import argparse
from dataclasses import dataclass
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np


READY_TOKEN = "IDRAGRA_COUPLING_READY"
CONTINUE_TOKEN = "IDRAGRA_COUPLING_CONTINUE"
DEFAULT_IDRAGRA_PARAMETERS_FILE = "idragra_parameters.txt"
DEFAULT_MODFLOW_PARAMETERS_FILE = "modflow_parameters.txt"
DEFAULT_IDRAGRA_INPUT_DIRECTORY = "spatial_data"
DEFAULT_DOMAIN_FILE = "domain.asc"
DEFAULT_COUPLING_DIRECTORY = "modflow_exchange"


class CouplerError(RuntimeError):
    pass


class CommandLineArgs(argparse.Namespace):
    working_dir: Path
    idragra_exe: Path
    idragra_args: list[str]


@dataclass(frozen=True)
class ModflowConfiguration:
    library: Path
    workspace: Path
    period_days: int = 7
    recharge_package: str | None = None
    well_package: str | None = None

    @classmethod
    def read(cls, path: Path) -> "ModflowConfiguration":
        values = _read_params_file(path)
        try:
            library = values["libmf6"]
            workspace = values["modflowworkspace"]
        except KeyError as exc:
            raise CouplerError(f"{path}: missing required setting {exc.args[0]}") from exc
        try:
            period_days = int(values["couplingperioddays"])
        except KeyError as exc:
            raise CouplerError(f"{path}: missing required setting CouplingPeriodDays") from exc
        except ValueError as exc:
            raise CouplerError(f"{path}: CouplingPeriodDays must be an integer") from exc
        if period_days <= 0:
            raise CouplerError(f"{path}: CouplingPeriodDays must be positive")
        return cls(
            library = _resolve_path_from(path.parent, library),
            workspace = _resolve_path_from(path.parent, workspace),
            period_days = period_days,
            recharge_package = values.get("rechargepackage"),
            well_package = values.get("wellpackage"),
        )


class AsciiGrid:
    def __init__(self, header: list[tuple[str, str]], data: np.ndarray):
        self.header = header
        self.data = np.asarray(data, dtype=float)

    @property
    def nodata(self) -> float:
        return float(dict((key.lower(), value) for key, value in self.header)["nodata_value"])

    @property
    def cellsize(self) -> float:
        return float(dict((key.lower(), value) for key, value in self.header)["cellsize"])

    @property
    def valid(self) -> np.ndarray:
        return self.data != self.nodata

    @classmethod
    def read(cls, path: Path) -> "AsciiGrid":
        with path.open("r", encoding="utf-8") as stream:
            header: list[tuple[str, str]] = []
            for _ in range(6):
                parts = stream.readline().split()
                if len(parts) < 2:
                    raise CouplerError(f"Invalid ESRI ASCII header in {path}")
                header.append((parts[0], parts[1]))
            data = np.loadtxt(stream, dtype=float)

        values = dict((key.lower(), value) for key, value in header)
        expected = (int(values["nrows"]), int(values["ncols"]))
        if data.shape != expected:
            raise CouplerError(f"{path} has shape {data.shape}; expected {expected}")
        return cls(header, data)

    def write(self, path: Path) -> None:
        temporary = path.with_suffix(path.suffix + ".tmp")
        with temporary.open("w", encoding="utf-8", newline="\n") as stream:
            for key, value in self.header:
                stream.write(f"{key} {value}\n")
            np.savetxt(stream, self.data, fmt="%.10g")
        temporary.replace(path)


@dataclass(frozen=True)
class GridLayout:
    header: list[tuple[str, str]]
    shape: tuple[int, int]
    cellsize: float

    @classmethod
    def from_grid(cls, grid: AsciiGrid) -> "GridLayout":
        return cls(grid.header.copy(), grid.data.shape, grid.cellsize)

    @property
    def size(self) -> int:
        return self.shape[0] * self.shape[1]


class IdrAgraProcess:
    def __init__(self, executable: Path, run_dir: Path, executable_args: list[str], coupling_dir: Path, period_days: int):
        self.coupling_dir = coupling_dir

        # Make coupling_dir and period_days available to IdrAgra as environmental variables
        environment = os.environ.copy()
        environment["IDRAGRA_COUPLING_DIR"] = str(coupling_dir)
        environment["IDRAGRA_COUPLING_DAYS"] = str(period_days)

        # Launches IdrAgra as a subprocess, redirecting its stdin and stdout to pipes for communication
        self.process = subprocess.Popen(
            [str(executable), *executable_args],
            cwd=run_dir,
            env=environment,
            stdin=subprocess.PIPE,    # Redirects anything written to self.process.stdin to IdrAgra's standard input
            stdout=subprocess.PIPE,   # Redirects all IdrAgra standard output to self.process.stdout for reading
            stderr=subprocess.STDOUT, # Adds IdrAgra's standard error to self.process.stdout as well
            text=True,
            bufsize=1,
            errors="replace",
        )

    def wait_for_idragra_flows(self) -> tuple[int, int, AsciiGrid, AsciiGrid]:
        if self.process.stdout is None:
            raise CouplerError("IdrAgra stdout is unavailable")
        while True:
            # Keep reading lines from IdrAgra's output
            line = self.process.stdout.readline()

            # If IdrAgra's stdout is closed, we receive an empty string (note that a purposefully empty line is received as "\n", not "")
            if line == "":
                code = self.process.poll()
                raise CouplerError(f"IdrAgra exited with code {code} before the next handshake")

            # Because Fortran's stdout no longer points to the console, it is up to python to print it to screen:
            print(line, end="")

            # Look for the handshake line that indicates IdrAgra is ready to exchange data
            stripped = line.strip()
            if not stripped.startswith(READY_TOKEN + " "):
                continue

            # When found, check it for errors and read the idragra-generated .asc files
            fields = stripped.split()
            if len(fields) != 3:
                raise CouplerError(f"Malformed IdrAgra handshake: {stripped}")
            period_index, period_days = int(fields[1]), int(fields[2])
            suffix = f"{period_index:06d}.asc"
            h_net_perc = AsciiGrid.read(self.coupling_dir / f"net_percolation_{suffix}")
            h_pumping = AsciiGrid.read(self.coupling_dir / f"pumping_{suffix}")
            return period_index, period_days, h_net_perc, h_pumping

    # Writes the updated water table to an .asc file and gives control back to IdrAgra
    def continue_idragra_sim(self, period_index: int, water_table: AsciiGrid) -> None:
        water_table.write(self.coupling_dir / f"water_table_{period_index:06d}.asc")
        if self.process.stdin is None:
            raise CouplerError("IdrAgra stdin is unavailable")
        self.process.stdin.write(f"{CONTINUE_TOKEN} {period_index}\n") # This is read by IdrAgra as if it were typed on the console (standard input)
        self.process.stdin.flush()

    # Waits for IdrAgra to finish; meanwhile keeps printing its stdout to the console
    def finish_simulation(self) -> None:
        if self.process.stdout is not None:
            for line in self.process.stdout:
                print(line, end="")
        code = self.process.wait()
        if code != 0:
            raise CouplerError(f"IdrAgra exited with code {code}")

    # Terminates the IdrAgra process if something went wrong elsewhere (e.g. MODFLOW failure)
    def terminate(self) -> None:
        if self.process.poll() is None:
            print("Forcibly terminating the IdrAgra process...")
            self.process.kill()
            self.process.wait()


class CouplingCallback:
    def __init__(
        self,
        idragra: IdrAgraProcess,
        callbacks: Any,
        rch_name: str | None,
        wel_name: str | None,
        grid_layout: GridLayout,
    ):
        self.idragra = idragra
        self.callbacks = callbacks
        self.rch_name = rch_name
        self.wel_name = wel_name
        self.grid_layout = grid_layout
        self.period_index = 0
        self.period_days = 0
        self.model: Any = None
        self.rch_package: Any = None
        self.wel_package: Any = None

    # This is the function that gets called whenever the MODFLOW api uses "callback(simulation, step)".
    # Here we configure it to do different things depending on said step.
    def __call__(self, simulation: Any, callback_step: Any) -> None:

        # Called once before the solving begins
        if callback_step == self.callbacks.initialize:
            self.model = simulation.get_model()
            _validate_modflow_grid(self.model, self.grid_layout)
            self.rch_package = _select_package(self.model, self.rch_name, "rch")
            self.wel_package = _select_package(self.model, self.wel_name, "wel")
            return

        # Called at the beginning of each stress period (import recharge and pumping data simulated by IdrAgra)
        if callback_step == self.callbacks.stress_period_start:
            period_index, period_days, h_net_perc, h_pumping = self.idragra.wait_for_idragra_flows()

            # Check that IdrAgra's period index (1-based) matches MODFLOW's stress period index (0-based)
            expected_index = int(simulation.kper) + 1
            if period_index != expected_index:
                raise CouplerError(f"IdrAgra period {period_index}; MODFLOW expected {expected_index}")

            recharge_m_per_day = h_net_perc.data / (1000.0 * period_days)
            pumping_m3_per_day = -(h_pumping.data / 1000.0) * self.grid_layout.cellsize**2 / period_days
            _set_mf_boundary_condition(self.rch_package, "recharge", h_net_perc, recharge_m_per_day, inactive_value=0.0)
            _set_mf_boundary_condition(self.wel_package, "q", h_pumping, pumping_m3_per_day, inactive_value=0.0)

            self.period_index = period_index
            self.period_days = period_days
            print(f"MODFLOW period {period_index}: importing IdrAgra flows for {period_days} day(s)")
            return

        # Called at the end of each stress period (gives control back to IdrAgra)
        if callback_step == self.callbacks.stress_period_end:
            heads = np.asarray(self.model.X, dtype=float).reshape(-1)
            top = np.asarray(self.model.dis.top.values, dtype=float).reshape(-1)
            if heads.size != self.grid_layout.size:
                raise CouplerError(
                    "Prototype requires one MODFLOW layer aligned with the complete IdrAgra grid; "
                    f"got {heads.size} heads and {self.grid_layout.size} IdrAgra cells"
                )
            if top.size == 1:
                top = np.full(heads.size, top.item())
            if top.size != heads.size:
                raise CouplerError(f"MODFLOW top has {top.size} values; expected {heads.size}")

            depth = np.maximum(top - heads, 0.0).reshape(self.grid_layout.shape)
            self.idragra.continue_idragra_sim(self.period_index, AsciiGrid(self.grid_layout.header, depth))



def main() -> int:

    # Use command line arguments to determine working directory and IdrAgra executable
    args = parse_command_line_args()
    run_dir = args.working_dir.resolve()
    executable = args.idragra_exe.resolve()

    # Strip leading "--" from the remaining arguments if present
    executable_args = args.idragra_args
    if executable_args[:1] == ["--"]:
        executable_args = executable_args[1:]

    try:
        # Ensure coupling is explicitly enabled in idragra_parameters.txt
        idragra_parameters = _parameter_file(run_dir, executable_args)
        if not _is_coupling_enabled(idragra_parameters):
            raise CouplerError(
                f"MODFLOW coupling is disabled in {idragra_parameters}; set UseModflowCoupling = T to enable it"
            )

        # Read domain.asc and store its information in GridLayout (data exchanged between the models must match this grid)
        domain_grid = AsciiGrid.read(_domain_file(run_dir, idragra_parameters))
        grid_layout = GridLayout.from_grid(domain_grid)

        # Read modflow_parameters.txt
        params_file_path = run_dir / DEFAULT_MODFLOW_PARAMETERS_FILE
        params = ModflowConfiguration.read(params_file_path)
    except (OSError, CouplerError) as exc:
        raise SystemExit(str(exc)) from exc

    # Ensure modflowapi is installed
    try:
        import modflowapi
    except ImportError as exc:
        raise SystemExit("Install the prototype dependencies with: pip install modflowapi numpy") from exc

    # Set up the coupling directory (it will store temporary files exchanged between IdrAgra and MODFLOW)
    coupling_dir = run_dir / DEFAULT_COUPLING_DIRECTORY
    coupling_dir.mkdir(parents=True, exist_ok=True)

    # Set up and launch the IdrAgra process
    idragra = IdrAgraProcess(
        executable.resolve(), run_dir, executable_args, coupling_dir.resolve(), params.period_days
    )

    # Set up the callback function that will be called by the MODFLOW API to exchange data with IdrAgra
    callback = CouplingCallback(
        idragra, modflowapi.Callbacks, params.recharge_package, params.well_package, grid_layout
    )

    try:
        modflowapi.run_simulation(str(params.library), str(params.workspace), callback, verbose=False)
        idragra.finish_simulation()
    finally:
        idragra.terminate()
    return 0



def parse_command_line_args() -> CommandLineArgs:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--working-dir", required=True, type=Path)
    parser.add_argument("--idragra-exe", required=True, type=Path)

    # Any additional argument is saved in "args.idragra_args" and is passed to IdrAgra to be interpreted as if it were on command line
    parser.add_argument(
        "idragra_args",
        nargs=argparse.REMAINDER,
        help="arguments passed to IdrAgra (place them after --, for example: -- -f custom_parameters.txt)",
    )
    return parser.parse_args(namespace=CommandLineArgs())

# Get the path to idragra_parameters.txt (or custom txt if specified in command line arguments)
def _parameter_file(run_dir: Path, executable_args: list[str]) -> Path:
    for index, argument in enumerate(executable_args):
        if argument.lower() in {"-f", "-filename"}:
            if index + 1 == len(executable_args):
                raise CouplerError(f"{argument} requires an IdrAgra parameter filename")
            return _resolve_path_from(run_dir, executable_args[index + 1])
    return (run_dir / DEFAULT_IDRAGRA_PARAMETERS_FILE).resolve()

# Locate IdrAgra's domain.asc
def _domain_file(run_dir: Path, idragra_parameters: Path) -> Path:
    values = _read_params_file(idragra_parameters)
    input_directory = values.get("inputpath", DEFAULT_IDRAGRA_INPUT_DIRECTORY).replace("\\", os.sep)
    return _resolve_path_from(run_dir, input_directory) / DEFAULT_DOMAIN_FILE

# Return an absolute path to a file given a base directory (if it isn't already absolute)
def _resolve_path_from(base_dir: Path, a_maybe_relative_path: str) -> Path:
    path = Path(a_maybe_relative_path)
    if not path.is_absolute():
        path = base_dir / path
    return path.resolve()

# Check whether "UseModflowCoupling = T" in idragra_parameters.txt
def _is_coupling_enabled(path: Path) -> bool:
    values = _read_params_file(path)
    usemodflowcoupling = values.get("usemodflowcoupling", "F")
    if usemodflowcoupling.strip().lower() == "t":
        return True
    return False

# Read an input text file in the "parameter = value" format and store its contents in the "values" dictionary
def _read_params_file(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open("r", encoding="utf-8-sig") as params_file:
        for line_number, raw_line in enumerate(params_file, 1):
            line = raw_line.partition("#")[0].strip() # Discard comments and leading/trailing whitespace
            if not line:
                continue
            if "=" not in line:
                raise CouplerError(f"{path}:{line_number}: expected 'name = value' format")
            parameter, value = (part.strip() for part in line.split("=", 1))
            if not parameter or not value:
                raise CouplerError(f"{path}:{line_number}: expected 'name = value' format")
            values[parameter.lower()] = value.strip("\"'") # Save the entry in the dictionary, using lowercase for the parameter name and stripping quotes if any
    return values


def _select_package(model: Any, name: str | None, package_type: str) -> Any:
    matches = {
        package_name: package
        for package_name, package in model.package_dict.items()
        if package.pkg_type.lower() == package_type
    }
    if name is not None:
        package = matches.get(name.lower())
        if package is None:
            available = ", ".join(matches) or "none"
            raise CouplerError(
                f"MODFLOW {package_type.upper()} package {name!r} was not found; available: {available}"
            )
        return package
    if len(matches) == 1:
        return next(iter(matches.values()))
    if not matches:
        raise CouplerError(f"MODFLOW model has no {package_type.upper()} package")
    raise CouplerError(
        f"MODFLOW model has multiple {package_type.upper()} packages ({', '.join(matches)}); "
        f"set {'RechargePackage' if package_type == 'rch' else 'WellPackage'}"
    )

# Make sure that the MODFLOW grid matches the IdrAgra domain.asc grid
def _validate_modflow_grid(model: Any, grid_layout: GridLayout) -> None:
    model_shape = tuple(int(value) for value in model.shape)
    expected_shape = (1, *grid_layout.shape)
    if model_shape != expected_shape:
        raise CouplerError(f"MODFLOW grid has shape {model_shape}; expected {expected_shape} from {DEFAULT_DOMAIN_FILE}")

    cell_area = np.asarray(model.dis.area.values, dtype=float).reshape(-1)
    if cell_area.size != grid_layout.size or not np.allclose(
        cell_area, grid_layout.cellsize**2, rtol=1e-9, atol=1e-9
    ):
        expected_area = grid_layout.cellsize**2
        raise CouplerError(f"MODFLOW cell areas must equal the {expected_area:g} m2 cell area from domain.asc")

    # Check that the MODFLOW IDOMAIN is positive for every cell in the complete IdrAgra grid rectangle
    # (i.e. MODFLOW simulates all cells, even if some are inactive in IdrAgra)
    idomain = np.asarray(model.dis.idomain.values).reshape(-1)
    if idomain.size != grid_layout.size or np.any(idomain <= 0):
        raise CouplerError("MODFLOW IDOMAIN must be positive for every cell in the complete IdrAgra grid rectangle")

# Currently used to set IdrAgra-prescribed recharge and pumping values to MODFLOW packages
def _set_mf_boundary_condition(
    package: Any, field: str, grid: AsciiGrid, values: np.ndarray, inactive_value: float | None = None
) -> None:
    stress_data = package.stress_period_data
    try:
        target = np.asarray(stress_data[field], dtype=float).reshape(-1)
    except (KeyError, TypeError, ValueError) as exc:
        raise CouplerError(f"MODFLOW package has no mutable '{field}' field") from exc

    flat = np.asarray(values, dtype=float).reshape(-1)
    valid = grid.valid.reshape(-1)
    if target.size != flat.size:
        raise CouplerError(
            f"MODFLOW {package.pkg_type.upper()} package has {target.size} entries; expected one for each of the "
            f"{flat.size} cells in the complete IdrAgra grid rectangle"
        )

    updated = target.copy()
    updated[valid] = flat[valid]

    # Cells that are inactive IdrAgra-side (e.g. outside the IdrAgra domain) receive the specified inactive_value
    if inactive_value is not None:
        updated[~valid] = inactive_value
    stress_data[field] = updated


if __name__ == "__main__":
    sys.exit(main())
