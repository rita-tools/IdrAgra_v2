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
            h_net_perc = AsciiGrid.read(self.coupling_dir / f"exchange_{suffix}")
            h_pumping = AsciiGrid.read(self.coupling_dir / f"pumping_{suffix}")
            if h_net_perc.data.shape != h_pumping.data.shape:
                raise CouplerError("Exchange and pumping grids have different shapes")
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
    def __init__(self, idragra: IdrAgraProcess, callbacks: Any, rch_name: str | None, wel_name: str | None):
        self.idragra = idragra
        self.callbacks = callbacks
        self.rch_name = rch_name
        self.wel_name = wel_name
        self.period_index = 0
        self.period_days = 0
        self.h_net_perc: AsciiGrid | None = None
        self.model: Any = None
        self.rch_package: Any = None
        self.wel_package: Any = None

    # This is the function that gets called whenever the MODFLOW api uses "callback(simulation, step)".
    # Here we configure it to do different things depending on said step.
    def __call__(self, simulation: Any, callback_step: Any) -> None:

        # Called once before the solving begins
        if callback_step == self.callbacks.initialize:
            self.model = simulation.get_model()
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
            pumping_m3_per_day = -(h_pumping.data / 1000.0) * h_net_perc.cellsize**2 / period_days
            _set_boundary_field(self.rch_package, "recharge", h_net_perc, recharge_m_per_day)
            _set_boundary_field(self.wel_package, "q", h_pumping, pumping_m3_per_day)

            self.period_index = period_index
            self.period_days = period_days
            self.h_net_perc = h_net_perc
            print(f"MODFLOW period {period_index}: importing IdrAgra flows for {period_days} day(s)")
            return

        # Called at the end of each stress period (gives control back to IdrAgra)
        if callback_step == self.callbacks.stress_period_end:

            # Safety check
            if self.h_net_perc is None:
                raise CouplerError("MODFLOW reached period end without IdrAgra flows - was callback not called at period start?")

            
            heads = np.asarray(self.model.X, dtype=float).reshape(-1)
            top = np.asarray(self.model.dis.top.values, dtype=float).reshape(-1)
            if heads.size != self.h_net_perc.data.size:
                raise CouplerError(
                    "Prototype requires one MODFLOW layer aligned with the complete IdrAgra grid; "
                    f"got {heads.size} heads and {self.h_net_perc.data.size} IdrAgra cells"
                )
            if top.size == 1:
                top = np.full(heads.size, top.item())
            if top.size != heads.size:
                raise CouplerError(f"MODFLOW top has {top.size} values; expected {heads.size}")

            depth = np.maximum(top - heads, 0.0).reshape(self.h_net_perc.data.shape)
            depth[~self.h_net_perc.valid] = self.h_net_perc.nodata
            self.idragra.continue_idragra_sim(self.period_index, AsciiGrid(self.h_net_perc.header, depth))
            self.h_net_perc = None



def main() -> int:

    # Use command line arguments to determine working directory and IdrAgra executable
    args = parse_command_line_args()
    run_dir = args.working_dir.resolve()
    executable = args.idragra_exe.resolve()

    # Strip leading "--" from the remaining arguments if present
    executable_args = args.idragra_args
    if executable_args[:1] == ["--"]:
        executable_args = executable_args[1:]

    # Ensure coupling is explicitly enabled in idragra_parameters.txt and read modflow_parameters.txt
    try:
        idragra_parameters = _parameter_file(run_dir, executable_args)
        if not _is_coupling_enabled(idragra_parameters):
            raise CouplerError(
                f"MODFLOW coupling is disabled in {idragra_parameters}; set UseModflowCoupling = T to enable it"
            )
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
        idragra, modflowapi.Callbacks, params.recharge_package, params.well_package
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


def _values_for_boundary(grid: AsciiGrid, values: np.ndarray, boundary_size: int) -> np.ndarray:
    flat = np.asarray(values, dtype=float).reshape(-1)
    valid = grid.valid.reshape(-1)
    if boundary_size == flat.size:
        return np.where(valid, flat, 0.0)
    if boundary_size == int(valid.sum()):
        return flat[valid]
    raise CouplerError(
        f"MODFLOW boundary has {boundary_size} entries, but the IdrAgra grid has "
        f"{flat.size} total and {int(valid.sum())} active cells"
    )


def _set_boundary_field(package: Any, field: str, grid: AsciiGrid, values: np.ndarray) -> None:
    stress_data = package.stress_period_data
    try:
        target = stress_data[field]
    except (KeyError, TypeError, ValueError) as exc:
        raise CouplerError(f"MODFLOW package has no mutable '{field}' field") from exc
    stress_data[field] = _values_for_boundary(grid, values, target.size)


if __name__ == "__main__":
    sys.exit(main())
