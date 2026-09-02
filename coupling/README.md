# IdrAgra v2 / MODFLOW API prototype

This protoype uses the official MF6 API to couple IdrAgra and MODFLOW simulations.
The two models run indipendently and exchange data every x days (specified in `modflow_parameters.txt`), waiting for each other before continuing. Information exchange is done through .asc files stored in a temporary coupling folder:
- exchange_NNNNNN.asc (cumulative net percolation through the bottom of IdrAgra's profile over the period [mm, positive downwards])
- pumping_NNNNNN.asc (cumulative uptake from irrigation wells over the period [mm, negative])
- water_table_NNNNNN.asc (MODFLOW-updated water table depth [m from soil surface, positive])
MODFLOW receives exchange and pumping, updates the water table, and feeds it back to IdrAgra.

It is assumed that the user has a pre-existing MODFLOW 6 simulation ready to use, which should:
- Have one stress period per coupling period, including the potentially shorter final period
- Represent groundwater as a 1-layer grid, aligned with IdrAgra's active domain
- Include RCH and WEL packages already containing one entry per complete grid cell or
  one entry per active IdrAgra cell, in IdrAgra row-major order.
- Have initial heads representing the same condition as IdrAgra's `waterdepth.asc`
- `exchange_NNNNNN.asc` is treated as net recharge: deep percolation minus
  capillary rise. `pumping_NNNNNN.asc` is converted to negative WEL flow.

## Configuration

Download the MODFLOW 6.7.0 windows zip package and exctract its libmf6.dll to a known location

Install the Python dependencies numpy and modflowapi:
    python -m pip install modflowapi numpy

In the IdrAgra working directory:
1. Set `UseModflowCoupling = T` in `idragra_parameters.txt` (also make sure that `CapillaryFlag = T`)
2. Copy `modflow_parameters_example.txt` to the working directory as `modflow_parameters.txt` and adjust its values.

## Run

NOTE: due to the stdin channel being needed for the exchange of signals between the models, the coupled run doesn't wait for human input before overwriting output folders. Proceed with caution!

The coupled run can exclusively be launched via the modflow_coupler.py script like so:
    python coupling/modflow_coupler.py --working-dir C:/path/to/idragra_project --idragra-exe C:/path/to/idragra_latest.exe

For this repository layout:
    python coupling/modflow_coupler.py --working-dir demo --idragra-exe release/idragra_latest.exe

Normal IdrAgra command line arguments are still usable by adding them after an empty "--", e.g.:
    python coupling/modflow_coupler.py --working-dir demo --idragra-exe release/idragra_latest.exe -- -verbose -f my_idragra_parameters.txt

## Additional notes

- If IdrAgra is set to run a warmup year (i.e. `InitialThetaFlag = F`, recommended) MODFLOW's configuration must include that year
- This prototype is currently only compatible with IdrAgra's "raster mode", i.e. that in which .asc file actually represent the spatial distribution of cells
- Note that if MODFLOW were to prescribe a very shallow / above-ground water table, IdrAgra internally enforces a minimum depth equal to its first layer's thickness