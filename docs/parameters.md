# IdrAgra_parameters.txt guide

`idragra_parameters.txt` is the main file which determines how the model behaves.

For example, it contains options and switches to set:
- the simulation mode (e.g. USE or NEED)
- the simulation start and end dates
- whether to simulate capillary uptake or not
- whether to run warmup or not
- the paths to input and output files
- which output files are printed and how frequently

and so on.

In this page we list all of the available settings. Note that most of them are technically optional, due to the model applying internal defaults, which are specified for each entry; settings missing a default are highligthed as required settings.

When launching the IdrAgra exe, the model looks for `idragra_parameters.txt` in that same directory. However, a different path and/or filename can be specified using the `-f` [command line argument](command_line_args.md).

## About the file structure

- Blank lines are ignored by the model and skipped.
- Anything that comes after "#" is a comment and is also ignored.
- Variables can appear in any order.
- The accepted structure is "VariableName = value".
- Both the variable name and its value are NOT case-sensitive
- Any number of spaces and tabs is accepted before and after the "=" symbol.
- It is not allowed to set multiple variables on the same line.
- Folder paths are relative to the directory from which IdrAgra is run (usually the exe's location), and must use `\\` as the delimiter.
- If a path includes spaces (not recommended), it must be enclosed in quotes (Path = "a path")
- Boolean variables accept only "T" (true) and "F" (false) as values.

{.parameter-list-group .parameter-list-first}
## Input and output paths

### `OutputPath`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `.\\sim_results\\`
:::

Folder in which simulation results are written.

### `InputPath`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `.\\spatial_data\\`
:::

Folder containing spatialized model inputs.

### `MeteoPath`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `.\\meteo_data\\`
:::

Folder containing meteorological station data.

### `MeteoFileName`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `weather_stations.dat`
:::

Root-level file listing the meteorological input files.

### `PhenoPath`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `.\\crop_series\\`
:::

Folder containing crop phenology parameter files.

### `PhenoFileRoot`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `pheno_`
:::

Common prefix used by phenology files associated with meteorological stations.

### `IrrMethPath`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `.\\irrmeth_data\\`
:::

Folder containing irrigation-method files.

### `IrrMethFileName`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `irrmethods.txt`
:::

File listing the irrigation-method input files.

### `WatSourPath`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `.\\watsour_data\\`
:::

Folder containing water-source input files.

{.parameter-list-group}
## Simulation mode, period and conditions

### `Mode`

:::{container} parameter-summary
**Type:** Integer  ·  **Default:** `2`
:::

Irrigation simulation mode [0...4].\
  0 = no irrigation  
  1 = irrigation consumptions  
  2 = satisfy field-capacity requirements  
  3 = fixed irrigation volumes  
  4 = fixed irrigation dates and volumes read from file

### `InitialThetaFlag`

:::{container} parameter-summary
**Type:** Boolean  ·  **Default:** `false`
:::

If true, the initial soil-moisture condition is read from the `InitialConditionPath\\InitialCondition.asc` file.\
If false, IdrAgra begins the simulation with all soils at field capacity but runs a warmup year reusing the first year of data.

{.parameter-dependent}
#### `InitialConditionPath`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `.\\spatial_data\\`
:::

{.parameter-dependent}
#### `InitialCondition`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `IC_thetaI`
:::

### `FinalThetaFlag`

:::{container} parameter-summary parameter-summary--required
**Type:** Boolean  ·  **Required:** no code default
:::

FinalThetaFlag: save the final soil-moisture state for possible reuse as an initial condition in a subsequent simulation [T/F].

{.parameter-dependent}
#### `FinalConditionPath`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `.\\sim_results\\`
:::

Final-condition output. Used when FinalThetaFlag = T.

{.parameter-dependent}
#### `FinalCondition`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `FC_thetaI`
:::

### `StartSimulation`

:::{container} parameter-summary
**Type:** Date  ·  **Default:** `date(29, 2, 1600, 2305507, 0)`
:::

StartSimulation / EndSimulation: inclusive simulation dates [dd/mm/yyyy].

### `EndSimulation`

:::{container} parameter-summary
**Type:** Date  ·  **Default:** `date(29, 2, 1600, 2305507, 0)`
:::

### `CapillaryFlag`

:::{container} parameter-summary
**Type:** Boolean  ·  **Default:** `false`
:::

CapillaryFlag: enable capillary-rise simulation [T/F].

### `SoilUseVarFlag`

:::{container} parameter-summary
**Type:** Boolean  ·  **Default:** `false`
:::

SoilUseVarFlag: allow soil use to vary between years [T/F].

{.parameter-list-group}
## Meteorology, land use, and sowing

### `MeteoStatTotNum`

:::{container} parameter-summary
**Type:** Integer  ·  **Default:** `1`
:::

MeteoStatTotNum: total number of meteorological stations.

### `MeteoStatWeightNum`

:::{container} parameter-summary
**Type:** Integer  ·  **Default:** `1`
:::

MeteoStatWeightNum: number of nearest stations used in the weighted  
meteorological/phenological calculations.

### `InterpolateTemperature`

:::{container} parameter-summary
**Type:** Boolean  ·  **Default:** `true`
:::

Whether simulation cells receive interpolated weather data (T, default&recommended) or the closest station's value (F, sometimes useful for rainfall).

### `InterpolateRain`

:::{container} parameter-summary
**Type:** Boolean  ·  **Default:** `true`
:::

### `InterpolateHumidity`

:::{container} parameter-summary
**Type:** Boolean  ·  **Default:** `true`
:::

### `InterpolateWind`

:::{container} parameter-summary
**Type:** Boolean  ·  **Default:** `true`
:::

### `InterpolateRadiation`

:::{container} parameter-summary
**Type:** Boolean  ·  **Default:** `true`
:::

### `SoilUsesNum`

:::{container} parameter-summary
**Type:** Integer  ·  **Default:** `1`
:::

SoilUsesNum: number of soil-use classes included in each phenological  
series.

### `SimulatedSoilUses`

:::{container} parameter-summary parameter-summary--required
**Type:** Integer array  ·  **Required:** no code default
:::

SimulatedSoilUses: IDs of the soil-use classes to simulate.

### `RandSowDaysSym`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `symmetric`
:::

RandSowDaysSym: shape of the sowing-date randomization window.

### `RandSowDaysWind`

:::{container} parameter-summary
**Type:** Integer  ·  **Default:** `0`  ·  **Unit:** days
:::

RandSowDaysWind: half-width or extent, in days, of the sowing-date  
randomization window.

### `Repeatable`

:::{container} parameter-summary
**Type:** Boolean  ·  **Default:** `true`
:::

Repeatable: use a repeatable random sequence, allowing identical results  
across equivalent runs [T/F].

### `Forecast_day`

:::{container} parameter-summary
**Type:** Integer  ·  **Default:** `5`
:::

Forecast_day: forecast horizon or reference day used by forecast-related  
calculations.

{.parameter-list-group}
## Periodic output

### `MonthlyFlag`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `monthly`
:::

MonthlyFlag: periodic-output schedule. In this demo, "monthly" requests  
monthly output. WeekDay and StartDate/EndDate/DeltaDate below are inactive  
unless a weekly or periodic schedule is selected.

### `WeekDay`

:::{container} parameter-summary
**Type:** String or Integer  ·  **Default:** `monday`
:::

WeekDay: used only when MonthlyFlag selects a weekly schedule.

### `StartDate`

:::{container} parameter-summary
**Type:** Integer  ·  **Default:** `10`
:::

StartDate / EndDate: used only when MonthlyFlag selects a periodic schedule;  
first and last day of year [1...366] included in that schedule.

### `EndDate`

:::{container} parameter-summary
**Type:** Integer  ·  **Default:** `100`
:::

### `DeltaDate`

:::{container} parameter-summary
**Type:** Integer  ·  **Default:** `30`  ·  **Unit:** d
:::

DeltaDate: used only with a periodic schedule; interval between outputs [d].

{.parameter-list-group}
## Soil-water, runoff, and layer parameters

### `01q_eva`

:::{container} parameter-summary
**Type:** Real  ·  **Default:** `0.575118`
:::

Lower and upper Ksat calibration breakpoints for the evaporative layer. IdrAgra does not calculate these percentiles: the supplied values are used directly to interpolate the irrigation-related percolation-booster coefficients. They affect irrigated simulations (Modes 1...4).
01q_eva: lower breakpoint; 09q_eva: upper breakpoint.

### `09q_eva`

:::{container} parameter-summary
**Type:** Real  ·  **Default:** `8.026400`
:::

### `01q_trasp`

:::{container} parameter-summary
**Type:** Real  ·  **Default:** `0.472116`
:::

Lower and upper Ksat calibration breakpoints for the transpirative layer.
As above, these are direct calibration inputs rather than statistics calculated by IdrAgra.
01q_trasp: lower breakpoint; 09q_trasp: upper breakpoint.

### `09q_trasp`

:::{container} parameter-summary
**Type:** Real  ·  **Default:** `7.706101`
:::

### `zEvap`

:::{container} parameter-summary
**Type:** Real  ·  **Default:** `0.10`  ·  **Unit:** m
:::

zEvap: evaporative-layer depth [m].

### `zRoot`

:::{container} parameter-summary
**Type:** Real  ·  **Default:** `0.90`  ·  **Unit:** m
:::

zRoot: transpirative/root-zone depth [m].

### `LambdaCN`

:::{container} parameter-summary
**Type:** Real  ·  **Default:** `0.2`  ·  **Unit:** -
:::

LambdaCN: initial-abstraction ratio used by the Curve Number runoff  
formulation [-].

### `lim_prec`

:::{container} parameter-summary
**Type:** Real  ·  **Default:** `5.0`  ·  **Unit:** mm
:::

lim_prec: minimum precipitation threshold used by rainfall/runoff  
calculations [mm].

### `h_maxpond`

:::{container} parameter-summary
**Type:** Real  ·  **Default:** `0.0D0`  ·  **Unit:** mm
:::

h_maxpond: general/out-of-season maximum surface-ponding depth [mm].
During the irrigation season, irrigation-method values override this value.

### `fc_ratio`

:::{container} parameter-summary
**Type:** Real  ·  **Default:** `1.0D0`  ·  **Unit:** -
:::

fc_ratio: field-capacity target multiplier [-]. Used only in Mode 2.

{.parameter-list-group}
## Irrigation season and source files

### `StartIrrSeason`

:::{container} parameter-summary
**Type:** Integer  ·  **Default:** `91`
:::

StartIrrSeason / EndIrrSeason: default irrigation-method season [1...366].
Values specified by an individual irrigation-method file override these.

### `EndIrrSeason`

:::{container} parameter-summary
**Type:** Integer  ·  **Default:** `304`
:::

### `WatSources_fn`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `watsources.txt`
:::

Water-source and irrigation-network input filenames. Most source/diversion files are used by Mode 1; sched_irr_fn is used only by Mode 4.

### `mon_sources_i_div_fn`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `monit_sources_i.txt`
:::

### `mon_sources_ii_div_fn`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `monit_sources_ii.txt`
:::

### `int_reuse_div_fn`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `int_reuse.txt`
:::

### `cr_sources_list_fn`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `cr_sources.txt`
:::

### `IrrDistr_list_fn`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `irr_districts.txt`
:::

### `sched_irr_fn`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `scheduled_irrigation.txt`
:::

{.parameter-list-group}
## Output controls

Decide which output files

### `prt_all`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `y`
:::

prt_all: master switch for the complete standard output set [y/n].
Individual switches below may be used to enable selected outputs when  
prt_all = n.

### `prt_annual`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `y`
:::

prt_annual: print annual summaries [y/n].

### `prt_yield`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `y`
:::

prt_yield: print crop-yield output [y/n].

### `prt_stp_irr`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `y`
:::

Time-step output switches [y/n].
irrigation

### `prt_stp_cap_rise`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `y`
:::

capillary rise

### `prt_stp_deep_perc`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `y`
:::

deep percolation

### `prt_stp_et_act`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `y`
:::

actual evapotranspiration

### `prt_stp_et_pot`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `y`
:::

potential evapotranspiration

{.parameter-list-group .parameter-list-last}
## DTx settings (deprecated)

### `DTxMode`

:::{container} parameter-summary
**Type:** String  ·  **Default:** `none`
:::

DTx support is incomplete. Keep all values in this block even when the mode is "none", because the current code does not provide safe defaults for every internal DTx dimension.

DTxMode: experimental DTx calculation mode [none/analysis/application].
  none        = disable DTx calculations  
  analysis    = accumulate samples and fit gamma-distribution parameters  
  application = export current-run accumulated deficit maps; despite its  
                name, it does NOT load or apply relationships from analysis

### `DTxNumXs`

:::{container} parameter-summary parameter-summary--required
**Type:** Integer  ·  **Required:** no code default
:::

DTxNumXs: number of DTx integration periods listed in DTx_x.

### `DTx_x`

:::{container} parameter-summary parameter-summary--required
**Type:** Integer array  ·  **Required:** no code default  ·  **Unit:** days
:::

DTx_x: integration periods in days. For example, DT10 represents the transpirative deficit accumulated over 10 days.

### `DTxDeltaDate`

:::{container} parameter-summary parameter-summary--required
**Type:** Integer  ·  **Required:** no code default  ·  **Unit:** days
:::

DTxDeltaDate: interval, in days, between DTx calculations.

### `DTxDelayDays`

:::{container} parameter-summary parameter-summary--required
**Type:** Integer  ·  **Required:** no code default  ·  **Unit:** days
:::

DTxDelayDays: delay, in days from the start of the year, before DTx  
calculations begin.

### `DTxMinCard`

:::{container} parameter-summary
**Type:** Integer  ·  **Default:** `0`
:::

DTxMinCard: minimum sample cardinality required for a valid estimate in analysis mode.
