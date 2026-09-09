# `IdrAgra_parameters.txt` guide

IdrAgra v2.3.0 compact feature demo

Anything that comes after "#" is a comment and is ignored by IdrAgra.  
Folder paths are relative to the directory from which IdrAgra is run.

## 0. Output controls

### `prt_all`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `y`
:::

prt_all: master switch for the complete standard output set [y/n].  
Individual switches below may be used to enable selected outputs when  
prt_all = n.

### `prt_annual`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `y`
:::

prt_annual: print annual summaries [y/n].

### `prt_yield`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `y`
:::

prt_yield: print crop-yield output [y/n].

### `prt_stp_irr`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `y`
:::

Time-step output switches [y/n].  
irrigation

### `prt_stp_cap_rise`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `y`
:::

capillary rise

### `prt_stp_deep_perc`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `y`
:::

deep percolation

### `prt_stp_et_act`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `y`
:::

actual evapotranspiration

### `prt_stp_et_pot`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `y`
:::

potential evapotranspiration

## 1. Input and output paths

### `OutputPath`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `.\\sim_results\\`
:::

OutputPath: folder in which simulation results are written.

### `InputPath`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `.\\spatial_data\\`
:::

InputPath: folder containing spatialized model inputs.

### `MeteoPath`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `.\\meteo_data\\`
:::

MeteoPath: folder containing meteorological station data.

### `MeteoFileName`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `weather_stations.dat`
:::

MeteoFileName: root-level file listing the meteorological input files.

### `PhenoPath`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `.\\crop_series\\`
:::

PhenoPath: folder containing crop phenology parameter files.

### `PhenoFileRoot`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `pheno_`
:::

PhenoFileRoot: common prefix used by phenology files associated with  
meteorological stations.

### `IrrMethPath`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `.\\irrmeth_data\\`
:::

IrrMethPath: folder containing irrigation-method files.

### `IrrMethFileName`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `irrmethods.txt`
:::

IrrMethFileName: file listing the irrigation-method input files.

### `WatSourPath`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `.\\watsour_data\\`
:::

WatSourPath: folder containing water-source input files.

## 2. Simulation mode, period and conditions

### `Mode`

:::{container} parameter-meta meta-keywords
**Type:** Integer  ·  **Default:** `2`
:::

Mode: irrigation simulation mode [0...4].  
  0 = no irrigation  
  1 = irrigation consumptions  
  2 = satisfy field-capacity requirements  
  3 = fixed irrigation volumes  
  4 = fixed irrigation dates and volumes read from file

### `InitialThetaFlag`

:::{container} parameter-meta meta-keywords
**Type:** Boolean  ·  **Default:** `false`
:::

InitialThetaFlag: read the initial soil-moisture state from an external  
file [T/F]. When F, IdrAgra generates the initial condition internally.

### `FinalThetaFlag`

:::{container} parameter-meta meta-required-soft
**Type:** Boolean  ·  **Required:** no code default
:::

FinalThetaFlag: save the final soil-moisture state for possible reuse as  
an initial condition in a subsequent simulation [T/F].

### `StartSimulation`

:::{container} parameter-meta meta-keywords
**Type:** Date  ·  **Default:** `date(29, 2, 1600, 2305507, 0)`
:::

StartSimulation / EndSimulation: inclusive simulation dates [dd/mm/yyyy].

### `EndSimulation`

:::{container} parameter-meta meta-keywords
**Type:** Date  ·  **Default:** `date(29, 2, 1600, 2305507, 0)`
:::

### `InitialConditionPath`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `.\\spatial_data\\`
:::

Initial-condition input. Used when InitialThetaFlag = T.

### `InitialCondition`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `IC_thetaI`
:::

### `FinalConditionPath`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `.\\sim_results\\`
:::

Final-condition output. Used when FinalThetaFlag = T.

### `FinalCondition`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `FC_thetaI`
:::

### `CapillaryFlag`

:::{container} parameter-meta meta-keywords
**Type:** Boolean  ·  **Default:** `false`
:::

CapillaryFlag: enable capillary-rise simulation [T/F].

### `SoilUseVarFlag`

:::{container} parameter-meta meta-keywords
**Type:** Boolean  ·  **Default:** `false`
:::

SoilUseVarFlag: allow soil use to vary between years [T/F].

## 3. Meteorology, land use, and sowing

### `MeteoStatTotNum`

:::{container} parameter-meta meta-keywords
**Type:** Integer  ·  **Default:** `1`
:::

MeteoStatTotNum: total number of meteorological stations.

### `MeteoStatWeightNum`

:::{container} parameter-meta meta-keywords
**Type:** Integer  ·  **Default:** `1`
:::

MeteoStatWeightNum: number of nearest stations used in the weighted  
meteorological/phenological calculations.

### `InterpolateTemperature`

:::{container} parameter-meta meta-keywords
**Type:** Boolean  ·  **Default:** `true`
:::

Whether simulation cells receive interpolated weather data (T, default&recommended) or the closest station's value (F, sometimes useful for rainfall).

### `InterpolateRain`

:::{container} parameter-meta meta-keywords
**Type:** Boolean  ·  **Default:** `true`
:::

### `InterpolateHumidity`

:::{container} parameter-meta meta-keywords
**Type:** Boolean  ·  **Default:** `true`
:::

### `InterpolateWind`

:::{container} parameter-meta meta-keywords
**Type:** Boolean  ·  **Default:** `true`
:::

### `InterpolateRadiation`

:::{container} parameter-meta meta-keywords
**Type:** Boolean  ·  **Default:** `true`
:::

### `SoilUsesNum`

:::{container} parameter-meta meta-keywords
**Type:** Integer  ·  **Default:** `1`
:::

SoilUsesNum: number of soil-use classes included in each phenological  
series.

### `SimulatedSoilUses`

:::{container} parameter-meta meta-required-soft
**Type:** Integer array  ·  **Required:** no code default
:::

SimulatedSoilUses: IDs of the soil-use classes to simulate.

### `RandSowDaysSym`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `symmetric`
:::

RandSowDaysSym: shape of the sowing-date randomization window.

### `RandSowDaysWind`

:::{container} parameter-meta meta-keywords
**Type:** Integer  ·  **Default:** `0`  ·  **Unit:** days
:::

RandSowDaysWind: half-width or extent, in days, of the sowing-date  
randomization window.

### `Repeatable`

:::{container} parameter-meta meta-keywords
**Type:** Boolean  ·  **Default:** `true`
:::

Repeatable: use a repeatable random sequence, allowing identical results  
across equivalent runs [T/F].

### `Forecast_day`

:::{container} parameter-meta meta-keywords
**Type:** Integer  ·  **Default:** `5`
:::

Forecast_day: forecast horizon or reference day used by forecast-related  
calculations.

## 4. Periodic output

### `MonthlyFlag`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `monthly`
:::

MonthlyFlag: periodic-output schedule. In this demo, "monthly" requests  
monthly output. WeekDay and StartDate/EndDate/DeltaDate below are inactive  
unless a weekly or periodic schedule is selected.

### `WeekDay`

:::{container} parameter-meta meta-keywords
**Type:** String or Integer  ·  **Default:** `monday`
:::

WeekDay: used only when MonthlyFlag selects a weekly schedule.

### `StartDate`

:::{container} parameter-meta meta-keywords
**Type:** Integer  ·  **Default:** `10`
:::

StartDate / EndDate: used only when MonthlyFlag selects a periodic schedule;  
first and last day of year [1...366] included in that schedule.

### `EndDate`

:::{container} parameter-meta meta-keywords
**Type:** Integer  ·  **Default:** `100`
:::

### `DeltaDate`

:::{container} parameter-meta meta-keywords
**Type:** Integer  ·  **Default:** `30`  ·  **Unit:** d
:::

DeltaDate: used only with a periodic schedule; interval between outputs [d].

## 5. Soil-water, runoff, and layer parameters

### `01q_eva`

:::{container} parameter-meta meta-keywords
**Type:** Real  ·  **Default:** `0.575118`
:::

Lower and upper Ksat calibration breakpoints for the evaporative layer.  
IdrAgra does not calculate these percentiles: the supplied values are used  
directly to interpolate the irrigation-related percolation-booster  
coefficients. They affect irrigated simulations (Modes 1...4).  
01q_eva: lower breakpoint; 09q_eva: upper breakpoint.

### `09q_eva`

:::{container} parameter-meta meta-keywords
**Type:** Real  ·  **Default:** `8.026400`
:::

### `01q_trasp`

:::{container} parameter-meta meta-keywords
**Type:** Real  ·  **Default:** `0.472116`
:::

Lower and upper Ksat calibration breakpoints for the transpirative layer.  
As above, these are direct calibration inputs rather than statistics  
calculated by IdrAgra.  
01q_trasp: lower breakpoint; 09q_trasp: upper breakpoint.

### `09q_trasp`

:::{container} parameter-meta meta-keywords
**Type:** Real  ·  **Default:** `7.706101`
:::

### `zEvap`

:::{container} parameter-meta meta-keywords
**Type:** Real  ·  **Default:** `0.10`  ·  **Unit:** m
:::

zEvap: evaporative-layer depth [m].

### `zRoot`

:::{container} parameter-meta meta-keywords
**Type:** Real  ·  **Default:** `0.90`  ·  **Unit:** m
:::

zRoot: transpirative/root-zone depth [m].

### `LambdaCN`

:::{container} parameter-meta meta-keywords
**Type:** Real  ·  **Default:** `0.2`  ·  **Unit:** -
:::

LambdaCN: initial-abstraction ratio used by the Curve Number runoff  
formulation [-].

### `lim_prec`

:::{container} parameter-meta meta-keywords
**Type:** Real  ·  **Default:** `5.0`  ·  **Unit:** mm
:::

lim_prec: minimum precipitation threshold used by rainfall/runoff  
calculations [mm].

### `h_maxpond`

:::{container} parameter-meta meta-keywords
**Type:** Real  ·  **Default:** `0.0D0`  ·  **Unit:** mm
:::

h_maxpond: general/out-of-season maximum surface-ponding depth [mm].  
During the irrigation season, irrigation-method values override this value.

### `fc_ratio`

:::{container} parameter-meta meta-keywords
**Type:** Real  ·  **Default:** `1.0D0`  ·  **Unit:** -
:::

fc_ratio: field-capacity target multiplier [-]. Used only in Mode 2;  
changing it has no effect in this demo while Mode = 1.

## 6. Irrigation season and source files

### `StartIrrSeason`

:::{container} parameter-meta meta-keywords
**Type:** Integer  ·  **Default:** `91`
:::

StartIrrSeason / EndIrrSeason: default irrigation-method season [1...366].  
Values specified by an individual irrigation-method file override these.

### `EndIrrSeason`

:::{container} parameter-meta meta-keywords
**Type:** Integer  ·  **Default:** `304`
:::

### `WatSources_fn`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `watsources.txt`
:::

Water-source and irrigation-network input filenames. Most source/diversion  
files are used by Mode 1; sched_irr_fn is used only by Mode 4.

### `mon_sources_i_div_fn`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `monit_sources_i.txt`
:::

### `mon_sources_ii_div_fn`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `monit_sources_ii.txt`
:::

### `int_reuse_div_fn`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `int_reuse.txt`
:::

### `cr_sources_list_fn`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `cr_sources.txt`
:::

### `IrrDistr_list_fn`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `irr_districts.txt`
:::

### `sched_irr_fn`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `scheduled_irrigation.txt`
:::

## 7. DTx settings (experimental; inactive in this demo)

### `DTxMode`

:::{container} parameter-meta meta-keywords
**Type:** String  ·  **Default:** `none`
:::

DTx support is incomplete. Keep all values in this block even when the mode  
is "none", because the current code does not provide safe defaults for every  
internal DTx dimension.

DTxMode: experimental DTx calculation mode [none/analysis/application].  
  none        = disable DTx calculations  
  analysis    = accumulate samples and fit gamma-distribution parameters  
  application = export current-run accumulated deficit maps; despite its  
                name, it does NOT load or apply relationships from analysis

### `DTxNumXs`

:::{container} parameter-meta meta-required-soft
**Type:** Integer  ·  **Required:** no code default
:::

DTxNumXs: number of DTx integration periods listed in DTx_x.

### `DTx_x`

:::{container} parameter-meta meta-required-soft
**Type:** Integer array  ·  **Required:** no code default  ·  **Unit:** days
:::

DTx_x: integration periods in days. For example, DT10 represents the  
transpirative deficit accumulated over 10 days.

### `DTxDeltaDate`

:::{container} parameter-meta meta-required-soft
**Type:** Integer  ·  **Required:** no code default  ·  **Unit:** days
:::

DTxDeltaDate: interval, in days, between DTx calculations.

### `DTxDelayDays`

:::{container} parameter-meta meta-required-soft
**Type:** Integer  ·  **Required:** no code default  ·  **Unit:** days
:::

DTxDelayDays: delay, in days from the start of the year, before DTx  
calculations begin.

### `DTxMinCard`

:::{container} parameter-meta meta-keywords
**Type:** Integer  ·  **Default:** `0`
:::

DTxMinCard: minimum sample cardinality required for a valid estimate in  
analysis mode. This two-year demo cannot normally satisfy a value of 3;  
the value is retained because DTx is disabled.
