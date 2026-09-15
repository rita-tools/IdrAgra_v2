# Project structure

In this section we look at how an IdrAgra working folder is structured.



## Input files

Input files make up the majority of the project's data. By default, they are organized in descriptive folders:

- **{ref}`spatial_data <parameter-inputpath>`**, containing .asc maps, i.e. anything that is spatialized at the {ref}`field <simulation-unit>` level;
- **{ref}`meteo_data <parameter-meteopath>`**, containing weather time series (one file per station);
- **{ref}`crop_series <parameter-phenopath>`**, containing crop parameters & time series, spatialized at the station level;
- **{ref}`irrmeth_data <parameter-irrmethpath>`**, containing each irrigation method's parameters;
- **{ref}`watsour_data <parameter-watsourpath>`**, containing water-source and scheduled-irrigation files;
- **{ref}`sim_results <parameter-outputpath>`**, containing generated outputs.

These are the default names. Each folder can be changed through its linked setting in `idragra_parameters.txt`.

...

(water-table-input-files)=
### Water table depth

If the user sets {ref}`CapillaryFlag <parameter-capflag>` to true, the model requires one or more .asc files describing the depth of the water table from the soil surface [m] to be placed in the {ref}`spatial_data <parameter-inputpath>` folder, following a strict naming convention:

- `waterdepth.asc` is read at the beginning of the simulation, and is always required.
- `waterdepth_yyyy_<doy>.asc` provides a recurring update on a day of year; `waterdepth_<year>_<doy>.asc` provides a year-specific update. Both forms are optional.

IdrAgra treats the received water table depths as static until they get updated by reading a new daily map. \
It is not necessary to provide a map for each day of the simulation - IdrAgra will simply continue the simulation using the last map it has read.

:::{container} llm-review-note
**LLM-authored draft — review required.** The practical inventory below was drafted from the 2025 manuals, the bundled demo, and the current file-reading code. Review it before relying on it as a complete input specification.
:::

(working-folder-layout)=
## Working-folder layout

IdrAgra locates each group of files through the path settings in `idragra_parameters.txt`. A typical project using the default names contains:

```text
project/
├── idragra_parameters.txt
├── weather_stations.dat
├── spatial_data/       # ESRI ASCII grids and rice_soilparam.txt
├── meteo_data/         # one daily weather file per station
├── crop_series/        # CanopyRes.dat and one phenology folder per station
├── irrmeth_data/       # irrigation-method list and method files
├── watsour_data/       # source, district, diversion, and schedule files
└── sim_results/        # generated outputs
```

Alternative names work when {ref}`InputPath <parameter-inputpath>`, {ref}`MeteoPath <parameter-meteopath>`, {ref}`PhenoPath <parameter-phenopath>`, {ref}`IrrMethPath <parameter-irrmethpath>`, {ref}`WatSourPath <parameter-watsourpath>`, and {ref}`OutputPath <parameter-outputpath>` point to the matching folders. When copying a parameter file between projects, copy the associated folder structure or update these paths.

(spatial-ascii-grids)=
## Spatial ASCII grids

Spatial inputs are ESRI ASCII grids. They must share alignment and resolution and cover the active domain. `domain.asc` defines which cells participate. The required filenames are:

| Group | Required root names |
|---|---|
| Domain | `domain` |
| Layer I water properties | `ThetaI_sat`, `ThetaI_FC`, `ThetaI_WP`, `ThetaI_r`, `Ksat_I`, `N_I` |
| Layer II water properties | `ThetaII_sat`, `ThetaII_FC`, `ThetaII_WP`, `ThetaII_r`, `Ksat_II`, `N_II` |
| Runoff | `slope`, `hydr_cond`, `hydr_group` |
| Land use and station weights | `soiluse`, `meteo_1`, ..., `meteo_N` |
| Irrigation, depending on mode | `irr_meth`, `appl_eff`, `irr_units`, `conv_eff` |
| Capillary rise, when enabled | `waterdepth`, `CapRisePar_a3`, `CapRisePar_a4`, `CapRisePar_b1`, `CapRisePar_b2`, `CapRisePar_b3`, `CapRisePar_b4` |
| Optional | `shapearea`, `irandom` |

Add the `.asc` extension to each root name shown above. Store these files in {ref}`InputPath <parameter-inputpath>`.

A conventional ESRI ASCII header is:

```text
ncols         100
nrows         80
xllcorner     500000
yllcorner     4990000
cellsize      250
NODATA_value  -9999
```

All grids must match the domain's row and column counts, lower-left origin, and cell size. NoData cells in a required parameter grid are excluded from the simulation domain.

:::{container} manual-code-divergence
**Validation gap to review.** IdrAgra does not currently warn when an input grid has a different `cellsize`. Check the resolution before running because IdrAgra does not resample input grids.
:::

:::{container} manual-code-divergence
**Manual/code divergence to review.** The installation manual describes `0` as the value outside `domain.asc`, but zero is treated as an active value. Use the grid's declared `NODATA_value` outside the domain.
:::

Enter slope values in `slope.asc` as percentages. Land-use IDs range from `1` to {ref}`SoilUsesNum <parameter-soilusesnum>`. Use irrigation-method ID `0` for no method, or a positive ID listed in the irrigation-method file.

### Time-varying maps

When {ref}`SoilUseVarFlag = T <parameter-soilusevarflag>`, yearly land-use and irrigation-method maps use:

```text
soiluse_<year>.asc
irr_meth_<year>.asc
appl_eff_<year>.asc   # read in Modes 2 and 4 when land use varies
```

The first simulation year must be present. The same naming convention is used when the maps are refreshed at the start of later years.

Water-table depth uses `waterdepth.asc` as its initial map. While {ref}`CapillaryFlag <parameter-capflag>` is enabled, IdrAgra looks for two optional update forms each day:

```text
waterdepth_yyyy_<doy>.asc   # recurring update in every year
waterdepth_<year>_<doy>.asc # update for one specific year
```

If both exist on the same day, the year-specific file takes precedence. A new map remains active until another update is found. Write the day of year as an unpadded integer, for example `waterdepth_2022_91.asc`.

:::{container} manual-code-divergence
**Manual/code divergence to review.** One map table in the installation manual uses `meth_eff_<year>.asc`. IdrAgra expects `appl_eff_<year>.asc`.
:::

(weather-input-files)=
## Weather inputs

{ref}`MeteoFileName <parameter-meteofilename>` points to a root-level station-list file. Its structure is:

```text
StatNum = 2
Table =
FileName X Y
station_west.dat 500250.0 5000750.0
station_east.dat 502000.0 5001000.0
EndTable =
```

`StatNum` must equal {ref}`MeteoStatTotNum <parameter-meteostattotnum>`, and the table must contain exactly that many station records. Coordinates must use the same projected coordinate system as the spatial grids.

Each station file begins with four header lines followed by one row per calendar day:

```text
Id station: 137, location: example
45.444092 102.0
01/01/2021 -> 31/12/2022
T_max T_min P_tot U_max U_min V_med RG_CORR
```

The second line contains latitude in degrees and elevation in metres. The daily columns are maximum and minimum air temperature, precipitation, maximum and minimum relative humidity, wind velocity, and solar radiation. Use a numeric filename stem as the station ID, or provide a positive ID in the first header line. IDs must be unique, and every station series must cover the requested simulation period.

The spatial weights are named with the lowercase root `meteo_1.asc`, ..., `meteo_N.asc`. Windows file lookup is normally case-insensitive, but projects intended to run on Linux should preserve this exact case even though the manual sometimes prints `Meteo_1.asc`.

:::{container} manual-code-divergence
**Manual/code divergence to review.** The installation manual restricts weather filenames to numeric station codes. Names such as `station_west.dat` also work when the first header line contains a positive, unique station ID.
:::

(selected-cell-input)=
## Selected-cell input

The presence of `cells.txt` in the working directory enables selected-cell outputs. The file contains a count followed by cell identifiers and one-based row and column indices:

```text
NCells = 2
Table =
ID Row Col
1 20 60
2 32 41
EndTable =
```

`NCells` must match the number of table rows. These are grid indices, not projected coordinates. The selected cells drive the annual daily-series and cell-information CSV outputs described on the [outputs page](outputs.md).

(crop-input-files)=
## Crop and phenology inputs

The folder selected by {ref}`PhenoPath <parameter-phenopath>` contains `CanopyRes.dat` plus one subfolder per station. Each station folder name is {ref}`PhenoFileRoot <parameter-phenofileroot>` followed by the weather filename without its extension. Include these case-sensitive filenames:

| File | Daily or static content used by IdrAgra |
|---|---|
| `CropId.dat` | crop slot active on each day |
| `Kcb.dat` | basal crop coefficient |
| `H.dat` | crop height |
| `Sr.dat` | root depth |
| `LAI.dat` | leaf-area index |
| `CNvalue.dat` | seasonal Curve Number state |
| `fc.dat` | vegetation cover fraction |
| `r_stress.dat` | stress-related crop parameter |
| `WPadj.dat` | yearly adjusted water productivity |
| `CropParam.dat` | crop classes and static yield, stress, interception, and root parameters |

`CanopyRes.dat` supplies a yearly canopy resistance. These files are normally generated by CropCoef rather than edited manually. The first header in the daily files determines the crop IDs and crop-cycle layout; those IDs must be compatible with {ref}`SoilUsesNum <parameter-soilusesnum>` and {ref}`SimulatedSoilUses <parameter-simulatedsoiluses>`.

:::{container} manual-code-divergence
**Manual/code divergence to review.** The installation manual lists `CropId.dat` as a debug output, but IdrAgra requires it for every station. The manual also lists `Ky.dat`; IdrAgra instead reads `kyT` and `ky1` through `ky4` from `CropParam.dat`.
:::

(soil-moisture-inputs)=
## Initial and final soil moisture

With {ref}`InitialThetaFlag = T <parameter-initialtheta>`, IdrAgra reads two files made from the {ref}`InitialCondition <parameter-initialcondition>` root in {ref}`InitialConditionPath <parameter-initialconditionpath>`:

```text
<InitialCondition>I.asc
<InitialCondition>II.asc
```

For example, `InitialCondition = IC_theta` selects `IC_thetaI.asc` and `IC_thetaII.asc`; omitting the setting uses those same default filenames. Missing values inside the active domain are filled from the corresponding field-capacity map.

With {ref}`InitialThetaFlag = F <parameter-initialtheta>`, IdrAgra starts from field capacity and runs its warm-up year. With {ref}`FinalThetaFlag = T <parameter-finaltheta>`, it writes both layers using {ref}`FinalCondition <parameter-finalcondition>` in {ref}`FinalConditionPath <parameter-finalconditionpath>`.

(cn-class)=
### Crop parameters

...
