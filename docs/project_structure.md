# Project structure

In this section we look at how an IdrAgra working folder is structured.



## Input files

Input files make up the majority of the project's data. By default, they are organized in descriptive folders:

- **{ref}`spatial_data <parameter-inputpath>`**, cotaining .asc maps, i.e. anything that is spatialized at the {ref}`field <simulation-unit>` level;
- **{ref}`meteo_data <parameter-meteopath>`**, containing weather time series (one file per station);
- **{ref}`crop_series <parameter-phenopath>`**, containing crop parameters & time series, spatialized at the station level;
- **{ref}`irrmeth_data <parameter-irrmethpath>`**, containing each irrigation method's parameters;
- **{ref}`watsour_data <parameter-watsourpath>`**, containing information about irrigation water sources (used in {ref}`simulation mode <parameter-mode>` 1 only)

Note: here we are using these folders' default names, but note that they may be modified in [idragra_parameters.txt](parameters.md).

...

(water-table-input-files)=
### Water table depth

If the user sets {ref}`CapillaryFlag <parameter-capflag>` to true, the model requires one or more .asc files describing the depth of the water table from the soil surface [m] to be placed in the {ref}`spatial_data <parameter-inputpath>` folder, following a strict naming convention:

- `waterdepth.asc` is read at the beginning of the simulation, and is always required.
- `waterdepth._yyyy_d.asc`, where yyyy and d are the simulation year and doy, are optional and used for subsequent updates.

IdrAgra treats the received water table depths as static until they get updated by reading a new daily map. \
It is not necessary to provide a map for each day of the simulation - IdrAgra will simply continue the simulation using the last map it has read.

(cn-class)=
### Crop parameters

...