# Model overview

IdrAgra is a spatially-distributed conceptual agro-hydrological model...


(simulation-unit)=
## Simulation unit

IdrAgra's domain is comprised of cells (sometimes referred to as "fields"), which are portions of land considered uniform in terms of weather, soil, landuse and management.\ 
It is at the cell level that the model performs most operations, such as solving the soil-crop-atmosphere balance, applying irrigation water, and calculating yield. 
Each cell is identified by its position in the simulation's domain (row, column)