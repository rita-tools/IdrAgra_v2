# Model overview

IdrAgra is a spatially-distributed conceptual agro-hydrological model...


(simulation-unit)=
## Simulation unit

IdrAgra's domain is comprised of cells (sometimes referred to as "fields"), which are portions of land considered uniform in terms of weather, soil, landuse and management.\ 

It is at the cell level that the model performs most operations, such as solving the soil-crop-atmosphere balance, applying irrigation water, and calculating yield. 
Each cell is identified by its position in the simulation's domain (row, column)

:::{container} llm-review-note
**LLM-authored draft — review required.** The sections below summarize the 2025 technical manual and were checked against the current Fortran implementation. They should be reviewed by a model maintainer before being treated as authoritative.
:::

## What the model represents

IdrAgra is a spatially distributed agro-hydrological model. It combines four closely connected parts:

1. **Crop phenology**, which supplies daily crop development and canopy/root parameters.
2. **Soil-crop water balance**, which accounts for precipitation, irrigation, interception, runoff, evaporation, transpiration, percolation, ponding, and optionally capillary rise.
3. **Irrigation**, which determines field applications from crop demand, available water, or a supplied calendar, depending on the selected mode.
4. **Crop yield**, which estimates potential biomass and the effects of water and heat stress on actual yield.

The simulation advances one calendar day at a time. Within each day, the current code solves the two-layer soil-water balance in 24 hourly steps, distributes daily meteorological and irrigation quantities across those steps, and accumulates the requested daily, periodic, and annual results.

## Spatial and temporal scale

The surface domain is represented by an ESRI ASCII grid. A valid grid cell is the smallest independent computational unit: it has one soil profile, land use, irrigation method, and set of interpolated meteorological and phenological inputs. Neighboring cells do not exchange soil water laterally in the current implementation.

The optional `shapearea.asc` map lets an exported polygon-based project retain the actual area represented by each computational point. When that file is absent, water-volume calculations use the square of the ASCII-grid cell size.

## Simulation modes

The technical manual groups irrigation runs into the conceptual **USE** and **NEED** families. IdrAgra provides five choices through {ref}`Mode <parameter-mode>`:

| Mode | Current executable behavior | Practical meaning |
|---:|---|---|
| 0 | No irrigation | Run crop and water-balance calculations without irrigation applications. |
| 1 | [USE](use_mode.md) | Route monitored and estimated source water through irrigation units, then distribute the available volume among eligible cells. Demand may remain unsatisfied. |
| 2 | [NEED, field-capacity target](need_modes.md) | Calculate the application needed to replenish the profile toward field capacity, adjusted by {ref}`fc_ratio <parameter-fc-ratio>`. |
| 3 | [NEED, fixed volume](need_modes.md) | Trigger irrigation from soil-water status and apply the irrigation method's fixed depth. |
| 4 | [Scheduled](scheduled_mode.md) | Read dated irrigation depths from the file selected by {ref}`sched_irr_fn <parameter-sched-irr-fn>` and apply them to irrigation units. |

Modes 2 and 3 calculate unconstrained field requirements and therefore belong to the manual's broader NEED concept. Mode 1 represents the manual's USE concept, where source availability, conveyance, and delivery order can limit irrigation.

:::{container} manual-code-divergence
**Manual/code divergence to review.** The technical manual mainly presents USE and NEED, while the model also supports no-irrigation Mode 0 and scheduled Mode 4. Use the mode numbering shown above when preparing `idragra_parameters.txt`.
:::

## Main data flow

For each simulated day, IdrAgra:

1. reads the next daily meteorological record and daily crop-parameter records;
2. updates yearly land-use and irrigation-method maps when configured, and updates the water-table map if a matching dated file exists;
3. interpolates meteorological and phenological information to cells;
4. determines irrigation according to the selected mode and irrigation season;
5. calculates interception and Curve Number runoff;
6. solves the evaporative- and transpirative-layer balances over 24 hourly substeps;
7. updates ponding, crop stress, biomass, and yield state; and
8. writes enabled cell, periodic, annual, and diagnostic outputs.

Meteorological variables can be inverse-distance weighted from several stations or taken from the nearest station, independently for temperature, precipitation, humidity, wind, and radiation. The spatial weight rasters are inputs generated outside the core executable; IdrAgra reads rather than calculates them.
