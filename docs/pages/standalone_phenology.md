# Standalone crop development: first implementation

IdrAgra reads raw crop definitions and advances development inside its daily
water-balance loop. CropCoef is no longer invoked, and generated station crop
tables are not read. This is a development implementation, with the explicit
modelling choices and limitations below; it is not numerically equivalent to
the old preprocessing workflow.

## Run and migrate

Build with `make` (or `make debug=1` after a clean build) and run the executable
from the `demo` directory. Raw demonstration inputs are in `demo/landuses`.
The old `demo/pheno` files, if present, are historical fixtures and are unused.

Replace `PhenoPath` and `PhenoFileRoot` in the simulation configuration with:

```text
CropRotationFile = landuses/soil_uses.txt
CropParameterPath = landuses/crop_parameters
CropTemperatureWindow = 2
CropCO2 = 0
```

Paths are relative to the process working directory. Crop input paths accept
forward slashes and may be quoted. The rotation file uses the existing format:

```text
Cr_ID Crop1 Crop2
1 winter_wheat.tab maize.tab
2 maize.tab *
3 * *
endTable
```

Define every ID from 1 to `SoilUsesNum`. An empty list represents bare soil;
zero-GDD entries are skipped. More than two rotation slots are supported within
the input line limit. Repeated crop filenames remain distinct slots.

Crop files use CropCoef's numeric parameter names and increasing GDD tables:
`GDD Kcb LAI Hc Sr [CN [fc [r_stress [Ky]]]]`. Missing entries use `*`.
Required scheduling fields are `SowingDate_min`, `SowingDelay_max`,
`HarvestDate_max`, `Tdaybase`, and `Tcutoff`. Missing curve knots are interpolated
in thermal time; completely absent required growth curves and malformed tables
are rejected. Missing cover fraction retains IdrAgra's computed-cover sentinel.

`CropCO2 = 0` keeps raw water productivity and uses canopy resistance 70 s/m.
A positive constant concentration activates CropCoef's water-productivity and
canopy-resistance formulas. **Annual CO2 time series are not implemented yet.**
No `CanopyRes.dat`, `WPadj.dat`, or `CropParam.dat` is required.

## Daily execution and storage

1. Prepare today's spatialized weather using the existing interpolation flags.
2. Advance persistent crop arrays: rotation slot, occurrence, cut, sowing date,
   GDD and vernalization. Generate today's crop-property matrices.
3. Update the root zone and solve the existing water balance and irrigation.
4. Accumulate crop production and finalize harvest/cut events after today's
   transpiration has been included.

State uses spatial arrays allocated once, with the first index in the inner
loop. Crop parameter tables are shared by rotation slot. Full temperature,
humidity and wind records are cached at stations, not for every cell. The
first implementation has not been benchmarked on large domains.

Thermal time uses the sine-wave GDD calculation and CropCoef's existing
vernalization/photoperiod factors. Properties are interpolated directly in GDD;
CN is a step function. Sowing has zero accumulated development on its first day.
`CropTemperatureWindow` specifies the number of days on either side of today
used to smooth mean temperature for sowing and vernalization. Raw daily Tmax
and Tmin still drive the GDD calculation.

## Rotations and variability

- The first crop is selected by the nearest remaining sowing deadline. Without
  crop warm-up, the simulation begins bare; an autumn crop is not reconstructed
  from soil-moisture initial-condition files.
- Subsequent crops follow the declared slot order. Suitable temperatures allow
  sowing within the window; unmet temperatures force sowing at its deadline.
- The preceding crop is terminated by the next deadline minus the required
  bare-soil interval. `CropsOverlap` is interpreted as a minimum interval, not
  simultaneous cultivation. Impossible schedules may establish a crop late;
  the event is reported as forced sowing.
- Positive sowing offsets delay eligibility; they do not translate an existing
  crop's growth curve or harvest date.
- Offset input precedence is annual `irandom_YEAR.asc`, static `irandom.asc`,
  then deterministic generation. Values must be in -365..365 on active cells.
  Generated offsets use `RandSowDaysWind`, `RandSowDaysSym`, `Repeatable`, and
  `RandSowDaysSeed`. Sampling does not depend on loop order. An offset is captured
  when the next occurrence is scheduled; later annual maps do not resample it.
- Annual land-use changes currently terminate the active crop and select a new
  sequence. This preserves the annual map replacement concept; redesigning
  multi-year land-use semantics remains separate work.
- Multiple cuts reset development and production accumulators while retaining
  the crop occurrence, cut count, and vernalization. Dedicated perennial
  dormancy, persistent woody roots and establishment-age models are not present.

## Kcb climate correction

The engine projects development from sowing using available weather to identify
mid-season and late-season intervals. It averages wind, daily minimum relative
humidity and crop height over those intervals and applies the FAO-56 Equation
70 correction to eligible coefficients. Corrections are interpolated into the
daily curve; the initial coefficient is unchanged. `adj_flag = 0` disables
correction for locally calibrated curves.

Stage boundaries are inferred once from the uncorrected Kcb knots. This needs
review for unusual curves, especially curves with no mid-season plateau.
A projected stage that is never reached receives no correction. Development
projection is repeated after a cut. Boundary lookahead uses the corresponding
calendar date in the nearest available weather year, clipped to the available
record when necessary. This is an explicit approximation, not observed future
weather. Forced early termination and atypical/perennial curves warrant further
scientific validation.

References: [FAO-56, Chapter 7, Equation 70](https://www.fao.org/4/X0490E/x0490e0c.htm)
and [Chapter 6, stage climate averages](https://www.fao.org/4/X0490E/x0490e0b.htm).

## Production and diagnostics

`crop_events.csv` records cell indices, land use, rotation slot, occurrence,
cut, GDD and sowing offset at sowing, harvest, forced termination, domain/map
changes and the end of a run. Living crops are reported as ongoing; ending a
simulation does not harvest them.

`crop_harvests.csv` contains one row per harvest/cut, with potential biomass,
potential yield, whole-cycle and stage water-stress factors, heat-stress factor,
and actual yield. Accumulators persist across calendar years. The formulas
follow v2's potential-transpiration/WP and HI approach, with zero-denominator
guards. First-pass heat sensitivity is at **45–75% of thermal development**,
replacing the old fraction of an already known calendar duration.

Early termination uses the supplied HI and reports its reason. It does not
infer silage yield or a maturity-dependent harvest index. Land-use/domain
termination records contain production accumulated before the transition day.

These CSVs replace annual crop-slot yield rasters and cell-parameter summaries.
They are currently always written; legacy `prt_yield`/individual annual yield
switches do not control them. Existing hydrological maps and daily control-cell
outputs remain available. Exporters and downstream yield consumers need updates.

## Initialization

The existing soil warm-up, selected with `InitialThetaFlag = F`, also advances
crops. Their biological state and unfinished production are retained when the
main run starts. Event dates are rebased across the weather rewind without
resetting development. Warm-up logs have a `warmup_` prefix; harvests containing
warm-up production are marked `contains_warmup`.

This is a minimal state handoff, not a redesign of warm-up. General mid-season
crop initialization, crop restart files and integrating warm-up into a single
regular simulation loop remain deferred.

## Verification

Run `python tests/run_standalone.py` with gfortran on PATH. The harness builds
into `.codex_tmp`, enables bounds checks and floating-point traps, and creates
isolated demo copies **without any generated phenology directory**. It checks
thermal equations, interpolation, calendar windows, a two-year USE simulation,
repeatable events, NEED-mode warm-up/cross-year crops, forced sowing without
irrigation, finite outputs, and rejection of malformed crop curves.

The initial integration is intended for review and experiments. Passing these
checks does not establish agronomic calibration or equivalence with CropCoef.
