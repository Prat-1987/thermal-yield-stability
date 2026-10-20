# July-August OpenET water-stress index

This folder is the **code-availability release** for the paper's native-OpenET-grid comparison of cold- and hot-dominant thermal-stability areas. It reproduces the July-August water-stress index (WSI) calculations. It does not contain field rasters, weather records, or results. The standalone script uses an input manifest so private directory names and machine-specific paths do not appear in the published code.

The implementation reproduces the analysis used for the paper. The name *PET* in that analysis denotes FAO-56 Penman-Monteith **reference evapotranspiration** (ET0), not a crop-specific potential ET estimate. The release keeps `pet_mm_day` in output column names for comparison with the original results.

## Files

| File | Purpose |
|---|---|
| `calculate_wsi.py` | Command-line analysis and documented equations |
| `example_config.json` | Anonymized three-field manifest; edit paths and crop years for another study |
| `requirements.txt` | Python dependencies |
| `test_calculate_wsi.py` | Data-free checks of fraction resampling and period aggregation |

## Install and run

Python 3.10+ is recommended. Run from this folder:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest test_calculate_wsi.py
python calculate_wsi.py \
  --config example_config.json \
  --input-root "REPLACE_WITH_INPUT_ROOT" \
  --output-root "REPLACE_WITH_OUTPUT_ROOT"
```

The script requires real input data; `example_config.json` is a **template**, not a downloadable demonstration dataset. Replace the `REPLACE_WITH_...` latitude values with numeric field latitudes, and replace the command-line input/output roots before running. The paths in the JSON manifest are illustrative relative paths: edit them to match your private input layout. Spaces in command-line paths must be quoted. The `id` values can be kept anonymous for publication; one `id` and one output folder are used per field. The configured crop years, not all years available in the weather CSV or OpenET directory, determine which years are analyzed. In the paper's analysis, only corn and soybean years were retained; wheat was excluded.

Expected private input layout for one field:

```text
REPLACE_WITH_INPUT_ROOT/
  Field_1/
    weather.csv
    dem.tif
    thermal_stability.tif
    openet/
      any_prefix_ensemble_et_2016_07_01_any_suffix.tif
      any_prefix_ensemble_et_2016_07_02_any_suffix.tif
      ...
```

The example manifest represents the original three-field crop-year selection. Its thermal-stability path is a placeholder: supply the **exact class raster used for each field in the original analysis**. In particular, the original Field 2 analysis used its corn-only long-term thermal-stability raster for the selected years of both crops. Changing that raster changes the class membership and WSI estimates.

## Required inputs and units

| Input | Required form | How used |
|---|---|---|
| OpenET | One single-band daily `ensemble_et_YYYY_MM_DD` GeoTIFF per available July-August date | Pixel ET in mm d-1, on its **native grid** |
| Thermal stability | Single-band georeferenced GeoTIFF; class 1 = cold, class 3 = hot | Fractional class coverage of each OpenET pixel |
| DEM | Single-band raster with elevation in metres | Mean valid elevation for field-level pressure and radiation |
| Weather CSV | `Year,DOY,SRAD,Tmax,Tmin,Rain,WindSp,Rmax,Rmin` | One daily series for the whole field |
| Manifest | JSON field entries with latitude, file paths, and `crop_by_year` | Explicit field selection and reproducible crop labels |

Weather CSV units **must already be**: `Tmax` and `Tmin` in degrees Celsius; `SRAD` in MJ m-2 d-1; `WindSp` in m s-1; `Rmax` and `Rmin` in percent; and `Rain` in mm d-1. The code makes no Kelvin-to-Celsius or W m-2-to-MJ m-2 d-1 conversion. Confirm those conversions upstream if starting from raw gridMET files. The field's latitude is used in the radiation calculation. A single field-level weather series is used for every OpenET pixel in that field. Longitude is not used by this script after the weather series has been selected. Rainfall characterizes the season; it does **not** enter the Penman-Monteith equation.

The source TIFFs used for the paper reported `scale=1`, `offset=0`, and no unit tag. Thus the mm d-1 interpretation is an **input assumption from the source workflow**, not something independently encoded in the GeoTIFF metadata. The script rejects non-unit raster scale/offset metadata rather than silently ignoring it. It removes NoData, non-finite, and negative ET values, but does not cap large positive ET or clamp WSI to 0-1.

The stored daily OpenET TIFFs for all three paper fields reported **EPSG:4326** and a pixel size of **0.00026949458523585647 degrees** in both axes. This is a geographic-grid resolution, so calling it exactly 30 m would overstate its precision. Their NoData tag was `-9999.0` (some missing pixels also read as non-finite), and their scale/offset were 1/0. This script records no assumed metre conversion for those pixel dimensions.

## Analysis definitions

1. **Native-grid classes.** On each field's native OpenET grid, create 2-m binary masks for thermal class 1 (cold) and class 3 (hot). Average-resample each binary mask to the OpenET grid with `rasterio.warp.reproject(..., Resampling.average)`. Source non-class and NoData cells contribute zero to the mask. A native OpenET pixel is cold-dominant when its cold fraction is at least 0.60, or hot-dominant when its hot fraction is at least 0.60. The ET raster is never resampled to the finer thermal/yield grid. If a field's daily OpenET grids are not identical, the script stops with an error.
2. **Daily ET0/PET.** Calculate saturation vapour pressure separately at maximum and minimum air temperature: `es(T) = 0.6108 exp(17.27 T/(T + 237.3))` kPa. Then `es = [es(Tmax) + es(Tmin)]/2` and `ea = [es(Tmin) RHmax/100 + es(Tmax) RHmin/100]/2`; RH is constrained to 0-100%. Vapour-pressure deficit is `max(es-ea, 0)`.
3. **Radiation and FAO-56.** Calculate extraterrestrial radiation using latitude and day of year. Calculate clear-sky radiation using field-mean DEM elevation, net shortwave radiation using albedo 0.23, and net longwave radiation with the FAO-56 Stefan-Boltzmann formulation. Net radiation is shortwave minus longwave. Daily soil heat flux `G` is set to zero. Reference ET is:

   ```text
   ET0 = [0.408 Delta (Rn - G) + gamma (900/(Tmean + 273)) u VPD]
         / [Delta + gamma (1 + 0.34 u)]
   ```

   Here `Delta` is the slope of the saturation-pressure curve and `gamma` is the elevation-dependent psychrometric constant. The code retains the original longwave-radiation details and constants in comments and equations.
4. **Wind-height caveat.** `WindSp` is passed directly as `u` in the FAO-56 equation. The original analysis did **not** convert nominal gridMET 10-m wind speed to the FAO-56 2-m wind convention. This is a documented simplifying assumption, applied consistently across field-years; do not describe it as a height-corrected FAO-56 calculation.
5. **Date matching and completeness.** ET and PET are matched by calendar date. If an OpenET raster is unavailable, its ET and corresponding PET are both omitted from the period sum; neither is imputed. For each field-year, a pixel is retained only when it has a valid ET value on every **sampled** OpenET date in that year. For the original Field 1 data, 2016 and 2019 had 61 included dates, and the other selected years had 62. Missing dates and missing/invalid weather are distinguished in `date_qc.csv`.
6. **Period WSI.** For each retained pixel, calculate `ET_total = sum(daily ET)` and `PET_total = sum(daily ET0)` over the matched July-August dates, followed by `WSI = 1 - ET_total/PET_total`. This is a ratio of period totals, **not** an average of daily `1 - ET/ET0`. The class-year WSI is the unweighted mean of retained native OpenET pixel WSI values. A positive `hot_minus_cold_wsi` means the hot-dominant region had higher calculated WSI.

The July-August rainfall total in output tables uses the **full available weather calendar**, including rainfall on a date lacking an OpenET raster. It is contextual, not part of either ET or PET sum. `required_days` instead records the number of sampled ET/PET dates.

## Output files

```text
outputs/
  all_fields_field_year_wsi.csv
  all_fields_hot_cold_wsi_gap.csv
  Field_1/
    daily_pet_weather.csv
    date_qc.csv
    pixel_year_wsi.csv
    field_year_wsi.csv
    hot_cold_wsi_gap.csv
  Field_2/
    ...
  Field_3/
    ...
```

`daily_pet_weather.csv` contains the selected growing-season weather dates and computed daily reference ET, including a date for which OpenET is absent. `date_qc.csv` gives one row per July-August calendar date and identifies missing OpenET, missing weather, invalid PET, no valid ET, or sampled dates. `pixel_year_wsi.csv` contains native-grid row/column, class fraction, ET and PET totals, WSI, and ET deficit for retained pixels. `field_year_wsi.csv` summarizes each class and crop-year. `hot_cold_wsi_gap.csv` gives the within-field-year difference in class-mean WSI and July-August rainfall. No statistical test or causal inference is implied by the gap.

## Reproducibility and limitations

- The original code was cross-checked against this standalone release for all 18 selected field-years (six per field): class-mean WSI differed by at most approximately `6e-9`. The two 61-day Field 1 years retained their original date counts. The tests cover fractional aggregation, complete-date period ratios, and a small end-to-end case without private data.
- The method compares class-dominant OpenET pixels at native spatial resolution. It does not interpret a 2-m thermal class pixel as an independent ET observation. Class means may share spatially correlated pixels and are descriptive.
- OpenET product selection and upstream extraction, original gridMET unit conversions, and creation of the thermal-stability raster are **outside** this WSI script. Record their provenance separately in the manuscript or data-availability materials.
- The output root may contain spatially resolved rows/columns; keep it private if the underlying field locations cannot be distributed.

## Methods references

- Allen, R. G., Pereira, L. S., Raes, D., and Smith, M. (1998). *Crop evapotranspiration: Guidelines for computing crop water requirements*. FAO Irrigation and Drainage Paper 56.
- Abatzoglou, J. T. (2013). Development of gridded surface meteorological data for ecological applications and modelling. *International Journal of Climatology*, 33, 121-131. Check the manuscript's reference style and source-data citation for the final paper.

