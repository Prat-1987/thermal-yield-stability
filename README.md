# Thermal and Yield Stability Mapping Pipeline

This repository contains reproducible processing code for three outputs used to
study persistent canopy thermal patterns and multi-year yield stability in
rainfed cropping systems:

1. one thermal map for each eligible field-year;
2. a long-term thermal-stability map; and
3. a matching multi-year yield-stability map.

The repository contains **code only**. The field-level thermal imagery and
yield-monitor rasters used in the study are proprietary and are not included.

## Repository contents

- `thermal_yield_stability_pipeline.py` — complete, self-contained analysis.
- `requirements.txt` — Python package requirements.
- `.gitignore` — prevents local rasters and generated results from being
  committed accidentally.

## Analysis workflow

### 1. Yearly thermal map

For each field-year:

1. Thermal GeoTIFFs are selected from the crop-specific growing-season window.
2. Each thermal raster is reprojected to that year's yield grid using bilinear
   resampling.
3. Thermal pixels are restricted to the positive-yield footprint.
4. Each thermal date is independently z-score normalized within the field:

   `z = (pixel temperature - field mean temperature) / field standard deviation`

5. For every pixel, the mean and population standard deviation (`ddof=0`) are
   calculated across the normalized thermal dates in that year.
6. Pixels are classified as:

   - Unstable when thermal standard deviation is greater than 1.3;
   - Cold when stable and in the bottom 20% of mean thermal z-score;
   - Medium when stable and in the middle 60%; or
   - Hot when stable and in the top 20%.

A field-year must contain at least three thermal images with valid overlap.

Yearly thermal-map codes are:

| Code | Class |
|---:|---|
| 0 | Unstable |
| 1 | Cold |
| 2 | Medium |
| 3 | Hot |
| 255 | NoData |

### 2. Long-term thermal-stability map

At least three eligible yearly thermal datasets are required. All individual
thermal dates from those years are independently z-score normalized within
their own field-year, aligned to the earliest eligible yield grid, and combined
into one stack.

The mean and population standard deviation are calculated for every pixel
across the full multi-year date stack. Classification then uses the same rules
as the yearly thermal maps:

- standard deviation > 1.3 = Unstable;
- stable bottom 20% = Cold;
- stable middle 60% = Medium; and
- stable top 20% = Hot.

As in the original analysis, every thermal date is one layer in the long-term
stack. Therefore, a year containing more thermal dates contributes more layers.

### 3. Multi-year yield-stability map

Only yield years that passed the thermal eligibility requirements are used.
Each yield raster is independently z-score normalized within its field-year,
aligned to the same reference grid, and stacked across years.

- standard deviation > 0.9 = Unstable;
- stable bottom 20% of mean yield z-score = Low;
- stable middle 60% = Medium; and
- stable top 20% = High.

To reproduce the comparison-friendly coding in the original analysis, High
yield uses code 1 and Low yield uses code 3. This aligns the expected
Cold/High and Hot/Low relationships:

| Code | Class |
|---:|---|
| 0 | Unstable |
| 1 | High yield |
| 2 | Medium yield |
| 3 | Low yield |
| 255 | NoData |

## Expected input structure

Create one directory per field. Dates must appear in thermal filenames as
`YYYY-MM-DD` or `YYYY_MM_DD`. Years must appear in yield filenames. Include
`corn`/`maize` or `soy` in each yield filename so the program can select the
correct growing-season window.

```text
INSERT_YOUR_INPUT_ROOT_PATH/
├── Field_001/
│   ├── THERMAL_Band1/
│   │   ├── thermal_2020-06-25.tif
│   │   ├── thermal_2020-07-10.tif
│   │   └── thermal_2020-08-05.tif
│   └── Yield/
│       └── corn_yield_2020.tif
└── Field_002/
    ├── THERMAL_Band1/
    └── Yield/
```

Accepted yield-directory names are `Yield`, `yield`, and `Yeild`. The thermal
directory name can be changed using `THERMAL_DIRECTORY_NAME` near the top of
the script.

## Growing-season windows

| Crop | Inclusive thermal window |
|---|---|
| Corn/maize | June 20–August 31 |
| Soybean | July 15–September 30 |

If a crop cannot be inferred from a yield filename, the script issues a warning
and applies the corn window, matching the original analysis.

## Installation

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On Windows, activate the environment with:

```powershell
.venv\Scripts\activate
```

## Running the pipeline

Either replace the path placeholders near the top of the Python file:

```python
INPUT_ROOT = Path("INSERT_YOUR_INPUT_ROOT_PATH")
OUTPUT_ROOT = Path("INSERT_YOUR_OUTPUT_ROOT_PATH")
```

or supply the paths at runtime:

```bash
python thermal_yield_stability_pipeline.py \
  --input-root "INSERT_YOUR_INPUT_ROOT_PATH" \
  --output-root "INSERT_YOUR_OUTPUT_ROOT_PATH"
```

Paths containing spaces should remain inside quotation marks.

## Outputs

```text
INSERT_YOUR_OUTPUT_ROOT_PATH/
├── Field_001/
│   ├── yearly_thermal/
│   │   └── 2020/
│   │       ├── yearly_thermal_mean_zscore.tif
│   │       ├── yearly_thermal_std_zscore.tif
│   │       ├── yearly_thermal_stability.tif
│   │       ├── yearly_thermal_stability_preview.png
│   │       └── thermal_files_used.csv
│   └── stability/
│       ├── all_crops/
│       ├── corn_only/
│       └── soybean_only/
│           ├── thermal_mean_zscore.tif
│           ├── thermal_std_zscore.tif
│           ├── thermal_stability.tif
│           ├── thermal_stability_preview.png
│           ├── thermal_files_used.csv
│           ├── yield_mean_zscore.tif
│           ├── yield_std_zscore.tif
│           ├── yield_stability.tif
│           ├── yield_stability_preview.png
│           └── yield_files_used.csv
├── yearly_thermal_summary_all_fields.csv
└── stability_summary_all_fields.csv
```

The `all_crops` group combines every eligible crop-year. Crop-specific
stability folders are created only when that crop independently has at least
three eligible years.

## Configurable parameters

The methodological constants are grouped near the top of the Python script:

- crop date windows;
- minimum thermal images per field-year;
- minimum field-years per stability analysis;
- stable-pixel percentile split;
- thermal instability threshold;
- yield instability threshold; and
- requested all-crop and crop-specific output groups.

Record any changes to these constants when reporting a sensitivity analysis.

## Reproducibility and input checks

- Thermal and yield inputs must be georeferenced GeoTIFFs with valid coordinate
  reference systems.
- Yield values less than or equal to zero are treated as NoData, following the
  original workflow.
- If multiple yield rasters are found for one field-year, the alphabetically
  first file is used and a warning is printed.
- Each yearly folder records the thermal files used for that map.
- Each long-term folder records the thermal dates and yield years used.
- The script does not calculate field area. Pixel totals should not be reported
  as hectares or acres unless converted using raster resolution in an
  appropriate projected coordinate system.

## Data availability statement



> The field-level thermal imagery and yield-monitor data analyzed in this study
> are subject to commercial and producer-privacy restrictions and therefore
> cannot be distributed publicly. The complete analysis code is publicly
> available and can be applied to georeferenced thermal and yield rasters with
> the documented directory structure. Requests concerning access to derived or
> aggregated data may be directed to the corresponding author and are subject
> to the original data-use agreements.



## Citation

