#!/usr/bin/env python3
"""Create yearly thermal maps and multi-year thermal/yield stability maps.

This publication version consolidates the relevant yearly and all-year
stability processing steps into one self-contained program. It contains no
study data and no private filesystem paths.

Expected input layout
---------------------

INPUT_ROOT/
    FIELD_NAME/
        THERMAL_Band1/
            ...YYYY-MM-DD....tif
            ...YYYY_MM_DD....tif
        Yield/                 # ``yield`` is also accepted
            ...YYYY...corn....tif
            ...YYYY...soy....tif

Important assumptions retained from the original analysis
---------------------------------------------------------

* Corn thermal window: June 20 through August 31.
* Soybean thermal window: July 15 through September 30.
* Thermal values are reprojected to the corresponding yield raster grid with
  bilinear resampling and clipped to the positive-yield footprint.
* Every thermal-date raster is z-score normalized within its field-year.
* Each yearly thermal map is calculated from at least three normalized dates.
* A year is eligible when it has at least three valid thermal dates.
* A stability map is eligible when it has at least three eligible years.
* Long-term thermal instability: standard deviation > 1.3.
* Long-term yield instability: standard deviation > 0.9.
* Stable thermal and yield pixels use a 20/60/20 percentile split.

The thresholds are constants near the top of this file so sensitivity analyses
can be run without modifying the processing functions.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import Resampling, reproject

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    MATPLOTLIB_AVAILABLE = True
except ImportError:  # GeoTIFF creation can still run without optional previews.
    MATPLOTLIB_AVAILABLE = False


# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Replace these placeholders, or supply --input-root and --output-root when
# running the script.  No proprietary study paths are included in this file.
INPUT_ROOT = Path("INSERT_YOUR_INPUT_ROOT_PATH")
OUTPUT_ROOT = Path("INSERT_YOUR_OUTPUT_ROOT_PATH")

THERMAL_DIRECTORY_NAME = "THERMAL_Band1"
YIELD_DIRECTORY_CANDIDATES = ("Yield", "yield")

# Crop-specific growing-season windows, written as MM-DD.
CORN_WINDOW = ("06-20", "08-31")
SOYBEAN_WINDOW = ("07-15", "09-30")

# Eligibility rules used for the publication workflow.
MIN_THERMAL_DATES_PER_YEAR = 3
MIN_YEARS_FOR_STABILITY = 3

# Percentile and variability thresholds retained from the existing code.
STABILITY_TAIL_PERCENT = 20   # Stable pixels: 20% / 60% / 20%.
THERMAL_STD_UNSTABLE = 1.3
YIELD_STD_UNSTABLE = 0.9

# By default, produce an all-crop result and crop-specific results whenever
# each group independently satisfies the minimum-year rule.
STABILITY_GROUPS = ("all_crops", "corn_only", "soybean_only")

NODATA_CLASS = 255
NODATA_FLOAT = -9999.0

# Class codes intentionally reproduce the comparison-friendly coding in the
# original all-year script.  Thermal code 1 (Cold) aligns with yield code 1
# (High), while thermal code 3 (Hot) aligns with yield code 3 (Low).
THERMAL_CLASS_NAMES = {0: "Unstable", 1: "Cold", 2: "Medium", 3: "Hot"}
YIELD_CLASS_NAMES = {0: "Unstable", 1: "High", 2: "Medium", 3: "Low"}

COLOR_UNSTABLE = "#FFFF00"
COLOR_COLD = "#0070FF"
COLOR_MEDIUM = "#70A800"
COLOR_HOT = "#A80000"

THERMAL_COLORS = [COLOR_UNSTABLE, COLOR_COLD, COLOR_MEDIUM, COLOR_HOT]
YIELD_COLORS = [COLOR_UNSTABLE, COLOR_COLD, COLOR_MEDIUM, COLOR_HOT]

DATE_PATTERN = re.compile(r"(20\d{2}|19\d{2})[-_](\d{2})[-_](\d{2})")
YEAR_PATTERN = re.compile(r"(20\d{2}|19\d{2})")
RASTER_SUFFIXES = {".tif", ".tiff"}


@dataclass(frozen=True)
class FieldYearRecord:
    """Files used for one field-year."""

    field_name: str
    year: int
    crop: str
    yield_path: Path
    thermal_paths: tuple[Path, ...]


@dataclass(frozen=True)
class RasterGrid:
    """Grid metadata needed for reprojection and GeoTIFF writing."""

    profile: dict
    transform: object
    crs: object
    height: int
    width: int


# =============================================================================
# FILE DISCOVERY AND BASIC RASTER OPERATIONS
# =============================================================================


def raster_files(directory: Path) -> list[Path]:
    """Return all GeoTIFF files immediately inside *directory*."""

    if not directory.is_dir():
        return []
    return sorted(
        path
        for path in directory.iterdir()
        if path.is_file() and path.suffix.lower() in RASTER_SUFFIXES
    )


def extract_year(path: Path) -> int | None:
    """Extract the first four-digit year from a filename."""

    match = YEAR_PATTERN.search(path.name)
    return int(match.group(1)) if match else None


def extract_date(path: Path) -> pd.Timestamp | None:
    """Extract YYYY-MM-DD or YYYY_MM_DD from a filename."""

    match = DATE_PATTERN.search(path.name)
    if not match:
        return None
    try:
        return pd.Timestamp(
            year=int(match.group(1)),
            month=int(match.group(2)),
            day=int(match.group(3)),
        )
    except ValueError:
        return None


def infer_crop(yield_path: Path) -> str:
    """Infer crop from the yield filename, reproducing the original rule."""

    name = yield_path.name.lower()
    if "soy" in name:
        return "soybean"
    if "corn" in name or "maize" in name:
        return "corn"

    print(
        f"WARNING: Crop not found in {yield_path.name!r}; assuming corn.",
        file=sys.stderr,
    )
    return "corn"


def crop_window(year: int, crop: str) -> tuple[pd.Timestamp, pd.Timestamp]:
    """Return the inclusive thermal-image date window for a crop-year."""

    start_mmdd, end_mmdd = SOYBEAN_WINDOW if crop == "soybean" else CORN_WINDOW
    start_month, start_day = map(int, start_mmdd.split("-"))
    end_month, end_day = map(int, end_mmdd.split("-"))
    return (
        pd.Timestamp(year=year, month=start_month, day=start_day),
        pd.Timestamp(year=year, month=end_month, day=end_day),
    )


def select_thermal_files(field_dir: Path, year: int, crop: str) -> tuple[Path, ...]:
    """Select dated thermal rasters inside the crop-specific date window."""

    start, end = crop_window(year, crop)
    selected: list[tuple[pd.Timestamp, Path]] = []
    for path in raster_files(field_dir / THERMAL_DIRECTORY_NAME):
        date = extract_date(path)
        if date is not None and date.year == year and start <= date <= end:
            selected.append((date, path))
    return tuple(path for _, path in sorted(selected))


def discover_records(field_dir: Path) -> list[FieldYearRecord]:
    """Discover one yield raster and matching thermal dates for each year."""

    yield_paths: list[Path] = []
    for directory_name in YIELD_DIRECTORY_CANDIDATES:
        yield_paths.extend(raster_files(field_dir / directory_name))

    # On case-insensitive filesystems, ``Yield`` and ``yield`` can point to the
    # same directory. Deduplicate by the underlying file identity rather than
    # by path spelling so one raster is never counted twice.
    unique_yield_paths: dict[tuple[int, int], Path] = {}
    for path in yield_paths:
        stat = path.stat()
        unique_yield_paths.setdefault((stat.st_dev, stat.st_ino), path)

    by_year: dict[int, list[Path]] = {}
    for path in sorted(unique_yield_paths.values()):
        year = extract_year(path)
        if year is not None:
            by_year.setdefault(year, []).append(path)

    records: list[FieldYearRecord] = []
    for year, candidates in sorted(by_year.items()):
        # This reproduces the prior scripts' deterministic first-file behavior.
        # A warning makes duplicate choices transparent for publication runs.
        if len(candidates) > 1:
            print(
                f"WARNING: {field_dir.name} {year} has {len(candidates)} yield "
                f"rasters; using {candidates[0].name!r}.",
                file=sys.stderr,
            )
        yield_path = candidates[0]
        crop = infer_crop(yield_path)
        thermal_paths = select_thermal_files(field_dir, year, crop)

        if len(thermal_paths) < MIN_THERMAL_DATES_PER_YEAR:
            print(
                f"SKIP {field_dir.name} {year}: found {len(thermal_paths)} "
                f"thermal date(s); need at least {MIN_THERMAL_DATES_PER_YEAR}."
            )
            continue

        records.append(
            FieldYearRecord(
                field_name=field_dir.name,
                year=year,
                crop=crop,
                yield_path=yield_path,
                thermal_paths=thermal_paths,
            )
        )
    return records


def read_yield(path: Path) -> tuple[np.ndarray, RasterGrid]:
    """Read yield and mask NoData, non-finite, zero, and negative values."""

    with rasterio.open(path) as source:
        array = source.read(1).astype("float32")
        nodata = source.nodata
        grid = RasterGrid(
            profile=source.profile.copy(),
            transform=source.transform,
            crs=source.crs,
            height=source.height,
            width=source.width,
        )

    if nodata is not None and np.isfinite(float(nodata)):
        array[array == nodata] = np.nan
    array[~np.isfinite(array) | (array <= 0)] = np.nan
    return array, grid


def same_grid(first: RasterGrid, second: RasterGrid) -> bool:
    """Return True when two raster grids have identical geometry."""

    return (
        first.height == second.height
        and first.width == second.width
        and first.transform == second.transform
        and first.crs == second.crs
    )


def reproject_thermal(path: Path, destination_grid: RasterGrid) -> np.ndarray:
    """Reproject one thermal raster to a yield grid using bilinear resampling."""

    with rasterio.open(path) as source:
        source_array = source.read(1).astype("float32")
        if source.nodata is not None and np.isfinite(float(source.nodata)):
            source_array[source_array == source.nodata] = np.nan

        destination = np.full(
            (destination_grid.height, destination_grid.width),
            np.nan,
            dtype="float32",
        )
        reproject(
            source=source_array,
            destination=destination,
            src_transform=source.transform,
            src_crs=source.crs,
            dst_transform=destination_grid.transform,
            dst_crs=destination_grid.crs,
            src_nodata=np.nan,
            dst_nodata=np.nan,
            resampling=Resampling.bilinear,
        )
    return destination


def reproject_array(
    array: np.ndarray,
    source_grid: RasterGrid,
    destination_grid: RasterGrid,
) -> np.ndarray:
    """Reproject an in-memory normalized array to a reference grid."""

    destination = np.full(
        (destination_grid.height, destination_grid.width),
        np.nan,
        dtype="float32",
    )
    reproject(
        source=array.astype("float32"),
        destination=destination,
        src_transform=source_grid.transform,
        src_crs=source_grid.crs,
        dst_transform=destination_grid.transform,
        dst_crs=destination_grid.crs,
        src_nodata=np.nan,
        dst_nodata=np.nan,
        resampling=Resampling.bilinear,
    )
    return destination


def zscore_within_field(array: np.ndarray, field_mask: np.ndarray) -> np.ndarray:
    """Calculate z = (x - field mean) / field standard deviation.

    NumPy's population standard deviation is used (``ddof=0``), matching the
    original scripts. Values outside the field mask remain NaN.
    """

    valid = np.isfinite(array) & field_mask
    normalized = np.full(array.shape, np.nan, dtype="float32")
    if not np.any(valid):
        return normalized

    mean = float(np.nanmean(array[valid]))
    standard_deviation = float(np.nanstd(array[valid], ddof=0))
    if not np.isfinite(standard_deviation) or standard_deviation == 0:
        return normalized

    normalized[valid] = (
        (array[valid] - mean) / standard_deviation
    ).astype("float32")
    return normalized


def write_float_tif(path: Path, array: np.ndarray, grid: RasterGrid) -> None:
    """Write a continuous float32 GeoTIFF."""

    profile = grid.profile.copy()
    profile.update(dtype="float32", count=1, nodata=NODATA_FLOAT, compress="lzw")
    output = np.where(np.isfinite(array), array, NODATA_FLOAT).astype("float32")
    with rasterio.open(path, "w", **profile) as destination:
        destination.write(output, 1)


def write_class_tif(path: Path, array: np.ndarray, grid: RasterGrid) -> None:
    """Write a categorical uint8 GeoTIFF."""

    profile = grid.profile.copy()
    profile.update(dtype="uint8", count=1, nodata=NODATA_CLASS, compress="lzw")
    with rasterio.open(path, "w", **profile) as destination:
        destination.write(array.astype("uint8"), 1)


def write_rows(path: Path, rows: Sequence[dict]) -> None:
    """Write dictionaries to CSV, including an empty file when needed."""

    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


# =============================================================================
# CLASSIFICATION AND MAP PREVIEWS
# =============================================================================


def classify_stability(
    mean_array: np.ndarray,
    standard_deviation_array: np.ndarray,
    reference_mask: np.ndarray,
    unstable_threshold: float,
    reverse_tail_codes: bool = False,
) -> tuple[np.ndarray, float, float]:
    """Classify unstable pixels, then split stable pixels 20/60/20.

    ``reverse_tail_codes`` reproduces the original yield-stability convention:
    high stable yield is code 1 and low stable yield is code 3. This makes the
    expected Cold/High and Hot/Low relationships share class codes.
    """

    valid = (
        np.isfinite(mean_array)
        & np.isfinite(standard_deviation_array)
        & reference_mask
    )
    classes = np.full(mean_array.shape, NODATA_CLASS, dtype="uint8")
    if not np.any(valid):
        return classes, np.nan, np.nan

    unstable = valid & (standard_deviation_array > unstable_threshold)
    stable = valid & ~unstable
    classes[unstable] = 0
    if not np.any(stable):
        return classes, np.nan, np.nan

    low_cut = float(
        np.nanpercentile(mean_array[stable], STABILITY_TAIL_PERCENT)
    )
    high_cut = float(
        np.nanpercentile(mean_array[stable], 100 - STABILITY_TAIL_PERCENT)
    )

    low_code, high_code = (3, 1) if reverse_tail_codes else (1, 3)
    classes[stable & (mean_array <= low_cut)] = low_code
    classes[stable & (mean_array > low_cut) & (mean_array < high_cut)] = 2
    classes[stable & (mean_array >= high_cut)] = high_code
    return classes, low_cut, high_cut


def save_class_preview(
    classes: np.ndarray,
    path: Path,
    title: str,
    names: dict[int, str],
    colors: Sequence[str],
    lower_cut: float,
    upper_cut: float,
    subtitle: str,
) -> None:
    """Save a publication-friendly categorical PNG preview."""

    if not MATPLOTLIB_AVAILABLE:
        print(
            "WARNING: matplotlib is unavailable; PNG previews were not created. "
            "Install the packages in requirements.txt to enable them.",
            file=sys.stderr,
        )
        return

    plot_array = classes.astype("float32")
    plot_array[classes == NODATA_CLASS] = np.nan
    cmap = ListedColormap(list(colors))
    if 0 in names:
        norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    else:
        norm = BoundaryNorm([0.5, 1.5, 2.5, 3.5], cmap.N)

    figure, axis = plt.subplots(figsize=(8, 7), facecolor="white")
    axis.imshow(plot_array, cmap=cmap, norm=norm, interpolation="nearest")
    axis.axis("off")
    cut_text = (
        f"lower cut={lower_cut:.3f}, upper cut={upper_cut:.3f}"
        if np.isfinite(lower_cut) and np.isfinite(upper_cut)
        else "class cuts unavailable"
    )
    axis.set_title(f"{title}\n{subtitle}\n{cut_text}", fontsize=11, weight="bold")
    axis.legend(
        handles=[
            mpatches.Patch(color=color, label=names[code])
            for code, color in zip(sorted(names), colors)
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=len(names),
        frameon=False,
        fontsize=9,
    )
    figure.tight_layout()
    figure.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)


# =============================================================================
# YEARLY THERMAL MAPS
# =============================================================================


def create_yearly_thermal_map(
    record: FieldYearRecord,
    field_output: Path,
) -> dict | None:
    """Create one seasonal thermal-stability map for a field-year.

    Every input date is independently z-scored within the positive-yield
    footprint. The per-pixel mean and population standard deviation are then
    calculated across dates. Pixels with standard deviation greater than 1.3
    are Unstable; the remaining pixels are split 20/60/20 by mean z-score.
    """

    yield_array, yield_grid = read_yield(record.yield_path)
    field_mask = np.isfinite(yield_array)
    output_dir = field_output / "yearly_thermal" / str(record.year)
    output_dir.mkdir(parents=True, exist_ok=True)

    layers: list[np.ndarray] = []
    files_used: list[dict] = []
    for thermal_path in record.thermal_paths:
        thermal = reproject_thermal(thermal_path, yield_grid)
        thermal[~field_mask] = np.nan
        normalized = zscore_within_field(thermal, field_mask)
        if not np.any(np.isfinite(normalized)):
            print(f"  SKIP invalid thermal/yield overlap: {thermal_path.name}")
            continue
        layers.append(normalized.astype("float32"))
        date = extract_date(thermal_path)
        files_used.append(
            {
                "field": record.field_name,
                "year": record.year,
                "crop": record.crop,
                "date": date.strftime("%Y-%m-%d") if date is not None else "",
                "thermal_source": str(thermal_path),
                "yield_mask_source": str(record.yield_path),
            }
        )

    if len(layers) < MIN_THERMAL_DATES_PER_YEAR:
        print(
            f"  EXCLUDE {record.year}: only {len(layers)} thermal date(s) "
            "had valid overlap."
        )
        return None

    mean_thermal, std_thermal = stack_mean_std(np.stack(layers, axis=0))
    classes, cold_cut, hot_cut = classify_stability(
        mean_thermal,
        std_thermal,
        field_mask,
        THERMAL_STD_UNSTABLE,
        reverse_tail_codes=False,
    )

    mean_path = output_dir / "yearly_thermal_mean_zscore.tif"
    std_path = output_dir / "yearly_thermal_std_zscore.tif"
    class_path = output_dir / "yearly_thermal_stability.tif"
    preview_path = output_dir / "yearly_thermal_stability_preview.png"
    files_used_path = output_dir / "thermal_files_used.csv"

    write_float_tif(mean_path, mean_thermal, yield_grid)
    write_float_tif(std_path, std_thermal, yield_grid)
    write_class_tif(class_path, classes, yield_grid)
    save_class_preview(
        classes=classes,
        path=preview_path,
        title=f"{record.field_name} | {record.year} | {record.crop}",
        names=THERMAL_CLASS_NAMES,
        colors=THERMAL_COLORS,
        lower_cut=cold_cut,
        upper_cut=hot_cut,
        subtitle=(
            f"Yearly thermal map from {len(layers)} dates; "
            f"Unstable when std > {THERMAL_STD_UNSTABLE}"
        ),
    )
    write_rows(files_used_path, files_used)

    return {
        "field": record.field_name,
        "year": record.year,
        "crop": record.crop,
        "thermal_dates_used": len(layers),
        "unstable_threshold": THERMAL_STD_UNSTABLE,
        "cold_cut": cold_cut,
        "hot_cut": hot_cut,
        "unstable_pixels": int(np.sum(classes == 0)),
        "cold_pixels": int(np.sum(classes == 1)),
        "medium_pixels": int(np.sum(classes == 2)),
        "hot_pixels": int(np.sum(classes == 3)),
        "mean_zscore_tif": str(mean_path),
        "std_zscore_tif": str(std_path),
        "yearly_thermal_map_tif": str(class_path),
        "preview_png": str(preview_path),
        "files_used_csv": str(files_used_path),
    }


# =============================================================================
# MULTI-YEAR STABILITY MAPS
# =============================================================================


def group_records(
    records: Sequence[FieldYearRecord], group_name: str
) -> list[FieldYearRecord]:
    """Select all-crop or crop-specific records."""

    if group_name == "corn_only":
        return [record for record in records if record.crop == "corn"]
    if group_name == "soybean_only":
        return [record for record in records if record.crop == "soybean"]
    return list(records)


def build_thermal_stack(
    records: Sequence[FieldYearRecord],
    reference_grid: RasterGrid,
    reference_mask: np.ndarray,
) -> tuple[np.ndarray | None, list[dict]]:
    """Build the all-date stack of within-field thermal z-scores."""

    layers: list[np.ndarray] = []
    rows: list[dict] = []
    for record in records:
        yield_array, year_grid = read_yield(record.yield_path)
        year_mask = np.isfinite(yield_array)
        for thermal_path in record.thermal_paths:
            thermal = reproject_thermal(thermal_path, year_grid)
            thermal[~year_mask] = np.nan
            normalized = zscore_within_field(thermal, year_mask)
            if not np.any(np.isfinite(normalized)):
                continue
            if not same_grid(year_grid, reference_grid):
                normalized = reproject_array(normalized, year_grid, reference_grid)
            normalized[~reference_mask] = np.nan
            if not np.any(np.isfinite(normalized)):
                continue

            layers.append(normalized.astype("float32"))
            date = extract_date(thermal_path)
            rows.append(
                {
                    "field": record.field_name,
                    "year": record.year,
                    "crop": record.crop,
                    "date": date.strftime("%Y-%m-%d") if date is not None else "",
                    "thermal_source": str(thermal_path),
                    "yield_mask_source": str(record.yield_path),
                }
            )
    return (np.stack(layers, axis=0) if layers else None), rows


def build_yield_stack(
    records: Sequence[FieldYearRecord],
    reference_grid: RasterGrid,
    reference_mask: np.ndarray,
) -> tuple[np.ndarray | None, list[dict]]:
    """Build the multi-year stack of within-field yield z-scores."""

    layers: list[np.ndarray] = []
    rows: list[dict] = []
    for record in records:
        yield_array, year_grid = read_yield(record.yield_path)
        year_mask = np.isfinite(yield_array)
        normalized = zscore_within_field(yield_array, year_mask)
        if not np.any(np.isfinite(normalized)):
            continue
        if not same_grid(year_grid, reference_grid):
            normalized = reproject_array(normalized, year_grid, reference_grid)
        normalized[~reference_mask] = np.nan
        if not np.any(np.isfinite(normalized)):
            continue

        layers.append(normalized.astype("float32"))
        rows.append(
            {
                "field": record.field_name,
                "year": record.year,
                "crop": record.crop,
                "yield_source": str(record.yield_path),
                "eligible_thermal_dates": len(record.thermal_paths),
            }
        )
    return (np.stack(layers, axis=0) if layers else None), rows


def stack_mean_std(stack: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Calculate per-pixel population mean and standard deviation."""

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean = np.nanmean(stack, axis=0).astype("float32")
        standard_deviation = np.nanstd(stack, axis=0, ddof=0).astype("float32")
    return mean, standard_deviation


def create_stability_maps(
    records: Sequence[FieldYearRecord],
    group_name: str,
    field_output: Path,
) -> list[dict]:
    """Create paired long-term thermal and yield stability products."""

    years = sorted({record.year for record in records})
    if len(years) < MIN_YEARS_FOR_STABILITY:
        print(
            f"  SKIP {group_name}: {len(years)} eligible year(s); "
            f"need at least {MIN_YEARS_FOR_STABILITY}."
        )
        return []

    reference_record = sorted(records, key=lambda record: record.year)[0]
    reference_yield, reference_grid = read_yield(reference_record.yield_path)
    reference_mask = np.isfinite(reference_yield)
    output_dir = field_output / "stability" / group_name
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []

    # ---- Long-term thermal stability ----
    thermal_stack, thermal_files = build_thermal_stack(
        records, reference_grid, reference_mask
    )
    if thermal_stack is not None:
        thermal_mean, thermal_std = stack_mean_std(thermal_stack)
        thermal_classes, cold_cut, hot_cut = classify_stability(
            thermal_mean,
            thermal_std,
            reference_mask,
            THERMAL_STD_UNSTABLE,
            reverse_tail_codes=False,
        )
        thermal_mean_path = output_dir / "thermal_mean_zscore.tif"
        thermal_std_path = output_dir / "thermal_std_zscore.tif"
        thermal_class_path = output_dir / "thermal_stability.tif"
        thermal_preview_path = output_dir / "thermal_stability_preview.png"
        write_float_tif(thermal_mean_path, thermal_mean, reference_grid)
        write_float_tif(thermal_std_path, thermal_std, reference_grid)
        write_class_tif(thermal_class_path, thermal_classes, reference_grid)
        save_class_preview(
            classes=thermal_classes,
            path=thermal_preview_path,
            title=f"{reference_record.field_name} | {group_name} | Thermal stability",
            names=THERMAL_CLASS_NAMES,
            colors=THERMAL_COLORS,
            lower_cut=cold_cut,
            upper_cut=hot_cut,
            subtitle=(
                f"{len(thermal_files)} dates across {len(years)} years; "
                f"Unstable when std > {THERMAL_STD_UNSTABLE}"
            ),
        )
        write_rows(output_dir / "thermal_files_used.csv", thermal_files)
        rows.append(
            {
                "field": reference_record.field_name,
                "group": group_name,
                "map_type": "thermal_stability",
                "years": ",".join(map(str, years)),
                "crops": ",".join(sorted({record.crop for record in records})),
                "layers_used": len(thermal_files),
                "reference_year": reference_record.year,
                "unstable_threshold": THERMAL_STD_UNSTABLE,
                "lower_percentile_cut": cold_cut,
                "upper_percentile_cut": hot_cut,
                "unstable_pixels": int(np.sum(thermal_classes == 0)),
                "class_1_pixels": int(np.sum(thermal_classes == 1)),
                "medium_pixels": int(np.sum(thermal_classes == 2)),
                "class_3_pixels": int(np.sum(thermal_classes == 3)),
                "output_tif": str(thermal_class_path),
            }
        )

    # ---- Multi-year yield stability, using exactly the same eligible years ----
    yield_stack, yield_files = build_yield_stack(
        records, reference_grid, reference_mask
    )
    if yield_stack is not None and yield_stack.shape[0] >= MIN_YEARS_FOR_STABILITY:
        yield_mean, yield_std = stack_mean_std(yield_stack)
        yield_classes, low_cut, high_cut = classify_stability(
            yield_mean,
            yield_std,
            reference_mask,
            YIELD_STD_UNSTABLE,
            reverse_tail_codes=True,
        )
        yield_mean_path = output_dir / "yield_mean_zscore.tif"
        yield_std_path = output_dir / "yield_std_zscore.tif"
        yield_class_path = output_dir / "yield_stability.tif"
        yield_preview_path = output_dir / "yield_stability_preview.png"
        write_float_tif(yield_mean_path, yield_mean, reference_grid)
        write_float_tif(yield_std_path, yield_std, reference_grid)
        write_class_tif(yield_class_path, yield_classes, reference_grid)
        save_class_preview(
            classes=yield_classes,
            path=yield_preview_path,
            title=f"{reference_record.field_name} | {group_name} | Yield stability",
            names=YIELD_CLASS_NAMES,
            colors=YIELD_COLORS,
            lower_cut=low_cut,
            upper_cut=high_cut,
            subtitle=(
                f"{len(yield_files)} years; Unstable when std > "
                f"{YIELD_STD_UNSTABLE}; codes align Cold/High and Hot/Low"
            ),
        )
        write_rows(output_dir / "yield_files_used.csv", yield_files)
        rows.append(
            {
                "field": reference_record.field_name,
                "group": group_name,
                "map_type": "yield_stability",
                "years": ",".join(map(str, years)),
                "crops": ",".join(sorted({record.crop for record in records})),
                "layers_used": len(yield_files),
                "reference_year": reference_record.year,
                "unstable_threshold": YIELD_STD_UNSTABLE,
                "lower_percentile_cut": low_cut,
                "upper_percentile_cut": high_cut,
                "unstable_pixels": int(np.sum(yield_classes == 0)),
                "class_1_pixels": int(np.sum(yield_classes == 1)),
                "medium_pixels": int(np.sum(yield_classes == 2)),
                "class_3_pixels": int(np.sum(yield_classes == 3)),
                "output_tif": str(yield_class_path),
            }
        )
    else:
        print(
            f"  SKIP {group_name} yield stability: fewer than "
            f"{MIN_YEARS_FOR_STABILITY} valid normalized yield layers."
        )

    return rows


# =============================================================================
# PROGRAM ENTRY POINT
# =============================================================================


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create yearly thermal maps and multi-year thermal/yield "
            "stability maps without distributing the source data."
        )
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=INPUT_ROOT,
        help="Root containing one folder per field (default: edit INPUT_ROOT).",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=OUTPUT_ROOT,
        help="Directory for generated outputs (default: edit OUTPUT_ROOT).",
    )
    return parser.parse_args()


def validate_paths(input_root: Path, output_root: Path) -> None:
    """Fail with a clear message when placeholders were not replaced."""

    if "INSERT_YOUR" in str(input_root) or "INSERT_YOUR" in str(output_root):
        raise ValueError(
            "Replace INPUT_ROOT and OUTPUT_ROOT at the top of the script, or "
            "pass --input-root and --output-root."
        )
    if not input_root.is_dir():
        raise FileNotFoundError(f"Input root does not exist: {input_root}")
    output_root.mkdir(parents=True, exist_ok=True)


def run_pipeline(input_root: Path, output_root: Path) -> None:
    """Run all requested products for every field under the input root."""

    validate_paths(input_root, output_root)
    field_dirs = sorted(
        path
        for path in input_root.iterdir()
        if path.is_dir() and path.resolve() != output_root.resolve()
    )
    print(f"Found {len(field_dirs)} candidate field folder(s) in {input_root}")

    yearly_rows: list[dict] = []
    stability_rows: list[dict] = []
    for field_dir in field_dirs:
        print(f"\n=== FIELD: {field_dir.name} ===")
        records = discover_records(field_dir)
        if not records:
            print("  No eligible field-years.")
            continue

        field_output = output_root / field_dir.name
        field_output.mkdir(parents=True, exist_ok=True)

        validated_records: list[FieldYearRecord] = []
        for record in records:
            print(
                f"  Yearly thermal map {record.year} ({record.crop}): "
                f"{len(record.thermal_paths)} thermal dates"
            )
            year_result = create_yearly_thermal_map(record, field_output)
            if year_result is not None:
                yearly_rows.append(year_result)
                validated_records.append(record)

        for group_name in STABILITY_GROUPS:
            selected = group_records(validated_records, group_name)
            print(f"  Stability group {group_name}: {len(selected)} eligible year(s)")
            stability_rows.extend(
                create_stability_maps(selected, group_name, field_output)
            )

    write_rows(output_root / "yearly_thermal_summary_all_fields.csv", yearly_rows)
    write_rows(output_root / "stability_summary_all_fields.csv", stability_rows)
    print(f"\nCompleted. Outputs: {output_root}")


def main() -> None:
    arguments = parse_arguments()
    try:
        run_pipeline(arguments.input_root.expanduser(), arguments.output_root.expanduser())
    except (ValueError, FileNotFoundError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2) from error


if __name__ == "__main__":
    main()
