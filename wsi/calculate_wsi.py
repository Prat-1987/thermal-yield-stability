#!/usr/bin/env python3
"""Reproduce July-August native-grid OpenET hot/cold WSI estimates.

All input paths come from a JSON manifest; no field data are bundled here.
The equations and aggregation follow the analysis scripts used for the paper.
See README.md for the scientific assumptions and output definitions.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject


FRACTION_THRESHOLD = 0.60
MONTHS = (7, 8)
WEATHER_COLUMNS = {"Year", "DOY", "SRAD", "Tmax", "Tmin", "Rain", "WindSp", "Rmax", "Rmin"}
DATE_PATTERN = re.compile(r"ensemble_et_(\d{4})_(\d{2})_(\d{2})(?:_|\.)")
CLASS_NAMES = ("Cold-dominant", "Hot-dominant")


def resolve_input(root: Path, value: str, field: str) -> Path:
    """Resolve a manifest path and report which field has a missing input."""
    path = root / value
    if not path.exists():
        raise FileNotFoundError(f"{field}: required input does not exist: {path}")
    return path


def read_config(path: Path, input_root: Path) -> list[dict]:
    with path.open(encoding="utf-8") as stream:
        config = json.load(stream)
    fields = config.get("fields")
    if not isinstance(fields, list) or not fields:
        raise ValueError("Configuration must contain a nonempty 'fields' list.")
    result = []
    seen = set()
    for entry in fields:
        field = str(entry.get("id", ""))
        if not re.fullmatch(r"[A-Za-z0-9_-]+", field) or field in seen:
            raise ValueError(f"Field ID must be unique and filename-safe: {field!r}")
        seen.add(field)
        try:
            latitude = float(entry["latitude_degrees"])
            crops = {int(k): str(v).lower() for k, v in entry["crop_by_year"].items()}
            paths = {key: resolve_input(input_root, entry[key], field) for key in
                     ("weather_csv", "dem_raster", "thermal_stability_raster", "openet_dir")}
        except KeyError as exc:
            raise ValueError(f"{field}: missing manifest key {exc.args[0]!r}") from exc
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"{field}: replace latitude_degrees with a numeric field latitude and "
                "check crop_by_year contains numeric years."
            ) from exc
        if not -90 < latitude < 90:
            raise ValueError(f"{field}: latitude_degrees must be between -90 and 90.")
        if not crops or any(c not in {"corn", "soybean"} for c in crops.values()):
            raise ValueError(f"{field}: crop_by_year must list only corn/soybean years.")
        result.append({"id": field, "latitude": latitude, "crops": crops, **paths})
    return result


def elevation_from_dem(path: Path) -> float:
    """Use the mean of valid DEM cells, not the elevation of an arbitrary point."""
    with rasterio.open(path) as src:
        if src.count != 1:
            raise ValueError(f"DEM must be single-band: {path}")
        dem = src.read(1, masked=True).astype("float32").filled(np.nan)
    elevation = float(np.nanmean(dem))
    if not np.isfinite(elevation):
        raise ValueError(f"DEM has no finite elevation values: {path}")
    return elevation


def extraterrestrial_radiation(doy: int, latitude_degrees: float) -> float:
    """FAO-56 daily extraterrestrial radiation, MJ m-2 d-1."""
    lat = math.radians(latitude_degrees)
    dr = 1.0 + 0.033 * math.cos(2.0 * math.pi * doy / 365.0)
    declination = 0.409 * math.sin(2.0 * math.pi * doy / 365.0 - 1.39)
    sunset = math.acos(-math.tan(lat) * math.tan(declination))
    return (24.0 * 60.0 / math.pi) * 0.0820 * dr * (
        sunset * math.sin(lat) * math.sin(declination)
        + math.cos(lat) * math.cos(declination) * math.sin(sunset)
    )


def saturation_pressure(temp_c):
    """Saturation vapour pressure at air temperature T, in kPa."""
    return 0.6108 * np.exp(17.27 * temp_c / (temp_c + 237.3))


def calculate_daily_pet(weather_path: Path, dem_path: Path, latitude: float) -> pd.DataFrame:
    """Compute the paper's FAO-56 PM ET0 (called PET in the original analysis).

    Input CSV values are assumed already in Celsius, MJ m-2 d-1, m s-1,
    percent RH, and mm d-1. No Kelvin or W m-2 conversion is performed.
    WindSp is used directly as u2: the original analysis did not convert
    nominal gridMET 10-m wind to the FAO-56 2-m convention. Daily G = 0.
    """
    weather = pd.read_csv(weather_path)
    missing = WEATHER_COLUMNS - set(weather.columns)
    if missing:
        raise ValueError(f"Weather CSV {weather_path} is missing columns: {sorted(missing)}")
    weather = weather.copy()
    weather["date"] = pd.to_datetime(weather["Year"].astype(str), format="%Y") + pd.to_timedelta(
        weather["DOY"].astype(int) - 1, unit="D"
    )
    if weather["date"].duplicated().any():
        raise ValueError(f"Weather CSV has duplicate dates: {weather_path}")
    for name in WEATHER_COLUMNS - {"Year", "DOY"}:
        weather[name] = pd.to_numeric(weather[name], errors="raise")
    z = elevation_from_dem(dem_path)
    tmax = weather["Tmax"].astype(float)
    tmin = weather["Tmin"].astype(float)
    tmean = (tmax + tmin) / 2.0
    rs = weather["SRAD"].astype(float)
    wind = weather["WindSp"].astype(float)
    rhmax = weather["Rmax"].astype(float).clip(0, 100)
    rhmin = weather["Rmin"].astype(float).clip(0, 100)

    # FAO-56 two-temperature and two-relative-humidity vapour pressures.
    e_tmax, e_tmin = saturation_pressure(tmax), saturation_pressure(tmin)
    es = (e_tmax + e_tmin) / 2.0
    ea = (e_tmin * rhmax / 100.0 + e_tmax * rhmin / 100.0) / 2.0
    vpd = (es - ea).clip(lower=0)
    slope = 4098.0 * saturation_pressure(tmean) / ((tmean + 237.3) ** 2)

    # Elevation enters both atmospheric pressure and clear-sky radiation.
    pressure = 101.3 * ((293.0 - 0.0065 * z) / 293.0) ** 5.26
    gamma = 0.000665 * pressure
    ra = weather["DOY"].apply(lambda d: extraterrestrial_radiation(int(d), latitude))
    rso = (0.75 + 2e-5 * z) * ra
    rns = (1.0 - 0.23) * rs
    rs_rso = (rs / rso).replace([np.inf, -np.inf], np.nan).clip(0.3, 1.0)
    rnl = (4.903e-9 * ((tmax + 273.16) ** 4 + (tmin + 273.16) ** 4) / 2.0
           * (0.34 - 0.14 * np.sqrt(ea.clip(lower=0))) * (1.35 * rs_rso - 0.35))
    rn = rns - rnl

    # FAO-56 Penman-Monteith reference ET with daily soil heat flux G = 0.
    pet = ((0.408 * slope * rn + gamma * (900.0 / (tmean + 273.0)) * wind * vpd)
           / (slope + gamma * (1.0 + 0.34 * wind))).clip(lower=0)
    return pd.DataFrame({"date": weather["date"], "year": weather["Year"].astype(int),
                         "doy": weather["DOY"].astype(int), "rain_mm_day": weather["Rain"],
                         "pet_mm_day": pet, "rn_mj_m2_day": rn, "elevation_m": z,
                         "wind_input_m_s": wind, "vpd_kpa": vpd})


def openet_files(directory: Path, years: set[int]) -> dict[pd.Timestamp, Path]:
    """Index July-August ensemble ET rasters and reject ambiguous dates."""
    files = {}
    for path in sorted(directory.rglob("*.tif")):
        match = DATE_PATTERN.search(path.name)
        if not match:
            continue
        try:
            date = pd.Timestamp("-".join(match.groups()))
        except ValueError as exc:
            raise ValueError(f"Invalid date in OpenET filename: {path}") from exc
        if date.year not in years or date.month not in MONTHS:
            continue
        if date in files:
            raise ValueError(f"Duplicate OpenET date {date.date()}: {files[date]} and {path}")
        files[date] = path
    if not files:
        raise FileNotFoundError(f"No July-August ensemble_et_YYYY_MM_DD TIFFs for configured years in {directory}")
    return files


def class_fractions(stability_path: Path, et_src) -> dict[str, np.ndarray]:
    """Average 2-m binary class masks onto the native OpenET grid.

    Source NoData is represented as zero in each mask, as in the original
    workflow. It therefore contributes to the 30-m denominator. Do not
    nearest-neighbour assign classes or interpolate the continuous ET raster.
    """
    with rasterio.open(stability_path) as src:
        if src.count != 1 or src.crs is None or et_src.crs is None:
            raise ValueError("Stability and OpenET rasters must be single-band and have a CRS.")
        classes = src.read(1)
        valid = np.ones(classes.shape, dtype=bool)
        if src.nodata is not None:
            valid &= classes != src.nodata
        fractions = {}
        for name, code in (("Cold-dominant", 1), ("Hot-dominant", 3)):
            binary = ((classes == code) & valid).astype("float32")
            fraction = np.zeros((et_src.height, et_src.width), dtype="float32")
            reproject(source=binary, destination=fraction,
                      src_transform=src.transform, src_crs=src.crs,
                      dst_transform=et_src.transform, dst_crs=et_src.crs,
                      src_nodata=None, dst_nodata=0, resampling=Resampling.average)
            fractions[name] = fraction
        return fractions


def raster_grid(src) -> tuple:
    """All dates within one field must share the native OpenET grid."""
    return str(src.crs), tuple(src.transform), src.height, src.width


def sample_field(field: dict, out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    name = field["id"]
    weather = calculate_daily_pet(field["weather_csv"], field["dem_raster"], field["latitude"])
    weather = weather[weather["year"].isin(field["crops"]) & weather["date"].dt.month.isin(MONTHS)]
    weather.sort_values("date").to_csv(out_dir / "daily_pet_weather.csv", index=False)
    weather_by_date = weather.set_index("date")
    files = openet_files(field["openet_dir"], set(field["crops"]))
    samples, qc = [], []
    reference_grid, fractions = None, None

    for year in sorted(field["crops"]):
        for date in pd.date_range(f"{year}-07-01", f"{year}-08-31"):
            path = files.get(date)
            record = {"field": name, "year": year, "crop": field["crops"][year],
                      "date": date.date().isoformat(), "openet_file": str(path) if path else "",
                      "cold_valid_pixels": 0, "hot_valid_pixels": 0}
            if path is None:
                record["status"] = "missing_openet"
                qc.append(record)
                continue
            if date not in weather_by_date.index:
                record["status"] = "missing_weather"
                qc.append(record)
                continue
            wx = weather_by_date.loc[date]
            pet = float(wx["pet_mm_day"])
            if not np.isfinite(pet) or pet <= 0:
                record["status"] = "invalid_pet"
                qc.append(record)
                continue
            with rasterio.open(path) as src:
                if src.count != 1 or src.crs is None:
                    raise ValueError(f"OpenET must be a georeferenced single-band TIFF: {path}")
                if not np.isclose(src.scales[0], 1.0) or not np.isclose(src.offsets[0], 0.0):
                    raise ValueError(f"OpenET raster has scale/offset metadata not used by this workflow: {path}")
                grid = raster_grid(src)
                if reference_grid is None:
                    reference_grid = grid
                    fractions = class_fractions(field["thermal_stability_raster"], src)
                    for cls in CLASS_NAMES:
                        if not np.any(fractions[cls] >= FRACTION_THRESHOLD):
                            raise ValueError(f"{name}: no {cls} pixels meet the {FRACTION_THRESHOLD} threshold on the OpenET grid")
                elif grid != reference_grid:
                    raise ValueError(f"{name}: OpenET grids differ; expected {reference_grid}, found {grid} in {path}")
                et = src.read(1).astype("float32")
                if src.nodata is not None:
                    et = np.where(et == src.nodata, np.nan, et)
                et = np.where(et < 0, np.nan, et)
            rows, cols = np.indices(et.shape)
            for cls in CLASS_NAMES:
                valid = (fractions[cls] >= FRACTION_THRESHOLD) & np.isfinite(et)
                record["cold_valid_pixels" if cls == "Cold-dominant" else "hot_valid_pixels"] = int(valid.sum())
                if not valid.any():
                    continue
                samples.append(pd.DataFrame({"field": name, "year": year, "crop": field["crops"][year],
                                             "date": date, "class": cls, "row": rows[valid], "col": cols[valid],
                                             "class_fraction": fractions[cls][valid], "et_mm_day": et[valid],
                                             "pet_mm_day": pet, "rain_mm_day": float(wx["rain_mm_day"])}))
            record["status"] = "sampled" if record["cold_valid_pixels"] or record["hot_valid_pixels"] else "no_valid_et"
            qc.append(record)
    pd.DataFrame(qc).to_csv(out_dir / "date_qc.csv", index=False)
    if not samples:
        raise ValueError(f"{name}: no valid cold/hot OpenET pixels were sampled; inspect date_qc.csv")
    return pd.concat(samples, ignore_index=True), weather


def summarize_field(samples: pd.DataFrame, weather: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Ratio of period totals per pixel, then equal-weight class mean.

    The denominator includes PET only for dates with an OpenET sample. As in
    the paper analysis, a pixel must have valid ET on all sampled dates in
    its field-year. This is not the mean of daily WSI values.
    """
    required = samples.groupby("year")["date"].nunique().rename("required_days")
    pixels = (samples.groupby(["field", "year", "crop", "class", "row", "col"], as_index=False)
              .agg(days_present=("date", "nunique"), class_fraction=("class_fraction", "mean"),
                   et_total_mm=("et_mm_day", "sum"), pet_total_mm=("pet_mm_day", "sum")))
    pixels = pixels.join(required, on="year")
    pixels = pixels[pixels["days_present"] == pixels["required_days"]].copy()
    if pixels.empty:
        raise ValueError("No pixels have valid OpenET ET on every sampled date; inspect date_qc.csv")
    pixels["et_pet"] = pixels["et_total_mm"] / pixels["pet_total_mm"]
    pixels["wsi"] = 1.0 - pixels["et_pet"]
    pixels["et_deficit_mm"] = pixels["pet_total_mm"] - pixels["et_total_mm"]

    summary = (pixels.groupby(["field", "year", "crop", "class"], as_index=False)
               .agg(n_pixels=("wsi", "size"), mean_class_fraction=("class_fraction", "mean"),
                    et_total_mm_mean=("et_total_mm", "mean"), pet_total_mm=("pet_total_mm", "first"),
                    et_pet_mean=("et_pet", "mean"), wsi_mean=("wsi", "mean"),
                    et_deficit_mm_mean=("et_deficit_mm", "mean"), required_days=("required_days", "first")))
    # Rainfall context uses all July-August weather dates, including a date
    # for which OpenET is missing. Rain was not an input to the PET equation.
    rain = weather.groupby("year")["rain_mm_day"].sum().rename("july_august_rain_mm")
    summary = summary.join(rain, on="year")
    return pixels, summary


def hot_cold_gap(summary: pd.DataFrame) -> pd.DataFrame:
    """Compare class means within each field-year; positive gap means hotter is drier."""
    index = ["field", "year", "crop"]
    wide = summary.pivot(index=index, columns="class", values="wsi_mean")
    for cls in CLASS_NAMES:
        if cls not in wide:
            wide[cls] = np.nan
    wide = wide.rename(columns={"Cold-dominant": "cold_wsi_mean", "Hot-dominant": "hot_wsi_mean"})
    wide["hot_minus_cold_wsi"] = wide["hot_wsi_mean"] - wide["cold_wsi_mean"]
    context = summary.drop_duplicates(index).set_index(index)[["july_august_rain_mm", "required_days"]]
    return wide.join(context).reset_index()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True, help="JSON field manifest")
    parser.add_argument("--input-root", type=Path, required=True, help="Base directory for paths in the manifest")
    parser.add_argument("--output-root", type=Path, required=True, help="Directory for generated CSV files")
    args = parser.parse_args()
    fields = read_config(args.config, args.input_root)
    args.output_root.mkdir(parents=True, exist_ok=True)
    combined, gaps = [], []
    for field in fields:
        out_dir = args.output_root / field["id"]
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"Analyzing {field['id']} ({len(field['crops'])} configured crop-years)", flush=True)
        samples, weather = sample_field(field, out_dir)
        pixels, summary = summarize_field(samples, weather)
        pixels.to_csv(out_dir / "pixel_year_wsi.csv", index=False)
        summary.to_csv(out_dir / "field_year_wsi.csv", index=False)
        gap = hot_cold_gap(summary)
        gap.to_csv(out_dir / "hot_cold_wsi_gap.csv", index=False)
        combined.append(summary)
        gaps.append(gap)
        print(summary[["year", "class", "n_pixels", "wsi_mean"]].to_string(index=False), flush=True)
    pd.concat(combined, ignore_index=True).to_csv(args.output_root / "all_fields_field_year_wsi.csv", index=False)
    pd.concat(gaps, ignore_index=True).to_csv(args.output_root / "all_fields_hot_cold_wsi_gap.csv", index=False)


if __name__ == "__main__":
    main()
