"""Small, data-free checks of the WSI definitions used in the paper."""

import unittest
import json
import tempfile
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import rasterio
from rasterio.io import MemoryFile
from rasterio.transform import from_origin

from calculate_wsi import class_fractions, hot_cold_gap, main, summarize_field


class WsiTests(unittest.TestCase):
    def test_period_ratio_and_complete_pixel_requirement(self):
        # Pixel (0, 2) lacks the second date and must not enter class means.
        dates = pd.to_datetime(["2016-07-01", "2016-07-02"])
        samples = pd.DataFrame([
            ("Cold-dominant", 0, 0, dates[0], 2.0),
            ("Cold-dominant", 0, 0, dates[1], 2.0),
            ("Hot-dominant", 0, 1, dates[0], 1.0),
            ("Hot-dominant", 0, 1, dates[1], 1.0),
            ("Hot-dominant", 0, 2, dates[0], 4.0),
        ], columns=["class", "row", "col", "date", "et_mm_day"])
        samples["field"] = "Field_1"
        samples["year"] = 2016
        samples["crop"] = "corn"
        samples["class_fraction"] = 0.75
        samples["pet_mm_day"] = 4.0
        samples["rain_mm_day"] = 0.0
        weather = pd.DataFrame({"year": [2016, 2016], "rain_mm_day": [2.0, 3.0]})
        pixels, summary = summarize_field(samples, weather)
        self.assertEqual(len(pixels), 2)
        values = summary.set_index("class")["wsi_mean"]
        self.assertAlmostEqual(values["Cold-dominant"], 0.5)
        self.assertAlmostEqual(values["Hot-dominant"], 0.75)
        self.assertAlmostEqual(hot_cold_gap(summary)["hot_minus_cold_wsi"].iloc[0], 0.25)
        self.assertEqual(summary["july_august_rain_mm"].iloc[0], 5.0)

    def test_binary_average_is_fractional_not_nearest_neighbor(self):
        # Four 1-m cells comprise one 2-m ET cell: 3 cold, 1 hot.
        source = np.array([[1, 1], [1, 3]], dtype="uint8")
        with MemoryFile() as source_memory, MemoryFile() as et_memory:
            with source_memory.open(driver="GTiff", height=2, width=2, count=1,
                                    dtype="uint8", crs="EPSG:32616",
                                    transform=from_origin(0, 2, 1, 1)) as src:
                src.write(source, 1)
            with et_memory.open(driver="GTiff", height=1, width=1, count=1,
                                dtype="float32", crs="EPSG:32616",
                                transform=from_origin(0, 2, 2, 2)):
                pass
            with et_memory.open() as et:
                fractions = class_fractions(source_memory.name, et)
                self.assertAlmostEqual(fractions["Cold-dominant"][0, 0], 0.75)
                self.assertAlmostEqual(fractions["Hot-dominant"][0, 0], 0.25)

    def test_cli_omits_pet_for_missing_openet_date(self):
        # The weather series has three days; OpenET has only days 1 and 3.
        # A two-pixel grid keeps one cold and one hot observation per date.
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            et_dir = root / "openet"
            et_dir.mkdir()
            transform = from_origin(0, 1, 1, 1)
            for path, values, dtype in (
                (root / "stability.tif", [[1, 3]], "uint8"),
                (root / "dem.tif", [[250.0, 250.0]], "float32"),
                (et_dir / "ensemble_et_2016_07_01_test.tif", [[2.0, 1.0]], "float32"),
                (et_dir / "ensemble_et_2016_07_03_test.tif", [[2.0, 1.0]], "float32"),
            ):
                with rasterio.open(path, "w", driver="GTiff", height=1, width=2,
                                   count=1, dtype=dtype, crs="EPSG:32616",
                                   transform=transform) as dst:
                    dst.write(np.asarray(values, dtype=dtype), 1)
            pd.DataFrame({"Year": [2016] * 3, "DOY": [183, 184, 185],
                          "SRAD": [20.0] * 3, "Tmax": [30.0] * 3,
                          "Tmin": [15.0] * 3, "Rain": [0.0, 5.0, 0.0],
                          "WindSp": [2.0] * 3, "Rmax": [90.0] * 3,
                          "Rmin": [40.0] * 3}).to_csv(root / "weather.csv", index=False)
            config = {"fields": [{"id": "TestField", "latitude_degrees": 40.0,
                                  "weather_csv": "weather.csv", "dem_raster": "dem.tif",
                                  "thermal_stability_raster": "stability.tif",
                                  "openet_dir": "openet", "crop_by_year": {"2016": "corn"}}]}
            (root / "config.json").write_text(json.dumps(config), encoding="utf-8")
            args = ["calculate_wsi.py", "--config", str(root / "config.json"),
                    "--input-root", str(root), "--output-root", str(root / "outputs")]
            with patch("sys.argv", args), redirect_stdout(StringIO()):
                main()
            output = root / "outputs" / "TestField"
            summary = pd.read_csv(output / "field_year_wsi.csv")
            pet = pd.read_csv(output / "daily_pet_weather.csv")
            qc = pd.read_csv(output / "date_qc.csv")
            self.assertEqual(set(summary["required_days"]), {2})
            self.assertEqual(len(pet), 3)
            self.assertEqual(qc.loc[qc.date == "2016-07-02", "status"].iloc[0], "missing_openet")
            matched_pet = pet.loc[pet.date.isin(["2016-07-01", "2016-07-03"]), "pet_mm_day"].sum()
            self.assertAlmostEqual(summary["pet_total_mm"].iloc[0], matched_pet)
            self.assertEqual(summary["july_august_rain_mm"].iloc[0], 5.0)


if __name__ == "__main__":
    unittest.main()
