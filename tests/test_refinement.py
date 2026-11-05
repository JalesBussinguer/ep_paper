import numpy as np
import pandas as pd
from PIL import Image
from shapely.geometry import box

from src import pixel_analysis as pa
from src import visual_validation
from src.gi_robustness import benjamini_hochberg, classify_fdr, cohen_kappa
from src.pixel_extraction import YEARS, decode_keys
from src.visualization import map_scale_length_km


def test_decode_keys_roundtrip():
    classes = np.array([[3, 4, 15, 15, 15, 21, 21]], dtype=np.uint8)
    key = sum(int(value) << (8 * index) for index, value in enumerate(classes[0]))
    decoded = decode_keys(np.array([key], dtype=np.uint64))
    assert decoded.shape == (1, len(YEARS))
    assert decoded.tolist() == classes.tolist()


def test_typology_categories():
    n, a = pa.NATURAL, pa.ANTHROPIC
    pre = np.array([n, n, n, a, a, n])
    det = np.array([a, n, n, a, n, n])
    post = np.array([a, a, n, a, n, n])
    missing = np.array([False, False, False, False, False, True])
    assert pa.typology(pre, det, post, missing).tolist() == [0, 1, 2, 3, 4, 6]


def test_cramers_v_bounds():
    assert pa.cramers_v([[50, 0], [0, 50]]) == 1.0
    assert pa.cramers_v([[25, 25], [25, 25]]) == 0.0


def test_bootstrap_draws_preserve_expected_total():
    rng = np.random.default_rng(0)
    matrix = np.tile(np.array([[1.0, 2.0]]), (50, 1))
    draws = pa.bootstrap_draws(matrix, rng, n_boot=20)
    assert draws.shape == (20, 2)
    assert np.allclose(draws, [[50.0, 100.0]])


def test_benjamini_hochberg_monotonic_and_bounded():
    p = np.array([0.001, 0.01, 0.03, 0.2, 0.9])
    adjusted = benjamini_hochberg(p)
    assert np.all(adjusted >= p)
    assert np.all(adjusted <= 1.0)
    assert list(adjusted) == sorted(adjusted)


def test_fdr_classification_and_kappa():
    z = np.array([5.0, -5.0, 0.1, 0.0])
    labels, _ = classify_fdr(z, None)
    assert labels[0].startswith("Hot") and labels[1].startswith("Cold")
    assert labels[2] == "Não significante"
    assert cohen_kappa([1, 0, 1, 0], [1, 0, 1, 0]) == 1.0


def test_map_scale_length_is_round_and_fits_map():
    assert map_scale_length_km(1_000_000) == 100
    assert map_scale_length_km(2_500_000) == 200


def test_visual_validation_sheets_keep_same_dimensions(monkeypatch, tmp_path):
    class RasterStub:
        crs = "EPSG:4326"

        def close(self):
            pass

    monkeypatch.setattr(visual_validation, "raster_paths", lambda: ["raster"] * len(visual_validation.config.MAPBIOMAS_YEARS))
    monkeypatch.setattr(visual_validation.rasterio, "open", lambda _: RasterStub())
    monkeypatch.setattr(visual_validation, "read_panel", lambda *_: np.zeros((2, 2), dtype=np.uint8))
    add_page_elements = visual_validation.add_page_map_elements
    page_crs = []
    monkeypatch.setattr(
        visual_validation,
        "add_page_map_elements",
        lambda fig, crs: (page_crs.append(crs), add_page_elements(fig, crs)),
    )

    geometry = box(0, 0, 1, 1)
    geometries = {1: geometry, 2: geometry, 3: geometry}
    sizes = []
    for sample_ids in ([1], [1, 2, 3]):
        rows = pd.DataFrame({"sample_id": sample_ids, "ano_prodes": [2020] * len(sample_ids)})
        path = tmp_path / f"sheet_{len(sample_ids)}.png"
        visual_validation.draw_sheet(rows, geometries, path, "Validação")
        with Image.open(path) as image:
            sizes.append(image.size)

    assert sizes[0] == sizes[1]
    assert visual_validation.PER_SHEET == 3
    assert visual_validation.SHEET_ROWS == visual_validation.SHEET_COLUMNS == 3
    assert len(page_crs) == 2
