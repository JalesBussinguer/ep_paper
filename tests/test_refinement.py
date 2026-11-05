import numpy as np
import pandas as pd

from src import pixel_analysis as pa
from src.gi_robustness import benjamini_hochberg, classify_fdr, cohen_kappa
from src.pixel_extraction import YEARS, decode_keys


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
