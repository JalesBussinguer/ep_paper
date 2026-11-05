import tempfile
import unittest
from pathlib import Path

import pandas as pd

from src.mapbiomas_event_transition_analysis import (
    INPUT_COLUMNS,
    TRANSITION_PAIRS,
    _aggregate_event_transitions,
    _main_outputs,
    _parse_class_codes,
)


class EventTransitionAnalysisTests(unittest.TestCase):
    def test_missing_class_is_not_converted_to_sem_dados(self):
        classes = pd.Series(["", "Sem dados", None], dtype="string")

        codes = _parse_class_codes(classes)

        self.assertTrue(codes.isna().iloc[0])
        self.assertEqual(int(codes.iloc[1]), 0)
        self.assertTrue(codes.isna().iloc[2])

    def test_transitions_use_valid_denominators_and_include_empty_cohort(self):
        rows = [
            [
                2020, "1", "Mun A", "MA", "Maranhao", 2.0,
                "3 - Formação Florestal", "15 - Pastagem", "15 - Pastagem",
            ],
            [
                2020, "1", "Mun A", "MA", "Maranhao", 1.0,
                "4 - Formação Savânica", "4 - Formação Savânica", "15 - Pastagem",
            ],
            [
                2021, "2", "Mun B", "TO", "Tocantins", 3.0,
                "15 - Pastagem", "21 - Mosaico de Agricultura e Pastagem", None,
            ],
            [
                2021, "2", "Mun B", "TO", "Tocantins", 4.0,
                "Sem dados", "Sem dados", "Sem dados",
            ],
            [
                2025, "3", "Mun C", "GO", "Goias", 5.0,
                "3 - Formação Florestal", "15 - Pastagem", None,
            ],
        ]

        with tempfile.TemporaryDirectory() as temporary_directory:
            input_path = Path(temporary_directory) / "events.csv"
            pd.DataFrame(rows, columns=INPUT_COLUMNS).to_csv(
                input_path,
                index=False,
                encoding="utf-8",
            )
            outputs = _main_outputs(_aggregate_event_transitions(input_path))

        summary, classes, groups, municipal, trajectories = outputs

        def summary_row(interval, year="Todos"):
            return summary.loc[
                (summary["intervalo"] == interval)
                & (summary["ano_prodes"] == year)
            ].iloc[0]

        self.assertEqual(
            summary_row(TRANSITION_PAIRS[0][0]).n_intersecoes_validas,
            5,
        )
        self.assertEqual(
            summary_row(TRANSITION_PAIRS[0][0]).n_classes_dominantes_alteradas,
            3,
        )
        self.assertEqual(
            summary_row(TRANSITION_PAIRS[1][0]).n_intersecoes_validas,
            3,
        )
        self.assertEqual(
            summary_row(TRANSITION_PAIRS[2][0]).n_intersecoes_validas,
            3,
        )
        empty_post_cohort = summary_row(TRANSITION_PAIRS[1][0], 2025)
        self.assertEqual(empty_post_cohort.n_intersecoes_validas, 0)
        self.assertEqual(empty_post_cohort.area_intersecoes_validas_km2, 0)

        for frame in (classes, groups):
            for _, cohort in frame.groupby(["intervalo", "ano_prodes"]):
                self.assertAlmostEqual(cohort["percentual_intersecoes"].sum(), 100)
                self.assertAlmostEqual(
                    cohort["percentual_area_intersecoes"].sum(),
                    100,
                )
        for _, cohort in municipal.groupby(
            ["CD_MUN", "ano_prodes", "intervalo"]
        ):
            self.assertAlmostEqual(cohort["percentual_intersecoes"].sum(), 100)
            self.assertAlmostEqual(
                cohort["percentual_area_intersecoes"].sum(),
                100,
            )
        self.assertEqual(
            trajectories.loc[trajectories["ano_prodes"] == "Todos",
                             "n_intersecoes"].sum(),
            3,
        )


if __name__ == "__main__":
    unittest.main()
