"""Regression checks for projection completeness and certificate rejection."""

import copy
from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest

from verify_split import (
    ORDERS, ZERO, add, band_column, build_system, check_identity, criterion,
    positive_circuits, project, scale, unit, verify,
)


class ProjectionTests(unittest.TestCase):
    def test_four_row_circuit_and_degenerate_supports(self):
        directions = [(-1, 1, 1), (0, -1, 0), (0, 0, -1), (1, 0, 0)]
        self.assertEqual(list(positive_circuits(directions)), [(tuple(directions), (1, 1, 1, 1))])
        self.assertEqual(list(positive_circuits([(1, 0, 0), (0, 1, 0)])), [])
        # A rank-one triple is not a circuit; its two opposite pairs are.
        rays = list(positive_circuits([(1, 0, 0), (2, 0, 0), (-1, 0, 0)]))
        self.assertEqual(len(rays), 2)
        self.assertTrue(all(len(support) == 2 for support, _ in rays))

    def test_all_rows_with_the_same_direction_and_zero_rows(self):
        rows = [add(unit(0), unit(37)), add(unit(1), unit(37)),
                add(unit(2), unit(38)), add(unit(3), scale(-1, unit(37)), scale(-1, unit(38))),
                unit(4)]
        projected, directions, circuits = project(rows)
        self.assertEqual(directions, 3)
        self.assertEqual(circuits, 1)
        self.assertEqual(set(projected), {
            add(unit(0), unit(2), unit(3))[:37],
            add(unit(1), unit(2), unit(3))[:37], unit(4)[:37],
        })

    def test_exact_weights_and_constant_entry(self):
        hypothesis = scale(2, unit(36))[:37]
        rho = unit(36)[:37]
        check_identity([hypothesis], rho, [[0, 1, 2]])
        for terms in ([[0, -1, 2]], [[0, 1, 3]], [[0, 1, 0]], [[1, 1, 2]]):
            with self.assertRaises(ValueError):
                check_identity([hypothesis], rho, terms)

    def test_grouped_formula_matches_individual_columns(self):
        # Compare affine row evaluation with the literal ungrouped H formulas.
        # Two copies per order test summation, all-equal values test ties,
        # and zero/one columns test boundary widths.
        for s_order in ORDERS:
            selected = []
            point = [Fraction(0)] * 40
            point[36] = 1
            numeric_selected = []
            for start, order, widths in [
                (0, (0, 1, 2), [Fraction(0), Fraction(1, 4), Fraction(3, 4)]),
                (3, (0, 1, 2), [Fraction(1, 3), Fraction(0), Fraction(0)]),
                (6, s_order, [Fraction(1, 5), Fraction(1, 5), Fraction(2, 5)]),
            ]:
                col, _ = band_column(start, order)
                selected.append(col)
                point[start:start + 3] = widths
                vals = [Fraction(0)] * 3
                for i, p in enumerate(order):
                    vals[p] = sum(widths[:i + 1])
                numeric_selected.append(vals)
            groups, numeric_groups, all_widths = [], [], []
            for i, order in enumerate(ORDERS):
                start = 9 + 4 * i
                col, widths = band_column(start, order, grouped=True)
                groups.append(col)
                all_widths.extend(widths)
                for bands in ([Fraction(1, 10), Fraction(2, 10), Fraction(3, 10), Fraction(4, 10)],
                              [Fraction(0), Fraction(0), Fraction(0), Fraction(1)]):
                    vals = [Fraction(0)] * 3
                    for j, p in enumerate(order):
                        vals[p] = sum(bands[:j + 1])
                    numeric_groups.append(vals)
                    for j, width in enumerate(bands):
                        point[start + j] += width
            q = [Fraction(1, 10), Fraction(1, 5), Fraction(3, 10)]
            point[33:36] = q
            rows = criterion(selected, groups, [unit(33 + p) for p in range(3)], add(*all_widths))
            actual = [sum(a * b for a, b in zip(row, point)) for row in rows]
            cols = numeric_selected + numeric_groups
            expected = list(q)
            for p in range(3):
                expected.append(q[p] - sum(c[p] for c in cols) + len(cols) - 1)
                expected.extend(c[p] - q[p] for c in numeric_selected)
            for a in range(3):
                for b in range(3):
                    if a != b:
                        expected.append(q[b] - q[a] + sum(max(c[a] - c[b], 0) for c in cols))
            for c in range(3):
                a, b = [p for p in range(3) if p != c]
                common = q[c] - q[a] - q[b] + sum(max(min(col[a], col[b]) - col[c], 0) for col in cols)
                expected.extend(common + max(col[a], col[b]) for col in numeric_selected)
            expected.append(sum(q) - sum(min(c) + max(c) for c in cols) + 2 * len(cols) - 2)
            self.assertEqual(actual, expected)

    def test_manuscript_counts_all_six_orders(self):
        for order in ORDERS:
            hypotheses, rows = build_system(order)
            projected, directions, circuits = project(rows)
            self.assertEqual((len(hypotheses), len(rows), directions, circuits, len(projected)),
                             (67, 87, 18, 160, 1893))


class CertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.path = Path(__file__).with_name('split_certificate.json')
        cls.certificate = json.loads(cls.path.read_text(encoding='utf-8'))

    def test_full_certificate(self):
        self.assertEqual(verify(self.path, verbose=False), 11358)

    def test_rejects_modified_or_incomplete_certificate(self):
        for mutation in ('weight', 'projection', 'order'):
            certificate = copy.deepcopy(self.certificate)
            case = certificate['cases']['123']
            if mutation == 'weight':
                case['identities'][0][0][1] += 1
            elif mutation == 'projection':
                case['projected'].pop()
                case['identities'].pop()
            else:
                del certificate['cases']['321']
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'bad.json'
                path.write_text(json.dumps(certificate), encoding='utf-8')
                with self.assertRaises(ValueError):
                    verify(path, verbose=False)


if __name__ == '__main__':
    unittest.main()
