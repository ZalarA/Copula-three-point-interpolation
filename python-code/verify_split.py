"""Exact Farkas projection verifier for the selected-coordinate merging lemma.

Reconstructs H1--H6, R1--R6 and C1--C4 from the manuscript, enumerates
every positive circuit, and verifies nonnegative rational certificates.
Verification uses only the Python standard library. See README.md.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from functools import reduce
from itertools import combinations, permutations, product
import json
from math import gcd
from pathlib import Path


ORDERS = tuple(permutations(range(3)))
SIZE = 40  # 36 data variables, constant, three eliminated variables
ZERO = (0,) * SIZE


def unit(i):
    return tuple(int(j == i) for j in range(SIZE))


ONE = unit(36)


def add(*rows):
    return tuple(map(sum, zip(*rows))) if rows else ZERO


def scale(n, row):
    return tuple(n * x for x in row)


def sub(a, b):
    return add(a, scale(-1, b))


def primitive(row):
    """Normalize by a positive factor only; never reverse an inequality."""
    divisor = reduce(gcd, row, 0)
    return tuple(x // divisor for x in row) if divisor else row


def unique(rows):
    return sorted({primitive(row) for row in rows if any(row)})


class Column:
    def __init__(self, values, order):
        self.values = values
        self.order = order
        self.rank = {p: i for i, p in enumerate(order)}

    def __getitem__(self, p):
        return self.values[p]

    def maximum(self, a, b):
        return self[max((a, b), key=self.rank.get)]

    def difference(self, a, b):
        return sub(self[a], self[b]) if self.rank[a] > self.rank[b] else ZERO

    def leakage(self, a, b, c):
        return self.difference(min((a, b), key=self.rank.get), c)


def band_column(start, order, grouped=False):
    """First three widths parameterize a unit column; groups use four totals."""
    widths = [unit(start + i) for i in range(3)]
    widths.append(unit(start + 3) if grouped else sub(ONE, add(*widths)))
    values = [ZERO] * 3
    for i, p in enumerate(order):
        values[p] = add(*widths[:i + 1])
    return Column(values, order), widths


def criterion(selected, groups, values, group_weight):
    """Slacks of H/R/C in the order 1,...,6, splitting outer minima.

    Each group column stores sums, not averages. group_weight is the
    sum of the four band totals over all six groups (the weighted r).
    """
    rows = list(values)  # H1 / R1 / C1 nonnegativity
    all_columns = selected + groups
    dimension = add(scale(len(selected), ONE), group_weight)
    for p in range(3):
        rows.append(add(values[p], scale(-1, add(*(c[p] for c in all_columns))),
                        dimension, scale(-1, ONE)))
        rows.extend(sub(c[p], values[p]) for c in selected)
    for a, b in permutations(range(3), 2):
        rows.append(add(values[b], scale(-1, values[a]),
                        *(c.difference(a, b) for c in all_columns)))
    for c in range(3):
        a, b = [p for p in range(3) if p != c]
        common = add(values[c], scale(-1, values[a]), scale(-1, values[b]),
                     *(col.leakage(a, b, c) for col in all_columns))
        rows.extend(add(common, col.maximum(a, b)) for col in selected)
    extremes = add(*(add(c[c.order[0]], c[c.order[2]]) for c in all_columns))
    rows.append(add(*values, scale(-1, extremes), scale(2, dimension), scale(-2, ONE)))
    return rows


def build_system(order):
    """Return hypothesis rows and simultaneous conclusion rows (slack >= 0)."""
    u, uw = band_column(0, (0, 1, 2))
    v, vw = band_column(3, (0, 1, 2))
    s, sw = band_column(6, order)
    groups, grid = [], uw + vw + sw
    grouped_widths = []
    for i, group_order in enumerate(ORDERS):
        col, widths = band_column(9 + 4 * i, group_order, grouped=True)
        groups.append(col)
        grid.extend(widths)
        grouped_widths.extend(widths)
    weight = add(*grouped_widths)
    q = [unit(33 + p) for p in range(3)]
    z = Column([unit(37 + p) for p in range(3)], (0, 1, 2))
    hypotheses = unique(grid + criterion([u, v, s], groups, q, weight))
    conclusions = unique(
        grid + criterion([z, s], groups, q, weight)
        + criterion([u, v], [], z.values, ZERO)
        + [z[0], sub(z[1], z[0]), sub(z[2], z[1]), sub(ONE, z[2])]
    )
    return [row[:37] for row in hypotheses], conclusions


def determinant(matrix):
    if not matrix:
        return 1
    return sum((-1) ** j * x * determinant([row[:j] + row[j + 1:] for row in matrix[1:]])
               for j, x in enumerate(matrix[0]))


def positive_circuits(directions):
    """All extreme multiplier rays on nonzero directions in R^3.

    A nonzero cofactor vector proves rank >= t-1; checking its cancellation
    against all three rows proves rank <= t-1. Strict positivity excludes
    smaller supports. This avoids any numerical rank decisions.
    """
    for t in range(2, 5):
        for support in combinations(directions, t):
            matrix = [[d[i] for d in support] for i in range(3)]
            for indices in combinations(range(3), t - 1):
                minor = [matrix[i] for i in indices]
                vector = [(-1) ** j * determinant([r[:j] + r[j + 1:] for r in minor])
                          for j in range(t)]
                if not any(vector):
                    continue
                if all(x < 0 for x in vector):
                    vector = [-x for x in vector]
                if (all(x > 0 for x in vector)
                        and all(sum(a * b for a, b in zip(row, vector)) == 0 for row in matrix)):
                    yield support, primitive(tuple(vector))
                break


def project(rows):
    """Eliminate z exactly, including every row sharing a direction."""
    by_direction = defaultdict(list)
    projected = set()
    for row in rows:
        direction = row[37:]
        if any(direction):
            by_direction[direction].append(row[:37])
        elif any(row[:37]):
            projected.add(primitive(row[:37]))
    circuits = list(positive_circuits(sorted(by_direction)))
    for support, weights in circuits:
        for selected in product(*(by_direction[d] for d in support)):
            rho = tuple(sum(w * row[i] for w, row in zip(weights, selected)) for i in range(37))
            if any(rho):
                projected.add(primitive(rho))
    return sorted(projected), len(by_direction), len(circuits)


def check_identity(hypotheses, rho, terms):
    """Check the entire rational identity, including its constant entry."""
    total = [Fraction(0)] * 37
    seen = set()
    for index, numerator, denominator in terms:
        if any(type(x) is not int for x in (index, numerator, denominator)):
            raise ValueError('Certificate terms must contain integers')
        if index in seen or not 0 <= index < len(hypotheses):
            raise ValueError('Repeated or invalid hypothesis index')
        seen.add(index)
        if numerator <= 0 or denominator <= 0:
            raise ValueError('Certificate weights must be strictly positive')
        weight = Fraction(numerator, denominator)
        for i, coefficient in enumerate(hypotheses[index]):
            total[i] += weight * coefficient
    if total != list(rho):
        raise ValueError('Rational identity does not match the projected row')


def verify(path, verbose=True):
    certificate = json.loads(Path(path).read_text(encoding='utf-8'))
    if certificate.get('format') != 'copula-split-v1':
        raise ValueError('Unknown certificate format')
    expected_orders = {''.join(str(p + 1) for p in order) for order in ORDERS}
    if set(certificate['cases']) != expected_orders:
        raise ValueError('Certificate must cover exactly all six orders of s')
    total = 0
    for order in ORDERS:
        key = ''.join(str(p + 1) for p in order)
        hypotheses, rows = build_system(order)
        projected, directions, circuits = project(rows)
        case = certificate['cases'][key]
        # Comparing complete row lists prevents omitted, reordered, or foreign
        # inequalities from silently changing which theorem is verified.
        if case['hypotheses'] != [list(row) for row in hypotheses]:
            raise ValueError(f'{key}: hypothesis rows differ')
        if case['projected'] != [list(row) for row in projected]:
            raise ValueError(f'{key}: projected row list is incomplete or differs')
        if len(case['identities']) != len(projected):
            raise ValueError(f'{key}: missing or extra identities')
        for i, (rho, terms) in enumerate(zip(projected, case['identities'])):
            try:
                check_identity(hypotheses, rho, terms)
            except (ValueError, TypeError) as error:
                raise ValueError(f'{key}, projected row {i}: {error}') from error
        total += len(projected)
        if verbose:
            print(f'{key}: {len(hypotheses)} hypotheses, {len(rows)} input rows, '
                  f'{directions} directions, {circuits} circuits, '
                  f'{len(projected)} exact identities OK', flush=True)
    if verbose:
        print(f'PASS: all {total} identities verified; Farkas implies simultaneous feasibility.')
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', nargs='?', type=Path,
                        default=Path(__file__).with_name('split_certificate.json'))
    args = parser.parse_args()
    try:
        verify(args.certificate)
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(1, f'FAIL: {error}\n')


if __name__ == '__main__':
    main()
