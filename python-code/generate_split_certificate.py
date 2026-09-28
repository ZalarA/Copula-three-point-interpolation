"""Generate rational witnesses; verify_split.py independently checks them.

SciPy is used only to find candidate supports. Every accepted witness is
solved again over Fraction and checked exactly before it is saved.
"""

import argparse
from fractions import Fraction
import json
from pathlib import Path

from verify_split import ORDERS, build_system, check_identity, project, verify


def exact_solution(hypotheses, support, target):
    """Solve H_support^T x = target by rational Gaussian elimination."""
    n = len(support)
    matrix = [[Fraction(hypotheses[j][i]) for j in support] + [Fraction(target[i])]
              for i in range(len(target))]
    pivot_row = 0
    pivots = []
    for col in range(n):
        pivot = next((r for r in range(pivot_row, len(matrix)) if matrix[r][col]), None)
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        divisor = matrix[pivot_row][col]
        matrix[pivot_row] = [x / divisor for x in matrix[pivot_row]]
        for r in range(len(matrix)):
            if r != pivot_row and matrix[r][col]:
                multiplier = matrix[r][col]
                matrix[r] = [a - multiplier * b for a, b in zip(matrix[r], matrix[pivot_row])]
        pivots.append(col)
        pivot_row += 1
    if any(not any(row[:n]) and row[n] for row in matrix):
        raise ValueError('Candidate support has no exact solution')
    solution = [Fraction(0)] * n
    for r, col in enumerate(pivots):
        solution[col] = matrix[r][n]
    if any(x < 0 for x in solution):
        raise ValueError('Candidate support yields a negative rational multiplier')
    return [[j, x.numerator, x.denominator] for j, x in zip(support, solution) if x]


def generate(path):
    try:
        import numpy as np
        from scipy.optimize import linprog
    except ImportError as error:
        raise RuntimeError('Generation needs SciPy: python -m pip install -r requirements.txt') from error

    certificate = {'format': 'copula-split-v1', 'cases': {}}
    for order in ORDERS:
        key = ''.join(str(p + 1) for p in order)
        hypotheses, rows = build_system(order)
        projected, directions, circuits = project(rows)
        print(f'{key}: {len(hypotheses)} hypotheses, {len(rows)} input rows, '
              f'{directions} directions, {circuits} circuits, {len(projected)} projections', flush=True)
        matrix = np.array(hypotheses, dtype=float).T
        direct = {row: [[j, 1, 1]] for j, row in enumerate(hypotheses)}
        identities = []
        for i, rho in enumerate(projected):
            if rho in direct:
                terms = direct[rho]
            else:
                result = linprog(np.ones(len(hypotheses)), A_eq=matrix,
                                 b_eq=np.array(rho, dtype=float), bounds=(0, None),
                                 method='highs-ds')
                if not result.success:
                    raise RuntimeError(f'{key}, row {i}: certificate search failed: {result.message}. '
                                       'This is not an exact infeasibility certificate.')
                support = [j for j, x in enumerate(result.x) if x > 1e-9]
                try:
                    terms = exact_solution(hypotheses, support, rho)
                except ValueError as error:
                    raise RuntimeError(f'{key}, row {i}: exact reconstruction failed: {error}') from error
            check_identity(hypotheses, rho, terms)
            identities.append(terms)
            if (i + 1) % 500 == 0:
                print(f'  {i + 1}/{len(projected)} exact witnesses constructed', flush=True)
        certificate['cases'][key] = {
            'hypotheses': hypotheses, 'projected': projected, 'identities': identities,
        }
    # Write only after all six orders have passed exact identity checks.
    path.write_text(json.dumps(certificate, separators=(',', ':')) + '\n', encoding='utf-8')
    verify(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', nargs='?', type=Path,
                        default=Path(__file__).with_name('split_certificate.json'))
    args = parser.parse_args()
    generate(args.output)


if __name__ == '__main__':
    main()
