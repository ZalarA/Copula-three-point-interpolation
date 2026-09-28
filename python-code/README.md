# Exact verification of the projection lemma

This is a reconstruction from the displayed formulas in
`TeX/multivariate_copulas_three_points_v2.tex`. The supplementary programs
and certificate mentioned in that manuscript were absent from this checkout;
these files provide a new implementation and a newly generated certificate.

The check proves that (H1)--(H6) and the coordinate-order constraints imply
simultaneous feasibility of (R1)--(R6), (C1)--(C4), and
`0 <= z1 <= z2 <= z3 <= 1`, for each of the six orders of `s`.
It covers arbitrary nonnegative band totals for all remaining columns,
including the weighted version of the selected-coordinate lemma.

## Run the exact verifier

Requires Python 3.10 or newer. From the repository root:

```sh
python python-code/verify_split.py
```

The included `split_certificate.json` is sufficient. Verification uses only
the standard library: no optimization solver, floating-point tolerances, or
third-party packages. An invalid or incomplete certificate causes a nonzero
exit status. A different certificate path can be supplied as the positional
argument.

For each order the verifier reconstructs:

| Quantity | Count |
| --- | ---: |
| Hypothesis and ordered-grid rows | 67 |
| Simultaneous input rows after deduplication | 87 |
| Nonzero coefficient directions in `z` | 18 |
| Positive circuits | 160 |
| Distinct nonzero projected inequalities | 1,893 |

Success ends with:

```text
PASS: all 11358 identities verified; Farkas implies simultaneous feasibility.
```

## What is verified

Every input row is represented as `h + b.z >= 0`. Cofactor determinants
enumerate all positive circuits on two, three, or four nonzero directions
`b`. Exact cancellation and strict positivity are checked. These are all
extreme rays of the nonnegative multiplier cone because there are three
eliminated variables. Every choice of input rows having the required
directions is included, as are all zero-direction rows. Positive gcd
normalization removes duplicates without reversing inequality signs.

For every projected row `rho`, the certificate stores nonnegative rational
weights giving `rho = sum(mu[j] * hypothesis[j])`. The verifier reconstructs
the complete projected row list independently of the certificate and checks
all 37 coefficients in every identity, including the constant. Thus every
Farkas obstruction is excluded and a common feasible `z` exists.

This is a universal implication check, not a numerical feasibility test for
one data set. It does not compute a particular `z`, construct the bivariate
copula, or verify the separate finite-dimensional base-case facet calculations.

## Regenerate the certificate

```sh
python -m pip install -r python-code/requirements.txt
python python-code/generate_split_certificate.py
```

Generation uses SciPy's linear programming solver to find candidate supports
for the multipliers. Each candidate is then solved again by rational Gaussian
elimination and checked exactly. Floating-point solver success alone is never
accepted as a proof. A failed search stops generation; it is not reported as
an exact disproof of the lemma. Once all cases pass, the generator writes the
certificate and invokes the independent verifier. An optional output path
can be supplied to preserve the included certificate.

## Representation

Each row has 40 integer entries, indexed from zero:

- `0:9`: three free ordered-band widths for each of `u`, `v`, `s`;
- `9:33`: four unrestricted nonnegative band totals per order of the other
  columns, with orders listed lexicographically;
- `33:36`: `q1`, `q2`, `q3`;
- `36`: constant term;
- `37:40`: `z1`, `z2`, `z3`.

For an individual selected column the fourth width is one minus the first
three. For grouped columns all four widths are independent nonnegative
totals; their sum is the group weight. Both the hypothesis system and the
conclusion system include the 36 ordered-grid constraints.

The JSON format is `copula-split-v1`, with one case per order of `s`.
Each case stores the sorted hypothesis and projected rows and an identity
for every projected row in that order. Each identity is a sparse list of
`[hypothesis_index, positive_numerator, positive_denominator]` entries.

## Regression checks

```sh
python -m unittest discover -s python-code -p "test_*.py" -v
```

Tests cover circuit degeneracies, repeated coefficient directions, the
constant entry, grouped versus literal coordinate formulas (including ties
and boundary coordinates), manuscript counts, the full certificate, and
rejection of altered weights or missing inequalities/cases.
