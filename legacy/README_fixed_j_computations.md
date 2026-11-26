# Fixed-\(j\) elliptic curves over number fields — computation prototype

**Status: first implementation, 13 October 2026.** The SageMath source code passes Python syntax compilation, but the Sage-specific test suite has **not yet been executed** because SageMath is unavailable in the current environment. An independent, Sage-free finite-count baseline **has been executed successfully**. The results below should be treated as research computations pending Sage execution and code review.

## Files

- `fixed_j_sage.py` — main SageMath/PARI-backed implementation.
- `test_fixed_j_sage.py` — assertions to run in Sage, including uniqueness of generic twist classes.
- `fixed_j_experiments.ipynb` — numerical examples, with **SageMath** kernel.
- `independent_example_tables.py` — a separate Python implementation of ideal-factor enumeration and numerical leading constants for three *specific* fields; requires `mpmath` and `sympy`.
- `verify_fixed_j_counts.py` — the existing independent finite ideal/element-count comparison; requires `sympy`.
- `fixed_j_independent_baseline.csv` — exact and predicted counts computed by the independent script.

## Quick start

Use SageMath 10.x, preferably 10.8 or later. Unzip the bundle and run from its directory:

```sh
sage -python test_fixed_j_sage.py
```

To experiment interactively, choose a SageMath kernel for `fixed_j_experiments.ipynb` or use:

```python
from sage.all import QQ, QuadraticField
from fixed_j_sage import (fixed_j_representatives, fixed_j_count,
                          fixed_j_minimal_curves, fixed_j_leading_constant)
K = QuadraticField(-5, 's')
j = K(QQ(3456)/29)
reps = fixed_j_representatives(K, j, ideal_norm_bound=2)
assert len(reps) == 4
assert fixed_j_count(K, j, ideal_norm_bound=10) == 24
smallest = fixed_j_minimal_curves(K, j)
assert smallest['min_height'] == 256
print([(r['A'],r['B']) for r in reps])
```

For **curve height** \(X\) rather than an ideal norm \(Y\), use `X=...` instead of `ideal_norm_bound=...`. Specify exactly one; the code uses an exact algebraic-real comparison (Sage `AA`) to determine the ideal-norm cutoff.

## Key functions

| Function | Role |
|---|---|
| `fixed_j_data(K,j,alpha=1,beta=1)` | Twist degree, normalizing fractional ideal, exact scale constant |
| `fixed_j_representatives(K,j,X=...)` | One explicit model of every class of bounded intrinsic height |
| `fixed_j_representatives(K,j,ideal_norm_bound=Y)` | Same enumeration indexed by ideal norm |
| `fixed_j_count(K,j,...)` | Exact count without constructing curves |
| `fixed_j_minimal_curves(K,j)` | Least intrinsic height and all minimizing classes |
| `intrinsic_height_direct(K,A,B)` | Independent direct height calculation via scaling-ideal valuations |
| `fixed_j_leading_constant(K,j)` | Numerically evaluate the leading constant in the height-ordered asymptotic |
| `fixed_j_leading_term(K,j,X=...)` | Height-ordered leading prediction |
| `fixed_j_comparison_table(K,j,Y_values)` | Actual counts, predicted main terms, relative errors |

In the generic case \(n=2\), a representative record contains a **parameter** \(t\in K^\times\), its model
\(y^2=x^3+(t^2/a)x+t^3/a\), the exact intrinsic height and the associated squarefree ideal.

For \(j=0\) and \(1728\), the same API uses \(n=6\) and \(4\), respectively. It constructs \(y^2=x^3+z\) and \(y^2=x^3+zx\).

**Caveat:** Models need not have integral coefficients, especially when the ideal class group is nontrivial. They represent the correct \(K\)-isomorphism classes. It is **incorrect** to compute intrinsic heights merely as the unnormalized maximum of their displayed coefficients.

## Algebraic algorithm

1. For generic \(j\), calculate \(a=4(1728-j)/(27j)\) and the fractional ideal
   \(\mathfrak N_a=\prod_{\mathfrak p}\mathfrak p^{\max(\lceil v_{\mathfrak p}(a)/2\rceil,\lceil v_{\mathfrak p}(a)/3\rceil)}\).
2. Enumerate integral \(n\)-free ideals \(\mathfrak b\) by norm (Sage `ideals_of_bdd_norm`).
3. For the generic case, find **all** ideal classes \([\mathfrak c]\) solving
   \([\mathfrak c]^2=[\mathfrak N_a\mathfrak b]^{-1}\); for the exceptional cases solve
   \([\mathfrak c]^n=[\mathfrak b]^{-1}\).
4. Construct a generator of the principal fractional ideal
   \(\mathfrak N_a\mathfrak b\mathfrak c^2\) (or \(\mathfrak b\mathfrak c^n\)).
5. Multiply by **all** representatives of \(\mathcal O_K^\times/(\mathcal O_K^\times)^n\).
6. Form the Weierstrass equations and attach the exact height from
   \(H(E)=(\Lambda_{K,j;\alpha,\beta}N\mathfrak b)^{12/n}\).

All ideal-class roots, including class-group \(n\)-torsion, are enumerated explicitly. The count-only routine multiplies the number of roots by the size of the unit quotient.

## Verified baseline (independent Python computation)

Counts at ideal-norm bound \(Y=10\):

| \(K\) | \(j\) | Exact count | Main term |
|---|---|---:|---:|
| \(\mathbb Q\) | \(6912/31\) | 14 | 12.15854 |
| \(\mathbb Q(\sqrt2)\) | \(6912/31\) | 20 | 17.37248 |
| \(\mathbb Q(\sqrt{-5})\) | \(3456/29\) | 24 | 15.14330 |

The \(\mathbb Q(\sqrt{-5})\) example starts at \(Y=2\), with exactly four classes of minimum intrinsic height \(256\). An integral set of representatives corresponds to parameters \(t=\pm2,\pm4\), giving models with \((A,B)=(2,\pm4),(8,\pm32)\).

The full baseline data at larger values of \(Y\) appears in the CSV file. The independent implementation also compares ideal enumeration with explicit quadratic-element square-class testing for small bounds.

## Important limitations and verification steps

1. **Unexecuted Sage code:** run `sage -python test_fixed_j_sage.py` and investigate any version-specific API differences. This step is required before using computations in a research manuscript.
2. **Exact inputs:** \(\alpha,\beta\) are required to be *positive rationals*; arbitrary real calibration is not implemented. The field \(K\) should be `QQ` or a Sage **absolute** number field with maximal order.
3. **Efficiency:** generating `K.ideals_of_bdd_norm(Y)` may be large. A later implementation could use an Euler-product/ideal-class dynamic program and construct models only when requested.
4. **Numerical constants:** \(\rho_K\), regulators and \(\zeta_K(n)\) are evaluated numerically. Leading terms are **approximations**, even when exact class enumeration is used.
5. **CM totals:** enumeration over all CM \(j\in K\) is a planned second phase; it is **not** included here. The current code supports each such \(j\) separately.
6. **Reproducibility:** for a publication, report Sage/PARI versions, tested proof flags, precision, and independent checks for at least one nontrivial class group.
