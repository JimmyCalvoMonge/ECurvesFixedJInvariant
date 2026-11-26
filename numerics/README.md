# Numerical verification

Code for the numerical section of "Counting elliptic curves with prescribed j-invariant over number fields".

```sh
python -m pip install cypari==2.5.7    # only dependency; bundles PARI 2.15.4
python fixed_j_numerics.py             # about 15 seconds; writes output/
```

| Table | Statement | What is compared |
|---|---|---|
| 1 | Theorem 1.1 | exact formula against direct enumeration of isomorphism classes (count, ideals, fibres, heights) |
| 1 | Theorem 1.4 | exact count at Y = 10^6 against the main term |
| 2 | Theorem 1.9, Algorithm 1 | least height and number of minimizing classes against direct enumeration |
| 3 | Theorem 1.6(ii) | X^{1/6}(1-p_0) against its limit C_{K,1728}/C_{K,0}, over Q, Q(i) and Q(sqrt -3), up to X = 10^24 |

Everything uses alpha = beta = 1. The direct enumeration uses only the definitions (scaling ideal,
height, isomorphism), never `n_a`, `D_{K,j}`, `Cl^n` or `q_n(K)`; both sides rely on PARI's field
arithmetic. `output/` holds the tables (`tab_*.tex`, pasted into the paper), the data (`*.csv`) and
`report.txt` (46 checks). The exit status is nonzero if a check fails.

As a sanity check, deliberately broken copies of the formulas do fail (of 46 checks): no gcd(n,w_K) 15,
no n^r 3, no |Cl[n]| 8, n_a with ceil(v/2) only 9, with ceil(v/3) only 9; with no coset condition
4 checks fail and the script then stops, since Algorithm 1 finds no curve.
