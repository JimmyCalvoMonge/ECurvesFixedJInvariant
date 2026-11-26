# ECurvesFixedJInvariant

Code for the paper *Counting elliptic curves with prescribed $j$-invariant over number fields*
by Adrian Barquero-Sanchez, Jimmy Calvo-Monge and Ingrid Salazar.

The paper counts the $K$-isomorphism classes of elliptic curves over a number field $K$ with a fixed
$j$-invariant and naive height at most $X$. This repository computes those counts with the paper's
algorithm and checks the main theorems numerically.

## The algorithm

Theorem 1.1 reduces the count to counting ideals:

$$\bigl|\mathcal E_{K,j}(X)\bigr| = q_n(K)\,|\mathrm{Cl}[n]|\cdot\Bigl|\bigl\lbrace\mathfrak b\ \text{an } n\text{-free integral ideal} : \mathrm N\mathfrak b\le X^{1/m}/D_{K,j},\ [\mathfrak n\,\mathfrak b]\in\mathrm{Cl}^n\bigr\rbrace\Bigr|,$$

where $(n,m) = (2,6)$ for $j\ne 0,1728$, $(4,3)$ for $j=1728$ and $(6,2)$ for $j=0$, and
$q_n(K)=\gcd(n,w_K)\,n^r$.

**Algorithm 2** (Section 7.3 of the paper) turns this into a finite computation:

1. Set $n$, $m$, the ideal $\mathfrak n$ and the constant $D_{K,j}$. For $j\ne0,1728$, put
   $a=4(1728-j)/(27j)$; then $\mathfrak n$ has exponents $\max(\lceil v_{\mathfrak p}(a)/2\rceil,\lceil v_{\mathfrak p}(a)/3\rceil)$.
2. Compute the class group $\mathrm{Cl}\cong\bigoplus_i\mathbb Z/e_i\mathbb Z$, the unit rank $r$ and $w_K$.
3. Count the $n$-free ideals $\mathfrak b$ with $\mathrm N\mathfrak b\le Y=X^{1/m}/D_{K,j}$ and $[\mathfrak n\,\mathfrak b]\in\mathrm{Cl}^n$.
4. Multiply by $\gcd(n,w_K)\,n^r\prod_i\gcd(n,e_i)$.

**Algorithm 1** (Section 6) finds the least height of a curve with invariant $j$ and an explicit curve attaining it.

In [`numerics/fixed_j_numerics.py`](numerics/fixed_j_numerics.py):

| Paper | Code |
|---|---|
| Algorithm 2, as written: loop over all ideals of norm at most $Y$ | `Family.count` |
| Algorithm 2, faster variant from the remark after it: build the $n$-free ideals from prime ideals, tracking their class | `free_counts` |
| Algorithm 1 with Remark 6.5: least height and all minimizing curves | `Family.least` |

The comparison $\mathrm N\mathfrak b\le Y$ is done exactly, as $D_{K,j}^m(\mathrm N\mathfrak b)^m\le X$ with rational arithmetic.
Class groups are certified with PARI's `bnfcertify`, so the results do not depend on GRH.

## The verifications

Each formula is compared with a direct enumeration of the curves. The enumeration lists the Weierstrass
models $y^2=x^3+Ax+B$ and computes the scaling ideal, the height and the isomorphism classes from their
definitions. It never uses $\mathfrak n$, $D_{K,j}$, $\mathrm{Cl}^n$ or $q_n(K)$.

The fields are $\mathbb Q$, $\mathbb Q(i)$, $\mathbb Q(\sqrt2)$, $\mathbb Q(\sqrt{-5})$ and $\mathbb Q(\sqrt{-23})$,
with $\mathbb Q(\sqrt{-3})$ and $\mathbb Q(\sqrt{-21})$ added where they are useful. Everything uses $\alpha=\beta=1$.

| Table in the paper | Theorem | What is checked |
|---|---|---|
| 1 | Theorem 1.1 | Algorithm 2 against the direct enumeration, for $j=1728/7$, $1728$ and $0$ |
| 1 | Theorem 1.4 | the exact count at $Y=10^6$ against the predicted leading constant (relative error below $10^{-4}$) |
| 2 | Theorem 1.9 | Algorithm 1: least height and number of minimizing classes, against the enumeration |
| 3 | Theorem 1.6(ii) | over $\mathbb Q$, $\mathbb Q(i)$, $\mathbb Q(\sqrt{-3})$: $X^{1/6}(1-p_0(X))$ approaches $C_{K,1728}/C_{K,0}$, up to $X=10^{24}$ |

The script runs 46 checks and exits with a nonzero status if any fails. To make sure the checks can
actually fail, the formula was broken on purpose in several ways (no $\gcd(n,w_K)$, no $n^r$, no $|\mathrm{Cl}[n]|$,
the wrong exponents in $\mathfrak n$, no class condition); each broken version fails between 3 and 15 checks.
See [`numerics/README.md`](numerics/README.md) for details.

## Running it

```sh
python -m venv venv && source venv/bin/activate
python -m pip install cypari==2.5.7      # the only dependency; it bundles PARI 2.15.4
python numerics/fixed_j_numerics.py      # about 15 seconds
```

It writes to `numerics/output/`:

- `tab_counts.tex`, `tab_least.tex`, `tab_cm.tex`: the three tables of the paper, as LaTeX (not tracked in git)
- `counts.csv`, `least.csv`, `cm.csv`: the same data as CSV
- `report.txt`: every check with its result

## Repository layout

- `numerics/`: the verification script, its README and its output.
- `legacy/`: the first prototype (SageMath code, a notebook and early independent checks). It is kept for reference only; the paper does not depend on it.
