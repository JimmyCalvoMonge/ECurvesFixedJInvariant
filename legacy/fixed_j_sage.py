"""Exact ideal-theoretic enumeration of fixed-j elliptic curves over number fields.

Research prototype for SageMath >= 10.x; execute with ``sage -python``.

Mathematical scope
------------------
* K = QQ or an absolute Sage number field with maximal order.
* j in K and positive *rational* height calibration alpha, beta.
* Elliptic curves are counted up to K-isomorphism, using the intrinsic
  (4,6)-weighted scaling-ideal height.
* Exact bounds use Sage's algebraic real field AA, not floating-point
  rounding, including mixed archimedean embeddings.
* Requires PARI-backed class groups, principal ideal tests and unit groups.
* Enumeration lists all integral ideals of norm <=Y; intended for small and
  medium bounds. The count-only path never constructs elliptic curves.

The local identity and ideal-class fibre theorem are proved in the accompanying
working manuscript. No reliance on a representative twist's raw coefficients.

IMPORTANT: This is a research implementation. Sage execution has not yet been
possible in the author's present container; see test_fixed_j_sage.py.
"""

from itertools import product
from math import gcd
from sage.all import (
    AA, QQ, QQbar, ZZ, RealField,
    EllipticCurve, pi, zeta, prod,
)


def _is_Q(K):
    return K == QQ


def _positive_root(x, n):
    """Positive real algebraic n-th root, exactly."""
    x = AA(x)
    if x < 0 or n <= 0:
        raise ValueError('Expected nonnegative x and n>0')
    return x ** (QQ(1) / ZZ(n))


def _validation(K, j, alpha, beta):
    alpha, beta = QQ(alpha), QQ(beta)
    if alpha <= 0 or beta <= 0:
        raise ValueError('alpha and beta must be positive rationals')
    j = K(j)
    return j, alpha, beta


def _na_rational(a):
    """Positive rational generator of the normalization ideal in Q."""
    a = QQ(a)
    out = QQ(1)
    for p, v in abs(a).factor():
        # A factorization of a positive rational returns prime powers.
        m = max((QQ(v) / 2).ceil(), (QQ(v) / 3).ceil())
        out *= QQ(p) ** m
    return out


def _na_number_field(K, a):
    """Fractional ideal N_a; works for v_p(a) negative as well as positive."""
    out = K.ideal(1)
    for prime, v in K.ideal(a).factor():
        m = max((QQ(v) / 2).ceil(), (QQ(v) / 3).ceil())
        out *= prime ** ZZ(m)
    return out


def _arch_constant(K, a, alpha, beta):
    """M = prod_{all embeddings sigma} max(alpha/|sigma(a)|^3,
       beta/|sigma(a)|^2), an exact positive algebraic real.

    Conjugate complex embeddings are BOTH present. This equals the product
    using one place with complex-place absolute values squared.
    """
    if _is_Q(K):
        x = AA(abs(QQ(a)))
        return max(AA(alpha) / x**3, AA(beta) / x**2)
    vals = [abs(embedding(a)) for embedding in K.embeddings(QQbar)]
    if len(vals) != K.degree():
        raise ArithmeticError('Failed to obtain all absolute embeddings')
    return prod([max(AA(alpha) / AA(v)**3, AA(beta) / AA(v)**2)
                 for v in vals], AA(1))


def fixed_j_data(K, j, alpha=1, beta=1):
    """Return twist degree n, height exponent m=12/n, Lambda, N_a.

    The exact height is (Lambda * Norm(b))^m, b an n-free ideal.
    Lambda is in AA, exact for rational alpha,beta and algebraic j.
    """
    j, alpha, beta = _validation(K, j, alpha, beta)
    d = ZZ(K.degree())
    if j == 0:
        return dict(K=K, j=j, alpha=alpha, beta=beta, n=6, m=2,
                    a=None, Na=None, Lambda=AA(beta)**(QQ(d)/2))
    if j == 1728:
        return dict(K=K, j=j, alpha=alpha, beta=beta, n=4, m=3,
                    a=None, Na=None, Lambda=AA(alpha)**(QQ(d)/3))
    a = K(4) * (K(1728) - j) / (K(27) * j)
    Na = _na_rational(a) if _is_Q(K) else _na_number_field(K, a)
    norm_Na = AA(Na) if _is_Q(K) else AA(Na.norm())
    M = _arch_constant(K, a, alpha, beta)
    Lambda = norm_Na * _positive_root(M, 6)
    return dict(K=K, j=j, alpha=alpha, beta=beta, n=2, m=6,
                a=a, Na=Na, Lambda=Lambda)


def _max_ideal_norm(data, X=None, ideal_norm_bound=None):
    if (X is None) == (ideal_norm_bound is None):
        raise ValueError('Specify exactly one of X or ideal_norm_bound')
    if ideal_norm_bound is not None:
        Y = ZZ(ideal_norm_bound)
        if Y < 0:
            raise ValueError('ideal_norm_bound must be nonnegative')
        return Y
    X = AA(X)
    if X < 0:
        raise ValueError('Height bound X must be nonnegative')
    # (Lambda * Norm(b))^m <= X iff Norm(b) <= X^(1/m)/Lambda.
    Y = _positive_root(X, data['m']) / data['Lambda']
    return ZZ(Y.floor())


def _norm_ideals(K, limit):
    """Yield (norm, ideal) with deterministic norm ordering."""
    if _is_Q(K):
        for n in range(1, ZZ(limit) + 1):
            yield ZZ(n), ZZ(n)  # integer n represents (n) in Z
        return
    by_norm = K.ideals_of_bdd_norm(ZZ(limit))
    for n in range(1, ZZ(limit) + 1):
        for I in by_norm.get(n, []):
            yield ZZ(n), I


def _n_free(K, I, n):
    if _is_Q(K):
        return all(exponent < n for p, exponent in ZZ(I).factor())
    return all(exponent < n for p, exponent in I.factor())


def _root_index(K, n):
    """Index all ideal classes c by c^n; includes every torsion root."""
    if _is_Q(K):
        return None, None
    G = K.class_group(proof=True)
    roots = {}
    for c in G:
        roots.setdefault(c**n, []).append(c)
    return G, roots


def _eligible_roots(K, I, data, G, roots):
    if _is_Q(K):
        return [None]
    if data['n'] == 2:
        target = (G(data['Na'] * I))**(-1)
    else:
        target = (G(I))**(-1)
    return roots.get(target, [])


def _unit_quotient_size(K, n):
    """Size of O_K^*/(O_K^*)^n without enumerating unit representatives."""
    if _is_Q(K):
        return gcd(int(n), 2)
    U = K.unit_group(proof=True)
    return ZZ(gcd(int(n), int(U.zeta_order()))) * ZZ(n)**ZZ(U.rank())


def _unit_representatives(K, n):
    """A complete system for O_K^*/(O_K^*)^n."""
    if _is_Q(K):
        return [QQ(1), QQ(-1)] if n % 2 == 0 else [QQ(1)]
    U = K.unit_group(proof=True)
    basis = list(U.gens_values())
    w = ZZ(U.zeta_order())
    torsion_count = ZZ(gcd(int(n), int(w)))
    root = basis[0]
    free_basis = basis[1:]
    output = []
    for e0 in range(torsion_count):
        for free_exp in product(range(n), repeat=len(free_basis)):
            unit = root**e0
            for g, e in zip(free_basis, free_exp):
                unit *= g**e
            output.append(K(unit))
    assert len(output) == torsion_count * n**len(free_basis)
    return output


def _twist_generator(K, I, data, root_class):
    """Choose a t or z generating N_a*b*c^n (or b*c^n)."""
    if _is_Q(K):
        return K(data['Na'] * I) if data['n'] == 2 else K(I)
    c = root_class.ideal()  # fixed fractional ideal in the root class
    F = ((data['Na'] * I) if data['n'] == 2 else I) * c**data['n']
    if not F.is_principal(proof=True):
        raise ArithmeticError('Class root predicted a nonprincipal ideal')
    gen = K(F.gens_reduced(proof=True)[0])
    if K.ideal(gen) != F:
        raise ArithmeticError('Principal ideal generator failed exact check')
    return gen


def _coefficients(data, parameter):
    K = data['K']
    parameter = K(parameter)
    if data['n'] == 2:
        return K(parameter**2 / data['a']), K(parameter**3 / data['a'])
    if data['n'] == 4:
        return parameter, K(0)
    return K(0), parameter


def fixed_j_count(K, j, X=None, ideal_norm_bound=None,
                  alpha=1, beta=1, return_details=False):
    """Count isomorphism classes EXACTLY, without constructing models.

    ``ideal_norm_bound=Y`` counts classes with norm(b)<=Y, useful for tables.
    ``X=X`` uses the genuine intrinsic-height cutoff (exact algebraic tests).
    """
    data = fixed_j_data(K, j, alpha, beta)
    Y = _max_ideal_norm(data, X, ideal_norm_bound)
    n = data['n']
    G, roots = _root_index(K, n)
    q = _unit_quotient_size(K, n)
    total, eligible_ideals = ZZ(0), ZZ(0)
    for norm, I in _norm_ideals(K, Y):
        if not _n_free(K, I, n):
            continue
        num_roots = len(_eligible_roots(K, I, data, G, roots))
        if num_roots:
            eligible_ideals += 1
            total += q * num_roots
    if return_details:
        return dict(count=total, eligible_ideals=eligible_ideals,
                    Y=Y, unit_multiplicity=q, n=n,
                    class_torsion=1 if G is None else len(next(iter(roots.values()))))
    return total


def fixed_j_representatives(K, j, X=None, ideal_norm_bound=None,
                            alpha=1, beta=1, verify=True):
    """List ONE explicit short Weierstrass model of EVERY class at the bound.

    Each record includes K-curve, (A,B), parameter (t or z), n-free ideal,
    ideal norm and *exact intrinsic height* in AA.

    Representatives may have nonintegral coefficients; they are genuine
    K-models. Integrality is not needed for intrinsic-height counting.
    """
    data = fixed_j_data(K, j, alpha, beta)
    Y = _max_ideal_norm(data, X, ideal_norm_bound)
    n = data['n']
    G, roots = _root_index(K, n)
    units = _unit_representatives(K, n)
    out = []
    for norm, I in _norm_ideals(K, Y):
        if not _n_free(K, I, n):
            continue
        candidates = _eligible_roots(K, I, data, G, roots)
        for c in candidates:
            generator = _twist_generator(K, I, data, c)
            for unit in units:
                t = generator * unit
                A, B = _coefficients(data, t)
                height = (data['Lambda'] * AA(norm))**data['m']
                if X is not None and height > AA(X):
                    # Guard against unexpected cutoff errors.
                    continue
                E = EllipticCurve(K, [K(0), K(0), K(0), A, B])
                if verify and E.j_invariant() != data['j']:
                    raise ArithmeticError('j-invariant verification failed')
                out.append(dict(curve=E, A=A, B=B, parameter=t,
                                n_free_ideal=I, ideal_norm=norm,
                                root_class=c, unit=unit,
                                height=height, n=n))
    if verify:
        count = fixed_j_count(K, j, ideal_norm_bound=Y,
                              alpha=alpha, beta=beta)
        assert len(out) == count, (len(out), count)
    return out



def intrinsic_height_direct(K, A, B, alpha=1, beta=1):
    """Independently calculate the intrinsic (4,6) height from coefficients.

    This routine intentionally does NOT use the squarefree-ideal formula.
    Instead it directly evaluates the scaling ideal
        I_(4,6)(A,B) = prod_P P^min(floor(v_P(A)/4),floor(v_P(B)/6))
    (with valuation of zero interpreted as +infinity) and the archimedean
    maximum. This is useful for mathematical verification of the reduction.
    Rational alpha,beta and an absolute number field K are required.
    """
    alpha, beta = QQ(alpha), QQ(beta)
    if alpha <= 0 or beta <= 0:
        raise ValueError('alpha,beta must be positive')
    A, B = K(A), K(B)
    if A == 0 and B == 0:
        raise ValueError('Singular coefficient pair (0,0)')
    if _is_Q(K):
        efA = {} if A == 0 else dict(abs(QQ(A)).factor())
        efB = {} if B == 0 else dict(abs(QQ(B)).factor())
        scale = QQ(1)
        for prime in set(efA) | set(efB):
            if A == 0:
                exponent = (QQ(efB.get(prime, 0)) / 6).floor()
            elif B == 0:
                exponent = (QQ(efA.get(prime, 0)) / 4).floor()
            else:
                exponent = min((QQ(efA.get(prime, 0)) / 4).floor(),
                               (QQ(efB.get(prime, 0)) / 6).floor())
            scale *= QQ(prime)**exponent
        numerator = max(AA(alpha)*AA(abs(A))**3,
                        AA(beta)*AA(abs(B))**2)
        return numerator / AA(scale)**12

    efA = {} if A == 0 else dict(K.ideal(A).factor())
    efB = {} if B == 0 else dict(K.ideal(B).factor())
    scale = K.ideal(1)
    for prime in set(efA) | set(efB):
        if A == 0:
            exponent = (QQ(efB.get(prime, 0)) / 6).floor()
        elif B == 0:
            exponent = (QQ(efA.get(prime, 0)) / 4).floor()
        else:
            exponent = min((QQ(efA.get(prime, 0)) / 4).floor(),
                           (QQ(efB.get(prime, 0)) / 6).floor())
        scale *= prime**ZZ(exponent)
    arch = AA(1)
    for embedding in K.embeddings(QQbar):
        vA, vB = AA(abs(embedding(A))), AA(abs(embedding(B)))
        arch *= max(AA(alpha)*vA**3, AA(beta)*vB**2)
    return arch / AA(scale.norm())**12

def fixed_j_minimal_curves(K, j, alpha=1, beta=1,
                            initial_norm_bound=1, max_norm_bound=65536):
    """Find least intrinsic height and every class realizing it.

    Searches increasing ideal norm bounds; stops at the first eligible norm.
    May raise if a user-imposed safety cap is reached.
    """
    b = max(ZZ(initial_norm_bound), ZZ(1))
    while b <= max_norm_bound:
        records = fixed_j_representatives(K, j, ideal_norm_bound=b,
                                         alpha=alpha, beta=beta)
        if records:
            min_norm = min(rec['ideal_norm'] for rec in records)
            minimum = [r for r in records if r['ideal_norm'] == min_norm]
            return dict(min_height=minimum[0]['height'],
                        min_ideal_norm=min_norm, representatives=minimum)
        b *= 2
    raise RuntimeError('Increase max_norm_bound to locate minimal classes')


def fixed_j_leading_constant(K, j, alpha=1, beta=1, prec=150):
    """Numerical leading C with N_{K,j}(X) ~ C X^(n/12)."""
    data = fixed_j_data(K, j, alpha, beta)
    R = RealField(ZZ(prec))
    n = data['n']
    if _is_Q(K):
        rho, zeta_n = R(1), R(zeta(R(n)))
        w, rank = 2, 0
    else:
        r1, r2 = K.signature()
        rank = r1+r2-1
        w = K.unit_group(proof=True).zeta_order()
        h = K.class_number(proof=True)
        regulator = R(K.regulator())
        rho = (R(2)**r1 * (2*R(pi))**r2 * h * regulator /
               (w * R(abs(K.discriminant())).sqrt()))
        zeta_n = R(K.zeta_function(prec=prec)(n).real())
    q = n**rank * gcd(n, int(w))
    return R(q) * rho / (zeta_n * R(data['Lambda']))


def fixed_j_leading_term(K, j, X=None, ideal_norm_bound=None,
                         alpha=1, beta=1, prec=150):
    """Numerical main term for height X or ideal-norm parameter Y."""
    data = fixed_j_data(K, j, alpha, beta)
    C = fixed_j_leading_constant(K, j, alpha, beta, prec)
    R = RealField(ZZ(prec))
    if (X is None) == (ideal_norm_bound is None):
        raise ValueError('Specify exactly one of X or ideal_norm_bound')
    if X is not None:
        return C * R(AA(X))**(R(data['n'])/12)
    # The Lambda factor cancels on substituting X=(Lambda*Y)^(12/n).
    return C * R(data['Lambda']) * R(ideal_norm_bound)


def fixed_j_comparison_table(K, j, Y_values, alpha=1, beta=1, prec=150):
    """Rows for a numerical comparison indexed by ideal norm bound Y."""
    rows = []
    data = fixed_j_data(K, j, alpha, beta)
    for Y in Y_values:
        exact = fixed_j_count(K, j, ideal_norm_bound=Y,
                              alpha=alpha, beta=beta)
        prediction = fixed_j_leading_term(K, j, ideal_norm_bound=Y,
                                           alpha=alpha, beta=beta, prec=prec)
        height = (data['Lambda'] * AA(Y))**data['m']
        rows.append(dict(Y=ZZ(Y), X=height, exact=exact,
                         main_term=prediction,
                         error=RealField(prec)(exact) - prediction,
                         relative_error=(RealField(prec)(exact) - prediction)
                                        / prediction))
    return rows
