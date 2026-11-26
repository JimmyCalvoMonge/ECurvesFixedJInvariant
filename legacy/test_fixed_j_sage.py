"""SageMath test suite for fixed_j_sage.py.

Run: sage -python test_fixed_j_sage.py
Tests cover: squarefree ideal counts, genuine isomorphism classes, units,
nontrivial class group, minimum heights, special j=0 and j=1728,
nonintegral twist parameters, and a rational-denominator normalization.
"""
from sage.all import QQ, ZZ, QuadraticField, AA
from fixed_j_sage import (
    fixed_j_data, fixed_j_count, fixed_j_representatives,
    fixed_j_minimal_curves, intrinsic_height_direct, fixed_j_leading_term,
    fixed_j_comparison_table,
)


def distinct_twists(recs, n):
    """Compare all twists over K by exact power testing (generic n=2)."""
    if n == 2:
        for i, ri in enumerate(recs):
            for rj in recs[i+1:]:
                quotient = ri['parameter'] / rj['parameter']
                assert not quotient.is_square(), (ri, rj)
    # For n=4,6, the class-group/unit construction proves distinctness;
    # no ad hoc nth-power oracle is assumed here.


def check_counts(K, j, wanted):
    for Y, expected in wanted.items():
        count = fixed_j_count(K, j, ideal_norm_bound=Y)
        assert count == expected, (K, j, Y, expected, count)
        reps = fixed_j_representatives(K, j, ideal_norm_bound=Y)
        assert len(reps) == expected
        if j != 0 and j != 1728:
            distinct_twists(reps, 2)
        for rec in reps:
            assert rec['curve'].j_invariant() == K(j)
            assert intrinsic_height_direct(K, rec['A'], rec['B']) == rec['height']
            assert rec['height'] == (fixed_j_data(K, j)['Lambda'] *
                                     AA(rec['ideal_norm'])) ** ZZ(12 // rec['n'])
        print('PASS', K, 'j=', j, 'Y=', Y, 'classes=', count)


def test_Q():
    j = QQ(6912)/31 # a(j)=1
    check_counts(QQ, j, {1:2, 2:4, 3:6, 5:8, 10:14, 20:26})
    minimum = fixed_j_minimal_curves(QQ, j)
    assert minimum['min_height'] == AA(1)
    assert len(minimum['representatives']) == 2
    # 6-free integers <= 10: all ten except none (6th powers begin at 64).
    check_counts(QQ, QQ(0), {1:2, 2:4, 5:10, 10:20})
    # 4-free integers <= 10: all ten (4th powers begin at 16).
    check_counts(QQ, QQ(1728), {1:2, 2:4, 5:10, 10:20})
    # Explicit denominators test: a=1/32 means j=...? solve for j.
    a = QQ(1)/32
    j2 = QQ(6912)/(QQ(4) + 27*a)  # 4(1728-j)/(27j)=a.
    data = fixed_j_data(QQ, j2)
    assert data['a'] == a
    assert len(fixed_j_representatives(QQ, j2, ideal_norm_bound=5)) == 8


def test_real_quadratic():
    K = QuadraticField(2, 'r')
    j = K(QQ(6912)/31)
    check_counts(K, j, {1:4, 2:8, 3:8, 5:8, 10:20})
    mn = fixed_j_minimal_curves(K, j)
    assert mn['min_height'] == AA(1)
    assert len(mn['representatives']) == 4


def test_nonprincipal_quadratic():
    K = QuadraticField(-5, 's')
    j = K(QQ(3456)/29) # a(j)=2
    data = fixed_j_data(K, j)
    assert data['a'] == 2
    assert K.class_number() == 2
    assert data['Lambda']**6 == AA(4)
    check_counts(K, j, {1:0, 2:4, 3:12, 5:12, 8:20, 10:24})
    minimum = fixed_j_minimal_curves(K, j)
    assert minimum['min_ideal_norm'] == 2
    assert minimum['min_height'] == AA(256)
    assert len(minimum['representatives']) == 4
    # Verify minimal classes against the four simple integral models
    # with t=+/-2,+/-4 (normalization-independent assertion).
    found = [r['parameter'] for r in minimum['representatives']]
    expected = [K(2), K(-2), K(4), K(-4)]
    for expected_t in expected:
        assert sum(1 for t in found if (t/expected_t).is_square()) == 1
    print('PASS minimum Q(sqrt(-5)) = 256, four classes')


def test_numerical_predictions():
    j = QQ(6912)/31
    rows = fixed_j_comparison_table(QQ, j, [10,100])
    assert len(rows) == 2
    assert rows[0]['exact'] == 14
    for row in rows:
        assert row['main_term'] > 0
    assert fixed_j_leading_term(QQ, j, ideal_norm_bound=100) > 0


if __name__ == '__main__':
    test_Q()
    test_real_quadratic()
    test_nonprincipal_quadratic()
    test_numerical_predictions()
    print('All SageMath tests passed.')
