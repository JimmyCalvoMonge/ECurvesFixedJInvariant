"""Sage-independent arithmetic baseline for the three proposed example fields.

Only handles the fixed examples in the research manuscript, not arbitrary K.
Uses (1) elementary integer squarefreeness and (2) independent ideal Euler
factor enumerators from verify_fixed_j_counts.py; analytic terms use mpmath.
"""
import csv
from math import sqrt, log, pi
from pathlib import Path
import mpmath as mp
from sympy import kronecker_symbol
from verify_fixed_j_counts import ideal_counts

mp.mp.dps = 50
DATA = Path(__file__).resolve().parent
BOUNDS = [1, 2, 3, 5, 10, 30, 100, 300, 1000]

def ideal_zeta(s, disc):
    conductor = abs(disc)
    chi = [int(kronecker_symbol(disc, n)) for n in range(conductor)]
    return mp.zeta(s) * mp.dirichlet(s, chi)

def main():
    rho = {
        'Q': mp.mpf(1),
        'Qsqrt2': 2*mp.log(1+mp.sqrt(2))/mp.sqrt(8),
        'Qsqrt-5': mp.pi/mp.sqrt(5),
    }
    zeta2 = {
        'Q': mp.zeta(2),
        'Qsqrt2': ideal_zeta(2, 8),
        'Qsqrt-5': ideal_zeta(2, -20),
    }
    q = {'Q': 2, 'Qsqrt2': 4, 'Qsqrt-5': 2}
    rows = []
    for Y in BOUNDS:
        # elementary squarefree positive integers
        def squarefree(n):
            k=2
            while k*k<=n:
                if n%(k*k)==0:return False
                k+=1
            return True
        exactQ = 2*sum(squarefree(n) for n in range(1,Y+1))
        d2 = ideal_counts(Y, 2)
        exact2 = 4*sum(number for (norm, cl),number in d2.items())
        dm5 = ideal_counts(Y, -5)
        exactm5 = 4*sum(number for (norm,cl),number in dm5.items() if cl==1)
        for field, j, exact in [('Q', '6912/31', exactQ),
                                 ('Qsqrt2', '6912/31', exact2),
                                 ('Qsqrt-5', '3456/29', exactm5)]:
            pred = q[field] * rho[field] * Y / zeta2[field]
            rows.append({
                'field': field, 'j': j, 'ideal_norm_bound_Y': Y,
                'height_bound_X': Y**6 if field!='Qsqrt-5' else 4*Y**6,
                'exact_count':exact, 'main_term': mp.nstr(pred, 14),
                'absolute_error':mp.nstr(exact-pred, 14),
                'relative_error':mp.nstr((exact-pred)/pred, 14),
            })
    dest = DATA/'fixed_j_independent_baseline.csv'
    with dest.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys())
        w.writeheader();w.writerows(rows)
    print(dest)
    for field in rho:
        print('\n',field,'leading ideal-norm coefficient:',mp.nstr(q[field]*rho[field]/zeta2[field],16))
        for row in rows:
            if row['field']==field:
                print('Y=',row['ideal_norm_bound_Y'], 'exact=',row['exact_count'],
                      'main=',row['main_term'],'relative=',row['relative_error'])

if __name__=='__main__':main()
