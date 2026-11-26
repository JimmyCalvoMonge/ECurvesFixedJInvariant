#!/usr/bin/env python3
"""Numerical verification for "Counting elliptic curves with prescribed
j-invariant over number fields".

    python -m pip install cypari==2.5.7     # only dependency; bundles PARI
    python fixed_j_numerics.py

Writes the three tables of the paper to ./output (tab_*.tex and *.csv) together
with report.txt; the exit status is nonzero if a check fails.  Here alpha=beta=1.

  Table 1  Theorem 1.1 (exact counts) and Theorem 1.4 (leading constants)
  Table 2  Theorem 1.9 (least heights)
  Table 3  Theorem 1.6 (CM curves over Q, Q(i), Q(sqrt-3))

The formulas are compared with an enumeration of isomorphism classes (Part 2)
that uses only the definitions: pairs (A,B), scaling ideal, height, and the
relation (A',B')=(u^4 A,u^6 B).
"""

import csv
import itertools
import math
import os
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache

import cypari
from cypari import pari

pari.allocatemem(2 * 10**9, silent=True)
pari('default(realprecision,100)')

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')

# ---------------------------------------------------------------------------
# Part 0.  PARI functions
# ---------------------------------------------------------------------------

# scaling ideal of (A,B), from the definition
pari('''scal(nf,A,B)=
  my(v=if(A!=0&&B!=0,concat(Vec(idealfactor(nf,A)[,1]),Vec(idealfactor(nf,B)[,1])),
          A==0,Vec(idealfactor(nf,B)[,1]),Vec(idealfactor(nf,A)[,1])),I=matid(poldegree(nf.pol)));
  foreach(Set(v),p,I=idealmul(nf,I,idealpow(nf,p,
    if(A==0,floor(idealval(nf,B,p)/6),B==0,floor(idealval(nf,A,p)/4),
       min(floor(idealval(nf,A,p)/4),floor(idealval(nf,B,p)/6))))));
  I'''.replace('\n', ' '))

# archimedean part of the height, and the height, from the definition
pari('''arch(nf,A,B)=my(x=nfeltembed(nf,A),y=nfeltembed(nf,B));
  prod(v=1,#x,my(h=max(abs(x[v])^3,abs(y[v])^2));if(v<=nf.r1,h,h^2),1.)'''.replace('\n', ' '))
pari('height(nf,A,B)=arch(nf,A,B)/idealnorm(nf,scal(nf,A,B))^12')

# nonzero z in the ideal I with |sigma_v(z)|<=R_v, up to sign
pari('''box(nf,I,R)=my(M=nf[5][1]*I,k=#R,e=1+1e-12,L=List(),G,V);
  G=sum(v=1,k,real(conj(M[v,]~)*M[v,])/R[v]^2);V=qfminim((G+G~)/2,k*e,,2)[3];
  for(i=1,#V,my(w=abs(M*V[,i]));
    if(vecmin(vector(k,v,w[v]<=R[v]*e)),listput(L,nfbasistoalg(nf,I*V[,i]))));
  Vec(L)'''.replace('\n', ' '))

# pairs (A,B) with scaling ideal c, invariant j and height <= X; modulo units one
# may take |sigma_v(A)|^3, |sigma_v(B)|^2 <= h_v
pari('''pairs(bnf,j,c,X)=my(nf=bnf.nf,T=X*idealnorm(nf,c)^12,k=nf.r1+nf.r2,L=List(),h,a,b);
  my(E=apply(u->nfeltembed(nf,u),bnf.fu));
  h=vector(k,v,T^(1./poldegree(nf.pol))*exp(6*sum(i=1,#E,abs(log(abs(E[i][v]))))));
  if(j==0,foreach(box(nf,idealpow(nf,c,6),vector(k,v,h[v]^(1/2))),z,
      listput(L,[0,z]);listput(L,[0,-z])),
     j==1728,foreach(box(nf,idealpow(nf,c,4),vector(k,v,h[v]^(1/3))),z,
      listput(L,[z,0]);listput(L,[-z,0])),
     a=4*(1728-j)/(27*j);
     foreach(box(nf,idealpow(nf,c,4),vector(k,v,h[v]^(1/3))),z0,foreach([z0,-z0],z,
      if(issquare(nfeltnorm(nf,a*z))&&nfeltispower(nf,a*z,2,&b),
       b=nfbasistoalg(nf,nfeltmul(nf,z,b));listput(L,[z,b]);listput(L,[z,-b])))));
  select(P->arch(nf,P[1],P[2])<=T*(1+1e-40)&&scal(nf,P[1],P[2])==c,Vec(L))'''.replace('\n', ' '))


def gp(code, *args):
    return pari(code)(*args)


def q(x):
    """Fraction or int -> PARI rational."""
    x = Fraction(x)
    return pari('%d/%d' % (x.numerator, x.denominator))


def iroot(x, k):
    """floor(x^(1/k)), exactly, for a nonnegative rational x."""
    t = int(Fraction(x).__floor__())
    return int(gp('(t,k)->sqrtnint(t,k)', t, k)) if t > 0 else 0


# ---------------------------------------------------------------------------
# Part 1.  The formulas of the paper
# ---------------------------------------------------------------------------

class Field:
    def __init__(self, pol, tex):
        self.tex = tex
        self.bnf = pari('bnfinit(%s,1)' % pol)
        assert int(gp('b->bnfcertify(b)', self.bnf)) == 1
        self.nf = gp('b->b.nf', self.bnf)
        self.d = int(gp('b->poldegree(b.pol)', self.bnf))
        r1, r2 = (int(x) for x in gp('b->b.sign', self.bnf))
        self.r = r1 + r2 - 1
        self.w = int(gp('b->b.tu[1]', self.bnf))
        self.h = int(gp('b->b.no', self.bnf))
        self.cyc = [int(x) for x in gp('b->b.cyc', self.bnf)]
        self.minkowski = float(gp('b->(4/Pi)^b.r2*(poldegree(b.pol))!/poldegree(b.pol)'
                                  '^poldegree(b.pol)*sqrt(abs(b.disc))', self.bnf))
        self.rho = gp('b->2.^b.r1*(2*Pi)^b.r2*b.no*b.reg/(b.tu[1]*sqrt(abs(b.disc)))', self.bnf)

    def q(self, n):                        # unit index q_n(K)
        return math.gcd(n, self.w) * n**self.r

    def tors(self, n):                     # |Cl[n]|
        return math.prod(math.gcd(n, e) for e in self.cyc)

    def slope(self, n):                    # q_n(K) rho_K / zeta_K(n)
        z = gp('(b,n)->if(poldegree(b.pol)==1,zeta(n),lfun(b.pol,n))', self.bnf, n)
        return float(self.q(n) * self.rho / z)

    def cls(self, I):                      # ideal class as a coordinate tuple
        return tuple(int(x) for x in gp('(b,I)->bnfisprincipal(b,I,0)', self.bnf, I))

    def coset(self, c, n):                 # its image in Cl/Cl^n
        return tuple(x % math.gcd(n, e) for x, e in zip(c, self.cyc))

    @lru_cache(None)
    def ideals(self, bound):
        L = gp('(nf,B)->ideallist(nf,B,4)', self.nf, bound)
        return [(N + 1, I) for N, row in enumerate(L) for I in row]

    @lru_cache(None)
    def class_reps(self):
        reps = {}
        for _, I in self.ideals(int(self.minkowski) + 1):
            reps.setdefault(self.cls(I), I)
        assert len(reps) == self.h
        return reps

    @lru_cache(None)
    def primes(self, Y):                   # sorted (norm, class) of primes of norm <= Y
        L = gp('(b,Y)->my(L=List());forprime(p=2,Y,foreach(idealprimedec(b,p,logint(Y,p)),P,'
               'listput(L,[p^P.f,bnfisprincipal(b,P,0)])));Vec(L)', self.bnf, Y)
        return sorted((int(N), tuple(int(x) for x in c)) for N, c in L)


class Family:
    """n, m, a, n_a, D^m and the coset C_j.

    Either j is rational or K has a single infinite place; in both cases the
    archimedean factor of D^6 is max(N^-3, N^-2) with N the absolute norm of a.
    """

    def __init__(self, K, j):
        self.K, self.j = K, q(j) if isinstance(j, (int, str)) else j
        one = gp('nf->matid(poldegree(nf.pol))', K.nf)
        if self.j == 0:
            self.n, self.m, self.Na, self.Dm = 6, 2, one, Fraction(1)
        elif self.j == 1728:
            self.n, self.m, self.Na, self.Dm = 4, 3, one, Fraction(1)
        else:
            self.n, self.m = 2, 6
            self.a = gp('(nf,j)->nfbasistoalg(nf,nfeltdiv(nf,4*(1728-j),27*j))', K.nf, self.j)
            self.Na = gp('(nf,a)->my(f=idealfactor(nf,a));if(#f~==0,matid(poldegree(nf.pol)),'
                         'idealfactorback(nf,f[,1],vector(#f~,i,'
                         'max(ceil(f[i,2]/2),ceil(f[i,2]/3)))~))', K.nf, self.a)
            N = Fraction(str(gp('(nf,I)->idealnorm(nf,I)', K.nf, self.Na)))
            A = abs(Fraction(str(gp('(nf,a)->nfeltnorm(nf,a)', K.nf, self.a))))
            self.Dm = N**6 * max(A**-3, A**-2)
        self.target = K.coset(tuple(-x for x in K.cls(self.Na)), self.n)
        self.mult = K.q(self.n) * K.tors(self.n)

    def max_norm(self, X):                 # largest N with D^m N^m <= X
        return iroot(Fraction(X) / self.Dm, self.m)

    def eligible(self, bound):             # n-free ideals of norm <= bound, class in C_j
        K = self.K
        free = '(nf,I,n)->my(f=idealfactor(nf,I));#f~==0||vecmax(f[,2])<n'
        return [(N, I) for N, I in K.ideals(bound) if K.coset(K.cls(I), self.n) == self.target
                and int(gp(free, K.nf, I, self.n))] if bound else []

    def count(self, X):
        """#E_{K,j}(X) by the exact formula."""
        return self.mult * len(self.eligible(self.max_norm(X)))

    def least(self):
        """Algorithm 1 with Remark 6.5: m_{K,j} and all minimizing pairs (A,B)."""
        K, n = self.K, self.n
        E = self.eligible(int(K.minkowski) + 1)
        m = E[0][0]
        out = []
        for I in [I for N, I in E if N == m]:
            for C in K.class_reps().values():
                J = gp('(nf,I,N,C,n)->idealmul(nf,idealmul(nf,I,N),idealpow(nf,C,n))',
                       K.nf, I, self.Na, C, n)
                if any(K.cls(J)):
                    continue
                g = gp('(b,J)->nfbasistoalg(b,bnfisprincipal(b,J)[2])', K.bnf, J)
                for ex in itertools.product(range(math.gcd(n, K.w)), *[range(n)] * K.r):
                    t = gp('(b,g,e)->my(u=concat([b.tu[2]],b.fu));for(i=1,#e,g=nfeltmul(b,g,'
                           'nfeltpow(b,u[i],e[i])));nfbasistoalg(b,g)', K.bnf, g, list(ex))
                    out.append((pari(0), t) if n == 6 else (t, pari(0)) if n == 4 else
                               (t**2 / self.a, t**3 / self.a))
        return m, out


def free_counts(K, n, Y):
    """{coset of Cl/Cl^n: number of n-free ideals of norm <= Y}, Algorithm 2."""
    cosets = list(itertools.product(*[range(math.gcd(n, e)) for e in K.cyc]))
    num = {c: i for i, c in enumerate(cosets)}
    add = [[num[K.coset([x + y for x, y in zip(a, b)], n)] for b in cosets] for a in cosets]
    P = [(N, num[K.coset(c, n)]) for N, c in K.primes(Y)] + [(Y + 1, 0)]
    cnt, stack = [0] * len(cosets), [(0, 1, 0)]
    while stack:
        i, N, c = stack.pop()
        cnt[c] += 1
        for k in range(i, len(P) - 1):
            p, g = P[k]
            if N * p > Y:
                break
            M, cc, nxt = N, c, P[k + 1][0]
            for _ in range(1, n):
                M *= p
                if M > Y:
                    break
                cc = add[cc][g]
                if M * nxt <= Y:
                    stack.append((k + 1, M, cc))
                else:
                    cnt[cc] += 1
    return dict(zip(cosets, cnt))


# ---------------------------------------------------------------------------
# Part 2.  Independent enumeration of isomorphism classes
# ---------------------------------------------------------------------------

def direct(K, j, X):
    """Classes of invariant j and height <= X, as (A, B, height, scaling ideal)."""
    j = Fraction(j)
    n = 6 if j == 0 else 4 if j == 1728 else 2
    out = []
    for c in K.class_reps().values():
        groups = defaultdict(list)
        for pr in gp('pairs', K.bnf, q(j), c, q(X)):
            A, B = pr[0], pr[1]
            H = gp('height', K.nf, A, B)
            z = B if n == 6 else A if n == 4 else B / A        # (A,B)~(A',B') iff z'/z in K^n
            g = groups[str(gp('h->precision(h,38)', H))[:30]]
            if not any(int(gp('(nf,x,n)->nfeltispower(nf,x,n)', K.nf, z / y[0], n)) for y in g):
                g.append((z, A, B, H, c))
        out += [y[1:] for g in groups.values() for y in g]
    return out


# ---------------------------------------------------------------------------
# Part 3.  The tables
# ---------------------------------------------------------------------------

LOG, FAILED = [], []


def say(s):
    print(s, flush=True)
    LOG.append(s)


def check(ok, what):
    say('[%s] %s' % ('ok' if ok else 'FAIL', what))
    if not ok:
        FAILED.append(what)


def write(name, spec, head, blocks, columns, rows):
    L = [r'\begin{tabular}{%s}' % spec, r'\toprule'] + head + [r'\midrule']
    L.append(('\n' + r'\midrule' + '\n').join(
        '\n'.join(' & '.join(r) + r'\\' for r in b) for b in blocks))
    with open(os.path.join(OUT, 'tab_%s.tex' % name), 'w') as f:
        f.write('\n'.join(L + [r'\bottomrule', r'\end{tabular}', '']))
    with open(os.path.join(OUT, name + '.csv'), 'w', newline='') as f:
        csv.writer(f).writerows([columns] + rows)


def M(x):
    return '$%s$' % x


Q, Qi = Field('y', r'\Q'), Field('y^2+1', r'\Q(i)')
Q2, Qm5 = Field('y^2-2', r'\Q(\sqrt2)'), Field('y^2+5', r'\Q(\sqrt{-5})')
Qm21, Qm23 = Field('y^2+21', r'\Q(\sqrt{-21})'), Field('y^2-y+6', r'\Q(\sqrt{-23})')
Qm3 = Field('y^2+y+1', r'\Q(\sqrt{-3})')


def table_counts():
    """Theorem 1.1 against enumeration; Theorem 1.4 at Y=10^6."""
    Y, rows, blocks = 10**6, [], []
    for K in [Q, Qi, Q2, Qm5, Qm23]:
        block = []
        # j=1728/7 has a=8/9, where both terms in the exponents of n_a matter
        for j, X in [('1728/7', 10**14), ('1728', 10**6), ('0', 10**4)]:
            F = Family(K, j)
            ideals = F.eligible(F.max_norm(X))
            E = direct(K, j, X)
            fibres, heights = Counter(), True
            for A, B, H, c in E:           # the ideal b([E]), and the height identity
                z = B / A if F.n == 2 else B if F.n == 6 else A
                b = gp('(nf,z,N,c,n)->idealdiv(nf,z,idealmul(nf,N,idealpow(nf,c,n)))',
                       K.nf, z, F.Na, c, F.n)
                N = int(gp('(nf,I)->idealnorm(nf,I)', K.nf, b))
                heights &= bool(gp('(x,y)->abs(x/y-1)<1e-60', H, q(F.Dm * N**F.m)))
                fibres[str(b)] += 1
            ok = heights and set(fibres) == {str(I) for _, I in ideals} \
                and set(fibres.values()) == {F.mult} and len(E) == F.mult * len(ideals)
            check(ok, 'Thm 1.1  %-14s j=%-8s X=%.0e: formula %d, enumeration %d'
                  % (K.tex, j, X, F.mult * len(ideals), len(E)))
            big, main = F.mult * free_counts(K, F.n, Y)[F.target], K.slope(F.n) * Y
            err = big / main - 1
            check(abs(err) < 1e-4, 'Thm 1.4  %-14s j=%-8s Y=1e6: exact %d, main term %.1f, '
                  'rel. error %+.1e' % (K.tex, j, big, main, err))
            e = math.floor(math.log10(abs(err)))
            rows.append([K.tex, j, K.q(F.n), K.tors(F.n), X, len(E), big, main, err])
            block.append([M(K.tex) if not block else '', M(j), M(K.q(F.n)), M(K.tors(F.n)),
                          M(len(E)), M(big), M('%.1f' % main),
                          M(r'%+.1f\cdot10^{%d}' % (err / 10.0**e, e))])
        blocks.append(block)
    write('counts', 'llccrrrr',
          [r'& & & & & \multicolumn{3}{c}{$Y=10^6$}\\', r'\cmidrule(l){6-8}',
           r'$K$ & $j$ & $q_n$ & $|\Cl[n]|$ & $\#\cE_{K,j}(X_0)$ & $\#\cE_{K,j}$ & main term '
           r'& rel.\ error\\'], blocks,
          ['field', 'j', 'q_n', 'Cl[n]', 'X', 'count_X', 'count_Y1e6', 'main_term', 'rel_error'], rows)


def model(K, A, B):
    """TeX for y^2=x^3+Ax+B; table_least prints a shortest one."""
    gen = r'\sqrt{-21}' if K is Qm21 else r'\sqrt2' if K is Q2 else r'\omega'

    def unit(c, t):                        # signed monomial c*t
        return '' if c == 0 else ('+' if c > 0 else '-') + (
            str(abs(c)) if abs(c) != 1 or not t else '') + t

    def term(x, t):
        a, b = (int(c) for c in gp('x->my(p=lift(x));[polcoeff(p,0),polcoeff(p,1)]', x))
        if a == 0 or b == 0:
            return unit(a, t) + unit(b, gen + t)
        s = 1 if a > 0 else -1
        return '%s(%d%s)%s' % ('+' if a > 0 else '-', abs(a), unit(s * b, gen), t)
    return 'y^2=x^3%s%s' % (term(A, 'x'), term(B, ''))


def table_least():
    """Theorem 1.9 / Algorithm 1 against enumeration up to the least height."""
    rows, blocks = [], []
    for block in [[(Q, '1728/7'), (Q2, '1728/7'), (Qm23, '1728/7'), (Qm23, '0')],
                  [(Qm5, '-3375'), (Qm5, '3456/29')],
                  [(Qm21, '6912/31'), (Qm21, '3456/29'), (Qm21, '6912/85'), (Qm21, '3456/83')]]:
        tex, last = [], None
        for K, j in block:
            F = Family(K, j)
            m, P = F.least()
            H = F.Dm * m**F.m
            ok = all(gp('(x,y)->abs(x/y-1)<1e-60', gp('height', K.nf, A, B), q(H))
                     and 6912 * A**3 == q(j) * (4 * A**3 + 27 * B**2) for A, B in P)
            E = direct(K, j, H)            # every class of height <= H
            ok &= len(E) == len(P) and all(gp('(x,y)->abs(x/y-1)<1e-60', e[2], q(H)) for e in E)
            check(ok and H.denominator == 1, 'Thm 1.9  %-14s j=%-8s m=%d least height %s, '
                  '%d classes; enumeration %d, none lower' % (K.tex, j, m, H, len(P), len(E)))
            rows.append([K.tex, j, m, H, len(P), '; '.join(
                '(%s, %s)' % tuple(gp('x->lift(x)', x) for x in p) for p in P)])
            A, B = min(P, key=lambda p: (len(model(K, *p)), model(K, *p).count('-')))
            cl = r'\times'.join(r'\Z/%d' % e for e in K.cyc) or '1'
            tex.append([M(K.tex) if K is not last else '', M(cl) if K is not last else '', M(j),
                        M(m), M(H), M(len(P)), M(model(K, A, B))])
            last = K
        blocks.append(tex)
    write('least', 'lccrrrl', [r'$K$ & $\Cl$ & $j$ & $m_{K,j}$ & least height & classes '
                               r'& a minimizing model\\'], blocks,
          ['field', 'j', 'm_Kj', 'least_height', 'classes', 'minimizing (A,B)'], rows)


CM = [0, 1728, -3375, 8000, -32768, 54000, 287496, -12288000, 16581375, -884736, -884736000,
      -147197952000, -262537412640768000]


def free_count_h1(K, n, Y):
    """Number of n-free ideals of norm <= Y for K=Q, Q(i) or Q(sqrt-3), by a Moebius sieve."""
    if Y < 1:
        return 0

    def ideals(R):                         # number of nonzero ideals of norm <= R
        if K is Q:
            return R
        if K is Qi:                        # a^2+b^2 <= R, over the 4 units
            s = math.isqrt(R)
            return (2 * sum(2 * math.isqrt(R - a * a) + 1 for a in range(1, s + 1)) + 2 * s) // 4
        G = 0                              # 4(x^2-xy+y^2) = u^2+3y^2 with u=2x-y, over the 6 units
        for y in range(-math.isqrt(4 * R // 3), math.isqrt(4 * R // 3) + 1):
            t = math.isqrt(4 * R - 3 * y * y)
            G += 2 * (t // 2) + 1 if y % 2 == 0 else 2 * ((t + 1) // 2)
        return (G - 1) // 6
    T = iroot(Y, n)
    P, total, stack = [N for N, _ in K.primes(max(T, 2))], 0, [(0, 1, 1)]
    while stack:
        i, N, s = stack.pop()
        total += s * ideals(Y // N**n)
        for k in range(i, len(P)):
            if N * P[k] > T:
                break
            stack.append((k + 1, N * P[k], -s))
    return total


def table_cm():
    """Theorem 1.6 over Q, Q(i), Q(sqrt-3), where the CM invariants are the 13 rational ones."""
    rows, blocks = [], []
    for K in [Q, Qi, Qm3]:
        F = [Family(K, j) for j in CM]
        check(all(free_count_h1(K, n, 10**4) == sum(free_counts(K, n, 10**4).values())
                  for n in (2, 4, 6)), 'sieve and Algorithm 2 agree over ' + K.tex)
        lim, block = K.slope(4) / K.slope(6), []
        for e in (6, 12, 18, 24):
            X = 10**e
            N = [f.mult * free_count_h1(K, f.n, f.max_norm(X)) for f in F]
            tot = sum(N)
            v = float(Fraction(tot - N[0], tot)) * 10**(e / 6)
            say('Thm 1.6  %-6s X=1e%d: CM %d = %d + %d + %d, X^{1/6}(1-p0)=%.4f, limit %.4f'
                % (K.tex, e, tot, N[0], N[1], sum(N[2:]), v, lim))
            rows.append([K.tex, X, tot, N[0], N[1], sum(N[2:]), N[0] / tot, v, lim])
            block.append([M(K.tex) if not block else '', M('10^{%d}' % e), M(tot), M(N[0]), M(N[1]),
                          M(sum(N[2:])), M('%.4f' % (N[0] / tot)), M('%.4f' % v)])
        check(abs(v / lim - 1) < 1e-3, 'Thm 1.6  %s: within 0.1%% of the limit at X=1e24' % K.tex)
        block.append(['', 'limit', '', '', '', '', M(1), M('%.4f' % lim)])
        blocks.append(block)
    write('cm', 'lrrrrrrr', [r'$K$ & $X$ & $\#\cE_K^{\CM}(X)$ & $j=0$ & $j=1728$ & other & $p_0(X)$ '
                             r'& $X^{1/6}(1-p_0(X))$\\'], blocks,
          ['field', 'X', 'CM', 'j=0', 'j=1728', 'other', 'p0', 'X^(1/6)(1-p0)', 'limit'], rows)


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    say('cypari %s, PARI %s' % (cypari.__version__, '.'.join(str(x) for x in pari('version()')[:3])))
    table_counts()
    table_least()
    table_cm()
    say('%d checks, %d failures' % (sum(s.startswith(('[ok]', '[FAIL]')) for s in LOG), len(FAILED)))
    with open(os.path.join(OUT, 'report.txt'), 'w') as f:
        f.write('\n'.join(LOG) + '\n')
    sys.exit(1 if FAILED else 0)
