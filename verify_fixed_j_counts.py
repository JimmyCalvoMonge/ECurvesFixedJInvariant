"""Independent finite checks for fixed-j twist counts.
Uses (1) exact enumeration of elements and square-class comparison,
(2) prime-ideal Euler factors. Requires sympy, no Sage/PARI.
This code checks the j != 0,1728 family with a(j)=1 over Q and Q(sqrt 2),
and a(j)=2 over Q(sqrt(-5)), with alpha=beta=1.
"""
from collections import defaultdict
from functools import lru_cache
from fractions import Fraction
from math import gcd, isqrt, sqrt
from sympy import factorint, kronecker_symbol, primerange


def is_rational_square(x):
    q=Fraction(x)
    if q < 0: return False
    return isqrt(q.numerator)**2 == q.numerator and isqrt(q.denominator)**2 == q.denominator


def is_square_quadratic(a,b,D):
    """Test whether a+b*sqrt(D) is a square in Q(sqrt(D)), exact fractions."""
    a,b=Fraction(a),Fraction(b)
    if b==0:
        return is_rational_square(a) or is_rational_square(a/Fraction(D))
    norm=a*a-D*b*b
    if not is_rational_square(norm): return False
    normroot=Fraction(isqrt(norm.numerator),isqrt(norm.denominator))
    return is_rational_square((a+normroot)/2) or is_rational_square((a-normroot)/2)


def equivalent_quadratic(x,y,D):
    """x/y square in K, x,y integer pairs."""
    a,b=x; c,d=y
    n=c*c-D*d*d
    return is_square_quadratic(Fraction(a*c-D*b*d,n),Fraction(b*c-a*d,n),D)


def pvaluation(n,p):
    count=0
    while n and n%p==0:
        n//=p; count+=1
    return count


def valuations_quadratic(t,D):
    """Return list (prime_norm, exponent) in the principal ideal (t).
    Q(sqrt D), D=2 or -5, discriminant 8 or -20 (integral basis 1,sqrt D).
    For unramified split odd primes use Hensel lifts of square roots mod p^e.
    """
    a,b=t; N=abs(a*a-D*b*b)
    if N==0: raise ValueError('zero')
    out=[]
    for p,e in factorint(N).items():
        disc=4*D
        if disc%p==0:
            out.append((p,e)); continue  # ramified, residue degree one
        k=int(kronecker_symbol(disc,p))
        if k==-1:
            assert e%2==0
            out.append((p*p,e//2));continue
        if k!=1: raise ValueError((p,k))
        roots=[r for r in range(p) if (r*r-D)%p==0]
        assert len(roots)==2
        v=0;root=roots[0];mod=p
        for depth in range(1,e+1):
            if depth>1:
                root=next(root+j*(mod//p) for j in range(p) if (root+j*(mod//p))**2%mod==D%mod)
            if (a+b*root)%mod==0: v+=1
            else: break
            mod*=p
        if v: out.append((p,v))
        if e-v: out.append((p,e-v))
    check=1
    for q,k in out:check*=q**k
    assert check==N,(t,D,N,out,check)
    return out


def bnorm_quadratic(t,D,a_parameter):
    """Norm of b(t) via local parity, independent direct local arithmetic."""
    if D not in [2,-5] or a_parameter not in [1,2]:raise ValueError
    vals=valuations_quadratic(t,D)
    norm=1
    # only ramified prime p2 divides a=2 in K(sqrt(-5)), v_p(2)=2, m=1
    m_at_2=int(D==-5 and a_parameter==2)
    if m_at_2 and 2 not in {q for q,_ in vals}: vals.append((2,0))
    for q,e in vals:
        m=m_at_2 if q==2 else 0
        eps=(e-m)%2
        norm*=q**eps
    return norm


def ideal_counts(Y,D):
    """Euler-product enumeration of squarefree *ideals*, by norm and class parity.
    D=2: class number 1, D=-5: class group of order two."""
    dp={(1,0):1}
    for p in primerange(2,Y+1):
        if p==2:
            opts=[(1,0,1),(2, 1 if D==-5 else 0,1)]
        elif D==-5 and p==5:
            opts=[(1,0,1),(5,0,1)]
        else:
            split=int(kronecker_symbol(4*D,p))
            if split==1:
                if D==2: cl=0
                else:
                    # principal iff p=a^2+5b^2 (norm of element)
                    cl=int(not any(u*u+5*v*v==p for v in range(isqrt(p//5)+1) for u in range(isqrt(p)+1)))
                opts=[(1,0,1),(p,cl,2)]
                if p*p<=Y:opts.append((p*p,0,1))
            elif split==-1:
                opts=[(1,0,1)]+([(p*p,0,1)] if p*p<=Y else [])
            else: raise ValueError((D,p,split))
        if len(opts)==1:continue
        ndp=defaultdict(int)
        for (n,cl),v in dp.items():
            for q,cp,m in opts:
                if n*q<=Y:ndp[(n*q,(cl+cp)%2)]+=v*m
        dp=ndp
    return dp


def brute_real_pid(Y):
    """Enumerate balanced t in Z[sqrt2] whose principal ideal is squarefree."""
    eps=1+sqrt(2)
    upper=eps**2
    bound=int(5*sqrt(Y)+12)
    selected=[]
    for a in range(-bound,bound+1):
        for b in range(-bound,bound+1):
            n=abs(a*a-2*b*b)
            if not (1<=n<=Y):continue
            ratio=abs(a+b*sqrt(2))/sqrt(n)
            if ratio < 1-1e-10 or ratio >=upper-1e-10: continue
            if all(e==1 for _,e in valuations_quadratic((a,b),2)):
                selected.append((a,b,n))
    # quotient by square equivalence test just for small Y
    reps=[]
    for a,b,n in selected:
        if not any(equivalent_quadratic((a,b), (x,y),2) for x,y,_ in reps):
            reps.append((a,b,n))
    return reps


def brute_imag_nonpid(Y):
    """For b squarefree with norm <=Y, all squareclasses have an integral t of norm <=8Y.
    This follows from the two class-group square roots, represented by O_K and p_2 (norm 2).
    Enumeration/dedup uses only direct element arithmetic and square testing.
    """
    alimit=isqrt(8*Y);blimit=isqrt(8*Y//5)
    classes=[]
    for a in range(-alimit,alimit+1):
        for b in range(-blimit,blimit+1):
            n=a*a+5*b*b
            if n<1 or n>8*Y:continue
            bn=bnorm_quadratic((a,b),-5,2)
            if bn>Y:continue
            if not any(equivalent_quadratic((a,b),(u,v),-5) for u,v,_,_ in classes):
                classes.append((a,b,n,bn))
    return classes


def main():
    print('Q: j=6912/31, a=1, alpha=beta=1')
    for Y in [1,2,3,5,10,20]:
        sq=sum(all(n%(p*p)!=0 for p in factorint(n)) for n in range(1,Y+1))
        print('Y=%3d height X=Y^6 count=%d' %(Y,2*sq))
    print('\nQ(sqrt 2): j=6912/31, a=1, rank 1 units, PID')
    for Y in [1,2,3,5,10,15]:
        dp=ideal_counts(Y,2)
        pred=4*sum(v for (norm,cl),v in dp.items())
        reps=brute_real_pid(Y)
        print('Y=%3d predicted=%3d direct=%3d'% (Y,pred,len(reps)))
        if len(reps)!=pred:print('DISAGREE',reps)
        assert len(reps)==pred
    print('\nQ(sqrt -5): j=3456/29, a=2, class group C2, N_a=p_2')
    for Y in [1,2,3,4,5,6,8,10,15]:
        dp=ideal_counts(Y,-5)
        elig=sum(v for (norm,cl),v in dp.items() if cl==1)
        predicted=4*elig
        direct=brute_imag_nonpid(Y)
        print('Y=%3d predicted=%3d direct=%3d (eligible squarefree ideals=%d)' % (Y,predicted,len(direct),elig))
        if predicted!=len(direct): print('DISAGREE',direct)
        assert predicted==len(direct)
    print('minimal representatives Q(sqrt-5), Y=2:',brute_imag_nonpid(2))
    print('All independent finite counts agree.')

if __name__=='__main__':main()
