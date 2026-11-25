"""Validated Fourier-Taylor calculation at finite R, using integer intervals.

All interval endpoints in the main recurrence share the scale 2**224.
The mathematical error estimate is in lag_bins_proof.md.
"""

from fractions import Fraction as F

import sys
from pathlib import Path

# Sibling module import that also works under python -I (script dir not on sys.path).
sys.path.insert(0, str(Path(__file__).resolve().parent))

from verify_lag_bins_exact import AMP, H, PI_UP, decimal_interval, log_lower


def require(condition, label):
    """Fail loudly, also under python -O (unlike assert)."""
    if not condition:
        raise SystemExit(f"CHECK FAILED: {label}")


S = 2**224
ZERO = (0, 0)
CZ = (ZERO, ZERO)
EDGES = tuple(map(F, ['1/4','1/2','3/4','1','3/2','2','3','4','6','10','16']))
WEIGHTS = tuple(F(n, 100) for n in (90,94,99,104,111,121,132,136,121,90,67))
U = 16


def iv(x):
    n, d = x.numerator*S, x.denominator
    return n//d, -((-n)//d)


def add(x, y):
    return x[0]+y[0], x[1]+y[1]


def neg(x):
    return -x[1], -x[0]


def product(x, y):
    ends = (x[0]*y[0], x[0]*y[1], x[1]*y[0], x[1]*y[1])
    return min(ends)//S, -((-max(ends))//S)


def cp(x, a):
    return product(x[0], a), product(x[1], a)


def plus(x, y):
    return add(x[0], y[0]), add(x[1], y[1])


def minus(x):
    return neg(x[0]), neg(x[1])


def times_i(x):
    return neg(x[1]), x[0]


def add_at(dic, k, val):
    dic[k] = plus(dic.get(k, CZ), val)


def atan_inv(n):
    # Enclose atan(1/n) by adjacent alternating Taylor sums.
    partial = sum(F((-1)**k, (2*k+1)*n**(2*k+1)) for k in range(75))
    return partial-F(1,151*n**151), partial


def pi_interval():
    a0,a1 = atan_inv(5)
    b0,b1 = atan_inv(239)
    return (iv(16*a0 - 4*b1)[0], iv(16*a1 - 4*b0)[1])


def phase_weights(h):
    prefactor = {0: (F(1),F(0)), 2: (-h/4,-h/4), -2: (-h/4,h/4)}
    rho = {}
    for k,(a,b) in prefactor.items():
        for j,v in ((-1,AMP/2),(0,F(1)),(1,AMP/2)):
            c,d = rho.get(k+j,(F(0),F(0)))
            rho[k+j] = c+v*a,d+v*b
    return {k:(iv(a),iv(b)) for k,(a,b) in rho.items()}


def expectation(dfield, rho):
    result = ZERO
    for k,(x,y) in rho.items():
        a,b = dfield.get(-k,CZ)
        result = add(result,add(product(a,x),neg(product(b,y))))
    return result


def apply_generator(dfield,cfield,delta,h,complete):
    """Coefficients of A acting on fields for the doublet and singleton."""
    out_d,out_c = {},{}
    wave = product(delta,iv(h/2))
    half_amp = iv(AMP/2)
    for k in dfield.keys() | cfield.keys():
        d,c = dfield.get(k,CZ),cfield.get(k,CZ)
        diff = plus(c,minus(d)) if not complete else minus(d)
        drift = cp(times_i(d),product(delta,iv(F(k))))
        transfer = cp(diff,delta)
        add_at(out_d,k,plus(drift,plus(transfer,minus(d))))
        side = cp(minus(d),half_amp)
        add_at(out_d,k-1,side)
        add_at(out_d,k+1,side)
        harmonic = cp(times_i(diff),wave)
        add_at(out_d,k+2,harmonic)
        add_at(out_d,k-2,minus(harmonic))
        if not complete:
            add_at(out_c,k,plus(cp(times_i(c),product(delta,iv(F(k)))),minus(transfer)))
            add_at(out_c,k+2,harmonic)
            add_at(out_c,k-2,minus(harmonic))
    return out_d,out_c


def abs_upper(field):
    return sum(max(-x[0],x[1])+max(-y[0],y[1]) for x,y in field.values())


def cdf(speed,h,complete):
    # V_n = U**(n+1) A**n (q,0) / (n+1)!.
    delta = product(pi_interval(),iv(F(2,speed)))
    rho = phase_weights(h)
    d = {0:(iv(F(U)),ZERO),-1:(iv(U*AMP/2),ZERO),1:(iv(U*AMP/2),ZERO)}
    c = {}
    projected = []
    for n in range(220):
        tail = max(abs_upper(d),abs_upper(c))
        # Semigroup contraction bounds the whole remaining Taylor error
        # by this next derivative on all requested 0<=u<=U.
        if tail*10**12 < S:
            break
        projected.append(expectation(d,rho))
        d,c = apply_generator(d,c,delta,h,complete)
        multiplier = iv(F(U,n+2))
        d={k:cp(v,multiplier) for k,v in d.items()}
        c={k:cp(v,multiplier) for k,v in c.items()}
    else:
        raise ArithmeticError('Taylor remainder did not reach the requested enclosure')

    bounds=[]
    for edge in EDGES:
        value = ZERO
        z = iv(edge/U)
        for coeff in reversed(projected):
            value = product(z,add(value,coeff))
        lo,hi = value
        bounds.append((F(lo-tail,S),F(hi+tail,S)))
    require(all(F(0)<lo<=hi<F(1) and hi-lo<F(1,10**10) for lo,hi in bounds), 'verify_moderate_speed:1')
    return bounds


def lower_bound(speed,complete):
    forward=cdf(speed,H,complete)
    reverse=cdf(speed,-H,complete)
    prev_f = prev_b = (F(0), F(0))
    total=F(0)
    for w,f,b in zip(WEIGHTS,forward,reverse):
        fmass=(f[0]-prev_f[1],f[1]-prev_f[0])
        bmass=(b[0]-prev_b[1],b[1]-prev_b[0])
        require(fmass[0]>0 and bmass[0]>0, 'verify_moderate_speed:2')
        score=log_lower(w)
        total+=min(score*v for v in fmass)-max((w-1)*v for v in bmass)
        prev_f,prev_b=f,b
    # The overflow bin and all non-adjacent-reporter pairs have weight 1.
    return F(speed,2)*total


if __name__=='__main__':
    for speed in (750,1000):
        for complete in (False,True):
            result=lower_bound(speed,complete)
            print('R =',speed,'complete observation =',complete,
                  'certified lower bound interval',decimal_interval(result,6),flush=True)
            require(result>4*PI_UP/5, 'verify_moderate_speed:3')
