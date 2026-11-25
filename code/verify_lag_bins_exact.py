"""Rational certificate for a nine-bin waiting-time witness.

The accompanying lag_bins_proof.md gives the error estimates used here.
This script uses only the Python standard library.
"""

from fractions import Fraction as F
from functools import cache
from math import comb, factorial


def require(condition, label):
    """Fail loudly, also under python -O (unlike assert)."""
    if not condition:
        raise SystemExit(f"CHECK FAILED: {label}")


H, AMP = F(4, 5), F(16, 17)
PI_UP = F(22, 7)
EDGES = tuple(map(F, ['1/4', '1/2', '1', '2', '3', '5', '8', '12']))
WEIGHTS = tuple(F(x, 100) for x in (90, 94, 101, 114, 132, 132, 103, 77, 56))


@cache
def cos_moment(n):
    return F(comb(n, n//2), 2**n) if n % 2 == 0 else F(0)


@cache
def frozen_moment(n, sign):
    # The event phase density is 2pq. Its sine component integrates to
    # zero against a function of cos(theta). sign=-1 is the forward arm.
    return sum(comb(n+1, j) * AMP**j *
               (cos_moment(j) + sign*H/2*(2*cos_moment(j+2)-cos_moment(j)))
               for j in range(n+2))


def frozen_survival(u, sign):
    # Taylor polynomials of odd/even degree sandwich exp(-x) for all x>=0.
    odd_order = 111
    lower = sum((-u)**n*frozen_moment(n, sign)/factorial(n)
                for n in range(odd_order+1))
    n = odd_order+1
    upper = lower + u**n*frozen_moment(n, sign)/factorial(n)
    require(F(0) < lower <= upper < F(1), 'verify_lag_bins_exact:1')
    require(upper-lower < F(1, 10**20), 'verify_lag_bins_exact:2')
    return lower, upper


def survival_bounds(speed, sign):
    bounds = []
    for u in EDGES:
        lo, hi = frozen_survival(u, sign)
        error = PI_UP*AMP*u*u/speed
        require(error < 1, 'verify_lag_bins_exact:3')
        # The integral of the time-dependent hazard differs from u*q(s)
        # by at most error. Use e^(-x)>=1-x and e^x<=1/(1-x).
        bounds.append((lo*(1-error), min(F(1), hi/(1-error))))
    return bounds


def log_lower(x):
    t = (x-1)/(x+1)
    n = 20
    polynomial = 2*sum(t**(2*k+1)/(2*k+1) for k in range(n))
    remainder = 2*abs(t)**(2*n+1)/((2*n+1)*(1-t*t))
    return polynomial - remainder


def expectation_range(values, survival):
    lo = hi = values[0]
    for last, nxt, (lower, upper) in zip(values, values[1:], survival):
        change = nxt-last
        lo += min(change*lower, change*upper)
        hi += max(change*lower, change*upper)
    return lo, hi


def bound(speed):
    flog = tuple(map(log_lower, WEIGHTS))
    bvalue = tuple(w-1 for w in WEIGHTS)
    f_lower, _ = expectation_range(flog, survival_bounds(speed, -1))
    _, b_upper = expectation_range(bvalue, survival_bounds(speed, 1))

    # Starting-phase distribution = 2pq. The Palm identity <qL>=1/R
    # gives a global bound on the chance of a hidden exit before the
    # always-on comparison reporter. We enclose p_max with 4/5.
    exit_upper = 2*F(4, 5)*2*PI_UP*(1+H)/speed
    require(min(flog) < 0 < max(flog), 'verify_lag_bins_exact:4')
    require(min(bvalue) < 0 < max(bvalue), 'verify_lag_bins_exact:5')
    correction = exit_upper*(max(flog)-min(flog)+max(bvalue)-min(bvalue))
    rate_lower = F(speed, 2)*(f_lower-b_upper-correction)
    return rate_lower, rate_lower/(4*PI_UP/5)


def decimal_interval(x, places):
    scale = 10**places
    lower = x.numerator*scale//x.denominator
    return f'[{lower}/{scale}, {lower+1}/{scale}]'


if __name__ == '__main__':
    for speed in (10_000, 100_000):
        rate, ratio = bound(speed)
        print(f'R = {speed}')
        print('  common lower bound on estimator, evaluated within', decimal_interval(rate, 5))
        print('  common lower bound on ratio, evaluated within', decimal_interval(ratio, 5))
        if speed == 100_000:
            require(ratio > 140, 'verify_lag_bins_exact:6')
        else:
            require(ratio > 1, 'verify_lag_bins_exact:7')
