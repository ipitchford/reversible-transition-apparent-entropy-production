"""Independent periodic backward-equation and fast-limit checks.

The backward equations integrate over future changes in the slow state and
the reporter hazard. They have no frozen-phase approximation.
"""

import math

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import roots_laguerre, xlogy


def require(condition, label):
    """Fail loudly, also under python -O (unlike assert)."""
    if not condition:
        raise SystemExit(f"CHECK FAILED: {label}")


H = 4 / 5
AMPLITUDE = 16 / 17
TAU = 2 * math.pi


def phase_functions(t, h, amplitude):
    theta = TAU * t
    drive = h * np.sin(2 * theta)
    a = TAU * (1 + drive)
    b = TAU * (1 - drive)
    q = 1 + amplitude * np.cos(theta)
    p = 0.5 + h * (np.sin(2 * theta) - np.cos(2 * theta)) / 4
    return a, b, q, p


def race_probabilities(speed, h, amplitude, tolerance):
    # u(t) is the probability of a visible event before an independent
    # Exp(R) deadline. Its backward equation is stable backwards in time.
    # The first two entries start in the active doublet and the singleton.
    # The third entry additionally requires the reporter to be the next
    # physical jump, as in the test with all six directed edges observed.
    def backward(t, u):
        a, b, q, _ = phase_functions(t, h, amplitude)
        total_active = speed * (1 + q) + b
        return (total_active * u[0] - b * u[1] - speed * q,
                (speed + a) * u[1] - a * u[0],
                total_active * u[2] - speed * q)

    def jacobian(t, u):
        a, b, q, _ = phase_functions(t, h, amplitude)
        total_active = speed * (1 + q) + b
        return np.array([[total_active, -b, 0.0],
                         [-a, speed + a, 0.0],
                         [0.0, 0.0, total_active]])

    # Start one extra period into the future. With terminal value zero,
    # missed probability on the used period is at most exp(-R).
    t = np.arange(1024) / 1024
    sol = solve_ivp(backward, (2.0, 0.0), np.zeros(3),
                    method='Radau', jac=jacobian, t_eval=t[::-1],
                    rtol=tolerance, atol=tolerance * 0.02)
    if not sol.success:
        raise RuntimeError(sol.message)
    active_race, _, direct_race = sol.y[:, ::-1]
    a, b, q, p = phase_functions(t, h, amplitude)
    event_weights = q * p
    visible_rate = speed * event_weights.mean()
    full_rate = visible_rate + np.mean((1 - p) * a + p * b)
    mass = np.average(active_race, weights=event_weights)
    full_mass = visible_rate / full_rate * np.average(direct_race, weights=event_weights)
    return mass, full_mass, visible_rate, full_rate


def bernoulli_kl(p, q):
    return xlogy(p, p / q) + xlogy(1 - p, (1 - p) / (1 - q))


def experiment(speed, tolerance, h=H, amplitude=AMPLITUDE):
    f, all_f, rate, all_rate = race_probabilities(speed, h, amplitude, tolerance)
    b, all_b, rev_rate, all_rev_rate = race_probabilities(speed, -h, amplitude, tolerance)
    require(abs(rate - rev_rate) < 1e-8, 'verify_numerically:1')
    require(abs(all_rate - all_rev_rate) < 1e-8, 'verify_numerically:2')
    return np.array([f, b, rate * bernoulli_kl(f, b),
                     all_f, all_b, all_rate * bernoulli_kl(all_f, all_b)])


def slope(angle_count, laguerre_order):
    # For g(u) = <p q^2 exp(-q u)>/<p q>, remove exp(-q_min u)
    # to keep the density ratio stable at the largest quadrature nodes.
    t = np.arange(angle_count) / angle_count
    _, _, q, pf = phase_functions(t, H, AMPLITUDE)
    pb = 1 - pf
    active = np.mean(pf * q)
    require(abs(active - 0.5) < 1e-14, 'verify_numerically:3')
    qmin = 1 - AMPLITUDE

    def scaled_integrand(u):
        rates = q*q * np.exp(-(q - qmin)*u) / active
        f, b = np.mean(pf*rates), np.mean(pb*rates)
        return f * np.log(f/b)

    nodes, weights = roots_laguerre(laguerre_order)
    gauss = sum(w * scaled_integrand(x/qmin) for x, w in zip(nodes, weights)) / qmin
    adaptive, error = quad(lambda u: math.exp(-qmin*u)*scaled_integrand(u),
                           0, np.inf, epsabs=2e-13, epsrel=2e-11, limit=200)
    return active * gauss, active * adaptive, active * error


def main():
    # Verify the entropy independently with the instantaneous flux formula.
    def entropy(t):
        a, b, _, p = phase_functions(t, H, AMPLITUDE)
        forward, reverse = (1-p)*a, p*b
        return (forward-reverse) * np.log(forward/reverse)

    sigma = quad(entropy, 0, 1, epsabs=1e-12, epsrel=1e-12)[0]
    require(abs(sigma - 4*np.pi/5) < 1e-12, 'verify_numerically:4')
    print(f'Entropy production (integrated flux): {sigma:.12g}')
    print('Finite-speed validation. Columns are R, forward/reverse race probabilities,')
    print('their one-observed-edge KL rate lower bound, and the same three quantities')
    print('when all physical jumps are observed. KL rates are numerical lower bounds on Eq. 14.')
    for speed in (100, 1_000, 100_000):
        results = experiment(speed, 2e-10)
        check = experiment(speed, 8e-9)
        require(max(abs(check[:2] - results[:2])) < 5e-8, 'verify_numerically:5')
        require(max(abs(check[3:5] - results[3:5])) < 5e-8, 'verify_numerically:6')
        print(f'{speed:6d} ' + ' '.join(f'{x:.11g}' for x in results))
        if speed == 100_000:
            require(min(results[2], results[5]) > 30*sigma, 'verify_numerically:7')
            require(abs(results[1] - results[0] - 17/600) < 4e-4, 'verify_numerically:8')

    print('Negative controls. Each forward/reverse probability difference should vanish.')
    for title, h, amp in [('undriven occupation', 0, AMPLITUDE),
                          ('constant reporter rate', H, 0)]:
        values = experiment(1000, 2e-10, h, amp)
        differences = [values[1] - values[0], values[4] - values[3]]
        require(max(abs(x) for x in differences) < 3e-9, 'verify_numerically:9')
        print(title, *(f'{x:.3g}' for x in differences))

    print('Exact-limit density slope K_F D(g_F||g_B)/R, quadrature values.')
    print('phase nodes, Laguerre order, Gauss-Laguerre, adaptive, adaptive error estimate')
    for nodes, order in [(256, 128), (1024, 256)]:
        gauss, adaptive, error = slope(nodes, order)
        print(nodes, order, f'{gauss:.12g}', f'{adaptive:.12g}', f'{error:.2g}')
    require(abs(gauss - adaptive) < 2e-8, 'verify_numerically:10')
    require(adaptive > (17/600)**2, 'verify_numerically:11')


if __name__ == '__main__':
    main()
