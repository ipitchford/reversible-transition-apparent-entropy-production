"""Check a finite-bin statistic directly from killed forward dynamics."""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import rel_entr

import sys
from pathlib import Path

# Sibling module import that also works under python -I (script dir not on sys.path).
sys.path.insert(0, str(Path(__file__).resolve().parent))

from verify_numerically import H, AMPLITUDE, phase_functions


def require(condition, label):
    """Fail loudly, also under python -O (unlike assert)."""
    if not condition:
        raise SystemExit(f"CHECK FAILED: {label}")


EDGES = np.array([0, 0.1, 0.3, 0.5, 0.8, 1.2, 1.7, 2.5, 3.5, 5, 7, 12, 20])


def histogram(speed, h, nodes, tolerance):
    # Average over the starting phase only after propagating the actual
    # event-conditioned state to every bin boundary.
    start_phase = np.arange(nodes) / nodes
    _, _, q0, p0 = phase_functions(start_phase, h, AMPLITUDE)
    event_weights = p0 * q0
    initial = np.concatenate([np.ones(nodes), np.zeros(nodes)])

    def killed_flow(u, state):
        active, inactive = np.split(state, 2)
        a, b, q, _ = phase_functions(start_phase + u/speed, h, AMPLITUDE)
        net_exchange = (b*active - a*inactive) / speed
        return np.concatenate([-net_exchange - q*active, net_exchange])

    solution = solve_ivp(killed_flow, (0, EDGES[-1]), initial, method='DOP853',
                         t_eval=EDGES, rtol=tolerance, atol=tolerance * 0.02)
    if not solution.success:
        raise RuntimeError(solution.message)
    active, inactive = np.split(solution.y, 2, axis=0)
    survival = np.average(active + inactive, axis=0, weights=event_weights)
    probabilities = np.append(-np.diff(survival), survival[-1])
    require(np.all(probabilities > 0), 'verify_histogram:1')
    require(abs(probabilities.sum() - 1) < 1e-12, 'verify_histogram:2')
    return probabilities


def lower_bound(speed, nodes, tolerance):
    f = histogram(speed, H, nodes, tolerance)
    b = histogram(speed, -H, nodes, tolerance)
    return speed/2 * np.sum(rel_entr(f, b))


def main():
    sigma = 4*np.pi/5
    print('Thirteen waiting-time bins, including a single overflow bin above R*Twait = 20.')
    print('R, one-edge Eq. 14 lower bound (numerical), ratio to physical entropy')
    for speed in (100, 500, 750, 1000, 100_000):
        coarse = lower_bound(speed, 64, 5e-9)
        refined = lower_bound(speed, 256, 3e-11)
        require(abs(coarse-refined) < 1e-6 * max(refined, 1.0), 'verify_histogram:3')
        print(f'{speed:6d} {refined:15.10g} {refined/sigma:15.10g}')
        if speed == 1000:
            require(refined > sigma, 'verify_histogram:4')


if __name__ == '__main__':
    main()
