"""Numerical check: the overshoot survives a constant cycle affinity (net pumped current).

The three-state model of Theorem 1 is modified by biasing the reporter (twin) edge,
k_{+ -} = R q exp(f/2), k_{- +} = R q exp(-f/2), with every other rate unchanged.
For f != 0 the triangle carries the constant cycle affinity f and a nonzero
period-averaged current. Twin symmetry is broken, so the directed marks are kept.

For each experiment E (forward, and reversed protocol k^B(s) = k^F(-s)) we compute by
killed forward dynamics the event-sampled law of (first mark, next mark, lag bin).
Data processing gives a lower bound on the HFB statistic, Eq. (14):
    S_1 >= K_F [ D(P_bin || Theta Q_bin) - D(nu_F || Theta nu_B) ].
Floating-point computation (illustrative, not a certificate); two phase meshes and
two tolerances are compared and must agree to a declared relative tolerance.
"""

import math
import sys

import numpy as np
from scipy.integrate import solve_ivp

H = 4 / 5
AMP = 16 / 17
TAU = 2 * math.pi
EDGES = np.array([0, 0.1, 0.3, 0.5, 0.8, 1.2, 1.7, 2.5, 3.5, 5, 7, 12, 20])
REL_TOL = 1e-5


def require(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)


def rates(t, R, f, sign):
    """Rate matrix entries at phase t; sign=+1 forward, -1 reversed protocol."""
    th = TAU * sign * np.asarray(t)
    drive = H * np.sin(2 * th)
    a = TAU * (1 + drive)
    b = TAU * (1 - drive)
    q = 1 + AMP * np.cos(th)
    # states 0, + (1), - (2)
    k = {
        (0, 1): a / 2, (0, 2): a / 2,
        (1, 0): b, (2, 0): b,
        (1, 2): R * q * math.exp(f / 2), (2, 1): R * q * math.exp(-f / 2),
    }
    return k


def periodic_state(R, f, sign, nodes):
    def rhs(t, p):
        k = rates(t, R, f, sign)
        dp = np.zeros(3)
        for (x, y), r in k.items():
            dp[x] -= p[x] * r
            dp[y] += p[x] * r
        return dp

    t_eval = 20 + np.arange(nodes) / nodes
    sol = solve_ivp(rhs, (0, 21), np.ones(3) / 3, method="Radau", t_eval=t_eval,
                    rtol=1e-11, atol=1e-13)
    require(sol.success, "periodic state solve")
    return sol.y  # shape (3, nodes), phases arange(nodes)/nodes


def mark_laws(R, f, sign, nodes, tol):
    """Return nu[i], law[i][j, bin] (conditional on first mark i), K, sigma, current."""
    s = np.arange(nodes) / nodes
    pi = periodic_state(R, f, sign, nodes)
    k0 = rates(s, R, f, sign)
    marks = [(1, 2), (2, 1)]
    intens = [pi[x] * k0[(x, y)] for (x, y) in marks]
    K = sum(float(np.mean(r)) for r in intens)
    nu = np.array([float(np.mean(r)) / K for r in intens])
    sigma = 0.0
    for (x, y), r in k0.items():
        sigma += float(np.mean(pi[x] * r * np.log(r / k0[(y, x)])))
    # cycle current 0->+ net flux, averaged
    current = float(np.mean(pi[0] * k0[(0, 1)] - pi[1] * k0[(1, 0)]))
    laws = []
    for idx, (x, y) in enumerate(marks):
        w = intens[idx] / intens[idx].sum()

        def rhs(u, z):
            p = z[:3 * nodes].reshape(3, nodes)
            k = rates(s + u / R, R, f, sign)
            dp = np.zeros_like(p)
            absorbed = np.zeros((2, nodes))
            for (a_, b_), r in k.items():
                flux = p[a_] * r / R
                dp[a_] -= flux
                if (a_, b_) in marks:
                    absorbed[marks.index((a_, b_))] += flux
                else:
                    dp[b_] += flux
            return np.concatenate([dp.ravel(), absorbed.ravel()])

        z0 = np.zeros(5 * nodes)
        z0[y * nodes:(y + 1) * nodes] = 1.0
        sol = solve_ivp(rhs, (0, EDGES[-1]), z0, method="DOP853", t_eval=EDGES,
                        rtol=tol, atol=tol * 0.01)
        require(sol.success, "killed flow solve")
        cum = sol.y[3 * nodes:].reshape(2, nodes, len(EDGES))
        cum = np.einsum("jnk,n->jk", cum, w)  # Palm average
        alive = 1 - cum[:, -1].sum()
        law = np.diff(cum, axis=1)  # (2 next marks, 12 bins)
        law = np.concatenate([law, np.zeros((2, 1))], axis=1)
        # overflow: unresolved beyond 20 lumped with next mark unknown -> put in own cell
        law = np.concatenate([law.ravel(), [alive]])
        require(np.all(law > -1e-12), "negative mass")
        require(abs(law.sum() - 1) < 1e-8, "mass")
        laws.append(np.clip(law, 1e-300, None))
    return nu, laws, K, sigma, current


def lower_bound(R, f, nodes, tol):
    nuF, lawF, KF, sigma, cur = mark_laws(R, f, +1, nodes, tol)
    nuB, lawB, KB, sigmaB, curB = mark_laws(R, f, -1, nodes, tol)
    # P(i, j, bin) = nuF(i) lawF[i][j, bin]; Theta Q(i, j, bin) = nuB(jbar) lawB[jbar][ibar, bin]
    # marks: index 0 = (+,-), 1 = (-,+); bar swaps 0 <-> 1.
    nb = len(EDGES)
    D = 0.0
    for i in range(2):
        P = nuF[i] * lawF[i]
        for j in range(2):
            ib, jb = 1 - i, 1 - j
            p = P[j * nb:(j + 1) * nb]
            q = nuB[jb] * lawB[jb][ib * nb:(ib + 1) * nb]
            D += float(np.sum(p * np.log(p / q)))
    # unresolved cell (no mark within u <= 20): compare forward and backward totals
    pu = sum(nuF[i] * lawF[i][-1] for i in range(2))
    qu = sum(nuB[i] * lawB[i][-1] for i in range(2))
    D += pu * math.log(pu / qu)
    Dnu = float(np.sum(nuF * np.log(nuF / nuB[::-1])))
    return KF * (D - Dnu), sigma, cur, KF, KB


def main():
    R = 1000
    print("Reporter bias f, net current J (0->+ per period), sigma, S_1 lower bound (13 lag bins x marks), ratio")
    for f in (0.0, 0.25, 0.5, 1.0):
        coarse = lower_bound(R, f, 64, 1e-9)
        fine = lower_bound(R, f, 128, 1e-11)
        require(abs(coarse[0] - fine[0]) < REL_TOL * abs(fine[0]) + 1e-8, f"mesh agreement at f={f}")
        S, sigma, cur, KF, KB = fine
        print(f"{f:5.2f} {cur: .6f} {sigma:.6f} {S:.6f} {S / sigma:.4f}")
        if f > 0:
            require(cur > 1e-3, "net current should be nonzero")
            require(S > sigma, f"overshoot at f={f}")
    print("PASS")


if __name__ == "__main__":
    sys.exit(main())
