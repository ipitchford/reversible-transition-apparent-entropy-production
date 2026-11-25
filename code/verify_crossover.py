"""Quadrature for the limiting phase experiment, independent of jump-path numerics."""
import numpy as np
from numpy.polynomial.hermite import hermgauss


def require(condition, label):
    """Fail loudly, also under python -O (unlike assert)."""
    if not condition:
        raise SystemExit(f"CHECK FAILED: {label}")


def divergence(scale, n_phase, n_normal):
    theta = np.arange(n_phase) * (2 * np.pi / n_phase)
    rate = 1 + (16 / 17) * np.cos(theta)
    rate_derivative = -(32 * np.pi / 17) * np.sin(theta)
    shift = scale**1.5 * rate_derivative / (np.sqrt(12) * rate**2)

    def phase_priors(angle):
        common = 1 + (16 / 17) * np.cos(angle)
        sine, cosine = np.sin(2 * angle), np.cos(2 * angle)
        return (common * (1 + 0.4 * (sine - cosine)),
                common * (1 + 0.4 * (sine + cosine)))

    eta, xi = phase_priors(theta)
    eta_other, xi_other = phase_priors(-theta)
    nodes, weights = hermgauss(n_normal)
    terms = []
    for node, weight in zip(nodes, weights):
        z = shift + np.sqrt(2) * node
        other_log_weight = -2 * z * shift
        log_num = np.logaddexp(np.log(eta), np.log(eta_other) + other_log_weight)
        log_den = np.logaddexp(np.log(xi), np.log(xi_other) + other_log_weight)
        terms.append(weight * np.mean(eta * (log_num - log_den)))
    return np.sum(terms) / np.sqrt(np.pi)


if __name__ == "__main__":
    d_rate = 2 * (1 - np.sqrt(21) / 5)
    d_phase = 1 - np.sqrt(17) / 5
    print("c, D(J_c(eta) || J_c(xi)), refinement difference")
    previous = 0
    for c in (0, 0.25, 0.5, 1, 2, 4, 8):
        coarse = divergence(c, 2048, 64)
        fine = divergence(c, 8192, 128)
        print(f"{c:4.2f}  {fine:.10f}  {abs(fine-coarse):.3g}")
        require(previous < fine <= d_phase + 1e-11, 'verify_crossover:1')
        require(abs(coarse - fine) < 1e-5, 'verify_crossover:2')
        if c == 0:
            require(abs(d_rate - fine) < 1e-12, 'verify_crossover:3')
        previous = fine
    print("Exact rate-only and full-phase endpoint values:", d_rate, d_phase)
