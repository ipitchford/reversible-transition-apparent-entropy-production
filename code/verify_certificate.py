"""Exact rational lower bounds for the three-state, unit-period example."""

from fractions import Fraction as F


def require(condition, label):
    """Fail loudly, also under python -O (unlike assert)."""
    if not condition:
        raise SystemExit(f"CHECK FAILED: {label}")


def main():
    contrast = F(4, 5)
    amplitude = F(16, 17)
    reporter_speed = 100_000
    pi_upper = F(22, 7)
    # Here sqrt(1-h^2) = 3/5 and sqrt(4-k^2) = 30/17.
    root_h, root_k = F(3, 5), F(30, 17)
    require(root_h**2 == 1 - contrast**2, 'verify_certificate:1')
    require(root_k**2 == 4 - amplitude**2, 'verify_certificate:2')
    poisson_parameter = (2 - root_k) / amplitude
    require(poisson_parameter == F(1, 4), 'verify_certificate:3')
    race_gap = contrast * poisson_parameter**2 / root_k
    require(race_gap == F(17, 600), 'verify_certificate:4')

    # Couple each physical experiment to its frozen-rate reporter until
    # the Exp(R) deadline. One slow exit or failure of the rate coupling
    # costs at most [max(b) + Lip(q)]/R, in each experiment separately.
    two_arm_error = 4 * pi_upper * (1 + contrast + amplitude) / reporter_speed
    remaining_gap = race_gap - two_arm_error
    require(remaining_gap > 0, 'verify_certificate:5')
    active_probability = F(1, 2)
    visible_rate = reporter_speed * active_probability
    lower = 2 * visible_rate * max(remaining_gap, F(0))**2  # explicit positive part
    entropy_upper = 2 * pi_upper * (1 - root_h)
    require(entropy_upper == F(88, 35), 'verify_certificate:6')
    require(lower > 30 * entropy_upper, 'verify_certificate:7')

    # The same test on adjacent reporter events in a fully observed jump
    # record is diluted by the fraction of all jumps that are reporters.
    slow_activity_upper = 2 * pi_upper * (1 - contrast**2 / 4)
    full_lower = lower * visible_rate / (visible_rate + slow_activity_upper)
    require(full_lower > 30 * entropy_upper, 'verify_certificate:8')
    for label, value in [
        ('limiting Bernoulli probability gap', race_gap),
        ('two-arm finite-speed error upper bound', two_arm_error),
        ('finite-speed gap lower bound', remaining_gap),
        ('entropy production upper bound', entropy_upper),
        ('one observed edge, estimator lower bound', lower),
        ('one observed edge, ratio lower bound', lower / entropy_upper),
        ('all edges observed, estimator lower bound', full_lower),
        ('all edges observed, ratio lower bound', full_lower / entropy_upper),
    ]:
        print(f'{label}: {value} = {float(value):.12g}')


if __name__ == '__main__':
    main()
