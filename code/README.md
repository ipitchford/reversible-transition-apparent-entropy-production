# Checks for "How a reversible transition amplifies apparent entropy production"

All scripts concern the three-state model of the paper (unit period, theta = 2 pi t, h = 4/5,
q = 1 + (16/17) cos theta, a = 2 pi (1 + h sin 2 theta), b = 2 pi (1 - h sin 2 theta), reporter rate R q).
Each script exits with a nonzero status if a tested inequality fails, also under `python -O`
(no `assert` statements are used). Run from this directory; `python3 -I -B -O <script>.py` works for all.
Expected output: `../results/results.txt`.

## Certificates (Python standard library only; exact integer/rational arithmetic)

| Script | Proves | Paper |
|---|---|---|
| `verify_certificate.py` | S_1 >= 78.336 and S_1/sigma >= 31.156 at R = 1e5 (one edge and all jumps), using pi < 22/7 | Theorem 1 |
| `verify_lag_bins_exact.py` | common lower bound > 21.62 (R = 1e4) and > 366.44 (R = 1e5); ratio > 145.74 | Proposition 2, Appendix C |
| `verify_moderate_speed.py` | validated Fourier-Taylor enclosure: > 2.746942 / 2.755619 at R = 750 and > 3.654098 / 3.661074 at R = 1000 (one edge / all jumps) | Proposition 2, Appendix C |

`verify_moderate_speed.py` imports constants and the logarithm enclosure from `verify_lag_bins_exact.py`; keep both
files together. `lag_bins_proof.md` records the error analysis behind both bin certificates.

## Floating-point illustrations (NumPy 2.3.5, SciPy 1.17.0)

| Script | Computes | Paper |
|---|---|---|
| `verify_numerically.py` | sigma by flux quadrature; race probabilities from periodic backward equations at R = 100, 1000, 1e5; negative controls (undriven, constant reporter); slope 0.0040707 of Theorem 4 | Theorems 1, 4 |
| `verify_histogram.py` | 13-bin direction-blind lower bound on S_1 at R = 100 ... 1e5 (ratio 0.76 at R = 500, 1.15 at 750) | Sec. 3 |
| `verify_crossover.py` | D(J_c(eta) ‖ J_c(xi)) for the critical Gaussian channel; endpoints D_q and 1 - sqrt(17)/5 | Theorem 10 |
| `verify_pumped_current.py` | biased reporter (cycle affinity f): mean current, sigma(f) and a binned lower bound on S_1 at R = 1000, f = 0 ... 1 (about 80 s) | Corollary 3 table |

Floating-point outputs are not enclosures; each script compares two resolutions or tolerances and declares its
tolerance. Values may differ in the last printed digits across platforms.

`../figure.py` (NumPy, Matplotlib) redraws `figure.pdf`.
