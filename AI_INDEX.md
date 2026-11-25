# AI index: How a reversible transition amplifies apparent entropy production (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) ("What
this is about") and the paper [manuscript.pdf](manuscript.pdf) (source [manuscript.tex](manuscript.tex) with its
`\input` files; Section 1 states the problem, why it matters, what was known, what is new and a roadmap). DOI
10.5281/zenodo.24529102 (concept 10.5281/zenodo.24529101). Status: unrefereed candidate; internal AI-assisted review
only. Follow-up to the Evidence Press release doi:10.5281/zenodo.23541298 (waiting-time-entropy-bound-counterexamples).

## Normalisation

Continuous-time Markov jump process on a finite connected graph with bidirected edges; rates k_xy(s) continuous,
T-periodic, bounded and bounded away from zero. Forward experiment F uses k; backward experiment B uses k^B(s) = k(-s)
(HFB Eq. (13)). sigma = physical entropy production rate (Eq. (1), HFB Eq. (2)), natural logarithms, units of k_B.
Visible marks: a reversal-closed set of directed edges; K_E = visible activity; mu_E = event-sampled (Palm)
phase-and-mark law; nu_E = mark marginal. S_1 = HFB's phase-blind estimator (HFB Eq. (14); their sigma-hat), Eq. (2).
sigma_L = phase-resolved rate (HFB Eq. (24)), sigma_L <= sigma. S_m = block estimator, Eq. (8); X = HFB's correction
(HFB Eq. (25)), Eq. (9). Example: three states {0, +, -}, unit period, theta = 2 pi t, h = 4/5,
a = 2 pi (1 + h sin 2theta), b = 2 pi (1 - h sin 2theta), q = 1 + (16/17) cos theta; reporter edge +<->- at rate R q.
HFB = Harunari, Fiore and Barato, Phys. Rev. E 110, 064126 (2024), doi:10.1103/PhysRevE.110.064126. Reference keys as
in [refs.bib](refs.bib).

## Claims and where they are proved

| ID | Claim | Location | Status / credit |
|---|---|---|---|
| 1 | Example: zero period-averaged current on every edge, frozen detailed balance, K_F = K_B = R/2, sigma = 4 pi/5; S_1 >= R [17/600 - 2C/R]_+^2 with C = 18 pi/5 + 32 pi/17; S_1 > 78 and S_1/sigma > 31 at R = 1e5; S_1/sigma unbounded | Theorem 1 | written proof; exact rational certificate (code/verify_certificate.py; pi < 22/7) |
| 2 | Certified lower bounds on S_1: R = 750: 2.74694 (one edge) / 2.75561 (every jump); R = 1000: 3.65409 / 3.66107; R = 1e4: 21.62; R = 1e5: 366.44, so S_1/sigma > 145 | Proposition 2, Appendix C | rational lag-bin certificate (verify_lag_bins_exact.py) and validated Fourier-Taylor enclosure (verify_moderate_speed.py); error analysis in code/lag_bins_proof.md |
| 3 | Constant cycle affinity f on the reporter: for R with R[17/600 - 2C/R]^2 > 4 pi/5 (e.g. any R >= 5400) there is f_0 > 0 with S_1(f) > sigma(f) for abs(f) < f_0; at R = 1e5 also S_1(f) > 30 sigma(f); also at R = 750 and 1000 via Proposition 2 | Corollary 3 | written proof (continuity); threshold 5400 checked in exact arithmetic (2.60 at R = 5400). The nonzero net current for f != 0 and the table at R = 1000 (f up to 1, ratio 1.10 to 1.53) are numerical (verify_pumped_current.py), not proved |
| 4 | High-activity slope: S_1/R -> (1/2) D(g_F ‖ g_B) = 0.00407070 (numerical quadrature) | Theorem 4 | written proof; slope numerical (verify_numerically.py) |
| 5 | Complete observation of all six directed jumps: S_1^all >= w R [17/600 - 2C/R]_+^2, > 78 at R = 1e5 | Section 4 | written proof |
| 6 | Block correction for every finite network: 0 <= S_m <= sigma_L + X/m <= sigma + X/m; m = 1 is HFB Eq. (26) | Theorem 5 | written proof; m = 1 credited to HFB |
| 7 | Safety criterion: S_1 - sigma <= X; if mu_F = Theta mu_B then X = 0 and S_m <= sigma for all m | Corollary 6 | written proof |
| 8 | In the example X/K_F = 1 - sqrt(17)/5 (X/R about 0.088); S_m >= (R/m)[17/600 - 2C/R]_+^2 | Eqs. (12)-(13) | written proof |
| 9 | Sharpness: in the example with m = 2k, R^(2/3)/k -> 0, k/R -> 0, m S_m / X -> 1; no universal coefficient below one | Theorem 7 | written proof; proved for the example only |
| 10 | Scaling regimes of S_m in the example (single waits to blocks spanning many periods; transition at m of order R^(2/3)) | Section 6 table, Lemma 8, Corollaries 9 and 12, Theorems 10 and 11, Appendices A-B | written proofs, example only; the general R^(2/3) expectation is stated as unproved; critical-channel values checked numerically (verify_crossover.py) |
| Open | Whether S_1 <= sigma throughout HFB's own three-state pump model; finite-sample estimation; general mechanism beyond the example | Section 7 | open |

## Evidence

| Check | Script | Scope |
|---|---|---|
| Theorem 1 finite inequalities at R = 1e5 (one edge and all jumps) | [code/verify_certificate.py](code/verify_certificate.py) | exact rational, standard library |
| Proposition 2 rows R = 1e4, 1e5 (nine lag bins) | [code/verify_lag_bins_exact.py](code/verify_lag_bins_exact.py) | exact rational, standard library |
| Proposition 2 rows R = 750, 1000 | [code/verify_moderate_speed.py](code/verify_moderate_speed.py) | integer-interval enclosure, standard library |
| sigma by quadrature; race probabilities; negative controls (undriven, constant reporter); Theorem 4 slope | [code/verify_numerically.py](code/verify_numerically.py) | floating point (NumPy/SciPy) |
| 13-bin direction-blind bound (ratio 0.76 at R = 500, 1.15 at R = 750) | [code/verify_histogram.py](code/verify_histogram.py) | floating point |
| Critical Gaussian channel (Theorem 10), endpoints | [code/verify_crossover.py](code/verify_crossover.py) | floating point |
| Corollary 3 table (biased reporter, R = 1000) | [code/verify_pumped_current.py](code/verify_pumped_current.py) | floating point; two phase meshes agree to relative 1e-5 |

Recorded output: [results/results.txt](results/results.txt). Every check exits non-zero on failure without `assert`,
also under `python -O`. Negative controls (release replay and CI): R = 1000 in verify_certificate.py, near-flat
weights in verify_lag_bins_exact.py and flat weights in verify_moderate_speed.py each exit 1 under -O.

Numbers quoted on the release page and in the audio: 4 pi/5 (about 2.5); more than 31 times; more than 145 times;
375 reporter events per period; under 0.008 nats of overshoot per event; ratio 1.10 to 1.53 in the pumped-current
table; R^(2/3). All come from the paper and results/results.txt.

## Limits

Model calculation, no experimental claim. The overshoot per event is small and needs a reporter several hundred times
faster than the drive. Scaling table and sharpness proved for the example only. The net current under a cycle affinity
is numerical. Finite-sample estimation is not analysed. Prior-art searches were of limited scope (see
[RESEARCH_GATES.json](RESEARCH_GATES.json)).

## Records

[REVIEW_SUMMARY.md](REVIEW_SUMMARY.md) (rating-free internal review record), [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md),
[ENVIRONMENT.txt](ENVIRONMENT.txt), [RESEARCH_GATES.json](RESEARCH_GATES.json), [refs.bib](refs.bib),
[CITATION.cff](CITATION.cff), [LICENSE.md](LICENSE.md), [MANIFEST.sha256](MANIFEST.sha256),
[.github/workflows/replay.yml](.github/workflows/replay.yml).
