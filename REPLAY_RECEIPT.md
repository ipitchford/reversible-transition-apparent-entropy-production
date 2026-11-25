# Replay receipt (release 0.1.0-candidate, 25 November 2026)

Input: the approved candidate (manuscript.tex with five `\input` files, figure.pdf, figure.py, manuscript.pdf, nine
files in code/ and results/results.txt), each SHA-256 checked against the frozen candidate list before use. After the
replay only the title-page DOI line was added to manuscript.tex and the PDF was rebuilt (pdflatex, three passes, 20
pages, as approved; word-level text comparison with the approved PDF shows only the added DOI line and a page number
moving within the extracted text). No script or recorded output changed.

## Runs (fresh copies, cwd = code/, macOS arm64)

| Command | Interpreter | Result |
|---|---|---|
| `python3 -I -B -O verify_certificate.py` | CPython 3.14.7 and 3.12.11 | exit 0; output identical to results/results.txt (S_1 >= 78.3368911365, ratio >= 31.1567180657) |
| `python3 -I -B -O verify_lag_bins_exact.py` | 3.14.7 and 3.12.11 | exit 0; identical (ratio within [14574516/100000, 14574517/100000] at R = 1e5) |
| `python3 -I -B -O verify_moderate_speed.py` | 3.14.7 and 3.12.11 | exit 0; identical (2.746942 and 2.755619 at R = 750) |
| `python3 -I -B -O verify_numerically.py` | 3.12.11, NumPy 2.3.5, SciPy 1.17.0 | exit 0; identical to the recorded output |
| `python3 -I -B -O verify_histogram.py` | same | exit 0; identical |
| `python3 -I -B -O verify_crossover.py` | same | exit 0; identical |
| `python3 -I -B -O verify_pumped_current.py` | same | exit 0; PASS; identical (about 65 s) |

## Negative controls (fresh copies, `python3 -I -B -O`)

- `reporter_speed = 100_000` changed to `1000` in verify_certificate.py: `CHECK FAILED: verify_certificate:5`, exit 1.
- Lag-bin weights in verify_lag_bins_exact.py changed to near-flat (99, 99, 100, 100, 101, 101, 100, 99, 99)/100:
  `CHECK FAILED: verify_lag_bins_exact:6`, exit 1.
- All weights in verify_moderate_speed.py changed to 1: `CHECK FAILED: verify_moderate_speed:3`, exit 1.

## Comparison policy

As in ENVIRONMENT.txt and README.md: the three certificates compare exactly with results/results.txt on every
platform; the four floating-point scripts compare token by token with numbers within a relative 1e-6 or an absolute
1e-9; interpreter versions are not compared. `.github/ci/compare_results.py` implements it.
