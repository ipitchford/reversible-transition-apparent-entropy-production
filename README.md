# How a reversible transition amplifies apparent entropy production

**A counterexample to the uncorrected waiting-time inequality of Harunari, Fiore and Barato, and the sharp correction
for blocks of waits.** Version 0.1.0-candidate. DOI [10.5281/zenodo.24529102](https://doi.org/10.5281/zenodo.24529102)
(all versions: [10.5281/zenodo.24529101](https://doi.org/10.5281/zenodo.24529101)). Release page:
https://evidencepress.org/releases/reversible-transition-apparent-entropy-production/

Unrefereed research candidate. Author: Anonymous. Text CC0-1.0; code MIT (see [LICENSE.md](LICENSE.md)).

## What this is about

Small driven systems, such as molecular motors, ion channels and colloidal beads in moving laser traps, turn part of
the energy they receive into heat. The rate at which they do so is their *entropy production*. Experiments usually see
only a few of the system's transitions and the waiting times between them, so physicists derive *lower bounds* on
dissipation from those waiting times: if the record looks different when played backwards, the system must be
dissipating at least a computable amount.

For systems driven by a repeating rhythm, Harunari, Fiore and Barato (Physical Review E, 2024) proposed a practical
estimator, here called S_1: run the experiment once with the rhythm and once with the rhythm reversed, compare the two
histograms of waiting times between visible jumps, and never record the phase of the rhythm. They proved that the
estimate minus a correction X is a valid lower bound, and asked whether the estimate alone is one.

An earlier Evidence Press release,
[Counterexamples to a conjectured waiting-time bound on entropy production in periodically driven Markov networks](https://doi.org/10.5281/zenodo.23541298),
refuted a related forward-only conjecture and gave numerical evidence on this question in one special case. It listed
estimators that use the reversed experiment as unaffected. That remains true of the phase-resolved estimator and of
Harunari, Fiore and Barato's corrected inequality, which stays valid; this release shows that it is not true of the
uncorrected estimate.

## What is new

1. **The uncorrected estimate is not a lower bound, and its overshoot is unbounded** (Theorem 1, Proposition 2). A
   three-state model whose currents all average to zero, with equal event rates in both experiments and fixed entropy
   production 4π/5 (about 2.5), has an estimate that grows in proportion to the speed R of one reversible transition.
   Exact rational arithmetic proves an estimate more than 31 times the true value at one speed, and a certified
   lag-bin bound more than 145 times. The smallest certified violation is at about 375 reporter events per period.
2. **Robustness.** Recording every microscopic jump does not help. A small constant force (cycle affinity) does not
   help either (Corollary 3); under it the model numerically pumps a net current, and the overshoot persists up to the
   largest force tested. The net current is a numerical observation, not a proved statement.
3. **The sharp correction for blocks** (Theorems 5 and 7, Corollary 6). For any finite driven network, using blocks of
   m consecutive waits, the overshoot is at most X/m, where X needs only the phase histogram of each kind of visible
   event, and this correction cannot be made smaller. At m = 1 this is Harunari, Fiore and Barato's corrected
   inequality. If events sample the phase symmetrically (X = 0), the uncorrected estimate is safe.
4. **How many waits recover the dissipation** (Section 6, appendices). In the example, the block estimate is
   determined in every scaling regime, with a transition near R^(2/3) waits, where the waits begin to reveal whether
   the drive is rising or falling.

**The idea.** A fast transition that produces no entropy (the same rate both ways) still reports which part of the
cycle the system is in, because its speed varies with the phase. Reversing the drive changes which speeds a typical
wait samples, even though the average number of events is the same. A single wait sees a mixture of exponential
waiting times, the two experiments mix them differently, and the estimator multiplies that small difference by a large
number of events.

## Limits

- A model calculation; no claim about a particular experiment.
- The overshoot per event is small (under 0.008 nats per event at the largest speed certified) and needs a reporter
  several hundred times faster than the drive.
- The scaling table and the sharpness statement are proved for the three-state example only; the block inequality
  holds for every finite network with continuous periodic rates bounded away from zero.
- Whether S_1 ≤ σ holds throughout Harunari, Fiore and Barato's own pump model remains open. Finite-sample estimation
  is not analysed.
- Practical advice: do not report the uncorrected estimate as a bound unless the visible events sample the phase
  symmetrically; use the block estimate minus X/m, or the phase-resolved estimator.

## Contents

| Path | What it is |
|---|---|
| [manuscript.pdf](manuscript.pdf) | The paper (20 pages). Source: [manuscript.tex](manuscript.tex) with the `\input` files `shortblocks.tex`, `critical.tex`, `macro.tex`, `longblocks.tex` and `bins.tex`; figure [figure.pdf](figure.pdf) drawn by [figure.py](figure.py) |
| [code/](code/) | Three certificates in exact integer and rational arithmetic and four floating-point illustrations; see [code/README.md](code/README.md) and [code/lag_bins_proof.md](code/lag_bins_proof.md) |
| [results/results.txt](results/results.txt) | Recorded output of every script |
| [AI_INDEX.md](AI_INDEX.md) | Machine-oriented map of claims, evidence and limits |
| [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md) | Internal AI-assisted review record (rating-free) |
| [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), [ENVIRONMENT.txt](ENVIRONMENT.txt) | Release replay and environment |
| [RESEARCH_GATES.json](RESEARCH_GATES.json) | Prior-art, source-reading and limitation record |
| [refs.bib](refs.bib), [CITATION.cff](CITATION.cff), [LICENSE.md](LICENSE.md), [MANIFEST.sha256](MANIFEST.sha256) | Bibliography, citation metadata, licences, checksums |

## Reproducing the checks

From `code/`, with Python 3.12 and the pinned packages in [requirements.txt](requirements.txt):

```
python3 -I -B -O verify_certificate.py      # standard library; exact rational
python3 -I -B -O verify_lag_bins_exact.py   # standard library; exact rational
python3 -I -B -O verify_moderate_speed.py   # standard library; integer intervals
python3 -I -B -O verify_numerically.py      # NumPy/SciPy
python3 -I -B -O verify_histogram.py
python3 -I -B -O verify_crossover.py
python3 -I -B -O verify_pumped_current.py   # about 80 seconds
```

Every check exits with a nonzero status if a tested inequality fails, also under `python -O` (no `assert` is used).

**Comparison policy** (declared before tagging; the same in CI, [ENVIRONMENT.txt](ENVIRONMENT.txt) and
[REPLAY_RECEIPT.md](REPLAY_RECEIPT.md)): the output of the three certificates (verify_certificate,
verify_lag_bins_exact, verify_moderate_speed) must equal its section of `results/results.txt` exactly. The four
floating-point scripts must match token by token, non-numeric tokens exactly and numbers within a relative 1e-6 or an
absolute 1e-9, because solver output may differ in the last printed digits across platforms; each script also enforces
its own declared tolerance. Interpreter versions are not compared.
[.github/ci/compare_results.py](.github/ci/compare_results.py) implements the comparison.

## Status and assurance

Written proofs with exact rational certificates. Internal AI-assisted review only: a five-role Claude review panel, a
replay of the exact rational certificates with negative controls, one revision round with extensions, then a re-review
and proof audit of the new proofs (which corrected the threshold in Corollary 3 to R ≥ 5400, checked exactly). No
Astra/Codex review took place. This is not peer review, independent reproduction or formal verification.

## Cite

Anonymous (2026). *How a reversible transition amplifies apparent entropy production*. Version 0.1.0-candidate.
Evidence Press. https://doi.org/10.5281/zenodo.24529102
