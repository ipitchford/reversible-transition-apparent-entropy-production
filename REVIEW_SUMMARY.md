# Review summary (internal, AI-assisted)

This file records the internal checking of this candidate before release, without any ratings. It is not peer review,
independent reproduction or formal verification. All review steps were carried out by AI systems (Claude) under the
Evidence Press publication protocol; no external referee has seen the paper.

## What was done (24–25 November 2026)

1. **Five reviewer roles** (editor, domain, perspective, methodology and devil's advocate), run as Claude subagents.
   The methodology role audited the proofs at argument level and reran the scripts. No role found an error in a proved
   statement; the objections concerned the checks under `python -O`, engagement with Harunari, Fiore and Barato's
   zero-current remark, the absence of a current-carrying case, the location of the certificate behind the factor 145,
   the novelty of the block inequality relative to their Eq. (26), the size of the effect, practical guidance, the
   framing of the earlier Evidence Press release, the scope of the scaling table and readability.
2. **Replay** of every script by the lead from a fresh copy, which found that the pass/fail checks relied on `assert`
   (and so would be skipped under `python -O`). All checks were converted to explicit failures. Three negative controls
   (R = 1000 in `verify_certificate.py`, near-flat weights in `verify_lag_bins_exact.py`, flat weights in
   `verify_moderate_speed.py`) each exit 1 under `-O`.
3. **One revision round** with extensions (listed below).
4. **Re-review and proof audit** of the new and rewritten proofs and of `verify_pumped_current.py`. The audit found one
   correctness item, a wrong threshold in Corollary 3 (the bound used there exceeds σ only for R above about 5282, not
   for every R above about 1216); the hypothesis now reads "for example any R ≥ 5400", checked in exact arithmetic with
   π < 22/7 (bound 2.60 at R = 5400, above 88/35). It also found two wording overclaims (the net current under a cycle
   affinity is numerical, not proved) and six minor wording and number items. All were applied and spot-checked against
   the rebuilt PDF.
5. No Astra/Codex review took place for this candidate.

## Main objections and how they were resolved

- **Checks under `python -O`.** Every `assert` replaced by an explicit failure; negative controls recorded.
- **Harunari, Fiore and Barato's zero-current remark.** The introduction quotes their Section VI, and Section 7
  discusses it; Corollary 6 identifies the phase asymmetry X as the quantity that controls the excess.
- **Current-carrying case (extension).** Corollary 3 proves that the violation persists under a small constant cycle
  affinity (and exceeds 30σ at R = 10^5); a floating-point table at R = 1000 (affinity up to 1, ratio 1.10 to 1.53) is
  produced by the new `verify_pumped_current.py`. The net current is reported as numerical.
- **The factor 145.** Proposition 2 and Appendix C put the lag-bin certificate and the moderate-speed enclosure in the
  paper.
- **Novelty of the block inequality.** Theorem 5 states that m = 1 is their Eq. (26); the new content is the extension
  to blocks, the safety criterion (Corollary 6, an extension) and the sharpness of the coefficient (Theorem 7).
- **Size of the effect.** New paragraphs give the excess in nats per event and a quasi-static reading; the smallest
  certified violation (R = 750) is stated.
- **Practical guidance and limits.** Section 7 states consequences for practice (X needs only one-dimensional phase
  histograms; plug-in bias; conservative variational estimators) and limits (no finite-sample analysis).
- **The earlier Evidence Press release.** Its list of unaffected estimators remains true of the phase-resolved
  estimator and of the corrected inequality, not of the uncorrected estimate; no wording implies that the corrected
  inequality fails.
- **Scope.** The scaling table and the sharpness statement are labelled as proved for the three-state example only.
- **Readability and literature.** Context-first introduction with a roadmap, purpose sentences for every section, notation
  defined before use; five further references added from Crossref records.

## What remains open or unchecked

- Whether the uncorrected estimate can exceed σ in Harunari, Fiore and Barato's own pump model (not scanned).
- Finite-sample estimation (bias and variance) is not analysed.
- A general mechanism theorem beyond the example; the general R^(2/3) expectation is not proved.
- Prior-art searches were of limited scope; this is not a proof of priority.
- No external specialist has checked the proofs.
