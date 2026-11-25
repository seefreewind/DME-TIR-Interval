# Release provenance

Version 1.0.0 packages the previously completed formal analysis and post hoc extensions. Scientific estimands, populations, calendars, seeds, interval and quantile algorithms are unchanged. Historical configuration files and pre-access addenda retain their original dates and authorization fields; those fields document the sequence of analysis, not the current code-publication authorization.

The folder names `phase2_ajo` and `phase2_ajo_revision` are retained solely to preserve the historical relative import/file structure. The public repository and software title are journal-neutral. No manuscript, submission package, individual record or original dataset is included.

## Engineering adaptation

`phase1/scripts/prepare.py` still verifies both ZIP SHA256 values, freshly extracts and byte-verifies each member, then regenerates the original normalized tables, eligibility and masks. Comparison against unavailable private historical tables is replaced by comparison against their full-file SHA256 values in `provenance/frozen_input_sha256.json`. No rows or scientific rules are changed. The standalone runner replaces manuscript/report-generation orchestration with code-only analysis, validation and a final digest comparison. Analysis engines are otherwise copied unchanged.

`provenance/final_table_sha256.json` records original final TSV digests, including private tables without redistributing their contents. The release checker compares all generated files for which a digest exists; files intentionally not regenerated are reported separately. Floating-point numerical variation across platforms can affect byte digests; the reproducibility report identifies failures instead of silently relaxing the criterion. The original reference environment uses NumPy 2.5.5.

## Interpretation

The four-week schedule is a denser reference, not daily visual-acuity truth. Record thinning under existing treatment paths is not an intervention on treatment management. Protocol AC fellow eyes are resampled together by participant. Post hoc extensions are labeled separately. AC selection weighting failed the prespecified qualification gate and remains diagnostic only; it must not be described as a qualified population correction.

The extension audit retains original independent cluster-expansion, model, weighting and CI checks. Two private historical-manifest loops are replaced by the public runner's final digest comparison; the audit's summary reports this substitution.

The public runner additionally checks NumPy/Matplotlib imports before starting a real-data run. Source tables in the synthetic fixture and public field mappings use LF line endings; their values are unchanged.
