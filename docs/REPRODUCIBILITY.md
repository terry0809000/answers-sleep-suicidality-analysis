# Provenance and reproduction scope

## Source

Packaged from the exact user-requested `ANSWERS_analysis_Rscript.R` (SHA-256 `ebc18012e8b0f7d6b97e833772f5cb8a2baec26e8fd0d3c0e1f8c15732c32fb8`). The original local script was not modified. Other audited scripts in the surrounding workspace were not substituted.

The root `renv.lock` retains package versions, repositories and dependency-resolution fields (descriptive package metadata is omitted) from the recorded `full_20260915_audited_v2` environment (R 4.5.1; 184 package entries). That environment came from a different audited copy of the analysis. It is a useful pinned dependency baseline, **not evidence that this portable copy has completed a fresh full analysis**. Setup must restore it before a real reproduction. `renv` is bootstrapped at version 1.2.4 from CRAN's version archive.

## Packaging changes

1. Removed machine-specific input/output paths and unconditional execution from the function module.
2. Added an explicit CLI with check, preparation and full-run modes. Functions are sourced into an isolated environment; existing globals used by the original functions are supplied there.
3. Changed the output lockfile snapshot's project argument from the raw-data directory to the code-project directory. This prevents environment management inside the data folder.
4. Added an empty-output-directory guard, dependency restoration, synthetic checks and ignored private directories.

No scientific estimator, hypothesis rule, scoring formula, seed, bootstrap total or cohort criterion was intentionally changed. The literal preregistration-question references in source comments are retained from the supplied script; the preregistration itself is not included or independently revalidated here.

## Verification

Run `Rscript --vanilla tests/run_tests.R`. It checks source isolation, frozen settings, independently expected synthetic scoring totals and missingness, signed jet lag, ideation severity versus count, decimal age boundaries, incompatible-input rejection, the six network specifications and script syntax. GitHub Actions repeats this data-free test on R 4.5.1. It does not access research data or fit models.

The full registered analysis is intentionally not run by CI. A successful synthetic test is not a validation of all estimators or reproduction of scientific results. Run `--mode check` and then `--mode full` with authorized data and the restored environment for that purpose. Floating-point results can vary across operating systems, BLAS implementations and compiled libraries even with fixed package versions. Compare numerical estimates, diagnostics and hypothesis decisions rather than assuming byte-identical PDFs or workbooks. The run records session and package information to support that comparison.

## Packaging validation record (22 October 2026)

- Synthetic tests and CLI help passed locally on R 4.5.1.
- An R expression-tree comparison with the requested source confirmed all 46 function definitions are unchanged except the snapshot project argument.
- A clean local dependency restore was attempted. Pure R and available binary packages installed in the staging environment, but compiled archived dependencies failed on this Windows machine, which has no Rtools45 toolchain. The restore did not complete. Setup now reports this prerequisite explicitly. No package versions were relaxed to hide that failure.
- Docker is provided as an alternative build recipe; it has not been built locally.
- No real participant-data analysis, private-output upload or fresh scientific-result reproduction was performed during packaging.
