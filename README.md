![ANSWERS — text-free conceptual cover illustrating sleep disturbance and psychological distress](assets/answers-cover.jpg)

# ANSWERS: sleep, psychological vulnerability and suicidal ideation

A portable research-analysis project built from `ANSWERS_analysis_Rscript.R`. It reconstructs registered scores, validates the expected release and cohorts, estimates the primary and sensitivity networks, and runs the separate classifier and rare-outcome descriptions. This is a runnable analysis bundle, not an installable CRAN package or clinical decision tool.

**Participant data are not distributed.** A complete reproduction requires authorized access to the exact raw release. The repository contains code and dependency metadata only; synthetic tests do not produce research findings.

## Quick start

Install **R 4.5.1**. Clone this repository, open a terminal in its directory, then run:

```sh
Rscript --vanilla scripts/setup.R
Rscript --vanilla tests/run_tests.R
Rscript --vanilla scripts/run_analysis.R --help
Rscript --vanilla scripts/run_analysis.R --mode check --input "/authorized/path/answers-dataset-0.2.0.csv"
Rscript --vanilla scripts/run_analysis.R --mode full --input "/authorized/path/answers-dataset-0.2.0.csv" --output-dir "outputs/private/run-001"
```

On Windows, quote the full path to `Rscript.exe` if R is not on PATH. Use PowerShell's `&` before a quoted executable path. `--mode prepare` writes only the two prepared participant-level datasets. `check` validates the raw release, reconstructed scores and cohort assertions without fitting models or writing data. The default mode is `check`. Paths passed on the command line are relative to the terminal's current directory; default paths are relative to the project root.

Setup restores the **existing recorded lockfile**, including transitive dependencies, into a project library. It does not resolve a fresh set of latest package versions. Internet access is required for installation. On Windows, install [Rtools45](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html) before setup: some pinned archived dependencies require compilation. Linux/macOS source builds also require appropriate compilers and development libraries. See [reproduction notes](docs/REPRODUCIBILITY.md).

## Container option

A Dockerfile supplies R 4.5.1 and a Linux compilation toolchain. Docker is optional; it is not installed by the setup script. From the repository root:

```sh
docker build -t answers-analysis .
docker run --rm answers-analysis --help
```

For a real run, mount an authorized input directory read-only at `/data` and a new empty local output directory at `/results`, then pass `--mode full --input /data/answers-dataset-0.2.0.csv --output-dir /results`. Do not copy research data into the image. `.dockerignore` excludes everything except the analysis, tests, setup and lockfile. The container build and full analysis were not executed during packaging; see the verification scope below.

## Analysis and outputs

- N1 is the primary network for H1a閳ユ弻1e and H2; N2閳ユ彉5 and the age sensitivity analysis retain their original roles.
- P1 is the separate classification analysis for H3; D1 describes rarer suicide-related outcomes.
- The original seed, single-core execution, scoring, cohort rules, model functions and resampling totals are retained. A full run performs the registered 1,000-draw stages and can take hours.
- A new, empty output directory is required for writing modes. The full run writes prepared datasets, result workbooks, figure PDFs, a private RDS object, session information and a run-specific lockfile.
- Prepared data and the RDS contain participant-level information, including IDs and predictions. Keep **all outputs private** until an appropriate disclosure review; ignored files must never be force-added.

Cross-sectional network edges describe conditional associations and do not establish causation or temporal ordering. Classification results are research estimates and are not a clinical deployment claim.

## Repository guide

| Path | Purpose |
|---|---|
| `R/answers_analysis.R` | Original analysis functions and frozen settings; sourcing does not run models |
| `scripts/setup.R` | Bootstrap pinned renv and restore recorded dependencies |
| `scripts/run_analysis.R` | Portable command-line entry point |
| `renv.lock` | Dependency versions from an existing local analysis record |
| `tests/run_tests.R` | Small synthetic regression checks requiring only base R |
| `docs/DATA.md` | Expected input and data-access boundary |
| `docs/REPRODUCIBILITY.md` | Provenance, changes and validation scope |

## Rights and attribution

No new software license is assigned by this packaging work. Obtain the relevant permissions before redistributing or reusing the code or data. Package dependencies retain their own licenses. This repository does not distribute the underlying preregistration or assert a public dataset license.

Cover: AI-generated conceptual raster artwork made with the built-in imagegen tool; it depicts no measured data or research findings. The generation prompt is recorded in `docs/COVER.md`.
