# Input contract and privacy

Obtain authorized access to `answers-dataset-0.2.0.csv` through the research project's approved process. No access URL or redistribution permission is inferred here. Do not upload the CSV, participant-level spreadsheets, saved R objects, fold manifests, predictions or credentials to GitHub.

The original raw-release checks require 971 rows, 294 columns, unique nonmissing `id`, unique `(id, visit)` and `visit == 1` throughout. The nine-column harmonized dataset is insufficient. Required item names and the released-summary comparisons are specified directly by `read_and_check_raw_data`, `score_registered_measures`, `validate_registered_scores` and `assemble_prepared_data` in `R/answers_analysis.R`.

The primary cohort is age >=18 and <26 (885 participants in the expected release); the age sensitivity cohort includes all observed ages (945). Further event-count assertions deliberately reject a mismatched release. Do not remove assertions to make another dataset pass.

Scoring requires every constituent item: there is no imputation or prorating. CES-D item 11 is excluded; items 4, 8, 12 and 16 are reversed. Signed social jet lag is converted from minutes to hours without taking its absolute value. The C-SSRS highest endorsed category is distinct from the count of endorsed items. These rules are covered by synthetic tests.

Keep authorized data outside the repository, or in the ignored `data/private/` directory. Output defaults to ignored `outputs/private/`. An external `--output-dir` is also supported; it remains your responsibility to secure it. Git ignore rules are an accident-prevention measure, not access control.
