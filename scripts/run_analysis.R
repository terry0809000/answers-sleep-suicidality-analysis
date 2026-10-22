# Run from any working directory; source functions into an isolated environment.
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."), winslash = "/", mustWork = TRUE)
usage <- function() cat("Usage: Rscript --vanilla scripts/run_analysis.R [--mode check|prepare|full] [--input PATH] [--output-dir PATH]\nDefault: check. Data and all outputs remain local.\n")
if ("--help" %in% args) { usage(); quit(status = 0) }
opts <- list(mode = "check", input = file.path(project, "data/private/answers-dataset-0.2.0.csv"), output = file.path(project, "outputs/private"))
i <- 1L
while (i <= length(args)) {
  key <- args[i]
  if (!key %in% c("--mode", "--input", "--output-dir") || i == length(args)) stop("Unknown or incomplete option: ", key, call. = FALSE)
  target <- switch(key, "--mode" = "mode", "--input" = "input", "--output-dir" = "output")
  opts[[target]] <- args[i + 1L]; i <- i + 2L
}
if (!opts$mode %in% c("check", "prepare", "full")) stop("Mode must be check, prepare or full.", call. = FALSE)
# --vanilla ignores .Rprofile; explicitly select the restored library.
.libPaths(c(file.path(project, ".bootstrap"), .libPaths()))
if (!requireNamespace("renv", quietly = TRUE)) stop("Run Rscript --vanilla scripts/setup.R first.", call. = FALSE)
renv::load(project, quiet = TRUE)
env <- new.env(parent = globalenv())
sys.source(file.path(project, "R/answers_analysis.R"), envir = env)
missing <- env$registered_packages[!vapply(env$registered_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Restore missing dependencies with scripts/setup.R: ", paste(missing, collapse = ", "), call. = FALSE)
if (!file.exists(opts$input)) stop("Authorized raw CSV not found. Supply --input PATH. See docs/DATA.md.", call. = FALSE)
env$PROJECT_DIR <- project
env$INPUT_FILE <- normalizePath(opts$input, winslash = "/", mustWork = TRUE)
env$input_file <- env$INPUT_FILE
env$data_dir <- dirname(env$INPUT_FILE)
env$OUTPUT_DIR <- normalizePath(opts$output, winslash = "/", mustWork = FALSE)
env$analysis_dir <- env$OUTPUT_DIR
env$output_file <- file.path(env$OUTPUT_DIR, "ANSWERS_preregistered_datasets.xlsx")
# Do not mix outputs from different runs or silently overwrite a previous analysis.
if (opts$mode != "check") {
  if (dir.exists(env$OUTPUT_DIR) && length(list.files(env$OUTPUT_DIR, all.files = TRUE, no.. = TRUE))) stop("Output directory must be empty; choose a new --output-dir.", call. = FALSE)
  dir.create(env$OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(env$OUTPUT_DIR, "objects"), showWarnings = FALSE)
}
if (opts$mode == "check") {
  invisible(env$prepare_answers_data(env$INPUT_FILE, write_output = FALSE))
  message("Input schema, released scores and cohort assertions passed. No models fitted or files written.")
} else if (opts$mode == "prepare") {
  invisible(env$prepare_answers_data(env$INPUT_FILE, env$output_file, write_output = TRUE))
  message("Prepared participant-level datasets saved locally: ", env$output_file)
} else {
  invisible(env$run_complete_analysis())
}
