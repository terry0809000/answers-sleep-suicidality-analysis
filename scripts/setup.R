# Bootstrap the version recorded in renv.lock, then restore the complete library.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."), winslash = "/", mustWork = TRUE)
if (getRversion() != "4.5.1") stop("Use R 4.5.1 to match the recorded analysis environment.", call. = FALSE)
if (.Platform$OS.type == "windows" && !nzchar(Sys.which("make")) && !dir.exists("C:/rtools45")) {
  stop("Rtools45 is required to compile archived pinned packages. Install it from https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html then rerun setup. Alternatively use the Dockerfile.", call. = FALSE)
}
options(repos = c(CRAN = "https://cloud.r-project.org"))
bootstrap <- file.path(project, ".bootstrap")
dir.create(bootstrap, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(bootstrap, .libPaths()))
renv_version <- "1.2.4"
if (!requireNamespace("renv", quietly = TRUE) || as.character(utils::packageVersion("renv")) != renv_version) {
  source_url <- sprintf("https://cloud.r-project.org/src/contrib/Archive/renv/renv_%s.tar.gz", renv_version)
  # Archive versions are immutable. renv has no compiled code.
  install.packages(source_url, repos = NULL, type = "source", lib = bootstrap)
}
renv::restore(project = project, lockfile = file.path(project, "renv.lock"), prompt = FALSE)
message("Environment restored. Use scripts/run_analysis.R to check, prepare or run.")
