# Synthetic examples only: no participant records, packages or model fitting.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
project <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."), winslash = "/", mustWork = TRUE)
e <- new.env(parent = globalenv())
sys.source(file.path(project, "R/answers_analysis.R"), envir = e)
stopifnot(!exists("codebook_results", envir = e, inherits = FALSE))
stopifnot(e$ANALYSIS_SEED == 20260817L, e$N_EDGE_BOOT == 1000L, e$N_CASE_BOOT == 1000L, e$N_CUSTOM_BOOT == 1000L, e$N_P1_BOOT == 1000L, e$N_CORES == 1L)
# Construct small synthetic item-level examples with independently known totals.
x <- as.data.frame(setNames(rep(list(rep(0, 4)), 20), paste0("cesd_", 1:20)))
for (prefix in c("isi", "gad")) for (j in 1:7) x[[paste0(prefix, "_", j)]] <- 0
for (j in 1:5) { x[[paste0("ddnsi_", j)]] <- 0; x[[paste0("cssrs_3mo_", j)]] <- 0 }
for (j in 1:15) x[[paste0("inq_", j)]] <- 1
for (j in c(6,8,13,15)) x[[paste0("upps_", j)]] <- 1
x$psqi_9 <- 1; x$socialjetlag <- c(-90,60,0,0)
x$cssrs_3mo_2[2] <- 1; x$cssrs_3mo_5[2] <- 1
x$cssrs_3mo_1[3] <- NA
x$cesd_11 <- NA # Excluded sleep item must not propagate NA.
x$cesd_1[4] <- NA; x$isi_1[4] <- NA
s <- e$score_registered_measures(x)
stopifnot(identical(s$cesd_no_sleep[1:3], rep(12,3)), is.na(s$cesd_no_sleep[4]), is.na(s$isi_total[4]), s$inq_perceivedburden[1] == 6, s$inq_thwartedbelong[1] == 45, s$upps_negativeurgency[1] == 16, s$psqi_c1[1] == 0)
stopifnot(s$socialjetlag_h[1] == -1.5, s$si_3mo_count[2] == 2, s$si_3mo_maxseverity[2] == 5, s$si_3mo_yn[2] == 1, is.na(s$si_3mo_maxseverity[3]), s$si_3mo_maxseverity[1] == 0)
stopifnot(e$same_with_na(c(1,NA),c(1,NA)), !e$same_with_na(c(1,NA),c(1,0)))
# Decimal age boundary: 25.9 is in primary; 26 is sensitivity only.
names <- c("id","visit","age","isi_1","isi_2","isi_3","isi_7","ddnsi_total","socialjetlag_h","cesd_no_sleep","gad_total","inq_perceivedburden","inq_thwartedbelong","upps_negativeurgency","si_3mo_maxseverity")
p <- as.data.frame(setNames(rep(list(rep(0,5)),length(names)),names))
p$id <- 1:5; p$age <- c(17.9,18,25.9,26,NA)
c <- e$define_analysis_cohorts(p)
stopifnot(identical(c$primary_18_25$id,c(2L,3L)), nrow(c$age_sensitivity)==4L)
# Raw release guard rejects incompatible CSV before scoring.
f <- tempfile(fileext=".csv"); write.csv(data.frame(id=1,visit=1),f,row.names=FALSE)
stopifnot(inherits(try(e$read_and_check_raw_data(f),silent=TRUE),"try-error")); unlink(f)
specs <- e$define_registered_networks()
stopifnot(length(specs)==6L)
for (path in c("R/answers_analysis.R","scripts/run_analysis.R","scripts/setup.R")) invisible(parse(file.path(project,path)))
cat("PASS: source isolation, frozen settings, scoring/missingness, signed jet lag, severity/count distinction, decimal-age cohorts, input rejection, six network specifications and script syntax.\n")
