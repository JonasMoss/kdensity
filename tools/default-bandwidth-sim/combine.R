#!/usr/bin/env Rscript
# Combine results/*.csv into summary tables.
#
#   Rscript combine.R [results_dir]
#
# Writes summary.csv (one row per kernel/start/target/n/selector) and prints
# the MISE ratio to the oracle, the median per-replication ISE ratio, the share
# of replications where a selector failed, and the mean selection time.
args <- commandArgs(trailingOnly = TRUE)
here <- tryCatch(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(), value = TRUE)))), error = function(e) ".")
dir <- if (length(args)) args[[1]] else file.path(here, "results")
files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
raw <- do.call(rbind, lapply(files, utils::read.csv))
key <- c("kernel", "start", "target", "n", "rep")
oracle <- raw[raw$selector == "oracle", c(key, "ise")]
names(oracle)[names(oracle) == "ise"] <- "oracle_ise"
raw <- merge(raw, oracle, by = key)
raw$ratio <- raw$ise / raw$oracle_ise
groups <- c("kernel", "start", "target", "n", "selector")
summary <- do.call(rbind, lapply(
  split(raw, raw[groups], drop = TRUE),
  function(d) {
    ok <- is.finite(d$ise)
    data.frame(
      d[1, groups],
      reps = nrow(d),
      failed = mean(!ok),
      mise = mean(d$ise[ok]),
      mise_ratio = mean(d$ise[ok]) / mean(d$oracle_ise[ok]),
      median_ratio = stats::median(d$ratio[ok]),
      median_h = stats::median(d$h, na.rm = TRUE),
      seconds = mean(d$seconds, na.rm = TRUE),
      hs_fallback = mean(d$hs_fallback),
      row.names = NULL
    )
  }
))
summary <- summary[do.call(order, summary[groups]), ]
utils::write.csv(summary, file.path(dirname(dir), "summary.csv"), row.names = FALSE)
cat("cells:", length(files), " rows:", nrow(raw), "\n\n")
wide <- stats::reshape(
  summary[summary$selector != "oracle",
          c("kernel", "start", "target", "n", "selector", "mise_ratio")],
  idvar = c("kernel", "start", "target", "n"),
  timevar = "selector", direction = "wide"
)
names(wide) <- sub("^mise_ratio\\.", "", names(wide))
print(format(wide, digits = 3), row.names = FALSE)
