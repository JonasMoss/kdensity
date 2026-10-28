#!/usr/bin/env Rscript
# Follow-up for the Gaussian-copula kernel: JH (a normal-reference rule on the
# probit scale) is 4-29x the oracle MISE on bimodal data. On the probit scale
# the Gaussian-copula kernel is a Gaussian kernel with sd s = sqrt(1 - rho^2),
# rho = 1 - h^2, so any Gaussian bandwidth s chosen for qnorm(x) maps exactly
# to h = sqrt(1 - sqrt(1 - s^2)). This script evaluates stats::bw.SJ and
# stats::bw.nrd0 on the probit scale, regenerating the simulation's data from
# its seeds and comparing with the stored oracle ISE in results/.
#
#   KDENSITY_LIB=... Rscript probit-rules.R [cores]

args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args)) as.integer(args[[1]]) else 4L
here <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE))))
setwd(here)
source("sim.R")

probit_h <- function(s) {
  s <- min(s, 1 - 1e-8)
  sqrt(1 - sqrt(1 - s^2))
}

files <- list.files("results", pattern = "^gcopula__.*\\.csv$", full.names = TRUE)
stored <- do.call(rbind, lapply(files, utils::read.csv))
stored <- stored[stored$selector %in% c("JH", "ucv_pkg", "oracle"), ]
keys <- unique(stored[, c("start", "target", "n", "rep")])

rows <- parallel::mclapply(seq_len(nrow(keys)), function(i) {
  key <- keys[i, ]
  target <- targets[[key$target]]
  set.seed(data_seed(key$target, key$n, key$rep))
  x <- target$r(key$n)
  z <- stats::qnorm(x)
  h <- c(
    SJ_probit = tryCatch(probit_h(stats::bw.SJ(z)), error = function(e) NA_real_),
    nrd0_probit = probit_h(stats::bw.nrd0(z))
  )
  data.frame(
    key[rep(1, length(h)), ],
    selector = names(h),
    h = unname(h),
    ise = vapply(h, function(b) {
      if (is.finite(b)) ise(x, b, "gcopula", key$start, target$d, c(0, 1)) else NA_real_
    }, numeric(1)),
    row.names = NULL
  )
}, mc.cores = cores)
new <- do.call(rbind, rows)

all <- rbind(stored[, c("start", "target", "n", "rep", "selector", "h", "ise")], new)
oracle <- all[all$selector == "oracle", c("start", "target", "n", "rep", "ise")]
names(oracle)[5] <- "oracle"
all <- merge(all[all$selector != "oracle", ], oracle)
ratio <- aggregate(cbind(ise, oracle) ~ start + target + n + selector, all, mean)
ratio$mise_ratio <- ratio$ise / ratio$oracle
wide <- stats::reshape(ratio[, c("start", "target", "n", "selector", "mise_ratio")],
  idvar = c("start", "target", "n"), timevar = "selector", direction = "wide")
names(wide) <- sub("^mise_ratio\\.", "", names(wide))
wide <- wide[do.call(order, wide[c("start", "target", "n")]), ]
utils::write.csv(wide, "probit-rules.csv", row.names = FALSE)
print(format(wide, digits = 3), row.names = FALSE)
cat("\ngeometric mean / max over target x n:\n")
for (start in unique(wide$start)) {
  d <- wide[wide$start == start, c("JH", "ucv_pkg", "SJ_probit", "nrd0_probit")]
  cat(start, "gmean:", format(exp(colMeans(log(d))), digits = 3),
    " max:", format(apply(d, 2, max), digits = 3), "\n")
}
