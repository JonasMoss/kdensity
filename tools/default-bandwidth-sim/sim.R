#!/usr/bin/env Rscript
# Default-bandwidth simulation for the asymmetric kernels (board TASK-4).
#
# Question: can the O(n^2) UCV defaults of the asymmetric kernels be replaced
# by the package's existing closed-form rules (GR for the gamma kernels, HS for
# the beta kernels, JH for the Gaussian-copula kernel), applied regardless of
# the parametric start, without a material loss in ISE? The biased kernels
# reuse the corrected kernels' rules (decision of 2026-10-28). GR needs a
# fallback when the fitted gamma shape is at most 3/2; capping the shape at 2
# or 3 (keeping the fitted mean) are the candidates.
#
# Usage:
#   Rscript sim.R list                 print every cell id, one per line
#   Rscript sim.R run <cell> <outdir>  run one cell, write <outdir>/<cell>.csv
#
# A cell is kernel x start x target x n x block. The CSV is written atomically
# (temporary file + rename) and is the cell's completion marker. Data for a
# given (target, n, rep) are identical across kernels and starts, so selector
# comparisons are paired.

library_path <- Sys.getenv("KDENSITY_LIB")
if (nzchar(library_path)) {
  .libPaths(c(library_path, .libPaths()))
}
suppressPackageStartupMessages(library(kdensity))
ns <- asNamespace("kdensity")

# ---------------------------------------------------------------------------
# Design.
# ---------------------------------------------------------------------------

targets <- list(
  gamma2 = list(
    domain = "pos",
    r = function(n) rgamma(n, 2),
    d = function(x) dgamma(x, 2)
  ),
  gamma8 = list(
    domain = "pos",
    r = function(n) rgamma(n, 8),
    d = function(x) dgamma(x, 8)
  ),
  exp = list(
    domain = "pos",
    r = function(n) rexp(n),
    d = function(x) dexp(x)
  ),
  lnorm05 = list(
    domain = "pos",
    r = function(n) rlnorm(n, 0, 0.5),
    d = function(x) dlnorm(x, 0, 0.5)
  ),
  lnorm1 = list(
    domain = "pos",
    r = function(n) rlnorm(n, 0, 1),
    d = function(x) dlnorm(x, 0, 1)
  ),
  weibull15 = list(
    domain = "pos",
    r = function(n) rweibull(n, 1.5),
    d = function(x) dweibull(x, 1.5)
  ),
  gammamix = list(
    domain = "pos",
    r = function(n) {
      ifelse(runif(n) < 0.5, rgamma(n, 2), rgamma(n, 10))
    },
    d = function(x) 0.5 * dgamma(x, 2) + 0.5 * dgamma(x, 10)
  ),
  beta25 = list(
    domain = "unit",
    r = function(n) rbeta(n, 2, 5),
    d = function(x) dbeta(x, 2, 5)
  ),
  beta0808 = list(
    domain = "unit",
    r = function(n) rbeta(n, 0.8, 0.8),
    d = function(x) dbeta(x, 0.8, 0.8)
  ),
  beta52 = list(
    domain = "unit",
    r = function(n) rbeta(n, 5, 2),
    d = function(x) dbeta(x, 5, 2)
  ),
  betamix = list(
    domain = "unit",
    r = function(n) {
      ifelse(runif(n) < 0.5, rbeta(n, 2, 8), rbeta(n, 8, 2))
    },
    d = function(x) 0.5 * dbeta(x, 2, 8) + 0.5 * dbeta(x, 8, 2)
  )
)

combos <- list(
  list(kernel = "gamma", start = "uniform", domain = "pos"),
  list(kernel = "gamma", start = "gamma", domain = "pos"),
  list(kernel = "gamma_biased", start = "uniform", domain = "pos"),
  list(kernel = "beta", start = "uniform", domain = "unit"),
  list(kernel = "beta", start = "beta", domain = "unit"),
  list(kernel = "beta_biased", start = "uniform", domain = "unit"),
  list(kernel = "gcopula", start = "uniform", domain = "unit"),
  list(kernel = "gcopula", start = "beta", domain = "unit")
)

# Replications per n, split into blocks of 50.
replications <- c("100" = 200L, "500" = 200L, "2000" = 100L)
block_size <- 50L

cell_ids <- function() {
  ids <- character()
  for (combo in combos) {
    for (target_name in names(targets)) {
      if (targets[[target_name]]$domain != combo$domain) next
      for (n in as.integer(names(replications))) {
        blocks <- replications[[as.character(n)]] %/% block_size
        for (block in seq_len(blocks)) {
          ids <- c(ids, sprintf(
            "%s__%s__%s__n%d__b%d",
            combo$kernel, combo$start, target_name, n, block
          ))
        }
      }
    }
  }
  ids
}

parse_cell <- function(id) {
  parts <- strsplit(id, "__", fixed = TRUE)[[1]]
  list(
    kernel = parts[[1]],
    start = parts[[2]],
    target = parts[[3]],
    n = as.integer(sub("^n", "", parts[[4]])),
    block = as.integer(sub("^b", "", parts[[5]]))
  )
}

# Stable seed for (target, n): the same data for every kernel and start.
data_seed <- function(target, n, rep) {
  codes <- utf8ToInt(paste0(target, "_", n))
  base <- sum(codes * seq_along(codes) * 7919) %% 100000L
  as.integer(base * 1000L + rep)
}

# ---------------------------------------------------------------------------
# Selectors.
# ---------------------------------------------------------------------------

# Gamma-reference bandwidth with the fitted shape capped from below; the rate
# keeps the fitted mean. cap = 0 is the package's GR when the shape exceeds
# 3/2 and NA otherwise.
gr_capped <- function(x, cap) {
  fit <- univariateML::mlgamma(x)
  shape <- max(unname(fit[["shape"]]), cap)
  if (shape <= 1.5) {
    return(NA_real_)
  }
  rate <- shape / mean(x)
  n <- length(x)
  log_v <- 0.5 * log(rate) + lgamma(shape - 0.5) - lgamma(shape)
  log_r <- 3 * log(rate) + lgamma(2 * shape - 3) -
    (2 * shape - 3) * log(2) - 2 * lgamma(shape) +
    log(shape - 1) + log(3 * shape - 4) - log(4)
  exp((2 / 5) * (log_v - log(2 * sqrt(pi)) - log_r - log(n)))
}

upper_cap <- function(kernel) {
  if (kernel %in% c("beta", "beta_biased")) {
    return(0.25 - 1e-6)
  }
  if (kernel == "gcopula") {
    return(0.5)
  }
  Inf
}

closed_form_rules <- function(x, kernel) {
  if (kernel %in% c("gamma", "gamma_biased")) {
    return(list(
      GR = gr_capped(x, 0),
      GRc2 = gr_capped(x, 2),
      GRc3 = gr_capped(x, 3)
    ))
  }
  if (kernel %in% c("beta", "beta_biased")) {
    hs_fallback <- FALSE
    hs <- withCallingHandlers(
      ns$compute_hs_bandwidth(x),
      warning = function(w) {
        hs_fallback <<- TRUE
        invokeRestart("muffleWarning")
      }
    )
    return(list(HS = hs, HS_fallback_used = hs_fallback))
  }
  list(JH = ns$bw_environment$JH(x))
}

# Centre of the wide search interval: the rule that would become the default.
rule_centre <- function(rules, kernel) {
  if (kernel %in% c("gamma", "gamma_biased")) {
    return(rules$GRc2)
  }
  if (kernel %in% c("beta", "beta_biased")) {
    return(rules$HS)
  }
  rules$JH
}

search_interval <- function(centre, kernel) {
  c(centre / 10, min(centre * 10, upper_cap(kernel)))
}

ise <- function(x, h, kernel, start, truth, support) {
  tryCatch({
    fit <- suppressWarnings(
      kdensity(x, kernel = kernel, start = start, bw = h)
    )
    integrate(
      function(t) (fit(t) - truth(t))^2,
      lower = support[[1]],
      upper = support[[2]],
      subdivisions = 2000L,
      rel.tol = 1e-6
    )$value
  }, error = function(e) NA_real_)
}

ucv_wide <- function(x, kernel, start, support, interval) {
  objective <- ns$make_ucv_objective(x, kernel, start, support)
  safe <- function(log_h) {
    value <- tryCatch(objective(exp(log_h)), error = function(e) NA_real_)
    if (is.finite(value)) value else .Machine$double.xmax
  }
  out <- stats::optimize(safe, log(interval), tol = 0.01)
  exp(out$minimum)
}

oracle <- function(x, kernel, start, truth, support, interval) {
  grid <- exp(seq(log(interval[[1]]), log(interval[[2]]), length.out = 25))
  values <- vapply(grid, function(h) {
    ise(x, h, kernel, start, truth, support)
  }, numeric(1))
  if (all(!is.finite(values))) {
    return(c(h = NA_real_, ise = NA_real_))
  }
  best <- which.min(values)
  lower <- grid[max(best - 1L, 1L)]
  upper <- grid[min(best + 1L, length(grid))]
  if (lower < upper) {
    refined <- stats::optimize(function(log_h) {
      value <- ise(x, exp(log_h), kernel, start, truth, support)
      if (is.finite(value)) value else .Machine$double.xmax
    }, log(c(lower, upper)), tol = 0.005)
    if (is.finite(refined$objective) && refined$objective < values[[best]]) {
      return(c(h = exp(refined$minimum), ise = refined$objective))
    }
  }
  c(h = grid[[best]], ise = values[[best]])
}

timed <- function(expr) {
  start <- proc.time()[["elapsed"]]
  value <- tryCatch(expr, error = function(e) NA_real_)
  list(value = value, seconds = proc.time()[["elapsed"]] - start)
}

# ---------------------------------------------------------------------------
# One replication.
# ---------------------------------------------------------------------------

run_rep <- function(cell, rep) {
  target <- targets[[cell$target]]
  support <- if (target$domain == "pos") c(0, Inf) else c(0, 1)
  set.seed(data_seed(cell$target, cell$n, rep))
  x <- target$r(cell$n)

  rules_time <- timed(closed_form_rules(x, cell$kernel))
  rules <- rules_time$value
  if (!is.list(rules)) {
    rules <- list(rule_failed = NA_real_)
  }
  centre <- rule_centre(rules, cell$kernel)
  if (length(centre) != 1L || !is.finite(centre)) {
    centre <- stats::bw.nrd0(x)
  }
  interval <- search_interval(centre, cell$kernel)

  bandwidths <- rules[setdiff(names(rules), "HS_fallback_used")]
  seconds <- as.list(rep(rules_time$seconds, length(bandwidths)))
  names(seconds) <- names(bandwidths)

  package_ucv <- timed(suppressWarnings(
    ns$bw_environment$ucv(x, cell$kernel, cell$start, support)
  ))
  bandwidths$ucv_pkg <- package_ucv$value
  seconds$ucv_pkg <- package_ucv$seconds

  wide <- timed(ucv_wide(x, cell$kernel, cell$start, support, interval))
  bandwidths$ucv_wide <- wide$value
  seconds$ucv_wide <- wide$seconds

  best <- oracle(x, cell$kernel, cell$start, target$d, support, interval)

  rows <- lapply(names(bandwidths), function(selector) {
    h <- bandwidths[[selector]]
    data.frame(
      kernel = cell$kernel,
      start = cell$start,
      target = cell$target,
      n = cell$n,
      rep = rep,
      selector = selector,
      h = h,
      ise = if (is.finite(h)) {
        ise(x, h, cell$kernel, cell$start, target$d, support)
      } else {
        NA_real_
      },
      seconds = seconds[[selector]],
      interval_lower = interval[[1]],
      interval_upper = interval[[2]],
      hs_fallback = isTRUE(rules$HS_fallback_used),
      stringsAsFactors = FALSE
    )
  })
  rows[[length(rows) + 1L]] <- data.frame(
    kernel = cell$kernel,
    start = cell$start,
    target = cell$target,
    n = cell$n,
    rep = rep,
    selector = "oracle",
    h = best[["h"]],
    ise = best[["ise"]],
    seconds = NA_real_,
    interval_lower = interval[[1]],
    interval_upper = interval[[2]],
    hs_fallback = isTRUE(rules$HS_fallback_used),
    stringsAsFactors = FALSE
  )
  do.call(rbind, rows)
}

run_cell <- function(id, outdir) {
  path <- file.path(outdir, paste0(id, ".csv"))
  if (file.exists(path)) {
    message("skip ", id)
    return(invisible(path))
  }
  cell <- parse_cell(id)
  reps <- (cell$block - 1L) * block_size + seq_len(block_size)
  result <- do.call(rbind, lapply(reps, function(rep) run_rep(cell, rep)))
  temporary <- paste0(path, ".tmp")
  utils::write.csv(result, temporary, row.names = FALSE)
  file.rename(temporary, path)
  message("done ", id)
  invisible(path)
}

# ---------------------------------------------------------------------------
# Command line.
# ---------------------------------------------------------------------------

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1L && args[[1]] == "list") {
    writeLines(cell_ids())
  } else if (length(args) == 3L && args[[1]] == "run") {
    dir.create(args[[3]], showWarnings = FALSE, recursive = TRUE)
    run_cell(args[[2]], args[[3]])
  } else {
    stop("usage: sim.R list | sim.R run <cell> <outdir>")
  }
}
