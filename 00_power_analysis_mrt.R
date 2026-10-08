# ============================================================================ #
# Power simulation – Study 2 (H4a): conflict-dependent distal prompt effect (MRT)
# Same structure as 00_power_analysis.R: model compiled ONCE, then
# ctStanFit(fit = FALSE) + stanoptimis(sm = sm) + rebuild_ctfit_from_manual()
# per replicate. Grid: N x availability x true_effect x rep.
# ============================================================================ #
library(ctsem)
library(tidyverse)
library(RhpcBLASctl)
library(tictoc)
library(yaml)
library(future.apply)
library(progressr)
handlers(global = TRUE)
handlers("txtprogressbar")

rstan::rstan_options(auto_write = TRUE)
blas_set_num_threads(1)
omp_set_num_threads(1)
ctsem_cores <- 8
plan(multisession, workers = 1)

source("R/simulate_mrt.R")
source("R/make_ct_model_mrt.R")
source("R/compile_model_mrt.R")
source("R/fit_once_mrt.R")
source("R/run_replicate_mrt.R")
source("R/rebuild_ctfit_from_manual.R")

# 1) Calibration (cached) -> muC is part of the model code, so the compiled
#    file name carries it. Delete model/scale_mrt.rds to recalibrate.
scale_mrt <- get_scale_mrt("model/scale_mrt.rds")
print(scale_mrt)

ctmodel <- make_ct_model_mrt(muC = scale_mrt[["muC"]])
path_compiled_model <- sprintf("model/ctmodel_mrt_muC%.4f.rds", scale_mrt[["muC"]])
sm <- compile_model_mrt(path_compiled_model, model = ctmodel, scale = scale_mrt)
###############################################################################

# default parameters (overwritten by config/power_analysis_mrt.yaml)
params <- list(Ns = 300, availability = 0.5, true_effect = -0.20, g0 = -0.10,
               n_reps = 2, priors = TRUE, days = 28, compliance = 0.80)
cfg <- "config/power_analysis_mrt.yaml"
if (file.exists(cfg)) {
  message("Using ", cfg)
  params <- modifyList(params, yaml::read_yaml(cfg))
} else message("No ", cfg, " found. Using default parameters.")
message(glue::glue(
  "Simulating with parameters:\n",
  "  Ns           = {toString(params$Ns)}\n",
  "  availability = {toString(params$availability)}\n",
  "  true_effect  = {toString(params$true_effect)}  (g0 = {params$g0})\n",
  "  n_reps       = {params$n_reps}, priors = {params$priors}"
))

set.seed(123)
tic()
pg <- tidyr::expand_grid(
  N = params$Ns,
  availability = params$availability,
  true_effect = params$true_effect,
  rep = seq_len(params$n_reps)
) |>
  dplyr::mutate(seed = sample.int(1e9, size = dplyr::n(), replace = FALSE))

progressr::with_progress({
  p <- progressr::progressor(along = seq_len(nrow(pg)))
  pg$out <- future.apply::future_lapply(
    X = seq_len(nrow(pg)),
    FUN = function(i) {
      on.exit(p(message = paste("Finished replicate", i)))
      options(mc.cores = ctsem_cores)
      run_replicate_mrt(
        N = pg$N[i],
        availability = pg$availability[i],
        true_effect = pg$true_effect[i],
        seed = pg$seed[i],
        scale = scale_mrt,
        g0 = params$g0,
        priors = isTRUE(params$priors),
        days = params$days,
        compliance = params$compliance
      )
    },
    future.stdout = FALSE,
    future.seed = TRUE
  )
})

runtime <- toc(quiet = TRUE)
runtime_sec <- runtime$toc - runtime$tic
pg <- structure(
  pg,
  params = params,
  scale = scale_mrt,
  runtime_sec = runtime_sec,
  runtime_readable = sprintf("%.2f minutes", runtime_sec / 60),
  timestamp = Sys.time()
)

dir.create("results/power_mrt", recursive = TRUE, showWarnings = FALSE)
saveRDS(pg, file = sprintf("results/power_mrt/power_results_raw.%s.rds",
                           format(attr(pg, "timestamp"), "%Y%m%d_%H%M%S")))
