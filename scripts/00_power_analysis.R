# ============================================================================ #
# Power simulation for continuous-time model (ctsem)
# Author: David Willinger
# Date: 2025-10-30
# Description:
#   Runs Monte Carlo simulations across sample sizes (Ns) and replicates
#   to estimate the power to detect a true effect (gamma_r1 = true_effect).
#   Results are saved as timestamped .rds files.
# ============================================================================ #
library(ctsem)
library(tidyverse)
library(RhpcBLASctl)
library(tictoc)
library(yaml)
library(future.apply)
library(progressr)
library(lme4)
handlers(global = TRUE)
handlers("txtprogressbar") 

rstan::rstan_options(auto_write = TRUE)
blas_set_num_threads(1)
omp_set_num_threads(1)
ctsem_cores <- 8
plan(multisession, workers = 1)

#plan(sequential)

source("R/simulate_ct_data.R")
source("R/make_ct_model.R")
source("R/fit_once.R")
source("R/run_replicate.R")
source("R/rebuild_ctfit_from_manual.R")
source("R/compile_model.R")

ctmodel <- make_ct_model2()
path_compiled_model <- "model/ctmodel_compiled3.rds"

sm <- compile_model(path_compiled_model)
###############################################################################

# default parameters
Ns <- c(10, 20)
n_reps <- 1
true_effect <- c(-0.2)

params <- tryCatch(
  {
    message("Using config/power_analysis.yaml:")
    params <- yaml::read_yaml("config/power_analysis.yaml")
    Ns <- params$Ns; n_reps <- params$n_reps; true_effect <- params$true_effect
  },
  warning = function(w) {
    message("No power_analysis.yaml found in `config`. Using default parameters.")
    list(Ns = Ns, n_reps = n_reps, true_effect = true_effect)
  },
  finally = {
    message(glue::glue(
      "Simulating with parameters:\n",
      "  Ns          = {toString(Ns)}\n",
      "  n_reps         = {n_reps}\n",
      "  true_effect = {toString(true_effect)}"
    ))
  }
)

set.seed(123)
tic()
pg <- tidyr::expand_grid(
  N = Ns,
  rep = seq_len(n_reps),
  true_effect = true_effect
) |>
  dplyr::mutate(seed = sample.int(1e9, size = dplyr::n(), replace = FALSE))

progressr::with_progress({
  p <- progressr::progressor(along = seq_len(nrow(pg)))
  pg$out <- future.apply::future_lapply(
    X = seq_len(nrow(pg)),
    FUN = function(i) {
      on.exit(p(message = paste("Finished replicate", i)))
      options(mc.cores = ctsem_cores)
      run_replicate(
        N = pg$N[i],
        true_effect = pg$true_effect[i],
        seed = pg$seed[i]
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
  runtime_sec = runtime_sec,
  runtime_readable = sprintf("%.2f minutes", runtime_sec / 60),
  timestamp = Sys.time()
)

dir.create("results/power", recursive = TRUE, showWarnings = FALSE)


saveRDS(pg, file = sprintf("results/power/power_results_raw.%s.rds", 
                       format(attr(pg, "timestamp"), "%Y%m%d_%H%M%S")))

########################################################
#########################################################

set.seed(125)
tic()
pg <- tidyr::expand_grid(
  N = c(100),                      # Power Check
  T_pts = c(80),                   # Länge der Studie
  true_effect = c(-0.2),          # False Positive vs. Echter Effekt
  time_var = c("Low", "High"),         # LMM-Killer 1
  miss_type = c("None", "MCAR","MNAR"),       # LMM-Killer 2
  rep = seq_len(15)                     # Z.B. 5 zum Testen, später 100+
) %>%
  dplyr::mutate(seed = sample.int(1e9, size = dplyr::n(), replace = FALSE))

progressr::with_progress({
  p <- progressr::progressor(along = seq_len(nrow(pg)))
  pg$out <- future.apply::future_lapply(
    X = seq_len(nrow(pg)),
    FUN = function(i) {
      on.exit(p(message = paste("Finished replicate", i)))
      options(mc.cores = ctsem_cores)
      run_replicate(
        N = pg$N[i],
        true_effect = pg$true_effect[i],
        seed = pg$seed[i],
        time_var = pg$time_var[i],
        T_pts = pg$T_pts[i],
        miss_type = pg$miss_type[i]
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
  runtime_sec = runtime_sec,
  runtime_readable = sprintf("%.2f minutes", runtime_sec / 60),
  timestamp = Sys.time()
)

results_flat <- pg |>
  select(out) |>
  unnest(out, names_repair = "minimal") 

(summary_table <- results_flat |>
  filter(ok == TRUE)|>
  group_by(model, time_variance, missing_type, N, T_points, true_effect) |>
  summarise(
    n_successful_reps = n(),
    mean_est = mean(est, na.rm = TRUE),
    
    # Bias: Schätzung minus Wahrheit
    abs_bias = mean(est, na.rm = TRUE) - first(true_effect),
    
    # Relativer Bias in Prozent (sehr wichtig für Paper!)
    rel_bias_pct = (mean(est, na.rm = TRUE) - first(true_effect)) / abs(first(true_effect)) * 100,
    
    # Power: Wie oft wurde der Effekt erkannt (0 nicht im CI)?
    power = mean(detected, na.rm = TRUE),
    
    # Precision: Standardabweichung der Schätzungen über die Replikationen
    sd_est = sd(est, na.rm = TRUE),
    
    .groups = "drop"
  ))

ggplot(summary_table %>% filter(true_effect != 0), 
       aes(x = factor(N), y = power, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(time_variance ~ missing_type) +
  geom_hline(yintercept = 0.8, linetype = "solid", color = "red") +
  labs(title = "Power-Vergleich: ctsem vs. nlme",
       subtitle = "Rote Linie = 80% Power-Schwelle",
       x = "Stichprobengröße (N)",
       y = "Power (Anteil signifikanter Ergebnisse)") +
  theme_minimal()


ggplot(summary_table %>% filter(true_effect != 0), 
       aes(x = factor(N), y = power, fill = model)) +
  geom_line(aes(group=model, color=model)) +
  facet_grid(time_variance ~ missing_type) +
  geom_hline(yintercept = 0.8, linetype = "solid", color = "red") +
  labs(title = "Power-Vergleich: ctsem vs. nlme",
       subtitle = "Rote Linie = 80% Power-Schwelle",
       x = "Stichprobengröße (N)",
       y = "Power (Anteil signifikanter Ergebnisse)") +
  theme_minimal()


