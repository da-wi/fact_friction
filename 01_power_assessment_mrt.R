library(tidyverse)

(result_file_path <- rev(dir("results/power_mrt/"))[1])
pg <- readRDS(paste0("results/power_mrt/", result_file_path))
attr(pg, "runtime_readable")

res <- pg |> select(-any_of(c("N", "availability", "true_effect", "seed"))) |> unnest(out)

(summary_table <- res |>
  group_by(N, availability, true_effect) |>
  summarise(
    n_reps     = n(),
    n_ok       = sum(ok),
    prompts_pp = mean(n_prompts) / first(N),
    obs_pp     = mean(n_obs) / first(N),
    mean_est   = mean(est[ok]),
    rel_bias_pct = (mean_est - first(true_effect)) / abs(first(true_effect)) * 100,
    sd_est     = sd(est[ok]),
    power      = mean(detected[ok]),
    power_lo   = pmax(0, power - 1.96 * sqrt(power * (1 - power) / n_ok)),  # MC uncertainty
    power_hi   = pmin(1, power + 1.96 * sqrt(power * (1 - power) / n_ok)),
    power_g0   = mean(detected0[ok]),
    .groups = "drop"
  ))

ggplot(summary_table, aes(x = factor(availability), y = power,
                          colour = factor(true_effect), group = factor(true_effect))) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = power_lo, ymax = power_hi), width = .05) +
  geom_hline(yintercept = 0.8, colour = "red") +
  facet_wrap(~ N, labeller = label_both) +
  labs(title = "Power H4a: conflict x prompt (distal, ctsem)",
       x = "Availability (share of prompts with news use)",
       y = "Power (95% CI excludes 0)", colour = "true g1 (SD)") +
  theme_minimal()

# Helper quantities for the proposal text
source("R/simulate_mrt.R")
cat(sprintf("Half-life derogation: %.1f h; g1 = -0.20 left at next daytime assessment: %.3f SD (%.1f VAS points at within-SD 15)\n",
            half_life_h(), effect_next(), effect_vas(effect_next())))
