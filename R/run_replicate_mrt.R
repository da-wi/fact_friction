run_replicate_mrt <- function(N, availability, true_effect, seed, scale,
                              g0 = -0.10, priors = TRUE, ...) {
  sim <- simulate_mrt(n_subjects = N, scale = scale, availability = availability,
                      g1_mean = true_effect, g0_mean = g0, seed = seed, ...)
  res <- fit_once_mrt(sim, scale = scale, priors = priors)
  res |> mutate(
    N = N, availability = availability, true_effect = true_effect, true_g0 = g0,
    seed = seed,
    n_obs     = sum(sim$TD1 == 0),          # answered + available assessments
    n_prompts = sum(sim$TD1 == 1),
    detected  = ok & (l95 > 0 | u95 < 0),    # 95% CI of g1 excludes 0
    detected0 = ok & (l95_0 > 0 | u95_0 < 0)
  )
}
