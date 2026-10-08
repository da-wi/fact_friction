run_replicate <- function(N, true_effect, seed) {
  #sim <- simulate_ct_data(
  #  n_subjects = N,
  #  gamma_r1_mean = true_effect,
  #  seed = seed
  #)
  sim <- generateSimData(n_subjects = N, eff_L_inter =true_effect, seed = seed)
  res <- fit_once(sim)
  res$N <- N
  res$seed <- seed
  res$detected <- with(res, ok & (l95 > 0 | u95 < 0))  # 95% CI excludes 0
  res |> select(ok, est, l95, u95, detected)
}


run_replicate <- function(N, T_pts, true_effect, time_var, miss_type, seed) {
  
  sim_data <- generateSimData(
    n_subjects = N, t_points = T_pts, eff_L_inter = true_effect, 
    time_variance = time_var, missing_type = miss_type, seed = seed
  )
  
  res <- fit_both_models(sim_data)
  
  # Tracking Variablen hinzufügen
  res <- res %>% mutate(
    N = N,
    n_obs = n_obs,
    T_points = T_pts,
    true_effect = true_effect,
    time_variance = time_var,
    missing_type = miss_type,
    seed = seed,
    detected = ok & (l95 > 0 | u95 < 0) # Ist 0 außerhalb des 95% CI?
  )
  
  return(res)
}
