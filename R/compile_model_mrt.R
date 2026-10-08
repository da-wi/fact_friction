# Same logic as compile_model(), but model and example data are passed in
# (compile_model() uses the global ctmodel and generateSimData()).
compile_model_mrt <- function(path_compiled_model, model, scale) {
  if (!file.exists(path_compiled_model)) {
    dir.create(dirname(path_compiled_model), recursive = TRUE, showWarnings = FALSE)
    ctfit <- ctStanFit(
      datalong = simulate_mrt(n_subjects = 10, scale = scale, days = 7, seed = 1),
      ctstanmodel = model,
      cores = 4,
      optimize = TRUE,
      fit = FALSE,
      saveCompile = TRUE,
      verbose = 0,
      priors = FALSE)
    sm <- rstan::stan_model(model_code = ctfit$stanmodeltext)
    saveRDS(sm, file = path_compiled_model)
  } else {
    sm <- readRDS(file = path_compiled_model)
  }
  if (!inherits(sm, "stanmodel")) stop(sprintf("%s is not a valid stan model.", path_compiled_model))
  message(sprintf("Model found at %s", path_compiled_model))
  invisible(sm)
}
