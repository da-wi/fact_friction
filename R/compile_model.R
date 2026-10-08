compile_model <- function(path_compiled_model = NULL) {
    if(is.null(path_compiled_model)) stop("Set path for model.")
  
    sm <- NULL
    if (!file.exists(path_compiled_model)) {
    # prepare model
    ctfit <- ctStanFit(
    #  datalong = simulate_ct_data(n_subjects = 10) |> select(id, time, Y1, Y2, Y3, Y4, TD1),
      datalong = generateSimData(n_subjects = 10),
      ctstanmodel = ctmodel,
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
  
  if (isClass(sm,"stanmodel") ) {
    message(sprintf("Model found at %s", path_compiled_model))
    invisible(sm)
  } else{
    stop("%s is not a valid stan model.", path_compiled_model)
  }
}
