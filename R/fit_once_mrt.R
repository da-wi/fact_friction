# Uses the globals `ctmodel`, `sm`, `ctsem_cores` exactly like myCtStanFit()
# (rebuild_ctfit_from_manual() also reads the global `sm`).
myCtStanFit_mrt <- function(sim_data, priors = TRUE) {
  ctfit <- ctStanFit(
    datalong = sim_data,
    ctstanmodel = ctmodel,
    cores = ctsem_cores,
    optimize = TRUE,
    fit = FALSE,
    saveCompile = TRUE,
    verbose = 0,
    priors = priors
  )
  optres <- ctsem:::stanoptimis(
    standata = ctfit$standata,
    sm = sm,
    init = NULL,
    tol = 1e-7,
    matsetup = data.frame(ctfit$ctstanmodel$modelmats$matsetup),
    verbose = 0,
    stochastic = FALSE,
    priors = priors,
    cores = ctfit$args$cores
  )
  rebuild_ctfit_from_manual(ctfit, optres)
}

safe_fit_mrt <- purrr::safely(myCtStanFit_mrt, otherwise = NULL)

#' Fit the MRT model once and return standardized effects
#' g1 (conflict x prompt) = gm1 * SD_C / SD_R ; g0 (prompt at mean conflict) = gm0 / SD_R
fit_once_mrt <- function(sim_data, scale, priors = TRUE) {
  fit <- safe_fit_mrt(sim_data, priors = priors)
  na_row <- tibble(ok = FALSE, est = NA_real_, l95 = NA_real_, u95 = NA_real_,
                   est0 = NA_real_, l95_0 = NA_real_, u95_0 = NA_real_)
  if (!is.null(fit$error)) { print(fit$error); return(na_row) }

  summ <- summary(fit$result)$popmeans |> as_tibble(rownames = "parameter")
  f1 <- scale[["C"]] / scale[["R"]]
  f0 <- 1 / scale[["R"]]
  g1 <- summ |> filter(parameter == "gm1")
  g0 <- summ |> filter(parameter == "gm0")
  if (nrow(g1) != 1 || nrow(g0) != 1) return(na_row)

  tibble(
    ok    = TRUE,
    est   = g1$mean * f1, l95   = g1$`2.5%` * f1, u95   = g1$`97.5%` * f1,
    est0  = g0$mean * f0, l95_0 = g0$`2.5%` * f0, u95_0 = g0$`97.5%` * f0
  )
}
