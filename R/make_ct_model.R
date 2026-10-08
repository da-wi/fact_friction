make_ct_model <- function() {
  ctModel(
    type = "ct",
    n.latent = 2,
    n.manifest = 4,
    n.TDpred = 1,
    LAMBDA = matrix(
      c(1, 0,
        #"loading_11", 0,
        0.8, 0,
        0, 1,
        0, 1.1),
      nrow = 4, byrow = TRUE
    ),
    TDPREDEFFECT = matrix(
      c("gamma_r0 + gamma_r1 * eta1",  # η_r equation gets TD1
        0),                            # η_c no TD1 effect
      nrow = 2, byrow = TRUE
    ),
    PARS = c("gamma_r0", "gamma_r1"),
    silent = TRUE
  )
}

make_ct_model2 <- function() {
  
  model_nonlinear <- ctModel(
    type = "ct",
    n.latent = 3,
    n.manifest = 4, # NEU: 4 Indikatoren
    n.TDpred = 1,
    latentNames = c("Load", "Agency", "ResetState"),
    manifestNames = c("load1", "load2", "agency1", "agency2"),
    TDpredNames = "mrt_trigger",
    
    DRIFT = matrix(c(
      "d_LL", "d_LA", "eff_L_base + eff_L_inter * Load", 
      "d_AL", "d_AA", "eff_A_base + eff_A_inter * Load",
      0,      0,      "decay_Reset"
    ), 3, 3, byrow=TRUE),
    
    TDPREDEFFECT = matrix(c(0, 0, 1), 3, 1),
    
    PARS = matrix(c("eff_L_base", "eff_L_inter", "eff_A_base", "eff_A_inter"), ncol = 1),
    
    # NEU: Das Measurement Model (LAMBDA)
    # Spalte 1 = Load, Spalte 2 = Agency, Spalte 3 = ResetState
    LAMBDA = matrix(c(
      1,           0,           0,   # load1 misst Load
      "lambda_L2", 0,           0,   # load2 misst Load (geschätzt)
      0,           1,           0,   # agency1 misst Agency
      0,           "lambda_A2", 0    # agency2 misst Agency (geschätzt)
    ), 4, 3, byrow = TRUE),
    
    CINT = matrix(c("cint_L||FALSE", "cint_A||FALSE", 0), 3, 1),
    
    DIFFUSION = matrix(c(
      "diff_L", 0, 0,
      "diff_LA", "diff_A", 0,
      0, 0, 0
    ), 3, 3, byrow = TRUE),
    
    # MANIFESTMEANS auf 0 setzen (wie in der Simulation)
    MANIFESTMEANS = matrix(0, nrow = 4, ncol = 1),
    
    # NEU: Individuelle Fehlervarianzen für die 4 Indikatoren
    MANIFESTVAR = matrix(c(
      "err_L1", 0, 0, 0,
      0, "err_L2", 0, 0,
      0, 0, "err_A1", 0,
      0, 0, 0, "err_A2"
    ), 4, 4, byrow = TRUE),
    
    T0MEANS = matrix(c("t0_L", "t0_A", 0), 3, 1),
    T0VAR = diag(0.1, 3)
  )
  
  # NA-Zeilen aus der Parameter-Tabelle putzen (Wichtig für Stan!)
  model_nonlinear$pars <- model_nonlinear$pars[!is.na(model_nonlinear$pars$matrix), ]
  
  return(model_nonlinear)
}
