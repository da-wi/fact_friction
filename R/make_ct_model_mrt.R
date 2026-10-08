# Analysis model for the MRT (impulse model M1): conflict-dependent prompt effect
# as state-dependent TDPREDEFFECT on derogation R, centred at mean conflict muC.
#   gm0 = impulse at average conflict, gm1 = conflict x prompt (raw latent units)
# NOTE: muC is written into the Stan code -> compiled model is specific to muC.
make_ct_model_mrt <- function(muC) {
  td <- sprintf("(gm0 + gm1 * (C - %.4f))", muC)
  ctModel(type = "ct",
    n.latent = 2, latentNames = c("C", "R"),
    n.manifest = 4, manifestNames = c("Y1", "Y2", "Y3", "Y4"),
    n.TDpred = 1, TDpredNames = "TD1",
    LAMBDA = matrix(c(0, 1,
                      0, "lR2",
                      1, 0,
                      "lC2", 0), 4, 2, byrow = TRUE),
    DRIFT = matrix(c("aC",  "bRC",
                     "bCR", "aR"), 2, 2, byrow = TRUE),
    PARS = matrix(c("gm0", "gm1"), ncol = 1),
    TDPREDEFFECT = matrix(c(0, td), 2, 1),
    CINT = matrix(c("cC||TRUE", "cR||TRUE"), 2, 1),
    DIFFUSION = matrix(c("dC", 0,
                         "dCR", "dR"), 2, 2, byrow = TRUE),
    MANIFESTMEANS = matrix(0, 4, 1),
    MANIFESTVAR = matrix(c("e1", 0, 0, 0,
                           0, "e2", 0, 0,
                           0, 0, "e3", 0,
                           0, 0, 0, "e4"), 4, 4, byrow = TRUE),
    T0MEANS = matrix(c("t0C", "t0R"), 2, 1),
    T0VAR = matrix(c("t0vC", 0,
                     "t0vCR", "t0vR"), 2, 2, byrow = TRUE),
    silent = TRUE)
}
