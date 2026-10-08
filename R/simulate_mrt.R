# ============================================================================ #
# Data-generating process – "Fact or Friction", Study 2 (H4a)
# Conflict-dependent, DISTAL effect of micro-randomized accuracy prompts on
# source derogation (MRT design).
#
#  - 28 days x 4 prompts (3 random times in blocks 08-11, 12-15, 16-19 h,
#    >= 2 h apart; evening prompt answered 20:00-22:30), 80% compliance
#  - availability: ratings + randomization only at prompts with news use
#  - randomization p = .50 at available prompts
#  - prompt = instantaneous impulse on derogation R at t + eps (after the
#    ratings at t): R <- R + (g0 + g1 * C)  -> distal test at t+1, ...
#  - g0, g1 are STANDARDIZED (within-person SD units, calibrated once without
#    prompts via calibrate_sd_mrt()); g1 = -0.20: a prompt at conflict +1 SD
#    lowers latent derogation by 0.20 within-person SD (on top of g0).
#
# Time unit: 1 unit = TIME_UNIT_H hours (3.5 h = mean daytime gap).
# Output columns: id, time, Y1, Y2 (derogation), Y3, Y4 (conflict), TD1 (prompt)
# ============================================================================ #

TIME_UNIT_H <- 3.5

make_schedule <- function(days = 28) {
  t <- numeric(0)
  for (d in 0:(days - 1)) {
    repeat {                                   # three daytime prompts >= 2 h apart
      x <- c(runif(1, 8, 11), runif(1, 12, 15), runif(1, 16, 19))
      if (all(diff(x) >= 2)) break
    }
    ev <- runif(1, 20, 22.5)                   # evening prompt answered 20:00-22:30
    t <- c(t, d * 24 + c(x, ev))
  }
  t / TIME_UNIT_H
}

#' @param scale named vector c(C = SD_C, R = SD_R, muC = mean C) from
#'   calibrate_sd_mrt(); passed explicitly (no global) so future workers get it.
simulate_mrt <- function(
    n_subjects   = 300,
    scale,
    days         = 28,
    compliance   = 0.80,
    availability = 0.50,
    p_prompt     = 0.50,
    g0_mean = -0.10, g0_sd = 0.05,
    g1_mean = -0.20, g1_sd = 0.05,
    demand  = 0,             # optional proximal demand artefact on R items at t (not modelled)
    a_r = -0.45, a_c = -0.35, beta_c_to_r = 0.15, beta_r_to_c = -0.55,
    b_r_mean = 0.00, b_c_mean = 0.10, b_sd = 0.10,
    sigma_r = 0.20, sigma_c = 0.20, sigma_common = 0.10,
    me_sd = 0.25, steps_per_unit = 50L, eps = 0.01, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  out <- vector("list", n_subjects)
  for (i in seq_len(n_subjects)) {
    g1 <- rnorm(1, g1_mean, g1_sd) * scale[["R"]] / scale[["C"]]
    g0 <- rnorm(1, g0_mean, g0_sd) * scale[["R"]] - g1 * scale[["muC"]]
    b_r <- rnorm(1, b_r_mean, b_sd); b_c <- rnorm(1, b_c_mean, b_sd)
    tt <- make_schedule(days)
    tt <- tt[runif(length(tt)) < compliance]
    avail  <- runif(length(tt)) < availability
    prompt <- avail & (runif(length(tt)) < p_prompt)
    C <- rnorm(1, 0.6, 0.6); R <- rnorm(1, -0.2, 0.6)
    rows <- vector("list", 2 * length(tt)); k <- 0; tprev <- 0
    for (j in seq_along(tt)) {
      Delta <- tt[j] - tprev
      if (Delta > 0) {                         # Euler-Maruyama from tprev to tt[j]
        n_steps <- max(1L, ceiling(steps_per_unit * Delta)); dt <- Delta / n_steps
        sq <- sqrt(dt)
        zc <- rnorm(n_steps, 0, sq) * sigma_common
        ec <- rnorm(n_steps, 0, sq) * sigma_c
        er <- rnorm(n_steps, 0, sq) * sigma_r
        for (s in seq_len(n_steps)) {
          dC <- a_c * C + beta_r_to_c * R + b_c
          dR <- a_r * R + beta_c_to_r * C + b_r
          C <- C + dC * dt + ec[s] + zc[s]
          R <- R + dR * dt + er[s] + zc[s]
        }
      }
      tprev <- tt[j]
      if (avail[j]) {                          # ratings at t (before prompt)
        art <- if (prompt[j]) demand else 0
        e <- rnorm(4, 0, me_sd)
        k <- k + 1
        rows[[k]] <- c(i, tt[j], R + art + e[1], 0.8 * R + art + e[2],
                       C + e[3], 1.1 * C + e[4], 0)
      }
      if (prompt[j]) {                         # impulse at t + eps (distal coding)
        tprev <- tt[j] + eps
        R <- R + g0 + g1 * C
        k <- k + 1
        rows[[k]] <- c(i, tt[j] + eps, NA, NA, NA, NA, 1)
      }
    }
    out[[i]] <- do.call(rbind, rows[seq_len(k)])
  }
  d <- as.data.frame(do.call(rbind, out))
  names(d) <- c("id", "time", "Y1", "Y2", "Y3", "Y4", "TD1")
  attr(d, "truth") <- list(g0 = g0_mean, g1 = g1_mean, scale = scale)
  d
}

# Within-person SDs of the latent processes without prompts (no measurement error)
calibrate_sd_mrt <- function(n = 80, seed = 99, ...) {
  d <- simulate_mrt(n_subjects = n, scale = c(C = 1, R = 1, muC = 0),
                    availability = 1, p_prompt = 0, me_sd = 0, seed = seed, ...)
  sdw <- function(x, id) sqrt(mean(tapply(x, id, var)))
  c(C = sdw(d$Y3, d$id), R = sdw(d$Y1, d$id), muC = mean(d$Y3))
}

# Calibrate once and cache (muC is hard-coded into the compiled model!)
get_scale_mrt <- function(path = "model/scale_mrt.rds", ...) {
  if (file.exists(path)) return(readRDS(path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  sc <- calibrate_sd_mrt(...)
  saveRDS(sc, path)
  sc
}

# ---- helper quantities for the proposal text --------------------------------
half_life_h <- function(a_r = -0.45) log(2) / abs(a_r) * TIME_UNIT_H
effect_next <- function(g1 = -0.20, a_r = -0.45, gap_h = 3.5) g1 * exp(a_r * gap_h / TIME_UNIT_H)
effect_vas  <- function(effect_sd, vas_sd_within = 15) effect_sd * vas_sd_within
