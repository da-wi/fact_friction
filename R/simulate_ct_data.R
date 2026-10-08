#' Simulate continuous-time longitudinal data with coupled latent processes
#'
#' @description
#' Generates synthetic longitudinal data from a simple continuous-time
#' dynamic system with two latent variables ("R" and "C") that evolve
#' continuously over irregular time intervals. At each measurement
#' occasion, subjects may receive an impulsive "nudge" that instantaneously
#' alters the latent state of the R process based on both a random intercept
#' and a random slope term. The latent trajectories are then converted into
#' observed indicators with measurement noise.
#'
#' @details
#' The system consists of two coupled latent processes (`eta_r`, `eta_c`)
#' evolving according to linear stochastic differential equations:
#'
#' \deqn{
#'   d(eta_r)/dt = a_r * eta_r + beta_c_to_r * eta_c + b_r + \varepsilon_r, \\
#'   d(eta_c)/dt = a_c * eta_c + beta_r_to_c * eta_r + b_c + \varepsilon_c.
#' }
#'
#' Random Gaussian perturbations are added at each integration step to
#' simulate process noise (`sigma_r`, `sigma_c`, `sigma_common`). Measurement
#' error is added to yield observed variables `Y1`–`Y4`. The observation
#' intervals follow a log-normal distribution with mean `dt_mean` and
#' log-SD `dt_logsd`. A binary time-dependent covariate `TD1` indicates
#' occasions when the impulsive "nudge" occurs.
#'
#' @param n_subjects Integer. Number of subjects to simulate.
#' @param Tpoints Integer. Number of measurement occasions per subject (≥ 2).
#' @param dt_mean Numeric. Mean interval length between measurements (on the original scale).
#' @param dt_logsd Numeric. Standard deviation of log-transformed intervals.
#' @param gamma_r1_mean Numeric. Mean random-slope effect for the impulsive mapping.
#' @param gamma_r1_sd Numeric. SD of the random-slope effect for impulses.
#' @param gamma_r0_mean Numeric. Mean random-intercept effect for impulses.
#' @param gamma_r0_sd Numeric. SD of the random-intercept effect for impulses.
#' @param p_nudge Numeric (0–1). Probability of receiving a nudge at each time point.
#' @param seed Optional integer. Random seed for reproducibility.
#'
#' @return
#' A tibble with simulated data containing the following columns:
#' \describe{
#'   \item{id}{Subject identifier.}
#'   \item{time}{Continuous measurement time in arbitrary units.}
#'   \item{eta_r, eta_c}{Latent process states at each observation.}
#'   \item{TD1}{Binary time-dependent covariate (nudge indicator).}
#'   \item{Y1–Y4}{Observed indicators with measurement error.}
#' }
#'
#' @examples
#' sim_data <- simulate_ct_data(n_subjects = 10, Tpoints = 50, seed = 123)
#' dplyr::glimpse(sim_data)
#'
#' @export
simulate_ct_data <- function(
    n_subjects    = 120,
    Tpoints       = 89,      # number of measurement occasions per subject
    dt_mean       = 1.0,     # mean Δt on original scale
    dt_logsd      = 0.20,    # sd of log(Δt)
    gamma_r1_mean = 0.2,
    gamma_r1_sd   = 0.05,
    gamma_r0_mean =  0,
    gamma_r0_sd   = 0.05,
    p_nudge       = 0.30,
    seed          = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  # ---- hardcoded dynamics (do not expose as args) ----
  a_r <- -0.45
  a_c <- -0.35
  beta_c_to_r <-  0.15
  beta_r_to_c <- -0.55
  b_r <- 0.00
  b_c <- 0.10
  sigma_r <- 0.20
  sigma_c <- 0.20
  sigma_common <- 0.10
  
  # measurement
  me_sd <- 0.25
  
  # integration control
  steps_per_unit <- 50L
  
  # convert (mean on original scale, log-sd) -> (meanlog, sdlog) for rlnorm
  meanlog <- log(dt_mean) - 0.5 * dt_logsd^2
  sdlog   <- dt_logsd
  
  # subject-level heterogeneity for impulse mapping
  gamma_r0_i <- rnorm(n_subjects, gamma_r0_mean, gamma_r0_sd)
  gamma_r1_i <- rnorm(n_subjects, gamma_r1_mean, gamma_r1_sd)
  
  out_list <- vector("list", n_subjects)
  
  for (id in seq_len(n_subjects)) {
    if (Tpoints < 2) stop("Tpoints must be >= 2")
    
    # irregular times per subject
    dts   <- rlnorm(Tpoints - 1L, meanlog = meanlog, sdlog = sdlog)
    times <- c(0, cumsum(dts))
    
    # initial latent states
    eta_r_state <- rnorm(1, mean = -0.2, sd = 0.6)
    eta_c_state <- rnorm(1, mean =  0.6, sd = 0.6)
    
    # TD schedule per occasion
    nudge_sched <- rbinom(Tpoints, size = 1, prob = p_nudge)
    
    gamma_r0 <- gamma_r0_i[id]
    gamma_r1 <- gamma_r1_i[id]
    
    eta_r_out <- numeric(Tpoints)
    eta_c_out <- numeric(Tpoints)
    
    for (t_idx in seq_len(Tpoints)) {
      
      # --- state-dependent impulse at the observation time ---
      if (nudge_sched[t_idx] == 1) {
        eta_r_state <- eta_r_state + (gamma_r0 + gamma_r1 * eta_r_state)
      }
      
      # record latent states at the observation
      eta_r_out[t_idx] <- eta_r_state
      eta_c_out[t_idx] <- eta_c_state
      
      # --- continuous evolution to next observation (no TD term here) ---
      if (t_idx < Tpoints) {
        Delta   <- dts[t_idx]
        n_steps <- max(1L, ceiling(steps_per_unit * Delta))
        dt_sub  <- Delta / n_steps
        
        for (s in seq_len(n_steps)) {
          d_eta_r <- a_r * eta_r_state + beta_c_to_r * eta_c_state + b_r
          d_eta_c <- a_c * eta_c_state + beta_r_to_c * eta_r_state + b_c
          
          z_common <- rnorm(1, 0, sqrt(dt_sub)) * sigma_common
          z_r      <- rnorm(1, 0, sqrt(dt_sub)) * sigma_r
          z_c      <- rnorm(1, 0, sqrt(dt_sub)) * sigma_c
          
          eta_r_state <- eta_r_state + d_eta_r * dt_sub + z_r + z_common
          eta_c_state <- eta_c_state + d_eta_c * dt_sub + z_c + z_common
        }
      }
    }
    
    out_list[[id]] <- tibble::tibble(
      id    = id,
      time  = times,
      eta_r = eta_r_out,
      eta_c = eta_c_out,
      TD1   = nudge_sched
    )
  }
  
  sim_latent <- dplyr::bind_rows(out_list)
  
  # measurement (fixed loadings to match DGP)
  sim_latent |>
    dplyr::mutate(
      Y1 = eta_r + rnorm(dplyr::n(), 0, me_sd),
      Y2 = 0.8 * eta_r + rnorm(dplyr::n(), 0, me_sd),
      Y3 = eta_c + rnorm(dplyr::n(), 0, me_sd),
      Y4 = 1.1 * eta_c + rnorm(dplyr::n(), 0, me_sd)
    )
}


generateSimData <- function(n_subjects = 120, t_points = 70, eff_L_inter = -0.2, 
                            seed = NULL, 
                            # Alle Parameter als Default-Werte definieren:
                            d_LL = -0.5, d_LA = 0.1, d_AL = 0.1, d_AA = -0.5,
                            cint_L = 0.5, cint_A = 0.5,
                            eff_L_base = -0.05, eff_A_base = 0.2, eff_A_inter = 0.1,
                            decay_Reset = -0.7, diff_L = 1, diff_A = 1,
                            dt_sim = 0.05) {
  
  if(!is.null(seed)) set.seed(seed)
  datalist <- list() # Initialisierung innerhalb der Funktion
  
  for(i in 1:n_subjects) {
    # 1. Zeitstruktur
    intervals <- exp(rnorm(t_points-1, log(1), 0.2))
    times <- c(0, cumsum(intervals))
    mrt_trigger <- rbinom(t_points, 1, 0.3)
    
    # 2. Startwerte (Load=1, Agency=1, Reset=0)
    latent_state <- c(1.0, 1.0, 0) 
    
    subject_data <- matrix(NA, nrow=t_points, ncol=5)
    colnames(subject_data) <- c("id", "time", "mental_load", "agency_control", "mrt_trigger")
    
    for(t_idx in 1:t_points) {
      # MRT-Impuls (Pulse)
      if(mrt_trigger[t_idx] == 1) {
        latent_state[3] <- latent_state[3] + 1.0 
      }
      
      # Messung (Manifeste Variablen)
      subject_data[t_idx,] <- c(
        i, times[t_idx], 
        latent_state[1] + rnorm(1, 0, 0.1), 
        latent_state[2] + rnorm(1, 0, 0.1), 
        mrt_trigger[t_idx]
      )
      
      # 3. Zeit-Integration (Zwischen den Messpunkten)
      if(t_idx < t_points) {
        wait <- times[t_idx+1] - times[t_idx]
        steps <- ceiling(wait / dt_sim)
        actual_dt <- wait / steps
        
        for(s in 1:steps) {
          # Ableitungen berechnen
          dL <- (d_LL * latent_state[1] + d_LA * latent_state[2] + cint_L) +
            (eff_L_base + eff_L_inter * latent_state[1]) * latent_state[3]
          
          dA <- (d_AL * latent_state[1] + d_AA * latent_state[2] + cint_A) +
            (eff_A_base + eff_A_inter * latent_state[1]) * latent_state[3]
          
          dR <- decay_Reset * latent_state[3]
          
          # Update mit skaliertem Rauschen (Diffusion * sqrt(dt))
          latent_state[1] <- latent_state[1] + dL * actual_dt + rnorm(1, 0, diff_L * sqrt(actual_dt))
          latent_state[2] <- latent_state[2] + dA * actual_dt + rnorm(1, 0, diff_A * sqrt(actual_dt))
          latent_state[3] <- latent_state[3] + dR * actual_dt
        }
      }
    }
    datalist[[i]] <- as.data.frame(subject_data)
  }
  
  df_sim <- do.call(rbind, datalist)
  return(df_sim) 
}


generateSimData <- function(n_subjects = 120, t_points = 70, eff_L_inter = -0.2, 
                            time_variance = "High", missing_type = "None",
                            seed = NULL, 
                            d_LL = -0.5, d_LA = 0.1, d_AL = 0.1, d_AA = -0.5,
                            cint_L = 0.5, cint_A = 0.5,
                            eff_L_base = -0.05, eff_A_base = 0.2, eff_A_inter = 0.1,
                            decay_Reset = -0.7, diff_L = 1.0, diff_A = 1.0,
                            dt_sim = 0.05) {
  
  if(!is.null(seed)) set.seed(seed)
  datalist <- list() 
  
  lambda_L2 <- 0.9   # Loading für den 2. Indikator von Load
  lambda_A2 <- 0.9   # Loading für den 2. Indikator von Agency
  me_sd <- 0.4       # Messfehler SD (0.4^2 = 0.16 Varianz -> Reliabilität ca .86)
  
  for(i in 1:n_subjects) {
    if(time_variance == "Low") {
      intervals <- rnorm(t_points - 1, mean = 1.0, sd = 0.05)
    } else {
      intervals <- exp(rnorm(t_points - 1, log(1), 0.4))
      night_indices <- seq(5, t_points - 1, by = 5)
      intervals[night_indices] <- 3.0
    }
    
    times <- c(0, cumsum(intervals))
    mrt_trigger <- rbinom(t_points, 1, 0.3)
    latent_state <- c(1.0, 1.0, 0) 
    
    # Matrix um 2 Spalten erweitert
    subject_data <- matrix(NA, nrow=t_points, ncol=7)
    colnames(subject_data) <- c("id", "time", "load1", "load2", "agency1", "agency2", "mrt_trigger")
    
    for(t_idx in 1:t_points) {
      if(mrt_trigger[t_idx] == 1) latent_state[3] <- latent_state[3] + 1.0 
      
      current_load <- latent_state[1]
      p_miss <- 0
      
      if(missing_type == "MCAR") {
        p_miss <- 0.20
      } else if(missing_type == "MNAR") {
        logit_p <- -3.5 + 2.0 * current_load
        p_miss <- 1 / (1 + exp(-logit_p))
      }
      
      is_missing <- rbinom(1, size = 1, prob = p_miss) == 1
      
      if(is_missing) {
        m_load1 <- NA; m_load2 <- NA
        m_agency1 <- NA; m_agency2 <- NA
      } else {
        # NEU: 2 Indikatoren pro Konstrukt generieren
        m_load1   <- current_load * 1.0       + rnorm(1, 0, me_sd)
        m_load2   <- current_load * lambda_L2 + rnorm(1, 0, me_sd)
        
        m_agency1 <- latent_state[2] * 1.0       + rnorm(1, 0, me_sd)
        m_agency2 <- latent_state[2] * lambda_A2 + rnorm(1, 0, me_sd)
      }
      
      subject_data[t_idx,] <- c(i, times[t_idx], 
                                m_load1, m_load2, m_agency1, m_agency2, 
                                mrt_trigger[t_idx])
      
      if(t_idx < t_points) {
        wait <- times[t_idx+1] - times[t_idx]
        steps <- ceiling(wait / dt_sim)
        actual_dt <- wait / steps
        for(s in 1:steps) {
          dL <- (d_LL * latent_state[1] + d_LA * latent_state[2] + cint_L) +
            (eff_L_base + eff_L_inter * latent_state[1]) * latent_state[3]
          dA <- (d_AL * latent_state[1] + d_AA * latent_state[2] + cint_A) +
            (eff_A_base + eff_A_inter * latent_state[1]) * latent_state[3]
          dR <- decay_Reset * latent_state[3]
          
          latent_state[1] <- latent_state[1] + dL * actual_dt + rnorm(1, 0, diff_L * sqrt(actual_dt))
          latent_state[2] <- latent_state[2] + dA * actual_dt + rnorm(1, 0, diff_A * sqrt(actual_dt))
          latent_state[3] <- latent_state[3] + dR * actual_dt
        }
      }
    }
    datalist[[i]] <- as.data.frame(subject_data)
  }
  return(do.call(rbind, datalist)) 
}
