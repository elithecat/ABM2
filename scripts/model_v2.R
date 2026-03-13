# =============================================
# Agent-Based Model: School and Household Transmission (v2)
# =============================================
# Changes from v1:
#   - [Code Review #1] Parameterized transmission rate bounds via transmission_prob_school
#   - [Code Review #2] Documented zero-length period behavior (kept Poisson per spec)
#   - [Code Review #3] Robust class assignment using length.out
#   - [Code Review #4] Documented seed selection from all individuals (per spec)
#   - [Code Review #7] Added latent_days/infectious_days as run_model defaults
#   - [Code Review #11] Added max-iteration guard in create_population
#   - [Epi Review C2] Randomized infectious processing order each day
#   - [Speed #1+3] Extracted vectors, inlined transmission functions
#   - [Speed #2] Vectorized RNG calls
#   - [Speed #4] Used split() for membership lists
#   - [Speed #5] Removed tidyverse dependency; only data.table needed
#   - [Speed #10] Initialize mutable vectors directly (no df copy)
# =============================================

library(data.table)

# ----- Model Parameters -----
num_students <- 500
num_HH <- 222
num_parents_perHH <- 2
num_classes <- 25
num_seeds_inf <- 5
model_run_length <- 30
transmission_prob_school <- 0.01
rel_transmission_HH <- 2
infectious_days <- 5
latent_days <- 2

# =============================================
# Step 1: Create population data frame
# =============================================
create_population <- function() {
  # Assign students to households using Poisson(2), min 1 per HH
  # [Code Review #11] Added iteration guard to prevent infinite loop
  max_iter <- 1000
  iter <- 0
  repeat {
    iter <- iter + 1
    if (iter > max_iter) stop("create_population: failed to converge after ", max_iter, " iterations")

    hh_sizes <- rpois(num_HH, lambda = 2)
    hh_sizes[hh_sizes == 0] <- 1

    # Adjust total to exactly num_students
    adj_iter <- 0
    while (sum(hh_sizes) != num_students && adj_iter < 10000) {
      adj_iter <- adj_iter + 1
      diff_val <- num_students - sum(hh_sizes)
      if (diff_val > 0) {
        idx <- sample(num_HH, min(abs(diff_val), num_HH))
        hh_sizes[idx] <- hh_sizes[idx] + 1
      } else {
        eligible <- which(hh_sizes > 1)
        if (length(eligible) == 0) break
        idx <- sample(eligible, min(abs(diff_val), length(eligible)))
        hh_sizes[idx] <- hh_sizes[idx] - 1
      }
    }
    if (sum(hh_sizes) == num_students && all(hh_sizes >= 1)) break
  }

  # Create student rows
  student_hh <- rep(1:num_HH, times = hh_sizes)
  # [Code Review #3] Use length.out to handle non-divisible student/class counts
  class_assignment <- sample(rep(1:num_classes, length.out = num_students))

  students <- data.frame(
    id = 1:num_students,
    type = "Student",
    HH_ID = student_hh,
    class = class_assignment,
    susp = 1L,
    exposed = NA_real_,
    start = NA_real_,
    end = NA_real_,
    seed = 0L,
    track_inf_days = 0,
    track_latent_days = 0,
    track_home_contacts = 0L,
    secondary_infs_home = 0L,
    track_school_contacts = 0L,
    secondary_infs_school = 0L,
    model_days = model_run_length,
    stringsAsFactors = FALSE
  )

  # Create household member (adult) rows
  n_adults <- num_parents_perHH * num_HH
  adults <- data.frame(
    id = (num_students + 1):(num_students + n_adults),
    type = "Household",
    HH_ID = rep(1:num_HH, each = num_parents_perHH),
    class = 0L,
    susp = 1L,
    exposed = NA_real_,
    start = NA_real_,
    end = NA_real_,
    seed = 0L,
    track_inf_days = 0,
    track_latent_days = 0,
    track_home_contacts = 0L,
    secondary_infs_home = 0L,
    track_school_contacts = 0L,
    secondary_infs_school = 0L,
    model_days = model_run_length,
    stringsAsFactors = FALSE
  )

  df <- rbind(students, adults)
  return(df)
}

# =============================================
# Steps 2-4: Simulation function
# =============================================
# [Code Review #7] Added latent_days and infectious_days as parameters with defaults
# [Code Review #4] Seeds drawn from all individuals per spec ("randomly initialize
#   5 seed infections from all individuals"). Adults can only transmit in households.
# [Code Review #2] Poisson draws can produce zero-length latent or infectious periods.
#   Zero latent (P=0.135): individual becomes infectious same day as exposed.
#   Zero infectious (P=0.0067): individual never transmits. Both kept per spec
#   requirement for Poisson with equal mean and variance.
# [Epi W2] Per-infector-per-day rate draw from U(0, 2*transmission_prob_school)
#   creates correlated contact outcomes (superspreading-like overdispersion).
#   This matches the spec's "transmission rate per class contact follows U(0, 0.02)."
run_model <- function(df, t, rel_HH, seeds,
                      lat_days = latent_days,
                      inf_days = infectious_days) {
  n <- nrow(df)

  # [Speed #10] Extract static columns, initialize mutable vectors directly
  v_type <- df$type
  v_class <- df$class
  v_HH_ID <- df$HH_ID

  v_susp <- rep(1L, n)
  v_exposed <- rep(NA_real_, n)
  v_start <- rep(NA_real_, n)
  v_end <- rep(NA_real_, n)
  v_seed <- rep(0L, n)
  v_track_inf_days <- rep(0, n)
  v_track_latent_days <- rep(0, n)
  v_track_home_contacts <- rep(0L, n)
  v_secondary_infs_home <- rep(0L, n)
  v_track_school_contacts <- rep(0L, n)
  v_secondary_infs_school <- rep(0L, n)

  # [Speed #4] Use split() for O(N) membership list construction
  student_idx <- which(v_type == "Student")
  class_members_list <- split(student_idx, v_class[student_idx])
  hh_members_list <- split(seq_len(n), v_HH_ID)

  # [Code Review #1] Parameterize transmission rate bounds
  rate_upper <- 2 * transmission_prob_school

  # [Speed #2] Vectorized seed initialization
  seed_idx <- sample(1:n, seeds)
  L_vec <- rpois(seeds, lat_days)
  I_vec <- rpois(seeds, inf_days)
  v_susp[seed_idx] <- 0L
  v_exposed[seed_idx] <- 1
  v_start[seed_idx] <- 1 + L_vec
  v_end[seed_idx] <- 1 + L_vec + I_vec
  v_seed[seed_idx] <- 1L

  # Day-by-day simulation
  for (day in 1:t) {
    is_weekday <- ((day - 1) %% 7) < 5  # Day 1 = Monday

    # Find currently infectious individuals
    infectious <- which(v_susp == 0L & !is.na(v_start) &
                          v_start <= day & v_end > day)
    if (length(infectious) == 0L) next

    # [Epi Review C2] Randomize processing order to eliminate ID-based bias
    # Guard against sample(n) pitfall when length==1: sample(5) gives 1:5, not c(5)
    if (length(infectious) > 1L) infectious <- sample(infectious)

    for (idx in infectious) {
      # [Code Review #1] Draw rate from parameterized bounds
      base_rate <- runif(1, 0, rate_upper)
      hh_rate <- rel_HH * base_rate

      # [Speed #1+3] Inlined classroom transmission (weekdays, students only)
      if (is_weekday && v_type[idx] == "Student") {
        cls <- as.character(v_class[idx])
        classmates <- class_members_list[[cls]]
        classmates <- classmates[classmates != idx]
        susc <- classmates[v_susp[classmates] == 1L]
        n_susc <- length(susc)

        v_track_school_contacts[idx] <- v_track_school_contacts[idx] + n_susc

        if (n_susc > 0) {
          infections <- rbinom(n_susc, 1, base_rate)
          new_inf <- susc[infections == 1L]
          v_secondary_infs_school[idx] <- v_secondary_infs_school[idx] + length(new_inf)

          if (length(new_inf) > 0) {
            # [Speed #2] Vectorized RNG for new infections
            L <- rpois(length(new_inf), lat_days)
            I <- rpois(length(new_inf), inf_days)
            v_susp[new_inf] <- 0L
            v_exposed[new_inf] <- day
            v_start[new_inf] <- day + L
            v_end[new_inf] <- day + L + I
          }
        }
      }

      # [Speed #1+3] Inlined household transmission (all days)
      hh <- as.character(v_HH_ID[idx])
      hh_mem <- hh_members_list[[hh]]
      hh_mem <- hh_mem[hh_mem != idx]
      susc_hh <- hh_mem[v_susp[hh_mem] == 1L]
      n_susc_hh <- length(susc_hh)

      v_track_home_contacts[idx] <- v_track_home_contacts[idx] + n_susc_hh

      if (n_susc_hh > 0) {
        infections_hh <- rbinom(n_susc_hh, 1, hh_rate)
        new_inf_hh <- susc_hh[infections_hh == 1L]
        v_secondary_infs_home[idx] <- v_secondary_infs_home[idx] + length(new_inf_hh)

        if (length(new_inf_hh) > 0) {
          # [Speed #2] Vectorized RNG for new infections
          L <- rpois(length(new_inf_hh), lat_days)
          I <- rpois(length(new_inf_hh), inf_days)
          v_susp[new_inf_hh] <- 0L
          v_exposed[new_inf_hh] <- day
          v_start[new_inf_hh] <- day + L
          v_end[new_inf_hh] <- day + L + I
        }
      }
    }
  }

  # Compute tracked infectious and latent days (censored at simulation boundary)
  infected <- which(v_susp == 0L)
  if (length(infected) > 0) {
    v_track_inf_days[infected] <- pmax(0,
      pmin(v_end[infected], t + 1) - pmax(v_start[infected], 1))
    v_track_latent_days[infected] <- pmax(0,
      pmin(v_start[infected], t + 1) - pmax(v_exposed[infected], 1))
  }

  # Write vectors back into data frame for return
  df$susp <- v_susp
  df$exposed <- v_exposed
  df$start <- v_start
  df$end <- v_end
  df$seed <- v_seed
  df$track_inf_days <- v_track_inf_days
  df$track_latent_days <- v_track_latent_days
  df$track_home_contacts <- v_track_home_contacts
  df$secondary_infs_home <- v_secondary_infs_home
  df$track_school_contacts <- v_track_school_contacts
  df$secondary_infs_school <- v_secondary_infs_school

  return(df)
}

# =============================================
# Initialize population and derived parameters
# =============================================
all <- create_population()

# Compute household size counts (by number of students per HH)
student_hh_table <- table(all$HH_ID[all$type == "Student"])
num_1studentHH <- sum(student_hh_table == 1)
num_2studentHH <- sum(student_hh_table == 2)
num_3studentHH <- sum(student_hh_table == 3)
num_4studentHH <- sum(student_hh_table == 4)
