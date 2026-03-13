# =============================================
# Agent-Based Model: School and Household Transmission
# =============================================

library(tidyverse)
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
  repeat {
    hh_sizes <- rpois(num_HH, lambda = 2)
    hh_sizes[hh_sizes == 0] <- 1

    # Adjust total to exactly num_students
    while (sum(hh_sizes) != num_students) {
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
  class_assignment <- sample(rep(1:num_classes, each = num_students %/% num_classes))

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
# Step 2: Classroom transmission subfunction
# =============================================
transmit_classroom <- function(idx, df, class_members_list, base_rate) {
  cls <- as.character(df$class[idx])
  classmates <- class_members_list[[cls]]
  classmates <- classmates[classmates != idx]
  susc <- classmates[df$susp[classmates] == 1L]
  n_susc <- length(susc)

  new_infected <- integer(0)
  df$track_school_contacts[idx] <- df$track_school_contacts[idx] + n_susc

  if (n_susc > 0) {
    infections <- rbinom(n_susc, 1, base_rate)
    new_infected <- susc[infections == 1]
    df$secondary_infs_school[idx] <- df$secondary_infs_school[idx] + length(new_infected)
  }

  return(list(df = df, new_infected = new_infected))
}

# =============================================
# Step 3: Household transmission subfunction
# =============================================
transmit_household <- function(idx, df, hh_members_list, hh_rate) {
  hh <- as.character(df$HH_ID[idx])
  hh_mem <- hh_members_list[[hh]]
  hh_mem <- hh_mem[hh_mem != idx]
  susc <- hh_mem[df$susp[hh_mem] == 1L]
  n_susc <- length(susc)

  new_infected <- integer(0)
  df$track_home_contacts[idx] <- df$track_home_contacts[idx] + n_susc

  if (n_susc > 0) {
    infections <- rbinom(n_susc, 1, hh_rate)
    new_infected <- susc[infections == 1]
    df$secondary_infs_home[idx] <- df$secondary_infs_home[idx] + length(new_infected)
  }

  return(list(df = df, new_infected = new_infected))
}

# =============================================
# Step 4: Simulation function
# =============================================
run_model <- function(df, t, rel_HH, seeds) {
  n <- nrow(df)

  # Reset infection-related columns
  df$susp <- 1L
  df$exposed <- NA_real_
  df$start <- NA_real_
  df$end <- NA_real_
  df$seed <- 0L
  df$track_inf_days <- 0
  df$track_latent_days <- 0
  df$track_home_contacts <- 0L
  df$secondary_infs_home <- 0L
  df$track_school_contacts <- 0L
  df$secondary_infs_school <- 0L

  # Pre-compute group memberships for speed
  student_classes <- unique(df$class[df$type == "Student"])
  class_members_list <- lapply(student_classes, function(c) which(df$class == c & df$type == "Student"))
  names(class_members_list) <- as.character(student_classes)

  all_hh <- unique(df$HH_ID)
  hh_members_list <- lapply(all_hh, function(h) which(df$HH_ID == h))
  names(hh_members_list) <- as.character(all_hh)

  # Seed infections
  seed_idx <- sample(1:n, seeds)
  for (s in seed_idx) {
    L <- rpois(1, latent_days)
    I <- rpois(1, infectious_days)
    df$susp[s] <- 0L
    df$exposed[s] <- 1
    df$start[s] <- 1 + L
    df$end[s] <- 1 + L + I
    df$seed[s] <- 1L
  }

  # Day-by-day simulation
  for (day in 1:t) {
    is_weekday <- ((day - 1) %% 7) < 5  # Day 1 = Monday

    # Find currently infectious individuals
    infectious <- which(df$susp == 0L & !is.na(df$start) &
                          df$start <= day & df$end > day)

    if (length(infectious) == 0) next

    for (idx in infectious) {
      # Draw base transmission rate ~ Uniform(0, 0.02)
      base_rate <- runif(1, 0, 0.02)
      hh_rate <- rel_HH * base_rate

      # Classroom transmission (weekdays only, students only)
      if (is_weekday && df$type[idx] == "Student") {
        result <- transmit_classroom(idx, df, class_members_list, base_rate)
        df <- result$df

        for (new_id in result$new_infected) {
          L <- rpois(1, latent_days)
          I <- rpois(1, infectious_days)
          df$susp[new_id] <- 0L
          df$exposed[new_id] <- day
          df$start[new_id] <- day + L
          df$end[new_id] <- day + L + I
        }
      }

      # Household transmission (all days)
      result_hh <- transmit_household(idx, df, hh_members_list, hh_rate)
      df <- result_hh$df

      for (new_id in result_hh$new_infected) {
        L <- rpois(1, latent_days)
        I <- rpois(1, infectious_days)
        df$susp[new_id] <- 0L
        df$exposed[new_id] <- day
        df$start[new_id] <- day + L
        df$end[new_id] <- day + L + I
      }
    }
  }

  # Compute tracked infectious and latent days (censored at simulation end)
  infected <- which(df$susp == 0L)
  if (length(infected) > 0) {
    df$track_inf_days[infected] <- pmax(0,
      pmin(df$end[infected], t + 1) - pmax(df$start[infected], 1))
    df$track_latent_days[infected] <- pmax(0,
      pmin(df$start[infected], t + 1) - pmax(df$exposed[infected], 1))
  }

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
