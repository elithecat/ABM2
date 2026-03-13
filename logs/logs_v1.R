# =============================================
# Log: Model v1 Generation Process
# =============================================
# Date: 2026-03-13
# Model: Agent-Based School/Household Transmission Model
# Output: /scripts/model_v1.R
# =============================================

# ----- Step 1: Input Review -----
# Read specification from GPT_structured_command.md
# Key parameters extracted:
#   - 500 students, 222 households, 2 adults per household (944 total individuals)
#   - 25 classes of 20 students each
#   - Students assigned to HH via Poisson(lambda=2), min 1 per HH
#   - Transmission rate per classroom contact: Uniform(0, 0.02), mean = 0.01
#   - Household transmission rate: 2x classroom rate
#   - Incubation period: Poisson(2 days)
#   - Infectious period: Poisson(5 days)
#   - 5 seed infections, 30-day simulation
#   - School mixing Mon-Fri only; household mixing all 7 days

# ----- Step 2: Test Script Analysis -----
# Read test_script.R to identify required outputs and validation criteria
# Required global variables:
#   all, model_run_length, rel_transmission_HH, num_seeds_inf,
#   infectious_days, latent_days, transmission_prob_school,
#   num_1studentHH, num_2studentHH, num_3studentHH, num_4studentHH,
#   num_HH, num_students, num_classes, num_parents_perHH
# Required function:
#   run_model(df, t, rel_HH, seeds) -> data.frame
# Required output columns:
#   id, type, HH_ID, class, susp, exposed, start, end, seed,
#   track_inf_days, track_latent_days, track_home_contacts,
#   secondary_infs_home, track_school_contacts, secondary_infs_school, model_days
# Validation: meets.threshold (relative_difference < 0.015) should be TRUE

# ----- Step 3: Population Initialization (create_population) -----
# Created function to generate population data frame:
#   - Draw HH sizes from Poisson(2), enforce min 1, adjust total to 500
#   - Assign students to 25 classes of exactly 20
#   - Create 2 adults per HH with class = 0 (single unique value for type grouping)
#   - Initialize all individuals as susceptible (susp = 1)
# Output: data.frame with 944 rows (500 students + 444 adults)

# ----- Step 4: Classroom Transmission Subfunction -----
# transmit_classroom(idx, df, class_members_list, base_rate):
#   - Finds susceptible classmates of infectious student idx
#   - Draws Bernoulli(base_rate) for each susceptible contact
#   - Tracks: track_school_contacts, secondary_infs_school
#   - Returns updated df and vector of newly infected IDs
# Transmission rate: base_rate ~ Uniform(0, 0.02), mean = 0.01
# Only called on weekdays (Mon-Fri) for students

# ----- Step 5: Household Transmission Subfunction -----
# transmit_household(idx, df, hh_members_list, hh_rate):
#   - Finds susceptible household members of infectious individual idx
#   - Draws Bernoulli(hh_rate) for each susceptible contact
#   - hh_rate = rel_HH * base_rate (2x classroom rate)
#   - Tracks: track_home_contacts, secondary_infs_home
#   - Returns updated df and vector of newly infected IDs
# Called on all days (Mon-Sun) for all individual types

# ----- Step 6: Simulation Function (run_model) -----
# run_model(df, t, rel_HH, seeds):
#   - Resets all infection columns (enables reuse of same df across runs)
#   - Pre-computes class and household membership lists for performance
#   - Seeds 5 random infections on day 1 (exposed=1, start=1+Pois(2), end=start+Pois(5))
#   - Simulates day-by-day for t=30 days:
#     * Day 1 = Monday; weekday = ((day-1) %% 7) < 5
#     * For each infectious individual (start <= day < end):
#       - Draw base_rate ~ Uniform(0, 0.02)
#       - If weekday & student: classroom transmission at base_rate
#       - Household transmission at rel_HH * base_rate (all days)
#       - New infections: exposed=day, start=day+Pois(2), end=start+Pois(5)
#   - Computes track_inf_days and track_latent_days (censored at simulation boundary)
#     * track_inf_days = max(0, min(end, t+1) - max(start, 1))
#     * track_latent_days = max(0, min(start, t+1) - max(exposed, 1))
#   - Returns complete data.frame of all individuals

# ----- Step 7: Derived Parameters -----
# After population creation, computed:
#   num_1studentHH, num_2studentHH, num_3studentHH, num_4studentHH
# These are counts of HHs with exactly 1, 2, 3, or 4 students

# ----- Step 8: Validation Design -----
# Test script runs 10,000 iterations via foreach/doMC parallel
# For each run, computes summary statistics by type (Student, Household)
# Key validation targets (must be within 1.5% relative difference):
#   - attack_rate_student ≈ 0.01 (mean of Uniform(0,0.02))
#   - attack_rate_HH ≈ 0.02 (2x classroom rate)
#   - rel_HH_obs ≈ 2.0
#   - track_inf_days ≈ 5 (Poisson mean, after censoring adjustment)
#   - track_latent_days ≈ 2 (Poisson mean, after censoring adjustment)
#   - num_indv: 500 (Student), 444 (Household)
#   - num_hh: 222
#   - num_class: 25 (Student), 1 (Household)
#   - avg_seeds_total: 5
#   - model_days: 30

# ----- Step 9: Output -----
# Model saved to: /scripts/model_v1.R
# Log saved to: /logs/logs_v1.R
# All files committed and pushed to GitHub
