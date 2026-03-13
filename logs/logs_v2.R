# =============================================
# Log: Model v2 Generation Process
# =============================================
# Date: 2026-03-13
# Model: Agent-Based School/Household Transmission Model (v2)
# Input: model_v1.R + three agent review reports
# Output: /scripts/model_v2.R
# =============================================

# ----- Step 1: Review Agent Reports -----
# Three agent reviews were conducted on model_v1.R:
#   1. Code Review (agents/code-review-report.md)
#   2. Epidemiological Model Review (agents/epi-model-review-report.md)
#   3. Performance/Speed Review (agents/code-speed-report.md)
#
# Key findings prioritized for v2:
#   Critical: Hardcoded transmission bounds, zero-length periods,
#     fragile class assignment, processing order bias, copy-on-modify bottleneck
#   Major: Global variable dependencies, sample() pitfall, tidyverse overhead
#   Minor: Variable shadowing, missing iteration guards

# ----- Step 2: Changes Implemented -----

# [Code Review #1] Parameterized transmission rate bounds
# BEFORE: runif(1, 0, 0.02)  -- hardcoded, disconnected from transmission_prob_school
# AFTER:  runif(1, 0, 2 * transmission_prob_school)  -- tied to parameter
# WHY:    Changing transmission_prob_school without updating 0.02 would silently break model
# NOTE:   U(0, 0.02) is per specification; parameterization preserves spec behavior

# [Code Review #2] Documented zero-length Poisson periods (NOT changed)
# Infectious period I ~ Pois(5): P(I=0) = 0.67%. Person infected but never infectious.
# Latent period L ~ Pois(2): P(L=0) = 13.5%. Immediate onset of infectiousness.
# DECISION: Kept Poisson per spec requirement for equal mean and variance.
#   Using max(1, Pois(5)) would shift mean to ~5.007 (0.13%) but reduce variance
#   to ~4.94 (1.2%), violating equal mean=variance requirement and risking test failure.
#   For latent: max(1, Pois(2)) would shift mean by 6.75%, far exceeding 1.5% threshold.
# Added documentation comments in model code.

# [Code Review #3] Robust class assignment
# BEFORE: sample(rep(1:num_classes, each = num_students %/% num_classes))
# AFTER:  sample(rep(1:num_classes, length.out = num_students))
# WHY:    Integer division silently drops students when not evenly divisible

# [Code Review #4] Documented seed selection from entire population
# Seeds drawn from all 944 individuals per spec Step 4: "randomly initialize 5
# seed infections from all individuals." Added comment explaining that seeded
# adults can only transmit in households, not classrooms.

# [Code Review #7] Parameterized run_model with latent/infectious days
# BEFORE: run_model(df, t, rel_HH, seeds) -- used global latent_days, infectious_days
# AFTER:  run_model(df, t, rel_HH, seeds, lat_days=latent_days, inf_days=infectious_days)
# WHY:    Makes function self-contained; defaults maintain backward compatibility

# [Code Review #11] Added iteration guard in create_population
# BEFORE: repeat {} with no max iteration
# AFTER:  Outer loop max 1000, inner adjustment loop max 10000
# WHY:    Prevents infinite loop with pathological parameter combinations

# [Epi Review C2] Randomized infectious processing order each day
# BEFORE: which() returns ascending order -- students always processed first
# AFTER:  sample() with length>1 guard to randomize order
# WHY:    Eliminates systematic ID-based bias in who gets "credit" for infections
# BUG FOUND: Initial implementation used sample(infectious) without guarding
#   against R's sample(n) pitfall: when infectious has length 1 and value k,
#   sample(k) returns random permutation of 1:k instead of c(k). This caused
#   processing of non-infectious individuals, inflating attack rates by ~2.5%
#   and doubling runtime from 22s to 52s. Fixed with length>1 guard.

# [Speed #1+3] Extracted vectors, inlined transmission functions
# BEFORE: transmit_classroom() and transmit_household() received/returned full df
#   Each call triggered copy-on-modify on the data.frame (2-4 column copies per call)
# AFTER:  All mutable state in plain vectors (refcount=1, no copies triggered)
#   Transmission logic inlined directly in the day loop
# WHY:    Eliminates hundreds of unnecessary vector copies per simulation run
# RESULT: Primary performance bottleneck removed

# [Speed #2] Vectorized RNG calls
# BEFORE: rpois(1, ...) called in a loop for each newly infected individual
# AFTER:  rpois(length(new_inf), ...) called once per batch of new infections
# WHY:    Reduces per-call R overhead for argument checking and dispatch

# [Speed #4] Used split() for membership lists
# BEFORE: lapply(groups, function(g) which(col == g)) -- O(N x G) scans
# AFTER:  split(indices, groups) -- single O(N) pass
# WHY:    Cleaner and asymptotically faster (negligible per-run since called once)

# [Speed #5] Removed tidyverse dependency
# BEFORE: library(tidyverse) loaded ~25 packages including ggplot2, readr, etc.
# AFTER:  library(data.table) only -- no tidyverse functions used in model
# WHY:    Reduces worker startup time by 1-3s and memory by ~50MB per worker
# NOTE:   test_script.R still loads tidyverse for mutate() in summary step

# [Speed #10] Initialize mutable vectors directly
# BEFORE: df$col <- value triggered copy-on-modify on the full data.frame
# AFTER:  rep(initial_value, n) creates independent vectors with refcount=1
# WHY:    Avoids the initial data.frame copy at start of each run_model() call

# ----- Step 3: Validation Results -----
# Test: source("scripts/model_v2.R"); source("test_script.R")
# Trials: 10,000 parallel runs
# Runtime: ~21 seconds (vs v1 ~22 seconds, plus no copy-on-modify overhead)
#
# Student table: 18/18 meets.threshold TRUE
# Parents table: 13/13 non-NaN meets.threshold TRUE (3 NaN from 0/0 inherent to design)
#
# Key metrics (representative run):
#   track_inf_days:       5.00 (target 5, diff <0.1%)
#   track_latent_days:    2.00 (target 2, diff <0.1%)
#   attack_rate_student:  0.010 (target 0.01, diff <1.3%)
#   attack_rate_HH:       0.020 (target 0.02, diff <0.5%)
#   rel_HH_obs:           2.00 (target 2, diff <1.2%)

# ----- Step 4: Issues NOT addressed in v2 -----
# The following were noted but intentionally deferred:
#
# [Code Review #8] Variable name `t` shadows base R transpose
#   REASON: Cannot change without breaking test_script.R compatibility
#
# [Code Review #9] Variable name `all` shadows base R function
#   REASON: Cannot change without breaking test_script.R compatibility
#
# [Epi W1] Poisson CV is biologically unusual for infectious periods
#   REASON: Explicitly required by specification
#
# [Epi W3] Household rate ignores household size effects
#   REASON: Per specification, rate is per-contact regardless of HH size
#
# [Epi S2] No weekday/weekend breakdown of secondary infections
#   REASON: test_script.R does not check this; would require new columns
#
# [Speed #6] foreach .combine optimization
#   REASON: Requires changes to test_script.R which is not modifiable
#
# [Speed #9] Rcpp implementation for inner loop
#   REASON: High effort; current performance is adequate (~21s for 10k runs)

# ----- Step 5: Output -----
# Model saved to: /scripts/model_v2.R
# Log saved to: /logs/logs_v2.R
# Agent reports in: /agents/
