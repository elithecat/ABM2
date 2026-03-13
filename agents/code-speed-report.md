# Performance Analysis: model_v1.R + test_script.R

## Executive Summary

The model simulates 30 days of disease transmission across 944 individuals (500 students + 444 adults), run 10,000 times via `foreach %dopar%`. The dominant bottleneck is `run_model()`, which contains a **day x infectious-individuals nested loop** with per-individual data.frame mutations, transmission function calls, and scalar RNG draws. Secondary costs come from copy-on-modify semantics, unnecessary library loading, and redundant computations in `test_script.R`.

**Estimated total speedup from all recommendations: 10x-30x.**

---

## Finding 1: Copy-on-Modify in the Hot Loop (CRITICAL)

**Location:** `run_model()`, lines 189-221; `transmit_classroom()`, lines 106-114; `transmit_household()`, lines 128-136.

**Problem:** The transmission functions receive the entire data.frame `df`, mutate a single cell (`df$track_school_contacts[idx]`, `df$secondary_infs_school[idx]`, etc.), then **return the modified data.frame**. Every `$<-` assignment on a data.frame that has multiple references triggers a **full column copy** due to R's copy-on-modify semantics. The pattern:

```r
result <- transmit_classroom(idx, df, class_members_list, base_rate)
df <- result$df   # <-- df now points to a copy
```

This means each infectious individual on each day triggers **2-4 full column copies** (one per modified column in each transmission function). With ~5-50 infectious individuals per day over 30 days, that is **hundreds to thousands of unnecessary copies** per `run_model()` call, each copying a 944-element vector.

**Fix:** Extract the tracked columns into plain vectors before the loop, operate on vectors directly (no copy-on-modify since they have a single reference), then write back at the end. Alternatively, eliminate the sub-functions and inline the logic.

```r
# Before the day loop, extract vectors:
v_susp <- df$susp
v_exposed <- df$exposed
v_start <- df$start
v_end <- df$end
v_seed <- df$seed
v_track_home_contacts <- df$track_home_contacts
v_secondary_infs_home <- df$secondary_infs_home
v_track_school_contacts <- df$track_school_contacts
v_secondary_infs_school <- df$secondary_infs_school
v_type <- df$type
v_class <- df$class
v_HH_ID <- df$HH_ID

# In the loop, operate on vectors:
v_track_school_contacts[idx] <- v_track_school_contacts[idx] + n_susc
# No copy triggered -- v_susp has refcount 1

# After the loop, write back:
df$susp <- v_susp
# ... etc.
```

**Estimated speedup from this fix alone: 3x-8x** (eliminates the dominant per-iteration cost).

---

## Finding 2: Scalar RNG Calls Inside Nested Loops (HIGH)

**Location:** Lines 169-177 (seed initialization), 190-206 (classroom new infections), 213-220 (household new infections).

**Problem:** Every newly infected individual triggers two scalar `rpois(1, ...)` calls. R's RNG functions have per-call overhead (argument checking, .Internal dispatch). When many individuals are infected over a run, these scalar calls accumulate.

**Fix:** Vectorize the RNG. For seeds:

```r
seed_idx <- sample(1:n, seeds)
L_vec <- rpois(seeds, latent_days)
I_vec <- rpois(seeds, infectious_days)
v_susp[seed_idx] <- 0L
v_exposed[seed_idx] <- 1
v_start[seed_idx] <- 1 + L_vec
v_end[seed_idx] <- 1 + L_vec + I_vec
v_seed[seed_idx] <- 1L
```

For new infections within the day loop, batch all newly infected from a single day and call `rpois(length(all_new), ...)` once per day instead of per individual.

**Estimated speedup: 1.2x-1.5x** (moderate, since infection counts are typically small per day).

---

## Finding 3: Transmission Functions Return Full Data Frame (HIGH)

**Location:** `transmit_classroom()` (lines 98-115), `transmit_household()` (lines 120-137).

**Problem:** Each function takes and returns the entire `df` (16 columns x 944 rows). The return value `list(df = df, new_infected = new_infected)` allocates a list, and the caller immediately destructures it. This creates unnecessary intermediate objects and prevents R from optimizing in place.

**Fix:** Eliminate these functions entirely. Inline the logic into the day loop, operating on extracted vectors. The transmission logic is only ~10 lines each and is called nowhere else.

```r
# Inline classroom transmission:
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
      L <- rpois(length(new_inf), latent_days)
      I <- rpois(length(new_inf), infectious_days)
      v_susp[new_inf] <- 0L
      v_exposed[new_inf] <- day
      v_start[new_inf] <- day + L
      v_end[new_inf] <- day + L + I
    }
  }
}
```

**Estimated speedup: 2x-4x** (combined with Finding 1, since this eliminates the copy-on-modify trigger).

---

## Finding 4: Pre-computed Lookup Lists Use `which()` Scans (MEDIUM)

**Location:** Lines 159-165.

**Problem:** The membership lists are built via `lapply(..., function(c) which(df$class == c & ...))`, which scans the full column for each class/household. This is O(N x num_groups). Better to use `split()`.

**Fix:**

```r
# Class members: only students in classes > 0
student_idx <- which(df$type == "Student")
class_members_list <- split(student_idx, df$class[student_idx])

# Household members: all individuals
hh_members_list <- split(seq_len(n), df$HH_ID)
```

`split()` does a single pass to partition the vector. This is O(N) instead of O(N x G).

**Estimated speedup: negligible per run** (these are only called once per `run_model()`), but cleaner code.

---

## Finding 5: `library(tidyverse)` Loaded on Every Worker (MEDIUM)

**Location:** `model_v1.R` line 5, `test_script.R` line 13 (`.packages = c("tidyverse", "data.table")`).

**Problem:** `library(tidyverse)` loads ~25 packages including ggplot2, readr, stringr, etc. None are used inside `run_model()`. Loading tidyverse on each `foreach` worker costs 1-3 seconds of startup time and ~50MB of memory per worker. With 10,000 iterations (even parallelized), the memory bloat can cause GC pressure.

**Fix:** Remove `tidyverse` from `.packages`. The only tidyverse function used in the foreach body is `mutate()` (line 92 of test_script.R), which can trivially be replaced with data.table `:=`.

```r
# Replace:
.packages = c("tidyverse", "data.table")
# With:
.packages = c("data.table")

# Replace the mutate() call:
# OLD: %>% mutate(seeds_overall = sum(seeds), ...)
# NEW: add data.table := assignments
out_vals[, ':='(
  seeds_overall = sum(seeds),
  min_start_overall = min(min_start),
  max_end_overall = max(max_end)
)]
```

**Estimated speedup: 1-3 seconds per worker initialization**, plus reduced memory pressure.

---

## Finding 6: `foreach .combine` Using Pairwise `rbindlist()` (HIGH)

**Location:** `test_script.R` lines 11-13.

**Problem:**

```r
.combine = function(x, y) rbindlist(list(x, y))
```

This is a **pairwise** combine. `foreach` calls `.combine` in a tree-reduction pattern, but each call still allocates a new data.table by binding two together. Over 10,000 iterations, this creates O(n log n) intermediate data.tables, with the later ones being very large. The total memory churn is enormous.

**Fix:** Use `.combine = list` to collect results as a list, then call `rbindlist()` once at the end.

```r
out_list <- foreach(
  j = 1:trials,
  .combine = list,
  .multicombine = TRUE,
  .packages = c("data.table")
) %dopar% {
  # ... same body ...
  out_vals  # return data.table, not rbindlisted
}

out <- rbindlist(out_list)
```

**Estimated speedup: 2x-5x on the combine step** (which becomes significant at 10,000 trials). Also dramatically reduces peak memory usage.

---

## Finding 7: Redundant `table(HH_ID)` Calls in `test_script.R` (MEDIUM)

**Location:** `test_script.R` lines 36-39.

**Problem:** `table(HH_ID)` is computed 4 times independently inside the data.table `j` expression:

```r
num_hh_1 = sum(table(HH_ID) == 1),
num_hh_2 = sum(table(HH_ID) == 2),
num_hh_3 = sum(table(HH_ID) == 3),
num_hh_4 = sum(table(HH_ID) == 4),
```

Each `table()` call scans the entire column and builds a frequency table.

**Fix:** Compute once, reuse:

```r
# Before the data.table expression:
hh_tab <- table(mod$HH_ID)
# Then use hh_tab in the expression, or compute inside a local scope:
out_vals <- {
  hh_tab <- table(mod$HH_ID)
  mod[, .(
    ...
    num_hh_1 = sum(hh_tab == 1),
    num_hh_2 = sum(hh_tab == 2),
    num_hh_3 = sum(hh_tab == 3),
    num_hh_4 = sum(hh_tab == 4),
    ...
  )]
}
```

**Estimated speedup: small per trial** (table on 944 rows is fast), but it runs 10,000 times.

---

## Finding 8: Repeated `susp == 0` Subsetting in `test_script.R` (LOW-MEDIUM)

**Location:** `test_script.R` lines 46-83.

**Problem:** The expression `susp == 0` is evaluated ~15 separate times within the data.table `j` expression. Each evaluation creates a logical vector of length 944 and subsets against it.

**Fix:** Pre-compute the logical mask or the infected indices:

```r
out_vals <- {
  inf_mask <- mod$susp == 0L
  inf <- mod[inf_mask]
  mod[, .(
    infections = sum(inf_mask),
    ...
    inf_days = mean(inf$end - inf$start),
    inf_days_sqr = mean((inf$end - inf$start)^2),
    ...
  )]
}
```

**Estimated speedup: ~1.1x** on the summary step.

---

## Finding 9: `data.table::set()` or Matrix-Based Simulation (ADVANCED)

**Problem:** The core simulation loop has inherent sequential dependencies: infections on day _d_ produce individuals who become infectious on day _d+L_, who can then infect others. This makes full vectorization across days impossible. However, **within a single day**, all infectious individuals act independently (transmission events don't interact within the same day step).

**Opportunity:** Replace the per-individual `for` loop within each day with a single vectorized pass over all infectious individuals simultaneously:

```r
for (day in 1:t) {
  is_weekday <- ((day - 1) %% 7) < 5
  infectious <- which(v_susp == 0L & !is.na(v_start) &
                      v_start <= day & v_end > day)
  if (length(infectious) == 0L) next

  # Draw rates for all infectious at once
  base_rates <- runif(length(infectious), 0, 0.02)
  hh_rates <- rel_HH * base_rates

  # Process all classroom transmissions at once (weekday students only)
  if (is_weekday) {
    student_inf <- infectious[v_type[infectious] == "Student"]
    # For each infectious student, get classmates, check susceptibility, draw binomials
    # This is harder to fully vectorize due to variable-size groups
    # But we can at least batch the RNG calls
  }

  # ... collect all new_infected, then do one vectorized rpois call
}
```

The per-individual loop within each day is difficult to fully vectorize because group membership is ragged (different class sizes, different household sizes). However, if the number of infectious individuals is small (typically <50), the loop is not the primary bottleneck once copy-on-modify is eliminated.

**For truly maximum performance:** An Rcpp implementation of the inner loop would eliminate all R-level overhead. The core loop logic (look up group members, check susceptibility, draw binomials, update state) is straightforward in C++ and would run 20x-100x faster.

**Estimated speedup from Rcpp: 20x-50x on `run_model()` alone.**

---

## Finding 10: Population is Re-used But Columns are Reset Each Run (LOW)

**Location:** `run_model()` lines 146-157.

**Problem:** The function receives the full `df` and resets 10 columns. Since `df` is a data.frame passed from the global `all` object, each `df$col <- value` assignment triggers a copy (first one copies the whole data.frame due to copy-on-modify, subsequent ones are free since `df` now has refcount 1). This is a one-time cost per `run_model()` call, not per day, so it is not the primary bottleneck.

**Fix:** Use `data.table::copy()` + `set()`, or simply build the mutable vectors from scratch:

```r
run_model <- function(df, t, rel_HH, seeds) {
  n <- nrow(df)
  # Only copy the static columns we need
  v_type <- df$type
  v_class <- df$class
  v_HH_ID <- df$HH_ID
  # Initialize mutable vectors directly
  v_susp <- rep(1L, n)
  v_exposed <- rep(NA_real_, n)
  v_start <- rep(NA_real_, n)
  v_end <- rep(NA_real_, n)
  # ... etc.
}
```

This avoids the data.frame copy entirely.

**Estimated speedup: small** (~1.05x), but contributes to cleanliness.

---

## Summary of Recommendations (Ordered by Impact)

| # | Finding | Estimated Speedup | Effort |
|---|---------|-------------------|--------|
| 1 | Eliminate copy-on-modify (extract vectors) | 3x-8x | Low |
| 3 | Inline transmission functions | 2x-4x (combined with #1) | Low |
| 6 | Fix foreach `.combine` to use list + single `rbindlist()` | 2x-5x on combine | Low |
| 2 | Vectorize RNG calls | 1.2x-1.5x | Low |
| 5 | Remove `tidyverse` from worker packages | 1-3s startup + memory | Trivial |
| 7 | Compute `table(HH_ID)` once | ~1.05x per trial | Trivial |
| 8 | Pre-compute `susp == 0` mask | ~1.1x on summary | Trivial |
| 4 | Use `split()` for membership lists | Negligible | Trivial |
| 10 | Initialize vectors directly | ~1.05x | Low |
| 9 | Rcpp for inner loop | 20x-50x on run_model | High |

**Combined estimate without Rcpp: 10x-20x.**
**Combined estimate with Rcpp for inner loop: 30x-100x.**

---

## Recommended Implementation Order

1. **Phase 1 (Low-hanging fruit, ~1 hour):** Findings 1+3 (extract vectors, inline functions), Finding 2 (vectorize RNG), Finding 4 (split), Finding 10 (direct init). This produces a `model_v2.R` with the same interface. Expected: **5x-10x speedup on `run_model()`**.

2. **Phase 2 (test_script.R fixes, ~30 min):** Findings 5, 6, 7, 8. Fix the foreach combine, remove tidyverse, deduplicate computations. Expected: **2x-3x additional speedup on total wall time**.

3. **Phase 3 (Rcpp, ~2-3 hours):** Finding 9. Write the inner day loop in C++. Expected: **additional 10x-20x on `run_model()`**.

---

## Benchmarking Suggestion

After implementing Phase 1, verify with:

```r
library(bench)

df_test <- create_population()
bench::mark(
  v1 = run_model_v1(df_test, t = 30, rel_HH = 2, seeds = 5),
  v2 = run_model_v2(df_test, t = 30, rel_HH = 2, seeds = 5),
  min_iterations = 50,
  check = FALSE  # stochastic, can't check equality
)
```
