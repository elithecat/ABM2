# Code Review Report: model_v1.R

## Summary

This is an agent-based model simulating disease transmission through school classrooms and households. The model creates a population of students (assigned to households and classes) and household adults, seeds initial infections, then simulates day-by-day transmission with latent and infectious periods drawn from Poisson distributions. The test script (`test_script.R`) validates that empirical output statistics converge to theoretical parameter targets within a 1.5% relative tolerance over 10,000 trials.

The model has several significant bugs affecting correctness of the simulation and its validation, plus a number of lesser issues.

---

## Critical Issues

### 1. Transmission probability is drawn from Uniform(0, 0.02) instead of using `transmission_prob_school`

- **Location**: `model_v1.R`, line 191
- **Description**: On each infectious day, for each infectious individual, a new `base_rate` is drawn from `runif(1, 0, 0.02)`. The parameter `transmission_prob_school` (set to 0.01 on line 15) is never used inside the simulation. The mean of `Uniform(0, 0.02)` is 0.01, which happens to equal the parameter, so this will converge in expectation -- but it introduces significant extra variance that is not documented or intentional. The household rate `hh_rate = rel_HH * base_rate` suffers the same issue: each infectious individual gets a different per-day transmission rate, creating individual-level overdispersion that may or may not be desired.
- **Impact**: The model does not implement the stated transmission probability. The test script checks that `attack_rate_student` converges to `transmission_prob_school`, which it will in the mean, but the variance structure of the epidemic is wrong. If anyone changes `transmission_prob_school` without also updating the hardcoded `0.02` on line 191, the model silently breaks.
- **Suggestion**: Replace `runif(1, 0, 0.02)` with the parameter `transmission_prob_school` (passed into `run_model` or referenced from the global scope). If overdispersion is intentionally desired, document it and parameterize the bounds.

### 2. Infectious period of zero is possible, producing individuals who are exposed but never infectious

- **Location**: `model_v1.R`, lines 171, 201, 215
- **Description**: `rpois(1, infectious_days)` with `infectious_days = 5` can return 0. When `I = 0`, `end = start`, so the individual is never infectious (`start <= day & end > day` is never true). Similarly, `rpois(1, latent_days)` can return 0, meaning the individual becomes infectious on the same day they are exposed.
- **Impact**: A zero infectious period means the individual is counted as infected (`susp == 0`) but contributes zero infectious days and zero transmission opportunities. This is epidemiologically questionable and will bias contact tracking and secondary infection statistics downward.
- **Suggestion**: Use `max(1, rpois(1, infectious_days))` or a shifted Poisson `1 + rpois(1, infectious_days - 1)` to ensure at least one infectious day, depending on the epidemiological intent.

### 3. Class assignment silently drops students when `num_students` is not divisible by `num_classes`

- **Location**: `model_v1.R`, line 47
- **Description**: `sample(rep(1:num_classes, each = num_students %/% num_classes))` generates `num_classes * (num_students %/% num_classes)` elements. With `num_students = 500` and `num_classes = 25`, this is `25 * 20 = 500`, which works. But if parameters are changed such that `num_students` is not evenly divisible by `num_classes`, the vector will be shorter than `num_students`, and R's recycling rule will silently recycle values during assignment.
- **Impact**: Incorrect class sizes would invalidate contact-tracking statistics and attack rate calculations. With current parameters (500/25 = 20) this does not trigger, but the code is fragile to parameter changes.
- **Suggestion**: Use a method that handles remainders, e.g., `sample(rep(1:num_classes, length.out = num_students))`.

### 4. Seed infections can select the same individual twice

- **Location**: `model_v1.R`, line 168
- **Description**: `sample(1:n, seeds)` draws without replacement by default in R when `seeds < n`, so this is actually fine for the default case. However, seeds can be adults who will never transmit in classrooms. If an adult is seeded, they can only transmit within households, reducing the expected school transmission.
- **Impact**: Seeds landing on household adults means those seeds can only transmit within households, not classrooms. This may be intentional but is not documented.
- **Suggestion**: If seeds should be students only, restrict to `sample(which(df$type == "Student"), seeds)`. If seeding adults is intended, document it.

---

## Major Concerns

### 5. `track_inf_days` and `track_latent_days` censoring logic is fragile

- **Location**: `model_v1.R`, lines 227-230; `test_script.R`, lines 56-67
- **Description**: The model computes `track_inf_days` as the number of infectious days that fall within `[1, t]` (censored). The test script then "un-censors" them by adding back the adjustment. This correction logic is convoluted and potentially error-prone. While numerically correct in the cases traced, the split between model-side censoring and test-side un-censoring is fragile and error-prone.
- **Suggestion**: Either store the raw (uncensored) `end - start` and `start - exposed` directly in the data frame, or do all censoring/un-censoring in one place.

### 6. Stale susceptibility check during within-day transmission

- **Location**: `model_v1.R`, lines 196-220
- **Description**: When multiple infectious individuals transmit within the same day, the `df` dataframe is updated in-place after each transmission event. The `infectious` vector is computed once at the start of the day. Newly exposed individuals whose latent period is 0 will NOT be in the `infectious` vector for the current day. Same-day transmission chains are impossible.
- **Impact**: The generation interval has a minimum of 1 day even when the latent period is 0. This may slightly underestimate transmission speed but is a common modeling simplification.

---

## Minor Suggestions

### 7. Global variable dependency makes `run_model` impure

- **Location**: `model_v1.R`, lines 170, 201, 215
- **Description**: `run_model` references global variables `latent_days` and `infectious_days` directly instead of taking them as parameters.
- **Suggestion**: Add `latent_days` and `infectious_days` as parameters to `run_model`.

### 8. The variable name `t` shadows R's built-in transpose function

- **Location**: `model_v1.R`, line 142
- **Description**: Using `t` as a parameter name shadows the base R function `t()`.
- **Suggestion**: Rename to `n_days` or `sim_days`.

### 9. The variable name `all` shadows R's built-in function

- **Location**: `model_v1.R`, line 239
- **Description**: `all` is a base R function. Assigning to it shadows the built-in.
- **Suggestion**: Use `pop` or `population` instead.

### 10. `create_population` relies entirely on global parameters

- **Location**: `model_v1.R`, line 23
- **Description**: The function takes no arguments but depends on `num_HH`, `num_students`, `num_classes`, `num_parents_perHH`, and `model_run_length` from the global environment.
- **Suggestion**: Pass these as function arguments for reusability and testability.

### 11. Potential infinite loop in `create_population`

- **Location**: `model_v1.R`, lines 25-43
- **Description**: The `repeat` loop with the inner `while` loop adjusts household sizes. With reasonable parameters this converges quickly, but there is no safeguard (e.g., max iteration count).
- **Suggestion**: Add a maximum iteration guard to prevent hanging.

---

## Positive Observations

- The pre-computation of class and household membership lists (lines 159-165) is a good optimization that avoids repeated subsetting inside the inner loop.
- The weekday/weekend distinction (line 181) is a thoughtful modeling detail.
- The separation of classroom and household transmission into distinct subfunctions improves readability.
- The test script's approach of validating the mean and variance of Poisson-distributed quantities against theoretical targets is a sound verification strategy.
- The censoring logic, while complex, is mathematically correct for computing observable infectious days within the simulation window.

---

## Recommendations (Prioritized)

1. **Fix the transmission rate** (Critical #1): Replace `runif(1, 0, 0.02)` with the actual parameter or parameterize the bounds.
2. **Prevent zero-duration infectious periods** (Critical #2): Add `max(1, ...)` guards or document intent.
3. **Make class assignment robust** (Critical #3): Handle non-divisible student/class counts.
4. **Clarify seed population** (Critical #4): Document or restrict seed selection.
5. **Simplify tracking variables** (Major #5): Store raw durations directly.
6. **Parameterize `run_model` fully** (Minor #7): Pass all parameters as explicit arguments.
7. **Rename shadowed variables** (Minor #8-9): Rename `t` and `all`.
