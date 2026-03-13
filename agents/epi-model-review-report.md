# Epidemiological Model Review Report

## File: `/workspace/scripts/model_v1.R`

**Specification:** `/workspace/GPT_structured_command.md`
**Validation tests:** `/workspace/test_script.R`

---

## 1. Model Summary

This is a discrete-time, agent-based SIR-like model simulating respiratory disease transmission among students and adults in a school-household setting.

**Population:** 500 students assigned to 25 classrooms of 20, allocated across 222 households (Poisson-distributed sibling counts), each household containing exactly 2 adults. Total population: 944 individuals.

**Compartments:** Susceptible -> Exposed (latent) -> Infectious -> Recovered (implicitly, via end-of-infectiousness). No re-infection.

**Transmission:** Per-contact daily probability drawn from Uniform(0, 0.02) for classroom contacts; household transmission is 2x classroom. Classroom mixing occurs on weekdays only; household mixing occurs all 7 days.

**Natural history:** Latent period ~ Poisson(mean=2 days); infectious period ~ Poisson(mean=5 days).

**Seeding:** 5 randomly chosen individuals from the entire population are seeded as exposed on day 1.

**Duration:** 30-day simulation, designed to be run 1,000+ times.

---

## 2. Findings by Severity

### CRITICAL [RED]

#### C1. Poisson distributions allow zero-length latent and infectious periods

The Poisson distribution with mean 2 has P(X=0) = e^(-2) = 0.135. The Poisson with mean 5 has P(X=0) = e^(-5) = 0.0067. This means:

- Approximately 13.5% of infections will have a zero-day latent period, becoming infectious on the same day they are exposed. For a seed infection exposed on day 1 with L=0, `start = 1 + 0 = 1`, so they are immediately infectious on day 1. This is epidemiologically aggressive but not necessarily wrong.

- More critically, approximately 0.67% of infections will have a zero-day infectious period. When I=0, `end = start`, so the individual is never infectious. This is biologically implausible for most pathogens -- an infected person who never becomes infectious and never transmits. These individuals are removed from the susceptible pool but contribute nothing to transmission.

**Suggested fix:** Use `rpois(1, latent_days) + 1` or `max(1, rpois(1, infectious_days))` to ensure minimum 1-day infectious period. Alternatively, use a shifted Poisson or a different discrete distribution.

#### C2. Same-day processing order bias in household transmission

The model iterates over infectious individuals in a `for` loop within each day. When individual A infects individual B on day `d`, B's susceptibility status is immediately updated. When another infectious individual in the same household is processed later, B is no longer susceptible. Since `which()` returns indices in ascending order, individuals with lower IDs (predominantly students) are systematically processed first.

**Impact:** This creates an order-dependent bias that slightly favors lower-ID infectious individuals getting "credit" for secondary infections and slightly reduces total unique infections. With low per-contact transmission probabilities (mean 0.01), the practical effect is small.

**Suggested fix:** Randomize the order of `infectious` each day: `infectious <- sample(infectious)`. Or collect all transmission events first and resolve conflicts afterward.

---

### WARNING [YELLOW]

#### W1. Poisson distribution for latent and infectious periods: equal mean and variance is biologically unusual

The specification explicitly requires Poisson distributions with "equal mean and variance." While mathematically consistent, for most real pathogens the infectious period has a coefficient of variation much lower than Poisson implies. For Poisson(5), CV = 1/sqrt(5) = 0.45, producing infectious periods commonly ranging from 1 to 9 days.

Real respiratory pathogens typically show less variation in infectious duration. A gamma or log-normal distribution with lower CV would be more realistic. However, since the specification explicitly requests Poisson, the implementation is correct per requirements.

#### W2. Uniform(0, 0.02) transmission drawn per infector per day creates correlated contact outcomes

Each infectious individual draws a single `base_rate ~ U(0, 0.02)` each day, applied to ALL their contacts. This creates correlated transmission across contacts of the same infector on the same day: a "highly infectious day" means high probability for all contacts simultaneously. This produces superspreading-like overdispersion.

The mean transmission probability is E[U(0,0.02)] = 0.01, matching the specification's "average transmission probability of 0.01 per day." The variance in secondary infections will be higher than if rates were drawn independently per contact.

#### W3. Household transmission rate applies uniformly regardless of household size

The household transmission rate is `2 * base_rate` per-contact regardless of household size. Larger households (5-8 members) will have more expected secondary infections per infectious individual. This is frequency-dependent transmission, which tends to overestimate transmission in large households compared to reality where contact time per person decreases with size.

#### W4. No recovery/removal compartment tracking

After `end`, individuals are no longer infectious but remain with `susp == 0`. This implicitly creates a Recovered state. There is no explicit tracking of daily incidence or prevalence curves -- only final-state data frames are output.

#### W5. Population structure is fixed across all runs

The call to `create_population()` creates the structure once. All simulation runs share identical household and classroom assignments. This means the model estimates transmission dynamics conditional on a fixed population structure, not averaging over structure uncertainty.

---

### SUGGESTION [GREEN]

#### S1. Missing output of average transmission rates (Specification step 4)

The specification requires: "Print out the average transmission rates in classroom vs. in households, and the actual number of infections in classroom vs. in households over weekdays vs. weekends." The model code does not print these values.

#### S2. No explicit weekday/weekend breakdown of secondary infections

The specification asks to track secondary infections over weekdays vs. weekends separately. The model tracks totals but does not break them down by day type. Since classroom transmission only occurs on weekdays, `secondary_infs_school` is implicitly weekday-only. But `secondary_infs_home` combines weekday and weekend household infections.

---

## 3. Detailed Analysis

### Transmission Dynamics Verification

Expected per-contact transmission probability:
- Classroom: E[U(0, 0.02)] = 0.01 per susceptible classmate per day
- Household: 2 * E[U(0, 0.02)] = 0.02 per susceptible household member per day

With 20 students per class, an infectious student contacts 19 susceptibles (maximum). Expected classroom secondary infections per day = 19 * 0.01 = 0.19. Over ~5 infectious days and ~4.3 weekdays per week, expected classroom infections per index case ≈ 0.68.

For households with mean ~4.25 members (2.25 students + 2 adults), an infectious person contacts ~3.25 susceptibles. Expected household secondary infections per day = 3.25 * 0.02 = 0.065. Over 5 infectious days ≈ 0.325.

These are low transmission rates, with R0 approximately 1.0 or slightly below. This is plausible for modeling early-stage transmission dynamics of a pathogen with moderate transmissibility.

### Seed Infection Initialization

Seeds are exposed on day 1 with normal Poisson-drawn latent and infectious periods. When L=0 (13.5% probability), the seed is infectious from day 1. Five seeds from 944 individuals is a 0.53% seeding rate, appropriate for modeling an introduction event. All seeds are introduced simultaneously, which is a common and acceptable simplification.

### Population Structure

- 222 households for 500 students with Poisson(2): E[sum] = 444, so the adjustment algorithm adds ~56 students, skewing the distribution rightward.
- 25 classes of exactly 20 students is a clean simplification.
- 2 adults per household (no single-parent, no multigenerational) is acceptable for a first-pass model.

### Conservation Check

No individuals are added or removed. Total population remains 944 throughout. All transitions are discrete state changes. The model is numerically stable with no risk of negative compartments.

---

## 4. Verification Checks

The following tests should be run to validate the model:

1. **Zero transmission test:** Set base_rate = 0. Verify only seed infections appear.
2. **No household transmission:** Set `rel_transmission_HH = 0`. Verify only classroom transmission occurs on weekdays.
3. **All weekends:** Force `is_weekday = FALSE`. Verify no classroom transmission but household continues.
4. **Single seed:** Set `seeds = 1`. Trace secondary infections manually.
5. **Conservation:** Verify `sum(susp == 1) + sum(susp == 0) == 944` at end of every run.
6. **Poisson parameter recovery:** Verify mean and variance of `(end - start)` and `(start - exposed)` match Poisson parameters.
7. **Attack rate recovery:** Verify ratios converge to 0.01 (classroom) and 0.02 (household).
8. **Sensitivity to zero-length periods:** Count fraction of infections with I=0 and assess impact.
9. **Order-dependence test:** Shuffle `infectious` vector and compare distributions.

---

## 5. Overall Assessment

The model is a competent implementation of the specification. The core transmission mechanics are correctly coded: classroom mixing on weekdays, household mixing on all days, binomial draws with appropriate rates, and proper state tracking. The population structure generation is sound.

The most consequential issues are:

1. **Zero-length infectious periods** (C1) -- a small but systematic bias. For Poisson(5), only 0.67% of cases are affected, so practical impact is minimal, but it represents a biological implausibility.

2. **Processing order bias** (C2) -- within-day ordering creates a subtle deterministic bias. Randomizing processing order would eliminate this.

3. **Correlated transmission draws** (W2) -- one rate per infector per day creates overdispersion. Whether this is a feature or bug depends on modeling intent.

The model is appropriate for its stated purpose of simulating early-phase school-household transmission with low transmission probabilities. The 30-day simulation window is sufficient to observe initial wave dynamics given the ~7-day serial interval (2-day latent + 5-day infectious).
