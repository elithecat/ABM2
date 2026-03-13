library(foreach)
library(tictoc)
# number of trials to run
trials = 1e4 ## changed for checking if any errors

# repeat model runs
tic("Model Starts")
doMC::registerDoMC(cores = parallel::detectCores())
out = foreach(
  j = 1:trials,
  .combine = function(x, y)
    rbindlist(list(x, y)),
  .packages = c("tidyverse", "data.table")
) %dopar% {
  # run model
  mod = data.table(
    run_model(
      all,
      t = model_run_length,
      rel_HH = rel_transmission_HH,
      seeds = num_seeds_inf
    )
  )
  
  # take sum with data.table
  # but this is easy to rewrite with base R or tidyverse
  
  
  
  out_vals = mod[, .(
    infections = sum(susp != 1),
    
    #Household characteristics
    num_indv = .N,
    num_hh = length(unique(HH_ID)),
    num_hh_1 = sum(table(HH_ID) == 1),
    num_hh_2 = sum(table(HH_ID) == 2),
    num_hh_3 = sum(table(HH_ID) == 3),
    num_hh_4 = sum(table(HH_ID) == 4),
    
    #Class characteristics
    num_students = sum(.N[type == "Student"]),
    num_class = length(unique(class)),
    
    # test inf days as set
    inf_days = mean(end[susp == 0] - start[susp == 0]),
    inf_days_sqr = mean((end[susp == 0] - start[susp == 0]) ^
                          2),
    
    # test latent days as set
    latent_days = mean(start[susp == 0] - exposed[susp == 0]),
    latent_days_sqr = mean((start[susp == 0] - exposed[susp ==
                                                         0]) ^ 2),
    
    # test inf_days as implemented
    track_inf_days = mean(track_inf_days[susp == 0] + pmax(
      pmin(end[susp == 0] - model_run_length - 1, end[susp == 0] - start[susp == 0]), 0
    )),
    track_inf_days_sqr = mean((track_inf_days[susp == 0] + pmax(
      pmin(end[susp == 0] - model_run_length - 1, end[susp == 0] - start[susp == 0]), 0
    )) ^ 2),
    
    # test latent_days as implemented
    track_latent_days = mean(track_latent_days[susp == 0] + pmax(start[susp == 0] - model_run_length - 1, 0)),
    track_latent_days_sqr = mean((
      track_latent_days[susp == 0] + pmax(start[susp == 0] - model_run_length - 1, 0)
    ) ^ 2),
    
    # track secondary HH infections
    track_home_contacts = sum(track_home_contacts[susp ==
                                                    0]),
    secondary_infs_home = sum(secondary_infs_home[susp ==
                                                    0]),
    
    # track secondary student infections
    track_school_contacts = sum(track_school_contacts[susp ==
                                                        0]),
    secondary_infs_school = sum(secondary_infs_school[susp ==
                                                        0]),
    
    # check the start/end_time
    min_start = min(start[susp == 0]),
    max_end = max(end[susp == 0]),
    
    #Track number of seeds
    seeds = sum(seed),
    
    #Track number of model days
    model_days = mean(model_days)
  ), by = c("type")][, id := j] %>%
    mutate(
      seeds_overall = sum(seeds),
      min_start_overall = min(min_start),
      max_end_overall = max(max_end)
    )
  
}
toc(log = TRUE) 

saveRDS(out, "out.rds")
out <- readRDS("out.rds")

# summarize averages
summary = out[, .(
  mean_infs = mean(infections),
  inf_days = weighted.mean(inf_days, w = infections),
  inf_days_var = weighted.mean(inf_days_sqr, w = infections) - weighted.mean(inf_days, w = infections)^2,
  latent_days = weighted.mean(latent_days, w = infections),
  latent_days_var = weighted.mean(latent_days_sqr, w = infections) - weighted.mean(latent_days, w = infections)^2,
  track_inf_days = weighted.mean(track_inf_days, w = infections),
  track_inf_days_var = weighted.mean(track_inf_days_sqr, w = infections) - weighted.mean(track_inf_days, w = infections)^2,
  track_latent_days = weighted.mean(track_latent_days, w = infections),
  track_latent_days_var = weighted.mean(track_latent_days_sqr, w = infections) - weighted.mean(track_latent_days, w = infections)^2,
  attack_rate_student = sum(secondary_infs_school) / sum(track_school_contacts),
  attack_rate_HH = sum(secondary_infs_home) / sum(track_home_contacts),
  min_start = min(min_start_overall),
  max_end = max(max_end_overall),
  frac_hh_1 = sum(num_hh_1) / sum(num_hh),
  frac_hh_2 = sum(num_hh_2) / sum(num_hh),
  frac_hh_3 = sum(num_hh_3) / sum(num_hh),
  frac_hh_4 = sum(num_hh_4) / sum(num_hh),
  avg_class = sum(num_indv) / sum(num_class),
  num_hh = mean(num_hh),
  num_indv = mean(num_indv),
  num_class = mean(num_class),
  avg_seeds_total = mean(seeds_overall),
  model_days = mean(model_days)
), by = c("type")] %>%
  mutate(rel_HH_obs = attack_rate_HH / attack_rate_student)

#Add parameter targets
target_student <- data.frame(
  type = "Target",
  mean_infs = "N/A (model output)",
  inf_days = infectious_days,
  inf_days_var = infectious_days,
  latent_days = latent_days,
  latent_days_var = latent_days,
  track_inf_days = infectious_days,
  track_inf_days_var = infectious_days,
  track_latent_days = latent_days,
  track_latent_days_var = latent_days,
  attack_rate_student = transmission_prob_school,
  attack_rate_HH = transmission_prob_school * rel_transmission_HH,
  min_start = 1,
  max_end = "N/A (model_output)",
  frac_hh_1 = num_1studentHH / num_HH,
  frac_hh_2 = num_2studentHH / num_HH,
  frac_hh_3 = num_3studentHH / num_HH,
  frac_hh_4 = num_4studentHH / num_HH,
  avg_class = num_students / num_classes,
  num_hh = num_HH,
  num_indv = num_students,
  num_class = num_classes,
  avg_seeds_total = num_seeds_inf,
  model_days = model_run_length,
  rel_HH_obs = rel_transmission_HH
)

target_parents <- data.frame(
  type = "Target",
  mean_infs = "N/A (model output)",
  inf_days = infectious_days,
  inf_days_var = infectious_days,
  latent_days = latent_days,
  latent_days_var = latent_days,
  track_inf_days = infectious_days,
  track_inf_days_var = infectious_days,
  track_latent_days = latent_days,
  track_latent_days_var = latent_days,
  attack_rate_student = "N/A",
  attack_rate_HH = transmission_prob_school * rel_transmission_HH,
  min_start = 1,
  max_end = "N/A (model_output)",
  frac_hh_1 = ifelse(num_parents_perHH == 1, 1, 0),
  frac_hh_2 = ifelse(num_parents_perHH == 2, 1, 0),
  frac_hh_3 = ifelse(num_parents_perHH == 3, 1, 0),
  frac_hh_4 = ifelse(num_parents_perHH == 4, 1, 0),
  avg_class = num_parents_perHH * num_HH,
  num_hh = num_HH,
  num_indv = num_parents_perHH * num_HH,
  num_class = 1,
  avg_seeds_total = num_seeds_inf,
  model_days = model_run_length,
  rel_HH_obs = rel_transmission_HH
)

student_table <-
  rbind(summary[type == "Student"], target_student) %>%
  #Drop model output and lower-level tracking variables
  select(!c(
    mean_infs,
    inf_days,
    inf_days_var,
    latent_days,
    latent_days_var,
    max_end
  )) %>%
  t() %>%
  as.data.frame()
colnames(student_table) <- student_table["type",]
student_table <- student_table %>%
  filter(row.names(.) != "type") %>%
  mutate(
    Student = as.numeric(Student),
    Target = as.numeric(Target),
    relative_difference = abs(Student - Target) / Student
  ) %>%
  mutate(meets.threshold = relative_difference < 0.015)

parents_table <-
  rbind(summary[type == "Household"], target_parents) %>%
  #Drop model output and lower-level tracking variables
  select(
    !c(
      mean_infs,
      inf_days,
      inf_days_var,
      latent_days,
      latent_days_var,
      max_end,
      attack_rate_student,
      rel_HH_obs
    )
  ) %>%
  t() %>%
  as.data.frame()
colnames(parents_table) <- parents_table["type",]
parents_table <- parents_table %>%
  filter(row.names(.) != "type") %>%
  mutate(
    Household = as.numeric(Household),
    Target = as.numeric(Target),
    relative_difference = abs(Household - Target) / Household
  ) %>%
  mutate(meets.threshold = relative_difference < 0.015)

# Want meets.threshold to be TRUE in all cases

## output (added)
print("Model Success!")
log_data <- tic.log(format = FALSE)
log_data[[length(log_data)]]$callback_msg