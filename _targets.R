# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)

# Set target options:
tar_option_set(packages = c("tidySEM", "kableExtra", "bain", "lavaan", "metafor"))

list(
  tar_target(
    name = res_rma,
    command = rma(yi = kuiper2013$beta,
                  vi = kuiper2013$vi)
  )
  , tar_target(
    name = res_pbf_t1,
    command = {
      set.seed(1)
      pbf(yi = kuiper2013$beta, vi = kuiper2013$vi, ni = kuiper2013$n,
          hypothesis = "y = 0; y > 0; y < 0")}
  )
  
  , tar_target(
    name = res_t2,
    command = {
      
      # Load lavaan package for SEM
      library(lavaan)
      
      # Specify SEM-model for latent variable regression
      model_nl <- "
fam =~ fam_1 + fam_2 + fam_3
con =~ sepa_soc_1 + sepa_soc_2 + sepa_soc_3 + sepa_soc_4 + sepa_soc_5 +
       sepa_eco_1 + sepa_eco_2 + sepa_eco_3 + sepa_eco_4 + sepa_eco_5
con ~ beta * fam"
      
      # Estimate the model in lavaan
      results_nl <- sem(model = model_nl, data = synthetic_nl)
      # Test that the effect labeled 'beta' is positive
      set.seed(1)
      bain(results_nl,
           hypothesis = "beta > .1",
           standardize = TRUE)
      
      
    }
  )
  
  , tar_target(
    name = res_pbf_t3,
    command = {
      model_nl <- "
fam =~ fam_1 + fam_2 + fam_3
con =~ sepa_soc_1 + sepa_soc_2 + sepa_soc_3 + sepa_soc_4 + sepa_soc_5 +
       sepa_eco_1 + sepa_eco_2 + sepa_eco_3 + sepa_eco_4 + sepa_eco_5
con ~ beta * fam"
      
      # Estimate the model in lavaan
      results_nl <- sem(model = model_nl, data = synthetic_nl)
      # Specify the models for synthetic_dk and synthetic_us
      model_dk <- "
fam =~ fam_1 + fam_2 + fam_3
con =~
sepa_soc_1 + sepa_soc_2 + sepa_soc_3 + sepa_soc_4 + sepa_soc_5 +
sepa_eco_1 + sepa_eco_2 + sepa_eco_3 + sepa_eco_4 + sepa_eco_5
con ~ beta * fam"
      model_us <- "
fam =~ fam_1 + fam_2 + fam_3
con =~
secs_soc_1 + secs_soc_2 + secs_soc_3 + secs_soc_4 + secs_soc_5 +
secs_soc_6 + secs_soc_7 +
secs_eco_1 + secs_eco_2 + secs_eco_3 + secs_eco_4 + secs_eco_5
con ~ beta * fam"
      
      # Estimate the model in lavaan
      results_dk <- sem(model = model_dk, data = synthetic_dk)
      results_us <- sem(model = model_us, data = synthetic_us)
      
      # Bind the models into a list
      results <- list(results_nl, results_dk, results_us)
      # Test the hypothesis that the effect size labeled 'beta' is positive
      set.seed(1)
      pbf(results, hypothesis = "beta > .1", standardize = TRUE)
    }
  )
  , tar_target(
    name = res_pbf4,
    command = {
      
      library(lavaan)
      model_nl <- "
fam =~ fam_1 + fam_2 + fam_3
con =~ sepa_soc_1 + sepa_soc_2 + sepa_soc_3 + sepa_soc_4 + sepa_soc_5 +
       sepa_eco_1 + sepa_eco_2 + sepa_eco_3 + sepa_eco_4 + sepa_eco_5
con ~ beta * fam"
      # Add the additional predictor to the model, label the effect beta2
      # Add the additional predictor to the model, label the effect beta2
      model_nl2 <- c(model_nl, "group =~ grp_1 + grp_2 + grp_3
                         con ~ beta2 * group")
      
      # Estimate the model in lavaan
      results_nl2 <- sem(model = model_nl2, data = synthetic_nl)
      
      model_dk <- "fam =~ fam_1 + fam_2 + fam_3
con =~
sepa_soc_1 + sepa_soc_2 + sepa_soc_3 + sepa_soc_4 + sepa_soc_5 +
sepa_eco_1 + sepa_eco_2 + sepa_eco_3 + sepa_eco_4 + sepa_eco_5
con ~ beta * fam"
      model_us <- "fam =~ fam_1 + fam_2 + fam_3
con =~
secs_soc_1 + secs_soc_2 + secs_soc_3 + secs_soc_4 + secs_soc_5 +
secs_soc_6 + secs_soc_7 +
secs_eco_1 + secs_eco_2 + secs_eco_3 + secs_eco_4 + secs_eco_5
con ~ beta * fam"
      
      # Estimate the model in lavaan
      results_dk <- sem(model = model_dk, data = synthetic_dk)
      results_us <- sem(model = model_us, data = synthetic_us)
      
      set.seed(1)
      bf_nl2 <- bain(results_nl2,
                     hypothesis = "beta > .1; beta2 < .1", 
                     standardize = TRUE)
      bf_dk <- bain(results_dk, hypothesis = "beta > .1", standardize = TRUE)
      bf_us <- bain(results_us, hypothesis = "beta > .1", standardize = TRUE)
      
      # Bind bain objects into a list
      bfs <- list(bf_nl2, bf_dk, bf_us)
      
      # Call pbf on that list
      pbf(bfs)
  
    }
)

, tar_target(name = res_pbf_t5,
             command = {
               # Create mean scale scores
               synthetic_nl <- data.frame(
                 family = rowMeans(synthetic_nl[c("fam_1", "fam_2", "fam_3")]),
                 conservative = rowMeans(synthetic_nl[c("sepa_soc_1", "sepa_soc_2", "sepa_soc_3",
                                                        "sepa_soc_4", "sepa_soc_5", "sepa_eco_1",
                                                        "sepa_eco_2", "sepa_eco_3", "sepa_eco_4",
                                                        "sepa_eco_5")]))
               synthetic_dk <- data.frame(
                 family = rowMeans(synthetic_dk[c("fam_1", "fam_2", "fam_3")]),
                 conservative = rowMeans(synthetic_dk[c("sepa_soc_1", "sepa_soc_2", "sepa_soc_3",
                                                        "sepa_soc_4", "sepa_soc_5", "sepa_eco_1",
                                                        "sepa_eco_2", "sepa_eco_3", "sepa_eco_4",
                                                        "sepa_eco_5")]))
               
               synthetic_us <- data.frame(
                 family = rowMeans(synthetic_us[c("fam_1", "fam_2", "fam_3")]),
                 conservative = rowMeans(synthetic_us[c("secs_soc_1", "secs_soc_2", "secs_soc_3",
                                                        "secs_soc_4", "secs_soc_5", "secs_soc_6",
                                                        "secs_soc_7", "secs_eco_1", "secs_eco_2",
                                                        "secs_eco_3", "secs_eco_4", "secs_eco_5")]))
               
               # synthetic_nl: Conduct t-test using Cohen's D
               synthetic_nl$group <- cut(synthetic_nl$conservative, breaks = 2,
                                         labels = c("liberal", "conservative"))
               sample_sizes <- table(synthetic_nl$group)
               sds <- tapply(synthetic_nl$family, synthetic_nl$group, sd)
               pooled_sd <- sqrt(sum((sample_sizes - 1) * sds) / (sum(sample_sizes) - 2))
               NL_est <- diff(tapply(synthetic_nl$family, synthetic_nl$group, mean)) / pooled_sd
               NL_var <- (sum(sample_sizes) / prod(sample_sizes)) +
                 (NL_est^2 / (2*sum(sample_sizes)))
               
               # synthetic_dk: Conduct bivariate regression
               DK_fit <- lm(conservative ~ family, data = synthetic_dk)
               DK_est <- coef(DK_fit)["family"]
               DK_var <- vcov(DK_fit)["family", "family"]
               
               # synthetic_us: Correlation coefficient
               US_est <- cor(synthetic_us)[1, 2]
               US_var <- (1 - US_est^2)^2 / (nrow(synthetic_us) - 1)
               
               # Name the estimates so hypotheses will be the same
               names(NL_est) <- names(DK_est) <- names(US_est) <- "parameter"
               
               # Use bain.default() to obtain BF for the central hypothesis 
               NL_bain <- bain(x = NL_est, 
                               Sigma = matrix(NL_var, 1, 1),
                               n = nrow(synthetic_nl),
                               hypothesis = "parameter > 0",
                               joint_parameters = 1)
               DK_bain <- bain(x = DK_est,
                               Sigma = matrix(DK_var, 1, 1),
                               n = nrow(synthetic_dk),
                               hypothesis = "parameter > 0",
                               joint_parameters = 1)
               US_bain <- bain(x = US_est,
                               Sigma = matrix(US_var, 1, 1),
                               n = nrow(synthetic_us),
                               hypothesis = "parameter > 0",
                               joint_parameters = 1)
               
               # Aggregate evidence using pbf()
               pbf(list(US_bain, DK_bain, NL_bain))
             })

, tarchetypes::tar_render(manuscript, "tutorial.Rmd", cue = tar_cue(mode = "always"))
)
