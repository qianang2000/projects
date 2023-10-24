#' The purpose of this file is to run preliminary analyses on the project data.
#'
#' Author: Hasan Abu-Amara
#' Date: 12-3-2022

library(R2jags)
library(coda)

# Set seed for replication.
set.seed(2398)

# Load the data.
data_dir = "H:/BIOSTAT 682/Project/Final_Data/"
air_pollutant_county_info_merged = read.csv(paste0(data_dir, "air_pollutant_and_covar_county_data.csv"), stringsAsFactors = FALSE)

# Subset data by year
air_pollutant_county_info_merged_2018 = subset(air_pollutant_county_info_merged, Year == 2018)
air_pollutant_county_info_merged_2019 = subset(air_pollutant_county_info_merged, Year == 2019)

# Make the data.
air_pollutant_county_info_merged_2018_jags = as.list(air_pollutant_county_info_merged_2018)
air_pollutant_county_info_merged_2018_jags[["N"]] = nrow(air_pollutant_county_info_merged_2018)

# Make the model.
copd_so2_model_2018_unadjusted_noninf  = function()
{
  # Likelihood
  for (i in 1:N) {
    COPD_deaths_age_adjusted[i] ~ dnorm(mu[i], tau)

    mu[i] <- beta[1] + So2_measurement[i] * beta[2]
  }

  # Priors
  beta[1] ~ dunif(-1e6, 1e6)
  beta[2] ~ dunif(-1e6, 1e6)
  tau ~ dgamma(1e-6, 1e-6)

  # Parameters of interest
  sigma2 <- 1 / tau
  mean_so2_IQR <- beta[2] * 1.121
}

copd_so2_model_2018_unadjusted_noninf_params = c("beta", "tau", "sigma2", "mean_so2_IQR")

copd_so2_model_2018_unadjusted_noninf_out = jags(data = air_pollutant_county_info_merged_2018_jags,
                                                parameters.to.save = copd_so2_model_2018_unadjusted_noninf_params,
                                                model.file = copd_so2_model_2018_unadjusted_noninf,
                                                n.iter = 100000,
                                                n.burnin = 50000,
                                                n.chains = 5,
                                                n.thin = 1,
                                                jags.seed = 30239)

print(copd_so2_model_2018_unadjusted_noninf_out)

# Adjusted
copd_so2_model_2018_adjusted_noninf  = function()
{
  # Likelihood
  for (i in 1:N) {
    COPD_deaths_age_adjusted[i] ~ dnorm(mu[i], tau)

    mu[i] <- beta[1] + So2_measurement[i] * beta[2] + percent_female[i] * beta[3]
    + percent_nonwhite[i] * beta[4] + percent_obesity[i] * beta[5] +
      percent_poverty[i] * beta[6] + smoking_crude[i] * beta[7]
  }

  # Priors
  beta[1] ~ dunif(-1e6, 1e6)
  beta[2] ~ dunif(-1e6, 1e6)
  beta[3] ~ dunif(-1e6, 1e6)
  beta[4] ~ dunif(-1e6, 1e6)
  beta[5] ~ dunif(-1e6, 1e6)
  beta[6] ~ dunif(-1e6, 1e6)
  beta[7] ~ dunif(-1e6, 1e6)
  tau ~ dgamma(1e-6, 1e-6)

  # Parameters of interest
  sigma2 <- 1 / tau
  mean_so2_IQR <- beta[2] * 1.121
}

copd_so2_model_2018_adjusted_noninf_params = c("beta", "tau", "sigma2", "mean_so2_IQR")

copd_so2_model_2018_adjusted_noninf_out = jags(data = air_pollutant_county_info_merged_2018_jags,
                                              parameters.to.save = copd_so2_model_2018_adjusted_noninf_params,
                                              model.file = copd_so2_model_2018_adjusted_noninf,
                                              n.iter = 100000,
                                              n.burnin = 50000,
                                              n.chains = 5,
                                              n.thin = 1,
                                              jags.seed = 30239)
print(copd_so2_model_2018_adjusted_noninf_out)

attach.jags(copd_so2_model_2018_adjusted_noninf_out)
sum(mean_so2_IQR > 0) / length(mean_so2_IQR)
detach.jags()

summary(air_pollutant_county_info_merged_2018$So2_measurement)
IQR(air_pollutant_county_info_merged_2018$So2_measurement)


##### SO2 - 2019 ####

# Make the data.
air_pollutant_county_info_merged_2019_jags = as.list(air_pollutant_county_info_merged_2019)
air_pollutant_county_info_merged_2019_jags[["N"]] = nrow(air_pollutant_county_info_merged_2019)

# Make the model.
copd_so2_model_2019_unadjusted_noninf  = function()
{
  # Likelihood
  for (i in 1:N) {
    COPD_deaths_age_adjusted[i] ~ dnorm(mu[i], tau)

    mu[i] <- beta[1] + So2_measurement[i] * beta[2]
  }

  # Priors
  beta[1] ~ dunif(-1e6, 1e6)
  beta[2] ~ dunif(-1e6, 1e6)
  tau ~ dgamma(1e-6, 1e-6)

  # Parameters of interest
  sigma2 <- 1 / tau
  mean_so2_IQR <- beta[2] * 1.1880
}

copd_so2_model_2019_unadjusted_noninf_params = c("beta", "tau", "sigma2", "mean_so2_IQR")

copd_so2_model_2019_unadjusted_noninf_out = jags(data = air_pollutant_county_info_merged_2019_jags,
                                                parameters.to.save = copd_so2_model_2019_unadjusted_noninf_params,
                                                model.file = copd_so2_model_2019_unadjusted_noninf,
                                                n.iter = 100000,
                                                n.burnin = 50000,
                                                n.chains = 5,
                                                n.thin = 1,
                                                jags.seed = 30239)

print(copd_so2_model_2019_unadjusted_noninf_out)

# Adjusted
copd_so2_model_2019_adjusted_noninf  = function()
{
  # Likelihood
  for (i in 1:N) {
    COPD_deaths_age_adjusted[i] ~ dnorm(mu[i], tau)

    mu[i] <- beta[1] + So2_measurement[i] * beta[2] + percent_female[i] * beta[3]
    + percent_nonwhite[i] * beta[4] + percent_obesity[i] * beta[5] +
      percent_poverty[i] * beta[6] + smoking_crude[i] * beta[7]
  }

  # Priors
  beta[1] ~ dunif(-1e6, 1e6)
  beta[2] ~ dunif(-1e6, 1e6)
  beta[3] ~ dunif(-1e6, 1e6)
  beta[4] ~ dunif(-1e6, 1e6)
  beta[5] ~ dunif(-1e6, 1e6)
  beta[6] ~ dunif(-1e6, 1e6)
  beta[7] ~ dunif(-1e6, 1e6)
  tau ~ dgamma(1e-6, 1e-6)

  # Parameters of interest
  sigma2 <- 1 / tau
  mean_so2_IQR <- beta[2] * 1.1880
}

copd_so2_model_2019_adjusted_noninf_params = c("beta", "tau", "sigma2", "mean_so2_IQR")

copd_so2_model_2019_adjusted_noninf_out = jags(data = air_pollutant_county_info_merged_2019_jags,
                                              parameters.to.save = copd_so2_model_2019_adjusted_noninf_params,
                                              model.file = copd_so2_model_2019_adjusted_noninf,
                                              n.iter = 100000,
                                              n.burnin = 50000,
                                              n.chains = 5,
                                              n.thin = 1,
                                              jags.seed = 30239)
print(copd_so2_model_2019_adjusted_noninf_out)

attach.jags(copd_so2_model_2019_adjusted_noninf_out)
sum(mean_so2_IQR > 0) / length(mean_so2_IQR)
detach.jags()

summary(air_pollutant_county_info_merged_2019$So2_measurement)
IQR(air_pollutant_county_info_merged_2019$So2_measurement)

# Hm. Let's try Zellner's G prior.
# Use frequentist methods to estimate the variance.

copd_so2_2018_adj_lm = lm(COPD_deaths_age_adjusted ~ So2_measurement + percent_female + percent_nonwhite + percent_obesity + percent_poverty + smoking_crude, air_pollutant_county_info_merged_2018)
summary(copd_so2_2018_adj_lm)
qqnorm(resid(copd_so2_2018_adj_lm))
qqline(resid(copd_so2_2018_adj_lm))

# Get the things needed for the Zellner's g
sigma2_mle_2018 = (summary(copd_so2_2018_adj_lm)$sigma)^2
sigma2_mle_2018_inv = 1 / sigma2_mle_2018
var_needed_vec = c("So2_measurement",
                   "percent_female",
                   "percent_nonwhite",
                   "percent_obesity",
                   "percent_poverty",
                   "smoking_crude")

X = cbind(rep(1, times = nrow(air_pollutant_county_info_merged_2018)), as.matrix(air_pollutant_county_info_merged_2018[, var_needed_vec]))
XtX = crossprod(X)

# Add this to the data.
air_pollutant_county_info_merged_2018_jags_zellner = air_pollutant_county_info_merged_2018_jags
air_pollutant_county_info_merged_2018_jags_zellner[["XtX"]] = XtX
air_pollutant_county_info_merged_2018_jags_zellner[["sigma2_mle_2018_inv"]] = sigma2_mle_2018_inv

copd_so2_model_2018_adjusted_zellner = function()
{
  # Likelihood
  for (i in 1:N) {
    COPD_deaths_age_adjusted[i] ~ dnorm(mu[i], sigma2_mle_2018_inv)

    mu[i] <- beta[1] + So2_measurement[i] * beta[2] + percent_female[i] * beta[3]
    + percent_nonwhite[i] * beta[4] + percent_obesity[i] * beta[5] +
      percent_poverty[i] * beta[6] + smoking_crude[i] * beta[7]
  }

  # Priors
  zeros <- c(0,0,0,0,0,0,0)
  prec_beta <- sigma2_mle_2018_inv * XtX / g
  beta[1:7] ~ dmnorm(zeros, prec_beta)
  g ~ dunif(0, 1e6)

  # Parameters of interest
  mean_so2_IQR <- beta[2] * 1.121
}

copd_so2_model_2018_adjusted_zellner_params = c("beta", "g", "mean_so2_IQR")
copd_so2_model_2018_adjusted_zellner_out = jags(data = air_pollutant_county_info_merged_2018_jags_zellner,
                                               parameters.to.save = copd_so2_model_2018_adjusted_zellner_params,
                                               model.file = copd_so2_model_2018_adjusted_zellner,
                                               n.iter = 100000,
                                               n.burnin = 50000,
                                               n.chains = 5,
                                               n.thin = 1,
                                               jags.seed = 30239)

print(copd_so2_model_2018_adjusted_zellner_out)


###### Asthma Prevalence #####
#### 2018 ####
# Make the model.
asthma_so2_model_2018_unadjusted_noninf  = function()
{
  # Likelihood
  for (i in 1:N) {
    prevalence_asthma_adults_age_adjusted[i] ~ dnorm(mu[i], tau)

    mu[i] <- beta[1] + So2_measurement[i] * beta[2]
  }

  # Priors
  beta[1] ~ dunif(-1e6, 1e6)
  beta[2] ~ dunif(-1e6, 1e6)
  tau ~ dgamma(1e-6, 1e-6)

  # Parameters of interest
  sigma2 <- 1 / tau
}

asthma_so2_model_2018_unadjusted_noninf_params = c("beta", "tau", "sigma2")

asthma_so2_model_2018_unadjusted_noninf_out = jags(data = air_pollutant_county_info_merged_2018_jags,
                                                  parameters.to.save = asthma_so2_model_2018_unadjusted_noninf_params,
                                                  model.file = asthma_so2_model_2018_unadjusted_noninf,
                                                  n.iter = 100000,
                                                  n.burnin = 50000,
                                                  n.chains = 5,
                                                  n.thin = 1,
                                                  jags.seed = 30239)

print(asthma_so2_model_2018_unadjusted_noninf_out)

# 2019 - model
asthma_so2_model_2019_unadjusted_noninf_out = jags(data = air_pollutant_county_info_merged_2019_jags,
                                                   parameters.to.save = asthma_so2_model_2018_unadjusted_noninf_params,
                                                   model.file = asthma_so2_model_2018_unadjusted_noninf,
                                                   n.iter = 100000,
                                                   n.burnin = 50000,
                                                   n.chains = 5,
                                                   n.thin = 1,
                                                   jags.seed = 30239)

print(asthma_so2_model_2019_unadjusted_noninf_out)

# Adjusted
asthma_so2_model_2018_adjusted_noninf  = function()
{
  # Likelihood
  for (i in 1:N) {
    prevalence_asthma_adults_age_adjusted[i] ~ dnorm(mu[i], tau)

    mu[i] <- beta[1] + So2_measurement[i] * beta[2] + percent_female[i] * beta[3]
    + percent_nonwhite[i] * beta[4] + percent_obesity[i] * beta[5] +
      percent_poverty[i] * beta[6] + smoking_crude[i] * beta[7]
  }

  # Priors
  beta[1] ~ dunif(-1e6, 1e6)
  beta[2] ~ dunif(-1e6, 1e6)
  beta[3] ~ dunif(-1e6, 1e6)
  beta[4] ~ dunif(-1e6, 1e6)
  beta[5] ~ dunif(-1e6, 1e6)
  beta[6] ~ dunif(-1e6, 1e6)
  beta[7] ~ dunif(-1e6, 1e6)
  tau ~ dgamma(1e-6, 1e-6)

  # Parameters of interest
  sigma2 <- 1 / tau
}

asthma_so2_model_2018_adjusted_noninf_params = c("beta", "tau", "sigma2")

asthma_so2_model_2018_adjusted_noninf_out = jags(data = air_pollutant_county_info_merged_2018_jags,
                                                parameters.to.save = asthma_so2_model_2018_adjusted_noninf_params,
                                                model.file = asthma_so2_model_2018_adjusted_noninf,
                                                n.iter = 1000000,
                                                n.burnin = 500000,
                                                n.chains = 1,
                                                n.thin = 1,
                                                jags.seed = 30239)
print(asthma_so2_model_2018_adjusted_noninf_out)
geweke.diag(asthma_so2_model_2018_adjusted_noninf_out)


# 2019 - model
asthma_so2_model_2019_adjusted_noninf_out = jags(data = air_pollutant_county_info_merged_2019_jags,
                                                 parameters.to.save = asthma_so2_model_2018_adjusted_noninf_params,
                                                 model.file = asthma_so2_model_2018_adjusted_noninf,
                                                 n.iter = 1000000,
                                                 n.burnin = 500000,
                                                 n.chains = 1,
                                                 n.thin = 1,
                                                 jags.seed = 30239)
print(asthma_so2_model_2019_adjusted_noninf_out)

