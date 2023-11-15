# Introduction
With the urbanization in modern society, trucks, households and factories can emit complex mixtures of ambient air pollutants, many of which can cause severe problems to human health, like respiratory diseases, cancer, and lower lung function. According to one study in Europe[1], 14% of pediatric asthma cases and 15% of all pediatric asthma exacerbations could be attributed to air pollution. Our main goal of the project is to find the relationship between health conditions of residents and outdoor air quality in 101 counties among 44 states of the US from 2018 to 2019, where we used the prevalence of asthma and mortality rate of COPD to assess the vulnerability towards residents in each county because they can be detected quickly when people get exposed to too much air pollutants compared to other diseases like cancer, which have important windows for exposures occurring decades before diagnosis. We first established the unadjusted models with respect to each pollutant. But many other risk factors (like tobacco smoke exposure, genomic influences and stress) also play a role in explaining the differences of prevalence, so we add some other covariates into the adjusted models. Then we’ll do some analysis and discussions based on the results obtained from the MCMC output using R2jags.
# Model
We conduct an ecological study in the United States across the years of 2018 and 2019, and mainly collected the data from two sources. The first one is outdoor air quality data from United States Environmental Protection Agency (EPA)[2], and the other one is national environmental public health data from Centers for Disease Control and Prevention (CDC)[3]. The variables we choose for the air quality are daily mean PM 2.5 concentration ($\mathrm{\mu g/m^3}$), daily max 1-hour $\mathrm{NO_2}$ and $\mathrm{SO_2}$ concentration (ppb) and daily max 8-hour $\mathrm{CO}$ concentration (ppm), then compute the mean value of these four covariates in by year. We separate our data by year. For each year and each air pollutant, we fit an unadjusted model, in which we make a linear regression of age-adjusted COPD mortality rate or age-adjusted asthma prevalence among adults on the air pollutant. We then make an adjusted model where we add several covariates to the model. We adjust for percent of females in the county, the percent of nonwhite individuals in the county, the percent of individuals living in poverty, the percent of obese individuals in the county, and the percent of smokers in the county. These variables were selected based on a literature review[1,4-6].

For the unadjusted model, we fit the following model for both age-adjusted COPD mortality and age-adjusted prevalence of asthma among adults:

$$
[Y_{ij}|\mu_{ij},\sigma_j^2]\sim\mathcal{N}(\mu_{ij},\sigma_j^2)
$$

$$
\mu_{ij}=\beta_{0j}+\beta_{1j}\cdot \text{pollutant}_{ij}
$$

$$
\pi(\boldsymbol{\beta}_j)\propto 1
$$

$$
\pi(\sigma_j^{-2})\propto \sigma_j^2
$$

where $j\in\set{1,2}$ represents the year 2018 ($j=1$) or 2019 ($j=2$) and $i$ is the county index. Pollutant is replaced by one of the following pollutants: $\mathrm{CO}$, $\mathrm{NO_2}$, $\mathrm{SO_2}$, or PM 2.5. This means that a separate model was used for each pollutant.

For the adjusted model, we consider pollutants, percentage of female, nonwhite, obesity, poverty and smoking when computing $\mu_{ij}$. For each pollutant, we calculate the posterior slope and its 95% credible interval (CrI). We additionally calculate the posterior probability of the slope being in the specific direction found from the estimate for the adjusted results. For example, if we find that the posterior mean of $\beta_{11} | Y,X$ was positive for $\mathrm{CO}$, then we calculate $\text{Pr}(\beta_{11}>0|Y,X)$. We considered a finding significant if this posterior probability was at least 90%.

# Data Analysis
The Table 1 in the Appendix shows the posteriors of all the parameters in our regression model. We used 5 chains, 100K iterations and 50K burn-in for COPD. Asthma models used 1M iterations, 500K burn-in, and 5 chains.

We will highlight the results we found that had at least a 90% posterior probability in the direction of the posterior mean. In 2019, we found that given this data, a 1 unit increase in $\mathrm{CO}$ concentration was associated with a 12.164 unit higher mean age-adjusted COPD mortality rate (95% CrI: -2.536, 26.906; posterior probability 0.9478). In both 2018, we found that $\mathrm{NO_2}$ was associated with age-adjusted COPD mortality rate. Given this data, the posterior mean estimate of the slope associated with a 1 unit increase in $\mathrm{NO_2}$ was -0.379 (95% CrI: -0.718, -0.038; posterior probability: 0.9850). In 2019, we found that SO2 was associated with age-adjusted COPD mortality rate, with a posterior mean slope of -1.407 (95% CrI: -3.324, 0.502; posterior probability 0.9262). Finally, we observed that in 2019 that PM 2.5 was associated with age-adjusted asthma prevalence among adults, with a posterior mean slope of -0.111 (95% CrI: -0.211, -0.012; posterior probability 0.9856).

According to Table 1, the mean values of posterior parameters of $\mathrm{NO_2}$, $\mathrm{SO_2}$ and PM 2.5 are all close to 0 compared to $\mathrm{CO}$ in the model whose response variable is the prevalence of asthma, whether it’s adjusted or not, so these three air pollutants don’t have significant effect to asthma based on our model. If we also take COPD into account, we find that the absolute value of posterior mean of $\mathrm{CO}$ is greater than the other three pollutants in both models, implying it has the most significant effect to the prevalence of asthma and mortality of COPD. Among all the confounders in the adjusted models we chose, we know that the posterior parameter of smoking rate is the greatest in age-adjusted COPD mortality model, and the percent of female is the most significant in age-adjusted prevalence of asthma model. From the posterior means of the adjusted models in Table 1, all the parameters relevant to smoking are bolded, implying that this factor has a great chance to deteriorate human health in terms of asthma and COPD.

# Discussion
We found that $\mathrm{CO}$ was positively associated with age-adjusted COPD mortality in 2019, $\mathrm{NO_2}$ was negatively associated with age-adjusted COPD mortality in 2018, that $\mathrm{SO_2}$ was negatively associated with age-adjusted COPD mortality in 2019, and that PM 2.5 was inversely associated with age-adjusted asthma prevalence among adults in 2019. To our knowledge, we are the first to report inverse associations between these pollutants and COPD mortality rates and asthma prevalence. There are many explanations for why we observe these associations.

One explanation is that given the 95% credible interval for these pollutants capture zero, it is probable that the true association for these pollutants is positive and that we did not have enough data to find this. However, given the posterior probability calculated for these associations are at least 90% and that we utilized flat priors, this explanation is unlikely. Another explanation is we are observing the ecological fallacy. A higher concentration of pollutants at the county level may probably result in behavior change at the individual level to avoid activities outside that may aggravate the condition of the subjects, even though the pollutants may truly be harmful at an individual level. Another explanation is the need to consider multiple pollutants within the same model. Given that these pollutants constitute a mixture in the atmosphere, an increase in one pollutant may represent a decrease in a more harmful pollutant at the county level. So the inverse association reflects the minor harm of one pollutant being outweighed by removing a more harmful pollutant. A solution to this would be to examine this data using a model with multiple pollutants at the same time.

According to the research[6], PM 2.5 has the greatest effect on respiratory health, but there’s no implication from our output. Most studies consider that PM 2.5 at or below 12 $\mathrm{\mu g/m^3}$ is considered healthy with little to no risk from exposure. In our data set, just 6 out of 202 observations of PM 2.5 are slightly over 12 $\mathrm{\mu g/m^3}$, indicating that almost no counties are affected by PM 2.5. Thus this covariate is outweighed by other risk factors. This may be attributed to the better urban planning and the use of renewable energy. We may address this problem by checking some developing countries, where air pollution is severe because of rapid urbanization.


From the data analysis part, we know that the effect of $\mathrm{NO_2}$, $\mathrm{SO_2}$ and PM 2.5 are mainly outweighed by smoking rate in COPD model. In a study about COPD[5], we know that active smoking remains the main factor of COPD over the decades, which is consistent with our results. The effects of the above three air pollutants are mainly outweighed by sex in asthma model. Though we don’t know if sex is the most significant factor, we do find that the asthma rate will increase if the percent of female grows. This result is reasonable because females were 62% more likely to develop asthma than males according to American Lung Association.

# References
[1] Stern, J., Pier, J., Litonjua, A. A. (2020, February). Asthma epidemiology and risk factors. In Seminars in immunopathology (Vol. 42, No. 1, pp. 5-15). Springer Berlin Heidelberg.

[2] Outdoor air quality data from United States Environmental Protection Agency (EPA)[online]: https://www.epa.gov/outdoor-air-quality-data.

[3] National environmental public health data from CDC[online]: https://ephtracking.cdc.gov/DataExplorer/?query= 51ED8370-BE00-4813-A4F8-AE641EF61672&fips=26161&G5=9999.

[4] Subbarao, P., Mandhane, P. J., Sears, M. R. (2009). Asthma: epidemiology, etiology and risk factors. Cmaj, 181(9), E181-E190.

[5] Raherison, C., Girodet, P. O. (2009). Epidemiology of COPD. European respiratory review, 18(114), 213-221.

[6] Ruvuna, L., Sood, A. (2020). Epidemiology of chronic obstructive pulmonary disease. Clinics in Chest Medicine, 41(3), 315-327.
