library(dplyr)
library(knitr)
library(tibble)
library(survival)
library(readr)
library(splines)
library(ggplot2)
library(survminer)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(condSURV)

df = read.csv('desktop/statistics/memo3/cross-verified-database.csv',encoding="latin1")

df = df %>% filter(un_region != "")
df = df %>% mutate(lifespan = death - birth)
dx = df %>% filter(birth >= 1500)
dx = dx %>% select(birth, lifespan, gender, un_region, level1_main_occ)
dx = dx %>% filter(gender %in% c("Female", "Male"))
dx = dx %>% filter(!(level1_main_occ %in% c("Missing", "Other")))
censor_year = 2022
dx = dx %>% mutate(clifespan = ifelse(is.na(lifespan), censor_year - birth, lifespan))
dx = dx %>% mutate(died = is.finite(dx$lifespan))
dx = dx %>% select(-lifespan)
dx = dx[complete.cases(dx),]
dx = dx %>% mutate(era = floor((birth - 1500) / 100))
dx = dx %>% filter(clifespan <= 100)
dx <- dx %>% mutate(gender = recode(gender, `Male` = 0, `Female` = 1))
dx <- dx %>% mutate(un_region = recode(un_region, `Europe` = 0, `America` = 1, 'Asia' = 2,
                                       'Africa' = 3, 'Oceania' = 4))
dx <- dx %>% mutate(level1_main_occ = recode(level1_main_occ, `Culture` = 0, `Leadership` = 1, 
                                             'Discovery/Science' = 2, 'Sports/Games' = 3))

s1 <- survfit(Surv(clifespan, died) ~ 1, data = dx)
str(s1)

survfit2(Surv(clifespan, died) ~ dx$era, data = dx) %>% 
  ggsurvfit() +
  labs(
    x = "Ages",
    y = "Overall survival probability"
  )
sf = survfit(Surv(clifespan, died) ~ level1_main_occ, dx)
plot(sf, xlab="Age", ylab="Proportion alive", col=1:6, xlim=c(0, 100))
legend(x="topright", c("Culture", "Leadership", "Discovery/Science", "Sports/Games"),
       col=1:4, lty=1)


