library(gdata)
library(coxme)
library(multcomp)
library(survival)
library(dplyr)

data <- read.csv("PATH/TO/FileS4_FitnessAndLifespan_MacroVsNoura.csv", header = TRUE)

### ANALYZE THE LIFESPAN DATA ###
lifespan <- subset(data, select = c("Treatment", "Female_spp", "Cross2male.20h.", "ID", "death_day", "censored"))
names(lifespan) <- c("strain.treatment", "female", "male", "id", "death.age", "censored")
lifespan$strain.treatment <- as.factor(lifespan$strain.treatment)
lifespan$dead <- 1

## Median Lifespan
span.summary <- as.data.frame(summarise(group_by(lifespan, strain.treatment), median(death.age, na.rm = TRUE), sd(death.age, na.rm = TRUE), sum(!is.na(death.age))))
names(span.summary) <- c("treatment", "median", "sd", "total_n")
# treatment         median   sd          total_n
# A_macro_con       10       3.491395    34
# B_noura_con        9       2.780648    43
# C_noura_hetero     2       3.876624    88
# D_macro_hetero    10       3.755245    66


## Estimation of treatment effects using planned comparison
model1 <- lm(death.age ~ strain.treatment, data = lifespan)
summary(model1)
# Do a planned comparison of treatment factors: Is there an effect of heterospecific mating on lifespan?
summary(glht(model1, linfct = mcp(strain.treatment = c("A_macro_con - D_macro_hetero = 0",
                                                       "B_noura_con - C_noura_hetero = 0"))))


# Do a planned comparison of treatment factors: Does lifespan differ between species?
summary(glht(model1, linfct = mcp(strain.treatment = c("A_macro_con - B_noura_con = 0"))))



## Cox Proportional Hazards treatment effects using random effects model within each strain
## CAUTION: This model may be over-fitting the data.
## CAUTION: It assumes that each "id" represents a *biological* replicate. However, if they are *technical* replicates, then this model is NOT correct and should not be used.
model2 <- coxme(Surv(death.age, dead) ~ strain.treatment + (1|id), data = lifespan)
summary(model2)
summary(glht(model2, linfct = mcp(strain.treatment = c("A_macro_con - D_macro_hetero = 0",
                                                       "A_macro_con - B_noura_con = 0",
                                                       "A_macro_con - C_noura_hetero = 0",
                                                       "B_noura_con - C_noura_hetero = 0"))))




### ANALYZE THE FECUNDITY DATA ###
offspring <- data[, c(1, 7:15)]
names(offspring) <- c("strain.treatment", "d0", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "fecundity")
offspring$strain.treatment <- as.factor(offspring$strain.treatment)

## Mean Fecundity
fecundity.summary <- as.data.frame(summarise(group_by(offspring, strain.treatment), mean(fecundity, na.rm = TRUE), sd(fecundity, na.rm = TRUE), sum(!is.na(fecundity))))
names(fecundity.summary) <- c("treatment", "mean", "sd", "total_n")
#        treatment      mean        sd total_n
# 1    A_macro_con 135.73529  80.46731      34
# 2    B_noura_con 225.37209 147.28463      43
# 3 C_noura_hetero  64.45977 125.89307      87
# 4 D_macro_hetero 124.62687  89.02721      67


## Estimation of treatment effects using planned comparison
model1 <- lm(fecundity ~ strain.treatment, data = offspring)
summary(model1)
# Do a planned comparison of treatment factors: Is there an effect of heterospecific mating on offspring?
summary(glht(model1, linfct = mcp(strain.treatment = c("A_macro_con - D_macro_hetero = 0",
                                                       "B_noura_con - C_noura_hetero = 0"))))


# Do a planned comparison of treatment factors: Does offspring differ between species?
summary(glht(model1, linfct = mcp(strain.treatment = c("A_macro_con - B_noura_con = 0"))))
