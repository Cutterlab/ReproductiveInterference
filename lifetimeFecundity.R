library(dplyr)
library(lme4)

#Read in female data
female <- read.table("PATH/TO/FileS1_Female_Lifetime_Fecundity.txt", sep = "\t", header = TRUE)

#Censor females who did not reproduce
female$total    <- rowSums(female[,3:10], na.rm = TRUE)
female.filtered <- subset(female, total != 0)

#Is there a difference in female fecundity between species?
model1 <- glm(total ~ Species, family = poisson(link = "log"), data = female.filtered)
summary(model1)

model2 <- glm(d1 ~ Species, family = poisson(link = "log"), data = female.filtered)
summary(model2)


fm <- c(subset(female.filtered, Species == "C. macrosperma")$total)
fn <- c(subset(female.filtered, Species == "C. nouraguensis")$total)
t.test(fm, fn)


#Is there a difference in female fecundity on day 1?
fm.d1 <- c(subset(female.filtered, Species == "C. macrosperma")$d1)
fn.d1 <- c(subset(female.filtered, Species == "C. nouraguensis")$d1)
t.test(fm.d1, fn.d1)
 


###
#Read in male data
male <- read.table("PATH/TO/FileS2_Male_Lifetime_Fecundity.txt", header = TRUE, sep = "\t")

male$total.d0 <- rowSums(male[, 5:15], na.rm = TRUE) 
male$total.d1 <- rowSums(male[, 16:26], na.rm = TRUE)
male$total.d2 <- rowSums(male[, 27:37], na.rm = TRUE)
male$total.d3 <- rowSums(male[, 38:48], na.rm = TRUE)
male$total.d4 <- rowSums(male[, 49:59], na.rm = TRUE)
male$total.d5 <- rowSums(male[, 60:70], na.rm = TRUE)
male$total.d6 <- rowSums(male[, 71:81], na.rm = TRUE)
male$total.d7 <- rowSums(male[, 82:92], na.rm = TRUE)
male$total.d8 <- rowSums(male[, 93:103], na.rm = TRUE)
male$total.d9 <- rowSums(male[, 104:114], na.rm = TRUE)

male$total.all <- rowSums(male[, 5:114], na.rm = TRUE)

#censor males who were lost during the experiment
censor <- c();
for (i in 1:length(male$notes)) {
  if (male$notes[i] == "") {
    out <- 0
  }
  else {
    out <- 1
  }
  censor <- rbind(censor, out)
}
male$censor <- censor

#Is there a difference is male fecundity between species?
model3 <- glmer(total.all ~ species + (1 | experiment), family = poisson(link = "log"), data = male)
summary(model3)

model4 <- glmer(total.all ~ species + (1 | experiment), family = poisson(link = "log"), data = subset(male, censor == 0))
summary(model4)

model5 <- glmer(total.d1 ~ species + (1 | experiment), family = poisson(link = "log"), data = male)
summary(model5)


mm <- c(subset(male, species == "C. macrosperma JU1857")$total.all)
mn <- c(subset(male, species == "C. nouraguensis JU1823")$total.all)
t.test(mm, mn)



#Is there a difference in male fecundity on day 0?
mm.d1 <- c(subset(male, species == "C. macrosperma JU1857")$total.d0)
mn.d1 <- c(subset(male, species == "C. nouraguensis JU1823")$total.d0)
t.test(mm.d1, mn.d1)
