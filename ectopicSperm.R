data <- read.table("PATH/TO/FileS5_SpermInvasion.txt", header = TRUE)

censored.data <- subset(data, censor == 0)
censored.data$ectopic <- rowSums(censored.data[, 12:13])

time <- c(0.5, 1, 1.5, 2, 4, 6)
prop.ectopic <- c(length(subset(censored.data, ectopic != 0 & time == 0.5)$ectopic) / length(subset(censored.data, time == 0.5)$ectopic),
                  length(subset(censored.data, ectopic != 0 & time == 1)$ectopic) / length(subset(censored.data, time == 1)$ectopic),
                  length(subset(censored.data, ectopic != 0 & time == 1.5)$ectopic) / length(subset(censored.data, time == 1.5)$ectopic),
                  length(subset(censored.data, ectopic != 0 & time == 2)$ectopic) / length(subset(censored.data, time == 2)$ectopic),
                  length(subset(censored.data, ectopic != 0 & time == 4)$ectopic) / length(subset(censored.data, time == 4)$ectopic),
                  length(subset(censored.data, ectopic != 0 & time == 6)$ectopic) / length(subset(censored.data, time == 6)$ectopic))


#Is there a significant effect of time on proportion of females with ectopic sperm?
model1 <- lm(prop.ectopic ~ time)
summary(model1)
