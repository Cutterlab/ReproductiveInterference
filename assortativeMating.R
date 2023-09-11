#read in assortative mating data file
data <- read.table("PATH/TO/FileS3_assortMating.txt", header = TRUE, sep = "\t")

data$treatment <- interaction(substr(data$focal.f, start = 1, 5), substr(data$foreign.m, start = 1, stop = 4), sep = ".")
data$freq.wrong <- data$numf.wrong / data$numf

#Random sampling loop that asks:
#Given that species x comprises p percent of the population, how frequently do species x females mate with species y males?

#set-up the population
#N = total number of females and males
N <- 602

#x = frequency of macro
x <- seq(0.1, 0.9, by = 0.05)
#y = frequency of nour
y <- seq(0.9, 0.1, by = -0.05)

#draw random pairs
focal <- c();
#macro is rare
for (j in 1:9) {
  #create female and male vectors
  female <- c(rep("macro", N/2 * x[j]), rep("nour", N/2 * y[j]))
  male   <- c(rep("macro", N/2 * x[j]), rep("nour", N/2 * y[j]))
  
  #draw random mating pairs
  nf <- sample(data$numf, size = 500, replace = TRUE)
  
  pairs <- c();
  for (i in 1:500) {
    out <- c();
    
    out$female <- as.data.frame(sample(female, size = nf[i], replace = FALSE))
    out$male   <- as.data.frame(sample(male, size = nf[i], replace = FALSE))
    
    out <- as.data.frame(out)
    names(out) <- c("female", "male")
    out$combo  <- interaction(out$female, out$male, sep = ".")
    
    z <- table(out$combo)
    
    pairs <- rbind(pairs, data.frame(macro.macro = z[1],
                                     macro.nour  = z[3],
                                     row.names = NULL))
    
  }
  
  focal <- rbind(focal, data.frame(ratio      = x[j],
                                   rareF.comM = pairs$macro.nour,
                                   total      = rowSums(pairs)))
}

#nour is rare
for (j in 9:length(x)) {
  #create female and male vectors
  female <- c(rep("macro", N/2 * x[j]), rep("nour", N/2 * y[j]))
  male   <- c(rep("macro", N/2 * x[j]), rep("nour", N/2 * y[j]))
  
  #draw random mating pairs
  nf <- sample(data$numf, size = 500, replace = TRUE)
  
  pairs <- c();
  for (i in 1:500) {
    out <- c();
    
    out$female <- as.data.frame(sample(female, size = nf[i], replace = FALSE))
    out$male   <- as.data.frame(sample(male, size = nf[i], replace = FALSE))
    
    out <- as.data.frame(out)
    names(out) <- c("female", "male")
    out$combo  <- interaction(out$female, out$male, sep = ".")
    
    z <- table(out$combo)
    
    pairs <- rbind(pairs, data.frame(nour.macro  = z[2],
                                     nour.nour   = z[4],
                                     row.names = NULL))
    
  }
  
  focal <- rbind(focal, data.frame(ratio      = x[j],
                                   rareF.comM = pairs$nour.macro,
                                   total      = rowSums(pairs)))
}

focal$ratio2 <- c(rep(0.1, 500), rep(0.15, 500),rep(0.2, 500), rep(0.25, 500), rep(0.3, 500), rep(0.35, 500), rep(0.4, 500), rep(0.45, 500), rep(0.5, 1000),
                  rep(0.45, 500), rep(0.4, 500), rep(0.35, 500), rep(0.3, 500), rep(0.25, 500), rep(0.2, 500), rep(0.15, 500), rep(0.1, 500))
focal <- focal[order(focal$ratio2, decreasing = FALSE), ]

#calculate 95% confidence interval
z <- seq(1, 9000, by = 1000)

conf <- c();
for (i in 1:length(z)) {
  #calculate mean incorrect and sample size
  k <- mean(focal[z[i]:(z[i] + 999), 2], na.rm = TRUE)
  n <- mean(focal[z[i]:(z[i] + 999), 3], na.rm = TRUE)
  
  #exact binomial test
  test <- binom.test(as.integer(k), as.integer(n), alternative = "two.sided", conf.level = 0.95)
  
  
  conf <- rbind(conf, data.frame(ratio            = focal[z[i], 1],
                                 mean.incorrect = test$estimate[1],
                                 lowerbound     = test$conf.int[1],
                                 upperbound     = test$conf.int[2]))
}

plot(conf$ratio, conf$mean.incorrect, type = "o", las = 1, pch = 19, xlim = c(0.1, 0.5), ylim = c(0, 1), xaxt = "n", bty = "l",
     xlab = "Fraction More Rare Female in Population", ylab = "Fraction of Females with Heterospecific Sperm")
axis(side = 1, at = c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5))

polygon(x = c(0.1, conf$ratio, 0.5, 0.5, rev(conf$ratio), 0.1, 0.1),
        y = c(conf$lowerbound[1], conf$lowerbound, conf$lowerbound[5], conf$upperbound[5], rev(conf$upperbound), conf$upperbound[1], conf$lowerbound[1]),
        col = alpha("gray50", 0.25), border = NA)

points(subset(data, treatment == "macro.nour")$ratio, subset(data, treatment == "macro.nour")$freq.wrong, pch = 17, col = alpha("#2171B5", 0.8), cex = 1.5)
points(subset(data, treatment == "noura.macr")$ratio, subset(data, treatment == "noura.macr")$freq.wrong, pch = 15, col = alpha("#CB181D", 0.8), cex = 1.5)

legend("topright", c("female x male", "macro x nour", "nour x macro"), pch = c(0, 17, 15), col = c("white", "#2171B5", "#CB181D"), bty = "n")


#difference between crosses
#0.1 fraction
t.test(data[7:9, 9], data[16:18, 9])
#t = 7.3909, df = 3.7396, p-value = 0.002317
#mean of x (m x n) mean of y (n x m) 
#0.5023588 0.2318459

#0.25 fraction
t.test(data[4:6, 9], data[13:15, 9])
#t = 0.51303, df = 2.8601, p-value = 0.6449
#mean of x (m x n) mean of y (n x m)
#0.4051923 0.3264007

#0.50 fraction
t.test(data[1:3, 9], data[10:12, 9])
#t = 0.52837, df = 3.1485, p-value = 0.6322
#mean of x (m x n) mean of y (n x m)
#0.1352918 0.1052362 #read in assortative mating data file
data <- read.table("Dropbox/ReproductiveInterferenceSharedFiles/data_files/working_for_R/Exp4_assortMating.txt", header = TRUE, sep = "\t")

data$treatment <- interaction(substr(data$focal.f, start = 1, 5), substr(data$foreign.m, start = 1, stop = 4), sep = ".")
data$freq.wrong <- data$numf.wrong / data$numf

#Random sampling loop that asks:
#Given that species x comprises p percent of the population, how frequently do species x females mate with species y males?

#set-up the population
#N = total number of females and males
N <- 602

#x = frequency of macro
x <- seq(0.1, 0.9, by = 0.05)
#y = frequency of nour
y <- seq(0.9, 0.1, by = -0.05)

#draw random pairs
focal <- c();
#macro is rare
for (j in 1:9) {
  #create female and male vectors
  female <- c(rep("macro", N/2 * x[j]), rep("nour", N/2 * y[j]))
  male   <- c(rep("macro", N/2 * x[j]), rep("nour", N/2 * y[j]))
  
  #draw random mating pairs
  nf <- sample(data$numf, size = 500, replace = TRUE)
  
  pairs <- c();
  for (i in 1:500) {
    out <- c();
    
    out$female <- as.data.frame(sample(female, size = nf[i], replace = FALSE))
    out$male   <- as.data.frame(sample(male, size = nf[i], replace = FALSE))
    
    out <- as.data.frame(out)
    names(out) <- c("female", "male")
    out$combo  <- interaction(out$female, out$male, sep = ".")
    
    z <- table(out$combo)
    
    pairs <- rbind(pairs, data.frame(macro.macro = z[1],
                                     macro.nour  = z[3],
                                     row.names = NULL))
    
  }
  
  focal <- rbind(focal, data.frame(ratio      = x[j],
                                   rareF.comM = pairs$macro.nour,
                                   total      = rowSums(pairs)))
}

#nour is rare
for (j in 9:length(x)) {
  #create female and male vectors
  female <- c(rep("macro", N/2 * x[j]), rep("nour", N/2 * y[j]))
  male   <- c(rep("macro", N/2 * x[j]), rep("nour", N/2 * y[j]))
  
  #draw random mating pairs
  nf <- sample(data$numf, size = 500, replace = TRUE)
  
  pairs <- c();
  for (i in 1:500) {
    out <- c();
    
    out$female <- as.data.frame(sample(female, size = nf[i], replace = FALSE))
    out$male   <- as.data.frame(sample(male, size = nf[i], replace = FALSE))
    
    out <- as.data.frame(out)
    names(out) <- c("female", "male")
    out$combo  <- interaction(out$female, out$male, sep = ".")
    
    z <- table(out$combo)
    
    pairs <- rbind(pairs, data.frame(nour.macro  = z[2],
                                     nour.nour   = z[4],
                                     row.names = NULL))
    
  }
  
  focal <- rbind(focal, data.frame(ratio      = x[j],
                                   rareF.comM = pairs$nour.macro,
                                   total      = rowSums(pairs)))
}

focal$ratio2 <- c(rep(0.1, 500), rep(0.15, 500),rep(0.2, 500), rep(0.25, 500), rep(0.3, 500), rep(0.35, 500), rep(0.4, 500), rep(0.45, 500), rep(0.5, 1000),
                  rep(0.45, 500), rep(0.4, 500), rep(0.35, 500), rep(0.3, 500), rep(0.25, 500), rep(0.2, 500), rep(0.15, 500), rep(0.1, 500))
focal <- focal[order(focal$ratio2, decreasing = FALSE), ]

#calculate 95% confidence interval
z <- seq(1, 9000, by = 1000)

conf <- c();
for (i in 1:length(z)) {
  #calculate mean incorrect and sample size
  k <- mean(focal[z[i]:(z[i] + 999), 2], na.rm = TRUE)
  n <- mean(focal[z[i]:(z[i] + 999), 3], na.rm = TRUE)
  
  #exact binomial test
  test <- binom.test(as.integer(k), as.integer(n), alternative = "two.sided", conf.level = 0.95)
  
  
  conf <- rbind(conf, data.frame(ratio            = focal[z[i], 1],
                                 mean.incorrect = test$estimate[1],
                                 lowerbound     = test$conf.int[1],
                                 upperbound     = test$conf.int[2]))
}



#difference between crosses
#0.1 fraction
t.test(data[7:9, 9], data[16:18, 9])

#0.25 fraction
t.test(data[4:6, 9], data[13:15, 9])

#0.50 fraction
t.test(data[1:3, 9], data[10:12, 9])
