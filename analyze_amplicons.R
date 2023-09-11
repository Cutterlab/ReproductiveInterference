library(dplyr)
library(tidyverse)

#NOTES:
#Amplicon length is 378bp and both primers are shared between species
#There are 312 identical bases between species (82% pairwise identity)


#length distribution of raw reads
setwd("~/Documents/UofT/experiments/Rebecca_ReproductiveInterference/RS_amplicon_round3/reads/")
lengths <- 
  list.files(pattern = "*.txt") %>% 
  map_df(~read.table(., header = FALSE, sep = "\t"))
names(lengths) <- c("ID", "bp")

plot(density(lengths$bp), las = 1)
hist(subset(lengths, bp <= 50)$bp)


#select the best BLAST hit based on e-value and bit score
setwd("~/Documents/UofT/experiments/Rebecca_ReproductiveInterference/RS_amplicon_round2/RS_BestBlast/")
allfiles = list.files(pattern = "*.txt")

final <- c();
identical <- c();
for (f in 1:length(allfiles)) {
  round <- read.table(allfiles[f], header = FALSE, sep = "\t")
  names(round) <- c("ID", "amplicon", "identity", "alignment.length", "mismatches", "gap.opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit.score")
  round$replicate <- substr(allfiles[f], start = 1, stop = 9)
  
  best.round <- c();
  for(i in 1:length(round[duplicated(round$ID), 1])){
    A <- subset(round, ID == round[duplicated(round$ID), 1][i])
    B <- subset(A, evalue == min(evalue) & max(bit.score))
    C <- length(B$ID) > 1
    
    if (C == FALSE){
      best.round <- rbind(best.round, B) 
    }
    else {
      identical <- rbind(identical, B)
    }
  }
  
  D <- round[!duplicated(round$ID), ]
  
  final <- rbind(final, D, best.round)
}

write.table(final, "FINAL_BestBlast.txt", row.names = FALSE, quote = FALSE, sep = "\t")

#data check
amplicons <- read.table("Documents/UofT/experiments/Rebecca_ReproductiveInterference/RS_amplicon_round2/FINAL_BestBlast.txt", header = TRUE, sep = "\t")
qc <- as.data.frame(summarise(group_by(amplicons, amplicon), mean(identity), sd(identity), mean(alignment.length), sd(alignment.length), min(alignment.length), mean(mismatches), sd(mismatches)))
#total reads: 4,443,793
#macrosperma: 2,962,567
#nouraguensis: 1,481,226

amplicons$replicate2 <- substr(amplicons$replicate, start = 1, stop = 6)
amplicons$replicate2[amplicons$replicate2 == "Ctl_R1"] <- "CTL"
amplicons$replicate2[amplicons$replicate2 == "Ctl_R2"] <- "CTL"

#tabulate allele counts and calculate frequencies
freq.table <- c();
for (i in 1:length(levels(as.factor(amplicons$replicate2)))) {
  A <- amplicons[amplicons$replicate2 == levels(as.factor(amplicons$replicate2))[i], ];
  
  B <- as.data.frame(table(A$amplicon))
  
  freq.table <- rbind(freq.table, data.frame(ID    = A$replicate2[1],
                                             MACRO = B[1, 2],
                                             NOUR  = B[2, 2]))
}

freq.table$TOTAL   <- rowSums(freq.table[, 2:3])
freq.table$F.MACRO <- freq.table$MACRO / freq.table$TOTAL
freq.table$F.NOUR  <- freq.table$NOUR / freq.table$TOTAL
freq.table$DAY     <- c(-1, rep(11, 15), rep(24, 15), rep(4, 15))
freq.table$TREATMENT <- c("CTL", substr(freq.table[2:31, 1], start = 6, stop = 6), substr(freq.table[32:46, 1], start = 5, stop = 5))
freq.table$REP <- c("CTL", substr(freq.table[2:31, 1], start = 4, stop = 5), substr(freq.table[32:46, 1], start = 3, stop = 4))


ids  <- c("D1R1C", "D1R2C", "D1R3C", "D1R4C", "D1R5M", "D1R1M", "D1R2M", "D1R3M", "D1R4M", "D1R5M", "D1R1N", "D1R2N", "D1R3N", "D1R4N", "D1R5N")
wrms <- c(rep(500, 5), rep(800, 5), rep(200, 5))
trts <- c(rep("C", 5), rep("M", 5), rep("N", 5))
rps  <- c(rep(c("R1", "R2", "R3", "R4", "R5"), 3))
for (j in 1:length(ids)) {
  freq.table <- rbind(freq.table, data.frame(ID = ids[j],
                                             MACRO = wrms[j],
                                             NOUR  = (1000 - wrms[j]),
                                             TOTAL = 1000,
                                             F.MACRO = (wrms[j] / 1000),
                                             F.NOUR  = ((1000 - wrms[j]) / 1000),
                                             DAY = 1,
                                             TREATMENT = trts[j],
                                             REP = rps[j]))
}


#estimate competition coefficient favoring macrosperma (positive slope) or favoring nouraguensis (negative slope)
competition.CTL <- glm(cbind(MACRO, NOUR) ~ DAY, family = binomial(link = "logit"), data = subset(freq.table, TREATMENT == "C"))
summary(competition.CTL)

competition.M <- glm(cbind(MACRO, NOUR) ~ DAY, family = binomial(link = "logit"), data = subset(freq.table, TREATMENT == "M"))
summary(competition.M)

competition.N <- glm(cbind(MACRO, NOUR) ~ DAY, family = binomial(link = "logit"), data = subset(freq.table, TREATMENT == "N"))
summary(competition.N)
