# counts <- readLines("totCounts.txt")
# freqs <- readLines("refFreqs.txt")

counts <- readLines("align1to978_totCounts.txt")
freqs <- readLines("align1to978_refFreqs.txt")




length(counts)

freqH <- strsplit(freqs[1], " ")
countH <- strsplit(counts[1], " ")

freqTab <- do.call(rbind, lapply(strsplit(freqs[-1], " "), as.numeric))
colnames(freqTab) <- gsub("_.*", "", freqH[[1]])
countTab <- do.call(rbind, lapply(strsplit(counts[-1], " "), as.numeric))
colnames(countTab) <- gsub("_.*", "", countH[[1]])


totalCounts <- c(countTab)
hist(totalCounts, col = "gray")

totCountByLine <- colSums(countTab)
totCountByLine / max(totCountByLine)


avgCount <- rowMeans(countTab)
mean(avgCount)
avgCount > 20

txtdensity(avgCount)

# dens <- list()
# for(i in 1:ncol(countTab)) dens[[i]] <- density(countTab[, i])
# maxY <- max(sapply(dens, function(x) max(x$y)))


# histL <- list()
# for(i in 1:ncol(countTab)) histL[[i]] <- hist(countTab[, i], plot = FALSE)
# maxY <- max(sapply(histL, function(x) max(x$counts)))

# grays <- c("#b2b2b2b2", "#808080b2", "#333333b2")
# pdf("cntDist.pdf", width = 14)
# # par(mfrow(1, 3))
# plot(histL[[1]], xlab = "Site Read Depth", ylab = "density", ylim = c(0, maxY), main = "Site frequency", col = grays[1]) 
# for(i in 2:ncol(countTab)) plot(histL[[i]], add = TRUE, col = grays[i])
# legend("topright", legend = colnames(countTab), fill = grays)
# dev.off()
# system("scp cntDist.pdf Bender:~/alfalfaBioinfo/")


inRange20 <- countTab >= 20
allInRange20 <- apply(inRange20, 1, all)

sum(allInRange20 & avgCount > 10 & avgCount < 50)

keep <- allInRange20 & avgCount < 50
sum(keep)


nSites20 <- sum(allInRange20)
print(nSites20)



# inRange <- countTab >= 30 & countTab <= 70
# allInRange <- apply(inRange, 1, all)

# nSites <- sum(allInRange)
# print(nSites)

alFreq <- freqTab[keep, ]
Vaf <- var(alFreq)
# Vaf <- crossprod(alFreq)

pdf("VarAlFreq.pdf")
image(Vaf[,ncol(Vaf):1])
heatmap(Vaf)
dev.off()
system("scp VarAlFreq.pdf Bender:~/alfalfaBioinfo/")


expfr <- rowMeans(alFreq[, c("Chilean", "Mfalcata")])
expfrOri <- rowMeans(alFreq[, c("ChileanOriginal", "MfalcataOriginal")])

alFreqExp <- cbind(alFreq, Echmf = expfr, EchmfOri = expfrOri)




cor(alFreqExp)

cor(alFreq)

cor(alFreq[, c("Mfalcata", "MfalcataOriginal")])



library(txtplot)
txtplot(alFreq[,3], expfr)

# pdf("CHxMfFreqOE.pdf")
# plot(expfr, alFreq[,3], col = "#0000001A", xlab = "Expected Allele Frequency", ylab = "Observed Allele Frequency") # 4D 33 1A 0D
# dev.off()
# system("scp CHxMfFreqOE.pdf Bender:~/Dropbox/robbinsLabUpdates/Jan28_2020/")


png("CHxMfFreqOE_seedIncrease.png", 1080, 1080)
par(mar = c(5, 5, 2, 2))
r = cor(alFreqExp[, "Echmf"], alFreqExp[,"ChileanxMfalcata"])
plot(alFreqExp[, "Echmf"], alFreqExp[,"ChileanxMfalcata"], col = "#0000000D", xlab = "Expected Allele Frequency", ylab = "Observed Allele Frequency", cex = 2, cex.lab = 2, cex.axis = 2) # 4D 33 1A 0D
text(0.1, 0.9, labels = paste0("r = ", round(r, 2)), cex = 3)
text(0.8, 0.1, labels = paste0("# sites = ", nSites), cex = 3)
dev.off()
system("scp CHxMfFreqOE_seedIncrease.png Bender:~/alfalfaBioinfo/")


png("CHxMfFreqOE_original.png", 1080, 1080)
par(mar = c(5, 5, 2, 2))
r = cor(alFreqExp[, "EchmfOri"], alFreqExp[,"ChileanxMfalcata"])
plot(alFreqExp[, "EchmfOri"], alFreqExp[,"ChileanxMfalcata"], col = "#0000000D", xlab = "Expected Allele Frequency", ylab = "Observed Allele Frequency", cex = 2, cex.lab = 2, cex.axis = 2) # 4D 33 1A 0D
text(0.1, 0.9, labels = paste0("r = ", round(r, 2)), cex = 3)
text(0.8, 0.1, labels = paste0("# sites = ", nSites), cex = 3)
dev.off()
system("scp CHxMfFreqOE_original.png Bender:~/alfalfaBioinfo/")


png("ChileanOriVSSeedInc.png", 1080, 1080)
par(mar = c(5, 5, 2, 2))
r = cor(alFreqExp[,"ChileanOriginal"], alFreqExp[, "Chilean"])
plot(alFreqExp[,"ChileanOriginal"], alFreqExp[, "Chilean"], col = "#0000000D", xlab = "Original Allele Frequency", ylab = "Seed Increase Allele Frequency", cex = 2, cex.lab = 2, cex.axis = 2) # 4D 33 1A 0D
text(0.1, 0.9, labels = paste0("r = ", round(r, 2)), cex = 3)
text(0.8, 0.1, labels = paste0("# sites = ", nSites), cex = 3)
dev.off()
system("scp ChileanOriVSSeedInc.png Bender:~/alfalfaBioinfo/")


png("MfalcataOriVSSeedInc.png", 1080, 1080)
par(mar = c(5, 5, 2, 2))
r = cor(alFreqExp[,"MfalcataOriginal"], alFreqExp[, "Mfalcata"])
plot(alFreqExp[,"MfalcataOriginal"], alFreqExp[, "Mfalcata"], col = "#0000000D", xlab = "Original Allele Frequency", ylab = "Seed Increase Allele Frequency", cex = 2, cex.lab = 2, cex.axis = 2) # 4D 33 1A 0D
text(0.1, 0.9, labels = paste0("r = ", round(r, 2)), cex = 3)
text(0.8, 0.1, labels = paste0("# sites = ", nSites), cex = 3)
dev.off()
system("scp MfalcataOriVSSeedInc.png Bender:~/alfalfaBioinfo/")



txtdensity(alFreq[,1])
txtdensity(alFreq[,2])


pdf("CHxMfdensAf.pdf")
plot(density(alFreq[,1]), xlab = "Allele Frequency", main = "")
lines(density(alFreq[,2]), lty = 2)
legend("topleft", legend = c("Chilean", "M. falcata"), lty = 1:2)
dev.off()
system("scp CHxMfdensAf.pdf Bender:~/alfalfaBioinfo/")

head(countTab)

head(totC)
head(counts)
head(freqs)

