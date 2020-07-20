# load("avgCountless50.RData")
pardir <- "/workdir/ns722/alfalfaVariantCalls/"
pardir <- "~/alfalfaBioinfo/dataFiles/"

info <- read.table(paste0(pardir, "alignAll_namePos_cnt20-125_mean20-75_tot0-10000.txt"), comment.char = "", header = TRUE)
freq <- read.table(paste0(pardir, "alignAll_refFreqs_cnt20-125_mean20-75_tot0-10000.txt"), header = TRUE)
cnt <- read.table(paste0(pardir, "alignAll_totCounts_cnt20-125_mean20-75_tot0-10000.txt"), header = TRUE)


# this was built in Dropbox/pairwiseFst/codis.R, varified using results from Weir and Hill 2002

# 4/23/2020: calc dominance by assuming HW equil, then use freq of exp hets in place of counts? Would have to assume number of indiviudals???


calcFstLin <- function(altab, scale = FALSE){
	one <- function(n) matrix(1, n, 1)
	onet <- function(n) matrix(1, 1, n)
	
	r <- nrow(altab)
	m <- ncol(altab)

	ni <- rowSums(altab)
	nic <- ni - ni^2 / sum(ni)
	nc <- 1/(r-1) * sum(nic)

	piu <- altab / rowSums(altab)
	cnts <- colSums(altab)
	pbar <- 1/sum(ni) * cnts

	num <- sum(nic) * ((ni / (ni-1)) %x% onet(m) * (piu * (1-piu))) %*% one(m)
	den <- (onet(m) %*% t(ni * (piu - one(r) %x% t(pbar))^2 + nic * (piu * (1-piu))) %*% one(r))[[1]]
	Bi <- c(1 - num / den)
	# print(Bi)
	Bii <- 1 - sum(nic) *  2 * tcrossprod(piu, 1-piu) / (2*den)
	Bii <- Bii[upper.tri(Bii)]
	list(Bi = Bi, Bii = Bii)
}

# save(info, freq, cnt, file = "tempCalls.RData")
# load("tempCalls.RData")
dim(cnt)
# colnames(info) <- c("CHROM", "POS", "ID", "REF", "ALT")
# names(info) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL")


head(info)
head(freq)
head(cnt)

freq <- as.matrix(freq)
cnt <- as.matrix(cnt)

colnames(freq) <- gsub("_sorted", "", colnames(freq))
colnames(cnt) <- gsub("_sorted", "", colnames(cnt))


goodQual <- info$QUAL > 20
info <- info[goodQual, ]
cnt <- cnt[goodQual, ]
freq <- freq[goodQual, ]


library(txtplot)
parents <- c("African", "Chilean", "Flemish", "Indian", "Ladak", "Mfalcata", "Mvaria", "Peruvian", "Turkistan")

useChilean <- "ChileanOriginal"
useMfalcata <- "MfalcataOriginal"

colnames(freq)[colnames(freq) %in% "Chilean"] <- "ChileanSeedInc"
colnames(freq)[colnames(freq) %in% "Mfalcata"] <- "MfalcataSeedInc"
colnames(cnt)[colnames(cnt) %in% "Chilean"] <- "ChileanSeedInc"
colnames(cnt)[colnames(cnt) %in% "Mfalcata"] <- "MfalcataSeedInc"

colnames(freq)[colnames(freq) %in% useChilean] <- "Chilean"
colnames(freq)[colnames(freq) %in% useMfalcata] <- "Mfalcata"
colnames(cnt)[colnames(cnt) %in% useChilean] <- "Chilean"
colnames(cnt)[colnames(cnt) %in% useMfalcata] <- "Mfalcata"


parCnt <- cnt[, parents]
parFreq <- freq[, parents]

freqmu <- rowMeans(parFreq)
all1or0 <- freqmu < 1 & freqmu > 0 
sum(all1or0) / nrow(parFreq)

cnt <- cnt[all1or0, ]
freq <- freq[all1or0, ]
info <- info[all1or0, ]
parCnt <- parCnt[all1or0, ]
parFreq <- parFreq[all1or0, ]

# colnames(cnt)

meanCount <- rowMeans(cnt)
meanFreq <- rowMeans(freq) # Note, this is not the actual population mean because it weights all populations equally. But all theoretical populations have similar Ne?
txtdensity(meanFreq)
sum(meanFreq < 0.9)
sum(meanFreq > 0.1)
# meanFreq < 0.95

hybridFreq <- matrix(NA, nrow(parFreq), choose(ncol(parFreq), 2))
hybridExpCnt <- matrix(NA, nrow(parCnt), choose(ncol(parCnt), 2))
hybridNames <- NULL
index <- 1
for(i in 1:(ncol(parFreq)-1)){
	for(j in (i+1):ncol(parFreq)){
		hybridFreq[, index] <- (parFreq[, i] + parFreq[, j]) / 2
		hybridExpCnt[, index] <- (cnt[, i] + cnt[, j]) / 2
		hybridNames <- c(hybridNames, paste0(colnames(parFreq)[i], "x", colnames(parFreq)[j], "Est"))
		index <- index + 1
	}
}

colnames(hybridFreq) <- hybridNames
colnames(hybridExpCnt) <- hybridNames

head(hybridFreq)
head(hybridExpCnt)

allFreq <- cbind(freq, hybridFreq) 
allCnt <- cbind(cnt, hybridExpCnt) 

refCnt <- allFreq * allCnt
altCnt <- allCnt - refCnt


useFreq <- meanFreq < 0.9 & meanFreq > 0.1
sum(useFreq)

useNames <- info[useFreq, ]
a1 <- refCnt[useFreq, ]
a2 <- altCnt[useFreq, ]

BiL <- list()
BiiL <- list()
for(i in 1:nrow(a1)){
	altabi <- cbind(a1[i, ], a2[i, ])
	fst <- calcFstLin(altabi)
	BiL[[i]] <- fst[["Bi"]]
	BiiL[[i]] <- fst[["Bii"]]
}



BiTab <- do.call(rbind, BiL)
BiiTab <- do.call(rbind, BiiL)

B <- diag(colMeans(BiTab))
B[upper.tri(B)] <- colMeans(BiiTab)
B[lower.tri(B)] <- t(B)[lower.tri(B)]



colnames(B) <- rownames(B) <- colnames(a1)

library(viridis)
image(B[, nrow(B):1], col = magma(50))

pdf("pairwiseFstHeatmap.pdf")
par(oma = c(2, 2, 2, 2))
# heatmap(B, col = magma(50), cexRow = 0.5, cexCol = 0.5)
heatmap(B,  Colv = "Rowv", symm = TRUE, col = magma(50), cexRow = 0.5, cexCol = 0.5)
dev.off()

write.table(B, "pairwiseFstAlfalfa.txt", sep = "\t")



# make covariance based on allele frequencies, scale diagonal to 1
keepFreq <- allFreq[useFreq, ]
dim(keepFreq)
Kfreq <- cov(keepFreq)
Kfreq <- Kfreq / mean(diag(Kfreq))


pdf("covAlleleFreqHeatmap.pdf")
par(oma = c(2, 2, 2, 2))
# heatmap(B, col = magma(50), cexRow = 0.5, cexCol = 0.5)
heatmap(Kfreq,  Colv = "Rowv", symm = TRUE, col = magma(50), cexRow = 0.5, cexCol = 0.5)
dev.off()

write.table(Kfreq, "covAlleleFreqAlfalfa.txt", sep = "\t")


# hetHW <- function(p) 2 * p*(1-p)

# for(i in 1:nrow(freq)){
# 	freq[i, ] * 
# }



# p1 <- freq[21, c("Chilean")]
# p2 <- freq[21, c("Mfalcata")]
# pbar <- meanFreq[21]
# getThetaiiLarge <- function(p1, p2) {
# 	1 - {(p1 * (1-p2) + p2 * (1-p1)) / (2 * pbar * (1-pbar))}
# }


# getThetaiiLarge <- function(p1, p2) {
# 	1 - {(p1 * (1-p2) + p2 * (1-p1)) / (2 * pbar * (1-pbar))}
# 	1 - {(p1 * (1-p2) + p2 * (1-p1) + (1-p1) * p2 + p2 * (1-p1)) / (2 * pbar * (1-pbar) + 2 * (1-pbar) * pbar)}
# }





# getThetaiiLarge <- function(p1, p2, pbar) {
# 	if(length(p1) != length(p2) | length(p1) != length(pbar)) stop("all populations must have all alleles present! if missing in population, use zero to indicate.")
# 	calcLastFreq <- if(sum(p1) < 1 | sum(p2) < 1 | sum(pbar) < 1) TRUE else FALSE

# 	num <- 0
# 	den <- 0
# 	for (i in 1:length(p1)){
# 		num <- num + p1[i] * (1-p2[i]) + p2[i] * (1-p1[i])
# 		den <- den + 2 * pbar[i] * (1-pbar[i])
# 	}

# 	if(calcLastFreq){
# 		p1l <- 1 - sum(p1)
# 		p2l <- 1 - sum(p2)
# 		pbarl <- 1 - sum(pbar)

# 		num <- num + p1l * (1-p2l) + p2l * (1-p1l)
# 		den <- den + 2 * pbarl * (1-pbarl)

# 	}
# 	1 - {num / den}
# }




# # want Fst or 1- Fst? Fst = var shared. 

# r <- ncol(parHybridFreq)
# n <- rep(100, r)
# nbar <- mean(n) # avg sample size
# nc <- (r * nbar - sum(n^2) / (r * nbar)) / (r - 1)
# # pbar <- n %*% t(parHybridFreq) / (r * nbar)
# pbar <- parHybridFreq %*% n / (r * nbar)

# s2 <- (parHybridFreq - matrix(1, 1, r) %x% pbar)^2 %*% n / ((r-1) * nbar)
# hbar <- hetHW(parHybridFreq) %*% n / (r * nbar)
# # pairwise Fst


# a <- nbar / nc * (s2 - 1/(nbar-1) * (pbar * (1 - pbar) - (r-1)/r * s2 - 1/4 * hbar))
# dim(a)
# txtdensity(a)
# b <- nbar / (nbar - 1) * (pbar*(1-pbar) - (r-1) / r * s2 - (2*nbar - 1) / (4 * nbar) * hbar)
# dim(b)
# txtdensity(b)
# cee <- 1/2 * hbar

# abc <- a + b + cee
# theta = a / a + b + cee

# txtdensity(theta)


# for(i in unique(info)){

# }

# meanVar <- meanFreq * (1 - meanFreq)
# txtdensity(meanVar)

# for(i in 1:ncol(freq)){
# 	for(j in (i+1)):ncol(freq)){
# 		freq[,i]

# 	}
# }


# head(keepInfoTab)
# head(keepCountTab)
# head(keepFreqTab)


# colMeans(keepCountTab)
# colMeans(keepFreqTab)
# avgFreq <- rowMeans(keepFreqTab)


# for(i in seq(0, 1, 0.05)) print(paste(i, length(avgFreq[avgFreq <= i]), sep = "  "))
# txtdensity(rowMeans(keepFreqTab))
# keepByFreq <- avgFreq >= 0.05 & avgFreq <= 0.95




# maxpos <- c(tapply(keepInfoTab[, "POS"], keepInfoTab[,"CHROM"], max))
# class(maxpos) <- "numeric"
# mPerbp <- c(table(keepInfoTab[,1])) / maxpos

# maxpos["Super-Scaffold_596"]
# summary(mPerbp)

# txtdensity(mPerbp[mPerbp < 0.001])


# dim(keepInfoTab)
# dim(keepCountTab)
# dim(keepFreqTab)




# diallelParents <- 

# combn()






# "MfalcataOriginal" 
# "ChileanOriginal"  

# "NY1430"
# "NY1517"
# "NY1518"
# "OneidaVR"
# "Regen"
# "SW315LH"
# "Ezra"
# "Vernal" 

# "AfricanxMfalcata"
# "AfricanxPeruvian"
# "ChileanxMfalcata"
# "IndianxLadak"
# "MfalcataxPeruvian"


# "African"           
# "AfricanxMfalcata"
# "AfricanxPeruvian"
# "ChileanOriginal"  
# "Chilean"           
# "ChileanxMfalcata"
# "Ezra"
# "Flemish"          
# "Indian"            
# "IndianxLadak"
# "Ladak"
# "MfalcataOriginal" 
# "Mfalcata"          
# "MfalcataxPeruvian"
# "Mvaria"
# "NY1430"
# "NY1517"
# "NY1518"
# "OneidaVR"
# "Peruvian"         
# "Regen"
# "SW315LH"
# "Turkistan"
# "Vernal" 