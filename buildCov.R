load("avgCountless50.RData")
library(txtplot)

useChilean <- "ChileanOriginal"
useMfalcata <- "MfalcataOriginal"

head(keepInfoTab)
head(keepCountTab)
head(keepFreqTab)


colMeans(keepCountTab)
colMeans(keepFreqTab)
avgFreq <- rowMeans(keepFreqTab)


for(i in seq(0, 1, 0.05)) print(paste(i, length(avgFreq[avgFreq <= i]), sep = "  "))
txtdensity(rowMeans(keepFreqTab))
keepByFreq <- avgFreq >= 0.05 & avgFreq <= 0.95


colnames(keepFreqTab)[colnames(keepFreqTab) %in% "Chilean"] <- "ChileanSeedInc"
colnames(keepFreqTab)[colnames(keepFreqTab) %in% "Mfalcata"] <- "MfalcataSeedInc"
colnames(keepCountTab)[colnames(keepCountTab) %in% "Chilean"] <- "ChileanSeedInc"
colnames(keepCountTab)[colnames(keepCountTab) %in% "Mfalcata"] <- "MfalcataSeedInc"

maxpos <- c(tapply(keepInfoTab[, "POS"], keepInfoTab[,"CHROM"], max))
class(maxpos) <- "numeric"
mPerbp <- c(table(keepInfoTab[,1])) / maxpos

maxpos["Super-Scaffold_596"]
summary(mPerbp)

txtdensity(mPerbp[mPerbp < 0.001])


dim(keepInfoTab)
dim(keepCountTab)
dim(keepFreqTab)




diallelParents <- 

combn()



parents <- c("African", "Chilean", "Flemish", "Indian", "Ladak", "Mfalcata", "Mvaria", "Peruvian", "Turkistan")



"MfalcataOriginal" 
"ChileanOriginal"  

"NY1430"
"NY1517"
"NY1518"
"OneidaVR"
"Regen"
"SW315LH"
"Ezra"
"Vernal" 

"AfricanxMfalcata"
"AfricanxPeruvian"
"ChileanxMfalcata"
"IndianxLadak"
"MfalcataxPeruvian"


"African"           
"AfricanxMfalcata"
"AfricanxPeruvian"
"ChileanOriginal"  
"Chilean"           
"ChileanxMfalcata"
"Ezra"
"Flemish"          
"Indian"            
"IndianxLadak"
"Ladak"
"MfalcataOriginal" 
"Mfalcata"          
"MfalcataxPeruvian"
"Mvaria"
"NY1430"
"NY1517"
"NY1518"
"OneidaVR"
"Peruvian"         
"Regen"
"SW315LH"
"Turkistan"
"Vernal" 