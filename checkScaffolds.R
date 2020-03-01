# scaffs <- read.table("scaffoldLengthInfo.txt", nrows = 5085)
scaffs <- read.table("scaffoldLengthInfo.txt")
names(scaffs) <- c("scaffold", "length")

scaffs$num <- as.numeric(gsub("obj", "-1", gsub(".*_", "", scaffs[["scaffold"]])))
# scaffs$num <- gsub(".*_", "", gsub("_obj|>tig", "", scaffs[["scaffold"]]))


library(txtplot)
txtdensity(scaffs$length[scaffs$num < 1000 & scaffs$num > 0])


sum(scaffs$length[scaffs$num < 1000 & scaffs$num > 0])  / 1e9
sum(scaffs$length[grep("100000", scaffs$num)]) / 1e6
sum(scaffs$length[grep("200000", scaffs$num)])  / 1e6


cat alfalfa3a.2/3a.2-bionano+unscaffolded-assembly.fa | awk '{if(NR%2==0) print length($1)}' > scaffoldLengths.txt
paste -d" " scaffoldNames.txt scaffoldLengths.txt > scaffoldLengthInfo.txt



Scaffolds <- nrow(scaffs)


MeanScaffold = mean(scaffs$length[scaffs$num > 0]),
ScaffoldN50 = quantile(scaffs$length[scaffs$num > 0], 0.5)[[1]],
ScaffoldN90 = quantile(scaffs$length[scaffs$num > 0], 0.1)[[1]],

# # from Connor:
# Scaffolds       10169
# Max Scaffold    16445029
# Min Scaffold    4
# Mean Scaffold   1956842
# Scaffold N50    2207780
# Scaffold N90    962372
# Total Scaffold Length   2698044700

genstats <- c(Scaffolds = nrow(scaffs),
MaxScaffold = max(scaffs$length),
MinScaffold = min(scaffs$length),
MeanScaffold = mean(scaffs$length),
ScaffoldN50 = quantile(scaffs$length, 0.5)[[1]],
ScaffoldN90 = quantile(scaffs$length, 0.1)[[1]],
TotalScaffoldLength = sum(scaffs$length))
round(matrix(genstats, ncol = 1, dimnames = list(names(genstats), "")))
# Scaffolds                10169
# MaxScaffold           16445029
# MinScaffold                  4
# MeanScaffold            265321
# ScaffoldN50              54857
# ScaffoldN90              10888
# TotalScaffoldLength 2698044700


genstats2 <- c(Scaffolds = nrow(scaffs),
MaxScaffold = max(scaffs$length),
MinScaffold = min(scaffs$length),
MeanScaffold = mean(scaffs$length[scaffs$num > 0]),
ScaffoldN50 = quantile(scaffs$length[scaffs$num > 0], 0.5)[[1]],
ScaffoldN90 = quantile(scaffs$length[scaffs$num > 0], 0.1)[[1]],
TotalScaffoldLength = sum(scaffs$length))
round(matrix(genstats2, ncol = 1, dimnames = list(names(genstats2), "")))

# Scaffolds                10169
# MaxScaffold           16445029
# MinScaffold                  4
# MeanScaffold           1505257
# ScaffoldN50             857229
# ScaffoldN90             136195
# TotalScaffoldLength 2698044700




genstats3 <- c(Scaffolds = nrow(scaffs),
MaxScaffold = max(scaffs$length),
MinScaffold = min(scaffs$length),
MeanScaffold = mean(scaffs$length[scaffs$num > 0 & scaffs$num < 1000]),
ScaffoldN50 = quantile(scaffs$length[scaffs$num > 0 & scaffs$num < 1000], 0.5)[[1]],
ScaffoldN90 = quantile(scaffs$length[scaffs$num > 0 & scaffs$num < 1000], 0.1)[[1]],
TotalScaffoldLength = sum(scaffs$length))
round(matrix(genstats3, ncol = 1, dimnames = list(names(genstats3), "")))

Scaffolds                10169
MaxScaffold           16445029
MinScaffold                  4
MeanScaffold           2092487
ScaffoldN50            1506560
ScaffoldN90             451995
TotalScaffoldLength 2698044700

meanL <- NULL
for(i in 1:nrow(scaffs)) {
	meanL[[i]] <- mean(scaffs[1:i, "length"])
}
which.min(abs(meanL - 1956842))