novaSeq <- read.table("alfalfaNovaSeq/alfalfaNovaSeqMd5sums.txt")
md5FromAdrien <- read.table("md5sum.source")

if(all(md5FromAdrien[[1]] == novaSeq[[1]])) cat("md5sums match!\n") else cat("Shoot dang! you gotta problem with these ones feller:\n", novaSeq[md5FromAdrien[[1]] != novaSeq[[1]], 2])

Robbins <- read.table("Robbins-NS-8397_2020_01_14/Robbins-NS-8397_2020_01_14Md5sums.txt")
md5FromAdrien <- read.table("md5sum.source")

if(all(md5FromAdrien[[1]] == Robbins[[1]])) cat("md5sums match!\n") else cat("Shoot dang! you gotta problem with these ones feller:\n", Robbins[md5FromAdrien[[1]] != Robbins[[1]], 2])

