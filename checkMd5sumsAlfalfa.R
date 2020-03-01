args <- commandArgs(TRUE)
print(args)

if(length(args)){
	md5calc <- args[[1]]
	md5True <- args[[2]]

	novaSeq <- read.table(md5calc)
	md5FromAdrien <- read.table(md5True)

	if(all(md5FromAdrien[[1]] == novaSeq[[1]])) cat("md5sums match!\n") else cat("Shoot dang! you gotta problem with these ones feller:\n", novaSeq[md5FromAdrien[[1]] != novaSeq[[1]], 2])
} else {
	stop("please provide the following:\npath to calculated md5sums\npath to true md5sums from Adrien\n")
}
