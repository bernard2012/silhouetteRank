#!/usr/bin/Rscript
c1<-commandArgs(TRUE)[1]
c2<-commandArgs(TRUE)[2]
library(qvalue)
x<-read.table(c1, header=FALSE)
#print(length(x[,1]))

do_qval <- function(yy){
	out <- tryCatch({
		qvalue(x[,1])
	}, error=function(cond){
		qvalue(x[,1], pi0=1)
	})
	return(out)
}

y <- do_qval(x[,1])
write.table(y$qvalues, file=c2, sep="\n", col.names=FALSE, row.names=FALSE)

