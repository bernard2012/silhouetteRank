library(eva)
f_name<- commandArgs(trailingOnly = T)[1]
x<-t(t(read.table(paste0(f_name), header=F)))
y<-gpdFit(x, nextremes=250, method="mle")
write.table(y$par.ests, file=paste0("par.", f_name), sep="\t", col.names=F)

