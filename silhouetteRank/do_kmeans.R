freq_file<-commandArgs(trailingOnly=T)[1]
par_seed <-commandArgs(trailingOnly=T)[2]
par_k <-commandArgs(trailingOnly=T)[3]
nstart<-commandArgs(trailingOnly=T)[4]
centroid_file<-commandArgs(trailingOnly=T)[5]
kmeans_file<-commandArgs(trailingOnly=T)[6]

par_k<-as.integer(par_k)
par_seed<-as.integer(par_seed)
nstart<-as.integer(nstart)

if(par_seed!=-1 & par_seed>0){
set.seed(par_seed)
}

xx<-read.table(freq_file, sep=" ", header=F)
y<-c(); for(i in seq(1, dim(xx)[1])){y<-append(y, rep(xx[i,2], xx[i,1]))}
kk<-kmeans(y, par_k, nstart=nstart, iter.max=300)
write.table(kk$cluster, file=kmeans_file, sep=" ", quote=F, col.names=F, row.names=T)
write.table(kk$centers, file=centroid_file, sep=" ", quote=F, col.names=F, row.names=T)
