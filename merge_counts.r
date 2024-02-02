library(data.table)
args = commandArgs(T)
path = "~/fcounts/"
res <- {}
for (i in 1:(length(args)-1)) {
    sample <- args[i]
    filename <- paste0(path,sample,".sam.tsv")
    dat <- fread(filename,header = T,skip = 1)
    dat <- dat[,c(1,7)]
    names(dat) <- c("Geneid",sample)
    res <- cbind(res,dat)
}
res <- res[,-seq(1,ncol(res),2)[-1],with=FALSE]
outname <- tail(args,1)
write.table(res,outname, sep="\t",row.names = FALSE,quote = FALSE)
