suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(DESeq2))
suppressMessages(library(BiocParallel))
suppressMessages(library(pheatmap))
suppressMessages(library(Glimma))
suppressMessages(library(crayon))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(aplot))
library(tidyr)

args = commandArgs(T)
path = "D:\\Y2022\\Skinner2020R\\sRNASeq\\AASRA\\results\\VZ"

group_name <- args[1]
file_name <- paste0(path,args[2])
#case_num <- args[3]
#control_num <- args[4]
#tmp <- file.path(getwd(),group_name) 
tmp = "D:\\Y2022\\Skinner2020R\\sRNASeq\\AASRA\\results\\VZ"
cat(yellow(paste("###",group_name,"###","DESeq2 analysis start\n",sep = " ")))
cts <- fread("D:\\Y2022\\Skinner2020R\\sRNASeq\\AASRA\\rawcounts\\All2Groups.cts",header = T)

cts <-cts[which(cts$Geneid != ""),]
cts.row.name <- cts$Geneid

CR <- cts[ ,c("14C10-3-2-10_S39","14C10-3-3-9_S49","14C20-3-28-10_S76","14C20-3-28-8_S74","AC6-3-1-5_S4","AC9-3-4-8_S16","14C10-3-2-8_S37","14C10-3-5-13_S61","14C20-3-28-11_S85","14C20-3-28-9_S75","AC6-3-1-6_S13","AC9-3-4-9_S25","14C10-3-2-9_S38","14C17-3-4-10_S51","14C20-3-28-5_S63","AC2-3-3-10_S2","AC6-3-6-8_S14","AC9-3-5-8_S27","14C10-3-3-10_S50","14C17-3-5-14_S62","14C20-3-28-6_S64","AC2-3-3-9_S1","AC6-3-6-9_S15","AC9-3-5-9_S28","14C10-3-3-8_S40","14C17-3-5-9_S52","14C20-3-28-7_S73","AC6-3-1-4_S3","AC9-3-4-10_S26")] 

#AA <- cts[ ,c("AA13_S10","AA1_S1","AA26_S21","AA32_S26","AA3_S3","AA45_S37","AA5_S2","AA14_S11","AA20_S16","AA28_S22","AA35_S27","AA40_S32","AA46_S38","AA6_S4","AA15_S12","AA21_S17","AA29_S23","AA36_S28","AA41_S33","AA47_S39","AA7_S5","AA16_S13","AA22_S18","AA2_S2","AA37_S29","AA42_S34","AA48_S40","AA8_S6","AA17_S14","AA23_S19","AA30_S24","AA38_S30","AA43_S35","AA49_S41","AA9_S7","AA19_S15","AA24_S20","AA31_S25","AA39_S31","AA44_S36","AA4_S1")] 
#AA <- AA[ ,-c("AA31_S25","AA16_S13","AA7_S5", "AA14_S11")]

#DD <- cts [,c("14D10-3-5-10_S65","14D10-3-5-11_S66","14D10-3-5-6_S53","14D10-3-5-7_S54","14D10-3-5-8_S55","14D10-3-5-9_S56","14D11-3-6-10_S67","14D11-3-6-11_S68","14D11-3-6-12_S77","14D11-3-6-13_S72","14D12-3-7-7_S78","14D12-3-7-8_S79","14D13-3-8-10_S86","14D13-3-8-9_S80","14D15-3-23-2_S87","14D15-3-23-3_S71","14D15-3-23-4_S89","14D15-3-23-5_S90","14D15-3-23-6_S91","14D15-3-25-3_S92","14D15-3-25-4_S93","14D15-3-25-5_S94","14D15-3-25-6_S95","14D15-3-25-7_S96","14D15-3-25-8_S84","14D5-3-21-10_S31","14D5-3-21-6_S18","14D5-3-21-7_S19","14D5-3-21-8_S29","14D5-3-21-9_S30","14D5-3-3-10_S7","14D5-3-3-8_S5","14D5-3-3-9_S6","14D5-3-4-4_S8","14D5-3-4-5_S17","14D8-3-10-11_S42","14D8-3-10-12_S43","14D8-3-10-13_S44","14D8-3-10-6_S32","14D8-3-10-7_S41")] 

#JF <- cts[ ,c("JF10_S51","JF11_S52","JF12_S53","JF13_S54","JF14_S55","JF15_S56","JF16_S57","JF17_S58","JF18_S59","JF19_S60","JF1_S42","JF20_S61","JF21_S62","JF23_S64","JF24_S65","JF25_S66","JF26_S67","JF27_S68","JF28_S69","JF29_S70","JF2_S43","JF30_S71","JF31_S72","JF32_S73","JF33_S74","JF34_S3","JF35_S75","JF36_S76","JF37_S77","JF3_S44","JF40_S79","JF4_S45","JF5_S46","JF6_S47","JF8_S49","JF9_S50")] 
#JF <- cts[ ,-c("JF38_S4", "JF7_S48", "JF39_S78", "JF22_S63", "JF41_S80")]

VZ <- cts [,c("14V1-3-1-3_S9","14V1-3-1-4_S10","14V1-3-4-5_S11","14V1-3-4-6_S12","14V1-3-4-7_S20","14V2-3-3-7_S21","14V2-3-3-8_S22","14V3-3-2-5_S33","14V3-3-2-6_S23","14V3-3-8-6_S24","14V3-3-8-7_S34","14V6-3-6-6_S35","14V6-3-6-7_S36","14V6-3-7-2_S45","15V14-3-14-6_S46","15V14-3-14-7_S47","15V14-3-14-8_S48","15V14-3-14-9_S57","15V14-3-15-10_S58","15V14-3-15-11_S59","15V14-3-15-7_S60","15V14-3-15-8_S69","15V14-3-15-9_S70","15V15-3-17-9_S88","15V25-3-20-5_S81","15V25-3-20-6_S82","15V25-3-20-7_S83")] 


case_num <- length(unique(VZ))
control_num <- length(unique(CR))
cts <- cbind(VZ,CR)

dir.create(group_name)
cts <- as.matrix(cts)
rownames(cts) <- cts.row.name
colData <- data.frame(condition=1:length(colnames(cts)))
rownames(colData) <- colnames(cts)
colData[1:as.numeric(case_num),] <- "Case"
colData[(as.numeric(case_num)+1):(as.numeric(case_num)+as.numeric(control_num)),] <- "Control"
all(rownames(colData) == colnames(cts))
colData$condition <- factor(colData$condition, levels = c('Control','Case'))
dds <- DESeqDataSetFromMatrix(cts, colData, ~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref ="Control")
register(MulticoreParam(20))
dds <- DESeq(dds,parallel = TRUE)
res <- results(dds, contrast=c('condition','Case','Control'))
res <- res[order(res$pvalue),]
res_data <-data.frame(res)
res_data$id <- rownames(res_data)
res_data$id <- gsub("_ncrna$","",res_data$id)
#anno <- fread("human_smallRNA_reference.oneline.anno", header = T)
#res_data$Gene = anno[res_data$id,on= "ID"][,2]
rownames(res_data) <- 1:dim(res_data)[1]
#DESeq2 result
write.table(as.data.table(res_data), file.path(tmp, "DESeq2res.tsv"), sep="\t",row.names = FALSE,quote = FALSE)
cat(blue(paste("###",group_name,"###","DESeq2 result output is done\n",sep = " ")))
normalized_dds = data.frame(counts(dds, normalized=TRUE))
write.table(normalized_dds, file.path(tmp, "normalized_dds.tsv"),sep="\t",quote = FALSE)
cat(blue(paste("###",group_name,"###","Normalized counts output is done\n",sep = " ")))
#P value & adjust P value histogram
res_histogram <- as.data.table(res_data)[which(!is.na(res_data$padj)),]
h1 <- ggplot(res_histogram) +
  geom_histogram(aes(x = pvalue,y=..count../length(res_histogram$pvalue)),fill = "#2f5688" ,stat="bin", 
                 colour = "#e9ecef") +
  scale_y_continuous("P value Frequency", expand = c(0,0)) + 
  scale_x_continuous("Heterogeneity", expand = c(0,0)) +
  theme_classic() +
  theme(legend.position = 'none')
h2 <- ggplot(res_histogram) +
  geom_histogram(aes(x = padj,y=..count../length(res_histogram$padj)), fill = "#CC0000",stat="bin", 
                 colour = "#e9ecef") +
  scale_y_continuous("adjP value Frequency", expand = c(0,0)) + 
  scale_x_continuous("P value & adjust P value Heterogeneity", expand = c(0,0)) +
  theme_classic() +
  theme(legend.position = 'none')
pdf(file.path(tmp, "Histogram.pdf"))
h2%>%insert_top(h1+xlab(NULL),height = 1)
dev.off()
cat(blue(paste("###",group_name,"###","Histogram is done\n",sep = " ")))
#PCA plot
dds <- estimateSizeFactors(dds)
vsd <- vst(dds)
pdf(file.path(tmp, "PCA.pdf"))
plotPCA(vsd, intgroup=c("condition"))
dev.off()
cat(blue(paste("###",group_name,"###","PCA plot is done\n",sep = " ")))
#filterNumRej plot
pdf(file.path(tmp, "filterNumRej.pdf"))
plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
dev.off()
cat(blue(paste("###",group_name,"###","filterNumRej plot is done\n",sep = " ")))
#Disp plot
pdf(file.path(tmp, "Disp.pdf"))
plotDispEsts(dds)
dev.off()
cat(blue(paste("###",group_name,"###","Disp plot is done\n",sep = " ")))
#rawcounts boxplot
pdf(file.path(tmp, "boxplot.pdf"))
par(mar=c(8,5,2,2))
n.sample=ncol(counts(dds))
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
boxplot(log10(assays(dds)[["cooks"]]),col = cols, range=0, las=2,main="expression raw counts")
dev.off()
cat(blue(paste("###",group_name,"###","rawcounts boxplot is done\n",sep = " ")))
#Top & Bottom heatmap
data_heat=subset(res_data,select = c("log2FoldChange", "padj","id"))
data_heat<-na.omit(data_heat)
results_top = data_heat[order(data_heat$padj),]
results_top = results_top[which(results_top$log2FoldChange>0),][1:51,]
results_top <- data.table(results_top)
results_bottom = data_heat[order(data_heat$padj),]
results_bottom = results_bottom[which(results_bottom$log2FoldChange<0),][1:51,]
results_bottom <- data.table(results_bottom)
data_top =subset(results_top,select = "log2FoldChange")
data_top<-t(data_top)
colnames(data_top) <- results_top$id
pdf(file.path(tmp,"TopmRNAHeatMap.pdf"),width = 25)
pheatmap(data_top,
         show_rownames = T,
         show_colnames = T,
         scale = "none",
         cluster_rows = F,
         cellwidth = 9, cellheight = 15,
         main = "Top 50 RNAS",
         fontsize_col = 6
)
dev.off()
cat(blue(paste("###",group_name,"###","Top heatmap is done\n",sep = " ")))
data_bottom =subset(results_bottom,select = "log2FoldChange")
data_bottom<-t(data_bottom)
colnames(data_bottom) <- results_bottom$id
pdf(file.path(tmp,"BottomRNAHeatMap.pdf"),width = 25,)
pheatmap(data_bottom,
         show_rownames = T,
         show_colnames = T,
         scale = "none",
         cluster_rows = F,
         cellwidth = 9, cellheight = 15,
         main = "Bottom 50 RNAs",
         fontsize_col = 6
)
dev.off()
cat(blue(paste("###",group_name,"###","Bottom heatmap is done\n",sep = " ")))
#Glimma
status <- as.numeric(res$padj < 0.1)
anno <- data.frame(GeneID=rownames(res), symbol = res_data$id)
glMDPlot(res, status=status, counts=counts(dds,normalized=TRUE),
         groups=dds$condition, transform=FALSE,
         samples=colnames(dds), anno=anno,
         path=tmp, folder="glimma", launch=FALSE)
cat(blue(paste("###",group_name,"###","Glimma is done\n",sep = " ")))
#MA plot
pdf(file.path(tmp,"MA.pdf"))
plotMA(res, main="MA plot")
dev.off()
cat(blue(paste("###",group_name,"###","MA plot is done\n",sep = " ")))
#valcano plot
res_data$change <- as.factor(
  ifelse(
    res_data$padj<0.05 & abs(res_data$log2FoldChange)>2,
    ifelse(res_data$log2FoldChange>2, "Up-regulated", "Down-regulated"),
    "Not-significant"
  )
)
#anno <- fread("human_smallRNA_reference.oneline.anno", header = T)
#resdata$id = anno[res_data$id,on= "ID"][,2]
valcano <- ggplot(data=res_data, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
  geom_point(alpha=0.8, size=1) + 
  theme_bw(base_size=15) +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) +
  ylab("-log10(Adjust p-value)")+ 
  scale_color_manual(name="", values=c("#CC0000", "#2f5688", "#BBBBBB"), limits=c("Up-regulated", "Down-regulated","Not-significant")) + 
  geom_vline(xintercept=c(-2, 2), lty=2, col="black", lwd=0.5) + 
  geom_hline(yintercept=1.3, lty=2, col="black", lwd=0.5)
pdf(file.path(tmp,"Valcano.pdf"))
valcano
dev.off()
cat(blue(paste("###",group_name,"###","Valcano plot is done\n",sep = " ")))
cat(yellow(paste("###",group_name,"###","DESeq2 analysis finished\n",sep = " ")))

ls *cts|grep -v ^d| sed 's/.cts//g'> groupname
ls *cts > filename
paste groupname filename>config
cat config |while read id;do arr=${id};group=${arr[0]};file=${arr[1]};echo "Rscript DESeq2.r $group $file 2 2">>DECMD.sh;done
bash DECMD.sh
