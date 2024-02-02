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
path = "D:\\Y2022\\Skinner2020R\\Atrazine\\results"

group_name <- args[1]
file_name <- paste0(path,args[2])
#case_num <- args[3]
#control_num <- args[4]
#tmp <- file.path(getwd(),group_name) 
tmp = "D:\\Y2022\\Skinner2020R\\Atrazine\\results"
cat(yellow(paste("###",group_name,"###","DESeq2 analysis start\n",sep = " ")))
cts <- fread("D:\\Y2022\\Skinner2020R\\lRNASeq\\All.cts",header = T)

cts <-cts[which(cts$Geneid != ""),]
cts.row.name <- cts$Geneid

CR <- cts[ ,c("14C10_3_2_10","14C10_3_2_8","14C10_3_2_9","14C10_3_3_10","14C10_3_3_8","14C10_3_3_9","14C10_3_5_13","14C17_3_4_10","14C17_3_5_14","14C17_3_5_9","14C20_3_28_10","14C20_3_28_11","14C20_3_28_5","14C20_3_28_6","14C20_3_28_7","14C20_3_28_8","14C20_3_28_9","AC2_3_3_10","AC2_3_3_9","AC6_3_1_4","AC6_3_1_5","AC6_3_1_6","AC6_3_6_8","AC6_3_6_9","AC9_3_4_10","AC9_3_4_8","AC9_3_4_9","AC9_3_5_8","AC9_3_5_9")] 
CR <- CR[,-c("14C20_3_28_11", "14C10_3_3_9", "14C20_3_28_7")]

AA <- cts[ ,c("AA2_3_11_5","AA2_3_11_6","AA2_3_11_7","AA2_3_11_8","AA2_3_11_9","AA2_3_6_8","AA2_3_6_9","AA2_3_9_10","AA2_3_9_11","AA2_3_9_6","AA2_3_9_7","AA2_3_9_8","AA3_3_10_2","AA3_3_1_5","AA3_3_1_6","AA3_3_2_8","AA3_3_8_11","AA3_3_8_12","AA3_3_8_13","AA4_3_4_3","AA4_3_4_5","AA4_3_4_6","AA8_3_2_10","AA8_3_2_7","AA8_3_2_8","AA8_3_2_9","AA9_3_12_6","AA9_3_12_7","AA9_3_12_8","AA9_3_12_9","AA9_3_5_10","AA9_3_5_7","AA9_3_5_8","AA9_3_5_9","AA9_3_7_4","AA9_3_7_5","AA9_3_7_6","AA9_3_7_7","AA9_3_7_8","AA9_3_7_9")] 
AA <- AA[ ,-c("AA2_3_11_9", "AA2_3_11_8", "AA4_3_4_5")]

#DD <- cts [,c("14D10_3_5_10","14D10_3_5_11","14D10_3_5_7","14D10_3_5_8","14D10_3_5_9","14D11_3_6_10","14D11_3_6_11","14D11_3_6_12","14D12_3_7_7","14D12_3_7_8","14D13_3_8_10","14D13_3_8_9","14D15_3_23_2","14D15_3_23_3","14D15_3_23_4","14D15_3_23_5","14D15_3_23_6","14D15_3_25_3","14D15_3_25_4","14D15_3_25_5","14D15_3_25_6","14D15_3_25_7","14D15_3_25_8","14D5_3_21_6","14D5_3_21_7","14D5_3_3_10","14D5_3_3_8","14D5_3_3_9","14D5_3_4_4","14D5_3_4_5","14D8_3_10_11","14D8_3_10_12","14D8_3_10_13","14D8_3_10_7")] 
#DD <- DD[ ,-c("14D15_3_25_6", "14D12_3_7_7")]

#JF <- cts [,c("DJA0_3_12_12","DJA0_3_12_13","DJA0_3_12_14","DJL2_3_11_6","DJL2_3_11_7","DJL2_3_11_9","DJL2_3_1_10","DJL2_3_1_11","DJL2_3_7_7","DJL2_3_7_9","DJT1_3_6_10","DJT1_3_6_7","DJT1_3_6_8","DJT1_3_6_9","DJW0_3_5_6","DJW0_3_5_7","DJW0_3_5_8","DJW0_3_5_9","DJZ2_3_10_6","DJZ2_3_10_7","DJZ2_3_10_8","DJZ2_3_2_3","DJZ2_3_3_15","DJZ2_3_9_10","DJZ2_3_9_11","DJZ2_3_9_12")] 

#VZ <- cts [,c("14V1_3_1_3","14V1_3_1_4","14V1_3_4_5","14V1_3_4_6","14V1_3_4_7","14V2_3_3_7","14V2_3_3_8","14V3_3_2_5","14V3_3_2_6","14V3_3_8_6","14V3_3_8_7","14V6_3_6_6","14V6_3_6_7","14V6_3_7_2","15V14_3_14_6","15V14_3_14_7","15V14_3_14_8","15V14_3_14_9","15V14_3_15_10","15V14_3_15_11","15V14_3_15_7","15V14_3_15_8","15V14_3_15_9","15V15_3_17_9","15V25_3_20_5","15V25_3_20_6","15V25_3_20_7")] 
#VZ <- VZ[ ,-c("14V1_3_4_6")]

case_num <- length(unique(AA))
control_num <- length(unique(CR))
cts <- cbind(AA,CR)

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
