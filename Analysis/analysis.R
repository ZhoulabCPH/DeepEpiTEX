###Meth:consistent up and down across five TEX subsets
meth<-read.table("Methyl_m.txt")

meth_up_gene<-read.table("gene_up.txt")
meth_down_gene<-read.table("gene_down.txt")

meth<-t(meth)
meth<-as.data.frame(meth)

meth_up<-meth[rownames(meth)%in%meth_up_gene[,1],]
meth_down<-meth[rownames(meth)%in%meth_down_gene[,1],]




label<-read.table("Methyl_l.txt")
rownames(label)<-label[,1]

label<-label[order(label$V2),]
meth_up<-meth_up[,rownames(label)]
meth_down<-meth_down[,rownames(label)]



annotation_col = data.frame( cluster = as.factor(label[,2]))
rownames(annotation_col) = rownames(label)
ann_colors = list(cluster = c('1' = "#0B5F99", 
                              "2" = "#86B5CE",
                              "3" = "#F2D5C1",
                              "4" = "#DE9B7D", 
                              "5" = "#B22023"))

blue <- rgb(63,117,179,maxColorValue = 255)
red <- rgb(214,47,36,maxColorValue = 255)
white <- rgb(255,255,255,maxColorValue = 255)


jieguo<-meth_up
jieguo<-meth_down

linshi <- apply(jieguo,1,scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(jieguo)
hist(linshi)
linshi[linshi>2] <- 2
linshi[linshi<(-2)] <- c(-2)
linshi<-as.data.frame(linshi)


out<- pheatmap(linshi,fontsize=6,cutree_col = 5,
               color  = colorRampPalette(c(blue,white,red))(100),
               annotation_col = annotation_col,
               annotation_colors = ann_colors,
               clustering_method = "ward.D",
               # border_color = "grey60",
               cluster_cols = F, cluster_rows = T,
               show_rownames = T, show_colnames = T,
               border=FALSE
)



meth_up<-t(meth_up)
meth_up<-as.data.frame(meth_up)
label<-label[rownames(meth_up),]
meth_up<-cbind(label[,2],meth_up)
colnames(meth_up)[1]<-c("label")

data<-aggregate(apply(meth_up[,2:605],2,as.numeric),by=list(meth_up$label),mean)
rownames(data)<-data[,1]
data<-data[,-1]

jieguo<-t(data)
jieguo<-as.data.frame(jieguo)

jieguo = jieguo[apply(jieguo, 1, function(x) sd(x)!=0),]   
linshi <- apply(jieguo,1,scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(jieguo)
hist(linshi)
linshi[linshi>2] <- 2
linshi[linshi<(-2)] <- c(-2)
linshi<-as.data.frame(linshi)

# linshi<-linshi[1:500,c(1:100,1400:1500,2600:2700)]
# library(pheatmap)
out<- pheatmap(linshi,fontsize=6,cutree_col = 5,
               color  = colorRampPalette(c(blue,white,red))(100),
               annotation_col = annotation_col,
               annotation_colors = ann_colors,
               clustering_method = "ward.D",
               # border_color = "grey60",
               cluster_cols = F, cluster_rows = T,
               show_rownames = T, show_colnames = T,
               border=FALSE
)



meth_down<-t(meth_down)
meth_down<-as.data.frame(meth_down)
label<-label[rownames(meth_down),]
meth_down<-cbind(label[,2],meth_down)
colnames(meth_down)[1]<-c("label")

data<-aggregate(apply(meth_down[,2:882],2,as.numeric),by=list(meth_down$label),mean)
rownames(data)<-data[,1]
data<-data[,-1]

jieguo<-t(data)
jieguo<-as.data.frame(jieguo)

jieguo = jieguo[apply(jieguo, 1, function(x) sd(x)!=0),]   
linshi <- apply(jieguo,1,scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(jieguo)
hist(linshi)
linshi[linshi>2] <- 2
linshi[linshi<(-2)] <- c(-2)
linshi<-as.data.frame(linshi)


out<- pheatmap(linshi,fontsize=6,cutree_col = 5,
               color  = colorRampPalette(c(blue,white,red))(100),
               annotation_col = annotation_col,
               annotation_colors = ann_colors,
               clustering_method = "ward.D",
               # border_color = "grey60",
               cluster_cols = F, cluster_rows = T,
               show_rownames = T, show_colnames = T,
               border=FALSE
)


####
library(GSEABase)
 GSE190113_mRNA_l 
data<-read.table("GSE190113_mRNA_m.txt")
label<-read.table("GSE190113_mRNA_l.txt")

#TEX-specific 
gs <- read.csv("耗竭分类pathway.csv",stringsAsFactors = F,row.names = 1)
gs<-gs[1:4,1:200]
gs <- t(gs)
gs <- as.data.frame(gs)
gs <- as.list(gs)
gs <- lapply(gs, function(x) x[which(x != "")])

data<-t(data)
data <- as.matrix(data)   ###
library(GSVA)
jieguo <- gsva(data,gs, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)

write.table(jieguo,"GSE190113-通路富集得分.txt",quote = F,sep = "\t")



label<-label[order(label$DNN_label),]
jieguo<-jieguo[,rownames(label)]

annotation_col = data.frame( cluster = as.factor(label[,3]))
rownames(annotation_col) = rownames(label)
ann_colors = list(cluster = c('1' = "#0B5F99", 
                              "2" = "#86B5CE",
                              "3" = "#F2D5C1",
                              "4" = "#DE9B7D", 
                              "5" = "#B22023"))


blue <- rgb(63,117,179,maxColorValue = 255)
red <- rgb(214,47,36,maxColorValue = 255)
white <- rgb(255,255,255,maxColorValue = 255)


linshi <- apply(jieguo,1,scale)



linshi <- t(linshi)
colnames(linshi) <- colnames(jieguo)
hist(linshi)
linshi[linshi>2] <- 2                
linshi[linshi<(-2)] <- c(-2)


library(pheatmap)
out<- pheatmap(linshi,fontsize=6,cutree_col = 5,
               color  = colorRampPalette(c(blue,white,red))(100),
               annotation_col = annotation_col,
               annotation_colors = ann_colors,
               clustering_method = "ward.D",
               # border_color = "grey60",
               cluster_cols = F, cluster_rows = T,
               show_rownames = T, show_colnames = T,
               border=FALSE
)




library(GSEABase)
data<-read.table("GSE190826_mRNA_m.txt")
label<-read.table("GSE190826_mRNA_l.txt")

#TEX-specific 
gs <- read.csv("D:/Check/zicheng_check/check/data/Fig1_data/耗竭分类pathway.csv",stringsAsFactors = F,row.names = 1)
gs<-gs[1:4,1:200]
gs <- t(gs)
gs <- as.data.frame(gs)
gs <- as.list(gs)
gs <- lapply(gs, function(x) x[which(x != "")])

data<-t(data)
data <- as.matrix(data)   
library(GSVA)
jieguo <- gsva(data,gs, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)




label<-label[order(label$DNN_label),]
jieguo<-jieguo[,rownames(label)]

annotation_col = data.frame( cluster = as.factor(label[,3]))
rownames(annotation_col) = rownames(label)
ann_colors = list(cluster = c('1' = "#0B5F99", 
                              "2" = "#86B5CE",
                              "3" = "#F2D5C1",
                              "4" = "#DE9B7D", 
                              "5" = "#B22023"))


red <- rgb(255,0,0,maxColorValue = 255)
blue <- rgb(0,0,255,maxColorValue = 255)
white <- rgb(255,255,255,maxColorValue = 255)
##############最终颜色
blue <- rgb(63,117,179,maxColorValue = 255)
red <- rgb(214,47,36,maxColorValue = 255)
white <- rgb(255,255,255,maxColorValue = 255)

linshi <- apply(jieguo,1,scale)



linshi <- t(linshi)
colnames(linshi) <- colnames(jieguo)
hist(linshi)
linshi[linshi>2] <- 2               
linshi[linshi<(-2)] <- c(-2)


library(pheatmap)
out<- pheatmap(linshi,fontsize=6,cutree_col = 5,
               color  = colorRampPalette(c(blue,white,red))(100),
               annotation_col = annotation_col,
               annotation_colors = ann_colors,
               clustering_method = "ward.D",
               # border_color = "grey60",
               cluster_cols = F, cluster_rows = T,
               show_rownames = T, show_colnames = T,
               border=FALSE
)


### GSEA 
#pathway-activation
a <- GSEABase::getGmt("D:/MSigDB-pathway/Go/GOBP_T_CELL_ACTIVATION.v2023.1.Hs.gmt")
T_acti<-a[["GOBP_T_CELL_ACTIVATION"]]@geneIds
T_acti<-as.data.frame(T_acti)
T_acti[,2]<-"GOBP_T_CELL_ACTIVATION"
T_acti<-T_acti[,c(2,1)]
colnames(T_acti)<-c("term","gene")

#pathway-PROLIFERATION 
a <- GSEABase::getGmt("D:/MSigDB-pathway/Go/GOBP_T_CELL_PROLIFERATION.v2023.1.Hs.gmt")
T_pro<-a[["GOBP_T_CELL_PROLIFERATION"]]@geneIds
T_pro<-as.data.frame(T_pro)
T_pro[,2]<-"GOBP_T_CELL_PROLIFERATION"
T_pro<-T_pro[,c(2,1)]
colnames(T_pro)<-c("term","gene")


#pathway-DIFFERENTIATION 
a <- GSEABase::getGmt("D:/MSigDB-pathway/Go/GOBP_T_CELL_DIFFERENTIATION.v2023.1.Hs.gmt")
T_diff<-a[["GOBP_T_CELL_DIFFERENTIATION"]]@geneIds
T_diff<-as.data.frame(T_diff)
T_diff[,2]<-"GOBP_T_CELL_DIFFERENTIATION"
T_diff<-T_diff[,c(2,1)]
colnames(T_diff)<-c("term","gene")


data_113<-cbind(gene_113,fc_113)
colnames(data_113)<-c("gene","logFC")

data_826<-cbind(gene_826,fc_826)
colnames(data_826)<-c("gene","logFC")


data_113 <- data_113[order(data_113$logFC,decreasing = T),]
genelist <- data_113$logFC
names(genelist) <- data_113$gene


data_826 <- data_826[order(data_826$logFC,decreasing = T),]
genelist <- data_826$logFC
names(genelist) <- data_826$gene


gsea.re1<- GSEA(genelist,  
                TERM2GENE = T_pro,  
                pvalueCutoff = 0.05, 
                pAdjustMethod = 'BH') 
g1<-as.data.frame(gsea.re1)

library(ggsci)
col_gsea1<-pal_simpsons()(16)
gseaplot2(gsea.re1,geneSetID = rownames(g1)[1:1],
          title = "T cell activation",
          color = col_gsea1[1:1],
          base_size = 14,
          rel_heights = c(1, 0.2, 0.4),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line"
)

gseaplot2(gsea.re1,geneSetID = rownames(g1)[2:2],
          title = "T cell differentiation",
          color = col_gsea1[2:2],
          base_size = 14,
          rel_heights = c(1, 0.2, 0.4),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line"
)

gseaplot2(gsea.re1,geneSetID = rownames(g1)[3:3],
          title = "T cell proliferation",
          color = col_gsea1[3:3],
          base_size = 14,
          rel_heights = c(1, 0.2, 0.4),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line"
)



T_cell<-rbind(T_acti,T_pro,T_diff)
gsea.re1<- GSEA(genelist,  
                TERM2GENE = T_cell,  
                pvalueCutoff = 1,  
                pAdjustMethod = 'BH',
                nPermSimple = 10000
)  
g1<-as.data.frame(gsea.re1)
num2=3
gseaplot2(gsea.re1,geneSetID = rownames(g1)[1:num2],
          title = "",
          color = col_gsea1[1:num2],
          base_size = 14,
          rel_heights = c(1, 0.2, 0.4),
          subplots = 1:3,
          pvalue_table = FALSE,
          ES_geom = "line"
)