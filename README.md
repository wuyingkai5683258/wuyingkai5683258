#55235+77298+12021（96）

library(tidyverse)
library(GEOquery)
library(tidyverse)
library(GEOquery)
library(limma) 
library(affy)
library(stringr)
gset = getGEO('GSE55235', destdir=".", AnnotGPL = T, getGPL = T)
class(gset)
gset[[1]]
gset2 = getGEO('GSE12021', destdir=".", AnnotGPL = T, getGPL = T)
class(gset2)
gset2[[1]]
gset3 = getGEO('GSE77298', destdir=".", AnnotGPL = T, getGPL = T)
class(gset3)
gset3[[1]]

plf1<-gset[[1]]@annotation
plf2<-gset2[[1]]@annotation
plf3<-gset3[[1]]@annotation
GPL_data<- getGEO(filename ="GPL96.annot.gz", AnnotGPL = T)
GPL_data_11 <- Table(GPL_data)
GPL_data1<- getGEO(filename ="GPL96.annot.gz", AnnotGPL = T)
GPL_data_22 <- Table(GPL_data1)
GPL_data2<- getGEO(filename ="GPL570.annot.gz", AnnotGPL = T)
GPL_data_33 <- Table(GPL_data2)
raw_geneid<-fData(gset[[1]])
raw_geneid<-raw_geneid[,c(1,3)]
raw_geneid2<-fData(gset2[[1]])
raw_geneid2<-raw_geneid2[,c(1,3)]
raw_geneid3<-fData(gset3[[1]])
raw_geneid3<-raw_geneid3[,c(1,3)]




exp <- exprs(gset[[1]])
probe_name<-rownames(exp)
exp2 <- exprs(gset2[[1]])
probe_name2<-rownames(exp2)
exp3 <- exprs(gset3[[1]])
probe_name3<-rownames(exp3)


raw_geneid$`Gene symbol`<-data.frame(sapply(raw_geneid$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]
raw_geneid2$`Gene symbol`<-data.frame(sapply(raw_geneid2$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]
raw_geneid3$`Gene symbol`<-data.frame(sapply(raw_geneid3$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]

loc<-match(GPL_data_11[,1],probe_name)
probe_exp<-exp[loc,]
raw_geneid<-(as.matrix(raw_geneid[,"Gene symbol"]))
index<-which(!is.na(raw_geneid))
geneid<-raw_geneid[index]
exp_matrix<-probe_exp[index,]
geneidfactor<-factor(geneid)
gene_exp_matrix<-apply(exp_matrix,2,function(x) tapply(x,geneidfactor,mean))
rownames(gene_exp_matrix)<-levels(geneidfactor)
gene_exp_matrix<-gene_exp_matrix[,-(11:20)]

loc2<-match(GPL_data_22[,1],probe_name2)
probe_exp2<-exp2[loc2,]
loc2
raw_geneid2<-(as.matrix(raw_geneid2[,"Gene symbol"]))
index2<-which(!is.na(raw_geneid2))
geneid2<-raw_geneid2[index2]
exp_matrix2<-probe_exp2[index2,]
geneidfactor2<-factor(geneid2)
gene_exp_matrix2<-apply(exp_matrix2,2,function(x) tapply(x,geneidfactor2,mean))
rownames(gene_exp_matrix2)<-levels(geneidfactor2)
gene_exp_matrix2<-gene_exp_matrix2[,-c(6,7,13,16:21,25)]

loc3<-match(GPL_data_33[,1],probe_name3)
loc3
probe_name3
probe_exp3<-exp3[loc3,]
raw_geneid3<-(as.matrix(raw_geneid3[,"Gene symbol"]))
index3<-which(!is.na(raw_geneid3))
geneid3<-raw_geneid3[index3]
exp_matrix3<-probe_exp3[index3,]
geneidfactor3<-factor(geneid3)
gene_exp_matrix3<-apply(exp_matrix3,2,function(x) tapply(x,geneidfactor3,mean))
rownames(gene_exp_matrix3)<-levels(geneidfactor3)



geo_exp_1=as.data.frame(gene_exp_matrix)
geo_exp_2=as.data.frame(gene_exp_matrix2)
geo_exp_3=as.data.frame(gene_exp_matrix3)
sameSample=intersect(rownames(geo_exp_1), rownames(geo_exp_2))
sameSample=as.data.frame(sameSample)
sameSample=intersect(rownames(geo_exp_3), sameSample$sameSample)
gene_exp1=geo_exp_1[sameSample,,drop=F]
gene_exp2=geo_exp_2[sameSample,,drop=F]
gene_exp3=geo_exp_3[sameSample,,drop=F]
bindgeo=cbind(gene_exp1,gene_exp2,gene_exp3)

pdata <- pData(gset[[1]])
pdata <- pdata[-(11:20),]
pdata2 <- pData(gset2[[1]])
pdata2<-pdata2[-c(6,7,13,16:21,25),]
pdata3 <- pData(gset3[[1]])
group_list <- ifelse(str_detect(pdata$title,"healthy joint"), "NC",
                     "RA")

group_list
group_list = factor(group_list,
                    levels = c("NC","RA"))
group_list
pdata$group=group_list

group_list2 <- ifelse(str_detect(pdata2$title,"Normal"), "NC",
                     "RA")

group_list2
group_list2 = factor(group_list2,
                     levels = c("NC","RA"))
group_list2
pdata2$group=group_list2

group_list3 <- ifelse(str_detect(pdata3$title,"HC"), "NC",
                      "RA")

group_list3
group_list3 = factor(group_list3,
                     levels = c("NC","RA"))
group_list3
pdata3$group=group_list3

group1<-(as.matrix(pdata[,"group"]))
row.names(group1)=rownames(pdata)
colnames(group1)="group"
group2<-(as.matrix(pdata2[,"group"]))
row.names(group2)=rownames(pdata2)
colnames(group2)="group"
group3<-(as.matrix(pdata3[,"group"]))
row.names(group3)=rownames(pdata3)
colnames(group3)="group"
talgroup=as.data.frame(rbind(group1,group2,group3))
talgroup_list=factor(talgroup$group,levels = c("RA","NC"))
write.csv(talgroup,file = "group.csv")

boxplot(bindgeo,outline=T, notch=T,col=talgroup_list, las=2)
dev.off()
bindgeo_normal=normalizeBetweenArrays(bindgeo)
boxplot(bindgeo_normal,outline=T, notch=T,col=talgroup_list, las=2)
range(bindgeo_normal)
bindgeo_normal <- log2(bindgeo_normal+1)
bindgeo_normal <-as.data.frame(bindgeo_normal)
bindgeo_normal=na.omit(bindgeo_normal)
write.csv(bindgeo_normal,file = "bindgeo_exp.csv")
range(bindgeo_normal)
dev.off()
#chayi fenxi
expFile="bindgeo_exp.csv"          #?????ļ?
#??ȡ?????ļ????????????ļ?????
expr_data=read.csv(expFile,sep=",",header=T,row.names = 1,check.names=F)
group<-read.csv("group.csv",header = T,row.names = 1,sep = ",")
design <- model.matrix(~0+factor(group$group))
colnames(design) <- levels(factor(group$group))
rownames(design) <- colnames(expr_data)
contrast.matrix <- makeContrasts(RA-NC,levels = design) 
fit <- lmFit(expr_data,design) #非线性最小二乘法
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)
DEG$regulate <- ifelse(DEG$adj.P.Val> 0.05, "unchanged",
                       ifelse(DEG$logFC > 1, "up-regulated",
                              ifelse(DEG$logFC < -1, "down-regulated", "unchanged")))
write.table(DEG,"DEG.csv",row.names=T,col.names=T,sep=",")
job <- "test"
write.table(table(DEG$regulate),file = paste0(job,"_","DEG_result_1_005.txt"),
            sep = "\t",quote = F,row.names = T,col.names = T)
write.table(data.frame(gene_symbol=rownames(DEG),DEG),file = paste0(job,"_","DEG_result.txt"),
            sep = "\t",quote = F,row.names = F,col.names = T)
DE_1_0.05 <- DEG[DEG$adj.P.Val<0.05&abs(DEG$logFC)>1,]
upGene_1_0.05 <- DE_1_0.05[DE_1_0.05$regulate == "up-regulated",]
downGene_1_0.05 <- DE_1_0.05[DE_1_0.05$regulate == "down-regulated",]
write.csv(upGene_1_0.05,paste0(job,"_","upGene_1_005.csv"))
write.csv(downGene_1_0.05,paste0(job,"_","downGene_1_005.csv"))
#CHAYI FUJI
library(dplyr) # 用于数据处理
library(gt) # 制作表格
#记得去GEO_all里面第一列加上Gene
data<-read.csv("DEG.csv", header = TRUE, sep = ",", dec = ".", quote = "\"", fill = TRUE, comment.char = "")
#data<-read.csv("GEO_diff.csv", header = TRUE, sep = ",", dec = ".", quote = "\"", fill = TRUE, comment.char = "")
data<- data%>%
  
  mutate(expression = case_when(logFC>= 1 &adj.P.Val< 0.05 ~ "up", # 上调
                                
                                logFC<= -1 &adj.P.Val< 0.05 ~ "down", # 下调
                                
                                TRUE ~ "Unchanged")) # 不变
data<-na.omit(data)
library(ggplot2)
library(ggrepel)
volc_plot<- ggplot(data, aes(logFC, -log10(P.Value))) +xlim(-9,9)+ylim(0,18)+# 将p值进行-log10转化
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1.2,1.2), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = expression),
             size =2.5, 
             alpha =0.5) +
  theme_bw(base_size = 12)+
  ggsci::scale_color_jama() +
  theme(panel.grid = element_blank(),
        legend.position = 'right')

volc_plot
top_20 <-bind_rows(
  
  data%>%
    
    filter(expression=='up') %>%
    
    arrange(P.Value, desc(abs(logFC))) %>%
    
    head(10),
  
  data%>%
    
    filter(expression=='down') %>%
    
    arrange(P.Value, desc(abs(logFC))) %>%
    
    head(10)
  
)
write.table(top_20,file="top_20.txt",sep="\t",quote=F,col.names=T)
options(ggrepel.max.overlaps = Inf)#不让他有重叠
top_20 %>% gt()

volc_plot2 <- volc_plot+
  
  geom_label_repel(data = top_20,
                   
                   aes(logFC, -log10(P.Value), label = gene),
                   
                   size = 2.5)

volc_plot2
ggsave(filename="volplot20.png",
       width=10,
       height=8,
       units="in",
       dpi=300)



#前**位基因热图
library(ggplot2)
library(pheatmap)
data<-read.csv("DEG.csv", header = TRUE, sep = ",", dec = ".", quote = "\"", fill = TRUE, comment.char = "")
#data<-read.csv("GEO_diff.csv", header = TRUE, sep = ",", dec = ".", quote = "\"", fill = TRUE, comment.char = "")
data<- data%>%
  
  mutate(expression = case_when(logFC>= 1 &adj.P.Val< 0.05 ~ "up", # 上调
                                
                                logFC<= -1 &adj.P.Val< 0.05 ~ "down", # 下调
                                
                                TRUE ~ "Unchanged")) # 不变
top_100 <-bind_rows(
  
  data%>%
    
    filter(expression=='up') %>%
    
    arrange(P.Value, desc(abs(logFC))) %>%
    
    head(50),
  
  data%>%
    
    filter(expression=='down') %>%
    
    arrange(P.Value, desc(abs(logFC))) %>%
    
    head(50)
  
)
write.table(top_100,file="top_100.csv",sep=",",quote=F,col.names=T,row.names = T)
data=read.csv("top_100.csv",header=TRUE,row.names=1,check.names = FALSE)  
expr_data<-read.csv("bindgeo_exp.csv", header = TRUE, sep = ",", dec = ".", quote = "\"", fill = TRUE, comment.char = "",row.names =1)
colnames(data)
rownames(data)
DEG_gene_expr100 <- expr_data[rownames(data),]

pheatmap(DEG_gene_expr100)
p <- pheatmap(DEG_gene_expr100,scale="row",
              border="white",  # 设置边框为白色
              cluster_cols = T, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #去掉横、纵坐标id
              show_colnames = F,
              annotation_col = group)
p
group <-read.csv(file="group.csv",sep=",",header=TRUE,row.names=1,check.names=FALSE)
dev.off()
#
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stats)
library(data.table)
library(dplyr)
DEG_data1 <- read.csv("DEGS.csv",sep=",",header = T)
##差异基因筛选，这里选择abs(logFC) > 2,FDR < 0.05的值作为差异基因。也可以根据自己的情况进行设定

#记得把GEO_diff.csv.里面的id更改位SYMBOL,Gene名转化为GeneID.并在前面加上一列
gene.df <- bitr(DEG_data1$gene, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                toType = c("ENTREZID", "SYMBOL"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = org.Hs.eg.db) #Orgdb是指对应的注释包是哪个
write.table(gene.df,file="gene.df.csv",sep=",",quote=F,col.names=T)
gene <- gene.df$ENTREZID
GO_all <-ego_ALL <- enrichGO(gene = gene,#我们上面定义了
                             OrgDb=org.Hs.eg.db,
                             keyType = "ENTREZID",
                             ont = "ALL",#富集的GO类型
                             pAdjustMethod = "BH",#这个不用管，一般都用的BH
                             minGSSize = 1,
                             pvalueCutoff = 0.05,#P值可以取0.05
                             qvalueCutoff = 0.05,
                             readable = TRUE)
GO_all
GO_result <- data.frame(GO_all)
go_enrichment_pathway <- GO_result %>% group_by(ONTOLOGY) %>% top_n(n = 5, wt = -p.adjust)
write.table(go_enrichment_pathway,file="go_enrichment_pathway.csv",sep=",",quote=F,row.names = F) 
PLOTYELLOW<-ggplot(go_enrichment_pathway, aes(x=reorder(Description, Count), y=Count)) +
  geom_point(aes(size=Count,color=-log10(p.adjust))) +
  scale_size_continuous(range=c(2, 6)) +
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  coord_flip() +  #让柱状图变为纵向
  theme_minimal() +
  scale_color_gradient(low = "pink",high ="red")+
  labs(color=expression(-log10(p.adjust),size="Count"), 
       x="Gene Ratio",y="Gene_Number",title="GO Enrichment")+
  theme_bw()
ggsave(filename="GO.png",
       width=10,
       height=5,
       units="in",
       dpi=300)
ggsave("GO of DEGS.pdf",width=10,height=4)# 保存为pdf格式
data <- read.csv("GO_enrichment_pathway.csv",header=TRUE)
GO_term_order=factor(as.integer(rownames(data)),labels=data$Description)
ggplot(data=data, aes(x=GO_term_order,y=Count, fill=ONTOLOGY)) + geom_bar(stat="identity", width=0.8) + coord_flip() +theme_bw()



ek <- enrichKEGG(gene = gene.df$ENTREZID, #需要分析的基因的EntrezID
                 organism = "hsa",  #人
                 pvalueCutoff =0.05, #设置pvalue界值
                 qvalueCutoff = 0.05) #设置qvalue界值(FDR校正后的p值）
write.table(ek,file="eRICHk.txt",sep="\t",quote=F,row.names = F)   
ek2 = setReadable(ek, #前面分析的结果
                  OrgDb = "org.Hs.eg.db", #人类数据库
                  keyType = "ENTREZID") #要转换的基因类型
head(ek2@result$geneID)
write.table(ek2,file="ek2.txt",sep="\t",quote=F,row.names = F)  
library("enrichplot")
pdf(file="ek_barplot-YELLOW.pdf",width = 7,height = 3) 
barplot(ek, x = "GeneRatio", color = "p.adjust", #默认参数
        showCategory =5) #只显示前10
dev.off()
library(tidyverse)
ek.rt = read.table("eRICHk.txt",header=TRUE,sep="\t",quote = "")  #读取第1部分enrichKEGG分析输出的文件ek。
ek.rt <- separate(data=ek.rt, col=GeneRatio, into = c("GR1", "GR2"), sep = "/") #劈分GeneRatio为2列（GR1、GR2）
ek.rt <- separate(data=ek.rt, col=BgRatio, into = c("BR1", "BR2"), sep = "/") #劈分BgRatio为2列（BR1、BR2）
ek.rt <- mutate(ek.rt, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) #计算Enrichment Factor
ek.rt10 <- ek.rt %>% filter(row_number() >= 1,row_number() <= 10)
p <- ggplot(ek.rt10,aes(enrichment_factor, fct_reorder(factor(Description), enrichment_factor))) + 
  geom_point(aes(size=Count,color=-1*log10(pvalue))) +
  scale_color_gradient(low="pink",high ="red") + 
  labs(color=expression(-log[10](p_value)),size="Count",
       x="Enrichment Factor",y="KEGG term",title="KEGG enrichment") + 
  theme_bw()
p
ggsave("er.rt10_KEGG.pdf",width=8,height=4)# 保存为pdf格式
ggsave("er.rt10_KEGG.png",width=8,height=4,dpi=300)# 保存为png格式
#DO富集分析
rt=read.table("gene.df.csv",sep=",",check.names=F,header=T,row.names = 1) 
library("clusterProfiler")

library("org.Hs.eg.db")

library("enrichplot")

library("ggplot2")

##在ID转换前，将要分析的基因名字+logFC保存为symbol.text文件。

##ID转换
genes=as.vector(rt[,1])

entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)

entrezIDs <- as.character(entrezIDs)

out=cbind(rt,entrezID=entrezIDs)

write.table(out,file="id.txt",sep="\t",quote=F,row.names=F)

##读取ID转换后文件

rt=read.table("id.txt",sep="\t",header=T,check.names=F)

rt=rt[is.na(rt[,"entrezID"])==F,]
gene=rt$entrezID

library(DOSE) #加载需要使用的包

erich.do<-DOSE::enrichDO(gene=gene,
                         
                         ont = "DO",
                         
                         pvalueCutoff = 0.05,
                         
                         qvalueCutoff = 0.05,
                         
                         readable = T) #DO代码，其中P值和Q值可以自己设置，自己按照需求设置，我这里设为0.05
DO<-barplot(erich.do)
ggsave("DO.png",width=8,height=4,dpi=300)
#gsva
library(ggplot2)
library(ComplexHeatmap) 
library(clusterProfiler) 
library(GSVA) 
library(GSEABase)
library(dplyr) 
library(data.table) 
library(tidyverse)
COAD <- read.csv("bindgeo_exp.csv", check.names = FALSE, header = TRUE, row.names = 1)
#COAD <- log2(COAD+1)
gene_set <- getGmt("c2.cp.kegg.v7.4.symbols.gmt")
group_list <-  read.csv("group.csv", check.names = FALSE, header = TRUE, row.names = 1)
table(group_list)
annotation <- data.frame(group_list)
rownames(annotation) <- colnames(COAD)
head(annotation)

gsva_result<- gsva(as.matrix(COAD), gene_set, method = "gsva",min.sz=1,
                   max.sz=Inf,kcdf="Gaussian")
library(dendextend)
library(circlize)
library(RColorBrewer)
colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
values <- seq(-0.8, 0.8, length.out = 101)[-101]
col_fun = colorRamp2(values, colors)
top_annotation<-HeatmapAnnotation(df=annotation,col=list(group=c("Normal"="blue","OA"="red")))
Heatmap(gsva_result, name = "GSVA", col = col_fun,cluster_rows = T,cluster_columns = F,show_row_names = T,
        show_column_names = F,column_split = annotation$group,)
annotation <- annotation[order(annotation$group == "RA", decreasing = TRUE), , drop = FALSE]


exp="RA"
ctr="NC"
design <- model.matrix(~0+factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_result)
contrast.matrix <- makeContrasts(contrasts=paste0(exp,'-',ctr),  #"exp/ctrl"
                                 levels = design)

fit1 <- lmFit(gsva_result,design)                 #拟合模型
fit2 <- contrasts.fit(fit1, contrast.matrix) #统计检验
efit <- eBayes(fit2)                         #修正

summary(decideTests(efit,lfc=1, p.value=0.05)) #统计查看差异结果
tempOutput <- topTable(efit, coef=paste0(exp,'-',ctr), n=Inf)
degs <- na.omit(tempOutput) 
keep <- rownames(degs[degs$adj.P.Val<0.25& degs$P.Value<0.05, ])
length(keep)
dat <- gsva_result[keep[1:20],] #选取前50进行展示
pheatmap(dat,cluster_rows = T,cluster_cols = F,
         color=colorRampPalette(c("navy","white","firebrick3"))(100),
         show_colnames = F,border_color = NA,scale = "row",show_rownames =T,
         annotation_col = group_list,fontsize=5.0,name = "GSVA")


#gsea
#oa
library(clusterProfiler)
library(org.Hs.eg.db)
KEGG_database <- 'hsa'
####4??GSEA????####
library(tidyverse)
GSEoa= read.table("DEG.csv",sep = ",",header = T)
library(dplyr)
GSEoa<- GSEoa %>% filter(abs(logFC) > 0) %>% filter(adj.P.Val < 0.05)
names(GSEoa)[1] <- "SYMBOL"
GSEoa1 = GSEoa[,c(1,2)]


genename <- GSEoa1$SYMBOL
En_id <- mget(genename, 
              org.Hs.egSYMBOL2EG, 
              ifnotfound=NA)
En_id <- as.character(En_id)
a=data.frame(ENTREZID=En_id,logFC=GSEoa1$logFC)

a$logFC<-sort(a$logFC,decreasing = T)

geneList = a[,2]
names(geneList) = as.character(a[,1])
geneList
library(stats)
GSEA_KEGG <- gseKEGG(geneList, organism = 'hsa', pvalueCutoff = 0.05)
ridgeplot(GSEA_KEGG) 
gseaplot2(GSEA_KEGG,1)
gseaplot2(GSEA_KEGG, title = GSEA_KEGG$Description[1], geneSetID = 1)
gseaplot2(GSEA_KEGG,1:5)
gseaplot2(GSEA_KEGG,geneSetID = 1,
          title = "",
          color = "green",
          base_size = 12,
          rel_heights = c(1.5, 0.5, 1),
          subplots = 1:3,
          pvalue_table = T,
          ES_geom = "line"
)

#ssgsea
library(GSVA)
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(pheatmap)
GSE39582_series_tumor<-read.csv("bindgeo_exp.csv",header = T,row.names = 1,sep = ",")
library(tidyverse)
cellMarker1 <- read.delim("1.txt", header = F, sep = "\t") #
cellMarker1 <- cellMarker1 %>% column_to_rownames("V1") %>% t()
a <- cellMarker1
a[1:5,1:5]
a <- a[1:nrow(a), ]
set <- colnames(a)
geneSet <- list()
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  geneSet[[i]] <-x
}
library(GSVA)
GSE39582_series_tumor=as.matrix(GSE39582_series_tumor)
gsva_matrix1 <- gsva(GSE39582_series_tumor, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
res <- gsva_matrix1
pheatmap(res, show_colnames = F)
group_list <- read.csv("group.csv",header = T,row.names = 1,sep = ",")
annotation <- data.frame(group_list)
rownames(annotation) <- colnames(res)
resm<-res
for(i in colnames(res)){
  resm[,i]<-(res[,i]-min(res[,1]))/max(res[,i]-min(res[,i]))
}
pheatmap(res,
         show_colnames = F,
         annotation_col = annotation,
         fontsize = 10, 
)
dt<-resm%>%t()%>%as.data.frame()%>%
  rownames_to_column("sample")%>%
  gather(key = cell_type,
         value = value,-sample)
head(dt)
dtt<-dt%>%
  group_by(sample)%>%
  mutate(proportion=round(value/sum(value),3))
head(dtt)
dtt$cell_type<-factor(dtt$cell_type,levels = unique(rownames(res)))
mytheme<-theme(axis.title = element_text(size=12),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               plot.title = element_text(size = 13,
                                         hjust = 0.5,
                                         face = "bold"),
               legend.text = element_text(size = 10),
               legend.position = "bottom")
library(paletteer)
d_palettes<-palettes_d_names
col<-paletteer_d("khroma::smoothrainbow",n=28)
p<-ggplot(dtt,
          aes(x=cell_type,y=proportion,fill=cell_type))+
  geom_boxplot(color="black",alpha=0.6,outlier.shape = 21,outlier.size = 1.2)+
  scale_fill_manual(values = col)+
  labs(x="cell type",y="proportion")+
  theme_bw()+mytheme
p
write.table(dtt,file="dtt.csv",sep=",",quote=T,row.names=T)
rest<-t(res)
write.table(rest,file="rest.csv",sep=",",quote=T,row.names=T)
#然后加上分组和列
dtt<-read.csv("dtt.csv",header = T,row.names = 1,sep = ",")
colnames(dtt)
library(tidyverse)
p1<-ggplot(dtt,
           aes(x=cell_type,y=proportion,fill=group))+
  geom_boxplot(color="black",alpha=0.6,outlier.shape = 21,outlier.size = 1.2)+
  scale_fill_manual(values = c("#4979b6","#d9352a"))+
  labs(x="cell type",y="proportion")+
  theme_bw()+mytheme+theme(axis.text.x = element_text(angle=45))
p1
library(ggsignif)
pvalues <- sapply(dtt$cell_type, function(x) {
  res <- wilcox.test(as.numeric(proportion) ~ group, data = subset(dtt, cell_type == x)) #两组，wilcox.test或t.test；多组，kruskal.test或aov(one-way ANOVA test)
  res$p.value
})
pv <- data.frame(gene = dtt$cell_type, pvalue = pvalues)

pv$sigcode <- cut(pv$pvalue, c(0,0.0001, 0.001, 0.01, 0.05, 1), 
                  labels=c('****','***', '**', '*', 'ns'))

p.box <- ggplot(dtt, aes(x=cell_type, y=proportion, color=group, fill=group)) +
  geom_boxplot(alpha = .5) + #半透明
  theme_classic() + #或theme_bw()
  scale_fill_brewer(palette = "Set1") + #按类填充颜色
  scale_color_brewer(palette = "Set1") + #按类给边框着色
  
  theme(axis.text.x = element_text(colour="black", size = 11,
                                   #名太挤，旋转45度
                                   angle = 90, hjust = .5, vjust = .5)) +
  geom_text(aes(x=gene, y=max(dtt$proportion) * 1.1,
                label = pv$sigcode),
            data=pv, 
            inherit.aes=F)
p.box
#wgcna
library(WGCNA)
exprMat <-"bindgeo_exp.csv"
dataExpr<-read.table(exprMat, sep=',', row.names=1, header=T,
                     quote="", comment="", check.names=F)
str(dataExpr)
my_mad <- function(x){mad(x,na.rm = TRUE)} #mad时将缺失值去掉
m.mad <- apply(dataExpr,2,my_mad)

#筛选时不设置mad最小值，直接使用前75%的基因或者探针。
datExpr0 <- dataExpr[,which(m.mad > quantile(m.mad, probs=seq(0, 1, 0.25))[2])]
write.table(datExpr0,"datExpr0.csv",row.names=T,col.names=T,sep=",")
#dataExprVar
## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(datExpr0))
datExpr0<-dataExpr
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "1.sampleClustering.pdf", width = 15, height = 8)
par(cex = 0.6)
par(mar = c(0,6,6,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,
     cex.axis = 1.5, cex.main = 2)
### Plot a line to show the cut
#abline(h = 7000, col = "red")##剪切高度不确定，故无红线
dev.off()
pdf("2_sample clutering_delete_outliers.pdf", width = 6, height = 8)
par(cex = 0.6)
par(mar = c(0,6,6,0))
cutHeight = 130  ### cut高度
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2, 
     cex.axis = 1.5,cex.main = 2) +
  abline(h = 130, col = "red")    ###'@"h = 1500"参数为你需要过滤的参数高度
dev.off()
#记得把原来表达矩阵的去掉的cut的样本去掉
#dataExpr<-dataExpr[-(3),]
#datExpr0<-dataExpr
traitData = read.csv("sampleinfo.csv",row.names=1)
# ##'@导入txt格式
# traitData = read.table("TraitData.txt",row.names=1,header=T,comment.char = "",check.names=F)
head(traitData)
allTraits = traitData
dim(allTraits)
names(allTraits)
# 形成一个类似于表达数据的数据框架
fpkmSamples = rownames(datExpr0)
traitSamples =rownames(allTraits)
traitRows = match(fpkmSamples, traitSamples)
datTraits = allTraits[traitRows,]
datTraits=allTraits
rownames(datTraits)
#记得把样本的该去掉的去掉
#datTraits=datTraits[-(4:5),]
rownames(datTraits)
datTraits
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
pdf(file="3_Sample_dendrogram_and_trait_heatmap.pdf",width=8 ,height= 6)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap",cex.colorLabels = 1.5, cex.dendroLabels = 1, cex.rowText = 2)
dev.off()
enableWGCNAThreads()
# 设置soft-thresholding powers的数量
powers = c(1:30)
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf(file="4_软阈值选择.pdf",width=12, height = 8)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
sft = pickSoftThreshold(datExpr0, powerVector = powers)
sft$powerEstimate
softPower =8
adjacency = adjacency(datExpr0, power = softPower)
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file="4_Gene clustering on TOM-based dissimilarity.pdf",width=6,height=4)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()
minModuleSize =30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
pdf(file="5_Dynamic Tree Cut.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
#sizeGrWindow(7, 6)
pdf(file="6_Clustering of module eigengenes.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25######剪切高度可修改
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs
table(mergedColors)

#输出所有modules
color<-unique(mergedColors)
for (i  in 1:length(color)) {
  y=t(assign(paste(color[i],"expr",sep = "."),datExpr0[mergedColors==color[i]]))
  write.csv(y,paste('6',color[i],"csv",sep = "."),quote = F)
}
#save.image(file = "module_splitted.RData")

##'@输出merge模块图形
pdf(file="7_merged dynamic.pdf", width = 8, height = 8)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
nGenes = ncol(datExpr0);
nSamples = nrow(datExpr0);
# 重新计算带有颜色标签的模块
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# 通过相关值对每个关联进行颜色编码
pdf(file="5.Module-trait relationships.pdf",width=6,height=6)
# 展示模块与表型数据的相关系数和 P值
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
trait_a = as.data.frame(datTraits$RA)
datTraits
names(trait_a) = "trait_a"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, trait_a, use = "p"));#和性状的关联
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(trait_a), sep="");
names(GSPvalue) = paste("p.GS.", names(trait_a), sep="")
module = "salmon"
column = match(module, modNames);
moduleGenes = moduleColors==module;
pdf(file="7.Module membership vs. gene significance.pdf",width=6,height=6);
par(mfrow = c(1,1));
verboseScatterplot(
  abs(geneModuleMembership[moduleGenes, column]),
  abs(geneTraitSignificance[moduleGenes, 1]),
  xlab = paste("Module Membership in", module, "module"),
  ylab = "Gene significance for proliferating",
  main = paste("Module membership vs. gene significance\n"),
  abline = TRUE,
  pch = 21,
  cex.main = 1.2,
  cex.lab = 1.2,
  cex.axis = 1.2,
  col = "black",
  bg = module
)
dev.off()
module = "salmon"
column = match(module, modNames)
moduleGenes = moduleColors==module
blue_module<-as.data.frame(dimnames(data.frame(datExpr0))[[2]][moduleGenes])
names(blue_module)="genename"
MM<-abs(geneModuleMembership[moduleGenes,column])
GS<-abs(geneTraitSignificance[moduleGenes, 1])
blue_MMGS<-as.data.frame(cbind(MM,GS))
rownames(blue_MMGS)=blue_module$genename
hub_b<-abs(blue_MMGS$MM)>0.8&abs(blue_MMGS$GS)>0.5
table(hub_b)
blue_hub_b<-subset(blue_MMGS, abs(blue_MMGS$MM)>0.8)
write.csv(blue_hub_b, "hubgene_MMGS_123.csv")
#vein
library (VennDiagram) 
dat <- read.csv('ragene.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE, check.names = FALSE,row.names = 1)

#以2个分组为例
#指定统计的分组列，并设置作图颜色、字体样式等
venn_list <- list(WGCNA_hubgenes = dat$ASPM, DEGs = dat$PLIN1)
fill.col<-c("#0055aa","#c40003")
venn.diagram(venn_list, 
             scaled=F,
             filename = "veinDEGS.tiff",
             cex=2,
             cat.cex=1.5,
             cat.pos=c(3.-3),
             fill=fill.col)
inter <- get.venn.partitions(venn_list)
#LASSO
VEIN=read.table("VEIN.csv",sep=",",check.names=F,header=T,row.names = 1) 
colnames(VEIN)
rownames(VEIN)
exp=read.table("bindgeo_exp.csv",sep=",",check.names=F,header=T,row.names = 1)
DEG_gene_exprvein1<-exp[rownames(VEIN),]
DEG_gene_exprvein1<-t(DEG_gene_exprvein1)
write.table(DEG_gene_exprvein1,file="DEG_gene_exprvein1.csv",sep=",",quote=T,col.names=T)
data<- read.csv("DEG_gene_exprvein1.csv", header = T, sep=",",row.names =1)
library(tidyverse)

library(glmnet)
library(sigFeature)

library(e1071)

library(caret)

library(randomForest) 
train <- read.csv("DEG_gene_exprvein1.csv",row.names = 1,
                  
                  as.is = F) #后面svmRFE函数要求group的类型为factor

dim(train)

train[1:4,1:4]

# 转为lasso需要的格式

x <- as.matrix(train[,-15])

(y <- ifelse(train$group == "NC", 0,1)) #把分组信息换成01

fit = glmnet(x, y, family = "binomial", alpha = 1, lambda = NULL)

# 绘制LASSO回归曲线图

pdf("1A_lasso.pdf", width = 30, height = 15)

plot(fit, xvar = "dev", label = TRUE)

dev.off()

#绘制LASSO回归10折交叉验证图

cvfit = cv.glmnet(x, y,
                  
                  nfold=10, #例文描述：10-fold cross-validation
                  
                  family = "binomial", type.measure = "class")

pdf("2cvfit.pdf")

plot(cvfit)

dev.off()

#查看最佳lambda

cvfit$lambda.min

# 获取LASSO选出来的特征

myCoefs <- coef(cvfit, s="lambda.min")

lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]

(lasso_fea <- lasso_fea[-1])

# 把lasso找到的特征保存到文件

write.csv(lasso_fea,"3feature_lasso.csv")
rm(list=ls())
#SVM-REF算法输入数据
predcm<-read.csv("DEG_gene_exprvein1.csv",row.names = 1,
                 
                 as.is = F)
control <- rfeControl(functions = caretFuncs, method = "cv", number = 5)  #cv 交叉验证次数5

# 执行SVM-RFE算法
results <- rfe(predcm[,1:14],   #1-8列为预测变量
               predcm[,15],    #9列为诊断变量
               sizes = c(1:14),  
               rfeControl = control,
               method = "svmRadial") # method = "svmRadial" specifies that the SVM model should use a radial kernel
# 结果分析
print(results)
# 列出选择的变量集
predictors(results)
# 绘制结果
pdf("6B_svm-accuracy.pdf",width = 5,height = 5)
plot(results, type=c("g", "o"))
dev.off()
p<-predictors(results)
write.csv(p,"feature_svm.csv")
#随街森林
VEIN=read.table("VEIN.csv",sep=",",check.names=F,header=T,row.names = 1) 
exp=read.table("bindgeo_exp.csv",sep=",",check.names=F,header=T,row.names = 1)
DEG_gene_exprvein1<-exp[rownames(VEIN),]
DEG_gene_exprvein1<-t(DEG_gene_exprvein1)
write.table(DEG_gene_exprvein1,file="DEG_gene_exprvein1.csv",sep=",",quote=T,col.names=T)
data<- read.csv("DEG_gene_exprvein1.csv", header = T, sep=",",row.names =1)
data$group<- ifelse(data$group == "NC", 0,1)
set.seed(11)
library(caret)
index <- createDataPartition(
  data$group,
  p = 0.7,
  list = FALSE
)
train <- data[index, ]
test <- data[-index, ]
set.seed(123)
x<-train %>% dplyr::select(-group)
y<-factor(train$group)
rf<-randomForest(x, y,data=train,importance=TRUE)
plot(margin(rf, train$group), main = '观测值被判断正确的概率图')
test_predict <- predict(rf, test)
compare_test <- table(test$group, test_predict, dnn = c('Actual', 'Predicted'))
compare_test
importance_otu <- data.frame(importance(rf), check.names = FALSE)
head(importance_otu)
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
#根据表格输出前30变量
rf
varImpPlot(rf, n.var = min(30, nrow(rf$importance)), main = 'Top 30 - variable importance')
write.table(importance_otu, 'importance_otu.txt', sep = '\t', col.names = NA, quote = FALSE)
set.seed(123)
otu_train.cv <- replicate(5, rfcv(train[-ncol(train)], train$group, cv.fold = 10,step = 1.5), simplify = FALSE)
otu_train.cv
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)
head(otu_train.cv.mean, 37)
library(ggplot2)
library(splines)  #用于在 geom_smooth() 中添加拟合线，或者使用 geom_line() 替代 geom_smooth() 绘制普通折线
p <- ggplot(otu_train.cv, aes(otus, value)) +
  geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = '',x = 'Number of OTUs', y = 'Cross-validation error')
p
p + geom_vline(xintercept =7)
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
importance_otu.select <- importance_otu[1:7, ]
importance_otu.select
write.table(importance_otu.select, 'importance_otu.select.csv', sep = ',', col.names = NA, quote = FALSE)
#
data<- read.csv("DEG_gene_exprvein1.csv", header = T, sep=",",row.names =1)
data$group<- ifelse(data$group == "NC", 0,1)
set.seed(11)
library(caret)
index <- createDataPartition(
  data$group,
  p = 0.7,
  list = FALSE
)
train <- data[index, ]
test <- data[-index, ]
set.seed(123)
x<-train %>% dplyr::select(-group)
y<-factor(train$group)
rf<-randomForest(x, y,data=train,importance=TRUE)
plot(margin(rf, train$group), main = '观测值被判断正确的概率图')
test_predict <- predict(rf, test)
compare_test <- table(test$group, test_predict, dnn = c('Actual', 'Predicted'))
hist(treesize(rf))
plot(rf,main="ERROR &TREES")
legend("top",
       legend=colnames(rf$err.rate),
       lty = 1:3,
       col = 1:3,
       horiz = T)

library (VennDiagram) 
dat <- read.csv('123 hub.csv', header = TRUE, sep = ',', stringsAsFactors = FALSE, check.names = FALSE,row.names = 1)

#以2个分组为例
#指定统计的分组列，并设置作图颜色、字体样式等
venn_list <- list(RandomForest_hubgenes = dat$RRM2, LASSO_hubgenes= dat$DLGAP5,SVM_hubgenes=dat$TOP2A)
fill.col<-c("#0055aa","#c40003",'green')
venn.diagram(venn_list, 
             scaled=F,
             filename = "veinDEGS.tiff",
             cex=2,
             cat.cex=1.5,
             cat.pos=c(3.-3),
             fill=fill.col)
inter <- get.venn.partitions(venn_list)
#RRM2 DLGAP5 KIF11
library(tidyverse)
library(GEOquery)
library(tidyverse)
library(GEOquery)
library(limma) 
library(affy)
library(stringr)
gset <- getGEO('GSE55457', destdir=".",
               AnnotGPL = T,     ## 注释文件
               getGPL = T)  
exp<-exprs(gset[[1]])
exp<-log2(exp+1)
cli<-pData(gset[[1]])	## 获取临床信息
GPL<-fData(gset[[1]])
gpl<-GPL[,c(1,3)]        
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"///"))[1]),stringsAsFactors=F)[,1]
gpl$`Gene symbol`<-data.frame(sapply(gpl$`Gene symbol`,function(x)unlist(strsplit(x,"-"))[1]),stringsAsFactors=F)[,1]
exp<-as.data.frame(exp)
exp$ID<-rownames(exp)
exp_symbol<-merge(exp,gpl,by="ID")
exp_symbol<-na.omit(exp_symbol)
#写入table后去更改table，将symbol放在第一位，然后分列-，去除-，然后筛选重复项后直接去除，得出的表达矩阵就是需要的矩阵，不要求平均值
table(duplicated(exp_symbol$`Gene symbol`))
exp_unique<-avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$`Gene symbol`)
write.table(exp_unique,"exp55457.csv",sep=",")
library(rms)
mydata<-read.csv("exp yuan.csv", header = TRUE, sep = ",")
str(mydata)
mydata$group<-ifelse(mydata$group =="RA",1,0)  
ddist <- datadist(mydata)
options(datadist='ddist')
fit1<-lrm(group~RRM2+DLGAP5+KIF11,data=mydata,x=TRUE,y=TRUE)
Nomo_LR<-nomogram(fit1,
                  fun=plogis,
                  lp=F,
                  fun.at = c(0.001,0.01,0.1,0.5,0.9,0.99),
                  funlabel = "Risk of RA")

plot(Nomo_LR)
cal1 <- calibrate(fit1, method='boot', B=1000) 
plot(cal1,
     xlim = c(0,1),
     xlab = "Predicted Probability",
     ylab = "Observed Probability",
     legend = FALSE,
     subtitles = FALSE)
abline(0,1,col = "black",lty = 2,lwd = 2)
lines(cal1[,c("predy","calibrated.orig")], type = "l",lwd = 2,col="red",pch =16)
lines(cal1[,c("predy","calibrated.corrected")], type = "l",lwd = 2,col="green",pch =16)
legend(0.55,0.35,
       c("Apparent","Ideal","Bias-corrected"),
       lty = c(2,1,1),
       lwd = c(2,1,1),
       col = c("black","red","green"),
       bty = "n") # "o"为加边框

library(rmda)
library(ggplot2)
library(rms)
library(caret)
complex<-decision_curve(group ~RRM2+DLGAP5+KIF11, 
                        data = mydata,family = binomial(link ='logit'),
                        thresholds = seq(0,1, by = 0.01),
                        confidence.intervals= 0.95)
plot_decision_curve(complex,
                    curve.names=c("Model"),#图例上每条曲线的名字
                    cost.benefit.axis =FALSE,# cost.benefit.axis是另外附加的一条横坐标轴，损失收益比，默认值是TRUE
                    col= c("blue","red","black"),#曲线颜色
                    confidence.intervals=FALSE,# 设置是否画出曲线的置信区间
                    standardize = FALSE)# 设置是否对净受益率（NB）使用患病率进行校正

library(multipleROC)
rm(list = ls())
df <- as.data.frame(mydata)
p1 <- multipleROC(group~RRM2,data=df)
#p2 <- multipleROC(group~ERBB4,data=df)
p3 <- multipleROC(group~DLGAP5,data=df)
p4<-multipleROC(group~KIF11,data=df)
p5<-multipleROC(group~LMCD1,data=df)
#p6<-multipleROC(group~CSF3R,data=df)
plot_ROC(list(p1,p3,p4),
         show.points = T, 
         show.eta = F, 
         show.sens = F, 
         show.AUC = T, 
         facet = F )
library(pROC)
library(ggplot2)
pred_f_training<-predict(fit1,mydata)
#下方参数中Death需改为你的研究的结局变量名

modelroc <- roc(mydata$group,pred_f_training)
#绘制ROC
auc(modelroc)# AUC
ci(modelroc) #AUC95%CI
# Area under the curve: 0.7851
# 95% CI: 0.7335-0.8368 (DeLong)
plot(modelroc,col="red",#颜色
     legacy.axes=T,#y轴格式更改
     print.auc=TRUE,#显示AUC面积
     print.thres=TRUE,#添加截点和95%CI
     grid=c(0.2,0.2),grid.col=c("blue","yellow"))#网格线设置

library(rms)
mydata<-read.csv("exp yuan.csv", header = TRUE, sep = ",")
str(mydata)
mydata$group<-ifelse(mydata$group =="RA",1,0)  
ddist <- datadist(mydata)
options(datadist='ddist')
fit1<-lrm(group~RRM2+DLGAP5+KIF11,data=mydata,x=TRUE,y=TRUE)
Nomo_LR<-nomogram(fit1,
                  fun=plogis,
                  lp=F,
                  fun.at = c(0.001,0.01,0.1,0.3,0.5,0.7,0.9,0.99),
                  funlabel = "Risk of RA")

plot(Nomo_LR)
cal1 <- calibrate(fit1, method='boot', B=1000) 
plot(cal1,
     xlim = c(0,1),
     xlab = "Predicted Probability",
     ylab = "Observed Probability",
     legend = FALSE,
     subtitles = FALSE)
abline(0,1,col = "black",lty = 2,lwd = 2)
lines(cal1[,c("predy","calibrated.orig")], type = "l",lwd = 2,col="red",pch =16)
lines(cal1[,c("predy","calibrated.corrected")], type = "l",lwd = 2,col="green",pch =16)
legend(0.55,0.35,
       c("Apparent","Ideal","Bias-corrected"),
       lty = c(2,1,1),
       lwd = c(2,1,1),
       col = c("black","red","green"),
       bty = "n") # "o"为加边框

library(rmda)
library(ggplot2)
library(rms)
library(caret)
complex<-decision_curve(group ~RRM2+DLGAP5+KIF11, 
                        data = mydata,family = binomial(link ='logit'),
                        thresholds = seq(0,1, by = 0.01),
                        confidence.intervals= 0.95)
plot_decision_curve(complex,
                    curve.names=c("Model"),#图例上每条曲线的名字
                    cost.benefit.axis =FALSE,# cost.benefit.axis是另外附加的一条横坐标轴，损失收益比，默认值是TRUE
                    col= c("blue","red","black"),#曲线颜色
                    confidence.intervals=FALSE,# 设置是否画出曲线的置信区间
                    standardize = FALSE)# 设置是否对净受益率（NB）使用患病率进行校正

library(multipleROC)
rm(list = ls())
df <- as.data.frame(mydata)
p1 <- multipleROC(group~RRM2,data=df)
#p2 <- multipleROC(group~ERBB4,data=df)
p3 <- multipleROC(group~DLGAP5,data=df)
p4<-multipleROC(group~KIF11,data=df)
p5<-multipleROC(group~LMCD1,data=df)
#p6<-multipleROC(group~CSF3R,data=df)
plot_ROC(list(p1,p3,p4),
         show.points = T, 
         show.eta = F, 
         show.sens = F, 
         show.AUC = T, 
         facet = F )
library(pROC)
library(ggplot2)
pred_f_training<-predict(fit1,mydata)
#下方参数中Death需改为你的研究的结局变量名

modelroc <- roc(mydata$group,pred_f_training)
#绘制ROC
auc(modelroc)# AUC
ci(modelroc) #AUC95%CI
# Area under the curve: 0.7851
# 95% CI: 0.7335-0.8368 (DeLong)
plot(modelroc,col="red",#颜色
     legacy.axes=T,#y轴格式更改
     print.auc=TRUE,#显示AUC面积
     print.thres=TRUE,#添加截点和95%CI
     grid=c(0.2,0.2),grid.col=c("blue","yellow"))#网格线设置

#mianyi ixangguan
library(ggcorrplot)
library(tidyr)
exp<-read.csv("bindgeo_exp.csv",row.names = 1)
re<-read.csv("res.csv",row.names = 1)
mygene <- c("RRM2","DLGAP5",'KIF11')  #定义你的目的基因
nc = t(rbind(re,exp[mygene,]))  ;#将你的目的基因匹配到表达矩阵---行名匹配--注意大小写
m = rcorr(nc)$r[1:nrow(re),(ncol(nc)-length(mygene)+1):ncol(nc)]

##计算p值
p = rcorr(nc)$P[1:nrow(re),(ncol(nc)-length(mygene)+1):ncol(nc)]
head(p)


tmp <- matrix(case_when(as.vector(p) < 0.01 ~ "**",
                        as.vector(p) < 0.05 ~ "*",
                        TRUE ~ ""), nrow = nrow(p))

##绘制热图
library(pheatmap)
p1 <- pheatmap(t(m),
               display_numbers =t(tmp),
               angle_col =45,
               color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
               border_color = "white",
               cellwidth = 20, 
               cellheight = 20,
               width = 7, 
               height=9.1,
               treeheight_col = 0,
               treeheight_row = 0)
##ENRICHR
library(enrichR)
dbs <- listEnrichrDbs() ###列出164个库
dbs[1:4,1:4]
#   geneCoverage genesPerTerm               libraryName
# 1        13362          275       Genome_Browser_PWMs
# 2        27884         1284  TRANSFAC_and_JASPAR_PWMs
# 3         6002           77 Transcription_Factor_PPIs
# 4        47172         1370                 ChEA_2013
#                                                       link
# 1 http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/
# 2                 http://jaspar.genereg.net/html/DOWNLOAD/
# 3                                                         
# 4           http://amp.pharm.mssm.edu/lib/cheadownload.jsp
###从中选择你要富集的库
dbs$libraryName ###查看库名
dbs <- c("DsigDB")###这里我选择GO库的3个process
library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
symbol=read.csv("all.diffsig.csv",header = T,row.names = 1)
symbol=symbol$gene
#BiocManager::install("org.Hs.eg.db")
###id转换
df <- bitr(unique(symbol$V1), fromType = "ENSEMBL",
           toType = c("SYMBOL","ENTREZID"),
           OrgDb = org.Hs.eg.db)
enrichr<- enrichr(symbol, dbs)
write.csv(enrichr,file="enrichr.csv")
###有点久，因为要联网
