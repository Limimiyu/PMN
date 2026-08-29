library(Seurat)
library(monocle)
library(harmony)
library(celldex)
library(dplyr)
library(RColorBrewer)
library(ggrepel)
library(scales)
library(ggpubr)
library(msigdbr)
library(AUCell)
library(irGSEA)
library(UCell)
library(GSVA)
library(pheatmap)
library(cowplot)
library(CellChat)
library(slingshot)
library(aplot)
library(forcats)
library(tidydr)
library(magrittr)
library(RColorBrewer)
library(NMF)
library(ggalluvial)
library(psych)
library(qgraph)
library(igraph)
library(clusterProfiler)
library(SCENIC)
library(GENIE3)
library(RcisTarget)
library(SCopeLoomR)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)
library(ggrepel)
library(tibble)
library(tidydr)
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
colors = c( '#1F78B4', '#FFFFB3', '#E31A1C')
cols <- hue_pal()((20))
##拟时颜色
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
#提取dimplot中图的颜色
col <- c("#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#D2B48C","#7B68EE","#DC143C","#800080","#A0522D","#D2691E","#87CEEB","#40E0D0")

p1 <- DimPlot(wt, group.by = "celltype_assign")
x<-ggplot_build(p1)
info = data.frame(colour = x$data[[1]]$colour, group = x$data[[1]]$group)
info <- unique((arrange(info, group)))
cols <- as.character(info$colour)


assays <- c("HC1","HC2", "WXT","MN1","MN2", "MN3","MN4","LC_220402")
dir=paste0("data5/practice/HE/",assays)
samples_names=c("HC1","HC2", "WXT","MN1","MN2", "MN3","MN4","LC_220402")
scRNAlist <- list()##创建空list
pbmc.function <- function(sample_dir, sample_name, group) {
  counts <- Read10X(data.dir = sample_dir)
  scRNA <- CreateSeuratObject(counts, project = sample_name, min.cells = 3, min.features = 200)
  scRNA <- RenameCells(scRNA, add.cell.id = sample_name)
  scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
  scRNA <- subset(scRNA, nFeature_RNA > 600 & nFeature_RNA < 5000 & percent.mt < 20)
  scRNA$group <- group
  return(scRNA)
}

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0 )

for (i in 1:length(dir)) {
  group <- ifelse(i <= 3, "healthy", "PMN")
  scRNAlist[[i]] <- pbmc.function(dir[i], samples_names[i], group)
}


pbmc.list <- lapply(X = scRNAlist, FUN = function(x) {
       x <- NormalizeData(x)
       x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
   })

features <- SelectIntegrationFeatures(object.list = pbmc.list)

anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = features)
pbmc <- IntegrateData(anchorset = anchors)

DefaultAssay(pbmc) <- "integrated"
pbmc <- ScaleData(pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc, npcs = 30, verbose = FALSE)
ElbowPlot(pbmc,ndims = 20)
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:15)
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.4)
DefaultAssay(pbmc) <- "RNA"
DimPlot(pbmc,reduction = "umap",label = TRUE)
DimPlot(pbmc,reduction = "umap",group.by = "orig.ident")
DefaultAssay(pbmc) <- "RNA"
FeaturePlot(pbmc,features = c("CD3E","CD4","CD8A","CD8B","CD79A","NKG7","CD14","FCGR3A","LYZ","FCER1A","LILRA4","PPBPp","MKI67","DEFA4"))
marker.plot <- c("CD3E","CD4","CD8A","CD8B","CD79A","NKG7","GNLY","FCN1","FCGR3A","CLEC10A","FCER1A","IRF8","LILRA4","MKI67","PPBP","DEFA4")
pbmc@active.ident <-factor(pbmc@active.ident,levels = c("0","1","3","4","8","19","5","9","12","6","7","2","13","16","10","14","17","15","11","18"))
DotPlot(pbmc,features = marker.plot,cols = c("lightgrey","red"))+theme(axis.text.x = element_text(face = "bold",angle = 30,hjust = 1))+coord_flip()
pbmc <- RenameIdents(pbmc,"0"="CD4+ T cell","1"="CD4+ T cell","3"="CD8+ T cell", "4"="γδT","8"="CD8+ T cell","19"="CD4+ T cell","5"="B","9"="B","12"= "B","6"="NK", "7"="NK","2"="Mono","13"="Mono","16"="Mono","10"="CD14- Mono","14"="cDC","17"="pDC","15"="MKI67+ cell","11"="PLT","18"="Neu")




##验证批次效应
embedding <- Embeddings(pbmc, reduction = "pca")  
embedding_pca <- Embeddings(pbmc, reduction = "pca")[, 1:10]
lisi_scores <- compute_lisi(
  X = embedding_pca,                
  meta_data = pbmc@meta.data,      
  label_colnames = c("orig.ident", "celltype"), 
  perplexity = 30                  
)
pbmc@meta.data$iLISI <- lisi_scores$bLISI
pbmc@meta.data$cLISI <- lisi_scores$cell_type_LISI
meta <- pbmc@meta.data
p1 <- ggplot(meta,aes(x=orig.ident,y=iLISI,fill=orig.ident))+geom_boxplot(outlier.shape = NA)+theme_classic()+theme(axis.text = element_text(size = 8,colour = "black"))+NoLegend()
p2 <- ggplot(meta,aes(x=orig.ident,y=cLISI,fill=orig.ident))+geom_boxplot(outlier.shape = NA)+theme_classic()+theme(axis.text = element_text(size = 8,colour = "black"))+NoLegend()
(p1+p2)




##基因表达热图
data.features <- FetchData(Mono, cells = colnames(Mono), vars = marker.plot, slot = "data")
p <- DotPlot(B,features = marker.plot)
dotplot_data<-p$data
heatmap_data<-dotplot_data%>%
  select(features.plot,id,avg.exp.scaled)%>%
  pivot_wider(names_from=features.plot,values_from=avg.exp.scaled)
heatmap_data=column_to_rownames(heatmap_data,var="id")

# Generate the heatmap
pheatmap(heatmap_data,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),          show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_rows = F,cluster_cols = F,
         scale = "column" ,cellheight = 10,cellwidth = 10,fontsize = 10, border_color = "white")




table(CD8$orig.ident)#查看各组细胞数
prop.table(table(Idents(CD8)))
table(Idents(CD8), CD8$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(CD8), CD8$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- data.frame(Cellratio)
Cellratio$Var2 <- as.character(Cellratio$Var2)
Cellratio$Var2[Cellratio$Var2 == "LC_220402"] <- "MN5"
Cellratio$Var2[Cellratio$Var2 == "WXT"] <- "HC3"
Cellratio$Var2 <- as.factor(Cellratio$Var2)
colnames(Cellratio)[1] <- "Celltype"

ggplot(Cellratio)+
  geom_ <- (aes(x = factor(Var2,levels=c("HC1","HC2","HC3","MN1","MN2","MN3","MN4","MN5")),y=Freq,fill=Celltype),stat = "identity",width = 0.8,show.legend = TRUE,position = "stack",color="black")+
  theme_classic()+
  xlab(label = "sample")+
  ylab(label = "Freq")+
  scale_fill_manual(values = cols)+
  guides(guide_legend(override.aes = list(size=50)))+
  theme(legend.text = element_text(face = "bold",colour = "#222222",size = 18))+
  theme(axis.title.x = element_text(size = 15,face = "bold",color="#222222" ))+
  theme(axis.text.x = element_text(size = 15,face = "bold",color = "#222222",angle = 45,hjust = 1,vjust = 1))+
  theme(axis.title.y = element_text(face = "bold",size = 15,color ="#222222" ))+
  theme(axis.text.y = element_text(face = "bold",size = 18,colour = "#222222"))+
  labs(title = "γδT")+
  theme(title  = element_text(size = 15,face = "bold",colour = "#222222"))


library(reshape2)
cellper <- dcast(Cellratio,Var2~Celltype, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]


###添加分组信息
sample <- c("HC1","HC2","MN1","MN2","MN3","MN4","MN5","HC3")
group <- c("healthy","healthy","PMN","PMN","PMN","PMN","PMN","healthy")
samples <- data.frame(sample, group)#创建数据框

rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']
cellper$group <- samples[rownames(cellper),'group']

###作图展示
pplist = list()
sce_groups = c("Naive CD8+ T",'GATA3+ Tcm','MAIT cell','Tem-1','CMC1+ Tem-2','Exhaust CD8+ T','CD8+ Temra-1','CD8+ Temra-2','NEAT1+ CD8+ T')
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)

for(i in 1:length(sce_groups)){
  cellper_  = cellper %>% select(one_of(c('sample','group',sce_groups[i])))#选择一组数据
  colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
  cellper_$percent = as.numeric(cellper_$percent)#数值型数据
  
  if (i==6) {
    pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
      geom_boxplot(shape = 21,aes(fill=group),width = 0.25) + 
      stat_summary(fun=mean, geom="point", color="grey60") +
      theme_cowplot() +
      theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 15),
            legend.title = element_text(size = 10),plot.title = element_text(size = 15,face = 'plain')) + 
      labs(title = sce_groups[i],y='Percentage')+
      theme(axis.text.y = element_text(size = 15))+
      theme(axis.title.y = element_text(size = 15,face = "bold"))+
      theme(axis.title.x = element_text(size = 15,face = "bold"))+
      theme(axis.text.x  = element_text(size = 15))
  }
  else {
    pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
      geom_boxplot(shape = 21,aes(fill=group),width = 0.25,show.legend = FALSE) + 
      stat_summary(fun=mean, geom="point", color="grey60") +
      theme_cowplot() +
      theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),plot.title = element_text(size = 15,face = 'plain')) + 
      labs(title =  sce_groups[i],y='Percentage')+
      theme(axis.text.y = element_text(size = 15))+
      theme(axis.title.y = element_text(size = 15,face = "bold"))+
      theme(axis.title.x = element_text(size = 15,face = "bold"))+
      theme(axis.text.x  = element_text(size = 15))
  }
  
  
  
  
  ###组间t检验分析
  my_comparisons <- list( c("PMN", "healthy") )
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 4,method = "t.test",)
  pplist[[i]] = pp1
}

plot_grid(pplist[[1]],
          pplist[[2]],
          pplist[[3]],
          pplist[[4]],
          pplist[[5]],
          pplist[[6]],
          pplist[[7]],
          pplist[[8]],
          pplist[[9]])
 


cellper.1 <- melt(cellper,id=c("sample","group"))
colnames(cellper.1)[3] <- "celltype"
colnames(cellper.1)[4] <- "percentage"
ggplot(cellper.1, aes(x = celltype, y = percentage, fill = group)) +
  geom_boxplot(shape = 21, aes(fill = group), width = 0.5, position = position_dodge(0.8)) +
  scale_y_sqrt() +
  theme_cowplot() +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15), 
        legend.text = element_text(size = 13), legend.title = element_text(size = 15),
        plot.title = element_text(size = 15, face = 'plain'),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1, vjust = 1),legend.position = "top") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", hide.ns = TRUE, size = 8) +
  scale_fill_manual(values = c('#6699CC', "#CC3333"))


##95%confidence
median_diff <- function(data, indices) {
  d <- data[indices, ]
  pmn_med <- median(d$percentage[d$group == "PMN"])
  hc_med <- median(d$percentage[d$group == "healthy"])
  return(pmn_med - hc_med)
}
Mono.result <- cellper.1 %>%
  group_by(celltype) %>%
  do({
    dat <- .
    diff_obs <- median(dat$percentage[dat$group == "PMN"]) -
      median(dat$percentage[dat$group == "healthy"])
    boot_res <- boot(dat, statistic = median_diff, R = 1000)
    ci <- boot.ci(boot_res, type = "perc")$percent[4:5]  # 95% 百分位 CI
    data.frame(
      celltype = unique(dat$celltype),
      median_diff = diff_obs,
      CI_lower = ci[1],
      CI_upper = ci[2],
      p_value = t.test(percentage ~ group, data = dat)$p.value
    )
  }) %>%
  ungroup()
print(Mono.result)


##gsva
expr <- as.data.frame(Mono@assays$RNA@data)
meta <- Mono@meta.data[,c("orig.ident","seurat_clusters")] 
m_df = msigdbr(species = "Homo sapiens", category = "C5",subcategory = "BP")
msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)
expr=as.matrix(expr) 
kegg <- gsva(expr, msigdbr_list, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
meta <- meta %>%arrange(meta$seurat_clusters)
data <- kegg[,rownames(meta)]
group <- factor(meta[,"seurat_clusters"],ordered = F)
data1 <-NULL
for(i in 0:(length(unique(group))-1)){
  ind <-which(group==i)
  dat <- apply(data[,ind], 1, mean)
  data1 <-cbind(data1,dat)
}
colnames(data1) <-c("C0","C1","C2","C3","C4","C5","C6","C7")
result<- t(scale(t(data1)))
result[result>2]=2
which(rownames(result) == "BIOCARTA_ASBCELL_PATHWAY")  
annotation_col=data.frame(celltype=factor(c(rep("Naive CD4+ T",1),rep("CD4+ Tcm-1",1),rep("CD4+ Tcm-2",1),rep("CD4+ Tem",1),rep("Th17",1),rep("Treg",1),rep("RPS26+ CD4+ T",1),rep("ISG+ CD4+ T",1),rep("NEAT1+ CD4+ T",1)),levels = c("Naive CD4+ T","CD4+ Tcm-1","CD4+ Tcm-2","CD4+ Tem","Th17","Treg","RPS26+ CD4+ T","ISG+ CD4+ T","NEAT1+ CD4+ T")))
row.names(annotation_col) <- colnames(result)
pheatmap(result[c(1705,1726,1738,1812),],
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T,
         show_colnames = T,
         color =colorRampPalette(c("blue", "white","red"))(100),
         cellwidth = 10, cellheight = 15,
         fontsize = 10 )




##AUcell评分
human_KEGG = msigdbr(species = "Homo sapiens",category = "C5",subcategory = "BP") %>%dplyr::select(gs_name,gene_symbol)
human_KEGG_Set = human_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
sce2 <- AddModuleScore(CD4,features = TCR_features,ctrl = 100,name = "TCR pathway")
colnames(sce2@meta.data)[9] <- "Chemokine_score"

cells_rankings <- AUCell_buildRankings(sce2@assays$RNA@data,plotStats = FALSE)/cells_rankings <- AUCell_buildRankings(sce2@assays$RNA@data,splitByBlocks=TRUE)
cells_AUC <- AUCell_calcAUC(human_KEGG_Set, cells_rankings,aucMaxRank=nrow(cells_rankings)*0.1)
geneSet <- "KEGG_CHEMOKINE_SIGNALING_PATHWAY"
AUCell_auc <- as.numeric(getAUC(cells_AUC)[geneSet, ])
sce2$AUCell <- AUCell_auc

names(MHC.regulation) <- "MHC.regulation"
cells_AUC.Antigen <- AUCell_calcAUC(MHC.regulation, cells_rankings,aucMaxRank=nrow(cells_rankings)*0.1)
AUCell_auc.Antigen <- as.numeric(getAUC(cells_AUC.Antigen))
B$Antigen <- AUCell_auc.Antigen


##基因评分
sce_groups <- levels(gene_expression$celltype)
pplist <- vector("list", length(sce_groups))
# 对每个sce_groups中的元素进行迭代
for(i in 1:length(sce_groups)){
  gene.expression <- subset(gene_expression, gene_expression$celltype == sce_groups[i])
  gene.expression <- gene.expression[, 2:65]
  gene.expression <- as.data.frame(gene.expression)
  # 将数据框架的每一列转换为数值类型
  for(j in 1:length(colnames(gene.expression))){
    gene.expression[, j] <- as.numeric(gene.expression[, j])
  }
  # 添加一行为列的平均数
  gene.expression[nrow(gene.expression) + 1, ] <- colMeans(gene.expression)
  gene.expression.sub <- gene.expression[nrow(gene.expression), ]
  # 设置行名
  rownames(gene.expression.sub) <- sce_groups[i]
  # 转置数据框架
  gene.expression.sub <- t(gene.expression.sub)
  # 将结果赋值给pplist的第i个元素
  pplist[[i]] <- gene.expression.sub
}
gene.expression <- cbind(pplist[[1]],pplist[[2]],pplist[[3]],pplist[[4]],pplist[[5]],pplist[[6]],pplist[[7]],pplist[[8]])
gene.expression <- t(gene.expression)
pheatmap(gene.expression,gaps_col  = c(14,18,36,44,64),cluster_rows = FALSE,cluster_cols = FALSE, color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100), border_color = "white",scale = "column",cellwidth = 10, cellheight = 10,fontsize = 8,fontface_col= "italic")




##基因评分条形图
CD8.score.sub <- meta.CD8.1%>%group_by(celltype)%>%dplyr::summarise_each(funs = mean)
CD8.score.sub$celltype <- factor(CD8.score.sub$celltype,levels = c("NEAT1+ CD8+ T","CD8+ Temra-2","CD8+ Temra-1","Exhaust CD8+ T","CMC1+ Tem-2","Tem-1","MAIT cell","GATA3+ Tcm","Naive CD8+ T"))
CD8.score.sub$exhaust.score <- as.numeric(CD8.score.sub$exhaust.score)
CD8.score.sub$above_average <- CD8.score.sub$exhaust > median(CD8.score.sub$exhaust)
ggplot(Mono.score.sub) +
  geom_bar(aes(x = chemokine, y = celltype, fill = above_average), stat = "identity", width = 0.8,show.legend = FALSE) +scale_fill_manual(values = c("FALSE" = "grey", "TRUE" = "blue")) + geom_vline(aes(xintercept = median(chemokine)), linetype = "dashed") + theme_classic()+theme(axis.text.y = element_blank(),axis.title.y = element_blank())+theme(axis.title.x = element_text(size = 20))








###monocle3
Idents(Mono) <- Mono$site
Mono <- subset(Mono,idents="tumor")
Idents(Mono) <- Mono$cell.type
mat <- GetAssayData(Mono, slot = "data")
cellInfo <- Mono@meta.data
geneInfo <- data.frame(gene_short_name = rownames(mat), row.names = rownames(mat))
cds <- new_cell_data_set(expression_data = mat,
                         cell_metadata = cellInfo,
                         gene_metadata = geneInfo)
cds <- preprocess_cds(cds, num_dim = 50)
reducedDim(cds, "UMAP") <- Embeddings(Mono, "umap")
cds = cluster_cells(cds, cluster_method = 'leiden',k = 20,reduction_method = "UMAP",resolution = 0.00001)
cds = learn_graph(cds, use_partition=T, verbose=T, learn_graph_control=list(
  minimal_branch_len=30
))

start = c("Mono_CD14_DLEU2")
closest_vertex = cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex = as.matrix(closest_vertex[colnames(cds), ])
root_pr_nodes = igraph::V(principal_graph(cds)[["UMAP"]])$name

flag = closest_vertex[as.character(colData(cds)$cell.type) %in% start,]
flag = as.numeric(names(which.max(table( flag ))))
root_pr_nodes = root_pr_nodes[flag]
cds = order_cells(cds, root_pr_nodes=root_pr_nodes)
p = plot_cells(cds,
               color_cells_by = "pseudotime",
               label_cell_groups=F,
               label_groups_by_cluster=F,
               label_roots=F,
               label_leaves=F,
               label_branch_points=F,
               cell_size=0.5,
               group_label_size=7,
               rasterize=F)




###slingshot应用
sce <- as.SingleCellExperiment(CD8,assay = "RNA")
sce_slingshot1 <- slingshot(sce,reducedDim="UMAP",clusterLabels = CD8$celltype,approx_points=150,start.clus="Naive CD8+ T")
SlingshotDataSet(sce_slingshot1)
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

# 使用自定义颜色向量直接替换
cell_colors <- colour[sce_slingshot1$celltype]
celltype_levels <- levels(sce_slingshot1$celltype)
named_colors <- setNames(colour, celltype_levels)
cell_colors <- named_colors[as.character(sce_slingshot1$celltype)]
plot(reducedDims(sce_slingshot1)$UMAP, col = cell_colors, pch=16, asp = 1, cex = 0.8)
lines(SlingshotDataSet(sce_slingshot1), lwd=2, col='black')

celltype_label <- CD8@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(celltype = CD8@meta.data$celltype) %>%
  group_by(celltype) %>%
  summarise(UMAP1 = median(UMAP_1),
            UMAP2 = median(UMAP_2))

for (i in 1:8) {
  text(celltype_label$celltype[i], x=celltype_label$UMAP1[i]-1, y=celltype_label$UMAP2[i])
}


##粗略点
lin1 <- getLineages(sce_slingshot1, 
                    clusterLabels = "celltype")
plot(reducedDims(sce_slingshot1)$UMAP,col = brewer.pal(10,'Paired')[sce$Cluster],pch=16,asp=1)
lines(SlingshotDataSet(lin1), lwd=2,col = 'black',type = 'lineages')         












###细胞通讯  cellchat
##第一：分组跑cellchat
healthy.object <- subset(pbmc,group=="healthy")
healthy.input <- GetAssayData(healthy.object,assay="RNA",slot = "data")  ##取出data
healthy.meta <- healthy.object@meta.data[,c("cell_type","group")]
healthy.meta$cell_type %<>% as.vector()
###     cellchat对象构建，需要两个数据##
healthy.cellchat <- createCellChat(object=healthy.input)
healthy.cellchat <- addMeta(healthy.cellchat,meta=healthy.meta)
healthy.cellchat <- setIdent(healthy.cellchat,ident.use="cell_type")

groupSize <- as.numeric(table(healthy.cellchat@idents))##各种类型的细胞数量

##cellchat提供的人的配受体数据库
healthy.cellchat@DB <- CellChatDB.human
healthy.cellchat <- subsetData(healthy.cellchat,features = NULL)
future::plan("multisession",workers=10)
healthy.cellchat <- identifyOverExpressedGenes(healthy.cellchat)##寻找一下高变基因
healthy.cellchat <- identifyOverExpressedInteractions(healthy.cellchat)##寻找高变基因的相互作用（通路）
healthy.cellchat <- projectData(healthy.cellchat,PPI.human)##将基因投射到PPI（蛋白蛋白网络互作），原因计算基因是不发挥作用的，需要翻译成蛋白执行

##cellchat分析##
##通过计算与每个信号通路相关的所有配体-受体相互作用的通信概率来推断信号通路水平上的通信概率
healthy.cellchat <- computeCommunProb(healthy.cellchat,raw.use = T)##raw.use默认为T，表示使用原始矩阵，也就是我们一开始提取的标准化表达矩阵；如果用的是F，那就是预处理之后的投射到ppi的数据，可能会引入一些非生物学因素
##过滤掉小于10个细胞的胞间通讯网络，通讯中的细胞很少没有意义
healthy.cellchat <- filterCommunication(healthy.cellchat,min.cells = 10)
##通过汇总所有相关的配体/受体，计算信号通路水平上的通信概率
healthy.cellchat <- computeCommunProbPathway(healthy.cellchat)
healthy.cellchat <- aggregateNet(healthy.cellchat)##计算聚合网络
##"netP"表示推断的信号通路的细胞间通信网络
healthy.cellchat <- netAnalysis_computeCentrality(healthy.cellchat,slot.name = "netP")

##结果数据查看/保存
group1.net <- subsetCommunication(healthy.cellchat)
write.csv(group1.net,file = "healthy.csv")
saveRDS(healthy.cellchat,"healthy.rds")
##结果可视化
##互作网络
par(mfrow=c(1,2))
netVisual_circle(healthy.cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge = F,arrow.size = 0.05,
                 title.name = "Number of interactions")
netVisual_circle(healthy.cellchat@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge = F,arrow.size = 0.05,
                 title.name = "Number of interactions")

##展示每个亚群作为source的信号传递
mat <- healthy.cellchat@net$weight
par(mfrow=c(4,4),mar=c(1,1,1,1))
for(i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,arrow.width = 0.2,arrow.size = 0.1,weight.scale = T,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}

#自定义，把某一个细胞单独拎出来看一下
par(mfrow=c(2,2),mar=c(1,1,1,1))
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[1,] <- mat[1,]
  netVisual_circle(mat2,vertex.weight = groupSize,arrow.width = 0.2,arrow.size = 0.1,weight.scale = T,edge.weight.max = max(mat),title.name = rownames(mat)[1])

##pathway选择感兴趣的通路进行可视化
  ##层级图
  healthy.cellchat@netP$pathways
  pathway.show <- group1.net$pathway_name
  pathway.show <- "MIF"##自定义
  levels(healthy.cellchat@idents)
  vertxt.receiver <- seq(8,12)##自定义
  netVisual_aggregate(healthy.cellchat,signaling = pathway.show,vertex.receiver = vertxt.receiver,layout = "hierarchy")
  ##在层次图中，实体圆和空心园分别表示源和目标
  ##线越粗，互作信号越强
  ##左图中间的target是我们选定的靶细胞
  ##右图是选中的靶细胞之外的另外一组放在中间看互作
  
##circle plot
  par(mfrow=c(1,1))
  netVisual_aggregate(healthy.cellchat,signaling = pathway.show,layout = "circle")
  
##热图
  netVisual_heatmap(healthy.cellchat,signaling = pathway.show,color.heatmap= "Reds")
  
##计算受配体对整个信号通路的贡献，并可视化单个配受体对介导的细胞-细胞通讯
  netAnalysis_contribution(healthy.cellchat,signaling = pathway.show)
##可视化由单个配体-受体对介导的细胞间通讯
  pairLR.CCL <- extractEnrichedLR(healthy.cellchat,signaling = pathway.show,geneLR.return = FALSE )
  ##提取对这个通路贡献最大的配体受体对来展示（也可选择其他的配受体对来展示）
  LR.show <- pairLR.CCL[1,]
  netVisual_individual(healthy.cellchat,signaling = pathway.show,pairLR.use = LR.show,layout = "circle")
  netVisual_individual(healthy.cellchat,signaling = pathway.show,pairLR.use = LR.show,layout = "chord")

  
##指定受体细胞和配体细胞
  ##这里指定B细胞，查看B配体和其他细胞受体
  netVisual_bubble(healthy.cellchat,targets.use  = "B",remove.isolate = FALSE,font.size = 14)

##选择淋巴细胞对应髓系细胞，source和target.use的选择
  netVisual_bubble(healthy.cellchat,targets.use  = "B",remove.isolate = FALSE,font.size = 14)###自己改
##参与某条信号通路的所有基因在细胞群中的表达情况展示（小提琴图和气泡图）
  plotGeneExpression(healthy.cellchat,signaling = "MIF")
  ##定义颜色渐变
  colors <- brewer.pal(8,"Oranges")
  plotGeneExpression(healthy.cellchat,signaling = "MIF",type = "dot",col=colors)

##计算和可视化网络中心性评分
  ##主要主要sender和receiver
  healthy.cellchat <- netAnalysis_computeCentrality(healthy.cellchat,slot.name = "netP")
  netAnalysis_signalingRole_network(healthy.cellchat,signaling = pathway.show,width = 15,height = 6,font.size = 10)
  ##使用散点图在2D空间中可视化主要的发送者（源）和接收者（目标）
  netAnalysis_signalingRole_scatter(healthy.cellchat)  ##ALL
  netAnalysis_signalingRole_scatter(healthy.cellchat,signaling = "MIF")  ##
  ##识别对某些细胞类群的输出和输入信号贡献最大的信号
  netAnalysis_signalingRole_heatmap(healthy.cellchat,pattern = "outgoing")
  netAnalysis_signalingRole_heatmap(healthy.cellchat,pattern = "incoming")
  netAnalysis_signalingRole_heatmap(healthy.cellchat,signaling = "MIF",pattern = "outgoing")
  netAnalysis_signalingRole_heatmap(healthy.cellchat,signaling = "MIF",pattern = "incoming")
  
##分泌细胞的信号流出通讯模式的鉴定和可视化（细胞通讯的聚类）
  selectK(healthy.cellchat,pattern = "outgoing")
  nPatterns=2
  healthy.cellchat <- identifyCommunicationPatterns(healthy.cellchat,pattern = "outgoing",k = 2 ,width = 5,height = 9,font.size = 6)
  
  
##组别之间的可视化
  ##合并两组结果
  object.list <- list(healthy=healthy.cellchat,PMN=PMN.cellchat)
  cellchat <- mergeCellChat(object.list,add.names = names(object.list))
  
  ##细胞间互作次数bar图
  ##比较两组互作数目
  gg1 <- compareInteractions(cellchat,show.legend = F,group = c(1,2),color.use  = c('#6699CC','#CC3333'))
  gg2 <- compareInteractions(cellchat,show.legend = F,group = c(1,2),measure = "weight",color.use  = c('#6699CC','#CC3333') )
  gg1+gg2
  ##细胞间互作次数网络图
  ##红色为第二个数据集相比与第一个数据集增加的信号，即为PMN组相比healthy组互作次数和互作强度增加，线越粗表示差异越大
  ##两个数据集之间的细胞-细胞通讯网络中的交互和交互程度的差异数量可以使用圆形图来可视化，其中红色表示第二个数据集相比于第一个数据集增加的信号
  par(mfrow=c(1,2),xpd=TRUE)
  netVisual_diffInteraction(cellchat,weight.scale = T)
  netVisual_diffInteraction(cellchat,weight.scale = T,measure = "weight")
  
  ##数量与强度差异的热图（主要看互作强度，看质不看量）
  par(mfrow=c(1,1))
  h1 <- netVisual_heatmap(cellchat)
  h2 <- netVisual_heatmap(cellchat,measure = "weight")
  h1+h2
  
  ###观察两组的source和target的区别（scatter）
  num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  gg <- list()
  for (i in 1:length(object.list)) {
    gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
  }

  
  ##保守和特异性信号通路的识别与可视化
  gg1 <- rankNet(cellchat,mode="comparison",stacked = T,do.stat = TRUE)
  gg2 <- rankNet(cellchat,mode="comparison",stacked = F,do.stat = TRUE )
  gg1+gg2
  
  ##识别通讯中的配受体贡献度
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  
  ##组别间的cellchat差异分析
  diff.count <- cellchat@net$PMN$count - cellchat@net$healthy$count
  pheatmap(diff.count,treeheight_row = "0",treeheight_col = "0",cluster_rows = T,cluster_cols = T)  ##结果与组别之间画的圈线图一致
  
  
  



  
  
  

          








##pyscenic
  loom <- open_loom('out_SCENIC.loom') 
  regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat)
  regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
  regulonAucThresholds <- get_regulon_thresholds(loom)
  tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])
  embeddings <- get_embeddings(loom)  
  close_loom(loom)
  sub_regulonAUC <- regulonAUC[,match(colnames(B),colnames(regulonAUC))]
  dim(sub_regulonAUC)
  identical(colnames(sub_regulonAUC), colnames(B))  ##确认是否一致
  cellClusters <- data.frame(row.names = colnames(B),seurat_clusters = as.character(B$celltype))
  cellTypes <- data.frame(row.names = colnames(B),celltype = B$celltype)
  head(cellTypes)
  head(cellClusters)
  sub_regulonAUC[1:4,1:4]
  
  ###TF活性均值
  # 看看不同单细胞亚群的转录因子活性平均值
  # Split the cells by cluster:
  selectedResolution <- "celltype"
  cellsPerGroup <- split(rownames(cellTypes),cellTypes[,selectedResolution])
  
  # 去除extened regulons
  sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
  dim(sub_regulonAUC)
  
  # Calculate average expression:
  regulonActivity_byGroup <- sapply(cellsPerGroup,function(cells) rowMeans(getAUC(sub_regulonAUC)[,cells]))
  
  # Scale expression. 
  # Scale函数是对列进行归一化，所以要把regulonActivity_byGroup转置成细胞为行，基因为列
  regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),center = T, scale=T)) 
  # 同一个regulon在不同cluster的scale处理
  dim(regulonActivity_byGroup_Scaled)
  regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
  

##以下为可视化的两种方式
#①
Heatmap(
  regulonActivity_byGroup_Scaled,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = TRUE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

##②
rss=regulonActivity_byGroup_Scaled
head(rss)
df = do.call(rbind,
             lapply(1:ncol(rss), function(i){
               dat= data.frame(
                 path  = rownames(rss),
                 cluster =   colnames(rss)[i],
                 sd.1 = rss[,i],
                 sd.2 = apply(rss[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster) 
n = rss[top5$path,] 
#rownames(rowcn) = rownames(n)
pheatmap(n,
         annotation_row = rowcn,
         show_rownames = T,cluster_cols  = FALSE,cluster_rows = FALSE)




##③

CD8_rss <- as.data.frame(rss)#rss特异性TF结果
#需要作图的细胞类型
celltype <- c("CMC1+ Tem-2","CD8+ Temra-1","CD8+ Temra-2")
rssRanklist <- list()
for(i in 1:length(celltype)) {
  
  data_rank_plot <- cbind(as.data.frame(rownames(regulonActivity_byGroup_Scaled)),
                          as.data.frame(regulonActivity_byGroup_Scaled[,celltype[i]]))#提取数据
  
  colnames(data_rank_plot) <- c("TF", "celltype")
  data_rank_plot=na.omit(data_rank_plot)#去除NA
  data_rank_plot <- data_rank_plot[order(data_rank_plot$celltype,decreasing=T),]#降序排列
  data_rank_plot$rank <- seq(1, nrow(data_rank_plot))#添加排序
  
  p <- ggplot(data_rank_plot, aes(x=rank, y=celltype)) + 
    geom_point(size=3, shape=16, color="#1F77B4",alpha =0.4)+
    geom_point(data = data_rank_plot[1:6,],
               size=3, color='#DC050C')+ #选择前6个标记，自行按照需求选择
    theme_bw()+
    theme(axis.title = element_text(colour = 'black', size = 20),
          axis.text = element_text(colour = 'black', size = 10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),title = element_text(size = 20))+
    labs(x='Regulons Rank', y='Specificity Score',title =celltype[i])+
    geom_text_repel(data= data_rank_plot[1:6,],
                    aes(label=TF), color="black", size=10, fontface="italic", 
                    arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 1,
                    point.padding = 1, segment.color = 'black', 
                    segment.size = 1, force = 1, max.iter = 3e3,max.overlaps = Inf)
  rssRanklist[[i]] <- p
}

plot_grid(rssRanklist[[1]])












CD4@active.ident <-factor(CD4@active.ident,levels = c("0","2","1","3","4","5","7","8","6","9"))
marker.plot <- c("CCR7","LEF1","SELL","ANXA1","VIM","GNLY","GZMK","KLRB1","RORA","RORC","CCR6","RPS26","ISG15","IFI44L","IFI6","FOXP3","NEAT1")
CD4 <- RenameIdents(CD4,"0"="Naive CD4+ T","2"="Naive CD4+ T","1"="CD4+ Tcm-1", "3"="CD4+ Tcm-2","4"="CD4+ Tem","5"="Th17","7"="RPS26+ CD4+ T","8"= "ISG+ CD4+ T","6"="Treg","9"="NEAT1+ CD4+ T")
DotPlot(CD4, features = marker.plot )+ scale_colour_gradientn(colours = colors,na.value="transparent")+theme(axis.title.y=element_text(size=8),axis.text.y =element_text(size = 10),axis.text.x = element_text(size = 10,angle = 90,vjust = 0.5),strip.background =element_rect(colour = "white")) +labs(title = "CD4")

CD8@active.ident <-factor(CD8@active.ident,levels = c("0","8","5","4","3","2","9","1","6","7"))
marker.plot <- c("CCR7","LEF1","SELL","GATA3","SLC4A10","GZMK","CMC1","CCL4","LAYN","CTLA4","GNLY","GZMB","HOPX","KLRD1","NEAT1")
CD8 <- RenameIdents(CD8,"0"="Naive CD8+ T","8"="Naive CD8+ T","5"="GATA3+ Tcm","4"="MAIT cell","3"="CD8+ Tem-1","2"="CMC1+ Tem-2","9"="Exhaust CD8+ T","1"="CD8+ Temra-1", "6"="CD8+ Temra-2","7"= "NEAT1+ CD8+ T")
DotPlot(CD8, features = marker.plot )+ scale_colour_gradientn(colours = colors,na.value="transparent")+theme(axis.title.y=element_text(size=8),axis.text.y =element_text(size = 10),axis.text.x = element_text(size = 10,angle = 90,vjust = 0.5),strip.background =element_rect(colour = "white"))+labs(title = "CD8") 



B@active.ident <-factor(B@active.ident,levels = c("0","3","5","7","1","2","4","6"))
marker.plot <- c("CD27","CD38","IGHM","IGHD","CD19","MME","SOX4","ID3","CXCR5","ITGAX","ISG15","IFI44L","IFIT3","FCRL5","MZB1","IGHG4")
B <- RenameIdents(B,"0"="Naive B cell","3"="Transitional B cell","5"="DN1 B cell", "7"="ISG+ B cell","1"="IgM Memory B cell","2"="Memory B cell","4"="FCRL5+ B cell","6"="Plasma cell")
DotPlot(B, features = marker.plot )+ scale_colour_gradientn(colours = colors,na.value="transparent")+theme(axis.title.y=element_text(size=8),axis.text.y =element_text(size = 10),axis.text.x = element_text(size = 10,angle = 90,vjust = 0.5),strip.background =element_rect(colour = "white")) 


gama@active.ident <- factor(gama@active.ident,levels = c("0","4","1","2","5","3"))
marker.plot <- c("GZMK","DUSP2","IL7R","FCGR3A","CD8A","NKG7","GZMH","GZMB","GNLY","KLRC1","NEAT1","CCR7","LEF1")
gama <- RenameIdents(gama,"0"="Memory γδT cell","4"="Memory γδT cell","1"="KLRC2+ γδT cell","2"="KLRC1+ γδT cell","5"="NEAT1+ γδT cell","3"="Naive γδT cell")
DotPlot(gama, features = marker.plot )+ scale_colour_gradientn(colours = colors,na.value="transparent")+theme(axis.title.y=element_text(size=8),axis.text.y =element_text(size = 10),axis.text.x = element_text(size = 10,angle = 90,vjust = 0.5),strip.background =element_rect(colour = "white")) +labs(title = "γδT")



NK@active.ident <- factor(NK@active.ident,levels = c("3","0","7","5","6","1","2","4"))
marker.plot <- c("NFKBIA","FOS","CD160","IFI44L","IFI6","ISG15","IL32","CX3CR1","FCGR3A","NCAM1","CD3E","NKG7")
NK <- RenameIdents(NK,"3"="NFKBIA+ NK","0"="CD160hi NK","7"="ISG+ NK","5"="IL32+ NK","6"="CD56bright  CD16hi NK","1"="NKT","2"="NKT","4"="NKT")
DotPlot(NK, features = marker.plot )+ scale_colour_gradientn(colours = colors,na.value="transparent")+theme(axis.title.y=element_text(size=8),axis.text.y =element_text(size = 10),axis.text.x = element_text(size = 10,angle = 90,vjust = 0.5),strip.background =element_rect(colour = "white")) +labs(title = "NK")


Mono@active.ident <- factor(Mono@active.ident,levels = c("2","1","0","6","4","7","3","5"))
marker.plot <- c("PADI4","FOS","KLF6","LGALS2","PPBP","PF4","IFI44L","ISG15","IFIT3","TMEM176A","TMEM176B","CD14","FCGR3A")
Mono <- RenameIdents(Mono,"2"="CD14+ PADI4+ Mono","1"="CD14+ FOS+ Mono","0"="CD14+ LGALS2+ Mono","6"="CD14+ PPBP+ Mono","4"="CD14+ ISG+ Mono","7"="CD14+ Mono DC","3"="Intermediate Mono","5"="CD16+ Mono")
DotPlot(Mono, features = marker.plot )+ scale_colour_gradientn(colours = colors,na.value="transparent")+theme(axis.title.y=element_text(size=8),axis.text.y =element_text(size = 10),axis.text.x = element_text(size = 10,angle = 90,vjust = 0.5),strip.background =element_rect(colour = "white")) +labs(title = "Mono")


DC@active.ident <- factor(DC@active.ident,levels = c("0","1","3","2"))
marker.plot <- c("FCER1A","CD1C","FSCN1","SERPINA1","CCL5","IRF8")
DC <- RenameIdents(DC,"0"="cDC-1","1"="cDC-2","3"="cDC-3","2"="pDC")
DotPlot(DC, features = marker.plot )+ scale_colour_gradientn(colours = colors,na.value="transparent")+theme(axis.title.y=element_text(size=8),axis.text.y =element_text(size = 10),axis.text.x = element_text(size = 10,angle = 90,vjust = 0.5),strip.background =element_rect(colour = "white")) +labs(title = "DC")




pbmc.list.CD4 <- list()
pbmc.list.CD8 <- list()
pbmc.list.T <- list()
for (i in seq_along(ns.2)) {
  ns.2[[i]]$barcode <- rownames(ns.2[[i]]@meta.data)
  
  # 只选择需要合并的列
  meta_to_merge <- ns@meta.data[, c("barcode", "celltype")]  # 替换为你需要的列名
  
  ns.2[[i]]@meta.data <- merge(
    ns.2[[i]]@meta.data,
    meta_to_merge,
    by = "barcode",
    all.x = TRUE
  )
  rownames(ns.2[[i]]@meta.data) <- ns.2[[i]]$barcode
  active.ident <- as.character(ns.2[[i]]@meta.data$celltype)
  names(active.ident) <- rownames(ns.2[[i]]@meta.data)
  ns.2[[i]]@active.ident <- as.factor(active.ident)
  
  pbmc.list.CD4[[i]] <- list(
    full = ns.2[[i]],
    CD4 = subset(ns.2[[i]], idents = "CD4+ T")
  )
  pbmc.list.CD8[[i]] <- list(
    full = ns.2[[i]],
    CD8 = subset(ns.2[[i]], idents = "CD8+ T")
  )
  pbmc.list.T[[i]] <- list(
    full = ns.2[[i]],
    CD4 = subset(ns.2[[i]], idents = c("CD4+ T","CD8+ T"))
  )
}









##作图相关
##1.基因评分的小提琴图
 ggplot(metadata,aes(x=seurat_clusters,y = Chemokine.score,fill=seurat_clusters))+
  geom_violin(show.legend = FALSE)+
  theme_classic()+
  stat_compare_means(method = "anova",label.x = 1,label.y = 0.65)+        
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.",size=8)+
  geom_boxplot(width=0.5,show.legend = FALSE,outlier.colour = NA)+
  geom_hline(yintercept = mean(metadata$Chemokine.score),linetype="dashed")+
  theme(axis.text.x = element_text(size = 5,face = "bold"))+
  theme(axis.text.y = element_text(size = 10,face = "bold"))+
  theme(axis.title.y = element_text(size = 15,face = "bold"))


  



##美化umap图
CD4[['cellType']] = CD4@active.ident
umap = CD4@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cellType = CD4@meta.data$cellType)
umap$cellType <- as.character(umap$cellType)
umap$cellType[umap$cellType=="Naive CD4+ T"] <- "1: Naive CD4+ T"
umap$cellType[umap$cellType=="CD4+ Tcm-1"] <- "2: CD4+ Tcm-1"
umap$cellType[umap$cellType=="CD4+ Tcm-2"] <- "3: CD4+ Tcm-2"
umap$cellType[umap$cellType=="CD4+ Tem"] <- "4: CD4+ Tem"
umap$cellType[umap$cellType=="Th17"] <- "5: Th17"
umap$cellType[umap$cellType=="RPS26+ CD4+ T"] <- "6: RPS26+ CD4+ T"
umap$cellType[umap$cellType=="ISG+ CD4+ T"] <- "7: ISG+ CD4+ T"
umap$cellType[umap$cellType=="Treg"] <- "8: Treg"
umap$cellType[umap$cellType=="NEAT1+ CD4+ T"] <- "9: NEAT1+ CD4+ T"

umap$cellType <- factor(umap$cellType,levels = c("1: Naive CD4+ T","2: CD4+ Tcm-1","3: CD4+ Tcm-2","4: CD4+ Tem","5: Th17",
                                                 "6: RPS26+ CD4+ T","7: ISG+ CD4+ T","8: Treg","9: NEAT1+ CD4+ T"))
cell_type_med <- umap %>%
  group_by(cellType) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )
cell_type_med$ident <- c(1:9)
Theme2<-theme(title = element_text(size = 14,face = "bold"),
              panel.background = element_blank(),panel.border = element_blank(),panel.grid=element_blank(), 
              axis.title = element_text(color='black',size=16),axis.ticks.length = unit(0.4,"lines"),
              axis.ticks = element_blank(),axis.line = element_blank(),axis.text=element_blank(),
              legend.title = element_text(size = 14,face = "plain"),legend.text=element_text(size=14),
              legend.key=element_blank(),legend.key.height=unit(0.5,'cm'),
              legend.position = "right",aspect.ratio = 1,plot.title = element_text(hjust = 0.5)) 
ggplot(umap,aes(x=UMAP_1,y=UMAP_2))+geom_point(aes(color=cellType),size=0.5)+scale_color_manual(values = c(colour))+
  geom_label(aes(label=ident),cell_type_med,nudge_x=0,alpha=0,size=6,label.size = NA)+labs(x=" ",y=" ")+theme_bw()+Theme2+
  guides(colour = guide_legend(override.aes = list(size=3.5),ncol = 1))+
  geom_segment(aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2),xend = min(umap$UMAP_1) +3.5, yend = min(umap$UMAP_2)),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap$UMAP_1),y = min(umap$UMAP_2),xend = min(umap$UMAP_1),yend = min(umap$UMAP_2) + 3.5),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+
  annotate("text", x = min(umap$UMAP_1) +2, y = min(umap$UMAP_2) -1, label = "UMAP_1",color="black",size = 4, fontface="bold" ) +
  annotate("text", x = min(umap$UMAP_1) -1, y = min(umap$UMAP_2) + 2, label = "UMAP_2",color="black",size = 4, fontface="bold" ,angle=90)+
  ggtitle("CD4")









###cellphonedb
write.table(as.matrix(pbmc3k@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F) 
meta_data <- cbind(rownames(pbmc3k@meta.data), pbmc3k@meta.data[,'celltype', drop=F]) 
meta_data <- as.matrix(meta_data) 
meta_data[is.na(meta_data)] = "Unkown" # 细胞类型中不能有NA 
write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)


all_pval = read.table("pvalues.txt", header=T, stringsAsFactors = F, sep='\t', comment.char = '', check.names=F)
all_means = read.table('means.txt', header=T, stringsAsFactors = F, sep='\t', comment.char = '', check.names=F)
#挑选需要的列，前面10列和作图没有关系，interacting_pair列是所有的受配体对
intr_pairs = all_pval$interacting_pair
all_pval = all_pval[,-c(1,3:11)]
all_means = all_means[,-c(1,3:11)]
selected_celltype = c("CD4+ T|CD8+ Temra-2","B|CD8+ Temra-2","Plasma cell|CD8+ Temra-2", "Mono|CD8+ Temra-2","cDC|CD8+ Temra-2","pDC|CD8+ Temra-2","Neu|CD8+ Temra-2")
sig_pairs <- all_pval
sig_pairs <- sig_pairs[which(rowSums(sig_pairs<=0.05)!=0), ]
View(sig_pairs)
sub_sig_pairs <- sig_pairs
dim(sub_sig_pairs)
selected_pairs = sub_sig_pairs$interacting_pair
sel_pval = all_pval[match(selected_pairs, intr_pairs), selected_celltype]#将上述需要呈现的受配体选出来
sel_means = all_means[match(selected_pairs, intr_pairs), selected_celltype]
View(sel_pval)
df_names = expand.grid(selected_pairs, selected_celltype)
pval = unlist(sel_pval)
pval[pval==0] = 0.00001
plot.data = cbind(df_names,pval)
pr = unlist(as.data.frame(sel_means))
plot.data = cbind(plot.data,pr)
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
plot.data$clusters <- gsub('[|]', '_', plot.data$clusters)
my_palette <- colorRampPalette(c("darkblue","yellow","red"))(n=1000)
# 指定要筛选的pair值
selected_pairs <- c("BTLA_TNFRSF14", "ICAM1_SPN", "ICAM1_ITGAL", "CD99_PILRA","CD40_CD40LG")

# 使用%in%操作符筛选出包含这些pair值的行
plot.data1 <- plot.data[plot.data$pair %in% selected_pairs, ]
ggplot(plot.data1, aes(x=clusters, y=pair)) +
  geom_point(aes(size=-log10(pvalue + 0.001), color=mean)) +
  scale_size_continuous(range=c(3, 7), breaks=c(0, 1.0, 2.0)) +  # 调整点的大小范围
  scale_color_gradientn('Mean expression', colors=my_palette, limits=c(0.3, 1.5), na.value="transparent") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=10, colour="black"),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0, size=10),
        axis.text.y = element_text(size=10, colour="black"),
        axis.title=element_blank(),
        panel.border = element_rect(size=0.7, linetype="solid", colour="black"))


all_means <-  read.delim("CD8.merge.cellphonedb/means.txt", check.names = FALSE)
all_pval <-  read.delim("CD8.merge.cellphonedb/pvalues.txt", check.names = FALSE)
plot_cpdb4(
  scdata = CD8.merge,interaction = "BTLA-TNFRSF14",
  cell_type1 = ".",
  cell_type2 = "CD8+ Temra-2",
  celltype_key = "celltype", # column name where the cell ids are located in the metadata
  means = all_means,
  pvals = all_pval,
  deconvoluted = decon_stat,desiredInteractions = list(c("pDC","CD8+ Temra-2"),c("pDC","CD4+ T"))
)






