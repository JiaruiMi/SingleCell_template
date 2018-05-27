#=============================================================================================
#
#            single-cell RNA-seq analysis using different packages (methods)
#
#=============================================================================================

#=============================================================================================
#                                      "Seurat" package
#=============================================================================================


### Load packages，加载数据前需要将文件夹中的三个文件分别命名为“matrix.mtx", "barcodes.tsv", "genes.tsv"
library(Seurat) # the seurat version now I am using is 2.3.1
library(dplyr)
library(Matrix)



### set working directory
setwd('/Users/mijiarui/Nature_Biotechnology_Paper/Pancreas_1_mRNA_GSM2830058_P5')
list.files("/Users/mijiarui/Nature_Biotechnology_Paper/Pancreas_1_mRNA_GSM2830058_P5")

### # Load the Pancreas_1_mRNA_GSM2830058_P5 dataset
pancreas_1.data <- Read10X(data.dir = "/Users/mijiarui/Nature_Biotechnology_Paper/Pancreas_1_mRNA_GSM2830058_P5")

### Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = pancreas_1.data))
dense.size

sparse.size <- object.size(x = pancreas_1.data)
sparse.size

dense.size / sparse.size

# Initialize the Seurat object with the raw (non-normalized data).  Keep all genes expressed in >= 3 cells (~0.1% of the data). 
# Keep all cells with at least 200 detected genes

pancreas_1 <- CreateSeuratObject(raw.data = pancreas_1.data, min.cells = 3, min.genes = 200, 
                                 project = "10X_Pancreas_1")
pancreas_1

###################################### Quality control ###############################################
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pancreas_1@data), value = TRUE)
percent.mito <- Matrix::colSums(pancreas_1@raw.data[mito.genes, ]) / Matrix::colSums(pancreas_1@raw.data)


# AddMetaData: adds columns to object@meta.data, and is a great place to stash QC stats
pancreas_1 <- AddMetaData(object = pancreas_1, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pancreas_1, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# GenePlot is typically used to visualize gene-gene relationships, but can be used for anything calculated by the object,
# i.e. columns in object@meta.data, PC scores etc.  Since there is a rare subset of cells with an outlier level of high 
# mitochondrial percentage and also low UMI content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = pancreas_1, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pancreas_1, gene1 = "nUMI", gene2 = "nGene")


# We filter out cells that have unique gene counts over 2,500 or less than 200 Note that low.thresholds and high.thresholds are 
# used to define a 'gate'.  -Inf and Inf should be used if you don't want a lower or upper
# threshold. "可以看到这里选择的QC标准是 200~2500基因范围内，以及线粒体基因表达占比小于5%的才保留。"
pancreas_1 <- FilterCells(object = pancreas_1, subset.names = c("nGene", "percent.mito"), 
                          low.thresholds = c(200, -Inf), high.thresholds = c(Inf, 0.05))


###################################### Normalization ###############################################
# "这里默认根据细胞测序文库大小进行normalization，简单的做一个log转换即可。"
summary(pancreas_1@raw.data[,1])
pancreas_1 <- NormalizeData(object = pancreas_1, normalization.method = "LogNormalize", 
                            scale.factor = 10000)

summary(pancreas_1@data[,1])


################### Detection of variable genes across the single cells ###################
par(mfrow = c(1, 1))
pancreas_1 <- FindVariableGenes(object = pancreas_1, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pancreas_1@var.genes)


################### Scaling the data and removing unwanted sources of variation ###################
# "需要去除那些technical noise,batch effects, or even biological sources of variation (cell cycle stage)"
pancreas_1 <- ScaleData(object = pancreas_1, vars.to.regress = c("nUMI", "percent.mito"))
summary(pancreas_1@scale.data[,1])

################### PCA分析(Principal component analysis) ###################
pancreas_1 <- RunPCA(object = pancreas_1, pc.genes = pancreas_1@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)


## 对PCA分析结果可以进行一系列的可视化： PrintPCA, VizPCA, PCAPlot, and PCHeatmap
par(mar = c(5,5,3,2))
VizPCA(object = pancreas_1, pcs.use = 1:2)
PCAPlot(object = pancreas_1, dim.1 = 1, dim.2 = 2)

# ProjectPCA scores each gene in the dataset (including genes not included in the PCA) based on their correlation 
# with the calculated components. Though we don't use this further here, it can be used to identify markers that 
# are strongly correlated with cellular heterogeneity, but may not have passed through variable gene selection. 
# The results of the projected PCA can be explored by setting use.full=T in the functions above
pancreas_1 <- ProjectPCA(object = pancreas_1, do.print = FALSE)

## 最重要的就是 PCHeatmap 函数了
PCHeatmap(object = pancreas_1, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = pancreas_1, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

######################### 找到有统计学显著性的主成分 #########################
# 主成分分析结束后需要确定哪些主成分所代表的基因可以进入下游分析，这里可以使用JackStraw做重抽样分析。
# 可以用JackStrawPlot可视化看看哪些主成分可以进行下游分析。这一步很耗时
pancreas_1 <- JackStraw(object = pancreas_1, num.replicate = 100) 
JackStrawPlot(object = pancreas_1, PCs = 1:12)

# 当然，也可以用最经典的碎石图来确定主成分。这一步很耗时。
PCElbowPlot(object = pancreas_1)


# 这个确定主成分是非常有挑战性的: - The first is more supervised, exploring PCs to determine relevant sources of heterogeneity, 
# and could be used in conjunction with GSEA for example. - The second implements a statistical test based on a random null model, 
# but is time-consuming for large datasets, and may not return a clear PC cutoff. The third is a heuristic that 
# is commonly used, and can be calculated instantly.

###################################### Cluster the cells ######################################
# save.SNN = T saves the SNN so that the clustering algorithm can be rerun using the same graph
# but with a different resolution value (see docs for full details)
pancreas_1 <- FindClusters(object = pancreas_1, reduction.type = "pca", dims.use = 1:15, resolution = 0.6, print.output = 0, save.SNN = TRUE)



# A useful feature in Seurat v2.0 is the ability to recall the parameters that were used in the latest function calls 
# for commonly used functions. For FindClusters, we provide the function PrintFindClustersParams to print a nicely 
# formatted formatted summary of the parameters that were chosen.
PrintFindClustersParams(object = pancreas_1)
# While we do provide function-specific printing functions, the more general function to 
# print calculation parameters is PrintCalcParams(). 


################### Run Non-linear dimensional reduction (tSNE) ###################
# 同样也是一个函数，这个结果也可以像PCA分析一下挑选合适的PC进行下游分析。
# 这一步很耗时，可以保存该对象，便于重复，以及分享交流 （出图时间很长)
pancreas_1 <- RunTSNE(object = pancreas_1, dims.use = 1:15, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters 
TSNEPlot(object = pancreas_1, do.label=T)
?TSNEPlot
save(pancreas_1, file = "pancreas_1.rData")
pancreas_1.markers %>% group_by(cluster)


################### Finding differentially expressed genes (cluster biomarkers) ###################
# 差异分析在seurat包里面被封装成了函数：FindMarkers，有一系列参数可以选择，然后又4种找差异基因的算法：

# 方法一： ROC test (“roc”)
# 方法二： t-test (“t”)
# 方法三： LRT test based on zero-inflated data (“bimod”, default)
# 方法四： LRT test based on tobit-censoring models (“tobit”)


# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = pancreas_1, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))

# find markers of each cluster
for (i in 0:11) {
  cluster.markers <- FindMarkers(object = pancreas_1, ident.1 = i, min.pct = 0.25)
  print(paste("marker of cluster: ",i ))
  print(x = head(x = cluster.markers, n = 50))
}

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = pancreas_1, ident.1 = 5, ident.2 = c(0,3), min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))


# find markers for every cluster compared to all remaining cells, report only the positive ones (this step takes some time!)
pancreas_1.markers <- FindAllMarkers(object = pancreas_1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
pancreas_1.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
pancreas_1.markers $cluster

# 值得注意的是： The ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect).
cluster1.markers <- FindMarkers(object = pancreas_1, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)

# 同时，该包提供了一系列可视化方法来检查差异分析的结果的可靠性：

# VlnPlot (shows expression probability distributions across clusters)
# FeaturePlot (visualizes gene expression on a tSNE or PCA plot) are our most commonly used visualizations
# JoyPlot, CellPlot, and DotPlot

VlnPlot(object = pancreas_1, features.plot = c("ins", "gcga","gcgb", 'sst2'))

# you can plot raw UMI counts as well
VlnPlot(object = pancreas_1, features.plot = c("ins", "gcga"), use.raw = TRUE, y.log = TRUE)


FeaturePlot(object = pancreas_1, 
            features.plot = c("ins", "gcga","gcgb", "sst2" ,'TMEM269','zgc:158463',"sox9b", 'cftr',
                              'her15.1','notch1a','try','ela2l','yap1','taz',
                              'ca2', 'nkx6.1', 'hnf1ba','gip', 'lats2','lats1', 'agrn','ptpn21','amotl2a',
                              'fev','gck','prdx4','spcs1','mgst3b','mknk2b','foxo1a','sik2b','ythdf2','igf2bp2a', 'fto',
                              'mettl3','mettl14','wtap','ythdf1', 'ythdf3','ythdc1', 'alkbh5',
                              'jag1a','jag1b','crb3a','crb3b','anxa4'), 
            cols.use = c("grey", "Red"), reduction.use = "tsne")

# DoHeatmap generates an expression heatmap for given cells and genes. In this case, we are plotting the top 20 markers 
# (or all markers if less than 20) for each cluster.
pancreas_1.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
DoHeatmap(object = pancreas_1, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)







#=============================================================================================
#                            "Monocle2" package (Monocle 2.8.0)
#=============================================================================================
# 在nature methods杂志发表的文章，更新为monocle2版本并且更换了主页，功能也不仅仅是差异分析那么简单。还包括pseudotime,
# clustering分析，而且还可以进行基于转录本的差异分析，其算法是BEAM (used in branch analysis) and Census (the core of relative2abs)，
# 也单独发表了文章。
# 用了4个公共的数据来测试说明其软件的用法和优点。
# the HSMM data set, GSE52529 (ref. 1);
# the lung data set, GSE52583 (ref. 8);
# the Paul et al. data set ;
# the Olsson data set9, synapse ID syn4975060.
# 也是有着非常详细的使用教程 , 读取表达矩阵和分组信息，需要理解其定义好的一些S4对象。

###################################### S4 对象 ######################################
# 主要是基于 CellDataSet 对象来进行下游分析，继承自ExpressionSet对象，也是常见的3个组成：
# exprs, a numeric matrix of expression values, where rows are genes, and columns are cells
# phenoData, an AnnotatedDataFrame object, where rows are cells, and columns are cell attributes (such as cell type, culture condition, day captured, etc.)
# featureData, an AnnotatedDataFrame object, where rows are features (e.g. genes), and columns are gene attributes, such as biotype, gc content, etc.
# 创建对象的时候需要指定引入的表达矩阵的方法，monocle2推荐用基于转录本的counts矩阵，同时也是默认的参数 
# expressionFamily=negbinomial.size() ，如果是其它RPKM/TMP等等，需要找到对应的参数

### Load packages: 我们使用的测试数据集是来自于HSMMSingleCell
library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)
library(HSMMSingleCell)
library(monocle)
library(M3Drop)
data(HSMM_expr_matrix) ## RPKM 矩阵,271个细胞，47192个基因
data(HSMM_gene_annotation)
data(HSMM_sample_sheet)
HSMM_expr_matrix[1:10,1:5]
dim(HSMM_expr_matrix)
head(HSMM_gene_annotation)
head(HSMM_sample_sheet)


###################################### 构建S4对象，CellDataSet ######################################
# 主要是读取表达矩阵和样本描述信息，这里介绍两种方式，一种是读取基于 subjunc+featureCounts 分析后的reads counts矩阵，
# 一种是读取 tophat+cufflinks 得到的RPKM表达矩阵
# 如何读取外部数据，请参考“https://mp.weixin.qq.com/s/zCfDkxbVTxjFQ5QAIULYjA”

# 在这里我们读取HSMMSingleCell包中的测试数据
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)


# First create a CellDataSet from the relative expression levels

## 这里仅仅是针对rpkm表达矩阵的读取
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),   
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit=0.1,
                       expressionFamily=tobit(Lower=0.1))

# Next, use it to estimate RNA counts
rpc_matrix <- relative2abs(HSMM)
rpc_matrix[1:10,1:5] 

## rpkm格式的表达值需要转换成reads counts之后才可以进行下游分析！

# Now, make a new CellDataSet using the RNA counts
HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit=0.5,
                       expressionFamily=negbinomial.size())

## 下面的分析，都基于内置数据构建的S4对象，HSMM

###################################### 过滤低质量细胞和未检测到的基因 ######################################
# 基于基因的过滤
##  这里只是把 基因挑选出来，并没有对S4对象进行过滤操作。 这个 detectGenes 函数还计算了 每个细胞里面表达的基因数量。

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
### Warning: Deprecated, use tibble::rownames_to_column() instead.

HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))
## 对每个基因都检查一下在多少个细胞里面是有表达量的。
## 只留下至少在10个细胞里面有表达量的那些基因，做后续分析
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
length(expressed_genes) ## 只剩下了14224个基因
print(head(pData(HSMM))) 

# 基于样本表达量进行过滤
## 这里选择的是通过不同时间点取样的细胞来进行分组查看，把超过2个sd 的那些样本的临界值挑选出来，下一步过滤的时候使用。
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) +
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) -
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))
table(pData(HSMM)$Hours)
qplot(Total_mRNAs, data = pData(HSMM), color = Hours, geom = "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

# 执行过滤并可视化检查一下
## 上面已经根据基因表达情况以及样本的总测序数据选择好了阈值，下面就可以可视化并且对比检验一下执行过滤与否的区别。
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound & 
               pData(HSMM)$Total_mRNAs < upper_bound]                                 
HSMM <- detectGenes(HSMM, min_expr = 0.1)
L <- log(exprs(HSMM[expressed_genes,]))
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
qplot(value, geom="density", data=melted_dens_df) +  stat_function(fun = dnorm, size=0.5, color='red') + 
  xlab("Standardized log(FPKM)") +
  ylab("Density")


########################################### 聚类 ###########################################
# 根据指定基因对单细胞转录组表达矩阵进行分类
## 下面这个代码只适用于这个测试数据， 主要是生物学背景知识，用MYF5基因和ANPEP基因来对细胞进行分类，可以区分Myoblast和Fibroblast。
## 如果是自己的数据，建议多读读paper看看如何选取合适的基因，或者干脆跳过这个代码。
## 根据基因名字找到其在表达矩阵的ID，这里是ENSEMBL数据库的ID
MYF5_id <- row.names(subset(fData(HSMM), gene_short_name == "MYF5"))
ANPEP_id <- row.names(subset(fData(HSMM), gene_short_name == "ANPEP"))
## 这里选取的基因取决于自己的单细胞实验设计
cth <- newCellTypeHierarchy()

cth <- addCellType(cth, "Myoblast", classify_func = function(x) { x[MYF5_id,] >= 1 })
cth <- addCellType(cth, "Fibroblast", classify_func = function(x){ x[MYF5_id,] < 1 & x[ANPEP_id,] > 1 })

HSMM <- classifyCells(HSMM, cth, 0.1)
HSMM <- clusterCells(HSMM)
pData(HSMM)$CellType
### Warning: Deprecated, use tibble::rownames_to_column() instead.

### Warning: Deprecated, use tibble::rownames_to_column() instead.
### 这个时候的HSMM已经被改变了，增加了属性。
table(pData(HSMM)$CellType)

pie <- ggplot(pData(HSMM), aes(x = factor(1), fill = factor(CellType))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

### 可以看到还有很大一部分细胞仅仅是根据这两个基因的表达量是无法成功的归类的。这个是很正常的，
### 因为单细胞转录组测序里面的mRNA捕获率不够好。 通过这个步骤成功的给HSMM这个S4对象增加了一个属性，就是CellType，
### 在下面的分析中会用得着。

## 无监督聚类: 这里需要安装最新版R包才可以使用里面的一些函数，因为上面的步骤基于指定基因的表达量进行细胞分组会漏掉很多信息，
## 所以需要更好的聚类方式。

disp_table <- dispersionTable(HSMM)
head(disp_table)

## 只有满足 条件的10198个基因才能进入聚类分析
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)


## 这里看看基因的表达量和基因的变异度之间的关系
## 处在灰色阴影区域的基因会被抛弃掉，不进入聚类分析。
plot_pc_variance_explained(HSMM, return_all = F) # norm_method = 'log',

HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 6, 
                        reduction_method = 'tSNE', verbose = T) 
HSMM <- clusterCells(HSMM, num_clusters=2)
## Distance cutoff calculated to 1.072748
## 这里先用tSNE的聚类方法处理HSMM数据集，并可视化展示
plot_cell_clusters(HSMM, 1, 2, color="CellType", markers=c("MYF5", "ANPEP"))
## 可以看到并不能把细胞类型完全区分开，这个是完全有可能的，因为虽然是同一种细胞，但是有着不同的培养条件。
head(pData(HSMM))
head(fData(HSMM))
## 所以这里也区分一下 培养基， a high-mitogen growth medium (GM) to a low-mitogen differentiation medium (DM). 
plot_cell_clusters(HSMM, 1, 2, color="Media")



## 因为我们假设就2种细胞类型，所以在做聚类的时候可以把这个参数添加进去，这样可以去除无关变量的干扰。
HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 2, reduction_method = 'tSNE', 
                        residualModelFormulaStr="~Media + num_genes_expressed", verbose = T) #
HSMM <- clusterCells(HSMM, num_clusters=2)
## Distance cutoff calculated to 1.284778
plot_cell_clusters(HSMM, 1, 2, color="CellType") 


plot_cell_clusters(HSMM, 1, 2, color="Cluster") + facet_wrap(~CellType)



# 半监督聚类
## 这里的差异分析非常耗时
marker_diff <- markerDiffTable(HSMM[expressed_genes,], 
                               cth, 
                               residualModelFormulaStr="~Media + num_genes_expressed",
                               cores=1)
head(marker_diff)

## 就是对每个基因增加了pval和qval两列信息，挑选出那些在不同media培养条件下显著差异表达的基因，310个，
candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.01))

## 计算这310个基因在不同的celltype的specificity值
marker_spec <- calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)
head(selectTopMarkers(marker_spec, 3)) 

semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
plot_ordering_genes(HSMM)

## 重新挑选基因，只用黑色高亮的基因来进行聚类。
plot_pc_variance_explained(HSMM, return_all = F) # norm_method = 'log',

HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 2, reduction_method = 'tSNE', 
                        residualModelFormulaStr="~Media + num_genes_expressed", verbose = T) 
HSMM <- clusterCells(HSMM, num_clusters=2) 
## Distance cutoff calculated to 1.02776
plot_cell_clusters(HSMM, 1, 2, color="CellType")

HSMM <- clusterCells(HSMM,
                     num_clusters=2, 
                     frequency_thresh=0.1,
                     cell_type_hierarchy=cth)
## Distance cutoff calculated to 1.02776
plot_cell_clusters(HSMM, 1, 2, color="CellType", markers = c("MYF5", "ANPEP"))

pie <- ggplot(pData(HSMM), aes(x = factor(1), fill = factor(CellType))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

########################################### Pseudotime分析 ###########################################
# 主要目的是：Constructing Single Cell Trajectories

# 发育过程中细胞状态是不断变化的，monocle包利用算法学习所有基因的表达模式来把每个细胞安排到各各自的发展轨迹。 
# 在大多数生物学过程中，参与的细胞通常不是同步发展的，只有单细胞转录组技术才能把处于该过程中各个中间状态的细胞分离开来，
# 而monocle包里面的pseudotime分析方法正是要探究这些。
# choose genes that define a cell’s progress
# reduce data dimensionality
# order cells along the trajectory

# 其中第一个步骤挑选合适的基因有3种策略，分别是：
# Ordering based on genes that differ between clusters
# Selecting genes with high dispersion across cells
# Ordering cells using known marker genes

# 无监督的Pseudotime分析
HSMM_myo <- HSMM[,pData(HSMM)$CellType == "Myoblast"]   
HSMM_myo <- estimateDispersions(HSMM_myo)
## Warning: Deprecated, use tibble::rownames_to_column() instead.
## Removing 143 outliers

# 使用不同的策略会给出不同的fData(State)
## 策略1：  Ordering based on genes that differ between clusters
if(F){
  diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,],
                                        fullModelFormulaStr="~Media")
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
}
## 策略2：Selecting genes with high dispersion across cells
disp_table <- dispersionTable(HSMM_myo)
ordering_genes <- subset(disp_table, 
                         mean_expression >= 0.5 & 
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id

HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
plot_ordering_genes(HSMM_myo)
## Warning: Transformation introduced infinite values in continuous y-axis

## 挑选变异度大的基因，如图所示
HSMM_myo <- reduceDimension(HSMM_myo, max_components=2)
HSMM_myo <- orderCells(HSMM_myo)
## 排序好的细胞可以直接按照发育顺序可视化
plot_cell_trajectory(HSMM_myo, color_by="State")
table(pData(HSMM_myo)$State)

########################################### 直接做差异分析 ###########################################
# 前面的聚类分析和Pseudotime分析都需要取基因子集，就已经利用过差异分析方法来挑选那些有着显著表达差异的基因。
# 如果对所有的基因来检验，非常耗时。
marker_genes <- row.names(subset(fData(HSMM_myo), 
                                 gene_short_name %in% c("MEF2C", "MEF2D", "MYF5", 
                                                        "ANPEP", "PDGFRA","MYOG", 
                                                        "TPM1",  "TPM2",  "MYH2", 
                                                        "MYH3",  "NCAM1", "TNNT1", 
                                                        "TNNT2", "TNNC1", "CDK1", 
                                                        "CDK2",  "CCNB1", "CCNB2", 
                                                        "CCND1", "CCNA1", "ID1")))

diff_test_res <- differentialGeneTest(HSMM_myo[marker_genes,], 
                                      fullModelFormulaStr="~Media")
# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)
sig_genes[,c("gene_short_name", "pval", "qval")]

###  还可以挑选其中几个基因来可视化看看它们是如何在不同组差异表达的。这个画图函数自己都可以写。

MYOG_ID1 <- HSMM_myo[row.names(subset(fData(HSMM_myo), 
                                      gene_short_name %in% c("MYOG", "CCNB2"))),]
plot_genes_jitter(MYOG_ID1, grouping="Media", ncol=2)

### 这样就可以测试某些基因，是否能区分细胞群体的不同类型及状态

to_be_tested <- row.names(subset(fData(HSMM), 
                                 gene_short_name %in% c("UBC", "NCAM1", "ANPEP"))) 
cds_subset <- HSMM[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr="~CellType")
diff_test_res[,c("gene_short_name", "pval", "qval")] 
plot_genes_jitter(cds_subset, grouping="CellType", color_by="CellType", 
                  nrow=1, ncol=NULL, plot_trend=TRUE)


full_model_fits <- fitModel(cds_subset, modelFormulaStr="~CellType")
reduced_model_fits <- fitModel(cds_subset, modelFormulaStr="~1")
diff_test_res <- compareModels(full_model_fits, reduced_model_fits)
diff_test_res

plot_genes_in_pseudotime(cds_subset, color_by="Hours")


########################################### 算法 ###########################################
# Monocole还提出了好几个算法：
## dpFeature: Selecting features from dense cell clusters
## Reversed graph embedding
## DRTree: Dimensionality Reduction via Learning a Tree
## DDRTree: discriminative dimensionality reduction via learning a tree
## Census: a normalization method to convert of single-cell mRNA transcript to relative transcript counts.
## BEAM : to test for branch-dependent gene expression by formulating the problem as a contrast between two negative binomial GLMs.
## Branch time point detection algorithm


#=============================================================================================
#                          Using "Monocle2" to analyze data from "Seurat"
#=============================================================================================
library(monocle)
pancreas.1 <- importCDS(pancreas_1)
exprs <- exprs(pancreas.1);dim(exprs)
pancreas_1@var.genes
pData(pancreas.1)[,5]
## 下面的分析，都基于内置数据构建的S4对象

d <- exprs[which(rownames(exprs) %in% pancreas_1@var.genes), ];dim(d)
d <- d[!duplicated(rownames(d)), ]
cellLabels <- as.numeric(pData(pancreas.1)[,5]); table(cellLabels)
geneNames <- rownames(exprs)

## run monocle
colnames(d) <- 1:ncol(d)
geneNames <- rownames(d)
rownames(d) <- 1:nrow(d)
dim(d)
pd <- data.frame(timepoint = cellLabels)
pd <- new("AnnotatedDataFrame", data=pd)
fd <- data.frame(gene_short_name = geneNames)
fd <- new("AnnotatedDataFrame", data=fd)

dCellData <- newCellDataSet(d, phenoData = pd, featureData = fd, expressionFamily = tobit())
dCellData <- setOrderingFilter(dCellData, which(geneNames %in% pancreas_1@var.genes))
dCellData <- estimateSizeFactors(dCellData)
dCellDataSet <- reduceDimension(dCellData, max_components = 12)
dCellDataSet <- orderCells(dCellDataSet, reverse = FALSE)
pData(dCellDataSet)$cluster <- factor(pData(dCellDataSet)$timepoint)
plot_cell_trajectory(dCellDataSet, color_by = 'cluster')

# Store the ordering
pseudotime_monocle <-
  data.frame(
    Timepoint = phenoData(dCellDataSet)$timepoint,
    pseudotime = phenoData(dCellDataSet)$Pseudotime,
    State = phenoData(dCellDataSet)$State
  )
rownames(pseudotime_monocle) <- 1:ncol(d)
pseudotime_order_monocle <-
  rownames(pseudotime_monocle[order(pseudotime_monocle$pseudotime), ])

# We can again compare the inferred pseudotime to the known sampling timepoints.
library(ggbeeswarm);
library(ggthemes)
pancreas.1$pseudotime_monocle <- pseudotime_monocle$pseudotime
pData(pancreas.1)$res.0.6 <- factor(as.numeric(pData(pancreas.1)$res.0.6))
ggplot(as.data.frame(pData(pancreas.1)), 
       aes(x = -pseudotime_monocle, 
           y = res.0.6, colour = res.0.6)) +
  geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("monocle pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by monocle pseudotime")




#=============================================================================================
#                          比较不同的对单细胞转录组数据normalization方法
#=============================================================================================
# 使用CPM去除文库大小影响
# 之所以需要normalization，就是因为测序的各个细胞样品的总量不一样，所以测序数据量不一样，就是文库大小不同，
# 这个因素是肯定需要去除。最简单的就是counts per million (CPM)，所有样本的所有基因的表达量都乘以各自的文库
# reads总数再除以一百万即可。(一般miRNA-seq数据结果喜欢用这个) 代码如下：

calc_cpm <-
  function (expr_mat, spikes = NULL) 
  {
    norm_factor <- colSums(expr_mat[-spikes, ])
    return(t(t(expr_mat)/norm_factor)) * 10^6
  }
# 但是CPM方法有一个很严重的缺陷，那些高表达并且在细胞群体表达差异很大的基因会严重影响那些低表达基因。

# RPKM, FPKM and TPM去除基因或者转录本长度影响
# 最常见的是下面3个：

# RPKM - Reads Per Kilobase Million (for single-end sequencing)
# FPKM - Fragments Per Kilobase Million (same as RPKM but for paired-end sequencing, makes sure that paired 
# ends 
# mapped to the same fragment are not counted twice)
# TPM - Transcripts Per Kilobase Million (same as RPKM, but the order of normalizations is reversed - length 
# first and sequencing depth second)
# 这些normalization方法并不适合单细胞转录组测序数据，因为有一些scRNA-seq建库方法具有3端偏好性，一般是没办法测
# 全长转录本的，所以转录本的长度跟表达量不是完全的成比例。
# 对于这样的数据，需要重新转换成 reads counts 才能做下游分析。

# 适用于bulk RNA-seq的normalization方法
# 比较流行的有：

# DESeq的size factor (SF)
# relative log expression(RLE)
# upperquartile (UQ)
# weighted trimmed mean of M-values(TMM)
# 这些适用于 bulk RNA-seq data 的normalization方法可能并不适合 single-cell RNA-seq data ，因为它们的基本假设是
# 有问题的。

# 特意为single-cell RNA-seq data 开发的normalization方法

# LSF (Lun Sum Factors)
# scran package implements a variant on CPM specialized for single-cell data
# 而scater包把这些normalization方法都包装到了normaliseExprs函数里面，可以直接调用。

# 并且通过plotPCA函数来可视化这些normalization的好坏。

## 工作环境
library(scRNA.seq.funcs)
library(scater)
library(scran)
library(SingleCellExperiment)
set.seed(1234567)

## 加载测试数据
## 这里选取的是芝加哥大学Yoav Gilad lab实验的Tung et al 2017的单细胞测序文章的数据
options(stringsAsFactors = FALSE)
set.seed(1234567)
## 这里直接读取过滤好的数据，是一个SCESet对象，适用于scater包的
## http://www.biotrainee.com/jmzeng/scRNA/tung_umi.rds
setwd('/Users/mijiarui/Nature_Biotechnology_Paper/Testing_dataset')
umi <- readRDS("tung_umi.rds")

## 如果没有这个rds对象，就自己把read counts的表达矩阵读进去，变成这个适用于scater包的SCESet对象，代码如下；
## Load the data and annotations:
molecules <- read.table("molecules.txt", sep = "\t")
anno <- read.table("annotation.txt", sep = "\t", header = TRUE)

umi.qc <- umi[fData(umi)$use, pData(umi)$use] 
## counts(umi) 和  exprs(umi) 这里是不一样的。
## 前面的过滤信息，这里直接用就好了。
endog_genes <- !fData(umi.qc)$is_feature_control
dim(exprs( umi.qc[endog_genes, ]))
## 可以看到是过滤后的654个单细胞的13997个基因的表达矩阵。
umi.qc
toSingleCellExperiment(umi.qc)


########################################### Raw ###########################################
# 先看看原始的表达值的分布情况，这里本来应该是对每一个样本画boxplot的，但是这里的样本数量太多了，
# 这样的可视化效果很差， 就用PCA的方式，看看这表达矩阵是否可以把样本区分开，只有那些区分度非常好的normalization
# 方法才是最优的。不过scater包提供了一个plotRLE函数，可以画出类似于样本boxplot的效果。

plotPCA(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual",
  exprs_values = "log2_counts"
)

########################################### CPM ###########################################
# scater默认对表达矩阵做了cpm转换，所以可以直接提取里面的信息
plotPCA(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual",
  exprs_values = "exprs"
)




























#=============================================================================================
#                          比较不同单细胞转录组数据寻找features方法
#=============================================================================================
# 挑选到的跟feature相关的基因集，有点类似于在某些组间差异表达的基因集，都需要后续功能注释。
# 背景介绍
# 单细胞转录组测序的确可以一次性对所有细胞都检测到上千个基因的表达，但是，大多数情况下，只有其中的少部分基因是有生物学意义的，
# 比如可以区分不同的细胞类型，或者分化发育相关的基因，或者细胞应对外界刺激的。而且大多数基因之所以在不同的细胞里面表达有差异，
# 其实是技术限制，背景噪音。这些技术限制，包括批次效应，都会阻碍我们发现那些真正的有生物学意义的基因。所以做 feature selection 
# 分析来去除那些技术噪音相关基因，可以显著的提高信噪比，降低后续分析的复杂度。
# Feature Selection (feature selection这部最建议使用的方法)
# M3Drop

# 包的安装
library('ROCR')
library(scRNA.seq.funcs)
library(matrixStats)
library(M3Drop)
library(RColorBrewer)
set.seed(1)

# 加载测试数据:这里选取的是Usoskin et al 文章单细胞测序文章的数据，包含4种细胞类型：

# NP = non-peptidergic nociceptors
# PEP = peptidergic nociceptors
# NF = neurofilament containing
# TH = tyrosine hydroxylase containing neurons.
# 对应着25334个基因在 622 个细胞里面的表达量
setwd('/Users/mijiarui/Nature_Biotechnology_Paper/Testing_dataset')
usoskin1 <- readRDS("usoskin1.rds")
dim(usoskin1)
table(colnames(usoskin1))

############################ 用M3Drop对表达矩阵进行一站式过滤 #############################
# 这个M3Drop的M3DropCleanData函数会自动过滤那些表达基因数量很少的细胞，过滤低表达基因，然后把reads counts的表达矩阵转换为 
# counts per million (CPM) 。过滤之后，只剩下 15708个基因在532个细胞的表达了。
uso_list <- M3Drop::M3DropCleanData(
  usoskin1,
  labels = colnames(usoskin1),
  min_detected_genes = 2000,
  is.counts = TRUE
)
expr_matrix <- uso_list$data # Normalized & filtered expression matrix
dim(expr_matrix)

celltype_labs <- uso_list$labels # filtered cell-type labels
cell_colors <- brewer.pal(max(3,length(unique(celltype_labs))), "Set3")

################################ 寻找highly variable genes (HVG) #################################
# 那些在样本群体里面表达量变异比较大的基因可能是真正的生物学现象，也有可能是技术误差，而且变异程度总是跟基因的表达量成正相关。
# 如下图所示：
plot(rowMeans(expr_matrix),rowVars(expr_matrix),log='xy')
# Brennecke et al.提出了算法来矫正这一影响，这个方法也被包装成了Brennecke_getVariableGenes(counts, spikes) 函数，
# 但是这个数据并没有ERCC spike-in，所以直接对整个表达矩阵处理即可。
Brennecke_HVG <- M3Drop::BrenneckeGetVariableGenes(
  expr_matrix,
  fdr = 0.01,
  minBiolDisp = 0.5
)


################################ 探究 High Dropout Genes ################################
# 另外一个寻找HVGs是查看它们是否有非常多的0表达量情况，这种多0表达的情况叫做dropout rate，通常单细胞转录组表达矩阵里面过半数的
# 基因都是0表达的。因为单细胞里面的mRNA很多无法被反转录，这种情况可以用Michaelis-Menten等式来模拟，如下图所示：
K = 49
S_sim = 10^seq(from=-3, to=4, by=0.05) # range of expression values
MM = 1-S_sim/(K+S_sim)
plot(S_sim, MM, type="l", lwd=3, xlab="Expression", ylab="Dropout Rate", xlim=c(1,1000))
S1 = 10; P1 = 1-S1/(K+S1) # Expression & dropouts for cells in condition 1
S2 = 750; P2 = 1-S2/(K+S2) # Expression & dropouts for cells in condition 2
points(c(S1,S2),c(P1,P2), pch=16, col="grey85", cex=3)
mix = 0.5; # proportion of cells in condition 1
points(S1*mix+S2*(1-mix), P1*mix+P2*(1-mix), pch=16, col="grey35", cex=3)
# 用来M3Drop包的M3DropFeatureSelection函数来挑选那些显著偏离了Michaelis-Menten曲线的基因，这里的阈值取1% FDR.
# 但是这个函数M3DropFeatureSelection依赖于正确的M3Drop包版本，下面就不运行了。
M3Drop_genes <- M3Drop::M3DropFeatureSelection(
  expr_matrix,
  mt_method = "fdr",
  mt_threshold = 0.01
)
title(main = "Usoskin")
M3Drop_genes <- M3Drop_genes$Gene

################################ Depth-Adjusted Negative Binomial (DANB) ################################
# 下面这个 NBumiConvertToInteger 也依赖于正确的M3Drop包版本，下面就不运行了
usoskin_int <- NBumiConvertToInteger(usoskin1)
DANB_fit <- NBumiFitModel(usoskin_int) # DANB is fit to the raw count matrix
# Perform DANB feature selection
DropFS <- NBumiFeatureSelectionCombinedDrop(DANB_fit)
DANB_genes <- names(DropFS[1:1500])

################################ 基因表达相关性 ################################
# 这个表达矩阵很大，所以计算所有基因之间的相关性耗时很长，为节约时间，也不运行了。
cor_mat <- cor(t(expr_matrix), method="spearman") #Gene-gene correlations
diag(cor_mat) <- rep(0, times=nrow(expr_matrix))
score <- apply(cor_mat, 1, function(x) {max(abs(x))}) #Correlation of highest magnitude
names(score) <- rownames(expr_matrix);
score <- score[order(-score)]
Cor_genes = names(score[1:1500])

################################ PCA挑选 ################################
# PCA 速度还行，挑选第1，2主成分的前1500个基因
pca <- prcomp(log(expr_matrix+1)/log(2)); 
# PCA is typically performed on log-transformed expression data
plot(pca$rotation[,1], pca$rotation[,2], pch=16, col=cell_colors[as.factor(celltype_labs)]) # plot projection
score <- rowSums(abs(pca$x[,c(1,2)])) 
# calculate loadings for components 1 and 2
names(score) <- rownames(expr_matrix)
score <- score[order(-score)]
PCA_genes = names(score[1:1500])


################################ 检查挑选的基因集的效果 ################################
# 热图+聚类可以看看基因是否在各个细胞类型差异表达，并且把细胞类型比较好的分开。这个热图非常耗时，如无必要，请不要运行这个代码
M3Drop::M3DropExpressionHeatmap(
  PCA_genes , ## 或者 M3Drop_genes 等其它方法挑到的基因
  uso_list$data,
  cell_labels = uso_list$labels
)

################################ 挑选的基因集跟DEseq得到的差异基因列表看交集 ################################
# 载入用DEseq得到的差异基因列，跟前面得到的M3Drop_genes比较一下。
# Load DE genes
# http://www.biotrainee.com/jmzeng/scRNA/DESeq_table.rds
DESeq_table <- readRDS("DESeq_table.rds")
DE_genes = unique(DESeq_table$Gene)

# Calculate precision
# sum(M3Drop_genes %in% DE_genes)/length(M3Drop_genes)
sum(PCA_genes %in% DE_genes)/length(PCA_genes)

#=============================================================================================
#                         比较不同的对单细胞转录组数据寻找差异基因的方法
#=============================================================================================
# 背景介绍
# 如果是bulk RNA-seq，那么现在最流行的就是DESeq2 和 edgeR啦，而且有很多经过了RT-qPCR 验证过的真实测序数据可以来评价不同的差异
# 基因算法的表现。对单细胞测序数据来说，通常需要先聚类之后把细胞群体进行分组，然后来比较不同的组的差异表达情况。当然，
# 也有不少单细胞测序实验设计本身就有时间点，不同个体来源，不同培养条件这样的分组！同时还有不少方法是不需要预先分类的，
# 因为分类本身就会引入偏差。跟bulk RNA-seq不一样的地方是，scRNA-seq通常涉及到的样本数量更多。这时候可以使用非参检验算法，
# 比如Kolmogorov-Smirnov test (KS-test) 等等。下面用一个测试数据来评价一下不同的算法的表现。处理同样的表达矩阵得到差异结果跟
# 已知的差异结果进行比较看看overlap怎么样。评价指标主要是：

# 召回率，True positive rate (TPR), TP/(TP + FN)
# 准确率，False positive rate (FPR), FP/(FP+TP)
# receiver-operating-characteristic curve (ROC)
# area under this curve (AUC)

library(ROCR)
library(edgeR)
library(DESeq2)

library(scater)
library(scde)
library(BPSC)
library(MAST)
library(monocle) 

# Differential Expression (在差异基因筛选这步优化的方法)
# Small number of cells and few groups : scde
# Replicates with batch effects : mixture/linear models
# Balanced batches: edgeR or MAST
# Large datasets: Kruskal-Wallis test (all groups at once), or Wilcox-test (compare 2-groups at a time).

# 加载测试数据
# 这里选取的是芝加哥大学Yoav Gilad lab实验的Tung et al 2017的单细胞测序文章的数据
# 我已经把测试数据保存为rdata数据格式了，直接加载。
setwd('/Users/mijiarui/Nature_Biotechnology_Paper/Testing_dataset')
load(file = 'scRNAseq_DEG_input.Rdata')
dim(counts)
dim(norm)
dim(DE)
dim(notDE)
table(group)

# 可以看到这里需要选择的测试数据来源于2个人，每个人都有288个细胞的表达数据。
# 就是要对它们进行差异比较，而已知的1083个基因是确定显著差异的，另外10897个基因是确定不显著的。(首先，
# 我们要假定这个是金标准！！！)
# 但是总共却是16026个基因，所以有一些基因是不确定显著与否的。

################################ Kolmogorov-Smirnov test ################################
# KS检验有两个弊端，首先是它假设基因表达量是连续的，如果有很多细胞表达量一致，比如都是0，表现就很差。
# 其次它对大样本量太敏感了，可能其实差异并不大，但是样本数量很多，也会被认为是显著差异。
pVals <- apply(norm, 1, function(x) {
  ks.test(x[group =="NA19101"], 
          x[group=="NA19239"])$p.value
})
# multiple testing correction
pVals <- p.adjust(pVals, method = "fdr")
sigDE <- names(pVals)[pVals < 0.05]
length(sigDE) # Number of KS-DE genes
sum(GroundTruth$DE %in% sigDE)  # Number of KS-DE genes that are true DE genes
sum(GroundTruth$notDE %in% sigDE) # Number of KS-DE genes that are not DE genes

tp <- sum(GroundTruth$DE %in% sigDE)
fp <- sum(GroundTruth$notDE %in% sigDE)
tn <- sum(GroundTruth$notDE %in% names(pVals)[pVals >= 0.05])
fn <- sum(GroundTruth$DE %in% names(pVals)[pVals >= 0.05])
tpr <- tp/(tp + fn)
fpr <- fp/(fp + tn)
cat(c(tpr, fpr))

ks_pVals=pVals

# 可以看到KS检验判断的显著差异基因实在是太多了，高达5095个。所以它能找回来792个真正的差异基因。
# 但是却找到了3190个假阳性。所以计算得到召回率73.46%，但是准确率只有29.44%，这个表现不佳。


# 再看看ROC和RUC
# Only consider genes for which we know the ground truth
pVals <- pVals[names(pVals) %in% GroundTruth$DE | 
                 names(pVals) %in% GroundTruth$notDE] 
truth <- rep(1, times = length(pVals));
truth[names(pVals) %in% GroundTruth$DE] = 0;
pred <- ROCR::prediction(pVals, truth)
perf <- ROCR::performance(pred, "tpr", "fpr")
ROCR::plot(perf)

aucObj <- ROCR::performance(pred, "auc")
aucObj@y.values[[1]] # AUC

# 把这两个评价分析包装成函数，后面可以直接使用！
DE_Quality_AUC <- function(pVals) {
  pVals <- pVals[names(pVals) %in% GroundTruth$DE | 
                   names(pVals) %in% GroundTruth$notDE]
  truth <- rep(1, times = length(pVals));
  truth[names(pVals) %in% GroundTruth$DE] = 0;
  pred <- ROCR::prediction(pVals, truth)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  ROCR::plot(perf)
  aucObj <- ROCR::performance(pred, "auc")
  return(aucObj@y.values[[1]])
}

DE_Quality_rate <- function(sigDE) {
  (length(sigDE) )
  # Number of KS-DE genes
  ( sum(GroundTruth$DE %in% sigDE) )
  # Number of KS-DE genes that are true DE genes
  (sum(GroundTruth$notDE %in% sigDE))
  tp <- sum(GroundTruth$DE %in% sigDE)
  fp <- sum(GroundTruth$notDE %in% sigDE)
  tn <- sum(GroundTruth$notDE %in% names(pVals)[pVals >= 0.05])
  fn <- sum(GroundTruth$DE %in% names(pVals)[pVals >= 0.05])
  tpr <- tp/(tp + fn)
  fpr <- fp/(fp + tn)
  cat(c(tpr, fpr))
}


##################################### Wilcox/Mann-Whitney-U Test #####################################
# 也是一种非参检验，通常比较两个组数据的median的差异。
pVals <- apply(norm, 1, function(x) {
  wilcox.test(x[group =="NA19101"], 
              x[group=="NA19239"])$p.value
})
# multiple testing correction
pVals <- p.adjust(pVals, method = "fdr")
sigDE <- names(pVals)[pVals < 0.05]
Wilcox_pVals=pVals

DE_Quality_rate(sigDE)

DE_Quality_AUC(pVals) 

# 召回率是81.9%，准确率是31.9%，这个表现不佳。

##################################### edgeR #####################################
# edgeR包在bulk RNA-seq测序领域应用很广泛，基于负二项分布模型，应用了 generalized linear model (GLM) 算法
library(edgeR)
dge <- DGEList(counts=counts, norm.factors = rep(1, length(counts[1,])), group=group)
group_edgeR <- factor(group)
design <- model.matrix(~group_edgeR)
dge <- estimateDisp(dge, design = design, trend.method="none")
fit <- glmFit(dge, design)
res <- glmLRT(fit)
pVals <- res$table[,4]
names(pVals) <- rownames(res$table)
pVals <- p.adjust(pVals, method = "fdr")
sigDE <- names(pVals)[pVals < 0.05]
edgeR_pVals=pVals

DE_Quality_rate(sigDE)

DE_Quality_AUC(pVals)

# 召回率是86.9%，准确率是39.4%，表现越来越好了。

##################################### Monocle #####################################
library(monocle)
pd <- data.frame(group=group, batch=batch)
rownames(pd) <- colnames(counts)
pd <- new("AnnotatedDataFrame", data = pd)
## 针对于基于read counts的表达矩阵
Obj <- newCellDataSet(as.matrix(counts), phenoData=pd, 
                      expressionFamily=negbinomial.size()) 
Obj <- estimateSizeFactors(Obj)
Obj <- estimateDispersions(Obj)
res <- differentialGeneTest(Obj,fullModelFormulaStr="~group")
pVals <- res[,3]
names(pVals) <- rownames(res)
pVals <- p.adjust(pVals, method = "fdr")
sigDE <- names(pVals)[pVals < 0.05]
monocle_pVals=pVals

DE_Quality_rate(sigDE)

DE_Quality_AUC(pVals)

# monocle做差异分析的耗时非常夸张，召回率是84.2%，准确率是38.1%

##################################### MAST #####################################
# MAST基于 zero-inflated negative binomial 分布模型
library(MAST)
log_counts <- log(counts+1)/log(2)
fData = data.frame(names=rownames(log_counts))
rownames(fData) = rownames(log_counts);
cData = data.frame(cond=group)
rownames(cData) = colnames(log_counts)
obj <- FromMatrix(as.matrix(log_counts), cData, fData)
colData(obj)$cngeneson <- scale(colSums(assay(obj)>0))
cond <- factor(colData(obj)$cond)
# Model expression as function of condition & number of detected genes
zlmCond <- zlm.SingleCellAssay(~cond + cngeneson, obj) 
summaryCond <- summary(zlmCond, doLRT="condNA19101")
summaryDt <- summaryCond$datatable
summaryDt <- as.data.frame(summaryDt)
pVals <- unlist(summaryDt[summaryDt$component == "H",4]) # H = hurdle model
names(pVals) <- unlist(summaryDt[summaryDt$component == "H",1])
pVals <- p.adjust(pVals, method = "fdr")
sigDE <- names(pVals)[pVals < 0.05]
MAST_pVals=pVals

DE_Quality_rate(sigDE)

DE_Quality_AUC(pVals)

# 召回率是82.8%，准确率是34.9.%


##################################### BPSC #####################################
# 这个用的是 Poisson-Beta 分布模型 (非常耗时)
library(BPSC)
bpsc_data <- norm[,batch=="NA19101.r1" | batch=="NA19239.r1"]
bpsc_group = group[batch=="NA19101.r1" | batch=="NA19239.r1"]
control_cells <- which(bpsc_group == "NA19101")
design <- model.matrix(~bpsc_group)
coef=2 # group label
res=BPglm(data=bpsc_data, controlIds=control_cells, design=design, coef=coef, 
          estIntPar=FALSE, useParallel = FALSE)
pVals = res$PVAL
pVals <- p.adjust(pVals, method = "fdr")
sigDE <- names(pVals)[pVals < 0.05]
BPSC_pVals=pVals

DE_Quality_rate(sigDE)

DE_Quality_AUC(pVals)

# 召回率是64.8%，准确率是30.7.%


##################################### SCDE #####################################
# SCDE是第一个特意针对单细胞转录组测序数据的差异分析而设计的，用贝叶斯统计方法把表达矩阵拟合到 
# zero-inflated negative binomial 分布模型里面。
library(scde)
cnts <- apply(
  counts,
  2,
  function(x) {
    storage.mode(x) <- 'integer'
    return(x)
  }
)
names(group) <- 1:length(group)
colnames(cnts) <- 1:length(group)
o.ifm <- scde::scde.error.models(
  counts = cnts,
  groups = group,
  n.cores = 1,
  threshold.segmentation = TRUE,
  save.crossfit.plots = FALSE,
  save.model.plots = FALSE,
  verbose = 0,
  min.size.entries = 2
)
priors <- scde::scde.expression.prior(
  models = o.ifm,
  counts = cnts,
  length.out = 400,
  show.plot = FALSE
)
resSCDE <- scde::scde.expression.difference(
  o.ifm,
  cnts,
  priors,
  groups = group,
  n.randomizations = 100,
  n.cores = 1,
  verbose = 0
)
# Convert Z-scores into 2-tailed p-values
pVals <- pnorm(abs(resSCDE$cZ), lower.tail = FALSE) * 2
pVals <- p.adjust(pVals, method = "fdr")
sigDE <- names(pVals)[pVals < 0.05]
SCDE_pVals=pVals

DE_Quality_rate(sigDE)

DE_Quality_AUC(pVals)

save(SCDE_pVals,BPSC_pVals,MAST_pVals,monocle_pVals,edgeR_pVals,Wilcox_pVals,ks_pVals,file = 'DEG_results.Rdata')

#=============================================================================================
#                                         统计学基础
#=============================================================================================

######################### 负二项分布 negative binomial model #########################
# 这个是被应用的最广泛的转录组表达数据分布模型。但是对单细胞转录组测序数据来说，因为有很高的dropout情况，
# 导致模型失准，所以就提出来了zero-inflated negative binomial models
set.seed(1)
hist(rnbinom(1000, mu=10, size=100), col="grey50", xlab="Read Counts", main="Negative Binomial")

######################### zero-inflated negative binomial models #########################
# 就是在原始的负二项分布数据里面随机挑选一些低表达量基因，给它们人为赋值为0表达量值。
d = 0.5;
counts <- rnbinom(1000, mu=10, size=100);
counts[runif(1000) < d] = 0;
hist(counts, col="grey50", xlab="Read Counts", main="Zero-inflated NB")

######################### Poisson-Beta distribution #########################
a = 0.1
b = 0.1
g = 100
lambdas = rbeta(1000, a, b)
counts = sapply(g*lambdas, function(l) {rpois(1, lambda=l)})
hist(counts, col="grey50", xlab="Read Counts", main="Poisson-Beta")

#=============================================================================================
#                          比较不同的对单细胞转录组数据聚类的方法
#=============================================================================================
# 背景介绍
# 聚类之前必须要对表达矩阵进行normalization，而且要去除一些批次效应等外部因素。通过对表达矩阵的聚类，可以把细胞群体分成不同的状态，
# 解释为什么会有不同的群体。不过从计算的角度来说，聚类还是蛮复杂的，各个细胞并没有预先标记好，而且也没办法事先知道可以聚多少类。
# 尤其是在单细胞转录组数据里面有很高的噪音，基因非常多，意味着的维度很高。对这样的高维数据，需要首先进行降维，可以选择PCA或者t-SNE
# 方法。聚类的话，一般都是无监督聚类方法，比如：hierarchical clustering, k-means clustering and graph-based clustering。
# 算法略微有一点复杂，略过吧。

# 这里主要比较6个常见的单细胞转录组数据的聚类包：

# SINCERA
# pcaReduce
# SC3
# tSNE + k-means
# SEURAT
# SNN-Cliq

# 加载代码如下：
library(pcaMethods)
library(pcaReduce)
library(SC3)
library(scater)
library(pheatmap)
set.seed(1234567)

################################ 加载测试数据 ################################

# 这里选取的是数据，加载了这个scater包的SCESet对象，包含着一个22431 features, 268 samples 的表达矩阵。供10已知的种细胞类型，
# 这样聚类的时候就可以跟这个已知信息做对比，看看聚类效果如何。可以直接用plotPCA来简单PCA并且可视化。
deng <- readRDS("deng-reads.rds")
deng
dim(counts(deng ))
head(colData(deng))
head(rowData(deng))
table(colData(deng)$cell_type2)
plotPCA(deng, colour_by = "cell_type2")
### 可以看到简单的PCA也是可以区分部分细胞类型的，只不过在某些细胞相似性很高的群体区分力度不够，
### 所以需要开发新的算法来解决这个聚类的问题。


#####################################  SC聚类  ##################################### 
deng <- sc3_estimate_k(deng)
metadata(deng)$sc3$k_estimation
plotPCA(deng, colour_by = "cell_type1")
## 准备 SCESet对象 数据给 SC3方法，先预测能聚多少个类。
## 这里是并行计算，所以速度还可以
deng <- sc3(deng, ks = 10, biology = TRUE)
deng
head(rowData(deng))
## 可以看到SC3方法处理后的SCESet对象的基因信息增加了5列，比较重要的是sc3_gene_filter信息，决定着该基因是否拿去聚类，
## 因为基因太多了，需要挑选
table(rowData(deng)$sc3_gene_filter)
### 只有一半的基因被挑选去聚类了

## 后面是一些可视化
sc3_plot_consensus(deng, k = 10, show_pdata = "cell_type2")
sc3_plot_silhouette(deng, k = 10)
sc3_plot_expression(deng, k = 10, show_pdata = "cell_type2")
sc3_plot_markers(deng, k = 10, show_pdata = "cell_type2")
plotPCA(deng, colour_by = "sc3_10_clusters")

## 还支持shiny的交互式聚类，暂时不显示
sc3_interactive(deng)
# 很明显可以看到SC3聚类的效果要好于普通的PCA


#####################################  pcaReduce  #####################################
# use the same gene filter as in SC3
input <- logcounts(deng[rowData(deng)$sc3_gene_filter, ])
# run pcaReduce 1 time creating hierarchies from 1 to 30 clusters
pca.red <- PCAreduce(t(input), nbt = 1, q = 30, method = 'S')[[1]]
##  这里对2~30种类别的情况都分别对样本进行分组。
## 我们这里取只有10组的时候，这些样本是如何分组的信息来可视化。
colData(deng)$pcaReduce <- as.character(pca.red[,32 - 10])
plotPCA(deng, colour_by = "pcaReduce")


#####################################  tSNE + kmeans  #####################################
# scater包包装了 Rtsne 和 ggplot2 来做tSNE并且可视化。
deng <- plotTSNE(deng, rand_seed = 1, return_SCE = TRUE)
## 上面的tSNE的结果，下面用kmeans的方法进行聚类，假定是8类细胞类型。
colData(deng)$tSNE_kmeans <- as.character(kmeans(deng@reducedDims$TSNE, centers = 8)$clust)
plotTSNE(deng, rand_seed = 1, colour_by = "tSNE_kmeans")


#####################################  SINCERA  #####################################  
# perform gene-by-gene per-sample z-score transformation
dat <- apply(input, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
# hierarchical clustering
dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)
hc <- hclust(dd, method = "average")
num.singleton <- 0
kk <- 1
for (i in 2:dim(dat)[2]) {
  clusters <- cutree(hc, k = i)
  clustersizes <- as.data.frame(table(clusters))
  singleton.clusters <- which(clustersizes$Freq < 2)
  if (length(singleton.clusters) <= num.singleton) {
    kk <- i
  } else {
    break;
  }
}
cat(kk)
pheatmap(
  t(dat),
  cluster_cols = hc,
  cutree_cols = kk,
  kmeans_k = 100,
  show_rownames = FALSE
)


#=============================================================================================
#                                      "Scater" package
#=============================================================================================
# scater 这个R包很强大，是McCarthy et al. 2017 发表的，包含的功能有：
# Automated computation of QC metrics
# Transcript quantification from read data with pseudo-alignment
# Data format standardisation
# Rich visualizations for exploratory analysis
# Seamless integration into the Bioconductor universe
# Simple normalisation methods

# 请深刻理解Scater的工作流程
# S4对象

# 主要是基于 SCESet 对象来进行下游分析，跟ExpressionSet对象类似，也是常见的3个组成：
##  exprs, a numeric matrix of expression values, where rows are features, and columns are cells

## phenoData, an AnnotatedDataFrame object, where rows are cells, and columns are cell attributes 
## (such as cell type, culture condition, day captured, etc.)

## featureData, an AnnotatedDataFrame object, where rows are features (e.g. genes), and columns are feature attributes, 
## such as biotype, gc content, etc.

## 主要就是读取scRNA上游分析处理得到的表达矩阵，加上每个样本的描述信息，形成矩阵之后。对样本进行过滤，然后对基因进行过滤。
## 针对过滤后的表达矩阵进行各种分类的可视化。

##################################### 加载测试数据和构建SCESet对象 #####################################
suppressPackageStartupMessages(library(scater))
data("sc_example_counts")
data("sc_example_cell_info") 

example_sce <- SingleCellExperiment(
  assays = list(counts = sc_example_counts), 
  colData = sc_example_cell_info
)

## 初识SCESet对象
example_sce
counts(example_sce)
rowData(example_sce)
colData(example_sce)

### 增加slot
exprs(example_sce) <- log2(calculateCPM(example_sce, use_size_factors = FALSE) + 1)
example_sce # 在assay slot增加了logcounts这个slot

### 对feature进行过滤
keep_feature <- rowSums(exprs(example_sce) > 0) > 0
example_sce <- example_sce[keep_feature,]
example_sce # 从原来2000个gene过滤到只剩1973个gene

example_sce <- calculateQCMetrics(example_sce, 
                                  feature_controls = list(eg = 1:40))

example_sce

##################################### shiny和可视化 #####################################

## Scater的shiny真的非常好用，所有的可视化都集中在了 scater_gui 这个函数产生的shiny网页里面：
## plotScater: a plot method exists for SingleCellExperiment objects, which gives an overview of expression across cells.
## plotQC: various methods are available for producing QC diagnostic plots.
## plotPCA: produce a principal components plot for the cells.
## plotTSNE: produce a t-distributed stochastic neighbour embedding (reduced dimension) plot for the cells.
## plotDiffusionMap: produce a diffusion map (reduced dimension) plot for the cells.
## plotMDS: produce a multi-dimensional scaling plot for the cells.
## plotReducedDim: plot a reduced-dimension representation of the cells.
## plotExpression: plot expression levels for a defined set of features.
## plotPlatePosition: plot cells in their position on a plate, coloured by cell metadata and QC metrics or feature expression level.
## plotColData: plot cell metadata and QC metrics.
## plotRowData: plot feature metadata and QC metrics.
scater_gui(example_sce)
library(destiny)

##可以充分的探索自己的数据，随便看一个可视化函数的结果：
plotExpression(example_sce, rownames(example_sce)[1:6],
               x = "Mutation_Status", exprs_values = "exprs", 
               colour = "Treatment")
colData(example_sce)


##################################### 详细的QC #####################################
# 做QC要结合上面的可视化步骤，所以没办法自动化，只能先可视化，肉眼分辨一下哪些样本或者基因数据是需要舍弃的。
## 初始化一些可视化参数
library(knitr)
opts_chunk$set(fig.align = 'center', fig.width = 6, fig.height = 5, dev = 'png')
library(ggplot2)
theme_set(theme_bw(12))

## quickstart -- loading data, 加载数据
suppressPackageStartupMessages(library(scater))
data("sc_example_counts")
data("sc_example_cell_info")

## quickstart -- make sce, 构建sce对象
gene_df <- DataFrame(Gene = rownames(sc_example_counts))
rownames(gene_df) <- gene_df$Gene
example_sce <- SingleCellExperiment(assays = list(counts = sc_example_counts), 
                                    colData = sc_example_cell_info, 
                                    rowData = gene_df)  # 这次相当于强行给feature data加了一列

example_sce <- normalise(example_sce)

## quickstart -- add-exprs
exprs(example_sce) <- log2(calculateCPM(example_sce, use_size_factors = FALSE) + 1)

## filter no-exprs
keep_feature <- rowSums(exprs(example_sce) > 0) > 0
example_sce <- example_sce[keep_feature,]
example_sceset <- calculateQCMetrics(example_sce, feature_controls = list(eg = 1:40)) 
colnames(colData(example_sceset))

## 首先是基于样本的过滤，用 colData(object) 可以查看各个样本统计情况
## 然后是基于基因的过滤，用 rowData(object) 可以查看各个基因统计情况

##################################### scater一站式过滤低质量样本 #####################################
# scater包自己提供了一个基于PCA的QC标准，不需要自己根据文库大小，覆盖的基因数量，外源的ERCC spike-ins 
# 含量以及线粒体DNA含量来进行人工过滤。
# 默认的筛选条件如下：
## pct_counts_top100features
## total_features
## pct_counts_feature_controls
## n_detected_feature_controls
## log10_counts_endogenous_features
## log10_counts_feature_controls

#一站式QC函数如下：
dat_pca <- scater::plotPCA(dat_qc,
                           size_by = "total_features", 
                           shape_by = "use",
                           pca_data_input = "pdata",
                           detect_outliers = TRUE,
                           return_SCESet = TRUE)

# 过滤只是它最基本的工具，它作为单细胞转录组3大R包，功能肯定是非常全面的，比如前面我们讲解的normalization，DEG, 
# features selection，cluster，它都手到擒来，只不过是包装的是其它R包的函数。



#=============================================================================================
#                                   Pseudotime analysis
#=============================================================================================
setwd('/Users/mijiarui/Nature_Biotechnology_Paper/Testing_dataset')
library(ggbeeswarm)
library(ggthemes)
# Let us take a first look at the Deng data, without yet applying sophisticated pseudotime methods. As the plot below shows, 
# simple PCA does a very good job of displaying the structure in these data. It is only once we reach the blast cell types 
# (“earlyblast”, “midblast”, “lateblast”) that PCA struggles to separate the distinct cell types.

deng_SCE <- readRDS("deng-reads.rds")
deng_SCE$cell_type2 <- factor(
  deng_SCE$cell_type2,
  levels = c("zy", "early2cell", "mid2cell", "late2cell",
             "4cell", "8cell", "16cell", "earlyblast",
             "midblast", "lateblast")
)
cellLabels <- deng_SCE$cell_type2
deng <- counts(deng_SCE)
colnames(deng) <- cellLabels
deng_SCE <- runPCA(deng_SCE, colour_by = "cell_type2", 
                   return_SCE = TRUE)

## PCA, here, provides a useful baseline for assessing different pseudotime methods. For a very naive pseudotime we can 
## just take the co-ordinates of the first principal component.
deng_SCE$PC1 <- reducedDim(deng_SCE, "PCA")[,1]
ggplot(as.data.frame(colData(deng_SCE)), aes(x = PC1, y = cell_type2, 
                                             colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("First principal component") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")

##################################### Monocle #####################################
# Monocle skips the clustering stage of TSCAN and directly builds a minimum spanning tree on a reduced dimension representation 
# of the cells to connect all cells. Monocle then identifies the longest path in this tree to determine pseudotime. If the data 
# contains diverging trajectories (i.e. one cell type differentiates into two different cell-types), monocle can identify these. 
# Each of the resulting forked paths is defined as a separate cell state.

# Unfortunately, Monocle does not work when all the genes are used, so we must carry out feature selection. First, we use M3Drop:
m3dGenes <- as.character(M3DropFeatureSelection(deng)$Gene)
d <- deng[which(rownames(deng) %in% m3dGenes), ]
d <- d[!duplicated(rownames(d)), ]

# Now run monocle
colnames(d) <- 1:ncol(d)
geneNames <- rownames(d)
rownames(d) <- 1:nrow(d)
pd <- data.frame(timepoint = cellLabels)
pd <- new("AnnotatedDataFrame", data=pd)
fd <- data.frame(gene_short_name = geneNames)
fd <- new("AnnotatedDataFrame", data=fd)

dCellData <- newCellDataSet(d, phenoData = pd, featureData = fd, expressionFamily = tobit())
dCellData <- setOrderingFilter(dCellData, which(geneNames %in% m3dGenes))
dCellData <- estimateSizeFactors(dCellData)
dCellDataSet <- reduceDimension(dCellData, pseudo_expr = 1)
dCellDataSet <- orderCells(dCellDataSet, reverse = FALSE)
plot_cell_trajectory(dCellDataSet)



##################################### TSCAN #####################################
# TSCAN combines clustering with pseudotime analysis. First it clusters the cells using mclust, which is based on a mixture of 
# normal distributions. Then it builds a minimum spanning tree to connect the clusters. The branch of this tree that connects 
# the largest number of clusters is the main branch which is used to determine pseudotime.
library(TSCAN)

procdeng <- TSCAN::preprocess(deng)
colnames(procdeng) <- 1:ncol(deng)
dengclust <- TSCAN::exprmclust(procdeng, clusternum = 10)
TSCAN::plotmclust(dengclust)

dengorderTSCAN <- TSCAN::TSCANorder(dengclust, orderonly = FALSE)
pseudotime_order_tscan <- as.character(dengorderTSCAN$sample_name)
deng_SCE$pseudotime_order_tscan <- NA
deng_SCE$pseudotime_order_tscan[as.numeric(dengorderTSCAN$sample_name)] <- 
  dengorderTSCAN$Pseudotime

# Frustratingly, TSCAN only provides pseudotime values for 221 of 268 cells, silently returning missing values for 
# non-assigned cells.
cellLabels[dengclust$clusterid == 10]

ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_order_tscan, 
           y = cell_type2, colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("TSCAN pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by TSCAN pseudotime")


##################################### Diffusion maps #####################################
# Diffusion maps were introduced by Ronald Coifman and Stephane Lafon, and the underlying idea is to assume that the data 
# are samples from a diffusion process. The method infers the low-dimensional manifold by estimating the eigenvalues and 
# eigenvectors for the diffusion operator related to the data.

# Haghverdi et al have applied the diffusion maps concept to the analysis of single-cell RNA-seq data to create an R package 
# called destiny.

# We will take the rank order of cells in the first diffusion map component as “diffusion map pseudotime” here.
deng <- logcounts(deng_SCE)
colnames(deng) <- cellLabels
dm <- DiffusionMap(t(deng))

tmp <- data.frame(DC1 = eigenvectors(dm)[,1],
                  DC2 = eigenvectors(dm)[,2],
                  Timepoint = deng_SCE$cell_type2)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
  geom_point() + scale_color_tableau() + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()

# Like the other methods, using the first diffusion map component from destiny as pseudotime does a good job at ordering the 
# early time-points (if we take high values as “earlier” in developement), but it is unable to distinguish the later ones.
deng_SCE$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])
ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_diffusionmap, 
           y = cell_type2, colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("Diffusion map pseudotime (first diffusion map component)") +
  ylab("Timepoint") +
  ggtitle("Cells ordered by diffusion map pseudotime")


#=============================================================================================
#                         在单细胞转录组表达矩阵里面去除细胞周期影响
#=============================================================================================
# 早在2015年发表在Nat. Biotechnol文章就提出了 scLVM (single-cell latent variable model)来在单细胞转录组数据
# 里面去除细胞周期影响 但是 scLVM 仅仅是考虑 细胞周期直接相关基因，而且没有考虑细胞类型，其实不同类型的细胞
# 哪怕是在同一个时间点的细胞周期状态，它们的细胞周期相关基因表达也是不同的。更重要的是，还有很多非直接细胞
# 周期相关基因也需要考虑。

# 所以作者开发了ccRemover来去除单细胞转录组数据里面去除细胞周期影响，发表于2016年的SR
# 核心代码就一句话：可以看到主要就是使用ccRemover函数处理我们的单细胞表达矩阵哦

# xhat <- ccRemover(dat, bar=FALSE)
# xhat <- ccRemover(dat, cutoff = 3, max_it = 4, nboot = 200, ntop = 10, bar=FALSE)
# 可以简单处理，也可以根据理解，加上一系列的参数。 表达矩阵需要是归一化的。 下面我们就具体讲解。

##################################### 表达矩阵的前期处理 #####################################
# 首先安装并且加载包和测试数据
library(ccRemover)
data(t.cell_data)
## 表达矩阵如下；
head(t.cell_data[,1:5])

## 基因表达量在样本的平均值的汇总
summary(apply(t.cell_data,1, mean))

## 样本的所有的基因的表达量的平均值的汇总
summary(apply(t.cell_data,2, mean))

## 保留每个基因的所有样本表达量平均值以备后用
mean_gene_exp <- rowMeans(t.cell_data)
# 每个基因减去其平均值后的表达量矩阵
t_cell_data_cen <- t.cell_data - mean_gene_exp
## 接近于 0 了
summary(apply(t_cell_data_cen,1,mean))

gene_names <- rownames(t_cell_data_cen)
head(gene_names)

## gene_indexer函数会根据包内置的细胞周期相关基因来判断我们的表达矩阵的基因是否属于
cell_cycle_gene_indices <- gene_indexer(gene_names, species = "mouse", 
                                        name_type = "symbols" )
## Invalid name type input. Switching to NULLNo name format input.
## Checking to see if match can be found:
## Best guess is  symbol  IDs
## 751  matches out of a possible  7073
if_cc <- rep(FALSE,nrow(t_cell_data_cen)) 
if_cc[cell_cycle_gene_indices] <- TRUE
summary(if_cc)

## 构造表达矩阵以及其对应的基因是否属于细胞周期相关基因集
dat <- list(x=t_cell_data_cen, if_cc=if_cc)
xhat <- ccRemover(dat, bar=FALSE)

xhat <- xhat + mean_gene_exp
# 最后可以把基因的平均值加回去

# 高级参数
xhat <- ccRemover(dat, cutoff = 3, max_it = 4, nboot = 200, ntop = 10, bar=FALSE)

# 建议直接看英文
# The ‘cutoff’ is used to determine which of the effects are cell-cycle effects. The default and recommended 
# value is 3, which roughly corresponds to a p-value of 0.01. For data sets which have very low levels of 
# cell-cycle activity this value can be lowered to increase the detection of cell-cycle effects. 
# Example 3 in the original manuscript was a case where a lower value of the cutoff was necessary.

# The ‘max.it’ value is the maximum number of iterations of the method. ccRemover will stop whenever it detects 
# no more significant effects present in the data or it reaches its maximum number of iterations. 
# The default value is 4 but we have found that for many data sets the cell-cycle effect will be effectively 
# removed after 1 or 2 iterations.

# The ‘nboot’ value corresponds to the number of bootstrap repetitions carried out to test the significance 
# of the components. Please refer to the methods section original manuscript for a detailed description of 
# the estimation process. The default value is 200 and we have found this work effectively for most data sets.

# The ‘ntop’ parameter determines the number of principal components which are to be tested as cell-cycle 
# effects upon each iteration. We have found that the default value of 10 works effectively for most data sets. 
# However, for extremely diverse data sets within which there are expected to be many sources of variation 
# this value could be raised so that all elements of the cell-cycle effect can be identified and removed.


#=============================================================================================
#                         了解scRNAseq包内置的单细胞转录组数据
#=============================================================================================

library(scRNAseq)
## ----- Load Example Data -----
data(fluidigm)
# Set assay to RSEM estimated counts
assay(fluidigm) = assays(fluidigm)$rsem_counts
# List all qc fields (accessible via colData())
metadata(fluidigm)$which_qc
## 可以看到，默认可以做QC的项目非常多。
## 值得注意的是 RALIGN 代表比对率， NREADS 代表文库大小。


# Joint distribution of "biological condition"" and "coverage type""
table(colData(fluidigm)$Coverage_Type,
      colData(fluidigm)$Biological_Condition)

# 这里面的表达矩阵是由 RSEM (Li and Dewey 2011) 软件根据 hg38 RefSeq transcriptome 得到的，总是130个文库，
# 每个细胞测了两次，测序深度不一样。

########### 可以用scran包的cyclone函数进行单细胞转录组的细胞周期状态推断 ###########
# 代码如下，需要注意是针对人类的测序数据，然后还有一些ID的转换。
library(scRNAseq)
library(RColorBrewer)
library(scran)
library(scater)
library(org.Hs.eg.db)
all.counts=assays(fluidigm)$rsem_counts
dim(all.counts)
head(rownames(all.counts))
sce <- SingleCellExperiment(list(counts=all.counts)) 
ensembl <- mapIds(org.Hs.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")
mm.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
assigned <- cyclone(sce, pairs=mm.pairs, gene.names=ensembl)
head(assigned$scores)
plot(assigned$score$G1, assigned$score$G2M, 
     xlab="G1 score", ylab="G2/M score", pch=16)

table(assigned$phases)
pheatmap::pheatmap(t(assigned$scores))







