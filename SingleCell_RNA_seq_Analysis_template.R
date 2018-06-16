#=============================================================================================
#
#            single-cell RNA-seq analysis using different packages (methods)
#
#=============================================================================================
# 现在R里处理single-cell RNA sequencing的packages主要分成两派：SingleCellExperiment派和Seurat派。
# SingleCellExperiment派使用SingleCellExperiment R包里定义的S4 object，主打分析包有s cater (scater，数据QC和初步visualisation)，SC3（SC3，聚类分析）……（待补充）
# Seurat派使用Seurat包里定义的S4 object储存数据，主打分析包有Seurat（Seurat，聚类分析等全套内容），monocle（Monocle，pseudotime analysis等）……（待补充）
# 有一个R package可以将Seurat S4 object和SingleCellExperiment S4 object相互转换（这个名字现在想不起来了，在推特上看过，待查证）。
#=============================================================================================
#                                      "Seurat" package
#=============================================================================================
# 牛津大学的Rahul Satija等开发的Seurat，最早公布在Nature biotechnology, 2015，文章是； Spatial reconstruction of single-cell 
# gene expression data , 在2017年进行了非常大的改动，所以重新在biorxiv发表了文章在 Integrated analysis of single cell transcriptomic 
# data across conditions, technologies, and species 。 功能涵盖了scRNA-seq的QC、过滤、标准化、批次效应、PCA、tNSE亚群聚类分析、
# 差异基因分析、亚群特异性标志物鉴定等等等。给初学者提供了一个2,700 PBMC scRNA-seq dataset from 10X genomics的数据实战指导；这里的测试
# 数据是经由Illumina NextSeq 500测到的2,700 single cells 表达矩阵。

# 新版的Seurat把对象的建立分解成7个步骤
# 1， 导入Expression matrix
# 2,  对象的初始化，使用CreateSeuratObject函数，结果会存储在object@raw.data中，一般输入进去的对象是下机处理后的原始矩阵数据
#     可以是counts，或是某些软件校正后的数值，比如FPKM，TPM等等，但是绝对不能是log转换后的数据；Seurat可以倒入任何方法得到的
#     矩阵数据，其中Read10X是专门指读入10X的数据；我们更建议使用counts data，因为我们可以使用gene scaling和差异基因分析，并
#     服从负二项分布
# 3， 对细胞进行过滤，select/remove cells based on QC metrics
# 4,  对数据进行归一化normalization，结果会存储在object@data中，是经过类似CPM然后再log转换后的数据，一般这部分数据是用来进行
#     数据可视化的，比如violinplot或者featureplot，并进行差异基因分析，HVG的寻找等等，并且作为ScaleData的输入。
# 5， 寻找差异基因，结果会存储在object@hvg.info和object@var.gene中
# 6， 对结果进行scaling，目的是scale和回归去除混杂因素，结果存储在object@scaled.data中；经过了对测序文库大小的校正，因为根据0
#     进行了scale，所以结果有正值也有负值。如果校正了数据变异的特定的来源，比如cell-cycle，那么scaled residuals from the model
#     也会被存储进来。只有经过scale的数据才能进行数据降维的操作。
# 7， 得到最终当然Seurat对象

# Seurat的slot比较复杂，建议?Seurat来查看slot的帮助文档


############################################### 载入数据 ###############################################

### Load packages，加载数据前需要将文件夹中的三个文件分别命名为“matrix.mtx", "barcodes.tsv", "genes.tsv"
### 注意，许多scRNA-seq的pipeline返回的已经是一个稀疏矩阵了(sparse matrix)，比如CellRanger返回的mtx格式文件
library(Seurat) # the seurat version now I am using is 2.3.1
library(dplyr)
library(Matrix)
library(ggplot2)



### set working directory
setwd('/Users/mijiarui/Nature_Biotechnology_Paper/Pancreas_3_endo_mRNA_GSM3032164_P7endo')
### 根据表达矩阵构建seurat对象：加载数据前需要将文件夹中的三个文件分别命名为“matrix.mtx", "barcodes.tsv", "genes.tsv"，需要准备好3个输入文件
list.files("/Users/mijiarui/Nature_Biotechnology_Paper/Pancreas_3_endo_mRNA_GSM3032164_P7endo") # 通过这个函数查看读入文件夹内包含的文件
### 这三个文件，必须按照上述命名要求命名，后面Read10X()会自动识别。

### 做完这一步，建议回到对一个文件夹，进入终端(terminal)，使用head命令查看matrix.mtx
#### matrix.mtx is a MatrixMarket file; see http://math.nist.gov/MatrixMarket/formats.html
#### only non-zero entries are stored in the file 这种格式的矩阵文件只记录非0的entry
#### comments start with a %, like LaTeX
#### the first line indicates the total number of rows, columns, and entries 第一行显示行数和列数，以及entry数，和dim()是一样的
#### the following lines provide a row and column number and the value at that coordinate 之后的每一行给出的是对应行和对应列以及对应的表达量

### Load the Pancreas_1_mRNA_GSM2830058_P5 dataset
pancreas_1.data <- Read10X(data.dir = "/Users/mijiarui/Nature_Biotechnology_Paper/Pancreas_3_endo_mRNA_GSM3032164_P7endo")


### 载入数据后一定注意查看下object!
pancreas_1.data[1:6,1:6] # 注意这一步仅仅是载入数据，还没有构建Seurat对象; 上面一步相当于就是把表达矩阵读进来了，metadata还没用整合进来
class(pancreas_1.data)    # https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/dgTMatrix-class.html
dim(pancreas_1.data)

str(pancreas_1.data) # 对于复杂的组学数据，使用str()函数可以更清楚的了解数据的组成，并且方便后续数据操作
pancreas_1.data@Dim    # 查看一下矩阵的维度。斑马鱼dataset注释的基因是28283个
length(pancreas_1.data@Dimnames[[1]]) # 每一个gene都是用gene_symbol来表示的，Dimnames下的第一个对象是基因名
length(pancreas_1.data@Dimnames[[2]]) # 每一个细胞，我们用对应的barcode来表示，Dimnames下的第二个对象是cell barcode

summary(colSums(pancreas_1.data)) # summary of total expression per single cell
# check how many genes have at least one transcript in each cell
at_least_one <- apply(pancreas_1.data, 2, function(x) sum(x>0))
hist(at_least_one, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")

hist(colSums(pancreas_1.data),
     breaks = 100, main = "Expression sum per cell",
     xlab = "Sum expression")

# 我们后续要进行数据过滤，保留基因（比如至少在3个细胞当中表达的gene）和细胞（至少有200个被探测到的细胞）
# manually check the number of genes detected in three or more cells
# a lot of genes are not detected in 3 or more cells
tmp <- apply(pancreas_1.data, 1, function(x) sum(x>0))
table(tmp>=3)

# all cells have at least 200 detected genes
keep <- tmp>=3
tmp <- pancreas_1.data[keep,]
at_least_one <- apply(tmp, 2, function(x) sum(x>0)>=200)
summary(at_least_one)
tmp <- tmp[,at_least_one]
dim(tmp)

### Examine the memory savings between regular and sparse matrices
### 这一步是可选项，分别查看密集矩阵和稀疏矩阵所占用的空间，大致可以判断矩阵当中0的数量，使用object.size()函数
dense.size <- object.size(x = as.matrix(x = pancreas_1.data))
dense.size # 密集矩阵所占用的空间

sparse.size <- object.size(x = pancreas_1.data)
sparse.size # 稀疏矩阵所占用的空间

dense.size / sparse.size # 稀疏矩阵比密集矩阵压缩了近20倍的空间(在这个数据集中)

############################################### 构建Seurat对象 ###############################################
# Initialize the Seurat object with the raw (non-normalized data).  Keep all genes expressed in >= 3 cells (~0.1% of the data). 
# Keep all cells with at least 200 detected genes
# 注意Seurat设定了自己的object，就叫做Seurat object
pancreas_1 <- CreateSeuratObject(raw.data = pancreas_1.data, min.cells = 3, min.genes = 200,  # 注意过滤条件
                                 project = "10X_Pancreas_1")
pancreas_1 # 构建对象后常规查看一下对象的内容，我们看到一部分基因和细胞被过滤掉了；统计结果和上面人工计算的结果稍有出入
class(pancreas_1)
slotNames(pancreas_1)
str(pancreas_1)  # 还是强烈建议用str()函数来查看完整信息，pancreas_1@raw.data就相当于Read10X读入的pancreas_1.data
rownames(pancreas_1@raw.data); rownames(pancreas_1@data) # 这两句代码都可以用来查看gene_symbol
## Seurat的tutorial指出Seurat会自动帮我们计算nGene和nUMI。
## nUMI的计算方法是num.mol <- colSums(object.raw.data),每一个转录本对应一个UMI
## nGene的计算方法是num.genes <- colSums(object.raw.data > is.expr) ，其中is.expr是等于0的一个变量

###################################### Quality control ###############################################
## 非常常见的QC步骤是根据线粒体gene的转录本在所有转录本中所占的比例；这是基于一个生物学原理，即如果细胞破坏裂解了，处在
## 细胞质中的RNA会损失，然后线粒体中的RNA因为包绕在线粒体的细胞器膜中，并没有收到损伤，所以得到相对完整的保留。
mito.genes <- grep(pattern = "^mt-", x = rownames(x = pancreas_1@data), value = TRUE) # 挑取mt-开头的gene_symbol，是线粒体基因，可以用于归一化
mito.genes; length(mito.genes) # 一共有13个基因是线粒体基因
## 当然在这个数据集中，mito.genes为空
percent.mito <- Matrix::colSums(pancreas_1@raw.data[mito.genes, ]) / Matrix::colSums(pancreas_1@raw.data) # 计算每一个样本中线粒体基因在所有counts中的比例
head(percent.mito)
summary(percent.mito) # 我们看一下线粒体基因所占的比例的分布

# check out the meta data
head(pancreas_1@meta.data) # 本质是一个矩阵
# AddMetaData: adds columns to object@meta.data, and is a great place to stash QC stats。可以将计算得到的percent.mito追加到Seurat对象的metadata里面
str(pancreas_1)
pancreas_1 <- AddMetaData(object = pancreas_1, metadata = percent.mito, col.name = "percent.mito")
str(pancreas_1) # 比较前后的差别，发现pancreas_1对象下的meta.data增加了percent.mito这个元素，其中nGene, nUMI和percent.mito都是向量
class(pancreas_1@meta.data$nGene); class(pancreas_1@meta.data$nUMI);class(pancreas_1@meta.data$percent.mito)
head(pancreas_1@meta.data)
VlnPlot(object = pancreas_1, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
ggplot(pancreas_1@meta.data, aes(nUMI, percent.mito)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = 0.09, linetype = "dashed", colour = "red")
## 我们查看一下平均的nGene, nUMI, percent.mito和中位数nGene, nUMI, percent.mito
apply(pancreas_1@meta.data[,c('nGene','nUMI','percent.mito')], 2, function(x) round(mean(x),4))  # 统计各参数的平均值
apply(pancreas_1@meta.data[,c('nGene','nUMI','percent.mito')], 2, function(x) round(median(x),4))   # 统计各参数的中位数

# GenePlot is typically used to visualize gene-gene relationships, but can be used for anything calculated by the Seurat object,
# i.e. columns in object@meta.data, PC scores etc.  For the PBMC dataset (not this one),since there is a rare subset of cells with an outlier level of high 
# mitochondrial percentage and also low UMI content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = pancreas_1, gene1 = "nUMI", gene2 = "percent.mito") # GenePlot是专门用来画二维散点图的
GenePlot(object = pancreas_1, gene1 = "nUMI", gene2 = "nGene")
## Left: There are some clear outliers in mitochondrial RNA vs. the poly-A selected RNA. Right: The more unique molecules 
## captured, the more genes that are probed.

GenePlot(object = pancreas_1, gene1 = "her15.1", gene2 = "epcam", 
         cex.use = 1, col.use = 'red', pch.use = '.')   # 在这边也可以check两个gene之间表达的相关性
GenePlot(object = pancreas_1, gene1 = "ins", gene2 = "sst2", 
         cex.use = 1, col.use = 'blue', pch.use = '.') # 我们使用rownames(pancreas_1@raw.data)来挑选gene


# 下一步是对细胞的筛选，根据每个细胞检测到的gene数和每个细胞的线粒体转录本的比例
# We filter out cells that have unique gene counts over 2,500 or less than 200 Note that low.thresholds and high.thresholds are 
# used to define a 'gate'.  -Inf and Inf should be used if you don't want a lower or upper
# threshold. "可以看到这里选择的QC标准是 200~2500基因范围内，以及线粒体基因表达占比小于5%（这个值是人为设定的，也可以是9%等等）的才保留。"

# manual check; I already know all cells have >200 genes
table(pancreas_1@meta.data$percent.mito < 0.09 & pancreas_1@meta.data$nGene<2500)

pancreas_1 <- FilterCells(object = pancreas_1, subset.names = c("nGene", "percent.mito"),    # 使用FilterCells()函数来过滤细胞
                          low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.09))  
pancreas_1   # 我们看到手工check的结果和使用FilterCells()是一样的。

## 针对这个数据集，我使用nGene和线粒体的比例作为过滤条件，过滤低于200个gene和高于2500个gene的细胞，同时过滤线粒体比例超过9%的细胞

###################################### Normalization ###############################################
# 数据的normalization的目的是让细胞之间可以比较，Seurat当中内嵌的数据归一化方法是log normalization
# 就是将基因除以这个细胞的total expression，然后scale by 10000，再log转换，很类似于CPM后在log转换
# 这种归一化方法是到目前位置Seurat唯一内嵌的方法，但是document里面也明确指出不久后更多的方法会被囊括进来。

# "这里默认根据细胞测序文库大小进行normalization，简单的做一个log转换即可。" 我们来解读下这句话的含义：
# 出自hemberg的注解：“After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ 
# a global-scaling normalization method LogNormalize that normalizes the gene expression measurements for each cell by the total 
# expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result”
# 其意思是，先针对测序文库的大小进行normalization，然后乘以一个scaling factor(非常类似于CPM，不过默认值是10000), 然后在对这个数值进行log转换
# 注意，这种类似于CPM的数据校正方法，并没有对基因的长度进行校正；为什么我们的scaling factor是10000？这是与单细胞数据的测序量有关系的，测到的
# UMI就是在10000这个数量级

## Normalization前
pancreas_1@raw.data[,1]
summary(pancreas_1@raw.data[,1])
hist(colSums(pancreas_1@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")

## Normalization
pancreas_1 <- NormalizeData(object = pancreas_1, normalization.method = "LogNormalize",  # 执行完这个函数，会在pancreas_1这个对象中增加归一化后的信息
                            scale.factor = 10000)

## Normalization后
str(pancreas_1)  # 增加了NormalizedData
summary(pancreas_1@data[,1]) # 比较一下normalize前后的数据分布，normalization后变成了正态分布
hist(colSums(pancreas_1@data),
     breaks = 100,
     main = "Total expression before normalisation",
     xlab = "Sum of expression")

############################# Detection of variable genes across the single cells ###########################
## 寻找HVG，我们认为表达量在一定程度以上，同时样本(细胞)之间有一定差异的基因，是可以用于后续分析(cluster, cell type identification)的
## 其余的基因，会成为噪声，需要滤除
length(x = pancreas_1@var.genes)  # 在执行FindVariableGenes()之前，var.gene这个slot还是空着的

## FindVariableGenes()计算每个gene在不同细胞的平均表达量和散度(dispersion，z-score)
## The FindVariableGenes() function calculates the average expression and dispersion for each gene, places these genes into 
## bins, and then calculates a z-score for dispersion within each bin. I interpret that as take each gene, get the average 
## expression and variance of the gene across the 2,638 cells, categorise genes into bins (default is 20) based on their 
## expression and variance, and finally normalise the variance in each bin. 
## 以下参数是针对UMI data，对应不同的数据集，需要使用不同的参数。
pancreas_1@var.genes  # the variable genes slot is empty before the analysis
par(mfrow = c(1,1))
pancreas_1 <- FindVariableGenes(object = pancreas_1, mean.function = ExpMean,   # 这一步稍微有点费时，尤其是出图加text
                                dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.8)   # y.cutoff = 0.8属于条件比较宽的
## 给出的plot中，x轴是pancreas_1@hvg.info$gene.mean，y轴是pancreas_1@hvg.info$gene.dispersion.scaled.
## 同样，执行完上述函数，会在对象中增加var.genes这个对象
length(x = pancreas_1@var.genes)
head(pancreas_1@var.genes)   # vector of variable genes

# mean and variance of genes are stored pancreas_1@hvg.info
head(pancreas_1@hvg.info)

############################# Scaling the data and removing unwanted sources of variation #########################
# "需要去除那些technical noise,batch effects, or even biological sources of variation (cell cycle stage)"
# To mitigate the effect of confounding factors, Seurat constructs linear models to predict gene expression based on user-defined 
# variables. The scaled z-scored residuals of these models are stored in the scale.data slot, and are used for dimensionality 
# reduction and clustering. Seurat can regress out cell-cell variation in gene expression driven by batch, cell alignment rate 
# (as provided by Drop-seq tools for Drop-seq data), the number of detected molecules, mitochondrial gene expression and cell cycle. 
# Here we regress on the number of detected molecules per cell/percentage of mitochondria genes.
## Seurat使用linear model，基于用户定义的变量来预测基因表达量，从而去除混杂因素，这些混杂因素包括batch effect和cell cycle stage。
## 校正混杂因素(可以认为设定，比如batch)，同时对数据进行scale(数值减去平均值，然后除以对应的标准差，得到z-score)
## 注意，这一步也是为了后续PCA和tSNE做的准备工作，因为PCA过程中必须要有centering和scaling两步，减少outlier的影响
str(pancreas_1)
pancreas_1@scale.data   # slot is empty before running ScaleData()
pancreas_1 <- ScaleData(object = pancreas_1, vars.to.regress = c("nUMI", "percent.mito")) # 这一步比较费时(好几分钟，第一步regression，第二步scale data matrix)
str(pancreas_1) # scale.data从Null变成了校正后的结果，存放了实际表达量与线性模型的值的差别，这可以用来进行数据降维和聚类
summary(pancreas_1@scale.data[,1])
class(pancreas_1@scale.data)
pancreas_1@scale.data[1:6,1:6]

############################## 线性降维Linear dimensionality reduction ##############################
## 现在降维和可视化一般都是两步法（其实是三步法），第一步是上面select HVG，第二步PCA，将高维数据降成十几维，滤除差异贡献度很小的维度
## 以进一步去噪声，然后使用tSNE; 在进行PCA分析的时候，一般都是挑选HVG，同时需要进行centering和scaling两步的，比较好的是上面的ScaleData
## 已经帮我们计算得到了这部分var.genes。具体原理可以见StatQuest对PCA的更新版视频（2018-04）
## 即使我们挑选了一批gene出来，但是相当一部分这些gene在single cell的层面是not interesting的，因为它们单个可解释的变异度很小。
## 使用PCA的一个好处就是综合很多gene，得到解释变异度最大的PC。RunPCA默认针对@var.gene进行PCA分析，当然也可以私人订制。
pancreas_1 <- RunPCA(object = pancreas_1, pc.genes = pancreas_1@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5) # 增加了RunPCA这个对象
PrintPCAParams(pancreas_1)
PrintPCA(object = pancreas_1, pcs.print = 1:2, genes.print = 5, use.full = FALSE)  # 显示对各个PC贡献度最大的正负向gene。
### 每个PC实际给出了10个gene，应该是PC正负向贡献度最大的各5个gene(判断)

## 对PCA分析结果可以进行一系列的可视化： PrintPCA, VizPCA, PCAPlot, and PCHeatmap
par(mar = c(5,5,3,2))
VizPCA(object = pancreas_1, pcs.use = 1:2) # 通过VizPCA，也可以确认上述的判断；visualise top genes associated with principal components
PCAPlot(object = pancreas_1, dim.1 = 1, dim.2 = 2) # cells are coloured by their identity class according to pancreas_1@ident.
pancreas_1@ident

# PCA是基于most variable genes，ProjectPCA对每个gene进行打分，这是基于gene与calculated component的correlation。
# 这一步的好处是，能够找到那些非variable gene，但是与细胞的异质性关系密切的gene。
# ProjectPCA scores each gene in the dataset (including genes not included in the PCA) based on their correlation 
# with the calculated components. Though we don't use this further here, it can be used to identify markers that 
# are strongly correlated with cellular heterogeneity, but may not have passed through variable gene selection. 
# The results of the projected PCA can be explored by setting use.full=T in the functions above
pancreas_1 <- ProjectPCA(object = pancreas_1, do.print = T)

## 最重要的就是 PCHeatmap 函数了，默认是绘制第一个PC，显示在指定的cells.use的细胞中的30个gene。如果我们人为设定cells.use
## 为一个确定的数值，函数会自动寻找最extreme的细胞用于画图，在大型数据集中会大大加快绘图速度。
## 绘制第一张图会费一点功夫，而且会有warning message，直接忽略。
PCHeatmap(object = pancreas_1, pc.use = 1, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = pancreas_1, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = pancreas_1, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = F)
### 需要指出的是，使用PCHeatmap也可以帮助决定使用多少个PC，这被称作是supervised analysis。该函数可以展示细胞和基因的极端值，并且
### 能够帮助排除哪些被ribosomal/mitochondrial gene抑或是cell cycle gene驱动的PC


## Seurat出了提供RunPCA以及后面提到RunTSNE以外，还有RunICA和RunDiffusionMap这些常用的scRNA-seq数据降维方法。


## 数据降维的每一步操作都会存储到object@dr这个slot里面，比如RunPCA之后就会出现object@dr$pca，其中的内容也非常丰富：
### cell.embeddings：每个细胞的坐标信息
### gene.loadings：每个参与PCA的基因，在各个PC上的贡献，即loading score
### gene.loadings.full：投射到所有基因（不仅仅是用于PCA的基因）在各个PC上的贡献，即loading score
### sdev：每个PC在总变异上的贡献值
### key：cell.embeddings和gene.loading的名字，比如PC
### jackstraw: jackstraw procedure的结果
### misc：任何其他的information

## 作者定义了一些函数来调取这些内容，比如：GetCellEmbeddings, GetGeneLoadings, GetDimReduction
head(x = GetCellEmbeddings(object = pancreas_1, reduction.type = "pca", dims.use = 1:5))
head(x = GetGeneLoadings(object = pancreas_1, reduction.type = "pca", dims.use = 1:5))
# We also provide shortcut functions for common dimensional reduction
# techniques like PCA PCAEmbed and PCALoad() will pull the PCA cell
# embeddings and gene loadings respectively
head(x = GetDimReduction(object = pancreas_1, reduction.type = "pca", slot = "sdev"))


############################################## MDS降维，无需求别run ##################################
# 我们还可以计算一些Seurat当中没有定义的数据降维方法，比如MDS, 
####
# Before running MDS, we first calculate a distance matrix between all pairs of cells.  Here we use a simple euclidean distance 
# metric on all genes, using object@scale.data as input
d <- dist(x = t(x = pancreas_1@scale.data))
# Run the MDS procedure, k determines the number of dimensions
mds <- cmdscale(d = d, k = 2)  # 很费时
# cmdscale returns the cell embeddings, we first label the columns to ensure downstream consistency
colnames(x = mds) <- paste0("MDS", 1:2)
# We will now store this as a new dimensional reduction called 'mds'
pancreas_1 <- SetDimReduction(object = pancreas_1, reduction.type = "mds", slot = "cell.embeddings", 
                        new.data = mds)
pancreas_1 <- SetDimReduction(object = pancreas_1, reduction.type = "mds", slot = "key", 
                        new.data = "MDS")

# We can now use this as you would any other dimensional reduction in all downstream functions (similar to PCAPlot, but 
# generalized for any reduction)
DimPlot(object = pancreas_1, reduction.use = "mds", pt.size = 0.5)


# If you wold like to observe genes that are strongly correlated with the first MDS coordinate (similar to ProjectPCA, but 
# generalized for any reduction):
pbmc <- ProjectDim(object = pancreas_1, reduction.type = "mds")

# Display the results as a heatmap (similar to PCHeatmap, but generalized for any dimensional reduction)
DimHeatmap(object = pancreas_1, reduction.type = "mds", dim.use = 1, cells.use = 500, 
           use.full = TRUE, do.balanced = TRUE, label.columns = FALSE, remove.key = TRUE)

# Explore how the first MDS dimension is distributed across clusters
VlnPlot(object = pancreas_1, features.plot = "MDS1", x.lab.rot = TRUE)

# See how the first MDS dimension is correlated with the first PC dimension
GenePlot(object = pancreas_1, gene1 = "MDS1", gene2 = "PC1")
################################### 找到有统计学显著性的主成分 #####################################
# 主成分分析结束后需要确定哪些主成分所代表的基因可以进入下游分析，这里可以使用JackStraw做重抽样分析(默认每次重抽样1%的数据)。
# Seurat randomly permutes a subset of the data (1% by default) and reruns PCA, constructing a null distribution of gene scores by 
# repeating this procedure. We identify significant PCs as those who have a strong enrichment of low p-value genes:
# JackStraw返回的结果是每个gene在每个PC当中的p-value
# 可以用JackStrawPlot可视化看看哪些主成分可以进行下游分析。这一步很耗时
system.time(
pancreas_1 <- JackStraw(object = pancreas_1, num.replicate = 100) 
)
JackStrawPlot(object = pancreas_1, PCs = 1:12)
## 图形化展示12个主成分，有意义的主成分是那些能够富集较多low p-value gene的主成分，solid curve above the dashed line
## The JackStrawPlot function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform 
## distribution (dashed line). Significant PCs will show a strong enrichment of genes with low p-values (solid curve above the dashed 
## line). In this case it appears that PCs 1-8 are significant.


# 当然，也可以用最经典的碎石图来确定主成分。这一步很耗时。A more ad hoc method for determining which PCs to use is to look at 
# a plot of the standard deviations of the principle components and draw your cutoff where there is a clear elbow in the graph. 
# This can be done with PCElbowPlot.
PCElbowPlot(object = pancreas_1)  # 对于endo的数据集，我们挑选15个PC


# 这个确定主成分是非常有挑战性的: - The first is more supervised, exploring PCs to determine relevant sources of heterogeneity, 
# and could be used in conjunction with GSEA for example. - The second implements a statistical test based on a random null model, 
# but is time-consuming for large datasets, and may not return a clear PC cutoff. The third is a heuristic that 
# is commonly used, and can be calculated instantly.

############################################### Cluster the cells #################################################
# Seurat的clustering是基于graph-based clustering approach inspired by SNN-Cliq和PhenoGraph；这个算法相当的technical
# save.SNN = T saves the SNN so that the clustering algorithm can be rerun using the same graph
# but with a different resolution value (see docs for full details)
str(pancreas_1)
## 在这个算法中，resolution参数能够校正granularity of the clustering, 其数值越高能够分出的cluster就越多；在3000个细胞左右的数据集中
### resolution设定在0.6-1.2左右往往能给出不错的结果。cluster信息是保存在@ident这个slot里面，还记得PCA颜色么！！
pancreas_1 <- FindClusters(object = pancreas_1, reduction.type = "pca", dims.use = 1:15, resolution = 0.6, print.output = 0, save.SNN = TRUE)
str(pancreas_1) # 增加了FindClusters这个元素
pancreas_1@ident # The clusters are saved in the object@ident slot(向量).具体来看每一个细胞被分到哪一个cluster
table(pancreas_1@ident) # 统计每个cluster有多少个细胞

# A useful feature in Seurat v2.0 is the ability to recall the parameters that were used in the latest function calls 
# for commonly used functions. For FindClusters, we provide the function PrintFindClustersParams to print a nicely 
# formatted formatted summary of the parameters that were chosen.
PrintFindClustersParams(object = pancreas_1)
## While we do provide function-specific printing functions, the more general function to 
## print calculation parameters is PrintCalcParams(). 我们以RunPCA为例子
PrintCalcParams(object = pancreas_1, calculation = 'RunPCA')


####################################### Run Non-linear dimensional reduction (tSNE) ##################################
# 同样也是一个函数，这个结果也可以像PCA分析一下挑选合适的PC进行下游分析。
# 这一步很耗时，可以保存该对象，便于重复，以及分享交流 （出图时间很长)
# 我们在进行tSNE计算的时候，倾向使用pc，这样噪声比较小；当然，在函数中也可以使用genes.use这个参数来使用z-scaled的gene表达来计算
# 而且作者建议使用上面clustering相同的PC
pancreas_1 <- RunTSNE(object = pancreas_1, dims.use = 1:15, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters 
TSNEPlot(object = pancreas_1, do.label=T, do.identify = T)
save(pancreas_1, file = "pancreas_1_endo.rData")



############################### Finding differentially expressed genes (cluster biomarkers) ###########################
# 差异分析在seurat包里面被封装成了函数：FindMarkers(); Seurat是基于非参数检验当中的Wilcoxon rank sum test进行差异基因表达
## 分析的，其表示方法为"wilcox"。之前的默认方法是“bimod”(Likelihood-ratio test for single cell gene expression)。
## “roc” : Standard AUC classifier
## “t” : Student’s t-test
## “tobit” : Tobit-test for differential gene expression
## “poisson” : Likelihood ratio test assuming an underlying negative binomial distribution. Use only for UMI-based datasets
## “negbinom” : Likelihood ratio test assuming an underlying negative binomial distribution. Use only for UMI-based datasets
## “MAST” : GLM-framework that treates cellular detection rate as a covariate
## “DESeq2” : DE based on a model using the negative binomial distribution 
## 在函数中使用test.use参数来指定进行差异基因表达的方法

# 默认设置是给出positive和negative的markers，one cluster compared against genes in all other cells
# p_val：未校正的p值
# avg_logFC：正值表示gene在第一组中的表达量高于第二组中
# pct.1：gene在第一组中被检测到的比例
# pct.2: gene在第二组中被检测到的比例
# p_val_adj：基于bonferroni校正的多重假设检验的校正后的p值

# 如果ident.2这个参数没有明确指定是哪个组的话，FindMarkers函数会检测ident.1与其他所有细胞的差异表达gene

# 为了加快计算速度，Seurat允许我们在进行差异基因分析的时候，先对gene或者细胞进行过滤，可以用于过滤的参数有：
## min.pct：minimal percentage, gene在细胞中的最小检出率
## logfc.threshold：A组相比于B组logFC比值在某个值以上的
## min.diff.pct: gene的检出率在两个cluster的细胞中，相差不超过的值，如25%，则用0.25表示
## max.cells.per.ident：subsample each group to a maximum of 多少个细胞，这是针对大的cluster或者computationally intensive
##                      DE
## 注意，我们限定了min.pct, logfc.threshold和min.diff.pct都会加快计算速度，但是会丢失一部分被filter掉的gene


# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = pancreas_1, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))

# find markers of each cluster
for (i in 0:13) {
  cluster.markers <- FindMarkers(object = pancreas_1, ident.1 = i, min.pct = 0.25)
  print(paste("marker of cluster: ",i ))
  print(x = head(x = cluster.markers, n = 20))
}

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = pancreas_1, ident.1 = 5, ident.2 = c(0,3), min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))


# find markers for every cluster compared to all remaining cells, report only the positive ones (this step takes some time!)
pancreas_1.markers <- FindAllMarkers(object = pancreas_1, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
head(pancreas_1.markers)
pancreas_1.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
# A tibble: 16 x 6
# Groups:   cluster [8]
pancreas_1.markers $cluster

# FindMarkers，针对某一个cluster，有一系列参数可以选择，然后又4种找差异基因的算法，使用test.use参数进行设置
# 方法一： ROC test (“roc”)
# 方法二： t-test (“t”)
# 方法三： LRT test based on zero-inflated data (“bimod”, default)
# 方法四： LRT test based on tobit-censoring models (“tobit”)，这个比较费时
# 值得注意的是： The ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect).
# 所以可以用来衡量找到的marker是否可靠，我们来比较一下使用以上四种方法找到的cluster1的marker
cluster0.markers_bimod <- FindMarkers(object = pancreas_1, ident.1 = 0, thresh.use = 0.25, test.use = "bimod", only.pos = TRUE)
cluster13.markers_roc <- FindMarkers(object = pancreas_1, ident.1 = 13, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)
cluster0.markers_t <- FindMarkers(object = pancreas_1, ident.1 = 0, thresh.use = 0.25, test.use = "t", only.pos = TRUE)
cluster0.markers_tobit <- FindMarkers(object = pancreas_1, ident.1 = 0, thresh.use = 0.25, test.use = "tobit", only.pos = TRUE)

# identical set of genes
dim(cluster1.markers_bimod); head(cluster1.markers_bimod)
dim(cluster1.markers_roc); head(cluster1.markers_bimod)
dim(cluster1.markers_t); head(cluster1.markers_bimod)
dim(cluster1.markers_tobit); head(cluster1.markers_bimod)

# the rankings of the genes are quite similar between the methods
my_gene <- row.names(cluster1.markers_bimod)
a <- 1:length(my_gene)
b <- match(my_gene, row.names(cluster1.markers_roc))
c <- match(my_gene, row.names(cluster1.markers_t))
d <- match(my_gene, row.names(cluster1.markers_tobit))

# bimod vs. bimod
cor(a, a, method = "spearman")
# bimod vs. roc
cor(a, b, method = "spearman")
# bimod vs. t
cor(a, c, method = "spearman")
# bimod vs. tobit
cor(a, d, method = "spearman")

par(mfrow=c(2,2))
barplot(a, main = 'bimod')
barplot(b, main = 'roc')
barplot(c, main = 't')
barplot(d, main = 'tobit')

# 同时，该包提供了一系列可视化方法来检查差异分析的结果的可靠性：
# VlnPlot (shows expression probability distributions across clusters)
# FeaturePlot (visualizes gene expression on a tSNE or PCA plot) are our most commonly used visualizations
# RidgePlot, CellPlot, and DotPlot 来试一下
## RidgePlot 是来自于ggridges这个包，在每个cluster看单个gene的表达丰度的分布distribution
## Dotplot中点的大小是这个cluster中表达这个gene的比例，颜色代表的是表达量

VlnPlot(object = pancreas_1, features.plot = c("ins", "gcga","gcgb", 'sst2','try','her15.1'))
RidgePlot(object = pancreas_1, features.plot = c("ins", "gcga",'sst2'))
CellPlot(pancreas_1,pancreas_1@cell.names[1],pancreas_1@cell.names[2],do.ident = FALSE)
DotPlot(pancreas_1,genes.plot = rownames(cluster1.markers_bimod)[1:10], plot.legend = T)

# you can plot raw UMI counts as well，使用use.raw = T来指定
VlnPlot(object = pancreas_1, features.plot = c('nf2a','nf2b','mitfb','cdc42'), use.raw = TRUE, y.log = TRUE)

# In tSNE plot
FeaturePlot(object = pancreas_1, 
            features.plot = c('nf2a'), 
            cols.use = c("grey", "Red"), reduction.use = "tsne",
            no.legend = F)  # 加上了表征表达量的图例 

FeaturePlot(object = pancreas_1, 
            features.plot = c('nf2a'), 
            cols.use = c("grey", "Red"), reduction.use = "tsne",
            no.legend = F,
            min.cutoff = 0.2,    # 为了增加对比度，可以设定min和max cutoff
            max.cutoff = 1)  


# Calculate gene-specific contrast levels based on quantiles of non-zero expression. Particularly useful when plotting 
# multiple markers 图例的颜色是根据不同的gene的上下10%来界定的。
FeaturePlot(object = pancreas_1, features.plot = c("ins", "gcga"), no.legend = FALSE, 
            min.cutoff = "q10", max.cutoff = "q90", dark.theme = T)

# 在一张tSNE上同时展现两个gene，这个非常有用，可以展示比如bihormonal gene
FeaturePlot(object = pancreas_1, features.plot = c("ins", "gcga"), 
            cols.use = c("grey", "red", "blue", "green"), 
            overlay = TRUE, no.legend = FALSE)


FeaturePlot(object = pancreas_1,          # 显示某个cluster的marker gene，很不错的方法
            features.plot = head(row.names(cluster7.markers_roc),9), 
            cols.use = c("grey", "Red"), reduction.use = "tsne")


FeaturePlot(object = pancreas_1, 
            features.plot = c('cdh5'), 
            cols.use = c("grey", "Red"), reduction.use = "tsne",
            do.hover = T, data.hover = c("ident", "PC1", "nGene",'percent.mito', 'ClusterNames_0.6', 'res.0.6', 'res.1'),
            do.identify = TRUE, dark.theme = T)
pancreas_1@meta.data

# DoHeatmap generates an expression heatmap for given cells and genes. In this case, we are plotting the top 20 markers 
# (or all markers if less than 20) for each cluster.
head(pancreas_1.markers); dim(pancreas_1.markers)
pancreas_1.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10 # 显示每个cluster的marker gene
head(top10) # A tibble: 6 x 6; Groups:cluster [1]
class(top10)

# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
DoHeatmap(object = pancreas_1, genes.use = top10$gene, order.by.ident = TRUE, slim.col.label = TRUE, remove.key = TRUE)

DoHeatmap(object = SubsetData(object = pancreas_1, max.cells.per.ident = 50), 
          genes.use = c('ins','gcga','sst2','ela2l','epcam'), 
          slim.col.label = TRUE, group.label.rot = TRUE)

###################################### Assigning cell type identity to clusters ####################################
# https://mp.weixin.qq.com/s/QZD1tvCgZVa5PQtbjvrkrg
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
new.cluster.ids <- c("Delta cells",
                     "Beta cells",
                     "Smooth Muscle cells",
                     "Alpha cells",
                     "Acinar cells",
                     "Endothelial cells",
                     "Centroacinar cells",
                     "Cell with both endocrine and exocrine markers",
                     "Unknown cells - exocrine probably",
                     "Unknown cells - endocrine probably",
                     "Unknown cells - exocrine probably",
                     "Blood cells - type 1",
                     "Ductal cells - type 2",
                     "Blood cells - type 2")

pancreas_1@ident <- plyr::mapvalues(x = pancreas_1@ident,
                              from = current.cluster.ids,
                              to = new.cluster.ids)   # 忽略warning message

TSNEPlot(object = pancreas_1, do.label = TRUE, pt.size = 0.5)

###################################### Further subdivisions within cell types ######################################
# https://mp.weixin.qq.com/s/QZD1tvCgZVa5PQtbjvrkrg
# Seurat provides the StashIdent() function for keeping cluster IDs; this is useful for testing various parameters and 
# comparing the clusters. For example, adjusting the parameters may lead to the CD4 T cells subdividing into two groups.
# stash cluster identities for later
# 我们在进行转换TSNEplot当中颜色的标记时，是需要控制object@indent的，所以第一步我们是先保存目前的identity到meta.data下的新的column
pancreas_1 <- StashIdent(object = pancreas_1, save.name = "ClusterNames_0.6")

pancreas_1 <- FindClusters(object = pancreas_1,
                     reduction.type = "pca",
                     dims.use = 1:15,
                     resolution = 1,   # 修改了resolution，变成了17个cluster
                     print.output = FALSE)

# plot two tSNE plots side by side, and colour points based on different criteria
plot1 <- TSNEPlot(object = pancreas_1,
                  do.return = TRUE,
                  no.legend = TRUE,
                  do.label = TRUE)

plot2 <- TSNEPlot(object = pancreas_1,
                  do.return = TRUE,
                  group.by = "ClusterNames_0.6",
                  no.legend = TRUE,
                  do.label = TRUE)

plot_grid(plot1, plot2)

#================================================================================================================
#
#          "Monocle2" package (Monocle 2.8.0)  生信技能树，重点Monocle对象的构建(RPKM或者counts data)
#                                      后续分析基于内置的RPKM数据进行
#
#================================================================================================================
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

# 这个数据集中包含了两种细胞类型：
# In the myoblast experiment, the culture contains fibroblasts that came from the original muscle biopsy used 
# to establish the primary cell culture. Myoblasts express some key genes that fibroblasts don't. Selecting 
# only the genes that express, for example, sufficiently high levels of MYF5 excludes the fibroblasts. Likewise, 
# fibroblasts express high levels of ANPEP (CD13), while myoblasts tend to express few if any transcripts of 
# this gene.


###################################### 构建S4对象，CellDataSet ######################################
# 读取矩阵的推荐层级：raw counts(UMI) > relative counts (FPKM/TPM) > raw counts(without UMI)
# 针对上述推荐层级，raw counts(UMI)在导入CellDataSet对象之前千万不要normalization或者把它变成FPKM/TPM的data(因为UMI raw counts是最优输入)。
# 主要是读取表达矩阵和样本描述信息，这里介绍两种方式，一种是读取基于 subjunc+featureCounts 分析后的reads counts矩阵，
# 一种是读取 tophat+cufflinks 得到的RPKM表达矩阵
# 如何读取外部数据，请参考“https://mp.weixin.qq.com/s/zCfDkxbVTxjFQ5QAIULYjA”
# 读取上游分析的输出文件(比如表达矩阵和样本注释都是txt格式的文件)
setwd('/Users/mijiarui/Nature_Biotechnology_Paper/Testing_dataset')
library(monocle)
library(scater, quietly = TRUE)
library(knitr)
options(stringsAsFactors = FALSE)

# 这个文件是表达矩阵，包括线粒体基因和 ERCC spike-ins 的表达量，可以用来做质控
molecules <- read.table("molecules.txt", sep = "\t")
dim(molecules)

## 这个文件是表达矩阵涉及到的所有样本的描述信息，包括样本来源于哪个细胞，以及哪个批次。
anno <- read.table("annotation.txt", sep = "\t", header = TRUE)  # anno是phenoData，存储的是样品的信息
dim(anno)
rownames(anno)=colnames(molecules)  # 为了让phenoNames is the same between assayData and phenoData，必不可少(重要)

## 以下是为了在featureData当中增加gene_short_name一列
library(org.Hs.eg.db)  
eg2symbol=toTable(org.Hs.egSYMBOL)
eg2ensembl=toTable(org.Hs.egENSEMBL)
egid=eg2ensembl[ match(rownames(molecules),eg2ensembl$ensembl_id),'gene_id']
symbol=eg2symbol[match( egid ,eg2symbol$gene_id),'symbol']
gene_annotation = data.frame(ensembl=rownames(molecules),   # gene_annotation是featureData，存储了基因的ensembl_id, gene_symbol和entrez_id
                             gene_short_name=symbol,
                             egid=egid)
rownames(gene_annotation)=rownames(molecules)  # 为了让featureNames is the same between assayData and featureData，必不可少(重要)

# 真正开始构建CellDataSet对象了：
pd <- new("AnnotatedDataFrame", data = anno)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
# tung <- newCellDataSet(as.matrix(molecules), phenoData = pd, featureData = fd)
tung <- newCellDataSet(as(as.matrix(molecules), "sparseMatrix"),   
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit=0.5,
                       expressionFamily=negbinomial.size())

tung
## sparseMatrix对于超大的数据集甚至是中等数据集都可以可以节省不少空间并提高运行速度，所以：To work with your data in a sparse format, 
## simply provide it to Monocle as a sparse matrix from the Matrix package；也就是使用Matrix函数将普通的输入密集矩阵转变为稀疏矩阵。
## 很多pipeline，例如CellRanger下机的数据本身就是稀疏矩阵了(数据格式为MTX)，在构建CellDataSet对象的时候，直接导入即可，别先使用as.matrix()
## 先转成密集矩阵，在用Matrix包再转成稀疏矩阵。针对CellRanger的输出，请进一步参考Monocle官网教程：
## http://cole-trapnell-lab.github.io/monocle-release/docs/#recommended-analysis-protocol

## 针对expressionFamily一共有4种参数可供选择，negbinomial.size()，negbinomial()，tobit()，gaussianff()。其中第一种和第三种比较常用。
## negbinomial()会比negbinomial.size()稍稍更加准确一点，但是运行速度慢了不少，所以不建议使用。tobit()是专门针对FPKM和RPMK的。gaussianff()
## 是针对已经normalization后正态分布的数据。正如之前所说Monocle是不建议先normalize在构建CellDataSet对象的，而且使用这个方法，后续部分
## Monocle feature不好使。另外值得注意的是，对于相对表达量数据，如RPKM/TPM，我们在进行处理的时候，会将其转换成transcript counts(通过
## 使用函数relative2abs())，转换产生的绝对counts就服从负二项分布了，可以使用negbinomial.size()，并且比相对值使用tobit()来得更好。


# 在这里我们读取HSMMSingleCell包中的测试数据，或者使用内置数据个构建S4对象：
getwd()
pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)

# First create a CellDataSet from the relative expression levels

## 这里仅仅是针对rpkm表达矩阵的读取，我们选取的统计学模型是tobit，这种方法没有把rpkm转换成counts的方法好
## Tobits are truncated normal distributions. Using tobit() will tell Monocle to log-transform your data where 
## appropriate. Do not transform it yourself.
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),   
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit=0.1,
                       expressionFamily=tobit(Lower=0.1))  # 注意因为读入的是RPKM，所以统计学模型使用的是tobit

# Next, use it to estimate RNA counts. RPC的含义是mRNAs per cell; rpkm格式的表达值需要转换成reads counts之后才可以进行下游分析
# Monocle 2 includes an algorithm called Census which performs this conversion. Census算法就是用来解决这个转换问题的。
rpc_matrix <- relative2abs(HSMM)    # 从校正后的相对值转换为绝对值
rpc_matrix[1:10,1:5] 

# Now, make a new CellDataSet using the RNA counts，既然转变成了counts data，则服从了负二项分布的规律
# 统计学模型相应的进行改变：Negative binomial distribution with fixed variance (which is automatically 
# calculated by Monocle). Recommended for most users. 同时lowerDetectionLimit也要相应更改。
HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit=0.5,  #  we have changed the value of lowerDetectionLimit to reflect the new scale of expression
                       expressionFamily=negbinomial.size())   # 注意，将RPKM或者TPM改成了counts以后，统计学模型也就相应改成了负二项分布

## 下面的分析，都基于内置数据构建的S4对象，HSMM

# 必要的矫正（量化因子和散度）：Size factors help us normalize for differences in mRNA recovered across cells, 
# and "dispersion" values will help us perform differential expression analysis later.
# 这两个estimate都是针对负二项分布的数据才有用。
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
### Warning: Deprecated, use tibble::rownames_to_column() instead.
### 至此，我们的HSMM构建完成，可以进行后续的分析了。

###################################### 过滤低质量细胞和未检测到的基因 ######################################
# 因为有可能有空的well，dead cell或者doublet的存在，会影响后续的分析，比如pseudotime。
# 基于基因的过滤
##  这里只是把基因挑选出来，并没有对S4对象进行过滤操作。 这个detectGenes函数还计算了每个细胞里面表达的基因数量。有点类似于CalculateQCMetrics
HSMM <- detectGenes(HSMM, min_expr = 0.1)  # a gene is “expressed” if there is at least one count since we set min_expr = 0.1
print(head(fData(HSMM)))    # num_cells_expressed，某个gene在多少个细胞当中有表达

## 对每个基因都检查一下在多少个细胞里面是有表达量的。
## 只留下至少在10个细胞里面有表达量的那些基因，做后续分析，当然如果你感兴趣的是一些稀有细胞，这个参数的设定就要琢磨一下了
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
length(expressed_genes) ## 只剩下了14224个基因

## The HSMM dataset included with this package has scoring columns built in. 说白了就是写质控数据反应这个细胞的测序
## 结果能不能在后续进行使用，对细胞(样本)进行过滤。
print(head(pData(HSMM))) 

valid_cells <- row.names(subset(pData(HSMM),     # 根据pData当中的参数，对细胞进行过滤
                                Cells.in.Well == 1 &
                                  Control == FALSE &
                                  Clump == FALSE &
                                  Debris == FALSE &
                                  Mapped.Fragments > 1000000))
HSMM <- HSMM[,valid_cells]

# 基于样本表达量进行过滤，以下代码说白了就是眼见为实，看看不同细胞的mRNA的分布情况
## 这里选择的是通过不同时间点取样的细胞来进行分组查看，把超过2个sd 的那些样本的临界值挑选出来，下一步过滤的时候使用。
HSMM; str(HSMM)
pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) +    # 用来在后面的图上画线，设定比较高的阈值是防止出现双细胞或者多细胞的情况出现
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) -    # 用来在后面的图上画线
                     2*sd(log10(pData(HSMM)$Total_mRNAs)))
table(pData(HSMM)$Hours)
qplot(Total_mRNAs, data = pData(HSMM), color = Hours, geom = "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

# 执行过滤并可视化检查一下；我们将不满足要求的细胞过滤掉以后，对每个细胞的表达量进行log转换后，理论上应该服从正态分布。
## 上面已经根据基因表达情况以及样本的总测序数据选择好了阈值，下面就可以可视化并且对比检验一下执行过滤与否的区别。
HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound & 
               pData(HSMM)$Total_mRNAs < upper_bound]                                 
HSMM <- detectGenes(HSMM, min_expr = 0.1)  # a gene is “expressed” if there is at least one count since we set min_expr = 0.1
HSMM

# Log-transform each value in the expression matrix.
L <- log(exprs(HSMM[expressed_genes,]))
# Standardize each gene, so that they are all on the same scale, then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
# Plot the distribution of the standardized gene expression values.
qplot(value, geom="density", data=melted_dens_df) +  stat_function(fun = dnorm, size=0.5, color='red') + 
  xlab("Standardized log(FPKM)") +
  ylab("Density")


########################################### 聚类 ###########################################
# 根据指定基因对单细胞转录组表达矩阵进行分类(Classify cells with known marker genes)
# leverage your knowledge of key marker genes to quickly and easily classify your cells by type:
## 下面这个代码只适用于这个测试数据， 主要是生物学背景知识，用MYF5基因和ANPEP基因来对细胞进行分类，可以区分Myoblast和Fibroblast。
## 如果是自己的数据，建议多读读paper看看如何选取合适的基因，或者干脆跳过这个代码。
## 根据基因名字找到其在表达矩阵的ID，这里是ENSEMBL数据库的ID（MYF5_id和ANPEP_id）
## For example, you could provide a function for each of several cell types. These functions accept as input the expression data for 
## each cell, and return TRUE to tell Monocle that a cell meets the criteria defined by the function. So you could have one function 
## that returns TRUE for cells that express myoblast-specific genes, another function for fibroblast-specific genes, etc. Here's an 
## example of such a set of "gating" functions
MYF5_id <- row.names(subset(fData(HSMM), gene_short_name == "MYF5"))
ANPEP_id <- row.names(subset(fData(HSMM), gene_short_name == "ANPEP"))
## 这里选取的基因取决于自己的单细胞实验设计

# 根据特征基因的表达量来给细胞赋予cell type属性，深度理解函数newCellTypeHierarchy()
# Monocle provides a simple system for tagging cells based on the expression of marker genes of your choosing. 
# You simply provide a set of functions that Monocle can use to annotate each cell. For example, you could 
# provide a function for each of several cell types. These functions accept as input the expression data for 
# each cell, and return TRUE to tell Monocle that a cell meets the criteria defined by the function. So you 
# could have one function that returns TRUE for cells that express myoblast-specific genes, another function 
# for fibroblast-specific genes, etc. Here's an example of such a set of "gating" functions:
# The functions are organized into a small data structure called a CellTypeHierarchy, that Monocle uses to 
# classify the cells. You first initialize a new CellTypeHierarchy object, then register your gating functions 
# within it. Once the data structure is set up, you can use it to classify all the cells in the experiment:

## 首先使用newCellTypeHierarchy()进行初始化
cth <- newCellTypeHierarchy() 

## 然后添加对细胞进行分类的内容
cth <- addCellType(cth, "Myoblast", classify_func = function(x) { x[MYF5_id,] >= 1 })
cth <- addCellType(cth, "Fibroblast", classify_func = function(x){ x[MYF5_id,] < 1 & x[ANPEP_id,] > 1 })

## 使用classifyCells()函数对每个细胞进行有监督的分类
HSMM <- classifyCells(HSMM, cth, 0.1) ## 这个时候的HSMM已经被改变了，增加了属性(phenoData中的CellType)。
HSMM 
## The function classifyCells applies each gating function to each cell, classifies the cells according to the 
## gating functions, and returns the CellDataSet with a new column, CellType in its pData table. We can now 
## count how many cells of each type there are in the experiment.
pData(HSMM)$CellType
table(pData(HSMM)$CellType)

## 可视化细胞分类和组成：
pie <- ggplot(pData(HSMM), aes(x = factor(1), fill = factor(CellType))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

### 可以看到还有很大一部分细胞仅仅是根据这两个基因的表达量是无法成功的归类的。这个是很正常的，
### 因为单细胞转录组测序里面的mRNA捕获率不够好。 通过这个步骤成功的给HSMM这个S4对象增加了一个属性，就是CellType，
### 在下面的分析中会用得着。注意Unknown意味这一个条件都没有满足；ambiguous意味着满足多个分类指标
### Note that many cells are marked "Unknown". This is common, largely because of the low rate of mRNA capture 
### in most single-cell RNA-Seq experiments. A cell might express a few MYF5 mRNAs, but we weren't lucky enough 
### to capture one of them. When a cell doesn't meet any of the criteria specified in your classification 
### functions, it's marked "Unknown". If it meets multiple functions' criteria, it's marked "Ambiguous". You 
### could just exclude such cells, but you'd be throwing out a lot of your data. In this case, we'd lose more 
### than half of the cells! 针对这些Unknown的细胞，可以考虑使用无监督的分类方法，基于的函数是clusterCells。它基于
### 的方法是即使只检测到了MYF5，但是表达了许多其它的myoblast基因


## 无监督聚类: 这里需要安装最新版R包才可以使用里面的一些函数，因为上面的步骤基于指定基因的表达量进行细胞分组会漏掉很多信息，
## 使用无监督聚类，相当于要对一部分unknown cell进行impute，使用的函数是clusterCells()。它是基于整体的基因表达情况来帮助聚类的。
## 比如某一些细胞缺乏MYF5，但是它们表达相当一部分myoblast特征性基因，那么仍然可以将这些细胞归类到myoblast当中去。
## 当然设计到聚类，就一定有降维去噪声的步骤。
## 所以需要更好的聚类方式。在进行无监督的分类的时候，筛选高表达的HVG可以尽可能的提高信噪比
## 首先我们挑选表达量比较高的基因:
## The first step is determining a subset of genes to use for clustering; this is because not all genes are informative, such as those 
## that are lowly expressed. The approach is to select gene based on their average expression and variability across cells. 
## The dispersionTable() function calculates the mean and dispersion values.
disp_table <- dispersionTable(HSMM)   # dispersionTable()计算平均值和散度
head(disp_table)
## 只有满足条件的10198个基因才能进入聚类分析，我们使用setOrderingFilter函数来标记那部分后续用clusterCell聚类的基因；
## 使用plot_ordering_genes函数来图形化展示基因的表达量和离散度，高表达的基因都可以纳入考虑，红色线是在不同的平均基因
## 表达量下离散度的估计值。
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)  # We will select genes, which have a mean expression >= 0.1, to use in the clustering step.
table(disp_table$mean_expression>=0.1) # 显示True的那部分是用来clustering的gene
 
## setOrderingFilter() function allows us to indicate which genes we want to use for clustering.
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id) # setOrderingFilter这一步可以将用于cluster的gene给标记出来
plot_ordering_genes(HSMM)   # 绘制散点图(dispersion vs average expression)
## 这里看看基因的表达量和基因的变异度之间的关系
## 处在灰色阴影区域的基因会被抛弃掉，不进入聚类分析。

# 聚类分析之前，首先对数据进行PCA降维和去噪，当然PCA需要对log转换后的数据进行分析
# HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(HSMM, return_all = F) # norm_method = 'log', 这一步往往比较耗时。
HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 6,     # num_dim参数决定我们要纳入几个PC，如果不指定的话，默认值是50个PC
                        reduction_method = 'tSNE', verbose = T) 
HSMM <- clusterCells(HSMM, num_clusters = 2)
## 这里先用tSNE的聚类方法处理HSMM数据集，并可视化展示；这个函数目前有问题，尤其是加入markers这个参数
plot_cell_clusters(HSMM, 1, 2, color_by  = 'CellType', markers = c("MYF5", "ANPEP"))  # Monocle默认使用tSNE来进行cluster的可视化,但是这个函数有点问题
### 似乎去除markers = c("MYF5", "ANPEP")就正确了
## 可以看到并不能把细胞类型完全区分开，这个是完全有可能的，因为虽然是同一种细胞，但是有着不同的培养条件。
head(pData(HSMM)) # 在这里另外一个问题是pdata中的Cluster，只有一个，这是不符合常理的，需要debug一下，已经在github上提问。
head(fData(HSMM))
## 所以这里也区分一下 培养基， a high-mitogen growth medium (GM) to a low-mitogen differentiation medium (DM). 
plot_cell_clusters(HSMM, 1, 2, color="Media")

## Remove batch effect
## 因为我们假设就2种细胞类型，所以在做聚类的时候可以把这个参数添加进去(去除medium带来的混杂因素)，这样可以去除无关变量的干扰。
## Monocle allows us to subtract the effects of "uninteresting" sources of variation to reduce their 
## impact on the clustering. You can do this with the residualModelFormulaStr argument to clusterCells 
## and several other Monocle functions. This argument accepts an R model formula string specifying 
## the effects you want to subtract prior to clustering.
HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 2, reduction_method = 'tSNE', 
                        residualModelFormulaStr="~Media + num_genes_expressed", verbose = T) #
HSMM <- clusterCells(HSMM, num_clusters=2)
## Distance cutoff calculated to 1.284778
pData(HSMM)
plot_cell_clusters(HSMM, 1, 2, color="CellType", cell_size = 4) # cell_size只能对应一个确定的数值，alpha还是不能用

## Now that we've accounted for some unwanted sources of variation, we're ready to take another crack at 
## classifying the cells by unsupervised clustering:
HSMM <- clusterCells(HSMM, num_clusters = 2)
plot_cell_clusters(HSMM, 1, 2, color = 'Cluster') + facet_wrap(~CellType)   # pData中的Cluster存在问题，需要debug一下


## 半监督聚类，也就是利用一部分marker genes来进行分析。正如之前所说，如果我们只单单用一个marker来指示细胞的类型的
## 时候，某些细胞不能很好的区分，因此我们需要用那些与指明的marker gene有相同表达差异的gene(co-vary)，构建一个大的gene list
## 这里的差异分析非常耗时
marker_diff <- markerDiffTable(HSMM[expressed_genes,], 
                               cth, 
                               residualModelFormulaStr="~Media + num_genes_expressed",
                               cores=1)
head(marker_diff)
### The function markerDiffTable takes a CellDataSet and a CellTypeHierarchy and classifies all the cells into 
### types according to your provided functions. It then removes all the "Unknown" and "Ambiguous" functions 
### before identifying genes that are differentially expressed between the types. Often it's best to pick the 
### top 10 or 20 genes that are most specific for each cell type. This ensures that the clustering genes aren't 
### dominated by markers for one cell type. You generally want a balanced panel of markers for each type if 
### possible. Monocle provides a handy function for ranking genes by how restricted their expression is for 
### each type.

## 就是对每个基因增加了pval和qval两列信息，挑选出那些在不同media培养条件下显著差异表达的基因，310个，
candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.01))
candidate_clustering_genes; length(candidate_clustering_genes)

## 计算这310个基因在不同的celltype的specificity值，一共620行，每个gene在两种细胞类型中都有各自的specificity
marker_spec <- calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)
head(selectTopMarkers(marker_spec, 3)) # 这句代码返回myoblast和fibroblast的3个特征基因
### 这个列表中的specificity的值是介于0到1之间，越接近于1提示在某一种细胞类型中越显得特异，可以用它来定义某个已知
### 细胞类型的marker genes或者用这些marker来定义新的细胞类型。

# To cluster the cells, we'll choose the top 500 markers for each of these cell types:（这句话没有理解，为什么是每种细胞类型500个gene？）
semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
plot_ordering_genes(HSMM)
## 重新挑选基因，只用黑色高亮的基因来进行聚类。

## 这里我们使用更少的基因（但是是特征基因，去噪非常好）来进行cluster
plot_pc_variance_explained(HSMM, return_all = F) # norm_method = 'log', 计算PCA，并通过碎石图查看变异解释度较大的PC，进一步去噪
HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 3, 
                        norm_method = 'log',
                        reduction_method = 'tSNE', 
                        residualModelFormulaStr="~Media + num_genes_expressed", verbose = T) 
HSMM <- clusterCells(HSMM, num_clusters=2) 
## Distance cutoff calculated to 1.02776
plot_cell_clusters(HSMM, 1, 2, color="CellType")

## 这个先不做，impute这个功能有点问题。
## Impute cell type，理论上这一步应该会把unknown和ambiguous的归类到myoblast或者fibroblast
## If you provide clusterCells with a the CellTypeHierarcy, Monocle will use it classify whole clusters, rather than just individual 
## cells. Essentially, clusterCells works exactly as before, except after the clusters are built, it counts the frequency of each cell 
## type in each cluster. When a cluster is composed of more than a certain percentage (in this case, 10%) of a certain type, all the 
## cells in the cluster are set to that type. If a cluster is composed of more than one cell type, the whole thing is marked "Ambiguous".
## If there's no cell type thats above the threshold, the cluster is marked "Unknown". Thus, Monocle helps you impute the type of each 
## cell even in the presence of missing marker data.
HSMM <- clusterCells(HSMM,
                     num_clusters=2, 
                     frequency_thresh=0.1,
                     cell_type_hierarchy=cth)
## Distance cutoff calculated to 1.02776
plot_cell_clusters(HSMM, 1, 2, color_by ="CellType",markers = c("MYF5", "ANPEP"))  # 建议查看一下CellType
table(pData(HSMM)$CellType) # CellType又出问题了，需要debug

pie <- ggplot(pData(HSMM), aes(x = factor(1), fill = factor(CellType))) +
  geom_bar(width = 1)
pie + coord_polar(theta = "y") + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank())

########################################### Pseudotime分析 ###########################################
# 主要目的是：Constructing Single Cell Trajectories

# 发育过程中细胞状态是不断变化的，monocle包利用算法学习所有基因的表达模式来把每个细胞安排到各自的发展轨迹。 
# 在大多数生物学过程中，参与的细胞通常不是同步发展的，其中很多细胞是处在transition state。只有单细胞转录组技术才能把处于该过程中
# 各个中间状态的细胞分离开来。Rather than purifying cells into discrete states experimentally, Monocle uses an algorithm to learn 
# the sequence of gene expression changes each cell must go through as part of a dynamic biological process. Monocle可以把细胞都安置
# 在这个trajectory特定的位置上面。使用Monocle的differential analysis toolkit，我们可以观察某些gene在这个trajectory上的表达变化。
# Monocle构建trajectory的算法称为reversed graph embedding。

# 什么是Pseudotime：
# Pseudotime是衡量单个细胞通过细胞分化等过程所发生的变化的指标。在许多生物过程中，细胞不能以完美的同步化方式进行。在诸如细胞分化的过程的
## 单细胞表达研究中，捕获的细胞可能在分化方面广泛分布于不同的阶段。也就是说，在同一时间捕获的细胞群体中，一些细胞可能远在前方，而另一些
# 细胞甚至可能尚未开始该过程。当你想了解随着细胞从一个状态转换到另一个状态时发生的变化顺序，这种异步会产生重大问题。追踪同时捕获的细胞表达
# 产生非常压缩的基因动力学，并且该基因表达的明显变异性将非常高。通过按照学习得到的轨迹(trajectory)的进度对每个细胞进行排序，Monocle可以缓解由于不同步
# 导致的问题。 Monocle并不追踪基因的表达随时间变化的变化，而是追踪基因表达沿学习得到的轨迹(trajectory)进展的函数，我们称之为pseudotime。
# pseudotime是一个抽象的进展单位：它就是沿着最短路径测量的一个细胞和轨迹开始之间的距离。轨迹的总长度是根据细胞在从起始状态移动到结束状态
# 时经历的转录变化总量来定义的。

# 在pseudotime中我们可以得到的两个重要的结果：
# Finding Genes that Change as a Function of Pseudotime
# Analyzing Branches in Single-Cell Trajectories

# 而monocle包里面的pseudotime分析方法正是要探究这些。整个workflow分为step1， step2，step3三步
# Step1: choose genes that define a cell’s progress: feature selection: Monocle会同时关注一些表达量很低，但是vary in interesting的非noisy genes
#        使用这些gene来order cell
# Step2: reduce data dimensionality of the data：数据降维方法是reversed graph embedding
# Step3: order cells along the trajectory (in pseudotime)：我们首先将表达数据投射到一个低维的空间当中去。Monocle假定trajectory是一个树状结构，
#        树的一端是根，另一端是叶子。Monocle的职责就是找到最佳的树型结构，这个任务称为manifold learning。一个细胞的pseudotime value是细胞
#        travel back回去的距离。

# 其中第一个步骤挑选合适的基因有3种策略，分别是：
# Ordering based on genes that differ between clusters
# Selecting genes with high dispersion across cells
# Ordering cells using known marker genes

###### Trajectory step 1: choose genes that define a cell's progress  
# 核心思想是找到set of genes increase or decrease in expression as a function of process。
# 当然我们最好不要依赖已有的知识，因为这会带来偏倚。一种方法是我们isolate在process的开始和结尾处的细胞，然后找差异表达基因
# 在这里我们使用的是简化版的筛选gene，但是Monocle明确建议使用一种称之为"dpFeature"的方法来筛选gene（见后，搜索dpFeature）

# 无监督的Pseudotime分析, 先不要执行上方impute cell type的环节
HSMM_myo <- HSMM[,pData(HSMM)$CellType == "Myoblast"]   
HSMM_myo <- estimateDispersions(HSMM_myo)
## Warning: Deprecated, use tibble::rownames_to_column() instead.
## Removing 143 outliers

# 使用不同的策略会给出不同的fData(State)：先不执行策略1，而执行策略2
## 策略1：  Ordering based on genes that differ between clusters
### find all genes that are differentially expressed in response to the switch from growth medium to differentiation medium:
### 这是一种简单的处理方法，就是比较process最早和最晚的细胞的表达基因，作为这个基因集
diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,],     # 这一步比较费时
                                      fullModelFormulaStr="~Media")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
### Choosing genes based on differential analysis of time points is often highly effective, but what if we don't have time series data? 
### If the cells are asynchronously moving through our biological process (as is usually the case), Monocle can often reconstruct their 
### trajectory from a single population captured all at the same time. 


## 策略2：Selecting genes with high dispersion across cells
disp_table <- dispersionTable(HSMM_myo)
ordering_genes <- subset(disp_table, 
                         mean_expression >= 0.5 & 
                           dispersion_empirical >= 1 * dispersion_fit)$gene_id


## 跳过策略1:Once we have a list of gene ids to be used for ordering, we need to set them in the HSMM object, because 
## the next several functions will depend on them.
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
plot_ordering_genes(HSMM_myo)
## Warning: Transformation introduced infinite values in continuous y-axis


##### Trajectory step 2: reduce data dimensionality
## 挑选变异度大的基因，如图上图所示，数据降维到可以可视化(1-2维)
HSMM_myo <- reduceDimension(HSMM_myo, max_components=2, method = 'DDRTree')


##### Trajectory step 3: order cells along the trajectory：万事具备，下面就是用orderCells函数来给细胞排排序
HSMM_myo <- orderCells(HSMM_myo)
pData(HSMM_myo)
## 排序好的细胞可以直接按照发育顺序可视化
plot_cell_trajectory(HSMM_myo, color_by="Hours")

plot_cell_trajectory(HSMM_myo, color_by = "State")

## "State" is just Monocle's term for the segment of the tree. The function below is handy for identifying the State which 
## contains most of the cells from time zero. 
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
HSMM_myo <- orderCells(HSMM_myo, root_state = GM_state(HSMM_myo))
plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime")

## If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. Sometimes it 
## can be handy to "facet" the trajectory plot so it's easier to see where each of the states are located: 分面，具体看每个state在
## pseudotime的哪一支。
plot_cell_trajectory(HSMM_myo, color_by = "State") +
  facet_wrap(~State, nrow = 1)


## 如果你的数据不是在多个时间点采样怎么办？
## And if you don't have a timeseries, you might need to set the root based on where certain marker genes are expressed, using your 
## biological knowledge of the system. For example, in this experiment, a highly proliferative population of progenitor cells are 
## generating two types of post-mitotic cells. So the root should have cells that express high levels of proliferation markers. 
## We can use the jitter plot to pick figure out which state corresponds to rapid proliferation:
blast_genes <- row.names(subset(fData(HSMM_myo),
                                gene_short_name %in% c("CCNB2", "MYOD1", "MYOG")))
plot_genes_jitter(HSMM_myo[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1) # a gene is “expressed” if there is at least one count since we set min_expr = 0.1



## 把gene映射到pseudotime上面，x轴是pseudotime
HSMM_expressed_genes <-  row.names(subset(fData(HSMM_myo),
                                          num_cells_expressed >= 10))
HSMM_filtered <- HSMM_myo[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("CDK1", "MEF2C", "MYH3")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "Hours")


##### Alternative choices for ordering genes -- for STEP 3 --Ordering based on genes that differ between clusters
## Ordering based on genes that differ between clusters(筛选gene用于order cell，我们选用dpFeature)
## dpFeature的第一步是选择那些至少在5%的细胞当中有表达的gene
HSMM_myo <- detectGenes(HSMM_myo, min_expr = 0.1) # a gene is “expressed” if there is at least one count since we set min_expr = 0.1
fData(HSMM_myo)$use_for_ordering <-
  fData(HSMM_myo)$num_cells_expressed > 0.05 * ncol(HSMM_myo)  # 这个5%好像有点太小了

## 然后我们使用PCA进行降维，通过scree Plot(碎石图), 选择变异解释度较大的PC
plot_pc_variance_explained(HSMM_myo, return_all = F)

## 然后基于筛选得到的PCA进一步数据降维，使用tSNE投射到一个二维空间；返回结果perplexity is too large。。。
HSMM_myo <- reduceDimension(HSMM_myo,
                            max_components = 2,
                            norm_method = 'log',
                            num_dim = 3,
                            reduction_method = 'tSNE',
                            verbose = T)

## This time without specifying the number of clusters; we will use thresholds on the cell’s local density (rho) and nearest distance 
## (delta) to determine the number of clusters.即此处不指定cluster数，而是通过p和delta来计算cluster
## 使用density peak clustering来发现二维空间中的tSNE，这个density peak算法是基于每个细胞周围的密度p和到距离较远的细胞的最短距离delta来决定的
## (不就是tSNE的算法核心么)，我们需要设定p和delta的阈值，来定义哪些细胞处在比较高的密度以内和距离范围以外。默认的clusterCell选择95%的p和
## delta来定义阈值。我们还可以定义cluster的数量。默认参数往往就挺好使的了。
HSMM_myo <- clusterCells(HSMM_myo, verbose = F)

### clutering结束后，看看结果
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)')
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)')

### decision plot，选择最佳p和delta
plot_rho_delta(HSMM_myo, rho_threshold = 2, delta_threshold = 4 )

### 然后我们可以重新计算clustering，使用优化后的p和delta两个参数，使用(skip_rho_sigma = T) 跳过计算 Ρ, Σ的步骤
HSMM_myo <- clusterCells(HSMM_myo,
                         rho_threshold = 2,
                         delta_threshold = 4,
                         skip_rho_sigma = T,
                         verbose = F)

### 看一下最终clustering的结果
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Cluster)')
plot_cell_clusters(HSMM_myo, color_by = 'as.factor(Hours)')

##### Alternative choices for ordering genes -- for STEP 3 --Ordering cells using known marker genes
## 使用非监督聚类方法有助于我们无偏移地进行分析，然而有时计算机会执迷于一些与你所探讨的生物学问题无关的gene，比如细胞周期相关gene，这个时候
## 使用半监督聚类方法就能很好的弥补了。方法是，首先使用CellTypeHierarchy来定义一些gene，然后通过这些gene找到co-vary的gene，最后基于这些所有
## 的gene进行非监督的聚类，所以非监督和半监督的聚类方法的差别仅仅在于挑选的是哪些gene进行ordering。
## 比如myoblast首先逃离细胞周期，经过一系列的变化，最终表达参与肌细胞收缩的关键基因。因此我们可以标记细胞周期基因cyclin B2(CCNB2)和myotubes
## gene, 比如myosin heavy chain(MYH3)

CCNB2_id <-
  row.names(subset(fData(HSMM_myo), gene_short_name == "CCNB2"))
MYH3_id <-
  row.names(subset(fData(HSMM_myo), gene_short_name == "MYH3"))

cth <- newCellTypeHierarchy()

cth <- addCellType(cth,
                   "Cycling myoblast",
                   classify_func = function(x) { x[CCNB2_id,] >= 1 })

cth <- addCellType(cth,
                   "Myotube",
                   classify_func = function(x) { x[MYH3_id,] >= 1 })

cth <- addCellType(cth,
                   "Reserve cell",
                   classify_func =
                     function(x) { x[MYH3_id,] == 0 & x[CCNB2_id,] == 0 })

HSMM_myo <- classifyCells(HSMM_myo, cth)


### 找到co-vary的gene
marker_diff <- markerDiffTable(HSMM_myo[HSMM_expressed_genes,],
                               cth,
                               cores = 1)
#semisup_clustering_genes <-
#row.names(subset(marker_diff, qval < 0.05))
semisup_clustering_genes <-
  row.names(marker_diff)[order(marker_diff$qval)][1:1000]

### 使用前1000个gene进行ordering
HSMM_myo <- setOrderingFilter(HSMM_myo, semisup_clustering_genes)
#plot_ordering_genes(HSMM_myo)
HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2,
                            method = 'DDRTree', norm_method = 'log')
HSMM_myo <- orderCells(HSMM_myo)
HSMM_myo <- orderCells(HSMM_myo, root_state = GM_state(HSMM_myo))
plot_cell_trajectory(HSMM_myo, color_by = "CellType") +
  theme(legend.position = "right")

### 我们通过观察某些标志gene的分布，可以大致评估这个trajectory是不是make sense
###  In this experiment, one of the branches corresponds to cells that successfully fuse to form myotubes, and the other to those that 
### fail to fully differentiate
HSMM_filtered <- HSMM_myo[HSMM_expressed_genes,]

my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("CDK1", "MEF2C", "MYH3")))

cds_subset <- HSMM_filtered[my_genes,]
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "Hours",
                               ncol = 1)

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

to_be_tested <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% c("MYH3", "MEF2C", "CCNB2", "TNNT1")))
cds_subset <- HSMM_myo[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res[,c("gene_short_name", "pval", "qval")]

plot_genes_in_pseudotime(cds_subset, color_by="Hours")


diff_test_res <- differentialGeneTest(HSMM_myo[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(HSMM_myo[sig_gene_names,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T)

to_be_tested <-
  row.names(subset(fData(HSMM),
                   gene_short_name %in% c("TPM1", "MYH3", "CCNB2", "GAPDH")))

cds_subset <- HSMM[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~CellType + Hours",
                                      reducedModelFormulaStr = "~Hours")
diff_test_res[,c("gene_short_name", "pval", "qval")]
plot_genes_jitter(cds_subset,
                  grouping = "Hours", color_by = "CellType", plot_trend = TRUE) +
  facet_wrap( ~ feature_label, scales= "free_y")

############################### Analyzing Branches in Single-Cell Trajectories ###############################
lung <- load_lung()
plot_cell_trajectory(lung, color_by = "Time")

BEAM_res <- BEAM(lung, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

plot_genes_branched_heatmap(lung[row.names(subset(BEAM_res,
                                                  qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)

lung_genes <- row.names(subset(fData(lung),
                               gene_short_name %in% c("Ccnd2", "Sftpb", "Pdpn")))
plot_genes_branched_pseudotime(lung[lung_genes,],
                               branch_point = 1,
                               color_by = "Time",
                               ncol = 1)
#################################################### 算法 ####################################################
# Monocole还提出了好几个算法：
## dpFeature: Selecting features from dense cell clusters
## Reversed graph embedding
## DRTree: Dimensionality Reduction via Learning a Tree
## DDRTree: discriminative dimensionality reduction via learning a tree
## Census: a normalization method to convert of single-cell mRNA transcript to relative transcript counts.
## BEAM : to test for branch-dependent gene expression by formulating the problem as a contrast between two negative binomial GLMs.
## Branch time point detection algorithm




#==================================================================================================================
#
#              "Monocle2" package (Monocle 2.8.0)  Cole Trapnell online tutorial（不要轻易run）
#
#==================================================================================================================
# Monocle的特点：
# In many single cell studies, individual cells are executing through a gene expression program in an unsynchronized manner. 
# In effect, each cell is a snapshot of the transcriptional program under study. 
# Pseudotime：
# The package Monocle provides tools for analyzing single-cell expression experiments. Monocle introduced the strategy of ordering 
# single cells in pseudotime, placing them along a trajectory corresponding to a biological process such as cell differentiation by 
# taking advantage of individual cell's asynchronous progression of those processes. Monocle orders cells by learning an explicit 
# principal graph from the single cell genomics data with advanced machine learning techniques (Reversed Graph Embedding), which 
# robustly and accurately resolves complicated biological processes.
# Monocle was originally developed to analyze dynamic biological processes such as cell differentiation


# Monocle的历史：主要是针对单细胞转录组测序数据开发的，用来找不同细胞类型或者不同细胞状态的差异表达基因。分析起始是表达矩阵，说的更直白一点
# Monocle是针对动态的细胞变化，比如细胞分化，从而进行开发的。
# 作者推荐用比较老旧的Tophat+Cufflinks流程(因为Cole Trapnell开发的)，或者RSEM, eXpress,Sailfish,等等。需要的是基于转录本的表达矩阵，
# 我一般用subjunc+featureCounts 来获取表达矩阵。
# 分为两个版本：1.0版发表于2014年的Nature Biotechnology，还用的是tophat+cufflinks组合来计算表达量， 就不过多介绍了。
# 2.0版本：发表于2017年的Nature Methods，功能也不仅仅是差异分析那么简单。还包括pseudotime,clustering分析，而且还可以进行基于转录本的
# 差异分析，并且对与超大的数据集响应更佳，其算法是BEAM (used in branch analysis) and Census (the core of relative2abs)，也单独发表了文章。
# 用了4个公共的数据来测试说明其软件的用法和优点。

# the HSMM data set, GSE52529 (ref. 1);
# the lung data set, GSE52583 (ref. 8);
# the Paul et al. data set ;
# the Olsson data set9, synapse ID syn4975060.

# 读取表达矩阵和分组信息，需要理解其定义好的一些S4对象。
# 还提出了好几个算法：

# dpFeature: Selecting features from dense cell clusters
# Reversed graph embedding
# DRTree: Dimensionality Reduction via Learning a Tree
# DDRTree: discriminative dimensionality reduction via learning a tree
# Census: a normalization method to convert of single-cell mRNA transcript to relative transcript counts.
# BEAM : to test for branch-dependent gene expression by formulating the problem as a contrast between two negative binomial GLMs.
# Branch time point detection algorithm :

library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)
library(HSMMSingleCell)
library(monocle)
library(M3Drop)

################################ STEP 0: 使用流程解析 ################################

# 载入表达矩阵并转化为CellDataSet对象
# 对表达矩阵进行基于基因和样本的过滤并可视化
# 无监督的聚类，发现细胞的亚型
# pseudotime分析，single-cell trajectories(细胞从一种状态向另一种状态的改变)
# 差异分析，从而更深入理解不同的细胞亚型和细胞状态

################################ STEP 1: The CellDataSet class ################################
# Monocle倾向于加载absolute transcript counts(from UMI experiments)，或者经过矫正后的数据(FPKM or TPM).可以接受
# 来自与10X Genomics公司的CellRange产生的数据。有一点非常值得注意的是，Although Monocle can be used with raw 
# read counts, these are not directly proportional to expression values unless you normalize them by length, 
# so some Monocle functions could produce nonsense results. If you don't have UMI counts, We recommend you 
# load up FPKM or TPM values instead of raw read counts.也就是Monocle后续的计算对原始的reads数不是很友好，不建议
# 加载原始reads。
# FPKM或者TPM的矫正后的定量值来自于Cufflinks。当然在后续的统计模型的处理时，Monocle默认使用对应于UMI的绝对定量
# 结果(基于转录本的counts矩阵)，并采用负二项分布的统计模型进行计算和下游数据的处理。如果使用FPKM或者TPM，则建议修改模型的设置。
# 在加载数据时，除了常规三张表之间的匹配关系意外，one of the columns of the featureData should be named "gene_short_name".

# 非常重要的一点，不要自说自话自己normalize数据：if you do have UMI data, you should not normalize it yourself prior 
# to creating your CellDataSet. You should also not try to convert the UMI counts to relative abundances (by 
# converting it to FPKM/TPM data). You should not use relative2abs() as discussed below in the section on 
# Converting TPM to mRNA Counts. Monocle will do all needed normalization steps internally. Normalizing it 
# yourself risks breaking some of Monocle's key steps.

# Monocle最建议加载的数据是transcript count data，尤其是基于UMI序列的绝对定量

## Monocle的S4对象--CellDataSet 
## 主要是基于 CellDataSet 对象来进行下游分析，继承自ExpressionSet对象，也是常见的3个组成：

## exprs, a numeric matrix of expression values, where rows are genes, and columns are cells
## phenoData, an AnnotatedDataFrame object, where rows are cells, and columns are cell attributes (such as cell type, 
## culture condition, day captured, etc.)
## featureData, an AnnotatedDataFrame object, where rows are features (e.g. genes), and columns are gene attributes, such as biotype, 
## gc content, etc.

# 从头创建对象，代码如下：创建对象的时候需要指定引入的表达矩阵的方法，monocle2推荐用基于转录本的counts矩阵，同时也是默认的参数 
# expressionFamily=negbinomial.size() ，如果是其它RPKM/TMP等等，需要找到对应的参数。
# do not run
# HSMM_expr_matrix <- read.table("fpkm_matrix.txt")
# HSMM_sample_sheet <- read.delim("cell_sample_sheet.txt")
# HSMM_gene_annotation <- read.delim("gene_annotations.txt")
# pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
# fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
# HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix), phenoData = pd, featureData = fd)


pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                       phenoData = pd, featureData = fd)

# FPKM/TPM values are generally log-normally distributed, while UMIs or read counts are better modeled with 
# the negative binomial. To work with count data, specify the negative binomial distribution as the 
# expressionFamily argument to newCellDataSet。FPKM/TPM服从log转换后的正态分布；count data服从负二项分布。
# 所以在统计学模型选择的时候，建议参考：http://cole-trapnell-lab.github.io/monocle-release/docs/#getting-started-with-monocle
# Do not run
# HSMM <- newCellDataSet(count_matrix,
#                       phenoData = pd,
#                       featureData = fd,
#                       expressionFamily=negbinomial.size())
#### 针对非常大的数据集，尤其是有上万个细胞，建议使用稀疏矩阵：理论基础是，大部分细胞中的基因数是0
# Using sparse matrices can help you work with huge datasets on a typical computer. We generally recommend the 
# use of sparseMatrices for most users, as it speeds up many computations even for more modestly sized datasets.
# To work with your data in a sparse format, simply provide it to Monocle as a sparse matrix from the Matrix package
# The output from a number of RNA-Seq pipelines, including CellRanger, is already in a sparseMatrix format (e.g. MTX). 
# If so, you should just pass it directly to newCellDataSet without first converting it to a dense matrix (via as.matrix(), because that may exceed your available memeory.



############################ Importing & exporting data with other packages ############################
# Monocle is able to convert Seurat objects from the package "Seurat" and SCESets from the package "scater" into 
# CellDataSet objects that Monocle can use. It's also worth noting that the function will also work with SCESets 
# from "Scran". To convert from either a Seurat object or a SCESet to a CellDataSet, execute the function 
# importCDS() as shown: Monocle可以兼容Seurat的对象，抑或是来自Scater或者Scran的对象(SingleCellExperiment),使用函数importCDS()即可解决

# Where 'data_to_be_imported' can either be a Seurat object
# or an SCESet.
importCDS(data_to_be_imported)

# We can set the parameter 'import_all' to TRUE if we'd like to import all the slots from our Seurat object or SingleCellExperiment.
# (Default is FALSE or only keep minimal dataset)；因为一个对象往往有不止一个存储槽，所以可以添加参数指定转换的存储槽。默认是最小化转化的。
importCDS(data_to_be_imported, import_all = TRUE)

# Monocle can also export data from CellDataSets to the "Seurat" and "scater" packages through the function exportCDS():
# Monocle同样也可以把数据倒入到Seurat或者Scater对象当中去。
lung <- load_lung()

# To convert to Seurat object
lung_seurat <- exportCDS(lung, 'Seurat')

# To convert to SCESet
lung_SCESet <- exportCDS(lung, 'Scater')


############################ Estimate size factors and dispersions ############################
# Size factors help us normalize for differences in mRNA recovered across cells, and "dispersion" values will 
# help us perform differential expression analysis later.
# estimateSizeFactors() and estimateDispersions() will only work, and are only needed, if you are working with 
# a CellDataSet with a negbinomial() or negbinomial.size() expression family.
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)




#===============================================================================================================
#
#                    "Monocle2" package (Monocle 2.8.0)  Dave Tang's Analysis (UMI data)
#
#===============================================================================================================
## 首先我们来了解一下这个数据集
## I will use an open dataset provided by 10x Genomics, which was generated from Peripheral Blood Mononuclear Cells (PBMCs) from a 
## healthy donor. PBMCs are primary cells, which are cells isolated directly from tissue using enzymatic or mechanical methods, with 
## relatively small amounts of RNA (around 1pg RNA per cell). The data was sequenced on an Illumina HiSeq 4000 with approximately 92,000 
## reads per 8,381 cells. This dataset was processed with Cell Ranger 2.0.1 with the parameter –cells=10000.
## 我们看到使用10X Genomics，我们每个细胞得到的平均reads数是92000个。


########################################### 读入数据，构建CellDataSet对象 ######################################

setwd('/Users/mijiarui/Monocle')
library(monocle)
library(cellrangerRkit)  # 我们使用cellrangerRkit包来将10X的下机后的单细胞数据读入R当中
packageVersion("cellrangerRkit")

# 注意，根据Dave Tang的博客，我们需要建立一个叫pbmc8k的文件夹，再建立一个叫outs的子文件夹，然后把tar的内容存放进去，这样load_cellranger_matrix
# 才能读取
my_dir <- "/Users/mijiarui/Monocle/pbmc8k"
# load data
gbm <- load_cellranger_matrix(my_dir)
class(gbm)  # 这个对象类型，也是基于ExpressionSet来进行开发的，所以可以使用exprs(), pData()和fData()来读取对象中的内容  

# check out the expression data using exprs()
# 33,694 genes across 8,381 single cells
dim(exprs(gbm))
exprs(gbm)[1:5, 1:5]   # 可以看出CellRanger的矩阵输出默认是稀疏矩阵sparse matrix; 我们这边的是UMI data，相当于counts数据，符合负二项分布

# check out the phenotypic data
dim(pData(gbm))

# check out the feature information
dim(fData(gbm))
## this table will be useful for matching gene IDs to symbols
head(fData(gbm))

## 正式构建CellDataSet对象，注意Monocle expects that the gene symbol column in the feature data is called gene_short_name 
## or else you will get a warning. 所以现对fData进行少许微调：

# rename gene symbol column
my_feat <- fData(gbm)  # 我们把featureData存到另一个变量当中去，然后修改列名，确保gene symbol的列名是gene_short_name
names(my_feat) <- c('id', 'gene_short_name')  # 如果不这么修改，在构建CellDataSet这个对象的时候，会有warning message

# no warning
my_cds <- newCellDataSet(exprs(gbm),
                         phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
                         featureData = new("AnnotatedDataFrame", data = my_feat),   # featureData输入的是经过列名修改过的my_feat
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())

my_cds
str(my_cds)
# see ?CellDataSet for more information
slotNames(my_cds)

########################################### Normalization and QC ######################################
# 计算 normalisation and variance estimation两步，这对后续差异基因的分析是必不可少的。可能会比较耗时。
my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds)

## exploratory data analysis
### detectGenes()函数会帮助我们计算一个gene在多少个细胞中表达，和一个细胞表达多少个gene，分别添加到fData和pData当中去
### 相当于一个很初步的统计
my_cds <- detectGenes(my_cds, min_expr = 0.1) # 只要有1个counts就会被认作这个gene是表达的，因为我们这边设定的最小识别阈值是0.1
head(fData(my_cds))  # num_cells_expressed列
summary(fData(my_cds)$num_cells_expressed)
head(pData(my_cds))  # num_genes_expressed列
summary(pData(my_cds)$num_genes_expressed)


### The num_cells_expressed column is a tally of the number of cells expressing a particular gene (a gene is “expressed” if there is at 
### least one count since we set min_expr = 0.1). For example if we tally ENSG00000239945 across all cells, we should get 1, as indicated 
### in the table above.
class(exprs(my_cds['ENSG00000239945',])) # 这个exprs的还是以稀疏矩阵的格式存在
sum(exprs(my_cds['ENSG00000239945',]))
sum(exprs(my_cds['ENSG00000238009',]))

### The number of genes expressed (num_genes_expressed) per cell is stored in phenoData. Note that if a gene has 10 UMIs or 1 UMI, 
### it is still tallied as 1.
head(pData(my_cds))
# sanity check for AAACCTGAGCATCATC-1
sum(exprs(my_cds)[,"AAACCTGAGCATCATC-1"]>0)
summary(pData(my_cds)$num_genes_expressed)
# standardise to Z-distribution，我们对每个细胞能够检测到的gene数的统计结果进行可视化
x <- pData(my_cds)$num_genes_expressed
x_1 <- (x - mean(x)) / sd(x)
summary(x_1)

### Histogram of genes detected for cells
library(ggplot2)
library(cowplot) # # I like the default theme of cowplot
df <- data.frame(x = x_1)
ggplot(df, aes(x)) +       # Only a few cells have an abnormal number of expressed genes, such as the ones that are 5 standard deviations away from the mean.
  geom_histogram(bins = 50) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')

### We can also add our own metadata, such as the UMI count, to the phenoData.
# add a UMI column into phenoData
# use Matrix because data is in the sparse format
pData(my_cds)$UMI <- Matrix::colSums(exprs(my_cds))  # 矩阵的纵向求和，得到每个细胞的UMI�数。在矩阵中每一个entry是一个细胞当中的一个gene的UMI个数
head(pData(my_cds))  # 在pData当中增加UMI一列
ggplot(pData(my_cds), aes(num_genes_expressed, UMI)) + geom_point()  # 我们可以滤除那些有很高UMI/gene的细胞，因为它们有可能是doublet

####################################### Clustering cells without marker genes #######################################
# However, since I don’t know of any specific markers for this dataset I will use the unsupervised approach. The first step is 
# determining a subset of genes to use for clustering; this is because not all genes are informative, such as those that are lowly 
# expressed. The approach is to select gene based on their average expression and variability across cells. The dispersionTable() 
# function calculates the mean and dispersion values (for each gene).
# 一提到聚类，脑子里第一个问题就是，筛选有效基因进行第一层次的数据降维；大部分gene是不informative的
# 我们一般会先绘制一个dispersion ～ mean expression的散点图，重点关注表达量高的gene或者HVG
disp_table <- dispersionTable(my_cds) # 得到一个矩阵，为啥只有19434行？(总基因数是33694行)
head(disp_table); dim(disp_table)

## We will select genes, which have a mean expression >= 0.1, to use in the clustering step. The setOrderingFilter() function allows us 
## to indicate which genes we want to use for clustering. The plot_ordering_genes() function plots mean expression against the empirical 
## dispersion and highlights the set of genes (as black dots) that will be used for clustering.
table(disp_table$mean_expression>=0.1) # 我们的筛选标准是mean expression >= 0.1，我们先来看一下统计结果，只有4012个gene用于cluster
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1) # 我们把满足

fData(my_cds)
my_cds <- setOrderingFilter(my_cds, unsup_clustering_genes$gene_id) # 这一步会在fData(my_cds)中增加use_for_ordering这个变量，存储T/F的结果
fData(my_cds)
plot_ordering_genes(my_cds)   # plot_ordering_genes默认将细胞以dispersion为纵轴，mean expression为横轴绘制散点图；用于后续cluster的gene黑色高亮

## It will take a lot of computational time to cluster 8,381 cells even when we have subsetted the total number of genes down to 4,012. 
## The usual approach is to perform dimension reduction (typically Principal Component Analysis [PCA]) and use the principal components 
## that explain the majority of the variance in the data. The plot_pc_variance_explained() function plots the percentage of variance 
## explained by each component based on a PCA performed on the normalised expression data. From the plot we can decide the number of 
## components to use in our downstream analyses.
## 目前单细胞测序的数据降维默认采用PCA+tSNE的方法
plot_pc_variance_explained(my_cds, return_all = FALSE)   # 如函数名，这个函数绘制碎石图，方便用户选择PC的个数用于后续聚类
### 这一步比较耗时，The goal is to identify the elbow in this plot, which can be a bit subjective.


## 开始进行聚类，使用clusterCells()函数
## The clusterCells() function performs the clustering on a CellDataSet; however, you need to input the number of requested clusters. 
## 我们先看看使用5个PC的结果；如果num_dim没有事先指定，默认值是50个PC。在这里我们只选择前5的PC进行降维去噪。
my_cds <- reduceDimension(my_cds, max_components = 2, num_dim = 5,         # 这一步比较耗时
                          reduction_method = 'tSNE', verbose = TRUE) # 我们在进程中可以看到remove noise by PCA和reduce dimension by tSNE

# perform unsupervised clustering requesting 15-1 clusters，我们最高选取15个cluster，计算机会从1-15进行计算，找到当中cluster数最合理的
my_cds <- clusterCells(my_cds, num_clusters = 15) # 这一步的结果与Dave的结果略有出入：4.797959

# cluster information is in the Cluster column，Cluster信息存储在pData当中，一共有14个cluster
head(pData(my_cds)); pData(my_cds)$Cluster

# store cluster info for comparing later
my_cluster_dim_5 <- pData(my_cds)$Cluster

# Visualization可视化
plot_cell_clusters(my_cds)

## 我们再看看使用10个PC的结果
my_cds <- reduceDimension(my_cds, max_components = 2, num_dim = 10,
                          reduction_method = 'tSNE', verbose = TRUE)
my_cds <- clusterCells(my_cds, num_clusters = 15)  # 结果还是和Dave的略有差别的，好在还是分出了14群细胞
table(pData(my_cds)$Cluster)
my_cluster_dim_10 <- pData(my_cds)$Cluster
plot_cell_clusters(my_cds)

## 怎样评估选取多少个PC比较合理呢：We’ll use the Adjusted Rand Index (ARI) to compare the clustering results. The adjustedRand() function 
## in the 'clues' package calculates the Rand Index (RI) and different variations of the ARI (HA, MA, and FM).
library(clues)
adjustedRand(as.numeric(my_cluster_dim_5), as.numeric(my_cluster_dim_10))
### The ARI (HA) suggests that the clustering of cells is quite different when including different numbers of principal components. 
### I also compared clustering results with num_dim = 4 and num_dim = 5 and while the ARI is higher than the comparison between 
### num_dim = 5 and num_dim = 10, it still implies that the cell clusters are quite different. 也就是我们在这边比较使用5个PC和10个PC进行
### cluster的结果，HA值的高低反映了两者的差别，这个值越大提示clustering得到的结果差别越大。作者另外还计算了使用4个PC进行cluster的结果，我
### 在这边不进行罗列，最后的比较4个PC和10个PC的结果得到的HA值�比5个PC和10个PC的HA值更大。

# For the rest of the analysis, I will be using num_dim = 10 and num_clusters = 15. 为了方便可重复性，我们选择PC是10个，最多cluster是15个
# include 10 dimensions，得到了14个细胞群。


######################################### Differential expression 和marker gene ##############################################
# One of the main interests of single cell transcriptomics is the identification of novel markers for cell types. We can use 
# differentialGeneTest(), which is quite a flexible function, to find marker genes by assessing each gene’s expression level across the 
# cells. The model formulae used in the function can use any column in the phenoData table, including columns that we manually add.
# 使用differentialGeneTest()这个函数。marker gene的筛选可以针对任何pData当中的column，类似于time，medium或者很多其它自己加入的属性  。

# I’ll create a vector that indicates whether a single was clustered in “Cluster 1” or not to identify genes differentially expressed in 
# cluster 1 versus the other clusters. Dave在这边将cluster1的细胞都标记出来，其余都是非cluster1，然后计算cluster1和非cluster1的差异基因

# create vector of no's
my_vector <- rep('no', nrow(pData(my_cds)))

# change status to yes if the cell was in cluster 1
table(pData(my_cds)$Cluster) # 我们看到一共有14个cluster, 其中cluster1的有1135个
my_vector[pData(my_cds)$Cluster == 1] <- rep('yes', sum(pData(my_cds)$Cluster == 1))
table(my_vector)  # 这样就将对应cluster==1的细胞的数值由默认的no改成了yes

# add vector to phenoData，放入test这一列
pData(my_cds)$test <- my_vector
head(pData(my_cds))

# I’ll perform the differential expression analysis on the subset of genes where the mean expression (across cells) was >= 0.1
length(unsup_clustering_genes$gene_id)  # 之前挑选出来用于cluster的gene(平均表达量大于0.1），共4012个，的dispersion， mean expression信息
de_cluster_one <- differentialGeneTest(my_cds[unsup_clustering_genes$gene_id,],
                                      fullModelFormulaStr = '~test',     # 这个formular提示哪两组进行比较
                                      cores = 8)
dim(de_cluster_one) # 得到一个矩阵，增加了pval和qval两列
de_cluster_one

# order by q-value
library(dplyr)
de_cluster_one %>% arrange(qval) %>% head()  # dplyr包中的arrange函数，对qval进行排序，默认从小到大，如果要降序，使用desc(qval)

# Let’s check out CTSS (ENSG00000163131), which is the most statistically significant gene, using plot_genes_fitter().
# 注意Dave的测试自己有一点问题，返回结果不是4012个gene，所以他的结果中最significant的是ENSG00000163131
plot_genes_jitter(my_cds['ENSG00000163131',], grouping = "Cluster") # 根据cluster来看jitterplot的结果
## CTSS is specific for cluster 1 but also for clusters 4 and 6. We count data hence the expression is on a discrete scale.
## 在我的数据集中CTSS高表达于3，5，8三个cluster，而非cluster1，但是为了后续分析的一致性，我们先使用这个gene

# We can use plot_cell_clusters() to highlight clusters 1, 4, and 6. （相应的修改成3，5，8）
# add another column to the phenotypic data
# where TRUE means cluster membership in clusters 1, 4, or 6 （相应的修改成3，5，8）
pData(my_cds)$my_colour <- pData(my_cds)$Cluster == 3 | pData(my_cds)$Cluster == 5 | pData(my_cds)$Cluster == 8 # 增加my_colour一列
table(pData(my_cds)$my_colour)
plot_cell_clusters(my_cds, color_by = 'my_colour')  # CTSS is overexpressed in clusters 3, 5, and 8, which are highlighted in turquoise.
## 这个结果和Dave的还是比较一致的。

############################################### Constructing Single Cell Trajectories ###############################################
# The “unit” used for the trajectory is pseudotime; a cell at the beginning of the trajectory, i.e. starting state, will have a lower 
# pseudotime than cells at the end of the trajectory, i.e. end state.
# The trajectory analysis consists of three stages:

# STEP1: Choose genes that define progress
# STEP2: Reduce the dimensionality of the data
# STEP3: Order cells in pseudotime

# For the trajectory analysis in this post, I will use only a subset of all cells (cluster 1, 4, and 6) and genes that are expressed in 
# at least 10 cells. 为了方便计算，Dave只筛选了在cluster1,4,6中的细胞，并且只纳入在大于等于10个细胞以上有表达的基因进行分析。当然在我的这个
# 数据里面就是3，5，8三个cluster，上面已经给出了。也就是说我们只分析这三个cluster当中细胞的pseudotime
expressed_genes <- row.names(subset(fData(my_cds), num_cells_expressed >= 10))

# my_colour is TRUE if a cell belongs to either cluster 1, 4, or 6 （对应修改成3，5，8）
my_cds_subset <- my_cds[expressed_genes, pData(my_cds)$my_colour]   # 我们把满足条件的细胞和基因给挑选出来了。

# 15,446 genes and 2,060 cells
my_cds_subset  # 稍有出入，这是和一开始的cluster的distance结果有出入是相关的

## We’ll use the recommended approach of ordering based on genes that differ between clusters. First, we’ll perform another subset of 
## genes, keeping only genes expressed in greater than 5% of cells. 进一步筛选在5%以上的细胞有表达的gene
my_cds_subset <- detectGenes(my_cds_subset, min_expr = 0.1)
fData(my_cds_subset)$use_for_ordering <- fData(my_cds_subset)$num_cells_expressed > 0.05 * ncol(my_cds_subset)

# how many genes are used?
table(fData(my_cds_subset)$use_for_ordering)

# PCA分析，并描绘碎石图，人为判断有多少个PC可以纳入后续的分析
plot_pc_variance_explained(my_cds_subset, return_all = FALSE)   # 耗时
## There’s the same warning message as before when running plot_pc_variance_explained(), so I’ll just use num_dim = 10 and perform 
## clustering as before (see section “Clustering cells without marker genes”) but this time without specifying the number of clusters; 
## we will use thresholds on the cell’s local density (rho) and nearest distance (delta) to determine the number of clusters. 
## 在这次更新之后（2.8.0，作者的版本是2.4.0），并没有出现warning message，在这里作者选择了10个PC。但是这次我们并不指定cluster的数目，而是
## 指定细胞局部的density(rho)和nearest distance(delta)来确定cluster的数目（特地用在了pseudotime分析上面）。
my_cds_subset <- reduceDimension(my_cds_subset,        # 耗时
                                 max_components = 2,
                                 norm_method = 'log',
                                 num_dim = 10,
                                 reduction_method = 'tSNE',
                                 verbose = TRUE)
my_cds_subset <- clusterCells(my_cds_subset, verbose = FALSE)  # Dave的结果是3.551792
plot_rho_delta(my_cds_subset, rho_threshold = 2, delta_threshold = 10)   # Scatter plot of rho versus delta.

# We’ll use rho = 2 and delta = 10 to cluster the cells again.
my_cds_subset <- clusterCells(my_cds_subset,
                              rho_threshold = 2,
                              delta_threshold = 10,
                              skip_rho_sigma = T,     # 别再计算rho和delta了
                              verbose = FALSE)
table(pData(my_cds_subset)$Cluster)    # Dave的结果是4个cluster
plot_cell_clusters(my_cds_subset)

# Now we’ll perform the differential gene expression analysis as before but across all cell clusters.
clustering_DEG_genes <- differentialGeneTest(my_cds_subset,      # 这一步相当耗时，是pairwise的计算差异表达然后进行multiple-testing的校正？
                                             fullModelFormulaStr = '~Cluster',
                                             cores = 8)
dim(clustering_DEG_genes)   
head(clustering_DEG_genes)
clustering_DEG_genes %>% arrange(qval) %>% head() # 我们还是根据qval升序，然后显示前6个；和Dave的结果差不多哟！

## We’ll use the top 1,000 most significantly differentially expressed genes as the set of ordering genes and perform the dimension 
## reduction and the trajectory analysis (using the orderCells() function). 对trajectory的分析也是基于差异基因的筛选的，而且是across all
## clusters，我们用这些gene进行降维和trajectory的分析
my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
my_cds_subset <- setOrderingFilter(my_cds_subset, ordering_genes = my_ordering_genes)
my_cds_subset <- reduceDimension(my_cds_subset, method = 'DDRTree')  # 耗时

# the warnings were for use of deprecated code
my_cds_subset <- orderCells(my_cds_subset)  # 使用orderCells函数，在pData中增加了Pseudotime一项？
pData(my_cds_subset)

plot_cell_trajectory(my_cds_subset, color_by = "Cluster")

################################### Finding Genes that Change as a Function of Pseudotime ###################################
# Once we have a trajectory, we can use differentialGeneTest() to find genes that have an expression pattern that varies according 
# to pseudotime. 也就是说trajectory构建好了以后，就可以使用的是differentialGeneTest()函数，来找到随着pseudotime发生变化的gene。

# pseudotime is now a column in the phenotypic data as well as the cell state
head(pData(my_cds_subset))

my_pseudotime_de <- differentialGeneTest(my_cds_subset,
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         cores = 8)
my_pseudotime_de %>% arrange(qval) %>% head()  # 显示pseudotime发现的差异基因，根据qval进行从小到大进行排序。


# save the top 6 genes(因为我这边用了head函数); 尤其是为啥select函数无法正确在dplyr包中使用
a <- my_pseudotime_de %>% arrange(qval) %>% head() 
a[,'id'] -> my_pseudotime_gene
plot_genes_in_pseudotime(my_cds_subset[my_pseudotime_gene,])   # State在Dave的结果里是5个，重新check一下。


################################### Clustering Genes by Pseudotemporal Expression Pattern ###################################
# 说白了就是根据pseudotime的结果，对gene进行cluster
# Monocle provides functionality for clustering genes according to their pseudotime value. The clustering analysis and plotting are 
# done in one step using the plot_pseudotime_heatmap() function; to save the clustering results, use return_heatmap = TRUE.
# cluster the top 50 genes that vary as a function of pseudotime
b <- my_pseudotime_de %>% arrange(qval) %>% head(50)
b[,'id']-> gene_to_cluster

my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[gene_to_cluster,],   # 结果和Dave的有点出入，HLA-B和HLA-C的位置
                                                 num_clusters = 3,   # 这边我们指定了number of cluster
                                                 cores = 8,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)

# The columns of the heatmap are pseudotime values binned into 100 bins. Here’s the source code:
newdata <- data.frame(Pseudotime = seq(min(pData(my_cds_subset)$Pseudotime), 
                                       max(pData(my_cds_subset)$Pseudotime), length.out = 100))
# Since we saved the results in my_pseudotime_cluster, we can extract the genes for each cluster.
# hierarchical clustering was used to cluster the genes
# we can cut the dendrogram to form the same 3 clusters as plot_pseudotime_heatmap
my_cluster <- cutree(my_pseudotime_cluster$tree_row, 3)
my_cluster

# genes in cluster 1, 这个结果和Dave的差不多
my_pseudotime_de[names(my_cluster[my_cluster == 1]),"gene_short_name"]

# genes in cluster 2
my_pseudotime_de[names(my_cluster[my_cluster == 2]),"gene_short_name"]

# genes in cluster 3
my_pseudotime_de[names(my_cluster[my_cluster == 3]),"gene_short_name"]


################################### Analyzing Branches in Single-Cell Trajectories ###################################
# In our trajectory, we have two branches, which represents cells that have alternative gene expression patterns. These represent cells 
# that have supposedly gone through different developmental paths. Monocle provides functions that allows you to identify the genes that 
# differ at a particular branch point. Here is the trajectory again.
plot_cell_trajectory(my_cds_subset, color_by = "Cluster")

# The BEAM() function takes a CellDataSet that has been ordered with orderCells() and a branch point in the trajectory. A table of genes 
# is returned with significance values that indicate whether genes have expression patterns that are branch dependent.
# warnings not shown
BEAM_res <- BEAM(my_cds_subset, branch_point = 1, cores = 8)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

# check out the results
head(BEAM_res)  

table(BEAM_res$qval < 1e-4)   # Dave的结果是15230，216

my_branched_heatmap <- plot_genes_branched_heatmap(my_cds_subset[row.names(subset(BEAM_res, qval < 1e-4)),],
                                                   branch_point = 1,
                                                   num_clusters = 4,
                                                   cores = 8,
                                                   use_gene_short_name = TRUE,
                                                   show_rownames = TRUE,
                                                   return_heatmap = TRUE)
## The heatmap shows how some genes are over-expressed or under-expressed depending on the trajectory path.

# We can return genes that belong to specific clusters that were identified by BEAM().
head(my_branched_heatmap$annotation_row)
dim(my_branched_heatmap$annotation_row)
table(my_branched_heatmap$annotation_row$Cluster)

my_row <- my_branched_heatmap$annotation_row
my_row <- data.frame(cluster = my_row$Cluster,
                     gene = row.names(my_row),
                     stringsAsFactors = FALSE)

head(my_row[my_row$cluster == 3,'gene'])

my_gene <- row.names(subset(fData(my_cds_subset),
                            gene_short_name %in% head(my_row[my_row$cluster == 3,'gene'])))

# plot genes that are expressed in a branch dependent manner
plot_genes_branched_pseudotime(my_cds_subset[my_gene,],
                               branch_point = 1,
                               ncol = 1)
## The trend lines show the expression pattern of genes along the paths formed by branch point 1.



#===============================================================================================================
#
#               "Seurat" package (Seurat 2.3.1)  Dave Tang: Merge two 10X single cell datasets
#
#===============================================================================================================
# 首先需要按要求建立文件夹，并下载两个数据集，请follow这个link：https://davetang.org/muse/2018/01/24/merging-two-10x-single-cell-datasets/
# 加载Seurat
library(Seurat)
library(ggplot2)
packageVersion('Seurat')

# 加载数据，使用Read10X()函数可以读入10X的数据，你所要做的就是指定matrix, genes, barcodes的文件路径，注意，是到GRCh38文件夹下
## 加载第一个数据集pbmc4k
pbmc4k.data <- Read10X(data.dir = "/Users/mijiarui/Seurat_SatijaLab/pbmc4k/filtered_gene_bc_matrices/GRCh38")
pbmc4k <- CreateSeuratObject(raw.data = pbmc4k.data, project = "PBMC4K")
## 加载第二个数据集pbmc8k
pbmc8k.data <- Read10X(data.dir = "/Users/mijiarui/Seurat_SatijaLab/pbmc8k/filtered_gene_bc_matrices/GRCh38")
pbmc8k <- CreateSeuratObject(raw.data = pbmc8k.data, project = "PBMC8K")

## 查看一下两个数据集，你应该看一下基因的数据是不是一致
pbmc4k
pbmc8k
4340 + 8381 # 两组数据加起来的细胞数时12721

## merge 两组数据
pbmc<- MergeSeurat(object1 = pbmc4k,
                   object2 = pbmc8k,
                   add.cell.id1 = "4K", 
                   add.cell.id2 = "8K",
                   project = "PBMC12K")
pbmc  # 再次查看新的merge过后的dataset，也是12721个细胞

## 之后是一些质量控制，我们过滤掉mitochondria基因比列比较高的细胞
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@raw.data), value = TRUE)
length(mito.genes) # 一共有13个mitochondria gene
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ]) / Matrix::colSums(pbmc@raw.data)
str(pbmc) # 目前的metadata只有nGene，nUMI和orig.ident
pbmc <- AddMetaData(object = pbmc,
                    metadata = percent.mito,
                    col.name = "percent.mito")
str(pbmc) # metadata增加了percent.mito

## cell QC，根据nGene，nUMI和percent.mito进行可视化
VlnPlot(object = pbmc,
        features.plot = c("nGene", "nUMI", "percent.mito"),
        nCol = 3)

## Dave在这边过滤掉了percent.mito高于9%的细胞：We’ll filter out cells with 9% or higher mitochondrial content (a completely arbitrary threshold that I chose).
###  可视化一下
ggplot(pbmc@meta.data, aes(nUMI, percent.mito)) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = 0.09, linetype = "dashed", colour = "red")

### 进行数据过滤，我们筛选gene在200-2500的细胞，同时线粒体基因的比例不超过9%
pbmc <- FilterCells(object = pbmc,
                    subset.names = c("nGene", "percent.mito"),
                    low.thresholds = c(200, -Inf),
                    high.thresholds = c(2500, 0.09))


## 然后我们需要进行PCA，但是PCA之前我们首先要选择进行PCA的gene，去除不必要的噪声
## Setting the y.cutoff parameter to 2 identifies genes that are more than two standard deviations away from the average dispersion 
## within a bin. Dave在这边设定的阈值是0.8(dispersion在0.8个标准差以外的gene），得到4000多个gene。
pbmc <- FindVariableGenes(object = pbmc,        # 结果存储在pbmc@var.genes里面
                          mean.function = ExpMean,
                          dispersion.function = LogVMR,
                          x.low.cutoff = 0.0125,
                          x.high.cutoff = 3,
                          y.cutoff = 0.8)
length(pbmc@var.genes)

## 然后对gene的表达量，根据UMI和线粒体基因的比例进行校正，有点费时，更要命的是非常耗内存，所以运行不下去，不过分析流程是这样的
## 也尝试过只校正nUMI仍然有内存溢出的现象。
pbmc <- ScaleData(object = pbmc,
                  vars.to.regress = c("nUMI", "percent.mito"))

## 终于得到干净的数据，可以进行PCA分析了
pbmc <- RunPCA(object = pbmc,
               pc.genes = pbmc@var.genes,
               do.print = TRUE,
               pcs.print = 1:5,
               genes.print = 5)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)   # 我们看到来自两个数据集的细胞重叠在一起，说明batch effect很小。

## 然后我们run tSNE
pbmc <- RunTSNE(object = pbmc,
                dims.use = 1:10,
                do.fast = TRUE)
TSNEPlot(object = pbmc, do.label = TRUE)  # 我们看到两团细胞还是overlap在一起的

## 最后我们将两组数据分开，看一下tSNE，不难发现分布非常一致，但是8k的确实比4k的细胞数多
my_4k_cell <- grep("^4K", row.names(pbmc@meta.data), value = TRUE)
my_8k_cell <- grep("^8K", row.names(pbmc@meta.data), value = TRUE)
my_4k_tsne <- TSNEPlot(object = pbmc, do.label = TRUE, cells.use = my_4k_cell, do.return = TRUE)
my_8k_tsne <- TSNEPlot(object = pbmc, do.label = TRUE, cells.use = my_8k_cell, do.return = TRUE)
plot_grid(
  my_4k_tsne,
  my_8k_tsne
)


#=============================================================================================
#
#         "Monocle2" package (Monocle 2.8.0)  Olsson Dataset Analysis
#
#=============================================================================================
## 研究背景：
## In this tutorial, we demonstrate how to use reversed graph embedding (RGE) in Monocle 2 to resolve the complicated haematopoiesis process 
## which involves that contains two major branch points. In regards to the branch points, HSC bifurcates into either megakaryote/erythroid 
## progenitor (MEP) or granulocyte/monocyte progenitor (GMP) and the bifurcation from GMP into either granulocyte and monocyte progenitor. 
## The reconstructed developmental trajectory is learned in four dimensions but can be visualized in two dimensions. With the trajectory 
## learned, we are able to identify genes showing a significant bifurcation pattern during each lineage bifurcation through BEAM. 
## We also provided a multi-way heatmap and multi-way kinetic curves to visualize important marker genes over the differentiation process. 
## Some additional analyses are included, which are used for the tutorial on analyzing the Paul dataset.



################################### 第一步，设定全局参数，加载包 ###################################
## turn warnings off by executing 'options(warn = -1)'
## execute helper function to identify the root cell
rm(list = ls()) # clear the environment 
options(warn=-1) # turn off warning message globally 

# helper function to set the root state correctly 
get_correct_root_state <- function(cds){
  T0_counts <- table(pData(cds)$State, pData(cds)$Type)[,"Lsk"]
  as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))])
}

suppressMessages(library(monocle)) # suppress the warnings
suppressMessages(library(stringr))
suppressMessages(library(plyr))
suppressMessages(library(netbiov))


######################################## 第二步，构建CDS对象(full dataset) ########################################
## read the expression matrix in FPKM values (Monocle可以读入FPKM矫正后的数据)
## create a CDS
## convert the FPKM values to relative census counts (但是在进行下游分析的时候，一般会把FPKM转换成census counts)
## create another CDS storing the relative census counts

# reading the exprs data and create a cell dataset:
## 按照常规构建好exprs，pd(pData或者colData, 在这里是sample_sheet)和fd(fData或者rowData, 在这里是gene_ann)
setwd('/Users/mijiarui/Nature_Biotechnology_Paper/simpleSingleCell/monocle2-rge-paper/Supplementary_scripts/Jupyter_notebooks')
hta_exprs <- read.csv("./Olsson_RSEM_SingleCellRNASeq.csv",row.names=1)
sample_sheet <- data.frame(groups = str_split_fixed(colnames(hta_exprs), "\\.+", 3), row.names = colnames(hta_exprs))
## 在这里正则表达式将列名根据"."号进行分割，新构建的sample_sheet的行名就是hta_exprs的列名，建议查看一下sample_sheet
str_split_fixed(colnames(hta_exprs), "\\.+", 3)
head(sample_sheet)
tail(sample_sheet)
gene_ann <- data.frame(gene_short_name = row.names(hta_exprs), row.names = row.names(hta_exprs)) # 在构建fd的时候必须要有gene_short_name这一列
pd <- new("AnnotatedDataFrame",data=sample_sheet)
fd <- new("AnnotatedDataFrame",data=gene_ann)

tpm_mat <- hta_exprs # 将读入的表达矩阵进行手工的CPM矫正，按照一列(一个样本/细胞)矫正，然后乘以1e6
tpm_mat <- apply(tpm_mat, 2, function(x) x / sum(x) * 1e6)

# 表达矩阵，样本注释和基因注释数据准备就绪，就可以构建CellDataSet(CDS)对象了，使用newCellDataSet函数，一开始构建的对象输入的是矫正后的数值
URMM_all_std <- newCellDataSet(as.matrix(tpm_mat),phenoData = pd,featureData =fd,
                               expressionFamily = tobit(),  # 没理解为什么用negbinomial.size()
                               lowerDetectionLimit=1)

# set up the experimental type for each cell，建议在每部执行之前查看pData(URMM_all_std)所包含的内容，说白了，这边根据实验的设计
# 追加了分组信息
pData(URMM_all_std)
pData(URMM_all_std)[, 'Type'] <- as.character(pData(URMM_all_std)[, 'groups.1']) #WT cells
pData(URMM_all_std)
pData(URMM_all_std)[453:593, 'Type'] <- paste(as.character(pData(URMM_all_std)[453:593, 'groups.1']), '_knockout', sep = '') #KO cells
pData(URMM_all_std)
pData(URMM_all_std)[594:640, 'Type'] <- paste(pData(URMM_all_std)[594:640, 'groups.1'], 
                                              pData(URMM_all_std)[594:640, 'groups.2'], 'knockout', sep = '_') #double KO cells
table(pData(URMM_all_std)$Type)

# 将FPKM/TPM这些矫正表达量转换成绝对数值RPC(mRNA per cell)是建议这么做的，采用的算法叫做census，建议将FPKM/TPM转换成RPC，这一步
# 在构建CellDataSet之前完成
# run Census to get the transcript counts；在得到FPKM/TPM矩阵后，使用relative2abs (Transform relative expression values into absolute transcript...)
# if you have FPKM/TPM data, you can still use negative binomial if you first convert your relative expression values to 
# transcript counts using relative2abs(). This often leads to much more accurate results than using tobit()
URMM_all_abs_list <- relative2abs(URMM_all_std, t_estimate = estimate_t(URMM_all_std), return_all = T, method = 'num_genes') # convert to RPC
str(URMM_all_abs_list)
URMM_all_abs <- newCellDataSet(as(URMM_all_abs_list$norm_cds, 'sparseMatrix'),               # 构建CDS, 表达矩阵是RPC
                               phenoData = new("AnnotatedDataFrame",data=pData(URMM_all_std)),
                               featureData = new("AnnotatedDataFrame",data=fData(URMM_all_std)),
                               expressionFamily = negbinomial.size(),
                               lowerDetectionLimit=1)

## 在完成CellDataSet的构建之后，计算量化因子(size factor)和散度(dispersion)都非常必要，这样可以矫正样本测序深度之间的差别，并且
## 对后续差异基因的筛选是有帮助的
URMM_all_abs <- estimateSizeFactors(URMM_all_abs)
suppressMessages(URMM_all_abs <- estimateDispersions(URMM_all_abs)) # suppress the warning produced from glm fit 

# The setOrderingFilter function marks genes that will be used for clustering in subsequent calls to clusterCells
URMM_all_abs <- setOrderingFilter(URMM_all_abs, row.names(fData(URMM_all_abs)))
colnames(fData(URMM_all_abs)) # 在fData多出来一列，叫做use_for_ordering
table(fData(URMM_all_abs)[,2]) # 当然，查看了一下这一列的内容，发现所有基因都用于ordering了。



########################################## 第三步，构建CDS对象 (WT) ##########################################
# Prepare the cds for the WT dataset
# read the annotation data for cell clustering on the wild-type(WT) data

# read data from figure 1b, data collected from the Nature website 
fig1b <- read.csv("./fig1b.txt",row.names=1, sep = '\t')
dim(fig1b)

# match up the column name in fig1b to the colnames in URMM_all_fig1b
# note that you should not run this mutliple times
URMM_all_fig1b <- URMM_all_abs[, pData(URMM_all_abs)$Type %in% c('Lsk', 'Cmp', 'Gmp', 'LK')]

fig1b_names <- colnames(fig1b)
match_id <- which(str_split_fixed(colnames(fig1b), "\\.+", 2)[, 2] %in% colnames(URMM_all_fig1b) == T) # 理解这句代码需要将下面一行代码看懂
## 相当于是正则表达式，使用"."为分隔符，然后从最左边的"."开始，将字段分割成2份，取第二份
head(colnames(fig1b)); head(str_split_fixed(colnames(fig1b), "\\.+", 2)); head(str_split_fixed(colnames(fig1b), "\\.+", 2)[, 2])


fig1b_names[match_id] <- str_split_fixed(colnames(fig1b), "\\.+", 2)[match_id, 2]
no_match_id <- which(str_split_fixed(colnames(fig1b), "\\.+", 2)[, 2] %in% colnames(URMM_all_fig1b) == F)
fig1b_names[no_match_id] <- str_split_fixed(colnames(fig1b), "\\.\\.", 2)[no_match_id, 2]
colnames(fig1b)[2:383] <- fig1b_names[2:383]

#set up the color for each experiment
cols <- c("Lsk" = "#edf8fb", "Cmp" = "#ccece6", "Gmp" = "#99d8c9", "GG1" = "#66c2a4", "IG2" = "#41ae76", "Irf8" = "#238b45", "LK" = "#005824",
          "Irf8_knockout" = "#fc8d59", "Gfi1_Irf8_knockout" = "#636363", "Gfi1_knockout" = "#dd1c77")

# assign clusters to each cell based on the clustering in the original study
pData(URMM_all_fig1b)$cluster <- 0 # 给pData中添加cluster信息
cluster_assignments <- as.numeric(fig1b[1, 2:383]) # 在fig2b矩阵中第一行的行名对应的时column_clusters-flat，相当于给每个细胞的分组信息
cluster_assignments
## 然后给每个cluster赋予其对应的cluster名称
cluster_assignments <- revalue(as.factor(cluster_assignments), c("1" = "HSCP-1", "2" = "HSCP-2", "3" = "Meg", "4" = "Eryth",
                                                                 "5" = "Multi-Lin", "6" = "MDP", "7" = "Mono", "8" = "Gran", "9" = "Myelocyte"))
cluster_assignments
names(cluster_assignments) <- colnames(fig1b[1, 2:383])
class(cluster_assignments); cluster_assignments # 这样得到一个cluster_assignment的向量，给向量赋予一个name方便取值
pData(URMM_all_fig1b)$cluster <- cluster_assignments[row.names(pData(URMM_all_fig1b))]

URMM_all_fig1b <- estimateSizeFactors(URMM_all_fig1b)
suppressMessages(URMM_all_fig1b <- estimateDispersions(URMM_all_fig1b))

########################################## 第四步，单细胞trajectory的理念的理解 ##########################################
# 在发育过程中，细胞从一个功能状态转变为另一个功能状态。不同状态的细胞呈现出动态改变的蛋白表达和代谢谱，而这一切都是基于部分
# 转录重整，即某一些基因的激活和某一些基因的沉默。Monocle introduced the strategy of using RNA-Seq for single cell trajectory 
# analysis. Rather than purifying cells into discrete states experimentally, Monocle uses an algorithm to learn the sequence 
# of gene expression changes each cell must go through as part of a dynamic biological process. Once it has learned the 
# overall "trajectory" of gene expression changes, Monocle can place each cell at its proper position in the trajectory. 
# You can then use Monocle's differential analysis toolkit to find genes regulated over the course of the trajectory
# If there are multiple outcome for the process, Monocle will reconstruct a "branched" trajectory. These branches correspond 
# to cellular "decisions", and Monocle provides powerful tools for identifying the genes affected by them and involved in 
# making them. 
# Monocle是基于一个称之为reversed graph embedding的机器学习算法来构建细胞的trajectory

# 我们来理解一些pseudotime: 在很多生物学过程中，细胞的演化并不是呈现齐头并进的模式。尤其是在细胞分化的过程中，捕获得到的细胞
# 广泛分布在不同的发育阶段，有的已经到了发育的终末期，有的在发育的前期或者刚开始。通过将这些细胞按顺序排列在经过学习的trajectory
# 上面，monocle可以很好的解决从一个细胞状态向另一个细胞状态发生改变的过程中调控的顺序变化。Monocle并不是建立这种变化和时间的函数
# 而是建立这种变化和progress of trajectory之间的函数，我们称之为pseudotime。Pseudotime的概念非常抽象，它是指一个细胞到trajectory的
# 起点之间的距离，找到这样一个最短的距离。而这个trajectory的全长是细胞从起始状态到终末状态转录本变化的总差异。

# Monocle的分析流程分为三大步：
## Step 1: Monocle非常重要的特点是在进行pseudotime的时候需要有feature selection一步，不能纳入所有基因
##        choosing genes that define progress 选择用于建立progress的基因集(发生表达上调或者下调的基因)。这一步可以
##        理解为feature seletion的一部分。在单细胞分析中，部分低表达gene可能会带来分析的噪声，但也可能蕴含了与细胞
##        状态相关的重要信息。Monocle善于分析那些variable(但是并非转录噪声)的基因。在这部分基因集选取的时候，我们既可以
##        采用完全非监督的方法，或者使用半监督的方法(整合之前的先验知识)来构建monocle的trajectory
##        在选择这个基因集的时候，尽量减少使用先验知识可以最大程度的减少偏倚，在这里monocle非常建议使用一种非监督的机器学习
##        方法，称为dpFeature。To use dpFeature, we first select superset of feature genes as genes expressed in at least 5% of all the cells.
## Step 2: reducing the dimensionality of the data 对数据进行降维处理，这一步是在得到feature selection的基因集基础之上的
##        数据降维的方法成为reverse graph embedding
## Step 3: ordering the cells in pseudotime 将细胞按照顺序排列在pseudotime上:
##        将表达数据投射到低维空间后，Monocle就可以来学习并产生这个trajectory了。Monocle假设trajectory是一个树型结构
##        一端是根，另一端是叶子。Monocle的任务就是找到最佳的树�型结构。这个方法称为manifold learning。Pseudotime vlaue
##        是从细胞的位置返回到根的距离。


# Running dpFeature for selecting ordering gene
## use clustering genes from the original paper to cluster cells
## cluster cells into five groups and select top 1000 genes as the ordering genes


#1. set ordering genes for the fig1b，将这部分gene赋值给一个对象，后续的很多分析都是要基于这个对象的。
URMM_all_fig1b <- setOrderingFilter(URMM_all_fig1b, ordering_genes = row.names(fig1b))
plot_ordering_genes(URMM_all_fig1b)
## By selecting only the high loading PCs, we effectively only focus on the more interesting biological variations.
URMM_pc_variance <- plot_pc_variance_explained(URMM_all_fig1b, return_all = T, norm_method = 'log') # 碎石图，通过肉眼看哪些PC比较重要
URMM_pc_variance  # 碎石图，通过肉眼看哪些PC比较重要，从而达到降维和去噪的作用。

#2. run reduceDimension with tSNE as the reduction_method
## We will then run reduceDimension with t-SNE as the reduction method on those top PCs and project them further down to two dimensions.
## 在使用PCA去噪的基础之上，使用tSNE进一步降维，并投射到一个二维的空间当中
set.seed(2017)
URMM_all_fig1b <- reduceDimension(URMM_all_fig1b, max_components=2, norm_method = 'log', reduction_method = 'tSNE', num_dim = 12,  verbose = F)

#3. initial run of clusterCells
## Then we can run density peak clustering to identify the clusters on the 2-D t-SNE space. The densityPeak algorithm 
## clusters cells based on each cell's local density (Ρ) and the nearest distance (Δ) of a cell to another cell with higher 
## distance. We can set a threshold for the Ρ, Δ and define any cell with a higher local density and distance than the 
## thresholds as the density peaks. Those peaks are then used to define the clusters for all cells. By default, clusterCells 
## choose 95% of Ρ and Δ to define the thresholds. We can also set a number of clusters (n) we want to cluster. In this 
## setting, we will find the top n cells with high Δ with Δ among the top 50% range. The default setting often gives good 
## clustering. 根据之前的提示和要求，将细胞分成5个cluster
URMM_all_fig1b <- clusterCells(URMM_all_fig1b, verbose = F, num_clusters = 5)

#4. check the clusters，图形化展示聚类结果
options(repr.plot.width=4, repr.plot.height=3)
plot_cell_clusters(URMM_all_fig1b, color_by = 'as.factor(Cluster)') + theme (legend.position="left", legend.title=element_blank())# show_density = F,
plot_cell_clusters(URMM_all_fig1b, color_by = 'cluster') + theme (legend.position="left", legend.title=element_blank())
plot_cell_clusters(URMM_all_fig1b, color_by = 'Type') + theme (legend.position="left", legend.title=element_blank())
plot_rho_delta(URMM_all_fig1b) # We also provide the decision plot for users to check the Ρ, Δ for each cell and decide the threshold for defining the cell clusters.

URMM_all_fig1b@expressionFamily <- negbinomial.size()
pData(URMM_all_fig1b)$Cluster <- factor(pData(URMM_all_fig1b)$Cluster)
## 当我们确定cluster的结果是有意义的，下面我们可以进行差异基因的分析来看看到底是哪些基因可以用来对细胞进行分群；
## 下面这一步非常耗时
URMM_clustering_DEG_genes <- differentialGeneTest(URMM_all_fig1b, fullModelFormulaStr = '~Cluster', cores = detectCores() - 2)
URMM_clustering_DEG_genes

# 然后我们取所有基因来作为trajectory的ordering genes. use all DEG gene from the clusters
URMM_ordering_genes <- row.names(URMM_clustering_DEG_genes)[order(URMM_clustering_DEG_genes$qval)]
## 下面几句代码有试验一下，选择前1000个为ordering gene
URMM_ordering_genes <- row.names(URMM_clustering_DEG_genes)[order(URMM_clustering_DEG_genes$qval)][1:1000]



########################################## 第五步，在野生型细胞中重构发育的trajectory ##########################################
# use the feature genes selected above to reconstruct the developmental trajectory
URMM_all_fig1b <- setOrderingFilter(URMM_all_fig1b, ordering_genes = c(URMM_ordering_genes))
URMM_all_fig1b <- reduceDimension(URMM_all_fig1b, verbose = F, scaling = T, max_components = 4, 
                                  maxIter = 100, norm_method = 'log',  lambda = 20 * ncol(URMM_all_fig1b)) 
URMM_all_fig1b <- orderCells(URMM_all_fig1b)
options(repr.plot.width=3, repr.plot.height=3)
plot_cell_trajectory(URMM_all_fig1b, color_by = 'Type')
options(repr.plot.width=8, repr.plot.height=8)
plot_cell_trajectory(URMM_all_fig1b, color_by = 'cluster') + facet_wrap(~cluster)
options(repr.plot.width=8, repr.plot.height=8)
plot_cell_trajectory(URMM_all_fig1b, color_by = 'cluster', x = 1, y = 3) + facet_wrap(~cluster)


########################################## 第六步，在所有细胞中重构发育的trajectory ##########################################
# use the same set of genes for the WT cells to reconstruct the developmental trajectory
# In agreement to the graphs shown above, two branch points are discovered
pData(URMM_all_abs)[colnames(URMM_all_fig1b), 'paper_cluster'] <- as.character(pData(URMM_all_fig1b)[, 'cluster'])

URMM_all_abs <- setOrderingFilter(URMM_all_abs, ordering_genes = URMM_ordering_genes)
URMM_all_abs <- reduceDimension(URMM_all_abs, verbose = F, scaling = T, maxIter = 100, norm_method = 'log', max_components = 4, 
                                param.gamma = 100, lambda = 14 * ncol(URMM_all_fig1b)) 
URMM_all_abs <- orderCells(URMM_all_abs)
options(repr.plot.width=6, repr.plot.height=5)
plot_cell_trajectory(URMM_all_abs, color_by = 'Type')
options(repr.plot.width=8, repr.plot.height=8)
plot_cell_trajectory(URMM_all_abs, color_by = 'Type') + facet_wrap(~paper_cluster)



########################################## 第七步，在所有细胞和野生型细胞中展示树型图 ##########################################
## Trajectories are reconstructed in 4 dimensions but can be visualized as a tree layout in two dimensions
## both the WT and full dataset include the similar branch points
## 根据细胞的种类和cluster赋予其相应的颜色，用于图形化展示。
### 所有细胞当中，根据type这个变量
table(pData(URMM_all_abs)$Type)  # Type这个分类变量下一共有9个因子
type_vec <- unique(pData(URMM_all_abs)$Type)
type_cols <- RColorBrewer::brewer.pal(9, name = 'Set1')
type_cols[6] <- "#6A3D9A"
names(type_cols) <- type_vec

### 野生型细胞当中，根据cluster这个变量
table(pData(URMM_all_fig1b)$cluster) # cluster这个分类变量下一共9个因子
cluster_vec <- unique(pData(URMM_all_fig1b)$cluster)
cluster_cols <- type_cols
cluster_cols[10] <- "#0600FC"   # 不懂这一步是干嘛的
names(cluster_cols) <- cluster_vec

## 绘图：一步法，传入表达矩阵即可 (结果与tutorial不一样)
options(repr.plot.width=6, repr.plot.height=4)
plot_complex_cell_trajectory(URMM_all_fig1b[, ], color_by = 'cluster', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3) + 
  facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = cluster_cols, name = "cluster")

options(repr.plot.width=6, repr.plot.height=5)
plot_complex_cell_trajectory(URMM_all_abs[, ], color_by = 'paper_cluster', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3) + 
  facet_wrap(~Type, nrow = 2) + scale_size(range = c(0.2, 0.2)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = cluster_cols, name = "cluster")


########################################## 第八步，将树型图展示在二维空间中 ##########################################
## 图所包含的内容是和上面一样的，只不过投射到一个二维的空间中(前两个dimension)
options(repr.plot.width=5, repr.plot.height=3)
plot_cell_trajectory(URMM_all_fig1b[, ], color_by = 'cluster', show_branch_points = T, theta = 120, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~Type, nrow = 1) + scale_size(range = c(0.2, 0.2)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = cluster_cols, name = "cluster") +theme (legend.position="right", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) + theme (legend.position="top", legend.title=element_blank()) #theme(axis.text.x = element_text(angle = 30, hjust = 1))

options(repr.plot.width=8, repr.plot.height=5)
plot_cell_trajectory(URMM_all_abs[, ], color_by = 'Type', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3, theta = 90) + facet_wrap(~Type, nrow = 2) + scale_size(range = c(0.2, 0.2)) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + scale_color_manual(values = type_cols, name = "Type")  + theme (legend.position="none", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25))


##################################### 第九步，在树型结构的每一个节点(或者叶子)观察细胞类型的分布 #####################################
## we can visualize the fraction of cell types over each state of the tree structure in a heatmap
## both heatmaps show that our method is able to put distinct cell types to different state (branch) of the learned tree structure
## 非常简单，实际就是构建一个matrix，matrix中每一个entry
state_cluster_stat <- table(pData(URMM_all_fig1b)[, c('State', 'cluster')])
state_cluster_stat # 一定要亲眼看看

state_cluster_stat <- apply(state_cluster_stat, 2, function(x) x / sum(x))
state_cluster_stat
state_cluster_stat_ordered <- t(state_cluster_stat)
state_cluster_stat

options(repr.plot.width=3, repr.plot.height=3)
pheatmap::pheatmap(state_cluster_stat_ordered, cluster_cols = F, cluster_rows = F, color = colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Greys'))(10))

state_cluster_stat <- table(pData(URMM_all_abs)[, c('State', 'paper_cluster')])

state_cluster_stat <- apply(state_cluster_stat, 2, function(x) x / sum(x))
state_cluster_stat_ordered <- t(state_cluster_stat)

options(repr.plot.width=3, repr.plot.height=3)
pheatmap::pheatmap(state_cluster_stat_ordered, cluster_cols = F, cluster_rows = F, 
                   color = colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Greys'))(10))


########################## Identify genes which can be used to define stemness or lineage score ##########################
# genes used to define stemness or lineage score are identified
# the selected genes will also be used for the Paul dataset analysis
fig1b_beam_genes_ery_meg <- BEAM(URMM_all_fig1b, branch_point = 1, verbose = F, cores = detectCores() - 2)
fig1b_beam_genes_ery_meg_ILRs <- calILRs(URMM_all_fig1b, branch_point = 1, verbose = F, cores = detectCores() - 2)

fig1b_beam_genes_ery_meg_ILRs_all <- calILRs(URMM_all_fig1b, branch_point = 1, verbose = F, cores = detectCores() - 2, return_all = T)
URMM_all_absbeam_genes_ery_meg_ILRs_all <- calILRs(URMM_all_abs, branch_point = 1, verbose = F, cores = detectCores() - 2, return_all = T)

initial_reference_val_ery <- rowMeans(fig1b_beam_genes_ery_meg_ILRs_all$str_branchA_expression_curve_matrix[, 25:30])
initial_reference_val_gmp <- rowMeans(fig1b_beam_genes_ery_meg_ILRs_all$str_branchB_expression_curve_matrix[, 25:30])

diff_reference_ery <- fig1b_beam_genes_ery_meg_ILRs_all$str_branchA_expression_curve_matrix[, 31:100] - matrix(rep(initial_reference_val_ery, 70), ncol = 70)
diff_reference_gmp <- fig1b_beam_genes_ery_meg_ILRs_all$str_branchB_expression_curve_matrix[, 31:55] - matrix(rep(initial_reference_val_gmp, 25), ncol = 25)

all_down <- apply(diff_reference_ery, 1, function(x) all(x < 0)) & apply(diff_reference_gmp, 1, function(x) all(x < 0))

URMM_all_fig1b_MEP <- URMM_all_fig1b[, pData(URMM_all_fig1b)$State %in% c(1, 5)]
URMM_all_fig1b_GMP <- URMM_all_fig1b[, pData(URMM_all_fig1b)$State %in% c(1:4)]
fig1b_pseudotime_MEP <- differentialGeneTest(URMM_all_fig1b_MEP, verbose = F, cores = detectCores() - 2)
fig1b_pseudotime_GMP <- differentialGeneTest(URMM_all_fig1b_GMP, verbose = F, cores = detectCores() - 2)


########################## Visualize the lineage score and stemness score on the tree ##########################
# color the tree by the lineage / stemness score to verify the continous transition of cell states
load('./valid_subset_GSE72857_cds2') #126 nodes in 10 dimensions
ery_meg_lineage_score <- rowMeans(fig1b_beam_genes_ery_meg_ILRs)[row.names(subset(fig1b_beam_genes_ery_meg, qval <0.01))]

positive_score_genes <- intersect(row.names(valid_subset_GSE72857_cds2), names(ery_meg_lineage_score[ery_meg_lineage_score > 0.5])) #positive genes (Ery/Meg lineage)
negtive_score_genes <- intersect(row.names(valid_subset_GSE72857_cds2), names(ery_meg_lineage_score[ery_meg_lineage_score < -0.5])) #negative genes (GMP lineage)

cell_ery_meg_lineage_score <- esApply(URMM_all_fig1b[c(positive_score_genes, negtive_score_genes), ], 2, function(x) mean(x[1:length(positive_score_genes)]) - mean(x[length(positive_score_genes):length(x)]))

pData(URMM_all_fig1b)$ery_meg_lineage_score <- -cell_ery_meg_lineage_score
options(repr.plot.width=4, repr.plot.height=4)
plot_complex_cell_trajectory(URMM_all_fig1b, color_by = 'ery_meg_lineage_score', show_branch_points = T, cell_link_size = 0.5, cell_size = 1) +
  scale_colour_gradient2() + scale_x_reverse() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + theme (legend.position="right", legend.title=element_blank()) +
  theme(legend.position="top", legend.title=element_blank())

pseudotime_degs <- intersect(row.names(subset(fig1b_pseudotime_MEP, qval < 0.01)), row.names(subset(fig1b_pseudotime_GMP, qval < 0.01)))
all_down_valid <- Reduce(intersect, list(row.names(valid_subset_GSE72857_cds2), names(all_down[all_down]), pseudotime_degs))

cell_stemness_score <- esApply(URMM_all_fig1b[all_down_valid, ], 2, function(x) mean(x))
pData(URMM_all_fig1b)$cell_stemness_score <- cell_stemness_score
options(repr.plot.width=4, repr.plot.height=4)
plot_complex_cell_trajectory(URMM_all_fig1b, color_by = 'cell_stemness_score') +  
  scale_colour_gradientn(colours = terrain.colors(10)) + theme(legend.position="top", legend.title=element_blank())




############################### Save the data for using in the Paul dataset ###############################
# all_down_valid: associate with the stemness score
# positive_score_genes, negtive_score_genes: associate with the lineage score
save(all_down_valid, positive_score_genes, negtive_score_genes, file = 'gene_set')



###############################  Create the multi-way kinetic curves as well as the heatmap ############################### 
# we can visualize the gene expression dynamics over fate commitment either with the multi-way kinetic curves or the multi-way heatmap
URMM_complex_tree_ery_meg_lineage_score_kinetic_curves_cols <- c("Ery/Meg" = "#B6A5D0", "Monocyte" = "#76C143", "Granulocyte" = "#3DC6F2")

options(repr.plot.width=6, repr.plot.height=6)
plot_multiple_branches_pseudotime(URMM_all_fig1b[c('Car1', 'Elane', 'Car2', 'Prtn3'),],
                                  branches=c(3, 4, 5), color_by = 'Branch',
                                  branches_name=c("Granulocyte", "Monocyte", "Ery/Meg"), nrow = 2, ncol = 2) +
  scale_color_manual(values = URMM_complex_tree_ery_meg_lineage_score_kinetic_curves_cols)


options(repr.plot.width=8, repr.plot.height=12)
plot_multiple_branches_heatmap(URMM_all_fig1b[unique(c(positive_score_genes, negtive_score_genes)),],
                               branches=c(3, 4, 5),
                               branches_name=c("Granulocyte", "Monocyte", "Ery/Meg"),
                               show_rownames=T,
                               num_clusters=4)







######################################## Show the regulatory network ########################################

load("./network_res")
options(repr.plot.width=4, repr.plot.height=4)
plot.netbiov(res) # create the hiearchical network















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
#                                    SingleCellExperiment
#=============================================================================================
# SingleCellExperiment是单细胞转录组数据的基础包，大部分其它单细胞转录组数据处理的包都依赖于它，就跟芯片数据
# 里面的ExpressionSet对象一样，需要拼了命的理解透，才有可能做好数据分析。











#========================================================================================================
#             比较不同的对单细胞转录组数据normalization方法, 特指测序文库大小的归一化 (生信技能树)
#========================================================================================================
# 使用CPM去除文库大小影响
# 之所以需要normalization，就是因为测序的各个细胞样品的总量不一样，所以测序数据量不一样，就是文库大小不同，
# 这个因素是肯定需要去除。最简单的就是counts per million (CPM)，所有样本的所有基因的表达量都除以各自的文库
# reads总数再乘以一百万即可。那些包装好的函数，比如在Seurat中，也是用的类似的方法，归一化，再用另一个函数去除混杂因素
# (一般miRNA-seq数据结果喜欢用这个) 代码如下：

## calc_cpm <-
##  function (expr_mat, spikes = NULL) 
##  {
##    norm_factor <- colSums(expr_mat[-spikes, ])
##    return(t(t(expr_mat)/norm_factor)) * 10^6
##  }
# 但是CPM方法有一个很严重的缺陷，那些高表达并且在细胞群体表达差异很大的基因会严重影响那些低表达基因。

# 对于同时校正基因/转录本长度的的方法：FPKM, RPKM, TPM
# RPKM, FPKM and TPM去除基因或者转录本长度影响
# 最常见的是下面3个：
# RPKM - Reads Per Kilobase Million (for single-end sequencing)
# FPKM - Fragments Per Kilobase Million (same as RPKM but for paired-end sequencing, makes sure that paired 
# ends 
# mapped to the same fragment are not counted twice)
# TPM - Transcripts Per Kilobase Million (same as RPKM, but the order of normalizations is reversed - length 
# first and sequencing depth second)
# 这些normalization方法并不适合单细胞转录组测序数据，因为有一些scRNA-seq建库方法具有3端偏好性，一般是没办法测
# 全长转录本的，所以转录本的长度跟表达量不是完全的成比例。对于UMI来定量的，肯定是不能用的。
# 对于这样的数据，需要重新转换成 reads counts 才能做下游分析。
# 需要注意的是，Monocle是可以接受RPKM，FPKM，TPM为输入的，不过也有一个转换为read counts的环节
# 另外，Smart-seq2是测定全长转录本的，是否可以使用上述方法呢？

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
# To illustrate cell QC, we consider a dataset of induced pluripotent stem cells generated from three 
# different individuals (Tung et al. 2017) in Yoav Gilad’s lab at the University of Chicago. The experiments 
# were carried out on the Fluidigm C1 platform and to facilitate the quantification both __unique molecular 
# identifiers (UMIs)__ and __ERCC spike-ins__ were used. The data files are located in the tung folder in your 
# working directory. These files are the copies of the original files made on the 15/03/16. We will use these 
# copies for reproducibility purposes.
options(stringsAsFactors = FALSE)
set.seed(1234567)
## 这里直接读取过滤好的数据，是一个SCESet对象，适用于scater包的
## http://www.biotrainee.com/jmzeng/scRNA/tung_umi.rds
setwd('/Users/mijiarui/Nature_Biotechnology_Paper/Testing_dataset')
umi <- readRDS("tung_umi.rds")   # 注意这个读入的对象是SCESet，而目前后续的操作都是建立在SingleCellExperiment这个对象上的
## We recently migrated to the SingleCellExperiment class, which uses the SummarizedExperiment class instead of the ExpressionSet class. 
## This provides a more modern interface that supports more data formats and is under active development. Try running updateSCESet or 
## toSingleCellExperiment on your SCESet object, which will convert this to a SingleCellExperiment object for you. This can be used in 
##  plotPCA as before. 
## "Scater"的开发者已经将SCESet全面升级为SingleCellExperiment，起最重要的变化是将其中的ExpressionSet Class改成了SummarizedExperiment class

umi  # 修改前，注意三个slot分别是counts, exprs(CPM后在log转换，现在标准的方法) 和log2_counts(仅仅log2转换)
str(umi)
umi <- updateSCESet(umi)
umi  # 修改后, 注意三个slot分别是counts，logcounts(CPM后在log转换，现在标准的方法)和log2_counts(仅仅log2转换)

umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
umi.qc
endog_genes <- !rowData(umi.qc)$is_feature_control
dim(exprs(umi.qc[endog_genes, ]))
umi.qc

################################################# 实践 #################################################
# 先总结，目前最最推荐的normalization方法(测序文库大小的归一化)是在CPM的基础上，在log转换，这个是scran的默认方法，也是最最推荐的方法。
# 一般经过CPM在log转换的方法，记录到SingleCellExperiment里面的logcounts这个slot当中
# 其他方法，看看就好，别较真。
# 先看看原始的表达值的分布情况，这里本来应该是对每一个样本画boxplot的，但是这里的样本数量太多了，这样的可视化效果很差， 就用PCA的方式，
# 看看这表达矩阵是否可以把样本区分开，只有那些区分度非常好的normalization方法才是最优的。
# 不过scater包提供了一个plotRLE函数，可以画出类似于样本boxplot的效果(归一化后能否柱子高度差不多)。

## Raw
plotPCA(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual",
  run_args = list(exprs_values = "log2_counts")
)


## CPM
## scater默认对表达矩阵做了cpm转换，所以可以直接提取里面的信息
plotPCA(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual",
  run_args = list(exprs_values = "logcounts")
)
# 还可以看看CPM(CPM后再log转换，现在标准的方法)和原始的log转换(log2_counts，仅仅log2转换)的表达矩阵的区别
plotRLE(
  umi.qc[endog_genes, ], 
  exprs_mats = list(Raw = "log2_counts", CPM = "logcounts"),
  exprs_logged = c(TRUE, TRUE),
  colour_by = "batch"
)

## scran: 类似于CPM然后log转换的方法，结果也自然存储到logcounts里面去，计算一个也叫size factor的东西
## scran package implements a variant on CPM specialized for single-cell data (L. Lun, Bach, and Marioni 2016). Briefly this method 
## deals with the problem of vary large numbers of zero values per cell by pooling cells together calculating a normalization factor 
## (similar to CPM) for the sum of each pool. Since each cell is found in many different pools, cell-specific factors can be 
## deconvoluted from the collection of pool-specific factors using linear algebra.
## For CPM normalisation we use scater’s calculateCPM() function. For RLE, UQ and TMM we use scater’s normaliseExprs() function. 
## For scran we use scran package to calculate size factors (it also operates on SingleCellExperiment class) and scater’s normalize() 
## to normalise the data. All these normalization functions save the results to the logcounts slot of the SCE object. 
## 这个scran package implements a variant on CPM specialized for single-cell data，所以需要特殊的代码
qclust <- quickCluster(umi.qc, min.size = 30) 
str(umi.qc)
umi.qc <- computeSumFactors(umi.qc, sizes = 15, clusters = qclust)
str(umi.qc)
umi.qc <- normalize(umi.qc)
str(umi.qc)
plotPCA(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual",
  run_args = list(exprs_values = 'logcounts')
)

plotRLE(
  umi.qc[endog_genes, ], 
  exprs_mats = list(Raw = "log2_counts", scran = "logcounts"),
  exprs_logged = c(TRUE, TRUE),
  colour_by = "batch"
)

summary(sizeFactors(umi.qc))

# TMM, RLE(size factor), Upperquantile是针对bulk RNA-seq开发的，对scRNA有可能不是很适用
## TMM 
## 需要用函数 normaliseExprs 来对SCESet对象里面的表达矩阵做TMM转换
library(edgeR)
umi.qc <- normaliseExprs(umi.qc, method = "TMM", feature_set = endog_genes) # 所有以normaliseExprs转换后的结果都会存入到normcounts中去
umi.qc   # 增加了一个新的slot: normcounts，存放的是经过TMM转换后再进行cpm转换的结果; 而logcounts存放的才是TMM的结果
normcounts(umi.qc)[1:10,1:3]
logcounts(umi.qc)[1:10,1:3]
plotPCASCE(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual",
  run_args = list(exprs_values = "logcounts")
)
plotRLE(
  umi.qc[endog_genes, ], 
  exprs_mats = list(Raw = "log2_counts", TMM = "logcounts"),
  exprs_logged = c(TRUE, TRUE),
  colour_by = "batch"
)


## Size-factor (RLE)
umi.qc <- normaliseExprs(
  umi.qc,
  method = "RLE", 
  feature_set = endog_genes
)
normcounts(umi.qc)[1:10,1:3]   # normcounts，存放的是经过RLE转换后再进行cpm转换的结果
logcounts(umi.qc)[1:10,1:3] # logcounts存放的才是RLE转换后的结果
plotPCA(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual",
  exprs_values = "logcounts"
)
plotRLE(
  umi.qc[endog_genes, ], 
  exprs_mats = list(Raw = "log2_counts", TMM = "logcounts"),
  exprs_logged = c(TRUE, TRUE),
  colour_by = "batch"
)


## Upperquantile
umi.qc <- normaliseExprs(
  umi.qc,
  method = "upperquartile", 
  feature_set = endog_genes,
  p = 0.99
)
normcounts(umi.qc)[1:10,1:3] # normcounts存放的是Upperquantile后CPM转换后的结果
logcounts(umi.qc)[1:10,1:3] # logcounts存放的是upperquantile后的结果
plotPCA(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual",
  run_args = list(exprs_values = "logcounts")
)
plotRLE(
  umi.qc[endog_genes, ], 
  exprs_mats = list(Raw = "log2_counts", TMM = "logcounts"),
  exprs_logged = c(TRUE, TRUE),
  colour_by = "batch"
)

## Downsampling
## 最后要介绍的这个去除文库大小差异的方法是从大的文库样本里面随机抽取部分reads使之文库大小缩减到跟其它文库一致。它的优点是抽样
## 过程中会造成一些基因表达量为0，这样人为创造了dropout情况，弥补了系统误差。但是有个很重要的缺点，就是每次抽样都是随机的，这样
## 结果无法重复，一般需要多次抽样保证结果的鲁棒性
## 抽样函数如下：

Down_Sample_Matrix <-
  function (expr_mat) {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
      prob <- min_lib_size/sum(x)
      return(unlist(lapply(x, function(y) {
        rbinom(1, y, prob)
      })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
  }
## 抽样后的counts矩阵赋值给SingleCellExperiment对象的新的属性。
logcounts(umi.qc) <- log2(Down_Sample_Matrix(counts(umi.qc)) + 1)  # 这一步比较耗时，对象存储到logcounts这个slot里面
plotPCA(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual",
  run_args = list(exprs_values = "logcounts")
)
plotRLE(
  umi.qc[endog_genes, ], 
  exprs_mats = list(Raw = "log2_counts", DownSample = "logcounts"),
  exprs_logged = c(TRUE, TRUE),
  colour_by = "batch"
)













## 我在这边使用University of Chicago的Tung数据集来进行演示，使用的是Fluidigm C1的平台，利用UMI和ERCC spike-in分别进行定量和数据归一化
## 如果没有这个rds对象，就自己把read counts的表达矩阵读进去，变成这个适用于scater包的SingleCellExperiment对象，代码如下；
## Load the data and annotations:
options(stringsAsFactors = FALSE)
set.seed(1234567)
setwd('/Users/mijiarui/Nature_Biotechnology_Paper/Testing_dataset')
molecules <- read.table("molecules.txt", sep = "\t")   # 这个文件是表达矩阵，包括线粒体基因和 ERCC spike-ins 的表达量，可以用来做质控
dim(molecules); head(molecules[ , 1:3])
table(grepl(pattern = '^ERCC', rownames(molecules))) # 查看下ERCC有多少个，在这个矩阵中线粒体基因也是以“ENSG”开头的，后面会指明

# 这个文件是表达矩阵涉及到的所有样本的描述信息，包括样本来源于哪个细胞，以及哪个批次。
anno <- read.table("annotation.txt", sep = "\t", header = TRUE)
anno[1:5,]

# The data consists of 3 individuals (biological replicates) and 3 replicates (technical replicates) and therefore has 9 batches in total.
table(anno$individual, anno$replicate)
umi <- SingleCellExperiment(
  assays = list(counts = as.matrix(molecules)), # 构建SingleCellExperiment对象的时候，表达矩阵可以有好几个slot，所以我们要用list，指明传入的是哪个slot
  colData = anno
)
# 构建好了SingleCellExperiment对象后一定记得看一眼
umi

# Remove genes that are not expressed in any cell:
keep_feature <- rowSums(counts(umi) > 0) > 0
umi <- umi[keep_feature, ]
dim(umi) # 经过这一步过滤，我们把在所有细胞中都不表达的gene给filter掉了

# Define control features (genes) - ERCC spike-ins and mitochondrial genes (provided by the authors):
isSpike(umi, "ERCC") <- grepl("^ERCC-", rownames(umi))
umi # 在spikeNames里面增加了ERCC一项
isSpike(umi, "MT") <- rownames(umi) %in% 
  c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888",
    "ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
    "ENSG00000198695", "ENSG00000198712", "ENSG00000198804",
    "ENSG00000198763", "ENSG00000228253", "ENSG00000198938",
    "ENSG00000198840")
umi # 在spikeNames里面又增加了MT一项

# Calculate the quality metrics:
colData(umi)
rowData(umi)
umi <- calculateQCMetrics(   # 这一步会在colData和rowData增加很多统计信息，这些信息对后续细胞的过滤非常重要
  umi,
  feature_controls = list(
    ERCC = isSpike(umi, "ERCC"), 
    MT = isSpike(umi, "MT")
  )
)
colData(umi)
rowData(umi)

########################################### Cell QC ###########################################
# Library size
# Next we consider the total number of RNA molecules detected per sample (if we were using read counts rather 
# than UMI counts this would be the total number of reads). Wells with few reads/molecules are likely to have 
# been broken or failed to capture a cell, and should thus be removed.
hist(
  umi$total_counts,
  breaks = 100
)
abline(v = 25000, col = "red")

filter_by_total_counts <- (umi$total_counts > 25000)
table(filter_by_total_counts)

# Detected genes
# In addition to ensuring sufficient sequencing depth for each sample, we also want to make sure that the 
# reads are distributed across the transcriptome. Thus, we count the total number of unique genes detected 
# in each sample.
hist(
  umi$total_features,
  breaks = 100
)
abline(v = 7000, col = "red")

filter_by_expr_features <- (umi$total_features > 7000)
table(filter_by_expr_features)

## From the plot we conclude that most cells have between 7,000-10,000 detected genes, which is normal for 
## high-depth scRNA-seq. However, this varies by experimental protocol and sequencing depth. For example, 
## droplet-based methods or samples with lower sequencing-depth typically detect fewer genes per cell. The most 
## notable feature in the above plot is the “heavy tail” on the left hand side of the distribution. If detection 
## rates were equal across the cells then the distribution should be approximately normal. Thus we remove those 
## cells in the tail of the distribution (fewer than 7,000 detected genes).

# ERCCs and MTs
# Another measure of cell quality is the ratio between ERCC spike-in RNAs and endogenous RNAs. This ratio can 
# be used to estimate the total amount of RNA in the captured cells. Cells with a high level of spike-in RNAs 
# had low starting amounts of RNA, likely due to the cell being dead or stressed which may result in the RNA 
# being degraded.
plotColData(
  umi, 
  x = "total_features",
  y = "pct_counts_MT",
  colour_by = "batch"
)

plotColData(
  umi,
  x = "total_features",
  y = "pct_counts_ERCC",
  colour_by =  "batch"
)
## The above analysis shows that majority of the cells from NA19098.r2 batch have a very high ERCC/Endo ratio. Indeed, 
## it has been shown by the authors that this batch contains cells of smaller size.
filter_by_ERCC <- umi$batch != "NA19098.r2"
table(filter_by_ERCC)
filter_by_MT <- umi$pct_counts_MT < 10
table(filter_by_MT)


########################################### Cell filtering ###########################################
# Now we can define a cell filter based on our previous analysis:
## Manuel cell filtering
umi$use <- (
    # sufficient features (genes)
    filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # sufficient endogenous RNA
    filter_by_ERCC &
    # remove cells with unusual number of reads in MT genes
    filter_by_MT
)
table(umi$use)  # 经过多步的探索（feature，counts，ERCC percentage，MT percentage）我们过滤掉一部分低质量的细胞

## Automatic cell filtering
# Another option available in scater is to conduct PCA on a set of QC metrics and then use automatic outlier 
# detection to identify potentially problematic cells.
# By default, the following metrics are used for PCA-based outlier detection:
  
# pct_counts_top_100_features
# total_features
# pct_counts_feature_controls
# n_detected_feature_controls
# log10_counts_endogenous_features
# log10_counts_feature_controls
# scater first creates a matrix where the rows represent cells and the columns represent the different QC 
# metrics. Here, the PCA plot provides a 2D representation of cells ordered by their quality metrics. 
# The outliers are then detected using methods from the mvoutlier package.

### automated 这边是有问题的，待定。。。
assay(umi, "logcounts") <- log2(counts(umi) + 1)  # 追加slot
umi <- plotPCA(
  umi,
  size_by = "total_features", 
  shape_by = "use",
  run_args = list(pca_data_input = "coldata", exprs_values = 'counts', detect_outliers = TRUE),
)
?plotPCA
table(umi$outlier)
table(colData(umi))
search()
# Compare the default, automatic and manual cell filters. Plot a Venn diagram of the outlier cells from these 
# filterings.
library(limma)
auto <- colnames(umi)[umi$outlier]
man <- colnames(umi)[!umi$use]
venn.diag <- vennCounts(
  cbind(colnames(umi) %in% auto,
        colnames(umi) %in% man)
)
vennDiagram(
  venn.diag,
  names = c("Automatic", "Manual"),
  circle.col = c("blue", "green")
)

## In addition to removing cells with poor quality, it is usually a good idea to exclude genes where we suspect 
## that technical artefacts may have skewed the results. Moreover, inspection of the gene expression profiles 
## may provide insights about how the experimental procedures could be improved.
## It is often instructive to consider the number of reads consumed by the top 50 expressed genes.
## The distributions are relatively flat indicating (but not guaranteeing!) good coverage of the full 
## transcriptome of these cells. However, there are several spike-ins in the top 15 genes which suggests a 
## greater dilution of the spike-ins may be preferrable if the experiment is to be repeated.
plotQC(umi, type = "highest-expression")

filter_genes <- apply(
  counts(umi[ , colData(umi)$use]), 
  1, 
  function(x) length(x[x > 1]) >= 2
)
rowData(umi)$use <- filter_genes
table(filter_genes)

# Dimensions of the QCed dataset (do not forget about the gene filter we defined above):
dim(umi[rowData(umi)$use, colData(umi)$use])
umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
logcounts(umi.qc) <- log2(calculateCPM(umi.qc, use_size_factors = FALSE) + 1)
# Save the data
saveRDS(umi, file = "umi.rds")

################################### Data visualization ###################################
umi <- readRDS("umi.rds")
umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
endog_genes <- !rowData(umi.qc)$is_feature_control

# Before QC：Without log-transformation:
plotPCA(
  umi[endog_genes, ],
  run_args = list(exprs_values = "counts"),
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual"
)


# Before QC: with log-transformation:
plotPCA(
  umi[endog_genes, ],
  run_args = list(exprs_values = "logcounts_raws"),
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual"
)

head(counts(umi))
head(logcounts(umi))
# After QC: without log-transformation:
plotPCA(
  as.data.frame(counts(umi.qc[endog_genes, ])),
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual"
)

plotPCA(
  umi.qc[endog_genes, ],
  ntop = 14214,
  exprs_values = "logcounts",
  colour_by = "batch",
  size_by = "total_features",
  shape_by = "individual"
)

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





#========================================================================================================
#                          比较不同的对单细胞转录组数据聚类的方法 （生信技能树）
#========================================================================================================
# 背景介绍
# 聚类之前必须要对表达矩阵进行normalization，而且要去除一些批次效应等外部因素。通过对表达矩阵的聚类，可以把细胞群体分成不同的状态，
# 解释为什么会有不同的群体。不过从计算的角度来说，聚类还是蛮复杂的，各个细胞并没有预先标记好，而且也没办法事先知道可以聚多少类。
# 尤其是在单细胞转录组数据里面有很高的噪音，基因非常多，意味着的维度很高。对这样的高维数据，需要首先进行降维去噪声，可以选择PCA或者t-SNE
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
pollen <- readRDS("pollen.rds") # 在生信技能树上下载的数据是SCESet对象，需要将其转换成SingleCellExperiment对象
pollen <- updateSCESet(pollen)
pollen # 注意一下表达矩阵的slot的名字
head(colData(pollen))
head(rowData(pollen))
table(colData(pollen)$cell_type1)
plotPCA(pollen, colour_by = "cell_type1")
### 可以看到简单的PCA也是可以区分部分细胞类型的，只不过在某些细胞相似性很高的群体区分力度不够，
### 所以需要开发新的算法来解决这个聚类的问题。


###################################################  SC3聚类  ################################################### 
## Vladimir Yu Kiselev 发明单细胞共识聚类SC3(single-cell consensus clustering)
## SC3聚类的原理:  single cell RNA-seq 分析中常用无监督的聚类来坚定细胞类型，各种聚类都需要调试很多参数，这篇文章利用了consensus matrix
## 这个矩阵进行二次分类，即多次不同参数的聚类两个细胞都被聚在同一个类中的频率，而使用不同参数的这个过程可以并行化。
## 首先创建基于相关分析（Spearman correlation）的距离矩阵，然后计算拉普拉斯算子的特征向量。接着在最初的d特征向量（d取不同的值，
## 从细胞总数的4%到7%）多次运用k-means算法。平均不同次运行的结果得到一个一致性矩阵（consensus matrix），从中可以看出两个细胞在
## 不同次运行中聚到一起的频率是多少。最后在一致性矩阵执行完整的层次聚类（hierarchical clustering with complete agglomeration），
## 得到k群（k clusters）。用于k-means和层次聚类的SC3的参数k，从2迭代到10。每一次运行SC3，都要计算剪影（silhouette），绘制一致性
## 矩阵，鉴定集群特异基因（cluster specific genes）。所有这三个方面可以帮助研究者根据经验得到最优的k和n。一旦确定出稳定的集群，
## 使用的程序将被迭代到每一个集群，用以揭示每个集群中细胞间的高变异基因，然后就可以使用这些变异基因来识别亚群。

## 准备SingleCellExperiment对象 数据给 SC3方法，先预测能聚多少个类。
pollen <- sc3_estimate_k(pollen) # 增加了metadata: sc3这个元素
str(pollen)
metadata(pollen)$sc3$k_estimation
plotPCA(pollen, colour_by = "cell_type1")
## Now we are ready to run SC3 (we also ask it to calculate biological properties of the clusters):
## 由于新的SingleCellExperiment对象pollen里面没有counts这个slot，要么我们赋予counts这个slot，要么我们在下面的参数中增加gene_filter = F
## 理论上，应该筛选基因来进行后续的聚类的，最好的方式是追加counts这一个slot
### pollen <- sc3(pollen, ks = 10, biology = TRUE, gene_filter = F) 不推荐
counts(pollen) <- logcounts(pollen)
pollen
pollen <- sc3(pollen, ks = 10, biology = TRUE) # 这里是并行计算，所以速度还可以
head(rowData(pollen))
## 可以看到SC3方法处理后的SingleCellExperiment对象的基因信息增加了5列，比较重要的是sc3_gene_filter信息，决定着该基因是否拿去聚类，
## 因为基因太多了，需要挑选
table(rowData(pollen)$sc3_gene_filter)
### 只有一半的基因被挑选去聚类了
colData(pollen) # 在coldata当中增加了一列，是细胞的cluster属性

## 后面是一些可视化
sc3_plot_consensus(pollen, k = 10, show_pdata = "cell_type1")
sc3_plot_silhouette(pollen, k = 10)
sc3_plot_expression(pollen, k = 10, show_pdata = "cell_type1")
sc3_plot_markers(pollen, k = 10, show_pdata = "cell_type1")
colData(pollen)
plotPCA(pollen, colour_by = "sc3_10_clusters")

## 还支持shiny的交互式聚类，暂时不显示
sc3_interactive(deng)
# 很明显可以看到SC3聚类的效果要好于普通的PCA


###########################################  pcaReduce  ###########################################
# use the same gene filter as in SC3
input <- logcounts(pollen[rowData(pollen)$sc3_gene_filter, ])
# run pcaReduce 1 time creating hierarchies from 1 to 30 clusters
pca.red <- PCAreduce(t(input), nbt = 1, q = 30, method = 'S')[[1]]
##  这里对2~30种类别的情况都分别对样本进行分组。
## 我们这里取只有10组的时候，这些样本是如何分组的信息来可视化。
colData(pollen)$pcaReduce <- as.character(pca.red[,32 - 10])
table(colData(pollen)$pcaReduce)  # pcaReduce相当于是cluster编号
plotPCA(pollen, colour_by = "pcaReduce")


#########################################  tSNE + kmeans  #########################################
# scater包包装了 Rtsne 和 ggplot2 来做tSNE并且可视化。
pollen <- runTSNE(pollen) # 增加了reducedDimNames(1): TSNE和colData(pollen)当中tSNE_kmeans这一列
plotTSNE(pollen, run_args = list(rand_seed = 1, return_SCESet = TRUE))
pollen; colData(pollen)
## 上面的tSNE的结果，下面用kmeans的方法进行聚类，假定是8类细胞类型。
colData(pollen)$tSNE_kmeans <- as.character(kmeans(pollen@reducedDims$TSNE, centers = 8)$clust)
plotTSNE(pollen, run_args = list(rand_seed = 1), colour_by = "tSNE_kmeans")


#####################################  SINCERA  #####################################  
# 没看懂，先跳过吧
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



#===============================================================================================================
#
#                          Velocyto on droplet data (Harvard University, 2017)
#
#===============================================================================================================

# velocyto.R安装有问题，我是根据这个帖子更改的，修改了某个文件后，不能用devtools下载github上的文件了
# https://github.com/velocyto-team/velocyto.R/issues/2
# 看一下这个文件的内容：~/.R/Makevars

