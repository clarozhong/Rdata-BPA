#### 数据获取 ----
# 加载必要库
suppressMessages({
  library(GEOquery)
  library(stringr)
  library(dplyr)
})

# 1. 数据下载 ####
setwd("D:\\lqhGroup\\R\\project\\Reproduce") # 设置工作路径

# 下载GEO数据
gset <- getGEO('GSE74986', destdir=".", AnnotGPL = F, getGPL = F)

# 获取表达矩阵
exp <- exprs(gset[[1]])  # 行为探针名，列为样本名的表达矩阵
exp[1:5, 1:5] # 查看前五行和前五列

# 提取临床信息
pdata <- pData(gset[[1]])
head(pdata)
# 提取分组信息
group_list <- str_extract(pdata$`source_name_ch1`, "HEALTHY|MODERATE_ASTHMA|SEVERE_ASTHMA")
pdata$group_list <- group_list

# 查看分组情况
table(group_list)

# 2. 筛选样本 - 只保留HEALTHY和SEVERE_ASTHMA ####
keep_samples <- pdata$group_list %in% c("HEALTHY", "SEVERE_ASTHMA")

# 筛选表达矩阵和临床信息
exp <- exp[, keep_samples]
pdata <- pdata[keep_samples, ]
group_list <- group_list[keep_samples]

# 检查筛选后的分组情况
table(group_list)

# 3. 保存数据 ####
save(gset, exp, group_list, file = "GSE74986allfile.RData")

# 输出信息
cat("数据筛选完成:\n")
cat("原始样本数:", ncol(exp), "\n")
cat("筛选后样本数:", ncol(exp_filtered), "\n")
cat("分组分布:\n")
print(table(group_list_filtered))


setwd("D:\\lqhGroup\\R\\project\\Reproduce")
rm(list = ls())
library(tidyverse)
library(ggvenn)

gene_module <- readRDS('wgcna\\combined_gene_list.rds')
DEGs <- readRDS('GSE74986_diff_genes.rds')
BPA <- readRDS('wgcna\\BPA.Rds')
ls()
head(DEGs)
head(BPA)
head(gene_module)
# 交集韦恩图 ---------------------------------------------------------------

ggvenn(
  list(gene_module = gene_module, DEGs = DEGs, BPA = BPA),
  # 下面要与列表中的命名一致
  c('gene_module', 'DEGs', 'BPA'),
  # 不展示比例
  show_percentage = F,
  fill_alpha = 0.5,
  stroke_color = NA,
  fill_color = c('#E64B35', '#4DBBD5', '#3C5488')
)

hub_gene <- Reduce(intersect, list(gene_module, DEGs, BPA))
hub_gene

saveRDS(hub_gene, file = 'hub_gene.Rds')

#2.数据检查####
##2.1箱线图查看####
exp_L=melt(exp)     #melt表示拆分数据
head(exp_L)
colnames(exp_L)=c('probe','sample','value') #设置拆分后列名
exp_L$group=rep(group_list,each=nrow(exp))  #获得其分组
# 按照group列进行排序
exp_L <- exp_L[order(exp_L$group), ]
head(exp_L)
#使用箱线图检查样本的数据分布，用ggplot2进行可视化
p=ggplot(exp_L,aes(x=sample,y=value,fill=group))+  
  geom_boxplot()+
  ggtitle("校正前数据分布")+  #标题
  theme(plot.title = element_text(hjust = 0.5)) +  #标题居中
  #theme(axis.text.x = element_blank()  #隐去x轴文字
  theme(axis.text.x=element_text(angle=90,hjust = 1,size=7)  #调整x轴名称文字属性
  )#绘制箱线图
print(p)
range(exp)
#[1]  4.385228 18.667093


##2.2PCA分布查看####
library(FactoMineR)
library(factoextra)
dat=as.data.frame(t(exp))
dat.pca <- PCA(dat, graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "point",
             col.ind = group_list,
             addEllipses = TRUE,
             palette = "jco", 
             legend.title = "Groups")
save(exp,group_list, file = "2outfile.RData")

setwd("D:\\lqhGroup\\R\\project\\Reproduce") 

# id转换 ####
# 下载数据库中的GPL平台soft文件获取对应关系
## 3.1 准备探针文件 ###
rm(list=ls())
load("GSE74986allfile.RData")

# 获取GPL平台信息
gpl <- getGEO(filename='GPL6480.soft.gz', destdir=".") 
colnames(Table(gpl))
head(Table(gpl)[,c(1,7)]) ## 选择对应的探针名与symbol名

probe2gene <- Table(gpl)[,c(1,7)]
probe2gene <- subset(probe2gene, `GENE_SYMBOL` != '') # 去除未匹配到的信息
ids <- probe2gene
head(ids)
colnames(ids) <- c("probe_id","symbol")

# 如果一个探针对应多个基因名（用///分隔），取第一个基因
ids$symbol <- sapply(strsplit(ids$symbol, "///"), function(x) x[1])

## 3.2 进行匹配 ###
library(dplyr)  # 确保加载dplyr包

# id匹配
exp <- as.data.frame(exp) %>% 
  mutate(probe_id = rownames(exp)) # 矩阵增加一列探针名

# 使用dplyr的inner_join和select
exp <- exp %>% 
  inner_join(ids, by = "probe_id") %>%
  dplyr::select(probe_id, symbol, everything())  # 明确指定使用dplyr的select

exp[1:5,1:5] # 查看前五
exp <- exp[!duplicated(exp$symbol),]  # 去除重复
rownames(exp) <- exp$symbol
exp <- exp[,-(1:2)] # 去掉第一列和第二列
exp[1:5,1:4]  # 这里得到处理完成后最终的表达矩阵

save(exp, group_list, file = "GSE74986.RData")

rm(list = ls())
setwd("D:\\lqhGroup\\R\\project\\Reproduce")
load("GSE74986.RData")

library(limma)
##4.1差异分析####
design = model.matrix(~group_list)
fit = lmFit(exp, design)
fit = eBayes(fit)
deg = topTable(fit, coef = 2, number = Inf)  # 提取所有基因的差异分析结果
colnames(deg)
## 设置 logFC 值和 P 值来标记上下调基因
logFC = 1
P.Value = 0.05
# 分组
deg$sig = ifelse((deg$P.Value < P.Value) & (deg$logFC < -logFC), "Down",
                 ifelse((deg$P.Value < P.Value) & (deg$logFC > logFC), "Up",
                        "Not-Sig"))
table(deg$sig)
head(deg)
save(deg, file = "GSE74986deg.RData")

# 提取差异基因（不区分上调或下调）
diff_genes <- rownames(deg[deg$sig %in% c("Up", "Down"), ])
# 统计上调和下调基因的数量
up_genes_count <- sum(deg$sig == "Up")
down_genes_count <- sum(deg$sig == "Down")

# 输出结果
cat("上调基因的数量: ", up_genes_count, "\n")
cat("下调基因的数量: ", down_genes_count, "\n")

# 保存为 rds 文件
saveRDS(diff_genes, file = "GSE74986_diff_genes.rds")
DEGs <- readRDS('GSE74986_diff_genes.Rds')
head(DEGs)
##4.3绘制火山图####

library(ggplot2)
threshold <- -log10(0.05)
p <- ggplot(deg, aes(x = logFC, y = -log10(P.Value))) +
  # 非显著基因（灰线以下）
  geom_point(
    data = subset(deg, -log10(P.Value) <= threshold),
    color = "#999999",
    alpha = 0.8, 
    size = 2
  ) +
  # 显著基因（灰线以上），自定义红蓝主导的渐变
  geom_point(
    data = subset(deg, -log10(P.Value) > threshold),
    aes(color = logFC),
    alpha = 0.8, 
    size = 2
  ) +
  scale_color_gradientn(
    name = "-Log10_q-value",
    # 减少黄色占比，红蓝占主导
    colours = c("#0072B2", "#0099E6", "#FFD700", "#FF6666", "#b02428"), 
    # 压缩黄色范围，扩大红蓝区间
    values = scales::rescale(c(-1.5, -0.5, 0, 0.5, 1.5), to = c(0, 1)), 
    limits = c(min(deg$logFC), max(deg$logFC))
  ) +
  xlab(expression('log'[2]*'FC')) +
  ylab(expression('-log'[10]*'(pvalue)')) +
  geom_hline(yintercept = threshold, lty = 2, col = "gray50", lwd = 0.8) +
  geom_vline(xintercept = c(-1, 1), lty = 2, col = "gray50", lwd = 0.8) + # 添加 logFC = 1 和 logFC = -1 的虚线
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key.width = unit(1, "cm")
  )

p
#热图
library(pheatmap)
group<-as.data.frame(group_list)
rownames(group)<-colnames(exp)
sig_row <- which(deg$sig %in% c("Up", "Down")) 
sig_rows <- deg[sig_row, ] 
gene<-rownames(sig_rows)  #差异表达基因
count1<-exp[gene,]  #提取差异基因的表达矩阵
color_vector <- colorRampPalette(c("#F08E59","#F7FCFD", "#76AECB"))(100) #设置渐变颜色
p1<-pheatmap(count1,
             cluster_rows = TRUE,    #是否行聚类
             cluster_cols = F,       #是否列聚类
             clustering_distance_cols = "correlation",  # 设置列聚类的距离度量
             clustering_method = "complete",  # 设置列聚类的方法，这里使用完全链接法
             treeheight_row = 20,  # 若不显示行聚类树，则为0
             legend = TRUE,
             color = color_vector,
             scale='row',
             show_rownames = F,
             show_colnames = F,
             annotation_col = group
)
# 确保group_list为因子并指定顺序（按实际情况调整levels顺序）
group_list <- factor(group_list, levels = c("HEALTHY", "SEVERE_ASTHMA")) # 替换为实际的组别名称

# 根据group_list排序样本
sample_order <- order(group_list)

# 对差异基因表达矩阵的列（样本）重新排序
count1_ordered <- count1[, sample_order]

# 调整分组注释数据框以匹配排序后的样本顺序
group_ordered <- group[sample_order, , drop = FALSE]
rownames(group_ordered) <- colnames(count1_ordered) # 确保行名与列名一致

# 绘制热图
p1 <- pheatmap(count1_ordered,
               cluster_rows = TRUE,
               cluster_cols = FALSE, # 关闭列聚类以保持排序
               treeheight_row = 20,
               legend = TRUE,
               color = color_vector,
               scale = 'row',
               show_rownames = F,
               show_colnames = F,
               annotation_col = group_ordered
)

#5.差异基因富集分析####
library(clusterProfiler)
library(ggplot2)
library (org.Hs.eg.db)
###
rm(list=ls())
load('GSE74986deg.RData')
table(deg$sig)
diff.genes<-rownames(subset(deg,sig !='Not-Sig'))

##5.1将Symbol ID转换为ENTREZ ID####
diff.df <- bitr(diff.genes,
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)
head(diff.df)

##5.2进行富集分析####
go.diff <- enrichGO(gene = diff.df$ENTREZID,
                    OrgDb = org.Hs.eg.db,
                    pAdjustMethod = 'BH',  #指定 p 值校正方法
                    pvalueCutoff =0.01,
                    qvalueCutoff = 0.05,
                    ont="all",   #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                    readable =T)

#KEGG
kk.diff <- enrichKEGG(gene =diff.df$ENTREZID,
                      organism = "hsa",
                      keyType = 'kegg',  #KEGG 富集
                      pAdjustMethod = 'BH',
                      pvalueCutoff =0.05)

View(go.diff)
View(kk.diff)
save(go.diff,kk.diff,file = 'enrich.re.RData')
#提取结果
go2 <- as.data.frame(go.diff@result)
kegg2 <- as.data.frame(kk.diff@result)

write.csv(go2, file = "GOenrich.csv", row.names = FALSE)
write.csv(kegg2, file = "KEGGenrich.csv", row.names = FALSE)

##5.3 GO富集基础可视化####
#网络图
library(enrichplot)
enrichplot::cnetplot(go.diff,circular=T,colorEdge =F)
enrichplot::cnetplot(kk.diff,circular=T,colorEdge =F)
#基础可视化（气泡图）
dotplot(go.diff, split="ONTOLOGY",showCategory = 10)+facet_grid(ONTOLOGY~., scale="free") #按类别分别提取前十的GO_Term
dotplot(kk.diff,showCategory = 15)  #展示前15的条目

setwd("D:\\lqhGroup\\R\\project\\Reproduce")
rm(list = ls())
load("GSE74986.RData")

library(tidyverse)
library(WGCNA)
dat_expr <- exp
dat_group <- data.frame(group = group_list)

# 替换元素
group_list <- ifelse(group_list == "HEALTHY", "HC", 
                     ifelse(group_list == "SEVERE_ASTHMA", "SA", group_list))

# 创建数据框
dat_group <- data.frame(group = group_list)

# 将 group 列转换为因子
dat_group$group <- factor(dat_group$group, levels = c('HC', 'SA'))

head(dat_group)
head(group_list)

# 数据处理 ---------------------------------------------------------------

# 计算各基因的方差
vars_all_gene <- apply(dat_expr, 1, var)
# 选择方差前n%的基因
gene0 <- names(vars_all_gene[vars_all_gene >= quantile(vars_all_gene, probs = 0.75)])
# 具体提取前百分之多少，以提取后基因数量在5000个左右为宜
dat_expr <- dat_expr[gene0, ] %>% t() %>% as.data.frame()

# 检查数据中是否存在缺失条目、权重低于阈值的条目和零方差的基因
res_gsg <- goodSamplesGenes(datExpr = dat_expr, verbose = 3)
# 去除不好的基因和样本
if(res_gsg$allOK) {
  print('all OK')
} else{
  if (sum(!res_gsg$goodGenes) == 0) {
    print('genes OK')
  } else{
    print(paste0('Removing gene(s): ', paste(colnames(dat_expr)[!res_gsg$goodGenes], collapse = ',')))
  }
  if (sum(!res_gsg$goodSamples) == 0) {
    print('samples OK')
  } else{
    print(paste0('Removing sample(s): ', paste(rownames(dat_expr)[!res_gsg$goodSamples], collapse = ',')))
  }
  good_genes <- colnames(expr)[res_gsg$goodGenes]
  good_samples <- rownames(expr)[res_gsg$goodSamples]
  dat_expr <- dat_expr[good_samples, good_genes]
  print('Bad gene(s) and sample(s) have been removed')
}


# 样本聚类 ---------------------------------------------------------------

# 对样本聚类
clust_sample_tree <- hclust(dist(dat_expr), method = 'average')

# 样本聚类树
pdf(file = 'wgcna/WGCNA_clust_sample_tree.pdf', width = 8, height = 6)
#源代码中width=8会产生较为严重的重叠现象，我改为了16
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(clust_sample_tree, main = 'Sample Clustering to Detect Outliers')
# 剪切线
# abline(h = 60, col = '#E64B35')
dev.off()

# 剪切样本聚类树
# cut_sample_clust <- cutreeStatic(clust_sample_tree,
#                                  # 剪切线
#                                  cutHeight = 150,
#                                  # 每个聚类所包含的最小样本数
#                                  minSize = 10)
# table(cut_sample_clust)
# dat_expr <- dat_expr[cut_sample_clust == 1, ]
# dat_group <- dat_group[cut_sample_clust == 1, ]
# 一般来说不对样本进行，除非有特别离谱的离群样本

# 构建分组信息
dat_trait <- model.matrix( ~ 0 + dat_group$group) %>% as.data.frame()
colnames(dat_trait) <- levels(dat_group$group)
rownames(dat_trait) <- dat_group$sample

# 重新对样本进行聚类
# 如果前面没有剪切样本的话，这里的聚类结果和前面是一致的
clust_sample_tree2 <- hclust(dist(dat_expr), method = 'average')

# 分组颜色
color_trait <- numbers2colors(dat_trait, colors = c('#4DBBD5', '#E64B35'))

# 样本聚类树
plotDendroAndColors(
  clust_sample_tree2,
  color_trait,
  groupLabels = colnames(dat_trait),
  main = 'Sample Dendrogram and Trait Heatmap'
)
dev.off()


# 软阈值 ---------------------------------------------------------------

# 允许多线程计算
enableWGCNAThreads()
# Mac系统不要跑上面这句，跑下面这句
# allowWGCNAThreads()

# 无标度拓扑拟合指数的软阈值
powers <- 1:20

# 软阈值的无标度拓扑分析
r2_cutoff <- 0.8
res_sft <- pickSoftThreshold(dat_expr,
                             # 软阈值
                             powerVector = powers,
                             # 无标度拟合指数
                             RsquaredCut = r2_cutoff,
                             # 详细程度
                             verbose = 5)

# 软阈值散点图
pdf(file = 'wgcna/WGCNA_softpower.pdf', width = 8, height = 6)
# 散点图两部分，排列方式的是一行两列
par(mfrow = c(1,2))
# 拟合指数与power值的散点图
plot(
  # 软阈值
  res_sft$fitIndices$Power,
  # 无缩放拓扑模型拟合
  - sign(res_sft$fitIndices$slope) * res_sft$fitIndices$SFT.R.sq,
  xlab = 'Soft Threshold (power)',
  ylab = 'Scale Free Topology Model Fit, signed R^2',
  type = 'n',
  main = 'Scale Independence'
)
text(
  # 软阈值
  res_sft$fitIndices$Power,
  # 无缩放拓扑模型拟合
  - sign(res_sft$fitIndices$slope) * res_sft$fitIndices$SFT.R.sq,
  labels = powers,
  col = '#E64B35'
)
abline(
  # 无标度拟合指数
  h = r2_cutoff, col = '#E64B35')
# 平均连通性与power值的散点图
plot(
  # 软阈值
  res_sft$fitIndices$Power,
  # 节点连接程度
  res_sft$fitIndices$mean.k.,
  xlab = 'Soft Threshold (power)',
  ylab = 'Mean Connectivity',
  type = 'n',
  main = 'Mean Connectivity'
)
text(
  # 软阈值
  res_sft$fitIndices$Power,
  # 节点连接程度
  res_sft$fitIndices$mean.k.,
  labels = powers,
  col = '#E64B35'
)
dev.off()


# 基因聚类 ---------------------------------------------------------------

# 最佳软阈值
soft_power <- res_sft$powerEstimate
print(soft_power)
soft_power <- 20
# 根据最佳软阈值计算网络邻接性
mat_adjacency <- adjacency(dat_expr, power = soft_power)

# 计算TOM矩阵（拓扑重叠矩阵，Topological Overlap Matrix）
mat_TOM <- TOMsimilarity(mat_adjacency)
# 计算TOM差异矩阵
mat_dissTOM <- 1 - mat_TOM
# TOM表示基因之间相似性，dissTOM表示基因之间的距离
# TOM越接近0 = 相似性越低 = dissTOM越接近1 = 距离越远

# 关闭多线程计算
disableWGCNAThreads()

# 对基因聚类
clust_gene_tree <- hclust(as.dist(mat_dissTOM), method = 'average')

# 基因聚类树
plot(
  clust_gene_tree,
  xlab = NA,
  sub = NA,
  main = 'Gene Clustering on TOM-Based Dissimilarity',
  labels = F
)
dev.off()

# 动态剪切模块识别
dynamic_mod <- cutreeDynamic(
  dendro = clust_gene_tree,
  distM = mat_dissTOM,
  # 树的分裂程度
  deepSplit = 2,
  # 将PAM算法应用于树型结构
  pamRespectsDendro = F,
  # 模块基因的最小数量
  minClusterSize = 50
)
table(dynamic_mod)

# 模块颜色
color_dynamic <- labels2colors(dynamic_mod)
table(color_dynamic)

# 基因聚类树
plotDendroAndColors(
  clust_gene_tree,
  color_dynamic,
  groupLabels = 'Dynamic Tree Cut',
  dendroLabels = F,
  hang = 0.03,
  addGuide = T,
  guideHang = 0.05,
  main = 'Gene Dendrogram and Module Colors'
)
dev.off()


# 相似模块聚类和合并 ---------------------------------------------------------------

# 计算模块的特征基因
list_me <- moduleEigengenes(dat_expr, colors = color_dynamic)

# 提取特征基因矩阵
dat_merged_me <- list_me$eigengenes

# 计算特征基因矩阵的相关性矩阵，即特征基因之间的差异
mat_diss_me <- 1 - WGCNA::cor(dat_merged_me)

# 对特征基因聚类
clust_me_tree <- hclust(as.dist(mat_diss_me), method = 'average')

# 特征基因聚类树
pdf(file = 'wgcna/WGCNA_eigengenes_cluster.pdf', width = 6, height = 5)
plot(clust_me_tree,
     main = 'Clustering of Module Eigengenes',
     xlab = NA,
     sub = NA)
# 合并相似模块的阈值
me_diss_threshold <- 0.4
abline(h = me_diss_threshold, col = '#E64B35')
dev.off()

# 合并相似模块
res_merge <- mergeCloseModules(dat_expr,
                               color_dynamic,
                               cutHeight = me_diss_threshold,
                               # 详细程度
                               verbose = 3)
color_merge <- res_merge$colors
dat_merge_me_new <- res_merge$newMEs

# 动态模块颜色的基因聚类树
pdf(file = 'wgcna/WGCNA_dendrogram.pdf', width = 6, height = 6)
plotDendroAndColors(
  clust_gene_tree,
  color_merge,
  groupLabels = 'Dynamic Tree Cut',
  dendroLabels = F,
  hang = 0.03,
  addGuide = T,
  guideHang = 0.05,
  main = 'Gene Dendrogram and Module Colors'
)
dev.off()

# 模块颜色
color_module <- color_merge
table(color_module)


# 模块基因热图 ---------------------------------------------------------------

# 计算相关性
mat_module_trait_cor <- WGCNA::cor(dat_merge_me_new,
                                   dat_trait,
                                   method = 'spearman',
                                   use = 'p')

# 计算p值
mat_module_trait_p <- corPvalueStudent(mat_module_trait_cor, nrow(dat_expr))

# 相关性和p值标注
mat_text <- paste0(
  'Correlation = ',
  round(mat_module_trait_cor, 2),
  '\npvalue = ',
  signif(mat_module_trait_p, 3)
)
dim(mat_text) <- dim(mat_module_trait_cor)

# 热图
pdf(file = 'wgcna/WGCNA_heatmap.pdf', width = 8, height = 8)
par(mar = c(5, 10, 3, 5))
labeledHeatmap(
  Matrix = mat_module_trait_cor,
  xLabels = colnames(dat_trait),
  xLabelsAngle = 0,
  yLabels = colnames(dat_merge_me_new),
  ySymbols = colnames(dat_merge_me_new),
  colors = colorRampPalette(colors = c('#4DBBD5', 'white', '#E64B35'))(50),
  textMatrix = mat_text,
  setStdMargins = F,
  zlim = c(-1, 1),
  main = 'Module-Trait Relationship'
)
dev.off()


# 模块基因 ---------------------------------------------------------------

# 挑选模块，注意不要挑选到'grey'模块
# 注：'grey60','lightgrey'等看着好像是'grey'模块的，但是都不是'grey'模块，都可以选
# 可以挑选多个模块基因，最后合并
select_module <- 'yellow'

# 批量选择模块基因
gene_module <- gene0[color_module %in% select_module]

saveRDS(gene_module, file = 'wgcna/gene_module_yellow.Rds')

# 计算模块与基因的相关性
mat_gene_module_cor <- WGCNA::cor(dat_expr, dat_merge_me_new, use = 'p')
mat_gene_module_p <- corPvalueStudent(mat_gene_module_cor, nrow(dat_expr))

# 计算分组与基因的相关性
mat_gene_trait_cor <- WGCNA::cor(dat_expr, dat_trait, use = 'p')
mat_gene_trait_p <- corPvalueStudent(mat_gene_trait_cor, nrow(dat_expr))

# 散点图
pdf(file = 'wgcna/WGCNA_scatterplot_yellow.pdf', width = 8, height = 8)
verboseScatterplot(
  abs(mat_gene_module_cor[gene_module, paste0('ME', select_module)]),
  abs(mat_gene_trait_cor[gene_module, 1]),
  xlab = paste0('Module Membership in ',select_module,' Module'),
  ylab = 'Gene Significance for Trait',
  main = 'Module membership vs Gene Significance\n',
  col = select_module
)
abline(h = 0.5, col = select_module)
abline(v = 0.8, col = select_module)
dev.off()

setwd("D:\\lqhGroup\\R\\project\\Reproduce\\wgcna")
#wgcna合并
# 读取两个RDS文件
genes_yellow <- readRDS("gene_module_yellow.Rds")
genes_salmon <- readRDS("gene_module_salmon.Rds")

# 合并并去除重复
combined_genes_unique <- unique(c(genes_yellow, genes_salmon))

# 保存合并后的结果
saveRDS(combined_genes_unique, file = "combined_gene_list.Rds")
#BPA靶点
# 创建基因列表
genes <- c("AR", "ESR1", "ESR1", "HTR6", "ESRRG", "ALOX5", "CA2", "CA4", "MAPT", 
           "RPS6KA3", "PTGS1", "BCL2L1", "BCL2", "ALOX15", "ALOX12", "RORC", 
           "HTR2B", "SLC6A2", "GABRG2", "GABRA1", "GABRB2", "PTGS2", "GABRB3", 
           "SHBG", "SRD5A1", "SRD5A2", "DAO", "ADRA2A", "ADRA2C", "ADRA2B", 
           "CHRM1", "SLC6A4", "TACR2", "SLC6A3", "CHRM3", "ADORA3", "LTA4H", 
           "PLEC", "DYRK1A", "MAPK1", "MAPK3", "IL1B", "IL6", "TNF", "PPARG", 
           "CAT", "SREBF1", "BAX", "FASN", "PGR", "STAR", "CYP11A1", "ACACA", 
           "AKT1", "CASP3", "DDIT3", "VTG1", "CASP9", "CYP19A1", "CYP19A1B", 
           "HSPA5", "CEBPA", "NOS2", "PCNA", "DNMT1", "SCD1", "SOD1", "ATF4", 
           "FABP4", "LPL", "PPARA", "CD36", "CDH1", "MMP9", "TP53", "HMOX1", 
           "S100G", "CASP8", "CCND1", "CYP17A1", "DNMT3A", "ESR2B", "FAS", 
           "IGF1", "SCD", "SOD2", "ATF6", "CYP1B1", "EIF2AK3", "GREB1", "INS1", 
           "KISS1", "LHB", "NFE2L2", "ODC1", "TRIB3", "CPT1A", "DDIT4", "FOS", 
           "GPER1", "GPX1", "HMGCR", "IL18", "MMP2", "NR1I2", "OCLN", "SLC2A4", 
           "THRB", "APOA1", "AQP4", "BAD", "CDH2", "CFD", "CYP1A1", "DHRS3", 
           "ESRRA", "FSHB", "GPAM", "GPT", "HSD3B1", "IL10", "JUN", "MTOR", 
           "PER2", "PRL", "SGK1", "SLC7A5", "SQLE", "TFF1", "TJP1")

# 去除重复基因（如果需要）
unique_genes <- unique(genes)

# 保存为RDS文件
saveRDS(unique_genes, file = "BPA.rds")

# 验证文件
read_genes <- readRDS("BPA.rds")
print(read_genes)

rm(list = ls())
library(tidyverse)
library(ggvenn)
setwd("D:\\lqhGroup\\R\\project\\Reproduce")
hub_gene_lasso <- readRDS('lasso/hub_gene_lasso.Rds')
hub_gene_rf <- readRDS('RF/hub_gene_rf.Rds')
hub_gene_xgb <- readRDS('XGB/hub_gene_xgb.Rds')

# 三种机器学习算法交集韦恩图 ---------------------------------------------------------------

ggvenn(
  list(LASSO = hub_gene_lasso, RF = hub_gene_rf, XGB = hub_gene_xgb),
  # 下面要与列表中的命名一致
  c('LASSO', 'RF', 'XGB'),
  # 不展示比例
  show_percentage = F,
  fill_alpha = 0.5,
  stroke_color = NA,
  fill_color = c('#E64B35', '#4DBBD5', '#3C5488')
)

key_gene <- Reduce(intersect, list(hub_gene_lasso, hub_gene_rf, hub_gene_xgb))
key_gene

saveRDS(key_gene, file = 'key_gene.Rds')

rm(list = ls())
library(tidyverse)
library(rms)
library(forestplot)
library(ggDCA)
setwd("D:\\lqhGroup\\R\\project\\Reproduce")
load("GSE74986.RData")
dat_expr_GSE74986 <- exp
group_list <- ifelse(group_list == "HEALTHY", "HC", 
                     ifelse(group_list == "SEVERE_ASTHMA", "SA", group_list))

# 创建数据框
dat_group_GSE74986 <- data.frame(group = group_list)

# 将 group 列转换为因子
dat_group_GSE74986$group <- factor(dat_group_GSE74986$group, levels = c('HC', 'SA'))
rm(exp, group_list)
load("GSE43696.RData")
ls()
dat_expr_GSE43696 <- exp
group_list <- ifelse(group_list == "Control", "HC", 
                     ifelse(group_list == "Severe Asthma", "SA", group_list))

# 创建数据框
dat_group_GSE43696 <- data.frame(group = group_list)

# 将 group 列转换为因子
dat_group_GSE43696$group <- factor(dat_group_GSE43696$group, levels = c('HC', 'SA'))
rm(exp, group_list)
key_genes <- readRDS('key_gene.Rds')
head(key_genes)
dat_expr_train <- dat_expr_GSE74986
dat_group_train <-dat_group_GSE74986
dat_expr_vali <- dat_expr_GSE43696
dat_group_vali <-dat_group_GSE43696
ls()
# 整理数据
dat_expr_train_key <- dat_expr_train[key_genes, ] %>%
  t() %>%                   # 转置使样本为行
  as.data.frame() %>%
  na.omit()                 # 删除含NA的样本

# 添加分组信息（确保样本顺序一致）
dat_expr_train_key$group <- ifelse(dat_group_train$group == 'SA', 1, 0)

# 验证数据维度
nrow(dat_expr_train_key) == nrow(dat_group_train)  # 应为TRUE

# 将数据打包
ddist <- datadist(dat_expr_train_key)
options(datadist = 'ddist')

# 4A 诊断列线图 --------------------------------------------------------------

# 用lrm方法再构建一个logistic回归模型
fit_log_lrm = lrm(as.formula(paste0(
  # 公式为所有筛选的基因与分组之间的关系
  'group ~ ', paste(key_genes, collapse = ' + ')
)), data = dat_expr_train_key, x = T, y = T)

# 诊断列线图
nomo <- nomogram(
  # 这个函数绘制诊断列线图要用lrm构建的回归模型
  fit_log_lrm,
  # 进行Logit转换
  fun = plogis,
  # 概率坐标轴刻度
  fun.at = c(0.01, 0.1, 0.5, 0.9, 0.99),
  # 显示预测值
  lp = T,
  funlabel = 'Risk'
)

pdf(file = '诊断模型/列线.pdf', width = 12, height = 10)
# 网络颜色
plot(nomo, col.grid = c('grey50', 'lightgrey'))
dev.off()
# 4B 诊断校准曲线 ---------------------------------------------------------------

# 校准分析
set.seed(2024)
dat_cal <- calibrate(
  # 这个函数绘制诊断校准曲线要用lrm构建的回归模型
  fit_log_lrm,
  # 抽样方法
  method = 'boot',
  # 抽样次数
  B = 500)

# 提取结果
dat_cal <- dat_cal[, 3:1] %>% as.data.frame()
colnames(dat_cal) <- c('bias_actual', 'apparent_actual', 'pre')

# 整理结果
dat_cal <- data.frame(
  actual = c(dat_cal$bias_actual, dat_cal$apparent_actual),
  pre = rep(dat_cal$pre, 2),
  group = c(rep('Bias Corrected', nrow(dat_cal)), rep('Apparent', nrow(dat_cal)))
)

# 校准曲线
p <- ggplot(dat_cal, aes(pre, actual, color = group)) +
  # 添加对角线
  geom_abline(
    # 斜率
    slope = 1,
    # 截距
    intercept = 0,
    color = 'grey',
    lty = 2
  ) +
  geom_line() +
  theme_bw() +
  theme(legend.title = element_blank()) +
  labs(x = 'Predicted Probability', y = 'Actual Probability') +
  # 坐标轴范围
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  # 固定xy轴比例
  coord_fixed() +
  scale_color_manual(values = c('#4DBBD5', '#E64B35'))

ggsave(file = '诊断模型/校准曲线.pdf', p, height = 8, width = 8)


# 4C 诊断DCA ---------------------------------------------------------------
# 1.16版本的R包data.table与可能与R包ggDCA不兼容，如果发生这种情况需要卸载并重新安装R包data.table
#remove.packages('data.table')
#devtools::install_version('data.table',version = '1.14.10')

# 表达矩阵的分组改为字符串
dat_expr_train_key$group <- as.character(dat_expr_train_key$group)

# DCA分析
p <- dca(fit_log_lrm) %>%
  # DCA曲线
  ggplot(aes(thresholds, NB, color = model)) +
  theme_bw() +
  scale_color_manual(
    values = c('#4DBBD5', '#E64B35', '#3C5488'),
    labels = c('Logistic Model', 'All', 'None')
  ) +
  scale_linetype_manual(values = rep(1, 3)) +
  # 不展示线型图例
  guides(lty = 'none') +
  labs(color = element_blank())

ggsave(file = '诊断模型/DCA.pdf', p, height = 6, width = 7)

rm(list = ls())
library(pROC)
library(xgboost)
library(tidyverse)
library(SHAPforxgboost)
setwd("D:\\lqhGroup\\R\\project\\Reproduce")
load("GSE74986.RData")
dat_group <- data.frame(group = group_list)
dat_expr <- exp
head(dat_group)
head(exp)
# 替换元素
group_list <- ifelse(group_list == "HEALTHY", "HC", 
                     ifelse(group_list == "SEVERE_ASTHMA", "SA", group_list))

# 创建数据框
dat_group <- data.frame(group = group_list)

# 将 group 列转换为因子
dat_group$group <- factor(dat_group$group, levels = c('HC', 'SA'))
head(dat_group)

key_gene <- readRDS('key_gene.Rds')

# 诊断ROC ---------------------------------------------------------------

# 以SCD为例
# 整理分组数据
dat_expr_HSPA5 <- dat_expr['SCD', ] %>% t() %>% as.data.frame()
colnames(dat_expr_HSPA5) <- 'SCD'
dat_expr_HSPA5$group <- dat_group$group

# ROC分析
roc_res_HSPA5 <- roc(group ~ SCD,
                   data = dat_expr_HSPA5,
                   # 计算AUC和置信区间
                   auc = T,
                   ci = T)

# 整理敏感性和特异性数据
dat_roc_plot_HSPA5 <- data.frame(specificity = roc_res_HSPA5$specificities,
                               sensitivity = roc_res_HSPA5$sensitivities) %>%
  dplyr::arrange(desc(specificity), sensitivity)

# ROC曲线
p <- ggplot(dat_roc_plot_HSPA5, aes(1 - specificity, sensitivity)) +
  geom_line(color = '#4DBBD5') +
  # 添加对角线
  geom_abline(
    # 斜率
    slope = 1,
    # 截距
    intercept = 0,
    color = 'grey',
    lty = 'dashed'
  ) +
  # 展示曲线下面积
  geom_area(fill = '#4DBBD5', alpha = 0.2) +
  # 展示AUC值和置信区间
  annotate(
    'text',
    label = paste0(
      'AUC = ',
      round(roc_res_HSPA5$auc, 3),
      ' (',
      round(roc_res_HSPA5$ci[1], 3),
      '-',
      round(roc_res_HSPA5$ci[3], 3),
      ')'
    ),
    x = 0.7,
    y = 0.1,
    color = '#4DBBD5'
  ) +
  theme_bw() +
  # 固定xy轴比例
  coord_fixed() +
  labs(x = '1-Specificity (FPR)', y = 'Sensitivity (TPR)')

ggsave(file = 'ROC_SCD.pdf', p, height = 6, width = 6)

rm(list = ls())
library(pROC)
library(xgboost)
library(tidyverse)
library(SHAPforxgboost)
setwd("D:\\lqhGroup\\R\\project\\Reproduce")
load("GSE74986.RData")
dat_group <- data.frame(group = group_list)
dat_expr <- exp
head(dat_group)
head(exp)
# 替换元素
group_list <- ifelse(group_list == "HEALTHY", "HC", 
                     ifelse(group_list == "SEVERE_ASTHMA", "SA", group_list))

# 创建数据框
dat_group <- data.frame(group = group_list)

# 将 group 列转换为因子
dat_group$group <- factor(dat_group$group, levels = c('HC', 'SA'))
head(dat_group)
# 定义目标基因和对应颜色
genes <- c("SGK1", "HSPA5", "SCD")
colors <- c("#E64B35", "#4DBBD5", "#00A087") # 红、蓝、绿配色

# 创建空数据框存储结果
roc_data <- data.frame()
auc_data <- data.frame()

# 循环处理每个基因
for (i in seq_along(genes)) {
  gene <- genes[i]
  
  # 提取基因表达数据
  dat_expr_gene <- dat_expr[gene, ] %>% t() %>% as.data.frame()
  colnames(dat_expr_gene) <- gene
  dat_expr_gene$group <- dat_group$group
  
  # ROC分析
  roc_res <- roc(group ~ get(gene), 
                 data = dat_expr_gene,
                 auc = T, ci = T)
  
  # 整理ROC曲线数据
  temp_roc <- data.frame(
    gene = gene,
    specificity = roc_res$specificities,
    sensitivity = roc_res$sensitivities
  )
  
  # 整理AUC数据
  temp_auc <- data.frame(
    gene = gene,
    auc = roc_res$auc,
    ci_lower = roc_res$ci[1],
    ci_upper = roc_res$ci[3],
    color = colors[i]
  )
  
  # 合并数据
  roc_data <- rbind(roc_data, temp_roc)
  auc_data <- rbind(auc_data, temp_auc)
}

# 数据排序（确保曲线绘制正确）
roc_data <- roc_data %>%
  group_by(gene) %>%
  dplyr::arrange(desc(specificity), sensitivity)

# 绘制ROC曲线
p <- ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity, color = gene)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
  scale_color_manual(values = colors) + # 自定义颜色映射
  
  # 添加AUC标注
  geom_text(
    data = auc_data,
    aes(x = 0.6, y = 0.15 - seq(0,0.2,length.out=3), # 纵向排列标注
        label = sprintf("AUC = %.2f (%.2f-%.2f)", auc, ci_lower, ci_upper),
        color = gene),
    hjust = 0, show.legend = FALSE
  ) +
  
  theme_bw() +
  coord_fixed() +
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)") +
  theme(
    legend.position = c(0.7, 0.2),
    legend.title = element_blank()
  )

# 保存图形
ggsave("Combined_ROC2.pdf", p, width = 7, height = 6)

rm(list = ls())
setwd("D:\\lqhGroup\\R\\project\\Reproduce")
library(CIBERSORT)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggcorrplot)
library(psych)
load("GSE74986.RData")
dat_expr <- exp
dat_group <- data.frame(group = group_list)
# 替换元素
group_list <- ifelse(group_list == "HEALTHY", "HC", 
                     ifelse(group_list == "SEVERE_ASTHMA", "SA", group_list))

# 创建数据框并添加样本名称
dat_group <- data.frame(
  group = group_list,
  sample = colnames(dat_expr)  # 新增样本名称列
)

# 将 group 列转换为因子
dat_group$group <- factor(dat_group$group, levels = c('HC', 'SA'))

# 转换sample列为因子（确保顺序一致）
key_gene <- readRDS('key_gene.Rds')


# CIBERSORT免疫浸润分析 ---------------------------------------------------------------

# CIBERSORT分析
dat_res_cibersort <- cibersort(
  # LM22免疫细胞
  LM22,
  as.matrix(dat_expr),
  # 排列次数
  perm = 100,
  # 当输入的表达谱为Count（整数）时选择F，其余情况均用T
  QN = T) %>%
  as.data.frame()

# 取免疫细胞浸润丰度数据
dat_res_cibersort <- dat_res_cibersort[, 1:22]
# 最后三列是p值、相关性和均方根误差，均不取

# 去除丰度均为零的细胞
dat_res_cibersort <- dat_res_cibersort[, colSums(dat_res_cibersort) > 0]


# CIBERSORT免疫浸润叠加柱状图 ---------------------------------------------------------------

# 整理叠加柱状图数据
dat_cibersort_stack <- dat_res_cibersort
dat_cibersort_stack$sample <- rownames(dat_cibersort_stack)
dat_cibersort_stack <- reshape2::melt(dat_cibersort_stack)
colnames(dat_cibersort_stack) <- c('sample', 'cell', 'score')
dat_cibersort_stack$sample <- factor(dat_cibersort_stack$sample,
                                     levels = unique(dat_cibersort_stack$sample))
dat_group$sample <- factor(dat_group$sample, levels = dat_group$sample)

# 叠加柱状图
p <- (
  ggplot(dat_cibersort_stack, aes(sample, score, fill = cell)) +
    geom_col(position = 'fill') +
    theme_classic() +
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(x = element_blank(), y = 'Relative Percentage', fill = 'Immune Cell') +
    scale_fill_manual(
      values = c(
        '#00A087',
        '#3C5488',
        '#F39B7F',
        '#8491B4',
        '#91D1C2',
        '#DC0000',
        '#7E6148',
        '#B09C85',
        '#0072B5',
        '#BC3C29',
        '#E18727',
        '#20854E',
        '#7876B1',
        '#6F99AD',
        '#FFDC91',
        '#EE4C97',
        '#3B4992',
        '#EE0000',
        '#008B45',
        '#631879',
        '#008280',
        '#BB0021'
      )
    )
) %>%
  aplot::insert_bottom(
    ggplot(dat_group, aes(sample, 'Group', fill = group)) +
      geom_tile() +
      theme_classic() +
      theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank()
      ) +
      labs(x = element_blank(), y = element_blank(), fill = 'Group') +
      scale_fill_manual(values = c('#4DBBD5', '#E64B35')),
    height = 0.03
  )

ggsave(file = '免疫/CIBERSORT_stackplot.pdf', p, width = 10, height = 6)


# CIBERSORT免疫浸润分组比较图 ---------------------------------------------------------------

# 分组比较图数据
long_dat_cibersort <- dat_res_cibersort
long_dat_cibersort$group <- dat_group$group
long_dat_cibersort <- reshape2::melt(long_dat_cibersort)
colnames(long_dat_cibersort) = c('group', 'cell', 'score')
long_dat_cibersort$group <- factor(long_dat_cibersort$group, levels = c('HC', 'SA'))

# 分组比较图
p <- ggplot(long_dat_cibersort, aes(cell, score, fill = group)) +
  # 箱线图，不显示离群值
  geom_boxplot(outlier.color = NA) +
  # 小提琴图
  # geom_violin(scale = 'width') +
  # 添加显著性标志
  stat_compare_means(
    aes(group = group),
    # 显著性标志为星号
    label = 'p.signif',
    # 统计方法
    method = 'wilcox.test',
    # 非配对样本
    paired = F,
    # 显著性的cutoff和symbol
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c('***', '**', '*', 'ns')
    )
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = element_blank(), y = 'Immune Infiltration', fill = 'Group') +
  scale_fill_manual(values = c('#4DBBD5', '#E64B35'))

ggsave(file = '免疫/CIBERSORT_boxplot.pdf', p, width = 12, height = 6)
# CIBERSORT免疫浸润相关性棒棒糖图 ---------------------------------------------------------------
dat_group_tumor <- dat_group[dat_group$group == 'SA', ]
dat_expr_tumor <- dat_expr[, dat_group_tumor$sample] %>% as.data.frame()

# 提取疾病样本的免疫细胞浸润丰度
dat_res_cibersort_tumor <- dat_res_cibersort[colnames(dat_expr_tumor), ]

# 计算相关性
# 以SIGLEC1为例
cor_res_cibersort_SIGLEC1 <- psych::corr.test(t(dat_expr_tumor['HSPA5', ]), 
                                              dat_res_cibersort_tumor,
                                              # 统计方法和p值校正方法
                                              method = 'spearman',
                                              adjust = 'BH')

# 整理相关性数据
cor_res_cibersort_SIGLEC1 <- data.frame(
  cell = colnames(dat_res_cibersort_tumor),
  cor = as.vector(cor_res_cibersort_SIGLEC1$r),
  pvalue = as.vector(cor_res_cibersort_SIGLEC1$p),
  row.names = colnames(dat_res_cibersort_tumor)
)
cor_res_cibersort_SIGLEC1 <- cor_res_cibersort_SIGLEC1 %>%
  # 加一列相关性的绝对值
  dplyr::mutate(cor_abs = abs(cor)) %>%
  dplyr::filter(pvalue < 0.05) %>%
  dplyr::arrange(desc(cor))
cor_res_cibersort_SIGLEC1$cell <- factor(cor_res_cibersort_SIGLEC1$cell, levels = cor_res_cibersort_SIGLEC1$cell)

# 相关性棒棒糖图
p <- ggplot(cor_res_cibersort_SIGLEC1,
            aes(cell, cor, color = pvalue, size = cor_abs)) +
  geom_segment(aes(xend = cell, yend = 0),
               color = 'black',
               lwd = 0) +
  geom_point(stat = 'identity') +
  # 水平线
  geom_hline(yintercept = 0) +
  theme_bw() +
  # xy轴颠倒
  coord_flip() +
  labs(x = element_blank(), y = 'Correlation', size = '|Correlation|') +
  scale_color_gradient(
    low = '#E64B35',
    high = '#4DBBD5',
    limits = c(0, 0.05)
  )

ggsave(file = '免疫/CIBERSORT_corr_lollipop_HSPA5.pdf', p, width = 6, height = 4)

library(clusterProfiler)
library(org.Hs.eg.db)  # 如果是人类数据
library(DOSE)
library(enrichplot)
library(tibble)
rm(list=ls())
load("GSE74986.RData")
# 1. 加载你的基因表达数据 (假设是dataframe格式, 行是基因，列是样本)
expr_data <- exp # 需要包含行名，即基因名
expr_data <- as.matrix(expr_data)
# 2. 选择你的目标基因
target_genes <- c("SGK1")

# 3. 计算相关性（Spearman相关性）
target_expr <- expr_data[target_genes, , drop = TRUE]  # 向量，长度=样本数
gene_cor <- apply(expr_data, 1, function(x) {
  cor(x, target_expr, method = "spearman", use = "pairwise.complete.obs")
})
# 4. 按相关性排序
gene_list <- sort(gene_cor, decreasing = TRUE)

# 5. 转换基因名为 ENTREZ ID
gene_df <- bitr(names(gene_list), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
gene_list <- gene_list[gene_df$SYMBOL]
names(gene_list) <- gene_df$ENTREZID  # 用 ENTREZ ID 进行 GSEA


# 6. 运行 GSEA
gsea_result <- gseGO(
  geneList = gene_list, 
  OrgDb = org.Hs.eg.db, 
  ont = "ALL",  # 可选 "MF", "CC"
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# 7. 结果可视化（只保留前 10 个显著结果） 
top_categories <- 10  
# 选择要显示的通路数量# 使用 `head()` 选取 p.adjust 最显著的前10个通路 
gsea_df <- gsea_result 
gsea_df@result <- gsea_result@result %>% arrange(p.adjust) %>% head(top_categories)
# 7. 结果可视化
dotplot(gsea_df, showCategory = 10)  # 气泡图
gseaplot2(gsea_df, geneSetID = 1, title = gsea_result$Description[1])  # GSEA 曲线


# ==============================================================================
# 0. 环境设置与结果目录初始化
# ==============================================================================
rm(list = ls())
setwd("~/ZYJ/BPA-TC")  # 请确保 .h5 文件在此目录下

result_dir <- "Analysis_Results_Asthma_BPA"
if(!dir.exists(result_dir)) {
  dir.create(result_dir)
  message(paste("已创建结果文件夹:", result_dir))
} else {
  message(paste("结果将保存至现有文件夹:", result_dir))
}

library(tidyverse)
library(Seurat)
library(harmony)
library(SingleR)
library(AUCell)
library(scTenifoldKnk)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(enrichplot)

# ==============================================================================
# 1. 读取 .h5 文件并构建 Seurat 对象 (GSE156285)
# ==============================================================================
message("正在读取 GSE156285 的 .h5 文件...")

# 获取当前目录下所有 .h5 文件（可根据实际文件名调整）
h5_files <- list.files(pattern = "*matrix\\.h5$")
if(length(h5_files) == 0) {
  stop("未找到 .h5 文件，请检查工作目录。")
}
message(paste("找到以下文件:", paste(h5_files, collapse = ", ")))

# 逐个读取并创建 Seurat 对象，同时添加分组信息
seu_list <- list()
for (f in h5_files) {
  # 从文件名提取样本标识（例如 "GSM4728360_IT_7741246983_matrix.h5" -> "IT_7741246983"）
  sample_id <- gsub("_matrix\\.h5$", "", f)
  
  # 读取 10X h5 文件
  counts <- Read10X_h5(f)
  
  # 创建 Seurat 对象，过滤低质量细胞（可选，后续还会统一过滤）
  seu <- CreateSeuratObject(counts = counts, project = sample_id, 
                            min.cells = 3, min.features = 200)
  
  # 根据文件名判断分组：IT（下鼻甲）或 P（息肉）
  if (grepl("_IT_", f) || grepl("_IT_matrix", f)) {
    seu$group <- "IT"
  } else if (grepl("_P_", f) || grepl("_P_matrix", f)) {
    seu$group <- "P"
  } else {
    seu$group <- "Unknown"  # 如果文件名不含 IT 或 P，可手动添加
  }
  
  # 存储到列表
  seu_list[[sample_id]] <- seu
}

# 合并所有样本（添加细胞前缀避免混淆）
sc_obj <- merge(seu_list[[1]], y = seu_list[-1], 
                add.cell.ids = names(seu_list))

message(paste("数据合并完成！总细胞数:", ncol(sc_obj)))

# 可选：如果基因名是 Ensembl ID，转换为 Symbol（检查前几个基因名）
if (grepl("^ENSG", rownames(sc_obj)[1])) {
  message("检测到基因名为 Ensembl ID，正在转换为 Symbol...")
  ensembl_ids <- rownames(sc_obj)
  symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, 
                    keytype = "ENSEMBL", column = "SYMBOL")
  # 处理 NA 和重复
  new_rownames <- ifelse(is.na(symbols), ensembl_ids, symbols)
  # 去重（保留表达量最高的基因）
  if (any(duplicated(new_rownames))) {
    # 简单合并：取每个基因的最大 counts（更严谨的做法可合并相同基因）
    sc_obj <- sc_obj[!duplicated(new_rownames), ]
    rownames(sc_obj) <- new_rownames[!duplicated(new_rownames)]
  } else {
    rownames(sc_obj) <- new_rownames
  }
  message("基因名转换完成。")
}

# ==============================================================================
# 2. 基础预处理与降维
# ==============================================================================
message("开始标准预处理流程...")
# 合并所有层为一个统一的层（解决多 layer 问题）
if (length(sc_obj@assays$RNA@layers) > 1) {
  sc_obj <- JoinLayers(sc_obj)
  message("已合并所有层为一个统一的层。")
}

# 计算线粒体基因比例
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
# ... 后续代码不变

# 计算线粒体基因比例
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")

# 可视化质控指标（可选，先查看分布再过滤）
# VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 质控过滤（阈值可根据小提琴图调整）
sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

# 标准化、找高变基因、缩放
sc_obj <- NormalizeData(sc_obj)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000)
sc_obj <- ScaleData(sc_obj)

# PCA
sc_obj <- RunPCA(sc_obj)
ElbowPlot(sc_obj)  # 查看主成分方差贡献，选择合适的主成分数（例如 1:15）

# UMAP 和聚类（使用前 15 个 PC）
sc_obj <- RunUMAP(sc_obj, dims = 1:15)
sc_obj <- FindNeighbors(sc_obj, dims = 1:15)
sc_obj <- FindClusters(sc_obj, resolution = 0.6)

# 保存 UMAP 聚类图
p_umap_cluster <- DimPlot(sc_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("Clusters (GSE156285)")
ggsave(file.path(result_dir, "01_UMAP_Clusters.pdf"), plot = p_umap_cluster, width = 6, height = 5)

# ==============================================================================
# 3. 细胞注释 (SingleR) & 自动寻找主要细胞
# ==============================================================================
message("正在加载参考数据集进行 SingleR 注释...")
ref_path <- "~/ZYJ/HRT/单细胞/HumanPrimaryCellAtlas_hpca.se_human.RData"
if(!file.exists(ref_path)) stop(paste("错误：找不到参考文件:", ref_path))

loaded_ref <- load(ref_path) 
ref <- get(loaded_ref[1]) 

testdata <- LayerData(sc_obj, assay = "RNA", layer = "data")
pred <- SingleR(test = testdata, ref = ref, labels = ref$label.main, 
                clusters = sc_obj$seurat_clusters)

celltype_df <- data.frame(ClusterID = rownames(pred), CellType = pred$labels)
sc_obj$celltype <- celltype_df$CellType[match(sc_obj$seurat_clusters, celltype_df$ClusterID)]

p_anno <- DimPlot(sc_obj, group.by = "celltype", label = TRUE, repel = TRUE) + 
  ggtitle("SingleR Annotation (GSE156285)")
ggsave(file.path(result_dir, "02_CellType_Annotation.pdf"), plot = p_anno, width = 8, height = 6)

# 自动寻找占比最高的主要细胞
cell_counts <- table(sc_obj$celltype)
main_celltype <- names(sort(cell_counts, decreasing = TRUE))[1]
message(paste("\n>>> 自动锁定主要细胞类型为:", main_celltype, "(细胞数:", max(cell_counts), ") <<<\n"))

# ==============================================================================
# 4. 关键靶点基因分析 (HSPA5 & HSPA5) - 纯 ggplot2 绕过版本报错
# ==============================================================================
message(">>> 分析哮喘关键靶标基因表达 (使用纯 ggplot2 绘图)...")

target_genes <- c("HSPA5", "HSPA5")
valid_genes <- target_genes[target_genes %in% rownames(sc_obj)]

if(length(valid_genes) > 0) {
  
  # 1. 安全提取表达矩阵 (使用 V5 规范的 layer="data")
  # 这一步绝对不会触发 slot 报错
  expr_mat <- GetAssayData(sc_obj, assay = "RNA", layer = "data")[valid_genes, , drop = FALSE]
  expr_df <- as.data.frame(t(as.matrix(expr_mat)))
  expr_df$celltype <- sc_obj$celltype
  
  # 2. 纯手工绘制 VlnPlot
  library(tidyr)
  df_long <- pivot_longer(expr_df, cols = all_of(valid_genes), names_to = "Gene", values_to = "Expression")
  
  p_vln <- ggplot(df_long, aes(x = celltype, y = Expression, fill = celltype)) +
    geom_violin(scale = "width", trim = TRUE) +
    facet_wrap(~ Gene, ncol = 1, scales = "free_y") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "none",
      strip.text = element_text(size = 12, face = "bold")
    ) +
    labs(x = "", y = "Expression Level")
  
  ggsave(file.path(result_dir, "03_TargetGenes_VlnPlot.pdf"), plot = p_vln, width = 8, height = 5)
  
  # 3. 纯手工绘制 FeaturePlot (UMAP 投射图)
  umap_df <- as.data.frame(Embeddings(sc_obj, "umap"))
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  plot_umap_df <- cbind(umap_df, expr_df)
  
  p_feat_list <- list()
  for(g in valid_genes) {
    p_feat_list[[g]] <- ggplot(plot_umap_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[g]])) +
      geom_point(size = 0.5, alpha = 0.8) +
      scale_color_gradient(low = "lightgrey", high = "#CC3333") +
      theme_classic() +
      ggtitle(g) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  }
  
  # 使用 ggpubr 拼图并保存
  p_feat <- ggarrange(plotlist = p_feat_list, ncol = length(valid_genes))
  ggsave(file.path(result_dir, "04_TargetGenes_FeaturePlot.pdf"), plot = p_feat, width = 10, height = 5)
  
  message(">>> 靶点基因表达分布图 (VlnPlot & FeaturePlot) 绘制并保存成功！")
  
} else {
  message("警告：在矩阵中未找到 HSPA5 或 HSPA5。请检查前几个基因名：", head(rownames(sc_obj), 5))
}



# ==============================================================================
# 5. 专项分析：Epithelial Cell 虚拟敲除 (极速版 - HSPA5)
# ==============================================================================
library(Seurat)
library(Matrix)

# 1. 提取上皮细胞
message(">>> 正在提取并准备 Epithelial_cells 数据...")
epi_sub <- subset(sc_obj, celltype == "Epithelial_cells")
epi_sub <- JoinLayers(epi_sub)

# ==========================================
# 【新增：核心提速模块 - 极速降维】
# ==========================================
# A. 细胞降采样 (最多保留 800 个细胞)
max_cells <- 800
if (ncol(epi_sub) > max_cells) {
  set.seed(42) # 保证每次抽取的细胞一样，结果可重复
  epi_sub <- subset(epi_sub, cells = sample(Cells(epi_sub), max_cells))
  message(paste(">>> 细胞数较多，已随机降采样至", max_cells, "个细胞..."))
}

# B. 基因过滤 (只保留 Top 2000 高变基因)
epi_sub <- FindVariableFeatures(epi_sub, selection.method = "vst", nfeatures = 2000)
hvgs <- VariableFeatures(epi_sub)

# ⚠️ 极其重要：强制将你的目标基因加入白名单，防止被当做非高变基因过滤掉！
target_genes <- c("SGK1", "HSPA5")
genes_to_keep <- unique(c(hvgs, target_genes))
# 确保这些基因确实在当前矩阵中
genes_to_keep <- genes_to_keep[genes_to_keep %in% rownames(epi_sub)] 

# 2. 提取过滤后的稀疏矩阵
raw_counts <- GetAssayData(epi_sub, assay = "RNA", layer = "counts")[genes_to_keep, ]
epi_mat <- as(raw_counts, "dgCMatrix")

message(">>> 过滤后矩阵大小: ", nrow(epi_mat), " 基因 x ", ncol(epi_mat), " 细胞")
# ==========================================

# 3. 运行 HSPA5 敲除
gKO <- "HSPA5"
message(">>> [任务启动] 正在上皮细胞中虚拟敲除: ", gKO, " (预计只需 2-5 分钟...)")

# 现在的 epi_mat 只有约 2000 个基因，速度会起飞
result_HSPA5 <- scTenifoldKnk(countMatrix = epi_mat, gKO = gKO, qc_minLSize = 10) # 适当调低 minLSize 防止严格质控报错
save(result_HSPA5, file = file.path(result_dir, paste0("Epi_KO_", gKO, "_Result.RData")))

# ==============================
# 4. 结果可视化 (HSPA5)
# ==============================
message(">>> 正在生成可视化图表...")
df_HSPA5 <- result_HSPA5$diffRegulation
df_HSPA5$log_pval <- -log10(df_HSPA5$p.adj)

# --- Top 20 基因条形图 ---
top_genes <- head(df_HSPA5[order(-df_HSPA5$FC), ], 20)

pdf(file = file.path(result_dir, paste0('Epi_KO_', gKO, '_top20.pdf')), width = 5, height = 4)
p1 <- ggplot(top_genes, aes(x = reorder(gene, FC), y = FC, fill = FC)) +
  geom_bar(stat = "identity") +
  scale_y_log10() +   
  scale_fill_gradient(low = "steelblue", high = "red") +
  labs(x = "", y = "log10(FC)", title = paste("Top 20 Regulated Genes -", gKO, "KO")) +
  theme_classic() +
  coord_flip()
print(p1)
dev.off()

# --- Z-score 火山图 ---
label_genes <- subset(df_HSPA5, abs(Z) > 6.5 & p.adj < 0.05)
top10_genes <- head(df_HSPA5[order(-df_HSPA5$FC), ], 10)

pdf(file = file.path(result_dir, paste0('Epi_KO_', gKO, '_zscore.pdf')), width = 5, height = 5)
p2 <- ggplot(df_HSPA5, aes(x = Z, y = log_pval)) +
  geom_point(alpha = 0.5, color = "black") +            
  geom_point(data = top10_genes, aes(x = Z, y = log_pval), color = "red", size = 2) +                 
  geom_vline(xintercept = c(1), linetype = "dashed", color = 'steelblue') +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = 'steelblue') +
  geom_text_repel(data = top10_genes, aes(label = gene), size = 3, max.overlaps = 50) +
  labs(title = paste("Differentially Regulated Genes (", gKO, "KO vs Ctrl)"),
       x = "Z-score", y = "-log10(P.adjust)") +
  theme_classic()
print(p2)
dev.off()

# ==============================
# ==============================================================================
# 富集分析与绘图 (修复版)
# ==============================================================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ggplot2)
library(enrichplot)

# 0. 准备数据 (假设 df_HSPA5 是之前的敲除结果)
# 如果你是单独运行这段，请确保 df_HSPA5 存在
df <- result_HSPA5$diffRegulation 
siggene <- subset(df, p.adj < 0.05)

# 如果差异基因太少，退而求其次选择 Z-score 较大的前 100 个
if(nrow(siggene) < 5){
  message(">>> 显著差异基因过少，正在提取 Z-score 前 100 的基因进行分析...")
  siggene <- head(df[order(-abs(df$Z)), ], 100)
}

# 1. 基因名转换
dat_genelist <- bitr(siggene$gene, 
                     fromType = 'SYMBOL', 
                     toType = 'ENTREZID', 
                     OrgDb = 'org.Hs.eg.db')

# 2. GO 分析 (放宽阈值以确保在小数据集下能运行)
res_go <- enrichGO(
  gene = dat_genelist$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.5, # 放宽以便筛选
  qvalueCutoff = 0.5,
  readable = TRUE
)
dat_res_go <- as.data.frame(res_go)

# 3. KEGG 分析
res_kegg <- enrichKEGG(
  gene = dat_genelist$ENTREZID,
  organism = 'hsa',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.5,
  qvalueCutoff = 0.5
) %>% setReadable(OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
dat_res_kegg <- as.data.frame(res_kegg)

# ==============================================================================
# 4. 关键词筛选与数据合并 (哮喘主题)
# ==============================================================================
# 定义哮喘相关关键词
keywords <- "asthma|interleukin|IL-|inflammation|bronchial|airway|eosinophil|mucus|epithelial"

# 提取显著结果
go_sig <- dat_res_go[dat_res_go$pvalue < 0.5, ]
kegg_sig <- dat_res_kegg[dat_res_kegg$pvalue < 0.5, ]

# 筛选关键词通路
go_sub <- go_sig[grepl(keywords, go_sig$Description, ignore.case = TRUE), ]
kegg_sub <- kegg_sig[grepl(keywords, kegg_sig$Description, ignore.case = TRUE), ]

# 扩充至每个类别 8 条通路 (保底机制)
target_n <- 8

process_enrich_data <- function(sub_df, sig_df, ontology_name) {
  if (nrow(sub_df) < target_n) {
    candidates <- sig_df[!sig_df$Description %in% sub_df$Description, ]
    # 排除 KEGG 中的癌症/疾病干扰 (非哮喘类)
    if(ontology_name == "KEGG"){
      candidates <- candidates[!grepl("cancer|disease|syndrome", candidates$Description, ignore.case = TRUE), ]
    }
    candidates <- candidates[order(candidates$p.adjust), ]
    needed <- target_n - nrow(sub_df)
    combined <- rbind(sub_df, head(candidates, needed))
  } else {
    combined <- head(sub_df, target_n)
  }
  if(nrow(combined) > 0) combined$ONTOLOGY <- ontology_name
  return(combined)
}

go_combined <- process_enrich_data(go_sub, go_sig, "GO")
kegg_combined <- process_enrich_data(kegg_sub, kegg_sig, "KEGG")

# 合并最终绘图数据
df3 <- rbind(go_combined, kegg_combined)

# ==============================================================================
# 5. 绘图 (PDF 生成)
# ==============================================================================
if (nrow(df3) > 0) {
  # 排序优化
  df3 <- df3 %>% 
    arrange(ONTOLOGY, desc(p.adjust)) %>%
    mutate(Description = factor(Description, levels = unique(Description)))
  
  # 计算分类侧栏位置
  cate_bar <- df3 %>% 
    group_by(ONTOLOGY) %>% 
    summarise(count = n(), .groups = 'drop') %>%
    mutate(
      ymax = cumsum(count) + 0.4,
      ymin = ymax - count,
      mid = (ymax + ymin) / 2
    )
  
  two_colors <- c("#E6E6AA", "#6A6AaF")
  annotate_bar_range <- c(-0.4, -0.2)
  
  pdf(file = paste0('Epi_KO_', gKO, '_Asthma_Enrich.pdf'), width = 8, height = 8)
  p <- ggplot(df3, aes(-log10(p.adjust), fill = ONTOLOGY)) +
    geom_col(aes(y = Description), width = 0.6, show.legend = FALSE) +
    geom_text(aes(x = 0.05, y = Description, label = Description), 
              hjust = 0, size = 3.5) +
    geom_text(aes(x = 0.05, y = Description, label = geneID, colour = ONTOLOGY),
              hjust = 0, vjust = 2.2, size = 2.8, fontface = 'italic', show.legend = FALSE) +
    # 左侧分类色块
    geom_rect(data = cate_bar, 
              mapping = aes(xmin = annotate_bar_range[1], xmax = annotate_bar_range[2],
                            ymin = ymin, ymax = ymax, fill = ONTOLOGY),
              inherit.aes = FALSE, show.legend = FALSE) +
    # 分类标签
    geom_text(data = cate_bar, 
              mapping = aes(x = sum(annotate_bar_range)/2, y = mid, label = ONTOLOGY),
              angle = 90, size = 4.5, color = "black") +
    # 坐标轴美化
    geom_segment(aes(x = 0, xend = max(-log10(p.adjust)) * 1.1, y = 0, yend = 0), color = "black") +
    geom_segment(aes(x = 0, xend = 0, y = 0, yend = nrow(df3) + 1), color = "black") +
    labs(x = "-log10(P.adjust)", y = NULL, title = paste(gKO, "Knockout - Asthma Related Pathways")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = two_colors) +
    scale_colour_manual(values = two_colors) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  print(p)
  dev.off()
  
  message("✅ PDF 已成功生成：Epi_KO_", gKO, "_Asthma_Enrich.pdf")
} else {
  message("⚠️ 警告：未找到符合条件的通路，未生成 PDF。")
}

save(dat_res_kegg, dat_res_go, file = paste0("Epi_KO_", gKO, "_EnrichData.RData"))


