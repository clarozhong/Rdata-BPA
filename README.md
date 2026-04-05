# Rdata-BPA
.
├── 数据获取与预处理
│   ├── GSE74986数据下载与筛选
│   ├── 探针ID转换（GPL6480平台）
│   └── 数据质量控制（箱线图、PCA）
├── WGCNA分析
│   ├── 软阈值筛选
│   ├── 模块识别与合并
│   └── 模块-性状关联分析
├── 差异表达分析
│   ├── limma差异分析
│   ├── 火山图与热图
│   └── GO/KEGG富集分析
├── 机器学习筛选
│   ├── LASSO回归
│   ├── 随机森林（RF）
│   └── XGBoost
├── 诊断模型构建
│   ├── Logistic回归列线图
│   ├── 校准曲线
│   └── 决策曲线分析（DCA）
├── 免疫浸润分析
│   └── CIBERSORT
├── 单细胞验证（GSE156285）
│   ├── Seurat标准流程
│   ├── SingleR细胞注释
│   └── 虚拟敲除分析（scTenifoldKnk）
└── 输出文件
    ├── RData/Rds格式中间数据
    ├── PDF/PNG格式图片
    └── CSV格式富集分析结果
    # 安装必要包
install.packages(c(
  "tidyverse", "ggplot2", "ggvenn", "pheatmap", 
  "reshape2", "ggpubr", "ggcorrplot", "psych", "rms",
  "forestplot", "pROC", "xgboost", "SHAPforxgboost"
))

# BiocManager安装包
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery", "limma", "WGCNA", "clusterProfiler", 
  "org.Hs.eg.db", "enrichplot", "DOSE", "FactoMineR",
  "factoextra", "SingleR", "Seurat", "harmony", "AUCell",
  "scTenifoldKnk", "ggDCA", "GSVA", "CIBERSORT"
))
