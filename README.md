# Rdata-BPA
.
├── Data Acquisition & Preprocessing
│   ├── GSE74986 data download and filtering
│   ├── Probe ID conversion (GPL6480 platform)
│   └── Quality control (boxplot, PCA)
├── WGCNA Analysis
│   ├── Soft-thresholding power selection
│   ├── Module identification and merging
│   └── Module-trait association analysis
├── Differential Expression Analysis
│   ├── limma differential analysis
│   ├── Volcano plot and heatmap
│   └── GO/KEGG enrichment analysis
├── Machine Learning Screening
│   ├── LASSO regression
│   ├── Random Forest (RF)
│   └── XGBoost
├── Diagnostic Model Construction
│   ├── Logistic regression nomogram
│   ├── Calibration curve
│   └── Decision Curve Analysis (DCA)
├── Immune Infiltration Analysis
│   └── CIBERSORT
├── Single-cell Validation (GSE156285)
│   ├── Seurat standard workflow
│   ├── SingleR cell annotation
│   └── Virtual knockout analysis (scTenifoldKnk)
└── Output Files
    ├── RData/Rds intermediate data
    ├── PDF/PNG figures
    └── CSV enrichment results
