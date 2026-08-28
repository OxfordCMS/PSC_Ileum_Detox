
# Evidence for bile acid mediated ileal barrier function deficits in primary sclerosing cholangitis (PSC_Ileum_Detox)

## Overview

This repository contains functions and rmarkdown files which contain code used to process and analyse the data that are published at [link when available]. The structure of the repository is as follows:

```
├── R
│   └── scfunctions.R
├── README.md
└── Rmd
    ├── bulk_rnaseq
    │   └── differential_expression.Rmd
    ├── proteomics
    │   ├── differential_proteomics.Rmd
    │   └── QC.Rmd
    └── scrnaseq
        ├── Caecum_Clustering.Rmd
        ├── Caecum_Epithelial_Cells.Rmd
        ├── Caecum_Immune_Cells.Rmd
        ├── Epithelium_Compare_Tissue_Differences.Rmd
        ├── Ileum_Clustering.Rmd
        ├── Ileum_Epithelial_Cells.Rmd
        ├── Ileum_Immune_Cells.Rmd
        ├── Ileum_Spatial.Rmd
        ├── Ileum_Stroma.Rmd
        ├── Ileum_T_Cells_Compass.Rmd
        ├── Immune_Cells_Compare_Tissue_Differences.Rmd
        ├── Mouse_Ileum.Rmd
        ├── QC.Rmd
        └── Spatial_Ileum_Epithelium.Rmd
```

The analyses are based on pre-processed scRNA-seq, LC-MS bulk proteomics and bulk RNA-seq - pre-processing was performed using commandline tools as described in the accompanying manuscript. Processed data used for the analyses and compatible with the code in this repo can be found at figshare in the collection [add link when ready].
