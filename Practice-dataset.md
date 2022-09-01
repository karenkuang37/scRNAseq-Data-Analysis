Practice dataset (1K)
================
Karen Kuang
2022-09-01

## Dataset overview

Moderate-size 1k peripheral blood mononuclear cells (PBMCs) from a
Healthy Donor (v3 chemistry). PBMCs are primary cells with relatively
small amounts of RNA.

Data retrieved from
<https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0>

### 1. Load Seurat object

``` r
# Load data using .RDS format

# Load data using .mtx file

# Load data using .loom files

# Load data using 10X CellRanger .HDF5 format
pbmc1k.data <- Read10X_h5("pbmc_1k_v3_raw_feature_bc_matrix.h5",
                          use.names = TRUE,
                          unique.features = TRUE)
str(pbmc1k.data)
```

    ## Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##   ..@ i       : int [1:3394796] 1814 4026 6625 17059 18002 20999 30468 33502 33498 7966 ...
    ##   ..@ p       : int [1:6794881] 0 0 0 0 0 0 0 0 0 8 ...
    ##   ..@ Dim     : int [1:2] 33538 6794880
    ##   ..@ Dimnames:List of 2
    ##   .. ..$ : chr [1:33538] "MIR1302-2HG" "FAM138A" "OR4F5" "AL627309.1" ...
    ##   .. ..$ : chr [1:6794880] "AAACCCAAGAAACACT-1" "AAACCCAAGAAACCAT-1" "AAACCCAAGAAACCCA-1" "AAACCCAAGAAACCCG-1" ...
    ##   ..@ x       : num [1:3394796] 1 1 1 1 1 1 1 1 1 1 ...
    ##   ..@ factors : list()

``` r
# Initialize the Seurat object with raw data
pbmc1k <- CreateSeuratObject(counts = pbmc1k.data, project = "PBMC1k", min.cells = 3, min.features = 200)
str(pbmc1k)
```

    ## Formal class 'Seurat' [package "SeuratObject"] with 13 slots
    ##   ..@ assays      :List of 1
    ##   .. ..$ RNA:Formal class 'Assay' [package "SeuratObject"] with 8 slots
    ##   .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##   .. .. .. .. .. ..@ i       : int [1:2504752] 21 26 28 29 30 39 45 48 49 54 ...
    ##   .. .. .. .. .. ..@ p       : int [1:1203] 0 2618 4423 5982 7207 9035 11081 12668 16088 19835 ...
    ##   .. .. .. .. .. ..@ Dim     : int [1:2] 15254 1202
    ##   .. .. .. .. .. ..@ Dimnames:List of 2
    ##   .. .. .. .. .. .. ..$ : chr [1:15254] "AL627309.1" "AL669831.5" "LINC00115" "FAM41C" ...
    ##   .. .. .. .. .. .. ..$ : chr [1:1202] "AAACCCAAGGAGAGTA-1" "AAACGCTTCAGCCCAG-1" "AAAGAACAGACGACTG-1" "AAAGAACCAATGGCAG-1" ...
    ##   .. .. .. .. .. ..@ x       : num [1:2504752] 1 1 1 1 2 1 1 2 2 1 ...
    ##   .. .. .. .. .. ..@ factors : list()
    ##   .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    ##   .. .. .. .. .. ..@ i       : int [1:2504752] 21 26 28 29 30 39 45 48 49 54 ...
    ##   .. .. .. .. .. ..@ p       : int [1:1203] 0 2618 4423 5982 7207 9035 11081 12668 16088 19835 ...
    ##   .. .. .. .. .. ..@ Dim     : int [1:2] 15254 1202
    ##   .. .. .. .. .. ..@ Dimnames:List of 2
    ##   .. .. .. .. .. .. ..$ : chr [1:15254] "AL627309.1" "AL669831.5" "LINC00115" "FAM41C" ...
    ##   .. .. .. .. .. .. ..$ : chr [1:1202] "AAACCCAAGGAGAGTA-1" "AAACGCTTCAGCCCAG-1" "AAAGAACAGACGACTG-1" "AAAGAACCAATGGCAG-1" ...
    ##   .. .. .. .. .. ..@ x       : num [1:2504752] 1 1 1 1 2 1 1 2 2 1 ...
    ##   .. .. .. .. .. ..@ factors : list()
    ##   .. .. .. ..@ scale.data   : num[0 , 0 ] 
    ##   .. .. .. ..@ key          : chr "rna_"
    ##   .. .. .. ..@ assay.orig   : NULL
    ##   .. .. .. ..@ var.features : logi(0) 
    ##   .. .. .. ..@ meta.features:'data.frame':   15254 obs. of  0 variables
    ##   .. .. .. ..@ misc         : list()
    ##   ..@ meta.data   :'data.frame': 1202 obs. of  3 variables:
    ##   .. ..$ orig.ident  : Factor w/ 1 level "PBMC1k": 1 1 1 1 1 1 1 1 1 1 ...
    ##   .. ..$ nCount_RNA  : num [1:1202] 8286 5509 4280 2754 6589 ...
    ##   .. ..$ nFeature_RNA: int [1:1202] 2618 1805 1559 1225 1828 2046 1587 3420 3747 1753 ...
    ##   ..@ active.assay: chr "RNA"
    ##   ..@ active.ident: Factor w/ 1 level "PBMC1k": 1 1 1 1 1 1 1 1 1 1 ...
    ##   .. ..- attr(*, "names")= chr [1:1202] "AAACCCAAGGAGAGTA-1" "AAACGCTTCAGCCCAG-1" "AAAGAACAGACGACTG-1" "AAAGAACCAATGGCAG-1" ...
    ##   ..@ graphs      : list()
    ##   ..@ neighbors   : list()
    ##   ..@ reductions  : list()
    ##   ..@ images      : list()
    ##   ..@ project.name: chr "PBMC1k"
    ##   ..@ misc        : list()
    ##   ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
    ##   .. ..$ : int [1:3] 4 1 1
    ##   ..@ commands    : list()
    ##   ..@ tools       : list()

``` r
pbmc1k
```

    ## An object of class Seurat 
    ## 15254 features across 1202 samples within 1 assay 
    ## Active assay: RNA (15254 features, 0 variable features)

``` r
# 15254 features across 1202 samples within 1 assay
```

### 2. QC

``` r
### % MT reads --- denotes the set of mitocondrial genes (low-quality / dying cells exhibiting extensive mt contamination)
pbmc1k[["percent.mt"]] <- PercentageFeatureSet(pbmc1k, pattern = "^MT-")
View(pbmc1k@meta.data)
```

    ## Warning in system2("/usr/bin/otool", c("-L", shQuote(DSO)), stdout = TRUE):
    ## running command ''/usr/bin/otool' -L '/Library/Frameworks/R.framework/Resources/
    ## modules/R_de.so'' had status 1

``` r
# Visualize QC metrics as a violin plot
VlnPlot(pbmc1k, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

![](Practice-dataset_files/figure-gfm/QC-1.png)<!-- -->

``` r
# Visualize feature-feature relationships with scatter plot
p1 <- FeatureScatter(pbmc1k,feature1 = "nCount_RNA", feature2 = "percent.mt")
# as shown, majority of transcripts have decently low % MT
p2 <- FeatureScatter(pbmc1k,feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
# x-axis = #transcripts, y-axis = #genes
# as shown, most cells follow the line which is nice. (good cells should have a good number of genes AND molecules detected)
p1+p2
```

    ## `geom_smooth()` using formula 'y ~ x'

![](Practice-dataset_files/figure-gfm/QC-2.png)<!-- -->

``` r
### ---------- FILTER -----------
pbmc1k <- subset(pbmc1k, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
# mt thresh selected according to plot, however as suggested in lit, a standardized 10% mtDNA should be applied to scRNA-seq QC of human samples?
# *Find tutorial on DubletFinder, could be useful

### ---------- NORMALIZATION -----------
pbmc1k <- NormalizeData(pbmc1k)
# pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
# Data is log-transformed, with a scale factor (10,000 by default)FEA

### ---------- FEATURE SELECTION ----------
pbmc1k <-FindVariableFeatures(pbmc1k, selection.method = "vst", nfeatures = 2000)
# number of features set to 2000 by default
# ID the top 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc1k),10)
#  [1] "GNLY"   "IGLC2"  "IGLC3"  "FCGR3A" "PPBP"   "CDKN1C" "GZMB"   "ITM2C"  "IGKC"   "S100A9"
p3 <- VariableFeaturePlot(pbmc1k)
p4 <- LabelPoints(plot = p3, points = top10, repel = TRUE)
```

    ## When using repel, set xnudge and ynudge to 0 for optimal results

``` r
p4
```

![](Practice-dataset_files/figure-gfm/QC-3.png)<!-- -->

``` r
# red dots refer to the 2000 highly variable features we indicated

### ---------- SCALING to reduce batch effect ---------
all.genes <- rownames(pbmc1k)
pbmc1k <- ScaleData(pbmc1k, features = all.genes)
```

    ## Centering and scaling data matrix

``` r
# this is where linear transformation is applied, to give equal weight so that highly-expressed genes don't dominate in downstream analysis

### ---------- LINEAR DIMENTIONALITY REDUCTION ------------
pbmc1k <- RunPCA(pbmc1k, features = VariableFeatures(object = pbmc1k))
```

    ## PC_ 1 
    ## Positive:  FCN1, LYZ, S100A9, FGL2, MNDA, CTSS, S100A8, SERPINA1, CST3, PSAP 
    ##     NCF2, VCAN, AIF1, CSTA, TYMP, KLF4, MPEG1, GRN, CPVL, CLEC7A 
    ##     CD68, LST1, MS4A6A, CD14, SLC7A7, LGALS1, CSF3R, CYBB, S100A12, CD36 
    ## Negative:  LTB, TRAC, CD3E, TRBC2, LCK, CD3D, IL32, IL7R, TCF7, CD3G 
    ##     CD69, ISG20, CD27, CD247, SPOCK2, ARL4C, CD7, CD2, LAT, BCL11B 
    ##     GZMM, TRBC1, CD6, NOSIP, RORA, RARRES3, SYNE2, ZAP70, CTSW, CCR7 
    ## PC_ 2 
    ## Positive:  CD3E, IL32, CD247, GZMM, LCK, CTSW, CD7, CD3D, GZMA, ANXA1 
    ##     NKG7, S100A4, CD3G, CST7, TRAC, PRF1, RARRES3, IL7R, LAT, KLRB1 
    ##     ARL4C, RORA, ZAP70, CD2, CCL5, TRBC1, KLRG1, SAMD3, SYNE2, ITGB2 
    ## Negative:  CD79A, MS4A1, IGHM, BANK1, HLA-DQA1, CD79B, LINC00926, IGHD, HLA-DQB1, CD22 
    ##     CD74, HLA-DRB1, HLA-DRA, TNFRSF13C, SPIB, HLA-DPA1, HLA-DQA2, HLA-DPB1, TCL1A, BCL11A 
    ##     IGKC, FCER2, VPREB3, RALGPS2, FAM129C, HLA-DRB5, MEF2C, FCRL1, HVCN1, FCRLA 
    ## PC_ 3 
    ## Positive:  CAVIN2, TUBB1, PF4, GNG11, PPBP, GP9, CD9, CLU, CMTM5, TREML1 
    ##     SPARC, HIST1H2AC, TMEM40, ACRBP, PTCRA, NRGN, ANKRD9, PRKAR2B, CTTN, MTURN 
    ##     ITGA2B, GFI1B, AC147651.1, GMPR, CLEC1B, MAP3K7CL, CLDN5, SMOX, C2orf88, TSC22D1 
    ## Negative:  MT-CO1, CYBA, FOSB, VIM, FOS, ITGB2, NEAT1, LSP1, HNRNPU, CALR 
    ##     LCP1, DUSP1, S100A10, PRELID1, CD74, LTB, KLF6, S100A6, SPCS1, PLAC8 
    ##     IFITM2, S100A4, SEC61B, ISG20, ZFP36L1, PPIB, PEBP1, LCK, ANXA1, APOBEC3G 
    ## PC_ 4 
    ## Positive:  GZMB, CLIC3, GNLY, NKG7, KLRD1, KLRF1, PRF1, SPON2, CST7, FGFBP2 
    ##     GZMA, ADGRG1, HOPX, CCL4, TRDC, TTC38, MATK, IL2RB, TBX21, APOBEC3G 
    ##     RHOC, CTSW, C12orf75, APMAP, FCGR3A, CD160, S1PR5, SH2D1B, CMC1, PTGDR 
    ## Negative:  LEF1, TRABD2A, MAL, TCF7, RCAN3, IL7R, CCR7, NOSIP, CD3D, LTB 
    ##     TRAC, CD3G, CAMK4, NELL2, PASK, CD27, BCL11B, CD5, SLC2A3, TSHZ2 
    ##     FHIT, EGR1, RGS10, RGCC, VIM, LINC01550, CD40LG, ADTRP, CD6, INPP4B 
    ## PC_ 5 
    ## Positive:  GNLY, KLRD1, FGFBP2, KLRF1, ADGRG1, PRF1, CST7, NKG7, CCL4, TRDC 
    ##     HOPX, GZMA, MS4A1, LINC00926, CD79B, SPON2, IGHD, TBX21, CD79A, MATK 
    ##     TTC38, IL2RB, CD160, CD22, MYOM2, S1PR5, SH2D1B, FCER2, TNFRSF13C, GZMH 
    ## Negative:  SERPINF1, LILRA4, SMPD3, CLEC4C, SCT, LRRC26, IL3RA, AL096865.1, PACSIN1, TPM2 
    ##     DNASE1L3, NOTCH4, TNFRSF21, GAS6, DERL3, CIB2, NRP1, RHEX, CUX2, MYBL2 
    ##     ITM2C, SCN9A, UGCG, PPP1R14B, PLD4, LINC00996, PPM1J, NAT8L, ZFAT, SERPINF2

``` r
print(pbmc1k[["pca"]], dims = 1:5, nfeatures = 5)
```

    ## PC_ 1 
    ## Positive:  FCN1, LYZ, S100A9, FGL2, MNDA 
    ## Negative:  LTB, TRAC, CD3E, TRBC2, LCK 
    ## PC_ 2 
    ## Positive:  CD3E, IL32, CD247, GZMM, LCK 
    ## Negative:  CD79A, MS4A1, IGHM, BANK1, HLA-DQA1 
    ## PC_ 3 
    ## Positive:  CAVIN2, TUBB1, PF4, GNG11, PPBP 
    ## Negative:  MT-CO1, CYBA, FOSB, VIM, FOS 
    ## PC_ 4 
    ## Positive:  GZMB, CLIC3, GNLY, NKG7, KLRD1 
    ## Negative:  LEF1, TRABD2A, MAL, TCF7, RCAN3 
    ## PC_ 5 
    ## Positive:  GNLY, KLRD1, FGFBP2, KLRF1, ADGRG1 
    ## Negative:  SERPINF1, LILRA4, SMPD3, CLEC4C, SCT

``` r
# this is where PCA is performed on the scaled data, default is to only use those previously defined variable features. there seem to be 3 clusters?
VizDimLoadings(pbmc1k, dims = 1:2, reduction = "pca")
```

![](Practice-dataset_files/figure-gfm/QC-4.png)<!-- -->

``` r
DimPlot(pbmc1k, reduction = "pca")
```

![](Practice-dataset_files/figure-gfm/QC-5.png)<!-- -->

``` r
# if we just look at 1 PC across 500 cells, 
DimHeatmap(pbmc1k, dims = 1, cells = 500, balanced = TRUE)
```

![](Practice-dataset_files/figure-gfm/QC-6.png)<!-- -->

``` r
### --------- DETERMINE DATA DIMENSIONS ---------
pbmc1k <- JackStraw(pbmc1k, num.replicate = 100)
pbmc1k <- ScoreJackStraw(pbmc1k, dims = 1:20)

JackStrawPlot(pbmc1k, dims = 1:15)
```

    ## Warning: Removed 21680 rows containing missing values (geom_point).

![](Practice-dataset_files/figure-gfm/QC%20continued-1.png)<!-- -->

``` r
### --------- DETERMINE DATA DIMENSIONS ---------
ElbowPlot(pbmc1k)
```

![](Practice-dataset_files/figure-gfm/QC%20continued.-1.png)<!-- -->

``` r
### --------- CLUSTERING ---------
# As plotted above, select top 10 PC to have captured the most variation in this dataset.
# Hence the dimentionality of pbmc1k is the first 10 PC.
pbmc1k <- FindNeighbors(pbmc1k, dims = 1:10)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
# this uses KNN method

# understanding resolution
pbmc1k <- FindClusters(pbmc1k, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 1134
    ## Number of edges: 33968
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9686
    ## Number of communities: 4
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 1134
    ## Number of edges: 33968
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9220
    ## Number of communities: 8
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 1134
    ## Number of edges: 33968
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8901
    ## Number of communities: 11
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 1134
    ## Number of edges: 33968
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8659
    ## Number of communities: 11
    ## Elapsed time: 0 seconds
    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 1134
    ## Number of edges: 33968
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8322
    ## Number of communities: 11
    ## Elapsed time: 0 seconds

``` r
# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

# Number of nodes: 1134
# Number of edges: 35235

# Running Louvain algorithm...
# Maximum modularity in 10 random starts: 0.8902
# Number of communities: 9

View(pbmc1k@meta.data)
```

    ## Warning in system2("/usr/bin/otool", c("-L", shQuote(DSO)), stdout = TRUE):
    ## running command ''/usr/bin/otool' -L '/Library/Frameworks/R.framework/Resources/
    ## modules/R_de.so'' had status 1

``` r
res0.1 <- DimPlot(pbmc1k, group.by = "RNA_snn_res.0.1", label = TRUE)
res0.3 <- DimPlot(pbmc1k, group.by = "RNA_snn_res.0.3", label = TRUE)
res0.5 <- DimPlot(pbmc1k, group.by = "RNA_snn_res.0.5", label = TRUE)
res0.7 <- DimPlot(pbmc1k, group.by = "RNA_snn_res.0.7", label = TRUE)
res1.0 <- DimPlot(pbmc1k, group.by = "RNA_snn_res.1", label = TRUE)


res0.1
```

![](Practice-dataset_files/figure-gfm/QC%20continued..-1.png)<!-- -->

``` r
res0.3
```

![](Practice-dataset_files/figure-gfm/QC%20continued..-2.png)<!-- -->

``` r
res0.5
```

![](Practice-dataset_files/figure-gfm/QC%20continued..-3.png)<!-- -->

``` r
# setting identity of clusters for cells. Here we picked 0.3 resolution. (default was 0.5)
Idents(pbmc1k) <- "RNA_snn_res.0.3"
Idents(pbmc1k)
```

    ## AAACCCAAGGAGAGTA-1 AAACGCTTCAGCCCAG-1 AAAGAACAGACGACTG-1 AAAGAACCAATGGCAG-1 
    ##                  0                  2                  4                  4 
    ## AAAGAACGTCTGCAAT-1 AAAGGATAGTAGACAT-1 AAAGGATCACCGGCTA-1 AAAGGATTCAGCTTGA-1 
    ##                  1                  2                  2                  0 
    ## AAAGGATTCCGTTTCG-1 AAAGGGCTCATGCCCT-1 AAAGGGCTCCGTAGGC-1 AAAGGTACAACTGCTA-1 
    ##                  0                  3                  4                  0 
    ## AAAGTCCAGCGGGTTA-1 AAAGTCCCACCAGCCA-1 AAAGTGATCGTACACA-1 AAATGGAAGCCGCTTG-1 
    ##                  3                  2                  6                  0 
    ## AAATGGACAATGCTCA-1 AAATGGAGTACCGCGT-1 AAATGGATCCTATTTG-1 AACAAAGGTGATGAAT-1 
    ##                  0                  1                  2                  1 
    ## AACAACCAGTAGTCCT-1 AACAACCCACGCTATA-1 AACAAGAGTTATAGAG-1 AACAGGGGTGGGAGAG-1 
    ##                  0                  4                  0                  4 
    ## AACCAACAGCTTGTTG-1 AACCCAACAACTGATC-1 AACCCAAGTGGGCTTC-1 AACCCAATCTTACCGC-1 
    ##                  3                  4                  2                  6 
    ## AACCTGACATCCTATT-1 AACCTTTGTTTCGGCG-1 AACGAAAAGGTTGGTG-1 AACGAAACAGCGTTTA-1 
    ##                  0                  6                  5                  0 
    ## AACGAAATCCATTTAC-1 AACGGGAGTCGCTCGA-1 AACGTCAAGACGCCCT-1 AAGAACAAGACCATTC-1 
    ##                  1                  0                  0                  4 
    ## AAGAACAAGCCTCAGC-1 AAGAACATCTTGCAAG-1 AAGACAACACTTCAGA-1 AAGACAACAGATCACT-1 
    ##                  0                  3                  0                  0 
    ## AAGACAATCCGCGAGT-1 AAGCATCAGGTCCGAA-1 AAGCATCCAGTCTTCC-1 AAGCGTTCACTGATTG-1 
    ##                  3                  3                  0                  0 
    ## AAGCGTTTCGCTATTT-1 AAGGTAAAGGAAGTAG-1 AAGGTAAAGGCCACCT-1 AAGGTAAGTCCTTAAG-1 
    ##                  2                  4                  1                  1 
    ## AAGTACCAGCGCCTTG-1 AAGTACCAGCTCCATA-1 AAGTACCCAAAGAACT-1 AAGTGAACATCAGCAT-1 
    ##                  3                  3                  3                  4 
    ## AAGTGAATCGGAAGGT-1 AAGTTCGAGGATACAT-1 AAGTTCGGTCAACACT-1 AATCGACGTGAGACGT-1 
    ##                  0                  3                  0                  0 
    ## AATGACCGTGTCATCA-1 AATGGAACACTCATAG-1 AATGGAATCTGATGGT-1 AATTCCTTCCATTGGA-1 
    ##                  0                  2                  0                  0 
    ## ACAAAGAAGCCTCAAT-1 ACAAAGAGTTCCATTT-1 ACAACCATCTATCCAT-1 ACAACCATCTGCCCTA-1 
    ##                  0                  1                  3                  1 
    ## ACAAGCTGTAGGTCAG-1 ACACCAAAGAGTCAGC-1 ACACGCGGTGTTGCCG-1 ACACGCGTCACCTCGT-1 
    ##                  1                  3                  0                  6 
    ## ACACGCGTCGCCGTGA-1 ACACTGATCACTGATG-1 ACACTGATCTTTCCGG-1 ACAGAAAAGTATTAGG-1 
    ##                  0                  0                  1                  0 
    ## ACAGAAACACCAGCGT-1 ACAGCCGTCGCTGTTC-1 ACATCCCCACAAATGA-1 ACATCCCGTAAGATTG-1 
    ##                  0                  2                  1                  0 
    ## ACATCGATCTAAGAAG-1 ACATCGATCTGAGAGG-1 ACATCGATCTTGGGCG-1 ACATGCAGTTGCTTGA-1 
    ##                  6                  4                  0                  0 
    ## ACATTTCGTGCCGAAA-1 ACATTTCGTTTCTTAC-1 ACCACAAAGGCCTTGC-1 ACCACAAGTCTGTGGC-1 
    ##                  3                  2                  0                  2 
    ## ACCATTTCATTGGGAG-1 ACCCAAAGTTGGGAAC-1 ACCCTTGGTATCGTAC-1 ACCCTTGTCATTTCGT-1 
    ##                  1                  5                  1                  4 
    ## ACCTACCCAACGGTAG-1 ACCTGAACATGAGAAT-1 ACCTGTCTCAAGTCTG-1 ACCTGTCTCTCAACGA-1 
    ##                  4                  1                  3                  7 
    ## ACGATCAGTCGTTCAA-1 ACGATCATCCACGGAC-1 ACGATCATCCGTTGGG-1 ACGCACGAGAGAGTGA-1 
    ##                  4                  4                  3                  4 
    ## ACGGAAGGTGAACCGA-1 ACGGTCGTCACTACGA-1 ACGGTTAAGTGCTCAT-1 ACGTAACAGGTGCTTT-1 
    ##                  5                  1                  2                  1 
    ## ACGTAACTCATATGGC-1 ACGTAACTCTTCGATT-1 ACGTACAAGCTGTTCA-1 ACGTACACAAAGCTCT-1 
    ##                  2                  0                  4                  0 
    ## ACGTACACAGACGATG-1 ACGTAGTCACGCCACA-1 ACGTAGTCATTAGGAA-1 ACGTTCCGTGGGTCAA-1 
    ##                  6                  2                  4                  0 
    ## ACTATCTGTGTGTACT-1 ACTATGGAGAAATCCA-1 ACTCCCAAGGTTGGTG-1 ACTCCCATCGTCCTCA-1 
    ##                  6                  2                  6                  1 
    ## ACTGATGAGCAGCGAT-1 ACTGATGCACCAGCCA-1 ACTGTGACAACTCGAT-1 ACTTATCTCTAAGGAA-1 
    ##                  2                  2                  3                  0 
    ## ACTTCGCAGAACCGCA-1 ACTTCGCGTACTCCGG-1 ACTTTCAAGTATGACA-1 ACTTTCAGTACGCGTC-1 
    ##                  1                  2                  4                  1 
    ## AGAAGCGGTCTAGATC-1 AGAAGCGTCTGCGTCT-1 AGACAAAGTAGCTTGT-1 AGACACTCAAGTCGTT-1 
    ##                  1                  1                  0                  1 
    ## AGACCATCATTACTCT-1 AGACTCATCAACGCTA-1 AGACTCATCAGAGCGA-1 AGAGAATAGAGGGTCT-1 
    ##                  3                  1                  0                  4 
    ## AGAGAATCACCCAACG-1 AGAGAGCAGTTGTCGT-1 AGAGCAGGTCATAAAG-1 AGAGCAGGTTGAGTCT-1 
    ##                  3                  2                  0                  3 
    ## AGAGCAGTCGCGCTGA-1 AGAGCCCGTCACTCGG-1 AGATAGAGTGACACGA-1 AGATCGTCAGGACTTT-1 
    ##                  4                  3                  0                  4 
    ## AGATGAACAGGCACAA-1 AGATGCTGTACGGTTT-1 AGATGCTGTGTGGTCC-1 AGCATCAGTCGGTGTC-1 
    ##                  1                  0                  3                  0 
    ## AGCATCAGTGGTCCGT-1 AGCCAATAGTCCCTAA-1 AGCCACGCAGGTTCGC-1 AGCCACGGTTCCGCGA-1 
    ##                  3                  1                  3                  0 
    ## AGCCACGTCCGAGAAG-1 AGCCAGCGTAATGCGG-1 AGCCAGCGTAGTTAGA-1 AGCGTATTCAGTGGGA-1 
    ##                  4                  3                  5                  2 
    ## AGCGTCGCATTGACCA-1 AGCTACACAGTCGGAA-1 AGCTACAGTCGTTTCC-1 AGGAAATTCACCATGA-1 
    ##                  1                  1                  4                  0 
    ## AGGAATACAAACTCGT-1 AGGAATATCATTACGG-1 AGGACGAAGACCTTTG-1 AGGACGAAGATTAGTG-1 
    ##                  3                  3                  6                  1 
    ## AGGACGATCACTGATG-1 AGGACGATCCAGTACA-1 AGGACTTTCCTACCAC-1 AGGATAAAGAACAAGG-1 
    ##                  3                  4                  4                  3 
    ## AGGCCACAGGGCCAAT-1 AGGCTGCAGGGCTTCC-1 AGGCTGCGTGTGTGTT-1 AGGGAGTGTGTCGCTG-1 
    ##                  2                  2                  2                  3 
    ## AGGGCTCCAATACGCT-1 AGGGTCCCAGGACATG-1 AGGGTCCCAGTAGTTC-1 AGGGTTTCAGCGTGCT-1 
    ##                  0                  0                  4                  0 
    ## AGGTAGGAGGAAGTAG-1 AGGTAGGCAAAGCGTG-1 AGGTCATTCAAGTCGT-1 AGGTGTTAGCAATTAG-1 
    ##                  2                  0                  0                  5 
    ## AGGTGTTCATTAAGCC-1 AGGTGTTGTCGATTAC-1 AGGTGTTTCAATGCAC-1 AGGTTACAGGTGCGAT-1 
    ##                  0                  1                  3                  3 
    ## AGGTTACGTCTGATAC-1 AGTAACCCAAGCTCTA-1 AGTAACCTCGTTGTAG-1 AGTACCAGTGGGTCAA-1 
    ##                  7                  0                  3                  0 
    ## AGTACCATCAGCTTGA-1 AGTAGCTAGGAGTACC-1 AGTAGCTCACTGCTTC-1 AGTAGCTCAGGCACTC-1 
    ##                  6                  2                  1                  0 
    ## AGTAGCTGTTAGAAAC-1 AGTAGTCGTATGACAA-1 AGTAGTCGTCGACTGC-1 AGTCAACAGTCATGCT-1 
    ##                  0                  1                  5                  0 
    ## AGTCAACGTACACGCC-1 AGTCACATCTCCCATG-1 AGTCTCCCAACCGTGC-1 AGTCTCCGTGAGTCAG-1 
    ##                  5                  1                  2                  4 
    ## AGTGCCGCATCTGGGC-1 AGTTAGCAGTCAGGGT-1 AGTTAGCGTACCATAC-1 AGTTAGCGTTCATCGA-1 
    ##                  4                  4                  2                  4 
    ## AGTTCCCCAAACAGGC-1 AGTTCCCGTTACACTG-1 AGTTCGAGTATGCGTT-1 ATACCGACACACGGTC-1 
    ##                  4                  3                  2                  1 
    ## ATACTTCGTACGACTT-1 ATACTTCGTTCGGTAT-1 ATACTTCTCAGCTGAT-1 ATACTTCTCGCCGTGA-1 
    ##                  0                  3                  1                  2 
    ## ATAGACCTCCGCTTAC-1 ATAGGCTGTGTTACAC-1 ATATCCTGTATGAGCG-1 ATCACAGAGACGTCCC-1 
    ##                  4                  0                  2                  0 
    ## ATCACAGCATGACGTT-1 ATCACGATCCGCAACG-1 ATCACTTGTGGCTGCT-1 ATCCACCCAGCAGTGA-1 
    ##                  0                  4                  0                  6 
    ## ATCCACCCATGTCTAG-1 ATCCACCGTTTGTTGG-1 ATCCATTAGACGCAGT-1 ATCCCTGCAATTTCTC-1 
    ##                  6                  0                  2                  4 
    ## ATCCGTCCAACTCCCT-1 ATCCTATGTCTGTGCG-1 ATCGATGCAACGGCTC-1 ATCGATGGTCAATCTG-1 
    ##                  2                  2                  2                  6 
    ## ATCGGATAGGAAAGGT-1 ATCGTAGCATCGGAGA-1 ATCGTCCAGCCTGAAG-1 ATCGTCCGTATGAGGC-1 
    ##                  1                  0                  0                  0 
    ## ATCGTCCTCGGAGCAA-1 ATCGTGAAGCAGATAT-1 ATCTCTACAAGAAACT-1 ATCTTCAAGAAGAGCA-1 
    ##                  7                  0                  0                  0 
    ## ATCTTCACATACTGAC-1 ATCTTCAGTCTCTCCA-1 ATCTTCATCTCATAGG-1 ATGAAAGTCCTTATGT-1 
    ##                  0                  0                  4                  1 
    ## ATGACCAAGGCTAACG-1 ATGACCACAATATCCG-1 ATGAGTCTCCAGTTCC-1 ATGATCGGTCGAACAG-1 
    ##                  0                  1                  2                  2 
    ## ATGCGATGTGCTCTTC-1 ATGCGATGTTGTAGCT-1 ATGGATCGTAGATTGA-1 ATGGATCTCGACCTAA-1 
    ##                  2                  1                  2                  1 
    ## ATGGGAGAGACGTCCC-1 ATGGGAGCACTCCGGA-1 ATGGGAGTCCACTGAA-1 ATGGGAGTCCCATACC-1 
    ##                  1                  0                  0                  3 
    ## ATGGTTGTCACGGAGA-1 ATGTCCCCAGGAATCG-1 ATGTCCCCATAGGTTC-1 ATGTCCCGTCCTACGG-1 
    ##                  0                  2                  2                  0 
    ## ATTACTCCAGCTTTGA-1 ATTACTCGTTGCAACT-1 ATTCACTAGAAGGCTC-1 ATTCACTTCTCCGCAT-1 
    ##                  1                  2                  4                  6 
    ## ATTCAGGCAAGTATCC-1 ATTCAGGCATAGCTGT-1 ATTCAGGTCTTCGTGC-1 ATTCCCGAGTCACTCA-1 
    ##                  0                  0                  1                  1 
    ## ATTCCTACATGAGATA-1 ATTCGTTCAGTCTTCC-1 ATTCGTTGTAACACCT-1 ATTCGTTGTGTAGCAG-1 
    ##                  0                  1                  3                  2 
    ## ATTCGTTGTGTTCAGT-1 ATTCGTTGTTCTTGTT-1 ATTCTACTCTCGTGAA-1 ATTCTTGGTACAATAG-1 
    ##                  1                  0                  3                  4 
    ## ATTTCACAGTGCAAAT-1 ATTTCACCAGTGGTGA-1 ATTTCACGTGGATACG-1 ATTTCTGCAAGGTCTT-1 
    ##                  0                  5                  0                  0 
    ## ATTTCTGGTATCCTCC-1 ATTTCTGTCGATGCAT-1 ATTTCTGTCGGCCAAC-1 CAAAGAACACATATCG-1 
    ##                  2                  0                  1                  1 
    ## CAAAGAACAGCTGTTA-1 CAAAGAAGTTCTCCCA-1 CAACAACTCGATACTG-1 CAACAGTCATCCGTTC-1 
    ##                  4                  0                  4                  3 
    ## CAACAGTTCGAGCACC-1 CAACCAAGTTAGAAAC-1 CAACCTCGTATCGGTT-1 CAACGATAGACGGAAA-1 
    ##                  3                  0                  1                  0 
    ## CAACGATAGCCTAGGA-1 CAACGATCAAAGCACG-1 CAACGATTCAAGCCCG-1 CAAGACTAGTGGAAAG-1 
    ##                  0                  1                  0                  4 
    ## CAAGACTTCAATCGGT-1 CAAGAGGTCTTCACGC-1 CAAGCTAGTGGCTGCT-1 CAAGGGAAGAGTCGAC-1 
    ##                  0                  0                  0                  1 
    ## CAATACGAGGCAGCTA-1 CAATCGAAGAGGATCC-1 CAATCGATCCGTACGG-1 CAATGACAGATGTTAG-1 
    ##                  4                  1                  1                  5 
    ## CAATGACTCAAGCCCG-1 CAATTTCAGATTGCGG-1 CACACAACAAGAAACT-1 CACACAACACGAGGAT-1 
    ##                  3                  2                  3                  0 
    ## CACAGATCACTCCGAG-1 CACAGATCATCATGAC-1 CACAGATGTTCCGTTC-1 CACAGGCCAACGATCT-1 
    ##                  0                  0                  7                  0 
    ## CACATGAGTTTCGACA-1 CACCAAAAGGACATCG-1 CACCGTTCACCCTAGG-1 CACGAATCACATAGCT-1 
    ##                  1                  0                  3                  1 
    ## CACGAATTCAAGTCGT-1 CACGGGTAGCATTTGC-1 CACGGGTTCCCGATCT-1 CACGTTCCAGCTTCCT-1 
    ##                  1                  0                  2                  0 
    ## CACGTTCTCTATCACT-1 CACTAAGGTGATAGTA-1 CACTGAAGTCTTTCTA-1 CACTGAATCCTGATAG-1 
    ##                  1                  1                  2                  5 
    ## CACTGGGAGGCGACAT-1 CACTGGGAGTTACTCG-1 CACTGGGGTCGCACAC-1 CACTGTCCACATATGC-1 
    ##                  0                  0                  0                  2 
    ## CACTTCGTCACCTACC-1 CAGAGCCGTGCCCGTA-1 CAGATACTCGTTCCCA-1 CAGATTGAGCTGAGCA-1 
    ##                  0                  3                  1                  1 
    ## CAGATTGGTGACTAAA-1 CAGCAATGTATCGTTG-1 CAGCAGCTCAGAGTTC-1 CAGCCAGGTGGGTTGA-1 
    ##                  2                  7                  4                  0 
    ## CAGCGTGTCCATCACC-1 CAGGCCAGTTCAATCG-1 CAGGCCATCACTCACC-1 CAGGCCATCCACTAGA-1 
    ##                  2                  1                  4                  0 
    ## CAGGGCTAGATTGGGC-1 CAGGGCTTCGACGAGA-1 CAGGTATGTCAGTCTA-1 CAGGTATGTGTAAATG-1 
    ##                  4                  2                  5                  2 
    ## CAGTTCCCAATCGTCA-1 CATAAGCGTAGCGCTC-1 CATAAGCTCCGTGTGG-1 CATACTTTCAACGAGG-1 
    ##                  4                  3                  0                  1 
    ## CATCAAGGTATGGTTC-1 CATCAAGGTCTAACTG-1 CATCAAGGTTGTCCCT-1 CATCCGTGTAGAGTTA-1 
    ##                  2                  2                  2                  3 
    ## CATCCGTGTATCGCTA-1 CATCGGGTCGCTACAA-1 CATCGTCTCCGTTGAA-1 CATGAGTCACCATTCC-1 
    ##                  3                  2                  4                  4 
    ## CATGCCTCAATTGCTG-1 CATGCCTGTTCCGCAG-1 CATGCGGTCGGTCGAC-1 CATTCTAGTAGCTTGT-1 
    ##                  2                  0                  3                  4 
    ## CATTCTAGTGCCGGTT-1 CATTCTATCAGTAGGG-1 CATTGAGAGGGACCAT-1 CATTTCATCTCTCTTC-1 
    ##                  5                  0                  0                  0 
    ## CCAAGCGAGCTAATCC-1 CCAAGCGAGGTTCATC-1 CCAATGAAGGTCACTT-1 CCAATTTCAAAGGCTG-1 
    ##                  3                  0                  2                  4 
    ## CCAATTTCATTCGATG-1 CCACACTCAGGACTTT-1 CCACCATAGGTCATCT-1 CCACGAGGTAATCAAG-1 
    ##                  4                  7                  1                  1 
    ## CCACGAGTCTTCCTAA-1 CCACGTTTCATCAGTG-1 CCACGTTTCGTGAGAG-1 CCACGTTTCTGCTGAA-1 
    ##                  0                  2                  1                  2 
    ## CCACTTGTCGCAGAGA-1 CCATAAGGTGAGGAAA-1 CCATCACCATTAAAGG-1 CCATCACTCTAACACG-1 
    ##                  7                  0                  2                  0 
    ## CCCGGAATCAAACCCA-1 CCCTAACCAGGTTCGC-1 CCCTAACGTAACTGCT-1 CCCTTAGGTTTCCAAG-1 
    ##                  3                  3                  0                  1 
    ## CCGAACGAGCTAAGTA-1 CCGAACGAGGAACGCT-1 CCGATCTGTGGAGGTT-1 CCGATGGGTTCTCGCT-1 
    ##                  1                  5                  0                  5 
    ## CCGATGGTCTGTCGCT-1 CCGCAAGCATTCAGGT-1 CCGCAAGGTACCAGAG-1 CCGGACAGTGATACCT-1 
    ##                  2                  3                  1                  0 
    ## CCGGGTACACTCCGAG-1 CCGGTAGGTAGTACGG-1 CCGTAGGGTTCCCACT-1 CCGTAGGTCATACAGC-1 
    ##                  0                  6                  6                  2 
    ## CCGTGAGAGAGGTCAC-1 CCGTGAGTCGCCTTTG-1 CCGTTCACAGAAGCTG-1 CCTAAGAAGACTCTAC-1 
    ##                  4                  2                  1                  4 
    ## CCTATCGAGACAGCGT-1 CCTCAACTCTTAGCCC-1 CCTCAGTAGCAAGGAA-1 CCTCAGTCATGATGCT-1 
    ##                  0                  3                  0                  3 
    ## CCTCAGTTCCGAGCTG-1 CCTCATGCACTATGTG-1 CCTCATGGTTTCGCTC-1 CCTCATGTCTCTCAAT-1 
    ##                  2                  6                  1                  0 
    ## CCTCCAAAGGCCCGTT-1 CCTCCAAAGTATGATG-1 CCTCCAACAACCCTCT-1 CCTCCAACACAATGAA-1 
    ##                  0                  2                  4                  2 
    ## CCTCCTCAGTCCCGAC-1 CCTCCTCCATTCAGGT-1 CCTCCTCTCAAGCGTT-1 CCTCTAGCAGACCTAT-1 
    ##                  2                  3                  3                  0 
    ## CCTCTCCTCGTTCCTG-1 CCTGCATAGTGGTGAC-1 CCTGCATCAACGGGTA-1 CCTGCATGTTCTTAGG-1 
    ##                  0                  5                  3                  1 
    ## CCTGTTGGTGTGCCTG-1 CCTGTTGTCAAGAGTA-1 CCTTCAGAGGAACGTC-1 CCTTCAGCACAGAGCA-1 
    ##                  0                  2                  5                  2 
    ## CGAAGGATCGACGCGT-1 CGAATTGGTCTGTTAG-1 CGAGAAGTCCAATGCA-1 CGAGGAAAGTCGCGAA-1 
    ##                  1                  3                  5                  3 
    ## CGAGGAACACTGATTG-1 CGAGGCTCAGATCACT-1 CGAGGCTCAGCATGCC-1 CGAGTGCGTTCCGCGA-1 
    ##                  2                  2                  0                  1 
    ## CGAGTTAGTCGAGTTT-1 CGATCGGCAATTGGTC-1 CGATCGGCAGCGTGAA-1 CGATGCGAGCAGAAAG-1 
    ##                  0                  0                  3                  3 
    ## CGATGCGTCGAACCAT-1 CGATGGCCACTCGATA-1 CGATGGCGTAGCGTCC-1 CGCAGGTGTGCTTCAA-1 
    ##                  2                  1                  1                  2 
    ## CGCATAAAGTGATGGC-1 CGCATGGAGGCCACTC-1 CGCGTGAGTACAATAG-1 CGGAACCGTATCACGT-1 
    ##                  0                  5                  2                  0 
    ## CGGAACCGTTTGCAGT-1 CGGACACAGGACAGCT-1 CGGGACTGTTTCACTT-1 CGGGTGTTCCGTGACG-1 
    ##                  2                  0                  6                  0 
    ## CGTAAGTGTGATAGTA-1 CGTAGTATCAGCCCAG-1 CGTCAAACAGAATGTA-1 CGTCCATGTTTGGAAA-1 
    ##                  1                  4                  6                  3 
    ## CGTGAATAGTCGGCCT-1 CGTGCTTAGCTTTGTG-1 CGTGCTTCAAGCGCTC-1 CGTGCTTGTTTAGACC-1 
    ##                  1                  0                  0                  1 
    ## CGTGTCTCACAGTACT-1 CGTGTCTCATGGGTTT-1 CGTGTCTGTCAGTTTG-1 CGTTAGATCAGCCTCT-1 
    ##                  1                  1                  0                  7 
    ## CGTTAGATCATTATCC-1 CGTTAGATCCAAGAGG-1 CGTTCTGTCAGCGCAC-1 CGTTCTGTCCCGAACG-1 
    ##                  0                  0                  3                  6 
    ## CTAACCCAGCCGCTTG-1 CTAACCCCACTTCAAG-1 CTAAGTGCACAATGCT-1 CTAAGTGGTCCACTCT-1 
    ##                  3                  1                  2                  4 
    ## CTACAGACAACAAAGT-1 CTACCTGTCTCGAACA-1 CTACGGGTCGTTCAGA-1 CTAGACAGTCTAACTG-1 
    ##                  2                  0                  4                  0 
    ## CTAGGTATCGTAGGAG-1 CTAGGTATCTAACGCA-1 CTATCCGAGCTGTTAC-1 CTCAACCGTGGATGAC-1 
    ##                  5                  2                  3                  2 
    ## CTCAACCGTGTGTCGC-1 CTCAAGAAGATTACCC-1 CTCAAGAGTCAAAGAT-1 CTCAATTTCGGAAGGT-1 
    ##                  2                  1                  6                  0 
    ## CTCACTGCACCTTCCA-1 CTCAGGGGTCATCGCG-1 CTCAGGGTCTCACCCA-1 CTCAGTCAGGGACAGG-1 
    ##                  3                  0                  7                  6 
    ## CTCATCGAGAGGCGTT-1 CTCATCGTCGTACCTC-1 CTCATTAAGGAAGTGA-1 CTCATTATCGTACACA-1 
    ##                  3                  2                  0                  5 
    ## CTCCAACCAATCTCTT-1 CTCCAACGTAGTGTGG-1 CTCCATGAGACAGCTG-1 CTCCATGGTGTTGACT-1 
    ##                  4                  1                  2                  0 
    ## CTCCATGTCATCGTAG-1 CTCCATGTCGTCGATA-1 CTCCCAATCAATCAGC-1 CTCCCTCAGCGCGTTC-1 
    ##                  3                  0                  6                  0 
    ## CTCCCTCCACCGGAAA-1 CTCCCTCGTTTACGAC-1 CTCCCTCTCATCGCTC-1 CTCCGATTCCTTTAGT-1 
    ##                  0                  4                  4                  4 
    ## CTCCGATTCTTTCTAG-1 CTCCTCCTCACCTCTG-1 CTCCTTTAGATGCTGG-1 CTCCTTTTCTACACTT-1 
    ##                  3                  3                  0                  3 
    ## CTCGAGGAGCCTATTG-1 CTCGAGGAGGAAAGGT-1 CTCTCAGGTTAAGCAA-1 CTCTCGATCCATAGGT-1 
    ##                  2                  4                  3                  2 
    ## CTCTGGTAGTAGGATT-1 CTGAATGAGCTCTATG-1 CTGAATGAGGCGAACT-1 CTGAGGCGTTCAGCGC-1 
    ##                  0                  0                  3                  2 
    ## CTGATCCGTCGCGTCA-1 CTGATCCTCCCTCTAG-1 CTGATCCTCTCTCCGA-1 CTGCAGGAGCCATTGT-1 
    ##                  1                  4                  5                  1 
    ## CTGCATCAGCGTCAAG-1 CTGCATCGTAAGCTCT-1 CTGCCATTCAGTGTTG-1 CTGCCATTCGGAACTT-1 
    ##                  0                  1                  1                  0 
    ## CTGCCTACAAGACTGG-1 CTGCCTATCTTGATTC-1 CTGCGAGAGATTAGAC-1 CTGCGAGAGTGCACTT-1 
    ##                  0                  2                  0                  0 
    ## CTGCGAGCAGATTCGT-1 CTGCGAGGTAGGTTTC-1 CTGCTCACAGCGTTTA-1 CTGGACGCAACGATTC-1 
    ##                  2                  0                  0                  1 
    ## CTGGACGCACACGGTC-1 CTGGACGGTTCCCAAA-1 CTGGACGTCTCATTTG-1 CTGGCAGGTTGCAAGG-1 
    ##                  0                  1                  0                  0 
    ## CTGGCAGTCGTCTCAC-1 CTGTCGTAGCTGGCCT-1 CTGTGAAAGGAACGCT-1 CTGTGAATCATCGCTC-1 
    ##                  1                  3                  0                  0 
    ## CTGTGGGGTCGTGATT-1 CTTAGGAGTCTGTGGC-1 CTTCAATGTTGGACTT-1 CTTCAATTCGAGAGCA-1 
    ##                  2                  2                  3                  5 
    ## CTTCGGTCACGTTCGG-1 CTTCGGTTCCAAGAGG-1 CTTCTCTTCCTACAAG-1 CTTGAGAAGAATTGCA-1 
    ##                  2                  2                  1                  3 
    ## CTTGATTAGGTAGGCT-1 GAACACTGTTGGACTT-1 GAACGTTGTCAGGAGT-1 GAACTGTCATCGAACT-1 
    ##                  1                  2                  0                  3 
    ## GAACTGTTCTCAAAGC-1 GAAGAATAGGTACAAT-1 GAAGCCCCAGGAACCA-1 GAAGCGATCAGTAGGG-1 
    ##                  1                  0                  0                  4 
    ## GAAGGACGTCTGCCTT-1 GAAGGGTTCGCAGATT-1 GAAGTAAGTCTAGATC-1 GAATCACCAGTGGTGA-1 
    ##                  0                  1                  1                  1 
    ## GACACGCCAGTCCGTG-1 GACACGCGTTCCAAAC-1 GACAGCCGTATCGTGT-1 GACATCAGTAACGATA-1 
    ##                  2                  2                  4                  0 
    ## GACCCTTCAGCGTATT-1 GACCCTTCATCATCTT-1 GACCTTCCACCTGCGA-1 GACCTTCGTCGAGCTC-1 
    ##                  7                  1                  0                  3 
    ## GACTATGCATCATTGG-1 GACTATGTCCAAGGGA-1 GACTCAATCTCCCAAC-1 GACTGATTCATCCTAT-1 
    ##                  0                  1                  0                  3 
    ## GACTTCCTCCTCTAGC-1 GAGAAATCATTGGATC-1 GAGACCCAGTAGGATT-1 GAGACTTAGTCACTCA-1 
    ##                  1                  1                  0                  6 
    ## GAGACTTCAGGCTATT-1 GAGAGGTTCATAGACC-1 GAGATGGCATGTAACC-1 GAGCCTGAGGCAGCTA-1 
    ##                  0                  3                  3                  5 
    ## GAGCCTGCAGTCTCTC-1 GAGCCTGGTCCACGCA-1 GAGCTGCGTCGGATTT-1 GAGGCAAAGTCACGAG-1 
    ##                  4                  1                  0                  1 
    ## GAGGGATCACTTGGCG-1 GAGGGTACAGTTTCAG-1 GAGGGTAGTCAACACT-1 GAGGGTATCTGTCCCA-1 
    ##                  3                  1                  3                  1 
    ## GAGTGAGGTCGGCTAC-1 GAGTGTTCACAAGCTT-1 GAGTGTTGTACGTAGG-1 GAGTTACAGACGACGT-1 
    ##                  2                  2                  4                  1 
    ## GAGTTGTAGCGTTACT-1 GAGTTGTAGTGCGTCC-1 GAGTTGTCACGGTGCT-1 GATAGCTGTCGGAAAC-1 
    ##                  0                  0                  4                  0 
    ## GATAGCTGTCTGTTAG-1 GATAGCTTCGTGCAGC-1 GATCAGTAGGAACATT-1 GATCAGTAGTGCGCTC-1 
    ##                  1                  2                  2                  3 
    ## GATCATGAGGCCTAAG-1 GATCCCTAGAGGATGA-1 GATCCCTGTAATTAGG-1 GATCGTACAGCGCGTT-1 
    ##                  1                  3                  3                  1 
    ## GATCGTATCTCTGGTC-1 GATGACTAGAGGTTAT-1 GATGACTTCGTGGCGT-1 GATGACTTCTCTAGGA-1 
    ##                  6                  2                  3                  2 
    ## GATGAGGCAAATGCGG-1 GATGATCGTATGGAAT-1 GATGATCGTGTAACGG-1 GATGCTAGTACAGGTG-1 
    ##                  3                  4                  2                  3 
    ## GATGGAGAGCCAGACA-1 GATGGAGCACTTGAAC-1 GATGGAGTCTGCGGCA-1 GATGGAGTCTTCCAGC-1 
    ##                  2                  1                  3                  0 
    ## GATTCTTGTACAGTCT-1 GATTGGTGTTGAGAGC-1 GATTTCTCAAGTGGGT-1 GATTTCTTCCTTCAGC-1 
    ##                  1                  1                  4                  2 
    ## GCAACCGTCATTTCCA-1 GCAACCGTCGAACGGA-1 GCACGGTAGTAGTCTC-1 GCACGGTCACCCAACG-1 
    ##                  4                  2                  1                  2 
    ## GCACGTGTCAGACATC-1 GCAGCCAAGAGGTCAC-1 GCAGCCACATACATCG-1 GCAGCCAGTTGCGTAT-1 
    ##                  6                  2                  2                  2 
    ## GCAGCCATCACCGGGT-1 GCAGCTGGTAAGTTAG-1 GCAGCTGGTGGGTCAA-1 GCAGCTGTCAAAGCCT-1 
    ##                  4                  0                  2                  5 
    ## GCAGCTGTCCTTCTGG-1 GCAGTTACACATTCTT-1 GCAGTTATCTCTCTTC-1 GCATCGGGTAACGATA-1 
    ##                  4                  0                  7                  0 
    ## GCCAACGTCTAACGCA-1 GCCAACGTCTGCATAG-1 GCCAGGTGTTAGAAAC-1 GCCAGTGAGTCCTACA-1 
    ##                  6                  0                  3                  2 
    ## GCCGATGAGTAACCTC-1 GCCGATGCATCTTTCA-1 GCCGTGACAACGACTT-1 GCCGTGACATGCCGGT-1 
    ##                  3                  2                  0                  4 
    ## GCCGTGAGTCAGTCCG-1 GCCTGTTCACAGTCAT-1 GCGAGAAAGTTGTAGA-1 GCGAGAATCCTCAGAA-1 
    ##                  5                  0                  7                  3 
    ## GCGGAAAGTGTGTCCG-1 GCGGATCAGCAACTTC-1 GCGTTTCAGTCTAGCT-1 GCTACCTCACGCACCA-1 
    ##                  4                  0                  1                  0 
    ## GCTACCTGTCCTACAA-1 GCTGAATAGCACTCTA-1 GCTGAATCACACGGAA-1 GCTGAATTCGTTAGTG-1 
    ##                  0                  2                  2                  0 
    ## GCTGCAGAGGGCAGAG-1 GCTGCAGGTGCTGATT-1 GCTGCAGGTTGACGGA-1 GCTGGGTAGTACAACA-1 
    ##                  2                  3                  0                  6 
    ## GCTTCACCAACACGAG-1 GCTTGGGAGAAGGTAG-1 GCTTTCGAGTGCAAAT-1 GCTTTCGGTTCGGGTC-1 
    ##                  5                  0                  0                  3 
    ## GCTTTCGTCCGGTAGC-1 GCTTTCGTCTGGGTCG-1 GGAACCCAGCCATTTG-1 GGAACCCTCTAATTCC-1 
    ##                  1                  0                  0                  2 
    ## GGAAGTGCAACCCTCT-1 GGAATCTCAAACTGCT-1 GGAATCTTCAATCTCT-1 GGAATGGAGGACTAAT-1 
    ##                  0                  3                  1                  0 
    ## GGAATGGCACGACGCT-1 GGACGTCGTTCAACGT-1 GGAGAACGTGCCCACA-1 GGAGGATGTATGAAAC-1 
    ##                  3                  0                  0                  0 
    ## GGAGGTAAGGCGTTGA-1 GGAGGTATCTGGTGCG-1 GGATCTAGTGCCTGCA-1 GGATCTATCACCTACC-1 
    ##                  0                  0                  5                  1 
    ## GGATCTATCATAAGGA-1 GGATGTTAGCCGGATA-1 GGCGTCACATTGTGCA-1 GGCTGTGCAAAGCTAA-1 
    ##                  4                  1                  4                  6 
    ## GGCTGTGCACGGTGAA-1 GGCTGTGGTCCACATA-1 GGCTTGGCACCGAATT-1 GGCTTGGGTAGATCGG-1 
    ##                  2                  0                  2                  2 
    ## GGCTTGGTCTTCTTCC-1 GGGACAATCCGATCTC-1 GGGACCTGTTCTTGTT-1 GGGACTCCACATATCG-1 
    ##                  0                  0                  1                  3 
    ## GGGACTCCACGGGCTT-1 GGGAGATCACGGGTAA-1 GGGAGATCATGCAGGA-1 GGGAGATTCACGGACC-1 
    ##                  0                  1                  0                  2 
    ## GGGAGTAAGCAGGCAT-1 GGGAGTATCACCTCGT-1 GGGAGTATCTTAGCCC-1 GGGATCCAGGGCAATC-1 
    ##                  0                  5                  0                  6 
    ## GGGATGAAGGGAGGGT-1 GGGATGAAGTCCCGAC-1 GGGCCATTCAGACCGC-1 GGGCTACAGACCAAAT-1 
    ##                  3                  1                  0                  3 
    ## GGGCTACAGGCAGTCA-1 GGGCTCATCGCAATTG-1 GGGTAGAAGTAGGGTC-1 GGGTCACCAGGAGGTT-1 
    ##                  3                  1                  3                  1 
    ## GGGTCACGTCCTTAAG-1 GGGTCACTCAATCTCT-1 GGGTCACTCACAAGGG-1 GGGTCTGGTTTGAAAG-1 
    ##                  0                  5                  0                  1 
    ## GGGTCTGTCCCGTTCA-1 GGGTGAACATAAGCAA-1 GGGTTATCAGCCATTA-1 GGGTTTACACAGAGCA-1 
    ##                  0                  1                  1                  0 
    ## GGGTTTACACTAAACC-1 GGTAATCCAAAGGAGA-1 GGTAATCCACAGCATT-1 GGTAATCCATACCACA-1 
    ##                  2                  4                  0                  0 
    ## GGTAATCTCAAGTTGC-1 GGTAGAGAGACCAGCA-1 GGTCACGCAATGGCCC-1 GGTCACGGTCGTGATT-1 
    ##                  0                  0                  1                  3 
    ## GGTCTGGTCAGACATC-1 GGTCTGGTCTCCCAAC-1 GGTGATTAGCTGTACT-1 GGTGATTCAGCTCGGT-1 
    ##                  2                  4                  2                  0 
    ## GGTGGCTGTGCGCTCA-1 GGTGTCGCAACGATCT-1 GGTGTCGGTGTTATCG-1 GGTGTTAGTGGATTTC-1 
    ##                  6                  0                  2                  2 
    ## GGTTCTCCAGGTTCCG-1 GGTTGTATCGAAGGAC-1 GTAACACAGAGCCGTA-1 GTAACACAGCTCCGAC-1 
    ##                  0                  2                  5                  7 
    ## GTAACACGTAGATGTA-1 GTAACACGTGGCGTAA-1 GTAAGTCTCGATTGGT-1 GTAATCGCAGTTTCGA-1 
    ##                  4                  0                  3                  4 
    ## GTAATCGGTTGCGAAG-1 GTAATGCCAATCGCGC-1 GTAATGCCACCAGTAT-1 GTACAACAGATAGCTA-1 
    ##                  6                  5                  0                  6 
    ## GTAGAAACATGGACAG-1 GTAGAAATCTGAGATC-1 GTAGAGGCAACTTCTT-1 GTAGATCCATCCCGTT-1 
    ##                  0                  0                  3                  1 
    ## GTAGATCTCCCGAAAT-1 GTAGCTATCTCATGGA-1 GTAGGAGTCCGTCACT-1 GTAGGTTCAGAAGCGT-1 
    ##                  2                  0                  6                  3 
    ## GTAGGTTGTGTTCGTA-1 GTAGTACCACCCTATC-1 GTAGTACTCTGTGTGA-1 GTATTTCAGAAGAGCA-1 
    ##                  4                  0                  5                  3 
    ## GTCAAACCAAGCACCC-1 GTCACGGAGCGGGTTA-1 GTCACGGCACCGTACG-1 GTCACTCTCTTGAACG-1 
    ##                  1                  3                  1                  1 
    ## GTCATGAAGTTTGGCT-1 GTCATGACATTGGATC-1 GTCATTTCAAGAGAGA-1 GTCATTTTCCTAGAGT-1 
    ##                  0                  2                  3                  4 
    ## GTCCCATCATGCCGAC-1 GTCCTCAAGTCACGAG-1 GTCGAATAGGAGTATT-1 GTCGAATCAGAACCGA-1 
    ##                  3                  2                  2                  0 
    ## GTCGAATGTCAAGCGA-1 GTCGCGAAGCAGCGAT-1 GTCGCGAGTTTACCTT-1 GTCTACCCAGGAAGTC-1 
    ##                  0                  4                  6                  3 
    ## GTCTAGAGTATGTCAC-1 GTCTAGATCGAGAGAC-1 GTCTGTCAGATACGAT-1 GTCTGTCTCGAGCCTG-1 
    ##                  5                  2                  1                  6 
    ## GTGAGCCCATCTTCGC-1 GTGAGCCGTGATCGTT-1 GTGAGCCTCTGGAAGG-1 GTGAGGAAGCATTGAA-1 
    ##                  1                  4                  0                  1 
    ## GTGAGGACACAACATC-1 GTGAGGACACTTTAGG-1 GTGAGTTTCCAATGCA-1 GTGATGTAGACATCAA-1 
    ##                  4                  3                  1                  1 
    ## GTGCGTGAGCTTCTAG-1 GTGCTGGCAGGCTACC-1 GTGCTGGTCAAACCTG-1 GTGCTTCAGACCATAA-1 
    ##                  0                  4                  3                  0 
    ## GTGCTTCCACTCATAG-1 GTGCTTCGTCGCTTGG-1 GTGGAAGAGCCAAGTG-1 GTGGCGTAGAATAACC-1 
    ##                  6                  0                  4                  4 
    ## GTGGCGTCAATCGAAA-1 GTGGGAAAGGTCGTAG-1 GTGGGAACATGAGTAA-1 GTGGTTAAGTAACCGG-1 
    ##                  4                  2                  1                  1 
    ## GTGTAACTCTTGGTGA-1 GTGTCCTAGTCGGCCT-1 GTGTCCTCAGTTGTTG-1 GTGTCCTGTGGCCTCA-1 
    ##                  1                  6                  2                  3 
    ## GTGTGATCACTACGGC-1 GTGTGATTCTTACGTT-1 GTGTGGCAGGTACTGG-1 GTGTGGCAGGTAGTCG-1 
    ##                  1                  6                  5                  4 
    ## GTGTGGCTCACCACAA-1 GTGTTAGAGACGGTCA-1 GTGTTCCAGGTGCCTC-1 GTGTTCCGTAGTTCCA-1 
    ##                  3                  5                  2                  5 
    ## GTTACAGAGTGTAGTA-1 GTTACAGCAAATCAGA-1 GTTACGAGTATGAGAT-1 GTTAGTGAGAACGTGC-1 
    ##                  3                  1                  0                  1 
    ## GTTAGTGCACCTCTAC-1 GTTATGGTCGTTCCTG-1 GTTCATTTCCGATGCG-1 GTTCCGTAGTTGCCCG-1 
    ##                  0                  1                  5                  0 
    ## GTTCCGTCACCCGTAG-1 GTTCTATCAAATGAGT-1 GTTCTATCACGTGTGC-1 GTTCTATTCTACACTT-1 
    ##                  4                  2                  0                  6 
    ## GTTGAACAGGAACATT-1 GTTGAACGTACCCGCA-1 GTTGCGGAGCTGCCTG-1 GTTGCTCAGTTTCGAC-1 
    ##                  6                  0                  3                  6 
    ## GTTGCTCTCATTTCCA-1 GTTGTAGCAGTTCTAG-1 GTTGTGATCTCGGTCT-1 GTTTGGAGTGTATACC-1 
    ##                  1                  4                  3                  1 
    ## TAACGACCACATATGC-1 TAACTTCGTTGTGGAG-1 TAAGCACAGAAGCTGC-1 TAAGCACTCACCATAG-1 
    ##                  3                  6                  1                  2 
    ## TAAGCCATCTGGCCTT-1 TAAGCGTTCGTTGCCT-1 TAAGCGTTCTGTCGTC-1 TAATTCCCACTCCGGA-1 
    ##                  2                  2                  1                  1 
    ## TAATTCCGTATTGACC-1 TAATTCCGTCAGTCCG-1 TACAACGCAAAGGTTA-1 TACAACGCATAGATCC-1 
    ##                  7                  1                  3                  3 
    ## TACACCCTCGACTCCT-1 TACAGGTTCCTTACCG-1 TACATTCGTGGCTTAT-1 TACCGGGTCCTCGATC-1 
    ##                  3                  1                  1                  5 
    ## TACCTCGCACCCAACG-1 TACCTGCAGTCTGCAT-1 TACCTGCTCCAAGCTA-1 TACGCTCTCTCGGTCT-1 
    ##                  3                  1                  2                  0 
    ## TACGCTCTCTTAAGGC-1 TACGGGCTCAACGAGG-1 TACGGGCTCTACCACC-1 TACGTCCCAAGTGATA-1 
    ##                  2                  3                  3                  0 
    ## TACTTCATCATGCCGG-1 TAGACTGAGATTGTGA-1 TAGACTGGTCTTGAAC-1 TAGACTGTCGCCGAGT-1 
    ##                  4                  2                  0                  0 
    ## TAGAGTCGTTGTGTAC-1 TAGATCGCAAATGGTA-1 TAGATCGCACACGCCA-1 TAGATCGTCACAACCA-1 
    ##                  0                  0                  3                  0 
    ## TAGCACAAGTACTCGT-1 TAGGAGGCAAGACCTT-1 TAGGAGGGTTCCAGGC-1 TAGGGTTAGAGAGCGG-1 
    ##                  0                  2                  2                  1 
    ## TAGGGTTAGGAACATT-1 TAGGGTTCAAGAAACT-1 TAGGTACGTGTCATCA-1 TAGGTACTCGAAGCAG-1 
    ##                  2                  2                  0                  6 
    ## TAGTGCACACTCGATA-1 TAGTGCACAGAGGAAA-1 TATACCTGTTCTCGTC-1 TATACCTTCTGTCCGT-1 
    ##                  0                  0                  4                  0 
    ## TATATCCCATGCGTGC-1 TATCGCCTCTCCCAAC-1 TATCTTGTCGGCTTGG-1 TATGTTCAGTATGGAT-1 
    ##                  2                  5                  0                  4 
    ## TATTGCTTCTCGAACA-1 TCAAGCATCCACTTTA-1 TCAATCTAGAGGTTTA-1 TCAATCTGTTTCCCAC-1 
    ##                  6                  1                  3                  2 
    ## TCACAAGGTAGTGGCA-1 TCACAAGGTGGATACG-1 TCACACCCACCAGCCA-1 TCACACCGTCTCACGG-1 
    ##                  3                  0                  0                  4 
    ## TCACACCGTTAAACCC-1 TCACACCTCTGGGAGA-1 TCACACCTCTTACCAT-1 TCACATTCAGAACTAA-1 
    ##                  0                  1                  4                  7 
    ## TCACATTCATGCCATA-1 TCACGCTAGGACATCG-1 TCACGCTGTCTGCCTT-1 TCACGCTTCACTAGCA-1 
    ##                  1                  4                  3                  3 
    ## TCACGCTTCGCTGACG-1 TCACGCTTCGGACTGC-1 TCACGGGCAACCGACC-1 TCACTATGTTTAGACC-1 
    ##                  5                  1                  0                  4 
    ## TCACTCGCAAAGCAAT-1 TCACTCGGTCAACCAT-1 TCACTCGTCGGTCATA-1 TCAGCAACAAATCAAG-1 
    ##                  0                  3                  0                  0 
    ## TCAGCAAGTATGACAA-1 TCAGGGCTCATGCGGC-1 TCAGGTACAATCGCAT-1 TCAGGTAGTATGGTAA-1 
    ##                  5                  3                  2                  1 
    ## TCAGTCCGTAAGATCA-1 TCAGTCCGTCTCTCAC-1 TCAGTTTAGGATATAC-1 TCATACTCAGGTACGA-1 
    ##                  0                  1                  0                  4 
    ## TCATACTTCAGAGCAG-1 TCATATCCATGCTGCG-1 TCATCATAGCACACAG-1 TCATCATGTTGCACGC-1 
    ##                  1                  5                  0                  0 
    ## TCATCCGTCCAACTGA-1 TCATGCCCAAAGCGTG-1 TCATGCCGTGTCTTAG-1 TCATGTTAGAGGCGGA-1 
    ##                  2                  4                  3                  2 
    ## TCATGTTGTGTGTGTT-1 TCATTACGTGCATGAG-1 TCATTGTGTGAGCAGT-1 TCCATCGGTTTAGTCG-1 
    ##                  3                  4                  3                  2 
    ## TCCCACACAGACCAAG-1 TCCCAGTAGACCGTTT-1 TCCCAGTTCGACATAC-1 TCCCATGTCGGCCCAA-1 
    ##                  0                  0                  0                  3 
    ## TCCGATCTCAAACCTG-1 TCCTCCCAGTCTAACC-1 TCCTCGATCTCTCGAC-1 TCCTCTTAGAGTGAAG-1 
    ##                  1                  0                  1                  0 
    ## TCCTCTTAGCCAAGGT-1 TCCTTCTCAACTCGTA-1 TCCTTTCTCGCTGTCT-1 TCGAAGTGTATCAGGG-1 
    ##                  0                  5                  5                  4 
    ## TCGACCTCAAAGACTA-1 TCGACCTCAAGTAGTA-1 TCGACGGTCGAAGTGG-1 TCGATTTGTACGCTTA-1 
    ##                  3                  2                  1                  1 
    ## TCGCAGGAGATTGACA-1 TCGCAGGGTAGGCTCC-1 TCGCAGGGTTCTAACG-1 TCGCTCACAAGATTGA-1 
    ##                  1                  3                  3                  3 
    ## TCGCTTGAGAGACAAG-1 TCGCTTGAGCATTTGC-1 TCGGGACGTAACGGTG-1 TCGGGACTCTAGATCG-1 
    ##                  4                  4                  0                  3 
    ## TCGGGCACAGTTACCA-1 TCGGGCATCCGATTAG-1 TCGGGCATCTGAACGT-1 TCGGGTGCAGTCTGGC-1 
    ##                  3                  1                  0                  0 
    ## TCGGTCTTCATTGCTT-1 TCGGTCTTCGACATAC-1 TCGTAGACAATGAAAC-1 TCGTGCTGTAAGGAGA-1 
    ##                  1                  1                  0                  3 
    ## TCGTGCTGTGATCATC-1 TCTACCGCAGTCGTTA-1 TCTACCGTCCACGGAC-1 TCTATACGTTATCCAG-1 
    ##                  5                  1                  3                  0 
    ## TCTATCAGTATGCTAC-1 TCTCACGAGGCTCCCA-1 TCTCCGAAGTCAATCC-1 TCTCCGATCTCTCTTC-1 
    ##                  0                  4                  0                  2 
    ## TCTCTGGCAAAGCAAT-1 TCTGCCAGTAACGCGA-1 TCTGCCAGTAGATGTA-1 TCTGCCATCCAGGACC-1 
    ##                  7                  3                  0                  5 
    ## TCTGGCTCATCCTGTC-1 TCTGTCGAGACTAGAT-1 TCTGTCGCACCAGCGT-1 TCTTAGTAGGCCTGCT-1 
    ##                  1                  0                  6                  3 
    ## TCTTAGTCAACACGAG-1 TCTTAGTGTTGCCAAT-1 TCTTAGTTCAGAATAG-1 TCTTGCGAGATGAAGG-1 
    ##                  0                  0                  4                  1 
    ## TCTTGCGCACTACTTT-1 TCTTGCGCAGAGGAAA-1 TCTTGCGGTCCAGCCA-1 TCTTTGAGTGGAGAAA-1 
    ##                  2                  1                  1                  2 
    ## TGAACGTCAAATGAAC-1 TGAACGTGTAACTAAG-1 TGAATGCGTGGACTAG-1 TGAATGCTCGACGAGA-1 
    ##                  4                  1                  2                  0 
    ## TGAATGCTCGCCGTGA-1 TGACAGTAGAGGTGCT-1 TGACAGTTCGAGATGG-1 TGACCCTTCTCATTGT-1 
    ##                  0                  3                  2                  0 
    ## TGACGCGGTGGCTTAT-1 TGACTCCAGTCTCCTC-1 TGACTCCCACCATAAC-1 TGAGACTAGGCTGAAC-1 
    ##                  0                  3                  0                  3 
    ## TGAGCATTCCATTGCC-1 TGAGCATTCCTTCACG-1 TGAGCATTCGATACTG-1 TGAGCATTCTGTCTCG-1 
    ##                  4                  1                  1                  0 
    ## TGAGCGCCACTATGTG-1 TGAGCGCCAGCATACT-1 TGAGCGCTCGTCTAAG-1 TGAGGAGTCCATCGTC-1 
    ##                  0                  0                  4                  2 
    ## TGAGGGATCAAACGTC-1 TGAGGTTCAAACTAGA-1 TGAGGTTCAGGACATG-1 TGAGGTTGTAGCTTTG-1 
    ##                  2                  0                  0                  0 
    ## TGAGTCACAGAACGCA-1 TGAGTCAGTAACCCTA-1 TGAGTCAGTCACTCTC-1 TGATCAGAGGATTCAA-1 
    ##                  0                  6                  0                  0 
    ## TGATGCACAGGGCTTC-1 TGATGGTAGACCCGCT-1 TGCACGGTCGAGTACT-1 TGCAGGCAGCGTTCAT-1 
    ##                  4                  6                  2                  1 
    ## TGCAGGCCAGGAATAT-1 TGCAGTAAGTAGTCCT-1 TGCAGTATCCACATAG-1 TGCATCCCACCCAACG-1 
    ##                  3                  2                  0                  2 
    ## TGCATGAAGCACCCAC-1 TGCGACGAGCGACAGT-1 TGCGATAAGAATCGCG-1 TGCGGCAAGTGGTGAC-1 
    ##                  3                  6                  4                  1 
    ## TGCGGCACAAAGACGC-1 TGCTCCACAATCTAGC-1 TGCTCGTAGACGAGCT-1 TGCTCGTTCCCATTTA-1 
    ##                  1                  4                  3                  2 
    ## TGCTGAACAGCACAAG-1 TGCTTCGCATAGCTGT-1 TGCTTGCGTGAACTAA-1 TGGAACTAGCGTCTCG-1 
    ##                  1                  0                  3                  0 
    ## TGGAGAGCAGTCAGTT-1 TGGAGGAAGTAGGAAG-1 TGGAGGAGTCTTTATC-1 TGGATCAGTTCTCCCA-1 
    ##                  0                  0                  5                  0 
    ## TGGATGTTCACGGGAA-1 TGGGAAGTCGATTGGT-1 TGGGAGAAGCAAGTGC-1 TGGGAGAAGTGGTGGT-1 
    ##                  4                  1                  2                  0 
    ## TGGGATTAGTCACGAG-1 TGGGATTGTACGATGG-1 TGGGCGTAGTGGAATT-1 TGGGCTGAGAGGGTCT-1 
    ##                  6                  0                  2                  1 
    ## TGGGCTGAGCGTTAGG-1 TGGGCTGCATCCGATA-1 TGGGCTGGTCGCGGTT-1 TGGGTTAAGGAACGAA-1 
    ##                  7                  2                  6                  0 
    ## TGGGTTAGTATTGAGA-1 TGGGTTAGTCCTTTGC-1 TGGGTTAGTGCACGCT-1 TGGGTTATCAGTAGGG-1 
    ##                  7                  5                  2                  0 
    ## TGGTACACAAAGGGTC-1 TGGTAGTCAGGCATGA-1 TGGTAGTCATCACAGT-1 TGGTAGTTCCCTCATG-1 
    ##                  3                  3                  3                  0 
    ## TGTAACGGTTGCGTAT-1 TGTACAGAGACTCGAG-1 TGTACAGAGTGAGTTA-1 TGTACAGGTGGATTTC-1 
    ##                  5                  0                  0                  3 
    ## TGTAGACTCCATTTAC-1 TGTCAGACAACTTGGT-1 TGTCAGACACTCACTC-1 TGTCCACCAACTGGTT-1 
    ##                  0                  1                  1                  4 
    ## TGTCCACTCGCTTACC-1 TGTCCCAAGCTCCATA-1 TGTCCCAAGGTCATAA-1 TGTCCCACATACTTTC-1 
    ##                  2                  1                  2                  3 
    ## TGTCCCAGTCAGGCAA-1 TGTCCCAGTCGGTGTC-1 TGTCCTGGTACTTGTG-1 TGTCCTGGTTGGCTAT-1 
    ##                  0                  0                  3                  2 
    ## TGTGATGAGCGACTGA-1 TGTGATGCATTCTCTA-1 TGTGCGGGTCGTCGGT-1 TGTGCGGTCACATACG-1 
    ##                  2                  1                  6                  0 
    ## TGTGGCGTCGTACCTC-1 TGTGTGATCTGACAGT-1 TGTTCATAGCTGCCTG-1 TGTTCATGTCCTGGGT-1 
    ##                  0                  3                  1                  4 
    ## TGTTCCGCAAGTATAG-1 TGTTGAGCACAACGAG-1 TGTTGAGTCTCGACCT-1 TGTTTGTGTGGAAATT-1 
    ##                  2                  0                  2                  0 
    ## TTAATCCAGCTACTAC-1 TTACAGGGTGATTCTG-1 TTACCGCGTCCGAAAG-1 TTACGCCCACGTCGTG-1 
    ##                  0                  4                  7                  0 
    ## TTACGCCTCAATGCAC-1 TTACGTTGTTTACGAC-1 TTACGTTTCTCGCTTG-1 TTACTGTCATACATCG-1 
    ##                  3                  1                  1                  2 
    ## TTACTGTGTAAGAACT-1 TTACTGTGTTAAGGGC-1 TTAGGGTGTCACATTG-1 TTAGTCTCACGTATAC-1 
    ##                  0                  2                  4                  3 
    ## TTATTGCCAAACGGCA-1 TTATTGCGTACGTAGG-1 TTATTGCGTCAACACT-1 TTATTGCGTGTCTTCC-1 
    ##                  2                  0                  6                  0 
    ## TTCACCGCAGCCCAGT-1 TTCACGCGTTAAGTCC-1 TTCACGCTCTATCCAT-1 TTCAGGAAGATTGGGC-1 
    ##                  2                  1                  1                  3 
    ## TTCAGGAGTCTACAAC-1 TTCATGTAGGTCATAA-1 TTCCACGGTCAACACT-1 TTCCGTGCACTCACTC-1 
    ##                  2                  1                  2                  3 
    ## TTCCGTGGTTCCAGGC-1 TTCCTAAGTGACGTCC-1 TTCCTCTGTACCTGTA-1 TTCCTTCAGGTAGTCG-1 
    ##                  2                  5                  0                  0 
    ## TTCCTTCAGTAGTCAA-1 TTCCTTCGTGACACGA-1 TTCGATTTCTGACCCT-1 TTCGCTGAGGGTTAGC-1 
    ##                  1                  0                  5                  4 
    ## TTCGCTGAGGTTCCGC-1 TTCGCTGCACTAGTAC-1 TTCGCTGTCGGACTGC-1 TTCTCTCGTATCGTAC-1 
    ##                  1                  2                  3                  1 
    ## TTCTTGAGTCTTGCTC-1 TTGACCCGTCTCGGGT-1 TTGACCCTCGTAGTGT-1 TTGATGGCATAAGATG-1 
    ##                  0                  4                  1                  2 
    ## TTGCATTGTGCCTTCT-1 TTGCCTGAGAGTGACC-1 TTGCCTGAGTGGAAAG-1 TTGCCTGGTAGCTGCC-1 
    ##                  5                  0                  0                  0 
    ## TTGCGTCCACGCCACA-1 TTGCTGCCATTGCCGG-1 TTGGGATCAACGGGTA-1 TTGGGATTCAGGAAGC-1 
    ##                  5                  7                  4                  0 
    ## TTGGGATTCTCTCTAA-1 TTGGGCGCACGGTGAA-1 TTGGGCGGTCGGAAAC-1 TTGGGTAGTGCTAGCC-1 
    ##                  0                  1                  2                  0 
    ## TTGGGTATCACCGACG-1 TTGGTTTCACTGGATT-1 TTGTGGATCTAAGAAG-1 TTGTGTTGTGTGTCCG-1 
    ##                  0                  4                  3                  3 
    ## TTGTTCACACTTGTGA-1 TTGTTCACAGTCGCTG-1 TTGTTCATCTTTACAC-1 TTGTTTGAGGTCGACA-1 
    ##                  5                  0                  0                  3 
    ## TTTACTGTCACGGGAA-1 TTTATGCGTAACATAG-1 TTTATGCGTTGATCGT-1 TTTCACATCGACGCTG-1 
    ##                  2                  1                  2                  5 
    ## TTTCACATCTCAGGCG-1 TTTCATGGTGCCTAAT-1 TTTCATGTCACTCACC-1 TTTCCTCCACAGAGCA-1 
    ##                  2                  2                  0                  2 
    ## TTTCCTCTCCTACACC-1 TTTCCTCTCTCTTGCG-1 TTTGATCTCTTTGGAG-1 TTTGGTTAGTAACCTC-1 
    ##                  1                  7                  1                  1 
    ## TTTGGTTGTAGAATAC-1 TTTGTTGCAATTAGGA-1 
    ##                  0                  2 
    ## Levels: 0 1 2 3 4 5 6 7

``` r
### ---------- NON_LINEAR dimensionality reduction ----------
pbmc1k <- RunUMAP(pbmc1k, dims = 1:10, label = TRUE)
```

    ## Warning: The following arguments are not used: label

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## 13:23:50 UMAP embedding parameters a = 0.9922 b = 1.112

    ## 13:23:50 Read 1134 rows and found 10 numeric columns

    ## 13:23:50 Using Annoy for neighbor search, n_neighbors = 30

    ## 13:23:50 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 13:23:50 Writing NN index file to temp file /var/folders/_6/hb24dn5x3sj71c6w0zmnlx640000gn/T//RtmpQAJ6ER/file1739755816c4
    ## 13:23:50 Searching Annoy index using 1 thread, search_k = 3000
    ## 13:23:50 Annoy recall = 100%
    ## 13:24:05 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 30
    ## 13:24:13 Initializing from normalized Laplacian + noise (using irlba)
    ## 13:24:13 Commencing optimization for 500 epochs, with 41998 positive edges
    ## 13:24:17 Optimization finished

``` r
DimPlot(pbmc1k, reduction = "umap")
```

![](Practice-dataset_files/figure-gfm/QC%20continued...-1.png)<!-- -->

``` r
saveRDS(pbmc1k, file = "pbmc1k.rds")
```
