---
output:
  BiocStyle::html_document
---

<!---
Because we split the analysis pipeline in different independent files,
to speed up processing it, here in the setup block we load libraries and
objects that were loaded or produced in the previously processed file,
and which are necessary in this file.
--->



# GSEA

Here we ran the GSEA analysis for the tumor dataset. Same as for the DE Analysis, we decided not to explain again every step of the pipeline, which in this case were all explained in detail in the FE Analysis section of the paired dataset.

Same as for the paired dataset, we decided to restrict our analysis to the KEGG, REACTOME and BIOCARTA genesets. We then built the incidence matrix and we only selected the genesets which are present in both our dataset and the incidence matrix.
Furthermore, we filtered out the genesets that contain less than 5 genes.

After that, we calculated the z-scores for the genesets as well as the adjusted p-values and focus on the ones with lower p-values.


```r
# collect gene sets
data(c2BroadSets)
gsc <- GeneSetCollection(c2BroadSets)
gsc <- gsc[c(grep("^KEGG", names(gsc)),grep("^REACTOME", names(gsc)), grep("^BIOCARTA", names(gsc)))]
gsc <- mapIdentifiers(gsc, AnnoOrEntrezIdentifier(metadata(coadse.tumor)$annotation))
# Know which gene is in which geneset
Im <- incidence(gsc)
# Discard genes that are not in our data
Im <- Im[, colnames(Im) %in% rownames(coadse.tumor)]
# Discard genes that are not in the dataset
coadse.tumor <- coadse.tumor[colnames(Im), ]
dge.tumor <- dge.tumor[colnames(Im),]
# Only select genesets of decent size
Im <- Im[rowSums(Im) >= 5, ]


# Z SCORES
# calc moderated test statistics
tGSgenes <- tt[match(colnames(Im), rownames(tt)), "t"]
# calc the z test statistic
zS <- sqrt(rowSums(Im)) * (as.vector(Im %*% tGSgenes)/rowSums(Im))
#Look at the first few geneset
rnkGS <- sort(abs(zS), decreasing = TRUE)
pv <- pmin(pnorm(zS), 1 - pnorm(zS))
sum(pv < 0.5)
```

```
[1] 812
```

```r
pvadj <- p.adjust(pv, method = "fdr")
DEgs <- names(pvadj)[which(pvadj < 0.1)]
length(DEgs)
```

```
[1] 22
```

```r
head(DEgs, n = 3)
```

```
[1] "KEGG_OXIDATIVE_PHOSPHORYLATION" "REACTOME_AXON_GUIDANCE"        
[3] "REACTOME_DIABETES_PATHWAYS"    
```



## Session Information

```r
sessionInfo()
```

```
R version 3.5.3 (2019-03-11)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Mojave 10.14.3

Matrix products: default
BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] it_IT.UTF-8/it_IT.UTF-8/it_IT.UTF-8/C/it_IT.UTF-8/it_IT.UTF-8

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] GOstats_2.48.0              EnhancedVolcano_1.0.1      
 [3] ggrepel_0.8.1               ggplot2_3.1.1              
 [5] calibrate_1.7.2             MASS_7.3-51.4              
 [7] sva_3.30.1                  genefilter_1.64.0          
 [9] mgcv_1.8-28                 nlme_3.1-139               
[11] Category_2.48.1             Matrix_1.2-17              
[13] GSVA_1.30.0                 GSVAdata_1.18.0            
[15] hgu95a.db_3.2.3             GSEABase_1.44.0            
[17] graph_1.60.0                org.Hs.eg.db_3.7.0         
[19] xtable_1.8-4                geneplotter_1.60.0         
[21] annotate_1.60.1             XML_3.98-1.19              
[23] AnnotationDbi_1.44.0        lattice_0.20-38            
[25] edgeR_3.24.3                limma_3.38.3               
[27] SummarizedExperiment_1.12.0 DelayedArray_0.8.0         
[29] BiocParallel_1.16.6         matrixStats_0.54.0         
[31] Biobase_2.42.0              GenomicRanges_1.34.0       
[33] GenomeInfoDb_1.18.2         IRanges_2.16.0             
[35] S4Vectors_0.20.1            BiocGenerics_0.28.0        
[37] knitr_1.22                  BiocStyle_2.10.0           

loaded via a namespace (and not attached):
 [1] bitops_1.0-6           bit64_0.9-7            RColorBrewer_1.1-2    
 [4] Rgraphviz_2.26.0       tools_3.5.3            R6_2.4.0              
 [7] DBI_1.0.0              lazyeval_0.2.2         colorspace_1.4-1      
[10] withr_2.1.2            tidyselect_0.2.5       bit_1.1-14            
[13] compiler_3.5.3         bookdown_0.9           scales_1.0.0          
[16] RBGL_1.58.2            stringr_1.4.0          digest_0.6.18         
[19] rmarkdown_1.12         AnnotationForge_1.24.0 XVector_0.22.0        
[22] pkgconfig_2.0.2        htmltools_0.3.6        rlang_0.3.4           
[25] RSQLite_2.1.1          shiny_1.3.2            dplyr_0.8.1           
[28] RCurl_1.95-4.12        magrittr_1.5           GO.db_3.7.0           
[31] GenomeInfoDbData_1.2.0 Rcpp_1.0.1             munsell_0.5.0         
[34] stringi_1.4.3          yaml_2.2.0             zlibbioc_1.28.0       
[37] plyr_1.8.4             grid_3.5.3             blob_1.1.1            
[40] promises_1.0.1         crayon_1.3.4           splines_3.5.3         
[43] locfit_1.5-9.1         pillar_1.3.1           codetools_0.2-16      
[46] glue_1.3.1             evaluate_0.13          BiocManager_1.30.4    
[49] httpuv_1.5.1           gtable_0.3.0           purrr_0.3.2           
[52] assertthat_0.2.1       xfun_0.6               mime_0.6              
[55] later_0.8.0            survival_2.44-1.1      tibble_2.1.1          
[58] shinythemes_1.1.2      memoise_1.1.0         
```
