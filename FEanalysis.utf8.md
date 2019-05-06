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



# Functional analysis

Here we do the functional analysis.

## Session information


```r
sessionInfo()
```

```
R version 3.5.3 (2019-03-11)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.2 LTS

Matrix products: default
BLAS: /usr/local/lib/R/lib/libRblas.so
LAPACK: /usr/local/lib/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=es_ES.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=es_ES.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=es_ES.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=es_ES.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] sva_3.30.1                  genefilter_1.64.0          
 [3] mgcv_1.8-27                 nlme_3.1-137               
 [5] geneplotter_1.60.0          annotate_1.60.1            
 [7] XML_3.98-1.19               AnnotationDbi_1.44.0       
 [9] lattice_0.20-38             edgeR_3.24.3               
[11] limma_3.38.3                SummarizedExperiment_1.12.0
[13] DelayedArray_0.8.0          BiocParallel_1.16.6        
[15] matrixStats_0.54.0          Biobase_2.42.0             
[17] GenomicRanges_1.34.0        GenomeInfoDb_1.18.2        
[19] IRanges_2.16.0              S4Vectors_0.20.1           
[21] BiocGenerics_0.28.0         knitr_1.22                 
[23] BiocStyle_2.10.0           

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.1             RColorBrewer_1.1-2     compiler_3.5.3        
 [4] BiocManager_1.30.4     XVector_0.22.0         bitops_1.0-6          
 [7] tools_3.5.3            zlibbioc_1.28.0        bit_1.1-14            
[10] digest_0.6.18          memoise_1.1.0          evaluate_0.13         
[13] RSQLite_2.1.1          Matrix_1.2-15          DBI_1.0.0             
[16] yaml_2.2.0             xfun_0.6               GenomeInfoDbData_1.2.0
[19] stringr_1.4.0          bit64_0.9-7            locfit_1.5-9.1        
[22] grid_3.5.3             survival_2.43-3        rmarkdown_1.12        
[25] bookdown_0.9           blob_1.1.1             magrittr_1.5          
[28] splines_3.5.3          htmltools_0.3.6        xtable_1.8-4          
[31] stringi_1.4.3          RCurl_1.95-4.12       
```
