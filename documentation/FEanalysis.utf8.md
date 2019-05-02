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
R version 3.5.0 (2018-04-23)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS: /opt/R/R-3.5.0/lib64/R/lib/libRblas.so
LAPACK: /opt/R/R-3.5.0/lib64/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF8       LC_NUMERIC=C             
 [3] LC_TIME=en_US.UTF8        LC_COLLATE=en_US.UTF8    
 [5] LC_MONETARY=en_US.UTF8    LC_MESSAGES=en_US.UTF8   
 [7] LC_PAPER=en_US.UTF8       LC_NAME=C                
 [9] LC_ADDRESS=C              LC_TELEPHONE=C           
[11] LC_MEASUREMENT=en_US.UTF8 LC_IDENTIFICATION=C      

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] sva_3.28.0                  genefilter_1.62.0          
 [3] mgcv_1.8-24                 nlme_3.1-137               
 [5] geneplotter_1.58.0          annotate_1.58.0            
 [7] XML_3.98-1.16               AnnotationDbi_1.42.1       
 [9] lattice_0.20-35             edgeR_3.22.5               
[11] limma_3.36.5                SummarizedExperiment_1.10.1
[13] DelayedArray_0.6.6          BiocParallel_1.14.2        
[15] matrixStats_0.54.0          Biobase_2.40.0             
[17] GenomicRanges_1.32.7        GenomeInfoDb_1.16.0        
[19] IRanges_2.14.12             S4Vectors_0.18.3           
[21] BiocGenerics_0.26.0         knitr_1.20                 
[23] BiocStyle_2.8.2            

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.19           RColorBrewer_1.1-2     compiler_3.5.0        
 [4] XVector_0.20.0         bitops_1.0-6           tools_3.5.0           
 [7] zlibbioc_1.26.0        bit_1.1-14             digest_0.6.18         
[10] memoise_1.1.0          evaluate_0.12          RSQLite_2.1.1         
[13] Matrix_1.2-14          DBI_1.0.0              yaml_2.2.0            
[16] xfun_0.3               GenomeInfoDbData_1.1.0 stringr_1.3.1         
[19] bit64_0.9-7            locfit_1.5-9.1         rprojroot_1.3-2       
[22] grid_3.5.0             survival_2.42-6        rmarkdown_1.10        
[25] bookdown_0.7           blob_1.1.1             magrittr_1.5          
[28] splines_3.5.0          backports_1.1.2        htmltools_0.3.6       
[31] xtable_1.8-3           stringi_1.2.4          RCurl_1.95-4.11       
```
