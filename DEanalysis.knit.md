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





# Differential expression


The following step is to perform a simple examination of expression changes and their associated p-values using the R/Bioconductor package [sva](http://bioconductor.org/packages/sva) as follows:


```r
library(sva)
tss <- substr(colnames(dge.filt), 6, 7)
mod <- model.matrix(~ coadse.filt$type + tss, colData(coadse.filt))
mod0 <- model.matrix(~tss, colData(coadse.filt))
pv <- f.pvalue(assays(coadse.filt)$logCPM, mod, mod0)
sum(p.adjust(pv, method="fdr") < 0.01)
```

```
[1] 7537
```

There are 7537 genes changing significantly
their expression at FDR < 1%. In Figure \@ref(fig:pdist) below we show the distribution of the
resulting p-values.

<div class="figure" style="text-align: center">
<img src="DEanalysis_files/figure-html/pdist-1.png" alt="Distribution of raw p-values for an F-test on every gene between tumor and normal samples." width="400px" />
<p class="caption">(\#fig:pdist)Distribution of raw p-values for an F-test on every gene between tumor and normal samples.</p>
</div>


Then we are going to estimate surrogate variables using the `sva()` function.


```r
sv <- sva(assays(coadse.filt)$logCPM, mod, mod0)
```

```
Number of significant surrogate variables is:  16 
Iteration (out of 5 ):1  2  3  4  5  
```

```r
sv$n
```

```
[1] 16
```

The SVA algorithm has found 16 surrogate variables. So, we are going to use them to
assess againt the extent of differential expression this time adjusting for these
surrogate variables.


```r
modsv <- cbind(mod, sv$sv)
mod0sv <- cbind(mod0, sv$sv)
pvsv <- f.pvalue(assays(coadse.filt)$logCPM, modsv, mod0sv)
sum(p.adjust(pvsv, method="fdr") < 0.01)
```

```
[1] 9796
```

We have increased the number of changing genes to 9796.
Figure \@ref(fig:psvdist) shows the resulting distribution of p-values.

<div class="figure" style="text-align: center">
<img src="DEanalysis_files/figure-html/psvdist-1.png" alt="Distribution of raw p-values for an F-test on every gene between tumor and normal samples, adjusting for surrogate variables estimated with SVA." width="400px" />
<p class="caption">(\#fig:psvdist)Distribution of raw p-values for an F-test on every gene between tumor and normal samples, adjusting for surrogate variables estimated with SVA.</p>
</div>

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
 [1] locfit_1.5-9.1         xfun_0.6               splines_3.5.3         
 [4] htmltools_0.3.6        yaml_2.2.0             blob_1.1.1            
 [7] survival_2.43-3        DBI_1.0.0              bit64_0.9-7           
[10] RColorBrewer_1.1-2     GenomeInfoDbData_1.2.0 stringr_1.4.0         
[13] zlibbioc_1.28.0        codetools_0.2-16       evaluate_0.13         
[16] memoise_1.1.0          highr_0.8              Rcpp_1.0.1            
[19] xtable_1.8-4           BiocManager_1.30.4     XVector_0.22.0        
[22] bit_1.1-14             digest_0.6.18          stringi_1.4.3         
[25] bookdown_0.9           grid_3.5.3             tools_3.5.3           
[28] bitops_1.0-6           magrittr_1.5           RCurl_1.95-4.12       
[31] RSQLite_2.1.1          Matrix_1.2-15          rmarkdown_1.12        
[34] compiler_3.5.3        
```
