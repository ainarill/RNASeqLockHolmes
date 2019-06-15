---
title: "QAanalysis_tumor"
output:
  BiocStyle::html_document
---

<!---
The following chunk of code, which should not be shown in the resulting document (echo=FALSE)
sets up global processing options, such as forcing 'knitr' to stop when an error
in the R code is encountered, caching of the results in the 'cache'
directory and asking 'knitr' to figure out automatically the dependencies among
code chunks to re-calculate cached results (autodep=TRUE).

Other options could be changing the name of the directory where figures end up
('figure' by default), etc. For a full account of 'knitr' options please consult
http://yihui.name/knitr/options

At the end of the chunk a 'cat()' call is made to dump a CSS file that gives
a better look-and-feel than the knitr default one. See the source css/ieo.css
and the resulting projectTemplate.html to understand where this is being dumpted.
--->




# Quality assessment

The Quality assessment for tumor data will be performed following the steps and the reasoning used in the Quality assesment for the paired data. For this reason, here we are only going to highlight the strategies that differ from the previous pipeline.

## Subset of samples: Tumor data

From the total of samples before starting to analyse the data, we proceed to do a subset that only includes the tumor data. The reason of doing this is because we are interested also in the differential expression analysis between the different stages of the Colon Adenocarcinoma Cancer. 

In order to create the subset, we looked for those samples that had "tumor" as their type and discarded those with "normal" type. We also discarded those samples that had no information about the stage of the tumor.

Although the data is categorized in diverse subgroups, we will work with the following ones: stage I, stage II, stage III and stage IV, being each one the sum of all the subcategories that share the same stage number:
 

```r
filter_tumor_mask <- rownames(colData(coadse)[colData(coadse)$type == "tumor",])
coadse.tumor <- coadse[, filter_tumor_mask]
coadse.tumor <- coadse.tumor[, !is.na(coadse.tumor$ajcc_pathologic_tumor_stage)]
coadse.tumor <- coadse.tumor[, coadse.tumor$ajcc_pathologic_tumor_stage != "[Not Available]"]
coadse.tumor$ajcc_pathologic_tumor_stage <- gsub(x = coadse.tumor$ajcc_pathologic_tumor_stage, pattern = "A", replacement = "")
coadse.tumor$ajcc_pathologic_tumor_stage <- gsub(x = coadse.tumor$ajcc_pathologic_tumor_stage, pattern = "B", replacement = "")
coadse.tumor$ajcc_pathologic_tumor_stage <- gsub(x = coadse.tumor$ajcc_pathologic_tumor_stage, pattern = "C", replacement = "")
table(coadse.tumor$ajcc_pathologic_tumor_stage)
```

```

  Stage I  Stage II Stage III  Stage IV 
       72       168       124        63 
```

## Sequencing depth among tumor samples

First of all, we want to collect an overview of the data that we are working with. Starting  examining the sequencing depth by plotting the total number of reads mapped to the genome per sample.


```
CM.4748 AA.3696 AA.3837 AA.A02H AA.A010 AA.A00E A6.2679 AA.3672 AA.3524 
    1.2     3.0     4.5     5.1     5.5     5.7     5.9     6.6     6.8 
AA.A01K AA.A004 AA.3947 AA.A01C CA.6717 D5.6536 AA.3976 AA.3506 AA.3841 
    7.3     7.4     7.6     8.0     8.0     8.1     8.3     8.6     8.6 
A6.2682 AA.A03F AA.3667 AA.3664 AA.A029 AA.3519 AA.3812 AA.3949 AA.3956 
    9.0     9.6     9.9    10.2    10.2    11.0    11.0    11.2    11.3 
AA.3842 AA.3848 AA.3875 AA.3930 AA.3955 AA.A01Q A6.3807 AA.3862 AA.3679 
   11.5    11.5    12.1    12.1    12.1    12.2    12.3    12.3    12.5 
A6.A5ZU AA.3684 AA.3542 AA.3850 AA.A00F A6.A56B AA.3560 AA.3680 AA.3715 
   12.6    12.7    12.8    12.8    13.0    13.1    13.2    13.2    13.2 
AA.A01F AA.A01T CM.4747 AA.3681 AZ.4313 AA.3495 AA.3531 AA.3521 AA.A01Z 
   13.4    13.4    13.4    13.5    13.5    13.6    13.6    13.7    13.7 
AA.3861 AA.3502 AA.3970 AA.A00N AA.3950 AA.A024 AA.3666 AA.3530 AA.3710 
   13.9    14.1    14.1    14.2    14.3    14.3    14.5    14.6    14.6 
AA.3517 AA.3821 AA.3972 AA.A00D AA.3971 AA.A00L A6.3808 AA.3544 AA.3548 
   14.7    14.7    14.7    14.7    14.9    14.9    15.0    15.0    15.0 
AA.3845 AA.3975 AA.3852 AA.A01G AZ.4315 AA.3552 AA.3939 AA.3982 AA.3815 
   15.0    15.0    15.1    15.1    15.1    15.2    15.2    15.2    15.3 
AA.A00R A6.A567 AA.3534 AA.3941 AA.3968 AA.3846 AA.3867 AA.3973 AA.A00K 
   15.3    15.4    15.4    15.4    15.4    15.5    15.5    15.6    15.7 
CM.6169 A6.2683 AA.A00J AA.A00U AA.A02E DM.A1DB AA.3811 AA.A00W AA.A03J 
   15.7    15.8    15.8    15.9    15.9    15.9    16.0    16.0    16.0 
AA.3833 AA.A017 AA.A00Q CM.4752 AA.A02F AA.3518 AA.3844 AA.3529 AA.3980 
   16.2    16.2    16.3    16.3    16.4    16.5    16.5    16.7    16.7 
AA.3520 AY.4070 AA.3522 AA.3554 AA.3860 AA.3543 AA.A00O AZ.4308 A6.2676 
   16.8    16.8    17.0    17.0    17.0    17.1    17.1    17.1    17.2 
AA.3516 AA.3989 AA.3858 AA.A00Z AY.4071 AA.3488 AA.3561 AA.3877 AA.A01V 
   17.2    17.2    17.3    17.3    17.3    17.4    17.4    17.4    17.4 
AA.3525 AA.3492 AA.3851 AA.A00A AA.3514 AA.3549 AA.3866 AA.3984 AA.3986 
   17.5    17.6    17.6    17.6    17.7    17.7    17.7    17.7    17.7 
AA.3994 AA.A022 AA.3870 AZ.4684 AA.3864 AA.3979 AA.A01X A6.2680 AA.3538 
   17.7    17.7    17.9    17.9    18.0    18.0    18.0    18.1    18.2 
AA.3556 AA.3695 AA.3814 AA.3509 AA.A01S CM.4746 AA.3869 AA.A01P AA.A02W 
   18.2    18.2    18.3    18.4    18.5    18.5    18.7    18.7    18.7 
AA.A02R CA.5256 G4.6304 AA.3510 AA.3527 AA.3952 AA.A01D AZ.4681 CM.5341 
   18.8    18.8    18.8    19.0    19.0    19.0    19.0    19.0    19.0 
AA.3553 A6.2685 AA.3562 AA.3855 CM.4750 AA.A01R AZ.4614 AA.3831 AA.A02J 
   19.1    19.2    19.2    19.2    19.2    19.4    19.6    19.7    19.7 
AA.A02O AA.3819 CK.4951 A6.4107 AA.3494 AA.3532 AA.3966 AA.3688 NH.A6GC 
   19.8    19.9    19.9    20.0    20.0    20.1    20.6    20.7    20.7 
AA.A01I AZ.4615 A6.2671 A6.2678 AA.3818 A6.2681 AA.3678 AA.3555 AA.3854 
   20.9    21.0    21.1    21.1    21.2    21.4    21.4    21.5    21.6 
4T.AA8H G4.6310 AA.3692 D5.6920 AA.3856 AA.3872 AA.3673 AA.3693 F4.6856 
   21.7    21.8    21.9    21.9    22.4    22.7    23.0    23.3    23.8 
G4.6306 F4.6806 D5.6924 CM.6677 CM.5861 AA.3697 CM.5864 D5.6533 CA.5796 
   23.8    24.3    24.5    25.1    25.4    25.6    25.8    26.3    26.6 
CM.5862 F4.6855 NH.A6GB G4.6323 A6.6141 CM.6166 F4.6808 G4.6320 AY.6196 
   26.8    26.9    27.0    27.4    28.0    28.0    28.2    28.3    29.0 
NH.A8F8 AZ.5407 A6.A566 D5.6929 DM.A288 CM.5860 D5.6927 CM.5349 G4.6309 
   29.0    29.1    29.2    29.3    29.4    29.8    29.8    29.9    29.9 
D5.5540 AZ.6598 CM.4743 F4.6805 CK.5912 G4.6298 D5.6930 G4.6315 A6.2675 
   30.0    30.2    30.2    30.3    30.4    30.6    30.8    30.8    31.0 
AA.3663 A6.6142 CM.6165 G4.6311 G4.6295 F4.6460 F4.6461 CM.6674 CK.5916 
   31.1    31.3    31.5    31.5    31.6    31.8    31.8    31.9    32.3 
G4.6321 AA.3662 DM.A280 G4.6293 CM.6167 G4.6588 AU.3779 CA.5255 QL.A97D 
   32.4    32.5    32.5    32.5    32.6    32.6    32.7    32.9    33.0 
A6.4105 AD.5900 AZ.4323 AZ.6607 D5.6530 F4.6463 AA.3655 D5.5537 F4.6807 
   33.2    33.2    33.2    33.2    33.2    33.2    33.5    33.5    33.5 
AA.3660 A6.5666 G4.6302 D5.6535 AZ.6599 RU.A8FL CK.5915 CK.6751 CM.6172 
   33.7    33.8    33.9    34.0    34.2    34.2    34.3    34.5    34.5 
4N.A93T AY.A8YK G4.6627 3L.AA1B AM.5821 AZ.6605 CM.5344 D5.6531 DM.A28K 
   34.6    34.6    34.6    34.7    34.7    34.7    34.9    34.9    34.9 
D5.5538 AD.6889 G4.6294 G4.6625 A6.6652 A6.6653 D5.6539 G4.6299 AD.A5EJ 
   35.3    35.4    35.4    35.5    35.6    35.6    35.6    35.7    36.3 
D5.6538 NH.A50V AZ.5403 CA.6719 CM.6171 DM.A28G G4.6626 A6.6137 D5.6540 
   36.3    36.3    36.4    36.4    36.5    36.6    36.6    36.7    36.8 
G4.6297 D5.6931 A6.6140 A6.6138 AZ.4616 DM.A28E CM.6675 CM.6680 A6.2686 
   36.9    37.0    37.1    37.3    37.4    37.5    37.6    37.7    37.8 
CM.6161 CM.4751 CM.6676 CK.6746 D5.5539 NH.A50T A6.5667 AM.5820 CM.5863 
   37.9    38.0    38.0    38.1    38.1    38.1    38.2    38.2    38.3 
G4.6586 A6.5660 AY.6386 D5.6926 DM.A282 AY.A54L AZ.6601 A6.5664 CK.4947 
   38.3    38.4    38.4    38.5    38.5    38.6    38.6    38.8    38.8 
AA.3526 CA.6716 A6.6654 A6.6782 DM.A0XD F4.6459 A6.5662 D5.6537 G4.6314 
   39.0    39.0    39.2    39.3    39.4    39.4    39.5    39.5    39.5 
D5.6932 DM.A1D4 G4.6628 DM.A28A AA.3496 AU.6004 D5.6532 AA.3685 CM.6163 
   39.6    39.6    39.6    39.7    39.8    39.8    39.8    39.9    39.9 
F4.6570 CM.6170 D5.6928 DM.A28C CM.5348 NH.A5IV AA.3511 AA.3712 F4.6809 
   39.9    40.0    40.0    40.1    40.2    40.2    40.3    40.4    40.4 
AA.A02K AZ.6606 AZ.6608 CK.4952 QG.A5YV QG.A5YX AY.6197 D5.6529 CM.6168 
   40.6    40.6    40.6    40.6    40.6    40.7    40.8    40.9    41.0 
DM.A28H CM.6162 CK.5914 AZ.6600 AA.3713 A6.6651 CA.5797 CM.6679 CA.6718 
   41.0    41.1    41.6    41.9    42.0    42.3    42.5    42.5    42.7 
CM.6678 AD.6548 AA.A02Y AD.A5EK AA.3675 CK.5913 D5.7000 DM.A1HA CA.6715 
   42.7    42.8    43.2    43.2    43.3    43.3    43.3    43.7    44.1 
D5.6922 QG.A5YW AY.A71X F4.6569 AZ.4682 AA.3489 D5.6898 CM.6164 DM.A1HB 
   44.2    44.3    44.4    44.5    44.6    45.2    45.2    45.3    45.9 
NH.A6GA AY.5543 CM.5868 G4.6322 D5.5541 CK.6748 QG.A5Z2 DM.A285 AD.6888 
   46.1    46.5    46.6    46.6    46.7    46.8    46.8    47.0    47.2 
F4.6854 G4.6303 DM.A28F AD.6899 CA.5254 CK.6747 DM.A1D7 DM.A28M A6.A565 
   47.3    47.9    48.0    48.2    48.2    48.3    48.4    48.5    48.6 
QG.A5Z1 D5.6541 AY.A69D DM.A1D9 WS.AB45 DM.A0XF CM.4744 NH.A50U DM.A1DA 
   48.9    49.0    49.1    49.4    49.4    50.6    50.7    50.7    50.9 
CK.4948 D5.6534 D5.6923 DM.A1D0 CK.4950 A6.6649 F4.6703 DM.A0X9 A6.5657 
   51.5    51.7    52.0    52.0    52.5    52.7    53.9    54.5    56.7 
G4.6307 DM.A1D6 A6.6648 F4.6704 
   56.9    62.7    63.1    67.5 
```

<div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/libsizestumor-1.png" alt="Library sizes in increasing order." width="600px" />
<p class="caption">(\#fig:libsizestumor)Library sizes in increasing order.</p>
</div>

Looking at the plot in figure \@ref(fig:libsizestumor), we can see two important details that we must take into consideration. 
The first one is that all four stages have different samples with different sequencing depth.
The other one is that we observe that within the considerable differences in the sequencing depth of some samples, we can distinguish some samples on the left part of the graph that present considerably low sequencing depths: this may generate problems in further analysis and that is why we then wanted to remove these samples as shown below:



```r
# Remove those samples with very low sequencing depth (<10)
coadse.tumor <- coadse.tumor[, dge.tumor$samples$lib.size >= 10*1e06]
dge.tumor <- DGEList(counts = assays(coadse.tumor)$counts, group = coadse.tumor$ajcc_pathologic_tumor_stage, genes = genes.tumor)
```


```r
# coverage_filtered plot
par(mfrow = c(1, 1), mar = c(4, 5, 1, 1))
ord.tumor <- order(dge.tumor$sample$lib.size)
ggplot(as.data.frame(dge.tumor$samples), aes(x = row.names(dge.tumor$samples), y = dge.tumor$samples$lib.size[ord.tumor]/1000000, fill = dge.tumor$samples$group[ord.tumor])) + geom_bar(stat = "identity") + theme_bw() + theme(legend.title = element_blank(), axis.text.x=element_blank()) + ylab("Milions of reads") + xlab("Samples")
```

<div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/libsizesfilt-1.png" alt="Library sizes in increasing order after the filtering" width="100%" />
<p class="caption">(\#fig:libsizesfilt)Library sizes in increasing order after the filtering</p>
</div>

In this second plot(Figure \@ref(fig:libsizesfilt)) we can see the subset of data that we are going to use after having removed the samples with less coverage depth. In total,we have lost 21 samples from our dataset.

## Gender among samples in tumor data


```r
ord2 <- order(dge.tumor$sample$lib.size)
table(coadse.tumor$gender)
```

```

FEMALE   MALE 
   194    212 
```

```r
barplot(dge.tumor$sample$lib.size[ord2]/1e+06, las = 1, ylab = "Millions of reads", xlab = "Samples",
        col = c("cyan", "orange")[coadse.tumor$gender[ord2]], border = NA)
legend("topleft", c("Female", "Male"), fill = c("cyan", "orange"), inset = 0.01)
```

<div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/barplotgender-1.png" alt="Library sizes with respect to gender" width="100%" />
<p class="caption">(\#fig:barplotgender)Library sizes with respect to gender</p>
</div>

In the figure \@ref(fig:barplotgender) we can see that, as happened with the paired analysis, not only we have a similar number of female and male samples, but their sequencing depth is quite well equilibrated as we cannot clearly identify any defined clustering event.

## Distribution of expression levels among samples in tumor data

Then, we first run a within sample normalization and we also explore the distribution of the expression levels through all the samples in terms of logarithmic CPM units.

<div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/disrRawExpAlTumorl-1.png" alt="Non-parametric density distribution of expression profiles per sample." width="800px" />
<p class="caption">(\#fig:disrRawExpAlTumorl)Non-parametric density distribution of expression profiles per sample.</p>
</div>

In figure \@ref(fig:disrRawExpAllTumor) we can observe a non-parametric density distribution of expression profiles per sample which follows a bimodal distribution.

<div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/disrRawExpCompTumor-1.png" alt="Non-parametric density distribution of expression profiles per sample." width="800px" />
<p class="caption">(\#fig:disrRawExpCompTumor)Non-parametric density distribution of expression profiles per sample.</p>
</div>

To provide a clearer visualization, we decided to divide the dataset into four smaller subsets, one for each stage. Figure \@ref(fig:distRawExpCompTumor) shows the expression levels in terms of logCPM for the different stages of tumor samples separately.

As we do not observe extreme differences in the distribution of the expression values in the different tumor stages, we decide to not exclude any sample from the dataset in this step.

## Distribution of expression levels among genes in tumor data

We now want to look if there is any gene which has very low expression to exclude them from the dataset. 

<div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/exprdisttumor-1.png" alt="Distribution of average expression level per gene." width="400px" />
<p class="caption">(\#fig:exprdisttumor)Distribution of average expression level per gene.</p>
</div>

In order to do so, we have calculated the average expression of each gene for all the samples and plotted their distribution in Figure \@ref(fig:exprdisttumor).

## Filtering of lowly-expressed genes in tumor data

Those genes that have very low counts should be removed prior to further analysis because a gene must be expressed at some minimal level before it is likely to be translated into a protein and also because genes statistically they are very unlikely to be assessed as significantly differential expressed. So, such genes can therefore be removed from the analysis without any loss of information.

Then, to filter the lowly-expressed genes we calculate first a minimum CPM cutoff value of expression; 10 reads per million in the sample with lowest depth.
We will keep also only those genes that have this minimum level of expression in at least as many samples as the smallest group of comparison.


```r
cpmcutoff.tumor <- round(10/min(dge.tumor$sample$lib.size/1e+06), digits = 1)
cpmcutoff.tumor 
```

```
[1] 1
```

```r
sprintf("Cutoff: %s", cpmcutoff.tumor)
```

```
[1] "Cutoff: 1"
```

```r
nsamplescutoff.tumor <- min(table(coadse.tumor$ajcc_pathologic_tumor_stage))
nsamplescutoff.tumor
```

```
[1] 60
```

```r
mask.tumor <- rowSums(cpm(dge.tumor) > cpmcutoff.tumor) >= nsamplescutoff.tumor
coadse.tumor.filt <- coadse.tumor[mask.tumor, ]
dge.tumor.filt <- dge.tumor[mask.tumor, ]
```

After the filtering of lowly-expressed genes, we have lost almost 6.500 genes.


```r
dim(coadse.tumor)
```

```
[1] 20115   406
```

```r
dim(coadse.tumor.filt)
```

```
[1] 13569   406
```
We can visually observe which genes have been left out from the datatset in Figure \@ref(fig:filter)


```r
par(mar = c(4, 5, 1, 1))
h <- hist(avgexp.tumor, xlab = expression("Expression level (" * log[2] * "CPM)"), main = "",
          las = 1, col = "grey", cex.axis = 1.2, cex.lab = 1.5)
x <- cut(rowMeans(assays(coadse.tumor.filt)$logCPM), breaks = h$breaks)
lines(h$mids, table(x), type = "h", lwd = 10, lend = 1, col = "darkred")
legend("topright", c("All genes", "Filtered genes"), fill = c("grey", "darkred"))
```

<div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/histfilt-1.png" alt="Distribution of average expression level per gene after the filtering" width="100%" />
<p class="caption">(\#fig:histfilt)Distribution of average expression level per gene after the filtering</p>
</div>

We can visually observe which genes have been left out from the datatset in Figure \@ref(fig:histfilt).

## Normalization tumor data

To make gene expression values comparable across samples, normalization step is needed to continue with our analysis.
We estimate a normalization factor for each library using the Trimmed Mean of M-Values and we apply it to our data.


```r
dge.tumor.filt <- calcNormFactors(dge.tumor.filt) 
assays(coadse.tumor.filt)$logCPM <- cpm(dge.tumor.filt, log=TRUE, normalized.lib.sizes=TRUE, prior.count=0.25)
```


```r
saveRDS(coadse.tumor.filt, file.path("results", "coadse.tumor.filt.rds"))
saveRDS(dge.tumor.filt, file.path("results", "dge.tumor.filt.rds"))
```


## MA plots tumor data

We now want to visualize the expression profiles of the normalized data for each tumor stage. 

For Stage I --> 
<div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/maPlotsTumor1-1.png" alt="MA-plots of the tumor samples in Stage I." width="600" />
<p class="caption">(\#fig:maPlotsTumor11)MA-plots of the tumor samples in Stage I.</p>
</div><div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/maPlotsTumor1-2.png" alt="MA-plots of the tumor samples in Stage I." width="600" />
<p class="caption">(\#fig:maPlotsTumor12)MA-plots of the tumor samples in Stage I.</p>
</div>

For Stage II --> 
<div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/maPlotsTumor2-1.png" alt="MA-plots of the tumor samples in Stage II." width="600" />
<p class="caption">(\#fig:maPlotsTumor21)MA-plots of the tumor samples in Stage II.</p>
</div><div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/maPlotsTumor2-2.png" alt="MA-plots of the tumor samples in Stage II." width="600" />
<p class="caption">(\#fig:maPlotsTumor22)MA-plots of the tumor samples in Stage II.</p>
</div><div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/maPlotsTumor2-3.png" alt="MA-plots of the tumor samples in Stage II." width="600" />
<p class="caption">(\#fig:maPlotsTumor23)MA-plots of the tumor samples in Stage II.</p>
</div><div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/maPlotsTumor2-4.png" alt="MA-plots of the tumor samples in Stage II." width="600" />
<p class="caption">(\#fig:maPlotsTumor24)MA-plots of the tumor samples in Stage II.</p>
</div>


For Stage III --> 
<div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/maPlotsTumor3-1.png" alt="MA-plots of the tumor samples in Stage III." width="600" />
<p class="caption">(\#fig:maPlotsTumor31)MA-plots of the tumor samples in Stage III.</p>
</div><div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/maPlotsTumor3-2.png" alt="MA-plots of the tumor samples in Stage III." width="600" />
<p class="caption">(\#fig:maPlotsTumor32)MA-plots of the tumor samples in Stage III.</p>
</div><div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/maPlotsTumor3-3.png" alt="MA-plots of the tumor samples in Stage III." width="600" />
<p class="caption">(\#fig:maPlotsTumor33)MA-plots of the tumor samples in Stage III.</p>
</div>

For Stage IV --> 
<div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/maPlotsTumor4-1.png" alt="MA-plots of the tumor samples in Stage IV." width="600" />
<p class="caption">(\#fig:maPlotsTumor41)MA-plots of the tumor samples in Stage IV.</p>
</div><div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/maPlotsTumor4-2.png" alt="MA-plots of the tumor samples in Stage IV." width="600" />
<p class="caption">(\#fig:maPlotsTumor42)MA-plots of the tumor samples in Stage IV.</p>
</div>

## Batch Identification for tumor data.

As normalization of the data can not always remove completely the batch effect, the next step will be indeed the search of potential surrogate of batch effect indicators derived from different elements of the TCGA barcode of each sample.


```r
tss.tumor <- substr(colnames(dge.tumor.filt), 6, 7)
table(data.frame(STAGE=coadse.tumor.filt$ajcc_pathologic_tumor_stage, TSS=tss.tumor))
```

```
           TSS
STAGE       3L 4N 4T A6 AA AD AM AU AY AZ CA CK CM D5 DM F4 G4 NH QG QL RU
  Stage I    1  0  0  2 31  3  0  1  3  3  0  4  9  5  0  4  2  0  1  1  0
  Stage II   0  0  1 12 63  2  1  1  2  5  8  4 12 13 15  6  9  2  1  0  0
  Stage III  0  1  0 15 37  2  0  0  3  4  1  4 10 11  8  5 10  3  3  0  1
  Stage IV   0  0  0  7 24  0  1  0  2  8  0  2  5  1  1  1  5  3  0  0  0
           TSS
STAGE       WS
  Stage I    0
  Stage II   1
  Stage III  0
  Stage IV   0
```

```r
center.tumor <- substr(colnames(dge.tumor.filt), 27, 28)
table(data.frame(STAGE=coadse.tumor.filt$ajcc_pathologic_tumor_stage, CENTER=center.tumor))
```

```
           CENTER
STAGE        07
  Stage I    70
  Stage II  158
  Stage III 118
  Stage IV   60
```

```r
samplevial.tumor <- substr(colnames(dge.tumor.filt), 14, 16)
table(data.frame(STAGE=coadse.tumor.filt$ajcc_pathologic_tumor_stage, SAMPLEVIAL=samplevial.tumor))
```

```
           SAMPLEVIAL
STAGE       01A 01B
  Stage I    70   0
  Stage II  158   0
  Stage III 117   1
  Stage IV   58   2
```

```r
portionanalyte.tumor <- substr(colnames(dge.tumor.filt), 18, 20)
table(data.frame(STAGE=coadse.tumor.filt$ajcc_pathologic_tumor_stage, PORTIONANALYTE=portionanalyte.tumor))
```

```
           PORTIONANALYTE
STAGE       01R 02R 03R 11R 12R 21R 23R 31R 33R 42R 43R 72R
  Stage I    28   6   0  31   2   1   1   0   0   0   1   0
  Stage II   68  12   0  59   4  14   0   0   0   1   0   0
  Stage III  33  14   0  51   0  18   0   2   0   0   0   0
  Stage IV   31   2   1  18   1   3   0   2   1   0   0   1
```

```r
plate.tumor <- substr(colnames(dge.tumor.filt), 22, 25)
table(data.frame(STAGE=coadse.tumor.filt$ajcc_pathologic_tumor_stage, PLATE=plate.tumor))
```

```
           PLATE
STAGE       0821 0826 0905 1022 1113 1410 1653 1723 1774 1839 1873 1928 A002
  Stage I      5    1    3   10    1    4   10    5    4    6    1    2    3
  Stage II    10    3   12   15    1   11   13   14   13   16    0   12    6
  Stage III    6    2    5    9    2    5    9   22    6   12    2    5    5
  Stage IV     1    0    7    9    0    6    6    5    6    8    0    1    0
           PLATE
STAGE       A00A A083 A089 A155 A16W A180 A28H A32Y A32Z A37K A41B
  Stage I      1    1    1    0    1    0    3    1    4    2    1
  Stage II     5    1    0    8    6    0    3    3    2    2    2
  Stage III    1    4    1    3    3    1    8    2    1    4    0
  Stage IV     2    0    1    0    1    0    1    1    1    1    3
```

Too look for potential surrogates of batch effect indicators we need the maximum information explained by the possible batch generator, that means no zeros along the different groups. 
In order to be able to identify a possible batch effect and looking the previous tables, we decided to continue our analysis focusing in 'plate' information.
Then we decided to see graphically how batch might be confounding the outcome of interest by doing a hierarchical clustering.


```r
logCPM <- cpm(dge.tumor.filt, log=TRUE, prior.count=3)
d <- as.dist(1-cor(logCPM, method="spearman"))
sampleClustering <- hclust(d)
batch <- as.integer(factor(plate.tumor))
sampleDendrogram <- as.dendrogram(sampleClustering, hang=0.1)
names(batch) <- colnames(coadse.tumor.filt)
outcome <- paste(substr(colnames(coadse.tumor.filt), 9, 12), as.character(coadse.tumor.filt$ajcc_pathologic_tumor_stage), sep="-")
names(outcome) <- colnames(coadse.tumor.filt)
sampleDendrogram <- dendrapply(sampleDendrogram,
                               function(x, batch, labels) {
                                 if (is.leaf(x)) {
                                   attr(x, "nodePar") <- list(lab.col=as.vector(batch[attr(x, "label")]))
                                   attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                 }
                                 x
                               }, batch, outcome)
plot(sampleDendrogram, main="Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch)), levels(factor(plate.tumor))), cex=0.75,  fill=sort(unique(batch)))
```

<div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/sampleClustering-1.png" alt="Hierarchical clustering of the samples." width="1400" />
<p class="caption">(\#fig:sampleClustering)Hierarchical clustering of the samples.</p>
</div>

As a result of this graph in figure \@ref(fig:sampleClustering), we concluded that the hierarchical clustering can become quite useless as a graphical tool when you have a huge quantity of samples, so then we decided to change our approach and try to identify graphically the batch effect with a MDS plot.


```r
col.stage <- c("blue","cyan", "orange", "purple")[factor(coadse.tumor.filt$ajcc_pathologic_tumor_stage)]
batch <- as.integer(factor(plate.tumor))
plotMDS(dge.tumor.filt, dim=c(3,4), col=col.stage, pch= batch, cex=2)
legend("topleft", legend=levels(factor(coadse.tumor.filt$ajcc_pathologic_tumor_stage)), col=col.stage ,pch=16, cex= 0.75)
legend("topright", legend=levels(factor(plate.tumor)),pch= batch, cex=0.75, ncol = 2)
```

<div class="figure" style="text-align: center">
<img src="QAanalysisTumor_files/figure-html/mds-1.png" alt=" MDS plot" width="100%" />
<p class="caption">(\#fig:mds) MDS plot</p>
</div>

With the MDSplot in figure \@ref(fig:mds) we are still not able to recognize any cluster of any group, so we can not conclude that 'plate' information could be a potential surrogate of batch effect.

## Session information


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
 [1] geneplotter_1.60.0          annotate_1.60.1            
 [3] XML_3.98-1.19               AnnotationDbi_1.44.0       
 [5] lattice_0.20-38             ggplot2_3.1.1              
 [7] edgeR_3.24.3                limma_3.38.3               
 [9] SummarizedExperiment_1.12.0 DelayedArray_0.8.0         
[11] BiocParallel_1.16.6         matrixStats_0.54.0         
[13] Biobase_2.42.0              GenomicRanges_1.34.0       
[15] GenomeInfoDb_1.18.2         IRanges_2.16.0             
[17] S4Vectors_0.20.1            BiocGenerics_0.28.0        
[19] knitr_1.22                  BiocStyle_2.10.0           

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.1             locfit_1.5-9.1         assertthat_0.2.1      
 [4] digest_0.6.18          R6_2.4.0               plyr_1.8.4            
 [7] RSQLite_2.1.1          evaluate_0.13          highr_0.8             
[10] pillar_1.3.1           zlibbioc_1.28.0        rlang_0.3.4           
[13] lazyeval_0.2.2         blob_1.1.1             Matrix_1.2-17         
[16] rmarkdown_1.12         labeling_0.3           stringr_1.4.0         
[19] RCurl_1.95-4.12        bit_1.1-14             munsell_0.5.0         
[22] compiler_3.5.3         xfun_0.6               pkgconfig_2.0.2       
[25] htmltools_0.3.6        tidyselect_0.2.5       tibble_2.1.1          
[28] GenomeInfoDbData_1.2.0 bookdown_0.9           codetools_0.2-16      
[31] crayon_1.3.4           dplyr_0.8.1            withr_2.1.2           
[34] bitops_1.0-6           grid_3.5.3             xtable_1.8-4          
[37] gtable_0.3.0           DBI_1.0.0              magrittr_1.5          
[40] scales_1.0.0           KernSmooth_2.23-15     stringi_1.4.3         
[43] XVector_0.22.0         RColorBrewer_1.1-2     tools_3.5.3           
[46] bit64_0.9-7            glue_1.3.1             purrr_0.3.2           
[49] yaml_2.2.0             colorspace_1.4-1       BiocManager_1.30.4    
[52] memoise_1.1.0         
```
