# GET THE PACKAGES
suppressPackageStartupMessages({
    library(edgeR)
    library(SummarizedExperiment)
    library(geneplotter)
    library(BiocStyle)
  }
)
# ------------------------------------------------------------------------------------------------------------

# GET THE DATA 
data <- readRDS("seCOAD.rds")
lclse<-data
dge <- DGEList(counts = assays(lclse)$counts, genes = as.data.frame(rowData(lclse)))
mask_remove_na_tumor <-!is.na(lclse$tumor_status)
mask_remove_na_tumor
lclse.test <- lclse[,mask_remove_na_tumor]

# ------------------------------------------------------------------------------------------------------------

# Having an overview 
ord <- order(dge$sample$lib.size)
barplot(dge$sample$lib.size[ord]/1e+06, las = 1, ylab = "Millions of reads", xlab = "Samples",
        col = c("cyan", "orange")[lclse$gender[ord]], border=NA )
legend("topleft", c("Female", "Male"), fill = c("cyan", "orange"), inset = 0.01)

#CPM scaling - within sample normalization
CPM <- t(t(dge$counts)/(dge$samples$lib.size/1e+06))
assays(lclse)$logCPM <- cpm(dge, log = TRUE, prior.count = 0.25) 
assays(lclse)$logCPM[1:3, 1:7]
TCGA.3L.AA1B.01A.11R.A37K.07 TCGA.4N.A93T.01A.11R.A37K.07 TCGA.4T.AA8H.01A.11R.A41B.07
1                     0.250769                      3.51782                   -0.6124583
2                    10.012995                      6.66845                    6.4718439
9                    -6.859130                     -6.85913                   -6.8591298
TCGA.5M.AAT4.01A.11R.A41B.07 TCGA.5M.AAT5.01A.21R.A41B.07 TCGA.5M.AAT6.01A.11R.A41B.07
1                    -1.603433                   -0.1085861                     1.564405
2                     7.316449                    7.2480187                     8.840315
9                    -6.859130                   -6.8591298                    -6.859130
TCGA.5M.AATA.01A.31R.A41B.07
1                    0.3822143
2                    9.1577534
9                   -6.8591298

# PLOTTING
#par(mfrow = c(1, 2), mar = c(4, 5, 1, 1))
#multidensity(as.list(as.data.frame(CPM)), xlab = "CPM", legend = NULL, main = "",
            # cex.axis = 1.2, cex.lab = 1.5, las = 1)
#multidensity(as.list(as.data.frame(assays(lclse)$logCPM)), xlab = "log2 CPM", legend = NULL,
            #main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)

#par(mfrow = c(1, 2), mar = c(4, 5, 1, 1)) 
#multidensity(as.list(as.data.frame(assays(lclse)$logCPM)), xlab = "log2 CPM", legend = NULL,
#                                                       main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)
#boxplot(assays(lclse)$logCPM, col = "gray", xlab = "Samples", ylab = expression(log[2] *
#                                                                                  "CPM"), cex.axis = 1.2, cex.lab = 1.5, las = 1)

# Filter out genes based on expression values

# Observe the distribution of expression values
avgexp <- rowMeans(assays(lclse)$logCPM)
hist(avgexp, xlab = expression(log[2] * "CPM"), main = "", las = 1, col = "gray")
abline(v = log(8.3), col = "red", lwd = 2)
cpmcutoff <- round(10/min(dge$sample$lib.size/1e+06), digits = 1)

# Eliminating from cutoff
nsamplescutoff <- min(table(lclse$gender))
mask <- rowSums(cpm(dge) > cpmcutoff) >= nsamplescutoff
lclse.filt <- lclse[mask, ]
dge.filt <- dge[mask, ]
dim(lclse.filt)
#--- DGE FILT NOW CONTAINS LESS GENES. 

# plotting
par(mar = c(4, 5, 1, 1))
h <- hist(avgexp, xlab = expression("Expression level (" * log[2] * "CPM)"), main = "",
          las = 1, col = "grey", cex.axis = 1.2, cex.lab = 1.5)
x <- cut(rowMeans(assays(lclse.filt)$logCPM), breaks = h$breaks)
lines(h$mids, table(x), type = "h", lwd = 10, lend = 1, col = "darkred")
legend("topright", c("All genes", "Filtered genes"), fill = c("grey", "darkred"))

# Quality assesment : MA plots
dge$samples$group <- lclse$gender
table(dge$samples$group)
FEMALE   MALE 
224    250 
dge.filt$samples$group <- lclse$gender
# REMOVING NAS
mask <- !is.na(dge$samples$group)
lclse.filt_na <- lclse.filt[, mask]
dge.filt_na <- dge.filt[,mask ]
# Plotting
plotSmear(dge.filt_na,lowess=TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2) 
abline(h = 0, col = "blue", lwd = 2)

# TMM NORMALIZATION 

#First Plot 
par(mfrow = c(1, 2))
mask_remove_na <- !is.na(dge$samples$group)
dge_na <- dge[,mask_remove_na ]
plotSmear(dge_na, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2) 
abline(h = 0, col = "blue", lwd = 2)

#Second plot
dge.filt_na_tmm <- calcNormFactors(dge.filt_na) 
head(dge.filt_na_tmm$samples$norm.factors)
[1] 1.1020263 1.0046975 1.0434021 0.9481915 1.0686638 0.9775405
head(dge.filt_na_tmm$samples$lib.size * dge.filt_na_tmm$samples$norm.factors)
[1] 38250959 34760643 22631088 20019183 33146612 16831736
plotSmear(dge.filt_na_tmm, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2) 
abline(h = 0, col = "blue", lwd = 2)


# MA PLOTS
summary(lclse.filt_na$gender)
FEMALE   MALE 
224    250 
summary(lclse.filt_na$race)
AMERICAN INDIAN OR ALASKA NATIVE                            ASIAN        BLACK OR AFRICAN AMERICAN 
1                               11                               51 
[Not Available]                            WHITE                             NA 
                             172                              201                               38 
                             
summary(lclse.filt_na$ethnicity)
HISPANIC OR LATINO        [Not Available]        [Not Evaluated] NOT HISPANIC OR LATINO 
1                    180                      3                    251 
[Unknown]                   NA
                     1                     38 
                     
# MA plots - only explorative
par(mfrow=c(1, 2), mar=c(4, 5, 3, 1)) 
for(i in 1:2){
  A <- rowMeans(assays(lclse.filt)$logCPM) ;
  M <- assays(lclse.filt)$logCPM[, i] - A 
  smoothScatter(A, M, main=colnames(lclse.filt)[i], las=1, cex.axis=1.2, cex.lab=1.5, cex.main=2) 
  abline(h=0, col="blue", lwd=2) ;
  lo <- lowess(M ~ A) ; 
  lines(lo$x, lo$y, col="red", lwd=2)
}

# Normalization

# Order values by samples keeping track of the original order
m <- assays(lclse.filt_na)$logCPM
originalOrder <- apply(m, 2, order)
m <- apply(m, 2, sort)
nsamples<- ncol(m)
ngenes<-nrow(m)

# calculate aritmetic mean and assign it.
m <- matrix(rep(floor(rowMeans(m)), nsamples), nrow = ngenes)

# re assign the original order
m <- sapply(1:ncol(m),
            function(x) m[order(originalOrder[, x]), x])
apply(m, 2, sd)


# visualization 
plotMDS(dge.filt_na_tmm, col = c("red", "blue")[as.integer(dge.filt$samples$group)], cex = 0.7) 
legend("topleft", c("female", "male"), fill = c("red", "blue"), inset = 0.05, cex = 0.7)

## In this second plot is needed to incorporate NAs in the legend (black ones in the plot)
plotMDS(dge.filt_na_tmm, col = c("red", "orange", "blue", "green", "cyan")[as.integer(lclse.filt_na$race)], cex = 0.7)
legend("topleft", na. ,levels(lclse.filt_na$race ), fill = c("red", "orange", "blue", "green", "cyan", "black"), inset = 0.05, cex = 0.7)

# eliminate outliers
maskbad <- colnames(lclse.filt_na) %in% c("AB45", "TCGA.F4.6703.01A.11R.1839.07")
table(maskbad)
dim(lclse.filt_na)
[1] 9615  474
dim(dge.filt_na_tmm)
[1] 9615  474
lclse.filt_bad <- lclse.filt_na[, !maskbad]
dge.filt_na_tmm_bad <- dge.filt_na_tmm[, !maskbad]
dim(lclse.filt_bad)
[1] 9615  472
dim(dge.filt_na_tmm_bad)
[1] 9615  472

# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------


# ------------------------------ BATCH EFFECT

head(colData(lclse))
lclse$bcr_patient_barcode

# all samples were sequenced at the same center
center <- substr(colnames(dge.filt_na_tmm), 27, 28)
table(center)

#Tss 
tss <- substr(colnames(dge.filt_na_tmm), 6, 7)
table(tss)
3L  4N  4T  A6  AA  AD  AM  AU  AY  AZ  CA  CK  CM  D5  DM  F4  G4  NH  QG  QL  RU  WS 
1   1   1  49 192  13   2   2  10  27  10  14  37  31  25  17  26   8   5   1   1   1 

patient <- substr(colnames(dge.filt_na_tmm), 9, 12)
table(patient)
2671 2675 2676 2678 2679 2680 2681 2682 2683 2685 2686 3488 3489 3492 3494 3495 3496 3502 3506 3509 
2    2    1    2    2    2    1    2    2    2    2    1    2    1    1    1    2    1    1    1 
3510 3511 3514 3516 3517 3518 3519 3520 3521 3522 3524 3525 3526 3527 3529 3530 3531 3532 3534 3538 
1    2    2    2    2    2    1    2    1    2    1    2    1    2    1    1    2    1    2    1 
3542 3543 3544 3548 3549 3552 3553 3554 3555 3556 3560 3561 3562 3655 3660 3662 3663 3664 3666 3667 
1    1    1    1    1    1    1    1    1    1    1    1    1    2    2    2    2    1    1    1 
3672 3673 3675 3678 3679 3680 3681 3684 3685 3688 3692 3693 3695 3696 3697 3710 3712 3713 3715 3779 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    2    1    2    2    1    1 
3807 3808 3811 3812 3814 3815 3818 3819 3821 3831 3833 3837 3841 3842 3844 3845 3846 3848 3850 3851 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
3852 3854 3855 3856 3858 3860 3861 3862 3864 3866 3867 3869 3870 3872 3875 3877 3930 3939 3941 3947 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
3949 3950 3952 3955 3956 3966 3968 3970 3971 3972 3973 3975 3976 3977 3979 3980 3982 3984 3986 3989 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
3994 4070 4071 4105 4107 4308 4313 4315 4323 4614 4615 4616 4681 4682 4684 4743 4744 4746 4747 4748 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
4750 4751 4752 4947 4948 4950 4951 4952 5254 5255 5256 5341 5344 5348 5349 5403 5407 5537 5538 5539 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
5540 5541 5543 5657 5660 5662 5664 5666 5667 5796 5797 5820 5821 5860 5861 5862 5863 5864 5868 5900 
1    1    1    1    1    2    1    1    2    1    1    1    1    1    1    1    1    1    1    1 
5912 5913 5914 5915 5916 6004 6137 6138 6140 6141 6142 6161 6162 6163 6164 6165 6166 6167 6168 6169 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
6170 6171 6172 6196 6197 6293 6294 6295 6297 6298 6299 6302 6303 6304 6306 6307 6309 6310 6311 6314 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
6315 6320 6321 6322 6323 6386 6459 6460 6461 6463 6529 6530 6531 6532 6533 6534 6535 6536 6537 6538 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
6539 6540 6541 6548 6569 6570 6586 6588 6598 6599 6600 6601 6603 6605 6606 6607 6608 6625 6626 6627 
1    1    1    1    1    1    1    1    2    2    2    2    2    2    1    1    1    1    1    1 
6628 6648 6649 6651 6652 6653 6654 6674 6675 6676 6677 6678 6679 6680 6703 6704 6715 6716 6717 6718 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    2    1    1    1    1 
6719 6746 6747 6748 6751 6782 6805 6806 6807 6808 6809 6854 6855 6856 6888 6889 6890 6895 6898 6899 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
6901 6920 6922 6923 6924 6926 6927 6928 6929 6930 6931 6932 6963 6964 6965 7000 A004 A00A A00D A00E 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
A00F A00J A00K A00L A00N A00O A00Q A00R A00U A00W A00Z A010 A017 A01C A01D A01F A01G A01I A01K A01P 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
A01Q A01R A01S A01T A01V A01X A01Z A022 A024 A029 A02E A02F A02H A02J A02K A02O A02R A02W A02Y A03F 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
A03J A0X9 A0XD A0XF A1D0 A1D4 A1D6 A1D7 A1D8 A1D9 A1DA A1DB A1HA A1HB A280 A282 A285 A288 A28A A28C 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
A28E A28F A28G A28H A28K A28M A50T A50U A50V A54L A565 A566 A567 A56B A5EJ A5EK A5IV A5YV A5YW A5YX 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 
A5Z1 A5Z2 A5ZU A69D A6GA A6GB A6GC A71X A8F8 A8FL A8YK A93T A97D AA1B AA8H AB45 
1    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1 

# Portion analyte
portionanalyte <- substr(colnames(dge.filt_na_tmm), 18, 20)
table(portionanalyte)
01R 02R 03R 11R 12R 21R 23R 31R 33R 42R 43R 72R 
214  37   1 170   7  36   1   4   1   1   1   1 

#Plate
plate <- substr(colnames(dge.filt_na_tmm), 22, 25)
table(plate)
0821 0826 0905 1022 1113 1410 1653 1672 1723 1774 1839 1873 1928 A002 A00A A083 A089 A155 A16W A180 
23    6   32   46    4   30   39    1   56   33   50    3   26   15   12    6    5   12   12    1 
A28H A32Y A32Z A37K A41B 
15    7   25    9    6 

#Sample vial
samplevial <- substr(colnames(dge.filt_na_tmm), 14, 16)
table(samplevial)
01A 01B 11A 
433   3  38 

# Testing 
table(data.frame(TYPE=lclse.filt_na$type, TSS=tss))
table(data.frame(TYPE=lclse.filt_na$type, plate=plate))
table(data.frame(TYPE=lclse.filt_na$type, vial = samplevial))
table(data.frame(TYPE=lclse.filt_na$type, analyte = portionanalyte))
# after having an overview of these, we decide to use 


logCPM <- cpm(dge.filt, log=TRUE, prior.count=3)
d <- as.dist(1-cor(logCPM, method="spearman"))
sampleClustering <- hclust(d)
batch <- as.integer(factor(tss))
sampleDendrogram <- as.dendrogram(sampleClustering, hang=0.1)
names(batch) <- colnames(lclse.filt_na)
outcome <- paste(substr(colnames(lclse.filt_na), 9, 12), as.character(lclse.filt_na$type), sep="-")
names(outcome) <- colnames(lclse.filt_na)
sampleDendrogram <- dendrapply(sampleDendrogram,
                               function(x, batch, labels) {
                                 if (is.leaf(x)) {
                                   attr(x, "nodePar") <- list(lab.col=as.vector(batch[attr(x, "label")]))
                                   attr(x, "label") <- as.vector(labels[attr(x, "label")])
                                 }
                                 x
                               }, batch, outcome)
plot(sampleDendrogram, main="Hierarchical clustering of samples")
legend("topright", paste("Batch", sort(unique(batch)), levels(factor(tss))), fill=sort(unique(batch)))


# --------------------------- Selecting the paired Samples for paired test analysis

# Found the ones that are potentially paired  - same patient 

# Find the paired samples and returns a list with the patient IDs that have a paired sample
find_paired_samples <- function(lclse){
  found <- c()
  paired <- c()
  for(barcode in colnames(lclse)){
    patient  <- substr(barcode, 9,12)
    if(is.element(patient,found)){ paired <- c(paired,patient) }
    else{ found <- c(found, patient) }
  }
  return(paired)
}

# Function that checks, given a vector of samples of the same patient, wether these 
# are both from tumor and non tumor samples. If they are from both, it returns TRUE. 
is_t_and_nt<- function(patient_samples){
  t_found <- FALSE; nt_found<- FALSE
  for (patient_sample in patient_samples){
    t_id <- substr(patient_sample, 14, 15)
    if(startsWith(t_id, "0")){  t_found<-TRUE }
    else if (startsWith(t_id, "1")){ nt_found<-TRUE }
  } 
  if (t_found && nt_found){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

# Filters out from the paired ones the ones that do not have both tumor and non tumor samples available
filter_ones_with_t_and_nt <- function(paired, lclse){
  filtered_paired <- c()
  for (patientcode in paired){
    patient_samples <- c()
    for(barcode in colnames(lclse)){
      patient  <- substr(barcode, 9,12)
      if(patient == patientcode){
        patient_samples <- c(patient_samples,barcode)
      }
    }
    if (is_t_and_nt(patient_samples)){
      filtered_paired <- c(filtered_paired,patient_samples)
    }
  }
  return(filtered_paired)
}



# Get the paired, tumor and non tumor samples
paired <- find_paired_samples(lclse)
filter_ones_with_t_and_nt(paired, lclse)

mask_paired <- is.element(colnames(lclse),filtered_paired)
lclse.paired <- lclse[,mask_paired]
dge.paired <- dge[,mask_paired]
table(lclse.paired$type)

lclse<- lclse.paired
dge<-dge.paired

# Having an overview 
ord <- order(dge$sample$lib.size)
barplot(dge$sample$lib.size[ord]/1e+06, las = 1, ylab = "Millions of reads", xlab = "Samples",
        col = c("cyan", "orange")[lclse$gender[ord]], border=NA )
legend("topleft", c("Female", "Male"), fill = c("cyan", "orange"), inset = 0.01)

#CPM scaling - within sample normalization
CPM <- t(t(dge$counts)/(dge$samples$lib.size/1e+06))
assays(lclse)$logCPM <- cpm(dge, log = TRUE, prior.count = 0.25) 
assays(lclse)$logCPM[1:3, 1:7]


