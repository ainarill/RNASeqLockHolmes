# ---------------- PACKAGES 
suppressPackageStartupMessages({
    library(edgeR)
    library(SummarizedExperiment)
    library(geneplotter)
    library(BiocStyle)
  }
)
# ------------------------------------------------------------------------------------------------------------

# GET THE DATA 
coadse <- readRDS("seCOAD.rds")
dge <- DGEList(counts = assays(coadse)$counts, genes = as.data.frame(rowData(coadse)))

# ------------------------------------------------------------------------------------------------------------
# Create a subset - paired 
# Only select the paired samples
df <- data.frame(rbind(table(substr(colnames(coadse),9,12), substr(colnames(coadse),14,15))))
df_filt <- df[df$X01 + df$X02 + df$X06> 0,]
df_paired <- df_filt[df_filt$X11>0,]
paired_mask <- substr(colnames(coadse),9,12) %in% rownames(df_paired)

coadse.paired<- coadse[,paired_mask]
dge.paired <- dge[,paired_mask]
table(coadse.paired$type)
normal  tumor 
41     46 
A6.2679
A6.2682

# ------------------------------------------------------------------------------------------------------------
# Analyze the sequencing depth
sampledepth <- round(dge.paired$sample$lib.size / 1e6, digits=1)
names(sampledepth) <- substr(colnames(coadse.paired), 6, 12)
sort(sampledepth)
A6.2679 A6.2682 AA.3531 AA.3517 AA.3534 A6.2683 A6.2684 AA.3518 A6.5665 AA.3520 AA.3522 AA.3516 
5.9     9.0    13.6    14.7    15.4    15.8    15.9    16.5    16.6    16.8    17.0    17.2 
AA.3525 A6.2684 AA.3514 A6.2680 AA.3527 A6.2685 A6.5659 A6.5659 A6.2671 A6.2678 AA.3496 AA.3697 
17.5    17.7    17.7    18.1    19.0    19.2    19.2    20.3    21.1    21.1    25.3    25.6 
AA.3660 AZ.6598 AA.3663 A6.2675 AA.3663 AZ.6600 AA.3662 AA.3662 A6.5665 AA.3713 A6.2684 A6.5662 
28.4    30.2    30.4    31.0    31.1    31.2    31.5    32.5    32.5    32.7    32.8    33.0 
A6.2675 AA.3655 AA.3660 AZ.6599 AA.3712 AZ.6605 AA.3489 AZ.6601 AZ.6603 A6.5665 AA.3525 AA.3655 
33.4    33.5    33.7    34.2    34.5    34.7    34.7    35.1    35.3    35.5    36.8    37.2 
A6.5667 AA.3534 A6.2686 AA.3511 A6.5667 A6.5659 AZ.6601 A6.2679 AZ.6598 A6.5662 AA.3496 A6.2684 
37.3    37.7    37.8    37.8    38.2    38.4    38.6    39.1    39.2    39.5    39.8    40.2 
AA.3511 AA.3712 AZ.6605 AA.3518 AZ.6599 AZ.6600 AA.3713 AA.3697 AA.3514 AZ.6603 AA.3489 A6.2682 
40.3    40.4    40.7    40.8    40.9    41.9    42.0    44.4    44.9    45.0    45.2    45.8 
A6.2678 A6.2671 F4.6704 AA.3516 A6.2685 A6.2683 A6.2686 A6.5659 A6.2680 AA.3531 AA.3522 AA.3517 
46.5    47.3    47.3    47.7    48.6    48.8    49.9    52.6    53.3    53.8    55.6    56.0 
AA.3527 AA.3520 F4.6704 
57.1    57.8    67.5 


paired_mask <- substr(colnames(coadse),9,12) %in% rownames(df_paired)
maskbad <- substr(colnames(coadse.paired$sample), 9, 12) %in% c("2679", "2682")


par(mfrow = c(1, 1), mar = c(4, 5, 1, 1))
ord <- order(dge.paired$sample$lib.size)
barplot(dge.paired$sample$lib.size[ord]/1e+06, las = 1, ylab = "Millions of reads", xlab = "Samples",
        col = c("cyan", "orange")[coadse.paired$type] )
legend("topleft", c("Normal", "Tumor"), fill = c("cyan", "orange"), inset = 0.01)

# TODO ! Possibly remove too poor sequencing depth 

# ------------------------------------------------------------------------------------------------------------
# Overview Gender vs Coverage
ord <- order(dge.paired$sample$lib.size)
barplot(dge.paired$sample$lib.size[ord]/1e+06, las = 1, ylab = "Millions of reads", xlab = "Samples",
        col = c("cyan", "orange")[coadse.paired$gender[ord]] )
legend("topleft", c("Female", "Male"), fill = c("cyan", "orange"), inset = 0.01)

# ------------------------------------------------------------------------------------------------------------
#CPM scaling - within sample normalization
CPM <- t(t(dge.paired$counts)/(dge.paired$samples$lib.size/1e+06))
assays(coadse.paired)$logCPM <- cpm(dge.paired, log = TRUE, prior.count = 0.25) 
#assays(coadse.paired)$logCPM[1:3, 1:7]
# visualizing
tumor_mask <- coadse.paired$type == "tumor"
non_tumor_mask <- coadse.paired$type == "normal"
coadse_tumor <- coadse.paired[,tumor_mask]
coadse_control <- coadse.paired[,non_tumor_mask]
# Plot alone
par(mfrow = c(1, 1), mar = c(4, 5, 1, 1))
multidensity(as.list(as.data.frame(assays(coadse.paired)$logCPM)), xlab = "log2 CPM", legend = NULL,
             main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)

# plot divided by tumor and control for a better overview
par(mfrow = c(1, 2), mar = c(4, 5, 1, 1))
multidensity(as.list(as.data.frame(assays(coadse_tumor)$logCPM)), xlab = "log2 CPM", legend = NULL,
           main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)
multidensity(as.list(as.data.frame(assays(coadse_control)$logCPM)), xlab = "log2 CPM", legend = NULL,
             main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)


# -------------------------------------------------------
# Distribution of expression among genes
avgexp <- rowMeans(assays(coadse.paired)$logCPM)
par(mfrow = c(1, 1), mar = c(4, 5, 1, 1))
hist(avgexp, xlab = expression(log[2] * "CPM"), main = "", las = 1, col = "gray")
abline(v = log(8.3), col = "red", lwd = 2)
cpmcutoff <- round(10/min(dge$sample$lib.size/1e+06), digits = 1)
# Filter out genes based on expression values



# Eliminating from cutoff
nsamplescutoff <- min(table(coadse$gender))
mask <- rowSums(cpm(dge) > cpmcutoff) >= nsamplescutoff
coadse.filt <- coadse[mask, ]
dge.filt <- dge[mask, ]
dim(coadse.filt)
#--- DGE FILT NOW CONTAINS LESS GENES. 

# plotting
par(mar = c(4, 5, 1, 1))
h <- hist(avgexp, xlab = expression("Expression level (" * log[2] * "CPM)"), main = "",
          las = 1, col = "grey", cex.axis = 1.2, cex.lab = 1.5)
x <- cut(rowMeans(assays(coadse.filt)$logCPM), breaks = h$breaks)
lines(h$mids, table(x), type = "h", lwd = 10, lend = 1, col = "darkred")
legend("topright", c("All genes", "Filtered genes"), fill = c("grey", "darkred"))

# Quality assesment : MA plots
dge$samples$group <- coadse$gender
table(dge$samples$group)
dge.filt$samples$group <- coadse$gender
# REMOVING NAS
mask <- !is.na(dge$samples$group)
coadse.filt_na <- coadse.filt[, mask]
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
head(dge.filt_na_tmm$samples$lib.size * dge.filt_na_tmm$samples$norm.factors)
plotSmear(dge.filt_na_tmm, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2) 
abline(h = 0, col = "blue", lwd = 2)


# MA PLOTS
summary(coadse.filt_na$gender)
summary(coadse.filt_na$race)
summary(coadse.filt_na$ethnicity)

# MA plots - only explorative
par(mfrow=c(1, 2), mar=c(4, 5, 3, 1)) 
for(i in 1:2){
  A <- rowMeans(assays(coadse.filt)$logCPM) ;
  M <- assays(coadse.filt)$logCPM[, i] - A 
  smoothScatter(A, M, main=colnames(coadse.filt)[i], las=1, cex.axis=1.2, cex.lab=1.5, cex.main=2) 
  abline(h=0, col="blue", lwd=2) ;
  lo <- lowess(M ~ A) ; 
  lines(lo$x, lo$y, col="red", lwd=2)
}

# Normalization

# Order values by samples keeping track of the original order
m <- assays(coadse.filt_na)$logCPM
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

plotMDS(dge.filt_na_tmm, col = c("red", "orange", "blue")[as.integer(coadse$concentration)], cex = 0.7)
legend("topleft", levels(coadse$concentration), fill = c("red", "orange", "blue"), inset = 0.05, cex = 0.7)
# missing: eliminate outliers

# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------


# ------------------------------ BATCH EFFECT

head(colData(coadse))
coadse$bcr_patient_barcode

# all samples were sequenced at teh same center
center <- substr(colnames(dge.filt_na_tmm), 27, 28)
table(center)

#Tss 
tss <- substr(colnames(dge.filt_na_tmm), 6, 7)
table(tss)

patient <- substr(colnames(dge.filt_na_tmm), 9, 12)
table(patient)

# Portion analyte
portionanalyte <- substr(colnames(dge.filt_na_tmm), 18, 20)
table(portionanalyte)

#Plate
plate <- substr(colnames(dge.filt_na_tmm), 22, 25)
table(plate)

#Sample vial
samplevial <- substr(colnames(dge.filt_na_tmm), 14, 16)
table(samplevial)

# Testing 
table(data.frame(TYPE=coadse.filt_na$type, TSS=tss))
table(data.frame(TYPE=coadse.filt_na$type, plate=plate))
table(data.frame(TYPE=coadse.filt_na$type, vial = samplevial))
table(data.frame(TYPE=coadse.filt_na$type, analyte = portionanalyte))
# after having an overview of these, we decide to use 


logCPM <- cpm(dge.filt, log=TRUE, prior.count=3)
d <- as.dist(1-cor(logCPM, method="spearman"))
sampleClustering <- hclust(d)
batch <- as.integer(factor(tss))
sampleDendrogram <- as.dendrogram(sampleClustering, hang=0.1)
names(batch) <- colnames(coadse.filt_na)
outcome <- paste(substr(colnames(coadse.filt_na), 9, 12), as.character(coadse.filt_na$type), sep="-")
names(outcome) <- colnames(coadse.filt_na)
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

resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}

