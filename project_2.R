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
data <- readRDS("/Users/luisasantus/Downloads/seCOAD.rds")
coadse<-data
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



# ------------------------------------------------------------------------------------------------------------
# Analyze the sequencing depth
# coverage_all
sampledepth <- round(dge.paired$sample$lib.size / 1e6, digits=1)
names(sampledepth) <- substr(colnames(coadse.paired), 6, 12)
sort(sampledepth)
par(mfrow = c(1, 1), mar = c(4, 5, 1, 1))
ord <- order(dge.paired$sample$lib.size)
barplot(dge.paired$sample$lib.size[ord]/1e+06, las = 1, ylab = "Millions of reads", xlab = "Samples",
        col = c("cyan", "orange")[coadse.paired$type] )
legend("topleft", c("Normal", "Tumor"), fill = c("cyan", "orange"), inset = 0.01)

# Remove the two samples with very low sequencing depth
mask_remove_low_coverage <- substr(colnames(coadse.paired),9,12) %in% c("2679", "2682")
coadse.paired <- coadse.paired[,!mask_remove_low_coverage]
dge.paired <- dge.paired[,!mask_remove_low_coverage]


# coverage_filtered plot
names(sampledepth) <- substr(colnames(coadse.paired), 6, 12)
sort(sampledepth)
par(mfrow = c(1, 1), mar = c(4, 5, 1, 1))
ord <- order(dge.paired$sample$lib.size)
barplot(dge.paired$sample$lib.size[ord]/1e+06, las = 1, ylab = "Millions of reads", xlab = "Samples",
        col = c("cyan", "orange")[coadse.paired$type] )
legend("topleft", c("Normal", "Tumor"), fill = c("cyan", "orange"), inset = 0.01)


# ------------------------------------------------------------------------------------------------------------
# Overview Gender vs Coverage
# The white columns in the plot correspond to the NAs - gender is not available for some samples!
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
             main = "All samples", cex.axis = 1.2, cex.lab = 1.5, las = 1)

# plot divided by tumor and control for a better overview
par(mfrow = c(1, 2), mar = c(4, 5, 1, 1))
multidensity(as.list(as.data.frame(assays(coadse_tumor)$logCPM)), xlab = "log2 CPM", legend = NULL,
           main = "Tumor samples", cex.axis = 1.2, cex.lab = 1.5, las = 1)
multidensity(as.list(as.data.frame(assays(coadse_control)$logCPM)), xlab = "log2 CPM", legend = NULL,
             main = "Control samples", cex.axis = 1.2, cex.lab = 1.5, las = 1)


# ------------------------------------------------------ FILTER OUT LOWLY EXPRESSED GENES-
# Distribution of expression among genes

# Visualize
avgexp <- rowMeans(assays(coadse.paired)$logCPM)
par(mfrow = c(1, 1), mar = c(4, 5, 1, 1))
hist(avgexp, xlab = expression(log[2] * "CPM"), main = "", las = 1, col = "gray")
abline(v = 0, col = "red", lwd = 2)
cpmcutoff <- round(10/min(dge.paired$sample$lib.size/1e+06), digits = 1)

# Eliminate
nsamplescutoff <- min(table(coadse.paired$gender))
mask <- rowSums(cpm(dge.paired) > cpmcutoff) >= nsamplescutoff
coadse.filt <- coadse.paired[mask, ]
dge.filt <- dge.paired[mask, ]
dim(coadse.filt)
# plotting
par(mar = c(4, 5, 1, 1))
h <- hist(avgexp, xlab = expression("Expression level (" * log[2] * "CPM)"), main = "",
          las = 1, col = "grey", cex.axis = 1.2, cex.lab = 1.5)
x <- cut(rowMeans(assays(coadse.filt)$logCPM), breaks = h$breaks)
lines(h$mids, table(x), type = "h", lwd = 10, lend = 1, col = "darkred")
legend("topright", c("All genes", "Filtered genes"), fill = c("grey", "darkred"))

# Save an unnormalized version 
saveRDS(coadse.filt, file.path(".", "coadse.filt.unnorm.rds"))
saveRDS(dge.filt, file.path(".", "dge.filt.unnorm.rds"))

# ------------------------------------------------------ MA PLOTS
# Quality assesment : MA plots

####### Exploration of gender 
dge.filt$samples$group <- coadse.filt$gender
table(dge.filt$samples$group)
dge.filt$samples$group <- coadse.filt$gender
# REMOVING NAS
mask <- !is.na(dge.filt$samples$group)
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

