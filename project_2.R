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
lclse<-data
dge <- DGEList(counts = assays(lclse)$counts, genes = as.data.frame(rowData(lclse)))

# ------------------------------------------------------------------------------------------------------------
# Create a subset - paired 
# Only select the paired samples
df <- data.frame(rbind(table(substr(colnames(lclse),9,12), substr(colnames(lclse),14,15))))
df_filt <- df[df$X01 + df$X02 + df$X06> 0,]
df_paired <- df_filt[df_filt$X11>0,]
paired_mask <- substr(colnames(lclse),9,12) %in% rownames(df_paired)

lclse.paired<- lclse[,paired_mask]
dge.paired <- dge[,paired_mask]
table(lclse.paired$type)


# ------------------------------------------------------------------------------------------------------------
# Analyze the sequencing depth
sampledepth <- round(dge.paired$sample$lib.size / 1e6, digits=1)
names(sampledepth) <- substr(colnames(lclse.paired), 6, 12)
sort(sampledepth)
ord <- order(dge.paired$sample$lib.size)
barplot(dge.paired$sample$lib.size[ord]/1e+06, las = 1, ylab = "Millions of reads", xlab = "Samples",
        col = c("cyan", "orange")[lclse.paired$type] )
legend("topleft", c("Normal", "Tumor"), fill = c("cyan", "orange"), inset = 0.01)

# TODO ! Possibly remove too poor sequencing depth 

# ------------------------------------------------------------------------------------------------------------
# Overview Gender vs Coverage
ord <- order(dge$sample$lib.size)
barplot(dge$sample$lib.size[ord]/1e+06, las = 1, ylab = "Millions of reads", xlab = "Samples",
        col = c("cyan", "orange")[lclse$gender[ord]], border=NA )
legend("topleft", c("Female", "Male"), fill = c("cyan", "orange"), inset = 0.01)

# ------------------------------------------------------------------------------------------------------------
#CPM scaling - within sample normalization
CPM <- t(t(dge.paired$counts)/(dge.paired$samples$lib.size/1e+06))
assays(lclse.paired)$logCPM <- cpm(dge.paired, log = TRUE, prior.count = 0.25) 
#assays(lclse.paired)$logCPM[1:3, 1:7]


# PLOTTING
tumor_mask <- lclse.paired$type == "tumor"
non_tumor_mask <- lclse.paired$type == "normal"
lclse_tumor <- lclse.paired[,tumor_mask]
lclse_control <- lclse.paired[,non_tumor_mask]
# Plot alone
par(mfrow = c(1, 1), mar = c(4, 5, 1, 1))
multidensity(as.list(as.data.frame(assays(lclse.paired)$logCPM)), xlab = "log2 CPM", legend = NULL,
             main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)

# plot divided by tumor and control for a better overview
par(mfrow = c(1, 2), mar = c(4, 5, 1, 1))
multidensity(as.list(as.data.frame(assays(lclse_tumor)$logCPM)), xlab = "log2 CPM", legend = NULL,
           main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)
multidensity(as.list(as.data.frame(assays(lclse_control)$logCPM)), xlab = "log2 CPM", legend = NULL,
             main = "", cex.axis = 1.2, cex.lab = 1.5, las = 1)




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
head(dge.filt_na_tmm$samples$lib.size * dge.filt_na_tmm$samples$norm.factors)
plotSmear(dge.filt_na_tmm, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2) 
abline(h = 0, col = "blue", lwd = 2)


# MA PLOTS
summary(lclse.filt_na$gender)
summary(lclse.filt_na$race)
summary(lclse.filt_na$ethnicity)

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

plotMDS(dge.filt_na_tmm, col = c("red", "orange", "blue")[as.integer(lclse$concentration)], cex = 0.7)
legend("topleft", levels(lclse$concentration), fill = c("red", "orange", "blue"), inset = 0.05, cex = 0.7)
# missing: eliminate outliers

# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------


# ------------------------------ BATCH EFFECT

head(colData(lclse))
lclse$bcr_patient_barcode

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

resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}

