library(readxl)
ds <- read_xlsx("TIO2+PTYR-human-MSS+MSIvsPD.XLSX")

# BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)

str(ds)
assay <- data.frame(ds[,5:16]) # només conté les mesures
row_data <- ds[,1] # conté la informació de les files que no son mesures del grups MSS i PD
col_data <- data.frame( # Conté la informació de les columnes de mesures
  sample_id = colnames(assay),
  group = c(rep('MSS',6), rep('PD',6)),
  replica = c(rep(c(1,2),6))
)

se <- SummarizedExperiment(
  assays = list(counts = assay),
  rowData = row_data,
  colData = col_data
)

# Es guarda l'objecte SE en binari
save(se, file = "dades_tumors.Rda")

# Es guarden les dades en format text
write.table(assay(se), file = "dades_assay.txt", sep = "\t", row.names = TRUE)

# Primer es comprova les dimensions de les dades obtingudes
dim(se)

# Es mira el tipus de cada columna i la seva distribució
str(assay(se))
summary(assay(se))

# Diferencies entre la primera i la segona repeticio (comparació de grups)
# S'extreuen les dades del SummarizedExperiment pel seu tractament, aquests objectes ja els
#tenim definits de la creació de l'objecte se, però a mode ilustratiu, s'extreuen les dades del se

counts <- assay(se)
groups <- colData(se)$group
repeticions <- colData(se)$replica

# Boxplots
library(ggplot2)
par(las = 2)
boxplot(counts, main = 'Boxplot de les mostres')

# Es passen els valors a escala logarítmica en base 10
par(las = 2)
boxplot(log10(counts), main = 'Boxplot de les mostres en log10')

# Es volen comparar mostres independents: MSS1 vs MSS2 
c_mss_1 <- counts[groups == "MSS" & repeticions == 1]
c_mss_2 <- counts[groups == "MSS" & repeticions == 2]
c_pd_1 <- counts[groups == "PD" & repeticions == 1]
c_pd_2 <- counts[groups == "PD" & repeticions == 2]

p_mss <- c()
for (i in 1:dim(c_mss_1)[2]){
  p_mss <- c(p_mss,t.test(c_mss_1[,i], c_mss_2[,i])$p.value)
}

p_pd <- c()
for (i in 1:dim(c_pd_1)[2]){
  p_pd <- c(p_pd,t.test(c_pd_1[,i], c_pd_2[,i])$p.value)
}

# Diferencies entre els dos grups MSS i PD per repetició si són diferents, general si son iguals
c_mss <- counts[groups == "MSS"]
c_pd <- counts[groups == "PD"]

p <- c()
for (i in 1:dim(c_mss)[1]){
  p <- c(p,t.test(c_mss[i,], c_pd[i,])$p.value)
}

# Es comproba quants registres rebutjen la H_0
length(which(p < 0.05))

# S'aplica una PCA
pcs <- prcomp(counts)

# Plot de la PCA
plot(pcs$rotation[,1], pcs$rotation[,2],
     main="Representació dels dos PCA",
     xlab = paste("PCA1", round(c(pcs$sdev^2 / sum(pcs$sdev^2))[1]*100,2), "%"),
     ylab = paste("PCA2", round(c(pcs$sdev^2 / sum(pcs$sdev^2))[2]*100,2), "%"),
     xlim = c(0, 0.6), ylim = c(-0.5, 0.8),
     col = as.factor(groups), pch = 19)
text(pcs$rotation[,1], pcs$rotation[,2], colnames(counts), cex=0.5, pos=3)

