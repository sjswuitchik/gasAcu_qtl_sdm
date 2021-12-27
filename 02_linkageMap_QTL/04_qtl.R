#### A5A12 - Klein 1 ####
## May 2021 redo with phased genotypes from LepMap3
library(qtl)
library(viridis)
library(scales)
setwd("~/Desktop/may2021_qtlRedo/klein/matrices/k1")
viridis <- viridis_pal(direction = 1, option = "D")
pal <- viridisLite::viridis(4, option = "D")
pdf('palette.pdf')
show_col(pal)
dev.off()
#### chr 1 ####
# load in data
chr1 <- read.cross(file = "k1_chr1.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr1 <- jittermap(chr1)
# drop null or duplicate markers
chr1 <- drop.nullmarkers(chr1)
chr1 <- drop.dupmarkers(chr1, verbose = T)
# check total markers after pruning
totmar(chr1)
# check for missing data
geno.image(chr1, reorder = T)
# fill in missing data (if any)
chr1 <- fill.geno(chr1, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr1, reorder = T)
# calculate genotype probabilities
chr1 <- calc.genoprob(chr1, step = 2, error.prob = 0.0001)
# simulate genos
chr1 <- sim.geno(chr1, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr1)
map <- pullMap(chr1)
write.csv(map, "chr1_map.csv", row.names = F)
# scan for QTL
chr1s1 <- scanone(chr1, method = "imp", pheno.col = c(2:5))
summary <- summary(chr1s1)
write.csv(summary, "chr1_summary.csv")
# run permutations
perms <- scanone(chr1, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr1_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr1)
pdf('chr1_scan.pdf')
plot(chr1s1, main = "Chr1 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr1s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr1s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr1s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr1s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr1s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr 2 ####
# load in data
chr2 <- read.cross(file = "k1_chr2.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr2 <- jittermap(chr2)
# drop null or duplicate markers
chr2 <- drop.nullmarkers(chr2)
chr2 <- drop.dupmarkers(chr2, verbose = T)
# check total markers after pruning
totmar(chr2)
# check for missing data
geno.image(chr2, reorder = T)
# fill in missing data (if any)
chr2 <- fill.geno(chr2, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr2, reorder = T)
# calculate genotype probabilities
chr2 <- calc.genoprob(chr2, step = 2, error.prob = 0.0001)
# simulate genos
chr2 <- sim.geno(chr2, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr2)
map <- pullMap(chr2)
write.csv(map, "chr2_map.csv", row.names = F)
# scan for QTL
chr2s1 <- scanone(chr2, method = "imp", pheno.col = c(2:5))
summary <- summary(chr2s1)
write.csv(summary, "chr2_summary.csv")
# run permutations
perms <- scanone(chr2, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr2_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr2)
pdf('chr2_scan.pdf')
plot(chr2s1, main = "chr2 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr2s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr2s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr2s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr2s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr2s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr3 - CTmin ####
# load in data
chr3 <- read.cross(file = "k1_chr3.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr3 <- jittermap(chr3)
# drop null or duplicate markers
chr3 <- drop.nullmarkers(chr3)
chr3 <- drop.dupmarkers(chr3, verbose = T)
# check total markers after pruning
totmar(chr3)
# check for missing data
geno.image(chr3, reorder = T)
# fill in missing data (if any)
chr3 <- fill.geno(chr3, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr3, reorder = T)
# calculate genotype probabilities
chr3 <- calc.genoprob(chr3, step = 2, error.prob = 0.0001)
# simulate genos
chr3 <- sim.geno(chr3, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr3)
map <- pullMap(chr3)
write.csv(map, "chr3_map.csv", row.names = F)
# scan for QTL
chr3s1 <- scanone(chr3, method = "imp", pheno.col = c(2:5))
summary <- summary(chr3s1)
write.csv(summary, "chr3_summary.csv")
# run permutations
perms <- scanone(chr3, pheno.col = c(2:5), n.perm = 5000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr3_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr3)
pdf('chr3_scan.pdf')
plot(chr3s1, main = "chr3 - May 2021 5k", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr3s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr3s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr3s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr3s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr3s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

# look at PxG plots for peaks
summary(chr3s1, threshold = 2.554995678, format = c("onepheno"), lodcolumn = 1, perms = perms, pvales = T)
pdf('chr3_ctmin_pxg.pdf')
plotPXG(chr3, marker = "628", pheno.col = 2, jitter = 1)
dev.off()
pdf('chr3_ctmin_eff.pdf')
effectplot(chr3, pheno.col = 2, mname1 = 628)
dev.off()
pdf("chr3_ctmin_effScan.pdf")
effectscan(chr3, pheno.col = 2, chr = 3, get.se = T, draw = T)
dev.off()
lodint(chr3s1, chr = 3, lodcolumn = 1, expandtomarkers = F)
q1 <- makeqtl(chr3, chr = 3, pos = 11.572037, qtl.name = "CTmin", what = "prob")
fitq1 <- fitqtl(chr3, pheno.col = 2, qtl = q1, method = "hk", model = "normal", get.ests = T)
summary(fitq1) # PVE = 28.6%, a = -0.5379 +/- 0.1933, d =  0.6212 +/- 0.2891, p = 0.003852923

#### chr4 ####
# load in data
chr4 <- read.cross(file = "k1_chr4.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr4 <- jittermap(chr4)
# drop null or duplicate markers
chr4 <- drop.nullmarkers(chr4)
chr4 <- drop.dupmarkers(chr4, verbose = T)
# check total markers after pruning
totmar(chr4)
# check for missing data
geno.image(chr4, reorder = T)
# fill in missing data (if any)
chr4 <- fill.geno(chr4, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr4, reorder = T)
# calculate genotype probabilities
chr4 <- calc.genoprob(chr4, step = 2, error.prob = 0.0001)
# simulate genos
chr4 <- sim.geno(chr4, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr4)
map <- pullMap(chr4)
write.csv(map, "chr4_map.csv", row.names = F)
# scan for QTL
chr4s1 <- scanone(chr4, method = "imp", pheno.col = c(2:5))
summary <- summary(chr4s1)
write.csv(summary, "chr4_summary.csv")
# run permutations
perms <- scanone(chr4, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr4_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr4)
pdf('chr4_scan.pdf')
plot(chr4s1, main = "chr4 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr4s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr4s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr4s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr4s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr4s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr5 ####
# load in data
chr5 <- read.cross(file = "k1_chr5.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr5 <- jittermap(chr5)
# drop null or duplicate markers
chr5 <- drop.nullmarkers(chr5)
chr5 <- drop.dupmarkers(chr5, verbose = T)
# check total markers after pruning
totmar(chr5)
# check for missing data
geno.image(chr5, reorder = T)
# fill in missing data (if any)
chr5 <- fill.geno(chr5, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr5, reorder = T)
# calculate genotype probabilities
chr5 <- calc.genoprob(chr5, step = 2, error.prob = 0.0001)
# simulate genos
chr5 <- sim.geno(chr5, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr5)
map <- pullMap(chr5)
write.csv(map, "chr5_map.csv", row.names = F)
# scan for QTL
chr5s1 <- scanone(chr5, method = "imp", pheno.col = c(2:5))
summary <- summary(chr5s1)
write.csv(summary, "chr5_summary.csv")
# run permutations
perms <- scanone(chr5, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr5_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr5)
pdf('chr5_scan.pdf')
plot(chr5s1, main = "chr5 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr5s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr5s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr5s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr5s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr5s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr6 ####
# load in data
chr6 <- read.cross(file = "k1_chr6.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr6 <- jittermap(chr6)
# drop null or duplicate markers
chr6 <- drop.nullmarkers(chr6)
chr6 <- drop.dupmarkers(chr6, verbose = T)
# check total markers after pruning
totmar(chr6)
# check for missing data
geno.image(chr6, reorder = T)
# fill in missing data (if any)
chr6 <- fill.geno(chr6, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr6, reorder = T)
# calculate genotype probabilities
chr6 <- calc.genoprob(chr6, step = 2, error.prob = 0.0001)
# simulate genos
chr6 <- sim.geno(chr6, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr6)
map <- pullMap(chr6)
write.csv(map, "chr6_map.csv", row.names = F)
# scan for QTL
chr6s1 <- scanone(chr6, method = "imp", pheno.col = c(2:5))
summary <- summary(chr6s1)
write.csv(summary, "chr6_summary.csv")
# run permutations
perms <- scanone(chr6, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr6_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr6)
pdf('chr6_scan.pdf')
plot(chr6s1, main = "chr6 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr6s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr6s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr6s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr6s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr6s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr7 ####
# load in data
chr7 <- read.cross(file = "k1_chr7.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr7 <- jittermap(chr7)
# drop null or duplicate markers
chr7 <- drop.nullmarkers(chr7)
chr7 <- drop.dupmarkers(chr7, verbose = T)
# check total markers after pruning
totmar(chr7)
# check for missing data
geno.image(chr7, reorder = T)
# fill in missing data (if any)
chr7 <- fill.geno(chr7, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr7, reorder = T)
# calculate genotype probabilities
chr7 <- calc.genoprob(chr7, step = 2, error.prob = 0.0001)
# simulate genos
chr7 <- sim.geno(chr7, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr7)
map <- pullMap(chr7)
write.csv(map, "chr7_map.csv", row.names = F)
# scan for QTL
chr7s1 <- scanone(chr7, method = "imp", pheno.col = c(2:5))
summary <- summary(chr7s1)
write.csv(summary, "chr7_summary.csv")
# run permutations
perms <- scanone(chr7, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr7_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr7)
pdf('chr7_scan.pdf')
plot(chr7s1, main = "chr7 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr7s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr7s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr7s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr7s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr7s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr8 ####
# load in data
chr8 <- read.cross(file = "k1_chr8.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr8 <- jittermap(chr8)
# drop null or duplicate markers
chr8 <- drop.nullmarkers(chr8)
chr8 <- drop.dupmarkers(chr8, verbose = T)
# check total markers after pruning
totmar(chr8)
# check for missing data
geno.image(chr8, reorder = T)
# fill in missing data (if any)
chr8 <- fill.geno(chr8, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr8, reorder = T)
# calculate genotype probabilities
chr8 <- calc.genoprob(chr8, step = 2, error.prob = 0.0001)
# simulate genos
chr8 <- sim.geno(chr8, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr8)
map <- pullMap(chr8)
write.csv(map, "chr8_map.csv", row.names = F)
# scan for QTL
chr8s1 <- scanone(chr8, method = "imp", pheno.col = c(2:5))
summary <- summary(chr8s1)
write.csv(summary, "chr8_summary.csv")
# run permutations
perms <- scanone(chr8, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr8_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr8)
pdf('chr8_scan.pdf')
plot(chr8s1, main = "chr8 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr8s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr8s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr8s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr8s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr8s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr9 ####
# load in data
chr9 <- read.cross(file = "k1_chr9.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr9 <- jittermap(chr9)
# drop null or duplicate markers
chr9 <- drop.nullmarkers(chr9)
chr9 <- drop.dupmarkers(chr9, verbose = T)
# check total markers after pruning
totmar(chr9)
# check for missing data
geno.image(chr9, reorder = T)
# fill in missing data (if any)
chr9 <- fill.geno(chr9, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr9, reorder = T)
# calculate genotype probabilities
chr9 <- calc.genoprob(chr9, step = 2, error.prob = 0.0001)
# simulate genos
chr9 <- sim.geno(chr9, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr9)
map <- pullMap(chr9)
write.csv(map, "chr9_map.csv", row.names = F)
# scan for QTL
chr9s1 <- scanone(chr9, method = "imp", pheno.col = c(2:5))
summary <- summary(chr9s1)
write.csv(summary, "chr9_summary.csv")
# run permutations
perms <- scanone(chr9, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr9_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr9)
pdf('chr9_scan.pdf')
plot(chr9s1, main = "chr9 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr9s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr9s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr9s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr9s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr9s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr10 ####
# load in data
chr10 <- read.cross(file = "k1_chr10.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr10 <- jittermap(chr10)
# drop null or duplicate markers
chr10 <- drop.nullmarkers(chr10)
chr10 <- drop.dupmarkers(chr10, verbose = T)
# check total markers after pruning
totmar(chr10)
# check for missing data
geno.image(chr10, reorder = T)
# fill in missing data (if any)
chr10 <- fill.geno(chr10, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr10, reorder = T)
# calculate genotype probabilities
chr10 <- calc.genoprob(chr10, step = 2, error.prob = 0.0001)
# simulate genos
chr10 <- sim.geno(chr10, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr10)
map <- pullMap(chr10)
write.csv(map, "chr10_map.csv", row.names = F)
# scan for QTL
chr10s1 <- scanone(chr10, method = "imp", pheno.col = c(2:5))
summary <- summary(chr10s1)
write.csv(summary, "chr10_summary.csv")
# run permutations
perms <- scanone(chr10, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr10_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr10)


pdf('chr10_scan.pdf')
plot(chr10s1, main = "chr10 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr10s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr10s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr10s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr10s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr10s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr11 ####
# load in data
chr11 <- read.cross(file = "k1_chr11.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr11 <- jittermap(chr11)
# drop null or duplicate markers
chr11 <- drop.nullmarkers(chr11)
chr11 <- drop.dupmarkers(chr11, verbose = T)
# check total markers after pruning
totmar(chr11)
# check for missing data
geno.image(chr11, reorder = T)
# fill in missing data (if any)
chr11 <- fill.geno(chr11, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr11, reorder = T)
# calculate genotype probabilities
chr11 <- calc.genoprob(chr11, step = 2, error.prob = 0.0001)
# simulate genos
chr11 <- sim.geno(chr11, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr11)
map <- pullMap(chr11)
write.csv(map, "chr11_map.csv", row.names = F)
# scan for QTL
chr11s1 <- scanone(chr11, method = "imp", pheno.col = c(2:5))
summary <- summary(chr11s1)
write.csv(summary, "chr11_summary.csv")
# run permutations
perms <- scanone(chr11, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr11_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr11)


pdf('chr11_scan.pdf')
plot(chr11s1, main = "chr11 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr11s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr11s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr11s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr11s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr11s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr12 ####
# load in data
chr12 <- read.cross(file = "k1_chr12.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr12 <- jittermap(chr12)
# drop null or duplicate markers
chr12 <- drop.nullmarkers(chr12)
chr12 <- drop.dupmarkers(chr12, verbose = T)
# check total markers after pruning
totmar(chr12)
# check for missing data
geno.image(chr12, reorder = T)
# fill in missing data (if any)
chr12 <- fill.geno(chr12, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr12, reorder = T)
# calculate genotype probabilities
chr12 <- calc.genoprob(chr12, step = 2, error.prob = 0.0001)
# simulate genos
chr12 <- sim.geno(chr12, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr12)
map <- pullMap(chr12)
write.csv(map, "chr12_map.csv", row.names = F)
# scan for QTL
chr12s1 <- scanone(chr12, method = "imp", pheno.col = c(2:5))
summary <- summary(chr12s1)
write.csv(summary, "chr12_summary.csv")
# run permutations
perms <- scanone(chr12, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr12_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr12)


pdf('chr12_scan.pdf')
plot(chr12s1, main = "chr12 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr12s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr12s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr12s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr12s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr12s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr13 ####
# load in data
chr13 <- read.cross(file = "k1_chr13.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr13 <- jittermap(chr13)
# drop null or duplicate markers
chr13 <- drop.nullmarkers(chr13)
chr13 <- drop.dupmarkers(chr13, verbose = T)
# check total markers after pruning
totmar(chr13)
# check for missing data
geno.image(chr13, reorder = T)
# fill in missing data (if any)
chr13 <- fill.geno(chr13, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr13, reorder = T)
# calculate genotype probabilities
chr13 <- calc.genoprob(chr13, step = 2, error.prob = 0.0001)
# simulate genos
chr13 <- sim.geno(chr13, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr13)
map <- pullMap(chr13)
write.csv(map, "chr13_map.csv", row.names = F)
# scan for QTL
chr13s1 <- scanone(chr13, method = "imp", pheno.col = c(2:5))
summary <- summary(chr13s1)
write.csv(summary, "chr13_summary.csv")
# run permutations
perms <- scanone(chr13, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr13_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr13)


pdf('chr13_scan.pdf')
plot(chr13s1, main = "chr13 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr13s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr13s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr13s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr13s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr13s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr14 ####
# load in data
chr14 <- read.cross(file = "k1_chr14.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr14 <- jittermap(chr14)
# drop null or duplicate markers
chr14 <- drop.nullmarkers(chr14)
chr14 <- drop.dupmarkers(chr14, verbose = T)
# check total markers after pruning
totmar(chr14)
# check for missing data
geno.image(chr14, reorder = T)
# fill in missing data (if any)
chr14 <- fill.geno(chr14, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr14, reorder = T)
# calculate genotype probabilities
chr14 <- calc.genoprob(chr14, step = 2, error.prob = 0.0001)
# simulate genos
chr14 <- sim.geno(chr14, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr14)
map <- pullMap(chr14)
write.csv(map, "chr14_map.csv", row.names = F)
# scan for QTL
chr14s1 <- scanone(chr14, method = "imp", pheno.col = c(2:5))
summary <- summary(chr14s1)
write.csv(summary, "chr14_summary.csv")
# run permutations
perms <- scanone(chr14, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr14_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr14)


pdf('chr14_scan.pdf')
plot(chr14s1, main = "chr14 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr14s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr14s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr14s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr14s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr14s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr15 ####
# load in data
chr15 <- read.cross(file = "k1_chr15.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr15 <- jittermap(chr15)
# drop null or duplicate markers
chr15 <- drop.nullmarkers(chr15)
chr15 <- drop.dupmarkers(chr15, verbose = T)
# check total markers after pruning
totmar(chr15)
# check for missing data
geno.image(chr15, reorder = T)
# fill in missing data (if any)
chr15 <- fill.geno(chr15, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr15, reorder = T)
# calculate genotype probabilities
chr15 <- calc.genoprob(chr15, step = 2, error.prob = 0.0001)
# simulate genos
chr15 <- sim.geno(chr15, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr15)
map <- pullMap(chr15)
write.csv(map, "chr15_map.csv", row.names = F)
# scan for QTL
chr15s1 <- scanone(chr15, method = "imp", pheno.col = c(2:5))
summary <- summary(chr15s1)
write.csv(summary, "chr15_summary.csv")
# run permutations
perms <- scanone(chr15, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr15_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr15)


pdf('chr15_scan.pdf')
plot(chr15s1, main = "chr15 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr15s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr15s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr15s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr15s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr15s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr16 ####
# load in data
chr16 <- read.cross(file = "k1_chr16.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr16 <- jittermap(chr16)
# drop null or duplicate markers
chr16 <- drop.nullmarkers(chr16)
chr16 <- drop.dupmarkers(chr16, verbose = T)
# check total markers after pruning
totmar(chr16)
# check for missing data
geno.image(chr16, reorder = T)
# fill in missing data (if any)
chr16 <- fill.geno(chr16, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr16, reorder = T)
# calculate genotype probabilities
chr16 <- calc.genoprob(chr16, step = 2, error.prob = 0.0001)
# simulate genos
chr16 <- sim.geno(chr16, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr16)
map <- pullMap(chr16)
write.csv(map, "chr16_map.csv", row.names = F)
# scan for QTL
chr16s1 <- scanone(chr16, method = "imp", pheno.col = c(2:5))
summary <- summary(chr16s1)
write.csv(summary, "chr16_summary.csv")
# run permutations
perms <- scanone(chr16, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr16_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr16)

pdf('chr16_scan.pdf')
plot(chr16s1, main = "chr16 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr16s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr16s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr16s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr16s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr16s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr17 - CBT ####
# load in data
chr17 <- read.cross(file = "k1_chr17.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr17 <- jittermap(chr17)
# drop null or duplicate markers
chr17 <- drop.nullmarkers(chr17)
chr17 <- drop.dupmarkers(chr17, verbose = T)
# check total markers after pruning
totmar(chr17)
# check for missing data
geno.image(chr17, reorder = T)
# fill in missing data (if any)
chr17 <- fill.geno(chr17, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr17, reorder = T)
# calculate genotype probabilities
chr17 <- calc.genoprob(chr17, step = 2, error.prob = 0.0001)
# simulate genos
chr17 <- sim.geno(chr17, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr17)
map <- pullMap(chr17)
write.csv(map, "chr17_map.csv", row.names = F)
# scan for QTL
chr17s1 <- scanone(chr17, method = "imp", pheno.col = c(2:5))
summary <- summary(chr17s1)
write.csv(summary, "chr17_summary.csv")
# run permutations
perms <- scanone(chr17, pheno.col = c(2:5), n.perm = 5000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr17_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr17)
pdf('chr17_scan.pdf')
plot(chr17s1, main = "chr17 - May 2021 5k", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr17s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr17s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr17s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr17s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr17s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

# look at PxG plots for peaks
summary(chr17s1, threshold = 2.7194388, format = c("onepheno"), lodcolumn = 4, perms = perms, pvales = T)
pdf('chr17_core_pxg.pdf')
plotPXG(chr17, marker = "75", pheno.col = 5, jitter = 1)
dev.off()
pdf('chr17_core_eff.pdf')
effectplot(chr17, pheno.col = 5, mname1 = 75)
dev.off()
pdf("chr17_core_effScan.pdf")
effectscan(chr17, pheno.col = 5, chr = 17, get.se = T, draw = T)
dev.off()
lodint(chr17s1, chr = 17, lodcolumn = 4, expandtomarkers = F)
q1 <- makeqtl(chr17, chr = 17, pos = 80.91639, qtl.name = "Core", what = "prob")
fitq1 <- fitqtl(chr17, pheno.col = 5, qtl = q1, method = "hk", model = "normal", get.ests = T)
summary(fitq1) # PVE = 33.4%, a = -1.0475 +/- 0.3158, d =  0.3819 +/- 0.4546, p = 0.001856709

#### chr18 ####
# load in data
chr18 <- read.cross(file = "k1_chr18.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr18 <- jittermap(chr18)
# drop null or duplicate markers
chr18 <- drop.nullmarkers(chr18)
chr18 <- drop.dupmarkers(chr18, verbose = T)
# check total markers after pruning
totmar(chr18)
# check for missing data
geno.image(chr18, reorder = T)
# fill in missing data (if any)
chr18 <- fill.geno(chr18, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr18, reorder = T)
# calculate genotype probabilities
chr18 <- calc.genoprob(chr18, step = 2, error.prob = 0.0001)
# simulate genos
chr18 <- sim.geno(chr18, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr18)
map <- pullMap(chr18)
write.csv(map, "chr18_map.csv", row.names = F)
# scan for QTL
chr18s1 <- scanone(chr18, method = "imp", pheno.col = c(2:5))
summary <- summary(chr18s1)
write.csv(summary, "chr18_summary.csv")
# run permutations
perms <- scanone(chr18, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr18_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr18)
pdf('chr18_scan.pdf')
plot(chr18s1, main = "chr18 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr18s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr18s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr18s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr18s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr18s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr19 ####
# load in data
chr19 <- read.cross(file = "k1_chr19.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr19 <- jittermap(chr19)
# drop null or duplicate markers
chr19 <- drop.nullmarkers(chr19)
chr19 <- drop.dupmarkers(chr19, verbose = T)
# check total markers after pruning
totmar(chr19)
# check for missing data
geno.image(chr19, reorder = T)
# fill in missing data (if any)
chr19 <- fill.geno(chr19, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr19, reorder = T)
# calculate genotype probabilities
chr19 <- calc.genoprob(chr19, step = 2, error.prob = 0.0001)
# simulate genos
chr19 <- sim.geno(chr19, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr19)
map <- pullMap(chr19)
write.csv(map, "chr19_map.csv", row.names = F)
# scan for QTL
chr19s1 <- scanone(chr19, method = "imp", pheno.col = c(2:5))
summary <- summary(chr19s1)
write.csv(summary, "chr19_summary.csv")
# run permutations
perms <- scanone(chr19, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr19_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr19)
pdf('chr19_scan.pdf')
plot(chr19s1, main = "chr19 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr19s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr19s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr19s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr19s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr19s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr20 ####
# load in data
chr20 <- read.cross(file = "k1_chr20.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr20 <- jittermap(chr20)
# drop null or duplicate markers
chr20 <- drop.nullmarkers(chr20)
chr20 <- drop.dupmarkers(chr20, verbose = T)
# check total markers after pruning
totmar(chr20)
# check for missing data
geno.image(chr20, reorder = T)
# fill in missing data (if any)
chr20 <- fill.geno(chr20, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr20, reorder = T)
# calculate genotype probabilities
chr20 <- calc.genoprob(chr20, step = 2, error.prob = 0.0001)
# simulate genos
chr20 <- sim.geno(chr20, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr20)
map <- pullMap(chr20)
write.csv(map, "chr20_map.csv", row.names = F)
# scan for QTL
chr20s1 <- scanone(chr20, method = "imp", pheno.col = c(2:5))
summary <- summary(chr20s1)
write.csv(summary, "chr20_summary.csv")
# run permutations
perms <- scanone(chr20, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr20_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr20)
pdf('chr20_scan.pdf')
plot(chr20s1, main = "chr20 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr20s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr20s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr20s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr20s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr20s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()

#### chr21 ####
# load in data
chr21 <- read.cross(file = "k1_chr21.csv", format = "csv", na.strings = c("-1","NA"), genotypes = c("AA","AB","BA", "BB"), alleles = c("A","B"), estimate.map = F)
# jitter markers at same position
chr21 <- jittermap(chr21)
# drop null or duplicate markers
chr21 <- drop.nullmarkers(chr21)
chr21 <- drop.dupmarkers(chr21, verbose = T)
# check total markers after pruning
totmar(chr21)
# check for missing data
geno.image(chr21, reorder = T)
# fill in missing data (if any)
chr21 <- fill.geno(chr21, method = c("argmax"), error.prob = 0.0001, map.function = "kosambi")
geno.image(chr21, reorder = T)
# calculate genotype probabilities
chr21 <- calc.genoprob(chr21, step = 2, error.prob = 0.0001)
# simulate genos
chr21 <- sim.geno(chr21, n.draws = 25, step = 0, error.prob = 0.0001, map.function = "kosambi", stepwidth = "fixed")
# summarize and write out map for linkage map viz
summaryMap(chr21)
map <- pullMap(chr21)
write.csv(map, "chr21_map.csv", row.names = F)
# scan for QTL
chr21s1 <- scanone(chr21, method = "imp", pheno.col = c(2:5))
summary <- summary(chr21s1)
write.csv(summary, "chr21_summary.csv")
# run permutations
perms <- scanone(chr21, pheno.col = c(2:5), n.perm = 1000, verbose = T, method = "imp")
perms_summary <- summary(perms)
write.csv(perms_summary, "chr21_perms.csv")
# plot chrom scan and write out as PDF
phenos <- phenames(chr21)
pdf('chr21_scan.pdf')
plot(chr21s1, main = "chr21 - May 2021", ylim = c(0,5), ylab = "LOD score")
for (i in 1:length(phenos)) plot(chr21s1, add = T, lodcolumn = i, col = pal[i])
add.threshold(out = chr21s1, perms = perms, alpha = 0.05, col = pal[1], lty = 2, lodcolumn = 1)
add.threshold(out = chr21s1, perms = perms, alpha = 0.05, col = pal[2], lty = 2, lodcolumn = 2)
add.threshold(out = chr21s1, perms = perms, alpha = 0.05, col = pal[3], lty = 2, lodcolumn = 3)
add.threshold(out = chr21s1, perms = perms, alpha = 0.05, col = pal[4], lty = 2, lodcolumn = 4)
dev.off()
