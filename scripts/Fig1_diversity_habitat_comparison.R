library('ape')
library('picante')

source('./scripts/spat_functions.R')

# first compute PD
tree = read.tree('./data/S1_FSP_97_DNA_RDPalign_tree.tree')
root_tip = 'FIB10_D58U3'
tree = root(tree, outgroup=root_tip, resolve.root=TRUE) ## root tree with ist sequence bc tree was not rooted
tree

abu = read.csv('./data/S1_FSP_97_DNA_name_abu_file.csv')
abu = as.data.frame(t(abu))                     ##transpose pgylogenetic data
names(abu) = as.matrix(abu[1,  ])               ## change the data as matrix
rownames_abu = as.matrix(rownames(abu))         ##change the row names as matrix
abu = abu[-1, ]                                 ## ----except first row
abu = apply(abu, 2, as.numeric)                 ## set all columns to numeric class
rownames(abu) = rownames_abu[-1,]               ## Adjust the row names

dim(abu)
sum(colnames(abu) %in% tree$tip.label)

# phylogenetic diversity 
PD_root = pd(abu, tree, include.root=TRUE)

# now read in habitat information
coords = read.csv('./data/Spatial_coordinates/treatment_forest_pasture_secForest_dd_coords.csv')

# make graphics
cor(PD_root$PD, PD_root$SR)

PD_root$habitat = coords$systems[match(rownames(PD_root), coords$sample)]

boxplot(PD_root$PD ~ PD_root$habitat)
boxplot(PD_root$SR ~ PD_root$habitat)

PDmod = lm(PD ~ habitat, data=PD_root)
SRmod = lm(SR ~ habitat, data=PD_root)

PD_avg = tapply(PD_root$PD, PD_root$habitat, mean)
SR_avg = tapply(PD_root$SR, PD_root$habitat, mean)

PD_SE = tapply(PD_root$PD, PD_root$habitat, function(x) sd(x) / sqrt(length(x)))
SR_SE = tapply(PD_root$SR, PD_root$habitat, function(x) sd(x) / sqrt(length(x)))

PD_CI = 1.96 * PD_SE
SR_CI = 1.96 * SR_SE

png('./figures/Fig1_div_hab_compare.png')
  par(mfrow=c(1,1))
  par(mar = c(4.1, 4.1, 4.1, 5.1))
  xlim = c(.75, 3.4)
  ylim = range(SR_avg + SR_CI, SR_avg - SR_CI)
  index = c(1,3,2)
  xlab = 'Habitat Type'
  ylab = 'Taxonomic Richness'
  SR_pch = c(2, 0, 1)
  PD_pch = c(17, 15, 19)
  col = c(1, 'grey50', 1)
  plot(1:3, SR_avg[index], xlim=xlim, ylim=ylim, xlab='', ylab='',
       type='n', axes=F, cex=1.25)
  addAxis(2)
  addAxis(1, lab=c('Forest', 'Secondary\nForest', 'Pasture'), at=1:3, cex.axis=1.5)
  mtext(side=2, ylab, padj=-2.5, cex=1.75)
  arrows(1:3, (SR_avg - SR_CI)[index], 1:3, (SR_avg + SR_CI)[index], 
         code=3, angle=90, length=.1, lwd=4)
  points(1:3, SR_avg[index], pch=SR_pch, cex=3, col=col, lwd=3)
  par(new=TRUE)
  x_incr = .2
  ylim = range(PD_avg + PD_CI, PD_avg - PD_CI)
  ylab = 'Phylogenetic Diversity'
  plot(1:3, PD_avg[index], xlim=xlim, ylim=ylim, xlab='', ylab='',
       type='n', axes=F, cex=1.25)
  arrows(1:3 + x_incr, (PD_avg - PD_CI)[index], 1:3 + x_incr, (PD_avg + PD_CI)[index], 
         code=3, angle=90, length=.1, lwd=4)
  points(1:3 + x_incr, PD_avg[index], pch=PD_pch, cex=3, col=col)
  addAxis(4)
  mtext(side=4, ylab, padj=3.25, cex=1.75)
dev.off()

## emphasize habitat comparison
png('./figures/Fig1_div_hab_compare_emphasis.png')
  par(mfrow=c(1,1))
  ## adjust margin widths to accomodate larger right margin
  par(mar = c(4.1, 4.1, 4.1, 5.1))
  xlim = c(.75, 3.4)
  ylim = range(SR_avg + SR_CI, SR_avg - SR_CI)
  index = c(1,3,2)
  xlab = 'Habitat Type'
  ylab = 'Taxonomic Richness'
  ## set point symbols
  SR_pch = c(2, 0, 1)
  PD_pch = c(17, 15, 19)
  col = 1
  plot(1:3, SR_avg[index], xlim=xlim, ylim=ylim, xlab='', ylab='',
       type='n', axes=F, cex=1.25)
  addAxis(2)
  addAxis(1, lab=c('Forest', 'Secondary\nForest', 'Pasture'), at=1:3, cex.axis=1.5)
  mtext(side=2, ylab, padj=-2.5, cex=1.75)
  arrows(1:3, (SR_avg - SR_CI)[index], 1:3, (SR_avg + SR_CI)[index], 
         code=3, angle=90, length=.1, lwd=4, col=col)
  points(1:3, SR_avg[index], pch=SR_pch, cex=3, col=col, lwd=3)
  ## overlay additional plot on existing plot
  par(new=TRUE)
  x_incr = .2
  col='grey50'
  ylim = range(PD_avg + PD_CI, PD_avg - PD_CI)
  ylab = 'Phylogenetic Diversity'
  plot(1:3, PD_avg[index], xlim=xlim, ylim=ylim, xlab='', ylab='',
       type='n', axes=F, cex=1.25)
  arrows(1:3 + x_incr, (PD_avg - PD_CI)[index], 1:3 + x_incr, (PD_avg + PD_CI)[index], 
         code=3, angle=90, length=.1, lwd=4, col=col)
  points(1:3 + x_incr, PD_avg[index], pch=PD_pch, cex=3, col=col)
  addAxis(4)
  mtext(side=4, ylab, padj=3.25, cex=1.75)
dev.off()

