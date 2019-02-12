##individual nmds of same size
library(vegan)

trt = read.csv('treatment_file3.csv') ## download data file
dat = read.csv('95%_FSP_S1_DNA.csv')
##set as data frame, transpose, data matrix and set as  a numaric varible

comm = as.data.frame(t(dat)) ## Make the test13 as data frame
names(comm) = as.matrix(comm[1,  ]) ## change the data as matrix
rownames_comm = as.matrix(rownames(comm)) ##change the row names as matrix
comm = comm[-1, ] ## ----except first row
comm = apply(comm, 2, as.numeric) ## set all data as numeric starting from row 2 ???
rownames(comm) = rownames_comm[-1,] ## Adjust the row names
#head(comm)
trt = subset(trt, subset = trt$sample %in% rownames(comm))
## check that rows of comm so that it matches coords
all(trt$sample == rownames(comm))
## if not true then select the same rows that match with the treatment files... subset coords to only samples we have community data for
comm = subset(comm, subset=trt$sample %in% rownames(comm))
## set as a factors
## run ANOSIM analysis 
anosim (comm, grouping=trt$treatment, permutations=999, distance="bray") ## did the anosim indirect ordination data (Bray-curtis)

##metanMDS analysis 
comm.nms = metaMDS(comm) ##metaMDS analysis

plot(comm.nms, display='sites', type='n', xlim=c(-2, 2), ylim=c(-2,2), cex=2) ##plot bank box
xlim=c(-2-2)
ylim=c(-2-2)

points= c(rep(17, 31), rep(22, 29), rep(1, 30)) ##repeat the selected points (3 time each)
bgs= c(rep("Black", 31), rep("grey70", 29), rep("white", 30)) 
points(comm.nms, pch=points, bg=bgs, cex=1.25) ##plot data points with specifed background 
##points(comm.nms, pch=points, cex=2) ## enter points
##
legend('topleft', c('Forest', 'Secondary Forest','Pasture'), pch=c(17, 22, 1), cex=.75) ## add ledgends on the top left (though not correct) cex is for the size of ledgend



#adding multiple figures

##individual nmds of same size
library(vegan)

trt = read.csv('treatment_file3.csv') ## download data file
dat = read.csv('97%_FSP_S1_DNA.csv')
##set as data frame, transpose, data matrix and set as  a numaric varible

comm = as.data.frame(t(dat)) ## Make the test13 as data frame
names(comm) = as.matrix(comm[1,  ]) ## change the data as matrix
rownames_comm = as.matrix(rownames(comm)) ##change the row names as matrix
comm = comm[-1, ] ## ----except first row
comm = apply(comm, 2, as.numeric) ## set all data as numeric starting from row 2 ???
rownames(comm) = rownames_comm[-1,] ## Adjust the row names
#head(comm)
trt = subset(trt, subset = trt$sample %in% rownames(comm))
## check that rows of comm so that it matches coords
all(trt$sample == rownames(comm))
## if not true then select the same rows that match with the treatment files... subset coords to only samples we have community data for
comm = subset(comm, subset=trt$sample %in% rownames(comm))
## set as a factors
## run ANOSIM analysis 
anosim (comm, grouping=trt$treatment, permutations=999, distance="bray") ## did the anosim indirect ordination data (Bray-curtis)

##metanMDS analysis 
comm.nms = metaMDS(comm) ##metaMDS analysis

par(mfrow=c(2,3))



plot(comm.nms, display='sites', type='n', xlim=c(-2, 2), ylim=c(-2,2), cex=2) ##plot bank box
points= c(rep(17, 31), rep(22, 29), rep(1, 30)) ##repeat the selected points (3 time each)
bgs= c(rep("Black", 31), rep("grey70", 29), rep("white", 30)) 
points(comm.nms, pch=points, bg=bgs, cex=1.40) ##plot data points with specifed background 
##points(comm.nms, pch=points, cex=2) ## enter points
##
legend('topleft', c('Forest', 'Secondary Forest','Pasture'), pch=c(17, 22, 1), cex=.75) ## add ledgends on the top left (though not correct) cex is for the size of ledgend




par(mfrow=c(2,2))



###
###fig1 97_DNA
#adding multiple figures

##individual nmds of same size
library(vegan)

trt = read.csv('treatment_file3.csv') ## download data file
dat1 = read.csv('97%_FSP_S1_DNA.csv')
##set as data frame, transpose, data matrix and set as  a numaric varible

comm1 = as.data.frame(t(dat1)) ## Make the test13 as data frame
names(comm1) = as.matrix(comm1[1,  ]) ## change the data as matrix
rownames_comm1 = as.matrix(rownames(comm1)) ##change the row names as matrix
comm1 = comm1[-1, ] ## ----except first row
comm1 = apply(comm1, 2, as.numeric) ## set all data as numeric starting from row 2 ???
rownames(comm1) = rownames_comm1[-1,] ## Adjust the row names
#head(comm)
trt = subset(trt, subset = trt$sample %in% rownames(comm1))
## check that rows of comm so that it matches coords
all(trt$sample == rownames(comm1))
## if not true then select the same rows that match with the treatment files... subset coords to only samples we have community data for
comm1 = subset(comm1, subset=trt$sample %in% rownames(comm1))
## set as a factors
## run ANOSIM analysis 
anosim (comm, grouping=trt$treatment, permutations=999, distance="bray") ## did the anosim indirect ordination data (Bray-curtis)

##metanMDS analysis 
comm.nms1 = metaMDS(comm1) ##metaMDS analysis

par(mfrow=c(2,3))



plot(comm.nms1, display='sites', type='n', xlim=c(-2, 2), ylim=c(-2,2),cex=2 ) ##plot bank box
points= c(rep(17, 31), rep(22, 29), rep(1, 30)) ##repeat the selected points (3 time each)
bgs= c(rep("Black", 31), rep("grey70", 29), rep("white", 30)) 
points(comm.nms1, pch=points, bg=bgs, cex=1.4) ##plot data points with specifed background 
##points(comm.nms, pch=points, cex=2) ## enter points
##
legend('topleft', c('Forest', 'Secondary Forest','Pasture'), pch=c(17, 22, 1), cex=.75) ## add ledgends on the top left (though not correct) cex is for the size of ledgend



###Fig2 # 95%_DNA adding multiple figures

##individual nmds of same size
library(vegan)

trt = read.csv('treatment_file3.csv') ## download data file
dat2 = read.csv('95%_FSP_S1_DNA.csv')
##set as data frame, transpose, data matrix and set as  a numaric varible

comm2 = as.data.frame(t(dat2)) ## Make the test13 as data frame
names(comm2) = as.matrix(comm2[1,  ]) ## change the data as matrix
rownames_comm2 = as.matrix(rownames(comm2)) ##change the row names as matrix
comm2 = comm2[-1, ] ## ----except first row
comm2 = apply(comm2, 2, as.numeric) ## set all data as numeric starting from row 2 ???
rownames(comm2) = rownames_comm2[-1,] ## Adjust the row names
#head(comm)
trt = subset(trt, subset = trt$sample %in% rownames(comm2))
## check that rows of comm so that it matches coords
all(trt$sample == rownames(comm2))
## if not true then select the same rows that match with the treatment files... subset coords to only samples we have community data for
comm2 = subset(comm2, subset=trt$sample %in% rownames(comm2))
## set as a factors
## run ANOSIM analysis 
anosim (comm, grouping=trt$treatment, permutations=999, distance="bray") ## did the anosim indirect ordination data (Bray-curtis)

##metanMDS analysis 
comm.nms2 = metaMDS(comm2) ##metaMDS analysis

par(mfrow=c(2,2))



plot(comm.nms2, display='sites', type='n', ylim=c(-2, 2), xlim=c(-2.2, 2.2)) ##plot bank box
points= c(rep(17, 31), rep(22, 29), rep(1, 30)) ##repeat the selected points (3 time each)
bgs= c(rep("Black", 31), rep("grey70", 29), rep("white", 30)) 
points(comm.nms2, pch=points, bg=bgs, cex=1) ##plot data points with specifed background 
##points(comm.nms, pch=points, cex=1) ## enter points
##
legend('topleft', c('Forest', 'Secondary Forest','Pasture'), pch=c(17, 22, 1), cex=.75) ## add ledgends on the top left (though not correct) cex is for the size of ledgend


###Fig3 90%_DNA #adding multiple figures

##individual nmds of same size
library(vegan)

trt = read.csv('treatment_file3.csv') ## download data file
dat3 = read.csv('90%_FSP_S1_DNA.csv')
##set as data frame, transpose, data matrix and set as  a numaric varible

comm3 = as.data.frame(t(dat3)) ## Make the test13 as data frame
names(comm3) = as.matrix(comm3[1,  ]) ## change the data as matrix
rownames_comm3 = as.matrix(rownames(comm3)) ##change the row names as matrix
comm3 = comm3[-1, ] ## ----except first row
comm3 = apply(comm3, 2, as.numeric) ## set all data as numeric starting from row 2 ???
rownames(comm3) = rownames_comm3[-1,] ## Adjust the row names
#head(comm)
trt = subset(trt, subset = trt$sample %in% rownames(comm3))
## check that rows of comm so that it matches coords
all(trt$sample == rownames(comm3))
## if not true then select the same rows that match with the treatment files... subset coords to only samples we have community data for
comm3 = subset(comm3, subset=trt$sample %in% rownames(comm3))
## set as a factors
## run ANOSIM analysis 
anosim (comm, grouping=trt$treatment, permutations=999, distance="bray") ## did the anosim indirect ordination data (Bray-curtis)

##metanMDS analysis 
comm.nms3 = metaMDS(comm3) ##metaMDS analysis

par(mfrow=c(2,3))

plot(comm.nms3, display='sites', type='n', xlim=c(-2, 2), ylim=c(-2,2),cex=2 ) ##plot bank box
points= c(rep(17, 31), rep(22, 29), rep(1, 30)) ##repeat the selected points (3 time each)
bgs= c(rep("Black", 31), rep("grey70", 29), rep("white", 30)) 
points(comm.nms3, pch=points, bg=bgs, cex=1) ##plot data points with specifed background 
##points(comm.nms, pch=points, cex=2) ## enter points
##
legend('topleft', c('Forest', 'Secondary Forest','Pasture'), pch=c(17, 22, 1), cex=.75) ## add ledgends on the top left (though not correct) cex is for the size of ledgend


###Fig4 99%_prot
#adding multiple figures

##individual nmds of same size
library(vegan)

trt = read.csv('treatment_file3.csv') ## download data file
dat4 = read.csv('99%_FSP_S1_prot.csv')
##set as data frame, transpose, data matrix and set as  a numaric varible

comm4 = as.data.frame(t(dat4)) ## Make the test13 as data frame
names(comm4) = as.matrix(comm4[1,  ]) ## change the data as matrix
rownames_comm4 = as.matrix(rownames(comm4)) ##change the row names as matrix
comm4 = comm4[-1, ] ## ----except first row
comm4 = apply(comm4, 2, as.numeric) ## set all data as numeric starting from row 2 ???
rownames(comm4) = rownames_comm4[-1,] ## Adjust the row names
#head(comm)
trt = subset(trt, subset = trt$sample %in% rownames(comm4))
## check that rows of comm so that it matches coords
all(trt$sample == rownames(comm4))
## if not true then select the same rows that match with the treatment files... subset coords to only samples we have community data for
comm4 = subset(comm4, subset=trt$sample %in% rownames(comm4))
## set as a factors
## run ANOSIM analysis 
anosim (comm, grouping=trt$treatment, permutations=999, distance="bray") ## did the anosim indirect ordination data (Bray-curtis)

##metanMDS analysis 
comm.nms4 = metaMDS(comm4) ##metaMDS analysis

par(mfrow=c(2,3))

plot(comm.nms4, display='sites', type='n', xlim=c(-2, 2), ylim=c(-2,2),cex=2 ) ##plot bank box
points= c(rep(17, 31), rep(22, 29), rep(1, 30)) ##repeat the selected points (3 time each)
bgs= c(rep("Black", 31), rep("grey70", 29), rep("white", 30)) 
points(comm.nms4, pch=points, bg=bgs, cex=1.4) ##plot data points with specifed background 
##points(comm.nms, pch=points, cex=2) ## enter points
##
legend('topleft', c('Forest', 'Secondary Forest','Pasture'), pch=c(17, 22, 1), cex=.75) ## add ledgends on the top left (though not correct) cex is for the size of ledgend



##Fig 5 97% prot
#adding multiple figures

##individual nmds of same size
library(vegan)

trt = read.csv('treatment_file3.csv') ## download data file
dat5 = read.csv('97%_FSP_S1_prot.csv')
##set as data frame, transpose, data matrix and set as  a numaric varible

comm5 = as.data.frame(t(dat5)) ## Make the test13 as data frame
names(comm5) = as.matrix(comm5[1,  ]) ## change the data as matrix
rownames_comm5 = as.matrix(rownames(comm5)) ##change the row names as matrix
comm5 = comm5[-1, ] ## ----except first row
comm5 = apply(comm5, 2, as.numeric) ## set all data as numeric starting from row 2 ???
rownames(comm5) = rownames_comm5[-1,] ## Adjust the row names
#head(comm)
trt = subset(trt, subset = trt$sample %in% rownames(comm5))
## check that rows of comm so that it matches coords
all(trt$sample == rownames(comm5))
## if not true then select the same rows that match with the treatment files... subset coords to only samples we have community data for
comm5 = subset(comm5, subset=trt$sample %in% rownames(comm5))
## set as a factors
## run ANOSIM analysis 
anosim (comm, grouping=trt$treatment, permutations=999, distance="bray") ## did the anosim indirect ordination data (Bray-curtis)

##metanMDS analysis 
comm.nms5 = metaMDS(comm5) ##metaMDS analysis

par(mfrow=c(2,2))

plot(comm.nms5, display='sites', type='n', xlim=c(-2, 2), ylim=c(-2,2),cex=2 ) ##plot bank box
points= c(rep(17, 31), rep(22, 29), rep(1, 30)) ##repeat the selected points (3 time each)
bgs= c(rep("Black", 31), rep("grey70", 29), rep("white", 30)) 
points(comm.nms5, pch=points, bg=bgs, cex=1) ##plot data points with specifed background 
##points(comm.nms, pch=points, cex=2) ## enter points
##
legend('topleft', c('Forest', 'Secondary Forest','Pasture'), pch=c(17, 22, 1), cex=.75) ## add ledgends on the top left (though not correct) cex is for the size of ledgend

##Fig 6 95% prot
#adding multiple figures

##individual nmds of same size
library(vegan)

trt = read.csv('treatment_file3.csv') ## download data file
dat6 = read.csv('95%_FSP_S1_prot.csv')
##set as data frame, transpose, data matrix and set as  a numaric varible

comm6 = as.data.frame(t(dat6)) ## Make the test13 as data frame
names(comm6) = as.matrix(comm6[1,  ]) ## change the data as matrix
rownames_comm6 = as.matrix(rownames(comm6)) ##change the row names as matrix
comm6 = comm6[-1, ] ## ----except first row
comm6 = apply(comm6, 2, as.numeric) ## set all data as numeric starting from row 2 ???
rownames(comm6) = rownames_comm6[-1,] ## Adjust the row names
#head(comm)
trt = subset(trt, subset = trt$sample %in% rownames(comm6))
## check that rows of comm so that it matches coords
all(trt$sample == rownames(comm6))
## if not true then select the same rows that match with the treatment files... subset coords to only samples we have community data for
comm6 = subset(comm6, subset=trt$sample %in% rownames(comm6))
## set as a factors
## run ANOSIM analysis 
anosim (comm, grouping=trt$treatment, permutations=999, distance="bray") ## did the anosim indirect ordination data (Bray-curtis)

##metanMDS analysis 
comm.nms6 = metaMDS(comm6) ##metaMDS analysis

par(mfrow=c(2,2))

plot(comm.nms6, display='sites', type='n', xlim=c(-2, 2), ylim=c(-2,2),cex=2 ) ##plot bank box
points= c(rep(17, 31), rep(22, 29), rep(1, 30)) ##repeat the selected points (3 time each)
bgs= c(rep("Black", 31), rep("grey70", 29), rep("white", 30)) 
points(comm.nms6, pch=points, bg=bgs, cex=1) ##plot data points with specifed background 
##points(comm.nms, pch=points, cex=2) ## enter points
##
legend('topleft', c('Forest', 'Secondary Forest','Pasture'), pch=c(17, 22, 1), cex=.75) ## add ledgends on the top left (though not correct) cex is for the size of ledgend


##Fig 6 95% prot
#adding multiple figures

##individual nmds of same size
library(vegan)

trt = read.csv('treatment_file3.csv') ## download data file
dat7 = read.csv('90%_FSP_S1_prot.csv')
##set as data frame, transpose, data matrix and set as  a numaric varible

comm7 = as.data.frame(t(dat7)) ## Make the test13 as data frame
names(comm7) = as.matrix(comm7[1,  ]) ## change the data as matrix
rownames_comm7 = as.matrix(rownames(comm7)) ##change the row names as matrix
comm7 = comm7[-1, ] ## ----except first row
comm7 = apply(comm7, 2, as.numeric) ## set all data as numeric starting from row 2 ???
rownames(comm7) = rownames_comm7[-1,] ## Adjust the row names
#head(comm)
trt = subset(trt, subset = trt$sample %in% rownames(comm7))
## check that rows of comm so that it matches coords
all(trt$sample == rownames(comm7))
## if not true then select the same rows that match with the treatment files... subset coords to only samples we have community data for
comm7 = subset(comm7, subset=trt$sample %in% rownames(comm6))
## set as a factors
## run ANOSIM analysis 
anosim (comm, grouping=trt$treatment, permutations=999, distance="bray") ## did the anosim indirect ordination data (Bray-curtis)

##metanMDS analysis 
comm.nms7 = metaMDS(comm7) ##metaMDS analysis

par(mfrow=c(2,3))

plot(comm.nms, display='sites', type='n', xlim=c(-2, 2), ylim=c(-2,2),cex=2 ) ##plot bank box
points= c(rep(17, 31), rep(22, 29), rep(1, 30)) ##repeat the selected points (3 time each)
bgs= c(rep("Black", 31), rep("grey70", 29), rep("white", 30)) 
points(comm.nms, pch=points, bg=bgs, cex=2) ##plot data points with specifed background 
##points(comm.nms, pch=points, cex=2) ## enter points
##
legend('topleft', c('Forest', 'Secondary Forest','Pasture'), pch=c(17, 22, 1), cex=.75) ## add ledgends on the top left (though not correct) cex is for the size of ledgend






par(mfrow=c(2,2))


plot(comm.nms2, display='sites', type='n', xlim=c(-2.5, 2.5), ylim=c(-2,2.5),cex=2 ) ##plot bank box
points= c(rep(17, 31), rep(22, 29), rep(1, 30)) ##repeat the selected points (3 time each)
bgs= c(rep("Black", 31), rep("grey70", 29), rep("white", 30)) 
points(comm.nms2, pch=points, bg=bgs, cex=2) ##plot data points with specifed background 

plot(comm.nms3, display='sites', type='n', xlim=c(-2, 2.5), ylim=c(-2,2.5),cex=2 ) ##plot bank box
points= c(rep(17, 31), rep(22, 29), rep(1, 30)) ##repeat the selected points (3 time each)
bgs= c(rep("Black", 31), rep("grey70", 29), rep("white", 30)) 
points(comm.nms3, pch=points, bg=bgs, cex=2) ##plot data points with specifed background 

##plot(comm.nms4, display='sites', type='n', xlim=c(-2, 2), ylim=c(-2,2),cex=2 ) ##plot bank box
##points= c(rep(17, 31), rep(22, 29), rep(1, 30)) ##repeat the selected points (3 time each)
##bgs= c(rep("Black", 31), rep("grey70", 29), rep("white", 30)) 
###points(comm.nms4, pch=points, bg=bgs, cex=2) ##plot data points with specifed background 

plot(comm.nms5, display='sites', type='n', xlim=c(-2, 2.5), ylim=c(-2,2.5),cex=2 ) ##plot bank box
points= c(rep(17, 31), rep(22, 29), rep(1, 30)) ##repeat the selected points (3 time each)
bgs= c(rep("Black", 31), rep("grey70", 29), rep("white", 30)) 
points(comm.nms5, pch=points, bg=bgs, cex=2) ##plot data points with specifed background 


plot(comm.nms6, display='sites', type='n', xlim=c(-2.5, 2,5), ylim=c(-2, 2),cex=2 ) ##plot bank box
points= c(rep(17, 31), rep(22, 29), rep(1, 30)) ##repeat the selected points (3 time each)
bgs= c(rep("Black", 31), rep("grey70", 29), rep("white", 30)) 
points(comm.nms6, pch=points, bg=bgs, cex=2) ##plot data points with specifed background 

##plot(comm.nms7, display='sites', type='n', xlim=c(-2, 2), ylim=c(-2,2),cex=2 ) ##plot bank box
##points= c(rep(17, 31), rep(22, 29), rep(1, 30)) ##repeat the selected points (3 time each)
##bgs= c(rep("Black", 31), rep("grey70", 29), rep("white", 30)) 
##points(comm.nms7, pch=points, bg=bgs, cex=2) ##plot data points with specifed background 





####


