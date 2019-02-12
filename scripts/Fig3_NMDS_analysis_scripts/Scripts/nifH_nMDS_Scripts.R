### step 1 processing the OTU for 454 arrA data for treatment differences for rda, ANOSIM, nMDS analysis


nifH_90_DNA = read.csv('90%_FSP_S1_DNA.csv') ## download data file

trt = read.csv('treatment_file2.csv') ## download data file
##this step is not required for this data. it is done to sum the rows with identical row names
for (i in 2:ncol(nifH_90_DNA)) {
  if (i == 2)
    nifH_90_DNA_aggr = aggregate(nifH_90_DNA[ , i], by=list(as.character(nifH_90_DNA[,1])), FUN='sum')
  else
    nifH_90_DNA_aggr = cbind(nifH_90_DNA_aggr, aggregate(nifH_90_DNA[ , i], by=list(as.character(nifH_90_DNA[,1])), FUN='sum')[,2]) 
} 

names(nifH_90_DNA_aggr) = names(nifH_90_DNA) ## change the data file names

nifH_90_DNA2 = nifH_90_DNA_aggr[rowSums(nifH_90_DNA_aggr[,-1]) > 0 , ] ## Select row having more than 0 values

nifH_90_DNA3 = as.data.frame(t(nifH_90_DNA2)) ## Make the test13 as data frame

names(nifH_90_DNA3) 

names(nifH_90_DNA3) = as.matrix(nifH_90_DNA3[1,  ]) ## change the data as matrix

rownames_nifH_90_DNA4 = as.matrix(rownames(nifH_90_DNA3)) ##change the row names as matrix

nifH_90_DNA3 = nifH_90_DNA3[-1, ] ## ----except first row

nifH_90_DNA3 = apply(nifH_90_DNA3, 2, as.numeric) ## set all data as numeric starting from row 2 ???

rownames(nifH_90_DNA3) = rownames_nifH_90_DNA4[-1,] ## Adjust the row names

head(nifH_90_DNA3)

anosim(nifH_90_DNA3, grouping=trt$treatment) ## did the anosim indirect ordination data (Bray-curtis)

rda.nifH_90_DNA3 = rda(nifH_90_DNA3 ~ trt$treatment) ## Did RDA using the direct ordination
rda.nifH_90_DNA3.anv = anova(rda.nifH_90_DNA3) ## RDA anova to see if the differences are significant or not
rda.nifH_90_DNA3.anv
rda.nifH_90_DNA3.anv$Var[1] / sum(rda.nifH_90_DNA3.anv$Var)  ## total variance explained by rda

nifH_90_DNA3.nms = metaMDS(nifH_90_DNA3) ##metaMDS analysis 
plot(nifH_90_DNA3.nms) ##plot MDS
plot(test3.nms, type='t') ## plot with sites
##pretty graphs
##shapes and color ledgened

plot(nifH_90_DNA3.nms, display='sites', type='n') ##plot bank box
points= c(rep(0, 30), rep(15, 29), rep(19, 30) ) ##repeat the selected points (3 time each)

points(nifH_90_DNA3.nms, pch=points) ## enter points
cols= c(rep("#66F0FF", 30), rep("#99F8FF", 29), rep("#FFCA99", 30)) ## made color files

points(nifH_90_DNA3.nms, pch=points, col=cols)
legend('lowleft', c('NP1T0', 'NP1G','NPGG','NP3G','NP3GG','NP8T0','NP8G', 'NP8GG'),
       pch=points[c(0, 15, 19, 17, 1, 2, 23, 5)], bty='n', col=cols,
       cex=.75) ## add ledgends on the top left (though not correct) cex is for the size of ledgend

plot(nifH_90_DNA3.nms, display='sites', type='t') ##plot with sites name box

plot(nifH_90_DNA3.nms, grouping='trt$treatment') ## plot with the reatment wise










### step 2 processing the OTU for 454 arrA data for treatment differences for rda, ANOSIM, nMDS analysis


nifH_95_DNA = read.csv('95%_FSP_S1_DNA.csv') ## download data file

trt = read.csv('treatment_file3.csv') ## download data file
##this step is not required for this data. it is done to sum the rows with identical row names
for (i in 2:ncol(nifH_95_DNA)) {
  if (i == 2)
    nifH_95_DNA_aggr = aggregate(nifH_95_DNA[ , i], by=list(as.character(nifH_95_DNA[,1])), FUN='sum')
  else
    nifH_95_DNA_aggr = cbind(nifH_95_DNA_aggr, aggregate(nifH_95_DNA[ , i], by=list(as.character(nifH_95_DNA[,1])), FUN='sum')[,2]) 
} 

names(nifH_95_DNA_aggr) = names(nifH_95_DNA) ## change the data file names

nifH_95_DNA2 = nifH_95_DNA_aggr[rowSums(nifH_95_DNA_aggr[,-1]) > 0 , ] ## Select row having more than 0 values

nifH_95_DNA3 = as.data.frame(t(nifH_95_DNA2)) ## Make the test13 as data frame

names(nifH_95_DNA3) 

names(nifH_95_DNA3) = as.matrix(nifH_95_DNA3[1,  ]) ## change the data as matrix

rownames_nifH_95_DNA4 = as.matrix(rownames(nifH_95_DNA3)) ##change the row names as matrix

nifH_95_DNA3 = nifH_95_DNA3[-1, ] ## ----except first row

nifH_95_DNA3 = apply(nifH_95_DNA3, 2, as.numeric) ## set all data as numeric starting from row 2 ???

rownames(nifH_95_DNA3) = rownames_nifH_95_DNA4[-1,] ## Adjust the row names

head(nifH_95_DNA3)

anosim(nifH_95_DNA3, grouping=trt$treatment) ## did the anosim indirect ordination data (Bray-curtis)

rda.nifH_95_DNA3 = rda(nifH_95_DNA3 ~ trt$treatment) ## Did RDA using the direct ordination
rda.nifH_95_DNA3.anv = anova(rda.nifH_95_DNA3) ## RDA anova to see if the differences are significant or not
rda.nifH_95_DNA3.anv
rda.nifH_95_DNA3.anv$Var[1] / sum(rda.nifH_95_DNA3.anv$Var)  ## total variance explained by rda

nifH_95_DNA3.nms = metaMDS(nifH_95_DNA3) ##metaMDS analysis 
plot(nifH_95_DNA3.nms) ##plot MDS
plot(test3.nms, type='t') ## plot with sites
##pretty graphs
##shapes and color ledgened

plot(nifH_95_DNA3.nms, display='sites', type='n') ##plot bank box
points= c(rep(0, 31), rep(15, 29), rep(19, 30) ) ##repeat the selected points (3 time each)

points(nifH_95_DNA3.nms, pch=points) ## enter points
cols= c(rep("#66F0FF", 31), rep("#99F8FF", 29), rep("#FFCA99", 30)) ## made color files

points(nifH_95_DNA3.nms, pch=points, col=cols)
legend('lowleft', c('Forest', 'SF','Pasture')),
       pch=points[c(0, 15, 19)], bty='n', col=cols,
       cex=.75) ## add ledgends on the top left (though not correct) cex is for the size of ledgend

plot(nifH_95_DNA3.nms, display='sites', type='t') ##plot with sites name box

plot(nifH_95_DNA3.nms, grouping='trt$treatment') ## plot with the reatment wise











### step 3 processing the OTU for 454 arrA data for treatment differences for rda, ANOSIM, nMDS analysis


nifH_97_DNA = read.csv('97%_FSP_S1_DNA.csv') ## download data file

trt = read.csv('treatment_file3.csv') ## download data file
##this step is not required for this data. it is done to sum the rows with identical row names
for (i in 2:ncol(nifH_97_DNA)) {
  if (i == 2)
    nifH_97_DNA_aggr = aggregate(nifH_97_DNA[ , i], by=list(as.character(nifH_97_DNA[,1])), FUN='sum')
  else
    nifH_97_DNA_aggr = cbind(nifH_97_DNA_aggr, aggregate(nifH_97_DNA[ , i], by=list(as.character(nifH_97_DNA[,1])), FUN='sum')[,2]) 
} 

names(nifH_97_DNA_aggr) = names(nifH_97_DNA) ## change the data file names

nifH_97_DNA2 = nifH_97_DNA_aggr[rowSums(nifH_97_DNA_aggr[,-1]) > 0 , ] ## Select row having more than 0 values

nifH_97_DNA3 = as.data.frame(t(nifH_97_DNA2)) ## Make the test13 as data frame

names(nifH_97_DNA3) 

names(nifH_97_DNA3) = as.matrix(nifH_97_DNA3[1,  ]) ## change the data as matrix

rownames_nifH_97_DNA4 = as.matrix(rownames(nifH_97_DNA3)) ##change the row names as matrix

nifH_97_DNA3 = nifH_97_DNA3[-1, ] ## ----except first row

nifH_97_DNA3 = apply(nifH_97_DNA3, 2, as.numeric) ## set all data as numeric starting from row 2 ???

rownames(nifH_97_DNA3) = rownames_nifH_97_DNA4[-1,] ## Adjust the row names

head(nifH_97_DNA3)

anosim(nifH_97_DNA3, grouping=trt$treatment) ## did the anosim indirect ordination data (Bray-curtis)

rda.nifH_97_DNA3 = rda(nifH_97_DNA3 ~ trt$treatment) ## Did RDA using the direct ordination
rda.nifH_97_DNA3.anv = anova(rda.nifH_97_DNA3) ## RDA anova to see if the differences are significant or not
rda.nifH_97_DNA3.anv
rda.nifH_97_DNA3.anv$Var[1] / sum(rda.nifH_97_DNA3.anv$Var)  ## total variance explained by rda

nifH_97_DNA3.nms = metaMDS(nifH_97_DNA3) ##metaMDS analysis 
plot(nifH_97_DNA3.nms) ##plot MDS
plot(test3.nms, type='t') ## plot with sites
##pretty graphs
##shapes and color ledgened

plot(nifH_97_DNA3.nms, display='sites', type='n') ##plot bank box
points= c(rep(0, 31), rep(15, 29), rep(19, 30) ) ##repeat the selected points (3 time each)

points(nifH_97_DNA3.nms, pch=points) ## enter points
cols= c(rep("#66F0FF", 31), rep("#99F8FF", 29), rep("#FFCA99", 30)) ## made color files

points(nifH_97_DNA3.nms, pch=points, col=cols)
legend('lowleft', c('Forest', 'SF','Pasture')),
pch=points[c(0, 15, 19)], bty='n', col=cols,
cex=.75) ## add ledgends on the top left (though not correct) cex is for the size of ledgend

plot(nifH_97_DNA3.nms, display='sites', type='t') ##plot with sites name box

plot(nifH_97_DNA3.nms, grouping='trt$treatment') ## plot with the reatment wise
