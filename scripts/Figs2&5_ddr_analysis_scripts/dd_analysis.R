## to analyze distance decay in Forest and Pasture

## load package and necessary functions ------------
library('vegan')
source('./Dist_Decay_analysis/spat_functions.R')

## get abundance file names -----------------------
files = dir('./Dist_Decay_analysis/abu_data') ## get list of files
comm_files = files[grep('FSP', files)] ## pull out file names with FSP

read.unifrac = function(file, sample_ids=NULL) {
  ## reads in unifrac distance matrix and optionally subsets 
  ## arguments
  ## file: file name with path to distance matrix
  ## sample_ids: which samples to keep
  unifrac_dat = read.table(file)
  if (!is.null(sample_ids)) {
    indices = match(sample_ids, rownames(unifrac_dat))
    indices = indices[!is.na(indices)]
    unifrac_dat = unifrac_dat[indices, indices]
  }  
  unifrac_dist = as.dist(unifrac_dat)
  return(list(dist=unifrac_dist, samples=rownames(unifrac_dat)))
}

## taxonomic and unifrac DDR ------------------------------------------------------------------------

calc_ddr = function(x, coords, metric, breaks, hmin=NA, hmax=NA) {
  ## computes the distance decay relationship
  ## Returns:
  ## list with first element the raw DDR values, second element the averaged values
  ## Arguments:
  ## x: either a site x sp matrix or a distance matrix
  ## coords: x y spatial coordinates
  ## metric: what similarity index to use
  ## breaks: the spatial breaks to average at
  ## hmin: minimum distance to examine, can be left blank
  ## hmax: maximum distance to examine, can be left blank
  if (class(x) == 'matrix')
    x = vegdist(x, method=metric)
  rawDDR = data.frame(gdist = as.vector(dist(coords)), 
                      csiml = as.vector(1 - x))
  avgDDR = vario(x, coords, breaks=breaks,
                 distance.metric='unifrac', hmin=hmin, hmax=hmax)$vario
  return(list(rawDDR, avgDDR))
}

load_ddr_data = function(file, metric) {
  ## extract marker and percent info from file
  split_name = strsplit(file, '_')             ## split up the file name into its components
  percent = sub('%', '', split_name[[1]][1])   ## identify the % similarity
  marker = sub('.csv', '', split_name[[1]][4]) ## identify the marker
  
  ## import data
  dat = read.csv(file.path('./Dist_Decay_analysis/abu_data/', file))
  
  ## reshape community data for analysis
  comm = t(dat[ , -1])
  colnames(comm) = as.character(dat[ , 1])
  
  coords = read.csv('./Dist_Decay_analysis/FNL_site_dat_coords.csv')
  ## subset coords to only samples we have community data for
  coords = subset(coords, subset=sample %in% rownames(comm))
  
  ## drop secondary forest
  #comm = comm[coords$systems != 'Secondary Forest', ]
  #coords = coords[coords$systems != 'Secondary Forest', ]
  
  ## check that rows of comm so that it matches coords
  #all(coords$sample == rownames(comm))
  
  habitat = unique(coords$systems)
  
  ## compute DDR for plotting----------------------------------------------
  ## compute the average dissimalarity for given spatial breaks 
  if (metric == 'unifrac') {
    unifrac_file = paste('./Dist_Decay_analysis/unifrac_dist/unifrac_dist', percent,
                         marker, 'abu.txt', sep='_') 
    unifrac = read.unifrac(unifrac_file)
    bact_dist = unifrac$dist
  }
  else
    bact_dist = vegdist(comm, method='bray')
  breaks = c(0.01, 0.5, 2, 25, 500, 1500, 15000)
  rawDDR = avgDDR = vector('list', length(habitat))
  names(rawDDR) = names(avgDDR) = habitat
  for (i in seq_along(habitat)) {
    true = coords$systems == habitat[i]                                        ## subset to only habitat of interest
    crds = coords[true, c('x','y')]                                            ## assign coordinates
    if (metric == 'unifrac') {
      unifrac = read.unifrac(unifrac_file, sample_ids= coords$sample[true])
      crds = crds[match(unifrac$samples, coords[true, 'sample']), ]
      cdist = unifrac$dist
    }
    else
      cdist = vegdist(comm[true, ], method = 'bray')
    DDR = calc_ddr(cdist, crds, metric, breaks, hmin=0, hmax=20000)
    rawDDR[[i]] = DDR[[1]]
    avgDDR[[i]] = DDR[[2]]
  }
  out = list(marker=marker, percent=percent, metric=metric, 
             comm=comm, bact_dist=bact_dist, coords=coords[ ,c('x', 'y')],
             habitat=coords$systems, rawDDR=rawDDR, avgDDR=avgDDR)
  return(out)
}

## compare statistical fits-----------------------------------------------------
compare_mods = function(ddr_data) {
  ## compares the r2 values for the linear, exponential, and pwr models
  ## of the DDR within each habitat type
  habitat_types = unique(ddr_data$habitat)
  r2 = matrix(NA, ncol=3, nrow=3)
  colnames(r2) = habitat_types
  rownames(r2) = c('lin','exp','pwr')
  for (i in seq_along(habitat_types)) {
    lin_fit = lm(csiml ~ gdist, data=ddr_data$rawDDR[[i]])
    exp_fit = lm(log(csiml) ~ gdist, data=ddr_data$rawDDR[[i]])
    pwr_fit = lm(log(csiml) ~ log(gdist), data=ddr_data$rawDDR[[i]])
    r2[ , i] = sapply(list(lin_fit, exp_fit, pwr_fit), function(x) summary(x)$r.squared)
  }
  return(r2)
}

bray_mods = lapply(comm_files, function(x) compare_mods(load_ddr_data(x, 'bray')))
unifrac_mods = lapply(comm_files, function(x) compare_mods(load_ddr_data(x, 'unifrac')))

bray_avg = unifrac_avg = matrix(0, 3, 3)
for(i in seq_along(bray_mods)) {
  bray_avg = bray_avg + bray_mods[[i]] / 6
  unifrac_avg = unifrac_avg + unifrac_mods[[i]] / 6
}

bray_avg
unifrac_avg
apply(bray_avg, 1, mean)
apply(unifrac_avg, 1, mean)
# take home is that the exponential tends to always be best choice

## make figures with exponential model fits-------------------------------------

make_figs = function(ddr_data, method, mod=NULL, add_avg=F,
                     path='./Dist_Decay_analysis/figs/') {
  if (method == 'threepanel') {
    fig_path = paste(path, ddr_data$metric, '_', ddr_data$percent, '_',
                     ddr_data$marker, '_DD_threepanel.png', sep='')
    png(file=fig_path, width=480 * 3)
    habitat_types = unique(ddr_data$habitat)
    if (ddr_data$metric == 'bray')
      ylim = c(.001, 1)
    else
      ylim = range(sapply(ddr_data$rawDDR, function(x) range(x$csiml)))
    xlim = range(sapply(ddr_data$rawDDR, function(x) range(x$gdist)))
    #(bottom, left, top, right)
    # c(5, 4, 4, 2) + 0.1
    par(mar=c(5, 5, 4, 1) + 0.1)
    par(mfrow=c(1,3))
    for (j in c('Forest', 'Secondary Forest', 'Pasture')) {
      i = match(j, habitat_types)
      ## assign coordinates
      ## set graphical x and y limits
      #xlim = c(0.001, 80000)
      #ylim = range(sapply(ddr_data$rawDDR, function(x) range(x[,2], na.rm=T)))
      ## plot result on log-log axes and convert dissimarlity to similarity
      if (ddr_data$metric == 'bray')
        metric_type = 'Bray'
      else
        metric_type = 'Unifrac'
      plot(csiml ~ gdist, data=ddr_data$rawDDR[[i]],log='xy', ylim=ylim, xlim=xlim,
           xlab='', ylab='', frame.plot=F, axes=F)
      axis(side=1, cex.axis=1.75, padj=.25, lwd=4)
      if (ddr_data$metric == 'bray')
        axis(side=2, cex.axis=1.75, padj= 0.1, lwd=4, at=c(.001, .01, .1, 1))
      else 
        axis(side=2, cex.axis=1.75, padj= 0.1, lwd=4)
      mtext(side=1, 'Geographic Distance (m)', cex=1.6, padj=2.25)
      mtext(side=2, paste('Community Similarity (', metric_type, ')', sep=''),
            cex=1.5, padj=-1.9)
      mtext(side=3, habitat_types[i], cex=2)
      ## add the averaged result
      if(!is.null(mod)) {
        if (mod == 'lin') {
          fit = lm(csiml ~ gdist, data=ddr_data$rawDDR[[i]])
          pred = predict(fit)
        }
        else if(mod == 'exp') {
          fit = lm(log(csiml) ~ gdist, data=ddr_data$rawDDR[[i]])
          pred = exp(predict(fit))
        }
        else if(mod == 'pwr') {
          fit = lm(log(csiml) ~ log(gdist), data=ddr_data$rawDDR[[i]])
          pred = exp(predict(fit))
        }
        gdist = ddr_data$rawDDR[[i]]$gdist
        lines(gdist[order(gdist)], pred[order(gdist)], lwd=3)
        r2 = round(summary(fit)$r.squared, 2)
        b0 = round(exp(coef(fit)[1]), 3)
        b1 = round(coef(fit)[2], 3)
        if (mod == 'pwr')
          mod_type = 'Power Model'
        if (mod == 'exp')
          mod_type = 'Exponential Model'
        if (mod == 'lin')
          mod_type = 'Linear Model'
        leg_txt = c(mod_type, paste('log(S) = ', b0, ' - ', abs(b1), '*log(D)', sep=''),
                    paste('R^2 =', r2[1]))
        legend('bottomleft', leg_txt, bty='n', cex=2)
      }
      if(add_avg)
        lines(I(1 - exp.var) ~ Dist, data=ddr_data$avgDDR[[i]], cex=1.5, lwd=5,
              type='o')
    }
    dev.off()
  }
}

lapply(comm_files, function(x)
                   make_figs(load_ddr_data(x, 'bray'), 'threepanel', 'pwr'))
lapply(comm_files, function(x) 
                   make_figs(load_ddr_data(x, 'unifrac'), 'threepanel', 'pwr'))

## ADONIS analysis ------------------------------------------------------------
## test if distance decay relationships differ between habitat types
if (metric == 'unifrac') 
  cdist = read.unifrac(unifrac_file, sample_ids=coords$samp)
if (metric == 'bray')
  cdist = vegdist(comm, method='bray') 
xy = as.matrix(coords[, c('x', 'y')])
gdist = dist(xy)
habitat = ifelse(coords$system == 'Pasture', 0, 1)
## Carry out permutation test for difference in avg similarity
print(paste('metric = ', metric, ', percent = ', percent,
            ', and marker = ', marker, sep=''))
## test if there is a difference in beta diversity between
## habitat types
## TO DO: write model outputs to text file
print(adonis(cdist ~ habitat))
## test if there is a relationship between beta diversity 
## spatial proximity, habitat, and their interaction
print(adonis(cdist ~ xy + habitat + xy * habitat))

## soil analysis --------------------------------------------------------------
soil = read.csv('./Env_data/FSP_Environmental_data.csv')
coords = read.csv('./Dist_Decay_analysis/FNL_site_dat_coords.csv')
sum(soil$sample %in% coords$sample)
missing_soil = coords$sample[!(coords$sample %in% soil$sample)]
coords[!(coords$sample %in% soil$sample), ]
#fill in missing soil samples
soil = read.csv('./Env_data//FSP_Environmental_data_missing_samples_added.csv')
rownames(soil) = soil[,1]
soil = soil[,-1]

## variation partitioning ----------------------------------------------------
pca = princomp(scale(soil))

plot(pca)
biplot(pca)
pca$loadings

## pca does not provide much resolution on the soil variables
## it does indicate that pastures tend to be higher in nutrients and organic matter

library(cluster)

for(i in seq_along(files)) {

  dat = load_ddr_data(files[i], 'unifrac')
  bact_dist = dat$bact_dist
  samples = rownames(as.matrix(bact_dist))
  
  land_dist = daisy(data.frame(habitat=dat$habitat))
  geog_dist = dist(dat$coords)
  
  ## check samples of soil matrix
  soil_tmp = soil[match(samples, rownames(soil)),]
  soil_dist = dist(scale(soil_tmp))
  
  ## convert to vectors
  bact_dist = as.vector(bact_dist)
  land_dist = as.vector(land_dist)
  soil_dist = as.vector(soil_dist)
  geog_dist = as.vector(geog_dist)
  varpart(bact_dist, land_dist, soil_dist, geog_dist)
  
  anova(rda(bact_dist, land_dist, cbind(soil_dist, geog_dist)))
  anova(rda(bact_dist, soil_dist, cbind(land_dist, geog_dist)))
  anova(rda(bact_dist, geog_dist, cbind(land_dist, soil_dist)))
}



## old code below -------------------------------------------------------------

  
    ## overlay pasture and forest
    fig_path = paste('./Dist_Decay_analysis/figs/', metric, '_', percent, '_', 
                     marker, '_DD_onepanel.png', sep='')
    png(file=fig_path, width=480)
    par(mfrow=c(1,1))
    col = c('darkgreen', 'red')
    for (i in seq_along(habitat)) {
      xlim = c(0.01, 80000)
      ylim = range(sapply(rawDDR, function(x) range(x[,2], na.rm=T)))
      if (i == 1) {
        plot(csiml ~ gdist, data=rawDDR[[i]] , log='xy', ylim=ylim, xlim=xlim,
             col=col[i], type='n', xlab='',ylab='', axes=F)
        addAxes()
        mtext(side=1, 'Geographic Distance (m)', cex=2, padj=2.5)
        mtext(side=2, 'Community Similarity (Bray-Curtis)', cex=2, padj=-2)
        legend('topright', c('Forest', 'Pasture'), col=col, lwd=5, bty='n', cex=2)
      }  
      points(I(1 - exp.var) ~ Dist, data=avgDDR[[i]], col=col[i], pch=19,
             cex=1.5, type='o', lwd=5)
    }
    dev.off()
    
    ## overlay pasture and forest with points
    fig_path = paste('./Dist_Decay_analysis/figs/', metric, '_', percent,
                     '_', marker, '_DD_onepanel_pts.png', sep='')
    png(file=fig_path, width=480)
    par(mfrow=c(1,1))
    col = c('darkgreen', 'red')
    for (i in seq_along(habitat)) {
      xlim = c(0.05, 80000)
      ylim = range(sapply(rawDDR, function(x) range(x[,2], na.rm=T)))
      if (i == 1) {
        plot(csiml ~ gdist, data=rawDDR[[i]], log='xy', ylim=ylim,
             xlim=xlim, col=col[i], type='n', xlab='',ylab='', axes=F)
        addAxes()
        mtext(side=1, 'Geographic Distance (m)', cex=2, padj=2.5)
        mtext(side=2, 'Community Similarity (Bray-Curtis)', cex=2, padj=-2)
        legend('bottomleft', c('Forest', 'Pasture'), col=col, lwd=5, bty='n', cex=2)
      }
      points(csiml ~ gdist, data=rawDDR[[i]] , col=col[i])
      points(I(1 - exp.var) ~ Dist, data=avgDDR[[i]], col=col[i], pch=19,
             cex=1.5, type='o', lwd=5)
    }
    dev.off()
    
    ## Mantel Test which tests for significant correlation between community turnover and geographic distance
    ## only test within a given habitat type
    #for (i in seq_along(habitat)) {
    #  true = coords$systems == habitat[i]  
    #  gdist = dist(coords[true, c('x','y')])
    #  csiml = 1 - vegdist(comm[true, ])
    #  csiml[gdist > 5000] = NA
    #  gdist[gdist > 5000] = NA
    #  x = mantel(csiml, gdist, na.rm = T, permutations=2000)
    #  print(x)
    #}
    
    ## Carry out analysis ------------------------------------------------------
    ## ADONIS ananlysis
    ## test if distance decay relationships differ between habitat types
    if (metric == 'unifrac') 
      cdist = read.unifrac(unifrac_file, sample_ids=coords$samp)
    if (metric == 'bray')
      cdist = vegdist(comm, method='bray') 
    xy = as.matrix(coords[, c('x', 'y')])
    gdist = dist(xy)
    habitat = ifelse(coords$system == 'Pasture', 0, 1)
    ## Carry out permutation test for difference in avg similarity
    print(paste('metric = ', metric, ', percent = ', percent,
                ', and marker = ', marker, sep=''))
    ## test if there is a difference in beta diversity between
    ## habitat types
    ## TO DO: write model outputs to text file
    print(adonis(cdist ~ habitat))
    ## test if there is a relationship between beta diversity 
    ## spatial proximity, habitat, and their interaction
    print(adonis(cdist ~ xy + habitat + xy * habitat))
  }
}


## potential residual analysis to clean up Pasture result...
i = 2
true = coords$systems == habitat[i] 
crds = coords[true, c('x','y')]
rda_mod = rda(comm[true, ] ~ coords[true, 'treatment'])
tst = residuals(rda_mod)
csiml = 1 - dist(tst)
v = vario(tst, crds, breaks=c(0.01, 0.5, 2, 25, 500, 1500, 15000),
          hmin=0, hmax=15000)$vario
xlim = range(dist(crds), v$Dist)
ylim = range(csiml , v$exp.var)
plot(dist(crds), csiml , main = habitat[i], ylim=ylim, xlim=xlim)
points(v$Dist, 1 - v$exp.var, col='dodgerblue', pch=19, cex=1.5, type='o', lwd=5)


summary(rda_mod)
plot(rda_mod)
rda_mso = mso(rda_mod, crds)
msoplot(rda_mso)

  

## ddR_ babur_test_run_scripts -------------------------------------------------------------

## to analyze distance decay in Forest and Pasture

## load package and necessary functions ------------
library('vegan')
source('./Dist_Decay_analysis/spat_functions.R')

## get abundance file names -----------------------
files = dir('./Dist_Decay_analysis/abu_data') ## get list of files
comm_files = files[grep('FSP', files)] ## pull out file names with FSP

read.unifrac = function(file, sample_ids=NULL) {
  ## reads in unifrac distance matrix and optionally subsets 
  ## arguments
  ## file: file name with path to distance matrix
  ## sample_ids: which samples to keep
  unifrac_dat = read.table(file)
  if (!is.null(sample_ids)) {
    indices = match(sample_ids, rownames(unifrac_dat))
    indices = indices[!is.na(indices)]
    unifrac_dat = unifrac_dat[indices, indices]
  }  
  unifrac_dist = as.dist(unifrac_dat)
  return(list(dist=unifrac_dist, samples=rownames(unifrac_dat)))
}

## taxonomic and unifrac DDR ------------------------------------------------------------------------

calc_ddr = function(x, coords, metric, breaks, hmin=NA, hmax=NA) {
  
  ## computes the distance decay relationship
  ## Returns:
  ## list with first element the raw DDR values, second element the averaged values
  ## Arguments:
  ## x: either a site x sp matrix or a distance matrix
  ## coords: x y spatial coordinates
  ## metric: what similarity index to use
  ## breaks: the spatial breaks to average at
  ## hmin: minimum distance to examine, can be left blank
  ## hmax: maximum distance to examine, can be left blank
  if (class(x) == 'matrix')
    x = vegdist(x, method=metric)
  rawDDR = data.frame(gdist = as.vector(dist(coords)), 
                      csiml = as.vector(1 - x))
  avgDDR = vario(x, coords, breaks=breaks,
                 distance.metric='unifrac', hmin=hmin, hmax=hmax)$vario
  return(list(rawDDR, avgDDR))
}

load_ddr_data = function(file, metric) {
  ## extract marker and percent info from file
  split_name = strsplit(file, '_')             ## split up the file name into its components
  percent = sub('%', '', split_name[[1]][1])   ## identify the % similarity
  marker = sub('.csv', '', split_name[[1]][4]) ## identify the marker
  
  ## import data
  dat = read.csv(file.path('./Dist_Decay_analysis/abu_data/', file))
  
  ## reshape community data for analysis
  comm = t(dat[ , -1])
  colnames(comm) = as.character(dat[ , 1])
  
  coords = read.csv('./Dist_Decay_analysis/FNL_site_dat_coords.csv')
  ## subset coords to only samples we have community data for
  coords = subset(coords, subset=sample %in% rownames(comm))
  
  ## drop secondary forest
  #comm = comm[coords$systems != 'Secondary Forest', ]
  #coords = coords[coords$systems != 'Secondary Forest', ]
  
  ## check that rows of comm so that it matches coords
  #all(coords$sample == rownames(comm))
  
  habitat = unique(coords$systems)
  
  ## compute DDR for plotting----------------------------------------------
  ## compute the average dissimalarity for given spatial breaks 
  if (metric == 'unifrac') {
    unifrac_file = paste('./Dist_Decay_analysis/unifrac_dist/unifrac_dist', percent,
                         marker, 'abu.txt', sep='_') 
  }
  breaks = c(0.01, 0.5, 2, 25, 500, 1500, 15000)
  rawDDR = avgDDR = vector('list', length(habitat))
  names(rawDDR) = names(avgDDR) = habitat
  for (i in seq_along(habitat)) {
    true = coords$systems == habitat[i]                                        ## subset to only habitat of interest
    crds = coords[true, c('x','y')]                                            ## assign coordinates
    if (metric == 'unifrac') {
      unifrac = read.unifrac(unifrac_file, sample_ids= coords$sample[true])
      crds = crds[match(unifrac$samples, coords[true, 'sample']), ]
      cdist = unifrac$dist
    }
    else
      cdist = vegdist(comm[true, ], method = 'bray')
    DDR = calc_ddr(cdist, crds, metric, breaks, hmin=0, hmax=20000)
    rawDDR[[i]] = DDR[[1]]
    avgDDR[[i]] = DDR[[2]]
  }
  out = list(marker=marker, percent=percent, metric=metric, 
             comm=comm, coords=coords[ ,c('x', 'y')], habitat=coords$systems,
             rawDDR=rawDDR, avgDDR=avgDDR)
  return(out)
}

## compare statistical fits-----------------------------------------------------
compare_mods = function(ddr_data) {
  ## compares the r2 values for the linear, exponential, and pwr models
  ## of the DDR within each habitat type
  habitat_types = unique(ddr_data$habitat)
  r2 = matrix(NA, ncol=3, nrow=3)
  colnames(r2) = habitat_types
  rownames(r2) = c('lin','exp','pwr')
  for (i in seq_along(habitat_types)) {
    lin_fit = lm(csiml ~ gdist, data=ddr_data$rawDDR[[i]])
    exp_fit = lm(log(csiml) ~ gdist, data=ddr_data$rawDDR[[i]])
    pwr_fit = lm(log(csiml) ~ log(gdist), data=ddr_data$rawDDR[[i]])
    r2[ , i] = sapply(list(lin_fit, exp_fit, pwr_fit), function(x) summary(x)$r.squared)
  }
  return(r2)
}

bray_mods = lapply(comm_files, function(x) compare_mods(load_ddr_data(x, 'bray')))
unifrac_mods = lapply(comm_files, function(x) compare_mods(load_ddr_data(x, 'unifrac')))

bray_avg = unifrac_avg = matrix(0, 3, 3)
for(i in seq_along(bray_mods)) {
  bray_avg = bray_avg + bray_mods[[i]] / 6
  unifrac_avg = unifrac_avg + unifrac_mods[[i]] / 6
}

bray_avg
unifrac_avg
apply(bray_avg, 1, mean)
apply(unifrac_avg, 1, mean)
# take home is that the exponential tends to always be best choice

## make figures with exponential model fits-------------------------------------

make_figs = function(ddr_data, method, mod=NULL, add_avg=F,
                     path='./Dist_Decay_analysis/figs/') {
  if (method == 'threepanel') {
    fig_path = paste(path, ddr_data$metric, '_', ddr_data$percent, '_',
                     ddr_data$marker, '_DD_threepanel.png', sep='')
    png(file=fig_path, width=480 * 3)
    habitat_types = unique(ddr_data$habitat)
    if (ddr_data$metric == 'bray')
      ylim = c(.001, 1)
    else
      ylim = range(sapply(ddr_data$rawDDR, function(x) range(x$csiml)))
    xlim = range(sapply(ddr_data$rawDDR, function(x) range(x$gdist)))
    #(bottom, left, top, right)
    # c(5, 4, 4, 2) + 0.1
    par(mar=c(5, 5, 4, 1) + 0.1)
    par(mfrow=c(1,3))
    for (j in c('Forest', 'Secondary Forest', 'Pasture')) {
      i = match(j, habitat_types)
      ## assign coordinates
      ## set graphical x and y limits
      #xlim = c(0.001, 80000)
      #ylim = range(sapply(ddr_data$rawDDR, function(x) range(x[,2], na.rm=T)))
      ## plot result on log-log axes and convert dissimarlity to similarity
      if (ddr_data$metric == 'bray')
        metric_type = 'Bray'
      else
        metric_type = 'Unifrac'
      plot(csiml ~ gdist, data=ddr_data$rawDDR[[i]],log='xy', ylim=ylim, xlim=xlim,
           xlab='', ylab='', frame.plot=F, axes=F)
      axis(side=1, cex.axis=1.75, padj=.25, lwd=4)
      if (ddr_data$metric == 'bray')
        axis(side=2, cex.axis=1.75, padj= 0.1, lwd=4, at=c(.001, .01, .1, 1))
      else 
        axis(side=2, cex.axis=1.75, padj= 0.1, lwd=4)
      mtext(side=1, 'Geographic Distance (m)', cex=1.6, padj=2.25)
      mtext(side=2, paste('Community Similarity (', metric_type, ')', sep=''),
            cex=1.5, padj=-1.9)
      mtext(side=3, habitat_types[i], cex=2)
      ## add the averaged result
      if(!is.null(mod)) {
        if (mod == 'lin') {
          fit = lm(csiml ~ gdist, data=ddr_data$rawDDR[[i]])
          pred = predict(fit)
        }
        else if(mod == 'exp') {
          fit = lm(log(csiml) ~ gdist, data=ddr_data$rawDDR[[i]])
          pred = exp(predict(fit))
        }
        else if(mod == 'pwr') {
          fit = lm(log(csiml) ~ log(gdist), data=ddr_data$rawDDR[[i]])
          pred = exp(predict(fit))
        }
        gdist = ddr_data$rawDDR[[i]]$gdist
        lines(gdist[order(gdist)], pred[order(gdist)], lwd=3)
        r2 = round(summary(fit)$r.squared, 2)
        b0 = round(exp(coef(fit)[1]), 3)
        b1 = round(coef(fit)[2], 3)
        if (mod == 'pwr')
          mod_type = 'Power Model'
        if (mod == 'exp')
          mod_type = 'Exponential Model'
        if (mod == 'lin')
          mod_type = 'Linear Model'
        leg_txt = c(mod_type, paste('log(S) = ', b0, ' - ', abs(b1), '*log(D)', sep=''),
                    paste('R^2 =', r2[1]))
        legend('bottomleft', leg_txt, bty='n', cex=2)
      }
      if(add_avg)
        lines(I(1 - exp.var) ~ Dist, data=ddr_data$avgDDR[[i]], cex=1.5, lwd=5,
              type='o')
    }
    dev.off()
  }
}

lapply(comm_files, function(x)
  make_figs(load_ddr_data(x, 'bray'), 'threepanel', 'pwr'))
lapply(comm_files, function(x) 
  make_figs(load_ddr_data(x, 'unifrac'), 'threepanel', 'pwr'))

## ADONIS analysis ------------------------------------------------------------
## test if distance decay relationships differ between habitat types
if (metric == 'unifrac') 
  cdist = read.unifrac(unifrac_file, sample_ids=coords$samp)
if (metric == 'bray')
  cdist = vegdist(comm, method='bray') 
xy = as.matrix(coords[, c('x', 'y')])
gdist = dist(xy)
habitat = ifelse(coords$system == 'Pasture', 0, 1)
## Carry out permutation test for difference in avg similarity
print(paste('metric = ', metric, ', percent = ', percent,
            ', and marker = ', marker, sep=''))
## test if there is a difference in beta diversity between
## habitat types
## TO DO: write model outputs to text file
print(adonis(cdist ~ habitat))
## test if there is a relationship between beta diversity 
## spatial proximity, habitat, and their interaction
print(adonis(cdist ~ xy + habitat + xy * habitat))








