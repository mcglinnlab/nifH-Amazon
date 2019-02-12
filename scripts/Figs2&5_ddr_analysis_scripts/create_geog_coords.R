
site_dat = read.csv('./Dist_Decay_analysis/treatment_forest_pasture_secForest_dd.csv')

## compute spacing between coordinates for secondary forest
#SF_dat = site_dat[site_dat$systems == 'Secondary Forest', ]
#SF_dist = sort(dist(SF_dat[!is.na(SF_dat$Easting), c('Easting', 'Northing')]))
#SF_dist 

## add SC forest where no sequence data exists
site_dat = rbind(site_dat, data.frame(sample='SCIA001', treatment='SCI', distance_class=3,
                                      transect='A', systems='Secondary Forest', Distance..m.=9791.317))
## add SC forest where no sequence data exists
site_dat = rbind(site_dat, data.frame(sample='SC93A001', treatment='SC93', distance_class=1,
                                      transect='A', systems='Secondary Forest', Distance..m.=.001))
## add pasture 
site_dat = rbind(site_dat, data.frame(sample='P72IIIA001', treatment='P72III', distance_class=1,
                                      transect='A', systems='Pasture', Distance..m.=.001))


## read in coordinate information
sample_info = read.csv('./Dist_Decay_analysis/sample_info.csv')
## rename 2nd column for merging
names(sample_info)[2] = 'sample'

## add column to sample_info
sample_info$pt_type = 'origin'

## add coordinate information to site_dat
site_dat = merge(site_dat, sample_info, by='sample', all.x=TRUE, all.y=TRUE)




create_coords = function(treatment, transect, pt_type, distance, x, y) {
  coords = matrix(NA, nrow=length(transect), ncol=2)
  trt_levels = unique(treatment)
  
  for (i in 1:length(transect)) {
    starting_index = which(transect == 'A' &
                             treatment == treatment[i] &
                             pt_type=='origin')
    x_start = x[starting_index]
    y_start = y[starting_index]
    if (transect[i] == 'A') {
      y_incr = 0
      if (distance[i] == 0.001)
        x_incr = 0
      else
        x_incr = distance[i]
    }
    if (transect[i] == 'B'){
      x_incr = 0
      y_incr = distance[i]
    }
    if (transect[i] == 'C'){
      y_incr = distance[i]
      x_incr = distance[i]
    }
    coords[i, ] = c(x_start + x_incr, y_new = y_start + y_incr)
  }
  return(coords)
}

coords = create_coords(site_dat$treatment, site_dat$transect, site_dat$Distance..m.,
                       site_dat$Easting, site_dat$Northing)
colnames(coords) = c('x', 'y')

site_dat = cbind(site_dat, coords)

write.csv(site_dat, file='./Dist_Decay_analysis/treatment_forest_pasture,secForest_dd_coords.csv',
          row.names=F)

