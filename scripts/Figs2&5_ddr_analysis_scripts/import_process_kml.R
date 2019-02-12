library(rgdal) 
library(maps)

sample_kml = readOGR('./Dist_Decay_analysis/ARMOallwaypoints.kml', 'ARMOallwaypoints') 

## check coordinates are reasonable
map('world', 'brazil')
points(sample_kml)

## project decimal degrees to UTM zone 20 L
sample_utm = spTransform(sample_kml, 
                         CRSobj=CRS("+proj=utm +zone=20 +units=m +datum=WGS84"))

## extract dataframe from sample_kml
map_data = sample_utm@data
map_data = data.frame(map_data, coordinates(sample_utm))
names(map_data) = c('Name', 'Date_of_sample', 'Easting', 'Northing', 'Elevation')

## pull out coordinate info for samples of interest
sample_lookup = read.csv('./Dist_Decay_analysis/sample_names_JR.csv')

sample_info = sample_lookup[as.character(sample_lookup$Sample.names.of.Google.map) != '', ]
## take only the SW coordinate at 001
sample_info = sample_info[grep('001', as.character(sample_info$Sample)), ]

## add a row for Secondary Forest SCA001 which does not have sequence data
sample_info = rbind(sample_info, data.frame(Sample= 'SCIA001', 
                                            Sample.names.of.Google.map='SCA001'))

## rename first column so that it matches the name used in map_data
names(sample_info)[2] = 'Name'

## now merge sample_info and map_data
sample_info = merge(sample_info, map_data, by='Name')

## correct name difference
sample_info$Sample = as.character(sample_info$Sample)
sample_info$Sample[sample_info$Name == 'ISC93A001V2'] =  'SC93A001'

## export sample_info as csv file
write.csv(sample_info, file='./Dist_Decay_analysis/sample_info.csv', row.names=F)

## this code only works when using decimal degrees not UTM
## exampine on map
plot(1:10, 1:10, type='n', xlim=c(-63,-62.7), ylim=c(-10.3, -10.1))
points(sample_kml)
points(sample_info$Longitude, sample_info$Latitude, col='red', pch=19)

## export lat/long projrection as kml file

coords = as.matrix(sample_info[, c('Longitude', 'Latitude')])
sample_info_sp = SpatialPointsDataFrame(coords, 
                                        data=sample_info)  

kmlPoints(sample_info_sp, kmlfile='./Dist_Decay_analysis/ddr_sample_points.kml')

