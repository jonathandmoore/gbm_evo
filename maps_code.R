setwd("X:/Daniel-zilberman/Projects/Jay-1001_methylomes")
data1=read.table(file="Maps/LFGBM_alleles_Del_UM_TEM_and_gbM.txt", sep="\t", header=TRUE)
data2= read.table(file="Maps/LFGBM_Southern_Sweden.txt", sep="\t", header=TRUE)
data3= read.table(file="Maps/FLC_UM_and_gbM.txt", sep="\t", header=TRUE)
data4= read.table(file="Maps/FLC_UM_and_gbM_Northern_Sweden.txt", sep="\t", header=TRUE)
data5= read.table(file="Maps/FLC_UM_and_gbM_48_to_52_degree_latitude.txt", sep="\t", header=TRUE)

library(ggmap)
library(mapproj)
register_google(key="AIzaSyBYZ2MkShIY-zsROZADNthYv0KZVvYpNyE")

our_maptype="satellite"
our_maptype="terrain-background"

map <- get_map(location = 'Europe', zoom = 4 , maptype=our_maptype)
ggmap(map) + geom_point(data=data1, aes(x=Longitude, y=Latitude, colour=Allele))
ggmap(map) + geom_point(data=data3, aes(x=Longitude, y=Latitude, colour=FLC_meth_status))
ggmap(map) + geom_point(data=data5, aes(x=Longitude, y=Latitude, colour=FLC_meth_status))

map <- get_map(location = 'Sweden', zoom = 5 , maptype=our_maptype)
ggmap(map) + geom_point(data=data2, aes(x=Longitude, y=Latitude, colour=Allele))
ggmap(map) + geom_point(data=data4, aes(x=Longitude, y=Latitude, colour=FLC_meth_status))




library(maps)
europe <- map_data("world")

eur_xlim=c(-10.5,41)
eur_ylim=c(37,67)

ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=eur_xlim, ylim= eur_ylim) + geom_point(data=data1, shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=Allele)) +scale_fill_manual(breaks=c("gbM", "UM", "Del", "TEm"), values=c("orange","white","grey","cyan"))

ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=c(10,20), ylim=c(55,60)) + geom_point(data=data2, shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=Allele)) +scale_fill_manual(breaks=c("gbM", "UM", "Del", "TEm"), values=c("orange","white","grey","cyan"))

ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=eur_xlim, ylim=eur_ylim) + geom_point(data=data3, shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=FLC_meth_status)) +scale_fill_manual(breaks=c("gbm", "UM"), values=c("orange","white"))

ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=c(10,20), ylim=c(60,65)) + geom_point(data=data4, shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=FLC_meth_status)) +scale_fill_manual(breaks=c("gbm", "UM"), values=c("orange","white"))

ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=eur_xlim, ylim=c(47,53)) + geom_point(data=data5, shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=FLC_meth_status)) +scale_fill_manual(breaks=c("gbm", "UM"), values=c("orange","white"))

#Version with UM on top (it was getting overplotted by gbM):
ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=eur_xlim, ylim=c(47,53)) + geom_point(data=data5[data5$FLC_meth_status=="gbm",], shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=FLC_meth_status))+ geom_point(data=data5[data5$FLC_meth_status=="UM",], shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=FLC_meth_status)) +scale_fill_manual(breaks=c("gbm", "UM"), values=c("orange","white"))

#LFGBM alleles Del UM TEM and gbM= full map
#LFGBM Southern Sweden= 55-60
#FLC UM and gbM= full map
#FLC UM and gbM Northern Sweden= 60-65
#FLC UM and gbM 48 to 52 degree latitude = 47-53


#I am coding gbM with orange color and UM with white for FLC. Both have black outlines. For LFGBM, gbM is orange, UM is white, Deletions are grey, and TEm is Turquoise, again all have black outline colors.

# Recode colours - white and grey are difficult to read, use black and red instead
ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=eur_xlim, ylim= eur_ylim) + geom_point(data=data1, shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=Allele)) +scale_fill_manual(breaks=c("gbM", "UM", "Del", "TEm"), values=c("orange","black","red","cyan"))

ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=c(10,20), ylim=c(55,60)) + geom_point(data=data2, shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=Allele)) +scale_fill_manual(breaks=c("gbM", "UM", "Del", "TEm"), values=c("orange","black","red","cyan"))

ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=eur_xlim, ylim=eur_ylim) + geom_point(data=data3, shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=FLC_meth_status)) +scale_fill_manual(breaks=c("gbm", "UM"), values=c("orange","black"))

ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=c(10,20), ylim=c(60,65)) + geom_point(data=data4, shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=FLC_meth_status)) +scale_fill_manual(breaks=c("gbm", "UM"), values=c("orange","black"))

ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=eur_xlim, ylim=c(47,53)) + geom_point(data=data5, shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=FLC_meth_status)) +scale_fill_manual(breaks=c("gbm", "UM"), values=c("orange","black"))

#Version with UM on top (it was getting overplotted by gbM):
ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=eur_xlim, ylim=c(47,53)) + geom_point(data=data5[data5$FLC_meth_status=="gbm",], shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=FLC_meth_status))+ geom_point(data=data5[data5$FLC_meth_status=="UM",], shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=FLC_meth_status)) +scale_fill_manual(breaks=c("gbm", "UM"), values=c("orange","black"))



#Make points transparent

# Recode colours - white and grey are difficult to read, use black and red instead
ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group, alpha=0.5), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=eur_xlim, ylim= eur_ylim) + geom_point(data=data1, shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=Allele), alpha=0.5) +scale_fill_manual(breaks=c("gbM", "UM", "Del", "TEm"), values=c("orange","black","red","cyan"))

ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=c(10,20), ylim=c(55,60)) + geom_point(data=data2, shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=Allele), alpha=0.5) +scale_fill_manual(breaks=c("gbM", "UM", "Del", "TEm"), values=c("orange","black","red","cyan"))

ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=eur_xlim, ylim=eur_ylim) + geom_point(data=data3, shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=FLC_meth_status), alpha=0.5) +scale_fill_manual(breaks=c("gbm", "UM"), values=c("orange","black"))

ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=c(10,20), ylim=c(60,65)) + geom_point(data=data4, shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=FLC_meth_status), alpha=0.5) +scale_fill_manual(breaks=c("gbm", "UM"), values=c("orange","black"))

ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=eur_xlim, ylim=c(47,53)) + geom_point(data=data5, shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=FLC_meth_status), alpha=0.5) +scale_fill_manual(breaks=c("gbm", "UM"), values=c("orange","black"))

#Version with UM on top (it was getting overplotted by gbM):
ggplot() + geom_polygon(data=europe, aes(x=long, y=lat, group=group), fill=NA, colour="black") + theme_minimal() + coord_fixed(xlim=eur_xlim, ylim=c(47,53)) + geom_point(data=data5[data5$FLC_meth_status=="gbm",], shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=FLC_meth_status), alpha=0.5)+ geom_point(data=data5[data5$FLC_meth_status=="UM",], shape=21, stroke=1, aes(x=Longitude, y=Latitude, fill=FLC_meth_status), alpha=0.5) +scale_fill_manual(breaks=c("gbm", "UM"), values=c("orange","black"))

