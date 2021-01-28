library(crawl)

data("northernFurSeal")
head(northernFurSeal)

colnames(northernFurSeal) <- c("time", "longitude", "latitude", "Argos_loc_class")

# Argos data with errors 1, 2, 3, A, B
northernFurSeal$Argos_loc_class <- factor(northernFurSeal$Argos_loc_class,
                                          levels=c("3", "2", "1","0","A"))

library(sp)
library(rgdal)

coordinates(northernFurSeal) = ~longitude+latitude

proj4string(northernFurSeal) <- CRS("+proj=longlat")

northernFurSeal <- spTransform(northernFurSeal, 
                               CRS(paste("+proj=aea +lat_1=30 +lat_2=70",
                                         "+lat_0=52 +lon_0=-170 +x_0=0 +y_0=0",
                                         "+ellps=GRS80 +datum=NAD83",
                                         "+units=m +no_defs"))
)

initial = list(a=c(coordinates(northernFurSeal)[1,1],0,
                   coordinates(northernFurSeal)[1,2],0),
               P=diag(c(10000^2,54000^2,10000^2,5400^2)))
