library(crawl)

data("northernFurSeal")
head(northernFurSeal)

colnames(northernFurSeal) <- c("Time", "longitude", "latitude", "Argos_loc_class")

# Argos data with errors 1, 2, 3, A, B
northernFurSeal$Argos_loc_class <- factor(northernFurSeal$Argos_loc_class,
                                          levels=c("3", "2", "1","0","A"))

northernFurSeal$Time <- as.POSIXct(strptime(northernFurSeal$Time, "%Y-%m-%d %H:%M:%S"))

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

fixPar = c(log(250), log(500), log(1500), rep(NA,3), NA)

displayPar(mov.model=~1,
           err.model=list(x=~Argos_loc_class-1),
           data=northernFurSeal,
           fixPar=fixPar)

constr=list(lower=c(rep(log(1500),2), rep(-Inf,2)),
            upper=rep(Inf,4))

ln.prior <- function(theta){
  -abs(theta[4]-3)/0.5
}

# fit

set.seed(123)
fit1 <- crwMLE(mov.model=~1, 
               err.model=list(x=~Argos_loc_class-1),
               data=northernFurSeal, 
               Time.name="Time",
               initial.state=initial,
               fixPar=fixPar, 
               constr=constr, 
               prior=ln.prior,
               control=list(maxit=30, trace=0,REPORT=1),
               initialSANN=list(maxit=200, trace=0, REPORT=1))

fit1

predTime <- seq(lubridate::ceiling_date(min(northernFurSeal$Time), unit = "days"), 
                lubridate::floor_date(max(northernFurSeal$Time), unit = "days"), 86400)


# predict 

predObj <- crwPredict(object.crwFit=fit1, 
                      predTime, 
                      speedEst=TRUE, 
                      flat=TRUE)

crwPredictPlot(predObj, "map")


set.seed(123)

simObj <- crwSimulator(fit1, 
                       predTime, 
                       method="IS", 
                       parIS=100, 
                       df=5, 
                       scale=18/20)


w <- simObj$thetaSampList[[1]][,1]
hist(w*100, main='Importance Sampling Weights', sub='More weights near 1 is desirable')


my.colors <-colorRampPalette(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a'))

iter <- 20
cols <- my.colors(iter)


crwPredictPlot(predObj, 'map')
for(i in 1:iter){
  samp <- crwPostIS(simObj)
  lines(samp$alpha.sim[,'mu.x'], samp$alpha.sim[,'mu.y'],col=cols[i]) 
}

