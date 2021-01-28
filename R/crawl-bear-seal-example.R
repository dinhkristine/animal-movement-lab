library(sp)
library(rgdal)
library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(crawl)
library(xts)

data("beardedSeals")
beardedSeals

beardedSeals %>% 
  group_by(deployid,date_time) %>% 
  filter(n()>1)


date_unique <- beardedSeals %>% 
  group_by(deployid) %>%
  do(unique_date = xts::make.time.unique(.$date_time,eps=1)) %>%
  tidyr::unnest(unique_date) %>%
  mutate(unique_posix = as.POSIXct(.$unique_date,origin='1970-01-01 00:00:00',tz='UTC')) %>%
  dplyr::arrange(deployid,unique_posix) %>% 
  dplyr::select(unique_posix)

beardedSeals <- beardedSeals %>% arrange(deployid,date_time) %>%
  bind_cols(date_unique)

beardedSeals <- beardedSeals %>%  
  dplyr::arrange(deployid,unique_posix)

library(doParallel)
library(argosfilter)

split_data <- split(beardedSeals,beardedSeals$deployid)

registerDoParallel(cores=2)
beardedSeals$filtered <- foreach(i = 1:length(split_data), .combine = c) %dopar% {
  argosfilter::sdafilter(
    lat=split_data[[i]]$latitude, 
    lon=split_data[[i]]$longitude, 
    dtime=split_data[[i]]$unique_posix,
    lc=split_data[[i]]$quality, 
    ang=-1,
    vmax=5)
}
stopImplicitCluster()

beardedSeals <- beardedSeals %>% 
  dplyr::filter(., filtered=="not" & !is.na(error_semimajor_axis)) %>%
  arrange(.,deployid,unique_posix)


beardedSeals <- as.data.frame(beardedSeals)
coordinates(beardedSeals) = ~longitude+latitude
proj4string(beardedSeals) = CRS("+proj=longlat +datum=WGS84")

beardedSeals <- spTransform(beardedSeals, CRS("+init=epsg:3571"))


# fit the model 

ids = unique(beardedSeals@data$deployid)      #define seal IDs

registerDoParallel(cores=2)
model_fits <-
  foreach(i = 1:length(ids)) %dopar% {
    id_data = subset(beardedSeals,deployid == ids[i])
    diag_data = model.matrix(
      ~ error_semimajor_axis + error_semiminor_axis + error_ellipse_orientation,
      model.frame( ~ ., id_data@data, na.action = na.pass)
    )[,-1]
    
    id_data@data = cbind(id_data@data, 
                         crawl::argosDiag2Cov(
                           diag_data[,1], 
                           diag_data[,2], 
                           diag_data[,3]))
    
    init = list(a = c(sp::coordinates(id_data)[1,1],0,
                      sp::coordinates(id_data)[1,2],0),
                P = diag(c(5000 ^ 2,10 * 3600 ^ 2, 
                           5000 ^ 2, 10 * 3600 ^ 2)))
    
    fit <- crawl::crwMLE(
      mov.model =  ~ 1,
      err.model = list(
        x =  ~ ln.sd.x - 1, 
        y =  ~ ln.sd.y - 1, 
        rho =  ~ error.corr
      ),
      data = id_data,
      Time.name = "unique_posix",
      initial.state = init,
      fixPar = c(1,1,NA,NA),
      theta = c(log(10), 3),
      initialSANN = list(maxit = 2500),
      control = list(REPORT = 10, trace = 1)
    )
    fit
  }
stopImplicitCluster()

names(model_fits) <- ids

print(model_fits)


# predict 

registerDoParallel(cores=2)
predData <- foreach(i = 1:length(model_fits), .combine = rbind) %dopar% {
  
  model_fits[[i]]$data$unique_posix <- lubridate::with_tz(
    model_fits[[i]]$data$unique_posix,"GMT")
  predTimes <- seq(
    lubridate::ceiling_date(min(model_fits[[i]]$data$unique_posix),"hour"),
    lubridate::floor_date(max(model_fits[[i]]$data$unique_posix),"hour"),
    "1 hour")
  tmp = crawl::crwPredict(model_fits[[i]], predTime=predTimes)
}
stopImplicitCluster()

predData$predTimes <- intToPOSIX(predData$TimeNum)

predData_sp <- predData
coordinates(predData_sp) <- ~mu.x+mu.y
proj4string(predData_sp) <- CRS("+init=epsg:3571")


# plot preds 

theme_map = function(base_size=9, base_family="")
{
  require(grid)
  theme_bw(base_size=base_size, base_family=base_family) %+replace%
    theme(axis.title.x=element_text(vjust=0),
          axis.title.y=element_text(angle=90, vjust=1.25),
          axis.text.y=element_text(angle=90),
          axis.ticks=element_line(colour="black", size=0.25),
          legend.background=element_rect(fill=NA, colour=NA),
          legend.direction="vertical",
          legend.key=element_rect(fill=NA, colour="white"),
          legend.text=element_text(),
          legend.title=element_text(face="bold", hjust=0),
          panel.border=element_rect(fill=NA, colour="black"),
          panel.grid.major=element_line(colour="grey92", size=0.3, linetype=1),
          panel.grid.minor=element_blank(),
          plot.title=element_text(vjust=1),
          strip.background=element_rect(fill="grey90", colour="black", size=0.3),
          strip.text=element_text()
    )
}

p1 <- ggplot(data=predData,aes(x=mu.x,y=mu.y)) + 
  geom_path(aes(colour=deployid)) + xlab("easting (meters)") +
  ylab("northing (meters)") + theme_map()
p1



