#### Packages ####
library(sp)
library(rgdal)
library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(crawl)
library(tidyverse)
library(magrittr)
library(xts)
library(doParallel)
library(argosfilter)
library(animation)


#### load data #### 

data("beardedSeals")

beardedSeals %>% glimpse() 


#### data processing ####

# Adjust Duplicate Times

beardedSeals %>% 
  group_by(deployid, date_time) %>% 
  filter(n() > 1)


# change duplicate times

date_unique <- beardedSeals %>% 
  group_by(deployid) %>%
  do(unique_date = make.time.unique(.$date_time,eps=1)) %>%
  unnest(unique_date) %>%
  mutate(unique_posix = as.POSIXct(.$unique_date,origin='1970-01-01 00:00:00',tz='UTC')) %>%
  arrange(deployid,unique_posix) %>% 
  select(unique_posix)

beardedSeals <- beardedSeals %>% 
  arrange(deployid, date_time) %>%
  bind_cols(date_unique)


# check again, there is no more dups

beardedSeals %>% 
  group_by(deployid, date_time) %>% 
  filter(n() > 1) %>% 
  select(deployid, date_time, unique_posix)


# Remove Obviously Erroneous Locations

beardedSeals <- beardedSeals %>%  
  arrange(deployid, unique_posix)

split_data <- beardedSeals %>% 
  split(.$deployid)

registerDoParallel(cores = 4)

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
  filter(filtered == "not" & !is.na(error_semimajor_axis)) %>%
  arrange(deployid, unique_posix)


# Convert to SpatialPointsDataFrame and Project

beardedSeals %<>% 
  as.data.frame() 

coordinates(beardedSeals) <- ~longitude+latitude

proj4string(beardedSeals) <- CRS("+proj=longlat +datum=WGS84")

beardedSeals %<>% 
  spTransform(., CRS("+init=epsg:3571"))


#### Fit the model #### 

ids <- unique(beardedSeals@data$deployid)      #define seal IDs

registerDoParallel(cores = 4)

model_fits <-
  foreach(i = 1:length(ids)) %dopar% {
    id_data = subset(beardedSeals,deployid == ids[i])
    diag_data = model.matrix( ~ error_semimajor_axis + 
                                error_semiminor_axis + 
                                error_ellipse_orientation,
      model.frame( ~ ., id_data@data, na.action = na.pass)
    )[,-1]
    
    id_data@data = cbind(id_data@data, 
                         crawl::argosDiag2Cov(
                           diag_data[,1], 
                           diag_data[,2], 
                           diag_data[,3]))
    
    init = list(a = c(sp::coordinates(id_data)[1,1],0,
                      sp::coordinates(id_data)[1,2],0),
                P = diag(c(5000 ^ 2, 10 * 3600 ^ 2, 
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
      fixPar = c(1, 1, NA, NA),
      theta = c(log(10), 3),
      initialSANN = list(maxit = 2500),
      control = list(REPORT = 10, trace = 1), time.scale = "days"
    )
    fit
  }

stopImplicitCluster()

names(model_fits) <- ids

print(model_fits)


#### Predict regularly-spaced locations ####

registerDoParallel(cores = 4)

predData <- foreach(i = 1:length(model_fits), .combine = rbind) %dopar% {
  
  model_fits[[i]]$data$unique_posix <- lubridate::with_tz(
    model_fits[[i]]$data$unique_posix,"GMT")
  predTimes <- seq(
    lubridate::ceiling_date(min(model_fits[[i]]$data$unique_posix), "hour"),
    lubridate::floor_date(max(model_fits[[i]]$data$unique_posix), "hour"),
    "1 day")
  tmp = crawl::crwPredict(model_fits[[i]], predTime=predTimes)
}

stopImplicitCluster()

predData$predTimes <- intToPOSIX(predData$TimeNum)


# change projection and change data type  

predData_sp <- predData

coordinates(predData_sp) <- ~mu.x + mu.y

proj4string(predData_sp) <- CRS("+init=epsg:3571")


#### Plot the output ####

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

p1 <- ggplot(predData, aes(x = mu.x, y = mu.y)) + 
  geom_path(aes(colour = deployid)) + 
  labs(x = "easting (meters)", 
       y = "northing (meters)") + 
  theme_map()

p1


#### simulation ####

i = 1

predTimes <- seq(
  lubridate::ceiling_date(min(model_fits[[i]]$data$unique_posix), "hour"),
  lubridate::floor_date(max(model_fits[[i]]$data$unique_posix), "hour"),
  "1 day")

simObj <- crwSimulator(model_fits[[i]],
                       predTimes,
                       method="IS",
                       parIS=100,
                       df=5,
                       scale=18/20)

my.colors <-colorRampPalette(c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a'))

iter <- 20
cols <- my.colors(iter)

predData %<>% na.omit()

crwPredictPlot(predData)

for(i in 1:iter){
  samp <- crwPostIS(simObj)
  lines(samp$alpha.sim[,'mu.x'], samp$alpha.sim[,'mu.y'],col=cols[i])
}


#### Animation of preds ####

# change class to be tidyverse friendly 

class(predData) <- c("bibliometrixDB", "data.frame")


# add a column date 

predData %<>% 
  mutate(date = as_date(date_time))


# list all dates , we will do by date for now 

all_dates <- predData %>% 
  .$date %>%  
  unique %>% 
  na.omit %>% 
  sort


# try plotting one time point at a time 

ggplot(predData %>% 
         filter(date == all_dates[[1]]), 
       aes(x = mu.x, y = mu.y)) + 
  geom_point(aes(colour = deployid)) + 
  labs(x = "easting (meters)", 
       y = "northing (meters)") +
  xlim(c(0, 1000000)) +
  ylim(c(-3200000, -1600000)) + 
  theme_map()


ggplot(predData %>% 
         filter(date == all_dates[[2]]), 
       aes(x = mu.x, y = mu.y)) + 
  geom_point(aes(colour = deployid)) + 
  labs(x = "easting (meters)", 
       y = "northing (meters)") + 
  xlim(c(0, 1000000)) +
  ylim(c(-3200000, -1600000)) + 
  theme_map()

ggplot(predData %>% 
         filter(date == all_dates[[3]]), 
       aes(x = mu.x, y = mu.y)) + 
  geom_point(aes(colour = deployid)) + 
  labs(x = "easting (meters)", 
       y = "northing (meters)") + 
  xlim(c(0, 1000000)) +
  ylim(c(-3200000, -1600000)) + 
  theme_map()

ggplot(predData %>% 
         filter(date == all_dates[[4]]), 
       aes(x = mu.x, y = mu.y)) + 
  geom_point(aes(colour = deployid)) + 
  labs(x = "easting (meters)", 
       y = "northing (meters)") + 
  xlim(c(0, 1000000)) +
  ylim(c(-3200000, -1600000)) + 
  theme_map()


## now try with simple animation 

desc <- c("Plot of Bearded Seal Example")

saveHTML({
  par(mar = c(4.1, 4.1, 0.1, 0.1))
  for (day in seq_along(all_dates)) {
    
    df <- predData %>% 
      filter(date == all_dates[[day]]) %>% 
      group_by(deployid, date) %>% 
      summarise(mu.x = mean(mu.x, na.rm = TRUE), 
                mu.y = mean(mu.y, na.rm = TRUE), 
                .groups = "drop")
    
    p <- ggplot(df, aes(x = mu.x, y = mu.y, colour = deployid)) +
      geom_point(alpha = 0.6, size = 3) +
      labs(title = all_dates[[day]], 
           x = "easting (meters)",
           y = "northing (meters)") +
      xlim(c(0, 1000000)) +
      ylim(c(-3200000, -1600000)) +
      theme_map() + 
      scale_colour_brewer(palette = "Set1") +
      guides(colour = guide_legend(nrow = 1)) + 
      theme(legend.position="bottom", legend.title = element_blank())
    print(p)
    ani.pause()
  }
}, img.name = "bearded_seal_plot", imgdir = "bearded_seal_dir", htmlfile = "bearded_seal.html", 
   title = "Daily Bearded Seals Paths", description = desc, interval = 0.2)


# animatrion with tail 

saveHTML({
  par(mar = c(4.1, 4.1, 0.1, 0.1))
  for (day in seq_along(all_dates)) {
    
    tail_len <- 3
    
    if(day < tail_len){
      days_before <- day
    } else {
      days_before <- day - tail_len
    }
    
    df <- predData %>% 
      filter(date %in% (all_dates[days_before:day])) %>% 
      group_by(deployid, date) %>% 
      summarise(mu.x = mean(mu.x, na.rm = TRUE), 
                mu.y = mean(mu.y, na.rm = TRUE), 
                .groups = "drop")
    
    df_today <- df %>% 
      filter(date == all_dates[[day]])
    
    p <- ggplot() +
      # geom_line(data = df, aes(x = mu.x, y = mu.y, colour = deployid, group = deployid), 
      #           color = "grey", alpha = .6, size = 3) +
      geom_smooth(data = df, aes(x = mu.x, y = mu.y, group = deployid),
                  method = "loess", formula = y ~ x, se = FALSE, alpha = .6, size = 3, color = "grey") +
      geom_point(data = df_today, aes(x = mu.x, y = mu.y, colour = deployid), alpha = 0.6, size = 3.2) +
      labs(title = all_dates[[day]], 
           x = "easting (meters)",
           y = "northing (meters)") +
      xlim(c(0, 1000000)) +
      ylim(c(-3200000, -1600000)) +
      theme_map() + 
      scale_colour_brewer(palette = "Set1") +
      guides(colour = guide_legend(nrow = 1)) + 
      theme(legend.position="bottom", legend.title = element_blank())
    print(p)
    ani.pause()
  }
}, img.name = "bearded_seal_with_tails_plot", 
   imgdir = "bearded_seal_with_tails_dir", 
   htmlfile = "bearded_seal_with_tails.html", 
   title = "Daily Bearded Seals Paths with Tails", 
   interval = 0.1)
 
