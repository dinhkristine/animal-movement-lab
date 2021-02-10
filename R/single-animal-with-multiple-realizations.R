#### packages ####

library(crawl)
library(tidyverse)
library(magrittr)
library(sp)
library(rgdal)
library(lubridate)


#### Load data #### 

data("northernFurSeal")

northernFurSeal %>% head

colnames(northernFurSeal) <- c("Time", "longitude", "latitude", "Argos_loc_class")

northernFurSeal$Argos_loc_class <- factor(northernFurSeal$Argos_loc_class,
                                          levels=c("3", "2", "1","0","A"))

# project the data 

coordinates(northernFurSeal) = ~longitude+latitude


# get projection info from the data 

proj4string(northernFurSeal) <- CRS("+proj=longlat")


# spatial transformation 


northernFurSeal <- spTransform(northernFurSeal, 
                               CRS(paste("+proj=aea +lat_1=30 +lat_2=70",
                                         "+lat_0=52 +lon_0=-170 +x_0=0 +y_0=0",
                                         "+ellps=GRS80 +datum=NAD83",
                                         "+units=m +no_defs")))


#### set initial parameters and priors for the model ####

initial <- list(a = c(coordinates(northernFurSeal)[1,1],0,
                      coordinates(northernFurSeal)[1,2],0),
                P = diag(c(10000^2,54000^2,10000^2,5400^2)))


# fix values in the model to 0 or 1 

fixPar <- c(log(250), log(500), log(1500), rep(NA,3), NA)

displayPar(mov.model=~1,
           err.model=list(x=~Argos_loc_class-1),
           data=northernFurSeal,
           fixPar=fixPar)

# parameters for upper and lower limits 

constr <- list(lower=c(rep(log(1500),2), rep(-Inf,2)),
               upper=rep(Inf,4))

# set a prior 

ln.prior <- function(theta){-abs(theta[4]-3)/0.5}


#### Fit the model ####

set.seed(123)

fit1 <- crwMLE(mov.model = ~ 1, 
               err.model = list(x = ~ Argos_loc_class - 1),
               data = northernFurSeal, 
               Time.name = "Time",
               initial.state = initial,
               fixPar = fixPar, 
               constr = constr, 
               prior = ln.prior,
               control = list(maxit = 30, trace = 0,REPORT = 1),
               initialSANN = list(maxit = 200, trace = 0, REPORT = 1))


# view model results 

fit1


#### Predict regularly-spaced locations ####

# define min and max times in the data

predTime <- seq(ceiling_date(min(northernFurSeal$Time), unit = "hours"), 
                floor_date(max(northernFurSeal$Time), unit = "hours"), 3600)


# predict using the time defined 

predObj <- crwPredict(object.crwFit = fit1, 
                      predTime, 
                      speedEst = TRUE, 
                      flat = TRUE)


# theme for the plot 

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


# plot the predictions

ggplot(predObj, aes(x = mu.x, y = mu.y)) + 
  geom_path() + 
  theme_map()


#### Simulation ####

set.seed(123)

simObj <- crwSimulator(fit1, 
                       predTime, 
                       method = "IS", 
                       parIS = 100, 
                       df = 5, 
                       scale = 18 / 20)


#### Examine the simulation ####

# important sampling weight distribution 

w <- simObj$thetaSampList[[1]][,1]

hist(w*100, 
     main = 'Importance Sampling Weights', 
     sub = 'More weights near 1 is desirable')


# number of independent samples 

round(100/(1+(sd(w)/mean(w))^2))


#### Plot simulation ####

# define colors 

my.colors <- colorRampPalette(c('#a6cee3','#1f78b4','#b2df8a',
                                '#33a02c','#fb9a99','#e31a1c',
                                '#fdbf6f','#ff7f00','#cab2d6',
                                '#6a3d9a'))


# set the number of tracks want to sample and define colors 

iter <- 5

cols <- my.colors(iter + 2)


# change pred class to dataframe 

class(predObj) <- c("bibliometrixDB", "data.frame")


# set loops to add simulation paths 

for(i in 1:iter){
  samp <- crwPostIS(simObj) 
  
  predObj %<>% 
    mutate(samp.x = samp$alpha.sim[,'mu.x'], 
           samp.y = samp$alpha.sim[,'mu.y'])
  
  colnames(predObj)[colnames(predObj) == "samp.x"] <- paste0("samp.x", i)
  colnames(predObj)[colnames(predObj) == "samp.y"] <- paste0("samp.y", i)
}


# gather data 

predObj.x <- predObj %>%
  select(Time, contains(".x"), -nu.x) %>% 
  gather(key, mu.x, -Time)

predObj.x$key %<>%
  str_remove(".x")

predObj.y <- predObj %>%
  select(Time, contains(".y"), -nu.x) %>% 
  gather(key, mu.y, -Time)

predObj.y$key %<>%
  str_remove(".y")

# combine data 

predObj <- predObj.x %>% 
  inner_join(predObj.y, by = c("Time", "key"))

ggplot(predObj, aes(x = mu.x, y = mu.y)) + 
  geom_path(aes(colour = key), size = 0.2) +  
  theme_map() + 
  theme(legend.position = "none") + 
  scale_colour_manual(values = cols)









