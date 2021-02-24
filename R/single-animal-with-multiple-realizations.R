#### packages ####

library(crawl)
library(tidyverse)
library(magrittr)
library(sp)
library(rgdal)
library(lubridate)
library(animation)
library(data.table)
library(ggmap)
library(ggblur)
# devtools::install_github("coolbutuseless/ggblur")


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
                      flat = TRUE, 
                      return.type = "flat")


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
  select(Time, contains(".x"), -nu.x, -se.nu.x) %>% 
  gather(key, mu.x, -Time, -se.mu.x)

predObj.x$key %<>%
  str_remove(".x")

predObj.y <- predObj %>%
  select(Time, contains(".y"), -nu.y, -se.nu.y) %>% 
  gather(key, mu.y, -Time, -se.mu.y)

predObj.y$key %<>%
  str_remove(".y")

# combine data 

predObj <- predObj.x %>% 
  inner_join(predObj.y, by = c("Time", "key"))

ggplot(predObj, aes(x = mu.x, y = mu.y)) + 
  geom_path(aes(colour = key), size = 0.2) +  
  theme_map() + 
  theme(legend.position = "none") + 
  scale_colour_manual(values = cols) + 
  xlim(c(-200000, 2350000)) + 
  ylim(c(-6.1e+05, 6.1e+05))


#### Animation ####

# add a column date 

predObj %<>% 
  mutate(date = as_date(Time))


# list all times

all_dates <- predObj %>% 
  .$date %>% 
  unique %>% 
  na.omit %>% 
  sort


# sort data by date time 

predObj %<>% 
  arrange(Time, key)


# try plotting one time point at a time 

ggplot(predObj %>% filter(date == all_dates[[1]]), 
       aes(x = mu.x, y = mu.y)) + 
  geom_point(aes(colour = key)) +
  scale_colour_manual(values = cols) + 
  xlim(c(-200000, 2350000)) + 
  ylim(c(-6.1e+05, 6.1e+05)) + 
  theme_map() + 
  theme(legend.position = "none")


#### simple animation ####

# filter data to 15:00:00 

predObj_sub <- predObj %>% 
  filter(Time %like% "15:00:00")


desc <- c("Plot of Animal with Multiple Realizations")

saveHTML({
  par(mar = c(4.1, 4.1, 0.1, 0.1))
  for (current_time in seq_along(all_dates)) {
    
    time_desc <- all_dates[[current_time]]
    
    df <- predObj_sub %>% 
      filter(date == time_desc)
    
    p <- ggplot(df, aes(x = mu.x, y = mu.y)) + 
      geom_point(aes(colour = key)) +
      scale_colour_manual(values = cols) + 
      xlim(c(-200000, 2350000)) + 
      ylim(c(-6.1e+05, 6.1e+05)) + 
      theme_map() + 
      theme(legend.position = "none") + 
      labs(title = time_desc)
    
    print(p)
    ani.pause()
  }
}, img.name = "multi_proj_plot", imgdir = "multi_proj_dir", htmlfile = "multi_proj.html", 
title = "Multiple Realization Animation", description = desc, interval = 0.2)


#### animatrion with tail #### 

saveHTML({
  par(mar = c(4.1, 4.1, 0.1, 0.1))
  for (current_time in seq_along(all_dates)) {
    
    tail_len <- 2
    
    if(current_time < tail_len){
      days_before <- current_time
    } else {
      days_before <- current_time - tail_len
    }
    
    df <- predObj_sub %>% 
      filter(date %in% (all_dates[days_before:current_time]))
    
    df_today <- df %>% 
      filter(date == all_dates[[current_time]])
    
    p <- ggplot() +
      geom_smooth(data = df, aes(x = mu.x, y = mu.y, group = key, colour = key),
                  method = "loess", formula = y ~ x, se = FALSE, alpha = .6, size = .5) +
      geom_point(data = df_today, aes(x = mu.x, y = mu.y, colour = key), alpha = 0.6, size = 1) +
      labs(title = all_dates[[current_time]]) +
      xlim(c(-200000, 2350000)) + 
      ylim(c(-6.1e+05, 6.1e+05)) + 
      theme_map() + 
      scale_colour_manual(values = cols) + 
      theme(legend.position="none")
    print(p)
    ani.pause()
  }
}, img.name = "multi_proj_plot_tail", imgdir = "multi_proj_dir_tail", htmlfile = "multi_proj_tail.html", 
title = "Multiple Realization Animation with Tails", description = desc, interval = 0.1)


#### Add map ####

predObj_sp <- predObj_sub

coordinates(predObj_sp) <- ~mu.x + mu.y

proj4string(predObj_sp) <- CRS(paste("+proj=aea +lat_1=30 +lat_2=70",
                                     "+lat_0=52 +lon_0=-170 +x_0=0 +y_0=0",
                                     "+ellps=GRS80 +datum=NAD83",
                                     "+units=m +no_defs"))

predObj_sp <- spTransform(predObj_sp, CRS("+proj=longlat"))


# mapbox <- c(predObj_sp@bbox[1] - 5, 
#             predObj_sp@bbox[2] - 1, 
#             predObj_sp@bbox[3] + 5, 
#             predObj_sp@bbox[4] + 1)

# map_frame <- get_map(location = mapbox, source = "stamen", maptype = "toner", zoom = 9)

predObj_sp_df <- as.data.frame(predObj_sp)

# ggmap(map_frame) +
#   geom_point(data = predObj_sp_df, aes(x = mu.x, y = mu.y, colour = key), alpha = 0.6, size = 1)

#### save map frame #### 
save(map_frame, file = "maps/single_animal_with_multiple_realizations.rda")


#### load saved map ####

load(file = "maps/single_animal_with_multiple_realizations.rda")


#### add animation with mapping ####

ani.options(interval = 0.2, ani.width = 1.5*480, ani.height = 1.5*480, verbose = F)

saveHTML({
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  for (current_time in seq_along(all_dates)) {
    
    tail_len <- 2
    
    if(current_time < tail_len){
      days_before <- current_time
    } else {
      days_before <- current_time - tail_len
    }
    
    df <- predObj_sp_df %>% 
      filter(date %in% (all_dates[days_before:current_time]))
    
    df_today <- df %>% 
      filter(date == all_dates[[current_time]])
    
    p <- ggmap(map_frame) +
      geom_smooth(data = df, aes(x = mu.x, y = mu.y, group = key, colour = key),
                  method = "loess", formula = y ~ x, se = FALSE, alpha = .6, size = .5) +
      geom_point(data = df_today, aes(x = mu.x, y = mu.y, colour = key), alpha = 0.6, size = 1) +
      labs(title = all_dates[[current_time]]) +
      theme_map() + 
      scale_colour_manual(values = cols) + 
      theme(legend.position="none")
    print(p)
    ani.pause()
  }
}, img.name = "multi_proj_plot_tail_map", 
imgdir = "multi_proj_dir_tail_map", 
htmlfile = "multi_proj_tail_map.html", 
title = "Multiple Realization Animation with Tails", 
description = desc)



#### create blurry points ####

predObj_sp_blur <- predObj %>% 
  group_by(date) %>% 
  summarise(move_range.x = abs(abs(max(mu.x, na.rm = T)) - abs(min(mu.x, na.rm = T))),
            move_range.y = abs(abs(max(mu.y, na.rm = T)) - abs(min(mu.y, na.rm = T))),
            move_range = move_range.x + move_range.y,
            .groups = "drop") %>% 
  left_join(predObj %>% filter(key == "mu"), by = "date")

ggplot(predObj_sp_blur[1:20,], aes(x = mu.x, y = mu.y)) + 
  geom_point_blur(aes(blur_size = move_range, size = move_range), 
                  blur_steps = 100, col = "#4c848a", alpha = .6) +  
  theme_map()


# create video with blurry points 

saveHTML({
  par(mar = c(4.1, 4.1, 0.1, 0.1))
  for (current_time in seq_along(all_dates)) {
    
    df <- predObj_sp_blur %>% 
      filter(date == all_dates[[current_time]])
    
    p <- ggplot(df, aes(x = mu.x, y = mu.y)) + 
      geom_point_blur(aes(blur_size = move_range, size = move_range), 
                      blur_steps = 100, col = "#4c848a", alpha = .6) +  
      theme_map() + 
      theme(legend.position="none") + 
      labs(title = all_dates[[current_time]]) + 
      xlim(c(-200000, 2350000)) + 
      ylim(c(-6.1e+05, 6.1e+05)) + 
      scale_blur_size(limits = c(4000, 150000)) +
      scale_size(limits = c(4000, 150000))
    print(p)
    ani.pause()
  }
}, img.name = "multi_proj_plot_blur", imgdir = "multi_proj_dir_blur", htmlfile = "multi_proj_blur.html", 
title = "Multiple Realization Animation Blurry", description = desc, interval = 0.1)


#### Blurry points with google map ####

# transform predObj to sp object to change coordinate system 

predObj_sp <- predObj

coordinates(predObj_sp) <- ~mu.x + mu.y

proj4string(predObj_sp) <- CRS(paste("+proj=aea +lat_1=30 +lat_2=70",
                                     "+lat_0=52 +lon_0=-170 +x_0=0 +y_0=0",
                                     "+ellps=GRS80 +datum=NAD83",
                                     "+units=m +no_defs"))

predObj_sp <- spTransform(predObj_sp, CRS("+proj=longlat"))


# tranform predObj_sp to dataframe 

predObj_df <- as.data.frame(predObj_sp)


# create a ggmap with blurry point for a few rows 

# predObj_df_blur <- predObj_df %>% 
#   group_by(date) %>% 
#   summarise(move_range.x = abs(abs(max(mu.x, na.rm = T)) - abs(min(mu.x, na.rm = T))),
#             move_range.y = abs(abs(max(mu.y, na.rm = T)) - abs(min(mu.y, na.rm = T))),
#             move_range = move_range.x + move_range.y,
#             .groups = "drop") %>% 
#   left_join(predObj_df %>% filter(key == "mu"), by = "date")

# ggmap(map_frame) +
#   geom_point_blur(data = predObj_df_blur[1:20,], 
#                   aes(x = mu.x, y = mu.y, blur_size = move_range, size = move_range), 
#                   blur_steps = 100, col = "#4c848a", alpha = .6) + 
#   theme(legend.position = "none")

predObj_df_blur <- predObj_df

ggplot(predObj_df_blur[1:20,] %>% filter(key == "mu"),
       aes(x = mu.x, y = mu.y)) +
  geom_point_blur(aes(blur_size = se.mu.x, size = se.mu.x), blur_steps = 100, col = "#4c848a", alpha = .6) + 
  theme(legend.position = "none")


# predObj_df_blur$move_range <- (predObj_df_blur$se.mu.x)^3

# put blurry point to video 

saveHTML({
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  for (current_time in seq_along(all_dates)) {
    
    lower_limit <- min(predObj_df_blur$se.mu.x)
    higher_limit <- max(predObj_df_blur$se.mu.x)
    
    df <- predObj_df_blur %>% 
      filter(date == all_dates[[current_time]] & key == "mu")
    
    p <- ggmap(map_frame) +
      geom_point_blur(data = df, 
                      aes(x = mu.x, y = mu.y, blur_size = se.mu.x, size = se.mu.x), 
                      blur_steps = 100, col = "#4c848a", alpha = .6) + 
      theme(legend.position = "none") + 
      labs(title = all_dates[[current_time]]) +
      scale_blur_size(limits = c(lower_limit, higher_limit)) + 
      scale_size(limits = c(lower_limit, higher_limit))
    print(p)
    ani.pause()
  }
}, img.name = "multi_proj_plot_blur_ggmap", 
imgdir = "multi_proj_dir_blur_ggmap", 
htmlfile = "multi_proj_blur_ggmap.html", 
title = "Multiple Realization Animation Blurry with Google Map", 
description = desc, interval = 0.1)
 




