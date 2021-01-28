#### Packages ####

library(anipaths)
library(tidyverse)
library(magrittr)


#### Load data ####

owl <- read_csv("data/Short-eared Owl, North America.csv")

owl %<>% 
  as.data.frame()

#### generate movement ####

animate_paths(paths = owl,
              delta.t = "week",
              coord = c("location-long", "location-lat"),
              Time.name = "timestamp",
              ID.name = "tag-local-identifier")

