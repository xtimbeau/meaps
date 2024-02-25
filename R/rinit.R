library(knitr)
opts_chunk$set(
  echo = FALSE,
  warning = FALSE,
  message = FALSE,
  fig.pos="htb", 
  out.extra="",
  dev="ragg_png",
  out.width="100%",
  fig.showtext=TRUE,
  cache=FALSE)

library(tidyverse)
library(ofce)
library(showtext)
library(markdown)
library(gt)
library(qs)
library(sf)
library(ggspatial)
library(stars)
library(glue)
library(gtsummary)
library(scales)
library(broom)
library(broom.mixed)
library(arrow)
library(vroom)
library(mapdeck)
options(
  ofce.base_family = "Open Sans",
  ofce.background_color = "transparent")

tkn <- Sys.getenv("mapbox_token")
mapdeck::set_token(tkn)

style <- "mapbox://styles/xtimbeau/ckyx5exex000r15n0rljbh8od"

showtext_auto()
options(ofce.background_color = "grey99")
options(ofce.base_family = "sans")
options(ofce.base_size = 9)

source("secrets/azure.R")
