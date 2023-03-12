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

library(PrettyCols)
showtext_auto()
options(ofce.background_color = "grey99")
options(ofce.base_family = "sans")
options(ofce.base_size = 9)