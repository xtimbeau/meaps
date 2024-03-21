library(tidyverse)
library(glue)
library(stars)
library(raster)
library(tmap)
library(r3035)
library(conflicted)
library(mapdeck)

conflicted::conflict_prefer_all("dplyr", quiet = TRUE)
GD <- Sys.getenv("GOOGLE_DRIVE")
rda <- glue("{GD}/DVFdata/rda")
tkn <- Sys.getenv("mapbox_token")
mapdeck::set_token(tkn)
style <- "mapbox://styles/xtimbeau/ckyx5exex000r15n0rljbh8od"
trim <- function(x, xm, xp) ifelse( x<= xm, xm, ifelse(x>= xp, xp, x))

# zones et données générales
c200 <- qs::qread(glue("{GD}/DVFdata/rda/c200.rda"))
iris15 <- qs::qread(glue("{GD}/DVFdata/rda/iris15.rda"))
uu758  <- iris15 %>% filter(UU2010=="00758") %>% st_union # Lyon
uu759  <- iris15 %>% filter(UU2010=="00759") %>% st_union # Marseille
uu851  <- iris15 %>% filter(UU2010=="00851") %>% st_union # Paris
c_uu851 <- uu851 %>% st_centroid
c_uu758 <- uu758 %>% st_centroid
c_uu759 <- uu759 %>% st_centroid
bb758 <- st_bbox((uu851[[1]]-c_uu851[[1]]+c_uu758[[1]]), crs=3035)
bb759 <- st_bbox((uu851[[1]]-c_uu851[[1]]+c_uu759[[1]]), crs=3035)
bb851 <- st_bbox((uu851[[1]]), crs=3035)

# Paris ------------------
access_paris <- qs::qread(glue("{GD}/DVFdata/AccessVilles/Access200_Paris.rda"))
iso_paris <- accessibility::iso2time(as(access_paris$emplois, "Raster") , seuils=c(100000)) |> 
  r2dt() |> 
  as_tibble() |> 
  select(idINS = idINS200, to100k) |>
  left_join(c200 |> as_tibble() |> select(idINS=idINS200, Ind, geometry), by="idINS") |> 
  mutate(x = idINS2point(idINS)[,1], y = idINS2point(idINS)[,2]) |>
  drop_na(Ind) 

# Lyon ----------------
access_lyon <- qs::qread(glue("{rda}/iso_transit_50_r5_Lyon.rda"))$EMP09 |> raster::aggregate(fact=4)
iso_lyon <- accessibility::iso2time(as(access_lyon, "Raster") , seuils=c(100000)) |> 
  r2dt() |> 
  as_tibble() |> 
  select(idINS = idINS200, to100k) |>
  left_join(c200 |> as_tibble() |> select(idINS=idINS200, Ind, geometry), by="idINS") |> 
  mutate(x = idINS2point(idINS)[,1], y = idINS2point(idINS)[,2]) |>
  drop_na(Ind) 

# Marseille ------------
access_marseille <- qs::qread(glue("{rda}/iso_transit_50_r5_Marseille.rda"))$EMP09 |> raster::aggregate(fact=4)
iso_marseille <- accessibility::iso2time(as(access_marseille, "Raster") , seuils=c(100000)) |>
  r2dt() |> 
  as_tibble() |> 
  select(idINS = idINS200, to100k) |>
  left_join(c200 |> as_tibble() |> select(idINS=idINS200, Ind, geometry), by="idINS") |> 
  mutate(x = idINS2point(idINS)[,1], y = idINS2point(idINS)[,2]) |>
  drop_na(Ind) 

# fonds de carte --------------

bb_paris <- iris15 %>% filter(UU2010=="00851") %>% st_union() %>% st_buffer(50000) %>% st_transform(4326) 
st_crs(bb_paris) <- st_crs("+proj=longlat +ellps=WGS84") 

paris.mb <- mapboxapi::get_static_tiles(
  location = bb_paris,
  zoom=8, 
  style_id = "ckjka0noe1eg819qrhuu1vigs", 
  username="xtimbeau",
  access_token=Sys.getenv("mapbox_token")) 

paris.mb <- stars::st_as_stars(paris.mb)[,,,1:3] |> st_transform(3035) |> st_rgb()

bb_lyon <- iris15 %>% filter(UU2010=="00758") %>% st_union() %>% st_buffer(50000) %>% st_transform(4326) 
st_crs(bb_lyon) <- st_crs("+proj=longlat +ellps=WGS84") 

lyon.mb <- mapboxapi::get_static_tiles(
  location = bb_lyon,
  zoom=8, 
  style_id = "ckjka0noe1eg819qrhuu1vigs", 
  username="xtimbeau",
  access_token=Sys.getenv("mapbox_token")) 

lyon.mb <- stars::st_as_stars(lyon.mb)[,,,1:3] |> st_transform(3035) |> st_rgb()

bb_marseille <- iris15 %>% filter(UU2010=="00759") %>% st_union() %>% st_buffer(50000) %>% st_transform(4326) 
st_crs(bb_marseille) <- st_crs("+proj=longlat +ellps=WGS84") 

marseille.mb <- mapboxapi::get_static_tiles(
  location = bb_marseille,
  zoom=8, 
  style_id = "ckjka0noe1eg819qrhuu1vigs", 
  username="xtimbeau",
  access_token=Sys.getenv("mapbox_token")) 

marseille.mb <- stars::st_as_stars(marseille.mb)[,,,1:3] |> st_transform(3035) |> st_rgb()

# ggplot -------------
stars_paris <- iso_paris |> st_drop_geometry() |> select(x, y, to100k) |> dt2r(resolution = 200) |> st_as_stars()
p <- ggplot()+geom_stars(data=paris.mb)+
  geom_stars(data=stars_paris)+
  scale_fill_distiller(name = "temps d'accès en TC (min)\n100k emplois", 
                       palette="YlGnBu",
                       na.value="transparent", 
                       limits = c(0,90),
                       oob = scales::oob_squish) +
  labs(title="Agglomération de Paris (uu00851)")+
  coord_sf(xlim=c(bb851$xmin, bb851$xmax), ylim=c(bb851$ymin, bb851$ymax))+
  ggspatial::annotation_scale() +
  theme_void()+
  theme(plot.title = element_text(size=9))

stars_lyon <- iso_lyon |> 
  st_drop_geometry() |>
  select(x, y, to100k) |> 
  dt2r(resolution = 200) |> 
  st_as_stars()
l <- ggplot()+geom_stars(data=lyon.mb)+
  geom_stars(data=stars_lyon)+
  scale_fill_distiller(name = "temps d'accès en TC (min)\n100k emplois", 
                       palette="YlGnBu",
                       na.value="transparent", 
                       limits = c(0,90),
                       oob = scales::oob_squish) +
  labs(title="Métropole de Lyon (uu00758)")+
  coord_sf(xlim=c(bb758$xmin, bb758$xmax), ylim=c(bb758$ymin, bb758$ymax))+
  ggspatial::annotation_scale() +
  theme_void()+
  theme(plot.title = element_text(size=9))

stars_marseille <- iso_marseille |> 
  st_drop_geometry() |>
  select(x, y, to100k) |> 
  dt2r(resolution = 200) |> 
  st_as_stars()
m <- ggplot()+geom_stars(data=marseille.mb)+
  geom_stars(data=stars_marseille)+
  scale_fill_distiller(name = "temps d'accès en TC (min)\n100k emplois", 
                       palette="YlGnBu",
                       na.value="transparent", 
                       limits = c(0,90),
                       oob = scales::oob_squish) +
  labs(title="Aix-Marseille-Provence (uu00759)")+
  coord_sf(xlim=c(bb759$xmin, bb759$xmax), ylim=c(bb759$ymin, bb759$ymax))+
  ggspatial::annotation_scale() +
  theme_void()+
  theme(plot.title = element_text(size=9))

library(patchwork)
access_plm <- p+l+m + 
  plot_layout(guides = 'collect') &
  theme(legend.position = "bottom", 
        legend.key.height = unit(6, "pt"),
        plot.margin = margin(b = 8, t = 8, l = 2, r = 2),
        legend.title = element_text(size=9),
        legend.text = element_text(size=6))
source("secrets/azure.R")
bd_write(access_plm)
access_plm <- bd_read("access_plm")
graph2png(access_plm, "access_plm", rep = "_presentations", ratio = 16/10 )
# mapdeck -------------
library(colourvalues)
pal <- grDevices::colorRamp(c("blue","green", "yellow"), bias =1)( (0:90)/90)
color_txxk <- function(x, min = 0, max = 90, palette = pal) {
  color_values(c(min, max, trim(x, min, max)), palette=palette) |> utils::tail(-2)
}
legend <- color_values(0:90, palette = pal, summary=TRUE, n_summaries=4)
le <- legend_element(as.integer(legend$summary_values), legend$summary_colours, "fill", "gradient", "minutes") |> 
  mapdeck_legend()

bind_rows(iso_paris, iso_lyon, iso_marseille) |> 
  st_as_sf() |>
  st_transform(4326) |>
  mutate(to100k = color_txxk(to100k)) |> 
  mapdeck(style = style,
        height = "100vh",
        width = "100%") |> 
  add_polygon(fill_colour = "to100k", 
              elevation = "Ind", 
              elevation_scale = 1.5, 
              legend = le)

bb_paris <- iris15 %>% filter(UU2010=="00851") %>% st_union() %>% st_buffer(50000) %>% st_transform(4326) 
st_crs(bb_paris) <- st_crs("+proj=longlat +ellps=WGS84") 

paris.mb <- mapboxapi::get_static_tiles(
  location = bb_paris,
  zoom=8, 
  style_id = "ckjka0noe1eg819qrhuu1vigs", 
  username="xtimbeau",
  access_token=Sys.getenv("mapbox_token")) 

paris.mb <- stars::st_as_stars(paris.mb)[,,,1:3] |> st_transform(3035) |> st_rgb()

p <- tm_shape(paris.mb, bbox = bb851)+tm_rgb()+
  tm_shape(iso_paris)+tm_raster(style="cont", palette="cividis", breaks=c(20,40,60,80,100))

ggplot() + geom_stars(data=paris.mb) + geom_stars(data=st_as_stars(iso_paris$to50k))

bb_lyon <- iris15 %>% filter(UU2010=="00758") %>% st_buffer(75000) %>% st_union() %>% st_transform(4326)
st_crs(bb_lyon) <- st_crs("+proj=longlat +ellps=WGS84")
lyon.mb <- cc_location(loc=bb_lyon, zoom = 9,
                       base_url = "https://api.mapbox.com/styles/v1/{username}/{style_id}/tiles/512/{zoom}/{x}/{y}")
maxs <- cellStats(lyon.mb, max)
lyon.mb <- projectRaster(from=lyon.mb, crs=st_crs(3035)$proj4string) # la projection fait un truc bizarre sur les entiers
lyon.mb <- lyon.mb/cellStats(lyon.mb, max)*maxs %>% as.integer # on remet tout comme avant mais en 3035
iso_lyon <- lload_DVF("iso_transit_50_r5_Lyon")
lyon <- iso2time(iso_lyon$EMP09, seuils=c(50000, 250000))
l <- tm_shape(lyon.mb, bbox = bb758)+tm_rgb()+tm_shape(lyon)+tm_raster(style="cont", palette=terrain.colors(20, rev=FALSE), breaks=c(20,40,60,80,100))

bb_marseille <- iris15 %>% filter(UU2010=="00759") %>% st_buffer(50000) %>% st_union() %>% st_transform(4326)
st_crs(bb_marseille) <- st_crs("+proj=longlat +ellps=WGS84")
marseille.mb <- cc_location(loc=bb_marseille, zoom = 9,
                            base_url = "https://api.mapbox.com/styles/v1/{username}/{style_id}/tiles/512/{zoom}/{x}/{y}")
maxs <- cellStats(marseille.mb, max)
marseille.mb <- projectRaster(from=marseille.mb, crs=st_crs(3035)$proj4string) # la projection fait un truc bizarre sur les entiers
marseille.mb <- marseille.mb/cellStats(marseille.mb, max)*maxs %>% as.integer # on remet tout comme avant mais en 3035
iso_marseille <- lload_DVF("iso_transit_50_r5_Marseille")
marseille <- iso2time(iso_marseille$EMP09, seuils=c(50000, 250000))
m <- tm_shape(marseille.mb, bbox = bb759)+tm_rgb()+tm_shape(marseille)+tm_raster(style="cont", palette=terrain.colors(20, rev=FALSE), breaks=c(20,40,60,80,100))

mmm <- tmap_arrange(p, l , m, nrow=3)
tmap_save(tm=mmm, "{DVFdata}/presentation/vv/plm transit.svg" %>% glue)
