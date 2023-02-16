#============================================================#
# MEAPS avec file d'attente et hétérogénité absorption/fuite #
# génération des graphiques et tableaux pour le .qmd         #
#============================================================#

# on utilise maintenant meaps (avec pl et odd ratio)
# le c++ est dans meaps_rcpp
# ggplot 3.4 fusille les binhex
# devtools::install_version("ggplot2", version = "3.3.6")
# init -------------------------
library(tidyverse)
library(conflicted)
library(furrr)
source("R/radiation functions.r")
Rcpp::sourceCpp("R/meaps_rcpp.cpp", echo=FALSE)
future::plan("multisession", workers = 8)
library(tictoc)
library(ggnewscale)
library(matrixStats)
library(patchwork)
library(ofce)
library(Rcpp)
library(progressr)
library(scales)
options(ofce.background_color = "grey97")
showtext::showtext_opts(dpi = 92)
showtext::showtext_auto()
conflict_prefer_all("dplyr", quiet = TRUE)
sysfonts::font_add_google('Nunito')
options(ofce.base_family = "Nunito")
options(ofce.base_size = 9)
plan("multisession", workers = 8)
Rcpp::sourceCpp("R/meaps2.cpp", echo=FALSE)
source("R/meaps2.r")

handlers(global = TRUE)
handlers("cli")

n <- 5000
k <- 5000
bins <- 1.2/0.05

## génération --------------
# on représente ici une agglomération centrale, plus des villages avec des emplois
# mais pas en nombre suffisant dans les villages
set.seed(1942) 
habc <- cbind(pos_cnorm(n=70*k/100, sigma = 0.2, centre = c(1, 1)), f=0.1, g = 1)
habv1 <- cbind(pos_cnorm(n=15*k/100, sigma = 0.1, centre = c(.6, .6)), f=0.1, g = 3)
habv12 <- cbind(pos_cnorm(n=15*k/100, sigma = 0.1, centre = c(.1, 0.1)), f=0.1, g = 3)
habv2 <- cbind(pos_cnorm(n=15*k/100, sigma = 0.1, centre = c(1.4, 1.4)), f=0.1, g = 2)
hab <- rbind(habc, habv2, habv1)
hab2 <- rbind(habc, habv2, habv12)
empc <- cbind(pos_cnorm(n=80/100*k, sigma = 0.05, centre = c(1, 1)), p=1, g=1)
empv1 <- cbind(pos_cnorm(n=5/100*k, sigma = 0.05, centre = c(0.6, 0.6)), p=1, g=3)
empv12 <- cbind(pos_cnorm(n=5/100*k, sigma = 0.05, centre = c(.1, .1)), p=1, g=3)
empv2 <- cbind(pos_cnorm(n=5/100*k, sigma = 0.05, centre = c(1.4, 1.4)), p=1, g=2)
emp <- rbind(empc, empv2, empv1)
emp2<- rbind(empc, empv2, empv12)


s1 <- make_tibs(emp, hab)
s2 <- make_tibs(emp2, hab2)
# save(s1, s2, file = "radiation/graphs/scenarios.rda")

dds <- rdist::cdist(cbind(s1$hgroupes$x, s1$hgroupes$y), cbind(s1$egroupes$x, s1$egroupes$y)) 
dds2 <- rdist::cdist(cbind(s2$hgroupes$x, s2$hgroupes$y), cbind(s2$egroupes$x, s2$egroupes$y)) 
rownames(dds) <- c("h1", "h2", "h3")
colnames(dds) <- c("e1", "e2", "e3")
rownames(dds2) <- c("h1", "h2", "h3")
colnames(dds2) <- c("e1", "e2", "e3")

dds |>
  knitr::kable(digits = 1)
dds2 |>
  knitr::kable(digits = 1)
save(dds, dds2, file="output/dds.rda")
meann <- function(n) function(x) ifelse(length(x)>n, mean(x), NA)
### cartes ----------------

coords <- coord_equal(xlim=c(0,2), ylim=c(0,2))
(gcarte_ss <- 
    (ggplot()+
       stat_binhex(data = as_tibble(s1$hab),
                   aes(x=x,y=y, fill=100*after_stat(density)), binwidth=0.05)+
       scale_fill_distiller(palette="Greens", direction=1, name = "densité\nd'habitants")+
       coords +
       geom_text(data = s1$hgroupes, aes(x=x, y=y, label = g_label), nudge_y = 0.4,  size = 2) +
       labs(title = "Habitants")+
       theme_void(base_size = 8)+ 
       theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
             plot.margin = margin(6,6,6,6),
             panel.background = element_rect(fill="grey97")))+
    (ggplot()+
       stat_binhex(data=as_tibble(s1$emp),
                   aes(x=x,y=y, fill=100*after_stat(density)), binwidth=0.05)+
       scale_fill_distiller(palette = "Oranges", direction=1, name = "densité\nd'emplois")+
       coords +
       geom_text(data = s1$egroupes, aes(x=x, y=y, label = g_label), size = 2, nudge_y = 0.4) +
       labs(title = "Emplois")+
       theme_void(base_size = 8)+ 
       theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
             plot.margin = margin(6,6,6,6),
             panel.background = element_rect(fill="grey97"))) + 
    plot_layout(guides = 'collect'))
(gcarte_ss2 <- 
    (ggplot()+
       stat_binhex(data = as_tibble(s2$hab),
                   aes(x=x,y=y, fill=100*after_stat(density)), binwidth=0.05)+
       scale_fill_distiller(palette="Greens", direction=1, name = "densité\nd'habitants")+
       coords +
       geom_text(data = s2$hgroupes, aes(x=x, y=y, label = g_label), nudge_y = 0.4,  size = 2) +
       labs(title = "Habitants")+
       theme_void(base_size = 8)+ 
       theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
             plot.margin = margin(6,6,6,6),
             panel.background = element_rect(fill="grey97")))+
    (ggplot()+
       stat_binhex(data=as_tibble(s2$emp),
                   aes(x=x,y=y, fill=100*after_stat(density)), binwidth=0.05)+
       scale_fill_distiller(palette = "Oranges", direction=1, name = "densité\nd'emplois")+
       coords +
       geom_text(data = s2$egroupes, aes(x=x, y=y, label = g_label), size = 2, nudge_y = 0.4) +
       labs(title = "Emplois")+
       theme_void(base_size = 8)+ 
       theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
             plot.margin = margin(6,6,6,6),
             panel.background = element_rect(fill="grey97"))) + 
    plot_layout(guides = 'collect'))
graph2png(gcarte_ss, rep="output", ratio = 2)
graph2png(gcarte_ss2, rep="output", ratio = 2)
save(gcarte_ss, file = "output/gcarte_ss.rda")
save(gcarte_ss2, file = "output/gcarte_ss2.rda")

# calculs ----------------------------

# bench::mark(v1 = meaps_cpp(s1$rk,f = s1$f, p = s1$p, shuf = 1:k))
# bench::mark(v2 = meaps_rcpp(s1$rk,emplois=rep(1, n), actifs = rep(1, k), odds = s1$p, f = s1$f, shuf = 1:k))
shufs <- do.call(rbind, purrr::map(1:100, ~sample.int(n,n)))
# la v2
# mm <- rmeaps(emp = emp, hab = hab, shuf = shufs, meaps_ver = 2)
tic();mmb <- rmeaps_bstp(s1, shufs, workers = 6); toc()
# la v3
tic();mmb2 <- rmeaps_bootmat(s1$rk, 
                             emplois = rep(1, k*0.9), 
                             actifs = rep(1, n), 
                             f = rep(0.1, n),
                             modds = matrix(1, nrow = n, ncol = k*0.9),
                             shuf = shufs, workers = 8); toc()
# la v4 qui doit être la finale
library(rmeaps)
tic(); mmnm <- rmeaps::meaps_bootstrap2(
  s1$rk, 
  emplois = rep(1, k*0.9), 
  actifs = rep(1, n), 
  f = rep(0.1, n),
  modds = matrix(1, nrow = n, ncol = k*0.9),
  shuf = shufs); toc()
# Les variations
#mm2 <- rmeaps(emp = emp2, hab = hab2, shuf = shufs, meaps_ver = 2)
mmb2 <- rmeaps_bstp(s2, shufs, workers = 6)
tic(); mmnm2 <- rmeaps::meaps_bootstrap2(
  s2$rk, 
  emplois = rep(1, k*0.9), 
  actifs = rep(1, n), 
  f = rep(0.1, n),
  modds = matrix(1, nrow = n, ncol = k*0.9),
  shuf = shufs); toc()


# save(mm, mm2, file = "radiation/graphs/mms.rda")

## matrice de flux ----------------
flux <- emp_flux(s1, mmb$emps)$s |>
  as_tibble() |>
  rename_all(~str_c("e",.x)) |> 
  mutate(gh = str_c("h", s1$hgroupes$g)) |>
  relocate(gh) |> 
  add_total() |> 
  rowwise() |> 
  mutate(total = sum(c_across(2:4))) |> 
  ungroup() |> 
  mutate(across(2:5, ~prettyNum(round(.x), format ="d", big.mark = " ")))

flux_new <- emp_flux(s1, mmnm)$s |>
  as_tibble() |>
  rename_all(~str_c("e",.x)) |> 
  mutate(gh = str_c("h", s1$hgroupes$g)) |>
  relocate(gh) |> 
  add_total() |> 
  rowwise() |> 
  mutate(total = sum(c_across(2:4))) |> 
  ungroup() |> 
  mutate(across(2:5, ~prettyNum(round(.x), format ="d", big.mark = " ")))

### flux2 ---------------
flux2 <- emp_flux(s2, mmb2$emps)$s |>
  as_tibble() |>
  rename_all(~str_c("e",.x)) |> 
  mutate(gh = str_c("h", s1$hgroupes$g)) |>
  relocate(gh) |> 
  add_total() |> 
  rowwise() |> 
  mutate(total = sum(c_across(2:4))) |> 
  ungroup() |> 
  mutate(across(2:5, ~prettyNum(round(.x), format ="d", big.mark = " ")))
knitr::kable(flux)
knitr::kable(flux2)

save(flux, flux2, file = "output/tblflux.rda")

(gdenshabg <- ggplot(mm$hab)+
    geom_density(aes(x=d, group=g, fill=factor(g), col=factor(g)), alpha = 0.5)+
    geom_density(data = mm2$hab, aes(x=d, group=g, fill=factor(g), col=factor(g)), alpha = 0.25, linetype ="dashed")+
    scico::scale_color_scico_d(
      palette ="tofino",
      begin=0.25,
      name="pôle d'habitation",
      aesthetics = c('color','fill'), 
      labels = c("h1", "h2", "h3"))+
    xlim(c(0,1.5))+ylim(c(0,5))+ylab(NULL)+xlab("distance")+
    theme_ofce(base_size=9, base_family = "Nunito") +theme(legend.position = "right"))
graph2png(gdenshabg, rep="/output", ratio = 16/9)
save(gdenshabg, file = "/output/gdenshabg.rda")

# flux |> gt() |> 
#   gt::fmt_integer(columns = 2:5,
#                   rows= everything(), sep_mark = " ") |> 
#   gt::summary_rows(columns = 2:4, fns = list(total = ~sum(.)), formatter = fmt_integer, sep_mark = " ") |> 
#   gt::cols_width(everything() ~ px(70))

## graphes à densité/distances moyennes -------------
### s1 -----------------
gdhab <- ggplot()+
  stat_summary_hex(data = mm$hab, aes(x=x, y=y, z=d), fun = meann(5), bins=25)+
  guides(fill=guide_colourbar("Distance\npar habitant"))+
  scale_fill_distiller(palette="Greens", direction=1)+
  xlab(NULL)+ylab(NULL)+labs(title="Distance moyenne parcourue par habitant")+
  geom_text(data = s1$hgroupes, aes(x=x, y=y, label = g_label), nudge_y = -0.2,  size = 2)+
  coords +
  theme_void(base_size = 8) +
  theme(legend.position="right", 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold",  margin = margin(b=4)),
        plot.margin = margin(2,2,2,2),
        panel.background = element_rect(fill="grey97"))
gdenshab <- ggplot(mm$hab)+
  geom_density(aes(x=d), fill="green", col=NA, alpha = 0.5)+
  xlim(c(0,1.5))+ylim(c(0,5))+ylab(NULL)+xlab(NULL)+
  theme_minimal()+
  theme(text=element_text(size=6))
gdemp <- ggplot()+
  stat_summary_hex(data = mm$emp, aes(x=x, y=y, z=d), fun = meann(5), bins=25)+
  guides(fill=guide_colourbar("Distance\npour un emploi"))+
  scale_fill_distiller(palette="Oranges", direction=1)+
  xlab(NULL)+ylab(NULL)+labs(title="Distance moyenne parcourue pour un emploi") +
  geom_text(data = s1$egroupes, aes(x=x, y=y, label = g_label), size = 2, nudge_y = -0.2) +
  coords +
  theme_void(base_size = 8)+ 
  theme(legend.position="right", 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
        plot.margin = margin(2,2,2,2),
        panel.background = element_rect(fill="grey97"))
gdensemp <- ggplot(mm$emp)+
  geom_density(aes(x=d), fill="orange", col=NA, alpha = 0.5, position="stack")+
  xlim(c(0,1.5))+
  scale_y_continuous(limits = c(0,10), oob =squish)+ylab(NULL)+xlab(NULL)+
  theme_minimal()+
  theme(text=element_text(size=6))

(gdistances <- 
    (gdhab+inset_element(gdenshab, 0.5, 0., 1, 0.33))+
    (gdemp+inset_element(gdensemp, 0.5, 0., 1, 0.33))+
    plot_layout(guides = "collect"))

graph2png(gdistances, rep="output", ratio = 2)
save(gdistances, file = "output/gdistances.rda")

### s2 ---------------------
gdhab <- ggplot()+
  stat_summary_hex(data = mm2$hab, aes(x=x, y=y, z=d), fun = meann(5), bins=25)+
  guides(fill=guide_colourbar("Distance\npar habitant"))+
  scale_fill_distiller(palette="Greens", direction=1)+
  xlab(NULL)+ylab(NULL)+labs(title="Distance moyenne parcourue par habitant")+
  geom_text(data = s2$hgroupes, aes(x=x, y=y, label = g_label), nudge_y = -0.2,  size = 2)+
  coords +
  theme_void(base_size = 8) +
  theme(legend.position="right", 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold",  margin = margin(b=4)),
        plot.margin = margin(2,2,2,2),
        panel.background = element_rect(fill="grey97"))
gdenshab <- ggplot(mm2$hab)+
  geom_density(aes(x=d), fill="green", col=NA, alpha = 0.5)+
  xlim(c(0,1.5))+ylim(c(0,5))+ylab(NULL)+xlab(NULL)+
  theme_minimal()+
  theme(text=element_text(size=6))
gdemp <- ggplot()+
  stat_summary_hex(data = mm2$emp, aes(x=x, y=y, z=d), fun = meann(5), bins=25)+
  guides(fill=guide_colourbar("Distance\npour un emploi"))+
  scale_fill_distiller(palette="Oranges", direction=1)+
  xlab(NULL)+ylab(NULL)+labs(title="Distance moyenne parcourue pour un emploi") +
  geom_text(data = s2$egroupes, aes(x=x, y=y, label = g_label), size = 2, nudge_y = -0.2) +
  coords +
  theme_void(base_size = 8)+ 
  theme(legend.position="right", 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
        plot.margin = margin(2,2,2,2),
        panel.background = element_rect(fill="grey97"))
gdensemp <- ggplot(mm2$emp)+
  geom_density(aes(x=d), fill="orange", col=NA, alpha = 0.5, position="stack")+
  xlim(c(0,1.5))+
  scale_y_continuous(limits = c(0,10), oob = squish)+ylab(NULL)+xlab(NULL)+
  theme_minimal()+
  theme(text=element_text(size=6))

(gdistances2 <- 
    (gdhab+inset_element(gdenshab, 0.5, 0., 1, 0.33))+
    (gdemp+inset_element(gdensemp, 0.5, 0., 1, 0.33))+
    plot_layout(guides = "collect"))

graph2png(gdistances2, rep="meaps-doc/output", ratio = 2)
save(gdistances2, file = "meaps-doc/output/gdistances2.rda")

## modèle gravitaire ----------------------

flux_grav <- function(s, delta) {
  f <-  exp(-s$dist/delta)
  f1 <- f/rowSums2(f)*(1-s$f)
  return(f1)
}

flux_grav_furness <- function(s, delta, tol=0.0001) {
  f <-  exp(-s$dist/delta) 
  err <- 1
  while(err>tol) {
    f1 <- f * (1-s$f)/rowSums2(f)
    f2 <- t(t(f1)/colSums2(f1))
    err <- sqrt(mean((rowSums2(f2)-(1-s$f))^2))
    f <- f2
  }
  return(f)
}

source("v2/annexes/f.normalisation.r")
score_grav <- function(emps, s, delta, tol = 0.0001, furness = FALSE) {
  if(furness)
    f <- flux_grav_furness(s, delta, tol)
  else
    f <- flux_grav(s, delta)
  sum((emp_flux(s1, emps)$s - emp_flux(s1, f)$s)^2)
}

kl_grav <- function(emps, s, delta, tol = 0.0001, furness=FALSE) {
  if(furness)
    f <- flux_grav_furness(s, delta, tol)
  else
    f <- flux_grav(s, delta)
  kullback_leibler(c(emp_flux(s1, f)$s), c(emp_flux(s1, emps)$s))
}

fr2 <- optim(.5, \(x) score_grav(mmb$emps,s1,x), lower = 0.001, upper = 10, method = "L-BFGS-B")
fkl <- optim(.5, \(x) kl_grav(mmb$emps,s1,x, furness=TRUE), lower = 0.001, upper = 10, method = "L-BFGS-B")
fkl2 <- optim(.5, \(x) kl_grav(mmb2$emps,s2,x, furness=TRUE), lower = 0.001, upper = 10, method = "L-BFGS-B")
fref <- flux_grav_furness(s1, fkl$par)
f2 <- flux_grav_furness(s2, fkl$par)
f22 <- flux_grav_furness(s2, fkl2$par)
r2kl(c(emp_flux(s1, fref)$s), c(emp_flux(s1, mm$meaps)$s))
fluxg <- emp_flux(s1, fref)$s |>
  as_tibble() |>
  rename_all(~str_c("e",.x)) |> 
  mutate(gh = str_c("h", s1$hgroupes$g)) |>
  relocate(gh) |> 
  add_total() |> 
  rowwise() |> 
  mutate(total = sum(c_across(2:4))) |> 
  ungroup() |> 
  mutate(across(2:5, ~prettyNum(round(.x), format ="d", big.mark = " ")))

fluxg2 <- emp_flux(s2, f2)$s |>
  as_tibble() |>
  rename_all(~str_c("e",.x)) |> 
  mutate(gh = str_c("h", s1$hgroupes$g)) |>
  relocate(gh) |> 
  add_total() |> 
  rowwise() |> 
  mutate(total = sum(c_across(2:4))) |> 
  ungroup() |> 
  mutate(across(2:5, ~prettyNum(round(.x), format ="d", big.mark = " ")))

fluxg22 <- emp_flux(s2, f22)$s |>
  as_tibble() |>
  rename_all(~str_c("e",.x)) |> 
  mutate(gh = str_c("h", s1$hgroupes$g)) |>
  relocate(gh) |> 
  add_total() |> 
  rowwise() |> 
  mutate(total = sum(c_across(2:4))) |> 
  ungroup() |> 
  mutate(across(2:5, ~prettyNum(round(.x), format ="d", big.mark = " ")))

save(fluxg, fluxg2, fluxg22, fkl, fkl2, file = "meaps-doc/output/flux_grav.srda")

## ergodicité ------------------------
# attention on utilise ici la version 1 de meaps
# ce qu'il y a en moins : pas de traitement des paquets, pas de traitement des odds
# pas de débordement
# mais dans ce cas ce n'est pas grave
if(fs::file_exists("meaps-doc/output/libres.raw.qs")) {
  libres.raw <- qs::qread("meaps-doc/output/libres.raw.qs") 
} else {
  plan("multisession", workers = 7)
  
  tic();
  libres.raw <- future_imap_dfr(1:500, ~{
    Rcpp::sourceCpp("radiation/meaps.cpp", showOutput = FALSE, echo = FALSE)
    shuf <- sample.int(k,k)
    raw <- meaps_cpp(s1$rk,f = s1$f, p = s1$p, shuf = shuf)
    colnames(raw$emps) <- str_c("e", 1:ncol(raw$emps))
    colnames(raw$dispo) <- str_c("e", 1:ncol(raw$dispo))
    tibrn <- tibble(emp = 1:ncol(raw$papn),
                    ehex = s1$emps$ehex[emp],
                    rangn = apply(raw$papn<0.01, FUN = which.max, MARGIN = 2)) |> 
      group_by(ehex) |> 
      summarize(rs = sum(rangn),
                rs2 = sum(as.numeric(rangn)^2),
                rn = n())
    tibe <- as_tibble(raw$emps)
    names(tibe) <- str_c("e", 1:ncol(raw$emps))
    tibe <- tibe |> 
      mutate(hab = 1:nrow(tibe)) |> 
      relocate(hab) |> 
      mutate(hhex = s1$habs$hhex[hab],
             hg = s1$habs$g[hab]) |> 
      pivot_longer(cols = starts_with("e"), values_to = 'pemp', names_to = "emp") |> 
      mutate(
        emp = as.numeric(str_sub(emp,2,-1)),
        ehex  = s1$emps$ehex[emp],
        eg = s1$emps$g[emp]) |>
      group_by(hhex, ehex) |>
      summarize(
        ps = sum(pemp),
        ps2 = sum(pemp^2),
        pn = n(), .groups = "drop") |>
      left_join(s1$hhex |> select(hhex, gh), by="hhex") |> 
      left_join(s1$ehex |> select(ehex, ge), by="ehex")
    tibo <- as_tibble(raw$dispo[-1,])
    names(tibo) <- str_c("e", 1:ncol(raw$dispo))
    tibo <- tibo |> 
      mutate(hab = 1:nrow(tibo)) |> 
      relocate(hab) |> 
      mutate(hhex = s1$hexhab[hab]) |> 
      pivot_longer(cols = starts_with("e"), values_to = 'libre', names_to = "emp") |> 
      mutate(
        emp = as.numeric(str_sub(emp,2,-1)),
        ehex  = s1$hexemp[emp]) |>
      group_by(hhex, ehex) |>
      summarize(ls = sum(libre),
                ls2 = sum(libre^2),
                ln = n(), 
                .groups = "drop")
    left_join(tibe, tibo, by = c('hhex', 'ehex')) |> 
      left_join(tibrn, by="ehex") |> 
      mutate(draw = .y) |> 
      relocate(hhex, ehex, draw)
  }, .options = furrr_options(seed = TRUE), .progress = TRUE)
  toc();
  
  qs::qsave(libres.raw, "meaps-doc/output/libres.raw.qs")
}

libres.raw <- qs::qread("meaps-doc/output/libres.raw.qs")

erg_libre <- libres.raw |>
  group_by(hhex, ehex) |> 
  arrange(draw) |> 
  mutate(
    rnc = cumsum(as.numeric(rn)), rm = cumsum(rs)/rnc, rsd = sqrt(cumsum(rs2)/rnc - rm^2),
    lnc=cumsum(as.numeric(ln)), lm=cumsum(ls)/lnc, lsd = sqrt(cumsum(ls2)/lnc -  lm^2),
    pnc=cumsum(as.numeric(pn)), pm=cumsum(ps)/pnc, psd = sqrt(cumsum(ps2)/pnc -  pm^2)) |> 
  ungroup()

# test 1 de tous les hab vers quelques emp
(gemploi_erg <- ggplot(erg_libre |> filter(ehex%in%sample(unique(ehex), 4), draw<100)) + 
    geom_line(aes(x=draw, y=ps/pn, group = hhex), col="white", alpha=0.1, size = 0.1) +
    geom_line(aes(x=draw, y=pm, group = hhex), col="green", alpha=0.1, size = 0.1) +
    scale_y_log10()+
    facet_wrap(vars(ehex))+
    theme_ofce()+
    xlab("Nombre de tirages")+
    ylab("part d'emploi affectée au carreau d'emploi")+
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill="black"),
          plot.background = element_rect(fill="white"),
          text = element_text(color = "black"),
          axis.text = element_text(color = "black"),
          strip.text = element_text(color = "black")))
graph2png(gemploi_erg, rep="meaps-doc/output", ratio=4/3)
save(gemploi_erg, file="meaps-doc/output/gemploi_erg.rda")
# test 2 vers tous les emp de quelques hab
gemploi_erg <- ggplot(erg_libre |> filter(hhex%in%sample(unique(hhex), 9))) + 
  geom_line(aes(x=draw, y=ps/pn, group = ehex), col="white", alpha=0.05, size=0.1) +
  geom_line(aes(x=draw, y=pm, group = ehex), col="green", alpha=0.1, size=0.2) +
  facet_wrap(vars(hhex))+
  scale_y_log10()+
  theme_ofce()+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="black"),
        plot.background = element_rect(fill="white"),
        text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black"))

erg_libre_emp <- erg_libre |> 
  group_by(ehex, hhex) |>
  ungroup() |> 
  filter(ehex%in%sample(unique(ehex), 3))

# test 3 pour quelques emp, état libre de i
ggplot(erg_libre_emp) + 
  geom_line(aes(x=draw, y=ls/ln, group = hhex), col="white", alpha=0.05, size=0.05) +
  geom_line(aes(x=draw, y=lm, group = hhex), col="green", alpha=0.1, size=0.05) +
  facet_wrap(vars(ehex))+
  theme_ofce()+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="black"),
        plot.background = element_rect(fill="white"),
        text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black"))

# test 4

### matrice de flux ---------------
### 
shufs <- purrr::map(1:500, ~sample.int(n,n))
plan("multisession", workers=16)
fluxs <- future_imap_dfr(shufs, ~{
  Rcpp::sourceCpp("radiation/meaps2.cpp", echo=FALSE)
  mm <- rmeaps(emp = emp, hab = hab, shuf=.x , meaps_ver = 2)
  emp_flux(s1, mm$meaps)$s |>
    as_tibble() |>
    rename_all(~str_c("e",.x)) |> 
    mutate(gh = str_c("h", s1$hgroupes$g)) |>
    relocate(gh) |> 
    mutate(draw = .y)
}, .options = furrr_options(seed = TRUE), .progress = TRUE)
save(fluxs, file="meaps-doc/output/fluxs.rda")
load("meaps-doc/output/fluxs.rda")
gfluxs <- ggplot(fluxs)+
  geom_boxplot(aes(y=flux, x=cat, color=gh),
               show.legend = FALSE)+
  scale_y_log10()+
  xlab(NULL)+ylab("Nombre de trajets")+
  theme_ofce()
graph2png(gfluxs, rep="meaps-doc/output")

fluxsq <- fluxs |>
  pivot_longer(cols = starts_with("e"), names_to = "ge", values_to = "flux") |> 
  group_by(ge, gh) |> 
  summarize(q5 = quantile(flux, 0.05),
            q50 = quantile(flux, 0.5),
            q95 = quantile(flux, 0.95)) |> 
  transmute(ge, gh,
            intervale_r = 
              str_c(round(q50), "<br>[", 
                    round(q5), "; ",
                    round(q95), "]")) |> 
  pivot_wider(id_cols = gh, names_from = ge, values_from = intervale_r)

save(fluxsq, file="meaps-doc/output/fluxsq.srda")

# test 5 rang n
rangns <- erg_libre |> group_by(ehex, draw) |> summarize(across(c(rm, rsd, ge), first))
rangns_m <- ggplot(rangns) + 
  geom_line(aes(x=draw, y=rm, group = ehex), col="white", alpha=0.3, size=0.2) +
  theme_ofce()+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="black", color = NA),
        plot.background = element_rect(fill="white", color = NA),
        text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black"))+
  xlab("Nombre de tirages cumulés")+ylab("rang moyen au moment du remplissage")

rangns_ec <- ggplot(rangns) + 
  geom_line(aes(x=draw, y=rsd, group = ehex), col="white", alpha=0.25, size=0.2) +
  theme_ofce()+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="black", color = NA),
        plot.background = element_rect(fill="white", color = NA),
        text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black"))+
  xlab("Nombre de tirages cumulés")+ylab("écart type du rang moyen au moment du remplissage")

library(patchwork)
g_rangns <- rangns_m + rangns_ec
graph2png(g_rangns, rep="meaps-doc/output", ratio = 16/7)

gnmsd <- ggplot(rangns |> filter(draw==500))+
  geom_point(aes(x=rm, y=rsd, color=ge)) +
  xlab("rang moyen au moment du remplissage")+
  ylab("écart type du rang moyen au moment du remplissage") + 
  scale_color_discrete(name = "Pôle d'emploi")+
  theme_ofce()+theme(legend.position="right")
graph2png(gnmsd, rep="radiation/svg")
save(gnmsd, file="meaps-doc/output/gmsd_erg.rda")


# martrice de dérivées -----------------------------
sr <- s1
n <- nrow(sr$habs)
k <- nrow(sr$emps)
gh <- unique(sr$habs$g) |> sort()
ge <- unique(sr$emps$g) |> sort()
groups <- cross2(ge,gh)

# groups <- list(h1g1 = list(1,1), h2g2 = list(2,2))
names(groups) <- map_chr(groups, ~str_c("e", .x[[1]], "-h", .x[[2]]))
ms_odd <- list_modify( 
  list(ref = matrix(1, nrow = n, ncol = k)),
  !!!imap(groups,  ~{
    m <- matrix(1, nrow = n, ncol = k) 
    m[which(sr$habs$g==.x[[2]]), which(sr$emp$g==.x[[1]])] <- 2
    m
  }))

shufs <- do.call(rbind, map(1:100, ~sample.int(n,n)))

plan("multisession", workers =6)

res <- future_imap(
  ms_odd, ~{
    sourceCpp("R/meaps_oddmatrix.cpp", showOutput = FALSE, echo = FALSE)
    mm <- meaps_boot2(sr$rk,
                      emplois = rep(1, n), actifs = rep(1, n),
                      f = sr$f,
                      modds = .x , shufs = shufs)
    flux <- emp_flux(sr, mm$emps, mm$emps2)
    ed <- mm$emps * sr$dist
    emp_i <- matrixStats::rowSums2(mm$emps)
    d_ind <- matrixStats::rowSums2(ed)/emp_i
    list(habs = bind_cols(sr$habs, tibble(d = d_ind), scn = .y),
         flux = flux$s,
         flux_ec = flux$ec,
         emps = mm$emps,
         empec = sqrt(mm$emps2- mm$emps^2))
  }, .options = furrr_options(seed=TRUE), .progress = TRUE) |> purrr::transpose()
tres <- bind_rows(res$habs)

mmflux <- map(res$flux[-1], ~.x - res$flux$ref)
save(mmflux, file="meaps-doc/output/mmflux.rda")
ee <- map(mmflux, as.vector) |> flatten_dbl() |> matrix(ncol=9, nrow=9) |> eigen() 
ee$values |> round(1)

bigflux <- bind_rows(
  imap(mmflux, ~{
    as_tibble(.x) |>
      mutate(h=str_c("h", 1:3),
             odd = .y) |>
      rename("e1"=`1`,"e2"=`2`, "e3"=`3`)})) |> 
  separate(odd, c("oe", "oh"), sep="-") |>
  mutate(oe = str_replace(oe, "e", "e")) |> 
  pivot_wider(id_cols = c(h, oh),
              names_from = oe, 
              names_glue = "{oe}_{.value}",
              values_from = c(e1,e2,e3)) 

bigflux <- bigflux |> 
  relocate(sort(names(bigflux))) |> 
  relocate(oh, h)

library(gt)
flux3x3 <- bigflux |> 
  group_by(oh) |> 
  gt() |>
  tab_spanner_delim(delim="_", split = "last") |>   
  cols_label(oh ="", h="") |> 
  tab_options(row_group.as_column = TRUE) |> 
  ofce::table_ofce() |> 
  fmt_number(columns = where(is.numeric),
             decimals = 0) |> 
  tab_style(
    style = cell_borders(
      side = "left", weight = px(1), color = "black"),
    locations = cells_body(
      columns = c(e1_e1, e2_e1, e3_e1))) |> 
  tab_style(
    style = cell_borders(
      side = "right", weight = px(1), color = "black"),
    locations = cells_body(
      columns = c(e3_e3))) |> 
  tab_style(
    style = cell_borders(
      side = c("bottom"), weight = px(1), color = "black"),
    locations = cells_body(
      rows = c(3,6,9)
    )) 

save(flux3x3, file="meaps-doc/output/flux3x3.rda" |> glue::glue())  

(gcarte_cfg2 <- 
    (ggplot()+
       stat_binhex(data = as_tibble(sr$hab),
                   aes(x=x,y=y, fill=100*after_stat(density)), binwidth=0.025)+
       scale_fill_distiller(palette="Greens", direction=1, name = "densité\nd'habitants")+
       coord_equal(xlim=c(0,1), ylim=c(0,1))+
       geom_text(data = sr$hgroupes, aes(x=x, y=y, label = g_label),
                 nudge_y = 0,  size = 2, hjust = 0.5, vjust = 0.5) +
       labs(title = "Habitants")+
       theme_void(base_size = 8)+ 
       theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
             plot.margin = margin(6,6,6,6),
             panel.background = element_rect(fill="grey97")))+
    (ggplot()+
       stat_binhex(data=as_tibble(sr$emp),
                   aes(x=x,y=y, fill=100*after_stat(density)), binwidth=0.025)+
       scale_fill_distiller(palette = "Oranges", direction=1, name = "densité\nd'emplois")+
       coord_equal(xlim=c(0,1), ylim=c(0,1))+
       geom_text(data = sr$egroupes, aes(x=x, y=y, label = g_label), size = 2, nudge_y = -0.2) +
       labs(title = "Emplois")+
       theme_void(base_size = 8)+ 
       theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
             plot.margin = margin(6,6,6,6),
             panel.background = element_rect(fill="grey97"))) + 
    plot_layout(guides = 'collect'))

tres <- tres |>
  separate(scn, into = c("e", "h"), sep="-") |> 
  mutate(h = replace_na(h, "ref"))
ggplot(data = tres |> filter(e!="ref") )+
  geom_density(aes(x=d, col = e, fill=e), alpha = 0.3, position="identity", show.legend = FALSE)+
  geom_density(data=tres |> filter(e=="ref") |> select(-e,-h), 
               aes(x=d),
               col = "gray",
               fill="gray",
               alpha = 0.3,
               position="identity")+
  facet_grid(rows = vars(h), cols = vars(e)) +
  ylab(NULL)+xlab("distances")+
  theme_ofce(base_size=9)+theme(legend.position = "right")