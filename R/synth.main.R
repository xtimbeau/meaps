#============================================================#
# MEAPS avec file d'attente et hétérogénité absorption/fuite #
# génération des graphiques et tableaux pour le .qmd         #
#============================================================#

# on utilise maintenant rmeaps
# dans un package à installer avec devtools::install_github("maxime2506/rmeaps")
# Il faut é&galement installer devtools::install_github("OFCE/OFCE")

# init -------------------------
library(tidyverse)
library(conflicted)
library(furrr)
library(tictoc)
library(ggnewscale)
library(matrixStats)
library(patchwork)
library(ofce)
library(Rcpp)
library(progressr)
library(scales)
library(rmeaps)
library(qs)
options(ofce.background_color = "grey99",
        ofce.basefamily = "Source Sans Pro")
showtext::showtext_opts(dpi = 200)
showtext::showtext_auto()
conflict_prefer_all("dplyr", quiet = TRUE)
options(ofce.base_family = "Nunito")
options(ofce.base_size = 9)

source("R/radiation functions.r")

# les anciennes versions hors package
if(FALSE) {
  Rcpp::sourceCpp("R/meaps_rcpp.cpp", echo=FALSE)
  future::plan("multisession", workers = 8)
  Rcpp::sourceCpp("R/meaps2.cpp", echo=FALSE)
  source("R/meaps2.r")
}
handlers(global = TRUE)
handlers("cli")

n <- 5000
k <- 4500
bins <- 1.2/0.05
binwidth <- 0.05
## génération --------------
# on représente ici une agglomération centrale, plus des villages avec des emplois
# mais pas en nombre suffisant dans les villages

set.seed(1942) 
s1 <- genere_3p(n = n, k = k, rayon = 0.35,
                part_h = 0.7, part_e = 0.7, beta = 1.6,
                d_cp2 = 0.75, d_cp3 = 0.75,
                binwidth = binwidth)
set.seed(1942) 
s2 <- genere_3p(n = n, k = k, rayon = 0.35,
                part_h = 0.7, part_e = 0.7, beta = 1.6,
                d_cp2 = 1, d_cp3 = 1,
                binwidth = binwidth)
qsave(s1, file = "output/s1.sqs")
qsave(s2, file = "output/s2.sqs")
s1 <- add_dist(s1)
s2 <- add_dist(s2)
dds <- rdist::cdist(cbind(s1$hgroupes$x, s1$hgroupes$y), cbind(s1$egroupes$x, s1$egroupes$y)) 
dds2 <- rdist::cdist(cbind(s2$hgroupes$x, s2$hgroupes$y), cbind(s2$egroupes$x, s2$egroupes$y)) 
rownames(dds) <- c("h1", "h2", "h3")
colnames(dds) <- c("e1", "e2", "e3")
rownames(dds2) <- c("h1", "h2", "h3")
colnames(dds2) <- c("e1", "e2", "e3")

dds |>
  knitr::kable(digits = 2)
dds2 |>
  knitr::kable(digits = 2)
save(dds, dds2, file="output/dds.rda")
meann <- function(n) function(x) ifelse(length(x)>n, mean(x), NA)
### cartes ----------------

coords <- coord_equal(xlim=c(-1,1), ylim=c(-1,1))
(gcarte_ss <- 
    (ggplot()+
       stat_binhex(data = s1$habs,
                   aes(x=x,y=y, fill=100*after_stat(density)), binwidth=binwidth)+
       scale_fill_distiller(palette="Greens", direction=1, name = "densité\nd'habitants",
                            breaks = c(1,2))+
       coords +
       geom_text(data = s1$hgroupes, aes(x=x, y=y, label = g_label), nudge_y = 0.3,  size = 2) +
       labs(title = "Habitants")+
       theme_ofce_void(base_size = 9)+ 
       theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
             plot.margin = margin(6,6,6,6),
             panel.background = element_rect(fill="grey97")))+
    (ggplot()+
       stat_binhex(data=as_tibble(s1$emps),
                   aes(x=x,y=y, fill=100*after_stat(density)), binwidth=binwidth)+
       scale_fill_distiller(palette = "Oranges", direction=1, name = "densité\nd'emplois", 
                            breaks = c(1,2))+
       coords +
       geom_text(data = s1$egroupes, aes(x=x, y=y, label = g_label), size = 2, nudge_y = 0.2) +
       labs(title = "Emplois")+
       theme_ofce_void(base_size = 9)+ 
       theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
             plot.margin = margin(6,6,6,6),
             panel.background = element_rect(fill="grey97"))) + 
    plot_layout(guides = 'collect'))
(gcarte_ss2 <- 
    (ggplot()+
       stat_binhex(data = as_tibble(s2$hab),
                   aes(x=x,y=y, fill=100*after_stat(density)), binwidth=binwidth)+
       scale_fill_distiller(palette="Greens", direction=1, name = "densité\nd'habitants", 
                            breaks = c(1,2))+
       coords +
       geom_text(data = s2$hgroupes, aes(x=x, y=y, label = g_label), nudge_y = 0.3,  size = 2) +
       labs(title = "Habitants")+
       theme_ofce_void(base_size = 8)+ 
       theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
             plot.margin = margin(6,6,6,6),
             panel.background = element_rect(fill="grey97")))+
    (ggplot()+
       stat_binhex(data=as_tibble(s2$emp),
                   aes(x=x,y=y, fill=100*after_stat(density)), binwidth=binwidth)+
       scale_fill_distiller(palette = "Oranges", direction=1, name = "densité\nd'emplois",
                            breaks = c(1,2))+
       coords +
       geom_text(data = s2$egroupes, aes(x=x, y=y, label = g_label), size = 2, nudge_y = 0.2) +
       labs(title = "Emplois")+
       theme_ofce_void(base_size = 9)+ 
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
# mm <- rmeaps(emp = emp, hab = hab, shuf = shufs, meaps_ver = 2)
tic();mmb <- rmeaps_multishuf(s1, s1$shufs, nthreads=8); toc()

# Les variations
#mm2 <- rmeaps(emp = emp2, hab = hab2, shuf = shufs, meaps_ver = 2)
tic();mmb2 <- rmeaps_multishuf(s2, s2$shufs, nthreads=8); toc()

qs::qsave(mmb, file = "output/mmb1.bigqs", preset = "archive")
qs::qsave(mmb2, file = "output/mmb2.bigqs")

## matrice de flux ----------------
flux <- emp_flux(s1, mmb$meaps)$s |>
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
flux2 <- emp_flux(s2, mmb2$meaps)$s |>
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

(gdenshabg <- ggplot(mmb$hab)+
    geom_density(aes(x=d, y=after_stat(count)/nrow(mmb$hab), group=g, fill=factor(g), col=factor(g)),
                 alpha = 0.5, position = "stack")+
    geom_density(data = mmb2$hab,
                 aes(x=d, y=after_stat(count)/nrow(mmb$hab), group=g, fill=factor(g), col=factor(g)), 
                 alpha = 0.2, linetype ="dashed", position="stack")+
    scale_color_brewer(
      palette ="Dark2",
      name="pôle d'habitation",
      aesthetics = c('color','fill'), 
      labels = c("h1", "h2", "h3"))+
    xlim(c(0,1))+ylab(NULL)+xlab("distance")+
    theme_ofce(base_size=9) +theme(legend.position = "right"))
graph2png(gdenshabg, rep="output", ratio = 16/10)

 # flux |> gt() |> 
#   gt::fmt_integer(columns = 2:5,
#                   rows= everything(), sep_mark = " ") |> 
#   gt::summary_rows(columns = 2:4, fns = list(total = ~sum(.)), formatter = fmt_integer, sep_mark = " ") |> 
#   gt::cols_width(everything() ~ px(70))

## graphes à densité/distances moyennes -------------
### s1 -----------------
gdhab <- ggplot()+
  stat_summary_hex(data = mmb$hab, aes(x=x, y=y, z=d), fun = meann(5), binwidth=binwidth)+
  guides(fill=guide_colourbar("Distance\npar habitant"))+
  scale_fill_distiller(palette="Greens", direction=1,
                       breaks = c(.2, .3, .4, .5))+
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
gdenshab <- ggplot(mmb$hab)+
  geom_density(aes(x=d), fill="green", col=NA, alpha = 0.5)+
  xlim(c(0,1.5))+ylim(c(0,10))+ylab(NULL)+xlab(NULL)+
  theme_minimal()+
  theme(text=element_text(size=6))
gdemp <- ggplot()+
  stat_summary_hex(data = mmb$emp, aes(x=x, y=y, z=d), fun = meann(5), binwidth=binwidth)+
  guides(fill=guide_colourbar("Distance\npour un emploi"))+
  scale_fill_distiller(palette="Oranges", direction=1,
                       breaks = c(.2, .3, .4, .5))+
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
gdensemp <- ggplot(mmb$emp)+
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
gdhab2 <- ggplot()+
  stat_summary_hex(data = mmb2$hab, aes(x=x, y=y, z=d), fun = meann(5), binwidth=binwidth)+
  guides(fill=guide_colourbar("Distance\npar habitant"))+
  scale_fill_distiller(palette="Greens", direction=1)+
  xlab(NULL)+ylab(NULL)+labs(title="Distance moyenne parcourue par habitant")+
  geom_text(data = s2$hgroupes, aes(x=x, y=y, label = g_label), nudge_y = 0.3,  size = 2)+
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
gdenshab2 <- ggplot(mmb2$hab)+
  geom_density(aes(x=d), fill="green", col=NA, alpha = 0.5)+
  xlim(c(0,1.5))+ylim(c(0,10))+ylab(NULL)+xlab(NULL)+
  theme_minimal()+
  theme(text=element_text(size=6))
gdemp2 <- ggplot()+
  stat_summary_hex(data = mmb2$emp, aes(x=x, y=y, z=d), fun = meann(5), binwidth=binwidth)+
  guides(fill=guide_colourbar("Distance\npour un emploi"))+
  scale_fill_distiller(palette="Oranges", direction=1)+
  xlab(NULL)+ylab(NULL)+labs(title="Distance moyenne parcourue pour un emploi") +
  geom_text(data = s2$egroupes, aes(x=x, y=y, label = g_label), size = 2, nudge_y = 0.3) +
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
gdensemp2 <- ggplot(mmb2$emp)+
  geom_density(aes(x=d), fill="orange", col=NA, alpha = 0.5, position="stack")+
  xlim(c(0,1.5))+
  scale_y_continuous(limits = c(0,10), oob = squish)+ylab(NULL)+xlab(NULL)+
  theme_minimal()+
  theme(text=element_text(size=6))

(gdistances2 <- 
    (gdhab2+inset_element(gdenshab2, 0.5, 0., 1, 0.33))+
    (gdemp2+inset_element(gdensemp2, 0.5, 0., 1, 0.33))+
    plot_layout(guides = "collect"))

graph2png(gdistances2, rep="output", ratio = 2)
save(gdistances2, file = "output/gdistances2.rda")

## variance de la matrice de flux ---------------
### 
shufs <- purrr::map(1:256, ~sample.int(n,n))
plan("multisession", workers=8)
modds <- matrix(1, ncol=ncol(s1$rk), nrow = nrow(s1$rk))
fluxs <- future_imap_dfr(shufs, ~{
  mm <- rmeaps::meaps_oneshuf(
    rkdist = s1$rk, 
    emplois = rep(1,nrow(s1$emp)), 
    actifs = rep(1,nrow(s1$hab)),
    modds = modds,
    f = s1$f,
    shuf = .x)
  
  emp_flux(s1, mm)$s |>
    as_tibble() |>
    rename_all(~str_c("e",.x)) |> 
    mutate(gh = str_c("h", s1$hgroupes$g)) |>
    relocate(gh) |> 
    mutate(draw = .y)
}, .options = furrr_options(seed = TRUE), .progress = TRUE)

save(fluxs, file="output/fluxs.rda")
load("output/fluxs.rda")
gfluxs <- ggplot(fluxs |> pivot_longer(cols=-c(gh, draw), names_to = "ge", values_to = "flux"))+
  geom_boxplot(aes(y=flux, x=interaction(gh, ge), color=gh),
               show.legend = FALSE)+
  scale_y_log10()+
  xlab(NULL)+ylab("Nombre de trajets")+
  theme_ofce()
graph2png(gfluxs, rep="output")

fluxsq <- fluxs |>
  pivot_longer(cols = starts_with("e"), names_to = "ge", values_to = "flux") |> 
  group_by(ge, gh) |> 
  summarize(q5 = quantile(flux, 0.025),
            q50 = quantile(flux, 0.5),
            q95 = quantile(flux, 0.975),
            te = sqrt(sd(flux))/mean(flux)) |> 
  transmute(ge, gh,
            intervale_r = 
              str_c(round(q50), "<br>[", 
                    round(q5), "; ",
                    round(q95), "]<br>",
                    "\u3B5\u2248", round(te*100,1), "%")) |> 
  pivot_wider(id_cols = gh, names_from = ge, values_from = intervale_r)

save(fluxsq, file="output/fluxsq.srda")


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

source("R/f.normalisation.r")
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
  kl(c(emp_flux(s1, f)$s), c(emp_flux(s1, emps)$s))
}

fr2 <- optim(.5, \(x) score_grav(mmb$meaps,s1,x), lower = 0.001, upper = 10, method = "L-BFGS-B")
fkl <- optim(.5, \(x) kl_grav(mmb$meaps,s1,x, furness=TRUE), lower = 0.001, upper = 10, method = "L-BFGS-B")
fkl2 <- optim(.5, \(x) kl_grav(mmb2$meaps,s2,x, furness=TRUE), lower = 0.001, upper = 10, method = "L-BFGS-B")
fref <- flux_grav_furness(s1, fkl$par)
f2 <- flux_grav_furness(s2, fkl$par)
f22 <- flux_grav_furness(s2, fkl2$par)
rmeaps::r2kl2(c(emp_flux(s1, fref)$s), c(emp_flux(s1, mmb$meaps)$s))
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

save(fluxg, fluxg2, fluxg22, fkl, fkl2, file = "output/flux_grav.srda")

# matrice de dérivées -----------------------------
sr <- s1
n <- nrow(sr$habs)
k <- nrow(sr$emps)
gh <- unique(sr$habs$g) |> sort()
ge <- unique(sr$emps$g) |> sort()
groups <- expand_grid(ge,gh) |> 
  mutate(l = str_c("h", gh,"-e", ge))
groups <- set_names(purrr::transpose(list(groups$gh, groups$ge)), groups$l)
# groups <- list(h1g1 = list(1,1), h2g2 = list(2,2))
names(groups) <- map_chr(groups, ~str_c("e", .x[[1]], "-h", .x[[2]]))
ms_odd <- list_modify( 
  list(ref = matrix(1, nrow = n, ncol = k)),
  !!!imap(groups,  ~{
    m <- matrix(1, nrow = n, ncol = k) 
    m[which(sr$habs$g==.x[[1]]), which(sr$emp$g==.x[[2]])] <- 2
    m
  }))

shufs <- do.call(rbind, map(1:128, ~sample.int(n,n)))

plan("multisession", workers=2)

res <- future_imap(
  ms_odd, ~{
    mm <- rmeaps::meaps_multishuf(sr$rk,
                      emplois = rep(1, k),
                      actifs = rep(1, n),
                      f = sr$f,
                      modds = .x ,
                      shuf = shufs,
                      nthreads = 4, 
                      progress=FALSE)
    flux <- emp_flux(sr, mm)
    ed <- mm * sr$dist
    emp_i <- matrixStats::rowSums2(mm)
    d_ind <- matrixStats::rowSums2(ed)/emp_i
    list(habs = bind_cols(sr$habs, tibble(d = d_ind), scn = .y),
         flux = flux$s,
         flux_ec = flux$ec,
         emps = mm,
         empec = 0)
  }, .options = furrr_options(seed=TRUE), .progress = TRUE) |> purrr::transpose()
tres <- bind_rows(res$habs)

mmflux <- map(res$flux[-1], ~.x - res$flux$ref)
save(mmflux, file="output/mmflux.rda")
ee <- map(mmflux, as.vector) |> flatten_dbl() |> matrix(ncol=9, nrow=9) |> eigen() 
save(ee, file="output/eigenvalues.srda")
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
  # ofce::table_ofce() |> 
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

save(flux3x3, file="output/flux3x3.rda" |> glue::glue())  

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

