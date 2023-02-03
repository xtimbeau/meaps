#============================================================#
# MEAPS avec file d'attente et hétérogénité absorption/fuite #
# génération des graphiques et tableaux pour le .qmd         #
#============================================================#

# on utilise maintenant meaps (avec pl et odd ratio)
# le c++ est dans meaps_rcpp

# init -------------------------
library(tidyverse)
library(conflicted)
library(furrr)
source("radiation/radiation functions.r")
Rcpp::sourceCpp("radiation/meaps_rcpp.cpp", echo=FALSE, showOutput = FALSE)
future::plan("multisession", workers = 8)
library(tictoc)
library(ggnewscale)
library(matrixStats)
library(patchwork)
library(ofce)
library(Rcpp)
options(ofce.background_color = "grey97")
showtext::showtext_opts(dpi = 200)
showtext::showtext_auto()
conflict_prefer_all("dplyr", quiet = TRUE)
sysfonts::font_add_google('Nunito')
options(ofce.base_family = "Nunito")
options(ofce.base_size = 9)
plan("multisession", workers = 8)
Rcpp::sourceCpp("radiation/meaps2.cpp", echo=FALSE, showOutput = FALSE)
source("radiation/meaps2.r")
n <- 5000
k <- 5000
bins <- 1.2/0.05

## génération --------------
# on représente ici une agglomération centrale, plus des villages avec des emplois
# mais pas en nombre suffisant dans les villages
set.seed(1942) 
habc <- cbind(pos_cnorm(n=70*n/100, sigma = 0.1, centre = c(0.8, 1)), f=0.1, g = 1)
habv1 <- cbind(pos_cnorm(n=15*n/100, sigma = 0.1, centre = c(1.2, 0.4)), f=0.1, g = 3)
habv12 <- cbind(pos_cnorm(n=15*n/100, sigma = 0.1, centre = c(1.6, 0.1)), f=0.1, g = 3)
habv2 <- cbind(pos_cnorm(n=15*n/100, sigma = 0.1, centre = c(0.2, 1.2)), f=0.1, g = 2)
hab <- rbind(habc, habv2, habv1)
hab2 <- rbind(habc, habv2, habv12)
empc <- cbind(pos_cnorm(n=80/100*k, sigma = 0.05, centre = c(0.8, 1)), p=1, g=1)
empv1 <- cbind(pos_cnorm(n=5/100*k, sigma = 0.05, centre = c(1.2, 0.4)), p=1, g=3)
empv12 <- cbind(pos_cnorm(n=5/100*k, sigma = 0.05, centre = c(1.6, 0.1)), p=1, g=3)
empv2 <- cbind(pos_cnorm(n=5/100*k, sigma = 0.05, centre = c(0.2, 1.2)), p=1, g=2)
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
save(dds, dds2, file="meaps-doc/datadoc/dds.rda")
meann <- function(n) function(x) ifelse(length(x)>n, mean(x), NA)
### cartes ----------------
(gcarte_ss <- 
   (ggplot()+
      stat_binhex(data = as_tibble(s1$hab),
                  aes(x=x, y=y, fill=100*after_stat(density)), binwidth=0.05)+
      scale_fill_distiller(palette="Greens", direction=1, name = "densité\nd'habitants")+
      coord_equal(xlim=c(0,1.5), ylim=c(0,1.5))+
      geom_text(data = s1$hgroupes, aes(x=x, y=y, label = g_label), nudge_y = -0.2,  size = 2) +
      labs(title = "Habitants")+
      theme_void(base_size = 8)+ 
      theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
            plot.margin = margin(6,6,6,6),
            panel.background = element_rect(fill="grey97")))+
   (ggplot()+
      stat_binhex(data=as_tibble(s1$emp),
                  aes(x=x,y=y, fill=100*after_stat(density)), binwidth=0.05)+
      scale_fill_distiller(palette = "Oranges", direction=1, name = "densité\nd'emplois")+
      coord_equal(xlim=c(0,1.5), ylim=c(0,1.5))+
      geom_text(data = s1$egroupes, aes(x=x, y=y, label = g_label), size = 2, nudge_y = -0.2) +
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
       coord_equal(xlim=c(0,1.5), ylim=c(0,1.5))+
       geom_text(data = s2$hgroupes, aes(x=x, y=y, label = g_label), nudge_y = -0.2,  size = 2) +
       labs(title = "Habitants")+
       theme_void(base_size = 8)+ 
       theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
             plot.margin = margin(6,6,6,6),
             panel.background = element_rect(fill="grey97")))+
    (ggplot()+
       stat_binhex(data=as_tibble(s2$emp),
                   aes(x=x,y=y, fill=100*after_stat(density)), binwidth=0.05)+
       scale_fill_distiller(palette = "Oranges", direction=1, name = "densité\nd'emplois")+
       coord_equal(xlim=c(0,1.5), ylim=c(0,1.5))+
       geom_text(data = s2$egroupes, aes(x=x, y=y, label = g_label), size = 2, nudge_y = -0.2) +
       labs(title = "Emplois")+
       theme_void(base_size = 8)+ 
       theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
             plot.margin = margin(6,6,6,6),
             panel.background = element_rect(fill="grey97"))) + 
    plot_layout(guides = 'collect'))
graph2png(gcarte_ss, rep="meaps-doc/svg", ratio = 2)
graph2png(gcarte_ss2, rep="meaps-doc/svg", ratio = 2)
save(gcarte_ss, file = "meaps-doc/graphs/gcarte_ss.rda")
save(gcarte_ss2, file = "meaps-doc/graphs/gcarte_ss2.rda")

# calculs ----------------------------

bench::mark(v1 = meaps_cpp(s1$rk,f = s1$f, p = s1$p, shuf = 1:k))
bench::mark(v2 = meaps_rcpp(s1$rk,emplois=rep(1, n), actifs = rep(1, k), odds = s1$p, f = s1$f, shuf = 1:k))

mm <- rmeaps(emp = emp, hab = hab, meaps_ver = 2)
mm2 <- rmeaps(emp = emp2, hab = hab2, meaps_ver = 2)
# save(mm, mm2, file = "radiation/graphs/mms.rda")

## matrice de flux ----------------
flux <- mm$meaps |>
  as_tibble() |>
  mutate(hab = 1:nrow(mm$meaps), 
         gh = factor(!!hab[,"g"])) |> 
  pivot_longer(cols = starts_with("V"), names_to = "emp", values_to = "pemp") |> 
  mutate(emp = str_sub(emp, 2,-1) |> as.integer(),
         ge = (!!emp[,"g"])[emp],
         gecroixgh = str_c(ge,"x",gh)) |> 
  group_by(gecroixgh) |> 
  summarize(pemp = sum(pemp), 
            nhab = n_distinct(hab),
            nemp = n_distinct(emp)) |> 
  mutate(ge = str_c("e", map_chr(str_split(gecroixgh, "x"), 1)),
         gh = str_c("h", map_chr(str_split(gecroixgh, "x"),2))) |> 
  pivot_wider(id_cols = gh, names_from = ge, values_from = pemp) |>
  add_total() |> 
  rowwise() |> 
  mutate(total = sum(c_across(2:4))) |> 
  ungroup() |> 
  mutate(across(2:5, ~prettyNum(round(.x), format ="d", big.mark = " ")))
### flux2 ---------------
flux2 <- mm2$meaps |>
  as_tibble() |>
  mutate(hab = 1:nrow(mm2$meaps), 
         gh = factor(!!hab2[,"g"])) |> 
  pivot_longer(cols = starts_with("V"), names_to = "emp", values_to = "pemp") |> 
  mutate(emp = str_sub(emp, 2,-1) |> as.integer(),
         ge = (!!emp2[,"g"])[emp],
         gecroixgh = str_c(ge,"x",gh)) |> 
  group_by(gecroixgh) |> 
  summarize(pemp = sum(pemp), 
            nhab = n_distinct(hab),
            nemp = n_distinct(emp)) |> 
  mutate(ge = str_c("e", map_chr(str_split(gecroixgh, "x"), 1)),
         gh = str_c("h", map_chr(str_split(gecroixgh, "x"),2))) |> 
  pivot_wider(id_cols = gh, names_from = ge, values_from = pemp) |>
  add_total() |> 
  rowwise() |> 
  mutate(total = sum(c_across(2:4))) |> 
  ungroup() |> 
  mutate(across(2:5, ~prettyNum(round(.x), format ="d", big.mark = " ")))
knitr::kable(flux)
knitr::kable(flux2)

save(flux, flux2, file = "meaps-doc/datadoc/tblflux.rda")

(gdenshabg <- ggplot(mm$hab)+
    geom_density(aes(x=d, group=g, fill=factor(g), col=factor(g)), alpha = 0.5)+
    geom_density(data = mm2$hab, aes(x=d, group=g, col=factor(g)), alpha = 0.5, linetype ="dashed")+
    scale_color_discrete(name="pôle d'habitation",
                         aesthetics = c('color','fill'), 
                         labels = c("h1", "h2", "h3"))+
    xlim(c(0,1.5))+ylim(c(0,10))+ylab(NULL)+xlab("distances")+
    theme_ofce(base_size=9, base_family = "Nunito", legend.position = "right"))
graph2png(gdenshabg, rep="meaps-doc/svg")
save(gdenshabg, file = "meaps-doc/graphs/gdenshabg.rda")


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
  coord_equal(xlim=c(0,1.5), ylim=c(0,1.5))+
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
  xlim(c(0,1.5))+ylim(c(0,8))+ylab(NULL)+xlab(NULL)+
  theme_minimal()+
  theme(text=element_text(size=6))
gdemp <- ggplot()+
  stat_summary_hex(data = mm$emp, aes(x=x, y=y, z=d), fun = meann(5), bins=25)+
  guides(fill=guide_colourbar("Distance\npour un emploi"))+
  scale_fill_distiller(palette="Oranges", direction=1)+
  xlab(NULL)+ylab(NULL)+labs(title="Distance moyenne parcourue pour un emploi") +
  geom_text(data = s1$egroupes, aes(x=x, y=y, label = g_label), size = 2, nudge_y = -0.2) +
  coord_equal(xlim=c(0,1.5), ylim=c(0,1.5))+
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
  xlim(c(0,1.5))+ylim(c(0,15))+ylab(NULL)+xlab(NULL)+
  theme_minimal()+
  theme(text=element_text(size=6))

(gdistances <- 
    (gdhab+inset_element(gdenshab, 0., 0., 0.3, 0.3))+
    (gdemp+inset_element(gdensemp, 0., 0., 0.3, 0.3))+
    plot_layout(guides = "collect"))

graph2png(gdistances, rep="radiation/svg", ratio = 2)
save(gdistances, file = "meaps-doc/graphs/gdistances.rda")

### s2 ---------------------
gdhab <- ggplot()+
  stat_summary_hex(data = mm2$hab, aes(x=x, y=y, z=d), fun = meann(5), bins=25)+
  guides(fill=guide_colourbar("Distance\npar habitant"))+
  scale_fill_distiller(palette="Greens", direction=1)+
  xlab(NULL)+ylab(NULL)+labs(title="Distance moyenne parcourue par habitant")+
  geom_text(data = s2$hgroupes, aes(x=x, y=y, label = g_label), nudge_y = -0.2,  size = 2)+
  coord_equal(xlim=c(0,1.5), ylim=c(0,1.5))+
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
  xlim(c(0,1.5))+ylim(c(0,8))+ylab(NULL)+xlab(NULL)+
  theme_minimal()+
  theme(text=element_text(size=6))
gdemp <- ggplot()+
  stat_summary_hex(data = mm2$emp, aes(x=x, y=y, z=d), fun = meann(5), bins=25)+
  guides(fill=guide_colourbar("Distance\npour un emploi"))+
  scale_fill_distiller(palette="Oranges", direction=1)+
  xlab(NULL)+ylab(NULL)+labs(title="Distance moyenne parcourue pour un emploi") +
  geom_text(data = s2$egroupes, aes(x=x, y=y, label = g_label), size = 2, nudge_y = -0.2) +
  coord_equal(xlim=c(0,1.5), ylim=c(0,1.5))+
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
  xlim(c(0,1.5))+ylim(c(0,15))+ylab(NULL)+xlab(NULL)+
  theme_minimal()+
  theme(text=element_text(size=6))

(gdistances2 <- 
    (gdhab+inset_element(gdenshab, 0., 0., 0.3, 0.3))+
    (gdemp+inset_element(gdensemp, 0., 0., 0.3, 0.3))+
    plot_layout(guides = "collect"))

graph2png(gdistances2, rep="meaps-doc/svg", ratio = 2)
save(gdistances2, file = "meaps-doc/graphs/gdistances2.rda")

## ergodicité ------------------------
# attention on utilise ici la version 1 de meaps
# ce qu'il y a en moins : pas de traitement des paquets, pas de traitement des odds
# pas de débordement
# mais dans ce cas ce n'est pas grave

plan("multisession", workers = 7)

tic();
libres.raw <- future_imap_dfr(1:500, ~{
  Rcpp::sourceCpp("radiation/meaps.cpp")
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

qs::qsave(libres.raw, "v2/libres.raw.qs")

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
    xlab("Nombre de tirages")+ylab("part d'emploi affectée au carreau d'emploi")+
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill="black"),
          plot.background = element_rect(fill="white"),
          text = element_text(color = "black"),
          axis.text = element_text(color = "black"),
          strip.text = element_text(color = "black")))
graph2png(gemploi_erg, rep="radiation/svg")
save(gemploi_erg, file="radiation/graphs/gemploi_erg.rda")
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
# matrice de flux
fluxs <- erg_libre |>
  group_by(ge, gh, draw) |> 
  summarize(flux = sum(ps)) |> 
  mutate( cat = str_c("h",gh,"->","e", ge))
gfluxs <- ggplot(fluxs)+
  geom_boxplot(aes(y=flux, x=cat, color=gh),
               show.legend = FALSE)+
  scale_y_log10()+
  xlab(NULL)+ylab("Nombre de trajets")+
  theme_ofce()
graph2png(gfluxs, rep="meaps-doc/svg")
save(fluxs, file="meaps-doc/graphs/fluxs.rda")

# test 5 rang n
rangns <- erg_libre |> group_by(ehex, draw) |> summarize(across(c(rm, rsd, ge), first))
ggplot(rangns) + 
  geom_line(aes(x=draw, y=rm, group = ehex), col="white", alpha=0.3, size=0.2) +
  theme_minimal()+theme(panel.grid = element_blank(),
                        panel.background = element_rect(fill="black"),
                        plot.background = element_rect(fill="white"),
                        text = element_text(color = "black"),
                        axis.text = element_text(color = "black"),
                        strip.text = element_text(color = "black"))+
  xlab("Nombre de tirages cumulés")+ylab("rang moyen au moment du remplissage")
ggplot(rangns) + 
  geom_line(aes(x=draw, y=rsd, group = ehex), col="white", alpha=0.25, size=0.2) +
  theme_minimal()+theme(panel.grid = element_blank(),
                        panel.background = element_rect(fill="black"),
                        plot.background = element_rect(fill="white"),
                        text = element_text(color = "black"),
                        axis.text = element_text(color = "black"),
                        strip.text = element_text(color = "black"))+
  xlab("Nombre de tirages cumulés")+ylab("écart type du rang moyen au moment du remplissage")
gnmsd <- ggplot(rangns |> filter(draw==500))+
  geom_point(aes(x=rm, y=rsd, color=ge)) +
  xlab("rang moyen au moment du remplissage")+
  ylab("écart type du rang moyen au moment du remplissage") + 
  scale_color_discrete(name = "Pôle d'emploi")+
  theme_ofce()+theme(legend.position="right")
graph2png(gnmsd, rep="meaps-doc/svg")
save(gnmsd, file="meaps-doc/graphs/gmsd_erg.rda")


## modification de la pa --------------------
# on le fait pour quelques tirages
sr <- s1
empg1 <- emp
empg1[emp[,"g"]==1, "p"] <- 0.5
empg1.1 <- emp
empg1.1[emp[,"g"]==1, "p"] <-  2
sg1 <- make_tibs(empg1, hab)
sg1.1 <- make_tibs(empg1.1, hab)

empg2 <- emp
empg2[emp[,"g"]==2, "p"] <- 0.5
sg2 <- make_tibs(empg2, hab)

empg2.1 <- emp
empg2.1[emp[,"g"]==2, "p"] <- 2
sg2.1 <- make_tibs(empg2.1, hab)

empg3 <- emp
empg3[emp[,"g"]==3, "p"] <- 0.5
sg3 <- make_tibs(empg3, hab)

empg3.1 <- emp
empg3.1[emp[,"g"]==3, "p"] <- 2
sg3.1 <- make_tibs(empg3.1, hab)

shufs <- do.call(rbind, map(1:100, ~sample.int(n,n)))
scns <- set_names(list(sr, sg1, sg1.1, sg2, sg2.1, sg3, sg3.1), 
                  c("sr", "sg1", "sg1.1", "sg2", "sg2.1", "sg3", "sg3.1"))
res.raw <- imap(scns, ~rmeaps_bstp(.x, shufs, workers = 8))

res <- imap_dfr(
  scns,
  ~{
    mm <- rmeaps_bstp(.x, shufs, workers = 8)
    emps <- mm$emps
    ed <- emps * .x$dist
    emp_i <- matrixStats::rowSums2(emps)
    d_ind <- matrixStats::rowSums2(ed)/emp_i
    bind_cols(.x$habs, tibble(d = d_ind), scn = .y)
  })
res <- res |> mutate(scn = factor(scn, c("sr", "sg1", "sg1.1", "sg2", "sg3")))

ggplot(data = res |> filter(scn %in% c("sr", "sg1", "sg1.1")))+
  geom_density(aes(x=d, group=interaction(scn, g), col=factor(g), fill=factor(g), linetype = scn), 
               alpha = 0.25)+
  scale_color_discrete(name="pôle d'habitation", 
                       aesthetics = c('color','fill'),
                       labels = c("h1", "h2", "h3"))+
  xlim(c(0,1.5))+ylim(c(0,10))+ylab(NULL)+xlab("distances")+
  theme_ofce(base_size=9)+theme(legend.position = "right")

ggplot(data = res |> filter(scn %in% c("sr", "sg2")))+
  geom_density(aes(x=d, group=interaction(scn, g), col=factor(g), fill=factor(g), linetype = scn),
               alpha = 0.5)+
  scale_color_discrete(name="pôle d'habitation", 
                       aesthetics = c('color','fill'),
                       labels = c("h1", "h2", "h3"))+
  xlim(c(0,1.5))+ylim(c(0,10))+ylab(NULL)+xlab("distances")+
  theme_ofce(base_size=9)+theme(legend.position = "right")

# res <- purrr::transpose(res)
# resm <- map(res, ~reduce(purrr::transpose(.x)$emps, `+`)/length(.x))
# habm <- map(res, ~reduce(purrr::transpose(.x)$hab, 
#                          function(t,t2) bind_rows(t,t2,.id="draw")) |>
#               group_by(hab) |> 
#               summarize(
#                 across(c(hhex,x,y,g,f), first),
#                 dm = mean(d),
#                 dsd = sd(d)))
# 
# flux_pa <- resm$pa |>
#   as_tibble() |>
#   mutate(hab = 1:n, 
#          gh = !!s3$hab$g) |> 
#   pivot_longer(cols = starts_with("V"), names_to = "emp", values_to = "pemp") |> 
#   mutate(emp = str_sub(emp, 2,-1) |> as.integer(),
#          ge = (!!s3$emp$g)[emp],
#          gecroixgh = str_c(ge,"x",gh)) |> 
#   group_by(gecroixgh) |> 
#   summarize(pemp = sum(pemp), 
#             ge = first(ge), 
#             gh = first(gh),
#             nhab = n_distinct(hab),
#             n_emp = n_distinct(emp)) |> 
#   mutate(ge = str_c("e", ge),
#          gh = str_c("h", gh)) |> 
#   pivot_wider(id_cols = gh, names_from = ge, values_from = pemp) |>
#   add_total() |> 
#   rowwise() |> 
#   mutate(total = sum(c_across(2:4))) |> 
#   ungroup() |> 
#   mutate(across(2:5, ~prettyNum(round(.x), format ="d", big.mark = " ")))
# flux_npa <- resm$npa |>
#   as_tibble() |>
#   mutate(hab = 1:n, 
#          gh = !!s1$hab$g) |> 
#   pivot_longer(cols = starts_with("V"), names_to = "emp", values_to = "pemp") |> 
#   mutate(emp = str_sub(emp, 2,-1) |> as.integer(),
#          ge = (!!s1$emp$g)[emp],
#          gecroixgh = str_c(ge,"x",gh)) |> 
#   group_by(gecroixgh) |> 
#   summarize(pemp = sum(pemp), 
#             ge = first(ge), 
#             gh = first(gh),
#             nhab = n_distinct(hab),
#             n_emp = n_distinct(emp)) |> 
#   mutate(ge = str_c("e", ge),
#          gh = str_c("h", gh)) |> 
#   pivot_wider(id_cols = gh, names_from = ge, values_from = pemp) |>
#   add_total() |> 
#   rowwise() |> 
#   mutate(total = sum(c_across(2:4))) |> 
#   ungroup() |> 
#   mutate(across(2:5, ~prettyNum(round(.x), format ="d", big.mark = " ")))
# 
# knitr::kable(flux_npa)
# knitr::kable(flux_pa)
# 
# gdenshabg3 <- ggplot(mapping=aes(x=dm, group=g, col=factor(g)))+
#   geom_density(data = habm$sr, mapping = aes(fill=factor(g), alpha = 0.5))+
#   geom_density(data = habm$g1, alpha = 0.5, linetype ="dashed")+
#   geom_density(data = habm$g1.1, alpha = 0.5, linetype ="dotted")+
#   scale_color_discrete(name="pôle d'habitation", 
#                        aesthetics = c('color','fill'),
#                        labels = c("h1", "h2", "h3"))+
#   xlim(c(0,1.5))+ylim(c(0,10))+ylab(NULL)+xlab("distances")+
#   theme_ofce(base_size=9)+theme(legend.position = "right")
# 
# gdenshabg3.1 <- ggplot(mapping=aes(x=dm, group=g, col=factor(g)))+
#   geom_density(data = habm$sr, mapping = aes(fill=factor(g), alpha = 0.5))+
#   geom_density(data = habm$g2, alpha = 0.5, linetype ="dashed")+
#   geom_density(data = habm$g3, alpha = 0.5, linetype ="dotted")+
#   scale_color_discrete(name="pôle d'habitation", 
#                        aesthetics = c('color','fill'),
#                        labels = c("h1", "h2", "h3"))+
#   xlim(c(0,1.5))+ylim(c(0,10))+ylab(NULL)+xlab("distances")+
#   theme_ofce(base_size=9)+theme(legend.position = "right")
# 
# graph2png(gdenshabg3, rep="meaps-doc/svg")
# save(gdenshabg3, file = "meaps-doc/graphs/gdenshabg2.rda")

bnd500 <- dplyr::count(erg_libre |> filter(draw==500)) |> pull(n)
ggplot(erg_libre |> filter(draw == 500)) + 
  geom_histogram(aes(x=lsd, y = log10(after_stat(count))), bins =1000)+ theme_ofce()

ggplot(erg_libre |> filter(draw < 50)) + 
  stat_ecdf(aes(x=lsd, group = factor(draw), color=draw), alpha = 0.2) + theme_ofce()

ggplot(erg_libre |> filter(draw < 50)) + 
  stat_ecdf(aes(x=psd, group = factor(draw), color=draw), alpha = 0.2) + theme_ofce()

fdata <- matrixStats::colSds(coco)/  matrixStats::colMeans2(coco)

bined_coco <- as_tibble(coco[-1,]) |> mutate(hex = hexs) |> pivot_longer(cols = starts_with("V"))
bined_coco |> group_by(hex, name) |> summarize(f = mean(value), sd = sd(value))


## modification de la pa config 2 --------------------
#
n <- 5000
hab <- cbind(pos_cnorm(n=n, sigma = 0.2, centre = c(0.5, 0.5)), f=0.1, g = 1)
hab[hab[,1]<0.5&hab[,2]<0.5, "g"] <- 1
hab[hab[,1]>=0.5&hab[,2]<0.5, "g"] <- 2
hab[hab[,1]>=0.5&hab[,2]>=0.5, "g"] <- 3
hab[hab[,1]<0.5&hab[,2]>=0.5, "g"] <- 4

emp <-do.call(
  rbind,
  imap(
    list(c(.2, .2), c(.2, .8), c(.8, .2), c(.8,.8)),
    ~cbind(pos_cnorm(n=22.5/100*n, sigma = 0.05, centre = .x), p=1, g=.y)))
sr <- make_tibs(emp, hab)
sgi <- imap(1:4, ~{
  empg <- emp
  empg[emp[,"g"]==.x, "p"] <- 0.5
  make_tibs(empg, hab)})
names(sgi) <- str_c("sg", 1:4)

shufs <- do.call(rbind, map(1:100, ~sample.int(n,n)))
scns <- list_modify(list(sr = sr), !!!sgi)
res <- imap(
  scns,
  ~{
    mm <- rmeaps_bstp(.x, shufs, workers = 6)
    flux <- emp_flux(.x, mm$emps, mm$emps2)
    ed <- mm$emps * .x$dist
    emp_i <- matrixStats::rowSums2(mm$emps)
    d_ind <- matrixStats::rowSums2(ed)/emp_i
    list(habs = bind_cols(.x$habs, tibble(d = d_ind), scn = .y),
         flux = flux$s,
         flux_ec = flux$ec,
         emps = mm$emps,
         empec = sqrt(mm$emps2- mm$emps^2))
  }) |> purrr::transpose()
tres <- bind_rows(res$habs)

map(res$flux[-1], ~round(.x-res$flux$sr))

(gcarte_cfg2 <- 
    (ggplot()+
       stat_binhex(data = as_tibble(sr$hab),
                   aes(x=x,y=y, fill=100*after_stat(density)), binwidth=0.05)+
       scale_fill_distiller(palette="Greens", direction=1, name = "densité\nd'habitants")+
       coord_equal(xlim=c(0,1), ylim=c(0,1))+
       geom_text(data = sr$hgroupes, aes(x=x, y=y, label = g_label), nudge_y = -0.2,  size = 2) +
       labs(title = "Habitants")+
       theme_void(base_size = 8)+ 
       theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
             plot.margin = margin(6,6,6,6),
             panel.background = element_rect(fill="grey97")))+
    (ggplot()+
       stat_binhex(data=as_tibble(sr$emp),
                   aes(x=x,y=y, fill=100*after_stat(density)), binwidth=0.05)+
       scale_fill_distiller(palette = "Oranges", direction=1, name = "densité\nd'emplois")+
       coord_equal(xlim=c(0,1), ylim=c(0,1))+
       geom_text(data = sr$egroupes, aes(x=x, y=y, label = g_label), size = 2, nudge_y = -0.2) +
       labs(title = "Emplois")+
       theme_void(base_size = 8)+ 
       theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
             plot.margin = margin(6,6,6,6),
             panel.background = element_rect(fill="grey97"))) + 
    plot_layout(guides = 'collect'))

ggplot(data = tres |> filter(scn!="sr") )+
  geom_density(aes(x=d, col = scn, fill=scn), alpha = 0.3, position="identity")+
  geom_density(data=tres |> rename(snc=scn) |> filter(snc=="sr"), 
               aes(x=d),
               col = "gray",
               fill="gray",
               alpha = 0.3,
               position="identity")+
  facet_wrap(vars(scn))+
  ylab(NULL)+xlab("distances")+
  theme_ofce(base_size=9)+theme(legend.position = "right")

## odds en matrice ---------------

sourceCpp("v2/annexes/meaps_oddmatrix.cpp")

# n <- 5000
# hab <- cbind(pos_cnorm(n=n, sigma = 0.15, centre = c(0.5, 0.5)), f=0.1, g = 1)
# hab[hab[,1]<0.5&hab[,2]<0.5, "g"] <- 1
# hab[hab[,1]>=0.5&hab[,2]<0.5, "g"] <- 2
# hab[hab[,1]>=0.5&hab[,2]>=0.5, "g"] <- 3
# hab[hab[,1]<0.5&hab[,2]>=0.5, "g"] <- 4
# 
# emp <-do.call(
#   rbind,
#   imap(
#     list(c(.2, .2), c(.2, .8), c(.8, .2), c(.8,.8)),
#     ~cbind(pos_cnorm(n=22.5/100*n, sigma = 0.05, centre = .x), p=1, g=.y)))


sr <- s1
gh <- unique(s1$habs$g) |> sort()
ge <- unique(s1$emps$g) |> sort()
ms_odd <- list_modify( 
  list(ref = matrix(1, nrow = n, ncol = n)),
  !!!imap(set_names(ge, str_c("e",ge, "-h1")),  ~{
    m <- matrix(1, nrow = n, ncol = n) 
    m[which(hab[,"g"]==1), which(emp[, "g"]==.x)] <- 2
    m
  }),
  !!!imap(set_names(ge, str_c("toush", "-e", ge)),  ~{
    m <- matrix(1, nrow = n, ncol = n) 
    m[,which(emp[,"g"]==.x)] <- 2
    m
  }))

shufs <- do.call(rbind, map(1:10, ~sample.int(n,n)))
plan("multisession", workers =6)
res <- furrr::future_imap(
  ms_odd, ~{
    sourceCpp("v2/annexes/meaps_oddmatrix.cpp")
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
  }, .options = furrr_options(seed=TRUE)) |> purrr::transpose()
tres <- bind_rows(res$habs)

map(map(res$flux[-1], ~.x-res$flux$ref), ~{
  etp1 <- cbind(.x, rowSums2(.x))
  rbind(etp1, colSums2(etp1)) |> round(1)
}) 
  

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

ggplot(data = tres |> filter(scn!="sr") )+
  geom_density(aes(x=d, col = scn, fill=scn), alpha = 0.3, position="identity")+
  geom_density(data=tres |> rename(snc=scn) |> filter(snc=="sr"), 
               aes(x=d),
               col = "gray",
               fill="gray",
               alpha = 0.3,
               position="identity")+
  facet_wrap(vars(scn))+
  ylab(NULL)+xlab("distances")+
  theme_ofce(base_size=9)+theme(legend.position = "right")

# Anciens ------------------------------------------
## cas uniforme en cercle ------------------

hab <- cbind(randompos(n=k, rayon = 0.2), f=0.1)
emp <- do.call(rbind, list(
  cbind(randompos(n=floor(n/4), rayon = 0.1, centre = c(0.2,0.2)), p=1),
  cbind(randompos(n=floor(n/4), rayon = 0.1, centre = c(0.2,0.8)), p=1),
  cbind(randompos(n=floor(n/4), rayon = 0.1, centre = c(0.8,0.2)), p=1),
  cbind(randompos(n=floor(n/4), rayon = 0.1, centre = c(0.8,0.8)), p=1)))
rownames(emp) <- 1:nrow(emp)
dist <- rdist::cdist(hab[,1:2], emp[,1:2])
dimnames(dist) <- list(rownames(hab), rownames(emp))
rkdist <- rowRanks(dist)
dimnames(rkdist) <- list(rownames(hab), rownames(emp))

# on accroit les proba dans un coin

emp2 <- emp
emp2[1:floor(n/4),"p"] <- 2

tic();res <- rmeaps(emp, hab); toc()
tic();res2 <- rmeaps(emp2, hab); toc()

tibhab <- as_tibble(hab) |> mutate(d_ind = res$d_ind) 
tibemp <- as_tibble(emp) |> mutate(file = res$file, efi = res$emplois, d_emp = res$d_emp)

(g1 <- ggplot()+
    stat_summary_hex(data = tibhab,
                     aes(x=x,y=y, z = d_ind), fun = mean, bins=50)+
    scale_fill_distiller(palette="Greens", direction=1)+
    new_scale_fill()+
    stat_summary_hex(data=tibemp,
                     aes(x=x,y=y, z=d_emp), alpha = 1, bins=25)+
    scale_fill_distiller(palette = "Oranges", direction=1)+
    coord_equal()
)

tibhab2 <- bind_rows(
  tibhab |> mutate(scenario = "uniforme"),
  tibhab |> mutate(scenario = "bas gauche", d_ind = res2$d_ind))
tibemp2 <- bind_rows(
  tibemp |> mutate(scenario = "uniforme"),
  tibemp |> mutate(scenario = "bas gauche", d_emp = res2$d_emp, file = res2$file, efi = res2$emplois))

(g3 <- ggplot()+
    stat_summary_hex(data = tibhab2,
                     aes(x=x,y=y, z = d_ind), fun = mean, bins=50)+
    scale_fill_distiller(palette="Greens", direction=1)+
    new_scale_fill()+
    stat_summary_hex(data=tibemp2,
                     aes(x=x,y=y, z=d_emp), fun = mean, alpha = 1, bins=25)+
    scale_fill_distiller(palette = "Oranges", direction=1)+
    coord_equal()+
    facet_wrap(vars(scenario)))

ggplot(tibhab2) +
  geom_histogram(aes(x=d_ind, 
                     fill = scenario),
                 position = "identity",
                 alpha=0.5, 
                 bins= 50)

data <- meaps_btsp(emp, hab)

(g2 <- ggplot()+
    stat_summary_hex(data = bind_rows(data$habs,
                                      data$all |> filter(h, mc==1)),
                     aes(x=x,y=y, z = d), fun = mean, bins=50)+
    scale_fill_distiller(palette="Greens", direction=1)+
    new_scale_fill()+
    stat_summary_hex(data=bind_rows( data$emps,
                                     data$all |> filter(!h, mc==1)),
                     aes(x=x,y=y, z=d), fun = mean, alpha = 1, bins=25)+
    scale_fill_distiller(palette = "Oranges", direction=1)+
    coord_equal()+
    facet_wrap(vars(mc)))

dist_tt <- as_tibble(t(dist_mc)) |>
  pivot_longer(cols = starts_with("V")) |>
  rename(tirage = name, distance = value)
ggplot(tibble(d = dist_mc[1,]))+
  geom_histogram(aes(x=d), bins = 100)

dist_tt |> 
  group_by(tirage) |>
  summarize(m = mean(distance), s = sd(distance), min = min(distance), max= max(distance))
ggplot(dist_tt) +
  geom_histogram(aes(x=distance, fill = tirage),
                 bins= 100, 
                 position = "identity", 
                 alpha=0.01, 
                 col=NA, 
                 show.legend = FALSE)

## cas gaussien --------------------

hab <- cbind(pos_cnorm(n=k, sigma = 0.2), f=0.1)
emp <- cbind(pos_cnorm(n=k/2, sigma = 0.05), p=1)
emp2 <- cbind(pos_cnorm(n=k/2, centre = c(0.8, 0.8), sigma = 0.05), p=1)

bench::mark(rmeaps(emp = emp, hab = hab), 
            rmeaps(emp = emp, hab = hab, rcpp = FALSE), check=FALSE )

meann <- function(n) function(x) ifelse(length(x)>n, mean(x), NA)

(g1 <- (ggplot()+
          stat_binhex(data = as_tibble(hab),
                      aes(x=x,y=y), bins=20)+
          scale_fill_distiller(palette="Greens", direction=1)+
          coord_equal(xlim=c(0,1.2), ylim=c(0,1.2)))+
    (ggplot()+
       stat_binhex(data=as_tibble(rbind(emp, emp2)),
                   aes(x=x,y=y), bins=20)+
       scale_fill_distiller(palette = "Oranges", direction=1)+
       coord_equal(xlim=c(0,1.2), ylim=c(0,1.2))))

mm <- rmeaps(emp = emp, hab = hab)
mm2 <- rmeaps(emp = rbind(emp,emp2), hab = hab)
gdhab <- ggplot()+
  stat_summary_hex(data = mm$hab, aes(x=x, y=y, z=d), fun = meann(10), bins=25)+
  guides(fill=guide_colourbar("Distance"))+
  scale_fill_distiller(palette="Greens", direction=1)+
  xlab(NULL)+ylab(NULL)+labs(title="Distance moyenne parcourue par habitant")+
  coord_equal(xlim=c(0,1.2), ylim=c(0,1.2))+ theme_ofce() +
  theme(legend.position="right", 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major.y = element_blank())
gdenshab <- ggplot(mm$hab)+
  geom_density(aes(x=d), 
               fill="green", 
               col=NA,
               alpha = 0.5)+
  xlim(c(0,0.8)) +
  ylab(NULL)+
  xlab(NULL)+
  theme(text=element_text(size=6))+theme_minimal()
gdemp <- ggplot()+
  stat_summary_hex(data = mm$emp,
                   aes(x=x, y=y, z=d),
                   fun = meann(10),
                   bins=25)+
  guides(fill=guide_colourbar("Distance"))+
  scale_fill_distiller(palette="Oranges", direction=1)+
  xlab(NULL)+ylab(NULL)+labs(title="Distance moyenne parcourue par emploi")+
  coord_equal(xlim=c(0,1.2), ylim=c(0,1.2))+ theme_ofce() + 
  theme(legend.position="right", 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major.y = element_blank())
gdensemp <- ggplot(mm$emp)+
  geom_density(aes(x=d),
               fill="Orange", 
               col=NA, alpha = 0.5)+ 
  xlim(c(0,0.8))+ylab(NULL)+xlab(NULL)+
  theme(text=element_text(size=6))+theme_minimal()

(gdistances <- 
    (gdhab+inset_element(gdenshab, 0.7, 0.7, 1, 1))+
    (gdemp+inset_element(gdensemp, 0.7, 0.7, 1, 1)))

gdhab2 <- ggplot()+
  stat_summary_hex(data = mm2$hab, aes(x=x, y=y, z=d), fun = meann(10), bins=25)+
  guides(fill=guide_colourbar("Distance"))+
  scale_fill_distiller(palette="Greens", direction=1)+
  xlab(NULL)+ylab(NULL)+
  coord_equal(xlim=c(0,1.2), ylim=c(0,1.2))+ 
  theme_ofce(legend.position="right")
gdenshab2 <- ggplot(mm2$hab)+
  geom_density(aes(x=d), fill="green", col=NA, alpha = 0.5)+ 
  xlim(c(0,0.8)) +ylab(NULL)+xlab(NULL)+
  theme(text=element_text(size=6))+theme_minimal()
gdemp2 <- ggplot()+
  stat_summary_hex(data = mm2$emp, aes(x=x, y=y, z=d), fun = meann(10), bins=25)+
  guides(fill=guide_colourbar("Distance"))+
  scale_fill_distiller(palette="Oranges", direction=1)+
  xlab(NULL)+ylab(NULL)+
  coord_equal(xlim=c(0,1.2), ylim=c(0,1.2))+
  theme_ofce() + theme(legend.position="right")
gdensemp2 <- ggplot(mm2$emp)+
  geom_density(aes(x=d), 
               fill="Orange", col=NA, alpha = 0.5)+ 
  xlim(c(0,0.8))+ylab(NULL)+xlab(NULL)+
  theme(text=element_text(size=6))+theme_minimal()

(gdistances2 <- 
    (gdhab2+inset_element(gdenshab2, 0.7, 0.7, 1, 1))+
    (gdemp2+inset_element(gdensemp2, 0.7, 0.7, 1, 1)))


(ggplot()+
    stat_summary_hex(data = mm2$hab,
                     aes(x=x,y=y, z = d), fun = mean, bins=20)+
    scale_fill_distiller(palette="Greens", direction=1)+
    coord_equal(xlim=c(0,1.5), ylim=c(0,1.5)))

## ergodicité
library(furrr)
plan("multisession", workers = 16)
dist <- rdist::cdist(hab[, 1:2], emp[,1:2])
rkdist <- matrixStats::rowRanks(dist)
f <- hab[, "f"]
p <- emp[, "p"]

bins <-floor( (range(hab[,"x"])/0.05)[[2]]- (range(hab[,"x"])/0.05)[[1]])

hexhab <- hexbin::hexbin(hab[,"x"], hab[,"y"], xbins=bins, IDs=TRUE)@cID
hexemp <- hexbin::hexbin(emp[,"x"], emp[,"y"], xbins=bins, IDs=TRUE)@cID

libres.raw <- future_imap_dfr(1:500, ~{
  Rcpp::sourceCpp("radiation/meaps.cpp")
  raw <- meaps_cpp(rkdist,f = f, p = p, shuf = sample.int(k,k))
  colnames(raw$emp_meaps) <- str_c("e", 1:ncol(raw$emp_meaps))
  colnames(raw$occup) <- str_c("e", 1:ncol(raw$occup))
  tibe <- as_tibble(raw$emp_meaps)
  names(tibe) <- str_c("e", 1:ncol(raw$emp_meaps))
  tibe <- tibe |> 
    mutate(hab = 1:nrow(tibe)) |> 
    relocate(hab) |> 
    mutate(hhex = hexhab[hab]) |> 
    pivot_longer(cols = starts_with("e"), values_to = 'pemp', names_to = "emp") |> 
    mutate(
      emp = as.numeric(str_sub(emp,2,-1)),
      ehex  = hexemp[emp]) |>
    group_by(hhex, ehex) |>
    summarize(ps = sum(pemp),
              ps2 = sum(pemp^2),
              pn = n(), .groups = "drop")
  tibo <- as_tibble(raw$occup[-1,])
  names(tibo) <- str_c("e", 1:ncol(raw$occup))
  tibo <- tibo |> 
    mutate(hab = 1:nrow(tibo)) |> 
    relocate(hab) |> 
    mutate(hhex = hexhab[hab]) |> 
    pivot_longer(cols = starts_with("e"), values_to = 'libre', names_to = "emp") |> 
    mutate(
      emp = as.numeric(str_sub(emp,2,-1)),
      ehex  = hexemp[emp]) |>
    group_by(hhex, ehex) |>
    summarize(ls = sum(libre),
              ls2 = sum(libre^2),
              ln = n(), .groups = "drop")
  left_join(tibe, tibo, by = c('hhex', 'ehex')) |> 
    mutate(draw = .y) |> 
    relocate(hhex, ehex, draw)
}, .options  =furrr_options(seed=TRUE))

erg_libre <- libres.raw |>
  group_by(hhex, ehex) |> 
  arrange(draw) |> 
  mutate(lm=cumsum(ls)/cumsum(ln), lsd = cumsum(ls2)/cumsum(ln) - 2*lm*cumsum(ls)/cumsum(ln) + lm^2,
         pm=cumsum(ps)/cumsum(pn), psd = cumsum(ps2)/cumsum(pn) - 2*pm*cumsum(ps)/cumsum(pn) + pm^2) |> 
  arrange(desc(psd)) |> 
  ungroup()

# test 1 de tous les hab vers quelques emp
ggplot(erg_libre |> filter(ehex%in%sample(unique(ehex), 9))) + 
  geom_line(aes(x=draw, y=pm, group = hhex), col="grey70", alpha=0.25) +
  facet_wrap(vars(ehex))+
  theme_minimal()

# test 2 vers tous les emp de quelques hab
ggplot(erg_libre |> filter(hhex%in%sample(unique(hhex), 9))) + 
  geom_line(aes(x=draw, y=pm, group = ehex), col="grey70", alpha=0.25) +
  facet_wrap(vars(hhex))+
  theme_minimal()

# test 3 pour quelques emp, état libre de i
ggplot(erg_libre |> filter(ehex%in%sample(unique(ehex), 9))) + 
  geom_line(aes(x=draw, y=lm, group = hhex), col="grey70", alpha=0.25) +
  facet_wrap(vars(ehex))+
  theme_minimal()



fdata <- matrixStats::colSds(coco)/  matrixStats::colMeans2(coco)

bined_coco <- as_tibble(coco[-1,]) |> mutate(hex = hexs) |> pivot_longer(cols = starts_with("V"))
bined_coco |> group_by(hex, name) |> summarize(f = mean(value), sd = sd(value))
