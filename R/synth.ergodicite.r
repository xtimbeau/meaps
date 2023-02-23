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

options(ofce.background_color = "grey99")
showtext::showtext_opts(dpi = 200)
showtext::showtext_auto()
conflict_prefer_all("dplyr", quiet = TRUE)
options(ofce.base_family = "Nunito")
options(ofce.base_size = 9)

source("R/radiation functions.r")

handlers(global = TRUE)
handlers("cli")

n <- 5000
k <- 4500
bins <- 1.2/0.05

load("output/scenarios.rda")

## ergodicité ------------------------

options(future.globals.maxSize=1024^3) # 1Gb
if(fs::file_exists("output/libres.raw.qs")&&FALSE) {
  libres.raw <- qs::qread("output/libres.raw.qs") 
} else {
  plan("multisession", workers = availableCores()-1)
  
  tic();
  libres.raw <- future_imap_dfr(1:500, ~{
    shuf <- sample.int(n,n)
    #raw <- meaps_cpp_v1(s1$rk, f = s1$f, p = s1$p, shuf = shuf)
    raw2 <- rmeaps::meaps_tension_alt(
      rkdist=s1$rk,
      emplois=rep(1,k),
      actifs=rep(1,n),
      f=s1$f,
      shuf = matrix(shuf, nrow=1),
      modds=matrix(1, nrow=n, ncol=k),
      nthreads=1)
    raw <- list(emps = raw2$flux)
    # raw$dispo <- raw2$tension
    colnames(raw$emps) <- str_c("e", 1:ncol(raw$emps))
    # colnames(raw$dispo) <- str_c("e", 1:ncol(raw$dispo))
    tibrn <- tibble(emp = 1:ncol(raw$emps),
                    ehex = s1$emps$ehex[emp],
                    rangn = raw2$tension) |> 
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
    # tibo <- as_tibble(raw$dispo[-1,])
    # names(tibo) <- str_c("e", 1:ncol(raw$dispo))
    # tibo <- tibo |> 
    #   mutate(hab = 1:nrow(tibo)) |> 
    #   relocate(hab) |> 
    #   mutate(hhex = s1$hexhab[hab]) |> 
    #   pivot_longer(cols = starts_with("e"), values_to = 'libre', names_to = "emp") |> 
    #   mutate(
    #     emp = as.numeric(str_sub(emp,2,-1)),
    #     ehex  = s1$hexemp[emp]) |>
    #   group_by(hhex, ehex) |>
    #   summarize(ls = sum(libre),
    #             ls2 = sum(libre^2),
    #             ln = n(), 
    #             .groups = "drop")
    # left_join(tibe, tibo, by = c('hhex', 'ehex')) |> 
    tibe |> 
      left_join(tibrn, by="ehex") |> 
      mutate(draw = .y) |> 
      relocate(hhex, ehex, draw)
  }, .options = furrr_options(seed = TRUE), .progress = TRUE)
  toc();
  
  qs::qsave(libres.raw, "output/libres.raw.qs")
}

libres.raw <- qs::qread("output/libres.raw.qs")

erg_libre <- libres.raw |>
  group_by(hhex, ehex) |> 
  arrange(draw) |> 
  mutate(
    rnc = cumsum(as.numeric(rn)), rm = cumsum(rs)/rnc, rsd = sqrt(cumsum(rs2)/rnc - rm^2),
    # lnc=cumsum(as.numeric(ln)), lm=cumsum(ls)/lnc, lsd = sqrt(cumsum(ls2)/lnc -  lm^2),
    pnc=cumsum(as.numeric(pn)), pm=cumsum(ps)/pnc, psd = sqrt(cumsum(ps2)/pnc -  pm^2)) |> 
  ungroup()

# test 1 de tous les hab vers quelques emp
(gemploi_erg <- ggplot(erg_libre |> filter(ehex%in%sample(unique(ehex), 4), draw<250)) + 
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
graph2png(gemploi_erg, rep="output", ratio=1)
save(gemploi_erg, file="output/gemploi_erg.rda")
# test 2 vers tous les emp de quelques hab
(gemploi_erg <- ggplot(erg_libre |> filter(hhex%in%sample(unique(hhex), 4))) + 
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
          strip.text = element_text(color = "black")))

erg_libre_emp <- erg_libre |> 
  group_by(ehex, hhex) |>
  ungroup() |> 
  filter(ehex%in%sample(unique(ehex), 3))

# test 3 pour quelques emp, état libre de i
# ggplot(erg_libre_emp) + 
#   geom_line(aes(x=draw, y=ls/ln, group = hhex), col="white", alpha=0.05, size=0.05) +
#   geom_line(aes(x=draw, y=lm, group = hhex), col="green", alpha=0.1, size=0.05) +
#   facet_wrap(vars(ehex))+
#   theme_ofce()+
#   theme(panel.grid = element_blank(),
#         panel.background = element_rect(fill="black"),
#         plot.background = element_rect(fill="white"),
#         text = element_text(color = "black"),
#         axis.text = element_text(color = "black"),
#         strip.text = element_text(color = "black"))
# 
# test 4


# test 5 rang n
rangns <- erg_libre |> group_by(ehex, draw) |> summarize(across(c(rm, rsd, ge), first))
rangns_m <- ggplot(rangns |> filter(draw<250)) + 
  geom_line(aes(x=draw, y=rm, group = ehex), col="white", alpha=0.3, size=0.2) +
  theme_ofce()+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill="black", color = NA),
        plot.background = element_rect(fill="white", color = NA),
        text = element_text(color = "black"),
        axis.text = element_text(color = "black"),
        strip.text = element_text(color = "black"))+
  xlab("Nombre de tirages cumulés")+ylab("rang moyen au moment du remplissage")

rangns_ec <- ggplot(rangns |> filter(draw <250)) + 
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

gnmsd <- ggplot(rangns |> filter(draw==max(draw)))+
  geom_point(aes(x=rm, y=rsd, color=ge)) +
  xlab("rang moyen au moment du remplissage")+
  ylab("écart type du rang moyen au moment du remplissage") + 
  scale_color_discrete(name = "Pôle d'emploi")+
  theme_ofce()+theme(legend.position="right")
graph2png(gnmsd, rep="output")
save(gnmsd, file="output/gmsd_erg.rda")
rangns <- s1$ehex |> 
  left_join(rangns |> filter(draw==max(draw)) |> select(ehex, rm, rsd), by="ehex")
carte_erg <- ggplot(rangns)+
  stat_summary_hex(aes(x=x,y=y, z=rm), binwidth=0.1)+
  scale_fill_distiller(palette="Greens", direction=1, name = "densité\nd'habitants")+
  coord_equal(xlim=c(0,2), ylim=c(0,2)) +
  geom_text(data = s1$hgroupes, aes(x=x, y=y, label = g_label), nudge_y = 0.3,  size = 2) +
  labs(title = "Habitants")+
  theme_ofce_void(base_size = 8)+ 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b=4)),
        plot.margin = margin(6,6,6,6),
        panel.background = element_rect(fill="grey97"))
