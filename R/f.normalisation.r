# Cette fonction permet d'agréger une matrice selon des regroupements lignes et colonnes.
# Elle sert notamment à agréger selon les communes.
emp_flux2 <- function(m, row_groups, col_groups) {
  icol <- map(unique(col_groups), function(gr) which(col_groups==gr))
  irow <- map(unique(row_groups), function(gr) which(row_groups==gr))
  m_red <- do.call(cbind, map(icol, function(col) matrixStats::rowSums2(m[,col, drop=FALSE])))
  m_red <- do.call(rbind, map(irow, function(row) matrixStats::colSums2(m_red[row, , drop=FALSE])))
  
  colnames(m_red) <- unique(col_groups)
  rownames(m_red) <- unique(row_groups)
  return(m_red)
}

# Cette fonction extrait des résultats de meaps les flux par communes,
# et calcul le odds ratio de ce flux avec celui de MOBPRO.
communaliser <- function(MOD, eps = 0.1) {
  # le paramètre esp est un incrément dans chaque case pour le calcul des or
  # de manière à éviter les divisions par zéro. Le biais introduit est d'atténuer
  # la force des or, ce qui n'est pas trop génant.
  liens <- MOD$emps
  dimnames(liens) <- dimnames(mat_rang)
  
  res <- emp_flux2(liens, 
                   les_actifs$COMMUNE,
                   les_emplois$DCLT) |> 
    as.data.table(keep.rownames = "COMMUNE")
  
  res <- melt(res,
              id.vars = "COMMUNE",
              variable.name = "DCLT", 
              value.name = "flux")
  res <- merge(res[flux > 0, ],
               mobpro, by = c("COMMUNE", "DCLT"),
               all = TRUE)
  res <- res[, flux := replace(flux, is.na(flux), 0)]
  res <- res[, mobpro := replace(mobpro, is.na(mobpro), 0)]
  
  res[,
      `:=`(marge_flux = sum(flux + eps),
           marge_mob = sum(mobpro + eps)),
      by = "COMMUNE"]
  res[,  `:=`(c_flux = (flux + eps) / (marge_flux - flux - eps), 
              c_mob = (mobpro + eps) / (marge_mob - mobpro - eps))]
  res[, 
      or := c_mob / c_flux][,
                            ':='(marge_flux = NULL, marge_mob = NULL, 
                                 c_flux = NULL, c_mob = NULL)]
  setorder(res, -mobpro)
  return(res[])
}

communaliser2 <- function(MOD, les_actifs, les_emplois, rangs, eps = 0.1) {
  # le paramètre esp est un incrément dans chaque case pour le calcul des or
  # de manière à éviter les divisions par zéro. Le biais introduit est d'atténuer
  # la force des or, ce qui n'est pas trop génant.
  liens <- MOD$emps
  dimnames(liens) <- dimnames(rangs)
  
  res <- emp_flux2(liens, 
                   les_actifs$COMMUNE,
                   les_emplois$DCLT) |> 
    as.data.table(keep.rownames = "COMMUNE")
  
  res <- melt(res,
              id.vars = "COMMUNE",
              variable.name = "DCLT", 
              value.name = "flux")
  res <- merge(res[flux > 0, ],
               mobpro, by = c("COMMUNE", "DCLT"),
               all = TRUE)
  res <- res[, flux := replace(flux, is.na(flux), 0)]
  res <- res[, mobpro := replace(mobpro, is.na(mobpro), 0)]
  
  res[,
      `:=`(marge_flux = sum(flux + eps),
           marge_mob = sum(mobpro + eps)),
      by = "COMMUNE"]
  res[,  `:=`(c_flux = (flux + eps) / (marge_flux - flux - eps), 
              c_mob = (mobpro + eps) / (marge_mob - mobpro - eps))]
  res[, 
      or := c_mob / c_flux][,
                            ':='(marge_flux = NULL, marge_mob = NULL, 
                                 c_flux = NULL, c_mob = NULL)]
  setorder(res, -mobpro)
  return(res[])
}


# Cette fonction prend un vecteur de odds de flux entre les communes
# et calcul la matrice des odds au carreau.
creer_matrice_odds <- function(x, or = "or") {
  odds <- expand_grid(les_actifs |> select(fromidINS = idINS, COMMUNE), 
                      les_emplois |> select(toidINS = idINS, DCLT)) |> 
    left_join(x |> select(COMMUNE,
                          DCLT, 
                          or = all_of(or)), 
              by = c("COMMUNE", "DCLT")) |> 
    mutate(or = replace(or, is.na(or), 1)) |> 
    pull(or)
  
  matrix(odds, nrow = N, ncol = K, byrow = TRUE)
}

NAsifie_distance <- function(dist, les_actifs, les_emplois, mobpro) {
  nas <- expand_grid(
    les_actifs |> select(fromidINS = idINS, COMMUNE), 
    les_emplois |> select(toidINS = idINS, DCLT)) |> 
    left_join(mobpro, by=c("COMMUNE", "DCLT")) |> 
    transmute(fromidINS, toidINS, mobpro = is.na(mobpro))
  nas <- matrix(
    nas$mobpro, 
    byrow=TRUE, 
    ncol = nrow(les_emplois), 
    nrow = nrow(les_actifs))
  dist[nas] <- NA
  return(dist)
}

# Un bête calcul de khi2,
khi2 <- function(est, ref) {
  sum((est - ref)^2/ref, na.rm=TRUE)
}

# Calcul de la béta-divergence Kullback-Leibler, soit l'entropie relative,
# avec une variation sur le traitement des 0
# ici on ajoute un petit bruit au zéro
kullback_leibler_2 <- function(nb_mod, nb_ref, seuil_collapse = .1, bruit = 1e-6) {
  tib <- tibble(x = nb_mod, y = nb_ref) |> 
    mutate(
      px = x / sum(x, na.rm=TRUE),
      py = y / sum(y, na.rm=TRUE)) |> 
    mutate(
      px = ifelse(px<=bruit/length(px[!is.na(px)]), 
                  bruit/length(px[!is.na(px)]),
                  px),
      px = px / sum(px, na.rm=TRUE),
      py = ifelse(py<=bruit/length(py[!is.na(py)]),
                  bruit/length(py[!is.na(py)]), 
                  py),
      py = py / sum(py, na.rm=TRUE)
    )
  sum(tib$py * log(tib$py / tib$px), na.rm=TRUE)
}

kl <- function(x,y) {
  yn <- y/sum(y)
  xn <- x/sum(x)
  sum( yn * log(yn/xn) )
}

# on ôte les 0
# 
# 
kullback_leibler_3 <- function(nb_mod, nb_ref, seuil_collapse = .1, bruit = 1e-6) {
  tib <- tibble(x = nb_mod, y = nb_ref) |>
    filter(x!=0, y!=0, !is.na(x), !is.na(y)) |> 
    mutate(
      px = x / sum(x),
      py = y / sum(y)) 
  sum(tib$py * log(tib$py / tib$px))
}

# Calcul de la béta-divergence Kullback-Leibler, soit l'entropie relative,
kullback_leibler <- function(nb_mod, nb_ref, seuil_collapse = 1, bruit = 1e-6) {
  tib <- tibble(x = nb_mod, y = nb_ref) |> 
    mutate(
      g=cur_group_id(),
      g = ifelse( pmin(nb_ref, nb_mod) >= seuil_collapse, g, 0)
    ) |>
    group_by(g) |> 
    summarise(x = sum(x, na.rm=TRUE), y = sum(y, na.rm=TRUE)) |> 
    ungroup() |> 
    mutate(px = x / sum(x, na.rm=TRUE), 
           py = y / sum(y, na.rm=TRUE))
  sum(tib$py * log(tib$py / tib$px), na.rm=TRUE)
}

r2kl <- function(nb_mod, nb_ref, seuil_collapse = .001, bruit=1e-6, version=1) {
  kl_f <- switch(version,
                 "1" = kullback_leibler,
                 "2" = kullback_leibler_2,
                 "3" = kullback_leibler_3)
  klref <- kl_f(mean(nb_ref, na.rm=TRUE), nb_ref, 
                seuil_collapse=seuil_collapse, bruit=bruit)
  kl <- kl_f(nb_mod, nb_ref, 
             seuil_collapse=seuil_collapse, bruit=bruit)
  return(1-kl/klref)
}
# évolution de r2kl
# on prend 3 en défaut, en comptant sur peu de 0 (les distances sont NAsifiées)
# on évalue l'impact potentiel du bruit en collapsant les petites cases
# 
r2kl2 <- function(nb_mod, nb_ref, seuil = .99, bruit=1e-6) {
  tib <- tibble(x=nb_mod, y=nb_ref) |> 
    filter(x!=0, y!=0)
  r2kln0 <- 1-kl(tib$x, tib$y)/kl(rep(mean(tib$y), nrow(tib)), tib$y)
  tib <- tib |> 
    arrange(desc(y)) |> 
    mutate(cum = cumsum(y)/sum(y)) |> 
    filter(cum<=seuil)
  r2klnb <- 1-kl(tib$x, tib$y)/kl(rep(mean(tib$y), nrow(tib)), tib$y)
  return(list(r2kl = r2kln0, r2kl_l = r2klnb))
}

shanon <- function(nb_ref, seuil_collapse = 1) {
  tib <- tibble(y = nb_ref) |> 
    mutate(
      g = cur_group_rows(),
      g = ifelse(nb_ref >= seuil_collapse, g, 0)
    ) |>
    group_by(g) |> 
    summarise( y = sum(y), na.rm=TRUE) |> 
    mutate(py = y / sum(y, na.rm=TRUE))
  
  sum(tib$py * log(tib$py), na.rm=TRUE)
}
# Compare différents modèles issus de meaps.
comparer_meaps_res <- function(...) {
  list_mod <- enquos(...)
  noms <- as.character(list_mod) |> str_replace("~", "")
  list_mod <- map(list_mod, rlang::eval_tidy)
  
  return(map2_dfr(list_mod, noms, ~ {
    setDT(.x)
    kl <- kullback_leibler_2(.x$flux, .x$mobpro)
    tibble(
      model = .y,
      khi2 = khi2(.x$flux + .1, .x$mobpro + .1), 
      kullback_leibler = kl,
      r2kl = r2kl2(.x$flux, .x$mobpro)$r2kl,
      dispersion_odds = sum(abs(log(.x$or))),
      delta_LR_entrant =
        .x[DCLT == 17300, sum(flux - mobpro, na.rm=TRUE)/sum(mobpro, na.rm=TRUE)],
      delta_LR_sortant = 
        .x[COMMUNE == 17300, sum(flux - mobpro, na.rm=TRUE)/sum(mobpro, na.rm=TRUE)]
    )}
  ))
}

comparer_meaps_bygrp <- function(data) {
  data |>
    mutate(nz = (mobpro!=0 & flux != 0)) |> 
    summarize(
    khi2 = khi2(flux + .1, mobpro + .1),
    kullback_leibler = kl(flux[nz], mobpro[nz]),
    r2kl = r2kl2(flux, mobpro)$r2kl,
    dispersion_odds = sum(abs(log(or[nz]))),
    delta_LR_entrant = sum(flux[DCLT == 17300] - mobpro[DCLT == 17300], na.rm=TRUE)/
      sum(mobpro[DCLT == 17300], na.rm=TRUE),
    delta_LR_sortant = sum(flux[COMMUNE == 17300] - mobpro[COMMUNE == 17300], na.rm=TRUE)/
      sum(mobpro[COMMUNE == 17300], na.rm=TRUE),
    nnz = sum(nz))
}