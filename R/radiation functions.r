# meaps <- function(rkdist, f, p, shuf = 1:nrow(rkdist), exact = FALSE) {
#   
#   n <- nrow(rkdist)
#   k <- ncol(rkdist)
#   
#   libre <- matrix(1, nrow=n+1, ncol=k)
#   emps <- matrix(0, nrow=n, ncol=k)
#   papn <- rep(k+1,k)
#   emp2i <- rep(0, k)
#   
#   fuite <- sum(f)
#   shuf[n+1] <- n+1
#   emp_max <- (n-fuite)/k
#   for (i in 1:n) {
#     i_o <- shuf[i]
#     i_o_plus_1 <- shuf[i+1]
#     picor <- p*libre[i_o,]
#     spic <- sum(picor)
#     if(spic>0)
#       ix <- -log(f[i_o])/spic
#     else
#       ix <- 1
#     xx <- ix
#     if(exact&spic>0)  xx <- uniroot(
#       function(x) prod(1-x*picor[rkdist[i_o,]]) - hab[i, "f"], 
#       interval = c(0, ix*10))$root
#     df <- abs(prod(1-xx*picor)-f[i_o])/f[i_o]
#     if(df>0.05) warning("pour {i}, l'erreur sur pf est de {df}" |> glue::glue())
#     picor_rk <- xx*picor[rkdist[i_o,]]
#     cpicor_rk <- cumprod(1-picor_rk)
#     cpicor_rk[2:k] <- cpicor_rk[1:(k-1)]
#     cpicor_rk[1] <- 1L
#     empi <- picor_rk*cpicor_rk
#     emps[i_o, ] <- empi[match(1:k,rkdist[i_o,])]
#     emp2i <- emp2i + emps[i_o, ]
#     libre[i_o_plus_1, ] <- emp2i<emp_max
#     papn <- xor(libre[i_o, ], libre[i_o_plus_1, ])
#   }
#   
#   return(list(
#     emp_meaps = emps,
#     occup = libre,
#     papn = papn
#   ))
# }

# fonction meaps à probabilité de disponibilité

meaps <- function(rkdist, f, p, shuf = 1:nrow(rkdist), exact = FALSE) {
  
  N <- nrow(rkdist)
  K <- ncol(rkdist)
  
  dispo <- matrix(1., nrow=N, ncol=K)
  emps <- matrix(0., nrow=N, ncol=K)
  papn <- matrix(0., nrow = N, ncol = K)
  emp2i <- rep(0., K)
  # on normalise p qui est maintenant un odd ratio relatif au moyen
  p <- p/mean(p) 
  fuite <- sum(f)
  emp_max <- (N-fuite)/K
  for (i in 1:N) {
    i_o <- shuf[i]
    dispo[i_o, ] <- pmax(0,1-emp2i/emp_max)
    papn[i,] <- dispo[i_o,]
    picor <- p*dispo[i_o, ]
    spicor <- sum(picor)
    if(spicor>0) {
      # ix est la chance moyenne
      ix <- -log(f[i_o])/spicor
      xx <- ix
      if(exact)  xx <- uniroot(
        function(x) prod(1-x/(1+x)*picor[rkdist[i_o,]]) - f[i_o], 
        interval = c(0, 100))$root
      df <- abs(prod(1-picor*xx/(1+xx))-f[i_o])/f[i_o]
      if(df>0.05) warning("pour {i}, l'erreur sur pf est de {df}" |> glue::glue())
      picor_rk <- xx/(1+xx)*picor[rkdist[i_o,]]
      cpicor_rk <- rep(1, K)
      for(j in 2:K) 
        cpicor_rk[j] <- cpicor_rk[j-1]*(1-picor_rk[j-1])
      empi <- picor_rk*cpicor_rk
      emps[i_o, ] <- empi[match(1:K,rkdist[i_o,])]
      emp2i <- emp2i + emps[i_o, ]
    }
  }
  
  return(emps)
}

meaps_summary <- function(emp, hab, dist, meaps, seuil_rang = 0.01) {
  ed <- meaps * dist
  emp_j <- matrixStats::colSums2(meaps)
  emp_i <- matrixStats::rowSums2(meaps)
  d_ind <- matrixStats::rowSums2(ed)/emp_i
  d_emp <- matrixStats::colSums2(ed)/emp_j
  # file <- matrixStats::colMeans2(meaps$dispo)
  # rangn <- apply((meaps$papn<seuil_rang),FUN=which.max, MARGIN = 2)
  
  return(list(
    hab = as_tibble(hab) |> mutate(d = d_ind, e_i = emp_i),
    emp = as_tibble(emp) |> mutate(d = d_emp, e_j = emp_j),
    meaps = meaps
  ))
}

rmeaps <- function(emp, hab, shuf = 1:nrow(hab), rcpp = TRUE, meaps_ver = 1) {
  k <- nrow(hab)
  n <- nrow(emp)
  ids <- rownames(hab)
  dist <- rdist::cdist(hab[,1:2], emp[,1:2])
  rkdist <- matrixStats::rowRanks(dist)
  if(meaps_ver==1) {
    if(rcpp)
      mm <- meaps_scpp(rkdist=rkdist, 
                       f = hab[, "f"], 
                       p = emp[, "p"], 
                       shuf = as.integer(shuf))
    else
      mm <- meaps(rkdist,
                  f = hab[, "f"], 
                  p = emp[, "p"], 
                  shuf = shuf)
  } 
  if(meaps_ver==2) {
    if(rcpp)
      mm <- meaps_rcpp(rkdist=rkdist, 
                       emplois = rep(1, n), 
                       actifs = rep(1, k),
                       odds =  emp[, "p"],
                       f = hab[, "f"],
                       shuf = as.integer(shuf))
    else
      mm <- meaps_odds(rkdist, 
                       emplois = rep(1, n), 
                       actifs = rep(1, k), 
                       odds = emp[, "p"],
                       f = hab[, "f"],
                       shuf = shuf)
  } 
  if(meaps_ver==4) {
    modds <- matrix(1, ncol=ncol(rkdist), nrow = nrow(rkdist))
    for (j in 1:ncol(rkdist)) modds[,j] <- emp[, "p"]
    mm <- rmeaps::meaps_oneshuf(
      rkdist=rkdist, 
      emplois = rep(1, n), 
      actifs = rep(1, k),
      modds =  modds,
      f = hab[, "f"],
      shuf = as.integer(shuf))
    mms <- meaps_summary(emp, hab, dist, mm)
    return(mms)
  }
}

pos_cunif <- function(n=100, centre = c(0.5, 0.5), rayon = 0.25) {
  rayons <- runif(n)^0.5*rayon
  angles <- runif(n)*2*pi
  res <- cbind(rayons*cos(angles),rayons*sin(angles))
  res[,1] <- res[,1]+centre[1]
  res[,2] <- res[,2]+centre[2]
  colnames(res) <- c("x", "y")
  rownames(res) <- 1:n
  return(res)
}

pos_cnorm  <- function(n=100, centre = c(0.5, 0.5), sigma = 0.1) {
  cbind(
    x = rnorm(n, mean=centre[1], sigma),
    y = rnorm(n, mean=centre[2], sigma))
}

# fonctions --------------------------
# 
add_total <- function(data) {
  un <- names(data)[[1]]
  tot <- data |>
    summarize(across(where(is.numeric), ~sum(.x))) |> 
    mutate({{un}}:= "total")
  bind_rows(data, tot)
}
make_tibs <- function(emp, hab) {
  dist <- rdist::cdist(hab[, 1:2], emp[,1:2])
  rkdist <- matrixStats::rowRanks(dist)
  f <- hab[, "f"]
  p <- emp[, "p"]
  hexhab <- hexbin::hexbin(hab[,"x"], hab[,"y"], xbins=bins, IDs=TRUE)@cID
  hexemp <- hexbin::hexbin(emp[,"x"], emp[,"y"], xbins=bins, IDs=TRUE)@cID
  habs <- hab |> 
    as_tibble() |> 
    mutate(hab = 1:nrow(hab),
           hhex = hexhab)
  hhex <- habs |> 
    group_by(hhex) |> 
    summarize(nh = n(),
              gh = names(table(g))[[1]],
              x = mean(x), y=mean(y))
  hgroupes <-  habs |> group_by(g) |> 
    summarize(x = mean(x),
              y= mean(y),
              pop = n())  |> 
    mutate(g_label = str_c("h", g ," " ,pop ," habitants"), size = 6/.pt)
  
  emps <- emp |> 
    as_tibble() |> 
    mutate(emp = 1:nrow(emp),
           ehex = hexemp)
  ehex <- emps |> 
    group_by(ehex) |> 
    summarize(ne = n(),
              ge = names(table(g))[[1]],
              x = mean(x), y=mean(y))
  egroupes <-  emps |> group_by(g) |> 
    summarize(x = mean(x),
              y= mean(y),
              pop = n()) |> 
    mutate(g_label = str_c("e", g ," " ,pop ," emplois"), size = 6/.pt)
  list(habs=habs, emps=emps,
       ehex = ehex, hhex=hhex,
       hexhab = hexhab, hexemp = hexemp,
       egroupes = egroupes, hgroupes = hgroupes, 
       dist = dist, rk = rkdist, f = f, p = p)
}

rmeaps_bstp <- function(scn, shufs, workers=1) {
  pl <- future::plan()
  sp <- split(1:nrow(shufs), ceiling(1:nrow(shufs)/max(1,(nrow(shufs)/workers))))
  future::plan("multisession", workers=min(length(sp), workers))
  res <- furrr::future_map(sp,~{
    Rcpp::sourceCpp("R/meaps2.cpp", showOutput = FALSE, echo = FALSE, verbose=FALSE)
    rr <- meaps_boot(scn$rk, 
                     rep(1,nrow(scn$emp)), 
                     rep(1,nrow(scn$hab)),
                     scn$p,
                     scn$f,
                     shufs[.x, , drop=FALSE]) # attention c'est divisé par le nombre de tirages
    rr <- map(rr, function(rrr) rrr * length(.x))
  }, .options = furrr::furrr_options(seed = TRUE))
  future::plan(pl)
  res <- purrr::transpose(res)
  res <- purrr::map(res, function(rr) reduce(rr, `+`)/nrow(shufs))
  return(res)
} 

rmeaps_multishuf <- function(scn, shufs, nthreads=0) {
  k <- nrow(hab)
  n <- nrow(emp)
  ids <- rownames(hab)
  dist <- rdist::cdist(hab[,1:2], emp[,1:2])
  rkdist <- matrixStats::rowRanks(dist)
  modds <- matrix(1, ncol=ncol(rkdist), nrow = nrow(rkdist))
  for (j in 1:ncol(rkdist)) 
    modds[,j] <- scn$p[[j]]
  rr <- rmeaps::meaps_multishuf(
    rkdist = scn$rk, 
    emplois = rep(1,nrow(scn$emp)), 
    actifs = rep(1,nrow(scn$hab)),
    modds = modds,
    f = scn$f,
    shuf = shufs,
    nthreads = nthreads) 
  mms <- meaps_summary(emp, hab, dist, rr)
  return(mms)
}

emp_flux <- function(s, emp, empec = NULL) {
  g_col <- unique(s$emps$g) |> sort()
  g_row <- unique(s$habs$g) |> sort()
  icol <- map(g_col, function(gr) which(s$emps$g==gr))
  irow <- map(g_row, function(gr) which(s$habs$g==gr))
  emp_red <- do.call(cbind, map(icol, function(col) matrixStats::rowSums2(emp[,col])))
  emp_red <- do.call(rbind, map(irow, function(row) matrixStats::colSums2(emp_red[row,])))
  dimnames(emp_red) <- list(g_row, g_col)
  if(!is.null(empec)) {
    emps2 <- empec^2
    emps2_red <- do.call(cbind, map(icol, function(col) matrixStats::rowSums2(emps2[,col])))
    emps2_red <- do.call(rbind, map(irow, function(row) matrixStats::colSums2(emps2_red[row,])))
    emps2_red <- sqrt(emps2_red)
    dimnames(emps2_red) <- list(g_row, g_col)
  } else {
    emps2_red <- NULL
  }
  return(list(s = emp_red, ec = emps2_red))
}