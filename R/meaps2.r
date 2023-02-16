meaps_odds <- function(rkdist, emplois, actifs, odds = 1, f, shuf = 1:nrow(rkdist), debordement = FALSE) {
  
  n <- nrow(rkdist)
  k <- ncol(rkdist)
  
  liaisons <- matrix(0, nrow = n, ncol = k)
  
  # recyclage forcé si besoin
  if (length(emplois) == 1) emplois <- rep(emplois, k)
  if (length(actifs) == 1) actifs <- rep(actifs, n)
  if (length(odds) == 1) odds <- rep(odds, k)
  if (length(f) == 1) f <- rep(f, n)
  
  # on décide de caler les emplois sur les actifs restants
  emplois <- emplois / sum(emplois) * sum(actifs *(1-f)) 
  
  for (i in shuf) {
    
    arrangement <- match(1:k, rkdist[i, ]) |> discard(is.na)
    la <- length(arrangement) 
    
    dispo <- (!is.na(rkdist[i,])) * emplois # la quantité d'emplois restants libres et accessibles pour i
    
    tot <- sum(dispo)
    if (tot <= 1e-1) next 
    
    p_ref <- 1 - f[i]^(1/tot)
    chances_ref <- p_ref / (1 - p_ref)
    chances_abs <- odds * chances_ref
    
    if (min(chances_abs) <= 0 | max(chances_abs) > 10000) cat(range(chances_abs), "\n")
    # 
    # fuite_f <- function(x) sum(dispo * log(1 + x * chances_abs)) + log(f[i]) 
    # fuite_derivee <- function(x) sum((dispo * chances_abs) / ( 1 + x * chances_abs))
    # 
    # xn = - log(f[i]) / (sum((dispo * odds))) / chances_ref 
    # xnp1 = xn - fuite_f(xn)/fuite_derivee(xn)
    # while (abs(xnp1 - xn) > 1e-9) {
    #   xn = xnp1
    #   xnp1 = xn - fuite_f(xn)/fuite_derivee(xn)
    # }
    
    # chances_abs <- xnp1 * odds * chances_ref
    p_abs <- chances_abs / (1 + chances_abs)
    
    passe <- (1 - p_abs[arrangement])^dispo[arrangement] |> cumprod()
    reste <- lag(passe, default = 1) - passe
    
    repartition <- actifs[i] * reste[rkdist[i,]]
    repartition <- replace(repartition, is.na(repartition), 0)
    
    if (debordement) {
      for (j in 1:(la-1)) if (emplois[arrangement[j]] < repartition[arrangement[j]]) {
        repartition[arrangement[j+1]] <- repartition[arrangement[j+1]] + repartition[arrangement[j]] - emplois[arrangement[j]]
        repartition[arrangement[j]] <- emplois[arrangement[j]]
      }
      # rq : on ne peut pas induire de décalage pour congestion sur le dernier de la file d'attente
      if (emplois[arrangement[la]] < repartition[arrangement[la]]) repartition[arrangement[la]] <- emplois[arrangement[la]] 
    }
    
    liaisons[i, ] <- repartition
    emplois <- emplois - repartition
    
  }
  
  return(liaisons)
}

rmeaps2_boot <- function(rkdist, emplois, actifs, odds = 1, f, shufs, workers=1, meaps_cpp = "meaps2.cpp") {
  pl <- future::plan()
  sp <- split(1:nrow(shufs), ceiling(1:nrow(shufs)/max(1,(nrow(shufs)/workers))))
  future::plan("multisession", workers=min(length(sp), workers))
  res <- furrr::future_map(sp,function(ss) {
    Rcpp::sourceCpp(meaps_cpp)
    rr <- meaps_boot(rkdist, emplois, actifs, odds = odds, f=f, shufs[ss, , drop=FALSE]) # attention c'est divisé par le nombre de tirages
    rr <- map(rr, function(rrr) rrr * length(ss))
  }, .options = furrr::furrr_options(seed = TRUE))
  future::plan(pl)
  res <- purrr::transpose(res)
  res <- purrr::map(res, function(rr) reduce(rr, `+`)/nrow(shufs))
  return(res)
} 

rmeaps_bootmat <- function(
    rkdist, emplois, actifs, 
    modds = matrix(1, nrow=nrow(rkdist), ncol=ncol(rkdist)), 
    f, shufs, workers=1) {
  pl <- future::plan()
  sp <- split(1:nrow(shufs), ceiling(1:nrow(shufs)/max(1,(nrow(shufs)/workers))))
  future::plan("multisession", workers=min(length(sp), workers))
  res <- furrr::future_map(sp,~{
    Rcpp::sourceCpp("R/meaps_oddmatrix.cpp")
    rr <- meaps_boot2(rkdist, emplois, actifs, modds = modds, f=f, shufs[.x, , drop=FALSE]) # attention c'est divisé par le nombre de tirages
    rr <- map(rr, function(rrr) rrr * length(.x))
  }, .options = furrr::furrr_options(seed = TRUE))
  future::plan(pl)
  res <- purrr::transpose(res)
  res <- purrr::map(res, function(rr) reduce(rr, `+`)/nrow(shufs))
  return(res)
} 


