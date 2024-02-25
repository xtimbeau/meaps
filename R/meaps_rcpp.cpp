#include <Rcpp.h>
#include <valarray>
#include <algorithm>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
// retourne plus de data que rcpp
List meaps_cpp(
    IntegerMatrix rkdist, 
    NumericVector f, 
    NumericVector p,
    IntegerVector shuf)
{
  int N = rkdist.nrow();
  int K = rkdist.ncol();
  
  if(f.length()!=N) {stop("erreur f n'a pas la bonne dimension");}
  if(p.length()!=K) {stop("erreur p n'a pas la bonne dimension");}
  if(shuf.length()!=N) {stop("erreur shuf n'a pas la bonne dimension");}
  
  int i_o;
  double fuite, emp_max, xx, spicor;
  NumericVector picor(K), empi(K), empi_inv(K), picor_rk(K), cpicor_rk(K), emp2i(K), disp(K), pn(K);
  IntegerVector ishuf(N);
  NumericMatrix dispo(N, K);
  IntegerVector rkinv(K), rk(K);
  NumericMatrix emps(N, K), papn(N,K);
  
  // init
  // décalage des index r->C
  ishuf = shuf-1L;
  
  // normalisation des p
  pn = p/mean(p);
  
  fuite = std::accumulate(f.begin(), f.end(), 0.0);
  emp_max = (N-fuite)/K;
  
  // boucle principale sur les lignes (habitants)  
  for (auto i = 0; i < N; i++) {
    // on met à jour les dispos
    disp = 1-emp2i/emp_max;
    disp[disp<0] = 0;
    // on suit la liste de priorité 
    i_o = ishuf[i];
    dispo(i_o, _) = disp;
    papn(i, _) = disp;
    // on met les probabilité à jour en fonction de dispo et des odd ratios
    picor = pn * disp;
    spicor = std::accumulate(picor.begin(), picor.end(), 0.0);
    // on calcule la modif de la fuite
    if(spicor>0) {
      xx = -log(f[i_o])/spicor;
      
      // on cumule les proba (attention xx a changé de signification, c'est un odd ratio)
      cpicor_rk[0L] = 1;
      rk = rkdist(i_o,_)-1L;
      picor_rk = as<NumericVector>(picor[rk])*xx/(1+xx);
      for (auto k = 1; k<K; k++) {
        cpicor_rk[k] = cpicor_rk[k-1]*(1- picor_rk[k-1]);
      }
      
      empi = cpicor_rk * picor_rk;
      
      // on calcule les emplois, remis dans l'ordre canonique
      for (auto k = 0; k<K; k++) {
        rkinv[rk[k]] = k;
      }
      empi_inv = empi[rkinv];
      emp2i = emp2i + empi_inv;
      emps(i_o, _) = empi_inv;
    }
  }
  
  // zou
  return List::create(
    _["emps"] = emps,
    _["dispo"] = dispo,
    _["papn"] = papn);
}


// [[Rcpp::export]]
NumericMatrix meaps_scpp(
    IntegerMatrix rkdist, 
    NumericVector f, 
    NumericVector p,
    IntegerVector shuf)
{
  int N = rkdist.nrow();
  int K = rkdist.ncol();
  
  if(f.length()!=N) {stop("erreur f n'a pas la bonne dimension");}
  if(p.length()!=K) {stop("erreur p n'a pas la bonne dimension");}
  if(shuf.length()!=N) {stop("erreur shuf n'a pas la bonne dimension");}
  
  int i_o;
  double fuite, emp_max, xx, spicor;
  NumericVector picor(K), empi(K), empi_inv(K), picor_rk(K), cpicor_rk(K), emp2i(K), disp(K), pn(K);
  IntegerVector ishuf(N);
  NumericMatrix dispo(N, K);
  IntegerVector rkinv(K), rk(K);
  NumericMatrix emps(N, K);
  
  // init
  // décalage des index r->C
  ishuf = shuf-1L;
  
  // normalisation des p
  pn = p/mean(p);
  
  fuite = std::accumulate(f.begin(), f.end(), 0.0);
  emp_max = (N-fuite)/K;
  
  // boucle principale sur les lignes (habitants)  
  for (auto i = 0; i < N; i++) {
    // on met à jour les dispos
    disp = 1-emp2i/emp_max;
    disp[disp<0] = 0;
    // on suit la liste de priorité 
    i_o = ishuf[i];
    dispo(i_o, _) = disp;
    
    // on met les probabilité à jour en fonction de dispo et des odd ratios
    picor = pn * disp;
    spicor = std::accumulate(picor.begin(), picor.end(), 0.0);
    // on calcule la modif de la fuite
    xx = -log(f[i_o])/spicor;
    
    // on cumule les proba (attention xx a changé de signification, c'est un odd ratio)
    cpicor_rk[0L] = 1;
    rk = rkdist(i_o,_)-1L;
    picor_rk = as<NumericVector>(picor[rk])*xx/(1+xx);
    for (auto k = 1; k<K; k++) {
      cpicor_rk[k] = cpicor_rk[k-1]*(1- picor_rk[k-1]);
    }
    
    empi = cpicor_rk * picor_rk;
    
    // on calcule les emplois, remis dans l'ordre canonique
    for (auto k = 0; k<K; k++) {
      rkinv[rk[k]] = k;
    }
    empi_inv = empi[rkinv];
    emp2i = emp2i + empi_inv;
    emps(i_o, _) = empi_inv;
  }
  
  // zou
  return emps;
}

// [[Rcpp::export]]
List meaps_btsp(
    List scenario,
    IntegerMatrix shufs)
{
  IntegerMatrix rkdist = as<Rcpp::IntegerMatrix>(scenario["rk"]);
  NumericMatrix dist = as<Rcpp::NumericMatrix>(scenario["dist"]);
  NumericVector f = as<Rcpp::NumericVector>(scenario["f"]);
  NumericVector p = as<Rcpp::NumericVector>(scenario["p"]);
  int N = rkdist.nrow();
  int K = rkdist.ncol();
  int Ns = shufs.nrow();
  IntegerVector shuf(N);
  NumericMatrix res(N,K), acc(N,K), acc2(N,K);
  
  for(auto s = 0; s < Ns; s++ ) {
    // on vérifie qu'on est pas trop long
    Rcpp::checkUserInterrupt();
    res = meaps_scpp(rkdist, f, p, shufs(s, Rcpp::_));
    for(auto ii=0; ii<N; ii++) {
      acc(ii,_) = acc(ii,_) + res(ii,_);
      acc2(ii,_) = acc2(ii,_) + res(ii,_)*res(ii,_); 
    }
  }
  
  acc = acc / Ns;
  acc2 = acc2 / Ns;
  
  return List::create(
    _["emps"] = acc,
    _["emps2"] = acc2
  );
}




// List meaps_btsp(
//     List scenario,
//     IntegerMatrix shufs,
//     int threads)
// {
//   IntegerMatrix rkdist = as<Rcpp::IntegerMatrix>(scenario["rk"]);
//   NumericMatrix dist = as<Rcpp::NumericMatrix>(scenario["dist"]);
//   NumericVector f = as<Rcpp::NumericVector>(scenario["f"]);
//   NumericVector p = as<Rcpp::NumericVector>(scenario["p"]);
//   int N = rkdist.nrow();
//   int K = rkdist.ncol();
//   int Ns = shufs.nrow();
//   int st = floor(Ns/threads);
//   IntegerVector shuf(N);
//   List res(threads), res2(threads);
//   NumericMatrix acc(N,K), acc2(N,K);
//   
//   omp_set_num_threads(threads);
//   
// #pragma omp parallel for
//   for(auto t=0; t<threads; t++) {
//     NumericMatrix local_res(N,K), local_acc(N,K), local_acc2(N,K);
//     for(auto s = t*st; s < std::max((t+1)*st, Ns); s++ ) {
//       local_res = meaps_rcpp(rkdist, f, p, shufs(s, Rcpp::_));
//       for(auto ii=0; ii<N; ii++) {
//         local_acc(ii,_) = local_acc(ii,_) + local_res(ii,_);
//         local_acc2(ii,_) = local_acc2(ii,_) + local_res(ii,_)*local_res(ii,_); 
//       }
//     }
//     res[t] = clone(local_acc);
//     res2[t] = clone(local_acc2);
//   }
//   
//   for(auto t=0; t<threads; t++) {
//     for(auto i=0; i<N; i++) {
//       NumericMatrix lr = res[t];
//       NumericMatrix lr2 = res2[t];
//       acc(i,_) = acc(i,_) + lr(i,_);
//       acc2(i, _) = acc2(i, _) + lr2(i,_);
//     }
//   }
//   
//   acc = acc / Ns;
//   acc2 = acc2 / Ns;
//   
//   return List::create(
//     _["emps"] = acc,
//     _["emps2"] = acc2
//   );
// }

// valarray<double> meaps_ts(
//     valarray<int> rkdist, 
//     valarray<double> f, 
//     valarray<double> p,
//     valarray<int> shuf)
// {
//   int N = f.size();
//   int K = p.size();
//   
//   int i_o;
//   double fuite, emp_max, xx, spicor;
//   valarray<double> picor(K), empi(K), empi_inv(K), picor_rk(K), cpicor_rk(K), emp2i(K), disp(K);
//   valarray<int> ishuf(N);
//   valarray<double> dispo(N*K), emps(N*K);
//   valarray<int> unak(K), rkinv(K);
//   
//   // init
//   // décalage des index r->C
//   ishuf = shuf-1L;
//   for(auto k=0; k<K; k++) {
//     unak = k+1;
//     }
//   
//   emp2i.fill(0.);
//   
//   fuite = std::accumulate(f.begin(), f.end(), 0.0);
//   emp_max = (N-fuite)/K;
//   
//   // boucle principale sur les lignes (habitants)  
//   for (int i = 0; i < N; i++) {
//     // on vérifie qu'on est pas trop long
//     if(i%1000 == 1) {Rcpp::checkUserInterrupt();}
//     // on met à jour les dispos
//     disp = 1-emp2i/emp_max;
//     disp[disp<0] = 0;
//     // on suit la liste de priorité 
//     i_o = ishuf[i];
//     dispo(i_o, _) = disp;
//     
//     // on met les probabilité à jour en fonction de dispo
//     picor = p * dispo(i_o, _);
//     spicor = std::accumulate(picor.begin(), picor.end(), 0.0);
//     // on calcule la modif de la fuite
//     xx = -log(f[i_o])/spicor;
//     
//     // on cumule les proba
//     cpicor_rk[0L] = 1;
//     picor_rk = as<NumericVector>(picor[rkdist(i_o, _) - 1L])*xx;
//     for (int k = 1; k<K; k++) {
//       cpicor_rk[k] = cpicor_rk[k-1]*(1- picor_rk[k-1]);
//     }
//     
//     empi = cpicor_rk * picor_rk;
//     
//     // for (int k = 1; k<K; k++) {
//     //   picor_rk[k] = xx*picor[rkdist(i_o,k)-1];
//     //   cpicor_rk[k] = cpicor_rk[k-1]*(1-picor_rk[k-1]);
//     //   empi[k] = cpicor_rk[k] * picor_rk[k];
//     // }
//     
//     // on calcule les emplois, remis dans l'ordre canonique
//     // on en profite pour remplir le résultat
//     rkinv = Rcpp::match(unak, rkdist(i_o,_))-1L;
//     empi_inv = empi[rkinv];
//     emp2i = emp2i + empi_inv;
//     emps(i_o, _) = empi_inv;
//     // for (int k = 0; k<K; k++) {
//     //   emp2i[k] = emp2i[k] + empi[rkinv(k)-1];
//     //   emps(i_o,k) = empi[rkinv(k)-1];
//     //   }
//   }
//   
//   // zou
//   return emps;
// }
