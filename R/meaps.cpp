#include <Rcpp.h>
#include <valarray>
#include <algorithm>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List meaps_cpp_v1(
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
  std::valarray<double> picor(K), empi(K), picor_rk(K), cpicor_rk(K), emp2i(K);
  std::valarray<int> ishuf(N);
  std::valarray<double> dispo(N*K);
  IntegerVector unak(K), rkinv(K);
  NumericMatrix out_emps(N, K);
  NumericMatrix out_papn(N, K);
  NumericMatrix out_dispo(N, K);
  
  // init
  // décalage des index r->C
  for(int i=0; i<N; i++) {
    ishuf[i] = shuf[i]-1;
  }
  
  for(int j=0; j<K; j++) {
    dispo[ishuf[0]*K+j] = 1;
    out_dispo(ishuf[0], j) = 1;
    unak[j] = j+1;
  }
  
  for (int j = 0; j < K; j++) {
    emp2i[j] = 0;
  }
  
  fuite = 0.;
  for (int i=0; i<N; i++) {
    fuite = fuite + f[i];
  }
  emp_max = (N-fuite)/K;
  
  // boucle principale sur les lignes (habitants)  
  for (int i = 0; i < N; i++) {
    // on vérifie qu'on est pas trop long
    if(i%1000 == 1) {Rcpp::checkUserInterrupt();}
    // on suit la liste de priorité 
    i_o = ishuf[i];
    for (auto k=0; k<K; k++) {
      dispo[i_o*K+k] = 1-emp2i[k]/emp_max;
      if(dispo[i_o*K+k]<0) {dispo[i_o*K+k]=0;}
    }
    // on met les probabilité à jour en fonction de dispo
    spicor = 0;
    for (int k = 0; k<K; k++) {
      picor[k] = p[k] * dispo[i_o*K+k];
      spicor = spicor + picor[k];
    }
    // on caclcule la modif de la fuite
    xx = -log(f[i_o])/spicor;
    
    // on cumule les proba
    picor_rk[0] = xx*picor[rkdist(i_o,0)-1];
    cpicor_rk[0] = 1;
    empi[0] = picor_rk[0];
    for (int k = 1; k<K; k++) {
      picor_rk[k] = xx*picor[rkdist(i_o,k)-1];
      cpicor_rk[k] = cpicor_rk[k-1]*(1-picor_rk[k-1]);
      empi[k] = cpicor_rk[k] * picor_rk[k];
    }
    
    // on calcule les emplois, remis dans l'ordre canonique
    // on en profite pour remplir le résultat
    rkinv = Rcpp::match(unak, rkdist(i_o, Rcpp::_));
    for (int k = 0; k<K; k++) {
      emp2i[k] = emp2i[k] + empi[rkinv(k)-1];
      out_papn(i, k) = dispo[i_o*K+k]; 
      out_dispo(i_o, k) = dispo[i_o*K+k];
      out_emps(i_o,k) = empi[rkinv(k)-1];
    }
  }
  
  // zou
  return List::create(
    _["emps"] = out_emps,
    _["dispo"] = out_dispo,
    _["papn"] = out_papn
  );
}

// [[Rcpp::export]]
NumericMatrix meaps_cpp_emp(
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
  std::valarray<double> picor(K), empi(K), picor_rk(K), cpicor_rk(K), emp2i(K);
  std::valarray<int> ishuf(N);
  std::valarray<double> dispo(N*K);
  IntegerVector unak(K), rkinv(K);
  NumericMatrix out_emps(N, K);
  // NumericMatrix out_papn(N, K);
  // NumericMatrix out_dispo(N, K);
  // 
  // init
  // décalage des index r->C
  for(int i=0; i<N; i++) {
    ishuf[i] = shuf[i]-1;
  }
  
  for(int j=0; j<K; j++) {
    dispo[ishuf[0]*K+j] = 1;
    // out_dispo(ishuf[0], j) = 1;
    unak[j] = j+1;
  }
  
  for (int j = 0; j < K; j++) {
    emp2i[j] = 0;
  }
  
  fuite = 0.;
  for (int i=0; i<N; i++) {
    fuite = fuite + f[i];
  }
  emp_max = (N-fuite)/K;
  
  // boucle principale sur les lignes (habitants)  
  for (int i = 0; i < N; i++) {
    // on vérifie qu'on est pas trop long
    if(i%1000 == 1) {Rcpp::checkUserInterrupt();}
    // on suit la liste de priorité 
    i_o = ishuf[i];
    for (auto k=0; k<K; k++) {
      dispo[i_o*K+k] = 1-emp2i[k]/emp_max;
      if(dispo[i_o*K+k]<0) {dispo[i_o*K+k]=0;}
    }
    // on met les probabilité à jour en fonction de dispo
    spicor = 0;
    for (int k = 0; k<K; k++) {
      picor[k] = p[k] * dispo[i_o*K+k];
      spicor = spicor + picor[k];
    }
    // on caclcule la modif de la fuite
    xx = -log(f[i_o])/spicor;
    
    // on cumule les proba
    picor_rk[0] = xx*picor[rkdist(i_o,0)-1];
    cpicor_rk[0] = 1;
    empi[0] = picor_rk[0];
    for (int k = 1; k<K; k++) {
      picor_rk[k] = xx*picor[rkdist(i_o,k)-1];
      cpicor_rk[k] = cpicor_rk[k-1]*(1-picor_rk[k-1]);
      empi[k] = cpicor_rk[k] * picor_rk[k];
    }
    
    // on calcule les emplois, remis dans l'ordre canonique
    // on en profite pour remplir le résultat
    rkinv = Rcpp::match(unak, rkdist(i_o, Rcpp::_));
    for (int k = 0; k<K; k++) {
      emp2i[k] = emp2i[k] + empi[rkinv(k)-1];
      // out_papn(i, k) = dispo[i_o*K+k]; 
      // out_dispo(i_o, k) = dispo[i_o*K+k];
      out_emps(i_o,k) = empi[rkinv(k)-1];
    }
  }
  
  // zou
  return out_emps;
}