#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

double ecart_fuite(
    double &x,
    NumericVector &dispo,
    NumericVector &c_abs,
    double &fuite) {
  int n = dispo.length();
  double res = 0.;
  
  for(auto i=0; i<n; i++) res += dispo[i] * log(1. + c_abs[i] * x);
  
  return (res + log(fuite));
}

double derivee_ef(
    double &x,
    NumericVector &dispo,
    NumericVector &c_abs,
    double &fuite) {
  int n = dispo.length();
  double res = 0.;
  
  for(auto i=0; i<n; i++) res += dispo[i] * c_abs[i] / (1. + c_abs[i] * x) ;
  
  return (res);
}

// [[Rcpp::export]]
NumericMatrix meaps_rcpp(
    IntegerMatrix rkdist, 
    NumericVector emplois,
    NumericVector actifs,
    NumericVector odds,
    NumericVector f, 
    IntegerVector shuf)
{
  int N = rkdist.nrow();
  int K = rkdist.ncol();
  int k;
  double tot, p_ref, c_ref, x_n, x_np1, temp, eps;
  NumericMatrix liaisons(N,K);
  NumericVector emp(K);
  IntegerVector rki(K), ishuf(N);
  LogicalVector nna_rki(K);
  bool need_norm;
  
  emp = emplois / sum(emplois) * sum(actifs * (1 - f));
  ishuf = shuf - 1L;
  need_norm = !is_true(all(odds == 1));
  
  for(auto i: ishuf) {
    
    // on vérifie qu'on est pas trop long
    if(i%1000 == 1) {Rcpp::checkUserInterrupt();}

    rki = rkdist(i, _);
    nna_rki = !is_na(rki);
    k = sum(!is_na(rki));
    IntegerVector arrangement(k);
    
    for(auto j=0; j<K; j++) {
      if(nna_rki[j]==TRUE)
        arrangement[rki[j]-1L] = j;
    }
    
    // arrangement = match(seq_len(K), rki);
    // arrangement = arrangement[!is_na(arrangement)];
    // arrangement = arrangement - 1L;
    // k = arrangement.length();
    
    NumericVector dispo = emp[arrangement];
    NumericVector c_abs(k), p_abs(k), passe(k), reste(k);
    
    tot = sum(dispo);
    
    // on teste que la probabilité de rester avant odd est supérieure à 10^-15, ce qui évite des dépassements de capacité
    
    if (abs(log10(f[i]))/tot<15) {
      p_ref = 1  - pow(f[i ], 1 / tot);
      c_ref = p_ref / (1 - p_ref);
      for (int j = 0; j < k; j++) c_abs[j] = odds[arrangement[j]] * c_ref;
      
      // normalisation des odds via Newton
      if (need_norm == TRUE) {
        temp = 0;
        for(int j = 0; j<k; j++) temp += dispo[j] * odds[arrangement[j]];
        x_n = - log(f[i]) / temp / c_ref;
        do {
          x_np1 = x_n - ecart_fuite(x_n, dispo, c_abs, f[i]) / derivee_ef(x_n, dispo, c_abs, f[i]); 
          eps = abs(x_np1 - x_n);
          x_n = x_np1;
        } while (eps > 1e-9);
        c_abs = c_abs * x_np1;
      }
      
      p_abs = c_abs / (1 + c_abs);
      
      passe[0] = pow(1 - p_abs[0], dispo[0]);
      for(auto j = 1; j < k; j++) passe[j] = passe[j-1] * pow(1 - p_abs[j], dispo[j]);
      
      
      reste[0] = 1 - passe[0];
      for(auto j = 1; j < k; j++) {
        reste[j] = passe[j - 1] - passe[j];
      }
      
      reste = reste * actifs[i];
      
      for(auto j = 0; j < k-1 ; j++) {
        // débordement des emplois disponibles
        if (emp[arrangement[j]] < reste[j]) {
          reste[j+1] = reste[j+1] + reste[j] - emp[arrangement[j]];
          reste[j] = emp[arrangement[j]];
        }
        emp[arrangement[j]] -= reste[j];
        liaisons(i, arrangement[j]) = reste[j];
      }
      
      if (emp[arrangement[k-1]] < reste[k-1]) {
        reste[k-1] = emp[arrangement[k-1]];
      } 
      
      emp[arrangement[k-1]] -= reste[k-1];
      liaisons(i, arrangement[k-1]) = reste[k-1];
      
      // for(auto j=0; j<K; j++) {
      //   
      //   if(Rcpp::traits::is_nan<REALSXP>(liaisons(i,j))) {
      //     cout << i << " " << arrangement[j] << " " << dispo[arrangement[j]] << " " << reste[j] << endl;
      //     cout << " " << emp[arrangement[j]];
      //   }
      // }
    }
  }
  
  return liaisons;
}


// [[Rcpp::export]]
List meaps_boot(
    IntegerMatrix &rkdist, 
    NumericVector &emplois,
    NumericVector &actifs,
    NumericVector &odds,
    NumericVector &f, 
    IntegerMatrix &shufs)
{
  int N = rkdist.nrow();
  int K = rkdist.ncol();
  int Ns = shufs.nrow();
  
  IntegerVector shuf(N);
  NumericMatrix res(N,K), acc(N,K), acc2(N,K);
  
  for(auto s = 0; s < Ns; s++ ) {
    res = meaps_rcpp(rkdist, emplois, actifs, odds, f, shufs(s, _));
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
