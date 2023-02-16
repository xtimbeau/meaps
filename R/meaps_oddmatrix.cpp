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
  
  for(int i=0; i<n; i++) res += dispo[i] * log(1. + c_abs[i] * x);
  
  return (res + log(fuite));
}

double derivee_ef(
    double &x,
    NumericVector &dispo,
    NumericVector &c_abs,
    double &fuite) {
  int n = dispo.length();
  double res = 0.;
  
  for(int i=0; i<n; i++) res += dispo[i] * c_abs[i] / (1. + c_abs[i] * x) ;
  
  return (res);
}

// [[Rcpp::export]]
NumericMatrix meaps_oddmat(
    IntegerMatrix rkdist, 
    NumericVector emplois,
    NumericVector actifs,
    NumericMatrix modds,
    NumericVector f, 
    IntegerVector shuf)
{
  int N = rkdist.nrow();
  int K = rkdist.ncol();
  int k;
  double tot, p_ref, c_ref, x_n, x_np1, temp, eps;
  NumericMatrix liaisons(N,K);
  NumericVector emp(K), odds(K);
  IntegerVector rki(K), ishuf(N);
  LogicalVector nna_rki(K);
  bool need_norm;
  
  emp = emplois / sum(emplois) * sum(actifs * (1 - f));
  ishuf = shuf - 1L;
  
  for(int i: ishuf) {
    
    // on vérifie si on est pas trop long
    if(i%1000 == 1) {Rcpp::checkUserInterrupt();}
    
    rki = rkdist(i, _);
    nna_rki = !is_na(rki);
    k = sum(!is_na(rki));
    IntegerVector arrangement(k);
    
    for(auto j=0; j<K; j++) {
      if(nna_rki[j]==TRUE)
        arrangement[rki[j]-1L] = j;
    }
    
    NumericVector dispo = emp[arrangement];
    NumericVector c_abs(k), p_abs(k), passe(k), reste(k);
    
    tot = sum(dispo);
    
    p_ref = 1  - pow(f[i ], 1 / tot);
    if (p_ref == 1)
      continue;
    c_ref = p_ref / (1 - p_ref);
    odds = modds(i, _);
    for (int j = 0; j < k; j++) c_abs[j] = odds[arrangement[j]] * c_ref;
    
    // normalisation des odds via Newton
    need_norm = !is_true(all(odds == 1));
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
    
    NumericVector infinite_abs(k);
    infinite_abs = is_infinite(c_abs);
    p_abs = c_abs / (1 + c_abs);
    for (auto j=0; j<k; j++) {
      if(infinite_abs[j]==TRUE)
        p_abs[j] = 1;
    }
    passe[0] = pow(1 - p_abs[0], dispo[0]);
    for(int j = 1; j < k; j++) passe[j] = passe[j-1] * pow(1 - p_abs[j], dispo[j]);
    
    reste[0] = 1 - passe[0];
    for(int j = 1; j < k; j++) {
      reste[j] = passe[j - 1] - passe[j];
    }
    
    reste = reste * actifs[i];
    // cas du débordement
    for(int j = 0; j < k-1 ; j++) {
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
  }
  return liaisons;
}


// [[Rcpp::export]]
NumericMatrix meaps_pabs(
    IntegerMatrix rkdist, 
    NumericVector emplois,
    NumericVector actifs,
    NumericMatrix modds,
    NumericVector f, 
    IntegerVector shuf)
{
  int N = rkdist.nrow();
  int K = rkdist.ncol();
  int k;
  double tot, p_ref, c_ref, x_n, x_np1, temp, eps;
  NumericMatrix p_abs(N,K);
  NumericVector emp(K), odds(K);
  IntegerVector rki(K), ishuf(N);
  LogicalVector nna_rki(K);
  bool need_norm;
  
  emp = emplois / sum(emplois) * sum(actifs * (1 - f));
  ishuf = shuf - 1L;
  
  for(int i: ishuf) {
    
    // on vérifie si on est pas trop long
    if(i%1000 == 1) {Rcpp::checkUserInterrupt();}
    
    rki = rkdist(i, _);
    nna_rki = !is_na(rki);
    k = sum(!is_na(rki));
    IntegerVector arrangement(k);
    
    for(auto j=0; j<K; j++) {
      if(nna_rki[j]==TRUE)
        arrangement[rki[j]-1L] = j;
    }
    
    NumericVector dispo = emp[arrangement];
    NumericVector c_abs(k);
    
    tot = sum(dispo);
    
    odds = modds(i, _);
    
    p_ref = 1  - pow(f[i ], 1 / tot);
    if (p_ref==1) {
  //    Rcout << "pref=1 i:" << i << " tot:" << tot <<"\n";
      continue;
      }
    c_ref = p_ref / (1 - p_ref);
    for (int j = 0; j < k; j++) c_abs[j] = odds[arrangement[j]] * c_ref;
    
    // normalisation des odds via Newton
    need_norm = !is_true(all(odds == 1));
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
    for(auto j=0; j<K; j++) {
      p_abs(i, j) = NA_REAL;
      if(nna_rki[j]==TRUE) {
        p_abs(i, j) = c_abs[rki[j]-1L] / (1 + c_abs[rki[j]-1L]);
      }
    }
  } // endfor
  return p_abs;
}

// [[Rcpp::export]]
NumericMatrix meaps_oddmat_2(
    IntegerMatrix rkdist, 
    NumericVector emplois,
    NumericVector actifs,
    NumericMatrix modds,
    NumericVector f, 
    IntegerVector shuf)
{
  int N = rkdist.nrow();
  int K = rkdist.ncol();
  int k;
  double tot, p_ref, c_ref, x_n, x_np1, temp, eps;
  NumericMatrix liaisons(N,K);
  NumericVector emp(K), odds(K);
  IntegerVector rki(K), ishuf(N);
  LogicalVector nna_rki(K);
  bool need_norm;
  // int underflow=0, maxu=5;
  
  emp = emplois / sum(emplois) * sum(actifs * (1 - f));
  ishuf = shuf - 1L;
  
  for(int i: ishuf) {
    
    // on vérifie si on est pas trop long
    if(i%1000 == 1) {Rcpp::checkUserInterrupt();}
    
    rki = rkdist(i, _);
    nna_rki = !is_na(rki);
    k = sum(!is_na(rki));
    IntegerVector arrangement(k);
    
    for(auto j=0; j<K; j++) {
      if(nna_rki[j]==TRUE)
        arrangement[rki[j]-1L] = j;
    }
    
    NumericVector dispo = emp[arrangement];
    NumericVector c_abs(k), p_abs(k), passe(k), reste(k);
    
    tot = sum(dispo);
    
    odds = modds(i, _);
    if(tot>0.1) {
      p_ref = 1  - pow(f[i ], 1 / tot);
      c_ref = p_ref / (1 - p_ref);
      if (p_ref==1) {
        c_ref = 1e+15;
        // Rcout << "pref est 1 " << i << " " << tot << " " << f[i] << "\n";
      }
      for (int j = 0; j < k; j++) { 
        c_abs[j] = odds[arrangement[j]] * c_ref;
        // if ((underflow<maxu) & (Rcpp::traits::is_nan<REALSXP>(c_abs[j])) ) {
        //   underflow += 1;
        //   Rcout << "c_abs " << c_abs[j] << " " << i << " " << j << "\n";
        //   Rcout << "  c_ref " << c_ref << "\n";
        //   Rcout << "  odd " << odds[arrangement[j]] << "\n";
        //   Rcout << "  p_ref " << p_ref << "\n";
        //   Rcout << "  f " << f[i] << "\n";
        //   Rcout << "  tot " << p_ref << "\n";
        // }
      }
      // normalisation des odds via Newton
      need_norm = !is_true(all(odds == 1));
      if (need_norm == TRUE) {
        temp = 0;
        for(int j = 0; j<k; j++) temp += dispo[j] * odds[arrangement[j]];
        
        x_n = - log(f[i]) / temp / c_ref;
        
        do {
          x_np1 = x_n - 
            ecart_fuite(x_n, dispo, c_abs, f[i]) / derivee_ef(x_n, dispo, c_abs, f[i]); 
          // if ( (underflow<maxu) & (Rcpp::traits::is_nan<REALSXP>(x_np1)) ) {
          //   Rcout << "x_np1 " << sum((c_abs)) << " " << i << "\n";
          //   Rcout << "  p_ref " << p_ref << "\n";
          //   underflow += 1;
          //   Rcout << "  x_np1 " << x_np1 << "\n";
          //   Rcout << "  x_n " << x_n << "\n";
          //   Rcout << "  dispo " << sum(dispo) << "\n";
          //   Rcout << "  f " << f[i] << "\n";
          // }
          eps = abs(x_np1 - x_n);
          x_n = x_np1;
          if(x_np1>1e+7) {
            x_np1 = 1e+7;
            eps = 0;
            // Rcout << "newton déconne " << i << "\n";
          }
        } while (eps > 1e-9);
        
        c_abs = c_abs * x_np1;
        for(auto j=0; j<k; j++) {
          // if ((underflow<maxu) & (Rcpp::traits::is_nan<REALSXP>(c_abs[j])) ) {
          //   underflow += 1;
          //   Rcout << "c_abs norm " << c_abs[j] << " " << i << " " << j << "\n";
          //   Rcout << "  c_ref " << c_ref << "\n";
          //   Rcout << "  odd " << odds[arrangement[j]] << "\n";
          //   Rcout << "  p_ref " << p_ref << "\n";
          //   Rcout << "  f " << f[i] << "\n";
          //   Rcout << "  tot " << tot << "\n";
          //   Rcout << "  x_np1 " << x_np1 << "\n";
          // }
        }
      }
      
      p_abs = c_abs / (1 + c_abs);
      // p_abs[is_infinite(c_abs)] = 1;
      passe[0] = pow(1 - p_abs[0], dispo[0]);
      
      for(int j = 1; j < k; j++) {
        passe[j] = passe[j-1] * pow(1 - p_abs[j], dispo[j]);
        // if ( (underflow<maxu) & (Rcpp::traits::is_nan<REALSXP>(passe[j])) ) {
        //   Rcout << "passe " << passe[j] << " " << i << " " << j << "\n";
        //   Rcout << "  p_ref " << p_ref << "\n";
        //   underflow += 1;
        //   Rcout << "  dispo " << dispo[j] << "\n";
        //   Rcout << "  c_abs " << c_abs[j] << "\n";
        //   Rcout << "  p_abs " << p_abs[j] << "\n\n";
        // }    
      }
      
      reste[0] = 1 - passe[0];
      for(int j = 1; j < k; j++) {
        reste[j] = passe[j - 1] - passe[j];
      }
      
      reste = reste * actifs[i];
      // cas du débordement
      for(int j = 0; j < k-1 ; j++) {
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
    }
  }
  return liaisons;
}

// [[Rcpp::export]]
List meaps_boot2(
    IntegerMatrix &rkdist, 
    NumericVector &emplois,
    NumericVector &actifs,
    NumericMatrix &modds,
    NumericVector &f, 
    IntegerMatrix &shufs)
{
  int N = rkdist.nrow();
  int K = rkdist.ncol();
  int Ns = shufs.nrow();
  
  IntegerVector shuf(N);
  NumericMatrix res(N,K), acc(N,K), acc2(N,K);
  
  for(int s = 0; s < Ns; s++ ) {
    res = meaps_oddmat(rkdist, emplois, actifs, modds, f, shufs(s, _));
    for(int ii=0; ii<N; ii++) {
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