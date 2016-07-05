#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double rootChoose(int n, int k, double root){
  double nomin = 1, denom = 1;
  
  // comp nomin
  for(int i = n; i > k; i--){
    nomin *= std::pow(i,1/root);
    }
  for(int i = n-k; i > 0; i--){
    denom *= std::pow(i,1/root);
    }
  return(nomin/denom);
  }


