#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix f_rescale_c(NumericVector x, NumericVector y, int npts) {
  int ind, n = x.size();
  NumericVector cumdiffs(n), steps(npts);
  NumericMatrix xyn(npts,2);
  double step, w1, w2, stepi, cumdiffi;
  for(int i = 0; i < n; i++){
    if(i < 1){
      cumdiffs[i] = sqrt(pow(x[i],2) + pow(y[i],2));
      } else {
      cumdiffs[i] = cumdiffs[i-1] + sqrt(pow(x[i] - x[i-1],2) + pow(y[i] - y[i-1],2));
      }
    }
  step = double(cumdiffs[n-1]) / double(npts-1);
  for(double i = 0; i < npts; i++){
    steps[i] = step * i;
    
    }
  for(int i = 0; i < npts; i++){
    ind = 0;
    for(int j = 0; j < n; j++){
      stepi = steps[i];
      cumdiffi = cumdiffs[j];
      if(stepi > cumdiffi) ind++;
      }
    if(i != (npts-1) && i != 0){
      w1 = std::abs(double(steps[i]) - double(cumdiffs[ind-1]));
      w2 = std::abs(double(steps[i]) - double(cumdiffs[ind]));
      xyn(i,0) = double(x[ind-1]) * w2/(w1+w2) + double(x[ind]) * w1/(w1+w2);
      xyn(i,1) = double(y[ind-1]) * w2/(w1+w2) + double(y[ind]) * w1/(w1+w2);        
      }      
    else if(i == 0){
      xyn(i,0) = 0.0;
      xyn(i,1) = 0.0;        
      }
    else {
      xyn(i,0) = double(x[n-1]);
      xyn(i,1) = double(y[n-1]);        
      }
    }
  return xyn;
  }

// [[Rcpp::export]]
NumericMatrix spatialRescale_c(GenericVector trs, int npts){
  int n = trs.size();
  NumericMatrix xyn(n * npts,2);
  for(int i = 0; i < n; i++){
    NumericMatrix tr = trs[i];    
    NumericMatrix rtr = f_rescale_c(tr(_,0), tr(_,1), npts);
    for(int j = 0; j < npts; j++){
      xyn(i * npts + j, 0) = rtr(j,0);
      xyn(i * npts + j, 1) = rtr(j,1);      
      }
    }
  return xyn;
  }


