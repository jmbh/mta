#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix distmat(NumericVector id, 
                      NumericVector x,
                      NumericVector y,
                      int n) {
  
  // declare variables
  NumericVector uni;
  int si;
  
  // get number of unique trajectories
  uni = unique(id);
  si = uni.size();
  
  //step 1: reorder 
  NumericMatrix mx(n,si);
  NumericMatrix my(n,si);
  
  for(int t = 0; t<si; t++) { // loop 
    for(int r = 0; r<n; r++) { // fill rows
      mx(r,t) =  x[t*n+r];
      my(r,t) =  y[t*n+r];
    }
  }
  
  //step 2: calc distance matrix
  NumericMatrix dist(si,si);
  
  for(int r=0; r<si; r++) { // loop rows
    
    for(int c=r; c<si; c++) { // loop cols (only half due to symmetry)
      
      double dummyres=0; //dummy value of final matrix entry
      NumericVector dummycalc(n); //dummy vector to collect differences for each 1:n      
      
      for(int v=0; v<n; v++) { // loop 1:n
        dummycalc[v] = sqrt( pow(mx(v,r)-mx(v,c),2) + pow(my(v,r)-my(v,c),2));
      } 
      
      for(int i=0; i<n; i++) { // sum up dummycalc vector
        dummyres += dummycalc[i];
      }
      
      dist(r,c) = dummyres;
      dist(c,r) = dummyres;
      
    }  
    
  }
  
  return dist;
  
}


