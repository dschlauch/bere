
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rcpp_ccorr(NumericMatrix x) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericMatrix resultMatrix(nrow,nrow);
  for (int i = 0; i < nrow; i++) {
    for (int j = i; j < nrow; j++) {
      double sumproduct = 0;
      for (int k = 0; k < ncol; k++){
        sumproduct+=x(i,k)*x(j,k);
      }
      resultMatrix(i,j) = sumproduct/(nrow-1);
      resultMatrix(j,i) = sumproduct/(nrow-1);
    }
  }
  return resultMatrix;
}