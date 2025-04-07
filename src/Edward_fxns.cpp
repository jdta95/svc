#include <RcppArmadillo.h>
#include "svc_fxns.h"

arma::mat Z_matrix(arma::mat X, arma::vec select)
{
  int number = arma::sum(select);
  
  arma::mat Z(X.n_rows, number);
  
  int col_index = 0;
  for(int i =0; i < X.n_cols; i++)
  {
    if(select(i)==1)
    {
      Z.col(col_index) = X.col(i);
      col_index++;
    }
  }
  
  return Z;
}
