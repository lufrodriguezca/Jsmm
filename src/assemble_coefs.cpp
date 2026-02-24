#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector assemble_coefs(NumericVector pqx, NumericVector rqx, NumericVector pqy, NumericVector rqy,
                             NumericVector C, String method, NumericVector XBeta_diffusion,
                             NumericVector XBeta_advection_1, NumericVector XBeta_advection_2,
                             NumericVector XBeta_mortality, int nele, NumericVector node_x, NumericVector node_y, NumericVector el_1){
  NumericVector data(nele*4*3);
  NumericVector XBeta_var (nele*3);

  for(int y=0; y < nele; y++){ //elements
    for(int x=0; x<4; x++){ //a,b1,b2,c
      if(x==0){
        XBeta_var = XBeta_diffusion;
      }
      if(x==1){
        XBeta_var = XBeta_advection_1;
      }
      if(x==2){
        XBeta_var = XBeta_advection_2;
      }
      if(x==3){
        XBeta_var = XBeta_mortality;
      }
      //A
      data[(x + y*4)*3 + 0] = pqy[y]*(XBeta_var[nele*2 + y] - XBeta_var[nele + y]) - (XBeta_var[y] - XBeta_var[nele + y])*rqy[y];
      //B
      data[(x + y*4)*3 + 1] = -( pqx[y]*(XBeta_var[nele*2 + y] - XBeta_var[nele + y]) - (XBeta_var[y] - XBeta_var[nele + y])*rqx[y]);
      //D
      data[(x + y*4)*3 + 2] = -(data[(x + y*4)*3 + 0]*node_x[el_1[y]-1] + data[(x + y*4)*3 + 1]*node_y[el_1[y]-1]+ C[y]*XBeta_var[y]);
    }
  };

  return data;
}

