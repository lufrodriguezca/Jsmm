#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector assemble_mass(NumericVector k, NumericVector ML, NumericVector detB, NumericVector vals, int nele, NumericVector xPos){
    int ind(1);
    for(int e = 0; e<nele; ++e){
        for(int j=0; j<9; ++j){
            if( j == 0 || j == 3 || j==6){
                ind = nele*0 + e;
            };
            if( j == 1 || j == 4 || j==7){
                ind = nele*1 + e;
            };
            if( j == 2 || j == 5 || j==8){
                ind = nele*2 + e;
            };
            vals[xPos[(e +nele*(j) )]-1] = vals[xPos[(e +nele*(j))]-1] + ML[j]*k[ind]*detB[e];
        };
    };
    return vals;
}

