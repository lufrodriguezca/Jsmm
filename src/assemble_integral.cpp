#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector assemble_integral(NumericVector k_v, NumericVector detB, int nnod, int ns, 
                                int nele, NumericVector el_v){
    NumericVector int_v(nnod*ns);
    
    for(int sp = 0; sp < ns; ++sp){
        for(int e = 0; e<nele; ++e){
            for(int i=0; i<3; ++i){
                int_v[nnod*sp + el_v[(i*nele + e)]-1] =  int_v[nnod*sp + el_v[(i*nele + e)]-1]  + k_v[ (nele*3)*sp + (i*nele + e)]*detB[e]/6;
            }
        }
    }
    return int_v;
}

