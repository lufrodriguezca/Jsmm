#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector assemble_stiffness(NumericVector k, NumericVector ML, NumericVector detB, NumericVector vals, 
                                int nele, NumericVector xPos, NumericVector coef, NumericVector Bj, NumericVector C,
                                NumericVector node_x, NumericVector node_y, NumericVector el_1, 
                                NumericVector pre11_01, NumericVector pre11_02, NumericVector pre11_03, NumericVector pre11_AC,
                                NumericVector pre22_01, NumericVector pre22_02, NumericVector pre22_03, NumericVector pre22_BC,
                                NumericVector preB1_01, NumericVector preB1_02, NumericVector preB1_03,
                                NumericVector preB2_01, NumericVector preB2_02, NumericVector preB2_03,
                                NumericVector Kxpp, NumericVector Kypp, NumericVector Kpp){
    int ind(1);
    double MLA11(1);
    double MLA22(1);
    double MLB1(1);
    double MLB2(1);
    double MLC(1);
    NumericVector c01(4);
    NumericVector c02(4);
    NumericVector c03(4);
    NumericVector A(4);
    NumericVector B(4);
            
    for(int e = 0; e<nele; ++e){
        for(int x=0; x<4; ++x){
            c01[x] =  ((coef[(x + e*4)*3 + 0]*Bj[e*4+0] + coef[(x + e*4)*3 + 1]*Bj[e*4+2])/(-C[e]));   
            c02[x] =  ((coef[(x + e*4)*3 + 0]*Bj[e*4+1] + coef[(x + e*4)*3 + 1]*Bj[e*4+3])/(-C[e]));  
            c03[x] =  (coef[(x + e*4)*3 + 0]*node_x[el_1[e]-1] + coef[(x + e*4)*3 + 1]*node_y[el_1[e]-1] +  coef[(x + e*4)*3 + 2])/(-C[e]);   
            A[x] = coef[(x + e*4)*3 + 0];
            B[x] = coef[(x + e*4)*3 + 1];
        }
        for(int j=0; j<9; ++j){
            MLA11 = c01[0]*pre11_01[9*e + j] + c02[0]*pre11_02[9*e + j] + c03[0]*pre11_03[9*e + j] - (A[0]/C[e])*pre11_AC[9*e + j];
            MLA22 = c01[0]*pre22_01[9*e + j] + c02[0]*pre22_02[9*e + j] + c03[0]*pre22_03[9*e + j] - (B[0]/C[e])*pre22_BC[9*e + j];
            MLB1 = c01[1]*preB1_01[9*e + j] + c02[1]*preB1_02[9*e + j] + c03[1]*preB1_03[9*e + j];
            MLB2 = c01[2]*preB2_01[9*e + j] + c02[2]*preB2_02[9*e + j] + c03[2]*preB2_03[9*e + j];
            MLC = c01[3]*Kxpp[j] + c02[3]*Kypp[j] + c03[3]*Kpp[j];
                    
            if( j == 0 || j == 3 || j==6){
                ind = nele*0 + e;
            };
            if( j == 1 || j == 4 || j==7){
                ind = nele*1 + e;
            };
            if( j == 2 || j == 5 || j==8){
                ind = nele*2 + e;
            };
            vals[xPos[(e +nele*(j) )]-1] = vals[xPos[(e +nele*(j))]-1] + (-MLA11 - MLA22 + MLB1 + MLB2 - MLC)*k[ind]*detB[e];
        }
    }
    return vals;
}

