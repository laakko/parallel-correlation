#include "cp.h"
#include <math.h>
#include <cmath>
#include <vector>

void correlate(int ny, int nx, const float* data, float* result) {


    double sum;
    constexpr int nb = 2;
    int na = (nx + nb - 1) / nb;
    int nab = na*nb;
    double *X = new double[nab*ny]; 


    // Loop input matrix
    for(int y=0; y<ny; ++y) {
        
        // Fill X with padding
        for(int i=0; i<nab; i++){
            X[i + y*nab] = 0.0;
        }

        sum = 0.0;
        int ynx = y*nx;
        const float *dp = &data[ynx];
        // Normalize: arithmetic mean of 0
        for(int x=0; x<nx; ++x) {
            sum += *dp++; 
        }

        sum /= nx;
        
        dp = &data[ynx];
        
        for(int x=0; x<nx; ++x) {
          X[x+y*nab] = *dp++ - sum;
        }

        sum = 0.0;

        // Normalize: for each row the sum of squares is 1
        for(int x=0; x<nx; ++x) {
            double element = X[x+y*nab];
            sum += element*element;    
        }
        
        sum = sqrt(sum);
                    
        for(int x=0; x<nx; ++x) {
            X[x+y*nab]  /= sum;
        } 
    
    }
    

    // Add correlations to the output matrix
    for(int y=0; y<ny; y++) {
        int ynab = y*nab; 
        for(int x=0; x<=y; x++) {    
            double sum1 = 0.0, sum2 = 0.0;
            int xnab = x*nab; 
            for(int i=0; i<na; i++) {
                sum1 += X[nb*i + xnab]*X[nb*i + ynab];
                sum2 += X[nb*i+1 + xnab]*X[nb*i+1 + ynab];
            }
            result[y+x*ny] = sum1 + sum2;
        }
    }
    
    delete[] X;   
}

   

