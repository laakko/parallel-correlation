#include "cp.h"
#include <math.h>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <malloc.h>


typedef double double4_t __attribute__ ((vector_size (4 * sizeof(double))));
static double4_t* double4_alloc(std::size_t n) {
    void* tmp = 0;
    if (posix_memalign(&tmp, sizeof(double4_t), sizeof(double4_t) * n)) {
        throw std::bad_alloc();
    }
    return (double4_t*)tmp;
}

constexpr double4_t d4zeros {
    0.0, 0.0, 0.0, 0.0, 
};


inline double sum4(double4_t v) {
    double final_sum = 0.0;
    for(int i=0; i<4; ++i) {
        final_sum += v[i];
    }
    return final_sum;
}



void correlate(int ny, int nx, const float* data, float* result) {


    double sum;
    constexpr int nb = 4;
    int na = (nx + nb - 1) / nb;
    double *X = new double[nx*ny]; 
    double4_t* X2 = double4_alloc(na*ny);

    // Loop input matrix
    for(int y=0; y<ny; ++y) {
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
            X[x+ynx] = *dp++ - sum;
        }

        sum = 0.0;

        // Normalize: for each row the sum of squares is 1
        for(int x=0; x<nx; ++x) {
            double element = X[x+ynx];
            sum += element*element;    
        }
        
        sum = sqrt(sum);
                    
        for(int x=0; x<nx; ++x) {
            X[x+ynx]  /= sum;
            

        }     

        // Copy X (double) into X2 (double4_t)
        for(int ka=0; ka < na; ++ka) {
            for(int kb=0; kb<nb; ++kb) {
                int i = ka*nb + kb;
                X2[y*na + ka][kb] = i < nx ? X[i+y*nx] : 0.0;
                }
            }
        
    }

    delete[] X;


    // Add correlations to the output matrix
    for(int y=0; y<ny; y++) {
        int yna = y*na;
        for(int x=0; x<=y; x++) {    
            int xna = x*na;
            double4_t finalsumV = d4zeros;
            for(int i=0; i<na; ++i) {
                 // 4 multiplications in parallel
                double4_t sum1 = X2[i+xna];
                double4_t sum2 = X2[i+yna];
                double4_t sum3 = sum1*sum2;
                finalsumV += sum3;
                
            }

            result[y+x*ny] = sum4(finalsumV);
        }
    }

    free(X2);  
}