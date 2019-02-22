#include "cp.h"
#include <math.h>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <tuple>
#include <algorithm>
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
    final_sum += v[0];
    final_sum += v[1];
    final_sum += v[2];
    final_sum += v[3];
    return final_sum;
}



void correlate(int ny, int nx, const float* data, float* result) {


    // Width padding
    constexpr int nb = 4;
    int na = (nx + nb - 1) / nb;

    // Height padding
    constexpr int nd = 4;
    int nc = (ny + nd - 1) / nd;
    int ncd = nc*nd;


    double *X = new double[nx*ny]; 
    double4_t* X2 = double4_alloc(na*ncd);

    double sum = 0.0;
    // Loop input matrix
    #pragma omp parallel for private(sum)
    for(int y=0; y<ny; ++y) {
        sum = 0.0;
        const unsigned int ynx = y*nx;
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
            X[x+ynx]  *= (1.0/sum);
            

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

    // Add padding
    #pragma omp parallel for collapse(3)
    for(int y=ny; y<ncd; ++y){
        for(int ka=0; ka < na; ++ka){
            for(int kb=0; kb < nb; ++kb) {
                X2[y*na + ka][kb] = 0.0;
            }
        }
    }


    std::vector<std::tuple<int,int,int>> rows(nc*nc);

    // Add correlations to the output matrix
    #pragma omp parallel for schedule(dynamic,2)
    for(int y=0; y<nc; ++y) {
        for(int x=0; x<=y; ++x) {   

            int yx = _pdep_u32(y, 0x55555555) | _pdep_u32(x, 0xAAAAAAAA);
            rows[x*nc + y] = std::make_tuple(yx,y,x);
            
            y = std::get<1>(rows[x*nc + y]);
            x = std::get<2>(rows[x*nc + y]);

            const unsigned int xnd = x*nd; 
            const unsigned int ynd = y*nd;
            double4_t finalsumV[nd][nd];
            
            finalsumV[0][0] = d4zeros;
            finalsumV[0][1] = d4zeros;
            finalsumV[0][2] = d4zeros;
            finalsumV[0][3] = d4zeros;
            finalsumV[1][0] = d4zeros;
            finalsumV[1][1] = d4zeros;
            finalsumV[1][2] = d4zeros;
            finalsumV[1][3] = d4zeros;
            finalsumV[2][0] = d4zeros;
            finalsumV[2][1] = d4zeros;
            finalsumV[2][2] = d4zeros; 
            finalsumV[2][3] = d4zeros; 
            finalsumV[3][0] = d4zeros;
            finalsumV[3][1] = d4zeros;
            finalsumV[3][2] = d4zeros;
            finalsumV[3][3] = d4zeros;
        



            /* Same as above looped
            // Initialize finalSumV
            for(int id = 0; id < nd; ++id){
                for(int jd=0; jd < nd; ++jd){
                    finalsumV[id][jd] = d4zeros;
                }
            }*/
    
            for(int i=0; i<na; ++i) {
                 // Do matrix multiplications in parallel & reuse data in registers

                constexpr int PF = 12;
                __builtin_prefetch(&X2[na*(xnd+0)+i +  PF]);
                __builtin_prefetch(&X2[na*(xnd+1)+i +  PF]);
                __builtin_prefetch(&X2[na*(xnd+2)+i +  PF]);
                __builtin_prefetch(&X2[na*(xnd+3)+i +  PF]);
                __builtin_prefetch(&X2[na*(ynd+0)+i +  PF]);
                __builtin_prefetch(&X2[na*(ynd+1)+i +  PF]);
                __builtin_prefetch(&X2[na*(ynd+2)+i +  PF]);
                __builtin_prefetch(&X2[na*(ynd+3)+i +  PF]);

            

                double4_t m1_0 = X2[na*(xnd+0)+i];
                double4_t m1_1 = X2[na*(xnd+1)+i];
                double4_t m1_2 = X2[na*(xnd+2)+i];
                double4_t m1_3 = X2[na*(xnd+3)+i];

                double4_t m2_0 = X2[na*(ynd+0)+i];
                double4_t m2_1 = X2[na*(ynd+1)+i];
                double4_t m2_2 = X2[na*(ynd+2)+i];
                double4_t m2_3 = X2[na*(ynd+3)+i];

                finalsumV[0][0] += m1_0*m2_0;
                finalsumV[0][1] += m1_0*m2_1;
                finalsumV[0][2] += m1_0*m2_2;
                finalsumV[0][3] += m1_0*m2_3;
                finalsumV[1][0] += m1_1*m2_0;
                finalsumV[1][1] += m1_1*m2_1;
                finalsumV[1][2] += m1_1*m2_2;
                finalsumV[1][3] += m1_1*m2_3;
                finalsumV[2][0] += m1_2*m2_0;
                finalsumV[2][1] += m1_2*m2_1;
                finalsumV[2][2] += m1_2*m2_2;
                finalsumV[2][3] += m1_2*m2_3;
                finalsumV[3][0] += m1_3*m2_0;
                finalsumV[3][1] += m1_3*m2_1;
                finalsumV[3][2] += m1_3*m2_2;
                finalsumV[3][3] += m1_3*m2_3;

                
                /* CP2C implementation
                double4_t sum1 = X2[i+xna];
                double4_t sum2 = X2[i+yna];
                double4_t sum3 = sum1*sum2;
                finalsumV += sum3; 
                */
                
            }

            
              for(int id=0;id<nd; ++id){
                for(int jd=0;jd<nd; ++jd){
                    int i = ynd + id;
                    int j = xnd + jd;
                    if(i<ny && j <= i){
                            result[i+j*ny] = sum4(finalsumV[jd][id]);       
                    }
                }   
            }


        }
    }

    std::sort(rows.begin(), rows.end());

    free(X2);  
    

}




