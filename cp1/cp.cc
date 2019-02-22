#include "cp.h"
#include <math.h>
#include <vector>

void correlate(int ny, int nx, const float* data, float* result) {

    double *X = new double[nx*ny]; 
    double sum;

    // Loop input matrix
    for(int y=0; y<ny; ++y) {
        sum = 0.0;

        // Normalize: arithmetic mean of 0
        for(int x=0; x<nx; ++x) {
            // Element at (x,y)
            sum += data[x + y*nx];
        }
        sum /= nx;
        for(int x=0; x<nx; ++x) {
            X[x+y*nx] = data[x+y*nx] - sum;
        }

        sum = 0.0;

        // Normalize: for each row the sum of squares is 1
        for(int x=0; x<nx; ++x) {
            // Element at (x,y)
            double element = X[x+y*nx];
            sum += element*element;    
        }
                    
        sum = sqrt(sum);
                    
        for(int x=0; x<nx; ++x) {
            X[x+y*nx]  /= sum;
        }     
    }
    
    // Add correlations to the output matrix
    for(int y=0; y<ny; y++) {
        for(int x=0; x<=y; x++) {
            sum = 0.0;
            for(int i=0; i<nx; ++i) {
                sum += X[i + x*nx]*X[i+ y*nx];
            }
            result[y+x*ny] = sum;
        }
    }
    
    delete[] X;      
}