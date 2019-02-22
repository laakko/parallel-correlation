# parallel-correlation
(*Assignments for a parallel programming course*)

Efficient parallel implementations of correlation-pair calculation.
This project can be used as a reference on how to efficiently calculate matrix multiplication by utilizing parallelism as much as possible.


input: `data` is a pointer to the input matrix, with `ny` rows and `nx` columns
output: correlation coefficients of the given input matrix


## Versions

### cp1
- Basic sequential implementation

### cp2a
- Instruction level parallelism utilized

### cp2c
- Parallelized with vector operations

### cp3a
- Fast solution with instruction level parallelism, vector operations and multithreading with OMP library


## Instructions
Navigate to folder of the desired implementation  
Run png correlations with `./pngcorrelate INPUT_IMG OUTPUT1 OUTPUT2`   
Run tests with `./cp-test` 
Run benchmarks with `./cp-benchmark Y X Iterations` 



## Results
(Best of 3 iterations)  
Benchmarked with input size `ny * nx = 4000 * 1000` on Intel Core i7-7500U:

| Version       | Running time  | 
| ------------- |:-------------:| 
| cp1           | 9.812s        |
| cp2a          | 6.082s        |
| cp2c          | 4.954s        |
| cp3a          | 1.123s        |

  
Benchmarked with input size `ny * nx = 3000 * 4000`(equal to 12 Megapixels image) on Intel Core i7-7500U:

| Version       | Running time  | 
| ------------- |:-------------:| 
| cp1           | 25.884s       |
| cp2a          | 15.969s       |
| cp2c          | 11.735s       |
| cp3a          |  2.902s       |