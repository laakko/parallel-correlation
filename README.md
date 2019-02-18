# parallel-correlation
(*Assignments for a parallel programming course*)

Efficient parallel implementations of correlation-pair calculation.
This project can be used as a reference on how to efficiently calculate matrix multiplication by utilizing parallelism as much as possible.


input: `data` is a pointer to the input matrix, with `ny` rows and `nx` columns
output: correlation coefficients of the given input matrix


## Versions

### cp1
- Sequential implementation

### cp2a
- Instruction level parallelism utilized

### cp2c
- Parallelized with vector operations

### cp3a
- Fast solution with instruction level parallelism, vector operations and multithreading with OMP library


## Results
Benchmarked with input size `nx * ny = 4000 * 1000` on Intel Core i7-7500U:

| Version       | Running time  | 
| ------------- |:-------------:| 
| cp1           |               |
| cp2a          |               |
| cp2c          |               |
| cp3a          |               |
