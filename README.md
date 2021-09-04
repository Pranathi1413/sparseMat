# sparseMat
This is an algorithm for representing sparse matrices in a particular format and multiplying them to get a product sparse matrix.

## Format of the Matrix
Say there is a matrix `A: R * C`. Every non-zero value `A[i][j]` in this matrix is represented as a pair `(position, value)` where `position` is basically `C * i + j`. The matrix is thus represented as a list of such pairs, of all non-zero values.

## Algorithm
Say, we have to do `A X B`. The algorithm involves transposing the second matrix `B` and sorting the list according to the `position`. Then we iterate the 2 lists corresponding to `A` and `B`, multiplying the values only when the row and column indices match, and update the counter accordingly when they don't. For more details, refer to the comments in the file `algo.cpp`.

## Plots
* For plotting the graphs, randomly generated matrices of different sparsity have been used, with dimensions `A: 400 * 500` and `B: 500 * 600`. 
* For the plots, both `A` and `B` are taken with equal sparsity at each point, which is uniformly increased.
* We see that the execution times are reduced to around 3% with sparse matrices with 5% non-zero values, when compared to the normal matrix multiplication with dense format. This reduction in execution time increase with increase in % of non-zero values as expected, at 13.5% of time with 15% sparsity, to 44.4% of time with 30% sparsity and to 91% of time at 40% non-zeroes.
* The second plot depicts the memory ratio in percentage of memory used in sparse format to that used in dense format, taking into account both operand and product matrices. This steadily increases from around 10% with 0.05 sparsity to almost 100% with 0.5 sparsity. 

## How to Run
Command line arguments are in the order: A's #rows, #columns, sparsity of A, B's #rows, #columns, sparsity of B

***Note: By sparsity, we mean the ratio of number of non-zero values to total number of values (R X C)***
``` 
g++ -o test algo.cpp
./test 400 500 0.1 500 600 0.1
```



##### Disclaimer: No sources or references have been used while coming up with this algorithm and its implementation.
