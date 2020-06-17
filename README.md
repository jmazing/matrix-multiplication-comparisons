A comparison of speed for matrix multiplication using different techniques and libraries

The file mm_comparisons.cpp is the file you want to execute to compare the times for a 1048576 * 1048576 using the following techniques:
- naive multiplication
- ikj multiplication
- bijk multiplication
- bikj multiplication
- eigen matrix multiplication
- cblas matrix multiplication

To compile mm_comparions.cpp:
Without optimizations:
g++ -std=c++11 -I /usr/local/include/eigen3/ mm_comparisons.cpp -o mm_comparisons -lopenblas

With optimizaitons:
g++ -std=c++11 -I /usr/local/include/eigen3/ mm_comparisons.cpp -o mm_comparisons -lopenblas -O3 -DEIGEN_USE_BLAS

To run mm_comparions:
./mm_comparisons
