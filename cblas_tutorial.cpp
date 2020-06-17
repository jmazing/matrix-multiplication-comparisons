#include <iostream>
#include <cblas.h>
#include <vector>

using std::vector;

void init_arr(double *arr, int size) {
    int i = 0;
    while(i < size) {
        arr[i] = i + 1;
        i++;
    }
}

void init_vec(vector<double> &vec) {
    int i = 0;
    while(i < vec.size()) {
        vec[i] = i + 1;
        i++;
    }
}

void zero_arr(double *arr, int size) {
    int i = 0;
    while(i < size) {
        arr[i] = 0.0;
        i++;
    }
}

void print_arr(double *arr, int size, int num_cols) {
    int i = 0, row = 0, col = 0;
    while(i < size) {
        if(i == num_cols) {
            row++;
            col = 0;
        }
        std::cout << "mat [" << row << ", " << col << "]: " << arr[i] << std::endl;
        col++;
        i++;
    }
}

void print_vec(vector<double> &C, int col_size) {
    int row = 0, col = 0;
    for(auto element: C) {
        if(col == col_size){
            col = 0;
            row++;
        }

        std::cout << "mat[ " << row << ", " << col << "]: " << element << std::endl;

        col++;
    }
}

int main(int argc, char const *argv[])
{
    // Multiply vector by scalar
    double x[] = {1.0, 2.0, 3.0};
    double coeff = 4.323;
    int x_length = 3;
    cblas_dscal(x_length, coeff, &x[0], 1);
    for(double e: x) {
        std::cout << e << std::endl;
    }

    std::cout << "-----------------------------------------------" << std::endl;

    vector<double> xv = {1.0, 2.0, 3.0};
    cblas_dscal(xv.size(), coeff, &xv[0], 1);
    for(double e: xv) {
        std::cout << e << std::endl;
    }

    std::cout << "-----------------------------------------------" << std::endl;

    int m = 2, k = 2, n = 2;
    double alpha = 1.0, beta = 0.0;
    double *A = new double[m*k];
    double *B = new double[k*n];
    double *C = new double[m*n];

    init_arr(A, m*n);
    init_arr(B, m*n);
    zero_arr(C, m*n);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, n, B, n, beta, C, n);

    print_arr(C, m*n, n);

    std::cout << "-----------------------------------------------" << std::endl;

    vector<double> Avec(m*n, 0);
    vector<double> Bvec(m*n, 0);
    vector<double> Cvec(m*n, 0);

    init_vec(Avec);
    init_vec(Bvec);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, &Avec[0], n, &Bvec[0], n, beta, &Cvec[0], n);

    print_vec(Cvec, n);

    return 0;
}
