#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <cblas.h>

using Eigen::MatrixXd;
using namespace std::chrono;
using std::vector;

static const int MATRIX_DIM = 1024;
static const int BLOCK_WIDTH = 16;

void init_matrix(vector<vector<double>> &mat) {
    double count = 1.0;
    for(int i = 0; i < MATRIX_DIM; i++) {
        for(int j = 0; j < MATRIX_DIM; j++) {
            mat[i][j] = count;
            count++;
        }
    }
}

void init_vec(vector<double> &vec) {
    for(int i = 0; i < vec.size(); i++) {
        vec[i] = i + 1;
    }
}

void zero_matrix(vector<vector<double>> &mat) {
    for(auto &row: mat) {
        std::fill(row.begin(), row.end(), 0);
    }

}

void init_EigenMatrix(MatrixXd &m) {
    double count = 1.0;
    for(int i = 0; i < MATRIX_DIM; i++) {
        for(int j = 0; j < MATRIX_DIM; j++) {
            m(i, j) = count;
            count++;
        }
    }
}

void print_matrix(vector<vector<double>> &mat) {
    for(int i = 0; i < MATRIX_DIM; i++) {
        for(int j = 0; j < MATRIX_DIM; j++) {
            std::cout << "mat[ " << i << ", " << j << "]: " << mat[i][j] << std::endl;
        }
    }
    std::cout << std::endl;
}

void print_row_matrix(vector<double> &C) {
    int row = 0, col = 0;
    for(auto element: C) {
        if(col == MATRIX_DIM){
            col = 0;
            row++;
        }

        std::cout << "mat[ " << row << ", " << col << "]: " << element << std::endl;

        col++;
    }
}

void naive_multiplication(vector<vector<double>> &result, vector<vector<double>> &A, vector<vector<double>> &B) {
    // i, j, k multiplication
    for(int i = 0; i < MATRIX_DIM; i++) {
        for(int j = 0; j < MATRIX_DIM; j++) {
            for(int k = 0; k < MATRIX_DIM; k++) {
                result[i][j] += A[i][k]  * B[k][j];
            }
        }
    }
}

void ikj_multiplication(vector<vector<double>> &result, vector<vector<double>> &A, vector<vector<double>> &B) {
    // i, k, j multiplication
    double r;
    for(int i = 0; i < MATRIX_DIM; i++) {
        for(int k = 0; k < MATRIX_DIM; k++) {
            r = A[i][k];
            for(int j = 0; j < MATRIX_DIM; j++) {
                result[i][j] += r  * B[k][j];
            }
        }
    }
}

void turn_to_row_major(vector<double> &row_major, vector<vector<double>> &mat) {
    for(int i = 0; i < MATRIX_DIM; i++) {
        for(int j = 0; j < MATRIX_DIM; j++) {
            row_major[i * MATRIX_DIM + j] = mat[i][j];
        }
    }
}

void cache_ijk_multiplication(vector<vector<double>> &result, vector<vector<double>> &A, vector<vector<double>> &B) {
    // move BLOCK_WIDTH*BLOCK_WIDTH block
    int i, j, k;
    for(i = 0; i < MATRIX_DIM; i += BLOCK_WIDTH) {
        for(j = 0; j < MATRIX_DIM; j += BLOCK_WIDTH) {
            for (k = 0; k < MATRIX_DIM; k += BLOCK_WIDTH) {
                // Get the indcies within a block
                for(int ib = i; ib < i + BLOCK_WIDTH; ib++) {
                    for(int jb = j; jb < j + BLOCK_WIDTH; jb++) {
                        for(int kb = k; kb < k + BLOCK_WIDTH; kb++) {
                            result[ib][jb] += A[ib][kb] * B[kb][jb];
                        }
                    }
                }
            }
        }
    }
}

void bijk(vector<vector<double>> &result, vector<vector<double>> &A, vector<vector<double>> &B) {
    int i, j, k, kk, jj;
    double sum;
    int en = BLOCK_WIDTH * (MATRIX_DIM/BLOCK_WIDTH);

    for(kk = 0; kk < en; kk += BLOCK_WIDTH) {
        for(jj = 0; jj < en; jj += BLOCK_WIDTH) {
            for(i = 0; i < MATRIX_DIM; i++) {
                for(j = jj; j < jj + BLOCK_WIDTH; j++) {
                    sum = result[i][j];
                    for(k = kk; k < kk + BLOCK_WIDTH; k++) {
                        sum += A[i][k] * B[k][j];
                    }
                    result[i][j] = sum;
                }
            }
        }
    }
}

void bikj(vector<vector<double>> &result, vector<vector<double>> &A, vector<vector<double>> &B) {
    int i, k, j, kk, jj;
    double r;
    int en = BLOCK_WIDTH * (MATRIX_DIM/BLOCK_WIDTH);

    for(kk = 0; kk < en; kk += BLOCK_WIDTH) {
        for(jj = 0; jj < en; jj += BLOCK_WIDTH) {
            for(i = 0; i < MATRIX_DIM; i++) {
                for(k = kk; k < kk + BLOCK_WIDTH; k++) {
                    r = A[i][k];
                    for(j = jj; j < jj + BLOCK_WIDTH; j++) {
                        result[i][j] += r * B[k][j];
                    }
                }
            }
        }
    }
}

void bijk_row(vector<double> &result, vector<double> &A, vector<double> &B) {
    int i, j, k, kk, jj;
    double sum;
    int en = BLOCK_WIDTH * (MATRIX_DIM/BLOCK_WIDTH);

    for(kk = 0; kk < en; kk += BLOCK_WIDTH) {
        for(jj = 0; jj < en; jj += BLOCK_WIDTH) {
            for(i = 0; i < MATRIX_DIM; i++) {
                for(j = jj; j < jj + BLOCK_WIDTH; j++) {
                    sum = result[i * MATRIX_DIM + j];
                    for(k = kk; k < kk + BLOCK_WIDTH; k++) {
                        sum += A[i * MATRIX_DIM + k] * B[k * MATRIX_DIM + j];
                    }
                    result[i * MATRIX_DIM + j] = sum;
                }
            }
        }
    }
}

void bikj_row(vector<double> &result, vector<double> &A, vector<double> &B) {
    int i, k, j, kk, jj;
    double r;
    int en = BLOCK_WIDTH * (MATRIX_DIM/BLOCK_WIDTH);

    for(kk = 0; kk < MATRIX_DIM; kk += BLOCK_WIDTH) {
        for(jj = 0; jj < MATRIX_DIM; jj += BLOCK_WIDTH) {
            for(i = 0; i < MATRIX_DIM; i++) {
                for(k = kk; k < kk + BLOCK_WIDTH; k++) {
                    r = A[i * MATRIX_DIM + k];
                    for(j = jj; j < jj + BLOCK_WIDTH; j++) {
                        result[i * MATRIX_DIM + j] += r * B[k * MATRIX_DIM + j];
                    }
                }
            }
        }
    }
}

void cache_row_blocking_multiplication(vector<double> &result, vector<double> &A, vector<double> &B) {
    // move BLOCK_WIDTH*BLOCK_WIDTH block
    for(int i = 0; i < MATRIX_DIM; i += BLOCK_WIDTH) {
        for(int j = 0; j < MATRIX_DIM; j += BLOCK_WIDTH) {
            for (int k = 0; k < MATRIX_DIM; k += BLOCK_WIDTH) {
                // Get the indcies within a block
                for(int ib = i; ib < i + BLOCK_WIDTH; ib++) {
                    for(int jb = j; jb < j + BLOCK_WIDTH; jb++) {
                        for(int kb = k; kb < k + BLOCK_WIDTH; kb++) {
                            result[ib * MATRIX_DIM + jb] += A[ib * MATRIX_DIM + kb] * B[kb * MATRIX_DIM + jb];
                        }
                    }
                }
            }
        }
    }
}

int main(int argc, char const *argv[])
{
    vector<vector<double>> A(MATRIX_DIM, vector<double> (MATRIX_DIM));
    vector<vector<double>> B(MATRIX_DIM, vector<double> (MATRIX_DIM)); 
    vector<vector<double>> C(MATRIX_DIM, vector<double> (MATRIX_DIM)); 
    init_matrix(A);
    init_matrix(B);

    // naive matrix multiplication
    auto start = high_resolution_clock::now();
    naive_multiplication(C, A, B);
    auto duration = duration_cast<milliseconds> (high_resolution_clock::now() - start);
    std::cout << "Naive multiplication took " << duration.count() <<  " milliseconds" << std::endl;
    zero_matrix(C);

    //ikj multiplication
    start = high_resolution_clock::now();
    ikj_multiplication(C, A, B);
    duration = duration_cast<milliseconds> (high_resolution_clock::now() - start);
    std::cout << "ikj multiplication took " << duration.count() <<  " milliseconds" << std::endl;
    zero_matrix(C);

    // cache ijk multiplication
    start = high_resolution_clock::now();
    cache_ijk_multiplication(C, A, B);
    duration = duration_cast<milliseconds> (high_resolution_clock::now() - start);
    std::cout << "cached multiplication took " << duration.count() <<  " milliseconds" << std::endl;
    zero_matrix(C);

    // cache bijk multiplication
    start = high_resolution_clock::now();
    bijk(C, A, B);
    duration = duration_cast<milliseconds> (high_resolution_clock::now() - start);
    std::cout << "bijk multiplication took " << duration.count() <<  " milliseconds" << std::endl;
    zero_matrix(C);

    // cache bikj multiplication
    start = high_resolution_clock::now();
    bikj(C, A, B);
    duration = duration_cast<milliseconds> (high_resolution_clock::now() - start);
    std::cout << "bikj multiplication took " << duration.count() <<  " milliseconds" << std::endl;
    zero_matrix(C);

    vector<double> A_row_major(MATRIX_DIM * MATRIX_DIM, 0);
    vector<double> B_row_major(MATRIX_DIM * MATRIX_DIM, 0);
    vector<double> C_row_major(MATRIX_DIM * MATRIX_DIM, 0);

    turn_to_row_major(A_row_major, A);
    turn_to_row_major(B_row_major, B);
    turn_to_row_major(C_row_major, C);

    // cache row blocking
    start = high_resolution_clock::now();
    cache_row_blocking_multiplication(C_row_major, A_row_major, B_row_major);
    duration = duration_cast<milliseconds> (high_resolution_clock::now() - start);
    std::cout << "cache row blocking multiplication took " << duration.count() <<  " milliseconds" << std::endl;
    C.clear();

    // bijk row blocking
    start = high_resolution_clock::now();
    bijk_row(C_row_major, A_row_major, B_row_major);
    duration = duration_cast<milliseconds> (high_resolution_clock::now() - start);
    std::cout << "bijk row blocking multiplication took " << duration.count() <<  " milliseconds" << std::endl;
    C.clear();

    // bikj row blocking
    start = high_resolution_clock::now();
    bikj_row(C_row_major, A_row_major, B_row_major);
    duration = duration_cast<milliseconds> (high_resolution_clock::now() - start);
    std::cout << "bikj row blocking multiplication took " << duration.count() <<  " milliseconds" << std::endl;
    C.clear();

    MatrixXd m(MATRIX_DIM, MATRIX_DIM);
    MatrixXd c(MATRIX_DIM, MATRIX_DIM);
    init_EigenMatrix(m);

    start = high_resolution_clock::now();
    c = m*m;
    duration = duration_cast<milliseconds> (high_resolution_clock::now() - start);
    std::cout << "Eigen Library matrix multiplication took " << duration.count() <<  " milliseconds" << std::endl;

    vector<double> Ablas(MATRIX_DIM*MATRIX_DIM, 0);
    vector<double> Bblas(MATRIX_DIM*MATRIX_DIM, 0);
    vector<double> Cblas(MATRIX_DIM*MATRIX_DIM, 0);

    init_vec(Ablas);
    init_vec(Bblas);

    start = high_resolution_clock::now();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, MATRIX_DIM, MATRIX_DIM, MATRIX_DIM, 1, 
                &Ablas[0], MATRIX_DIM, &Bblas[0], MATRIX_DIM, 0, &Cblas[0], MATRIX_DIM);
    duration = duration_cast<milliseconds> (high_resolution_clock::now() - start);
    std::cout << "cblas_dgemm matrix multiplication took " << duration.count() <<  " milliseconds" << std::endl;
    

}
