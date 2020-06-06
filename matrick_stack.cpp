#include <iostream>
#include <random>
#include <chrono>

using namespace std::chrono;

static const int _MATRIX_DIM = 400;

void init_matrix(double (&mat)[_MATRIX_DIM][_MATRIX_DIM]) {
    std::default_random_engine generator;
    std::normal_distribution<double> distribution{5.0, 1.0};
    
    for(int i = 0; i < _MATRIX_DIM; i++) {
        for(int j = 0; j < _MATRIX_DIM; j++) {
            mat[i][j] = distribution(generator);
        }
    }
}

void print_matrix(double (&mat)[_MATRIX_DIM][_MATRIX_DIM]) {
    for(int i = 0; i < _MATRIX_DIM; i++) {
        for(int j = 0; j < _MATRIX_DIM; j++) {
            std::cout << "mat[ " << i << ", " << j << "]: " << mat[i][j] << std::endl;
        }
    }
}

void naive_multiplication(double (&result)[_MATRIX_DIM][_MATRIX_DIM], 
                        double (&A)[_MATRIX_DIM][_MATRIX_DIM], 
                        double (&B)[_MATRIX_DIM][_MATRIX_DIM]) {
    for(int i = 0; i < _MATRIX_DIM; i++) {
        for(int j = 0; j < _MATRIX_DIM; j++) {
            for(int k = 0; k < _MATRIX_DIM; k++) {
                result[i][j] += A[i][k]  * B[k][j];
            }
        }
    }
}

void turn_to_row_major(double (&row_major)[_MATRIX_DIM * _MATRIX_DIM], 
                    double (&A)[_MATRIX_DIM][_MATRIX_DIM]) {
    int row_major_index = 0;
    for(int i = 0; i < _MATRIX_DIM; i++) {
        for(int j = 0; j < _MATRIX_DIM; j++) {
            row_major[row_major_index] = A[i][j];
            row_major_index++;
        }
    }
}

void turn_to_column_major(double (&column_major)[_MATRIX_DIM * _MATRIX_DIM], 
                    double (&B)[_MATRIX_DIM][_MATRIX_DIM]) {
    int column_major_index = 0;
    for(int i = 0; i < _MATRIX_DIM; i++) {
        for(int j = 0; j < _MATRIX_DIM; j++) {
            column_major[column_major_index] = B[j][i];
            column_major_index++;
        }
    }
}

void row_column_major_multiplication(double (&C_major)[_MATRIX_DIM * _MATRIX_DIM], 
                        double (&A_row_major)[_MATRIX_DIM * _MATRIX_DIM], 
                        double (&B_column_major)[_MATRIX_DIM * _MATRIX_DIM]) {
    int c_major_index = 0;
    for(int i = 0; i < _MATRIX_DIM; i++) {
        for(int j = 0; j < _MATRIX_DIM; j++) {
            for(int k = 0; k < _MATRIX_DIM; k++) {
                int row_offset = i * 3 + k;
                int col_offset = j * 3 + k;
                C_major[c_major_index] = A_row_major[row_offset] * B_column_major[col_offset];
            }
        }
    }
}


int main(int argc, char const *argv[])
{
    double A[_MATRIX_DIM][_MATRIX_DIM];
    double B[_MATRIX_DIM][_MATRIX_DIM];
    double C[_MATRIX_DIM][_MATRIX_DIM];

    init_matrix(A);
    init_matrix(B);

    auto start = high_resolution_clock::now();
    naive_multiplication(C, A, B);
    auto duration = duration_cast<milliseconds> (high_resolution_clock::now() - start);

    std::cout << "Milliseconds it took to perform naive multiplication is: " << duration.count() << std::endl;

    double A_row_major[_MATRIX_DIM * _MATRIX_DIM];
    double B_column_major[_MATRIX_DIM * _MATRIX_DIM];
    double C_major[_MATRIX_DIM * _MATRIX_DIM];

    turn_to_row_major(A_row_major, A);
    turn_to_column_major(B_column_major, B);

    start = high_resolution_clock::now();
    row_column_major_multiplication(C_major, A_row_major, B_column_major);
    duration = duration_cast<milliseconds> (high_resolution_clock::now() - start);

    std::cout << "Milliseconds it took to perform row_column_major multiplication is: " << duration.count() << std::endl;

    return 0;
}
