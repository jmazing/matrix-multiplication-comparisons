#include <iostream>
#include <random>
#include <chrono>
#include <vector>

using namespace std::chrono;
using std::vector;

static const int _MATRIX_DIM = 1000;

void init_matrix(vector<vector<double>> &mat) {
    std::default_random_engine generator;
    std::normal_distribution<double> distribution{5.0, 1.0};
    for(int i = 0; i < _MATRIX_DIM; i++) {
        for(int j = 0; j < _MATRIX_DIM; j++) {
            mat[i][j] = distribution(generator);
        }
    }
}

void print_matrix(vector<vector<double>> &mat) {
    for(int i = 0; i < _MATRIX_DIM; i++) {
        for(int j = 0; j < _MATRIX_DIM; j++) {
            std::cout << "mat[ " << i << ", " << j << "]: " << mat[i][j] << std::endl;
        }
    }
}

void naive_multiplication(vector<vector<double>> &result, 
                        vector<vector<double>> &A, 
                        vector<vector<double>> &B) {
    for(int i = 0; i < _MATRIX_DIM; i++) {
        for(int j = 0; j < _MATRIX_DIM; j++) {
            for(int k = 0; k < _MATRIX_DIM; k++) {
                result[i][j] += A[i][k]  * B[k][j];
            }
        }
    }
}

void turn_to_row_major(vector<double> &row_major, 
                    vector<vector<double>> &mat) {
    int row_major_index = 0;
    for(int i = 0; i < _MATRIX_DIM; i++) {
        for(int j = 0; j < _MATRIX_DIM; j++) {
            row_major[row_major_index] = mat[i][j];
            row_major_index++;
        }
    }
}

void turn_to_column_major(vector<double> &column_major, 
                        vector<vector<double>> &mat) {
    int column_major_index = 0;
    for(int i = 0; i < _MATRIX_DIM; i++) {
        for(int j = 0; j < _MATRIX_DIM; j++) {
            column_major[column_major_index] = mat[j][i];
            column_major_index++;
        }
    }
}

void row_column_major_multiplication(vector<double> &result, 
                                    vector<double> &row_major, 
                                    vector<double> &column_major) {
    int c_major_index = 0;
    for(int i = 0; i < _MATRIX_DIM; i++) {
        for(int j = 0; j < _MATRIX_DIM; j++) {
            for(int k = 0; k < _MATRIX_DIM; k++) {
                int row_offset = i * _MATRIX_DIM + k;
                int col_offset = j * _MATRIX_DIM + k;
                result[c_major_index] += row_major[row_offset] * column_major[col_offset];
            }
            c_major_index++;
        }
    }
}


int main(int argc, char const *argv[])
{
    vector<vector<double>> A(_MATRIX_DIM, vector<double> (_MATRIX_DIM));
    vector<vector<double>> B(_MATRIX_DIM, vector<double> (_MATRIX_DIM)); 
    vector<vector<double>> C(_MATRIX_DIM, vector<double> (_MATRIX_DIM)); 
    init_matrix(A);
    init_matrix(B);

    /*
    vector<vector<double>> A{{1, 2}, {3, 4}};
    vector<vector<double>> B{{5, 6}, {7, 8}};
    vector<vector<double>> C(_MATRIX_DIM, vector<double> (_MATRIX_DIM)); 
    */

    auto start = high_resolution_clock::now();
    naive_multiplication(C, A, B);
    auto duration = duration_cast<seconds> (high_resolution_clock::now() - start);

    std::cout << "Milliseconds it took to perform naive multiplication is: " << duration.count() << std::endl;

    vector<double> A_row_major(_MATRIX_DIM * _MATRIX_DIM, 0);
    vector<double> B_column_major(_MATRIX_DIM * _MATRIX_DIM, 0);
    vector<double> C_major(_MATRIX_DIM * _MATRIX_DIM, 0);

    turn_to_row_major(A_row_major, A);

    /*
    std::cout << "A row major is: " << std::endl;
    for(auto element: A_row_major){
        std::cout << element << std::endl;
    }
    */

    turn_to_column_major(B_column_major, B);

    /*
    std::cout << "B column major is: " << std::endl;
    for(auto element: B_column_major){
        std::cout << element << std::endl;
    }
    */

    start = high_resolution_clock::now();
    row_column_major_multiplication(C_major, A_row_major, B_column_major);
    duration = duration_cast<seconds> (high_resolution_clock::now() - start);

    std::cout << "Milliseconds it took to perform row_column_major multiplication is: " << duration.count() << std::endl;

    /*
    for(auto element: C_major){
        std::cout << element << std::endl;
    }
    */

    return 0;
}
