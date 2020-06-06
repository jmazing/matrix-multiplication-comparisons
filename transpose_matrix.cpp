#include <iostream>
#include <random>
#include <chrono>
#include <vector>

using namespace std::chrono;
using std::vector;

static const int N = 3;
static const int M = 5;

void init_matrix(vector<vector<double>> &mat) {
    
    double count = 0.0;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            mat[i][j] = count;
            count++;
        }
    }
}

void print_matrix(vector<vector<double>> &mat) {
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            std::cout << "mat[ " << i << ", " << j << "]: " << mat[i][j] << std::endl;
        }
    }
    std::cout << std::endl;
}

void naive_multiplication(vector<vector<double>> &result, 
                        vector<vector<double>> &A, 
                        vector<vector<double>> &B) {
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            for(int k = 0; k < N; k++) {
                result[i][j] += A[i][k]  * B[k][j];
            }
        }
    }
}

// Assuming the matrix is in square form
void transpose_square_matrix(vector<vector<double>> &mat) {
    int diagonal_index = N + 1;
    for(int i = 0; i < mat.size(); i++) {
        if(i % diagonal_index == 0)
            continue;
        int col = i % N;
        int row = i / N;
        double tmp = mat[row][col];
        mat[row][col] = mat[col][row];
        mat[col][row] = tmp;
    }
}

void transpose_matrix(vector<double> &contiguous_mat) {
    int index = 1;
    int row_index = -1;
    int saved_value = contiguous_mat[1];
    int moved = 0;
    int saved_index = -1;
    
    vector<bool> checked(N * M, false);
    // initialize first and last spot
    checked[0] = true;
    checked[N * M - 1] = true;

    for(int i = 1; i < checked.size(); i++) {
        if(checked[i] == true) {
            std::cout << "i is: " << i << std::endl;
            continue;
        }
        else {
            index = i;
            saved_value = contiguous_mat[index];
        }

        while(moved <= N * M) {
            saved_index = index;
            std::cout << "saved index is: " << saved_index << std::endl;

            row_index = index / M;

            std::cout << "row index is: " << row_index << std::endl;

            index = (index % M) * N + row_index;

            std::cout << "index is: " << index << std::endl;

            if(checked[saved_index] == true) {
                std::cout << "breaking out of while loop" << std::endl;
                moved = 0;
                break;
            }

            int tmp = contiguous_mat[index];
            std::cout << "saved value is: " << saved_value << std::endl;
            contiguous_mat[index] = saved_value;
            saved_value = tmp;
            moved++;
            checked[saved_index] = true;

            std:: cout << "tmp is: " << tmp << std::endl;
            for(auto element: contiguous_mat) {
                std::cout << element << " ";
            }
            std::cout << std::endl;

            for(auto element: checked) {
                std::cout << element << " ";
            }
            std::cout << std::endl;
        }
    }
}

int main(int argc, char const *argv[])
{
    /*
    vector<vector<double>> A(N, vector<double> (N));
    vector<vector<double>> B(N, vector<double> (N)); 
    vector<vector<double>> C(N, vector<double> (N)); 
    init_matrix(A);
    init_matrix(B);

    print_matrix(B);

    transpose_square_matrix(B);

    print_matrix(B);
    */

    vector<double> D {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
    transpose_matrix(D);
    for(auto element: D) {
        std::cout << element << " ";
    }
    std::cout << std::endl;


    /*
    auto start = high_resolution_clock::now();
    naive_multiplication(C, A, B);
    auto duration = duration_cast<seconds> (high_resolution_clock::now() - start);

    std::cout << "Milliseconds it took to perform naive multiplication is: " << duration.count() << std::endl;

    vector<double> A_row_major(N * N, 0);
    vector<double> B_column_major(N * N, 0);
    vector<double> C_major(N * N, 0);


    start = high_resolution_clock::now();

    duration = duration_cast<seconds> (high_resolution_clock::now() - start);

    std::cout << "Milliseconds it took to perform row_column_major multiplication is: " << duration.count() << std::endl;
    */

    return 0;
}
