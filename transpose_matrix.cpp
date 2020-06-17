#include <iostream>
#include <random>
#include <chrono>
#include <vector>

using namespace std::chrono;
using std::vector;

static const int N = 3;
static const int M = 4;

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

        saved_index = i;
        index = i;
        saved_value = contiguous_mat[index];
        while(moved <= N * M) {
            saved_index = index;
            row_index = index / M;
            index = (index % M) * N + row_index;
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
        }
    }
}

int main(int argc, char const *argv[])
{
    vector<double> A {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    transpose_matrix(A);
    for(auto element: A) {
        std::cout << element << " ";
    }
    std::cout << std::endl;
    return 0;
}
