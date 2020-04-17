#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <future>

using namespace std;

typedef vector<vector<int64_t>> MyMatrix;

bool isPowerOfTwo(int64_t number) {
    return (number & (number - 1)) == 0;
}

MyMatrix firstQuarter(MyMatrix matrix) {
    MyMatrix resMatrix;
    size_t size = matrix.size();
    for (int i = 0; i < size / 2; ++i) {
        vector<int64_t> resRow;
        for (int j = 0; j < size / 2; ++j) {
            resRow.push_back(matrix[i][j]);
        }
        resMatrix.push_back(resRow);
    }
    return resMatrix;
}

MyMatrix secondQuarter(MyMatrix matrix) {
    MyMatrix resMatrix;
    size_t size = matrix.size();
    for (int i = 0; i < size / 2; ++i) {
        vector<int64_t> resRow;
        for (int j = size / 2; j < size; ++j) {
            resRow.push_back(matrix[i][j]);
        }
        resMatrix.push_back(resRow);
    }
    return resMatrix;
}

MyMatrix thirdQuarter(MyMatrix matrix) {
    MyMatrix resMatrix;
    size_t size = matrix.size();
    for (int i = size / 2; i < size; ++i) {
        vector<int64_t> resRow;
        for (int j = 0; j < size / 2; ++j) {
            resRow.push_back(matrix[i][j]);
        }
        resMatrix.push_back(resRow);
    }
    return resMatrix;
}

MyMatrix fourthQuarter(MyMatrix matrix) {
    MyMatrix resMatrix;
    size_t size = matrix.size();
    for (int i = size / 2; i < size; ++i) {
        vector<int64_t> resRow;
        for (int j = size / 2; j < size; ++j) {
            resRow.push_back(matrix[i][j]);
        }
        resMatrix.push_back(resRow);
    }
    return resMatrix;
}

MyMatrix joinQuarters(const MyMatrix &first, const MyMatrix &second, const MyMatrix &third, const MyMatrix &fourth) {
    int64_t size = first.size() * 2;
    MyMatrix m;
    for (int i = 0; i < size; ++i) {
        vector<int64_t> row;
        for (int j = 0; j < size; ++j) {
            if (i < size / 2 && j < size / 2)row.emplace_back(first[i][j]);
            if (i < size / 2 && j >= size / 2)row.emplace_back(second[i][j - size / 2]);
            if (i >= size / 2 && j < size / 2)row.emplace_back(third[i - size / 2][j]);
            if (i >= size / 2 && j >= size / 2)row.emplace_back(fourth[i - size / 2][j - size / 2]);
        }
        m.emplace_back(row);
    }
    return m;

}

MyMatrix initializeSquareMatrix(int64_t size) {
    MyMatrix m;
    for (int i = 0; i < size; i++) {
        vector<int64_t> row;
        row.reserve(size);
        for (int j = 0; j < size; j++) {
            row.emplace_back(0);
        }
        m.emplace_back(row);
    }
    return m;
}

MyMatrix initializeSquareMatrixRandom(int64_t size) {
    MyMatrix m;
    std::random_device rd;
    std::mt19937 g(rd());
    auto distribution = std::uniform_int_distribution<int64_t>(0, 10);

    for (int i = 0; i < size; i++) {
        vector<int64_t> row;
        row.reserve(size);
        for (int j = 0; j < size; j++) {
            row.emplace_back(distribution(g));
        }
        m.emplace_back(row);

    }
    return m;
}

MyMatrix reformMatrix(MyMatrix A, int64_t newSize) {
    int64_t oldSize = A.size();
    MyMatrix m = initializeSquareMatrix(newSize);
    for (int i = 0; i < newSize; ++i) {
        for (int j = 0; j < newSize; ++j) {
            if (i < oldSize && j < oldSize)m[i][j] = A[i][j];
        }
    }
    return m;
}

MyMatrix strassen(int n, MyMatrix A, MyMatrix B, MyMatrix C);

void input(int n, MyMatrix p);

void output(int n, MyMatrix C);

void strassenThreadRoutine(int n, MyMatrix A, MyMatrix B, MyMatrix C, std::promise<MyMatrix> res) {
    res.set_value(strassen(n, A, B, C));
}

int main() {
    size_t matrixSize = 100;
    MyMatrix A = initializeSquareMatrixRandom(matrixSize), B = initializeSquareMatrixRandom(
            matrixSize), C;


    if (!isPowerOfTwo(matrixSize)) {
        matrixSize = pow(2, floor((log2(matrixSize))) + 1);
        A = reformMatrix(A, matrixSize);
        B = reformMatrix(B, matrixSize);
    }


    C = initializeSquareMatrix(matrixSize);

    auto start = std::chrono::steady_clock::now();
    std::promise<MyMatrix> res1Promise, res2Promise, res3Promise, res4Promise;

    std::future<MyMatrix> job1 = res1Promise.get_future();
    std::future<MyMatrix> job2 = res2Promise.get_future();
    std::future<MyMatrix> job3 = res3Promise.get_future();
    std::future<MyMatrix> job4 = res4Promise.get_future();


    std::thread calcThread1(strassenThreadRoutine, matrixSize / 2, firstQuarter(A), firstQuarter(B), C,
                            std::move(res1Promise));
    std::thread calcThread2(strassenThreadRoutine, matrixSize / 2, secondQuarter(A), secondQuarter(B), C,
                            std::move(res2Promise));
    std::thread calcThread3(strassenThreadRoutine, matrixSize / 2, thirdQuarter(A), thirdQuarter(B), C,
                            std::move(res3Promise));
    std::thread calcThread4(strassenThreadRoutine, matrixSize / 2, fourthQuarter(A), fourthQuarter(B), C,
                            std::move(res4Promise));

    //C = strassen(matrixSize, A, B, C);

    MyMatrix res1 = job1.get();
    MyMatrix res2 = job2.get();
    MyMatrix res3 = job3.get();
    MyMatrix res4 = job4.get();

    calcThread1.join();
    calcThread2.join();
    calcThread3.join();
    calcThread4.join();

    C = joinQuarters(res1, res2, res3, res4);

    auto end = std::chrono::steady_clock::now();
    auto durationMs = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);


    cout << "Actual size: " << matrixSize << endl;
    std::cout << "Calculated in " << durationMs.count() / 1000.0 << "s" << std::endl;


    //output(matrixSize, C);

    return 0;
}


void input(int n, MyMatrix p) {
    for (int i = 0; i < n; i++) {
        cout << "Please Input Line " << i + 1 << endl;
        for (int j = 0; j < n; j++) {
            cin >> p[i][j];
        }
    }
}

void output(int n, MyMatrix C) {
    cout << "The Output Matrix is :" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << C[i][j] << " ";
        }
        cout << endl;
    }
}

MyMatrix Matrix_Multiply(MyMatrix A, MyMatrix B) {
    MyMatrix C;
    for (int i = 0; i < 2; i++) {
        vector<int64_t> row;
        for (int j = 0; j < 2; j++) {
            int64_t tempo = 0;
            for (int t = 0; t < 2; t++) {
                tempo += A[i][t] * B[t][j];
            }
            row.emplace_back(tempo);
        }
        C.emplace_back(row);
    }
    return C;
}

MyMatrix Matrix_Add(int n, MyMatrix X, MyMatrix Y) {
    MyMatrix C;
    for (int i = 0; i < n; i++) {
        vector<int64_t> row;
        for (int j = 0; j < n; j++) {
            row.emplace_back(X[i][j] + Y[i][j]);
        }
        C.emplace_back(row);
    }
    return C;
}

MyMatrix Matrix_Sub(int n, MyMatrix X, MyMatrix Y) {
    MyMatrix C;
    for (int i = 0; i < n; i++) {
        vector<int64_t> row;
        for (int j = 0; j < n; j++) {
            row.emplace_back(X[i][j] - Y[i][j]);
        }
        C.emplace_back(row);
    }
    return C;
}

MyMatrix strassen(int n, MyMatrix A, MyMatrix B, MyMatrix C) {
    MyMatrix A11 = initializeSquareMatrix(n), A12 = initializeSquareMatrix(
            n), A21 = initializeSquareMatrix(n), A22 = initializeSquareMatrix(n);
    MyMatrix B11 = initializeSquareMatrix(n), B12 = initializeSquareMatrix(
            n), B21 = initializeSquareMatrix(n), B22 = initializeSquareMatrix(n);
    MyMatrix C11 = initializeSquareMatrix(n), C12 = initializeSquareMatrix(
            n), C21 = initializeSquareMatrix(n), C22 = initializeSquareMatrix(n);
    MyMatrix M1 = initializeSquareMatrix(n), M2 = initializeSquareMatrix(n), M3 = initializeSquareMatrix(
            n), M4 = initializeSquareMatrix(n), M5 = initializeSquareMatrix(n), M6 = initializeSquareMatrix(
            n), M7 = initializeSquareMatrix(n);
    MyMatrix AA = initializeSquareMatrix(n), BB = initializeSquareMatrix(n);

    if (n == 2) {
        C = Matrix_Multiply(A, B);
    } else {

        for (int i = 0; i < n / 2; i++) {
            for (int j = 0; j < n / 2; j++) {
                A11[i][j] = A[i][j];
                A12[i][j] = A[i][j + n / 2];
                A21[i][j] = A[i + n / 2][j];
                A22[i][j] = A[i + n / 2][j + n / 2];

                B11[i][j] = B[i][j];
                B12[i][j] = B[i][j + n / 2];
                B21[i][j] = B[i + n / 2][j];
                B22[i][j] = B[i + n / 2][j + n / 2];
            }
        }

        //Calculate M1 = (A0 + A3) × (B0 + B3)
        AA = Matrix_Add(n / 2, A11, A22);
        BB = Matrix_Add(n / 2, B11, B22);
        M1 = strassen(n / 2, AA, BB, M1);

        //Calculate M2 = (A2 + A3) × B0
        AA = Matrix_Add(n / 2, A21, A22);
        M2 = strassen(n / 2, AA, B11, M2);

        //Calculate M3 = A0 × (B1 - B3)
        BB = Matrix_Sub(n / 2, B12, B22);
        M3 = strassen(n / 2, A11, BB, M3);

        //Calculate M4 = A3 × (B2 - B0)
        BB = Matrix_Sub(n / 2, B21, B11);
        M4 = strassen(n / 2, A22, BB, M4);

        //Calculate M5 = (A0 + A1) × B3
        AA = Matrix_Add(n / 2, A11, A12);
        M5 = strassen(n / 2, AA, B22, M5);

        //Calculate M6 = (A2 - A0) × (B0 + B1)
        AA = Matrix_Sub(n / 2, A21, A11);
        BB = Matrix_Add(n / 2, B11, B12);
        M6 = strassen(n / 2, AA, BB, M6);

        //Calculate M7 = (A1 - A3) × (B2 + B3)
        AA = Matrix_Sub(n / 2, A12, A22);
        BB = Matrix_Add(n / 2, B21, B22);
        M7 = strassen(n / 2, AA, BB, M7);

        //Calculate C0 = M1 + M4 - M5 + M7
        AA = Matrix_Add(n / 2, M1, M4);
        BB = Matrix_Sub(n / 2, M7, M5);
        C11 = Matrix_Add(n / 2, AA, BB);

        //Calculate C1 = M3 + M5
        C12 = Matrix_Add(n / 2, M3, M5);

        //Calculate C2 = M2 + M4
        C21 = Matrix_Add(n / 2, M2, M4);

        //Calculate C3 = M1 - M2 + M3 + M6
        AA = Matrix_Sub(n / 2, M1, M2);
        BB = Matrix_Add(n / 2, M3, M6);
        C22 = Matrix_Add(n / 2, AA, BB);

        //Set the result to C[][N]
        for (int i = 0; i < n / 2; i++) {
            for (int j = 0; j < n / 2; j++) {
                C[i][j] = C11[i][j];
                C[i][j + n / 2] = C12[i][j];
                C[i + n / 2][j] = C21[i][j];
                C[i + n / 2][j + n / 2] = C22[i][j];
            }
        }
    }
    return C;
}