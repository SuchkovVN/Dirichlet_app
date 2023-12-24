#pragma once

#include <cmath>
#include <cstddef>
#include <iostream>
template <class T, class U>
void SeidelMethod(T A, U rhs, U x, const size_t& n, const size_t& maxIter, const double& eps) {
    // A is an object that represents a matrix A[i][j] of size n x n
    // rhs is a right-hand-side vector rhs[i] w/ size n
    // x is a vars vector x[i] w/ size n

    size_t s = 0;
    double curr_eps = eps;
    double diff;
    double prev_x;

    while ((s < maxIter) && (curr_eps >= eps)) {
        for (size_t i = 0; i < n; i++) {
            prev_x = x[i];
            x[i] = rhs[i];
            for (size_t j = 0; j < n; j++) {
                if (j != i) {
                    x[i] = x[i] - A[i][j] * x[j];
                }
            }
            x[i] = x[i] / A[i][i];

            diff = std::abs(x[i] - prev_x);
            curr_eps = diff < curr_eps ? diff : curr_eps;
        }
        s++;
    }

    std::cout << "Seidel Method: " << '\n' << "Iterations: " << s << '\n' << "Result eps: " << curr_eps << '\n';
}