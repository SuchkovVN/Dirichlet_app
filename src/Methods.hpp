#pragma once

#include <cmath>
#include <cstddef>
#include <iostream>
template <class T, class U>
void SeidelMethod(const T& A, const U& rhs, U& x, const size_t& n, const size_t& maxIter, const double& eps) {
    // A is an object that represents a matrix A[i][j] of size n x n
    // rhs is a right-hand-side vector rhs[i] w/ size n
    // x is a vars vector x[i] w/ size n

    size_t s = 0;
    double curr_eps = 0.l;
    double iter_eps = eps;
    double diff;
    double prev_x;

    while ((s < maxIter) && (iter_eps >= eps)) {
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
            curr_eps = diff > curr_eps ? diff : curr_eps;
        }
        s++;
        iter_eps = curr_eps;
        curr_eps = 0.l;
    }

    std::cout << "---Seidel Method: " << '\n' << "Iterations: " << s << '\n' << "Result eps (epsN): " << iter_eps << '\n';
}

template <class T, class U>
double residual_norm(const T& A, const U& v, const U& F, const size_t& n) {
    double res = 0.l;
    double r = 0.l;

    for (size_t i = 0; i < n; i++) { 
        for (size_t j = 0; j < n; j++) {
            r += A[i][j] * v[j];  
        }
        r += -F[i];
        r = r*r;
        res += r;
        r = 0.l;
    }

    return std::sqrt(res);
}
template <class T>
double matrix_dnorm(const T& A, const size_t& m, const size_t& n) {
    double res = 0.l;

    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            res += A[i][j] * A[i][j];
        }
    }

    return std::sqrt(res);
}

template <class T>
double find_abs_min(const T& v, const size_t& n) {
    double min = std::abs(v[0]);
    double tmp;
    for (size_t i = 0; i < n; i++) {
        tmp = std::abs(v[i]);
        min = tmp < min ? tmp : min;
    }

    return min;
}