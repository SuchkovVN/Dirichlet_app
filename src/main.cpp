#include "Data.hpp"
#include "Methods.hpp"
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>

int main(int argc, char* argv[]) {
    size_t N = 4, M = 4, maxIter = 100;
    double eps = 1e-6;
    if (argc == 5) {
        N = atoi(argv[1]);
        M = atoi(argv[2]);
        maxIter = atoi(argv[3]);
        eps = atof(argv[4]);
    }

    const double h = 2.l / N, k = 1.l / M, hsq = h * h, ksq = k * k;

    auto mu12 = [](const double& x) -> double {
        double res = (x - 0.5);
        res *= res;
        return -res;
    };

    auto mu34 = [](const double& x) -> double {
        double res = (x - 1);
        res *= res;
        res = 0.75 - res;
        return res;
    };

    auto f = [&](const size_t& i, const size_t& j) -> double {
        double res = -4.0l;
        if (i == 1 || i == N - 1) {
            res += -mu12(j * k) / (hsq);
        }

        if (j == 1 || j == M - 1) {
            res += -mu34(i * h) / (ksq);
        }

        return res;
    };

    size_t dims[2] = { (N - 1) * (M - 1), (N - 1) * (M - 1) };
    std::vector<std::vector<double>> A(dims[0]);

    for (size_t i = 0; i < dims[0]; i++) {
        A[i].resize(dims[1]);
        for (size_t j = 0; j < dims[1]; j++) {
            A[i][j] = 0.l;
        }
    }

    const size_t s = N - 1;
    for (size_t i = 0; i < dims[0]; i++) {
        A[i][i] = -2 * (1 / hsq + 1 / ksq);
        if (i % s != 0) {
            A[i][i - 1] = 1.l / hsq;
        }

        if (i > s - 1) {
            A[i][i - s] = 1.l / ksq;
        }

        if (i < (dims[0] - s)) {
            A[i][i + s] = 1.l / ksq;
        }

        if (i % s != s - 1) {
            A[i][i + 1] = 1.l / hsq;
        }
    }

    for (size_t i = 0; i < dims[0]; i++) {
        for (size_t j = 0; j < dims[1]; j++) {
            std::cout << A[i][j] << ' ';
        }
        std::cout << '\n';
    }

    double V[dims[0]];

    for (size_t i = 0; i < dims[0]; i++) {
        V[i] = 0.l;
    }
    double F[dims[0]];

    size_t i = 1, j = 1;
    for (size_t k = 0; k < dims[0]; k++) {
        F[k] = f(i, j);
        i++;
        if (i % s == 0) {
            i = 1;
            j++;
        }
        std::cout << F[k] << '\n';
    }

    SeidelMethod(A, F, V, dims[0], maxIter, eps);

    for (size_t i = 0; i < dims[0]; i++) {
        std::cout << V[i] << '\n';
    }

    std::cout << "Result vector is" << '\n';
    double X[M * N];
    i = 0;
    j = 0;
    for (size_t k = 0; k < M * N; k++) {
        if (i == 0 || i == N) {
            X[k] = mu12(j * k);
        } else if (j == 0 || j == M) {
            X[k] = mu34(i * h);
        } else {
            X[k] = V[(i - 1) + (j - 1) * (M - 1)];
        }
        i++;
        if (i % N == 0) {
            i = 0;
            j++;
        }
        std::cout << X[k] << '\n';
    }

    return 0;
}