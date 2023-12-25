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

    auto u = [](const double& x, const double& y) -> double {
        return 1 - (x - 1) * (x - 1) - (0.5 - y) * (0.5 - y);
    };

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

    double F[dims[0]];

    size_t d = 0;
    for (size_t j = 1; j < M; j++) {
        for (size_t i = 1; i < N; i++) {
            F[d] = f(i, j);
            d++;
        }
    }

    std::cout << "linear system:\n";
    for (size_t i = 0; i < dims[0]; i++) {
        for (size_t j = 0; j < dims[1]; j++) {
            std::cout << A[i][j] << ' ';
        }
        std::cout << " | " << F[i] <<'\n';
    }

    double V[dims[0]];

    for (size_t i = 0; i < dims[0]; i++) {
        V[i] = 0.l;
    }

    SeidelMethod(A, F, V, dims[0], maxIter, eps);

    std::cout << "Solution from Seidel:\n";
    for (size_t i = 0; i <dims[0]; i++) {
        std::cout << V[i] << '\n';
    }

    std::cout << "Result vector is" << '\n';
    double X[(N + 1)][(M + 1)];

    for (size_t j = 0; j < M + 1; j++) {
        for (size_t i = 0; i < N + 1; i++) {
            if (i == 0 || i == N) {
                X[i][j] = mu12(j * k);
            } else if (j == 0 || j == M) {
                X[i][j] = mu34(i * h);
            } else {
                X[i][j] = V[(i - 1) + (j - 1) * (N - 1)];
            }

            std::cout << X[i][j] << '\t';
        }
        std::cout << '\n';
    }

    std::cout << "Output diff\n";
    for (size_t j = 0; j < M + 1; j++) {
        for (size_t i = 0; i < N + 1; i++) {
            std::cout << std::abs(X[i][j] - u(i * h, j * k)) << '\t';
        }
        std::cout << '\n';
    }

    return 0;
}