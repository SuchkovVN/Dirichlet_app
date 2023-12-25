#include "Data.hpp"
#include "Methods.hpp"
#include "linalg/linalg.hpp"
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

    std::cout << "---Parameters\n"
              << "(N, M) = (" << N << ',' << ' ' << M << ")\n"
              << "Max iterations: " << maxIter << '\n'
              << "eps: " << eps << '\n';

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

    std::vector<double> F(dims[0]);

    size_t d = 0;
    for (size_t j = 1; j < M; j++) {
        for (size_t i = 1; i < N; i++) {
            F[d] = f(i, j);
            d++;
        }
    }

    std::cout << "---linear system:\n";
    for (size_t i = 0; i < dims[0]; i++) {
        for (size_t j = 0; j < dims[1]; j++) {
            std::cout << A[i][j] << ' ';
        }
        std::cout << " | " << F[i] << '\n';
    }

    std::vector<double> V(dims[0]);

    for (size_t i = 0; i < dims[0]; i++) {
        V[i] = 0.l;
    }

    SeidelMethod(A, F, V, dims[0], maxIter, eps);

    std::cout << "---Solution from Seidel:\n";
    for (size_t i = 0; i < dims[0]; i++) {
        std::cout << V[i] << '\n';
    }

    double residual = residual_norm(A, V, F, dims[0]);
    std::cout << "Residual 2-morm: " << residual << '\n';

    std::vector<std::vector<double>> vec(dims[0], std::vector<double>(dims[1]));
    std::vector<double> lambdas(dims[0]);
    double* D = lambdas.data();
    std::vector<double> a(dims[0] * dims[1]);
    double* aptr = a.data();
    for (size_t i = 0; i < dims[0]; i++) {
        for (size_t j = 0; j < dims[1]; j++) {
            a[i + j * dims[0]] = A[i][j];
        }
    }

    int it_max = 50, it_num = 0, rot_num = 0;
    std::vector<double> vv(dims[0] * dims[0]);
    double* vvptr = vv.data();
    jacobi_eigenvalue(dims[0], aptr, it_max, vvptr, D, it_num, rot_num);

    double minLambdaAbs = find_abs_min(lambdas, dims[0]);
    std::cout << "invA norm approx: " << 1 / minLambdaAbs << '\n';
    double topest = residual / minLambdaAbs;

    std::cout << "Top estimation: " << topest << '\n';

    std::cout << "---Result vector is" << '\n';
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

    double max_diff = 0.l;
    double diff = 0.l;
    std::cout << "---Output diff\n";
    for (size_t j = 0; j < M + 1; j++) {
        for (size_t i = 0; i < N + 1; i++) {
            diff = std::abs(X[i][j] - u(i * h, j * k));
            std::cout << diff << '\t';

            max_diff = diff > max_diff ? diff : max_diff;
        }
        std::cout << '\n';
    }

    std::cout << "Max diff: " << max_diff << '\n';

    std::cout << "---True solution:\n";
    for (size_t j = 0; j < M + 1; j++) {
        for (size_t i = 0; i < N + 1; i++) {
            std::cout << u(i * h, j * k) << '\t';
        }
        std::cout << '\n';
    }

    return 0;
}