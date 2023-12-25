#pragma once 

#include <iostream>
#include <vector>
#include <math.h>
#include <iomanip>
#include <stdexcept>

double getDeterminant(const std::vector<std::vector<double>> vect);

std::vector<std::vector<double>> getTranspose(const std::vector<std::vector<double>> matrix1);

std::vector<std::vector<double>> getCofactor(const std::vector<std::vector<double>> vect);

std::vector<std::vector<double>> getInverse(const std::vector<std::vector<double>> vect);

void jacobi(const size_t& n, std::vector<std::vector<double>> a, std::vector<double> d, std::vector<std::vector<double>> v);

void jacobi_eigenvalue(int n, double a[], int it_max, double v[], double d[], int& it_num, int& rot_num);
void r8mat_diag_get_vector(int n, double a[], double v[]);
void r8mat_identity(int n, double a[]);
double r8mat_is_eigen_right(int n, int k, double a[], double x[], double lambda[]);
double r8mat_norm_fro(int m, int n, double a[]);
void r8mat_print(int m, int n, double a[], std::string title);
void r8mat_print_some(int m, int n, double a[], int ilo, int jlo, int ihi, int jhi, std::string title);
void r8vec_print(int n, double a[], std::string title);
void timestamp();