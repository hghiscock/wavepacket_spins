#pragma once
#include <complex>

using namespace std;

void kron(complex<double>*, complex<double>*, complex<double>*, int, int, int);
void cpy_mat(complex<double>*, complex<double>*, int, int, int);
void cpy_mat(complex<double>*, complex<double>*, int, int);
void cpy_mat(complex<double>*, complex<double>*, int);
void add_mat(complex<double>*, complex<double>*, complex<double>*, double, int);
void add_mat(complex<double>*, complex<double>*, complex<double>*, complex<double>, int);
void print_mat(double*, int, int);
void print_mat(complex<double>*, int, int);
template <class T, class U> void print_mat_test(T, U);
void matmul(complex<double>*, complex<double>*, complex<double>*, int, int, int, int);
void matmul(complex<double>*, complex<double>*, complex<double>*, int);
int foo();
