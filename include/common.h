#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <complex>
#include <vector>
#include <tuple>
#include <iomanip>
#include <cmath>
#include <random>
#include <numeric>
#include <cstdint>
#include <set>
#include <map>
#include <chrono>
#include <omp.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>
#include <filesystem>
/* #include <sys/stat.h> */

using namespace std;

#define MKL_Complex16 complex<double>
#include <mkl.h>
#define my_int MKL_INT
#define my_cplx complex<double>
#define my_vec vector<double>
#define my_cplx_vec vector<complex<double>>
#define my_mat vector<vector<double>>
#define my_cplx_mat vector<vector<complex<double>>>
#define PI 3.14159265358979323846264338327950288419
#define cplxi std::complex<double>(0.0, 1.0)

#define ACCURATECUT 1.0e-9

void GetParaFromInput_double(const string &fname, const string &string_match, double &para);

void Normalize(const int &dim, double *wf);

void Normalize(const int &dim, complex<double> *wf);

double Vec_Dot(const int &dim, double *x, const int incx, double *y, const int incy);

my_cplx Vec_Dot(const int &dim, my_cplx *x, const int incx, my_cplx *y, const int incy);

void Vec_ax(const int &dim, const my_cplx &a, my_cplx *x, my_cplx *y);

void Vec_ax(const int &dim, const double &a, double *x, double *y);

void Vec_Cross(const double *vec1, const double *vec2, double *result);

void Vec_Cross(const my_cplx *vec1, const my_cplx *vec2, my_cplx *result);

void Mat_Read(const string &filename, vector<vector<double>> &matrix);

void Mat_Read(const string &filename, vector<vector<complex<double>>> &matrix);
//*****************************************


#endif