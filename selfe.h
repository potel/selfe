#include <iostream>
#include <cstdlib>

#ifdef __linux__
#define LINUX 1
#else
#define LINUX 0
#endif
#ifdef _WIN32
#define WINDOWS 1
#else
#define WINDOWS 0
#endif

#include <iostream>
#include <cstdlib>
#ifdef __linux__
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#else
#ifdef _WIN32
#include <gsl\gsl_multiroots.h>
#include <gsl\gsl_math.h>
#include <gsl\gsl_vector.h>
#include <gsl\gsl_sf.h>
#include <gsl\gsl_complex.h>
#include <gsl\gsl_complex_math.h>
#include <gsl\gsl_deriv.h>
#include <gsl\gsl_integration.h>
#include <gsl\gsl_sf_coupling.h>
#include <gsl\gsl_eigen.h>
#include <gsl\gsl_errno.h>
#endif
#endif
#include <cstdio>
#include <iostream>
#include <complex>
#include <valarray>
#include <vector>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <iomanip>
typedef complex<double> complejo;
typedef vector<double> vector_dbl;
typedef vector<double> vector_cmpx;
#define MAX_PTS 6000
#define MAX_POTS 10
#define MAX_ST 100
#define MAX_GAUSS 200
#define MIN_ENERGIA -1000
#define MAX_ENERGIA  1000
#define MAX_L  200
#define MAX_PARES  1000
#define EPSILON  1.e-7
// Constantes fisicas
#define PI 3.14159265358979323846264338327
#define E_CUADRADO 1.44
//#define AMU 931.49432
#define AMU 938.2720813
#define HC  197.3270533 // MeV*fm
#define HBAR  6.582119e-22 // MeV*s
#define E2HC	0.00729927
#define C 2.99792458e23 // velocidad de la luz en fm/s
#define NEULER 5000
complejo const I(0., 1.);
extern ofstream misc1;
extern ofstream misc2;
extern ofstream misc3;
extern ofstream misc4;
extern ofstream misc5;
extern ofstream misc6;
extern ofstream misc7;
extern ofstream misc8;
extern ofstream informe;


