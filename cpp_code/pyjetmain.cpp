#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fenv.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

extern void jetmain(std::vector<double>& ear, size_t ne, std::vector<double>& param,
                    std::vector<double>& photeng, std::vector<double>& photspec);

extern "C" void pyjetmain(double* ear, size_t ne, double* param, double* photeng,
                          double* photspec) {
    std::vector<double> earvec(ear, ear + ne);
    std::vector<double> paramvec(param, param + ne);
    std::vector<double> photengvec(photeng, photeng + ne);
    std::vector<double> photspecvec(photspec, photspec + ne);

    jetmain(earvec, ne, paramvec, photengvec, photspecvec);
}
