#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fenv.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#else
double omp_get_wtime() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;
}
#endif

namespace fs = std::filesystem;

void read_params(const std::string& path, std::vector<double>& pars);

extern void jetinterp(std::vector<double>& ear, std::vector<double>& energ,
                      std::vector<double>& phot, std::vector<double>& photar, size_t ne,
                      size_t newne);
extern void jetmain(std::vector<double>& ear, size_t ne, std::vector<double>& param,
                    std::vector<double>& photeng, std::vector<double>& photspec);

int main([[maybe_unused]] int argc, char* argv[]) {
    // set input path file relative to executable
    fs::path input_path = argv[0];
    input_path.replace_filename("Input/ip.dat");

    double start = omp_get_wtime();

    size_t npar = 28;
    size_t ne = 201;
    double emin = -10;
    double emax = 10;
    double einc = (emax - emin) / static_cast<double>(ne);

    std::vector<double> ebins(ne, 0.0);
    std::vector<double> param(npar, 0.0);
    std::vector<double> spec(ne - 1, 0.0);
    std::vector<double> dumarr(ne - 1, 0.0);

    for (size_t i = 0; i < ne; i++) {
        ebins[i] = std::pow(10, (emin + static_cast<double>(i) * einc));
    }

    read_params(input_path, param);

    jetmain(ebins, ne - 1, param, spec, dumarr);

    double end = omp_get_wtime();
    std::cout << "Total running time: " << end - start << " seconds\n";

    return EXIT_SUCCESS;
}    // ----------  end of function main  ----------

// Read input file
//
//  @param file		Parameters file
//
//  @return pars         Parameters
//
void read_params(const std::string& path, std::vector<double>& pars) {
    std::ifstream inFile;
    inFile.open(path);
    std::string line;
    size_t line_nb = 0;
    if (!inFile) {
        std::cerr << "Can't open input file: " << path << "\n";
        exit(1);
    }
    while (getline(inFile, line)) {
        // Remove whitespace from the beginning of the line
        line.erase(line.begin(), std::find_if(line.begin(), line.end(),
                                              [](unsigned char c) { return !std::isspace(c); }));
        if (line[0] == '#') {
            continue;
        } else {
            pars[line_nb] = atof(line.c_str());
            line_nb++;
        }
    }
    inFile.close();
}
