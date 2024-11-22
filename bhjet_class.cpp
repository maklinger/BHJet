#include "bhjet_class.h"

#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;


BhJetClass::BhJetClass() : params(28, 0.0) { // define constructor 

    //internal BHJet Parameters ----- 

    IsShock;   // Jet shock heating flag
    nz;        // Total number of zones
    nel;       // Number of elements for electron distribution
    syn_res;   // Resolution for synchrotron calculations
    com_res;   // Resolution for Compton calculations
    nsyn;      // Number of synchrotron frequency bins
    ncom;      // Number of Compton frequency bins
    npsw;      // Proton calculation switch
    
}


// Read parameters from a file (like .dat, normal)
void BhJetClass::read_params(const std::string& file) {
    std::ifstream inFile(file);
    if (!inFile) {
        throw std::runtime_error("Cannot open parameter file: " + file);
    }

    std::string line;
    int line_nb = 0;
    while (std::getline(inFile, line)) {
        line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](unsigned char c) { return !std::isspace(c); }));
        if (line.empty() || line[0] == '#') {
            continue;
        }
        params[line_nb] = std::atof(line.c_str());
        line_nb++;
    }
    inFile.close();

    update_internal_parameters();
}

// assigning the names and values 
void BhJetClass::update_internal_parameters() {
    Mbh = params[0];
    Eddlum = 1.25e38 * Mbh;
    Rg = params[0] * 1.989e30 / (2.998e10 * 2.998e10); // Example gravitational radius
    theta = params[1];
    dist = params[2];
    redsh = params[3];
    jetrat = params[4] * Eddlum;
    r_0 = params[5] * Rg;
    z_diss = params[6] * Rg;
    z_acc = params[7] * Rg;
    z_max = params[8] * Rg;
    t_e = params[9];
    f_nth = params[10];
    f_pl = params[11];
    pspec = params[12];
    f_heat = params[13];
    f_beta = params[14];
    f_sc = params[15];
    p_beta = params[16];
    sig_acc = params[17];
    l_disk = params[18];
    r_in = params[19] * Rg;
    r_out = params[20] * Rg;
    compar1 = params[21];
    compar2 = params[22];
    compar3 = params[23];
    compsw = params[24];
    velsw = params[25];
    infosw = static_cast<int>(params[26]);
    EBLsw = static_cast<int>(params[27]);
    zmin = 2.0 * Rg;
}

// Set individual parameter
void BhJetClass::set_parameter(int index, double value) {
    if (index < 0 || index >= params.size()) {
        throw std::out_of_range("Parameter index out of range");
    }
    params[index] = value;
    update_internal_parameters();
}


// Set all parameters
void BhJetClass::set_parameters(const std::vector<double>& new_params) {
    if (new_params.size() != params.size()) {
        throw std::invalid_argument("Parameter vector must have 28 elements");
    }
    params = new_params;
    update_internal_parameters();
}

// Get current parameters
std::vector<double> BhJetClass::get_parameters() const {
    return params;
}

void BhJetClass::run() {
    std::cout << "Running BHJet Model with the following parameters:\n";
    std::cout << "Mbh (Black hole mass): " << Mbh << '\n';
    std::cout << "theta (Viewing angle): " << theta << '\n';
    std::cout << "dist (Distance): " << dist << '\n';
    std::cout << "redsh (Redshift): " << redsh << '\n';
    std::cout << "jetrat (Jet power): " << jetrat << '\n';
    std::cout << "r_0 (Initial jet radius): " << r_0 << '\n';
    std::cout << "z_acc (Shock distance): " << z_acc << '\n';
    std::cout << "z_diss (Magnetic acceleration distance): " << z_diss << '\n';
    std::cout << "z_max (Maximum distance): " << z_max << '\n';
    std::cout << "t_e (Electron temperature): " << t_e << '\n';
    std::cout << "f_nth (% Nonthermal particles): " << f_nth << '\n';
    std::cout << "f_pl (Change in PL fraction): " << f_pl << '\n';
    std::cout << "pspec (Nonthermal slope): " << pspec << '\n';
    std::cout << "f_heat (Shock heating): " << f_heat << '\n';
    std::cout << "f_beta (Dynamic time scale): " << f_beta << '\n';
    std::cout << "f_sc (Particle acceleration time scale): " << f_sc << '\n';
    std::cout << "p_beta (Plasma beta): " << p_beta << '\n';
    std::cout << "sig_acc (Acceleration sigma): " << sig_acc << '\n';
    std::cout << "l_disk (Disk luminosity): " << l_disk << '\n';
    std::cout << "r_in (Disk inner radius): " << r_in << '\n';
    std::cout << "r_out (Disk outer radius): " << r_out << '\n';
    std::cout << "compar1: " << compar1 << '\n';
    std::cout << "compar2: " << compar2 << '\n';
    std::cout << "compar3: " << compar3 << '\n';
    std::cout << "compsw (Compton switch): " << compsw << '\n';
    std::cout << "velsw (Velocity switch): " << velsw << '\n';
    std::cout << "infosw (Info switch): " << infosw << '\n';
    std::cout << "EBLsw (EBL switch): " << EBLsw << '\n';
}