#ifndef PYBHJET_CLASS
#define PYBHJET_CLASS

#include <unordered_map>
#include <vector>
#include <string>

class BhJetClass {
public:
    BhJetClass(); //constructor, initializing class 

    // Read parameters from a file
    void read_params(const std::string& file);

    // Set individual parameter by index
    void set_parameter(int index, double value);

    // Set all parameters
    void set_parameters(const std::vector<double>& new_params);

    // Get current parameters
    std::vector<double> get_parameters() const;

    // Run the model
    void run();

private:
    // internal parameters for bhjet to run, not used in python ----- 
    bool IsShock;
    int nz, nel, syn_res, com_res, nsyn, ncom, npsw;

    // Parameter vector
    std::vector<double> params;

    // internal named variables -- to use in the run function 
    double Mbh, Eddlum, Rg, theta, dist, redsh, jetrat, zmin, r_0, h, z_acc, z_diss, z_max, t_e;
    double f_nth, f_pl, pspec, f_heat, f_beta, f_sc, p_beta, sig_acc, l_disk, r_in, r_out;
    double compar1, compar2, compar3, compsw, velsw;
    int infosw, EBLsw;

    // Update named variables from params
    void update_internal_parameters();

};

#endif