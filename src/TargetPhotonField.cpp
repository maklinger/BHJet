#include "TargetPhotonField.hpp"


namespace bhjet
{

std::vector<double> TargetPhotonField::get_observed_energy() {return observed_energy;}
std::vector<double> TargetPhotonField::get_observed_energy_flux() {return observed_energy_flux;}
std::vector<double> TargetPhotonField::get_observed_number_flux() {
    std::vector<double> number_flux(observed_energy_flux.size());
    for(size_t i=0; i < observed_energy_flux.size(); i++)
        number_flux[i] = observed_energy_flux[i] / observed_energy[i];
    return number_flux;
}


}