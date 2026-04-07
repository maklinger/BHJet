
#include "MLJet.hpp"
#include "iostream"
#include "filesystem"
#include "kariba/Thermal.hpp"

namespace karcst = kariba::constants;

namespace massloadedjet
{






double jet_mass_density(double nprotons, double nelectrons)
{
	return nprotons * karcst::pmgm + nelectrons * karcst::emgm;
}

double MassLoadedJet::specific_enthalpy_simple(double ge_av, double gp_av)
{
    return ( Ge * (ge_av - 1.) + Gp * (gp_av - 1.) * karcst::pmgm / (karcst::emgm * eta_e_) ) / 
            (1. + karcst::pmgm / karcst::emgm / eta_e_);
}

/// @brief number density of the jet in cm^-3
/// @param g bulk Lorentz factor of the jet
/// @param r radius of the jet segment in cm
/// @return number density of the jet in cm^-3
double MassLoadedJet::number_density_jet(double g, double r)
{
	double beta_j;
	beta_j = (g <= 5.) ? 
                std::sqrt( g * g - 1.) / g : 
                1. - 0.5 * std::pow(g, -2) ;
	return jet_mass_loading_parameters.lepdens * jet_dyn.gam0 * jet_dyn.beta0 / 
                (g * beta_j) * std::pow(r / jet_dyn.r0, -2);
}


void MassLoadedJet::reinit_grid_arrays_for_mass_loading()
{
    if (verbosity_level > 1)
        std::cout << "MassLoadedJet::reinit_grid_arrays_for_mass_loading" << std::endl;
    
    mu_grid = std::vector<double>(n_zones, 0.0);
    specific_enthalpy_grid = std::vector<double>(n_zones, 0.0);
    sigma_grid = std::vector<double>(n_zones, 0.0);
    mass_density_grid = std::vector<double>(n_zones, 0.0);
    electron_energy_density_grid = std::vector<double>(n_zones, 0.0);
    proton_energy_density_grid = std::vector<double>(n_zones, 0.0);
    electron_average_gamma_grid = std::vector<double>(n_zones, 0.0);
    proton_average_gamma_grid = std::vector<double>(n_zones, 0.0);
    reinit_grid_arrays();
}




/// @brief average gamma Lorentz factor of the particles
/// @param pspec the spectral index of the particles
/// @param gmin minimum Lorentz factor
/// @param gmax maximum Lorentz factor
/// @return average Lorentz factor
/// @exception std::invalid_argument if gmin or gmax are not positive
/// @exception std::invalid_argument if gmax is smaller than gmin
/// @exception std::invalid_argument if pspec exactly 1
double average_gamma(double pspec, double gmin, double gmax)
{
    constexpr double eps = 1e-12;

    if (gmin <= 0.0 || gmax <= 0.0) {
        throw std::invalid_argument("average_gamma: gmin and gmax must be positive.");
    }

    if (gmax < gmin) {
        throw std::invalid_argument("average_gamma: gmax must be >= gmin.");
    }

    if (std::abs(gmax - gmin) < eps) {
        return gmin;
    }

    if (std::abs(pspec - 1.0) < eps) {
        throw std::invalid_argument("average_gamma: pspec must not be 1.");
    }

    if (std::abs(pspec - 2.0) < eps) {
        return std::log(gmax / gmin) / (-1.0 / gmax + 1.0 / gmin);
    }

    return (1.0 - pspec) / (2.0 - pspec) *
           (std::pow(gmax, 2.0 - pspec) - std::pow(gmin, 2.0 - pspec)) /
           (std::pow(gmax, 1.0 - pspec) - std::pow(gmin, 1.0 - pspec));
}


/// @brief calculate the max particle energy at the dissipation region
/// @param[inout] ge_av average gamma of the electrons
/// @param[inout] gp_av average gamma of the protons
/// @param z the jet segment distance along the jet axis
/// @file MLJet.cpp
/// @exception std::runtime_error if the approximation fails @todo check if needs to be changed
void MassLoadedJet::approximate_max_particle_energy(double &ge_av, double &gp_av, double z)
{
	/*I approximately find the max energy of accelerated electrons (and protons).
	For the electrons I neglect any photon field and for protons, I only account for escape*/
	double theta = opening_angle_constant_ / jet_dyn.gamf;

	double Rapprox = jet_dyn.r0 + std::max(z - jet_dyn.h0,0.) * std::tan(theta);  //for max particle energy
	double Bapprox = jet_mass_loading_parameters.Bacc * std::pow(z/jet_dyn.acc, -1);
	//approximation for the max energy of electrons; I neglect any photon field so TODO think how to
	//improve it in case the max energy is constrained by ICS

    /// synchrotron losses timescale constant
	double Tsyn0 = 6. * karcst::pi * karcst::emerg * karcst::emerg / 
                    (karcst::sigtom * Bapprox * Bapprox * karcst::cee);
    /// accelerated electrons timescale constant
    double Tesc = factor_break_electrons * Rapprox / karcst::cee;
	/// acceleration electrons timescale constant
    double Tacc0 = 4.0 / (3. * factor_max_energy_electrons * karcst::charg * 
                    karcst::cee * Bapprox);
	/// constant for the calculation of the maximum gamma
    double D = Tsyn0 / Tesc * Tsyn0 / Tesc + 4. * Tsyn0 / Tacc0;
	/// maximum gamma of electrons/protons
    double gmax = ( - Tsyn0 / Tesc + std::sqrt(D)) / (2. * karcst::emerg);


    ge_av = average_gamma(index_injected_electrons, gamma_electron_min_dissipation_, gmax); 
	/// @todo I assume here that the electrons are not cooled at all but might have to compare with the final Eemax
		
    gp_av = 1.;
	if (gamma_proton_min_dissipation_ >= 1.){
		gmax = Tesc / Tacc0 / (karcst::pmgm * karcst::cee * karcst::cee);
		gp_av = average_gamma(index_injected_protons, gamma_proton_min_dissipation_, gmax);
	}
}





/// @brief calculate the jet parameters for the mass loading jet
/// @file MLJet.cpp
/// @exception std::runtime_error if the sigma magnetization at the 
/// dissipation region is not positive
/// @exception std::runtime_error if the approximation of the 
/// maximum particle energy fails due to average_gamma()
void MassLoadedJet::calc_mass_loading_jet()
{
	double dyn_fac = 0.0;
	double rho0 = 0.0;
	double theta_acc = 0.15/jet_dyn.gamf;
    
	//calculate the initial mass density
	dyn_fac = 2.*karcst::pi*jet_dyn.r0*jet_dyn.r0*jet_dyn.gam0*jet_dyn.beta0*karcst::cee;

    jet_mass_loading_parameters.sig0 = magnetization_jet_base_;
    jet_mass_loading_parameters.lepdens = jet_mass_loading_parameters.Nj / 
        (dyn_fac*(karcst::pmgm/jet_mass_loading_parameters.eta+karcst::emgm) *
        karcst::cee*karcst::cee*(1.+jet_mass_loading_parameters.sig0));
    jet_mass_loading_parameters.protdens = jet_mass_loading_parameters.lepdens / jet_mass_loading_parameters.eta;
    
    /// @todo should I make them methods of the class ?
    rho0 = jet_mass_density(jet_mass_loading_parameters.protdens,jet_mass_loading_parameters.lepdens);

    jet_mass_loading_parameters.h_0_specific_enthalpy = specific_enthalpy_simple(jet_mass_loading_parameters.av_gamma, 1.0);
    jet_mass_loading_parameters.mu_initial = jet_dyn.gam0 * 
        (jet_mass_loading_parameters.sig0+jet_mass_loading_parameters.h_0_specific_enthalpy + 1.); //equal to μ_0
    jet_mass_loading_parameters.sig_acc = jet_mass_loading_parameters.mu_initial / jet_dyn.gamf - 
        (jet_mass_loading_parameters.h_0_specific_enthalpy + 1.);


    if (jet_mass_loading_parameters.sig_acc < 0.0) {
        throw std::runtime_error(
            "MassLoadedJet::calc_mass_loading_jet: computed negative sig_acc = "
            + std::to_string(jet_mass_loading_parameters.sig_acc) + "\n" +
            "\tIncrease the initial magnetisation above: " + std::to_string(jet_mass_loading_parameters.sig0) + "\n" +
            "\t and/or decrease η below: " + std::to_string(jet_mass_loading_parameters.eta) + "\n" +
            "\t and/or decrease kBTe below: " + std::to_string(round((gamma_electron_min_dissipation_-1.)*karcst::emerg*karcst::erg/1.0e3)) + "\n" +
            "\t and/or decrease gacc below: " + std::to_string(jet_dyn.gamf) + 
            "\n\t##########################################################\n" +
            "\tYou might want to keep:\n" + 
            "\t the initial magnetisation between 5 and 30\n" +
            "\t the initial temperature between 20 and 1000 keV\n" +
            "\t the maximum Lorentz factor gacc between 3 and 10\n"
        );
    }

    jet_mass_loading_parameters.hdiss = jet_mass_loading_parameters.h_0_specific_enthalpy;
    jet_mass_loading_parameters.rdiss = jet_dyn.r0 + (z_dissipation_cm - jet_dyn.h0) * std::tan(theta_acc);
    jet_mass_loading_parameters.nacc = number_density_jet(jet_dyn.gamf,jet_mass_loading_parameters.rdiss);
    jet_mass_loading_parameters.rhoacc=jet_mass_density(
        jet_mass_loading_parameters.nacc/jet_mass_loading_parameters.eta,jet_mass_loading_parameters.nacc);//the mass density before mass loading
    jet_mass_loading_parameters.Bacc = 
        std::sqrt(4.*karcst::pi*jet_mass_loading_parameters.sig_acc*jet_mass_loading_parameters.rhoacc*karcst::cee*karcst::cee);

    //I keep the following to compare ge_av_diss with the one after loading and calculate tshift(need it?)

    /// @note the following is mainly a check to see if the mass loading jet is working
    try {
        approximate_max_particle_energy(gamma_electron_average_dissipation_,
            gamma_proton_average_dissipation_,z_dissipation_cm);
    } catch (const std::exception &e) {
        throw std::runtime_error(
            std::string("MassLoadedJet::calc_mass_loading_jet: approximate_max_particle_energy failed: ")
            + e.what()
        );
    }

    double rload_end = jet_dyn.r0 + (100.0 * z_dissipation_cm - jet_dyn.h0) * std::tan(theta_acc);
    double nload_end = number_density_jet(jet_dyn.gamf,rload_end);//before mass loading
    jet_mass_loading_parameters.rhoload_end = jet_mass_density(nload_end/jet_mass_loading_parameters.eta,nload_end);//before mass loading

	jet_mass_loading_parameters.bfield = 
            std::sqrt(4.*karcst::pi*jet_mass_loading_parameters.sig0*rho0*karcst::cee*karcst::cee);


}






/// @brief the main method to compute the jet dynamics for MassLoadedJet
/// @details it will initiali all the structrures and arrays, and then call
/// the calc_mass_loading_jet() method to compute the jet dynamics basic quantities
/// Finally, it will call the calc_zone_properties() method for every jet segment
void MassLoadedJet::compute_jet_dynamics()
{
    if (verbosity_level > 1)
        std::cout << "Computing MassLoadedJet dynamics" << std::endl;

    // clean up existing gsl interpolation
    if (verbosity_level > 2)
        std::cout << "cleaning up spline memory" << std::endl;
    if (spline_speed)
        gsl_spline_free(spline_speed);


    // init basic internal params
    eddington_luminosity = 1.25e38 * mass_bh;
    r_g = karcst::gconst * mass_bh * karcst::msun / (karcst::cee * karcst::cee);
    
    jet_dyn.min = z_jet_launching * r_g;
    jet_dyn.max = z_max_calculation * r_g;
    jet_dyn.h0 = 2. * r_initial * r_g + z_jet_launching * r_g;
    jet_dyn.r0 = r_initial * r_g;
    jet_dyn.acc = z_end_of_acceleration * r_g;
    z_dissipation_cm = z_dissipation * r_g;
    jet_dyn.beta0 = std::sqrt(4. / 3. * (4. / 3. - 1.) /
                         (4. / 3. + 1.)); // set initial jet speed for relativistic fluid, g=4/3
    jet_dyn.gam0 = 1. / 
                std::sqrt(1. - jet_dyn.beta0 * jet_dyn.beta0); // set corresponding lorentz factor
    jet_dyn.gamf = gamma_final;
    jet_dyn.Rg = r_g;

    kariba::Thermal dummy_elec(70); // use 70 bins for integral
    dummy_elec.set_temp_kev(electron_temperature_jet_base);
    dummy_elec.set_p();
    dummy_elec.set_norm(1.);
    dummy_elec.set_ndens();


    // jet_mass_loading_parameters.pbeta = plasma_beta_jet_base;
    jet_mass_loading_parameters.Nj = jet_power_eddington * eddington_luminosity;
   	jet_mass_loading_parameters.sig_acc = sigma_final;
	jet_mass_loading_parameters.av_gamma = dummy_elec.av_gamma();
   	jet_mass_loading_parameters.eta = eta_e_;



    try {
        calc_mass_loading_jet();
    } catch (const std::exception &e) {
        throw std::runtime_error(
            std::string("MassLoadedJet::compute_jet_dynamics: failed to compute mass loading jet: ")
            + e.what()
        );
    }


    // grid building
    // double dlgz = 0.1;
    double zmin = jet_dyn.min;
    double zmax = jet_dyn.max;
    double zdiss = z_dissipation*r_g;

    if (zdiss < zmin * std::pow(10, dlgz * 1.5)) {
        zdiss = zmin * std::pow(10, dlgz / 2.0);
        if (verbosity_level > 1)
            std::cout << "z_dissipation<z_jet_launching, set z_dissipation=z_jet_launching*10^(dlgz/2)" << std::endl;
        
    } else if (zdiss > zmax * std::pow(10, -dlgz * 1.5)) {
        zdiss = zmax * pow(10, -dlgz / 2.0);
        if (verbosity_level > 1)
            std::cout << "z_dissipation>z_max_calculation, set z_dissipation=z_max_calculation*10^(-dlgz/2)" << std::endl;
    }

    n_zones= static_cast<size_t>(std::log10(zmax / zmin) / dlgz);
    size_t N1 = static_cast<size_t>(std::log10(zdiss * std::pow(10, -dlgz / 2.0) / zmin) / dlgz);
    double dlgz1 = std::log10(zdiss * std::pow(10, -dlgz / 2.0) / zmin) / static_cast<double>(N1);
    size_t N2 = n_zones - N1 - 1;
    double dlgz2 = std::log10(zmax / (zdiss * std::pow(10, dlgz / 2.0))) / static_cast<double>(N2);
    

    /// @todo why not to add the following function in the constructor?
    reinit_grid_arrays_for_mass_loading();
    // std::cout<< "Nzones=" << n_zones << ", N1=" << N1 << std::endl;

    for (size_t i = 0; i < n_zones; i++)
    {
        // std::cout<< "zone " << i << std::endl;
        // calculate the size of the next zone (formerly "jetgrid")
        // fills z_min_grid and z_height_grid
        // calc_grid_next_zone(i, cut, zcut);
        if(i<N1) {
            z_min_grid[i] = zmin*pow(10, dlgz1 * i);
            z_center_grid[i] = zmin*pow(10, dlgz1 * (i+0.5));
            z_height_grid[i] = zmin*(pow(10, dlgz1 * (i+1)) - pow(10, dlgz1 * i));
        } else if(i==N1){
            z_min_grid[i] = zdiss*pow(10, -dlgz/2);
            z_center_grid[i] = zdiss;
            z_height_grid[i] = zdiss*(pow(10, dlgz/2) - pow(10, -dlgz/2));
        } else {
            z_min_grid[i] = zdiss*pow(10, dlgz/2 + dlgz2 * (i-N1-1));
            z_center_grid[i] = zdiss*pow(10, dlgz/2 +  dlgz2 * (i-N1-1+0.5));
            z_height_grid[i] = zdiss*(pow(10, dlgz/2 + dlgz2 * (i-N1-1+1)) - pow(10, dlgz/2 + dlgz2 * (i-N1-1)));
        }

        // calculate the parameters of the cell (formerly "bljetpars")
        calc_zone_properties(i);
    }


}








/// @brief number density of the jet in cm^-3
/// @param i index of the zone/jet segment
/// @return number density of the jet in cm^-3
double MassLoadedJet::number_density_jet(size_t i)
{
	return jet_mass_loading_parameters.lepdens * jet_dyn.gam0 * jet_dyn.beta0 / 
            beta_gamma_grid[i] * std::pow(radius_grid[i] / jet_dyn.r0, -2);
}




double internal_energy_density(double g_av, double n, double mass)
{
	//the mass in grams
	return (g_av - 1.) * n * mass * karcst::cee * karcst::cee;
}


/// @brief parameterisation of the sigma magnetization of the mass loading jet
/// from Chatterjee et al 2019
/// @param z the jet segment distance along the jet axis in cm
/// @param zdiss dissipation region in cm
/// @return the magnetization of the jet segment at z 
double sigma_loading(double z, double zdiss)
{
	double x = std::log10(z / zdiss);
    return std::pow(10., 0.621 * std::pow(x, 5) - 
            3.005 * std::pow(x, 4) + 
            4.599 * std::pow(x, 3) - 
            2.502 * std::pow(x, 2) + 
            0.2419 * x + 0.5631);
}

/// @brief  parameterisation of the gamma bulk Lorentz factor of the mass loading jet
/// from Chatterjee et al 2019
/// @param z the jet segment distance along the jet axis in cm
/// @param zdiss dissipation region in cm
/// @return the bulk Lorentz factor of the jet segment at z
double gamma_loading(double z, double zdiss)
{
	/*zdiss: where the mass loading inititates (in cm)*/
	double x = std::log10(z / zdiss);
    return std::pow(10.0,
            -0.2764 * std::pow(x, 6) + 
            1.412 * std::pow(x, 5) - 
            2.207 * std::pow(x, 4) + 
            0.8534 * std::pow(x, 3) + 
            0.2565 * std::pow(x, 2) - 
            0.07496 * x + 0.3941);
}

/// @brief  parameterisation of the specific enthalpy of the mass loading jet
/// from Chatterjee et al 2019
/// @param z the jet segment distance along the jet axis in cm
/// @param zdiss dissipation region in cm
/// @return the specific enthalpy of the jet segment at z
double specific_enthalpy_loading(double z, double zdiss)
{
	double x = std::log10(z / zdiss);
    return std::pow(10.,
			0.4668 * std::pow(x, 5) - 
            1.903 * std::pow(x, 4) + 
            1.1 * std::pow(x, 3) + 
            2.482 * std::pow(x, 2) - 
            1.171 * x - 1.826);
}

/// @brief  parameterisation of the mu of the mass loading jet
/// from Chatterjee et al 2019
/// @param z the jet segment distance along the jet axis in cm
/// @param zdiss dissipation region in cm
/// @return the mu of the jet segment at z
double mu_loading(double z, double zdiss)
{
	/*	We only call this function if the initial conditions lead mu to increase */
	double x = std::log10(z / zdiss);
	return std::pow(10.,
            0.2231 * std::pow(x, 5) - 
            0.7242 * std::pow(x, 4) + 
            0.4546 * std::pow(x, 3) + 
            0.104 * std::pow(x, 2) - 
            0.09267 * x + 1.031);
}



/// @brief the method to calculate the jet parameters for the mass loading jet for every jet segment
/// @param i the index of the jet segment
/// @file MLJet.cpp
/// @details the parametrization comes from Chatterjee et al 2019. 
/// It used to be the loadingjetpars() function
void MassLoadedJet::calc_zone_properties(size_t i){
// (double z,jet_dynpars &dyn,jet_enpars &en,particle_pars &prtcls,zone_pars &zone,double &t,
// 	int infosw,double loading_factor_,double heat,std::string outputConfiguration){
	
    /// jet opening angle
    double theta;
    /// average gamma of electrons
	double ge_av;
    /// average gamma of protons
    double gp_av;

	if (z_center_grid[i] < z_dissipation_cm) { 
        /// @note the pre-loading region

		gamma_grid[i] = std::pow(10.,log10(jet_dyn.gam0) + 
				std::log10(jet_dyn.gamf / jet_dyn.gam0) / 
                std::log10(z_dissipation_cm / jet_dyn.min) * 
                std::log10(z_center_grid[i] / jet_dyn.min));
		beta_grid[i] = (gamma_grid[i] <= 5.) ? 
                std::sqrt(gamma_grid[i] * gamma_grid[i] - 1.) / gamma_grid[i] : 
                1. - 0.5 * std::pow(gamma_grid[i], -2) ;
        beta_gamma_grid[i] = beta_grid[i] * gamma_grid[i];
        mu_grid[i] = jet_mass_loading_parameters.mu_initial;
		
        /// @note I can define the specific enthalpy before the zdiss here
        /// or I can calculate the ge_av and compute it later; but it will be the same
        specific_enthalpy_grid[i] = jet_mass_loading_parameters.h_0_specific_enthalpy;
		sigma_grid[i] = mu_grid[i] / gamma_grid[i] - (specific_enthalpy_grid[i] + 1.);

		theta = opening_angle_constant_ / gamma_grid[i];
		radius_grid[i] = jet_dyn.r0 + std::max(z_center_grid[i] - jet_dyn.h0, 0.) * 
                    std::tan(theta);

        electron_density_grid[i] = number_density_jet(i);
		proton_density_grid[i] = electron_density_grid[i] / jet_mass_loading_parameters.eta;
		mass_density_grid[i] = jet_mass_density(proton_density_grid[i], 
                                            electron_density_grid[i]);
		
        proton_energy_density_grid[i] = 
                    internal_energy_density(1, proton_density_grid[i], karcst::pmgm);
        electron_energy_density_grid[i] = specific_enthalpy_grid[i] * mass_density_grid[i] * 
                    karcst::cee * karcst::cee / Ge;
		
		ge_av = jet_mass_loading_parameters.av_gamma;
		gp_av = 1.;
		
	} else if (z_center_grid[i] >= z_dissipation_cm && z_center_grid[i] 
        <= 100. * z_dissipation_cm){

        /// @todo figure out that eta is changing after that point properly
        if (jet_mass_loading_parameters.eta > 10.0) {
            /// after mass loading I switch to 10
            jet_mass_loading_parameters.eta = 10.0;
        }

		double tanh2 = std::tanh(std::log10(z_center_grid[i] / z_dissipation_cm)) * 
                        std::tanh(std::log10(z_center_grid[i] / z_dissipation_cm));
		sigma_grid[i] = sigma_loading(z_center_grid[i], z_dissipation_cm) * 
                        jet_mass_loading_parameters.sig_acc / 
                        sigma_loading(z_dissipation_cm, z_dissipation_cm);
		
        if (heating_fraction_ > 0.0 && heating_fraction_ < 1.0) { 
            ///a fraction "heating" goes to h
            sigma_grid[i] *= (1. - tanh2 * heating_fraction_);
        }
		gamma_grid[i] = gamma_loading(z_center_grid[i], z_dissipation_cm) * 
                        jet_dyn.gamf / gamma_loading(z_dissipation_cm, z_dissipation_cm);
		specific_enthalpy_grid[i] = std::max(
                {
                    specific_enthalpy_loading(z_center_grid[i], z_dissipation_cm) * 
                            jet_mass_loading_parameters.hdiss / 
                            specific_enthalpy_loading(z_dissipation_cm, z_dissipation_cm),
                    specific_enthalpy_grid[i - 1],
                    jet_mass_loading_parameters.h_0_specific_enthalpy
                });
		if (heating_fraction_ > 0.0 && heating_fraction_ < 1.) { 
            ///a fraction "heating" comes from σ magnetization 
            specific_enthalpy_grid[i] += tanh2*heating_fraction_ * sigma_grid[i]; 
		}
        /// @note If h leads to an increas of mu, I'm trying to fix some extremes
		if (z_center_grid[i] < 9.0 * z_dissipation_cm) {
            /// @note force the first few segements to be constant
			mu_grid[i] = jet_mass_loading_parameters.mu_initial;
			if (gamma_grid[i] * (sigma_grid[i] + specific_enthalpy_grid[i] + 1.) > 
                    jet_mass_loading_parameters.mu_initial)
            { //when h from the profile forces mu to increase in the first segments:
				specific_enthalpy_grid[i] = std::max(
                    {
                        mu_grid[i]/gamma_grid[i] - (sigma_grid[i]+1.),
                        specific_enthalpy_grid[i - 1],
                        // zone.helec,
                        jet_mass_loading_parameters.h_0_specific_enthalpy
                    });
			}
		} else {
			if (gamma_grid[i] * (sigma_grid[i] + specific_enthalpy_grid[i] + 1.) > 
                    jet_mass_loading_parameters.mu_end) {
				//when h from the profile forces mu to increase beyound the 8*zdiss segments:
				mu_grid[i] = mu_loading(z_center_grid[i], z_dissipation_cm) * 
                        jet_mass_loading_parameters.mu_initial / mu_loading(z_dissipation_cm, z_dissipation_cm);
                specific_enthalpy_grid[i] = std::max(
                        {
                            mu_grid[i ] /gamma_grid[i] - (sigma_grid[i] + 1.),
                            specific_enthalpy_grid[i - 1],
                        });
				//in case the h is not from the mu equation, I have to fix gamma:
				if (gamma_grid[i] * (sigma_grid[i] + specific_enthalpy_grid[i] + 1.) > 
                        mu_grid[i]) {
                    gamma_grid[i] = mu_grid[i] / (sigma_grid[i] + specific_enthalpy_grid[i] + 1.);
                }
			} else { 
                /// @note this is the default case if nothing weird happens
				mu_grid[i] = gamma_grid[i] * (sigma_grid[i] + specific_enthalpy_grid[i] + 1.);
			}
			
		}
		theta = opening_angle_constant_ / gamma_grid[i];
		radius_grid[i] = jet_dyn.r0 + std::max(z_center_grid[i] - jet_dyn.h0, 0.) * std::tan(theta);

		/// @note the factor 2 comes from log10(100*zdiss/zdiss)
		mass_density_grid[i] = std::pow(10., std::log10(jet_mass_loading_parameters.rhoacc) + 
                    std::log10(loading_factor_ * jet_mass_loading_parameters.rhoload_end / 
                        jet_mass_loading_parameters.rhoacc) / 
                    2.0 * std::log10(z_center_grid[i] / z_dissipation_cm));
		
        electron_density_grid[i] = mass_density_grid[i]  / 
                    (karcst::pmgm / jet_mass_loading_parameters.eta + karcst::emgm);

        /// @note I assume η=10 for loaded matter
        proton_energy_density_grid[i] = electron_density_grid[i] / jet_mass_loading_parameters.eta; 
		
        /// @note the very first time, helec is from pre-shock region; beyond zdiss it's from the previous segment
		electron_energy_density_grid[i] = specific_enthalpy_grid[i - 1] * 
                    mass_density_grid[i] * karcst::cee * karcst::cee / Ge;
		proton_energy_density_grid[i] = std::max(1.e-50,
                    (specific_enthalpy_grid[i] * mass_density_grid[i] * karcst::cee * karcst::cee - 
                        Ge * electron_energy_density_grid[i]) / Gp);
		ge_av = std::max(
                electron_energy_density_grid[i] / (electron_density_grid[i] * karcst::emerg) + 1., 
                1.
            );
	    gp_av = std::max(
                proton_energy_density_grid[i] / (proton_energy_density_grid[i] * 
                    karcst::pmgm * karcst::cee * karcst::cee) + 1., 
                1.
            );
	    jet_mass_loading_parameters.mu_end = mu_grid[i];
	    jet_mass_loading_parameters.hdiss_end = specific_enthalpy_grid[i];
	} else {
        /// beyond the 100 dissipation region 
        double tanh2 = std::tanh(std::log10(z_center_grid[i] / z_dissipation_cm)) * 
                        std::tanh(std::log10(z_center_grid[i] / z_dissipation_cm));
        //The one degree polynomial fit of CLTM19 between 1e4-1e5rg
        sigma_grid[i] = std::pow(10.0,
                -0.09708 * std::log10(z_center_grid[i] / z_dissipation_cm) - 0.178) * 
                jet_mass_loading_parameters.sig_acc / sigma_loading(z_dissipation_cm, z_dissipation_cm);
		if (heating_fraction_ > 0.0 && heating_fraction_ < 1.) {
            sigma_grid[i] *= (1.0 - tanh2 * heating_fraction_);
            //a fraction "heating" goes to h specific enthalpy
        }
        specific_enthalpy_grid[i] = std::max(
            {
                specific_enthalpy_grid[i - 1],
	            std::pow(10.0, - 0.2451 * std::log10(z_center_grid[i] / z_dissipation_cm) - 0.5763) * 
                jet_mass_loading_parameters.hdiss_end / 
                specific_enthalpy_loading(100.0 * z_dissipation_cm, z_dissipation_cm)
            });
		if (heating_fraction_ > 0 && heating_fraction_ < 1.) {
            specific_enthalpy_grid[i] += tanh2 * heating_fraction_ * sigma_grid[i];
        }
		mu_grid[i] = jet_mass_loading_parameters.mu_end;
		gamma_grid[i] = std::max(
                        1.,
                        mu_grid[i] / (sigma_grid[i] + specific_enthalpy_grid[i] + 1.)
                    );

		theta = opening_angle_constant_ / gamma_grid[i];
		radius_grid[i] = jet_dyn.r0 + std::max(z_center_grid[i] - jet_dyn.h0, 0.) * std::tan(theta);
		/// @note based on CLTM19 the mass density has increased by a factor of loading_factor_
		mass_density_grid[i] = loading_factor_ * jet_mass_loading_parameters.rhoload_end * 
                    std::pow(z_center_grid[i] / (100. * z_dissipation_cm), -2);
		electron_density_grid[i] = mass_density_grid[i] / 
                    (karcst::pmgm / jet_mass_loading_parameters.eta + karcst::emgm);
        /// I assume η=1 for loaded matter @todo I'm not sure this is correct
		proton_energy_density_grid[i] = electron_density_grid[i] / jet_mass_loading_parameters.eta; 

		electron_energy_density_grid[i] = specific_enthalpy_grid[i] * 
                    mass_density_grid[i] * karcst::cee * karcst::cee / Ge;
		proton_energy_density_grid[i] = std::max(
                        1.e-50,
                        (specific_enthalpy_grid[i] * mass_density_grid[i] * 
                            karcst::cee * karcst::cee - Ge * electron_energy_density_grid[i]) / Gp
                    );
		ge_av = std::max(
                electron_energy_density_grid[i] / (electron_density_grid[i] * karcst::emerg) + 1.,
                1.
            );
        gp_av = std::max(
                proton_energy_density_grid[i] / 
                    (proton_energy_density_grid[i] * karcst::pmgm * karcst::cee * karcst::cee) + 1.,
                1.
            );
	}
	magnetic_field_grid[i] = std::sqrt(4. * karcst::pi * mass_density_grid[i] * 
                karcst::cee * karcst::cee * sigma_grid[i] * (1. + specific_enthalpy_grid[i]));
    beta_grid[i] = (gamma_grid[i] <= 5.0) ? 
                std::sqrt(gamma_grid[i] * gamma_grid[i] - 1.) / gamma_grid[i] : 
                1. - 0.5 * std::pow(gamma_grid[i], -2) ;
    beta_gamma_grid[i] = beta_grid[i] * gamma_grid[i];

    electron_average_gamma_grid[i] = ge_av;
    proton_average_gamma_grid[i] = gp_av;

}


void MassLoadedJet::print_MLJet_kinematics(std::ostream& os) const
{
 
    const auto& p = jet_mass_loading_parameters;

    const std::array<std::pair<std::string, double>, 18> params{{
        {"av_gamma",     p.av_gamma},
        // {"pbeta",        p.pbeta},
        {"Nj",           p.Nj},
        {"bfield",       p.bfield},
        {"lepdens",      p.lepdens},
        {"protdens",     p.protdens},
        {"eta",          p.eta},
        {"sig0",         p.sig0},
        {"sig_acc",      p.sig_acc},
        {"nacc",         p.nacc},
        {"h_0",          p.h_0_specific_enthalpy},
        {"hdiss",        p.hdiss},
        {"hdiss_end",    p.hdiss_end},
        {"Bacc",         p.Bacc},
        {"rhoacc",       p.rhoacc},
        {"rhoload_end",  p.rhoload_end},
        {"mu_initial",   p.mu_initial},
        {"mu_end",       p.mu_end},
        {"rdiss",        p.rdiss}
    }};

    std::size_t width = 0;
    for (const auto& [name, value] : params) {
        width = std::max(width, name.size());
    }
    width += 2;

    os << "JetMassLoadingParameters:\n";
    for (const auto& [name, value] : params) {
        os << "  "
           << std::left << std::setw(static_cast<int>(width)) << name
           << " : " << value << '\n';
    }

    os << "\nMLJet kinematics ... " << std::endl;
    os << std::left << 
        std::setw(10) << "z_center" << 
        std::setw(10) << "radius" << 
        std::setw(10) << "z_height" << 
        std::setw(10) << "gamma" << 
        std::setw(15) << "elec_density" << 
        std::setw(15) << "prot_density" << 
        // std::setw(15) << "mass_density" << 
        std::setw(17) << "e_energy_density" << 
        std::setw(17) << "p_energy_density" << 
        std::setw(12) << "h" << 
        std::setw(12) << "sigma" << 
        std::setw(12) << "B-field" << 
        std::setw(15) << "mu" << 
        std::setw(15) << "<elec_gamma>" << 
        std::setw(15) << "<prot_gamma>" <<
    std::endl;

    for (size_t i = 0; i < n_zones; i++)
    {
        std::cout<<std::left << 
            std::setw(10) << z_center_grid[i] / r_g << 
            std::setw(10) << radius_grid[i] / r_g << 
            std::setw(10) << z_height_grid[i] / r_g <<
            std::setw(10) << gamma_grid[i] <<
            std::setw(15) << electron_density_grid[i] <<
            std::setw(15) << proton_density_grid[i] <<
            // std::setw(15) << mass_density_grid[i] <<
            std::setw(17) << electron_energy_density_grid[i] <<
            std::setw(17) << proton_energy_density_grid[i] <<
            std::setw(12) << specific_enthalpy_grid[i] <<
            std::setw(12) << sigma_grid[i] <<
            std::setw(12) << magnetic_field_grid[i] <<
            std::setw(15) << mu_grid[i] << 
            std::setw(15) << electron_average_gamma_grid[i] <<
            std::setw(15) << proton_average_gamma_grid[i] <<
        std::endl;
    }

}


} // namespace mljet