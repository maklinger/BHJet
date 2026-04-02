
#include "MLJet.hpp"
#include "iostream"
#include "filesystem"
#include "kariba/Thermal.hpp"

namespace karcst = kariba::constants;

namespace mljet
{


void MLJet::compute_jet_dynamics()
{
    if (verbosity_level > 1)
        std::cout << "Computing MLJet dynamics" << std::endl;

    // clean up existing gsl interpolation
    if (verbosity_level > 2)
        std::cout << "cleaning up spline memory" << std::endl;
    if (spline_speed)
        gsl_spline_free(spline_speed);


    // init basic internal params
    eddington_luminosity = 1.25e38 * mass_bh;
    r_g = karcst::gconst * mass_bh * karcst::msun / (karcst::cee * karcst::cee * karcst::erg);
    
    jet_dyn.min = z_jet_launching * r_g;
    jet_dyn.max = z_max_calculation * r_g;
    jet_dyn.h0 = 2. * r_initial * r_g + z_jet_launching * r_g;
    jet_dyn.r0 = r_initial * r_g;
    jet_dyn.acc = z_end_of_acceleration * r_g;
    jet_dyn.beta0 = sqrt(4. / 3. * (4. / 3. - 1.) /
                         (4. / 3. + 1.)); // set initial jet speed for relativistic fluid, g=4/3
    jet_dyn.gam0 = 1. / sqrt(1. - (std::pow(jet_dyn.beta0,
                                            2.))); // set corresponding lorentz factor
    jet_dyn.gamf = gamma_final;
    jet_dyn.Rg = r_g;


    kariba::Thermal dummy_elec(70); // use 70 bins for integral
    dummy_elec.set_temp_kev(electron_temperature_jet_base);
    dummy_elec.set_p();
    dummy_elec.set_norm(1.);
    dummy_elec.set_ndens();



    jet_mass_loading_parameters.pbeta = plasma_beta_jet_base;
    jet_mass_loading_parameters.Nj = jet_power_eddington * eddington_luminosity;
   	jet_mass_loading_parameters.sig_acc = sigma_final;
	jet_mass_loading_parameters.av_gamma = dummy_elec.av_gamma();
	// ge_av_0 = jet_mass_loading_parameters.av_gamma;

    /// @todo add this in the free parameter list
   	jet_mass_loading_parameters.eta=eta_e;



    try {
        calc_mass_loading_jet();
    }
    catch (const std::exception &e) {
        throw std::runtime_error(
            std::string("MLJet::compute_jet_dynamics: failed to compute mass loading jet: ")
            + e.what()
        );
    }


    /// @todo the velocity profile should be calculated somewhere
    // fill the spline for the speed profile (formerly "velprof_mag")
    // calc_velocity_profile_magnetized_jet();

    // grid building
    // double dlgz = 0.1;
    double zmin = jet_dyn.min;
    double zmax = jet_dyn.max;
    double zdiss = z_dissipation*r_g;

    if (zdiss < zmin*pow(10, dlgz*1.5)) {
        zdiss = zmin*pow(10, dlgz/2);
        if (verbosity_level > 1)
            std::cout << "z_dissipation<z_jet_launching, set z_dissipation=z_jet_launching*10^(dlgz/2)" << std::endl;
        
    } else if (zdiss > zmax*pow(10, -dlgz*1.5)) {
        zdiss = zmax*pow(10, -dlgz/2);
        if (verbosity_level > 1)
            std::cout << "z_dissipation>z_max_calculation, set z_dissipation=z_max_calculation*10^(-dlgz/2)" << std::endl;
    }

    n_zones= static_cast<size_t>(std::log10(zmax/zmin)/dlgz);
    size_t N1 = static_cast<size_t>(std::log10(zdiss*pow(10, -dlgz/2)/zmin)/dlgz);
    double dlgz1 = std::log10(zdiss*pow(10, -dlgz/2)/zmin)/static_cast<double>(N1);
    size_t N2 = n_zones - N1 - 1;
    double dlgz2 = std::log10(zmax/(zdiss*pow(10, dlgz/2)))/static_cast<double>(N2);
    

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
        // bljetpars(z, jet_dyn, nozzle_ener, tshift, zone, spline_speed, acc_speed);
    }


}


/// @brief calculate the jet parameters for the mass loading jet
/// @file MLJet.cpp
/// @exception std::runtime_error if the sigma magnetization at the dissipation region is not positive
/// @exception std::runtime_error if the approximation of the maximum particle energy fails due to 
/// average_gamma()
void MLJet::calc_mass_loading_jet()
{
	double dyn_fac = 0.0;
	double rho0 = 0.0;
	double theta_acc = 0.15/jet_dyn.gamf;
    
	//calculate the initial mass density
	dyn_fac = 2.*karcst::pi*jet_dyn.r0*jet_dyn.r0*jet_dyn.gam0*jet_dyn.beta0*karcst::cee;

    jet_mass_loading_parameters.sig0 = jet_mass_loading_parameters.pbeta;
    jet_mass_loading_parameters.lepdens = jet_mass_loading_parameters.Nj / 
        (dyn_fac*(karcst::pmgm/jet_mass_loading_parameters.eta+karcst::emgm) *
        karcst::cee*karcst::cee*(1.+jet_mass_loading_parameters.sig0));
    jet_mass_loading_parameters.protdens = jet_mass_loading_parameters.lepdens / jet_mass_loading_parameters.eta;
    
    /// @todo should I make them methods of the class ?
    rho0 = jet_mass_density(jet_mass_loading_parameters.protdens,jet_mass_loading_parameters.lepdens);

    jet_mass_loading_parameters.h0 = specific_enthalpy_simple(jet_mass_loading_parameters.av_gamma, 1.0);
    jet_mass_loading_parameters.mu_initial = jet_dyn.gam0 * 
        (jet_mass_loading_parameters.sig0+jet_mass_loading_parameters.h0+1.); //equal to μ_0
    jet_mass_loading_parameters.sig_acc = jet_mass_loading_parameters.mu_initial / jet_dyn.gamf - 
        (jet_mass_loading_parameters.h0+1.);


    if (jet_mass_loading_parameters.sig_acc < 0.0) {
        throw std::runtime_error(
            "MLJet::calc_mass_loading_jet: computed negative sig_acc = "
            + std::to_string(jet_mass_loading_parameters.sig_acc) + "\n" +
            "\tIncrease the initial magnetisation above: " + std::to_string(jet_mass_loading_parameters.sig0) + "\n" +
            "\tand/or decrease η below: " + std::to_string(jet_mass_loading_parameters.eta) + "\n" +
            "\tand/or decrease kBTe below: " + std::to_string(round((gamma_electron_min_dissipation-1.)*karcst::emerg*karcst::erg/1.0e3)) + "\n" +
            "\tand/or decrease gacc below: " + std::to_string(jet_dyn.gamf) + 
            "\n##########################################################\n" +
            "\tYou might want to keep:\n" + 
            "\t the initial magnetisation between 5 and 30\n" +
            "\t the initial temperature between 20 and 1000 keV\n" +
            "\t the maximum Lorentz factor gacc between 3 and 10\n"
        );
    }

    jet_mass_loading_parameters.hdiss = jet_mass_loading_parameters.h0;
    jet_mass_loading_parameters.rdiss = jet_dyn.r0 + (z_dissipation - jet_dyn.h0) * std::tan(theta_acc);
    jet_mass_loading_parameters.nacc = number_density_jet(jet_mass_loading_parameters.lepdens,jet_dyn.gamf,jet_mass_loading_parameters.rdiss);
    jet_mass_loading_parameters.rhoacc=jet_mass_density(
        jet_mass_loading_parameters.nacc/jet_mass_loading_parameters.eta,jet_mass_loading_parameters.nacc);//the mass density before mass loading
    jet_mass_loading_parameters.Bacc = 
        std::sqrt(4.*karcst::pi*jet_mass_loading_parameters.sig_acc*jet_mass_loading_parameters.rhoacc*karcst::cee*karcst::cee);

    //I keep the following to compare ge_av_diss with the one after loading and calculate tshift(need it?)

    /// @todo do I need actually the average gamma below? 
    try {
        approximate_max_particle_energy(gamma_electron_average_dissipation,
            gamma_proton_average_dissipation,z_dissipation);
    } catch (const std::exception &e) {
        throw std::runtime_error(
            std::string("MLJet::calc_mass_loading_jet: approximate_max_particle_energy failed: ")
            + e.what()
        );
    }

    double rload_end = jet_dyn.r0 + (100.0 * z_dissipation - jet_dyn.h0) * std::tan(theta_acc);
    double nload_end = number_density_jet(jet_mass_loading_parameters.lepdens,jet_dyn.gamf,rload_end);//before mass loading
    jet_mass_loading_parameters.rhoload_end = jet_mass_density(nload_end/jet_mass_loading_parameters.eta,nload_end);//before mass loading

	jet_mass_loading_parameters.bfield = 
            std::sqrt(4.*karcst::pi*jet_mass_loading_parameters.sig0*rho0*karcst::cee*karcst::cee);


}



double jet_mass_density(double nprotons, double nelectrons)
{
	return nprotons * karcst::pmgm + nelectrons * karcst::emgm;
}

double MLJet::specific_enthalpy_simple(double ge_av, double gp_av)
{
    return ( Ge * (ge_av - 1.) + Gp * (gp_av - 1.) * karcst::pmgm / (karcst::emgm * eta_e) ) / 
            (1. + karcst::pmgm / karcst::emgm / eta_e);
}

double MLJet::number_density_jet(double n0, double g, double r)
{
	double beta_j;
	beta_j = (g<=5.) ? sqrt(g*g-1.)/g : 1.-0.5*pow(g,-2) ;
	return n0*jet_dyn.gam0*jet_dyn.beta0/(g*beta_j)*pow(r/jet_dyn.r0,-2);
}


/// @brief calculate the max particle energy at the dissipation region
/// @param[inout] ge_av average gamma of the electrons
/// @param[inout] gp_av average gamma of the protons
/// @param z the jet segment distance along the jet axis
/// @file MLJet.cpp
/// @exception std::runtime_error if the approximation fails @todo check if needs to be changed
void MLJet::approximate_max_particle_energy(double &ge_av, double &gp_av, double z)
{
	/*I approximately find the max energy of accelerated electrons (and protons).
	For the electrons I neglect any photon field and for protons, I only account for escape*/
	double theta = 0.15 / jet_dyn.gamf;
	double Rapprox = jet_dyn.r0 + std::max(z-jet_dyn.h0,0.)*std::tan(theta);  //for max particle energy
	double Bapprox = jet_mass_loading_parameters.Bacc * std::pow(z/jet_dyn.acc,-1);
	//approximation for the max energy of electrons; I neglect any photon field so TODO think how to
	//improve it in case the max energy is constrained by ICS
	// double Tacc0,Tesc,Tsyn0,D,gmax;//,escom,accon,syncon,b,c, gmax;

    /// synchrotron losses timescale constant
	double Tsyn0 = 6.*karcst::pi*karcst::emerg*karcst::emerg/(karcst::sigtom*Bapprox*Bapprox*karcst::cee);
    /// accelerated electrons timescale constant
    double Tesc = factor_break_electrons*Rapprox/karcst::cee;
	/// acceleration electrons timescale constant
    double Tacc0 = 4./(3.*factor_max_energy_electrons*karcst::charg*karcst::cee*Bapprox);
	/// constant for the calculation of the maximum gamma
    double D = std::pow(Tsyn0/Tesc,2)+4.*Tsyn0/Tacc0;
	/// maximum gamma of electrons/protons
    double gmax = (-Tsyn0/Tesc+std::sqrt(D))/(2.*karcst::emerg);
	
    ge_av = average_gamma(index_injected_electrons,gamma_electron_min_dissipation,gmax); 
	//TODO I assume here that the electrons are not cooled at all but might have to compare with the final Eemax
		
    gp_av=1.;
	if (gamma_proton_min_dissipation>=1.){
		gmax = Tesc/Tacc0/(karcst::pmgm*karcst::cee*karcst::cee);
		gp_av = average_gamma(index_injected_protons,gamma_proton_min_dissipation,gmax);
	}
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






void MLJet::reinit_grid_arrays_for_mass_loading()
{
    if (verbosity_level > 1)
        std::cout << "MLJet::reinit_grid_arrays_for_mass_loading" << std::endl;
    
    mu_grid = std::vector<double>(n_zones, 0.0);
    reinit_grid_arrays();
}








/// @todo the following function is the one that sets the jet parameters, must be in a function and fill in a vectro (?)

/*We parametrise Chatterjee et al 2019 mass loading paper; heat=Up/Ue */
/// It was the loadingjetpars() function
void MLJet::calc_zone_properties(size_t i){
// (double z,jet_dynpars &dyn,jet_enpars &en,particle_pars &prtcls,zone_pars &zone,double &t,
// 	int infosw,double load_frac,double heat,std::string outputConfiguration){
	
    double theta;
	double ge_av,gp_av;

	if (z_center_grid[i] < z_dissipation)
    { /// @note the pre-loading region

		gamma_grid[i] = std::pow(10.,log10(jet_dyn.gam0) + 
				std::log10(jet_dyn.gamf/jet_dyn.gam0)/std::log10(z_dissipation/jet_dyn.min) * std::log10(z_center_grid[i]/jet_dyn.min));
		mu_grid[i] = jet_mass_loading_parameters.mu_initial;
		
        /// @todo I need to define h that at jet base is h0
        /// @todo should I define the zone.helec as a grid vector? 
        specific_enthalpy_grid[i] = zone.helec;
		sigma_grid[i] = mu_grid[i] / gamma_grid[i] - (specific_enthalpy_grid[i] + 1.);

		theta = 0.15 / gamma_grid[i];
		radius_grid[i] = jet_dyn.r0 + std::max(z_center_grid[i] - jet_dyn.h0, 0.) * std::tan(theta);
		/// @todo should I define the zone.lepdens as a grid vector?
		/// @todo should I define the zone.protdens as a grid vector?
        zone.lepdens = number_density_jet(jet_mass_loading_parameters.lepdens, gamma_grid[i], radius_grid[i]);
		zone.protdens = zone.lepdens/jet_mass_loading_parameters.eta;
		mass_density_grid[i] = jet_mass_density(zone.protdens, zone.lepdens);
		
        zone.Up = internal_energy_density(1,zone.protdens,karcst::pmgm); //always zero assuming cold protons
		zone.Ue = specific_enthalpy_grid[i]*mass_density_grid[i]*karcst::cee*karcst::cee/Ge;
		
		ge_av = jet_mass_loading_parameters.av_gamma;
		gp_av = 1.;
		
	}else if (z_center_grid[i]>=z_dissipation && z_center_grid[i]<=100.*z_dissipation){		
		
        /// @todo figure out that eta is changing after that point properly
        if (en.eta>10.) en.eta=10.;
//		en.eta=1.;//after mass loading I switch to 10
		double tanh2=pow(tanh(log10(z_center_grid[i]/z_dissipation)),2);
		sigma_grid[i] = sigma_loading(z_center_grid[i],z_dissipation)*en.sig_acc/sigma_loading(z_dissipation,z_dissipation);
		
        /// @todo the heat parameter is not used so far, I need to pass it in the input
        if (heat>0 && heat<1.) sigma_grid[i] *=(1.-tanh2*heat);//a fraction "heat" goes to h
		gamma_grid[i] = gamma_loading(z_center_grid[i],z_dissipation)*jet_dyn.gamf/gamma_loading(z_dissipation,z_dissipation);
		specific_enthalpy_grid[i] = std::max({h_loading(z_center_grid[i],z_dissipation)*\
				en.hdiss/h_loading(z_dissipation,z_dissipation),zone.helec,en.h0});
		if (heat>0 && heat<1.) specific_enthalpy_grid[i]+=tanh2*heat*sigma_grid[i]; //a fraction "heat" comes from σ 
		//If h leads to an increas of mu:
		
		//I'm trying to fix some extremes
		if (z_center_grid[i]<9.*z_dissipation){
			mu_grid[i]=en.muload;//force the first few segements to be constant
			if (gamma_grid[i]*(sigma_grid[i]+specific_enthalpy_grid[i]+1.)>en.muload){
				//when h from the profile forces mu to increase in the first segments:
				specific_enthalpy_grid[i] = std::max({mu_grid[i]/gamma_grid[i] - (sigma_grid[i]+1.),zone.helec,en.h0});
			}
		}else{
			if (gamma_grid[i]*(sigma_grid[i]+specific_enthalpy_grid[i]+1.)>en.muend){
				//when h from the profile forces mu to increase beyound the 8*zdiss segments:
				mu_grid[i] = mu_loading(z_center_grid[i],z_dissipation)*en.muload/mu_loading(z_dissipation,z_dissipation);
                specific_enthalpy_grid[i] = std::max({mu_grid[i]/gamma_grid[i] - (sigma_grid[i]+1.),zone.helec});
				//in case the h is not from the mu equation, I have to fix gamma:
				if (gamma_grid[i]*(sigma_grid[i]+specific_enthalpy_grid[i]+1.)>mu_grid[i]) gamma_grid[i] = mu_grid[i]/(sigma_grid[i]+specific_enthalpy_grid[i]+1.);
			}else{
				mu_grid[i] = gamma_grid[i]*(sigma_grid[i]+specific_enthalpy_grid[i]+1.);
			}
			
		}
		theta = 0.15/gamma_grid[i];
		radius_grid[i] = jet_dyn.r0+std::max(z_center_grid[i]-jet_dyn.h0,0.)*tan(theta);

		//the factor 2 comes from log10(100*zdiss/zdiss)
		mass_density_grid[i] = pow(10.,log10(en.rhoacc)+log10(load_frac*en.rhoload_end/en.rhoacc)/2.*log10(z_center_grid[i]/z_dissipation));
		zone.lepdens = mass_density_grid[i]/(karcst::pmgm/en.eta+karcst::emgm);
		zone.protdens = zone.lepdens/en.eta; //I assume η=10 for loaded matter

		//the very first time, helec is from pre-shock region; beyond zdiss it's from the previous segment
		zone.Ue = zone.helec*mass_density_grid[i]*karcst::cee*karcst::cee/Ge;
		zone.Up = std::max(1.e-50,(specific_enthalpy_grid[i]*mass_density_grid[i]*karcst::cee*karcst::cee-Ge*zone.Ue)/Gp);
		ge_av = std::max(zone.Ue/(zone.lepdens*karcst::emerg)+1.,1.);
	    gp_av = std::max(zone.Up/(zone.protdens*karcst::pmgm*karcst::cee*karcst::cee)+1.,1.);
	    en.muend = mu_grid[i];
	    en.hdiss_end = specific_enthalpy_grid[i];
	}else{
		double tanh2=pow(tanh(log10(z_center_grid[i]/z_dissipation)),2);
        //The one degree polynomial fit of CLTM19 between 1e4-1e5rg
        sigma_grid[i] = pow(10.,-0.09708*log10(z_center_grid[i]/z_dissipation)-0.178)*en.sig_acc/sigma_loading(z_dissipation,z_dissipation);
		if (heat>0 && heat<1.) sigma_grid[i] *=(1.-tanh2*heat);//a fraction "heat" goes to h
    	 specific_enthalpy_grid[i] =std::max({zone.helec,\
	    	pow(10.,-0.2451*log10(z_center_grid[i]/z_dissipation)-0.5763)*en.hdiss_end/h_loading(100.*z_dissipation,z_dissipation)});
		if (heat>0 && heat<1.) specific_enthalpy_grid[i]+=tanh2*heat*sigma_grid[i];
		mu_grid[i] = en.muend;
		gamma_grid[i] = std::max(1.,mu_grid[i]/(sigma_grid[i]+specific_enthalpy_grid[i]+1.));

		theta = 0.15/gamma_grid[i];
		radius_grid[i] = jet_dyn.r0+std::max(z_center_grid[i]-jet_dyn.h0,0.)*tan(theta);
		//based on CLTM19 the mass density has increased by a factor of load_frac
		mass_density_grid[i] = load_frac*en.rhoload_end*pow(z_center_grid[i]/(100.*z_dissipation),-2);
		zone.lepdens = mass_density_grid[i]/(karcst::pmgm/en.eta+karcst::emgm);
		zone.protdens = zone.lepdens/en.eta; //I assume η=1 for loaded matter

		zone.Ue = zone.helec*mass_density_grid[i]*karcst::cee*karcst::cee/Ge;
		zone.Up = std::max(1.e-50,(specific_enthalpy_grid[i]*mass_density_grid[i]*karcst::cee*karcst::cee-Ge*zone.Ue)/Gp);
		ge_av = std::max(zone.Ue/(zone.lepdens*karcst::emerg)+1.,1.);
        gp_av = std::max(zone.Up/(zone.protdens*karcst::pmgm*karcst::cee*karcst::cee)+1.,1.);
	}
	zone.bfield = sqrt(4.*karcst::pi*mass_density_grid[i]*karcst::cee*karcst::cee*sigma_grid[i]*(1.+specific_enthalpy_grid[i]));
	zone.beta = pow(1.-1./gamma_grid[i]/gamma_grid[i],0.5);
	t=1.;

	if (infosw>=2) print_zone_kinematics(outputConfiguration,z_center_grid[i],jet_dyn.Rg,zone,mass_density_grid[i],specific_enthalpy_grid[i],sigma_grid[i],gp_av,ge_av);	
}



double internal_energy_density(double g_av,double n,double mass){
	//the mass in grams
	return (g_av - 1.) * n * mass * karcst::cee * karcst::cee;
}



} // namespace mljet