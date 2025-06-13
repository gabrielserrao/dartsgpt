#include "dartsflash/flash/flash.hpp"
#include "dartsflash/flash/px_flash.hpp"
#include "dartsflash/stability/stability.hpp"
#include "dartsflash/eos/helmholtz/cubic.hpp"
#include "dartsflash/eos/aq/jager.hpp"
#include "dartsflash/eos/aq/ziabakhsh.hpp"
#include "dartsflash/eos/vdwp/ballard.hpp"
#include "dartsflash/eos/solid/solid.hpp"
#include "dartsflash/global/global.hpp"

int test_purecomponent_ph();
int test_phflash_vapour_liquid();
int test_phflash_vapour_liquid_water();

struct Reference
{
	double pressure, h_s, T_tol{ 1e-1 }, T_ref;
	std::vector<double> composition;

	Reference(const double p, const double h_spec, const std::vector<double>& z, const double T_)
	: pressure(p), h_s(h_spec), T_ref(T_), composition(z) {}

	int test(PXFlash *flash, bool verbose, bool test_derivs=false)
	{
		if (verbose)
		{
			std::cout << "==================================\n";
			print("p, h_spec", std::vector<double>{pressure, h_s});
			print("z", composition);
		}
		int error = flash->evaluate(pressure, h_s, composition);
		if (error > 0)
		{
			print("Error in Flash", error);
			return error;
		}

		// Output and compare results
		double T_root = flash->getT();

		if (verbose)
		{
			std::cout << "\nResults:\n";
			print("nu", flash->getnu());
			print("X", flash->getx());
			print("T", T_root);
		}

		if (std::fabs(T_root-T_ref) > T_tol)
		{
			std::cout << "T and T_ref are not the same \n";
			print("T", T_root);
			print("T_ref", T_ref);
			error++;
		}

		if (test_derivs && this->test_derivatives(flash, verbose))
		{
			std::cout << "Error in PXFlash::test_derivatives()\n";
			error++;
		}
		return error;
	}

	int test_derivatives(Flash *flash, bool verbose)
	{
		int nc = static_cast<int>(composition.size());
		int error = 0;
		
		// Output results and derivatives
		FlashResults flash_results0 = flash->get_flash_results(true);
		int np = flash_results0.np;
		// int np_tot = flash_results0.np_tot;

		std::vector<double> dnudP, dnudT, dnudzk, dnudH, dxdP, dxdT, dxdzk, dxdH;
		flash_results0.dP_derivs(dnudP, dxdP);
		flash_results0.dT_derivs(dnudT, dxdT);
		flash_results0.dz_derivs(dnudzk, dxdzk);
		flash_results0.dX_derivs(EoS::Property::ENTHALPY, dnudT, dxdT, dnudH, dxdH);
		double d, dX_tol{ 1e-4 };

		// // Calculate numerical derivatives w.r.t. pressure
		// double dp = 1e-6 * pressure;
		// error += flash->evaluate(pressure+dp, h_s, composition);
		// FlashResults flash_results1 = flash->get_flash_results(false);

		// // Compare analytical and numerical
		// for (int jj = 0; jj < np; jj++)
		// {
		// 	int j = flash_results0.phase_idxs[jj];

		// 	double dnujdP_num = (flash_results1.nuj[j] - flash_results0.nuj[j])/dp;
		// 	// Use logarithmic scale to compare
		// 	d = std::log(std::fabs(dnujdP_num + 1e-15)) - std::log(std::fabs(dnudP[j] + 1e-15));
		// 	if (verbose || std::isnan(dnudP[j]) || (std::fabs(d) > dX_tol && (std::fabs(dnudP[j]) > 1e-8 && std::fabs(dnujdP_num) > 1e-8)))
		// 	{
		// 		print("phase", j);
		// 		print("dnuj/dP", std::vector<double>{dnudP[j], dnujdP_num, d});
		// 		error++;
		// 	}

		// 	for (int i = 0; i < nc; i++)
		// 	{
		// 		double dxdP_num = (flash_results1.Xij[j*nc+i] - flash_results0.Xij[j*nc+i])/dp;
		// 		// Use logarithmic scale to compare
		// 		d = std::log(std::fabs(dxdP_num + 1e-15)) - std::log(std::fabs(dxdP[j*nc+i] + 1e-15));
		// 		if (verbose || std::isnan(dxdP[j*nc+i]) || (std::fabs(d) > dX_tol && (std::fabs(dxdP[j*nc+i]) > 1e-8 && std::fabs(dxdP_num) > 1e-8)))
		// 		{
		// 			print("phase, comp", {j, i});
		// 			print("dXij/dP", std::vector<double>{dxdP[j*nc+i], dxdP_num, d});
		// 			error++;
		// 		}
		// 	}
		// }

		// Calculate numerical derivatives w.r.t. enthalpy
		double dH = 1e-6 * h_s;
		error += flash->evaluate(pressure, h_s + dH, composition);
		FlashResults flash_results1 = flash->get_flash_results(false);

		// Compare analytical and numerical
		for (int jj = 0; jj < np; jj++)
		{
			int j = flash_results0.phase_idxs[jj];

			double dnujdH_num = (flash_results1.nuj[j] - flash_results0.nuj[j])/dH;
			// Use logarithmic scale to compare
			d = std::log(std::fabs(dnujdH_num + 1e-15)) - std::log(std::fabs(dnudH[j] + 1e-15));
			if (verbose || std::isnan(dnudH[j]) || (std::fabs(d) > dX_tol && (std::fabs(dnudH[j]) > 1e-8 && std::fabs(dnujdH_num) > 1e-8)))
			{
				print("phase", j);
				print("dnuj/dH", std::vector<double>{dnudH[j], dnujdH_num, d});
				error++;
			}

			for (int i = 0; i < nc; i++)
			{
				double dxdH_num = (flash_results1.Xij[j*nc+i] - flash_results0.Xij[j*nc+i])/dH;
				// Use logarithmic scale to compare
				d = std::log(std::fabs(dxdH_num + 1e-15)) - std::log(std::fabs(dxdH[j*nc+i] + 1e-15));
				if (verbose || std::isnan(dxdH[j*nc+i]) || (std::fabs(d) > dX_tol && (std::fabs(dxdH[j*nc+i]) > 1e-8 && std::fabs(dxdH_num) > 1e-8)))
				{
					print("phase, comp", {j, i});
					print("dXij/dH", std::vector<double>{dxdH[j*nc+i], dxdH_num, d});
					error++;
				}
			}
		}

		// // Calculate numerical derivatives w.r.t. composition
		// std::vector<double> z(nc);
		// double nT_inv = 1./std::accumulate(composition.begin(), composition.end(), 0.);
		// std::transform(composition.begin(), composition.end(), z.begin(), [&nT_inv](double element) { return element *= nT_inv; });

		// double dz = 1e-6;
		// for (int k = 0; k < nc - 1; k++)
		// {
		// 	// Transform to +dz
		// 	double dzk = dz * z[k];
		// 	z[k] += dzk;
		// 	z[nc-1] -= dzk;
			
		// 	// Numerical derivative of lnphi w.r.t. zk
		// 	error += flash->evaluate(pressure, h_s, z);
		// 	flash_results1 = flash->get_flash_results(false);

		// 	// Compare analytical and numerical
		// 	for (int jj = 0; jj < np; jj++)
		// 	{
		// 		int j = flash_results0.phase_idxs[jj];

		// 		double dnujdzk_num = (flash_results1.nuj[j] - flash_results0.nuj[j])/dzk;
		// 		// Use logarithmic scale to compare
		// 		d = std::log(std::fabs(dnujdzk_num + 1e-15)) - std::log(std::fabs(dnudzk[k * np_tot + j] + 1e-15));
		// 		if (verbose || std::isnan(dnudzk[k * np_tot + j]) || (std::fabs(d) > dX_tol && (std::fabs(dnudzk[k * np_tot + j]) > 1e-8 && std::fabs(dnudzk[k * np_tot + j]) > 1e-8)))
		// 		{
		// 			print("phase, zk", {j, k});
		// 			print("dnuj/dzk", std::vector<double>{dnudzk[k * np_tot + j], dnujdzk_num, d});
		// 			error++;
		// 		}

		// 		for (int i = 0; i < nc-1; i++)
		// 		{
		// 			double dxdzk_num = (flash_results1.Xij[j*nc+i] - flash_results0.Xij[j*nc+i])/dzk;
		// 			// Use logarithmic scale to compare
		// 			d = std::log(std::fabs(dxdzk_num + 1e-15)) - std::log(std::fabs(dxdzk[k * np_tot * nc + j*nc+i] + 1e-15));
		// 			if (verbose || std::isnan(dxdzk[k * np_tot * nc + j*nc+i]) || (std::fabs(d) > dX_tol && (std::fabs(dxdzk[k * np_tot * nc + j*nc+i]) > 1e-8 && std::fabs(dxdzk[k * np_tot * nc + j*nc+i]) > 1e-8)))
		// 			{
		// 				print("phase, comp, zk", {j, i, k});
		// 				print("dXij/dzk", std::vector<double>{dxdzk[k * np_tot * nc + j*nc+i], dxdzk_num, d});
		// 				error++;
		// 			}
		// 		}	
		// 	}

		// 	// Return to original z
		// 	z[k] -= dzk;
		// 	z[nc-1] += dzk;
		// }
		return error;
	}
};

int main()
{
	int error_output = 0;
	
	error_output += test_purecomponent_ph();
	error_output += test_phflash_vapour_liquid();
	error_output += test_phflash_vapour_liquid_water();

    return error_output;
}

int test_purecomponent_ph()
{
	// Test pure component PH-flash
	// Test pure CO2 going from liquid(-like) to vapour(-like)
	// Test pure H2O going from ice to liquid to vapour
	bool verbose = false;
	std::cout << (verbose ? "TESTING 1C FLASH VAPOUR-LIQUID AT P-H\n" : "");
    int error_output = 0;

	std::vector<std::string> comp = {"CO2"};
	std::vector<double> z = {1.};
	CompData comp_data(comp);
	comp_data.Pc = {comp_data::Pc["CO2"]};
	comp_data.Tc = {comp_data::Tc["CO2"]};
	comp_data.ac = {comp_data::ac["CO2"]};
	comp_data.Mw = {comp_data::Mw["CO2"]};
	comp_data.kij = std::vector<double>(1, 0.);

	CubicEoS ceos(comp_data, CubicEoS::PR);

	FlashParams flash_params(comp_data);
	flash_params.verbose = verbose;

	flash_params.add_eos("CEOS", &ceos);
	flash_params.eos_params["CEOS"].root_order = {EoS::RootFlag::MAX, EoS::RootFlag::MIN};

	PXFlash flash(flash_params, PXFlash::StateSpecification::ENTHALPY);

	std::vector<double> pressure = linspace(10., 100., 10);
	double hmin = ceos.H(pressure.back(), 260., z, 0, true) * M_R;
	double hmax = ceos.H(pressure.back(), 350., z, 0, true) * M_R;
	std::vector<double> enthalpy = linspace(hmin, hmax, 20);

	for (double h: enthalpy)
	{
		for (double p: pressure)
		{
			int error = flash.evaluate(p, h);
			if (verbose || error) { print("Error in CO2 PHFlash", {p, h}); error_output += error; }
			// error_output += (std::fabs(flash.getH(flash.getT()) - T) < 1e-3) ? 0 : 1;
		}
	}

	// Define H2O
	comp = {"H2O"};
	comp_data = CompData(comp);
	comp_data.Pc = {comp_data::Pc["H2O"]};
	comp_data.Tc = {comp_data::Tc["H2O"]};
	comp_data.ac = {comp_data::ac["H2O"]};
	comp_data.Mw = {comp_data::Mw["H2O"]};
	comp_data.kij = std::vector<double>(1., 0.);

	flash_params = FlashParams(comp_data);
	flash_params.verbose = verbose;

	ceos = CubicEoS(comp_data, CubicEoS::PR);
	ceos.set_preferred_roots(0, 0.75, EoS::RootFlag::MAX);
	flash_params.add_eos("CEOS", &ceos);

	std::map<AQEoS::CompType, AQEoS::Model> evaluator_map = {
		{AQEoS::CompType::water, AQEoS::Model::Jager2003},
		{AQEoS::CompType::solute, AQEoS::Model::Ziabakhsh2012},
		{AQEoS::CompType::ion, AQEoS::Model::Jager2003}
	};
	AQEoS aq(comp_data, evaluator_map);
	aq.set_eos_range(0, std::vector<double>{0.6, 1.});
	flash_params.add_eos("AQ", &aq);

	flash_params.eos_order = {"AQ", "CEOS"};
	flash_params.T_min = 250.;
	flash_params.T_max = 500.;

	flash = PXFlash(flash_params, PXFlash::StateSpecification::ENTHALPY);

	pressure = linspace(10., 100., 10);
	hmin = aq.H(pressure.back(), 270., z, 0, true) * M_R;
	hmax = ceos.H(pressure.back(), 470., z, 0, true) * M_R;
	enthalpy = linspace(hmin, hmax, 20);

	for (double h: enthalpy)
	{
		for (double p: pressure)
		{
			int error = flash.evaluate(p, h);
			if (verbose || error) { print("Error in H2O PHFlash", {p, h}); error_output += error; }
			// error_output += (std::fabs(flash.getH(flash.getT()) - T) < 1e-3) ? 0 : 1;
		}
	}

	if (error_output > 0)
	{
		// std::cout << ref_string;
		print("Errors occurred in test_purecomponent_ph()", error_output);
	}
	else
	{
		print("No errors occurred in test_purecomponent_ph()", error_output);
	}
    return error_output;
}

int test_phflash_vapour_liquid()
{
	// Test C1/C4 mixture (Zhu, 2014), data from (Zhu, 2014)
	bool verbose = false;
	bool test_derivs = true;
	std::cout << (verbose ? "TESTING PH-FLASH WITH NP STABILITYFLASH FOR BINARY MIXTURE\n" : "");
	int error_output = 0;

	std::vector<std::string> comp = {"C1", "C4"};
	CompData comp_data = CompData(comp);
	comp_data.Pc = {46.0, 38.0};
	comp_data.Tc = {190.60, 425.20};
	comp_data.ac = {0.008, 0.193};
	comp_data.kij = std::vector<double>(2*2, 0.);
	comp_data.T_0 = 273.15;

	std::vector<double> z = {0.99, 0.01};

	FlashParams flash_params(comp_data);
	flash_params.split_variables = FlashParams::lnK;
    flash_params.split_switch_tol = 1e-4;
	flash_params.stability_variables = FlashParams::alpha;
	flash_params.verbose = verbose;

	CubicEoS pr(comp_data, CubicEoS::PR);
	flash_params.add_eos("PR", &pr);
    flash_params.eos_params["PR"].stability_switch_tol = 1e-3;
	flash_params.eos_params["PR"].initial_guesses = {InitialGuess::Yi::Wilson,
													 InitialGuess::Yi::Wilson13};  // pure H2O initial guess
	flash_params.eos_params["PR"].root_order = {EoS::RootFlag::MAX, EoS::RootFlag::MIN};
    flash_params.T_min = {100};
	flash_params.T_max = {900};
	flash_params.T_init = {400};

	PXFlash flash(flash_params, PXFlash::StateSpecification::ENTHALPY);

	std::vector<Reference> references = {
		Reference(50., -6500, z, 195.6676379),
	};

	for (Reference condition: references)
	{
		error_output += condition.test(&flash, verbose, test_derivs);
	}

	if (error_output > 0)
	{
		print("Errors occurred in test_phflash_vapour_liquid()", error_output);
	}
	else
	{
		print("No errors occurred in test_phflash_vapour_liquid", error_output);
	}
    return error_output;
}

int test_phflash_vapour_liquid_water()
{
	// Test Water/NWE mixture (Khan et al, 1992), data from (Li, 2018)
	bool verbose = false;
	bool test_derivs = true;
	std::cout << (verbose ? "TESTING PH-FLASH WITH NP STABILITYFLASH FOR WATER MIXTURE\n" : "");
	int error_output = 0;

	std::vector<std::string> comp = {"H2O", "PC1", "PC2", "PC3", "PC4"};
	CompData comp_data = CompData(comp);
	comp_data.Pc = {220.89, 48.82, 19.65, 10.20, 7.72};
	comp_data.Tc = {647.3, 305.556, 638.889, 788.889, 838.889};
	comp_data.ac = {0.344, 0.098, 0.535, 0.891, 1.085};
	comp_data.kij = std::vector<double>(5*5, 0.);
	comp_data.set_binary_coefficients(0, {0., 0.71918, 0.45996, 0.26773, 0.24166});
	
	comp_data.cpi = {comp_data::cpi["H2O"],
					{-3.5 / M_R, 0.005764 / M_R, 5.09E-7 / M_R, 0. },
					{-0.404 / M_R, 0.0006572 / M_R, 5.41E-8 / M_R, 0.},
					{-6.1 / M_R, 0.01093 / M_R, 1.41E-6 / M_R, 0.},
					{-4.5 / M_R, 0.008049 / M_R, 1.04E-6 / M_R, 0.}};
	comp_data.T_0 = 273.15;

	std::vector<double> z = {0.5, 0.15, 0.10, 0.10, 0.15};

	FlashParams flash_params(comp_data);
	flash_params.split_variables = FlashParams::nik;
    flash_params.split_switch_tol = 1e-3;
	flash_params.stability_variables = FlashParams::Y;
	flash_params.verbose = verbose;

	CubicEoS pr(comp_data, CubicEoS::PR);
	flash_params.add_eos("PR", &pr);
    flash_params.eos_params["PR"].stability_switch_tol = 1e-1;
	flash_params.eos_params["PR"].initial_guesses = {InitialGuess::Yi::Wilson,
													 InitialGuess::Yi::Wilson13};  // pure H2O initial guess
	flash_params.eos_params["PR"].root_order = {EoS::RootFlag::MAX, EoS::RootFlag::MIN};
    flash_params.T_min = {200};
	flash_params.T_max = {900};
	flash_params.T_init = {400};

	PXFlash flash(flash_params, PXFlash::StateSpecification::ENTHALPY);

	std::vector<Reference> references = {
		Reference(30., 0, z, 742.7160),
		Reference(60., 0, z, 782.646),
		Reference(90., 0, z, 828.752),
	};

	for (Reference condition: references)
	{
		error_output += condition.test(&flash, verbose, test_derivs);
	}

	if (error_output > 0)
	{
		print("Errors occurred in test_phflash_vapour_liquid_water()", error_output);
	}
	else
	{
		print("No errors occurred in test_phflash_vapour_liquid_water()", error_output);
	}
    return error_output;
}