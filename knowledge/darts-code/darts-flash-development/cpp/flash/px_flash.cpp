#include <iostream>
#include <vector>
#include <cassert>
#include <numeric>
#include <functional>

#include "dartsflash/maths/root_finding.hpp"
#include "dartsflash/flash/flash_results.hpp"
#include "dartsflash/flash/px_flash.hpp"

#include <Eigen/Dense>

PXFlash::PXFlash(FlashParams& flashparams, PXFlash::StateSpecification state_spec_) : Flash(flashparams)
{
	this->state_spec = state_spec_;
}

void PXFlash::init(double p_, double X_spec_)
{
	// Initialize PH-flash algorithm at P and Hspec
	p = p_;
	X_spec = X_spec_;

	this->T_a = this->flash_params.T_min;
	this->T_b = this->flash_params.T_max;
	this->T = this->flash_params.T_init;
	return;
}

void PXFlash::init(double p_, double X_spec_, std::vector<double>& z_)
{
	// Initialize PH-flash algorithm at P, Hspec and z
	this->init(p_, X_spec_);

    // Check if feed composition needs to be corrected for 0 values
    z.resize(flash_params.ns);
    for (int i = 0; i < flash_params.ns; i++)
    {
        z[i] = (z_[i] > flash_params.min_z) ? z_[i] : flash_params.min_z;
    }
    return;
}

int PXFlash::evaluate(double p_, double X_spec_)
{
	PXFlash::init(p_, X_spec_);
	this->z = {1.};

	RootFinding root;
   	auto f = std::bind(&PXFlash::obj_fun, this, std::placeholders::_1);
	// auto df = std::bind(&PXFlash::gradient, this, std::placeholders::_1);
	this->T = NAN;
   	// int output = root.brent_newton(f, df, this->T, T_a, T_b, flash_params.phflash_Htol, flash_params.phflash_Ttol);
	int output = root.brent(f, this->T, T_a, T_b, flash_params.phflash_Htol, flash_params.phflash_Ttol);
	this->T = root.getx();

   	if (output == 0)
	{
		return 0;
	}
	else if (output == -1)
	{
		// T_max - T_min < eps, H-Hspec > eps
		// Indicates that solution is oscillating between two phases and solution is 2P
		std::vector<TrialPhase> comps_a, comps_b;

		// Stable composition at point b equals latest Flash results
		FlashResults results = this->get_flash_results();
		double X_b = XX;
		int phase_b = std::distance(results.nuj.begin(), std::find_if(results.nuj.begin(), results.nuj.end(), [](double nu_j) { return nu_j > 0.; }));
		int eos_idx_b = results.eos_idx[phase_b];
		EoS::RootFlag root_b = results.root_type[phase_b];

		// Stable composition at point a
		int iter = 0;
		double X_a;
		int phase_a, eos_idx_a;
		EoS::RootFlag root_a;
		double t_incr = (X_b - X_spec > 0.) ? -flash_params.phflash_Ttol : flash_params.phflash_Ttol;
		while (true)
		{
			iter++;
			T_a += t_incr;
			this->obj_fun(T_a);
			X_a = XX;

			results = this->get_flash_results();
			phase_a = std::distance(results.nuj.begin(), std::find_if(results.nuj.begin(), results.nuj.end(), [](double nu_j) { return nu_j > 0.; }));
			eos_idx_a = results.eos_idx[phase_a];
			root_a = results.root_type[phase_a];
			if ((eos_idx_a != eos_idx_b) || // If PT flash gives same EoS
				(root_a != root_b))  // and same root, evaluate the other root to combine the two ref compositions
			{
				break;
			}
			if (iter > 100)
			{
				return 1;
			}
		}

		// Find phase fraction nu_a at solution of X = Xspec
		// X = nu_a * X_a + (1-nu_a) * X_b;
		double nu_a = (X_spec-X_b)/(X_a-X_b);
		this->ref_compositions = {TrialPhase{eos_idx_a, flash_params.eos_order[eos_idx_a], z}, TrialPhase{eos_idx_b, flash_params.eos_order[eos_idx_b], z}};
		ref_compositions[0].nu = nu_a;
		ref_compositions[0].root = root_a;
		ref_compositions[1].nu = 1.-nu_a;
		ref_compositions[1].root = root_b;

		return 0;
	}
	else if (output == 1)
   	{
      	if (flash_params.verbose)
		{
			std::cout << "ERROR in PXFlash\n";
        	print("p, X", std::vector<double>{p, X_spec}, 10);
		}
      	return 1;
   	}
   	else
   	{
      	return 0;
   	}
}

int PXFlash::evaluate(double p_, double X_spec_, std::vector<double>& z_)
{
	if (z_.size() == 1)
    {
        return PXFlash::evaluate(p_, X_spec_);
    }
	PXFlash::init(p_, X_spec_, z_);

	RootFinding root;
   	auto f = std::bind(&PXFlash::obj_fun, this, std::placeholders::_1);
	// auto df = std::bind(&PXFlash::gradient, this, std::placeholders::_1);
	this->T = NAN;
   	// int output = root.brent_newton(f, df, this->T, T_a, T_b, flash_params.phflash_Htol, flash_params.phflash_Ttol);
	int output = root.brent(f, this->T, T_a, T_b, flash_params.phflash_Htol, flash_params.phflash_Ttol);
	this->T = root.getx();

	if (output == 0)
	{
		if (flash_params.verbose)
        {
            print("PXFlash", "===============");
            print("p, X", std::vector<double>{p, X_spec}, 1, 10);
			print("z", z_, 1, 10);
            this->get_flash_results().print_results();
        }
		return 0;
	}
   	else if (output == 1)
   	{
      	if (flash_params.verbose)
		{
			std::cout << "ERROR in PXFlash\n";
        	print("p, X", std::vector<double>{p, X_spec}, 10);
			print("z", z_, 10);
		}
      	return 1;
   	}
   	else
   	{
      	return 0;
   	}
}

double PXFlash::obj_fun(double T_)
{
	// Perform PT-flash at P,T,z
	if (flash_params.verbose)
	{
		print("PT evaluation at P, T", std::vector<double>{p, T_});
	}
	error = (flash_params.nc == 1) ? Flash::evaluate(p, T_) : Flash::evaluate(p, T_, z);

    // In case of error, return 1
    if (error)
    {
		return NAN;
    }

	// Get data from PT-flash
    FlashResults results = Flash::get_flash_results();
	X = results.Xij;
	nu = results.nuj;
    np = results.np;

    // Calculate total enthalpy of mixture
    XX = 0.;
	int j = 0;
	for (double nuj: nu)
	{
		if (nuj > 0.)
		{
			// Ideal and residual enthalpy/entropy
			// X = Σj Xj
			if (state_spec == StateSpecification::ENTHALPY)
			{
				XX += flash_params.eos_params[flash_params.eos_order[results.eos_idx[j]]].eos->H(p, T_, results.nij, j*flash_params.ns, true) * M_R;
			}
 			else  // ENTROPY
			{
				XX += flash_params.eos_params[flash_params.eos_order[results.eos_idx[j]]].eos->S(p, T_, results.nij, j*flash_params.ns, true) * M_R;
			}
		}
		j++;
	}

	return XX - X_spec;
}

double PXFlash::gradient(double T_)
{
	// Calculate gradient of specified variable with respect to temperature to apply gradient update of objective function
	(void) T_;
	FlashResults results = Flash::get_flash_results(true);
	
	// Calculate derivative of flash at solution with P, T(H_spec) with respect to temperature
	std::vector<double> dnudT(np), dxdT(np * flash_params.nc);
	results.dT_derivs(dnudT, dxdT);

	// Calculate derivative of specified property (H, S) with respect to temperature
	EoS::Property prop = (state_spec == StateSpecification::ENTHALPY) ? EoS::Property::ENTHALPY : EoS::Property::ENTROPY;
	return results.dXdT(prop, dnudT, dxdT);
}
