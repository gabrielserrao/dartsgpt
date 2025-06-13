#include <iostream>
#include <vector>
#include <cassert>
#include <numeric>
#include <algorithm>

#include "dartsflash/flash/flash.hpp"
#include <Eigen/Dense>

NegativeFlash::NegativeFlash(FlashParams& flashparams, const std::vector<std::string>& eos_used, const std::vector<int>& initial_guesses_) 
	: Flash(flashparams) 
{
	assert((eos_used.size() == (initial_guesses_.size() + 1)) && 
			"Length of specified EoS incompatible with number of initial guesses\n");
	this->np = static_cast<int>(eos_used.size());
	this->eos = eos_used;
	this->eos_idxs.resize(np);
	for (int j = 0; j < np; j++)
	{
		eos_idxs[j] = std::distance(flash_params.eos_order.begin(), std::find(flash_params.eos_order.begin(), flash_params.eos_order.end(), eos[j]));
	}
	this->initial_guesses = initial_guesses_;

	this->sp_idxs.resize(np);
	std::iota(this->sp_idxs.begin(), this->sp_idxs.end(), 0);
}

int NegativeFlash::evaluate(double p_, double T_, std::vector<double>& z_)
{
	// Solve N-phase split with negative flash

	// Initialize flash
	// Initialize EoS at p, T, correct composition vector and generate initial guess for two-phase split
	this->init(p_, T_, z_);

	// Initialize vector with TrialPhase for each phase
	this->stationary_points.clear();
    for (int j = 0; j < np; j++)
    {
        stationary_points.push_back(TrialPhase(eos_idxs[j], eos[j], z));
    }

	// Perform negative flash with eos specified
	int split_output = this->run_split(sp_idxs);

	// Return output
	if (split_output == 0)
	{
		// If split output is 0, output is correct
		return 0;
	}
	else if (split_output == -1)
	{
		// If split output is -1, one or more negative phases present
		(void) this->check_negative();
		if (flash_params.verbose)
        {
            print("NegativeFlash", "===============");
            print("p, T", std::vector<double>{p, T});
			print("z", z_);
            print("nu", nu);
            print("X", X, np);
        }
		return 0;
	}
	else
	{
		// Else, error occurred in split. Check if negative phases are present
		if (this->check_negative())
		{
			return 0;
		}
		// If not, phase split output is incorrect
		if (flash_params.verbose)
		{
			print("ERROR in NegativeFlash", split_output);
        	print("p, T", std::vector<double>{p, T}, 10);
			print("z", z_, 10);
		}
		for (TrialPhase comp: ref_compositions)
		{
			comp.nu = NAN;
			comp.X = std::vector<double>(flash_params.ns, NAN);
		}
		return 1;
	}
}

std::vector<double> NegativeFlash::generate_lnK(std::vector<int>& sp_idxs_)
{
	(void) sp_idxs_;
	std::vector<double> lnk = flash_params.initial_guess.evaluate(this->initial_guesses);
	std::vector<double> lnK = {};
	for (int i: flash_params.eos_params.begin()->second.comp_idxs)
	{
		lnK.push_back(lnk[i]);
	}
	return lnK;
}

bool NegativeFlash::check_negative() {
	if (ref_compositions[0].nu > 1.)
	{
		ref_compositions[0].nu = 1.;
		ref_compositions[1].nu = 0.;
		ref_compositions[0].X = z;
		ref_compositions[1].X = std::vector<double>(flash_params.ns, NAN);
		return true;
	}
	else if (ref_compositions[1].nu > 1.)
	{
		ref_compositions[1].nu = 1.;
		ref_compositions[0].nu = 0.;
		ref_compositions[1].X = z;
		ref_compositions[0].X = std::vector<double>(flash_params.ns, NAN);
		return true;
	}
	return false;
}	
