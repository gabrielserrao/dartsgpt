//--------------------------------------------------------------------------
#ifndef OPENDARTS_FLASH_FLASH_FLASH_H
#define OPENDARTS_FLASH_FLASH_FLASH_H
//--------------------------------------------------------------------------

#include <vector>
#include "dartsflash/global/global.hpp"
#include "dartsflash/flash/flash_params.hpp"
#include "dartsflash/flash/flash_results.hpp"
#include "dartsflash/flash/trial_phase.hpp"

class Flash
{
public:
	enum class StateSpecification : int { TEMPERATURE = 0, ENTHALPY, ENTROPY };

protected:
	StateSpecification state_spec = StateSpecification::TEMPERATURE;
    int np;
	int total_ssi_flash_iter, total_ssi_stability_iter, total_newton_flash_iter, total_newton_stability_iter;

	double p, T;
	double gibbs;
	std::vector<double> z, nu, X, gpure, z_mid;
	std::vector<std::string> eos;
	std::vector<TrialPhase> ref_compositions, stationary_points;
	std::vector<int> sp_idxs;
	bool run_max_np;
	FlashParams flash_params;

public:
	Flash(FlashParams& flashparams);
	virtual ~Flash() = default;

	virtual int evaluate(double p_, double T_);  // single-component "flash"
	virtual int evaluate(double p_, double T_, std::vector<double>& z_);  // multicomponent flash

protected:
	virtual void init(double p_, double T_);
	virtual void init(double p_, double T_, std::vector<double>& z_);
	
	int run_stability(std::vector<TrialPhase>& ref_comps);
	int run_stability(std::vector<TrialPhase>& ref_comps, std::vector<int>& sp_idxs_);
	int run_split(std::vector<int>& sp_idxs_);
	int run_loop(std::vector<TrialPhase>& ref_comps);

	virtual std::vector<double> generate_lnK(std::vector<int>& sp_idxs_);
	bool compare_stationary_points(TrialPhase& stationary_point);

public:
	FlashResults get_flash_results(bool derivs = false) { return FlashResults(this->flash_params, this->p, this->T, this->z, this->ref_compositions, derivs); }

    int get_flash_total_ssi_iter(){return total_ssi_flash_iter;}
    int get_flash_total_newton_iter(){return total_newton_flash_iter;}
    int get_stability_total_ssi_iter(){return total_ssi_stability_iter;}
    int get_stability_total_newton_iter(){return total_ssi_stability_iter;}

	std::vector<TrialPhase> find_stationary_points(double p_, double T_, std::vector<double>& X_);
};

class NegativeFlash : public Flash
{
protected:
	std::vector<int> initial_guesses, eos_idxs;

public:
	NegativeFlash(FlashParams& flashparams, const std::vector<std::string>& eos_used, const std::vector<int>& initial_guesses_);

	virtual int evaluate(double p_, double T_, std::vector<double>& z_) override;

protected:
	virtual std::vector<double> generate_lnK(std::vector<int>& sp_idxs_) override;
	virtual bool check_negative();
};

//--------------------------------------------------------------------------
#endif // OPENDARTS_FLASH_FLASH_FLASH_H
//--------------------------------------------------------------------------
