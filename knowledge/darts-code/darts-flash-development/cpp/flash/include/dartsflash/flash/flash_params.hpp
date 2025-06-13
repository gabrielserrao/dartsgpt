//--------------------------------------------------------------------------
#ifndef OPENDARTS_FLASH_FLASH_FLASHPARAMS_H
#define OPENDARTS_FLASH_FLASH_FLASHPARAMS_H
//--------------------------------------------------------------------------

#include <vector>
#include <string>
#include <map>
#include "dartsflash/global/global.hpp"
#include "dartsflash/global/units.hpp"
#include "dartsflash/global/timer.hpp"
#include "dartsflash/eos/eos.hpp"
#include "dartsflash/flash/initial_guess.hpp"
#include "dartsflash/flash/trial_phase.hpp"

struct EoSParams
{
	EoSParams() {}
	EoSParams(CompData& comp_data);

	// EoS object
	EoS* eos;

	// Root selection and phase ordering
	EoS::RootFlag root_flag = EoS::RootFlag::STABLE;
	std::vector<EoS::RootFlag> root_order = {EoS::RootFlag::STABLE};
	std::vector<int> rich_phase_order = {};
	double rich_phase_composition = 0.5;
	
	// EoS related parameters
	double stability_tol{ 1e-20 };
	double stability_switch_tol = stability_tol;
	double stability_line_tol{ 1e-8 };
	int stability_max_iter{ 500 };
	int stability_line_iter{ 10 };

	// In-/exclude components from flash with constant salinity or (augmented) free-water flash
	int nc, ni, ns;
	std::vector<int> comp_idxs;
	void set_active_components(std::vector<int>& idxs);

	// Method to obtain initial guesses for this EoS
	InitialGuess initial_guess;
	std::vector<int> initial_guesses = {};
	std::vector<TrialPhase> evaluate_initial_guesses(int eos_idx, std::string eos_name, std::vector<TrialPhase>& ref_comps);
	bool use_gmix{ false };
};

struct FlashParams
{
	// Timers
	Timer timer;

	// Input/output units
	Units units; // default BAR, KELVIN, M3, KJ

	// Flash-related parameters
	double min_z{ 1e-13 };
	double y_pure{ 0.9 };
	double tpd_tol{ 1e-8 };
	double tpd_1p_tol{ 1e-4 };
	double tpd_close_to_boundary{ 1e-3 };
	double comp_tol{ 1e-4 };

	double rr2_tol{ 1e-12 };
	double rrn_tol{ 1e-14 };
	int rr_max_iter{ 100 };
	int rr_line_iter{ 10 };

	double split_tol{ 1e-15 };
	double split_switch_tol = split_tol;
	double split_line_tol{ 1e-8 };
	int split_max_iter{ 500 };
	int split_line_iter{ 100 };
	int split_negative_flash_iter = split_max_iter;
	double split_negative_flash_tol{ 1e-4 };
	enum SplitVars { nik = 0, lnK, lnK_chol };
	SplitVars split_variables = SplitVars::nik;
	bool modChol_split = true;

	enum StabilityVars { Y = 0, lnY, alpha };
	StabilityVars stability_variables = StabilityVars::Y;
	bool modChol_stab = true;

	// Ph-flash parameters. Temperature (K)
	double T_min{100};
	double T_max{1000};
	double T_init{300};

	double phflash_Htol{ 1e-3 };
	double phflash_Ttol{ 1e-8 };

	InitialGuess initial_guess;

	// verbose = flashpar::verbose::NONE;
	bool verbose = false;

	// EoS-related parameters
	int nc, ni, ns;
	CompData comp_data;
	std::unordered_map<std::string, EoSParams> eos_params = {};
	std::vector<std::string> eos_order;

	// Constructor
	FlashParams() {}
	FlashParams(CompData& comp_data);
	~FlashParams() = default;

	// Add EoS to EoSMap
	void add_eos(std::string name, EoS* eos);
	void init_eos(double p, double T);
	
	// EoS computations
	TrialPhase find_ref_comp(double p, double T, std::vector<double>& z);
	std::vector<std::string> find_pure_phase(double p, double T);
	std::vector<double> G_pure(double p, double T);
	std::vector<double> H_pure(double p, double T);

	// EoS property computations with flash results
	std::vector<double> prop_pure(EoS::Property prop, double p, double T);
	std::vector<double> prop_1p(EoS::Property prop, double p, double T, std::vector<double>& X, std::vector<double>& eos_idxs, std::vector<double>& roots);
	std::vector<double> prop_np(EoS::Property prop, double p, double T, std::vector<double>& X, std::vector<double>& eos_idxs, std::vector<double>& roots);
};

//--------------------------------------------------------------------------
#endif // OPENDARTS_FLASH_FLASH_FLASHPARAMS_H
//--------------------------------------------------------------------------
