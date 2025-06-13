#include <cmath>
#include <numeric>
#include "dartsflash/global/timer.hpp"
#include "dartsflash/flash/flash_params.hpp"

EoSParams::EoSParams(CompData& comp_data)
{
	this->initial_guess = InitialGuess(comp_data);

	this->comp_idxs.resize(comp_data.ns);
	std::iota(comp_idxs.begin(), comp_idxs.end(), 0);

	this->nc = comp_data.nc;
	this->ni = comp_data.ni;
	this->ns = comp_data.ns;
}

void EoSParams::set_active_components(std::vector<int>& idxs)
{
	this->comp_idxs = idxs;

	int NC = this->nc;
	nc = 0; ni = 0;
	for (int i: idxs)
	{
		if (i < NC)
		{
			nc++;
		}
		else
		{
			ni++;
		}
	}
	this->ns = nc + ni;
	return;
}

FlashParams::FlashParams(CompData& compdata) : FlashParams()
{
	this->comp_data = compdata;
	this->nc = compdata.nc;
	this->ni = compdata.ni;
	this->ns = nc + ni;

	this->initial_guess = InitialGuess(compdata);
	this->eos_order = {};
}

void FlashParams::add_eos(std::string name, EoS* eos)
{
	eos_params[name] = EoSParams(this->comp_data);
	eos_params[name].eos = eos->getCopy();
	
	eos_order.push_back(name);
	return;
}

void FlashParams::init_eos(double p, double T)
{
	// Initialise EoS component parameters at p, T
    this->timer.start(Timer::timer::EOS);
	for (auto& it: this->eos_params) 
	{
		it.second.eos->set_root_flag(it.second.root_flag);
		it.second.eos->init_PT(p, T);
		it.second.initial_guess.init(p, T);
	}
	this->initial_guess.init(p, T);
    this->timer.stop(Timer::timer::EOS);
}

TrialPhase FlashParams::find_ref_comp(double p, double T, std::vector<double>& z)
{
    // Find reference compositions - hypothetical single phase
    double gmin = NAN;
    int ref_eos_idx = 0;
    EoS::RootFlag ref_root = EoS::RootFlag::STABLE;
    for (size_t j = 0; j < this->eos_order.size(); j++)
    {
		EoSParams eosparams = this->eos_params[this->eos_order[j]];
        if (eosparams.eos->eos_in_range(z.begin()))
        {
            eosparams.eos->set_root_flag(static_cast<EoS::RootFlag>(eosparams.root_flag));
            double gr = eosparams.eos->Gr(p, T, z, 0, true);

            if (eosparams.eos->select_root(z.begin()) >= EoS::RootSelect::ACCEPT  // root should be selected or not
                        && (std::isnan(gmin) // gmin not initialized
                        || (gr < gmin))) // Gibbs energy of EoS is lower
            {
				ref_root = eosparams.eos->is_root_type();
                ref_eos_idx = static_cast<int>(j);
                gmin = gr;
            }
        }
    }
    TrialPhase ref_comp = TrialPhase(ref_eos_idx, this->eos_order[ref_eos_idx], z);
    ref_comp.root = ref_root;

    if (this->verbose)
    {
        ref_comp.print_point("Reference phase");
    }
    return ref_comp;
}

std::vector<std::string> FlashParams::find_pure_phase(double p, double T)
{
	std::vector<std::string> eos_pure(nc);
	std::vector<double> G_pure(nc, NAN);

	for (auto it: this->eos_params)
	{
        it.second.eos->set_root_flag(it.second.root_flag);

		it.second.eos->init_PT(p, T);
		std::vector<double> gpure = it.second.eos->get_gpure();

		for (int i = 0; i < nc; i++)
		{
			if (!std::isnan(gpure[i]) && (std::isnan(G_pure[i]) || gpure[i] < G_pure[i]))
			{
				G_pure[i] = gpure[i];
				eos_pure[i] = it.first;
			}
		}
	}
	return eos_pure;
}

std::vector<double> FlashParams::G_pure(double p, double T)
{
	std::vector<std::string> eos_pure = this->find_pure_phase(p, T);
	std::vector<double> Gpure(nc, NAN);

	for (int i = 0; i < nc; i++)
	{
		EoSParams eosparams = this->eos_params[eos_pure[i]];
        eosparams.eos->set_root_flag(eosparams.root_flag);

		eosparams.eos->init_PT(p, T);
		Gpure[i] = eosparams.eos->get_gpure()[i];
	}

	return Gpure;
}

std::vector<double> FlashParams::H_pure(double p, double T)
{
	std::vector<std::string> eos_pure = this->find_pure_phase(p, T);
	std::vector<double> H_pure(nc, NAN);
	std::vector<double> n(nc, 0.);

	for (int i = 0; i < nc; i++)
	{
		EoSParams eosparams = this->eos_params[eos_pure[i]];
        eosparams.eos->set_root_flag(eosparams.root_flag);

		n[i] = 1.;
		eosparams.eos->init_PT(p, T);
		eosparams.eos->solve_PT(n.begin());
		H_pure[i] = eosparams.eos->H(p, T, n, 0, true);

		n[i] = 0.;
	}

	return H_pure;
}

std::vector<TrialPhase> EoSParams::evaluate_initial_guesses(int eos_idx, std::string eos_name, std::vector<TrialPhase>& ref_comps)
{
	std::vector<TrialPhase> trial_comps = this->initial_guess.evaluate(eos_idx, eos_name, this->initial_guesses, ref_comps);

	for (TrialPhase trial_comp: trial_comps)
	{
		trial_comp.root = static_cast<EoS::RootFlag>(this->root_flag);
	}
	return trial_comps;
}

std::vector<double> FlashParams::prop_pure(EoS::Property prop, double p, double T)
{
	switch (prop)
	{
		case EoS::Property::GIBBS:
		{
			return this->G_pure(p, T);
		}
		default:
		// case EoS::Property::ENTHALPY:
		{
			return this->H_pure(p, T);
		}
	}
}

std::vector<double> FlashParams::prop_1p(EoS::Property prop, double p, double T, std::vector<double>& X, std::vector<double>& eos_idxs, std::vector<double>& roots)
{
	std::vector<double> result(eos_idxs.size(), NAN);
	
	int j = 0;
	for (double eos_idx: eos_idxs)
	{
		if (std::isnan(eos_idx))
		{
			return result;
		}

		this->eos_params[this->eos_order[static_cast<int>(eos_idx)]].eos->set_root_flag(static_cast<EoS::RootFlag>(roots[j]));
		switch (prop)
		{
			case EoS::Property::GIBBS:
			{
				result[j] = this->eos_params[this->eos_order[static_cast<int>(eos_idx)]].eos->G(p, T, X, 0, true);
				break;
			}
			case EoS::Property::ENTHALPY:
			{
				result[j] = this->eos_params[this->eos_order[static_cast<int>(eos_idx)]].eos->H(p, T, X, 0, true);
				break;
			}
			default:
			{
				result[j] = this->eos_params[this->eos_order[static_cast<int>(eos_idx)]].eos->S(p, T, X, 0, true);
			}
		}
		
		j++;
	}
	return result;
}

std::vector<double> FlashParams::prop_np(EoS::Property prop, double p, double T, std::vector<double>& X, std::vector<double>& eos_idxs, std::vector<double>& roots)
{
	std::vector<double> result(eos_idxs.size(), NAN);
	
	int j = 0;
	for (double eos_idx: eos_idxs)
	{
		if (std::isnan(eos_idx))
		{
			return result;
		}

		if (!std::isnan(X[j*ns]))
		{
			this->eos_params[this->eos_order[static_cast<int>(eos_idx)]].eos->set_root_flag(static_cast<EoS::RootFlag>(roots[j]));
			switch (prop)
			{
				case EoS::Property::GIBBS:
				{
					result[j] = this->eos_params[this->eos_order[static_cast<int>(eos_idx)]].eos->G(p, T, X, j*ns, true);
					break;
				}
				case EoS::Property::ENTHALPY:
				{
					result[j] = this->eos_params[this->eos_order[static_cast<int>(eos_idx)]].eos->H(p, T, X, j*ns, true);
					break;
				}
				default:
				{
					result[j] = this->eos_params[this->eos_order[static_cast<int>(eos_idx)]].eos->S(p, T, X, j*ns, true);
				}
			}
		}
		
		j++;
	}
	return result;
}
