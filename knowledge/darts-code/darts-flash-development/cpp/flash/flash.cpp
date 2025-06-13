#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <map>

#include "dartsflash/maths/maths.hpp"
#include "dartsflash/maths/geometry.hpp"
#include "dartsflash/flash/flash.hpp"
#include "dartsflash/flash/flash_params.hpp"
#include "dartsflash/flash/trial_phase.hpp"
#include "dartsflash/stability/stability.hpp"
#include "dartsflash/phase-split/twophasesplit.hpp"
#include "dartsflash/phase-split/multiphasesplit.hpp"

Flash::Flash(FlashParams& flashparams) {
	this->flash_params = flashparams;

    this->z.resize(flashparams.ns);
    this->gpure.resize(flashparams.ns);
    // this->eos.reserve(NP_MAX);
    this->nu.reserve(NP_MAX);
    this->X.reserve(NP_MAX*flashparams.ns);

    this->stationary_points.reserve(NP_MAX);
    this->ref_compositions.reserve(NP_MAX);
    this->sp_idxs.reserve(NP_MAX + 1);
}

void Flash::init(double p_, double T_)
{
    // Set iterations to zero
    total_ssi_flash_iter = total_ssi_stability_iter = total_newton_flash_iter = total_newton_stability_iter = 0;
    
    // Initialize EoS and InitialGuess at p, T
    this->p = p_; this->T = T_;
    this->flash_params.init_eos(p, T);

    return;
}

void Flash::init(double p_, double T_, std::vector<double>& z_)
{
    Flash::init(p_, T_);

    // Check if feed composition needs to be corrected for 0 values
    z.resize(flash_params.ns);
    for (int i = 0; i < flash_params.ns; i++)
    {
        z[i] = (z_[i] > flash_params.min_z) ? z_[i] : flash_params.min_z;
    }
    return;
}

int Flash::evaluate(double p_, double T_)
{
    // Find reference compositions - hypothetical single phase
    Flash::init(p_, T_);
    this->np = 1;
    z = {1.};
    this->ref_compositions = {this->flash_params.find_ref_comp(p_, T_, z)};

    // Set nu and X for Results
    this->ref_compositions[0].nu = 1.;
    this->ref_compositions[0].X = z;
    return 0;
}

int Flash::evaluate(double p_, double T_, std::vector<double>& z_)
{
    // Evaluate sequential stability + flash algorithm
    
    // Initialize flash
    // Initialize EoS at p, T and check if feed composition needs to be corrected
    if (z_.size() == 1)
    {
        return Flash::evaluate(p_, T_);
    }
    Flash::init(p_, T_, z_);

    // Calculate pure component Gibbs free energies
    this->gpure = flash_params.G_pure(p, T);

    // Perform stability and phase split loop starting from np = 1
    // Run stability test + flash loop over phases until np_max has been reached
    ref_compositions = {flash_params.find_ref_comp(p, T, z)};
    ref_compositions[0].nu = 1.;
    ref_compositions[0].X = z;
    // eos = {ref_compositions[0].eos_name};
    this->np = 1;
    int it = 1;
    while (true)
    {
        int output = this->run_loop(ref_compositions);
        if (output == -1)
        {
            // Output -1, all phases stable
            if (flash_params.verbose)
            {
                print("StabilityFlash", "===============");
                print("p, T", std::vector<double>{p, T}, 1, 10);
			    print("z", z_, 1, 10);
                this->get_flash_results().print_results();
            }
            return 0;
        }
        else if (output > 0 || it > 10)
        {
            // Else, error occurred in split
            if (flash_params.verbose)
            {
                print("ERROR in StabilityFlash", output);
    		    print("p, T", std::vector<double>{p, T}, 1, 10);
			    print("z", z_, 1, 10);
            }
            return output;
        }
        it++;
    }

    return 0;
}

int Flash::run_stability(std::vector<TrialPhase>& ref_comps)
{
    // Run stability test on feed X with initial guesses Y
    Stability stab(flash_params);
    stab.init(ref_comps);

    // Add each of the minima trivial stationary points (for np > 1, this corresponds to all phase compositions)
    this->stationary_points.clear();
    EoS* eos0 = flash_params.eos_params[ref_comps[0].eos_name].eos;
    if (ref_comps.size() > 1 || eos0->is_convex(ref_comps[0].Y.begin()) || !eos0->has_multiple_minima())
    {
        this->stationary_points = ref_comps;
        if (flash_params.verbose)
        {
            std::cout << "Ref compositions are stationary points\n";
        }
    }
    sp_idxs = std::vector<int>(stationary_points.size());
    std::iota(sp_idxs.begin(), sp_idxs.end(), 0);

    // Iterate over initial guesses in Y to run stability tests
    for (size_t j = 0; j < this->flash_params.eos_order.size(); j++)
    {
        std::string eosname = this->flash_params.eos_order[j];
        std::vector<TrialPhase> trial_comps = this->flash_params.eos_params[eosname].evaluate_initial_guesses(j, eosname, ref_comps);
        for (TrialPhase trial_comp: trial_comps)
        {
            int error = stab.run(trial_comp);

            if (error > 0)
            {
                return error;
            }

            // Get TPD value and check if it is already in the list
            if (trial_comp.is_preferred_root >= EoS::RootSelect::ACCEPT && !this->compare_stationary_points(trial_comp))
            {
                // No duplicate found: add stationary point to vector of stationary points
                stationary_points.push_back(trial_comp);
            }

            // Get number of iterations from stability
            this->total_ssi_stability_iter += stab.get_ssi_iter();
            this->total_newton_stability_iter += stab.get_newton_iter();
        }
    }
    return 0;
}

int Flash::run_stability(std::vector<TrialPhase>& ref_comps, std::vector<int>& sp_idxs_)
{
    // Run stability test on flash compositions
    Stability stab(flash_params);
    stab.init(ref_comps);

    // Run stability tests starting from each stationary point
    for (auto it = stationary_points.begin(); it != stationary_points.end(); it++)
    {
        int error = stab.run(*it);

        if (error > 0)
        {
            return error;
        }

        // // Get TPD value and check if it is already in the list
        // if (it->is_preferred_root < EoS::RootSelect::ACCEPT)
        // {
        //     // No duplicate found: add stationary point to vector of stationary points
        //     // it->print_point("UNACcEPTABLE");
        //     // stationary_points.erase(it);
        // }
    }

    for (int sp_idx: sp_idxs_)
    {
        (void) sp_idx;
    //     // it.first: eos_name; it.second: EoSParams
    //     std::vector<TrialPhase> trial_comps = it.second.evaluate_initial_guesses(it.first, ref_comps);
    //     for (TrialPhase trial_comp: trial_comps)
    //     {
    //         int error = stab.run(trial_comp);

    //         if (error > 0)
    //         {
    //             return error;
    //         }

    //         // Get TPD value and check if it is already in the list
    //         if (trial_comp.is_preferred_root >= EoS::RootSelect::ACCEPT && !this->compare_stationary_points(trial_comp))
    //         {
    //             // No duplicate found: add stationary point to vector of stationary points
    //             stationary_points.push_back(trial_comp);
    //         }

    //         // Get number of iterations from stability
    //         this->total_ssi_stability_iter += stab.get_ssi_iter();
    //         this->total_newton_stability_iter += stab.get_newton_iter();
    //     }
    }
    return 0;
}

int Flash::run_split(std::vector<int>& sp_idxs_)
{
    // Choose proper set of lnK
    std::vector<double> lnK = this->generate_lnK(sp_idxs_);
    if (std::isnan(lnK[0])) { return -1; }

    std::vector<TrialPhase> trial_comps{};
    for (int sp_idx: sp_idxs_)
    {
        trial_comps.push_back(stationary_points[sp_idx]);
    }

    // Initialize split object
    int split_output = 0;
    np = static_cast<int>(sp_idxs_.size());
    if (np == 2)
    {
        // Initialize PhaseSplit object for two phases
        TwoPhaseSplit split(flash_params);

        // Run split and return nu and x
        split_output = split.run(z, lnK, trial_comps);  // Run multiphase split algorithm at P, T, z with initial guess lnK
        if (split_output == 1)
        {
            // Find a composition that is in the "middle" of the stationary points for easier solution
            split_output = split.run(this->z_mid, lnK, trial_comps);

            // Use converged solution as initial guess for actual composition
            std::vector<double> lnk = split.get_lnk();
            split_output = split.run(this->z, lnk, trial_comps);
        }
        gibbs = split.get_gibbs();
        total_ssi_flash_iter += split.get_ssi_iter();
        total_newton_flash_iter += split.get_newton_iter();
    }
    else
    {
        // Initialize MultiPhaseSplit object for 3 or more phases
        MultiPhaseSplit split(flash_params, np);

        // Run split and return nu and x
        split_output = split.run(z, lnK, trial_comps);  // Run multiphase split algorithm at P, T, z with initial guess lnK
        if (split_output == 1)
        {
            // Find a composition that is in the "middle" of the stationary points for easier solution
            split_output = split.run(this->z_mid, lnK, trial_comps);

            // Use converged solution as initial guess for actual composition
            std::vector<double> lnk = split.get_lnk();
            split_output = split.run(this->z, lnk, trial_comps);
        }
        gibbs = split.get_gibbs();
        total_ssi_flash_iter += split.get_ssi_iter();
        total_newton_flash_iter += split.get_newton_iter();
    }

    // Determine output value
    if (split_output <= 0)
    {
        ref_compositions = trial_comps;
    }
    return split_output;
}

int Flash::run_loop(std::vector<TrialPhase>& ref_comps)
{
    // Perform stability test starting from each of the initial guesses
    this->ref_compositions = ref_comps;
    int stability_output = this->run_stability(ref_compositions);
    if (stability_output > 0)
    {
        // Error occurred in stability test
        return stability_output;
    }
    for (auto sp = stationary_points.begin(); sp != stationary_points.end(); sp++)
    {
        sp->ymin = this->flash_params.eos_params[sp->eos_name].eos->mix_min(sp->y, this->gpure, sp->gmin);
    }

    // Count number of stationary points with negative TPD and sort stationary points according to tpd
    int negative_tpds = 0;
    for (size_t j = sp_idxs.size(); j < stationary_points.size(); j++)
    {
        TrialPhase* sp = &stationary_points[j];
        sp->ymin = this->flash_params.eos_params[sp->eos_name].eos->mix_min(sp->y, this->gpure, sp->gmin);
        if (stationary_points[j].tpd < -flash_params.tpd_tol)  // negative TPD
        {
            negative_tpds++;
        }
        // Find location of tpd in sorted idxs
        std::vector<int>::iterator it;
        for (it = sp_idxs.begin(); it != sp_idxs.end(); it++)
        {
            if (stationary_points[*it].tpd > stationary_points[j].tpd)
            {
                sp_idxs.insert(it, static_cast<int>(j));
                break;
            }
        }
        if (it == sp_idxs.end())
        {
            sp_idxs.push_back(static_cast<int>(j));
        }
    }

    // If for EoS gmix is preferred over stationary point, evaluate local minimum
    for (TrialPhase& sp: stationary_points)
    {
        if (flash_params.eos_params[sp.eos_name].use_gmix)
        {
            sp.ymin = this->flash_params.eos_params[sp.eos_name].eos->mix_min(sp.y, this->gpure, sp.gmin);
        }
    }

    // If no negative TPDs, feed is stable, return -1
    if (negative_tpds == 0 || stationary_points[sp_idxs[0]].tpd > -std::fabs(this->flash_params.tpd_1p_tol))
    {
        return -1;
    }
    else
    {
        // Else not stable
        if (flash_params.verbose)
        {
            print("Unstable feed", "============");
            for (TrialPhase stationary_point : stationary_points) { stationary_point.print_point(); }
        }

        // Determine lnK-values from set of stationary points for phase split
        int tot_sp = static_cast<int>(sp_idxs.size());
        int split_output = 0;

        if (tot_sp == 1)
        {
            if (flash_params.verbose)
            {
                std::cout << "Unstable feed, but only 1 stationary point in ref_compositions!\n";
                print("p, T", std::vector<double>{p, T});
                print("z", z);
            }
            // Add trivial solution, even though it is non-convex
            stationary_points.push_back(flash_params.find_ref_comp(p, T, z));
            sp_idxs.push_back(tot_sp);
            tot_sp++;
        }
        if (tot_sp == 2)  // TWO STATIONARY POINTS <= 0 --> 2 PHASES: Select both stationary points
        {
            // Find if composition is within simplex of stationary point compositions
            std::vector<std::vector<double>> coords = {};
            for (TrialPhase& sp: stationary_points)
            {
                coords.push_back(sp.y);
            }
            this->z_mid = find_midpoint(coords);

            // Run split with two stationary points
            split_output = this->run_split(sp_idxs);
            return split_output;
        }
        else  // In case of 3 or more stationary points <= 0, find possible combinations and determine which one has lowest Gibbs energy
        {
            // Maximum number of phases is equal to number of stationary points, but not larger than NC
            int NP_max = std::min(tot_sp, flash_params.nc);
            bool stable = false;
            double gibbs_lowest = NAN;
            std::vector<TrialPhase> comps;

            for (int NP = NP_max; NP >= 2; NP--)
            {
                // Find combinations of stationary points of maximum length np
                Combinations c(tot_sp, NP);
                
                // Find for each combination of stationary points if simplex with stationary points as vertices contains the feed composition
                // If it does, combination may be a solution to the multiphase equilibrium
                for (std::vector<int>& combination: c.combinations)
                {
                    bool run_combination = false;

                    std::vector<int> idxs = {};
                    bool close_to_boundary = true;  // if minimum tpd is too close to zero, use of lnK becomes problematic
                    for (int idx: combination)
                    {
                        int sp_idx = sp_idxs[idx];
                        idxs.push_back(sp_idx);
                        TrialPhase* sp = &stationary_points[sp_idx];

                        // Check if combination contains negative tpd
                        if (sp->tpd < 0.)
                        {
                            run_combination = true;
                        }
                        // Check if combination contains positive gmin, if so do not run combination and break from loop
                        if (flash_params.eos_params[sp->eos_name].use_gmix && sp->gmin > 0.)
                        {
                            run_combination = false;
                            break;
                        }
                        if (sp->tpd < -flash_params.tpd_close_to_boundary)
                        {
                            close_to_boundary = false;
                        }
                    }

                    if (run_combination && NP == flash_params.nc)
                    {
                        // Find if composition is within simplex of stationary point compositions
                        std::vector<std::vector<double>> coords = {};
                        for (int sp_idx: idxs)
                        {
                            coords.push_back(stationary_points[sp_idx].y);
                        }
                        this->z_mid = find_midpoint(coords);
                        if (!is_in_simplex(z, coords))
                        {
                            // Or if it is within minimum Gibbs energy compositions
                            coords = {};
                            for (int sp_idx: idxs)
                            {
                                coords.push_back(stationary_points[sp_idx].ymin);
                            }
                            if (!is_in_simplex(z, coords))
                            {
                                run_combination = false;
                            }
                        }
                    }

                    if (run_combination)
                    {
                        // If too close to phase boundary, change split variables to nik
                        FlashParams::SplitVars split_vars = this->flash_params.split_variables;
                        this->flash_params.split_variables = (close_to_boundary) ? FlashParams::SplitVars::nik : split_vars;

                        // Run phase split
                        split_output = this->run_split(idxs);

                        // Change back the variables
                        this->flash_params.split_variables = split_vars;

                        if (split_output == 0 && (std::isnan(gibbs_lowest) || gibbs < gibbs_lowest))
                        {
                            stability_output = this->run_stability(ref_compositions, idxs);

                            // Count number of stationary points with negative TPD and sort stationary points according to tpd
                            negative_tpds = 0;
                            for (auto it = stationary_points.begin(); it != stationary_points.end(); it++)
                            {
                                if (it->tpd < -flash_params.tpd_tol && it->is_preferred_root)
                                {
                                    negative_tpds++;
                                }
                            }

                            if (negative_tpds == 0)
                            {
                                // NP-Split success, found stable equilibrium of NP phases
                                ref_comps = ref_compositions;
                                return -1;
                            }

                            if (flash_params.verbose)
                            {
                                print("Unstable flash", "============");
                                for (TrialPhase stationary_point : stationary_points) { stationary_point.print_point(); }
                            }
                            gibbs_lowest = gibbs;
                            comps = ref_compositions;
                            stable = true;
                        }
                    }
                }
                if (stable)
                {
                    // NP-Split success, found equilibrium of NP phases
                    ref_comps = comps;
                    return 0;
                }
            }
            std::cout << "Error occurred in Flash::run_loop(), no stable flashes have been found\n";
            print("p, T", std::vector<double>{p, T});
			print("z", z);
            return 1;
        }
    }
}

std::vector<double> Flash::generate_lnK(std::vector<int>& sp_idxs_)
{
    // Determine lnK initialization for phase split
    np = static_cast<int>(sp_idxs_.size());
    int nc = flash_params.nc;

    std::vector<double> lnY0(nc);
    int sp_idx = sp_idxs_[0];
    TrialPhase* sp = &stationary_points[sp_idxs_[0]];
    if (flash_params.eos_params[sp->eos_name].use_gmix)
    {
        for (int i = 0; i < nc; i++)
        {
            lnY0[i] = std::log(sp->ymin[i]);
        }
    }
    else
    {
        for (int i = 0; i < nc; i++)
        {
            lnY0[i] = std::log(sp->y[i]);
        }
    }
    // eos = {stationary_points[sp_idx].eos_name};
    std::vector<TrialPhase> trial_comps = {stationary_points[sp_idx]};

    std::vector<double> lnK((np-1)*nc);
    for (size_t j = 1; j < sp_idxs_.size(); j++)
    {
        sp = &stationary_points[sp_idxs_[j]];

        // Determine composition to initialize K-values with
        if (flash_params.eos_params[sp->eos_name].use_gmix)
        {
            for (int i = 0; i < nc; i++)
            {
                lnK[(j-1) * nc + i] = std::log(sp->ymin[i]) - lnY0[i];
            }
        }
        else
        {
            for (int i = 0; i < nc; i++)
            {
                lnK[(j-1) * nc + i] = std::log(sp->y[i]) - lnY0[i];
            }
        }
        // eos.push_back(sp->eos_name);
        trial_comps.push_back(*sp);
    }
    return lnK;
}

bool Flash::compare_stationary_points(TrialPhase& stationary_point)
{
    // Compare stationary point with entries in vector of stationary points to check if it is unique
    // Returns true if point is already in the list
    double tpd0 = stationary_point.tpd;
    double lntpd0 = std::log(std::fabs(tpd0));
    for (size_t j = 0; j < stationary_points.size(); j++)
    {
        double tpdj = stationary_points[j].tpd;
        // For small tpd difference (tpd < 1), compare absolute difference; for large tpd values, logarithmic scale is used to compare
        double tpd_diff = lntpd0 < 0. ? tpdj-tpd0 : lntpd0 - std::log(std::fabs(tpdj) + 1e-15);
        if (stationary_points[j].eos_name == stationary_point.eos_name // eos is the same
             && (std::fabs(tpd_diff) < flash_params.tpd_tol || // tpd is within tolerance
                (std::fabs(tpd0) < flash_params.tpd_tol && std::fabs(tpdj) < flash_params.tpd_tol))) // both tpds are within absolute tpd tolerance
        {
            // Similar TPD found; Check if composition is also the same
            if (compare_compositions(stationary_point.Y, stationary_points[j].Y, flash_params.comp_tol))
            {
                if (stationary_point.root != stationary_points[j].root)
                {
                    stationary_points[j].root = EoS::RootFlag::STABLE;
                }
                return true;
            }
        }
    }
    return false;
}

std::vector<TrialPhase> Flash::find_stationary_points(double p_, double T_, std::vector<double>& X_)
{
    // Initialize EoS and InitialGuess at p, T
    Flash::init(p_, T_);
    
    // Perform stability test starting from each of the initial guesses
    std::vector<TrialPhase> ref_comps = {this->flash_params.find_ref_comp(p, T, X_)};
    int stability_output = this->run_stability(ref_comps);
    if (stability_output > 0)
    {
        // Error occurred in stability test
        if (this->flash_params.verbose)
		{
            print("ERROR in find_stationary_points()", stability_output);
		    print("p", p, 10);
            print("T", T, 10);

            for (TrialPhase ref: ref_comps)
            {
                ref.print_point();
            }
        }
        return ref_comps;
    }
    return this->stationary_points;
}
