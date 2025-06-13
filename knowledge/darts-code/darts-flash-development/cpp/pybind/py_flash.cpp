// #ifdef PYBIND11_ENABLED
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

#include "dartsflash/flash/flash.hpp"
#include "dartsflash/flash/px_flash.hpp"
#include "dartsflash/flash/initial_guess.hpp"
#include "dartsflash/flash/trial_phase.hpp"

namespace py = pybind11;

template <class FlashBase = Flash> class PyFlash : public FlashBase {
public:
    /* Inherit the constructors */
    using FlashBase::FlashBase;

    /* Trampoline (need one for each virtual function) */
    int evaluate(double p_, double T_, std::vector<double>& z_) override {
        PYBIND11_OVERRIDE(
            int,  /* Return type */
            Flash, /* Parent class */
            evaluate,   /* Name of function in C++ (must match Python name) */
            p_, T_, z_    /* Argument(s) */
        );
    }
    int evaluate(double p_, double T_) override {
        PYBIND11_OVERRIDE(
            int, Flash, evaluate, p_, T_
        );
    }
};

void pybind_flash(py::module& m) 
{
    using namespace pybind11::literals;  // bring in '_a' literal

    py::class_<Flash, PyFlash<>> flash(m, "Flash", R"pbdoc(
            This is a base class for Flash.

            Each Flash child class overrides the methods for `evaluate(p, T, z)`.
            )pbdoc");
    
    flash.def(py::init<FlashParams&>(), R"pbdoc(
            This is the constructor of the Flash base class.

            :param flashparams: Flash parameters object
            :type flashparams: FlashParams
            )pbdoc", py::arg("flashparams"))
        
        .def("evaluate", py::overload_cast<double, double>(&Flash::evaluate), R"pbdoc(
            Evaluate single-component equilibrium at (P, T)

            :param p: Pressure
            :type p: double
            :param T: Temperature
            :type T: double
            :returns: Error output of flash procedure
            :rtype: int
            )pbdoc", "p"_a, "T"_a)
        .def("evaluate", py::overload_cast<double, double, std::vector<double>&>(&Flash::evaluate), R"pbdoc(
            Evaluate multicomponent flash at (P, T, z)

            :param p: Pressure
            :type p: double
            :param T: Temperature
            :type T: double
            :param z: Feed composition
            :type z: list
            :returns: Error output of flash procedure
            :rtype: int
            )pbdoc", "p"_a, "T"_a, "z"_a)
        
        .def("get_flash_results", &Flash::get_flash_results, R"pbdoc(
            :returns: Results of flash
            :rtype: FlashResults
            )pbdoc", py::arg("derivs")=false)
        // using array = py::array_t<double, pybind11::array::c_style | pybind11::array::forcecast>;

        .def("find_stationary_points", &Flash::find_stationary_points, R"pbdoc(
            Determine stationary points of TPD function at (multiphase) composition X

            :param p: Pressure
            :type p: double
            :param T: Temperature
            :type T: double
            :param X: List of reference compositions
            :type X: list
            
            :returns: Set of stationary points
            :rtype: list
            )pbdoc", "p"_a, "T"_a, "X"_a)
        ;

    py::enum_<Flash::StateSpecification>(flash, "StateSpecification", "State specification for P-based flash (PT/PH/PS)")
        .value("TEMPERATURE", Flash::StateSpecification::TEMPERATURE)
        .value("ENTHALPY", Flash::StateSpecification::ENTHALPY)
        .value("ENTROPY", Flash::StateSpecification::ENTROPY)
        .export_values()
        ;

    py::class_<FlashResults>(m, "FlashResults", R"pbdoc(
            This is a struct containing flash results and derivatives.
            )pbdoc")
        .def("print_results", &FlashResults::print_results)
        .def_readonly("pressure", &FlashResults::pressure)
        .def_readonly("temperature", &FlashResults::temperature)
        .def_readonly("np", &FlashResults::np)
        .def_readonly("nu", &FlashResults::nuj)
        .def_readonly("X", &FlashResults::Xij)
        .def_readonly("eos_idx", &FlashResults::eos_idx)
        .def_readonly("root_type", &FlashResults::root_type)
        ;

    py::class_<NegativeFlash, PyFlash<NegativeFlash>, Flash>(m, "NegativeFlash", R"pbdoc(
            This is the two-phase implementation of negative flash.

            It evaluates two-phase negative flash and returns single-phase composition if either phase fraction is negative.
            )pbdoc")
        .def(py::init<FlashParams&, const std::vector<std::string>&, const std::vector<int>&>(), R"pbdoc(
            This is the constructor of the NegativeFlash class.

            :param flashparams: Flash parameters object
            :type flashparams: FlashParams
            :param eos_used: List of EoS names
            :type eos_used: list
            :param initial_guesses: List of initial guesses for K-values
            :type initial_guesses: list
            )pbdoc", "flashparams"_a, "eos_used"_a, "initial_guesses"_a)
        ;

    py::class_<PXFlash, PyFlash<PXFlash>, Flash>(m, "PXFlash", R"pbdoc(
            This is the implementation of P-based flash (PH/PS).

            )pbdoc")

        .def(py::init<FlashParams&, Flash::StateSpecification>(), R"pbdoc(
            This is the constructor of the PXFlash class.

            :param flashparams: Flash parameters object
            :type flashparams: FlashParams
            :param state_spec: State specification for P-based flash (ENTHALPY/ENTROPY)
            :type state_spec: Flash.StateSpecification
            )pbdoc", "flashparams"_a, "state_spec"_a)
        
        .def("evaluate", py::overload_cast<double, double>(&PXFlash::evaluate), R"pbdoc(
            Evaluate single-component equilibrium at (P, X)

            :param p: Pressure
            :type p: double
            :param X: Enthalpy/entropy
            :type X: double
            :returns: Error output of flash procedure
            :rtype: int
            )pbdoc", "p"_a, "X"_a)
        .def("evaluate", py::overload_cast<double, double, std::vector<double>&>(&PXFlash::evaluate), R"pbdoc(
            Evaluate multicomponent flash at (P, X, z)

            :param p: Pressure
            :type p: double
            :param h: Enthalpy/entropy
            :type h: double
            :param z: Feed composition
            :type z: list
            :returns: Error output of flash procedure
            :rtype: int
            )pbdoc", "p"_a, "X"_a, "z"_a)

        .def("evaluate_PT", py::overload_cast<double, double>(&PXFlash::evaluate_PT), R"pbdoc(
            Evaluate single-component equilibrium at (P, T)

            :param p: Pressure
            :type p: double
            :param T: Temperature
            :type T: double
            :returns: Error output of flash procedure
            :rtype: int
            )pbdoc", "p"_a, "T"_a)
        .def("evaluate_PT", py::overload_cast<double, double, std::vector<double>&>(&PXFlash::evaluate_PT), R"pbdoc(
            Evaluate multicomponent flash at (P, T, z)

            :param p: Pressure
            :type p: double
            :param T: Temperature
            :type T: double
            :param z: Feed composition
            :type z: list
            :returns: Error output of flash procedure
            :rtype: int
            )pbdoc", "p"_a, "T"_a, "z"_a)
        ;

    py::class_<FlashParams> flash_params(m, "FlashParams", R"pbdoc(
            This is a class that contains all required parameters that can be passed to a flash algorithm.
            )pbdoc");

    flash_params.def(py::init<CompData&>(), R"pbdoc(
            :param comp_data: Component data object
            :type comp_data: CompData
            )pbdoc", "comp_data"_a)

        .def_readwrite("timer", &FlashParams::timer, "Timer object")
        .def_readwrite("units", &FlashParams::units, "Units object")

        // Flash-related parameters
        .def_readwrite("min_z", &FlashParams::min_z, "Minimum value for composition")
        .def_readwrite("y_pure", &FlashParams::y_pure, "Mole fraction for preferred EoS range")

        .def_readwrite("rr2_tol", &FlashParams::rr2_tol, "Tolerance for two-phase Rachford-Rice norm")
        .def_readwrite("rrn_tol", &FlashParams::rrn_tol, "Tolerance for N-phase Rachford-Rice norm")
        .def_readwrite("rr_max_iter", &FlashParams::rr_max_iter, "Maximum number of iterations for Rachford-Rice procedure")

        .def_readwrite("split_tol", &FlashParams::split_tol, "Tolerance for phase split norm")
        .def_readwrite("split_switch_tol", &FlashParams::split_switch_tol, "Tolerance for switch to Newton in phase split")
        .def_readwrite("split_line_tol", &FlashParams::split_line_tol, "Tolerance for line search in phase split")
        .def_readwrite("split_max_iter", &FlashParams::split_max_iter, "Maximum number of iterations for phase split")
        .def_readwrite("split_line_iter", &FlashParams::split_line_iter, "Maximum number of iterations for line search in phase split")
        .def_readwrite("split_negative_flash_iter", &FlashParams::split_negative_flash_iter, "Number of iterations to quit split if negative flash")
        .def_readwrite("split_variables", &FlashParams::split_variables, "Variables for phase split: 0) n_ik, 1) lnK, 2) lnK-chol")
        
        .def_readwrite("stability_variables", &FlashParams::stability_variables, "Variables for stability: 0) Y, 1) lnY, 2) alpha")
        .def_readwrite("tpd_tol", &FlashParams::tpd_tol, "Tolerance for comparing tpd; also used to determine limit for tpd that is considered stable")
        .def_readwrite("tpd_1p_tol", &FlashParams::tpd_1p_tol, "Tolerance to determine limit for tpd that is considered stable")
        .def_readwrite("tpd_close_to_boundary", &FlashParams::tpd_close_to_boundary, "Tolerance to check if too close to phase boundary and switch PhaseSplit variables")
        .def_readwrite("comp_tol", &FlashParams::comp_tol, "Tolerance for comparing compositions")

        .def_readwrite("T_min", &FlashParams::T_min, "Minimum temperature for PXFlash")
        .def_readwrite("T_max", &FlashParams::T_max, "Maximum temperature for PXFlash")
        .def_readwrite("T_init", &FlashParams::T_init, "Initial temperature for PXFlash")
        .def_readwrite("phflash_Htol", &FlashParams::phflash_Htol, "Tolerance for comparing enthalpies in PXFlash")
        .def_readwrite("phflash_Ttol", &FlashParams::phflash_Ttol, "Tolerance for comparing temperatures in PXFlash")
        
        .def_readwrite("verbose", &FlashParams::verbose, "Verbose level")

        // EoS-related parameters
        .def_readwrite("eos_params", &FlashParams::eos_params, "Map of EoSParams object associated with each EoS object")
        .def_readwrite("eos_order", &FlashParams::eos_order, "Order of EoS for output")
        .def("add_eos", &FlashParams::add_eos, R"pbdoc(
            Add EoSParams object to map. This function creates a copy of the EoS object inside the EoSParams struct.

            :param name: EoS name
            :type name: str
            :param eos: EoS object
            :type eos: EoS
            )pbdoc", "name"_a, "eos"_a)
        .def("init_eos", &FlashParams::init_eos, R"pbdoc(
            Initialize EoS parameters at (P,T)

            :param p: Pressure
            :type p: double
            :param T: Temperature
            :type T: double
            )pbdoc", "p"_a, "T"_a)
        .def("find_ref_comp", &FlashParams::find_ref_comp, R"pbdoc(
            Find stable phase at state (P,T,n)

            :param p: Pressure
            :type p: float
            :param T: Temperature
            :type T: float
            :param n: List of compositions
            :type n: list
            )pbdoc", "p"_a, "T"_a, "n"_a)
        .def("find_pure_phase", &FlashParams::find_pure_phase, R"pbdoc(
            Find pure phases at state (P,T)

            :param p: Pressure
            :type p: float
            :param T: Temperature
            :type T: float
            )pbdoc", "p"_a, "T"_a)
        .def("G_pure", &FlashParams::G_pure)
        .def("H_pure", &FlashParams::H_pure)
        .def("prop_pure", &FlashParams::prop_pure)
        .def("prop_1p", &FlashParams::prop_1p)
        .def("prop_np", &FlashParams::prop_np)
        ;

    py::enum_<FlashParams::SplitVars>(flash_params, "SplitVars", "Primary variables for phase split")
        .value("nik", FlashParams::SplitVars::nik)
        .value("lnK", FlashParams::SplitVars::lnK)
        .value("lnK_chol", FlashParams::SplitVars::lnK_chol)
        .export_values()
        ;

    py::enum_<FlashParams::StabilityVars>(flash_params, "StabilityVars", "Primary variables for stability")
        .value("Y", FlashParams::StabilityVars::Y)
        .value("lnY", FlashParams::StabilityVars::lnY)
        .value("alpha", FlashParams::StabilityVars::alpha)
        .export_values()
        ;

    py::class_<EoSParams> eos_params(m, "EoSParams", R"pbdoc(
            This is a class that contains all required flash/stability parameters for each EoS object.
            )pbdoc");
    
    eos_params.def_readwrite("eos", &EoSParams::eos, "EoS object")
        .def_readwrite("root_flag", &EoSParams::root_flag, "EoS::RootFlag enum for roots to be selected in stability test. 0) STABLE, 1) MIN, 2) MAX")
        .def_readwrite("root_order", &EoSParams::root_order, "Order of EoS::RootFlag types in FlashResults, default is STABLE (unordered)")
        .def_readwrite("rich_phase_order", &EoSParams::rich_phase_order, "Order of rich phases component idxs in FlashResults, default is empty (unordered)")
        .def_readwrite("rich_phase_composition", &EoSParams::rich_phase_composition, "Composition to be considered rich phase, default is 0.5")
        .def_readwrite("use_gmix", &EoSParams::use_gmix, "Flag to use minimum of gmix for phase split rather than stationary point")
        .def_readwrite("initial_guesses", &EoSParams::initial_guesses, "Set of initial guesses for EoS")
        .def_readwrite("stability_tol", &EoSParams::stability_tol, "Tolerance for stability norm")
        .def_readwrite("stability_switch_tol", &EoSParams::stability_switch_tol, "Tolerance for switch to Newton in stability")
        .def_readwrite("stability_line_tol", &EoSParams::stability_line_tol, "Tolerance for line search in stability")
        .def_readwrite("stability_max_iter", &EoSParams::stability_max_iter, "Maximum number of iterations for stability")
        .def_readwrite("stability_line_iter", &EoSParams::stability_line_iter, "Maximum number of iterations for line search in stability")
        
        .def("set_active_components", &EoSParams::set_active_components, R"pbdoc(
            :param idxs: List of active component indices for EoS object
            :type idxs: list
            )pbdoc")
        ;

    // Expose TrialPhase and InitialGuess classes
    py::class_<TrialPhase>(m, "TrialPhase", R"pbdoc(
            This is a struct containing all information for a trial phase
            )pbdoc")
        .def(py::init<int, std::string, std::vector<double>&>(), R"pbdoc(
            :param eos_idx: Index of EoS in FlashParams::eos_order
            :param eos_name: Name of EoS of trial phase
            :param Y: Trial phase composition Y
            )pbdoc", "eos_idx"_a, "eos_name"_a, "Y"_a)

        .def_readonly("Y", &TrialPhase::Y)
        .def_readonly("y", &TrialPhase::y)
        .def_readonly("tpd", &TrialPhase::tpd)
        .def_readonly("eos_name", &TrialPhase::eos_name)
        .def_readonly("eos_idx", &TrialPhase::eos_idx)
        .def_readonly("root_type", &TrialPhase::root)

        .def("print_point", &TrialPhase::print_point, R"pbdoc(
            Function to print out trial phase information
            )pbdoc")
        ;

    py::class_<InitialGuess> ig(m, "InitialGuess", R"pbdoc(
            This is a base class for providing initial guesses to Flash, Stability or PhaseSplit algorithms.

            It contains methods for evaluating Wilson's (1968) correlation for vapour-liquid, Henry's law coefficients (Sander, 2015)
            vapour-hydrate and `pure phase` ideal K-values.
            )pbdoc");
    
    ig.def(py::init<CompData&>(), R"pbdoc(
            This is the constructor for generating initial guesses for Stability algorithms

            :param comp_data: Component data
            :type comp_data: CompData
            )pbdoc", "comp_data"_a)
        
        .def("init", &InitialGuess::init, R"pbdoc(
            Initialize P and T inside InitialGuess object.

            :param p: Pressure
            :type p: double
            :param T: Temperature
            :type T: double
            )pbdoc", "p"_a, "T"_a)
        
        .def("wilson", &InitialGuess::k_wilson, R"pbdoc(        
            :returns: List of Wilson (1968) K-values at P-T
            :rtype: list
            )pbdoc")
        .def("henry", &InitialGuess::k_henry, R"pbdoc(        
            :returns: List of Henry's law K-values at P-T, Sander (2015)
            :rtype: list
            )pbdoc")
        .def("y_henry", &InitialGuess::y_henry, R"pbdoc(        
            :returns: List of Henry's law K-values at P-T, Sander (2015)
            :rtype: list
            )pbdoc")
        .def("vapour_sI", &InitialGuess::k_vapour_sI, R"pbdoc(        
            :returns: List of vapour-sI K-values at P-T
            :rtype: list
            )pbdoc")
        .def("vapour_sII", &InitialGuess::k_vapour_sII, R"pbdoc(        
            :returns: List of vapour-sII K-values at P-T
            :rtype: list
            )pbdoc")
        .def("y_pure", &InitialGuess::y_pure, R"pbdoc(        
            :param j: j'th component
            :type j: int
            :param pure: Mole fraction of pure component, default is 0.9
            :type pure: double
            :returns: List of pure component K-values at P-T-x
            :rtype: list
            )pbdoc", "j"_a)
        .def("set_ypure", &InitialGuess::set_ypure, "Set mole fraction for pure phase guess of component j", "j"_a, "pure"_a)
        ;

    py::enum_<InitialGuess::Ki>(ig, "Ki", "Correlation type for phase split initial guess (K)")
        .value("Wilson_VL", InitialGuess::Ki::Wilson_VL)
        .value("Wilson_LV", InitialGuess::Ki::Wilson_LV)
        .value("Henry_VA", InitialGuess::Ki::Henry_VA)
        .value("Henry_AV", InitialGuess::Ki::Henry_AV)
        .export_values()
        ;
    
    py::enum_<InitialGuess::Yi>(ig, "Yi", "Correlation type for stability test initial guess (Y)")
        .value("Wilson", InitialGuess::Yi::Wilson)
        .value("Wilson13", InitialGuess::Yi::Wilson13)
        .value("Henry", InitialGuess::Yi::Henry)
        .value("sI", InitialGuess::Yi::sI)
        .value("sII", InitialGuess::Yi::sII)
        .value("sH", InitialGuess::Yi::sH)
        .value("Pure", InitialGuess::Yi::Pure)
        .export_values()
        ;
};

// #endif //PYBIND11_ENABLED