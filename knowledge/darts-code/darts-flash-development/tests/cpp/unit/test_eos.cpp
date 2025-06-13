#include <algorithm>
#include <numeric>
#include "dartsflash/eos/eos.hpp"
#include "dartsflash/eos/ideal.hpp"
#include "dartsflash/eos/aq/aq.hpp"
#include "dartsflash/eos/aq/jager.hpp"
#include "dartsflash/eos/aq/ziabakhsh.hpp"
#include "dartsflash/eos/helmholtz/cubic.hpp"
// #include "dartsflash/eos/helmholtz/gerg.hpp"
#include "dartsflash/eos/iapws/iapws95.hpp"
#include "dartsflash/eos/iapws/iapws_ice.hpp"
#include "dartsflash/eos/solid/solid.hpp"
#include "dartsflash/eos/vdwp/ballard.hpp"
#include "dartsflash/eos/vdwp/munck.hpp"
#include "dartsflash/flash/flash_params.hpp"

int test_ideal();
int test_cubic();
int test_aq();
int test_iapws();
// int test_gerg();
int test_solid();
int test_vdwp();
int test_gmix();

int gmix_test(EoS* eos, double p, double T, std::vector<double>& x, std::vector<double>& gpure, bool verbose = false)
{
	// Test for finding local minimum of Gibbs energy of mixing surface
	double gmix = 0.;
	eos->init_PT(p, T);
	std::vector<double> xmin = eos->mix_min(x, gpure, gmix);

	if (std::isnan(xmin[0]))
	{
		std::cout << "Minimum of gmix not found\n";
		print("p, T", std::vector<double>{p, T});
		print("X", x);
		print("xmin", xmin);
		return 1;
	}
	else
	{
		if (verbose)
		{
			std::cout << "Local minimum of gmix\n";
			print("p, T", std::vector<double>{p, T});
			print("x", x);
			print("xmin", xmin);
		}
		return 0;
	}
}

int main() 
{
    int error_output = 0;

	error_output += test_ideal();
	error_output += test_cubic();
	error_output += test_iapws();
	error_output += test_aq();
	error_output += test_solid();
	error_output += test_vdwp();
	error_output += test_gmix();

    return error_output;
}

int test_ideal()
{
	// Test implementation of ideal gas enthalpy and derivative w.r.t. T
	int error_output = 0;
	double d;

	std::vector<std::string> comp = {"C1"};
	CompData comp_data(comp);
	comp_data.cpi = {comp_data::cpi["C1"]};

	IdealGas ig(comp_data);
	std::vector<double> x{ 1. };
	double T = 300.;
	double dT = 1e-5;
	double tol = 1e-5;

	// Enthalpy and heat capacity test
	double H0 = ig.hi(T, 0);
	double H1 = ig.hi(T + dT, 0);
	double CP = ig.cpi(T, 0);
	double CP_num = (H1-H0)/dT;
	d = std::log(std::fabs(CP + 1e-15)) - std::log(std::fabs(CP_num + 1e-15));
	if (std::fabs(d) > tol) { print("dH/dT != dH/dT", std::vector<double>{CP, CP_num, d}); error_output++; }

	if (error_output > 0)
	{
		print("Errors occurred in test_ideal()", error_output);
	}
	else
	{
		print("No errors occurred in test_ideal()", error_output);
	}
	return error_output;
}

int test_cubic()
{
	// Test Cubics (Peng-Robinson and Soave-Redlich-Kwong)
	int error, error_output = 0;
	bool verbose = false;
	double p{ 30. }, v, T{ 300. };

	// Pure CO2
	std::vector<std::string> comp = {"CO2"};
	CompData comp_data(comp);
	comp_data.Pc = {73.75};
    comp_data.Tc = {304.1};
    comp_data.ac = {0.239};
    comp_data.kij = {0.};
	std::vector<double> n = {1.};

	CubicEoS pr(comp_data, CubicEoS::PR);
    CubicEoS srk(comp_data, CubicEoS::SRK);

    // Consistency tests for cubics
	v = pr.V(p, T, n);
	error = 0;
	error += pr.dlnphi_test(p, T, n, 1e-3);
	error += pr.derivatives_test(v, T, n, 2e-4);
	error += pr.lnphi_test(p, T, n, 2e-4);
	error += pr.pressure_test(p, T, n, 1e-5);
	// error += pr.temperature_test(p, T, n, 1e-5);
	error += pr.composition_test(p, T, n, 1e-5);
	error += pr.pvt_test(p, T, n, 1e-5);
	error += pr.critical_point_test(n, 1e4);
	error += pr.properties_test(p, T, n, 1e-5);
	error += pr.mix_dT_test(T, n, 1e-4);
	if (error || verbose) { print("PR.test() CO2", error); error_output += error; }

	v = srk.V(p, T, n);
	error = 0;
	error += srk.dlnphi_test(p, T, n, 1e-3);
	error += srk.derivatives_test(v, T, n, 2e-4);
	error += srk.lnphi_test(p, T, n, 3e-4);
	error += srk.pressure_test(p, T, n, 1e-5);
	// error += srk.temperature_test(p, T, n, 1e-5);
	error += srk.composition_test(p, T, n, 1e-5);
	error += srk.pvt_test(p, T, n, 1e-5);
	error += srk.critical_point_test(n, 1e4);
	error += srk.properties_test(p, T, n, 1e-5);
	error += srk.mix_dT_test(T, n, 1e-4);
	if (error || verbose) { print("SRK.test() CO2", error); error_output += error; }

	// MY10 mixture
	//// Sour gas mixture, data from (Li, 2012)
	comp = std::vector<std::string>{"CO2", "N2", "H2S", "C1", "C2", "C3"};
    comp_data = CompData(comp);
	comp_data.Pc = {73.819, 33.9, 89.4, 45.992, 48.718, 42.462};
    comp_data.Tc = {304.211, 126.2, 373.2, 190.564, 305.322, 369.825};
    comp_data.ac = {0.225, 0.039, 0.081, 0.01141, 0.10574, 0.15813};
    comp_data.kij = std::vector<double>(6*6, 0.);
    comp_data.set_binary_coefficients(0, {0., -0.02, 0.12, 0.125, 0.135, 0.150});
    comp_data.set_binary_coefficients(1, {-0.02, 0., 0.2, 0.031, 0.042, 0.091});
    comp_data.set_binary_coefficients(2, {0.12, 0.2, 0., 0.1, 0.08, 0.08});
	n = {0.9, 0.03, 0.04, 0.06, 0.04, 0.03};

	pr = CubicEoS(comp_data, CubicEoS::PR);
    srk = CubicEoS(comp_data, CubicEoS::SRK);
	
    // // Consistency tests for cubics
	v = pr.V(p, T, n);
	error = 0;
	error += pr.dlnphi_test(p, T, n, 1e-3);
	error += pr.derivatives_test(v, T, n, 2e-4);
	error += pr.lnphi_test(p, T, n, 2e-4);
	error += pr.pressure_test(p, T, n, 1e-5);
	// error += pr.temperature_test(p, T, n, 1e-5);
	error += pr.composition_test(p, T, n, 1e-5);
	error += pr.pvt_test(p, T, n, 1e-5);
	error += pr.critical_point_test(n, 1e3);
	error += pr.properties_test(p, T, n, 3e-3);
	error += pr.mix_dT_test(T, n, 1e-4);
	if (error || verbose) { print("PR.test() MY10", error); error_output += error; }

	v = srk.V(p, T, n);
	error = 0;
	error += srk.dlnphi_test(p, T, n, 1e-3);
	error += srk.derivatives_test(v, T, n, 2e-4);
	error += srk.lnphi_test(p, T, n, 1e-3);
	error += srk.pressure_test(p, T, n, 1e-5);
	// error += srk.temperature_test(p, T, n, 1e-5);
	error += srk.composition_test(p, T, n, 1e-5);
	error += srk.pvt_test(p, T, n, 1e-5);
	error += srk.critical_point_test(n, 1e3);
	error += srk.properties_test(p, T, n, 4e-3);
	error += srk.mix_dT_test(T, n, 1e-4);
	if (error || verbose) { print("SRK.test() MY10", error); error_output += error; }

	// H2O-CO2 mixture
	comp = {"H2O", "CO2"};
    comp_data = CompData(comp);
    comp_data.Pc = {220.50, 73.75};
    comp_data.Tc = {647.14, 304.10};
    comp_data.ac = {0.328, 0.239};
    comp_data.kij = std::vector<double>(2*2, 0.);
	comp_data.set_binary_coefficients(0, {0., 0.19014});
	n = {1., 0.};

	pr = CubicEoS(comp_data, CubicEoS::PR);
	srk = CubicEoS(comp_data, CubicEoS::SRK);
	error = pr.pvt_test(1., 360., n, 1e-4);
	if (error || verbose) { print("PR.pvt_test()", error); error_output += error; }
	error = pr.critical_point_test(n, 1e0);
	if (error || verbose) { print("PR.critical_point_test()", error); error_output += error; }
	error = srk.pvt_test(1., 360., n, 1e-4);
	if (error || verbose) { print("SRK.pvt_test()", error); error_output += error; }
	error = srk.critical_point_test(n, 1e0);
	if (error || verbose) { print("SRK.critical_point_test()", error); error_output += error; }
	

	if (error_output > 0)
	{
		print("Errors occurred in test_cubic()", error_output);
	}
	else
	{
		print("No errors occurred in test_cubic()", error_output);
	}
	return error_output;
}

int test_aq()
{
    int error, error_output = 0;
	bool verbose = false;

	// Without ions
    std::vector<std::string> comp{"H2O", "CO2", "C1"};
	CompData comp_data(comp);

	std::vector<double> n{0.95, 0.045, 0.005};
	double p{20.}, T{300.};

	Ziabakhsh2012 zia(comp_data);
	error = zia.derivatives_test(p, T, n, 1e-5);
	if (error || verbose) { print("Ziabakhsh.derivatives_test()", error); error_output += error; }

	Jager2003 jag(comp_data);
	error = jag.derivatives_test(p, T, n, 1e-5);
	if (error || verbose) { print("Jager.derivatives_test()", error); error_output += error; }

	// Test AQEoS
	// Jager (2003) for water and solutes
	AQEoS aq(comp_data, AQEoS::Model::Jager2003);
	error = aq.dlnphi_test(p, T, n, 5e-3);
	if (error || verbose) { print("Jager.dlnphi_test()", error); error_output += error; }
	error = aq.properties_test(p, T, n, 2e-3);
	if (error || verbose) { print("Jager.properties_test()", error); error_output += error; }

	// Ziabakhsh (2012) for water and solutes
	aq = AQEoS(comp_data, AQEoS::Model::Ziabakhsh2012);
	error = aq.dlnphi_test(p, T, n, 1e-4);
	if (error || verbose) { print("Ziabakhsh.dlnphi_test()", error); error_output += error; }
	error = aq.properties_test(p, T, n, 2e-3);
	if (error || verbose) { print("Ziabakhsh.properties_test()", error); error_output += error; }

	// Test mixed AQEoS models: Jager (2003) for water and Ziabakhsh (2012) for solutes
	// Test passing evaluator_map with [CompType, Model]
	std::map<AQEoS::CompType, AQEoS::Model> evaluator_map = {
		{AQEoS::CompType::water, AQEoS::Model::Jager2003},
		{AQEoS::CompType::solute, AQEoS::Model::Ziabakhsh2012},
		{AQEoS::CompType::ion, AQEoS::Model::Jager2003}
	};
	aq = AQEoS(comp_data, evaluator_map);
	error = aq.dlnphi_test(p, T, n, 5e-3);
	if (error || verbose) { print("AQEoS.dlnphi_test()", error); error_output += error; }
	error = aq.properties_test(p, T, n, 2e-3);
	if (error || verbose) { print("AQEoS.properties_test()", error); error_output += error; }

	// Test passing map of evaluator pointers to AQEoS constructor
	std::map<AQEoS::Model, AQBase*> evaluators = {
		{AQEoS::Model::Jager2003, &jag},
		{AQEoS::Model::Ziabakhsh2012, &zia}
	};
	aq = AQEoS(comp_data, evaluator_map, evaluators);
	error = aq.dlnphi_test(p, T, n, 5e-3);
	if (error || verbose) { print("AQEoS.dlnphi_test()", error); error_output += error; }
	error = aq.properties_test(p, T, n, 2e-3);
	if (error || verbose) { print("AQEoS.properties_test()", error); error_output += error; }

	// With ions
	// Na+ and Cl-
	std::vector<std::string> ions = {"Na+", "Cl-"};
	comp_data = CompData(comp, ions);
	comp_data.charge = {1, -1};

	n = {0.9, 0.045, 0.005, 0.025, 0.025};

	zia = Ziabakhsh2012(comp_data);
	error = zia.derivatives_test(p, T, n, 1e-5);
	if (error || verbose) { print("Ziabakhsh.derivatives_test() NaCl", error); error_output += error; }

	jag = Jager2003(comp_data);
	error = jag.derivatives_test(p, T, n, 2e-2);
	if (error || verbose) { print("Jager.derivatives_test() NaCl", error); error_output += error; }

	// Test Jager2003 with ions
	aq = AQEoS(comp_data, AQEoS::Model::Jager2003);
	error = aq.dlnphi_test(p, T, n, 2.e-2);
	if (error || verbose) { print("Jager.dlnphi_test() NaCl", error); error_output += error; }
	error = aq.properties_test(p, T, n, 2e-3);
	if (error || verbose) { print("Jager.properties_test() NaCl", error); error_output += error; }

	// Test Ziabakhsh2012 with ions
	aq = AQEoS(comp_data, AQEoS::Model::Ziabakhsh2012);
	error = aq.dlnphi_test(p, T, n, 1.e-3);
	if (error || verbose) { print("Ziabakhsh.dlnphi_test() NaCl", error); error_output += error; }
	error = aq.properties_test(p, T, n, 2e-3);
	if (error || verbose) { print("Ziabakhsh.properties_test() NaCl", error); error_output += error; }

	// Test mixed evaluators
	evaluators[AQEoS::Model::Jager2003] = &jag;
	evaluators[AQEoS::Model::Ziabakhsh2012] = &zia;

	aq = AQEoS(comp_data, evaluator_map, evaluators);
	error = aq.dlnphi_test(p, T, n, 3.e-3);
	if (error || verbose) { print("AQEoS.dlnphi_test() NaCl", error); error_output += error; }
	error = aq.properties_test(p, T, n, 2e-3);
	if (error || verbose) { print("AQEoS.properties_test() NaCl", error); error_output += error; }

	// Ca2+ and Cl-
	ions = {"Ca2+", "Cl-"};
	comp_data = CompData(comp, ions);
	comp_data.charge = {2, -1};

	n = {0.9, 0.045, 0.005, 0.01667, 0.03333};

	zia = Ziabakhsh2012(comp_data);
	error = zia.derivatives_test(p, T, n, 1e-5);
	if (error || verbose) { print("Ziabakhsh.derivatives_test() CaCl2", error); error_output += error; }

	jag = Jager2003(comp_data);
	// error = jag.derivatives_test(p, T, n, 2.1e-2);
	// if (error || verbose) { print("Jager.derivatives_test() CaCl2", error); error_output += error; }

	// Test Jager2003 with ions
	aq = AQEoS(comp_data, AQEoS::Model::Jager2003);
	error = aq.dlnphi_test(p, T, n, 1.e-1);
	if (error || verbose) { print("Jager.dlnphi_test() CaCl2", error); error_output += error; }
	error = aq.properties_test(p, T, n, 2e-3);
	if (error || verbose) { print("Jager.properties_test() CaCl2", error); error_output += error; }

	// Test Ziabakhsh2012 with ions
	aq = AQEoS(comp_data, AQEoS::Model::Ziabakhsh2012);
	error = aq.dlnphi_test(p, T, n, 1.e-3);
	if (error || verbose) { print("Ziabakhsh.dlnphi_test() CaCl2", error); error_output += error; }
	error = aq.properties_test(p, T, n, 2e-3);
	if (error || verbose) { print("Ziabakhsh.properties_test() CaCl2", error); error_output += error; }

	// Test mixed evaluators
	evaluators[AQEoS::Model::Jager2003] = &jag;
	evaluators[AQEoS::Model::Ziabakhsh2012] = &zia;

	aq = AQEoS(comp_data, evaluator_map, evaluators);
	error = aq.dlnphi_test(p, T, n, 3.e-3);
	if (error || verbose) { print("AQEoS.dlnphi_test() CaCl2", error); error_output += error; }
	error = aq.properties_test(p, T, n, 2e-3);
	if (error || verbose) { print("AQEoS.properties_test() CaCl2", error); error_output += error; }

	if (error_output > 0)
	{
		print("Errors occurred in test_aq()", error_output);
	}
	else
	{
		print("No errors occurred in test_aq()", error_output);
	}
    return error_output;
}

int test_iapws()
{
	// Test IAPWS-95 EoS
	int error, error_output = 0;
	bool verbose = false;
	double p{ 30. }, V{ 1.807725808e-05 }, T{ 300. };

	// Pure H2O
	std::vector<std::string> comp = {"H2O"};
	CompData comp_data(comp);
	std::vector<double> n = {2.};

	// use IAPWS95 ideal gas properties
	bool iapws_ideal = true;
	IAPWS95 iapws95(comp_data, iapws_ideal);

	// Test of analytical derivatives
	error = iapws95.dlnphi_test(p, T, n, 2.2e-2);
	if (error || verbose) { print("IAPWS95.dlnphi_test()", error); error_output += error; }

    // Consistency tests for IAPWS-95
	error = iapws95.derivatives_test(2*V, T, n, 5e-4);
	if (error || verbose) { print("IAPWS95.derivatives_test()", {static_cast<double>(error), 2*V, T}); error_output += error; }
	double V2 = 2. * iapws95::Mw * 1e-3 / (0.6585693359 * iapws95::rhoc);
	error = iapws95.derivatives_test(V2, 640., n, 5e-2);  // problematic conditions
	if (error || verbose) { print("IAPWS95.derivatives_test()", {static_cast<double>(error), V2, 640.}); error_output += error; }
	error = iapws95.lnphi_test(p, T, n, 2e-4);
	if (error || verbose) { print("IAPWS95.lnphi_test()", error); error_output += error; }
	error = iapws95.pressure_test(p, T, n, 1e-5);
	if (error || verbose) { print("IAPWS95.pressure_test()", error); error_output += error; }
	error = iapws95.temperature_test(p, T, n, 1e-5);
	if (error || verbose) { print("IAPWS95.temperature_test()", error); error_output += error; }
	error = iapws95.composition_test(p, T, n, 1e-5);
	if (error || verbose) { print("IAPWS95.composition_test()", error); error_output += error; }
	error = iapws95.pvt_test(p, T, n, 6e-3);
	if (error || verbose) { print("IAPWS95.pvt_test()", error); error_output += error; }
	error = iapws95.properties_test(p, T, n, 3e-3);
	if (error || verbose) { print("IAPWS95.properties_test()", error); error_output += error; }
	error = iapws95.references_test(1e-5);
	if (error || verbose) { print("IAPWS95.references_test()", error); error_output += error; }

	// Test PT solver
	std::vector<double> pressure = logspace(1e-5, 1e4, 10);
	std::vector<double> temperature = {270., 300., 330., 473., 630., 640., 643., 645., 647., 648., 700., 1000.};

	for (double temp: temperature)
	{
		for (double pres: pressure)
		{
			double v = iapws95.V(pres, temp, n);
			if (v < 1.e-2)
			{
				error = iapws95.pvt_test(pres, temp, n, 1e-2);
				if (error || verbose) { print("IAPWS95.pvt_test() IAPWS-95 ref", {static_cast<double>(error), pres, temp}); error_output += error; }
				error = iapws95.properties_test(pres, temp, n, 2e-3);
				if (error || verbose) { print("IAPWS95.properties_test() IAPWS-95 ref", {static_cast<double>(error), pres, temp}); error_output += error; }
			}
		}
	}

	// use default EoS ideal gas properties
	iapws_ideal = false;
	iapws95 = IAPWS95(comp_data, iapws_ideal);

	// Test PT solver
	for (double temp: temperature)
	{
		for (double pres: pressure)
		{
			double v = iapws95.V(pres, temp, n);
			if (v < 1.e-2)
			{
				error = iapws95.properties_test(pres, temp, n, 1e-3);
				if (error || verbose) { print("IAPWS95.properties_test() EoS", {static_cast<double>(error), pres, temp}); error_output += error; }
			}
		}
	}

	if (error_output > 0)
	{
		print("Errors occurred in test_iapws()", error_output);
	}
	else
	{
		print("No errors occurred in test_iapws()", error_output);
	}
	return error_output;
}

int test_solid()
{
    int error, error_output = 0;
	bool verbose = false;
	double p{ 30.}, T{ 273.15 };
	std::vector<double> n{ 2. };

	// Test IAPWSIce EoS
	CompData comp_data({"H2O"});
	IAPWSIce ice(comp_data, true);

	// error = ice.references_test(1e-3);
	// if (error || verbose) { print("IAPWSIce.references_test()", error); error_output += error; }
	error = ice.derivatives_test(p, T, n, 1e-5);
	if (error || verbose) { print("IAPWSIce.derivatives_test()", error); error_output += error; }
	error = ice.dlnphi_test(p, T, n, 1e-4);
	if (error || verbose) { print("IAPWSIce.dlnphi_test()", error); error_output += error; }
	error = ice.properties_test(p, T, n, 1e-3);
	if (error || verbose) { print("IAPWSIce.properties_test()", error); error_output += error; }
	error = ice.pvt_test(p, T, n, 1e-5);
	if (error || verbose) { print("IAPWSIce.pvt_test()", error); error_output += error; }

	// Test for PureSolid EoS for pure component only
	comp_data = CompData({"NaCl"});
	PureSolid s(comp_data, "NaCl");
	error = s.dlnphi_test(p, T, n, 2e-4);
	if (error || verbose) { print("PureSolid.dlnphi_test() NaCl", error); error_output += error; }
	error = s.pvt_test(p, T, n, 1e-5);
	if (error || verbose) { print("PureSolid.pvt_test() NaCl", error); error_output += error; }
	error = s.properties_test(p, T, n, 5e-3);
	if (error || verbose) { print("PureSolid.properties_test() NaCl", error); error_output += error; }
	error = s.derivatives_test(p, T, n, 1e-5);
	if (error || verbose) { print("PureSolid.derivatives_test() NaCl", error); error_output += error; }
	
	// Test PureSolid EoS for Ice for vector of components
	comp_data = CompData({"H2O", "CO2", "C1"});
	n = {1., 0., 0.};
	s = PureSolid(comp_data, "Ice");
	error = s.dlnphi_test(p, T, n, 2e-4);
	if (error || verbose) { print("PureSolid.dlnphi_test() Ice", error); error_output += error; }
	error = s.pvt_test(p, T, n, 1e-5);
	if (error || verbose) { print("PureSolid.pvt_test() Ice", error); error_output += error; }
	error = s.properties_test(p, T, n, 1e-5);
	if (error || verbose) { print("PureSolid.properties_test() Ice", error); error_output += error; }
	error = s.derivatives_test(p, T, n, 2e-5);
	if (error || verbose) { print("PureSolid.derivatives_test() Ice", error); error_output += error; }

	if (error_output > 0)
	{
		print("Errors occurred in test_solid()", error_output);
	}
	else
	{
		print("No errors occurred in test_solid()", error_output);
	}
	return error_output;
}

int test_vdwp()
{
	int error, error_output = 0;
	bool verbose = false;

	double p{ 30. }, T{ 280. };

	// Single-component sI hydrates
	std::vector<std::string> comp{"H2O", "CO2"};
	CompData comp_data(comp);
	std::vector<double> n{0.9, 0.1};

	Ballard ballard(comp_data, "sI");
	error = ballard.derivatives_test(p, T, n, 2e-4);
	if (error || verbose) { print("Ballard.derivatives_test() sI 1C", error); error_output += error; }
	error = ballard.dlnphi_test(p, T, n, 2e-3);
	if (error || verbose) { print("Ballard.dlnphi_test() sI 1C", error); error_output += error; }
	error = ballard.properties_test(p, T, n, 8e-4);
	if (error || verbose) { print("Ballard.properties_test() sI 1C", error); error_output += error; }

	Munck munck(comp_data, "sI");
	error = munck.derivatives_test(p, T, n, 1e-5);
	if (error || verbose) { print("Munck.derivatives_test() sI 1C", error); error_output += error; }
	error = munck.dlnphi_test(p, T, n, 1e-3);
	if (error || verbose) { print("Munck.dlnphi_test() sI 1C", error); error_output += error; }
	error = munck.properties_test(p, T, n, 2e-4);
	if (error || verbose) { print("Munck.properties_test() sI 1C", error); error_output += error; }

	// Multi-component sI hydrates
	comp = {"H2O", "CO2", "C1"};
	comp_data = CompData(comp);
	n = {0.86, 0.1, 0.04};
	
	ballard = Ballard(comp_data, "sI");
	error = ballard.derivatives_test(p, T, n, 2e-4);
	if (error || verbose) { print("Ballard.derivatives_test() sI 2C", error); error_output += error; }
	error = ballard.dlnphi_test(p, T, n, 2e-3);
	if (error || verbose) { print("Ballard.dlnphi_test() sI 2C", error); error_output += error; }
	error = ballard.properties_test(p, T, n, 2e-3);
	if (error || verbose) { print("Ballard.properties_test() sI 2C", error); error_output += error; }

	munck = Munck(comp_data, "sI");
	error = munck.derivatives_test(p, T, n, 2e-5);
	if (error || verbose) { print("Munck.derivatives_test() sI 2C", error); error_output += error; }
	error = munck.dlnphi_test(p, T, n, 1e-3);
	if (error || verbose) { print("Munck.dlnphi_test() sI 2C", error); error_output += error; }
	error = munck.properties_test(p, T, n, 7e-4);
	if (error || verbose) { print("Munck.properties_test() sI 2C", error); error_output += error; }

	// Single-component sII hydrates
	comp = {"H2O", "CO2"};
	comp_data = CompData(comp);
	n = {0.9, 0.1};

	ballard = Ballard(comp_data, "sII");
	error = ballard.derivatives_test(p, T, n, 2e-4);
	if (error || verbose) { print("Ballard.derivatives_test() sII 1C", error); error_output += error; }
	error = ballard.dlnphi_test(p, T, n, 2e-3);
	if (error || verbose) { print("Ballard.dlnphi_test() sII 1C", error); error_output += error; }
	error = ballard.properties_test(p, T, n, 3e-4);
	if (error || verbose) { print("Ballard.properties_test() sII 1C", error); error_output += error; }

	munck = Munck(comp_data, "sII");
	error = munck.derivatives_test(p, T, n, 1e-5);
	if (error || verbose) { print("Munck.derivatives_test() sII 1C", error); error_output += error; }
	error = munck.dlnphi_test(p, T, n, 1e-3);
	if (error || verbose) { print("Munck.dlnphi_test() sII 1C", error); error_output += error; }
	error = munck.properties_test(p, T, n, 2e-4);
	if (error || verbose) { print("Munck.properties_test() sII 1C", error); error_output += error; }

	// Multi-component sII hydrates
	comp = {"H2O", "CO2", "C1"};
	comp_data = CompData(comp);
	n = {0.86, 0.1, 0.04};
	
	ballard = Ballard(comp_data, "sII");
	error = ballard.derivatives_test(p, T, n, 2e-4);
	if (error || verbose) { print("Ballard.derivatives_test() sII 2C", error); error_output += error; }
	error = ballard.dlnphi_test(p, T, n, 3e-3);
	if (error || verbose) { print("Ballard.dlnphi_test() sII 2C", error); error_output += error; }
	error = ballard.properties_test(p, T, n, 2e-3);
	if (error || verbose) { print("Ballard.properties_test() sII 2C", error); error_output += error; }

	munck = Munck(comp_data, "sII");
	error = munck.derivatives_test(p, T, n, 2e-5);
	if (error || verbose) { print("Munck.derivatives_test() sII 2C", error); error_output += error; }
	error = munck.dlnphi_test(p, T, n, 1e-3);
	if (error || verbose) { print("Munck.dlnphi_test() sII 2C", error); error_output += error; }
	error = munck.properties_test(p, T, n, 3e-3);
	if (error || verbose) { print("Munck.properties_test() sII 2C", error); error_output += error; }

	if (error_output > 0)
	{
		print("Errors occurred in test_vdwp()", error_output);
	}
	else
	{
		print("No errors occurred in test_vdwp()", error_output);
	}
	return error_output;
}

int test_gmix()
{
	// Test implementation of Gibbs energy of mixing
	const bool verbose = false;
	int error_output = 0;

	// Binary mixture H2O-CO2
	std::vector<std::string> comp = {"H2O", "CO2"};
	CompData comp_data(comp);
	comp_data.Pc = {220.50, 73.75};
    comp_data.Tc = {647.14, 304.10};
    comp_data.ac = {0.328, 0.239};
    comp_data.kij = std::vector<double>(2*2, 0.);
	comp_data.set_binary_coefficients(0, {0., 0.19014});
	comp_data.cpi = {comp_data::cpi["H2O"], comp_data::cpi["CO2"]};

	IdealGas ig(comp_data);
	CubicEoS ceos(comp_data, CubicEoS::PR);
	ceos.set_preferred_roots(0, 0.6, EoS::RootFlag::MAX);

	std::map<AQEoS::CompType, AQEoS::Model> evaluator_map = {
		{AQEoS::CompType::water, AQEoS::Model::Jager2003},
		{AQEoS::CompType::solute, AQEoS::Model::Ziabakhsh2012},
		{AQEoS::CompType::ion, AQEoS::Model::Jager2003}
	};
	AQEoS aq(comp_data, evaluator_map);
	aq.set_eos_range(0, std::vector<double>{0.6, 1.});

	Ballard si(comp_data, "sI");
	Ballard sii(comp_data, "sII");

	PureSolid ice(comp_data, "Ice");
	IAPWSIce iapws_ice(comp_data, true);

	FlashParams flashparams(comp_data);
	// flashparams.add_eos("IG", &ig);
	flashparams.add_eos("CEOS", &ceos);
	flashparams.add_eos("AQ", &aq);
	flashparams.add_eos("sI", &si);
	flashparams.add_eos("sII", &sii);
	flashparams.add_eos("Ice", &ice);
	// flashparams.add_eos("Ice", &iapws_ice);
	flashparams.verbose = verbose;

	std::vector<std::pair<std::vector<double>, std::vector<double>>> references = {
		{{1., 273.15}, {-5.038285578, -0.007120529782}},
		{{10., 273.15}, {-7.333254784, -0.07266156188}},
		{{100., 273.15}, {-9.564637647, -1.201389886}},
		{{1., 373.15}, {-0.008472870891, -0.002637647353}},
		{{10., 373.15}, {-2.260979589, -0.02638413411}},
		{{100., 373.15}, {-4.509159992, -0.2619219487}},
		{{1., 473.15}, {-0.004420302368, -0.001043064989}},
		{{10., 473.15}, {-0.04489612603, -0.01035300458}},
		{{100., 473.15}, {-1.87757965, -0.09494685447}},
	};

	std::vector<std::vector<double>> n = {
		{0.95, 0.05}, {0.05, 0.95}
	};

	for (auto ref: references)
	{
		double p = ref.first[0];
		double T = ref.first[1];
		std::vector<double> gpure = flashparams.G_pure(p, T);

		for (size_t i = 0; i < gpure.size(); i++)
		{
			if (std::fabs(gpure[i]-ref.second[i]) > 1e-6)
			{
				print("Different values for Gpure", std::vector<double>{p, T, gpure[i]-ref.second[i]});
				print("result", gpure);
				print("ref", ref.second);
				error_output++;
			}
		}

		// Test for finding local minimum of Gibbs energy of mixing
		TrialPhase trial;

		trial = flashparams.find_ref_comp(p, T, n[0]);
		error_output += gmix_test(flashparams.eos_params[trial.eos_name].eos, p, T, n[0], gpure, verbose);

		trial = flashparams.find_ref_comp(p, T, n[1]);
		error_output += gmix_test(flashparams.eos_params[trial.eos_name].eos, p, T, n[1], gpure, verbose);
	}

	// Ternary mixture H2O-CO2-C1
	comp = {"H2O", "CO2", "C1"};
	comp_data = CompData(comp);
	comp_data.Pc = {220.50, 73.75, 46.04};
	comp_data.Tc = {647.14, 304.10, 190.58};
	comp_data.ac = {0.328, 0.239, 0.012};
	comp_data.kij = std::vector<double>(3*3, 0.);
    comp_data.set_binary_coefficients(0, {0., 0.19014, 0.47893});
	comp_data.set_binary_coefficients(1, {0.19014, 0., 0.0936});
	comp_data.Mw = {18.015, 44.01, 16.043};
	comp_data.cpi = {comp_data::cpi["H2O"], comp_data::cpi["CO2"], comp_data::cpi["C1"]};

	ig = IdealGas(comp_data);
	ceos = CubicEoS(comp_data, CubicEoS::PR);
	ceos.set_preferred_roots(0, 0.6, EoS::RootFlag::MAX);

	aq = AQEoS(comp_data, evaluator_map);
	aq.set_eos_range(0, std::vector<double>{0.6, 1.});

	si = Ballard(comp_data, "sI");
	sii = Ballard(comp_data, "sII");

	ice = PureSolid(comp_data, "Ice");

	flashparams = FlashParams(comp_data);
	// flashparams.add_eos("IG", &ig);
	flashparams.add_eos("CEOS", &ceos);
	flashparams.add_eos("AQ", &aq);
	flashparams.add_eos("sI", &si);
	flashparams.add_eos("sII", &sii);
	// flashparams.add_eos("Ice", &ice);
	flashparams.verbose = verbose;

	references = {
		{{1., 273.15}, {-5.037806881, -0.007120529782, -0.00293389646}},
		{{10., 273.15}, {-7.333254784, -0.07266156188, -0.02922065073}},
		{{100., 273.15}, {-9.564637647, -1.201389886, -0.2713925062}},
		{{1., 373.15}, {-0.008472870891, -0.002637647353, -0.0009766608968}},
		{{10., 373.15}, {-2.260979589, -0.02638413411, -0.009632032391}},
		{{100., 373.15}, {-4.509159992, -0.2619219487, -0.08209063185}},
		{{1., 473.15}, {-0.004420302368, -0.001043064989, -0.0002830630341}},
		{{10., 473.15}, {-0.04489612603, -0.01035300458, -0.002754486846}},
		{{100., 473.15}, {-1.87757965, -0.09494685447, -0.02013809187}},
	};

	n = {
		{0.95, 0.025, 0.025}, {0.025, 0.95, 0.025}, {0.025, 0.025, 0.95}
	};

	for (auto ref: references)
	{
		double p = ref.first[0];
		double T = ref.first[1];
		std::vector<double> gpure = flashparams.G_pure(p, T);

		for (size_t i = 0; i < gpure.size(); i++)
		{
			if (std::fabs(gpure[i]-ref.second[i]) > 1e-6)
			{
				print("Different values for Gpure", std::vector<double>{p, T, gpure[i]-ref.second[i]});
				print("result", gpure);
				print("ref", ref.second);
				error_output++;
			}
		}

		// Test for finding local minimum of Gibbs energy of mixing
		TrialPhase trial;

		trial = flashparams.find_ref_comp(p, T, n[0]);
		error_output += gmix_test(flashparams.eos_params[trial.eos_name].eos, p, T, n[0], gpure, verbose);

		trial = flashparams.find_ref_comp(p, T, n[1]);
		error_output += gmix_test(flashparams.eos_params[trial.eos_name].eos, p, T, n[1], gpure, verbose);

		trial = flashparams.find_ref_comp(p, T, n[2]);
		error_output += gmix_test(flashparams.eos_params[trial.eos_name].eos, p, T, n[2], gpure, verbose);
	}

	if (error_output > 0)
	{
		print("Errors occurred in test_gmix()", error_output);
	}
	else
	{
		print("No errors occurred in test_gmix()", error_output);
	}
	return error_output;
}