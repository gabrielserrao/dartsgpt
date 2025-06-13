//--------------------------------------------------------------------------
#ifndef OPENDARTS_FLASH_FLASH_PXFLASH_H
#define OPENDARTS_FLASH_FLASH_PXFLASH_H
//--------------------------------------------------------------------------

#include <vector>
#include "dartsflash/global/global.hpp"
#include "dartsflash/flash/flash_params.hpp"
#include "dartsflash/flash/flash_results.hpp"
#include "dartsflash/flash/flash.hpp"
#include "dartsflash/flash/trial_phase.hpp"

class PXFlash : public Flash
{
protected:
	int total_ssi_flash_iter, total_ssi_stability_iter, total_newton_flash_iter, total_newton_stability_iter;
	int error = 0;
	double XX, X_spec;  // current and specified state specification
	double T_a, T_b;  // a and b values for Brent's method root finding

public:	
	PXFlash(FlashParams& flashparams, Flash::StateSpecification state_spec_);
	~PXFlash() = default;

	virtual int evaluate(double p_, double X_spec_) override;
	virtual int evaluate(double p_, double X_spec_, std::vector<double>& z_) override;

	int evaluate_PT(double p_, double T_) { return Flash::evaluate(p_, T_); }
	int evaluate_PT(double p_, double T_, std::vector<double>& z_) { return Flash::evaluate(p_, T_, z_); }

protected:
	virtual void init(double p_, double X_spec_) override;
	virtual void init(double p_, double X_spec_, std::vector<double>& z_) override;
	double obj_fun(double T_);
	double gradient(double T_);

public:
	std::vector<double>& getnu() { return this->nu; }
	std::vector<double>& getx() { return this->X; }
	double getT() { return this->T; }
	double getSpec(double T_test) { return this->obj_fun(T_test); }

};

//--------------------------------------------------------------------------
#endif // OPENDARTS_FLASH_FLASH_PXFLASH_H
//--------------------------------------------------------------------------
