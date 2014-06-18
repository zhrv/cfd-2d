#ifndef _EMPTYVISCOSITYMODEL_H_
#define _EMPTYVISCOSITYMODEL_H_

#include "viscosity_model.h"

class EmptyViscosityModel : public ViscosityModel
{
public:
	EmptyViscosityModel(void);
	~EmptyViscosityModel(void);

	void init( Grid * grid, double * ro, double *ru, double * rv, double * ro_m, double * u_m, double * v_m, Vector * gradU, Vector * gradV, double * Txx, double * Tyy, double * Txy, const double mu );
	double getMuT(const int iCell);
	void calcMuT( double * cTau );
	void done();

	void fprintParams(FILE * file);
};

#endif