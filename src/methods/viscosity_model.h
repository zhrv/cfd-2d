#ifndef _VISCOSITYMODEL_H_
#define _VISCOSITYMODEL_H_

#include "global.h"
#include "grid.h"

class ViscosityModel
{
public:
	virtual void init(Grid * grid, double * ro, double *ru, double * rv, double * ro_m, double * u_m, double * v_m, Vector * gradU, Vector * gradV, double * Txx, double * Tyy, double * Txy, const double mu) = 0;
	virtual double getMuT(const int iCell) = 0;
	virtual void calcMuT( double * cTau ) = 0;
	virtual void done() = 0;

	virtual void fprintParams(FILE * file) = 0;

protected:
	Grid * grid;
	
	double * ro;
	double * ru;
	double * rv;

	double * ro_m;
	double * u_m;
	double * v_m;

	Vector *gradU, *gradV;
	double *Txx, *Tyy, *Txy;

	double mu;
	double * muT;
};

#endif