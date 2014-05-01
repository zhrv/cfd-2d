#ifndef _VISCOSITYMODEL_H_
#define _VISCOSITYMODEL_H_

#include "global.h"
#include "grid.h"

class ViscosityModel
{
public:
	virtual double getMuT(int iCell) = 0;

protected:
	Grid grid;
	
	double * ro;
	double * ru;
	double * rv;

	double * muT;
	double * rk;
	double * reps;

	Vector *gradK, *gradEps;
};

#endif