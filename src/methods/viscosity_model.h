#ifndef _VISCOSITYMODEL_H_
#define _VISCOSITYMODEL_H_

#include "global.h"
#include "grid.h"

class ViscosityModel
{
public:
	virtual void init(Grid * grid, double * ro, double *ru, double * rv, double * Txx, double * Tyy, double * Txy) = 0;
	virtual double getMuT(int iCell) = 0;
	virtual void calcMuT( const double TAU ) = 0;
	virtual void done() = 0;

protected:
	Grid * grid;
	
	double * ro;
	double * ru;
	double * rv;
	double *Txx, *Tyy, *Txy;

	double * muT;
};

#endif