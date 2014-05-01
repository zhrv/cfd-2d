#ifndef _KEPSMODEL_H_
#define _KEPSMODEL_H_

#include "viscosity_model.h"

class KEpsModel :
	public ViscosityModel
{
public:
	KEpsModel(void);
	~KEpsModel(void);

	void init(Grid * grid, double * ro, double *ru, double * rv);
	double getMuT(int iCell);
	void calcMuT( const double TAU );
	void done();
private:
	double * rk;
	double * reps;

	Vector *gradK, *gradEps;
};

#endif