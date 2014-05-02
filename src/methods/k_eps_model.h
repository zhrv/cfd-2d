#ifndef _KEPSMODEL_H_
#define _KEPSMODEL_H_

#include "viscosity_model.h"

/**
 *	Турбулентные параметры
 */
struct KEpsParam
{
	double k;		//!< плотность
	double eps;		//!< давление
	
};

class KEpsModel :
	public ViscosityModel
{
public:
	KEpsModel(void);
	~KEpsModel(void);

	void init(Grid * grid, double * ro, double *ru, double * rv, double * Txx, double * Tyy, double * Txy);
	double getMuT(int iCell);
	void calcMuT( const double TAU );
	void done();

private:
	double * rk;
	double * reps;

	double * rk_int;
	double * reps_int;

	Vector *gradK, *gradEps;

	void calcGrad();
	void kEpsReconstruct( int iEdge, KEpsParam& pL, KEpsParam& pR );
	void kEpsConvertConsToPar( int iCell, KEpsParam&  par );
	void boundaryCond( int iEdge, KEpsParam& pL, KEpsParam& pR );
};

#endif