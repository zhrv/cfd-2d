#ifndef _KEPSMODEL_H_
#define _KEPSMODEL_H_

#include "viscosity_model.h"

/**
 *	Турбулентные параметры
 */
struct KEpsParam
{
	double r;		//!< плотность
	double u;		//!< первая компонента вектора скорости
	double v;		//!< вторая компонента вектора скорости

	double k;
	double eps;
	double muT;	
};

class KEpsModel :
	public ViscosityModel
{
public:
	KEpsModel(void);
	~KEpsModel(void);

	void init( Grid * grid, double * ro, double *ru, double * rv, double * ro_m, double * u_m, double * v_m, Vector * gradU, Vector * gradV, double * Txx, double * Tyy, double * Txy, const double mu );
	double getMuT(const int iCell);
	void calcMuT( const double TAU );
	void done();

private:
	static const double C_mu;
	static const double C_eps1;
	static const double C_eps2;
	static const double Sigma_K;
	static const double Sigma_Eps;

	static const double It_Start;
	static const double Lt_Start;
	
	double * rk;
	double * reps;

	double * rk_int;
	double * reps_int;

	Vector *gradK, *gradEps;

	void calcGrad();
	void kEpsReconstruct( int iEdge, KEpsParam& pL, KEpsParam& pR );
	void kEpsConvertConsToPar( int iCell, KEpsParam& par );
	void kEpsConvertParToCons( int iCell, KEpsParam& par );
	void startCond();
	void boundaryCond( int iEdge, KEpsParam& pL, KEpsParam& pR );
};

#endif