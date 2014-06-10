#ifndef _SAMODEL_H_
#define _SAMODEL_H_

#include "viscosity_model.h"

/**
 *	Турбулентные параметры
 */
struct SAParam
{
	double r;		//!< плотность
	double u;		//!< первая компонента вектора скорости
	double v;		//!< вторая компонента вектора скорости

	double nt;		//!< ню с волной
	double muT;		//!< тубулентная вязкость
};

class SAModel :
	public ViscosityModel
{
public:
	SAModel(void);
	~SAModel(void);

	void init(Grid * grid, double * ro, double *ru, double * rv, double * ro_m, double * u_m, double * v_m, Vector * gradU, Vector * gradV, double * Txx, double * Tyy, double * Txy, const double mu);
	double getMuT(const int iCell);
	void calcMuT( const double TAU );
	void done();

private:
	static const double K;
	static const double C_b1;
	static const double C_b2;
	static const double C_v1;
	static const double C_w1;
	static const double C_w2;
	static const double C_w3;
	static const double Sigma;
	
	static const double It_Start;
	static const double Lt_Start;

	double * rnt;

	double * rnt_int;

	Vector *gradNT;

	void calcGrad();
	void startCond();
	void SAReconstruct( int iEdge, SAParam& pL, SAParam& pR );
	void SAConvertConsToPar( int iCell, SAParam& par );
	void SAConvertParToCons( int iCell, SAParam& par );
	void boundaryCond( int iEdge, SAParam& pL, SAParam& pR );
};

#endif