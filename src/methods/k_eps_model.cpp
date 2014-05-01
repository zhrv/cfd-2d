#include "k_eps_model.h"


KEpsModel::KEpsModel(void)
{
}


KEpsModel::~KEpsModel(void)
{
}

void KEpsModel::init( Grid * grid, double * ro, double *ru, double * rv )
{
	this->grid = grid;
	this->ro = ro;
	this->ru = ru;
	this->rv = rv;

	this->muT = new double[grid->cCount];
	this->rk = new double[grid->cCount];
	this->reps = new double[grid->cCount];

	this->gradK = new Vector[grid->cCount];
	this->gradEps = new Vector[grid->cCount];

	memset(muT, 0, grid->cCount*sizeof(double));
	memset(rk, 0, grid->cCount*sizeof(double));
	memset(reps, 0, grid->cCount*sizeof(double));
}

inline double KEpsModel::getMuT( int iCell )
{
	return muT[iCell];
}

void KEpsModel::done()
{
	delete [] muT;
	delete [] rk;
	delete [] reps;

	delete [] gradK;
	delete [] gradEps;
}

void KEpsModel::calcMuT( const double TAU )
{
	
}

