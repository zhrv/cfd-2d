#include "k_eps_model.h"


KEpsModel::KEpsModel(void)
{
}


KEpsModel::~KEpsModel(void)
{
}

void KEpsModel::init( Grid * grid, double * ro, double *ru, double * rv, double * Txx, double * Tyy, double * Txy )
{
	this->grid = grid;
	this->ro = ro;
	this->ru = ru;
	this->rv = rv;

	this->Txx = Txx;
	this->Tyy = Tyy;
	this->Txy = Txy;

	this->muT = new double[grid->cCount];
	this->rk = new double[grid->cCount];
	this->reps = new double[grid->cCount];

	this->rk_int = new double[grid->cCount];
	this->reps_int = new double[grid->cCount];

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
	int nc = grid->cCount;
	int ne = grid->eCount;
	
	memset(rk_int, 0, nc*sizeof(double));
	memset(reps_int, 0, nc*sizeof(double));

	calcGrad();

	// to be continued...
}

void KEpsModel::calcGrad()
{
	int nc = grid->cCount;
	int ne = grid->eCount;

	memset(gradK, 0, nc * sizeof(Vector));
	memset(gradEps, 0, nc * sizeof(Vector));

	for (int iEdge = 0; iEdge < ne; iEdge++)
	{

		int c1 = grid->edges[iEdge].c1;
		int c2 = grid->edges[iEdge].c2;

        KEpsParam pL, pR;
		kEpsReconstruct(iEdge, pL, pR);

		Vector n = grid->edges[iEdge].n;
		double l = grid->edges[iEdge].l;

		gradK[c1].x += (pL.k+pR.k)/2*n.x*l;
		gradK[c1].y += (pL.k+pR.k)/2*n.y*l;
		gradEps[c1].x += (pL.eps+pR.eps)/2*n.x*l;
		gradEps[c1].y += (pL.eps+pR.eps)/2*n.y*l;

		if (c2 > -1)
		{
			gradK[c2].x -= (pL.k+pR.k)/2*n.x*l;
			gradK[c2].y -= (pL.k+pR.k)/2*n.y*l;
			gradEps[c2].x -= (pL.eps+pR.eps)/2*n.x*l;
			gradEps[c2].y -= (pL.eps+pR.eps)/2*n.y*l;
		}

	}

	for (int iCell = 0; iCell < nc; iCell++)
	{
		register double si = grid->cells[iCell].S;
		gradK[iCell].x /= si;
		gradK[iCell].y /= si;
		gradEps[iCell].x /= si;
		gradEps[iCell].y /= si;
	}
}

void KEpsModel::kEpsReconstruct( int iEdge, KEpsParam& pL, KEpsParam& pR )
{
	if (grid->edges[iEdge].type == Edge::TYPE_INNER) 
	{
		int c1	= grid->edges[iEdge].c1;
		int c2	= grid->edges[iEdge].c2;
		kEpsConvertConsToPar(c1, pL);
		kEpsConvertConsToPar(c2, pR);
	} else {
		int c1	= grid->edges[iEdge].c1;
		kEpsConvertConsToPar(c1, pL);
		boundaryCond(iEdge, pL, pR);
	}
}

void KEpsModel::kEpsConvertConsToPar( int iCell, KEpsParam& par )
{
	par.k = rk[iCell]/ro[iCell];
	par.eps = reps[iCell]/ro[iCell];
}

void KEpsModel::boundaryCond( int iEdge, KEpsParam& pL, KEpsParam& pR )
{
	// TODO: узнать, правильно ли это
	pR = pL;
}

