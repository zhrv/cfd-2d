#include "sa_model.h"

const double SAModel::K = 0.4187;
const double SAModel::C_b1 = 0.1355;
const double SAModel::C_b2 = 0.622;
const double SAModel::C_v1 = 7.1;
const double SAModel::C_w1 = C_b1/(K*K) + (1 + C_b2)/Sigma;
const double SAModel::C_w2 = 0.3;
const double SAModel::C_w3 = 2.0;
const double SAModel::Sigma = 2.0/3.0;


SAModel::SAModel(void)
{
}


SAModel::~SAModel(void)
{
}

void SAModel::init( Grid * grid, double * ro, double *ru, double * rv, Vector * gradU, Vector * gradV, double * Txx, double * Tyy, double * Txy, const double mu )
{
	this->grid = grid;
	this->ro = ro;
	this->ru = ru;
	this->rv = rv;

	this->gradU = gradU;
	this->gradV = gradV;

	this->Txx = Txx;
	this->Tyy = Tyy;
	this->Txy = Txy;

	this->mu = mu;

	this->muT = new double[grid->cCount];
	this->rnul = new double[grid->cCount];

	this->rnul_int = new double[grid->cCount];

	this->gradNUL = new Vector[grid->cCount];

	memset(muT, 0, grid->cCount*sizeof(double));
	memset(rnul, 0, grid->cCount*sizeof(double));
}

inline double SAModel::getMuT( const int iCell )
{
	return muT[iCell];
}

void SAModel::done()
{
	delete [] muT;
	delete [] rnul;

	delete [] gradNUL;
}

void SAModel::calcMuT( const double TAU )
{
	int nc = grid->cCount;
	int ne = grid->eCount;

	double * min_dist = grid->min_dist;
	
	memset(rnul_int, 0, nc*sizeof(double));

	calcGrad();

	for (int iEdge = 0; iEdge < ne; iEdge++)
	{
		int c1	= grid->edges[iEdge].c1;
		int c2	= grid->edges[iEdge].c2;

		Vector n = grid->edges[iEdge].n;
		double l = grid->edges[iEdge].l;

		SAParam pL, pR;
		SAReconstruct(iEdge, pL, pR);

		if (c2 <= -1)
			continue;

		double ro_m = (pL.r + pR.r) / 2.0;
		double u_m = (pL.u + pR.u) / 2.0;
		double v_m = (pL.v + pR.v) / 2.0;
		double nul_m = (pL.nul + pR.nul) / 2.0;
		double nu_m = ro_m / mu;

		Vector gradNUL_m;
		gradNUL_m.x = (gradNUL[c1].x + gradNUL[c2].x) / 2.0;
		gradNUL_m.y = (gradNUL[c1].y + gradNUL[c2].y) / 2.0;

		// Первая скобка
		register double tmp1 = ro_m * (u_m * n.x + v_m * n.y) * l;
		rnul_int[c1] -= nul_m * tmp1;
		rnul_int[c2] += nul_m * tmp1;

		// Вторая скобка
		tmp1 = (ro_m * (nu_m + nul_m) / Sigma) * ( gradNUL_m.x * u_m + gradNUL_m.y * v_m ) * l;
		rnul_int[c1] += tmp1;
		rnul_int[c2] -= tmp1;
	}

	for (int iCell = 0; iCell < nc; iCell++)
	{
		register double si = grid->cells[iCell].S;
		// TODO: проверить
		double nu = mu / ro[iCell];

		double Xi = rnul[iCell] / nu;
		double f_v1 = pow(Xi, 3) / ( pow(Xi, 3) + pow(C_v1, 3) );
		double f_v2 = 1 - Xi / (1 + Xi * f_v1);
		double Omega = (gradU[iCell].y + gradV[iCell].x) / 2.0; //TODO: переделать
		double Sl = fabs(Omega) + ( rnul[iCell] / pow(K * min_dist[iCell], 2) ) * f_v2;

		double r = rnul[iCell] / ( pow(K * min_dist[iCell], 2) * Sl );
		double g = r + C_w2 * (pow(r, 6) - r);
		double f_w = g * pow( (1 + pow(C_w3, 6))/(pow(g, 6) + pow(C_w3, 6)), 1/6 );
		double eps_v = C_w1 * f_w * pow(rnul[iCell] / min_dist[iCell], 2);
		double rPnul = C_b1 * Sl * rnul[iCell];
		// TODO: знаки, знаки!
		
		rnul_int[iCell] += si * C_b2 * pow(gradNUL[iCell].x + gradNUL[iCell].y, 2) / Sigma + si * rPnul - si * eps_v; //TODO: переделать
	}

	
	for (int iCell = 0; iCell < nc; iCell++)
	{
		// TODO: знаки, знаки!
		register double cfl = TAU/grid->cells[iCell].S;
		rnul[iCell] += cfl * rnul_int[iCell];

		double nu = mu / ro[iCell];
		double Xi = rnul[iCell] / nu;
		double f_v1 = pow(Xi, 3) / ( pow(Xi, 3) + pow(C_v1, 3) );
		
		muT[iCell] = rnul[iCell] * f_v1;
	}
	
}

void SAModel::calcGrad()
{
	int nc = grid->cCount;
	int ne = grid->eCount;

	memset(gradNUL, 0, nc * sizeof(Vector));

	for (int iEdge = 0; iEdge < ne; iEdge++)
	{

		int c1 = grid->edges[iEdge].c1;
		int c2 = grid->edges[iEdge].c2;

        SAParam pL, pR;
		// TODO: здесь нужен только nul
		SAReconstruct(iEdge, pL, pR);

		Vector n = grid->edges[iEdge].n;
		double l = grid->edges[iEdge].l;

		gradNUL[c1].x += (pL.nul+pR.nul)/2.0*n.x*l;
		gradNUL[c1].y += (pL.nul+pR.nul)/2.0*n.y*l;

		if (c2 > -1)
		{
			gradNUL[c2].x -= (pL.nul+pR.nul)/2.0*n.x*l;
			gradNUL[c2].y -= (pL.nul+pR.nul)/2.0*n.y*l;
		}

	}

	for (int iCell = 0; iCell < nc; iCell++)
	{
		register double si = grid->cells[iCell].S;
		gradNUL[iCell].x /= si;
		gradNUL[iCell].y /= si;
	}
}

void SAModel::SAReconstruct( int iEdge, SAParam& pL, SAParam& pR )
{
	if (grid->edges[iEdge].type == Edge::TYPE_INNER) 
	{
		int c1	= grid->edges[iEdge].c1;
		int c2	= grid->edges[iEdge].c2;
		SAConvertConsToPar(c1, pL);
		SAConvertConsToPar(c2, pR);
	} else {
		int c1	= grid->edges[iEdge].c1;
		SAConvertConsToPar(c1, pL);
		boundaryCond(iEdge, pL, pR);
	}
}

void SAModel::SAConvertConsToPar( int iCell, SAParam& par )
{
	par.r = ro[iCell];
	par.u = ru[iCell]/ro[iCell];
	par.v = rv[iCell]/ro[iCell];
	
	par.nul = rnul[iCell]/ro[iCell];

	par.muT = muT[iCell];
}

void SAModel::boundaryCond( int iEdge, SAParam& pL, SAParam& pR )
{
	// TODO: узнать, правильно ли это
	pR = pL;
}

