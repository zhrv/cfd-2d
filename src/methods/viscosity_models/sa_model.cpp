#include "sa_model.h"

const double SAModel::K = 0.4187;
const double SAModel::C_b1 = 0.1355;
const double SAModel::C_b2 = 0.622;
const double SAModel::C_v1 = 7.1;
const double SAModel::C_w1 = C_b1/(K*K) + (1 + C_b2)/Sigma;
const double SAModel::C_w2 = 0.3;
const double SAModel::C_w3 = 2.0;
const double SAModel::Sigma = 2.0/3.0;

const double SAModel::It_Start = 0.01;
const double SAModel::Lt_Start = 1.4e-5;


SAModel::SAModel(void)
{
} 


SAModel::~SAModel(void)
{
}

void SAModel::init( Grid * grid, double * ro, double *ru, double * rv, double * ro_m, double * u_m, double * v_m, Vector * gradU, Vector * gradV, double * Txx, double * Tyy, double * Txy, const double mu )
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

	this->ro_m = ro_m;
	this->u_m = u_m;
	this->v_m = v_m;

	this->muT = new double[grid->cCount]; 
	this->rnt = new double[grid->cCount];

	this->rnt_int = new double[grid->cCount];

	this->gradNT = new Vector[grid->cCount];

	/*
	memset(muT, 0, grid->cCount*sizeof(double));
	memset(rnt, 0, grid->cCount*sizeof(double));
	*/
	startCond();
}

void SAModel::startCond()
{
	int nc = grid->cCount;

	for (int iCell = 0; iCell < nc; iCell++)
	{
		SAParam par;
		SAConvertConsToPar(iCell, par);

		double q = sqrt( par.u * par.u + par.v * par.v );
		par.nt = sqrt( 3.0 / 2.0 ) * It_Start * Lt_Start * q;
		double hi = par.r * par.nt / mu;
		double f_v1 = pow(hi, 3.0) / ( pow(hi, 3.0) + pow(C_v1, 3.0));
		par.muT = par.r * par.nt * f_v1;

		SAConvertParToCons(iCell, par);
	}
}

inline double SAModel::getMuT( const int iCell )
{
	return muT[iCell];
}

void SAModel::done()
{
	delete [] muT;
	delete [] rnt;

	delete [] gradNT;
}

void SAModel::calcMuT( const double TAU )
{
	int nc = grid->cCount;
	int ne = grid->eCount;

	double * min_dist = grid->min_dist;
	
	memset(rnt_int, 0, nc*sizeof(double));

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

		ro_m[iEdge] = (pL.r + pR.r) / 2.0;
		u_m[iEdge] = (pL.u + pR.u) / 2.0;
		v_m[iEdge] = (pL.v + pR.v) / 2.0;
		double nt_m = (pL.nt + pR.nt) / 2.0;
		double nu_m = mu / ro_m[iEdge];

		Vector gradNT_m;
		gradNT_m.x = (gradNT[c1].x + gradNT[c2].x) / 2.0;
		gradNT_m.y = (gradNT[c1].y + gradNT[c2].y) / 2.0;

		register double tmp1 = ro_m[iEdge] * (u_m[iEdge] * n.x + v_m[iEdge] * n.y) * l;
		rnt_int[c1] -= nt_m * tmp1;
		rnt_int[c2] += nt_m * tmp1;

		tmp1 = ro_m[iEdge] * ((nu_m + nt_m) / Sigma) * ( gradNT_m.x * n.x + gradNT_m.y * n.y ) * l;
		rnt_int[c1] += tmp1;
		rnt_int[c2] -= tmp1;
	}

	for (int iCell = 0; iCell < nc; iCell++)
	{
		register double si = grid->cells[iCell].S;

		double nu = mu / ro[iCell];

		double hi = rnt[iCell] / nu;
		double f_v1 = pow(hi, 3.0) / ( pow(hi, 3.0) + pow(C_v1, 3.0) );
		double f_v2 = 1 - hi / (1 + hi * f_v1);
		double omega = (gradU[iCell].y - gradV[iCell].x) / 2.0;
		double S = sqrt(2.0) * fabs(omega) + ( rnt[iCell] / pow(K * min_dist[iCell], 2.0) ) * f_v2;

		double r = rnt[iCell] / ( pow(K * min_dist[iCell], 2.0) * S );
		double g = r + C_w2 * (pow(r, 6.0) - r);
		double f_w = g * pow( (1 + pow(C_w3, 6.0))/(pow(g, 6.0) + pow(C_w3, 6.0)), 1.0/6.0 );
		double eps = C_w1 * f_w * pow(rnt[iCell] / min_dist[iCell], 2.0);
		double rP = C_b1 * S * rnt[iCell];
		
		rnt_int[iCell] += si * ( C_b2 * pow(gradNT[iCell].x + gradNT[iCell].y, 2.0) / Sigma + rP - eps );
	}

	
	for (int iCell = 0; iCell < nc; iCell++)
	{
		register double cfl = TAU/grid->cells[iCell].S;
		rnt[iCell] += cfl * rnt_int[iCell];

		double nu = mu / ro[iCell];
		double hi = rnt[iCell] / nu;
		double f_v1 = pow(hi, 3.0) / ( pow(hi, 3.0) + pow(C_v1, 3.0) );
		
		muT[iCell] = rnt[iCell] * f_v1;
	}
	
}

void SAModel::calcGrad()
{
	int nc = grid->cCount;
	int ne = grid->eCount;

	memset(gradNT, 0, nc * sizeof(Vector));

	for (int iEdge = 0; iEdge < ne; iEdge++)
	{
		int c1 = grid->edges[iEdge].c1;
		int c2 = grid->edges[iEdge].c2;

        SAParam pL, pR;
		SAReconstruct(iEdge, pL, pR);

		Vector n = grid->edges[iEdge].n;
		double l = grid->edges[iEdge].l;

		gradNT[c1].x += (pL.nt+pR.nt)/2.0*n.x*l;
		gradNT[c1].y += (pL.nt+pR.nt)/2.0*n.y*l;

		if (c2 > -1)
		{
			gradNT[c2].x -= (pL.nt+pR.nt)/2.0*n.x*l;
			gradNT[c2].y -= (pL.nt+pR.nt)/2.0*n.y*l;
		}
	}

	for (int iCell = 0; iCell < nc; iCell++)
	{
		register double si = grid->cells[iCell].S;
		gradNT[iCell].x /= si;
		gradNT[iCell].y /= si;
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
	
	par.nt = rnt[iCell]/ro[iCell];

	par.muT = muT[iCell];
}

void SAModel::SAConvertParToCons( int iCell, SAParam& par )
{
	ro[iCell] = par.r;
	ru[iCell] = par.u * par.r;
	rv[iCell] = par.v * par.r;

	rnt[iCell] = par.nt * par.r;

	muT[iCell] = par.muT;
}

void SAModel::boundaryCond( int iEdge, SAParam& pL, SAParam& pR )
{
	pR = pL;
}

