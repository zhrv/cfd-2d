#include "k_eps_model.h"


const double KEpsModel::C_mu = 0.09;
const double KEpsModel::C_eps1 = 1.44;
const double KEpsModel::C_eps2 = 1.92;
const double KEpsModel::Sigma_K = 1.0;
const double KEpsModel::Sigma_Eps = 1.3;

const double KEpsModel::It_Start = 0.01;
const double KEpsModel::Lt_Start = 1.4e-5;


KEpsModel::KEpsModel(void)
{
}


KEpsModel::~KEpsModel(void)
{
}

void KEpsModel::init( Grid * grid, double * ro, double *ru, double * rv, double * ro_m, double * u_m, double * v_m, Vector * gradU, Vector * gradV, double * Txx, double * Tyy, double * Txy, const double mu )
{
	this->grid = grid;
	this->ro = ro;
	this->ru = ru;
	this->rv = rv;

	this->ro_m = ro_m;
	this->u_m = u_m;
	this->v_m = v_m;

	this->gradU = gradU;
	this->gradV = gradV;

	this->Txx = Txx;
	this->Tyy = Tyy;
	this->Txy = Txy;

	this->mu = mu;

	this->muT = new double[grid->cCount];
	this->rk = new double[grid->cCount];
	this->reps = new double[grid->cCount];

	this->rk_int = new double[grid->cCount];
	this->reps_int = new double[grid->cCount];

	this->rk_int_kinematic_flow = new double[grid->cCount];
	this->rk_int_turbulent_diffusion = new double[grid->cCount];
	this->rk_int_generation = new double[grid->cCount];
	this->rk_int_dissipation = new double[grid->cCount];

	this->reps_int_kinematic_flow = new double[grid->cCount];
	this->reps_int_turbulent_diffusion = new double[grid->cCount];
	this->reps_int_generation = new double[grid->cCount];
	this->reps_int_dissipation = new double[grid->cCount];

	this->gradK = new Vector[grid->cCount];
	this->gradEps = new Vector[grid->cCount];

	/*
	memset(muT, 0, grid->cCount*sizeof(double));
	memset(rk, 0, grid->cCount*sizeof(double));
	memset(reps, 0, grid->cCount*sizeof(double));
	*/

	this->step = 0;
	startCond();
}


void KEpsModel::startCond()
{
	int nc = grid->cCount;
	
	for (int iCell = 0; iCell < nc; iCell++)
	{
		KEpsParam par;
		kEpsConvertConsToPar(iCell, par);
		
		double q = sqrt( par.u * par.u + par.v * par.v );
		
		par.k = 3.0 / 2.0 * ( It_Start * q ) * ( It_Start * q );
		par.eps = pow(C_mu, 3.0 / 4.0) * pow(par.k, 3.0 / 2.0) / Lt_Start;
		par.muT = par.r * C_mu * ( par.k * par.k / par.eps );

		kEpsConvertParToCons(iCell, par);
	}
}

inline double KEpsModel::getMuT( const int iCell )
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

void KEpsModel::calcMuT( double * cTau )
{
	step++;
	
	int nc = grid->cCount;
	int ne = grid->eCount;

	for (int iTau = 1; iTau <= 10; iTau++)
	{
		memset(rk_int, 0, nc*sizeof(double));
		memset(reps_int, 0, nc*sizeof(double));

		memset(rk_int_kinematic_flow, 0, nc*sizeof(double));
		memset(rk_int_turbulent_diffusion, 0, nc*sizeof(double));
		memset(rk_int_generation, 0, nc*sizeof(double));
		memset(rk_int_dissipation, 0, nc*sizeof(double));

		memset(reps_int_kinematic_flow, 0, nc*sizeof(double));
		memset(reps_int_turbulent_diffusion, 0, nc*sizeof(double));
		memset(reps_int_generation, 0, nc*sizeof(double));
		memset(reps_int_dissipation, 0, nc*sizeof(double));

		calcGrad();

		for (int iEdge = 0; iEdge < ne; iEdge++)
		{
			int c1	= grid->edges[iEdge].c1;
			int c2	= grid->edges[iEdge].c2;

			Vector n = grid->edges[iEdge].n;
			double l = grid->edges[iEdge].l;

			KEpsParam pL, pR;
			kEpsReconstruct(iEdge, pL, pR);

			double ro_m = (pL.r + pR.r) / 2.0;
			double u_m = (pL.u + pR.u) / 2.0;
			double v_m = (pL.v + pR.v) / 2.0;
			double k_m = (pL.k + pR.k) / 2.0;
			double eps_m = (pL.eps + pR.eps) / 2.0;
			double muT_m = (pL.muT + pR.muT) / 2.0;

			double Txx_m, Tyy_m, Txy_m;
			Vector gradK_m, gradEps_m;
			if (c2 <= -1)
			{
				Txx_m = Txx[c1];
				Tyy_m = Tyy[c1];
				Txy_m = Txy[c1];

				gradK_m.x = gradK[c1].x;
				gradK_m.y = gradK[c1].y;

				gradEps_m.x = gradEps[c1].x;
				gradEps_m.y = gradEps[c1].y;
			}
			else
			{
				Txx_m = (Txx[c1] + Txx[c2]) / 2.0;
				Tyy_m = (Tyy[c1] + Tyy[c2]) / 2.0;
				Txy_m = (Txy[c1] + Txy[c2]) / 2.0;

				gradK_m.x = (gradK[c1].x + gradK[c2].x) / 2.0;
				gradK_m.y = (gradK[c1].y + gradK[c2].y) / 2.0;

				gradEps_m.x = (gradEps[c1].x + gradEps[c2].x) / 2.0;
				gradEps_m.y = (gradEps[c1].y + gradEps[c2].y) / 2.0;
			}

			// Первая скобка
			register double tmp1 = ro_m * (u_m * n.x + v_m * n.y) * l;
			rk_int_kinematic_flow[c1] -= k_m * tmp1;
			if (c2 > -1)
				rk_int_kinematic_flow[c2] += k_m * tmp1;

			reps_int_kinematic_flow[c1] -= eps_m * tmp1;
			if (c2 > -1)
				reps_int_kinematic_flow[c2] += eps_m * tmp1;

			// Вторая скобка
			tmp1 = (mu + muT_m / Sigma_K) * ( gradK_m.x * n.x + gradK_m.y * n.y ) * l;
			rk_int_turbulent_diffusion[c1] += tmp1;
			if (c2 > -1)
				rk_int_turbulent_diffusion[c2] -= tmp1;

			tmp1 = (mu + muT_m / Sigma_Eps) * ( gradEps_m.x * n.x + gradEps_m.y * n.y ) * l;
			reps_int_turbulent_diffusion[c1] += tmp1;
			if (c2 > -1)
				reps_int_turbulent_diffusion[c2] -= tmp1;
		}

		for (int iCell = 0; iCell < nc; iCell++)
		{
			register double si = grid->cells[iCell].S;
			// TODO: проверить
			double rPk = Txx[iCell] * gradU[iCell].x + Txy[iCell] * ( gradU[iCell].y + gradV[iCell].x ) + Tyy[iCell] * gradV[iCell].y;

			// несжимаемые потоки!
			rPk += 2.0 / 3.0 * ( ( ( rk[iCell] + muT[iCell] * gradU->x ) * gradU->x ) + ( ( rk[iCell] + muT[iCell] * gradV->y ) * gradV->y ) );

			// TODO: знаки, знаки!
			rk_int_generation[iCell] += si * rPk;
			reps_int_generation[iCell] -= si * rPk * C_eps1 * reps[iCell] / rk[iCell];

			rk_int_dissipation[iCell] -= si * reps[iCell];
			reps_int_dissipation[iCell] -= si * C_eps2 * reps[iCell] * reps[iCell] / rk[iCell];

			rk_int[iCell] = rk_int_kinematic_flow[iCell] + rk_int_turbulent_diffusion[iCell] + rk_int_generation[iCell] + rk_int_dissipation[iCell];
			reps_int[iCell] = rk_int_kinematic_flow[iCell] + rk_int_turbulent_diffusion[iCell] + rk_int_generation[iCell] + rk_int_dissipation[iCell];
		}

		for (int iCell = 0; iCell < nc; iCell++)
		{
			// TODO: знаки, знаки!
			register double cfl = cTau[iCell] / 10.0 / grid->cells[iCell].S;
			rk[iCell] += cfl * rk_int[iCell];
			reps[iCell] += cfl * reps_int[iCell];

			muT[iCell] = C_mu * rk[iCell] * rk[iCell] / reps[iCell];
		}

		// saveTurbulentParamsToFile(step, iTau);
	}
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

		gradK[c1].x += ( pL.k+pR.k ) / 2.0 * n.x * l;
		gradK[c1].y += ( pL.k+pR.k ) / 2.0 * n.y * l;
		gradEps[c1].x += ( pL.eps+pR.eps ) / 2.0 * n.x * l;
		gradEps[c1].y += ( pL.eps+pR.eps ) / 2.0 * n.y * l;

		if (c2 > -1)
		{
			gradK[c2].x -= ( pL.k+pR.k ) / 2.0 * n.x * l;
			gradK[c2].y -= ( pL.k+pR.k ) / 2.0 * n.y * l;
			gradEps[c2].x -= ( pL.eps+pR.eps ) / 2.0 * n.x * l;
			gradEps[c2].y -= ( pL.eps+pR.eps ) / 2.0 * n.y * l;
		}

	}

	for (int iCell = 0; iCell < nc; iCell++)
	{
#ifdef _DEBUG
		KEpsParam par;
		kEpsConvertConsToPar(iCell, par);
#endif

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
	par.r = ro[iCell];
	par.u = ru[iCell]/ro[iCell];
	par.v = rv[iCell]/ro[iCell];
	
	par.k = rk[iCell]/ro[iCell];
	par.eps = reps[iCell]/ro[iCell];

	par.muT = muT[iCell];
}

void KEpsModel::kEpsConvertParToCons( int iCell, KEpsParam& par )
{
	ro[iCell] = par.r;
	ru[iCell] = par.u * par.r;
	rv[iCell] = par.v * par.r;

	rk[iCell] = par.k * par.r;
	reps[iCell] = par.eps * par.r;

	muT[iCell] = par.muT;
}

void KEpsModel::boundaryCond( int iEdge, KEpsParam& pL, KEpsParam& pR )
{
	// TODO: узнать, правильно ли это
	pR = pL;
}

void KEpsModel::fprintParams(FILE * file)
{
	fprintf(file, "SCALARS Turb_K float 1\nLOOKUP_TABLE default\n", grid->cCount);
	for (int i = 0; i < grid->cCount; i++)
	{
		KEpsParam p;
		kEpsConvertConsToPar(i, p);
		fprintf(file, "%25.16f ", p.k);
		if (i+1 % 8 == 0 || i+1 == grid->cCount) fprintf(file, "\n");
	}

	fprintf(file, "SCALARS Turb_Eps float 1\nLOOKUP_TABLE default\n", grid->cCount);
	for (int i = 0; i < grid->cCount; i++)
	{
		KEpsParam p;
		kEpsConvertConsToPar(i, p);
		fprintf(file, "%25.16f ", p.eps);
		if (i+1 % 8 == 0 || i+1 == grid->cCount) fprintf(file, "\n");
	}
}

void KEpsModel::saveTurbulentParamsToFile(int step, int iTau)
{
	/*
	char fName[50];

	sprintf(fName, "res_turb_%05d%05d.vtk", step, iTau);
	FILE * fp = fopen(fName, "w");
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "GASDIN data file\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp, "POINTS %d float\n", grid->nCount);
	for (int i = 0; i < grid->nCount; i++)
	{
		fprintf(fp, "%f %f %f  ", grid->nodes[i].x,  grid->nodes[i].y, 0.0);
		if (i+1 % 8 == 0) fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	fprintf(fp, "CELLS %d %d\n", grid->cCount, 4*grid->cCount);
	for (int i = 0; i < grid->cCount; i++)
	{
		fprintf(fp, "3 %d %d %d\n", grid->cells[i].nodesInd[0], grid->cells[i].nodesInd[1], grid->cells[i].nodesInd[2]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "CELL_TYPES %d\n", grid->cCount);
	for (int i = 0; i < grid->cCount; i++) fprintf(fp, "5\n");
	fprintf(fp, "\n");

	fprintf(fp, "SCALARS MuT float 1\nLOOKUP_TABLE default\n", grid->cCount);
	for (int i = 0; i < grid->cCount; i++)
	{
		fprintf(fp, "%25.16f ", viscosityModel->getMuT(i));
		if (i+1 % 8 == 0 || i+1 == grid->cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS MuT float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%25.16f ", viscosityModel->getMuT(i));
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS MuT float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%25.16f ", viscosityModel->getMuT(i));
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS MuT float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%25.16f ", viscosityModel->getMuT(i));
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS MuT float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%25.16f ", viscosityModel->getMuT(i));
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS MuT float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%25.16f ", viscosityModel->getMuT(i));
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS MuT float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%25.16f ", viscosityModel->getMuT(i));
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS MuT float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%25.16f ", viscosityModel->getMuT(i));
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS MuT float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%25.16f ", viscosityModel->getMuT(i));
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS MuT float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%25.16f ", viscosityModel->getMuT(i));
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS MuT float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%25.16f ", viscosityModel->getMuT(i));
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}


	fclose(fp);
	printf("File '%s' saved...\n", fName);
	*/
}