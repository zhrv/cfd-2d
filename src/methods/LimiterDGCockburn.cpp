#include "LimiterDGCockburn.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

LimiterDGCockburn::LimiterDGCockburn(FEM_DG *mthd, Grid *grid, double **ro, double **ru, double **rv, double **re, int fCount)
{
	this->method = mthd;
	this->grid = grid;
	cellsCount = grid->cCount;
	funcCount = fCount;
	initLimiterParameters();

	fRO = ro;
	fRU = ru;
	fRV = rv;
	fRE = re;
}


LimiterDGCockburn::~LimiterDGCockburn()
{
	// TODO: очистка памяти
}



void LimiterDGCockburn::initLimiterParameters()
{
	limAlfa = new double**[cellsCount];
	limNeigh = new int**[cellsCount];
	limLm = new Vector*[cellsCount];
	limLmN = new Vector*[cellsCount];
	limPm = new Point*[cellsCount];

	fROlim = new double*[cellsCount];
	fRUlim = new double*[cellsCount];
	fRVlim = new double*[cellsCount];
	fRElim = new double*[cellsCount];

	for (int i = 0; i < cellsCount; i++)
	{
		limAlfa[i] = new double*[3];
		limAlfa[i][0] = new double[2];
		limAlfa[i][1] = new double[2];
		limAlfa[i][2] = new double[2];

		limNeigh[i] = new int*[3];
		limNeigh[i][0] = new int[2];
		limNeigh[i][1] = new int[2];
		limNeigh[i][2] = new int[2];

		limLm[i] = new Vector[3];
		limLmN[i] = new Vector[3];
		limPm[i] = new Point[3];

		fROlim[i] = new double[funcCount];
		fRUlim[i] = new double[funcCount];
		fRVlim[i] = new double[funcCount];
		fRElim[i] = new double[funcCount];

	}

	deltaS1 = new double[4];
	deltaS2 = new double[4];
	deltaU1 = new double*[3];
	deltaU2 = new double*[3];
	for (int m = 0; m < 3; m++)
	{
		deltaU1[m] = new double[4];
		deltaU2[m] = new double[4];
	}

	matrL = new double*[4];
	matrR = new double*[4];
	for (int i = 0; i < 4; i++)
	{
		matrL[i] = new double[4];
		matrR[i] = new double[4];
	}
	
	// находим угол и коэффициенты разложения
	for (int iCell = 0; iCell < cellsCount; iCell++)
	{
		int n0 = grid->cells[iCell].neigh[0];
		int n1 = grid->cells[iCell].neigh[1];
		int n2 = grid->cells[iCell].neigh[2];
		for (int m = 0; m < 3; m++)
		{
			int iEdge = __getEdgeByCells(iCell, grid->cells[iCell].neigh[m]);
			limPm[iCell][m].x = grid->edges[iEdge].c->x; 
			limPm[iCell][m].y = grid->edges[iEdge].c->y;
			limLm[iCell][m].x = grid->edges[iEdge].c->x - grid->cells[iCell].c.x;
			limLm[iCell][m].y = grid->edges[iEdge].c->y - grid->cells[iCell].c.y;
			double tmp = sqrt(limLm[iCell][m].x*limLm[iCell][m].x + limLm[iCell][m].y*limLm[iCell][m].y);
			limLmN[iCell][m].x = limLm[iCell][m].x / tmp;
			limLmN[iCell][m].y = limLm[iCell][m].y / tmp;
			choiseDirection(limNeigh[iCell][m][0], limNeigh[iCell][m][1], limAlfa[iCell][m][0], limAlfa[iCell][m][1], iCell, n0, n1, n2, limPm[iCell][m], m);
			if (limAlfa[iCell][m][0]<0.0 || limAlfa[iCell][m][1]<0.0) log("ERROR!!!\n");
		}
		int zhrv = 0;
	}
}

inline double MINMOD(double a, double b)
{
	//if (abs(a) < EPS  ||  abs(b) < EPS) return 0.0;
	if (a*b>0.0)
	{
		if (abs(a)<abs(b))
			return a;
		else
			return b;
	}
	else {
		return 0;
	}
}

inline double MINMOD_B(double a, double b)
{
	//if (abs(a) < MM_M*MM_DX*MM_DX) return a;
	return MINMOD(a, b);
}


void LimiterDGCockburn::getMatrLR(double** L, double** R, double u, double v, double c)
{
	double h = c*c / AGAM + 0.5*(u*u + v*v);

	R[0][0] = 1.0;		R[0][1] = 1.0;		R[0][2] = 0.0;	R[0][3] = 1.0;
	R[1][0] = u - c;		R[1][1] = u;		R[1][2] = 0.0;	R[1][3] = u + c;
	R[2][0] = v;		R[2][1] = v;		R[2][2] = 1.0;	R[2][3] = v;
	R[3][0] = h - u*c;	R[3][1] = 0.5*v*v;	R[3][2] = v;	R[3][3] = h + u*c;

	//R[0][0] = 1.0;		R[0][1] = 0.0;		R[0][2] = 0.0;	R[0][3] = 0.0;
	//R[1][0] = 0.0;		R[1][1] = 1.0;		R[1][2] = 0.0;	R[1][3] = 0.0;
	//R[2][0] = 0.0;		R[2][1] = 0.0;		R[2][2] = 1.0;	R[2][3] = 0.0;
	//R[3][0] = 0.0;		R[3][1] = 0.0;		R[3][2] = 0.0;	R[3][3] = 1.0;

	inverseMatr(R, L, 4);
}

void LimiterDGCockburn::run()
{
	//static const double GAM = 1.4;
	//static const double AGAM = GAM - 1.0;
	for (int i = 0; i < cellsCount; i++)
	{
		memcpy(fROlim[i], fRO[i], funcCount*sizeof(double));
		memcpy(fRUlim[i], fRU[i], funcCount*sizeof(double));
		memcpy(fRVlim[i], fRV[i], funcCount*sizeof(double));
		memcpy(fRElim[i], fRE[i], funcCount*sizeof(double));

		//memcpy(fROold[i], fRO[i], funcCount*sizeof(double));
		//memcpy(fRUold[i], fRU[i], funcCount*sizeof(double));
		//memcpy(fRVold[i], fRV[i], funcCount*sizeof(double));
		//memcpy(fREold[i], fRE[i], funcCount*sizeof(double));
	}



	//return;
	for (int iCell = 0; iCell < cellsCount; iCell++) // цикл по ячейкам
	{
		int n0 = grid->cells[iCell].neigh[0];
		int n1 = grid->cells[iCell].neigh[1];
		int n2 = grid->cells[iCell].neigh[2];
		if ((n0 < 0) || (n1 < 0) || (n2 < 0))
		{
			fROlim[iCell][1] = 0.0;
			fRUlim[iCell][1] = 0.0;
			fRVlim[iCell][1] = 0.0;
			fRElim[iCell][1] = 0.0;

			fROlim[iCell][2] = 0.0;
			fRUlim[iCell][2] = 0.0;
			fRVlim[iCell][2] = 0.0;
			fRElim[iCell][3] = 0.0;

			continue;
		}

		double ROc, RUc, RVc, REc;
		method->getFields(ROc, RUc, RVc, REc, iCell, grid->cells[iCell].c.x, grid->cells[iCell].c.y);
		for (int m = 0; m < 3; m++) // цикл по направлениям
		{
			// вычисляем приращения
			double ROm, RUm, RVm, REm;
			double RO1, RU1, RV1, RE1;
			double RO2, RU2, RV2, RE2;

			method->getFields(ROm, RUm, RVm, REm, iCell, limPm[iCell][m].x, limPm[iCell][m].y);
			method->getFields(RO1, RU1, RV1, RE1, limNeigh[iCell][m][0], grid->cells[limNeigh[iCell][m][0]].c.x, grid->cells[limNeigh[iCell][m][0]].c.y);
			method->getFields(RO2, RU2, RV2, RE2, limNeigh[iCell][m][1], grid->cells[limNeigh[iCell][m][1]].c.x, grid->cells[limNeigh[iCell][m][1]].c.y);

			deltaU1[m][0] = ROm - ROc;
			deltaU1[m][1] = RUm - RUc;
			deltaU1[m][2] = RVm - RVc;
			deltaU1[m][3] = REm - REc;

			deltaU2[m][0] = limAlfa[iCell][m][0] * (RO1 - ROc) + limAlfa[iCell][m][1] * (RO2 - ROc);
			deltaU2[m][1] = limAlfa[iCell][m][0] * (RU1 - RUc) + limAlfa[iCell][m][1] * (RU2 - RUc);
			deltaU2[m][2] = limAlfa[iCell][m][0] * (RV1 - RVc) + limAlfa[iCell][m][1] * (RV2 - RVc);
			deltaU2[m][3] = limAlfa[iCell][m][0] * (RE1 - REc) + limAlfa[iCell][m][1] * (RE2 - REc);

			// скорости по направлению
			double un, ut, c, utmp, vtmp;
			double nx = limLmN[iCell][m].x;
			double ny = limLmN[iCell][m].y;
			utmp = RUc / ROc;
			vtmp = RVc / ROc;
			un = utmp*nx + vtmp*ny;
			ut = utmp*ny - vtmp*nx;
			c = sqrt(GAM*AGAM*(REc / ROc - (utmp*utmp + vtmp*vtmp) / 2.0));

			// поворачиваем скорости
			utmp = deltaU1[m][1] * nx + deltaU1[m][2] * ny;
			vtmp = deltaU1[m][1] * ny - deltaU1[m][2] * nx;
			deltaU1[m][1] = utmp;
			deltaU1[m][2] = vtmp;

			utmp = deltaU2[m][1] * nx + deltaU2[m][2] * ny;
			vtmp = deltaU2[m][1] * ny - deltaU2[m][2] * nx;
			deltaU2[m][1] = utmp;
			deltaU2[m][2] = vtmp;

			// переходим к инвариантам
			getMatrLR(matrL, matrR, un, ut, c);
			memset(deltaS1, 0, 4 * sizeof(double));
			memset(deltaS2, 0, 4 * sizeof(double));
			for (int k = 0; k < 4; k++)
			{
				deltaS1[0] += matrL[0][k] * deltaU1[m][k];
				deltaS1[1] += matrL[1][k] * deltaU1[m][k];
				deltaS1[2] += matrL[2][k] * deltaU1[m][k];
				deltaS1[3] += matrL[3][k] * deltaU1[m][k];

				deltaS2[0] += matrL[0][k] * deltaU2[m][k];
				deltaS2[1] += matrL[1][k] * deltaU2[m][k];
				deltaS2[2] += matrL[2][k] * deltaU2[m][k];
				deltaS2[3] += matrL[3][k] * deltaU2[m][k];
			}

			// лимитируем инварианты
			deltaS1[0] = MINMOD_B(deltaS1[0], LIMITER_ALFA*deltaS2[0]);
			deltaS1[1] = MINMOD_B(deltaS1[1], LIMITER_ALFA*deltaS2[1]);
			deltaS1[2] = MINMOD_B(deltaS1[2], LIMITER_ALFA*deltaS2[2]);
			deltaS1[3] = MINMOD_B(deltaS1[3], LIMITER_ALFA*deltaS2[3]);


			// переходим к консервативным
			memset(deltaU1[m], 0, 4 * sizeof(double));
			for (int k = 0; k < 4; k++)
			{
				deltaU1[m][0] += matrR[0][k] * deltaS1[k];
				deltaU1[m][1] += matrR[1][k] * deltaS1[k];
				deltaU1[m][2] += matrR[2][k] * deltaS1[k];
				deltaU1[m][3] += matrR[3][k] * deltaS1[k];
			}

			// поворачиваем скорости обратно
			utmp = deltaU1[m][1] * nx + deltaU1[m][2] * ny;
			vtmp = deltaU1[m][1] * ny - deltaU1[m][2] * nx;
			deltaU1[m][1] = utmp;
			deltaU1[m][2] = vtmp;
		}

		double pos, neg;
		for (int k = 0; k < 4; k++) // цикл по компонентам {RO,RU,RV,RE}
		{

			double & d1n = deltaU1[0][k];
			double & d2n = deltaU1[1][k];
			double & d3n = deltaU1[2][k];
			// делаем поправки для приращений
			if (abs(deltaU1[0][k] + deltaU1[1][k] + deltaU1[2][k]) > EPS)
			{
				pos = 0; neg = 0;
				for (int m = 0; m < 3; m++)
				{
					pos += MAX(0.0, deltaU1[m][k]);
					neg += MAX(0.0, -deltaU1[m][k]);
				}
				double thetaP, thetaM;
				if (abs(pos)<EPS)
				{
					thetaP = 1.0;
					thetaM = 0.0;
				}
				else if (abs(neg)<EPS)
				{
					thetaP = 0.0;
					thetaM = 1.0;
				}
				else
				{
					thetaP = MIN(1.0, neg / pos);
					thetaM = MIN(1.0, pos / neg);
				}
				deltaU1[0][k] = thetaP * MAX(0.0, deltaU1[0][k]) - thetaM * MAX(0.0, -deltaU1[0][k]);
				deltaU1[1][k] = thetaP * MAX(0.0, deltaU1[1][k]) - thetaM * MAX(0.0, -deltaU1[1][k]);
				deltaU1[2][k] = thetaP * MAX(0.0, deltaU1[2][k]) - thetaM * MAX(0.0, -deltaU1[2][k]);
			}
			// вычисляем отлимитированные коэффициенты разложения
			double& xb0 = grid->cells[iCell].c.x;
			double& yb0 = grid->cells[iCell].c.y;
			double& xm1 = limPm[iCell][0].x;
			double& ym1 = limPm[iCell][0].y;
			double& xm2 = limPm[iCell][1].x;
			double& ym2 = limPm[iCell][1].y;
			double& xm3 = limPm[iCell][2].x;
			double& ym3 = limPm[iCell][2].y;

			double det = (xm2*ym3 - xm3*ym2) - (xm1*ym3 - xm3*ym1) + (xm1*ym2 - xm2*ym1);

			double a1 = (ym2 - ym3) / det;
			double b1 = (xm3 - xm2) / det;
			//double c1 = (xm2*ym3-xm3*ym2)/det;

			double a2 = (ym3 - ym1) / det;
			double b2 = (xm1 - xm3) / det;
			//double c2 = (xm3*ym1-xm1*ym3)/det;

			double a3 = (ym1 - ym2) / det;
			double b3 = (xm2 - xm1) / det;
			//double c3 = (xm1*ym2-xm2*ym1)/det;



			double *pRR, *pRC;
			switch (k)
			{
			case 0:
				pRR = fROlim[iCell];
				pRC = &ROc;
				break;
			case 1:
				pRR = fRUlim[iCell];
				pRC = &RUc;
				break;
			case 2:
				pRR = fRVlim[iCell];
				pRC = &RVc;
				break;
			case 3:
				pRR = fRElim[iCell];
				pRC = &REc;
				break;
			}

			pRR[1] = (a1*d1n + a2*d2n + a3*d3n)*grid->cells[iCell].HX;
			pRR[2] = (b1*d1n + b2*d2n + b3*d3n)*grid->cells[iCell].HY;
			//pRR[0] = (*pRC);

			//if (abs(pRR[0]) < EPS) pRR[0] = 0.0;
			if (abs(pRR[1]) < EPS) pRR[1] = 0.0;
			if (abs(pRR[2]) < EPS) pRR[2] = 0.0;
		}

	}

	for (int i = 0; i < cellsCount; i++)
	{
		memcpy(fRO[i], fROlim[i], funcCount*sizeof(double));
		memcpy(fRU[i], fRUlim[i], funcCount*sizeof(double));
		memcpy(fRV[i], fRVlim[i], funcCount*sizeof(double));
		memcpy(fRE[i], fRElim[i], funcCount*sizeof(double));
	}


	calcLimiter_II();
}

/*!
Вычисление вектора площади
*/
double LimiterDGCockburn::triSquare(Point p0, Point p1, Point p2)
{
	double s = (p1.x - p0.x)*(p2.y - p0.y) - (p2.x - p0.x)*(p1.y - p0.y);
	return (abs(s) < EPS) ? 0.0 : (p1.x - p0.x)*(p2.y - p0.y) - (p2.x - p0.x)*(p1.y - p0.y);
}




void LimiterDGCockburn::choiseDirection(int& nn1, int& nn2, double& a1, double& a2, int n0, int n1, int n2, int n3, Point pm, int mm)
{
	int nn[3] = { n1, n2, n3 };
	Point p[3];
	p[0] = grid->cells[n1].c;
	p[1] = grid->cells[n2].c;
	p[2] = grid->cells[n3].c;
	for (int m = 0; m < 3; m++)
	{
		nn1 = nn[m];
		nn2 = nn[(m + 1) % 3];
		if (nn1 >= 0 && nn2 >= 0)
		{
			double D  = (grid->cells[nn1].c.x - grid->cells[n0].c.x)*(grid->cells[nn2].c.y - grid->cells[n0].c.y) - (grid->cells[nn1].c.y - grid->cells[n0].c.y)*(grid->cells[nn2].c.x - grid->cells[n0].c.x);
			double D1 = (pm.x                 - grid->cells[n0].c.x)*(grid->cells[nn2].c.y - grid->cells[n0].c.y) - (pm.y                 - grid->cells[n0].c.y)*(grid->cells[nn2].c.x - grid->cells[n0].c.x);
			double D2 = (grid->cells[nn1].c.x - grid->cells[n0].c.x)*(pm.y                 - grid->cells[n0].c.y) - (grid->cells[nn1].c.y - grid->cells[n0].c.y)*(pm.x                 - grid->cells[n0].c.x);

			a1 = D1 / D;
			a2 = D2 / D;
			//_NORM_(a1, a2);
			if (a1 >= 0.0 && a2 >= 0.0) return;
		}
	}
	if (n1 < 0 || n2 < 0 || n3 < 0)
	{
		a1 = 0.0;
		a2 = 0.0;
		nn1 = 0;
		nn2 = 0;
		return;
	}
	else {
		log(" ERROR: choiseDirection()\n");
		log("N0 = %d, N1 = %d, N2 = %d, N3 = %d, pm.x = %f, pm.y = %f, m = %d\n", n0, n1, n2, n3, pm.x, pm.y, mm);
		log("p0.x = %f, p0.y = %f\n", grid->cells[n0].c.x, grid->cells[n0].c.y);
		log("p1.x = %f, p1.y = %f\n", grid->cells[n1].c.x, grid->cells[n1].c.y);
		log("p2.x = %f, p2.y = %f\n", grid->cells[n2].c.x, grid->cells[n2].c.y);
		log("p3.x = %f, p3.y = %f\n", grid->cells[n3].c.x, grid->cells[n3].c.y);
		log("\n\n");
		for (int m = 0; m < 3; m++)
		{
			nn1 = nn[m];
			nn2 = nn[(m + 1) % 3];
			log("nn1 = %d, nn2 = %d\n", nn1, nn2);
			if (nn1 >= 0 && nn2 >= 0)
			{
				double D = (grid->cells[nn1].c.x - grid->cells[n0].c.x)*(grid->cells[nn2].c.y - grid->cells[n0].c.y) - (grid->cells[nn1].c.y - grid->cells[n0].c.y)*(grid->cells[nn2].c.x - grid->cells[n0].c.x);
				double D1 = (pm.x - grid->cells[n0].c.x)*(grid->cells[nn2].c.y - grid->cells[n0].c.y) - (pm.y - grid->cells[n0].c.y)*(grid->cells[nn2].c.x - grid->cells[n0].c.x);
				double D2 = (grid->cells[nn1].c.x - grid->cells[n0].c.x)*(pm.y - grid->cells[n0].c.y) - (grid->cells[nn1].c.y - grid->cells[n0].c.y)*(pm.x - grid->cells[n0].c.x);

				log("D = %f, D1 = %f, D2 = %f\n", D, D1, D2);

				a1 = D1 / D; if (abs(a1)<EPS) a1 = 0.0;
				a2 = D2 / D; if (abs(a2)<EPS) a2 = 0.0;
				log("a1 = %f, a2 = %f\n", a1, a2);
				//if (a1>=0.0 && a2>=0.0) return;
			}
		}
		exit(1);
	}
}

int LimiterDGCockburn::__getEdgeByCells(int c1, int c2)
{
	for (int iEdge = 0; iEdge < grid->eCount; iEdge++)
	{
		Edge &edge = grid->edges[iEdge];
		if ((edge.c1 == c1 && edge.c2 == c2) || (edge.c1 == c2 && edge.c2 == c1)) return iEdge;
	}
}


void LimiterDGCockburn::calcLimiter_II()
{
	double xtt[9], ytt[9], t[3];
	double Rcur1, RUcur1, RVcur1, Ecur1;
	double Rcur2, RUcur2, RVcur2, Ecur2;
	double Rcur3, RUcur3, RVcur3, Ecur3;
	double sR, sRU, sRV, sE;
	double Pmin, Amin, Amax, Acur, fcur, fmax, ffcur, ffmax;
	int pointP;
	double gamma = GAM;



	for (int k = 0; k < cellsCount; k++)
	{
		double x1 = grid->getNode(grid->cells[k].nodesInd[0]).x; //nodes[cellNodes[k][0]].x;//xi1(int(point1(k)))
		double x2 = grid->getNode(grid->cells[k].nodesInd[1]).x;//xi1(int(point2(k)))
		double x3 = grid->getNode(grid->cells[k].nodesInd[2]).x;//xi1(int(point3(k)))

		double y1 = grid->getNode(grid->cells[k].nodesInd[0]).y;//yi1(int(point1(k)))
		double y2 = grid->getNode(grid->cells[k].nodesInd[1]).y;//yi1(int(point2(k)))
		double y3 = grid->getNode(grid->cells[k].nodesInd[2]).y;//yi1(int(point3(k)))

		double xc = (x1 + x2 + x3) / 3.;
		double yc = (y1 + y2 + y3) / 3.;

		//dx = fdx(k)
		//dy = fdy(k)

		method->getFields(Rcur1, RUcur1, RVcur1, Ecur1, k, x1, y1);
		method->getFields(Rcur2, RUcur2, RVcur2, Ecur2, k, x2, y2);
		method->getFields(Rcur3, RUcur3, RVcur3, Ecur3, k, x3, y3);

		if ((Ecur1 - (RUcur1*RUcur1 + RVcur1*RVcur1) / (2.*Rcur1) > 0.) && (Ecur2 - (RUcur2*RUcur2 + RVcur2*RVcur2) / (2.*Rcur2) > 0.) && (Ecur3 - (RUcur3*RUcur3 + RVcur3*RVcur3) / (2.*Rcur3) > 0.)){
			goto lbl_Lim_1;

		}
		else{
			t[0] = -sqrt(3.0 / 5.0);   //!-0.7745966
			t[1] = 0.;
			t[2] = -t[0];

			double feps = 1.e-13;

			for (int i = 0; i < 3; i++) {
				xtt[i] = (x2 - x1)*(t[i] + 1) / 2. + x1;
				ytt[i] = (y2 - y1)*(t[i] + 1) / 2. + y1;
			}
			for (int i = 3; i < 6; i++) {
				xtt[i] = (x3 - x2)*(t[i - 3] + 1) / 2. + x2;
				ytt[i] = (y3 - y2)*(t[i - 3] + 1) / 2. + y2;
			}

			for (int i = 6; i < 9; i++) {

				xtt[i] = (x1 - x3)*(t[i - 6] + 1) / 2. + x3;
				ytt[i] = (y1 - y3)*(t[i - 6] + 1) / 2. + y3;
			}

			double fRcur, Rmin;
			for (int j = 0; j < 9; j++)
			{
				fRcur = fRO[k][0] + fRO[k][1] * method->getF(1, k, xtt[j], ytt[j]) + fRO[k][2] * method->getF(2, k, xtt[j], ytt[j]);
				if (j == 0)
				{
					Rmin = fRcur;
				}
				else
				{
					if (fRcur <= Rmin) Rmin = fRcur;
				}
			}


			double alfR = (fRO[k][0] == Rmin) ? 1.0 : MIN((fRO[k][0] - feps) / (fRO[k][0] - Rmin), 1.);
			for (int i = 1; i < funcCount; i++) fRO[k][i] *= alfR;


			for (int j = 0; j < 9; j++)
			{
				double sR, sRU, sRV, sE;
				method->getFields(sR, sRU, sRV, sE, k, xtt[j], ytt[j]);

				double sep = sE / sR - ((sRU*sRU) / (sR*sR) + (sRV*sRV) / (sR*sR)) / 2.;

				double P = (gamma - 1.)*sR*sep;







				if (j == 0){
					Pmin = P;
					pointP = j;
				}
				else{
					if (P <= Pmin){
						Pmin = P;
						pointP = j;
					}
				}
			}


			if (Pmin >= feps){
				goto lbl_Lim_1;
			}
			else {
				Amin = 0.;
				Amax = 1.;
				Acur = (Amax - Amin) / 2. + Amin;

				sR = 0.;
				sRU = 0.;
				sRV = 0.;
				sE = 0.;


				for (int i = 1; i < funcCount; i++)
				{
					sR = sR + fRO[k][i] * method->getF(i, k, xtt[pointP], ytt[pointP]);//R(k,i)*fiI(i,xtt(pointP),xc,ytt(pointP),yc,dx,dy)
					sRU = sRU + fRU[k][i] * method->getF(i, k, xtt[pointP], ytt[pointP]);//RU(k,i)*fiI(i,xtt(pointP),xc,ytt(pointP),yc,dx,dy)
					sRV = sRV + fRV[k][i] * method->getF(i, k, xtt[pointP], ytt[pointP]);//RV(k,i)*fiI(i,xtt(pointP),xc,ytt(pointP),yc,dx,dy)
					sE = sE + fRE[k][i] * method->getF(i, k, xtt[pointP], ytt[pointP]);//E(k,i)*fiI(i,xtt(pointP),xc,ytt(pointP),yc,dx,dy)
				}

				fcur = fRE[k][0] + Acur*sE - ((fRU[k][0] + Acur*sRU)*(fRU[k][0] + Acur*sRU) + (fRV[k][0] + Acur*sRV)*(fRV[k][0] + Acur*sRV)) / (2.*(fRO[k][0] + Acur*sR));
				fmax = fRE[k][0] + Amax*sE - ((fRU[k][0] + Amax*sRU)*(fRU[k][0] + Amax*sRU) + (fRV[k][0] + Amax*sRV)*(fRV[k][0] + Amax*sRV)) / (2.*(fRO[k][0] + Amax*sR));
				ffmax = (gamma - 1.)*fmax - feps;
				ffcur = (gamma - 1.)*fcur - feps;

				double tfl = 0.;

				while (abs(ffmax)>feps) {
					if (tfl>10) goto lbl_Lim_2;
					if (ffmax*ffcur>0.)
						Amax = Acur;
					else
						Amin = Acur;
					Acur = (Amax - Amin) / 2. + Amin;
					for (int i = 1; i < funcCount; i++)
					{
						sR = sR + fRO[k][i] * method->getF(i, k, xtt[pointP], ytt[pointP]);//R(k,i)*fiI(i,xtt(pointP),xc,ytt(pointP),yc,dx,dy)
						sRU = sRU + fRU[k][i] * method->getF(i, k, xtt[pointP], ytt[pointP]);//RU(k,i)*fiI(i,xtt(pointP),xc,ytt(pointP),yc,dx,dy)
						sRV = sRV + fRV[k][i] * method->getF(i, k, xtt[pointP], ytt[pointP]);//RV(k,i)*fiI(i,xtt(pointP),xc,ytt(pointP),yc,dx,dy)
						sE = sE + fRE[k][i] * method->getF(i, k, xtt[pointP], ytt[pointP]);//E(k,i)*fiI(i,xtt(pointP),xc,ytt(pointP),yc,dx,dy)
					}

					fcur = fRE[k][0] + Acur*sE - ((fRU[k][0] + Acur*sRU)*(fRU[k][0] + Acur*sRU) + (fRV[k][0] + Acur*sRV)*(fRV[k][0] + Acur*sRV)) / (2.*(fRO[k][0] + Acur*sR));
					fmax = fRE[k][0] + Amax*sE - ((fRU[k][0] + Amax*sRU)*(fRU[k][0] + Amax*sRU) + (fRV[k][0] + Amax*sRV)*(fRV[k][0] + Amax*sRV)) / (2.*(fRO[k][0] + Amax*sR));
					ffmax = (gamma - 1.)*fmax - feps;
					ffcur = (gamma - 1.)*fcur - feps;


					tfl = tfl + 1;
				}


			lbl_Lim_2:
				for (int i = 1; i< funcCount; i++)
				{
					fRO[k][i] *= Acur;
					fRU[k][i] *= Acur;
					fRV[k][i] *= Acur;
					fRE[k][i] *= Acur;
				}
			}
		}
	lbl_Lim_1:;
	}

}


