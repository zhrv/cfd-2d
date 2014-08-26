#pragma once

#include "LimiterDG.h"
#include "grid.h"
#include "FEM_DG.h"

class LimiterDGCockburn : public LimiterDG
{
public:
	LimiterDGCockburn(FEM_DG *mthd, Grid *grid, double **ro, double **ru, double **rv, double **re, int fCount);
	virtual ~LimiterDGCockburn();
	virtual void run();
	virtual const char* getName() { return "Cockburn"; }
private:
	void initLimiterParameters();
	int __getEdgeByCells(int c1, int c2);
	void choiseDirection(int& nn1, int& nn2, double& a1, double& a2, int n0, int n1, int n2, int n3, Point pm, int m);
	void getMatrLR(double** L, double** R, double u, double v, double c);
	double triSquare(Point p0, Point p1, Point p2);
	void calcLimiter_II();
private:
	Grid *grid;
	int cellsCount;

	double*** limAlfa;
	int*** limNeigh;
	Vector** limLm;
	Vector** limLmN;
	Point** limPm;

	double* deltaS1;
	double* deltaS2;
	double** deltaU1;
	double** deltaU2;

	double** matrL;
	double** matrR;

	FEM_DG *method;

	int funcCount;

	double **fRO;
	double **fRU;
	double **fRV;
	double **fRE;

	double **fROlim;
	double **fRUlim;
	double **fRVlim;
	double **fRElim;

	const double LIMITER_ALFA = 1.9;

	const double GAM = 1.4;
	const double AGAM = GAM - 1.0;

	const double EPS = 1.0e-10;
};

