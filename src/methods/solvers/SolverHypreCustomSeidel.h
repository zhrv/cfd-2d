#pragma once
#include "SolverHypre.h"
class SolverHypreCustomSeidel :
	public SolverHypre
{
public:
	SolverHypreCustomSeidel();
	virtual ~SolverHypreCustomSeidel();

	virtual void init(int cellsCount, int blockDimension);
	virtual void zero();

	virtual void setX(double* x);

	virtual int solve(double eps, int& maxIter);
	virtual char* getName() { return (char*)"Seidel based on the HYPRE structures"; }

	virtual void addRightElement(int i, double* vectDim);
	virtual void setRightElement(int i, double* vectDim);
};

