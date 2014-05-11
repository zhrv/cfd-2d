#ifndef _MatrixSolver_
#define _MatrixSolver_

#include "CSR.h"

class MatrixSolver
{
public:

	static const int RESULT_OK					= 0x0000;
	static const int RESULT_ERR_ZERO_DIAG = 0x0001;
	static const int RESULT_ERR_MAX_ITER = 0x0002;

	virtual ~MatrixSolver();

	void init(int cellsCount, int blockDimension);
	
	void zero();
	
	void setMatrElement(int i, int j, double** matrDim);
	void setRightElement(int i, double* vectDim);
	void addMatrElement(int i, int j, double** matrDim);
	void addRightElement(int i, double* vectDim);
	void createMatrElement(int i, int j);

	virtual int solve(double eps, int& maxIter) = 0;

	void printToFile(const char* fileName);

	CSRMatrix	*a;
	int			 blockDim;
	double		*b;
	double		*x;
};

class SolverZeidel: public MatrixSolver 
{
	virtual int solve(double eps, int& maxIter);
};

class SolverJacobi: public MatrixSolver
{
public:
	SolverJacobi() {
		tempXAlloc = false;
	}
	~SolverJacobi() {
		if (tempXAlloc) {
			delete [] tempX;
			tempXAlloc = false;
		}
	}
	virtual int	solve(double eps, int& maxIter);
	double			*tempX;
	bool			tempXAlloc;
};

#endif