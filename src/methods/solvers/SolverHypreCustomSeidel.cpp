#include "SolverHypreCustomSeidel.h"
#include "global.h"
#include "_hypre_utilities.h"

SolverHypreCustomSeidel::SolverHypreCustomSeidel()
{
}


SolverHypreCustomSeidel::~SolverHypreCustomSeidel()
{
	HYPRE_IJMatrixDestroy(A);
	delete[] x;
	delete[] b;
	delete[] cols;
	delete[] values;
}


void SolverHypreCustomSeidel::setX(double* x)
{
	MatrixSolver::setX(x);
}

int SolverHypreCustomSeidel::solve(double eps, int& maxIter)
{
	int result = MatrixSolver::RESULT_OK;

	/* Assemble after setting the coefficients */
	HYPRE_IJMatrixAssemble(A);

	/* Get the parcsr matrix object to use */
	HYPRE_IJMatrixGetObject(A, (void**)&parcsr_A);


	double	aii;
	double	err = 1.0;
	int		step = 0;
	double	tmp;
	double	*values;
	int		*cols;
	int		size = 0;


	while (err > eps && step < maxIter)
	{
		step++;
		for (int i = ilower; i <= iupper; i++)
		{
			HYPRE_ParCSRMatrixGetRow(parcsr_A, i, (HYPRE_Int*)&size, (HYPRE_Int**)&cols, &values);
			HYPRE_ParCSRMatrixRestoreRow(parcsr_A, i, (HYPRE_Int*)&size, (HYPRE_Int**)&cols, &values);
			tmp = 0.0;
			aii = 0;
			for (int k = 0; k < size; k++) {
				if (cols[k] == i) {
					aii = values[k];
				}
				else {
					tmp += values[k] * x[cols[k]];
				}
			}
//			if (fabs(aii) <= eps*eps)
//			{
//				log((char*)"ZEIDEL_SOLVER: error: a[%d, %d] = 0\n", i, i);
//				return MatrixSolver::RESULT_ERR_ZERO_DIAG;
//			}
			x[i] = (-tmp + b[i]) / aii;
		}
		err = 0.0;
		for (int i = ilower; i <= iupper; i++)
		{
			HYPRE_ParCSRMatrixGetRow(parcsr_A, i, (HYPRE_Int*)&size, (HYPRE_Int**)&cols, &values);
			HYPRE_ParCSRMatrixRestoreRow(parcsr_A, i, (HYPRE_Int*)&size, (HYPRE_Int**)&cols, &values);
			tmp = 0.0;
			for (int k = 0; k < size; k++) {
				tmp += values[k] * x[cols[k]];
			}
			err += fabs(tmp - b[i]);
		}
		if (PRINT_LEVEL > 0) {
			printf("SEIDEL SOLVER: step = %5d\terr = %16.8e\n", step, err);
		}
	}
	if (step >= maxIter)
	{
		log((char*)"ZEIDEL_SOLVER: (warning) maximum iterations done (%d); error: %e\n", step, err);
		maxIter = step;
		return MatrixSolver::RESULT_ERR_MAX_ITER;
	}
	maxIter = step;
	return MatrixSolver::RESULT_OK;
}



void SolverHypreCustomSeidel::addRightElement(int i, double* vectDim)
{
	MatrixSolver::addRightElement(i, vectDim);
}

void SolverHypreCustomSeidel::setRightElement(int i, double* vectDim)
{
	MatrixSolver::setRightElement(i, vectDim);
}


void SolverHypreCustomSeidel::init(int cellsCount, int blockDimension)
{
	blockDim = blockDimension;
	int n = cellsCount*blockDim;
	ilower = 0;
	iupper = n - 1;
	local_size = iupper - ilower + 1;

	cols	= new HYPRE_Int[blockDim];
	values	= new double[local_size];
	x		= new double[local_size];
	b		= new double[local_size];
	
	memset(x, 0, local_size*sizeof(double));
	memset(b, 0, local_size*sizeof(double));

	HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
	HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(A);
}

void SolverHypreCustomSeidel::zero() {
	HYPRE_IJMatrixDestroy(A);

	HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
	HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(A);

	memset(x, 0, local_size*sizeof(double));
	memset(b, 0, local_size*sizeof(double));
}



void SolverHypreCustomSeidel::printToFile(const char* fileName)
{
	HYPRE_Int n = local_size;
	HYPRE_Int nc = n;
	double * x = new double[n];
	HYPRE_Int * cols = new HYPRE_Int[n];
	FILE * fp = fopen(fileName, "w");
	for (int i = 0; i < nc; i++) {
		cols[i] = ilower + i;
	}

	for (HYPRE_Int row = 0; row < n; row++) {
		HYPRE_IJMatrixGetValues(A, 1, &nc, &row, cols, x);
		for (HYPRE_Int i = 0; i < nc; i++) {
			fprintf(fp, "%25.16e  ", x[i]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n\n=============================================================================================\n\n\n");
	//HYPRE_IJVectorGetValues(bb, local_size, cols, x);
	for (int i = 0; i < n; i++) {
		fprintf(fp, "%25.16e  \n", b[i]);
	}

	fclose(fp);
	delete[] x;
	delete[] cols;
}

