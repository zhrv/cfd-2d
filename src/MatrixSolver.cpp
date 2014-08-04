#include "MatrixSolver.h"
#include "global.h"
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

void MatrixSolver::init(int cellsCount, int blockDimension)
{
	blockDim = blockDimension;
	int n = cellsCount*blockDim;
	a = new CSRMatrix(n);
	b = new double[n];
	x = new double[n];
}

void MatrixSolver::zero() {
	memset(x, 0, sizeof(double)*a->n);
	memset(b, 0, sizeof(double)*a->n);
	a->zero();
}

MatrixSolver::~MatrixSolver()
{
	delete a;
	delete[] b;
	delete[] x;
}

void MatrixSolver::setMatrElement(int i, int j, double** matrDim)
{
	for (int ii = 0; ii < blockDim; ++ii)
	{
		for (int jj = 0; jj < blockDim; ++jj)
		{
			a->set(ii+i*blockDim, jj+j*blockDim, matrDim[ii][jj]);
		}
	}
}

void MatrixSolver::setRightElement(int i, double* vectDim)
{
	for (int ii = 0; ii < blockDim; ++ii)
	{
		b[ii+i*blockDim] = vectDim[ii];
	}
}


void MatrixSolver::addMatrElement(int i, int j, double** matrDim)
{
	for (int ii = 0; ii < blockDim; ++ii)
	{
		for (int jj = 0; jj < blockDim; ++jj)
		{
			a->add(ii+i*blockDim, jj+j*blockDim, matrDim[ii][jj]);
		}
	}
}

void MatrixSolver::createMatrElement(int i, int j) {
	for (int ii = 0; ii < blockDim; ++ii)
	{
		for (int jj = 0; jj < blockDim; ++jj)
		{
			a->set(ii + i*blockDim, jj + j*blockDim, 0.0);
		}
	}

}

void MatrixSolver::addRightElement(int i, double* vectDim)
{
	for (int ii = 0; ii < blockDim; ii++)
	{
		b[ii+i*blockDim] += vectDim[ii];
	}
}

void MatrixSolver::printToFile(const char* fileName)
{
	a->printToFile(fileName);
}

int SolverZeidel::solve(double eps, int& maxIter)
{
	double	aii;
	double	err = 1.0;
	int		step = 0;
	double	tmp;
	//memset(x, 0, sizeof(double)*a->n);
	while(err > eps && step < maxIter)
	{
		step++;
		for (int i = 0; i < a->n; i++)
		{
			tmp = 0.0;
			aii = 0;
			for (int k = a->ia[i]; k < a->ia[i+1]; k++)
			{
				if (i == a->ja[k])
				{
					aii = a->a[k];
				} else {
					tmp += a->a[k]*x[a->ja[k]];
				}
			}
			if (fabs(aii) <= eps*eps) 
			{
				log("ZEIDEL_SOLVER: error: a[%d, %d] = 0\n", i, i);
				return MatrixSolver::RESULT_ERR_ZERO_DIAG;
			}
			x[i] = (-tmp+b[i])/aii;
		}
		err = 0.0;
		for (int i = 0; i < a->n; i++)
		{
			tmp = 0.0;
			for (int k = a->ia[i]; k < a->ia[i+1]; k++)
			{
				tmp += a->a[k]*x[a->ja[k]];
			}
			err += fabs(tmp-b[i]);
		}
		int qqqqq = 0; // ZHRV_WARN
	}
	if (step >= maxIter)
	{
		log("ZEIDEL_SOLVER: (warning) maximum iterations done (%d); error: %e\n", step, err);
		maxIter = step;
		return MatrixSolver::RESULT_ERR_MAX_ITER;
	}
	maxIter = step;
	return MatrixSolver::RESULT_OK;
}

int SolverJacobi::solve(double eps, int& maxIter)
{
	if (!tempXAlloc) {
		tempX = new double [a->n];
		tempXAlloc = true;
	}
	double	aii;
	double	err = 1.0;
	int		step = 0;
	double	tmp;
	//memset(x, 0, sizeof(double)*a->n);
	while(err > eps && step < maxIter)
	{
		step++;
		for (int i = 0; i < a->n; i++)
		{
			tmp = 0.0;
			aii = 0;
			for (int k = a->ia[i]; k < a->ia[i+1]; k++)
			{
				if (i == a->ja[k])	// i == j
				{
					aii = a->a[k];
				} else {
					tmp += a->a[k]*x[a->ja[k]];
				}
			}
			if (fabs(aii) <= eps*eps) 
			{
				log("JACOBI_SOLVER: error: a[%d, %d] = 0\n", i, i);
				return MatrixSolver::RESULT_ERR_ZERO_DIAG;
			}
			//x[i] = (-tmp+b[i])/aii;
			tempX[i] = (-tmp+b[i])/aii;
		}
		err = 0.0;
		for (int i = 0; i < a->n; i++)
		{
			tmp = 0.0;
			for (int k = a->ia[i]; k < a->ia[i+1]; k++)
			{
				x[a->ja[k]] = tempX[a->ja[k]];
				tmp += a->a[k]*x[a->ja[k]];
			}
			err += fabs(tmp-b[i]);
		}
		int qqqqq = 0; // ZHRV_WARN
	}
	if (step >= maxIter)
	{
		log("JACOBI_SOLVER: (warning) maximum iterations done (%d); error: %e\n", step, err);
		maxIter = step;
		return MatrixSolver::RESULT_ERR_MAX_ITER;
	}
	maxIter = step;
	return MatrixSolver::RESULT_OK;
}


int SolverHYPREBoomerAMG::solve(double eps, int& maxIter)
{
	HYPRE_IJMatrix A;
	HYPRE_ParCSRMatrix parcsr_A;
	HYPRE_IJVector bb;
	HYPRE_ParVector par_bb;
	HYPRE_IJVector xx;
	HYPRE_ParVector par_xx;

	HYPRE_Solver solver, precond;
	HYPRE_Int ilower, iupper, local_size;

	int solver_id = 0;

	ilower = 0;
	iupper = a->n-1;
	local_size = iupper - ilower+1;

	/* Create the matrix.
	Note that this is a square matrix, so we indicate the row partition
	size twice (since number of rows = number of cols) */
	HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);

	/* Choose a parallel csr format storage (see the User's Manual) */
	HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

	/* Initialize before setting coefficients */
	HYPRE_IJMatrixInitialize(A);

	double	aii;
	double	err = 1.0;
	int		step = 0;
	double	tmp;
	//memset(x, 0, sizeof(double)*a->n);
	{
		int nnz;
		double *values = new double[a->n];
		int *cols = new int[a->n];

		for (int i = 0; i < a->n; i++)
		{
			nnz = 0;
			for (int k = a->ia[i]; k < a->ia[i + 1]; k++)
			{
				cols[nnz] = a->ja[k];
				values[nnz] = a->a[k];
				nnz++;
			}
			HYPRE_IJMatrixSetValues(A, 1, &nnz, &i, cols, values);
		}

		/* Assemble after setting the coefficients */
		HYPRE_IJMatrixAssemble(A);

		delete[] values;
		delete[] cols;
	}

	/* Get the parcsr matrix object to use */
	HYPRE_IJMatrixGetObject(A, (void**)&parcsr_A);


	/* Create the rhs and solution */
	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &bb);
	HYPRE_IJVectorSetObjectType(bb, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(bb);

	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &xx);
	HYPRE_IJVectorSetObjectType(xx, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(xx);

	/* Set the rhs values to h^2 and the solution to zero */
	{
		int    *rows;

		rows = (int*)calloc(local_size, sizeof(int));

		for (int i = 0; i < local_size; i++)
		{
			rows[i] = ilower + i;
		}

		HYPRE_IJVectorSetValues(bb, local_size, rows, b);
		HYPRE_IJVectorSetValues(xx, local_size, rows, x);

		free(rows);
	}


	HYPRE_IJVectorAssemble(bb);
	/*  As with the matrix, for testing purposes, one may wish to read in a rhs:
	HYPRE_IJVectorRead( <filename>, MPI_COMM_WORLD,
	HYPRE_PARCSR, &b );
	as an alternative to the
	following sequence of HYPRE_IJVectors calls:
	Create, SetObjectType, Initialize, SetValues, and Assemble
	*/
	HYPRE_IJVectorGetObject(bb, (void **)&par_bb);

	HYPRE_IJVectorAssemble(xx);
	HYPRE_IJVectorGetObject(xx, (void **)&par_xx);


	/* Choose a solver and solve the system */

	/* AMG */
	if (solver_id == 0)
	{
		int num_iterations;
		double final_res_norm;

		/* Create solver */
		HYPRE_BoomerAMGCreate(&solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_BoomerAMGSetMinIter(solver, maxIter);
		HYPRE_BoomerAMGSetPrintLevel(solver, 3);  /* print solve info + parameters */
		HYPRE_BoomerAMGSetCoarsenType(solver, 6); /* Falgout coarsening */
		HYPRE_BoomerAMGSetRelaxType(solver, 3);   /* G-S/Jacobi hybrid relaxation */
		HYPRE_BoomerAMGSetNumSweeps(solver, 1);   /* Sweeeps on each level */
		HYPRE_BoomerAMGSetMaxLevels(solver, 20);  /* maximum number of levels */
		HYPRE_BoomerAMGSetTol(solver, eps);      /* conv. tolerance */

		/* Now setup and solve! */
		HYPRE_BoomerAMGSetup(solver, parcsr_A, par_bb, par_xx);
		HYPRE_BoomerAMGSolve(solver, parcsr_A, par_bb, par_xx);

		/* Run info - needed logging turned on */
		HYPRE_BoomerAMGGetNumIterations(solver, &maxIter);
		HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
		//if (myid == 0)
		{
			log("\n");
			log("Iterations = %d\n", maxIter);
			log("Final Relative Residual Norm = %e\n", final_res_norm);
			log("\n");
		}

		/* Destroy solver */
		HYPRE_BoomerAMGDestroy(solver);
	}
	/* PCG */
	else if (solver_id == 50)
	{
		int num_iterations;
		double final_res_norm;

		/* Create solver */
		HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_PCGSetMaxIter(solver, 1000); /* max iterations */
		HYPRE_PCGSetTol(solver, 1e-7); /* conv. tolerance */
		HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
		HYPRE_PCGSetPrintLevel(solver, 2); /* prints out the iteration info */
		HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

		/* Now setup and solve! */
		HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_bb, par_xx);
		HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_bb, par_xx);

		/* Run info - needed logging turned on */
		HYPRE_PCGGetNumIterations(solver, &num_iterations);
		HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
		//if (myid == 0)
		{
			log("\n");
			log("Iterations = %d\n", num_iterations);
			log("Final Relative Residual Norm = %e\n", final_res_norm);
			log("\n");
		}

		/* Destroy solver */
		HYPRE_ParCSRPCGDestroy(solver);
	}

	{
		int *rows = (int*)calloc(a->n, sizeof(int));
		for (int i = 0; i < a->n; i++)
			rows[i] = ilower + i;

		/* get the local solution */
		HYPRE_IJVectorGetValues(xx, a->n, rows, x);
		
		delete[] rows;
	}


	return MatrixSolver::RESULT_OK;
}