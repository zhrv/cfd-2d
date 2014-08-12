#include "MatrixSolver.h"
#include "global.h"
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"


void SolverHYPREBoomerAMG::initMatrVectors()
{
	/* Create the matrix.
	Note that this is a square matrix, so we indicate the row partition
	size twice (since number of rows = number of cols) */
	HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);

	/* Choose a parallel csr format storage (see the User's Manual) */
	HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

	/* Initialize before setting coefficients */
	HYPRE_IJMatrixInitialize(A);

	/* Create the rhs and solution */
	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &bb);
	HYPRE_IJVectorSetObjectType(bb, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(bb);

	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &xx);
	HYPRE_IJVectorSetObjectType(xx, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(xx);
}

void SolverHYPREBoomerAMG::init(int cellsCount, int blockDimension)
{
	blockDim = blockDimension;
	int n = cellsCount*blockDim;
	ilower = 0;
	iupper = n-1;
	local_size = iupper - ilower + 1;

	cols = new int[blockDim];
	x = new double[n];
	
	initMatrVectors();

}

void SolverHYPREBoomerAMG::zero() {
	HYPRE_IJMatrixDestroy(A);
	HYPRE_IJVectorDestroy(bb);
	HYPRE_IJVectorDestroy(xx);
	
	initMatrVectors();
}

SolverHYPREBoomerAMG::~SolverHYPREBoomerAMG()
{
	HYPRE_IJMatrixDestroy(A);
	HYPRE_IJVectorDestroy(bb);
	HYPRE_IJVectorDestroy(xx);
	delete[] x;
}

void SolverHYPREBoomerAMG::setMatrElement(int i, int j, double** matrDim)
{
	int row, nnz;

	for (int ii = 0; ii < blockDim; ++ii)
	{
		row = ii + i*blockDim;
		for (int jj = 0; jj < blockDim; ++jj)
		{
			x[jj] = matrDim[ii][jj];
			cols[jj] = jj + j*blockDim;
		}
		HYPRE_IJMatrixSetValues(A, 1, &blockDim, &row, cols, x);
	}
	
	
}

void SolverHYPREBoomerAMG::setRightElement(int i, double* vectDim)
{
	for (int ii = 0; ii < blockDim; ++ii)
	{
		cols[ii] = ii + i*blockDim;
	}
	HYPRE_IJVectorSetValues(bb, blockDim, cols, vectDim);
}


void SolverHYPREBoomerAMG::addMatrElement(int i, int j, double** matrDim)
{
	int row, nnz;

	for (int ii = 0; ii < blockDim; ++ii)
	{
		row = ii + i*blockDim;
		for (int jj = 0; jj < blockDim; ++jj)
		{
			x[jj] = matrDim[ii][jj];
			cols[jj] = jj + j*blockDim;
		}
		HYPRE_IJMatrixAddToValues(A, 1, &blockDim, &row, cols, x);
	}

}

void SolverHYPREBoomerAMG::createMatrElement(int i, int j) {
}

void SolverHYPREBoomerAMG::addRightElement(int i, double* vectDim)
{
	for (int ii = 0; ii < blockDim; ++ii)
	{
		cols[ii] = ii + i*blockDim;
	}
	HYPRE_IJVectorAddToValues(bb, blockDim, cols, vectDim);
}

int SolverHYPREBoomerAMG::solve(double eps, int& maxIter)
{
	/* Set the solution to zero */
	{
		int    *rows;

		rows = (int*)calloc(local_size, sizeof(int));

		for (int i = 0; i < local_size; i++)
		{
			rows[i] = ilower + i;
			x[i] = 0.0;
		}

		
		HYPRE_IJVectorSetValues(xx, local_size, rows, x);

		free(rows);
	}


	/* Assemble after setting the coefficients */
	HYPRE_IJMatrixAssemble(A);
	HYPRE_IJVectorAssemble(bb);
	HYPRE_IJVectorAssemble(xx);


	/* Get the parcsr matrix object to use */
	HYPRE_IJMatrixGetObject(A, (void**)&parcsr_A);
	HYPRE_IJVectorGetObject(bb, (void **)&par_bb);
	HYPRE_IJVectorGetObject(xx, (void **)&par_xx);

	////printToFile("matr.txt");
	//{
	//	int *rows = (int*)calloc(local_size, sizeof(int));
	//	for (int i = 0; i < local_size; i++)
	//		rows[i] = ilower + i;

	//	/* get the local solution */
	//	HYPRE_IJVectorGetValues(bb, local_size, rows, x);

	//	delete[] rows;
	//}



	/* Choose a solver and solve the system */

	/* AMG */
	if (solver_id == 0)
	{
		double final_res_norm;

		/* Create solver */
		HYPRE_BoomerAMGCreate(&solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_BoomerAMGSetMaxIter(solver, maxIter);
		HYPRE_BoomerAMGSetPrintLevel(solver, 2);  /* print solve info + parameters */
		HYPRE_BoomerAMGSetCoarsenType(solver, 6); /* Falgout coarsening */
		HYPRE_BoomerAMGSetRelaxType(solver, 3);   /* G-S/Jacobi hybrid relaxation */
		HYPRE_BoomerAMGSetNumSweeps(solver, 1);   /* Sweeeps on each level */
		HYPRE_BoomerAMGSetMaxLevels(solver, 20);  /* maximum number of levels */
		HYPRE_BoomerAMGSetTol(solver, eps);       /* conv. tolerance */

		/* Now setup and solve! */
		HYPRE_BoomerAMGSetup(solver, parcsr_A, par_bb, par_xx);
		HYPRE_BoomerAMGSolve(solver, parcsr_A, par_bb, par_xx);

		/* Run info - needed logging turned on */
		HYPRE_BoomerAMGGetNumIterations(solver, &maxIter);
		HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
		//if (myid == 0)
		/*{
		log("\n");
		log("Iterations = %d\n", maxIter);
		log("Final Relative Residual Norm = %e\n", final_res_norm);
		log("\n");
		}*/

		/* Destroy solver */
		HYPRE_BoomerAMGDestroy(solver);
	}
	/* PCG */
	else if (solver_id == 50)
	{
		double final_res_norm;

		/* Create solver */
		HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_PCGSetMaxIter(solver, maxIter); /* max iterations */
		HYPRE_PCGSetTol(solver, eps); /* conv. tolerance */
		HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
		HYPRE_PCGSetPrintLevel(solver, 1); /* prints out the iteration info */
		HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

		/* Now setup and solve! */
		HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_bb, par_xx);
		HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_bb, par_xx);

		/* Run info - needed logging turned on */
		HYPRE_PCGGetNumIterations(solver, &maxIter);
		HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
		//if (myid == 0)
		{
			log("\n");
			log("Iterations = %d\n", maxIter);
			log("Final Relative Residual Norm = %e\n", final_res_norm);
			log("\n");
		}

		/* Destroy solver */
		HYPRE_ParCSRPCGDestroy(solver);
	}
	/* PCG with AMG preconditioner */
	else if (solver_id == 1)
	{
		int num_iterations;
		double final_res_norm;

		/* Create solver */
		HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_PCGSetMaxIter(solver, maxIter); /* max iterations */
		HYPRE_PCGSetTol(solver, 1e-7); /* conv. tolerance */
		HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
		HYPRE_PCGSetPrintLevel(solver, 1); /* print solve info */
		HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

		/* Now set up the AMG preconditioner and specify any parameters */
		HYPRE_BoomerAMGCreate(&precond);
		HYPRE_BoomerAMGSetPrintLevel(precond, 1); /* print amg solution info */
		HYPRE_BoomerAMGSetCoarsenType(precond, 6);
		HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
		HYPRE_BoomerAMGSetNumSweeps(precond, 1);
		HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
		HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */

		/* Set the PCG preconditioner */
		HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
			(HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);

		/* Now setup and solve! */
		HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_bb, par_xx);
		HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_bb, par_xx);

		/* Run info - needed logging turned on */
		HYPRE_PCGGetNumIterations(solver, &maxIter);
		HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
		//if (myid == 0)
		{
			printf("\n");
			printf("Iterations = %d\n", maxIter);
			printf("Final Relative Residual Norm = %e\n", final_res_norm);
			printf("\n");
		}

		/* Destroy solver and preconditioner */
		HYPRE_ParCSRPCGDestroy(solver);
		HYPRE_BoomerAMGDestroy(precond);
	}
	/* PCG with Parasails Preconditioner */
	else if (solver_id == 8)
	{
		int    num_iterations;
		double final_res_norm;

		int      sai_max_levels = 1;
		double   sai_threshold = 0.1;
		double   sai_filter = 0.05;
		int      sai_sym = 1;

		/* Create solver */
		HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_PCGSetMaxIter(solver, maxIter); /* max iterations */
		HYPRE_PCGSetTol(solver, 1e-7); /* conv. tolerance */
		HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
		HYPRE_PCGSetPrintLevel(solver, 2); /* print solve info */
		HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

		/* Now set up the ParaSails preconditioner and specify any parameters */
		HYPRE_ParaSailsCreate(MPI_COMM_WORLD, &precond);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_ParaSailsSetParams(precond, sai_threshold, sai_max_levels);
		HYPRE_ParaSailsSetFilter(precond, sai_filter);
		HYPRE_ParaSailsSetSym(precond, sai_sym);
		HYPRE_ParaSailsSetLogging(precond, 3);

		/* Set the PCG preconditioner */
		HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_ParaSailsSolve,
			(HYPRE_PtrToSolverFcn)HYPRE_ParaSailsSetup, precond);

		/* Now setup and solve! */
		HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_bb, par_xx);
		HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_bb, par_xx);


		/* Run info - needed logging turned on */
		HYPRE_PCGGetNumIterations(solver, &maxIter);
		HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
		//if (myid == 0)
		{
			printf("\n");
			printf("Iterations = %d\n", maxIter);
			printf("Final Relative Residual Norm = %e\n", final_res_norm);
			printf("\n");
		}

		/* Destory solver and preconditioner */
		HYPRE_ParCSRPCGDestroy(solver);
		HYPRE_ParaSailsDestroy(precond);
	}
	/* Flexible GMRES with  AMG Preconditioner */
	else if (solver_id == 61)
	{
		int    num_iterations;
		double final_res_norm;
		int    restart = 30;
		int    modify = 1;


		/* Create solver */
		HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, &solver);

		/* Set some parameters (See Reference Manual for more parameters) */
		HYPRE_FlexGMRESSetKDim(solver, restart);
		HYPRE_FlexGMRESSetMaxIter(solver, maxIter); /* max iterations */
		HYPRE_FlexGMRESSetTol(solver, 1e-7); /* conv. tolerance */
		HYPRE_FlexGMRESSetPrintLevel(solver, 2); /* print solve info */
		HYPRE_FlexGMRESSetLogging(solver, 1); /* needed to get run info later */


		/* Now set up the AMG preconditioner and specify any parameters */
		HYPRE_BoomerAMGCreate(&precond);
		HYPRE_BoomerAMGSetPrintLevel(precond, 1); /* print amg solution info */
		HYPRE_BoomerAMGSetCoarsenType(precond, 6);
		HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
		HYPRE_BoomerAMGSetNumSweeps(precond, 1);
		HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero */
		HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */

		/* Set the FlexGMRES preconditioner */
		//HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
		//	(HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);




		/* Now setup and solve! */
		HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_bb, par_xx);
		HYPRE_ParCSRFlexGMRESSolve(solver, parcsr_A, par_bb, par_xx);

		/* Run info - needed logging turned on */
		HYPRE_FlexGMRESGetNumIterations(solver, &maxIter);
		HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
		//if (myid == 0)
		{
			printf("\n");
			printf("Iterations = %d\n", maxIter);
			printf("Final Relative Residual Norm = %e\n", final_res_norm);
			printf("\n");
		}

		/* Destory solver and preconditioner */
		HYPRE_ParCSRFlexGMRESDestroy(solver);
		HYPRE_BoomerAMGDestroy(precond);

	}
	else
	{
		printf("Invalid solver id specified.\n");
	}

	{
		int *rows = (int*)calloc(local_size, sizeof(int));
		for (int i = 0; i < local_size; i++)
			rows[i] = ilower + i;

		/* get the local solution */
		HYPRE_IJVectorGetValues(xx, local_size, rows, x);

		delete[] rows;
	}


	return MatrixSolver::RESULT_OK;
}


void SolverHYPREBoomerAMG::printToFile(const char* fileName)
{
	int n = local_size;
	int nc = 1;
	double * x = new double[n];
	int * cols = new int[n];
	FILE * fp = fopen(fileName, "w");
	for (int i = 0; i < n; i++) {
		cols[i] = ilower + i;
	}
	 
	for (int row = 0; row < n; row++) {
		HYPRE_IJMatrixGetValues(A, 1, &nc, &row, &row, x);
		for (int i = 0; i < nc; i++) {
			fprintf(fp, "%25.16e  ", x[i]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n\n=============================================================================================\n\n\n");
	HYPRE_IJVectorGetValues(xx, local_size, cols, x);
	for (int i = 0; i < n; i++) {
		fprintf(fp, "%25.16e  ", x[i]);
	}
	
	fclose(fp);
	delete[] x;
	delete[] cols;
}