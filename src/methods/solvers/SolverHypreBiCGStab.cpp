#include "SolverHypreBiCGStab.h"
#include "global.h"
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

int SolverHypreBiCGStab::solve(double eps, int& maxIter)
{
    int result = MatrixSolver::RESULT_OK;

    /* Set the solution to zero */
    int    *rows;

    rows = (int*)calloc(local_size, sizeof(int));

    for (int i = 0; i < local_size; i++)
    {
        rows[i] = ilower + i;
        x[i] = 0.0;
    }
    free(rows);

    HYPRE_IJVectorSetValues(xx, local_size, (HYPRE_Int*)rows, x);


    /* Assemble after setting the coefficients */
    HYPRE_IJMatrixAssemble(A);
    HYPRE_IJVectorAssemble(bb);
    HYPRE_IJVectorAssemble(xx);



    /* Get the parcsr matrix object to use */
    HYPRE_IJMatrixGetObject(A, (void**)&parcsr_A);
    HYPRE_IJVectorGetObject(bb, (void **)&par_bb);
    HYPRE_IJVectorGetObject(xx, (void **)&par_xx);


    /* Choose a solver and solve the system */

    /* BiCGSTAB */
    {
        double final_res_norm;

        /* Create solver */
        HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &solver);

        /* Set some parameters (See Reference Manual for more parameters) */
        HYPRE_BiCGSTABSetMaxIter(solver,		maxIter);			/* max iterations */
        HYPRE_BiCGSTABSetTol(solver,			eps);				/* conv. tolerance */
        HYPRE_BiCGSTABSetAbsoluteTol(solver,	0.0);
        HYPRE_BiCGSTABSetPrintLevel(solver,	PRINT_LEVEL);		/* print solve info */
        HYPRE_BiCGSTABSetLogging(solver,		1);					/* needed to get run info later */


        /* Now setup and solve! */
        HYPRE_ParCSRBiCGSTABSetup(solver, parcsr_A, par_bb, par_xx);
        HYPRE_ParCSRBiCGSTABSolve(solver, parcsr_A, par_bb, par_xx);


        /* Run info - needed logging turned on */
        int initMaxIter = maxIter;
        HYPRE_BiCGSTABGetNumIterations(solver, (HYPRE_Int*)&maxIter);
        HYPRE_BiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
        if (initMaxIter <= maxIter) {
            result |= MatrixSolver::RESULT_ERR_MAX_ITER;
        }
        if (final_res_norm >= eps || !std::isfinite(final_res_norm)) {
            result |= MatrixSolver::RESULT_ERR_CONVERG;
        }

        /* Destory solver */
        HYPRE_ParCSRBiCGSTABDestroy(solver);
    }

    if (result == MatrixSolver::RESULT_OK) {
        int *rows = (int*)calloc(local_size, sizeof(int));
        for (int i = 0; i < local_size; i++)
            rows[i] = ilower + i;

        /* get the local solution */
        HYPRE_IJVectorGetValues(xx, local_size, (HYPRE_Int*)rows, x);

        delete[] rows;
    }


    return result;
}


void SolverHypreBiCGStab::setParameter(const char* name, int val)
{
    if (strcmp(name, "KRYLOV_DIM") == 0) {
        KRYLOV_DIM = val;
    }
    else {
        SolverHypre::setParameter(name, val);
    }
}

