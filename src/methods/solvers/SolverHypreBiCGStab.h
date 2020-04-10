#pragma once
#include "SolverHypre.h"
class SolverHypreBiCGStab  : public SolverHypre
{
public:
    virtual int solve(double eps, int& maxIter);
    virtual char* getName() { return (char*)"HYPRE BiCGSTAB"; }
    void setParameter(const char* name, int val);
private:
    int KRYLOV_DIM = 30;
};

