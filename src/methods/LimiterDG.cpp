#include "LimiterDG.h"
#include <cstring>
#include "LimiterDGCockburn.h"

LimiterDG::~LimiterDG()
{
}

LimiterDG* LimiterDG::create(const char* limiterName, FEM_DG_IMPLICIT* solver)
{
	if (strcmp(limiterName, "Cockburn") == 0) {
		return new LimiterDGCockburn(solver, &(solver->grid), solver->ro, solver->ru, solver->rv, solver->re, solver->BASE_FUNC_COUNT);
	}
	else {
		return NULL;
	}
}
