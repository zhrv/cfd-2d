#pragma once

#include "fem_dg.h"

class LimiterDG
{
public:
	virtual ~LimiterDG();
	virtual void run() = 0;
	virtual const char* getName() = 0;
	static LimiterDG* create(const char *limiterName, FEM_DG *solver);
};

