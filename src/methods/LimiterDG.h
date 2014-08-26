#pragma once

#include "FEM_DG.h"

class LimiterDG
{
public:
	LimiterDG();
	virtual ~LimiterDG();
	virtual void run() = 0;
	virtual const char* getName() = 0;
	static LimiterDG* create(const char* limiterName, FEM_DG* solver);
};

