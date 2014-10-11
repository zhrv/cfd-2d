#pragma once

#include "grid.h"
#include <cstring>

class MeshLoader
{
public:
	virtual void load() = 0;
	MeshLoader* create(const char* meshTypeName, Grid* grid);
};

