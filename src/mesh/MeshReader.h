#pragma once

#include "grid.h"

class MeshReader
{
public:
	static const int TYPE_BERKLEY_TRI	= 1;
	static const int TYPE_SALOME		= 2;

	virtual void read(Grid*) = 0;
	static MeshReader* create(int type, char* fileName);
	static int getType(char* name);
};

