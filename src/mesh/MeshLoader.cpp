#include "MeshLoader.h"
#include "grid.h"

MeshLoader* MeshLoader::create(const char* meshTypeName, Grid* grid)
{
	if (strcmp("berkeley_triangle", meshTypeName) == 0) {
		return NULL;
	}
}
