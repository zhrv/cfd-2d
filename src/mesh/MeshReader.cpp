#include "MeshReader.h"
#include "MeshReaderBerkleyTriangle.h"
#include "MeshReaderSalomeUnv.h"
#include "MeshReaderGmsh.h"



MeshReader* MeshReader::create(int type, char* fileName)
{
	switch (type) {
        case TYPE_BERKLEY_TRI:
            return new MeshReaderBerkleyTriangle(fileName);
            break;
        case TYPE_SALOME:
            return new MeshReaderSalomeUnv(fileName);
            break;
        case TYPE_GMSH:
            return new MeshReaderGmsh(fileName);
            break;
	}
}

int MeshReader::getType(char* name) {
	if (strcmp(name, "berkeley_triangle") == 0) return MeshReader::TYPE_BERKLEY_TRI;
    if (strcmp(name, "salome_unv") == 0) return MeshReader::TYPE_SALOME;
    if (strcmp(name, "gmsh") == 0) return MeshReader::TYPE_GMSH;
	throw Exception("Wrong mesh type.", Exception::TYPE_MESH_WRONG_NAME);
}
