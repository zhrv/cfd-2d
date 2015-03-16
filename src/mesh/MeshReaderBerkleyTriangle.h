#pragma once
#include "MeshReader.h"
class MeshReaderBerkleyTriangle :
	public MeshReader
{
private:
	char* fileName;
public:
	MeshReaderBerkleyTriangle(char* fName): fileName(fName) {}
	~MeshReaderBerkleyTriangle() { delete[] fileName; }
	
	virtual void read(Grid*);
};

