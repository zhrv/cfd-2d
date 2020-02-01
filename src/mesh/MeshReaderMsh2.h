#pragma once
#include "MeshReader.h"
#include <vector>
#include <string>
#include <map>
#include <set>

using namespace std;


class MeshReaderMsh2 :
        public MeshReader
{
private:
    char* fileName;

    char* getLineUpper (ifstream &stream);

public:
    MeshReaderMsh2(char* fName) : fileName(fName) {}
    ~MeshReaderMsh2() { delete[] fileName; }

    virtual void read(Grid*);

};



