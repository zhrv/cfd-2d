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
    void addNeigh(Cell &c, int iNeigh);

public:
    MeshReaderMsh2(char* fName) : fileName(fName) {}
    ~MeshReaderMsh2() { delete[] fileName; }

    virtual void read(Grid*);

    struct element {
        int ind;
        int type;

        static const int TYPE_EDGE = 11;
        static const int TYPE_CELL = 41;

        element() : ind(0), type(0){}
        element(int i, int t) : ind(i), type(t) {}
    };

    typedef vector<string> string_list;
    typedef vector<int> indexes;
    typedef vector<indexes> index_list;
    typedef map<string, indexes> bnd_map;
    typedef set<int> ind_set;
    typedef vector<element> ele_map;

};



