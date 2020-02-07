#pragma once
#include "MeshReader.h"
#include <vector>
#include <string>
#include <map>
#include <set>

using namespace std;


class MeshReaderMsh2 : public MeshReader
{
public:
    struct element {
        int ind;
        int type;

        static const int TYPE_EDGE = 1;
        static const int TYPE_CELL = 2;

        element() : ind(0), type(0){}
        element(int i, int t) : ind(i), type(t) {}
    };

    typedef vector<string> string_list;
    typedef vector<int> indexes;
    typedef vector<indexes> index_list;
    typedef map<string, indexes> bnd_map;
    typedef set<int> ind_set;
    typedef vector<element> ele_map;
    typedef std::pair<int, int> edge_t;
    typedef std::map<edge_t, indexes> edge_cells_t;

private:
    char* fileName;
    index_list edges, cells;
    bnd_map bounds;
    vector<Point> points;
    ele_map elements;

    int find_edge(int n1, int n2);

public:
    MeshReaderMsh2(char* fName) : fileName(fName) {}
    ~MeshReaderMsh2() { delete[] fileName; }

    virtual void read(Grid*);
};



