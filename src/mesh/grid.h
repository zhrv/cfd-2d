#ifndef _GRID_H_
#define _GRID_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>
#include "global.h"
#include "bnd_cond.h"

const unsigned int CELL_FLAG_GOOD	= 0x000000;
const unsigned int CELL_FLAG_BAD	= 0x000001;
const unsigned int CELL_FLAG_LIM	= 0x000002;



class Cell 
{
public:
	Cell(): nCount(0), eCount(0), nodesInd(NULL), edgesInd(NULL) {};
	~Cell();

	int   nCount;
	int   eCount;
	int*  nodesInd;
	int*  edgesInd;
	int	  neigh[3];
	int	  type;
	char  typeName[64];
	double  S;
	Point c; // центр ячейки
	double HX;
	double HY;
	unsigned int flag;

	friend class Grid;
};

class Edge 
{
public:
	Edge(): cCount(0), c(nullptr) {};
	~Edge();

	int      n1;        // узел в начале
	int      n2;        // узел в конце
	int      c1;        // ячейка слева
	int      c2;        // ячейка справа, нормаль из с1 в с2
	Vector   n;         // нормаль к грани
	double     l;         // длина грани 
	double     cnl1;       // сумма длин проекций векторов из центра грани до центра смежной ячейки
	double     cnl2;       // сумма длин проекций векторов из центра грани до центра смежной ячейки
	int      cCount;    // количество точек на грани
	Point*   c;         // точки на грани
	int      type;      // тип грани (внутр., гранич.)
	char	typeName[64];
	CFDBoundary * bnd = NULL;
	friend class Grid;
public:
	static const int TYPE_INNER		= 0;  
	static const int TYPE_INLET		= 1;  
	static const int TYPE_OUTLET	= 2;  
	static const int TYPE_WALL		= 3;  

	static const int TYPE_NAMED = 0x100;
};

class Grid 
{
public:
	Grid(): nodes(NULL), cells(NULL), edges(NULL),
		nCount(0), cCount(0), eCount(0) {};
	~Grid();

	void initFromFiles(char* fName);
	inline Point& getNode(int i) { return nodes[i]; };
	inline Cell&  getCell(int i) { return cells[i]; };
	inline Edge&  getEdge(int i) { return edges[i]; };
	void replaceEdges(int if1, int if2);
	virtual void reorderEdges();
	int findEdge(int n1, int n2);

	Point* nodes;
	Cell*  cells;
	Edge*  edges;
	int    nCount;	//количество nodes.
	int    cCount;  //количество ячеек(cells).
	int    eCount;  //количество ребер(edges).
	int    nCountEx;	//количество nodes с фиктивными ячйками.
	int    cCountEx;  //количество ячеек(cells) с фиктивными ячйками.
	int    eCountEx;  //количество ребер(edges) с фиктивными ячйками.

	std::vector<int>					recvCount;
	std::vector<int>					recvShift;
	std::vector< std::vector<int> >		sendInd;

    void readMeshFiles();
    void procMeshFiles();
	void saveMeshInfo();
};

inline double _max_(double a, double b)
{
	if (a>b) { return a; } else {return b; }
}

inline double _min_(double a, double b)
{
	if (a<b) { return a; } else {return b; }
}

inline double _max_(double a, double b, double c)
{
	return _max_(a, _max_(b, c));
}

inline double _min_(double a, double b, double c)
{
	return _min_(a, _min_(b, c));
}

inline double _sqr_(double a)
{
	return a*a;
}

#endif