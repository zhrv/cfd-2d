#ifndef _GRID_H_
#define _GRID_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "global.h"

class Cell 
{
public:
	Cell(): nodesInd(NULL), edgesInd(NULL), nCount(0), eCount(0) {};
	~Cell();

	int   nCount;
	int   eCount;
	int*  nodesInd;
	int*  edgesInd;
	int	  neigh[3];
	int	  type;
	double  S;
	Point c; // центр €чейки
	double HX;
	double HY;
	friend class Grid;
};

class Edge 
{
public:
	Edge(): c(NULL), cCount(0) {};
	~Edge();

	int      n1;        // узел в начале
	int      n2;        // узел в конце
	int      c1;        // €чейка слева
	int      c2;        // €чейка справа, нормаль из с1 в с2
	Vector   n;         // нормаль к грани
	double     l;         // длина грани 
	double     cnl1;       // сумма длин проекций векторов из центра грани до центра смежной €чейки
	double     cnl2;       // сумма длин проекций векторов из центра грани до центра смежной €чейки
	int      cCount;    // количество точек на грани
	Point*   c;         // точки на грани
	int      type;      // тип грани (внутр., гранич.)
	friend class Grid;
public:
	static const int TYPE_INNER		= 0;  
	static const int TYPE_INLET		= 1;  
	static const int TYPE_OUTLET	= 2;  
	static const int TYPE_WALL		= 3;  
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
	int    cCount;  //количество €чеек(cells).
	int    eCount;  //количество ребер(edges).
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