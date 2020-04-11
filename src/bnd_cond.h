#ifndef _BND_COND_H_
#define _BND_COND_H_

#include <cstdlib>
#include "global.h"
#include "tinyxml.h"

#include <vector>

class Edge;
class Grid;
struct CFDBoundary
{
	static CFDBoundary* create(TiXmlNode* bNode, Grid * g);
	virtual void run(int iEdgs, Param& pL, Param& pR) = 0;

	CFDBoundary() : par(NULL) {}
	~CFDBoundary() { if (par) delete[] par; par = NULL; }


	int			edgeType;
	char		name[64];
	double *	par;
	int			parCount;
	Grid *		g;

	static const char* TYPE_INLET;
	static const char* TYPE_OUTLET;
    static const char* TYPE_PRESSURE;
	static const char* TYPE_WALL_SLIP;
    static const char* TYPE_WALL_NO_SLIP;

};
typedef std::vector< CFDBoundary* > CFDBoundaries;

struct CFDBndInlet : public CFDBoundary
{
	virtual void run(int iEdge, Param& pL, Param& pR);
};

struct CFDBndOutlet : public CFDBoundary
{
	virtual void run(int iEdge, Param& pL, Param& pR);
};

struct CFDBndWallSlip : public CFDBoundary
{
	virtual void run(int iEdge, Param& pL, Param& pR);
};

struct CFDBndWallNoSlip : public CFDBoundary
{
    virtual void run(int iEdge, Param& pL, Param& pR);
};

struct CFDBndPressure : public CFDBoundary
{
    virtual void run(int iEdge, Param& pL, Param& pR);
};


#endif
