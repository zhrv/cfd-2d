#include "bnd_cond.h"
#include "grid.h"

const char* CFDBoundary::TYPE_INLET = "BOUND_INLET";
const char* CFDBoundary::TYPE_OUTLET = "BOUND_OUTLET";
const char* CFDBoundary::TYPE_WALL_SLIP = "BOUND_WALL_SLIP";
const char* CFDBoundary::TYPE_WALL_NO_SLIP = "BOUND_WALL_NO_SLIP";


CFDBoundary* CFDBoundary::create(TiXmlNode* bNode, Grid * g)
{
	CFDBoundary * b = 0;
	TiXmlNode* node = 0;
	TiXmlNode* node1 = 0;

	const char * type = bNode->FirstChild("type")->ToElement()->GetText();

	if (strcmp(type, CFDBoundary::TYPE_INLET) == 0) {
		b = new CFDBndInlet();
		bNode->ToElement()->Attribute("edgeType", &b->edgeType);
		b->parCount = 4;
		b->par = new double[4];
		node = bNode->FirstChild("parameters");

		node1 = node->FirstChild("Vx");
		if (!node1) throw Exception((char*)"Parameter 'Vx' isn't specified  for BOUND_INLET.", Exception::TYPE_BOUND_NOPAR);
		node1->ToElement()->Attribute("value", &b->par[0]);

		node1 = node->FirstChild("Vy");
		if (!node1) throw Exception((char*)"Parameter 'Vy' isn't specified  for BOUND_INLET.", Exception::TYPE_BOUND_NOPAR);
		node1->ToElement()->Attribute("value", &b->par[1]);

		node1 = node->FirstChild("T");
		if (!node1) throw Exception((char*)"Parameter 'T' isn't specified  for BOUND_INLET.", Exception::TYPE_BOUND_NOPAR);
		node1->ToElement()->Attribute("value", &b->par[2]);

		node1 = node->FirstChild("P");
		if (!node1) throw Exception((char*)"Parameter 'P' isn't specified  for BOUND_INLET.", Exception::TYPE_BOUND_NOPAR);
		node1->ToElement()->Attribute("value", &b->par[3]);

	}

	if (strcmp(type, CFDBoundary::TYPE_OUTLET) == 0) {
		b = new CFDBndOutlet();
		bNode->ToElement()->Attribute("edgeType", &b->edgeType);
		b->parCount = 0;
		b->par = NULL;
	}

	if (strcmp(type, CFDBoundary::TYPE_WALL_NO_SLIP) == 0) {
		b = new CFDBndWallSlip();
		bNode->ToElement()->Attribute("edgeType", &b->edgeType);
		b->parCount = 0;
		b->par = NULL;
	}

	if (strcmp(type, CFDBoundary::TYPE_WALL_SLIP) == 0) {
		b = new CFDBndWallSlip();
		bNode->ToElement()->Attribute("edgeType", &b->edgeType);
		b->parCount = 0;
		b->par = NULL;
	}

	if (!b) {
		throw Exception((char*)"Unknown boundary type '%s' specified.", Exception::TYPE_BOUND_UNKNOWN);
	}

	const char * name = bNode->FirstChild("name")->ToElement()->GetText();
	strcpy(b->name, name);
	b->g = g;
	return b;
}


void CFDBndInlet::run(int iEdge, Param& pL, Param& pR)
{
	pR.u = par[0];
	pR.v = par[1];
	pR.T = par[2];
	pR.p = par[3];
}

void CFDBndOutlet::run(int iEdge, Param& pL, Param& pR) // @todo исправить - необходимо учитывать шаг сетки
{
	pR = pL;
}

void CFDBndWallNoSlip::run(int iEdge, Param& pL, Param& pR)
{
	Edge &e = g->edges[iEdge];
	pR = pL;
	double Un = pL.u*e.n.x + pL.v*e.n.y;
	Vector V;
	V.x = e.n.x*Un*2.0;
	V.y = e.n.y*Un*2.0;
	pR.u = pL.u - V.x;
	pR.v = pL.v - V.y;
}

void CFDBndWallSlip::run(int iEdge, Param& pL, Param& pR)
{
	Edge &e = g->edges[iEdge];
	pR = pL;
	double Un = pL.u*e.n.x + pL.v*e.n.y;
	Vector V;
	V.x = e.n.x*Un*2.0;
	V.y = e.n.y*Un*2.0;
	pR.u = pL.u - V.x;
	pR.v = pL.v - V.y;
}

