#include "bnd_cond.h"
#include "grid.h"

const char* CFDBoundary::TYPE_INLET         = "BOUND_INLET";
const char* CFDBoundary::TYPE_OUTLET        = "BOUND_OUTLET";
const char* CFDBoundary::TYPE_PRESSURE      = "BOUND_PRESSURE";
const char* CFDBoundary::TYPE_WALL_SLIP     = "BOUND_WALL_SLIP";
const char* CFDBoundary::TYPE_WALL_NO_SLIP  = "BOUND_WALL_NO_SLIP";

const int CFDBoundary::TYPE_ID_INLET            = 0;
const int CFDBoundary::TYPE_ID_OUTLET           = 1;
const int CFDBoundary::TYPE_ID_PRESSURE         = 2;
const int CFDBoundary::TYPE_ID_WALL_SLIP        = 3;
const int CFDBoundary::TYPE_ID_WALL_NO_SLIP     = 4;


CFDBoundary* CFDBoundary::create(TiXmlNode* bNode, Grid * g)
{
	CFDBoundary * b = nullptr;
	TiXmlNode* node = nullptr;
	TiXmlNode* node1 = nullptr;

	const char * type = bNode->FirstChild("type")->ToElement()->GetText();

	if (strcmp(type, CFDBoundary::TYPE_INLET) == 0) {
		b = new CFDBndInlet();
		//bNode->ToElement()->Attribute("edgeType", &b->edgeType);
        b->edgeType = CFDBoundary::TYPE_ID_INLET;
		b->parCount = 4;
		b->par = new double[4];
		node = bNode->FirstChild("parameters");

		node1 = node->FirstChild("Vx");
		if (!node1) throw Exception("Parameter 'Vx' isn't specified  for BOUND_INLET.", Exception::TYPE_BOUND_NOPAR);
		node1->ToElement()->Attribute("value", &b->par[0]);

		node1 = node->FirstChild("Vy");
		if (!node1) throw Exception("Parameter 'Vy' isn't specified  for BOUND_INLET.", Exception::TYPE_BOUND_NOPAR);
		node1->ToElement()->Attribute("value", &b->par[1]);

		node1 = node->FirstChild("T");
		if (!node1) throw Exception("Parameter 'T' isn't specified  for BOUND_INLET.", Exception::TYPE_BOUND_NOPAR);
		node1->ToElement()->Attribute("value", &b->par[2]);

		node1 = node->FirstChild("P");
		if (!node1) throw Exception("Parameter 'P' isn't specified  for BOUND_INLET.", Exception::TYPE_BOUND_NOPAR);
		node1->ToElement()->Attribute("value", &b->par[3]);
	}

	if (strcmp(type, CFDBoundary::TYPE_OUTLET) == 0) {
		b = new CFDBndOutlet();
		//bNode->ToElement()->Attribute("edgeType", &b->edgeType);
        b->edgeType = CFDBoundary::TYPE_ID_OUTLET;
        b->parCount = 1;
        b->par = new double[1];
        node = bNode->FirstChild("parameters");

        node1 = node->FirstChild("P");
        if (!node1) throw Exception("Parameter 'P' isn't specified  for BOUND_INLET.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[0]);
	}

	if (strcmp(type, CFDBoundary::TYPE_WALL_NO_SLIP) == 0) {
		b = new CFDBndWallNoSlip();
		//bNode->ToElement()->Attribute("edgeType", &b->edgeType);
        b->edgeType = CFDBoundary::TYPE_ID_WALL_NO_SLIP;
        b->parCount = 1;
        b->par = new double[1];
        node = bNode->FirstChild("parameters");

        node1 = node->FirstChild("T");
        if (!node1) throw Exception("Parameter 'T' isn't specified  for BOUND_INLET.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[0]);
	}

    if (strcmp(type, CFDBoundary::TYPE_WALL_SLIP) == 0) {
        b = new CFDBndWallSlip();
        //bNode->ToElement()->Attribute("edgeType", &b->edgeType);
        b->edgeType = CFDBoundary::TYPE_ID_WALL_SLIP;
        b->parCount = 0;
        b->par = nullptr;
    }

    if (strcmp(type, CFDBoundary::TYPE_PRESSURE) == 0) {
        b = new CFDBndPressure();
        //bNode->ToElement()->Attribute("edgeType", &b->edgeType);
        b->edgeType = CFDBoundary::TYPE_ID_PRESSURE;
        b->parCount = 2;
        b->par = new double[2];

        node = bNode->FirstChild("parameters");

        node1 = node->FirstChild("T");
        if (!node1) throw Exception("Parameter 'T' isn't specified  for BOUND_PRESSURE.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[1]);

        node1 = node->FirstChild("P");
        if (!node1) throw Exception("Parameter 'P' isn't specified  for BOUND_PRESSURE.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[0]);
    }

    if (!b) {
		throw Exception("Unknown boundary type '%s' specified.", Exception::TYPE_BOUND_UNKNOWN);
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
//	double u2 = pL.U2();
//	if (u2 < pL.cz*pL.cz) {
//        double u = sqrt(u2);
//
//        pR.p = par[0];
//        pR.r = pow( pR.p / (pL.p / pow( pL.r, pL.gam )), 1.0 / pL.gam );
//
//        double v = ( u + 2.0*pL.cz/(pL.gam - 1.0) - 2.0 / ( pL.gam - 1.0 ) * sqrt( pL.gam * pR.p / pR.r ) ) / u;
//
//        pR.u = v * pL.u;
//        pR.v = v * pL.v;
//
//	}
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

void CFDBndPressure::run(int iEdge, Param& pL, Param& pR)
{
    double p = par[0];
    double T = par[1];
    Edge &e = g->edges[iEdge];
    pR = pL;
    double Un = pR.u*e.n.x + pR.v*e.n.y;
    if ( Un < 0.0 ) {
        pR.T = T;
        pR.p = p - 0.5 * pR.r * Un * Un;
    } else {
        if ( pR.U2() < pR.cz*pR.cz ) {
            pR.p = p;
        }
    }
}

