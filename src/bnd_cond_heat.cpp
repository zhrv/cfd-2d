//
// Created by appmath on 19.06.16.
//

#include <mesh/grid.h>
#include "bnd_cond_heat.h"


const char* HeatBoundary::TYPE_CONST        = "BOUND_HEAT_CONST";
const char* HeatBoundary::TYPE_FLUX_CONST   = "BOUND_HEAT_FLUX_CONST";
const char* HeatBoundary::TYPE_FLUX_ZERO    = "BOUND_HEAT_FLUX_ZERO";

CFDBoundary* HeatBoundary::create(TiXmlNode* bNode, Grid * g)
{
    CFDBoundary * b = 0;
    TiXmlNode* node = 0;
    TiXmlNode* node1 = 0;

    const char * type = bNode->FirstChild("type")->ToElement()->GetText();

    if (strcmp(type, HeatBoundary::TYPE_CONST) == 0) {
        b = new HeatBndConst();
        bNode->ToElement()->Attribute("edgeType", &b->edgeType);
        b->parCount = 1;
        b->par = new double[1];
        node = bNode->FirstChild("parameters");

        node1 = node->FirstChild("T");
        if (!node1) throw Exception((char*)"Parameter 'T' isn't specified  for BOUND_HEAT_CONST.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[0]);
    }

    if (strcmp(type, HeatBoundary::TYPE_FLUX_ZERO) == 0) {
        b = new HeatBndFluxZero();
        bNode->ToElement()->Attribute("edgeType", &b->edgeType);
        b->parCount = 0;
        b->par = NULL;
    }

    if (strcmp(type, HeatBoundary::TYPE_FLUX_CONST) == 0) {
        b = new HeatBndFluxConst();
        bNode->ToElement()->Attribute("edgeType", &b->edgeType);
        b->parCount = 1;
        b->par = new double[1];
        node = bNode->FirstChild("parameters");

        node1 = node->FirstChild("flux");
        if (!node1) throw Exception((char*)"Parameter 'flux' isn't specified  for BOUND_HEAT_FLUX_CONST.", Exception::TYPE_BOUND_NOPAR);
        node1->ToElement()->Attribute("value", &b->par[0]);
    }

    if (!b) {
        char msg[256];
        sprintf(msg, "Unknown boundary type '%s' specified.", type);
        throw Exception((char*)msg, Exception::TYPE_BOUND_UNKNOWN);
    }

    const char * name = bNode->FirstChild("name")->ToElement()->GetText();
    strcpy(b->name, name);
    b->g = g;
    return b;
}


void HeatBndConst::run(int iEdge, Param& pL, Param& pR)
{
    Edge &e = g->edges[iEdge];
    Cell &c = g->cells[e.c1];
    Vector &n = e.n;
    Vector ce = e.c[0];
    ce -= c.c;
    double h = 2.0*_LENGTH_(ce);
    pR.T = par[0];
    // TODO calc grad(T)
    pR.Qt[0] = pR.Qt[1] = (pR.T-pL.T)/(h*(n.x+n.y));

}

void HeatBndFluxConst::run(int iEdge, Param& pL, Param& pR)
{
    // TODO
    pR.T = pL.T;
    pR.Qt[0] = 0.0;
    pR.Qt[1] = 0.0;
}

void HeatBndFluxZero::run(int iEdge, Param& pL, Param& pR)
{
    pR.T = pL.T;
    pR.Qt[0] = 0.0;
    pR.Qt[1] = 0.0;
}

