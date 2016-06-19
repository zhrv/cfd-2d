//
// Created by appmath on 19.06.16.
//

#ifndef CFD_2D_BND_COND_HEAT_H
#define CFD_2D_BND_COND_HEAT_H

#include "bnd_cond.h"

struct HeatBoundary : public CFDBoundary
{
    static CFDBoundary* create(TiXmlNode* bNode, Grid * g);
    virtual void run(int iEdgs, Param& pL, Param& pR) = 0;

    HeatBoundary() : CFDBoundary() {}
    ~HeatBoundary() { if (par) delete[] par; par = NULL; }


    static const char* TYPE_CONST;
    static const char* TYPE_FLUX_CONST;
    static const char* TYPE_FLUX_ZERO;

};

//typedef std::vector< CFDBoundary* > CFDBoundaries;

struct HeatBndConst : public HeatBoundary
{
    virtual void run(int iEdgs, Param& pL, Param& pR);
};

struct HeatBndFluxConst : public HeatBoundary
{
    virtual void run(int iEdgs, Param& pL, Param& pR);
};

struct HeatBndFluxZero : public HeatBoundary
{
    virtual void run(int iEdgs, Param& pL, Param& pR);
};


#endif //CFD_2D_BND_COND_HEAT_H
