#pragma once

#include "method.h"

class FEM_DG : public Method
{
	friend class LimiterDG;
protected:
	const static int BASE_FUNC_COUNT = 3;

	double **ro;
	double **ru;
	double **rv;
	double **re;

public:
	virtual void getFields(double &fRO, double &fRU, double &fRV, double &fRE, int iCell, Point p) = 0;
	virtual void getFields(double &fRO, double &fRU, double &fRV, double &fRE, int iCell, double x, double y) = 0;

	virtual double getField(int fld, int iCell, Point p) = 0;
	virtual double getField(int fld, int iCell, double x, double y) = 0;

	virtual double getF(int id, int iCell, Point p) = 0;
	virtual double getF(int id, int iCell, double x, double y) = 0;


};