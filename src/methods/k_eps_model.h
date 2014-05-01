#ifndef _KEPSMODEL_H_
#define _KEPSMODEL_H_

#include "viscosity_model.h"

class KEpsModel :
	public ViscosityModel
{
public:
	KEpsModel(void);
	~KEpsModel(void);

	double getMuT(int iCell);
};

#endif