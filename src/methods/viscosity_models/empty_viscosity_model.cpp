#include "empty_viscosity_model.h"

EmptyViscosityModel::EmptyViscosityModel(void)
{
}

EmptyViscosityModel::~EmptyViscosityModel(void)
{
}

void EmptyViscosityModel::init( Grid * grid, double * ro, double *ru, double * rv, double * ro_m, double * u_m, double * v_m, Vector * gradU, Vector * gradV, double * Txx, double * Tyy, double * Txy, const double mu )
{
}

double EmptyViscosityModel::getMuT( const int iCell )
{
	return 0.0;
}

void EmptyViscosityModel::calcMuT( double * cTau )
{
}

void EmptyViscosityModel::done()
{
}
