#include "fem_rkdg.h"
#include "tinyxml.h"
#include <string>
#include "global.h"

void FEM_RKDG::init(char * xmlFileName)
{
	TiXmlDocument doc( xmlFileName );
	bool loadOkay = doc.LoadFile( TIXML_ENCODING_UTF8 );
	if (!loadOkay)
	{
		log("ERROR: %s\n", doc.ErrorDesc());
		exit(doc.ErrorId());
	}
	
	TiXmlNode* task = 0;
	TiXmlElement* el = 0;
	TiXmlNode* node0 = 0;
	TiXmlNode* node1 = 0;
	task = doc.FirstChild( "task" );


	node0 = task->FirstChild("control");
	node0->FirstChild("TAU")->ToElement()->Attribute("value", &TAU);
	node0->FirstChild("TMAX")->ToElement()->Attribute("value", &TMAX);
	node0->FirstChild("CFL")->ToElement()->Attribute("value", &CFL);
	node0->FirstChild("FILE_OUTPUT_STEP")->ToElement()->Attribute("value", &FILE_SAVE_STEP);
	node0->FirstChild("LOG_OUTPUT_STEP")->ToElement()->Attribute("value", &PRINT_STEP);

	// чтение параметров о МАТЕРИАЛАХ
	node0 = task->FirstChild("materials");
	node0->ToElement()->Attribute("count", &matCount);;
	materials = new Material[matCount];
	TiXmlNode* matNode = node0->FirstChild("material"); 
	for (int i = 0; i < matCount; i++)
	{
		Material & mat = materials[i];
		matNode->ToElement()->Attribute("id", &mat.id);
		node1 = matNode->FirstChild("name");
		el = node1->ToElement();
		mat.name = el->GetText();
		node1 = matNode->FirstChild("parameters");
		node1->FirstChild( "M"  )->ToElement()->Attribute( "value", &mat.M  );
		node1->FirstChild( "Cp" )->ToElement()->Attribute( "value", &mat.Cp );
		node1->FirstChild( "K"  )->ToElement()->Attribute( "value", &mat.K  );
		node1->FirstChild( "ML" )->ToElement()->Attribute( "value", &mat.ML );
		matNode = matNode->NextSibling("material");
	}

	// чтение параметров о РЕГИОНАХ
	node0 = task->FirstChild("regions");
	node0->ToElement()->Attribute("count", &regCount);
	regions = new Region[regCount];
	TiXmlNode* regNode = node0->FirstChild("region");
	for (int i = 0; i < regCount; i++)
	{
		Region & reg = regions[i];
		regNode->ToElement()->Attribute("id", &reg.id);
		regNode->FirstChild("material")->ToElement()->Attribute("id", &reg.matId);
		regNode->FirstChild("cell")->ToElement()->Attribute("type", &reg.cellType);
		
		node1 = regNode->FirstChild("parameters");
		node1->FirstChild( "Vx" )->ToElement()->Attribute( "value", &reg.par.u );
		node1->FirstChild( "Vy" )->ToElement()->Attribute( "value", &reg.par.v );
		node1->FirstChild( "T"  )->ToElement()->Attribute( "value", &reg.par.T );
		node1->FirstChild( "P"  )->ToElement()->Attribute( "value", &reg.par.p );
		
		Material& mat = materials[reg.matId];
		mat.URS(reg.par, 2);	// r=r(p,T)
		mat.URS(reg.par, 1);	// e=e(p,r)
		
		regNode = regNode->NextSibling("region");
	}

	// чтение параметров о ГРАНИЧНЫХ УСЛОВИЯХ
	node0 = task->FirstChild("boundaries");
	node0->ToElement()->Attribute("count", &bCount);
	boundaries = new Boundary[bCount];
	TiXmlNode* bNode = node0->FirstChild("boundCond");
	for (int i = 0; i < bCount; i++)
	{
		Boundary & b = boundaries[i];
		bNode->ToElement()->Attribute("edgeType", &b.edgeType);
		const char * str = bNode->FirstChild("type")->ToElement()->GetText();
		if (strcmp(str, "BOUND_WALL") == 0) 
		{
			b.parCount = 0;
			b.par = NULL;
			b.type = Boundary::BOUND_WALL;
		} else
		if (strcmp(str, "BOUND_OUTLET") == 0) 
		{
			b.parCount = 0;
			b.par = NULL;
			b.type = Boundary::BOUND_OUTLET;
		} else
		if (strcmp(str, "BOUND_INLET") == 0) 
		{
			b.parCount = 4;
			b.par = new double[4];
			b.type = Boundary::BOUND_INLET;

			node1 = bNode->FirstChild("parameters");
			node1->FirstChild( "T"  )->ToElement()->Attribute( "value", &b.par[0] );
			node1->FirstChild( "P"  )->ToElement()->Attribute( "value", &b.par[1] );
			node1->FirstChild( "Vx" )->ToElement()->Attribute( "value", &b.par[2] );
			node1->FirstChild( "Vy" )->ToElement()->Attribute( "value", &b.par[3] );
		} else {
			log("ERROR: unsupported boundary condition type '%s'", str);
			EXIT(1);
		}

		
		
		
		bNode = bNode->NextSibling("boundCond");
	}


	node0 = task->FirstChild("mesh");
	const char* fName = task->FirstChild("mesh")->FirstChild("name")->ToElement()->Attribute("value");
	grid.initFromFiles((char*)fName);


	ro		= new VECTOR[grid.cCount];
	ru		= new VECTOR[grid.cCount];
	rv		= new VECTOR[grid.cCount];
	re		= new VECTOR[grid.cCount];

	ro_old	= new VECTOR[grid.cCount];
	ru_old	= new VECTOR[grid.cCount];
	rv_old	= new VECTOR[grid.cCount];
	re_old	= new VECTOR[grid.cCount];

	ro_int	= new VECTOR[grid.cCount];
	ru_int	= new VECTOR[grid.cCount];
	rv_int	= new VECTOR[grid.cCount];
	re_int	= new VECTOR[grid.cCount];

	for (int i = 0; i < grid.cCount; i++)
	{
		ro[i].init(FUNC_COUNT);
		ru[i].init(FUNC_COUNT);
		rv[i].init(FUNC_COUNT);
		re[i].init(FUNC_COUNT);

		ro_old[i].init(FUNC_COUNT);
		ru_old[i].init(FUNC_COUNT);
		rv_old[i].init(FUNC_COUNT);
		re_old[i].init(FUNC_COUNT);

		ro_int[i].init(FUNC_COUNT);
		ru_int[i].init(FUNC_COUNT);
		rv_int[i].init(FUNC_COUNT);
		re_int[i].init(FUNC_COUNT);
	}

	for (int i = 0; i < grid.cCount; i++)
	{
		Region & reg = getRegion(i);
		convertParToCons(i, reg.par);
	}

	memcpy(ro_old, ro, grid.cCount*sizeof(double));
	memcpy(ru_old, ru, grid.cCount*sizeof(double));
	memcpy(rv_old, rv, grid.cCount*sizeof(double));
	memcpy(re_old, re, grid.cCount*sizeof(double));

	calcTimeStep();
	save(0);
}


void FEM_RKDG::calcTimeStep()
{
	double tau = 1.0e+20;
	for (int iCell = 0; iCell < grid.cCount; iCell++)
	{
		Param p;
		convertConsToPar(iCell, grid.cells[iCell].c, p);
		double tmp = grid.cells[iCell].S/_max_(abs(p.u)+p.cz, abs(p.v)+p.cz);
		if (tmp < tau) tau = tmp;
	}
	tau  *= CFL;
	TAU = _min_(TAU, tau);
	printf("\n\nTime step TAU = %e.\n\n", TAU);
}

void FEM_RKDG::calcGP() {}   //TODO:

void FEM_RKDG::calcMatr() {} //TODO:


void FEM_RKDG::run() 
{


	int nc = grid.cCount;
	int ne = grid.eCount;

	double			t		= 0.0;
	unsigned int	step	= 0;
	while (t < TMAX) 
	{
		t += TAU; 
		step++;
		
		calcLimiters();
		calcConvectionVol();
		calcConvectionSurf();
		calcDiffusionVol();
		calcDiffusionSurf();
		calcNewFields();

		
		if (step % FILE_SAVE_STEP == 0)
		{
			save(step);
		}
		if (step % PRINT_STEP == 0)
		{
			log("step: %d\t\ttime step: %.16f\n", step, t);
		}
	}

}

void FEM_RKDG::save(int step)
{
	char fName[50];
	sprintf(fName, "res_%010d.dat", step);
	FILE * fp = fopen(fName, "w");
	//fprintf(fp, "# x y ro p u v e M\n");
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, grid.cells[i].c, p);
		fprintf(fp, "%25.16f %25.16f %25.16f %25.16f %25.16f %25.16f %25.16f %25.16f\n",	grid.cells[i].c.x,
																							grid.cells[i].c.y,
																							p.r,								//	плотность
																							p.p,								//	давление
																							p.u,								//	скорость
																							p.v,								//
																							p.e,								//	внутренняя энергия
																							sqrt(_sqr_(p.u)+_sqr_(p.v))/p.cz);	//	число Маха
	}
	fclose(fp);
	// выводим градиенты примитивных переменных
	printf("File '%s' saved...\n", fName);
	//sprintf(fName, "grad_%010d.dat", step);
	//fp = fopen(fName, "w");
	////fprintf(fp, "# x y ro p u v e M\n");
	//for (int i = 0; i < grid.cCount; i++)
	//{
	//	fprintf(fp, "%25.16f\t%25.16f\t\t%25.16f\t%25.16f\t\t%25.16f\t%25.16f\t\t%25.16f\t%25.16f\n",	primGrad[ID_R][i].x,
	//																									primGrad[ID_R][i].y,
	//																									primGrad[ID_P][i].x,
	//																									primGrad[ID_P][i].y,
	//																									primGrad[ID_U][i].x,
	//																									primGrad[ID_U][i].y,
	//																									primGrad[ID_V][i].x,
	//																									primGrad[ID_V][i].y	);
	//}
	//fclose(fp);
	//printf("File '%s' saved...\n", fName);
}



void FEM_RKDG::calcFlux(double& fr, double& fu, double& fv, double& fe, Param pL, Param pR, Vector n, double GAM)
{
	{	// GODUNOV FLUX
		double RI, EI, PI, UI, VI, UN;
		double unl = pL.u*n.x+pL.v*n.y;
		double unr = pR.u*n.x+pR.v*n.y;
		rim(RI, EI, PI, UN, UI, VI,  pL.r, pL.p, unl, pL.u, pL.v,  pR.r, pR.p, unr, pR.u, pR.v, n, GAM);
			
		fr = RI*UN;
		fu = fr*UI+PI*n.x;
		fv = fr*VI+PI*n.y;
		fe = (RI*(EI+0.5*(UI*UI+VI*VI))+PI)*UN;
	}
	//{	// LAX-FRIEDRIX FLUX
	//	double unl = pL.u*n.x+pL.v*n.y;
	//	double unr = pR.u*n.x+pR.v*n.y;
	//	double rol, rul, rvl, rel,  ror, rur, rvr, rer;
	//	double alpha = _max_(fabs(unl)+sqrt(GAM*pL.p/pL.r), fabs(unr)+sqrt(GAM*pR.p/pR.r));
	//	pL.getToCons(rol, rul, rvl, rel);
	//	pR.getToCons(ror, rur, rvr, rer);
	//	double frl = rol*unl;
	//	double frr = ror*unr;
	//	fr = 0.5*(frr+frl								- alpha*(ror-rol));
	//	fu = 0.5*(frr*pR.u+frl*pL.u + (pR.p+pL.p)*n.x	- alpha*(rur-rul));
	//	fv = 0.5*(frr*pR.v+frl*pL.v + (pR.p+pL.p)*n.y	- alpha*(rvr-rvl));
	//	fe = 0.5*((rer+pR.p)*unr + (rel+pL.p)*unl		- alpha*(rer-rel));
	//}
}


void FEM_RKDG::calcConvectionSurf()
{
	double fr, fu, fv, fe;
	for (int i = 0; i < grid.eCount; i++)
	{
		Edge& edge = grid.edges[i];
		Vector n = edge.n;
		if (edge.type == Edge::TYPE_INNER)
		{
			int c1 = edge.c1;
			int c2 = edge.c2;
			VECTOR intRO1(FUNC_COUNT), intRU1(FUNC_COUNT), intRV1(FUNC_COUNT), intRE1(FUNC_COUNT);
			VECTOR intRO2(FUNC_COUNT), intRU2(FUNC_COUNT), intRV2(FUNC_COUNT), intRE2(FUNC_COUNT);
			for (int iGP = 0; iGP < CELL_GP_COUNT; iGP++)
			{
				Point& pt = edgeGP[i][iGP];
				Param pL, pR;
				convertConsToPar(c1, pt, pL);
				convertConsToPar(c2, pt, pR);
				double __GAM = 1.4; // TODO: сделать правильное вычисление показателя адиабаты
				calcFlux(fr, fu, fv, fe, pL, pR, n, __GAM);
				
				VECTOR tmp(FUNC_COUNT), tmp1(FUNC_COUNT), tmp2(FUNC_COUNT);
				tmp1     = getF(c1, pt);
				tmp1    *= edgeGW[iGP];
				tmp2     = getF(c2, pt);
				tmp2    *= edgeGW[iGP];
				
				tmp     = tmp1;
				tmp    *= fr;
				intRO1 += tmp;
				tmp     = tmp2;
				tmp    *= fr;
				intRO2 += tmp;

				tmp     = tmp1;
				tmp    *= fu;
				intRU1 += tmp;
				tmp     = tmp2;
				tmp    *= fu;
				intRU2 += tmp;

				tmp     = tmp1;
				tmp    *= fv;
				intRV1 += tmp;
				tmp     = tmp2;
				tmp    *= fv;
				intRV2 += tmp;

				tmp     = tmp1;
				tmp    *= fe;
				intRE1 += tmp;
				tmp     = tmp2;
				tmp    *= fe;
				intRE2 += tmp;
			}
			intRO1 *= edgeJ[i];
			intRO2 *= edgeJ[i];
			intRU1 *= edgeJ[i];
			intRU2 *= edgeJ[i];
			intRV1 *= edgeJ[i];
			intRV2 *= edgeJ[i];
			intRE1 *= edgeJ[i];
			intRE2 *= edgeJ[i];

			ro_int[c1] -= intRO1;
			ru_int[c1] -= intRU1;
			rv_int[c1] -= intRV1;
			re_int[c1] -= intRE1;

			ro_int[c2] += intRO1;
			ru_int[c2] += intRU1;
			rv_int[c2] += intRV1;
			re_int[c2] += intRE1;
		} else {
			int c1 = edge.c1;
			VECTOR intRO1(FUNC_COUNT), intRU1(FUNC_COUNT), intRV1(FUNC_COUNT), intRE1(FUNC_COUNT);
			for (int iGP = 0; iGP < CELL_GP_COUNT; iGP++)
			{
				Point& pt = edgeGP[i][iGP];
				Param pL, pR;
				convertConsToPar(c1, pt, pL);
				boundaryCond(i, pt, pL, pR);
				double __GAM = 1.4; // TODO: сделать правильное вычисление показателя адиабаты
				calcFlux(fr, fu, fv, fe, pL, pR, n, __GAM);
				
				VECTOR tmp(FUNC_COUNT), tmp1(FUNC_COUNT), tmp2(FUNC_COUNT);
				tmp1     = getF(c1, pt);
				tmp1    *= edgeGW[iGP];
				
				tmp     = tmp1;
				tmp    *= fr;
				intRO1 += tmp;

				tmp     = tmp1;
				tmp    *= fu;
				intRU1 += tmp;

				tmp     = tmp1;
				tmp    *= fv;
				intRV1 += tmp;

				tmp     = tmp1;
				tmp    *= fe;
				intRE1 += tmp;
			}
			intRO1 *= edgeJ[i];
			intRU1 *= edgeJ[i];
			intRV1 *= edgeJ[i];
			intRE1 *= edgeJ[i];

			ro_int[c1] -= intRO1;
			ru_int[c1] -= intRU1;
			rv_int[c1] -= intRV1;
			re_int[c1] -= intRE1;
		}
	}
}


void FEM_RKDG::calcConvectionVol()
{
	for (int i = 0; i < grid.cCount; i++)
	{
		VECTOR intRO(FUNC_COUNT), intRU(FUNC_COUNT), intRV(FUNC_COUNT), intRE(FUNC_COUNT);
		for (int iGP = 0; iGP < CELL_GP_COUNT; iGP++)
		{
			Point& pt = cellGP[i][iGP];
			Param par;
			convertConsToPar(i, pt, par);
			
			double fr = par.r*par.u;
			double fu = fr*par.u+par.p;
			double fv = fr*par.v;
			double fe = par.u*(par.r*par.E+par.p);
			
			double gr = par.r*par.v;
			double gu = gr*par.u;
			double gv = gr*par.v+par.p;
			double ge = par.v*(par.r*par.E+par.p);
			
			VECTOR tmpInt(FUNC_COUNT), tmpInt1(FUNC_COUNT);
			
			
			// RO
			tmpInt = 0.0;
			tmpInt1  = getDFDX(i, pt); 
			tmpInt1 *= fr;
			tmpInt += tmpInt1;
			
			tmpInt1  = getDFDY(i, pt); 
			tmpInt1 *= gr;
			tmpInt += tmpInt1;
			
			tmpInt *= cellGW[iGP];
			intRO  += tmpInt;

			// RU
			tmpInt = 0.0;
			tmpInt1  = getDFDX(i, pt); 
			tmpInt1 *= fu;
			tmpInt += tmpInt1;

			tmpInt1  = getDFDY(i, pt); 
			tmpInt1 *= gu;
			tmpInt += tmpInt1;

			tmpInt *= cellGW[iGP];
			intRU  += tmpInt;

			// RV
			tmpInt = 0.0;
			tmpInt1  = getDFDX(i, pt); 
			tmpInt1 *= fv;
			tmpInt += tmpInt1;

			tmpInt1  = getDFDY(i, pt); 
			tmpInt1 *= gv;
			tmpInt += tmpInt1;

			tmpInt *= cellGW[iGP];
			intRV  += tmpInt;

			// RE
			tmpInt = 0.0;
			tmpInt1  = getDFDX(i, pt); 
			tmpInt1 *= fe;
			tmpInt += tmpInt1;

			tmpInt1  = getDFDY(i, pt); 
			tmpInt1 *= ge;
			tmpInt += tmpInt1;

			tmpInt *= cellGW[iGP];
			intRE  += tmpInt;
		}
		intRO *= cellJ[i];
		intRU *= cellJ[i];
		intRV *= cellJ[i];
		intRE *= cellJ[i];

		ro_int[i] += intRO;
		ru_int[i] += intRU;
		rv_int[i] += intRV;
		re_int[i] += intRE;
	}
}

void FEM_RKDG::calcDiffusionVol() {} //TODO:

void FEM_RKDG::calcDiffusionSurf() {} //TODO:

void FEM_RKDG::calcNewFields() 
{
	VECTOR tmp(FUNC_COUNT);

	for (int i = 0; i < grid.cCount; i++)
	{
		tmp = ro_int[i];
		tmp *= TAU;
		ro[i] += tmp;

		tmp = ru_int[i];
		tmp *= TAU;
		ru[i] += tmp;

		tmp = rv_int[i];
		tmp *= TAU;
		rv[i] += tmp;

		tmp = re_int[i];
		tmp *= TAU;
		re[i] += tmp;
	}
}

void FEM_RKDG::calcLimiters() {}

void FEM_RKDG::reconstruct(int iEdge, Point pt, Param& pL, Param& pR)
{
	if (grid.edges[iEdge].type == Edge::TYPE_INNER) 
	{
		int c1	= grid.edges[iEdge].c1;
		int c2	= grid.edges[iEdge].c2;
		convertConsToPar(c1, pt, pL);
		convertConsToPar(c2, pt, pR);
	} else {
		int c1	= grid.edges[iEdge].c1;
		convertConsToPar(c1, pt, pL);
		boundaryCond(iEdge, pt, pL, pR);
	}
}


void FEM_RKDG::boundaryCond(int iEdge, Point pt, Param& pL, Param& pR)
{
	int iBound = -1;
	for (int i = 0; i < bCount; i++)
	{
		if (grid.edges[iEdge].type == boundaries[i].edgeType) 
		{
			iBound = i;
			break;
		}
	}
	if (iBound < 0)
	{
		log("ERROR (boundary condition): unknown edge type of edge %d...\n", iEdge);
		EXIT(1);
	}
	Boundary& b = boundaries[iBound];
	int c1	= grid.edges[iEdge].c1;
	Material& m = getMaterial(c1);
	switch (b.type)
	{
	case Boundary::BOUND_INLET:
		pR.T  = b.par[0];		//!< температура
		pR.p  = b.par[1];		//!< давление
		pR.u  = b.par[2];		//!< первая компонента вектора скорости
		pR.v  = b.par[3];		//!< вторая компонента вектора скорости
		
		m.URS(pR, 2);
		m.URS(pR, 1);
		break;
	
	case Boundary::BOUND_OUTLET:
		pR = pL;
		break;
	
	case Boundary::BOUND_WALL:
		pR = pL;
		double Un = pL.u*grid.edges[iEdge].n.x+pL.v*grid.edges[iEdge].n.y;
		Vector V;
		V.x = grid.edges[iEdge].n.x*Un*2.0; 
		V.y = grid.edges[iEdge].n.y*Un*2.0;
		pR.u = pL.u-V.x;
		pR.v = pL.v-V.y;
		break;
	}
}


void FEM_RKDG::done()
{
	delete[] ro;
	delete[] ru;
	delete[] rv;
	delete[] re;

	delete[] ro_old;
	delete[] ru_old;
	delete[] rv_old;
	delete[] re_old;

	delete[] ro_int;
	delete[] ru_int;
	delete[] rv_int;
	delete[] re_int;
}




Region & FEM_RKDG::getRegionByCellType(int type)
{
	for (int i = 0; i < regCount; i++)
	{
		if (regions[i].cellType == type) return regions[i];
	}
	log("ERROR: unknown cell type %d...\n", type);
	EXIT(1);
}


Region   &	FEM_RKDG::getRegion	(int iCell)
{
	return getRegionByCellType( grid.cells[iCell].type );
}

Material &	FEM_RKDG::getMaterial	(int iCell)
{
	Region & reg = getRegion(iCell);
	return materials[reg.matId];
}


void FEM_RKDG::convertParToCons(int iCell, Param & par)
{
	ro[iCell][0] = par.r;
	ru[iCell][0] = par.r*par.u;
	rv[iCell][0] = par.r*par.v;
	re[iCell][0] = par.r*(par.e+0.5*(par.u*par.u+par.v*par.v));

	for (int i = 1; i < FUNC_COUNT; i++)
	{
		ro[iCell][i] = 0.0;
		ru[iCell][i] = 0.0;
		rv[iCell][i] = 0.0;
		re[iCell][i] = 0.0;
	}
}

void FEM_RKDG::convertConsToPar(int iCell, Point pt, Param & par)
{
	par.r = getRO(iCell, pt);
	par.u = getRU(iCell, pt)/par.r;
	par.v = getRV(iCell, pt)/par.r;
	par.E = getRE(iCell, pt)/par.r;
	par.e = par.E-0.5*(par.u*par.u+par.v*par.v);
	Material& mat = getMaterial(iCell);
	mat.URS(par, 0);
}






VECTOR FEM_RKDG::getF(int iCell, Point pt)
{
	if (FUNC_COUNT == 3 )
	{
		VECTOR v(3);
		v[0] = 1.0;
		v[1] = (pt.x - grid.cells[iCell].c.x)/grid.cells[iCell].HX;
		v[2] = (pt.y - grid.cells[iCell].c.y)/grid.cells[iCell].HY;
		return v;
	}
	else if (FUNC_COUNT == 6)
	{
		VECTOR v(6);
		v[0] = 1.0;
		v[1] = (pt.x - grid.cells[iCell].c.x)/grid.cells[iCell].HX;
		v[2] = (pt.y - grid.cells[iCell].c.y)/grid.cells[iCell].HY;
		v[3] = (pt.x - grid.cells[iCell].c.x)*(pt.x - grid.cells[iCell].c.x)/(grid.cells[iCell].HX*grid.cells[iCell].HX);
		v[4] = (pt.y - grid.cells[iCell].c.y)*(pt.y - grid.cells[iCell].c.y)/(grid.cells[iCell].HY*grid.cells[iCell].HY);
		v[5] = (pt.x - grid.cells[iCell].c.x)*(pt.y - grid.cells[iCell].c.y)/(grid.cells[iCell].HX*grid.cells[iCell].HY);
		return v;
	} else {
		log("ERROR: wrong basic functions count!\n");
		EXIT(1);
	}
}

VECTOR FEM_RKDG::getDFDX(int iCell, Point pt)
{
	if (FUNC_COUNT == 3 )
	{
		VECTOR v(3);
		v[0] = 0.0;
		v[1] = 1.0/grid.cells[iCell].HX;
		v[2] = 0.0;
		return v;
	}
	else if (FUNC_COUNT == 6)
	{
		VECTOR v(6);
		v[0] = 0.0;
		v[1] = 1.0/grid.cells[iCell].HX;
		v[2] = 0.0;
		v[3] = 2.0*(pt.x - grid.cells[iCell].c.x)/(grid.cells[iCell].HX*grid.cells[iCell].HX);
		v[4] = 0.0;
		v[5] = (pt.y - grid.cells[iCell].c.y)/(grid.cells[iCell].HX*grid.cells[iCell].HY);
		return v;
	} else {
		log("ERROR: wrong basic functions count!\n");
		EXIT(1);
	}
}

VECTOR FEM_RKDG::getDFDY(int iCell, Point pt)
{
	if (FUNC_COUNT == 3 )
	{
		VECTOR v(3);
		v[0] = 0.0;
		v[1] = 0.0;
		v[2] = 1.0/grid.cells[iCell].HY;
		return v;
	}
	else if (FUNC_COUNT == 6)
	{
		VECTOR v(6);
		v[0] = 0.0;
		v[1] = 0.0;
		v[2] = 1.0/grid.cells[iCell].HY;
		v[3] = 0.0;
		v[4] = 2.0*(pt.y - grid.cells[iCell].c.y)/(grid.cells[iCell].HY*grid.cells[iCell].HY);
		v[5] = (pt.x - grid.cells[iCell].c.x)/(grid.cells[iCell].HX*grid.cells[iCell].HY);
		return v;
	} else {
		log("ERROR: wrong basic functions count!\n");
		EXIT(1);
	}
}
