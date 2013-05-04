#include "fvm_tvd.h"
#include "tinyxml.h"
#include <string>
#include "global.h"

void FVM_TVD::init(char * xmlFileName)
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


	ro		= new double[grid.cCount];
	ru		= new double[grid.cCount];
	rv		= new double[grid.cCount];
	re		= new double[grid.cCount];

	ro_old	= new double[grid.cCount];
	ru_old	= new double[grid.cCount];
	rv_old	= new double[grid.cCount];
	re_old	= new double[grid.cCount];

	ro_int	= new double[grid.cCount];
	ru_int	= new double[grid.cCount];
	rv_int	= new double[grid.cCount];
	re_int	= new double[grid.cCount];


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


void FVM_TVD::calcTimeStep()
{
	double tau = 1.0e+20;
	for (int iCell = 0; iCell < grid.cCount; iCell++)
	{
		Param p;
		convertConsToPar(iCell, p);
		double tmp = grid.cells[iCell].S/_max_(abs(p.u)+p.cz, abs(p.v)+p.cz);
		if (tmp < tau) tau = tmp;
	}
	tau  *= CFL;
	TAU = _min_(TAU, tau);
	printf("\n\nTime step TAU = %e.\n\n", TAU);
}

void FVM_TVD::calcGrad() 
{
	// ... сюда нужно вставить код Бариновой Марии
}

void FVM_TVD::run() 
{


	int nc = grid.cCount;
	int ne = grid.eCount;

	double			t		= 0.0;
	unsigned int	step	= 0;
	while (t < TMAX) 
	{
		t += TAU; 
		step++;
		memcpy(ro, ro_old, nc*sizeof(double));
		memcpy(ru, ru_old, nc*sizeof(double));
		memcpy(rv, rv_old, nc*sizeof(double));
		memcpy(re, re_old, nc*sizeof(double));

		// первый подшаг метода Р.-К.
		memset(ro_int, 0, nc*sizeof(double));
		memset(ru_int, 0, nc*sizeof(double));
		memset(rv_int, 0, nc*sizeof(double));
		memset(re_int, 0, nc*sizeof(double));
		calcGrad();
		for (int iEdge = 0; iEdge < ne; iEdge++)
		{
			double fr, fu, fv, fe;
			int c1	= grid.edges[iEdge].c1;
			int c2	= grid.edges[iEdge].c2;
			Vector n	= grid.edges[iEdge].n;
			double l	= grid.edges[iEdge].l*0.5;
			Param pL, pR;
			fr = 0.0;
			fu = 0.0;
			fv = 0.0;
			fe = 0.0;
			for (int iGP = 1; iGP < grid.edges[iEdge].cCount; iGP++) 
			{
				double fr1, fu1, fv1, fe1;
				reconstruct(iEdge, pL, pR, grid.edges[iEdge].c[iGP]);
				double __GAM = 1.4; // TODO: сделать правильное вычисление показателя адиабаты
				calcFlux(fr1, fu1, fv1, fe1, pL, pR, n, __GAM);
				fr += fr1;
				fu += fu1;
				fv += fv1;
				fe += fe1;

			}
			ro_int[c1] -= fr*l;
			ru_int[c1] -= fu*l;
			rv_int[c1] -= fv*l;
			re_int[c1] -= fe*l;
			if (c2 > -1) 
			{
				ro_int[c2] += fr*l;
				ru_int[c2] += fu*l;
				rv_int[c2] += fv*l;
				re_int[c2] += fe*l;
			}

		}
		for (int iCell = 0; iCell < nc; iCell++)
		{
			register double cfl = TAU/grid.cells[iCell].S;
			ro[iCell] += cfl*ro_int[iCell];
			ru[iCell] += cfl*ru_int[iCell];
			rv[iCell] += cfl*rv_int[iCell];
			re[iCell] += cfl*re_int[iCell];
		}

		// второй подшаг метода Р.-К.
		memset(ro_int, 0, nc*sizeof(double));
		memset(ru_int, 0, nc*sizeof(double));
		memset(rv_int, 0, nc*sizeof(double));
		memset(re_int, 0, nc*sizeof(double));
		calcGrad();
		for (int iEdge = 0; iEdge < ne; iEdge++)
		{
			double fr, fu, fv, fe;
			int c1	= grid.edges[iEdge].c1;
			int c2	= grid.edges[iEdge].c2;
			Vector n	= grid.edges[iEdge].n;
			double l	= grid.edges[iEdge].l*0.5;
			Param pL, pR;
			fr = 0.0;
			fu = 0.0;
			fv = 0.0;
			fe = 0.0;
			for (int iGP = 1; iGP < grid.edges[iEdge].cCount; iGP++) 
			{
				double fr1, fu1, fv1, fe1;
				reconstruct(iEdge, pL, pR, grid.edges[iEdge].c[iGP]);
				double __GAM = 1.4; // TODO: сделать правильное вычисление показателя адиабаты
				calcFlux(fr1, fu1, fv1, fe1, pL, pR, n, __GAM);
				fr += fr1;
				fu += fu1;
				fv += fv1;
				fe += fe1;

			}
			ro_int[c1] -= fr*l;
			ru_int[c1] -= fu*l;
			rv_int[c1] -= fv*l;
			re_int[c1] -= fe*l;
			if (c2 > -1) 
			{
				ro_int[c2] += fr*l;
				ru_int[c2] += fu*l;
				rv_int[c2] += fv*l;
				re_int[c2] += fe*l;
			}

		}
		for (int iCell = 0; iCell < nc; iCell++)
		{
			register double cfl = TAU/grid.cells[iCell].S;
			ro[iCell] += cfl*ro_int[iCell];
			ru[iCell] += cfl*ru_int[iCell];
			rv[iCell] += cfl*rv_int[iCell];
			re[iCell] += cfl*re_int[iCell];
		}

		// полусумма: формула (4.10) из icase-1997-65.pdf
		for (int iCell = 0; iCell < nc; iCell++)
		{
			ro[iCell] = 0.5*(ro_old[iCell]+ro[iCell]);
			ru[iCell] = 0.5*(ru_old[iCell]+ru[iCell]);
			rv[iCell] = 0.5*(rv_old[iCell]+rv[iCell]);
			re[iCell] = 0.5*(re_old[iCell]+re[iCell]);
		}




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

void FVM_TVD::save(int step)
{
	char fName[50];
	//sprintf(fName, "res_%010d.dat", step);
	//FILE * fp = fopen(fName, "w");
	////fprintf(fp, "# x y ro p u v e M\n");
	//for (int i = 0; i < grid.cCount; i++)
	//{
	//	Param p;
	//	convertConsToPar(i, p);
	//	fprintf(fp, "%25.16f %25.16f %25.16f %25.16f %25.16f %25.16f %25.16f %25.16f\n",	grid.cells[i].c.x,
	//																						grid.cells[i].c.y,
	//																						p.r,								//	плотность
	//																						p.p,								//	давление
	//																						p.u,								//	скорость
	//																						p.v,								//
	//																						p.e,								//	внутренняя энергия
	//																						sqrt(_sqr_(p.u)+_sqr_(p.v))/p.cz);	//	число Маха
	//}
	//fclose(fp);
	//// выводим градиенты примитивных переменных
	//printf("File '%s' saved...\n", fName);
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

	sprintf(fName, "res_%010d.vtk", step);
	FILE * fp = fopen(fName, "w");
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "GASDIN data file\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp, "POINTS %d float\n", grid.nCount);
	for (int i = 0; i < grid.nCount; i++)
	{
		fprintf(fp, "%f %f %f  ", grid.nodes[i].x,  grid.nodes[i].y, 0.0);
		if (i+1 % 8 == 0) fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	fprintf(fp, "CELLS %d %d\n", grid.cCount, 4*grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "3 %d %d %d\n", grid.cells[i].nodesInd[0], grid.cells[i].nodesInd[1], grid.cells[i].nodesInd[2]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "CELL_TYPES %d\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++) fprintf(fp, "5\n");
	fprintf(fp, "\n");

	fprintf(fp, "CELL_DATA %d\nSCALARS Density float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.r);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Pressure float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.p);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Temperature float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.T);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS MachNumber float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", (p.u*p.u+p.v*p.v)/p.cz);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "VECTORS Velosity float\n");
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f %25.16f %25.16f ", p.u, p.v, 0.0);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}


	fclose(fp);
	printf("File '%s' saved...\n", fName);


}

void FVM_TVD::calcFlux(double& fr, double& fu, double& fv, double& fe, Param pL, Param pR, Vector n, double GAM)
{
	{	// GODUNOV FLUX
		double RI, EI, PI, UI, VI, WI, UN, UT;
		//double unl = pL.u*n.x+pL.v*n.y;
		//double unr = pR.u*n.x+pR.v*n.y;
		//rim(RI, EI, PI, UN, UI, VI,  pL.r, pL.p, unl, pL.u, pL.v,  pR.r, pR.p, unr, pR.u, pR.v, n, GAM);
			
		double unl = pL.u*n.x+pL.v*n.y;
		double unr = pR.u*n.x+pR.v*n.y;
		double utl = pL.u*n.y-pL.v*n.x;
		double utr = pR.u*n.y-pR.v*n.x;
		rim_orig(RI, EI, PI, UN, UT, WI,  pL.r, pL.p, unl, utl, 0,  pR.r, pR.p, unr, utr, 0, GAM);
		
		UI = UN*n.x+UT*n.y;
		VI = UN*n.y-UT*n.x;

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


void FVM_TVD::reconstruct(int iEdge, Param& pL, Param& pR, Point p)
{
	// ... сюда нужно вставить код Бариновой Марии
	if (grid.edges[iEdge].type == Edge::TYPE_INNER) 
	{
		int c1	= grid.edges[iEdge].c1;
		int c2	= grid.edges[iEdge].c2;
		convertConsToPar(c1, pL);
		convertConsToPar(c2, pR);
	} else {
		int c1	= grid.edges[iEdge].c1;
		convertConsToPar(c1, pL);
		boundaryCond(iEdge, pL, pR);
	}
}


void FVM_TVD::boundaryCond(int iEdge, Param& pL, Param& pR)
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


void FVM_TVD::done()
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




Region & FVM_TVD::getRegionByCellType(int type)
{
	for (int i = 0; i < regCount; i++)
	{
		if (regions[i].cellType == type) return regions[i];
	}
	log("ERROR: unknown cell type %d...\n", type);
	EXIT(1);
}


Region   &	FVM_TVD::getRegion	(int iCell)
{
	return getRegionByCellType( grid.cells[iCell].type );
}

Material &	FVM_TVD::getMaterial	(int iCell)
{
	Region & reg = getRegion(iCell);
	return materials[reg.matId];
}


void FVM_TVD::convertParToCons(int iCell, Param & par)
{
	ro[iCell] = par.r;
	ru[iCell] = par.r*par.u;
	rv[iCell] = par.r*par.v;
	re[iCell] = par.r*(par.e+0.5*(par.u*par.u+par.v*par.v));
}

void FVM_TVD::convertConsToPar(int iCell, Param & par)
{
	par.r = ro[iCell];
	par.u = ru[iCell]/ro[iCell];
	par.v = rv[iCell]/ro[iCell];
	par.E = re[iCell]/ro[iCell];
	par.e = par.E-0.5*(par.u*par.u+par.v*par.v);
	Material& mat = getMaterial(iCell);
	mat.URS(par, 0);
}



