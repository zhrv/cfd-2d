#include "fvm_tvd_implicit.h"
#include "tinyxml.h"
#include <string>
#include "global.h"

void FVM_TVD_IMPLICIT::init(char * xmlFileName)
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
/*
	ro_old	= new double[grid.cCount];
	ru_old	= new double[grid.cCount];
	rv_old	= new double[grid.cCount];
	re_old	= new double[grid.cCount];

	ro_int	= new double[grid.cCount];
	ru_int	= new double[grid.cCount];
	rv_int	= new double[grid.cCount];
	re_int	= new double[grid.cCount];
*/

	for (int i = 0; i < grid.cCount; i++)
	{
		Region & reg = getRegion(i);
		convertParToCons(i, reg.par);
	}
/*
	memcpy(ro_old, ro, grid.cCount*sizeof(double));
	memcpy(ru_old, ru, grid.cCount*sizeof(double));
	memcpy(rv_old, rv, grid.cCount*sizeof(double));
	memcpy(re_old, re, grid.cCount*sizeof(double));
*/
	//заполнение массива cellsEdges значениями.
	cellsEdges = new int* [grid.cCount];
	for (int iCells = 0; iCells < grid.cCount; ++iCells)
	{
		cellsEdges[iCells] = new int [3];
		for (int j = 0; j < 3; ++j)
			cellsEdges[iCells][j] = -1;
	}
	for (int iEdge = 0; iEdge < grid.eCount; ++iEdge)
	{
		assert(grid.edges[iEdge].c1 < grid.cCount);
		if (grid.edges[iEdge].c1 >= 0)
		{	
			bool edgeExist = false;
			for (int j = 0; j < 3; ++j)
			{
				if (cellsEdges[grid.edges[iEdge].c1][j] == iEdge) 
				{
						edgeExist = true;
						break;
				}
			}
			if (!edgeExist)
			{
				for (int j = 0; j < 3; ++j)	
				{
					if (cellsEdges[grid.edges[iEdge].c1][j] == -1)
					{
						cellsEdges[grid.edges[iEdge].c1][j] = iEdge;
						break;
					}
				}
			}
		}
		assert(grid.edges[iEdge].c2 < grid.cCount);
		if (grid.edges[iEdge].c2 >= 0)
		{	
			bool edgeExist = false;
			for (int j = 0; j < 3; ++j)
			{
				if (cellsEdges[grid.edges[iEdge].c2][j] == iEdge) 
				{
						edgeExist = true;
						break;
				}
			}
			if (!edgeExist)
			{
				for (int j = 0; j < 3; ++j)	
				{
					if (cellsEdges[grid.edges[iEdge].c2][j] == -1)
					{
						cellsEdges[grid.edges[iEdge].c2][j] = iEdge;
						break;
					}
				}
			}
		}
	}

	calcTimeStep();
	save(0);
}

double **FVM_TVD_IMPLICIT::allocMtx4()
{
	double		**tempMtx4 = new double* [4];
	for (int i = 0; i < 4; ++i) tempMtx4[i] = new double [4];
	return tempMtx4;
}
void   FVM_TVD_IMPLICIT::freeMtx4(double **mtx4)
{
	for (int i = 0; i < 4; ++i)
		delete [] mtx4[i];
	delete [] mtx4;
}

void FVM_TVD_IMPLICIT::calcTimeStep()
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

void FVM_TVD_IMPLICIT::eigenValues(double** dst4, double c, double u, double nx, double v, double ny)
{
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			dst4[i][j] = 0.0;
	double qn = u*nx + v*ny;
	dst4[0][0] = qn - c;
	dst4[1][1] = qn;
	dst4[2][2] = qn + c;
	dst4[3][3] = qn;
}

void FVM_TVD_IMPLICIT::rightEigenVector(double **dst4, double c, double u, double nx, double v, double ny, double H)
{
	double qn = u*nx + v*ny;
	double q2 = u*u + v*v;
	double lx = -ny; 
	double ly = nx;
	double ql = u*lx + v*ly;
	dst4[0][0] = 1;			dst4[0][1] = 1;		dst4[0][2] = 1;			dst4[0][3] = 0;
	dst4[1][0] = u-c*nx;	dst4[1][1] = u;		dst4[1][2] = u+c*nx;	dst4[1][3] = lx;
	dst4[2][0] = v-c*ny;	dst4[2][1] = v;		dst4[2][2] = v+c*ny;	dst4[2][3] = ly;
	dst4[3][0] = H-qn*c;	dst4[3][1] = q2/2;	dst4[3][2] = H+qn*c;	dst4[3][3] = ql;
}

void FVM_TVD_IMPLICIT::leftEigenVector(double **dst4, double c, double GAM, double u, double nx, double v, double ny)
{
	double qn = u*nx + v*ny;
	double q2 = u*u + v*v;
	double lx = -ny; 
	double ly = nx;
	double ql = u*lx + v*ly;
	double g1 = GAM-1;
	double c2 = c*c;
	dst4[0][0] = 0.5*(0.5*q2*g1/c2 + qn/c);	dst4[0][1] = -0.5*(g1*u/c2 + nx/c);	dst4[0][2] = -0.5*(g1*v/c2 + ny/c);	dst4[0][3] = 0.5*g1/c2;
	dst4[1][0] = 1-0.5*q2*g1/c2;			dst4[1][1] = g1*u/c2;				dst4[1][2] = g1*v/c2;				dst4[1][3] = -g1/c2;
	dst4[2][0] = 0.5*(0.5*q2*g1/c2 - qn/c);	dst4[2][1] = -0.5*(g1*u/c2 - nx/c);	dst4[2][2] = -0.5*(g1*v/c2 - ny/c);	dst4[2][3] = 0.5*g1/c2;
	dst4[3][0] = -ql;						dst4[3][1] = lx;					dst4[3][2] = ly;					dst4[3][3] = 0;	
}

void FVM_TVD_IMPLICIT::multMtx4(double **dst4, double **srcA4, double **srcB4)
{
	double sum;
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
		{
			sum = 0;
			for (int k = 0; k < 4; ++k)
				sum += srcA4[i][k] * srcB4[k][j];
			dst4[i][j] = sum;
		}
}

void FVM_TVD_IMPLICIT::clearMtx4(double **mtx4)
{
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			mtx4[i][j] = 0;
}

void FVM_TVD_IMPLICIT::calcAP(double **dst4, double **rightEgnVecl4, double **egnVal4, double **leftEgnVecl4)
{
	double		**tempMtx4 = allocMtx4();
	
	//egnVal4 +.
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			dst4[i][j] = egnVal4[i][j] > 0 ? egnVal4[i][j] : 0;
	
	multMtx4(tempMtx4, rightEgnVecl4, dst4);
	multMtx4(dst4, tempMtx4, leftEgnVecl4);
	freeMtx4(tempMtx4);
}

void FVM_TVD_IMPLICIT::calcAM(double **dst4, double **rightEgnVecl4, double **egnVal4, double **leftEgnVecl4)
{
	double		**tempMtx4 = allocMtx4();
	
	//egnVal4 -.
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
			dst4[i][j] = egnVal4[i][j] < 0 ? egnVal4[i][j] : 0;
	
	multMtx4(tempMtx4, rightEgnVecl4, dst4);
	multMtx4(dst4, tempMtx4, leftEgnVecl4);
	freeMtx4(tempMtx4);
}

void FVM_TVD_IMPLICIT::calcRoeAverage(Param& average, Param pL, Param pR, double GAM)
{
	double WI;
	roe_orig(average.r, average.e, average.p, average.u, average.v, WI,
		pL.r, pL.p, pL.u, pL.v, 0.0,
		pR.r, pR.p, pR.u, pR.v, 0.0, GAM);
	average.cz = sqrt(GAM*average.p/average.r);
	average.E = average.e + 0.5*(average.u * average.u + average.v * average.v);
}

void FVM_TVD_IMPLICIT::reconstruct(int iCell, Param& cell, Param neighbor[3])
{
	assert(iCell >= 0 && iCell < grid.cCount);
	for (int j = 0; j < 3; ++j)
	{
		int iEdge = cellsEdges[iCell][j];
		if (grid.edges[iEdge].type == Edge::TYPE_INNER)
		{
			int c1	= grid.edges[iEdge].c1;
			int c2	= grid.edges[iEdge].c2;
			if (c1 != iCell)	std::swap(c1, c2);
			assert(c1 == iCell);
			convertConsToPar(c1, cell);
			convertConsToPar(c2, neighbor[j]);
		} else {
			int c1	= grid.edges[iEdge].c1;
			assert(c1 == iCell);
			convertConsToPar(c1, cell);
			boundaryCond(iEdge, cell, neighbor[j]);
		}
	}
}

void FVM_TVD_IMPLICIT::run() 
{
	int						nc = grid.cCount; // количество ячеек.
	int						ne = grid.eCount; // количество ребер.
	double					t = 0.0;
	unsigned int			step = 0;

	MatrixSolver			*solverMtx = new SolverZeidel();
	double					**eigenMtx4, **rEigenVector4, **lEigenVector4, **A1mtx4, **A2mtx4, **mtx4_1,**mtx4_2;
	double					*right4;

	right4 = new double [4];
	eigenMtx4 = allocMtx4(); 
	rEigenVector4 = allocMtx4(); 
	lEigenVector4 = allocMtx4(); 
	A1mtx4 = allocMtx4(); 
	A2mtx4 = allocMtx4();
	mtx4_1 = allocMtx4();
	mtx4_2 = allocMtx4();
/*
	memcpy(ro, ro_old, nc*sizeof(double));
	memcpy(ru, ru_old, nc*sizeof(double));
	memcpy(rv, rv_old, nc*sizeof(double));
	memcpy(re, re_old, nc*sizeof(double));
*/
	solverMtx->init(nc, 4);

	while (t < TMAX)
	{
		t += TAU;
		++step;

		//заполнение матрицы.
		for (int iCell = 0; iCell < nc; ++iCell)
		{
			Param		average, cell, neighbor[3];
			double		__GAM = 1.4; // TODO: сделать правильное вычисление показателя адиабаты
			
			reconstruct(iCell, cell, neighbor);

			clearMtx4(mtx4_1);
			for (int neighborIndex = 0; neighborIndex < 3; ++neighborIndex)
			{
				int numberOfEdge = cellsEdges[iCell][neighborIndex];
				double l = grid.edges[numberOfEdge].l;
				Vector n = grid.edges[numberOfEdge].n;
				
				calcRoeAverage(average, cell, neighbor[neighborIndex], __GAM);
				double H = average.E + average.p/average.r;

				eigenValues(eigenMtx4, average.cz, average.u, n.x, average.v, 0.0);
				rightEigenVector(rEigenVector4, average.cz, average.u, n.x, average.v, 0.0, H);
				leftEigenVector(lEigenVector4, average.cz, __GAM, average.u, n.x, average.v, 0.0);
				calcAP(A1mtx4, rEigenVector4, eigenMtx4, lEigenVector4);

				eigenValues(eigenMtx4, average.cz, average.u, 0.0, average.v, n.y);
				rightEigenVector(rEigenVector4, average.cz, average.u, 0.0, average.v, n.y, H);
				leftEigenVector(lEigenVector4, average.cz, __GAM, average.u, 0.0, average.v, n.y);
				calcAP(A2mtx4, rEigenVector4, eigenMtx4, lEigenVector4);

				for (int i = 0; i < 4; ++i)
				{
					for (int j = 0; j < 4; ++j)
					{
						mtx4_1[i][j] += (A1mtx4[i][j] + A2mtx4[i][j])*l;
					}
				}

				eigenValues(eigenMtx4, average.cz, average.u, n.x, average.v, 0.0);
				rightEigenVector(rEigenVector4, average.cz, average.u, n.x, average.v, 0.0, H);
				leftEigenVector(lEigenVector4, average.cz, __GAM, average.u, n.x, average.v, 0.0);
				calcAM(A1mtx4, rEigenVector4, eigenMtx4, lEigenVector4);

				eigenValues(eigenMtx4, average.cz, average.u, 0.0, average.v, n.y);
				rightEigenVector(rEigenVector4, average.cz, average.u, 0.0, average.v, n.y, H);
				leftEigenVector(lEigenVector4, average.cz, __GAM, average.u, 0.0, average.v, n.y);
				calcAM(A2mtx4, rEigenVector4, eigenMtx4, lEigenVector4);

				for (int i = 0; i < 4; ++i)
				{
					for (int j = 0; j < 4; ++j)
					{
						mtx4_2[i][j] = (A1mtx4[i][j] + A2mtx4[i][j])*l;
					}
				}
				int c1 = grid.edges[numberOfEdge].c1;
				int c2 = grid.edges[numberOfEdge].c2;
				int neight = (c1 == iCell)? c2 : c1; //номер соседней ячейки.
				solverMtx->setMatrElement(iCell, 
					neight,
					mtx4_2);	//du~(i,j,n+1)
				
				double		fr, fu, fv, fe;
				calcFlux(fr, fu, fv, fe, cell, neighbor[neighborIndex], n, __GAM);

				right4[0] -= l*fr;
				right4[1] -= l*fu;
				right4[2] -= l*fv;
				right4[3] -= l*fe;
				solverMtx->setRightElement(iCell, right4); //F(i,j,n)
			}

			for (int i = 0; i < 4; ++i)
			{
				for (int j = 0; j < 4; ++j)
				{
					if (i == j) mtx4_1[i][j] += grid.cells[iCell].S/TAU;
				}
			}
			solverMtx->setMatrElement(iCell, iCell, mtx4_1);	//du(i, n+1)
		}

		int maxIter = 1000;
		const double eps = 1.0e-5;
		
		solverMtx->solve(eps, maxIter);
		for (int cellIndex = 0, ind = 0; cellIndex < nc; cellIndex++, ind += 4)
		{
			ro[cellIndex] += solverMtx->x[ind+0];
			ru[cellIndex] += solverMtx->x[ind+1];
			rv[cellIndex] += solverMtx->x[ind+2];
			re[cellIndex] += solverMtx->x[ind+3];			
		}
		solverMtx->zero();

		if (step % FILE_SAVE_STEP == 0)
		{
			save(step);
		}
		if (step % PRINT_STEP == 0)
		{
			log("step: %d\t\ttime step: %.16f\n", step, t);
		}
	}

	freeMtx4(eigenMtx4); 
	freeMtx4(rEigenVector4); 
	freeMtx4(lEigenVector4); 
	freeMtx4(A1mtx4); 
	freeMtx4(A2mtx4);
	freeMtx4(mtx4_1);
	freeMtx4(mtx4_2);
	delete [] right4;
	delete solverMtx;
}

void FVM_TVD_IMPLICIT::save(int step)
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

void FVM_TVD_IMPLICIT::calcFlux(double& fr, double& fu, double& fv, double& fe, Param pL, Param pR, Vector n, double GAM)
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

void FVM_TVD_IMPLICIT::boundaryCond(int iEdge, Param& pL, Param& pR)
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


void FVM_TVD_IMPLICIT::done()
{
	delete[] ro;
	delete[] ru;
	delete[] rv;
	delete[] re;
/*
	delete[] ro_old;
	delete[] ru_old;
	delete[] rv_old;
	delete[] re_old;

	delete[] ro_int;
	delete[] ru_int;
	delete[] rv_int;
	delete[] re_int;
*/
	for (int iCell = 0; iCell < grid.cCount; ++iCell)
		delete cellsEdges[iCell];
	delete [] cellsEdges;
}




Region & FVM_TVD_IMPLICIT::getRegionByCellType(int type)
{
	for (int i = 0; i < regCount; i++)
	{
		if (regions[i].cellType == type) return regions[i];
	}
	log("ERROR: unknown cell type %d...\n", type);
	EXIT(1);
}


Region   &	FVM_TVD_IMPLICIT::getRegion	(int iCell)
{
	return getRegionByCellType( grid.cells[iCell].type );
}

Material &	FVM_TVD_IMPLICIT::getMaterial	(int iCell)
{
	Region & reg = getRegion(iCell);
	return materials[reg.matId];
}


void FVM_TVD_IMPLICIT::convertParToCons(int iCell, Param & par)
{
	ro[iCell] = par.r;
	ru[iCell] = par.r*par.u;
	rv[iCell] = par.r*par.v;
	re[iCell] = par.r*(par.e+0.5*(par.u*par.u+par.v*par.v));
}

void FVM_TVD_IMPLICIT::convertConsToPar(int iCell, Param & par)
{
	par.r = ro[iCell];
	par.u = ru[iCell]/ro[iCell];
	par.v = rv[iCell]/ro[iCell];
	par.E = re[iCell]/ro[iCell];
	par.e = par.E-0.5*(par.u*par.u+par.v*par.v);
	Material& mat = getMaterial(iCell);
	mat.URS(par, 0);
}



