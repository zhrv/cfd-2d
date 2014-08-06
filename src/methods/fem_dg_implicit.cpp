#include "fem_dg_implicit.h"
#include "tinyxml.h"
#include "global.h"
#include <ctime>

#define POW_2(x) ((x)*(x))

FEM_DG_IMPLICIT::FEM_DG_IMPLICIT()
{
}


FEM_DG_IMPLICIT::~FEM_DG_IMPLICIT()
{
}


void FEM_DG_IMPLICIT::init(char * xmlFileName)
{
	TiXmlDocument doc(xmlFileName);
	bool loadOkay = doc.LoadFile(TIXML_ENCODING_UTF8);
	if (!loadOkay)
	{
		log("ERROR: %s\n", doc.ErrorDesc());
		exit(doc.ErrorId());
	}

	TiXmlNode* task = 0;
	TiXmlElement* el = 0;
	TiXmlNode* node0 = 0;
	TiXmlNode* node1 = 0;
	task = doc.FirstChild("task");

	int steadyVal = 1;
	node0 = task->FirstChild("control");
	node0->FirstChild("STEADY")->ToElement()->Attribute("value", &steadyVal);
	node0->FirstChild("TAU")->ToElement()->Attribute("value", &TAU);
	node0->FirstChild("TMAX")->ToElement()->Attribute("value", &TMAX);
	node0->FirstChild("STEP_MAX")->ToElement()->Attribute("value", &STEP_MAX);
	node0->FirstChild("FILE_OUTPUT_STEP")->ToElement()->Attribute("value", &FILE_SAVE_STEP);
	node0->FirstChild("LOG_OUTPUT_STEP")->ToElement()->Attribute("value", &PRINT_STEP);
	const char * flxStr = node0->FirstChild("FLUX")->ToElement()->Attribute("value");
	if (strcmp(flxStr, "GODUNOV") == 0) {
		FLUX = FLUX_GODUNOV;
	}
	else if (strcmp(flxStr, "LAX") == 0) {
		FLUX = FLUX_LAX;
	}
	else {
		FLUX = FLUX_GODUNOV;
	}

	if (steadyVal == 0) {
		STEADY = false;
	}
	else {
		STEADY = true;
		node1 = node0->FirstChild("CFL");
		node1->FirstChild("start")->ToElement()->Attribute("value", &CFL);
		node1->FirstChild("scale")->ToElement()->Attribute("value", &scaleCFL);
		node1->FirstChild("max")->ToElement()->Attribute("value", &maxCFL);
		node1->FirstChild("step")->ToElement()->Attribute("value", &stepCFL);
		node1->FirstChild("max_limited_cells")->ToElement()->Attribute("value", &maxLimCells);
	}

	// сглаживание невязок
	int smUsing = 1;
	node0 = task->FirstChild("smoothing");
	node0->FirstChild("using")->ToElement()->Attribute("value", &smUsing);
	node0->FirstChild("coefficient")->ToElement()->Attribute("value", &SMOOTHING_PAR);
	SMOOTHING = (smUsing == 1);

	// чтение параметров о ПРЕДЕЛЬНЫХ ЗНАЧЕНИЯХ
	node0 = task->FirstChild("limits");
	node0->FirstChild("ro")->ToElement()->Attribute("min", &limitRmin);
	node0->FirstChild("ro")->ToElement()->Attribute("max", &limitRmax);
	node0->FirstChild("p")->ToElement()->Attribute("min", &limitPmin);
	node0->FirstChild("p")->ToElement()->Attribute("max", &limitPmax);
	node0->FirstChild("u")->ToElement()->Attribute("max", &limitUmax);

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
		node1->FirstChild("M")->ToElement()->Attribute("value", &mat.M);
		node1->FirstChild("Cp")->ToElement()->Attribute("value", &mat.Cp);
		node1->FirstChild("K")->ToElement()->Attribute("value", &mat.K);
		node1->FirstChild("ML")->ToElement()->Attribute("value", &mat.ML);
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
		node1->FirstChild("Vx")->ToElement()->Attribute("value", &reg.par.u);
		node1->FirstChild("Vy")->ToElement()->Attribute("value", &reg.par.v);
		node1->FirstChild("T")->ToElement()->Attribute("value", &reg.par.T);
		node1->FirstChild("P")->ToElement()->Attribute("value", &reg.par.p);

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
		}
		else
		if (strcmp(str, "BOUND_OUTLET") == 0)
		{
			b.parCount = 0;
			b.par = NULL;
			b.type = Boundary::BOUND_OUTLET;
		}
		else
		if (strcmp(str, "BOUND_INLET") == 0)
		{
			b.parCount = 4;
			b.par = new double[4];
			b.type = Boundary::BOUND_INLET;

			node1 = bNode->FirstChild("parameters");
			node1->FirstChild("T")->ToElement()->Attribute("value", &b.par[0]);
			node1->FirstChild("P")->ToElement()->Attribute("value", &b.par[1]);
			node1->FirstChild("Vx")->ToElement()->Attribute("value", &b.par[2]);
			node1->FirstChild("Vy")->ToElement()->Attribute("value", &b.par[3]);
		}
		else {
			log("ERROR: unsupported boundary condition type '%s'", str);
			EXIT(1);
		}

		bNode = bNode->NextSibling("boundCond");
	}


	node0 = task->FirstChild("mesh");
	const char* fName = task->FirstChild("mesh")->FirstChild("name")->ToElement()->Attribute("value");
	grid.initFromFiles((char*)fName);

	memAlloc();

	// вычисляем узлы и коэффициенты квадратур
	// для ячеек
	for (int i = 0; i < grid.cCount; i++) {

		double a = 1.0 / 3.0;
		double b = 1.0 / 5.0;
		double c = 3.0 / 5.0;
		double x1 = grid.nodes[grid.cells[i].nodesInd[0]].x;
		double y1 = grid.nodes[grid.cells[i].nodesInd[0]].y;
		double x2 = grid.nodes[grid.cells[i].nodesInd[1]].x;
		double y2 = grid.nodes[grid.cells[i].nodesInd[1]].y;
		double x3 = grid.nodes[grid.cells[i].nodesInd[2]].x;
		double y3 = grid.nodes[grid.cells[i].nodesInd[2]].y;
		double a1 = x1 - x3;
		double a2 = y1 - y3;
		double b1 = x2 - x3;
		double b2 = y2 - y3;
		double c1 = x3;
		double c2 = y3;

		cellGW[i][0] = -27.0 / 48.0;
		cellGP[i][0].x = a1*a + b1*a + c1;
		cellGP[i][0].y = a2*a + b2*a + c2;

		cellGW[i][1] = 25.0 / 48.0;
		cellGP[i][1].x = a1*c + b1*b + c1;
		cellGP[i][1].y = a2*c + b2*b + c2;

		cellGW[i][2] = 25.0 / 48.0;
		cellGP[i][2].x = a1*b + b1*c + c1;
		cellGP[i][2].y = a2*b + b2*c + c2;

		cellGW[i][3] = 25.0 / 48.0;
		cellGP[i][3].x = a1*b + b1*b + c1;
		cellGP[i][3].y = a2*b + b2*b + c2;


		cellJ[i] = 0.5*fabs(a1*b2 - a2*b1);


		//double a = 1.0/6.0;
		//double b = 2.0/3.0;
		//double x1 = nodes[cellNodes[iCell][0]].x;
		//double y1 = nodes[cellNodes[iCell][0]].y;
		//double x2 = nodes[cellNodes[iCell][1]].x;
		//double y2 = nodes[cellNodes[iCell][1]].y;
		//double x3 = nodes[cellNodes[iCell][2]].x;
		//double y3 = nodes[cellNodes[iCell][2]].y;
		//double a1 = x1-x3;
		//double a2 = y1-y3;
		//double b1 = x2-x3;
		//double b2 = y2-y3;
		//double c1 = x3;
		//double c2 = y3;

		//cellWGP[0] = 1.0/6.0; 
		//cellGP[iCell][0].x = a1*a+b1*a+c1;
		//cellGP[iCell][0].y = a2*a+b2*a+c2;

		//cellWGP[1] = 1.0/6.0; 
		//cellGP[iCell][1].x = a1*a+b1*b+c1;
		//cellGP[iCell][1].y = a2*a+b2*b+c2;

		//cellWGP[2] = 1.0/6.0; 
		//cellGP[iCell][2].x = a1*b+b1*a+c1;
		//cellGP[iCell][2].y = a2*b+b2*a+c2;

		//cellJ[iCell] = a1*b2-a2*b1;


		//double a = 0.5;
		//double b = 0.0;
		//double x1 = nodes[cellNodes[iCell][0]].x;
		//double y1 = nodes[cellNodes[iCell][0]].y;
		//double x2 = nodes[cellNodes[iCell][1]].x;
		//double y2 = nodes[cellNodes[iCell][1]].y;
		//double x3 = nodes[cellNodes[iCell][2]].x;
		//double y3 = nodes[cellNodes[iCell][2]].y;
		//double a1 = x1-x3;
		//double a2 = y1-y3;
		//double b1 = x2-x3;
		//double b2 = y2-y3;
		//double c1 = x3;
		//double c2 = y3;

		//cellWGP[0] = 1.0/3.0; 
		//cellGP[iCell][0].x = a1*a+b1*a+c1;
		//cellGP[iCell][0].y = a2*a+b2*a+c2;

		//cellWGP[1] = 1.0/3.0; 
		//cellGP[iCell][1].x = a1*a+b1*b+c1;
		//cellGP[iCell][1].y = a2*a+b2*b+c2;

		//cellWGP[2] = 1.0/3.0; 
		//cellGP[iCell][2].x = a1*b+b1*a+c1;
		//cellGP[iCell][2].y = a2*b+b2*a+c2;

		//cellJ[iCell] = (a1*b2-a2*b1)*0.5;
	}

	// для ребер
	for (int i = 0; i < grid.eCount; i++) {
		double gp1 = -3.0 / 5.0;
		double gp2 = 0.0;
		double gp3 = 3.0 / 5.0;
		double x1 = grid.nodes[grid.edges[i].n1].x;
		double y1 = grid.nodes[grid.edges[i].n1].y;
		double x2 = grid.nodes[grid.edges[i].n2].x;
		double y2 = grid.nodes[grid.edges[i].n2].y;

		edgeGW[i][0] = 5.0 / 9.0;
		edgeGP[i][0].x = ((x1 + x2) + gp1*(x2 - x1)) / 2.0;
		edgeGP[i][0].y = ((y1 + y2) + gp1*(y2 - y1)) / 2.0;

		edgeGW[i][1] = 8.0 / 9.0;
		edgeGP[i][1].x = (x1 + x2) / 2.0;
		edgeGP[i][1].y = (y1 + y2) / 2.0;

		edgeGW[i][2] = 5.0 / 9.0;
		edgeGP[i][2].x = ((x1 + x2) + gp3*(x2 - x1)) / 2.0;
		edgeGP[i][2].y = ((y1 + y2) + gp3*(y2 - y1)) / 2.0;

		edgeJ[i] = sqrt(POW_2(x2 - x1) + POW_2(y2 - y1))*0.5;

	}

	// вычисляем матрицу масс
	for (int iCell = 0; iCell < grid.cCount; iCell++) {
		double **	A		= matrA[iCell];
		double **	invA	= matrInvA[iCell];
		int			N		= BASE_FUNC_COUNT;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				A[i][j] = 0.0;
				for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++)
				{
					A[i][j] += cellGW[iCell][iGP] * getF(i, iCell, cellGP[iCell][iGP].x, cellGP[iCell][iGP].y)*getF(j, iCell, cellGP[iCell][iGP].x, cellGP[iCell][iGP].y);
				}
				A[i][j] *= cellJ[iCell];
			}
		}

		inverseMatr(A, invA, N);
	}


	tmpArr = new double[grid.cCount];
	tmpArrInt = new int[grid.cCount];

	for (int i = 0; i < grid.cCount; i++)
	{
		Region & reg = getRegion(i);
		convertParToCons(i, reg.par);
	}

	calcTimeStep();

	solverMtx = new SolverHYPREBoomerAMG();
	solverMtx->init(grid.cCount, 4*BASE_FUNC_COUNT);

	save(0);
}


void FEM_DG_IMPLICIT::memAlloc()
{
	int n = grid.cCount;
	cTau = new double[n];

	ro = new double*[n];
	ru = new double*[n];
	rv = new double*[n];
	re = new double*[n];

	cellGP = new Point*[n];
	edgeGP = new Point*[n];

	cellGW = new double*[n];
	edgeGW = new double*[n];
	cellJ  = new double[n];
	edgeJ  = new double[n];

	matrA		= new double**[n];
	matrInvA	= new double**[n];

	for (int i = 0; i < grid.cCount; i++) {
		ro[i] = new double[BASE_FUNC_COUNT];
		ru[i] = new double[BASE_FUNC_COUNT];
		rv[i] = new double[BASE_FUNC_COUNT];
		re[i] = new double[BASE_FUNC_COUNT];

		cellGP[i] = new Point[GP_CELL_COUNT];
		edgeGP[i] = new Point[GP_EDGE_COUNT];

		cellGW[i] = new double[GP_CELL_COUNT];
		edgeGW[i] = new double[GP_EDGE_COUNT];

		matrA[i]	= new double*[BASE_FUNC_COUNT];
		matrInvA[i]	= new double*[BASE_FUNC_COUNT];
		for (int j = 0; j < BASE_FUNC_COUNT; j++) {
			matrA[i][j] = new double[BASE_FUNC_COUNT];
			matrInvA[i][j] = new double[BASE_FUNC_COUNT];
		}
	}

	fields = new double**[4];
	fields[FIELD_RO] = ro;
	fields[FIELD_RU] = ru;
	fields[FIELD_RV] = rv;
	fields[FIELD_RE] = re;
}


Region & FEM_DG_IMPLICIT::getRegionByCellType(int type)
{
	for (int i = 0; i < regCount; i++)
	{
		if (regions[i].cellType == type) return regions[i];
	}
	log("ERROR: unknown cell type %d...\n", type);
	EXIT(1);
}


Region   &	FEM_DG_IMPLICIT::getRegion(int iCell)
{
	return getRegionByCellType(grid.cells[iCell].type);
}

Material &	FEM_DG_IMPLICIT::getMaterial(int iCell)
{
	Region & reg = getRegion(iCell);
	return materials[reg.matId];
}


void FEM_DG_IMPLICIT::convertParToCons(int iCell, Param & par)
{
	//ro[iCell] = par.r;
	//ru[iCell] = par.r*par.u;
	//rv[iCell] = par.r*par.v;
	//re[iCell] = par.r*(par.e + 0.5*(par.u*par.u + par.v*par.v));
}

void FEM_DG_IMPLICIT::convertConsToPar(int iCell, Param & par)
{
	double fRO = getField(FIELD_RO, iCell, grid.cells[iCell].c);
	double fRU = getField(FIELD_RU, iCell, grid.cells[iCell].c);
	double fRV = getField(FIELD_RV, iCell, grid.cells[iCell].c);
	double fRE = getField(FIELD_RE, iCell, grid.cells[iCell].c);

	par.r = fRO;
	par.u = fRU / fRO;
	par.v = fRV / fRO;
	par.E = fRE / fRO;
	par.e = par.E - 0.5*(par.u*par.u + par.v*par.v);
	Material& mat = getMaterial(iCell);
	mat.URS(par, 0);
}

double FEM_DG_IMPLICIT::getField(int fldId, int iCell, double x, double y)
{
	double *fld = fields[fldId][iCell];
	double result = 0.0;
	for (int i = 0; i < BASE_FUNC_COUNT; i++) {
		result += fld[i] * getF(i, iCell, x, y);
	}
	return result;
}

double FEM_DG_IMPLICIT::getField(int fldId, int iCell, Point p)
{
	return getField(fldId, iCell, p.x, p.y);
}

double FEM_DG_IMPLICIT::getF(int id, int iCell, Point p)
{
	return getF(id, iCell, p.x, p.y);
}

double FEM_DG_IMPLICIT::getF(int id, int iCell, double x, double y)
{
	switch (id) {
	case 0:
		return 1.0;
		break;
	case 1:
		return (x - grid.cells[iCell].c.x) / grid.cells[iCell].HX;
		break;
	case 2:
		return (y - grid.cells[iCell].c.y) / grid.cells[iCell].HY;
		break;
	case 3:
		return (x - grid.cells[iCell].c.x)*(x - grid.cells[iCell].c.x) / grid.cells[iCell].HX / grid.cells[iCell].HX;
		break;
	case 4:
		return (y - grid.cells[iCell].c.y)*(y - grid.cells[iCell].c.y) / grid.cells[iCell].HY / grid.cells[iCell].HY;
		break;
	case 5:
		return (x - grid.cells[iCell].c.x)*(y - grid.cells[iCell].c.y) / grid.cells[iCell].HX / grid.cells[iCell].HY;
		break;
	default:
		return 0.0;
		break;
	}
}

double FEM_DG_IMPLICIT::getDfDx(int id, int iCell, Point p)
{
	return getDfDx(id, iCell, p.x, p.y);
}

double FEM_DG_IMPLICIT::getDfDx(int id, int iCell, double x, double y)
{
	switch (id) {
	case 0:
		return 0.0;
		break;
	case 1:
		return 1.0 / grid.cells[iCell].HX;
		break;
	case 2:
		return 0.0;
		break;
	case 3:
		return 2.0*(x - grid.cells[iCell].c.x) / grid.cells[iCell].HX / grid.cells[iCell].HX;
		break;
	case 4:
		return 0.0;
		break;
	case 5:
		return (y - grid.cells[iCell].c.y) / grid.cells[iCell].HX / grid.cells[iCell].HY;
		break;
	default:
		return 0.0;
		break;
	}
}

double FEM_DG_IMPLICIT::getDfDy(int id, int iCell, Point p)
{
	return getDfDy(id, iCell, p.x, p.y);
}

double FEM_DG_IMPLICIT::getDfDy(int id, int iCell, double x, double y)
{
	switch (id) {
	case 0:
		return 0.0;
		break;
	case 1:
		return 0.0;
		break;
	case 2:
		return 1.0 / grid.cells[iCell].HY;
		break;
	case 3:
		return 0.0;
		break;
	case 4:
		return 2.0*(y - grid.cells[iCell].c.y) / grid.cells[iCell].HY / grid.cells[iCell].HY;
		break;
	case 5:
		return (x - grid.cells[iCell].c.x) / grid.cells[iCell].HX / grid.cells[iCell].HY;
		break;
	default:
		return 0.0;
		break;
	}
}

void FEM_DG_IMPLICIT::calcTimeStep()
{

}

void FEM_DG_IMPLICIT::save(int step)
{
	char fName[50];

	sprintf(fName, "res_%010d.vtk", step);
	FILE * fp = fopen(fName, "w");
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "GASDIN data file\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

	fprintf(fp, "POINTS %d float\n", grid.nCount);
	for (int i = 0; i < grid.nCount; i++)
	{
		fprintf(fp, "%f %f %f  ", grid.nodes[i].x, grid.nodes[i].y, 0.0);
		if (i + 1 % 8 == 0) fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	fprintf(fp, "CELLS %d %d\n", grid.cCount, 4 * grid.cCount);
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
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Pressure float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.p);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Temperature float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.T);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS MachNumber float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", sqrt(p.u*p.u + p.v*p.v) / p.cz);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "VECTORS Velosity float\n");
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f %25.16f %25.16f ", p.u, p.v, 0.0);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Total_energy float 1\nLOOKUP_TABLE default\n");
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.E);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Total_pressure float 1\nLOOKUP_TABLE default\n");
	for (int i = 0; i < grid.cCount; i++)
	{
		Material &mat = getMaterial(i);
		double gam = mat.getGamma();
		double agam = gam - 1.0;
		Param p;
		convertConsToPar(i, p);
		double M2 = (p.u*p.u + p.v*p.v) / (p.cz*p.cz);
		fprintf(fp, "%25.16e ", p.p*::pow(1.0 + 0.5*M2*agam, gam / agam));
		if ((i + 1) % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS TAU float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%25.16f ", cTau[i]);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fclose(fp);
	printf("File '%s' saved...\n", fName);

}

int FEM_DG_IMPLICIT::getLimitedCellsCount() {
	int n = 0;
	for (int iCell = 0; iCell < grid.cCount; iCell++) {
		if ((grid.cells[iCell].flag & CELL_FLAG_LIM) > 0) n++;
	}
	return n;
}

void FEM_DG_IMPLICIT::remediateLimCells()
{

}


void FEM_DG_IMPLICIT::decCFL()
{
	if (CFL > 0.1) {
		CFL *= 0.75;
		if (CFL < 0.1) CFL = 0.1;
		log(" CFL Number has been decreased : %f \n", CFL);
	}
}

void FEM_DG_IMPLICIT::incCFL()
{
	if (CFL < maxCFL) {
		CFL *= scaleCFL;
		if (CFL > maxCFL) CFL = maxCFL;
		log(" CFL Number has been increased : %f \n", CFL);
	}
}

void FEM_DG_IMPLICIT::calcIntegral()
{

}

void FEM_DG_IMPLICIT::calcMatrWithTau()
{

}

void FEM_DG_IMPLICIT::calcMatrFlux()
{

}

void FEM_DG_IMPLICIT::calcRHS()
{

}


void FEM_DG_IMPLICIT::run()
{
	double __GAM = 1.4;
	double t = 0.0;
	int step = 0;
	while (t < TMAX && step < STEP_MAX) {
		long timeStart, timeEnd;
		timeStart = clock();

		if (STEADY) {
			calcTimeStep();
		}
		else {
			t + TAU;
		}

		int solverErr = 0;
		solverMtx->zero();

		/* Заполняем элементы матрицы */
		calcIntegral();			// вычисляем интеграл от(dF / dU)*deltaU*dFi / dx
		calcMatrWithTau();		// вычисляем матрицы перед производной по времени
		calcMatrFlux();			// Вычисляем потоковые величины 


		/* Заполняем правую часть */
		calcRHS();				// Вычисляем столбец правых членов


		/* Решаем СЛАУ */
		int maxIter = 20;
		const double eps = 1.0e-7;

		solverErr = solverMtx->solve(eps, maxIter);

		if (solverErr == MatrixSolver::RESULT_OK) {

			//if (SMOOTHING) smoothingDelta(solverMtx->x);

			for (int cellIndex = 0, ind = 0; cellIndex < grid.cCount; cellIndex++)
			{
				//if (cellIsLim(cellIndex))	continue;
				for (int iFld = 0; iFld < 4; iFld++) {
					for (int iF = 0; iF < BASE_FUNC_COUNT; iF++) {
						fields[iFld][cellIndex][iF] += solverMtx->x[ind++];
					}
				}

				Param par;
				convertConsToPar(cellIndex, par);
				if (par.r > limitRmax)			{ par.r = limitRmax;	setCellFlagLim(cellIndex); }
				if (par.r < limitRmin)			{ par.r = limitRmin;	setCellFlagLim(cellIndex); }
				if (fabs(par.u) > limitUmax)		{ par.u = limitUmax*par.u / fabs(par.u);	setCellFlagLim(cellIndex); }
				if (fabs(par.v) > limitUmax)		{ par.v = limitUmax*par.v / fabs(par.v);	setCellFlagLim(cellIndex); }
				if (par.p > limitPmax)			{ par.p = limitPmax;	setCellFlagLim(cellIndex); }
				if (par.p < limitPmin)			{ par.p = limitPmin;	setCellFlagLim(cellIndex); }
				if (cellIsLim(cellIndex)) 		{ par.e = par.p / ((__GAM - 1)*par.r); convertParToCons(cellIndex, par); }
			}
			remediateLimCells();
			int limCells = getLimitedCellsCount();
			if (STEADY && (limCells >= maxLimCells)) decCFL();

			timeEnd = clock();

			if (step % FILE_SAVE_STEP == 0)
			{
				save(step);
			}
			if (step % PRINT_STEP == 0)
			{
				//calcLiftForce();
				if (!STEADY) {

					log("step: %d\ttime step: %.16f\tmax iter: %d\tlim: %d \ttime: %d ms\n", step, t, maxIter, limCells, timeEnd - timeStart);
				}
				else {
					log("step: %d\tmax iter: %d\tlim: %d ttime: %d ms\n", step, maxIter, limCells, timeEnd - timeStart);
				}
			}

			if (STEADY && (step % stepCFL == 0)) incCFL();
		}
		else {
			if (STEADY) {
				//decCFL();
			}
			else {
				solverErr = 0;
			}
		}
	}
}


void FEM_DG_IMPLICIT::done()
{
}


