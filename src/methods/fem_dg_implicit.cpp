#include "fem_dg_implicit.h"
#include "tinyxml.h"
#include "global.h"
#include <ctime>
#include <cfloat>
#include "LimiterDG.h"
#include "MatrixSolver.h"

#define POW_2(x) ((x)*(x))


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
	node0->FirstChild("SIGMA")->ToElement()->Attribute("value", &SIGMA);
	node0->FirstChild("TAU")->ToElement()->Attribute("value", &TAU);
	node0->FirstChild("TMAX")->ToElement()->Attribute("value", &TMAX);
	node0->FirstChild("STEP_MAX")->ToElement()->Attribute("value", &STEP_MAX);
	node0->FirstChild("FILE_OUTPUT_STEP")->ToElement()->Attribute("value", &FILE_SAVE_STEP);
	node0->FirstChild("LOG_OUTPUT_STEP")->ToElement()->Attribute("value", &PRINT_STEP);

	const char * limiterName = node0->FirstChild("LIMITER")->ToElement()->Attribute("value");

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
	double maxP = 0.0;
	double maxR = 0.0;
	double maxT = 0.0;
	double maxU = 0.0;

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
		reg.par.ML = mat.ML;

		regNode = regNode->NextSibling("region");

		if (reg.par.p > maxP) maxP = reg.par.p;
		if (reg.par.r > maxR) maxR = reg.par.r;
		if (reg.par.T > maxT) maxT = reg.par.T;
		if (reg.par.U2() > maxU) maxU = reg.par.U2();
	}

	maxU = sqrt(maxU);

	// параметры обезразмеривания
	L_ = 1.0;
	R_ = maxR;					// характерная плотность = начальная плотность
	P_ = maxP;					// характерное давление = начальное давление
	T_ = maxT;					// характерная температура = начальная температура
	U_ = sqrt(P_ / R_);		// характерная скорость = sqrt( P_ / R_ )
	E_ = POW_2(U_);			// характерная энергия  = U_**2
	CV_ = POW_2(U_) / T_;	// характерная теплоёмкость  = U_**2 / T_
	TIME_ = L_ / U_;			// характерное время
	MU_ = R_ * U_ * L_;		// характерная вязкость = R_ * U_ * L_
	KP_ = R_ * POW_2(U_) * U_ * L_ / T_;	// коэффициент теплопроводности = R_ * U_**3 * L_ / T_
	CV_ = POW_2(U_) / T_;	// характерная теплоёмкость  = U_**2 / T_


	// Обезразмеривание всех параметров
	TAU /= TIME_;
	TMAX /= TIME_;

	for (int i = 0; i < regCount; i++) {
		Region &reg = regions[i];
		Param &par = reg.par;

		Material& mat = materials[reg.matId];
		mat.URS(reg.par, 2);	// r=r(p,T)
		mat.URS(reg.par, 1);	// e=e(p,r)

		par.p /= P_;
		par.u /= U_;
		par.v /= U_;
		par.T /= T_;

		par.r /= R_;
		par.e /= E_;
		par.cz /= U_;
		par.ML /= MU_;
		par.E = par.e + par.U2()*0.5;
	}

	Material::gR *= R_ * T_ / P_;	// Газовая постоянная
	for (int i = 0; i < matCount; i++)
	{
		Material & mat = materials[i];
		mat.Cp /= CV_;
		mat.K /= KP_;
		mat.ML /= MU_;
	}



	// чтение параметров о ГРАНИЧНЫХ УСЛОВИЯХ
    node0 = task->FirstChild("boundaries");
    TiXmlNode* bNode = node0->FirstChild("boundCond");
    while (bNode != NULL)
    {
        int edgeType;
        bNode->ToElement()->Attribute("edgeType", &edgeType);

        CFDBoundary * b;

        try {
            b = CFDBoundary::create(bNode, &grid);
        }
        catch (Exception e) {
            log("ERROR: %s\n", e.getMessage());
            exit(e.getType());
        }

        // TODO: !!!! ??????? ??? ????????????????
        if (b->parCount == 4) {
            b->par[0] /= U_;
            b->par[1] /= U_;
            b->par[2] /= T_;
            b->par[3] /= P_;
        }

        boundaries.push_back(b);

        bNode = bNode->NextSibling("boundCond");
    }

    bCount = boundaries.size();
//	node0 = task->FirstChild("boundaries");
//	node0->ToElement()->Attribute("count", &bCount);
//	//boundaries = new CFDBoundaries[bCount];
//	TiXmlNode* bNode = node0->FirstChild("boundCond");
//	for (int i = 0; i < bCount; i++)
//	{
//		CFDBoundary & b = boundaries[i];
//		bNode->ToElement()->Attribute("edgeType", &b.edgeType);
//		const char * str = bNode->FirstChild("type")->ToElement()->GetText();
//		if (strcmp(str, "BOUND_WALL") == 0)
//		{
//			b.parCount = 0;
//			b.par = NULL;
//			b.type = Boundary::BOUND_WALL;
//		}
//		else
//		if (strcmp(str, "BOUND_OUTLET") == 0)
//		{
//			b.parCount = 0;
//			b.par = NULL;
//			b.type = Boundary::BOUND_OUTLET;
//		}
//		else
//		if (strcmp(str, "BOUND_INLET") == 0)
//		{
//			b.parCount = 4;
//			b.par = new double[4];
//			b.type = Boundary::BOUND_INLET;
//
//			node1 = bNode->FirstChild("parameters");
//			node1->FirstChild("T")->ToElement()->Attribute("value", &b.par[0]); b.par[0] /= T_;
//			node1->FirstChild("P")->ToElement()->Attribute("value", &b.par[1]); b.par[1] /= P_;
//			node1->FirstChild("Vx")->ToElement()->Attribute("value", &b.par[2]); b.par[2] /= U_;
//			node1->FirstChild("Vy")->ToElement()->Attribute("value", &b.par[3]); b.par[3] /= U_;
//		}
//		else {
//			log("ERROR: unsupported boundary condition type '%s'", str);
//			EXIT(1);
//		}
//
//		bNode = bNode->NextSibling("boundCond");
//	}

	
	node0 = task->FirstChild("mesh");
	const char* fName = task->FirstChild("mesh")->FirstChild("name")->ToElement()->Attribute("value");
	grid.initFromFiles((char*)fName);

	memAlloc();

    /* определение ГУ для каждого ребра. */
    for (int iEdge = 0; iEdge < grid.eCount; iEdge++) {
        Edge & e = grid.edges[iEdge];
        if (e.type == Edge::TYPE_INNER) {
            e.bnd = NULL;
            continue;
        }
        if (e.type == Edge::TYPE_NAMED) {
            int iBound = -1;
            for (int i = 0; i < bCount; i++)
            {
                if (strcmp(e.typeName, boundaries[i]->name) == 0)
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

            e.bnd = boundaries[iBound];
        }
        else {
            int iBound = -1;
            for (int i = 0; i < bCount; i++)
            {
                if (e.type == boundaries[i]->edgeType)
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
            e.bnd = boundaries[iBound];
            strcpy(e.typeName, boundaries[iBound]->name);
        }
    }

	// инициализация лимитера
	limiter = LimiterDG::create(limiterName, this);
	if (limiter != NULL) {
		log("Limiter: %s\n", limiter->getName());
	}
	else {
		log("Without limiter\n");
	}

	// инициализация решателя
	node0 = task->FirstChild("solver");
	const char * solverName = node0->ToElement()->Attribute("type");
	//solverMtx = new SolverZeidel();
	solverMtx = MatrixSolver::create(solverName);
	solverMtx->init(grid.cCount, MATR_DIM);
	//solverMtx->init(&(grid), FIELD_COUNT, BASE_FUNC_COUNT, 1);
	log("Solver type: %s.\n", solverMtx->getName());

	node0->FirstChild("iterations")->ToElement()->Attribute("value", &SOLVER_ITER);
	node0->FirstChild("epsilon")->ToElement()->Attribute("value", &SOLVER_EPS);
	if (strstr(solverName, "HYPRE") != NULL) {
		node1 = node0->FirstChild("hypre");

		int tmp = 0;
		node1->FirstChild("printLevel")->ToElement()->Attribute("value", &tmp);
		solverMtx->setParameter("PRINT_LEVEL", tmp);
		node1->FirstChild("krylovDim")->ToElement()->Attribute("value", &tmp);
		solverMtx->setParameter("KRYLOV_DIM", tmp);
	}


	

	calcGaussPar();

	calcMassMatr();


	for (int i = 0; i < grid.cCount; i++)
	{
		Region & reg = getRegion(i);
		convertParToCons(i, reg.par);
        memset(tau_xx[i], 0, sizeof(double)*BASE_FUNC_COUNT);
        memset(tau_xy[i], 0, sizeof(double)*BASE_FUNC_COUNT);
        memset(tau_yy[i], 0, sizeof(double)*BASE_FUNC_COUNT);
	}

	calcTimeStep();
	//log("TAU_MIN = %25.16e\n", TAU_MIN);

	save(0);
}

void FEM_DG_IMPLICIT::calcMassMatr()
{
	for (int iCell = 0; iCell < grid.cCount; iCell++) {
		double **	A = matrA[iCell];
		for (int i = 0; i < BASE_FUNC_COUNT; i++) {
			for (int j = 0; j < BASE_FUNC_COUNT; j++) {
				A[i][j] = 0.0;
				for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
					A[i][j] += cellGW[iCell][iGP] * getF(i, iCell, cellGP[iCell][iGP].x, cellGP[iCell][iGP].y)
						                          * getF(j, iCell, cellGP[iCell][iGP].x, cellGP[iCell][iGP].y);
				}
				A[i][j] *= cellJ[iCell];
			}
		}
	}
}

void FEM_DG_IMPLICIT::calcGaussPar()
{
	
	// для ячеек
	if (GP_CELL_COUNT == 4) {
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

			//cellJ[i] = 0.5*fabs(a1*b2 - a2*b1);
			cellJ[i] = 0.5*(a1*b2 - a2*b1);
			if (cellJ[i] <= 0) {
				int zhrv = 0;
			}
		}
	}
	else if (GP_CELL_COUNT == 3) {
		for (int i = 0; i < grid.cCount; i++) {
			double a = 1.0 / 6.0;
			double b = 2.0 / 3.0;
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

			cellGW[i][0] = 1.0 / 6.0;
			cellGP[i][0].x = a1*a + b1*a + c1;
			cellGP[i][0].y = a2*a + b2*a + c2;

			cellGW[i][1] = 1.0 / 6.0;
			cellGP[i][1].x = a1*a + b1*b + c1;
			cellGP[i][1].y = a2*a + b2*b + c2;

			cellGW[i][2] = 1.0 / 6.0;
			cellGP[i][2].x = a1*b + b1*a + c1;
			cellGP[i][2].y = a2*b + b2*a + c2;

			cellJ[i] = a1*b2 - a2*b1;


		}
	}


	// для ребер
	if (GP_EDGE_COUNT == 3) {
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
	}
	else if (GP_EDGE_COUNT == 2) {
		for (int i = 0; i < grid.eCount; i++) {
			double gp1 = -1.0 / sqrt(3.0);
			double gp2 = 1.0 / sqrt(3.0);
			double x1 = grid.nodes[grid.edges[i].n1].x;
			double y1 = grid.nodes[grid.edges[i].n1].y;
			double x2 = grid.nodes[grid.edges[i].n2].x;
			double y2 = grid.nodes[grid.edges[i].n2].y;

			edgeGW[i][0] = 1.0;
			edgeGP[i][0].x = (x1 + x2) / 2.0 + gp1*(x2 - x1) / 2.0;
			edgeGP[i][0].y = (y1 + y2) / 2.0 + gp1*(y2 - y1) / 2.0;

			edgeGW[i][1] = 1.0;
			edgeGP[i][1].x = (x1 + x2) / 2.0 + gp2*(x2 - x1) / 2.0;
			edgeGP[i][1].y = (y1 + y2) / 2.0 + gp2*(y2 - y1) / 2.0;

			edgeJ[i] = sqrt(POW_2(x2 - x1) + POW_2(y2 - y1))*0.5;
		}
	}
}

double** allocMtx(int N)
{
	double **m;
	m = new double*[N];
	for (int i = 0; i < N; i++)
		m[i] = new double[N];

	return m;
}

void freeMtx(double** m, int N)
{
	for (int i = 0; i < N; i++)
		delete[] m[i];
	delete[] m;
}

void FEM_DG_IMPLICIT::memAlloc()
{
	int n = grid.cCount;
	cTau = new double[n];

	ro = new double*[n];
	ru = new double*[n];
	rv = new double*[n];
	re = new double*[n];
    tau_xx = new double*[n];
    tau_xy = new double*[n];
    tau_yy = new double*[n];

	cellGP = new Point*[n];
	cellGW = new double*[n];
	cellJ = new double[n];

	edgeGW = new double*[grid.eCount];
	edgeJ  = new double[grid.eCount];
	edgeGP = new Point*[grid.eCount];

	matrA		= new double**[n];
	matrInvA	= new double**[n];

	for (int i = 0; i < n; i++) {
		ro[i] = new double[BASE_FUNC_COUNT];
		ru[i] = new double[BASE_FUNC_COUNT];
		rv[i] = new double[BASE_FUNC_COUNT];
		re[i] = new double[BASE_FUNC_COUNT];
        tau_xx[i] = new double[BASE_FUNC_COUNT];
        tau_xy[i] = new double[BASE_FUNC_COUNT];
        tau_yy[i] = new double[BASE_FUNC_COUNT];

		cellGP[i] = new Point[GP_CELL_COUNT];
		cellGW[i] = new double[GP_CELL_COUNT];

		matrA[i] = allocMtx(BASE_FUNC_COUNT);
		matrInvA[i] = allocMtx(BASE_FUNC_COUNT);
	}

	for (int i = 0; i < grid.eCount; i++) {
		edgeGP[i] = new Point[GP_EDGE_COUNT];
		edgeGW[i] = new double[GP_EDGE_COUNT];
	}

	tmpArr = new double[n];
	tmpArr1 = new double[n];
	tmpArr2 = new double[n];
	tmpCFL = new double[n];
	tmpArrInt = new int[n];

	fields = new double**[FIELD_COUNT_EXT];
	fields[FIELD_RO] = ro;
	fields[FIELD_RU] = ru;
	fields[FIELD_RV] = rv;
	fields[FIELD_RE] = re;
	fields[FIELD_TAU_XX] = tau_xx;
	fields[FIELD_TAU_XY] = tau_xy;
	fields[FIELD_TAU_YY] = tau_yy;

	matrSmall = allocMtx(BASE_FUNC_COUNT);
	matrSmall2 = allocMtx(BASE_FUNC_COUNT);
	matrBig = allocMtx(MATR_DIM);
	matrBig2 = allocMtx(MATR_DIM);

}

void FEM_DG_IMPLICIT::memFree() 
{
	int n = grid.cCount; 
	
	for (int i = 0; i < n; i++) {
		delete[] ro[i];
		delete[] ru[i];
		delete[] rv[i];
		delete[] re[i];
        delete[] tau_xx[i];
        delete[] tau_xy[i];
        delete[] tau_yy[i];

		delete[] cellGP[i];
		delete[] cellGW[i];

		freeMtx(matrA[i], BASE_FUNC_COUNT);
		freeMtx(matrInvA[i], BASE_FUNC_COUNT);
	}
	delete[] cTau;

	delete[] ro;
	delete[] ru;
	delete[] rv;
	delete[] re;
    delete[] tau_xx;
    delete[] tau_xy;
    delete[] tau_yy;

	delete[] cellGP;
	delete[] cellGW;
	delete[] cellJ;

	delete[] matrA;
	delete[] matrInvA;

	for (int i = 0; i < grid.eCount; i++) {
		delete[] edgeGP[i];
		delete[] edgeGW[i];
	}
	delete[] edgeGW;
	delete[] edgeJ;
	delete[] edgeGP;

	delete[] tmpArr;
	delete[] tmpArr1;
	delete[] tmpArr2;
	delete[] tmpCFL;

	delete[] fields;

	freeMtx(matrSmall, BASE_FUNC_COUNT);
	freeMtx(matrSmall2, BASE_FUNC_COUNT);
	freeMtx(matrBig, MATR_DIM);
	freeMtx(matrBig2, MATR_DIM);

	delete[] tmpArr;
	delete[] tmpArr1;
	delete[] tmpArr2;
	delete[] tmpArrInt;


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
	memset(ro[iCell], 0, sizeof(double)*BASE_FUNC_COUNT);
	memset(ru[iCell], 0, sizeof(double)*BASE_FUNC_COUNT);
	memset(rv[iCell], 0, sizeof(double)*BASE_FUNC_COUNT);
	memset(re[iCell], 0, sizeof(double)*BASE_FUNC_COUNT);
	ro[iCell][0] = par.r;
	ru[iCell][0] = par.r*par.u;
	rv[iCell][0] = par.r*par.v;
	re[iCell][0] = par.r*(par.e + 0.5*(par.u*par.u + par.v*par.v));
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

void FEM_DG_IMPLICIT::getFields(double &fRO, double &fRU, double &fRV, double &fRE, int iCell, double x, double y)
{
	fRO = getField(FIELD_RO, iCell, x, y);
	fRU = getField(FIELD_RU, iCell, x, y);
	fRV = getField(FIELD_RV, iCell, x, y);
	fRE = getField(FIELD_RE, iCell, x, y);
}

void FEM_DG_IMPLICIT::getFields(double &fRO, double &fRU, double &fRV, double &fRE, int iCell, Point p)
{
	getFields(fRO, fRU, fRV, fRE, iCell, p.x, p.y);
}

double FEM_DG_IMPLICIT::getF(int id, int iCell, Point p)
{
	return getF(id, iCell, p.x, p.y);
}

double FEM_DG_IMPLICIT::getF(int id, int iCell, double x, double y)
{
	Point &c = grid.cells[iCell].c;
	double &hx = grid.cells[iCell].HX;
	double &hy = grid.cells[iCell].HY;

	switch (id) {
	case 0:
		return 1.0;
		break;
	case 1:
		return (x - c.x) / hx;
		break;
	case 2:
		return (y - c.y) / hy;
		break;
	case 3:
		return (x - c.x)*(x - c.x) / hx / hx;
		break;
	case 4:
		return (y - c.y)*(y - c.y) / hy / hy;
		break;
	case 5:
		return (x - c.x)*(y - c.y) / hx / hy;
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
	Point &c = grid.cells[iCell].c;
	double &hx = grid.cells[iCell].HX;
	double &hy = grid.cells[iCell].HY;

	switch (id) {
	case 0:
		return 0.0;
		break;
	case 1:
		return 1.0 / hx;
		break;
	case 2:
		return 0.0;
		break;
	case 3:
		return 2.0*(x - c.x) / hx / hx;
		break;
	case 4:
		return 0.0;
		break;
	case 5:
		return (y - c.y) / hx / hy;
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
	Point &c = grid.cells[iCell].c;
	double &hx = grid.cells[iCell].HX;
	double &hy = grid.cells[iCell].HY;

	switch (id) {
	case 0:
		return 0.0;
		break;
	case 1:
		return 0.0;
		break;
	case 2:
		return 1.0 / hy;
		break;
	case 3:
		return 0.0;
		break;
	case 4:
		return 2.0*(y - c.y) / hy / hy;
		break;
	case 5:
		return (x - c.x) / hx / hy;
		break;
	default:
		return 0.0;
		break;
	}
}

void FEM_DG_IMPLICIT::calcTimeStep()
{
	double tauMin = 1.0e+20;
	for (int iCell = 0; iCell < grid.cCount; iCell++)
	{
		if (STEADY) {
			Param p;
			convertConsToPar(iCell, p);
			double tmp = 0.5*grid.cells[iCell].S/(sqrt(p.u*p.u+p.v*p.v)+p.cz+FLT_EPSILON);
			cTau[iCell] = _min_(CFL*tmp, TAU);
			if (cTau[iCell] < tauMin) tauMin = cTau[iCell];
		} else {
			cTau[iCell] = TAU;
		}
	}
	if (STEADY)	TAU_MIN = tauMin;
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
		fprintf(fp, "%f %f %f  ", grid.nodes[i].x*L_, grid.nodes[i].y*L_, 0.0);
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
		fprintf(fp, "%25.16f ", p.r*R_);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Pressure float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.p*P_);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Temperature float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		Material &mat = getMaterial(i);
		mat.URS(p, 1);
		fprintf(fp, "%25.16f ", p.T*T_);
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
		fprintf(fp, "%25.16f %25.16f %25.16f ", p.u*U_, p.v*U_, 0.0);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Total_energy float 1\nLOOKUP_TABLE default\n");
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.E*E_);
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
		fprintf(fp, "%25.16e ", p.p*::pow(1.0 + 0.5*M2*agam, gam / agam)*P_);
		if ((i + 1) % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS TAU float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%25.16f ", cTau[i]*TIME_);
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
	for (int iCell = 0; iCell < grid.cCount; iCell++)
	{
		if (cellIsLim(iCell))
		{
			// пересчитываем по соседям			
			double sRO = 0.0;
			double sRU = 0.0;
			double sRV = 0.0;
			double sRE = 0.0;
			double S = 0.0;
			for (int i = 0; i < grid.cells[iCell].eCount; i++)
			{
				int		iEdge = grid.cells[iCell].edgesInd[i];
				int		j = grid.edges[iEdge].c2;
				if (j == iCell)	{
					//std::swap(j, grid.edges[iEdge].c1); // так нужно еще нормаль поворчивать тогда
					j = grid.edges[iEdge].c1;
				}
				if (j >= 0) {
					double  s = grid.cells[j].S;
					S += s;
					sRO += getField(FIELD_RO, j, grid.cells[j].c) * s;
					sRU += getField(FIELD_RO, j, grid.cells[j].c) * s;
					sRV += getField(FIELD_RO, j, grid.cells[j].c) * s;
					sRE += getField(FIELD_RO, j, grid.cells[j].c) * s;
				}
			}
			if (S >= TAU*TAU) {
				memset(ro[iCell], 0, sizeof(double)*BASE_FUNC_COUNT);
				memset(ru[iCell], 0, sizeof(double)*BASE_FUNC_COUNT);
				memset(rv[iCell], 0, sizeof(double)*BASE_FUNC_COUNT);
				memset(re[iCell], 0, sizeof(double)*BASE_FUNC_COUNT);

				ro[iCell][0] = sRO / S;
				ru[iCell][0] = sRU / S;
				rv[iCell][0] = sRV / S;
				re[iCell][0] = sRE / S;
			}
			// после 0x20 итераций пробуем вернуть ячейку в счет
			grid.cells[iCell].flag += 0x010000;
			if (grid.cells[iCell].flag & 0x200000) grid.cells[iCell].flag &= 0x001110;
		}
	}
}


void FEM_DG_IMPLICIT::decCFL()
{
	if (CFL > 0.01) {
		CFL *= 0.75;
		if (CFL < 0.01) CFL = 0.01;
		log(" CFL Number has been decreased : %25.25e \n", CFL);
	}
	//log(" <<< CFL Number has been decreased : %25.25e \n", CFL);
}

void FEM_DG_IMPLICIT::incCFL()
{
	if (CFL < maxCFL) {
		CFL *= scaleCFL;
		if (CFL > maxCFL) CFL = maxCFL;
		log(" CFL Number has been increased : %25.25e \n", CFL);
	}
}

double **FEM_DG_IMPLICIT::allocMtx4()
{
	double		**tempMtx4 = new double*[4];
	for (int i = 0; i < 4; ++i) tempMtx4[i] = new double[4];
	return tempMtx4;
}
void   FEM_DG_IMPLICIT::freeMtx4(double **mtx4)
{
	for (int i = 0; i < 4; ++i)
		delete[] mtx4[i];
	delete[] mtx4;
}
void FEM_DG_IMPLICIT::multMtx4(double **dst4, double **srcA4, double **srcB4)
{
	double sum;
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			sum = 0;
			for (int k = 0; k < 4; ++k)
				sum += srcA4[i][k] * srcB4[k][j];
			dst4[i][j] = sum;
		}
	}
}
void FEM_DG_IMPLICIT::clearMtx4(double **mtx4)
{
	for (int i = 0; i < 4; ++i)
	for (int j = 0; j < 4; ++j)
		mtx4[i][j] = 0;
}

void FEM_DG_IMPLICIT::multMtxToVal(double **dst, double x, int N)
{
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			dst[i][j] *= x;
		}
	}
}

void FEM_DG_IMPLICIT::fillMtx(double** dst, double x, int N)
{
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			dst[i][j] = x;
		}
	}
}

void FEM_DG_IMPLICIT::eigenValues(double** dst4, double c, double u, double nx, double v, double ny)
{
	clearMtx4(dst4);
	double qn = u*nx + v*ny;
	dst4[0][0] = qn - c;
	dst4[1][1] = qn;
	dst4[2][2] = qn + c;
	dst4[3][3] = qn;
}
void FEM_DG_IMPLICIT::rightEigenVector(double **dst4, double c, double u, double nx, double v, double ny, double H)
{
	double qn = u*nx + v*ny;
	double q2 = u*u + v*v;
	double lx = -ny;
	double ly = nx;
	double ql = u*lx + v*ly;
	dst4[0][0] = 1.0;			dst4[0][1] = 1.0;		dst4[0][2] = 1.0;			dst4[0][3] = 0.0;
	dst4[1][0] = u - c*nx;	dst4[1][1] = u;		dst4[1][2] = u + c*nx;	dst4[1][3] = lx;
	dst4[2][0] = v - c*ny;	dst4[2][1] = v;		dst4[2][2] = v + c*ny;	dst4[2][3] = ly;
	dst4[3][0] = H - qn*c;	dst4[3][1] = q2 / 2;	dst4[3][2] = H + qn*c;	dst4[3][3] = ql;
}
void FEM_DG_IMPLICIT::leftEigenVector(double **dst4, double c, double GAM, double u, double nx, double v, double ny)
{
	double qn = u*nx + v*ny;
	double q2 = u*u + v*v;
	double lx = -ny;
	double ly = nx;
	double ql = u*lx + v*ly;
	double g1 = GAM - 1.0;
	double c2 = c*c;
	dst4[0][0] = 0.5*(0.5*q2*g1 / c2 + qn / c);	dst4[0][1] = -0.5*(g1*u / c2 + nx / c);	dst4[0][2] = -0.5*(g1*v / c2 + ny / c);	dst4[0][3] = 0.5*g1 / c2;
	dst4[1][0] = 1 - 0.5*q2*g1 / c2;			dst4[1][1] = g1*u / c2;				dst4[1][2] = g1*v / c2;				dst4[1][3] = -g1 / c2;
	dst4[2][0] = 0.5*(0.5*q2*g1 / c2 - qn / c);	dst4[2][1] = -0.5*(g1*u / c2 - nx / c);	dst4[2][2] = -0.5*(g1*v / c2 - ny / c);	dst4[2][3] = 0.5*g1 / c2;
	dst4[3][0] = -ql;						dst4[3][1] = lx;					dst4[3][2] = ly;					dst4[3][3] = 0;
}

void FEM_DG_IMPLICIT::calcAP(double **dst4, double **rightEgnVecl4, double **egnVal4, double **leftEgnVecl4)
{
	double		**tempMtx4 = allocMtx4();

	//egnVal4 +.
	for (int i = 0; i < 4; ++i)
	for (int j = 0; j < 4; ++j)
		dst4[i][j] = egnVal4[i][j] > 0 ? egnVal4[i][j] : 0;

	multMtx4(tempMtx4, rightEgnVecl4, dst4);
	multMtx4(dst4, tempMtx4, leftEgnVecl4);
	//multMtxToVal(dst4, 0.5, 4);
	freeMtx4(tempMtx4);
}
void FEM_DG_IMPLICIT::calcAM(double **dst4, double **rightEgnVecl4, double **egnVal4, double **leftEgnVecl4)
{
	double		**tempMtx4 = allocMtx4();

	//egnVal4 -.
	for (int i = 0; i < 4; ++i)
	for (int j = 0; j < 4; ++j)
		dst4[i][j] = egnVal4[i][j] < 0 ? egnVal4[i][j] : 0;

	multMtx4(tempMtx4, rightEgnVecl4, dst4);
	multMtx4(dst4, tempMtx4, leftEgnVecl4);
	//multMtxToVal(dst4, 0.5, 4);
	freeMtx4(tempMtx4);
}

void FEM_DG_IMPLICIT::_calcA(double **dst4, double **rightEgnVecl4, double **egnVal4, double **leftEgnVecl4)
{
	double		**tempMtx4 = allocMtx4();

	multMtx4(tempMtx4, rightEgnVecl4, egnVal4);
	multMtx4(dst4, tempMtx4, leftEgnVecl4);
	freeMtx4(tempMtx4);
}

void FEM_DG_IMPLICIT::calcA(double **dst4, double c, double GAM, double u, double nx, double v, double ny, double H)
{
	double **rightEgnVecl4 = allocMtx4();
	double **egnVal4 = allocMtx4();
	double **leftEgnVecl4 = allocMtx4();

	eigenValues(egnVal4, c, u, nx, v, ny);
	rightEigenVector(rightEgnVecl4, c, u, nx, v, ny, H);
	leftEigenVector(leftEgnVecl4, c, GAM, u, nx, v, ny);
	_calcA(dst4, rightEgnVecl4, egnVal4, leftEgnVecl4);


	freeMtx4(rightEgnVecl4);
	freeMtx4(egnVal4);
	freeMtx4(leftEgnVecl4);
}

void FEM_DG_IMPLICIT::calcAx(double **dst4, double c, double GAM, double u, double v, double H)
{
	calcA(dst4, c, GAM, u, 1.0, v, 0.0, H);
}

void FEM_DG_IMPLICIT::calcAy(double **dst4, double c, double GAM, double u, double v, double H)
{
	calcA(dst4, c, GAM, u, 0.0, v, 1.0, H);
}

void FEM_DG_IMPLICIT::calcAx_(double **dst4, Param par, double GAM)
{
	double AGAM = GAM - 1.0;
	dst4[0][0] = 0.0;
	dst4[0][1] = 1.0;
	dst4[0][2] = 0.0;
	dst4[0][3] = 0.0;

	dst4[1][0] = -POW_2(par.u) + 0.5*par.U2()*AGAM;
	dst4[1][1] = 2.0*par.u - par.u*AGAM;
	dst4[1][2] = -par.v*AGAM;
	dst4[1][3] = AGAM;

	dst4[2][0] = -par.u*par.v;
	dst4[2][1] = par.v;
	dst4[2][2] = par.u;
	dst4[2][3] = 0.0;

	dst4[3][0] = par.U2()*par.u*AGAM - par.E*par.u*GAM;
	dst4[3][1] = -POW_2(par.u)*AGAM + par.E*GAM - 0.5*par.U2()*AGAM;
	dst4[3][2] = par.u*par.v*AGAM;
	dst4[3][3] = par.u*GAM;
}

void FEM_DG_IMPLICIT::calcAy_(double **dst4, Param par, double GAM)
{
	double AGAM = GAM - 1.0;
	dst4[0][0] = 0.0;
	dst4[0][1] = 0.0;
	dst4[0][2] = 1.0;
	dst4[0][3] = 0.0;

	dst4[1][0] = -par.u*par.v;
	dst4[1][1] = par.v;
	dst4[1][2] = par.u;
	dst4[1][3] = 0.0;

	dst4[2][0] = -POW_2(par.u) + (par.r*par.E + 0.5*par.U2())*AGAM;
	dst4[2][1] = 2.0*par.u + (par.r*par.E - par.u)*AGAM;
	dst4[2][2] = -par.v*AGAM;
	dst4[2][3] = AGAM;

	dst4[3][0] = par.U2()*par.u*AGAM - par.E*par.u*GAM;
	dst4[3][1] = -POW_2(par.u)*AGAM + par.E*GAM - 0.5*par.U2()*AGAM;
	dst4[3][2] = par.u*par.v*AGAM;
	dst4[3][3] = par.u*GAM;
}

void FEM_DG_IMPLICIT::calcRoeAverage(Param& average, Param pL, Param pR, double GAM, Vector n)
{
	double WI, UN, UT;
	//double unl = pL.u*n.x+pL.v*n.y;
	//double unr = pR.u*n.x+pR.v*n.y;
	//double utl = pL.u*n.y-pL.v*n.x;
	//double utr = pR.u*n.y-pR.v*n.x;
	//rim_orig(average.r, average.e, average.p, UN, UT, WI, pL.r, pL.p, unl, utl, 0, pR.r, pR.p, unr, utr, 0, GAM);
	//	
	//double UI = UN*n.x+UT*n.y;
	//double VI = UN*n.y-UT*n.x;

	//average.u = UI;
	//average.v = VI;

	roe_orig(average.r, average.e, average.p, average.u, average.v, WI,
		pL.r, pL.p, pL.u, pL.v, 0.0,
		pR.r, pR.p, pR.u, pR.v, 0.0, GAM);

	average.cz = sqrt(GAM*average.p / average.r);
	average.E = average.e + 0.5*(average.u*average.u + average.v*average.v);
}


void FEM_DG_IMPLICIT::consToPar(double fRO, double fRU, double fRV, double fRE, Param& par)
{
	par.r = fRO;
	par.u = fRU / fRO;
	par.v = fRV / fRO;
	par.E = fRE / fRO;
	par.e = par.E - par.U2()*0.5;
}

void FEM_DG_IMPLICIT::calcMatrWithTau()
{
	for (int iCell = 0; iCell < grid.cCount; iCell++) {

		fillMtx(matrBig, 0.0, MATR_DIM);
		

		for (int i = 0; i < BASE_FUNC_COUNT; i++) {
			for (int j = 0; j < BASE_FUNC_COUNT; j++) {
				matrSmall[i][j] = matrA[iCell][i][j] / cTau[iCell];
			}
		}

		for (int ii = 0; ii < FIELD_COUNT; ii++) {
			addSmallMatrToBigMatr(matrBig, matrSmall, ii, ii);
		}
		//for (int ii = 0; ii < FIELD_COUNT; ii++) {
		//	addSmallMatrToBigMatr(matrBig, matrA[iCell], ii, ii);
		//}

		solverMtx->addMatrElement(iCell, iCell, matrBig);

	}
}

void FEM_DG_IMPLICIT::calcIntegral()
{
	double fRO, fRU, fRV, fRE;
	double FR, FU, FV, FE;
	double **mx = allocMtx7();
	double **my = allocMtx7();

	for (int iCell = 0; iCell < grid.cCount; iCell++) {
		
		fillMtx(matrBig, 0.0, MATR_DIM);
		
		for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
			Point& p = cellGP[iCell][iGP];
			double w = cellGW[iCell][iGP];
			getFields(fRO, fRU, fRV, fRE, iCell, p);
			Param par;
			consToPar(fRO, fRU, fRV, fRE, par);
			Material& mat = getMaterial(iCell);
			mat.URS(par, 0); // p=p(r,e)
			double H = par.E + par.p / par.r;
			calcA(mx, par.cz, getGAM(iCell), par.u, 1.0, par.v, 0.0, H);
			calcA(my, par.cz, getGAM(iCell), par.u, 0.0, par.v, 1.0, H);
			//calcAx_(my, par, getGAM(iCell));
			for (int i = 0; i < FIELD_COUNT; i++) {
				for (int j = 0; j < FIELD_COUNT; j++) {
					for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
						for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
							matrSmall[ii][jj] = -mx[i][j] * getDfDx(ii, iCell, p);
							matrSmall[ii][jj] -= my[i][j] * getDfDy(ii, iCell, p);
							matrSmall[ii][jj] *= getF(jj, iCell, p);
							matrSmall[ii][jj] *= w;
						}
					}
					addSmallMatrToBigMatr(matrBig, matrSmall, i, j);
				}
			}

		}
		
		multMtxToVal(matrBig, cellJ[iCell] * SIGMA, MATR_DIM);
		
		solverMtx->addMatrElement(iCell, iCell, matrBig);
	}

	freeMtx7(mx);
	freeMtx7(my);

}

void FEM_DG_IMPLICIT::calcMatrFlux()
{
	double **eigenMtx4, **rEigenVector4, **lEigenVector4;
	double **Amtx4P, **Amtx4M, **mtx4;
	
	mtx4 = allocMtx4();
	eigenMtx4 = allocMtx4();
	rEigenVector4 = allocMtx4();
	lEigenVector4 = allocMtx4();
	Amtx4P = allocMtx4();
	Amtx4M = allocMtx4();


	double fRO1, fRU1, fRV1, fRE1, fRO2, fRU2, fRV2, fRE2;
	Param par1, par2;

	//int mSize = BASE_FUNC_COUNT * 4;

	for (int iEdge = 0; iEdge < grid.eCount; iEdge++) {
		Edge& edge = grid.edges[iEdge];
		Vector&	n = grid.edges[iEdge].n;
		int c1 = edge.c1;
		int c2 = edge.c2;
		
		if (c2 >= 0) {
			
			fillMtx(matrBig, 0.0, MATR_DIM); 
			fillMtx(matrBig2, 0.0, MATR_DIM); 
			
			for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
				Point& p = edgeGP[iEdge][iGP];
				double w = edgeGW[iEdge][iGP];
				getFields(fRO1, fRU1, fRV1, fRE1, c1, p);
				consToPar(fRO1, fRU1, fRV1, fRE1, par1);
				Material& mat1 = getMaterial(c1);
				mat1.URS(par1, 0); // p=p(r,e)
				
				getFields(fRO2, fRU2, fRV2, fRE2, c2, p); 
				consToPar(fRO2, fRU2, fRV2, fRE2, par2);
				Material& mat2 = getMaterial(c2);
				mat2.URS(par2, 0); // p=p(r,e)

				Param average;
				calcRoeAverage(average, par1, par2, getGAM(c1), n);


				double H = average.E + average.p / average.r;
				eigenValues(eigenMtx4, average.cz, average.u, n.x, average.v, n.y);
				rightEigenVector(rEigenVector4, average.cz, average.u, n.x, average.v, n.y, H);
				leftEigenVector(lEigenVector4, average.cz, getGAM(c1), average.u, n.x, average.v, n.y);
				calcAP(Amtx4P, rEigenVector4, eigenMtx4, lEigenVector4);
				calcAM(Amtx4M, rEigenVector4, eigenMtx4, lEigenVector4);
				for (int i = 0; i < FIELD_COUNT; i++) {
					for (int j = 0; j < FIELD_COUNT; j++) {
						for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
							for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
								matrSmall[ii][jj]  = Amtx4P[i][j] * getF(jj, c1, p) * getF(ii, c1, p) * w;
								matrSmall2[ii][jj] = Amtx4M[i][j] * getF(jj, c2, p) * getF(ii, c1, p) * w;
							}
						}
						addSmallMatrToBigMatr(matrBig,  matrSmall, i, j);
						addSmallMatrToBigMatr(matrBig2, matrSmall2, i, j);
					}
				}

			}
			
			multMtxToVal(matrBig, edgeJ[iEdge] * SIGMA, MATR_DIM);
			multMtxToVal(matrBig2, edgeJ[iEdge] * SIGMA, MATR_DIM);

			solverMtx->addMatrElement(c1, c1, matrBig);
			solverMtx->addMatrElement(c1, c2, matrBig2);


			
			fillMtx(matrBig, 0.0, MATR_DIM);
			fillMtx(matrBig2, 0.0, MATR_DIM);

			for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
				Point& p = edgeGP[iEdge][iGP];
				double w = edgeGW[iEdge][iGP];
				getFields(fRO1, fRU1, fRV1, fRE1, c1, p);
				consToPar(fRO1, fRU1, fRV1, fRE1, par1);
				Material& mat1 = getMaterial(c1);
				mat1.URS(par1, 0); // p=p(r,e)

				getFields(fRO2, fRU2, fRV2, fRE2, c2, p);
				consToPar(fRO2, fRU2, fRV2, fRE2, par2);
				Material& mat2 = getMaterial(c2);
				mat2.URS(par2, 0); // p=p(r,e)

				Param average;
				calcRoeAverage(average, par1, par2, getGAM(c1), n);
				double H = average.E + average.p / average.r;

				eigenValues(eigenMtx4, average.cz, average.u, -n.x, average.v, -n.y);
				rightEigenVector(rEigenVector4, average.cz, average.u, -n.x, average.v, -n.y, H);
				leftEigenVector(lEigenVector4, average.cz, getGAM(c1), average.u, -n.x, average.v, -n.y);
				calcAP(Amtx4P, rEigenVector4, eigenMtx4, lEigenVector4);
				calcAM(Amtx4M, rEigenVector4, eigenMtx4, lEigenVector4);
				for (int i = 0; i < FIELD_COUNT; i++) {
					for (int j = 0; j < FIELD_COUNT; j++) {
						for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
							for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
								matrSmall[ii][jj]  = Amtx4P[i][j] * getF(jj, c2, p) * getF(ii, c2, p) * w;
								matrSmall2[ii][jj] = Amtx4M[i][j] * getF(jj, c1, p) * getF(ii, c2, p) * w;
							}
						}
						addSmallMatrToBigMatr(matrBig, matrSmall, i, j);
						addSmallMatrToBigMatr(matrBig2, matrSmall2, i, j);
					}
				}

			}

			multMtxToVal(matrBig, edgeJ[iEdge] * SIGMA, MATR_DIM);
			multMtxToVal(matrBig2, edgeJ[iEdge] * SIGMA, MATR_DIM);

			solverMtx->addMatrElement(c2, c2, matrBig);
			solverMtx->addMatrElement(c2, c1, matrBig2);

		}
		else {

			fillMtx(matrBig, 0.0, MATR_DIM);

			for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
				Point& p = edgeGP[iEdge][iGP];
				double w = edgeGW[iEdge][iGP];
				getFields(fRO1, fRU1, fRV1, fRE1, c1, p);
				consToPar(fRO1, fRU1, fRV1, fRE1, par1);
				Material& mat1 = getMaterial(c1);
                mat1.URS(par1, 0); // p=p(r,e)
                mat1.URS(par1, 1); // T=T(e)

				boundaryCond(iEdge, par1, par2);

				Param average;
				calcRoeAverage(average, par1, par2, getGAM(c1), n);

				double H = average.E + average.p / average.r;

				eigenValues(eigenMtx4, average.cz, average.u, n.x, average.v, n.y);
				rightEigenVector(rEigenVector4, average.cz, average.u, n.x, average.v, n.y, H);
				leftEigenVector(lEigenVector4, average.cz, getGAM(c1), average.u, n.x, average.v, n.y);
				calcAP(Amtx4P, rEigenVector4, eigenMtx4, lEigenVector4);
				for (int i = 0; i < FIELD_COUNT; i++) {
					for (int j = 0; j < FIELD_COUNT; j++) {
						for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
							for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
								matrSmall[ii][jj] = Amtx4P[i][j] * getF(jj, c1, p) * getF(ii, c1, p) * w;
							}
						}
						addSmallMatrToBigMatr(matrBig, matrSmall, i, j);
					}
				}

			}
			multMtxToVal(matrBig, edgeJ[iEdge] * SIGMA, MATR_DIM);

			solverMtx->addMatrElement(c1, c1, matrBig);

		}
	}
	freeMtx4(eigenMtx4);
	freeMtx4(rEigenVector4);
	freeMtx4(lEigenVector4);
	freeMtx4(Amtx4P);
	freeMtx4(Amtx4M);
	freeMtx4(mtx4);
}

void FEM_DG_IMPLICIT::calcRHS()
{
	/* volume integral */

	//const int arrSize = BASE_FUNC_COUNT * 4;


	for (int iCell = 0; iCell < grid.cCount; iCell++) {
		memset(tmpArr, 0, sizeof(double)*MATR_DIM);
		for (int iBF = 0; iBF < BASE_FUNC_COUNT; iBF++) {
			double sRO = 0.0;
			double sRU = 0.0;
			double sRV = 0.0;
			double sRE = 0.0;
			for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
				Point& p = cellGP[iCell][iGP];
				double w = cellGW[iCell][iGP];
				double fRO, fRU, fRV, fRE;
				getFields(fRO, fRU, fRV, fRE, iCell, p);
				Param par;
				consToPar(fRO, fRU, fRV, fRE, par);
				Material& mat = getMaterial(iCell);
				mat.URS(par, 0); // p=p(r,e)

				double FR = par.r*par.u;
				double FU = FR*par.u + par.p;
				double FV = FR*par.v;
				double FE = (fRE + par.p)*par.u;

				double GR = par.r*par.v;
				double GU = GR*par.u;
				double GV = GR*par.v + par.p;
				double GE = (fRE + par.p)*par.v;

				double dFdx = getDfDx(iBF, iCell, p) * w;
				double dFdy = getDfDy(iBF, iCell, p) * w;

				sRO += (FR*dFdx + GR*dFdy);
				sRU += (FU*dFdx + GU*dFdy);
				sRV += (FV*dFdx + GV*dFdy);
				sRE += (FE*dFdx + GE*dFdy);
			}
			sRO *= cellJ[iCell];
			sRU *= cellJ[iCell];
			sRV *= cellJ[iCell];
			sRE *= cellJ[iCell];
			
			int shift = 0;
			tmpArr[shift + iBF] = sRO; shift += BASE_FUNC_COUNT;
			tmpArr[shift + iBF] = sRU; shift += BASE_FUNC_COUNT;
			tmpArr[shift + iBF] = sRV; shift += BASE_FUNC_COUNT;
			tmpArr[shift + iBF] = sRE; //shift += BASE_FUNC_COUNT;
		}

		solverMtx->addRightElement(iCell, tmpArr);
	}

	/* surf integral */

	memset(tmpCFL, 0, grid.cCount*sizeof(double));
	
	for (int iEdge = 0; iEdge < grid.eCount; iEdge++) {
		memset(tmpArr1, 0, sizeof(double)*MATR_DIM);
		memset(tmpArr2, 0, sizeof(double)*MATR_DIM);

		int c1 = grid.edges[iEdge].c1;
		int c2 = grid.edges[iEdge].c2;
		if (c2 >= 0) {
			for (int iBF = 0; iBF < BASE_FUNC_COUNT; iBF++) {
				double sRO1 = 0.0;
				double sRU1 = 0.0;
				double sRV1 = 0.0;
				double sRE1 = 0.0;
				double sRO2 = 0.0;
				double sRU2 = 0.0;
				double sRV2 = 0.0;
				double sRE2 = 0.0;
				for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
					double fRO, fRU, fRV, fRE;
					double FR, FU, FV, FE;

					Point& p = edgeGP[iEdge][iGP];
					double w = edgeGW[iEdge][iGP];

					getFields(fRO, fRU, fRV, fRE, c1, p);
					Param par1;
					consToPar(fRO, fRU, fRV, fRE, par1);
					Material& mat1 = getMaterial(c1);
					mat1.URS(par1, 0); // p=p(r,e)

					getFields(fRO, fRU, fRV, fRE, c2, p);
					Param par2;
					consToPar(fRO, fRU, fRV, fRE, par2);
					Material& mat2 = getMaterial(c2);
					mat2.URS(par2, 0); // p=p(r,e)

					calcFlux(FR, FU, FV, FE, par1, par2, grid.edges[iEdge].n, getGAM(c1));
					
					// вычисляем спектральный радиус для вычисления шага по времени
					if (STEADY) {
						double u1 = sqrt(par1.U2()) + par1.cz;
						double u2 = sqrt(par2.U2()) + par2.cz;
						double lambda = _max_(u1, u2);
						lambda *= w;
						lambda *= edgeJ[iEdge];
						tmpCFL[c1] += lambda;
						tmpCFL[c2] += lambda;
					}

					double cGP1 = w * getF(iBF, c1, p);
					double cGP2 = w * getF(iBF, c2, p);

					sRO1 += FR*cGP1;
					sRU1 += FU*cGP1;
					sRV1 += FV*cGP1;
					sRE1 += FE*cGP1;

					sRO2 += FR*cGP2;
					sRU2 += FU*cGP2;
					sRV2 += FV*cGP2;
					sRE2 += FE*cGP2;
				}

				sRO1 *= edgeJ[iEdge];
				sRU1 *= edgeJ[iEdge];
				sRV1 *= edgeJ[iEdge];
				sRE1 *= edgeJ[iEdge];

				int shift = 0;
				tmpArr1[shift + iBF] = -sRO1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = -sRU1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = -sRV1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = -sRE1; //shift += BASE_FUNC_COUNT;

				sRO2 *= edgeJ[iEdge];
				sRU2 *= edgeJ[iEdge];
				sRV2 *= edgeJ[iEdge];
				sRE2 *= edgeJ[iEdge];

				shift = 0;
				tmpArr2[shift + iBF] = sRO2; shift += BASE_FUNC_COUNT;
				tmpArr2[shift + iBF] = sRU2; shift += BASE_FUNC_COUNT;
				tmpArr2[shift + iBF] = sRV2; shift += BASE_FUNC_COUNT;
				tmpArr2[shift + iBF] = sRE2; //shift += BASE_FUNC_COUNT;
			}
			
			solverMtx->addRightElement(c1, tmpArr1);
			solverMtx->addRightElement(c2, tmpArr2);
		}
		else {
			for (int iBF = 0; iBF < BASE_FUNC_COUNT; iBF++) {
				double sRO1 = 0.0;
				double sRU1 = 0.0;
				double sRV1 = 0.0;
				double sRE1 = 0.0;
				for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
					double fRO, fRU, fRV, fRE;
					double FR, FU, FV, FE;

					Point& p = edgeGP[iEdge][iGP];
					double w = edgeGW[iEdge][iGP];

					getFields(fRO, fRU, fRV, fRE, c1, p);
					Param par1;
					consToPar(fRO, fRU, fRV, fRE, par1);
					Material& mat = getMaterial(c1);
                    mat.URS(par1, 0); // p=p(r,e)
                    mat.URS(par1, 1);

					Param par2;
					boundaryCond(iEdge, par1, par2);

					calcFlux(FR, FU, FV, FE, par1, par2, grid.edges[iEdge].n, getGAM(c1));

					// вычисляем спектральный радиус для вычисления шага по времени
					if (STEADY) {
						double u1 = sqrt(par1.u*par1.u + par1.v*par1.v) + par1.cz;
						double u2 = sqrt(par2.u*par2.u + par2.v*par2.v) + par2.cz;
						double lambda = _max_(u1, u2);
						lambda *= w;
						lambda *= edgeJ[iEdge];
						tmpCFL[c1] += lambda;
					}

					double cGP1 = w * getF(iBF, c1, p);

					sRO1 += FR*cGP1;
					sRU1 += FU*cGP1;
					sRV1 += FV*cGP1;
					sRE1 += FE*cGP1;

				}

				sRO1 *= edgeJ[iEdge];
				sRU1 *= edgeJ[iEdge];
				sRV1 *= edgeJ[iEdge];
				sRE1 *= edgeJ[iEdge];

				int shift = 0;
				tmpArr1[shift + iBF] = -sRO1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = -sRU1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = -sRV1; shift += BASE_FUNC_COUNT;
				tmpArr1[shift + iBF] = -sRE1; //shift += BASE_FUNC_COUNT;

			}

			solverMtx->addRightElement(c1, tmpArr1);
		}
	}
}


void FEM_DG_IMPLICIT::run()
{
	//double __GAM = 1.4;
	
	//// инициализируем портрет матрицы
	//log("Matrix structure initialization:\n");
	//CSRMatrix::DELTA = 65536;
	//for (int iEdge = 0; iEdge < grid.eCount; iEdge++) {
	//	int		c1 = grid.edges[iEdge].c1;
	//	int		c2 = grid.edges[iEdge].c2;
	//	solverMtx->createMatrElement(c1, c1);
	//	if (c2 >= 0){
	//		solverMtx->createMatrElement(c1, c2);
	//		solverMtx->createMatrElement(c2, c2);
	//		solverMtx->createMatrElement(c2, c1);
	//	}
	//	if (iEdge % 100 == 0) {
	//		log("\tfor edge: %d\n", iEdge);
	//	}
	//}
	//solverMtx->initCSR();
	//log("\tcomplete...\n");

	int solverErr = 0;
	double t = 0.0;
	int step = 0;
	long totalCalcTime = 0;
	while (t < TMAX && step < STEP_MAX) {
		long timeStart, timeEnd;
		timeStart = clock();

		if (!solverErr) step++;
		
		if (STEADY) {
			calcTimeStep();
		}
		else {
			t += TAU;
		}

		long tmStart = clock();
		solverErr = 0;
		solverMtx->zero();

		/* Заполняем правую часть */
		calcRHS();
        /* Заполняем правую часть от диффузионных членов уравнений */
		calcDiffusionRHS();

		/* Вычисляем шаги по времени в ячейках по насчитанным ранее значениям спектра */
		if (STEADY) {
			double minTau = DBL_MAX;
			for (int iCell = 0; iCell < grid.cCount; iCell++)
			{
				//Param p;
				//convertConsToPar(iCell, p);
				cTau[iCell] = CFL*grid.cells[iCell].S / (tmpCFL[iCell] + 1.0e-100);
				if (cTau[iCell] < minTau) minTau = cTau[iCell];
			}
			//log("MIN_TAU = %25.15e\n", minTau); 
		}

		/* Заполняем элементы матрицы */
		calcMatrWithTau();		// вычисляем матрицы перед производной по времени
		calcIntegral();			// вычисляем интеграл от(dF / dU)*deltaU*dFi / dx
		calcMatrFlux();			// Вычисляем потоковые величины

		/* Заполняем элементы матрицы от диффузионных членов уравнений */
        calcMatrTensor();			//!< Вычисляем матрицы перед компонентами тензора вязких напряжений
        calcDiffusionIntegral(); 	//!< Вычисляем интеграл от (dH / dU)*dFi / dx и (dG / dU)*dFi / dx
        calcMatrDiffusionFlux();	//!< Вычисляем потоковые величины от диффузионных членов
        //calcTensorIntegral();		//!< Вычисляем интеграл от (dG / dU)*dFi / dx
        //calcMatrTensorFlux();          //!< Вычисляем потоковые величины от градиента полей

		/* Решаем СЛАУ */
		int maxIter = SOLVER_ITER;
		double eps = SOLVER_EPS;

		solverErr = solverMtx->solve(eps, maxIter);

		if (solverErr == MatrixSolver::RESULT_OK) {

			//if (SMOOTHING) smoothingDelta(solverMtx->x);

			for (int cellIndex = 0, ind = 0; cellIndex < grid.cCount; cellIndex++)
			{
				Cell &cell = grid.cells[cellIndex];

				//if (cellIsLim(cellIndex))	continue;
				for (int iFld = 0; iFld < FIELD_COUNT_EXT; iFld++) {
					for (int iF = 0; iF < BASE_FUNC_COUNT; iF++) {
						fields[iFld][cellIndex][iF] += solverMtx->x[ind++];
					}
				}
			}
			
			
			if (limiter != NULL) {
				limiter->run();
			}

			for (int cellIndex = 0, ind = 0; cellIndex < grid.cCount; cellIndex++)
			{
				Cell &cell = grid.cells[cellIndex];

				Param par;
				convertConsToPar(cellIndex, par);
				if (par.r > limitRmax)			{ par.r = limitRmax;	setCellFlagLim(cellIndex); }
				if (par.r < limitRmin)			{ par.r = limitRmin;	setCellFlagLim(cellIndex); }
				if (fabs(par.u) > limitUmax)		{ par.u = limitUmax*par.u / fabs(par.u);	setCellFlagLim(cellIndex); }
				if (fabs(par.v) > limitUmax)		{ par.v = limitUmax*par.v / fabs(par.v);	setCellFlagLim(cellIndex); }
				if (par.p > limitPmax)			{ par.p = limitPmax;	setCellFlagLim(cellIndex); }
				if (par.p < limitPmin)			{ par.p = limitPmin;	setCellFlagLim(cellIndex); }
				if (cellIsLim(cellIndex)) 		{ par.e = par.p / ((getGAM(cellIndex) - 1)*par.r); convertParToCons(cellIndex, par); }

				double fRO, fRU, fRV, fRE;
				for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {

					fRO = getField(FIELD_RO, cellIndex, cellGP[cellIndex][iGP]);
					fRU = getField(FIELD_RU, cellIndex, cellGP[cellIndex][iGP]);
					fRV = getField(FIELD_RV, cellIndex, cellGP[cellIndex][iGP]);
					fRE = getField(FIELD_RE, cellIndex, cellGP[cellIndex][iGP]);
					consToPar(fRO, fRU, fRV, fRE, par);
					if (par.r > limitRmax)			{ par.r = limitRmax;	setCellFlagLim(cellIndex); }
					if (par.r < limitRmin)			{ par.r = limitRmin;	setCellFlagLim(cellIndex); }
					if (fabs(par.u) > limitUmax)		{ par.u = limitUmax*par.u / fabs(par.u);	setCellFlagLim(cellIndex); }
					if (fabs(par.v) > limitUmax)		{ par.v = limitUmax*par.v / fabs(par.v);	setCellFlagLim(cellIndex); }
					if (par.p > limitPmax)			{ par.p = limitPmax;	setCellFlagLim(cellIndex); }
					if (par.p < limitPmin)			{ par.p = limitPmin;	setCellFlagLim(cellIndex); }
					if (cellIsLim(cellIndex)) 		{ par.e = par.p / ((getGAM(cellIndex) - 1)*par.r); convertParToCons(cellIndex, par); }
				}

				for (int iEdge = 0; iEdge < cell.eCount; iEdge++) {
					int edgInd = cell.edgesInd[iEdge];
					for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
						fRO = getField(FIELD_RO, cellIndex, edgeGP[edgInd][iGP]);
						fRU = getField(FIELD_RU, cellIndex, edgeGP[edgInd][iGP]);
						fRV = getField(FIELD_RV, cellIndex, edgeGP[edgInd][iGP]);
						fRE = getField(FIELD_RE, cellIndex, edgeGP[edgInd][iGP]);
						consToPar(fRO, fRU, fRV, fRE, par);
						if (par.r > limitRmax)			{ par.r = limitRmax;	setCellFlagLim(cellIndex); }
						if (par.r < limitRmin)			{ par.r = limitRmin;	setCellFlagLim(cellIndex); }
						if (fabs(par.u) > limitUmax)		{ par.u = limitUmax*par.u / fabs(par.u);	setCellFlagLim(cellIndex); }
						if (fabs(par.v) > limitUmax)		{ par.v = limitUmax*par.v / fabs(par.v);	setCellFlagLim(cellIndex); }
						if (par.p > limitPmax)			{ par.p = limitPmax;	setCellFlagLim(cellIndex); }
						if (par.p < limitPmin)			{ par.p = limitPmin;	setCellFlagLim(cellIndex); }
						if (cellIsLim(cellIndex)) 		{ par.e = par.p / ((getGAM(cellIndex) - 1)*par.r); convertParToCons(cellIndex, par); }
					}
				}


			}
			remediateLimCells();

			int limCells = getLimitedCellsCount();
			if (STEADY && (limCells >= maxLimCells)) decCFL();

			timeEnd = clock();
			totalCalcTime += (timeEnd - timeStart);
			if (step % FILE_SAVE_STEP == 0)
			{
				save(step);
			}
			if (step % PRINT_STEP == 0)
			{
				calcLiftForce();
				if (!STEADY) {

					log("step: %6d  time step: %.16f\tmax iter: %5d\tlim: %4d\tlift force (Fx, Fy) = (%.16f, %.16f)\ttime: %6d ms\ttotal calc time: %ld\n", step, t, maxIter, limCells, Fx, Fy, timeEnd - timeStart, totalCalcTime);
				}
				else {
					log("step: %6d  max iter: %5d\tlim: %4d\tlift force (Fx, Fy) = (%.16f, %.16f)\ttime: %6d ms\ttotal calc time: %ld\n", step, maxIter, limCells, Fx, Fy, timeEnd - timeStart, totalCalcTime);
				}
			}

			if (STEADY && (step % stepCFL == 0)) incCFL();
		}
		else {
			if (solverErr & MatrixSolver::RESULT_ERR_CONVERG) {
				log("Solver error: residual condition not considered.\n");
			}
			if (solverErr & MatrixSolver::RESULT_ERR_MAX_ITER) {
				log("Solver error: max iterations done.\n");
			}
			if (STEADY) {
				decCFL();
			}
			else {
				solverErr = 0;
			}
		}
	}
}

void FEM_DG_IMPLICIT::calcFlux(double& fr, double& fu, double& fv, double& fe, Param pL, Param pR, Vector n, double GAM)
{
	if (FLUX == FLUX_GODUNOV) {	// GODUNOV FLUX
		double RI, EI, PI, UI, VI, WI, UN, UT;
		//double unl = pL.u*n.x+pL.v*n.y;
		//double unr = pR.u*n.x+pR.v*n.y;
		//rim(RI, EI, PI, UN, UI, VI,  pL.r, pL.p, unl, pL.u, pL.v,  pR.r, pR.p, unr, pR.u, pR.v, n, GAM);

		double unl = pL.u*n.x + pL.v*n.y;
		double unr = pR.u*n.x + pR.v*n.y;
		double utl = pL.u*n.y - pL.v*n.x;
		double utr = pR.u*n.y - pR.v*n.x;
		rim_orig(RI, EI, PI, UN, UT, WI, pL.r, pL.p, unl, utl, 0, pR.r, pR.p, unr, utr, 0, GAM);

		UI = UN*n.x + UT*n.y;
		VI = UN*n.y - UT*n.x;

		fr = RI*UN;
		fu = fr*UI + PI*n.x;
		fv = fr*VI + PI*n.y;
		fe = (RI*(EI + 0.5*(UI*UI + VI*VI)) + PI)*UN;
	}
	if (FLUX == FLUX_LAX) {	// LAX-FRIEDRIX FLUX
		double unl = pL.u*n.x + pL.v*n.y;
		double unr = pR.u*n.x + pR.v*n.y;
		double rol, rul, rvl, rel, ror, rur, rvr, rer;
		double alpha = _max_(fabs(unl) + sqrt(GAM*pL.p / pL.r), fabs(unr) + sqrt(GAM*pR.p / pR.r));
		//pL.getToCons(rol, rul, rvl, rel);
		//pR.getToCons(ror, rur, rvr, rer);
		rol = pL.r;
		rul = pL.r*pL.u;
		rvl = pL.r*pL.v;
		rel = pL.p / (GAM - 1.0) + 0.5*pL.r*(pL.u*pL.u + pL.v*pL.v);
		ror = pR.r;
		rur = pR.r*pR.u;
		rvr = pR.r*pR.v;
		rer = pR.p / (GAM - 1.0) + 0.5*pR.r*(pR.u*pR.u + pR.v*pR.v);
		double frl = rol*unl;
		double frr = ror*unr;
		fr = 0.5*(frr + frl - alpha*(ror - rol));
		fu = 0.5*(frr*pR.u + frl*pL.u + (pR.p + pL.p)*n.x - alpha*(rur - rul));
		fv = 0.5*(frr*pR.v + frl*pL.v + (pR.p + pL.p)*n.y - alpha*(rvr - rvl));
		fe = 0.5*((rer + pR.p)*unr + (rel + pL.p)*unl - alpha*(rer - rel));
	}
}

void FEM_DG_IMPLICIT::boundaryCond(int iEdge, Param& pL, Param& pR)
{
    Edge &edge = grid.edges[iEdge];
    int c1 = edge.c1;
    Material& m = getMaterial(c1);
    if (edge.bnd) {
        edge.bnd->run(iEdge, pL, pR);
        m.URS(pR, 2);
        m.URS(pR, 1);
        pR.E = pR.e + 0.5*(pR.U2());
        return;
    }
    else {
        char msg[128];
        sprintf(msg, "Not defined boundary condition for edge %d\n", iEdge);
        throw Exception(msg, Exception::TYPE_BOUND_UNKNOWN);
    }
//	int iBound = -1;
//	for (int i = 0; i < bCount; i++)
//	{
//		if (grid.edges[iEdge].type == boundaries[i].edgeType)
//		{
//			iBound = i;
//			break;
//		}
//	}
//	if (iBound < 0)
//	{
//		log("ERROR (boundary condition): unknown edge type of edge %d...\n", iEdge);
//		EXIT(1);
//	}
//	Boundary& b = boundaries[iBound];
//    CFDBoundary *b = grid.edges[iEdge].bnd;
//	int c1 = grid.edges[iEdge].c1;
//	Material& m = getMaterial(c1);
//	switch (b->name)
//	{
//	    case CFDBoundary::TYPE_INLET:
//		pR.T = b->par[0];		//!< температура
//		pR.p = b->par[1];		//!< давление
//		pR.u = b->par[2];		//!< первая компонента вектора скорости
//		pR.v = b->par[3];		//!< вторая компонента вектора скорости
//
//		m.URS(pR, 2);
//		m.URS(pR, 1);
//		break;
//
//	case Boundary::BOUND_OUTLET:
//		pR = pL;
//		break;
//
//	case Boundary::BOUND_WALL:
//		pR = pL;
//		double Un = pL.u*grid.edges[iEdge].n.x + pL.v*grid.edges[iEdge].n.y;
//		Vector V;
//		V.x = grid.edges[iEdge].n.x*Un*2.0;
//		V.y = grid.edges[iEdge].n.y*Un*2.0;
//		pR.u = pL.u - V.x;
//		pR.v = pL.v - V.y;
//		break;
//	}
}

void FEM_DG_IMPLICIT::addSmallMatrToBigMatr(double **mB, double **mS, int i, int j)
{

	int ii = i * BASE_FUNC_COUNT;
	int jj = j * BASE_FUNC_COUNT;
	for (int i1 = 0; i1 < BASE_FUNC_COUNT; i1++) {
		for (int j1 = 0; j1 < BASE_FUNC_COUNT; j1++) {
			mB[ii + i1][jj + j1] += mS[i1][j1];
		}
	}
}

void FEM_DG_IMPLICIT::done()
{
}

void FEM_DG_IMPLICIT::calcLiftForce()
{
	const double width = 1.0; // предполагаемая ширина профиля по z.
	Param		 par;
	Fx = Fy = 0.0;
	for (int iEdge = 0; iEdge < grid.eCount; ++iEdge)
	{
		if (grid.edges[iEdge].type == Edge::TYPE_WALL)
		{
			int			cellIndex = grid.edges[iEdge].c1;
			double		nx = grid.edges[iEdge].n.x;
			double		ny = grid.edges[iEdge].n.y;
			if (cellIndex < 0) {
				cellIndex = grid.edges[iEdge].c2;
				nx *= -1.0;
				ny *= -1.0;
			}
			convertConsToPar(cellIndex, par);
			Fx += par.p*P_*width*grid.edges[iEdge].l*nx;
			Fy += par.p*P_*width*grid.edges[iEdge].l*ny;
		}
	}
}

double **FEM_DG_IMPLICIT::allocMtx7() {
    double		**tempMtx7 = new double*[7];
    for (int i = 0; i < 7; ++i) tempMtx7[i] = new double[7];
    return tempMtx7;
}

void FEM_DG_IMPLICIT::freeMtx7(double **mtx7) {
    for (int i = 0; i < 7; ++i)
        delete[] mtx7[i];
    delete[] mtx7;
}

void FEM_DG_IMPLICIT::multMtx7(double **dst7, double **srcA7, double **srcB7) {
    double sum;
    for (int i = 0; i < 7; ++i)
    {
        for (int j = 0; j < 7; ++j)
        {
            sum = 0;
            for (int k = 0; k < 7; ++k)
                sum += srcA7[i][k] * srcB7[k][j];
            dst7[i][j] = sum;
        }
    }
}

void FEM_DG_IMPLICIT::clearMtx7(double **mtx7) {
    for (int i = 0; i < 7; ++i)
        for (int j = 0; j < 7; ++j)
            mtx7[i][j] = 0;
}

void FEM_DG_IMPLICIT::getTensorComponents(double &fTAU_XX, double &fTAU_XY, double &fTAU_YY, int iCell, Point p) {
    getTensorComponents(fTAU_XX, fTAU_XY, fTAU_YY, iCell, p.x, p.y);
}

void
FEM_DG_IMPLICIT::getTensorComponents(double &fTAU_XX, double &fTAU_XY, double &fTAU_YY, int iCell, double x, double y) {
    fTAU_XX = getField(FIELD_TAU_XX, iCell, x, y);
    fTAU_XY = getField(FIELD_TAU_XY, iCell, x, y);
    fTAU_YY = getField(FIELD_TAU_YY, iCell, x, y);
}

void FEM_DG_IMPLICIT::calcMatrTensor() {
    for (int iCell = 0; iCell < grid.cCount; iCell++) {

        fillMtx(matrBig, 0.0, MATR_DIM);


        for (int i = 0; i < BASE_FUNC_COUNT; i++) {
            for (int j = 0; j < BASE_FUNC_COUNT; j++) {
                matrSmall[i][j] = matrA[iCell][i][j];
            }
        }

        for (int ii = FIELD_COUNT; ii < FIELD_COUNT_EXT; ii++) {
            addSmallMatrToBigMatr(matrBig, matrSmall, ii, ii);
        }

        solverMtx->addMatrElement(iCell, iCell, matrBig);

    }
}

void FEM_DG_IMPLICIT::calcDiffusionIntegral() {
    double fRO, fRU, fRV, fRE, fTAU_XX, fTAU_XY, fTAU_YY;
    double FR, FU, FV, FE;
    double **mx = allocMtx7();
    double **my = allocMtx7();

    for (int iCell = 0; iCell < grid.cCount; iCell++) {

        fillMtx(matrBig, 0.0, MATR_DIM);

        for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
            Point& p = cellGP[iCell][iGP];
            double w = cellGW[iCell][iGP];
            getFields(fRO, fRU, fRV, fRE, iCell, p);
            Param par;
            consToPar(fRO, fRU, fRV, fRE, par);
            getTensorComponents(fTAU_XX, fTAU_XY, fTAU_YY, iCell, p);
            Material& mat = getMaterial(iCell);
            mat.URS(par, 0); // p=p(r,e)
            mat.getML(par);
            double H = par.E + par.p / par.r;

            calcJ(mx, par.r, par.u, 1.0, par.v, 0.0, par.ML, fTAU_XX, fTAU_XY, fTAU_YY);
            calcJ(my, par.r, par.u, 0.0, par.v, 1.0, par.ML, fTAU_XX, fTAU_XY, fTAU_YY);
            for (int i = 0; i < FIELD_COUNT_EXT; i++) {
                for (int j = 0; j < FIELD_COUNT_EXT; j++) {
                    for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
                        for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
                            matrSmall[ii][jj] = mx[i][j] * getDfDx(ii, iCell, p);
                            matrSmall[ii][jj] += my[i][j] * getDfDy(ii, iCell, p);
                            matrSmall[ii][jj] *= getF(jj, iCell, p); // basis function
                            matrSmall[ii][jj] *= w;
                        }
                    }
                    addSmallMatrToBigMatr(matrBig, matrSmall, i, j);
                }
            }

        }

        multMtxToVal(matrBig, cellJ[iCell] * SIGMA, MATR_DIM);

        solverMtx->addMatrElement(iCell, iCell, matrBig);
    }

    freeMtx7(mx);
    freeMtx7(my);
}

void FEM_DG_IMPLICIT::calcMatrDiffusionFlux() {
    double  **Amtx7P, **Amtx7M;

    Amtx7P = allocMtx7();
    Amtx7M = allocMtx7();


    double fRO1, fRU1, fRV1, fRE1, fRO2, fRU2, fRV2, fRE2, fTAU_XX1, fTAU_XY1, fTAU_YY1, fTAU_XX2, fTAU_XY2, fTAU_YY2;
    Param par1, par2;

    for (int iEdge = 0; iEdge < grid.eCount; iEdge++) {
        Edge& edge = grid.edges[iEdge];
        Vector&	n = grid.edges[iEdge].n;
        int c1 = edge.c1;
        int c2 = edge.c2;

        if (c2 >= 0) {

            fillMtx(matrBig, 0.0, MATR_DIM);
            fillMtx(matrBig2, 0.0, MATR_DIM);

            for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
                Point& p = edgeGP[iEdge][iGP];
                double w = edgeGW[iEdge][iGP];
                getFields(fRO1, fRU1, fRV1, fRE1, c1, p);
                consToPar(fRO1, fRU1, fRV1, fRE1, par1);
                Material& mat1 = getMaterial(c1);
                mat1.URS(par1, 0); // p=p(r,e)
                mat1.getML(par1);
                getTensorComponents(fTAU_XX1, fTAU_XY1, fTAU_YY1, c1, p);

                getFields(fRO2, fRU2, fRV2, fRE2, c2, p);
                consToPar(fRO2, fRU2, fRV2, fRE2, par2);
                Material& mat2 = getMaterial(c2);
                mat2.URS(par2, 0); // p=p(r,e)
                mat2.getML(par2);
                getTensorComponents(fTAU_XX2, fTAU_XY2, fTAU_YY2, c2, p);

                calcJ(Amtx7P, par1.r, par1.u, n.x, par1.v, n.y, par1.ML, fTAU_XX1, fTAU_XY1, fTAU_YY1);
                calcJ(Amtx7M, par2.r, par2.u, n.x, par2.v, n.y, par2.ML, fTAU_XX2, fTAU_XY2, fTAU_YY2);
                for (int i = 0; i < FIELD_COUNT_EXT; i++) {
                    for (int j = 0; j < FIELD_COUNT_EXT; j++) {
                        for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
                            for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
                                matrSmall[ii][jj]  = Amtx7P[i][j] * getF(jj, c1, p) * getF(ii, c1, p) * w;
                                matrSmall2[ii][jj] = Amtx7M[i][j] * getF(jj, c2, p) * getF(ii, c1, p) * w;
                            }
                        }
                        addSmallMatrToBigMatr(matrBig,  matrSmall, i, j);
                        addSmallMatrToBigMatr(matrBig2, matrSmall2, i, j);
                    }
                }

            }

            multMtxToVal(matrBig, edgeJ[iEdge] * SIGMA, MATR_DIM);
            multMtxToVal(matrBig2, edgeJ[iEdge] * SIGMA, MATR_DIM);

            solverMtx->addMatrElement(c1, c1, matrBig);
            solverMtx->addMatrElement(c1, c2, matrBig2);



            fillMtx(matrBig, 0.0, MATR_DIM);
            fillMtx(matrBig2, 0.0, MATR_DIM);

            for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
                Point& p = edgeGP[iEdge][iGP];
                double w = edgeGW[iEdge][iGP];
                getFields(fRO1, fRU1, fRV1, fRE1, c1, p);
                consToPar(fRO1, fRU1, fRV1, fRE1, par1);
                Material& mat1 = getMaterial(c1);
                mat1.URS(par1, 0); // p=p(r,e)
                mat1.getML(par1);
                getTensorComponents(fTAU_XX1, fTAU_XY1, fTAU_YY1, c1, p);

                getFields(fRO2, fRU2, fRV2, fRE2, c2, p);
                consToPar(fRO2, fRU2, fRV2, fRE2, par2);
                Material& mat2 = getMaterial(c2);
                mat2.URS(par2, 0); // p=p(r,e)
                mat2.getML(par2);
                getTensorComponents(fTAU_XX1, fTAU_XY1, fTAU_YY1, c1, p);

                calcJ(Amtx7P, par1.r, par1.u, -n.x, par1.v, -n.y, par1.ML, fTAU_XX1, fTAU_XY1, fTAU_YY1);
                calcJ(Amtx7M, par2.r, par2.u, -n.x, par2.v, -n.y, par2.ML, fTAU_XX2, fTAU_XY2, fTAU_YY2);
                for (int i = 0; i < FIELD_COUNT_EXT; i++) {
                    for (int j = 0; j < FIELD_COUNT_EXT; j++) {
                        for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
                            for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
                                matrSmall[ii][jj]  = Amtx7P[i][j] * getF(jj, c2, p) * getF(ii, c2, p) * w;
                                matrSmall2[ii][jj] = Amtx7M[i][j] * getF(jj, c1, p) * getF(ii, c2, p) * w;
                            }
                        }
                        addSmallMatrToBigMatr(matrBig, matrSmall, i, j);
                        addSmallMatrToBigMatr(matrBig2, matrSmall2, i, j);
                    }
                }

            }

            multMtxToVal(matrBig, edgeJ[iEdge] * SIGMA, MATR_DIM);
            multMtxToVal(matrBig2, edgeJ[iEdge] * SIGMA, MATR_DIM);

            solverMtx->addMatrElement(c2, c2, matrBig);
            solverMtx->addMatrElement(c2, c1, matrBig2);

        }
    }

    freeMtx7(Amtx7P);
    freeMtx7(Amtx7M);
}

void FEM_DG_IMPLICIT::calcTensorIntegral() {
    double fRO, fRU, fRV, fRE, fTAU_XX, fTAU_XY, fTAU_YY;
    double FR, FU, FV, FE;
    double **mx = allocMtx7();
    double **my = allocMtx7();

    for (int iCell = 0; iCell < grid.cCount; iCell++) {

        fillMtx(matrBig, 0.0, MATR_DIM);

        for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
            Point& p = cellGP[iCell][iGP];
            double w = cellGW[iCell][iGP];
            getFields(fRO, fRU, fRV, fRE, iCell, p);
            Param par;
            consToPar(fRO, fRU, fRV, fRE, par);
            getTensorComponents(fTAU_XX, fTAU_XY, fTAU_YY, iCell, p);
            Material& mat = getMaterial(iCell);
            mat.URS(par, 0); // p=p(r,e)
            mat.getML(par);

            calcJ(mx, par.r, par.u, 1.0, par.v, 0.0, par.ML, fTAU_XX, fTAU_XY, fTAU_YY);
            calcJ(my, par.r, par.u, 0.0, par.v, 1.0, par.ML, fTAU_XX, fTAU_XY, fTAU_YY);
            for (int i = FIELD_COUNT; i < FIELD_COUNT_EXT; i++) {
                for (int j = 0; j < FIELD_COUNT_EXT; j++) {
                    for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
                        for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
                            matrSmall[ii][jj] = mx[i][j] * getDfDx(ii, iCell, p);
                            matrSmall[ii][jj] += my[i][j] * getDfDy(ii, iCell, p);
                            matrSmall[ii][jj] *= getF(jj, iCell, p);
                            matrSmall[ii][jj] *= w;
                        }
                    }
                    addSmallMatrToBigMatr(matrBig, matrSmall, i, j);
                }
            }

        }

        multMtxToVal(matrBig, cellJ[iCell] * SIGMA, MATR_DIM);

        solverMtx->addMatrElement(iCell, iCell, matrBig);
    }

    freeMtx7(mx);
    freeMtx7(my);
}

void FEM_DG_IMPLICIT::calcMatrTensorFlux() {
    double  **Amtx7P, **Amtx7M;//, **mtx7;

    //mtx7 = allocMtx7();
    Amtx7P = allocMtx7();
    Amtx7M = allocMtx7();


    double fRO1, fRU1, fRV1, fRE1, fRO2, fRU2, fRV2, fRE2, fTAU_XX1, fTAU_XY1, fTAU_YY1, fTAU_XX2, fTAU_XY2, fTAU_YY2;
    Param par1, par2;

    for (int iEdge = 0; iEdge < grid.eCount; iEdge++) {
        Edge& edge = grid.edges[iEdge];
        Vector&	n = grid.edges[iEdge].n;
        int c1 = edge.c1;
        int c2 = edge.c2;

        if (c2 >= 0) {

            fillMtx(matrBig, 0.0, MATR_DIM);
            fillMtx(matrBig2, 0.0, MATR_DIM);

            for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
                Point& p = edgeGP[iEdge][iGP];
                double w = edgeGW[iEdge][iGP];
                getFields(fRO1, fRU1, fRV1, fRE1, c1, p);
                consToPar(fRO1, fRU1, fRV1, fRE1, par1);
                Material& mat1 = getMaterial(c1);
                mat1.URS(par1, 0); // p=p(r,e)
                mat1.getML(par1);
                getTensorComponents(fTAU_XX1, fTAU_XY1, fTAU_YY1, c1, p);

                getFields(fRO2, fRU2, fRV2, fRE2, c2, p);
                consToPar(fRO2, fRU2, fRV2, fRE2, par2);
                Material& mat2 = getMaterial(c2);
                mat2.URS(par2, 0); // p=p(r,e)
                mat2.getML(par2);
                getTensorComponents(fTAU_XX2, fTAU_XY2, fTAU_YY2, c2, p);

                calcJ(Amtx7P, par1.r, par1.u, n.x, par1.v, n.y, par1.ML, fTAU_XX1, fTAU_XY1, fTAU_YY1);
                calcJ(Amtx7M, par2.r, par2.u, n.x, par2.v, n.y, par2.ML, fTAU_XX2, fTAU_XY2, fTAU_YY2);
                for (int i = FIELD_COUNT; i < FIELD_COUNT_EXT; i++) {
                    for (int j = 0; j < FIELD_COUNT_EXT; j++) {
                        for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
                            for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
                                matrSmall[ii][jj]  = 0.5 * Amtx7P[i][j] * getF(jj, c1, p) * getF(ii, c1, p) * w;
                                matrSmall2[ii][jj] = 0.5 * Amtx7M[i][j] * getF(jj, c2, p) * getF(ii, c1, p) * w;
                            }
                        }
                        addSmallMatrToBigMatr(matrBig,  matrSmall, i, j);
                        addSmallMatrToBigMatr(matrBig2, matrSmall2, i, j);
                    }
                }

            }

            multMtxToVal(matrBig, edgeJ[iEdge] * SIGMA, MATR_DIM);
            multMtxToVal(matrBig2, edgeJ[iEdge] * SIGMA, MATR_DIM);

            solverMtx->addMatrElement(c1, c1, matrBig);
            solverMtx->addMatrElement(c1, c2, matrBig2);



            fillMtx(matrBig, 0.0, MATR_DIM);
            fillMtx(matrBig2, 0.0, MATR_DIM);

            for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
                Point& p = edgeGP[iEdge][iGP];
                double w = edgeGW[iEdge][iGP];
                getFields(fRO1, fRU1, fRV1, fRE1, c1, p);
                consToPar(fRO1, fRU1, fRV1, fRE1, par1);
                Material& mat1 = getMaterial(c1);
                mat1.URS(par1, 0); // p=p(r,e)
                mat1.getML(par1);
                getTensorComponents(fTAU_XX1, fTAU_XY1, fTAU_YY1, c1, p);

                getFields(fRO2, fRU2, fRV2, fRE2, c2, p);
                consToPar(fRO2, fRU2, fRV2, fRE2, par2);
                Material& mat2 = getMaterial(c2);
                mat2.URS(par2, 0); // p=p(r,e)
                mat2.getML(par2);
                getTensorComponents(fTAU_XX1, fTAU_XY1, fTAU_YY1, c1, p);

                calcJ(Amtx7P, par1.r, par1.u, -n.x, par1.v, -n.y, par1.ML, fTAU_XX1, fTAU_XY1, fTAU_YY1);
                calcJ(Amtx7M, par2.r, par2.u, -n.x, par2.v, -n.y, par2.ML, fTAU_XX2, fTAU_XY2, fTAU_YY2);
                for (int i = FIELD_COUNT; i < FIELD_COUNT_EXT; i++) {
                    for (int j = 0; j < FIELD_COUNT_EXT; j++) {
                        for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
                            for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
                                matrSmall[ii][jj]  = 0.5 * Amtx7P[i][j] * getF(jj, c2, p) * getF(ii, c2, p) * w;
                                matrSmall2[ii][jj] = 0.5 * Amtx7M[i][j] * getF(jj, c1, p) * getF(ii, c2, p) * w;
                            }
                        }
                        addSmallMatrToBigMatr(matrBig, matrSmall, i, j);
                        addSmallMatrToBigMatr(matrBig2, matrSmall2, i, j);
                    }
                }

            }

            multMtxToVal(matrBig, edgeJ[iEdge] * SIGMA, MATR_DIM);
            multMtxToVal(matrBig2, edgeJ[iEdge] * SIGMA, MATR_DIM);

            solverMtx->addMatrElement(c2, c2, matrBig);
            solverMtx->addMatrElement(c2, c1, matrBig2);

        }
//        else {
//
//            fillMtx(matrBig, 0.0, MATR_DIM);
//
//            for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
//                Point& p = edgeGP[iEdge][iGP];
//                double w = edgeGW[iEdge][iGP];
//                getFields(fRO1, fRU1, fRV1, fRE1, c1, p);
//                consToPar(fRO1, fRU1, fRV1, fRE1, par1);
//                Material& mat1 = getMaterial(c1);
//                mat1.URS(par1, 0); // p=p(r,e)
//                mat1.getML(par1);
//
//                boundaryCond(iEdge, par1, par2);
//
////                Param average;
////                calcRoeAverage(average, par1, par2, getGAM(c1), n);
////
////                double H = average.E + average.p / average.r;
//
////                eigenValues(eigenMtx4, average.cz, average.u, n.x, average.v, n.y);
////                rightEigenVector(rEigenVector4, average.cz, average.u, n.x, average.v, n.y, H);
////                leftEigenVector(lEigenVector4, average.cz, getGAM(c1), average.u, n.x, average.v, n.y);
////                calcAP(Amtx4P, rEigenVector4, eigenMtx4, lEigenVector4);
//                calcJ(Amtx7P, par2.r, par2.u, n.x, par2.v, n.y, par2.ML, fTAU_XX1, fTAU_XY1, fTAU_YY1);
//                //calcJ(Amtx7M, par2.r, par2.u, n.x, par2.v, n.y, par2.ML, fTAU_XX2, fTAU_XY2, fTAU_YY2);
//                for (int i = 0; i < FIELD_COUNT; i++) {
//                    for (int j = 0; j < FIELD_COUNT; j++) {
//                        for (int ii = 0; ii < BASE_FUNC_COUNT; ii++) {
//                            for (int jj = 0; jj < BASE_FUNC_COUNT; jj++) {
//                                matrSmall[ii][jj] = Amtx7P[i][j] * getF(ii, c1, p) * getF(jj, c1, p) * w;
//                            }
//                        }
//                        addSmallMatrToBigMatr(matrBig, matrSmall, i, j);
//                    }
//                }
//
//            }
//            multMtxToVal(matrBig, edgeJ[iEdge] * SIGMA, MATR_DIM);
//
//            solverMtx->addMatrElement(c1, c1, matrBig);
//
//        }
    }

    freeMtx7(Amtx7P);
    freeMtx7(Amtx7M);
//    freeMtx4(mtx7);
}

void FEM_DG_IMPLICIT::calcJ(double **dst7, double r, double u, double nx, double v, double ny, double mu, double tau_xx, double tau_xy,
                            double tau_yy) {
    dst7[0][0] = 0.0; dst7[0][1] = 0.0; dst7[0][2] = 0.0; dst7[0][3] = 0.0; dst7[0][4] = 0.0; dst7[0][5] = 0.0; dst7[0][6] = 0.0;
    dst7[1][0] = 0.0; dst7[1][1] = 0.0; dst7[1][2] = 0.0; dst7[1][3] = 0.0; dst7[1][4] = nx; dst7[1][5] = ny; dst7[1][6] = 0.0;
    dst7[2][0] = 0.0; dst7[2][1] = 0.0; dst7[2][2] = 0.0; dst7[2][3] = 0.0; dst7[2][4] = 0.0; dst7[2][5] = nx; dst7[2][6] = ny;
    dst7[3][0] = -(nx*(tau_xx*u + tau_xy*v) + ny*(tau_xy*u + tau_yy*v)) / r; dst7[3][1] = (tau_xx*nx + tau_xy*ny)/r; dst7[3][2] = (tau_xy*nx + tau_yy*ny) / r; dst7[3][3] = 0.0; dst7[3][4] = u*nx; dst7[3][5] = v*nx + u*ny; dst7[3][6] = v*ny;
    dst7[4][0] = -mu*(4.*u*nx - 2.*v*ny) / 3. / r; dst7[4][1] = 4.*mu*nx / 3. / r; dst7[4][2] = -2.*mu*ny / 3. / r; dst7[4][3] = 0.0; dst7[4][4] = 0.0; dst7[4][5] = 0.0; dst7[4][6] = 0.0;
    dst7[5][0] = -mu*((v - 2.*u / 3.)*nx + (u - 2.*v / 3.)*ny) / r; dst7[5][1] = mu*(-2.*nx / 3. + ny) / r; dst7[5][2] = mu*(nx - 2.*ny / 3.) / r; dst7[5][3] = 0.0; dst7[5][4] = 0.0; dst7[5][5] = 0.0; dst7[5][6] = 0.0;
    dst7[6][0] = -mu*(-2.*u*nx + 4.*v*ny) / 3. / r; dst7[6][1] = -2.*mu*nx / 3. / r; dst7[6][2] = 4.*mu*ny / 3. / r; dst7[6][3] = 0.0; dst7[6][4] = 0.0; dst7[6][5] = 0.0; dst7[6][6] = 0.0;
}

void FEM_DG_IMPLICIT::calcDiffusionRHS() {

    /* volume integral */

    for (int iCell = 0; iCell < grid.cCount; iCell++) {
        memset(tmpArr, 0, sizeof(double)*MATR_DIM);
        for (int iBF = 0; iBF < BASE_FUNC_COUNT; iBF++) {
            double s1 = 0.0;
            double s2 = 0.0;
            double s3 = 0.0;
            double s4 = 0.0;
            double s5 = 0.0;
            double s6 = 0.0;
            double s7 = 0.0;
            for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
                Point& p = cellGP[iCell][iGP];
                double w = cellGW[iCell][iGP];
                double fRO, fRU, fRV, fRE, fTAU_XX, fTAU_XY, fTAU_YY;
                getFields(fRO, fRU, fRV, fRE, iCell, p);
                Param par;
                consToPar(fRO, fRU, fRV, fRE, par);
                getTensorComponents(fTAU_XX, fTAU_XY, fTAU_YY, iCell, p);
                Material& mat = getMaterial(iCell);
                mat.URS(par, 0); // p=p(r,e)
                mat.getML(par);

                double F1 = 0.;
                double F2 = fTAU_XX;
                double F3 = fTAU_XY;
                double F4 = fTAU_XX*par.u + fTAU_XY*par.v;
                double F5 = par.ML*(4.*par.u / 3.);
                double F6 = par.ML*(par.v - 2.*par.u / 3.);
                double F7 = par.ML*(-2.*par.u / 3.);

                double G1 = 0.;
                double G2 = fTAU_XY;
                double G3 = fTAU_YY;
                double G4 = fTAU_XY*par.u + fTAU_YY*par.v;
                double G5 = par.ML*(-2.*par.v / 3.);
                double G6 = par.ML*(par.u - 2.*par.v / 3.);
                double G7 = par.ML*(4.*par.v / 3.);

                double dFdx = getDfDx(iBF, iCell, p) * w;
                double dFdy = getDfDy(iBF, iCell, p) * w;

                s1 += (F1*dFdx + G1*dFdy);
                s2 += (F2*dFdx + G2*dFdy);
                s3 += (F3*dFdx + G3*dFdy);
                s4 += (F4*dFdx + G4*dFdy);
                s5 += (F5*dFdx + G5*dFdy);
                s6 += (F6*dFdx + G6*dFdy);
                s7 += (F7*dFdx + G7*dFdy);
            }
            s1 *= cellJ[iCell];
            s2 *= cellJ[iCell];
            s3 *= cellJ[iCell];
            s4 *= cellJ[iCell];
            s5 *= cellJ[iCell];
            s6 *= cellJ[iCell];
            s7 *= cellJ[iCell];

            int shift = 0;
            tmpArr[shift + iBF] = -s1; shift += BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s2; shift += BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s3; shift += BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s4; shift += BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s5; shift += BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s6; shift += BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s7; //shift += BASE_FUNC_COUNT;
        }

        solverMtx->addRightElement(iCell, tmpArr);
    }
    // W^m contribution
    for (int iCell = 0; iCell < grid.cCount; iCell++) {
        memset(tmpArr, 0, sizeof(double)*MATR_DIM);
        for (int iBF = 0; iBF < BASE_FUNC_COUNT; iBF++) {
            double s5 = 0.0;
            double s6 = 0.0;
            double s7 = 0.0;
            for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
                Point& p = cellGP[iCell][iGP];
                double w = cellGW[iCell][iGP];
                double fTAU_XX, fTAU_XY, fTAU_YY;
                getTensorComponents(fTAU_XX, fTAU_XY, fTAU_YY, iCell, p);

                double cGP = w*getF(iBF, iCell, p);

                s5 += fTAU_XX*cGP;
                s6 += fTAU_XY*cGP;
                s7 += fTAU_YY*cGP;
            }
            s5 *= cellJ[iCell];
            s6 *= cellJ[iCell];
            s7 *= cellJ[iCell];

            int shift = FIELD_COUNT*BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s5; shift += BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s6; shift += BASE_FUNC_COUNT;
            tmpArr[shift + iBF] = -s7; //shift += BASE_FUNC_COUNT;
        }

        solverMtx->addRightElement(iCell, tmpArr);
    }

    /* surf integral */

    memset(tmpCFL, 0, grid.cCount*sizeof(double));

    for (int iEdge = 0; iEdge < grid.eCount; iEdge++) {
        memset(tmpArr1, 0, sizeof(double)*MATR_DIM);
        memset(tmpArr2, 0, sizeof(double)*MATR_DIM);

        int c1 = grid.edges[iEdge].c1;
        int c2 = grid.edges[iEdge].c2;
        Vector n = grid.edges[iEdge].n;
        if (c2 >= 0) {
            for (int iBF = 0; iBF < BASE_FUNC_COUNT; iBF++) {
                //double s11 = 0.0;
                double s21 = 0.0;
                double s31 = 0.0;
                double s41 = 0.0;
                double s51 = 0.0;
                double s61 = 0.0;
                double s71 = 0.0;

                //double s12 = 0.0;
                double s22 = 0.0;
                double s32 = 0.0;
                double s42 = 0.0;
                double s52 = 0.0;
                double s62 = 0.0;
                double s72 = 0.0;
                for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
                    double fRO, fRU, fRV, fRE, fTAU_XX1, fTAU_XY1, fTAU_YY1, fTAU_XX2, fTAU_XY2, fTAU_YY2;
                    double FS2, FS3, FS4, FS5, FS6, FS7;

                    Point& p = edgeGP[iEdge][iGP];
                    double w = edgeGW[iEdge][iGP];

                    getFields(fRO, fRU, fRV, fRE, c1, p);
                    Param par1;
                    consToPar(fRO, fRU, fRV, fRE, par1);
                    Material& mat1 = getMaterial(c1);
                    mat1.URS(par1, 0); // p=p(r,e)
                    mat1.getML(par1);
                    getTensorComponents(fTAU_XX1, fTAU_XY1, fTAU_YY1, c1, p);

                    getFields(fRO, fRU, fRV, fRE, c2, p);
                    Param par2;
                    consToPar(fRO, fRU, fRV, fRE, par2);
                    Material& mat2 = getMaterial(c2);
                    mat2.URS(par2, 0); // p=p(r,e)
                    mat2.getML(par2);
                    getTensorComponents(fTAU_XX2, fTAU_XY2, fTAU_YY2, c2, p);

                    FS2 = 0.5*((fTAU_XX1 + fTAU_XX2)*n.x + (fTAU_XY1 + fTAU_XY2)*n.y);
                    FS3 = 0.5*((fTAU_XY1 + fTAU_XY2)*n.x + (fTAU_YY1 + fTAU_YY2)*n.y);
                    FS4 = 0.5*((fTAU_XX1*par1.u + fTAU_XY1*par1.v + fTAU_XX2*par2.u + fTAU_XY2*par2.v)*n.x + (fTAU_XY1*par1.u + fTAU_YY1*par1.v + fTAU_XY2*par2.u + fTAU_YY2*par2.v)*n.y);
                    FS5 = 0.5*(par1.ML + par2.ML)*(4.*(par1.u + par2.u) / 3.*n.x - 2.*(par1.v + par2.v) / 3.*n.y);
                    FS6 = 0.5*(par1.ML + par2.ML)*((par1.v - 2.*par1.u / 3. + par2.v - 2.*par2.u / 3.)*n.x + (par1.u - 2.*par1.v / 3. + par2.u - 2.*par2.v / 3.)*n.y);
                    FS7 = 0.5*(par1.ML + par2.ML)*(-2.*(par1.u + par2.u) / 3.*n.x + 4.*(par1.v + par2.v) / 3.*n.y);

                    double cGP1 = w * getF(iBF, c1, p);
                    double cGP2 = w * getF(iBF, c2, p);

                    //s11 += 0.;
                    s21 += FS2*cGP1;
                    s31 += FS3*cGP1;
                    s41 += FS4*cGP1;
                    s51 += FS5*cGP1;
                    s61 += FS6*cGP1;
                    s71 += FS7*cGP1;

                    //s12 += 0.;
                    s22 += FS2*cGP2;
                    s32 += FS3*cGP2;
                    s42 += FS4*cGP2;
                    s52 += FS5*cGP2;
                    s62 += FS6*cGP2;
                    s72 += FS7*cGP2;
                }

                //s11 *= edgeJ[iEdge];
                s21 *= edgeJ[iEdge];
                s31 *= edgeJ[iEdge];
                s41 *= edgeJ[iEdge];
                s51 *= edgeJ[iEdge];
                s61 *= edgeJ[iEdge];
                s71 *= edgeJ[iEdge];

                int shift = BASE_FUNC_COUNT;
                tmpArr1[shift + iBF] = s21; shift += BASE_FUNC_COUNT;
                tmpArr1[shift + iBF] = s31; shift += BASE_FUNC_COUNT;
                tmpArr1[shift + iBF] = s41; shift += BASE_FUNC_COUNT;
                tmpArr1[shift + iBF] = s51; shift += BASE_FUNC_COUNT;
                tmpArr1[shift + iBF] = s61; shift += BASE_FUNC_COUNT;
                tmpArr1[shift + iBF] = s71; //shift += BASE_FUNC_COUNT;

                s22 *= edgeJ[iEdge];
                s32 *= edgeJ[iEdge];
                s42 *= edgeJ[iEdge];
                s52 *= edgeJ[iEdge];
                s62 *= edgeJ[iEdge];
                s72 *= edgeJ[iEdge];

                shift = BASE_FUNC_COUNT;
                tmpArr2[shift + iBF] = -s22; shift += BASE_FUNC_COUNT;
                tmpArr2[shift + iBF] = -s32; shift += BASE_FUNC_COUNT;
                tmpArr2[shift + iBF] = -s42; shift += BASE_FUNC_COUNT;
                tmpArr2[shift + iBF] = -s52; shift += BASE_FUNC_COUNT;
                tmpArr2[shift + iBF] = -s62; shift += BASE_FUNC_COUNT;
                tmpArr2[shift + iBF] = -s72; //shift += BASE_FUNC_COUNT;
            }

            solverMtx->addRightElement(c1, tmpArr1);
            solverMtx->addRightElement(c2, tmpArr2);
        }
        else {
            if (strcmp(grid.edges[iEdge].typeName, CFDBoundary::TYPE_WALL_NO_SLIP) == 0){
                for (int iBF = 0; iBF < BASE_FUNC_COUNT; iBF++) {
                    //double sRO1 = 0.0;
                    double sRU1 = 0.0;
                    double sRV1 = 0.0;
                    double sRE1 = 0.0;

                    double fRO, fRU, fRV, fRE;
                    double FR, FU, FV, FE;
                    double fDL = grid.edges[iEdge].l;
                    Point  pC  = grid.edges[iEdge].c[0];

                    getFields(fRO, fRU, fRV, fRE, c1, pC);
                    Param par1;
                    consToPar(fRO, fRU, fRV, fRE, par1);
                    Material& mat = getMaterial(c1);
                    mat.getML(par1);

                    // вектор скорости в ячейке
                    Vector vV( par1.u, par1.v);
                    //Определяем заданную скорость вращения на регионе
                    Vector vOmega(0., 0.);
                    Vector vOmegaR = vOmega; //vector_prod( vOmega, pC ); TODO: в двумерном пространстве векторное произведение вырождается в скаляр
                    Vector vVwall = vOmegaR;
                    vV -= vVwall;

                    // значение скорости по нормали
                    double fVn = scalar_prod( vV, n );

                    // нормальная компонента скорости
                    Vector vVn = n;  vVn *= fVn;

                    // тангенциальная компонента скорости
                    Vector vVt = vV;  vVt -= vVn;

                    double fMU = par1.ML;			// молекулярная вязкость
                    double fKP = 0.0;			// коэффициент теплопроводности TODO: добавить в код

                    Vector vRP  = grid.edges[iEdge].c[0];
                    vRP -= grid.cells[c1].c;			// расстояние от центра ячейки до центра грани
                    double fLP_ = scalar_prod(vRP, n);		// от точки P' до центра грани

                    Vector vRP_ = n;  vRP_ *= - fLP_;  vRP_ += vRP;

                    fLP_ = 1.0 / fLP_;

                    Vector vTauN;

                    vTauN.x = - fMU * vVt.x * fLP_;
                    vTauN.y = - fMU * vVt.y * fLP_;

                    double fQV = scalar_prod( vTauN, vVwall );			// работа вязких сил

                    double fQT = 0.0;		// тепловой поток на границе

                    sRU1 += vTauN.x*fDL;
                    sRV1 += vTauN.y*fDL;
                    sRE1 += (fQV + fQT) * fDL;

                    int shift = BASE_FUNC_COUNT;
                    tmpArr1[shift + iBF] = sRU1; shift += BASE_FUNC_COUNT;
                    tmpArr1[shift + iBF] = sRV1; shift += BASE_FUNC_COUNT;
                    tmpArr1[shift + iBF] = sRE1; //shift += BASE_FUNC_COUNT;

                }
            }

            solverMtx->addRightElement(c1, tmpArr1);
        }
    }
}

