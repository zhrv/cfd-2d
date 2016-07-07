#include "heat_dg_implicit.h"
#include "tinyxml.h"
#include "global.h"
#include <ctime>
#include <cfloat>
#include "LimiterDG.h"
#include "MatrixSolver.h"
#include "bnd_cond_heat.h"
#include "MeshReader.h"

#define POW_2(x) ((x)*(x))


void HEAT_DG_IMPLICIT::init(char * xmlFileName)
{
    TiXmlDocument doc(xmlFileName);
    bool loadOkay = doc.LoadFile(TIXML_ENCODING_UTF8);
    if (!loadOkay)
    {
        log((char*)"ERROR: %s\n", doc.ErrorDesc());
        exit(doc.ErrorId());
    }

    TiXmlNode* task = 0;
    TiXmlElement* el = 0;
    TiXmlNode* node0 = 0;
    TiXmlNode* node1 = 0;
    task = doc.FirstChild("task");

    int steadyVal = 1;
    node0 = task->FirstChild("control");
//    node0->FirstChild("STEADY")->ToElement()->Attribute("value", &steadyVal);
//    node0->FirstChild("SIGMA")->ToElement()->Attribute("value", &SIGMA);
    node0->FirstChild("TAU")->ToElement()->Attribute("value", &TAU);
    node0->FirstChild("TMAX")->ToElement()->Attribute("value", &TMAX);
    node0->FirstChild("STEP_MAX")->ToElement()->Attribute("value", &STEP_MAX);
    node0->FirstChild("FILE_OUTPUT_STEP")->ToElement()->Attribute("value", &FILE_SAVE_STEP);
    node0->FirstChild("LOG_OUTPUT_STEP")->ToElement()->Attribute("value", &PRINT_STEP);

    STEADY = false;
//    const char * limiterName = node0->FirstChild("LIMITER")->ToElement()->Attribute("value");


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
//        node1->FirstChild("K")->ToElement()->Attribute("value", &mat.K);
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

        reg.name = regNode->FirstChild("name")->ToElement()->GetText();

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

        int zhrv = 0;
    }

    Material::gR *= R_ * T_ / P_;	// Газовая постоянная
    for (int i = 0; i < matCount; i++)
    {
        Material & mat = materials[i];
        mat.Cp /= CV_;
        mat.K /= KP_;
        mat.ML /= MU_;
    }



//    // чтение параметров о ГРАНИЧНЫХ УСЛОВИЯХ
//    node0 = task->FirstChild("boundaries");
//    node0->ToElement()->Attribute("count", &bCount);
//    boundaries = new Boundary[bCount];
//    TiXmlNode* bNode = node0->FirstChild("boundCond");
//    for (int i = 0; i < bCount; i++)
//    {
//        Boundary & b = boundaries[i];
//        bNode->ToElement()->Attribute("edgeType", &b.edgeType);
//        const char * str = bNode->FirstChild("type")->ToElement()->GetText();
//        if (strcmp(str, "BOUND_WALL") == 0)
//        {
//            b.parCount = 0;
//            b.par = NULL;
//            b.type = Boundary::BOUND_WALL;
//        }
//        else
//        if (strcmp(str, "BOUND_OUTLET") == 0)
//        {
//            b.parCount = 0;
//            b.par = NULL;
//            b.type = Boundary::BOUND_OUTLET;
//        }
//        else
//        if (strcmp(str, "BOUND_INLET") == 0)
//        {
//            b.parCount = 4;
//            b.par = new double[4];
//            b.type = Boundary::BOUND_INLET;
//
//            node1 = bNode->FirstChild("parameters");
//            node1->FirstChild("T")->ToElement()->Attribute("value", &b.par[0]); b.par[0] /= T_;
//            node1->FirstChild("P")->ToElement()->Attribute("value", &b.par[1]); b.par[1] /= P_;
//            node1->FirstChild("Vx")->ToElement()->Attribute("value", &b.par[2]); b.par[2] /= U_;
//            node1->FirstChild("Vy")->ToElement()->Attribute("value", &b.par[3]); b.par[3] /= U_;
//        }
//        else {
//            log((char*)"ERROR: unsupported boundary condition type '%s'", str);
//            EXIT(1);
//        }
//
//        bNode = bNode->NextSibling("boundCond");
//    }

    /* Чтение параметров ГУ */
    node0 = task->FirstChild("boundaries");
    TiXmlNode* bNode = node0->FirstChild("boundCond");
    while (bNode != NULL)
    {
        int edgeType;
        bNode->ToElement()->Attribute("edgeType", &edgeType);

        CFDBoundary * b;

        try {
            b = HeatBoundary::create(bNode, &grid);
        }
        catch (Exception e) {
            log((char*)"ERROR: %s\n", e.getMessage());
            exit(e.getType());
        }

        boundaries.push_back(b);

        bNode = bNode->NextSibling("boundCond");
    }

    bCount = boundaries.size();

    /* Чтение данных сетки. */
    node0 = task->FirstChild("mesh");
    const char* fName = node0->FirstChild("name")->ToElement()->Attribute("value");
    const char* tName = node0->FirstChild("filesType")->ToElement()->Attribute("value");
    MeshReader* mr = MeshReader::create(MeshReader::getType((char*)tName), (char*)fName);
    mr->read(&grid);

    /* Определение ГУ для каждого ребра. */
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
                log((char*)"ERROR (boundary condition): unknown edge type of edge %d...\n", iEdge);
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
                log((char*)"ERROR (boundary condition): unknown edge type of edge %d...\n", iEdge);
                EXIT(1);
            }

            e.bnd = boundaries[iBound];
        }
    }


    memAlloc();


    // инициализация решателя
    node0 = task->FirstChild("solver");
    const char * solverName = node0->ToElement()->Attribute("type");
    //solverMtx = new SolverZeidel();
    solverMtx = MatrixSolver::create(solverName);
    solverMtx->init(grid.cCount, MATR_DIM);
    log((char*)"Solver type: %s.\n", solverMtx->getName());

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
        Cell & c = grid.cells[i];
        Region & reg = getRegion(c.typeName);
        convertParToCons(i, reg.par);
    }

    calcTimeStep();
    //log("TAU_MIN = %25.16e\n", TAU_MIN);

    save(0);
}

void HEAT_DG_IMPLICIT::calcMassMatr()
{
    for (int iCell = 0; iCell < grid.cCount; iCell++) {
        double **	A  = matrA[iCell];
        double **	Bx = matrBx[iCell];
        double **	By = matrBy[iCell];
        for (int i = 0; i < BASE_FUNC_COUNT; i++) {
            for (int j = 0; j < BASE_FUNC_COUNT; j++) {
                A[i][j]  = 0.0;
                Bx[i][j] = 0.0;
                By[i][j] = 0.0;
                for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
                    A[i][j] += cellGW[iCell][iGP] * getF(i, iCell, cellGP[iCell][iGP].x, cellGP[iCell][iGP].y)
                               * getF(j, iCell, cellGP[iCell][iGP].x, cellGP[iCell][iGP].y);
                    Bx[i][j] += cellGW[iCell][iGP] * getF(i, iCell, cellGP[iCell][iGP].x, cellGP[iCell][iGP].y)
                               * getDfDx(j, iCell, cellGP[iCell][iGP].x, cellGP[iCell][iGP].y);
                    By[i][j] += cellGW[iCell][iGP] * getF(i, iCell, cellGP[iCell][iGP].x, cellGP[iCell][iGP].y)
                               * getDfDy(j, iCell, cellGP[iCell][iGP].x, cellGP[iCell][iGP].y);
                }
                A[i][j] *= cellJ[iCell];
            }
        }
    }
}

void HEAT_DG_IMPLICIT::calcGaussPar()
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

void HEAT_DG_IMPLICIT::memAlloc()
{
    int n = grid.cCount;
    cTau = new double[n];

    u = new double*[n];
    qx = new double*[n];
    qy = new double*[n];

    cellGP = new Point*[n];
    cellGW = new double*[n];
    cellJ = new double[n];

    edgeGW = new double*[grid.eCount];
    edgeJ  = new double[grid.eCount];
    edgeGP = new Point*[grid.eCount];

    matrA		= new double**[n];
    matrBx		= new double**[n];
    matrBy		= new double**[n];
    matrInvA	= new double**[n];

    for (int i = 0; i < n; i++) {
        u[i] = new double[BASE_FUNC_COUNT];
        qx[i] = new double[BASE_FUNC_COUNT];
        qy[i] = new double[BASE_FUNC_COUNT];

        cellGP[i] = new Point[GP_CELL_COUNT];
        cellGW[i] = new double[GP_CELL_COUNT];

        matrA[i]  = allocMtx(BASE_FUNC_COUNT);
        matrBx[i] = allocMtx(BASE_FUNC_COUNT);
        matrBy[i] = allocMtx(BASE_FUNC_COUNT);
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

    fields = new double**[3];
    fields[FIELD_U]  = u;
    fields[FIELD_QX] = qx;
    fields[FIELD_QY] = qy;

    matrSmall = allocMtx(BASE_FUNC_COUNT);
    matrSmall1 = allocMtx(BASE_FUNC_COUNT);
    matrSmall2 = allocMtx(BASE_FUNC_COUNT);
    matrBig = allocMtx(MATR_DIM);
    matrBig2 = allocMtx(MATR_DIM);

}

void HEAT_DG_IMPLICIT::memFree()
{
    int n = grid.cCount;

    for (int i = 0; i < n; i++) {
        delete[] u[i];
        delete[] qx[i];
        delete[] qy[i];

        delete[] cellGP[i];
        delete[] cellGW[i];

        freeMtx(matrA[i], BASE_FUNC_COUNT);
        freeMtx(matrBx[i], BASE_FUNC_COUNT);
        freeMtx(matrBy[i], BASE_FUNC_COUNT);
        freeMtx(matrInvA[i], BASE_FUNC_COUNT);
    }
    delete[] cTau;

    delete[] u;
    delete[] qx;
    delete[] qy;

    delete[] cellGP;
    delete[] cellGW;
    delete[] cellJ;

    delete[] matrA;
    delete[] matrBx;
    delete[] matrBy;
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

Region & HEAT_DG_IMPLICIT::getRegionByCellType(int type)
{
    for (int i = 0; i < regCount; i++)
    {
        if (regions[i].cellType == type) return regions[i];
    }
    log((char*)"ERROR: unknown cell type %d...\n", type);
    EXIT(1);
}

Region & HEAT_DG_IMPLICIT::getRegionByName(char* name)
{
    for (int i = 0; i < regCount; i++)
    {
        if (strcmp(regions[i].name.c_str(), name) == 0) return regions[i];
    }
    log((char*)"ERROR: unknown cell name '%d'...\n", name);
    EXIT(1);
}

Region   &	HEAT_DG_IMPLICIT::getRegion(int iCell)
{
    return getRegionByName(grid.cells[iCell].typeName);
    //return getRegionByCellType(grid.cells[iCell].type);
}

Region & HEAT_DG_IMPLICIT::getRegion(char * name)
{
    return getRegionByName(name);
}

Material &	HEAT_DG_IMPLICIT::getMaterial(int iCell)
{
    Region & reg = getRegion(iCell);
    return materials[reg.matId];
}

void HEAT_DG_IMPLICIT::convertParToCons(int iCell, Param & par)
{
    memset(u[iCell], 0, sizeof(double)*BASE_FUNC_COUNT);
    memset(qx[iCell], 0, sizeof(double)*BASE_FUNC_COUNT);
    memset(qy[iCell], 0, sizeof(double)*BASE_FUNC_COUNT);
    u[iCell][0] = par.T;
}

void HEAT_DG_IMPLICIT::convertConsToPar(int iCell, Param & par) {
    double fT = getField(FIELD_U, iCell, grid.cells[iCell].c);

    par.T = fT;
    par.u = 0;
    par.v = 0;
    par.E = 0;
    par.e = 0;
    par.r = 0;
    par.p = 0;
    par.cz = 0;
    par.ML = 0;
}



double HEAT_DG_IMPLICIT::getField(int fldId, int iCell, double x, double y)
{
    double *fld = fields[fldId][iCell];
    double result = 0.0;
    for (int i = 0; i < BASE_FUNC_COUNT; i++) {
        result += fld[i] * getF(i, iCell, x, y);
    }
    return result;
}

double HEAT_DG_IMPLICIT::getField(int fldId, int iCell, Point p)
{
    return getField(fldId, iCell, p.x, p.y);
}

void HEAT_DG_IMPLICIT::getFields(double &fRO, double &fRU, double &fRV, double &fRE, int iCell, double x, double y)
{
    fRO = getField(FIELD_U, iCell, x, y);
    fRU = getField(FIELD_QX, iCell, x, y);
    fRV = getField(FIELD_QY, iCell, x, y);
}

void HEAT_DG_IMPLICIT::getFields(double &fRO, double &fRU, double &fRV, double &fRE, int iCell, Point p)
{
    getFields(fRO, fRU, fRV, fRE, iCell, p.x, p.y);
}

double HEAT_DG_IMPLICIT::getF(int id, int iCell, Point p)
{
    return getF(id, iCell, p.x, p.y);
}

double HEAT_DG_IMPLICIT::getF(int id, int iCell, double x, double y)
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

double HEAT_DG_IMPLICIT::getDfDx(int id, int iCell, Point p)
{
    return getDfDx(id, iCell, p.x, p.y);
}

double HEAT_DG_IMPLICIT::getDfDx(int id, int iCell, double x, double y)
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

double HEAT_DG_IMPLICIT::getDfDy(int id, int iCell, Point p)
{
    return getDfDy(id, iCell, p.x, p.y);
}

double HEAT_DG_IMPLICIT::getDfDy(int id, int iCell, double x, double y)
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

void HEAT_DG_IMPLICIT::calcTimeStep()
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

void HEAT_DG_IMPLICIT::save(int step)
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

    fprintf(fp, "CELL_DATA %d\nSCALARS Temperature float 1\nLOOKUP_TABLE default\n", grid.cCount);
    for (int i = 0; i < grid.cCount; i++)
    {
        Param p;
        convertConsToPar(i, p);
        fprintf(fp, "%25.16f ", p.T*T_);
        if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
    }


    fprintf(fp, "SCALARS TAU float 1\nLOOKUP_TABLE default\n");
    for (int i = 0; i < grid.cCount; i++)
    {
        fprintf(fp, "%25.16f ", cTau[i]*TIME_);
        if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
    }

    fclose(fp);
    printf("File '%s' saved...\n", fName);

}

//int HEAT_DG_IMPLICIT::getLimitedCellsCount() {
//    int n = 0;
//    for (int iCell = 0; iCell < grid.cCount; iCell++) {
//        if ((grid.cells[iCell].flag & CELL_FLAG_LIM) > 0) n++;
//    }
//    return n;
//}
//
//void HEAT_DG_IMPLICIT::remediateLimCells()
//{
////    for (int iCell = 0; iCell < grid.cCount; iCell++)
////    {
////        if (cellIsLim(iCell))
////        {
////            // пересчитываем по соседям
////            double sRO = 0.0;
////            double sRU = 0.0;
////            double sRV = 0.0;
////            double sRE = 0.0;
////            double S = 0.0;
////            for (int i = 0; i < grid.cells[iCell].eCount; i++)
////            {
////                int		iEdge = grid.cells[iCell].edgesInd[i];
////                int		j = grid.edges[iEdge].c2;
////                if (j == iCell)	{
////                    //std::swap(j, grid.edges[iEdge].c1); // так нужно еще нормаль поворчивать тогда
////                    j = grid.edges[iEdge].c1;
////                }
////                if (j >= 0) {
////                    double  s = grid.cells[j].S;
////                    S += s;
////                    sRO += getField(FIELD_RO, j, grid.cells[j].c) * s;
////                    sRU += getField(FIELD_RO, j, grid.cells[j].c) * s;
////                    sRV += getField(FIELD_RO, j, grid.cells[j].c) * s;
////                    sRE += getField(FIELD_RO, j, grid.cells[j].c) * s;
////                }
////            }
////            if (S >= TAU*TAU) {
////                memset(ro[iCell], 0, sizeof(double)*BASE_FUNC_COUNT);
////                memset(ru[iCell], 0, sizeof(double)*BASE_FUNC_COUNT);
////                memset(rv[iCell], 0, sizeof(double)*BASE_FUNC_COUNT);
////                memset(re[iCell], 0, sizeof(double)*BASE_FUNC_COUNT);
////
////                ro[iCell][0] = sRO / S;
////                ru[iCell][0] = sRU / S;
////                rv[iCell][0] = sRV / S;
////                re[iCell][0] = sRE / S;
////            }
////            // после 0x20 итераций пробуем вернуть ячейку в счет
////            grid.cells[iCell].flag += 0x010000;
////            if (grid.cells[iCell].flag & 0x200000) grid.cells[iCell].flag &= 0x001110;
////        }
////    }
//}


void HEAT_DG_IMPLICIT::decCFL()
{
    if (CFL > 0.01) {
        CFL *= 0.75;
        if (CFL < 0.01) CFL = 0.01;
        log((char*)" CFL Number has been decreased : %25.25e \n", CFL);
    }
    //log(" <<< CFL Number has been decreased : %25.25e \n", CFL);
}

void HEAT_DG_IMPLICIT::incCFL()
{
    if (CFL < maxCFL) {
        CFL *= scaleCFL;
        if (CFL > maxCFL) CFL = maxCFL;
        log((char*)" CFL Number has been increased : %25.25e \n", CFL);
    }
}

//double **HEAT_DG_IMPLICIT::allocMtx4()
//{
//    double		**tempMtx4 = new double*[4];
//    for (int i = 0; i < 4; ++i) tempMtx4[i] = new double[4];
//    return tempMtx4;
//}
//void   HEAT_DG_IMPLICIT::freeMtx4(double **mtx4)
//{
//    for (int i = 0; i < 4; ++i)
//        delete[] mtx4[i];
//    delete[] mtx4;
//}
//void HEAT_DG_IMPLICIT::multMtx4(double **dst4, double **srcA4, double **srcB4)
//{
//    double sum;
//    for (int i = 0; i < 4; ++i)
//    {
//        for (int j = 0; j < 4; ++j)
//        {
//            sum = 0;
//            for (int k = 0; k < 4; ++k)
//                sum += srcA4[i][k] * srcB4[k][j];
//            dst4[i][j] = sum;
//        }
//    }
//}
//void HEAT_DG_IMPLICIT::clearMtx4(double **mtx4)
//{
//    for (int i = 0; i < 4; ++i)
//        for (int j = 0; j < 4; ++j)
//            mtx4[i][j] = 0;
//}

void HEAT_DG_IMPLICIT::multMtxToVal(double **dst, double x, int N)
{
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            dst[i][j] *= x;
        }
    }
}

void HEAT_DG_IMPLICIT::multMtxToVec(double *dst, double **mtx, double *vec, int N)
{
    for (int i = 0; i < N; ++i)
    {
        dst[i] = 0.0;
        for (int j = 0; j < N; ++j)
        {
            dst[i] += mtx[i][j]*vec[j];
        }
    }
}

void HEAT_DG_IMPLICIT::fillMtx(double** dst, double x, int N)
{
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            dst[i][j] = x;
        }
    }
}


void HEAT_DG_IMPLICIT::calcMatrWithTau()
{
//    for (int iCell = 0; iCell < grid.cCount; iCell++) {
//
//        fillMtx(matrBig, 0.0, MATR_DIM);
//
//
//        //for (int ii = 0; ii < FIELD_COUNT; ii++) {
//        //	addSmallMatrToBigMatr(matrBig, matrA[iCell], ii, ii);
//        //}
//
//        solverMtx->addMatrElement(iCell, iCell, matrBig);
//
//    }
}

void HEAT_DG_IMPLICIT::calcIntegral()
{
    for (int iCell = 0; iCell < grid.cCount; iCell++) {

        fillMtx(matrBig, 0.0, MATR_DIM);

        for (int i = 0; i < BASE_FUNC_COUNT; i++) {
            for (int j = 0; j < BASE_FUNC_COUNT; j++) {
                matrSmall[i][j]  = matrA[iCell][i][j] / cTau[iCell];
                matrSmall1[i][j] = -matrBx[iCell][i][j];
                matrSmall2[i][j] = -matrBy[iCell][i][j];
            }
        }

        addSmallMatrToBigMatr(matrBig, matrSmall, 0, 0);
        for (int ii = 1; ii < FIELD_COUNT; ii++) {
            addSmallMatrToBigMatr(matrBig, matrA[iCell], ii, ii);
        }

        addSmallMatrToBigMatr(matrBig, matrSmall1, 0, 1);
        addSmallMatrToBigMatr(matrBig, matrSmall2, 0, 2);

        addSmallMatrToBigMatr(matrBig, matrSmall1, 1, 0);
        addSmallMatrToBigMatr(matrBig, matrSmall2, 2, 0);

        solverMtx->addMatrElement(iCell, iCell, matrBig);
    }
}

void HEAT_DG_IMPLICIT::calcMatrFlux()
{
    double** mtxII = matrSmall1;
    double** mtxKI = matrSmall2;
    double** mtxTmp = matrSmall;

    for (int iEdge = 0; iEdge < grid.eCount; iEdge++) {
        Edge& edge = grid.edges[iEdge];
        Vector&	n = grid.edges[iEdge].n;
        int c1 = edge.c1;
        int c2 = edge.c2;

        if (c2 >= 0) {

            for (int i = 0; i < BASE_FUNC_COUNT; i++) {
                for (int j = 0; j < BASE_FUNC_COUNT; j++) {
                    mtxII[i][j] = 0.0;
                    mtxKI[i][j] = 0.0;
                    for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
                        mtxII[i][j] += edgeGW[iEdge][iGP] * getF(i, c1, edgeGP[iEdge][iGP].x, edgeGP[iEdge][iGP].y)
                                                          * getF(j, c1, edgeGP[iEdge][iGP].x, edgeGP[iEdge][iGP].y);
                        mtxKI[i][j] += edgeGW[iEdge][iGP] * getF(i, c2, edgeGP[iEdge][iGP].x, edgeGP[iEdge][iGP].y)
                                                          * getF(j, c1, edgeGP[iEdge][iGP].x, edgeGP[iEdge][iGP].y);

                    }
                    mtxII[i][j] *= edgeJ[iEdge];
                    mtxKI[i][j] *= edgeJ[iEdge];
                }
            }

            fillMtx(matrBig, 0.0, MATR_DIM);
            fillMtx(matrBig2, 0.0, MATR_DIM);

            copyMtx(mtxTmp, mtxII, BASE_FUNC_COUNT);
            multMtxToVal(mtxTmp, n.x*0.5, BASE_FUNC_COUNT);
            addSmallMatrToBigMatr(matrBig, mtxTmp, 0, 1);
            addSmallMatrToBigMatr(matrBig, mtxTmp, 1, 0);
            copyMtx(mtxTmp, mtxII, BASE_FUNC_COUNT);
            multMtxToVal(mtxTmp, n.y*0.5, BASE_FUNC_COUNT);
            addSmallMatrToBigMatr(matrBig, mtxTmp, 0, 2);
            addSmallMatrToBigMatr(matrBig, mtxTmp, 2, 0);

            copyMtx(mtxTmp, mtxKI, BASE_FUNC_COUNT);
            multMtxToVal(mtxTmp, n.x*0.5, BASE_FUNC_COUNT);
            addSmallMatrToBigMatr(matrBig2, mtxTmp, 0, 1);
            addSmallMatrToBigMatr(matrBig2, mtxTmp, 1, 0);
            copyMtx(mtxTmp, mtxKI, BASE_FUNC_COUNT);
            multMtxToVal(mtxTmp, n.y*0.5, BASE_FUNC_COUNT);
            addSmallMatrToBigMatr(matrBig2, mtxTmp, 0, 2);
            addSmallMatrToBigMatr(matrBig2, mtxTmp, 2, 0);

            solverMtx->addMatrElement(c1, c1, matrBig);
            solverMtx->addMatrElement(c1, c2, matrBig2);


            fillMtx(matrBig, 0.0, MATR_DIM);
            fillMtx(matrBig2, 0.0, MATR_DIM);

            copyMtx(mtxTmp, mtxII, BASE_FUNC_COUNT);
            multMtxToVal(mtxTmp, -n.x*0.5, BASE_FUNC_COUNT);
            addSmallMatrToBigMatr(matrBig, mtxTmp, 0, 1);
            addSmallMatrToBigMatr(matrBig, mtxTmp, 1, 0);
            copyMtx(mtxTmp, mtxII, BASE_FUNC_COUNT);
            multMtxToVal(mtxTmp, -n.y*0.5, BASE_FUNC_COUNT);
            addSmallMatrToBigMatr(matrBig, mtxTmp, 0, 2);
            addSmallMatrToBigMatr(matrBig, mtxTmp, 2, 0);

            copyMtx(mtxTmp, mtxKI, BASE_FUNC_COUNT);
            multMtxToVal(mtxTmp, -n.x*0.5, BASE_FUNC_COUNT);
            addSmallMatrToBigMatr(matrBig2, mtxTmp, 0, 1);
            addSmallMatrToBigMatr(matrBig2, mtxTmp, 1, 0);
            copyMtx(mtxTmp, mtxKI, BASE_FUNC_COUNT);
            multMtxToVal(mtxTmp, -n.y*0.5, BASE_FUNC_COUNT);
            addSmallMatrToBigMatr(matrBig2, mtxTmp, 0, 2);
            addSmallMatrToBigMatr(matrBig2, mtxTmp, 2, 0);

            solverMtx->addMatrElement(c2, c2, matrBig);
            solverMtx->addMatrElement(c2, c1, matrBig2);
        }
    }

}

void HEAT_DG_IMPLICIT::calcRHS()
{

    for (int iCell = 0; iCell < grid.cCount; iCell++) {
        memset(tmpArr, 0, sizeof(double)*MATR_DIM);
        multMtxToVec(tmpArr, matrA[iCell], u[iCell], BASE_FUNC_COUNT);
        for (int i = 0; i < BASE_FUNC_COUNT; i++) {
            tmpArr[i] /= cTau[iCell];
        }
        solverMtx->addRightElement(iCell, tmpArr);
    }

    for (int iEdge = 0; iEdge < grid.eCount; iEdge++) {
        Edge &edge = grid.edges[iEdge];
        Vector &n = grid.edges[iEdge].n;
        int c1 = edge.c1;
        int c2 = edge.c2;

        if (c2 < 0) {
            memset(tmpArr, 0, sizeof(double)*MATR_DIM);

            for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
                Param p1, p2;
                double fU  = getField(FIELD_U,  c1, edgeGP[iEdge][iGP]);
                double fQx = getField(FIELD_QX, c1, edgeGP[iEdge][iGP]);
                double fQy = getField(FIELD_QY, c1, edgeGP[iEdge][iGP]);

                p1.T = fU;
                p1.Qt[0] = fQx;
                p1.Qt[1] = fQy;


                edge.bnd->run(iEdge, p1, p2);

                fU  = 0.5*(p1.T+p2.T);
                fQx = 0.5*(p1.Qt[0]+p2.Qt[0]);
                fQy = 0.5*(p1.Qt[1]+p2.Qt[1]);

                for (int m = 0; m < BASE_FUNC_COUNT; m++) {
                    tmpArr[ m ]                    += (fQx*n.x+fQy*n.y)*getF(m, c1, edgeGP[iEdge][iGP]);
                    tmpArr[ m + BASE_FUNC_COUNT]   += fU*n.x*getF(m, c1, edgeGP[iEdge][iGP]);
                    tmpArr[ m + BASE_FUNC_COUNT*2] += fU*n.y*getF(m, c1, edgeGP[iEdge][iGP]);
                }
            }

            solverMtx->addRightElement(c1, tmpArr);
        }
    }

}


void HEAT_DG_IMPLICIT::run()
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

        /* Заполняем элементы матрицы */
//        calcMatrWithTau();		// вычисляем матрицы перед производной по времени
        calcIntegral();			// вычисляем интеграл от(dF / dU)*deltaU*dFi / dx
        calcMatrFlux();			// Вычисляем потоковые величины

        /* Задаем начальное приближение */
        for (int cellIndex = 0, ind = 0; cellIndex < grid.cCount; cellIndex++)
        {
            for (int iFld = 0; iFld < FIELD_COUNT; iFld++) {
                for (int iF = 0; iF < BASE_FUNC_COUNT; iF++) {
                    solverMtx->x[ind++] = fields[iFld][cellIndex][iF];
                }
            }
        }

        /* Решаем СЛАУ */
        int maxIter = SOLVER_ITER;
        double eps = SOLVER_EPS;

        solverErr = solverMtx->solve(eps, maxIter);

        solverMtx->printToFile("matr.txt");
        if (solverErr == MatrixSolver::RESULT_OK) {

            //if (SMOOTHING) smoothingDelta(solverMtx->x);

            for (int cellIndex = 0, ind = 0; cellIndex < grid.cCount; cellIndex++)
            {
                Cell &cell = grid.cells[cellIndex];

                //if (cellIsLim(cellIndex))	continue;
                for (int iFld = 0; iFld < FIELD_COUNT; iFld++) {
                    for (int iF = 0; iF < BASE_FUNC_COUNT; iF++) {
                        fields[iFld][cellIndex][iF] = solverMtx->x[ind++];
                    }
                }
            }


//            if (limiter != NULL) {
//                limiter->run();
//            }

//            for (int cellIndex = 0, ind = 0; cellIndex < grid.cCount; cellIndex++)
//            {
//                Cell &cell = grid.cells[cellIndex];
//
//                Param par;
//                convertConsToPar(cellIndex, par);
//                if (par.r > limitRmax)			{ par.r = limitRmax;	setCellFlagLim(cellIndex); }
//                if (par.r < limitRmin)			{ par.r = limitRmin;	setCellFlagLim(cellIndex); }
//                if (fabs(par.u) > limitUmax)		{ par.u = limitUmax*par.u / fabs(par.u);	setCellFlagLim(cellIndex); }
//                if (fabs(par.v) > limitUmax)		{ par.v = limitUmax*par.v / fabs(par.v);	setCellFlagLim(cellIndex); }
//                if (par.p > limitPmax)			{ par.p = limitPmax;	setCellFlagLim(cellIndex); }
//                if (par.p < limitPmin)			{ par.p = limitPmin;	setCellFlagLim(cellIndex); }
//                if (cellIsLim(cellIndex)) 		{ par.e = par.p / ((getGAM(cellIndex) - 1)*par.r); convertParToCons(cellIndex, par); }
//
//                double fRO, fRU, fRV, fRE;
//                for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
//
//                    fRO = getField(FIELD_RO, cellIndex, cellGP[cellIndex][iGP]);
//                    fRU = getField(FIELD_RU, cellIndex, cellGP[cellIndex][iGP]);
//                    fRV = getField(FIELD_RV, cellIndex, cellGP[cellIndex][iGP]);
//                    fRE = getField(FIELD_RE, cellIndex, cellGP[cellIndex][iGP]);
//                    consToPar(fRO, fRU, fRV, fRE, par);
//                    if (par.r > limitRmax)			{ par.r = limitRmax;	setCellFlagLim(cellIndex); }
//                    if (par.r < limitRmin)			{ par.r = limitRmin;	setCellFlagLim(cellIndex); }
//                    if (fabs(par.u) > limitUmax)		{ par.u = limitUmax*par.u / fabs(par.u);	setCellFlagLim(cellIndex); }
//                    if (fabs(par.v) > limitUmax)		{ par.v = limitUmax*par.v / fabs(par.v);	setCellFlagLim(cellIndex); }
//                    if (par.p > limitPmax)			{ par.p = limitPmax;	setCellFlagLim(cellIndex); }
//                    if (par.p < limitPmin)			{ par.p = limitPmin;	setCellFlagLim(cellIndex); }
//                    if (cellIsLim(cellIndex)) 		{ par.e = par.p / ((getGAM(cellIndex) - 1)*par.r); convertParToCons(cellIndex, par); }
//                }
//
//                for (int iEdge = 0; iEdge < cell.eCount; iEdge++) {
//                    int edgInd = cell.edgesInd[iEdge];
//                    for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
//                        fRO = getField(FIELD_RO, cellIndex, edgeGP[edgInd][iGP]);
//                        fRU = getField(FIELD_RU, cellIndex, edgeGP[edgInd][iGP]);
//                        fRV = getField(FIELD_RV, cellIndex, edgeGP[edgInd][iGP]);
//                        fRE = getField(FIELD_RE, cellIndex, edgeGP[edgInd][iGP]);
//                        consToPar(fRO, fRU, fRV, fRE, par);
//                        if (par.r > limitRmax)			{ par.r = limitRmax;	setCellFlagLim(cellIndex); }
//                        if (par.r < limitRmin)			{ par.r = limitRmin;	setCellFlagLim(cellIndex); }
//                        if (fabs(par.u) > limitUmax)		{ par.u = limitUmax*par.u / fabs(par.u);	setCellFlagLim(cellIndex); }
//                        if (fabs(par.v) > limitUmax)		{ par.v = limitUmax*par.v / fabs(par.v);	setCellFlagLim(cellIndex); }
//                        if (par.p > limitPmax)			{ par.p = limitPmax;	setCellFlagLim(cellIndex); }
//                        if (par.p < limitPmin)			{ par.p = limitPmin;	setCellFlagLim(cellIndex); }
//                        if (cellIsLim(cellIndex)) 		{ par.e = par.p / ((getGAM(cellIndex) - 1)*par.r); convertParToCons(cellIndex, par); }
//                    }
//                }
//
//
//            }
//            remediateLimCells();
//
//            int limCells = getLimitedCellsCount();
//            if (STEADY && (limCells >= maxLimCells)) decCFL();

            timeEnd = clock();
            totalCalcTime += (timeEnd - timeStart);
            if (step % FILE_SAVE_STEP == 0)
            {
                save(step);
            }
            if (step % PRINT_STEP == 0)
            {
                if (!STEADY) {

                    log((char*)"step: %6d  time step: %.16f\tmax iter: %5d\tlim: %4d\tlift force (Fx, Fy) = (%.16f, %.16f)\ttime: %6d ms\ttotal calc time: %ld\n", step, t, maxIter, limCells, Fx, Fy, timeEnd - timeStart, totalCalcTime);
                }
                else {
                    log((char*)"step: %6d  max iter: %5d\tlim: %4d\tlift force (Fx, Fy) = (%.16f, %.16f)\ttime: %6d ms\ttotal calc time: %ld\n", step, maxIter, limCells, Fx, Fy, timeEnd - timeStart, totalCalcTime);
                }
            }

            if (STEADY && (step % stepCFL == 0)) incCFL();
        }
        else {
            if (solverErr & MatrixSolver::RESULT_ERR_CONVERG) {
                log((char*)"Solver error: residual condition not considered.\n");
            }
            if (solverErr & MatrixSolver::RESULT_ERR_MAX_ITER) {
                log((char*)"Solver error: max iterations done.\n");
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

void HEAT_DG_IMPLICIT::addSmallMatrToBigMatr(double **mB, double **mS, int i, int j)
{

    int ii = i * BASE_FUNC_COUNT;
    int jj = j * BASE_FUNC_COUNT;
    for (int i1 = 0; i1 < BASE_FUNC_COUNT; i1++) {
        for (int j1 = 0; j1 < BASE_FUNC_COUNT; j1++) {
            mB[ii + i1][jj + j1] += mS[i1][j1];
        }
    }
}

void HEAT_DG_IMPLICIT::done()
{
}

