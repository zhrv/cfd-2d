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

#define _PI_ 3.14159265358979323846

//double phi_init(double x, double y)
//{
//    return sin(_PI_*x)*sin(_PI_*y);
//}
//
//double phi_init_dx(double x, double y)
//{
//    return _PI_*cos(_PI_*x)*sin(_PI_*y);
//}
//
//double phi_init_dy(double x, double y)
//{
//    return _PI_*sin(_PI_*x)*cos(_PI_*y);
//}
//
//double phi_exac(double t, double x, double y)
//{
//    return exp(-_PI_*_PI_*t*2.0)*sin(_PI_*x)*sin(_PI_*y);
//}


double phi_init(double x, double y)
{
    return sin(2.0*_PI_*(x+y-1));
}

double phi_init_dx(double x, double y)
{
    return 2.0*_PI_*cos(2.0*_PI_*(x+y-1));
}

double phi_init_dy(double x, double y)
{
    return 2.0*_PI_*cos(2.0*_PI_*(x+y-1));
}

double phi_exac(double t, double x, double y)
{
    return exp(-t*2.0)*sin(2.0*_PI_*(x+y-1)-2.0*t);
}

double phi_exac_dx(double t, double x, double y)
{
    return 2.0*_PI_*exp(-t*2.0)*cos(2.0*_PI_*(x+y-1)-2.0*t);
}

double phi_exac_dy(double t, double x, double y)
{
    return 2.0*_PI_*exp(-t*2.0)*cos(2.0*_PI_*(x+y-1)-2.0*t);
}


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
        node1->FirstChild("K")->ToElement()->Attribute("value", &mat.K);
        node1->FirstChild("ML")->ToElement()->Attribute("value", &mat.ML);
        node1->FirstChild("rho")->ToElement()->Attribute("value", &mat.rho);
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
        reg.par.r = mat.rho;
        //mat.URS(reg.par, 3);
        reg.par.Qt[0] = reg.par.Qt[1] = 0.0;

        regNode = regNode->NextSibling("region");

        if (reg.par.p > maxP) maxP = reg.par.p;
        if (reg.par.r > maxR) maxR = reg.par.r;
        if (reg.par.T > maxT) maxT = reg.par.T;
        if (reg.par.U2() > maxU) maxU = reg.par.U2();
    }

    maxU = sqrt(maxU);

    //TODO
//    maxR = 1.0;
//    maxP = 1.0;
//    maxT = 1.0;

    // параметры обезразмеривания
    L_ = 150.0;
    R_ = maxR;					// характерная плотность = начальная плотность
    P_ = maxP;					// характерное давление = начальное давление 
    T_ = maxT;					// характерная температура = начальная температура
    U_ = maxU;
    //P_ = U_*U_*R_;
    U_ = sqrt(P_ / R_);		// характерная скорость = sqrt( P_ / R_ )
    E_ = POW_2(U_);			// характерная энергия  = U_**2
    CV_ = POW_2(U_) / T_;	// характерная теплоёмкость  = U_**2 / T_
    TIME_ = L_ / U_;			// характерное время
    MU_ = R_ * U_ * L_;		// характерная вязкость = R_ * U_ * L_
    //KP_ = R_ * POW_2(U_) * U_ * L_ / T_;	// коэффициент теплопроводности = R_ * U_**3 * L_ / T_
    KP_ = U_ * L_;	// коэффициент теплопроводности = R_ * U_**3 * L_ / T_
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

        par.r  /= R_;
        par.e  /= E_;
        par.cz /= U_;
        par.ML /= MU_;
        par.E   = par.e + par.U2()*0.5;

        int zhrv = 0;
    }

    Material::gR *= R_ * T_ / P_;	// Газовая постоянная
    for (int i = 0; i < matCount; i++)
    {
        Material & mat = materials[i];
        mat.Cp /= CV_;
        mat.K  /= KP_;
        mat.ML /= MU_;
    }



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
        Param &par = reg.par;
        convertParToCons(i, par);
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
                    Bx[i][j] += cellGW[iCell][iGP] * getF(j, iCell, cellGP[iCell][iGP].x, cellGP[iCell][iGP].y)
                               * getDfDx(i, iCell, cellGP[iCell][iGP].x, cellGP[iCell][iGP].y);
                    By[i][j] += cellGW[iCell][iGP] * getF(j, iCell, cellGP[iCell][iGP].x, cellGP[iCell][iGP].y)
                               * getDfDy(i, iCell, cellGP[iCell][iGP].x, cellGP[iCell][iGP].y);
                }
                A[i][j]  *= cellJ[iCell];
                Bx[i][j] *= cellJ[iCell];
                By[i][j] *= cellJ[iCell];
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
    qx[iCell][0] = par.Qt[0];
    qy[iCell][1] = par.Qt[1];
}

void HEAT_DG_IMPLICIT::convertConsToPar(int iCell, Param & par) {
    double fT = getField(FIELD_U, iCell, grid.cells[iCell].c);
    double fQx = getField(FIELD_QX, iCell, grid.cells[iCell].c);
    double fQy = getField(FIELD_QY, iCell, grid.cells[iCell].c);

    par.T = fT;
    par.Qt[0] = fQx;
    par.Qt[1] = fQy;
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

void HEAT_DG_IMPLICIT::getFields(double &fT, double &fQX, double &fQY, int iCell, double x, double y)
{
    fT  = getField(FIELD_U, iCell, x, y);
    fQX = getField(FIELD_QX, iCell, x, y);
    fQY = getField(FIELD_QY, iCell, x, y);
}

void HEAT_DG_IMPLICIT::getFields(double &fT, double &fQX, double &fQY, int iCell, Point p)
{
    getFields(fT, fQX, fQY, iCell, p.x, p.y);
}

double HEAT_DG_IMPLICIT::getF(int id, int iCell, Point p)
{
    return getF(id, iCell, p.x, p.y);
}

double HEAT_DG_IMPLICIT::getF(int id, int iCell, double x, double y)
{
    Point  &c  = grid.cells[iCell].c;
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
    for (int iCell = 0; iCell < grid.cCount; iCell++)
    {
        cTau[iCell] = TAU;
    }
}

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


void HEAT_DG_IMPLICIT::calcIntegral()
{
    for (int iCell = 0; iCell < grid.cCount; iCell++) {

        double k = getMaterial(iCell).K;


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

        for (int i = 0; i < BASE_FUNC_COUNT; i++) {
            for (int j = 0; j < BASE_FUNC_COUNT; j++) {
                matrSmall1[i][j] = -k*matrBx[iCell][i][j];
                matrSmall2[i][j] = -k*matrBy[iCell][i][j];
            }
        }
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

            Material &mat1 = getMaterial(c1);
            Material &mat2 = getMaterial(c2);

            double k1 = mat1.K;
            double k2 = mat2.K;

            /* C1 */

            for (int i = 0; i < BASE_FUNC_COUNT; i++) {
                for (int j = 0; j < BASE_FUNC_COUNT; j++) {
                    mtxII[i][j] = 0.0;
                    mtxKI[i][j] = 0.0;
                    for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
                        mtxII[i][j] += edgeGW[iEdge][iGP] * getF(i, c1, edgeGP[iEdge][iGP].x, edgeGP[iEdge][iGP].y)
                                       * getF(j, c1, edgeGP[iEdge][iGP].x, edgeGP[iEdge][iGP].y);
                        mtxKI[i][j] += edgeGW[iEdge][iGP] * getF(i, c1, edgeGP[iEdge][iGP].x, edgeGP[iEdge][iGP].y)
                                       * getF(j, c2, edgeGP[iEdge][iGP].x, edgeGP[iEdge][iGP].y);

                    }
                    mtxII[i][j] *= edgeJ[iEdge];
                    mtxKI[i][j] *= edgeJ[iEdge];
                }
            }

            fillMtx(matrBig, 0.0, MATR_DIM);
            fillMtx(matrBig2, 0.0, MATR_DIM);

//            if (n.x+n.y>=0) {
                copyMtx(mtxTmp, mtxII, BASE_FUNC_COUNT);
                multMtxToVal(mtxTmp, n.x, BASE_FUNC_COUNT);
                multMtxToVal(mtxTmp, k1, BASE_FUNC_COUNT);
                addSmallMatrToBigMatr(matrBig, mtxTmp, 1, 0);

                copyMtx(mtxTmp, mtxII, BASE_FUNC_COUNT);
                multMtxToVal(mtxTmp, n.y, BASE_FUNC_COUNT);
                multMtxToVal(mtxTmp, k1, BASE_FUNC_COUNT);
                addSmallMatrToBigMatr(matrBig, mtxTmp, 2, 0);


                copyMtx(mtxTmp, mtxKI, BASE_FUNC_COUNT);
                multMtxToVal(mtxTmp, n.x, BASE_FUNC_COUNT);
                addSmallMatrToBigMatr(matrBig2, mtxTmp, 0, 1);

                copyMtx(mtxTmp, mtxKI, BASE_FUNC_COUNT);
                multMtxToVal(mtxTmp, n.y, BASE_FUNC_COUNT);
                addSmallMatrToBigMatr(matrBig2, mtxTmp, 0, 2);
//            }
//            else {
//                copyMtx(mtxTmp, mtxII, BASE_FUNC_COUNT);
//                multMtxToVal(mtxTmp, n.x, BASE_FUNC_COUNT);
//                addSmallMatrToBigMatr(matrBig, mtxTmp, 0, 1);
//
//                copyMtx(mtxTmp, mtxII, BASE_FUNC_COUNT);
//                multMtxToVal(mtxTmp, n.y, BASE_FUNC_COUNT);
//                addSmallMatrToBigMatr(matrBig, mtxTmp, 0, 2);
//
//
//                copyMtx(mtxTmp, mtxKI, BASE_FUNC_COUNT);
//                multMtxToVal(mtxTmp, n.x, BASE_FUNC_COUNT);
//                multMtxToVal(mtxTmp, k2, BASE_FUNC_COUNT);
//                addSmallMatrToBigMatr(matrBig2, mtxTmp, 1, 0);
//
//                copyMtx(mtxTmp, mtxKI, BASE_FUNC_COUNT);
//                multMtxToVal(mtxTmp, n.y, BASE_FUNC_COUNT);
//                multMtxToVal(mtxTmp, k2, BASE_FUNC_COUNT);
//                addSmallMatrToBigMatr(matrBig2, mtxTmp, 2, 0);
//            }

            solverMtx->addMatrElement(c1, c1, matrBig);
            solverMtx->addMatrElement(c1, c2, matrBig2);


            /* C2 */

            for (int i = 0; i < BASE_FUNC_COUNT; i++) {
                for (int j = 0; j < BASE_FUNC_COUNT; j++) {
                    mtxII[i][j] = 0.0;
                    mtxKI[i][j] = 0.0;
                    for (int iGP = 0; iGP <  GP_EDGE_COUNT; iGP++) {
                        mtxII[i][j] += edgeGW[iEdge][iGP] * getF(i, c2, edgeGP[iEdge][iGP].x, edgeGP[iEdge][iGP].y)
                                       * getF(j, c2, edgeGP[iEdge][iGP].x, edgeGP[iEdge][iGP].y);
                        mtxKI[i][j] += edgeGW[iEdge][iGP] * getF(i, c2, edgeGP[iEdge][iGP].x, edgeGP[iEdge][iGP].y)
                                       * getF(j, c1, edgeGP[iEdge][iGP].x, edgeGP[iEdge][iGP].y);

                    }
                    mtxII[i][j] *= edgeJ[iEdge];
                    mtxKI[i][j] *= edgeJ[iEdge];
                }
            }
            fillMtx(matrBig, 0.0, MATR_DIM);
            fillMtx(matrBig2, 0.0, MATR_DIM);

//            if (-n.x-n.y>=0) {
//                copyMtx(mtxTmp, mtxII, BASE_FUNC_COUNT);
//                multMtxToVal(mtxTmp, -n.x, BASE_FUNC_COUNT);
//                multMtxToVal(mtxTmp, k2, BASE_FUNC_COUNT);
//                addSmallMatrToBigMatr(matrBig, mtxTmp, 1, 0);
//
//                copyMtx(mtxTmp, mtxII, BASE_FUNC_COUNT);
//                multMtxToVal(mtxTmp, -n.y, BASE_FUNC_COUNT);
//                multMtxToVal(mtxTmp, k2, BASE_FUNC_COUNT);
//                addSmallMatrToBigMatr(matrBig, mtxTmp, 2, 0);
//
//
//                copyMtx(mtxTmp, mtxKI, BASE_FUNC_COUNT);
//                multMtxToVal(mtxTmp, -n.x, BASE_FUNC_COUNT);
//                addSmallMatrToBigMatr(matrBig2, mtxTmp, 0, 1);
//
//                copyMtx(mtxTmp, mtxKI, BASE_FUNC_COUNT);
//                multMtxToVal(mtxTmp, -n.y, BASE_FUNC_COUNT);
//                addSmallMatrToBigMatr(matrBig2, mtxTmp, 0, 2);
//            }
//            else {
                copyMtx(mtxTmp, mtxII, BASE_FUNC_COUNT);
                multMtxToVal(mtxTmp, -n.x, BASE_FUNC_COUNT);
                addSmallMatrToBigMatr(matrBig, mtxTmp, 0, 1);

                copyMtx(mtxTmp, mtxII, BASE_FUNC_COUNT);
                multMtxToVal(mtxTmp, -n.y, BASE_FUNC_COUNT);
                addSmallMatrToBigMatr(matrBig, mtxTmp, 0, 2);


                copyMtx(mtxTmp, mtxKI, BASE_FUNC_COUNT);
                multMtxToVal(mtxTmp, -n.x, BASE_FUNC_COUNT);
                multMtxToVal(mtxTmp, k1, BASE_FUNC_COUNT);
                addSmallMatrToBigMatr(matrBig2, mtxTmp, 1, 0);

                copyMtx(mtxTmp, mtxKI, BASE_FUNC_COUNT);
                multMtxToVal(mtxTmp, -n.y, BASE_FUNC_COUNT);
                multMtxToVal(mtxTmp, k1, BASE_FUNC_COUNT);
                addSmallMatrToBigMatr(matrBig2, mtxTmp, 2, 0);
//            }

            solverMtx->addMatrElement(c2, c2, matrBig);
            solverMtx->addMatrElement(c2, c1, matrBig2);
        }
    }

}

void HEAT_DG_IMPLICIT::calcRHS()
{

    for (int iCell = 0; iCell < grid.cCount; iCell++) {
        Region &reg = getRegion(iCell);


        memset(tmpArr, 0, sizeof(double)*MATR_DIM);
        multMtxToVec(tmpArr, matrA[iCell], u[iCell], BASE_FUNC_COUNT);
        for (int i = 0; i < BASE_FUNC_COUNT; i++) {
            tmpArr[i] /= cTau[iCell];
        }
        solverMtx->addRightElement(iCell, tmpArr);

//        // CONVECTION
//        memset(tmpArr, 0, sizeof(double)*MATR_DIM);
//        for (int m = 0; m < BASE_FUNC_COUNT; m++) {
//            for (int iGP = 0; iGP < GP_CELL_COUNT; iGP++) {
//                tmpArr[m] += cellGW[iCell][iGP] * (getField(FIELD_U, iCell, cellGP[iCell][iGP])) *
//                          (reg.par.u*getDfDx(m, iCell, cellGP[iCell][iGP]) + reg.par.v*getDfDy(m, iCell, cellGP[iCell][iGP]));
//            }
//            tmpArr[m] *= cellJ[iCell];
//        }
//        solverMtx->addRightElement(iCell, tmpArr);

    }

    for (int iEdge = 0; iEdge < grid.eCount; iEdge++) {
        Edge &edge = grid.edges[iEdge];
        Vector &n = grid.edges[iEdge].n;
        int c1 = edge.c1;
        int c2 = edge.c2;

        Region   &reg1 = getRegion(c1);
        Material &mat1 = getMaterial(c1);
        double k1 = mat1.K;

        if (c2 < 0) {
            memset(tmpArr, 0, sizeof(double)*MATR_DIM);

            for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
                Param p1, p2;
                double fU, fQx, fQy;
                getFields(fU, fQx, fQy, c1, edgeGP[iEdge][iGP]);

                p1.T     = fU;
                p1.Qt[0] = fQx;
                p1.Qt[1] = fQy;


                edge.bnd->run(iEdge, p1, p2);

                if (n.x+n.y >= 0) {
                    fU  = p1.T;
                    fQx = p2.Qt[0];
                    fQy = p2.Qt[1];
                }
                else {
                    fU  = p2.T;
                    fQx = p1.Qt[0];
                    fQy = p1.Qt[1];
                }

                for (int m = 0; m < BASE_FUNC_COUNT; m++) {
                    double fw = getF(m, c1, edgeGP[iEdge][iGP])*edgeGW[iEdge][iGP];
                    tmpArr[ m ]                    -= (fQx*n.x+fQy*n.y)*fw;
                    tmpArr[ m + BASE_FUNC_COUNT]   -= k1*fU*n.x*fw;
                    tmpArr[ m + BASE_FUNC_COUNT*2] -= k1*fU*n.y*fw;
                }
            }

            for (int i = 0; i < MATR_DIM; i++) {
                tmpArr[i] *= edgeJ[iEdge];
            }

            solverMtx->addRightElement(c1, tmpArr);
        }


//        // CONVECTION
//
//        if (c2 >= 0) {
//            // C1
//            memset(tmpArr, 0, sizeof(double)*MATR_DIM);
//            for (int m = 0; m < BASE_FUNC_COUNT; m++) {
//                for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
//                    double vn = reg1.par.u*n.x+reg1.par.v*n.y;
//                    double flux = ((vn >=0) ? getField(FIELD_U, c1, edgeGP[iEdge][iGP]) : getField(FIELD_U, c2, edgeGP[iEdge][iGP])) * vn;
//                    tmpArr[m] -= edgeGW[iEdge][iGP]*flux*getF(m, c1, edgeGP[iEdge][iGP]);
//                }
//                tmpArr[m] *= edgeJ[iEdge];
//            }
//            solverMtx->addRightElement(c1, tmpArr);
//
//            // C2
//            Region &reg2 = getRegion(c2);
//            memset(tmpArr, 0, sizeof(double)*MATR_DIM);
//            for (int m = 0; m < BASE_FUNC_COUNT; m++) {
//                for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
//                    double vn = reg2.par.u*n.x+reg2.par.v*n.y;
//                    double flux = ((vn >=0) ? getField(FIELD_U, c1, edgeGP[iEdge][iGP]) : getField(FIELD_U, c2, edgeGP[iEdge][iGP])) * vn;
//                    tmpArr[m] += edgeGW[iEdge][iGP]*flux*getF(m, c2, edgeGP[iEdge][iGP]);
//                }
//                tmpArr[m] *= edgeJ[iEdge];
//            }
//            solverMtx->addRightElement(c2, tmpArr);
//        }
//        else {
//            // C1
//            memset(tmpArr, 0, sizeof(double)*MATR_DIM);
//            for (int m = 0; m < BASE_FUNC_COUNT; m++) {
//                for (int iGP = 0; iGP < GP_EDGE_COUNT; iGP++) {
//                    double fU1 = getField(FIELD_U, c1, edgeGP[iEdge][iGP]);
//                    double fU2 = phi_exac(time, edgeGP[iEdge][iGP].x, edgeGP[iEdge][iGP].y);
//
//                    double vn = reg1.par.u*n.x+reg1.par.v*n.y;
//                    double flux = ((vn >=0) ? fU1 : fU2) * vn;
//                    tmpArr[m] += edgeGW[iEdge][iGP]*flux*getF(m, c1, edgeGP[iEdge][iGP]);
//                }
//                tmpArr[m] *= edgeJ[iEdge];
//            }
//            solverMtx->addRightElement(c1, tmpArr);
//        }
    }

}


void HEAT_DG_IMPLICIT::run()
{
    int solverErr = 0;
    time = 0.0;
    int step = 0;
    long totalCalcTime = 0;
    while (time < TMAX && step < STEP_MAX) {
        long timeStart, timeEnd;
        timeStart = clock();

        solverMtx->zero();

        /* Заполняем правую часть */
        calcRHS();

        /* Заполняем элементы матрицы */
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

        if (solverErr == MatrixSolver::RESULT_OK) {

            step++;
            time += TAU;

            for (int cellIndex = 0, ind = 0; cellIndex < grid.cCount; cellIndex++)
            {
                for (int iFld = 0; iFld < FIELD_COUNT; iFld++) {
                    for (int iF = 0; iF < BASE_FUNC_COUNT; iF++) {
                        fields[iFld][cellIndex][iF] = solverMtx->x[ind++];
                    }
                }
            }

            timeEnd = clock();
            totalCalcTime += (timeEnd - timeStart);
            if (step % FILE_SAVE_STEP == 0)
            {
                save(step);
            }
            if (step % PRINT_STEP == 0)
            {
                log((char*)"step: %6d  time step: %.16f\tmax iter: %5d\ttime: %6d ms\ttotal calc time: %ld\n", step, time, maxIter, timeEnd - timeStart, totalCalcTime);
            }
        }
        else {
            if (solverErr & MatrixSolver::RESULT_ERR_CONVERG) {
                log((char*)"Solver error: residual condition not considered.\n");
            }
            if (solverErr & MatrixSolver::RESULT_ERR_MAX_ITER) {
                log((char*)"Solver error: max iterations done.\n");
            }
            solverErr = 0;
        }
    }
}

void HEAT_DG_IMPLICIT::done()
{
    //memFree();
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
        fprintf(fp, "%16.8f ", p.T*T_);
        if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
    }


    fprintf(fp, "SCALARS TAU float 1\nLOOKUP_TABLE default\n");
    for (int i = 0; i < grid.cCount; i++)
    {
        fprintf(fp, "%16.8f ", cTau[i]*TIME_);
        if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
    }

    fprintf(fp, "SCALARS CellID float 1\nLOOKUP_TABLE default\n");
    for (int i = 0; i < grid.cCount; i++)
    {
        fprintf(fp, "%d ", i);
        if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
    }

    fclose(fp);
    printf("File '%s' saved...\n", fName);

}



