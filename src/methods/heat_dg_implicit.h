#pragma once
#include "method.h"
#include "MatrixSolver.h"

class HEAT_DG_IMPLICIT : public Method
{

protected:
    const static int BASE_FUNC_COUNT = 3;

    double **u;
    double **qx;
    double **qy;

public:
    void init(char * xmlFileName);
    void run();
    void done();

    virtual void getFields(double &fRO, double &fRU, double &fRV, double &fRE, int iCell, Point p);
    virtual void getFields(double &fRO, double &fRU, double &fRV, double &fRE, int iCell, double x, double y);

    virtual double getField(int fld, int iCell, Point p);
    virtual double getField(int fld, int iCell, double x, double y);

    virtual double getF(int id, int iCell, Point p);
    virtual double getF(int id, int iCell, double x, double y);

    void memAlloc();
    void memFree();

    void calcMassMatr(); //!< вычисляем матрицу масс
    void calcGaussPar(); //!< вычисляем узлы и коэффициенты квадратур

    void calcTimeStep();

    Region & getRegionByCellType(int type);
    Region & getRegionByName(char* name);
    Region & getRegion(char * name);

    Region   &	getRegion(int iCell);
    Material &	getMaterial(int iCell);

    inline void convertParToCons(int iCell, Param & par);

    inline void convertConsToPar(int iCell, Param & par);


    inline double getDfDx(int id, int iCell, Point p);
    inline double getDfDx(int id, int iCell, double x, double y);

    inline double getDfDy(int id, int iCell, Point p);
    inline double getDfDy(int id, int iCell, double x, double y);

    void save(int step);

//    inline void setCellFlagLim(int iCell){ grid.cells[iCell].flag |= CELL_FLAG_LIM; }
//    inline bool cellIsLim(int iCell)		{ return (grid.cells[iCell].flag & CELL_FLAG_LIM) > 0; }
//    int getLimitedCellsCount();
//    void remediateLimCells();

    void incCFL();
    void decCFL();

    void calcIntegral();		//!< вычисляем интеграл от(dF / dU)*deltaU*dFi / dx
    void calcMatrWithTau();		//!< вычисляем матрицы перед производной по времени
    void calcMatrFlux();		//!< Вычисляем потоковые величины
    void calcRHS();				//!< Вычисляем столбец правых членов

//    void calcLiftForce();

//    double** allocMtx4();
//    void freeMtx4(double **mtx4);
//    void multMtx4(double **dst4, double **srcA4, double **srcB4);
//    void clearMtx4(double **mtx4);
    void multMtxToVal(double **dst, double x, int N);
    void multMtxToVec(double *dst, double **mtx, double *vec, int N);
    void fillMtx(double** dst, double x, int N);

//    void eigenValues(double** dst4, double c, double u, double nx, double v, double ny);
//    void rightEigenVector(double **dst4, double c, double u, double nx, double v, double ny, double H);
//    void leftEigenVector(double **dst4, double c, double GAM, double u, double nx, double v, double ny);
//    void calcAP(double **dst4, double **rightEgnVecl4, double **egnVal4, double **leftEgnVecl4);
//    void calcAM(double **dst4, double **rightEgnVecl4, double **egnVal4, double **leftEgnVecl4);
//    void calcRoeAverage(Param& average, Param pL, Param pR, double GAM, Vector n);
//    void _calcA(double **dst4, double **rightEgnVecl4, double **egnVal4, double **leftEgnVecl4);
//    void calcA(double **dst4, double c, double GAM, double u, double nx, double v, double ny, double H);
//    void calcAx(double **dst4, double c, double GAM, double u, double v, double H);
//    void calcAy(double **dst4, double c, double GAM, double u, double v, double H);
//    void calcAx_(double **dst4, Param par, double GAM);
//    void calcAy_(double **dst4, Param par, double GAM);
//    void consToPar(double fRO, double fRU, double fRV, double fRE, Param& par);
//
//    void calcFlux(double& fr, double& fu, double& fv, double& fe, Param pL, Param pR, Vector n, double GAM);
//    void boundaryCond(int iEdge, Param& pL, Param& pR);

    void addSmallMatrToBigMatr(double **mB, double **mS, int i, int j);

private:
    double			TMAX;
    int				STEP_MAX;
    double			TAU;
    double			CFL;
    double			scaleCFL;
    double			maxCFL;
    int				stepCFL;
    double			maxLimCells;
    int				FILE_SAVE_STEP;
    int				PRINT_STEP;

    double			TAU_MIN;

    bool			STEADY;	// false - нестационарное течение, true - стационанрное течение.
    double			*cTau;  // локальный шаг по времени в ячейке.
    bool			SMOOTHING;
    double			SMOOTHING_PAR;
    int				FLUX;

    double			SIGMA;

    double			SOLVER_EPS		= 1.e-7;
    int				SOLVER_ITER		= 50;

    int				matCount;
    int				regCount;
    int				bCount;
    Material	   *materials;
    Region		   *regions;
    CFDBoundaries	boundaries;


    double			***fields;

    double			*tmpArr;
    double			*tmpArr1;
    double			*tmpArr2;
    int				*tmpArrInt;

    //! градиенты.
    Vector			*gradR;
    Vector			*gradP;
    Vector			*gradU;
    Vector			*gradV;

    //! лимиты
    double			limitRmin;
    double			limitRmax;
    double			limitPmin;
    double			limitPmax;
    double			limitUmax;

    int				limCells;

    //! подъемная сила.
    double			Fx;
    double			Fy;

    //! узлы и коэффициенты квадратур
    Point			**cellGP;
    Point			**edgeGP;
    double			**cellGW;
    double			*cellJ;
    double			**edgeGW;
    double			*edgeJ;

    // матрица масс
    double			***matrA;
    double			***matrBx;
    double			***matrBy;
    double			***matrInvA;

    MatrixSolver	*solverMtx;

    double			**matrBig;
    double			**matrSmall;
    double			**matrBig2;
    double			**matrSmall1;
    double			**matrSmall2;

    double			*tmpCFL;
    double			pCFL;


    // параметры обезразмеривания
    double			L_		= 1.0;
    double			U_		= 1.0;
    double			R_		= 1.0;
    double			P_		= 1.0;
    double			T_		= 1.0;
    double			E_		= 1.0;
    double			CV_ = 1.0;
    double			MU_ = 1.0;
    double			KP_ = 1.0;
    double			TIME_ = 1.0;

protected:
    const static int FLUX_GODUNOV = 0;
    const static int FLUX_LAX = 1;


    const static int GP_CELL_COUNT = 3;
    const static int GP_EDGE_COUNT = 2;

    const static int FIELD_COUNT = 3;
    const static int FIELD_U  = 0;
    const static int FIELD_QX = 1;
    const static int FIELD_QY = 1;

    const static int MATR_DIM = FIELD_COUNT * BASE_FUNC_COUNT;

    inline double getGAM(int iCell) { return getMaterial(iCell).getGamma(); } // TODO: сделать
};
