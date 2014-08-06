#pragma once
#include "method.h"
#include "MatrixSolver.h"


class FEM_DG_IMPLICIT :	public Method
{
public:
	FEM_DG_IMPLICIT();
	~FEM_DG_IMPLICIT();
public:
	void init(char * xmlFileName);
	void run();
	void done();

private:
	void memAlloc();
	void memFree();

	void calcTimeStep();
	
	Region & getRegionByCellType(int type);

	Region   &	getRegion(int iCell);
	Material &	getMaterial(int iCell);

	/**
	*	Преобразование примитивных переменных в консервативные
	*/
	void convertParToCons(int iCell, Param & par);

	/**
	*	Преобразование консервативных переменных в примитивные
	*/
	void convertConsToPar(int iCell, Param & par);

	double getField(int fld, int iCell, Point p);
	double getField(int fld, int iCell, double x, double y);

	double getF(int id, int iCell, Point p);
	double getF(int id, int iCell, double x, double y);

	double getDfDx(int id, int iCell, Point p);
	double getDfDx(int id, int iCell, double x, double y);

	double getDfDy(int id, int iCell, Point p);
	double getDfDy(int id, int iCell, double x, double y);

	void save(int step);

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


	int				matCount;
	int				regCount;
	int				bCount;
	Material	   *materials;
	Region		   *regions;
	Boundary	   *boundaries;

	//! консервативные переменные на текущем временном слое.
	double			**ro;
	double			**ru;
	double			**rv;
	double			**re;

	double			***fields;

	double			*tmpArr;
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
	double			***matrInvA;

	MatrixSolver	*solverMtx;

protected:
	const static int FLUX_GODUNOV = 0;
	const static int FLUX_LAX = 1;

	const static int BASE_FUNC_COUNT = 3;
	const static int GP_CELL_COUNT = 3;
	const static int GP_EDGE_COUNT = 2;

	const static int FIELD_RO = 0;
	const static int FIELD_RU = 1;
	const static int FIELD_RV = 2;
	const static int FIELD_RE = 3;
};

