#include "global.h"

#ifndef _FEM_RKDG_H_
#define _FEM_RKDG_H_

#include <method.h>


class FEM_RKDG: public Method {
public:
	virtual void init(char * xmlFileName);
	virtual void run();
	virtual void done();
protected:
	Region & getRegionByCellType(int type);

	Region   &	getRegion	(int iCell);
	Material &	getMaterial	(int iCell);
	
	/**
	 *	Преобразование примитивных переменных в консервативные
	 */
	void convertParToCons(int iCell, Param & par);
	
	/**
	 *	Преобразование консервативных переменных в примитивные
	 */
	void convertConsToPar(int iCell, Point pt, Param & par);
	
	/**
	 *	Вычисление параметров справа и слева от границы ячейки
	 */
	void reconstruct(int iEdge, Point pt, Param& pL, Param& pR);

	/**
	 *	Вычисление параметров с внешней стороны от границы ячейки согласно граничным условиям
	 */
	void boundaryCond(int iEdge, Point pt, Param& pL, Param& pR);

	/**
	 *	Вычисление шага по времени по значению CFL, если значение TAU из XML 
	 *	меньше вычисленного, то используется значение, заданное в XML
	 */
	void calcTimeStep();

	/**
	 *	Запись значений газодинамических параметров в файл
	 */
	void save(int);

	/**
	 *	Вычисление численного потока
	 */
	void calcFlux(double& fr, double& fu, double& fv, double& fe, Param pL, Param pR, Vector n, double GAM);

	void calcLimiters();

	void calcConvectionVol();

	void calcConvectionSurf();
	
	void calcDiffusionVol();

	void calcDiffusionSurf();
	
	void calcNewFields();

	void calcGP();

	void calcMatr();
	
	/**
	 *	Получение значения вектора базовых функций
	 */
	VECTOR getF(int iCell, Point pt);

	VECTOR getDFDX(int iCell, Point pt);

	VECTOR getDFDY(int iCell, Point pt);

	inline double getRO(int iCell, Point pt) { return VECTOR::SCALAR_PROD(ro[iCell], getF(iCell, pt)); }
	inline double getRU(int iCell, Point pt) { return VECTOR::SCALAR_PROD(ru[iCell], getF(iCell, pt)); }
	inline double getRV(int iCell, Point pt) { return VECTOR::SCALAR_PROD(rv[iCell], getF(iCell, pt)); }
	inline double getRE(int iCell, Point pt) { return VECTOR::SCALAR_PROD(re[iCell], getF(iCell, pt)); }

private:
	double TMAX;
	double TAU;
	double CFL;
	int FILE_SAVE_STEP;
	int PRINT_STEP;
	int EDGE_GP_COUNT;
	int CELL_GP_COUNT;
	int FUNC_COUNT;

	int				matCount;
	int				regCount;
	int				bCount;
	Material	*	materials;
	Region		*	regions;
	Boundary	*	boundaries;

	//! консервативные переменные на текущем временном слое (коэф. разлож. по базисным функциям)
	VECTOR * ro;			 
	VECTOR * ru;			
	VECTOR * rv;			
	VECTOR * re;			

	//! консервативные переменные на предыдущем временном слое (коэф. разлож. по базисным функциям)
	VECTOR * ro_old;
	VECTOR * ru_old;
	VECTOR * rv_old;
	VECTOR * re_old;

	//! правые части системы уравнений разностной схемы
	VECTOR * ro_int;
	VECTOR * ru_int;
	VECTOR * rv_int;
	VECTOR * re_int;
	
	Point ** edgeGP;
	Point ** cellGP;

	double * edgeGW;
	double * cellGW;

	double * edgeJ;
	double * cellJ;

	MATRIX * cellA;
	MATRIX * cellInvA;

};

#endif