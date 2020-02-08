#ifndef _FVM_TVD_H_
#define _FVM_TVD_H_

#include <method.h>
#include "bnd_cond.h"

class FVM_TVD: public Method
{
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
	void convertConsToPar(int iCell, Param & par);

    /**
     *  Вычисление одного шага по времени
     */
	void singleTimeStep();

	/**
	 *	Вычисление параметров справа и слева от границы ячейки
	 */
	void reconstruct(int iFace, Param& pL, Param& pR, Point p);

	/**
	 *	Вычисление параметров с внешней стороны от границы ячейки согласно граничным условиям
	 */
	void boundaryCond(int iFace, Param& pL, Param& pR);

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

	/**
	 *	Вычисление градиентов в ячейках
	 */
	void calcGrad();

	void setCellFlagLim(int iCell)	{ grid.cells[iCell].flag |= CELL_FLAG_LIM; }
	bool cellIsLim(int iCell)		{ return (grid.cells[iCell].flag & CELL_FLAG_LIM) > 0; }

	void remediateLimCells();
	Region & getRegionByName(char* name);
	Region & getRegion(char * name);

    void procExchangeFields();
    void procExchangeGrads();
private:
	double TMAX;
	double TAU;
	double CFL;
	int STEP_MAX;
	int FILE_SAVE_STEP;
	int PRINT_STEP;

	bool			STEADY;		// false - нестационарное течение, true - стационанрное течение.
	double		*	cTau;		// локальный шаг по времени в ячейке

	int				matCount;
	int				regCount;
	int				bCount;
	Material	*	materials;
	Region		*	regions;
	CFDBoundaries	boundaries;

	//! консервативные переменные на текущем временном слое
	double * ro;			 
	double * ru;			
	double * rv;			
	double * re;			

	//! консервативные переменные на предыдущем временном слое
	double * ro_old;
	double * ru_old;
	double * rv_old;
	double * re_old;

	//! правые части системы уравнений разностной схемы
	double * ro_int;
	double * ru_int;
	double * rv_int;
	double * re_int;
	//! градиенты
	Vector *gradR;
	Vector *gradP;
	Vector *gradU;
    Vector *gradV;
    Vector *gradT;

	//! лимиты
	double limitRmin;
	double limitRmax;
	double limitPmin;
	double limitPmax;
	double limitUmax;

};

#endif