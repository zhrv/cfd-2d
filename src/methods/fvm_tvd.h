#ifndef _FVM_TVD_H_
#define _FVM_TVD_H_

#include <method.h>
#include "k_eps_model.h"

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
	 *	Вычисление параметров справа и слева от границы ячейки
	 */
	void reconstruct(int iFace, Param& pL, Param& pR);

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


	//_нач
	void calcGrad(Vector *gradU, Vector *gradV);
	void calcTensor(double l, double m, Vector *gradU, Vector *gradV, double *Txx, double *Tyy, double *Txy);
	void calcDiffFlux(double& fu_diff, double& fv_diff, double& fe_diff, Vector n, int c1, int c2);
	//_кон

private:
	double TMAX;
	double TAU;
	double CFL;
	int FILE_SAVE_STEP;
	int PRINT_STEP;

	int				matCount;
	int				regCount;
	int				bCount;
	Material	*	materials;
	Region		*	regions;
	Boundary	*	boundaries;

	ViscosityModel * viscosityModel;

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
	
	Vector *gradU, *gradV;
	double *Txx, *Tyy, *Txy;

	// TODO: Вынести
	double lambda, mu;
};

#endif