#ifndef _FVM_TVD_H_
#define _FVM_TVD_H_

#include <method.h>
#include "viscosity_model.h"

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
	 *	ѕреобразование примитивных переменных в консервативные
	 */
	void convertParToCons(int iCell, Param & par);
	
	/**
	 *	ѕреобразование консервативных переменных в примитивные
	 */
	void convertConsToPar(int iCell, Param & par);
	
	/**
	 *	¬ычисление параметров справа и слева от границы €чейки
	 */
	void reconstruct(int iFace, Param& pL, Param& pR);
	void reconstruct(int iFace, Param& pL, Param& pR, Point p);

	/**
	 *	¬ычисление параметров с внешней стороны от границы €чейки согласно граничным услови€м
	 */
	void boundaryCond(int iFace, Param& pL, Param& pR);

	/**
	 *	¬ычисление шага по времени по значению CFL, если значение TAU из XML 
	 *	меньше вычисленного, то используетс€ значение, заданное в XML
	 */
	void calcTimeStep();

	/**
	 *	«апись значений газодинамических параметров в файл
	 */
	void save(int);

	/**
	 *	¬ычисление численного потока
	 */
	void calcFlux(double& fr, double& fu, double& fv, double& fe, double& ro_m, double& u_m, double& v_m, Param pL, Param pR, Vector n, double GAM);


	//_нач
	void calcGrad(Vector *gradR, Vector *gradP, Vector *gradU, Vector *gradV);
	void calcTensor(double l, double m, Vector *gradU, Vector *gradV, double *Txx, double *Tyy, double *Txy);
	void calcDiffFlux(double& fu_diff, double& fv_diff, double& fe_diff, Param pL, Param pR, Vector n, int c1, int c2);
	void chooseTurbulenceModel( const char * turbModelStr );
	//_кон

	void setCellFlagLim(int iCell)	{ grid.cells[iCell].flag |= CELL_FLAG_LIM; }
	bool cellIsLim(int iCell)		{ return (grid.cells[iCell].flag & CELL_FLAG_LIM) > 0; }

	void remediateLimCells();

private:
	double TMAX;
	double TAU;
	double CFL;
	int STEP_MAX;
	int FILE_SAVE_STEP;
	int PRINT_STEP;

	bool			STEADY;		// false - нестационарное течение, true - стационанрное течение.
	double		*	cTau;		// локальный шаг по времени в €чейке

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

	//! средние значени€ консервативных переменных на рЄбрах
	double * ro_m;
	double * u_m;
	double * v_m;
	
	//! градиенты
	Vector *gradR, *gradP, *gradU, *gradV;

	//! тензоры
	double *Txx, *Tyy, *Txy;

	// TODO: ¬ынести
	double lambda, mu;

	//! лимиты
	double limitRmin;
	double limitRmax;
	double limitPmin;
	double limitPmax;
	double limitUmax;
};

#endif