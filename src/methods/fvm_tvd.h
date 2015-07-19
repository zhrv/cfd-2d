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
	 *	�������������� ����������� ���������� � ��������������
	 */
	void convertParToCons(int iCell, Param & par);
	
	/**
	 *	�������������� �������������� ���������� � �����������
	 */
	void convertConsToPar(int iCell, Param & par);
	
	/**
	 *	���������� ���������� ������ � ����� �� ������� ������
	 */
	void reconstruct(int iFace, Param& pL, Param& pR, Point p);

	/**
	 *	���������� ���������� � ������� ������� �� ������� ������ �������� ��������� ��������
	 */
	void boundaryCond(int iFace, Param& pL, Param& pR);

	/**
	 *	���������� ���� �� ������� �� �������� CFL, ���� �������� TAU �� XML 
	 *	������ ������������, �� ������������ ��������, �������� � XML
	 */
	void calcTimeStep();

	/**
	 *	������ �������� ���������������� ���������� � ����
	 */
	void save(int);

	/**
	 *	���������� ���������� ������
	 */
	void calcFlux(double& fr, double& fu, double& fv, double& fe, Param pL, Param pR, Vector n, double GAM);

	/**
	 *	���������� ���������� � �������
	 */
	void calcGrad();

	void setCellFlagLim(int iCell)	{ grid.cells[iCell].flag |= CELL_FLAG_LIM; }
	bool cellIsLim(int iCell)		{ return (grid.cells[iCell].flag & CELL_FLAG_LIM) > 0; }

	void remediateLimCells();
	Region & getRegionByName(char* name);
	Region & getRegion(char * name);

private:
	double TMAX;
	double TAU;
	double CFL;
	int STEP_MAX;
	int FILE_SAVE_STEP;
	int PRINT_STEP;

	bool			STEADY;		// false - �������������� �������, true - ������������� �������.
	double		*	cTau;		// ��������� ��� �� ������� � ������

	int				matCount;
	int				regCount;
	int				bCount;
	Material	*	materials;
	Region		*	regions;
	CFDBoundaries	boundaries;

	//! �������������� ���������� �� ������� ��������� ����
	double * ro;			 
	double * ru;			
	double * rv;			
	double * re;			

	//! �������������� ���������� �� ���������� ��������� ����
	double * ro_old;
	double * ru_old;
	double * rv_old;
	double * re_old;

	//! ������ ����� ������� ��������� ���������� �����
	double * ro_int;
	double * ru_int;
	double * rv_int;
	double * re_int;
	//! ���������
	Vector *gradR;
	Vector *gradP;
	Vector *gradU;
	Vector *gradV;
	
	//! ������
	double limitRmin;
	double limitRmax;
	double limitPmin;
	double limitPmax;
	double limitUmax;

};

#endif