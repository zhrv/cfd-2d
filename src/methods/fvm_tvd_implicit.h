#ifndef _FVM_TVD_IMPLICIT_H_
#define _FVM_TVD_IMPLICIT_H_

#include <method.h>
#include "MatrixSolver.h"
#include <algorithm>

class FVM_TVD_IMPLICIT: public Method
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

private:
	double **allocMtx4();
	void freeMtx4(double **mtx4);
	void multMtx4(double **dst4, double **srcA4, double **srcB4);
	void clearMtx4(double **mtx4);
	/**
	 * ��������� � mtx4 ����������� �������� ������� ����� ������� ������.
	 */
	void eigenValues(double **dst4, double c, double u, double nx, double v, double ny);
	void rightEigenVector(double **dst4, double c, double u, double nx, double v, double ny, double H);
	void leftEigenVector(double **dst4, double c, double GAM, double u, double nx, double v, double ny);
	/**
	 * ��������� � mtx4 A+.
	*/
	void calcAP(double **dst4, double **rightEgnVecl4, double **egnVal4, double **leftEgnVecl4);
	/**
	 * ��������� � mtx4 A-.
	*/
	void calcAM(double **dst4, double **rightEgnVecl4, double **egnVal4, double **leftEgnVecl4);
	void calcRoeAverage(Param& average, Param pL, Param pR, double GAM);
	void reconstruct(int iCell, Param& cell, Param neighbor[3]);
private:
	//������ ����������� 3 x grid.cCount. �� �������� ����� ����������, ����� 3-�������� ������, ����� ������� ����� �����.
	int **cellsEdges;

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

	//! �������������� ���������� �� ������� ��������� ����
	double * ro;			 
	double * ru;			
	double * rv;			
	double * re;			
/*
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
*/
};

#endif