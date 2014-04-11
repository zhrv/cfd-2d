#ifndef _FVM_TVD_IMPLICIT_H_
#define _FVM_TVD_IMPLICIT_H_

#include <method.h>
#include <algorithm>

//TODO:  ѕереписать! Ќе хватает пам€ти. ћатрица очень разр€жена. »спользовать CSR способ хранени€.
struct ConservativeSystem: public MATRIX
{
	ConservativeSystem(int an, int dimension): MATRIX(dimension*an) { conservativeDim = an; dim = dimension; conservativeCount = dim*an*dim*an; }
	ConservativeSystem(const ConservativeSystem& m) : MATRIX(m) 
	{ 
		conservativeDim = m.conservativeDim; conservativeCount = m.conservativeCount; dim = m.dim;
		for (int i = 0; i < n; ++i) rightElem[i] = m.rightElem[i]; 
	}
	~ConservativeSystem() { if (rightElem) delete [] rightElem; }
	void clear()
	{
		memset(elem, 0x0, sizeof(double)*count);
		memset(rightElem, 0x0, sizeof(double)*n);
	}
	void setMatrElement(int i, int j, double **mtxDim) 
	{
		for (int ii = 0; ii < dim; ++ii)
			for (int jj = 0; jj < dim; ++jj)
				elem[ii+i*dim, ii+j*dim] = mtxDim[ii][jj];
	}
	void addMatrElement(int i, int j, double **mtxDim) 
	{
		for (int ii = 0; ii < dim; ++ii)
			for (int jj = 0; jj < dim; ++jj)
				elem[ii+i*dim, ii+j*dim] += mtxDim[ii][jj];
	}
	void setRightElement(int i, double *rightDim) 
	{
		for (int ii = 0; ii < dim; ii++)
			rightElem[ii+i*dim] = rightDim[ii];
	}
	void addRightElement(int i, double *rightDim) 
	{
		for (int ii = 0; ii < dim; ii++)
			rightElem[ii+i*dim] += rightDim[ii];
	}
	//! решает систему методом «ейдел€. решение записываетс€ в вектор rightElem или result,
	//! если result передан как аргумент, то следует заполнить result первым приближением.
	void solveZeidel(double eps = 1.0e-5, int maxIt = 1000, double *result = 0)
	{
		bool resultNull = result == 0? true : false;
		if (resultNull) { result = new double [n]; memset(result, 0, sizeof(double)*n); }
		double *previous = new double [n];
		int it = 0;
		while (1)
		{
			for (int i = 0; i < n; ++i) previous[i] = result[i];
			for (int i = 0; i < n; i++)
			{
				double var = 0;
				for (int j = 0; j < i; j++)
					var += (elem[i*n + j] * result[j]);
				for (int j = i + 1; j < n; j++)
					var += (elem[i*n + j] * previous[j]);
				result[i] = (rightElem[i] - var) / elem[i*n + i];
			}
			++it;
			bool next = false;
			if (it < maxIt)
			{
				for (int i = 0; i < n; ++i)
				{
					if (abs(result[i] - previous[i]) >= eps) 
					{
						next = true;
						break;
					}
				}
			}
			if (!next) break;
		}
		if (resultNull)
		{
			for (int i = 0; i < n; ++i) rightElem[i] = result[i];
		}
		delete [] previous;
		if (resultNull) delete [] result;
	}

	int		conservativeDim;
	int		conservativeCount;
	int		dim;
	double *rightElem;
};

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
	void calcFlux(double& fr, double& fu, double& fv, double& fe, Param pL, Param pR, Vector n, double GAM);

private:
	double **allocMtx4();
	void freeMtx4(double **mtx4);
	/**
	 * «аполн€ет в mtx4 собственные значени€ матрицы якоби функции потока.
	 */
	void eigenValues(double **dst4, double c, double u, double nx, double v, double ny);
	void rightEigenVector(double **dst4, double c, double u, double nx, double v, double ny, double H);
	void leftEigenVector(double **dst4, double c, double GAM, double u, double nx, double v, double ny);
	void mult(double **dst4, double **srcA4, double **srcB4);
	/**
	 * «аполн€ет в mtx4 A+.
	*/
	void calcAP(double **dst4, double **rightEgnVecl4, double **egnVal4, double **leftEgnVecl4);
	/**
	 * «аполн€ет в mtx4 A-.
	*/
	void calcAM(double **dst4, double **rightEgnVecl4, double **egnVal4, double **leftEgnVecl4);
	void calcRoeAverage(double& fr, double& fu, double& fv, double& fe, Param pL, Param pR, Vector n, double GAM);
	void reconstruct(int iCell, Param& cell, Param neighbor[3]);
private:
	//массив размерность 3 x grid.cCount. по которому можно определить, кака€ 3-угольна€ €чейка, какие нормера ребер имеет.
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
};

#endif