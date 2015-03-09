#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
	#include <direct.h>
	#define MK_DIR(name) _mkdir("mesh");
#else
	#include <sys/stat.h>
	#define MK_DIR(name) mkdir(name, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif



// ! универсальная газовая постоянная
//static double gR = 8.314472;	// м2 кг с-2 К-1 Моль-1

/**
 *	Точка на плоскости
 */
struct Point
{
	double x;
	double y;
};


/**
 *	Двумерный вектор
 */
typedef Point Vector;

/**
 *	Вектор размерности n
 */
struct VECTOR
{
	VECTOR(int an = 0): n(an) { if (n) { elem = new double[n]; memset(elem, 0, sizeof(double)*n); } else elem = NULL; }
	VECTOR(const VECTOR& v): n(v.n) { elem = new double[n]; memcpy(elem, v.elem, n*sizeof(double)); }
	~VECTOR() { if (elem) delete[] elem; n = 0; }

	void init(int an) { n = an; if (elem) delete[] elem; elem = new double[n]; }

	double&		operator []		(int i)					{ return elem[i]; }
	
	double&		operator ()		(int i)					{ return elem[i]; }
	
	void		operator =		(const VECTOR& v)		{	
		if (elem) delete[] elem; 
		n = v.n; 
		elem = new double[n]; 
		memcpy(elem, v.elem, n*sizeof(double)); 
	}

	void		operator =		(const double& x)		{ for (int i = 0; i < n; i++) elem[i] = x; }

	void		operator +=		(const VECTOR& v)		{ for (int i = 0; i < n; i++) elem[i] += v.elem[i]; }	

	void		operator +=		(const double& x)		{ for (int i = 0; i < n; i++) elem[i] += x;	}	

	void		operator -=		(const VECTOR& v)		{ for (int i = 0; i < n; i++) elem[i] -= v.elem[i];	}		
	
	void		operator -=		(const double& x)		{ for (int i = 0; i < n; i++) elem[i] -= x; }

	void		operator *=		(const double& x)		{ for (int i = 0; i < n; i++) elem[i] *= x;	}		

	static double SCALAR_PROD(const VECTOR& v1, const VECTOR& v2) {
		double s = 0;
		for (int i = 0; i < v1.n; i++)
		{
			s += v1.elem[i]*v2.elem[i];
		}
		return s;
	}

	double *	elem;
	int			n;
};



/**
 * Матрица размерности nxn
 */
struct MATRIX
{
	MATRIX(int an): n(an), count(an*an) { elem = new double[n*n]; memset(elem, 0, sizeof(double)*count); }
	MATRIX(const MATRIX& m): n(m.n), count(m.n*m.n) { if (elem) delete[] elem; elem = new double[count]; memcpy(elem, m.elem, count*sizeof(double)); }
	~MATRIX() { if (elem) delete[] elem; n = 0; }

	void operator = (const MATRIX& m) {
		if (elem) delete[] elem; 
		n = m.n;
		elem = new double[count]; 
		memcpy(elem, m.elem, count*sizeof(double));
	}

	MATRIX&		operator  =		(const double& x)		{ for (int i = 0; i < count; i++) elem[i] = x;	}	

	MATRIX&		operator +=		(const MATRIX& v)		{ for (int i = 0; i < count; i++) elem[i] += v.elem[i]; }	

	MATRIX&		operator +=		(const double& x)		{ for (int i = 0; i < count; i++) elem[i] += x;	}	

	MATRIX&		operator -=		(const MATRIX& v)		{ for (int i = 0; i < count; i++) elem[i] -= v.elem[i];	}		
	
	MATRIX&		operator -=		(const double& x)		{ for (int i = 0; i < count; i++) elem[i] -= x; }

	MATRIX&		operator *=		(const double& x)		{ for (int i = 0; i < count; i++) elem[i] *= x; }		


	double *	elem;
	int			n;
	int			count;
};

/**
 *	Газодинамические параметры
 */
struct Param
{
	double r;		//!< плотность
	double p;		//!< давление
	double e;		//!< внутренняя энергия
	double E;		//!< полная энергия
	double u;		//!< первая компонента вектора скорости
	double v;		//!< вторая компонента вектора скорости
	double cz;		//!< скорость звука
	double T;		//!< температура
	double ML;		//!< динамическая вязкость
	
	inline double U2() { return  u*u + v*v; }
};

/**
 *	Параметры региона геометрии задачи
 */
struct Region
{
	int		id;
	int		matId;
	int		cellType;		//!< тип ячейки
	Param	par;			//!< параметры региона
};

struct Material
{
	const char*	name;
	int			id;
	
	double		M;		//!< молярная масса
	double		Cp;		//!< теплоемкость при постоянном давлении
	double		ML;		//!< динамическая вязкость
	double		K;		//!< коэффициент теплопроводности

	void URS(Param &par, int flag);
	inline double getGamma() { return Cp/(Cp-gR/M); }

	static double gR;
};

struct Boundary
{
	Boundary(): par(NULL) {}
	~Boundary() { if (par) delete[] par; par = NULL; }
	int			edgeType;
	int			type;
	double *	par;
	int			parCount;

	static const int BOUND_INLET	= 1;
	static const int BOUND_OUTLET	= 2;
	static const int BOUND_WALL		= 3;
};

extern FILE * hLog;


extern void log(char * format, ...);
extern void EXIT(int err);
extern void inverseMatr(double** a_src, double **am, int N);



/**
 *	Решение задачи о распаде произвольного разрыва
 *	(с) ИПМ им. М.В.Келдыша РАН, Тишкин, Никишин, Змитренко
 *
 *	c==========================================================
 *	C    Nikichine
 *	C    module for tube.for /Zmitrenko/
 *	c==========================================================
 */
extern void rim(double& RI, double& EI, double& PI, double& UN, double& UI, double& VI,
         double RB, double PB, double UNB, double UB, double VB,
         double RE, double PE, double UNE, double UE, double VE, Vector n, double GAM);


extern void rim_orig(double& RI, double& EI, double& PI, double& UI, double& VI, double& WI,
         double RB, double PB, double UB, double VB, double WB,
         double RE, double PE, double UE, double VE, double WE, double GAM);

extern void roe_orig(double& RI, double& EI, double& PI, double& UI, double& VI, double& WI,
         double RB, double PB, double UB, double VB, double WB,
         double RE, double PE, double UE, double VE, double WE, double GAM);

#endif