#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "mpi.h"
#include <math.h>

#ifdef _WIN32
	#include <direct.h>
	#define MK_DIR(name) _mkdir("mesh");
#else
	#include <sys/stat.h>
	#define MK_DIR(name) mkdir("mesh", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
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

	Point() : x(0.0), y(0.0) {}
	Point(double ax, double ay) : x(ax), y(ay) {}
	inline void operator *= (double q) { x *= q; y *= q; }
	inline void operator /= (double q) { x /= q; y /= q; }
	
	inline void operator  = (double q) { x = q;  y = q; }
	inline void operator  = (Point p) { x = p.x;  y = p.y; }
	
	inline void operator += (Point p) { x += p.x; y += p.y; }
	inline void operator -= (Point p) { x -= p.x; y -= p.y; }
	
	inline void operator += (double q) { x += q; y += q; }
	inline void operator -= (double q) { x -= q; y -= q; }
	
	friend inline double operator * (Point p1, Point p2) { return p1.x*p2.x + p1.y*p2.y; }
	friend inline Point operator + (Point p1, Point p2) { return Point(p1.x + p2.x, p1.y + p2.y); }
	friend inline Point operator - (Point p1, Point p2) { return Point(p1.x - p2.x, p1.y - p2.y); }
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

	void		operator *=		(const double& x)		{ for (int i = 0; i < n; i++) elem[i] *= x; }

	void		operator *=		(double** matr)		
	{ 
		double * tmp = new double[n];
		for (int i = 0; i < n; i++) {
			tmp[i] = 0.0;
			for (int j = 0; j < n; j++) {
				tmp[i] += matr[i][j] * elem[j];
			}
		}
		delete[] elem;
		elem = tmp;
	}

	void zero()
	{
		if (elem != NULL && n != 0) {
			memset(elem, 0, n*sizeof(double));
		}
	}

	void abs()
	{
		for (int i = 0; i < n; i++) {
			elem[i] = fabs(elem[i]);
		}
	}

	double norm()
	{
		double s = 0.0;
		for (int i = 0; i < n; i++) {
			s += elem[i] * elem[i];
		}
		return sqrt(s);
	}

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
	inline double magU() { return sqrt(U2()); }
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
	std::string	name;		//!< имя региона
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

class Exception
{
public:
	static const int TYPE_BOUND_NOPAR = 101;
	static const int TYPE_BOUND_UNKNOWN = 102;

	static const int TYPE_MESH_WRONG_NAME = 201;
	static const int TYPE_MESH_UNV_UNKNOWN_ELEMENT = 221;
	static const int TYPE_MESH_UNV_NOT_DEFINED_BND_EDGE = 222;

	static const int FILE_OPENING_ERROR = 301;

	Exception(char* msg, int t) : message(msg), type(t) {}
	char* getMessage() { return message; }
	int getType() { return type; }

private:
	char* message;
	int type;
};


struct Parallel
{
	static void init(int* argc, char*** argv);
	static void done();

	static bool isRoot() { return (procId == 0); }
	static void barrier() { MPI_Barrier(MPI_COMM_WORLD); }

	static void send(int pid, int tag, int n, double* data);
	static void send(int pid, int tag, int n, int* data);
	static void send(int pid, int tag, int n, VECTOR* data);

	static void recv(int pid, int tag, int n, double* data);
	static void recv(int pid, int tag, int n, int* data);
	static void recv(int pid, int tag, int n, VECTOR* data);

	//static void bcast(int tag, int n, double* data);
	//static void bcast(int tag, int n, int* data);
	//static void bcast(int tag, int n, VECTOR* data);

	static int procCount;
	static int procId;
	static double * buf;
};


extern FILE * hLog;

inline double scalar_prod(Vector a, Vector b) { return a.x*b.x + a.y*b.y; }

extern void log(char * format, ...);
extern void EXIT(int err);
extern void inverseMatr(double** a_src, double **am, int N);

extern double** allocMtx(int N);

extern void freeMtx(double** m, int N);



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