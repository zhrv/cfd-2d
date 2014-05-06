#include "CSR.h"

int CSRMatrix::DELTA = 1024;

CSRMatrix::CSRMatrix(int N)
{
	n  = N;
	_na = 0;
	na = 0;
	a  = 0;
	ja = 0;
	ia = (int*)malloc(sizeof(int)*(n+1));
	memset(ia, 0L, sizeof(int)*(n+1));
}

CSRMatrix::~CSRMatrix()
{
	_na = 0;
	na = 0;
	n = 0;
	free(a);
	free(ia);
	free(ja);
	a  = 0;
	ia = 0;
	ja = 0;
}

void CSRMatrix::zero() 
{
	memset(a, 0, sizeof(double)*na);
	//na = 0;
	//if (a != 0) free(a);
	//if (ja != 0) free(ja);
	//a  = 0;
	//ja = 0;
	//memset(ia, 0L, sizeof(int)*(n+1));
}

void CSRMatrix::set(int i, int j, double aa)
{
	for (int k = ia[i]; k < ia[i+1]; k++)
	{
		if (ja[k] == j)
		{
			a[k] = aa;
			return;
		}
	}
	if (na == _na) {
		double	* ta = (double*)malloc((na + CSRMatrix::DELTA)*sizeof(double));
		int		* tj = (int*)malloc((na + CSRMatrix::DELTA)*sizeof(int));
		memcpy(ta, a, na*sizeof(double));
		memcpy(tj, ja, na*sizeof(int));
		free(a);	free(ja);
		a = ta;		ja = tj;
		_na += CSRMatrix::DELTA;
	}
	int ii = ia[i];
	memcpy(&a[ii + 1], &a[ii], (na - ii)*sizeof(double));
	memcpy(&ja[ii+1], &ja[ii], (na-ii)*sizeof(int));
	a[ii]  = aa;
	ja[ii] = j;
	na++;
	for (int ii = i+1; ii < n+1; ii++)
	{
		ia[ii]++;
	}
}

double CSRMatrix::get(int i, int j)
{
	for (int k = ia[i]; k < ia[i+1]; k++)
	{
		if (ja[k] == j) return a[k];
	}
	return 0.0;
}

void CSRMatrix::add(int i, int j, double aa) {
	set(i, j, get(i,j)+aa);
}

//void CSRMatrix::printToFile(const char *fileName) {
//	FILE * fp = fopen(fileName, "w");
//	for (int i = 0; i < n; i++) {
//		for (int j = 0; j < n; j++) {
//			fprintf(fp, "%16.8e ", get(i, j));
//		}
//		fprintf(fp, "\n");
//	}
//	fclose(fp);
//}
//

void CSRMatrix::printToFile(const char *fileName) {
	FILE * fp = fopen(fileName, "w");
	for (int i = 0; i < n; i++) {
		for (int k = ia[i]; k < ia[i + 1]; k++)
		{
			fprintf(fp, "{{%d: %16.8e}}  ", ja[k], a[k]);    //if (ja[k] == j) return a[k];
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}
