#include "CSR.h"

CSRMatrix::CSRMatrix(int N)
{
	n  = N;
	na = 0;
	a  = 0;
	ja = 0;
	ia = (int*)malloc(sizeof(int)*(n+1));
	memset(ia, 0L, sizeof(int)*(n+1));
}

CSRMatrix::~CSRMatrix()
{
	na = 0;
	n  = 0;
	free(a);
	free(ia);
	free(ja);
	a  = 0;
	ia = 0;
	ja = 0;
}

void CSRMatrix::zero() 
{
	na = 0;
	if (a != 0) free(a);
	if (ja != 0) free(ja);
	a  = 0;
	ja = 0;
	memset(ia, 0L, sizeof(int)*(n+1));
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
	int ii = ia[i];
	double	* ta = (double*) malloc((na+1)*sizeof(double));
	int		* tj = (int*)    malloc((na+1)*sizeof(int));
	memcpy(ta, a,  na*sizeof(double));
	memcpy(tj, ja, na*sizeof(int));
	free(a);	free(ja);
	a = ta;		ja = tj;
	memcpy(&a[ii+1],  &a[ii],  (na-ii)*sizeof(double));
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

void CSRMatrix::printToFile(const std::string& fileName) {
	FILE * fp = fopen(fileName.c_str(), "w");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fprintf(fp, "%16.8e ", get(i,j));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}
