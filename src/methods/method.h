#ifndef _METHOD_H_
#define _METHOD_H_

#include "grid.h"
#include "global.h"
class Method
{
public:
	virtual void init(char * xmlFileName) = 0;
	virtual void run() = 0;
	virtual void done() = 0;

	void exchange(double* field)
	{
		for (int p = 0; p < Parallel::procCount; p++) {
			if (p < Parallel::procId) {

			}
			else if (p > Parallel::procId) {
			}
		}
	}

	void exchange(int* field)
	{

	}

	void exchange(VECTOR* field)
	{

	}

protected:
	Grid grid;

	double * dBuf;
	int    * iBuf;
	VECTOR * vBuf;
};

#endif