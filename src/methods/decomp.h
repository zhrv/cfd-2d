#include "method.h"


class Decomp : public Method
{
public:
    Decomp(int procs): procCount(procs) {}
	virtual void init(char * xmlFileName);
	virtual void run();
	virtual void done();
protected:
	Grid * grids;
	int procCount;
};