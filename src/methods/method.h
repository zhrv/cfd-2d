#ifndef _METHOD_H_
#define _METHOD_H_

#include "grid.h"

class Method
{
public:
	virtual void init(char * xmlFileName) = 0;
	virtual void run() = 0;
	virtual void done() = 0;
protected:
	Grid grid;
};

#endif