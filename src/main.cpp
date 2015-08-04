#include <stdlib.h>
#include <stdio.h>
#include <solver.h>
#include "global.h"
#include "mpi.h"
#include <float.h>

#include "SolverZeidel.h"

int main(int argc, char** argv)
{
	//CSRMatrix mtr(4);
	//mtr.set(0, 0, 1);
	//mtr.set(0, 1, 1);
	//mtr.add(0, 0, 1);
	//mtr.set(1, 1, 1);
	//mtr.add(3, 3, 4);
	//mtr.add(2, 3, 9);
	//mtr.set(3, 1, 10);
	//mtr.assemble();

#ifdef _DEBUG
	_controlfp(~(_MCW_EM & (~_EM_INEXACT) & (~_EM_UNDERFLOW)), _MCW_EM);
#endif

	Parallel::init(&argc, &argv);

	hLog = fopen("task.log", "w"); // открываем файл для записи лога; вместо printf(...) необходимо использовать log(...)
	
	Method * method = Solver::initMethod( "task.xml" ); 
	
	Solver :: runMethod		( method );
	Solver :: destroyMethod	( method );

	fclose(hLog);

	Parallel::done();

	return 0;
}