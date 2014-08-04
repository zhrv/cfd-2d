#include <stdlib.h>
#include <stdio.h>
#include <solver.h>
#include "global.h"
#include "mpi.h"

int main(int argc, char** argv)
{
	int myid, num_procs;

	/* Initialize MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	hLog = fopen("task.log", "w"); // открываем файл для записи лога; вместо printf(...) необходимо использовать log(...)
	
	Method * method = Solver::initMethod( "task.xml" ); 
	
	Solver :: runMethod		( method );
	Solver :: destroyMethod	( method );

	fclose(hLog);

	return 0;
}