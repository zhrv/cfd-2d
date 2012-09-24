#include <stdlib.h>
#include <stdio.h>
#include <solver.h>
#include "global.h"

int main(int argc, char** argv)
{

	hLog = fopen("task.log", "w"); // открываем файл для записи лога; вместо printf(...) необходимо использовать log(...)
	
	Method * method = Solver::initMethod( "task.xml" ); 
	
	Solver :: runMethod		( method );
	Solver :: destroyMethod	( method );

	fclose(hLog);

	return 0;
}