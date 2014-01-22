#include <stdlib.h>
#include <stdio.h>
#include <solver.h>
#include "global.h"
/*void ZhWENO(double Um, double Up, double *UU) 
{
double BETA[3];
double ALPHA[3];
double EPS = 0.0001;
//REAL*8 BETA(0:2), ALPHA(0:2), EPS, AS, MINU,MAXU
//COMMON /ACCUR/ EPS
//REAL*8, EXTERNAL :: GW,AMM

BETA[0] = (13.0/12.0)*((UU[2]-2.0*UU[3]+UU[4])*(UU[2]-2.0*UU[3]+UU[4]))+
   0.25*((3.0*UU[2]-4.0*UU[3]+UU[4])*(3.0*UU[2]-4.0*UU[3]+UU[4]));

BETA[1] = (13.0/12.0)*((UU[1]-2.0*UU[2]+UU[3])*(UU[1]-2.0*UU[2]+UU[3]))+
   0.25*((UU[1]-UU[3])*(UU[1]-UU[3]));

BETA[2] = (13.0/12.0)*((UU[0]-2.0*UU[1]+UU[2])*(UU[0]-2.0*UU[1]+UU[2]))+
   0.25*((UU[0]-4.0*UU[1]+3.0*UU[2])*(UU[0]-4.0*UU[1]+3.0*UU[2]));

BETA[2] = (13.0/12.0)*((UU[0]-2.0*UU[1]+UU[2])*(UU[0]-2.0*UU[1]+UU[2]))+
   0.25*((UU[0]-4.0*UU[1]+3.0*UU[2])*(UU[0]-4.0*UU[1]+3.0*UU[2]));
ALPHA[0] = 0.30/((EPS+BETA[0])*(EPS+BETA[0]));
ALPHA[1] = 0.60/((EPS+BETA[1])*(EPS+BETA[1]));
ALPHA[2] = 0.10/((EPS+BETA[2])*(EPS+BETA[2]));

Um = (ALPHA[0]*(2.0*UU[2]+5.0*UU[3]-UU[4])+
  ALPHA[1]*(-UU[1]+5.0*UU[2]+2.0*UU[3])+
  ALPHA[2]*(2.0*UU[0]-7.0*UU[1]+11.0*UU[2]) )/
 ((ALPHA[0]+ALPHA[1]+ALPHA[2])*6.0);
printf("beta[0] = %lf\n",BETA[0]);
printf("beta[1] = %lf\n",BETA[1]);
printf("beta[2] = %lf\n",BETA[2]);
printf("ALPHA[0] = %lf\n",ALPHA[0]);
printf("ALPHA[1] = %lf\n",ALPHA[1]);
printf("ALPHA[2] = %lf\n",ALPHA[2]);

BETA[0] = (13.0/12.0)*((UU[3]-2.0*UU[4]+UU[5])*((UU[3]-2.0*UU[4]+UU[5])))+
     0.25*((3.0*UU[3]-4.0*UU[4]+UU[5])*(3.0*UU[3]-4.0*UU[4]+UU[5]));
BETA[1] = (13.0/12.0)*((UU[2]-2.0*UU[3]+UU[4])*(UU[2]-2.0*UU[3]+UU[4]))+
     0.25*((UU[2]-UU[4])*(UU[2]-UU[4]));
BETA[2] = (13.0/12.0)*((UU[1]-2.0*UU[2]+UU[3])*(UU[1]-2.0*UU[2]+UU[3]))+
     0.25*((UU[1]-4.0*UU[2]+3.0*UU[3])*(UU[1]-4.0*UU[2]+3.0*UU[3]));
ALPHA[0] = 0.1/((EPS+BETA[0])*(EPS+BETA[0]));
ALPHA[1] = 0.6/((EPS+BETA[1])*(EPS+BETA[1]));
ALPHA[2] = 0.3/((EPS+BETA[2])*(EPS+BETA[2]));

Up = (ALPHA[0]*(11.0*UU[3]-7.0*UU[4]+2.0*UU[5])+
 ALPHA[1]*(2.0*UU[2]+5.0*UU[3]-UU[4])+
 ALPHA[2]*(-UU[1]+5.0*UU[2]+2.0*UU[3]) )/
 ((ALPHA[0]+ALPHA[1]+ALPHA[2])*6.0);
printf("beta[0] = %lf\n",BETA[0]);
printf("beta[1] = %lf\n",BETA[1]);
printf("beta[2] = %lf\n",BETA[2]);
printf("ALPHA[0] = %lf\n",ALPHA[0]);
printf("ALPHA[1] = %lf\n",ALPHA[1]);
printf("ALPHA[2] = %lf\n",ALPHA[2]);
printf("um = %lf\n",Um);
printf("up = %lf\n",Up);

}*/

int main(int argc, char** argv)
{
	double a, uum, uup;
	uum = 0.0;
	uup = 0.0;
	double u[6];
	u[0] = 1.0;
	u[1] = 2.0;
	u[2] = 3.0;
	u[3] = 4.0;
	u[4] = 5.0;
	u[5] = 6.0;
	//ZhWENO(uum, uup, u); 
	//printf("uum = %lf\n",uum);
	//printf("uup = %lf\n",uup);
	hLog = fopen("task.log", "w"); // открываем файл для записи лога; вместо printf(...) необходимо использовать log(...)
	
	Method * method = Solver::initMethod( "task.xml" ); 
	
	Solver :: runMethod		( method );
	Solver :: destroyMethod	( method );

	fclose(hLog);

	return 0;
}