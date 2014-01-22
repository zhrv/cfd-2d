#include "grid.h"

Cell::~Cell() 
{
	delete[] nodesInd;
	delete[] edgesInd;
}

Edge::~Edge() 
{
	delete[] c;
}

Grid::~Grid() 
{
	delete[] nodes;
	delete[] cells;
	delete[] edges;
}

void Grid::replaceEdges(int if1, int if2) 
{
	Edge tmp = edges[if1];
	for (int i = 0; i < 0; i++) 
	{

	}
}

void Grid::reorderEdges() 
{

}

int Grid::findEdge(int n1, int n2)
{
	for (int iEdge = 0; iEdge < eCount; iEdge++)
	{
		if ((edges[iEdge].n1 == n1 && edges[iEdge].n2 == n2) || (edges[iEdge].n1 == n2 && edges[iEdge].n2 == n1))
		{
			return iEdge;
		}
	}
	return -1;
}
int Grid::location_node (int possible_tri,double x0,double y0)
{
	double x1,y1,x2,y2,x3,y3; //координаты вершины фиксированного треугольника
	int num_tri;//искомый треугольник
	double x_center = 0.0;
	double y_center = 0.0;
	//printf("possible_tri= %d x0= %lf y0= %lf\n",possible_tri,x0,y0);
	int sum = 0;
			 num_tri = 0;
			 //printf("x_center = %d\n",cCount);
		 for (int j = 0; j < cCount; j++){
			x1 = nodes[cells[j].nodesInd[0]].x;
			x2 = nodes[cells[j].nodesInd[1]].x;
			x3 = nodes[cells[j].nodesInd[2]].x;
			y1 = nodes[cells[j].nodesInd[0]].y;
			y2 = nodes[cells[j].nodesInd[1]].y;
			y3 = nodes[cells[j].nodesInd[2]].y;
			//printf("x_center = %lf\n",x1);
			x_center = cells[j].c.x;
		    y_center = cells[j].c.y;
			
			
			//определяем лежит ли точка пересеч серед перп внутри соответсвующего треугольника
			if (x1 != x2){ 	
				if (y_center <= (x_center*(y2-y1)/(x2-x1)+(y1-x1*(y2-y1)/(x2-x1)))) 
				{
					if (y0 <= (x0*(y2-y1)/(x2-x1)+(y1-x1*(y2-y1)/(x2-x1)))){
						//if (j==9) printf("vnutri11\n");
						//printf("%lf %lf %lf\n",y_center,y0,x0*(y2-y1)/(x2-x1)+(y1-x1/(y2-y1)/(x2-x1)));
						sum +=1;}
				}
				else
				{
					if (y0 >= (x0*(y2-y1)/(x2-x1)+(y1-x1*(y2-y1)/(x2-x1)))){
						//if (j==9) printf("vnutri12\n");
						//printf("%lf %lf %lf\n",y_center,y0,x0*(y2-y1)/(x2-x1)+(y1-x1/(y2-y1)/(x2-x1)));
						sum +=1;}
				}
			}
			if (x1 != x3){
				if (y_center <= (x_center*(y3-y1)/(x3-x1)+(y1-x1*(y3-y1)/(x3-x1)))) 
				{
					if (y0 <= (x0*(y3-y1)/(x3-x1)+(y1-x1*(y3-y1)/(x3-x1)))){
						//if (j==9) printf("vnutri21\n");
						sum +=1;}
				}
				else
				{
					if (y0 >= (x0*(y3-y1)/(x3-x1)+(y1-x1*(y3-y1)/(x3-x1)))){
						//if (j==9) printf("vnutri22\n");
						//printf("%lf %lf %lf\n",y_center,y0,x0*(y3-y1)/(x3-x1)+(y1-x1*(y3-y1)/(x3-x1)));
						sum +=1;}
				}
			}
			if (x2 != x3){
				if (y_center <= (x_center*(y2-y3)/(x2-x3)+(y3-x3*(y2-y3)/(x2-x3)))) 
				{
					//if (j==102) printf("vnutri331\n");
					if (y0 <= (x0*(y2-y3)/(x2-x3)+(y3-x3*(y2-y3)/(x2-x3)))){
						//if (j==102) printf("vnutri31\n");
						//printf("%lf %lf %lf\n",y_center,y0,x0*(y2-y3)/(x2-x3)+(y3-x3*(y2-y3)/(x2-x3)));
						sum +=1;}
				}
				else
				{
					if (y0 >= (x0*(y2-y3)/(x2-x3)+(y3-x3*(y2-y3)/(x2-x3)))){
						//if (j==102) printf("vnutri32\n");
						sum +=1;}
				}
			}
			if (x1 == x2){
				//if (j==21) printf("sfdd1 %d\n",sum);
				if (((x_center <= x1) && (x0 <= x1))||
					((x_center >= x1)&&(x0 >= x1))||
					((x_center == x1) && (x0 == x1)))
				{
					sum +=1;
					//printf("%d\n",sum);
				}
			}
			//if (j==9) printf("j = 9 sum = %d\n",sum);
			if (x2 == x3){
				//if (j==21) printf("sfdd2 %d\n",sum);
				if (((x_center <= x2) && (x0 <= x2))||
					((x_center >= x2)&&(x0 >= x2))||
					((x_center == x2) && (x0 == x2)))
				{
					sum +=1;
					//if (j==9) printf("vnutri423 %lf %lf\n",y_center,x2);
				}
			}
			if (x1 == x3){
				//if (j==21) {printf("sfdd3 %d\n",sum);
				//printf("x0 = %.10lf\n",x0);
				//printf("x1 = %.10lf\n",x1);
				//printf("x3 = %.10lf\n",x3);
				//printf("x_center = %.10lf\n",x_center);}
				if (((x_center <= x1) && (x0 <= x1))||
					((x_center >= x1)&&(x0 >= x1))||
					((x_center == x1)&&(x0 == x1)))
				{
					sum +=1;
					//if (j==21) printf("sfdd3 %d\n",sum);
					//printf("%d\n",sum);
				}
			}
			if (sum > 3) printf("%d\n",sum);
			if (sum == 3)  {
				//if ((pbound[i] != 0) printf("!!!!!!!!!!\n");
				if (num_tri != possible_tri) 
				{
					num_tri = j;
					//printf("ro = %d\n",num_tri);
					//printf("%d %d\n",sum_tmp,sum);
					//printf("%d (%lf;%lf) (%lf;%lf) (%lf;%lf) \n",j,x1,y1,x2,y2,x3,y3);
				}
				//printf("%d (%lf;%lf) (%lf;%lf) (%lf;%lf) (%lf;%lf) (%lf;%lf) \n",i,x1,y1,x2,y2,x3,y3,x_center,y_center,x0,y0);
				//printf("%d\n",pbound[i]);
			}
			sum = 0;
			//printf("%d %lf %lf\n",i,new_nodes_coord[1][i],new_nodes_coord[2][i]);
			}
	return num_tri;
}
int Grid::inArea (double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
	double k1, k2; //угловые коэффициенты прямых
	double b1 = 0.0, b2 = 0.0; //свободные члены прямых
	double x0, y0; //точка пересечения прямых
	if (x1 != x2)
	{
		k1 = (y2-y1)/(x2-x1);
		b1 = y1-k1*x1;
	}
	else 
	{
		k1 = 0.0;
	}

	if (x3 != x4)
	{
		k2 = (y4-y3)/(x4-x3);
		b2 = y3-k2*x3;
	}
	else 
	{
		k2 = 0.0;
	}

	if ((k1 != k2)&&(k1 != 0.0)&&(k2 != 0.0)) 
	{
		x0 = (b1-b2)/(k2-k1);
		y0 = k1*x0+b1;
		if ((((x0-x1)*(x0-x2)+(y0-y1)*(y0-y2))>0) || (((x0-x3)*(x0-x4)+(y0-y3)*(y0-y4))>0)) 
		{
			return 0;
		}
		else 
		{
			return 1;
		}
	}

	if ((k1 == 0.0)&&(k2 != 0.0))
	{
		if (y1 == y2)
		{
			y0 = y1;
			x0 = (y0-b2)/k2;
		}
		else
		{
			x0 = x1;
			y0 = k2*x0+b2;
		}
		if ((((x0-x1)*(x0-x2)+(y0-y1)*(y0-y2))>0) || (((x0-x3)*(x0-x4)+(y0-y3)*(y0-y4))>0)) 
		{
			return 0;
		}
		else 
		{
			return 1;
		}
	}

	if ((k2 == 0.0)&&(k1 != 0.0))
	{
		if (y3 == y4)
		{
			y0 = y3;
			x0 = (y0-b1)/k1;
		}
		else
		{
			x0 = x3;
			y0 = k1*x0+b1;
		}
		if ((((x0-x1)*(x0-x2)+(y0-y1)*(y0-y2))>0) || (((x0-x3)*(x0-x4)+(y0-y3)*(y0-y4))>0)) 
		{
			return 0;
		}
		else 
		{
			return 1;
		}
	}
	
	if ((k1 == k2)&&(k1 != 0.0))
	{
		return 0;
	}

	if ((k1 == 0.0)&&(k2==0.0))
	{
		if ((x1==x2) && (y3==y4))
		{
			x0 = x1;
			y0 = y3;
			if ((((x0-x1)*(x0-x2)+(y0-y1)*(y0-y2))>0) || (((x0-x3)*(x0-x4)+(y0-y3)*(y0-y4))>0)) 
			{
				return 0;
			}
			else 
			{
				return 1;
			}
		}
		if ((x3==x4) && (y1==y2))
		{
			x0 = x3;
			y0 = y1;
			if ((((x0-x1)*(x0-x2)+(y0-y1)*(y0-y2))>0) || (((x0-x3)*(x0-x4)+(y0-y3)*(y0-y4))>0)) 
			{
				return 0;
			}
			else 
			{
				return 1;
			}
		}
	}
}
void Grid::WENOforEdge (int node1, int node2, int cellL, int cellR, double* uwx, double* uwy)
{
    double x[2];
    double y[2];
	x[0] = nodes[node1].x;
    x[1] = nodes[node2].x;
    y[0] = nodes[node1].y;
    y[1] = nodes[node2].y;
	/*if (node1 == 1)
	{
		printf("x[0] = %lf y[0] = %lf\n",x[0],y[0]);
		printf("x[1] = %lf y[1] = %lf\n",x[1],y[1]);
	}*/
    double xtmp1, xtmp2, ytmp1, ytmp2;
    
    //cells[i].sideType = 0;
    xtmp1 =                (((y[1]+y[0])/2.0-cells[cellL].c.y)*(x[0]-x[1])*(y[0]-y[1])+cells[cellL].c.x*(y[0]-y[1])*(y[0]-y[1])+(x[1]+x[0])/2.0*(x[0]-x[1])*(x[0]-x[1]))/
                                     ((y[0]-y[1])*(y[0]-y[1])+(x[0]-x[1])*(x[0]-x[1]));
    if (x[1] != x[0]) {
                ytmp1 = (y[0]-y[1])/(x[0]-x[1])*xtmp1+cells[cellL].c.y-(y[0]-y[1])/(x[0]-x[1])*cells[cellL].c.x;}
            else
            {
                ytmp1 = cells[cellL].c.y;
            }

    xtmp2 =                (((y[1]+y[0])/2.0-cells[cellR].c.y)*(x[0]-x[1])*(y[0]-y[1])+cells[cellR].c.x*(y[0]-y[1])*(y[0]-y[1])+(x[1]+x[0])/2.0*(x[0]-x[1])*(x[0]-x[1]))/
                                     ((y[0]-y[1])*(y[0]-y[1])+(x[0]-x[1])*(x[0]-x[1]));
    if (x[1] != x[0]) {
                ytmp2 = (y[0]-y[1])/(x[0]-x[1])*xtmp2+cells[cellR].c.y-(y[0]-y[1])/(x[0]-x[1])*cells[cellR].c.x;}
            else
            {
                ytmp2 = cells[cellR].c.y;
            }
    uwx[2] = (xtmp1-xtmp2)/2.0+(x[0]+x[1])/2.0;
    uwy[2] = (ytmp1-ytmp2)/2.0+(y[0]+y[1])/2.0;
    //printf("xtmp1=%lf xtmp2=%lf x0=%lf x=%lf\n",xtmp1,xtmp2,(x[0]+x[1])/2.0,uwx[2]);
    uwx[1] = 3.0*uwx[2]-(x[1]+x[0]);
    uwy[1] = 3.0*uwy[2]-(y[1]+y[0]);
    uwx[0] = 3.0*uwx[1]-4*uwx[2]+x[1]+x[0];
    uwy[0] = 3.0*uwy[1]-4*uwy[2]+y[1]+y[0];
    uwx[3] = x[1]+x[0]-uwx[2];
    uwy[3] = y[1]+y[0]-uwy[2];
    uwx[4] = 3.0*uwx[3]-x[1]-x[0];
    uwy[4] = 3.0*uwy[3]-y[1]-y[0];
    uwx[5] = 3.0*uwx[4]-4.0*uwx[3]+x[1]+x[0];
    uwy[5] = 3.0*uwy[4]-4.0*uwy[3]+y[1]+y[0];
    //printf("x[0]=%lf x[1]=%lf x[2]=%lf x[3]=%lf x[4]=%lf x[5]=%lf\n",uwx[0],uwx[1],uwx[2],uwx[3],uwx[4],uwx[5]);
    //printf("y[0]=%lf y[1]=%lf y[2]=%lf y[3]=%lf y[4]=%lf y[5]=%lf\n",uwy[0],uwy[1],uwy[2],uwy[3],uwy[4],uwy[5]);

                
      /*cells[i].side[j][1].x = 3.0*cells[i].side[j][2].x-(x[2]+x[1]);  
        cells[i].side[j][1].y = 3.0*cells[i].side[j][2].y-(y[2]+y[1]);
        cells[i].side[j][0].x = 3.0*cells[i].side[j][1].x-4*cells[i].side[j][2].x+x[2]+x[1];
        cells[i].side[j][0].y = 3.0*cells[i].side[j][1].y-4*cells[i].side[j][2].y+y[2]+y[1];
        cells[i].side[j][3].x = x[2]+x[1]-cells[i].side[j][2].x;
        cells[i].side[j][3].y = y[2]+y[1]-cells[i].side[j][2].y;
        cells[i].side[j][4].x = 3.0*cells[i].side[j][3].x-x[2]-x[1];
        cells[i].side[j][4].y = 3.0*cells[i].side[j][3].y-y[2]-y[1];
        cells[i].side[j][5].x = 3.0*cells[i].side[j][4].x-4.0*cells[i].side[j][3].x+x[2]+x[1];
        cells[i].side[j][5].y = 3.0*cells[i].side[j][4].y-4.0*cells[i].side[j][3].y+y[2]+y[1];*/

}

void Grid::initFromFiles(char* fName) 
{
	char str[50];
	FILE *fp;
	int tmp; 
 
	//вспомогательные переменные
	/*double vtmpx, vtmpy;
	double x[3];
	double y[3];*/
	// читаем данные об УЗЛАХ
	sprintf(str, "%s.node", fName);
	fp = fopen(str, "r");
	fscanf(fp, "%d %d %d %d", &nCount, &tmp, &tmp, &tmp);
	nodes = new Point[nCount];
	for (int i = 0; i < nCount; i++) 
	{
		fscanf(fp, "%d %lf %lf %d", &tmp, &(nodes[i].x), &(nodes[i].y), &tmp);
	}
	fclose(fp);

	// читаем данные о ЯЧЕЙКАХ
	sprintf(str, "%s.ele", fName);
	fp = fopen(str, "r");
	fscanf(fp, "%d %d %d", &cCount, &tmp, &tmp);
	cells = new Cell[cCount];
	double RadiusMiddle = 0.0;
	double abc = 0.0;
	for (int i = 0; i < cCount; i++) 
	{
		cells[i].nCount = 3;
		cells[i].nodesInd = new int[cells[i].nCount];
		fscanf(fp, "%d %d %d %d %d", &tmp, &(cells[i].nodesInd[0]), &(cells[i].nodesInd[1]), &(cells[i].nodesInd[2]), &(cells[i].type));
		cells[i].nodesInd[0]--;
		cells[i].nodesInd[1]--;
		cells[i].nodesInd[2]--;
		cells[i].c.x = (nodes[cells[i].nodesInd[0]].x+nodes[cells[i].nodesInd[1]].x+nodes[cells[i].nodesInd[2]].x)/3.0;
		cells[i].c.y = (nodes[cells[i].nodesInd[0]].y+nodes[cells[i].nodesInd[1]].y+nodes[cells[i].nodesInd[2]].y)/3.0;
		cells[i].HX = _max_( fabs(nodes[cells[i].nodesInd[0]].x-nodes[cells[i].nodesInd[1]].x), 
			                 fabs(nodes[cells[i].nodesInd[1]].x-nodes[cells[i].nodesInd[2]].x),
							 fabs(nodes[cells[i].nodesInd[0]].x-nodes[cells[i].nodesInd[2]].x) );
		cells[i].HY = _max_( fabs(nodes[cells[i].nodesInd[0]].y-nodes[cells[i].nodesInd[1]].y), 
			                 fabs(nodes[cells[i].nodesInd[1]].y-nodes[cells[i].nodesInd[2]].y),
							 fabs(nodes[cells[i].nodesInd[0]].y-nodes[cells[i].nodesInd[2]].y) );
		cells[i].eCount = 3;
		cells[i].edgesInd = new int[cells[i].eCount];
		double x[3], y[3];
		x[0] = nodes[cells[i].nodesInd[0]].x;
		x[1] = nodes[cells[i].nodesInd[1]].x;
		x[2] = nodes[cells[i].nodesInd[2]].x;
		y[0] = nodes[cells[i].nodesInd[0]].y;
		y[1] = nodes[cells[i].nodesInd[1]].y;
		y[2] = nodes[cells[i].nodesInd[2]].y;
		double a, b, c, p;
		a = sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
		b = sqrt((x[2]-x[1])*(x[2]-x[1])+(y[2]-y[1])*(y[2]-y[1]));
		c = sqrt((x[0]-x[2])*(x[0]-x[2])+(y[0]-y[2])*(y[0]-y[2]));
		abc += a*b*c;
		p = 0.5*(a+b+c);
		cells[i].R = a*b*c/(4*sqrt(p*(p-a)*(p-b)*(p-c)));
		//printf("R[%d] = %lf\n",i,cells[i].R);
		RadiusMiddle += cells[i].R;
	}
	fclose(fp);
	RadiusMiddle /= cCount;
	abc /= cCount;
	printf("RadiusMiddle = %lf\n",RadiusMiddle);
	printf("abc = %lf\n",abc);

	int num_tri = 0;
	/*for (int i = 0; i < cCount; i++) 
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 6; k++)
	{
	  if (cells[i].sideType == 0) {
		  num_tri = location_node (1,cells[i].side[j][k].x,cells[i].side[j][k].y);
		  printf("%d ",i);
	  }
	}*/
	// формируем данные о РЕБРАХ
	sprintf(str, "%s.neigh", fName);
	fp = fopen(str, "r");
	fscanf(fp, "%d %d", &tmp, &tmp);
	int** neigh;
	neigh = new int*[cCount]; 
	for (int i = 0; i < cCount; i++) 
	{
		neigh[i] = new int[3];
		fscanf(fp, "%d %d %d %d", &tmp, &(neigh[i][0]), &(neigh[i][1]), &(neigh[i][2]));
		neigh[i][0]--;
		neigh[i][1]--;
		neigh[i][2]--;
		cells[i].neigh[0] = neigh[i][0];
		cells[i].neigh[1] = neigh[i][1];
		cells[i].neigh[2] = neigh[i][2];
	}
	fclose(fp);
	eCount = 0;
	for (int i = 0; i < cCount; i++) 
	{
		for (int j = 0; j < 3; j++) 
		{
			int p = neigh[i][j];
			if (p > -1) 
			{
				for (int k = 0; k < 3; k++) 
				{ // убираем у соседа номер этой ячейки, чтобы грань не повторялась
					if (neigh[p][k] == i) neigh[p][k] = -1;
				}
				eCount++;
			}
			if (p == -2) eCount++;
		}
	}
	edges = new Edge[eCount];

	int iEdge = 0;
	int * cfi = new int[cCount];
	for (int i = 0; i < cCount; i++) 
	{
		cfi[i] = 0;
	}
	// ::memset(cfi, 0, cCount*sizeof(int));
	double trix[6];
	double triy[6];
	for (int i = 0; i < cCount; i++) 
	{
		for (int j = 0; j < 3; j++) 
		{
			int p = neigh[i][j];
			if (p != -1) 
			{
				edges[iEdge].n1     = cells[i].nodesInd[(j+1)%3];
				edges[iEdge].n2     = cells[i].nodesInd[(j+2)%3];
				
				edges[iEdge].cCount = 3;
				edges[iEdge].c      = new Point[edges[iEdge].cCount];
				double _sqrt3 = 1.0/sqrt(3.0);
				// центр ребра
				edges[iEdge].c[0].x = (nodes[edges[iEdge].n1].x+nodes[edges[iEdge].n2].x)/2.0;
				edges[iEdge].c[0].y = (nodes[edges[iEdge].n1].y+nodes[edges[iEdge].n2].y)/2.0;
				// первая точка Гаусса
				edges[iEdge].c[1].x = (nodes[edges[iEdge].n1].x+nodes[edges[iEdge].n2].x)/2.0-_sqrt3*(nodes[edges[iEdge].n2].x-nodes[edges[iEdge].n1].x)/2.0;
				edges[iEdge].c[1].y = (nodes[edges[iEdge].n1].y+nodes[edges[iEdge].n2].y)/2.0-_sqrt3*(nodes[edges[iEdge].n2].y-nodes[edges[iEdge].n1].y)/2.0;
				// вторая точка Гаусса
				edges[iEdge].c[2].x = (nodes[edges[iEdge].n1].x+nodes[edges[iEdge].n2].x)/2.0+_sqrt3*(nodes[edges[iEdge].n2].x-nodes[edges[iEdge].n1].x)/2.0;
				edges[iEdge].c[2].y = (nodes[edges[iEdge].n1].y+nodes[edges[iEdge].n2].y)/2.0+_sqrt3*(nodes[edges[iEdge].n2].y-nodes[edges[iEdge].n1].y)/2.0;

				edges[iEdge].n.x    = nodes[edges[iEdge].n2].y-nodes[edges[iEdge].n1].y;
				edges[iEdge].n.y    = nodes[edges[iEdge].n1].x-nodes[edges[iEdge].n2].x;
				edges[iEdge].l      = sqrt(edges[iEdge].n.x*edges[iEdge].n.x+edges[iEdge].n.y*edges[iEdge].n.y);
				edges[iEdge].n.x    /= edges[iEdge].l;
				edges[iEdge].n.y    /= edges[iEdge].l;
				edges[iEdge].c1     = i;
				cells[i].edgesInd[cfi[i]] = iEdge;
				cfi[i]++;
				edges[iEdge].cnl1 = fabs(edges[iEdge].n.x*(edges[iEdge].c[0].x-cells[edges[iEdge].c1].c.x)+edges[iEdge].n.y*(edges[iEdge].c[0].y-cells[edges[iEdge].c1].c.y) );

				if (p > -1) 
				{

					edges[iEdge].c2 = p;
					cells[p].edgesInd[cfi[p]] = iEdge;
					cfi[p]++;
					edges[iEdge].cnl2 = fabs(edges[iEdge].n.x*(cells[edges[iEdge].c2].c.x-edges[iEdge].c[0].x)+edges[iEdge].n.y*(cells[edges[iEdge].c2].c.y-edges[iEdge].c[0].y) );
					edges[iEdge].type = Edge::TYPE_INNER;
				}
				if (p == -2) 
				{
					edges[iEdge].c2 = -1;
					edges[iEdge].cnl2 = 0;
					edges[iEdge].type = -1;
				}
				if (edges[iEdge].type != -1) {
					WENOforEdge (edges[iEdge].n1, edges[iEdge].n2, edges[iEdge].c1, edges[iEdge].c2, trix, triy);
				}
				edges[iEdge].otr[0].x = trix[0];
				edges[iEdge].otr[1].x = trix[1];
				edges[iEdge].otr[2].x = trix[2];
				edges[iEdge].otr[3].x = trix[3];
				edges[iEdge].otr[4].x = trix[4];
				edges[iEdge].otr[5].x = trix[5];
				
				edges[iEdge].otr[0].y = triy[0];
				edges[iEdge].otr[1].y = triy[1];
				edges[iEdge].otr[2].y = triy[2];
				edges[iEdge].otr[3].y = triy[3];
				edges[iEdge].otr[4].y = triy[4];
				edges[iEdge].otr[5].y = triy[5];
				iEdge++;
			}
		}
	}

	int * tri_tmp = new int[cCount];
	for (int i = 0; i < cCount; i++)
	{
		tri_tmp[i] = 0;
	}
	for (int i = 0; i < iEdge; i++) 
	{
		for (int k = 0; k < iEdge; k++) 
		{
				if (edges[k].type == -1) 
				{
					edges[k].wenoType = 1;
					if (edges[i].wenoType != 1)
						{
							edges[i].wenoType = inArea(nodes[edges[k].n1].x, nodes[edges[k].n1].y,nodes[edges[k].n2].x,nodes[edges[k].n2].y,edges[i].otr[0].x,edges[i].otr[0].y,edges[i].otr[5].x,edges[i].otr[5].y);
					        if (edges[i].wenoType == 1)
							{
								tri_tmp[edges[i].c1]--;
			                    tri_tmp[edges[i].c2]--;
							}
					    }
				}
		}
		/*if (i==23)
		{
			printf("wenoType = %d\n",edges[i].wenoType);
		}*/
		if ((edges[i].c1 > 0) && (edges[i].c2 > 0))
		{
			
			tri_tmp[edges[i].c1]++;
			tri_tmp[edges[i].c2]++;
		}
	}

	/*for (int i = 0; i < cCount; i++)
	{
		if (tri_tmp[i] == 1) printf("tri_tmp[%d] = %d\n",i,tri_tmp[i]);
	}*/

	//for (int i = 0; i < cCount; i++)
		for (int k = 0; k < iEdge; k++) 
			{
				if (edges[k].type == -1) 
				{
					edges[k].wenoType = 1;
				}
				/*if (edges[k].wenoType != 1)
				{		
					edges[k].wenoType = inArea(nodes[edges[k].n1].x, nodes[edges[k].n1].y,nodes[edges[k].n2].x,nodes[edges[k].n2].y,edges[i].otr[0].x,edges[i].otr[0].y,edges[i].otr[5].x,edges[i].otr[5].y);
				}*/
				/*if (k==23)
						{
							printf("wenoType = %d\n",edges[k].wenoType);
						}*/
				if ((edges[k].wenoType != 1) && ((tri_tmp[edges[k].c1] != 3) && (tri_tmp[edges[k].c2] != 3)))
				{edges[k].wenoType = 1;}
				/*if (k==23)
						{
							printf("wenoType = %d\n",edges[k].wenoType);
						}*/
			}
	FILE *f11=fopen("raspad.gp", "w");
	for (int i = 0; i < cCount; i++) 
	{  
		fprintf(f11,"%f\t%f\n",nodes[cells[i].nodesInd[0]].x,nodes[cells[i].nodesInd[0]].y);
		fprintf(f11,"%f\t%f\n",nodes[cells[i].nodesInd[1]].x,nodes[cells[i].nodesInd[1]].y);
		fprintf(f11,"%f\t%f\n",nodes[cells[i].nodesInd[2]].x,nodes[cells[i].nodesInd[2]].y);
		fprintf(f11,"%f\t%f\n",nodes[cells[i].nodesInd[0]].x,nodes[cells[i].nodesInd[0]].y);
		fprintf(f11,"\n");
	}
	for (int i = 0; i < iEdge; i++) 
	{ 	
		if (edges[i].wenoType != 1)
		{
			for (int j = 0; j < 6; j++)
			{
			  edges[i].seq_tri[j] = location_node (1,edges[i].otr[j].x,edges[i].otr[j].y);
			  //printf("%d ",edges[i].seq_tri[j]);
			}
			//printf("tri = %d\n",edges[i].n1);
			//printf("tri = %d\n",edges[i].n2);
			//printf("\n");
		
		fprintf(f11,"%f\t%f\n",edges[i].otr[0].x,edges[i].otr[0].y);
		fprintf(f11,"%f\t%f\n",edges[i].otr[1].x,edges[i].otr[1].y);
		fprintf(f11,"%f\t%f\n",edges[i].otr[2].x,edges[i].otr[2].y);
		fprintf(f11,"%f\t%f\n",edges[i].otr[3].x,edges[i].otr[3].y);
		fprintf(f11,"%f\t%f\n",edges[i].otr[4].x,edges[i].otr[4].y);
		fprintf(f11,"%f\t%f\n",edges[i].otr[5].x,edges[i].otr[5].y);
		fprintf(f11,"%f\t%f\n",edges[i].otr[0].x,edges[i].otr[0].y);
		fprintf(f11,"\n");
		}
	}

	/*for (int i = 0; i < cCount; i++) 
		for (int k = 0; k < iEdge; k++) 
			for (int j = 0; j < 3; j++) 
			{
				if ((cells[i].sideType == 0) && (edges[k].type == -1))
						{
							cells[i].sideType = inArea(nodes[edges[k].n1].x, nodes[edges[k].n1].y,nodes[edges[k].n2].x,nodes[edges[k].n2].y,cells[i].side[j][0].x,cells[i].side[j][0].y,cells[i].side[j][5].x,cells[i].side[j][5].y);
					    }
			}*/
	/*FILE *f11=fopen("raspad.gp", "w");
	for (int i = 0; i < cCount; i++) 
	{  
		fprintf(f11,"%f\t%f\n",nodes[cells[i].nodesInd[0]].x,nodes[cells[i].nodesInd[0]].y);
		fprintf(f11,"%f\t%f\n",nodes[cells[i].nodesInd[1]].x,nodes[cells[i].nodesInd[1]].y);
		fprintf(f11,"%f\t%f\n",nodes[cells[i].nodesInd[2]].x,nodes[cells[i].nodesInd[2]].y);
		fprintf(f11,"%f\t%f\n",nodes[cells[i].nodesInd[0]].x,nodes[cells[i].nodesInd[0]].y);
		fprintf(f11,"\n");
		if (cells[i].sideType == 0) 
			for (int j = 0; j < 3; j++)
			 {
		        if (cells[i].sideType == 0) {
				fprintf(f11,"%f\t%f\n",cells[i].side[j][0].x,cells[i].side[j][0].y);
				fprintf(f11,"%f\t%f\n",cells[i].side[j][1].x,cells[i].side[j][1].y);
				fprintf(f11,"%f\t%f\n",cells[i].side[j][2].x,cells[i].side[j][2].y);
				fprintf(f11,"%f\t%f\n",cells[i].side[j][3].x,cells[i].side[j][3].y);
				fprintf(f11,"%f\t%f\n",cells[i].side[j][4].x,cells[i].side[j][4].y);
				fprintf(f11,"%f\t%f\n",cells[i].side[j][5].x,cells[i].side[j][5].y);
				fprintf(f11,"%f\t%f\n",cells[i].side[j][0].x,cells[i].side[j][0].y);
				fprintf(f11,"\n");}

			  }
	}
	fprintf(f11,"%f\t%f\n",0.0000,-0.059641);
				fprintf(f11,"%f\t%f\n",1.0000,-0.059641);
				fprintf(f11,"%f\t%f\n",1.0000,0.059641);
				fprintf(f11,"%f\t%f\n",0.0000,0.059641);
				fprintf(f11,"%f\t%f\n",0.0000,-0.059641);
				fprintf(f11,"\n");
	fclose(f11);*/
	
	// чтение данных о граничных гранях
	sprintf(str, "%s.poly", fName);
	fp = fopen(str, "r");
	int bndCount;
	fscanf(fp, "%d %d %d %d", &tmp, &tmp, &tmp, &tmp);
	fscanf(fp, "%d %d", &bndCount, &tmp);
	for (int i = 0; i < bndCount; i++) 
	{
		int n, n1, n2, type;
		fscanf(fp, "%d %d %d %d", &n, &n1, &n2, &type);
		n1--;
		n2--;
		int iEdge = findEdge(n1, n2);
		if (iEdge >= 0) edges[iEdge].type = type;
	}
	fclose(fp);

	for (int i = 0; i < cCount; i++) 
	{
		double a = edges[cells[i].edgesInd[0]].l;
		double b = edges[cells[i].edgesInd[1]].l;
		double c = edges[cells[i].edgesInd[2]].l;
		double p = (a+b+c)/2.0;
		cells[i].S = sqrt(p*(p-a)*(p-b)*(p-c));
	}

	for (int i = 0; i < cCount; i++) 
	{
		delete[] neigh[i];
	}
	delete[] neigh;
	delete[] cfi;


}

