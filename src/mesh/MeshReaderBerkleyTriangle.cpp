#include "MeshReaderBerkleyTriangle.h"

void MeshReaderBerkleyTriangle::read(Grid* g)
{
	char str[50];
	FILE *fp;
	int tmp;

	

	// читаем данные об ”«Ћј’
	sprintf(str, "%s.node", fileName);
	fp = fopen(str, "r");
	fscanf(fp, "%d %d %d %d", &g->nCount, &tmp, &tmp, &tmp);
	g->nodes = new Point[g->nCount];
	for (int i = 0; i < g->nCount; i++)
	{
		fscanf(fp, "%d %lf %lf %d", &tmp, &(g->nodes[i].x), &(g->nodes[i].y), &tmp);
	}
	fclose(fp);

	// читаем данные о я„≈… ј’
	sprintf(str, "%s.ele", fileName);
	fp = fopen(str, "r");
	fscanf(fp, "%d %d %d", &g->cCount, &tmp, &tmp);
	g->cells = new Cell[g->cCount];
	for (int i = 0; i < g->cCount; i++)
	{
		g->cells[i].nCount = 3;
		g->cells[i].nodesInd = new int[g->cells[i].nCount];
		fscanf(fp, "%d %d %d %d %d", &tmp, &(g->cells[i].nodesInd[0]), &(g->cells[i].nodesInd[1]), &(g->cells[i].nodesInd[2]), &(g->cells[i].type));
		g->cells[i].nodesInd[0]--;
		g->cells[i].nodesInd[1]--;
		g->cells[i].nodesInd[2]--;
		g->cells[i].c.x = (g->nodes[g->cells[i].nodesInd[0]].x + g->nodes[g->cells[i].nodesInd[1]].x + g->nodes[g->cells[i].nodesInd[2]].x) / 3.0;
		g->cells[i].c.y = (g->nodes[g->cells[i].nodesInd[0]].y + g->nodes[g->cells[i].nodesInd[1]].y + g->nodes[g->cells[i].nodesInd[2]].y) / 3.0;
		g->cells[i].HX = _max_(fabs(g->nodes[g->cells[i].nodesInd[0]].x - g->nodes[g->cells[i].nodesInd[1]].x),
			fabs(g->nodes[g->cells[i].nodesInd[1]].x - g->nodes[g->cells[i].nodesInd[2]].x),
			fabs(g->nodes[g->cells[i].nodesInd[0]].x - g->nodes[g->cells[i].nodesInd[2]].x));
		g->cells[i].HY = _max_(fabs(g->nodes[g->cells[i].nodesInd[0]].y - g->nodes[g->cells[i].nodesInd[1]].y),
			fabs(g->nodes[g->cells[i].nodesInd[1]].y - g->nodes[g->cells[i].nodesInd[2]].y),
			fabs(g->nodes[g->cells[i].nodesInd[0]].y - g->nodes[g->cells[i].nodesInd[2]].y));
		g->cells[i].eCount = 3;
		g->cells[i].edgesInd = new int[g->cells[i].eCount];
	}
	fclose(fp);

	// формируем данные о –≈Ѕ–ј’
	sprintf(str, "%s.neigh", fileName);
	fp = fopen(str, "r");
	fscanf(fp, "%d %d", &tmp, &tmp);
	int** neigh;
	neigh = new int*[g->cCount];
	for (int i = 0; i < g->cCount; i++)
	{
		neigh[i] = new int[3];
		fscanf(fp, "%d %d %d %d", &tmp, &(neigh[i][0]), &(neigh[i][1]), &(neigh[i][2]));
		neigh[i][0]--;
		neigh[i][1]--;
		neigh[i][2]--;
		g->cells[i].neigh[0] = neigh[i][0];
		g->cells[i].neigh[1] = neigh[i][1];
		g->cells[i].neigh[2] = neigh[i][2];
	}
	fclose(fp);
	g->eCount = 0;
	for (int i = 0; i < g->cCount; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int p = neigh[i][j];
			if (p > -1)
			{
				for (int k = 0; k < 3; k++)
				{ // убираем у соседа номер этой €чейки, чтобы грань не повтор€лась
					if (neigh[p][k] == i) neigh[p][k] = -1;
				}
				g->eCount++;
			}
			if (p == -2) g->eCount++;
		}
	}
	g->edges = new Edge[g->eCount];

	int iEdge = 0;
	int * cfi = new int[g->cCount];
	for (int i = 0; i < g->cCount; i++)
	{
		cfi[i] = 0;
	}
	// ::memset(cfi, 0, cCount*sizeof(int));
	for (int i = 0; i < g->cCount; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int p = neigh[i][j];
			if (p != -1)
			{
				g->edges[iEdge].n1 = g->cells[i].nodesInd[(j + 1) % 3];
				g->edges[iEdge].n2 = g->cells[i].nodesInd[(j + 2) % 3];
				g->edges[iEdge].cCount = 1;
				g->edges[iEdge].c = new Point[g->edges[iEdge].cCount];
				g->edges[iEdge].c[0].x = (g->nodes[g->edges[iEdge].n1].x + g->nodes[g->edges[iEdge].n2].x) / 2.0;
				g->edges[iEdge].c[0].y = (g->nodes[g->edges[iEdge].n1].y + g->nodes[g->edges[iEdge].n2].y) / 2.0;
				g->edges[iEdge].n.x = g->nodes[g->edges[iEdge].n2].y - g->nodes[g->edges[iEdge].n1].y;
				g->edges[iEdge].n.y = g->nodes[g->edges[iEdge].n1].x - g->nodes[g->edges[iEdge].n2].x;
				g->edges[iEdge].l = sqrt(g->edges[iEdge].n.x*g->edges[iEdge].n.x + g->edges[iEdge].n.y*g->edges[iEdge].n.y);
				g->edges[iEdge].n.x /= g->edges[iEdge].l;
				g->edges[iEdge].n.y /= g->edges[iEdge].l;
				g->edges[iEdge].c1 = i;
				g->cells[i].edgesInd[cfi[i]] = iEdge;
				cfi[i]++;
				g->edges[iEdge].cnl1 = fabs(g->edges[iEdge].n.x*(g->edges[iEdge].c[0].x - g->cells[g->edges[iEdge].c1].c.x) + g->edges[iEdge].n.y*(g->edges[iEdge].c[0].y - g->cells[g->edges[iEdge].c1].c.y));

				if (p > -1)
				{

					g->edges[iEdge].c2 = p;
					g->cells[p].edgesInd[cfi[p]] = iEdge;
					cfi[p]++;
					g->edges[iEdge].cnl2 = fabs(g->edges[iEdge].n.x*(g->cells[g->edges[iEdge].c2].c.x - g->edges[iEdge].c[0].x) + g->edges[iEdge].n.y*(g->cells[g->edges[iEdge].c2].c.y - g->edges[iEdge].c[0].y));
					g->edges[iEdge].type = Edge::TYPE_INNER;
				}
				if (p == -2)
				{
					g->edges[iEdge].c2 = -1;
					g->edges[iEdge].cnl2 = 0;
					g->edges[iEdge].type = -1;
				}
				iEdge++;
			}
		}
	}

	// чтение данных о граничных гран€х
	sprintf(str, "%s.poly", fileName);
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
		int iEdge = g->findEdge(n1, n2);
		if (iEdge >= 0) g->edges[iEdge].type = type;
	}
	fclose(fp);

	for (int i = 0; i < g->cCount; i++)
	{
		double a = g->edges[g->cells[i].edgesInd[0]].l;
		double b = g->edges[g->cells[i].edgesInd[1]].l;
		double c = g->edges[g->cells[i].edgesInd[2]].l;
		double p = (a + b + c) / 2.0;
		g->cells[i].S = sqrt(p*(p - a)*(p - b)*(p - c));
	}

	for (int i = 0; i < g->cCount; i++)
	{
		delete[] neigh[i];
	}
	delete[] neigh;
	delete[] cfi;


}


