#include "grid.h"
#include <map>
#include <set>
#include <algorithm>

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

//Grid::initFromFiles: загрузка сетки из файла fName.
void Grid::initFromFiles(char* fName) 
{
	char str[50];
	FILE *fp;
	int tmp; 

	// читаем данные об ”«Ћј’
	sprintf(str, "%s.node", fName);
	fp = fopen(str, "r");
	if (!fp) {
		log("Can not open file '%s'\n", str);
		EXIT(1);
	}
	fscanf(fp, "%d %d %d %d", &nCount, &tmp, &tmp, &tmp);
	nodes = new Point[nCount];
	for (int i = 0; i < nCount; i++) 
	{
		fscanf(fp, "%d %lf %lf %d", &tmp, &(nodes[i].x), &(nodes[i].y), &tmp);
	}
	fclose(fp);

	// читаем данные о я„≈… ј’
	sprintf(str, "%s.ele", fName);
	fp = fopen(str, "r");
	if (!fp) {
		log("Can not open file '%s'\n", str);
		EXIT(1);
	}
	fscanf(fp, "%d %d %d", &cCount, &tmp, &tmp);
	cells = new Cell[cCount];
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
	}
	fclose(fp);

	// формируем данные о –≈Ѕ–ј’
	sprintf(str, "%s.neigh", fName);
	fp = fopen(str, "r");
	if (!fp) {
		log("Can not open file '%s'\n", str);
		EXIT(1);
	}
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
				{ // убираем у соседа номер этой €чейки, чтобы грань не повтор€лась
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
				// перва§ точка vаусса
				edges[iEdge].c[1].x = (nodes[edges[iEdge].n1].x+nodes[edges[iEdge].n2].x)/2.0-_sqrt3*(nodes[edges[iEdge].n2].x-nodes[edges[iEdge].n1].x)/2.0;
				edges[iEdge].c[1].y = (nodes[edges[iEdge].n1].y+nodes[edges[iEdge].n2].y)/2.0-_sqrt3*(nodes[edges[iEdge].n2].y-nodes[edges[iEdge].n1].y)/2.0;
				// втора§ точка vаусса
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
				iEdge++;
			}
		}
		
	}
	
	// чтение данных о граничных гран€х
	sprintf(str, "%s.poly", fName);
	fp = fopen(str, "r");
	if (!fp) {
		log("Can not open file '%s'\n", str);
		EXIT(1);
	}
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
		cells[i].flag = CELL_FLAG_GOOD;		
	}

	for (int i = 0; i < cCount; i++) 
	{
		cells[i].flag = 0;
	}
	
	for (int i = 0; i < cCount; i++) 
	{
		delete[] neigh[i];
	}
	delete[] neigh;
	delete[] cfi;


	for (int i = 0; i < eCount; i++) {

	}


}


std::map< std::pair<int, int>, int > nodeToEdge;
int findEdgeByNodes(int n1, int n2)
{
	if (n1 > n2) {
		int tmp = n1;
		n1 = n2;
		n2 = tmp;
	}
	return nodeToEdge[std::pair<int, int>(n1, n2)];
}


void Grid::readMeshFiles()
{
	int p = Parallel::procId;
	char fName[64];
	int tmp;
	sprintf(fName, "mesh/mesh.%04d.proc", p);
	FILE * fp = fopen(fName, "r");
	if (!fp) {
		log("Can not open file '%s'\n", fName);
		EXIT(1);
	}

	log("Reading mesh part structure:\n");

	// Nodes
	log("\t- nodes;\n");
	fscanf(fp, "%d %d", &nCount, &nCountEx);
	nodes = new Point[nCountEx];
	for (int i = 0; i < nCountEx; i++)
	{
		fscanf(fp, "%d %lf %lf", &tmp, &(nodes[i].x), &(nodes[i].y));
	}

	// Cells
	log("\t- cells;\n");
	fscanf(fp, "%d %d", &cCount, &cCountEx);
	cells = new Cell[cCountEx];
	for (int i = 0; i < cCountEx; i++)
	{
		cells[i].nCount = 3;
		cells[i].nodesInd = new int[cells[i].nCount];
		fscanf(fp, "%d %d %d %d %s", &tmp, &(cells[i].nodesInd[0]), &(cells[i].nodesInd[1]), &(cells[i].nodesInd[2]), &(cells[i].typeName));
		cells[i].c.x = (nodes[cells[i].nodesInd[0]].x + nodes[cells[i].nodesInd[1]].x + nodes[cells[i].nodesInd[2]].x) / 3.0;
		cells[i].c.y = (nodes[cells[i].nodesInd[0]].y + nodes[cells[i].nodesInd[1]].y + nodes[cells[i].nodesInd[2]].y) / 3.0;
		cells[i].HX = _max_(fabs(nodes[cells[i].nodesInd[0]].x - nodes[cells[i].nodesInd[1]].x),
			fabs(nodes[cells[i].nodesInd[1]].x - nodes[cells[i].nodesInd[2]].x),
			fabs(nodes[cells[i].nodesInd[0]].x - nodes[cells[i].nodesInd[2]].x));
		cells[i].HY = _max_(fabs(nodes[cells[i].nodesInd[0]].y - nodes[cells[i].nodesInd[1]].y),
			fabs(nodes[cells[i].nodesInd[1]].y - nodes[cells[i].nodesInd[2]].y),
			fabs(nodes[cells[i].nodesInd[0]].y - nodes[cells[i].nodesInd[2]].y));
		cells[i].eCount = 3;
		cells[i].edgesInd = new int[cells[i].eCount];
	}
	std::map<int, std::set<int> > node_cells;

	for (int i = 0; i < cCountEx; i++) {
		node_cells[cells[i].nodesInd[0]].insert(i);
		node_cells[cells[i].nodesInd[1]].insert(i);
		node_cells[cells[i].nodesInd[2]].insert(i);
	}

	i = 0;
	int** neigh = new int*[cCountEx];
	for (int i = 0; i < cCountEx; i++) {
		Cell & c = cells[i];
		neigh[i] = new int[3];

		for (int k = 0; k < 3; k++) {
			std::map<int, std::set<int> >::iterator out_it;
			std::vector<int> res(10);
			std::vector<int>::iterator it_res = set_intersection(
				node_cells[c.nodesInd[k % 3]].begin(), node_cells[c.nodesInd[k % 3]].end(),
				node_cells[c.nodesInd[(k + 1) % 3]].begin(), node_cells[c.nodesInd[(k + 1) % 3]].end(),
				res.begin());
			res.resize(it_res - res.begin());
			c.neigh[k] = -2;
			for (std::vector<int>::iterator rit = res.begin(); rit != res.end(); rit++) {
				if (*rit != i) {
					c.neigh[k] = *rit;
				}
			}
			neigh[i][k] = c.neigh[k];
		}
	}

	// Edges
	log("\t- edges;\n");
	fscanf(fp, "%d %d", &eCount, &eCountEx);
	edges = new Edge[eCountEx];
	for (int i = 0; i < eCountEx; i++)
	{
		fscanf(fp, "%d %d %d %d", &tmp, &(edges[i].n1), &(edges[i].n2), &(edges[i].type));
		if (edges[i].type != 0) {
			fscanf(fp, "%s", &(edges[i].typeName));
		}
		else {
			strcpy(edges[i].typeName, "");
		}
		int n1, n2;
		if (edges[i].n1 < edges[i].n2) {
			n1 = edges[i].n1;
			n2 = edges[i].n2;
		}
		else {
			n1 = edges[i].n2;
			n2 = edges[i].n1;
		}
		std::pair<int, int> nodes(n1, n2);
		nodeToEdge[nodes] = i;
	}

	int iEdge = 0;
	int * cfi = new int[cCountEx];
	for (int i = 0; i < cCountEx; i++)
	{
		cfi[i] = 0;
	}
	// ::memset(cfi, 0, cCount*sizeof(int));
	for (int i = 0; i < cCountEx; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int p = neigh[i][j];
			if (p != -1)
			{
				iEdge = findEdgeByNodes(cells[i].nodesInd[(j + 0) % 3], cells[i].nodesInd[(j + 1) % 3]);
				//g->edges[iEdge].cCount = 1;
				//g->edges[iEdge].c = new Point[g->edges[iEdge].cCount];
				//g->edges[iEdge].c[0].x = (g->nodes[g->edges[iEdge].n1].x + g->nodes[g->edges[iEdge].n2].x) / 2.0;
				//g->edges[iEdge].c[0].y = (g->nodes[g->edges[iEdge].n1].y + g->nodes[g->edges[iEdge].n2].y) / 2.0;
				edges[iEdge].cCount = 3;
				edges[iEdge].c = new Point[edges[iEdge].cCount];
				double _sqrt3 = 1.0 / sqrt(3.0);
				// центр ребра
				edges[iEdge].c[0].x = (nodes[edges[iEdge].n1].x + nodes[edges[iEdge].n2].x) / 2.0;
				edges[iEdge].c[0].y = (nodes[edges[iEdge].n1].y + nodes[edges[iEdge].n2].y) / 2.0;
				// перва€ точка √аусса
				edges[iEdge].c[1].x = (nodes[edges[iEdge].n1].x + nodes[edges[iEdge].n2].x) / 2.0 - _sqrt3*(nodes[edges[iEdge].n2].x - nodes[edges[iEdge].n1].x) / 2.0;
				edges[iEdge].c[1].y = (nodes[edges[iEdge].n1].y + nodes[edges[iEdge].n2].y) / 2.0 - _sqrt3*(nodes[edges[iEdge].n2].y - nodes[edges[iEdge].n1].y) / 2.0;
				// втора€ точка √аусса
				edges[iEdge].c[2].x = (nodes[edges[iEdge].n1].x + nodes[edges[iEdge].n2].x) / 2.0 + _sqrt3*(nodes[edges[iEdge].n2].x - nodes[edges[iEdge].n1].x) / 2.0;
				edges[iEdge].c[2].y = (nodes[edges[iEdge].n1].y + nodes[edges[iEdge].n2].y) / 2.0 + _sqrt3*(nodes[edges[iEdge].n2].y - nodes[edges[iEdge].n1].y) / 2.0;
				// нормаль и длина ребра
				edges[iEdge].n.x = nodes[edges[iEdge].n2].y - nodes[edges[iEdge].n1].y;
				edges[iEdge].n.y = nodes[edges[iEdge].n1].x - nodes[edges[iEdge].n2].x;
				edges[iEdge].l = sqrt(edges[iEdge].n.x*edges[iEdge].n.x + edges[iEdge].n.y*edges[iEdge].n.y);
				edges[iEdge].n.x /= edges[iEdge].l;
				edges[iEdge].n.y /= edges[iEdge].l;
				edges[iEdge].c1 = i;
				cells[i].edgesInd[cfi[i]] = iEdge;
				cfi[i]++;
				edges[iEdge].cnl1 = fabs(edges[iEdge].n.x*(edges[iEdge].c[0].x - cells[edges[iEdge].c1].c.x) + edges[iEdge].n.y*(edges[iEdge].c[0].y - cells[edges[iEdge].c1].c.y));

				// коррекци€ направлений нормалей
				Vector vc;
				vc.x = cells[i].c.x - edges[iEdge].c[0].x;
				vc.y = cells[i].c.y - edges[iEdge].c[0].y;
				if (scalar_prod(vc, edges[iEdge].n) > 0) {
					edges[iEdge].n.x *= -1;
					edges[iEdge].n.y *= -1;
				}

				if (p > -1)
				{

					edges[iEdge].c2 = p;
					cells[p].edgesInd[cfi[p]] = iEdge;
					cfi[p]++;
					edges[iEdge].cnl2 = fabs(edges[iEdge].n.x*(cells[edges[iEdge].c2].c.x - edges[iEdge].c[0].x) + edges[iEdge].n.y*(cells[edges[iEdge].c2].c.y - edges[iEdge].c[0].y));
					edges[iEdge].type = Edge::TYPE_INNER;
				}
				if (p == -2)
				{
					edges[iEdge].c2 = -1;
					edges[iEdge].cnl2 = 0;
					edges[iEdge].type = Edge::TYPE_NAMED;

					//int jEdge = find_edge(g->edges[iEdge].n1, g->edges[iEdge].n2);
					//if (jEdge != -1) {
					//	elements[edges[jEdge][0]] = element(iEdge, element::TYPE_EDGE);
					//}
					//else {
					//	throw Exception("Boundary edge #%d not defined in UNV file.", Exception::TYPE_MESH_UNV_NOT_DEFINED_BND_EDGE);
					//}
				}


			}
		}
	}
	
	for (int i = 0; i < cCountEx; i++)	{
		double a = edges[cells[i].edgesInd[0]].l;
		double b = edges[cells[i].edgesInd[1]].l;
		double c = edges[cells[i].edgesInd[2]].l;
		double p = (a + b + c) / 2.0;
		cells[i].S = sqrt(p*(p - a)*(p - b)*(p - c));

		delete[] neigh[i];
	}

	delete[] neigh;
	delete[] cfi;

	recvCount.clear();
	for (int i = 0; i < Parallel::procCount; i++) {
		fscanf(fp, "%d", &tmp);
		recvCount.push_back(tmp);
	}
	
	sendInd.clear();
	for (int i = 0; i < Parallel::procCount; i++) {
		int np;
		fscanf(fp, "%d %d", &tmp, &np);
		std::vector<int> ind;
		ind.clear();
		for (int j = 0; j < np; j++) {
			fscanf(fp, "%d", &tmp);
			ind.push_back(tmp);
		}
		sendInd.push_back(ind);
	}

	fclose(fp);


	log("  complete...\n");


}