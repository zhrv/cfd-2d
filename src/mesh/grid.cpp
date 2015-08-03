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

typedef std::vector<int> idx_t;
typedef std::vector<idx_t> idx2_t;
typedef std::pair<int, int> edge_t;
typedef std::map<edge_t, int> edge_idx_t;

edge_idx_t nodeToEdge;
int findEdgeByNodes(int n1, int n2)
{
	if (n1 > n2) {
		int tmp = n1;
		n1 = n2;
		n2 = tmp;
	}
	return nodeToEdge[edge_t(n1, n2)];
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
		nodeToEdge[edge_t(n1, n2)] = i;
	}

	idx2_t edge_cells;
	edge_cells.resize(eCountEx);
	for (int i = 0; i < cCountEx; i++) {
		for (int j = 0; j < 3; j++) {
			int n1 = cells[i].nodesInd[j % 3];
			int n2 = cells[i].nodesInd[(j + 1) % 3];
			int iEdge = findEdgeByNodes(n1,n2);
			edge_cells[iEdge].push_back(i);
			cells[i].edgesInd[j] = iEdge;
		}
	}

	for (int iEdge = 0; iEdge < eCountEx; iEdge++) {
		idx_t & ec = edge_cells[iEdge];
		Edge & e = edges[iEdge];
		if (ec.size() == 2) {
			e.c1 = ec[0];
			e.c2 = ec[1];
		}
		else if (ec.size() == 1) {
			e.c1 = ec[0];
			e.c2 = -1;
		}
		else {
			log("ERROR: edge #%d is not assigned to any cell. Sourse: &s; line: %d\n", iEdge, __FILE__, __LINE__);
			EXIT(1);
		}
		e.cCount = 3;
		e.c = new Point[e.cCount];
		double _sqrt3 = 1.0 / sqrt(3.0);
		// центр ребра
		e.c[0].x = (nodes[e.n1].x + nodes[e.n2].x) / 2.0;
		e.c[0].y = (nodes[e.n1].y + nodes[e.n2].y) / 2.0;
		// перва€ точка √аусса
		e.c[1].x = (nodes[e.n1].x + nodes[e.n2].x) / 2.0 - _sqrt3*(nodes[e.n2].x - nodes[e.n1].x) / 2.0;
		e.c[1].y = (nodes[e.n1].y + nodes[e.n2].y) / 2.0 - _sqrt3*(nodes[e.n2].y - nodes[e.n1].y) / 2.0;
		// втора€ точка √аусса
		e.c[2].x = (nodes[e.n1].x + nodes[e.n2].x) / 2.0 + _sqrt3*(nodes[e.n2].x - nodes[e.n1].x) / 2.0;
		e.c[2].y = (nodes[e.n1].y + nodes[e.n2].y) / 2.0 + _sqrt3*(nodes[e.n2].y - nodes[e.n1].y) / 2.0;
		// нормаль и длина ребра
		e.n.x = nodes[e.n2].y - nodes[e.n1].y;
		e.n.y = nodes[e.n1].x - nodes[e.n2].x;
		e.l = sqrt(e.n.x*e.n.x + e.n.y*e.n.y);
		e.n.x /= e.l;
		e.n.y /= e.l;
		e.cnl1 = fabs(e.n.x*(e.c[0].x - cells[e.c1].c.x) + e.n.y*(e.c[0].y - cells[e.c1].c.y));
		if (e.c2 > -1) {
			e.cnl2 = fabs(e.n.x*(cells[e.c2].c.x - e.c[0].x) + e.n.y*(cells[e.c2].c.y - e.c[0].y));
		}
		else {
			e.cnl2 = 0.0;
		}
		// коррекци€ направлений нормалей
		Vector vc;
		
		vc.x = cells[e.c1].c.x - e.c[0].x;
		vc.y = cells[e.c1].c.y - e.c[0].y;
		if (scalar_prod(vc, e.n) > 0) {
			e.n.x *= -1;
			e.n.y *= -1;
		}
	}

	
	for (int i = 0; i < cCountEx; i++)	{
		double a = edges[cells[i].edgesInd[0]].l;
		double b = edges[cells[i].edgesInd[1]].l;
		double c = edges[cells[i].edgesInd[2]].l;
		double p = (a + b + c) / 2.0;
		cells[i].S = sqrt(p*(p - a)*(p - b)*(p - c));

	}

	recvCount.clear();
	int sh = cCount;
	for (int i = 0; i < Parallel::procCount; i++) {
		fscanf(fp, "%d", &tmp);
		recvCount.push_back(tmp);
		recvShift.push_back(sh);
		sh += tmp;
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

void Grid::saveMeshInfo() {
	char name[64];
	FILE * fp;

	sprintf(name, "mesh.%04d.nodes.info", Parallel::procId);
	fp = fopen(name, "w");
	fprintf(fp, "COUNT: %6d    COUNT_EX: %6d\n", nCount, nCountEx);
	for (int i = 0; i < nCount; i++) {
		Point & p = nodes[i];
		fprintf(fp, "============================================================\n");
		fprintf(fp, "#%06d\n", i);
		fprintf(fp, "POINT: %25.15e %25.15e\n", p.x, p.y);
		fprintf(fp, "============================================================\n\n\n");
	}
	fprintf(fp, "\n********************* Extended *********************\n", cCount, cCountEx);
	for (int i = nCount; i < nCountEx; i++) {
		Point & p = nodes[i];
		fprintf(fp, "============================================================\n");
		fprintf(fp, "#%06d\n", i);
		fprintf(fp, "POINT: %25.15e %25.15e\n", p.x, p.y);
		fprintf(fp, "============================================================\n\n\n");
	}
	fclose(fp);

	sprintf(name, "mesh.%04d.cells.info", Parallel::procId);
	fp = fopen(name, "w");
	fprintf(fp, "COUNT: %6d    COUNT_EX: %6d\n", cCount, cCountEx);
	for (int i = 0; i < cCount; i++) {
		Cell & c = cells[i];
		fprintf(fp, "============================================================\n");
		fprintf(fp, "#%06d\n", i);
		fprintf(fp, "NODES: %6d %6d %6d\n", c.nodesInd[0], c.nodesInd[1], c.nodesInd[2]);
		fprintf(fp, "EDGES: %6d %6d %6d\n", c.edgesInd[0], c.edgesInd[1], c.edgesInd[2]);
		fprintf(fp, "CENTER: %25.15e %25.15e\n", c.c.x, c.c.y);
		fprintf(fp, "HX: %25.15e   HY: %25.15e\n", c.HX, c.HY);
		fprintf(fp, "SQUARE: %25.15e\n", c.S);
		fprintf(fp, "TYPE: %s\n", c.typeName);
		fprintf(fp, "============================================================\n\n\n");
	}
	fprintf(fp, "\n********************* Extended *********************\n", cCount, cCountEx);
	for (int i = cCount; i < cCountEx; i++) {
		Cell & c = cells[i];
		fprintf(fp, "============================================================\n");
		fprintf(fp, "#%06d\n", i);
		fprintf(fp, "NODES: %6d %6d %6d\n", c.nodesInd[0], c.nodesInd[1], c.nodesInd[2]);
		fprintf(fp, "EDGES: %6d %6d %6d\n", c.edgesInd[0], c.edgesInd[1], c.edgesInd[2]);
		fprintf(fp, "CENTER: %25.15e %25.15e\n", c.c.x, c.c.y);
		fprintf(fp, "HX: %25.15e   HY: %25.15e\n", c.HX, c.HY);
		fprintf(fp, "SQUARE: %25.15e\n", c.S);
		fprintf(fp, "TYPE: %s\n", c.typeName);
		fprintf(fp, "============================================================\n\n\n");
	}
	fclose(fp);

	sprintf(name, "mesh.%04d.edges.info", Parallel::procId);
	fp = fopen(name, "w");
	fprintf(fp, "COUNT: %6d    COUNT_EX: %6d\n", eCount, eCountEx);
	for (int i = 0; i < eCount; i++) {
		Edge & e = edges[i];
		fprintf(fp, "============================================================\n");
		fprintf(fp, "#%06d\n", i);
		fprintf(fp, "NODES: %6d %6d\n", e.n1, e.n2);
		fprintf(fp, "CELLS: %6d %6d\n", e.c1, e.c2);
		fprintf(fp, "CENTER: %25.15e %25.15e\n", e.c[0].x, e.c[0].y);
		fprintf(fp, "LENGTH: %25.15e\n", e.l);
		fprintf(fp, "NORMAL: %25.15e %25.15e\n", e.n.x, e.n.y);
		fprintf(fp, "TYPE: %s\n", e.typeName);
		fprintf(fp, "============================================================\n\n\n");
	}
	fprintf(fp, "\n********************* Extended *********************\n", cCount, cCountEx);
	for (int i = eCount; i < eCountEx; i++) {
	}
	fclose(fp);
}