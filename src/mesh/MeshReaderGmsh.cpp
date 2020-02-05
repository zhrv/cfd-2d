#include "MeshReaderGmsh.h"
#include <iostream>
#include <fstream> // подключаем файлы
#include <algorithm>
#include <cstring>

std::vector<int> MeshReaderGmsh::getIntersection(index_list &sets)
{
	std::vector <int> result;  // To store the reaultant set 
	int smallSetInd = 0;  // Initialize index of smallest set 
	int minSize = sets[0].size(); // Initialize size of smallest set 

	// sort all the sets, and also find the smallest set 
	for (int i = 1; i < sets.size(); i++)
	{
		// sort this set 
		std::sort(sets[i].begin(), sets[i].end());

		// update minSize, if needed 
		if (minSize > sets[i].size())
		{
			minSize = sets[i].size();
			smallSetInd = i;
		}
	}

	std::map<int, int> elementsMap;

	// Add all the elements of smallest set to a map, if already present, 
	// update the frequency 
	for (int i = 0; i < sets[smallSetInd].size(); i++)
	{
		if (elementsMap.find(sets[smallSetInd][i]) == elementsMap.end())
			elementsMap[sets[smallSetInd][i]] = 1;
		else
			elementsMap[sets[smallSetInd][i]]++;
	}

	// iterate through the map elements to see if they are present in 
	// remaining sets 
	std::map<int, int>::iterator it;
	for (it = elementsMap.begin(); it != elementsMap.end(); ++it)
	{
		int elem = it->first;
		int freq = it->second;

		bool bFound = true;

		// Iterate through all sets 
		for (int j = 0; j < sets.size(); j++)
		{
			// If this set is not the smallest set, then do binary search in it 
			if (j != smallSetInd)
			{
				// If the element is found in this set, then find its frequency 
				if (binary_search(sets[j].begin(), sets[j].end(), elem))
				{
					int lInd = lower_bound(sets[j].begin(), sets[j].end(), elem)
						- sets[j].begin();
					int rInd = upper_bound(sets[j].begin(), sets[j].end(), elem)
						- sets[j].begin();

					// Update the minimum frequency, if needed 
					if ((rInd - lInd) < freq)
						freq = rInd - lInd;
				}
				// If the element is not present in any set, then no need  
				// to proceed for this element. 
				else
				{
					bFound = false;
					break;
				}
			}
		}

		// If element was found in all sets, then add it to result 'freq' times 
		if (bFound)
		{
			for (int k = 0; k < freq; k++)
				result.push_back(elem);
		}
	}
	return result;
}

int MeshReaderGmsh::findEdge(int n1, int n2)
{
	std::set <int> s1 = { n1, n2};
	std::set <int> s2;
    for (int i = all_edges.size() - 1; i > 0; --i) {
        s2 = { all_edges[i][0], all_edges[i][1] };
        if (s1 == s2) {
            return i;
        }
        s2.clear();
    }
	return -1;
}

void MeshReaderGmsh::read(Grid* g)
{
	char str[64];
	long long int tmp, elem, cell, edge, type;
	double tmp1;
	int numPatches, numElements, retval;
	std::string line;

	// читаем сеточные данные
	std::sprintf(str, "%s.msh", fileName);
	std::ifstream file(str);
	string_list patches;
	while (getline(file, line)) {
		if (line[0] == '$') {
			if (line == "$PhysicalNames") {
				getline(file, line);
				sscanf(line.c_str(), "%d", &numPatches);
				patches.resize(numPatches);
				for (int i = 0; i < numPatches; i++) {
					getline(file, line);
					int dim, id;
					char name[128];
					sscanf(line.c_str(), "%d %d \"%[^\"]", &(dim), &(id), name);
					patches[--id] = name;
				}
			}
			else if (line == "$Nodes") {
				getline(file, line);
				sscanf(line.c_str(), "%d", &(g->nCount));
				g->nodes = new Point[g->nCount];
				for (int i = 0; i < g->nCount; i++) {
					getline(file, line);
					sscanf(line.c_str(), "%lld %lf %lf %lf", &tmp, &g->nodes[i].x, &g->nodes[i].y, &tmp1);
				}
			}
			else if (line == "$Elements") {
				getline(file, line);
				sscanf(line.c_str(), "%d", &numElements);
				for (int i = 0; i < numElements; i++) {
					getline(file, line);
					retval = sscanf(line.c_str(), "%lld %lld", &elem, &type);
					if (retval != 2) {
						log("Premature end of file");
						exit(2);
					}

					if (type == 1) {
						indexes v;
						v.clear();

						int node1, node2;
						int tagc, tag1, tag2;
						retval = sscanf(line.c_str(), "%lld %lld %d %d %d %d %d", &edge, &type, &tagc, &tag1, &tag2, &node1, &node2);
						if (retval != 7) {
							log("Premature end of file");
							exit(2);
						}
						v.push_back(type);
						v.push_back(tagc);
						v.push_back(tag1);
						v.push_back(tag2);
						v.push_back(node1);
						v.push_back(node2);
						edges.push_back(v);
					}

					if (type == 2) {
						indexes v;
						v.clear();
						int node1, node2, node3;
						int tagc, tag1, tag2;
						retval = sscanf(line.c_str(), "%lld %lld %d %d %d %d %d %d", &cell, &type, &tagc, &tag1, &tag2,
							&node1, &node2, &node3);
						if (retval != 8) {
							log("Premature end of file");
							exit(2);
						}
						v.push_back(type);
						v.push_back(tagc);
						v.push_back(tag1);
						v.push_back(tag2);
						v.push_back(node1);
						v.push_back(node2);
						v.push_back(node3);
						cells.push_back(v);
					}
				}
			}
		}
	}
	file.close();
	
	std::map<int, ind_set> node_cells;
	// читаем данные о ЯЧЕЙКАХ
	int iCell = 0;
	g->cCount = cells.size();
	g->cells = new Cell[g->cCount];
	for (index_list::iterator it = cells.begin(); it != cells.end(); it++)
	{
		g->cells[iCell].nCount = 3;
		g->cells[iCell].nodesInd = new int[g->cells[iCell].nCount];
		g->cells[iCell].nodesInd[0] = --(*it)[4];
		g->cells[iCell].nodesInd[1] = --(*it)[5];
		g->cells[iCell].nodesInd[2] = --(*it)[6];
		g->cells[iCell].type        =   (*it)[2];
		//g->cells[iCell].neigh = new int[3];

		node_cells[g->cells[iCell].nodesInd[0]].insert(iCell);
		node_cells[g->cells[iCell].nodesInd[1]].insert(iCell);
		node_cells[g->cells[iCell].nodesInd[2]].insert(iCell);
	
		//sprintf(g->cells[iCell].typeName, "%d", g->cells[iCell].type);
		strcpy(g->cells[iCell].typeName, patches[g->cells[iCell].type-1].c_str());
		g->cells[iCell].c.x = (g->nodes[g->cells[iCell].nodesInd[0]].x + g->nodes[g->cells[iCell].nodesInd[1]].x + g->nodes[g->cells[iCell].nodesInd[2]].x) / 3.;
		g->cells[iCell].c.y = (g->nodes[g->cells[iCell].nodesInd[0]].y + g->nodes[g->cells[iCell].nodesInd[1]].y + g->nodes[g->cells[iCell].nodesInd[2]].y) / 3.;
		
		g->cells[iCell].HX = _max_(_max_(g->nodes[g->cells[iCell].nodesInd[0]].x, g->nodes[g->cells[iCell].nodesInd[1]].x), g->nodes[g->cells[iCell].nodesInd[2]].x) -
			_min_(_min_(g->nodes[g->cells[iCell].nodesInd[0]].x, g->nodes[g->cells[iCell].nodesInd[1]].x), g->nodes[g->cells[iCell].nodesInd[2]].x);
		g->cells[iCell].HY = _max_(_max_(g->nodes[g->cells[iCell].nodesInd[0]].y, g->nodes[g->cells[iCell].nodesInd[1]].y), g->nodes[g->cells[iCell].nodesInd[2]].y) -
			_min_(_min_(g->nodes[g->cells[iCell].nodesInd[0]].y, g->nodes[g->cells[iCell].nodesInd[1]].y), g->nodes[g->cells[iCell].nodesInd[2]].y);
	
		g->cells[iCell].eCount = 3;
		g->cells[iCell].edgesInd = new int[g->cells[iCell].eCount];

		iCell++;
	}
	
	// определяем соседние ячейки
	int** neigh;
	neigh = new int*[g->cCount];
	for (int i = 0; i < g->cCount; i++)
	{
		neigh[i] = new int[3];
		Cell & c = g->cells[i];
		for (int k = 0; k < 3; k++) {
			index_list sets;
			
			sets.push_back(indexes(node_cells[c.nodesInd[(k + 1) % 3]].begin(), node_cells[c.nodesInd[(k + 1) % 3]].end()));
			sets.push_back(indexes(node_cells[c.nodesInd[(k + 2) % 3]].begin(), node_cells[c.nodesInd[(k + 2) % 3]].end()));

			indexes res = getIntersection(sets);
			sets.clear();

			c.neigh[k] = -2;
			for (indexes::iterator rit = res.begin(); rit != res.end(); rit++) {
				if (*rit != i) {
					c.neigh[k] = *rit;
				}
			}
			neigh[i][k] = c.neigh[k];
		}
		g->cells[i].neigh[0] = neigh[i][0];
		g->cells[i].neigh[1] = neigh[i][1];
		g->cells[i].neigh[2] = neigh[i][2];
	}
	
	// формируем данные о РЕБРАХ
	g->eCount = 0;
	for (int i = 0; i < g->cCount; i++)
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
				g->eCount++;
			}
			if (p == -2) g->eCount++;
		}
	}
	g->eCountEx = g->eCount;
	g->edges = new Edge[g->eCount];

	int iEdge = 0;
	int * cfi = new int[g->cCount];

	::memset(cfi, 0, g->cCount*sizeof(int));
	for (int i = 0; i < g->cCount; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int p = neigh[i][j];
			if (p != -1)
			{
				g->edges[iEdge].n1 = g->cells[i].nodesInd[(j + 1) % 3];
				g->edges[iEdge].n2 = g->cells[i].nodesInd[(j + 2) % 3];

				std::vector<int> v;
				v.clear();
				v.push_back(g->edges[iEdge].n1);
				v.push_back(g->edges[iEdge].n2);
				all_edges.push_back(v);

                g->edges[iEdge].cCount = 3;
                g->edges[iEdge].c = new Point[g->edges[iEdge].cCount];
                double _sqrt3 = 1.0 / sqrt(3.0);
                // центр ребра
                g->edges[iEdge].c[0].x = (g->nodes[g->edges[iEdge].n1].x + g->nodes[g->edges[iEdge].n2].x) / 2.0;
                g->edges[iEdge].c[0].y = (g->nodes[g->edges[iEdge].n1].y + g->nodes[g->edges[iEdge].n2].y) / 2.0;
                // первая точка Гаусса
                g->edges[iEdge].c[1].x = (g->nodes[g->edges[iEdge].n1].x + g->nodes[g->edges[iEdge].n2].x) / 2.0 - _sqrt3*(g->nodes[g->edges[iEdge].n2].x - g->nodes[g->edges[iEdge].n1].x) / 2.0;
                g->edges[iEdge].c[1].y = (g->nodes[g->edges[iEdge].n1].y + g->nodes[g->edges[iEdge].n2].y) / 2.0 - _sqrt3*(g->nodes[g->edges[iEdge].n2].y - g->nodes[g->edges[iEdge].n1].y) / 2.0;
                // вторая точка Гаусса
                g->edges[iEdge].c[2].x = (g->nodes[g->edges[iEdge].n1].x + g->nodes[g->edges[iEdge].n2].x) / 2.0 + _sqrt3*(g->nodes[g->edges[iEdge].n2].x - g->nodes[g->edges[iEdge].n1].x) / 2.0;
                g->edges[iEdge].c[2].y = (g->nodes[g->edges[iEdge].n1].y + g->nodes[g->edges[iEdge].n2].y) / 2.0 + _sqrt3*(g->nodes[g->edges[iEdge].n2].y - g->nodes[g->edges[iEdge].n1].y) / 2.0;
				
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
					sprintf(g->edges[iEdge].typeName, "%d", g->edges[iEdge].type);
				}
				if (p == -2)
				{
					g->edges[iEdge].c2 = -1;
					g->edges[iEdge].cnl2 = 0;
					g->edges[iEdge].type = Edge::TYPE_NAMED;
					// коррекция направлений нормалей
					Vector vc;
					vc.x = g->cells[i].c.x - g->edges[iEdge].c[0].x;
					vc.y = g->cells[i].c.y - g->edges[iEdge].c[0].y;
					if (scalar_prod(vc, g->edges[iEdge].n) > 0) {
						g->edges[iEdge].n.x *= -1.0;
						g->edges[iEdge].n.y *= -1.0;
					}
				}
				iEdge++;
			}
		}
	}

	for (index_list::iterator it = edges.begin(); it != edges.end(); it++)
	{
		int iEdge = findEdge(--(*it)[4], --(*it)[5]);
		if (iEdge != -1) {
			g->edges[iEdge].type = Edge::TYPE_NAMED;
			//sprintf(g->edges[iEdge].typeName, "%s", patches[(*it)[2] - 1].c_str());
            strcpy(g->edges[iEdge].typeName, patches[(*it)[2]-1].c_str());
		}
		else {
			throw Exception((char*)"Boundary edge #%d not defined in msh file.", Exception::TYPE_MESH_GMSH_NOT_DEFINED_BND_EDGE);
		}
	}

	/*for (int i = 0; i < g->eCount; i++) {
		printf("%d %s %d %d\n", g->edges[i].type, g->edges[i].typeName, g->edges[i].c1, g->edges[i].c2);
	}*/

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


