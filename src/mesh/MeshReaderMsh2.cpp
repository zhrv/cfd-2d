#include "MeshReaderMsh2.h"
#include <iostream>
#include <fstream> // подключаем файлы
#include <algorithm>
#include <cstring>


//int MeshReaderMsh2::findEdge(int n1, int n2)
//{
//    std::set <int> s1 = { n1, n2};
//    std::set <int> s2;
//    for (int i = all_edges.size() - 1; i > 0; --i) {
//        s2 = { all_edges[i][0], all_edges[i][1] };
//        if (s1 == s2) {
//            return i;
//        }
//        s2.clear();
//    }
//    return -1;
//}

int MeshReaderMsh2::find_edge(int n1, int n2)
{
    for (int i = 0; i < edges.size(); i++) {
        if (((edges[i][4] == n1) && (edges[i][5] == n2)) || ((edges[i][4] == n2) && (edges[i][5] == n1)) ){
            return i;
        }
    }
    return -1;
}


void MeshReaderMsh2::read(Grid* g)
{
    char str[64];
    long long int tmp, elem, cell, edge, type;
    double tmp1;
    int numPatches, numElements, retval;
    std::string line;

    // читаем сеточные данные
    std::sprintf(str, "%s.msh", fileName);
    log("Reading file '%s'\n", str);
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
                        retval = sscanf(line.c_str(), "%lld %lld %d %d %d %d %d", &edge, &type, &tagc, &tag1, &tag2,
                                        &node1, &node2);
                        if (retval != 7) {
                            log("Premature end of file");
                            exit(2);
                        }
                        v.push_back(type);
                        v.push_back(tagc);
                        v.push_back(tag1);
                        v.push_back(tag2);
                        v.push_back(node1-1);
                        v.push_back(node2-1);
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
                        v.push_back(node1-1);
                        v.push_back(node2-1);
                        v.push_back(node3-1);
                        cells.push_back(v);
                    }
                }
            }
        }
    }
    file.close();

    log("Building mesh structure:\n");

    std::map<int, ind_set> node_cells;
    // читаем данные о ЯЧЕЙКАХ
    g->cCount = cells.size();
    g->cells = new Cell[g->cCount];
    int** neigh;
    neigh = new int*[g->cCount];
    int i = 0;
    for (index_list::iterator it = cells.begin(); it != cells.end(); it++, i++) {
        indexes & ind = *it;
        node_cells[ind[4]].insert(i);
        node_cells[ind[5]].insert(i);
        node_cells[ind[6]].insert(i);
    }

    log("\t- cells;\n");
    i = 0;
    for (auto it = cells.begin(); it != cells.end(); it++, i++)
    {
        neigh[i] = new int[3];
        indexes & ind = *it;

//        elements[ind[0]] = element(i, element::TYPE_CELL);

        Cell & c = g->cells[i];
        c.nCount = 3;
        c.nodesInd = new int[g->cells[i].nCount];
        c.nodesInd[0] = ind[4];
        c.nodesInd[1] = ind[5];
        c.nodesInd[2] = ind[6];

        g->cells[i].type = cells[i][2];
        strcpy(g->cells[i].typeName, patches[g->cells[i].type-1].c_str());

        c.c.x = (g->nodes[g->cells[i].nodesInd[0]].x + g->nodes[g->cells[i].nodesInd[1]].x + g->nodes[g->cells[i].nodesInd[2]].x) / 3.0;
        c.c.y = (g->nodes[g->cells[i].nodesInd[0]].y + g->nodes[g->cells[i].nodesInd[1]].y + g->nodes[g->cells[i].nodesInd[2]].y) / 3.0;
        c.HX = _max_(fabs(g->nodes[g->cells[i].nodesInd[0]].x - g->nodes[g->cells[i].nodesInd[1]].x),
                               fabs(g->nodes[g->cells[i].nodesInd[1]].x - g->nodes[g->cells[i].nodesInd[2]].x),
                               fabs(g->nodes[g->cells[i].nodesInd[0]].x - g->nodes[g->cells[i].nodesInd[2]].x));
        c.HY = _max_(fabs(g->nodes[g->cells[i].nodesInd[0]].y - g->nodes[g->cells[i].nodesInd[1]].y),
                               fabs(g->nodes[g->cells[i].nodesInd[1]].y - g->nodes[g->cells[i].nodesInd[2]].y),
                               fabs(g->nodes[g->cells[i].nodesInd[0]].y - g->nodes[g->cells[i].nodesInd[2]].y));
        c.eCount = 3;
        c.edgesInd = new int[g->cells[i].eCount];

        for (int k = 0; k < 3; k++) {
            //map<int, ind_set>::iterator out_it;
            indexes res(50);
            auto it_res = set_intersection(
                    node_cells[c.nodesInd[k % 3]].begin(), node_cells[c.nodesInd[k % 3]].end(),
                    node_cells[c.nodesInd[(k + 1) % 3]].begin(), node_cells[c.nodesInd[(k + 1) % 3]].end(),
                    res.begin());
            res.resize(it_res-res.begin());
            c.neigh[k] = -2;
            for (auto rit = res.begin(); rit != res.end(); rit++) {
                if (*rit != i) {
                    c.neigh[k] = *rit;
                }
            }
            neigh[i][k] = c.neigh[k];
        }
    }


    log("\t- edges;\n");
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
                g->edges[iEdge].n1 = g->cells[i].nodesInd[(j + 0) % 3];
                g->edges[iEdge].n2 = g->cells[i].nodesInd[(j + 1) % 3];
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

                // коррекция направлений нормалей
                Vector vc;
                vc.x = g->cells[i].c.x - g->edges[iEdge].c[0].x;
                vc.y = g->cells[i].c.y - g->edges[iEdge].c[0].y;
                if (scalar_prod(vc, g->edges[iEdge].n) > 0) {
                    g->edges[iEdge].n.x *= -1;
                    g->edges[iEdge].n.y *= -1;
                }

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
                    g->edges[iEdge].type = Edge::TYPE_NAMED;

                    int jEdge = find_edge(g->edges[iEdge].n1, g->edges[iEdge].n2);
                    if (jEdge != -1) {
                        //elements[edges[jEdge][0]] = element(iEdge, element::TYPE_EDGE);
                        strcpy(g->edges[jEdge].typeName, patches[edges[jEdge][2]-1].c_str());
                    }
                    else {
                        throw Exception("Boundary edge #%d not defined in UNV file.", Exception::TYPE_MESH_UNV_NOT_DEFINED_BND_EDGE);
                    }
                }


                iEdge++;
            }
        }
    }

    for (int i = 0; i < g->cCount; i++)
    {
        double a = g->edges[g->cells[i].edgesInd[0]].l;
        double b = g->edges[g->cells[i].edgesInd[1]].l;
        double c = g->edges[g->cells[i].edgesInd[2]].l;
        double p = (a + b + c) / 2.0;
        g->cells[i].S = sqrt(p*(p - a)*(p - b)*(p - c));
    }

    for (int i = 0; i < g->eCount; i++) {
        Edge &e = g->edges[i];
        if (e.type == Edge::TYPE_NAMED and strcmp(e.typeName, "") == 0) {
            int kkk=0;
        }
    }
    for (int i = 0; i < g->cCount; i++)
    {
        delete[] neigh[i];
    }
    delete[] neigh;
    delete[] cfi;
}


