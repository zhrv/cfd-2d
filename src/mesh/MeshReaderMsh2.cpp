#include "MeshReaderMsh2.h"
#include <iostream>
#include <fstream> // подключаем файлы
#include <algorithm>
#include <cstring>



MeshReaderMsh2::edge_t make_edge(int n1, int n2) {
    if (n1 > n2) std::swap(n1, n2);
    MeshReaderMsh2::edge_t e(n1, n2);
    return e;
}

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
    int numPatches, numElements, retval, i;
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

    log("\t- cells;\n");
    g->cCount = cells.size();
    g->cells = new Cell[g->cCount];
    edge_cells_t e2c;
    i = 0;
    for (auto ind : cells)
    {
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
            e2c[make_edge(c.nodesInd[k % 3], c.nodesInd[(k + 1) % 3])].push_back(i);
        }

        ++i;
    }


    log("\t- edges;\n");
    g->eCount = e2c.size();
    g->edges = new Edge[g->eCount];
    int iEdge = 0;
    int * cfi = new int[g->cCount];
    for (int i = 0; i < g->cCount; i++)
    {
        cfi[i] = 0;
    }
    for (auto e : e2c) {
        edge_t e_nodes = e.first;
        indexes &c = e.second;

        g->edges[iEdge].n1 = e_nodes.first;
        g->edges[iEdge].n2 = e_nodes.second;
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
        if (c.size() != 2 and c.size() != 1) {
            throw Exception("Edge hasn'n nearest cells...", Exception::TYPE_MESH_GMSH_WRONG_EDGE);
        }
        i = c[0];
        g->edges[iEdge].c1 = i;
        g->cells[i].edgesInd[cfi[i]] = iEdge;
        cfi[i]++;
        g->edges[iEdge].cnl1 = fabs(
                g->edges[iEdge].n.x * (g->edges[iEdge].c[0].x - g->cells[g->edges[iEdge].c1].c.x) +
                g->edges[iEdge].n.y * (g->edges[iEdge].c[0].y - g->cells[g->edges[iEdge].c1].c.y));
        if (c.size() == 2) {
            i = c[1];
            g->edges[iEdge].c2 = i;
            g->cells[i].edgesInd[cfi[i]] = iEdge;
            cfi[i]++;
            g->edges[iEdge].cnl2 = fabs(
                    g->edges[iEdge].n.x * (g->edges[iEdge].c[0].x - g->cells[g->edges[iEdge].c2].c.x) +
                    g->edges[iEdge].n.y * (g->edges[iEdge].c[0].y - g->cells[g->edges[iEdge].c2].c.y));
            g->edges[iEdge].type = Edge::TYPE_INNER;
        }
        else {
            g->edges[iEdge].c2 = -1;
            g->edges[iEdge].cnl2 = 0;
            g->edges[iEdge].type = Edge::TYPE_NAMED;

            int jEdge = find_edge(g->edges[iEdge].n1, g->edges[iEdge].n2);
            if (jEdge == -1) {
                throw Exception("Boundary edge not defined...", Exception::TYPE_MESH_GMSH_WRONG_EDGE);
            }
            indexes &ind = edges[jEdge];
            strcpy(g->edges[iEdge].typeName, patches[ind[2]-1].c_str());
        }
        // коррекция направлений нормалей
        Vector vc;
        vc.x = g->cells[g->edges[iEdge].c1].c.x - g->edges[iEdge].c[0].x;
        vc.y = g->cells[g->edges[iEdge].c1].c.y - g->edges[iEdge].c[0].y;
        if (scalar_prod(vc, g->edges[iEdge].n) > 0) {
            g->edges[iEdge].n.x *= -1;
            g->edges[iEdge].n.y *= -1;
        }


        ++iEdge;
    }

    log("\t- neighbors;\n");
    for (int iCell = 0; iCell < g->cCount; iCell++) {
        Cell &cell = g->cells[iCell];
        for (int k = 0; k < 3; k++) {
            Edge &edge = g->edges[cell.edgesInd[k]];
            cell.neigh[k] = (edge.c1 == iCell ? edge.c2 : edge.c1);
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
}


