
#include "MeshReaderMsh2.h"
#include <fstream>
#include <algorithm>
#include <cstring>
#include <ctime>

/*
 * Read a line from a file. Obtained from:
 * http://stackoverflow.com/questions/314401/
 * how-to-read-a-line-from-the-console-in-c/314422#314422
 *
 * Using this avoids a dependence on IEEE Std 1003.1-2008 (``POSIX.1'') for the
 * getline function.
 */
char* MeshReaderMsh2::getLineUpper(ifstream &stream)
{
    char               *line = new char[1024], *linep = line;
    size_t              lenmax = 1024, len = lenmax;
    int                 c;

    for (;;) {
        c = stream.get();
        c = toupper (c);
        if (c == EOF && linep == line) {
            delete[] linep;
            return nullptr;
        }

        if (--len == 0) {
            char               *linen;

            len = lenmax;
            lenmax *= 2;

            linen = new char[lenmax];
            memcpy(linen, linep, sizeof(char)*len);

            line = linen + (line - linep);
            linep = linen;
        }
        if ((*line++ = c) == '\n')
            break;
    }
    *line = '\0';
    return linep;
}

typedef struct {
    char  name[128];
    int   id;
    int   dim;
} patch_t;

typedef struct {
    long long id;
    double x, y, z;
} node_t;

typedef struct {
    long long id, type, n[3];
} cell_t;

typedef struct {
    long long id, type, n[2];
} edge_t;


void addNeigh(Cell &c, int j)
{
    if (j == c.neigh[0] or j == c.neigh[1] or j == c.neigh[2]) return;
    for (int i = 0; i < 3; i++) {
        if (c.neigh[i] == -1) {
            c.neigh[i] == j;
        }
    }

}


void MeshReaderMsh2::read(Grid * g)
{
    int num_vertices = 0,
        num_patches = 0,
        num_trees = 0,
        num_trees_real = 0,
        tree;
    vector<patch_t> patches;
    patch_t patch;
    vector<node_t> nodes;
    node_t node;
    vector<edge_t> edges;
    edge_t edge;
    vector<cell_t> cells;
    cell_t cell;
    double          *vert;
    long long tmp_ll;
    double tmp_d;
    int tmp_i;
    streampos fpos;

    ifstream fin(fileName);
    if (!fin.is_open()) {
        throw Exception("Mesh file opening error.", Exception::FILE_OPENING_ERROR);
    }

    while (1) {
        char *line = getLineUpper(fin);

        if (line == NULL) {
            break;
        }

        if (line[0] == '$') {
            if (strstr(line, "$PHYSICALNAMES")) {
                delete[] line;
                line = getLineUpper(fin);
                sscanf(line, "%d", &num_patches);
                for (int i = 0; i < num_patches; i++) {
                    delete[] line;
                    line = getLineUpper(fin);
                    sscanf(line, "%d %d \"%[^\"]", &(patch.dim), &(patch.id), patch.name);
                    patches.push_back(patch);
                }
            }
            else if (strstr(line, "$NODES")) {
                delete[] line;
                line = getLineUpper(fin);
                sscanf(line, "%d", &g->nCount);
                g->nodes = new Point[g->nCount];
                for (int i = 0; i < g->nCount; i++) {
                    delete[] line;
                    line = getLineUpper(fin);
                    int retval = sscanf (line, "%lld %lf %lf %lf", &tmp_ll, &(g->nodes[i].x), &(g->nodes[i].y), &tmp_d);
                    if (retval != 4) {
                        delete[] line;
                        throw Exception("Premature end of file", Exception::FILE_OPENING_ERROR);
                    }
                }
            }
            else if(strstr(line, "$ELEMENTS")) {
                fpos = fin.tellg();
                delete[] line;
                long long id, type, el_count;
                line = getLineUpper(fin);
                sscanf(line, "%d", &el_count);
                for (int i = 0; i < el_count; i++) {
                    line = getLineUpper(fin);
                    int retval = sscanf (line, "%lld %lld", &id, &type);
                    if (retval != 2) {
                        delete[] line;
                        throw Exception("Premature end of file", Exception::FILE_OPENING_ERROR);
                    }

                    if (type == 2) { // 3-node triangle
                        long long tagc, tag2;
                        retval = sscanf (line, "%lld %lld %lld %lld %lld %lld %lld %lld", &cell.id, &type, &tagc, &cell.type, &tag2,
                                         &(cell.n[0]), &(cell.n[1]), &(cell.n[2]) );
                        if (retval != 8) {
                            delete[] line;
                            throw Exception("Premature end of file", Exception::FILE_OPENING_ERROR);
                        }
                        --cell.n[0];
                        --cell.n[1];
                        --cell.n[2];
                        cells.push_back(cell);
                    }
                    else if (type == 1) { // 2-node line
                        long long tagc, tag2;
                        retval = sscanf (line, "%lld %lld %lld %lld %lld %lld %lld", &edge.id, &type, &tagc, &edge.type, &tag2,
                                         &(edge.n[0]), &(edge.n[1]) );
                        if (retval != 7) {
                            delete[] line;
                            throw Exception("Premature end of file", Exception::FILE_OPENING_ERROR);
                        }
                        --edge.n[0];
                        --edge.n[1];
                        if (edge.n[0] > edge.n[1]) {
                            tmp_i = edge.n[0];
                            edge.n[0] = edge.n[1];
                            edge.n[1] = tmp_i;
                        }
                        edges.push_back(edge);
                    }
                }
            }
        }

        delete[] line;
    }

    log("Building mesh structure:\n");

    // nodes
    log("\t- nodes;\n");
    g->nCount = nodes.size();
    g->nodes = new Point[g->nCount];
    int i = 0;
    for (auto it = nodes.begin(); it != nodes.end(); it++, i++) {
        Point & p = g->nodes[i];
        p.x = it->x;
        p.y = it->y;
    }

    log("\t- cells;\n");
    g->cCount = cells.size();
    g->cells = new Cell[g->cCount];

    map<int, ind_set> node_cells;
    map<pair<int, int>, vector<int>> edges_cells;
    i = 0;
    for (auto it = cells.begin(); it != cells.end(); it++, i++) {
        Cell & c = g->cells[i];
        c.nCount = 3;
        c.nodesInd = new int[g->cells[i].nCount];
        c.nodesInd[0] = it->n[0];
        c.nodesInd[1] = it->n[1];
        c.nodesInd[2] = it->n[2];

        g->cells[i].type = it->type;
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

        for (int j = 0; j < 3; j++) {
            int n1 = c.nodesInd[j];
            int n2 = c.nodesInd[(j+1)%3];
            if (n1 > n2) {
                tmp_i = n1;
                n1 = n2;
                n2 = tmp_i;
            }
            edges_cells[pair<int,int>(n1, n2)].push_back(i);
        }
    }


}