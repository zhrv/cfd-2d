#include "decomp.h"
#include "tinyxml.h"
#include <string>
#include "global.h"
#include "metis.h"
#include <string.h>
#include <vector>
#include <map>
#include <string>
#include "MeshReader.h"

//using namespace std;

typedef std::map<int, int> gl_map;
typedef std::vector<int> indexes;
typedef std::vector<indexes> exch_map;

struct ProcMesh
{
	int cCount, cCountEx;
	int eCount, eCountEx;
	int nCount, nCountEx;

	indexes gCells, gEdges, gNodes;
	gl_map lCells, lEdges, lNodes;

	indexes recvCount;
	exch_map sendInd;
};

const char FLAG_IN = 1;
const char FLAG_EX = 2;


void Decomp::init(char * xmlFileName)
{
	TiXmlDocument doc( xmlFileName );
	bool loadOkay = doc.LoadFile( TIXML_ENCODING_UTF8 );
	if (!loadOkay)
	{
		log("ERROR: %s\n", doc.ErrorDesc());
		exit(doc.ErrorId());
	}
	
	TiXmlNode* task = 0;
	TiXmlElement* el = 0;
	TiXmlNode* node0 = 0;
	TiXmlNode* node1 = 0;
	task = doc.FirstChild( "task" );


//	node0 = task->FirstChild("decomp");
//	node0->FirstChild("processors")->ToElement()->Attribute("value", &procCount);

	grids = new Grid[procCount];

	//const char* fName = task->FirstChild("mesh")->FirstChild("name")->ToElement()->Attribute("value");
	//grid.initFromFiles((char*)fName);

	/* Чтение данных сетки. */
	node0 = task->FirstChild("mesh");
	const char* fName = node0->FirstChild("name")->ToElement()->Attribute("value");
	const char* tName = node0->FirstChild("filesType")->ToElement()->Attribute("value");
	MeshReader* mr = MeshReader::create(MeshReader::getType((char*)tName), (char*)fName);
	mr->read(&grid);
}


void Decomp::run()
{
	char fName[128];

	log("Decomposition to %d processors started...\n", procCount);
	log(" Initial grid info:\n");
	log("  cells count:     %d\n", grid.cCount);
	log("  edges count:     %d\n", grid.eCount);
	log("  nodes count:     %d\n", grid.nCount);

	// декомпозиция области средствами METIS
	int n = grid.cCount;
	int m = grid.eCount;
	int ncon = 1;
	idx_t objval;
	idx_t * part   = new idx_t[n];
	if (procCount > 1) {
		idx_t * xadj = new idx_t[n + 1];
		idx_t * adjncy = new idx_t[2 * m];
		int jj = 0;
		xadj[0] = 0;
		for (int i = 0; i < n; i++)
		{
			Cell & cell = grid.cells[i];
			for (int j = 0; j < 3; j++)
			{
				if (cell.neigh[j] >= 0)
				{
					adjncy[jj] = cell.neigh[j];
					jj++;
				}
			}
			xadj[i + 1] = jj;
		}
		// TODO: подключить х64 версию METIS
		METIS_PartGraphRecursive(&n, &ncon, xadj, adjncy, NULL, NULL, NULL, &procCount, NULL, NULL, NULL, &objval, part);
		delete[] xadj;
		delete[] adjncy;
	}
	else {
		memset(part, 0, sizeof(idx_t)*n);
	}

	FILE * fp = fopen("parts.vtk", "w");
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "GASDIN data file\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp, "POINTS %d float\n", grid.nCount);
	for (int i = 0; i < grid.nCount; i++)
	{
		fprintf(fp, "%f %f %f  ", grid.nodes[i].x,  grid.nodes[i].y, 0.0);
		if (i+1 % 8 == 0) fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	fprintf(fp, "CELLS %d %d\n", grid.cCount, 4*grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "3 %d %d %d\n", grid.cells[i].nodesInd[0], grid.cells[i].nodesInd[1], grid.cells[i].nodesInd[2]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "CELL_TYPES %d\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++) fprintf(fp, "5\n");
	fprintf(fp, "\n");

	fprintf(fp, "CELL_DATA %d\nSCALARS Proc int 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%d ", part[i]);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fclose(fp);

	log(" Written file 'parts.vtk'...\n");

	// формирование файлов c данными сетки для каждого процессора
	MK_DIR('mesh'); // @todo: добавить удаление старого содержимого
	MK_DIR('vtk_data');  
	int * nProc = new int[procCount];
	memset(nProc, 0, procCount*sizeof(int));
	for (int i = 0; i < n; i++)
	{
		nProc[part[i]]++;
	}

	std::vector<ProcMesh> procMesh;
	for (int p = 0; p < procCount; p++)	{
		ProcMesh pm;
		int np = nProc[p];
		int ep = 0;
		//int * cellGlobalIdx = NULL;
		//int * cellFlag = NULL;
		indexes procCellEx;

		for (int i = 0; i < n; i++)
		{
			if (part[i] == p)
			{
				pm.gCells.push_back(i);
			}
		}
		for (int i = 0; i < pm.gCells.size(); i++)
		{
			Cell& c = grid.cells[pm.gCells[i]];
			for (int j = 0; j < 3; j++)
			{
				if (c.neigh[j] < 0) continue;
				if ((part[c.neigh[j]] != p)) procCellEx.push_back(c.neigh[j]);
			}
		}

		pm.cCount = pm.gCells.size();
		pm.cCountEx = pm.cCount + procCellEx.size();
		for (int i = 0; i < procCellEx.size(); i++) 
		{
			pm.gCells.push_back(procCellEx[i]);
		}
		procCellEx.clear();

		// сортировка фиктивных ячеек в порядке убывания процессора
		int exCnt = pm.cCountEx - pm.cCount;

		pm.recvCount.resize(procCount);
		pm.sendInd.resize(procCount);
		for (int i = 0; i < procCount; i++) {
			pm.recvCount[i] = 0;
			pm.sendInd[p].clear();
		}
		for (int ii = 0; ii < exCnt; ii++) {
			for (int j = 0; j < exCnt - 1 - ii; j++) {
				int i = pm.cCount + j;
				if (part[pm.gCells[i]] > part[pm.gCells[i + 1]]) {
					int tmp = pm.gCells[i];
					pm.gCells[i] = pm.gCells[i+1];
					pm.gCells[i+1] = tmp;
				}
			}
		}
		for (int i = 0; i < pm.cCountEx; i++) {
			pm.lCells[pm.gCells[i]] = i;
		}
		for (int i = pm.cCount; i < pm.cCountEx; i++) {
			pm.recvCount[part[pm.gCells[i]]]++;
		}


		char * edgeFlg = new char[grid.eCount];
		memset(edgeFlg, 0, grid.eCount*sizeof(char));
		for (int i = 0; i < pm.cCount; i++)
		{
			for (int j = 0; j < grid.cells[pm.gCells[i]].eCount; j++)
			{
				edgeFlg[grid.cells[pm.gCells[i]].edgesInd[j]] |= FLAG_IN;
			}
		}
		for (int i = pm.cCount; i < pm.cCountEx; i++)
		{
			for (int j = 0; j < grid.cells[pm.gCells[i]].eCount; j++)
			{
				edgeFlg[grid.cells[pm.gCells[i]].edgesInd[j]] |= FLAG_EX;
			}
		}

		for (int i = 0; i < grid.eCount; i++)
		{
			if (edgeFlg[i] & FLAG_IN) pm.gEdges.push_back(i);
		}
		pm.eCount = pm.gEdges.size();
		for (int i = 0; i < grid.eCount; i++)
		{
			if (((edgeFlg[i] & FLAG_EX) == FLAG_EX) && ((edgeFlg[i] & FLAG_IN) == 0)) pm.gEdges.push_back(i);
		}
		pm.eCountEx = pm.gEdges.size();
		delete[] edgeFlg;
		for (int i = 0; i < pm.eCountEx; i++) {
			pm.lEdges[pm.gEdges[i]] = i;
		}

		char * nodeFlg = new char[grid.nCount];
		memset(nodeFlg, 0, grid.nCount*sizeof(char));
		for (int i = 0; i < pm.cCount; i++)
		{
			for (int j = 0; j < grid.cells[pm.gCells[i]].nCount; j++)
			{
				nodeFlg[grid.cells[pm.gCells[i]].nodesInd[j]] |= FLAG_IN;
			}
		}
		for (int i = pm.cCount; i < pm.cCountEx; i++)
		{
			for (int j = 0; j < grid.cells[pm.gCells[i]].nCount; j++)
			{
				nodeFlg[grid.cells[pm.gCells[i]].nodesInd[j]] |= FLAG_EX;
			}
		}

		for (int i = 0; i < grid.nCount; i++) {
			if (nodeFlg[i] & FLAG_IN) pm.gNodes.push_back(i);
		}
		pm.nCount = pm.gNodes.size();
		for (int i = 0; i < grid.nCount; i++)
		{
			if (((nodeFlg[i] & FLAG_EX) == FLAG_EX) && ((nodeFlg[i] & FLAG_IN) == 0)) pm.gNodes.push_back(i);
		}
		pm.nCountEx = pm.gNodes.size();
		delete[] nodeFlg;
		for (int i = 0; i < pm.nCountEx; i++) {
			pm.lNodes[pm.gNodes[i]] = i;
		}


		procMesh.push_back(pm);
	}

	for (int p = 0; p < procCount; p++) {
		ProcMesh & pm = procMesh[p];

		for (int i = pm.cCount; i < pm.cCountEx; i++) {
			int giCell = pm.gCells[i];
			ProcMesh & pn = procMesh[part[giCell]];
			pn.sendInd[p].push_back(pn.lCells[giCell]);
		}
	}


	for (int p = 0; p < procCount; p++) {
		
		sprintf(fName, "mesh/mesh.%04d.proc", p);
		fp = fopen(fName, "w");
		ProcMesh & pm = procMesh[p];
		fprintf(fp, "%d %d\n", pm.nCount, pm.nCountEx);
		for (int i = 0; i < pm.nCountEx; i++) {
			int iNode = pm.gNodes[i];
			fprintf(fp, "%6d %25.15e %25.15e\n", i, grid.nodes[iNode].x, grid.nodes[iNode].y);
		}

		fprintf(fp, "\n%d %d\n", pm.cCount, pm.cCountEx);
		for (int i = 0; i < pm.cCountEx; i++) {
			Cell & c = grid.cells[pm.gCells[i]];
			fprintf(fp, "%6d %6d %6d %6d %s\n", i, pm.lNodes[c.nodesInd[0]], pm.lNodes[c.nodesInd[1]], pm.lNodes[c.nodesInd[2]], c.typeName);
		}

		fprintf(fp, "\n%d %d\n", pm.eCount, pm.eCountEx);
		for (int i = 0; i < pm.eCountEx; i++) {
			Edge & e = grid.edges[pm.gEdges[i]];
			fprintf(fp, "%6d %6d %6d %6d %s\n", i, pm.lNodes[e.n1], pm.lNodes[e.n2], e.type, e.type != 0 ? e.typeName : "");
		}
		
		fprintf(fp, "\n");
		for (int i = 0; i < procCount; i++) {
			fprintf(fp, "%d%s", pm.recvCount[i], (i+1 % 8 == 0 || i == procCount - 1) ? "\n" : " ");
		}
		fprintf(fp, "\n");

		for (int i = 0; i < procCount; i++) {
			int n = pm.sendInd[i].size();
			fprintf(fp, "%d %d\n", i, n);
			for (int j = 0; j < n; j++) {
				fprintf(fp, "%d%s", pm.sendInd[i][j], (j+1 % 8 == 0 || j == n - 1) ? "\n" : " ");
			}
		}

		fclose(fp);

	}
}

void Decomp::done()
{
	delete[] grids;
}
