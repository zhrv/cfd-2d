#include "decomp.h"
#include "tinyxml.h"
#include <string>
#include "global.h"
#include "metis.h"
#include <string.h>
#include <vector>

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


	node0 = task->FirstChild("decomp");
	node0->FirstChild("processors")->ToElement()->Attribute("value", &procCount);

	grids = new Grid[procCount];

	const char* fName = task->FirstChild("mesh")->FirstChild("name")->ToElement()->Attribute("value");
	grid.initFromFiles((char*)fName);
}


void Decomp::run()
{
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
	idx_t * xadj   = new idx_t[n+1];
	idx_t * adjncy = new idx_t[2*m];
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
		xadj[i+1] = jj;
	}
	// TODO: подключить х64 версию METIS
	//METIS_PartGraphRecursive(&n, &ncon, xadj, adjncy, NULL, NULL, NULL, &procCount, NULL, NULL, NULL, &objval, part);
	delete[] xadj;
	delete[] adjncy;


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
	int * nProc = new int[procCount];
	memset(nProc, 0, procCount*sizeof(int));
	for (int i = 0; i < n; i++)
	{
		nProc[part[i]]++;
	}

	for (int p = 0; p < procCount; p++)
	{
		int np = nProc[p];
		int ep = 0;
		int * cellGlobalIdx = NULL;
		int * cellFlag = NULL;
		std::vector<int> procCell;
		std::vector<int> procCellEx;

		for (int i = 0; i < n; i++)
		{
			if (part[i] == p)
			{
				procCell.push_back(i);
			}
		}
		for (int i = 0; i < procCell.size(); i++)
		{
			Cell& c = grid.cells[procCell[i]];
			for (int j = 0; j < 3; j++)
			{
				if (c.neigh[j] < 0) continue;
				if ((part[c.neigh[j]] != p)) procCellEx.push_back(c.neigh[j]);
			}
		}

		int cellCount = procCell.size();
		int cellCountEx = cellCount + procCellEx.size();
		for (int i = 0; i < procCellEx.size(); i++) 
		{
			procCell.push_back(procCellEx[i]);
		}
		procCellEx.clear();

		char FLG_EDGE_IN = 1;
		char FLG_EDGE_EX = 2;
		char * edgeFlg = new char[grid.eCount];
		memset(edgeFlg, 0, grid.eCount*sizeof(char));
		for (int i = 0; i < cellCount; i++)
		{
			for (int j = 0; j < grid.cells[procCell[i]].eCount; j++)
			{
				edgeFlg[grid.cells[procCell[i]].edgesInd[j]] |= FLG_EDGE_IN;
			}
		}
		for (int i = cellCount; i < cellCountEx; i++)
		{
			for (int j = 0; j < grid.cells[procCell[i]].eCount; j++)
			{
				edgeFlg[grid.cells[procCell[i]].edgesInd[j]] |= FLG_EDGE_EX;
			}
		}

		std::vector<int> edgeProc;
		for (int i = 0; i < grid.eCount; i++)
		{
			if (edgeFlg[i] & FLG_EDGE_IN) edgeProc.push_back(i);
		}
		int edgeCount = edgeProc.size();
		for (int i = 0; i < grid.eCount; i++)
		{
			if ( ((edgeFlg[i] & FLG_EDGE_EX) == FLG_EDGE_EX) && ((edgeFlg[i] & FLG_EDGE_IN) == 0) ) edgeProc.push_back(i);
		}
		int edgeCountEx = edgeProc.size();

		char fName[20];
		sprintf(fName, "mesh/mesh.%04d.proc", p);
		fp = fopen(fName, "w");
		fprintf(fp, "%d\n", np);
		fclose(fp);
	}
}

void Decomp::done()
{
	delete[] grids;
}
