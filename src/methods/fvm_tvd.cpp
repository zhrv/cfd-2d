#include "fvm_tvd.h"
#include "tinyxml.h"
#include <string>
#include "global.h"
#include "MeshReader.h"

void FVM_TVD::init(char * xmlFileName)
{
	STEADY = true;
	
	
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

	int steadyVal = 1;
	node0 = task->FirstChild("control");
	node0->FirstChild("STEADY")->ToElement()->Attribute("value", &steadyVal);
	node0->FirstChild("TAU")->ToElement()->Attribute("value", &TAU);
	node0->FirstChild("TMAX")->ToElement()->Attribute("value", &TMAX);
	node0->FirstChild("CFL")->ToElement()->Attribute("value", &CFL);
	node0->FirstChild("STEP_MAX")->ToElement()->Attribute("value", &STEP_MAX);
	node0->FirstChild("FILE_OUTPUT_STEP")->ToElement()->Attribute("value", &FILE_SAVE_STEP);
	node0->FirstChild("LOG_OUTPUT_STEP")->ToElement()->Attribute("value", &PRINT_STEP);

	if (steadyVal == 0) {
		STEADY = false;
	} else {
		STEADY = true;
	}

	// чтение параметров о ПРЕДЕЛЬНЫХ ЗНАЧЕНИЯХ
	node0 = task->FirstChild("limits");
	node0->FirstChild("ro")->ToElement()->Attribute("min", &limitRmin);
	node0->FirstChild("ro")->ToElement()->Attribute("max", &limitRmax);
	node0->FirstChild("p")->ToElement()->Attribute( "min", &limitPmin);
	node0->FirstChild("p")->ToElement()->Attribute( "max", &limitPmax);
	node0->FirstChild("u")->ToElement()->Attribute( "max", &limitUmax);

	// чтение параметров о МАТЕРИАЛАХ
	node0 = task->FirstChild("materials");
	node0->ToElement()->Attribute("count", &matCount);;
	materials = new Material[matCount];
	TiXmlNode* matNode = node0->FirstChild("material"); 
	for (int i = 0; i < matCount; i++)
	{
		Material & mat = materials[i];
		matNode->ToElement()->Attribute("id", &mat.id);
		node1 = matNode->FirstChild("name");
		el = node1->ToElement();
		mat.name = el->GetText();
		node1 = matNode->FirstChild("parameters");
		node1->FirstChild( "M"  )->ToElement()->Attribute( "value", &mat.M  );
		node1->FirstChild( "Cp" )->ToElement()->Attribute( "value", &mat.Cp );
		node1->FirstChild( "K"  )->ToElement()->Attribute( "value", &mat.K  );
		node1->FirstChild( "ML" )->ToElement()->Attribute( "value", &mat.ML );
		matNode = matNode->NextSibling("material");
	}

	// чтение параметров о РЕГИОНАХ
	node0 = task->FirstChild("regions");
	node0->ToElement()->Attribute("count", &regCount);
	regions = new Region[regCount];
	TiXmlNode* regNode = node0->FirstChild("region");
	for (int i = 0; i < regCount; i++)
	{
		Region & reg = regions[i];
		regNode->ToElement()->Attribute("id", &reg.id);
		regNode->FirstChild("material")->ToElement()->Attribute("id", &reg.matId);
		regNode->FirstChild("cell")->ToElement()->Attribute("type", &reg.cellType);
		
		reg.name = regNode->FirstChild("name")->ToElement()->GetText();

		node1 = regNode->FirstChild("parameters");
		node1->FirstChild( "Vx" )->ToElement()->Attribute( "value", &reg.par.u );
		node1->FirstChild( "Vy" )->ToElement()->Attribute( "value", &reg.par.v );
		node1->FirstChild( "T"  )->ToElement()->Attribute( "value", &reg.par.T );
		node1->FirstChild( "P"  )->ToElement()->Attribute( "value", &reg.par.p );
		
		Material& mat = materials[reg.matId];
		mat.URS(reg.par, 2);	// r=r(p,T)
		mat.URS(reg.par, 1);	// e=e(p,r)
		
		regNode = regNode->NextSibling("region");
	}


	/* Чтение параметров ГУ */
	node0 = task->FirstChild("boundaries");
	TiXmlNode* bNode = node0->FirstChild("boundCond");
	while (bNode != NULL)
	{
		int edgeType;
		bNode->ToElement()->Attribute("edgeType", &edgeType);

		CFDBoundary * b;

		try {
			b = CFDBoundary::create(bNode, &grid);
		}
		catch (Exception e) {
			log("ERROR: %s\n", e.getMessage());
			exit(e.getType());
		}

		boundaries.push_back(b);

		bNode = bNode->NextSibling("boundCond");
	}

	bCount = boundaries.size();

	/* Чтение данных сетки. */
//	node0 = task->FirstChild("mesh");
//	const char* fName = node0->FirstChild("name")->ToElement()->Attribute("value");
//	const char* tName = node0->FirstChild("filesType")->ToElement()->Attribute("value");
//	MeshReader* mr = MeshReader::create(MeshReader::getType((char*)tName), (char*)fName);
//	mr->read(&grid);

    grid.readMeshFiles();

	/* Определение ГУ для каждого ребра. */
	for (int iEdge = 0; iEdge < grid.eCount; iEdge++) {
		Edge & e = grid.edges[iEdge];
		if (e.type == Edge::TYPE_INNER) {
			e.bnd = NULL;
			continue;
		}
		if (e.type == Edge::TYPE_NAMED) {
			int iBound = -1;
			for (int i = 0; i < bCount; i++)
			{
				if (strcmp(e.typeName, boundaries[i]->name) == 0)
				{
					iBound = i;
					break;
				}
			}
			if (iBound < 0)
			{
				log("ERROR (boundary condition): unknown edge type of edge %d...\n", iEdge);
				EXIT(1);
			}

			e.bnd = boundaries[iBound];
		}
		else {
			int iBound = -1;
			for (int i = 0; i < bCount; i++)
			{
				if (e.type == boundaries[i]->edgeType)
				{
					iBound = i;
					break;
				}
			}
			if (iBound < 0)
			{
				log("ERROR (boundary condition): unknown edge type of edge %d...\n", iEdge);
				EXIT(1);
			}

			e.bnd = boundaries[iBound];
		}
	}

    Parallel::buf = new double[100*grid.cCountEx];
    dBuf = new double[grid.cCountEx];
    iBuf = new int[grid.cCountEx];
    vBuf = new VECTOR[grid.cCountEx];
    pBuf = new Point[grid.cCountEx];

	cTau    = new double[grid.cCountEx];

	ro		= new double[grid.cCountEx];
	ru		= new double[grid.cCountEx];
	rv		= new double[grid.cCountEx];
	re		= new double[grid.cCountEx];

	ro_old	= new double[grid.cCountEx];
	ru_old	= new double[grid.cCountEx];
	rv_old	= new double[grid.cCountEx];
	re_old	= new double[grid.cCountEx];

	ro_int	= new double[grid.cCountEx];
	ru_int	= new double[grid.cCountEx];
	rv_int	= new double[grid.cCountEx];
	re_int	= new double[grid.cCountEx];

	gradR		= new Vector[grid.cCountEx];
	gradP		= new Vector[grid.cCountEx];
	gradU		= new Vector[grid.cCountEx];
    gradV		= new Vector[grid.cCountEx];
    gradT		= new Vector[grid.cCountEx];

	for (int i = 0; i < grid.cCount; i++)
	{
		Cell & c = grid.cells[i];
		Region & reg = getRegion(c.typeName);
		convertParToCons(i, reg.par);
	}
	procExchangeFields();

	memcpy(ro_old, ro, grid.cCountEx*sizeof(double));
	memcpy(ru_old, ru, grid.cCountEx*sizeof(double));
	memcpy(rv_old, rv, grid.cCountEx*sizeof(double));
	memcpy(re_old, re, grid.cCountEx*sizeof(double));

	calcTimeStep();
	save(0);
}


void FVM_TVD::calcTimeStep()
{
	if (STEADY) {
		for (int iCell = 0; iCell < grid.cCount; iCell++)
		{
			Param p;
			convertConsToPar(iCell, p);
			cTau[iCell] = CFL * grid.cells[iCell].S / _max_(abs(p.u) + p.cz, abs(p.v) + p.cz);
		}
	}
	else {
//		for (int iCell = 0; iCell < grid.cCount; iCell++)
//		{
//			Param p;
//			convertConsToPar(iCell, p);
//			if (TAU > CFL* grid.cells[iCell].S / _max_(abs(p.u) + p.cz, abs(p.v) + p.cz)) {
//				TAU = CFL * grid.cells[iCell].S / _max_(abs(p.u) + p.cz, abs(p.v) + p.cz);
//			}
//		}
		for (int iCell = 0; iCell < grid.cCount; iCell++) {
			cTau[iCell] = TAU;
		}
		log("time step: %25.16E\n", TAU);
	}

	exchange(cTau);
}

void FVM_TVD::calcGrad() 
{
	int nc = grid.cCountEx;
	int ne = grid.eCount;
	

	memset(gradR, 0, nc*sizeof(Vector));
	memset(gradP, 0, nc*sizeof(Vector));
	memset(gradU, 0, nc*sizeof(Vector));
    memset(gradV, 0, nc*sizeof(Vector));
    memset(gradT, 0, nc*sizeof(Vector));

    for (int iEdge = 0; iEdge < ne; iEdge++)
	{
			
		int c1	= grid.edges[iEdge].c1;
		int c2	= grid.edges[iEdge].c2;
			
		Param pL, pR;
		convertConsToPar(c1, pL);
		if (c2 > -1) {
			convertConsToPar(c2, pR);
		} else {
			boundaryCond(iEdge, pL, pR);
		}

		Vector n	= grid.edges[iEdge].n;
		double l	= grid.edges[iEdge].l;
			
			
			
		gradR[c1].x += (pL.r+pR.r)/2*n.x*l;
		gradR[c1].y += (pL.r+pR.r)/2*n.y*l;
		gradP[c1].x += (pL.p+pR.p)/2*n.x*l;
		gradP[c1].y += (pL.p+pR.p)/2*n.y*l;
		gradU[c1].x += (pL.u+pR.u)/2*n.x*l;
		gradU[c1].y += (pL.u+pR.u)/2*n.y*l;
        gradV[c1].x += (pL.v+pR.v)/2*n.x*l;
        gradV[c1].y += (pL.v+pR.v)/2*n.y*l;
        gradT[c1].x += (pL.T+pR.T)/2*n.x*l;
        gradT[c1].y += (pL.T+pR.T)/2*n.y*l;
		if (c2 > -1)
		{
			gradR[c2].x -= (pL.r+pR.r)/2*n.x*l;
			gradR[c2].y -= (pL.r+pR.r)/2*n.y*l;
			gradP[c2].x -= (pL.p+pR.p)/2*n.x*l;
			gradP[c2].y -= (pL.p+pR.p)/2*n.y*l;
			gradU[c2].x -= (pL.u+pR.u)/2*n.x*l;
			gradU[c2].y -= (pL.u+pR.u)/2*n.y*l;
			gradV[c2].x -= (pL.v+pR.v)/2*n.x*l;
			gradV[c2].y -= (pL.v+pR.v)/2*n.y*l;
            gradT[c2].x -= (pL.T+pR.T)/2*n.x*l;
            gradT[c2].y -= (pL.T+pR.T)/2*n.y*l;
		}

	}
	for (int iCell = 0; iCell < nc; iCell++)
	{
		double si = grid.cells[iCell].S;
		gradR[iCell] /= si;
		gradP[iCell] /= si;
		gradU[iCell] /= si;
        gradV[iCell] /= si;
        gradT[iCell] /= si;
	}

    procExchangeGrads();
}

void FVM_TVD::singleTimeStep()
{
    memset(ro_int, 0, grid.cCountEx*sizeof(double));
    memset(ru_int, 0, grid.cCountEx*sizeof(double));
    memset(rv_int, 0, grid.cCountEx*sizeof(double));
    memset(re_int, 0, grid.cCountEx*sizeof(double));
    for (int iEdge = 0; iEdge < grid.eCount; iEdge++)
    {
        Edge &edge = grid.edges[iEdge];
        double fr, fu, fv, fe;
        int c1	= edge.c1;
        int c2	= edge.c2;
        Vector n	= edge.n;
        double l	= edge.l*0.5;
        Material &mat = getMaterial(c1);
        double gam = mat.getGamma();
        double mu = mat.ML;
        double lambda = 0.;
        double kt = mat.K;
        Param pL, pR;
        fr = 0.0;
        fu = 0.0;
        fv = 0.0;
        fe = 0.0;
        for (int iGP = 1; iGP < edge.cCount; iGP++)
        {
            double fr1, fu1, fv1, fe1;
            reconstruct(iEdge, pL, pR, edge.c[iGP]);
            calcFlux(fr1, fu1, fv1, fe1, pL, pR, n, gam);
            fr += fr1;
            fu += fu1;
            fv += fv1;
            fe += fe1;
        }

        if (c2 > -1) {
            reconstruct(iEdge, pL, pR, edge.c[0]);
            double rl = sqrt(pL.r);
            double rr = sqrt(pR.r);

            double u_  = (rl*pL.u+rr*pR.u)/(rl+rr);
            double v_  = (rl*pL.v+rr*pR.v)/(rl+rr);
            double u_x = (gradU[c1].x+gradU[c2].x)*0.5;
            double v_x = (gradV[c1].x+gradV[c2].x)*0.5;
            double t_x = (gradT[c1].x+gradT[c2].x)*0.5;
            double u_y = (gradU[c1].y+gradU[c2].y)*0.5;
            double v_y = (gradV[c1].y+gradV[c2].y)*0.5;
            double t_y = (gradT[c1].y+gradT[c2].y)*0.5;

            double txx = (lambda-(2./3.)*mu)*(u_x+v_y)-2.*mu*u_x;
            double tyy = (lambda-(2./3.)*mu)*(u_x+v_y)-2.*mu*v_y;
            double txy = mu*(u_y+v_x);

            fu -= txx*n.x+txy*n.y;
            fv -= txy*n.x+tyy*n.y;
            fe -= (u_*txx+v_*txy+kt*t_x)*n.x+(u_*txy+v_*tyy+kt*t_y)*n.y;
        }
        else {
            if (instanceof<CFDBndWallNoSlip>(edge.bnd)) {
                reconstruct(iEdge, pL, pR, edge.c[0]);
                double Tbnd = edge.bnd->par[0];
                double vn = pL.u*n.x+pL.v*n.y;
                Vector Vn = n;
                Vn *= vn;
                Vector Vt = Vn;
                Vt -= Vn;

                double txx = mu*Vt.x/edge.cnl1;
                double tyy = mu*Vt.y/edge.cnl1;
                double Q   = kt*(Tbnd-pL.T)/edge.cnl1;

                fu -= txx;
                fv -= tyy;
                fe -= Q;
            }
        }

        ro_int[c1] -= fr*l;
        ru_int[c1] -= fu*l;
        rv_int[c1] -= fv*l;
        re_int[c1] -= fe*l;
        if (c2 > -1)
        {
            ro_int[c2] += fr*l;
            ru_int[c2] += fu*l;
            rv_int[c2] += fv*l;
            re_int[c2] += fe*l;
        }

    }

    for (int iCell = 0; iCell < grid.cCount; iCell++)
    {
        if (cellIsLim(iCell)) continue;
        double cfl = cTau[iCell]/grid.cells[iCell].S;
        ro[iCell] += cfl*ro_int[iCell];
        ru[iCell] += cfl*ru_int[iCell];
        rv[iCell] += cfl*rv_int[iCell];
        re[iCell] += cfl*re_int[iCell];
    }

    procExchangeFields();
}

void FVM_TVD::run() 
{
    double			t		= 0.0;
	int	            step	= 0;
	while (t <= TMAX && step <= STEP_MAX)
	{
		if (!STEADY) {
			t += TAU; 
		} else {
			calcTimeStep();
		}
		step++;
		memcpy(ro_old, ro, grid.cCount*sizeof(double));
		memcpy(ru_old, ru, grid.cCount*sizeof(double));
		memcpy(rv_old, rv, grid.cCount*sizeof(double));
		memcpy(re_old, re, grid.cCount*sizeof(double));

		// первый подшаг метода Р.-К.
		calcGrad();
        singleTimeStep();

		// второй подшаг метода Р.-К.
		calcGrad();
		singleTimeStep();

		// полусумма: формула (4.10) из icase-1997-65.pdf
		for (int iCell = 0; iCell < grid.cCount; iCell++)
		{
			if (cellIsLim(iCell)) continue;

			ro[iCell] = 0.5*(ro_old[iCell]+ro[iCell]);
			ru[iCell] = 0.5*(ru_old[iCell]+ru[iCell]);
			rv[iCell] = 0.5*(rv_old[iCell]+rv[iCell]);
			re[iCell] = 0.5*(re_old[iCell]+re[iCell]);

			Param par;
			convertConsToPar(iCell, par);
			if (par.r < limitRmin)			{ par.r = limitRmin; setCellFlagLim(iCell); }
			if (par.r > limitRmax)			{ par.r = limitRmax; setCellFlagLim(iCell); }
			if (par.p < limitPmin)			{ par.p = limitPmin; setCellFlagLim(iCell); }
			if (par.p > limitPmax)			{ par.p = limitPmax; setCellFlagLim(iCell); }
			if (fabs(par.u) > limitUmax)	{ par.u = limitUmax; setCellFlagLim(iCell); }
			if (fabs(par.v) > limitUmax)	{ par.v = limitUmax; setCellFlagLim(iCell); }
		}

        procExchangeFields();

		remediateLimCells();

        procExchangeFields();

		if (step % FILE_SAVE_STEP == 0)
		{
            Parallel::barrier();

            save(step);
		}
		if (step % PRINT_STEP == 0)
		{
			log("step: %d\t\ttime step: %.16f\n", step, t);
		}
	}

}

void FVM_TVD::remediateLimCells()
{
	for (int iCell = 0; iCell < grid.cCount; iCell++) 
	{
		if (cellIsLim(iCell)) 
		{
			// пересчитываем по соседям
			double sRO = 0.0;
			double sRU = 0.0;
			double sRV = 0.0;
			double sRE = 0.0;
			double S   = 0.0;
			for (int i = 0; i < grid.cells[iCell].eCount; i++)
			{
				int iEdge = grid.cells[iCell].edgesInd[i];
				int j = grid.edges[iEdge].c2;
				if (j >= 0) {
					double s = grid.cells[j].S;
					S += s;
					sRO += ro[j] * s;
					sRU += ru[j] * s;
					sRV += rv[j] * s;
					sRE += re[j] * s;
				} 
			}
			ro[iCell] = sRO/S;
			ru[iCell] = sRU/S;
			rv[iCell] = sRV/S;
			re[iCell] = sRE/S;

			// после 0x20 итераций пробуем вернуть ячейку в счет
			grid.cells[iCell].flag += 0x010000;
			if (grid.cells[iCell].flag & 0x200000) grid.cells[iCell].flag &= 0x001110;
		}
	}
}

void FVM_TVD::save(int step)
{
    char fName[50];
    FILE * fp;
    if (Parallel::isRoot()) {
        sprintf(fName, "res_%010d.pvtu", step);
        fp = fopen(fName, "w");
        fprintf(fp, "<?xml version=\"1.0\"?>\n");
        fprintf(fp, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
        fprintf(fp, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
        fprintf(fp, "    <PPoints><PDataArray type=\"Float32\" NumberOfComponents = \"3\" /></PPoints>\n");
        fprintf(fp, "    <PCellData Scalars=\"Density, Pressure, Total_pressure, Temperature, Mach, Proc, TAU\" Vectors=\"Velocity\">\n");
        fprintf(fp, "      <PDataArray type=\"Float32\" Name=\"Density\" format=\"ascii\" />\n");
        fprintf(fp, "      <PDataArray type=\"Float32\" Name=\"Pressure\" format=\"ascii\" />\n");
        fprintf(fp, "      <PDataArray type=\"Float32\" Name=\"Total_pressure\" format=\"ascii\" />\n");
        fprintf(fp, "      <PDataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\" />\n");
        fprintf(fp, "      <PDataArray type=\"Float32\" Name=\"Mach\" format=\"ascii\" />\n");
        fprintf(fp, "      <PDataArray type=\"Float32\" Name=\"Velocity\" format=\"ascii\" NumberOfComponents=\"3\"/>\n");
        fprintf(fp, "      <PDataArray type=\"Int32\" Name=\"Proc\" format=\"ascii\" />\n");
        fprintf(fp, "      <PDataArray type=\"Float32\" Name=\"TAU\" format=\"ascii\" />\n");
        fprintf(fp, "    </PCellData>\n");

        for (int pid = 0; pid < Parallel::procCount; pid++) {
            fprintf(fp, "    <Piece Source=\"res_%010d.%04d.vtu\"/>\n", step, pid);
        }
        fprintf(fp, "  </PUnstructuredGrid>\n");

        fprintf(fp, "</VTKFile>\n");
        fclose(fp);
        log("File '%s' saved...\n", fName);
    }
    sprintf(fName, "res_%010d.%04d.vtu", step, Parallel::procId);
    fp = fopen(fName, "w");
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(fp, "  <UnstructuredGrid GhostLevel=\"0\">\n");
    fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", grid.nCount, grid.cCount);
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    fprintf(fp, "          ");
    for (int i = 0; i < grid.nCount; i++) {
        fprintf(fp, "%25.15f %25.15f %25.15f ", grid.nodes[i].x, grid.nodes[i].y, 0.0);
    }
    fprintf(fp, "\n");
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Points>\n");
    fprintf(fp, "      <Cells>\n");
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
    fprintf(fp, "          ");
    for (int i = 0; i < grid.cCount; i++) {
        fprintf(fp, "%d ", (i + 1) * 3);
    }
    fprintf(fp, "\n");
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
    fprintf(fp, "          ");
    for (int i = 0; i < grid.cCount; i++) {
        fprintf(fp, "%d %d %d ", grid.cells[i].nodesInd[0], grid.cells[i].nodesInd[1], grid.cells[i].nodesInd[2]);
    }
    fprintf(fp, "\n");
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
    fprintf(fp, "          ");
    for (int i = 0; i < grid.cCount; i++) {
        fprintf(fp, "%d ", 5);
    }
    fprintf(fp, "\n");
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </Cells>\n");


    fprintf(fp, "      <CellData Scalars=\"Density, Pressure, Total_pressure, Temperature, Mach, Proc, TAU\" Vectors=\"Velocity\">\n");

    // плотность
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Density\" format=\"ascii\">\n");
    fprintf(fp, "          ");
    for (int i = 0; i < grid.cCount; i++) {
        Param p;
        convertConsToPar(i, p);
        fprintf(fp, "%25.16f ", p.r);
        if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");

    // давление
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Pressure\" format=\"ascii\">\n");
    fprintf(fp, "          ");
    for (int i = 0; i < grid.cCount; i++) {
        Param p;
        convertConsToPar(i, p);
        fprintf(fp, "%25.16f ", p.p);
        if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");

    // полное давление
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Total_pressure\" format=\"ascii\">\n");
    fprintf(fp, "          ");
    for (int i = 0; i < grid.cCount; i++) {
        Material &mat = getMaterial(i);
        double gam = mat.getGamma();
        double agam = gam - 1.0;
        Param p;
        convertConsToPar(i, p);
        double M2 = (p.u*p.u+p.v*p.v)/(gam*p.p/p.r);
        fprintf(fp, "%f ", p.p*::pow(1.0+0.5*M2*agam, gam/(gam-1.0)) );
        if ((i+1) % 8 == 0  ||  i+1 == grid.cCount) fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");

    // температура
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">\n");
    fprintf(fp, "          ");
    for (int i = 0; i < grid.cCount; i++) {
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.T);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");

    // число Маха
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Mach\" format=\"ascii\">\n");
    fprintf(fp, "          ");
    for (int i = 0; i < grid.cCount; i++) {
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", sqrt(p.u*p.u+p.v*p.v)/p.cz);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");

    // скорость
    fprintf(fp,
            "        <DataArray type=\"Float32\" Name=\"Velocity\" format=\"ascii\" NumberOfComponents=\"3\">\n");
    fprintf(fp, "          ");
    for (int i = 0; i < grid.cCount; i++) {
        Param p;
        convertConsToPar(i, p);
        fprintf(fp, "%25.16f %25.16f %25.16f ", p.u, p.v, 0.0);
        if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");

    // номер процессора
    fprintf(fp, "        <DataArray type=\"Int32\" Name=\"Proc\" format=\"ascii\">\n");
    fprintf(fp, "          ");
    for (int i = 0; i < grid.cCount; i++) {
        fprintf(fp, "%d ", Parallel::procId);
        if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");

    // шаг по времени
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"TAU\" format=\"ascii\">\n");
    fprintf(fp, "          ");
    for (int i = 0; i < grid.cCount; i++) {
        fprintf(fp, "%25.16f ", cTau[i]);
        if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
    }
    fprintf(fp, "        </DataArray>\n");
    fprintf(fp, "      </CellData>\n");


    fprintf(fp, "    </Piece>\n");


    fprintf(fp, "  </UnstructuredGrid>\n");

    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
}

void FVM_TVD::calcFlux(double& fr, double& fu, double& fv, double& fe, Param pL, Param pR, Vector n, double GAM)
{
	{	// GODUNOV FLUX
		double RI, EI, PI, UI, VI, WI, UN, UT;
		double unl = pL.u*n.x+pL.v*n.y;
		double unr = pR.u*n.x+pR.v*n.y;
		double utl = pL.u*n.y-pL.v*n.x;
		double utr = pR.u*n.y-pR.v*n.x;
		rim_orig(RI, EI, PI, UN, UT, WI,  pL.r, pL.p, unl, utl, 0,  pR.r, pR.p, unr, utr, 0, GAM);
		
		UI = UN*n.x+UT*n.y;
		VI = UN*n.y-UT*n.x;
		fr = RI*UN;
		fu = fr*UI+PI*n.x;
		fv = fr*VI+PI*n.y;
		fe = (RI*(EI+0.5*(UI*UI+VI*VI))+PI)*UN;
	}
// 	{	// LAX-FRIEDRIX FLUX
// 		double unl = pL.u*n.x+pL.v*n.y;
// 		double unr = pR.u*n.x+pR.v*n.y;
// 		double rol, rul, rvl, rel,  ror, rur, rvr, rer;
// 		double alpha = _max_(fabs(unl)+sqrt(GAM*pL.p/pL.r), fabs(unr)+sqrt(GAM*pR.p/pR.r));
// 		rol = pL.r;
// 		rul = pL.r*pL.u;
// 		rvl = pL.r*pL.v;
// 		rel = pL.r*pL.E;
// 		ror = pR.r;
// 		rur = pR.r*pR.u;
// 		rvr = pR.r*pR.v;
// 		rer = pR.r*pR.E;
// 		double frl = rol*unl;
// 		double frr = ror*unr;
// 		fr = 0.5*(frr+frl								- alpha*(ror-rol));
// 		fu = 0.5*(frr*pR.u+frl*pL.u + (pR.p+pL.p)*n.x	- alpha*(rur-rul));
// 		fv = 0.5*(frr*pR.v+frl*pL.v + (pR.p+pL.p)*n.y	- alpha*(rvr-rvl));
// 		fe = 0.5*((rer+pR.p)*unr + (rel+pL.p)*unl		- alpha*(rer-rel));
// 	}
}


void FVM_TVD::reconstruct(int iEdge, Param& pL, Param& pR, Point p)
{
	if (grid.edges[iEdge].type == Edge::TYPE_INNER) 
	{
		int c1	= grid.edges[iEdge].c1;
		int c2	= grid.edges[iEdge].c2;
		convertConsToPar(c1, pL);
		convertConsToPar(c2, pR);
		Point &PE = p;
		Point P1 = grid.cells[c1].c;
		Point P2 = grid.cells[c2].c;
		Vector DL1;
		Vector DL2;
		DL1.x=PE.x-P1.x;
		DL1.y=PE.y-P1.y;
		DL2.x=PE.x-P2.x;
		DL2.y=PE.y-P2.y;
		pL.r+=gradR[c1].x*DL1.x+gradR[c1].y*DL1.y;
		pL.p+=gradP[c1].x*DL1.x+gradP[c1].y*DL1.y;
		pL.u+=gradU[c1].x*DL1.x+gradU[c1].y*DL1.y;
		pL.v+=gradV[c1].x*DL1.x+gradV[c1].y*DL1.y;
		pR.r+=gradR[c2].x*DL2.x+gradR[c2].y*DL2.y;
		pR.p+=gradP[c2].x*DL2.x+gradP[c2].y*DL2.y;
		pR.u+=gradU[c2].x*DL2.x+gradU[c2].y*DL2.y;
		pR.v+=gradV[c2].x*DL2.x+gradV[c2].y*DL2.y;
	} else {
		int c1 = grid.edges[iEdge].c1;
		convertConsToPar(c1, pL);
		//return;
		Point &PE = p;
		Point P1 = grid.cells[c1].c;
		Vector DL1;
		DL1.x = PE.x - P1.x;
		DL1.y = PE.y - P1.y;
		pL.r += gradR[c1].x*DL1.x + gradR[c1].y*DL1.y;
		pL.p += gradP[c1].x*DL1.x + gradP[c1].y*DL1.y;
		pL.u += gradU[c1].x*DL1.x + gradU[c1].y*DL1.y;
		pL.v += gradV[c1].x*DL1.x + gradV[c1].y*DL1.y;
		boundaryCond(iEdge, pL, pR);
	}
}


void FVM_TVD::boundaryCond(int iEdge, Param& pL, Param& pR)
{
	Edge &edge = grid.edges[iEdge];
	int c1 = edge.c1;
	Material& m = getMaterial(c1);
	if (edge.bnd) {
		edge.bnd->run(iEdge, pL, pR);
		if (instanceof<CFDBndOutlet>(edge.bnd)) {
		    m.URS(pR, 4);
		}
		else {
            m.URS(pR, 2);
            m.URS(pR, 1);
        }
		pR.E = pR.e + 0.5*(pR.U2());
		return;
	}
	else {
		char msg[128];
		sprintf(msg, "Not defined boundary condition for edge %d\n", iEdge);
		throw Exception(msg, Exception::TYPE_BOUND_UNKNOWN);
	}
}


void FVM_TVD::done()
{
	delete[] ro;
	delete[] ru;
	delete[] rv;
	delete[] re;

	delete[] ro_old;
	delete[] ru_old;
	delete[] rv_old;
	delete[] re_old;

	delete[] ro_int;
	delete[] ru_int;
	delete[] rv_int;
	delete[] re_int;

	delete[] gradR;
	delete[] gradP;
	delete[] gradU;
    delete[] gradV;
    delete[] gradT;
}




Region & FVM_TVD::getRegionByCellType(int type)
{
	for (int i = 0; i < regCount; i++)
	{
		if (regions[i].cellType == type) return regions[i];
	}
	log("ERROR: unknown cell type %d...\n", type);
	EXIT(1);
}


Region   &	FVM_TVD::getRegion	(int iCell)
{
	return getRegionByCellType( grid.cells[iCell].type );
}

Region & FVM_TVD::getRegionByName(char* name)
{
	for (int i = 0; i < regCount; i++)
	{
		if (strcmp(regions[i].name.c_str(), name) == 0) return regions[i];
	}
	log("ERROR: unknown cell name '%d'...\n", name);
	EXIT(1);
}

Region & FVM_TVD::getRegion(char * name)
{
	return getRegionByName(name);
}

Material &	FVM_TVD::getMaterial(int iCell)
{
	Region & reg = getRegion(grid.cells[iCell].typeName);
	return materials[reg.matId];
}


void FVM_TVD::convertParToCons(int iCell, Param & par)
{
	ro[iCell] = par.r;
	ru[iCell] = par.r*par.u;
	rv[iCell] = par.r*par.v;
	re[iCell] = par.r*(par.e+0.5*(par.u*par.u+par.v*par.v));
}

void FVM_TVD::convertConsToPar(int iCell, Param & par)
{
	par.r = ro[iCell];
	par.u = ru[iCell]/ro[iCell];
	par.v = rv[iCell]/ro[iCell];
	par.E = re[iCell]/ro[iCell];
	par.e = par.E-0.5*(par.u*par.u+par.v*par.v);
	getMaterial(iCell).URS(par, 3);
}

void FVM_TVD::procExchangeGrads()
{
    exchange(gradR);
    exchange(gradP);
    exchange(gradU);
    exchange(gradV);
    exchange(gradT);
}


void FVM_TVD::procExchangeFields()
{
    exchange(ro);
    exchange(ru);
    exchange(rv);
    exchange(re);
}



