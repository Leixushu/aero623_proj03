void read_data(string fname, MESH *mesh, INPUT_TEMP *temp)
{
	ifstream input(fname);
	int i, j, k;
	int nNode, nElemTot, Dim;
	input >> nNode >> nElemTot >> Dim;
	double **node_pos;
	node_pos = new double *[nNode];
	/*-------------------------------------------------------------------------------*/
	for (i = 0; i < nNode; i++)
	{
		node_pos[i] = new double[2];
		input >> node_pos[i][0] >> node_pos[i][1];
	}
	/*-------------------------------------------------------------------------------*/
	int nBGroup;
	input >> nBGroup;
	BGROUP *BGroup;
	BGroup = new BGROUP[nBGroup];
	// Use a temporary list "B_temp" to store the boundary nodes before the establishment of the global edge list
	// "B_temp" is a 3D list, nBGroup * nBFace * nf
	int ***B_temp;
	B_temp = new int**[nBGroup];
	for (i = 0; i < nBGroup; i++)
	{
		input >> BGroup[i].nBFace >> BGroup[i].nf >> BGroup[i].Title;
		B_temp[i] = new int*[BGroup[i].nBFace];
		for (j = 0; j < BGroup[i].nBFace; j++)
		{
			B_temp[i][j] = new int[BGroup[i].nf]; // store the nf nodes per face
			for (k = 0; k < BGroup[i].nf; k++)
			{
				input >> B_temp[i][j][k];
				B_temp[i][j][k] -= 1;
			}
		}
	}
	/*-------------------------------------------------------------------------------*/
	int curTot = 0;
	int nEGroup = 0;
	vector<EGROUP> EGroup;
	CELL *E_global;
	E_global = new CELL[nElemTot];
	// Use a temporary list "edge_temp" to store the edges before the establishment of the global edge list
	EDGE_TEMP **edge_temp;
	edge_temp = new EDGE_TEMP*[nElemTot];
	int q, nq;
	while (curTot < nElemTot)
	{
		EGROUP new_EGroup;
		input >> new_EGroup.nElem >> new_EGroup.Order >> new_EGroup.Basis;
		new_EGroup.E = new int[new_EGroup.nElem];
		q = new_EGroup.Order;
		nq = int(0.5*(q + 1)*(q + 2));		
		for (i = 0; i < new_EGroup.nElem; i++)
		{
			E_global[curTot].q = q;
			E_global[curTot].nodes = new int[nq];
			for (j = 0; j < nq; j++)
			{
				input >> E_global[curTot].nodes[j];
				E_global[curTot].nodes[j] -= 1;
			}
			edge_temp[curTot] = new EDGE_TEMP[3];
			for (j = 0; j < 3; j++)
			{
				edge_temp[curTot][j].nodes = new int[q+1];
				edge_temp[curTot][j].t = curTot;
				edge_temp[curTot][j].e = j;
			}
			if (q == 1)
			{
				edge_temp[curTot][0].nodes[0] = E_global[curTot].nodes[1];
				edge_temp[curTot][0].nodes[1] = E_global[curTot].nodes[2];
				edge_temp[curTot][1].nodes[0] = E_global[curTot].nodes[2];
				edge_temp[curTot][1].nodes[1] = E_global[curTot].nodes[0];
				edge_temp[curTot][2].nodes[0] = E_global[curTot].nodes[0];
				edge_temp[curTot][2].nodes[1] = E_global[curTot].nodes[1];
			}
			if (q == 3)
			{
				edge_temp[curTot][0].nodes[0] = E_global[curTot].nodes[3];
				edge_temp[curTot][0].nodes[1] = E_global[curTot].nodes[6];
				edge_temp[curTot][0].nodes[2] = E_global[curTot].nodes[8];
				edge_temp[curTot][0].nodes[3] = E_global[curTot].nodes[9];
				edge_temp[curTot][1].nodes[0] = E_global[curTot].nodes[9];
				edge_temp[curTot][1].nodes[1] = E_global[curTot].nodes[7];
				edge_temp[curTot][1].nodes[2] = E_global[curTot].nodes[4];
				edge_temp[curTot][1].nodes[3] = E_global[curTot].nodes[0];
				edge_temp[curTot][2].nodes[0] = E_global[curTot].nodes[0];
				edge_temp[curTot][2].nodes[1] = E_global[curTot].nodes[1];
				edge_temp[curTot][2].nodes[2] = E_global[curTot].nodes[2];
				edge_temp[curTot][2].nodes[3] = E_global[curTot].nodes[3];
			}
			new_EGroup.E[i] = curTot;
			curTot += 1;
		}
		EGroup.push_back(new_EGroup);
		nEGroup += 1;
	}
	input.close();
	mesh->nNode = nNode;
	mesh->node_pos = node_pos;
	mesh->nElemTot = nElemTot;
	mesh->nEGroup = nEGroup;
	mesh->EGroup = EGroup;
	mesh->nBGroup = nBGroup;
	mesh->BGroup = BGroup;
	mesh->Elems = E_global;
	mesh->nEdge = 0;
	// mesh->Edges = no edges established currently;
	temp->B_temp = B_temp;
	temp->edge_temp = edge_temp;
}