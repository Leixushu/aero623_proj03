/* How to use gsl library
https://blog.csdn.net/raodotcong/article/details/8998379
https://zoujiemeng.github.io/post/2018/03/gsl-vs-linux.html
https://www.cnblogs.com/flipped/p/9314461.html

copy .dll files from ../../gsl-cmake-x64/bin/Release/ to working directory

properties -> VC++ Directories -> Include Directories
Add ../gsl-my-files/header/

properties -> VC++ Directories -> Library Directories
Add ../gsl-my-files/lib/

properties -> Linker -> Input -> Additional Dependencies
Add gsl.lib; gslcblas.lib;
*/

#include "gsl/gsl_spmatrix.h"

void create_adjacent_cell(int t1, int e1, int t2, int e2, MESH *mesh)
{
	ADJ_CELL new_info_1, new_info_2;
	new_info_1.loc_edge = e1;
	new_info_1.adj_edge = e2;
	new_info_1.adj_tri = t2;
	mesh->Elems[t1].adj_cell.push_back(new_info_1);
	new_info_2.loc_edge = e2;
	new_info_2.adj_edge = e1;
	new_info_2.adj_tri = t1;
	mesh->Elems[t2].adj_cell.push_back(new_info_2);
}

void check(EDGE_TEMP edge, gsl_spmatrix *conn_vertex, MESH *mesh)
{
	int n0, n1, t, e;
	n0 = edge.nodes[0];
	int len = _msize(edge.nodes)/sizeof(edge.nodes[0]);
	n1 = edge.nodes[len-1];
	t = edge.t;
	e = edge.e;
	int gloal_edge_no;
	if (gsl_spmatrix_get(conn_vertex, n0, n1) == 0)
	{
		gloal_edge_no = mesh->nEdge;
		EDGE new_edge;
		new_edge.nodes = new int[len];
		for (int i = 0; i < len; i++)
			new_edge.nodes[i] = edge.nodes[i];
		new_edge.t1 = t;
		new_edge.e1 = e;
		new_edge.t2 = -1; //will be updated when visit again
		new_edge.e2 = -1;
		mesh->Edges.push_back(new_edge);
		mesh->nEdge = mesh->nEdge + 1;
		gsl_spmatrix_set(conn_vertex, n0, n1, gloal_edge_no + 1);
		gsl_spmatrix_set(conn_vertex, n1, n0, gloal_edge_no + 1);
		mesh->Elems[t].edge[e] = gloal_edge_no;
		//the index in sparse matrix conn_vertex starts from 1
	}
	else
	{
		//if the edge already in the sparse matrix
		int t1, e1, t2, e2;
		t2 = t;
		e2 = e;
		//decide which edge already in the list
		gloal_edge_no = gsl_spmatrix_get(conn_vertex, n0, n1) - 1;
		mesh->Elems[t2].edge[e2] = gloal_edge_no;
		t1 = mesh->Edges[gloal_edge_no].t1;
		e1 = mesh->Edges[gloal_edge_no].e1;
		//store the new edge
		mesh->Edges[gloal_edge_no].t2 = t2;
		mesh->Edges[gloal_edge_no].e2 = e2;
		mesh->Edges[gloal_edge_no].type = "Interior";
		create_adjacent_cell(t1, e1, t2, e2, mesh);
	}
}

void find_adjecent_cell(MESH *mesh, INPUT_TEMP temp, gsl_spmatrix *conn_vertex)
{
	int nNode = mesh->nNode;
	int nElemTot = mesh->nElemTot;
	EDGE_TEMP **edge_temp = temp.edge_temp;
	int i, j;
	for (i = 0; i < nElemTot; i++)
		for (j = 0; j < 3; j++)
			check(edge_temp[i][j], conn_vertex, mesh);
}

/*-------------------------------------------------------------------------------*/
void find_mass_matrix(MESH *mesh)                 // p is the order of solution p = 0,1,2
{
	int i, j, k;
	int p = mesh->p;
	int np = (p + 1)*(p + 2) / 2;                 // # of Lagrange pts, np, i.e. # of unknowns per cell
	for (i = 0; i < mesh->nEGroup; i++)           // loop over all EGroups (linear and curved)
	{
		EGROUP EGroup = mesh->EGroup[i];          // for the i th EGroup
		int nElem = EGroup.nElem;
		int q = EGroup.Order;                     // define q, which is the geo order, q = 1 (linear) or q = 3 (curved)
		int nq = (q + 1)*(q + 2) / 2;             // 3 nodes for linear cell, 10 nodes per curved cell
		QUADRATURE_2D quad;
		select_2d_quad_pts(p, q, &quad);
		double **node_pos;                        // define local node_pos
		node_pos = new double*[nq];
		for (j = 0; j < nq; j++)
			node_pos[j] = new double[2];         // x and y coordinates
		for (j = 0; j < nElem; j++)              // loop over Elem within the i th elem group
		{
			int global_elem_no = EGroup.E[j];
			mesh->Elems[global_elem_no].M = gsl_matrix_alloc(np, np); // allocated the mass matrix
			for (k = 0; k < nq; k++)                                  // traverse all nq Lagrange nodes in elem E[j] and find their positions, which will be used when computing J
			{
				int global_node_no = mesh->Elems[global_elem_no].nodes[k];
				node_pos[k][0] = mesh->node_pos[global_node_no][0];
				node_pos[k][1] = mesh->node_pos[global_node_no][1];
			}
			// ---------------------------------------- DEBUG ----------------------------------------
			// cout << np << " "<< quad.n << endl;
			// ---------------------------------------- DEBUG ----------------------------------------
			gsl_matrix *Phi = gsl_matrix_alloc(quad.n, np); // quad.n row * np column
			gsl_matrix *Diag = gsl_matrix_alloc(quad.n, quad.n);
			gsl_matrix_set_identity(Diag);
			for (k = 0; k < quad.n; k++)                    // traverse all quad pts and construct matrix Phi and Diag
			{
				gsl_matrix *J = gsl_matrix_alloc(2, 2);
				find_Jacobian(q, quad.x[k], quad.y[k], node_pos, J);          // find the Jacobian matrix at the k th quad pt in elem E[j]
				double J_det = find_det(J);
				gsl_vector *psi = gsl_vector_alloc(np);
				find_2d_basis_fcn(p, quad.x[k], quad.y[k], psi); // compute the basis fcn at those quad pts
				/*// ---------------------------------------- DEBUG ----------------------------------------
			    for (int ib = 0;ib < psi->size;ib++)
				cout << psi->data[ib] << " ";
				cout << endl;
				// ---------------------------------------- DEBUG ----------------------------------------*/																		  
				gsl_matrix_set_row(Phi, k, psi);
				gsl_matrix_set(Diag, k, k, quad.w[k] * J_det);
				// free memory
				gsl_matrix_free(J);
				gsl_vector_free(psi);
			}
			gsl_matrix *Temp = gsl_matrix_alloc(np, quad.n); // np row * quad.m column
			gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Phi, Diag, 0.0, Temp);
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Temp, Phi, 0.0, mesh->Elems[global_elem_no].M);
			// free memory
			gsl_matrix_free(Phi);
			gsl_matrix_free(Diag);
			gsl_matrix_free(Temp);
		}
		// free memory
		for (j = 0; j < nq; j++)
			delete[] node_pos[j];
		delete[] node_pos;
		quad.clear();
	}
}

/*-------------------------------------------------------------------------------*/
void furnish_boundary(MESH *mesh, INPUT_TEMP temp, gsl_spmatrix *conn_vertex)
{
	int i, j, n0, n1, global_edge_no;
	string Title;
	int ***B_temp = temp.B_temp;
	for (i = 0;i < mesh->nBGroup;i++) // 4 boundary groups
	{
		Title = mesh->BGroup[i].Title;
		mesh->BGroup[i].B = new int[mesh->BGroup[i].nBFace]; // allocate nBFace integers in B
		for (j = 0;j < mesh->BGroup[i].nBFace;j++)
		{
			n0 = B_temp[i][j][0];
			int len = _msize(B_temp[i][j])/sizeof(B_temp[i][j][0]);
			n1 = B_temp[i][j][len-1];
			global_edge_no = gsl_spmatrix_get(conn_vertex, n0, n1) - 1;
			mesh->Edges[global_edge_no].type = Title;
			mesh->BGroup[i].B[j] = global_edge_no;
		}
	}
}

/*-------------------------------------------------------------------------------*/
void furnish_edge(MESH *mesh)
{
	int i, A, B;
	double xA, yA, xB, yB, vec_len;
	for (i = 0;i < mesh->nEdge;i++)
	{
		int n = _msize(mesh->Edges[i].nodes) / sizeof(mesh->Edges[i].nodes[0]); //n = nodes_per_edge
		// _msize is used for dynamic arrary (defined by "new"), sizeof for static arrary
		A = mesh->Edges[i].nodes[0]; // first node on the edge
		B = mesh->Edges[i].nodes[n-1]; // last node on the edge
		xA = mesh->node_pos[A][0];
		yA = mesh->node_pos[A][1];
		xB = mesh->node_pos[B][0];
		yB = mesh->node_pos[B][1];
		//find length
		vec_len = sqrt((xA - xB)*(xA - xB) + (yA - yB)*(yA - yB));
		mesh->Edges[i].length = vec_len;                           // not accurate for curved edges
		//find normal vector, based on the first and last node
		mesh->Edges[i].norm_vec[0] = (yB - yA) / vec_len;          // not accurate for curved edges
		mesh->Edges[i].norm_vec[1] = (xA - xB) / vec_len;
	}
}

/*
void furnish_cell(MESH *mesh)
{
	int i, j, k, global_edge_no, global_node_no;
	double len[3], S, vertex_pos[3][2];
	for (i = 0;i < mesh->nElemTot;i++)
	{
		for (j = 0;j < 3;j++)
		{
			global_edge_no = mesh->Elems[i].edge[j];
			len[j] = mesh->Edges[global_edge_no].length;
		}
		S = 0.5*(len[0] + len[1] + len[2]);
		mesh->Elems[i].A = sqrt(S*(S - len[0])*(S - len[1])*(S - len[2])); // not accurate for curved elem
	}
}
*/

void furnish_cell(MESH *mesh) // used the mass matrix to find the cell area
{
	for (int i = 0;i < mesh->nElemTot;i++)
	{
		gsl_matrix *M = mesh->Elems[i].M;
		double sum = 0;
		for (int j = 0;j < M->size1;j++)
			for (int k = 0;k < M->size2;k++)
				sum += gsl_matrix_get(M, j, k);
		mesh->Elems[i].A = sum;
	}
}

void find_local_node_pos(int global_elem_no, MESH *mesh, double **node_pos)
{
	int q = mesh->Elems[global_elem_no].q; // geometry order, for linear elems q = 1, for  curved elems q = 3
	int nq = (q + 1)*(q + 2) / 2;
	// pick the node_pos for current elem
	for (int j = 0; j < nq; j++)
	{
		int global_node_no = mesh->Elems[global_elem_no].nodes[j];
		node_pos[j][0] = mesh->node_pos[global_node_no][0];
		node_pos[j][1] = mesh->node_pos[global_node_no][1];
	}
}

void trans_quad_1d_to_2d(QUADRATURE_1D quad_1d, int e, double **quad_2d_x, int rev)
{
	for (int i = 0; i < quad_1d.n; i++)
	{
		double s;
		if (rev == 0) // counter - clockwise
			s = quad_1d.x[i];
		else // clockwise
			s = 1.0 - quad_1d.x[i];
		switch (e)
		{
		case 0:
			quad_2d_x[i][0] = 1.0 - s; // xi
			quad_2d_x[i][1] = s;       // eta
			break;
		case 1:
			quad_2d_x[i][0] = 0;
			quad_2d_x[i][1] = 1.0 - s;
			break;
		case 2:
			quad_2d_x[i][0] = s;
			quad_2d_x[i][1] = 0;
			break;
		}
	}
}

void find_J_edge_norm_vec(int e, int n, int q, double **quad_2d_x, double **node_pos, double *J_edge, double**norm_vec)
{
	double dxi_dsig, deta_dsig;
	switch (e)
	{
	case 0:
		dxi_dsig = -1;
		deta_dsig = 1;
		break;
	case 1:
		dxi_dsig = 0;
		deta_dsig = -1;
		break;
	case 2:
		dxi_dsig = 1;
		deta_dsig = 0;
		break;
	}
	for (int i = 0;i < n;i++)
	{
		double xi = quad_2d_x[i][0];
		double eta = quad_2d_x[i][1];
		gsl_matrix *J = gsl_matrix_alloc(2, 2);
		find_Jacobian(q, xi, eta, node_pos, J);
		double tan_vec[2];
		tan_vec[0] = gsl_matrix_get(J, 0, 0) * dxi_dsig + gsl_matrix_get(J, 0, 1) * deta_dsig;
		tan_vec[1] = gsl_matrix_get(J, 1, 0) * dxi_dsig + gsl_matrix_get(J, 1, 1) * deta_dsig;
		J_edge[i] = sqrt(tan_vec[0] * tan_vec[0] + tan_vec[1] * tan_vec[1]);
		norm_vec[i][0] = tan_vec[1] / J_edge[i];
		norm_vec[i][1] = -tan_vec[0] / J_edge[i];
		// free memory
		gsl_matrix_free(J);
	}
}

void fix_edge_length(MESH *mesh)
{
	for (int i = 0;i < mesh->nEdge;i++)
	{
		if (mesh->Edges[i].type == "Bottom")
		{
			EDGE edge = mesh->Edges[i];
			int t = edge.t1;
			double **node_pos; // node_pos for t elem, delete after use
			int q = mesh->Elems[t].q;
			int nq = (q + 1)*(q + 2) / 2;
			node_pos = new double*[nq];
			for (int j = 0; j < nq; j++)
				node_pos[j] = new double[2];
			find_local_node_pos(t, mesh, node_pos);
			/* ---------------------------------------------------------------------- */
			int e = edge.e1;
			QUADRATURE_1D quad_1d;
			int p = mesh->p;
			select_1d_quad_pts(p, q, &quad_1d);
			// mapping the 1d quad pts position to 2d position (xi,eta) based on the local edge number
			double **quad_2d_x;
			quad_2d_x = new double*[quad_1d.n]; // remember to delete after use
			for (int i = 0; i < quad_1d.n; i++)
				quad_2d_x[i] = new double[2];
			trans_quad_1d_to_2d(quad_1d, e, quad_2d_x, 0);
			/* ---------------------------------------------------------------------- */
			double *J_edge;
			J_edge = new double[quad_1d.n];
			double **norm_vec; // not used here
			norm_vec = new double*[quad_1d.n];
			for(int i =0;i < quad_1d.n; i++)
				norm_vec[i] = new double[2]; // remember to delete after use
			find_J_edge_norm_vec(e, quad_1d.n, q, quad_2d_x, node_pos, J_edge, norm_vec);
			/* ---------------------------------------------------------------------- */
			double delta_l = 0;
			for (int i = 0; i < quad_1d.n; i++)
				delta_l += J_edge[i] * quad_1d.w[i];
			mesh->Edges[i].length = delta_l;
			// free memory
			for (int j = 0; j < nq; j++)
				delete[] node_pos[j];
			delete[] node_pos;
			quad_1d.clear();
			for (int i = 0; i < quad_1d.n; i++)
			{
				delete[] quad_2d_x[i];
				delete[] norm_vec[i];
			}	
			delete[] J_edge;
		}
	}
}

void pre_calculation(MESH *mesh, INPUT_TEMP temp)
{
	gsl_spmatrix *conn_vertex = gsl_spmatrix_alloc(mesh->nNode, mesh->nNode);
	find_adjecent_cell(mesh, temp, conn_vertex);
	find_mass_matrix(mesh);
	furnish_boundary(mesh, temp, conn_vertex);
	furnish_edge(mesh);
	furnish_cell(mesh);
	// fix_edge_length(mesh);
}