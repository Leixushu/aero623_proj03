#include "stdafx.h"
#include "read_data.cpp"
#include "find_basis_fcns_and_grad.cpp"
#include "basic_matrix_operation.cpp"
#include "select_quad_pts.cpp"
#include "pre_calculation.cpp"
#include "initialization.cpp"
#include "fluxes.cpp"
#include "time_iteration.cpp"

void free_stream_test(string case_name, int p)
{
	// task_2_a
	cout << "task_2_a" << endl;
	MESH mesh;
	INPUT_TEMP temp;
	read_data("mesh/" + case_name + ".gri", &mesh, &temp);
	mesh.p = p;
	pre_calculation(&mesh, temp);
	int np = (p + 1)*(p + 2) / 2;
	gsl_matrix *state = gsl_matrix_alloc(mesh.nElemTot * np, 4);
	initialization(&mesh, state);
	int nTime_step = 1000;
	string fname = "data/task_2_ab/free_stream_test";
	double *R_L_inf;
	R_L_inf = new double[nTime_step];
	time_iteration(&mesh, state, nTime_step, fname, R_L_inf);
	save_history(R_L_inf, nTime_step, fname + "_history" + ".txt");
}

void FEM(string case_name, int p)
{
	MESH mesh;
	INPUT_TEMP temp;
	read_data("mesh/" + case_name + ".gri", &mesh, &temp);
	mesh.p = p;
	pre_calculation(&mesh, temp);
	int np = (p + 1)*(p + 2) / 2;
	gsl_matrix *state = gsl_matrix_alloc(mesh.nElemTot * np, 4);
	initialization(&mesh, state);
	long long int nTime_step = 1000000;
	string fname = "data/task_2_c/" + case_name + "_p" + to_string(p) + "/" + case_name;
	double *R_L_inf;
	R_L_inf = new double[nTime_step];
	time_iteration(&mesh, state, nTime_step, fname, R_L_inf);
}

void mapping_from_ref_to_global(double **node_pos, double **sample_ref, double **sample_global, int q, int ns)
{
	int nq = (q + 1)*(q + 2) / 2;
	// store the node_pos of elem into a (2 by nq) gsl_matrix
	gsl_matrix *NODE_POS = gsl_matrix_alloc(2, nq);
	for (int i = 0; i < nq; i++)
	{
		gsl_matrix_set(NODE_POS, 0, i, node_pos[i][0]);
		gsl_matrix_set(NODE_POS, 1, i, node_pos[i][1]);
	}
	// loop over all sample points
	for (int i = 0;i < ns;i++)
	{
		double xi = sample_ref[i][0];
		double eta = sample_ref[i][1];
		gsl_vector *psi = gsl_vector_alloc(nq);
		find_2d_basis_fcn(q, xi, eta, psi); // psi: nq*1 array
		gsl_vector_view SAMPLE_GLOBAL = gsl_vector_view_array(sample_global[i], 2);
		// matrix - vector multiplcation
		gsl_blas_dgemv(CblasNoTrans, 1.0, NODE_POS, psi, 0.0, &SAMPLE_GLOBAL.vector);
		gsl_vector_free(psi);
	}
}

void post_processing_mach(string case_name, int p)
{
	// load mesh
	MESH mesh;
	INPUT_TEMP temp;
	read_data("mesh/" + case_name + ".gri", &mesh, &temp);
	mesh.p = p;
	pre_calculation(&mesh, temp);
	int np = (p + 1)*(p + 2) / 2;
	
	// load state
	gsl_matrix *state = gsl_matrix_alloc(mesh.nElemTot * np, 4);
	string fname = "data/task_2_c/" + case_name + "_p" + to_string(p) + "/" + case_name + "_state_final.txt";
	load_state(&mesh, fname, state);
	
	// load single right elem
	MESH mesh_s;
	INPUT_TEMP temp_s;
	read_data("mesh/single_right_tri5.gri", &mesh_s, &temp_s);
	mesh.p = 0;
	pre_calculation(&mesh_s, temp_s);
	
	// find position of sample pts in ref space 
	int ns = mesh_s.nNode; // number of sampling pts
	double **sample_ref;
	sample_ref = new double*[ns]; // remember to delete after use
	for (int i = 0; i < ns; i++)
	{
		sample_ref[i] = new double[2];
		sample_ref[i][0] = mesh_s.node_pos[i][0];
		sample_ref[i][1] = mesh_s.node_pos[i][1];
	}

	// pre-define position of sample pts in global space 
	double **sample_global;
	sample_global = new double*[ns]; // remember to delete after use
	for (int i = 0; i < ns; i++)
		sample_global[i] = new double[2];
	
	// first output trianglulation result, which is independent of elems
	ofstream output;
	fname = "data/task_2_c/" + case_name + "_p" + to_string(p) + "/" + case_name + "_triangulation.txt";
	output.open(fname);
	output << mesh_s.nElemTot << " " << 3 << endl;
	for (int i = 0;i < mesh_s.nElemTot;i++)
	{
		for (int j = 0;j < 3;j++)
			output << mesh_s.Elems[i].nodes[j] << " ";
		output << endl;
	}
	output.close();

	// pre-define the states @ sample pts
	double **U_at_sample_pts;
	U_at_sample_pts = new double*[ns]; // remember to delete after use
	for (int i = 0; i < ns; i++)
		U_at_sample_pts[i] = new double[4];

	// output data for Mach plot
	fname = "data/task_2_c/" + case_name + "_p" + to_string(p) + "/" + case_name + "_Mach_number.txt";
	output.open(fname);

	output << mesh.nElemTot << " " << ns << endl;
	// output x,y coordinates of sample pts in global space and their Mach number
	for (int i = 0; i < mesh.nElemTot; i++) // traverse all elems
	{
		int q = mesh.Elems[i].q;
		int nq = (q + 1)*(q + 2) / 2;
		// first find the position of node pts in global space for the ith elem
		double **node_pos; // delete after use
		node_pos = new double*[nq]; // nq*2 array, node_pos of the element
		for (int j = 0; j < nq; j++)
			node_pos[j] = new double[2];
		find_local_node_pos(i, &mesh, node_pos);

		// mapping the position of sample pts from the ref space into global space
		mapping_from_ref_to_global(node_pos, sample_ref, sample_global, q, ns);	

        // pick the (address of) submatrix of "state"
		gsl_matrix_view STATE_ELEM = gsl_matrix_submatrix(state, i * np, 0, np, 4);
		// find states @ sample pts
		find_U_at_quad(&STATE_ELEM.matrix, sample_ref, p, ns, U_at_sample_pts);
		// compute Mach number @ those sample pts
		double U[4];
		double *M;
		M = new double[ns]; // remember to delete
		for (int j = 0;j < ns;j++)
		{
			for (int k = 0; k < 4; k++)
				U[k] = U_at_sample_pts[j][k];
			double rho = U[0];
			double u = U[1] / U[0];
			double v = U[2] / U[0];
			double E = U[3] / U[0];
			double p = (gamma - 1.0)*(rho*E - 0.5*rho*(u*u + v*v));
			double a = sqrt(gamma*p / rho);
			M[j] = sqrt(u*u + v*v) / a;
		}
		// fout to file
		for (int j = 0; j < ns; j++)
			output << sample_global[j][0] << " " << sample_global[j][1] << " " << M[j] << endl;
		// free memory
		for (int j = 0; j < nq; j++)
			delete[] node_pos[j];
		delete[] node_pos;
		delete[] M;
	}
	output.close();

	// free memory
	for (int i = 0; i < ns; i++)
	{
		delete[] sample_ref[i];
		delete[] sample_global[i];
		delete[] U_at_sample_pts[i];
	}
	delete[] sample_ref;
	delete[] sample_global;
	delete[] U_at_sample_pts;
	gsl_matrix_free(state);
}

void post_processing_coef(string case_name, int p)
{
	// load mesh
	MESH mesh;
	INPUT_TEMP temp;
	read_data("mesh/" + case_name + ".gri", &mesh, &temp);
	mesh.p = p;
	int np = (p + 1)*(p + 2) / 2;
	pre_calculation(&mesh, temp);

	printf("%5d  ", mesh.nElemTot * np);

	// load state
	gsl_matrix *state = gsl_matrix_alloc(mesh.nElemTot * np, 4);
	string fname = "data/task_2_c/" + case_name + "_p" + to_string(p) + "/" + case_name + "_state_final.txt";
	load_state(&mesh, fname, state);

	int q = 3; // for bottom edges, q = 3 always hold true
	int nq = (q + 1)*(q + 2) / 2;
	
	QUADRATURE_1D quad_1d;
	select_1d_quad_pts(p, q, &quad_1d); 

	// pre-define position of sample pts in global space 
	double **quad_global;
	quad_global = new double*[quad_1d.n]; // remember to delete after use
	for (int i = 0; i < quad_1d.n; i++)
		quad_global[i] = new double[2];


	fname = "data/task_2_c/" + case_name + "_p" + to_string(p) + "/" + case_name + "_cp_distribution.txt";
	ofstream output;
	output.open(fname);
	
	// loop over the bottom boundary and sum over "cd_Tot","cl_Tot"
	double cd_Tot = 0, cl_Tot = 0;
	for (int i = 0;i < mesh.nEdge;i++)
	{
		if (mesh.Edges[i].type == "Bottom")
		{
			int e = mesh.Edges[i].e1; // e is the local edge #
			int t = mesh.Edges[i].t1; // t is the adjacent elem #
			// figure out the state for adjacent element "t"
			gsl_matrix_view STATE_ELEM = gsl_matrix_submatrix(state, t * np, 0, np, 4);

			// figure out the node pos for adjacent element "e"
			double **node_pos; // node_pos for t1 elem, delete after use
			node_pos = new double*[nq];
			for (int j = 0; j < nq; j++)
				node_pos[j] = new double[2];
			find_local_node_pos(t, &mesh, node_pos);
		
			// mapping the 1d quad pts position to 2d position (xi,eta) based on the local edge number
			double **quad_2d_x;
			quad_2d_x = new double*[quad_1d.n]; // remember to delete after use
			for (int j = 0; j < quad_1d.n; j++)
				quad_2d_x[j] = new double[2];
			trans_quad_1d_to_2d(quad_1d, e, quad_2d_x, 0);

			// mapping the position of sample pts from the ref space into global space
			mapping_from_ref_to_global(node_pos, quad_2d_x, quad_global, q, quad_1d.n);

			// find the left state on the edge
			double **U_edge;
			U_edge = new double*[quad_1d.n]; // remember to delete after use
			for (int j = 0; j < quad_1d.n; j++)
				U_edge[j] = new double[4];
			find_U_at_quad(&STATE_ELEM.matrix, quad_2d_x, p, quad_1d.n, U_edge);  // find the "left" states on the edge
			
			// find the normal vector and J_edge at each quad pt, for curved edges they are NOT const
			double *J_edge;
			J_edge = new double[quad_1d.n];
			
			double **n;
			n = new double*[quad_1d.n];
			for (int j = 0;j < quad_1d.n; j++)
				n[j] = new double[2]; // remember to delete after use

			// when geo order q > 1, curved elem, always true for bottom elems
			find_J_edge_norm_vec(e, quad_1d.n, q, quad_2d_x, node_pos, J_edge, n);
			
			// loop over quad pts inside an elem
			double cd_Elem = 0, cl_Elem = 0;
			for (int j = 0; j < quad_1d.n; j++)
			{
				PROP prop;
				Euler_flux(U_edge[j], &prop);
				double v_b[2];
				
				for (int k = 0;k < 2;k++)
					v_b[k] = prop.v[k] - (prop.v[0] * n[j][0] + prop.v[1] * n[j][1]) * n[j][k];
				double p_b = (gamma - 1)*(prop.rho*prop.E - 0.5*prop.rho*(v_b[0] * v_b[0] + v_b[1] * v_b[1]));
				
				// double p_b = prop.p;

				cl_Elem += (p_b - p_inf) * n[j][1] * quad_1d.w[j] * J_edge[j];
				cd_Elem += (p_b - p_inf) * n[j][0] * quad_1d.w[j] * J_edge[j];
				// figure out quad pt position "x_i" and pressure coef @ that pt "cp_i"
				double x_i = quad_global[j][0];
				double cp_i = (p_b - p_inf) / (0.5*gamma*p_inf*M_inf*M_inf);
				output << x_i << " " << cp_i << endl;
			}
			cl_Tot += cl_Elem;
			cd_Tot += cd_Elem;

			// free memory
			for (int j = 0; j < nq; j++)
				delete[] node_pos[j];
			delete[] node_pos;

			for (int j = 0; j < quad_1d.n; j++)
			{
				delete[] quad_2d_x[j];
				delete[] U_edge[j];
				delete[] n[j];
			}
			delete[] quad_2d_x;
			delete[] U_edge;
			delete[] J_edge;
			delete[] n;
		}
	}
	output.close();
	
	cl_Tot = cl_Tot / (0.5 * gamma * p_inf * M_inf*M_inf*h);
	cd_Tot = cd_Tot / (0.5 * gamma * p_inf * M_inf*M_inf*h);
	printf("%.16f %.16f ", cl_Tot, cd_Tot);

	// free memory
	quad_1d.clear();
	for (int i = 0; i < quad_1d.n; i++)
		delete[] quad_global[i];
	delete[] quad_global;

	/*-------------------------------------------------------------------------------*/

	// loop over all elems and sum over Es
	double T_t = 1 + ((gamma - 1) / 2) * M_inf * M_inf;
	double p_t = pow(T_t, gamma / (gamma - 1));
	double rho_t = p_t / (R*T_t);
	double s_t = p_t / pow(rho_t, gamma);
	double sums_Tot = 0, sumA_Tot = 0;

	for (int i = 0;i < mesh.nElemTot;i++)
	{

		q = mesh.Elems[i].q;
		nq = (q + 1)*(q + 2) / 2;

		// figure state for current element i
		gsl_matrix_view STATE_ELEM = gsl_matrix_submatrix(state, i * np, 0, np, 4);

		// figure out the node pos for current element "i"
		double **node_pos; // node_pos for t1 elem, delete after use
		node_pos = new double*[nq];
		for (int j = 0; j < nq; j++)
			node_pos[j] = new double[2];
		find_local_node_pos(i, &mesh, node_pos);

		QUADRATURE_2D quad;
		select_2d_quad_pts(p, q, &quad);

		gsl_matrix *J = gsl_matrix_alloc(2, 2);
		// gsl_matrix *J_inv = gsl_matrix_alloc(2, 2);
		double J_det;
		
		// for linear elems (q=1), the Jacobian matrix is const
		find_Jacobian(q, 0.0, 0.0, node_pos, J);
		// find_inv(J, J_inv);
		J_det = find_det(J);
		
		
		// loop over quad pts inside an elem
		double sums_Elem = 0;
		for (int j = 0;j < quad.n;j++)
		{
			double xi = quad.x[j];
			double eta = quad.y[j];
			// for curved elems (q>1), we need to update J, J_inv and J_det @ each quad pt
			if (q > 1)
			{
				find_Jacobian(q, xi, eta, node_pos, J);
				// find_inv(J, J_inv);
				J_det = find_det(J);
			}
			// find states @ current quad point
			gsl_vector *psi = gsl_vector_alloc(np);
			find_2d_basis_fcn(p, xi, eta, psi); // psi: np*1 array
			double state_quad[4];
			gsl_vector_view STATE_QUAD = gsl_vector_view_array(state_quad, 4);
			gsl_blas_dgemv(CblasTrans, 1.0, &STATE_ELEM.matrix, psi, 0.0, &STATE_QUAD.vector);
			// find (Euler) flux @ current quad point
			PROP prop;
			Euler_flux(state_quad, &prop);
			double s = prop.p / pow(prop.rho, gamma);
			sums_Elem += (s / s_t - 1)*(s / s_t - 1) * quad.w[j] * J_det;
			// free memory
			gsl_vector_free(psi);
		}
		sums_Tot += sums_Elem;
		sumA_Tot += mesh.Elems[i].A;
		// free memory
		for (int j = 0; j < nq; j++)
			delete[] node_pos[j];
		delete[] node_pos;
		quad.clear();
		gsl_matrix_free(J);
		// gsl_matrix_free(J_inv);
	}
	double Es = sqrt(sums_Tot / sumA_Tot);
	printf("%.16f\n", Es);

	// free memory
	gsl_matrix_free(state);
	
}



int main()
{
	// free_stream_test("bump0_curved",2); // 2nd order solution
	/*for (int i = 0; i < 3; i++)
	{
		FEM("bump0_curved", i);
		FEM("bump1_curved", i);
		FEM("bump2_curved", i);
	}*/

	cout << "    bump    | p | DOF |        cl        |        cd        |        ES        |" << endl;
	for (int p = 0; p < 3; p++)
	{
		for (int i = 0; i < 3; i++)
		{
			string mesh_name = "bump" + to_string(i) + "_curved";
			printf("bump%d_curved  %d ", i, p);
			// post_processing_mach(mesh_name , p);
			post_processing_coef(mesh_name, p);
		}
	}
		
	system("pause");
    return 0;
}

