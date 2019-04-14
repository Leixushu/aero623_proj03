void interior_contribution(gsl_matrix*state, gsl_matrix*R, double **node_pos, int p, int q) // note state and R are np*4 matirx for single elem
{
	QUADRATURE_2D quad;
	select_2d_quad_pts(p, q, &quad);
	int np = (p + 1)*(p + 2) / 2;
	gsl_matrix *J = gsl_matrix_alloc(2, 2);
	gsl_matrix *J_inv = gsl_matrix_alloc(2, 2);
	double J_det;
	// for linear elems (q=1), the Jacobian matrix is const
	find_Jacobian(q, 0.0, 0.0, node_pos, J);
	find_inv(J, J_inv);
	J_det = find_det(J);
	// loop over all quad points
	gsl_matrix *sum = gsl_matrix_alloc(np, 4);
	gsl_matrix_set_zero(sum);
	for (int i = 0;i < quad.n;i++)
	{
		double xi = quad.x[i];
		double eta = quad.y[i];
		// for curved elems (q>1), we need to update J, J_inv and J_det @ each quad pt
		if (q > 1)
		{
			find_Jacobian(q, xi, eta, node_pos, J);
			find_inv(J, J_inv);
			J_det = find_det(J);
		}
		// find states @ current quad point
		gsl_vector *psi = gsl_vector_alloc(np);
		find_2d_basis_fcn(p, xi, eta, psi); // psi: np*1 array
		double state_quad[4];
		gsl_vector_view STATE_QUAD = gsl_vector_view_array(state_quad, 4);
		gsl_blas_dgemv(CblasTrans, 1.0, state, psi, 0.0, &STATE_QUAD.vector);
		// find (Euler) flux @ current quad point
		PROP prop;
		Euler_flux(state_quad, &prop);
		gsl_matrix *F = gsl_matrix_alloc(2, 4);
		for (int j = 0; j < 2; j++)
			for (int k = 0; k < 4; k++)
				gsl_matrix_set(F, j, k, prop.F[j][k]);
		gsl_matrix *grad_psi = gsl_matrix_alloc(np, 2);
		find_2d_grad_basis_fcn(p, xi, eta, grad_psi); // grad_psi: np*2 matrix
		gsl_matrix *temp1 = gsl_matrix_alloc(np, 2);               // np row * 2 column
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, grad_psi, J_inv, 0.0, temp1);
		gsl_matrix *temp2 = gsl_matrix_alloc(np, 4);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -J_det*quad.w[i], temp1, F, 0.0, temp2);
		gsl_matrix_add(sum, temp2);
		// free memory
		gsl_vector_free(psi);
		gsl_matrix_free(F);
		gsl_matrix_free(grad_psi);
		gsl_matrix_free(temp1);
		gsl_matrix_free(temp2);
	}
	gsl_matrix_add(R, sum);
	// free memory
	quad.clear();
	gsl_matrix_free(J);
	gsl_matrix_free(J_inv);
	gsl_matrix_free(sum);
}

void find_U_at_quad(gsl_matrix*state, double**quad_2d_x, int p, int quad_n, double **U_edge)
{
	for (int i = 0; i < quad_n; i++)
	{
		double xi = quad_2d_x[i][0];
		double eta = quad_2d_x[i][1];
		int np = (p + 1)*(p + 2) / 2;
		gsl_vector *psi = gsl_vector_alloc(np);
		find_2d_basis_fcn(p, xi, eta, psi); // psi: np*1 array
		double state_quad[4];
		gsl_vector_view STATE_QUAD = gsl_vector_view_array(state_quad, 4);
		gsl_blas_dgemv(CblasTrans, 1.0, state, psi, 0.0, &STATE_QUAD.vector);
		for (int j = 0; j < 4; j++) // or we may replace state_quad with U_edge[i]
			U_edge[i][j] = state_quad[j];
		gsl_vector_free(psi);
	}
}

void edge_contribution_interior(EDGE *edge, gsl_matrix *U1, gsl_matrix *U2, gsl_matrix *R1, gsl_matrix *R2, int p)
{
	int q = 1; // note that all interior edges are considered as linear, even if the elem is q=3 on both side
	int e1 = edge->e1;
	int e2 = edge->e2;
	QUADRATURE_1D quad_1d;
	select_1d_quad_pts(p, q, &quad_1d);
	// mapping the 1d quad pts position to 2d position (xi,eta) based on the local edge number
	double **quad_2d_x1, **quad_2d_x2;
	quad_2d_x1 = new double*[quad_1d.n]; // remember to delete after use
	quad_2d_x2 = new double*[quad_1d.n];
	for (int i = 0; i < quad_1d.n; i++)
	{
		quad_2d_x1[i] = new double[2];
		quad_2d_x2[i] = new double[2];
	}
	trans_quad_1d_to_2d(quad_1d, e1, quad_2d_x1, 0);
	trans_quad_1d_to_2d(quad_1d, e2, quad_2d_x2, 1); 
	// find the left and right state on the edge
	double **U_edge1, **U_edge2;
	U_edge1 = new double*[quad_1d.n]; // remember to delete after use
	U_edge2 = new double*[quad_1d.n];
	for (int i = 0; i < quad_1d.n; i++)
	{
		U_edge1[i] = new double[4];
		U_edge2[i] = new double[4];
	}
	find_U_at_quad(U1, quad_2d_x1, p, quad_1d.n, U_edge1);  // find the "left" states on the edge
	find_U_at_quad(U2, quad_2d_x2, p, quad_1d.n, U_edge2);  // find the "right" states on the edge
	double J_edge = edge->length;
	double norm_vec[] = { edge->norm_vec[0], edge->norm_vec[1] };
	// find Roe flux @ each quad pt
	double **F;
	F = new double*[quad_1d.n]; // remember to delete after use
	for (int i = 0; i < quad_1d.n; i++)
	{
		F[i] = new double[4];
		double s;
		Roe_flux(U_edge1[i], U_edge2[i], norm_vec, F[i], &s);
		if (s > edge->s)
			edge->s = s;
	}
	int np = (p + 1)*(p + 2) / 2;
	gsl_matrix *sum1 = gsl_matrix_alloc(np, 4);
	gsl_matrix_set_zero(sum1);
	gsl_matrix *sum2 = gsl_matrix_alloc(np, 4);
	gsl_matrix_set_zero(sum2);
	// update residual to R1 R2
	for (int i = 0; i < quad_1d.n; i++)
	{
		double xi1 = quad_2d_x1[i][0];
		double eta1 = quad_2d_x1[i][1];
		gsl_vector *psi1 = gsl_vector_alloc(np);
		find_2d_basis_fcn(p, xi1, eta1, psi1); // psi: np*1 vector
		gsl_matrix *PSI1 = gsl_matrix_alloc(np, 1);
		gsl_matrix_set_col(PSI1, 0, psi1);
		double xi2 = quad_2d_x2[i][0];
		double eta2 = quad_2d_x2[i][1];
		gsl_vector *psi2 = gsl_vector_alloc(np);
		find_2d_basis_fcn(p, xi2, eta2, psi2); // psi: np*1 vector
		gsl_matrix *PSI2 = gsl_matrix_alloc(np, 1);
		gsl_matrix_set_col(PSI2, 0, psi2);
		gsl_matrix_view FLUX = gsl_matrix_view_array(F[i], 1, 4);
		gsl_matrix *temp1 = gsl_matrix_alloc(np, 4);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, J_edge * quad_1d.w[i], PSI1, &FLUX.matrix, 0.0, temp1);
		gsl_matrix_add(sum1, temp1);
		gsl_matrix *temp2 = gsl_matrix_alloc(np, 4);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -J_edge * quad_1d.w[i], PSI2, &FLUX.matrix, 0.0, temp2);
		gsl_matrix_add(sum2, temp2);
		// free memory
		gsl_vector_free(psi1);
		gsl_matrix_free(PSI1);
		gsl_vector_free(psi2);
		gsl_matrix_free(PSI2);
		gsl_matrix_free(temp1);
		gsl_matrix_free(temp2);
	}
	gsl_matrix_add(R1, sum1);
	gsl_matrix_add(R2, sum2);
	// free memory
	for (int i = 0; i < quad_1d.n; i++)
	{	
		delete[] quad_2d_x1[i];
		delete[] quad_2d_x2[i];
		delete[] U_edge1[i];
		delete[] U_edge2[i];
		delete[] F[i];
	}
	delete[] quad_2d_x1;
	delete[] quad_2d_x2;
	delete[] U_edge1;
	delete[] U_edge2;
	delete[] F;
	quad_1d.clear();
	gsl_matrix_free(sum1);
	gsl_matrix_free(sum2);
}

void edge_contribution_boundary(EDGE *edge, gsl_matrix *U, gsl_matrix *R, double **node_pos, int p, int q)
{
	int e = edge->e1;
	QUADRATURE_1D quad_1d;
	select_1d_quad_pts(p, q, &quad_1d);
	// mapping the 1d quad pts position to 2d position (xi,eta) based on the local edge number
	double **quad_2d_x;
	quad_2d_x = new double*[quad_1d.n]; // remember to delete after use
	for (int i = 0; i < quad_1d.n; i++)
		quad_2d_x[i] = new double[2];
	trans_quad_1d_to_2d(quad_1d, e, quad_2d_x, 0);
	// find the left state on the edge
	double **U_edge;
	U_edge = new double*[quad_1d.n]; // remember to delete after use
	for (int i = 0; i < quad_1d.n; i++)
		U_edge[i] = new double[4];
	find_U_at_quad(U, quad_2d_x, p, quad_1d.n, U_edge);  // find the "left" states on the edge
	// find the normal vector and J_edge at each quad pt, for curved edges they are NOT const
	double *J_edge;
	J_edge = new double[quad_1d.n];
	double **norm_vec; // not used here
	norm_vec = new double*[quad_1d.n];
	for (int i = 0;i < quad_1d.n; i++)
		norm_vec[i] = new double[2]; // remember to delete after use
	// cout << "q = " << q << endl;
	if (q > 1) // when geo order q > 1, curved elem
		find_J_edge_norm_vec(e, quad_1d.n, q, quad_2d_x, node_pos, J_edge, norm_vec);
	else
	{
		for (int i = 0;i < quad_1d.n;i++)
		{
			J_edge[i] = edge->length;
			norm_vec[i][0] = edge->norm_vec[0];
			norm_vec[i][1] = edge->norm_vec[1];
		}
	}
	// find boundary flux @ each quad pt
	double **F;
	F = new double*[quad_1d.n]; // remember to delete after use
	for (int i = 0; i < quad_1d.n; i++)
	{
		F[i] = new double[4];
		double s;
		/* // ----- free stream test start -----
		double U_inf[4];
		define_U(0.5, 0, gamma, U_inf);
		Roe_flux(U_edge[i], U_inf, norm_vec[i], F[i], &s);
		// ----- free stream test end ----- */

		// ----- boundary condition start -----
		if (edge->type == "Left") // inflow
		{
			double T_t = 1 + ((gamma - 1) / 2) * M_inf * M_inf;
			double p_t = pow(T_t, gamma / (gamma - 1));
			inlet_flux(U_edge[i], norm_vec[i], T_t, p_t, alpha, F[i], &s);
		}
		if (edge->type == "Top" || edge->type == "Bottom")
			wall_flux(U_edge[i], norm_vec[i], F[i], &s);
		if (edge->type == "Right")
			outflow_flux(U_edge[i], norm_vec[i], p_inf, F[i], &s);
		// ----- boundary condition end ----- */
		if (s > edge->s)
			edge->s = s;
	}
	// update residual to R1
	int np = (p + 1)*(p + 2) / 2;
	gsl_matrix *sum = gsl_matrix_alloc(np, 4);
	gsl_matrix_set_zero(sum);
	for (int i = 0; i < quad_1d.n; i++)
	{
		double xi = quad_2d_x[i][0];
		double eta = quad_2d_x[i][1];
		gsl_vector *psi = gsl_vector_alloc(np);
		find_2d_basis_fcn(p, xi, eta, psi); // psi: np*1 vector
		gsl_matrix *PSI = gsl_matrix_alloc(np, 1);
		gsl_matrix_set_col(PSI, 0, psi);
		gsl_matrix_view FLUX = gsl_matrix_view_array(F[i], 1, 4);
		gsl_matrix *temp = gsl_matrix_alloc(np, 4);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, J_edge[i] * quad_1d.w[i], PSI, &FLUX.matrix, 0.0, temp);
		gsl_matrix_add(sum, temp);
		// free memory
		gsl_vector_free(psi);
		gsl_matrix_free(PSI);
		gsl_matrix_free(temp);
	}
	gsl_matrix_add(R, sum);
	// free memory
	quad_1d.clear();
	for (int i = 0; i < quad_1d.n; i++)
	{
		delete[] quad_2d_x[i];
		delete[] U_edge[i];
		delete[] norm_vec[i];
		delete[] F[i];
	}
	delete[] quad_2d_x;
	delete[] U_edge;
	delete[] J_edge;
	delete[] norm_vec;
	delete[] F;
	gsl_matrix_free(sum);
}

void refresh_R(MESH *mesh, gsl_matrix *state, gsl_matrix *R)
{
	// initialize R to zero
	gsl_matrix_set_zero(R);
	int i, j, k;
	int p = mesh->p;  // solution order
	int np = (p + 1)*(p + 2) / 2;
	// loop over elem interiors
	for (i = 0; i < mesh->nElemTot; i++)
	{
		int q = mesh->Elems[i].q; // geometry order, for linear elems q = 1, for  curved elems q = 3
		int nq = (q + 1)*(q + 2) / 2;
		gsl_matrix_view STATE_ELEM = gsl_matrix_submatrix(state, i * np, 0, np, 4); // pick the (address of) submatrix of "state"
		gsl_matrix_view R_ELEM = gsl_matrix_submatrix(R, i * np, 0, np, 4);         // pick the (address of) submatrix of residual vector "R"
		double **node_pos; // define node_pos for current elem, delete after use
		node_pos = new double*[nq];
		for (int j = 0; j < nq; j++)
			node_pos[j] = new double[2];
		find_local_node_pos(i, mesh, node_pos);	
		// update the residual of current elem, meanwhile obtain s_max
		interior_contribution(&STATE_ELEM.matrix, &R_ELEM.matrix, node_pos, p, mesh->Elems[i].q);
		// delete node_pos for current elem
		for (int j = 0;j < nq;j++)
			delete[] node_pos[j];
		delete[] node_pos;
	}
	// loop over edges
	for (i = 0; i < mesh->nEdge; i++)
	{
		EDGE edge = mesh->Edges[i];
		if (edge.type == "Interior") // for interior edges
		{
			// For the left elem t1
			int t1 = edge.t1;
			gsl_matrix_view STATE_ELEM1 = gsl_matrix_submatrix(state, t1 * np, 0, np, 4); // pick the (address of) submatrix of "state"
			gsl_matrix_view R_ELEM1 = gsl_matrix_submatrix(R, t1 * np, 0, np, 4);
			int q1 = mesh->Elems[t1].q;
			// For the right elem t2
			int t2 = edge.t2;
			gsl_matrix_view STATE_ELEM2 = gsl_matrix_submatrix(state, t2 * np, 0, np, 4); // pick the (address of) submatrix of "state"
			gsl_matrix_view R_ELEM2 = gsl_matrix_submatrix(R, t2 * np, 0, np, 4);
			// find the node_pos for right elem
			int q2 = mesh->Elems[t2].q;
			edge_contribution_interior(&(mesh->Edges[i]), &STATE_ELEM1.matrix, &STATE_ELEM2.matrix, &R_ELEM1.matrix, &R_ELEM2.matrix, p);
		}
		else
		{
			int t = edge.t1;
			gsl_matrix_view STATE_ELEM = gsl_matrix_submatrix(state, t * np, 0, np, 4); // pick the (address of) submatrix of "state"
			gsl_matrix_view R_ELEM = gsl_matrix_submatrix(R, t * np, 0, np, 4);
			// pick the node_pos for current elem
			double **node_pos; // node_pos for t1 elem, delete after use
			int q = mesh->Elems[t].q;
			int nq = (q + 1)*(q + 2) / 2;
			node_pos = new double*[nq];
			for (j = 0; j < nq; j++)
				node_pos[j] = new double[2];
			find_local_node_pos(t, mesh, node_pos);
			edge_contribution_boundary(&(mesh->Edges[i]), &STATE_ELEM.matrix, &R_ELEM.matrix, node_pos, p, q);
			// delete node_pos for current elem
			for (j = 0; j < nq; j++)
				delete[] node_pos[j];
			delete[] node_pos;
		}
	}
}

void compute_time_step(MESH *mesh, double CFL)
{
	int i, j, global_edge_no;
	double SUM, s, delta_l;
	for (i = 0; i < mesh->nElemTot; i++)
	{
		SUM = 0;
		for (j = 0; j < 3; j++)
		{
			global_edge_no = mesh->Elems[i].edge[j];
			s = mesh->Edges[global_edge_no].s;
			delta_l = mesh->Edges[global_edge_no].length;
			SUM += s*delta_l;
		}
		mesh->Elems[i].dt = 2 * CFL * mesh->Elems[i].A / SUM;
		// the edge length for curved edges needs to be updated (if using the max wave speed)
	}
}

double RK4_intermediate_steps(MESH *mesh, gsl_matrix *state_in, gsl_matrix *state_out, gsl_matrix *R, bool timestep_update, double timestep_coef)
{
	int i, j, k;
	int p = mesh->p;
	int np = (p + 1)*(p + 2) / 2;
	refresh_R(mesh, state_in, R);
	if (timestep_update == 1)
		compute_time_step(mesh, CFL_ref/(p + 1.0));
	for (i = 0; i < mesh->nElemTot; i++)
	{
		gsl_matrix *M_inv = gsl_matrix_alloc(np, np);
		find_inv(mesh->Elems[i].M, M_inv);                                                // find inv of mass matrix for each elem
		gsl_matrix_view STATE_IN_ELEM = gsl_matrix_submatrix(state_in, i * np, 0, np, 4);       // pick the (address of) submatrix of "state"
		gsl_matrix_view STATE_OUT_ELEM = gsl_matrix_submatrix(state_out, i * np, 0, np, 4); // pick the (address of) submatrix of "state_FE"
		gsl_matrix_view R_ELEM = gsl_matrix_submatrix(R, i * np, 0, np, 4);               // pick the (address of) submatrix of residual vector "R"
		gsl_matrix *term_1 = gsl_matrix_alloc(np, 4);                                     // the 1^{st} term
		gsl_matrix_memcpy(term_1, &STATE_IN_ELEM.matrix);
		gsl_matrix *term_2 = gsl_matrix_alloc(np,4);                                      // the 2^{nd} term
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -timestep_coef*mesh->Elems[i].dt, M_inv, &R_ELEM.matrix, 0.0, term_2);
		gsl_matrix *sum = gsl_matrix_alloc(np, 4);
		gsl_matrix_set_zero(sum);                                                         // initialize sum as zero
		gsl_matrix_add(sum, term_1);                                                      // add two terms into sum
		gsl_matrix_add(sum, term_2);
		gsl_matrix_memcpy(&STATE_OUT_ELEM.matrix, sum);
		// free memory
		gsl_matrix_free(M_inv);
		gsl_matrix_free(term_1);
		gsl_matrix_free(term_2);
		gsl_matrix_free(sum);
	}
	double min_out, max_out, max_R;
	gsl_matrix_minmax(R, &min_out, &max_out);
	min_out = abs(min_out);
	max_out = abs(max_out);
	if (min_out > max_out)
		max_R = min_out;
	else
		max_R = max_out;
	return max_R;
}

void RK4_final_step(MESH *mesh, gsl_matrix *state, gsl_matrix *R0, gsl_matrix *R1, gsl_matrix *R2, gsl_matrix *R3)
{
	// cout << "RK4_final_step" << endl;
	int i, j;
	int p = mesh->p;
	int np = (p + 1)*(p + 2) / 2;
	for (i = 0; i < mesh->nElemTot; i++)
	{
		gsl_matrix *M_inv = gsl_matrix_alloc(np, np);
		find_inv(mesh->Elems[i].M, M_inv);                                                // find inv of mass matrix for each elem
		gsl_matrix_view STATE_ELEM = gsl_matrix_submatrix(state, i * np, 0, np, 4);       // find the submatrix of state

		gsl_matrix_view R0_ELEM = gsl_matrix_submatrix(R0, i * np, 0, np, 4);         // find the residual of elem i, np*4 matrix
		gsl_matrix_view R1_ELEM = gsl_matrix_submatrix(R1, i * np, 0, np, 4);
		gsl_matrix_view R2_ELEM = gsl_matrix_submatrix(R2, i * np, 0, np, 4);
		gsl_matrix_view R3_ELEM = gsl_matrix_submatrix(R3, i * np, 0, np, 4);
		
		gsl_matrix_scale(&R1_ELEM.matrix, 2.0);
		gsl_matrix_scale(&R2_ELEM.matrix, 2.0);

		gsl_matrix *sum = gsl_matrix_alloc(np, 4);
		gsl_matrix_set_zero(sum);                                                         // initialize sum as zero
		gsl_matrix_add(sum, &R0_ELEM.matrix);                                                      // add two terms into sum
		gsl_matrix_add(sum, &R1_ELEM.matrix);
		gsl_matrix_add(sum, &R2_ELEM.matrix);
		gsl_matrix_add(sum, &R3_ELEM.matrix);

		gsl_matrix *inc = gsl_matrix_alloc(np, 4);                                      // the increasement of state of elem i
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0/6.0*mesh->Elems[i].dt, M_inv, sum, 0.0, inc);
		
		gsl_matrix_add(&STATE_ELEM.matrix, inc);
		// free memory
		gsl_matrix_free(M_inv);
		gsl_matrix_free(sum);
		gsl_matrix_free(inc);
	}
}

void auto_save(MESH *mesh, gsl_matrix *state, string fname)
{
	int p = mesh->p;
	int np = (p + 1)*(p + 2) / 2;
	ofstream output;
	output.open(fname);
	output.precision(16);
	for (int i = 0; i < mesh->nElemTot * np; i++)
	{
		for (int j = 0; j < 4; j++)
			output << gsl_matrix_get(state,i,j) << " ";
		output << endl;
	}
	output.close();
}

void save_history(double *R_L_inf, int nTime_step_curr, string fname)
{
	ofstream output;
	output.open(fname);
	output.precision(16);
	output << nTime_step_curr << endl;
	for (int i = 0; i < nTime_step_curr; i++)
		output << R_L_inf[i] << endl;
	output.close();
}

void time_iteration(MESH *mesh, gsl_matrix *state, int nTime_step, string fname, double *R_L_inf)
{
	int i,j;
	int p = mesh->p;
	int np = (p + 1)*(p + 2) / 2;
	gsl_matrix *state0 = gsl_matrix_alloc(mesh->nElemTot * np, 4);
	gsl_matrix *state1 = gsl_matrix_alloc(mesh->nElemTot * np, 4);
	gsl_matrix *state2 = gsl_matrix_alloc(mesh->nElemTot * np, 4);
	gsl_matrix *R0 = gsl_matrix_alloc(mesh->nElemTot * np, 4);
	gsl_matrix *R1 = gsl_matrix_alloc(mesh->nElemTot * np, 4);
	gsl_matrix *R2 = gsl_matrix_alloc(mesh->nElemTot * np, 4);
	gsl_matrix *R3 = gsl_matrix_alloc(mesh->nElemTot * np, 4);

	for (i = 0;i < nTime_step;i++)
	{
		// MESH *mesh, gsl_matrix *state_in, gsl_matrix *state_out, gsl_matrix *R, bool timestep_update, double timestep_coef)

		// cout << "KR4_STEP1" << endl;
		R_L_inf[i] = RK4_intermediate_steps(mesh, state, state0, R0, 1, 0.5);
		// cout << "KR4_STEP2" << endl;
		RK4_intermediate_steps(mesh, state0, state1, R1, 0, 0.5);
		// cout << "KR4_STEP3" << endl;
		RK4_intermediate_steps(mesh, state1, state2, R2, 0, 1.0);
		// cout << "KR4_STEP4" << endl;
		refresh_R(mesh, state2, R3);
		// cout << "KR4_FINAL_STEP" << endl;
		RK4_final_step(mesh, state, R0, R1, R2, R3);
		if (R_L_inf[i] < 1.0e-7)
		{
			cout << "iteration = " << i << ", |R|_L_inf = " << R_L_inf[i] << endl;
			auto_save(mesh, state, fname + "_state_final" + ".txt");
			save_history(R_L_inf, i, fname + "_history" + ".txt");
			break;
		}
		if (i % 500 == 0)
		{
			cout << "iteration = " << i << ", |R|_L_inf = " << R_L_inf[i] << endl;
			auto_save(mesh, state, fname + "_state_iter_" + to_string(i) + ".txt");
		}	
	}
	gsl_matrix_free(state0);
	gsl_matrix_free(state1);
	gsl_matrix_free(state2);
	gsl_matrix_free(R0);
	gsl_matrix_free(R1);
	gsl_matrix_free(R2);
	gsl_matrix_free(R3);
}