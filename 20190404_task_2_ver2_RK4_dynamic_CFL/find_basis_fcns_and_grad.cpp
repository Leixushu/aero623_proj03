/*
Function: int gsl_blas_sgemm (CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, float alpha, const gsl_matrix_float * A, const gsl_matrix_float * B, float beta, gsl_matrix_float * C)
Function: int gsl_blas_dgemm(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, double alpha, const gsl_matrix * A, const gsl_matrix * B, double beta, gsl_matrix * C)
Function : int gsl_blas_cgemm(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, const gsl_complex_float alpha, const gsl_matrix_complex_float * A, const gsl_matrix_complex_float * B, const gsl_complex_float beta, gsl_matrix_complex_float * C)
Function : int gsl_blas_zgemm(CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, const gsl_complex alpha, const gsl_matrix_complex * A, const gsl_matrix_complex * B, const gsl_complex beta, gsl_matrix_complex * C)
These functions compute the matrix - matrix product and sum C = \alpha op(A) op(B) + \beta C where op(A) = A, A^T, A^H for TransA = CblasNoTrans, CblasTrans, CblasConjTrans and similarly for the parameter TransB.
*/

void find_2d_basis_fcn(int p, double xi, double eta, gsl_vector *psi)
{
	// gsl_vector *psi = gsl_vector_alloc(np);
	int N, i, j;
	int np = (p + 1)*(p + 2) / 2;
	gsl_vector_set_zero(psi);
	switch (p)
	{
	case 0:
		gsl_vector_set(psi, 0, 1.0000000000000000);
	break;
	case 1:
	{
		double m[] = { 
			1.0000000000000000, -1.0000000000000000, -1.0000000000000000,
			0.0000000000000000, 1.0000000000000000, 0.0000000000000000,
			-0.0000000000000000, -0.0000000000000000, 1.0000000000000000 };
		double basis[] = { 1.0, xi, eta };
		gsl_matrix_view M = gsl_matrix_view_array(m, np, np);
		gsl_vector_view BASIS = gsl_vector_view_array(basis, np);
		gsl_blas_dgemv(CblasNoTrans, 1.0, &M.matrix, &BASIS.vector, 0.0, psi);
	}
	break;
	case 2:
	{
		double m[] = { 
			1.0000000000000000, -3.0000000000000000, 2.0000000000000000, -3.0000000000000000, 4.0000000000000000, 2.0000000000000000,
			0.0000000000000000, 4.0000000000000000, -4.0000000000000000, 0.0000000000000000, -4.0000000000000000, 0.0000000000000000,
			0.0000000000000000, -1.0000000000000000, 2.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
			0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 4.0000000000000000, -4.0000000000000000, -4.0000000000000000,
			-0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, 4.0000000000000000, 0.0000000000000000,
			0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -1.0000000000000000, 0.0000000000000000, 2.0000000000000000 };
		double basis[] = { 1.0, xi, xi*xi, eta, xi*eta, eta*eta };
		gsl_matrix_view M = gsl_matrix_view_array(m, np, np);
		gsl_vector_view BASIS = gsl_vector_view_array(basis, np);
		gsl_blas_dgemv(CblasNoTrans, 1.0, &M.matrix, &BASIS.vector, 0.0, psi);
	}
	break;
	case 3:
	{
		double m[] = { 
			1.0000000000000000, -5.4999999999999982, 8.9999999999999947, -4.4999999999999964, -5.4999999999999982, 18.0000000000000000, -13.5000000000000000, 8.9999999999999947, -13.5000000000000000, -4.4999999999999964,
			0.0000000000000000, 8.9999999999999964, -22.4999999999999858, 13.4999999999999911, 0.0000000000000000, -22.5000000000000000, 27.0000000000000000, 0.0000000000000000, 13.5000000000000000, 0.0000000000000000,
			0.0000000000000000, -4.4999999999999964, 17.9999999999999858, -13.4999999999999911, 0.0000000000000000, 4.5000000000000000, -13.5000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
			0.0000000000000000, 0.9999999999999987, -4.4999999999999947, 4.4999999999999964, 0.0000000000000000, 0.0000000000000005, 0.0000000000000013, 0.0000000000000000, -0.0000000000000005, 0.0000000000000000,
			-0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, 8.9999999999999964, -22.5000000000000000, 13.5000000000000000, -22.4999999999999858, 27.0000000000000000, 13.4999999999999911,
			0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 27.0000000000000000, -27.0000000000000000, -0.0000000000000000, -27.0000000000000000, 0.0000000000000000,
			0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -4.5000000000000000, 13.5000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
			0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -4.4999999999999964, 4.5000000000000000, 0.0000000000000000, 17.9999999999999858, -13.5000000000000000, -13.4999999999999911,
			0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -4.5000000000000000, 0.0000000000000000, 0.0000000000000000, 13.5000000000000000, 0.0000000000000000,
			0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.9999999999999987, 0.0000000000000005, -0.0000000000000006, -4.4999999999999947, -0.0000000000000013, 4.4999999999999964 };
		double basis[] = { 1.0, xi, xi*xi, xi*xi*xi, eta, xi*eta,  xi*xi*eta, eta*eta, xi*eta*eta, eta*eta*eta };
		gsl_matrix_view M = gsl_matrix_view_array(m, np, np);
		gsl_vector_view BASIS = gsl_vector_view_array(basis, np);
		gsl_blas_dgemv(CblasNoTrans, 1.0, &M.matrix, &BASIS.vector, 0.0, psi);
	}
	break;
	}
}

void find_2d_grad_basis_fcn(int p, double xi, double eta, gsl_matrix *grad_psi)
{
	// gsl_matrix *grad_psi = gsl_matrix_alloc(np, 2);
	int N, i, j;
	int np = (p + 1)*(p + 2) / 2;
	gsl_vector *grad_psi_xi = gsl_vector_alloc(np);
	gsl_vector *grad_psi_eta = gsl_vector_alloc(np);
	gsl_vector_set_zero(grad_psi_xi);
	gsl_vector_set_zero(grad_psi_eta);
	switch (p)
	{
	case 0:
	{
		// psi[0]       = 1;
		gsl_vector_set(grad_psi_xi, 0, 0.0000000000000000); // partial psi[0] / partial xi;
		gsl_vector_set(grad_psi_eta, 0, 0.0000000000000000); // partial psi[0] / partial eta;
	}
	break;
	case 1:
	{
		double m[] = { 
			1.0000000000000000, -1.0000000000000000, -1.0000000000000000,
			0.0000000000000000, 1.0000000000000000, 0.0000000000000000,
			-0.0000000000000000, -0.0000000000000000, 1.0000000000000000 };
		// basis[3]                = { 1, xi, eta};
		double basis_partial_xi[]  = { 0.0,  1.0,  0.0 }; // partial basis / partial xi;
		double basis_partial_eta[] = { 0.0,  0.0,  1.0 }; // partial basis / partial eta;
		gsl_matrix_view M = gsl_matrix_view_array(m, np, np);
		gsl_vector_view BASIS_PARTIAL_XI = gsl_vector_view_array(basis_partial_xi, np);
		gsl_vector_view BASIS_PARTIAL_ETA = gsl_vector_view_array(basis_partial_eta, np);
		gsl_blas_dgemv(CblasNoTrans, 1.0, &M.matrix, &BASIS_PARTIAL_XI.vector, 0.0, grad_psi_xi);
		gsl_blas_dgemv(CblasNoTrans, 1.0, &M.matrix, &BASIS_PARTIAL_ETA.vector, 0.0, grad_psi_eta);
	}
	break;
	case 2:
	{
		double m[] = { 
			1.0000000000000000, -3.0000000000000000, 2.0000000000000000, -3.0000000000000000, 4.0000000000000000, 2.0000000000000000,
			0.0000000000000000, 4.0000000000000000, -4.0000000000000000, 0.0000000000000000, -4.0000000000000000, 0.0000000000000000,
			0.0000000000000000, -1.0000000000000000, 2.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
			0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 4.0000000000000000, -4.0000000000000000, -4.0000000000000000,
			-0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, 4.0000000000000000, 0.0000000000000000,
			0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -1.0000000000000000, 0.0000000000000000, 2.0000000000000000 };
		//  basis[6]               = { 1, xi,  xi*xi, eta, xi*eta, eta*eta };
		double basis_partial_xi[]  = { 0.0,  1.0, 2.0 * xi,   0.0,    eta, 0.0 };       // partial basis / partial xi;
		double basis_partial_eta[] = { 0.0,  0.0,      0.0,   1.0,     xi, 2.0 * eta }; // partial basis / partial eta;
		gsl_matrix_view M = gsl_matrix_view_array(m, np, np);
		gsl_vector_view BASIS_PARTIAL_XI = gsl_vector_view_array(basis_partial_xi, np);
		gsl_vector_view BASIS_PARTIAL_ETA = gsl_vector_view_array(basis_partial_eta, np);
		gsl_blas_dgemv(CblasNoTrans, 1.0, &M.matrix, &BASIS_PARTIAL_XI.vector, 0.0, grad_psi_xi);
		gsl_blas_dgemv(CblasNoTrans, 1.0, &M.matrix, &BASIS_PARTIAL_ETA.vector, 0.0, grad_psi_eta);
	}
	break;
	case 3:
	{
		double m[] = { 
			1.0000000000000000, -5.4999999999999982, 8.9999999999999947, -4.4999999999999964, -5.4999999999999982, 18.0000000000000000, -13.5000000000000000, 8.9999999999999947, -13.5000000000000000, -4.4999999999999964,
			0.0000000000000000, 8.9999999999999964, -22.4999999999999858, 13.4999999999999911, 0.0000000000000000, -22.5000000000000000, 27.0000000000000000, 0.0000000000000000, 13.5000000000000000, 0.0000000000000000,
			0.0000000000000000, -4.4999999999999964, 17.9999999999999858, -13.4999999999999911, 0.0000000000000000, 4.5000000000000000, -13.5000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
			0.0000000000000000, 0.9999999999999987, -4.4999999999999947, 4.4999999999999964, 0.0000000000000000, 0.0000000000000005, 0.0000000000000013, 0.0000000000000000, -0.0000000000000005, 0.0000000000000000,
			-0.0000000000000000, -0.0000000000000000, -0.0000000000000000, -0.0000000000000000, 8.9999999999999964, -22.5000000000000000, 13.5000000000000000, -22.4999999999999858, 27.0000000000000000, 13.4999999999999911,
			0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 27.0000000000000000, -27.0000000000000000, -0.0000000000000000, -27.0000000000000000, 0.0000000000000000,
			0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -4.5000000000000000, 13.5000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
			0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -4.4999999999999964, 4.5000000000000000, 0.0000000000000000, 17.9999999999999858, -13.5000000000000000, -13.4999999999999911,
			0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, -4.5000000000000000, 0.0000000000000000, 0.0000000000000000, 13.5000000000000000, 0.0000000000000000,
			0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.9999999999999987, 0.0000000000000005, -0.0000000000000006, -4.4999999999999947, -0.0000000000000013, 4.4999999999999964 };
		//  basis[10]              = { 1, xi, xi*xi, xi*xi*xi, eta, xi*eta,  xi*xi*eta, eta*eta, xi*eta*eta, eta*eta*eta};
		double basis_partial_xi[]  = { 0.0,  1.0,  2.0 * xi,  3.0 * xi*xi,   0.0,    eta,   2.0 * xi*eta,         0.0,        eta*eta, 0.0 }; // partial basis / partial xi;
		double basis_partial_eta[] = { 0.0,  0.0,       0.0,          0.0,   1.0,     xi,          xi*xi,   2.0 * eta,   2.0 * xi*eta, 3.0 * eta*eta }; // partial basis / partial eta;
		gsl_matrix_view M = gsl_matrix_view_array(m, np, np);
		gsl_vector_view BASIS_PARTIAL_XI = gsl_vector_view_array(basis_partial_xi, np);
		gsl_vector_view BASIS_PARTIAL_ETA = gsl_vector_view_array(basis_partial_eta, np);
		gsl_blas_dgemv(CblasNoTrans, 1.0, &M.matrix, &BASIS_PARTIAL_XI.vector, 0.0, grad_psi_xi);
		gsl_blas_dgemv(CblasNoTrans, 1.0, &M.matrix, &BASIS_PARTIAL_ETA.vector, 0.0, grad_psi_eta);
	}
	break;
	}
	gsl_matrix_set_col(grad_psi, 0, grad_psi_xi);
	gsl_matrix_set_col(grad_psi, 1, grad_psi_eta);
	
	// free memory
	gsl_vector_free(grad_psi_xi);
	gsl_vector_free(grad_psi_eta);
}