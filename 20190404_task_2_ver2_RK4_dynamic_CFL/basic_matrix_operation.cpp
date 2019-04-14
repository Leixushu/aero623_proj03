void find_Jacobian(int q, double xi, double eta, double **node_pos, gsl_matrix *J) // node_pos is the nodes position for a single elem
{
	// basis fcns correspond to geo order
	int nq = (q + 1)*(q + 2) / 2;
	int i, j, k; 
	double J_jk;
	gsl_matrix *grad_psi = gsl_matrix_alloc(nq, 2);
    find_2d_grad_basis_fcn(q, xi, eta, grad_psi); // nq * 2 matrix, 1st column: partial / partial xi, 2nd column: partial / partial eta
	gsl_matrix_set_zero(J);
	for (int i = 0; i < nq; i++)                                             // i th node
	{
		for (j = 0; j < 2; j++)                                              // j th row of J
		{
			for (k = 0; k < 2; k++)                                          // k th column of J
			{
				J_jk = gsl_matrix_get(J, j, k) + node_pos[i][j] * gsl_matrix_get(grad_psi,i,k);
				gsl_matrix_set(J, j, k, J_jk);
			}
		}
	}
	// free memory
	gsl_matrix_free(grad_psi);
}

void find_inv(gsl_matrix *A, gsl_matrix *A_inv)
{
	int n = A->size1;
	gsl_matrix *A_new = gsl_matrix_alloc(n, n);    // A_new is a n by n matrix
	gsl_matrix_memcpy(A_new, A);                   // copy A to A_new                                          
	int s;
	gsl_permutation *p = gsl_permutation_alloc(n);
	gsl_linalg_LU_decomp(A_new, p, &s);            // A_new will also be modified
	gsl_linalg_LU_invert(A_new, p, A_inv);
	// free memory
	gsl_matrix_free(A_new);
	gsl_permutation_free(p);
}

double find_det(gsl_matrix *J)
{
	gsl_matrix *J_new = gsl_matrix_alloc(2, 2);
	gsl_matrix_memcpy(J_new, J);                                           // copy J to J_new
	int n = J_new->size1;                                                  // J_new is a n by n matrix
	int s;
	gsl_permutation *p = gsl_permutation_alloc(n);
	gsl_linalg_LU_decomp(J_new, p, &s);                                    // J_new will also be modified
	double J_det = gsl_linalg_LU_det(J_new, s);
	// free memory
	gsl_matrix_free(J_new);
	gsl_permutation_free(p);
	return J_det;
}