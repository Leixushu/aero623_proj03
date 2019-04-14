void define_U(double M, double alpha, double gamma, double U[4])
{
	U[0] = 1;
	U[1] = M*cos(alpha);
	U[2] = M*sin(alpha);
	U[3] = 1 / (gamma - 1) + 0.5*M*M*gamma;
	// U[3] = 1 / ((gamma - 1)*gamma) + 0.5*M*M; // <--- AE 523 Proj 2 ?
}

void initialization(MESH *mesh, gsl_matrix *state)
{
	double U_inf[4];
	define_U(M_inf, alpha, gamma, U_inf);
	gsl_vector_view U_INF = gsl_vector_view_array(U_inf, 4);
	int p = mesh->p;
	int np = (p + 1)*(p + 2) / 2;
	for (int i = 0;i < mesh->nElemTot * np;i++)
		gsl_matrix_set_row(state, i, &U_INF.vector);
}

void load_state(MESH *mesh, string fname, gsl_matrix *state)
{
	ifstream input(fname);
	double s;
	int p = mesh->p;
	int np = (p + 1)*(p + 2) / 2;
	for (int i = 0; i < mesh->nElemTot * np; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			input >> s;
			gsl_matrix_set(state, i, j, s);
		}		
	}
}