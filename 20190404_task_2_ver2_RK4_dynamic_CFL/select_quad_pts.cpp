#include "quad1d.c"
#include "quad2d.c"

// note that # of quad pts should be more than 2p + 1

void select_1d_quad_pts(int p, int q, QUADRATURE_1D *quad)
{
	// p = 0,1,2; q = 1,3; The order is determined by 2*p + 2*(q-1) + 1
	int i;
	switch (2*p + 2*(q-1) + 1)
	{
	case 1:
	{
		quad->n = n1_1D;
		quad->x = new double[quad->n];
		quad->w = new double[quad->n];
		for (i = 0; i < quad->n; i++)
		{
			quad->x[i] = x1_1D[i];
			quad->w[i] = w1_1D[i];
		}
	}
	break;
	case 3:
	{
		quad->n = n3_1D;
		quad->x = new double[quad->n];
		quad->w = new double[quad->n];
		for (i = 0; i < quad->n; i++)
		{
			quad->x[i] = x3_1D[i];
			quad->w[i] = w3_1D[i];
		}
	}
	break;
	case 5:
	{
		quad->n = n5_1D;
		quad->x = new double[quad->n];
		quad->w = new double[quad->n];
		for (i = 0; i < quad->n; i++)
		{
			quad->x[i] = x5_1D[i];
			quad->w[i] = w5_1D[i];
		}
	}
	break;
	case 7:
	{
		quad->n = n7_1D;
		quad->x = new double[quad->n];
		quad->w = new double[quad->n];
		for (i = 0; i < quad->n; i++)
		{
			quad->x[i] = x7_1D[i];
			quad->w[i] = w7_1D[i];
		}
	}
	break;
	case 9:
	{
		quad->n = n9_1D;
		quad->x = new double[quad->n];
		quad->w = new double[quad->n];
		for (i = 0; i < quad->n; i++)
		{
			quad->x[i] = x9_1D[i];
			quad->w[i] = w9_1D[i];
		}
	}
	break;
	}
}

void select_2d_quad_pts(int p, int q, QUADRATURE_2D *quad)
{
	int i;
	switch (2 * p + 2 * (q - 1) + 1)
	{
	case 1:
	{
		quad->n = n1_2D;
		quad->x = new double[quad->n];
		quad->y = new double[quad->n];
		quad->w = new double[quad->n];
		for (i = 0; i < quad->n; i++)
		{
			quad->x[i] = x1_2D[2 * i];
			quad->y[i] = x1_2D[2 * i + 1];
			quad->w[i] = w1_2D[i];
		}
	}
	break;
	case 3:
	{
		quad->n = n3_2D;
		quad->x = new double[quad->n];
		quad->y = new double[quad->n];
		quad->w = new double[quad->n];
		for (i = 0; i < quad->n; i++)
		{
			quad->x[i] = x3_2D[2 * i];
			quad->y[i] = x3_2D[2 * i + 1];
			quad->w[i] = w3_2D[i];
		}
	}
	break;
	case 5:
	{
		quad->n = n5_2D;
		quad->x = new double[quad->n];
		quad->y = new double[quad->n];
		quad->w = new double[quad->n];
		for (i = 0; i < quad->n; i++)
		{
			quad->x[i] = x5_2D[2 * i];
			quad->y[i] = x5_2D[2 * i + 1];
			quad->w[i] = w5_2D[i];
		}
	}
	break;
	case 7:
	{
		quad->n = n7_2D;
		quad->x = new double[quad->n];
		quad->y = new double[quad->n];
		quad->w = new double[quad->n];
		for (i = 0; i < quad->n; i++)
		{
			quad->x[i] = x7_2D[2 * i];
			quad->y[i] = x7_2D[2 * i + 1];
			quad->w[i] = w7_2D[i];
		}
	}
	break;
	case 9:
	{
		quad->n = n9_2D;
		quad->x = new double[quad->n];
		quad->y = new double[quad->n];
		quad->w = new double[quad->n];
		for (i = 0; i < quad->n; i++)
		{
			quad->x[i] = x9_2D[2 * i];
			quad->y[i] = x9_2D[2 * i + 1];
			quad->w[i] = w9_2D[i];
		}
	}
	break;
	}
}