#pragma once
#include <string>
#include <vector>

class EDGE
{
public:
	// Data Members 
	int nNodes;
	int *nodes;
	int t1;
	int e1;
	int t2;             // if not exist -1
	int e2;             // if not exist -1
	double norm_vec[2]; // only for linear elem, when curved, set as -1 when pre-computing
	double length;      // only for linear elem, when curved, set as -1 when pre-computing
	double s;           // wave speed
	string type;        // Left / Right / Bottom / Top / Interior
};

class ADJ_CELL
{
public:
	// Data Members
	int loc_edge;
	int adj_edge;
	int adj_tri;
};

class CELL
{
public:
	// Data Members 
	int q;                     // geometry order
	int *nodes;
	int edge[3];
	vector<ADJ_CELL> adj_cell;
	gsl_matrix *M;             // mass matrix
	// double s_max;              // max wave speed inside a cell
	double dt;                 // local time step
	double A;
	// double P;
};

class EGROUP
{
public:
	// Data Members 
	int nElem;
	int Order;
	string Basis;
	int *E;
};

class BGROUP
{
public:
	// Data Members 
	int nBFace;   // number of faces on that boundary
	int nf;
	string Title;
	int *B;       // the address of B[0]
};

class MESH
{
public:
	// Data Members 
	int nNode;
	double **node_pos;
	int nElemTot;
	int nEGroup;
	vector<EGROUP> EGroup;
	int nBGroup;
	BGROUP *BGroup;
	CELL *Elems;
	int nEdge;
	vector<EDGE> Edges;
	int p; // order of solution
};

class EDGE_TEMP
{
public:
	// Data Members
	int *nodes;  // nodes list
	int t;
	int e;
};

class INPUT_TEMP
{
public:
	// Data Members
	int ***B_temp;       // 3D array nBGroup * nBFace * nf, where nf = number of nodes per face
	EDGE_TEMP **edge_temp; // 2D array nElemTot*3
};

class QUADRATURE_1D
{
public:
	// Data Members
	int n;
	double *x;
	double *w;
	void clear()
	{
		delete[] x;
		delete[] w;
	}
};

class QUADRATURE_2D
{
public:
	// Data Members
	int n;
	double *x;
	double *y;
	double *w;
	void clear()
	{
		delete[] x;
		delete[] y;
		delete[] w;
	}
};

class PROP               // key properties of the state
{
public:
	double F[2][4];      // Euler flux
	double rho;          // density
	double v[2];         // velocity
	double c;            // speed of sound
	double E;
	double p;
	double H;
};