// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently

#pragma once

#include "targetver.h"
#include <stdio.h>
#include <tchar.h>

#include <cstdlib>  // system("pause") need this!
#include <iostream> // cin cout need this!
#include <fstream>  // file input and output
#include <cmath>
#include <algorithm>

// TODO: reference additional headers your program requires here
using namespace std;
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "my_class.h"
const double R = 1.0;
const double gamma = 1.4;
const double M_inf = 0.5;
const double alpha = 0.0;
const double p_inf = 1.0;
const double CFL_ref = 0.85;
const double h = 0.0625;