/*
** Header-file for Main
*/

#ifndef MAIN_H
#define MAIN_H

#define SMALL    1e-7
#define AREA(x)  (1.398 + 0.347*tanh(0.8*x - 4))

typedef struct
{
	double length;
	int    im;
	
	char   scheme;
	double CFL;
	double epsilon;
	double kappa;

	double gamma;
	double R;

	double M_start;
	double p_start;
	double rho_start;
	double T_start;
	double a_start;
	double u_start;

	double u_exit;

	double T_0;
	double a_0;
	double p_0;
	double rho_0;
} tData;

typedef struct
{
	int      im;

	double   timeStep;

	double   *Q1, *Q2, *Q3;
	double   *E1, *E2, *E3;
	double   *H2;

	double   *x;
	double   *A;
} tResult;

#endif
