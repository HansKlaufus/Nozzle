#include <stdio.h>
#include <math.h>

#include "main.h"
#include "roe.h"
#include "schemes.h"

/*
** Function Constant
**    Calculates the left and right conservative variables
**    for use in an approximate riemann solver.
**
** In:       FILE          log      = pointer to log file
**           tResult       Result   = structure containing results
**           int           i        = number of node left of inerface of interest
** Out:      tConservative left     = structure containing left conservative variables
**           tConservative right    = structure containing right conservative variables
** Return:   0 on success, -1 on failure
**
** Author:   J.L. Klaufus
*/

int Constant(FILE *log, tResult *Result, int i, tConservative *left, tConservative *right)
{
	int ret;

	ret = 0;

	left->Q1  = Result->Q1[i];
	right->Q1 = Result->Q1[i+1];

	left->Q2  = Result->Q2[i];
	right->Q2 = Result->Q2[i+1];

	left->Q3  = Result->Q3[i];
	right->Q3 = Result->Q3[i+1];

	/* Write report */
	if (log)
	{
		fprintf(log, "\n******* FUNCTION CONSTANT *******\n");

		fprintf(log, "               Q1         Q2         Q3\n");
		fprintf(log, "Left : %10.4f %10.4f %10.4f\n", left->Q1, left->Q2, left->Q3);
		fprintf(log, "Right: %10.4f %10.4f %10.4f\n", right->Q1, right->Q2, right->Q3);

		fprintf(log, "\n******************************\n");
	}

	return ret;
}


/*
** Function Muscl
**    Calculates the left and right conservative variables
**    for use in an approximate riemann solver.
**
** In:       FILE          log      = pointer to log file
**           tData         Data     = structure containing all data
**           tResult       Result   = structure containing results
**           int           i        = number of node left of inerface of interest
** Out:      tConservative left     = structure containing left conservative variables
**           tConservative right    = structure containing right conservative variables
** Return:   0 on success, -1 on failure
**
** Author:   J.L. Klaufus
*/

int Muscl(FILE *log, tData *Data, tResult *Result, int i, tConservative *left, tConservative *right)
{
	int ret;

	double gamma;
	double kappa;

	double Q1[4], Q2[4], Q3[4];
	double rho0, rho1, rho2, rho3;
	double u0, u1, u2, u3;
	double p0, p1, p2, p3;

	double r;
	double psi_l, psi_r;

	ret = 0;

	gamma     = Data->gamma;
	kappa     = Data->kappa;

	if (i==0)
	{
		/* Result->Q1[i-1] does not exist; use linear extrapolation */
		rho2  = Result->Q1[i+1];
		rho1  = Result->Q1[i  ];
		rho0  = 2*rho1 - rho2;

		u2    = Result->Q2[i+1]/Result->Q1[i+1];
		u1    = Result->Q2[i  ]/Result->Q1[i  ];
		u0    = 2*u1 - u2;

		p2    = (Result->Q3[i+1] - 0.5*rho2*u2*u2)*(gamma-1);
		p1    = (Result->Q3[i  ] - 0.5*rho1*u1*u1)*(gamma-1);
		p0    = 2*p1 - p2;

		Q1[0] = rho0;
		Q1[1] = Result->Q1[i  ]; 
		Q1[2] = Result->Q1[i+1];
		Q1[3] = Result->Q1[i+2];

		Q2[0] = rho0*u0;
		Q2[1] = Result->Q2[i  ]; 
		Q2[2] = Result->Q2[i+1];
		Q2[3] = Result->Q2[i+2];

		Q3[0] = 0.5*rho0*u0*u0 + p0/(gamma-1);
		Q3[1] = Result->Q3[i  ]; 
		Q3[2] = Result->Q3[i+1];
		Q3[3] = Result->Q3[i+2];
	}
	else if (i==Data->im-1)
	{
		/* Result->Q1[i+2] does not exist; use linear extrapolation */
		rho1  = Result->Q1[i  ];
		rho2  = Result->Q1[i+1];
		rho3  = 2*rho2 - rho1;

		u1    = Result->Q2[i  ]/Result->Q1[i  ];
		u2    = Result->Q2[i+1]/Result->Q1[i+1];
		u3    = 2*u2 - u1;

		p1    = (Result->Q3[i  ] - 0.5*rho1*u1*u1)*(gamma-1);
		p2    = (Result->Q3[i+1] - 0.5*rho2*u2*u2)*(gamma-1);
		p3    = 2*p2 - p1;

		Q1[0] = Result->Q1[i-1]; 
		Q1[1] = Result->Q1[i  ]; 
		Q1[2] = Result->Q1[i+1];
		Q1[3] = rho3;

		Q2[0] = Result->Q2[i-1]; 
		Q2[1] = Result->Q2[i  ]; 
		Q2[2] = Result->Q2[i+1];
		Q2[3] = rho3*u3;

		Q3[0] = Result->Q3[i-1]; 
		Q3[1] = Result->Q3[i  ]; 
		Q3[2] = Result->Q3[i+1];
		Q3[3] = 0.5*rho3*u3*u3 + p3/(gamma-1);
	}
	else
	{
		/* Everything OK; copy all */
		Q1[0] = Result->Q1[i-1]; 
		Q1[1] = Result->Q1[i  ]; 
		Q1[2] = Result->Q1[i+1];
		Q1[3] = Result->Q1[i+2];

		Q2[0] = Result->Q2[i-1]; 
		Q2[1] = Result->Q2[i  ]; 
		Q2[2] = Result->Q2[i+1];
		Q2[3] = Result->Q2[i+2];

		Q3[0] = Result->Q3[i-1]; 
		Q3[1] = Result->Q3[i  ]; 
		Q3[2] = Result->Q3[i+1];
		Q3[3] = Result->Q3[i+2];
	}

	/*
	** Find Q1
	*/
	if (Q1[1]-Q1[0] < SMALL)
	{
		left->Q1  = Q1[1];
	}
	else
	{
		r         = (Q1[2] - Q1[1])/(Q1[1] - Q1[0]);
		psi_l     = VanLeer(r);
		left->Q1  = Q1[1] + 0.5*psi_l*(Q1[1] - Q1[0]);
	}

	if (Q1[2]-Q1[1] < SMALL)
	{
		right->Q1 = Q1[2];
	}
	else
	{
		r         = (Q1[3] - Q1[2])/(Q1[2] - Q1[1]);
		psi_r     = VanLeer(1/r);
		right->Q1 = Q1[2] - 0.5*psi_r*(Q1[3] - Q1[2]);
	}

	/*
	** Find Q2
	*/
	if (Q2[1]-Q2[0] < SMALL)
	{
		left->Q2  = Q2[1];
	}
	else
	{
		r         = (Q2[2] - Q2[1])/(Q2[1] - Q2[0]);
		psi_l     = VanLeer(r);
		left->Q2  = Q2[1] + 0.5*psi_l*(Q2[1] - Q2[0]);
	}

	if (Q2[2]-Q2[1] < SMALL)
	{
		right->Q2 = Q2[2];
	}
	else
	{
		r         = (Q2[3] - Q2[2])/(Q2[2] - Q2[1]);
		psi_r     = VanLeer(1/r);
		right->Q2 = Q2[2] - 0.5*psi_r*(Q2[3] - Q2[2]);
	}

	/*
	** Find Q3
	*/
	if (Q3[1]-Q3[0] < SMALL)
	{
		left->Q3  = Q3[1];
	}
	else
	{
		r         = (Q3[2] - Q3[1])/(Q3[1] - Q3[0]);
		psi_l     = VanLeer(r);
		left->Q3  = Q3[1] + 0.5*psi_l*(Q3[1] - Q3[0]);
	}

	if (Q3[2]-Q3[1] < SMALL)
	{
		right->Q3 = Q3[2];
	}
	else
	{
		r         = (Q3[3] - Q3[2])/(Q3[2] - Q3[1]);
		psi_r     = VanLeer(1/r);
		right->Q3 = Q3[2] - 0.5*psi_r*(Q3[3] - Q3[2]);
	}

	/* Write report */
	if (log)
	{
		fprintf(log, "\n******* FUNCTION MUSCL *******\n");

		fprintf(log, "               Q1         Q2         Q3\n");
		fprintf(log, "Left : %10.4f %10.4f %10.4f\n", left->Q1, left->Q2, left->Q3);
		fprintf(log, "Right: %10.4f %10.4f %10.4f\n", right->Q1, right->Q2, right->Q3);

		fprintf(log, "\n******************************\n");
	}

	return ret;
}



/*
** Function VanLeer
**   Functions as a limiter to the MUSCL-scheme.
** 
**   In:      double r      = the ratio of the consecutive gradients.
**   Out:     -
**   Return:  double psi    = limiter function
**
**   Author:  J.L. Klaufus
*/

double VanLeer(double r)
{
	double psi;

	psi = (r + fabs(r))/(fabs(r) + 1);
	//psi = (r + fabs(r))/(r*r + 1);

	return psi;
}


/*
** Function VanAlbada
**   Functions as a limiter to the MUSCL-scheme.
** 
**   In:      double r      = the ratio of the consecutive gradients.
**   Out:     -
**   Return:  double psi    = limiter function
**
**   Author:  J.L. Klaufus
*/

double VanAlbada(double r)
{
	double psi;

	psi = (r*r + r)/(r*r + 1);

	return psi;
}


/*
** Function KappaScheme
**   Functions as a limiter to the MUSCL-scheme.
** 
**   In:      double r      = the ratio of the consecutive gradients.
**            double kappa  = value for kappa
**   Out:     -
**   Return:  double psi    = limiter function
**
**   Author:  J.L. Klaufus
*/

double KappaScheme(double r, double kappa)
{
	double phi;
	double psi;

	phi = (2*r)/(r*r + 1);
	psi = ((1-kappa)/2 + r*(1+kappa)/2)*phi;

	return psi;
}

