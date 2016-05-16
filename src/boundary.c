/*
** Function Boundary
**   Updates exit boundary using linear extrapolation.
**
** In:
** Out:
** Return:  0 = success / -1 = failure
**
** Author:  J.L. Klaufus
*/

#include <stdio.h>
#include <math.h>

#include "main.h"
#include "boundary.h"

int Boundary(FILE *log, tData *Data, tResult *Result)
{
	int ret;

	int    im, ii;

	double gamma;
	double X1,   X2,   X3;
	double A1,   A2,   A3;
	double u1,   u2,   u3;
	double rho1, rho2, rho3;
	double p1,   p2,   p3;
	double Et1,  Et2;

	ret = 0;

	/*printf("Updating boundary...\n");*/

	/* Get values */
	im       = Data->im-1;
	gamma    = Data->gamma;

	/* Node before previous node */
	ii   = im-2;
	X1   = Result->x[ii];
	A1   = Result->A[ii];
	rho1 = Result->Q1[ii]/A1;
	u1   = Result->Q2[ii]/Result->Q1[ii];
	Et1  = Result->Q3[ii]/A1;
	p1   = (Et1-0.5*rho1*u1*u1)*(gamma-1);

	/* Previous node */
	ii   = im-1;
	X2   = Result->x[ii];
	A2   = Result->A[ii];
	rho2 = Result->Q1[ii]/A2;
	u2   = Result->Q2[ii]/Result->Q1[ii];
	Et2  = Result->Q3[ii]/A2;
	p2   = (Et2-0.5*rho2*u2*u2)*(gamma-1);

	/* Fixed exit values */
	X3   = Result->x[im];
	A3   = Result->A[im];
	u3   = Data->u_exit;

	/* Use linear extrapolation to calculate density and pressure */
	rho3 = rho2 + (rho2-rho1)/(X2-X1)*(X3-X2);
	p3   = p2   + (p2-p1)/(X2-X1)*(X3-X2);

	/* Store values */
	Result->Q1[im] = rho3*A3;
	Result->Q2[im] = rho3*A3*u3;
	Result->Q3[im] = (0.5*rho3*u3*u3 + p3/(gamma-1))*A3;

	/* Write report */
	if (log)
	{
		fprintf(log, "\n***** FUNCTION BOUNDARY *****\n\n");

		fprintf(log, "  X1   = %10.3f\n", X1);
		fprintf(log, "  X2   = %10.3f\n", X2);
		fprintf(log, "  X3   = %10.3f\n", X3);

		fprintf(log, "  A1   = %10.3f\n", Result->A[im-2]);
		fprintf(log, "  A2   = %10.3f\n", Result->A[im-1]);
		fprintf(log, "  A3   = %10.3f\n", Result->A[im]);

		fprintf(log, "  u1   = %10.3f\n", u1);
		fprintf(log, "  u2   = %10.3f\n", u2);
		fprintf(log, "  u3   = %10.3f\n", u3);

		fprintf(log, "  rho1 = %10.3f\n", rho1);
		fprintf(log, "  rho2 = %10.3f\n", rho2);
		fprintf(log, "  rho3 = %10.3f\n", rho3);

		fprintf(log, "  p1   = %10.3f\n", p1);
		fprintf(log, "  p2   = %10.3f\n", p2);
		fprintf(log, "  p3   = %10.3f\n", p3);

		fprintf(log, "  Et1  = %10.3f\n", Et1);
		fprintf(log, "  Et2  = %10.3f\n", Et2);

		fprintf(log, "  Q1   = %10.3f\n", Result->Q1[im]);
		fprintf(log, "  Q2   = %10.3f\n", Result->Q2[im]);
		fprintf(log, "  Q3   = %10.3f\n", Result->Q3[im]);

		fprintf(log, "\n*****************************\n\n");
	}

	return ret;
}
