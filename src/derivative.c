#include <stdio.h>
#include <math.h>

#include "main.h"
#include "derivative.h"

/*
** Function Derivative
**   Calculates the derivative of the area function.
**
** In:      log    = name of logfile
**          i      = node number
**          Result = structure containing results
** Out:     -
** Return:  dA_dx  = derivative of area function at node i
**
** Author: J.L. Klaufus
*/

double Derivative(FILE *log, int i, tResult *Result)
{
	double dA_dx, dA_dx_old;
	double X1, X2;
	double A1, A2;

	if (i < Result->im-1)
	{
		/* Find x-coordinates */
		X1 = Result->x[i];
		X2 = Result->x[i+1];

		/* Find areas */
		A1 = AREA(X1);
		A2 = AREA(X2);

		/* Compute derivative */
		dA_dx = (A2-A1)/(X2-X1);

		dA_dx_old = 0;
		while (fabs(dA_dx - dA_dx_old) > SMALL)
		{
			X2 = X1 + (X2-X1)/10;
			A2 = AREA(X2);

			dA_dx_old = dA_dx;
			dA_dx     = (A2-A1)/(X2-X1);
		}
	}
	else
	{
		/* Find x-coordinates */
		X1 = Result->x[i-1];
		X2 = Result->x[i];

		/* Find areas */
		A1 = AREA(X1);
		A2 = AREA(X2);

		/* Compute derivative */
		dA_dx = (A2-A1)/(X2-X1);

		dA_dx_old = 0;
		while (fabs(dA_dx - dA_dx_old) > SMALL)
		{
			X1 = X2 - (X2-X1)/10;
			A1 = AREA(X1);

			dA_dx_old = dA_dx;
			dA_dx     = (A2-A1)/(X2-X1);
		}
	}

	/*
	X1 = Result->x[i];
	dA_dx = 0.347*0.8/cosh(0.8*X1-4);
	*/

		/* Write report */
	if (log)
	{
		fprintf(log, "\n***** FUNCTION DERIVATIVE *****\n\n");

		fprintf(log, "  dA_dx = %10.3f\n", dA_dx);

		fprintf(log, "\n*******************************\n\n");
	}

	return dA_dx;
}

