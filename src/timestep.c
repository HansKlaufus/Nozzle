/*
** Function TimeStep
** Calculates the appropiate timestep
**
** In:       tData   Data   = structure containing all data
** Out:      tResult Result = structure containing results
** Return:   0 on success, -1 on failure
**
** Author:   J.L. Klaufus
*/

#include <stdio.h>
#include <math.h>

#include "main.h"
#include "timestep.h"

int TimeStep(FILE *log, tData *Data, tResult *Result)
{
	int ret;
	int i;

	double CFL;
	double X1, X2;

	double gamma;
	double A, rho, u, Et, p, a;

	double localTimeStep;

	/*printf("Calculating timestep...\n");*/

	ret = 0;

	for (i=1; i<Result->im-1; i++)
	{
		gamma = Data->gamma;
		CFL   = Data->CFL;

		/* Solve for primitives */
		X1    = Result->x[i];
		X2    = Result->x[i+1];
		A     = Result->A[i];
		rho   = Result->Q1[i]/A;
		u     = Result->Q2[i]/Result->Q1[i];
		Et    = Result->Q3[i]/A;
		p     = (Et-0.5*rho*u*u)*(gamma-1);
		a     = sqrt(gamma*p/rho);

		/* Calculate timestep */
		if (a < SMALL)
		{
			fprintf(stderr, "ERROR in function TimeStep: a = %10.3f\n", a);
			ret = -1;
		}
		else
		{
			/* Calculate deltaT */
			localTimeStep = CFL*(X2 - X1)/(fabs(u)+a);

			if ((i==1) || (localTimeStep < Result->timeStep))
				Result->timeStep = localTimeStep;
		}
	}

	if (log)
	{
		fprintf(log, "\n***** FUNCTION TIMESTEP *****\n\n");

		fprintf(log, "Timestep = %f\n", Result->timeStep);

		fprintf(log, "\n*****************************\n\n");
	}

	return ret;
}

