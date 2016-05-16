/*
** Function CalcEH
**   Calculates the vectors E ans H.
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
#include "derivative.h"
#include "eh.h"

int CalcEH(FILE *log, tData *Data, tResult *Result)
{
	int ret;
	int i;

	double gamma;
	double A, rho, u, p, Et;

	/*printf("Calculating vectors E and H...\n");*/

	ret = 0;

	for (i=0; i<Result->im; i++)
	{
		/* Solve for primitives */
		gamma = Data->gamma;
		A     = Result->A[i];
		rho   = Result->Q1[i]/A;
		u     = Result->Q2[i]/Result->Q1[i];
		Et    = Result->Q3[i]/A;
		p     = (Et-0.5*rho*u*u)*(gamma-1);

		/* Calculate the vectors */
		Result->E1[i] = rho*u*A;
		Result->E2[i] = (rho*u*u + p)*A;
		Result->E3[i] = u*(Et + p)*A;
		Result->H2[i] = p*Derivative(&(*log), i, &(*Result));
	}

	if (log)
	{
		fprintf(log, "\n***** FUNCTION CALCEH *****\n\n");

		fprintf(log, "  I         E1         E2         E3         H2\n");
		for (i=0; i<Result->im; i++)
			fprintf(log, "%3d %10.4f %10.4f %10.4f %10.4f\n", i, Result->E1[i], Result->E2[i], Result->E3[i], Result->H2[i]);

		fprintf(log, "\n*****************************\n\n");
	}

	return ret;
}
