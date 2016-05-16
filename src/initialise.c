/*
** Function Init
** Initialises boundary values
**
** In:       tData   Data   = structure containing all data
** Out:      tResult Result = structure containing all results
** Return:   0 on success, -1 on failure
**
** Author:   J.L. Klaufus
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "main.h"
#include "initialise.h"

int Init(FILE *log, tData *Data, tResult *Result)
{
	int ret;
	int i;

	double deltaX;

	int    im;
	double length;

	double gamma;
	double R;

	double M_start;
	double p_start;
	double rho_start;
	double T_start;
	double a_start;
	double u_start;

	double rho, u, p;

	printf("Initialising...\n");

	ret = 0;

	/* Get handy variables */
	im        = Data->im;
	length    = Data->length;

	gamma     = Data->gamma;
	R         = Data->R;

	M_start   = Data->M_start;
	p_start   = Data->p_start;
	rho_start = Data->rho_start;

	/* Calculate start conditions and store */
	a_start = sqrt(gamma*p_start/rho_start);
	u_start = M_start*a_start;
	T_start = (rho_start*u_start*u_start/2 + p_start/(gamma-1))*(gamma-1)/R;

	Data->T_start = T_start;
	Data->a_start = a_start;
	Data->u_start = u_start;
	
	/* Initialise field */
	for (i=0; i<im; i++)
	{
		/* Calculate x co-ordinates; use uniform distribution */
		deltaX = length/(im-1);
		Result->x[i] = i*deltaX;

		/* Calculate areas */
		Result->A[i] = AREA(Result->x[i]);

		/* Quess starting conditions */
		rho = rho_start;
		p   = p_start;

		if (i != im-1)
			u = u_start;
		else
			u = Data->u_exit;

		/* Initialise Q-vector */
		Result->Q1[i] = rho*Result->A[i];
		Result->Q2[i] = rho*u*Result->A[i];
		Result->Q3[i] = (rho*u*u/2 + p/(gamma-1))*Result->A[i];
	}

	/* Write report */
	if (log)
	{
		fprintf(log, "\n***** FUNCTION INIT *****\n\n");

		fprintf(log, " T_start   = %10.3f\n", Data->T_start);
		fprintf(log, " a_start   = %10.3f\n", Data->a_start);
		fprintf(log, " u_start   = %10.3f\n\n", Data->u_start);

		fprintf(log, "  I         Q1         Q2         Q3\n");
		for (i=0; i<im; i++)
			fprintf(log, "%3d %10.3f %10.3f %10.3f\n", i, Result->Q1[i], Result->Q2[i], Result->Q3[i]);

		fprintf(log, "\n*************************\n\n");
	}

	return ret;
}

