/*
** Function MacCormack
**    Uses the MacCormack finite difference equations for
**    calculating the flow characteristics in a quasi-
**    onedimesional flow.
**
**    CONSERVATION FORM FOR SHOCK CAPTURING
**
** In:       tData Data     = structure containing all data
** Out:      tResult Result = structure containing results
** Return:   0 on success, -1 on failure
**
** Author:   J.L. Klaufus
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "main.h"
#include "av.h"
#include "derivative.h"
#include "maccormack.h"

int MacCormack(FILE *log, tData *Data, tResult *Result, double *residual)
{
	int ret;

	int i, im;

	double A, rho, e, p, u;
	double gamma;
	double timeStep;
	double deltaX;
	double tau;
	double rhoBefore, rhoAfter;

	tAV    AV;

	double *Q1_b=NULL;
	double *Q2_b=NULL;
	double *Q3_b=NULL;

	double *Q1_bb=NULL;
	double *Q2_bb=NULL;
	double *Q3_bb=NULL;

	double *E1_b=NULL;
	double *E2_b=NULL;
	double *E3_b=NULL;

	double *H2_b=NULL;

	/*printf("Solving using MacCormack...\n");*/

	ret = 0;
	im  = Data->im;

	AV.D1 = AV.D2 = AV.D3 = 0;

	gamma    = Data->gamma;
	timeStep = Result->timeStep;

	/* Allocate memory for temporary arrays */
	Q1_b  = (double*)malloc(Data->im*sizeof(double));
	Q2_b  = (double*)malloc(Data->im*sizeof(double));
	Q3_b  = (double*)malloc(Data->im*sizeof(double));

	Q1_bb = (double*)malloc(Data->im*sizeof(double));
	Q2_bb = (double*)malloc(Data->im*sizeof(double));
	Q3_bb = (double*)malloc(Data->im*sizeof(double));

	E1_b  = (double*)malloc(Data->im*sizeof(double));
	E2_b  = (double*)malloc(Data->im*sizeof(double));
	E3_b  = (double*)malloc(Data->im*sizeof(double));

	H2_b  = (double*)malloc(Data->im*sizeof(double));

	if ((Q1_b == NULL)  || (Q2_b == NULL)  || (Q3_b == NULL) ||
	    (E1_b == NULL)  || (E2_b == NULL)  || (E3_b == NULL) ||
	    (H2_b == NULL)  ||
	    (Q1_bb == NULL) || (Q2_bb == NULL) || (Q3_bb == NULL))
	{
		fprintf(stderr, "ERROR in function MacCormack: Could not allocate memory.\n");
		ret = -1;
	}
	else
	{
		/*
		** Predictor step; using forward differencing
		**   Therefore: highest array-count: im-1
		**              highest loop-count : im-2.
		*/
		for(i=0; i<im-1 && ret!=-1; i++)
		{
			deltaX = Result->x[i+1] - Result->x[i];
			tau    = Result->timeStep/deltaX;

			/* Calculate artificial viscosity */
			ret = CalcAV(&(*log), gamma, Data->epsilon, i, im, Result->A[i], Result->Q1, Result->Q2, Result->Q3, &AV);

			/* Calculate Q-bar; defined for [0, im-2] */
			Q1_b[i] = Result->Q1[i] - tau*(Result->E1[i+1]-Result->E1[i]) + tau*AV.D1;
			Q2_b[i] = Result->Q2[i] - tau*(Result->E2[i+1]-Result->E2[i]) + tau*AV.D2 + timeStep*Result->H2[i];
			Q3_b[i] = Result->Q3[i] - tau*(Result->E3[i+1]-Result->E3[i]) + tau*AV.D3;

			/* Get primitives */
			A   = Result->A[i];
			rho = Q1_b[i]/A;
			u   = Q2_b[i]/Q1_b[i];
			e   = Q3_b[i]/A;
			p   = (e-0.5*rho*u*u)*(gamma-1);

			/* Calculate E-bar and H-bar; only defined for [0, im-2] */
			E1_b[i]   = rho*u*A;
			E2_b[i]   = (rho*u*u+p)*A;
			E3_b[i]   = u*(e+p)*A;
			H2_b[i]   = p*Derivative(&(*log), i, &(*Result))/Result->A[i];
		}

		/*
		** Corrector step; using backward differencing
		**   Only change inner field: i=1 to i=im-2
		**   Necessary for i=1      : E_b[0]
		*/
		*residual = 0;
		for(i=1; i<im-1 && ret!=-1; i++)
		{
			deltaX = Result->x[i] - Result->x[i-1];
			tau    = Result->timeStep/deltaX;

			/* Calculate artificial viscosity */
			ret = CalcAV(&(*log), gamma, Data->epsilon, i, im-1, Result->A[i], Q1_b, Q2_b, Q3_b, &AV);

			/* Calculate Q-double-bar
			**    i=1   : E1_b[0]    needed
			**    i=im-2: E1_b[im-2] needed
			*/
			Q1_bb[i] = Result->Q1[i] - tau*(E1_b[i]-E1_b[i-1]) + tau*AV.D1;
			Q2_bb[i] = Result->Q2[i] - tau*(E2_b[i]-E2_b[i-1]) + tau*AV.D2 + timeStep*H2_b[i];
			Q3_bb[i] = Result->Q3[i] - tau*(E3_b[i]-E3_b[i-1]) + tau*AV.D3;

			/* Use density for residual calculation */
			rhoBefore = Result->Q1[i]/Result->A[i];
			
			/* Calculate Q at the new timestep */
			Result->Q1[i] = 0.5*(Q1_b[i] + Q1_bb[i]);
			Result->Q2[i] = 0.5*(Q2_b[i] + Q2_bb[i]);
			Result->Q3[i] = 0.5*(Q3_b[i] + Q3_bb[i]);

			/* Calculate the residual */
			rhoAfter    = Result->Q1[i]/Result->A[i];
			*residual  += pow((rhoAfter-rhoBefore)/timeStep, 2);
		}
	}

	/* Deallocate memory */
	if (Q1_b)
		free(Q1_b);

	if (Q2_b)
		free(Q2_b);

	if (Q3_b)
		free(Q3_b);

	if (Q1_bb)
		free(Q1_bb);

	if (Q2_bb)
		free(Q2_bb);

	if (Q3_bb)
		free(Q3_bb);

	if (E1_b)
		free(E1_b);

	if (E2_b)
		free(E2_b);

	if (E3_b)
		free(E3_b);

	if (H2_b)
		free(H2_b);

	/* Write report */
	if (log)
	{
		fprintf(log, "\n***** FUNCTION MACCORMACK *****\n\n");

		if (ret != -1)
		{
			fprintf(log, "  I         Q1         Q2         Q3\n");

			for(i=1; i<Result->im-1; i++)
				fprintf(log, "%3d %10.4f %10.4f %10.4f\n", i, Result->Q1[i], Result->Q2[i], Result->Q3[i]);
		}
		else
		{
			fprintf(log, "  Function MacCormack NOT succesfully ended.\n");
		}

		fprintf(log, "\n*********************************\n\n");
	}

	return ret;
}

