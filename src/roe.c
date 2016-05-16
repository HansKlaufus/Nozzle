/*
** Function Roe
**    Uses Roe's approximate Riemann solver for calculating the 
**    flow characteristics in a quasi-onedimensional flow.
**
** In:       FILE    log      = pointer to log file
**           tData   Data     = structure containing all data
** Out:      tResult Result   = structure containing results
**           double  residual = residual for convergence testing
** Return:   0 on success, -1 on failure
**
** Author:   J.L. Klaufus
*/

#include <stdio.h>
#include <math.h>

#include "main.h"
#include "roe.h"
#include "schemes.h"

int Roe(FILE *log, tData *Data, tResult *Result, double *residual)
{
	int ret;
	int i, im;

	double gamma;
	double epsilon;
	double timeStep, tau;
	double rhoBefore, rhoAfter;

	tConservative left, right;

	double rho_l, rho_r;
	double u_l, u_r;
	double Et_l, Et_r;
	double p_l, p_r;
	double H_l, H_r;
	double R;
	double rho_tilde, u_tilde, H_tilde, a_tilde;
	double rhoDelta, uDelta, pDelta;
	double alpha_1, alpha_2, alpha_3;
	double lambda_tilde_1, lambda_tilde_2, lambda_tilde_3;
	double R1_tilde[3], R2_tilde[3], R3_tilde[3];
	double E_l[3], E_r[3];
	double E_tilde_right[3], E_tilde_left[3];


	ret = 0;

	im        = Data->im;
	gamma     = Data->gamma;
	epsilon   = Data->epsilon;
	timeStep  = Result->timeStep;

	*residual = 0;
	for(i=0; i<im-1; i++)
	{
		if (Data->scheme == 'R')
		{
			/* First order approximation */
			if (ret!=-1)
				ret = Constant(&(*log), &(*Result), i, &left, &right);
		}
		else if (Data->scheme == 'M')
		{
			/* Use the MUSCL-scheme */
			if (ret!=-1)
				ret = Muscl(&(*log), &(*Data), &(*Result), i, &left, &right);
		}
		else
		{
			fprintf(stderr, "ERROR in function ROE: Unknown scheme...\n");
			fprintf(log   , "ERROR in function ROE: Unknown scheme...\n");
			ret = -1;
		}


		/* Decode the primitives */
		rho_l  = left.Q1;
		rho_r  = right.Q1;

		u_l    = left.Q2/left.Q1;
		u_r    = right.Q2/right.Q1;

		Et_l   = left.Q3;
		Et_r   = right.Q3;

		p_l    = (Et_l - 0.5*rho_l*u_l*u_l)*(gamma-1);
		p_r    = (Et_r - 0.5*rho_r*u_r*u_r)*(gamma-1);

		H_l    = (Et_l + p_l)/rho_l;
		H_r    = (Et_r + p_r)/rho_r;

		/* Calculate Roe averaged values */
		R         = sqrt(rho_r/rho_l);
		rho_tilde = R*rho_l;
		u_tilde   = (u_l + R*u_r)/(1+R);
		H_tilde   = (H_l + R*H_r)/(1+R);
		a_tilde   = sqrt((gamma-1)*(H_tilde - 0.5*u_tilde*u_tilde));

		/* Calculate deltas */
		rhoDelta = rho_r - rho_l;
		uDelta   = u_r   - u_l;
		pDelta   = p_r   - p_l;

		/* Calculate wave strengths */
		alpha_1  = (pDelta - rho_tilde*a_tilde*uDelta)/(2*a_tilde*a_tilde);
		alpha_2  = rhoDelta - pDelta/(a_tilde*a_tilde);
		alpha_3  = (pDelta + rho_tilde*a_tilde*uDelta)/(2*a_tilde*a_tilde);

		/* Calculate averaged eigenvalues */
		lambda_tilde_1 = fabs(u_tilde - a_tilde);
		lambda_tilde_2 = fabs(u_tilde);
		lambda_tilde_3 = fabs(u_tilde + a_tilde);

		/* Entropy fix by Harten and Hyman */
		if (lambda_tilde_1 < epsilon)
			lambda_tilde_1 = 0.5*(lambda_tilde_1/epsilon + epsilon);

		if (lambda_tilde_2 < epsilon)
			lambda_tilde_2 = 0.5*(lambda_tilde_2/epsilon + epsilon);

		if (lambda_tilde_3 < epsilon)
			lambda_tilde_3 = 0.5*(lambda_tilde_3/epsilon + epsilon);

		/* Calculate averaged eigenvectors */
		R1_tilde[0] = 1.0;
		R1_tilde[1] = u_tilde - a_tilde;
		R1_tilde[2] = H_tilde - u_tilde*a_tilde;

		R2_tilde[0] = 1.0;
		R2_tilde[1] = u_tilde;
		R2_tilde[2] = 0.5*u_tilde*u_tilde;

		R3_tilde[0] = 1.0;
		R3_tilde[1] = u_tilde + a_tilde;
		R3_tilde[2] = H_tilde + u_tilde*a_tilde;

		/* Calculate averaged flux through interface right of current node */
		E_l[0] = Result->E1[i];
		E_l[1] = Result->E2[i];
		E_l[2] = Result->E3[i];

		E_r[0] = Result->E1[i+1];
		E_r[1] = Result->E2[i+1];
		E_r[2] = Result->E3[i+1];

		E_tilde_right[0] = 0.5*(E_l[0]+E_r[0]) - 
		                   0.5*(alpha_1*lambda_tilde_1*R1_tilde[0] + 
		                        alpha_2*lambda_tilde_2*R2_tilde[0] + 
		                        alpha_3*lambda_tilde_3*R3_tilde[0]);

		E_tilde_right[1] = 0.5*(E_l[1]+E_r[1]) - 
		                   0.5*(alpha_1*lambda_tilde_1*R1_tilde[1] + 
		                        alpha_2*lambda_tilde_2*R2_tilde[1] + 
		                        alpha_3*lambda_tilde_3*R3_tilde[1]);

		E_tilde_right[2] = 0.5*(E_l[2]+E_r[2]) - 
		                   0.5*(alpha_1*lambda_tilde_1*R1_tilde[2] + 
		                        alpha_2*lambda_tilde_2*R2_tilde[2] + 
		                        alpha_3*lambda_tilde_3*R3_tilde[2]);

		/* Calculate new Q-vector      */
		/* Only in inner field, so i>0 */
		if (i>0)
		{
			/* Use density for residual calculation */
			rhoBefore = Result->Q1[i]/Result->A[i];
			
			tau = timeStep/(Result->x[i+1]-Result->x[i]);
			Result->Q1[i] += -tau*(E_tilde_right[0] - E_tilde_left[0]);
			Result->Q2[i] += -tau*(E_tilde_right[1] - E_tilde_left[1]) + timeStep*Result->H2[i];
			Result->Q3[i] += -tau*(E_tilde_right[2] - E_tilde_left[2]);

			/* Calculate the residual */
			rhoAfter    = Result->Q1[i]/Result->A[i];
			*residual  += pow((rhoAfter-rhoBefore)/timeStep, 2);
		}

		/* Store E_tilde_right as E_tilde_left for next node */
		E_tilde_left[0] = E_tilde_right[0];
		E_tilde_left[1] = E_tilde_right[1];
		E_tilde_left[2] = E_tilde_right[2];
	}


	/* Write report */
	if (log)
	{
		fprintf(log, "\n***** FUNCTION ROE *****\n\n");

		if (ret != -1)
		{
			fprintf(log, "   I         Q1         Q2         Q3\n");
			for (i=0; i<im; i++)
				fprintf(log, " %3d %10.4f %10.4f %10.4f\n", i, Result->Q1[i], Result->Q2[i], Result->Q3[i]);
		}
		else
		{
			fprintf(log, " Function Roe NOT succesfully ended.\n");
		}

		fprintf(log, "\n*********************************\n\n");
	}

	return ret;
}

