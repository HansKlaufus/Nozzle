#include <stdio.h>
#include <math.h>

#include "main.h"
#include "av.h"

/*
** Function CalcAV
**    Creates the source term using artificial viscosity at a 
**    given node.
**
** In:       FILE    log     = logfile
**           double  gamma   = constant
**           double  epsilon = constant factor for artificial viscosity
**           int     i       = current node number
**           int     im      = last node number
**           double  A       = area at current node
**           double  Q1      = first unknowns of vector Q
**           double  Q2      = second unknowns of vector Q
**           double  Q3      = third unknowns of vector Q
** Out:      tAV     AV      = structure containing results S1, S2 and S3.
** Return:   int     ret     = 0 on success / -1 on failure
**
** Author:   J.L. Klaufus
*/

int CalcAV(FILE *log, double gamma, double epsilon, int i, int im, double A, double *Q1, double *Q2, double* Q3, tAV *AV)
{
	int    ret;
	int    ii;
	double rho_prev, rho_cur, rho_next;
	double u_prev,   u_cur,   u_next;
	double p_prev,   p_cur,   p_next;
	double a_prev,   a_cur,   a_next;
	double p_term;
	double u_plus,   u_min;
	double av_plus,  av_min;
	double Q1_prev,  Q1_cur,  Q1_next;
	double Q2_prev,  Q2_cur,  Q2_next;
	double Q3_prev,  Q3_cur,  Q3_next;

	ret = 0;

	if (A < SMALL)
	{
		fprintf(stderr, "ERROR in function AV: A = %10.4f\n", A);
		ret = -1;
	}
	else
	{
		/*
		** Find the pressures
		*/

		/* Current node */
		ii      = i;
		rho_cur = Q1[ii]/A;
		u_cur   = Q2[ii]/Q1[ii];
		p_cur   = (Q3[ii]/A - 0.5*rho_cur*u_cur*u_cur)*(gamma-1);
		a_cur   = sqrt(p_cur*gamma/rho_cur);

		if (i==im-1)
		{
			/* Previous node */
			ii       = i-1;
			rho_prev = Q1[ii]/A;
			u_prev   = Q2[ii]/Q1[ii];
			p_prev   = (Q3[ii]/A - 0.5*rho_prev*u_prev*u_prev)*(gamma-1);

			if (u_prev < 0) u_prev = u_cur;
			if (p_prev < 0) p_prev = p_cur;

			a_prev   = sqrt(p_prev*gamma/rho_prev);

			/* Next node */
			/*   Does not exist; use linear interpolation */
			rho_next = 2*rho_cur - rho_prev;
			u_next   = 2*u_cur - u_prev;
			p_next   = 2*p_cur - p_prev;

			if (u_next < 0) u_next = u_cur;
			if (p_next < 0) p_next = p_cur;

			a_next   = sqrt(p_next*gamma/rho_next);
		}
		else if (i==0)
		{
			/* Next node */
			ii       = i+1;
			rho_next = Q1[ii]/A;
			u_next   = Q2[ii]/Q1[ii];
			p_next   = (Q3[ii]/A - 0.5*rho_next*u_next*u_next)*(gamma-1);

			if (u_next < 0) u_next = u_cur;
			if (p_next < 0) p_next = p_cur;

			a_next   = sqrt(p_next*gamma/rho_next);

			/* Previous node */
			/*   Does not exist; use linear interpolation */
			rho_prev = 2*rho_cur - rho_next;
			u_prev   = 2*u_cur - u_next;
			p_prev   = 2*p_cur - p_next;

			if (u_prev < 0) u_prev = u_cur;
			if (p_prev < 0) p_prev = p_cur;

			a_prev   = sqrt(p_prev*gamma/rho_prev);
		}
		else
		{
			/* Previous node */
			ii       = i-1;
			rho_prev = Q1[ii]/A;
			u_prev   = Q2[ii]/Q1[ii];
			p_prev   = (Q3[ii]/A - 0.5*rho_prev*u_prev*u_prev)*(gamma-1);

			if (u_prev < 0) u_prev = u_cur;
			if (p_prev < 0) p_prev = p_cur;

			a_prev   = sqrt(p_prev*gamma/rho_prev);

			/* Next node */
			ii       = i+1;
			rho_next = Q1[ii]/A;
			u_next   = Q2[ii]/Q1[ii];
			p_next   = (Q3[ii]/A - 0.5*rho_next*u_next*u_next)*(gamma-1);

			if (u_next < 0) u_next = u_cur;
			if (p_next < 0) p_next = p_cur;

			a_next   = sqrt(p_next*gamma/rho_next);
		}

		if (p_next < 0)
			fprintf(stderr, "ERROR: P<0; i_next = %d\n", i+1);

		/* Calculate the pressure term */
		if (p_next+2*p_cur+p_prev < 0)
		{
			fprintf(stderr, "ERROR in function AV: p_term division by zero\n");
			fprintf(log   , "ERROR in function AV: p_term division by zero\n");
			fprintf(log   , "  i = %d\n", i);
			fprintf(log   , "  P1 = %10.4f\n", p_prev);
			fprintf(log   , "  P2 = %10.4f\n", p_cur);
			fprintf(log   , "  P3 = %10.4f\n", p_next);

			p_term =  0;
			ret    = -1;
		}
		else
			p_term = fabs(p_next-2*p_cur+p_prev)/(p_next+2*p_cur+p_prev);
	
		/* Calculate velocity dependent terms */
		u_plus = (fabs(u_cur) + a_cur + fabs(u_next) + a_next)/2;
		u_min  = (fabs(u_cur) + a_cur + fabs(u_prev) + a_prev)/2;

		/* Calculate the artificial terms */
		av_plus = epsilon*u_plus*p_term;
		av_min  = epsilon*u_min*p_term;

		/* Calculate the correct Q-terms */
		Q1_prev = rho_prev*A;
		Q1_cur  = rho_cur*A;
		Q1_next = rho_next*A;

		Q2_prev = rho_prev*u_prev*A;
		Q2_cur  = rho_cur*u_cur*A;
		Q2_next = rho_next*u_next*A;

		Q3_prev = (0.5*rho_prev*u_prev*u_prev + p_prev/(gamma-1))*A;
		Q3_cur  = (0.5*rho_cur*u_cur*u_cur    + p_cur/(gamma-1))*A;
		Q3_next = (0.5*rho_next*u_next*u_next + p_next/(gamma-1))*A;
		
		/* Calculate the D-terms */
		AV->D1 = av_plus*(Q1_next-Q1_cur) - av_min*(Q1_cur-Q1_prev);
		AV->D2 = av_plus*(Q2_next-Q2_cur) - av_min*(Q2_cur-Q2_prev);
		AV->D3 = av_plus*(Q3_next-Q3_cur) - av_min*(Q3_cur-Q3_prev);

		/* Andersson */
		/*
		AV->D1 = epsilon*p_term*(Q1_next-2*Q1_cur+Q1_prev);
		AV->D2 = epsilon*p_term*(Q2_next-2*Q2_cur+Q1_prev);
		AV->D3 = epsilon*p_term*(Q3_next-2*Q3_cur+Q1_prev);
		*/

		/* Write report */
		if (log)
		{
			fprintf(log, "\n***** FUNCTION AV *****\n");

			fprintf(log, "  I        = %10d\n", i);

			fprintf(log, "  Rho_prev = %10.3f\n", rho_prev);
			fprintf(log, "  Rho_cur  = %10.3f\n", rho_cur);
			fprintf(log, "  Rho_next = %10.3f\n", rho_next);

			fprintf(log, "  P_prev   = %10.3f\n", p_prev);
				fprintf(log, "  P_cur    = %10.3f\n", p_cur);
			fprintf(log, "  P_next   = %10.3f\n", p_next);

			fprintf(log, "  U_prev   = %10.3f\n", u_prev);
			fprintf(log, "  U_cur    = %10.3f\n", u_cur);
			fprintf(log, "  U_next   = %10.3f\n", u_next);

			fprintf(log, "  A_prev   = %10.3f\n", a_prev);
			fprintf(log, "  A_cur    = %10.3f\n", a_cur);
			fprintf(log, "  A_next   = %10.3f\n", a_next);

			fprintf(log, "  U_plus   = %10.3f\n", u_plus);
			fprintf(log, "  U_min    = %10.3f\n", u_min);

			fprintf(log, "  AV_plus  = %10.3f\n", av_plus);
			fprintf(log, "  AV_min   = %10.3f\n", av_min);

			fprintf(log, "  P_term   = %10.3f\n", p_term);

			fprintf(log, "  D1       = %10.3f\n", AV->D1);
			fprintf(log, "  D2       = %10.3f\n", AV->D2);
			fprintf(log, "  D3       = %10.3f\n", AV->D3);

			fprintf(log, "\n***********************\n");
		}
	}

	return ret;
}

