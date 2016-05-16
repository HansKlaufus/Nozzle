/*
** Function ReadData.
** Reads data from data-file 'nozzle.in'
**
** In:       -
** Out:      Data = structure containing all data
** Return:   0 on success, -1 on failure
**
** Datafile: nozzle.in or user defined filename
**
** Author:   J.L. Klaufus
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "main.h"
#include "data.h"

int ReadData(FILE *log, char *dataFileName, tData *Data)
{
	FILE   *dataFile;
	int    ret;
	char   eol;

	double gamma, R;
	double M_start, p_start, rho_start;
	double u_exit;
	double length;
	char   scheme;
	double CFL;
	double epsilon, kappa;
	int    im;

	printf("Reading data...\n");

	ret = 0;

	dataFile = fopen(dataFileName, "r");
	if (dataFile)
	{
		fscanf(dataFile, "%lf %lf", &gamma,  &R);
		fscanf(dataFile, "%lf %lf %lf", &M_start, &p_start, &rho_start);
		fscanf(dataFile, "%lf", &u_exit);
		fscanf(dataFile, "%lf", &length);
		fscanf(dataFile, "%c %c", &eol, &scheme);
		fscanf(dataFile, "%lf %lf %lf", &CFL, &epsilon, &kappa);
		fscanf(dataFile, "%d", &im);

		fclose(dataFile);

		Data->gamma     = gamma;
		Data->R         = R;
		Data->M_start   = M_start;
		Data->p_start   = p_start;
		Data->rho_start = rho_start;
		Data->u_exit    = u_exit;
		Data->length    = length;
		Data->scheme    = scheme;
		Data->CFL       = CFL;
		Data->epsilon   = epsilon;
		Data->kappa     = kappa;
		Data->im        = im;

		/* Write report */
		if (log)
		{
			fprintf(log, "\n***** FUNCTION READDATA *****\n\n");

			fprintf(log, "   gamma     = %10.3f\n", Data->gamma);
			fprintf(log, "   R         = %10.3f\n", Data->R);
			fprintf(log, "   M_start   = %10.3f\n", Data->M_start);
			fprintf(log, "   p_start   = %10.3f\n", Data->p_start);
			fprintf(log, "   rho_start = %10.3f\n", Data->rho_start);
			fprintf(log, "   u_exit    = %10.3f\n", Data->u_exit);
			fprintf(log, "   length    = %10.3f\n", Data->length);
			fprintf(log, "   scheme    = %c\n", Data->scheme);
			fprintf(log, "   CFL       = %10.3f\n", Data->CFL);
			fprintf(log, "   epsilon   = %10.3f\n", Data->epsilon);
			fprintf(log, "   kappa     = %10.3f\n", Data->kappa);
			fprintf(log, "   im        = %10d\n", Data->im);

			fprintf(log, "\n*****************************\n\n");
		}
	}
	else
	{
		fprintf(stderr, "ERROR in function ReadData: Could not open '%s'.\n", dataFileName);
		ret = -1;
	}

	return ret;
}


/*
** Function WriteData.
** Writes data to files
**
** In:       tResult Result = structure containing all results
**
** Out:      -
**
** Return:   0 on success, -1 on failure
**
** Author:   J.L. Klaufus
*/

int WriteData(FILE *log, tData *Data, tResult *Result)
{
	int   ret;

	printf("Writing data...\n");

	ret = 0;

	ret = WriteGNUData(&(*log), &(*Data), &(*Result));

	if (log)
	{
		fprintf(log, "\n***** FUNCTION WRITEDATA *****\n\n");

		if (ret != -1)
			fprintf(log, "Data succesfully written to file.\n");
		else
			fprintf(log, "Data NOT succesfully written to file.\n");

		fprintf(log, "\n*****************************\n\n");
	}

	return ret;
}


/*
** Function WriteGNUData.
** Writes data in an ASCII format suitable for the visualisation
** package GNUPlot.
**
** In:       tResult Result  = structure containing all results.
**
** Out:      -
**
** Return:   0 on success, -1 on failure
**
** Datfiles: nozzle.gnu: Data file containing flow characteristics.
**
** Author:   J.L. Klaufus
*/

int WriteGNUData(FILE *log, tData *Data, tResult *Result)
{
	FILE   *dataFile = NULL;

	int    ret;
	int    i;

	double x, A;
	double gamma, R;
	double rho, u, e, T, p, a, M;

	printf("Creating GNUPlot datafile...\n"); 

	ret = 0;

	R     = Data->R;
	gamma = Data->gamma;

	/*
	** Write data for GNUPlot
	*/
	dataFile = fopen("nozzle.gnu", "w");
	if (dataFile)
	{
		fprintf(dataFile, "# I          x          A        rho          u          T          p          M\n");
		for (i=0; i<Result->im; i++)
		{
			/* Solve for primitives */
			x     = Result->x[i];
			A     = Result->A[i];
			rho   = Result->Q1[i]/A;
			u     = Result->Q2[i]/Result->Q1[i];
			e     = Result->Q3[i]/A;
			p     = (e-0.5*rho*u*u)*(gamma-1);
			T     = e*(gamma-1)/R;
			a     = sqrt(gamma*p/rho);
			M     = u/a;

			fprintf(dataFile, "%3d %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n", i, x, A, rho, u, T, p, M);
		}
	}
	else
	{
		fprintf(stderr, "Could not open output file.\n");
		ret = -1;
	}

	/* Close the files */
	if (dataFile)
		fclose(dataFile);

	/* Write report */
	if (log)
	{
		fprintf(log, "\n***** FUNCTION WRITEGNUDATA *****\n\n");

		if (ret != -1)
			fprintf(log, "GNU data succesfully written to file.\n");
		else
			fprintf(log, "GNU data NOT succesfully written to file.\n");

		fprintf(log, "\n*********************************\n\n");
	}

	return ret;
}

