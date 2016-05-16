/*
** Function InitMem
** Initialises all arrays in structure Result
**
** In:       tData Data = structure containing all data
** Out:      -
** Return:   0 on success, -1 on failure
**
** Author:   J.L. Klaufus
*/

#include <stdio.h>
#include <stdlib.h>

#include "main.h"
#include "memory.h"

int InitMem(FILE *log, tData *Data, tResult *Result)
{
	int ret;

	printf("Allocating memory...\n");

	ret = 0;

	Result->im = Data->im;

	Result->x   = (double*)malloc(Result->im*sizeof(double));
	Result->A   = (double*)malloc(Result->im*sizeof(double));

	Result->Q1  = (double*)malloc(Result->im*sizeof(double));
	Result->Q2  = (double*)malloc(Result->im*sizeof(double));
	Result->Q3  = (double*)malloc(Result->im*sizeof(double));

	Result->E1  = (double*)malloc(Result->im*sizeof(double));
	Result->E2  = (double*)malloc(Result->im*sizeof(double));
	Result->E3  = (double*)malloc(Result->im*sizeof(double));

	Result->H2  = (double*)malloc(Result->im*sizeof(double));

	if((Result->x == NULL)  || (Result->A == NULL)  ||
	   (Result->Q1 == NULL) || (Result->Q2 == NULL) || (Result->Q3 == NULL) ||
	   (Result->E1 == NULL) || (Result->E2 == NULL) || (Result->E3 == NULL) ||
	   (Result->H2 == NULL))
	{
		if (log)
			fprintf(log, "ERROR in function InitMem: could not allocate memory...\n");

		ret = -1;
	}

	if (log)
	{
		fprintf(log, "\n***** FUNCTION INITMEM *****\n\n");

		if (ret == 0)
			fprintf(log, "Function InitMem succesfully ended.\n");
		else
			fprintf(log, "Function InitMem NOT succesfully ended.\n");

		fprintf(log, "\n****************************\n\n");
	}

	return ret;
}

/*
** Function FreeMem
** Deallocates all memory
**
** In:       tData   Data     = structure containing data
**           tResult Result   = structure containing results
**           
** Out:      -
**
** Return:   0 on success, -1 on failure
**
** Author:   J.L. Klaufus
*/

int FreeMem(tResult *Result)
{
	int ret = 0;

	printf("Deallocating memory...\n");

	if (Result->x)
		free(Result->x);

	if (Result->A)
		free(Result->A);

	if (Result->Q1)
		free(Result->Q1);

	if (Result->Q2)
		free(Result->Q2);

	if (Result->Q3)
		free(Result->Q3);

	if (Result->E1)
		free(Result->E1);

	if (Result->E2)
		free(Result->E2);

	if (Result->E3)
		free(Result->E3);

	if (Result->H2)
		free(Result->H2);

	return ret;
}

