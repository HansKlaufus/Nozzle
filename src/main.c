/*
** Function Main
** Main routine of program Nozzle
**
** Author:   J.L. Klaufus
*/

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "main.h"
#include "boundary.h"
#include "data.h"
#include "eh.h"
#include "initialise.h"
#include "maccormack.h"
#include "memory.h"
#include "roe.h"
#include "timestep.h"

int main(int argc, char *argv[])
{
	int    ret;
	int    i;
	int    debug;
	int    down;

	double residual, normResidual, oldResidual;

	time_t t1, t2;

	FILE   *logFile      = NULL;
	FILE   *residualFile = NULL;
	char   dataFileName[50];

	tData   Data;
	tResult Result;

	printf("\nStarting program Main...\n");

	ret        = 0;
	debug      = 0;
	strcpy(dataFileName, "nozzle.in");

	/* get  commandline arguments */
	for(i=1; i<argc; i++)
	{
		if (strcmp(argv[i], "-l") == 0)
		{
			/* Turn logging on */
			logFile = fopen("nozzle.log", "w");
			debug   = 1;
		}
		else if (strcmp(argv[i], "-f") == 0)
		{
			/* Use different datafile */
			strcpy(dataFileName, argv[++i]);
		}
		else
		{
			printf("\nUnknown commandline option: '%s'\n", argv[i]);
			printf("Use : nozzle [-l] [-p G|V|B] [-f FILENAME]\n");
			ret = -1;
		}
	}

	/* Check if logFile was opened succesfully */
	if (logFile == NULL && debug == 1)
	{
		fprintf(stderr, "ERROR in function Main: Could not open logFile: nozzle.log.\n");
		ret = -1;
	}
	
	/* Open file for residual and check for success */
	residualFile = fopen("residual.gnu", "w");
	if (residualFile == NULL)
	{
		fprintf(stderr, "ERROR in function Main: Could not open residualFile: residual.log.\n");
		ret = -1;
	}
	else
		fprintf(residualFile, "#   I   Residual\n");

	if (ret != -1)
	{
		/* Read data from file */
		if (ret != -1)
			ret = ReadData(logFile, dataFileName,  &Data);

		/* Allocate memory */
		if (ret != -1)
			ret = InitMem(logFile, &Data, &Result);

		/* Set start time */
		t1 = time(&t1);

		/* Initialise */
		if (ret != -1)
			ret = Init(logFile, &Data, &Result);

		i            = 0;
		down         = 0;
		oldResidual  = 0;
		normResidual = 0;
		residual     = SMALL+1;
		while ((residual > SMALL) && (ret != -1))
		{
			i++;

			/* Calculate E and H vectors */
			if (ret != -1)
				ret = CalcEH(logFile, &Data, &Result);

			/* Calculate timestep */
			if (ret != -1)
				ret = TimeStep(logFile, &Data, &Result);

			/* Solve */
			if (ret != -1)
			{
				oldResidual = residual;

				if (Data.scheme == 'C')
					ret = MacCormack(logFile, &Data, &Result, &residual);
				else if (Data.scheme == 'R')
					ret = Roe(logFile, &Data, &Result, &residual);
				else if (Data.scheme == 'M')
					ret = Roe(logFile, &Data, &Result, &residual);
				else
				{
					fprintf(stderr, "ERROR in function Main: UNKNOWN scheme type...\n");
					fprintf(logFile,"ERROR in function Main: UNKNOWN scheme type...\n");
					ret = -1;
				}
			}

			/* Update boundaries */
			if (ret != -1)
				ret = Boundary(logFile, &Data, &Result);

			/* Normalise residual */
			if (i==1)
				normResidual = residual;

			residual /= normResidual;

			/* Write the residual */
			if ((residual < oldResidual) && (i>1))
				down++;

			fprintf(stderr,       "I = %d Residual = %10.7f [DECREASING = %d%%]\n", i, residual, (int)((float)(100*down)/i));

			/*
			if ((residual > oldResidual) && (i>1))
			{
				fprintf(stderr,       "I = %d Residual = %10.7f INCREASING\n", i, residual);
			}
			else
			{
				fprintf(stderr,       "I = %d Residual = %10.7f DECREASING\n", i, residual);
			}
			*/

			fprintf(residualFile, "%5d %10.7f\n", i, residual);
		}
		printf("Iterations  : %d\n", i);

		/* Set end time */
		t2 = time(&t2);
		printf("Calculation time = %d sec.\n", (int) (t2-t1));

		/* Write the data to outputfile */
		//if (ret != -1)
			ret = WriteData(logFile, &Data, &Result);

		/* Free allocated memory */
		if (ret != -1)
			ret = FreeMem(&Result);

	}

	/* Close files if opened */
	if (logFile)
	{
		fclose(logFile);

		if (ret == -1)
			printf("\nERRORS occurred, please read LOG.\n");

		printf("Log can be found in nozzle.log\n");
	}

	if (residualFile)
		fclose(residualFile);


	printf("Done.\n\n");

	return ret;
}

