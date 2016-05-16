/*
** Header-file for Roe
*/

#ifndef ROE_H
#define ROE_H

typedef struct
{
	double Q1;
	double Q2;
	double Q3;
} tConservative;

int Roe(FILE*, tData*, tResult*, double*);

#endif
