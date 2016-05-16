/*
** Header-file for AV
*/

#ifndef AV_H
#define AV_H

typedef struct
{
	double D1;
	double D2;
	double D3;
} tAV;

int CalcAV(FILE*, double, double, int, int, double, double*, double*, double*, tAV*);

#endif
