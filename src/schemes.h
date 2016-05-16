/*
** Header-file for Schemes
*/

#ifndef SCHEMES_H
#define SCHEMES_H

int    Constant(FILE*, tResult*, int, tConservative*, tConservative*);
int    Muscl(FILE*, tData*, tResult*, int, tConservative*, tConservative*);
double VanLeer(double);
double VanAlbada(double);
double KappaScheme(double, double);

#endif
