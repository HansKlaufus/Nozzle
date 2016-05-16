/*
** Header-file for Data
*/

#ifndef DATA_H
#define DATA_H

int ReadData(FILE*, char*, tData*);
int WriteData(FILE*, tData*, tResult*);
int WriteVigieData(FILE*, tResult*);
int WriteGNUData(FILE*, tData*, tResult*);

#endif
