/*
    writer.h
    a header file of matrix writer
*/
#ifndef WRITER_H
#define WRITER_H
#include <stdio.h>

extern void lda_write(FILE *tp, double **theta, int nclass, int ndoc);
void write_matrix(FILE *fp, double **matrix, int rows, int cols);
#endif
