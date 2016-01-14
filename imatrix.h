/*
    imatrix.h
    a header file of int matrix
*/
#ifndef IMATRIX_H
#define IMATRIX_H
#include <stdio.h>

extern int **load_n_wz(char *filename, int *rows, int *cols);
extern int **imatrix(int rows, int cols);
extern void free_imatrix(int **matrix, int rows);
#endif
