/*
 * writer.h
 */
#ifndef WRITER_H
#define WRITER_H
#include <stdio.h>

extern void write_matrix(FILE *fp, double *mat, int n_rows, int n_cols);
extern void write_2d_matrix(FILE *fp, double **mat, int n_rows, int n_cols);

#endif
