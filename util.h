/*
 * util.h
 */
#ifndef UTIL_H
#define UTIL_H

extern double **dmatrix(char *filename, int *n_samples, int *n_features);
extern double **matrix(int n_samples, int n_features);
extern void free_matrix(double **matrix, int n_samples);
extern void print_matrix(double *d, int n_rows, int n_cols);
extern void print_vector(double *d, int len);

#endif
