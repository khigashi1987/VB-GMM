/*
 * lb.h
 */
#ifndef LB_H
#define LB_H
#include <stdlib.h>

extern double lower_bound(double alpha0, double *alpha, double beta0, double *beta, double nyu0, double *nyu, double *W0, double **W, double *r, int n_samples, int n_class, int n_features);

#endif
