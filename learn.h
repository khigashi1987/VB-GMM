/*
 * learn.h
 */
#ifndef VBGMM_LEARN_H
#define VBGMM_LEARN_H
#include <stdlib.h>

extern void vbgmm_learn(double **data, double *m, double **l, double *r, int n_class, int n_samples, int n_features, int n_iter, char *model);

#endif
