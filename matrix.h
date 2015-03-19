/*
 * matrix.h
 */
#ifndef MATRIX_H
#define MATRIX_H

extern double determinant(double *mat, int n_rows);
extern void invert(double *mat, double *res, int n_rows);
extern void get_row(double *mat, double *res, int n_rows, int n_cols, int row_index);
extern void get_col(double *mat, double *res, int n_rows, int n_cols, int col_index);
extern double vec_max(double *v, int len);
extern void vec_subtract(double *v1, double *v2, double *res, int len);
extern void matrix_add(double *m1, double *m2, double *res, int n_rows, int n_cols);
extern void prod_mat_vec(double *mat, double *vec, double *res, int n_rows, int n_cols);
extern double dot_to_scalar(double *v1, double *v2, int len);
extern void dot_to_mat(double *v1, double *v2, double *res, int len);

#endif
