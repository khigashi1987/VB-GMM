/*
 * matrix.c
 */
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include "matrix.h"

double determinant(double *mat, int n_rows){
    int s;
    double det;
    int i,j;
    gsl_matrix *m = gsl_matrix_alloc(n_rows,n_rows);
    for(i = 0;i < n_rows;i++)
        for(j = 0;j < n_rows;j++)
            gsl_matrix_set(m,i,j,mat[i*n_rows+j]);
    gsl_permutation *p = gsl_permutation_alloc(n_rows);
    gsl_linalg_LU_decomp(m,p,&s);
    det = gsl_linalg_LU_det(m,s);
    gsl_permutation_free(p);
    gsl_matrix_free(m);
    return det;
}

void invert(double *mat, double *res, int n_rows){
    int s;
    int i,j;
    gsl_matrix *gmat = gsl_matrix_alloc(n_rows,n_rows);
    for(i = 0;i < n_rows;i++)
        for(j = 0;j < n_rows;j++)
            gsl_matrix_set(gmat,i,j,mat[i*n_rows+j]);
    gsl_matrix *gmat_inv = gsl_matrix_alloc(n_rows,n_rows);
    gsl_permutation *p = gsl_permutation_alloc(n_rows);
    gsl_linalg_LU_decomp(gmat,p,&s);
    gsl_linalg_LU_invert(gmat,p,gmat_inv);
    for(i = 0;i < n_rows;i++){
        for(j = 0;j < n_rows;j++){
            res[i*n_rows+j] = gsl_matrix_get(gmat_inv,i,j);
        }
    }
    gsl_permutation_free(p);
    gsl_matrix_free(gmat);
    gsl_matrix_free(gmat_inv);
}

void get_row(double *mat, double *res, int n_rows, int n_cols, int row_index){
    int i;
    for(i = 0;i < n_cols;i++){
        res[i] = mat[row_index*n_cols+i];
    }
}
void get_col(double *mat, double *res, int n_rows, int n_cols, int col_index){
    int i;
    for(i = 0;i < n_rows;i++){
        res[i] = mat[i*n_cols+col_index];
    }
}

double vec_max(double *v, int len){
    int i;
    double max = v[0];
    for(i = 1;i < len;i++){
        if(v[i] > max) max = v[i];
    }
    return max;
}

void vec_subtract(double *v1, double *v2, double *res, int len){
    int i;
    for(i = 0;i < len;i++){
        res[i] = v1[i] - v2[i];
    }
}

void matrix_add(double *m1, double *m2, double *res, int n_rows, int n_cols){
    int i,j;
    for(i = 0;i < n_rows;i++){
        for(j = 0;j < n_cols;j++){
            res[i*n_cols+j] = m1[i*n_cols+j] + m2[i*n_cols+j];
        }
    }
}

void prod_mat_vec(double *mat, double *vec, double *res, int n_rows, int n_cols){
    int i,j;
    for(i = 0;i < n_rows;i++)
        res[i] = 0.0;
    for(i = 0;i < n_rows;i++){
        for(j = 0;j < n_rows;j++){
            res[i] += mat[i*n_cols+j] * vec[j];
        }
    }
}

double dot_to_scalar(double *v1, double *v2, int len){
    int i;
    double res = 0.0;
    for(i = 0;i < len;i++){
        res += v1[i] * v2[i];
    }
    return res;
}

void dot_to_mat(double *v1, double *v2, double *res, int len){
    int i,j;
    for(i = 0;i < len;i++){
        for(j = 0;j < len;j++){
            res[i*len+j] = v1[i] * v2[j];
        }
    }
}
