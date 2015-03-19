/*
 * writer.c
 */
#include <stdio.h>
#include "writer.h"

void write_matrix(FILE *fp, double *mat, int n_rows, int n_cols){
    int i,j;
    for(i = 0;i < n_rows;i++)
        for(j = 0;j < n_cols;j++)
            fprintf(fp,"%.7e%s",mat[i*n_cols+j],(j == n_cols-1)? "\n":"\t");
}

void write_2d_matrix(FILE *fp, double **mat, int n_rows, int n_cols){
    int i,j;
    for(i = 0;i < n_rows;i++)
        for(j = 0;j < n_cols;j++)
            fprintf(fp,"%.7e%s",mat[i][j],(j == n_cols-1)? "\n":"\t");
}
