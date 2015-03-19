/*
 * util.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "util.h"
#define BUFFSIZE 65536
static int file_lines(FILE *fp);
static int line_fields(char *line);
static int isspaces(char *s);

double **dmatrix(char *filename, int *n_samples, int *n_features){
    double **d;
    int n,m,ni,mi;
    FILE *fp;
    char line[BUFFSIZE];
    char First_Line = 1;
    
    if((fp = fopen(filename,"r")) == NULL)
        return NULL;
    n = file_lines(fp);
    if((d = (double **)calloc(n,sizeof(double *))) == NULL)
        return NULL;
    
    ni = 0;
    while(fgets(line,sizeof(line),fp)){
        int i;
        char *sp,*lp = line;
        
        if(isspaces(line))
            continue;
        
        if(First_Line){
            m = line_fields(line);
            for(i = 0;i < n;i++){
                if((d[i] = (double *)calloc(m,sizeof(double))) == NULL)
                    return NULL;
            }
            First_Line = 0;
        }
        
        mi = 0;
        while(*lp){
            if((sp = strpbrk(lp," \t\n")) == NULL)
                break;
            *sp = '\0';
            d[ni][mi] = atof(lp);
            lp = sp + 1;
            mi++;
        }
        ni++;
    }
    fclose(fp);
    *n_samples = n;
    *n_features = m;
    return d;
}

static int file_lines(FILE *fp){
    int n = 0;
    char buf[BUFFSIZE];
    
    while(fgets(buf,sizeof(buf),fp)){
        if(!isspaces(buf))
            n++;
    }
    rewind(fp);
    return n;
}

static int line_fields(char *line){
    int m = 0;
    char *sp,*lp = line;
    
    while(*lp){
        if((sp = strpbrk(lp," \t\n")) == NULL)
            break;
        lp = sp + 1;
        m++;
    }
    return m;
}

static int isspaces(char *s){
    char *c = s;
    while(*c){
        if(!isspace(*c))
            return 0;
        c++;
    }
    return 1;
}

double **matrix(int n_samples, int n_features){
    double **d;
    int i;
    
    if((d = (double **)calloc(n_samples,sizeof(double *))) == NULL)
        return NULL;
    
    for(i = 0;i < n_samples;i++){
        if((d[i] = (double *)calloc(n_features,sizeof(double))) == NULL)
            return NULL;
    }
    return d;
}

void free_matrix(double **matrix, int n_samples){
    int i;
    
    for(i = 0;i < n_samples;i++){
        free(matrix[i]);
    }
    free(matrix);
}

void print_matrix(double *d, int n_rows, int n_cols){
    int i,j;
    for(i = 0;i < n_rows;i++){
        for(j = 0;j < n_cols;j++){
            printf("%.8f\t",d[i*n_cols+j]);
        }
        printf("\n");
    }
}

void print_vector(double *d, int len){
    int i;
    for(i = 0;i < len;i++)
        printf("%.8f\t",d[i]);
    printf("\n");
}
