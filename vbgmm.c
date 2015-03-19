/*
 * vbgmm.c
 * Variational Bayesian Gaussian Mixture Model
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "vbgmm.h"
#include "util.h"
#include "writer.h"
#include "learn.h"

int main(int argc, char **argv){
    double **data;
    int n_class = CLASS_DEFAULT;
    int n_iter  = MAX_ITER_DEFAULT;
    char c;
    int i;
    double *m;
    double **l;
    double *r;
    FILE *mp, *lp, *rp; // for m (mean of cluster prior), lambda (precision: mean of wishart) and r (response)
    int n_samples,n_features;
    
    while((c = getopt(argc,argv,"K:I:h")) != -1){
        switch(c){
            case 'K': n_class = atoi(optarg); break;
            case 'I': n_iter = atoi(optarg); break;
            case 'h': usage(); break;
            default : usage(); break;
        }
    }
    if(!(argc - optind == 2))
        usage();
    
    /* open data */
    if((data = dmatrix(argv[optind], &n_samples, &n_features)) == NULL){
        fprintf(stderr,"vbgmm:: cannot open training data.\n");
        exit(1);
    }
    /* allocate parameters */
    if((m = (double *)calloc(n_class*n_features,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm:: cannot allocate m.\n");
        exit(1);
    }
    if((l = (double **)calloc(n_class,sizeof(double *))) == NULL){
        fprintf(stderr,"vbgmm:: cannot allocate lambda (top level).\n");
        exit(1);
    }
    for(i = 0;i < n_class;i++){
        if((l[i] = (double *)calloc(n_features*n_features,sizeof(double))) == NULL){
            fprintf(stderr,"vbgmm:: cannot allocate lambda (second level).\n");
            exit(1);
        }
    }
    if((r = (double *)calloc(n_samples*n_class,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm:: cannot allocate r.\n");
        exit(1);
    }
    /* open model output*/
    char m_name[strlen(argv[optind+1])+2];
    char l_name[strlen(argv[optind+1])+2];
    char r_name[strlen(argv[optind+1])+2];
    sprintf(m_name,"%s.m",argv[optind+1]);
    sprintf(l_name,"%s.l",argv[optind+1]);
    sprintf(r_name,"%s.r",argv[optind+1]);
    if(((mp = fopen(m_name,"w")) == NULL)
    || ((lp = fopen(l_name,"w")) == NULL)
    || ((rp = fopen(r_name,"w")) == NULL)){
        fprintf(stderr,"vbgmm:: cannot open model outputs.\n");
        exit(1);
    }
    
    vbgmm_learn(data, m, l, r, n_class, n_samples, n_features, n_iter, argv[optind+1]);
    
    write_matrix(mp,m,n_class,n_features);
    write_2d_matrix(lp,l,n_class,(n_features*n_features));
    write_matrix(rp,r,n_samples,n_class);
    
    free_matrix(data,n_samples);
    free(m);
    for(i = 0;i < n_class;i++){
        free(l[i]);
    }
    free(l);
    free(r);
    fclose(mp);
    fclose(lp);
    fclose(rp);
    
    exit(0);
}

void usage(){
    printf("usage: %s -K classes -I iterations train model\n","vbgmm");
    exit(0);
}
