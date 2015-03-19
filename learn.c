/*
 * learn.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_psi.h>
#include "learn.h"
#include "util.h"
#include "matrix.h"
#include "lb.h"

#define FIX_SEED

void vbgmm_learn(double **data, double *m, double **l, double *r, int n_class, int n_samples, int n_features, int n_iter, char *model){
    int s;
    int i,j,t;
    int k,d,n;
    double alpha0,beta0,nyu0;
    double *m0, *W0;
    double *alpha,*beta,*nyu;
    double **W;
    double *ln_lambda, *ln_pi, *ln_rho, *N_k;
    double *x_k, *W0_inv, *Wk_inv;
    double **S_k;
    double *xn_m, *wk_xn_m, *mk, *m0k, *xkk, *xn_xk, *xn_xk_mat, *xk_m0, *xk_m0_mat, *NkSk;
    double alpha_sum, quad, ln_rho_max, r_sum;
    double lowerbound;
    FILE *lbp;
    char lb_name[strlen(model)+3];
    
    
    // initialze parameters
    if((alpha = (double *)calloc(n_class,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate alpha.\n");
        return;
    }
    for(i = 0;i < n_class;i++)
        alpha[i] = 0.001 + (double)n_samples / (double)n_features;
    alpha0 = 0.001;
    
    if((beta = (double *)calloc(n_class,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate beta.\n");
        return;
    }
    for(i = 0;i < n_class;i++)
        beta[i] = 25.0 + (double)n_samples / (double)n_features;
    beta0 = 25.0;
    
    if((nyu = (double *)calloc(n_class,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate nyu.\n");
        return;
    }
    for(i = 0;i < n_class;i++)
        nyu[i] = (double)n_features + (double)n_samples / (double)n_features;
    nyu0 = (double)n_features;
    
    if((W = (double **)calloc(n_class,sizeof(double *))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate W (top-level).\n");
        return;
    }
    for(i = 0;i < n_class;i++){
        if((W[i] = (double *)calloc(n_features*n_features,sizeof(double))) == NULL){
            fprintf(stderr,"vbgmm_learn:: cannot allocate W (second-level).\n");
            return;
        }
        for(j = 0;j < (n_features*n_features);j += (n_features+1)){
            W[i][j] = 1.0; // identity matrix
        }
    }
    if((W0 = (double *)calloc(n_features*n_features,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate W0.\n");
        return;
    }
    for(i = 0;i < (n_features*n_features);i += (n_features+1)){
        W0[i] = 1.0;
    }
    
    const gsl_rng_type *T;
    gsl_rng *rn;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    rn = gsl_rng_alloc(T);
    #ifdef FIX_SEED
        gsl_rng_set(rn,1024);
    #else
        gsl_rng_set(rn,(unsigned)time(NULL));
    #endif
    for(i = 0;i < n_class;i++){
        for(j = 0;j < n_features;j++){
            m[i*n_features+j] = gsl_ran_ugaussian(rn);
        }
    }
    if((m0 = (double *)calloc(n_class*n_features,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate m0.\n");
        return;
    }
    for(i = 0;i < n_class;i++){
        for(j = 0;j < n_features;j++){
            m0[i*n_features+j] = m[i*n_features+j];
        }
    }
    gsl_rng_free(rn);
    
    // initialze buffers
    if((ln_lambda = (double *)calloc(n_class,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate ln_lambda.\n");
        return;
    }
    if((ln_pi = (double *)calloc(n_class,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate ln_pi.\n");
        return;
    }
    if((ln_rho = (double *)calloc(n_class,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate ln_rho.\n");
        return;
    }
    if((N_k = (double *)calloc(n_class,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate N_k.\n");
        return;
    }
    if((x_k = (double *)calloc(n_class*n_features,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate x_k.\n");
        return;
    }
    
    if((W0_inv = (double *)calloc(n_features*n_features,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate W0_inv.\n");
        return;
    }
    invert(W0,W0_inv,n_features);
    
    if((Wk_inv = (double *)calloc(n_features*n_features,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate Wk_inv.\n");
        return;
    }
    if((S_k = (double **)calloc(n_class,sizeof(double *))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate S_k (top-level).\n");
        return;
    }
    for(i = 0;i < n_class;i++){
        if((S_k[i] = (double *)calloc(n_features*n_features,sizeof(double))) == NULL){
            fprintf(stderr,"vbgmm_learn:: cannot allocate S_k (second-level).\n");
            return;
        }
    }
    if((xn_m = (double *)calloc(n_features,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate xn_m.\n");
        return;
    }
    if((wk_xn_m = (double *)calloc(n_features,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate wk_xn_m.\n");
        return;
    }
    if((mk = (double *)calloc(n_features,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate mk.\n");
        return;
    }
    if((m0k = (double *)calloc(n_features,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate m0k.\n");
        return;
    }
    if((xkk = (double *)calloc(n_features,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate xkk.\n");
        return;
    }
    if((xn_xk = (double *)calloc(n_features,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate xn_xk.\n");
        return;
    }
    if((xn_xk_mat = (double *)calloc(n_features*n_features,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate xn_xk_mat.\n");
        return;
    }
    if((xk_m0 = (double *)calloc(n_features,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate xk_m0.\n");
        return;
    }
    if((xk_m0_mat = (double *)calloc(n_features*n_features,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate xk_m0_mat.\n");
        return;
    }
    if((NkSk = (double *)calloc(n_features*n_features,sizeof(double))) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot allocate NkSk.\n");
        return;
    }
    
    sprintf(lb_name,"%s.lb",model);
    if((lbp = fopen(lb_name,"w")) == NULL){
        fprintf(stderr,"vbgmm_learn:: cannot open lb output.\n");
        return;
    }
    
    printf("Number of samples     = %d\n",n_samples);
    printf("Number of features    = %d\n",n_features);
    printf("Number of clusters    = %d\n",n_class);
    printf("Number of iteration   = %d\n",n_iter);
    
    /*
     * learn main
     */
    for(t = 0;t < n_iter;t++){
        //printf("iteration %2d/%3d..\n",t+1,n_iter);
        printf("iteration %2d/%3d..\r",t+1,n_iter);
        fflush(stdout);
        /*
         * VB-E step
         */
        for(k = 0;k < n_class;k++){
            ln_lambda[k] = 0.0;
            for(d = 0;d < n_features;d++){
                ln_lambda[k] += gsl_sf_psi((nyu[k] + 1.0 - (double)(d+1)) / 2.0);
            }
            ln_lambda[k] += (double)n_features * log(2.0);
            ln_lambda[k] += log(determinant(W[k],n_features));
        }
        
        alpha_sum = 0.0;
        for(k = 0;k < n_class;k++){
            alpha_sum += alpha[k];
        }
        for(k = 0;k < n_class;k++){
            ln_pi[k] = gsl_sf_psi(alpha[k]) - gsl_sf_psi(alpha_sum);
        }
        
        for(n = 0;n < n_samples;n++){
            for(k = 0;k < n_class;k++){
                get_row(m,mk,n_class,n_features,k);
                vec_subtract(data[n],mk,xn_m,n_features);
                prod_mat_vec(W[k],xn_m,wk_xn_m,n_features,n_features);
                quad = dot_to_scalar(xn_m,wk_xn_m,n_features);
                ln_rho[k] = ln_pi[k] + ln_lambda[k]/2.0 - (double)n_features/(2.0*beta[k]) - nyu[k]*quad/2.0;
            }
            ln_rho_max = vec_max(ln_rho,n_class);
            r_sum = 0.0;
            for(k = 0;k < n_class;k++){
                r[n*n_class+k] = exp(ln_rho[k] - ln_rho_max);
                r_sum += r[n*n_class+k];
            }
            for(k = 0;k < n_class;k++){
                r[n*n_class+k] /= r_sum;
                if(r[n*n_class+k] < 1.0e-32)
                    r[n*n_class+k] = 1.0e-32;
            }
        }
        
        /*
         * VB-M step
         */
        for(k = 0;k < n_class;k++){
            N_k[k] = 0.0;
        }
        for(k = 0;k < n_class;k++){
            for(n = 0;n < n_samples;n++){
                N_k[k] += r[n*n_class+k];
            }
        }
        
        for(k = 0;k < n_class;k++){
            for(d = 0;d < n_features;d++){
                x_k[k*n_features+d] = 0.0;
                for(n = 0;n < n_samples;n++){
                    x_k[k*n_features+d] += r[n*n_class+k] * data[n][d];
                }
                x_k[k*n_features+d] /= N_k[k];
            }
        }
        
        for(k = 0;k < n_class;k++){
            get_row(x_k,xkk,n_class,n_features,k);
            for(d = 0;d < n_features*n_features;d++)
                S_k[k][d] = 0.0;
            for(n = 0;n < n_samples;n++){
                vec_subtract(data[n],xkk,xn_xk,n_features);
                dot_to_mat(xn_xk,xn_xk,xn_xk_mat,n_features);
                for(d = 0;d < n_features*n_features;d++)
                    xn_xk_mat[d] *= r[n*n_class+k];
                matrix_add(S_k[k],xn_xk_mat,S_k[k],n_features,n_features);
            }
            for(d = 0;d < n_features*n_features;d++)
                S_k[k][d] /= N_k[k];
        }
        
        for(k = 0;k < n_class;k++){
            alpha[k] = alpha0 + N_k[k];
            beta[k] = beta0 + N_k[k];
            nyu[k] = nyu0 + N_k[k];
            
            get_row(x_k,xkk,n_class,n_features,k);
            for(d = 0;d < n_features;d++){
                m[k*n_features+d] = (beta0*m0[k*n_features+d] + N_k[k]*xkk[d]) / beta[k];
            }
            
            for(d = 0;d < n_features*n_features;d++){
                Wk_inv[d] = 0.0;
            }
            for(d = 0;d < n_features*n_features;d++){
                NkSk[d] = N_k[k] * S_k[k][d];
            }
            get_row(m0,m0k,n_class,n_features,k);
            vec_subtract(xkk,m0k,xk_m0,n_features);
            dot_to_mat(xk_m0,xk_m0,xk_m0_mat,n_features);
            for(d = 0;d < n_features*n_features;d++)
                xk_m0_mat[d] *= (beta0*N_k[k] / (beta0+N_k[k]));
            matrix_add(xk_m0_mat,NkSk,Wk_inv,n_features,n_features);
            matrix_add(W0_inv,Wk_inv,Wk_inv,n_features,n_features);
            invert(Wk_inv,W[k],n_features);
        }
        lowerbound = lower_bound(alpha0,alpha,beta0,beta,nyu0,nyu,W0,W,r,n_samples,n_class,n_features);
        fprintf(lbp,"%.8f\n",lowerbound);
        //printf("\tLower_bound = %.8f\n",lowerbound);
    }
    printf("\n");
    fclose(lbp);
    // calculate posterior lambda
    for(k = 0;k < n_class;k++)
        for(d = 0;d < n_features*n_features;d++)
            l[k][d] = W[k][d]*nyu[k];
}
