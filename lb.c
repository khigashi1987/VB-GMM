/*
 * lb.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lb.h"
#include "matrix.h"

double ln_C0(double alpha0, int n_class);
double ln_C(double *alpha, int n_class);
double ln_B(double *W, double nyu, int n_features);

double lower_bound(double alpha0, double *alpha, double beta0, double *beta, double nyu0, double *nyu, double *W0, double **W, double *r, int n_samples, int n_class, int n_features){
    double lb = 0.0;
    double ln_beta = 0.0;
    double ln_BB = 0.0;
    double r_ln_r = 0.0;
    double result;
    int i,j;
    for(i = 0;i < n_class;i++){
        ln_beta += (log(beta0) - log(beta[i]));
        ln_BB += (ln_B(W0,nyu0,n_features) - ln_B(W[i],nyu[i],n_features));
    }
    for(i = 0;i < n_samples;i++){
        for(j = 0;j < n_class;j++){
            if(r[i*n_class+j] == 0.0){
                continue;
            }else{
                r_ln_r += r[i*n_class+j] * log(r[i*n_class+j]);
            }
        }
    }
    result = ln_C0(alpha0,n_class) - ln_C(alpha,n_class) + (double)n_features*ln_beta / 2.0 + ln_BB - r_ln_r - (double)n_features*(double)n_samples*log(2.0*M_PI)/2.0;
    /*
    printf("ln_C0 = %.8f\n",ln_C0(alpha0,n_class));
    printf("ln_C = %.8f\n",ln_C(alpha,n_class));
    for(i = 0;i < n_class;i++)
        printf("%.8f ",alpha[i]);
    printf("\n");
    printf("ln_BB = %.8f\n",ln_BB);
    printf("r_ln_r = %.8f\n",r_ln_r);
    */
    return result;
}

double ln_C0(double alpha0, int n_class){
    int i;
    double alpha_sum,result;
    alpha_sum = alpha0 * n_class;
    result = lgamma(alpha_sum);
    for(i = 0;i < n_class;i++){
        result -= lgamma(alpha0);
    }
    return result;
}

double ln_C(double *alpha, int n_class){
    int i;
    double alpha_sum = 0.0;
    double result;
    for(i = 0;i < n_class;i++){
        alpha_sum += alpha[i];
    }
    result = lgamma(alpha_sum);
    for(i = 0;i < n_class;i++){
        result -= lgamma(alpha[i]);
    }
    return result;
}

double ln_B(double *W, double nyu, int n_features){
    int i;
    double result;
    result = - nyu * log(determinant(W,n_features)) / 2.0 - nyu * (double)n_features * log(2.0) / 2.0 - (double)n_features * (double)(n_features - 1)  * log(M_PI) / 4.0;
    for(i = 0;i < n_features;i++){
        result -= lgamma((nyu + 1.0 - (double)(i+1)) / 2.0);
    }
    return result;
}
