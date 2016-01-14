/*
    learn.c
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "MT.h"
#include "learn.h"
#include "feature.h"
#include "imatrix.h"
#include "dmatrix.h"
#include "util.h"

void ldapf_learn(document *data, double alpha, double beta, int ndocs, int nclass, int nlex, int dlenmax, int nparticle, int ess, int rejuvenation, int **n_zw, double **theta){
    document *dp;
    int *n_z;
    int word_index;
    int word_num;
    double *left;
    double *right;
    double *p_z;
    double sum_p_z;
    double *cum_sum_p_z;
    double sampling;
    double **temp_theta;
    int ***particles_z;
    int **particles_n_mz;
    int ***particles_n_zw;
    int **particles_n_z;
    double *weights;
    double sum_weights;
    double sum_square_weights;
    double effective_sample_size;
    int z;
    int word_counter;
    int m, w, p, i, j, k, r;
    int random_word;
    int random_word_index;
    int random_word_num;
    int random_word_num_random;
    int previous_z;
    
    // count data length
    for(dp = data, ndocs = 0;(dp->len) != -1;dp++, ndocs++)
        ;
    // initialize buffers
    if((n_z = calloc(nclass,sizeof(int))) == NULL){
        fprintf(stderr,"ldapf_learn:: cannot allocate n_z.\n");
        return;
    }
    for(k = 0;k < nclass;k++)
        for(i = 0;i < nlex;i++)
            n_z[k] += n_zw[k][i];
    if((left = calloc(nclass,sizeof(double))) == NULL){
        fprintf(stderr,"ldapf_learn:: cannot allocate left.\n");
        return;
    }
    if((right = calloc(nclass,sizeof(double))) == NULL){
        fprintf(stderr,"ldapf_learn:: cannot allocate right.\n");
        return;
    }
    if((p_z = calloc(nclass,sizeof(double))) == NULL){
        fprintf(stderr,"ldapf_learn:: cannot allocate p_z.\n");
        return;
    }
    if((cum_sum_p_z = calloc((nclass+1),sizeof(double))) == NULL){
        fprintf(stderr,"ldapf_learn:: cannot allocate cum_sum_p_z.\n");
        return;
    }
    if((temp_theta = dmatrix(ndocs, nclass)) == NULL){
        fprintf(stderr,"ldapf_learn:: cannot allocate temp_theta.\n");
        exit(1);
    }
    
    printf("Number of new documents      = %d\n",ndocs);
    printf("Number of words              = %d\n",nlex);
    printf("Number of latent classes     = %d\n",nclass);
    printf("Number of particles          = %d\n",nparticle);
    printf("Rejuvenation steps           = %d\n",rejuvenation);
    printf("Effective sample size        = %d\n",ess);
    
    init_genrand(time(NULL));
    // each document
    for(dp = data, m = 0;(dp->len) != -1;dp++,m++){
        printf("document %2d / %3d..\n",m+1,ndocs);
        if(((particles_z = calloc(nparticle,sizeof(int **))) == NULL)
        || ((particles_n_mz = calloc(nparticle,sizeof(int *))) == NULL)
        || ((particles_n_zw = calloc(nparticle,sizeof(int **))) == NULL)
        || ((particles_n_z = calloc(nparticle,sizeof(int *))) == NULL)
        || ((weights = calloc(nparticle,sizeof(double))) == NULL)){
            fprintf(stderr,"cannot allocate particles.\n");
            return;
        }
        for(p = 0;p < nparticle;p++){
            if((particles_z[p] = calloc((dp->len),sizeof(int *))) == NULL){
                fprintf(stderr,"cannot allocate particles_z.\n");
                return;
            }
            if((particles_n_mz[p] = calloc(nclass,sizeof(int))) == NULL){
                fprintf(stderr,"cannot allocate particles_n_mz.\n");
                return;
            }
            if((particles_n_zw[p] = calloc(nclass,sizeof(int *))) == NULL){
                fprintf(stderr,"cannot allocate particles_n_zw.\n");
                return;
            }
            for(k = 0;k < nclass;k++){
                if((particles_n_zw[p][k] = calloc(nlex,sizeof(int))) == NULL){
                    fprintf(stderr,"cannot allocate particles_n_zw.\n");
                    return;
                }
            }
            for(k = 0;k < nclass;k++)
                for(i = 0;i < nlex;i++)
                    particles_n_zw[p][k][i] = n_zw[k][i];
            if((particles_n_z[p] = calloc(nclass,sizeof(int))) == NULL){
                fprintf(stderr,"cannot allocate particles_n_z.\n");
                return;
            }
            for(k = 0;k < nclass;k++)
                particles_n_z[p][k] = n_z[k];
        }
        // initialize weights
        for(p = 0;p < nparticle;p++)
            weights[p] = 1.0 / (double)nparticle;
        // each word in each document
        word_counter = 0;
        for(w = 0;w < (dp->len);w++){
            printf("\tword %2d / %3d..\r",w+1,(dp->len)); fflush(stdout);
            word_index = dp->id[w];
            word_num = dp->cnt[w];
            for(p = 0;p < nparticle;p++){
                if((particles_z[p][w] = calloc(word_num,sizeof(int))) == NULL){
                    fprintf(stderr,"cannot allocate particles_z.\n");
                    return;
                }
            }
            // each individual word
            for(i = 0;i < word_num;i++){
                // each particle
                for(p = 0;p < nparticle;p++){
                    // p_z left
                    for(k = 0;k < nclass;k++){
                        left[k] = (double)particles_n_mz[p][k] + alpha;
                        left[k] /= ((double)word_counter + (double)nclass * alpha);
                    }
                    // p_z right
                    for(k = 0;k < nclass;k++){
                        right[k] = (double)particles_n_zw[p][k][word_index] + beta;
                        right[k] /= ((double)particles_n_z[p][k] + (double)nlex * beta);
                    }
                    // conditional distribution p_z
                    sum_p_z = 0.0;
                    for(k = 0;k < nclass;k++){
                        p_z[k] = left[k] * right[k];
                        sum_p_z += left[k] * right[k];
                    }
                    for(k = 0;k < nclass;k++){
                        p_z[k] /= sum_p_z; // normalize to obtain probabilities
                    }
                    // random sampling from p_z
                    cum_sum_p_z[0] = 0.0;
                    for(k = 0;k < nclass;k++){
                        cum_sum_p_z[k+1] = cum_sum_p_z[k] + p_z[k];
                    }
                    sampling = genrand_real3();
                    for(k = 0;k < nclass;k++){
                        if((sampling >= cum_sum_p_z[k]) && (sampling < cum_sum_p_z[k+1])){
                            z = k;
                            break;
                        }
                    }
                    // update buffers
                    particles_n_z[p][z] += 1;
                    particles_n_mz[p][z] += 1;
                    particles_n_zw[p][z][word_index] += 1;
                    particles_z[p][w][i] = z;
                    // update weights
                    weights[p] *= sum_p_z;
                }
                // normalize weights
                sum_weights = 0.0;
                for(p = 0;p < nparticle;p++)
                    sum_weights += weights[p];
                for(p = 0;p < nparticle;p++)
                    weights[p] /= sum_weights;
                for(p = 0;p < nparticle;p++)
                    temp_theta[m][particles_z[p][w][i]] += weights[p];
                // check effective sample size
                sum_square_weights = 0.0;
                for(p = 0;p < nparticle;p++)
                    sum_square_weights += weights[p] * weights[p];
                effective_sample_size = 1.0 / sum_square_weights;
                // rejuvenation.
                if(word_counter > rejuvenation && effective_sample_size < (double)ess){
                    for(r = 0;r < rejuvenation;r++){
                        if(w == 0) break;
                        random_word = genrand_int32() % w; // random selection of previous words
                        random_word_index = dp->id[random_word];
                        random_word_num = dp->cnt[random_word];
                        if(random_word_num == 0) continue;
                        random_word_num_random = genrand_int32() % random_word_num;
                        for(p = 0;p < nparticle;p++){
                            previous_z = particles_z[p][random_word][random_word_num_random];
                            // remove assignment
                            particles_n_z[p][previous_z] -= 1;
                            particles_n_mz[p][previous_z] -= 1;
                            particles_n_zw[p][previous_z][random_word_index] -= 1;
                            //p_z left
                            for(k = 0;k < nclass;k++){
                                left[k] = (double)particles_n_mz[p][k] + alpha;
                                left[k] /= ((double)word_counter + (double)nclass * alpha);
                            }
                            // p_z right
                            for(k = 0;k < nclass;k++){
                                right[k] = (double)particles_n_zw[p][k][random_word_index] + beta;
                                right[k] /= ((double)particles_n_z[p][k] + (double)nlex * beta);
                            }
                            // conditional distribution p_z
                            sum_p_z = 0.0;
                            for(k = 0;k < nclass;k++){
                                p_z[k] = left[k] * right[k];
                                sum_p_z += left[k] * right[k];
                            }
                            for(k = 0;k < nclass;k++){
                                p_z[k] /= sum_p_z; // normalize to obtain probabilities
                            }
                            // random sampling from p_z
                            cum_sum_p_z[0] = 0.0;
                            for(k = 0;k < nclass;k++){
                                cum_sum_p_z[k+1] = cum_sum_p_z[k] + p_z[k];
                            }
                            sampling = genrand_real3();
                            for(k = 0;k < nclass;k++){
                                if((sampling >= cum_sum_p_z[k]) && (sampling < cum_sum_p_z[k+1])){
                                    z = k;
                                    break;
                                }
                            }
                            // update buffers
                            particles_n_z[p][z] += 1;
                            particles_n_mz[p][z] += 1;
                            particles_n_zw[p][z][random_word_index] += 1;
                            particles_z[p][random_word][random_word_num_random] = z;
                        }
                    }
                    for(p = 0;p < nparticle;p++)
                        weights[p] = 1.0 / (double)nparticle;
                }
                word_counter += 1;
            }
        }
        printf("\n");
        // free
        for(p = 0;p < nparticle;p++){
            for(w = 0;w < (dp->len);w++)
                free(particles_z[p][w]);
            free(particles_z[p]);
            free(particles_n_mz[p]);
            for(k = 0;k < nclass;k++){
                free(particles_n_zw[p][k]);
            }
            free(particles_n_zw[p]);
            free(particles_n_z[p]);
        }
        free(particles_z);
        free(particles_n_mz);
        free(particles_n_zw);
        free(particles_n_z);
        free(weights);
    }
    
    normalize_matrix_row(theta, temp_theta, ndocs, nclass);
    
    free(n_z);
    free(left);
    free(right);
    free(p_z);
    free(cum_sum_p_z);
    free_dmatrix(temp_theta, ndocs);
    
    return;
}
