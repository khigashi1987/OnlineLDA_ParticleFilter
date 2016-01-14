/*
    learn.h
*/
#ifndef LEARN_H
#define LEARN_H
#include <stdlib.h>
#include "feature.h"

#define RANDOM ((double)rand()/(double)RAND_MAX)

extern void ldapf_learn(document *data, double alpha, double beta, int ndocs, int nclass, int nlex, int dlenmax, int nparticle, int ess, int rejuvenation, int **n_zw, double **theta);
#endif
