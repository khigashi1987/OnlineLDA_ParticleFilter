/*
    ldapf.h
    header file for LDA-Particle Filter
*/
#ifndef LDAPF_H
#define LDAPF_H
#include <stdlib.h>
#include "feature.h"

#define NPARTICLE_DEFAULT 100
#define ESS_DEFAULT 10
#define REJUVENATION_DEFAULT 10
#define ALPHA_DEFAULT 0.1
#define BETA_DEFAULT 0.1

void usage(void);

#endif
