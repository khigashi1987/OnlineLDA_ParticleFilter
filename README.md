# OnlineLDA_ParticleFilter
Topic prediction of new documents with pre-learned LDA model by Particle filter posterior estimation.

15.01.14 HIGASHI Koichi

The C implementation for online Latent Dirichlet Allocation by particle filter posterior estimation.
This tool is intended for use only in prediction of topic distributions for new documents and needs LDA model outputs learned in advance by LDA Batch-learning method.
The Mersenne twister are used for a pseudorandom number generator. C codes are taken from http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c

usage: ldapf [-P particles] [-E effective sample size] [-R rejuvenation steps] [-A alpha] [-B beta]  test model

-P (int) The number of particles.  
-E (int) Effective sample size threshold. If ESS falls below this threshold, "Rejuvenation" steps will be executed.  
-R (int) The number of "Rejuvenation" sampling.   
-A (float) Alpha parameter of pre-learned LDA model.  
-B (float) Beta parameter of pre-learned LDA model.  
test (file) New documents to be predicted.  
model (model prefix) Pre-learned model.  
