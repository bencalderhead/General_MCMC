#ifndef GMCMC_CHAIN_H
#define GMCMC_CHAIN_H

#include <gmcmc/model.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct __gmcmc_chain_st gmcmc_chain;

gmcmc_error gmcmc_chain_create(gmcmc_chain **, const gmcmc_model *);
gmcmc_error gmcmc_chain_create_copy(gmcmc_chain **, const gmcmc_chain *);
gmcmc_error gmcmc_chain_copy(gmcmc_chain *, const gmcmc_chain *);
void gmcmc_destroy(gmcmc_chain *);


#ifdef __cplusplus
}
#endif

#endif
