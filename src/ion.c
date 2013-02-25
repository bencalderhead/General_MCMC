typedef struct __ion_model_st * ion_model;

int ion_evaluate_mh(const ion_model, markov_chain);
int ion_proposal_mh(const ion_model, const markov_chain, double *, double *);
