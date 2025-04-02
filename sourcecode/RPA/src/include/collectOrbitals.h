#ifndef COLLECTORBIT
#define COLLECTORBIT

#include "isddft.h"
#include "main.h"

/**
 * @brief Collect all psis, epsilons and occupations in different bandcomms, which are used for composing initial guess of Sternheimer eq.s
 *
 * @param pSPARC pointer to SPARC_OBJ
 * @param pRPA pointer to RPA_OBJ
 */
void collect_allXorb_allLambdas_gamma(SPARC_OBJ *pSPARC, RPA_OBJ *pRPA);

void collect_allXorb_allLambdas_kpt(SPARC_OBJ *pSPARC, RPA_OBJ *pRPA);

void send_recv_allXorb_allLambdas_kPq(SPARC_OBJ *pSPARC, RPA_OBJ *pRPA, int qptIndex);

#endif