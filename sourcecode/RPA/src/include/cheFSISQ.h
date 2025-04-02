#ifndef CHEFSIRPASQ
#define CHEFSIRPASQ

#include "isddft.h"
#include "main.h"

/**
 * @brief Make CheFSI on every nu^0.5 chi0 nu^0.5 (qpt, omega) operator to get its RPA energy term by SQ
 *
 * @param pSPARC pointer to SPARC_OBJ
 * @param pRPA pointer to RPA_OBJ
 * @param qptIndex the index of qpt of the current nu chi0 operator
 * @param omegaIndex the index of omega of the current nu chi0 operator
 */
void subspace_iteration_RPASQ(SPARC_OBJ *pSPARC, RPA_OBJ *pRPA, int qptIndex, int omegaIndex);

#endif