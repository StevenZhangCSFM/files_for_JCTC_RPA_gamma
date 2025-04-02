#ifndef CHEBFILTERRPA
#define CHEBFILTERRPA

#include "isddft.h"
#include "main.h"

/**
 * @brief Initialize the Delta V vectors to make CheFSI. The initialized vectors are pointed by pRPA->deltaVs_phi
 *
 * @param pSPARC pointer to SPARC_OBJ
 * @param pRPA pointer to RPA_OBJ
 */
void initialize_deltaVs(SPARC_OBJ *pSPARC, RPA_OBJ *pRPA);

/**
 * @brief Make CheFSI on every nu chi0(qpt, omega) operator to get its RPA energy term
 *
 * @param pSPARC pointer to SPARC_OBJ
 * @param pRPA pointer to RPA_OBJ
 * @param qptIndex the index of qpt of the current nu chi0 operator
 * @param omegaIndex the index of omega of the current nu chi0 operator
 */
void subspace_iteration_RPA(SPARC_OBJ *pSPARC, RPA_OBJ *pRPA, int qptIndex, int omegaIndex);

/**
 * @brief Compute the RPA energy term of the nu chi0(qpt, omega) operator from its eigenvalues
 *
 * @param RRnuChi0Eigs eigenvalues of the nu chi0(qpt, omega) operator
 * @param nuChi0Neig amount of solved eigenvalues
 * @param omegaMesh01 meshs in interval 0~1, mapped from 0~infty, Gauss-Legendre integral 
 * @param qptOmegaWeight the weight of the current nu chi0(qpt, omega), the weight of qpt multiplying the weight of omega
 */
double compute_ErpaTerm(double *RRnuChi0Eigs, int nuChi0Neig, double omegaMesh01, double qptOmegaWeight);

void initialize_timing_counter(RPA_OBJ *pRPA);

void finalize_timing_counter(RPA_OBJ *pRPA, int qptIndex, int omegaIndex, int flagNoDmcomm);

#endif