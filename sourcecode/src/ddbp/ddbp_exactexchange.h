/***
 * @file    exactexchange_ddbp.h
 * @brief   This file contains the function declarations for Exact Exchange.
 *
 * @authors Xin Jing <xjing30@gatech.edu>
 *          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
 * 
 * Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech.
 */
 

#ifndef DDBP_EXACTEXCHANGE_H
#define DDBP_EXACTEXCHANGE_H

#include "isddft.h"

///////////////////////////////////////////////////////////////
///////////////  Initialization and finalization  /////////////
///////////////////////////////////////////////////////////////

/**
 * @brief   Initialization of all variables for exact exchange DDBP
 */
void init_exx_ddbp(SPARC_OBJ *pSPARC);

/**
 * @brief   Initialization of DDBP_info for exact exchange DDBP
 */
void init_exx_ddbp_info(SPARC_OBJ *pSPARC);

/**
 * @brief   free variables for exact exchange DDBP
 */
void free_exx_ddbp(SPARC_OBJ *pSPARC);


///////////////////////////////////////////////////////////////
//////////////  DDBP + Non-ACE kernel functions ///////////////
///////////////////////////////////////////////////////////////

/**
 * @brief   Compute projected Vx
 */
void Calculate_Vx_DDBP_Hamiltonian(SPARC_OBJ *pSPARC, int count, int kpt, int spn_i);

/**
 * @brief   Compute dense matrix for projected Vx without ACE
 */
void Calculate_Vx_DDBP_Hamiltonian_NonACE(SPARC_OBJ *pSPARC, int count, int kpt, int spn_i);

/**
 * @brief   Compute dense matrix for projected Vx without ACE
 */
void Calculate_H_DDBP_Vx_block(SPARC_OBJ *pSPARC);

/**
 * @brief   Compute VxVi locally
 */
void Calculate_VxVi(SPARC_OBJ *pSPARC);

/**
 * @brief   calculate all blocks of projected Vx in H_DDBP 
 */
void Calculate_VTVxV(DDBP_INFO *DDBP_info, double *VxViES);

/**
 * @brief   Projected exact exchange operator (V'VxV) times a vector. 
 */
void DDBP_Hamiltonian_ExactExchange_vectors_mult_kernel(
    DDBP_HAMILTONIAN *H_DDBP, const DDBP_ARRAY *X, 
    DDBP_ARRAY *HX, MPI_Comm comm);

/**
 * @brief   Projected exact exchange operator (V'VxV) times a vector without ACE
 */
void DDBP_Hamiltonian_ExactExchange_vectors_mult_kernel_NonACE(
    DDBP_HAMILTONIAN *H_DDBP, const DDBP_ARRAY *X, 
    DDBP_ARRAY *HX, MPI_Comm comm);

/**
 * @brief   evaluate exact exchange energy in case of DDBP
 */
void evaluate_exact_exchange_energy_DDBP(SPARC_OBJ *pSPARC);

/**
 * @brief   evaluate exact exchange energy in case of DDBP+nonACE
 */
void evaluate_exact_exchange_energy_DDBP_nonACE(SPARC_OBJ *pSPARC);


///////////////////////////////////////////////////////////////
//////////////////  DDBP + ACE kernel functions ///////////////
///////////////////////////////////////////////////////////////

/**
 * @brief   Create ACE operator in DDBP basis
 */
void ACE_operator_DDBP(SPARC_OBJ *pSPARC);

/**
 * @brief   Compute operator for projected Vx with ACE
 */
void Calculate_Vx_DDBP_Hamiltonian_ACE(SPARC_OBJ *pSPARC, int count, int kpt, int spn_i);

/**
 * @brief   Projected exact exchange operator (V'VxV) times a vector without ACE
 */
void DDBP_Hamiltonian_ExactExchange_vectors_mult_kernel_ACE(
    DDBP_HAMILTONIAN *H_DDBP, const DDBP_ARRAY *X, 
    DDBP_ARRAY *HX, MPI_Comm comm);

/**
 * @brief   Evaluate exact exchange energy in case of DDBP-ACE
 */
void evaluate_exact_exchange_energy_DDBP_ACE(SPARC_OBJ *pSPARC);


///////////////////////////////////////////////////////////////
//////////////////  Preparation and tools /////////////////////
///////////////////////////////////////////////////////////////

/**
 * @brief   update the usefock flag everywhere
 */
void update_usefock(SPARC_OBJ *pSPARC);

/**
 * @brief   Copy basis in other basiscomm for all elements
 *
 *          The default comm is DDBP_info->blacscomm
 */
void copy_ddbp_basis_xorb(SPARC_OBJ *pSPARC);

/**
 * @brief   Extract ACE oeprator for all elements
 */
void ExtractElementACE(SPARC_OBJ *pSPARC);

/**
 * @brief   Gather ddbp basis in other basiscomm
 *
 *          The default comm is DDBP_info->blacscomm
 */
void gather_ddbp_basis(DDBP_INFO *DDBP_info, 
    DDBP_ELEM *E_k, double *v, double *v_full, MPI_Comm comm);

/**
 * @brief   Gather element-wise vectors from all other elements 
 *          and rearrange into band-wise distribution
 *
 *          The default comm is DDBP_info->blacscomm
 */
void gather_ESvec_and_E2B(DDBP_INFO *DDBP_info, int ncol, double **vec_loc, double *vec_BS);

/**
 * @brief   gather full psi for each bandcomm
 */
void gather_psi_and_E2B(SPARC_OBJ *pSPARC);

/**
 * @brief   Copy occupations for DDBP
 */
void copy_occ_DDBP(SPARC_OBJ *pSPARC);

/**
 * @brief   rearrange array from element-wise to band-wise
 */
void ElementWise2BandWise_local(DDBP_INFO *DDBP_info, int ncol, double *x_in, double *x_out);

/**
 * @brief   rearrange local array from band-wise to element-wise
 */
void BandWise2ElementWise_local(DDBP_INFO *DDBP_info, int ncol, double *x_in, double *x_out);

/**
 * @brief   Find local Nd for each elements and DMvertices
 */
int find_local_Nde(DDBP_INFO *DDBP_info, int *coords, double *meshes, int *DMvertices);

/**
 * @brief   Transfer basis-wise array to element-wise array
 */
void BasisWise2ElementWise(DDBP_INFO *DDBP_info, double *x_in, double *x_out);

/**
 * @brief   Rearrange Column-Wise array with Multiple Elements
 */
void RearrangeColumnWiseMultipleElements(DDBP_INFO *DDBP_info, double *x_in, double *x_out);

/**
 * @brief   Poisson's solver for all pairs of local basis for DDBP
 */
void poissonSolveDDBP(SPARC_OBJ *pSPARC, 
    DDBP_ELEM *E_k, int ncol, double *v, double *v2, double *solg);

/**
 * @brief   Global vector (Ndg x nALB) times inner DDBP orbitals
 */
void GlobalvecTimesX(
    DDBP_INFO *DDBP_info, DDBP_ELEM *E_k, int k, 
    const double alpha, double *v, DDBP_ARRAY *xorb,
    const double beta, double *vX);

/**
 * @brief   Compute Vx times local basis
 */
void computeVxViDDBP(DDBP_INFO *DDBP_info, 
    double *solX, double *VxVi, int element_index, int basis_index);



/**
 * @brief   Compute Vx times local basis non-blocking way
 */
void computeVxViDDBP_nonblocking(DDBP_INFO *DDBP_info, 
    double *solX, double *VxVi, int element_index, int basis_index, MPI_Request *request);


/**
 * @brief   Get rank in basiscomm given basis index
 */
int basis_index_to_basiscomm_rank(int nALB, int nproc, int basis_index);

/**
 * @brief   Set up havx exchange information 
 */
void setup_havxX_DDBP_Array(DDBP_INFO *DDBP_info, DDBP_ARRAY *X, MPI_Comm comm);

/**
 * @brief   copy havxX_info
 */
void copy_havxX_info(const havxX_t *havxX_src, havxX_t *havxX_dest);

/**
 * @brief   free havxX_info
 */
void free_havxX_info(DDBP_ARRAY *X);

/**
 * @brief   Gather components from all elemnts for a DDBP array
 */
void DDBP_element_allelements_allgatherv(const DDBP_ARRAY *X, double *buff, MPI_Comm comm);

/**
 * @brief   Reduce and extract local VxViX
 */
void ReduceExtractVxViX(DDBP_INFO *DDBP_info, double *VxViX_g, double **VxViX, int opt);

/**
 * @brief   Gather vecotrs horizontally
 *
 * Note:    The default comm is blacscomm. 
 *          Generalization of gather_ddbp_basis function
 */
void gather_horizontally(
    int N1, int ncol_loc, int start_index, int ncol_tot, int np, 
    double *v, double *v_full, MPI_Comm comm);

/**
 * @brief   Change Veff to PBE in hybrid DDBP basis calculation
 * 
 * Note:    Always use Recover_hybrid_Veff together
 */
void Change_PBE_Veff(SPARC_OBJ *pSPARC);

/**
 * @brief   Recover hybrid Veff in hybrid DDBP basis calculation
 * 
 * Note:    Always use Change_PBE_Veff together
 */
void Recover_hybrid_Veff(SPARC_OBJ *pSPARC);

/**
 * @brief   Update the Vx part in H_DDBP after basis changed
 */
void Update_Vx_part_of_DDBP_Hamiltonian(SPARC_OBJ *pSPARC, int count, int kpt, int spn_i);


///////////////////////////////////////////////////////////////
//////////  Functions for parallelly solving poisson //////////
///////////////////////////////////////////////////////////////

/**
 * @brief   Solving poissons equations parallelly
 */
void solve_poisson_parallelly(SPARC_OBJ *pSPARC, DDBP_ELEM *E_k, double *sol);

/**
 * @brief   Distribute poissons solutions 
 */
void distribute_solutions(DDBP_INFO *DDBP_info, 
    DDBP_EXX_POIS_PARAL_D *pois_paral_d, double *sol, int n_basis);

/**
 * @brief   Solve poisson's equations construted by all pairs of basis
 */
void solve_allpair_poissons_equations(SPARC_OBJ *pSPARC, 
    DDBP_ELEM *E_k, DDBP_EXX_POIS_PARAL_G *pois_paral_g, 
    double *v, double *recvbuf, double *sol);

/**
 * @brief   Solve poisson's equations construted by distinct pairs of basis
 */
void solve_halfpair_poissons_equations(SPARC_OBJ *pSPARC, 
    DDBP_ELEM *E_k, double *v, double *sol);

/**
 * @brief   Gather half basis for solving poissons' equations
 */
void gather_half_basis(DDBP_EXX_POIS_PARAL_G *pois_paral_g, 
    int Nde, int ncol, double *v, double *recvbuf);

/**
 * @brief   Create communicator for parallelly solving poissons
 */
void Create_communicator_for_parallel_poisson(DDBP_INFO *DDBP_info);

/**
 * @brief   Free info for parallelly solving poissons
 */
void free_pois_paral_info(DDBP_INFO *DDBP_info);


void calculate_QQG(DDBP_INFO *DDBP_info, DDBP_ARRAY *X, double ***QQG);

void cwkp(double *A, double *B, int nrowA, int nrowB, int ncol, double *AkpB);

void fsp(double *A, double *B, int nrow, int ncolA, int ncolB, double *AfspB);

void gather_vertically(DDBP_INFO *DDBP_info, 
    int ncol, double **vec_loc, double *vec_full, int option);

void computeVxViN2method(DDBP_INFO *DDBP_info, 
    double *v_full_all, double *sol_es, double **QQGk, double *VxVi, int element_index, int basis_index, double **vfsps_save);
    
#endif