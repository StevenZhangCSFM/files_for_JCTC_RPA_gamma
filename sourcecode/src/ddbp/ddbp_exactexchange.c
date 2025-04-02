/***
 * @file    ddbp_exactexchange.c
 * @brief   This file contains the functions for DDBP + Exact Exchange.
 *
 * @authors Xin Jing <xjing30@gatech.edu>
 *          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
 * 
 * Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech.
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include <complex.h>
#include <limits.h>
/** BLAS and LAPACK routines */
#ifdef USE_MKL
    // #define MKL_Complex16 double complex
    #include <mkl.h>
#else
    #include <cblas.h>
    #include <lapacke.h>
#endif
/** ScaLAPACK routines */
#ifdef USE_MKL
    #include "blacs.h"     // Cblacs_*
    #include <mkl_blacs.h>
    #include <mkl_pblas.h>
    #include <mkl_scalapack.h>
#endif
#ifdef USE_SCALAPACK
    #include "blacs.h"     // Cblacs_*
    #include "scalapack.h" // ScaLAPACK functions
#endif
#ifdef USE_FFTW
    #include <fftw3.h>
#endif

#include "exactExchange.h"
#include "exactExchangeInitialization.h"
#include "ddbp_exactexchange.h"
#include "ddbp_tools.h"
#include "ddbp_paral.h"
#include "linearAlgebra.h"
#include "ddbp_matvec.h"
#include "ddbp_basis.h"
#include "exchangeCorrelation.h"
#include "electronicGroundState.h"
#include "tools.h"

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))


#define TEMP_TOL (1e-14)


///////////////////////////////////////////////////////////////
///////////////  Initialization and finalization  /////////////
///////////////////////////////////////////////////////////////

/**
 * @brief   Initialization of all variables for exact exchange DDBP
 */
void init_exx_ddbp(SPARC_OBJ *pSPARC)
{
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;

///////////////////////////////////
// remove later
DDBP_info->H_DDBP.t_ddbp_Vx = 0;
DDBP_info->H_DDBP.c_ddbp_Vx = 0;
///////////////////////////////////

    int Nstates = pSPARC->Nstates;
    int Nkpts = pSPARC->Nkpts_kptcomm;
    int Nspin = pSPARC->Nspin_spincomm;
    int nelem = DDBP_info->n_elem_elemcomm;
    int Nband = DDBP_info->n_band_bandcomm; // #bands (KS orbitals)
    int ncol = DDBP_info->n_basis_basiscomm;
    int Nd_global = DDBP_info->Nd_global;

    DDBP_info->occ_outer = (double *) calloc(Nstates, sizeof(double));
    DDBP_info->psi_full = (double *) calloc(Nd_global*Nband, sizeof(double));
    assert(DDBP_info->occ_outer != NULL && DDBP_info->psi_full != NULL);

    if (pSPARC->ExxAcc == 0) {
        // set up DDBP Kohn-Sham orbitals from previous outer loop    
        DDBP_info->xorb_outer = malloc(Nspin * sizeof(*DDBP_info->xorb_outer));
        assert(DDBP_info->xorb_outer != NULL);
        for(int spn_i = 0; spn_i < Nspin; spn_i++) {
            DDBP_info->xorb_outer[spn_i] = malloc(Nkpts * sizeof(DDBP_ARRAY));
            assert(DDBP_info->xorb[spn_i] != NULL);
            for (int kpt = 0; kpt < Nkpts; kpt++) {
                DDBP_ARRAY *X_ks = &DDBP_info->xorb[spn_i][kpt];
                DDBP_ARRAY *Xouter_ks = &DDBP_info->xorb_outer[spn_i][kpt];
                duplicate_DDBP_Array_template(X_ks, Xouter_ks);
            }
        }

        DDBP_info->VxVi = (double *) calloc(Nd_global * ncol * nelem, sizeof(double));
        assert(DDBP_info->VxVi != NULL);
        
        // set up DDBP basis from previous outer loop
        for (int k = 0; k < nelem; k++) {
            DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
            E_k->v_outer = (double *) calloc(E_k->nd_d * ncol, sizeof(double));
            assert(E_k->v_outer != NULL);
        }
    
        // space for VxVi in H_DDBP_Ek
        for (int k = 0; k < nelem; k++) {
            DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
            DDBP_HAMILT_ERBLKS *H_DDBP_Ek = &E_k->H_DDBP_Ek;
            DDBP_EXX *Exx_DDBP = &H_DDBP_Ek->Exx_DDBP;
            Exx_DDBP->Ne_tot = DDBP_info->Ne_tot;
            Exx_DDBP->nALB = E_k->nALB;
            Exx_DDBP->h_ij = (double **) calloc(Exx_DDBP->Ne_tot, sizeof(double *));
            assert(Exx_DDBP->h_ij != NULL);
            for (int i = 0; i < Exx_DDBP->Ne_tot; i++) {
                Exx_DDBP->h_ij[i] = (double *) calloc(Exx_DDBP->nALB * Exx_DDBP->nALB, sizeof(double));
                assert(Exx_DDBP->h_ij[i] != NULL);
            }
        }
    } else {
        DDBP_info->Xi = (double **) calloc(nelem, sizeof(double *));
        assert(DDBP_info->Xi != NULL);
        for (int k = 0; k < nelem; k++) {
            DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
            // ACE operator in elements
            int Nde = E_k->nd_d;
            DDBP_info->Xi[k] = (double *) calloc(Nde * Nband, sizeof(double));
            assert(DDBP_info->Xi[k] != NULL);

            // allocate projected ACE operator
            DDBP_HAMILT_ERBLKS *H_DDBP_Ek = &E_k->H_DDBP_Ek;
            DDBP_EXX *Exx_DDBP = &H_DDBP_Ek->Exx_DDBP;
            Exx_DDBP->nALB = E_k->nALB;
            Exx_DDBP->ncol = Nstates;
            Exx_DDBP->xi_bar = (double *) calloc(Nstates*Exx_DDBP->nALB,sizeof(double));
            assert(Exx_DDBP->xi_bar != NULL);
            
            // ACE operator in extended elements
            int Nde_ex = E_k->nd_ex_d;
            SPARC_OBJ *ESPRC_k = E_k->ESPRC;
            ESPRC_k->Xi = (double *) calloc(Nde_ex * Nstates, sizeof(double));
            assert(ESPRC_k->Xi != NULL);
            assert(ESPRC_k->Nd_d_dmcomm == Nde_ex);
            
            // ACE operator in kptcomm_topo
            int Nde_ex_topo = E_k->nd_ex_d_topo;
            ESPRC_k->Xi_kptcomm_topo = (double *) calloc(Nde_ex_topo * Nstates, sizeof(double));
            assert(ESPRC_k->Xi_kptcomm_topo);
        }
    }

    // create n_basis info 
    MPI_Comm blacscomm = DDBP_info->blacscomm;
    int blacscomm_rank, blacscomm_size;
    MPI_Comm_rank(blacscomm, &blacscomm_rank);
    MPI_Comm_size(blacscomm, &blacscomm_size);

    int Nb_k = DDBP_info->nALB_tot / DDBP_info->Ne_tot;
    int nsubcomm = DDBP_info->npbasis;
    int basis_end_index;

    DDBP_info->basis_start_index_list = (int *) calloc(nsubcomm, sizeof(int));
    DDBP_info->n_basis_list = (int *) calloc(nsubcomm, sizeof(int));    
    assert(DDBP_info->basis_start_index_list != NULL &&
        DDBP_info->n_basis_list != NULL);
    
    int blacscomm_NZ_size = 0;    
    for (int i = 0; i < nsubcomm; i++) {
        assign_task_BLCYC_fashion(
            Nb_k, DDBP_info->npbasis, i,
            DDBP_info->basis_start_index_list + i, 
            &basis_end_index, DDBP_info->n_basis_list + i
        );
        if (*(DDBP_info->n_basis_list + i) > 0) blacscomm_NZ_size++;
    }
    DDBP_info->pois_paral_g.size = blacscomm_NZ_size;

#ifdef DEBUG
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank) printf("blacscomm size %d, proc with non-zero basis %d\n", blacscomm_size, blacscomm_NZ_size);
#endif


    // DDBP_info->psol_save_flag = 1;
    assert(DDBP_info->psol_save_flag == 1 || DDBP_info->psol_save_flag == 0);
    if (DDBP_info->psol_save_flag == 1) {
        // initialize communicator for parallelly solving poissons
        Create_communicator_for_parallel_poisson(DDBP_info);
        // allocate memory space for poissons' equations' solutions
        DDBP_info->psol = (double **) calloc(nelem, sizeof(double*));
        assert(DDBP_info->psol != NULL);
        for (int k = 0; k < nelem; k++) {
            DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
            int nALB = E_k->nALB;
            size_t psol_size = (size_t)Nd_global*ncol*nALB;
            DDBP_info->psol[k] = (double *) calloc(psol_size, sizeof(double));
            assert(DDBP_info->psol[k] != NULL);
        }

#ifdef DEBUG
    if (!rank) printf("\n== Memory space for saving Poisson's solutions ==\n");
    double memory_usage;
    char mem_str[32];
    int nALB = DDBP_info->nALB_tot / DDBP_info->Ne_tot;
    memory_usage = (double) Nd_global*nALB*nALB*DDBP_info->Ne_tot;
    formatBytes(memory_usage, 32, mem_str);
    if (!rank) printf("Total memory: %s\n", mem_str);
    memory_usage = (double) Nd_global*nALB*nALB;
    formatBytes(memory_usage, 32, mem_str);
    if (!rank) printf("Memory per element: %s\n", mem_str);
    memory_usage = (double) Nd_global*nALB*ncol*nelem;
    formatBytes(memory_usage, 32, mem_str);
    if (!rank) printf("Memory per processor: %s\n", mem_str);
#endif

    }

    if (DDBP_info->hybrid_memory_flag == 1) {
        // allocate memory spapce for vfsps
        DDBP_info->vfsps = (double****) malloc(sizeof(double) * nelem);
        for (int k = 0; k < nelem; k++) {
            DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
            int nALB = E_k->nALB;
            DDBP_info->vfsps[k] = (double ***) malloc(sizeof(double) * nALB);
            for (int i = 0; i < nALB; i++) {
                int Netot = DDBP_info->Ne_tot;
                DDBP_info->vfsps[k][i] = (double **) malloc(sizeof(double) * Netot);
                int nbasis = DDBP_info->n_basis_basiscomm;
                for (int n = 0; n < Netot; n++) {
                    int Nde = DDBP_info->Ndelist[n]; 
                    DDBP_info->vfsps[k][i][n] = (double *) malloc(sizeof(double) * Nde*nbasis*nALB);
                }
            }
        }
    }

#if defined(USE_MKL) || defined(USE_SCALAPACK)
    // set up descriptor for poisson's equaions' solution times X
    DDBP_info->desc_solX = (int **) malloc(nelem * sizeof(int*));
    assert(DDBP_info->desc_solX != NULL);
    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        SPARC_OBJ *ESPRC_k = E_k->ESPRC;
        DDBP_info->desc_solX[k] = malloc(9 * sizeof(int));
        assert(DDBP_info->desc_solX[k] != NULL);
        int ZERO = 0, mb, nb, llda, info;
        int Nrows = pSPARC->Nd; // global #rows
        int Ncols = Nstates; // global #cols
        mb = max(1, Nrows);
        nb = (Ncols - 1) / DDBP_info->npband + 1; // equal to ceil(Ncols/npband)
        // set up descriptor for storage of KS orbitals in ictxt_bandcomm (original)
        llda = max(1, Nrows);
        if (DDBP_info->bandcomm_index != -1 && DDBP_info->bandcomm != MPI_COMM_NULL) {
            descinit_(DDBP_info->desc_solX[k], &Nrows, &Ncols, &mb, &nb,
                &ZERO, &ZERO, &DDBP_info->ictxt_elemcomm, &llda, &info);
            assert(info == 0);
        } else {
            for (int i = 0; i < 9; i++) DDBP_info->desc_solX[k][i] = -1;
        }

        // set up ScaLAPACK descriptors for poisson's equaions' solution 
        Nrows = pSPARC->Nd; // global #rows
        Ncols = E_k->nALB; // global #cols
        mb = max(1, Nrows);
        nb = (E_k->nALB - 1) / DDBP_info->npbasis + 1; // equal to ceil(nALB/npbasis)
        // set up descriptor for storage of extended basis funcs in ictxt_blacscomm (original)
        llda = max(1, Nrows);
        if (DDBP_info->basiscomm_index != -1 && DDBP_info->dmcomm != MPI_COMM_NULL) {
            descinit_(E_k->desc_sol, &Nrows, &Ncols,
                  &mb, &nb, &ZERO, &ZERO, &ESPRC_k->ictxt_blacs, &llda, &info);
        } else {
            for (int i = 0; i < 9; i++) E_k->desc_sol[i] = -1;
        }
    }

    if (pSPARC->ExxAcc == 1) {
        int ZERO = 0, info;
        int mb, nb;
        mb = nb = Nstates;
        int nprow, npcol, myrow, mycol;
        // get coord of each process in block cyclic topology context
        Cblacs_gridinfo( DDBP_info->ictxt_elemcomm, &nprow, &npcol, &myrow, &mycol );

        if (DDBP_info->basiscomm_index != -1 && DDBP_info->dmcomm != MPI_COMM_NULL) {
            DDBP_info->nrows_M = numroc_( &Nstates, &mb, &myrow, &ZERO, &nprow);
            DDBP_info->ncols_M = numroc_( &Nstates, &nb, &mycol, &ZERO, &npcol);
            // checking only 1 proc save M
            int elemcomm_rank;
            MPI_Comm_rank(DDBP_info->elemcomm, &elemcomm_rank);
            if (!elemcomm_rank) 
                assert(DDBP_info->nrows_M == Nstates && DDBP_info->ncols_M == Nstates);
            else
                assert(DDBP_info->nrows_M * DDBP_info->ncols_M == 0);
            int llda = max(1, DDBP_info->nrows_M);
            descinit_(DDBP_info->desc_M, &Nstates, &Nstates, 
                &mb, &nb, &ZERO, &ZERO, &DDBP_info->ictxt_elemcomm, &llda, &info);
            assert(info == 0);

            // in case that col distribution is different from band parallelization
            DDBP_info->band_start_index_M = 0;
            for (int i = 0; i < mycol; i++) {
                int ncols = numroc_( &Nstates, &nb, &i, &ZERO, &npcol);
                DDBP_info->band_start_index_M += ncols;
            }
        } else {
            DDBP_info->nrows_M = DDBP_info->ncols_M = 0;
            for (int i = 0; i < 9; i++) DDBP_info->desc_M[i] = -1;
            DDBP_info->band_start_index_M = -1;
        }
    }
#endif

    pSPARC->Nkpts_shift = 1;
    // only for gamma point
    assert(pSPARC->Nkpts_shift == 1);

    allocate_singularity_removal_const(pSPARC);
    
    assert(pSPARC->ExxMethod == 0);
    if (pSPARC->ExxMethod == 0) {
        compute_pois_fft_const(pSPARC);
    }

    // initialize Eexx
    pSPARC->Eexx = 0;        
}

/**
 * @brief   Initialization of DDBP_info for exact exchange DDBP
 */
void init_exx_ddbp_info(SPARC_OBJ *pSPARC) 
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    // copy usefock flag and hybrid mixing coefficient
    DDBP_info->H_DDBP.usefock = pSPARC->usefock;
    DDBP_info->H_DDBP.exx_frac = pSPARC->exx_frac;
    DDBP_info->H_DDBP.Ne_tot = DDBP_info->Ne_tot;
    DDBP_info->H_DDBP.ExxAcc = pSPARC->ExxAcc;
    DDBP_info->H_DDBP.Nstates = pSPARC->Nstates;
    DDBP_info->H_DDBP.dV = pSPARC->dV;
    DDBP_info->psol_save_flag = pSPARC->DDBP_exx_psol_save_flag;
    DDBP_info->hybrid_memory_flag = pSPARC->DDBP_hybrid_memory_flag;
    if (pSPARC->usefock == 0) return;

    int Nstates = pSPARC->Nstates;
    int Nkpts = pSPARC->Nkpts_kptcomm;
    int Nspin = pSPARC->Nspin_spincomm;
    int nelem = DDBP_info->n_elem_elemcomm;
    DDBP_info->Nx_global = pSPARC->Nx;
    DDBP_info->Ny_global = pSPARC->Ny;
    DDBP_info->Nz_global = pSPARC->Nz;
    DDBP_info->Nd_global = pSPARC->Nd;
    DDBP_info->dV = pSPARC->dV;

    int nelem_global = DDBP_info->Ne_tot;
    DDBP_info->DMvert_list = (int **) calloc(nelem_global, sizeof(int *));
    assert(DDBP_info->DMvert_list != NULL);
    for (int index = 0; index < nelem_global; index++) {
        DDBP_info->DMvert_list[index] = (int *) calloc(6,sizeof(int));
        assert(DDBP_info->DMvert_list[index] != NULL);
    }

    // find out # of FD grid points in all elements
    int E_dims[3];
    double meshes[3];
    E_dims[0] = DDBP_info->Nex;
    E_dims[1] = DDBP_info->Ney;
    E_dims[2] = DDBP_info->Nez;
    meshes[0] = pSPARC->delta_x;
    meshes[1] = pSPARC->delta_y;
    meshes[2] = pSPARC->delta_z;
    DDBP_info->Ndelist = (int *) calloc(nelem_global, sizeof(int));
    DDBP_info->Ndedisp = (int *) calloc(nelem_global, sizeof(int));
    assert(DDBP_info->Ndelist != NULL && DDBP_info->Ndedisp != NULL);

    // go over all element index
    DDBP_info->Ndemax = 0;
    DDBP_info->Ndedisp[0] = 0;
    for (int index = 0; index < nelem_global; index++) {
        int coords[3];
        // find the global coordinates of the element
        DDBP_Index_Cart(E_dims, index, coords);
        DDBP_info->Ndelist[index] = find_local_Nde(DDBP_info, coords, meshes, DDBP_info->DMvert_list[index]);
        DDBP_info->Ndemax = max(DDBP_info->Ndemax, DDBP_info->Ndelist[index]);
        if (index != (nelem_global-1))
            DDBP_info->Ndedisp[index+1] = DDBP_info->Ndedisp[index] + DDBP_info->Ndelist[index];
    }

    // find list of number of elements 
    DDBP_info->Nelem_list = (int *) calloc(DDBP_info->npelem, sizeof(int));
    DDBP_info->elemstart_list = (int *) calloc(DDBP_info->npelem, sizeof(int));
    int elem_end_index;
    for (int i = 0; i < DDBP_info->npelem; i++) {
        assign_task_subcomm(
            nelem_global, DDBP_info->npelem, i,
            DDBP_info->elemstart_list+i, &elem_end_index, DDBP_info->Nelem_list+i
        );
    }
}

/**
 * @brief   free variables for exact exchange DDBP
 */
void free_exx_ddbp(SPARC_OBJ *pSPARC)
{
    free(pSPARC->k1_hf);
    free(pSPARC->k2_hf);
    free(pSPARC->k3_hf);
    free(pSPARC->kpthf_ind);
    free(pSPARC->kpthf_ind_red);
    free(pSPARC->kpthfred2kpthf);
    free(pSPARC->kpthf_pn);
    free(pSPARC->kpts_hf_red_list);
    
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    int Nkpts = pSPARC->Nkpts_kptcomm;
    int Nspin = pSPARC->Nspin_spincomm;
    int nelem = DDBP_info->n_elem_elemcomm;
    int nelem_global = DDBP_info->Ne_tot;

    for (int index = 0; index < nelem_global; index++) 
        free(DDBP_info->DMvert_list[index]);
    free(DDBP_info->DMvert_list);

    free(DDBP_info->Ndelist);
    free(DDBP_info->Ndedisp);
    free(DDBP_info->Nelem_list);
    free(DDBP_info->elemstart_list);
    free(DDBP_info->basis_start_index_list);
    free(DDBP_info->n_basis_list);

    free(DDBP_info->occ_outer);
    free(DDBP_info->psi_full);
    
    if (pSPARC->ExxAcc == 0) {
        for(int spn_i = 0; spn_i < Nspin; spn_i++) {
            for (int kpt = 0; kpt < Nkpts; kpt++) {
                DDBP_ARRAY *X_ks = &DDBP_info->xorb_outer[spn_i][kpt];
                delete_DDBP_Array(X_ks);
            }
            free(DDBP_info->xorb_outer[spn_i]);
        }
        free(DDBP_info->xorb_outer);
        free(DDBP_info->VxVi);
        for (int k = 0; k < nelem; k++) {
            DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
            free(E_k->v_outer);
        }

        for (int k = 0; k < nelem; k++) {
            DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
            DDBP_HAMILT_ERBLKS *H_DDBP_Ek = &E_k->H_DDBP_Ek;
            DDBP_EXX *Exx_DDBP = &H_DDBP_Ek->Exx_DDBP;
            for (int i = 0; i < Exx_DDBP->Ne_tot; i++) {
                free(Exx_DDBP->h_ij[i]);
            }
            free(Exx_DDBP->h_ij);
        }
    } else {
        for (int k = 0; k < nelem; k++) {
            DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
            free(DDBP_info->Xi[k]);
            DDBP_HAMILT_ERBLKS *H_DDBP_Ek = &E_k->H_DDBP_Ek;
            DDBP_EXX *Exx_DDBP = &H_DDBP_Ek->Exx_DDBP;
            free(Exx_DDBP->xi_bar);
            SPARC_OBJ *ESPRC_k = E_k->ESPRC;
            free(ESPRC_k->Xi);
            free(ESPRC_k->Xi_kptcomm_topo);
        }
        free(DDBP_info->Xi);
    }

    if (DDBP_info->psol_save_flag == 1 && 0) {
        for (int k = 0; k < nelem; k++) 
            free(DDBP_info->psol[k]);
        free(DDBP_info->psol);
        free_pois_paral_info(DDBP_info);
    }

    if (DDBP_info->hybrid_memory_flag == 1) {                
        for (int k = 0; k < nelem; k++) {            
            DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
            int nALB = E_k->nALB;
            for (int i = 0; i < nALB; i++) {
                int Netot = DDBP_info->Ne_tot;                
                for (int n = 0; n < Netot; n++) {
                    free(DDBP_info->vfsps[k][i][n]);
                }
                free(DDBP_info->vfsps[k][i]);
            }
            free(DDBP_info->vfsps[k]);
        }
        free(DDBP_info->vfsps);
    }

#if defined(USE_MKL) || defined(USE_SCALAPACK)
    for (int k = 0; k < nelem; k++) {
        free(DDBP_info->desc_solX[k]);
    }
    free(DDBP_info->desc_solX);
#endif

    if (pSPARC->isGammaPoint == 1) {
        free(pSPARC->pois_const);
    }
}


///////////////////////////////////////////////////////////////
//////////////  DDBP + Non-ACE kernel functions ///////////////
///////////////////////////////////////////////////////////////

/**
 * @brief   Compute projected Vx
 */
void Calculate_Vx_DDBP_Hamiltonian(SPARC_OBJ *pSPARC, int count, int kpt, int spn_i)
{
#ifdef DEBUG
    int rank;
    double t1, t2;
    t1 = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank) printf("\n== Construting H_DDBP of Vx operator == \n");
#endif

    if (pSPARC->ExxAcc == 0)
        Calculate_Vx_DDBP_Hamiltonian_NonACE(pSPARC, count, kpt, spn_i);
    else
        Calculate_Vx_DDBP_Hamiltonian_ACE(pSPARC, count, kpt, spn_i);

#ifdef DEBUG
    t2 = MPI_Wtime();
    if (!rank) printf("== Total time for Construting H_DDBP of Vx operator: %.3f ms\n", (t2-t1)*1e3);
#endif
}


/**
 * @brief   Compute dense matrix for projected Vx without ACE
 */
void Calculate_Vx_DDBP_Hamiltonian_NonACE(SPARC_OBJ *pSPARC, int count, int kpt, int spn_i)
{
    // solve poisson's equations and find local VxVi
    Calculate_VxVi(pSPARC);
    // calculate the dense part of H_DDBP_Vx
    Calculate_H_DDBP_Vx_block(pSPARC);
}

/**
 * @brief   Compute dense matrix for projected Vx without ACE
 */
void Calculate_H_DDBP_Vx_block(SPARC_OBJ *pSPARC)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    int nelem = DDBP_info->n_elem_elemcomm;
    int Ndg = DDBP_info->Nd_global;
    int nbasis = DDBP_info->n_basis_basiscomm;    
    double t1, t2;
    
    t1 = MPI_Wtime();
    double *VxVi_ = (double *) calloc(Ndg * nbasis * nelem, sizeof(double));
    assert(VxVi_ != NULL);
    // rearrange the x_in for communication
    RearrangeColumnWiseMultipleElements(DDBP_info, DDBP_info->VxVi, VxVi_);    

    int Nde_tot_loc = 0;
    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        Nde_tot_loc += E_k->nd_d;
    }
    int Ne_tot = DDBP_info->Ne_tot;
    double *VxViES = (double *) calloc(Nde_tot_loc*nbasis*Ne_tot, sizeof(double));
    assert(VxViES != NULL);
    // communication for changing to element-wise
    BasisWise2ElementWise(DDBP_info, VxVi_, VxViES);
    free(VxVi_);
    t2 = MPI_Wtime();
#ifdef DEBUG
    if (!rank) printf("== Time for rearrange and communication of VxVi    : %.3f ms\n", (t2-t1)*1e3);
#endif

    // compute V^TVxV
    t1 = MPI_Wtime();
    Calculate_VTVxV(DDBP_info, VxViES);
    t2 = MPI_Wtime();
#ifdef DEBUG
    if (!rank) printf("== Time for calculating VTVxV                      : %.3f ms\n", (t2-t1)*1e3);
#endif
    free(VxViES);
}


/**
 * @brief   Compute VxVi locally
 */
void Calculate_VxVi(SPARC_OBJ *pSPARC)
{
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    if (DDBP_info->elemcomm_index == -1 || DDBP_info->basiscomm_index == -1 ||
        DDBP_info->dmcomm == MPI_COMM_NULL) {
        return;
    }
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm blacscomm = DDBP_info->blacscomm;
    int nelem = DDBP_info->n_elem_elemcomm;
    int Ndg = DDBP_info->Nd_global;
    int nbasis = DDBP_info->n_basis_basiscomm;
    int Nband = DDBP_info->n_band_bandcomm;
    int Nstates = pSPARC->Nstates;
    double t1, t2, t_gather, t_pois, t_poisX, t_VxVi;
    t_gather = t_pois = t_poisX = t_VxVi = 0.0;

    double *solX = (double *) calloc(Ndg * Nband, sizeof(double));
    assert(solX != NULL);

    // 3 conditions for skipping solving poissons' equations. 
    int PoissonSkipFlag = (pSPARC->fock_err < pSPARC->DDBP_fock_tol_updatebasis)
                          && (DDBP_info->psol_save_flag == 1) && (DDBP_info->skip_pois_flag == 1);
    assert(PoissonSkipFlag == 1 || PoissonSkipFlag == 0);
    // first outer loop after fock_err drops below tol, we still need to update Vx, i.e. solve pois
    if (pSPARC->fock_err < pSPARC->DDBP_fock_tol_updatebasis) DDBP_info->skip_pois_flag = 1;

    double *VxVi = DDBP_info->VxVi;
    for (int k = 0; k < nelem; k++) {
        t1 = MPI_Wtime();
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        int Nde = E_k->nd_d;
        int nALB = E_k->nALB;
        double *v_full = (double *) calloc(Nde * nALB, sizeof(double));
        assert(v_full != NULL);
        // gather current DDBP basis
        // gather_ddbp_basis(DDBP_info, E_k, E_k->v, v_full, blacscomm);
        gather_horizontally(
            Nde, nbasis, DDBP_info->basis_start_index, nALB, 
            DDBP_info->npbasis, E_k->v, v_full, blacscomm);
        t2 = MPI_Wtime();
        t_gather += (t2 - t1);

        // creating requests for non-blocking reducing
        MPI_Request *request_reduce = (MPI_Request *) calloc(nALB, sizeof(MPI_Request));
        assert(request_reduce != NULL);

// TODO: Actually we could multiply all poisson's solution with X together
// TODO: by stacking solutions vertically and it should be faster. But the
// TODO: memory required is nALB times more than it needs now. 

        for (int i = 0; i < nALB; i++) {
            double *sol;
            if (DDBP_info->psol_save_flag == 1)
                sol = DDBP_info->psol[k] + i*Ndg*nbasis;
            else
                sol = (double *) calloc(Ndg * nbasis, sizeof(double));
            assert(sol != NULL);

            if (PoissonSkipFlag == 0) {
                t1 = MPI_Wtime();
                double *v_full_i = v_full + i * Nde;
                // solve poisson's equation of v_full_i * v_outer 
                poissonSolveDDBP(pSPARC, E_k, nbasis, E_k->v_outer, v_full_i, sol);
                t2 = MPI_Wtime();
                t_pois += (t2 - t1);
            }

            t1 = MPI_Wtime();
            // solution times X to get poisson's equation of v_full_i times psi_outer            
            DDBP_ARRAY *xorb_outer = &DDBP_info->xorb_outer[0][0];
            GlobalvecTimesX(DDBP_info, E_k, k, 1.0, sol, xorb_outer, 0.0, solX);
            t2 = MPI_Wtime();
            t_poisX += (t2 - t1);
            if (DDBP_info->psol_save_flag == 0) free(sol);
            
            t1 = MPI_Wtime();
            // compute VxVi and put into correct places
            // computeVxViDDBP(DDBP_info, solX, VxVi + k*Ndg*nbasis, k, i);
            computeVxViDDBP_nonblocking(DDBP_info, solX, VxVi + k*Ndg*nbasis, k, i, request_reduce+i);
            t2 = MPI_Wtime();
            t_VxVi += (t2 - t1);            
        }
        free(v_full);

        // for non-blocking reduce
        t1 = MPI_Wtime();
        MPI_Waitall(nALB, request_reduce, MPI_STATUS_IGNORE);
        free(request_reduce);
        t2 = MPI_Wtime();
        t_VxVi += (t2 - t1);
    }
    free(solX);

#ifdef DEBUG
    if (!rank) printf("== Time for gathering ddbp basis accross blacscomm : %.3f ms\n", t_gather*1e3);
    if (PoissonSkipFlag == 0) {
        if (!rank) printf("== Time for solving poisson's equations            : %.3f ms\n", t_pois*1e3);
    } else {
        if (!rank) printf("== Poisson's equations are not solved and solution retrived\n");
    }
    if (!rank) printf("== Time for poisson's equation times DDBP orbitals : %.3f ms\n", t_poisX*1e3);
    if (!rank) printf("== Time for calculating VxVi and reduce            : %.3f ms\n", t_VxVi*1e3);
#endif    
}


/**
 * @brief   calculate all blocks of projected Vx in H_DDBP 
 */
void Calculate_VTVxV(DDBP_INFO *DDBP_info, double *VxViES)
{
    MPI_Comm comm = DDBP_info->blacscomm;
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    int nelem = DDBP_info->n_elem_elemcomm;
    int npelem = DDBP_info->npelem;
    int nbasis = DDBP_info->n_basis_basiscomm;

    double *VxVi = VxViES;
    for (int r = 0; r < npelem; r++) {
        int nelem_src = DDBP_info->Nelem_list[r];
        int elem_start_index = DDBP_info->elemstart_list[r];
        
        for (int j = 0; j < nelem; j++) {
            DDBP_ELEM *E_k = &DDBP_info->elem_list[j];
            int Nde = E_k->nd_d;
            int nALB = E_k->nALB;
            double *v = E_k->v;
            DDBP_HAMILT_ERBLKS *H_DDBP_Ek = &E_k->H_DDBP_Ek;

            for (int i = 0; i < nelem_src; i++) {
                int global_i = i + elem_start_index;
                double *h_ij = H_DDBP_Ek->Exx_DDBP.h_ij[global_i];
            
            #if defined(USE_MKL) || defined(USE_SCALAPACK)
                int deschij[9];
                int ZERO = 0, info;
                int mb = max(1, nALB);
                int nb = max(1, nALB);
                int llda = max(1, nALB);
                double beta = 0.0;

                descinit_(deschij, &nALB, &nALB, &mb, &nb, &ZERO, &ZERO, &DDBP_info->ictxt_blacs, &llda, &info);

                pdgemm_subcomm("T", "N", nALB, nALB, Nde, 1.0, v, E_k->desc_v, VxVi, E_k->desc_v,
                            beta, h_ij, deschij, DDBP_info->blacscomm, best_max_nproc(Nde, nALB, "pdgemm"));
            #else
                assert(0);
            #endif

                if (DDBP_info->npdm > 1) { // in the current implementation, this won't happen
                    MPI_Allreduce(MPI_IN_PLACE, h_ij, nALB*nALB, MPI_DOUBLE, MPI_SUM, DDBP_info->dmcomm);
                }

                if (DDBP_info->npbasis > 1)
                    MPI_Bcast(h_ij, nALB*nALB, MPI_DOUBLE, 0, DDBP_info->blacscomm);

                VxVi += Nde * nbasis;
            }
        }
    }
}


/**
 * @brief   Projected exact exchange operator (V'VxV) times a vector. 
 */
void DDBP_Hamiltonian_ExactExchange_vectors_mult_kernel(
    DDBP_HAMILTONIAN *H_DDBP, const DDBP_ARRAY *X, 
    DDBP_ARRAY *HX, MPI_Comm comm)
{    
    double t1, t2;
    t1 = MPI_Wtime();
    if (H_DDBP->ExxAcc == 0) 
        DDBP_Hamiltonian_ExactExchange_vectors_mult_kernel_NonACE(H_DDBP, X, HX, comm);
    else
        DDBP_Hamiltonian_ExactExchange_vectors_mult_kernel_ACE(H_DDBP, X, HX, comm);    
    t2 = MPI_Wtime();
    H_DDBP->t_ddbp_Vx += (t2 - t1);
    H_DDBP->c_ddbp_Vx ++;
}

/**
 * @brief   Projected exact exchange operator (V'VxV) times a vector without ACE
 */
void DDBP_Hamiltonian_ExactExchange_vectors_mult_kernel_NonACE(
    DDBP_HAMILTONIAN *H_DDBP, const DDBP_ARRAY *X, 
    DDBP_ARRAY *HX, MPI_Comm comm)
{
    const havxX_t *havxX_info = &X->havxX_info;
    // in case that LANCZOS_DDBP proc with havxX_flag = 0 will enter here.
    if (havxX_info->havxX_flag == 0) return;
    assert(havxX_info->havxX_flag == 1);
    int ncol = X->ncol;
    int nrows_total = havxX_info->nrows_total;
    double *buff = (double *) calloc(ncol*nrows_total, sizeof(double));
    assert(buff != NULL);

    DDBP_element_allelements_allgatherv(X, buff, comm);

    int Ne_tot = H_DDBP->Ne_tot;
    int blksz = havxX_info->blksz;
    int nelem = X->nelem;
    // warning: only correct when nrows are the same 
    // for all elements (in current implementation)
    int nALB = X->nrows[0];
    double alpha = H_DDBP->exx_frac;
    double beta = 1.0;

    // TODO: Assemble the matrices and change to use only 1 dgemm
    for (int i = 0; i < nelem; i++) {
        double *hx_i = HX->array[i];
        for (int j = 0; j < Ne_tot; j++) {
            double *xj = buff + j*blksz;
            double *h_ij = H_DDBP->H_DDBP_Ek_list[i]->Exx_DDBP.h_ij[j];
            adaptive_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nALB, ncol, nALB,
                alpha, h_ij, nALB, xj, nALB, beta, hx_i, nALB);
        }
    }    

    free(buff);
}


/**
 * @brief   evaluate exact exchange energy in case of DDBP
 */
void evaluate_exact_exchange_energy_DDBP(SPARC_OBJ *pSPARC)
{
    double t1, t2;
    t1 = MPI_Wtime();
    if (pSPARC->ExxAcc == 0)
        evaluate_exact_exchange_energy_DDBP_nonACE(pSPARC);
    else
        evaluate_exact_exchange_energy_DDBP_ACE(pSPARC);
    t2 = MPI_Wtime();

#ifdef DEBUG
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(!rank) {
        printf("\n== Time for calculating exact exchange enrgy : %.3f ms\n", (t2-t1)*1e3);
        printf("\n Exact exchange energy %.6f.\n", pSPARC->Eexx);
    }
#endif  
}


/**
 * @brief   evaluate exact exchange energy in case of DDBP+nonACE
 */
void evaluate_exact_exchange_energy_DDBP_nonACE(SPARC_OBJ *pSPARC)
{
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    // if (DDBP_info->elemcomm_index == -1) return;
    if (DDBP_info->elemcomm_index == -1 || DDBP_info->basiscomm_index == -1) return;
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm blacscomm = DDBP_info->blacscomm;
    int nelem = DDBP_info->n_elem_elemcomm;
    int Ndg = DDBP_info->Nd_global;
    int nbasis = DDBP_info->n_basis_basiscomm;
    int Nband = DDBP_info->n_band_bandcomm;
    int Nstates = pSPARC->Nstates;
    int spn_i = 0, kpt = 0;
    double t1, t2, t_gather, t_pois, t_poisX, t_VxVi;
    t_gather = t_pois = t_poisX = t_VxVi = 0.0;

    double *VxViX_g = (double *) calloc(Ndg * Nband, sizeof(double));    
    assert(VxViX_g != NULL);

    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        DDBP_ARRAY *xorb = &DDBP_info->xorb[spn_i][kpt];
        double *VxVi = DDBP_info->VxVi + k*Ndg*nbasis;
        GlobalvecTimesX(DDBP_info, E_k, k, 1.0, VxVi, xorb, 1.0, VxViX_g);
    }
    
    double **VxViX = (double **) calloc(nelem, sizeof(double *));
    assert(VxViX != NULL);
    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        int Nde = E_k->nd_d;
        VxViX[k] = (double *) calloc(Nde*Nband, sizeof(double));
        assert(VxViX[k] != NULL);
    }
    ReduceExtractVxViX(DDBP_info, VxViX_g, VxViX, 1);
    free(VxViX_g);
    
    // calculate exact exchange energy
    double Eexx = 0;
    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        int Nde = E_k->nd_d;
        double *psi_k = DDBP_info->psi[spn_i][kpt][k];
        double *VxViX_k = VxViX[k];
        for (int n = 0; n < Nband; n++) {
            double occ = pSPARC->occ[n + DDBP_info->band_start_index];
            double *psi_kn = psi_k + n*Nde;
            double *VxViX_kn = VxViX_k + n*Nde;
            for (int i = 0; i < Nde; i++) {
                Eexx += occ * psi_kn[i] * VxViX_kn[i];
            }
        }
    }

    MPI_Comm bandcomm = DDBP_info->bandcomm;
    int size_bandcomm, rank_bandcomm;
    MPI_Comm_size(bandcomm, &size_bandcomm);
    MPI_Comm_rank(bandcomm, &rank_bandcomm);
    if (size_bandcomm > 1) {
        MPI_Allreduce(MPI_IN_PLACE, &Eexx, 1, MPI_DOUBLE, MPI_SUM, bandcomm);
    }
    
    int size_blacscomm, rank_blacscomm;
    MPI_Comm_size(blacscomm, &size_blacscomm);
    MPI_Comm_rank(blacscomm, &rank_blacscomm);
    if (size_blacscomm > 1) {        
        MPI_Allreduce(MPI_IN_PLACE, &Eexx, 1, MPI_DOUBLE, MPI_SUM, blacscomm);
    }    
    pSPARC->Eexx = Eexx * pSPARC->exx_frac;
    
    for (int k = 0; k < nelem; k++) 
        free(VxViX[k]);
    free(VxViX);
}



///////////////////////////////////////////////////////////////
//////////////////  DDBP + ACE kernel functions ///////////////
///////////////////////////////////////////////////////////////

/**
 * @brief   Create ACE operator in DDBP basis
 */
void ACE_operator_DDBP_old(SPARC_OBJ *pSPARC)
{
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    if (DDBP_info->elemcomm_index == -1 || DDBP_info->basiscomm_index == -1 ||
        DDBP_info->dmcomm == MPI_COMM_NULL) {
        return;
    }
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm blacscomm = DDBP_info->blacscomm;
    int nelem = DDBP_info->n_elem_elemcomm;
    int Ndg = DDBP_info->Nd_global;
    int nbasis = DDBP_info->n_basis_basiscomm;
    int Nband = DDBP_info->n_band_bandcomm;
    int Nstates = pSPARC->Nstates;
    int spn_i = 0, kpt = 0;
    double t1, t2, t_gather, t_pois, t_poisX, t_VxVi, t_VxViX, t_VxViXReduce, t_M, t_Choltrsm;
    t_gather = t_pois = t_poisX = t_VxVi = t_VxViX = t_VxViXReduce = t_M = t_Choltrsm = 0.0;

    double *VxViX_g = (double *) calloc(Ndg * Nband, sizeof(double));    
    assert(VxViX_g != NULL);

    // 3 conditions for skipping solving poissons' equations. 
    int PoissonSkipFlag = (pSPARC->fock_err < pSPARC->DDBP_fock_tol_updatebasis)
                          && (DDBP_info->psol_save_flag == 1) && (DDBP_info->skip_pois_flag == 1);
    assert(PoissonSkipFlag == 1 || PoissonSkipFlag == 0);
    // first outer loop after fock_err drops below tol, we still need to update Vx, i.e. solve pois
    if (pSPARC->fock_err < pSPARC->DDBP_fock_tol_updatebasis) DDBP_info->skip_pois_flag = 1;

    // if (!rank) printf("fock_err %e, DDBP_fock_tol_updatebasis %e, psol_save_flag %d, skip_pois_flag %d, PoissonSkipFlag %d\n",
        // pSPARC->fock_err, pSPARC->DDBP_fock_tol_updatebasis, DDBP_info->psol_save_flag, DDBP_info->skip_pois_flag, PoissonSkipFlag);

    for (int k = 0; k < nelem; k++) {
        t1 = MPI_Wtime();
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        int Nde = E_k->nd_d;
        int nALB = E_k->nALB;
        double *v_full = (double *) calloc(Nde * nALB, sizeof(double));
        assert(v_full != NULL);
        // gather current DDBP basis
        // gather_ddbp_basis(DDBP_info, E_k, E_k->v, v_full, blacscomm);
        gather_horizontally(
            Nde, nbasis, DDBP_info->basis_start_index, nALB, 
            DDBP_info->npbasis, E_k->v, v_full, blacscomm);
        t2 = MPI_Wtime();
        t_gather += (t2 - t1);

        double *VxVi = (double *) calloc(Ndg * nbasis, sizeof(double));
        assert(VxVi != NULL);
        // creating requests for non-blocking reducing
        MPI_Request *request_reduce = (MPI_Request *) calloc(nALB, sizeof(MPI_Request));
        assert(request_reduce != NULL);

        // The storage of total solutions is O(Nde*Ne*nALB*nALB*Ne) = O(Nde*Ne^2*nALB^2)
        // For large systems e.g. 100 elements, each element needs around 10Gb for solutions
        if (DDBP_info->psol_save_flag == 1) {
            if (PoissonSkipFlag == 0) {
                t1 = MPI_Wtime();
                solve_poisson_parallelly(pSPARC, E_k, DDBP_info->psol[k]);
                t2 = MPI_Wtime();
                t_pois += (t2 - t1);
            }
            
            for (int i = 0; i < nALB; i++) {
                double *sol = DDBP_info->psol[k] + i*Ndg*nbasis;
                double *solX = (double *) calloc(Ndg * Nband, sizeof(double));
                assert(solX != NULL);
                t1 = MPI_Wtime();
                // solution times X to get poisson's equation of v_full_i times psi_outer
                DDBP_ARRAY *xorb = &DDBP_info->xorb[spn_i][kpt];
                GlobalvecTimesX(DDBP_info, E_k, k, 1.0, sol, xorb, 0.0, solX);
                t2 = MPI_Wtime();
                t_poisX += (t2 - t1);

                t1 = MPI_Wtime();
                // compute VxVi and put into correct places
                // printf("rank %d ind %d enter\n", rank, i);
                // computeVxViDDBP(DDBP_info, solX, VxVi, k, i);
                computeVxViDDBP_nonblocking(DDBP_info, solX, VxVi, k, i, request_reduce + i);
                t2 = MPI_Wtime();
                t_VxVi += (t2 - t1);
                free(solX);
            }
        } else {
            for (int i = 0; i < nALB; i++) {
                double *sol = (double *) calloc(Ndg * nbasis, sizeof(double));
                assert(sol != NULL);

                t1 = MPI_Wtime();
                double *v_full_i = v_full + i * Nde;
                // solve poisson's equation of v_full_i * v
                poissonSolveDDBP(pSPARC, E_k, nbasis, E_k->v, v_full_i, sol);
                t2 = MPI_Wtime();
                t_pois += (t2 - t1);

                double *solX = (double *) calloc(Ndg * Nband, sizeof(double));
                assert(solX != NULL);
                t1 = MPI_Wtime();
                // solution times X to get poisson's equation of v_full_i times psi_outer
                DDBP_ARRAY *xorb = &DDBP_info->xorb[spn_i][kpt];
                GlobalvecTimesX(DDBP_info, E_k, k, 1.0, sol, xorb, 0.0, solX);
                t2 = MPI_Wtime();
                t_poisX += (t2 - t1);
                free(sol);

                t1 = MPI_Wtime();
                // compute VxVi and put into correct places
                // computeVxViDDBP(DDBP_info, solX, VxVi, k, i);
                computeVxViDDBP_nonblocking(DDBP_info, solX, VxVi, k, i, request_reduce + i);
                t2 = MPI_Wtime();
                t_VxVi += (t2 - t1);
                free(solX);
            }
        }
        free(v_full);

        // for non-blocking reduce
        t1 = MPI_Wtime();
        MPI_Waitall(nALB, request_reduce, MPI_STATUS_IGNORE);
        free(request_reduce);
        t2 = MPI_Wtime();
        t_VxVi += (t2 - t1);

        DDBP_ARRAY *xorb = &DDBP_info->xorb[spn_i][kpt];
        t1 = MPI_Wtime();
        GlobalvecTimesX(DDBP_info, E_k, k, 1.0, VxVi, xorb, 1.0, VxViX_g);
        t2 = MPI_Wtime();
        t_VxViX += (t2 - t1);
        free(VxVi);
    }
    
    t1 = MPI_Wtime();
    ReduceExtractVxViX(DDBP_info, VxViX_g, NULL, 2);
    free(VxViX_g);
    t2 = MPI_Wtime();
    t_VxViXReduce += (t2 - t1);

    t1 = MPI_Wtime();
    int nrows_M, ncols_M;
    nrows_M = DDBP_info->nrows_M;
    ncols_M = DDBP_info->ncols_M;
    double *M = (double *) calloc(nrows_M*ncols_M, sizeof(double));
    assert(M != NULL);

    // compute local M = psi'*VxViX_g
    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        int Nde = E_k->nd_d;
        double *psi_k = DDBP_info->psi[spn_i][kpt][k];
        int *desc_psi_k = DDBP_info->desc_psi[k];
        double *Xi_k = DDBP_info->Xi[k];
        double alpha = 1.0, beta = 0.0;
        int ONE = 1;
        // M = xk^T * yk + M
        // generally, pdgemm_subcomm is not faster than pdgemm here.
        pdgemm_subcomm(
            "T", "N", Nstates, Nstates, Nde, 1.0, psi_k, desc_psi_k, 
            Xi_k, desc_psi_k, 1.0, M, DDBP_info->desc_M,
            blacscomm, best_max_nproc(Nstates, Nde, "pdgemm")
        );
        // pdgemm_(
        //     "T", "N", &Nstates, &Nstates, &Nde, 
        //     &alpha, psi_k, &ONE, &ONE, desc_psi_k, 
        //     Xi_k, &ONE, &ONE, desc_psi_k, 
        //     &beta, M, &ONE, &ONE, DDBP_info->desc_M);
    }
    t2 = MPI_Wtime();
    t_M += (t2 - t1);

    t1 = MPI_Wtime();
    MPI_Comm bandcomm = DDBP_info->bandcomm;
    int rank_bandcomm, size_bandcomm;
    MPI_Comm_rank(bandcomm, &rank_bandcomm);
    MPI_Comm_size(bandcomm, &size_bandcomm);

    if (size_bandcomm > 1)
        MPI_Allreduce(MPI_IN_PLACE, M, 
            nrows_M*ncols_M, MPI_DOUBLE, MPI_SUM, bandcomm);

    for (int i = 0; i < nrows_M*ncols_M; i++) M[i] *= -1.0;

#if defined(USE_MKL) || defined(USE_SCALAPACK)
    int ONE = 1, info;
    double alpha = 1.0;
    // pdpotrf_("U", &Nstates, M, &ONE, &ONE, DDBP_info->desc_M, &info);
    if (nrows_M*ncols_M > 0) {
        info = LAPACKE_dpotrf (LAPACK_COL_MAJOR, 'U', Nstates, M, Nstates);
    }

#ifdef DEBUG
    if (!rank) printf("info = %d for LAPACKE_dpotrf.\n", info);
#endif
    assert(info == 0);
    
    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        int *desc_psi_k = DDBP_info->desc_psi[k];
        int Nde = E_k->nd_d;
        double *Xi_k = DDBP_info->Xi[k];
        pdtrsm_("R", "U", "N", "N", &Nde, &Nstates, &alpha, 
                M, &ONE, &ONE, DDBP_info->desc_M, 
                Xi_k, &ONE, &ONE, desc_psi_k);
    }
#endif 
    free(M);
    t2 = MPI_Wtime();
    t_Choltrsm += (t2 - t1);

#ifdef DEBUG
    if (!rank) printf("\n== Creating ACE operator == \n");
    if (!rank) printf("== Time for gathering ddbp basis accross blacscomm : %.3f ms\n", t_gather*1e3);
    if (PoissonSkipFlag == 0) {
        if (!rank) printf("== Time for solving poisson's equations            : %.3f ms\n", t_pois*1e3);
    } else {
        if (!rank) printf("== Poisson's equations are not solved and solution retrived\n");
    }
    if (!rank) printf("== Time for poisson's equation times DDBP orbitals : %.3f ms\n", t_poisX*1e3);
    if (!rank) printf("== Time for calculating VxVi and reduce            : %.3f ms\n", t_VxVi*1e3);
    if (!rank) printf("== Time for calculating VxVi times DDBP orbitals   : %.3f ms\n", t_VxViX*1e3);
    if (!rank) printf("== Time for reducing VxVi and extract components   : %.3f ms\n", t_VxViXReduce*1e3);
    if (!rank) printf("== Time for creating M                             : %.3f ms\n", t_M*1e3);
    if (!rank) printf("== Time for choleskey and apply triag              : %.3f ms\n", t_Choltrsm*1e3);
#endif
}


/**
 * @brief   Compute operator for projected Vx with ACE
 */
void Calculate_Vx_DDBP_Hamiltonian_ACE(SPARC_OBJ *pSPARC, int count, int kpt, int spn_i)
{
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    MPI_Comm blacscomm = DDBP_info->blacscomm;

    int nelem = DDBP_info->n_elem_elemcomm;
    int Nxg = DDBP_info->Nx_global;
    int Nyg = DDBP_info->Ny_global;
    int Nzg = DDBP_info->Nz_global;
    int Ndg = DDBP_info->Nd_global;
    int Nstates = pSPARC->Nstates;

    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        int Nde = E_k->nd_d;
        int nALB = E_k->nALB;
        DDBP_HAMILT_ERBLKS *H_DDBP_Ek = &E_k->H_DDBP_Ek;
        DDBP_EXX *Exx_DDBP = &H_DDBP_Ek->Exx_DDBP;
        double *Xi_k = DDBP_info->Xi[k];
        int *desc_psi_k = DDBP_info->desc_psi[k];
        double *Xi_bar_k = Exx_DDBP->xi_bar;
        assert(Nstates == Exx_DDBP->ncol);
        double *v = E_k->v;
        
        // xi_bar = v^T \xi    
    #if defined(USE_MKL) || defined(USE_SCALAPACK)
        int desc_xibar[9];
        int ZERO = 0, info;
        int mb = max(1, nALB);
        int nb = max(1, Nstates);
        int llda = max(1, nALB);        

        descinit_(desc_xibar, &nALB, &Nstates, &mb, &nb, &ZERO, &ZERO, &DDBP_info->ictxt_blacs, &llda, &info);
        assert(info == 0);
        
        pdgemm_subcomm("T", "N", nALB, Nstates, Nde, 1.0, v, E_k->desc_v, Xi_k, desc_psi_k,
                    0.0, Xi_bar_k, desc_xibar, blacscomm, best_max_nproc(Nde, Nstates, "pdgemm"));
    #else
        assert(0);
    #endif    

        if (DDBP_info->npdm > 1) { // in the current implementation, this won't happen
            MPI_Allreduce(MPI_IN_PLACE, Xi_bar_k, nALB*nALB, MPI_DOUBLE, MPI_SUM, DDBP_info->dmcomm);
        }

        if (DDBP_info->npbasis > 1)
            MPI_Bcast(Xi_bar_k, nALB*Nstates, MPI_DOUBLE, 0, blacscomm); 
    }        
}


/**
 * @brief   Projected exact exchange operator (V'VxV) times a vector without ACE
 */
void DDBP_Hamiltonian_ExactExchange_vectors_mult_kernel_ACE(
    DDBP_HAMILTONIAN *H_DDBP, const DDBP_ARRAY *X, 
    DDBP_ARRAY *HX, MPI_Comm comm)
{
    int ncol = X->ncol;
    if(ncol == 0) return;
    
    int rank, size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    int nelem = X->nelem;
    int Nstates = H_DDBP->Nstates;
    double alpha = 1.0;
    double beta = 1.0;

    double *xibarTx = (double *) calloc(Nstates*ncol, sizeof(double));
    assert(xibarTx != NULL);

    for (int i = 0; i < nelem; i++) {
        assert(Nstates == H_DDBP->H_DDBP_Ek_list[i]->Exx_DDBP.ncol);
        int nALB = H_DDBP->H_DDBP_Ek_list[i]->Exx_DDBP.nALB;
        double *x_i = X->array[i];
        double *xi_bar = H_DDBP->H_DDBP_Ek_list[i]->Exx_DDBP.xi_bar;
        assert(nALB == X->nrows[i]);

        adaptive_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, Nstates, ncol, nALB,
                alpha, xi_bar, nALB, x_i, nALB, beta, xibarTx, Nstates);
    }

    if (size > 1)
        MPI_Allreduce(MPI_IN_PLACE, xibarTx, Nstates*ncol, MPI_DOUBLE, MPI_SUM, comm);

    alpha = -1.0*H_DDBP->exx_frac;
    for (int i = 0; i < nelem; i++) {
        int nALB = X->nrows[i];
        double *xi_bar = H_DDBP->H_DDBP_Ek_list[i]->Exx_DDBP.xi_bar;
        double *hx_i = HX->array[i];
        adaptive_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nALB, ncol, Nstates,
                alpha, xi_bar, nALB, xibarTx, Nstates, beta, hx_i, nALB);        
    }    
    free(xibarTx);    
}


/**
 * @brief   Evaluate exact exchange energy in case of DDBP-ACE
 */
void evaluate_exact_exchange_energy_DDBP_ACE(SPARC_OBJ *pSPARC)
{
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    // if (DDBP_info->elemcomm_index == -1) return;
    if (DDBP_info->elemcomm_index == -1 || DDBP_info->basiscomm_index == -1) return;
    MPI_Comm blacscomm = DDBP_info->blacscomm;
    int blacscomm_size, blacscomm_rank;
    MPI_Comm_size(blacscomm, &blacscomm_size);
    MPI_Comm_rank(blacscomm, &blacscomm_rank);

    int nelem = DDBP_info->n_elem_elemcomm;    
    int Nband = DDBP_info->n_band_bandcomm;
    int Nstates = pSPARC->Nstates;
    int spn_i, kpt;
    spn_i = kpt = 0;

    // Xi^T psi
    int nrow_M = DDBP_info->nrows_M;
    int ncol_M = DDBP_info->ncols_M;
    assert(nrow_M == Nstates);
    double *XiTpsi = (double *) calloc(nrow_M*ncol_M, sizeof(double));
    assert(XiTpsi != NULL);

    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        int Nde = E_k->nd_d;
        double *Xi_k = DDBP_info->Xi[k];
        double *psi_k = DDBP_info->psi[spn_i][kpt][k];
        int *desc_psi_k = DDBP_info->desc_psi[k];        
        double alpha = 1.0, beta = 1.0;
        int ONE = 1;

        pdgemm_subcomm(
            "T", "N", Nstates, Nstates, Nde, 1.0, Xi_k, desc_psi_k, 
            psi_k, desc_psi_k, 1.0, XiTpsi, DDBP_info->desc_M,
            blacscomm, best_max_nproc(Nstates, Nde, "pdgemm")
        );
    }

    MPI_Comm bandcomm = DDBP_info->bandcomm;
    int rank_bandcomm, size_bandcomm;
    MPI_Comm_rank(bandcomm, &rank_bandcomm);
    MPI_Comm_size(bandcomm, &size_bandcomm);

    if (size_bandcomm > 1)
        MPI_Allreduce(MPI_IN_PLACE, XiTpsi, nrow_M*ncol_M, MPI_DOUBLE, MPI_SUM, bandcomm);

    double Eexx = 0;
    for (int n = 0; n < ncol_M; n++) {
        double occ = DDBP_info->occ_outer[n+DDBP_info->band_start_index_M];
        for (int k = 0; k < nrow_M; k++) {
            Eexx -= occ * XiTpsi[k + n*nrow_M] * XiTpsi[k + n*nrow_M];
        }
    }

    // In current implementation, no row parallelization over Xi^Tpsi
    // otherwise, another row reduce is required
    if (blacscomm_size > 1)
        MPI_Allreduce(MPI_IN_PLACE, &Eexx, 1, MPI_DOUBLE, MPI_SUM, blacscomm);
    
    pSPARC->Eexx = Eexx * pSPARC->exx_frac;
    free(XiTpsi);
}


///////////////////////////////////////////////////////////////
//////////////////  Preparation and tools /////////////////////
///////////////////////////////////////////////////////////////

/**
 * @brief   update the usefock flag everywhere
 */
void update_usefock(SPARC_OBJ *pSPARC)
{
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    DDBP_info->H_DDBP.usefock = pSPARC->usefock;

    int nelem = DDBP_info->n_elem_elemcomm;
    if (pSPARC->ExxAcc == 1) {
        for (int k = 0; k < nelem; k++) {
            DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
            SPARC_OBJ *ESPRC_k = E_k->ESPRC;
            ESPRC_k->usefock = pSPARC->usefock;
        }
    } else {
        // use PBE basis here, thus don't update usefock
    }

    // forced skip_pois_flag = 0 for every outer loop
    DDBP_info->skip_pois_flag = 0;
    DDBP_info->vfsps_saved_flag = 0;
}


/**
 * @brief   Copy basis in other basiscomm for all elements
 *
 *          The default comm is DDBP_info->blacscomm
 */
void copy_ddbp_basis_xorb(SPARC_OBJ *pSPARC)
{
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    int nelem = DDBP_info->n_elem_elemcomm;
    int ncol = DDBP_info->n_basis_basiscomm;
    int Nstates = pSPARC->Nstates;
    int Nkpts = pSPARC->Nkpts_kptcomm;
    int Nspin = pSPARC->Nspin_spincomm;

    // copy ddbp basis
    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        int Nd_d = E_k->nd_d;
        memcpy(E_k->v_outer, E_k->v, Nd_d*ncol*sizeof(double));
    }
    
    // copy xorb
    for(int spn_i = 0; spn_i < Nspin; spn_i++) {
        for (int kpt = 0; kpt < Nkpts; kpt++) {
            DDBP_ARRAY *X_ks = &DDBP_info->xorb[spn_i][kpt];
            DDBP_ARRAY *Xouter_ks = &DDBP_info->xorb_outer[spn_i][kpt];
            void copy_DDBP_Array(double a, const DDBP_ARRAY *X, DDBP_ARRAY *Y);
            copy_DDBP_Array(1.0, X_ks, Xouter_ks);
        }
    }
}


/**
 * @brief   Extract ACE oeprator for all elements
 */
void ExtractElementACE(SPARC_OBJ *pSPARC)
{
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    // if (DDBP_info->elemcomm_index == -1) return;
    if (DDBP_info->elemcomm_index == -1 || DDBP_info->basiscomm_index == -1) return;

    MPI_Comm comm = DDBP_info->bandcomm;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int Nband = DDBP_info->n_band_bandcomm;
    int band_start_index = DDBP_info->band_start_index;
    int nelem = DDBP_info->n_elem_elemcomm;
    int Nxg = DDBP_info->Nx_global;
    int Nyg = DDBP_info->Ny_global;
    int Nzg = DDBP_info->Nz_global;
    int Ndg = DDBP_info->Nd_global;
    int Nstates = pSPARC->Nstates;

    double *Xi_full = (double *) calloc(Nband*Ndg, sizeof(double));
    assert(Xi_full != NULL);

    // collect full Xi
    // TODO: No need to collect all Xi. 
    // TODO: Implement a clever way to collect required parts of Xi
    gather_ESvec_and_E2B(DDBP_info, Nband, DDBP_info->Xi, Xi_full);

    int gridsizes[3] = {Nxg, Nyg, Nzg};
    int BCs[3] = {DDBP_info->BCx, DDBP_info->BCy, DDBP_info->BCz};

    // collect for Xi in dmcomm
    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        SPARC_OBJ *ESPRC_k = E_k->ESPRC;
        int Nde_ex = E_k->nd_ex_d;
        assert(Nde_ex == E_k->nd_ex);
        double *Xi_temp = (double *) calloc(Nde_ex*Nband, sizeof(double));
        for (int n = 0; n < Nband; n++) {
            double *Xi_full_n = Xi_full + n*Ndg;  
            extract_local_from_global(gridsizes, BCs, Xi_full_n, E_k->DMVert_ex,
                E_k->is_ex, E_k->js_ex, E_k->ks_ex, Xi_temp + n*Nde_ex);
        }
        // collect all bands of local Xi
        gather_horizontally(
            Nde_ex, Nband, band_start_index, Nstates, DDBP_info->npband, 
            Xi_temp, ESPRC_k->Xi, DDBP_info->blacscomm);
        free(Xi_temp);
    }
    free(Xi_full);

    // collect for Xi in kptcomm_topo
    // gather all bands of Xi
    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        SPARC_OBJ *ESPRC_k = E_k->ESPRC;
        int Nde_ex = E_k->nd_ex_d;

        int Nx_ex = E_k->nx_ex_d;
        assert(Nx_ex == E_k->nx_ex);
        int Ny_ex = E_k->ny_ex_d;
        assert(Ny_ex == E_k->ny_ex);
        int Nz_ex = E_k->nz_ex_d;
        assert(Nz_ex == E_k->nz_ex);

        int nd_ex_d_topo = E_k->nd_ex_d_topo;
        int nx_ex_d_topo = E_k->nx_ex_d_topo;
        int ny_ex_d_topo = E_k->ny_ex_d_topo;
        int nz_ex_d_topo = E_k->nz_ex_d_topo;

        //  only correct when no parallization in dmcomm
        for (int n = 0; n < Nstates; n++) {
            double *Xi_dmcomm_n = ESPRC_k->Xi + n*Nde_ex;
            double *Xi_kptcomm_topo_n = ESPRC_k->Xi_kptcomm_topo + n*nd_ex_d_topo;
            extract_subgrid(Xi_dmcomm_n, Nx_ex, Ny_ex, Nz_ex,
                E_k->DMVert_ex_topo[0],E_k->DMVert_ex_topo[1],E_k->DMVert_ex_topo[2],
                E_k->DMVert_ex_topo[3],E_k->DMVert_ex_topo[4],E_k->DMVert_ex_topo[5],
                Xi_kptcomm_topo_n, nx_ex_d_topo, ny_ex_d_topo, nz_ex_d_topo,
                0, nx_ex_d_topo-1, 0, ny_ex_d_topo-1, 0, nz_ex_d_topo-1);
        }        
    }
}


/**
 * @brief   Gather ddbp basis in other basiscomm
 *
 *          The default comm is DDBP_info->blacscomm
 */
void gather_ddbp_basis(DDBP_INFO *DDBP_info, 
    DDBP_ELEM *E_k, double *v, double *v_full, MPI_Comm comm)
{
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    int Nde = E_k->nd_d;
    int nbasis = DDBP_info->n_basis_basiscomm;
    int start_index = DDBP_info->basis_start_index;
    int shift = Nde * start_index;
    
    memcpy(v_full + shift, v, Nde*nbasis*sizeof(double));
    int *recvcounts, *displs;
    recvcounts = (int *) calloc(size, sizeof(int));
    displs = (int *) calloc(size, sizeof(int));
    
    // find n_basis_basiscomm in each basiscomm
    int Nb_k = DDBP_info->nALB_tot / DDBP_info->Ne_tot;
    int basis_start_index, basis_end_index, n_basis_basiscomm;

    displs[0] = 0;
    for (int i = 0; i < size; i++) {
        assign_task_BLCYC_fashion(
            Nb_k, DDBP_info->npbasis, i,
            &basis_start_index, &basis_end_index, &n_basis_basiscomm
        );
        if (i == rank) assert(nbasis == n_basis_basiscomm);
        recvcounts[i] = n_basis_basiscomm * Nde;
        if (i < size-1)
            displs[i+1] = displs[i] + recvcounts[i];
    }

    // TODO: change to NON-blocking
    MPI_Allgatherv(MPI_IN_PLACE, 1, MPI_DOUBLE, v_full, recvcounts, displs, MPI_DOUBLE, comm);
    free(recvcounts);
    free(displs);
}


/**
 * @brief   Gather element-wise vectors from all other elements 
 *          and rearrange into band-wise distribution
 *
 *          The default comm is DDBP_info->bandcomm
 */
void gather_ESvec_and_E2B(DDBP_INFO *DDBP_info, int ncol, double **vec_loc, double *vec_BS)
{
    MPI_Comm comm = DDBP_info->bandcomm;
    assert(comm != MPI_COMM_NULL);
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    assert(size == DDBP_info->npelem);

    int Nd_global = DDBP_info->Nd_global;
    int nelem = DDBP_info->n_elem_elemcomm;
    int spn_i, kpt;
    spn_i = kpt = 0;

    double *vec_full = (double *) calloc(ncol * Nd_global, sizeof(double));
    assert(vec_full != NULL);

    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        int global_index = k + DDBP_info->elem_start_index;
        int Nd_d = E_k->nd_d;
        memcpy(vec_full + DDBP_info->Ndedisp[global_index]*ncol, 
               vec_loc[k], Nd_d*ncol*sizeof(double));
    }

    int *recvcounts, *displs;
    recvcounts = (int *) calloc(size, sizeof(int));
    displs = (int *) calloc(size, sizeof(int));
    assert(recvcounts != NULL && displs != NULL);
    
    displs[0] = 0;
    for (int i = 0; i < size; i++) {
        int Nd_d_total = 0;
        for (int k = 0; k < DDBP_info->Nelem_list[i]; k++) {
            int global_index = k + DDBP_info->elemstart_list[i];
            Nd_d_total += DDBP_info->Ndelist[global_index];
        }
        recvcounts[i] = Nd_d_total * ncol;
        if (i < size - 1)
            displs[i+1] = displs[i] + recvcounts[i];
    }

    // TODO: change to non-blocking
    MPI_Allgatherv(MPI_IN_PLACE, 1, MPI_DOUBLE, vec_full, recvcounts, displs, MPI_DOUBLE, comm);

    free(recvcounts);
    free(displs);

    // rearrange vec_full into band-wise partition
    ElementWise2BandWise_local(DDBP_info, ncol, vec_full, vec_BS);
    free(vec_full);
}


/**
 * @brief   gather full psi for each bandcomm
 */
void gather_psi_and_E2B(SPARC_OBJ *pSPARC)
{
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    // if (DDBP_info->elemcomm_index == -1) return;
    if (DDBP_info->elemcomm_index == -1 || DDBP_info->basiscomm_index == -1) return;
    int Nband = DDBP_info->n_band_bandcomm;
    gather_ESvec_and_E2B(
            DDBP_info, Nband, 
            DDBP_info->psi[0][0], DDBP_info->psi_full);
}


/**
 * @brief   Copy occupations for DDBP
 */
void copy_occ_DDBP(SPARC_OBJ *pSPARC)
{
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    int Nstates = pSPARC->Nstates;
    // copy occpuations
    memcpy(DDBP_info->occ_outer, pSPARC->occ, Nstates*sizeof(double));
}

/**
 * @brief   rearrange array from element-wise to band-wise
 */
void ElementWise2BandWise_local(DDBP_INFO *DDBP_info, int ncol, double *x_in, double *x_out)
{
    int npelem = DDBP_info->npelem;
    int Nxg = DDBP_info->Nx_global;
    int Nyg = DDBP_info->Ny_global;
    int Nzg = DDBP_info->Nz_global;
    int Ndg = DDBP_info->Nd_global;

    double *x_in_ = x_in;
    for (int i = 0; i < npelem; i++) {
        for (int k = 0; k < DDBP_info->Nelem_list[i]; k++) {
            int global_index = DDBP_info->elemstart_list[i] + k;
            int *DMvertices = DDBP_info->DMvert_list[global_index];
            int DMnx = DMvertices[1] - DMvertices[0] + 1;
            int DMny = DMvertices[3] - DMvertices[2] + 1;
            int DMnz = DMvertices[5] - DMvertices[4] + 1;
            int DMnd = DMnx * DMny * DMnz;
            assert(DDBP_info->Ndelist[global_index] == DMnd);

            for (int n = 0; n < ncol; n++) {
                double *x_out_ = x_out + n*Ndg;
                extract_subgrid(
                    x_in_, DMnx, DMny, DMnz,
                    0, DMnx-1, 0, DMny-1, 0, DMnz-1,
                    x_out_, Nxg, Nyg, Nzg,
                    DMvertices[0], DMvertices[1],DMvertices[2], 
                    DMvertices[3],DMvertices[4], DMvertices[5]
                );
                x_in_ += DMnd;
            }
        }
    }
}

/**
 * @brief   rearrange local array from band-wise to element-wise
 */
void BandWise2ElementWise_local(DDBP_INFO *DDBP_info, int ncol, double *x_in, double *x_out)
{
    int npelem = DDBP_info->npelem;
    int Nxg = DDBP_info->Nx_global;
    int Nyg = DDBP_info->Ny_global;
    int Nzg = DDBP_info->Nz_global;
    int Ndg = DDBP_info->Nd_global;

    double *x_out_ = x_out;
    for (int i = 0; i < npelem; i++) {
        for (int k = 0; k < DDBP_info->Nelem_list[i]; k++) {
            int global_index = DDBP_info->elemstart_list[i] + k;
            int *DMvertices = DDBP_info->DMvert_list[global_index];
            int DMnx = DMvertices[1] - DMvertices[0] + 1;
            int DMny = DMvertices[3] - DMvertices[2] + 1;
            int DMnz = DMvertices[5] - DMvertices[4] + 1;
            int DMnd = DMnx * DMny * DMnz;
            for (int n = 0; n < ncol; n++) {
                double *x_in_ = x_in + n * Ndg;
                extract_subgrid(
                    x_in_, Nxg, Nyg, Nzg,
                    DMvertices[0], DMvertices[1],DMvertices[2], 
                    DMvertices[3],DMvertices[4], DMvertices[5],
                    x_out_, DMnx, DMny, DMnz,
                    0, DMnx-1, 0, DMny-1, 0, DMnz-1
                );
                x_out_ += DMnd;
            }
        }
    }
}

/**
 * @brief   Find local Nd for each elements and DMvertices
 */
int find_local_Nde(DDBP_INFO *DDBP_info, int *coords, double *meshes, int *DMvertices)
{
    int ii_E = coords[0];
    int jj_E = coords[1];
    int kk_E = coords[2];

    double xs = ii_E * DDBP_info->Lex;
    double xe = (ii_E + 1) * DDBP_info->Lex;
    double ys = jj_E * DDBP_info->Ley;
    double ye = (jj_E + 1) * DDBP_info->Ley;
    double zs = kk_E * DDBP_info->Lez;
    double ze = (kk_E + 1) * DDBP_info->Lez;

    double dx = meshes[0];
    double dy = meshes[1];
    double dz = meshes[2];

    DMvertices[0] = ceil(xs / dx - TEMP_TOL);
    DMvertices[1] = ceil(xe / dx - TEMP_TOL) - 1;
    DMvertices[2] = ceil(ys / dy - TEMP_TOL);
    DMvertices[3] = ceil(ye / dy - TEMP_TOL) - 1;
    DMvertices[4] = ceil(zs / dz - TEMP_TOL);
    DMvertices[5] = ceil(ze / dz - TEMP_TOL) - 1;

    int nx = DMvertices[1] - DMvertices[0] + 1;
    int ny = DMvertices[3] - DMvertices[2] + 1;
    int nz = DMvertices[5] - DMvertices[4] + 1;
    int nd = nx * ny * nz;
    return nd;
}


/**
 * @brief   Transfer basis-wise array to element-wise array
 */
void BasisWise2ElementWise(DDBP_INFO *DDBP_info, double *x_in, double *x_out)
{
    MPI_Comm comm = DDBP_info->bandcomm;
    assert(comm != MPI_COMM_NULL);
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    assert(size == DDBP_info->npelem);

    int Ndg = DDBP_info->Nd_global;
    int nbasis = DDBP_info->n_basis_basiscomm;
    int nelem_loc = DDBP_info->n_elem_elemcomm;

    // find alltoallv parameters
    int *scounts, *sdispls, *rcounts, *rdispls;
    scounts = (int *) calloc(size, sizeof(int));
    sdispls = (int *) calloc(size, sizeof(int));
    rcounts = (int *) calloc(size, sizeof(int));
    rdispls = (int *) calloc(size, sizeof(int));
    assert(scounts != NULL && sdispls != NULL && 
           rcounts != NULL && rdispls != NULL);
    
    sdispls[0] = 0;
    for (int i = 0; i < size; i++) {
        int nelem_dest = DDBP_info->Nelem_list[i];
        int elem_start_index = DDBP_info->elemstart_list[i];
        int Nde_tot = 0;
        for (int k = 0; k < nelem_dest; k++) {
            int global_index = elem_start_index + k;
            Nde_tot += DDBP_info->Ndelist[global_index];
        }
        scounts[i] = Nde_tot * nelem_loc * nbasis;
        if (i < size-1) 
            sdispls[i+1] = sdispls[i] + scounts[i];
    }

    int Nde_tot_loc = 0;
    for (int k = 0; k < nelem_loc; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        Nde_tot_loc += E_k->nd_d;
    }
    rdispls[0] = 0;
    for (int i = 0; i < size; i++) {
        int nelem_src = DDBP_info->Nelem_list[i];
        rcounts[i] = Nde_tot_loc * nelem_src * nbasis;
        if (i < size-1)
            rdispls[i+1] = rdispls[i] + rcounts[i];
    }

#ifdef DEBUG
    int rcounts_tot = 0;
    for (int i = 0; i < size; i++) {
        rcounts_tot += rcounts[i];
    }
    assert(rcounts_tot == Nde_tot_loc*nbasis*DDBP_info->Ne_tot);
#endif

    // TODO: change to non-blocking
    MPI_Alltoallv(x_in, scounts, sdispls, MPI_DOUBLE, 
                  x_out, rcounts, rdispls, MPI_DOUBLE, comm);
    
    free(scounts);
    free(sdispls);
    free(rcounts);
    free(rdispls);
}

/**
 * @brief   Rearrange Column-Wise array with Multiple Elements
 */
void RearrangeColumnWiseMultipleElements(DDBP_INFO *DDBP_info, double *x_in, double *x_out)
{
    int npelem = DDBP_info->npelem;
    int nbasis = DDBP_info->n_basis_basiscomm;
    int nelem = DDBP_info->n_elem_elemcomm;
    int Nxg = DDBP_info->Nx_global;
    int Nyg = DDBP_info->Ny_global;
    int Nzg = DDBP_info->Nz_global;
    int Ndg = DDBP_info->Nd_global;

    double *x_out_ = x_out;
    for (int i = 0; i < npelem; i++) {
        int start_inex = DDBP_info->elemstart_list[i];
        for (int j = 0; j < DDBP_info->Nelem_list[i]; j++) {
            int global_index = j + start_inex;
            int *DMVertices = DDBP_info->DMvert_list[global_index];
            int DMnx = DMVertices[1] - DMVertices[0] + 1;
            int DMny = DMVertices[3] - DMVertices[2] + 1;
            int DMnz = DMVertices[5] - DMVertices[4] + 1;
            int DMnd = DDBP_info->Ndelist[global_index];
            assert(DMnd == DDBP_info->Ndelist[global_index]);

            for (int k = 0; k < nelem; k++) {
                // extract blocks 
                for (int col = 0; col < nbasis; col++) {
                    double *x_in_ = x_in + col*Ndg + k*Ndg*nbasis;
                    extract_subgrid(
                        x_in_, Nxg, Nyg, Nzg,
                        DMVertices[0], DMVertices[1], DMVertices[2],
                        DMVertices[3], DMVertices[4], DMVertices[5],
                        x_out_, DMnx, DMny, DMnz,
                        0, DMnx-1, 0, DMny-1, 0, DMnz-1
                    );
                    x_out_ += DMnd;
                }
            }
        }
    }
}

/**
 * @brief   Poisson's solver for all pairs of local basis for DDBP
 */
void poissonSolveDDBP(SPARC_OBJ *pSPARC, 
    DDBP_ELEM *E_k, int ncol, double *v, double *v2, double *solg)
{
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    assert(ncol == DDBP_info->n_basis_basiscomm);
    if (ncol == 0) return;

    int Nde = E_k->nd_d;
    int Ndg = pSPARC->Nd;
    int Nxg = pSPARC->Nx;
    int Nyg = pSPARC->Ny;
    int Nzg = pSPARC->Nz;

    double *rhs = (double *) calloc(Nde, sizeof(double));
    double *rhsg = (double *) calloc(ncol*Ndg, sizeof(double));    
    memset(rhsg, 0, sizeof(double)* ncol*Ndg);
    assert(rhs != NULL && rhsg != NULL);

    for (int n = 0; n < ncol; n++) {
        double *v_n = v + n * Nde;
        for (int i = 0; i < Nde; i++) {
            rhs[i] = v_n[i] * v2[i];
        }        
        
        int DMnx = E_k->nx_d;
        int DMny = E_k->ny_d;
        int DMnz = E_k->nz_d;
        double *rhsg_ = rhsg + n *Ndg;
        // construct the rhs global vector 
        extract_subgrid(
            rhs, DMnx, DMny, DMnz,
            0, DMnx-1, 0, DMny-1, 0, DMnz-1, 
            rhsg_, Nxg, Nyg, Nzg,
            E_k->is, E_k->ie, E_k->js, E_k->je, E_k->ks, E_k->ke
        );
    }

    // solve poisson's equation by FFT    
    pois_fft(pSPARC, rhsg, pSPARC->pois_const, ncol, solg);
    free(rhs);
    free(rhsg);
}


/**
 * @brief   Global vector (Ndg x nALB) times inner DDBP orbitals
 */
void GlobalvecTimesX(
    DDBP_INFO *DDBP_info, DDBP_ELEM *E_k, int k, 
    const double alpha, double *v, DDBP_ARRAY *xorb,
    const double beta, double *vX)
{
    int spn_i, kpt;
    spn_i = kpt = 0;

    double *xorb_k = xorb->array[k];
    int *descX = DDBP_info->desc_xorb;
    int *desc_solX = DDBP_info->desc_solX[k];

    // pdgemm_subcomm will be much slower than pdgemm for tall matrix x flat matrix
    // when processors in the communicator less than a threshold. Generally, the larger
    // the matrix, the smaller the threshold. We hardcoded 8 here. it works for systems 
    // upto 400 atoms systems. 

    int nproc_rowcomm;
    MPI_Comm_size(DDBP_info->rowcomm, &nproc_rowcomm);
    int proc_active = (DDBP_info->is_active_rowcomm) && (DDBP_info->elemcomm_index >= 0);

    if (nproc_rowcomm >= 8) {
        DDBP_Element_Basis_Coeff_Multiply(
            alpha, v, E_k->desc_sol, xorb_k, descX,
            beta, vX, desc_solX, DDBP_info->rowcomm, proc_active
        );
    } else {
        // simply use pdgemm 
        int ONE = 1;
        int M = E_k->desc_sol[2], N = descX[3], K = E_k->desc_sol[3];
        assert (M == desc_solX[2] && N == desc_solX[3] && K == descX[2]);

#if defined(USE_MKL) || defined(USE_SCALAPACK)        
        // TODO: Careful about this descriptor. It's only correct in current implementation
        int desc[9];
        int mb = max(1, M);
        int nb = (K - 1) / DDBP_info->npbasis + 1;
        int ZERO = 0;
        // set up descriptor for storage of extended basis funcs in ictxt_blacscomm (original)
        int llda = max(1, M);
        int info;
        descinit_(desc, &M, &K,
                &mb, &nb, &ZERO, &ZERO, &DDBP_info->ictxt_elemcomm, &llda, &info);
        assert(info == 0);

        pdgemm_("N", "N", &M, &N, &K, &alpha, 
                v, &ONE, &ONE, desc, xorb_k, 
                &ONE, &ONE, descX, &beta, vX, &ONE, &ONE, 
                desc_solX);
#endif
    }
}


/**
 * @brief   Compute Vx times local basis
 */
void computeVxViDDBP(DDBP_INFO *DDBP_info, 
    double *solX, double *VxVi, int element_index, int basis_index)
{
    MPI_Comm comm = DDBP_info->blacscomm;
    int rank, size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    int Nband = DDBP_info->n_band_bandcomm;
    int Ndg = DDBP_info->Nd_global;
    int nALB = DDBP_info->elem_list[element_index].nALB;
    double dV = DDBP_info->dV;

    double *VxVil = (double *) calloc(Ndg, sizeof(double));
    memset(VxVil, 0, sizeof(double)* Ndg);
    assert(VxVil != NULL);

    for (int n = 0; n < Nband; n++) {
        int global_band_index = n + DDBP_info->band_start_index;
        double occ = DDBP_info->occ_outer[global_band_index];
        double *psi_n = DDBP_info->psi_full + n * Ndg;
        double *solX_ = solX + n * Ndg;
        for (int i = 0; i < Ndg; i++) {
            VxVil[i] -= occ*psi_n[i]*solX_[i]/dV;
        }
    }
    
    // find destination for reduce of VxVi
    int dest = basis_index_to_basiscomm_rank(nALB, size, basis_index);
    int shift = basis_index - DDBP_info->basis_start_index;
    if (dest == rank) assert(shift >= 0);
    double *VxVi_ = VxVi + shift * Ndg;

    MPI_Reduce(VxVil, VxVi_, Ndg, MPI_DOUBLE, MPI_SUM, dest, comm);
    free(VxVil);
}



/**
 * @brief   Compute Vx times local basis non-blocking way
 */
void computeVxViDDBP_nonblocking(DDBP_INFO *DDBP_info, 
    double *solX, double *VxVi, int element_index, int basis_index, MPI_Request *request)
{
    MPI_Comm comm = DDBP_info->blacscomm;
    int rank, size;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    int Nband = DDBP_info->n_band_bandcomm;
    int Ndg = DDBP_info->Nd_global;
    int nALB = DDBP_info->elem_list[element_index].nALB;
    double dV = DDBP_info->dV;

    double *VxVil = (double *) calloc(Ndg, sizeof(double));
    memset(VxVil, 0, sizeof(double)* Ndg);
    assert(VxVil != NULL);

    for (int n = 0; n < Nband; n++) {
        int global_band_index = n + DDBP_info->band_start_index;
        double occ = DDBP_info->occ_outer[global_band_index];
        double *psi_n = DDBP_info->psi_full + n * Ndg;
        double *solX_ = solX + n * Ndg;
        for (int i = 0; i < Ndg; i++) {
            VxVil[i] -= occ*psi_n[i]*solX_[i]/dV;
        }
    }
    
    // find destination for reduce of VxVi
    int dest = basis_index_to_basiscomm_rank(nALB, size, basis_index);
    int shift = basis_index - DDBP_info->basis_start_index;
    if (dest == rank) assert(shift >= 0);
    double *VxVi_ = VxVi + shift * Ndg;

    MPI_Ireduce(VxVil, VxVi_, Ndg, MPI_DOUBLE, MPI_SUM, dest, comm, request);
    free(VxVil);
}

/**
 * @brief   Get rank in basiscomm given basis index
 */
int basis_index_to_basiscomm_rank(int nALB, int nproc, int basis_index)
{   
    int basis_start_index, basis_end_index, n_basis_basiscomm;
    int rank = -1;

    for (int i = 0; i < nproc; i++) {
        assign_task_BLCYC_fashion(
            nALB, nproc, i,
            &basis_start_index, &basis_end_index, &n_basis_basiscomm
        );
        if (basis_end_index >= basis_index) {
            rank = i;
            break;
        }
    }    
    return rank;
}

/**
 * @brief   Set up havx exchange information 
 */
void setup_havxX_DDBP_Array(DDBP_INFO *DDBP_info, DDBP_ARRAY *X, MPI_Comm comm)
{
    havxX_t *havxX_info = &X->havxX_info;
    if (DDBP_info->H_DDBP.usefock == 0 || DDBP_info->elemcomm_index == -1) {
        havxX_info->havxX_flag = 0;
        return;
    } else {
        havxX_info->havxX_flag = 1;
    }
    assert(comm != MPI_COMM_NULL);
    int size;
    MPI_Comm_size(comm, &size);
    assert(size == DDBP_info->npelem);

    int nALB = DDBP_info->nALB_tot/DDBP_info->Ne_tot;
    int ncol = X->ncol;
    havxX_info->nrows_total = DDBP_info->nALB_tot;
    havxX_info->blksz = nALB * ncol;
    havxX_info->nelem = DDBP_info->n_elem_elemcomm;
    havxX_info->size = size;
    int nelem = havxX_info->nelem;

    havxX_info->global_indexes = (int *) calloc(nelem, sizeof(int));
    assert(havxX_info->global_indexes != NULL);
    for (int i = 0; i < nelem; i++)
        havxX_info->global_indexes[i] = i + DDBP_info->elem_start_index;
    
    havxX_info->recvcounts = (int *) calloc(size, sizeof(int));
    havxX_info->displs = (int *) calloc(size, sizeof(int));
    assert(havxX_info->recvcounts != NULL && havxX_info->displs != NULL);
    
    havxX_info->displs[0] = 0;
    for (int i = 0; i < size; i++) {
        havxX_info->recvcounts[i] = havxX_info->blksz * DDBP_info->Nelem_list[i];
        if (i < size - 1)
            havxX_info->displs[i+1] = havxX_info->displs[i] + havxX_info->recvcounts[i];
    }
}

/**
 * @brief   copy havxX_info
 */
void copy_havxX_info(const havxX_t *havxX_src, havxX_t *havxX_dest)
{
    havxX_dest->havxX_flag = havxX_src->havxX_flag;
    if (havxX_src->havxX_flag == 0) return;
    
    havxX_dest->nrows_total = havxX_src->nrows_total;
    havxX_dest->blksz = havxX_src->blksz;
    int nelem = havxX_dest->nelem = havxX_src->nelem;
    int size = havxX_dest->size = havxX_src->size;

    havxX_dest->global_indexes = (int *) calloc(nelem, sizeof(int));
    havxX_dest->recvcounts = (int *) calloc(size, sizeof(int));
    havxX_dest->displs = (int *) calloc(size, sizeof(int));
    assert(havxX_dest->global_indexes != NULL 
        && havxX_dest->recvcounts != NULL
        && havxX_dest->displs != NULL);

    for (int i = 0; i < nelem; i++) {
        havxX_dest->global_indexes[i] = havxX_src->global_indexes[i];
    }

    for (int i = 0; i < size; i++) {
        havxX_dest->recvcounts[i] = havxX_src->recvcounts[i];
        havxX_dest->displs[i] = havxX_src->displs[i];
    }
}

/**
 * @brief   free havxX_info
 */
void free_havxX_info(DDBP_ARRAY *X)
{
    havxX_t *havxX_info = &X->havxX_info;
    if (havxX_info->havxX_flag == 0) return;

    free(havxX_info->global_indexes);
    free(havxX_info->recvcounts);
    free(havxX_info->displs);
}

/**
 * @brief   Gather components from all elemnts for a DDBP array
 */
void DDBP_element_allelements_allgatherv(const DDBP_ARRAY *X, double *buff, MPI_Comm comm)
{
    const havxX_t *havxX_info = &X->havxX_info;
    assert(havxX_info->havxX_flag == 1);

    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    assert(size == havxX_info->size);

    int nelem = havxX_info->nelem;
    assert(nelem == X->nelem);
    for (int k = 0; k < nelem; k++) {
        int nrow = X->nrows[k];
        int ncol = X->ncol;
        int blksz = havxX_info->blksz;
        assert(nrow * ncol == blksz);
        int shift = havxX_info->global_indexes[k]*blksz;

        memcpy(buff+shift, X->array[k], blksz*sizeof(double));
    }

    // TODO: change to non-blocking
    int *recvcounts = havxX_info->recvcounts;
    int *displs = havxX_info->displs;
    if (size > 1)
        MPI_Allgatherv(MPI_IN_PLACE, 1, MPI_DOUBLE, 
                buff, recvcounts, displs, MPI_DOUBLE, comm);
}



/**
 * @brief   Reduce and extract local VxViX
 */
void ReduceExtractVxViX(DDBP_INFO *DDBP_info, double *VxViX_g, double **VxViX, int opt)
{
    // if opt == 1, use VxViX as output
    // if opt == 2, use E_k->Xi as output
    assert(opt == 1 || opt == 2);

    MPI_Comm comm = DDBP_info->bandcomm;
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    int Nband = DDBP_info->n_band_bandcomm;
    int Nxg = DDBP_info->Nx_global;
    int Nyg = DDBP_info->Ny_global;
    int Nzg = DDBP_info->Nz_global;
    int Ndg = DDBP_info->Nd_global;
    int nelem = DDBP_info->n_elem_elemcomm;

    if (Nband == 0) return;
    if (size > 1)
        MPI_Allreduce(MPI_IN_PLACE, VxViX_g, Ndg*Nband, MPI_DOUBLE, MPI_SUM, comm);

    for (int k = 0; k < nelem; k++) {        
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        int DMnx = E_k->nx_d;
        int DMny = E_k->ny_d;
        int DMnz = E_k->nz_d;
        int DMnd = E_k->nd_d;
        double *VxViX_k = (opt == 1) ? VxViX[k] : DDBP_info->Xi[k];

        for (int n = 0; n < Nband; n++) {
            double *VxViX_g_n = VxViX_g + n*Ndg;
            double *VxViX_k_n = VxViX_k + n*DMnd;
            // extract the local VxViX            
            extract_subgrid(
                VxViX_g_n, Nxg, Nyg, Nzg,
                E_k->is, E_k->ie, E_k->js, 
                E_k->je, E_k->ks, E_k->ke,
                VxViX_k_n, DMnx, DMny, DMnz,
                0, DMnx-1, 0, DMny-1, 0, DMnz-1
            );
        }
    }
}



/**
 * @brief   Gather vecotrs horizontally
 *
 * Note:    The default comm is blacscomm. 
 *          Generalization of gather_ddbp_basis function
 */
void gather_horizontally(
    int N1, int ncol_loc, int start_index, int ncol_tot, int np, 
    double *v, double *v_full, MPI_Comm comm)
{
    if (comm == MPI_COMM_NULL) return;
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    int shift = N1 * start_index;    
    memcpy(v_full + shift, v, N1*ncol_loc*sizeof(double));
    int *recvcounts, *displs;
    recvcounts = (int *) calloc(size, sizeof(int));
    displs = (int *) calloc(size, sizeof(int));
    
    // find ncol_loc in each comm    
    int basis_start_index, basis_end_index, ncol;

    displs[0] = 0;
    for (int i = 0; i < size; i++) {
        assign_task_BLCYC_fashion(
            ncol_tot, np, i,
            &basis_start_index, &basis_end_index, &ncol
        );
        if (i == rank) assert(ncol_loc == ncol);
        recvcounts[i] = ncol * N1;
        if (i < size-1)
            displs[i+1] = displs[i] + recvcounts[i];
    }

    // TODO: change to NON-blocking
    MPI_Allgatherv(MPI_IN_PLACE, 1, MPI_DOUBLE, v_full, recvcounts, displs, MPI_DOUBLE, comm);

    free(recvcounts);
    free(displs);
}


/**
 * @brief   Change Veff to PBE in hybrid DDBP basis calculation
 * 
 * Note:    Always use Recover_hybrid_Veff together
 */
void Change_PBE_Veff(SPARC_OBJ *pSPARC)
{
    pSPARC->usefock --;
    // save the PBE0 Veff
    size_t Veff_size = pSPARC->Nd_d_dmcomm * pSPARC->Nspin;
    if (pSPARC->MixingVariable == 1) { // potential mixing
        pSPARC->Veff_loc_dmcomm_temp = (double *) calloc(Veff_size, sizeof(double));
        assert(pSPARC->Veff_loc_dmcomm_temp != NULL);
        memcpy(pSPARC->Veff_loc_dmcomm_temp, pSPARC->Veff_loc_dmcomm, Veff_size*sizeof(double));
    }
    Calculate_Vxc(pSPARC);
    Calculate_Veff_loc_dmcomm_phi(pSPARC);
    int DMnd = pSPARC->Nd_d;
    for (int i = 0; i < pSPARC->Nspin; i++)
        Transfer_Veff_loc(pSPARC, pSPARC->Veff_loc_dmcomm_phi + i*DMnd, pSPARC->Veff_loc_dmcomm + i*pSPARC->Nd_d_dmcomm);
    pSPARC->usefock ++;
}


/**
 * @brief   Recover hybrid Veff in hybrid DDBP basis calculation
 * 
 * Note:    Always use Change_PBE_Veff together
 */
void Recover_hybrid_Veff(SPARC_OBJ *pSPARC)
{   
    // Note: In case of potential mixing, the Veff is the reuslt after mixing, thus we 
    // need to recover it exactly after finding DDBP basis. And e_xc etc. will be updated
    // again later in scf_loop, thus no worries about it. In case of density mixing, we 
    // need to update the e_xc etc. which have been changed when calculating PBE Vxc. 

    if (pSPARC->MixingVariable == 1) { // potential mixing
        size_t Veff_size = pSPARC->Nd_d_dmcomm * pSPARC->Nspin;
        memcpy(pSPARC->Veff_loc_dmcomm, pSPARC->Veff_loc_dmcomm_temp, Veff_size*sizeof(double));
        setup_ddbp_element_Veff(pSPARC);
        free(pSPARC->Veff_loc_dmcomm_temp);
    } else {
        Calculate_Vxc(pSPARC);
        Calculate_Veff_loc_dmcomm_phi(pSPARC);
        int DMnd = pSPARC->Nd_d;
        for (int i = 0; i < pSPARC->Nspin; i++)
            Transfer_Veff_loc(pSPARC, pSPARC->Veff_loc_dmcomm_phi + i*DMnd, pSPARC->Veff_loc_dmcomm + i*pSPARC->Nd_d_dmcomm);
        setup_ddbp_element_Veff(pSPARC);
    }
}

/**
 * @brief   Update the Vx part in H_DDBP after basis changed
 */
void Update_Vx_part_of_DDBP_Hamiltonian(SPARC_OBJ *pSPARC, int count, int kpt, int spn_i)
{
#ifdef DEBUG
    int rank;
    double t1, t2;
    t1 = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank) printf("\n== Update the Vx part in H_DDBP after basis changed == \n");
#endif

    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    if (DDBP_info->elemcomm_index == -1 || DDBP_info->basiscomm_index == -1 ||
        DDBP_info->dmcomm == MPI_COMM_NULL) {
        return;
    }
    
    if (pSPARC->ExxAcc == 0)
        Calculate_H_DDBP_Vx_block(pSPARC);
    else
        Calculate_Vx_DDBP_Hamiltonian_ACE(pSPARC, count, kpt, spn_i);

#ifdef DEBUG
    t2 = MPI_Wtime();
    if (!rank) printf("== Total time for Update the Vx: %.3f ms\n", (t2-t1)*1e3);
#endif
}


///////////////////////////////////////////////////////////////
//////////  Functions for parallelly solving poisson //////////
///////////////////////////////////////////////////////////////

/**
 * @brief   Solving poissons equations parallelly
 */
void solve_poisson_parallelly(SPARC_OBJ *pSPARC, DDBP_ELEM *E_k, double *sol)
{
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    DDBP_EXX_POIS_PARAL_G *pois_paral_g = &DDBP_info->pois_paral_g;
    DDBP_EXX_POIS_PARAL_D *pois_paral_d = &DDBP_info->pois_paral_d;
    int nbasis = DDBP_info->n_basis_basiscomm;
    int nelem = DDBP_info->n_elem_elemcomm;
    int Nd_global = DDBP_info->Nd_global;
    if (nbasis == 0) return;

    double t1, t2, t_cal, t_com;
    t_cal = t_com = 0;
    t1 = MPI_Wtime();
    // first solve locally
    solve_halfpair_poissons_equations(pSPARC, E_k, E_k->v, sol);
    t2 = MPI_Wtime();
    t_cal += (t2 - t1);

    // gather half other basis and solve all and them distribution solutions         
    int Nde = E_k->nd_d;
    int rcount_tot = pois_paral_g->rcount_tot * Nde;
    double *recvbuf = (double *) calloc(rcount_tot, sizeof(double));
    t1 = MPI_Wtime();
    gather_half_basis(pois_paral_g, Nde, nbasis, E_k->v, recvbuf);
    t2 = MPI_Wtime();
    t_com += (t2 - t1);

    t1 = MPI_Wtime();
    // solve all pair of poissons' equations by recvbuff and local v
    solve_allpair_poissons_equations(pSPARC, 
        E_k, pois_paral_g, E_k->v, recvbuf, sol);
    free(recvbuf);
    t2 = MPI_Wtime();
    t_cal += (t2 - t1);
    
    t1 = MPI_Wtime();
    // communication
    distribute_solutions(DDBP_info, pois_paral_d, sol, nbasis);
    t2 = MPI_Wtime();
    t_com += (t2 - t1);

#ifdef DEBUG
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank) printf("Solving half of poissons' equations took %.3f ms, (calc %.3f + com %.3f)\n",
        (t_cal + t_com)*1e3, t_cal*1e3, t_com*1e3);
#endif
}

/**
 * @brief   Distribute poissons solutions 
 */
void distribute_solutions(DDBP_INFO *DDBP_info, 
    DDBP_EXX_POIS_PARAL_D *pois_paral_d, double *sol, int n_basis)
{
    size_t scount_tot = pois_paral_d->scount_tot;
    size_t rcount_tot = pois_paral_d->rcount_tot;
    double *send_buff = (double*) calloc(scount_tot, sizeof(double));
    double *recv_buff = (double*) calloc(rcount_tot, sizeof(double));
    assert(send_buff != NULL && recv_buff != NULL);
    
    int Ndg = DDBP_info->Nd_global;
    int outdegree = pois_paral_d->outdegree;
    double *send_buff_ = send_buff;
    for (int i = 0; i < outdegree; i++) {
        int dest = pois_paral_d->destinations[i];
        int n_basis_dest = DDBP_info->n_basis_list[dest];
        int basis_start_index = DDBP_info->basis_start_index_list[dest];
        double *sol_out = sol + n_basis * basis_start_index * Ndg;
        size_t blksz = sizeof(double) * n_basis_dest * n_basis * Ndg;
        memcpy(send_buff_, sol_out, blksz);
        send_buff_ += n_basis_dest * n_basis * Ndg;
    }
    
    MPI_Comm dist_graph_comm = pois_paral_d->dist_graph_comm;
    if (dist_graph_comm != MPI_COMM_NULL) {
        MPI_Neighbor_alltoallv( send_buff, pois_paral_d->scount, pois_paral_d->sdispls, MPI_DOUBLE, 
                                recv_buff, pois_paral_d->rcount, pois_paral_d->rdispls, MPI_DOUBLE, 
                                dist_graph_comm );
    }

    // assemble received solutions
    int indegree = pois_paral_d->indegree;
    double *recv_buff_ = recv_buff;
    for (int i = 0; i < indegree; i++) {
        int src = pois_paral_d->sources[i];
        int n_basis_src = DDBP_info->n_basis_list[src];
        int basis_start_index = DDBP_info->basis_start_index_list[src];

        for (int n = 0; n < n_basis; n++) {
            double *sol_nX = sol + n * Ndg;
            for (int m = 0; m < n_basis_src; m++) {
                double *sol_nm = sol_nX + (basis_start_index + m) * n_basis * Ndg;
                memcpy(sol_nm, recv_buff_, Ndg*sizeof(double));
                recv_buff_ += Ndg;
            }
        }
    }

    free(send_buff);
    free(recv_buff);
}

/**
 * @brief   Solve poisson's equations construted by all pairs of basis
 */
void solve_allpair_poissons_equations(SPARC_OBJ *pSPARC, 
    DDBP_ELEM *E_k, DDBP_EXX_POIS_PARAL_G *pois_paral_g, 
    double *v, double *recvbuf, double *sol)
{
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    int Nde = E_k->nd_d;
    int Ndg = pSPARC->Nd;
    int Nxg = pSPARC->Nx;
    int Nyg = pSPARC->Ny;
    int Nzg = pSPARC->Nz;
    int basis_start_index = DDBP_info->basis_start_index;
    int n_basis = DDBP_info->n_basis_basiscomm;

    double *rhs = (double *) calloc(Nde, sizeof(double));
    double *rhsg = (double *) calloc(Ndg*n_basis, sizeof(double));
    assert(rhs != NULL && rhsg != NULL);

    int indegree = pois_paral_g->indegree;
    double *recvbuf_ = recvbuf;
    for (int i = 0; i < indegree; i++) {
        int src = pois_paral_g->sources[i];
        int rcount = pois_paral_g->rcount[i];
        int basis_start_index = DDBP_info->basis_start_index_list[src];
        for (int n = 0; n < rcount; n++) {
            int n_global = n + basis_start_index;
            
            for (int m = 0; m < n_basis; m++) {
                double *v_m = v + m * Nde;
                for (int j = 0; j < Nde; j++) {
                    rhs[j] = v_m[j] * recvbuf_[j];
                }
                int DMnx = E_k->nx_d;
                int DMny = E_k->ny_d;
                int DMnz = E_k->nz_d;
                double *rhsg_m = rhsg + m *Ndg;
                // construct the rhs global vector 
                extract_subgrid(
                    rhs, DMnx, DMny, DMnz,
                    0, DMnx-1, 0, DMny-1, 0, DMnz-1, 
                    rhsg_m, Nxg, Nyg, Nzg,
                    E_k->is, E_k->ie, E_k->js, E_k->je, E_k->ks, E_k->ke
                );
            }
            recvbuf_ += Nde;

            double *sol_Xn = sol + (basis_start_index + n) * n_basis * Ndg;            
            pois_fft(pSPARC, rhsg, pSPARC->pois_const, n_basis, sol_Xn);;
        }
    }
    free(rhs);
    free(rhsg);
}

/**
 * @brief   Solve poisson's equations construted by distinct pairs of basis
 */
void solve_halfpair_poissons_equations(SPARC_OBJ *pSPARC, 
    DDBP_ELEM *E_k, double *v, double *sol)
{
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    int Nde = E_k->nd_d;
    int Ndg = pSPARC->Nd;
    int Nxg = pSPARC->Nx;
    int Nyg = pSPARC->Ny;
    int Nzg = pSPARC->Nz;
    int basis_start_index = DDBP_info->basis_start_index;
    int n_basis = DDBP_info->n_basis_basiscomm;

    double *rhs = (double *) calloc(Nde, sizeof(double));
    assert(rhs != NULL);

    for (int n = 0; n < n_basis; n++) {
        double *v_n = v + n *Nde;
        double *rhsg = (double *) calloc((n_basis-n)*Ndg, sizeof(double));
        assert(rhsg != NULL);
        memset(rhsg, 0, sizeof(double)* (n_basis-n)*Ndg);

        for (int m = n; m < n_basis; m++) {
            double *v_m = v + m * Nde;
            for (int i = 0; i < Nde; i++) {
                rhs[i] = v_n[i] * v_m[i];
            }

            int DMnx = E_k->nx_d;
            int DMny = E_k->ny_d;
            int DMnz = E_k->nz_d;
            double *rhsg_ = rhsg + (m-n) *Ndg;
            // construct the rhs global vector 
            extract_subgrid(
                rhs, DMnx, DMny, DMnz,
                0, DMnx-1, 0, DMny-1, 0, DMnz-1, 
                rhsg_, Nxg, Nyg, Nzg,
                E_k->is, E_k->ie, E_k->js, E_k->je, E_k->ks, E_k->ke
            );
        }
        
        double *sol_Xn = sol + (basis_start_index + n) * n_basis * Ndg;
        double *sol_nn = sol_Xn + n * Ndg;        
        pois_fft(pSPARC, rhsg, pSPARC->pois_const, n_basis-n, sol_nn);
        free(rhsg);

        // copy off-diagonal solutions
        for (int m = n + 1; m < n_basis; m++) {
            double *sol_mn = sol_Xn + m * Ndg;
            double *sol_nm = sol + (basis_start_index + m) * n_basis * Ndg + n * Ndg;
            memcpy(sol_nm, sol_mn, sizeof(double)*Ndg);
        }
    }
    free(rhs);
}

/**
 * @brief   Gather half basis for solving poissons' equations
 */
void gather_half_basis(DDBP_EXX_POIS_PARAL_G *pois_paral_g, 
    int Nde, int ncol, double *v, double *recvbuf)
{
    int indegree = pois_paral_g->indegree;
    int *rcounts, *displs;
    rcounts = (int *) calloc(indegree, sizeof(int));
    displs = (int *) calloc(indegree, sizeof(int));
    assert(rcounts != NULL && displs != NULL);

    for (int i = 0; i < indegree; i++) {
        rcounts[i] = pois_paral_g->rcount[i] * Nde;
        displs[i] = pois_paral_g->displs[i] * Nde;
    }

    MPI_Comm dist_graph_comm = pois_paral_g->dist_graph_comm;

    if (dist_graph_comm != MPI_COMM_NULL) {
        MPI_Neighbor_allgatherv(v, ncol*Nde, MPI_DOUBLE, 
                        recvbuf, rcounts, displs, MPI_DOUBLE, dist_graph_comm);
    }

    free(rcounts);
    free(displs);
}

/**
 * @brief   Create communicator for parallelly solving poissons
 */
void Create_communicator_for_parallel_poisson(DDBP_INFO *DDBP_info)
{
    if (DDBP_info->elemcomm_index == -1 || 
        DDBP_info->basiscomm_index == -1 || 
        DDBP_info->dmcomm == MPI_COMM_NULL) {
        DDBP_EXX_POIS_PARAL_G *pois_paral_g = &DDBP_info->pois_paral_g;
        pois_paral_g->dist_graph_comm = MPI_COMM_NULL;
        DDBP_EXX_POIS_PARAL_D *pois_paral_d = &DDBP_info->pois_paral_d;
        pois_paral_d->dist_graph_comm = MPI_COMM_NULL;
        return;
    }

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    MPI_Comm blacscomm = DDBP_info->blacscomm;
    int blacscomm_rank, blacscomm_size;
    MPI_Comm_rank(blacscomm, &blacscomm_rank);
    MPI_Comm_size(blacscomm, &blacscomm_size);

    int Nb_k = DDBP_info->nALB_tot / DDBP_info->Ne_tot;
    int nsubcomm = DDBP_info->npbasis;
    int n_basis = DDBP_info->n_basis_basiscomm;
    int Ndg = DDBP_info->Nd_global;

    int basiscomm_index = DDBP_info->basiscomm_index;
    int color = (basiscomm_index == -1 || n_basis == 0) ? MPI_UNDEFINED : 0;
    MPI_Comm_split(blacscomm, color, basiscomm_index, &DDBP_info->blacscomm_NZ);

    DDBP_EXX_POIS_PARAL_G *pois_paral_g = &DDBP_info->pois_paral_g;
    int blacscomm_NZ_size = pois_paral_g->size;
    int indegree = pois_paral_g->indegree = (blacscomm_NZ_size+1)/2 - blacscomm_NZ_size%2;
    int outdegree = pois_paral_g->outdegree = pois_paral_g->indegree;
    pois_paral_g->sources = (int *) calloc(indegree, sizeof(int));
    pois_paral_g->destinations = (int *) calloc(outdegree, sizeof(int));
    assert(pois_paral_g->sources != NULL && pois_paral_g->destinations != NULL);
    
    for (int i = 0; i < indegree; i++) {
        pois_paral_g->sources[i] = (blacscomm_rank + i + 1) % blacscomm_NZ_size;
    }
    
    // order doesn't matter here
    for (int i = 0; i < outdegree; i++) {
        pois_paral_g->destinations[i] = (blacscomm_rank - outdegree + i + blacscomm_NZ_size) % blacscomm_NZ_size;
    }

    if (DDBP_info->blacscomm_NZ != MPI_COMM_NULL) {
        MPI_Dist_graph_create_adjacent(DDBP_info->blacscomm_NZ, 
            indegree, pois_paral_g->sources, MPI_UNWEIGHTED,
            outdegree, pois_paral_g->destinations, MPI_UNWEIGHTED,
            MPI_INFO_NULL, 0, &pois_paral_g->dist_graph_comm);
    } else {
        pois_paral_g->dist_graph_comm = MPI_COMM_NULL;
    }
    
    
    // create rcount and displs for neighbor gather
    pois_paral_g->rcount = (int *) calloc(indegree, sizeof(int));
    pois_paral_g->displs = (int *) calloc(indegree, sizeof(int));
    pois_paral_g->rcount_tot = 0;
    pois_paral_g->displs[0] = 0;
    for (int i = 0; i < indegree; i++) {
        pois_paral_g->rcount[i] = DDBP_info->n_basis_list[pois_paral_g->sources[i]];
        pois_paral_g->rcount_tot += pois_paral_g->rcount[i];
        if (i < indegree-1)
            pois_paral_g->displs[i+1] = pois_paral_g->displs[i] + pois_paral_g->rcount[i];
    }

    // creating communicator for distributing
    DDBP_EXX_POIS_PARAL_D *pois_paral_d = &DDBP_info->pois_paral_d;
    pois_paral_d->size = blacscomm_NZ_size;
    indegree = pois_paral_d->indegree = blacscomm_NZ_size - 1 - pois_paral_g->indegree;
    outdegree = pois_paral_d->outdegree = pois_paral_d->indegree;
    pois_paral_d->sources = (int *) calloc(indegree, sizeof(int));
    pois_paral_d->destinations = (int *) calloc(outdegree, sizeof(int));
    assert(pois_paral_d->sources != NULL && pois_paral_d->destinations != NULL);

    for (int i = 0; i < indegree; i++) {
        pois_paral_d->sources[i] = (blacscomm_rank - indegree + i + blacscomm_NZ_size) % blacscomm_NZ_size;
    }

    for (int i = 0; i < outdegree; i++) {
        pois_paral_d->destinations[i] = (blacscomm_rank + i + 1) % blacscomm_NZ_size;
    }

    if (DDBP_info->blacscomm_NZ != MPI_COMM_NULL) {
        MPI_Dist_graph_create_adjacent(DDBP_info->blacscomm_NZ, 
            indegree, pois_paral_d->sources, MPI_UNWEIGHTED,
            outdegree, pois_paral_d->destinations, MPI_UNWEIGHTED,
            MPI_INFO_NULL, 0, &pois_paral_d->dist_graph_comm);
    } else {
        pois_paral_d->dist_graph_comm = MPI_COMM_NULL;
    }
    
    // creating required arrays for distributing
    pois_paral_d->scount = (int *) calloc(outdegree, sizeof(int));
    pois_paral_d->sdispls = (int *) calloc(outdegree, sizeof(int));
    assert(pois_paral_d->scount != NULL && pois_paral_d->sdispls != NULL);
    pois_paral_d->scount_tot = 0;
    pois_paral_d->sdispls[0] = 0;
    for (int i = 0; i < outdegree; i++) {
        int dest = pois_paral_d->destinations[i];
        pois_paral_d->scount[i] = DDBP_info->n_basis_list[dest] * n_basis * Ndg;
        pois_paral_d->scount_tot += pois_paral_d->scount[i];
        if (i < outdegree-1)
            pois_paral_d->sdispls[i+1] = pois_paral_d->sdispls[i] + pois_paral_d->scount[i];
    }

    pois_paral_d->rcount = (int *) calloc(indegree, sizeof(int));
    pois_paral_d->rdispls = (int *) calloc(indegree, sizeof(int));
    assert(pois_paral_d->rcount != NULL && pois_paral_d->rdispls != NULL);
    pois_paral_d->rcount_tot = 0;
    pois_paral_d->rdispls[0] = 0;
    for (int i = 0; i < indegree; i++) {
        int src = pois_paral_d->sources[i];
        pois_paral_d->rcount[i] = DDBP_info->n_basis_list[src] * n_basis * Ndg;
        pois_paral_d->rcount_tot += pois_paral_d->rcount[i];
        if (i < indegree-1)
            pois_paral_d->rdispls[i+1] = pois_paral_d->rdispls[i] + pois_paral_d->rcount[i];
    }
}

/**
 * @brief   Free info for parallelly solving poissons
 */
void free_pois_paral_info(DDBP_INFO *DDBP_info)
{
    if (DDBP_info->elemcomm_index == -1 || 
        DDBP_info->basiscomm_index == -1 || 
        DDBP_info->dmcomm == MPI_COMM_NULL) {
        return;
    }

    if (DDBP_info->blacscomm_NZ != MPI_COMM_NULL)
        MPI_Comm_free(&DDBP_info->blacscomm_NZ);

    // free pois_paral_g structure
    DDBP_EXX_POIS_PARAL_G *pois_paral_g = &DDBP_info->pois_paral_g;
    if (pois_paral_g->dist_graph_comm != MPI_COMM_NULL)
        MPI_Comm_free(&pois_paral_g->dist_graph_comm);
    free(pois_paral_g->rcount);
    free(pois_paral_g->displs);
    free(pois_paral_g->sources);
    free(pois_paral_g->destinations);

    // free pois_paral_d structure
    DDBP_EXX_POIS_PARAL_D *pois_paral_d = &DDBP_info->pois_paral_d;
    if (pois_paral_d->dist_graph_comm != MPI_COMM_NULL)
        MPI_Comm_free(&pois_paral_d->dist_graph_comm);
    free(pois_paral_d->sources);
    free(pois_paral_d->destinations);
    free(pois_paral_d->scount);
    free(pois_paral_d->sdispls);
    free(pois_paral_d->rcount);
    free(pois_paral_d->rdispls);
}




/**
 * @brief   Create ACE operator in DDBP basis
 */
void ACE_operator_DDBP(SPARC_OBJ *pSPARC)
{
    DDBP_INFO *DDBP_info = pSPARC->DDBP_info;
    if (DDBP_info->elemcomm_index == -1 || DDBP_info->basiscomm_index == -1 ||
        DDBP_info->dmcomm == MPI_COMM_NULL) {
        return;
    }
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm blacscomm = DDBP_info->blacscomm;
    int nelem = DDBP_info->n_elem_elemcomm;
    int Ndg = DDBP_info->Nd_global;
    int nbasis = DDBP_info->n_basis_basiscomm;
    int Nband = DDBP_info->n_band_bandcomm;
    int Nstates = pSPARC->Nstates;
    int Netot = DDBP_info->Ne_tot;
    int spn_i = 0, kpt = 0;
    double t1, t2, t_QQG, t_gather_v, t_pois, t_reorg, t_VxViN2, t_VxViX, t_VxViXReduce, t_M, t_Choltrsm;
    t_QQG = t_gather_v = t_pois = t_reorg = t_VxViN2 = t_VxViX = t_VxViXReduce = t_M = t_Choltrsm = 0.0;

    double *VxViX_g = (double *) calloc(Ndg * Nband, sizeof(double));    
    assert(VxViX_g != NULL);

    // 3 conditions for skipping solving poissons' equations. 
    int saveFspvFlag = (pSPARC->fock_err < pSPARC->DDBP_fock_tol_updatebasis) && (DDBP_info->hybrid_memory_flag == 1);
    assert(saveFspvFlag == 1 || saveFspvFlag == 0);

#ifdef DEBUG
    if (!rank) printf("fock_err %e, DDBP_fock_tol_updatebasis %e, saveFspvFlag %d, vfsps_saved_flag %d\n",
        pSPARC->fock_err, pSPARC->DDBP_fock_tol_updatebasis, saveFspvFlag, DDBP_info->vfsps_saved_flag); fflush(stdout);
#endif 

    // here assuming nALB is the same for all elements!
    int nALB_ = DDBP_info->nALB_tot / DDBP_info->Ne_tot;
    double ***QQG = (double ***) calloc(nelem, sizeof(double **));
    assert(QQG != NULL);
    for (int k = 0; k <nelem; k++) {
        QQG[k] = (double **) calloc(Netot, sizeof(double *));
        for (int j = 0; j < Netot; j++) {
            QQG[k][j] = (double *) calloc(nALB_*nALB_, sizeof(double));
            assert(QQG[k][j] != NULL);
        }
        assert(QQG[k] != NULL);
    }

    t1 = MPI_Wtime();
    calculate_QQG(DDBP_info, &DDBP_info->xorb[spn_i][kpt], QQG);
    t2 = MPI_Wtime();
    t_QQG += (t2 - t1);

    t1 = MPI_Wtime();
    double **v_full = (double **) calloc(nelem, sizeof(double*));
    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        int Nde = E_k->nd_d;
        int nALB = E_k->nALB;
        assert(nALB == nALB_);

        v_full[k] = (double *) calloc(Nde*nALB, sizeof(double));
        assert(v_full[k] != NULL);

        gather_horizontally(
            Nde, nbasis, DDBP_info->basis_start_index, nALB, 
            DDBP_info->npbasis, E_k->v, v_full[k], blacscomm);
    }

    // here assume nALB is the same for all elements
    double *v_full_all = (double *) calloc(Ndg*nALB_, sizeof(double));
    assert(v_full_all != NULL);

    gather_vertically(DDBP_info, nALB_, v_full, v_full_all, 0);
    t2 = MPI_Wtime();
    t_gather_v += (t2 - t1);

    for (int k = 0; k < nelem; k++) {
        // t1 = MPI_Wtime();
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        int Nde = E_k->nd_d;
        int nALB = E_k->nALB;

        double *VxVi = (double *) calloc(Ndg * nbasis, sizeof(double));
        assert(VxVi != NULL);

        // The storage of total solutions is O(Nde*Ne*nALB*nALB*Ne) = O(Nde*Ne^2*nALB^2)
        // For large systems e.g. 100 elements, each element needs around 10Gb for solutions
        if (DDBP_info->vfsps_saved_flag == 1) {
            for (int i = 0; i < nALB; i++) {
                t1 = MPI_Wtime();
                computeVxViN2method(DDBP_info, NULL, NULL, QQG[k], VxVi, k, i, DDBP_info->vfsps[k][i]);
                t2 = MPI_Wtime();
                t_VxViN2 += (t2 - t1); 
            }
        } else {
            for (int i = 0; i < nALB; i++) {
                // storage of solution and element-wise arrangement solution
                double *sol = (double *) calloc(Ndg * nbasis, sizeof(double));
                double *sol_es = (double *) calloc(Ndg * nbasis, sizeof(double));
                assert(sol != NULL && sol_es != NULL);

                t1 = MPI_Wtime();
                double *v_full_i = v_full[k] + i * Nde;
                // solve poisson's equation of v_full_i * v
                poissonSolveDDBP(pSPARC, E_k, nbasis, E_k->v, v_full_i, sol);
                t2 = MPI_Wtime();
                t_pois += (t2 - t1);

                t1 = MPI_Wtime();
                BandWise2ElementWise_local(DDBP_info, nbasis, sol, sol_es);
                t2 = MPI_Wtime();
                t_reorg += (t2 - t1);
                free(sol);

                t1 = MPI_Wtime();
                double **vfsps = (DDBP_info->hybrid_memory_flag) ? DDBP_info->vfsps[k][i]: NULL;
                computeVxViN2method(DDBP_info, v_full_all, sol_es, QQG[k], VxVi, k, i, vfsps);
                t2 = MPI_Wtime();
                t_VxViN2 += (t2 - t1);
                free(sol_es);
            }
        }
        free(v_full[k]);

        DDBP_ARRAY *xorb = &DDBP_info->xorb[spn_i][kpt];
        t1 = MPI_Wtime();
        GlobalvecTimesX(DDBP_info, E_k, k, 1.0, VxVi, xorb, 1.0, VxViX_g);
        t2 = MPI_Wtime();
        t_VxViX += (t2 - t1);
        free(VxVi);
    }

    if (saveFspvFlag) DDBP_info->vfsps_saved_flag++;

    // free memory for QQG
    for (int k = 0; k <nelem; k++) {
        for (int j = 0; j < Netot; j++) {
            free(QQG[k][j]);
        }
        free(QQG[k]);
    }
    free(QQG);
    free(v_full);
    free(v_full_all);

    t1 = MPI_Wtime();
    ReduceExtractVxViX(DDBP_info, VxViX_g, NULL, 2);
    free(VxViX_g);
    t2 = MPI_Wtime();
    t_VxViXReduce += (t2 - t1);

    t1 = MPI_Wtime();
    int nrows_M, ncols_M;
    nrows_M = DDBP_info->nrows_M;
    ncols_M = DDBP_info->ncols_M;
    double *M = (double *) calloc(nrows_M*ncols_M, sizeof(double));
    assert(M != NULL);

    // compute local M = psi'*VxViX_g
    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        int Nde = E_k->nd_d;
        double *psi_k = DDBP_info->psi[spn_i][kpt][k];
        int *desc_psi_k = DDBP_info->desc_psi[k];
        double *Xi_k = DDBP_info->Xi[k];
        double alpha = 1.0, beta = 0.0;
        int ONE = 1;
        // M = xk^T * yk + M
        // generally, pdgemm_subcomm is not faster than pdgemm here.
        pdgemm_subcomm(
            "T", "N", Nstates, Nstates, Nde, 1.0, psi_k, desc_psi_k, 
            Xi_k, desc_psi_k, 1.0, M, DDBP_info->desc_M,
            blacscomm, best_max_nproc(Nstates, Nde, "pdgemm")
        );
        // pdgemm_(
        //     "T", "N", &Nstates, &Nstates, &Nde, 
        //     &alpha, psi_k, &ONE, &ONE, desc_psi_k, 
        //     Xi_k, &ONE, &ONE, desc_psi_k, 
        //     &beta, M, &ONE, &ONE, DDBP_info->desc_M);
    }
    t2 = MPI_Wtime();
    t_M += (t2 - t1);

    t1 = MPI_Wtime();
    MPI_Comm bandcomm = DDBP_info->bandcomm;
    int rank_bandcomm, size_bandcomm;
    MPI_Comm_rank(bandcomm, &rank_bandcomm);
    MPI_Comm_size(bandcomm, &size_bandcomm);

    if (size_bandcomm > 1)
        MPI_Allreduce(MPI_IN_PLACE, M, 
            nrows_M*ncols_M, MPI_DOUBLE, MPI_SUM, bandcomm);

    for (int i = 0; i < nrows_M*ncols_M; i++) M[i] *= -1.0;

#if defined(USE_MKL) || defined(USE_SCALAPACK)
    int ONE = 1, info = 0;
    double alpha = 1.0;
    // pdpotrf_("U", &Nstates, M, &ONE, &ONE, DDBP_info->desc_M, &info);
    if (nrows_M*ncols_M > 0) {
        info = LAPACKE_dpotrf (LAPACK_COL_MAJOR, 'U', Nstates, M, Nstates);
    }

#ifdef DEBUG
    if (!rank) printf("info = %d for LAPACKE_dpotrf.\n", info);
    fflush(stdout);
#endif
    
    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        int *desc_psi_k = DDBP_info->desc_psi[k];
        int Nde = E_k->nd_d;
        double *Xi_k = DDBP_info->Xi[k];
        pdtrsm_("R", "U", "N", "N", &Nde, &Nstates, &alpha, 
                M, &ONE, &ONE, DDBP_info->desc_M, 
                Xi_k, &ONE, &ONE, desc_psi_k);
    }
#endif 
    free(M);
    t2 = MPI_Wtime();
    t_Choltrsm += (t2 - t1);

#ifdef DEBUG
    if (!rank) printf("\n== Creating ACE operator == \n");
    if (!rank) printf("== Time for calculating QcwkpQG                    : %.3f ms\n", t_QQG*1e3);
    if (!rank) printf("== Time for gathering ddbp basis accross bandcomm  : %.3f ms\n", t_gather_v*1e3);
    if (!rank) printf("== Time for solving poisson's equations            : %.3f ms\n", t_pois*1e3);        
    if (!rank) printf("== Time for reorganizing poisson's solutions       : %.3f ms\n", t_reorg*1e3);
    if (!rank) printf("== Time for calculating VxVi and reduce            : %.3f ms\n", t_VxViN2*1e3);
    if (!rank) printf("== Time for calculating VxVi times DDBP orbitals   : %.3f ms\n", t_VxViX*1e3);
    if (!rank) printf("== Time for reducing VxVi and extract components   : %.3f ms\n", t_VxViXReduce*1e3);
    if (!rank) printf("== Time for creating M                             : %.3f ms\n", t_M*1e3);
    if (!rank) printf("== Time for choleskey and apply triag              : %.3f ms\n", t_Choltrsm*1e3);
#endif
}



void calculate_QQG(DDBP_INFO *DDBP_info, DDBP_ARRAY *X, double ***QQG)
{
    const havxX_t *havxX_info = &X->havxX_info;
    // in case that LANCZOS_DDBP proc with havxX_flag = 0 will enter here.
    if (havxX_info->havxX_flag == 0) return;
    assert(havxX_info->havxX_flag == 1);

    MPI_Comm colcomm = DDBP_info->bandcomm;
    MPI_Comm rowcomm = DDBP_info->rowcomm;
    int rowcomm_size, rowcomm_rank;
    MPI_Comm_size(rowcomm, &rowcomm_size);
    MPI_Comm_rank(rowcomm, &rowcomm_rank);
    
    int nelem = DDBP_info->n_elem_elemcomm;
    assert(nelem == X->nelem);
    int netot = DDBP_info->Ne_tot;
    // warning: only correct when nrows are the same 
    // for all elements (in current implementation)
    int blksz = havxX_info->blksz;
    int ncol = X->ncol;
    if (ncol == 0) {
        for (int k = 0; k < nelem; k++) {
            int nrows = X->nrows[k];
            int nrows2 = nrows*nrows;
            for (int j = 0; j < netot; j++) {
                if (rowcomm_size > 1) {
                    MPI_Allreduce(MPI_IN_PLACE, QQG[k][j], nrows2, MPI_DOUBLE, MPI_SUM, rowcomm);
                }
            }
        }
        return;
    }

    int nrows_total = havxX_info->nrows_total;
    double *buff = (double *) calloc(ncol*nrows_total, sizeof(double));
    assert(buff != NULL);
    
    // first gather all DDBP coefficients
    DDBP_element_allelements_allgatherv(X, buff, colcomm);
    
    double *occ_start = DDBP_info->occ_outer + DDBP_info->band_start_index;

    for (int k = 0; k < nelem; k++) {
        int nrows = X->nrows[k];
        int nrows2 = nrows*nrows;
        for (int j = 0; j < netot; j++) {
            double *xj = buff + j*blksz; 
            double *QQ = (double *) calloc(nrows2*ncol, sizeof(double));
            assert(QQ != NULL);
            cwkp(X->array[k], xj, nrows, nrows, ncol, QQ);
            // QQG = QQ * G
            // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
            //     nrows2, 1, ncol, 1.0, QQ, nrows2, occ_start, ncol,
            //     0.0, QQG[k][j], nrows2);
            cblas_dgemv(CblasColMajor, CblasNoTrans, nrows2, 
                ncol, 1.0, QQ, nrows2, 
                occ_start, 1, 0.0, QQG[k][j], 1);
            
            free(QQ);

            if (rowcomm_size > 1) {
                MPI_Allreduce(MPI_IN_PLACE, QQG[k][j], nrows2, MPI_DOUBLE, MPI_SUM, rowcomm);
                // replace with MPI_Iallreduce
                // int MPI_Iallreduce(const void* sendbuf, void* recvbuf, int count,
                //     MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                //     MPI_Request *request)
            }        
        }
    }

    free(buff);
}


void cwkp(double *A, double *B, int nrowA, int nrowB, int ncol, double *AkpB) 
{
#define A(i,j) A[i+j*nrowA]
#define B(i,j) B[i+j*nrowB]
#define AkpB(i,j) AkpB[i+j*nrowAkpB]

    int nrowAkpB = nrowA * nrowB;
    if (nrowAkpB == 0) return;
    for (int n = 0; n < ncol; n++) {
        int count = 0;
        // kron(A(:,n),B(:,n))
        for (int i = 0; i < nrowA; i++) {
            for (int j = 0; j < nrowB; j++) {                
                AkpB(count,n) = A(i,n) * B(j,n);
                count ++;
            }
        }
    }

#undef A
#undef B
#undef AkpB
}

void fsp(double *A, double *B, int nrow, int ncolA, int ncolB, double *AfspB)
{   
    int count = 0, idx_a = 0;
    for (int j = 0; j < ncolA; j++) {        
        int idx_b = 0;
        for (int k = 0; k < ncolB; k++) {            
            for (int i = 0; i < nrow; i++) {
                AfspB[count++] = A[idx_a+i] * B[idx_b+i];
            }
            idx_b += nrow;
        }
        idx_a += nrow;
    }
}


void fsp_omp(double *A, double *B, int nrow, int ncolA, int ncolB, double *AfspB) {
    int count = 0;
    // Parallelize the outer loop
    #pragma omp parallel for collapse(2) reduction(+:count)
    for (int j = 0; j < ncolA; j++) {
        for (int k = 0; k < ncolB; k++) {
            int idx_a = j * nrow;
            int idx_b = k * nrow;
            int idx_r_col = k + j * ncolB;
            for (int i = 0; i < nrow; i++) {
                int local_count = idx_r_col * nrow + i; // Compute index without contention
                AfspB[local_count] = A[idx_a + i] * B[idx_b + i];
            }
        }
    }
}

// replace gather_ESvec_and_E2B later 
void gather_vertically(DDBP_INFO *DDBP_info, 
    int ncol, double **vec_loc, double *vec_full, int option)
{
    // option 0, using vec_loc, option 1, using E_k->v (gather v)
    assert(option == 0 || option == 1);
    MPI_Comm comm = DDBP_info->bandcomm;
    assert(comm != MPI_COMM_NULL);
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    assert(size == DDBP_info->npelem);

    int Nd_global = DDBP_info->Nd_global;
    int nelem = DDBP_info->n_elem_elemcomm;
    int spn_i, kpt;
    spn_i = kpt = 0;

    for (int k = 0; k < nelem; k++) {
        DDBP_ELEM *E_k = &DDBP_info->elem_list[k];
        int global_index = k + DDBP_info->elem_start_index;
        int Nd_d = E_k->nd_d;
        double *src = (option == 0) ? vec_loc[k] : E_k->v;
        memcpy(vec_full + DDBP_info->Ndedisp[global_index]*ncol, 
               src, Nd_d*ncol*sizeof(double));
    }

    int *recvcounts, *displs;
    recvcounts = (int *) calloc(size, sizeof(int));
    displs = (int *) calloc(size, sizeof(int));
    assert(recvcounts != NULL && displs != NULL);
    
    displs[0] = 0;
    for (int i = 0; i < size; i++) {
        int Nd_d_total = 0;
        for (int k = 0; k < DDBP_info->Nelem_list[i]; k++) {
            int global_index = k + DDBP_info->elemstart_list[i];
            Nd_d_total += DDBP_info->Ndelist[global_index];
        }
        recvcounts[i] = Nd_d_total * ncol;
        if (i < size - 1)
            displs[i+1] = displs[i] + recvcounts[i];
    }

    // TODO: change to non-blocking
    MPI_Allgatherv(MPI_IN_PLACE, 1, MPI_DOUBLE, vec_full, recvcounts, displs, MPI_DOUBLE, comm);

    free(recvcounts);
    free(displs);
}


void computeVxViN2method(DDBP_INFO *DDBP_info, 
    double *v_full_all, double *sol_es, double **QQGk, double *VxVi, int element_index, int basis_index, double **vfsps_save)
{
    int Netot = DDBP_info->Ne_tot;
    int nbasis = DDBP_info->n_basis_basiscomm;
    // if (nbasis == 0) return;
    DDBP_ELEM *E_k = &DDBP_info->elem_list[element_index];
    int nALB = E_k->nALB;
    int Ndg = DDBP_info->Nd_global;
    int Nxg = DDBP_info->Nx_global;
    int Nyg = DDBP_info->Ny_global;
    int Nzg = DDBP_info->Nz_global;
    double dV = DDBP_info->dV;
    MPI_Comm blacscomm = DDBP_info->blacscomm;
    int rank, size;
    MPI_Comm_rank(blacscomm, &rank);
    MPI_Comm_size(blacscomm, &size);

    double *VxVil = (double *) calloc(Ndg*Netot, sizeof(double));    
    assert(VxVil != NULL);

    double t1, t2, t_gather, t_dgemv, t_reduce, t_VxVi, t_fsp;
    t_gather = t_dgemv = t_reduce = t_VxVi = t_fsp = 0;

    double *vi = v_full_all;
    double *si = sol_es;
    double *VxVil_ = VxVil;
    int Ndemax = DDBP_info->Ndemax;    

    for (int n = 0; n < Netot; n++) {
        int Nde = DDBP_info->Ndelist[n]; 

        t1 = MPI_Wtime();
        // compute local part of vfsps
        double *vfsps;
        if (DDBP_info->hybrid_memory_flag) {
            vfsps = vfsps_save[n];
        } else {
            vfsps = (double *) malloc(Nde*nbasis*nALB * sizeof(double));
        }
        assert(vfsps != NULL);

        if (DDBP_info->vfsps_saved_flag == 0) {
            fsp_omp(si, vi, Nde, nbasis, nALB, vfsps);
        }
        t2 = MPI_Wtime();
        t_fsp += (t2 - t1);

        t1 = MPI_Wtime();
        if (nbasis > 0) {
            // local dgemm of vfsps and QQG
            int nbasisdisp = DDBP_info->basis_start_index_list[rank];
            double *QQGkn_loc = QQGk[n] + nbasisdisp*nALB;
            // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
            //         Nde, 1, nALB*nbasis, -1.0/dV, vfsps, Nde, 
            //         QQGkn_loc, nALB*nbasis, 0.0, VxVil_, Nde);
            cblas_dgemv(CblasColMajor, CblasNoTrans, Nde, 
                nALB*nbasis, -1.0/dV, vfsps, Nde, 
                QQGkn_loc, 1, 0.0, VxVil_, 1);
        }

        t2 = MPI_Wtime();
        t_dgemv += (t2 - t1);

        if (!DDBP_info->hybrid_memory_flag) free(vfsps);
        VxVil_ += Nde;
        vi += nALB * Nde;
        si += nbasis * Nde;
    }

    t1 = MPI_Wtime();
    // find destination for reduce of VxVi
    int dest = basis_index_to_basiscomm_rank(nALB, size, basis_index);
    // reducing 
    if (dest == rank)
        MPI_Reduce(MPI_IN_PLACE, VxVil, Ndg, MPI_DOUBLE, MPI_SUM, dest, blacscomm);
    else
        MPI_Reduce(VxVil, VxVil, Ndg, MPI_DOUBLE, MPI_SUM, dest, blacscomm);

    t2 = MPI_Wtime();
    t_reduce += (t2 - t1);

    int shift = basis_index - DDBP_info->basis_start_index;
    if (dest == rank) assert(shift >= 0);
    
    t1 = MPI_Wtime();
    VxVil_ = VxVil;
    for (int n = 0; n < Netot; n++) {
        int Nde = DDBP_info->Ndelist[n];
        double *VxVi_ = VxVi + shift * Ndg;
        if (dest == rank) {
            // convert to column wise and save into correct position
            int *DMVertices = DDBP_info->DMvert_list[n];
            int DMnx = DMVertices[1] - DMVertices[0] + 1;
            int DMny = DMVertices[3] - DMVertices[2] + 1;
            int DMnz = DMVertices[5] - DMVertices[4] + 1;
            int DMnd = DDBP_info->Ndelist[n];
            assert(DMnd == DDBP_info->Ndelist[n]);
            extract_subgrid(
                VxVil_, DMnx, DMny, DMnz,
                0, DMnx-1, 0, DMny-1, 0, DMnz-1,
                VxVi_, Nxg, Nyg, Nzg,
                DMVertices[0], DMVertices[1], DMVertices[2],
                DMVertices[3], DMVertices[4], DMVertices[5]
            );
        }

        VxVil_ += Nde;
    }
    free(VxVil);

    t2 = MPI_Wtime();
    t_VxVi += (t2 - t1);

    // int grank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &grank);
    // if (!grank) printf("Time decomposition of computeVxViN2method:\n"
    //                   "t_fsp = %.3f ms\n"                      
    //                   "t_dgemv = %.3f ms\n"
    //                   "t_reduce = %.3f ms\n"
    //                   "t_VxVi = %.3f ms\n",
    //                   t_fsp*1e3, t_dgemv*1e3, t_reduce*1e3, t_VxVi*1e3);
}

