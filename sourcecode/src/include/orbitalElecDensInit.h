/**
 * @file    orbitalElecDensInit.h
 * @brief   This file declares the functions for electron density initialization.
 *
 * @authors Qimen Xu <qimenxu@gatech.edu>
 *          Abhiraj Sharma <asharma424@gatech.edu>
 *          Phanish Suryanarayana <phanish.suryanarayana@ce.gatech.edu>
 * 
 * Copyright (c) 2020 Material Physics & Mechanics Group, Georgia Tech.
 */

#ifndef ORBITALELECDENSINIT_H
#define ORBITALELECDENSINIT_H

#include "isddft.h"

/**
 * @brief Read vector(s) from Cube file(s) and distribute in comm.
 */
void read_cube_and_dist_vec(
    SPARC_OBJ *pSPARC, char filenames[3][L_STRING+L_PSD], double *data_dist, int n,
    int DMverts[6], MPI_Comm comm);

/**
 * @brief   Initialze electron density.
 */
void Init_electronDensity(SPARC_OBJ *pSPARC);

/*
 @ brief: function to perform charge extrapolation to provide better rho_guess for future relax/MD steps
*/
void elecDensExtrapolation(SPARC_OBJ *pSPARC);

/**
 * @brief   initialize Kohn-Sham orbitals.
 */
void Init_orbital(SPARC_OBJ *pSPARC);

#endif // ORBITALELECDENSINIT_H

