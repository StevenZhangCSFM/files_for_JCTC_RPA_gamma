#ifndef STERNEQKPT
#define STERNEQKPT

#include "isddft.h"
#include "main.h"

/**
 * @brief Solve all Sternheimer eq.s of trial vector dVs (from its nuChi0Eigscomm), k-points (from its kptcomm) and bands psis (from its bandcomm) assigned to this processor
 *
 * @param pSPARC pointer to SPARC_OBJ
 * @param pRPA pointer to RPA_OBJ
 * @param qptIndex the index of the q-point
 * @param omegaIndex the index of omega of the current nu chi0 operator
 * @param nuChi0EigsAmount amount of solved eigenvalues in this nuChi0Eigscomm, at here it is the amount of DV vectors
 * @param DVs (input) deltaVs in Sternheimer eq. (Hscf - epsilon_n + i omega) = - deltaV*psi_n
 * @param printFlag flag to print mid variables
 */
void sternheimer_eq_kpt(SPARC_OBJ *pSPARC, RPA_OBJ *pRPA, int qptIndex, int omegaIndex, int nuChi0EigsAmount, double _Complex *DVs, int printFlag);

void sternheimer_solver_kpt(SPARC_OBJ *pSPARC, int spn_i, int kPqIndex, int kptIndex, double epsilon, double omega, int flagPQ,
    int flagCOCGinitial, double _Complex *allXorb_thekPq, double *allLambdas_thekPq, int NusedOrbitalForInitGuess, 
    int flagPreconditioner, double *poisC_invE_precond, double Icoeff,
    double _Complex *deltaPsisPlusOmega, double _Complex *deltaPsisMinusOmega, double _Complex *deltaVs_kpt, double _Complex *psi_kpt, double bandWeight, double _Complex *deltaRhos, int nuChi0EigsAmounts,
    double sternRelativeResTol, int SternBlockSize, int printFlag, char *outputName, double *timeRecorder);

/**
 * @brief obtain the product of the LHS of the Sternheimer equation and any set of complex vectors X, 
 *        LHS at here is (K-S Hamiltonian - epsilon_n + Q, if (flagPQ), |psi_n><psi_n| + i omega)X
 *
 * @param pSPARC pointer to SPARC_OBJ
 * @param spn_i spin of the input orbital psi_n
 * @param kPq index of bloch vector k+q in the complete k-point list
 * @param epsilon (input) epsilon_n in Sternheimer eq. (Hscf - epsilon_n + i omega) = - deltaV*psi_n
 * @param omega (input) omega in Sternheimer eq. (Hscf - epsilon_n + i omega) = - deltaV*psi_n
 * @param flagPQ the flag to use linear operator P (|psi_n><psi_n| / dV) and Q (|psi_n><psi_n|)
 * @param X (input) the vectors to be left multiplied by LHS
 * @param lhsX (output) the product of the LHS of the Sternheimer equation and complex vectors X
 * @param nuChi0EigsAmounts the amounts of vectors in X
 */
void Sternheimer_lhs_kpt(SPARC_OBJ* pSPARC, int spn_i, int kPq, double epsilon, double omega, double _Complex *X, double _Complex *lhsX, int nuChi0EigsAmounts);

/**
 * @brief compute the initial guess solutions of the Sternheimer equations based on the eigenpairs of Hamiltonian
 *
 * @param pSPARC pointer to SPARC_OBJ
 * @param epsilon (input) epsilon_n in Sternheimer eq. (Hscf - epsilon_n +- i omega) = - deltaV*psi_n
 * @param omega (input) omega in Sternheimer eq. (Hscf - epsilon_n +- i omega) = - deltaV*psi_n
 * @param flagPQ the flag to use linear operator P (|psi_n><psi_n| / dV) and Q (|psi_n><psi_n|)
 * @param flagCOCGinitial the flag to use the initial guess solutions of Sternheimer equations. If it is closed, the initial guesses are set as zero vectors
 * @param NusedOrbital (input) the number of first several orbitals to be used to generate the initial guess solutions. The more orbitals used, the better the initial guesses.
 * @param allXorb_thekPq all of the solved orbitals of the block vector k+q from the Kohn-Sham DFT calculation
 * @param allLambdas_thekPq all of the solved eigenvalues of the block vector k+q from the Kohn-Sham DFT calculation
 * @param SternheimerRhs (input) all RHS of the Sternheimer equation, -psi_n .* [dV1 dV2 ... dVn]
 * @param deltaPsisPlusOmega (output) the initial guesses for Sternheimer equations having +i omega
 * @param deltaPsisMinusOmega (output) the initial guesses for Sternheimer equations having -i omega
 * @param nuChi0EigsAmounts the amounts of vectors in RHS
 */
void set_initial_guess_deltaPsis_kpt(SPARC_OBJ *pSPARC, double epsilon, double omega, int flagPQ, int flagCOCGinitial, int NusedOrbital, double _Complex *allXorb_thekPq, double *allLambdas_thekPq,
    double _Complex *SternheimerRhs, double _Complex *deltaPsisPlusOmega, double _Complex *deltaPsisMinusOmega, int nuChi0EigsAmounts);

#endif