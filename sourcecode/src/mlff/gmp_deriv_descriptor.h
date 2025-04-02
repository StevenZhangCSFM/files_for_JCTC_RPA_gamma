#ifndef GMP_DERIV_DESCRIPTOR_H
#define GMP_DERIV_DESCRIPTOR_H

#include "mlff_types.h"
#include "isddft.h"
#include "solid_harmonics.h"

void build_gmpObj(DescriptorObj *gmp_str, NeighList *nlist, int nmcsh, double *atompos, int **params_i, double **params_d, double **atom_gaussian, int *ngaussians, int *element_index_to_order, int *atom_type_to_indices, int *atom_indices_p);

void free_GMP(DescriptorObj *gmp_str);

double get_scaled_basis(int order, double sigma);

void create_logspace_list(double radial_min, double radial_max, int nsig, double *sigmas);

void read_GMP_inpt(MLFF_Obj *mlff_str, SPARC_OBJ *pSPARC);

void calc_mcsh_type_1(
    int ii,                            // Index of the current atom
    int m,                             // Index of the current descriptor
    int i,                             // Atom index in the atom_indices array
    double xi, double yi, double zi,   // Coordinates of the current atom
    NeighList *nlist,            // Neighbor list structure
    double *atompos,             // Array of atom positions
    int *atom_indices,           // Array mapping atom indices
    int *element_index_to_order, // Mapping from element index to order
    int *atom_type_to_indices,   // Mapping from atom types to indices
    double L1, double L2, double L3,   // Lattice parameters
    int cell_typ,                      // Cell type indicator
    double twist,                      // Twist parameter for non-orthogonal cells
    double A,                          // Parameter A for the Gaussian function
    double alpha,                      // Parameter alpha for the Gaussian function
    int *ngaussians,             // Number of Gaussians per element
    double **atom_gaussian,            // Gaussian parameters for each element
    int mcsh_order, int group_index, // MCSH order and group index
    double group_coefficient,          // Coefficient for the group
    DescriptorObj *gmp_str,             // Structure to store results
    double *sum_square
);

void calc_mcsh_type_2(
    int ii,                            // Index of the current atom
    int m,                             // Index of the current descriptor
    int i,                             // Atom index in the atom_indices array
    double xi, double yi, double zi,   // Coordinates of the current atom
    NeighList *nlist,            // Neighbor list structure
    double *atompos,             // Array of atom positions
    int *atom_indices,           // Array mapping atom indices
    int *element_index_to_order, // Mapping from element index to order
    int *atom_type_to_indices,   // Mapping from atom types to indices
    double L1, double L2, double L3,   // Lattice parameters
    int cell_typ,                      // Cell type indicator
    double twist,                      // Twist parameter for non-orthogonal cells
    double A,                          // Parameter A for the Gaussian function
    double alpha,                      // Parameter alpha for the Gaussian function
    int *ngaussians,             // Number of Gaussians per element
    double **atom_gaussian,            // Gaussian parameters for each element
    int mcsh_order, int group_index, // MCSH order and group index
    double group_coefficient,          // Coefficient for the group
    DescriptorObj *gmp_str,             // Structure to store results
    double *sum_square
);

void calc_mcsh_type_3(
    int ii,                            // Index of the current atom
    int m,                             // Index of the current descriptor
    int i,                             // Atom index in the atom_indices array
    double xi, double yi, double zi,   // Coordinates of the current atom
    NeighList *nlist,            // Neighbor list structure
    double *atompos,             // Array of atom positions
    int *atom_indices,           // Array mapping atom indices
    int *element_index_to_order, // Mapping from element index to order
    int *atom_type_to_indices,   // Mapping from atom types to indices
    double L1, double L2, double L3,   // Lattice parameters
    int cell_typ,                      // Cell type indicator
    double twist,                      // Twist parameter for non-orthogonal cells
    double A,                          // Parameter A for the Gaussian function
    double alpha,                      // Parameter alpha for the Gaussian function
    int *ngaussians,             // Number of Gaussians per element
    double **atom_gaussian,            // Gaussian parameters for each element
    int mcsh_order, int group_index, // MCSH order and group index
    double group_coefficient,          // Coefficient for the group
    DescriptorObj *gmp_str,             // Structure to store results
    double *sum_square
);

#endif