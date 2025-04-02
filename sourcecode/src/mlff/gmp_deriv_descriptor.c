#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>
#include "tools.h"
#include "helper_gmp.h"
#include "solid_harmonics.h"
#include "gmp_deriv_descriptor.h"
#include "tools_mlff.h"
#include "ddbp_tools.h"
#include "soap_descriptor.h"
#include "sparc_mlff_interface.h"

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

#define L_ATOMTYP 8
/*
delete_GMPObj function frees the memory allocated for dynamic arrays in GMPObj structure

[Input]
1. gmp_str: pointer to GMPObj structure

[Output]
1. pointer to deallocated GMPObj structure
*/

void free_GMP(DescriptorObj *gmp_str){

    for (int i = 0; i < gmp_str->natom_domain; i++){
        free(gmp_str->X3[i]);
        int uniq_natms = gmp_str->unique_Nneighbors[i];
        if (uniq_natms < gmp_str->natom) {
            uniq_natms ++;
        }
        for (int j = 0; j < uniq_natms; j++){
            free(gmp_str->dX3_dX[i][j]);
            free(gmp_str->dX3_dY[i][j]);
            free(gmp_str->dX3_dZ[i][j]);
        }
        for(int j=0; j < 6; j++){
            free(gmp_str->dX3_dF[i][j]);
        }
        free(gmp_str->dX3_dX[i]);
        free(gmp_str->dX3_dY[i]);
        free(gmp_str->dX3_dZ[i]);
        free(gmp_str->dX3_dF[i]);
    }
    free(gmp_str->X3);
    free(gmp_str->dX3_dF);
    free(gmp_str->dX3_dX);
    free(gmp_str->dX3_dY);
    free(gmp_str->dX3_dZ);

}


/*
build_gmpObj function calculates the GMP descriptor and its derivatives for each atom on each processor

[Input]
1. nlist: pointer to Neighlist structure
2. nmcsh: total number of descriptors/atom
3. atompos: pointer to atomic positions
4. integer GMP params: mcsh order, square (1) or square root (0), solid (1) or surface (0) harmonic
5. double GMP params: gaussian probe sigma, 1 (unexplained, never referenced - will change later), sigma manipulation, sigma manipulation, sigma manipulation
6. atom_gaussian: pointer to primitive gaussian parameters
7. ngaussians: pointer to number of primitives
8. element_index_to_order: param to convert atomic number to order in list of gaussian parameters
9. atom_type_to_indices: param to convert atomic number to index in list of gaussian parameters
10. atom_indices: param to convert index in nlist to atomic number

[Output]
1. gmp_str: pointer to GMPObj structure
*/

void build_gmpObj(DescriptorObj *gmp_str, NeighList *nlist, int nmcsh, double *atompos, int **params_i, double **params_d, double** atom_gaussian, int* ngaussians, int* element_index_to_order, int* atom_type_to_indices, int* atom_indices){

    double xi, yi, zi;
    int cell_typ = nlist->cell_typ;
    double twist = nlist->twist;
    double L1 = nlist->cell_len[0];
	double L2 = nlist->cell_len[1];
	double L3 = nlist->cell_len[2];
    // Reshape stress indices
    int index[6] = {0,1,2,3,4,5};
	reshape_stress(cell_typ,nlist->BC,index);
    for (int ii=0; ii < nlist->natom_domain; ++ii) {
        // i is the index of atom being calculated - this logic is a bit more complicated than it needs to be.
        //Index just loops from 0 - natoms - 1. Just be consistent and everything should work.    
        int i = nlist->atom_idx_domain[ii];
        xi = atompos[3*i];
        yi = atompos[3*i+1];
        zi = atompos[3*i+2];        
        for (int m = 0; m < nmcsh; ++m) {
            int mcsh_order = params_i[m][0];
            int num_groups = get_num_groups(mcsh_order);
            // Reminder, A, alpha terms are from probe gaussian.
            double A = params_d[m][2], alpha = params_d[m][3];
            double sum_square = 0.0;
            for (int group_index = 1; group_index < (num_groups+1); ++group_index){
                double group_coefficient = get_group_coefficients(mcsh_order, group_index);
                int mcsh_type = get_mcsh_type(mcsh_order, group_index);
                switch (mcsh_type){
                    case 1:
                        calc_mcsh_type_1(ii, m, i, xi, yi, zi, nlist, atompos, atom_indices, element_index_to_order, atom_type_to_indices, L1, L2, L3, cell_typ, twist, A, alpha, ngaussians, atom_gaussian, mcsh_order, group_index, group_coefficient, gmp_str, &sum_square);
                        break;
                    case 2:
                        calc_mcsh_type_2(ii, m, i, xi, yi, zi, nlist, atompos, atom_indices, element_index_to_order, atom_type_to_indices, L1, L2, L3, cell_typ, twist, A, alpha, ngaussians, atom_gaussian, mcsh_order, group_index, group_coefficient, gmp_str, &sum_square);
                        break;
                    case 3:
                        calc_mcsh_type_3(ii, m, i, xi, yi, zi, nlist, atompos, atom_indices, element_index_to_order, atom_type_to_indices, L1, L2, L3, cell_typ, twist, A, alpha, ngaussians, atom_gaussian, mcsh_order, group_index, group_coefficient, gmp_str, &sum_square);
                        break;
                    default:
                        printf("MCSH type not implemented\n");
                        exit(1);
                }
            }
            gmp_str->X3[ii][m] = sum_square;
        }
    }
}

// The final result depends on units, but this shouldn't since we don't introduce units to the gaussian integrals until we plug in
//alpha here. Everything is just an integral over all space.
double get_scaled_basis(int order, double sigma){
    double alpha = 1/(2*sigma*sigma);
    double _const = pow((M_PI / (2 * alpha)), 1.5);
    double temp;
    if (order == 0){
        temp = _const;
        return 1/sqrt(temp);
    } else if (order == 1){
        temp = 0.75 * pow( alpha, -1 ) * _const;
        return 1/sqrt(temp);
    } else if (order == 2){
        temp = 5.625 * pow( alpha, -2 ) * _const;
        return 1/sqrt(temp);
    } else if (order == 3){
        temp = 147.65625 * pow( alpha, -3 ) * _const;
        return 1/sqrt(temp);
    } else if (order == 4){
        temp = 9302.34375 * pow( alpha, -4 ) * _const;
        return 1/sqrt(temp);
    } else if (order == 5){
        temp = 1151165.0390625 * pow( alpha, -5 ) * _const;
        return 1/sqrt(temp);
    } else if (order == 6){
        temp = 246924900.87890625 * pow( alpha, -6 ) * _const;
        return 1/sqrt(temp);
    } else if (order == 7){
        temp = 84263122424.92676 * pow( alpha, -7 ) * _const;
        return 1/sqrt(temp);
    } else if (order == 8){
        temp = 42974192436712.65 * pow( alpha, -8 ) * _const;
        return 1/sqrt(temp);
    } else if (order == 9){
        temp = 3.163831150894136e+16 * pow( alpha, -9 ) * _const;
        return 1/sqrt(temp);
    } else{
        printf("Order not implemented\n");
        exit(1);
    }
}

void read_GMP_inpt(MLFF_Obj *mlff_str, SPARC_OBJ *pSPARC){

	char *str            = malloc(L_STRING * sizeof(char));
    int *ngaussians;
	
	ngaussians = (int*) malloc(sizeof(int)*pSPARC->Ntypes);
    int nsig, nang, Ntypes, _tmp;
	nsig = pSPARC->N_max_SOAP;
	nang = pSPARC->L_max_SOAP+1;
    Ntypes = pSPARC->Ntypes;
    double **atom_gaussians_p = (double**) malloc(sizeof(double*)*Ntypes);

	FILE *input_fp;
	// TODO: Update how to handle transfer learning
	/*
	if (pSPARC->transfer_learning) {
		input_fp = fopen("pretrained/GMP_input.in", "r");
	} else {
		input_fp = fopen(pSPARC->hnl_file_name, "r");
	}
	*/
	input_fp = fopen(pSPARC->hnl_file_name, "r");
	if (NULL == input_fp) {
		printf("GMP input file cannot be opened\n");
		exit(1);
	}

	/* TODO: Find way to read ASE_SORT
	else if (strcmpi(str,"ASE_SORT:") == 0) {
		for (int i=0; i<pSPARC->n_atom; i++){
			fscanf(input_fp,"%d",&mlff_str->ase_ordering[i]);
		}
		fscanf(input_fp, "%*[^\n]\n");
	// Read max number of prmitives (could use built in function, but just an input for now) - int
	}
	*/

    for (int i = 0; i < Ntypes; i++){
        rewind(input_fp);
        while (!feof(input_fp)){
            int count = fscanf(input_fp,"%s",str);
            if (count < 0) continue;  // for some specific cases
            
            // enable commenting with '#'
            if (str[0] == '#' || str[0] == '!' || str[0] == '\n'|| strcmpi(str,"undefined") == 0) {
                fscanf(input_fp, "%*[^\n]\n"); // skip current line
                continue;
            }
            // Using char* may require doctoring
            if (strcmpi(str, &pSPARC->atomType[L_ATOMTYP*i]) == 0){
                fscanf(input_fp,"%d %d",&_tmp,&ngaussians[i]);
                fscanf(input_fp, "%*[^\n]\n");
                atom_gaussians_p[i] = (double*) malloc(ngaussians[i]*2*sizeof(double));
                for (int j=0; j<ngaussians[i]; j++){
                    fscanf(input_fp,"%lf %lf",&atom_gaussians_p[i][j*2],&atom_gaussians_p[i][j*2+1]);
                }
            }
		}
    }
    fclose(input_fp);

    // Store cutoff
    double rcut = pSPARC->rcut_SOAP;
    // Get descriptor vector size from nsig and nang probes
    int nmcsh = nsig*nang;
    
    // Create integer and double parameter arrays. D1 = nmcsh, D2 = 6 for double, 3 for int
    // Allocate memory for GMP double params and populate 2D array
    mlff_str->params_d = (double**) malloc(sizeof(double*)*nmcsh);
    for (int i = 0; i < nmcsh; i++) mlff_str->params_d[i] = (double*) malloc(sizeof(double)*6); // Hardcoding 6 is fine - this does not change
    for (int i = 0; i < nang; i++){
        for (int j = 0; j < nsig; j++){
            mlff_str->params_d[i*nsig+j][0] = mlff_str->sigmas[j];
            mlff_str->params_d[i*nsig+j][1] = 1.0;
            // Normalization - change depending on test or production
			/*
			if (pSPARC->gmp_tests == 1) {
            	mlff_str->params_d[i*nsig+j][2] = pow(1/(sigmas[j]*sqrt(2.0*M_PI)),3);
			} else {
            	mlff_str->params_d[i*nsig+j][2] = get_scaled_basis(i, sigmas[j]);
			}
			*/
			//mlff_str->params_d[i*nsig+j][2] = get_scaled_basis(i, mlff_str->sigmas[j]);
            mlff_str->params_d[i*nsig+j][2] = pow(1/(mlff_str->sigmas[j]*sqrt(2.0*M_PI)),3);
            // Alpha
            mlff_str->params_d[i*nsig+j][3] = 1/(2*pow(mlff_str->sigmas[j],2));
            mlff_str->params_d[i*nsig+j][4] = rcut;
            mlff_str->params_d[i*nsig+j][5] = 1.0;
        }
    }

    // Allocate memory for GMP integer params and populate 2D array
	// 0: angular order, 1: squared (1) or square root (0), 2: solid (1) or surface (0)
    mlff_str->params_i = (int**) malloc(sizeof(int*)*nmcsh);
    for (int i=0; i<nmcsh; i++) mlff_str->params_i[i] = (int*) malloc(sizeof(int)*3);  // Hardcoding 3 is fine - this does not change
    for (int i = 0; i < nang; i++){
        for (int j = 0; j < nsig; j++){
            mlff_str->params_i[i*nsig+j][0] = i;
            mlff_str->params_i[i*nsig+j][1] = 1;
            mlff_str->params_i[i*nsig+j][2] = 1;
        }
    }

    // Allocate memory for gaussian primitive parameters and populate array. 2D array determined by Ntypes and max number of primitives
    // Primitives read from file which at current must be in working directory
    mlff_str->atom_gaussians_p = (double**) malloc(sizeof(double*)*Ntypes);
    for (int i=0; i < Ntypes; i++) mlff_str->atom_gaussians_p[i] = (double*) malloc(sizeof(double)*ngaussians[i]*2);
    for (int i=0; i < Ntypes; i++){
    	for (int j=0; j < 2*ngaussians[i]; j++){
    		mlff_str->atom_gaussians_p[i][j] = atom_gaussians_p[i][j];
    	}
    }

    mlff_str->ngaussians = (int*) malloc(sizeof(int)*Ntypes);
    for (int i = 0; i < Ntypes; i++){
    	mlff_str->ngaussians[i] = ngaussians[i]; 
    }

	free(str);
	free(ngaussians);
	for (int i=0; i<Ntypes; i++) free(atom_gaussians_p[i]);
	free(atom_gaussians_p);
}

void create_logspace_list(double radial_min, double radial_max, int nsig, double *sigmas){
    if (nsig <= 0 || radial_min <= 0 || radial_max <= 0 || radial_min > radial_max) {
        fprintf(stderr, "Invalid input parameters to logspace: %d %lf %lf\n", nsig, radial_min, radial_max);
        exit(1);
    }

    // Calculate the factor between consecutive values
    double log_min = log10(radial_min);
    double log_max = log10(radial_max);
    double step = (log_max - log_min) / (nsig - 1);

    for (int i = 0; i < nsig; i++) {
        sigmas[i] = pow(10, log_min + i * step);
    }
}

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
) {
    SolidGMPFunction mcsh_solid_function = get_solid_mcsh_function(mcsh_order, group_index);
    double sum_desc = 0.0;
    double *sum_dmiu_dxj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double *sum_dmiu_dyj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double *sum_dmiu_dzj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double *sum_dmiu_dFj = (double *) calloc(6, sizeof(double));
    double m_desc[1], deriv[3], xj, yj, zj, dx = 0.0, dy = 0.0, dz = 0.0, r0_sqr;
    int _count = 0;
    int *_indices = (int*) malloc(sizeof(int)*nlist->Nneighbors[ii]);
    int _tmp = -1, idx_neigh, elem_typ;
    // -1 accounts for the fact that atom i in the fundamental unit cell is not included. However, it's mirror images are included.
    //for (int j = -1; j < nneigh[ii]; ++j) {
    for (int j = -1; j < nlist->Nneighbors[ii]; ++j) {
        int neigh_atom_element_order;
        if (j < 0) {
            xj = xi;
            yj = yi;
            zj = zi;
            // gives the atomic number of the atom for a given index of atom - correponds to index in .ion file. 
            int neigh_atom_element_index = atom_indices[i];
            // Gives the order a given element appears in the .ion file and the primitive gaussians arrays
            neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
            _tmp = i;
        } else {
            idx_neigh = (nlist->neighborList + ii)->array[j];
            elem_typ = (nlist->neighborAtmTyp + ii)->array[j];
            xj = atompos[3*idx_neigh] + L1 * (nlist->neighborList_imgX + ii)->array[j];
            yj = atompos[3*idx_neigh+1] + L2 * (nlist->neighborList_imgY + ii)->array[j];
            zj = atompos[3*idx_neigh+2] + L3 * (nlist->neighborList_imgZ + ii)->array[j];
            int neigh_atom_element_index = atom_type_to_indices[elem_typ];
            neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
            //This assumes neighbor list is constructed sequentially, which it is. Basically, store deriv wrt self at 0,
            //then starting from 1, store the deriv of X wrt pos of neighbor 1 (idx in nlist of [0]). Keep adding to pos 1 until
            //a new atom is encountered. Then, start adding to pos 2, etc. Anytime the other occurences of the self atom are encountered,
            //set the index to 0. Then, next time a new atom is encountered continue incrementing. 
            if (idx_neigh == i) {
                _indices[j] = 0;
            } else if (_tmp < 0) {
                perror("Error in GMP descriptor calculation: indices loop broken\n");
            } else if (_tmp != idx_neigh) {
                _count += 1;
                _tmp = idx_neigh;
                _indices[j] = _count;
            } else {
                _indices[j] = _count;
            }
        }
        if (cell_typ > 20 && cell_typ < 30) { // add for nonorthogonal
            dx = xj*cos(yj+twist*zj) - xi*cos(yi+twist*zi);
            dy = xj*sin(yj+twist*zj) - xi*sin(yi+twist*zi);
            dz = zj - zi;
        } else if (cell_typ > 10 && cell_typ < 20){
            double dx1 = xj - xi;
            double dx2 = yj - yi;
            double dx3 = zj - zi;
            dx = nlist->LatUVec[0] * dx1 + nlist->LatUVec[3] * dx2 + nlist->LatUVec[6] * dx3;
            dy = nlist->LatUVec[1] * dx1 + nlist->LatUVec[4] * dx2 + nlist->LatUVec[7] * dx3;
            dz = nlist->LatUVec[2] * dx1 + nlist->LatUVec[5] * dx2 + nlist->LatUVec[8] * dx3;
        } else if (cell_typ == 0) {
            dx = xj - xi;
            dy = yj - yi;
            dz = zj - zi;
        }
        r0_sqr = dx*dx + dy*dy + dz*dz;
        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
            mcsh_solid_function(dx, dy, dz, r0_sqr, A, B, alpha, beta, m_desc, deriv);
            sum_desc += m_desc[0];
            sum_dmiu_dxj[j+1] += deriv[0];
            sum_dmiu_dyj[j+1] += deriv[1];
            sum_dmiu_dzj[j+1] += deriv[2];
        }
        sum_dmiu_dFj[0] += sum_dmiu_dxj[j+1] * dx;
        sum_dmiu_dFj[1] += sum_dmiu_dxj[j+1] * dy;
        sum_dmiu_dFj[2] += sum_dmiu_dxj[j+1] * dz;
        sum_dmiu_dFj[3] += sum_dmiu_dyj[j+1] * dy;
        sum_dmiu_dFj[4] += sum_dmiu_dyj[j+1] * dz;
        sum_dmiu_dFj[5] += sum_dmiu_dzj[j+1] * dz;
    }
    *sum_square += group_coefficient * sum_desc * sum_desc;
    for (int j = 0; j < 6; j++) {
        sum_dmiu_dFj[j] *= sum_desc * group_coefficient * 2.0;
        gmp_str->dX3_dF[ii][j][m] += sum_dmiu_dFj[j];
    }
    double dmdx, dmdy, dmdz;
    for (int j = -1; j < nlist->Nneighbors[ii]; ++j) {
        dmdx = (sum_desc * sum_dmiu_dxj[j+1]) * group_coefficient * 2.0;
        dmdy = (sum_desc * sum_dmiu_dyj[j+1]) * group_coefficient * 2.0;
        dmdz = (sum_desc * sum_dmiu_dzj[j+1]) * group_coefficient * 2.0;
        if (j < 0) {
            gmp_str->dX3_dX[ii][0][m] += dmdx;
            gmp_str->dX3_dY[ii][0][m] += dmdy;
            gmp_str->dX3_dZ[ii][0][m] += dmdz;
        }
        else {
            gmp_str->dX3_dX[ii][_indices[j]][m] += dmdx;
            gmp_str->dX3_dY[ii][_indices[j]][m] += dmdy;
            gmp_str->dX3_dZ[ii][_indices[j]][m] += dmdz;
        }
        gmp_str->dX3_dX[ii][0][m] -= dmdx;
        gmp_str->dX3_dY[ii][0][m] -= dmdy;
        gmp_str->dX3_dZ[ii][0][m] -= dmdz;						
    }
    free(sum_dmiu_dxj);
    free(sum_dmiu_dyj);
    free(sum_dmiu_dzj);
    free(sum_dmiu_dFj);
    free(_indices);
}

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
) {
    SolidGMPFunction mcsh_solid_function = get_solid_mcsh_function(mcsh_order, group_index);
    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0;
    double* sum_dmiu1_dxj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu2_dxj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu3_dxj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu1_dyj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu2_dyj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu3_dyj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu1_dzj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu2_dzj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu3_dzj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu1_dFj = (double *) calloc(6, sizeof(double));
    double* sum_dmiu2_dFj = (double *) calloc(6, sizeof(double));
    double* sum_dmiu3_dFj = (double *) calloc(6, sizeof(double));
    double miu[3], deriv[9], xj, yj, zj, dx = 0.0, dy = 0.0, dz = 0.0, r0_sqr;
    int _count = 0;
    int *_indices = (int*) malloc(sizeof(int)*nlist->Nneighbors[ii]);
    int _tmp = -1, idx_neigh, elem_typ;
    // -1 accounts for the fact that atom i in the fundamental unit cell is not included. However, it's mirror images are included.
    for (int j = -1; j < nlist->Nneighbors[ii]; ++j) {
        int neigh_atom_element_order;
        if (j < 0) {
            xj = xi;
            yj = yi;
            zj = zi;
            // gives the atomic number of the atom for a given index of atom - correponds to index in .ion file. 
            int neigh_atom_element_index = atom_indices[i];
            // Gives the order a given element appears in the .ion file and the primitive gaussians arrays
            neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
            _tmp = i;
        } else {
            idx_neigh = (nlist->neighborList + ii)->array[j];
            elem_typ = (nlist->neighborAtmTyp + ii)->array[j];
            xj = atompos[3*idx_neigh] + L1 * (nlist->neighborList_imgX + ii)->array[j];
            yj = atompos[3*idx_neigh+1] + L2 * (nlist->neighborList_imgY + ii)->array[j];
            zj = atompos[3*idx_neigh+2] + L3 * (nlist->neighborList_imgZ + ii)->array[j];
            int neigh_atom_element_index = atom_type_to_indices[elem_typ];
            neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
            if (idx_neigh == i) {
                _indices[j] = 0;
            } else if (_tmp < 0) {
                perror("Error in GMP descriptor calculation: indices loop broken\n");
            } else if (_tmp != idx_neigh) {
                _count += 1;
                _tmp = idx_neigh;
                _indices[j] = _count;
            } else {
                _indices[j] = _count;
            }
        }
        if (cell_typ > 20 && cell_typ < 30) { // add for nonorthogonal
            dx = xj*cos(yj+twist*zj) - xi*cos(yi+twist*zi);
            dy = xj*sin(yj+twist*zj) - xi*sin(yi+twist*zi);
            dz = zj - zi;
        } else if (cell_typ > 10 && cell_typ < 20){
            double dx1 = xj - xi;
            double dx2 = yj - yi;
            double dx3 = zj - zi;
            dx = nlist->LatUVec[0] * dx1 + nlist->LatUVec[3] * dx2 + nlist->LatUVec[6] * dx3;
            dy = nlist->LatUVec[1] * dx1 + nlist->LatUVec[4] * dx2 + nlist->LatUVec[7] * dx3;
            dz = nlist->LatUVec[2] * dx1 + nlist->LatUVec[5] * dx2 + nlist->LatUVec[8] * dx3;
        } else if (cell_typ == 0) {
            dx = xj - xi;
            dy = yj - yi;
            dz = zj - zi;
        }
        r0_sqr = dx*dx + dy*dy + dz*dz;
        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
            mcsh_solid_function(dx, dy, dz, r0_sqr, A, B, alpha, beta, miu, deriv);
            sum_miu1 += miu[0];
            sum_miu2 += miu[1];
            sum_miu3 += miu[2];
            sum_dmiu1_dxj[j+1] += deriv[0];
            sum_dmiu1_dyj[j+1] += deriv[1];
            sum_dmiu1_dzj[j+1] += deriv[2];
            sum_dmiu2_dxj[j+1] += deriv[3];
            sum_dmiu2_dyj[j+1] += deriv[4];
            sum_dmiu2_dzj[j+1] += deriv[5];
            sum_dmiu3_dxj[j+1] += deriv[6];
            sum_dmiu3_dyj[j+1] += deriv[7];
            sum_dmiu3_dzj[j+1] += deriv[8];
        }
        sum_dmiu1_dFj[0] += sum_dmiu1_dxj[j+1] * dx;
        sum_dmiu1_dFj[1] += sum_dmiu1_dxj[j+1] * dy;
        sum_dmiu1_dFj[2] += sum_dmiu1_dxj[j+1] * dz;
        sum_dmiu1_dFj[3] += sum_dmiu1_dyj[j+1] * dy;
        sum_dmiu1_dFj[4] += sum_dmiu1_dyj[j+1] * dz;
        sum_dmiu1_dFj[5] += sum_dmiu1_dzj[j+1] * dz;
        sum_dmiu2_dFj[0] += sum_dmiu2_dxj[j+1] * dx;
        sum_dmiu2_dFj[1] += sum_dmiu2_dxj[j+1] * dy;
        sum_dmiu2_dFj[2] += sum_dmiu2_dxj[j+1] * dz;
        sum_dmiu2_dFj[3] += sum_dmiu2_dyj[j+1] * dy;
        sum_dmiu2_dFj[4] += sum_dmiu2_dyj[j+1] * dz;
        sum_dmiu2_dFj[5] += sum_dmiu2_dzj[j+1] * dz;
        sum_dmiu3_dFj[0] += sum_dmiu3_dxj[j+1] * dx;
        sum_dmiu3_dFj[1] += sum_dmiu3_dxj[j+1] * dy;
        sum_dmiu3_dFj[2] += sum_dmiu3_dxj[j+1] * dz;
        sum_dmiu3_dFj[3] += sum_dmiu3_dyj[j+1] * dy;
        sum_dmiu3_dFj[4] += sum_dmiu3_dyj[j+1] * dz;
        sum_dmiu3_dFj[5] += sum_dmiu3_dzj[j+1] * dz;
    }
    *sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3);
    for (int j = 0; j < 6; j++) {
        sum_dmiu1_dFj[j] *= sum_miu1 * group_coefficient * 2.0;
        sum_dmiu2_dFj[j] *= sum_miu2 * group_coefficient * 2.0;
        sum_dmiu3_dFj[j] *= sum_miu3 * group_coefficient * 2.0;
        gmp_str->dX3_dF[ii][j][m] += sum_dmiu1_dFj[j];
        gmp_str->dX3_dF[ii][j][m] += sum_dmiu2_dFj[j];
        gmp_str->dX3_dF[ii][j][m] += sum_dmiu3_dFj[j];
    }
    double dmdx, dmdy, dmdz;
    for (int j = -1; j < nlist->Nneighbors[ii]; ++j) {
        dmdx = (sum_miu1 * sum_dmiu1_dxj[j+1] + sum_miu2 * sum_dmiu2_dxj[j+1] + sum_miu3 * sum_dmiu3_dxj[j+1]) * group_coefficient * 2.0;
        dmdy = (sum_miu1 * sum_dmiu1_dyj[j+1] + sum_miu2 * sum_dmiu2_dyj[j+1] + sum_miu3 * sum_dmiu3_dyj[j+1]) * group_coefficient * 2.0;
        dmdz = (sum_miu1 * sum_dmiu1_dzj[j+1] + sum_miu2 * sum_dmiu2_dzj[j+1] + sum_miu3 * sum_dmiu3_dzj[j+1]) * group_coefficient * 2.0;
        if (j < 0) {
            gmp_str->dX3_dX[ii][0][m] += dmdx;
            gmp_str->dX3_dY[ii][0][m] += dmdy;
            gmp_str->dX3_dZ[ii][0][m] += dmdz;
        }
        else {
            gmp_str->dX3_dX[ii][_indices[j]][m] += dmdx;
            gmp_str->dX3_dY[ii][_indices[j]][m] += dmdy;
            gmp_str->dX3_dZ[ii][_indices[j]][m] += dmdz;
        }

        gmp_str->dX3_dX[ii][0][m] -= dmdx;
        gmp_str->dX3_dY[ii][0][m] -= dmdy;
        gmp_str->dX3_dZ[ii][0][m] -= dmdz;
    }
    free(sum_dmiu1_dxj);
    free(sum_dmiu1_dyj);
    free(sum_dmiu1_dzj);
    free(sum_dmiu2_dxj);
    free(sum_dmiu2_dyj);
    free(sum_dmiu2_dzj);
    free(sum_dmiu3_dxj);
    free(sum_dmiu3_dyj);
    free(sum_dmiu3_dzj);
    free(sum_dmiu1_dFj);
    free(sum_dmiu2_dFj);
    free(sum_dmiu3_dFj);
    free(_indices);
}

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
) {
    SolidGMPFunction mcsh_solid_function = get_solid_mcsh_function(mcsh_order, group_index);
    double sum_miu1 = 0.0, sum_miu2 = 0.0, sum_miu3 = 0.0, sum_miu4 = 0.0, sum_miu5 = 0.0, sum_miu6 = 0.0;
    double* sum_dmiu1_dxj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu2_dxj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu3_dxj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu4_dxj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu5_dxj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu6_dxj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu1_dyj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu2_dyj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu3_dyj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu4_dyj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu5_dyj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu6_dyj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu1_dzj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu2_dzj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu3_dzj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu4_dzj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu5_dzj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu6_dzj = (double *) calloc(nlist->Nneighbors[ii]+1, sizeof(double));
    double* sum_dmiu1_dFj = (double *) calloc(6, sizeof(double));
    double* sum_dmiu2_dFj = (double *) calloc(6, sizeof(double));
    double* sum_dmiu3_dFj = (double *) calloc(6, sizeof(double));
    double* sum_dmiu4_dFj = (double *) calloc(6, sizeof(double));
    double* sum_dmiu5_dFj = (double *) calloc(6, sizeof(double));
    double* sum_dmiu6_dFj = (double *) calloc(6, sizeof(double));                    
    double miu[6], deriv[18], xj, yj, zj, dx = 0.0, dy = 0.0, dz = 0.0, r0_sqr;
    int _count = 0;
    int *_indices = (int*) malloc(sizeof(int)*nlist->Nneighbors[ii]);
    int _tmp = -1, idx_neigh, elem_typ;
    // -1 accounts for the fact that atom i in the fundamental unit cell is not included. However, it's mirror images are included.
    for (int j = -1; j < nlist->Nneighbors[ii]; ++j) {
        int neigh_atom_element_order;
        if (j < 0) {
            xj = xi;
            yj = yi;
            zj = zi;
            // gives the atomic number of the atom for a given index of atom - correponds to index in .ion file. 
            int neigh_atom_element_index = atom_indices[i];
            // Gives the order a given element appears in the .ion file and the primitive gaussians arrays
            neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
            _tmp = i;
        } else {
            idx_neigh = (nlist->neighborList + ii)->array[j];
            elem_typ = (nlist->neighborAtmTyp + ii)->array[j];
            xj = atompos[3*idx_neigh] + L1 * (nlist->neighborList_imgX + ii)->array[j];
            yj = atompos[3*idx_neigh+1] + L2 * (nlist->neighborList_imgY + ii)->array[j];
            zj = atompos[3*idx_neigh+2] + L3 * (nlist->neighborList_imgZ + ii)->array[j];
            int neigh_atom_element_index = atom_type_to_indices[elem_typ];
            neigh_atom_element_order = element_index_to_order[neigh_atom_element_index];
            if (idx_neigh == i) {
                _indices[j] = 0;
            } else if (_tmp < 0) {
                perror("Error in GMP descriptor calculation: indices loop broken\n");
            } else if (_tmp != idx_neigh) {
                _count += 1;
                _tmp = idx_neigh;
                _indices[j] = _count;
            } else {
                _indices[j] = _count;
            }
        }
        if (cell_typ > 20 && cell_typ < 30) { // add for nonorthogonal
            dx = xj*cos(yj+twist*zj) - xi*cos(yi+twist*zi);
            dy = xj*sin(yj+twist*zj) - xi*sin(yi+twist*zi);
            dz = zj - zi;
        } else if (cell_typ > 10 && cell_typ < 20){
            double dx1 = xj - xi;
            double dx2 = yj - yi;
            double dx3 = zj - zi;
            dx = nlist->LatUVec[0] * dx1 + nlist->LatUVec[3] * dx2 + nlist->LatUVec[6] * dx3;
            dy = nlist->LatUVec[1] * dx1 + nlist->LatUVec[4] * dx2 + nlist->LatUVec[7] * dx3;
            dz = nlist->LatUVec[2] * dx1 + nlist->LatUVec[5] * dx2 + nlist->LatUVec[8] * dx3;
        } else if (cell_typ == 0) {
            dx = xj - xi;
            dy = yj - yi;
            dz = zj - zi;
        }
        r0_sqr = dx*dx + dy*dy + dz*dz;
        for (int g = 0; g < ngaussians[neigh_atom_element_order]; ++g){
            double B = atom_gaussian[neigh_atom_element_order][g*2], beta = atom_gaussian[neigh_atom_element_order][g*2+1];
            mcsh_solid_function(dx, dy, dz, r0_sqr, A, B, alpha, beta, miu, deriv);
            sum_miu1 += miu[0];
            sum_miu2 += miu[1];
            sum_miu3 += miu[2];
            sum_miu4 += miu[3];
            sum_miu5 += miu[4];
            sum_miu6 += miu[5];
            sum_dmiu1_dxj[j+1] += deriv[0];
            sum_dmiu1_dyj[j+1] += deriv[1];
            sum_dmiu1_dzj[j+1] += deriv[2];
            sum_dmiu2_dxj[j+1] += deriv[3];
            sum_dmiu2_dyj[j+1] += deriv[4];
            sum_dmiu2_dzj[j+1] += deriv[5];
            sum_dmiu3_dxj[j+1] += deriv[6];
            sum_dmiu3_dyj[j+1] += deriv[7];
            sum_dmiu3_dzj[j+1] += deriv[8];
            sum_dmiu4_dxj[j+1] += deriv[9];
            sum_dmiu4_dyj[j+1] += deriv[10];
            sum_dmiu4_dzj[j+1] += deriv[11];
            sum_dmiu5_dxj[j+1] += deriv[12];
            sum_dmiu5_dyj[j+1] += deriv[13];
            sum_dmiu5_dzj[j+1] += deriv[14];
            sum_dmiu6_dxj[j+1] += deriv[15];
            sum_dmiu6_dyj[j+1] += deriv[16];
            sum_dmiu6_dzj[j+1] += deriv[17];
        }
        sum_dmiu1_dFj[0] += sum_dmiu1_dxj[j+1] * dx;
        sum_dmiu1_dFj[1] += sum_dmiu1_dxj[j+1] * dy;
        sum_dmiu1_dFj[2] += sum_dmiu1_dxj[j+1] * dz;
        sum_dmiu1_dFj[3] += sum_dmiu1_dyj[j+1] * dy;
        sum_dmiu1_dFj[4] += sum_dmiu1_dyj[j+1] * dz;
        sum_dmiu1_dFj[5] += sum_dmiu1_dzj[j+1] * dz;
        sum_dmiu2_dFj[0] += sum_dmiu2_dxj[j+1] * dx;
        sum_dmiu2_dFj[1] += sum_dmiu2_dxj[j+1] * dy;
        sum_dmiu2_dFj[2] += sum_dmiu2_dxj[j+1] * dz;
        sum_dmiu2_dFj[3] += sum_dmiu2_dyj[j+1] * dy;
        sum_dmiu2_dFj[4] += sum_dmiu2_dyj[j+1] * dz;
        sum_dmiu2_dFj[5] += sum_dmiu2_dzj[j+1] * dz;
        sum_dmiu3_dFj[0] += sum_dmiu3_dxj[j+1] * dx;
        sum_dmiu3_dFj[1] += sum_dmiu3_dxj[j+1] * dy;
        sum_dmiu3_dFj[2] += sum_dmiu3_dxj[j+1] * dz;
        sum_dmiu3_dFj[3] += sum_dmiu3_dyj[j+1] * dy;
        sum_dmiu3_dFj[4] += sum_dmiu3_dyj[j+1] * dz;
        sum_dmiu3_dFj[5] += sum_dmiu3_dzj[j+1] * dz;
        sum_dmiu4_dFj[0] += sum_dmiu4_dxj[j+1] * dx;
        sum_dmiu4_dFj[1] += sum_dmiu4_dxj[j+1] * dy;
        sum_dmiu4_dFj[2] += sum_dmiu4_dxj[j+1] * dz;
        sum_dmiu4_dFj[3] += sum_dmiu4_dyj[j+1] * dy;
        sum_dmiu4_dFj[4] += sum_dmiu4_dyj[j+1] * dz;
        sum_dmiu4_dFj[5] += sum_dmiu4_dzj[j+1] * dz;
        sum_dmiu5_dFj[0] += sum_dmiu5_dxj[j+1] * dx;
        sum_dmiu5_dFj[1] += sum_dmiu5_dxj[j+1] * dy;
        sum_dmiu5_dFj[2] += sum_dmiu5_dxj[j+1] * dz;
        sum_dmiu5_dFj[3] += sum_dmiu5_dyj[j+1] * dy;
        sum_dmiu5_dFj[4] += sum_dmiu5_dyj[j+1] * dz;
        sum_dmiu5_dFj[5] += sum_dmiu5_dzj[j+1] * dz;
        sum_dmiu6_dFj[0] += sum_dmiu6_dxj[j+1] * dx;
        sum_dmiu6_dFj[1] += sum_dmiu6_dxj[j+1] * dy;
        sum_dmiu6_dFj[2] += sum_dmiu6_dxj[j+1] * dz;
        sum_dmiu6_dFj[3] += sum_dmiu6_dyj[j+1] * dy;
        sum_dmiu6_dFj[4] += sum_dmiu6_dyj[j+1] * dz;
        sum_dmiu6_dFj[5] += sum_dmiu6_dzj[j+1] * dz;
    }
    *sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                    sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6);
    for (int j = 0; j < 6; j++) {
        sum_dmiu1_dFj[j] *= sum_miu1 * group_coefficient * 2.0;
        sum_dmiu2_dFj[j] *= sum_miu2 * group_coefficient * 2.0;
        sum_dmiu3_dFj[j] *= sum_miu3 * group_coefficient * 2.0;
        sum_dmiu4_dFj[j] *= sum_miu4 * group_coefficient * 2.0;
        sum_dmiu5_dFj[j] *= sum_miu5 * group_coefficient * 2.0;
        sum_dmiu6_dFj[j] *= sum_miu6 * group_coefficient * 2.0;
        gmp_str->dX3_dF[ii][j][m] += sum_dmiu1_dFj[j];
        gmp_str->dX3_dF[ii][j][m] += sum_dmiu2_dFj[j];
        gmp_str->dX3_dF[ii][j][m] += sum_dmiu3_dFj[j];
        gmp_str->dX3_dF[ii][j][m] += sum_dmiu4_dFj[j];
        gmp_str->dX3_dF[ii][j][m] += sum_dmiu5_dFj[j];
        gmp_str->dX3_dF[ii][j][m] += sum_dmiu6_dFj[j];
    }
    double dmdx, dmdy, dmdz;
    for (int j = -1; j < nlist->Nneighbors[ii]; ++j) {
        dmdx = (sum_miu1 * sum_dmiu1_dxj[j+1] + sum_miu2 * sum_dmiu2_dxj[j+1] +
                sum_miu3 * sum_dmiu3_dxj[j+1] + sum_miu4 * sum_dmiu4_dxj[j+1] +
                sum_miu5 * sum_dmiu5_dxj[j+1] + sum_miu6 * sum_dmiu6_dxj[j+1]) * group_coefficient * 2.0;

        dmdy = (sum_miu1 * sum_dmiu1_dyj[j+1] + sum_miu2 * sum_dmiu2_dyj[j+1] +
                sum_miu3 * sum_dmiu3_dyj[j+1] + sum_miu4 * sum_dmiu4_dyj[j+1] +
                sum_miu5 * sum_dmiu5_dyj[j+1] + sum_miu6 * sum_dmiu6_dyj[j+1]) * group_coefficient * 2.0;

        dmdz = (sum_miu1 * sum_dmiu1_dzj[j+1] + sum_miu2 * sum_dmiu2_dzj[j+1] +
                sum_miu3 * sum_dmiu3_dzj[j+1] + sum_miu4 * sum_dmiu4_dzj[j+1] +
                sum_miu5 * sum_dmiu5_dzj[j+1] + sum_miu6 * sum_dmiu6_dzj[j+1]) * group_coefficient * 2.0;
        if (j < 0) {
            gmp_str->dX3_dX[ii][0][m] += dmdx;
            gmp_str->dX3_dY[ii][0][m] += dmdy;
            gmp_str->dX3_dZ[ii][0][m] += dmdz;
        }
        else {
            gmp_str->dX3_dX[ii][_indices[j]][m] += dmdx;
            gmp_str->dX3_dY[ii][_indices[j]][m] += dmdy;
            gmp_str->dX3_dZ[ii][_indices[j]][m] += dmdz;
        }

        gmp_str->dX3_dX[ii][0][m] -= dmdx;
        gmp_str->dX3_dY[ii][0][m] -= dmdy;
        gmp_str->dX3_dZ[ii][0][m] -= dmdz;
    }
    free(sum_dmiu1_dxj);
    free(sum_dmiu1_dyj);
    free(sum_dmiu1_dzj);
    free(sum_dmiu2_dxj);
    free(sum_dmiu2_dyj);
    free(sum_dmiu2_dzj);
    free(sum_dmiu3_dxj);
    free(sum_dmiu3_dyj);
    free(sum_dmiu3_dzj);
    free(sum_dmiu4_dxj);
    free(sum_dmiu4_dyj);
    free(sum_dmiu4_dzj);
    free(sum_dmiu5_dxj);
    free(sum_dmiu5_dyj);
    free(sum_dmiu5_dzj);
    free(sum_dmiu6_dxj);
    free(sum_dmiu6_dyj);
    free(sum_dmiu6_dzj);
    free(sum_dmiu1_dFj);
    free(sum_dmiu2_dFj);
    free(sum_dmiu3_dFj);
    free(sum_dmiu4_dFj);
    free(sum_dmiu5_dFj);
    free(sum_dmiu6_dFj);
    free(_indices);
}