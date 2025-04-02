#ifndef BESSELFUNCTIONS_H
#define BESSELFUNCTIONS_H

#include "mlff_types.h"
#include "isddft.h"


double spherical_bessel(int order, double z);
void modifiedsphericalbesselfirst_l_01(double x, double *vals);
double modifiedsphericalbesselfirst(int n, double x);
double derivative_modifiedsphericalbesselfirst(int n, double x);
void roots_sph_bessel_notinuse(int order, double *roots);

#endif
