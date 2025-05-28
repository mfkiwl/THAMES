/**
@file StrainEnergy.h
@brief Defines the strainenergy vector.

The strainenergy vector stores the elastic strain energy associated with
each microstructure phase, which is then communicated to the GEM3K library
as a modification to the standard Gibbs energy of formation of those phases.
This is therefore a way to modify the usual thermodynamic calculations to
account for elastic strain energy.
*/
#ifndef SRC_THAMESLIB_STRAINENERGY_H_
#define SRC_THAMESLIB_STRAINENERGY_H_

#include <math.h>
#include <string>
#include <vector>

using namespace std;

extern vector<double> strainenergy;

#endif // SRC_THAMESLIB_STRAINENERGY_H_
