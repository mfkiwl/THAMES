/**
@file Solution.h
@brief Declare the Solution class.

@note Perhaps this is really a derived class from the ChemicalSystem.
@todo Investigate whether to make this a derived class.

*/
#ifndef SOLUTION_H
#define SOLUTION_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>

#include "../GEMS3K-standalone/Resources/include/GEMS3K/node.h"
// #include <GEMS3K/io_arrays.h>

#include "global.h"

using namespace std;

/**
@class Solution
@brief Declares the Solution class.

This class keeps track of the composition of the aqueous solution, and
has a lot of member variables that enable direct communication with the
GEM3K library.

The solution is kept separate from the ChemicalSystem because we
often want to just know the driving force (saturation index) for
dissolution or precipitation of solid phases.  Therefore we may
want to set the upper limit of solid phases to a very low number
approximating zero.

@note This class's members consist almost entirely of a subset of the
members of the ChemicalSystem class.

@todo Find out why we should not make the ChemicalSystem object a member
of this class and then access the members that way.
*/
class Solution {

private:

TNode *node_;                       /**< Pointer to the GEM3K Tnode object */

bool jsonFormat_;                   /**< True if GEM input data files are in JSON format */

double T_;                          /**< System-wide temperature [K] */
double P_;                          /**< System-wide pressure [Pa] */
double Vs_;                         /**< System total volume [m<sup>3</sup>] */
double Ms_;                         /**< System total mass [kg] */
double Gs_;                         /**< System total Gibbs energy [J] */
double Hs_;                         /**< System total enthalpy [J] */
double ionicStrength_;              /**< Solution ionic strength [mol/kgw] */
double pH_;                         /**< Solution pH */
double pe_;                         /**< Solution pe */
double Eh_;                         /**< Solution Eh [volts] */

long int nodeHandle_;               /**< integer flag used to identify a node */
long int nodeStatus_;               /**< integer flag used to identify node's status */
long int iterDone_;                 /**< number of iterations performed in the most
                                                recent GEM calculation on the node */

int timesGEMFailed_;                /**< number of consecutive times a call to GEM_run
                                                has failed */

int maxGEMFails_;                   /**< maximum allowed number of consecutive GEM_run fails */

unsigned int numICs_;                /**< Number of independent components (IC) */
unsigned int numDCs_;                /**< Number of dependent components (DC) */
unsigned int numGEMPhases_;             /**< Number of GEM phases in the CSD */
unsigned int numSolutionPhases_;     /**< Number of GEM solution phases in the CSD;
                                              solution phases are non-stoichiometric */

double *ICMoles_;                   /**< List of number of moles of each IC in system */
double *ICResiduals_;               /**< List of errors in IC moles for mass balance */
double *ICChemicalPotential_;                 /**< List of chemical potentials of each IC, in
                                            the GEM dual solution */
double *DCMoles_;                   /**< List of moles of each DC */
double *DCActivityCoeff_;           /**< List of activity coefficients for each DC */

double *GEMPhaseMoles_;                /**< List of moles of each phase in the system */
double *GEMPhaseMass_;                 /**< List of mass of each phase in the system [kg] */
double *GEMPhaseVolume_;               /**< List of volume of each phase in the system
                                            [m<sup>3</sup> */

double *pGEMPhaseStoich_;               /**< List of amount of moles of each IC in a
                                                GEM CSD phase (pointer form) */
/**
@brief Solid stoichiometry list for communicating with GEM-IPM.

@warning Not sure how this variable is used
*/
double *solidStoich_;  

double *carrier_;                   /**< List of moles of carrier (solvent) in
                                            multicomponent asymmetric phases */
double *surfaceArea_;               /**< List of specific surface area of each
                                            phase, in m<sup>2</sup>/kg */
double *DCUpperLimit_;              /**< List of upper bound on moles of each DC */
double *DCLowerLimit_;              /**< List of lower bound on moles of each DC,
                                            generally non-zero for numerical
                                            stability of the thermodynamic
                                            calculations */
vector<string> ICName_;             /**< Names of ICs in the GEM CSD */
vector<string> DCName_;             /**< Names of DCs in the GEM CSD */
vector<string> GEMPhaseName_;          /**< Names of phases in the GEM CSD */

/**
@brief Saturation index of each phase in the GEM CSD.

The departure from equilibrium between a given solid phase and an aqueous solution
is characterized by the saturation index, `SI_`, which is defined as the activity
product for the dissolution reaction divided by the equilibrium value of that activity
product.  This variable stores the current SI for each solid phase in the GEM CSD.
*/
vector<double> SI_;

double crystalStrain_;                  /**< Assigned strain from FE model */

bool verbose_;                      /**< Flag to determine verbose output */
bool debug_;                        /**< Flag to determine debugging output */

public:

/**
@brief Constructor.

The class members can all be assigned by reading GEM3K input files,
which are passed to the constructor.

@param dchFileName is the name of the GEM DCH file
@param dbrFileName is the name of the GEM DBR (data bridge) file
@param verbose is true if verbose output should be produced
@param debug is true if debugging output should be produced
*/

Solution (const string &dchFileName,
          const string &dbrFileName,
          const bool verbose,
          const bool debug = false);

/**
@brief Destructor.

A destructor is needed to free the memory allocated by the constructor.
*/
~Solution ();

/**
@brief Determine whether GEM input data files are in JSON format

@param masterFileName is the lst name containing the others
@return true if JSON is indicated, false otherwise
*/
bool isInputFormatJSON(const char *masterFileName);

/**
@brief Get the name of the three JSON input files for GEMS data

@param masterFileName is the lst name containing the others
@param dchName is the name of the DCH data file
@param ipmName is the name of the IPM data file
@param dbrName is the name of the DBR data file
*/
void getJSONFiles(const char *masterFileName,
                  string &dchName,
                  string &ipmName,
                  string &dbrName);

/**
@brief Calculate the thermodynanmic equilibrium state of the solution.

This method calculates the thermodynamic equilibrium state of the solution
under the assumption that no solids are allowed to precipitate or dissolve.
As a result, GEM3K calculates and stores the saturation indices with respect
to each dependent component in the GEM chemical system definition (CSD).

@todo Find out how GEM knows not to precipitate or dissolve in this version.

@param isFirst is true if this is the first equilibrium calculation, false otherwise
*/
void calculateState (bool isFirst);

/**
@brief Get the list of multicomponent phase volumes.

@return a pointer to the array of phase volumes [m<sup>3</sup>]
*/
double *getGEMPhaseVolume () const
{
  return GEMPhaseVolume_; 
}

/**
@brief Get the volume of one multicomponent phase.

@param idx is the index of the phase in the array of phase volumes
@return the volume [m<sup>3</sup>]
*/
double getGEMPhaseVolume (const unsigned int idx)
{
  if (idx >= numGEMPhases_) {
      cout << "idx beyond the range of numGEMPhases_" << endl;
      cerr << "idx beyond the range of numGEMPhases_" << endl;
      exit(1);
  }
  return GEMPhaseVolume_[idx];
}

/**
@brief Get the total number of phases in the thermodynamic system.

@return the total number of defined phases
*/
int getNumGEMPhases (void)
{
    return numGEMPhases_;
}

/**
@brief Get the name of a phase specified by its index in the phase array.

@param idx is the index of the phase in the array of phase names
@return the phase name
*/
string &getGEMPhaseName (const unsigned int idx)
{
    if (idx >= GEMPhaseName_.size()) {
        cout << "idx beyond the range of GEMPhaseName_" << endl;
        cerr << "idx beyond the range of GEMPhaseName_" << endl;
        exit(1);
    }
    return (string &)GEMPhaseName_[idx];
}

/**
@brief Get the list of multicomponent phase moles.

@return a pointer to the array of phase moles
*/
double *getGEMPhaseMoles (void) const
{
    return GEMPhaseMoles_;
}

/**
@brief Get the list of all multicomponent phase names.

@return a pointer to the array of phase names
*/
vector<string> getGEMPhaseName (void) const
{
    return GEMPhaseName_;
}

/**
@brief Get the number of iterations needed for equilibration by GEM.

@return the number of iterations executed by GEM
*/
long int getIterdone (void);

/**
@brief Set the number of moles of a given independent component (IC) in the solution.

@param idx is the index of the IC in the array of ICs
@param val is the number of moles to set for this component
*/
void setICMoles (const unsigned int idx,
                 const double val)
{
    if (idx >= numICs_) {
        cout << "index beyond the range of ICnum, "
             << "so exit the program." << endl;
        cerr << "index beyond the range of ICnum, "
             << "so exit the program." << endl;
        exit(1);
    } else {
        ICMoles_[idx] = val;
    }
    return;
}

/**
@brief Set the number of moles of all independent components (ICs) in the solution.

@param val is an array of mole values, one for each IC
*/
void setICMoles (vector<double> val)
{
    for (int i = 0; i < numICs_; i++) {
      ICMoles_[i] = val[i];
    }

    return;
}

/**
@brief Get the list of all independent component (IC) moles in the solution.

@return a pointer to the array of IC moles
*/
double *getICMoles (void)
{
    return ICMoles_;
}

/**
@brief Get the number of moles of an independent component (IC) specified by its index.

@todo Perform bounds checking and throw an error if out of bounds.

@param idx is the index of the IC to query
@return the number of moles of that IC in the system
*/
double getICMoles (const unsigned int idx)
{
    return ICMoles_[idx];
}

/**
@brief Set the saturation index of each dependent component.
*/
void setSI (void)
{
    SI_.clear();
    double *Falp;
    Falp = (node_->ppmm())->Falp;
    
    for (int i = 0; i < numGEMPhases_; i++) {
        double si = pow(10,Falp[i]);
        SI_.push_back(si);
    }
    return;
}

/**
@brief Get the list of all saturation indices.

@return the vector of saturation indices, one for each solid phase
*/
vector<double> getSI (void)
{
    return SI_;
}

/**
@brief Get the saturation index of a phase based on its index in the array.

@todo Perform bounds checking and throw an error if out of bounds.

@param GEMPhaseId is the position of the queried phase in the array of SIs
@return the saturation index of that phase
*/
double getSI (int GEMPhaseId)
{
    return SI_[GEMPhaseId];
}
   
/**
@brief Get a pointer to the GEM node doing the equilibrium calculations.

@return the GEM node doing the calculations
*/
TNode *getNode (void)
{
    return node_;
}  
 
/**
@brief Calculate the site strain based on local crystallization pressure.

The calculation of local crystallization pressure happens in this method,
based on the saturation index of the phase.  Therefore, the method is
currently restricted to problems of sulfate attack, and it probably can
be generalized somewhat.

Assuming that the saturation index \f$\beta\f$ is known for the phase,
then the crystallization pressure should be the difference in the Laplace pressure
between the large pore entrance, \f$r_{cp}\f$, and the size of the average gel
pore, \f$r_{gp}\f$.  This pressure difference is

\f[
    p_c = 2 \gamma \left[ \frac{1}{r_{gp} - \delta} - \frac{1}{r_{cp} - \delta} \right]
\f]

where we subtract \f$\delta\f$, the liquid film thickness, so that the terms in
the denominator are the radii of curvature of the actual crystal in these two
locations.

In THAMES, we <i>assume</i> that the largest pore entrance to the gel porosity
is about half the size of a lattice site.  The usual dimension of a lattice site
in THAMES is 1 micrometer--- <b>although this is not necessary</b>--- so we assume
that \f$r_{cp} = 500\f$ nm.

The Thompson-Freundlich relation tells us the size of a crystal that is in equilibrium
with a supersaturated solution with saturation index, relative to that crystal,
of \f$\beta\f$.  The condition of equilibrium is that the chemical potential of
formula units in the crystal be equal to the chemical potential of dissolved formula
units in the solution:

\f[
    V_c \gamma \kappa_{cr} = R T \ln \frac{Q}{K} = R T \ln \beta
\f]

where \f$V_c\f$ is the stress-free molar volume of the crystal, \f$\gamma\f$ is the
surface energy of the crystal-liquid interface, \f$\kappa_{cr}\f$ is the mean
curvature of the crystal in equilibrium with the solution, <i>R</i> and <i>T</i>
are the gas constant and absolute temperature, respectively, and <i>Q</i> and
<i>K</i> are the activity product and equilibrium constant of the crystal dissociation
reaction.  Therefore, the mean curvature is

\f[
    \kappa_{cr} = \frac{R T}{V_c \gamma} \ln \beta
\f]

And since \f$\kappa_{cr} \equiv 2/r_{cr}\f$, we have \f$r_{cr} = 2/\kappa_{cr}\f$
But \f$r_{gp}\f$ is the radius of curvature of the gel pore, not the crystal, so
we must add the liquid film thickness: \f$r_{gp} = r_{cr} + \delta\f$.

With the crystallization pressure, <i>p</i><sub>c</sub>, calculated, we must
now estimate the local strain.  Ordinarily, this could be done only with detailed
knowledge of the local pore structure.  However, to make the calculation tractable,
we adopt the assumptions of poromechanics, in which case the strain is estimated
by a volume-averaged equation,

\f[
    \epsilon_x \approx \left( \frac{1}{3 K_p} - \frac{1}{3 K_s} \right)
    \left[ \phi p_c + pl - p_{atm} \right]
\f]

where \f$\epsilon_x\f$ is the linear strain (strain is assumed to be isotropic),
<i>K</i><sub>p</sub> and <i>K</i><sub>s</sub> are the volume averaged bulk and
shear moduli of the porous body, \f$\phi\f$ is the local porosity, <i>p</i><sub>l</sub>
is the hydrostatic pressure in the liquid film, and <i>p</i><sub>atm</sub> is
atmospheric pressure.

@todo Consider renaming this function to calcCrystallizationStrain.

@param SI is the saturation index of the growing crystal
@param poreVolFrac is the local porosity
@param Kp is the effective bulk modulus of the porous body
@param Ks is the effective shear modulus of the porous body
@return the calculated crystallization strain
*/
double calculateCrystalStrain (double SI,
                               double poreVolFrac,
                               double Kp,
                               double Ks);

/**
@brief Make sure that the IC moles are all greater than e-10
to ensure stability of the IPM calculation.

The exception is for charge, which should always be set to zero
for these kinds of simulations.
*/
void checkICMoles (void)
{
    for (unsigned int i = 0; i < numICs_; i++) {
        if (ICMoles_[i] < 1.0e-9) ICMoles_[i] = 1.0e-9;
        if (i == numICs_ - 1) ICMoles_[i] = 0.0;
    }
    return;
}

/**
@brief Set the jsonFormat_ flag

@param jsonFormat is true if input files are in JSON format
*/
void setJSONformat (const bool jsonFormat)
{
    jsonFormat_ = jsonFormat;
    return;
}

/**
@brief Get the jsonformat_ flag

@return the jsonformat flag
*/
bool getJSONFormat (void) const
{
    return jsonFormat_;
}

/**
@brief Set the number of times in a row that GEM_run has failed.

Calculations of equilibrium state sometimes fail to converge.  It is
possible that some small tweak in the composition could help the
algorithm converge, so we allow up to `maxGEMFails_' attempts for
convergence, each time tweaking the composition in some way
before giving up and throwing an exception.

@param nTimes is the number of consecutive times GEM_run has failed
*/
void setTimesGEMFailed (const int nTimes)
{
    timesGEMFailed_ = (nTimes >= 0) ? nTimes : 0;
    return;
}

/**
@brief Get the number of consecutive times a call to GEM_run has failed.

@return the number of consecutive times a call to GEM_run has failed.
*/
int getTimesGEMFailed (void) const
{
    return timesGEMFailed_;
}


/**
@brief Set the verbose flag

@param isverbose is true if verbose output should be produced
*/
void setVerbose (const bool isverbose)
{
    verbose_ = isverbose;
    return;
}

/**
@brief Get the verbose flag

@return the verbose flag
*/
bool getVerbose (void) const
{
    return verbose_;
}

/**
@brief Set the debug flag

@param isdebug is true if debug output should be produced
*/
void setDebug (const bool isdebug)
{
    debug_ = isdebug;
    return;
}

/**
@brief Get the debug flag

@return the debug flag
*/
bool getDebug (void) const
{
    return debug_;
}

};      // End of the Solution class

#endif
