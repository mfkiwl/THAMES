/**
@file KineticModel.h
@brief Declaration of the KineticModel class.

@section Introduction
In THAMES, the `KineticModel` class can be perceived as the engine that calculates
the kinetic changes in the system during a given time increment.  The primary
kinetic aspect that is calculated is the extent of dissolution of mineral phases in
the original clinker.
As of March 2017, the kinetic model implements the Parrot and Killoh (PK) model of
1984 [1]---in the same way as described by Lothenbach and
Winnefeld [2]---for cement clinker phases.

The PK model accounts only for dissolution of the four major clinker minerals:
alite (impure C<sub>3</sub>S), belite (impure C<sub>2</sub>S), aluminate (C<sub>3</sub>A),
and ferrite (C<sub>4</sub>AF).  The model provides mathematical equations for the rates
of three broad rate-controlling phenomena:

    -# nucleation and growth,
    -# early-age diffusion, and
    -# late-age diffusion.
    
For each clinker mineral <i>i</i>, the three rate equations are

@f{eqnarray*}
R_{ng}(i) &=& \frac{A K_1(i)}{N_1(i)} \left( 1 - \alpha \right) \left( - \ln \left( 1 - \alpha \right) \right)^{1 - N_1(i)} \\
R_{de}(i) &=& \frac{K_2(i) \left( 1 - \alpha \right)^{2/3}}{1 - \left( 1 - \alpha \right)^{1/3}} \\
R_{dl}(i) &=& K_3(i) \left( 1 - \alpha \right)^{N_3(i)}
@f}

where <i>A</i> is the overall surface area of the cement powder (cm<sup>2</sup>/kg)
and \f$\alpha\f$ is the overall degree of hydration on a mass basis.
<i>K</i><sub>1</sub>, <i>K</i><sub>2</sub>, <i>K</i><sub>3</sub>, <i>N</i><sub>1</sub>,
and <i>N</i><sub>3</sub> are constants defined for each clinker mineral.  The values
of these constants used by Parrot and Killoh in Ref. [1] are shown in the
table.  In any particular time interval, the predicted rate of dissolution of a
clinker mineral

@f[
R(i) = \min (R_{ng}(i),R_{de}(i),R_{dl}(i)) \cdot f(\text{RH}) \cdot g(w/c)
@f]

where \f$f(\text{RH})\f$ and \f$g(w/c)\f$ account empirically for the influences of
relative humidity and water-cement mass ratio (w/c), respectively, according to

@f{eqnarray*}
f(\text{RH}) &=& \left( \frac{ \text{RH} - 0.55}{0.45} \right)^4 \\
g(w/c) &=&
\begin{cases}
1 & \text{if}\ \alpha \le 1.333\, w/c \\
(1 + 4.444 w/c - 3.333 \alpha)^4 & \text{otherwise}
\end{cases}
@f}

The new degree of hydration at the end of the time interval is calculated according
to the difference equation

@f[
\alpha(t+\Delta t) = \alpha(t) + R(t) \Delta t
@f]

<table>
<caption id="multi_row">Empirical constants used by Parrot and Killoh</caption>
<tr><th>Parameter               <th>Alite           <th>Belite          <th>Aluminate           <th>Ferrite
<tr><td><i>K</i><sub>1</sub>    <td>1.5             <td>0.5             <td>1.0                 <td>0.37
<tr><td><i>N</i><sub>1</sub>    <td>0.7             <td>1.0             <td>0.85                <td>0.7
<tr><td><i>K</i><sub>2</sub>    <td>0.05            <td>0.006           <td>0.05                <td>0.015
<tr><td><i>K</i><sub>3</sub>    <td>1.1             <td>0.2             <td>1.0                 <td>0.4
<tr><td><i>N</i><sub>3</sub>    <td>3.3             <td>5.0             <td>3.2                 <td>3.7
\end{center}
</table>

@section References

    -# Parrot, L.J., Killoh, D.C., Prediction of cement hydration, British Ceramic
        Proceedings 35 (1984) 41-53.
    -# Lothenbach, B., Winnefeld, F., Thermodynamic modelling of the hydration of
        portland cement, Cement and Concrete Research 36 (2006) 209--226.

*/

#ifndef KINETICSH
#define KINETICSH

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <ctime>
#include "ChemicalSystem.h"
#include "Lattice.h"
#include "global.h"

using namespace std;

/**
@struct KineticData
@brief Stores data about each phase possible in the system for ease of parsing the input files.

In THAMES, phases are identified either thermodynamically--- in the 
GEM data repository--- or microstructurally.  A microstructural phase can
be one, or a combination of more than on, thermodynamically defined phase.

The structure is defined to make it easier to parse the input file for the
kinetic model. It is not used elsewhere in the code.  In fact,
the same members are identified as class variables in the `KineticModel` class.

Most of the members have self-evident meanings:
    - `name` is the name of the microstructure phase
    - `microPhaseId` is the integer id for the phase in the microstructure
    - `GEMPhaseId` is the integer id of the same phase in the GEM Chemical System
        Definition (CSD) 
    - `DCId` is the integer id of the GEM dependent component making up the phase
    - `type` is a string specifying whether the phase is under kinetic control
        or thermodynamic control
    - `scaledMass` is the number of grams of the phase per 100 grams of solid
    - `k1` is the Parrot and Killoh <i>K</i><sub>1</sub> parameter for the phase
    - `k2` is the Parrot and Killoh <i>K</i><sub>2</sub> parameter for the phase
    - `k3` is the Parrot and Killoh <i>K</i><sub>3</sub> parameter for the phase
    - `n1` is the Parrot and Killoh <i>N</i><sub>1</sub> parameter for the phase
    - `n3` is the Parrot and Killoh <i>N</i><sub>3</sub> parameter for the phase
    - `Ea` is the activation energy [J/mol/K]
    - `critDOH` is the critical degree of hydration used in the equation for
        calculating the influence of w/c ratio.
    - `RdId` is a vector of GEM CSD independent component (IC) ids
    - `RdVal` is a vector of the Rd values for each IC
*/

#ifndef KINETICDATASTRUCT
#define KINETICDATASTRUCT
struct KineticData {
    string name;            /**< Name of the microstructure phase */
    int microPhaseId;       /**< Integer id of the microstructure phase */
    int GEMPhaseId;         /**< Integer id of the phase in the GEM CSD */
    int DCId;               /**< Integer id of the DC making up the phase */
    string type;            /**< Specifies kinetic or thermodynamic control */
    double scaledMass;        /**< Mass percent on a total solids basis */
    double k1;                /**< Parrot and Killoh <i>K</i><sub>1</sub> parameter */
    double k2;                /**< Parrot and Killoh <i>K</i><sub>2</sub> parameter */
    double k3;                /**< Parrot and Killoh <i>K</i><sub>3</sub> parameter */
    double n1;                /**< Parrot and Killoh <i>N</i><sub>1</sub> parameter */
    double n3;                /**< Parrot and Killoh <i>N</i><sub>3</sub> parameter */
    double Ea;                /**< Apparent activation energy [J/mol/K] */
    double critDOH;           /**< Critical degree of hydration for w/c effect */
    vector<int> RdId;         /**< Vector of IC ids of the partitioned components in the phase */
    vector<double> RdVal;     /**< Vector of Rd values for each IC */
};
#endif

/**
@class KineticModel
@brief Handles the kinetically controlled phase dissolution and growth.

For now, the Parrot and Killoh model [1] is used to empirically estimate the
mass fraction of each clinker phase that dissolves in a unit time.  Eventually
this can be expanded to handle other kinetically controlled phases outside the
Parrot and Killoh model, such as the growth of portlandite or C--S--H.
*/

class KineticModel {

protected:
    
int numPhases_;             /**< Total number of phases in the kinetic model */
ChemicalSystem *chemSys_;   /**< Pointer to the ChemicalSystem object for this simulation */
Solution *solut_;           /**< Pointer to the aqueous phase object for this simulation */
Lattice *lattice_;          /**< Pointer to the lattice object holding the microstructure */
double wcRatio_;            /**< water-cement mass ratio */
double initSolidMass_;      /**< initial total mass of solids [g] */
double blaine_;             /**< Blaine fineness [m<sup>2</sup>/kg] */
double refBlaine_;          /**< Reference Blaine value for PK model,
                                    usually 385.0 m<sup>2</sup>/kg */
double blaineFactor_;       /**< `blaine_`/`refBlaine_` */
double temperature_;        /**< Temperature [K] */
double refT_;               /**< Reference temperature for PK model [K] */
double sulfateAttackTime_;  /**< Time at which sulfate attack simulation starts [days] */
double leachTime_;          /**< Time at which leaching simulation starts [days] */

/**
@brief List of target sodium concentrations [mol/kgw] at different time steps.

This member may be obsolete.  Values for this quantity are read from a file.

@warning May be obsolete.  Seems to not be used.
@todo Remove if determined to be obsolete.
*/
vector<double> NaTarget_;

/**
@brief List of target potassium concentrations [mol/kgw] at different time steps.

This member may be obsolete.  Values for this quantity are read from a file.

@warning May be obsolete.  Seems to not be used.
@todo Remove if determined to be obsolete.
*/
vector<double> KTarget_;

/**
@brief List of target magnesium concentrations [mol/kgw] at different time steps.

This member may be obsolete.  Values for this quantity are read from a file.

@warning May be obsolete.  Seems to not be used.
@todo Remove if determined to be obsolete.
*/
vector<double> MgTarget_;

/**
@brief List of target sulfate concentrations [mol/kgw] at different time steps.

This member may be obsolete.  Values for this quantity are read from a file.

@warning May be obsolete.  Seems to not be used.
@todo Remove if determined to be obsolete.
*/
vector<double> SO4Target_;

vector<string> name_;           /**< List of names of phases in the kinetic model */
vector<bool> isKinetic_;        /**< List of ids of phases that are kinetically controlled */
vector<bool> isPK_;             /**< List of ids of phases that are PK-model controlled */
vector<bool> isThermo_;         /**< List of ids of phases that are thermodynamically
                                        controlled */
vector<bool> isSoluble_;        /**< List of ids of phases that are instantly dissolved */
vector<int> microPhaseId_;      /**< List of microstructure ids that are in kinetic model */
vector<int> DCId_;              /**< List of DC ids from the ChemicalSystem object */
vector<int> GEMPhaseId_;        /**< List of phase ids from the ChemicalSystem object */
vector<vector<int> > RdICId_;   /**< List of IC ids for each phase */
vector<double> k1_;             /**< List of Parrot and Killoh <i>K</i><sub>1</sub> values */
vector<double> k2_;             /**< List of Parrot and Killoh <i>K</i><sub>2</sub> values */
vector<double> k3_;             /**< List of Parrot and Killoh <i>K</i><sub>3</sub> values */
vector<double> n1_;             /**< List of Parrot and Killoh <i>N</i><sub>1</sub> values */
vector<double> n3_;             /**< List of Parrot and Killoh <i>N</i><sub>3</sub> values */
vector<vector<double> > Rd_;    /**< List of Rd values for each IC in each kinetic phase */
vector<double> scaledMass_;     /**< List of phase mass percents, total solids basis */
vector<double> initScaledMass_; /**< List of initial phase mass percents */
vector<double> activationEnergy_;  /**< List of apparent activation energies for each
                                    kinetic phase [J/mol/K] */
vector<double> critDOH_;        /**< List of critical degrees of hydration for w/c
                                        effect in the Parrot and Killoh model */
vector<double> degreeOfHydration_; /**< Degree of hydration of each kinetic (clinker) phase */

double pfk1_,pfk2_,pfk3_;       /**< Pozzolanic factors for the Parrott and Killoh parameters */

bool verbose_;                  /**< Flag for verbose output */
bool warning_;                  /**< Flag for warnining output */
bool debug_;                    /**< Flag for debugging output */

public:
    
/**
@brief Default constructor.

This constructor is not used in THAMES.  It just establishes default values for
all the member variables.

@note NOT USED.
*/
KineticModel ();

/**
@brief Overloaded constructor.

This constructor is the one invoked by THAMES.  It can only be called once the
various other objects for the simulation are allocated and constructed.

@param cs is a pointer to the ChemicalSystem object for the simulation
@param solut is a pointer to the aqeuous solution object for the simulation
@param lattic is a pointer to the Lattice object holding the microstructure
@param fileName is the name of the XML file with the input for the kinetic model
@param verbose is true if verbose output should be produced
@param warning is false if suppressing warning output
@param debug is true if debugging output should be produced
*/
KineticModel (ChemicalSystem *cs,
              Solution *solut,
              Lattice *lattic,
              const string &fileName,
              const bool verbose,
              const bool warning,
              const bool debug);
     
/**
@brief Master method controlling the parsing of XML input to the kinetic model.

@param docName is the name of the (purported) XML input file
*/
void parseDoc (const string &docName);

/**
@brief Parse the input data for one phase in the XML input file.

Note that this method uses the libxml library, so this must be included.

@param doc is a libxml pointer to the document head
@param cur is a libxml pointer to the current node being parsed
@param numEntry is the number of solid entries in the XML file, will be incremented
@param kineticData is a reference to the KineticData structure for temporarily storing
            the input parameters.
*/
void parsePhase (xmlDocPtr doc,
                 xmlNodePtr cur,
                 int &numEntry,
                 KineticData &kineticData);

/**
@brief Parse the kinetic data for one phase in the XML input file.

Note that this method uses the libxml library, so this must be included.
This method parses the complex grouping of data associated with kinetic
control of a given phase, including the Parrot and Killoh parameters.

@param doc is a libxml pointer to the document head
@param cur is a libxml pointer to the current node being parsed
@param kineticData is a reference to the KineticData structure for temporarily storing
            the input parameters.
*/
void parseKineticData (xmlDocPtr doc,
                       xmlNodePtr cur,
                       KineticData &kineticData);

/**
@brief Parse the Rd data (impurity partitioning) for one phase in the XML input file.

Note that this method uses the libxml library, so this must be included.
This method parses the data about partitioning of impurities within the
kinetically controlled phases.

@param doc is a libxml pointer to the document head
@param cur is a libxml pointer to the current node being parsed
@param kineticData is a reference to the KineticData structure for temporarily storing
            the input parameters
*/
void parseRdData (xmlDocPtr doc,
                  xmlNodePtr cur,
                  KineticData &kineticData);
    
/**
@brief Compute normalized initial microstructure phase masses

Given the initial masses of all phases in the microstructure,
this method scales them to 100 grams of solid.  In the process,
this method also sets the initial moles of water in the
chemical system definition.
*/

void getPhaseMasses (void);

/**
@brief Set the initial total microstructure volume

This method computes the sums of all microstructure phase volumes
and assigns the total.
*/
void setInitialTotalVolume(void);

/**
@brief Set the initial moles of water

This method uses the mass of water and the molar mass
of water as defined in the CSD to calculate the moles
of water in the system.
*/
void setInitialWaterMoles(void);

/**
@brief Get the list of all GEM DC ids in the ChemicalSystem object.

@note NOT USED.

@return the list of all GEM DC ids in the ChemicalSystem object.
*/
vector<int> getDCId () const
{
    return DCId_;
}
  
/**
@brief Get the GEM DC id (in the ChemicalSystem object) for a given kinetic phase.

@note NOT USED.

@param i is the id of the kinetic phase in the kinetic model
@return the GEM phase id of the kinetic phase
*/
int getDCId (const unsigned int i) const
{
    try { return DCId_.at(i); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getDCId",
                           "DCId_",DCId_.size(),i);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the list of all microstructure ids in the KineticModel.

Usually this will be all of the ChemicalSystem microstructure phases
except for VOID and H2O.

@return the list of all microstructure ids.
*/
vector<int> getMicroPhaseId () const
{
    return microPhaseId_;
}
  
/**
@brief Get a microstructure id in the KineticModel.

@param idx is the the index of all the phases in the kinetic model
@return the microstructure id of that phase
*/
int getMicroPhaseId (const int idx)
{
    try { return microPhaseId_.at(idx); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getMicroPhaseId",
                           "microPhaseId_",microPhaseId_.size(),idx);
        ex.printException();
        exit(1);
    }
}
  
/**
@brief Set the total number of phases in the kinetic model.

@note NOT USED.

@param numphases is the total number of phases in the kinetic model
*/
void setNumPhases (const unsigned int numphases)
{
    numPhases_ = numphases;
}

/**
@brief Get the total number of phases in the kinetic model.

@note NOT USED.

@return the total number of phases in the kinetic model
*/
int getNumPhases () const
{
    return numPhases_;
}

/**
@brief Get the kinetic status of every microstructure phase

@note NOT USED.

@return the list of kinetic status for all microstructure phases
*/
vector<bool> isKinetic () const
{
    return isKinetic_;
}

/**
@brief Get the thermo status of every microstructure phase

@note NOT USED.

@return the list of thermo status for all microstructure phases
*/
vector<bool> isThermo () const
{
    return isThermo_;
}

/**
@brief Get the readily soluble status of every microstructure phase

@note NOT USED.

@return the list of readily soluble status for all microstructure phases
*/
vector<bool> isSoluble () const
{
    return isSoluble_;
}

/**
@brief Get the kinetic status of a microstructure phase by its id

@note NOT USED.

@param idx is a microstructure id number
@return true if it is kinetically controlled
*/
bool isKinetic (const unsigned int idx)
{
    try { return isKinetic_.at(idx); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","isKinetic",
                           "isKinetic_",isKinetic_.size(),idx);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the kinetic status of a microstructure phase by its name

@note NOT USED.

@param micname is the name of a microstructure phase
@return true if it is kinetically controlled
*/
bool isKinetic (const string micname)
{
    int idx;
    try {
        idx = chemSys_->getMicroPhaseId(micname);
        return isKinetic_.at(idx);
    }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","isKinetic",
                           "isKinetic_",isKinetic_.size(),idx);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the Parrott and Killoh model status of a microstructure phase by its id

@param idx is a microstructure id number
@return true if it is kinetically controlled by PK model
*/
bool isPK (const unsigned int idx)
{
    try { return isPK_.at(idx); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","isPK",
                           "isPK_",isPK_.size(),idx);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the Parrott and Killoh model status of a microstructure phase by its name

@param micname is the name of a microstructure phase
@return true if it is kinetically controlled by PK model
*/
bool isPK (const string micname)
{
    int idx;
    try {
        idx = chemSys_->getMicroPhaseId(micname);
        return isPK_.at(idx);
    }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","isPK",
                           "isPK_",isPK_.size(),idx);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the thermo status of a microstructure phase by its id

@note NOT USED.

@param idx is a microstructure id number
@return true if it is thermodynamically controlled
*/
bool isThermo (const unsigned int idx)
{
    try { return isThermo_.at(idx); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","isThermo",
                           "isThermo_",isThermo_.size(),idx);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the thermo status of a microstructure phase by its name

@note NOT USED.

@param micname is the name of a microstructure phase
@return true if it is thermodynamically controlled
*/
bool isThermo (const string micname)
{
    int idx;
    try {
        idx = chemSys_->getMicroPhaseId(micname);
        return isThermo_.at(idx);
    }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","isThermo",
                           "isThermo_",isThermo_.size(),idx);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the readily soluble status of a microstructure phase by its id

@note NOT USED.

@param idx is a microstructure id number
@return true if it is readily soluble
*/
bool isSoluble (const unsigned int idx)
{
    try { return isSoluble_.at(idx); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","isSoluble",
                           "isSoluble_",isSoluble_.size(),idx);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the soluble status of a microstructure phase by its name

@note NOT USED.

@param micname is the name of a microstructure phase
@return true if it is readily soluble
*/
bool isSoluble (const string micname)
{
    int idx;
    try {
        idx = chemSys_->getMicroPhaseId(micname);
        return isSoluble_.at(idx);
    }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","isSoluble",
                           "isSoluble_",isSoluble_.size(),idx);
        ex.printException();
        exit(1);
    }
}

/**
@brief Set the w/c mass ratio of the system for the kinetic model equations.

@note NOT USED.

@param wcr is the w/c ratio to set
*/
void setWcRatio (double wcr)
{
    wcRatio_ = wcr;
}

/**
@brief Get the w/c mass ratio of the system used by the kinetic model equations.

@return the w/c mass ratio
*/
double getWcRatio () const
{
    return wcRatio_;
}

/**
@brief Set the initial total solid mass

@note NOT USED

@param ism is the initial solid mass [g]
*/
void setInitSolidMass (double ism)
{
    initSolidMass_ = ism;
}

/**
@brief Get the initial solid mass

@note NOT USED

@return the initial solid mass [g]
*/
double getInitSolidMass () const
{
    return initSolidMass_;
}

/**
@brief Set the Blaine fineness of the cement.

@note NOT USED.

@param bval is the Blaine fineness [m<sup>2</sup>/kg]
*/
void setBlaine (double bval)
{
    blaine_ = bval;
}

/**
@brief Get the Blaine fineness of the cement.

@note NOT USED.

@return the Blaine fineness [m<sup>2</sup>/kg]
*/
double getBlaine () const
{
    return blaine_;
}

/**
@brief Set the reference Blaine fineness parameter for the Parrot and Killoh model.

The value set in the Parrot and Killoh model is 385 m<sup>2</sup>/kg, and there
is no particular reason to change it.

@note NOT USED.

@param rbval is the reference Blaine fineness [m<sup>2</sup>/kg]
*/
void setRefBlaine (double rbval)
{
    refBlaine_ = rbval;
}

/**
@brief Get the reference Blaine fineness parameter for the Parrot and Killoh model.

@note NOT USED.

@return the reference Blaine fineness [m<sup>2</sup>/kg]
*/
double getRefBlaine () const
{
    return refBlaine_;
}

/**
@brief Set the ratio of the true Blaine fineness to the model reference value.

@note NOT USED.

@param bfact is the ratio of the actual Blaine fineness to the reference value
*/
void setBlaineFactor (double bfact)
{
    blaineFactor_ = bfact;
}

/**
@brief Get the ratio of the true Blaine fineness to the model reference value.

@note NOT USED.

@return the ratio of the actual Blaine fineness to the reference value
*/
double getBlaineFactor () const
{
    return blaineFactor_;
}

/**
@brief Set the absolute temperature.

@note NOT USED.

@param tval is the absolute temperature [K]
*/
void setTemperature (double tval)
{
    temperature_ = tval;
}

/**
@brief Get the absolute temperature.

@note NOT USED.

@return the absolute temperature [K]
*/
double getTemperature () const
{
    return temperature_;
}

/**
@brief Set the model reference temperature.

@note NOT USED.

@param rtval is the reference temperature [K]
*/
void setRefT (double rtval)
{
    refT_ = rtval;
}

/**
@brief Get the model reference temperature.

@note NOT USED.

@return the reference temperature [K]
*/
double getRefT () const
{
    return refT_;
}

/**
@brief Initialize a kinetic data structure 
*/
void initKineticData(struct KineticData &kd)
{
    kd.name.clear();
    kd.microPhaseId = 0;
    kd.GEMPhaseId = 0;
    kd.DCId = 0;
    kd.type.clear();
    kd.scaledMass = 0.0;
    kd.k1 = 0.0;
    kd.k2 = 0.0;
    kd.k3 = 0.0;
    kd.n1 = 0.0;
    kd.n3 = 0.0;
    kd.Ea = 0.0;
    kd.critDOH = 0.0;
    kd.RdId.clear();
    kd.RdVal.clear();
    return;
}

/**
@brief Get the ChemicalSystem object for the simulation used by the kinetic model.

@note NOT USED.

@return a pointer to the ChemicalSystem object
*/
ChemicalSystem *getChemSys () const
{
    return chemSys_;
}

/**
@brief Set the simulation time at which to begin external sulfate attack.

@param sattacktime is the simulation time to begin sulfate attack [days]
*/
void setSulfateAttackTime (double sattacktime)
{
    sulfateAttackTime_ = sattacktime;
}

/**
@brief Get the simulation time at which to begin external sulfate attack.

@note NOT USED.

@return the simulation time to begin sulfate attack [days]
*/
double getSulfateAttackTime (void) const
{
    return sulfateAttackTime_;
}

/**
@brief Set the simulation time at which to begin leaching.

@param leachtime is the simulation time to begin leaching [days]
*/
void setLeachTime (double leachtime)
{
    leachTime_ = leachtime;
}

/**
@brief Get the simulation time at which to begin leaching.

@note NOT USED.

@return the simulation time to begin leaching [days]
*/
double getLeachTime (void) const
{
    return leachTime_;
}

/**
@brief Get the list of phase names used by the kinetic model.

@note NOT USED.

@return the vector of names of phases in the kinetic model
*/
vector<string> getName () const
{
    return name_;
}

/**
@brief Get the name of phase with a given index in the kinetic model.

@param i is the index of the phase in the kinetic model
@return the name of the phase with index i
*/
string getName (const unsigned int i) const
{
    try { return name_.at(i); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getName","name_",
                           name_.size(),i);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the list of <i>K</i><sub>1</sub> values for clinker phases in the PK model.

@note NOT USED.

@return the vector of <i>K</i><sub>1</sub> values for clinker phases in the PK model
*/
vector<double> getK1 () const
{
    return k1_;
}

/**
@brief Get the <i>K</i><sub>1</sub> value for a particular clinker phase in the PK model.

@note NOT USED.

@param i is the index of the clinker phase in the kinetic model
@return the <i>K</i><sub>1</sub> value for the clinker phase in the PK model
*/
double getK1 (const unsigned int i) const
{
    try { return k1_.at(i); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getK1","k1_",k1_.size(),i);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the list of <i>K</i><sub>2</sub> values for clinker phases in the PK model.

@note NOT USED.

@return the vector of <i>K</i><sub>2</sub> values for clinker phases in the PK model
*/
vector<double> getK2 () const
{
    return k2_;
}

/**
@brief Get the <i>K</i><sub>2</sub> value for a particular clinker phase in the PK model.

@note NOT USED.

@param i is the index of the clinker phase in the kinetic model
@return the <i>K</i><sub>2</sub> value for the clinker phase in the PK model
*/
double getK2 (const unsigned int i) const
{
    try { return k2_.at(i); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getK2","k2_",k2_.size(),i);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the list of <i>K</i><sub>3</sub> values for clinker phases in the PK model.

@note NOT USED.

@return the vector of <i>K</i><sub>3</sub> values for clinker phases in the PK model
*/
vector<double> getK3 () const
{
    return k3_;
}

/**
@brief Get the <i>K</i><sub>3</sub> value for a particular clinker phase in the PK model.

@note NOT USED.

@param i is the index of the clinker phase in the kinetic model
@return the <i>K</i><sub>3</sub> value for the clinker phase in the PK model
*/
double getK3 (const unsigned int i) const
{
    try { return k3_.at(i); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getK3","k3_",k3_.size(),i);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the list of <i>N</i><sub>1</sub> values for clinker phases in the PK model.

@note NOT USED.

@return the vector of <i>N</i><sub>1</sub> values for clinker phases in the PK model
*/
vector<double> getN1 () const
{
    return n1_;
}

/**
@brief Get the <i>N</i><sub>1</sub> value for a particular clinker phase in the PK model.

@note NOT USED.

@param i is the index of the clinker phase in the kinetic model
@return the <i>N</i><sub>1</sub> value for the clinker phase in the PK model
*/
double getN1 (const unsigned int i) const
{
    try { return n1_.at(i); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getN1","n1_",n1_.size(),i);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the list of <i>N</i><sub>3</sub> values for clinker phases in the PK model.

@note NOT USED.

@return the vector of <i>N</i><sub>3</sub> values for clinker phases in the PK model
*/
vector<double> getN3 () const
{
    return n3_;
}

/**
@brief Get the <i>N</i><sub>3</sub> value for a particular clinker phase in the PK model.

@note NOT USED.

@param i is the index of the clinker phase in the kinetic model
@return the <i>N</i><sub>3</sub> value for the clinker phase in the PK model
*/
double getN3 (const unsigned int i) const
{
    try { return n3_.at(i); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getN3","n3_",n3_.size(),i);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the list of scaled masses for the phases in the kinetic model.

The scaled mass of a phase is its mass percent on a total solids basis.

@note NOT USED.

@return the vector of scaled masses [percent solids]
*/
vector<double> getScaledMass () const
{
    return scaledMass_;
}

/**
@brief Get the scaled mass of a particular clinker phase in the kinetic model.

The scaled mass of a phase is its mass percent on a total solids basis.

@note NOT USED.

@param i is the index of the phase in the kinetic model
@return the scaled mass of the phase [percent solids]
*/
double getScaledMass (const unsigned int i) const
{
    try { return scaledMass_.at(i); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getScaledMass",
                           "scaledMass_",scaledMass_.size(),i);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the list of scaled <i>initial</i> masses for the phases in the kinetic model.

The scaled mass of a phase is its mass percent on a total solids basis.

@return the vector of initial scaled masses [percent solids]
*/
vector<double> getInitScaledMass () const
{
    return initScaledMass_;
}

/**
@brief Get the initial scaled mass of a particular clinker phase in the kinetic model.

The scaled mass of a phase is its mass percent on a total solids basis.

@note NOT USED.

@param i is the index of the phase in the kinetic model
@return the initial scaled mass of the phase [percent solids]
*/
double getInitScaledMass (const unsigned int i) const
{
    try { return initScaledMass_.at(i); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getInitScaledMass",
                           "initScaledMass_",initScaledMass_.size(),i);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the list of activation energies for the phases in the kinetic model.

@note NOT USED.

@return the vector of activation energies [J/mol/K]
*/
vector<double> getActivationEnergy () const
{
    return activationEnergy_;
}

/**
@brief Get the activation energy of a clinker phase in the kinetic model.

@note NOT USED.

@param i is the index of the phase in the kinetic model
@return the activation energy of the phase [J/mol/K]
*/
double getActivationEnergy (const unsigned int i) const
{
    try { return activationEnergy_.at(i); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getActivationEnergy",
                        "activationEnergy_",activationEnergy_.size(),i);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the list of critical degrees of hydration for w/c effects in the kinetic model.

@note NOT USED.

@return the vector of critical degrees of hydration for the phases
*/
vector<double> getCritDOH () const
{
    return critDOH_;
}

/**
@brief Get the critical degree of hydration for a particular clinker phase in the kinetic model.

@note NOT USED.

@param i is the index of the phase in the kinetic model
@return the critical degree of hydration of the phase
*/
double getCritDOH (const unsigned int i) const
{
    try { return critDOH_.at(i); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getCritDOH",
                           "critDOH_",critDOH_.size(),i);
        ex.printException();
        exit(1);
    }
}

/**
@brief Set the "pozzolanic factor" for the Parrott and Killoh k1 parameter

@param pk1 is the value to set for the k1 parameter's pozzolanic factor
*/
void setPfk1 (const double pk1)
{
    pfk1_ = pk1;
}

/**
@brief Get the "pozzolanic factor" for the Parrott and Killoh k1 parameter

@return the "pozzolanic factor" for the Parrott and Killoh k1 parameter
*/
double getPfk1 () const
{
    return pfk1_;
}

/**
@brief Set the "pozzolanic factor" for the Parrott and Killoh k2 parameter

@param pk2 is the value to set for the k2 parameter's pozzolanic factor
*/
void setPfk2 (const double pk2)
{
    pfk2_ = pk2;
}

/**
@brief Get the "pozzolanic factor" for the Parrott and Killoh k2 parameter

@return the "pozzolanic factor" for the Parrott and Killoh k2 parameter
*/
double getPfk2 () const
{
    return pfk2_;
}

/**
@brief Set the "pozzolanic factor" for the Parrott and Killoh k3 parameter

@param pk3 is the value to set for the k3 parameter's pozzolanic factor
*/
void setPfk3 (const double pk3)
{
    pfk3_ = pk3;
}

/**
@brief Get the "pozzolanic factor" for the Parrott and Killoh k3 parameter

@return the "pozzolanic factor" for the Parrott and Killoh k3 parameter
*/
double getPfk3 () const
{
    return pfk3_;
}

/**
@brief Get the list of degrees of hydration for clinker phases in the kinetic model.

@note NOT USED.

@return the vector of degrees of hydration of the clinker phases
*/
vector<double> getDegreeOfHydration () const
{
    return degreeOfHydration_;
}

/**
@brief Get the degree of hydration for a particular clinker phase in the kinetic model.

@note NOT USED.

@param i is the index of the phase in the kinetic model
@return the degree of hydration of the phase
*/
double getDegreeOfHydration (const unsigned int i) const
{
    try { return degreeOfHydration_.at(i); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getDegreeOfHydration","degreeOfHydration_",
                           degreeOfHydration_.size(),i);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the list of independent components (ICs) for each phase in the kinetic model.

@note NOT USED.

@return the vector of lists of ICs for the phases in the kinetic model
*/
vector<vector<int> > getRdICId () const
{
    return RdICId_;
}

/**
@brief Get the list of Rd values for the independent components (ICs) for each phase in the kinetic model.

The Rd values for each phase are the partionings of impurities in the clinker phases.

@note NOT USED.

@return the vector of Rd values of the ICs for the phases in the kinetic model
*/
vector<vector<double> > getRd () const
{
    return Rd_;
}

/**
@brief Get the ids of all ICs of a particular phase in the kinetic model.

@note NOT USED.

@param idx is the index of the phase in the kinetic model
@return the list of IC ids in that phase.
*/
vector<int> getRdICId (const unsigned int idx)
{
    try { return RdICId_.at(idx); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getRdICId",
                           "RdICId_",RdICId_.size(),idx);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the list of Rd values for the independent components (ICs) in a phase in the kinetic model.

The Rd values for each phase are the partionings of impurities in the clinker phases.

@note NOT USED.

@param idx is the index of the phase in the kinetic model
@return the vector of Rd values of the ICs for the phase specified by idx
*/
vector<double> getRd (const unsigned int idx)
{
    try { return Rd_.at(idx); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getRd",
                           "Rd_",Rd_.size(),idx);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the id of one IC of a particular phase in the kinetic model.

@note NOT USED.

@param idx1 is the index of the phase in the kinetic model
@param idx2 is the element location of the IC for that phase
@return the IC id value
*/
unsigned int getRdICId (const unsigned int idx1,
                        const unsigned int idx2)
{
    try { return RdICId_.at(idx1).at(idx2); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getRdICId",
                           "RdICId_",RdICId_.size(),idx1);
        ex.printException();
        exit(1);
    }
}

/**
@brief Get the Rd of one IC of a particular phase in the kinetic model.

The Rd values for each phase are the partionings of impurities in the clinker phases.

@note NOT USED.

@param idx1 is the index of the phase in the kinetic model
@param idx2 is the element location of the IC for that phase
@return the IC id value
*/
double getRd (const unsigned int idx1,
              const unsigned int idx2)
{
    try { return Rd_.at(idx1).at(idx2); }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","getRd",
                           "Rd_",Rd_.size(),idx1);
        ex.printException();
        exit(1);
    }
}
  
/**
@brief Master method for implementing one kinetic time step.

In a given time step, a certain number of moles of each clinker phase will dissolve,
and the instantly soluble phases will dissolve in the first time step.  This function
determines the number of moles of each phase to dissolve, based on the time interval
being simulated.  It then calculates the number of IC moles to promote to the thermodynamic
system from those phases (which are outside the thermodynamic system because they
are kinetically controlled), based on the stoichiometry.  Those IC moles are then
added to the thermodynamic system, and the moles and mass of each kinetically controlled
phase are changed accordingly.

@remark This method is very long and several parts are hard-coded when they
should be made more general.

@todo Split this method into more convenient chunks
@todo Make the methods more general, less hardwiring of parameters
@todo Make the local variable names more descriptive

@param timestep is the time interval to simulate [days]
@param temperature is the absolute temperature during this step [K]
@param isFirst is true if this is the first time step of the simulation, false otherwise
*/
void calculateKineticStep (const double timestep,
                           const double temperature,
                           bool isFirst);
     
/**
@brief Determine the change in moles of a given kinetically controlled phase.

This method is called by the `calculateKineticStep` method, to calculate the change
in moles of a given kinetically controlled phase during a given time interval.

@note NOT USED.

@todo Generalize the rate equation for other phases more than it is.
@todo Change the variable names to be more descriptive

@param pid is the id of the microstructure phase to change
@param k is the effective rate constant in the rate equation
@param gamma is the exponent for the driving force term, (SI - 1)
@param timestep is the time interval to simulate [days]
*/
void calculatePhaseChange (int pid,
                           double k,
                           double gamma,
                           double timestep);
 
/**
@brief Set up the number of moles of dependent components in the kinetic phases.

This method loops over the <i>kinetically</i> controlled phases in the kinetic
model, gets the DC stoichiometry of each phase, and determines the number of moles
of each DC component based on the number of moles of the kinetically controlled phases.
*/
void setKineticDCMoles ();

/**
@brief Set the number of moles of dependent components to zero.

This method loops over the <i>kinetically</i> controlled phases in the kinetic
model, and sets the number of moles of each DC component of that phase to zero.
*/
void zeroKineticDCMoles ();

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
bool getVerbose () const
{
    return verbose_;
}

/**
@brief Set the warning flag

@param iswarning is true if verbose output should be produced
*/
void setWarning (const bool iswarning)
{
    warning_ = iswarning;
    return;
}

/**
@brief Get the warning flag

@return the warning flag
*/
bool getWarning () const
{
    return warning_;
}

/**
@brief Set the debug flag

@param isdebug is true if debugging output should be produced
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
bool getDebug () const
{
    return debug_;
}

};      // End of KineticModel class

#endif
