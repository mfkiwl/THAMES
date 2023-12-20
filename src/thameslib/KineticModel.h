/**
@file KineticModel.h
@brief Declaration of the KineticModel class.

@section Introduction

In THAMES, the `KineticModel` class can be perceived as the engine that calculates
the kinetic changes in the system during a given time increment.  The primary
kinetic aspect that is calculated is the extent of dissolution of mineral phases in
the original cement.  This is just the base class.  It is not used.

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
// #include "KineticController.h"
#include "Lattice.h"
#include "KineticData.h"
#include "global.h"

using namespace std;

/**
@class KineticModel
@brief Base class for all kinetic models

THAMES allows some flexibility in defining different types of kinetic models.
*/

class KineticModel {

protected:
    
int numPhases_;             /**< Total number of phases in the kinetic model */
ChemicalSystem *chemSys_;   /**< Pointer to the ChemicalSystem object for this simulation */
Solution *solut_;           /**< Pointer to the aqueous phase object for this simulation */
Lattice *lattice_;          /**< Pointer to the lattice object holding the microstructure */
double initSolidMass_;      /**< initial total mass of solids controlled by this model [g] */
double temperature_;        /**< Temperature [K] */
double refT_;               /**< Reference temperature for PK model [K] */
double sulfateAttackTime_;  /**< Time at which sulfate attack simulation starts [days] */
double leachTime_;          /**< Time at which leaching simulation starts [days] */

string name_;                   /**< Name of phase controlled by this kinetic model */
int microPhaseId_;              /**< Microstructure id controlled by this model */
int DCId_;                      /**< List of DC ids from the ChemicalSystem object */
int GEMPhaseId_;                /**< List of phase ids from the ChemicalSystem object */
int waterId_;                   /**< DC index for liquid water */
int ICNum_;                     /**< Number of ICs in chemical system */
int DCNum_;                     /**< Number of DCs in chemical system */
int GEMPhaseNum_;               /**< Number of GEM phases in chemical system */
vector<string> ICName_;         /**< Names of ICs */
vector<string> DCName_;         /**< Names of DCs */
double scaledMass_;             /**< Phase mass percent, total solids basis */
double initScaledMass_;         /**< Initial phase mass percents */
double activationEnergy_;       /**< Apparent activation energy for the reaction [J/mol/K] */
double specificSurfaceArea_;    /**< Specific surface area (m2/kg) */
double refSpecificSurfaceArea_; /**< Reference specific surface area (m2/kg) */
double ssaFactor_;              /**< Reference specific surface area (m2/kg) */
double degreeOfReaction_;       /**< Degree of reaction of this component (mass
basis) */
double lossOnIgnition_;         /**< Loss on ignition of this component (ignited
mass basis) */
bool verbose_;                  /**< Flag for verbose output */
bool warning_;                  /**< Flag for warnining output */

public:
    
/**
@brief Default constructor.

This constructor is not used in THAMES.  It just establishes default values for
all the member variables.

@note NOT USED.
*/
KineticModel ();

     
/**
@brief Destructor does nothing.
*/
virtual ~KineticModel() {};

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
@brief Set the multiplicative adjustment to clinker phase rate constant due to
pozzolans

@param pfk is the specific surface area [m<sup>2</sup>/kg]
*/
virtual void setPfk (double pfk)
{
    // Most components do not have a pfk variable, generically return
    return;
}

/**
@brief Get the type of kinetic model

@return a string indicating the model type
*/
virtual string getType () const
{
    return (GenericType);
}

/**
@brief Set the specific surface area

@note NOT USED.

@param sval is the specific surface area [m<sup>2</sup>/kg]
*/
void setSpecificSurfaceArea (double sval)
{
    specificSurfaceArea_ = sval;
}

/**
@brief Get the specific surface area.

@note NOT USED.

@return the specific surface area fineness [m<sup>2</sup>/kg]
*/
double getSpecificSurfaceArea () const
{
    return specificSurfaceArea_;
}

/**
@brief Set the reference specific surface area.

The value set in the Parrot and Killoh model is 385 m<sup>2</sup>/kg, and there
is no particular reason to change it.

@note NOT USED.

@param rsval is the reference specific surface area [m<sup>2</sup>/kg]
*/
void setRefSpecificSurfaceArea (const double rsval)
{
    refSpecificSurfaceArea_ = rsval;
}

/**
@brief Get the reference specific surface area.

@note NOT USED.

@return the reference specific surface area [m<sup>2</sup>/kg]
*/
double getRefSpecificSurfaceArea () const
{
    return refSpecificSurfaceArea_;
}

/**
@brief Set the ratio of the true specific surface area to the model reference value.

@note NOT USED.

@param sfact is the ratio of the actual specific surface area to the reference value
*/
void setSsaFactor (const double sfact)
{
    ssaFactor_ = sfact;
}

/**
@brief Get the ratio of the true specific surface area to the model reference value.

@note NOT USED.

@return the ratio of the actual specific surface area to the reference value
*/
double getSsaFactor () const
{
    return ssaFactor_;
}

/**
@brief Set the degree of reaction of this component (mass basis)

@note NOT USED.

@param dor is the degree of reaction to set
*/
virtual void setDegreeOfReaction (const double dor)
{
    degreeOfReaction_ = dor;
}

/**
@brief Get the degree of reaction of this component (mass basis)

@note NOT USED.

@return the degree of reaction of this component
*/
virtual double getDegreeOfReaction () const
{
    return degreeOfReaction_;
}

/**
@brief Set the loss on ignition of this component (ignited mass basis)

@note NOT USED.

@param loi is the degree of reaction to set
*/
virtual void setLossOnIgnition (const double loi)
{
    lossOnIgnition_ = loi;
}

/**
@brief Get the loss on ignition of this component (ignited mass basis)

@note NOT USED.

@return the loss on ignition of this component
*/
virtual double getLossOnIgnition () const
{
    return lossOnIgnition_;
}

/**
@brief Compute normalized initial microstructure phase masses

Given the initial masses of all phases in the microstructure,
this method scales them to 100 grams of solid.  In the process,
this method also sets the initial moles of water in the
chemical system definition.
*/

void getPhaseMasses (void);

/**
@brief Get the microstructure id in the KineticModel.


@return the list of all microstructure ids.
*/
int getMicroPhaseId () const
{
    return microPhaseId_;
}
  
/**
@brief Set the DC index for liquid water

@note NOT USED.

@param waterid is the index value to use
*/
void setWaterId (const int waterid)
{
    waterId_ = waterid;
}

/**
@brief Get the DC index for liquid water

@return the DC index for liquid water
*/
int getWaterId () const
{
    return waterId_;
}
  
/**
@brief Set the number of ICs

@note NOT USED.

@param icnum is the number of ICs to specify
*/
void setICNum (const int icnum)
{
    ICNum_ = icnum;
}

/**
@brief Get the number of ICs

@return the number of ICs
*/
int getICNum () const
{
    return ICNum_;
}
  
/**
@brief Set the number of DCs

@note NOT USED.

@param dcnum is number of DCs to specify
*/
void setDCNum (const int dcnum)
{
    DCNum_ = dcnum;
}

/**
@brief Get the number of DCs

@return the number of DCs
*/
int getDCNum () const
{
    return DCNum_;
}
  
/**
@brief Set the number of GEM phases 

@note NOT USED.

@param gpnum is number of GEM phases to specify
*/
void setGEMPhaseNum (const int gpnum)
{
    GEMPhaseNum_ = gpnum;
}

/**
@brief Get the number of GEM phases

@return the number of GEM phases
*/
int getGEMPhaseNum () const
{
    return GEMPhaseNum_;
}
  
/**
@brief Set the IC names

@note NOT USED.

@param icname is the list of IC names
*/
void setICName (vector<string> icname)
{
    ICName_ = icname;
}

/**
@brief Get the IC names

@return the list of IC names
*/
vector<string> getICName () const
{
    return ICName_;
}
  
/**
@brief Set the DC names

@note NOT USED.

@param dcname is the list of DC names
*/
void setDCName (vector<string> dcname)
{
    DCName_ = dcname;
}

/**
@brief Get the DC names

@return the list of DC names
*/
vector<string> getDCName () const
{
    return DCName_;
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
string getName () const
{
    return name_;
}

/**
@brief Get the list of activation energies for the phases in the kinetic model.

@note NOT USED.

@return the vector of activation energies [J/mol/K]
*/
double getActivationEnergy () const
{
    return activationEnergy_;
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
@brief Get the scaled mass of the phase in the kinetic model.

The scaled mass of a phase is its mass percent on a total solids basis.

@note NOT USED.

@return the vector of scaled masses [percent solids]
*/
double getScaledMass () const
{
    return scaledMass_;
}

/**
@brief Get the <i>initial</i> mass of the phase in the kinetic model.

The scaled mass of a phase is its mass percent on a total solids basis.

@return the vector of initial scaled masses [percent solids]
*/
double getInitScaledMass () const
{
    return initScaledMass_;
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

This is now a pure virtual function.

@remark This method is very long and several parts are hard-coded when they
should be made more general.

@todo Split this method into more convenient chunks
@todo Make the methods more general, less hardwiring of parameters
@todo Make the local variable names more descriptive

@param timestep is the time interval to simulate [days]
@param temperature is the absolute temperature during this step [K]
@param isFirst is true if this is the first time step of the simulation, false otherwise
@param rh is the internal relative humidity
@param dICMoles is the vector of moles of each IC change due to kinetics
@param dsolutICMoles is the vector of moles of each IC change in solution due to
kinetics
@param DCMoles is the vector of moles of each DC
@param GEMPhaseMoles is the vector of moles of each phase in GEMS
*/
virtual void calculateKineticStep (const double timestep,
                                   const double temperature,
                                   bool isFirst,
                                   double rh,
                                   vector<double> &dICMoles,
                                   vector<double> &dsolutICMoles,
                                   vector<double> &DCMoles,
                                   vector<double> &GEMPhaseMoles) = 0;
     
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

};      // End of KineticModel class

#endif
