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
#include "KineticController.h"
#include "Lattice.h"
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

vector<string> name_;           /**< List of names of phases in the kinetic model */
vector<int> microPhaseId_;      /**< List of microstructure ids that are controlled by this model */
vector<int> DCId_;              /**< List of DC ids from the ChemicalSystem object */
vector<int> GEMPhaseId_;        /**< List of phase ids from the ChemicalSystem object */
vector<vector<int> > RdICId_;   /**< List of IC ids for each phase */
vector<vector<double> > Rd_;    /**< List of Rd values for each IC in each kinetic phase */
vector<double> scaledMass_;     /**< List of phase mass percents, total solids basis */
vector<double> initScaledMass_; /**< List of initial phase mass percents */
vector<double> activationEnergy_;  /**< List of apparent activation energies for each
                                    kinetic phase [J/mol/K] */
vector<bool> isKinetic_;
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
*/
virtual void calculateKineticStep (const double timestep,
                                   const double temperature,
                                   bool isFirst) = 0;
     
/**
@brief Calculate the change in phase moles

Pure virtual function must be defined for each child class

@param microPhaseId is the identity of the phase that is changing
@param rateconst is the rate constant (mol/m2/s)
@param dfexp is the exponent on the driving force term (1 - beta)
@param timestep is the time interval over which to calculate the change
*/
virtual void calculatePhaseChange (const int microPhaseId,
                           double rateconst,
                           double dfexp,
                           double timestep) = 0;

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
