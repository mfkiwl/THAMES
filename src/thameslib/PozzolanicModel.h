/**
@file PozzolanicModel.h
@brief Declaration of the PozzolanicModel class.

@section Introduction
This class implements a dissolution equation like Dove's for 
amorphous silicates [1].

@section References

    -# PM Dove, N Han, AF Wallace, JJ De Yoreo, Kinetics of amorphous silica dissolution and the paradox of the silica polymorphs, Proceedings of the National Academy of Sciences USA, 105 (2008) 9903–9908.

*/

#ifndef POZZOLANICMODELH
#define POZZOLANICMODELH

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <ctime>
#include "KineticController.h"
#include "KineticModel.h"
#include "ChemicalSystem.h"
#include "Lattice.h"
#include "KineticData.h"
#include "global.h"

using namespace std;

/**
@class PozzolanicModel
@brief Handles the Parrot and Killoh(1984) kinetic model of clinker ractions

The Parrot and Killoh model [1] is used to empirically estimate the
mass fraction of each clinker phase that dissolves in a unit time.  Eventually
this can be expanded to handle other kinetically controlled phases outside the
Parrot and Killoh model, such as the growth of portlandite or C--S--H.
*/

class PozzolanicModel : public KineticModel {

protected:
    
double rateconst_;               /**< Rate constant for reaction (mol/m2/s) */
double siexp_;                   /**< Exponent on saturation index (unitless) */
double dfexp_;                   /**< Exponent on driving force (unitless) */
double ohexp_;                   /**< Exponent on OH ion activity (unitless) */
double sio2_;                    /**< Mass fraction of SiO2 (unitless) */
double al2o3_;                   /**< Mass fraction of Al2O3 (unitless) */
double cao_;                     /**< Mass fraction of CaO (unitless) */
double loi_;                     /**< Loss on ignition (unitless) */

public:
    
/**
@brief Default constructor.

This constructor is not used in THAMES.  It just establishes default values for
all the member variables.

@note NOT USED.
*/
PozzolanicModel ();

/**
@brief Overloaded constructor.

This constructor is the one invoked by THAMES.  It can only be called once the
various other objects for the simulation are allocated and constructed.

@param cs is a pointer to the ChemicalSystem object for the simulation
@param solut is a pointer to the aqeuous solution object for the simulation
@param lattice is a pointer to the Lattice object holding the microstructure
@param kineticData is the collection of kinetic parameters already stored
@param verbose is true if verbose output should be produced
@param warning is false if suppressing warning output
*/
PozzolanicModel (ChemicalSystem *cs,
                 Solution *solut,
                 Lattice *lattice,
                 struct KineticData &kineticData,
                 const bool verbose,
                 const bool warning);
     
/**
@brief Set the rate constant

@note NOT USED.

@param rc is the rate constant value to use
*/
void setRateconstant (const double rc)
{
    rateconst_ = max(rc,0.0);
}

/**
@brief Get the rate constant

@note NOT USED.

@return the rate constant
*/
double getRateconst () const
{
    return rateconst_;
}

/**
@brief Set the exponent on the saturation index

@note NOT USED.

@param siexp is the exponent value to use
*/
void setSiexp (const double siexp)
{
    siexp_ = max(siexp,0.0);
}

/**
@brief Get the exponent on the saturation index

@note NOT USED.

@return the exponent on the saturation index
*/
double getSiexp () const
{
    return siexp_;
}

/**
@brief Set the exponent on the driving force

@note NOT USED.

@param dfexp is the exponent value to use
*/
void setDfexp (const double dfexp)
{
    dfexp_ = max(dfexp,0.0);
}

/**
@brief Get the exponent on the driving force

@note NOT USED.

@return the exponent on the driving force
*/
double getDfexp () const
{
    return dfexp_;
}

/**
@brief Set the exponent on the hydroxyl ion activity

@note NOT USED.

@param ohexp is the exponent value to use
*/
void setOhexp (const double ohexp)
{
    ohexp_ = max(ohexp,0.0);
}

/**
@brief Get the exponent on the hydroxyl ion activity

@note NOT USED.

@return the exponent on the hydroxyl ion activity
*/
double getOhexp () const
{
    return ohexp_;
}

/**
@brief Set the mass fraction of SiO2

@note NOT USED.

@param sio2 is the mass fraction of SiO2 
*/
void setSio2 (const double sio2)
{
    sio2_ = max(sio2,0.0);
    if (sio2_ > 1.0) sio2_ = 1.0;
}

/**
@brief Get the SiO2 mass fraction

@note NOT USED.

@return the SiO2 mass fraction 
*/
double getSio2 () const
{
    return sio2_;
}

/**
@brief Set the mass fraction of Al2O3

@note NOT USED.

@param al2o3 is the mass fraction of Al2O3
*/
void setAl2o3 (const double al2o3)
{
    al2o3_ = max(al2o3,0.0);
    if (al2o3_ > 1.0) al2o3_ = 1.0;
}

/**
@brief Get the Al2O3 mass fraction

@note NOT USED.

@return the Al2O3 mass fraction 
*/
double getAl2o3 () const
{
    return al2o3_;
}

/**
@brief Set the mass fraction of CaO

@note NOT USED.

@param cao is the mass fraction of CaO
*/
void setCao (const double cao)
{
    cao_ = max(cao,0.0);
    if (cao_ > 1.0) cao_ = 1.0;
}

/**
@brief Get the CaO mass fraction

@note NOT USED.

@return the CaO mass fraction 
*/
double getCao () const
{
    return cao_;
}

/**
@brief Set the loss on ignition

@note NOT USED.

@param loi is the loss on ignition
*/
void setLoi (const double loi)
{
    loi_ = max(loi,0.0);
    if (loi_ > 1.0) loi_ = 1.0;
}

/**
@brief Get the loss on ignition

@note NOT USED.

@return the loss on ignition
*/
double getLoi () const
{
    return loi_;
}

/**
@brief Master method for implementing one kinetic time step.

Overloaded from base class to handle pozzolanic materials

@todo Split this method into more convenient chunks
@todo Make the methods more general, less hardwiring of parameters
@todo Make the local variable names more descriptive

@param timestep is the time interval to simulate [days]
@param temperature is the absolute temperature during this step [K]
@param isFirst is true if this is the first time step of the simulation, false otherwise
@param doTweak is true if trying to recover from failed convergence
@param rh is the internal relative humidity
@param ICMoles is the vector of moles of each IC
@param solutICMoles is the vector of moles of each IC in solution
@param DCMoles is the vector of moles of each DC
@param GEMPhaseMoles is the vector of moles of each phase in GEMS
*/
virtual void calculateKineticStep (const double timestep,
                                   const double temperature,
                                   bool isFirst,
                                   bool doTweak,
                                   double rh,
                                   vector<double> &ICMoles,
                                   vector<double> &solutICMoles,
                                   vector<double> &DCMoles,
                                   vector<double> &GEMPhaseMoles);

};      // End of PozzolanicModel class

#endif
