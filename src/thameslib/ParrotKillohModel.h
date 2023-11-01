/**
@file ParrotKillohModel.h
@brief Declaration of the ParrotKillohModel class.

@section Introduction
This class implements the Parrot and Killoh (PK) model of
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

#ifndef PARROTKILLOHH
#define PARROTKILLOHH

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
#include "global.h"

using namespace std;

/**
@class ParrotKillohModel
@brief Handles the Parrot and Killoh(1984) kinetic model of clinker ractions

The Parrot and Killoh model [1] is used to empirically estimate the
mass fraction of each clinker phase that dissolves in a unit time.  Eventually
this can be expanded to handle other kinetically controlled phases outside the
Parrot and Killoh model, such as the growth of portlandite or C--S--H.
*/

class ParrotKillohModel : public KineticModel {

protected:
    
double wcRatio_;            /**< water-cement mass ratio */
double blaine_;             /**< Blaine fineness [m<sup>2</sup>/kg] */
double refBlaine_;          /**< Reference Blaine value for PK model,
                                    usually 385.0 m<sup>2</sup>/kg */
double blaineFactor_;       /**< `blaine_`/`refBlaine_` */

vector<double> k1_;             /**< List of Parrot and Killoh <i>K</i><sub>1</sub> values */
vector<double> k2_;             /**< List of Parrot and Killoh <i>K</i><sub>2</sub> values */
vector<double> k3_;             /**< List of Parrot and Killoh <i>K</i><sub>3</sub> values */
vector<double> n1_;             /**< List of Parrot and Killoh <i>N</i><sub>1</sub> values */
vector<double> n3_;             /**< List of Parrot and Killoh <i>N</i><sub>3</sub> values */
vector<double> critDOH_;        /**< List of critical degrees of hydration for w/c
                                        effect in the Parrot and Killoh model */
vector<double> degreeOfHydration_; /**< Degree of hydration of each kinetic (clinker) phase */

public:
    
/**
@brief Default constructor.

This constructor is not used in THAMES.  It just establishes default values for
all the member variables.

@note NOT USED.
*/
ParrotKillohModel ();

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
*/
ParrotKillohModel (ChemicalSystem *cs,
                   Solution *solut,
                   Lattice *lattic,
                   const string &fileName,
                   const bool verbose,
                   const bool warning);
     
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
@brief Master method for implementing one kinetic time step.

Overloaded from base class to handle Parrot and Killoh model.

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
 
};      // End of ParrotKillohModel class

#endif
