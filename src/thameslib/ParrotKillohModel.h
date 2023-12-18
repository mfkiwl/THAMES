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

#ifndef PARROTKILLOHMODELH
#define PARROTKILLOHMODELH

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <ctime>
#include "KineticModel.h"
// #include "KineticController.h"
#include "ChemicalSystem.h"
#include "Lattice.h"
#include "KineticData.h"
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
    
double wsRatio_;            /**< water-solid mass ratio */

double k1_;             /**< List of Parrot and Killoh <i>K</i><sub>1</sub> values */
double k2_;             /**< List of Parrot and Killoh <i>K</i><sub>2</sub> values */
double k3_;             /**< List of Parrot and Killoh <i>K</i><sub>3</sub> values */
double n1_;             /**< List of Parrot and Killoh <i>N</i><sub>1</sub> values */
double n3_;             /**< List of Parrot and Killoh <i>N</i><sub>3</sub> values */
double critDOH_;        /**< List of critical degrees of hydration for w/c
                                        effect in the Parrot and Killoh model */
double pfk_;            /**< Multiplicative factor for k's to account for

effects of pozzolanic additions */
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
@param lattice is a pointer to the Lattice object holding the microstructure
@param kineticData is the collection of kinetic parameters already stored
@param verbose is true if verbose output should be produced
@param warning is false if suppressing warning output
*/
ParrotKillohModel (ChemicalSystem *cs,
                   Solution *solut,
                   Lattice *lattice,
                   struct KineticData &kineticData,
                   const bool verbose,
                   const bool warning);
     

/**
@brief Get the type of kinetic model

@return a string indicating the model type
*/
string getType () const
{
    return (ParrotKillohType);
}

/**
@brief Set the w/s mass ratio of the system for the kinetic model equations.

@note NOT USED.

@param wsr is the w/s mass ratio to set
*/
void setWsRatio (double wsr)
{
    wsRatio_ = wsr;
}

/**
@brief Get the w/s mass ratio of the system used by the kinetic model equations.

@return the w/s mass ratio
*/
double getWsRatio () const
{
    return wsRatio_;
}

/**
@brief Get the list of <i>K</i><sub>1</sub> values for clinker phases in the PK model.

@note NOT USED.

@return the <i>K</i><sub>1</sub> value for the phase in the PK model
*/
double getK1 () const
{
    return k1_;
}

/**
@brief Get the list of <i>K</i><sub>2</sub> values for clinker phases in the PK model.

@note NOT USED.

@return the <i>K</i><sub>2</sub> value for the phase in the PK model
*/
double getK2 () const
{
    return k2_;
}

/**
@brief Get the list of <i>K</i><sub>3</sub> values for clinker phases in the PK model.

@note NOT USED.

@return the <i>K</i><sub>2</sub> value for the phase in the PK model
*/
double getK3 () const
{
    return k3_;
}

/**
@brief Set the multiplicative factor for k's due to pozzolanic effects

@param pfk is the multiplicative factor to set.
*/
void setPfk (const double pfk)
{
    pfk_ = pfk;
    if (pfk_ < 1.0e-5) pfk_ = 1.0e-5;
}

/**
@brief Get the multiplicative factor for k's due to pozzolanic effects

@note NOT USED.

@return the pfk_ value for the phase in the PK model
*/
double getPfk () const
{
    return pfk_;
}

/**
@brief Get the list of <i>N</i><sub>1</sub> values for clinker phases in the PK model.

@note NOT USED.

@return the <i>N</i><sub>2</sub> value for the phase in the PK model
*/
double getN1 () const
{
    return n1_;
}

/**
@brief Get the list of <i>N</i><sub>3</sub> values for clinker phases in the PK model.

@note NOT USED.

@return the <i>N</i><sub>3</sub> value for the phase in the PK model
*/
double getN3 () const
{
    return n3_;
}

/**
@brief Get the critical degrees of hydration for w/c effects in the kinetic model.

@note NOT USED.

@return the critical degree of hydration for this phase
*/
double getCritDOH () const
{
    return critDOH_;
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
@param rh is the internal relative humidity
@param dICMoles is the vector of moles of each IC changed by kinetics
@param dsolutICMoles is the vector of moles of each IC in solution changed by
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
                                   vector<double> &GEMPhaseMoles);

};      // End of ParrotKillohModel class

#endif
