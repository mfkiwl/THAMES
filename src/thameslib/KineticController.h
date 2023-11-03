/**
@file KineticController.h
@brief Declaration of the KineticController class.

@section Introduction
The `KineticController` class keeps track of all the
different kinetic models that govern the rate of hydration.

*/

#ifndef KINETICCONTROLLERH
#define KINETICCONTROLLERH

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <ctime>
#include "ChemicalSystem.h"
#include "Lattice.h"
#include "KineticModel.h"
#include "ParrotKillohModel.h"
#include "PozzolanicModel.h"
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
    - `rateconst` is a generic rate constant for dissolution or growth (flux units)
    - `siexp` is the exponent on the saturation index in the rate equation
    - `dfexp` is the exponent on the driving force term in the rate equation
    - `ohexp` is the exponent on the hydroxyl ion activity in the rate equation
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
    double scaledMass;      /**< Mass percent on a total solids basis */
    double ssa;             /**< Specific surface area [m2/kg] */
    double refssa;          /**< Reference specific surface area [m2/kg] */
    double temperature;       /**< Temperature [K] */
    double reftemperature;    /**< Reference temperature [K] */
    double k1;                /**< Parrot and Killoh <i>K</i><sub>1</sub> parameter */
    double k2;                /**< Parrot and Killoh <i>K</i><sub>2</sub> parameter */
    double k3;                /**< Parrot and Killoh <i>K</i><sub>3</sub> parameter */
    double n1;                /**< Parrot and Killoh <i>N</i><sub>1</sub> parameter */
    double n3;                /**< Parrot and Killoh <i>N</i><sub>3</sub> parameter */
    double rateconst;         /**< Generic rate constant [mol/m2/s] */
    double siexp;             /**< Exponent on saturation index [unitless] */
    double dfexp;             /**< Exponent on driving force [unitless] */
    double ohexp;             /**< Exponent on OH ion activity [unitless] */
    double Ea;                /**< Apparent activation energy [J/mol/K] */
    double critDOH;           /**< Critical degree of hydration for w/c effect */
    vector<int> RdId;         /**< Vector of IC ids of the partitioned components in the phase */
    vector<double> RdVal;     /**< Vector of Rd values for each IC */
};
#endif

/**
@class KineticController
@brief Manages kinetic models

THAMES allows some flexibility in defining different types of kinetic models.
*/

class KineticController {

private:
    
int numPhases_;             /**< Total number of phases in the kinetic model */
ChemicalSystem *chemSys_;   /**< Pointer to the ChemicalSystem object for this simulation */
Solution *solut_;           /**< Pointer to the aqueous phase object for this simulation */
Lattice *lattice_;          /**< Pointer to the lattice object holding the microstructure */
vector<KineticModel * > phaseKineticModel_;  /***< Kinetic model for each phase */
double temperature_;        /**< Temperature [K] */
double refT_;               /**< Reference temperature for PK model [K] */
double sulfateAttackTime_;  /**< Time at which sulfate attack simulation starts [days] */
double leachTime_;          /**< Time at which leaching simulation starts [days] */

vector<string> name_;           /**< List of names of phases in the kinetic model */
vector<string> ICName_;
vector<string> DCName_;
vector<int> microPhaseId_;      /**< List of microstructure ids that are in kinetic model */
vector<bool> isKinetic_;
int waterId_;                   /**< DC index for liquid water */
int ICNum_;                     /**< Number of ICs in chemical system */
int DCNum_;                     /**< Number of DCs in chemical system */
int GEMPhaseNum_;               /**< Number of GEM phases in chemical system */
bool verbose_;                  /**< Flag for verbose output */
bool warning_;                  /**< Flag for warnining output */

public:
    
/**
@brief Default constructor.

This constructor is not used in THAMES.  It just establishes default values for
all the member variables.

@note NOT USED.
*/
KineticController ();

/**
@brief Overloaded constructor.

This constructor is used in THAMES.

@param cs is a pointer to the ChemicalSystem object for the simulation
@param solut is a pointer to the aqeuous solution object for the simulation
@param lattic is a pointer to the Lattice object holding the microstructure
@param fileName is the name of the XML file with the input for the kinetic model
@param verbose is true if verbose output should be produced
@param warning is false if suppressing warning output
*/
KineticController::KineticController (ChemicalSystem *cs,
                                      Solution *solut,
                                      Lattice *lattice,
                                      const string &fileName,
                                      const bool verbose,
                                      const bool warning);
     
/**
@brief Destructor does nothing.
*/
virtual ~KineticController() {};

/**
@brief Initialize the kinetic data structure
*/
void initKineticData()
{
    name = "";
    microPhaseId = GEMPhaseId = DCId = 0; 
    type = "thermo";
    scaledMass = 0.0;
    ssa = refssa = 385.0;
    temperature = reftemperature = 293.15;
    k1 = k2 = k3 = n1 = n3 = 1.0;
    rateconst = 0.0;
    siexp = dfexp = ohexp = 0.0;
    Ea = 0.0;
    critDOH = 0.0;
    RdId.clear();
    RdVal.clear();
}

/**
@brief Master method controlling the parsing of XML input to the kinetic model.

@param docName is the name of the (purported) XML input file
*/
void parseDoc (const string &docName);

/**
@brief Parse the input data for one phase in the XML input file.

This method uses the libxml library, so this must be included.

@param doc is a libxml pointer to the document head
@param cur is a libxml pointer to the current node being parsed
@param numEntry is the number of solid entries in the XML file, will be incremented
@param kineticData is a structure that can store all kinetic data for a phase
@param kineticData is a reference to the KineticData structure for temporarily storing
            the input parameters.
*/
void parsePhase (xmlDocPtr doc,
                 xmlNodePtr cur,
                 int &numEntry,
                 KineticData &kineticData);

/**
@brief Parse the kinetic data for one phase in the XML input file.

This method uses the libxml library, so this must be included.

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

This method uses the libxml library, so this must be included.

@param doc is a libxml pointer to the document head
@param cur is a libxml pointer to the current node being parsed
@param kineticData is a reference to the KineticData structure for temporarily storing
            the input parameters
*/
void parseRdData (xmlDocPtr doc,
                  xmlNodePtr cur,
                  KineticData &kineticData);


/**
@brief Make a kinetic model for a given phase

This method uses the libxml library, so this must be included.

@param doc is a libxml pointer to the document head
@param cur is a libxml pointer to the current node being parsed
@param numEntry is the number of solid entries in the XML file, will be incremented
@param kineticData is a reference to the KineticData structure for temporarily storing
            the input parameters.
*/
void makeModel (xmlDocPtr doc,
                xmlNodePtr cur,
                int &numEntry,
                KineticData &kineticData);

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
        EOBException ex("KineticController","getName","name_",
                           name_.size(),i);
        ex.printException();
        exit(1);
    }
}
/**
@brief Initialize a kinetic data structure 
*/
void initKineticData(struct KineticData &kd)
{
    kd.name = "";
    kd.microPhaseId = 0;
    kd.GEMPhaseId.clear();
    kd.DCId.clear();
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
    kd.ssa = 385.0;         
    kd.refssa = 385.0;       
    kd.temperature = 296.0;     
    kd.reftemperature = 296.0;  
    kd.rateconst = 0.0;       
    kd.siexp = 0.0;          
    kd.dfexp = 0.0;         
    kd.ohexp = 0.0;        

    return;
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
void calculateKineticStep (const double timestep,
                                   const double temperature,
                                   bool isFirst);
     
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

};      // End of KineticController class

#endif
