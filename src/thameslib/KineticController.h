/**
@file KineticController.h
@brief Declaration of the KineticController class.

@section Introduction
The `KineticController` class keeps track of all the
different kinetic models that govern the rate of hydration.

*/

#ifndef KINETICCONTROLLERH
#define KINETICCONTROLLERH

#include "ChemicalSystem.h"
#include "KineticData.h"
#include "KineticModel.h"
#include "Lattice.h"
#include "ParrotKillohModel.h"
#include "PozzolanicModel.h"
#include "StandardKineticModel.h"
#include "global.h"
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

/**
@class KineticController
@brief Manages kinetic models

THAMES allows some flexibility in defining different types of kinetic models.
*/

class KineticController {

private:
  int numPhases_; /**< Total number of phases in the kinetic model */
  ChemicalSystem *
      chemSys_; /**< Pointer to the ChemicalSystem object for this simulation */
  Solution
      *solut_; /**< Pointer to the aqueous phase object for this simulation */
  Lattice *
      lattice_; /**< Pointer to the lattice object holding the microstructure */
  vector<KineticModel *>
      phaseKineticModel_;    /***< Kinetic model for each phase */
  double temperature_;       /**< Temperature [K] */
  double refT_;              /**< Reference temperature for PK model [K] */
  double sulfateAttackTime_; /**< Time at which sulfate attack simulation starts
                                [days] */
  double leachTime_; /**< Time at which leaching simulation starts [days] */

  vector<string> name_; /**< List of names of phases in the kinetic model */
  vector<string> ICName_;
  vector<string> DCName_;
  vector<int> microPhaseId_; /**< List of microstructure ids that are in kinetic
                                model */
  vector<double> initScaledMass_;      /**< List of initial scaled masses */
  vector<double> scaledMass_;          /**< List of scaled masses */
  vector<double> specificSurfaceArea_; /**< List of specific surface areas */
  vector<double>
      refSpecificSurfaceArea_; /**< List of reference specific surface areas */
  vector<bool> isKinetic_;
  int waterId_;       /**< DC index for liquid water */
  int ICNum_;         /**< Number of ICs in chemical system */
  int DCNum_;         /**< Number of DCs in chemical system */
  int numPK_;         /**< Number of kinetically Parrott/Killoh components */
  int numPozzolanic_; /**< Number of kinetically pozzolanic components */
  int numStandard_;   /**< Number of kinetically standard components */
  int GEMPhaseNum_;   /**< Number of GEM phases in chemical system */
  bool verbose_;      /**< Flag for verbose output */
  bool warning_;      /**< Flag for warnining output */

public:
  /**
  @brief Default constructor.

  This constructor is not used in THAMES.  It just establishes default values
  for all the member variables.

  @note NOT USED.
  */
  KineticController();

  /**
  @brief Overloaded constructor.

  This constructor is used in THAMES.

  @param cs is a pointer to the ChemicalSystem object for the simulation
  @param solut is a pointer to the aqeuous solution object for the simulation
  @param lattic is a pointer to the Lattice object holding the microstructure
  @param fileName is the name of the XML file with the input for the kinetic
  model
  @param verbose is true if verbose output should be produced
  @param warning is false if suppressing warning output
  */
  KineticController(ChemicalSystem *cs, Solution *solut, Lattice *lattice,
                    const string &fileName, const bool verbose,
                    const bool warning);

  /**
  @brief Destructor does nothing.
  */
  virtual ~KineticController(){};

  /**
  @brief Initialize the kinetic data structure
  */
  void initKineticData(struct KineticData &kineticData) {
    kineticData.name = "";
    kineticData.microPhaseId = kineticData.GEMPhaseId = kineticData.DCId = 0;
    kineticData.type = "thermo";
    kineticData.scaledMass = 0.0;
    kineticData.specificSurfaceArea = kineticData.refSpecificSurfaceArea =
        385.0;
    kineticData.temperature = kineticData.reftemperature = 293.15;
    kineticData.k1 = kineticData.k2 = kineticData.k3 = 1.0;
    kineticData.n1 = kineticData.n3 = 1.0;
    kineticData.critDOR = 0.0;
    kineticData.dissolutionRateConst = 0.0;
    kineticData.diffusionRateConstEarly = 0.0;
    kineticData.diffusionRateConstLate = 0.0;
    kineticData.siexp = kineticData.dfexp = 0.0;
    kineticData.dorexp = kineticData.ohexp = 0.0;
    kineticData.dissolvedUnits = 1.0;
    kineticData.activationEnergy = 0.0;
    kineticData.loi = kineticData.sio2 = kineticData.al2o3 = kineticData.cao =
        0.0;
  }

  /**
  @brief Master method controlling the parsing of XML input to the kinetic
  model.

  @param docName is the name of the (purported) XML input file
  */
  void parseDoc(const string &docName);

  /**
  @brief Parse the input data for one phase in the XML input file.

  This method uses the libxml library, so this must be included.

  @param doc is a libxml pointer to the document head
  @param cur is a libxml pointer to the current node being parsed
  @param numEntry is the number of solid entries in the XML file, will be
  incremented
  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void parseMicroPhase(xmlDocPtr doc, xmlNodePtr cur, int &numEntry,
                       struct KineticData &kineticData);

  /**
  @brief Parse the kinetic data for one phase in the XML input file.

  This method uses the libxml library, so this must be included.

  @param doc is a libxml pointer to the document head
  @param cur is a libxml pointer to the current node being parsed
  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void parseKineticData(xmlDocPtr doc, xmlNodePtr cur,
                        struct KineticData &kineticData);

  /**
  @brief Parse the kinetic data for the Parrot-Killoh kinetic model.

  This method uses the libxml library, so this must be included.

  @param doc is a libxml pointer to the document head
  @param cur is a libxml pointer to the current node being parsed
  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void parseKineticDataForParrotKilloh(xmlDocPtr doc, xmlNodePtr cur,
                                       struct KineticData &kineticData);

  /**
  @brief Parse the kinetic data for the standard kinetic model.

  This method uses the libxml library, so this must be included.

  @param doc is a libxml pointer to the document head
  @param cur is a libxml pointer to the current node being parsed
  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void parseKineticDataForStandard(xmlDocPtr doc, xmlNodePtr cur,
                                   struct KineticData &kineticData);

  /**
  @brief Parse the kinetic data for the pozzolanic kinetic model.

  This method uses the libxml library, so this must be included.

  @param doc is a libxml pointer to the document head
  @param cur is a libxml pointer to the current node being parsed
  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void parseKineticDataForPozzolanic(xmlDocPtr doc, xmlNodePtr cur,
                                     struct KineticData &kineticData);

  /**
  @brief Get the scaled mass of the phase in the kinetic model.

  The scaled mass of a phase is its mass percent on a total solids basis.

  @return the vector of scaled masses [percent solids]
  */
  vector<double> getScaledMass() const { return scaledMass_; }

  /**
  @brief Get the scaled mass of one phase

  The scaled mass of a phase is its mass percent on a total solids basis.

  @note NOT USED.

  @param midx is the microstructure id of the phase to query
  @return the vector of scaled masses [percent solids]
  */

  double getScaledMass(const int midx) {
    try {
      return scaledMass_.at(midx);
    } catch (out_of_range &oor) {
      EOBException ex("KineticController", "getScaledMass", "scaledMass_",
                      scaledMass_.size(), midx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the <i>initial</i> mass of the microstructure phases

  The scaled mass of a phase is its mass percent on a total solids basis.

  @return the vector of initial scaled masses [percent solids]
  */
  vector<double> getInitScaledMass() const { return initScaledMass_; }

  /**
  @brief Get the <i>initial</i> scaled mass of one microstructure phase

  The scaled mass of a phase is its mass percent on a total solids basis.

  @note NOT USED.

  @param midx is the microstructure id of the phase to query
  @return the initial scaled mass [percent solids]
  */

  double getInitScaledMass(const int midx) {
    try {
      return initScaledMass_.at(midx);
    } catch (out_of_range &oor) {
      EOBException ex("KineticController", "getInitScaledMass",
                      "initScaledMass_", initScaledMass_.size(), midx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Compute normalized initial microstructure phase masses

  @note NOT USED

  Given the initial masses of all phases in the microstructure,
  this method scales them to 100 grams of solid.  In the process,
  this method also sets the initial moles of water in the
  chemical system definition.
  */

  void calcPhaseMasses(void);

  /**
  @brief Get sum of phase masses

  */

  double getSolidMass(void);

  /**
  @brief Get the list of specific surface areas

  @return the vector of specific surface areas [m2/kg]
  */
  vector<double> getSpecificSurfaceArea() const { return specificSurfaceArea_; }

  /**
  @brief Get the specific surface area of one phase


  @param midx is the microstructure id of the phase to query
  @return the specific surface area [m2/kg]
  */
  double getSpecificSurfaceArea(const int midx) {
    try {
      return specificSurfaceArea_.at(midx);
    } catch (out_of_range &oor) {
      EOBException ex("KineticController", "getSpecificSurfaceArea",
                      "specificSurfaceArea_", specificSurfaceArea_.size(),
                      midx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the list of reference specific surface areas

  @return the vector of reference specific surface areas [m2/kg]
  */
  vector<double> getRefSpecificSurfaceArea() const {
    return refSpecificSurfaceArea_;
  }

  /**
  @brief Get the reference specific surface area of one phase


  @param midx is the microstructure id of the phase to query
  @return the reference specific surface area [m2/kg]
  */
  double getRefSpecificSurfaceArea(const int midx) {
    try {
      return refSpecificSurfaceArea_.at(midx);
    } catch (out_of_range &oor) {
      EOBException ex("KineticController", "getRefSpecificSurfaceArea",
                      "refSpecificSurfaceArea_", refSpecificSurfaceArea_.size(),
                      midx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Make a kinetic model for a given phase

  This method uses the libxml library, so this must be included.

  @param doc is a libxml pointer to the document head
  @param cur is a libxml pointer to the current node being parsed
  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void makeModel(xmlDocPtr doc, xmlNodePtr cur,
                 struct KineticData &kineticData);

  /**
  @brief Get the ChemicalSystem object for the simulation used by the kinetic
  model.

  @note NOT USED.

  @return a pointer to the ChemicalSystem object
  */
  ChemicalSystem *getChemSys() const { return chemSys_; }

  /**
  @brief Set the simulation time at which to begin external sulfate attack.

  @param sattacktime is the simulation time to begin sulfate attack [days]
  */
  void setSulfateAttackTime(double sattacktime) {
    sulfateAttackTime_ = sattacktime;
  }

  /**
  @brief Get the simulation time at which to begin external sulfate attack.

  @note NOT USED.

  @return the simulation time to begin sulfate attack [days]
  */
  double getSulfateAttackTime(void) const { return sulfateAttackTime_; }

  /**
  @brief Set the simulation time at which to begin leaching.

  @param leachtime is the simulation time to begin leaching [days]
  */
  void setLeachTime(double leachtime) { leachTime_ = leachtime; }

  /**
  @brief Get the simulation time at which to begin leaching.

  @note NOT USED.

  @return the simulation time to begin leaching [days]
  */
  double getLeachTime(void) const { return leachTime_; }

  /**
  @brief Get the list of phase names used by the kinetic model.

  @note NOT USED.

  @return the vector of names of phases in the kinetic model
  */
  vector<string> getName() const { return name_; }

  /**
  @brief Get the name of phase with a given index in the kinetic model.

  @param i is the index of the phase in the kinetic model
  @return the name of the phase with index i
  */
  string getName(const unsigned int i) const {
    try {
      return name_.at(i);
    } catch (out_of_range &oor) {
      EOBException ex("KineticController", "getName", "name_", name_.size(), i);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Set kinetic model DC moles

  */
  void setKineticDCMoles() {
    for (int i = 0; i < phaseKineticModel_.size(); ++i) {
      phaseKineticModel_[i]->setKineticDCMoles();
    }
    return;
  }

  /**
  @brief Zero kinetic model DC moles

  */
  void zeroKineticDCMoles() {
    for (int i = 0; i < phaseKineticModel_.size(); ++i) {
      phaseKineticModel_[i]->zeroKineticDCMoles();
    }
    return;
  }

  /**
  @brief Set the effect of pozzolans on Parrot-Killoh kinetics

  @param timestep is the time interval to simulate [days]
  @param temperature is the absolute temperature during this step [K]
  @param isFirst is true if this is the first time step of the simulation, false
  otherwise
  */
  void setPozzEffectOnPK(void);

  /**
  @brief Master method for implementing one kinetic time step.

  In a given time step, a certain number of moles of each clinker phase will
  dissolve, and the instantly soluble phases will dissolve in the first time
  step.  This function determines the number of moles of each phase to dissolve,
  based on the time interval being simulated.  It then calculates the number of
  IC moles to promote to the thermodynamic system from those phases (which are
  outside the thermodynamic system because they are kinetically controlled),
  based on the stoichiometry.  Those IC moles are then added to the
  thermodynamic system, and the moles and mass of each kinetically controlled
  phase are changed accordingly.

  This is now a pure virtual function.

  @remark This method is very long and several parts are hard-coded when they
  should be made more general.

  @todo Split this method into more convenient chunks
  @todo Make the methods more general, less hardwiring of parameters
  @todo Make the local variable names more descriptive

  @param timestep is the time interval to simulate [days]
  @param temperature is the absolute temperature during this step [K]
  @param isFirst is true if this is the first time step of the simulation, false
  otherwise
  */
  void calculateDissolutionEvents(const double timestep,
                                  const double temperature, bool isFirst);

  /**
  @brief Set the verbose flag

  @param isverbose is true if verbose output should be produced
  */
  void setVerbose(const bool isverbose) {
    verbose_ = isverbose;
    return;
  }

  /**
  @brief Get the verbose flag

  @return the verbose flag
  */
  bool getVerbose() const { return verbose_; }

  /**
  @brief Set the warning flag

  @param iswarning is true if verbose output should be produced
  */
  void setWarning(const bool iswarning) {
    warning_ = iswarning;
    return;
  }

  /**
  @brief Get the warning flag

  @return the warning flag
  */
  bool getWarning() const { return warning_; }

  /**
  @brief Set the number of kinetically Parrot-Killoh phases

  @param numpk is the number of Parrot-Killot phases
  */
  void setNumPK(const int numpk) {
    if (numpk >= 0) {
      numPK_ = numpk;
    } else {
      numPK_ = 0;
    }
    return;
  }

  /**
  @brief Get the number of kinetically Parrot-Killoh phases

  @return the number of kinetically Parrot-Killoh phases
  */
  int getNumPK() const { return numPK_; }
  int numpk = 0;

  /**
  @brief Set the number of kinetically pozzolanic phases

  @param numpk is the number of kinetically pozzolanic phases
  */
  void setNumPozzolanic(const int numpozzolanic) {
    if (numpozzolanic >= 0) {
      numPozzolanic_ = numpozzolanic;
    } else {
      numPozzolanic_ = 0;
    }
    return;
  }

  /**
  @brief Get the number of kinetically pozzolanic phases

  @return the number of kinetically pozzolanic phases
  */
  int getNumPozzolanic() const { return numPozzolanic_; }

  /**
  @brief Set the number of kinetically standard phases

  @param numpk is the number of kinetically standard phases
  */
  void setNumStandard(const int numstandard) {
    if (numstandard >= 0) {
      numStandard_ = numstandard;
    } else {
      numStandard_ = 0;
    }
    return;
  }

  /**
  @brief Get the number of kinetically standard phases

  @return the number of kinetically standard phases
  */
  int getNumStandard() const { return numStandard_; }

}; // End of KineticController class

#endif
