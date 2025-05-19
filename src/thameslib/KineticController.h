/**
@file KineticController.h
@brief Declaration of the KineticController class.

@section Introduction
The `KineticController` class keeps track of all the
different kinetic models that govern the rate of hydration.

*/

#ifndef KINETICCONTROLLERH
#define KINETICCONTROLLERH

#include "../Resources/include/nlohmann/json.hpp"
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
using json = nlohmann::json;

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
  Lattice *
      lattice_; /**< Pointer to the lattice object holding the microstructure */
  vector<KineticModel *>
      phaseKineticModel_;    /***< Kinetic model for each phase */
  double temperature_;       /**< Temperature [K] */
  double refT_;              /**< Reference temperature for PK model [K] */
  double sulfateAttackTime_; /**< Time at which sulfate attack simulation starts
                                [hours] */
  double leachTime_; /**< Time at which leaching simulation starts [hours] */

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
  // int waterId_;     /**< DC index for liquid water */
  int ICNum_;       /**< Number of ICs in chemical system */
  int DCNum_;       /**< Number of DCs in chemical system */
  int GEMPhaseNum_; /**< Number of GEM phases in chemical system */
  bool verbose_;    /**< Flag for verbose output */
  bool warning_;    /**< Flag for warnining output */

  vector<double> ICMoles_;
  vector<double> DCMoles_;
  vector<double> DCMolesIni_;
  vector<double> ICMolesTot_;
  vector<double>
      scaledMassIni_; /**< List of scaled masses before a given time step*/

  vector<int> impurityDCID_;
  vector<double> impurity_K2O_;
  vector<double> impurity_Na2O_;
  vector<double> impurity_Per_;
  vector<double> impurity_SO3_;

  int pKMsize_;

  double initScaledCementMass_;
  double hydTimeIni_;
  double hyd_time_;
  int waterDCId_; /**< coresp to DCName = "H2O@" */
  double beginAttackTime_;

  vector<double> surfaceAreaIni_;

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
  @param lattice is a pointer to the Lattice object holding the microstructure
  @param jsonFileName is the name of the JSON file with the input for the
  kinetic model
  @param verbose is true if verbose output should be produced
  @param warning is false if suppressing warning output
  */
  KineticController(ChemicalSystem *cs, Lattice *lattice,
                    const string &jsonFileName, const bool verbose,
                    const bool warning);

  /**
  @brief Destructor does nothing.
  */
  virtual ~KineticController();

  /**
  @brief Initialize the kinetic data structure
  */
  void initKineticData(struct KineticData &kineticData) {
    kineticData.name = "";
    kineticData.microPhaseId = kineticData.GEMPhaseId = kineticData.DCId = 0;
    kineticData.type = "thermo";
    kineticData.scaledMass = 0.0;
    kineticData.surfaceAreaMultiplier = 1.0;
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
  @brief Master method controlling the parsing of JSON input to the kinetic
  model.

  @param docName is the name of the (purported) JSON input file
  */
  void parseDoc(const string &docName);

  /**
  @brief Parse the input data for one phase in the JSON input file.

  @param cdi is an iterator over the JSON data
  @param numEntry is the number of solid entries in the JSON file, will be
  incremented
  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void parseMicroPhases(const json::iterator cdi, int &numEntry,
                        struct KineticData &kineticData);

  /**
  @brief Parse the kinetic data for one phase in the JSON input file.

  @todo Need error checking for what to do if a required entry is not present

  @param p is an iterator over the JSON data
  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void parseKineticData(const json::iterator p,
                        struct KineticData &kineticData);

  /**
  @brief Parse the kinetic data for the Parrot-Killoh kinetic model.

  @todo Need error checking for what to do if a required entry is not present

  @param pp is an iterator over the JSON data
  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void parseKineticDataForParrotKilloh(const json::iterator pp,
                                       struct KineticData &kineticData);

  /**
  @brief Parse the kinetic data for the standard kinetic model.

  @todo Need error checking for what to do if a required entry is not present

  @param pp is an iterator over the JSON data
  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void parseKineticDataForStandard(const json::iterator pp,
                                   struct KineticData &kineticData);

  /**
  @brief Parse the kinetic data for the pozzolanic kinetic model.

  @todo Need error checking for what to do if a required entry is not present

  @param pp is an iterator over the JSON data
  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void parseKineticDataForPozzolanic(const json::iterator pp,
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

  double getScaledMass(const int midx) { return scaledMass_[midx]; }

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
  double getInitScaledMass(const int midx) { return initScaledMass_[midx]; }

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
    return specificSurfaceArea_[midx];
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
    return refSpecificSurfaceArea_[midx];
  }

  /**
  @brief Make a kinetic model for a given phase

  @param kineticData is a reference to the KineticData structure for temporarily
  storing the input parameters.
  */
  void makeModel(struct KineticData &kineticData);

  /**
  @brief Get the ChemicalSystem object for the simulation used by the kinetic
  model.

  @note NOT USED.

  @return a pointer to the ChemicalSystem object
  */
  ChemicalSystem *getChemSys() const { return chemSys_; }

  /**
  @brief Set the simulation time at which to begin external sulfate attack.

  @param sattacktime is the simulation time to begin sulfate attack [hours]
  */
  void setSulfateAttackTime(double sattacktime) {
    sulfateAttackTime_ = sattacktime;
  }

  /**
  @brief Get the simulation time at which to begin external sulfate attack.

  @note NOT USED.

  @return the simulation time to begin sulfate attack [hours]
  */
  double getSulfateAttackTime(void) const { return sulfateAttackTime_; }

  /**
  @brief Set the simulation time at which to begin leaching.

  @param leachtime is the simulation time to begin leaching [hours]
  */
  void setLeachTime(double leachtime) { leachTime_ = leachtime; }

  /**
  @brief Get the simulation time at which to begin leaching.

  @note NOT USED.

  @return the simulation time to begin leaching [hours]
  */
  double getLeachTime(void) const { return leachTime_; }

  /**
  @brief Get the list of phase names used by the kinetic model.

  @note NOT USED.

  @return the vector of names of phases in the kinetic model
  */
  // vector<string> getName() const { return name_; }

  /**
  @brief Get the name of phase with a given index in the kinetic model.

  @param i is the index of the phase in the kinetic model
  @return the name of the phase with index i
  */
  // string getName(const unsigned int i) const { return name_[i]; }

  /**
  @brief Set kinetic model DC moles

  */
  // void setKineticDCMoles() {
  //   int size = phaseKineticModel_.size();
  //   for (int i = 0; i < size; ++i) {
  //     phaseKineticModel_[i]->setKineticDCMoles();
  //   }
  //   return;
  // }

  /**
  @brief Zero kinetic model DC moles

  */
  // void zeroKineticDCMoles() {
  //   int size = phaseKineticModel_.size();
  //   for (int i = 0; i < size; ++i) {
  //     phaseKineticModel_[i]->zeroKineticDCMoles();
  //   }
  //   return;
  // }

  /**
  @brief Set the effect of pozzolans on Parrot-Killoh kinetics

  @param timestep is the time interval to simulate [hours]
  @param temperature is the absolute temperature during this step [K]
  @param isFirst is true if this is the first time step of the simulation, false
  otherwise
  */
  void setPozzEffectOnPK(void);

  /**
  @brief Master method for implementing phase kinetics

  In a given time step, a certain number of moles of each phase may
  dissolve or precipitate.  This function determines the number of moles of each
  phase to change, based on the time interval being simulated.

  @remark This method is very long and several parts are hard-coded when they
  should be made more general.

  @todo Split this method into more convenient chunks
  @todo Make the methods more general, less hardwiring of parameters
  @todo Make the local variable names more descriptive

  @param timestep is the time interval to simulate [hours]
  @param temperature is the absolute temperature during this step [K]
  @param isFirst is true if this is the first time step of the simulation, false
  otherwise
  */
  void calculateKineticStep(double time, const double timestep, int cyc);

  void updateKineticStep(int cyc, int pId, double scaledMass, double timestep);

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

  vector<double> getICMoles(void) { return ICMoles_; }

  vector<double> getDCMoles(void) { return DCMoles_; }

  vector<bool> getIsKinetic(void) { return isKinetic_; }

  void setHydTimeIni(double val) { hydTimeIni_ = val; }

  void setIniAttackTime(const double val) { beginAttackTime_ = val; }

}; // End of KineticController class

#endif
