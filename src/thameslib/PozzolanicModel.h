/**
@file PozzolanicModel.h
@brief Declaration of the PozzolanicModel class.

@section Introduction
This class implements a dissolution equation like Dove's for
amorphous silicates [1].

@section References

    -# PM Dove, N Han, AF Wallace, JJ De Yoreo, Kinetics of amorphous silica
dissolution and the paradox of the silica polymorphs, Proceedings of the
National Academy of Sciences USA, 105 (2008) 9903–9908.

*/

#ifndef POZZOLANICMODELH
#define POZZOLANICMODELH

#include "ChemicalSystem.h"
#include "KineticController.h"
#include "KineticData.h"
#include "KineticModel.h"
#include "Lattice.h"
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
@class PozzolanicModel
@brief Handles the kinetic model of pozzolanic materials

*/

class PozzolanicModel : public KineticModel {

protected:
  double dissolutionRateConst_;    /**< Rate constand for dissolution
                                      (mol/m2/s) */
  double diffusionRateConstEarly_; /**< Rate constant for early-age diffusion
                                      (mol/m2/s) */
  double diffusionRateConstLate_;  /**< Rate constant for later-age diffusion
                                      (mol/m2/s) */
  /**
  @brief Number of dissolved DC units per unit dissolution reaction
  */
  double dissolvedUnits_;
  double siexp_;  /**< Exponent on saturation index (unitless) */
  double dfexp_;  /**< Exponent on driving force (unitless) */
  double dorexp_; /**< Exponent on degree of reaction (unitless) */
  double ohexp_;  /**< Exponent on OH ion activity (unitless) */
  double sio2_;   /**< Mass fraction of SiO2 (unitless) */
  double al2o3_;  /**< Mass fraction of Al2O3 (unitless) */
  double cao_;    /**< Mass fraction of CaO (unitless) */

public:
  /**
  @brief Default constructor.

  This constructor is not used in THAMES.  It just establishes default values
  for all the member variables.

  @note NOT USED.
  */
  PozzolanicModel();

  /**
  @brief Overloaded constructor.

  This constructor is the one invoked by THAMES.  It can only be called once the
  various other objects for the simulation are allocated and constructed.

  @param cs is a pointer to the ChemicalSystem object for the simulation
  @param lattice is a pointer to the Lattice object holding the microstructure
  @param kineticData is the collection of kinetic parameters already stored
  @param verbose is true if verbose output should be produced
  @param warning is false if suppressing warning output
  */
  PozzolanicModel(ChemicalSystem *cs, Lattice *lattice,
                  struct KineticData &kineticData, const bool verbose,
                  const bool warning);

  /**
  @brief Get the type of kinetic model

  @return a string indicating the model type
  */
  string getType() const { return (PozzolanicType); }

  /**
  @brief Set the dissolution rate constant

  @note NOT USED.

  @param rc is the rate constant value to use
  */
  void setDissolutionRateConst(const double rc) {
    dissolutionRateConst_ = max(rc, 0.0);
  }

  /**
  @brief Get the dissolution rate constant

  @note NOT USED.

  @return the dissolution rate constant
  */
  double getDissolutionRateConst() const { return dissolutionRateConst_; }

  /**
  @brief Set the early-age diffusion rate constant

  @note NOT USED.

  @param rc is the early-age diffusion rate constant value to use
  */
  void setDiffusionRateConstEarly(const double rc) {
    diffusionRateConstEarly_ = max(rc, 0.0);
  }

  /**
  @brief Get the dissolution rate constant

  @note NOT USED.

  @return the early-age diffusion rate constant
  */
  double getDiffusionRateConstEarly() const { return diffusionRateConstEarly_; }

  /**
  @brief Set the later-age diffusion rate constant

  @note NOT USED.

  @param rc is the later-age diffusion rate constant value to use
  */
  void setDiffusionRateConstLate(const double rc) {
    diffusionRateConstLate_ = max(rc, 0.0);
  }

  /**
  @brief Get the later-age diffusion rate constant

  @note NOT USED.

  @return the later-age diffusion rate constant
  */
  double getDiffusionRateConstLate() const { return diffusionRateConstLate_; }

  /**
  @brief Set the exponent on the degree of reaction in the diffusion rate
  equation

  @note NOT USED.

  @param dorexp is the exponent value to use
  */
  void setDorexp(const double dorexp) { dorexp_ = max(dorexp, 0.0); }

  /**
  @brief Get the exponent on the degree of reaction in the diffusion rate
  equation

  @note NOT USED.

  @return the exponent on the degree of reaction
  */
  double getDorexp() const { return dorexp_; }

  /**
  @brief Set the number of dissolved DC units per unit dissolution

  @note NOT USED.

  @param dissolvedUnits is the value to set
  */
  void setDissolvedUnits(const double dissolvedUnits) {
    dissolvedUnits_ = max(dissolvedUnits, 1.0);
  }

  /**
  @brief Get the number of dissolved DC units per unit dissolution

  @note NOT USED.

  @return the number of dissolved DC units
  */
  double getDissolvedUnits() const { return dissolvedUnits_; }

  /**
  @brief Set the exponent on the saturation index

  @note NOT USED.

  @param siexp is the exponent value to use
  */
  void setSiexp(const double siexp) { siexp_ = max(siexp, 0.0); }

  /**
  @brief Get the exponent on the saturation index

  @note NOT USED.

  @return the exponent on the saturation index
  */
  double getSiexp() const { return siexp_; }

  /**
  @brief Set the exponent on the driving force

  @note NOT USED.

  @param dfexp is the exponent value to use
  */
  void setDfexp(const double dfexp) { dfexp_ = max(dfexp, 0.0); }

  /**
  @brief Get the exponent on the driving force

  @note NOT USED.

  @return the exponent on the driving force
  */
  double getDfexp() const { return dfexp_; }

  /**
  @brief Set the exponent on the hydroxyl ion activity

  @note NOT USED.

  @param ohexp is the exponent value to use
  */
  void setOhexp(const double ohexp) { ohexp_ = max(ohexp, 0.0); }

  /**
  @brief Get the exponent on the hydroxyl ion activity

  @note NOT USED.

  @return the exponent on the hydroxyl ion activity
  */
  double getOhexp() const { return ohexp_; }

  /**
  @brief Set the mass fraction of SiO2

  @note NOT USED.

  @param sio2 is the mass fraction of SiO2
  */
  void setSio2(const double sio2) {
    sio2_ = max(sio2, 0.0);
    if (sio2_ > 1.0)
      sio2_ = 1.0;
  }

  /**
  @brief Get the SiO2 mass fraction

  @note NOT USED.

  @return the SiO2 mass fraction
  */
  double getSio2() const { return sio2_; }

  /**
  @brief Set the mass fraction of Al2O3

  @note NOT USED.

  @param al2o3 is the mass fraction of Al2O3
  */
  void setAl2o3(const double al2o3) {
    al2o3_ = max(al2o3, 0.0);
    if (al2o3_ > 1.0)
      al2o3_ = 1.0;
  }

  /**
  @brief Get the Al2O3 mass fraction

  @note NOT USED.

  @return the Al2O3 mass fraction
  */
  double getAl2o3() const { return al2o3_; }

  /**
  @brief Set the mass fraction of CaO

  @note NOT USED.

  @param cao is the mass fraction of CaO
  */
  void setCao(const double cao) {
    cao_ = max(cao, 0.0);
    if (cao_ > 1.0)
      cao_ = 1.0;
  }

  /**
  @brief Get the CaO mass fraction

  @note NOT USED.

  @return the CaO mass fraction
  */
  double getCao() const { return cao_; }

  /**
  @brief Master method for implementing one kinetic time step.

  Overloaded from base class to handle pozzolanic materials

  @todo Split this method into more convenient chunks
  @todo Make the methods more general, less hardwiring of parameters
  @todo Make the local variable names more descriptive

  @param timestep is the time interval to simulate [days]
  @param temperature is the absolute temperature during this step [K]
  @param isFirst is true if this is the first time step of the simulation, false
  otherwise
  @param rh is the internal relative humidity
  @param dICMoles is the vector of moles of each IC changed by kinetics
  @param DCMoles is the vector of moles of each DC
  @param GEMPhaseMoles is the vector of moles of each phase in GEMS
  */
  virtual void calculateKineticEvent(const double timestep,
                                     const double temperature, bool isFirst,
                                     double rh, vector<double> &dICMoles,
                                     vector<double> &DCMoles,
                                     vector<double> &GEMPhaseMoles);

}; // End of PozzolanicModel class

#endif
