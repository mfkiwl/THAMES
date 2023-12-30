/**
@file  KineticModel.cc
@brief Method definitions for the KineticModel class.

*/
#include "KineticModel.h"

KineticModel::KineticModel() {
  ///
  /// Clear out the vectors so they can be populated with values from the
  /// XML input file
  ///

  numPhases_ = 0;
  name_ = "";
  microPhaseId_ = 0;
  DCId_ = 0;
  GEMPhaseId_ = 0;
  activationEnergy_ = 0.0;
  scaledMass_ = 0.0;
  initScaledMass_ = 0.0;
  waterId_ = 1;
  ICNum_ = 0;
  DCNum_ = 0;
  ICName_.clear();
  DCName_.clear();
  GEMPhaseNum_ = 0;
  specificSurfaceArea_ = refSpecificSurfaceArea_ = 0.0;
  degreeOfReaction_ = lossOnIgnition_ = 0.0;

  ///
  /// The default is to not have sulfate attack or leaching, so we set the
  /// default time for initiating these simulations to an absurdly large value:
  /// 10 billion days or 27 million years
  ///

  sulfateAttackTime_ = 1.0e10;
  leachTime_ = 1.0e10;

  return;
}

void KineticModel::setKineticDCMoles() {

#ifdef DEBUG
  cout << "KineticModel::setKineMicDCmoles" << endl;
  cout.flush();
#endif

  try {
    int waterId = chemSys_->getDCId("H2O@");
    double waterMoles = chemSys_->getDCMoles(waterId);
    double waterMolarMass = chemSys_->getDCMolarMass(waterId);
    double waterMass = waterMoles * waterMolarMass;
#ifdef DEBUG
    cout << "KineticModel::setKineticDCmoles        "
         << chemSys_->getDCName(waterId) << ": Mass = " << waterMass
         << ", Molar mass = " << waterMolarMass << endl;
    cout.flush();
#endif
    if (chemSys_->getDCMolarMass(DCId_) <= 0.0) {
      throw FloatException("KineticModel", "setKineticDCmoles",
                           "Divide by zero error");
    }
#ifdef DEBUG
    cout << "KineticModel::setKineticDCmoles        Clinker phase " << name_
         << ": Mass = " << scaledMass_
         << ", Molar mass = " << chemSys_->getDCMolarMass(DCId_) << endl;
#endif
    chemSys_->setDCMoles(DCId_,
                         (scaledMass_ / chemSys_->getDCMolarMass(DCId_)));
  } catch (EOBException eex) {
    eex.printException();
    exit(1);
  } catch (FloatException fex) {
    fex.printException();
    exit(1);
  } catch (out_of_range &oor) {
    EOBException ex("KineticModel", "setKineticDCMoles", oor.what(), 0, 0);
    ex.printException();
    exit(1);
  }
  return;
}

void KineticModel::zeroKineticDCMoles() {
  try {
    chemSys_->setDCMoles(DCId_, 0.0);
  } catch (EOBException eex) {
    eex.printException();
    exit(0);
  }
  return;
}
