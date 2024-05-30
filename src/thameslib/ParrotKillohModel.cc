/**
@file  ParrotKillohModel.cc
@brief Method definitions for the ParrotKillohModel class.

*/
#include "ParrotKillohModel.h"

ParrotKillohModel::ParrotKillohModel() {
  ///
  /// Default value for w/c ratio in PK model is 0.45
  ///

  wsRatio_ = 0.45;

  ///
  /// Default value for specific surface area in PK model is 385
  /// m<sup>2</sup>/kg
  ///

  specificSurfaceArea_ = 385.0;
  refSpecificSurfaceArea_ = 385.0; // reference SSA (m2/kg)

  ///
  /// Default temperature in the PK model is 20 C (or 293 K)
  ///

  temperature_ = 293.15; // default temperature (K)
  refT_ = 293.15;        // default temperature (K)

  ///
  /// Clear out the vectors so they can be populated with values from the
  /// XML input file
  ///

  name_ = "";
  microPhaseId_ = 2;
  DCId_ = 2;
  GEMPhaseId_ = 2;
  k1_ = 1.0;
  k2_ = 1.0;
  k3_ = 1.0;
  pfk_ = 1.0;
  n1_ = 1.0;
  n3_ = 1.0;
  activationEnergy_ = 0.0;
  scaledMass_ = 0.0;
  initScaledMass_ = 0.0;
  initScaledMoles_ = 0.0;
  critDOR_ = 100.0;
  degreeOfReaction_ = 0.0;
  ICNum_ = 0;
  ICName_.clear();
  DCNum_ = 0;
  DCName_.clear();
  GEMPhaseNum_ = 0;

  ///
  /// The default is to not have sulfate attack or leaching, so we set the
  /// default time for initiating these simulations to an absurdly large value:
  /// 10 billion days or 27 million years
  ///

  sulfateAttackTime_ = 1.0e10;
  leachTime_ = 1.0e10;

  return;
}

ParrotKillohModel::ParrotKillohModel(ChemicalSystem *cs, Lattice *lattice,
                                     struct KineticData &kineticData,
                                     const bool verbose, const bool warning) {
  // Set the verbose and warning flags

  verbose_ = verbose;
  warning_ = warning;
#ifdef DEBUG
  verbose_ = true;
  warning_ = true;
  cout << "ParrotKillohModel::ParrotKillohModel Constructor" << endl;
  cout.flush();
#else
  verbose_ = verbose;
  warning_ = warning;
#endif

  chemSys_ = cs;
  lattice_ = lattice;

  ///
  /// Default value for specific surface area in PK model is 385
  /// m<sup>2</sup>/kg
  ///

  specificSurfaceArea_ = kineticData.specificSurfaceArea;
  refSpecificSurfaceArea_ = kineticData.refSpecificSurfaceArea;
  ssaFactor_ = specificSurfaceArea_ / refSpecificSurfaceArea_;

  ///
  /// Default value for w/c ratio in PK model is 0.45
  ///

  wsRatio_ = lattice_->getWsratio();

  ///
  /// Default initial solid mass is 100 g
  ///

  initSolidMass_ = 100.0;

  ///
  /// Default temperature in the PK model is 20 C (or 293 K)
  ///

  lossOnIgnition_ = kineticData.loi;
  temperature_ = kineticData.temperature;
  refT_ = kineticData.reftemperature;

  name_ = kineticData.name;
  microPhaseId_ = kineticData.microPhaseId;
  DCId_ = kineticData.DCId;
  GEMPhaseId_ = kineticData.GEMPhaseId;
  ;
  k1_ = kineticData.k1;
  k2_ = kineticData.k2;
  k3_ = kineticData.k3;

  /// Default to no pozzolanic influence on clinker phases

  pfk_ = 1.0;

  n1_ = kineticData.n1;
  n3_ = kineticData.n3;
  activationEnergy_ = kineticData.activationEnergy;
  scaledMass_ = kineticData.scaledMass;
  initScaledMass_ = kineticData.scaledMass;
  initScaledMoles_ = initScaledMass_ / (chemSys_->getDCMolarMass(DCId_));
  critDOR_ = kineticData.critDOR;
  degreeOfReaction_ = 0.0;

  waterId_ = chemSys_->getDCId("H2O@");
  ICNum_ = chemSys_->getNumICs();
  ICName_ = chemSys_->getICName();
  DCNum_ = chemSys_->getNumDCs();
  DCName_ = chemSys_->getDCName();
  GEMPhaseNum_ = chemSys_->getNumGEMPhases();

  /// The default is to not have sulfate attack or leaching, so we set the
  /// default time for initiating these simulations to an absurdly large value:
  /// 10 billion days or 27 million years
  ///

  sulfateAttackTime_ = 1.0e10;
  leachTime_ = 1.0e10;

  return;
}

void ParrotKillohModel::calculateKineticEvent(const double timestep,
                                              const double temperature,
                                              bool isFirst, double rh,
                                              vector<double> &dICMoles,
                                              vector<double> &DCMoles,
                                              vector<double> &GEMPhaseMoles) {
  ///
  /// Initialize local variables
  ///

  double T = temperature;
  double wsFactor = 1.0;
  double rhFactor = 1.0;
  double DOR = getDegreeOfReaction();
  double cDOR = 0.0;
  double arrhenius = 1.0;

  double ngrate = 1.0e-10;   // Nucleation and growth rate
  double hsrate = 1.0e-10;   // Hydration shell rate
  double diffrate = 1.0e-10; // Diffusion rate
  double rate = 1.0e-10;     // Selected rate

  double newDOR = DOR; // Updated value of doh
  double massDissolved = 0.0;
  double DCMolarMass = chemSys_->getDCMolarMass(DCId_);
  double scaledMoles, scaledMolesDissolved;

  // Is this the first time? If so, then for the ParrotKilloh model
  // we will add all the IC moles for this phase at the beginning

  if (isFirst) {
    chemSys_->setMicroPhaseMass(microPhaseId_, initScaledMass_);
    chemSys_->setMicroPhaseMassDissolved(microPhaseId_, 0.0);

    /// Convert mass and mass dissolved to moles and moles dissolved

    scaledMoles = initScaledMass_ / DCMolarMass;
    scaledMolesDissolved = massDissolved / DCMolarMass;

    cout << "dICMoles in ParrotKilloh Model:" << endl;
    for (int ii = 0; ii < dICMoles.size(); ++ii) {
      dICMoles[ii] +=
          ((initScaledMass_ / DCMolarMass) * chemSys_->getDCStoich(DCId_, ii));
    }
  }

  ///
  /// Determine if this is a normal step or a necessary
  /// tweak from a failed GEM_run call
  ///

  try {
    static int conc_index = 0;
    int microPhaseId, DCId, ICId;
    double molarMass;
    GEMPhaseMoles.resize(GEMPhaseNum_, 0.0);
    string icn;

    // @todo BULLARD PLACEHOLDER
    // Still need to implement constant gas phase composition
    // Will involve equilibrating gas with aqueous solution
    //
    // First step each iteration is to equilibrate gas phase
    // with the electrolyte, while forbidding anything new
    // from precipitating.

    // if (verbose_) {
    // cout << "ParrotKillohModel::calculateKineticEvent for " << name_
    // << endl;
    // cout.flush();
    // }

    vector<double> impurityRelease;
    impurityRelease.clear();
    impurityRelease.resize(chemSys_->getNumMicroImpurities(), 0.0);

    // RH factor is the same for all clinker phases
    // rh variable passed to this function pre-calculated by KineticController
    // class

    rh = rh > 0.55 ? rh : 0.551;
    rhFactor = pow(((rh - 0.55) / 0.45), 4.0);

    if (initScaledMass_ > 0.0) {
      DOR = (initScaledMass_ - scaledMass_) / (initScaledMass_);
      DOR = min(DOR, 0.99); // prevents DOR from prematurely
                            // stopping PK calculations
      // if (verbose_) {
      // cout << "~~~~>DOR for " << name_ << " = " << DOR << endl;
      // cout.flush();
      // }
    } else {
      throw FloatException("ParrotKillohModel", "calculateKineticEvent",
                           "initScaledMass_ = 0.0");
    }

    cDOR = 1.333 * wsRatio_;
    wsFactor = 1.0;
    if (DOR > cDOR) {
      wsFactor += ((4.444 * wsRatio_) - (3.333 * DOR));
      wsFactor = pow(wsFactor, 4.0);
    }

    /*
    wsFactor = 1.0 + (3.333 *
           (pow(((critDOR_[i] * wsRatio_) - DOR),4.0)));
    */

    arrhenius =
        exp((activationEnergy_ / GASCONSTANT) * ((1.0 / refT_) - (1.0 / T)));

    if (DOR < 1.0) {

      // Normal Parrott and Killoh implementation here

      if (fabs(n1_) > 0.0) {
        ngrate =
            (k1_ / n1_) * (1.0 - DOR) * pow((-log(1.0 - DOR)), (1.0 - n1_));
        ngrate *= (ssaFactor_); // only used for the N+G rate

        if (ngrate < 1.0e-10)
          ngrate = 1.0e-10;
      } else {
        throw FloatException("ParrotKillohModel", "calculateKineticEvent",
                             "n1_ = 0.0");
      }

      hsrate = k3_ * pow((1.0 - DOR), n3_);
      if (hsrate < 1.0e-10)
        hsrate = 1.0e-10;

      if (DOR > 0.0) {
        diffrate = (k2_ * pow((1.0 - DOR), (2.0 / 3.0))) /
                   (1.0 - pow((1.0 - DOR), (1.0 / 3.0)));
        if (diffrate < 1.0e-10)
          diffrate = 1.0e-10;
      } else {
        diffrate = 1.0e9;
      }

      rate = (ngrate < hsrate) ? ngrate : hsrate;
      if (diffrate < rate)
        rate = diffrate;
      rate *= (pfk_ * wsFactor * rhFactor * arrhenius);
      newDOR = DOR + (rate * timestep);

      setDegreeOfReaction(newDOR);

      // if (verbose_) {
      // cout << "PK model for " << name_ << ", ngrate = " << ngrate
      // << ", hsrate = " << hsrate << ", diffrate = " << diffrate
      // << ", rhFactor = " << rhFactor << ", wsFactor = " << wsFactor
      // << ", RATE = " << rate << ", timestep " << timestep
      // << ", oldDOR = " << DOR << ", new DOR = " << newDOR << endl;
      // cout.flush();
      // }

      /// @note This where we can figure out the volume dissolved
      /// and link it back to the current volume to see how many
      /// voxels need to dissolve
      ///

      scaledMass_ = initScaledMass_ * (1.0 - newDOR);
      massDissolved = (newDOR - DOR) * initScaledMass_;

      chemSys_->setMicroPhaseMass(microPhaseId_, scaledMass_);
      chemSys_->setMicroPhaseMassDissolved(microPhaseId_, massDissolved);

      /// Convert mass and mass dissolved to moles and moles dissolved

      scaledMoles = scaledMass_ / DCMolarMass;
      scaledMolesDissolved = massDissolved / DCMolarMass;

      chemSys_->setDCLowerLimit(DCId_, (scaledMoles - scaledMolesDissolved));

      // if (verbose_) {
      // cout << "ParrotKillohModel::calculateKineticEvent "
      // << "Original scaled mass = " << initScaledMass_
      // << ", dissolved scaled mass = " << massDissolved << endl;
      // cout << "New scaled mass = "
      // << chemSys_->getMicroPhaseMass(microPhaseId_)
      // << " and new volume = "
      // << chemSys_->getMicroPhaseVolume(microPhaseId_) << endl;
      // cout << "DC limits: [" << chemSys_->getDCLowerLimit(DCId_) << ","
      // << chemSys_->getDCUpperLimit(DCId_) << "]" << endl;
      // cout.flush();
      // cout << endl;
      // cout << "IC composition:" << endl;
      // for (int i = 0; i < chemSys_->getNumICs(); ++i) {
      // cout << chemSys_->getICName(i) << ": " << chemSys_->getICMoles(i)
      // << endl;
      // }
      // cout.flush();
      // }

      /// @note impurityRelease index values are assumed to
      /// be uniquely associated with particular chemical
      /// elements

      /// @todo Make this more general so that any indexing
      /// can be used

      /// @todo Allow any IC element to be an impurity, not
      /// necessarily just the ones hard-coded here

      impurityRelease[0] = (massDissolved * chemSys_->getK2o(microPhaseId_));
      impurityRelease[1] = (massDissolved * chemSys_->getNa2o(microPhaseId_));
      impurityRelease[2] = (massDissolved * chemSys_->getMgo(microPhaseId_));
      impurityRelease[3] = (massDissolved * chemSys_->getSo3(microPhaseId_));

      for (int ii = 0; ii < dICMoles.size(); ii++) {

        // Handling the actual clinker phase as a DC lower limits now,
        // because we set the IC moles of the phase at the very beginning
        // all at once.

        // Therefore, only need to do IC moles for impurities dissolved within
        // the clinker phase

        if (ICName_[ii] == "O") {
          // Dissolved K2O in this phase
          if (chemSys_->isIC("K")) {
            icn = "K";
            molarMass = 2.0 * chemSys_->getICMolarMass(icn);
            icn = "O";
            molarMass += chemSys_->getICMolarMass(icn);
            dICMoles[ii] += (impurityRelease[0] / molarMass);
          }
          // Dissolved Na2O in this phase
          if (chemSys_->isIC("Na")) {
            icn = "Na";
            molarMass = 2.0 * chemSys_->getICMolarMass(icn);
            icn = "O";
            molarMass += chemSys_->getICMolarMass(icn);
            dICMoles[ii] += (impurityRelease[1] / molarMass);
          }
          // Dissolved MgO in this phase
          if (chemSys_->isIC("Mg")) {
            icn = "Mg";
            molarMass = chemSys_->getICMolarMass(icn);
            icn = "O";
            molarMass += chemSys_->getICMolarMass(icn);
            dICMoles[ii] += (impurityRelease[2] / molarMass);
          }
          // Dissolved SO3 in this phase
          if (chemSys_->isIC("S")) {
            icn = "S";
            molarMass = chemSys_->getICMolarMass(icn);
            icn = "O";
            molarMass += (3.0 * chemSys_->getICMolarMass(icn));
            dICMoles[ii] += (3.0 * (impurityRelease[3] / molarMass));
          }
        } else if (ICName_[ii] == "S") {
          // Dissolved SO3  in this phase
          icn = "S";
          molarMass = chemSys_->getICMolarMass(icn);
          icn = "O";
          molarMass += (3.0 * chemSys_->getICMolarMass(icn));
          dICMoles[ii] += (impurityRelease[3] / molarMass);
        } else if (ICName_[ii] == "K") {
          // Dissolved K2O in this phase
          icn = "K";
          molarMass = 2.0 * chemSys_->getICMolarMass(icn);
          icn = "O";
          molarMass += chemSys_->getICMolarMass(icn);
          dICMoles[ii] += (2.0 * (impurityRelease[0] / molarMass));
        } else if (ICName_[ii] == "Na") {
          // Dissolved Na2O in this phase
          icn = "Na";
          molarMass = 2.0 * chemSys_->getICMolarMass(icn);
          icn = "O";
          molarMass += chemSys_->getICMolarMass(icn);
          dICMoles[ii] += (2.0 * (impurityRelease[1] / molarMass));
        } else if (ICName_[ii] == "Mg") {
          // Dissolved MgO in this phase
          icn = "Mg";
          molarMass = chemSys_->getICMolarMass(icn);
          icn = "O";
          molarMass += chemSys_->getICMolarMass(icn);
          dICMoles[ii] += (impurityRelease[2] / molarMass);
        }
      }

    } else {
      throw DataException("ParrotKillohModel", "calculateKineticEvent",
                          "DOR >= 1.0");
    }

  } // End of try block

  catch (EOBException eex) {
    eex.printException();
    exit(1);
  } catch (DataException dex) {
    dex.printException();
    exit(1);
  } catch (FloatException fex) {
    fex.printException();
    exit(1);
  } catch (out_of_range &oor) {
    EOBException ex("ParrotKillohModel", "calculateKineticEvent", oor.what(), 0,
                    0);
    ex.printException();
    exit(1);
  }

  return;
}
