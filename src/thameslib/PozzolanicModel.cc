/**
@file  PozzolanicModel.cc
@brief Method definitions for the PozzolanicModel class.

*/
#include "PozzolanicModel.h"

PozzolanicModel::PozzolanicModel() {

  ///
  /// Default value for specific surface area is 385 m<sup>2</sup>/kg
  ///

  specificSurfaceArea_ = 385.0;
  refSpecificSurfaceArea_ = 385.0; // reference specific surface area (m2/kg)
  ssaFactor_ = 1.0;

  ///
  /// Default temperature in the PK model is 20 C (or 293 K)
  ///

  temperature_ = 293.15; // default temperature (K)
  refT_ = 293.15;        // default temperature (K)

  ///
  /// Default values for the rate constants
  ///

  dissolutionRateConst_ = 0.0;
  diffusionRateConstEarly_ = 0.0;
  diffusionRateConstLate_ = 0.0;

  ///
  /// Default value for the exponents in the rate equation
  ///

  siexp_ = 1.0;
  dfexp_ = 1.0;
  ohexp_ = 0.0;
  sio2_ = 1.0;
  al2o3_ = cao_ = 0.0;
  lossOnIgnition_ = 0.0;

  name_ = "";
  microPhaseId_ = 2;
  DCId_ = 2;
  GEMPhaseId_ = 2;
  activationEnergy_ = 0.0;
  scaledMass_ = 0.0;
  initScaledMass_ = 0.0;

  ///
  /// The default is to not have sulfate attack or leaching, so we set the
  /// default time for initiating these simulations to an absurdly large value:
  /// 10 billion days or 27 million years
  ///

  sulfateAttackTime_ = 1.0e10;
  leachTime_ = 1.0e10;

  return;
}

PozzolanicModel::PozzolanicModel(ChemicalSystem *cs, Lattice *lattice,
                                 struct KineticData &kineticData,
                                 const bool verbose, const bool warning) {

  // Set the verbose and warning flags

  verbose_ = verbose;
  warning_ = warning;
#ifdef DEBUG
  verbose_ = true;
  warning_ = true;
  cout << "PozzolanicModel::PozzolanicModel Constructor sio2 value = "
       << kineticData.sio2 << endl;
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
  setSio2(kineticData.sio2);
  setAl2o3(kineticData.al2o3);
  setCao(kineticData.cao);
  setDissolutionRateConst(kineticData.dissolutionRateConst);
  setDiffusionRateConstEarly(kineticData.diffusionRateConstEarly);
  setDiffusionRateConstLate(kineticData.diffusionRateConstLate);
  setDissolvedUnits(kineticData.dissolvedUnits);
  setSiexp(kineticData.siexp);
  setDfexp(kineticData.dfexp);
  setOhexp(kineticData.ohexp);
  lossOnIgnition_ = kineticData.loi;

  ///
  /// Default initial solid mass is 100 g
  ///

  initSolidMass_ = 100.0;

  temperature_ = kineticData.temperature;
  refT_ = kineticData.reftemperature;

  name_ = kineticData.name;
  microPhaseId_ = kineticData.microPhaseId;
  DCId_ = kineticData.DCId;
  GEMPhaseId_ = kineticData.GEMPhaseId;
  activationEnergy_ = kineticData.activationEnergy;
  scaledMass_ = kineticData.scaledMass;
  initScaledMass_ = kineticData.scaledMass;

  ///
  /// The default is to not have sulfate attack or leaching, so we set the
  /// default time for initiating these simulations to an absurdly large value:
  /// 10 billion days or 27 million years
  ///

  sulfateAttackTime_ = 1.0e10;
  leachTime_ = 1.0e10;

  return;
}

void PozzolanicModel::calculateKineticEvent(const double timestep,
                                            const double temperature,
                                            bool isFirst, double rh,
                                            vector<double> &dICMoles,
                                            vector<double> &DCMoles,
                                            vector<double> &GEMPhaseMoles) {
  ///
  /// Initialize local variables
  ///

  double T = temperature;
  double arrhenius = 1.0;

  double dissrate = 1.0e9; // Nucleation and growth rate
  double diffrate = 1.0e9; // Diffusion rate

  double rate = 1.0e-10; // Selected rate

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

    for (int ii = 0; ii < dICMoles.size(); ++ii) {
      dICMoles[ii] +=
          ((initScaledMass_ / DCMolarMass) * chemSys_->getDCStoich(DCId_, ii));
    }
  }

  ///
  /// Determine if this is a normal step or a necessary
  /// tweak from a failed GEM_run call
  ///

  static double hyd_time = 0.0;
  hyd_time = hyd_time + timestep;

  try {
    static int conc_index = 0;
    int ICNum = chemSys_->getNumICs();
    int DCNum = chemSys_->getNumDCs();
    int numGEMPhases = chemSys_->getNumGEMPhases();
    int DCId, ICId;
    double molarMass;
    GEMPhaseMoles.clear();
    GEMPhaseMoles.resize(numGEMPhases, 0.0);
    string icn;

    vector<string> ICName;
    ICName.clear();
    ICName.resize(ICNum, " ");

    int waterId = chemSys_->getDCId("H2O@");

    // Each component now has its own kinetic model and we
    // just want to know the *change* in IC moles caused by
    // this component's dissolution or growth.

    for (int i = 0; i < ICNum; i++) {
      ICName[i] = chemSys_->getICName(i);
    }

    for (int i = 0; i < numGEMPhases; i++) {
      GEMPhaseMoles[i] = chemSys_->getGEMPhaseMoles(i);
    }

    if (hyd_time < leachTime_ && hyd_time < sulfateAttackTime_) {

      // @todo BULLARD PLACEHOLDER
      // Still need to implement constant gas phase composition
      // Will involve equilibrating gas with aqueous solution
      //
      // First step each iteration is to equilibrate gas phase
      // with the electrolyte, while forbidding anything new
      // from precipitating.

      vector<double> impurityRelease;
      impurityRelease.clear();
      impurityRelease.resize(chemSys_->getNumMicroImpurities(), 0.0);

      // RH factor is the same for all clinker phases
      double vfvoid = lattice_->getVolumefraction(VOIDID);
      double vfh2o = lattice_->getVolumefraction(ELECTROLYTEID);

      /// This is a big kluge for internal relative humidity
      /// @note Using new gel and interhydrate pore size distribution model
      ///       which is currently contained in the Lattice object.
      ///
      /// Surface tension of water is gamma = 0.072 J/m2
      /// Molar volume of water is Vm = 1.8e-5 m3/mole
      /// The Kelvin equation is
      ///    p/p0 = exp (-4 gamma Vm / d R T) = exp (-6.23527e-7 / (d T))
      ///
      ///    where d is the pore diameter in meters and T is absolute
      ///    temperature

      /// Assume a zero contact angle for now.
      /// @todo revisit the contact angle issue

      double critporediam = lattice_->getLargestSaturatedPore(); // in nm
      critporediam *= 1.0e-9;                                    // in m
      double rh = exp(-6.23527e-7 / critporediam / T);

      /// Assume a zero contact angle for now.
      /// @todo revisit the contact angle issue

      rh = rh > 0.55 ? rh : 0.551;
      double rhFactor = rh;
      // rhFactor = pow(((rh - 0.55)/0.45),4.0);

      double DOR = 0.0;
      double newDOR = 0.0;

      if (initScaledMass_ > 0.0) {
        DOR = (initScaledMass_ - scaledMass_) / (initScaledMass_);
        // prevent DOR from prematurely stopping PK calculations
        DOR = min(DOR, 0.99);
      } else {
        throw FloatException("PozzolanicModel", "calculateKineticEvent",
                             "initScaledMass_ = 0.0");
      }

      arrhenius =
          exp((activationEnergy_ / GASCONSTANT) * ((1.0 / refT_) - (1.0 / T)));

      if (DOR < 1.0) {

        /// BULLARD placeholder
        /// @note playing with different base rate constants here

        // double baserateconst = 5.6e-3 * arrhenius;  // mol m-2 s-1
        //

        double baserateconst = dissolutionRateConst_;

        /// @note The following influence alkali and alkali earth cations
        /// was asserted by Dove and Crerar (1990) but only at near-neutral pH

        double ca = chemSys_->getDCConcentration("Ca+2");
        double kca = 4.0e-7; // mol m-2 s-1 ads.
                             // rate const for Ca (guess)
        double Kca = 10.0;   // adsorption equilibrium
                             // constant is a guess
        double na = chemSys_->getDCConcentration("Na+");
        double kna = 6.35e-7; // mol m-2 s-1 ads. rate
                              // const from Dove and Crerar
        double Kna = 58.3;    // adsorption equilibrium
                              // constant from Dove and Crerar
        double k = chemSys_->getDCConcentration("K+");
        double kk = 5.6e-7; // mol m-2 s-1 ads. rate
                            // const from Dove and Crerar
        double Kk = 46.6;   // adsorption equilibrium constant
                            // from Dove and Crerar

        // Langmuir adsorption isotherms assumed to be additive

        baserateconst += (kca * Kca * ca / (1.0 + (Kca * ca)));
        baserateconst += (kna * Kna * na / (1.0 + (Kna * na)));
        baserateconst += (kk * Kk * k / (1.0 + (Kk * k)));

        double ohActivity = chemSys_->getDCActivity("OH-");
        double area = (specificSurfaceArea_ / 1000.0) * scaledMass_; // m2

        // Saturation index

        double saturationIndex = chemSys_->getMicroPhaseSI(microPhaseId_);

        // activity of water
        double waterActivity =
            chemSys_->getDCActivity(chemSys_->getDCId("H2O@"));

        // This equation basically implements the Dove and Crerar rate
        // equation for quartz.  Needs to be calibrated for silica fume, but
        // hopefully the BET area and LOI will help do that.

        if (saturationIndex < 1.0) {
          dissrate = baserateconst * rhFactor * pow(ohActivity, ohexp_) * area *
                     pow(waterActivity, 2.0) *
                     (1.0 - (lossOnIgnition_ / 100.0)) *
                     (sio2_)*pow((1.0 - pow(saturationIndex, siexp_)), dfexp_);
        } else {
          dissrate = -baserateconst * rhFactor * pow(ohActivity, ohexp_) *
                     area * pow(waterActivity, 2.0) *
                     (1.0 - (lossOnIgnition_ / 100.0)) *
                     (sio2_)*pow((pow(saturationIndex, siexp_) - 1.0), dfexp_);
        }

        /// Assume steady-state diffusion, with the surface being
        /// at equilibrium and the bulk being at the current
        /// saturation index.
        ///
        /// Also assume a particular, fixed boundary layer thickness
        /// through which diffusion occurs, like one micrometer

        double boundaryLayer = 1.0;

        double average_cdiff = 1.0e9;
        if (DOR > 0.0) {
          diffrate = (diffusionRateConstEarly_ * ssaFactor_ * (1.0 - DOR));
          /// Below is very rough approximation to chemical potential gradient
          /// Would be better if we knew the equilibrium constant of
          /// the dissociation reaction.  We would need to raise
          /// it to the power 1/dissolvedUnits and then multiply
          /// it by average_cdiff.
          if (saturationIndex < 1.0) {
            average_cdiff =
                (1.0 - pow(saturationIndex, (1.0 / dissolvedUnits_)));
            diffrate *= (average_cdiff) / boundaryLayer;
            if (abs(diffrate) < 1.0e-10)
              diffrate = 1.0e-10;
          } else {
            average_cdiff =
                -(pow(saturationIndex, (1.0 / dissolvedUnits_)) - 1.0);
            diffrate *= (average_cdiff) / boundaryLayer;
            if (abs(diffrate) < 1.0e-10)
              diffrate = -1.0e-10;
          }
        } else if (saturationIndex < 1.0) {
          diffrate = 1.0e9;
        } else {
          diffrate = -1.0e9;
        }

        cout << "PozzolanicModel::calculateKineticEvent for " << name_ << endl;
        cout << "  dissrate = " << dissrate << endl;
        cout << "    (baserateconst = " << baserateconst << ")" << endl;
        cout << "    (ssaFactor = " << ssaFactor_ << ")" << endl;
        cout << "    (DOR = " << DOR << ")" << endl;
        cout << "    (rhFactor = " << rhFactor << ")" << endl;
        cout << "    (arrhenius = " << arrhenius << ")" << endl;
        cout << "    (ohactivity = " << ohActivity << ")" << endl;
        cout << "    (waterterm = " << pow(waterActivity, 2.0) << ")" << endl;
        cout << "    (area = " << area << ")" << endl;
        cout << "    (saturationIndex = " << saturationIndex << ")" << endl;
        cout << "    (LOI = " << lossOnIgnition_ << ")" << endl;
        cout << "    (sio2 = " << sio2_ << ")" << endl;
        cout << "  diffrate = " << diffrate << endl;
        cout << "    (diffusionRateConstEarly = " << diffusionRateConstEarly_
             << ")" << endl;
        cout << "    (dissolvedUnits = " << dissolvedUnits_ << ")" << endl;
        cout << "    (average_cdiff = " << average_cdiff << ")" << endl;
        cout.flush();

        rate = dissrate;
        if (abs(diffrate) < abs(rate))
          rate = diffrate;
        rate *= (rhFactor * arrhenius);
        newDOR = DOR + (rate * timestep);

        cout << "  Final rate = " << rate << endl;
        cout.flush();

        scaledMass_ = max(initScaledMass_ * (1.0 - newDOR), 0.0);
        massDissolved = max((newDOR - DOR) * initScaledMass_, 0.0);

        chemSys_->setMicroPhaseMass(microPhaseId_, scaledMass_);
        chemSys_->setMicroPhaseMassDissolved(microPhaseId_, massDissolved);

        /// Convert mass and mass dissolved to moles and moles dissolved

        scaledMoles = scaledMass_ / DCMolarMass;
        scaledMolesDissolved = massDissolved / DCMolarMass;

        chemSys_->setDCLowerLimit(DCId_, (scaledMoles - scaledMolesDissolved));

        if (verbose_) {
          cout << "PozzolanicModel::calculateKineticEvent Original scaled "
                  "mass = "
               << initScaledMass_ << " and new scaled mass = "
               << chemSys_->getMicroPhaseMass(microPhaseId_)
               << " and new volume = "
               << chemSys_->getMicroPhaseVolume(microPhaseId_) << endl;
          cout.flush();
        }

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

          if (ICName[ii] == "O") {
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
          } else if (ICName[ii] == "S") {
            // Dissolved SO3  in this phase
            icn = "S";
            molarMass = chemSys_->getICMolarMass(icn);
            icn = "O";
            molarMass += (3.0 * chemSys_->getICMolarMass(icn));
            dICMoles[ii] += (impurityRelease[3] / molarMass);
          } else if (ICName[ii] == "K") {
            // Dissolved K2O in this phase
            icn = "K";
            molarMass = 2.0 * chemSys_->getICMolarMass(icn);
            icn = "O";
            molarMass += chemSys_->getICMolarMass(icn);
            dICMoles[ii] += (2.0 * (impurityRelease[0] / molarMass));
          } else if (ICName[ii] == "Na") {
            // Dissolved Na2O in this phase
            icn = "Na";
            molarMass = 2.0 * chemSys_->getICMolarMass(icn);
            icn = "O";
            molarMass += chemSys_->getICMolarMass(icn);
            dICMoles[ii] += (2.0 * (impurityRelease[1] / molarMass));
          } else if (ICName[ii] == "Mg") {
            // Dissolved MgO in this phase
            icn = "Mg";
            molarMass = chemSys_->getICMolarMass(icn);
            icn = "O";
            molarMass += chemSys_->getICMolarMass(icn);
            dICMoles[ii] += (impurityRelease[2] / molarMass);
          }
        }

      } else {
        throw DataException("PozzolanicModel", "calculateKineticEvent",
                            "DOR >= 1.0");
      }

    } // End of normal hydration block
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
    EOBException ex("PozzolanicModel", "calculateKineticEvent", oor.what(), 0,
                    0);
    ex.printException();
    exit(1);
  }

  return;
}
