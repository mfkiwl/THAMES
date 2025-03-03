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

  double critporediam = lattice_->getLargestSaturatedPore(); // in nm
  critporediam *= 1.0e-9;                                    // in m
  rh_ = exp(-6.23527e-7 / critporediam / temperature_);
  rh_ = rh_ > 0.55 ? rh_ : 0.551;
  rhFactor_ = rh_;

  arrhenius_ = exp((activationEnergy_ / GASCONSTANT) *
                   ((1.0 / refT_) - (1.0 / temperature_)));

  ///
  /// The default is to not have sulfate attack or leaching, so we set the
  /// default time for initiating these simulations to an absurdly large value:
  /// 10 billion hours or 114,000 years
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

  refT_ = kineticData.reftemperature;

  modelName_ = "PozzolanicModel";
  name_ = kineticData.name;
  microPhaseId_ = kineticData.microPhaseId;
  DCId_ = kineticData.DCId;
  GEMPhaseId_ = kineticData.GEMPhaseId;
  activationEnergy_ = kineticData.activationEnergy;
  scaledMass_ = kineticData.scaledMass;
  initScaledMass_ = kineticData.scaledMass;

  temperature_ = lattice_->getTemperature();
  double critporediam = lattice_->getLargestSaturatedPore(); // in nm
  critporediam *= 1.0e-9;                                    // in m
  rh_ = exp(-6.23527e-7 / critporediam / temperature_);
  rh_ = rh_ > 0.55 ? rh_ : 0.551;
  rhFactor_ = rh_;

  arrhenius_ = exp((activationEnergy_ / GASCONSTANT) *
                   ((1.0 / refT_) - (1.0 / temperature_)));

  ///
  /// The default is to not have sulfate attack or leaching, so we set the
  /// default time for initiating these simulations to an absurdly large value:
  /// 10 billion hours or 114,000 years
  ///

  sulfateAttackTime_ = 1.0e10;
  leachTime_ = 1.0e10;

  return;
}

void PozzolanicModel::calculateKineticStep(const double timestep,
                                           double &scaledMass,
                                           double &massDissolved, int cyc,
                                           double totalDOR) {

  ///
  /// Initialize local variables
  ///

  // double arrhenius; // = 1.0;

  double dissrate = 1.0e9; // Nucleation and growth rate
  double diffrate = 1.0e9; // Diffusion rate

  double rate = 1.0e-10; // Selected rate

  // double massDissolved = 0.0;

  double DOR, newDOR;

  ///
  /// Determine if this is a normal step or a necessary
  /// tweak from a failed GEM_run call
  ///

  try {

    // Each component now has its own kinetic model and we
    // just want to know the *change* in IC moles caused by
    // this component's dissolution or growth.

    // if (timestep < leachTime_ && timestep < sulfateAttackTime_) {

    // @todo BULLARD PLACEHOLDER
    // Still need to implement constant gas phase composition
    // Will involve equilibrating gas with aqueous solution
    //
    // First step each iteration is to equilibrate gas phase
    // with the electrolyte, while forbidding anything new
    // from precipitating.

    // RH factor is the same for all clinker phases
    // double vfvoid = lattice_->getVolumeFraction(VOIDID);
    // double vfh2o = lattice_->getVolumeFraction(ELECTROLYTEID);

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

    if (initScaledMass_ > 0.0) {
      DOR = (initScaledMass_ - scaledMass_) / (initScaledMass_);
      // prevent DOR from prematurely stopping PK calculations
      DOR = min(DOR, 0.99);
    } else {
      throw FloatException("PozzolanicModel", "calculateKineticStep",
                           "initScaledMass_ = 0.0");
    }

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
      // double area = (specificSurfaceArea_ / 1000.0) * scaledMass_; // m2

      // JWB BEWARE: The new definition of area is truly a geometric calculation
      // made on the microstructure. It does not catch BET surface area
      // if that ends up being important.
      // This area value has units of m2 per 100 g of initial solid
      double area = lattice_->getSurfaceArea(microPhaseId_);

      // Saturation index , but be sure that there is only one GEM Phase
      /// @note Assumes there is only one phase in this microstructure
      /// component
      /// @todo Generalize to multiple phases in a component (how?)

      // double saturationIndex = solut_->getSI(GEMPhaseId_);
      double saturationIndex = chemSys_->getMicroPhaseSI(microPhaseId_);

      // activity of water
      double waterActivity = chemSys_->getDCActivity(chemSys_->getDCId("H2O@"));

      // This equation basically implements the Dove and Crerar rate
      // equation for quartz.  Needs to be calibrated for silica fume, but
      // hopefully the BET area and LOI will help do that.

      // baserateconst_ has units of mol/m2/h
      // area has units of m2 of phase per 100 g of total solid
      // Therefore dissrate has units of mol of phase per 100 g of all solid
      // per h
      if (saturationIndex < 1.0) {
        dissrate = baserateconst * rhFactor_ * pow(ohActivity, ohexp_) * area *
                   pow(waterActivity, 2.0) * (1.0 - (lossOnIgnition_ / 100.0)) *
                   (sio2_)*pow((1.0 - pow(saturationIndex, siexp_)), dfexp_);
      } else {
        dissrate = -baserateconst * rhFactor_ * pow(ohActivity, ohexp_) * area *
                   pow(waterActivity, 2.0) * (1.0 - (lossOnIgnition_ / 100.0)) *
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
          average_cdiff = (1.0 - pow(saturationIndex, (1.0 / dissolvedUnits_)));
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

      // dissrate has units of mol of phase per 100 g of all solid
      // per h
      /// @todo JWB Check to make sure that diffrate has same units as dissrate
      rate = dissrate;
      if (abs(diffrate) < abs(rate))
        rate = diffrate;
      int rate_ini = rate;
      rate *= (rhFactor_ * arrhenius_);

      // Mass dissolved has units of g of phase per 100 g of all initial solid
      massDissolved = rate * timestep * chemSys_->getDCMolarMass(DCId_); //

      if (verbose_) {
        cout
            << "    PozzolanicModel::calculateKineticStep rate/massDissolved : "
            << rate << " / " << massDissolved << endl;
      }

      scaledMass_ = max(scaledMass_ - massDissolved, 0.0); //

      newDOR = (initScaledMass_ - scaledMass_) / initScaledMass_; //

      scaledMass = scaledMass_;

      if (verbose_) {
        cout << "  ****************** PZM_hT = " << timestep
             << "\tcyc = " << cyc << "\tmicroPhaseId_ = " << microPhaseId_
             << "    microPhase = " << name_
             << "\tGEMPhaseIndex = " << GEMPhaseId_ << " ******************"
             << endl;
        cout << "   PZM_hT   " << "rhFacto_r: " << rhFactor_
             << "\tarrhenius_: " << arrhenius_
             << "\tsaturationIndex: " << saturationIndex
             << "\twaterActivity: " << waterActivity << endl;
        cout << "   PZM_hT   " << "dissrate: " << dissrate
             << "\tdiffrate: " << diffrate << "\trate_ini: " << rate_ini
             << "\trate: " << rate << endl;
        cout << "   PZM_hT   " << "DOR: " << DOR << "\tnewDOR: " << newDOR
             << "\tinitScaledMass_: " << initScaledMass_
             << "\tscaledMass_: " << scaledMass_
             << "\tmassDissolved: " << massDissolved << endl;
        cout.flush();
      }

    } else {
      throw DataException("PozzolanicModel", "calculateKineticStep",
                          "DOR >= 1.0");
    }

    //} // End of normal hydration block
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
    EOBException ex("PozzolanicModel", "calculateKineticStep", oor.what(), 0,
                    0);
    ex.printException();
    exit(1);
  }

  return;
}
