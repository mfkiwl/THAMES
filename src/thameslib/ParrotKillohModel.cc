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
  critDOH_ = 100.0;
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
  wcRatio_ = lattice_->getWcratio();

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

  k1_ = kineticData.k1;
  k2_ = kineticData.k2;
  k3_ = kineticData.k3;

  /// Default to no pozzolanic influence on clinker phases

  pfk_ = 1.0;

  n1_ = kineticData.n1;
  n3_ = kineticData.n3;

  HLK_ = kineticData.HLK;

  activationEnergy_ = kineticData.activationEnergy;
  scaledMass_ = kineticData.scaledMass;
  initScaledMass_ = kineticData.scaledMass;
  critDOH_ = kineticData.critDOH;
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

//void ParrotKillohModel::calculateKineticStep(...) initial
/*
void ParrotKillohModel::calculateKineticStep (const double timestep, const double temperature,
                                             double rh, double &scaledMass,
                                             double &massDissolved, int cyc, double totalDOR) {

    ///
    /// Initialize local variables
    ///

    double T = temperature;
    double DOR, newDOR;

    double wsFactor = 1.0;
    double rhFactor = 1.0;
    double cDOR = 0.0;
    double arrhenius = 1.0;
    double ngrate = 1.0e-10;   // Nucleation and growth rate
    double hsrate = 1.0e-10;   // Hydration shell rate
    double diffrate = 1.0e-10; // Diffusion rate
    double rate = 1.0e-10;     // Selected rate

    ///
    /// Determine if this is a normal step or a necessary
    /// tweak from a failed GEM_run call
    ///

    try {

        // @todo BULLARD PLACEHOLDER
        // Still need to implement constant gas phase composition
        // Will involve equilibrating gas with aqueous solution
        //
        // First step each iteration is to equilibrate gas phase
        // with the electrolyte, while forbidding anything new
        // from precipitating.

        if (verbose_) {
            cout << "ParrotKillohModel::calculateKineticStep for " << name_ << endl;
            cout.flush();
        }

        // RH factor is the same for all clinker phases

        /// This is a big kluge for internal relative humidity
        /// @note Using new gel and interhydrate pore size distribution model
        ///       which is currently contained in the Lattice object.
        ///
        /// Surface tension of water is gamma = 0.072 J/m2
        /// Molar volume of water is Vm = 1.8e-5 m3/mole
        /// The Kelvin equation is
        ///    p/p0 = exp (-4 gamma Vm / d R T) = exp (-6.23527e-7 / (d T))
        ///
        ///    where d is the pore diameter in meters and T is absolute temperature

        /// Assume a zero contact angle for now.
        /// @todo revisit the contact angle issue

        //rh = rh > 0.55 ? rh : 0.551;
        rhFactor = pow(((rh - 0.55) / 0.45), 4.0);
        if (initScaledMass_ > 0.0) {
            DOR = (initScaledMass_ - scaledMass_) / initScaledMass_;
            DOR = min(DOR, 0.99); // prevents DOR from prematurely stopping PK calculations

            if (verbose_) {
                cout << "~~~~>DOR for " << name_ << " = " << DOR << endl;
                cout.flush();
            }
        } else {
            throw FloatException("ParrotKillohModel", "calculateKineticStep",
                                 "initScaledMass_ = 0.0");
        }
        cDOR = 1.333 * wsRatio_;
        wsFactor = 1.0;
        if (DOR > cDOR) {
            wsFactor += ((4.444 * wsRatio_) - (3.333 * DOR));
            wsFactor = pow(wsFactor, 4.0);
        }

        //    wsFactor = 1.0 + (3.333 *
        //      (pow(((critDOH_[i] * wsRatio_) - DOH),4.0)));


        arrhenius =
            exp((activationEnergy_ / GASCONSTANT) * ((1.0 / refT_) - (1.0 / T)));

        //cout << "PKM model for " << name_ << "    k1 = " << k1_ << "    k2 = " << k2_ << "    k3 = " << k3_
        //     << "    pfk = " << pfk_ << "    n1 = " << n1_ << "    n3 = " << n3_
        //     << "    Ea = " << activationEnergy_ << "    ssaFactor_ = " << ssaFactor_ << "    DOR = " << DOR << endl;
        //cout.flush();

        if (DOR < 1.0) {
            // Normal Parrott and Killoh implementation here

            if (fabs(n1_) > 0.0) {
                ngrate =
                    (k1_ / n1_) * (1.0 - DOR) * pow((-log(1.0 - DOR)), (1.0 - n1_));
                ngrate *= (ssaFactor_); // only used for the N+G rate

                if (ngrate < 1.0e-10)
                    ngrate = 1.0e-10;
            } else {
                throw FloatException("ParrotKillohModel", "calculateKineticStep",
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

            //cout << "ParrotKillohModel::calculateKineticStep for " << name_ << "\tDCId_: " << DCId_
            //     << "\tmicroPhaseId_: " << microPhaseId_ << endl;
            //cout << "  dissrate = " << rate << endl;
            //cout << "    (DOR = " << DOR << ")" << endl;
            //cout << "    (rhFactor = " << rhFactor << ")" << endl;
            //cout << "    (arrhenius = " << arrhenius << ")" << endl;
            //cout << "    (LOI = " << lossOnIgnition_ << ")" << endl;
            //cout.flush();

            double rate_ini = rate;

            rate *= (pfk_ * wsFactor * rhFactor * arrhenius);

            /// @note This where we can figure out the volume dissolved
            /// and link it back to the current volume to see how many
            /// voxels need to dissolve
            ///

            /// @note This all depends on concept of degree of
            /// hydration as defined by the PK model

            /// @todo Make this independent of PK model

            //double prod = rate * timestep;
            //newDOR = DOR + prod;
            scaledMass_ = initScaledMass_ * (1.0 - newDOR);
            //massDissolved = (newDOR - DOR) * initScaledMass_;
            massDissolved = initScaledMass_ * prod;
            scaledMass = scaledMass_;

            cout << "****************** PKM_hT = " << timestep << "    cyc = " << cyc
                 << "    microPhaseId_ = " << microPhaseId_ << "    microPhase = " << name_
                 << "    GEMPhaseIndex = " << GEMPhaseId_ << " ******************" << endl;
            cout << "PKM_hT   " << "pfk_: " << pfk_ << "\twsFactor: " << wsFactor
                 << "\trhFactor: " << rhFactor << "\tarrhenius: " << arrhenius << endl;
            cout << "PKM_hT   " << "HLK_: " << HLK_ << "\twcRatio_: " << wcRatio_ << "\twsRatio_: " << wsRatio_ << endl;
            cout << "PKM_hT   " << "ngrate: " << ngrate << "\thsrate: " << hsrate << "\tdiffrate: " << diffrate
                 << "\trate_ini: " << rate_ini << "\trate: " << rate << endl;
            cout << "PKM_hT   " << "DOR: " << DOR << "\tnewDOR: " << newDOR << "\ttotalDOR: " << totalDOR
                 << "\tinitScaledMass_: " << initScaledMass_ << "\tscaledMass_: " << scaledMass_
                 << "\tmassDissolved: " << massDissolved << endl;
            cout.flush();

            if (verbose_) {
                cout << "ParrotKillohModel::calculateKineticStep "
                     << "Original scaled mass = " << initScaledMass_
                     << ", dissolved scaled mass = " << massDissolved << endl;
                cout.flush();
            }

        } else {
            throw DataException("ParrotKillohModel", "calculateKineticStep",
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
        EOBException ex("ParrotKillohModel", "calculateKineticStep", oor.what(), 0, 0);
        ex.printException();
        exit(1);
    }

    return ;
}
*/

void ParrotKillohModel::calculateKineticStep (const double timestep, const double temperature,
                                             double rh, double &scaledMass,
                                             double &massDissolved, int cyc, double totalDOR) {

    ///
    /// Initialize local variables
    ///

    double T = temperature;
    double DOR, newDOR;

    double wsFactor = 1.0;
    double rhFactor = 1.0;
    double arrhenius = 1.0;
    double ngrate = 1.0e-10;   // Nucleation and growth rate
    double hsrate = 1.0e-10;   // Hydration shell rate
    double diffrate = 1.0e-10; // Diffusion rate
    double rate = 1.0e-10;     // Selected rate

    ///
    /// Determine if this is a normal step or a necessary
    /// tweak from a failed GEM_run call
    ///

    try {

        // @todo BULLARD PLACEHOLDER
        // Still need to implement constant gas phase composition
        // Will involve equilibrating gas with aqueous solution
        //
        // First step each iteration is to equilibrate gas phase
        // with the electrolyte, while forbidding anything new
        // from precipitating.

        if (verbose_) {
            cout << "ParrotKillohModel::calculateKineticStep for " << name_ << endl;
            cout.flush();
        }

        // RH factor is the same for all clinker phases

        /// This is a big kluge for internal relative humidity
        /// @note Using new gel and interhydrate pore size distribution model
        ///       which is currently contained in the Lattice object.
        ///
        /// Surface tension of water is gamma = 0.072 J/m2
        /// Molar volume of water is Vm = 1.8e-5 m3/mole
        /// The Kelvin equation is
        ///    p/p0 = exp (-4 gamma Vm / d R T) = exp (-6.23527e-7 / (d T))
        ///
        ///    where d is the pore diameter in meters and T is absolute temperature

        /// Assume a zero contact angle for now.
        /// @todo revisit the contact angle issue

        if (initScaledMass_ > 0.0) { //DOR @ t-1
            DOR = (initScaledMass_ - scaledMass_) / initScaledMass_;
            DOR = min(DOR, 0.99); // prevents DOR from prematurely stopping PK calculations

            if (verbose_) {
                cout << "~~~~>DOR for " << name_ << " = " << DOR << endl;
                cout.flush();
            }
        } else {
            throw FloatException("ParrotKillohModel", "calculateKineticStep",
                                 "initScaledMass_ = 0.0");
        }

        //cout << "PKM model for " << name_ << "    k1 = " << k1_ << "    k2 = " << k2_ << "    k3 = " << k3_
        //     << "    pfk = " << pfk_ << "    n1 = " << n1_ << "    n3 = " << n3_
        //     << "    Ea = " << activationEnergy_ << "    ssaFactor_ = " << ssaFactor_ << "    DOR = " << DOR << endl;
        //cout.flush();

        if (DOR < 1.0) {
            // Normal Parrott and Killoh implementation here

            if (fabs(n1_) > 0.0) {
                ngrate =
                    (k1_ / n1_) * (1.0 - DOR) * pow((-log(1.0 - DOR)), (1.0 - n1_));
                ngrate *= (ssaFactor_); // only used for the N+G rate

                if (ngrate < 1.0e-10)
                    ngrate = 1.0e-10;
            } else {
                throw FloatException("ParrotKillohModel", "calculateKineticStep",
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

            rhFactor = pow(((rh - 0.55) / 0.45), 4.0);
            arrhenius =
                exp((activationEnergy_ / GASCONSTANT) * ((1.0 / refT_) - (1.0 / T)));

            double rate_ini = rate;

            rate *= (pfk_ * rhFactor * arrhenius); // rate is R @ t-1

            double prod = rate * timestep;
            newDOR = DOR + prod;
            wsFactor = 1;

            if (totalDOR > HLK_ * wcRatio_){
                wsFactor = pow((1 + 3.333 * (HLK_ * wcRatio_ - newDOR)),4);
                //prod = timestep * rate * pow((1 + 3.333 * (HLK_ * wcRatio_ - newDOR)),4);
                prod = timestep * rate * wsFactor;
                newDOR = DOR + prod;
            }

            scaledMass_ = initScaledMass_ * (1.0 - newDOR);
            //massDissolved = (newDOR - DOR) * initScaledMass_;
            massDissolved = initScaledMass_ * prod;
            scaledMass = scaledMass_;

            cout << "****************** PKM_hT = " << timestep << "    cyc = " << cyc
                 << "    microPhaseId_ = " << microPhaseId_ << "    microPhase = " << name_
                 << "    GEMPhaseIndex = " << GEMPhaseId_ << " ******************" << endl;
            cout << "PKM_hT   " << "pfk_: " << pfk_ << "    k1 = " << k1_ << "    n1 = " << n1_
                 << "    k2 = " << k2_ << "    k3 = " << k3_ << "    n3 = " << n3_ << endl;
            cout << "PKM_hT   " << "HLK_: " << HLK_ << "    Ea = " << activationEnergy_ << endl;
            cout << "PKM_hT   " << "specificSurfaceArea_ = " << specificSurfaceArea_
                 << "    refSpecificSurfaceArea_ = " << refSpecificSurfaceArea_ << "    ssaFactor_ = " << ssaFactor_ << endl;
            cout << "PKM_hT   " << "wcRatio_: " << wcRatio_ << "\twsRatio_: " << wsRatio_ << endl;
            cout << "PKM_hT   " << "ngrate: " << ngrate << "\thsrate: " << hsrate << "\tdiffrate: " << diffrate
                 << "\trate_ini: " << rate_ini << "\trate: " << rate << endl;
            cout << "PKM_hT   " << "wsFactor: " << wsFactor << "\trhFactor: " << rhFactor << "\tarrhenius: " << arrhenius << endl;
            cout << "PKM_hT   " << "DOR: " << DOR << "\tnewDOR: " << newDOR << "\ttotalDOR: " << totalDOR
                 << "\tinitScaledMass_: " << initScaledMass_ << "\tscaledMass_: " << scaledMass_
                 << "\tmassDissolved: " << massDissolved << endl;
            cout.flush();

            if (verbose_) {
                cout << "ParrotKillohModel::calculateKineticStep "
                     << "Original scaled mass = " << initScaledMass_
                     << ", dissolved scaled mass = " << massDissolved << endl;
                cout.flush();
            }

        } else {
            throw DataException("ParrotKillohModel", "calculateKineticStep",
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
        EOBException ex("ParrotKillohModel", "calculateKineticStep", oor.what(), 0, 0);
        ex.printException();
        exit(1);
    }

    return ;
}


