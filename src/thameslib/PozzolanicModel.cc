/**
@file  PozzolanicModel.cc
@brief Method definitions for the PozzolanicModel class.

*/
#include "PozzolanicModel.h"

PozzolanicModel::PozzolanicModel ()
{
    ///
    /// Default value for specific surface area is 385 m<sup>2</sup>/kg
    ///

    specificSurfaceArea_ = 385.0;
    refSpecificSurfaceArea_ = 385.0;     // reference specific surface area (m2/kg)
    ssaFactor_ = 1.0;

    ///
    /// Default temperature in the PK model is 20 C (or 293 K)
    ///

    temperature_ = 293.15;  // default temperature (K)
    refT_ = 293.15;         // default temperature (K)

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
    /// The default is to not have sulfate attack or leaching, so we set the default
    /// time for initiating these simulations to an absurdly large value: 10 billion
    /// days or 27 million years
    ///

    sulfateAttackTime_ = 1.0e10;
    leachTime_ = 1.0e10;


    return;
}

PozzolanicModel::PozzolanicModel (ChemicalSystem *cs,
                 Solution *solut,
                 Lattice *lattice,
                 struct KineticData &kineticData,
                 const bool verbose,
                 const bool warning)
{

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
    solut_ = solut;
    lattice_ = lattice;


    ///
    /// Default value for specific surface area in PK model is 385 m<sup>2</sup>/kg
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
    setSiexp(kineticData.siexp);
    setDfexp(kineticData.dfexp);
    setOhexp(kineticData.ohexp);
    lossOnIgnition_ = kineticData.loi;

    ///
    /// Default initial solid mass is 100 g
    ///
    
    initSolidMass_ = 100.0;

    ///
    /// Default temperature in the PK model is 20 C (or 293 K)
    ///

    temperature_ = kineticData.temperature;
    refT_ = kineticData.reftemperature;

    name_ = kineticData.name;
    microPhaseId_ = kineticData.microPhaseId;
    DCId_ = kineticData.DCId;
    GEMPhaseId_ = kineticData.GEMPhaseId;;
    activationEnergy_ = kineticData.activationEnergy;
    scaledMass_ = kineticData.scaledMass;
    initScaledMass_ = kineticData.scaledMass;

    ///
    /// The default is to not have sulfate attack or leaching, so we set the default
    /// time for initiating these simulations to an absurdly large value: 10 billion
    /// days or 27 million years
    ///

    sulfateAttackTime_ = 1.0e10;
    leachTime_ = 1.0e10;
    
                                   
    return;
}

void PozzolanicModel::calculateKineticStep (const double timestep,
                                            const double temperature,
                                            bool isFirst,
                                            double rh,
                                            vector<double> &dICMoles,
                                            vector<double> &dsolutICMoles,
                                            vector<double> &DCMoles,
                                            vector<double> &GEMPhaseMoles)
{
    ///
    /// Initialize local variables
    ///

    double T = temperature;
    double arrhenius = 1.0;

    double dissrate = 1.0e9;          // Nucleation and growth rate
    double hsrate = 1.0e9;          // Hydration shell rate
    double diffrate = 1.0e9;        // Diffusion rate

    double rate = 1.0e-10;            // Selected rate

    double massDissolved = 0.0;

    ///
    /// Determine if this is a normal step or a necessary
    /// tweak from a failed GEM_run call
    ///
    
    static double hyd_time = 0.0;
    hyd_time = hyd_time + timestep;

    #ifdef DEBUG
        cout << "PozzolanicModel::calculateKineticStep Hydration Time = "
             << hyd_time << endl;
        cout.flush();
    #endif

    try {
        static int conc_index = 0;     
        int ICNum = chemSys_->getNumICs();
        int DCNum = chemSys_->getNumDCs();
        int numGEMPhases = chemSys_->getNumGEMPhases();
        int DCId,ICId;
        double molarMass;
        dICMoles.clear();
        dICMoles.resize(ICNum,0.0);
        dsolutICMoles.clear();
        dsolutICMoles.resize(ICNum,0.0);
        GEMPhaseMoles.clear();
        GEMPhaseMoles.resize(numGEMPhases,0.0);
        string icn;

        vector<string> ICName;
        ICName.clear();
        ICName.resize(ICNum," ");

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
           
          #ifdef DEBUG
             cout << "PozzolanicModel::calculateKineticStep "
                  << "Looping over kinetically controlled phases.  " << endl;
             cout.flush();
          #endif

          vector<double> impurityRelease;
          impurityRelease.clear();
          impurityRelease.resize(chemSys_->getNumMicroImpurities(),0.0);

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
          ///    where d is the pore diameter in meters and T is absolute temperature
          
          /// Assume a zero contact angle for now.
          /// @todo revisit the contact angle issue
          
          double critporediam = lattice_->getLargestSaturatedPore(); // in nm
          critporediam *= 1.0e-9;                                    // in m
          double rh = exp(-6.23527e-7/critporediam/T);

          /// Assume a zero contact angle for now.
          /// @todo revisit the contact angle issue
          
          rh = rh > 0.55 ? rh : 0.551;
          double rhFactor = rh;
          // rhFactor = pow(((rh - 0.55)/0.45),4.0);
          
          int thisGEMPhase = chemSys_->getMicroPhaseToGEMPhase(microPhaseId_,0);
          double saturationIndex = solut_->getSI(thisGEMPhase);
          double DOR = 0.0;
          double newDOR = 0.0;

          if (initScaledMass_ > 0.0) {
              DOR = (initScaledMass_ - scaledMass_) /
                        (initScaledMass_);
              DOR = min(DOR,0.99);  // prevents DOR from prematurely stopping PK calculations
          } else {
              throw FloatException("PozzolanicModel","calculateKineticStep",
                             "initScaledMass_ = 0.0");
          }

          arrhenius = exp((activationEnergy_/GASCONSTANT)*((1.0/refT_) - (1.0/T)));

          if (DOR < 1.0) {
                    
              /// @todo Modify dissolution rate equation for hydroxyl activity
              /// and loi

              if (fabs(dissolutionRateConst_) > 0.0) {
                
                /// BULLARD placeholder
                /// @note playing with different base rate constants here
               
                // double baserateconst = 5.6e-3 * arrhenius;  // mol m-2 s-1
                                                                       //

                double baserateconst = dissolutionRateConst_ * arrhenius;

                double ca = chemSys_->getDCConcentration("Ca+2");
                double kca = 4.0e-7; // mol m-2 s-1 ads. rate const for Ca (guess)
                double Kca = 10.0; // adsorption equilibrium constant is a guess
                double na = chemSys_->getDCConcentration("Na+");
                double kna = 6.35e-7; // mol m-2 s-1 ads. rate const from Dove and Crerar
                double Kna = 58.3; // adsorption equilibrium constant from Dove and Crerar
                double k = chemSys_->getDCConcentration("K+");
                double kk = 5.6e-7; // mol m-2 s-1 ads. rate const from Dove and Crerar
                double Kk = 46.6; // adsorption equilibrium constant from Dove and Crerar

                // Langmuir adsorption isotherms assumed to be additive
               
                baserateconst += (kca * Kca * ca / (1.0 + (Kca * ca)));
                baserateconst += (kna * Kna * na / (1.0 + (Kna * na)));
                baserateconst += (kk * Kk * k / (1.0 + (Kk * k)));

                double ohActivity = chemSys_->getDCActivity("OH-");
                double area = (specificSurfaceArea_ / 1000.0) * scaledMass_; // m2
                double molarmass = chemSys_->getDCMolarMass("Amor-Sl"); // g mol-1
                int GEMPhaseIndex = chemSys_->getMicroPhaseToGEMPhase(microPhaseId_,0);

                // Saturation index , but be sure that there is only one GEM Phase
                /// @note Assumes there is only one phase in this microstructure
                /// component
                /// @todo Generalize to multiple phases in a component (how?)
               
                double saturationIndex = solut_->getSI(GEMPhaseIndex);

                // activity of water
                double waterActivity = chemSys_->getDCActivity(chemSys_->getDCId("H2O@"));

                // Make sure sio2 is given in mole fraction silica fume
                // This equation basically implements the Dove and Crerar rate equation
                // for quartz.  Needs to be calibrated for silica fume, but hopefully
                // the BET area and LOI will help do that.
               
                dissrate = baserateconst * rhFactor * pow(ohActivity,ohexp_) * area * pow(waterActivity,2.0)
                       * (1.0 - (lossOnIgnition_/100.0)) * (sio2_)
                       * pow((1.0 - pow(saturationIndex,siexp_)),dfexp_); 
                            
                 cout << "PozzolanicModel::calculateKineticStep for " << name_ << endl;
                 cout << "  dissrate = " << dissrate << endl;
                 cout << "    (baserateconst = " << baserateconst << ")"  << endl;
                 cout << "    (rhFactor = " << rhFactor << ")" << endl;
                 cout << "    (arrhenius = " << arrhenius << ")" << endl;
                 cout << "    (ohactivity = " << ohActivity << ")"  << endl;
                 cout << "    (waterterm = " << pow(waterActivity,2.0) << ")"  << endl;
                 cout << "    (area = " << area << ")"  << endl;
                 cout << "    (saturationIndex = " << saturationIndex << ")" << endl;
                 cout << "    (LOI = " << lossOnIgnition_ << ")" << endl;
                 cout << "    (sio2 = " << sio2_ << ")" << endl;
                 cout.flush();
                 if (dissrate < 1.0e-10) dissrate = 1.0e-10;
              } else {
                  throw FloatException("PozzolanicModel","calculateKineticStep",
                                       "diffusionRateConst_ = 0.0");
              }
        
              /// @note Need to get some more constants in here
              /// for diffusion coefficient, thickness, etc.
             
              hsrate = diffusionRateConstLate_ * ssaFactor_
                     * pow((1.0 - DOR),dorexp_) * pow((1.0 - pow(saturationIndex,siexp_)),dfexp_);
              if (hsrate < 1.0e-10) hsrate = 1.0e-10;

              /// @note Temporarily eliminate early age diffusion model
              if (DOR > 0.0) {
                  // diffrate = (diffusionRateConstEarly_
                  //          * ssaFactor_ * pow((1.0 - pow(saturationIndex,siexp_)),dfexp_)
                  //          * pow((1.0 - pow(DOR,dorexp_)),(2.0/3.0))) /
                  //                (1.0 - pow((1.0 - pow(DOR,dorexp_)),(1.0/3.0)));
                  // if (diffrate < 1.0e-10) diffrate = 1.0e-10;
                  diffrate = 1.0e9;
              } else {
                  diffrate = 1.0e9;
              }

              rate = (dissrate < hsrate) ? dissrate : hsrate;
              if (diffrate < rate) rate = diffrate;
              rate *= (rhFactor * arrhenius);
              newDOR = DOR + (rate * timestep);

              cout << "  hsrate = " << hsrate << endl;
              cout << "    (baserateconst = " << diffusionRateConstLate_ << ")"  << endl;
              cout << "    (dor = " << DOR << ")"  << endl;
              cout << "  diffrate = " << diffrate << endl;
              cout << "  RATE = " << rate << endl;
              cout << "  timestep " << timestep << endl;
              cout << "  oldDOR = " << DOR << ", new DOR = "
                   << newDOR << endl;
              cout.flush();


              /// ... to here
                    
              scaledMass_ = max(initScaledMass_ * (1.0 - newDOR),0.0);
              massDissolved = max((newDOR - DOR) * initScaledMass_,0.0);

                       
              chemSys_->setMicroPhaseMass(microPhaseId_,scaledMass_);
              chemSys_->setMicroPhaseMassDissolved(microPhaseId_,massDissolved);

              #ifdef DEBUG
                  cout << "PozzolanicModel::calculateKineticStep Original scaled mass = " << initScaledMass_
                       << " and new scaled mass = "
                       << chemSys_->getMicroPhaseMass(microPhaseId_)
                       << " and new volume = "
                       << chemSys_->getMicroPhaseVolume(microPhaseId_) << endl;
                  cout.flush();
              #endif

              /// @note impurityRelease index values are assumed to
              /// be uniquely associated with particular chemical
              /// elements
              
              /// @todo Make this more general so that any indexing
              /// can be used
               
              /// @todo Allow any IC element to be an impurity, not
              /// necessarily just the ones hard-coded here
                   
              impurityRelease[0] = (massDissolved *
                      chemSys_->getK2o(microPhaseId_));
              impurityRelease[1] = (massDissolved *
                      chemSys_->getNa2o(microPhaseId_));
              impurityRelease[2] = (massDissolved *
                      chemSys_->getMgo(microPhaseId_));
              impurityRelease[3] = (massDissolved *
                      chemSys_->getSo3(microPhaseId_));

              for (int ii = 0; ii < dICMoles.size(); ii++) {

                  dICMoles[ii] += ((massDissolved
                                  / chemSys_->getDCMolarMass(DCId_))
                                  * chemSys_->getDCStoich(DCId_,ii));
                      
                  if (ICName[ii] == "O") {
                    // Dissolved K2O in this phase
                    if (chemSys_->isIC("K")) {
                        icn = "K";
                        molarMass = 2.0 * chemSys_->getICMolarMass(icn);
                        icn = "O";
                        molarMass += chemSys_->getICMolarMass(icn);
                        dICMoles[ii] += (impurityRelease[0]/molarMass);
                    }
                    // Dissolved Na2O in this phase
                    if (chemSys_->isIC("Na")) {
                        icn = "Na";
                        molarMass = 2.0 * chemSys_->getICMolarMass(icn);
                        icn = "O";
                        molarMass += chemSys_->getICMolarMass(icn);
                        dICMoles[ii] += (impurityRelease[1]/molarMass);
                    }
                    // Dissolved MgO in this phase
                    if (chemSys_->isIC("Mg")) {
                        icn = "Mg";
                        molarMass = chemSys_->getICMolarMass(icn);
                        icn = "O";
                        molarMass += chemSys_->getICMolarMass(icn);
                        dICMoles[ii] += (impurityRelease[2]/molarMass);
                    }
                    // Dissolved SO3 in this phase
                    if (chemSys_->isIC("S")) {
                        icn = "S";
                        molarMass = chemSys_->getICMolarMass(icn);
                        icn = "O";
                        molarMass += (3.0 * chemSys_->getICMolarMass(icn));
                        dICMoles[ii] += (3.0 * (impurityRelease[3]/molarMass));
                    }
                  } else if (ICName[ii] == "S") {
                      // Dissolved SO3  in this phase
                      icn = "S";
                      molarMass = chemSys_->getICMolarMass(icn);
                      icn = "O";
                      molarMass += (3.0 * chemSys_->getICMolarMass(icn));
                      dICMoles[ii] += (impurityRelease[3]/molarMass);
                  } else if (ICName[ii] == "K") {
                      // Dissolved K2O in this phase
                      icn = "K";
                      molarMass = 2.0 * chemSys_->getICMolarMass(icn);
                      icn = "O";
                      molarMass += chemSys_->getICMolarMass(icn);
                      dICMoles[ii] += (2.0 * (impurityRelease[0]/molarMass));
                  } else if (ICName[ii] == "Na") {
                      // Dissolved Na2O in this phase
                      icn = "Na";
                      molarMass = 2.0 * chemSys_->getICMolarMass(icn);
                      icn = "O";
                      molarMass += chemSys_->getICMolarMass(icn);
                      dICMoles[ii] += (2.0 * (impurityRelease[1]/molarMass));
                  } else if (ICName[ii] == "Mg") {
                      // Dissolved MgO in this phase
                      icn = "Mg";
                      molarMass = chemSys_->getICMolarMass(icn);
                      icn = "O";
                      molarMass += chemSys_->getICMolarMass(icn);
                      dICMoles[ii] += (impurityRelease[2]/molarMass);
                  }
    
              }

          } else {
              throw DataException("PozzolanicModel","calculateKineticStep",
                            "DOR >= 1.0");
          }   

          #ifdef DEBUG
            cout << "PozzolanicModel::calculateKineticStep ICmoles after dissolving:" << endl;
            for (int i = 0; i < ICNum; i++) {
              cout << "    " << ICName[i] << ": " << dICMoles[i] << " mol" << endl;
            }
            cout.flush();
          #endif

        } // End of normal hydration block
    }     // End of try block


    catch (EOBException eex) {
        eex.printException();
        exit(1);
    }
    catch (DataException dex) { 
        dex.printException();
        exit(1);
    }
    catch (FloatException fex) {
        fex.printException();
        exit(1);
    }
    catch (out_of_range &oor) {
        EOBException ex("PozzolanicModel","calculateKineticStep",
                           oor.what(),0,0);
        ex.printException();
        exit(1);
    }

	
    return;
}
