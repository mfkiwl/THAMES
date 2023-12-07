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
    /// Default value for the rate constant
    ///
    
    rateconst_ = 0.0;

    ///
    /// Default value for the exponents in the rate equation
    ///
    
    siexp_ = 1.0;
    dfexp_ = 1.0;
    ohexp_ = 0.0;
    sio2_ = 1.0;
    al2o3_ = cao_ = 0.0;
    loi_ = 0.0;

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
        cout << "ParrotKillohModel::ParrotKillohModel Constructor" << endl;
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

    specificSurfaceArea_ = kineticData.ssa;
    refSpecificSurfaceArea_ = kineticData.refssa;
    ssaFactor_ = specificSurfaceArea_ / refSpecificSurfaceArea_;
    setSio2(kineticData.sio2);
    setAl2o3(kineticData.al2o3);
    setCao(kineticData.cao);
    setLoi(kineticData.loi);
    setRateconstant(kineticData.rateconst);
    setSiexp(kineticData.siexp);
    setDfexp(kineticData.dfexp);
    setOhexp(kineticData.ohexp);

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
                                            bool doTweak,
                                            double rh,
                                            vector<double> &ICMoles,
                                            vector<double> &solutICMoles,
                                            vector<double> &DCMoles,
                                            vector<double> &GEMPhaseMoles)
{
    ///
    /// Initialize local variables
    ///

    double T = temperature;
    double arrhenius = 1.0;

    double rate = 1.0e-10;            // Selected rate

    double massDissolved = 0.0;

    ///
    /// Determine if this is a normal step or a necessary
    /// tweak from a failed GEM_run call
    ///
    
    static double hyd_time = 0.0;
    if (!doTweak) hyd_time = hyd_time + timestep;

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
        vector<double> ICMoles,solutICMoles,DCMoles,GEMPhaseMoles;
        ICMoles.clear();
        ICMoles.resize(ICNum,0.0);
        solutICMoles.clear();
        solutICMoles.resize(ICNum,0.0);
        DCMoles.clear();
        DCMoles.resize(DCNum,0.0);
        GEMPhaseMoles.clear();
        GEMPhaseMoles.resize(numGEMPhases,0.0);
        string icn;

        vector<string> ICName;
        ICName.clear();
        ICName.resize(ICNum," ");

        int waterId = chemSys_->getDCId("H2O@");

        #ifdef DEBUG
            cout << "PozzolanicModel::calculateKineticStep ICmoles before dissolving:" << endl;
            cout << "PozzolanicModel::calculateKineticStep DC moles of water = " << chemSys_->getDCMoles(waterId);
            cout.flush();
            for (int i = 0; i < ICNum; i++) {
              if (isFirst) {
                  ICMoles[i] = 1.0e-9;
              } else {
                  ICMoles[i] = chemSys_->getICMoles(i);
              }
              ICName[i] = chemSys_->getICName(i);
              cout << "PozzolanicModel::calculateKineticStep     " << ICName[i] << ": " << ICMoles[i] << " mol" << endl;
            }
        #else 
            for (int i = 0; i < ICNum; i++) {
              ICMoles[i] = chemSys_->getICMoles(i);
              if (isFirst) {
                  ICMoles[i] = 1.0e-9;
              } else {
                  ICMoles[i] = chemSys_->getICMoles(i);
              }
              ICName[i] = chemSys_->getICName(i);
            }
        #endif

        for (int i = 0; i < DCNum; i++) {
          DCMoles[i] = chemSys_->getDCMoles(i);
        }
        for (int i = 0; i < numGEMPhases; i++) {
          GEMPhaseMoles[i] = chemSys_->getGEMPhaseMoles(i);
        }

        solutICMoles = chemSys_->getSolution();

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
          
          double DOH = 0.0;
          double newDOH = 0.0;

          if (initScaledMass_ > 0.0) {
              DOH = (initScaledMass_ - scaledMass_) /
                        (initScaledMass_);
              DOH = min(DOH,0.99);  // prevents DOH from prematurely stopping PK calculations
          } else {
              throw FloatException("PozzolanicModel","calculateKineticStep",
                             "initScaledMass_ = 0.0");
          }

          arrhenius = exp((activationEnergy_/GASCONSTANT)*((1.0/refT_) - (1.0/T)));

          if (DOH < 1.0 && !doTweak) {
                    
              /// @todo Make a rate equation here ...

              newDOH = 0.5;

              /// ... to here
                    
              scaledMass_ = initScaledMass_ * (1.0 - newDOH);
              massDissolved = (newDOH - DOH) * initScaledMass_;

                       
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

              for (int ii = 0; ii < ICMoles.size(); ii++) {

                  /// @todo BULLARD PLACEHOLDER
                  /// Special case for Silica-amorph here to account for SiO2 < 100%
                     

                  ICMoles[ii] += ((massDissolved
                                  / chemSys_->getDCMolarMass(DCId_))
                                  * chemSys_->getDCStoich(DCId_,ii));
                      
                  if (ICName[ii] == "O") {
                    // Dissolved K2O in this phase
                    if (chemSys_->isIC("K")) {
                        icn = "K";
                        molarMass = 2.0 * chemSys_->getICMolarMass(icn);
                        icn = "O";
                        molarMass += chemSys_->getICMolarMass(icn);
                        ICMoles[ii] += (impurityRelease[0]/molarMass);
                    }
                    // Dissolved Na2O in this phase
                    if (chemSys_->isIC("Na")) {
                        icn = "Na";
                        molarMass = 2.0 * chemSys_->getICMolarMass(icn);
                        icn = "O";
                        molarMass += chemSys_->getICMolarMass(icn);
                        ICMoles[ii] += (impurityRelease[1]/molarMass);
                    }
                    // Dissolved MgO in this phase
                    if (chemSys_->isIC("Mg")) {
                        icn = "Mg";
                        molarMass = chemSys_->getICMolarMass(icn);
                        icn = "O";
                        molarMass += chemSys_->getICMolarMass(icn);
                        ICMoles[ii] += (impurityRelease[2]/molarMass);
                    }
                    // Dissolved SO3 in this phase
                    if (chemSys_->isIC("S")) {
                        icn = "S";
                        molarMass = chemSys_->getICMolarMass(icn);
                        icn = "O";
                        molarMass += (3.0 * chemSys_->getICMolarMass(icn));
                        ICMoles[ii] += (3.0 * (impurityRelease[3]/molarMass));
                    }
                  } else if (ICName[ii] == "S") {
                      // Dissolved SO3  in this phase
                      icn = "S";
                      molarMass = chemSys_->getICMolarMass(icn);
                      icn = "O";
                      molarMass += (3.0 * chemSys_->getICMolarMass(icn));
                      ICMoles[ii] += (impurityRelease[3]/molarMass);
                  } else if (ICName[ii] == "K") {
                      // Dissolved K2O in this phase
                      icn = "K";
                      molarMass = 2.0 * chemSys_->getICMolarMass(icn);
                      icn = "O";
                      molarMass += chemSys_->getICMolarMass(icn);
                      ICMoles[ii] += (2.0 * (impurityRelease[0]/molarMass));
                  } else if (ICName[ii] == "Na") {
                      // Dissolved Na2O in this phase
                      icn = "Na";
                      molarMass = 2.0 * chemSys_->getICMolarMass(icn);
                      icn = "O";
                      molarMass += chemSys_->getICMolarMass(icn);
                      ICMoles[ii] += (2.0 * (impurityRelease[1]/molarMass));
                  } else if (ICName[ii] == "Mg") {
                      // Dissolved MgO in this phase
                      icn = "Mg";
                      molarMass = chemSys_->getICMolarMass(icn);
                      icn = "O";
                      molarMass += chemSys_->getICMolarMass(icn);
                      ICMoles[ii] += (impurityRelease[2]/molarMass);
                  }
    
              }

          } else if (DOH < 1.0) {

              ///
              /// We will just tweak the icmoles a bit to try to
              /// cure a previous failed convergence with GEM_run
              ///
             
              for (int ii = 0; ii < ICMoles.size(); ii++) {
                  if (ICName[ii] != "H" && ICName[ii] != "O" && ICName[ii] != "Zz") {
                      ICMoles[ii] *= 1.01;
                  }
              }

          } else {
              throw DataException("PozzolanicModel","calculateKineticStep",
                            "DOH >= 1.0");
          }   

          #ifdef DEBUG
            if (!doTweak) {
                cout << "PozzolanicModel::calculateKineticStep ICmoles after dissolving:" << endl;
            } else {
                cout << "PozzolanicModel::calculateKineticStep ICmoles after tweaking:" << endl;
            }
            for (int i = 0; i < ICNum; i++) {
              cout << "    " << ICName[i] << ": " << ICMoles[i] << " mol" << endl;
            }
            cout.flush();
          #endif

          if (doTweak) {
              for (int ii = 0; ii < ICMoles.size(); ii++) {
                  chemSys_->setICMoles(ii,ICMoles[ii]);
              }
              return;
          }

          for (int ii = 0; ii < ICMoles.size(); ii++) {
              chemSys_->setICMoles(ii,ICMoles[ii]);
              #ifdef DEBUG
                  cout << "PozzolanicModel::calculateKineticStep ICmoles of " << ICName[ii]
                       << " is: " << ICMoles[ii] << endl;
                  cout.flush();
              #endif
          }

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
