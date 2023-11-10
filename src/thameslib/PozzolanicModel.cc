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
    microPhaseId = 2;
    DCId_ = 2;
    GEMPhaseId_ = 2;
    RdICId_.clear();
    Rd_.clear();
    activationEnergy_ = 0.0;
    scaledMass_ = 0.0;
    initScaledMass_ = 0.0;
    degreeOfHydration_ = 0.0;

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
:chemSys_(cs),solut_(solut),lattice_(lattice)
{
    // Set the verbose and warning flags
   
    verbose_ = verbose;
    warning_ = warning;
    #ifdef DEBUG
        verbose_ = true;
        warning_ = true;
        cout << "PozzolanicModel::PozzolanicModel Constructor" << endl;
        cout.flush();
    #else
        verbose_ = verbose;
        warning_ = warning;
    #endif

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
    setRateconst(kineticData.rateconst);
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
    Rd_ = kineticData.Rd;
    RdICId_ = kineticData.RdICId;
    activationEnergy_ = kineticData.activationEnergy;
    scaledMass_ = kineticData.scaledMass;
    initScaledMass_ kineticData.initScaledMass;
    degreeOfHydration_ = 0.0;

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
    
    bool doTweak = (chemSys_->getTimesGEMFailed() > 0) ? true : false;

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
        int microPhaseId,DCId,ICId;
        double molarMass,Rd;
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

        if (isFirst) {  // Beginning of special first-time setup tasks
            
            // Set the proper amount of water for the total solid mass
            // and the water-solid ratio
            
            // Determine total initial solid mass
            double solidMass = 0.0;
            vector<double> initSolidMasses = getInitScaledMass();
            for (int i = 0; i < initSolidMasses.size(); i++) {
                solidMass += initSolidMasses[i];
            }
            double waterMass = solidMass * getWcRatio();
            double waterMolarMass = chemSys_->getDCMolarMass(waterId);
            double waterMoles = waterMass / waterMolarMass;
            if (verbose_) {
                cout << "PozzolanicModel::calculateKineticStep *** Initial solid mass = "
                     << solidMass << endl;
                cout << "PozzolanicModel::calculateKineticStep *** w/s ratio = "
                     << getWcRatio() << endl;
                cout << "PozzolanicModel::calculateKineticStep *** Initial water mass = "
                     << waterMass << endl;
                cout << "PozzolanicModel::calculateKineticStep *** Initial water moles = "
                     << waterMoles << endl;
                cout << "PozzolanicModel::calculateKineticStep ***" << endl;
                cout.flush();
            }

            if (waterMass <= 0.0) {
                throw FloatException("Controller","calculateState",
                                     "Divide by zero error");
            }

            // Print out the initial volumes of microstructure phases
           
            microPhaseId = 0;
            double psMass,psVolume;
            double volume = 0.0;
            if (verbose_) {
                cout << "PozzolanicModel::calculateKineticStep Initial MICROSTRUCTURE phase amount:" << endl;
                psMass = chemSys_->getMicroPhaseMass(microPhaseId_);
                psVolume = chemSys_->getMicroPhaseVolume(microPhaseId_);
                cout << "PozzolanicModel::calculateKineticStep     " << chemSys_->getMicroPhaseName(microPhaseId_)
                     << " (" << microPhaseId_ << "): mass = " << psMass
                     << ", vol = " << psVolume << endl;
                cout.flush();
            }

            for (int i = 0; i < ICMoles.size(); i++) {
                if (ICName[i] == "H") {
                    #ifdef DEBUG
                        cout << "PozzolanicModel::calculateKineticStep Previous IC moles for H is: "
                             << ICMoles[i] << endl;
                        cout.flush();
                    #endif
                    ICMoles[i] = (2.0 * waterMoles);
                    solut_->setICMoles(i,ICMoles[i]);
                    chemSys_->setDCMoles(waterId,waterMoles);
                    #ifdef DEBUG
                    if (verbose_) {
                        cout << "PozzolanicModel::calculateKineticStep New ICmoles for H is: "
                             << ICMoles[i] << endl;
                        cout.flush();
                    }
                    #endif
                }
                if (ICName[i] == "O") {
                    #ifdef DEBUG
                        cout << "PozzolanicModel::calculateKineticStep Previous IC moles for O is: "
                             << ICMoles[i] << endl;
                        cout.flush();
                    #endif
                    ICMoles[i] = waterMoles;
                    solut_->setICMoles(i,ICMoles[i]);
                    #ifdef DEBUG
                        cout << "PozzolanicModel::calculateKineticStep New ICmoles for O is: "
                             << ICMoles[i] << endl;
                        cout.flush();
                    #endif
               }
            }
            double wmv = chemSys_->getNode()->DC_V0(chemSys_->getDCId("H2O@"),
                                                    chemSys_->getP(),
                                                    chemSys_->getTemperature());
            chemSys_->setGEMPhaseMass(chemSys_->getGEMPhaseId("aq_gen"),waterMass);
            chemSys_->setGEMPhaseVolume(chemSys_->getGEMPhaseId("aq_gen"),wmv/waterMoles);


            // Modify initial pore solution composition if desired
            // input units are mol/kgw
            // watermass is in units of grams, not kg
            
            double kgWaterMass = waterMass / 1000.0;
            
            map<int,double> isComp = chemSys_->getInitialSolutionComposition();
            map<int,double>::iterator p = isComp.begin();

            while (p != isComp.end()) {
                if (verbose_) {
                    cout << "PozzolanicModel::calculateKineticStep "
                         << "modifying initial pore solution" << endl;
                    cout.flush();
                }
                #ifdef DEBUG
                    cout << "PozzolanicModel::calculateKineticStep "
                         << "--->Adding " << p->second
                         << " mol/kgw of "
                         << ICName[p->first] << " to initial solution." << endl;
                    cout.flush();
                #endif
                ICMoles[p->first] += (p->second * kgWaterMass);
                p++;
            }

            // @todo BULLARD PLACEHOLDER
            // Still need to implement constant gas phase composition
           
            // @todo BULLARD PLACEHOLDER
            // While we are still using PK model for clinker phases, we
            // must scan for and account for pozzolanic reactive components
            // that will alter the hydration rate of clinker
            // phases, presumably by making a denser hydration shell
            // with slower transport

            // @todo BULLARD PLACEHOLDER
            // Currently silica fume is the only pozzolanically reactive
            // component
            
            // @todo Find a way to make this general to all pozzolans
              
        }   // End of special first-time tasks

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
          rhFactor = pow(((rh - 0.55)/0.45),4.0);
          

          if (initScaledMass > 0.0) {
              DOH = (initScaledMass_ - scaledMass_) /
                        (initScaledMass_);
              DOH = min(DOH,0.99);  // prevents DOH from prematurely stopping PK calculations
          } else {
              throw FloatException("PozzolanicModel","calculateKineticStep",
                             "initScaledMass_ = 0.0");
          }

          cDOH = 1.333 * wcRatio_;
          wcFactor = 1.0;
          if (DOH > cDOH) {
              wcFactor += ((4.444 * wcRatio_) -
                           (3.333 * DOH));
              wcFactor = pow(wcFactor,4.0);
          }

          /*
          wcFactor = 1.0 + (3.333 *
                 (pow(((critDOH_[i] * wcRatio_) - DOH),4.0)));
          */

          arrhenius = exp((activationEnergy_/GASCONSTANT)*((1.0/refT_) - (1.0/T)));

          if (DOH < 1.0 && !doTweak) {
                    
              // @todo BULLARD PLACEHOLDER
              // Handle silica fume as a special case here:
              // Only for silica fume, we let k1 = BET surface area in m2/g,
              // k2 = LOI in percent by solid mass.
              // critDOH = silica content in PERCENT BY MASS of silica fume
              //
              // This is a TOTAL KLUGE right now from here...

              cDOH = 1.333 * wcRatio_;
              wcFactor = 1.0;
              if (DOH > cDOH) {
                  wcFactor += ((4.444 * wcRatio_) -
                               (3.333 * DOH));
                  wcFactor = pow(wcFactor,4.0);
              }

              // Normal Parrott and Killoh implementation here
                    
              // @todo BULLARD PLACEHOLDER pfk1_, pfk2_, pfk3_ are
              // pozzolanic modification factors for the Parrott and Killoh
              // rate constants k1, k2, and k3.
        
              if (fabs(n1_) > 0.0) {
                  ngrate = (k1_/n1_) * (1.0 - DOH)
                                 * pow((-log(1.0 - DOH)),(1.0 - n1_));
                  ngrate *= (ssaFactor_);  // only used for the N+G rate
              
                  if (ngrate < 1.0e-10) ngrate = 1.0e-10;
              } else {
                  throw FloatException("PozzolanicModel","calculateKineticStep",
                                         "n1_ = 0.0");
              }
        
              hsrate = k3_ * pow((1.0 - DOH),n3_[i]);
              if (hsrate < 1.0e-10) hsrate = 1.0e-10;

              if (DOH > 0.0) {
                  diffrate = (k2_ * pow((1.0 - DOH),(2.0/3.0))) /
                                       (1.0 - pow((1.0 - DOH),(1.0/3.0)));
                  if (diffrate < 1.0e-10) diffrate = 1.0e-10;
              } else {
                  diffrate = 1.0e9;
              }

              rate = (ngrate < hsrate) ? ngrate : hsrate;
              if (diffrate < rate) rate = diffrate;
              rate *= (wcFactor * rhFactor * arrhenius);
              newDOH = DOH + (rate * timestep);

              cout << "PK model for " << name_
                   << ", ngrate = " << ngrate
                   << ", hsrate = " << hsrate
                   << ", diffrate = " << diffrate
                   << ", rhFactor = " << rhFactor
                   << ", wcFactor = " << wcFactor
                   << ", RATE = " << rate
                   << ", timestep " << timestep
                   << ", oldDOH = " << DOH << ", new DOH = "
                   << newDOH << endl;
              cout.flush();

              /// @note This where we can figure out the volume dissolved
              /// and link it back to the current volume to see how many
              /// voxels need to dissolve
              ///

              /// @note This all depends on concept of degree of
              /// hydration as defined by the PK model 

              /// @todo Make this independent of PK model
                    
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
