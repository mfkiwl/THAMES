/**
@file  ParrotKillohModel.cc
@brief Method definitions for the ParrotKillohModel class.

*/
#include "ParrotKillohModel.h"

ParrotKillohModel::ParrotKillohModel ()
{
    ///
    /// Default value for w/c ratio in PK model is 0.45
    ///

    wcRatio_ = 0.45;

    ///
    /// Default value for specific surface area in PK model is 385 m<sup>2</sup>/kg
    ///

    specificSurfaceArea_ = 385.0;
    refSpecificSurfaceArea_ = 385.0;     // reference SSA (m2/kg)

    ///
    /// Default temperature in the PK model is 20 C (or 293 K)
    ///

    temperature_ = 293.15;  // default temperature (K)
    refT_ = 293.15;         // default temperature (K)

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
    n1_ = 1.0;
    n3_ = 1.0;
    activationEnergy_ = 0.0;
    scaledMass_ = 0.0;
    initScaledMass_ = 0.0;
    critDOH_ = 100.0;
    degreeOfHydration_ = 0.0;
    ICNum_ = 0;
    ICName_.clear();
    DCNum_ = 0;
    DCName_.clear();
    GEMPhaseNum_ = 0;
    
    ///
    /// The default is to not have sulfate attack or leaching, so we set the default
    /// time for initiating these simulations to an absurdly large value: 10 billion
    /// days or 27 million years
    ///

    sulfateAttackTime_ = 1.0e10;
    leachTime_ = 1.0e10;


    return;
}

ParrotKillohModel::ParrotKillohModel (ChemicalSystem *cs,
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

    specificSurfaceArea_ = kineticData.specificSurfaceArea;
    refSpecificSurfaceArea_ = kineticData.refSpecificSurfaceArea;
    ssaFactor_ = specificSurfaceArea_ / refSpecificSurfaceArea_;

    ///
    /// Default value for w/c ratio in PK model is 0.45
    ///

    wcRatio_ = lattice_->getWsratio();

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
    k1_ = kineticData.k1;
    k2_ = kineticData.k2;
    k3_ = kineticData.k3;
    n1_ = kineticData.n1;
    n3_ = kineticData.n3;
    activationEnergy_ = kineticData.activationEnergy;
    scaledMass_ = kineticData.scaledMass;
    initScaledMass_ = kineticData.scaledMass;
    critDOH_ = kineticData.critDOH;
    degreeOfHydration_ = 0.0;
    
    waterId_ = chemSys_->getDCId("H2O@");
    ICNum_ = chemSys_->getNumICs();
    ICName_ = chemSys_->getICName();
    DCNum_ = chemSys_->getNumDCs();
    DCName_ = chemSys_->getDCName();
    GEMPhaseNum_ = chemSys_->getNumGEMPhases();

    /// The default is to not have sulfate attack or leaching, so we set the default
    /// time for initiating these simulations to an absurdly large value: 10 billion
    /// days or 27 million years
    ///

    sulfateAttackTime_ = 1.0e10;
    leachTime_ = 1.0e10;
    
                                   
    return;
}

void ParrotKillohModel::calculateKineticStep (const double timestep,
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
    double wcFactor = 1.0;
    double rhFactor = 1.0;
    double DOH = 0.0;
    double cDOH = 0.0;
    double arrhenius = 1.0;

    double ngrate = 1.0e-10;          // Nucleation and growth rate
    double hsrate = 1.0e-10;          // Hydration shell rate
    double diffrate = 1.0e-10;        // Diffusion rate
    double rate = 1.0e-10;            // Selected rate

    double newDOH = 0.0;              // Updated value of doh
    double massDissolved = 0.0;

    ///
    /// Determine if this is a normal step or a necessary
    /// tweak from a failed GEM_run call
    ///
    
    try {
        static int conc_index = 0;     
        int microPhaseId,DCId,ICId;
        double molarMass;
        GEMPhaseMoles.resize(GEMPhaseNum_,0.0);
        string icn;

        // @todo BULLARD PLACEHOLDER
        // Still need to implement constant gas phase composition
        // Will involve equilibrating gas with aqueous solution
        //
        // First step each iteration is to equilibrate gas phase
        // with the electrolyte, while forbidding anything new
        // from precipitating.
           
        #ifdef DEBUG
           cout << "ParrotKillohModel::calculateKineticStep for "
                << name_ << endl;
           cout.flush();
        #endif

        vector<double> impurityRelease;
        impurityRelease.clear();
        impurityRelease.resize(chemSys_->getNumMicroImpurities(),0.0);

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
          
        rh = rh > 0.55 ? rh : 0.551;
        rhFactor = pow(((rh - 0.55)/0.45),4.0);
          
        if (initScaledMass_ > 0.0) {
            DOH = (initScaledMass_ - scaledMass_) /
                      (initScaledMass_);
            DOH = min(DOH,0.99);  // prevents DOH from prematurely
                                  // stopping PK calculations
            #ifdef DEBUG
              cout << "~~~~>DOH for " << name_ << " = " << DOH << endl;
              cout.flush();
            #endif // DEBUG
        } else {
            throw FloatException("ParrotKillohModel","calculateKineticStep",
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

        cout << "PK model for " << name_ << endl;
        cout << "    k1 = " << k1_ << endl;
        cout << "    k2 = " << k2_ << endl;
        cout << "    k3 = " << k3_ << endl;
        cout << "    n1 = " << n1_ << endl;
        cout << "    n3 = " << n3_ << endl;
        cout << "    Ea = " << activationEnergy_ << endl;
        cout.flush();

        if (DOH < 1.0) {
                    
            // Normal Parrott and Killoh implementation here
                    
            if (fabs(n1_) > 0.0) {
                ngrate = (k1_/n1_) * (1.0 - DOH)
                               * pow((-log(1.0 - DOH)),(1.0 - n1_));
                ngrate *= (ssaFactor_);  // only used for the N+G rate
            
                if (ngrate < 1.0e-10) ngrate = 1.0e-10;
            } else {
                throw FloatException("ParrotKillohModel","calculateKineticStep",
                                       "n1_ = 0.0");
            }
        
            hsrate = k3_ * pow((1.0 - DOH),n3_);
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
                cout << "ParrotKillohModel::calculateKineticStep "
                     << "Original scaled mass = " << initScaledMass_
                     << ", dissolved scaled mass = " << massDissolved << endl;
                cout << "New scaled mass = "
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

#ifdef DEBUG
            cout << "!!!Before dissolving some " << name_ << endl;
            for (int ii = 0; ii < dICMoles.size(); ii++) {
              cout << "dIC moles of " << ICName_[ii] << " = " << dICMoles[ii] << endl;
            }
            cout.flush();
#endif
            for (int ii = 0; ii < dICMoles.size(); ii++) {

                dICMoles[ii] += ((massDissolved
                                / chemSys_->getDCMolarMass(DCId_))
                                * chemSys_->getDCStoich(DCId_,ii));
                      
                if (ICName_[ii] == "O") {
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
                } else if (ICName_[ii] == "S") {
                    // Dissolved SO3  in this phase
                    icn = "S";
                    molarMass = chemSys_->getICMolarMass(icn);
                    icn = "O";
                    molarMass += (3.0 * chemSys_->getICMolarMass(icn));
                    dICMoles[ii] += (impurityRelease[3]/molarMass);
                } else if (ICName_[ii] == "K") {
                    // Dissolved K2O in this phase
                    icn = "K";
                    molarMass = 2.0 * chemSys_->getICMolarMass(icn);
                    icn = "O";
                    molarMass += chemSys_->getICMolarMass(icn);
                    dICMoles[ii] += (2.0 * (impurityRelease[0]/molarMass));
                } else if (ICName_[ii] == "Na") {
                    // Dissolved Na2O in this phase
                    icn = "Na";
                    molarMass = 2.0 * chemSys_->getICMolarMass(icn);
                    icn = "O";
                    molarMass += chemSys_->getICMolarMass(icn);
                    dICMoles[ii] += (2.0 * (impurityRelease[1]/molarMass));
                } else if (ICName_[ii] == "Mg") {
                    // Dissolved MgO in this phase
                    icn = "Mg";
                    molarMass = chemSys_->getICMolarMass(icn);
                    icn = "O";
                    molarMass += chemSys_->getICMolarMass(icn);
                    dICMoles[ii] += (impurityRelease[2]/molarMass);
                }
    
            }

        } else {
            throw DataException("ParrotKillohModel","calculateKineticStep",
                          "DOH >= 1.0");
        }   

//        #ifdef DEBUG
//          if (!doTweak) {
//              cout << "ParrotKillohModel::calculateKineticStep ICmoles after dissolving:" << endl;
//          } else {
//              cout << "ParrotKillohModel::calculateKineticStep ICmoles after tweaking:" << endl;
//          }
//          for (int i = 0; i < ICNum_; i++) {
//            cout << "    " << ICName_[i] << ": " << ICMoles[i] << " mol" << endl;
//          }
//          cout.flush();
//        #endif

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
        EOBException ex("ParrotKillohModel","calculateKineticStep",
                           oor.what(),0,0);
        ex.printException();
        exit(1);
    }

	
    return;
}
