/**
@file  KineticModel.cc
@brief Method definitions for the KineticModel class.

*/
#include "KineticModel.h"

KineticModel::KineticModel ()
{
    ///
    /// Clear out the vectors so they can be populated with values from the
    /// XML input file
    ///

    numPhases_ = 0;
    name_ = "";
    microPhaseId_ = 0;
    DCId_ = 0;
    GEMPhaseId_ = 0;
    RdICId_.clear();
    Rd_.clear();
    activationEnergy_ = 0.0;
    scaledMass_ = 0.0;
    initScaledMass_ = 0.0;
    waterId_ = 1;
    ICNum_ = 0;
    DCNum_ = 0;
    ICName_.clear();
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

void KineticModel::calculatePhaseChange (double k,
                                         double gamma,
                                         double timestep)
{
    ///
    /// Initialize local variables
    ///

    double DCChange = 0.0;
    vector<int> GEMPhaseIndex;
    GEMPhaseIndex.clear();
    GEMPhaseIndex = chemSys_->getMicroPhaseToGEMPhase(microPhaseId_);

    ///
    /// Should the next block be only for verbose output
    ///

    #ifdef DEBUG
        cout << "ParrotKillohModel::calculatePhaseChange Microphase "
             << chemSys_->getMicroPhaseName(microPhaseId_) << " contains phases: "
             << endl;
        for (int i = 0; i < GEMPhaseIndex.size(); i++) {
            cout << "ParrotKillohModel::calculatePhaseChange            "
                 << chemSys_->getGEMPhaseName(i) << endl;
        }    
        cout.flush();
    #endif

    double GEMPhaseMoles = 0.0;
    double *GEMPhaseMolesArray = chemSys_->getGEMPhaseMoles();
    for (int i = 0; i < GEMPhaseIndex.size(); i++) {
      GEMPhaseMoles += GEMPhaseMolesArray[i];
    }

    vector<double> GEMPhaseFrac;
    GEMPhaseFrac.clear();
    GEMPhaseFrac.resize(GEMPhaseIndex.size(), 0.0);

    for (int i = 0; i < GEMPhaseIndex.size(); i++) {
      double A = 0.0; // A is the surface area
      GEMPhaseFrac[i] = chemSys_->getGEMPhaseMoles(GEMPhaseIndex[i]) / GEMPhaseMoles;
      A = lattice_->getSurfaceArea(microPhaseId) * GEMPhaseFrac[i];
      double SI = 0.0; // SI is the saturation index for the phase i
      SI = solut_->getSI(GEMPhaseIndex[i]);
      
      vector<int> DCIndex;
      DCIndex.clear();
      DCIndex = chemSys_->getGEMPhaseDCMembers(GEMPhaseIndex[i]);
      #ifdef DEBUG
          cout << "ParrotKillohModel::calculatePhaseChange Phase "
               << chemSys_->getGEMPhaseName(GEMPhaseIndex[i]) 
               << " contains DC members: " << endl;
      #endif
      for (int j = 0; j < DCIndex.size(); j++) {
        chemSys_->getDCName(DCIndex[j]);
      }
      

      vector<double> DCFrac;
      DCFrac.clear();
      DCFrac.resize(DCIndex.size(), 0.0);
      double DCMoles = 0.0;
      for (int j = 0; j < DCIndex.size(); j++) {
        DCMoles += chemSys_->getDCMoles(DCIndex[j]);
      }
      
      double surfaceArea = 0.0;
      for (int j = 0; j < DCIndex.size(); j++) {
        DCFrac[j] = chemSys_->getDCMoles(DCIndex[j]) / DCMoles;
        surfaceArea = A * DCFrac[j];
        DCChange = k * surfaceArea * pow((SI - 1),gamma);

        ///
        /// Convert time step from days to seconds, because DCChange
        /// has units of mol/s
        ///

        DCChange = DCChange * timestep * 24.0 * 60.0 * 60.0; // DCChange unit: moles/s

        double newDCMoles = chemSys_->getDCMoles(DCIndex[j]) + DCChange;
        chemSys_->setDCUpperLimit(DCIndex[j],newDCMoles);
        chemSys_->setDCLowerLimit(DCIndex[j],newDCMoles);

      }
    }

    return;
}

void KineticModel::getPhaseMasses(void)
{
    int microPhaseId;
    double pscaledMass = 0.0;

    if (microPhaseId != VOIDID && microPhaseId != ELECTROLYTEID) {
        scaledMass_ = chemSys_->getMicroPhaseMass(microPhaseId);
        initScaledMass_ = scaledMass_;

        // Setting the phase mass will also automatically calculate the phase volume
            
        #ifdef DEBUG
            cout << "KineticModel::getPhaseMasses Kinetic model reads solid micphase mass of "
                 << chemSys_->getMicroPhaseName(microPhaseId)
                 << " as " << initScaledMass_ << endl;
            cout.flush();
        #endif
    }

    return;
}

void KineticModel::setKineticDCMoles ()
{

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
            throw FloatException("KineticModel","setKineticDCmoles",
                                 "Divide by zero error");
        }
        #ifdef DEBUG
            cout << "KineticModel::setKineticDCmoles        Clinker phase "
                 << name_ << ": Mass = " << scaledMass_
                 << ", Molar mass = " << chemSys_->getDCMolarMass(DCId)
                 << endl;
        #endif
        chemSys_->setDCMoles(DCId_,(scaledMass_
                             / chemSys_->getDCMolarMass(DCId_)));
    }
    catch (EOBException eex) {
        eex.printException();
        exit(1);
    }
    catch (FloatException fex) {
        fex.printException();
        exit(1);
    }
    catch (out_of_range &oor) {
        EOBException ex("KineticModel","setKineticDCMoles",
                           oor.what(),0,0);
        ex.printException();
        exit(1);
    }
    return;
}

void KineticModel::zeroKineticDCMoles ()
{
    try {
        chemSys_->setDCMoles(DCId_,0.0);
    }
    catch (EOBException eex) {
        eex.printException();
        exit(0);
    }
    return;
}
