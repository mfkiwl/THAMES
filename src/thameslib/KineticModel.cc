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
    name_.clear();
    microPhaseId_.clear();
    DCId_.clear();
    GEMPhaseId_.clear();
    RdICId_.clear();
    Rd_.clear();
    activationEnergy_.clear();
    scaledMass_.clear();
    initScaledMass_.clear();
    isKinetic_.clear();
    
    ///
    /// The default is to not have sulfate attack or leaching, so we set the default
    /// time for initiating these simulations to an absurdly large value: 10 billion
    /// days or 27 million years
    ///

    sulfateAttackTime_ = 1.0e10;
    leachTime_ = 1.0e10;

    return;
}

void KineticModel::getPhaseMasses(void)
{
    int microPhaseId;
    double pscaledMass = 0.0;

    wcRatio_ = lattice_->getWsratio();

    for (int i = 0; i < microPhaseId_.size(); i++) {
        microPhaseId = microPhaseId_[i];
        if (microPhaseId != VOIDID && microPhaseId != ELECTROLYTEID) {
            pscaledMass = chemSys_->getMicroPhaseMass(microPhaseId);
            scaledMass_[i] = pscaledMass;
            initScaledMass_[i] = pscaledMass;

            // Setting the phase mass will also automatically calculate the phase volume
            
            #ifdef DEBUG
                cout << "KineticModel::getPhaseMasses Kinetic model reads solid micphase mass of "
                     << chemSys_->getMicroPhaseName(microPhaseId)
                     << " as " << initScaledMass_[i] << endl;
                cout.flush();
            #endif
        }
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
        for (int i = 0; i < microPhaseId_.size(); i++) {
            if (isKinetic(i)) {
                if (chemSys_->getDCMolarMass(DCId_[i]) <= 0.0) {
                    throw FloatException("KineticModel","setKineticDCmoles",
                                         "Divide by zero error");
                }
                #ifdef DEBUG
                    cout << "KineticModel::setKineticDCmoles        Clinker phase "
                         << name_[i] << ": Mass = " << scaledMass_[i]
                         << ", Molar mass = " << chemSys_->getDCMolarMass(DCId_[i])
                         << endl;
                #endif
                chemSys_->setDCMoles(DCId_[i],(scaledMass_[i]
                                 / chemSys_->getDCMolarMass(DCId_[i])));
            } // END OF IF KINETIC
        }
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
        for (int i = 0; i < microPhaseId_.size(); i++) {
            if (isKinetic(i)) {
                chemSys_->setDCMoles(DCId_[i],0.0);
            }
        }
    }
    catch (EOBException eex) {
        eex.printException();
        exit(0);
    }
    return;
}
