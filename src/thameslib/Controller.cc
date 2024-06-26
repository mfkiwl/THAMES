/**
@file Controller.cc
@brief Definition of Controller class methods
*/

#include "Controller.h"

Controller::Controller(Lattice *msh, KineticController *kc, ChemicalSystem *cs,
                       ThermalStrain *thmstr, const int simtype,
                       const string &parfilename, const string &jobname,
                       const bool verbose, const bool warning)
    : lattice_(msh), kineticController_(kc), chemSys_(cs),
      sim_type_(simtype), thermalstr_(thmstr), jobroot_(jobname) {
  unsigned int i;
  double tvalue, pvalue;
  string buff;
  vector<double> phases;
  const string imgfreqstr = "Image_frequency:";
  const string outtimestr = "OutTime:";
  const string calctimestr = "CalcTime:";

#ifdef DEBUG
  verbose_ = true;
  warning_ = true;
#else
  verbose_ = verbose;
  warning_ = warning;
#endif

  ///
  /// Set default values for all parameters prior to any customization
  ///

  ///
  /// Setting the default times for outputting images, and for initiating the
  /// simulations for leaching or external sulfate attack.
  ///
  /// All times are given in days, and the leaching and sulfate attack times are
  /// set to very high values so that they usually won't happen
  ///

  imgfreq_ = 7.0;
  leach_time_ = 1.0e10;
  sattack_time_ = 1.0e10;

  damagecount_ = 0;

  ///
  /// Load up the pointers to the `ChemicalSystem` object and `Lattice` object
  ///

  chemSys_ = lattice_->getChemSys();
  lattice_->setJobroot(jobroot_);

  ///
  /// Output the class codes for the solution and for DC components.
  /// Output the header for the microstructure phase stats file
  /// Output header for the file tracking pH
  /// Output header for the file tracking the C-S-H composition and Ca/Si ratios
  /// Output header for the file tracking the IC moles in the system
  ///

  try {
    string outfilename = jobroot_ + "_Solution.csv";
    ofstream out(outfilename.c_str());
    if (!out) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    char cc;
    out << "Time(d)";
    for (int i = 0; i < chemSys_->getNumDCs(); i++) {
      cc = chemSys_->getDCClassCode(i);
      if (cc == 'S' || cc == 'T' || cc == 'W') {
        out << "," << chemSys_->getDCName(i);
      }
    }
    out << endl;
    out.close();

    outfilename = jobroot_ + "_DCVolumes.csv";
    ofstream out1(outfilename.c_str());
    if (!out1) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }

    out1 << "Time(d)";
    for (int i = 0; i < chemSys_->getNumDCs(); i++) {
      cc = chemSys_->getDCClassCode(i);
      if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M' || cc == 'W') {
        out1 << "," << chemSys_->getDCName(i);
      }
    }
    out1 << endl;
    out1.close();

    outfilename = jobroot_ + "_Microstructure.csv";
    ofstream out2(outfilename.c_str());
    if (!out2) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }

    out2 << "Time(d)";
    for (int i = 0; i < chemSys_->getNumMicroPhases(); i++) {
      out2 << "," << chemSys_->getMicroPhaseName(i);
    }
    out2 << ",Total Volume (m3),Chemical Shrinkage (m3)";
    out2 << endl;
    out2.close();

    outfilename = jobroot_ + "_pH.csv";
    ofstream out3(outfilename.c_str());
    if (!out3) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    out3 << "Time(d),pH" << endl;
    out3.close();

    outfilename = jobroot_ + "_CSH.csv";
    ofstream out4(outfilename.c_str());
    if (!out4) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    out4 << "Time(d)";
    for (int i = 0; i < chemSys_->getNumICs(); i++) {
      out4 << "," << chemSys_->getICName(i);
    }
    out4 << ",Ca/Si" << endl;
    out4.close();

    outfilename = jobroot_ + "_CSratio_solid.csv";
    ofstream out5(outfilename.c_str());
    if (!out5) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    out5 << "Time(d),C/S in solid" << endl;
    out5.close();

    outfilename = jobroot_ + "_Enthalpy.csv";
    ofstream out6(outfilename.c_str());
    if (!out6) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    out6 << "Time(d),Enthalpy(J)" << endl;
    out6.close();

  } catch (FileException fex) {
    throw fex;
  }

  ///
  /// Open and read the Controller parameter file
  ///

  string xmlext = ".xml";
  size_t foundxml;

  try {
    foundxml = parfilename.find(xmlext);

    time_.clear();
    if (foundxml != string::npos) {
      parseDoc(parfilename);
    } else {
      cout << "Parameter file must be XML" << endl;
      throw FileException("Controller", "Controller", parfilename,
                          "NOT XML FORMAT");
    }
  } catch (FileException fex) {
    throw fex;
  }

}

void Controller::doCycle(const string &statfilename, int choice) {
    unsigned int i;
    int time_index;
    //double next_stat_time = statfreq_;
    RestoreSystem iniLattice;

    ///
    /// This block arbitrarily sets the leaching initiation time to 100 days if
    /// the leaching module is to be run, or sets the sulfate attack initiation
    /// time to 100 days if the sulfate attack module is to be run
    ///
    /// @todo Think about generalizing this more, or allowing combinations of more
    /// than one
    ///

    if (choice == LEACHING) {
        leach_time_ = 100.0;
    } else if (choice == SULFATE_ATTACK) {
        sattack_time_ = 100.0;
    }

    /*
  kineticController_->setSattack_time(sattack_time_);
  kineticController_->setLeach_time(leach_time_);
  */

    chemSys_->setSulfateAttackTime(sattack_time_);
    chemSys_->setLeachTime(leach_time_);
    lattice_->setSattack_time(sattack_time_);
    lattice_->setLeach_time(leach_time_);

    // Initialize the list of all interfaces in the lattice

    if (verbose_) {
        cout << "Controller::doCycle Entering Lattice::findInterfaces" << endl;
    }

    lattice_->findInterfaces();

    if (verbose_) {
        cout << "Controller::Entering Main time loop" << endl;
    }

    double timestep = 0.0;
    bool capwater = true; // True if some capillary water is available
    time_index = 0;

    ///
    /// Output a file that directly links the microstructure ids to their
    /// rgb color.  This is only for easier image processing after the simulation
    /// is finished so we don't have to read the xml file
    ///

    lattice_->writeMicroColors(jobroot_);


    ///
    /// Write the initial microstructure image and its png image
    ///

    lattice_->writeLattice(0.0, sim_type_, jobroot_);
    lattice_->writeLatticePNG(0.0, sim_type_, jobroot_);

    int timesGEMFailed_loc = 0;

    // init to 0 all DC moles corresponding to the kinetic controlled microphases
    //      these DCmoles will be updated by KineticController::calculateKineticStep and
    //      passedd to GEM together the other DC moles in the stystem (ChemicalSystem::calculateState)
    int numMicPh = chemSys_->getNumMicroPhases();
    //cout << "numMicPh : " << numMicPh << endl;


    int DCId;
    for (int i = FIRST_SOLID; i < numMicPh; i++){
        if(chemSys_->isKinetic(i)){
            DCId = chemSys_->getMicroPhaseDCMembers(i, 0);
            //chemSys_->setDCMoles(DCId,0.0); //coment if DCLowerLimit in kineticControllerStep/GEM_from_MT
            chemSys_->setIsDCKinetic(DCId,true);
        }
    }
    cout << endl << "***   numGEMPhases_  = " <<  chemSys_->getNumGEMPhases() << endl;
    cout << "***   numDCs_        = " <<  chemSys_->getNumDCs() << endl;
    cout << "***   numICs_        = " <<  chemSys_->getNumICs() << endl;

    //cout << "Starting with a pore solution without dissolved DCs  => all microPhaseSI_ = 0" << endl;
    //init to 0 all microPhaseSI_
    //chemSys_->setZeroMicroPhaseSI();

    bool writeICsDCs = true;
    if (writeICsDCs) writeTxtOutputFiles_onlyICsDCs(0); // to check the total ICs

    cout << endl << "     ===== START SIMULATION =====" << endl ;

    int cyc;
    for (i = 0; (i < time_.size()) && (capwater); ++i) { //main computation cycle loop

        ///
        /// Do not advance the time step if GEM_run failed the last time
        ///

        bool isFirst = (i == 0) ? true : false;

        if (chemSys_->getTimesGEMFailed() > 0)
            i -= 1;//if(i > 69)DCUpperLimit_[i] = 0;

        //cout << "Time = " << time_[i] << endl;
        cyc = i + 1;
        cout << endl << "Controller::doCycle     i/cyc/Time(i.e. time_[i]) : "
             << i << " / " << cyc << " / " << time_[i] << endl;

        //time_t lt10 = time(NULL);
        //struct tm *time10;
        //time10 = localtime(&lt10);
        //cout << asctime(time10);

        timestep = (i > 0) ? (time_[i] - time_[i - 1]) : (time_[i]);

        ///
        /// Assume that only capillary pore water is chemically reactive,
        /// while water in nanopores is chemically inert.
        ///
        ///
        /// This is the main step of the cycle; the calculateState method
        /// runs all the major steps of a computational cycle
        ///

        try {

            timesGEMFailed_loc = calculateState(time_[i], timestep, isFirst, cyc);

        } catch (GEMException gex) {
            lattice_->writeLattice(time_[i], sim_type_, jobroot_);
            lattice_->writeLatticePNG(time_[i], sim_type_, jobroot_);
            throw gex;
        }

        ///
        /// Once the change in state is determined, propagate the consequences
        /// to the 3D microstructure only if the GEM_run calculation succeeded.
        /// Otherwise we will need to tweak the IC moles so just return from
        /// function without doing anything else
        ///

        if (timesGEMFailed_loc > 0) {
            // Skip the remainder of this iteration and go to the next iteration
            if (warning_) {
                cout << "Controller::doCycle WARNING: Previous call to "
                     << "GEM_run failed, so I will" << endl
                     << "Controller::doCycle not update the microstructure "
                     << "or do anything else" << endl
                     << "controller::doCycle during this time step" << endl;
                cout.flush();
            }
            cout << "i/cyc/time[i]/getTimesGEMFailed: " << i << " / " << cyc << " / "
                 << time_[i] << "\t" << timesGEMFailed_loc << endl;
            continue;
        }

        if (verbose_) {
            cout << "Controller::doCycle Entering Lattice::changeMicrostructure"
                 << endl;
            cout.flush();
        }

        ///
        /// Next function can encounter EOB exceptions within but they
        /// are caught there and the program will then exit from within
        /// this function rather than throwing an exception itself
        ///

        try {
            //set iniLattice i.e. a copy of initial lattice/system configuration (including all DCs values)
            //from ChemicalSystem:
            iniLattice.DCMoles = kineticController_->getDCMoles();

            //from Lattice:
            iniLattice.count = lattice_->getCount();
            iniLattice.site.clear();
            RestoreSite site_l; // only one declaration
            int dimLatticeSite = lattice_->getNumsites(); // only one declaration
            for(int i = 0; i < dimLatticeSite; i++){
                site_l.microPhaseId = (lattice_->getSite(i))->getMicroPhaseId();
                site_l.dissolution = (lattice_->getSite(i))->getDissolutionPhases();
                site_l.growth = (lattice_->getSite(i))->getGrowthPhases();
                site_l.wmc = (lattice_->getSite(i))->getWmc();
                site_l.wmc0 = (lattice_->getSite(i))->getWmc0();
                site_l.visit = 0;
                iniLattice.site.push_back(site_l);
            }
            iniLattice.interface.clear();
            RestoreInterface interface_l; //only one declaration
            int dimLatticeInterface = lattice_->getInterfaceSize();  //only one declaration
            for(int i = 0; i < dimLatticeInterface; i++){
                interface_l.microPhaseId = lattice_->getInterface(i).getMicroPhaseId();
                interface_l.growthSites = lattice_->getInterface(i).getGrowthSites();
                interface_l.dissolutionSites = lattice_->getInterface(i).getDissolutionSites();
                iniLattice.interface.push_back(interface_l);
            }

            /*
            int sizeG, sizeD;
            string nameFileG;
            string nameFileD;
            for(int i = 0; i < dimLatticeInterface; i++){
                nameFileG = "contrGrow_0_" + to_string(i) + ".dat";
                nameFileD = "contrDiss_0_" + to_string(i) + ".dat";
                ofstream outG(nameFileG.c_str());
                ofstream outD(nameFileD.c_str());
                sizeG = iniLattice.interface[i].growthSites.size();
                sizeD = iniLattice.interface[i].dissolutionSites.size();
                outG << "interface_[" << i << "]growthSites.size() = " << sizeG
                     << "\t& microPhaseId_ = " << iniLattice.interface[i].microPhaseId << endl;
                for(int j = 0; j < sizeG; j++){
                    outG << iniLattice.interface[i].growthSites[j].getId() << "\t"
                         << iniLattice.interface[i].growthSites[j].getAffinity() << "\t"
                         << iniLattice.interface[i].growthSites[j].getVerbose() << "\t"
                         << iniLattice.interface[i].growthSites[j].getProb() << "\t"
                         << iniLattice.interface[i].growthSites[j].getProbIni() << endl;
                }
                outG.close();
                outD << "interfacsolidMasse_[" << i << "]dissolutionSites.size() = "
                     << sizeD << "\t& microPhaseId_ = "
                     << iniLattice.interface[i].microPhaseId << endl;
                for(int j = 0; j < sizeD; j++){
                    outD << iniLattice.interface[i].dissolutionSites[j].getId() << "\t"
                         << iniLattice.interface[i].dissolutionSites[j].getAffinity() << "\t"
                         << iniLattice.interface[i].dissolutionSites[j].getVerbose() << "\t"
                         << iniLattice.interface[i].dissolutionSites[j].getProb() << "\t"
                         << iniLattice.interface[i].dissolutionSites[j].getProbIni() << endl;
                }
                outD.close();
            }
            */

            int numDiff = 1000, phDiff = 1000;
            string nameDiff = "testDiff";
            int changeLattice;
            int whileCount = 0;
            changeLattice = lattice_->changeMicrostructure(time_[i], sim_type_, isFirst, capwater, numDiff,
                                                           phDiff, nameDiff, whileCount);

            //if error from changeMicrostructure (not all the voxels given by GEM can
            //be switched to a new DCId):
            //  - comeback to the initial system configuration contained by iniLattice
            //  - re-run GEM with restrictions impossed by the system configuration (DC distribution
            //      on the lattice sites) i.e.
            //      the primal solution must contained a number of moles corresponding to
            //      "numDiff" lattice sites for the microphase "phDiff"
            //
            //restore initial system:
            while (changeLattice == 0){ //} // - for many phases or more loops for the same phase!
            //if (changeLattice == 0) {
                DCId = chemSys_->getMicroPhaseDCMembers(phDiff,0);
                //cout << endl << "phDiff/DCId = " << phDiff << " / " << DCId
                //     << "\tvolMol: " << chemSys_->getDCMolarVolume(phDiff) << " / "
                //     << chemSys_->getDCMolarVolume(DCId) << endl << endl;cout.flush();
                double volMolDiff = chemSys_->getDCMolarVolume(DCId); // m3/mol
                double molarMassDiff = chemSys_->getDCMolarMass(DCId); // g/mol

                double vfracDiff = ((double)numDiff) / ((double)dimLatticeSite);

                double microPhaseMassDiff = vfracDiff * molarMassDiff / volMolDiff / 1.0e6; // g/cm3

                double scaledMassDiff = microPhaseMassDiff * 100.0 / lattice_->getInitSolidMass();

                double numMolesDiff = scaledMassDiff / molarMassDiff;

                //for ChemicalSystem:
                unsigned int numDCs = chemSys_->getNumDCs();
                for(int i = 0; i < numDCs; i++){
                    chemSys_->setDCMoles(i,iniLattice.DCMoles[i]);
                }

                //for Lattice:
                lattice_->setCount(iniLattice.count);

                for(int i = 0; i < dimLatticeSite; i++){
                    (lattice_->getSite(i))->setMicroPhaseId(iniLattice.site[i].microPhaseId);
                    (lattice_->getSite(i))->setDissolutionSite(iniLattice.site[i].dissolution);
                    (lattice_->getSite(i))->setGrowthPhases(iniLattice.site[i].growth);
                    (lattice_->getSite(i))->setWmc(iniLattice.site[i].wmc);
                    (lattice_->getSite(i))->setWmc0(iniLattice.site[i].wmc0);
                    (lattice_->getSite(i))->setVisit(iniLattice.site[i].visit); //or 0!
                }

                for(int i = 0; i < dimLatticeInterface; i++){
                    lattice_->setInterfaceMicroPhaseId(i,iniLattice.interface[i].microPhaseId); // same as before!
                    lattice_->setGrowthSites(i,iniLattice.interface[i].growthSites);
                    lattice_->setDissolutionSites(i,iniLattice.interface[i].dissolutionSites);
                }

                chemSys_->setDCLowerLimit(DCId,numMolesDiff);
                int timesGEMFailed_recall = chemSys_->calculateState(time_[i], isFirst, cyc, true);
                chemSys_->setDCLowerLimit(DCId,0);
                if (timesGEMFailed_recall > 0) {
                    // Skip the remainder of this iteration and go to the next iteration
                    if (warning_) {
                        cout << "Controller::doCycle WARNING: Previous call to "
                             << "GEM_run failed, so I will" << endl
                             << "Controller::doCycle not update the microstructure "
                             << "or do anything else" << endl
                             << "controller::doCycle during this time step" << endl;
                        cout.flush();
                    }
                    cout << "i/time[i]/getTimesGEMFailed_recall: " << i << "\t" << time_[i] << "\t" << timesGEMFailed_recall << endl;
                    continue;
                }

                if (verbose_) {
                    cout << endl << "*Returned from ChemicalSystem::calculateState" << endl;
                    cout << "*called by function Controller::calculateState" << endl;
                    cout << "*timesGEMFailed = " << timesGEMFailed_recall << endl;
                    cout.flush();
                }

                numDiff = 5000, phDiff = 5000;
                nameDiff = "testDiff_recall";

                whileCount++;
                changeLattice = lattice_->changeMicrostructure(time_[i], sim_type_, isFirst, capwater, numDiff, phDiff, nameDiff, whileCount);
                //if(changeLattice != 1 ) {cout << endl << " end changeLattice recall changeLattice = " << changeLattice << endl; exit(0);}
            }


        } catch (DataException dex) {
            lattice_->writeLattice(time_[i], sim_type_, jobroot_);
            lattice_->writeLatticePNG(time_[i], sim_type_, jobroot_);
            throw dex;
        } catch (EOBException ex) {
            lattice_->writeLattice(time_[i], sim_type_, jobroot_);
            lattice_->writeLatticePNG(time_[i], sim_type_, jobroot_);
            throw ex;
        } catch (MicrostructureException mex) {
            lattice_->writeLattice(time_[i], sim_type_, jobroot_);
            lattice_->writeLatticePNG(time_[i], sim_type_, jobroot_);
            throw mex;
        }

        // write output .dat files
        if (writeICsDCs) writeTxtOutputFiles_onlyICsDCs(time_[i]);
        writeTxtOutputFiles(time_[i]);

        ///
        /// Calculate the pore size distribution and saturation
        ///

        lattice_->calculatePoreSizeDistribution();

        ///
        /// Check if there is any capillary pore water remaining.  If not then
        /// we ASSUME hydration has stopped.
        ///
        /// @todo Generalize this idea to allow nanopore water to react by taking
        /// into account its lower chemical potential.

        if (verbose_) {
            cout << "Controller::doCycle Returned from Lattice::changeMicrostructure" << endl;
            cout.flush();
        }

        if ((time_[i] >= output_time_[time_index]) &&
            (time_index < output_time_.size())) {
            if (verbose_) {
                cout << "Controller::doCycle Writing lattice at time_[" << i
                     << "] = " << time_[i] << ", output_time_[" << time_index
                     << "] = " << output_time_[time_index] << endl;
            }
            lattice_->writeLattice(time_[i], sim_type_, jobroot_);
            lattice_->writeLatticePNG(time_[i], sim_type_, jobroot_);
            lattice_->writePoreSizeDistribution(time_[i], sim_type_, jobroot_);

            // lattice_->CheckPoint(jobroot_);
            time_index++;
        }

        double watervolume = chemSys_->getMicroPhaseVolume(ELECTROLYTEID);

        if (watervolume < 2.0e-18) { // Units in m3, so this is about two voxels,
            // we will stop hydration
            if (warning_) {
                cout << "Controller::doCycle WARNING: System is out of capillary pore "
                        "water."
                     << endl;
                cout << "Controller::doCycle          This version of code assumes "
                        "that only capillary"
                     << endl;
                cout << "Controller::doCycle          water is chemically reactive, so "
                        "the system is"
                     << endl;
                cout << "Controller::doCycle          is assumed to be incapable of "
                        "further hydration."
                     << endl;
                cout.flush();
            }
        }

        ///
        /// The following block executes only for sulfate attack simulations
        ///

        if (time_[i] >= sattack_time_) {

            cout << endl << " Controller::doCycle - for sulfate attack, check conditions for addDissolutionSites & coordination sphere " << endl;
            cout << " program stops " << endl;
            exit(1);


            if (verbose_) {
                cout << "Controller::doCycle Sulfate attack module" << endl;
                cout.flush();
            }
            map<int, vector<double>> expansion;
            expansion = lattice_->getExpansion();

            ifstream instopexp("stopexp.dat");
            if (!instopexp) {
                if (verbose_)
                    cout << "keep expanding." << endl;
            } else {
                expansion.clear();
                cout << "expansion has been stopped due to the percolation of damage."
                     << endl;
            }
            cout.flush();

            ///
            /// Stop FM temporarily
            /// @todo What is this?  The following if block will never be run if
            /// uncommented!
            ///

            /*
      expansion.clear();
      */

            if (expansion.size() > 1) {

                damagecount_ = 0;
                double poreintroduce = 0.5;

                if (verbose_) {
                    cout << "Controller::doCycle Sulfate attack module writing " << endl;
                    cout << "Controller::doCycle lattice at time_[" << i
                         << "] = " << time_[i] << ", " << endl;
                    cout << "controller::doCycle output_time_[" << time_index
                         << "] = " << output_time_[time_index] << endl;
                    cout.flush();
                }

                lattice_->writeLattice(time_[i], sim_type_, jobroot_);
                lattice_->writeLatticePNG(time_[i], sim_type_, jobroot_);
                string ofileName(jobroot_);
                ostringstream ostr1, ostr2;
                ostr1 << (int)(time_[i] * 100);
                ostr2 << setprecision(3) << chemSys_->getTemperature();
                string timestr(ostr1.str());
                string tempstr(ostr2.str());
                ofileName = ofileName + "." + timestr + "." + tempstr + ".img";

                ///
                /// In the sulfate attack algorithm, calculate the stress and strain
                /// distributions
                ///

                thermalstr_->setEigen();
                for (map<int, vector<double>>::iterator it = expansion.begin();
                     it != expansion.end(); it++) {

                    int expindex = it->first;
                    vector<double> expanval = it->second;
                    vector<int> expcoordin = lattice_->getExpansionCoordin(expindex);
                    thermalstr_->setEigen(expindex, expanval[0], expanval[1], expanval[2],
                                          0.0, 0.0, 0.0);
                    thermalstr_->setExp(expindex, expcoordin);

                    ///
                    /// Set expansion site to be damaged if there is one, as determined by
                    /// the setEigen function returning every site above damage stress
                    /// threshold
                    ///
                    Site *ste;
                    ste = lattice_->getSite(expindex);
                    /*
          lattice_->dWaterchange(poreintroduce);
          */

                    double dwmcval = poreintroduce;
                    lattice_->dWmc(expindex, dwmcval);
                    for (int j = 0; j < ste->nbSize(2); j++) {
                        Site *stenb = ste->nb(j);
                        stenb->dWmc(dwmcval);
                    }
                }

                ///
                /// Calculate the stress-free strain (thermal strain) state in the
                /// microstructure, and then write the displacement field
                ///

                thermalstr_->Calc(time_[i], ofileName, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

                // thermalstr_ -> writeStress(jobroot_,time_[i],0); //write strxx
                // thermalstr_ -> writeStrainEngy(jobroot_,time_[i]);
                thermalstr_->writeDisp(jobroot_, time_[i]);

                ///
                /// Get the true volume of each voxel after FEM calculation
                ///

                double otruevolume = 0.0;
                double truevolume = 0.0;
                for (int ii = 0; ii < lattice_->getNumsites(); ii++) {
                    Site *ste;
                    ste = lattice_->getSite(i);
                    otruevolume = ste->getTrueVolume();
                    for (int j = 0; j < 3; j++) {
                        truevolume += thermalstr_->getEleStrain(ii, j);
                    }
                    truevolume = otruevolume * (1 + truevolume);
                    ste->setTrueVolume(truevolume);
                }

                for (int index = 0; index < lattice_->getNumsites(); index++) {
                    Site *ste;
                    ste = lattice_->getSite(index);
                    int pid = ste->getMicroPhaseId();

                    if (ste->IsDamage()) {
                        damagecount_++;
                    }

                    if ((ste->IsDamage())) {
                        double strxx, stryy, strzz;
                        strxx = stryy = strzz = 0.0;
                        strxx = thermalstr_->getEleStress(index, 0);
                        stryy = thermalstr_->getEleStress(index, 1);
                        strzz = thermalstr_->getEleStress(index, 2);
                        if ((strxx >= 1.0) || (stryy >= 1.0) || (strzz >= 1.0)) {
                            vector<double> damageexp;
                            damageexp.clear();
                            double poreincrease = 0.2;
                            damageexp.resize(3, (1.0 / 3.0 * poreincrease));
                            vector<double> damageexpo;
                            damageexpo.clear();
                            damageexpo = lattice_->getExpansion(index);
                            for (int i = 0; i < 3; i++) {
                                damageexp[i] += damageexpo[i];
                            }
                            lattice_->setExpansion(index, damageexp);
                            lattice_->dWaterchange(poreincrease);
                            /// JWB: This next line must be from some earlier version
                            /// The called method does not exist any longer
                            ///
                            /// @todo: Determine whether it is necessary to add this back
                            /// in for crystallization pressure calculations
                            //
                            // ste->setVolume(VOIDID,(ste->getVolume(VOIDID) + poreincrease));
                            //
                        }
                    }

                    ///
                    /// The next block gets the stress in each voxel that does NOT
                    /// contain a clinker phase (C3S, C2S, C3A, or C4AF), then determine
                    /// if the voxel should be damaged as a result

                    /// Prefer to make this independent of whether or not there is C4AF in
                    /// the phase definitions.  What if this is a white cement or
                    /// something?
                    ///
                    /// @note Associating the last clinker phase with id 5 is a kluge
                    /// @todo Give each phase a calcstress property or something like that
                    ///       that can be checked instead of hardwiring phase ids

                    if (pid > 5) {
                        double strxx, stryy, strzz;
                        strxx = stryy = strzz = 0.0;
                        strxx = thermalstr_->getEleStress(index, 0);
                        stryy = thermalstr_->getEleStress(index, 1);
                        strzz = thermalstr_->getEleStress(index, 2);
                        if ((strxx >= thermalstr_->getTstrength(index)) ||
                            (stryy >= thermalstr_->getTstrength(index)) ||
                            (strzz >= thermalstr_->getTstrength(index))) {
                            // if (verbose_) cout << "Phase " << pid << " is damaged." <<
                            // endl;
                            if (!ste->IsDamage()) {
                                // if (verbose_) cout << " it has not been damaged before." <<
                                // endl;
                                ste->setDamage();
                                damagecount_++;
                                // lattice_->dWaterchange(poreintroduce);

                                double dwmcval = poreintroduce;
                                lattice_->dWmc(index, dwmcval);
                                for (int j = 0; j < ste->nbSize(1); j++) { //NN_NNN?
                                    Site *stenb = ste->nb(j);
                                    stenb->dWmc(dwmcval);
                                    if ((stenb->getWmc() > 0.0) &&
                                        (stenb->getMicroPhaseId() != ELECTROLYTEID) &&
                                        (stenb->getMicroPhaseId() != VOIDID)) {
                                        lattice_->addDissolutionSite(stenb,
                                                                        stenb->getMicroPhaseId());
                                    }
                                }
                                for (int j = ste->nbSize(1); j < ste->nbSize(2); j++) {
                                    Site *stenb = ste->nb(j);
                                    stenb->dWmc(dwmcval);
                                }

                                /*
                vector<double> damageexp;
                damageexp.clear();
                double poreindamage = 0.6;
                damageexp.resize(3,(1.0 / 3.0 * poreindamage));
                lattice_->setExpansion(index,damageexp);
                vector<int> coordin;
                coordin.clear();
                coordin.resize(3,0);
                coordin[0] = ste->getX();
                coordin[1] = ste->getY();
                coordin[2] = ste->getZ();
                lattice_->setExpansionCoordin(index,coordin);
                lattice_->dWaterchange(poreindamage);
                ste->setVolume(VOIDID,poreindamage);
                */
                            }
                        }
                    }

                } // End of loop over all voxels

                if (verbose_) {
                    cout << "Controller::doCycle sulfate attack module Time = "
                         << time_[i] << " damagecount_ is: " << damagecount_ << endl;
                    cout.flush();
                }
                ofstream outdamage("damage.dat");
                outdamage << damagecount_;
                outdamage.close();

                string damagejobroot = jobroot_ + ".damage";
                lattice_->writeDamageLattice(time_[i], damagejobroot);
                lattice_->writeDamageLatticePNG(time_[i], damagejobroot);
                // to see whether new damage is generated
            }
        }
    }

    ///
    /// Write the final lattice state to an ASCII file and to a PNG file for
    /// visualization
    ///

    lattice_->writeLattice(time_[i - 1], sim_type_, jobroot_);
    lattice_->writeLatticePNG(time_[i - 1], sim_type_, jobroot_);

    return;
}

int Controller::calculateState(double time, double dt, bool isFirst, int cyc) {

    int timesGEMFailed = 0;

    try {

        if (isFirst) {

            double T = chemSys_->getTemperature();
            lattice_->setTemperature(T);
        }

        ///
        /// We must pass some vectors to the `calculateKineticStep` method that
        /// will hold the amounts of impurity elements released from the clinker
        /// phases.  These values do not go in the `ChemicalSystem` object, but will
        /// still need to be processed afterward.
        ///

        //vector<double> impurityrelease;
        //impurityrelease.clear();
        //impurityrelease.resize(chemSys_->getNumMicroImpurities(), 0.0);

        ///
        /// Get the number of moles of each IC dissolved from kinetically controlled
        /// phases
        ///

        double T = lattice_->getTemperature();


        /// The thermodynamic calculation returns the saturation index of phases,
        /// which is needed for calculations of driving force for dissolution
        /// or growth.
        ///
        /// 2024-05-29:  At the moment, if a microstructure phase is defined
        /// to be one or more GEM phases, the SI of the microstructure phase
        /// is calculated as the mole-weighted average of the SIs of the
        /// constituent GEM CSD phases.  Can't think of a better way to do this
        /// except to prohibit users from defining mixtures of CSD phases
        /// as microstructure phases.

        chemSys_->setMicroPhaseSI(cyc);

        kineticController_->calculateKineticStep(dt, T, cyc);
        //kineticController_->calculateKineticEvents(dt, T, isFirst);

        //if (time >= sattack_time_) {// for sulfate attack iterations}

        ///
        /// Now that the method is done determining the change in moles of each IC,
        /// launch a thermodynamic calculation to determine new equilibrium state
        ///
        /// The `ChemicalSystem` object provides an interface for these calculations
        ///

        try {
            timesGEMFailed = chemSys_->calculateState(time, isFirst, cyc, false);
            if (verbose_) {
                cout << "*Returned from ChemicalSystem::calculateState" << endl;
                cout << "*called by function Controller::calculateState" << endl;
                cout << "*timesGEMFailed = " << timesGEMFailed << endl;
                cout.flush();
            }
        } catch (GEMException gex) {
            gex.printException();
            exit(1);
        }

        if (verbose_) {
            cout << "Done!" << endl;
            cout.flush();
        }

        ///
        /// The thermodynamic calculation returns the saturation index of phases,
        /// which is needed for calculations of driving force for dissolution
        /// or growth.  Assign this to the lattice in case crystallization pressures
        /// should be calculated.
        ///

        if (timesGEMFailed > 0) return timesGEMFailed;

        /*
        - all these must be done after the lattice update -
        try {
            double aveSI = 0.0;
            double moles = 0.0;
            vector<int> microPhaseMembers;
            for (int i = 0; i < chemSys_->getNumMicroPhases(); ++i) {
                string pname = chemSys_->getMicroPhaseName(i);
                int newMicroPhaseId =
                    chemSys_->getMicroPhaseId(chemSys_->getMicroPhaseName(i));
                aveSI = moles = 0.0;
                microPhaseMembers = chemSys_->getMicroPhaseMembers(newMicroPhaseId);
                for (int ii = 0; ii < microPhaseMembers.size(); ++ii) {
                    int newGEMPhaseId = microPhaseMembers[ii];
                    aveSI += ((solut_->getSI(newGEMPhaseId) *
                               chemSys_->getGEMPhaseMoles(newGEMPhaseId)));
                    moles += chemSys_->getGEMPhaseMoles(newGEMPhaseId);
                }
                if (moles > 0.0) {
                    aveSI = aveSI / moles;
                } else {
                    aveSI = aveSI / (static_cast<double>(microPhaseMembers.size()));
                }
                lattice_->setSI(newMicroPhaseId, aveSI);
            }
        } catch (EOBException eex) {
            eex.printException();
            exit(1);
        }

        if (isFirst) {
            chemSys_->writeChemSys();
        }
        */

        // Output to files the solution composition data, phase data, DC data,
        // microstructure data, pH, and C-S-H composition and Ca/Si ratio
        // DONE AFTER LATTICE UPDATE

    } catch (FileException fex) {
        fex.printException();
        exit(1);
    } catch (FloatException flex) {
        flex.printException();
        exit(1);
    }
    return timesGEMFailed;
}

void Controller::writeTxtOutputFiles (double time){
    ///
    /// Set the kinetic DC moles.  This adds the clinker components to the DC
    /// moles.
    ///
    /// @todo Find out what this is and why it needs to be done
    ///
    ///
    int numICs = chemSys_->getNumICs();
    int numDCs = chemSys_->getNumDCs();
    int i, j;

    // Output to files the solution composition data, phase data, DC data,
    // microstructure data, pH, and C-S-H composition and Ca/Si ratio

    string outfilename = jobroot_ + "_Solution.csv";
    ofstream out3(outfilename.c_str(), ios::app);
    if (!out3) {
        throw FileException("Controller", "calculateState", outfilename,
                            "Could not append");
    }

    out3 << setprecision(5) << time;
    char cc;
    for (i = 0; i < numDCs; i++) {
        cc = chemSys_->getDCClassCode(i);
        if (cc == 'S' || cc == 'T' || cc == 'W') {
            out3 << "," << (chemSys_->getNode())->Get_cDC((long int)i); // molality
        }
    }
    out3 << endl;
    out3.close();

    outfilename = jobroot_ + "_DCVolumes.csv";
    ofstream out4(outfilename.c_str(), ios::app);
    if (!out4) {
        throw FileException("Controller", "calculateState", outfilename,
                            "Could not append");
    }

    out4 << setprecision(5) << time;
    for (i = 0; i < numDCs; i++) {
        if (chemSys_->getDCMolarMass(i) > 0.0) {
            cc = chemSys_->getDCClassCode(i);
            if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M' || cc == 'W') {
                string dcname = chemSys_->getDCName(i);
                double V0 =
                    chemSys_->getDCMoles(dcname) * chemSys_->getDCMolarVolume(dcname);
                out4 << "," << V0;
                if (verbose_) {
                    cout << "Controller::calculateState    DC = "
                         << chemSys_->getDCName(i)
                         << ", moles = " << chemSys_->getDCMoles(i)
                         << ", molar mass = " << chemSys_->getDCMolarMass(i) << endl;
                    cout.flush();
                }
            }
        } else {
            string msg = "Divide by zero error for DC " + chemSys_->getDCName(i);
            out4.close();
            throw FloatException("Controller", "calculateState", msg);
        }
    }
    out4 << endl;
    out4.close();

    outfilename = jobroot_ + "_Microstructure.csv";
    ofstream out5(outfilename.c_str(), ios::app);
    if (!out5) {
        throw FileException("Controller", "calculateState", outfilename,
                            "Could not append");
    }

    out5 << setprecision(5) << time;
    for (i = 0; i < chemSys_->getNumMicroPhases(); i++) {
        out5 << "," << (lattice_->getVolumefraction(i));
    }
    double micvol = lattice_->getMicrostructurevolume();
    double initmicvol = lattice_->getInitialmicrostructurevolume();
    out5 << "," << micvol << "," << (initmicvol - micvol);
    out5 << endl;
    out5.close();

    outfilename = jobroot_ + "_pH.csv";
    ofstream out6(outfilename.c_str(), ios::app);
    if (!out6) {
        throw FileException("Controller", "calculateState", outfilename,
                            "Could not append");
    }
    out6 << setprecision(5) << time;
    out6 << "," << (chemSys_->getPH()) << endl;
    out6.close();

    chemSys_->setGEMPhaseStoich();

    double *CSHcomp;
    try {
        CSHcomp =
            chemSys_->getPGEMPhaseStoich(chemSys_->getGEMPhaseId(CSHGEMName));
    } catch (EOBException eex) {
        eex.printException();
        exit(1);
    }
    if (verbose_) {
        cout << "Done!" << endl;
        cout.flush();
    }
    //double CaMoles = 0.0, SiMoles = 0.0, CaSiRatio = 0.0;
    outfilename = jobroot_ + "_CSH.csv";
    ofstream out7(outfilename.c_str(), ios::app);
    if (!out7) {
        throw FileException("Controller", "calculateState", outfilename,
                            "Could not append");
    }
    out7 << setprecision(5) << time;
    for (i = 0; i < numICs; i++) {
        out7 << "," << CSHcomp[i];
        //if (chemSys_->getICName(i) == "Ca") {
        //    CaMoles = CSHcomp[i];
        //}
        //if (chemSys_->getICName(i) == "Si") {
        //    SiMoles = CSHcomp[i];
        //}
    }
    int id_Ca = chemSys_->getICId("Ca");
    int id_Si = chemSys_->getICId("Si");
    double CaMoles = CSHcomp[id_Ca];
    double SiMoles = CSHcomp[id_Si];
    if (CaMoles < 1.0e-16)
        CaMoles = 1.0e-16;
    if (SiMoles < 1.0e-16)
        SiMoles = 1.0e-16;
    double CaSiRatio = CaMoles / SiMoles;
    out7 << "," << CaSiRatio << endl;
    out7.close();

    //chemSys_->setGEMPhaseStoich();
    double *phaseRecord;
    int ICIndex;
    CaMoles = SiMoles = 0.0;
    outfilename = jobroot_ + "_CSratio_solid.csv";
    ofstream out8(outfilename.c_str(), ios::app);
    if (!out8) {
        throw FileException("Controller", "calculateState", outfilename,
                            "Could not append");
    }
    out8 << setprecision(5) << time;
    for (i = 0; i < chemSys_->getNumGEMPhases(); i++) {
        cc = chemSys_->getGEMPhaseClassCode(i);
        if (cc == 's') {
            phaseRecord = chemSys_->getPGEMPhaseStoich(i);
            //ICIndex = chemSys_->getICId("Ca");
            //CaMoles += phaseRecord[ICIndex];
            //ICIndex = chemSys_->getICId("Si");
            //SiMoles += phaseRecord[ICIndex];
            CaMoles += phaseRecord[id_Ca];
            SiMoles += phaseRecord[id_Si];
        }
    }
    if (SiMoles != 0) {
        CaSiRatio = CaMoles / SiMoles;
        out8 << "," << CaSiRatio << endl;
    } else {
        out8 << ",Si_moles is ZERO" << endl;
    }
    out8.close();

    outfilename = jobroot_ + "_Enthalpy.csv";
    ofstream out10(outfilename.c_str(), ios::app);
    if (!out10) {
        throw FileException("Controller", "calculateState", outfilename,
                            "Could not append");
    }
    if (verbose_) {
        cout << "Writing Enthalpy values...";
        cout.flush();
    }

    double enth = 0.0;
    for (i = 0; i < numDCs; i++) {
        enth += (chemSys_->getDCEnthalpy(i));
    }

    out10 << setprecision(5) << time;
    out10 << "," << enth << endl;
    if (verbose_) {
        cout << "Done!" << endl;
        cout.flush();
    }
    out10.close();

}

void Controller::writeTxtOutputFiles_onlyICsDCs(double time){

    int i, j;
    int numICs = chemSys_->getNumICs();
    int numDCs = chemSys_->getNumDCs();
    vector<double> ICMoles;
    ICMoles.resize(numICs,0.0);
    vector<double> DCMoles;
    DCMoles.resize(numDCs,0.0);
    for (i = 0; i < numDCs; i++){
        DCMoles[i] = chemSys_->getDCMoles(i);
    }

    string outfilenameIC = jobroot_ + "_icmoles.csv";
    string outfilenameDC = jobroot_ + "_dcmoles.csv";
    if(time < 1.e-10){
        ofstream out0IC(outfilenameIC.c_str());
        out0IC << "Time(d)";
        for (i = 0; i < numICs; i++) {
            out0IC << "," << chemSys_->getICName(i);
        }
        out0IC << endl;
        out0IC.close();

        ofstream out0DC(outfilenameDC.c_str());
        out0DC << "Time(d)";
        for (i = 0; i < numDCs; i++) {
            out0DC << "," << chemSys_->getDCName(i);
        }
        out0DC << endl;
        out0DC.close();
    }   

    vector <int> impurityDCID;
    impurityDCID.clear();
    impurityDCID.push_back(chemSys_->getDCId("K2O"));
    impurityDCID.push_back(chemSys_->getDCId("Na2O"));
    impurityDCID.push_back(chemSys_->getDCId("Per"));
    impurityDCID.push_back(chemSys_->getDCId("SO3"));

    double scMass, molMass;
    int mPhId;
    double massImpurity, totMassImpurity;

    //cout << endl << "getIsDCKinetic: " << endl;
    for(j = 0; j < numDCs; j++){
        if(chemSys_->getIsDCKinetic(j)){
            //molMass = chemSys_->getDCMolarMass(j);
            mPhId = chemSys_->getDC_to_MPhID(j);
            scMass = chemSys_->getMicroPhaseMass(mPhId);

            totMassImpurity = 0;

            massImpurity = scMass * chemSys_->getK2o(mPhId);
            totMassImpurity += massImpurity;
            DCMoles[impurityDCID[0]] += massImpurity / chemSys_->getDCMolarMass("K2O");

            massImpurity = scMass * chemSys_->getNa2o(mPhId);
            totMassImpurity += massImpurity;
            DCMoles[impurityDCID[1]] += massImpurity / chemSys_->getDCMolarMass("Na2O");

            massImpurity = scMass * chemSys_->getMgo(mPhId);
            totMassImpurity += massImpurity;
            DCMoles[impurityDCID[2]] += massImpurity / chemSys_->getDCMolarMass("Per");//MgO

            massImpurity = scMass * chemSys_->getSo3(mPhId);
            totMassImpurity += massImpurity;
            DCMoles[impurityDCID[3]] += massImpurity / chemSys_->getDCMolarMass("SO3");

            DCMoles[j] = (scMass - totMassImpurity) / chemSys_->getDCMolarMass(j);

            //DCMoles[j] = scMass / molMass;
            //DCMoles[impurityDCID[0]] += scMass * chemSys_->getK2o(mPhId) / chemSys_->getDCMolarMass("K2O");
            //DCMoles[impurityDCID[1]] += scMass * chemSys_->getNa2o(mPhId) / chemSys_->getDCMolarMass("Na2O");
            //DCMoles[impurityDCID[2]] += scMass * chemSys_->getMgo(mPhId) / chemSys_->getDCMolarMass("Per");//MgO
            //DCMoles[impurityDCID[3]] += scMass * chemSys_->getSo3(mPhId) / chemSys_->getDCMolarMass("SO3");

            //cout << j << "  scMass : " << scMass << "  mPhId :" << mPhId << "  molMass : " << molMass << endl;
            //cout << "      impurityDCID [0/1/2/3]   : " << impurityDCID[0] << " / " << impurityDCID[1] << endl; cout.flush();
            //                                    //  " / " << impurityDCID[2] << " / " << impurityDCID[3] << endl; cout.flush();
            //cout << "      mass(%) K2O/Na2O/MgO/SO3 : " << chemSys_->getK2o(mPhId) << " / " << chemSys_->getNa2o(mPhId) << endl; cout.flush();
            //                                    //  " / " << chemSys_->getMgo(mPhId) << " / " << chemSys_->getSo3(mPhId) << endl; cout.flush();
            //cout << "      molMass K2O/Na2O/MgO/SO3 : " << chemSys_->getDCMolarMass("K2O") << " / " << chemSys_->getDCMolarMass("Na2O") << endl; cout.flush();
            //                                    //  " / " << chemSys_->getDCMolarMass("MgO") << " / " << chemSys_->getDCMolarMass("SO3") << endl; cout.flush();
        }
    }

    for(j = 0; j < numDCs; j++){
        for(i = 0; i < numICs; i++){
            ICMoles[i] += DCMoles[j]* chemSys_->getDCStoich(j,i);
        }
    }

    ofstream out1(outfilenameIC.c_str(), ios::app);
    if (!out1) {
        throw FileException("Controller", "writeTxtOutputFiles_onlyICsDCs", outfilenameIC,
                            "Could not append");
    }

    out1 << setprecision(5) << time;
    for (i = 0; i < numICs; i++) {
        out1 << "," << ICMoles[i];
    }
    out1 << endl;
    out1.close();

    ofstream out2(outfilenameDC.c_str(), ios::app);
    if (!out2) {
        throw FileException("Controller", "writeTxtOutputFiles_onlyICsDCs", outfilenameDC,
                            "Could not append");
    }

    out2 << setprecision(5) << time;
    for (i = 0; i < numDCs; i++) {
        out2 << "," << DCMoles[i];
    }
    out2 << endl;
    out2.close();

}




void Controller::parseDoc(const string &docName) {
  xmlDocPtr doc;
  xmlNodePtr cur;
  xmlChar *key;
  cout.flush();
  double testtime;
  doc = xmlParseFile(docName.c_str());

  // check if the xml file is valid
  /// @note This block requires the schema file to be local

  try {
    string paramxsd = "parameters.xsd";
    if (!is_xml_valid(doc, paramxsd.c_str())) {
      cout << "XML file is NOT valid" << endl;
      throw FileException("Controller", "parseDoc", docName, "XML not valid");
    } else {
      int boolval;
      if (doc == NULL) {
        throw FileException("Controller", "parseDoc", docName,
                            "XML not parsed successfully");
      }

      cur = xmlDocGetRootElement(doc);

      if (cur == NULL) {
        throw FileException("Controller", "parseDoc", docName,
                            "Empty parameter XML document");
      }

      cur = cur->xmlChildrenNode;
      while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"calctime"))) {
          key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
          string st((char *)key);
          from_string(testtime, st);
          time_.push_back(testtime);
          xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"outtime"))) {
          key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
          string st((char *)key);
          from_string(testtime, st);
          output_time_.push_back(testtime);
          xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"image_frequency"))) {
          key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
          string st((char *)key);
          from_string(imgfreq_, st);
          xmlFree(key);
        }

        cur = cur->next;
      }
    }
    if (doc != NULL)
      xmlFreeDoc(doc);
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }
  return;
}
