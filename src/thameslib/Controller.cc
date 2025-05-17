/**
@file Controller.cc
@brief Definition of Controller class methods
*/

#include "Controller.h"

Controller::Controller(Lattice *msh, KineticController *kc, ChemicalSystem *cs,
                       ThermalStrain *thmstr, const int simtype,
                       const string &jsonFileName, const string &jobname,
                       const bool verbose, const bool warning, const bool xyz)
    : lattice_(msh), kineticController_(kc), chemSys_(cs), simType_(simtype),
      thermalstr_(thmstr), jobRoot_(jobname) {

  xyz_ = xyz;

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
  /// All times are given in hours, and the leaching and sulfate attack times
  /// are set to very high values so that they usually won't happen
  ///

  imgFreq_ = 168.0; // hours = 7 days
  leachTime_ = 1.0e10;
  sulfateAttackTime_ = 1.0e10;

  damageCount_ = 0;

  ///
  /// Load up the pointers to the `ChemicalSystem` object and `Lattice` object
  ///

  chemSys_ = lattice_->getChemSys();
  lattice_->setJobRoot(jobRoot_);

  ///
  /// Output the class codes for the solution and for DC components.
  /// Output the header for the microstructure phase stats file
  /// Output header for the file tracking pH
  /// Output header for the file tracking the C-S-H composition and Ca/Si ratios
  /// Output header for the file tracking the IC moles in the system
  ///

  ofstream outfs;
  string outfilename;

  try {
    outfilename = jobRoot_ + "_Solution.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    char cc;
    outfs << "Time(h)";
    for (int i = 0; i < chemSys_->getNumDCs(); i++) {
      cc = chemSys_->getDCClassCode(i);
      if (cc == 'S' || cc == 'T' || cc == 'W') {
        outfs << "," << chemSys_->getDCName(i);
      }
    }
    outfs << endl;
    outfs.close();

    outfilename = jobRoot_ + "_DCVolumes.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }

    outfs << "Time(h)";
    for (int i = 0; i < chemSys_->getNumDCs(); i++) {
      cc = chemSys_->getDCClassCode(i);
      if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M' || cc == 'W') {
        outfs << "," << chemSys_->getDCName(i) << "(m3/100g)";
      }
    }
    outfs << endl;
    outfs.close();

    outfilename = jobRoot_ + "_SurfaceAreas.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }

    outfs << "Time(h)";
    // JWB: Start with micro phase id 1 to avoid the Void phase
    for (int i = 1; i < chemSys_->getNumMicroPhases(); i++) {
      int dcid = chemSys_->getMicroPhaseDCMembers(i, 0);
      cc = chemSys_->getDCClassCode(dcid);
      if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M') {
        outfs << ",A_" << chemSys_->getMicroPhaseName(i) << "(m2/100g)";
      }
    }
    outfs << endl;
    outfs.close();

    outfilename = jobRoot_ + "_SI.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }

    outfs << "Time(h)";
    // JWB: Start with micro phase id 1 to avoid the Void phase
    for (int i = 1; i < chemSys_->getNumMicroPhases(); i++) {
      int dcid = chemSys_->getMicroPhaseDCMembers(i, 0);
      cc = chemSys_->getDCClassCode(dcid);
      if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M') {
        outfs << ",SI_" << chemSys_->getMicroPhaseName(i);
      }
    }
    outfs << endl;
    outfs.close();

    outfilename = jobRoot_ + "_CSH.csv";
    outfs.open(outfilename.c_str(), ios::app);
    if (!outfs) {
      throw FileException("Controller", "calculateState", outfilename,
                          "Could not append");
    }
    outfs << "Time(h)";
    for (int i = 0; i < numICs_; i++) {
      outfs << "," << chemSys_->getICName(i);
    }
    outfs << ",Ca/Si" << endl;
    outfs.close();

    outfilename = jobRoot_ + "_CSratio_solid.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "calculateState", outfilename,
                          "Could not append");
    }

    outfs << "Time(h),Ca/Si Ratio" << endl;

    outfs.close();

    outfilename = jobRoot_ + "_Microstructure.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }

    outfs << "Time(h)";
    for (int i = 0; i < chemSys_->getNumMicroPhases(); i++) {
      outfs << "," << chemSys_->getMicroPhaseName(i);
    }
    outfs << ",Total Volume (m3/100g),Chemical Shrinkage (m3/100g)";
    outfs << endl;
    outfs.close();

    outfilename = jobRoot_ + "_pH.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    outfs << "Time(h),pH" << endl;
    outfs.close();

    outfilename = jobRoot_ + "_Enthalpy.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    outfs << "Time(h),Enthalpy(J/100g)" << endl;
    outfs.close();

  } catch (FileException fex) {
    throw fex;
  }

  numDCs_ = chemSys_->getNumDCs();
  numICs_ = chemSys_->getNumICs();
  numMicroPhases_ = chemSys_->getNumMicroPhases();
  numGEMPhases_ = chemSys_->getNumGEMPhases();

  cout << endl << "***   numGEMPhases_   = " << numGEMPhases_ << endl;
  cout << "***   numMicroPhases_ = " << numMicroPhases_ << endl;
  cout << "***   numDCs_         = " << numDCs_ << endl;
  cout << "***   numICs_         = " << numICs_ << endl;

  temperature_ = chemSys_->getTemperature();
  lattice_->setTemperature(temperature_);
  presure_ = chemSys_->getP();
  waterDCId_ = chemSys_->getDCId("H2O@");
  waterMolarMass_ = chemSys_->getDCMolarMass("H2O@");
  numSites_ = lattice_->getNumSites();
  initMicroVolume_ = chemSys_->getInitMicroVolume();

  ///
  /// Open and read the Controller parameter file
  ///

  string jsonext = ".json";
  size_t foundjson;

  try {
    foundjson = jsonFileName.find(jsonext);

    if (foundjson != string::npos) {
      parseDoc(jsonFileName);
    } else {
      cout << "Parameter file must be JSON" << endl;
      throw FileException("Controller", "Controller", jsonFileName,
                          "NOT JSON FORMAT");
    }
  } catch (FileException fex) {
    throw fex;
  }

  for (int i = 0; i < time_.size() - 1; i++) {
    if (abs(time_[i] - time_[i + 1]) <= 1.0e-6) {
      time_.erase(time_.begin() + i);
    }
  }

  int time_Size = time_.size();
  int outputTime_Size = outputTime_.size();

  cout << endl << "Controller::Controller(...) :" << endl;

  if (time_.size() == 0) {
    cout << endl
         << endl
         << "Controller::Controller error : a final time and "
         << "at least one output time value must be present in "
         << "the simulation parameters file!" << endl;
    cout << endl
         << "check and modify the simulation parameters file and run thames "
            "again"
         << endl;
    cout << endl << "end program" << endl;
    // exit(0);
    throw FileException(
        "Controller", "Controller", "simparams.json",
        "a final time and at least one "
        "outtime value must be present into simparams.json file!");
  }

  cout << "   final time_.size()          = " << time_Size << endl;

  outfilename = jobRoot_ + "-times_used.json";
  outfs.open(outfilename.c_str());
  if (!outfs) {
    throw FileException("Controller", "Controller", outfilename,
                        "Could not append");
  }

  outfs << "{" << endl;
  outfs << "  \"time_parameters\": {" << endl;
  outfs << "    \"calctimes\": [" << endl;
  // outfs << "        ";
  int j = 0;
  for (int i = 0; i < time_Size; i++) {
    j++;
    if (i < (time_Size - 1)) {
      if (j == 1) {
        outfs << "        " << fixed << time_[i] << ", ";
      } else if (j < 7) {
        outfs << time_[i] << ", ";
      } else {
        j = 0;
        outfs << time_[i] << "," << endl;
      }
    } else {
      if (j == 1) {
        outfs << "        " << time_[i] << endl;
      } else {
        outfs << time_[i] << endl;
      }
    }
  }
  outfs << "    ]," << endl;
  outfs << "    \"outtimes\": [" << endl;
  j = 0;
  for (int i = 0; i < outputTime_Size; i++) {
    j++;
    if (i < (outputTime_Size - 1)) {
      if (j == 1) {
        outfs << "        " << outputTime_[i] << ", ";
      } else if (j < 7) {
        outfs << outputTime_[i] << ", ";
      } else {
        j = 0;
        outfs << outputTime_[i] << "," << endl;
      }
    } else {
      outfs << outputTime_[i] << endl;
    }
  }
  outfs << "    ]" << endl;
  outfs << "  }" << endl;
  outfs << "}" << endl;
  outfs.close();

  cout << "   => new time values (calctimes & outtimes) have been used and "
          "saved in the file :"
       << endl;
  cout << "         " << outfilename << endl;
}

void Controller::doCycle(double elemTimeInterval) {
  unsigned int i;
  int time_index;
  RestoreSystem iniLattice;

  ///
  /// This block arbitrarily sets the leaching initiation time to 100 days if
  /// the leaching module is to be run, or sets the sulfate attack initiation
  /// time to 100 days if the sulfate attack module is to be run
  ///
  /// @todo Think about generalizing this more, or allowing combinations of more
  /// than one
  ///

  if (simType_ == LEACHING) {
    leachTime_ = 2400.0; // 2400 hours = 100 days
  } else if (simType_ == SULFATE_ATTACK) {
    sulfateAttackTime_ = 2400.0; // 2400 hours = 100 days
  }

  /*
  kineticController_->setSattack_time(sulfateAttackTime_);
  kineticController_->setLeach_time(leachTime_);
  */

  chemSys_->setSulfateAttackTime(sulfateAttackTime_);
  chemSys_->setLeachTime(leachTime_);
  lattice_->setSulfateAttackTime(sulfateAttackTime_);
  lattice_->setLeachTime(leachTime_);

  // Initialize the list of all interfaces in the lattice

  cout << endl
       << "Controller::doCycle(...) Entering Lattice::findInterfaces()" << endl;

  lattice_->findInterfaces();

  // lattice_->checkSite(8);
  // cout << endl << " exit controller" << endl;// exit(0);

  cout << endl << "Controller::doCycle(...) Entering Main time loop" << endl;

  static double timestep = 0.0;
  bool capwater = true; // True if some capillary water is available
  time_index = 0;

  ///
  /// Output a file that directly links the microstructure ids to their
  /// rgb color.  This is only for easier image processing after the simulation
  /// is finished so we don't have to read the json file
  ///

  lattice_->writeMicroColors();

  ///
  /// Write the initial microstructure image and its png image
  ///

  lattice_->writeLattice(0.0);
  lattice_->writeLatticePNG(0.0);
  if (xyz_)
    lattice_->appendXYZ(0.0);
  int timesGEMFailed_loc = 0;

  // init to 0 all DC moles corresponding to the kinetic controlled microphases
  //      these DCmoles will be updated by
  //      KineticController::calculateKineticStep and passedd to GEM together
  //      the other DC moles in the stystem (ChemicalSystem::calculateState)
  int numMicPh = chemSys_->getNumMicroPhases();
  // cout << "numMicPh : " << numMicPh << endl;

  int DCId;
  for (int i = FIRST_SOLID; i < numMicPh; i++) {
    if (chemSys_->isKinetic(i)) {
      DCId = chemSys_->getMicroPhaseDCMembers(i, 0);
      // chemSys_->setDCMoles(DCId,0.0); //coment if DCLowerLimit in
      // kineticControllerStep/GEM_from_MT
      chemSys_->setIsDCKinetic(DCId, true);
    }
  }
  cout << endl
       << "***   numGEMPhases_  = " << chemSys_->getNumGEMPhases() << endl;
  cout << "***   numDCs_        = " << chemSys_->getNumDCs() << endl;
  cout << "***   numICs_        = " << chemSys_->getNumICs() << endl;

  // cout << "Starting with a pore solution without dissolved DCs  => all
  // microPhaseSI_ = 0" << endl; init to 0 all microPhaseSI_
  // chemSys_->setZeroMicroPhaseSI();

  bool writeICsDCs = true;
  if (writeICsDCs)
    writeTxtOutputFiles_onlyICsDCs(0); // to check the total ICs

  // variables used in DCLowerLimit computation
  double volMolDiff, molarMassDiff, vfracDiff, massDissolved;
  double microPhaseMassDiff, scaledMassDiff, numMolesDiff;
  int numDCs = chemSys_->getNumDCs();
  int timesGEMFailed_recall;

  cout << endl << endl << "     ===== START SIMULATION =====" << endl;

  int cyc = 0;
  int timeSize = time_.size();

  double minTime, rNum;
  double timeTemp;

  // JWB 2024-03-18: Florin had this set to 1.e-7
  // JWB 2024-03-18: I increased it due to going from days to hours for time
  // units
  double deltaTime = elemTimeInterval; // 1.e-6; // 1.e-5;
  double delta2Time_0 = 2.0 * deltaTime;
  double delta2Time;
  int numIntervals = 0, numMaxIntervals = 1;
  int numGenMax = 3000;
  int numGen = 0;

  int numTotGen = 0;
  double nextTimeStep;
  int fracNum = 10;
  double fracNextTimeStep;
  double timeZero;

  int lastGoodI = 0;
  double lastGoodTime = 0.0;

  // Main computation cycle
  for (i = 0; (i < timeSize) && (capwater); ++i) {

    ///
    /// Do not advance the time step if GEM_run failed the last time
    ///

    bool isFirst = (i == 0) ? true : false;

    cyc++;

    // new
    if (timesGEMFailed_loc > 0) {
      i--;
      time_[i] += (0.1 * (time_[i + 1] - time_[i]));
      timestep = time_[i] - lastGoodTime;

      cout << endl
           << endl
           << endl
           << "##### Controller::doCycle  GEMFailed => add next timestep "
              "before to START NEW CYCLE   "
              "i/cyc/time_[i]/lastGoodI/lastGoodTime/timestep: "
           << i << " / " << cyc << " / " << time_[i] << " / " << lastGoodI
           << " / " << lastGoodTime << " / " << timestep << " #####" << endl;

    } else {

      timestep = (i > 0) ? (time_[i] - time_[i - 1]) : time_[i];
      if (i == 0) {
        lastGoodTime = 0;
        lastGoodI = 0;
        cout << endl
             << endl
             << endl
             << "##### Controller::doCycle  START NEW CYCLE   "
                "i/cyc/time_[i]/timestep: "
             << i << " / " << cyc << " / " << time_[i] << " / " << timestep
             << " #####" << endl;
      } else {
        lastGoodTime = time_[i - 1];
        lastGoodI = i - 1;
        cout << endl
             << endl
             << endl
             << "##### Controller::doCycle  START NEW CYCLE   "
                "i/cyc/time_[i]/time_[i-1]/timestep: "
             << i << " / " << cyc << " / " << time_[i] << " / " << time_[i - 1]
             << " / " << timestep << " #####" << endl;
      }
    }

    /// Assume that only capillary pore water is chemically reactive,
    /// while water in nanopores is chemically inert.
    ///
    ///
    /// This is the main step of the cycle; the calculateState method
    /// runs all the major steps of a computational cycle

    try {

      chemSys_->initDCLowerLimit(0.0);
      timesGEMFailed_loc = calculateState(time_[i], timestep, isFirst, cyc);

    } catch (GEMException gex) {
      lattice_->writeLattice(time_[i]);
      lattice_->writeLatticePNG(time_[i]);
      if (xyz_)
        lattice_->appendXYZ(time_[i]);
      throw gex;
    }

    ///
    /// Once the change in state is determined, propagate the consequences
    /// to the 3D microstructure only if the GEM_run calculation succeeded.
    /// Otherwise we adjust the time step and try again
    ///

    if (timesGEMFailed_loc > 0) {
      cout << endl
           << "  Controller::doCycle first GEM_run failed "
              "i/cyc/time[i]/getTimesGEMFailed_loc : "
           << i << " / " << cyc << " / " << time_[i] << " / "
           << timesGEMFailed_loc << endl;

      //**********************

      cout << endl
           << "  Controller::doCycle - PRBL_0      i/cyc/time_[i]/timestep     "
              "   : "
           << i << " / " << cyc << " / " << time_[i] << " / " << timestep
           << "   =>   WAIT..." << endl;
      cout.flush();

      numTotGen = 0;
      nextTimeStep = time_[i + 1] - time_[i];
      fracNextTimeStep = nextTimeStep / fracNum;

      for (int indFracNum = 0; indFracNum < fracNum; indFracNum++) {
        timeZero = time_[i] + (((double)indFracNum) * fracNextTimeStep);
        minTime = timeZero - deltaTime;
        numGen = 0;
        numIntervals = 0;
        delta2Time = delta2Time_0;
        while (timesGEMFailed_loc > 0) {
          if (numGen % numGenMax == 0) {
            if (numGen > 0) {
              numIntervals++;
              delta2Time = delta2Time * 10.0;
              minTime = timeZero - (0.5 * delta2Time);
            }
            if (numIntervals == numMaxIntervals) {
              // cout << "      for
              // cyc/indFracNum/delta2Time/numGen/timeZero/minTime : " << cyc
              //      << " / " << indFracNum << " / " << delta2Time << " / " <<
              //      numGen << " / "
              //      << timeZero << " / " << minTime << endl;
              // cout << "         =>   numIncreaseInterval = " <<
              // numMaxIntervals
              //      << " (max val) => change indFracNum (next timeZero)!!!" <<
              //      endl;
              // cout.flush();
              break;
            }
          }

          rNum = lattice_->callRNG();
          numGen++;
          numTotGen++;
          timeTemp = minTime + (rNum * delta2Time);
          timestep = timeTemp - lastGoodTime;
          chemSys_->initDCLowerLimit(0.0);
          timesGEMFailed_loc = calculateState(timeTemp, timestep, isFirst, cyc);
          if (timesGEMFailed_loc == 0) {
            time_[i] = timeTemp;
            cout << "  Controller::doCycle - PRBL_0 solved for "
                    "i/cyc/time_[i]/timestep/numTotGen : "
                 << i << " / " << cyc << " / " << time_[i] << " / " << timestep
                 << " / " << numTotGen << endl;
            cout.flush();
            break;
          }
        }
        if (timesGEMFailed_loc == 0) {
          delta2Time = delta2Time_0;
          numIntervals = 0;
          break;
        }
      }

      if (timesGEMFailed_loc > 0) {
        cout << "  Controller::doCycle - PRBL_0 not solved for "
                "i/cyc/time_[i]/lastGoodI/lastGoodTime/timestep: "
             << i << " / " << cyc << " / " << time_[i] << " / " << lastGoodI
             << " / " << lastGoodTime << " / " << timestep << endl;
        continue;
      }
      //**********************

    } else {

      cout << endl
           << "  Controller::doCycle first GEM_run OK "
              "i/cyc/time[i]/getTimesGEMFailed_loc: "
           << i << " / " << cyc << " / " << time_[i] << " / "
           << timesGEMFailed_loc << endl;
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
      // set iniLattice i.e. a copy of initial lattice/system configuration
      // (including RNG state and all DCs values)
      // from ChemicalSystem:
      iniLattice.DCMoles = kineticController_->getDCMoles();

      // from Lattice:
      iniLattice.count = lattice_->getCount();
      iniLattice.growthInterfaceSize = lattice_->getGrowthInterfaceSize();
      iniLattice.dissolutionInterfaceSize =
          lattice_->getDissolutionInterfaceSize();
      iniLattice.site.clear();
      RestoreSite site_l; // only one declaration
      for (int i = 0; i < numSites_; i++) {
        site_l.microPhaseId = (lattice_->getSite(i))->getMicroPhaseId();
        site_l.growth = (lattice_->getSite(i))->getGrowthPhases();
        site_l.wmc = (lattice_->getSite(i))->getWmc();
        site_l.wmc0 = (lattice_->getSite(i))->getWmc0();
        site_l.visit = 0;
        site_l.inGrowInterfacePos =
            (lattice_->getSite(i))->getInGrowInterfacePosVector();
        site_l.inDissInterfacePos =
            (lattice_->getSite(i))->getInDissInterfacePos();
        iniLattice.site.push_back(site_l);
      }
      iniLattice.interface.clear();
      RestoreInterface interface_l; // only one declaration
      int dimLatticeInterface =
          lattice_->getInterfaceSize(); // only one declaration
      for (int i = 0; i < dimLatticeInterface; i++) {
        interface_l.microPhaseId = lattice_->getInterface(i).getMicroPhaseId();
        interface_l.growthSites = lattice_->getInterface(i).getGrowthSites();
        interface_l.dissolutionSites =
            lattice_->getInterface(i).getDissolutionSites();
        iniLattice.interface.push_back(interface_l);
      }
      iniLattice.numRNGcall_0 = lattice_->getNumRNGcall_0();
      iniLattice.numRNGcallLONGMAX = lattice_->getNumRNGcallLONGMAX();
      iniLattice.lastRNG = lattice_->getLastRNG();

      int phId = 0;
      int changeLattice = -100;
      int whileCount = 0;

      vector<int> numSitesNotAvailable;
      numSitesNotAvailable.clear();
      vector<int> vectPhIdDiff;
      vectPhIdDiff.clear();
      vector<string> vectPhNameDiff;
      vectPhNameDiff.clear();

      changeLattice = lattice_->changeMicrostructure(
          time_[i], simType_, capwater, numSitesNotAvailable, vectPhIdDiff,
          vectPhNameDiff, whileCount, cyc);

      // if not all the voxels requested by KM/GEM for a certain microphase
      //  phDiff (DCId)can be dissolved because of the system configuration
      //  (lattice):
      //   - comeback to the initial system configuration : iniLattice
      //   - re-run GEM with restrictions impossed by the system configuration
      //       (DC distribution on the lattice sites) i.e. the primal solution
      //        must contain a number of moles corresponding to
      //        numSitesNotAvailable ("numDiff" in
      //        lattice_->changeMicrostructure" - sites that cannot be
      //        dissolved) lattice sites for the microphase phDiff

      if (changeLattice == 0) {
        timesGEMFailed_recall = -1;
        bool testDiff = true;
        while (changeLattice == 0) { // - for many phases!
          whileCount++;

          cout << endl
               << "  Controller::doCycle - cyc = " << cyc
               << " :  changeLattice = " << changeLattice
               << "  =>  whileCount = " << whileCount << endl;

          while (timesGEMFailed_recall != 0) {

            // reset for ChemicalSystem:
            for (int i = 0; i < numDCs; i++) {
              chemSys_->setDCMoles(i, iniLattice.DCMoles[i]);
            }

            // reset for Lattice:
            lattice_->setCount(iniLattice.count);
            for (int i = 0; i < numSites_; i++) {
              (lattice_->getSite(i))
                  ->setMicroPhaseId(iniLattice.site[i].microPhaseId);
              (lattice_->getSite(i))
                  ->setGrowthPhases(iniLattice.site[i].growth);
              (lattice_->getSite(i))->setWmc(iniLattice.site[i].wmc);
              (lattice_->getSite(i))->setWmc0(iniLattice.site[i].wmc0);
              (lattice_->getSite(i))
                  ->setVisit(iniLattice.site[i].visit); // or 0!
              (lattice_->getSite(i))
                  ->setInGrowInterfacePosVector(
                      iniLattice.site[i].inGrowInterfacePos);
              (lattice_->getSite(i))
                  ->setInDissInterfacePos(
                      iniLattice.site[i].inDissInterfacePos);
            }
            for (int i = 0; i < dimLatticeInterface; i++) {
              lattice_->setInterfaceMicroPhaseId(
                  i, iniLattice.interface[i].microPhaseId); // same as before!
              lattice_->setGrowthSites(i, iniLattice.interface[i].growthSites);
              lattice_->setDissolutionSites(
                  i, iniLattice.interface[i].dissolutionSites);
            }
            lattice_->setGrowthInterfaceSize(iniLattice.growthInterfaceSize);
            lattice_->setDissolutionInterfaceSize(
                iniLattice.dissolutionInterfaceSize);
            lattice_->resetRNG(iniLattice.numRNGcall_0,
                               iniLattice.numRNGcallLONGMAX, iniLattice.lastRNG,
                               cyc, whileCount);

            cout << "  Controller::doCycle - cyc = " << cyc
                 << " :  reset system OK & GEM_run recall for "
                    "i/whileCount/numSitesNotAvailable.size() = "
                 << i << " / " << whileCount << " / "
                 << numSitesNotAvailable.size() << endl;
            cout << "  Controller::doCycle - cyc = " << cyc
                 << " :  reset DCLowerLimits :" << endl;

            for (int i = 0; i < numSitesNotAvailable.size(); i++) {

              phId = vectPhIdDiff[i];

              DCId = chemSys_->getMicroPhaseDCMembers(phId, 0);

              volMolDiff = chemSys_->getDCMolarVolume(DCId);  // m3/mol
              molarMassDiff = chemSys_->getDCMolarMass(DCId); // g/mol

              vfracDiff =
                  ((double)numSitesNotAvailable[i]) / ((double)numSites_);

              microPhaseMassDiff =
                  vfracDiff * molarMassDiff / volMolDiff / 1.0e6; // g/cm3

              scaledMassDiff =
                  microPhaseMassDiff * 100.0 / lattice_->getInitSolidMass();

              if (chemSys_->isKinetic(phId)) {
                massDissolved = kineticController_->updateKineticStep(
                    cyc, phId, scaledMassDiff);
              } else {

                cout << "    Controller::doCycle - not a KM phase - for cyc = "
                     << cyc << " & phaseId = " << phId << " ["
                     << chemSys_->getMicroPhaseName(phId) << " / DCId:" << DCId
                     << "]" << endl;

                numMolesDiff = scaledMassDiff / molarMassDiff;

                cout << "      DCMoles_/keepNumDCMoles : "
                     << chemSys_->getDCMoles(DCId) << " / " << numMolesDiff
                     << endl;

                chemSys_->setDCLowerLimit(DCId, numMolesDiff);
              }
            }

            cout << endl
                 << "  Controller::doCycle - cyc = " << cyc
                 << " :  #  i#/ "
                    "phName/phId/count/dissInterfaceSize/numSitesNotAvailable"
                    "/DCId/DCMoles/DCLowerLimit :"
                 << endl;
            for (int i = 0; i < numSitesNotAvailable.size(); i++) {
              phId = vectPhIdDiff[i];
              DCId = chemSys_->getMicroPhaseDCMembers(phId, 0);
              cout << "                        cyc = " << cyc << " :  #"
                   << setw(3) << right << i << "#/ " << setw(15) << left
                   << vectPhNameDiff[i] << "   " << setw(5) << right << phId
                   << "   " << setw(9) << right << lattice_->getCount(phId)
                   << "   " << setw(9) << right
                   << lattice_->getDissolutionInterfaceSize(phId) << "   "
                   << setw(9) << right << numSitesNotAvailable[i] << "   "
                   << setw(5) << right << DCId << "   "
                   << chemSys_->getDCMoles(DCId) << "   "
                   << chemSys_->getDCLowerLimit(DCId) << endl;
            }
            cout << endl;

            timesGEMFailed_recall =
                chemSys_->calculateState(time_[i], isFirst, cyc);

            cout << endl
                 << "  Controller::doCycle - cyc = " << cyc
                 << " :  i/time[i]/getTimesGEMFailed_recall = " << i << " / "
                 << time_[i] << " / " << timesGEMFailed_recall << endl;

            if (timesGEMFailed_recall > 0) {
              cout << "  Controller::doCycle - GEM_run failed for whileCount = "
                   << whileCount << endl;
              cout.flush();
              timesGEMFailed_loc = timesGEMFailed_recall;
              testDiff = false;
              for (int iii = 0; iii < numSitesNotAvailable.size(); iii++) {
                phId = vectPhIdDiff[iii];
                if (lattice_->getCount(phId) > numSitesNotAvailable[iii]) {
                  numSitesNotAvailable[iii]++;
                  cout << "  Controller::doCycle - for i/cyc/phId/iii = " << i
                       << " / " << cyc << " / " << phId << " / " << iii
                       << "   =>   numSitesNotAvailable[iii] = "
                       << numSitesNotAvailable[iii] << endl;
                  cout.flush();
                  testDiff = true;
                  break;
                }
              }
              if (!testDiff) {
                cout
                    << "Controller::doCycle - do not update the microstructure "
                    << endl;
                cout.flush();
                break;
              }

            } else {
              cout << "  Controller::doCycle - cyc = " << cyc
                   << " :  GEM_run OK for whileCount = " << whileCount << endl;
              testDiff = true;
            }
          }

          if (testDiff) {

            numSitesNotAvailable.clear();
            vectPhIdDiff.clear();
            vectPhNameDiff.clear();
            changeLattice = lattice_->changeMicrostructure(
                time_[i], simType_, capwater, numSitesNotAvailable,
                vectPhIdDiff, vectPhNameDiff, whileCount, cyc);
            cout << endl
                 << "  Controller::doCycle - cyc = " << cyc
                 << "  &  whileCount = " << whileCount
                 << "  :  timesGEMFailed_recall = " << timesGEMFailed_recall
                 << "  &  changeLattice = " << changeLattice << endl;
            if (timesGEMFailed_recall == 0 && changeLattice == 0) {
              cout << "     => redo adjustment for cyc/whileCount = " << cyc
                   << " / " << whileCount << endl;
              timesGEMFailed_recall = -1;
              testDiff = true;
            }
          } else {
            break;
          }
        } // while (changeLattice == 0) { // - for many phases!

        if (timesGEMFailed_loc > 0) {
          continue;
        } else {
          kineticController_->setHydTimeIni(time_[i]);
          cout << endl
               << "Controller::doCycle => normal end after reset system for "
                  "cyc = "
               << cyc << " (i = " << i << ")" << endl;
        }
      } else {
        kineticController_->setHydTimeIni(time_[i]);

        cout << endl
             << "Controller::doCycle => normal end - cyc = " << cyc
             << " (i = " << i << ")" << endl;
      }

    } catch (DataException dex) {
      dex.printException();
      lattice_->writeLattice(time_[i]);
      lattice_->writeLatticePNG(time_[i]);
      if (xyz_)
        lattice_->appendXYZ(time_[i]);
      throw dex;
    } catch (EOBException ex) {
      ex.printException();
      lattice_->writeLattice(time_[i]);
      lattice_->writeLatticePNG(time_[i]);
      if (xyz_)
        lattice_->appendXYZ(time_[i]);
      throw ex;
    } catch (MicrostructureException mex) {
      mex.printException();
      lattice_->writeLattice(time_[i]);
      lattice_->writeLatticePNG(time_[i]);
      if (xyz_)
        lattice_->appendXYZ(time_[i]);
      throw mex;
    }

    // write output .txt files
    writeTxtOutputFiles(time_[i]);
    if (writeICsDCs)
      writeTxtOutputFiles_onlyICsDCs(time_[i]);

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
      cout << "Controller::doCycle Returned from Lattice::changeMicrostructure"
           << endl;
      cout.flush();
    }

    if ((time_[i] >= outputTime_[time_index]) &&
        (time_index < outputTime_.size())) {
      if (verbose_) {
        cout << "Controller::doCycle Writing lattice at time_[" << i
             << "] = " << time_[i] << ", outputTime_[" << time_index
             << "] = " << outputTime_[time_index] << endl;
      }

      lattice_->writeLattice(time_[i]);
      lattice_->writeLatticePNG(time_[i]);

      if (xyz_)
        lattice_->appendXYZ(time_[i]);

      lattice_->writePoreSizeDistribution(time_[i]);

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

    if (time_[i] >= sulfateAttackTime_) {

      cout << endl
           << " Controller::doCycle - for sulfate attack, check conditions for "
              "addDissolutionSites & coordination sphere "
           << endl;
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

      // expansion.clear();

      if (expansion.size() > 1) {

        damageCount_ = 0;
        double poreintroduce = 0.5;

        if (verbose_) {
          cout << "Controller::doCycle Sulfate attack module writing " << endl;
          cout << "Controller::doCycle lattice at time_[" << i
               << "] = " << time_[i] << ", " << endl;
          cout << "controller::doCycle outputTime_[" << time_index
               << "] = " << outputTime_[time_index] << endl;
          cout.flush();
        }

        lattice_->writeLattice(time_[i]);
        lattice_->writeLatticePNG(time_[i]);
        if (xyz_)
          lattice_->appendXYZ(time_[i]);
        string ofileName(jobRoot_);
        ostringstream ostr1, ostr2;
        ostr1 << static_cast<int>(time_[i] * 10.0); // tenths of an hour
        ostr2 << setprecision(3) << temperature_;
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
          lattice_->dWaterChange(poreintroduce);
          */

          double dwmcval = poreintroduce;
          lattice_->dWmc(expindex, dwmcval);
          for (int j = 0; j < NN_NNN; j++) { // ste->nbSize(2)
            Site *stenb = ste->nb(j);
            stenb->dWmc(dwmcval);
          }
        }

        ///
        /// Calculate the stress-free strain (thermal strain) state in the
        /// microstructure, and then write the displacement field
        ///

        thermalstr_->Calc(time_[i], ofileName, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

        // thermalstr_ -> writeStress(jobRoot_,time_[i],0); //write strxx
        // thermalstr_ -> writeStrainEngy(jobRoot_,time_[i]);
        thermalstr_->writeDisp(jobRoot_, time_[i]);

        ///
        /// Get the true volume of each voxel after FEM calculation
        ///

        double otruevolume = 0.0;
        double truevolume = 0.0;
        for (int ii = 0; ii < numSites_; ii++) {
          Site *ste;
          ste = lattice_->getSite(i);
          otruevolume = ste->getTrueVolume();
          for (int j = 0; j < 3; j++) {
            truevolume += thermalstr_->getEleStrain(ii, j);
          }
          truevolume = otruevolume * (1 + truevolume);
          ste->setTrueVolume(truevolume);
        }

        for (int index = 0; index < numSites_; index++) {
          Site *ste;
          ste = lattice_->getSite(index);
          int pid = ste->getMicroPhaseId();

          if (ste->IsDamage()) {
            damageCount_++;
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
              lattice_->dWaterChange(poreincrease);
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
                damageCount_++;
                // lattice_->dWaterChange(poreintroduce);

                double dwmcval = poreintroduce;
                lattice_->dWmc(index, dwmcval);
                for (int j = 0; j < NUM_NEAREST_NEIGHBORS; j++) { // NN_NNN?
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
                lattice_->dWaterChange(poreindamage);
                ste->setVolume(VOIDID,poreindamage);
                */
              }
            }
          }

        } // End of loop over all voxels

        if (verbose_) {
          cout << "Controller::doCycle sulfate attack module Time = "
               << time_[i] << " damageCount_ is: " << damageCount_ << endl;
          cout.flush();
        }
        ofstream outdamage("damage.dat");
        outdamage << damageCount_;
        outdamage.close();

        lattice_->writeDamageLattice(time_[i]);
        lattice_->writeDamageLatticePNG(time_[i]);
        // to see whether new damage is generated
      }
    }
  }

  ///
  /// Write the final lattice state to an ASCII file and to a PNG file for
  /// visualization
  ///

  lattice_->writeLattice(time_[i - 1]);
  lattice_->writeLatticePNG(time_[i - 1]);
  if (xyz_)
    lattice_->appendXYZ(time_[i - 1]);

  return;
}

int Controller::calculateState(double time, double dt, bool isFirst, int cyc) {

  int timesGEMFailed = 0;

  try {

    ///
    /// We must pass some vectors to the `calculateKineticStep` method that
    /// will hold the amounts of impurity elements released from the clinker
    /// phases.  These values do not go in the `ChemicalSystem` object, but will
    /// still need to be processed afterward.
    ///

    // vector<double> impurityrelease;
    // impurityrelease.clear();
    // impurityrelease.resize(chemSys_->getNumMicroImpurities(), 0.0);

    ///
    /// Get the number of moles of each IC dissolved from kinetically controlled
    /// phases
    ///

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

    chemSys_->setMicroPhaseSI();

    kineticController_->calculateKineticStep(dt, cyc);

    // if (time >= sulfateAttackTime_) {// for sulfate attack iterations}

    ///
    /// Now that the method is done determining the change in moles of each IC,
    /// launch a thermodynamic calculation to determine new equilibrium state
    ///
    /// The `ChemicalSystem` object provides an interface for these calculations
    ///

    try {
      timesGEMFailed = chemSys_->calculateState(time, isFirst, cyc);
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

    ///
    /// The thermodynamic calculation returns the saturation index of phases,
    /// which is needed for calculations of driving force for dissolution
    /// or growth.  Assign this to the lattice in case crystallization pressures
    /// should be calculated.
    ///

    if (timesGEMFailed > 0)
      return timesGEMFailed;

  } catch (FileException fex) {
    fex.printException();
    exit(1);
  } catch (FloatException flex) {
    flex.printException();
    exit(1);
  }
  return timesGEMFailed;
}

void Controller::writeTxtOutputFiles(double time) {
  ///
  /// Set the kinetic DC moles.  This adds the clinker components to the DC
  /// moles.
  ///
  /// @todo Find out what this is and why it needs to be done
  ///
  ///
  // int numICs = chemSys_->getNumICs();
  // int numDCs = chemSys_->getNumDCs();
  int i, j;

  // Output to files the solution composition data, phase data, DC data,
  // microstructure data, pH, and C-S-H composition and Ca/Si ratio

  string outfilename = jobRoot_ + "_Solution.csv";
  ofstream outfs(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "calculateState", outfilename,
                        "Could not append");
  }

  outfs << setprecision(5) << time;
  char cc;
  for (i = 0; i < numDCs_; i++) {
    cc = chemSys_->getDCClassCode(i);
    if (cc == 'S' || cc == 'T' || cc == 'W') {
      outfs << "," << (chemSys_->getNode())->Get_cDC((long int)i); // molality
    }
  }
  outfs << endl;
  outfs.close();

  outfilename = jobRoot_ + "_DCVolumes.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "calculateState", outfilename,
                        "Could not append");
  }

  outfs << setprecision(5) << time;
  for (i = 0; i < numDCs_; i++) {
    if (chemSys_->getDCMolarMass(i) > 0.0) {
      cc = chemSys_->getDCClassCode(i);
      if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M' || cc == 'W') {
        string dcname = chemSys_->getDCName(i);
        double V0 =
            chemSys_->getDCMoles(dcname) * chemSys_->getDCMolarVolume(dcname);
        outfs << "," << V0;
      }
    } else {
      string msg = "Divide by zero error for DC " + chemSys_->getDCName(i);
      outfs.close();
      throw FloatException("Controller", "calculateState", msg);
    }
  }
  outfs << endl;
  outfs.close();

  outfilename = jobRoot_ + "_SurfaceAreas.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "Controller", outfilename,
                        "Could not append");
  }
  outfilename = jobRoot_ + "_SI.csv";
  ofstream outfs01;
  outfs01.open(outfilename.c_str(), ios::app);
  if (!outfs01) {
    throw FileException("Controller", "Controller", outfilename,
                        "Could not append");
  }

  outfs << setprecision(5) << time;
  outfs01 << setprecision(5) << time;
  // JWB: Start with micro phase id 1 to avoid the Void phase
  for (int i = 1; i < chemSys_->getNumMicroPhases(); i++) {
    int dcid = chemSys_->getMicroPhaseDCMembers(i, 0);
    cc = chemSys_->getDCClassCode(dcid);
    if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M') {
      outfs << "," << lattice_->getSurfaceArea(i);
      ;
      outfs01 << "," << chemSys_->getMicroPhaseSI(i);
    }
  }
  outfs << endl;
  outfs01 << endl;
  outfs.close();
  outfs01.close();

  outfilename = jobRoot_ + "_Microstructure.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "calculateState", outfilename,
                        "Could not append");
  }

  outfs << setprecision(5) << time;
  for (i = 0; i < numMicroPhases_; i++) {
    outfs << "," << (lattice_->getVolumeFraction(i));
  }
  double micvol = lattice_->getMicrostructureVolume();
  double initmicvol = lattice_->getInitialMicrostructureVolume();
  outfs << "," << micvol << "," << (initmicvol - micvol);
  outfs << endl;
  outfs.close();

  outfilename = jobRoot_ + "_pH.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "calculateState", outfilename,
                        "Could not append");
  }
  outfs << setprecision(5) << time;
  outfs << "," << (chemSys_->getPH()) << endl;
  outfs.close();

  chemSys_->setGEMPhaseStoich();

  double *CSHcomp;
  try {
    CSHcomp = chemSys_->getPGEMPhaseStoich(chemSys_->getGEMPhaseId(CSHGEMName));
  } catch (EOBException eex) {
    eex.printException();
    exit(1);
  }
  if (verbose_) {
    cout << "Done!" << endl;
    cout.flush();
  }

  outfilename = jobRoot_ + "_CSH.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "calculateState", outfilename,
                        "Could not append");
  }
  outfs << setprecision(5) << time;
  for (i = 0; i < numICs_; i++) {
    outfs << "," << CSHcomp[i];
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
  outfs << "," << CaSiRatio << endl;
  outfs.close();

  double *phaseRecord;
  CaMoles = SiMoles = 0.0;
  outfilename = jobRoot_ + "_CSratio_solid.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "calculateState", outfilename,
                        "Could not append");
  }
  outfs << setprecision(5) << time;
  for (i = 0; i < numGEMPhases_; i++) {
    cc = chemSys_->getGEMPhaseClassCode(i);
    if (cc == 's') {
      phaseRecord = chemSys_->getPGEMPhaseStoich(i);
      // ICIndex = chemSys_->getICId("Ca");
      // CaMoles += phaseRecord[ICIndex];
      // ICIndex = chemSys_->getICId("Si");
      // SiMoles += phaseRecord[ICIndex];
      CaMoles += phaseRecord[id_Ca];
      SiMoles += phaseRecord[id_Si];
    }
  }
  if (SiMoles != 0) {
    CaSiRatio = CaMoles / SiMoles;
    outfs << "," << CaSiRatio << endl;
  } else {
    outfs << ",Si_moles is ZERO" << endl;
  }
  outfs.close();

  outfilename = jobRoot_ + "_Enthalpy.csv";
  outfs.open(outfilename.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "calculateState", outfilename,
                        "Could not append");
  }

  double enth = 0.0;
  for (i = 0; i < numDCs_; i++) {
    enth += (chemSys_->getDCEnthalpy(i));
  }

  outfs << setprecision(5) << time;
  outfs << "," << enth << endl;
  outfs.close();
}

void Controller::writeTxtOutputFiles_onlyICsDCs(double time) {

  int i, j;
  vector<double> ICMoles;
  ICMoles.resize(numICs_, 0.0);
  vector<double> DCMoles;
  DCMoles.resize(numDCs_, 0.0);
  for (i = 0; i < numDCs_; i++) {
    DCMoles[i] = chemSys_->getDCMoles(i);
  }

  string outfilenameIC = jobRoot_ + "_icmoles.csv";
  string outfilenameDC = jobRoot_ + "_dcmoles.csv";
  ofstream outfs;
  if (time < 1.e-10) {
    outfs.open(outfilenameIC.c_str());
    if (!outfs) {
      throw FileException("Controller", "writeTxtOutputFiles_onlyICsDCs",
                          outfilenameIC, "Could not append");
    }
    outfs << "Time(h)";
    for (i = 0; i < numICs_; i++) {
      outfs << "," << chemSys_->getICName(i);
    }
    outfs << endl;
    outfs.close();

    outfs.open(outfilenameDC.c_str());
    if (!outfs) {
      throw FileException("Controller", "writeTxtOutputFiles_onlyICsDCs",
                          outfilenameDC, "Could not append");
    }
    outfs << "Time(h)";
    for (i = 0; i < numDCs_; i++) {
      outfs << "," << chemSys_->getDCName(i);
    }
    outfs << endl;
    outfs.close();
  }

  vector<int> impurityDCID;
  impurityDCID.clear();
  impurityDCID.push_back(chemSys_->getDCId("K2O"));
  impurityDCID.push_back(chemSys_->getDCId("Na2O"));
  impurityDCID.push_back(chemSys_->getDCId("Per"));
  impurityDCID.push_back(chemSys_->getDCId("SO3"));

  double scMass;
  int mPhId;
  double massImpurity, totMassImpurity;

  // cout << endl << "getIsDCKinetic: " << endl;
  for (j = 0; j < numDCs_; j++) {
    if (chemSys_->getIsDCKinetic(j)) {
      // molMass = chemSys_->getDCMolarMass(j);
      mPhId = chemSys_->getDC_to_MPhID(j);
      scMass = chemSys_->getMicroPhaseMass(mPhId);

      totMassImpurity = 0;

      massImpurity = scMass * chemSys_->getK2o(mPhId);
      totMassImpurity += massImpurity;
      DCMoles[impurityDCID[0]] +=
          massImpurity / chemSys_->getDCMolarMass("K2O");

      massImpurity = scMass * chemSys_->getNa2o(mPhId);
      totMassImpurity += massImpurity;
      DCMoles[impurityDCID[1]] +=
          massImpurity / chemSys_->getDCMolarMass("Na2O");

      massImpurity = scMass * chemSys_->getMgo(mPhId);
      totMassImpurity += massImpurity;
      DCMoles[impurityDCID[2]] +=
          massImpurity / chemSys_->getDCMolarMass("Per"); // MgO

      massImpurity = scMass * chemSys_->getSo3(mPhId);
      totMassImpurity += massImpurity;
      DCMoles[impurityDCID[3]] +=
          massImpurity / chemSys_->getDCMolarMass("SO3");

      DCMoles[j] = (scMass - totMassImpurity) / chemSys_->getDCMolarMass(j);
    }
  }

  for (j = 0; j < numDCs_; j++) {
    for (i = 0; i < numICs_; i++) {
      ICMoles[i] += DCMoles[j] * chemSys_->getDCStoich(j, i);
    }
  }

  outfs.open(outfilenameIC.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "writeTxtOutputFiles_onlyICsDCs",
                        outfilenameIC, "Could not append");
  }

  outfs << setprecision(5) << time;
  for (i = 0; i < numICs_; i++) {
    outfs << "," << ICMoles[i];
  }
  outfs << endl;
  outfs.close();

  outfs.open(outfilenameDC.c_str(), ios::app);
  if (!outfs) {
    throw FileException("Controller", "writeTxtOutputFiles_onlyICsDCs",
                        outfilenameDC, "Could not append");
  }

  outfs << setprecision(5) << time;
  for (i = 0; i < numDCs_; i++) {
    outfs << "," << DCMoles[i];
  }
  outfs << endl;
  outfs.close();
}

void Controller::parseDoc(const string &docName) {

  /// check if the JSON file exists

  ifstream f(docName.c_str());
  if (!f.is_open()) {
    cout << "JSON output times file not found" << endl;
    throw FileException("Controller", "parseDoc", docName, "File not found");
  }

  /// Parse the JSON file all at once

  json data = json::parse(f);
  f.close();

  try {

    /// Get an iterator to the root node of the JSON file
    /// @todo Add a better JSON validity check.

    json::iterator it = data.find("time_parameters");
    json::iterator cdi = it.value().find("finaltime");

    if (cdi == it.value().end() || it == data.end()) {
      throw FileException("Controller", "parseDoc", docName, "Empty JSON file");
    }

    // Input times are conventionally in days
    // Immediately convert to hours within model
    double finalTime = cdi.value();
    finalTime *= (H_PER_DAY);

    // Input times are conventionally in days
    // Immediately convert to hours within model
    cdi = it.value().find("outtimes");
    double testTime = 0.0;
    int outtimenum = cdi.value().size();
    for (int i = 0; i < outtimenum; ++i) {
      testTime = cdi.value()[i];
      testTime *= (H_PER_DAY);
      outputTime_.push_back(testTime);
    }

    // Now populate the calculation times on a natural log scale
    time_.clear();
    testTime = 0.0;
    int j = 0;
    bool done = false;
    while (testTime < finalTime) {
      testTime += (0.1 * (testTime + 0.024));
      done = false;
      if (j < outputTime_.size()) {
        while (j < outputTime_.size() && !done) {
          if (testTime >= outputTime_[j]) {
            time_.push_back(outputTime_[j]);
            j++;
          } else if (testTime < finalTime) {
            time_.push_back(testTime);
            done = true;
          }
        }
      } else if (testTime < finalTime) {
        time_.push_back(testTime);
      } else {
        time_.push_back(finalTime);
      }
    }

  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  return;
}
