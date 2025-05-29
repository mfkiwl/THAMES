/**
@file Controller.cc
@brief Definition of Controller class methods
*/

#include "Controller.h"

Controller::Controller(Lattice *msh, KineticController *kc, ChemicalSystem *cs,
                       ThermalStrain *thmstr, const int simtype,
                       const std::string &jsonFileName,
                       const std::string &jobname, const bool verbose,
                       const bool warning, const bool xyz) {

  xyz_ = xyz;
  simType_ = simtype;
  chemSys_ = cs;
  kineticController_ = kc;
  lattice_ = msh;
  thermalstr_ = thmstr;
  jobRoot_ = jobname;

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

  // imgFreq_ = 168.0; // hours = 7 days
  imgFreq_ = -1;

  leachTime_ = 1.0e10;
  sulfateAttackTime_ = 1.0e10;

  beginAttackTime_ = -1.0;
  endAttackTime_ = -1.0;
  attackTimeInterval_ = -1.0;

  oldDamageCount_ = 0;
  allDamageCount_ = 0;

  isParrotKilloh_ = chemSys_->getIsParrotKilloh();
  sizePK_ = isParrotKilloh_.size();

  numDCs_ = chemSys_->getNumDCs();
  numICs_ = chemSys_->getNumICs();
  numMicroPhases_ = chemSys_->getNumMicroPhases();
  numGEMPhases_ = chemSys_->getNumGEMPhases();

  ///
  /// Load up the pointers to the `ChemicalSystem` object and `Lattice` object
  ///

  // chemSys_ = lattice_->getChemSys();
  lattice_->setJobRoot(jobRoot_);

  ///
  /// Output the class codes for the solution and for DC components.
  /// Output the header for the microstructure phase stats file
  /// Output header for the file tracking pH
  /// Output header for the file tracking the C-S-H composition and Ca/Si ratios
  /// Output header for the file tracking the IC moles in the system
  ///

  ofstream outfs;
  std::string outfilename;

  try {
    outfilename = jobRoot_ + "_Solution.csv";
    outfs.open(outfilename.c_str());
    if (!outfs) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    char cc;
    outfs << "Time(h)";
    for (int i = 0; i < numDCs_; i++) {
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
    for (int i = 0; i < numDCs_; i++) {
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
    for (int i = 1; i < numMicroPhases_; i++) {
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
    for (int i = 0; i < chemSys_->getNumICs(); i++) {
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

  temperature_ = chemSys_->getTemperature();
  lattice_->setTemperature(temperature_);
  presure_ = chemSys_->getP();
  waterDCId_ = chemSys_->getDCId("H2O@");
  waterMolarMass_ = chemSys_->getDCMolarMass("H2O@");
  numSites_ = lattice_->getNumSites();
  initMicroVolume_ = chemSys_->getInitMicroVolume();

  cout << endl << "Controller::Controller(...) :" << endl;

  cout << endl
       << "   numGEMPhases_   = " << setw(3) << right << numGEMPhases_ << endl;
  cout << "   numMicroPhases_ = " << setw(3) << right << numMicroPhases_
       << endl;
  cout << "   numDCs_         = " << setw(3) << right << numDCs_ << endl;
  cout << "   numICs_         = " << setw(3) << right << numICs_ << endl;
  cout << "   waterDCId_      = " << setw(3) << right << waterDCId_
       << " (waterDCName = \"" << chemSys_->getDCName(waterDCId_) << "\")"
       << endl;

  ///
  /// Output a file that directly links the microstructure ids to their
  /// rgb color.  This is only for easier image processing after the simulation
  /// is finished so we don't have to read the xml file
  ///

  lattice_->writeMicroColors();

  ///
  /// Write the initial microstructure image and its png image
  ///

  cout << endl
       << "Controller::Controller(...) - write initial "
          "microstructure files (writeLattice(0.0), etc)"
       << endl;

  TimeStruct resolvedTime = getResolvedTime(0.0);
  lattice_->writeLattice(0.0, resolvedTime);
  lattice_->writeLatticePNG(0.0, resolvedTime);
  if (xyz_)
    lattice_->appendXYZ(0.0);

  ///
  /// Open and read the Controller parameter file
  ///

  std::string jsonext = ".json";
  size_t foundjson;

  try {
    foundjson = jsonFileName.find(jsonext);

    if (foundjson != std::string::npos) {
      parseDoc(jsonFileName);
    } else {
      cout << "Parameter file must be JSON" << endl;
      throw FileException("Controller", "Controller", jsonFileName,
                          "NOT JSON FORMAT");
    }
  } catch (FileException fex) {
    throw fex;
  }

  for (int i = 0; i < static_cast<int>(time_.size() - 1); i++) {
    if (abs(time_[i] - time_[i + 1]) <= 1.0e-6) {
      time_.erase(time_.begin() + i);
    }
  }

  int time_Size = time_.size();
  int outputTime_Size = outputTime_.size();

  if (time_Size == 0) {
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

  outfilename = jobRoot_ + "-times_used.json";
  outfs.open(outfilename.c_str());
  if (!outfs) {
    throw FileException("Controller", "Controller", outfilename,
                        "Could not append");
  }

  outfs << "{" << endl;
  outfs << "  \"time_parameters\": {" << endl;
  outfs << "    \"calctimes\": [" << endl;

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

  if (attack_) {
    outfs << "    ]," << endl;
    outfs << "    \"beginattacktime\": " << beginAttackTime_ << "," << endl;
    outfs << "    \"endattacktime\": " << endAttackTime_ << "," << endl;
    outfs << "    \"attacktimeinterval\": " << attackTimeInterval_ << endl;
    outfs << "  }" << endl;
    outfs << "}" << endl;
  } else {
    outfs << "    ]," << endl;
    outfs << "    \"beginattacktime\": -1," << endl;
    outfs << "    \"endattacktime\": -1," << endl;
    outfs << "    \"attacktimeinterval\": -1" << endl;
    outfs << "  }" << endl;
    outfs << "}" << endl;
  }

  outfs.close();

  imgFreq_ *= 24.0; // !!!!!

  cout << endl
       << "   => new time values (calctime & outtime) have been used and "
          "writen as :"
       << endl;
  cout << "         " << outfilename << endl;

  leachTime_ = beginAttackTime_; // not default leachTime_ = 1.0e10
  sulfateAttackTime_ =
      beginAttackTime_; // not default sulfateAttackTime_ = 1.0e10

  chemSys_->setIniAttackTime(sulfateAttackTime_);
  lattice_->setSulfateAttackTime(sulfateAttackTime_);
  lattice_->setLeachTime(leachTime_);

  kineticController_->setIniAttackTime(sulfateAttackTime_);

  if (simType_ == SULFATE_ATTACK) {
    cout << endl << "   => attack = " << attack_ << endl;
    cout << "   parameters in hours:" << endl;
    cout << "     -> beginattacktime = " << setw(7) << right
         << (int)beginAttackTime_ << endl;
    cout << "     -> endattacktime = " << setw(7) << right
         << (int)endAttackTime_ << endl;
    cout << "     -> attacktimeinterval = " << setw(7) << right
         << (int)attackTimeInterval_ << endl;

    // if (simType_ == SULFATE_ATTACK)
    lattice_->createGrowingVectSA();
  }
}

void Controller::doCycle(double elemTimeInterval) {
  int i;
  int time_index;
  RestoreSystem iniLattice;

  ///
  /// This block arbitrarily sets the leaching initiation time to 100 days if
  /// the leaching module is to be run, or sets the sulfate attack initiation
  /// time to 100 days if the sulfate attack module is to be run
  ///
  /// @todo Think about generalizing this more, or allowing combinations of more
  /// than one (LEACHING/SULFATE_ATTACK/...)
  ///

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

  int timesGEMFailed_loc = 0;

  // init to 0 all DC moles corresponding to the kinetic controlled microphases
  //      these DCmoles will be updated by
  //      KineticController::calculateKineticStep and passedd to GEM together
  //      the other DC moles in the stystem (ChemicalSystem::calculateState)
  // int numMicPh = chemSys_->getNumMicroPhases();
  // cout << "numMicPh : " << numMicPh << endl;

  int DCId;
  for (int i = FIRST_SOLID; i < numMicroPhases_; i++) {
    if (chemSys_->isKinetic(i)) {
      DCId = chemSys_->getMicroPhaseDCMembers(i, 0);
      // chemSys_->setDCMoles(DCId,0.0); //coment if DCLowerLimit in
      // kineticControllerStep/GEM_from_MT
      chemSys_->setIsDCKinetic(DCId, true);
    }
  }
  cout << endl
       << "   numGEMPhases_  = " << setw(3) << right
       << chemSys_->getNumGEMPhases() << endl;
  cout << "   numDCs_        = " << setw(3) << right << chemSys_->getNumDCs()
       << endl;
  cout << "   numICs_        = " << setw(3) << right << chemSys_->getNumICs()
       << endl;

  // cout << "Starting with a pore solution without dissolved DCs  => all
  // microPhaseSI_ = 0" << endl; init to 0 all microPhaseSI_
  // chemSys_->setZeroMicroPhaseSI();

  bool writeICsDCs = true;
  if (writeICsDCs)
    writeTxtOutputFiles_onlyICsDCs(0); // to check the total ICs

  // variables used in DCLowerLimit computation
  double volMolDiff, molarMassDiff, vfracDiff;
  double microPhaseMassDiff, scaledMassDiff, numMolesDiff;
  // int numDCs = chemSys_->getNumDCs();
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

  double thrTimeToWriteLattice = 0.0167; // threshold ~ 1 minute

  // Main computation cycle
  for (i = 0; (i < timeSize) && (capwater); ++i) {

    TimeStruct resolvedTime = getResolvedTime(time_[i]);
    ///
    /// Do not advance the time step if GEM_run failed the last time
    ///

    bool isFirst = (i == 0) ? true : false;

    cyc++;

    // new

    if (time_[i] < beginAttackTime_) {
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
               << " (time in hours) #####" << endl;
        } else {
          lastGoodTime = time_[i - 1];
          lastGoodI = i - 1;
          cout << endl
               << endl
               << endl
               << "##### Controller::doCycle  START NEW CYCLE   "
                  "i/cyc/time_[i]/time_[i-1]/timestep: "
               << i << " / " << cyc << " / " << time_[i] << " / "
               << time_[i - 1] << " / " << timestep << " (time in hours) #####"
               << endl;
        }
      }
    } else {
      timestep = time_[i] - time_[i - 1];
      cout << endl
           << endl
           << endl
           << "##### Controller::doCycle  START NEW CYCLE - SA   "
              "i/cyc/time_[i]/time_[i-1]/timestep: "
           << i << " / " << cyc << " / " << time_[i] << " / " << time_[i - 1]
           << " / " << timestep << " (time in hours) #####" << endl;

      if (timesGEMFailed_loc > 0) {
        cout << endl
             << ">>>>> Controller::doCycle - SA - i/cyc/time_[i] : " << i
             << " / " << cyc << " / " << time_[i] << endl;
        cout << endl
             << ">>>>> \"normal exit\" (for now!) under sulfate attack "
                "conditions"
             << endl;
        // exit(0);
        bool is_Error = false;
        throw MicrostructureException(
            "Controller", "doCycle", "SA - first GEM_run failed (0)", is_Error);
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
      lattice_->writeLattice(time_[i], resolvedTime);
      lattice_->writeLatticePNG(time_[i], resolvedTime);
      if (xyz_)
        lattice_->appendXYZ(time_[i]);
      throw gex;
    }

    ///
    /// Once the change in state is determined, propagate the consequences
    /// to the 3D microstructure only if the GEM_run calculation succeeded.
    /// Otherwise we adjust the time step and try again
    ///

    if (time_[i] < beginAttackTime_) {
      if (timesGEMFailed_loc > 0) {
        cout << endl
             << "  Controller::doCycle first GEM_run failed "
                "i/cyc/time[i]/getTimesGEMFailed_loc : "
             << i << " / " << cyc << " / " << time_[i] << " / "
             << timesGEMFailed_loc << endl;

        //**********************

        cout << endl
             << "  Controller::doCycle - PRBL_0      i/cyc/time_[i]/timestep   "
                "  "
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
                //      << " / " << indFracNum << " / " << delta2Time << " / "
                //      << numGen << " / "
                //      << timeZero << " / " << minTime << endl;
                // cout << "         =>   numIncreaseInterval = " <<
                // numMaxIntervals
                //      << " (max val) => change indFracNum (next timeZero)!!!"
                //      << endl;
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
            timesGEMFailed_loc =
                calculateState(timeTemp, timestep, isFirst, cyc);
            if (timesGEMFailed_loc == 0) {
              time_[i] = timeTemp;
              cout << "  Controller::doCycle - PRBL_0 solved for "
                      "i/cyc/time_[i]/timestep/numTotGen : "
                   << i << " / " << cyc << " / " << time_[i] << " / "
                   << timestep << " / " << numTotGen << endl;
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
    } else {
      if (timesGEMFailed_loc > 0) {
        cout << endl
             << "  Controller::doCycle => SA - first GEM_run failed "
                "i/cyc/time[i]/getTimesGEMFailed_loc : "
             << i << " / " << cyc << " / " << time_[i] << " / "
             << timesGEMFailed_loc << endl;

        cout << endl
             << ">>>>> \"normal exit\" (for now!) under sulfate attack "
                "conditions"
             << endl;
        // exit(0);
        bool is_Error = false;
        throw MicrostructureException(
            "Controller", "doCycle", "SA - first GEM_run failed (1)", is_Error);
      }
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
      for (int ij = 0; ij < numSites_; ij++) {
        site_l.microPhaseId = (lattice_->getSite(ij))->getMicroPhaseId();
        site_l.growth = (lattice_->getSite(ij))->getGrowthPhases();
        site_l.wmc = (lattice_->getSite(ij))->getWmc();
        site_l.wmc0 = (lattice_->getSite(ij))->getWmc0();
        site_l.visit = 0;
        site_l.inGrowInterfacePos =
            (lattice_->getSite(ij))->getInGrowInterfacePosVector();
        site_l.inDissInterfacePos =
            (lattice_->getSite(ij))->getInDissInterfacePos();
        iniLattice.site.push_back(site_l);
      }
      iniLattice.interface.clear();
      RestoreInterface interface_l; // only one declaration
      int dimLatticeInterface =
          lattice_->getInterfaceSize(); // only one declaration
      for (int ij = 0; ij < dimLatticeInterface; ij++) {
        interface_l.microPhaseId = lattice_->getInterface(ij).getMicroPhaseId();
        interface_l.growthSites = lattice_->getInterface(ij).getGrowthSites();
        interface_l.dissolutionSites =
            lattice_->getInterface(ij).getDissolutionSites();
        iniLattice.interface.push_back(interface_l);
      }
      iniLattice.numRNGcall_0 = lattice_->getNumRNGcall_0();
      iniLattice.numRNGcallLONGMAX = lattice_->getNumRNGcallLONGMAX();
      iniLattice.lastRNG = lattice_->getLastRNG();

      int phId = 0;
      int changeLattice = -100;
      int whileCount = 0;

      std::vector<int> numSitesNotAvailable;
      numSitesNotAvailable.clear();
      std::vector<int> vectPhIdDiff;
      vectPhIdDiff.clear();
      std::vector<std::string> vectPhNameDiff;
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
        int numSitesNotAvailableSize;
        while (changeLattice == 0) { // - for many phases!
          whileCount++;
          numSitesNotAvailableSize = numSitesNotAvailable.size();
          cout << endl
               << "  Controller::doCycle - cyc = " << cyc
               << " :  changeLattice = " << changeLattice
               << "  =>  whileCount = " << whileCount << endl;

          while (timesGEMFailed_recall != 0) {

            // reset for ChemicalSystem:
            for (int ij = 0; ij < numDCs_; ij++) {
              chemSys_->setDCMoles(ij, iniLattice.DCMoles[ij]);
            }

            // reset for Lattice:
            lattice_->setCount(iniLattice.count);
            for (int ij = 0; ij < numSites_; ij++) {
              (lattice_->getSite(ij))
                  ->setMicroPhaseId(iniLattice.site[ij].microPhaseId);
              (lattice_->getSite(ij))
                  ->setGrowthPhases(iniLattice.site[ij].growth);
              (lattice_->getSite(ij))->setWmc(iniLattice.site[ij].wmc);
              (lattice_->getSite(ij))->setWmc0(iniLattice.site[ij].wmc0);
              (lattice_->getSite(ij))
                  ->setVisit(iniLattice.site[ij].visit); // or 0!
              (lattice_->getSite(ij))
                  ->setInGrowInterfacePosVector(
                      iniLattice.site[ij].inGrowInterfacePos);
              (lattice_->getSite(ij))
                  ->setInDissInterfacePos(
                      iniLattice.site[ij].inDissInterfacePos);
            }
            for (int ij = 0; ij < dimLatticeInterface; ij++) {
              lattice_->setInterfaceMicroPhaseId(
                  ij, iniLattice.interface[ij].microPhaseId); // same as before!
              lattice_->setGrowthSites(ij,
                                       iniLattice.interface[ij].growthSites);
              lattice_->setDissolutionSites(
                  ij, iniLattice.interface[ij].dissolutionSites);
            }
            lattice_->setGrowthInterfaceSize(iniLattice.growthInterfaceSize);
            lattice_->setDissolutionInterfaceSize(
                iniLattice.dissolutionInterfaceSize);
            lattice_->resetRNG(iniLattice.numRNGcall_0,
                               iniLattice.numRNGcallLONGMAX,
                               iniLattice.lastRNG);

            cout << "  Controller::doCycle - cyc = " << cyc
                 << " :  reset system OK & GEM_run recall for "
                    "i/whileCount/numSitesNotAvailable.size() = "
                 << i << " / " << whileCount << " / "
                 << numSitesNotAvailableSize << endl;
            cout << "  Controller::doCycle - cyc = " << cyc
                 << " :  reset DCLowerLimits :" << endl;

            chemSys_->initDCLowerLimit(0); // check!
            for (int ij = 0; ij < numSitesNotAvailableSize; ij++) {

              phId = vectPhIdDiff[ij];

              DCId = chemSys_->getMicroPhaseDCMembers(phId, 0);

              volMolDiff = chemSys_->getDCMolarVolume(DCId);  // m3/mol
              molarMassDiff = chemSys_->getDCMolarMass(DCId); // g/mol

              vfracDiff =
                  ((double)numSitesNotAvailable[ij]) / ((double)numSites_);

              microPhaseMassDiff =
                  vfracDiff * molarMassDiff / volMolDiff / 1.0e6; // g/cm3

              scaledMassDiff =
                  microPhaseMassDiff * 100.0 / lattice_->getInitSolidMass();

              if (chemSys_->isKinetic(phId) && (time_[i] < beginAttackTime_)) {
                kineticController_->updateKineticStep(cyc, phId, scaledMassDiff,
                                                      timestep);
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
            for (int ij = 0; ij < numSitesNotAvailableSize; ij++) {
              phId = vectPhIdDiff[ij];
              DCId = chemSys_->getMicroPhaseDCMembers(phId, 0);
              cout << "                        cyc = " << cyc << " :  #"
                   << setw(3) << right << ij << "#/ " << setw(15) << left
                   << vectPhNameDiff[ij] << "   " << setw(5) << right << phId
                   << "   " << setw(9) << right << lattice_->getCount(phId)
                   << "   " << setw(9) << right
                   << lattice_->getDissolutionInterfaceSize(phId) << "   "
                   << setw(9) << right << numSitesNotAvailable[ij] << "   "
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
              for (int iii = 0; iii < numSitesNotAvailableSize; iii++) {
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
             << "Controller::doCycle - hydration & lattice update => normal "
                "end - cyc = "
             << cyc << " (i = " << i << ")" << endl;
      }

    } catch (DataException dex) {
      dex.printException();
      lattice_->writeLattice(time_[i], resolvedTime);
      lattice_->writeLatticePNG(time_[i], resolvedTime);
      if (xyz_)
        lattice_->appendXYZ(time_[i]);
      throw dex;
    } catch (EOBException ex) {
      ex.printException();
      lattice_->writeLattice(time_[i], resolvedTime);
      lattice_->writeLatticePNG(time_[i], resolvedTime);
      if (xyz_)
        lattice_->appendXYZ(time_[i]);
      throw ex;
    } catch (MicrostructureException mex) {
      cout << endl
           << "Controller::doCycle MicroEx from Lattice::changeMicrostructure "
              "- cyc = "
           << cyc << endl;
      mex.printException();
      lattice_->writeLattice(time_[i], resolvedTime);
      lattice_->writeLatticePNG(time_[i], resolvedTime);
      if (xyz_)
        lattice_->appendXYZ(time_[i]);

      // write output .txt files
      writeTxtOutputFiles(time_[i]);
      if (writeICsDCs)
        writeTxtOutputFiles_onlyICsDCs(time_[i]);

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

    // thrTimeToWriteLattice threshold ~ 1 minute i.e 0.0167 hours
    if ((time_index < (int)(outputTime_.size())) &&
        ((time_[i] >= outputTime_[time_index]) ||
         (abs(time_[i] - outputTime_[time_index]) < thrTimeToWriteLattice))) {

      double writeTime = time_[i];
      if (abs(time_[i] - outputTime_[time_index]) < thrTimeToWriteLattice)
        writeTime = outputTime_[time_index];

      // if (verbose_)
      cout << endl
           << "Controller::doCycle - write microstructure files at time_[" << i
           << "] = " << time_[i] << ", outputTime_[" << time_index
           << "] = " << outputTime_[time_index] << ", writeTime = " << writeTime
           << endl;
      //

      lattice_->writeLattice(time_[i], resolvedTime);
      lattice_->writeLatticePNG(time_[i], resolvedTime);

      if (xyz_)
        lattice_->appendXYZ(writeTime);

      lattice_->writePoreSizeDistribution(time_[i], resolvedTime);

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

      // cout << endl
      //      << " Controller::doCycle - for sulfate attack, check conditions
      //      for "
      //         "addDissolutionSites & coordination sphere "
      //      << endl;

      if (verbose_) {
        cout << "Controller::doCycle Sulfate attack module" << endl;
        cout.flush();
      }
      std::map<int, std::vector<double>> expansion;
      expansion = lattice_->getExpansion();

      ifstream instopexp("stopexp.dat"); // check!
      if (!instopexp) {
        // if (verbose_)
        cout << endl
             << "Controller::doCycle - sulfate attack module : cyc = " << cyc
             << "   =>   expansion.size() = " << expansion.size() << endl;
      } else {
        expansion.clear();
        cout << endl
             << "Controller::doCycle - sulfate attack module : cyc = " << cyc
             << "   =>   expansion has been stopped due to the "
                "percolation of damage"
             << endl;
      }
      cout.flush();

      ///
      /// Stop FM temporarily
      /// @todo What is this?  The following if block will never be run if
      /// uncommented!
      ///

      // expansion.clear();

      if (expansion.size() > 0) { // check! : expansion.size() > 0 !!!

        // double strxx, stryy, strzz;
        std::vector<double> locEleStress;
        double locTstrength;

        // aliId => aliId_ etc
        // int aliId = chemSys_->getMicroPhaseId_SA("Alite");
        // int belId = chemSys_->getMicroPhaseId_SA("Belite");
        // int aluId = chemSys_->getMicroPhaseId_SA("Aluminate");
        // int ferId = chemSys_->getMicroPhaseId_SA("Ferrite");
        // std::vector<int> isParrotKilloh_ = chemSys_->getIsParrotKilloh();
        // int sizePK_ = isParrotKilloh_.size();
        bool notPKPhase = true;

        int oldDamageCount = 0;
        int newDamageCount = 0;
        // double poreintroduce = 0.5;

        if (verbose_) {
          cout << "Controller::doCycle Sulfate attack module writing " << endl;
          cout << "Controller::doCycle lattice at time_[" << i
               << "] = " << time_[i] << ", " << endl;
          cout << "controller::doCycle outputTime_[" << time_index
               << "] = " << outputTime_[time_index] << endl;
          cout.flush();
        }

        std::string ofileName(jobRoot_);

        ostringstream ostrT;
        ostrT << setprecision(3) << temperature_;
        std::string tempstr(ostrT.str());

        int days, hours, mins;
        double hours_dbl;
        days = floor(time_[i] / 24);
        hours_dbl = time_[i] - (days * 24);
        hours = floor(hours_dbl);
        mins = floor((hours_dbl - hours) * 60);

        ostringstream ostrD, ostrH, ostrM;
        ostrD << setfill('0') << setw(4) << days;
        std::string timestrD(ostrD.str());
        ostrH << setfill('0') << setw(2) << hours;
        std::string timestrH(ostrH.str());
        ostrM << setfill('0') << setw(2) << mins;
        std::string timestrM(ostrM.str());

        std::string timeString =
            timestrD + "d" + timestrH + "h" + timestrM + "m";

        ofileName = ofileName + "." + timeString + "." + tempstr + "_SA.img";

        ///
        /// In the sulfate attack algorithm, calculate the stress and strain
        /// distributions
        ///

        thermalstr_->setEigen();

        int expindex;
        std::vector<double> expanval;
        std::vector<int> expcoordin;
        for (std::map<int, std::vector<double>>::iterator it =
                 expansion.begin();
             it != expansion.end(); it++) {

          expindex = it->first;
          expanval = it->second;
          // std::vector<int> expcoordin =
          // lattice_->getExpansionCoordin(expindex);
          expcoordin = lattice_->getSite(expindex)->getXYZ();
          thermalstr_->setEigen(expindex, expanval[0], expanval[1], expanval[2],
                                0.0, 0.0, 0.0);
          thermalstr_->setExpansionCoord(expindex, expcoordin);
        }

        ///
        /// Calculate the stress-free strain (thermal strain) state in the
        /// microstructure, and then write the displacement field
        ///

        thermalstr_->Calc(time_[i], ofileName, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

        // thermalstr_ -> writeStress(jobRoot_,time_[i],0); //write strxx
        // thermalstr_ -> writeStrainEngy(jobRoot_,time_[i]);

        // thermalstr_->writeDisp(jobRoot_, time_[i]);
        thermalstr_->writeDisp(jobRoot_, timeString);

        ///
        /// Get the true volume of each voxel after FEM calculation
        ///

        //// truevolume not used!!!
        // double otruevolume = 0.0;
        // double truevolume = 0.0;
        // for (int ii = 0; ii < numSites_; ii++) {
        //   Site *ste;
        //   ste = lattice_->getSite(i); // check! i -> ii !!!
        //   otruevolume = ste->getTrueVolume();
        //   for (int j = 0; j < 3; j++) {
        //     truevolume += thermalstr_->getEleStrain(ii, j);
        //   }
        //   truevolume = otruevolume * (1 + truevolume);
        //   ste->setTrueVolume(truevolume);
        // }

        // double dwmcval = poreintroduce;
        double poreincrease = 0.2;
        double damageexp_0 = 1.0 / 3.0 * poreincrease;
        std::vector<double> damageexp(3, damageexp_0);
        std::vector<double> damageexpo;
        Site *ste; // , *stenb;
        int pid;

        for (int index = 0; index < numSites_; index++) {

          ste = lattice_->getSite(index);
          pid = ste->getMicroPhaseId();

          if ((ste->IsDamage())) {

            oldDamageCount++;
            // cout << endl << "SA-test: oldDamageCount_ = " << oldDamageCount_
            //      << "  index = " << index << "  pid = " << pid << endl;

            // double strxx, stryy, strzz;
            // strxx = stryy = strzz = 0.0;
            // strxx = thermalstr_->getEleStress(index, 0);
            // stryy = thermalstr_->getEleStress(index, 1);
            // strzz = thermalstr_->getEleStress(index, 2);
            locEleStress = thermalstr_->getEleStressMod(index);
            // strxx = locEleStress[0];
            // stryy = locEleStress[1];
            // strzz = locEleStress[2];
            // if ((strxx >= 1.0) || (stryy >= 1.0) || (strzz >= 1.0)) {
            if ((locEleStress[0] >= 1.0) || (locEleStress[1] >= 1.0) ||
                (locEleStress[2] >= 1.0)) {
              damageexpo = lattice_->getExpansion(index);
              damageexpo[0] += damageexp_0;
              damageexpo[1] += damageexp_0;
              damageexpo[2] += damageexp_0;
              lattice_->setExpansion(index, damageexpo);

              ///
              /// Set expansion site to be damaged if there is one, as
              /// determined by the setEigen function returning every site above
              /// damage stress threshold
              ///

              // lattice_->dWaterChange(poreintroduce);

              // lattice_->setWmc0(index,dwmcval); // *
              // lattice_->dWmc(index, dwmcval);
              // for (int j = 0; j < NN_NNN; j++) { //ste->nbSize(2)
              //   stenb = ste->nb(j);
              //   stenb->dWmc(dwmcval);
              // }
              // need interface update???

              // lattice_->dWaterChange(poreincrease);
              /// JWB: This next line must be from some earlier version
              /// The called method does not exist any longer
              ///
              /// @todo: Determine whether it is necessary to add this back
              /// in for crystallization pressure calculations
              //
              // ste->setVolume(VOIDID,(ste->getVolume(VOIDID) + poreincrease));
              //
            }
          } else {

            ///
            /// The next block gets the stress in each voxel that does NOT
            /// contain a clinker phase (C3S, C2S, C3A, or C4AF), then determine
            /// if the voxel should be damaged as a result

            /// Prefer to make this independent of whether or not there is C4AF
            /// in the phase definitions.  What if this is a white cement or
            /// something?
            ///
            /// @note Associating the last clinker phase with id 5 is a kluge
            /// @todo Give each phase a calcstress property or something like
            /// that
            ///       that can be checked instead of hardwiring phase ids

            notPKPhase = true;
            for (int i = 0; i < sizePK_; i++) {
              if (pid == isParrotKilloh_[i]) {
                notPKPhase = false;
                break;
              }
            }
            if ((pid > ELECTROLYTEID) && notPKPhase &&
                (chemSys_->isPorous(pid) || chemSys_->isWeak(pid))) {
              // double strxx, stryy, strzz;
              // strxx = stryy = strzz = 0.0;
              // strxx = thermalstr_->getEleStress(index, 0);
              // stryy = thermalstr_->getEleStress(index, 1);
              // strzz = thermalstr_->getEleStress(index, 2);
              locEleStress = thermalstr_->getEleStressMod(index);
              // strxx = locEleStress[0];
              // stryy = locEleStress[1];
              // strzz = locEleStress[2];
              // if ((strxx >= thermalstr_->getTstrength(index)) ||
              //     (stryy >= thermalstr_->getTstrength(index)) ||
              //     (strzz >= thermalstr_->getTstrength(index))) {
              locTstrength = thermalstr_->getTstrength(index);
              // if ((strxx >= locTstrength) ||
              //     (stryy >= locTstrength) ||
              //     (strzz >= locTstrength)) {
              if ((locEleStress[0] >= locTstrength) ||
                  (locEleStress[1] >= locTstrength) ||
                  (locEleStress[2] >= locTstrength)) {

                newDamageCount++;
                ste->setDamage();
                lattice_->setExpansion(index, damageexp);
                // lattice_->dWaterChange(poreintroduce);

                // lattice_->setWmc0(index,dwmcval); // *
                // lattice_->dWmc(index, dwmcval);
                // for (int j = 0; j < NN_NNN; j++) {
                //   stenb = ste->nb(j);
                //   stenb->dWmc(dwmcval);
                //   if ((stenb->getWmc() > 0.0) &&
                //       (stenb->getMicroPhaseId() > ELECTROLYTEID)) {
                //     lattice_->addDissolutionSite(stenb,
                //                                  stenb->getMicroPhaseId());
                //   }
                // }

                // std::vector<double> damageexp;
                // damageexp.clear();
                // double poreindamage = 0.6;
                // damageexp.resize(3,(1.0 / 3.0 * poreindamage));
                // lattice_->setExpansion(index,damageexp);
                // std::vector<int> coordin;
                // coordin.clear();
                // coordin.resize(3,0);
                // coordin[0] = ste->getX();
                // coordin[1] = ste->getY();
                // coordin[2] = ste->getZ();
                // lattice_->setExpansionCoordin(index,coordin);
                // lattice_->dWaterChange(poreindamage);
                // ste->setVolume(VOIDID,poreindamage);
                // }
              }
            }
          }
        } // End of loop over all voxels
        if (oldDamageCount_ != oldDamageCount) {
          cout << endl
               << "Controller::doCycle SA -  error : oldDamageCount_ != "
                  "oldDamageCount"
                  " <-> cyc/oldDamageCount_/oldDamageCount = "
               << cyc << " / " << oldDamageCount_ << " / " << oldDamageCount
               << endl;
          cout << endl
               << "         newDamageCount/allDamageCount_ = " << newDamageCount
               << " / " << allDamageCount_ << endl;
          cout << endl << "exit" << endl;
          exit(0);
        } else {
          cout << endl
               << "Controller::doCycle SA - "
                  "cyc/newDamageCount_/oldDamageCount_/allDamageCount_ = "
               << cyc << " / " << newDamageCount << " / " << oldDamageCount_;

          allDamageCount_ = newDamageCount + oldDamageCount_;
          oldDamageCount_ = allDamageCount_;

          cout << " / " << allDamageCount_ << endl;
          cout.flush();
        }

        // ofstream outdamage("damage.dat");
        // outdamage.close();

        lattice_->writeDamageLattice(time_[i], resolvedTime);
        lattice_->writeDamageLatticePNG(time_[i], resolvedTime);
        // to see whether new damage is generated
      }

      cout << endl
           << "Controller::doCycle - sulfate attack module => normal end - cyc "
              "= "
           << cyc << " (i = " << i << ")" << endl;
    }
  }

  ///
  /// Write the final lattice state to an ASCII file and to a PNG file for
  /// visualization
  ///

  TimeStruct timebefore = getResolvedTime(time_[i - 1]);
  lattice_->writeLattice(time_[i - 1], timebefore);
  lattice_->writeLatticePNG(time_[i - 1], timebefore);

  // if (xyz_ && (time_[i - 1] < sattack_time_))... ?
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

    // std::vector<double> impurityrelease;
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

    // chemSys_->setMicroPhaseSI();

    // chemSys_->initDCLowerLimit(0); // check!
    kineticController_->calculateKineticStep(time, dt, cyc); // check!

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
  int i;

  // Output to files the solution composition data, phase data, DC data,
  // microstructure data, pH, and C-S-H composition and Ca/Si ratio

  std::string outfilename = jobRoot_ + "_Solution.csv";
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
        std::string dcname = chemSys_->getDCName(i);
        double V0 =
            chemSys_->getDCMoles(dcname) * chemSys_->getDCMolarVolume(dcname);
        outfs << "," << V0;
      }
    } else {
      std::string msg = "Divide by zero error for DC " + chemSys_->getDCName(i);
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
  std::vector<double> ICMoles;
  ICMoles.resize(numICs_, 0.0);
  std::vector<double> DCMoles;
  DCMoles.resize(numDCs_, 0.0);
  for (i = 0; i < numDCs_; i++) {
    DCMoles[i] = chemSys_->getDCMoles(i);
  }

  std::string outfilenameIC = jobRoot_ + "_icmoles.csv";
  std::string outfilenameDC = jobRoot_ + "_dcmoles.csv";
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

  std::vector<int> impurityDCID;
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

void Controller::parseDoc(const std::string &docName) {

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

    // There may be times associated with chemical attack
    // Next three blocks search for this
    cdi = it.value().find("beginttacktime");
    if (cdi != it.value().end()) {
      beginAttackTime_ = cdi.value();
      beginAttackTime_ *= (H_PER_DAY);
    }

    // Input times are conventionally in days
    // Immediately convert to hours within model
    cdi = it.value().find("endattacktime");
    if (cdi != it.value().end()) {
      endAttackTime_ = cdi.value();
      endAttackTime_ *= (H_PER_DAY);
    }

    // Input times are conventionally in days
    // Immediately convert to hours within model
    cdi = it.value().find("attacktimeinterval");
    if (cdi != it.value().end()) {
      attackTimeInterval_ = cdi.value();
      attackTimeInterval_ *= (H_PER_DAY);
    }

    // Done searching for chemical attack times
    attack_ = false;
    bool errorAttack = false;
    if (simType_ == SULFATE_ATTACK || simType_ == LEACHING) {
      if (beginAttackTime_ >= 0) {
        if (endAttackTime_ > beginAttackTime_) {
          if (attackTimeInterval_ > 0) {
            attack_ = true;
          } else {
            errorAttack = true;
          }
        } else {
          errorAttack = true;
        }
      } else {
        errorAttack = true;
      }
    }

    if (attack_ && errorAttack) {
      cout << endl << endl << "************" << endl;
      cout << endl
           << "=> you decided to simulate a leaching or a sulfate attack "
              "but your time parameters are not set accordingly!"
           << endl;
      cout << endl
           << "=> to do it, the specific controll parameters "
              "must fulfill some additional conditions:"
           << endl;
      cout << "  -> beginattacktime >= 0" << endl;
      cout << "  -> endattacktime > beginattacktime" << endl;
      cout << "  -> attacktimeinterval > 0" << endl;

      cout << endl
           << "=> by default all these variables are set to -1.0" << endl;
      cout << endl << "=> for the current simulation their values are:" << endl;
      cout << "  -> beginattacktime = " << beginAttackTime_ << endl;
      cout << "  -> endattacktime = " << endAttackTime_ << endl;
      cout << "  -> attacktimeinterval = " << attackTimeInterval_ << endl;

      cout << endl
           << "=> before to restart the program, please modify their values "
              "into parameters.json file)"
           << endl;

      cout << endl << endl << "STOP" << endl;
      cout << endl << endl << "************" << endl;
      throw DataException("Controller", "Controller",
                          "leaching or sulfate attack time parameters setting");

    } else if (attack_) {

      cout << endl << "   => you decided to simulate a ";
      if (simType_ == LEACHING) {
        leachTime_ = beginAttackTime_; // not default leachTime_ = 1.0e10
        cout << "leaching ";
      } else if (simType_ == SULFATE_ATTACK) {
        sulfateAttackTime_ =
            beginAttackTime_; // not default sulfateAttackTime_ = 1.0e10
        cout << "sulfate attack ";
      }

      cout << "using these time parameters (in days):" << endl;
      cout << "     -> beginattacktime = " << setw(5) << right
           << (int)beginAttackTime_ << endl;
      cout << "     -> endattacktime = " << setw(5) << right
           << (int)endAttackTime_ << endl;
      cout << "     -> attacktimeinterval = " << setw(5) << right
           << (int)attackTimeInterval_ << endl;

      double tp = beginAttackTime_;
      int tempSize;
      int i, j;

      tempSize = outputTime_.size() - 1;
      for (i = 0; i < tempSize; i++) {
        if (outputTime_[i] > tp) {
          int last = tempSize;
          outputTime_.erase(outputTime_.begin() + i,
                            outputTime_.begin() + last);
          break;
        }
      }

      tempSize = time_.size() - 1;
      for (i = 0; i < tempSize; i++) {
        if (time_[i] > tp) {
          int last = tempSize;
          time_.erase(time_.begin() + i, time_.begin() + last);
          break;
        }
      }

      while (tp <= endAttackTime_) {
        // time_.push_back(tp);
        outputTime_.push_back(tp);
        tp += attackTimeInterval_;
      }

      int time_Size = outputTime_.size();
      for (i = 0; i < time_Size - 1; i++) {
        for (j = i + 1; j < time_Size; j++) {
          if (outputTime_[i] > outputTime_[j]) {
            tp = outputTime_[i];
            outputTime_[i] = outputTime_[j];
            outputTime_[j] = tp;
          }
        }
      }

      time_Size = outputTime_.size() - 1;
      for (i = 0; i < time_Size; i++) {
        if (abs(outputTime_[i] - outputTime_[i + 1]) <= 1.0e-9) {
          outputTime_.erase(outputTime_.begin() + i);
        }
      }

    } else if (simType_ == HYDRATION) {
      cout << endl << "   => you decided to simulate a hydration" << endl;
      beginAttackTime_ = 1.e10;
      endAttackTime_ = 1.e10;
      attackTimeInterval_ = 1.e10;
    }

    //
    // Now populate the calculation times on a natural log scale
    // and fold in the output times in order
    time_.clear();
    testTime = 0.0;
    int j = 0;
    bool done = false;
    int outputTimeSize = static_cast<int>(outputTime_.size());
    while (testTime < finalTime) {
      testTime += (0.1 * (testTime + 0.024));
      done = false;
      if (j < outputTimeSize) {
        while (j < outputTimeSize && !done) {
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
