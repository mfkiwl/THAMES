/**
@file Controller.cc
@brief Definition of Controller class methods
*/

#include "Controller.h"

Controller::Controller(Lattice *msh, KineticController *kc, ChemicalSystem *cs,
                       Solution *solut, ThermalStrain *thmstr,
                       const int simtype, const string &parfilename,
                       const string &jobname, const bool verbose,
                       const bool warning)
    : lattice_(msh), kineticController_(kc), chemSys_(cs), solut_(solut),
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
    ofstream out(outfilename.c_str(), ios::app);
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
    ofstream out1(outfilename.c_str(), ios::app);
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
    ofstream out2(outfilename.c_str(), ios::app);
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
    ofstream out3(outfilename.c_str(), ios::app);
    if (!out3) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    out3 << "Time(d),pH" << endl;
    out3.close();

    outfilename = jobroot_ + "_CSH.csv";
    ofstream out4(outfilename.c_str(), ios::app);
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
    ofstream out5(outfilename.c_str(), ios::app);
    if (!out5) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    out5 << "Time(d),C/S in solid" << endl;
    out5.close();

    outfilename = jobroot_ + "_icmoles.csv";
    ofstream out6(outfilename.c_str(), ios::app);
    if (!out6) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    out6 << "Time(d)";
    for (int i = 0; i < chemSys_->getNumICs(); i++) {
      out6 << "," << chemSys_->getICName(i);
    }
    out6 << endl;
    out6.close();

    outfilename = jobroot_ + "_Enthalpy.csv";
    ofstream out7(outfilename.c_str(), ios::app);
    if (!out7) {
      throw FileException("Controller", "Controller", outfilename,
                          "Could not append");
    }
    out7 << "Time(d),Enthalpy(J)" << endl;
    out7.close();

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
  double next_stat_time = statfreq_;

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

#ifdef DEBUG
  cout << "Controller::doCycle Entering Lattice::findInterfaces" << endl;
  cout.flush();
#endif

  lattice_->findInterfaces();

#ifdef DEBUG
  cout << "Controller::doCycle Returned from Lattice::findInterfaces" << endl;
#endif

  ///
  /// The next for loop is the main computation cycle loop, iterating over
  /// the array of cycle times and executing all the tasks of
  /// a computational cycle.
  ///
  /// If the time is greater than or equal to one of the output times just
  /// defined, then output the microstructure to an ASCII file and a PNG file
  ///

  double timestep = 0.0;
  bool capwater = true; // True if some capillary water is available
  time_index = 0;

  for (i = 0; (i < time_.size()) && (capwater); ++i) {

    ///
    /// Do not advance the time step if GEM_run failed the last time
    ///

    bool isFirst = (i == 0) ? true : false;

    if (chemSys_->getTimesGEMFailed() > 0)
      i -= 1;

    cout << "Time = " << time_[i] << endl;
    if (time_index < output_time_.size()) {
#ifdef DEBUG
      cout << "Controller::doCycle Next output time = "
           << output_time_[time_index] << endl;
      cout.flush();
#endif
    } else {
      int lasttime = time_.size() - 1;
#ifdef DEBUG
      cout << "Controller::doCycle Next output time = " << time_[lasttime]
           << endl;
      cout.flush();
#endif
    }

    time_t lt10 = time(NULL);
    struct tm *time10;
    time10 = localtime(&lt10);
    cout << asctime(time10);

    timestep = (i > 0) ? (time_[i] - time_[i - 1]) : (time_[i]);

    ///
    /// Assume that only capillary pore water is chemically reactive,
    /// while water in nanopores is chemically inert.
    ///
    ///
    /// This is the main step of the cycle; the calculateState method
    /// runs all the major steps of a computational cycle
    ///

#ifdef DEBUG
    cout << "Controller::doCycle Entering Controller::calculateState "
         << "with isFirst = " << isFirst << endl;
    cout.flush();
#endif
    try {
      calculateState(time_[i], timestep, isFirst);
    } catch (GEMException gex) {
      lattice_->writeLattice(time_[i], sim_type_, jobroot_);
      lattice_->writeLatticePNG(time_[i], sim_type_, jobroot_);
      throw gex;
    }

#ifdef DEBUG
    cout << "Controller::doCycle Returned from Controller::calculateState("
         << time_[i] << "," << timestep << "," << isFirst << ")" << endl;
    cout.flush();
#endif

    ///
    /// Once the change in state is determined, propagate the consequences
    /// to the 3D microstructure only if the GEM_run calculation succeeded.
    /// Otherwise we will need to tweak the IC moles so just return from
    /// function without doing anything else
    ///

    if (chemSys_->getTimesGEMFailed() > 0) {
      // Skip the remainder of this iteration and go to the next iteration
      if (warning_) {
        cout << "Controller::doCycle  WARNING: Previous call to "
             << "GEM_run failed, so I will" << endl
             << "Controller::doCycle  not update the microstructure "
             << "or do anything else" << endl
             << "controller::doCycle  during this time step" << endl;
        cout.flush();
      }
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
    ////

    try {
      lattice_->changeMicrostructure(time_[i], sim_type_, isFirst, capwater);
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
#ifdef DEBUG
      cout << "Controller::doCycle Returned from writing lattice" << endl;
#endif
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

#ifdef DEBUG
      cout << "Controller::doCycle expansion.size() = " << expansion.size()
           << endl;
#endif
      if (expansion.size() > 1) {

#ifdef DEBUG
        cout << "Controller::doCycle time_ is: " << time_[i]
             << ". expansion.size() is: " << expansion.size() << endl;
#endif

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

#ifdef DEBUG
        cout << "Controller::doCycle sulfate attack module entering "
             << "ThermalStrain:Calc" << endl;
        cout.flush();
#endif

        thermalstr_->Calc(time_[i], ofileName, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

#ifdef DEBUG
        cout << "Controller::doCycle sulfate attack module returned from "
             << "ThermalStrain:Calc" << endl;
        cout.flush();
#endif

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
                for (int j = 0; j < ste->nbSize(1); j++) {
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

void Controller::calculateState(double time, double dt, bool isFirst) {
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

    vector<double> impurityrelease;
    impurityrelease.clear();
    impurityrelease.resize(chemSys_->getNumMicroImpurities(), 0.0);

    ///
    /// Get the number of moles of each IC dissolved from kinetically controlled
    /// phases
    ///

    double T = lattice_->getTemperature();
#ifdef DEBUG
    cout << "Controller::calculateState Entering "
            "KineticController::calculateKineticStep"
         << endl;
    cout.flush();
#endif

    kineticController_->calculateKineticStep(dt, T, isFirst);

#ifdef DEBUG
    cout << "Controller::calculateState Returned from "
            "KineticController::calculateKineticStep"
         << endl;
    cout.flush();
#endif

    ///
    /// The next block only operates for sulfate attack iterations
    /// It determines how many IC moles of hydrogen, oxygen, sulfur, etc to add
    ///

    if (time >= sattack_time_) {

#ifdef DEBUG
      cout << "Controller::calculateState waterchange_ in Lattice is: "
           << lattice_->getWaterchange() << endl;
      cout.flush();
#endif

      double addwatervol = lattice_->getWaterchange() /
                           lattice_->getNumsites() *
                           chemSys_->getInitMicroVolume();

      ///
      /// Get the molar volume of water from the GEM node
      ///

      double water_v0 = chemSys_->getNode()->DC_V0(
          chemSys_->getMicroPhaseDCMembers(ELECTROLYTEID, 0), chemSys_->getP(),
          chemSys_->getTemperature());
      double addwatermol = addwatervol / water_v0;

#ifdef DEBUG
      if (verbose_) {
        cout << "Controller::calculateState Molar volume of water is: "
             << water_v0 << endl;
        cout << "Controller::calculateState The moles of water added into the "
                "system are: "
             << addwatermol << endl;
        cout.flush();
      }
#endif
      for (int i = 0; i < chemSys_->getNumICs(); i++) {
        if (chemSys_->getICName(i) == "H") {
#ifdef DEBUG
          cout << "Controller::calculateState previous IC moles for H is: "
               << chemSys_->getICMoles(i) << endl;
          cout.flush();
#endif
          chemSys_->setICMoles(i, (chemSys_->getICMoles(i) + 2 * addwatermol));
#ifdef DEBUG
          cout << "Controller::calculateState new ICmoles for H is: "
               << chemSys_->getICMoles(i) << endl;
          cout.flush();
#endif
        }
        if (chemSys_->getICName(i) == "O") {
#ifdef DEBUG
          cout << "Controller::calculateState previous IC moles for O is: "
               << chemSys_->getICMoles(i) << endl;
          cout.flush();
#endif
          chemSys_->setICMoles(i, (chemSys_->getICMoles(i) + addwatermol));
#ifdef DEBUG
          cout << "Controller::calculateState new ICmoles for O is: "
               << chemSys_->getICMoles(i) << endl;
          cout.flush();
#endif
        }
      }
    }

    ///
    /// Now that the method is done determining the change in moles of each IC,
    /// launch a thermodynamic calculation to determine new equilibrium state
    ///
    /// The `ChemicalSystem` object provides an interface for these calculations
    ///

    int timesGEMFailed = 0;

    if (verbose_)
      cout << "Going to launch thermodynamic calculation now with isFirst = "
           << isFirst << endl;
    ;
    try {
      timesGEMFailed = chemSys_->calculateState(time, isFirst);
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

    if (chemSys_->getTimesGEMFailed() > 0) {
      return;
    }

    try {
      double aveSI = 0.0;
      double moles = 0.0;
      vector<int> microPhaseMembers;
      for (int i = 0; i < chemSys_->getNumMicroPhases(); ++i) {
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
          aveSI = aveSI / (double)(microPhaseMembers.size());
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

    ///
    /// Set the kinetic DC moles.  This adds the clinker components to the DC
    /// moles.
    ///
    /// @todo Find out what this is and why it needs to be done
    ///

#ifdef DEBUG
    cout << "Controller::calculateState Entering setKineticDCMoles" << endl;
    cout.flush();
#endif
    kineticController_->setKineticDCMoles();
#ifdef DEBUG
    cout << "Controller::calculateState Returned from setKineticDCMoles"
         << endl;
    cout.flush();
#endif

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
    for (int i = 0; i < chemSys_->getNumDCs(); i++) {
      cc = chemSys_->getDCClassCode(i);
      if (cc == 'S' || cc == 'T' || cc == 'W') {
        out3 << "," << (chemSys_->getNode())->Get_cDC((long int)i);
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
#ifdef DEBUG
    cout << "Controller::calculateState Writing DC volumes file at time = "
         << time << endl;
    cout.flush();
#endif
    for (int i = 0; i < chemSys_->getNumDCs(); i++) {
      if (chemSys_->getDCMolarMass(i) > 0.0) {
        cc = chemSys_->getDCClassCode(i);
        if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M' || cc == 'W') {
          string dcname = chemSys_->getDCName(i);
          double V0 =
              chemSys_->getDCMoles(dcname) * chemSys_->getDCMolarVolume(dcname);
          out4 << "," << V0;
#ifdef DEBUG
          cout << "Controller::calculateState    DC = "
               << chemSys_->getDCName(i)
               << ", moles = " << chemSys_->getDCMoles(i)
               << ", molar mass = " << chemSys_->getDCMolarMass(i) << endl;
          cout.flush();
#endif
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

#ifdef DEBUG
    cout << "Controller::calculateState Writing microstructure "
         << "volume fractions at time " << time << endl;
    cout.flush();
#endif

    out5 << setprecision(5) << time;
    for (int i = 0; i < chemSys_->getNumMicroPhases(); i++) {
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
#ifdef DEBUG
    cout << "Controller::calculateState Writing pH values";
    cout.flush();
#endif
    out6 << setprecision(5) << time;
    out6 << "," << (chemSys_->getPH()) << endl;
    out6.close();

#ifdef DEBUG
    cout << "Controller::calculateState Entering "
            "ChemicalSystem::setGEMPhaseStoich"
         << endl;
    cout.flush();
#endif

    chemSys_->setGEMPhaseStoich();

#ifdef DEBUG
    cout << "Controller::calculateState Returned from "
            "ChemicalSystem::setGEMPhaseStoich"
         << endl;
    cout.flush();
#endif

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
    double CaMoles = 0.0, SiMoles = 0.0, CaSiRatio = 0.0;
    outfilename = jobroot_ + "_CSH.csv";
    ofstream out7(outfilename.c_str(), ios::app);
    if (!out7) {
      throw FileException("Controller", "calculateState", outfilename,
                          "Could not append");
    }
    out7 << setprecision(5) << time;
    for (unsigned int i = 0; i < chemSys_->getNumICs(); i++) {
      out7 << "," << CSHcomp[i];
      if (chemSys_->getICName(i) == "Ca") {
        CaMoles = CSHcomp[i];
      }
      if (chemSys_->getICName(i) == "Si") {
        SiMoles = CSHcomp[i];
      }
    }
    if (CaMoles < 1.0e-16)
      CaMoles = 1.0e-16;
    if (SiMoles < 1.0e-16)
      SiMoles = 1.0e-16;
    CaSiRatio = CaMoles / SiMoles;
    out7 << "," << CaSiRatio << endl;
    out7.close();

    chemSys_->setGEMPhaseStoich();
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
    for (int i = 0; i < chemSys_->getNumGEMPhases(); i++) {
      cc = chemSys_->getGEMPhaseClassCode(i);
      if (cc == 's') {
        phaseRecord = chemSys_->getPGEMPhaseStoich(i);
        ICIndex = chemSys_->getICId("Ca");
        CaMoles += phaseRecord[ICIndex];
        ICIndex = chemSys_->getICId("Si");
        SiMoles += phaseRecord[ICIndex];
      }
    }
    if (SiMoles != 0) {
      CaSiRatio = CaMoles / SiMoles;
      out8 << "," << CaSiRatio << endl;
    } else {
      out8 << ",Si_moles is ZERO" << endl;
    }
    out8.close();

    outfilename = jobroot_ + "_icmoles.csv";
    ofstream out9(outfilename.c_str(), ios::app);
    if (!out9) {
      throw FileException("Controller", "calculateState", outfilename,
                          "Could not append");
    }
    out9 << setprecision(5) << time;
    for (int i = 0; i < chemSys_->getNumICs(); i++) {
      out9 << "," << chemSys_->getICMoles(i);
    }
    out9 << endl;
    out9.close();

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
    for (int i = 0; i < chemSys_->getNumDCs(); i++) {
      enth += (chemSys_->getDCEnthalpy(i));
    }

    out10 << setprecision(5) << time;
    out10 << "," << enth << endl;
    if (verbose_) {
      cout << "Done!" << endl;
      cout.flush();
    }
    out10.close();

    ///
    /// Now that the end of the iteration is reached, zero out the kinetic
    /// DC moles in preparation for the next iteration
    ///
    /// @todo Find out what this is doing and why it is needed
    ///

    kineticController_->zeroKineticDCMoles();
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  } catch (FloatException flex) {
    flex.printException();
    exit(1);
  }
  return;
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
