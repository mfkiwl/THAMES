/**
@file Controller.cc
@brief Definition of Controller class methods
*/

#include "Controller.h"

Controller::Controller (Lattice *msh,
                        KineticModel *km,
                        ChemicalSystem *cs,
                        Solution *solut,
                        ThermalStrain *thmstr,
                        const int simtype,
                        const string &parfilename,
                        const string &jobname,
                        const bool verbose,
                        const bool warning,
                        const bool debug)
    :lattice_(msh),
     kineticmodel_(km),
     chemsys_(cs),
     solut_(solut),
     sim_type_(simtype),
     thermalstr_(thmstr),
     jobroot_(jobname),
     verbose_(verbose),
     warning_(warning),
     debug_(debug)
{
  unsigned int i;
  double tvalue,pvalue;
  string buff;
  vector<double> phases;
  const string imgfreqstr = "Image_frequency:";
  const string outtimestr = "OutTime:";
  const string calctimestr = "CalcTime:";

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

  chemsys_ = lattice_->getChemsys();
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
    ofstream out(outfilename.c_str(),ios::app);
    if (!out) {
        throw FileException("Controller","Controller",outfilename,"Could not append");
    }
    char cc;
    out << "Time(d)";
    for (int i = 0; i < chemsys_->getDCnum(); i++) {
        cc = chemsys_->getDCclasscode(i); 
        if (cc == 'S' || cc == 'T' || cc == 'W') {
            out << "," << chemsys_->getDCname(i);
        }
    }  
    out << endl;
    out.close();

    outfilename = jobroot_ + "_DCVolumes.csv";
    ofstream out1(outfilename.c_str(),ios::app);
    if (!out1) {
        throw FileException("Controller","Controller",outfilename,"Could not append");
    }

    out1 << "Time(d)";
    for (int i = 0; i < chemsys_->getDCnum(); i++) {
        cc = chemsys_->getDCclasscode(i); 
        if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M' || cc == 'W') {
            out1 << "," << chemsys_->getDCname(i);
        }
    }   
    out1 << endl;
    out1.close();

    outfilename = jobroot_ + "_Microstructure.csv";
    ofstream out2(outfilename.c_str(),ios::app);
    if (!out2) {
        throw FileException("Controller","Controller",outfilename,"Could not append");
    }

    out2 << "Time(d)";
    for (int i = 0; i < chemsys_->getMicphasenum(); i++) {
        out2 << "," << chemsys_->getMicphasename(i);
    }   
    out2 << endl;
    out2.close();

    outfilename = jobroot_ + "_pH.csv";
    ofstream out3(outfilename.c_str(),ios::app);
    if (!out3) {
        throw FileException("Controller","Controller",outfilename,"Could not append");
    }
    out3 << "Time(d),pH" << endl;
    out3.close();

    outfilename = jobroot_ + "_CSH.csv";
    ofstream out4(outfilename.c_str(),ios::app);
    if (!out4) {
        throw FileException("Controller","Controller",outfilename,"Could not append");
    }
    out4 << "Time(d)";
    for (int i = 0; i < chemsys_->getICnum(); i++){
        out4 << "," << chemsys_->getICname(i);
    }
    out4 << ",Ca/Si" << endl;
    out4.close();

    outfilename = jobroot_ + "_CSratio_solid.csv";
    ofstream out5(outfilename.c_str(),ios::app);
    if(!out5) {
        throw FileException("Controller","Controller",outfilename,"Could not append");
    }
    out5 << "Time(d),C/S in solid" << endl;
    out5.close();
    
    outfilename = jobroot_ + "_icmoles.csv";
    ofstream out6(outfilename.c_str(),ios::app);
    if(!out6) {
        throw FileException("Controller","Controller",outfilename,"Could not append"); 
    }
    out6 << "Time(d)";
    for (int i = 0; i < chemsys_->getICnum(); i++){
        out6 << "," << chemsys_->getICname(i);
    }
    out6 << endl;
    out6.close();
  }
  catch (FileException fex) {
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
        throw FileException("Controller","Controller",parfilename,"NOT XML FORMAT");
    }
  }
  catch (FileException fex) {
    throw fex;
  }
}

void Controller::doCycle (const string &statfilename,
                          int choice)
{    
  unsigned int i;
  int time_index;
  double next_stat_time = statfreq_;
    
  ///
  /// This block arbitrarily sets the leaching initiation time to 100 days if the leaching
  /// module is to be run, or sets the sulfate attack initiation time to 100 days if the
  /// sulfate attack module is to be run
  ///
  /// @todo Think about generalizing this more, or allowing combinations of more than one
  ///
 
  if (choice == LEACHING) {
    leach_time_ = 100.0;
  } else if (choice == SULFATE_ATTACK) {
    sattack_time_ = 100.0;
  } 

  /*
  kineticmodel_->setSattack_time(sattack_time_);
  kineticmodel_->setLeach_time(leach_time_);
  */

  chemsys_->setSattack_time(sattack_time_);
  chemsys_->setLeach_time(leach_time_);
  lattice_->setSattack_time(sattack_time_);
  lattice_->setLeach_time(leach_time_);
    
  // Initialize the list of all interfaces in the lattice

  if (verbose_) cout << "Going into Lattice::FindInterfaces..." << endl;
  lattice_->findInterfaces();
  if (verbose_) cout << "...Done!" << endl;
    
  ///
  /// The next for loop is the main computation cycle loop, iterating over
  /// the array of cycle times and executing all the tasks of
  /// a computational cycle.
  ///
  /// If the time is greater than or equal to one of the output times just
  /// defined, then output the microstructure to an ASCII file and a PNG file
  ///

  double timestep = 0.0;
  bool capwater = true;      // True if some capillary water is available
  time_index = 0;

  for (i = 0; (i < time_.size()) && (capwater); ++i) {

    cout << "Time = " << time_[i] << endl;
    cout << "Next output time = " << output_time_[time_index] << endl;

    time_t lt10 = time(NULL);
    struct tm *time10;
    time10 = localtime(&lt10);
    cout << asctime(time10);

    timestep = (i > 0) ? (time_[i] - time_[i-1]) : (time_[i]);
    bool isfirst = (i == 0) ? true : false;

    /// 
    /// Assume that only capillary pore water is chemically reactive,
    /// while water in nanopores is chemically inert.
    ///
    /// @todo the following is a kluge to force the simulation to
    /// stop when the capillary pore water is exhausted.
    ///
   
    // if (isfirst) {
    //     cout << "Going into Controller::initializeState" << endl;
    //     initializeState(time_[i]);
    // }

    ///
    /// This is the main step of the cycle; the calculateState method
    /// runs all fo the major steps of a computational cycle
    ///

    if (verbose_) cout << "Going into Controller::calculateState with isfirst = "
                       << isfirst << endl;
    try {
        calculateState(time_[i],timestep,isfirst);
    }
    catch (GEMException gex) {
        lattice_->writeLattice(time_[i],sim_type_,jobroot_);
        lattice_->writeLatticePNG(time_[i],sim_type_,jobroot_);
        throw gex;
    }

    if (verbose_) {
        cout << "*Returned from Controller::calculateState(" << time_[i] << ","
             << timestep << "," << isfirst << ")" << endl;
        cout << "*called by Controller::doCycle" << endl;
        cout.flush();
    }

    ///
    /// Once the change in state is determined, propagate the consequences
    /// to the 3D microstructure
    ///

    if (verbose_) {
        cout << "Going into Lattice::changeMicrostructure" << endl;
        cout.flush();
    }

    ///
    /// Next function can encounter EOB exceptions within but they
    /// are caught there and the program will then exit from within
    /// this function rather than throwing an exception itself
    ////
    
    try {
        lattice_->changeMicrostructure(time_[i],sim_type_,isfirst,capwater);
    }
    catch (DataException dex) {
        lattice_->writeLattice(time_[i],sim_type_,jobroot_);
        lattice_->writeLatticePNG(time_[i],sim_type_,jobroot_);
        throw dex;
    }
    catch (EOBException ex) {
        lattice_->writeLattice(time_[i],sim_type_,jobroot_);
        lattice_->writeLatticePNG(time_[i],sim_type_,jobroot_);
        throw ex;
    }
    catch (MicrostructureException mex) {
        lattice_->writeLattice(time_[i],sim_type_,jobroot_);
        lattice_->writeLatticePNG(time_[i],sim_type_,jobroot_);
        throw mex;
    }

    ///
    /// Check if there is any capillary pore water remaining.  If not then
    /// we ASSUME hydration has stopped.
    ///
    /// @todo Generalize this idea to allow nanopore water to react by taking
    /// into account its lower chemical potential.
   
    if ((time_[i] >= output_time_[time_index]) && (time_index < output_time_.size())) {
        if (verbose_) {
            cout << "Writing lattice now... time_[" << i << "] = "
                 << time_[i] << ", output_time_[" << time_index << "] = "
                 << output_time_[time_index] << endl;
        }
        lattice_->writeLattice(time_[i],sim_type_,jobroot_);
        lattice_->writeLatticePNG(time_[i],sim_type_,jobroot_);
        // lattice_->CheckPoint(jobroot_);
        time_index++;
        if (verbose_) cout << "...Done!" << endl;
    }
    
    if (!capwater) {  // We will stop hydration
        if (warning_) {
            cout << "WARNING: System is out of capillary pore water." << endl;
            cout << "         This version of code assumes that only capillary" << endl;
            cout << "         water is chemically reactive, so the system is" << endl;
            cout << "         is assumed to be incapable of further hydration." << endl;
            cout.flush();
        }
    } 

    ///
    /// The following block executes only for sulfate attack simulations
    ///

    if (time_[i] >= sattack_time_) {
   
      map<int, vector<double> > expansion;
      expansion = lattice_->getExpansion();
      
      ifstream instopexp("stopexp.dat");
      if (!instopexp) {
        if (verbose_) cout << "keep expanding." << endl;
      } else {
        expansion.clear();
        cout << "expansion has been stopped due to the percolation of damage." << endl;
      }

      ///
      /// Stop FM temporarily
      /// @todo What is this?  The following if block will never be run if uncommented!
      ///
 
      /*
      expansion.clear();
      */

      if (verbose_) cout << "In Controller, expansion.size() = " << expansion.size() << endl;
      if (expansion.size() > 1) {
          if (verbose_) {
              cout << "time_ is: " << time_[i] << ". expansion.size() is: " << expansion.size()
                   << " now create new microstructure..." << endl; 
          }
    
        damagecount_ = 0;  
        double poreintroduce = 0.5;
    
        lattice_->writeLattice(time_[i],sim_type_,jobroot_);
        lattice_->writeLatticePNG(time_[i],sim_type_,jobroot_);
        string ofname(jobroot_);
        ostringstream ostr1,ostr2;
        ostr1 << (int)(time_[i] * 100);
        ostr2 << setprecision(3) << chemsys_->getTemperature();
        string timestr(ostr1.str());
        string tempstr(ostr2.str());
        ofname = ofname + "." + timestr + "." + tempstr + ".img";
    
        ///
        /// In the sulfate attack algorithm, calculate the stress and strain distributions
        ///
 
        thermalstr_->setEigen();
        for (map<int, vector<double> >::iterator it = expansion.begin(); 
             it != expansion.end(); it++) {

          int expindex = it->first;
          vector<double> expanval = it->second;
          vector<int> expcoordin = lattice_->getExpansionCoordin(expindex);
          thermalstr_->setEigen(expindex,expanval[0],expanval[1],expanval[2],0.0,0.0,0.0);
          thermalstr_->setExp(expindex,expcoordin);

          ///
          /// Set expansion site to be damaged if there is one, as determined by
          /// the setEigen function returning every site above damage stress threshold
          ///
          Site *ste;
          ste = lattice_->getSite(expindex);
          /*
          lattice_->dWaterchange(poreintroduce);
          */
      
          double dwmcval = poreintroduce;
          lattice_->dWmc(expindex,dwmcval);
          for (int j = 0; j < ste->nbSize(2); j++) {
            Site *stenb = ste->nb(j);
            stenb->dWmc(dwmcval);
          }
        }
    
        ///
        /// Calculate the stress-free strain (thermal strain) state in the microstructure,
        /// and then write the displacement field
        ///
 
        if (verbose_) cout << "Entering ThermalStrain calculation..." << endl;
        thermalstr_->Calc(time_[i],ofname,0.0,0.0,0.0,0.0,0.0,0.0);
        //thermalstr_ -> writeStress(jobroot_,time_[i],0); //write strxx
        //thermalstr_ -> writeStrainEngy(jobroot_,time_[i]);
        thermalstr_->writeDisp(jobroot_,time_[i]);

        ///
        /// Get the true volume of each voxel after FEM calculation
        ///

        double otruevolume = 0.0;
        double truevolume = 0.0;
        for (int ii = 0; ii < lattice_->getNumsites(); ii++) {
          Site *ste;
          ste = lattice_->getSite(i);
          otruevolume = ste->getTruevolume();
          for (int j = 0; j < 3; j++) {
            truevolume += thermalstr_->getEleStrain(ii,j);
          }
          truevolume = otruevolume * (1 + truevolume);
          ste->setTruevolume(truevolume);
        }

        for (int index = 0; index < lattice_->getNumsites(); index++) {
          Site *ste;
          ste = lattice_->getSite(index);
          int pid = ste->getPhaseId();

          if (ste->IsDamage()) {
              damagecount_++;
          }

          /*
          if ((pid == DAMAGEID)) {
            double strxx, stryy, strzz;
            strxx = stryy = strzz = 0.0;
            strxx = thermalstr_->getEleStress(index,0);
            stryy = thermalstr_->getEleStress(index,1);
            strzz = thermalstr_->getEleStress(index,2);
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
              lattice_->setExpansion(index,damageexp);
              lattice_->dWaterchange(poreincrease);
              ste->setVolume(VOIDID,(ste->getVolume(VOIDID) + poreincrease));
            }
          }
          */

          ///
          /// The next block gets the stress in each voxel that does NOT
          /// contain a clinker phase (C3S, C2S, C3A, or C4AF), then determine
          /// if the voxel should be damaged as a result

          if (pid > C4AFID) {
            double strxx, stryy, strzz;
            strxx = stryy = strzz = 0.0;
            strxx = thermalstr_->getEleStress(index,0);
            stryy = thermalstr_->getEleStress(index,1);
            strzz = thermalstr_->getEleStress(index,2);
            if ((strxx >= thermalstr_->getTstrength(index)) 
              || (stryy >= thermalstr_->getTstrength(index)) 
              || (strzz >= thermalstr_->getTstrength(index))) {
              // if (verbose_) cout << "Phase " << pid << " is damaged." << endl;
              if (!ste->IsDamage()) {
                // if (verbose_) cout << " it has not been damaged before." << endl;
                ste->setDamage();
                damagecount_++;
                //lattice_->dWaterchange(poreintroduce);
        
                double dwmcval = poreintroduce;
                lattice_->dWmc(index,dwmcval);
                for (int j = 0; j < ste->nbSize(1); j++) {
                  Site *stenb = ste->nb(j);
                  stenb->dWmc(dwmcval);
                  if ((stenb->getWmc() > 0.0) && (stenb->getPhaseId() != WATERID) 
                    && (stenb->getPhaseId() != VOIDID)) {
                    lattice_->addDissolutionSite(stenb,stenb->getPhaseId());
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

        }   // End of loop over all voxels

        cout << "Time = " << time_[i] << " damagecount_ is: " << damagecount_ << endl;
        ofstream outdamage("damage.dat");
        outdamage << damagecount_;
        outdamage.close();
    
        string damagejobroot = jobroot_ + ".damage";
        lattice_->writeDamageLattice(time_[i],damagejobroot);
        lattice_->writeDamageLatticePNG(time_[i],damagejobroot);
        // to see whether new damage is generated
      }
    }
  }

  ///
  /// Write the lattice state to an ASCII file and to a PNG file for visualization
  ///

  // lattice_->writeLattice(jobroot_);
  // lattice_->writeLatticePNG(jobroot_);
    
  return;
}

void Controller::initializeState (double time)
{

    // Set up the amount of water in the system
    kineticmodel_->initializeMoles();      

    // Set the temperature
    double T = chemsys_->getTemperature();
    lattice_->setTemperature(T);

    ///
    /// Now that the method is done determining the change in moles of each IC,
    /// launch a thermodynamic calculation to determine new equilibrium state
    /// 
    /// The `ChemicalSystem` object provides an interface for these calculations
    /// 

    if (verbose_) {
        cout << "Going to launch thermodynamic calculation now... ";
        cout.flush();
    }
    try {
        chemsys_->calculateState(time,true);
        if (verbose_) {
            cout << "*Returned from ChemicalSystem::calculateState" << endl;
            cout << "*called by function Controller::initializeState" << endl;
            cout.flush();
        }
    }
    catch (GEMException gex) {
        throw gex;
    }

    if (verbose_) {
        cout << "Done!" << endl;
        cout.flush();
    }
    
  return;
}

void Controller::calculateState (double time,
                                 double dt,
                                 bool isfirst) 
{
  try {

    if (isfirst) {

        double T = chemsys_->getTemperature();
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
    impurityrelease.resize(chemsys_->getMicimpuritynum(),0.0);
       
    if (verbose_) {
        cout << "Before KineticModel::calculateKineticStep, print out "
             << "ICmoles for solution..." << endl;
    }
    solut_->getICmoles();
 
    ///
    /// Get the number of moles of each IC dissolved from kinetically controlled phases
    ///

    double T = lattice_->getTemperature();
    if (verbose_) cout << "Going into KineticModel::calculateKineticStep now... " << endl;
    kineticmodel_->calculateKineticStep(dt,T,isfirst);
    if (verbose_) cout << "Done!" << endl;

    ///
    /// The next block only operates for sulfate attack iterations
    /// It determines how many IC moles of hydrogen, oxygen, sulfur, etc to add
    ///

    if (time >= sattack_time_) {            

      if (verbose_) cout << "waterchange_ in Lattice is: "
                         << lattice_->getWaterchange() << endl;

      double addwatervol = lattice_->getWaterchange()
          / lattice_->getNumsites() * chemsys_->getMictotinitvolume();

      ///
      /// Get the molar volume of water from the GEM node
      ///

      double water_v0 = chemsys_->getNode()->DC_V0(chemsys_->getMic2DC(WATERID,0),
                                            chemsys_->getP(),chemsys_->getTemperature());
      double addwatermol = addwatervol / water_v0;

      if (verbose_) {
          cout << "Molar volume of water is: " << water_v0 << endl;
          cout << "The moles of water added into the system are: " << addwatermol << endl;
      }
      for (int i = 0; i < chemsys_->getICnum(); i++) {
        if (chemsys_->getICname(i) == "H") {
          if (verbose_) cout << "previous IC moles for H is: " << chemsys_->getICmoles(i) << endl;
          chemsys_->setICmoles(i,(chemsys_->getICmoles(i) + 2 * addwatermol));
          if (verbose_) cout << "new ICmoles for H is: " << chemsys_->getICmoles(i) << endl;
        }
        if (chemsys_->getICname(i) == "O") {
          if (verbose_) cout << "previous IC moles for O is: " << chemsys_->getICmoles(i) << endl;
          chemsys_->setICmoles(i,(chemsys_->getICmoles(i) + addwatermol));
          if (verbose_) cout << "new ICmoles for O is: " << chemsys_->getICmoles(i) << endl;
        }
      }
    }
        
    ///
    /// Now that the method is done determining the change in moles of each IC,
    /// launch a thermodynamic calculation to determine new equilibrium state
    /// 
    /// The `ChemicalSystem` object provides an interface for these calculations
    /// 

    if (verbose_) cout << "Going to launch thermodynamic calculation now with isfirst = "
                       << isfirst << endl;;
    try {
        chemsys_->calculateState(time,isfirst);
        if (verbose_) {
            cout << "*Returned from ChemicalSystem::calculateState" << endl;
            cout << "*called by function Controller::calculateState" << endl;
            cout.flush();
        }
    }
    catch (GEMException gex) {
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

    try {
        double aveSI = 0.0;
        double moles = 0.0;
        vector<int> micphasemembers;
        for (int i = 0; i < chemsys_->getMicphasenum(); ++i) {
            int newmicid = chemsys_->getMicid(chemsys_->getMicphasename(i));
            if (verbose_) {
                cout << "Mic phase " << chemsys_->getMicphasename(i)
                     << " is associated with GEMS phases ";
                cout.flush();
            }
            aveSI = moles = 0.0;
            micphasemembers = chemsys_->getMicphasemembers(newmicid);
            for (int ii = 0; ii < micphasemembers.size(); ++ii) {
                int newpid = micphasemembers[ii];
                if (verbose_) {
                    if (ii == 0) {
                        cout << chemsys_->getPhasename(newpid);
                    } else {
                        cout << " , " << chemsys_->getPhasename(newpid);
                    }
                }
                aveSI += (solut_->getSI(newpid) * chemsys_->getPhasemoles(newpid));
                moles += chemsys_->getPhasemoles(newpid);
            }
            if (moles > 0.0) {
                aveSI = aveSI / moles;
            } else {
                aveSI = aveSI / (double)(micphasemembers.size());
            }
            if (verbose_) {
                cout << endl;
                cout.flush();
            }
            lattice_->setSI(newmicid,aveSI);    
        }
    }
    catch (EOBException eex) {
        eex.printException();
        exit(1);
    }
    
    if (isfirst) {        
      string csfilename("ChemSysOutput.dat");
      ofstream out2(csfilename.c_str(),ios::app);
      chemsys_->writeChemSys(out2);
      out2.close();
    }
        
    ///
    /// Set the kinetic DC moles.  This adds the clinker components to the DC moles.
    ///
    /// @todo Find out what this is and why it needs to be done
    ///

    if (verbose_) {
        cout << "About to enter setKineticDCmoles..." << endl;
        cout.flush();
    }
    kineticmodel_->setKineticDCmoles();
    if (verbose_) {
        cout << "Finished with setKineticDCmoles..." << endl;
        cout.flush();
    }
        
    // Output to files the solution composition data, phase data, DC data,
    // microstructure data, pH, and C-S-H composition and Ca/Si ratio

    string outfilename = jobroot_ + "_Solution.csv";
    ofstream out3(outfilename.c_str(),ios::app);
    if (!out3) {
      throw FileException("Controller","calculateState",
                         outfilename,"Could not append");
    }

    out3 << setprecision(5) << time;
    char cc;
    for (int i = 0; i < chemsys_->getDCnum(); i++) {
      cc = chemsys_->getDCclasscode(i); 
      if (cc == 'S' || cc == 'T' || cc == 'W') {
        out3 << "," << (chemsys_->getNode())->Get_cDC((long int) i);
      }
    }   
    out3 << endl;
    out3.close();
        
    outfilename = jobroot_ + "_DCVolumes.csv";
    ofstream out4(outfilename.c_str(),ios::app);
    if (!out4) {
      throw FileException("Controller","calculateState",
                        outfilename,"Could not append");
    }

    out4 << setprecision(5) << time;
    if (verbose_) cout << "Writing DC volumes file at time = " << time << endl;
    for (int i = 0; i < chemsys_->getDCnum(); i++) {
      if (chemsys_->getDCmolarmass(i) > 0.0) {
        cc = chemsys_->getDCclasscode(i); 
        if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M' || cc == 'W') {
            string dcname = chemsys_->getDCname(i);
            double V0 = chemsys_->getDCmoles(dcname) * chemsys_->getDCmolarvolume(dcname);
            out4 << "," << V0;
            if (verbose_) {
                cout << "    DC = " << chemsys_->getDCname(i)
                    << ", moles = " << chemsys_->getDCmoles(i) << ", molar mass = "
                    << chemsys_->getDCmolarmass(i) << endl;
            }
        }
      } else {
        string msg = "Divide by zero error for DC " + chemsys_->getDCname(i);
        out4.close();
        throw FloatException("Controller","calculateState",msg);
      }
    }   
    out4 << endl;
    out4.close();
        
    outfilename = jobroot_ + "_Microstructure.csv";
    ofstream out5(outfilename.c_str(),ios::app);
    if (!out5) {
      throw FileException("Controller","calculateState",outfilename,
                        "Could not append");
    }
    if (verbose_) cout << "Writing microstructure volume fractions at time " << time << endl;
    out5 << setprecision(5) << time;
    for (int i = 0; i < chemsys_->getMicphasenum(); i++) {
      out5 << "," << (lattice_->getVolumefraction(i));
    }   
    out5 << endl;
    out5.close();
        
    outfilename = jobroot_ + "_pH.csv";
    ofstream out6(outfilename.c_str(),ios::app);
    if (!out6) {
      throw FileException("Controller","calculateState",outfilename,
                        "Could not append");
    }
    if (verbose_) {
        cout << "Writing pH values...";
        cout.flush();
    }
    out6 << setprecision(5) << time;
    out6 << "," << (chemsys_->getPH()) << endl;
    if (verbose_) {
        cout << "Done!" << endl;
        cout.flush();
    }
    out6.close();
        
    if (verbose_) {
        cout << "Call to chemsys_setPhasestoich() ... " << endl;
        cout.flush();
    }
    chemsys_->setPhasestoich();
    if (verbose_) {
        cout << "Done!" << endl;
        cout << "Call to chemsys_getPhasestoich() ... " << endl;
        cout.flush();
    }
    double *CSHcomp;
    try {
        CSHcomp = chemsys_->getPhasestoich(chemsys_->getPhaseid("CSHQ"));
    }
    catch (EOBException eex) {
        eex.printException();
        exit(1);
    }
    if (verbose_) {
        cout << "Done!" << endl;
        cout.flush();
    }
    double Ca_moles = 0.0, Si_moles = 0.0, Ca_Si_Ratio = 0.0;
    outfilename = jobroot_ + "_CSH.csv";
    ofstream out7(outfilename.c_str(),ios::app);
    if(!out7){
      throw FileException("Controller","calculateState",
                        outfilename,"Could not append");
    }
    out7 << setprecision(5) << time;
    for (unsigned int i = 0; i < chemsys_->getICnum(); i++){
      out7 << "," << CSHcomp[i];
      if (chemsys_->getICname(i) == "Ca"){
        Ca_moles = CSHcomp[i];
      }
      if (chemsys_->getICname(i) == "Si"){
        Si_moles = CSHcomp[i];
      }
    }
    if (Ca_moles < 1.0e-16) Ca_moles = 1.0e-16;
    if (Si_moles < 1.0e-16) Si_moles = 1.0e-16;
    Ca_Si_Ratio = Ca_moles/Si_moles;
    out7 << "," << Ca_Si_Ratio << endl;
    out7.close();
        
    chemsys_->setPhasestoich();
    double *phaserecord;
    int ICindex;
    Ca_moles = Si_moles = 0.0;
    double C_S_ratio = 0.0;
    outfilename = jobroot_ + "_CSratio_solid.csv";
    ofstream out8(outfilename.c_str(),ios::app);
    if(!out8){
      throw FileException("Controller","calculateState",
                        outfilename,"Could not append");
    }
    out8 << setprecision(5) << time;
    for (int i = 0; i < chemsys_->getPhasenum(); i++) {
      cc = chemsys_->getPhaseclasscode(i);
      if (cc == 's'){
        phaserecord = chemsys_->getPhasestoich(i);
        ICindex = chemsys_->getICid("Ca");
        Ca_moles += phaserecord[ICindex];
        ICindex = chemsys_->getICid("Si");
        Si_moles += phaserecord[ICindex];
      }
    } 
    if (Si_moles != 0) {
      C_S_ratio = Ca_moles/Si_moles;
      out8 << "," << C_S_ratio << endl;
    }else{
      out8 << ",Si_moles is ZERO" << endl;
    }
    out8.close();
        
    outfilename = jobroot_ + "_icmoles.csv";
    ofstream out9(outfilename.c_str(),ios::app);
    if(!out9){
      throw FileException("Controller","calculateState",
                        outfilename,"Could not append");
    }
    out9 << setprecision(5) << time;
    for (int i = 0; i < chemsys_->getICnum(); i++) {
      out9 << "," << chemsys_->getICmoles(i);
    }
    out9 << endl;
    out9.close();

    /*
    ofstream out10("initicmoles.dat");
    for (int i = 0; i < chemsys_->getICnum(); i++) {
      out10 << chemsys_->getICmoles(i) << " " ;
    }
    out10 << endl;
    out10.close();*/
        
    ///
    /// Now that the end of the iteration is reached, zero out the kinetic
    /// DC moles in preparation for the next iteration
    ///
    /// @todo Find out what this is doing and why it is needed
    ///

    kineticmodel_->zeroKineticDCmoles();
  }
  catch (FileException fex) {
    fex.printException();
    exit(1);
  }
  catch (FloatException flex) {
    flex.printException();
    exit(1);
  }
  return;

}

void Controller::parseDoc (const string &docname) 
{
    xmlDocPtr doc;
    xmlNodePtr cur;
    xmlChar *key;
    cout.flush();
    double testtime;
    doc = xmlParseFile(docname.c_str());

    // check if the xml file is valid
    try {
        string paramxsd = xsd_files_path;
        paramxsd += "/parameters.xsd";
        if(!is_xml_valid(doc,paramxsd.c_str())) {
            cout << "XML file is NOT valid" << endl;
            throw FileException("Controller","parseDoc",
                                docname,"XML not valid");
        } else {
            int boolval;
	    if (doc == NULL ) {
                throw FileException("Controller","parseDoc",
                                    docname,"XML not parsed successfully");
            }
  	
            cur = xmlDocGetRootElement(doc);
	
            if (cur == NULL) {
                throw FileException("Controller","parseDoc",
                                    docname,"Empty parameter XML document");
            }


            cur = cur->xmlChildrenNode;
            while (cur != NULL) {
                if ((!xmlStrcmp(cur->name, (const xmlChar *)"calctime"))) {
                  key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                  string st((char *)key);
                  from_string(testtime,st);
                  time_.push_back(testtime);
                  xmlFree(key);
                }
                if ((!xmlStrcmp(cur->name, (const xmlChar *)"outtime"))) {
                  key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                  string st((char *)key);
                  from_string(testtime,st);
                  output_time_.push_back(testtime);
                  xmlFree(key);
                }
                if ((!xmlStrcmp(cur->name, (const xmlChar *)"image_frequency"))) {
                  key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
                  string st((char *)key);
                  from_string(imgfreq_,st);
                  xmlFree(key);
                }

                cur = cur->next;
            }

        }
        if (doc != NULL) xmlFreeDoc(doc);
    }
    catch (FileException fex) {
        fex.printException();
        exit(1);
    }
    return;
}
