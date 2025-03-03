/**
@file thames.cc

*/

#include "thames.h"
#include "version.h"
#include <limits>

/**
@brief The main block for running THAMES.

@return 0 on successful completion, non-zero otherwise
*/
int main(int argc, char **argv) {
  //
  // Set up the strainenergy vector.  We allow no more than 156 phases,
  // but this can be changed below.
  //

  checkargs(argc, argv);

  strainenergy.clear();
  strainenergy.resize(156, 0.0);

  int choice, simtype;
  string buff = "";
  ChemicalSystem *ChemSys = NULL;
  Lattice *Mic = NULL;
  ThermalStrain *ThermalStrainSolver = NULL;
  AppliedStrain *AppliedStrainSolver = NULL;
  KineticController *KController = NULL;
  Controller *Ctrl = NULL;
  RanGen *RNG = NULL;

  bool stopProgram = false;

  //
  // Main menu where user decides what kind of simulation this will be.
  //

  cout << "Enter choice: " << endl;
  cout << "  " << QUIT_PROGRAM << ") Exit program " << endl;
  cout << "  " << HYDRATION << ") Hydration " << endl;
  cout << "  " << LEACHING << ") Leaching " << endl;
  cout << "  " << SULFATE_ATTACK << ") Sulfate attack " << endl;
  cin >> choice;
  cout << endl << "choice = " << choice << endl;

  // cout << "epsilon for double : \t" << numeric_limits<double>::epsilon() <<
  // endl; cout << "epsilon for int : \t" << numeric_limits<int>::epsilon() <<
  // endl; cout << "epsilon for float : \t" << numeric_limits<float>::epsilon()
  // << endl;

  if (choice <= QUIT_PROGRAM || choice > SULFATE_ATTACK) {

    cout << "Exiting program now." << endl << endl;
    exit(1);
  }

  int seedRNG = -25943; // -142234;
  //cin >> seedRNG;
  cout << endl << "The RNG seed is                 : seedRNG = " << seedRNG << endl;

  double elemTimeInterval = 1.e-7;
  //cin >> elemTimeInterval;
  cout << "The elementary time interval is : elemTimeInterval = "
       << setprecision(3) << elemTimeInterval << endl;

  cout << scientific << setprecision(15) << endl;

  time_t lt = time(NULL);
  struct tm *inittime;
  inittime = localtime(&lt);
  cout << asctime(inittime);
  clock_t starttime = clock();
  simtype = choice;

  //
  // User must provide the name of the GEM chemical system definition (CSD) file
  // for the aqueous solution
  //

  // Read the newline character.  Wish there was a better way!
  getline(cin, buff);

  //
  // User must provide the name of the GEM CSD for the whole system
  //

  cout << endl << "What is the name of the GEM input file? " << endl;
  getline(cin, buff);
  const string geminput_filename(buff);
  cout << "geminput_filename : " << geminput_filename << endl;
  cout.flush();
  //
  // User must provide the name of the GEM data bridge (DBR) file
  // for the whole system
  //

  cout << endl << "What is the name of the GEM DBR file? " << endl;
  getline(cin, buff);
  const string geminput_dbrname(buff);
  cout << "geminput_dbrname  : " << geminput_dbrname << endl;

  //
  // User must provide the name of the file specifying the microstructre
  // phase data
  //

  cout << endl
       << "What is the name of the microstructure phase definition file? "
       << endl;
  getline(cin, buff);
  const string pi_filename(buff);
  const string cement_filename(buff);
  cout << "pi_filename       : " << pi_filename << endl;
  cout << "cement_filename   : " << cement_filename << endl;

  //
  // Create the ChemicalSystem object
  //

  try {
    ChemSys = new ChemicalSystem(geminput_filename, geminput_dbrname,
                                 pi_filename, VERBOSE, WARNING);
  } catch (bad_alloc &ba) {
    cout << "Bad memory allocation in ChemicalSystem constructor: " << ba.what()
         << endl;
    stopProgram = true;
  } catch (FileException fex) {
    fex.printException();
    stopProgram = true;
  } catch (GEMException gex) {
    gex.printException();
    stopProgram = true;
  } catch (DataException dex) {
    dex.printException();
    stopProgram = true;
  }
  if (stopProgram) {
    errorDelDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                        AppliedStrainSolver, KController, Ctrl, starttime, lt);
  }

  //
  // Create the random number generator
  //

  RNG = new RanGen(seedRNG);


  //
  // User must specifiy the file containing the 3D microstructure itself
  //

  cout << endl << "What is the name of the MICROSTRUCTURE file? " << endl;
  getline(cin, buff);
  const string mic_filename(buff);
  cout << "mic_filename      : " << mic_filename << endl;

  //
  // Create the Lattice object to hold the microstructure
  //

  try {
    Mic = new Lattice(ChemSys, RNG, seedRNG, mic_filename, VERBOSE, WARNING);
    cout << endl <<"Lattice creation done... " << endl;
    cout << "X size of lattice is " << Mic->getXDim() << endl;
    cout << "Y size of lattice is " << Mic->getYDim() << endl;
    cout << "Z size of lattice is " << Mic->getZDim() << endl;
    cout << "Total number of sites is " << Mic->getNumsites() << endl;
  } catch (bad_alloc &ba) {
    cout << "Bad memory allocation in Lattice constructor: " << ba.what()
         << endl;
    stopProgram = true;
  } catch (FileException ex) {
    ex.printException();
    stopProgram = true;
  } catch (GEMException ex) {
    ex.printException();
    stopProgram = true;
  }
  if (stopProgram) {
    errorDelDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                        AppliedStrainSolver, KController, Ctrl, starttime, lt);
  }

  if (choice == SULFATE_ATTACK) {

    //
    // This block is executed only if simulating external sulfate attack,
    // in which case we need information about the elastic moduli of the
    // constituent phases, and will need to include a finite element solver
    //

    cout << endl << "What is the name of the elastic modulus file?" << endl;
    buff = "";
    // cin >> buff;  // C++ >> operator does not allow spaces
    getline(cin, buff);
    const string phasemod_fileName(buff);
    cout << "phasemod_fileName : " << phasemod_fileName << endl;

    //
    // Create the ThermalStrain FE solver, which handles phase transformation
    // misfit
    //

    try {
      ThermalStrainSolver =
          new ThermalStrain(Mic->getXDim(), Mic->getYDim(), Mic->getZDim(),
                            (Mic->getNumsites() + 2),
                            ChemSys->getNumMicroPhases(), 1, VERBOSE, WARNING);
      cout << "ThermalStrain object creation done... " << endl;
      ThermalStrainSolver->setPhasemodfileName(phasemod_fileName);
    } catch (bad_alloc &ba) {
      cout << "Bad memory allocation in ThermalStrain constructor: "
           << ba.what() << endl;
      stopProgram = true;
    } catch (FileException ex) {
      ex.printException();
      stopProgram = true;
    } catch (GEMException ex) {
      ex.printException();
      stopProgram = true;
    }
    if (stopProgram) {
      errorDelDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                          AppliedStrainSolver, KController, Ctrl, starttime, lt);
    }

    int nx, ny, nz;
    nx = ny = nz = 3;
    int ns = nx * ny * nz;

    //
    // Create the AppliedStrain FE solver, which handles applied external
    // strains
    //

    try {
      AppliedStrainSolver = new AppliedStrain(
          nx, ny, nz, ns, ChemSys->getNumMicroPhases(), 1, VERBOSE, WARNING);
      AppliedStrainSolver->setPhasemodfileName(phasemod_fileName);
    } catch (bad_alloc &ba) {
      cout << "Bad memory allocation in AppliedStrain constructor: "
           << ba.what() << endl;
      stopProgram = true;
    } catch (FileException ex) {
      ex.printException();
      stopProgram = true;
    } catch (GEMException ex) {
      ex.printException();
      stopProgram = true;
    }
    if (stopProgram) {
      errorDelDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                          AppliedStrainSolver, KController, Ctrl, starttime, lt);
    }

    Mic->setFEsolver(AppliedStrainSolver);
  }

  string jobroot, par_filename, statfilename;
  if (VERBOSE) {
    cout << "About to enter KineticController constructor" << endl;
    cout.flush();
  }

  //
  // Create the KineticController object
  //

  try {
    KController =
        new KineticController(ChemSys, Mic, cement_filename, VERBOSE, WARNING);
  } catch (bad_alloc &ba) {
    cout << "Bad memory allocation in KineticController constructor: "
         << ba.what() << endl;
    stopProgram = true;
  } catch (FileException ex) {
    ex.printException();
    stopProgram = true;
  } catch (GEMException ex) {
    ex.printException();
    stopProgram = true;
  } catch (FloatException ex) {
    ex.printException();
    stopProgram = true;
  } catch (DataException ex) {
    ex.printException();
    stopProgram = true;
  }
  if (stopProgram) {
    errorDelDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                        AppliedStrainSolver, KController, Ctrl, starttime, lt);
  }

  if (VERBOSE) {
    cout << "Finished constructing KineticController KController" << endl;
    cout.flush();
  }
  cout << endl << "What is the name of the simulation parameter file? " << endl;
  getline(cin, buff);
  par_filename.assign(buff);
  cout << "par_filename      : " << par_filename << endl;

  cout << endl << "What is the root name of output files?" << endl;
  getline(cin, buff);
  jobroot.assign(buff);
  cout << "files root name   : " << jobroot << endl;

  buff = "mkdir -p Result";
  int resCallSystem = system(buff.c_str());
  if (resCallSystem == -1) {
    // handle the error
    cout << endl << endl << "    thames.cc - error : resCallSystem = -1" << endl;
    errorDelDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                        AppliedStrainSolver, KController, Ctrl, starttime, lt);
  }
  jobroot = "Result/" + jobroot;
  cout << "jobroot           : " << jobroot << endl;

  statfilename = jobroot + ".stats";

  //
  // Create the Controller object to direct flow of the program
  //

  try {
    Ctrl = new Controller(Mic, KController, ChemSys, ThermalStrainSolver,
                          simtype, par_filename, jobroot, VERBOSE, WARNING);
  } catch (bad_alloc &ba) {
    cout << "Bad memory allocation in Controller constructor: " << ba.what()
         << endl;
    stopProgram = true;
  } catch (FileException ex) {
    ex.printException();
    stopProgram = true;
  } catch (GEMException ex) {
    ex.printException();
    stopProgram = true;
  }
  if (stopProgram) {
    errorDelDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                        AppliedStrainSolver, KController, Ctrl, starttime, lt);
  }

  //
  // Write a formatted output of the simulation parameters for later reference
  //

  writeReport(jobroot, inittime, mic_filename, par_filename, geminput_filename,
              ChemSys, Ctrl);

  //
  // Launch the main controller to run the simulation
  //

  // if (VERBOSE) {
  cout << endl << "Going into Controller::doCycle" << endl;
  cout.flush();
  //}

  try {

    Ctrl->doCycle(statfilename, choice, elemTimeInterval);

  } catch (GEMException gex) {
    gex.printException();
    stopProgram = true;
  } catch (DataException dex) {
    dex.printException();
    stopProgram = true;
  } catch (EOBException ex) {
    ex.printException();
    stopProgram = true;
  } catch (MicrostructureException mex) {
    mex.printException();
    stopProgram = true;
  }
  if (stopProgram) {
    errorDelDynAllocMem(ChemSys, Mic, RNG, ThermalStrainSolver,
                        AppliedStrainSolver, KController, Ctrl, starttime, lt);
  }

  //
  // Simulation is finished.  Record and output the timing data.
  //

  cout << endl << "THAMES END" << endl;
  timeCount(starttime, lt);

  //
  // Delete the dynamically allocated memory
  //

  if (Ctrl) {
    delete Ctrl;
  }
  if (KController) {
    delete KController;
  }
  if (AppliedStrainSolver) {
    delete AppliedStrainSolver;
  }
  if (ThermalStrainSolver) {
    delete ThermalStrainSolver;
  }
  if (Mic) {
    delete Mic;
  }
  if (ChemSys) {
    delete ChemSys;
  }
  if (RNG) {
    delete RNG;
  }

  return 0;
}

void errorDelDynAllocMem (ChemicalSystem *ChemSys, Lattice *Mic, RanGen *RNG,
                          ThermalStrain *ThermalStrainSolver,
                          AppliedStrain *AppliedStrainSolver,
                          KineticController *KController,
                          Controller *Ctrl, clock_t st_time, time_t lt) {
  if (Ctrl) {
    delete Ctrl;
  }
  if (KController) {
    delete KController;
  }
  if (AppliedStrainSolver) {
    delete AppliedStrainSolver;
  }
  if (ThermalStrainSolver) {
    delete ThermalStrainSolver;
  }
  if (Mic) {
    delete Mic;
  }
  if (ChemSys) {
    delete ChemSys;
  }
  if (RNG) {
    delete RNG;
  }
  cout << endl << "STOP Program" << endl;
  timeCount(st_time, lt);
  exit(1);
}

void timeCount(clock_t time_, time_t lt_) {

  time_t lt1 = time(NULL);
  struct tm *inittime1;
  inittime1 = localtime(&lt1);
  cout << endl << asctime(inittime1);
  clock_t endtime = clock();

  double elapsedtime = (double)(endtime - time_) / CLOCKS_PER_SEC;
  double ltD = difftime(lt1, lt_);
  cout << endl << "Total time = " << ltD << " seconds" << endl;
  cout << endl
       << "Total time with clock = " << elapsedtime << " seconds" << endl;
}

void printHelp(void) {
  cout << endl;
  cout << "Usage: \"thames [--verbose|-v] [--help|-h]\"" << endl;
  cout << "        --verbose [-v]      Produce verbose output" << endl;
  cout << "        --suppress [-s]     Suppress warning messages" << endl;
  cout << "        --help [-h]         Print this help message" << endl;
  cout << endl;

  return;
}

void checkargs(int argc, char **argv) {

  // Many of the variables here are defined in the getopts.h system header file
  // Can define more options here if we want

  const char *const short_opts = "vsh";
  const option long_opts[] = {{"verbose", no_argument, nullptr, 'v'},
                              {"suppress", no_argument, nullptr, 's'},
                              {"help", no_argument, nullptr, 'h'},
                              {nullptr, no_argument, nullptr, 0}};

  VERBOSE = false;
  WARNING = true;

  while (true) {

    const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

    if (-1 == opt)
      break; // breaks out of this while loop

    switch (opt) {
    case 'v':
      VERBOSE = true; // Verbose defined in thameslib global.h
      cout << "**Will produce verbose output**" << endl;
      break;
    case 's':
      WARNING = false; // Verbose defined in thameslib global.h
      cout << "**Will suppress warning messages**" << endl;
      break;
    case 'h': // -h or --help
    case '?': // Unrecognized option
    default:
      printHelp();
      break;
    }
  }
}

void writeReport(const string &jobroot, struct tm *itime,
                 const string &mfileName, const string &parfilename,
                 const string &csname, ChemicalSystem *csys, Controller *ctr) {
  string statname = jobroot + ".stats";
  string jfilename = jobroot + ".report";
  ofstream out(jfilename.c_str());
  if (!out.is_open()) {
    if (WARNING)
      cout << "WARNING:  Could not open report file" << endl;
    return;
  }

  //
  // Write the time the job was executed
  //

  out << "THAMES job " << jobroot;
  out << " initialized on " << asctime(itime) << endl;
  out << endl;
  out << "INPUT FILES USED:" << endl;
  out << "   Microstructure file name: " << mfileName << endl;
  out << "        GEM input file name: " << csname << endl;
  out << endl;
  out << "OUTPUT FILES GENERATED:" << endl;
  out << "     Global phase fractions: " << statname << endl;
  out << endl;
  csys->writeChemSys();
  out << endl;

  out.close();
  return;
}
