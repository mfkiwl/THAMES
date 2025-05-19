/**
@file ChemicalSystem.cc
@brief Method definitions for the ChemicalSystem base class
*/

#include "ChemicalSystem.h"

string CSHMicroName("");
string MonocarbMicroName("");
string HydrotalcMicroName("");
string AFTMicroName("");
string MonosulfMicroName("");

ChemicalSystem::ChemicalSystem(const string &GEMfilename,
                               const string &jsonFileName, const bool verbose,
                               const bool warning) {
  int i, j;
  double *amat;
  string exmsg;
  long int gemflag = 0;

  double *icmolarmass, *dcmolarmass;
  char *cc;

  nodeStatus_ = NEED_GEM_AIA;
  nodeHandle_ = 0;
  iterDone_ = 0;
  timesGEMFailed_ = 0;
  maxGEMFails_ = 100000000; // 1000; // = 3;

  sulfateAttackTime_ = 1.0e10;
  leachTime_ = 1.0e10;
  beginAttackTime_ = 1.0e10;

  ///  The constructor initializes all the members to default values,
  ///  then launches the initial thermodynamic calculation, and sets
  ///  up the correspondences between GEM CSD phases and microstructure
  ///  phases.
  ///
  ///  All members are initialized to default values first, and all
  ///  vectors and maps are cleared.  The thermodynamic variables
  ///  are set to be consistent with neutral water at STP
  ///

  Eh_ = 0.0;
  T_ = 298.0;    // Default temperature [K]
  P_ = 101325.0; // Default pressure in [Pa]
  Vs_ = 1.0;
  Gs_ = 0.0;
  Ms_ = 0.0;
  verbose_ = verbose;
  jsonFormat_ = false;
  warning_ = warning;
  numMicroPhases_ = numGEMPhases_ = numSolutionPhases_ = 0;
  initScaledCementMass_ = scaledCementMass_ = 0;
  numMicroImpurities_ = 4;
  microPhaseName_.clear();
  stressPhaseName_.clear();
  porousPhaseName_.clear();
  microPhaseId_.clear();
  isKinetic_.clear();
  isParrotKilloh_.clear();
  isParrotKillohSize_ = 0;
  ICName_.clear();
  DCName_.clear();
  GEMPhaseName_.clear();
  numICs_ = numDCs_ = numGEMPhases_ = 0;
  microPhaseMembers_.clear();
  microPhaseMemberVolumeFraction_.clear();
  microPhaseDCMembers_.clear();
  growthTemplate_.clear();
  affinity_.clear();
  contactAngle_.clear();
  microPhasePorosity_.clear();
  poreSizeDistribution_.clear();
  k2o_.clear();
  na2o_.clear();
  mgo_.clear();
  so3_.clear();
  // RdICId_.clear();
  // Rd_.clear();
  GEMPhaseStoich_.clear();
  GEMPhaseDCMembers_.clear();
  microPhaseDCPorosities_.clear();
  microPhaseIdLookup_.clear();
  ICIdLookup_.clear();
  DCIdLookup_.clear();
  GEMPhaseIdLookup_.clear();
  microPhaseToGEMPhase_.clear();
  DCStoich_.clear();
  DCCharge_.clear();
  ICClassCode_.clear();
  DCClassCode_.clear();
  GEMPhaseClassCode_.clear();

  GEMPhaseName_.clear();
  ICIdLookup_.clear();
  DCIdLookup_.clear();
  GEMPhaseIdLookup_.clear();
  microPhaseVolume_.clear();
  microPhaseMass_.clear();
  microPhasePorosity_.clear();
  microPhaseMassDissolved_.clear();
  initialSolutionComposition_.clear();
  fixedSolutionComposition_.clear();
  gasSolidRatio_ = 0.0;
  initialGasComposition_.clear();
  fixedGasComposition_.clear();
  cementComponent_.clear();

  color_.clear();
  colorN_.clear(); // used in initColorMap() and output files
  initColorMap();

  elasticModuli_.clear();
  initElasticModuliMap();

  SI_.clear();

  node_ = new TNode();

  ///
  /// Initialize the thermodynamic system for both hydrates and solution
  /// in order to initialize GEMPhaseVolume_
  ///

  char *cGEMfilename = (char *)GEMfilename.c_str();
  // if (verbose_) {
  cout << endl
       << "ChemicalSystem::Going into GEM_init (1) to read CSD file "
       << cGEMfilename << endl; // *-dat.lst
  // }

  /// Find out if the input data are in json format or in key-value format

  try {
    jsonFormat_ = isInputFormatJSON(cGEMfilename);

    ///
    /// GEM_init initializes the IPM and DCH data structures
    /// This function will read the IPM, DCH, and one or more DBRs
    /// The argument is type const char *, and is the name of the
    /// data.lst file with the names of the IPM,
    /// DCH, and root DBR file.
    ///
    /// Return values are :
    ///    0 if successful, and the node_ object will hold the data
    ///    1 if input file(s) were not found or are corrupt
    ///   -1 if internal memory allocation error occurred
    ///

    // if (jsonFormat_) {
    //     if (verbose) {
    //         cout << "Detected JSON input file format for "
    //              << "ChemicalSystem GEM data files" << endl;
    //         cout.flush();
    //     }
    //     getJSONFiles(cGEMfilename,json_dch,json_ipm,json_dbr);
    //     gemflag = node_->GEM_init(json_dch,json_ipm,json_dbr);
    // } else {
    //     if (verbose) {
    //         cout << "Detected key-value input file format for "
    //              << "ChemicalSystem GEM data files" << endl;
    //         cout.flush();
    //     }
    gemflag = node_->GEM_init(cGEMfilename);
    // }
  } catch (FileException fex) {
    throw fex;
  }

  if (gemflag == 1) {
    exmsg = "Bad return from GEM_init: " + GEMfilename + " missing or corrupt";
    throw GEMException("ChemicalSystem", "ChemicalSystem", exmsg);
  }
  if (gemflag == -1) {
    exmsg = "Bad return from GEM_init: internal memory allocation error";
    throw GEMException("ChemicalSystem", "ChemicalSystem", exmsg);
  }

  ///
  /// Determine the number of possible ICs, DCs, and phases from the
  /// GEM CSD input that was read by GEM-IPM during initialization
  ///

  numICs_ = (int)((node_->pCSD())->nIC);
  numDCs_ = (int)((node_->pCSD())->nDC);
  numGEMPhases_ = (int)((node_->pCSD())->nPH);
  numSolutionPhases_ = (int)((node_->pCSD())->nPS);

  ///
  /// Knowing the dimensions, allocate the memory for all the arrays that
  /// must be created to store thermodynamic calculation results and communicate
  /// them to the microstructure
  ///

  try {
    exmsg = "ICMoles_";
    ICMoles_ = new double[numICs_];
    exmsg = "ICResiduals_";
    ICResiduals_ = new double[numICs_];
    exmsg = "ICChemicalPotential_";
    ICChemicalPotential_ = new double[numICs_];
    exmsg = "DCMoles_";
    DCMoles_ = new double[numDCs_];
    exmsg = "DCH0_";
    DCH0_ = new double[numDCs_];
    exmsg = "DCActivityCoeff_";
    DCActivityCoeff_ = new double[numDCs_];
    exmsg = "GEMPhaseMoles_";
    GEMPhaseMoles_ = new double[numGEMPhases_];
    exmsg = "solutPhaseMoles_";
    solutPhaseMoles_ = new double[numGEMPhases_];
    exmsg = "prevGEMPhaseMoles_";
    prevGEMPhaseMoles_ = new double[numGEMPhases_];
    exmsg = "GEMPhaseVolume_";
    GEMPhaseVolume_ = new double[numGEMPhases_];
    exmsg = "solutPhaseVolume_";
    solutPhaseVolume_ = new double[numGEMPhases_];
    exmsg = "prevGEMPhaseVolume_";
    prevGEMPhaseVolume_ = new double[numGEMPhases_];
    exmsg = "GEMPhaseMass_";
    GEMPhaseMass_ = new double[numGEMPhases_];
    exmsg = "solutPhaseMass_";
    solutPhaseMass_ = new double[numGEMPhases_];
    exmsg = "prevGEMPhaseMass_";
    prevGEMPhaseMass_ = new double[numGEMPhases_];
    exmsg = "specificSurfaceArea_";
    specificSurfaceArea_ = new double[numGEMPhases_];
    exmsg = "carrier_";
    carrier_ = new double[numSolutionPhases_];
    exmsg = "DCUpperLimit_";
    DCUpperLimit_ = new double[numDCs_];
    exmsg = "DCLowerLimit_";
    DCLowerLimit_ = new double[numDCs_];
    exmsg = "pGEMPhaseStoich_";
    pGEMPhaseStoich_ = new double[numGEMPhases_ * numICs_];
    exmsg = "pSolidStoich_";
    pSolidStoich_ = new double[numICs_];
    exmsg = "pSolutPhaseStoich_";
    pSolutPhaseStoich_ = new double[numGEMPhases_ * numICs_];
    exmsg = "pSolutSolidStoich_";
    pSolutSolidStoich_ = new double[numICs_];
  } catch (bad_alloc &ba) {
    cout << endl << "Bad_alloc Exception Thrown:" << endl;
    cout << "    Details:" << endl;
    cout << "    Offending function ChemicalSystem::ChemicalSystem" << endl;
    cout << "    Error in allocating memory for array " << exmsg << endl;
    cerr << endl << "Bad_alloc Exception Thrown:" << endl;
    cerr << "    Details:" << endl;
    cerr << "    Offending function ChemicalSystem::ChemicalSystem" << endl;
    cerr << "    Error in allocating memory for array " << exmsg << endl;
    exit(1);
  }

  ///
  /// Attempt to run GEM with auto initial approximation (AIA)
  ///
  /// This starts the thermodynamic calculation and returns the results,
  /// including the ionic strength, pH, IC chemical potentials, DC moles, phase
  /// moles, phase volumes, and other results of the calculation.  All of these
  /// parameters are loaded into the THAMES vectors that keep track of these
  /// things, since, they were passed to the GEM calculation by reference.
  ///
  /// The argument is false if we wamt to use activity coefficients and
  /// speciation from a previous GEM_run, but is true if we want to use the
  /// activity coefficients and speciation stored in a DBR memory structure read
  /// from a DBR file
  ///
  /// Possible return values for nodeStatus_:
  ///    0 (NO_GEM_SOLVER): No GEM recalculation needed for node
  ///    1 (NEED_GEM_AIA) : Need GEM calc with LPP (auto initial approx, AIA)
  ///    2 (OK_GEM_AIA)   : OK after GEM calc with LPP AIA
  ///    3 (BAD_GEM_AIA)  : Not fully trusworthy result after calc with LPP AIA
  ///    4 (ERR_GEM_AIA)  : Failure (no result) in GEM calc with LPP AIA
  ///    5 (NEED_GEM_SIA) : Need GEM calc with no-LPP (smart initial approx,
  ///    SIA) 6 (OK_GEM_SIA)   : OK after GEM calc with SIA 7 (BAD_GEM_SIA)  :
  ///    Not fully trusworthy result after calc with SIA 8 (ERR_GEM_SIA)  :
  ///    Failure (no result) in GEM calc with SIA 9 (T_ERROR_GEM ) : Terminal
  ///    error (e.g., memory corruption). Need restart
  ///

  (node_->pCNode())->NodeStatusCH = NEED_GEM_AIA;
  if (verbose_) {
    cout << "ChemicalSystem::Constructor: Entering GEM_run (1) with node "
            "status = "
         << nodeStatus_ << endl;
    cout.flush();
  }

  /// All ICs and DCs, etc. are zeroed out right now

  nodeStatus_ = node_->GEM_run(true);
  if (verbose_) {
    cout << "Done! nodeStatus is " << nodeStatus_ << endl;
    cout.flush();
  }

  if (!(nodeStatus_ == OK_GEM_AIA || nodeStatus_ == OK_GEM_SIA)) {
    bool dothrow = false;
    cerr << "ERROR: Call to GEM_run in ChemicalSystem constructor had an "
            "issue..."
         << endl;
    cerr << "       nodeStatus_ = " << nodeStatus_;
    switch (nodeStatus_) {
    case NEED_GEM_AIA:
      exmsg = " Need GEM calc with auto initial approx (AIA)";
      cout << exmsg << endl;
      dothrow = true;
      break;
    case BAD_GEM_AIA:
      exmsg = " Untrustworthy result with auto initial approx (AIA)",
      cout << exmsg << endl;
      dothrow = true;
      break;
    case ERR_GEM_AIA:
      exmsg = " Failed result with auto initial approx (AIA)";
      cout << exmsg << endl;
      dothrow = true;
      break;
    case NEED_GEM_SIA:
      exmsg = " Need GEM calc with smart initial approx (SIA)";
      dothrow = true;
      break;
    case BAD_GEM_SIA:
      exmsg = " Untrustworthy result with smart initial approx (SIA)";
      cout << exmsg << endl;
      dothrow = true;
      break;
    case ERR_GEM_SIA:
      exmsg = " Failed result with smart initial approx (SIA)";
      cout << exmsg << endl;
      dothrow = true;
      break;
    case T_ERROR_GEM:
      exmsg = " Terminal GEM error; need restart";
      cout << exmsg << endl;
      dothrow = true;
      break;
    case NO_GEM_SOLVER:
      exmsg = " No GEM recalculation needed for node";
      cout << exmsg << endl;
      dothrow = false;
      break;
    }
    if (dothrow) {
      throw GEMException("ChemicalSystem", "Constructor", exmsg);
    }
  }

  if (verbose_) {
    cout << "ChemicalSystem::Constructor: Entering GEM_restore_MT (1) ... "
         << endl;
    cout.flush();
  }

  ///
  /// Next call passes (copies) the GEMS3K input data from the DBR structure.
  /// This is useful after the GEM_init and GEM_run() calls to initialize the
  /// arrays that keep the chemical data for all the nodes (one node in most
  /// cases)
  ///
  /// This function returns nothing and appears unable of throwing exceptions
  /// @todo Check carefully whether this function can throw an exception
  ///

  node_->GEM_restore_MT(nodeHandle_, nodeStatus_, T_, P_, Vs_, Ms_,
                        &ICMoles_[0], &DCUpperLimit_[0], &DCLowerLimit_[0],
                        &specificSurfaceArea_[0]);

  if (verbose_) {
    cout << "Done!" << endl;
    cout << "ChemicalSystem::Constructor: Entering GEM_to_MT (1) ... " << endl;
    cout.flush();
  }

  ///
  /// Next call retrieves the GEMIPM chemical speciation calculation
  /// results from the DBR structure instance into memory provided by
  /// the THAMES code.  The dimensions and ordering of the arrays must
  /// correspond to those in currently existing DCH memory structure
  ///

  node_->GEM_to_MT(nodeHandle_, nodeStatus_, iterDone_, Vs_, Ms_, Gs_, Hs_,
                   ionicStrength_, pH_, pe_, Eh_, &ICResiduals_[0],
                   &ICChemicalPotential_[0], &DCMoles_[0], &DCActivityCoeff_[0],
                   &GEMPhaseMoles_[0], &GEMPhaseVolume_[0], &GEMPhaseMass_[0],
                   &pGEMPhaseStoich_[0], &carrier_[0], &specificSurfaceArea_[0],
                   &pSolidStoich_[0]);

  /// At this point all the IC moles and DC moles have the values that
  /// are loaded in the <bIC> vector of the DBR file.

  /// The results of the thermodynamic calculation are now known, and
  /// the constructor can cast them into appropriate units and set up
  /// the data structure to make correspondences between GEM and microstructure
  ///
  /// Convert all IC and DC molar masses from kg/mol to g/mol
  ///

  ICMolarMass_.resize(numICs_, 0.0);
  icmolarmass = (node_->pCSD())->ICmm;
  for (i = 0; i < numICs_; i++) {
    // Convert to g per mole
    ICMolarMass_[i] = (1000.0 * static_cast<double>(*icmolarmass));
    icmolarmass++;
  }
  DCMolarMass_.resize(numDCs_, 0.0);
  dcmolarmass = (node_->pCSD())->DCmm;
  for (i = 0; i < numDCs_; i++) {
    // Convert to g per mole
    DCMolarMass_[i] = (1000.0 * static_cast<double>(*dcmolarmass));
    dcmolarmass++;
  }

  /// Get the molar enthalpy of each DC component

  for (i = 0; i < numDCs_; i++) {
    DCH0_[i] = node_->DC_H0(i, P_, T_);
  }

  int maxICnameLength = node_->getMaxICnameLength();
  int maxDCnameLength = node_->getMaxDCnameLength();
  int maxPHnameLength = node_->getMaxPHnameLength();

  string string1;
  for (i = 0; i < numICs_; i++) {
    string1 = node_->xCH_to_IC_name(i);
    if ((int)(string1.length()) >= maxICnameLength)
      string1.resize(maxICnameLength);
    if (verbose_) {
      cout << "IC number " << i << " is " << string1 << endl;
    }
    ICName_.push_back(string1);
    ICIdLookup_.insert(make_pair(string1, i));
  }
  for (i = 0; i < numDCs_; i++) {
    string1 = node_->xCH_to_DC_name(i);
    if ((int)(string1.length()) >= maxDCnameLength)
      string1.resize(maxDCnameLength);
    if (verbose_) {
      cout << "DC number " << i << " is " << string1 << endl;
    }
    DCName_.push_back(string1);
    DCIdLookup_.insert(make_pair(string1, i));
  }
  for (int i = 0; i < numGEMPhases_; i++) {
    string1 = node_->xCH_to_Ph_name(i);
    if ((int)(string1.length()) >= maxPHnameLength)
      string1.resize(maxPHnameLength);
    if (verbose_) {
      cout << "PH number " << i << " is " << string1 << endl;
    }
    GEMPhaseName_.push_back(string1);
    GEMPhaseIdLookup_.insert(make_pair(string1, i));
  }

  ///
  /// Set up the stoichiometry matrix for dependent components (DCs) in terms
  /// of independent components (ICs).  This is the GEM CSD A matrix
  ///

  vector<double> scplaceholder;
  scplaceholder.clear();
  scplaceholder.resize(numICs_, 0);
  DCStoich_.resize(numDCs_, scplaceholder);
  DCCharge_.resize(numDCs_, 0.0);
  amat = (node_->pCSD())->A;
  for (i = 0; i < numDCs_; i++) {
    for (j = 0; j < numICs_; j++) {
      DCStoich_[i][j] = static_cast<double>(*amat);
      amat++;
    }
    DCCharge_[i] = static_cast<double>(DCStoich_[i][numICs_ - 1]);
  }

  ///
  /// Set up the stoichiometry and molar masses of the GEM CSD phases
  ///

  setPGEMPhaseStoich();
  setGEMPhaseStoich();

  ///
  /// Normally we can call setGEMPhaseMass() to read in the GEM phase masses
  /// from the GEM CSD, but this first time we have to do it manually
  /// because the units are converted from kg to g and this will mess
  /// up the assignment of the prevGEMPhaseMass_ values, which would still be
  /// in kg if we called the setGEMPhaseMass() function.
  ///

  for (int i = 0; i < numGEMPhases_; i++) {
    GEMPhaseMass_[i] =
        static_cast<double>(node_->Ph_Mass(i) * 1000.0); // in g, not kg
    prevGEMPhaseMass_[i] =
        static_cast<double>(node_->Ph_Mass(i) * 1000.0); // in g, not kg
  }

  setGEMPhaseVolume();
  setGEMPhaseMolarMass();

  ///
  /// Set up the class codes for ICs, DCs, and phases, based on the type of
  /// component they are.
  ///
  ///   IC class codes:
  ///       e = element
  ///       h = hydrogen
  ///       o = oxygen
  ///       a = Nit
  ///       z = Zz (charge)
  ///
  ///   DC class codes:
  ///       G = gas
  ///       I,M,J,O = solid components
  ///       W = water (solvent)
  ///       T = H^+
  ///       S = solute in electrolyte
  ///
  ///   Phase class codes:
  ///       a = aqueous solution
  ///       g = gas mixture
  ///       s = solid

  ICClassCode_.resize(numICs_, ' ');
  cc = (node_->pCSD())->ccIC;
  for (i = 0; i < numICs_; i++) {
    ICClassCode_[i] = *cc;
    cc++;
  }

  DCClassCode_.resize(numDCs_, ' ');
  cc = (node_->pCSD())->ccDC;
  for (i = 0; i < numDCs_; i++) {
    DCClassCode_[i] = *cc;
    cc++;
  }

  GEMPhaseClassCode_.resize(numGEMPhases_, ' ');
  cc = (node_->pCSD())->ccPH;
  for (int i = 0; i < numGEMPhases_; i++) {
    GEMPhaseClassCode_[i] = *cc;
    cc++;
  }

  if (verbose_) {
    cout << "ChemicalSystem::Constructor after GEM_run" << endl;
    cout << "    Number of ICs = " << getNumICs() << endl;
    for (int i = 0; i < getNumICs(); ++i) {
      cout << "        " << getICName(i) << ": " << getICMoles(i) << " moles"
           << endl;
    }
    // cout << endl;
    // cout << "    Number of DCs = " << chemSys_->getNumDCs() << endl;
    // for (int i = 0; i < chemSys_->getNumDCs(); ++i) {
    // cout << "        " << chemSys_->getDCName(i) << ": "
    // << chemSys_->getDCMoles(i) << " moles" << endl;
    // }
    cout << endl;
    cout.flush();
  }

  ///
  /// Begin parsing the chemistry input JSON file
  ///

  string msg;
  string jsonext = ".json";
  size_t foundjson = jsonFileName.find(jsonext);
  try {
    if (foundjson != string::npos) {
      parseDoc(jsonFileName);
      microPhaseVolume_.resize(numMicroPhases_, 0.0);
      microPhaseMass_.resize(numMicroPhases_, 0.0);
      microPhasePorosity_.resize(numMicroPhases_, 0.0);
      if (verbose_) {
        cout << " Setting microPhaseMass size to " << numMicroPhases_ << endl;
        cout.flush();
      }
      microPhaseMassDissolved_.resize(numMicroPhases_, 0.0);
    } else {
      msg = "Not a JSON file";
      throw FileException("ChemicalSystem", "ChemicalSystem", jsonFileName,
                          msg);
    }
  } catch (FileException fex) {
    throw fex;
  } catch (DataException dex) {
    throw dex;
  }

  ///
  /// Set up the main map that correlates microstructure phases with GEM CSD
  /// phases
  ///

  initMicroVolume_ = 0.0;
  // microPhaseBelongsToCement_.resize(numMicroPhases_,false);
  for (int i = 0; i < numMicroPhases_; i++) {
    microPhaseToGEMPhase_.insert(
        make_pair(static_cast<int>(i), microPhaseMembers_[i]));
  }

  ///
  /// Set up the vector of saturation indices for each GEM phase.
  /// This determines the driving force for growth and is also used in
  /// calculations of the crystallization pressure during external sulfate
  /// attack

  // setSI();

  isDCKinetic_.resize(numDCs_, false);
  DC_to_MPhID_.resize(numDCs_, -1);

  microPhaseSI_.clear();
  microPhaseSI_.resize(numMicroPhases_, 0.0);

  waterDCId_ = getDCId("H2O@");
  waterMollarMass_ = getDCMolarMass(waterDCId_);
  waterMollarVol_ = getDCMolarVolume(waterDCId_);

  aliteDCId_ = getDCId(AliteDCName);
  beliteDCId_ = getDCId(BeliteDCName);
  aluminateDCId_ = getDCId(AluminateDCName);
  ferriteDCId_ = getDCId(FerriteDCName);

  // for (int i = FIRST_SOLID; i < numMicroPhases_; i++) {
  //   cout << endl << "   affinity for i = " << i << endl;
  //   for (int j = FIRST_SOLID; j < numMicroPhases_; j++) {
  //     cout << "  " << j << "/" << affinity_[i][j] << "/" <<
  //     contactAngle_[i][j];
  //   }
  //   cout << endl;
  // }
  // cout << " exit chemsys" << endl; exit(0);

  // checkChemSys();

  writeSatElectrolyteGasConditions();
}

bool ChemicalSystem::isInputFormatJSON(const char *masterFileName) {
  ifstream in(masterFileName);
  if (!in) {
    throw FileException("ChemicalSystem", "isInputFormatJSON", masterFileName,
                        "Could not open");
  }

  string filetypeflag;
  if (in.peek() != EOF) {
    in >> filetypeflag;
  } else {
    throw FileException("ChemicalSystem", "isInputFormatJSON", masterFileName,
                        "Bad or corrupt format");
  }
  in.close();
  if (filetypeflag == "-j") {
    return true;
  }
  return false;
}

void ChemicalSystem::getJSONFiles(const char *masterFileName, string &dchName,
                                  string &ipmName, string &dbrName) {
  ifstream in(masterFileName);
  if (!in) {
    throw FileException("ChemicalSystem", "getJSONFiles", masterFileName,
                        "Could not open");
  }

  string filetypeflag;
  if (in.peek() != EOF) {
    in >> filetypeflag;
  } else {
    throw FileException("ChemicalSystem", "getJSONFiles", masterFileName,
                        "Bad or corrupt format");
  }
  if (in.peek() != EOF) {
    in >> dchName;
  } else {
    throw FileException("ChemicalSystem", "getJSONFiles", masterFileName,
                        "Bad or corrupt format");
  }
  if (in.peek() != EOF) {
    in >> ipmName;
  } else {
    throw FileException("ChemicalSystem", "getJSONFiles", masterFileName,
                        "Bad or corrupt format");
  }
  if (in.peek() != EOF) {
    in >> dbrName;
  } else {
    throw FileException("ChemicalSystem", "getJSONFiles", masterFileName,
                        "Bad or corrupt format");
  }

  /// Assuming that all the file names were read correctly, strip their quote
  /// marks

  dchName.erase(remove(dchName.begin(), dchName.end(), '"'), dchName.end());
  ipmName.erase(remove(ipmName.begin(), ipmName.end(), '"'), ipmName.end());
  dbrName.erase(remove(dbrName.begin(), dbrName.end(), '"'), dbrName.end());
}

vector<double> ChemicalSystem::getSolution(void) {
  ///
  /// Get IC moles for solution, which fully characterizes the
  /// composition of the aqueous solution
  ///

  double waterMass = DCMoles_[getDCId("H2O@")] * DCMolarMass_[getDCId("H2O@")];

  if (verbose_) {
    cout << "water mass is: " << waterMass << endl;
  }
  vector<double> tempICMoles;
  tempICMoles.clear();
  tempICMoles.resize(numICs_, 0.0);
  cout << endl;
  for (int i = 0; i < numDCs_; i++) {
    char cc = getDCClassCode(i);
    if (cc == 'S' || cc == 'T' || cc == 'W') {
      // cout << "tempICMoles " << i <<  "    DCName_: " << DCName_[i] << endl;
      for (int j = 0; j < (numICs_ - 1); j++) {
        tempICMoles[j] += DCMoles_[i] * DCStoich_[i][j];
      }
    }
  }
  cout << endl;

  /*
  for (int i = 0; i < numDCs_; i++) {
    char cc = getDCClassCode(i);
    if (cc == 'S' || cc == 'T') {
      double moles = node_->Get_cDC(i) * waterMass * 1.0e-3;
      for (int j = 0; j < (numICs_ - 1); j++) {
        tempICMoles[j] += moles * DCStoich_[i][j];
      }
    }
  }

  // Treat H2O separately

  double waterMoles = DCMoles_[getDCId("H2O@")];
  for (int j = 0; j < numICs_; j++) {
    if (ICName_[j] == "H")
      tempICMoles[j] += waterMoles * 2;
    if (ICName_[j] == "O")
      tempICMoles[j] += waterMoles;
  }
  */

  return tempICMoles;
}

void ChemicalSystem::parseDoc(const string &docName) {
  // Need to open the docName and scan it somehow for
  // phase names and id numbers

  PhaseData phaseData;

  /// Test for JSON existence and non-emptiness
  /// @todo Add a better JSON validity check.

  ifstream f(docName.c_str());
  if (!f.is_open()) {
    cout << "JSON parameter file not found" << endl;
    throw FileException("Controller", "parseDoc", docName, "File not found");
  }

  json data = json::parse(f);
  f.close();

  // Get an iterator to the root node of the JSON file
  json::iterator it = data.find("microstructure");
  json::iterator cdi = it.value().find("numentries");
  cdi = it.value().find("phases");

  map<string, int> phaseids;
  phaseids.clear();

  parseMicroPhaseNames(cdi, phaseids);

  ///
  /// Now start the iterator at the beginning and scan properly
  ///
  /// The file contains information about each microstructure
  /// phase that is defined, including the list of GEM CSD phases that
  /// are to be associated with that phase, the phase's internal porosity,
  /// dissolved impurities, and visualization properties.
  ///

  if (verbose_) {
    cout << "Back to the beginning to mine data from JSON file" << endl;
    cout.flush();
  }

  // Make "it" a json interator over the environment data object
  // to mine out the environmental (curing) input data
  it = data.find("environment");

  // Find saturation state
  cdi = it.value().find("saturated");
  int satstate = cdi.value();
  isSaturated_ = (satstate != 0) ? true : false;

  // See if electrolyte composition is specified
  cdi = it.value().find("electrolyte_conditions");
  if (cdi != it.value().end()) {
    try {
      parseSolutionComp(cdi);
    } catch (DataException dex) {
      throw dex;
    }
  }

  // See if gas composition is specified
  cdi = it.value().find("gas_conditions");

  if (cdi != it.value().end()) {
    try {
      parseGasComp(cdi);
    } catch (DataException dex) {
      throw dex;
    }
  } else {
  }

  // Now make "it" a json interator over the "microstructure" object
  // to mine out the microstructure phase input data
  it = data.find("microstructure");

  // Find number of phase entries
  cdi = it.value().find("numentries");
  int testnumEntries = cdi.value();

  cdi = it.value().find("phases");

  try {
    parseMicroPhases(cdi, testnumEntries, phaseids, phaseData);
  } catch (FileException fex) {
    // fex.printException();
    throw fex;
    cout << endl;
  } catch (GEMException gex) {
    gex.printException();
    cout << endl;
  }

  return;
}

void ChemicalSystem::parseSolutionComp(const json::iterator cdi) {

  // Clear the associative map to initialize it
  fixedSolutionComposition_.clear();
  initialSolutionComposition_.clear();
  attackSolutionComposition_.clear();

  string testName, testCondition;
  double testConc;
  int testDCId;

  json::iterator p = cdi.value()[0].find("DCName");
  // cout << "JSON iterator okay so far" << endl;
  // cout.flush();
  for (int i = 0; i < (int)(cdi.value().size()); ++i) {
    p = cdi.value()[i].find("DCname");
    // testName = p.value();
    testName = p.value().get<string>();
    testDCId = getDCId(testName);
    // Only add the data if this is a solution component
    if (DCClassCode_[testDCId] == 'S' || DCClassCode_[testDCId] == 'T') {
      p = cdi.value()[i].find("condition");
      // testCondition = p.value();
      testCondition = p.value().get<string>();
      p = cdi.value()[i].find("concentration");
      testConc = p.value();
      if (testCondition == "initial") {
        initialSolutionComposition_.insert(make_pair(testDCId, testConc));
      } else if (testCondition == "fixed") {
        fixedSolutionComposition_.insert(make_pair(testDCId, testConc));
      } else if (testCondition == "attack") {
        attackSolutionComposition_.insert(make_pair(testDCId, testConc));
      }
    }
  }

  // At this point all the composition constraints have been read
  // Make sure they have charge balance

  double totcharge = 0.0;
  map<int, double>::iterator it = initialSolutionComposition_.begin();
  while (it != initialSolutionComposition_.end()) {
    totcharge += ((it->second) * (DCCharge_[it->first]));
    // cout << DCName_[it->first]
    //      << ": Total initial charge so far = " << totcharge << endl;
    it++;
  }
  if (abs(totcharge) > 1.0e-9) {
    throw DataException("ChemicalSystem", "parseSolutionComp",
                        "Initial electrolyte charge imbalance");
  }

  totcharge = 0.0;
  it = fixedSolutionComposition_.begin();
  while (it != fixedSolutionComposition_.end()) {
    totcharge += ((it->second) * (DCCharge_[it->first]));
    // cout << DCName_[it->first] << ": Total fixed charge so far = " <<
    // totcharge
    //      << endl;
    it++;
  }
  if (abs(totcharge) > 1.0e-9) {
    throw DataException("ChemicalSystem", "parseSolutionComp",
                        "Fixed electrolyte charge imbalance");
  }

  totcharge = 0.0;
  it = attackSolutionComposition_.begin();
  while (it != attackSolutionComposition_.end()) {
    totcharge += ((it->second) * (DCCharge_[it->first]));
    // cout << DCName_[it->first]
    //      << ": Total initial charge so far = " << totcharge << endl;
    it++;
  }
  if (abs(totcharge) > 1.0e-9) {
    throw DataException("ChemicalSystem", "parseSolutionComp",
                        "Attack electrolyte charge imbalance");
  }

  return;
}

void ChemicalSystem::parseGasComp(const json::iterator cdi) {
  // Clear the associative map to initialize it

  fixedGasComposition_.clear();
  initialGasComposition_.clear();
  attackGasComposition_.clear();

  string DCName, testCondition;
  double DCMoles;
  int DCId;

  json::iterator p = cdi.value().begin();
  for (int i = 0; i < (int)(cdi.value().size()); ++i) {
    p = cdi.value()[i].find("name");
    DCName = p.value();
    DCId = getDCId(DCName);
    // Only add the data if this is a gas component
    if (DCClassCode_[DCId] == 'G' || DCClassCode_[DCId] == 'V' ||
        DCClassCode_[DCId] == 'H' || DCClassCode_[DCId] == 'N') {
      p = cdi.value()[i].find("condition");
      testCondition = p.value();
      p = cdi.value()[i].find("moles");
      DCMoles = p.value();
      if (testCondition == "initial") {
        initialGasComposition_.insert(make_pair(DCId, DCMoles));
      } else if (testCondition == "fixed") {
        fixedGasComposition_.insert(make_pair(DCId, DCMoles));
      } else if (testCondition == "attack") {
        attackGasComposition_.insert(make_pair(DCId, DCMoles));
      }
    }
  }

  // At this point all the composition constraints have been read
  // Make sure they have charge balance

  double totcharge = 0.0;
  map<int, double>::iterator it = initialGasComposition_.begin();
  while (it != initialGasComposition_.end()) {
    totcharge += ((it->second) * (DCCharge_[it->first]));
    // cout << DCName_[it->first]
    //      << ": Total initial charge so far = " << totcharge << endl;
    it++;
  }
  if (abs(totcharge) > 1.0e-9) {
    throw DataException("ChemicalSystem", "parseGasComp",
                        "Initial gas charge imbalance");
  }

  totcharge = 0.0;
  it = fixedGasComposition_.begin();
  while (it != fixedGasComposition_.end()) {
    totcharge += ((it->second) * (DCCharge_[it->first]));
    // cout << DCName_[it->first] << ": Total fixed charge so far = " <<
    // totcharge
    //      << endl;
    it++;
  }
  if (abs(totcharge) > 1.0e-9) {
    throw DataException("ChemicalSystem", "parseGasComp",
                        "Fixed gas charge imbalance");
  }

  totcharge = 0.0;
  it = attackGasComposition_.begin();
  while (it != attackGasComposition_.end()) {
    totcharge += ((it->second) * (DCCharge_[it->first]));
    // cout << DCName_[it->first] << ": Total fixed charge so far = " <<
    // totcharge
    //      << endl;
    it++;
  }
  if (abs(totcharge) > 1.0e-9) {
    throw DataException("ChemicalSystem", "parseGasComp",
                        "Attack gas charge imbalance");
  }

  return;
}

void ChemicalSystem::parseMicroPhaseNames(const json::iterator cdi,
                                          map<string, int> &phaseids) {

  json::iterator cii = cdi.value()[0].find("thamesname");
  string pname;
  int pid, cemcompval;
  bool cemComp;

  for (int pnum = 0; pnum < (int)(cdi.value().size()); ++pnum) {
    cii = cdi.value()[pnum].find("thamesname");
    pname = cii.value();
    cii = cdi.value()[pnum].find("id");
    pid = cii.value();
    if (verbose_) {
      cout << "    Phase name = " << pname << " with id " << pid << endl;
      cout.flush();
    }
    cii = cdi.value()[pnum].find("cement_component");
    cemcompval = cii.value();
    cemComp = (cemcompval == 1) ? true : false;
    phaseids.insert(make_pair(pname, pid));
    cementComponent_.push_back(cemComp);
  }
  return;
}

void ChemicalSystem::parseMicroPhases(const json::iterator cdi, int numEntries,
                                      map<string, int> phaseids,
                                      PhaseData &phaseData) {

  json::iterator p;
  for (int i = 0; i < (int)(cdi.value().size()); ++i) {

    phaseData.growthTemplate.clear();
    phaseData.affinity.clear();
    phaseData.contactAngle.clear();

    /// @note The affinity vector is always the same length, one entry for every
    /// microstructure phase, and the default value is zero (contactanglevalue =
    /// 180); for self-affinity (the growing phase and the template are the
    /// same) the default value is 1 (contactanglevalue = 0). Both, affinity and
    /// self-affinity can be modified supplying, in the chemistry.json file, the
    /// desired values for the contact angle

    // f(th) = (2 - 3 cos(th) + (cos(th))^3)/4
    // f(th) : [0, 1]
    // f(th) = 0 when th = 0 (perfect wetting)
    // f(th) = 1 when th = 180° (dewetting)
    // affinity(th) = 1 - f(th)
    // affinity(th) = 1 when th = 0 (perfect wetting)
    // affinity(th) = 0 when th = 180° (dewetting)
    //
    phaseData.contactAngle.resize(numEntries, 180);
    phaseData.affinity.resize(numEntries, 0.0);
    phaseData.k2o = 0.0;
    phaseData.na2o = 0.0;
    phaseData.mgo = 0.0;
    phaseData.so3 = 0.0;
    phaseData.red = 0;
    phaseData.green = 0;
    phaseData.blue = 0;
    phaseData.gray = 0;

    phaseData.DCId.clear();
    phaseData.DCName.clear();
    phaseData.GEMPhaseId.clear();
    phaseData.GEMPhaseName.clear();
    phaseData.microPhaseDCPorosities.clear();
    //  phaseData.RdId.clear();
    //  phaseData.RdVal.clear();
    phaseData.stressCalc = 0;
    phaseData.weak = 0;

    /// @note This parsing ignores the kinetic data portion for each
    /// phase.  The kinetic data parsing is handled by the KineticModel class

    p = cdi.value()[i].find("id");
    phaseData.id = p.value();
    p = cdi.value()[i].find("thamesname");
    if (p != cdi.value()[i].end()) {
      phaseData.thamesName = p.value();
    } else {
      p = cdi.value()[i].find("microphasename");
      if (p != cdi.value()[i].end()) {
        phaseData.thamesName = p.value();
      } else {
        throw DataException("ChemicalSystem", "parseMicroPhases",
                            "Microstructure phase name not found");
      }
    }
    p = cdi.value()[i].find("gemphase_data");
    if (p != cdi.value()[i].end()) {
      try {
        parseGEMPhaseData(p, phaseData);
      } catch (DataException dex) {
        throw dex;
      }
    } else if (phaseData.id != VOIDID) {
      throw DataException("ChemicalSystem", "parseMicroPhases",
                          "No GEM phase data for a microstructure phase");
    }
    p = cdi.value()[i].find("poresize_distribution");
    if (p != cdi.value()[i].end()) {
      parsePoreSizeDistribution(p, phaseData);
    }
    p = cdi.value()[i].find("stresscalc");
    if (p != cdi.value()[i].end()) {
      phaseData.stressCalc = p.value();
    }
    p = cdi.value()[i].find("weak");
    if (p != cdi.value()[i].end()) {
      // Weak means the phase can be damaged by stress
      phaseData.weak = p.value();
    }
    p = cdi.value()[i].find("display_data");
    if (p != cdi.value()[i].end()) {
      try {
        parseDisplayData(p, phaseData);
      } catch (DataException dex) {
        throw dex;
      }
    }
    p = cdi.value()[i].find("impurity_data");
    if (p != cdi.value()[i].end()) {
      try {
        parseImpurityData(p, phaseData);
      } catch (DataException dex) {
        throw dex;
      }
    }
    p = cdi.value()[i].find("interface_data");
    if (p != cdi.value()[i].end()) {
      try {
        parseInterfaceData(p, phaseids, phaseData);
      } catch (DataException dex) {
        throw dex;
      }
    }
    //    // Impurity partitioning data
    //    p = cdi.value()[pnum].find("Rd");
    //    if (p != cdi.value().end()) {
    //
    //      ///
    //      /// The data about partitioning of impurities among the hydration
    //      /// product phases as they grow are grouped within a complex
    //      /// field in the input
    //      /// file, so we have a special method to parse it.
    //      ///
    //      try {
    //        parseRdData(p, phaseData);
    //      }
    //    }

    // We have now found all the findable data for
    // this microstructure phase.  Time to process it

    if (phaseData.stressCalc > 0) {
      stressPhaseName_.push_back(phaseData.thamesName);
      stressPhaseId_.push_back(phaseData.id);
    }
    if (phaseData.weak > 0) {
      weakPhaseName_.push_back(phaseData.thamesName);
      weakPhaseId_.push_back(phaseData.id);
    }

    // Figure out if this is a porous phase or not
    // The criterion is that at least one of the constituent GEM phases
    // shall have at least one DC with porosity between 0.0 and 1.0

    bool done = false;
    for (int j = 0; j < (int)(phaseData.microPhaseDCPorosities.size()) && !done;
         ++j) {
      if (phaseData.microPhaseDCPorosities[j] > 0.0 &&
          phaseData.microPhaseDCPorosities[j] < 1.0) {
        porousPhaseName_.push_back(phaseData.thamesName);
        porousPhaseId_.push_back(phaseData.id);
        done = true;
      }
    }

    // Next is a lookup table that will return
    // the ordered list of DC porosities for a given microstructure phase id

    microPhaseDCPorosities_.insert(
        make_pair(phaseData.id, phaseData.microPhaseDCPorosities));

    phaseData.microPhaseDCPorosities.clear();

    microPhaseName_.push_back(phaseData.thamesName);
    microPhaseId_.push_back(phaseData.id);
    microPhaseIdLookup_.insert(make_pair(phaseData.thamesName, phaseData.id));

    double cs;
    for (int j = FIRST_SOLID; j < numEntries; j++) {
      if (j == phaseData.id && abs(phaseData.contactAngle[j] - 180.) < 1.e16) {
        phaseData.contactAngle[j] = 0;
      }
      cs = cos(phaseData.contactAngle[j] * Pi / 180.);
      phaseData.affinity[j] = 1.0 - (2. - 3 * cs + pow(cs, 3)) / 4.;
    }
    affinity_.push_back(phaseData.affinity);
    contactAngle_.push_back(phaseData.contactAngle);

    /// Growth template is based on positive affinities only

    growthTemplate_.push_back(calcGrowthTemplate(phaseData.affinity));

    poreSizeDistribution_.push_back(phaseData.poreSizeDist);
    if (verbose_) {
      cout << "Pushed pore size distribution data for phase "
           << phaseData.thamesName << endl;
      cout << "This phase distribution has " << phaseData.poreSizeDist.size()
           << " entries" << endl;
      cout << "Have now registered " << poreSizeDistribution_.size() << " PSDs"
           << endl;
      cout.flush();
    }
    phaseData.poreSizeDist.clear();
    grayscale_.push_back(phaseData.gray);
    color_.push_back(phaseData.colors);
    k2o_.push_back(phaseData.k2o);
    na2o_.push_back(phaseData.na2o);
    mgo_.push_back(phaseData.mgo);
    so3_.push_back(phaseData.so3);
    //  RdICId_.push_back(phaseData.RdId);
    //  Rd_.push_back(phaseData.RdVal);
    microPhaseMembers_.insert(make_pair(phaseData.id, phaseData.GEMPhaseId));
    microPhaseDCMembers_.insert(make_pair(phaseData.id, phaseData.DCId));

    numMicroPhases_++;

    if (verbose_) {
      cout << "Parsed phase " << phaseData.thamesName << endl;
      cout.flush();
    }
  }
  return;
}

void ChemicalSystem::parsePoreSizeDistribution(const json::iterator p,
                                               PhaseData &phaseData) {

  if (verbose_) {
    cout << "Reading Pore Size Distribution:" << endl;
    cout.flush();
  }

  double sum;
  struct PoreSizeVolume datarow;
  json::iterator pp;

  phaseData.poreSizeDist.clear();

  sum = 0.0;
  for (int i = 0; i < (int)(p.value().size()); ++i) {
    pp = p.value()[i].find("diameter");
    if (pp != p.value()[i].end()) {
      datarow.diam = pp.value();
      pp = p.value()[i].find("volumefraction");
      if (pp != p.value()[i].end()) {
        datarow.volfrac = pp.value();
      } else {
        cout << "      WARNING: Malformed data row for pore size distribution "
             << "(no volume fraction)... assuming zero" << endl;
        cout.flush();
        datarow.volfrac = 0.0;
      }
      phaseData.poreSizeDist.push_back(datarow);
      sum += datarow.volfrac;
      if (verbose_) {
        // i++;
        // cout << "---> " << "   i = " << i << "   " << datarow.diam << " , "
        // << datarow.volfrac << endl; cout.flush();
        cout << "---> " << datarow.diam << " , " << datarow.volfrac << endl;
        cout.flush();
      }
    } else {
      cout << "      WARNING: Malformed data row for pore size distribution "
           << "(no diameter)... skipping" << endl;
      cout.flush();
    }
  }

  if (verbose_) {
    cout << "<---- sum = " << sum << "   phaseData.poreSizeDist.size() : "
         << phaseData.poreSizeDist.size() << endl;
    cout.flush();
  }

  return;
}

void ChemicalSystem::parseGEMPhaseData(const json::iterator p,
                                       PhaseData &phaseData) {
  int GEMPhaseId = 0;
  string mypstr;
  phaseData.GEMPhaseDCMembers.clear();
  phaseData.microPhaseDCPorosities.clear();
  // bool scrapeWaterDCs = false;

  json::iterator pp;

  for (int i = 0; i < (int)(p.value().size()); ++i) {
    pp = p.value()[i].find("gemphasename");
    if (pp != p.value()[i].end()) {
      mypstr = pp.value().get<string>();
      phaseData.GEMPhaseName.push_back(const_cast<char *>(mypstr.c_str()));
      // mypstr = pp.value();
      // phaseData.GEMPhaseName.push_back(mypstr);
      // if (mypstr == WaterGEMName) {
      //   scrapeWaterDCs = true;
      // } else {
      //   scrapeWaterDCs = false;
      // }
      // cout << endl
      //      << "GEM Phase name = " << mypstr // << endl;
      //      << ", scrapeWaterDCs = " << scrapeWaterDCs << endl;
      // cout.flush();
      //  Assign the global microstructure phase name associated with CSH
      if (mypstr == CSHGEMName) {
        CSHMicroName = phaseData.thamesName;
      }
      if (mypstr == MonocarbGEMName) {
        MonocarbMicroName = phaseData.thamesName;
      }
      if (mypstr == HydrotalcGEMName) {
        HydrotalcMicroName = phaseData.thamesName;
      }
      if (mypstr == AFTGEMName) {
        AFTMicroName = phaseData.thamesName;
      }
      // if (mypstr == MonosulfGEMName) {
      //   AFTMicroName = phaseData.thamesName;
      // }
      GEMPhaseId = getGEMPhaseId(mypstr);
      phaseData.GEMPhaseId.push_back(GEMPhaseId);
    } else {
      throw DataException("ChemicalSystem", "parseGEMPhaseData",
                          "No GEM phase name given");
    }

    pp = p.value()[i].find("gemdc");
    if (pp != p.value()[i].end()) {
      try {
        parseGEMPhaseDCData(pp, phaseData);
      } catch (DataException dex) {
        throw dex;
      }
    } else {
      throw DataException("ChemicalSystem", "parseGEMPhaseData",
                          "No GEM DCs specified");
    }
    GEMPhaseDCMembers_.insert(
        make_pair(GEMPhaseId, phaseData.GEMPhaseDCMembers));
  }
  return;
}

void ChemicalSystem::parseGEMPhaseDCData(const json::iterator pp,
                                         PhaseData &phaseData) {
  string mydcstr;
  // char *mydchar = const_cast<char *>(str.c_str());

  int dcid = 0;
  double porosity = 0.0;
  bool scrapeWaterDCs = false;
  phaseData.GEMPhaseDCMembers.clear();

  json::iterator ppp = pp.value().begin();
  for (int i = 0; i < (int)(pp.value().size()); ++i) {
    ppp = pp.value()[i].find("gemdcname");
    if (ppp != pp.value()[i].end()) {
      mydcstr = ppp.value().get<string>();
      phaseData.DCName.push_back(const_cast<char *>(mydcstr.c_str()));
      // cout << "Phase " << phaseData.thamesName << " found DC " << mydcstr
      //     << endl;
      // if (mydcstr == AFTDCName) {
      //   AFTMicroName = phaseData.thamesName;
      // }
      if (mydcstr == MonosulfDCName) {
        MonosulfMicroName = phaseData.thamesName;
      }
      if (!scrapeWaterDCs || mydcstr == WaterDCName) {
        // Only for aqueous solution, we keep only the water DC,
        // not all the dissolved components
        dcid = getDCId(mydcstr);
        phaseData.DCId.push_back(dcid);
        phaseData.GEMPhaseDCMembers.push_back(dcid);
        // cout << "    GEM DC id = " << dcid // << endl;
        //     << ", scrapeWaterDCs = " << scrapeWaterDCs << endl;
      }
      // Make certain that there will be a porosity associated
      // with this DC
      if (phaseData.microPhaseDCPorosities.size() < phaseData.DCName.size()) {
        phaseData.microPhaseDCPorosities.push_back(0.0);
      }
    } else {
      throw DataException("ChemicalSystem", "parseGEMPhaseDCData",
                          "No DC name given");
    }

    ppp = pp.value()[i].find("gemdcporosity");
    if (ppp != pp.value()[i].end()) {
      porosity = ppp.value();
      // Make certain that there will be a porosity associated
      // with this DC
      if (phaseData.microPhaseDCPorosities.size() < phaseData.DCName.size()) {
        phaseData.microPhaseDCPorosities.push_back(porosity);
      } else {
        phaseData.microPhaseDCPorosities.at(
            phaseData.microPhaseDCPorosities.size() - 1) = porosity;
      }
    }
  }
  return;
}

void ChemicalSystem::parseDisplayData(const json::iterator p,
                                      PhaseData &phaseData) {

  int red, green, blue, gray;

  red = green = blue = gray = 0;
  bool rgbBool = false;

  json::iterator pp = p.value().find("red");
  if (pp != p.value().end()) {
    red = pp.value();
    rgbBool = true;
  }
  pp = p.value().find("green");
  if (pp != p.value().end()) {
    green = pp.value();
    rgbBool = true;
  }
  pp = p.value().find("blue");
  if (pp != p.value().end()) {
    blue = pp.value();
    rgbBool = true;
  }
  pp = p.value().find("gray");
  if (pp != p.value().end()) {
    gray = pp.value();
    rgbBool = true;
  }

  phaseData.colors.clear();
  phaseData.colors.push_back(red);
  phaseData.colors.push_back(green);
  phaseData.colors.push_back(blue);
  phaseData.gray = gray;

  float rgbfColor;

  if (rgbBool) {
    map<string, elemColor>::iterator it = colorN_.find(phaseData.thamesName);
    if (it != colorN_.end()) {
      colorN_[phaseData.thamesName].rgb[0] = red;
      colorN_[phaseData.thamesName].rgb[1] = green;
      colorN_[phaseData.thamesName].rgb[2] = blue;
      colorN_[phaseData.thamesName].gray = gray;

      colorN_[phaseData.thamesName].rgbf[0] =
          (float)(colorN_[phaseData.thamesName].rgb[0]) / 255.0;
      colorN_[phaseData.thamesName].rgbf[1] =
          (float)(colorN_[phaseData.thamesName].rgb[1]) / 255.0;
      colorN_[phaseData.thamesName].rgbf[2] =
          (float)(colorN_[phaseData.thamesName].rgb[2]) / 255.0;
      colorN_[phaseData.thamesName].grayf =
          (float)(colorN_[phaseData.thamesName].gray) / 255.0;

    } else {
      colorN_[phaseData.thamesName].colorId = colorN_.size();
      colorN_[phaseData.thamesName].altName = phaseData.thamesName;
      colorN_[phaseData.thamesName].rgb.push_back(red);
      colorN_[phaseData.thamesName].rgb.push_back(green);
      colorN_[phaseData.thamesName].rgb.push_back(blue);
      colorN_[phaseData.thamesName].gray = gray;

      rgbfColor = colorN_[phaseData.thamesName].rgb[0] / 255.0;
      colorN_[phaseData.thamesName].rgbf.push_back(rgbfColor);

      rgbfColor = colorN_[phaseData.thamesName].rgb[1] / 255.0;
      colorN_[phaseData.thamesName].rgbf.push_back(rgbfColor);

      rgbfColor = colorN_[phaseData.thamesName].rgb[2] / 255.0;
      colorN_[phaseData.thamesName].rgbf.push_back(rgbfColor);

      rgbfColor = colorN_[phaseData.thamesName].gray / 255.0;
      colorN_[phaseData.thamesName].grayf = rgbfColor;
    }
    // colorN_[phaseData.thamesName].rgbf[0] =
    //     (float)(colorN_[phaseData.thamesName].rgb[0]) / 255.0;
    // colorN_[phaseData.thamesName].rgbf[1] =
    //     (float)(colorN_[phaseData.thamesName].rgb[1]) / 255.0;
    // colorN_[phaseData.thamesName].rgbf[2] =
    //     (float)(colorN_[phaseData.thamesName].rgb[2]) / 255.0;
    // colorN_[phaseData.thamesName].grayf =
    //     (float)(colorN_[phaseData.thamesName].gray) / 255.0;
  }

  return;
}

void ChemicalSystem::parseImpurityData(const json::iterator p,
                                       PhaseData &phaseData) {

  json::iterator pp = p.value().begin();
  while (pp != p.value().end()) {
    if (pp.key() == "k2ocoeff") {
      phaseData.k2o = pp.value();
    }
    if (pp.key() == "na2ocoeff") {
      phaseData.na2o = pp.value();
    }
    if (pp.key() == "mgocoeff") {
      phaseData.mgo = pp.value();
    }
    if (pp.key() == "so3coeff") {
      phaseData.so3 = pp.value();
    }
    ++pp;
  }
  return;
}

void ChemicalSystem::parseInterfaceData(const json::iterator p,
                                        map<string, int> &phaseids,
                                        PhaseData &phaseData) {

  json::iterator pp = p.value().find("affinity");
  if (pp != p.value().end()) {
    try {
      parseAffinityData(pp, phaseids, phaseData);
    } catch (DataException dex) {
      throw dex;
    }
  }
  return;
}

void ChemicalSystem::parseAffinityData(const json::iterator pp,
                                       map<string, int> &phaseids,
                                       PhaseData &phaseData) {
  int testaftyid;
  double testangleval = 90.0;
  // map<string, int>::iterator it = phaseids.begin();
  string mystr("Null");

  json::iterator ppp;
  for (int i = 0; i < (int)(pp.value().size()); ++i) {
    ppp = pp.value()[i].find("affinityphase");
    if (ppp != pp.value()[i].end()) {
      mystr = ppp.value();
      map<string, int>::iterator it = phaseids.find(mystr);
      if (it != phaseids.end()) {
        testaftyid = it->second;
        ppp = pp.value()[i].find("contactanglevalue");
        if (ppp != pp.value()[i].end()) {
          testangleval = ppp.value();
          phaseData.contactAngle[testaftyid] = testangleval;
        }
      }
    } else {
      cout << "      WARNING: No valid substrate given... skipping" << endl;
      cout.flush();
    }
  }

  return;
}

// void ChemicalSystem::parseRdData(const json::iterator p,
//                                  struct PhaseData &phaseData) {
//   int RdId;
//   double RdVal;
//
//   json::iterator pp = p.value().begin();
//   while (pp != p.value().end()) {
//     if (pp.key() == "Rdelement") {
//       string st(pp.value());
//       RdId = getICId(st);
//       phaseData.RdId.push_back(RdId);
//     }
//
//     if (pp.key() == "Rdvalue") {
//       RdVal = pp.value();
//       phaseData.RdVal.push_back(RdVal);
//     }
//     ++pp;
//   }
// }

ChemicalSystem::ChemicalSystem(const ChemicalSystem &obj) {

  ///
  /// This is a straightforward copy constructor
  /// Each member and vector is copied into the newly constructed object
  ///

  numMicroPhases_ = obj.getNumMicroPhases();
  numMicroImpurities_ = obj.getNumMicroImpurities();
  numICs_ = obj.getNumICs();
  numDCs_ = obj.getNumDCs();
  numGEMPhases_ = obj.getNumGEMPhases();
  numSolutionPhases_ = obj.getNumSolutionPhases();
  microPhaseName_ = obj.getMicroPhaseName();
  ICName_ = obj.getICName();
  DCName_ = obj.getDCName();
  GEMPhaseName_ = obj.getGEMPhaseName();
  GEMPhaseDCMembers_ = obj.getGEMPhaseDCMembers();
  microPhaseId_ = obj.getMicroPhaseId();
  isKinetic_ = obj.getIsKinetic();
  ICMoles_ = obj.getICMoles();
  DCMoles_ = obj.getDCMoles();
  ICMolarMass_ = obj.getICMolarMass();
  DCMolarMass_ = obj.getDCMolarMass();
  growthTemplate_ = obj.getGrowthTemplate();
  affinity_ = obj.getAffinity();
  contactAngle_ = obj.getContactAngle();
  microPhaseMembers_ = obj.getMicroPhaseMembers();
  microPhaseMemberVolumeFraction_ = obj.getMicroPhaseMemberVolumeFraction();
  microPhaseDCMembers_ = obj.getMicroPhaseDCMembers();
  microPhasePorosity_ = obj.getMicroPhasePorosity();
  poreSizeDistribution_ = obj.getPoreSizeDistribution();
  k2o_ = obj.getK2o();
  na2o_ = obj.getNa2o();
  mgo_ = obj.getMgo();
  so3_ = obj.getSo3();
  // RdICId_ = obj.getRdICId();
  // Rd_ = obj.getRd();
  grayscale_ = obj.getGrayscale();
  color_ = obj.getColor();
  microPhaseIdLookup_ = obj.getMicroPhaseIdLookup();
  ICIdLookup_ = obj.getICIdLookup();
  DCIdLookup_ = obj.getDCIdLookup();
  GEMPhaseIdLookup_ = obj.getGEMPhaseIdLookup();
  microPhaseToGEMPhase_ = obj.getMicroPhaseToGEMPhase();
  ICClassCode_ = obj.getICClassCode();
  DCClassCode_ = obj.getDCClassCode();
  GEMPhaseClassCode_ = obj.getGEMPhaseClassCode();
  DCStoich_ = obj.getDCStoich();
  DCCharge_ = obj.getDCCharge();
  pGEMPhaseStoich_ = obj.getPGEMPhaseStoich();
  GEMPhaseStoich_ = obj.getGEMPhaseStoich();
  ICResiduals_ = obj.getICResiduals();
  ICChemicalPotential_ = obj.getICChemicalPotential();
  DCActivityCoeff_ = obj.getDCActivityCoeff();
  GEMPhaseMoles_ = obj.getGEMPhaseMoles();
  prevGEMPhaseMoles_ = obj.getPrevGEMPhaseMoles();
  GEMPhaseMass_ = obj.getGEMPhaseMass();
  prevGEMPhaseMass_ = obj.getPrevGEMPhaseMass();
  GEMPhaseVolume_ = obj.getGEMPhaseVolume();
  prevGEMPhaseVolume_ = obj.getPrevGEMPhaseVolume();
  carrier_ = obj.getCarrier();
  specificSurfaceArea_ = obj.getSpecificSurfaceArea();
  DCLowerLimit_ = obj.getDCLowerLimit();
  DCUpperLimit_ = obj.getDCUpperLimit();
  /*
  node_ = obj.getNode();
  */
  T_ = obj.getTemperature();
  P_ = obj.getP();
  Vs_ = obj.getVs();
  Ms_ = obj.getMs();
  pH_ = obj.getPH();
  pe_ = obj.getPe();
  Eh_ = obj.getEh();
  ionicStrength_ = obj.getIonicStrength();
  Gs_ = obj.getGs();
  Hs_ = obj.getHs();
  nodeHandle_ = obj.getNodeHandle();
  nodeStatus_ = obj.getNodeStatus();
  iterDone_ = obj.getIterDone();
  microPhaseVolume_ = obj.getMicroPhaseVolume();
  microVolume_ = obj.getMicroVolume();
  initMicroVolume_ = obj.getInitMicroVolume();
  microPhaseMass_ = obj.getMicroPhaseMass();
  microPhaseMassDissolved_ = obj.getMicroPhaseMassDissolved();
  microVoidVolume_ = obj.getMicroVoidVolume();
  verbose_ = obj.getVerbose();
  warning_ = obj.getWarning();
}

ChemicalSystem::~ChemicalSystem(void) {
  ///
  /// Clear out the maps
  ///

  microPhaseIdLookup_.clear();
  DCIdLookup_.clear();
  ICIdLookup_.clear();
  GEMPhaseIdLookup_.clear();
  microPhaseToGEMPhase_.clear();

  ///
  /// Clear out the vectors
  ///

  microPhaseName_.clear();
  ICName_.clear();
  DCName_.clear();
  GEMPhaseName_.clear();
  microPhaseId_.clear();
  DCStoich_.clear();
  DCCharge_.clear();
  growthTemplate_.clear();
  affinity_.clear();
  contactAngle_.clear();
  microPhaseMembers_.clear();
  microPhaseMemberVolumeFraction_.clear();
  microPhaseMass_.clear();
  microPhaseMassDissolved_.clear();
  microPhaseDCMembers_.clear();
  microPhasePorosity_.clear();
  poreSizeDistribution_.clear();
  k2o_.clear();
  na2o_.clear();
  mgo_.clear();
  so3_.clear();
  // RdICId_.clear();
  // Rd_.clear();
  grayscale_.clear();
  color_.clear();
  ICClassCode_.clear();
  DCClassCode_.clear();
  GEMPhaseClassCode_.clear();
  GEMPhaseStoich_.clear();

  ///
  /// Free up the dynamically allocated memory
  ///

  delete[] DCLowerLimit_;
  delete[] DCUpperLimit_;
  delete[] specificSurfaceArea_;
  delete[] prevGEMPhaseMass_;
  delete[] prevGEMPhaseVolume_;
  delete[] prevGEMPhaseMoles_;
  delete[] GEMPhaseMass_;
  delete[] GEMPhaseVolume_;
  delete[] carrier_;
  delete[] GEMPhaseMoles_;
  delete[] DCActivityCoeff_;
  delete[] DCMoles_;
  delete[] DCH0_;
  delete[] ICChemicalPotential_;
  delete[] ICResiduals_;
  delete[] ICMoles_;
  delete[] pGEMPhaseStoich_;
  delete[] pSolidStoich_;
  delete[] pSolutSolidStoich_;
  delete[] solutPhaseMoles_;
  delete[] solutPhaseVolume_;
  delete[] solutPhaseMass_;
  delete[] pSolutPhaseStoich_;

  delete node_;
  node_ = NULL;
}

void ChemicalSystem::getPGEMPhaseStoich(void) {
  double *arout = new double[numICs_];
  for (int i = 0; i < numGEMPhases_; i++) {
    arout = node_->Ph_BC(i, arout);
    for (int j = 0; j < numICs_; j++) {
      pGEMPhaseStoich_[(i * numICs_) + j] = arout[j];
    }
  }
  delete[] arout;
}

void ChemicalSystem::getGEMPhaseStoich(void) {
  double minval = 0.0;
  GEMPhaseStoich_.clear();
  vector<double> vplace;
  vplace.clear();
  vplace.resize(numICs_, 0.0);
  GEMPhaseStoich_.resize(numGEMPhases_, vplace);
  int indexval, oval;
  for (int i = 0; i < numGEMPhases_; i++) {
    if (GEMPhaseName_[i] == "aq_gen") {

      ///
      /// Normalize to one mole of oxygen
      ///

      oval = (i * numICs_) + getICId("O");
      if (pGEMPhaseStoich_[oval] > 0.0) {
        for (int j = 0; j < numICs_; j++) {
          indexval = (i * numICs_) + j;
          GEMPhaseStoich_[i][j] =
              (pGEMPhaseStoich_[indexval] / pGEMPhaseStoich_[oval]);
        }
      }
    } else {
      for (int j = 0; j < numICs_; j++) {
        indexval = (i * numICs_) + j;
        GEMPhaseStoich_[i][j] = (pGEMPhaseStoich_[indexval] / minval);
      }
    }
  }
}

/*
void ChemicalSystem::writeDb(ostream &stream) {
  unsigned int i;

  ///
  /// Make the header
  ///

  stream << "--------------------------------------------------------" << endl;
  stream << "CONTENTS OF PHASE DATABASE:" << endl;
  stream << endl;

  ///
  /// Format one line at a time
  ///

  for (i = 0; i < numMicroPhases_; i++) {
    writeMember(i, stream);
  }

  stream << endl;
  stream << "--------------------------------------------------------" << endl;
  stream << endl;
}

void ChemicalSystem::writeMember(const unsigned int i, ostream &stream) {

  // unsigned int idnum;
  if (i >= numMicroPhases_) {
    throw EOBException("ChemicalSystem", "writeMember", "microPhaseName_",
                       numMicroPhases_, i);
  }

  ///
  /// Format the output for one phase
  ///

  stream << "------------------------------------------------------" << endl;
  stream << "DATA FOR MATERIAL " << i << ":" << endl;
  stream << "       Name = " << microPhaseName_[i] << endl;
  stream << "         Id = " << microPhaseId_[i] << endl;
  stream << "   Porosity = " << microPhasePorosity_[i] << endl;
  stream << "------------------------------------------------------" << endl;
}
*/

void ChemicalSystem::writeChemSys(ofstream &out) {
  int j;

  ///
  /// First we will list details for the ICs
  ///

  out << "Report on the Material Database" << endl;
  out << "-------------- ChemicalSystem -----------------" << endl << endl;
  out << endl << "List of Independent Components :" << endl;
  out << "numICs_ = " << numICs_ << endl;
  out << "IC_Id) ICName_[i] / ICClassCode_[i] / ICMolarMass_[i]" << endl;
  for (int i = 0; i < numICs_; i++) {
    out << i << ") Name: " << ICName_[i] << endl;
    out << "        classcode: " << ICClassCode_[i] << endl;
    out << "       molar mass: " << ICMolarMass_[i] << endl << endl;
  }

  out << endl << "List of Dependent Components:" << endl;
  out << "numDCs_ = " << numDCs_ << endl;
  out << "DC_Id) DCName_[i] / DCClassCode_[i] / DCMolarMass_[i] / "
         "node_->DCtoPh_DBR(i)"
      << endl;
  for (int i = 0; i < numDCs_; i++) {
    out << endl << i << ") Name: " << DCName_[i] << endl;
    out << "        classcode: " << DCClassCode_[i] << endl;
    out << "       molar mass: " << DCMolarMass_[i] << endl;
    out << "        DBR index: " << node_->DCtoPh_DBR(i)
        << "\t(GEM Phase : " << GEMPhaseName_[node_->DCtoPh_DBR(i)] << ")"
        << endl;
    if (DC_to_MPhID_[i] != -1) {
      out << "     microPhaseId: " << DC_to_MPhID_[i]
          << "\t(THAMES Phase : " << microPhaseName_[DC_to_MPhID_[i]] << ")"
          << endl;
    } else {
      out << "     microPhaseId: " << DC_to_MPhID_[i] << "\t- no THAMES Phase"
          << endl;
    }
  }

  out << endl << "List of GEM Phases:" << endl;
  out << "numGEMPhases_ = " << numGEMPhases_ << endl;
  out << "GEMPhase_Id) GEMPhaseName_[i] / GEMPhaseClassCode_[i]" << endl;
  for (int i = 0; i < numGEMPhases_; i++) {
    out << endl << i << ") Name: " << GEMPhaseName_[i] << endl;
    out << "        classcode: " << GEMPhaseClassCode_[i] << endl;
  }

  vector<int> compDC;
  out << endl << "List of Microstructure Phases (THAMES Phases):" << endl;
  out << "numMicroPhases_ = " << numMicroPhases_ << endl;
  out << "microPhase_Id) microPhaseName_[i] / microPhaseId_[i] / "
         "randomGrowth_[i] / affinity_[i][j] / growthTemplate_[i][j] / "
         "microPhasePorosity_[i]"
      << endl;
  for (int i = 0; i < numMicroPhases_; i++) {
    out << endl << i << ") Name: " << microPhaseName_[i] << endl;
    out << "                     id: " << microPhaseId_[i] << endl;
    out << "               affinity: " << endl;
    for (j = 0; j < (int)(affinity_[i].size()); j++) {
      if (affinity_[i][j] != 0) {
        out << "                  affinity to " << j << ": " << affinity_[i][j]
            << "   " << microPhaseName_[j] << endl;
      } else {
        out << "                  affinity to " << j << "t: " << affinity_[i][j]
            << endl;
      }
    }
    out << "        growthTemplate:";
    if (growthTemplate_[i].size() != 0) {
      for (j = 0; j < (int)(growthTemplate_[i].size()); j++) {
        out << " " << growthTemplate_[i][j] << "("
            << microPhaseName_[growthTemplate_[i][j]] << ")";
      }
    } else {
      out << "  - no templates";
    }
    out << endl;
    out << "              porosity: " << microPhasePorosity_[i] << endl;
    out << "            impurities: " << endl;
    out << "                  k2o_[i]: " << k2o_[i] << endl;
    out << "                 na2o_[i]: " << na2o_[i] << endl;
    out << "                  mgo_[i]: " << mgo_[i] << endl;
    out << "                  so3_[i]: " << so3_[i] << endl;
    compDC = getMicroPhaseDCMembers(i);
    out << "        DCs components:";
    if (compDC.size() != 0) {
      for (j = 0; j < (int)(compDC.size()); j++) {
        out << " " << compDC[j] << "(" << DCName_[compDC[j]] << ")";
      }
    } else {
      out << "  - no DC components";
    }
    out << endl;
  }

  out << endl << "        equivalence microPhaseId/DCId/GEMPhaseId:" << endl;
  for (int i = ELECTROLYTEID; i < numMicroPhases_; ++i) {
    // if (i >= 1) {
    string pname = getMicroPhaseName(i);
    int DCId = getMicroPhaseDCMembers(i, 0);
    int indDBR = node_->DCtoPh_DBR(DCId);

    out << endl << "   " << i << "\tpname: " << pname << endl;
    out << "          GEMPhaseId: " << indDBR
        << "\tGEMPhaseName_: " << GEMPhaseName_[indDBR] << endl;
    out << "          DCId      : " << DCId << "\tDCName_: " << DCName_[DCId]
        << endl;
    // }
  }

  return;
}

void ChemicalSystem::setMicroPhaseMass(const unsigned int idx,
                                       const double val) {
  try {
    microPhaseMass_.at(idx) = val;
  } catch (out_of_range &oor) {
    EOBException ex("ChemicalSystem", "setMicroPhaseMass", "microPhaseMass_",
                    microPhaseMass_.size(), idx);
    ex.printException();
    exit(1);
  }
  int DCId = 0;
  if (idx == ELECTROLYTEID) {
    DCId = getDCId("H2O@");
  } else if (idx != VOIDID) {
    DCId = getMicroPhaseDCMembers(idx, 0);
  }
  if (idx != VOIDID) {
    double v0 = node_->DC_V0(DCId, P_, T_);
    double dcmm = getDCMolarMass(DCId);
    if (dcmm < 1.0e-9) {
      FloatException fex("ChemicalSystem", "setMicroPhaseMass",
                         "Divide by zero (dcmm)");
      fex.printException();
      exit(1);
    }
    setMicroPhaseVolume(idx, (val * v0 / dcmm));
  }

  return;
}

void ChemicalSystem::calcMicroPhasePorosity(const unsigned int idx) {

  // No point in doing this if the phase is not present

  double porosity = 0.0;

  // Get the first GEM phase for this microstructure phase

  /// @todo Generalize calcMicroPhasePorosity so that we
  /// don't just check the first GEM Phase listed as part
  /// of this microstructure phase
  unsigned int gemphaseid = getMicroPhaseToGEMPhase(idx, 0);

  // Find all the DC ids for this GEM phase
  vector<int> DClist = getGEMPhaseDCMembers(gemphaseid);
  // Get the porosities for each DC in the microstructure phase
  int DClistSize = DClist.size();

  vector<double> DCporosities = getMicroPhaseDCPorosities(idx);

  if (verbose_) {
    cout << "ChemicalSystem::calcMicroPhasePorosity for "
         << getMicroPhaseName(idx) << endl;
    cout << "    This phase's GEM phase id = " << gemphaseid
         << " and volume = " << getMicroPhaseVolume(idx) << " m3" << endl;
    cout << "    The DC members for this GEM phase are:" << endl;
    for (int ii = 0; ii < DClistSize; ++ii) {
      cout << "        " << DClist[ii] << " (" << getDCName(DClist[ii]) << ")"
           << endl;
      cout.flush();
    }
  }

  /// @todo Do we need to check that DClist and DCporosities are the same
  /// size?
  /// @todo Do we need to check that DClist and DCporosities are in the same
  /// order?
  /// @note I think this is already guaranteed when we parse the phase data

  // Loop over them one by one and get each one's volume and
  // subvoxel porosity.  Create a volume weighted average
  // of porosities

  /// @todo Will this work for capillary porosity?  Conc will be molal

  double conc = 0.0; // Temporary variable for holding concentrations
                     // For solid phases this will be mole fraction

  int DCId = 0;
  double vol = 0.0;
  double sumvol = 0.0;
  double weightedporosities = 0.0;

  // JWB (2023-Apr-22) This algorithm works ONLY if there are no
  // solid solutions that are kinetic phases.  But THAMES is not
  // set up for that in all sorts of ways and this is not the
  // biggest of them.

  if (DClistSize == 1) {
    DCId = DClist[0];
    porosity = DCporosities[0];
    conc = 1.0;
    vol = conc * getDCMolarVolume(DCId);
    if (verbose_) {
      cout << "    " << getDCName(DCId) << " (" << DCId
           << ") concentration = " << conc << endl;
      cout << "    " << getDCName(DCId) << " (" << DCId
           << ") porosity = " << porosity << endl;
      cout << "    " << getDCName(DCId) << " (" << DCId
           << ") molarvolume = " << getDCMolarVolume(DCId) << endl;
      cout << "****" << endl;
      cout.flush();
    }
  } else {
    // int DClistSize = DClist.size();
    for (int i = 0; i < DClistSize; ++i) {
      DCId = DClist[i];
      conc = getDCConcentration(DCId);
      porosity = DCporosities[i];
      vol = conc * getDCMolarVolume(DCId);
      weightedporosities += (vol * porosity);
      sumvol += vol;
      if (verbose_) {
        cout << "    " << getDCName(DCId) << " (" << DCId
             << ") concentration = " << conc << endl;
        cout << "    " << getDCName(DCId) << " (" << DCId
             << ") porosity = " << porosity << endl;
        cout << "    " << getDCName(DCId) << " (" << DCId
             << ") molarvolume = " << getDCMolarVolume(DCId) << endl;
        cout << "****" << endl;
        cout.flush();
      }
    }
    if (sumvol > 0.0) {
      porosity = (weightedporosities / sumvol);
    } else {
      porosity = 0.0;
    }
    if (verbose_) {
      cout << "    " << getMicroPhaseName(idx)
           << " subvoxel porosity = " << porosity << endl;
      cout.flush();
    }
  }

  // setMicroPhasePorosity(idx, porosity);
  microPhasePorosity_[idx] = porosity;

  return;
}

int ChemicalSystem::calculateState(double time, bool isFirst = false,
                                   int cyc = 0) {
  // int status = 0;
  string msg; // , iniStrNodeStatus, endStrNodeStatus;

  // isFirst = true;

  // vector<double> oDCMoles;
  // oDCMoles.clear();
  // oDCMoles.resize(numDCs_, 0.0);
  // for (int i = 0; i < numDCs_; i++) {
  //   oDCMoles[i] = DCMoles_[i];
  // }

  // writeICMoles();
  // writeDCMoles();

  vector<double> microPhaseVolumes = getMicroPhaseVolume();

  nodeStatus_ = NEED_GEM_AIA;
  // int iniNodeStatus = nodeStatus_;
  // if (iniNodeStatus == 1) {
  //   iniStrNodeStatus = "NEED_GEM_AIA";
  // } else if (iniNodeStatus == 5) {
  //   iniStrNodeStatus = "NEED_GEM_SIA";
  // }
  // ALL ICs/DCs in the system are set to zero in Lattice constructor before
  // to call normalizePhaseMasses() DCs are updated in
  // Lattice::normalizePhaseMasses only the ICMoles_ that are less than 10^-9
  // after the first call of calculateKineticStep(...) are set to 10^-9

  if (isFirst) {
    checkICMoles();
  }

  //  cout << endl << "chemSys after checkICMoles for cyc = " << cyc << " :
  //  ICMoles_/ICName_" << endl; for(int i = 0; i < numICs_; i++){
  //      cout << i << "\t" << ICMoles_[i] << "\t" << ICName_[i] << endl;
  //  }
  //  writeDCMoles();
  //  exit(0);

  ///
  /// Next function loads the input data for the THAMES node into the
  /// instance of the DBR structure.  This call precedes the GEM_run call
  ///
  /// This function returns nothing and appears to be incapable of throwing
  /// an exception.
  ///
  /// @todo Check carefully if this function can throw an exception
  ///
  /// @note MT in the function name stands for "mass transport", which is
  /// the generic designation given to the code that couples to GEMS, THAMES
  /// in this case.
  ///

  bool doAttack = (time >= beginAttackTime_) ? true : false;

  // cout << endl << "  ----> ChemicalSystem::calculateState -
  // time/beginAttackTime_/doAttack : "
  //      << time << " / " << beginAttackTime_ << " / " << doAttack << endl;

  // Check and set chemical conditions on electrolyte and gas phase
  setElectrolyteComposition(isFirst, doAttack, cyc);
  setGasComposition(isFirst, doAttack);

  if (doAttack) {
    DCLowerLimit_[aliteDCId_] = DCMoles_[aliteDCId_];
    DCLowerLimit_[beliteDCId_] = DCMoles_[beliteDCId_];
    DCLowerLimit_[aluminateDCId_] = DCMoles_[aluminateDCId_];
    DCLowerLimit_[ferriteDCId_] = DCMoles_[ferriteDCId_];
    // cout << endl
    //      << "  ----> ChemicalSystem::calculateState -  DCLowerLimit_:"
    //      << endl;
    // for (int i = 0; i < numDCs_; i++) {
    //   cout << "        i/DCLowerLimit_[i]/DCName_[i] -0- : " << i
    //        << " / " << DCLowerLimit_[i]
    //        << " / " << DCName_[i] << endl;
    // }

    /*
    cout << endl // check!
         << "  ----> ChemicalSystem::calculateState -  ICMoles_ : " << endl;
    for (int i = 0; i < numICs_; i++) {
      cout << "        i/ICMoles_/ICName_[i] -i- : " << i << " / " <<
    ICMoles_[i]
           << " / " << ICName_[i] << endl;
    }
    checkICMoles();
    cout;
    for (int i = 0; i < numICs_; i++) {
      cout << "        i/ICMoles_/ICName_[i] -f- : " << i << " / " <<
    ICMoles_[i]
           << " / " << ICName_[i] << endl;
    }
    */
  }

  if (verbose_) {
    cout << endl
         << "ChemicalSystem::calculateState Entering GEM_from_MT cyc = " << cyc
         << endl;
    cout << "DCMoles:" << endl;
    for (int i = 0; i < numDCs_; ++i) {
      cout << "    " << DCName_[i] << ": " << DCMoles_[i] << " mol, ["
           << DCLowerLimit_[i] << ", " << DCUpperLimit_[i] << "] mol" << endl;
    }
    cout.flush();
    // cout << "ICMoles:" << endl;
    // for (int i = 0; i < numICs_; ++i) {
    //   cout << "    " << ICName_[i] << ": " << ICMoles_[i] << endl;
    // }
    // cout.flush();
  }

  node_->GEM_from_MT(nodeHandle_, nodeStatus_, T_, P_, Vs_, Ms_, ICMoles_,
                     DCUpperLimit_, DCLowerLimit_, specificSurfaceArea_,
                     DCMoles_);

  if (isFirst) {
    for (int i = 0; i < numICs_; i++) {
      ICMoles_[i] = 0.0;
    }
  }

  if (verbose_) {
    cout << "ChemicalSystem::calculateState Exiting GEM_from_MT" << endl;
    cout << "DCMoles:" << endl;
    for (int i = 0; i < numDCs_; ++i) {
      cout << "    " << DCName_[i] << ": " << DCMoles_[i] << ", ["
           << DCLowerLimit_[i] << ", " << DCUpperLimit_[i] << "]" << endl;
    }
    cout.flush();
  }

  /// For passing the current THAMES time and time step into the working
  /// instance of the DBR structure.
  ///
  /// This function returns nothing and appears to be incapable of throwing an
  /// exception
  ///
  ///
  /// @todo Check carefully if this function can throw an exception
  ///
  /// @note This function call seems to be optional.  It is not needed for
  /// the GEM data structures or calculations, but mostly for debugging
  /// purposes.
  ///
  /// @note This function used to be called GEM_set_MT in older versions of
  /// GEMS3K
  ///

  node_->GEM_from_MT_time(time, 1.0);

  ///
  /// Attempt to run GEM with automatic initial approximation (AIA)
  ///
  /// This starts the thermodynamic calculation and returns the results,
  /// including the ionic strength, pH, IC chemical potentials, DC moles,
  /// phase moles, phase volumes, and other results of the calculation.  All
  /// of these parameters are loaded into the THAMES vectors that keep track
  /// of these things, since, they were passed to the GEM calculation by
  /// reference.
  ///
  /// The argument is false if we want to use activity coefficients and
  /// speciation from a previous GEM_run, but is true if we want to use the
  /// activity coefficients and speciation stored in a DBR memory structure
  /// read from a DBR file
  ///
  /// Possible return values for nodeStatus_:
  ///    0 (NO_GEM_SOLVER): No GEM recalculation needed for node
  ///    1 (NEED_GEM_AIA) : Need GEM calc with LPP (auto initial approx, AIA)
  ///    2 (OK_GEM_AIA)   : OK after GEM calc with LPP AIA
  ///    3 (BAD_GEM_AIA)  : Not fully trusworthy result after calc with LPP AIA
  ///    4 (ERR_GEM_AIA)  : Failure (no result) in GEM calc with LPP AIA
  ///    5 (NEED_GEM_SIA) : Need GEM calc with no-LPP (smart initial approx,
  ///    SIA) 6 (OK_GEM_SIA)   : OK after GEM calc with SIA 7 (BAD_GEM_SIA)  :
  ///    Not fully trusworthy result after calc with SIA 8 (ERR_GEM_SIA)  :
  ///    Failure (no result) in GEM calc with SIA 9 (T_ERROR_GEM ) : Terminal
  ///    error (e.g., memory corruption). Need restart
  ///

  nodeStatus_ = node_->GEM_run(true);

  if (verbose_) {
    cout << "Done!  nodeStatus is " << nodeStatus_ << endl;
    cout.flush();
  }

  if (!(nodeStatus_ == OK_GEM_AIA || nodeStatus_ == OK_GEM_SIA)) {
    bool dothrow = false;
    if (verbose_) {
      cerr << endl
           << "  ChemicalSystem::calculateState - GEM_run ERROR: nodeStatus_ = "
           << nodeStatus_ << endl;
    }
    switch (nodeStatus_) {
    case NEED_GEM_AIA:
      msg = "    Need GEM calc with auto initial approx (AIA)";
      cerr << msg << endl;
      dothrow = false;
      break;
    case BAD_GEM_AIA:
      msg = "    Untrustworthy result with auto initial approx (AIA)",
      cerr << msg << endl;
      dothrow = false;
      break;
    case ERR_GEM_AIA:
      msg =
          "  ChemicalSystem::calculateState - Failed result with auto initial "
          "approx (AIA)";
      if (verbose_) {
        cerr << msg << ", GEMS failed " << timesGEMFailed_ << " times" << endl;
      }
      node_->GEM_print_ipm("IPM_dump.txt");
      timesGEMFailed_++;
      dothrow = (timesGEMFailed_ > maxGEMFails_) ? true : false;
      break;
    case NEED_GEM_SIA:
      msg = "    Need GEM calc with smart initial approx (SIA)";
      cerr << msg << endl;
      dothrow = false;
      break;
    case BAD_GEM_SIA:
      msg = "    Untrustworthy result with smart initial approx (SIA)";
      cerr << msg << endl;
      dothrow = false;
      break;
    case ERR_GEM_SIA:
      msg = "    ChemicalSystem::calculateState - Failed result with smart "
            "initial "
            "approx (SIA)";
      cerr << msg << ", GEMS failed " << timesGEMFailed_ << " times" << endl;
      node_->GEM_print_ipm("IPM_dump.txt");
      timesGEMFailed_++;
      dothrow = (timesGEMFailed_ > maxGEMFails_) ? true : false;
      break;
    case T_ERROR_GEM:
      msg = "    Terminal GEM error; need restart";
      cerr << msg << endl;
      dothrow = true;
      break;
    case NO_GEM_SOLVER:
      msg = "    No GEM recalculation needed for node";
      cerr << msg << endl;
      dothrow = false;
      break;
    }
    if (dothrow) {
      throw GEMException("ChemicalSystem", "calculateState", msg);
    }
  } else {
    // string finStrNodeStatus;
    // if (nodeStatus_ == 2) {
    //   finStrNodeStatus = "OK_GEM_AIA";
    // } else if (iniNodeStatus == 6) {
    //   finStrNodeStatus = "OK_GEM_SIA";
    // }
    // cout << endl <<"  ChemicalSystem::calculateState - for cyc = " << cyc <<
    // " : OK" << endl; cout << "    ChemicalSystem::calculateState - initial
    // nodeStatus_ = " << iniNodeStatus << " [" << iniStrNodeStatus << "]" <<
    // endl; cout << "    ChemicalSystem::calculateState - final nodeStatus_   =
    // " << nodeStatus_ << " [" << finStrNodeStatus << "]" << endl;

    cout << "  ChemicalSystem::calculateState - cyc = " << cyc
         << " : OK & GEM_run has failed " << timesGEMFailed_
         << " consecutive times before to find this solution" << endl;

    // cout << "    ChemicalSystem::calculateState - solution for kinetic
    // controlled phases:" << endl; for (int i = 0; i < numMicroPhases_; i++) {
    //   if (isKinetic_[i]) {
    //     cout << "      i = " << setw(3) << right << i << " : " << setw(15) <<
    //     left
    //          << microPhaseName_[i]
    //          << " => updated scaledMass (microPhaseMass_[i]) = " <<
    //          microPhaseMass_[i]
    //          << " , microPhaseMassDissolved_[i] = " <<
    //          microPhaseMassDissolved_[i]
    //          << " and volume = " << microPhaseVolume_[i] << endl;
    //   }
    // }

    timesGEMFailed_ = 0;
  }

  if (timesGEMFailed_ > 0) {
    // cout << "  ChemicalSystem::calculateState - GEM_run has failed "
    //      << timesGEMFailed_ << " consecutive times  cyc = " << cyc << endl;
    // if (timesGEMFailed_ % 10000  == 0) {
    //   cout << endl << "  ChemicalSystem::calculateState - test : GEM_run has
    //   failed "
    //        << timesGEMFailed_ << " consecutive times  cyc = " << cyc << endl;
    // }
    return timesGEMFailed_;
  }

  if (verbose_) {
    cout << "    Going into ChemicalSystem::calculateState::GEM_to_MT (2)... ";
    cout.flush();
  }

  ///
  /// Next call retrieves the GEMIPM chemical speciation calculation
  /// results from the DBR structure instance into memory provided by
  /// the THAMES code.  The dimensions and ordering of the arrays must
  /// correspond to those in currently existing DCH memory structure
  ///
  /// This function returns nothing and appears unable of throwing exceptions
  /// @todo Check carefully whether this function can throw an exception
  ///

  node_->GEM_to_MT(nodeHandle_, nodeStatus_, iterDone_, Vs_, Ms_, Gs_, Hs_,
                   ionicStrength_, pH_, pe_, Eh_, &ICResiduals_[0],
                   &ICChemicalPotential_[0], &DCMoles_[0], &DCActivityCoeff_[0],
                   &solutPhaseMoles_[0], &solutPhaseVolume_[0],
                   &solutPhaseMass_[0], &pSolutPhaseStoich_[0], &carrier_[0],
                   &specificSurfaceArea_[0], &pSolidStoich_[0]);

  if (verbose_) {
    cout << endl << "Done!" << endl;
    cout << "ChemicalSystem::calculateState Exiting GEM_from_MT cyc = " << cyc
         << endl;
    cout << "DCMoles:" << endl;
    for (int i = 0; i < numDCs_; ++i) {
      cout << "    " << DCName_[i] << ": " << DCMoles_[i] << ", ["
           << DCLowerLimit_[i] << ", " << DCUpperLimit_[i] << "]" << endl;
    }
    cout << "ICMoles:" << endl;
    for (int i = 0; i < numICs_; ++i) {
      cout << "    " << ICName_[i] << ": " << ICMoles_[i] << endl;
    }
    cout << "after GEM_to_MT...Ms_ = " << Ms_ << ", Hs_ = " << Hs_ << endl;
    cout.flush();
  }

  // writePhasemoles();

  microVolume_ = 0.0;
  setPGEMPhaseStoich(); // call getPGEMPhaseStoich() => pGEMPhaseStoich_[i] //
                        // check! number of moles all ICs in all GEM CSD phases.
  setGEMPhaseStoich();  // call getGEMPhaseStoich() => GEMPhaseStoich_[i][j]  //
                        // check!
  setGEMPhaseMass();    // => GEMPhaseMass_[i]
  setGEMPhaseVolume();  // => GEMPhaseVolume_[i]
  setGEMPhaseMolarMass(); // =>GEMPhaseMolarMass_[pidx]  // check!

  if (verbose_) {
    cout << endl
         << "    ~~~~>After calculateState, "
         << "printing microPhaseVolumes" << endl;
    for (int i = 0; i < numMicroPhases_; ++i) { // microPhaseVolumes.size()
      cout << "    Phase name " << microPhaseName_[i]
           << ": volume = " << microPhaseVolumes[i] << endl;
      cout.flush();
    }
    // cout << "%%%%%%%%%% Printing GEM Masses and "
    //      << "Volumes in this Step %%%%%%%" << endl;

    // for (int myid = 0; myid < numGEMPhases_; myid++) {
    //   cout << "Mass and volume of GEM phase " << node_->pCSD()->PHNL[myid]
    //        << " = " << GEMPhaseMass_[myid] << " g and " <<
    //        GEMPhaseVolume_[myid]
    //        << " m3" << endl;
    // }
    // cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    //      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  }

  /// JWB 2023-07-29
  /// Assign the molar enthalpy of each DC

  for (int i = 0; i < numDCs_; i++) {
    DCH0_[i] = node_->DC_H0(i, P_, T_);
  }

  /// JWB 2023-04-14
  /// Account for subvoxel porosity within each microstructure phase here
  /// We used to do this in the Lattice class but it is cleaner here

  double phi; // local variable to store subvoxel volume fraction of pores of
              // a phase 0 <= phi <= 1

  int microPhaseMembersSize_i;
  for (int i = 1; i < numMicroPhases_; i++) {
    if (verbose_) {
      cout << "Setting microPhase amounts for " << i << " = "
           << microPhaseName_[i] << endl;
      cout.flush();
    }

    if (!isKinetic_[i]) {
      calcMicroPhasePorosity(i);
      // phi = getMicroPhasePorosity(i);
      phi = microPhasePorosity_[i];
      microPhaseMass_[i] = microPhaseVolume_[i] = 0.0;
      microPhaseMembersSize_i = microPhaseMembers_[i].size();

      // if (i == 1) cout << endl << "microPhaseMembersSize_i for i = 1 is : "
      //                  <<  microPhaseMembersSize_i << endl;

      for (int j = 0; j < microPhaseMembersSize_i; j++) {

        microPhaseMass_[i] += GEMPhaseMass_[microPhaseMembers_[i][j]];

        if (verbose_) {
          cout << "    Is NOT a KINETIC phase: is composed of "
               << GEMPhaseName_[microPhaseMembers_[i][j]]
               << " having mass = " << GEMPhaseMass_[microPhaseMembers_[i][j]]
               << " and volume = " << GEMPhaseVolume_[microPhaseMembers_[i][j]]
               << " and porosity = " << phi << endl;
          cout.flush();
        }

        /// Here is where the subvoxel porosity is included
        /// @todo Need a more disciplined approach to making sure
        /// these are solid phases and not capillary porosity

        if (i == ELECTROLYTEID) {
          microPhaseVolume_[i] += GEMPhaseVolume_[microPhaseMembers_[i][j]];
          // cout << endl << "   i = " << i << "   j = " << j
          //      << "   GEMPhaseName_ = " <<
          //      GEMPhaseName_[microPhaseMembers_[i][j]]
          //      << "   microPhaseVolume_[i](use GEMPhaseVolume_!!! ) = " <<
          //      microPhaseVolume_[i]
          //      << "   numMoles(use waterMollarVol_(waterDCId)!!!) = "
          //      << microPhaseVolume_[i]/waterMollarVol_ <<  endl;
          // cout << endl << "   DCMoles_[waterDCId_] = " <<
          // DCMoles_[waterDCId_]
          //      << "   waterMollarVol_ = " << waterMollarVol_ << endl;
        } else if (i != VOIDID && phi < 1.0) {
          microPhaseVolume_[i] +=
              (GEMPhaseVolume_[microPhaseMembers_[i][j]] / (1.0 - phi));
        } else {
          cout << "WARNING: A solid phase with porosity = 1.0?" << endl;
          cout.flush();
          microPhaseVolume_[i] +=
              (GEMPhaseVolume_[microPhaseMembers_[i][j]] / (0.001));
        }
      }
      microVolume_ += microPhaseVolume_[i];
      if (microPhaseVolume_[i] < 0) {
        cout << endl
             << "ChemicalSystem::calculateState error1 - microPhaseVolume_ < 0 "
                "for cyc = "
             << cyc << " :" << endl;
        for (int i = ELECTROLYTEID; i < numMicroPhases_; i++) {
          cout << "   " << i
               << " : phName/microPhaseVolume_ : " << microPhaseName_[i]
               << " / " << microPhaseVolume_[i] << endl;
        }
        cout << endl << "end program" << endl;

        throw GEMException("ChemicalSystem", "calculateState",
                           "error1 : microPhaseVolume_ < 0");
        // exit(0);
      }
    } else {
      calcMicroPhasePorosity(i);
      phi = getMicroPhasePorosity(i);
      if (verbose_) {
        cout << "    IS a KINETIC phase: is composed of "
             << GEMPhaseName_[microPhaseMembers_[i][0]]
             << "  having mass = " << microPhaseMass_[i]
             << "  and volume = " << microPhaseVolume_[i]
             << "  and internal porosity = " << phi << endl;
        cout.flush();
      }

      ///
      /// microPhaseMass and microPhaseVolume for kinetic phases are
      /// already set in KineticModel::calculateKineticStep
      ///

      microVolume_ += (microPhaseVolume_[i] / (1.0 - phi));
      if (microPhaseVolume_[i] < 0) {
        cout << endl
             << "ChemicalSystem::calculateState error2 - microPhaseVolume_ < 0 "
                "for cyc = "
             << cyc << " :" << endl;
        for (int i = ELECTROLYTEID; i < numMicroPhases_; i++) {
          cout << "   " << i
               << " : phName/microPhaseVolume_ : " << microPhaseName_[i]
               << " / " << microPhaseVolume_[i] << endl;
        }
        cout << endl << "end program" << endl;

        throw GEMException("ChemicalSystem", "calculateState",
                           "error2 : microPhaseVolume_ < 0");
        // exit(0);
      }
    }
  }

  scaledCementMass_ = 0;
  for (int i = 0; i < isParrotKillohSize_; i++) {
    scaledCementMass_ += microPhaseMass_[isParrotKilloh_[i]];
    // cout << "  ChemicalSystem::calculateState for cyc = " << cyc
    //      << "  i = " << i << "  isParrotKilloh_[i] = " << isParrotKilloh_[i]
    //      << "   "
    //      << microPhaseMass_[isParrotKilloh_[i]] << endl;
  }

  if (isSaturated_) { // System is saturated
    if (initMicroVolume_ > microVolume_) {
      // double water_molarv;
      double water_molesincr;
      // int wDCId = getDCId("H2O@");
      // water_molarv = node_->DC_V0(wDCId, P_, T_);
      water_molesincr = (initMicroVolume_ - microVolume_) / waterMollarVol_;
      if (verbose_) {
        cout << "System is saturated: wDCId = " << waterDCId_ << endl;
        cout << "    water_molarv = " << waterMollarVol_ << endl;
        cout << "    volume increase of water is: "
             << (initMicroVolume_ - microVolume_) << endl;
        cout << "    water_molesincr = " << water_molesincr << endl;
      }
      DCMoles_[waterDCId_] += water_molesincr;

      // double waterMolarMass = getDCMolarMass(wDCId);
      addWaterMassAndVolume(water_molesincr * waterMollarMass_,
                            initMicroVolume_ - microVolume_); // necessary

      cout << "  ChemicalSystem::calculateState - cyc = " << cyc
           << " : water_molesincr = " << water_molesincr << endl;
      // newMicroVolume_ = initMicroVolume_;
    }
  }

  // cout << endl << "ChemicalSystem::calculateState-2 waterVolume_ = "
  //      << microPhaseVolume_[1] << "   waterMoles_ = " << DCMoles_[waterDCId_]
  //      << endl;
  // cout << endl << "ChemicalSystem::calculateState initMicroVolume_ = " <<
  // initMicroVolume_ << endl; double waterDensity = waterMollarMass_ /
  // waterMollarVol_ / 1.0e6; // g/cm3 double waterTotMass_0 =
  // (microPhaseVolume_[1]/initMicroVolume_) * waterDensity * 100 /
  // initSolidMass_; double waterTotMoles_0 = waterTotMass_0 / waterMollarMass_;
  // // mol cout << endl << "ChemicalSystem::calculateState waterMolesCalc_M = "
  // << waterTotMoles_0 << endl; cout << endl << "ChemicalSystem::calculateState
  // waterMolesCalc_V = " <<
  //         microPhaseVolume_[1] / waterMollarVol_ << endl;

  if (verbose_) {
    cout << "GEM volume change = "
         << 100.0 * (microVolume_ - initMicroVolume_) / (initMicroVolume_)
         << " %" << endl;
    cout.flush();
  }

  if (initMicroVolume_ < microVolume_) {
    cout << "  ChemicalSystem::calculateState - cyc = " << cyc
         << " => initMicroVolume_ < microVolume_ : " << initMicroVolume_
         << " < " << microVolume_ << endl;
    // initMicroVolume_ = microVolume_;
    // newMicroVolume_ = microVolume_;
  }

  if (verbose_) {
    cout << "Leaving ChemicalSystem::calculateState now" << endl;
    cout.flush();
  }

  // setMicroPhaseSI(time);
  setMicroPhaseSI();

  return timesGEMFailed_;
}

void ChemicalSystem::setMicroPhaseSI() {

  microPhaseSI_.clear();
  microPhaseSI_.resize(numMicroPhases_, 0.0);

  try {
    double aveSI = 0.0;
    double moles = 0.0;
    double tmoles = 0.0;
    vector<int> microPhaseDCMembers;
    int sizeMicroPhaseDCMembers;
    string pname;
    int newDCId;

    // Query CSD node to set the SI of every microPhase
    // if (isFirst) {
    //} else {

    // setSI();
    for (int i = FIRST_SOLID; i < numMicroPhases_; ++i) {
      pname = microPhaseName_[i];
      aveSI = moles = 0.0;
      microPhaseDCMembers = getMicroPhaseDCMembers(i);
      sizeMicroPhaseDCMembers = microPhaseDCMembers.size();
      // cout << endl << "   " << i << "\tpname: " << pname
      //      << "\tmicroPhaseMembers.size: " << sizeMicroPhaseDCMembers << "
      //      : " << endl;
      for (int ii = 0; ii < sizeMicroPhaseDCMembers; ++ii) {
        newDCId = microPhaseDCMembers.at(ii);
        tmoles = DCMoles_[newDCId];
        // if ( microPhaseDCMembers.size() == 1) {
        //   vector<int> microPhaseMembers = getMicroPhaseMembers(i);
        //   int newGEMPhaseId = microPhaseMembers.at(0);
        //   cout << "       " << ii << "\tcyc/newDCId: " << cyc << " / "
        //        << newDCId << "\tnewDCName: " << DCName_[newDCId]
        //        << "\ttmoles: " << tmoles
        //        << "\tDC_a: " << node_->DC_a(newDCId)
        //        << "\tmono -> newGEMPhaseId: " << newGEMPhaseId
        //        << "\tSI_: " << SI_[newGEMPhaseId] << endl;
        // } else {
        //   cout << "       " << ii << "\tcyc/newDCId: " << cyc << " / "
        //        << newDCId << "\tnewDCName: " << DCName_[newDCId]
        //        << "\ttmoles: " << tmoles
        //        << "\tDC_a: " << node_->DC_a(newDCId) << endl;
        // }

        aveSI += (node_->DC_a(newDCId) * tmoles);
        moles += tmoles;
      }
      // cout << "          aveSI: " << aveSI << "\tmoles: " << moles << endl;
      if (moles > 0.0) {
        aveSI = aveSI / moles;
      }
      microPhaseSI_.at(i) = aveSI;
      // cout << "          pname = " << pname << "  =>     microPhaseSI_(cyc
      // = " << cyc
      //      << ") = " << microPhaseSI_[i] << "\tmoles: " << moles << endl;
    }
    //} //if (isFirst) {
    // cout << endl << "ChemicalSystem::setMicroPhaseSI end" << endl; exit(0);
  } catch (EOBException eex) {
    eex.printException();
    exit(1);
  }

  return;
}

void ChemicalSystem::setMicroPhaseSI(double time) {

  microPhaseSI_.clear();
  microPhaseSI_.resize(numMicroPhases_, 0.0);
  double sumSI = 0.0;

  try {
    double aveSI = 0.0;
    double moles = 0.0;
    vector<int> microPhaseDCMembers;
    string pname;

    // Query CSD node to set the SI of every microPhase
    // if (isFirst) {
    //} else {

    // setSI();

    //***********
    // solutPhaseMoles_
    double totMoles = 0;
    vector<int> microPhasePhMembers;
    int sizeMicroPhasePhMembers;
    int newPhId;

    if (time >= beginAttackTime_) {
      cout << endl
           << "          SI calc from GEM phases for mPhId = 15 - mPhName = "
           << microPhaseName_[15] << " : " << endl;
    }

    for (int i = FIRST_SOLID; i < numMicroPhases_; ++i) {
      pname = microPhaseName_[i];
      aveSI = totMoles = 0.0;
      microPhasePhMembers =
          microPhaseMembers_[i]; // getMicroPhaseMembers(i);
                                 // //microPhaseMembers_ vs GEMPhaseIdLookup_
      sizeMicroPhasePhMembers = (int)microPhasePhMembers.size();
      if ((time >= beginAttackTime_) && (i == 15)) {
        for (int ii = 0; ii < sizeMicroPhasePhMembers; ++ii) {
          newPhId = microPhaseMembers_[i][ii]; // microPhasePhMembers[i];
          moles = solutPhaseMoles_[newPhId];
          cout << "           ii/newPhId/moles/Ph_SatInd/SI/newPhName = " << ii
               << " / " << newPhId << " / " << moles << " / "
               << node_->Ph_SatInd(newPhId) << " / "
               << pow(10, node_->Ph_SatInd(newPhId)) << " / "
               << GEMPhaseName_[newPhId] << endl;
        }
      }
      for (int ii = 0; ii < sizeMicroPhasePhMembers; ++ii) {
        newPhId = microPhaseMembers_[i][ii]; // microPhasePhMembers[i];
        moles = solutPhaseMoles_[newPhId];
        aveSI += (pow(10, node_->Ph_SatInd(newPhId)) * moles);
        totMoles += moles;
      }
      if (totMoles > 0.0) {
        aveSI = aveSI / totMoles;
      }
      microPhaseSI_.at(i) = aveSI;
    }
    if (time >= beginAttackTime_)
      cout << "          aveSI for mPhId = 15 - mPhName = "
           << microPhaseName_[15] << " : " << microPhaseSI_.at(15) << endl;
    //***********

    // if (time >= beginAttackTime_) {
    //   setSI();
    //   cout << endl << "          SI calc from GEM DCs for mPhId = 15 -
    //   mPhName = "
    //        << microPhaseName_[15] << " : " << endl;
    // }

    /*
    for (int i = FIRST_SOLID; i < numMicroPhases_; ++i) {
      pname = microPhaseName_[i];
      aveSI = moles = 0.0;
      microPhaseDCMembers = getMicroPhaseDCMembers(i);
      sizeMicroPhaseDCMembers = microPhaseDCMembers.size();
      // cout << endl << "   " << i << "\tpname: " << pname
      //      << "\tmicroPhaseMembers.size: " << sizeMicroPhaseDCMembers << "
      //      : " << endl;
      if ((time >= beginAttackTime_) && (i == 15)) {
        double *Falp = (node_->ppmm())->Falp;
        for (int ii = 0; ii < sizeMicroPhaseDCMembers; ++ii) {
          newDCId = microPhaseDCMembers.at(ii);
          tmoles = DCMoles_[newDCId];
          if (tmoles > 1.e-15) sumSI += node_->DC_a(newDCId);
          cout << "           ii/newDCId/tmoles/Falp/SI(Falp)/SI(DC_a)/newDCName
    = "
               << ii << " / " << newDCId << " / " << tmoles << " / "
               << Falp[newDCId] <<   " / " << pow(10, Falp[newDCId]) << " / "
               << node_->DC_a(newDCId) << " / " << DCName_[newDCId] << endl;
        }
      }
      for (int ii = 0; ii < sizeMicroPhaseDCMembers; ++ii) {
        newDCId = microPhaseDCMembers.at(ii);
        tmoles = DCMoles_[newDCId];
        // if ( microPhaseDCMembers.size() == 1) {
        //   vector<int> microPhaseMembers = getMicroPhaseMembers(i);
        //   int newGEMPhaseId = microPhaseMembers.at(0);
        //   cout << "       " << ii << "\tcyc/newDCId: " << cyc << " / "
        //        << newDCId << "\tnewDCName: " << DCName_[newDCId]
        //        << "\ttmoles: " << tmoles
        //        << "\tDC_a: " << node_->DC_a(newDCId)
        //        << "\tmono -> newGEMPhaseId: " << newGEMPhaseId
        //        << "\tSI_: " << SI_[newGEMPhaseId] << endl;
        // } else {
        //   cout << "       " << ii << "\tcyc/newDCId: " << cyc << " / "
        //        << newDCId << "\tnewDCName: " << DCName_[newDCId]
        //        << "\ttmoles: " << tmoles
        //        << "\tDC_a: " << node_->DC_a(newDCId) << endl;
        // }

        aveSI += (node_->DC_a(newDCId) * tmoles);
        moles += tmoles;
      }
      // cout << "          aveSI: " << aveSI << "\tmoles: " << moles << endl;
      if (moles > 0.0) {
        aveSI = aveSI / moles;
      }
      microPhaseSI_.at(i) = aveSI;
      // cout << "          pname = " << pname << "  =>     microPhaseSI_(cyc
      // = " << cyc
      //      << ") = " << microPhaseSI_[i] << "\tmoles: " << moles << endl;
    }
    //} //if (isFirst) {
    // cout << endl << "ChemicalSystem::setMicroPhaseSI end" << endl; exit(0);
    */
  } catch (EOBException eex) {
    eex.printException();
    exit(1);
  }

  if (time >= beginAttackTime_)
    cout << "          aveSI for mPhId = 15 - mPhName = " << microPhaseName_[15]
         << " : " << microPhaseSI_.at(15) << "   &   sumSI = " << sumSI << endl;

  return;
}

void ChemicalSystem::initElasticModuliMap(void) {
  // (1): Table 1 from Haecker et al. / Cement and Concrete Research
  // 35(2005)1948–1960 (k,G) -> (E,n) : E = 9KG/(3K + G)  &  n = (3K - 2G)/(2(3K
  // + G)) (E,n) -> (k,G) : K = E/(3(1 - 2n)) &  G = E/(2(1 + n))

  // struct elMod {
  //   float K; // bulk modulus    : K(GPa)
  //   float G; // shear modulus   : G(GPa)
  //   float E; // Young's modulus : E(GPa)
  //   float n; // Poisson's ratio
  //   string altName;
  // };
  // map<string, elMod> elasticModuli_;

  // cement chemistry notation :
  //   C=CaO, S=SiO2,A=Al2O3, F=Fe2O3, H=H2O, K=K2O, N=Na2O, M=MgO
  //   s=SO3, c=CO2

  elasticModuli_["Void"].K = 0.0;
  elasticModuli_["Void"].G = 0.0;
  elasticModuli_["Void"].E = 0.0;
  elasticModuli_["Void"].n = 0.0;
  elasticModuli_["Void"].altName = "Empty";

  // (1)
  elasticModuli_["Electrolyte"].K = 2.2;
  elasticModuli_["Electrolyte"].G = 0.0;
  elasticModuli_["Electrolyte"].E = 0.0;
  elasticModuli_["Electrolyte"].n = 0.5;
  elasticModuli_["Electrolyte"].altName = "Porosity";

  // (1)
  elasticModuli_["Alite"].K = 105.2;
  elasticModuli_["Alite"].G = 44.80;
  elasticModuli_["Alite"].E = 117.6;
  elasticModuli_["Alite"].n = 0.314;
  elasticModuli_["Alite"].altName = "C3S";

  // (1) - same as Alite
  elasticModuli_["Belite"].K = 105.2;
  elasticModuli_["Belite"].G = 44.80;
  elasticModuli_["Belite"].E = 117.6;
  elasticModuli_["Belite"].n = 0.314;
  elasticModuli_["Belite"].altName = "C2S";

  // (1) - same as Alite
  elasticModuli_["Aluminate"].K = 105.2;
  elasticModuli_["Aluminate"].G = 44.80;
  elasticModuli_["Aluminate"].E = 117.6;
  elasticModuli_["Aluminate"].n = 0.314;
  elasticModuli_["Aluminate"].altName = "C3A";

  // (1) - same as Alite
  elasticModuli_["Ferrite"].K = 105.2;
  elasticModuli_["Ferrite"].G = 44.80;
  elasticModuli_["Ferrite"].E = 117.6;
  elasticModuli_["Ferrite"].n = 0.314;
  elasticModuli_["Ferrite"].altName = "C4AF";

  // (1)
  elasticModuli_["Gypsum"].K = 42.5;
  elasticModuli_["Gypsum"].G = 15.7;
  elasticModuli_["Gypsum"].E = 45.7;
  elasticModuli_["Gypsum"].n = 0.33;
  elasticModuli_["Gypsum"].altName = "CsH2, CaSO4*2H2O, Dihydrate";

  // (1)
  elasticModuli_["Bassanite"].K = 52.4;
  elasticModuli_["Bassanite"].G = 24.2;
  elasticModuli_["Bassanite"].E = 62.9;
  elasticModuli_["Bassanite"].n = 0.30;
  elasticModuli_["Bassanite"].altName = "CsH1/2, CaSO4*0.5H2O, Hemihydrate";

  // (1)
  elasticModuli_["Anhydrite"].K = 54.90;
  elasticModuli_["Anhydrite"].G = 29.30;
  elasticModuli_["Anhydrite"].E = 80.00;
  elasticModuli_["Anhydrite"].n = 0.275;
  elasticModuli_["Anhydrite"].altName = "Cs, CaSO4";

  // (1)
  elasticModuli_["Arcanite"].K = 31.90;
  elasticModuli_["Arcanite"].G = 17.40;
  elasticModuli_["Arcanite"].E = 44.20;
  elasticModuli_["Arcanite"].n = 0.269;
  elasticModuli_["Arcanite"].altName = "Ks, K2SO4";

  // (1)
  elasticModuli_["Thenardite"].K = 43.40;
  elasticModuli_["Thenardite"].G = 22.30;
  elasticModuli_["Thenardite"].E = 57.10;
  elasticModuli_["Thenardite"].n = 0.281;
  elasticModuli_["Thenardite"].altName = "Ns, Na2SO4";

  // (1)
  elasticModuli_["SilicaAm"].K = 36.50;
  elasticModuli_["SilicaAm"].G = 31.20;
  elasticModuli_["SilicaAm"].E = 72.80;
  elasticModuli_["SilicaAm"].n = 0.167;
  elasticModuli_["SilicaAm"].altName = "S, SiO2, Silica fume, Amorphous silica";

  // (1)
  elasticModuli_["Portlandite"].K = 40.00;
  elasticModuli_["Portlandite"].G = 16.00;
  elasticModuli_["Portlandite"].E = 42.30;
  elasticModuli_["Portlandite"].n = 0.324;
  elasticModuli_["Portlandite"].altName = "CH, Ca(OH)2";

  // (1)
  elasticModuli_["CSHQ"].K = 14.9;
  elasticModuli_["CSHQ"].G = 9.00;
  elasticModuli_["CSHQ"].E = 22.4;
  elasticModuli_["CSHQ"].n = 0.25;
  elasticModuli_["CSHQ"].altName = "CSH, C-S-H";

  // (1)
  elasticModuli_["Calcite"].K = 69.8;
  elasticModuli_["Calcite"].G = 30.4;
  elasticModuli_["Calcite"].E = 79.6;
  elasticModuli_["Calcite"].n = 0.31;
  elasticModuli_["Calcite"].altName =
      "Cc, CaCO3, Limestone, Aragonite, Vaterite";

  // (1) - same as CSHQ
  elasticModuli_["Hydrogarnet"].K = 14.9;
  elasticModuli_["Hydrogarnet"].G = 9.00;
  elasticModuli_["Hydrogarnet"].E = 22.4;
  elasticModuli_["Hydrogarnet"].n = 0.25;
  elasticModuli_["Hydrogarnet"].altName = "C3AH6, Hydrogarnet";

  // (1) - same as CSHQ
  elasticModuli_["AFt"].K = 14.9;
  elasticModuli_["AFt"].G = 9.00;
  elasticModuli_["AFt"].E = 22.4;
  elasticModuli_["AFt"].n = 0.25;
  elasticModuli_["AFt"].altName = "C6As3H32, Ettringite";

  // (1) - same as Portlandite
  elasticModuli_["AFm"].K = 40.00;
  elasticModuli_["AFm"].G = 16.00;
  elasticModuli_["AFm"].E = 42.30;
  elasticModuli_["AFm"].n = 0.234;
  elasticModuli_["AFm"].altName = "C4AsH12, AFm";

  // (1) - same as CSHQ
  elasticModuli_["FH3"].K = 14.9;
  elasticModuli_["FH3"].G = 9.00;
  elasticModuli_["FH3"].E = 22.4;
  elasticModuli_["FH3"].n = 0.25;
  elasticModuli_["FH3"].altName = "FH3, Iron hydroxide";

  // (1) - same as Portlandite
  elasticModuli_["CaCl2"].K = 40.00;
  elasticModuli_["CaCl2"].G = 16.00;
  elasticModuli_["CaCl2"].E = 42.30;
  elasticModuli_["CaCl2"].n = 0.324;
  elasticModuli_["CaCl2"].altName = "Calcium chloride";

  // (1) - same as CSHQ
  elasticModuli_["Hydrocalumite"].K = 14.9;
  elasticModuli_["Hydrocalumite"].G = 9.00;
  elasticModuli_["Hydrocalumite"].E = 22.4;
  elasticModuli_["Hydrocalumite"].n = 0.25;
  elasticModuli_["Hydrocalumite"].altName = "C3A(CaCl2)H10, Friedel salt";

  // (1) - same as CSHQ
  elasticModuli_["Stratlingite"].K = 14.9;
  elasticModuli_["Stratlingite"].G = 9.00;
  elasticModuli_["Stratlingite"].E = 22.4;
  elasticModuli_["Stratlingite"].n = 0.25;
  elasticModuli_["Stratlingite"].altName = "C2ASH8, Stratlingite";

  // (1) - same as AFm
  elasticModuli_["AFmc"].K = 40.00;
  elasticModuli_["AFmc"].G = 16.00;
  elasticModuli_["AFmc"].E = 42.30;
  elasticModuli_["AFmc"].n = 0.234;
  elasticModuli_["AFmc"].altName =
      "C4AcH11, AFm-c, Monoarbonate, Carboaluminate";

  // elasticModuli_["Monosulfate"].K = ;
  // elasticModuli_["Monosulfate"].G = ;
  // elasticModuli_["Monosulfate"].E = ;
  // elasticModuli_["Monosulfate"].n = ;
  // elasticModuli_["Monosulfate"].altName = "Sulfoaluminate";

  // gitCheckout375c4798 * - as Portlandite
  elasticModuli_["Brucite"].K = 40.00;
  elasticModuli_["Brucite"].G = 16.00;
  elasticModuli_["Brucite"].E = 42.30; // *
  elasticModuli_["Brucite"].n = 0.324; // *
  elasticModuli_["Brucite"].altName = "MH, Mg(OH)2";

  // gitCheckout375c4798 * - as Alite
  elasticModuli_["CaO"].K = 105.2;
  elasticModuli_["CaO"].G = 44.80;
  elasticModuli_["CaO"].E = 117.6; // *
  elasticModuli_["CaO"].n = 0.314; // *
  elasticModuli_["CaO"].altName = "C, Free lime";

  // gitCheckout375c4798 * - as CSHQ
  elasticModuli_["Hydrotalcite"].K = 14.9;
  elasticModuli_["Hydrotalcite"].G = 9.00;
  elasticModuli_["Hydrotalcite"].E = 22.4; // *
  elasticModuli_["Hydrotalcite"].n = 0.25; // *
  elasticModuli_["Hydrotalcite"].altName = "Hydrotalcite";

  // gitCheckout375c4798 * - as Alite
  elasticModuli_["Periclase"].K = 105.2;
  elasticModuli_["Periclase"].G = 44.80;
  elasticModuli_["Periclase"].E = 117.0; // *
  elasticModuli_["Periclase"].n = 0.314; // *
  elasticModuli_["Periclase"].altName = "M, MgO";

  // gitCheckout375c4798 * - as CSHQ
  elasticModuli_["Damage"].K = 14.9;
  elasticModuli_["Damage"].G = 9.00;
  elasticModuli_["Damage"].E = 22.4; // *
  elasticModuli_["Damage"].n = 0.25; // *
  elasticModuli_["Damage"].altName = "Damage";

  /*
  elasticModuli_["C2AS"].K = ;
  elasticModuli_["C2AS"].G = ;
  elasticModuli_["C2AS"].E = ;
  elasticModuli_["C2AS"].n = ;
  elasticModuli_["C2AS"].altName = "C2AS(am)";

  elasticModuli_["CA2S"].K = ;
  elasticModuli_["CA2S"].G = ;
  elasticModuli_["CA2S"].E = ;
  elasticModuli_["CA2S"].n = ;
  elasticModuli_["CA2S"].altName = "CA2S(am)";

  elasticModuli_["K6A2S"].K = ;
  elasticModuli_["K6A2S"].G = ;
  elasticModuli_["K6A2S"].E = ;
  elasticModuli_["K6A2S"].n = ;
  elasticModuli_["K6A2S"].altName = "K6A2S(am)";

  elasticModuli_["Diopside"].K = ;
  elasticModuli_[""].G = ;
  elasticModuli_[""].E = ;
  elasticModuli_[""].n = ;
  elasticModuli_["Diopside"].altName = "Pyroxene";

  elasticModuli_["Forsterite"].K = ;
  elasticModuli_["Forsterite"].G = ;
  elasticModuli_["Forsterite"].E = ;
  elasticModuli_["Forsterite"].n = ;
  elasticModuli_["Forsterite"].altName = "Olivine";

  elasticModuli_["Goethite"].K = ;
  elasticModuli_["Goethite"].G = ;
  elasticModuli_["Goethite"].E = ;
  elasticModuli_["Goethite"].n = ;
  elasticModuli_["Goethite"].altName = "Goethite";

  elasticModuli_["Mullite"].K = ;
  elasticModuli_["Mullite"].G = ;
  elasticModuli_["Mullite"].E = ;
  elasticModuli_["Mullite"].n = ;
  elasticModuli_["Mullite"].altName = "Mullite";

  elasticModuli_["Pyrite"].K = ;
  elasticModuli_["Pyrite"].G = ;
  elasticModuli_["Pyrite"].E = ;
  elasticModuli_["Pyrite"].n = ;
  elasticModuli_["Pyrite"].altName = "FeS";

  elasticModuli_["Pyrrhotite"].K = ;
  elasticModuli_["Pyrrhotite"].G = ;
  elasticModuli_["Pyrrhotite"].E = ;
  elasticModuli_["Pyrrhotite"].n = ;
  elasticModuli_["Pyrrhotite"].altName = "Pyrrhotite";

  elasticModuli_["Quartz"].K = ;
  elasticModuli_["Quartz"].G = ;
  elasticModuli_["Quartz"].E = ;
  elasticModuli_["Quartz"].n = ;
  elasticModuli_["Quartz"].altName = "Silica";

  elasticModuli_["Sodalite"].K = ;
  elasticModuli_["Sodalite"].G = ;
  elasticModuli_["Sodalite"].E = ;
  elasticModuli_["Sodalite"].n = ;
  elasticModuli_["Sodalite"].altName = "Sodalite";

  elasticModuli_["Syngenite"].K = ;
  elasticModuli_["Syngenite"].G = ;
  elasticModuli_["Syngenite"].E = ;
  elasticModuli_["Syngenite"].n = ;
  elasticModuli_["Syngenite"].altName = "Syngenite";

  elasticModuli_["Troilite"].K = ;
  elasticModuli_["Troilite"].G = ;
  elasticModuli_["Troilite"].E = ;
  elasticModuli_["Troilite"].n = ;
  elasticModuli_["Troilite"].altName = "Troilite";

  elasticModuli_["Zeolite"].K = ;
  elasticModuli_["Zeolite"].G = ;
  elasticModuli_["Zeolite"].E = ;
  elasticModuli_["Zeolite"].n = ;
  elasticModuli_["Zeolite"].altName = "Zeolite";

  elasticModuli_["Dolomite"].K = ;
  elasticModuli_["Dolomite"].G = ;
  elasticModuli_["Dolomite"].E = ;
  elasticModuli_["Dolomite"].n = ;
  elasticModuli_["Dolomite"].altName = "Dolomite";
  */
}

void ChemicalSystem::initColorMap(void) {
  // struct elemColor {
  //   int colorId;
  //   string altName;
  //   vector<int> rgb;
  //   int gray;
  // };

  colorN_["Void"].colorId = 0;
  colorN_["Void"].altName = "Empty";
  colorN_["Void"].rgb.push_back(0);
  colorN_["Void"].rgb.push_back(0);
  colorN_["Void"].rgb.push_back(0);
  colorN_["Void"].gray = 0;

  colorN_["Electrolyte"].colorId = 1;
  colorN_["Electrolyte"].altName = "Porosity";
  colorN_["Electrolyte"].rgb.push_back(0);
  colorN_["Electrolyte"].rgb.push_back(25);
  colorN_["Electrolyte"].rgb.push_back(25);
  colorN_["Electrolyte"].gray = 25;

  colorN_["AFm"].colorId = 2;
  colorN_["AFm"].altName = "AFm";
  colorN_["AFm"].rgb.push_back(244);
  colorN_["AFm"].rgb.push_back(70);
  colorN_["AFm"].rgb.push_back(203);
  colorN_["AFm"].gray = 106;

  colorN_["AFmc"].colorId = 3;
  colorN_["AFmc"].altName = "AFm-c, Carboaluminate";
  colorN_["AFmc"].rgb.push_back(250);
  colorN_["AFmc"].rgb.push_back(198);
  colorN_["AFmc"].rgb.push_back(220);
  colorN_["AFmc"].gray = 108;

  colorN_["AFt"].colorId = 4;
  colorN_["AFt"].altName = "Ettringite";
  colorN_["AFt"].rgb.push_back(127);
  colorN_["AFt"].rgb.push_back(0);
  colorN_["AFt"].rgb.push_back(255);
  colorN_["AFt"].gray = 113;

  colorN_["Alite"].colorId = 5;
  colorN_["Alite"].altName = "C3S";
  colorN_["Alite"].rgb.push_back(42);
  colorN_["Alite"].rgb.push_back(42);
  colorN_["Alite"].rgb.push_back(210);
  colorN_["Alite"].gray = 220;

  colorN_["Aluminate"].colorId = 6;
  colorN_["Aluminate"].altName = "C3A";
  colorN_["Aluminate"].rgb.push_back(178);
  colorN_["Aluminate"].rgb.push_back(178);
  colorN_["Aluminate"].rgb.push_back(178);
  colorN_["Aluminate"].gray = 195;

  colorN_["Anhydrite"].colorId = 7;
  colorN_["Anhydrite"].altName = "Anhydrite";
  colorN_["Anhydrite"].rgb.push_back(255);
  colorN_["Anhydrite"].rgb.push_back(255);
  colorN_["Anhydrite"].rgb.push_back(128);
  colorN_["Anhydrite"].gray = 140;

  colorN_["Arcanite"].colorId = 8;
  colorN_["Arcanite"].altName = "K2SO4";
  colorN_["Arcanite"].rgb.push_back(255);
  colorN_["Arcanite"].rgb.push_back(0);
  colorN_["Arcanite"].rgb.push_back(0);
  colorN_["Arcanite"].gray = 100;

  colorN_["Bassanite"].colorId = 9;
  colorN_["Bassanite"].altName = "Hemihydrate";
  colorN_["Bassanite"].rgb.push_back(255);
  colorN_["Bassanite"].rgb.push_back(240);
  colorN_["Bassanite"].rgb.push_back(86);
  colorN_["Bassanite"].gray = 140;

  colorN_["Belite"].colorId = 10;
  colorN_["Belite"].altName = "C2S";
  colorN_["Belite"].rgb.push_back(139);
  colorN_["Belite"].rgb.push_back(79);
  colorN_["Belite"].rgb.push_back(19);
  colorN_["Belite"].gray = 200;

  colorN_["C2AS"].colorId = 11;
  colorN_["C2AS"].altName = "C2AS(am)";
  colorN_["C2AS"].rgb.push_back(255);
  colorN_["C2AS"].rgb.push_back(165);
  colorN_["C2AS"].rgb.push_back(0);
  colorN_["C2AS"].gray = 110;

  colorN_["CA2S"].colorId = 12;
  colorN_["CA2S"].altName = "CA2S(am)";
  colorN_["CA2S"].rgb.push_back(255);
  colorN_["CA2S"].rgb.push_back(192);
  colorN_["CA2S"].rgb.push_back(65);
  colorN_["CA2S"].gray = 108;

  colorN_["Calcite"].colorId = 13;
  colorN_["Calcite"].altName = "Limestone, CaCO3, Aragonite, Vaterite";
  colorN_["Calcite"].rgb.push_back(0);
  colorN_["Calcite"].rgb.push_back(204);
  colorN_["Calcite"].rgb.push_back(0);
  colorN_["Calcite"].gray = 121;

  colorN_["CaO"].colorId = 14;
  colorN_["CaO"].altName = "Free lime";
  colorN_["CaO"].rgb.push_back(0);
  colorN_["CaO"].rgb.push_back(230);
  colorN_["CaO"].rgb.push_back(0);
  colorN_["CaO"].gray = 123;

  colorN_["CSHQ"].colorId = 15;
  colorN_["CSHQ"].altName = "CSH, C-S-H";
  colorN_["CSHQ"].rgb.push_back(245);
  colorN_["CSHQ"].rgb.push_back(222);
  colorN_["CSHQ"].rgb.push_back(179);
  colorN_["CSHQ"].gray = 159;

  colorN_["Diopside"].colorId = 16;
  colorN_["Diopside"].altName = "Pyroxene";
  colorN_["Diopside"].rgb.push_back(150);
  colorN_["Diopside"].rgb.push_back(150);
  colorN_["Diopside"].rgb.push_back(0);
  colorN_["Diopside"].gray = 157;

  colorN_["Ferrite"].colorId = 17;
  colorN_["Ferrite"].altName = "C4AF";
  colorN_["Ferrite"].rgb.push_back(253);
  colorN_["Ferrite"].rgb.push_back(253);
  colorN_["Ferrite"].rgb.push_back(253);
  colorN_["Ferrite"].gray = 245;

  colorN_["Forsterite"].colorId = 18;
  colorN_["Forsterite"].altName = "Olivine";
  colorN_["Forsterite"].rgb.push_back(127);
  colorN_["Forsterite"].rgb.push_back(127);
  colorN_["Forsterite"].rgb.push_back(0);
  colorN_["Forsterite"].gray = 165;

  colorN_["Goethite"].colorId = 19;
  colorN_["Goethite"].altName = "Goethite";
  colorN_["Goethite"].rgb.push_back(210);
  colorN_["Goethite"].rgb.push_back(210);
  colorN_["Goethite"].rgb.push_back(210);
  colorN_["Goethite"].gray = 200;

  colorN_["Gypsum"].colorId = 20;
  colorN_["Gypsum"].altName = "Dihydrate";
  colorN_["Gypsum"].rgb.push_back(255);
  colorN_["Gypsum"].rgb.push_back(255);
  colorN_["Gypsum"].rgb.push_back(0);
  colorN_["Gypsum"].gray = 112;

  colorN_["Hydrotalcite"].colorId = 21;
  colorN_["Hydrotalcite"].altName = "Hydrotalcite";
  colorN_["Hydrotalcite"].rgb.push_back(0);
  colorN_["Hydrotalcite"].rgb.push_back(200);
  colorN_["Hydrotalcite"].rgb.push_back(200);
  colorN_["Hydrotalcite"].gray = 109;

  colorN_["K6A2S"].colorId = 22;
  colorN_["K6A2S"].altName = "K6A2S(am)";
  colorN_["K6A2S"].rgb.push_back(255);
  colorN_["K6A2S"].rgb.push_back(170);
  colorN_["K6A2S"].rgb.push_back(128);
  colorN_["K6A2S"].gray = 115;

  colorN_["Monosulfate"].colorId = 23;
  colorN_["Monosulfate"].altName = "Sulfoaluminate";
  colorN_["Monosulfate"].rgb.push_back(152);
  colorN_["Monosulfate"].rgb.push_back(93);
  colorN_["Monosulfate"].rgb.push_back(175);
  colorN_["Monosulfate"].gray = 106;

  colorN_["Mullite"].colorId = 24;
  colorN_["Mullite"].altName = "Mullite";
  colorN_["Mullite"].rgb.push_back(178);
  colorN_["Mullite"].rgb.push_back(34);
  colorN_["Mullite"].rgb.push_back(34);
  colorN_["Mullite"].gray = 128;

  colorN_["Portlandite"].colorId = 25;
  colorN_["Portlandite"].altName = "CH, Ca(OH)2";
  colorN_["Portlandite"].rgb.push_back(7);
  colorN_["Portlandite"].rgb.push_back(72);
  colorN_["Portlandite"].rgb.push_back(142);
  colorN_["Portlandite"].gray = 186;

  colorN_["Pyrite"].colorId = 26;
  colorN_["Pyrite"].altName = "FeS";
  colorN_["Pyrite"].rgb.push_back(212);
  colorN_["Pyrite"].rgb.push_back(175);
  colorN_["Pyrite"].rgb.push_back(55);
  colorN_["Pyrite"].gray = 235;

  colorN_["Pyrrhotite"].colorId = 27;
  colorN_["Pyrrhotite"].altName = "Pyrrhotite";
  colorN_["Pyrrhotite"].rgb.push_back(49);
  colorN_["Pyrrhotite"].rgb.push_back(142);
  colorN_["Pyrrhotite"].rgb.push_back(7);
  colorN_["Pyrrhotite"].gray = 240;

  colorN_["Quartz"].colorId = 28;
  colorN_["Quartz"].altName = "Silica";
  colorN_["Quartz"].rgb.push_back(32);
  colorN_["Quartz"].rgb.push_back(141);
  colorN_["Quartz"].rgb.push_back(142);
  colorN_["Quartz"].gray = 100;

  colorN_["SilicaAm"].colorId = 29;
  colorN_["SilicaAm"].altName = "Amorphous silica";
  colorN_["SilicaAm"].rgb.push_back(40);
  colorN_["SilicaAm"].rgb.push_back(173);
  colorN_["SilicaAm"].rgb.push_back(175);
  colorN_["SilicaAm"].gray = 100;

  colorN_["Sodalite"].colorId = 30;
  colorN_["Sodalite"].altName = "Sodalite";
  colorN_["Sodalite"].rgb.push_back(100);
  colorN_["Sodalite"].rgb.push_back(100);
  colorN_["Sodalite"].rgb.push_back(255);
  colorN_["Sodalite"].gray = 100;

  colorN_["Syngenite"].colorId = 31;
  colorN_["Syngenite"].altName = "Syngenite";
  colorN_["Syngenite"].rgb.push_back(255);
  colorN_["Syngenite"].rgb.push_back(192);
  colorN_["Syngenite"].rgb.push_back(203);
  colorN_["Syngenite"].gray = 210;

  colorN_["Thenardite"].colorId = 32;
  colorN_["Thenardite"].altName = "Na2SO4";
  colorN_["Thenardite"].rgb.push_back(255);
  colorN_["Thenardite"].rgb.push_back(20);
  colorN_["Thenardite"].rgb.push_back(0);
  colorN_["Thenardite"].gray = 95;

  colorN_["Troilite"].colorId = 33;
  colorN_["Troilite"].altName = "Troilite";
  colorN_["Troilite"].rgb.push_back(255);
  colorN_["Troilite"].rgb.push_back(178);
  colorN_["Troilite"].rgb.push_back(220);
  colorN_["Troilite"].gray = 232;

  colorN_["Zeolite"].colorId = 34;
  colorN_["Zeolite"].altName = "Zeolite";
  colorN_["Zeolite"].rgb.push_back(130);
  colorN_["Zeolite"].rgb.push_back(130);
  colorN_["Zeolite"].rgb.push_back(255);
  colorN_["Zeolite"].gray = 232;

  colorN_["Dolomite"].colorId = 35;
  colorN_["Dolomite"].altName = "Dolomite";
  colorN_["Dolomite"].rgb.push_back(255);
  colorN_["Dolomite"].rgb.push_back(215);
  colorN_["Dolomite"].rgb.push_back(0);
  colorN_["Dolomite"].gray = 150;

  colorN_["Brucite"].colorId = 36;
  colorN_["Brucite"].altName = "Mg(OH)2";
  colorN_["Brucite"].rgb.push_back(26);
  colorN_["Brucite"].rgb.push_back(100);
  colorN_["Brucite"].rgb.push_back(26);
  colorN_["Brucite"].gray = 83;

  colorN_["Periclase"].colorId = 37;
  colorN_["Periclase"].altName = "MgO";
  colorN_["Periclase"].rgb.push_back(255);
  colorN_["Periclase"].rgb.push_back(128);
  colorN_["Periclase"].rgb.push_back(128);
  colorN_["Periclase"].gray = 220; // as Alite

  colorN_["DAMAGE"].colorId = 38; // always last microphase
  colorN_["DAMAGE"].altName = "DAMAGE";
  colorN_["DAMAGE"].rgb.push_back(255);
  colorN_["DAMAGE"].rgb.push_back(0);
  colorN_["DAMAGE"].rgb.push_back(255);
  colorN_["DAMAGE"].gray = 255;

  map<string, elemColor>::iterator it = colorN_.begin();

  while (it != colorN_.end()) {
    // cout << "   " << it->first << " red = " << (it->second).rgb[0] << endl;
    // cout.flush();
    (it->second).rgbf.push_back((float)((it->second).rgb[0]) / 255.0);
    (it->second).rgbf.push_back((float)((it->second).rgb[1]) / 255.0);
    (it->second).rgbf.push_back((float)((it->second).rgb[2]) / 255.0);
    (it->second).grayf = (float)((it->second).gray) / 255.0;
    ++it;
  }

  // cout << "Made it past iterator" << endl;
  // cout.flush();

  /*
  colorN_[""].colorId = ;
  colorN_[""].altName = "";
  colorN_[""].rgb.push_back(255);
  colorN_[""].rgb.push_back();
  colorN_[""].rgb.push_back();
  colorN_[""].gray = ;
  */
}

//*@******************************************
//*@******************************************

void ChemicalSystem::checkChemSys(void) {
  int i, j, size, size_sec;

  cout << "" << endl;
  cout << "numMicroPhases_ " << numMicroPhases_ << endl;
  cout << "nnumICs_/numDCs_ " << numICs_ << " / " << numDCs_ << endl;
  cout << "numGEMPhases_ " << numGEMPhases_ << endl;
  cout << "numSolutionPhases_ " << numSolutionPhases_ << endl;
  cout << " " << endl;
  cout << "****************" << endl;
  cout << "vectors" << endl;
  size = microPhaseName_.size();
  cout << "microPhaseName_.size " << microPhaseName_.size() << endl;
  for (i = 0; i < size; i++) {
    cout << "   microPhaseName_[" << i << "] " << microPhaseName_[i] << endl;
  }
  cout << " " << endl;
  size = ICName_.size();
  cout << "ICName_.size " << ICName_.size() << endl;
  for (i = 0; i < size; i++) {
    cout << "   ICName_[" << i << "] " << ICName_[i] << endl;
  }
  size = DCName_.size();
  cout << "DCName_.size " << DCName_.size() << endl;
  for (i = 0; i < size; i++) {
    cout << "   DCName_[" << i << "] " << DCName_[i] << endl;
  }
  cout << " " << endl;
  size = GEMPhaseName_.size();
  cout << "GEMPhaseName_.size " << GEMPhaseName_.size() << endl;
  for (i = 0; i < size; i++) {
    cout << "   GEMPhaseName_[" << i << "] " << GEMPhaseName_[i] << endl;
  }
  size = microPhaseId_.size();
  cout << "microPhaseId_.size " << microPhaseId_.size() << endl;
  for (i = 0; i < size; i++) {
    cout << "   microPhaseId_[" << i << "] " << microPhaseId_[i] << endl;
  }

  cout << " " << endl;
  cout << "****************" << endl;
  cout << "maps" << endl;
  vector<int> second;
  // map<int,vector<int> > microPhaseMembers_;
  size = microPhaseMembers_.size();
  cout << "microPhaseMembers_.size " << size << endl;
  for (i = 0; i < size; i++) {
    cout << " " << endl;
    second = getMicroPhaseMembers(i);
    size_sec = second.size();
    cout << "   i/microPhaseName_/second.size " << i << " "
         << microPhaseName_[i] << " " << size_sec << endl;
    for (j = 0; j < size_sec; j++) {
      cout << "      second/microPhaseName_ " << second[j] << " "
           << microPhaseName_[second[j]] << endl;
    }
  }
  second.clear();

  cout << " " << endl;
  size = microPhaseDCMembers_.size();
  cout << "microPhaseDCMembers_.size " << size << endl;
  for (i = 0; i < size; i++) {
    second = getMicroPhaseDCMembers(i);
    size_sec = second.size();
    cout << "   second.size " << size_sec << endl;
    for (j = 0; j < size_sec; j++) {
      cout << "      second " << second[j] << endl;
    }
  }
  second.clear();

  cout << " " << endl;
  size = microPhaseToGEMPhase_.size();
  cout << "microPhaseToGEMPhase_.size " << size << endl;
  for (i = 0; i < size; i++) {
    second = getMicroPhaseToGEMPhase(i);
    size_sec = second.size();
    cout << "   second.size " << size_sec << endl;
    for (j = 0; j < size_sec; j++) {
      cout << "      second " << second[j] << endl;
    }
  }

  cout << " " << endl;
  cout << "****************" << endl;
  cout << "maps" << endl;

  cout << "numMicroPhases_ " << numMicroPhases_ << endl;
  size = microPhaseIdLookup_.size();
  if (size == numMicroPhases_) {
    cout << "microPhaseIdLookup_.size() OK! " << size << endl;
  } else {
    cout << "error -> microPhaseIdLookup_.size() /= numMicroPhases_ : " << size
         << " / " << numMicroPhases_ << endl;
    cout << "STOP";
    exit(1);
  }
  for (i = 0; i < size; i++) {
    cout << "   " << microPhaseName_[i] << " "
         << getMicroPhaseIdLookup(microPhaseName_[i]) << endl;
  }

  cout << " " << endl;
  cout << "numICs_" << numICs_ << endl;
  // map<string,int> ICIdLookup_;
  size = ICIdLookup_.size();
  if (size == numICs_) {
    cout << "ICIdLookup_.size() OK! " << size << endl;
  } else {
    cout << "error -> ICIdLookup_.size() /= numICs_ : " << size << " / "
         << numICs_ << endl;
    cout << "STOP";
    exit(1);
  }
  for (i = 0; i < size; i++) {
    cout << "   " << ICName_[i] << " " << getICIdLookup(ICName_[i]) << endl;
  }

  cout << " " << endl;
  cout << "numDCs_" << numDCs_ << endl;
  // map<string,int> DCIdLookup_;
  size = DCIdLookup_.size();
  if (size == numDCs_) {
    cout << "DCIdLookup_.size() OK! " << size << endl;
  } else {
    cout << "error -> DCIdLookup_.size() /= numDCs_ : " << size << " / "
         << numDCs_ << endl;
    cout << "STOP";
    exit(1);
  }
  for (i = 0; i < size; i++) {
    cout << "   " << DCName_[i] << " " << getDCIdLookup(DCName_[i]) << endl;
  }

  cout << " " << endl;
  cout << "numGEMPhases_ " << numGEMPhases_ << endl;
  // map<string,int> GEMPhaseIdLookup_
  size = GEMPhaseIdLookup_.size();
  if (size == numGEMPhases_) {
    cout << "GEMPhaseIdLookup_.size() OK! " << size << endl;
  } else {
    cout << "error -> GEMPhaseIdLookup_.size() /= numGEMPhases_ : " << size
         << " / " << numGEMPhases_ << endl;
    cout << "STOP";
    exit(1);
  }
  for (i = 0; i < size; i++) {
    cout << "   " << GEMPhaseName_[i] << " "
         << getGEMPhaseIdLookup(GEMPhaseName_[i]) << endl;
  }
}

void ChemicalSystem::writeSatElectrolyteGasConditions(void) {
  int DCId;
  double DCconc;  // mol/kgw units
  double DCmoles; // mol units
  cout << endl << "ChemicalSystem::writeSatElectrolyteGasConditions" << endl;
  if (isSaturated_) {
    cout << "   - saturated : 1 " << endl;
  } else {
    cout << "   - saturated : 0 " << endl;
  }

  bool initial = false;
  bool fixed = false;

  if (initialSolutionComposition_.size() > 0)
    initial = true;
  if (fixedSolutionComposition_.size() > 0)
    fixed = true;

  if (initial || fixed) {
    if (initial && fixed) {
      cout << endl << "error : check the chemistry.json file!" << endl;
      cout << endl
           << "if \"electrolyte_conditions\" is present, \"condition\" "
              "can be only one of the next combinations:"
           << endl;
      cout << "\"initial\"" << endl;
      cout << "\"fixed\"" << endl;
      cout << "\"attack\"" << endl;
      cout << "\"initial\" & \"attack\"" << endl;
      cout << "\"fixed\" & \"attack\"" << endl;
      cout << endl << "stop" << endl;
      exit(0);
    } else if (initial) {
      cout << "   - hydration starts with a given electrolyte composition :"
           << endl;
      cout << "      condition - initial" << endl;
      map<int, double>::iterator it = initialSolutionComposition_.begin();
      while (it != initialSolutionComposition_.end()) {
        DCId = it->first;
        if (DCId != waterDCId_) {
          DCconc = it->second;
          cout << "        DCId = " << DCId << "   DCName = " << DCName_[DCId]
               << "   DCconc = " << DCconc << " mol/kgw" << endl;
        }
        it++;
      }
    } else { // fixed
      cout << "   - during hydration the electrolyte has a fixed composition :"
           << endl;
      cout << "      condition - fixed" << endl;
      map<int, double>::iterator it = fixedSolutionComposition_.begin();
      while (it != fixedSolutionComposition_.end()) {
        DCId = it->first;
        if (DCId != waterDCId_) {
          DCconc = it->second;
          cout << "        DCId = " << DCId << "   DCName = " << DCName_[DCId]
               << "   DCconc = " << DCconc << " mol/kgw" << endl;
        }
        it++;
      }
    }
  } else {
    cout << "   - hydration starts without any dissolved ionic species (pure "
            "water)"
         << endl;
  }
  if (attackSolutionComposition_.size() > 0) {
    cout << "   - after hydration, simulation continues with an attack - "
            "solution composition being :"
         << endl;
    cout << "      condition - attack" << endl;
    map<int, double>::iterator it = attackSolutionComposition_.begin();
    while (it != attackSolutionComposition_.end()) {
      DCId = it->first;
      if (DCId != waterDCId_) {
        DCconc = it->second;
        cout << "        DCId = " << DCId << "   DCName = " << DCName_[DCId]
             << "   DCconc = " << DCconc << " mol/kgw" << endl;
      }
      it++;
    }
  }

  /*
  if (initialSolutionComposition_.size() > 0) {
    cout << "   - simulation starts with given solution composition :" << endl;
    cout << "      electrolyte - initial" << endl;
    map<int, double>::iterator it = initialSolutionComposition_.begin();
    while (it != initialSolutionComposition_.end()) {
      DCId = it->first;
      if (DCId != waterDCId_) {
        DCconc = it->second;
        cout << "        DCId = " << DCId << "   DCName = " << DCName_[DCId]
             << "   DCconc = " << DCconc << " mol/kgw" << endl;
      }
      it++;
    }
  } else {
    cout << "   - simulation without initial solution composition." << endl;
  }
  if (fixedSolutionComposition_.size() > 0) {
    cout << "   - simulation starts with given solution composition :" << endl;
    cout << "      electrolyte - fixed" << endl;
    map<int, double>::iterator it = fixedSolutionComposition_.begin();
    while (it != fixedSolutionComposition_.end()) {
      DCId = it->first;
      if (DCId != waterDCId_) {
        DCconc = it->second;
        cout << "        DCId = " << DCId << "   DCName = " << DCName_[DCId]
             << "   DCconc = " << DCconc << " mol/kgw" << endl;
      }
      it++;
    }
  } else {
    cout << "   - simulation without fixed solution composition." << endl;
  }
  if (attackSolutionComposition_.size() > 0) {
    cout << "   - simulation continues with attack, solution composition being
  :" << endl; cout << "      electrolyte - attack" << endl; map<int,
  double>::iterator it = attackSolutionComposition_.begin(); while (it !=
  attackSolutionComposition_.end()) { DCId = it->first; if (DCId != waterDCId_)
  { DCconc = it->second; cout << "        DCId = " << DCId << "   DCName = " <<
  DCName_[DCId]
             << "   DCconc = " << DCconc << " mol/kgw" << endl;
      }
      it++;
    }
  }
  */

  initial = false;
  fixed = false;
  if (initialGasComposition_.size() > 0)
    initial = true;
  if (fixedGasComposition_.size() > 0)
    fixed = true;

  if (initial || fixed) {
    if (initial && fixed) {
      cout << endl << "error : check the chemistry.json file!" << endl;
      cout << endl
           << "if \"gas_conditions\" is present, \"condition\" "
              "can be only one of the next combinations:"
           << endl;
      cout << "\"initial\"" << endl;
      cout << "\"fixed\"" << endl;
      cout << "\"attack\"" << endl;
      cout << "\"initial\" & \"attack\"" << endl;
      cout << "\"fixed\" & \"attack\"" << endl;
      cout << endl << "stop" << endl;
      exit(0);
    } else if (initial) {
      cout << "   - hydration starts in an atmosphere having an initial gas "
              "composition :"
           << endl;
      cout << "      gas - initial" << endl;
      map<int, double>::iterator it = initialGasComposition_.begin();
      while (it != initialGasComposition_.end()) {
        DCId = it->first;
        DCmoles = it->second;
        cout << "        DCId = " << DCId << "   DCName = " << DCName_[DCId]
             << "   DCmoles = " << DCmoles << " mol" << endl;
        it++;
      }
    } else { // fixed
      cout << "   - during hydration the atmosphere has a fixed gas "
              "composition :"
           << endl;
      cout << "      gas - fixed" << endl;
      map<int, double>::iterator it = fixedGasComposition_.begin();
      while (it != fixedGasComposition_.end()) {
        DCId = it->first;
        DCmoles = it->second;
        cout << "        DCId = " << DCId << "   DCName = " << DCName_[DCId]
             << "   DCmoles = " << DCmoles << " mol" << endl;
        it++;
      }
    }
  } else {
    cout << "   - hydration takes place in an inert atmosphere" << endl;
  }
  if (attackGasComposition_.size() > 0) {
    cout << "   - after hydration, simulation continues with gas attack - gas "
            "composition being :"
         << endl;
    cout << "      condition - attack" << endl;
    map<int, double>::iterator it = attackGasComposition_.begin();
    while (it != attackGasComposition_.end()) {
      DCId = it->first;
      if (DCId != waterDCId_) {
        DCconc = it->second;
        cout << "        DCId = " << DCId << "   DCName = " << DCName_[DCId]
             << "   DCconc = " << DCconc << " mol/kgw" << endl;
      }
      it++;
    }
  }

  if (attackSolutionComposition_.size() == 0 &&
      attackGasComposition_.size() == 0) {
    cout << "   - after hydration, no attack (electyrolyte/gas) takes place"
         << endl;
  }

  /*
  if (initialGasComposition_.size() > 0) {
    cout << "   - simulation with given gas composition :" << endl;
    cout << "      gas - initial" << endl;
    map<int, double>::iterator it = initialGasComposition_.begin();
    while (it != initialGasComposition_.end()) {
      DCId = it->first;
      DCmoles = it->second;
      cout << "        DCId = " << DCId << "   DCName = " << DCName_[DCId]
           << "   DCmoles = " << DCmoles << " mol" << endl;
      it++;
    }
  } else if (fixedGasComposition_.size() > 0) {
    cout << "   - simulation with given gas composition :" << endl;
    cout << "      gas - fixed" << endl;
    map<int, double>::iterator it = fixedGasComposition_.begin();
    while (it != fixedGasComposition_.end()) {
      DCId = it->first;
      DCmoles = it->second;
      cout << "        DCId = " << DCId << "   DCName = " << DCName_[DCId]
           << "   DCmoles = " << DCmoles << " mol" << endl;
      it++;
    }
  } else {
    cout << "   - simulation without given gas composition." << endl;
  }
  */
}

void ChemicalSystem::setElectrolyteComposition(const bool isFirst,
                                               bool doAttack, int cyc) {
  int DCId;
  double DCconc = 0.0; // mol/kgw units
  double waterMoles = DCMoles_[waterDCId_];
  double waterMass = 0.001 * waterMoles * waterMollarMass_; // in kg

  // cout << "  ----> ChemicalSystem::setElectrolyteComposition -
  // doAttack/waterMass : "
  //      << doAttack << " / " << waterMass << endl;

  if (doAttack) {
    cout << endl
         << "  ----> ChemicalSystem::setElectrolyteComposition - "
            "doAttack/waterMass : "
         << doAttack << " / " << waterMass << endl;

    if (attackSolutionComposition_.size() > 0) {
      map<int, double>::iterator it = attackSolutionComposition_.begin();
      while (it != attackSolutionComposition_.end()) {
        DCId = it->first;
        if (DCId != waterDCId_) {
          DCconc = it->second;
          DCMoles_[DCId] = DCconc * waterMass;
        }
        cout << "    ----> ChemicalSystem::setElectrolyteComposition attack - "
                "waterMass/DCId/DCconc/DCMoles_ : "
             << waterMass << " / " << DCId << " / " << DCconc << " / "
             << DCMoles_[DCId] << endl;
        it++;
      }
    }
  } else {
    if (isFirst && initialSolutionComposition_.size() > 0) {
      map<int, double>::iterator it = initialSolutionComposition_.begin();
      while (it != initialSolutionComposition_.end()) {
        DCId = it->first;
        if (DCId != waterDCId_) {
          DCconc = it->second;
          DCMoles_[DCId] = DCconc * waterMass;
        }
        it++;
      }
    }
    if (fixedSolutionComposition_.size() > 0) {
      map<int, double>::iterator it = fixedSolutionComposition_.begin();
      while (it != fixedSolutionComposition_.end()) {
        DCId = it->first;
        if (DCId != waterDCId_) {
          DCconc = it->second;
          DCMoles_[DCId] = DCconc * waterMass;
        }
        it++;
      }
    }
  }
  return;
}

void ChemicalSystem::setGasComposition(const bool isFirst, bool doAttack) {
  int DCId;
  double DCmoles; // mole units

  if (doAttack) {
    if (attackGasComposition_.size() > 0) {
      map<int, double>::iterator it = attackGasComposition_.begin();
      while (it != attackGasComposition_.end()) {
        DCId = it->first;
        DCmoles = it->second;
        DCMoles_[DCId] = DCmoles;
        it++;
      }
      // } else {
      //   cout << endl << "ChemicalSystem::setGasComposition() error for cyc =
      //   "
      //        << cyc << "  : " << endl;
      //   cout << "doAttack is true but atackGasComposition_.size() = "
      //        << atackGasComposition_.size() << endl;
      //   cout << endl << "stop program!" << endl;
      //   exit(0);
    }
  } else {
    if (isFirst && initialGasComposition_.size() > 0) {
      map<int, double>::iterator it = initialGasComposition_.begin();
      while (it != initialGasComposition_.end()) {
        DCId = it->first;
        DCmoles = it->second;
        DCMoles_[DCId] = DCmoles;
        it++;
      }
    } else if (fixedGasComposition_.size() > 0) {
      map<int, double>::iterator it = fixedGasComposition_.begin();
      while (it != fixedGasComposition_.end()) {
        DCId = it->first;
        DCmoles = it->second;
        DCMoles_[DCId] = DCmoles;
        it++;
      }
    }
  }
  return;
}

double ChemicalSystem::calculateCrystalStrain(int growPhId, double poreVolFrac,
                                              double Kp, double Ks) {
  double crystalStrain = -1; // 0.0;

  /*
  calculateState(false);
  */

  ///
  /// Crystallization pressure only exists if the solution is supersaturated,
  /// which means that \f$\beta > 1\f$.
  ///

  if (microPhaseSI_[growPhId] > 1.0) {

    ///
    /// Estimate the hydrostatic pressure in the pore solution as 1 atmosphere
    /// Unit of pressure for this calculation is MPa
    ///

    double pl = 0.101;

    ///
    /// The assumed largest radius of capillary pores offering entrance
    /// into gel porosity of C-S-H or other nanoporous component.
    /// 500 nm is chosen because it is half the size of a typical lattice
    /// site in THAMES (although the lattice resolution can be varied in
    /// the future if desired, so we may want to revisit this assumption).
    /// Unit of length for this calculation is millimeters.
    ///

    double r = 5.0e-4;

    ///
    /// Thickness of liquid film separating the crystallizing solid and
    /// the pore walls, assumed to be 1 nm.
    /// Unit of length for this calculation is millimeters.
    ///

    double delta = 1.0e-6;

    ///
    /// Crystal-liquid surface energy, assumed to be 100 mJ/m<sup>2</sup>.
    /// Unit of surface energy for this calculation is N/mm.
    ///

    double gamma = 1.0e-4; // N/mm

    ///
    /// Stress-free molar volume of the growing crystal
    /// Units of molar volume for this calcualtion is mm<sup>3</sup>/mol.
    ///
    /// @note This could be loaded up directly from the GEM CSD, rather than
    ///       hardwiring it into the code here.
    ///

    double Vc = 7.070e5;

    ///
    /// The ideal gas constant, with units of (N mm)/(mol K)
    ///

    double Rg = 8.314e3; // gas constant; N.mm/mol.K

    if (verbose_) {
      cout << "ChemicalSystem::calculateCrystalStrain SI for this phase is: "
           << microPhaseSI_[growPhId] << endl;
      cout.flush();
    }

    ///
    /// Calculate the crystal mean curvature in equilibrium with the
    /// solution with this saturation index (Thompson-Freundlich effect)
    ///

    double kcr = Rg * T_ * log(microPhaseSI_[growPhId]) / (Vc * gamma);

    ///
    /// If the portion of the crystal near the wall is a hemispherical cap,
    /// then the mean curvature is 2/r, where r is the radius of curvature
    /// of the crystal.  Therefore, the smallest pore within which the crystal
    /// can fit is delta larger than this.
    ///

    double rcr = (2.0 / kcr) + delta;

    ///
    /// Crystallization pressure associated with this pore size
    ///

    double pa = 2.0 * gamma * (1.0 / (rcr - delta) - 1.0 / (r - delta));

    ///
    /// Strain can be retrieved from the stress via the effeictive elastic
    /// of the porous medium (poromechanics assumption)
    ///

    crystalStrain =
        (1.0 / (3.0 * Kp) - 1.0 / (3.0 * Ks)) * (poreVolFrac * pa + pl);

    if (verbose_) {
      cout << "ChemicalSystem::calculateCrystalStrain crystalStrain = "
           << crystalStrain << endl;
      cout.flush();
    }
  }

  return crystalStrain;
}
