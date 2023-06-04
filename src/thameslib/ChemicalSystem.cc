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

ChemicalSystem::ChemicalSystem (Solution *Solut,
                                const string &GEMfilename,
                                const string &GEMdbrname,
                                const string &interfaceFileName,
                                const bool verbose,
                                const bool warning) : solut_(Solut)
{
    unsigned int i,j;
    double *amat;
    string exmsg;
    long int gemflag = 0;

    double *icmolarmass,*dcmolarmass;
    char *cc;
    unsigned int ii,jj;
    int k;
    bool found = false;
    timesGEMFailed_ = 0;
    maxGEMFails_ = 3;
    
    ///  The constructor initializes all the members to default values,
    ///  then launches the initial thermodynamic calculation, and sets
    ///  up the correspondences between GEM CSD phases and microstructure
    ///  phases.
    ///
    ///  All members are initialized to default values first, and all
    ///  vectors and maps are cleared.  The thermodynamic variables
    ///  are set to be consistent with neutral water at STP
    ///
 
    verbose_ = verbose;
    jsonFormat_ = false;
    warning_ = warning;
    numMicroPhases_ = numGEMPhases_ = numSolutionPhases_ = 0;
    numMicroImpurities_ = 4;
    microPhaseName_.clear();
    stressPhaseName_.clear();
    porousPhaseName_.clear();
    microPhaseId_.clear();
    isKinetic_.clear();
    DCName_.clear();
    GEMPhaseName_.clear();
    numICs_ = numDCs_ = numGEMPhases_ = 0;
    microPhaseMembers_.clear();
    microPhaseMemberVolumeFraction_.clear();
    microPhaseDCMembers_.clear();
    randomGrowth_.clear();
    growthTemplate_.clear();
    affinity_.clear();
    microPhasePorosity_.clear();
    poreSizeDistribution_.clear();
    k2o_.clear();
    na2o_.clear();
    mgo_.clear();
    so3_.clear();
    color_.clear();
    GEMPhaseStoich_.clear();
    GEMPhaseDCMembers_.clear();
    microPhaseDCPorosities_.clear();
    microPhaseIdLookup_.clear();
    ICIdLookup_.clear();
    DCIdLookup_.clear();
    GEMPhaseIdLookup_.clear();
    microPhaseToGEMPhase_.clear();
    DCStoich_.clear();
    ICClassCode_.clear();
    DCClassCode_.clear();
    GEMPhaseClassCode_.clear();

    Eh_ = 0.0;
    T_ = 298.0;             // Default temperature [K]
    P_ = 101325.0;          // Default pressure in [Pa]
    Vs_ = Ms_ = 1.0;
    Gs_ = Ms_ = 0.0;
    nodeStatus_ = NEED_GEM_AIA;
    nodeHandle_ = iterDone_ = 0;
    sulfateAttackTime_ = 1.0e10;
    leachTime_ = 1.0e10; 
    ICName_.clear();
    DCName_.clear();
    GEMPhaseName_.clear();
    ICIdLookup_.clear();
    DCIdLookup_.clear();
    GEMPhaseIdLookup_.clear();
    microPhaseVolume_.clear();
    microPhaseMass_.clear();
    microPhasePorosity_.clear();
    microPhaseMassDissolved_.clear();
    initialSolutionComposition_.clear();
    gasSolidRatio_ = 0.0;
    gasComposition_.clear();

    SI_.clear();
      
    node_ = new TNode();
   
    ///
    /// Initialize the thermodynamic system for both hydrates and solution 
    /// in order to initialize GEMPhaseVolume_ 
    ///

    char *cGEMfilename = (char*)GEMfilename.c_str();
    char *cGEMdbrname = (char*)GEMdbrname.c_str();
    if (verbose_) {
        cout << "ChemicalSystem::Going into GEM_init (1) to read CSD file "
             << cGEMfilename << endl;
    }

    /// Find out if the input data are in json format or in key-value format
    
    try {
        string json_dch = "";
        string json_ipm = "";
        string json_dbr = "";
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
    }
    catch (FileException fex) {
        throw fex;
    }

    if (gemflag == 1) {
        exmsg = "Bad return from GEM_init: " + GEMfilename + " missing or corrupt";
        throw GEMException("ChemicalSystem","ChemicalSystem",exmsg);
    }
    if (gemflag == -1) {
        exmsg = "Bad return from GEM_init: internal memory allocation error";
        throw GEMException("ChemicalSystem","ChemicalSystem",exmsg);
    }
    
    /// 
    /// Determine the number of possible ICs, DCs, and phases from the
    /// GEM CSD input that was read by GEM-IPM during initialization
    ///
    
    numICs_ = (unsigned int)((node_->pCSD())->nIC);
    numDCs_ = (unsigned int)((node_->pCSD())->nDC);
    numGEMPhases_ = (unsigned int)((node_->pCSD())->nPH);
    numSolutionPhases_ = (unsigned int)((node_->pCSD())->nPS);
    
    ///
    /// Knowing the dimensions, allocate the memory for all the arrays that
    /// must be created to store thermodynamic calculation results and communicate
    /// them to the microstructure
    ///

    try {
      exmsg = "ICMoles_";
      ICMoles_ = new double [numICs_];
      exmsg = "ICResiduals_";
      ICResiduals_ = new double [numICs_];
      exmsg = "ICChemicalPotential_";
      ICChemicalPotential_ = new double [numICs_];
      exmsg = "DCMoles_";
      DCMoles_ = new double [numDCs_];
      exmsg = "DCActivityCoeff_";
      DCActivityCoeff_ = new double [numDCs_];
      exmsg = "GEMPhaseMoles_";
      GEMPhaseMoles_ = new double [numGEMPhases_];
      exmsg = "solutPhaseMoles_";
      solutPhaseMoles_ = new double [numGEMPhases_];
      exmsg = "prevGEMPhaseMoles_";
      prevGEMPhaseMoles_ = new double [numGEMPhases_];
      exmsg = "GEMPhaseVolume_";
      GEMPhaseVolume_ = new double [numGEMPhases_];
      exmsg = "solutPhaseVolume_";
      solutPhaseVolume_ = new double [numGEMPhases_];
      exmsg = "prevGEMPhaseVolume_";
      prevGEMPhaseVolume_ = new double [numGEMPhases_];
      exmsg = "GEMPhaseMass_";
      GEMPhaseMass_ = new double [numGEMPhases_];
      exmsg = "solutPhaseMass_";
      solutPhaseMass_ = new double [numGEMPhases_];
      exmsg = "prevGEMPhaseMass_";
      prevGEMPhaseMass_ = new double [numGEMPhases_];
      exmsg = "surfaceArea_";
      surfaceArea_ = new double [numGEMPhases_];
      exmsg = "carrier_";
      carrier_ = new double [numSolutionPhases_];
      exmsg = "DCUpperLimit_";
      DCUpperLimit_ = new double [numDCs_];
      exmsg = "DCLowerLimit_";
      DCLowerLimit_ = new double [numDCs_];
      exmsg = "pGEMPhaseStoich_";
      pGEMPhaseStoich_ = new double [numGEMPhases_ * numICs_];
      exmsg = "pSolidStoich_";
      pSolidStoich_ = new double [numICs_];
      exmsg = "pSolutPhaseStoich_";
      pSolutPhaseStoich_ = new double [numGEMPhases_ * numICs_];
      exmsg = "pSolutSolidStoich_";
      pSolutSolidStoich_ = new double [numICs_];
    }
    catch (bad_alloc& ba) {
      cout << endl << "Bad_alloc Exception Thrown:" << endl;
      cout << "    Details:" << endl;
      cout << "    Offending function ChemicalSystem::ChemicalSystem" << endl;
      cout << "    Error in allocating memory for array " << exmsg << endl;
      cerr << endl << "Bad_alloc Exception Thrown:" << endl;
      cerr << "    Details:" << endl;
      cerr << "    Offending function ChemicalSystem::ChemicalSystem" << endl;
      cerr << "    Error in allocating memory for array " << exmsg << endl;
      exit(0);
    }
    
    ///
    /// Attempt to run GEM with auto initial approximation (AIA)
    ///
    /// This starts the thermodynamic calculation and returns the results, including
    /// the ionic strength, pH, IC chemical potentials, DC moles, phase moles, phase
    /// volumes, and other results of the calculation.  All of these parameters
    /// are loaded into the THAMES vectors that keep track of these things, since,
    /// they were passed to the GEM calculation by reference.
    /// 
    /// The argument is false if we wamt to use activity coefficients and speciation
    /// from a previous GEM_run, but is true if we want to use the activity coefficients
    /// and speciation stored in a DBR memory structure read from a DBR file
    ///
    /// Possible return values for nodeStatus_:
    ///    0 (NO_GEM_SOLVER): No GEM recalculation needed for node
    ///    1 (NEED_GEM_AIA) : Need GEM calc with LPP (auto initial approx, AIA)
    ///    2 (OK_GEM_AIA)   : OK after GEM calc with LPP AIA
    ///    3 (BAD_GEM_AIA)  : Not fully trusworthy result after calc with LPP AIA
    ///    4 (ERR_GEM_AIA)  : Failure (no result) in GEM calc with LPP AIA
    ///    5 (NEED_GEM_SIA) : Need GEM calc with no-LPP (smart initial approx, SIA)
    ///    6 (OK_GEM_SIA)   : OK after GEM calc with SIA
    ///    7 (BAD_GEM_SIA)  : Not fully trusworthy result after calc with SIA
    ///    8 (ERR_GEM_SIA)  : Failure (no result) in GEM calc with SIA
    ///    9 (T_ERROR_GEM ) : Terminal error (e.g., memory corruption).  Need restart
    ///

    (node_->pCNode())->NodeStatusCH = NEED_GEM_AIA;
    if (verbose_) {
        cout << "ChemicalSystem::Constructor: Entering GEM_run (1) with node status = "
             << nodeStatus_ << endl;
        cout.flush();
    }
    nodeStatus_ = node_->GEM_run(true);
    if (verbose_) {
      cout << "Done! nodeStatus is " << nodeStatus_ << endl;
      cout.flush();
    }
    if (!(nodeStatus_ == OK_GEM_AIA || nodeStatus_ == OK_GEM_SIA)) {
        bool dothrow = false;
        cerr << "ERROR: Call to GEM_run in ChemicalSystem constructor had an issue..." << endl;
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
                exmsg =  " Need GEM calc with smart initial approx (SIA)";
                dothrow = true;
                break;
            case BAD_GEM_SIA:
                exmsg = " Untrustworthy result with smart initial approx (SIA)";
                cout << exmsg << endl;
                dothrow = true;
                break;
            case ERR_GEM_SIA:
                exmsg =  " Failed result with smart initial approx (SIA)";
                cout << exmsg << endl;
                dothrow = true;
                break;
            case T_ERROR_GEM:
                exmsg = " Terminal GEM error; need restart";
                cout << exmsg << endl;
                dothrow = true;
                break;
            case NO_GEM_SOLVER:
                exmsg =  " No GEM recalculation needed for node";
                cout << exmsg << endl;
                dothrow = false;
                break;
        }
        if (dothrow) {
            throw GEMException("ChemicalSystem","Constructor",exmsg);
        }
    }

    if (verbose_) {
        cout << "ChemicalSystem::Constructor: Entering GEM_restore_MT (1) ... " << endl;
        cout.flush();
    }

    ///
    /// Next call passes (copies) the GEMS3K input data from the DBR structure.
    /// This is useful after the GEM_init and GEM_run() calls to initialize the
    /// arrays that keep the chemical data for all the nodes (one node in most cases)
    ///
    /// This function returns nothing and appears unable of throwing exceptions
    /// @todo Check carefully whether this function can throw an exception
    ///
    

    node_->GEM_restore_MT(nodeHandle_,nodeStatus_,T_,P_,Vs_,
                          Ms_,&ICMoles_[0],&DCUpperLimit_[0],
                          &DCLowerLimit_[0],&surfaceArea_[0]);

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
    /// This function returns nothing and appears unable of throwing exceptions
    /// @todo Check carefully whether this function can throw an exception
    ///

    node_->GEM_to_MT(nodeHandle_,nodeStatus_,iterDone_,Vs_,
        Ms_,Gs_,Hs_,ionicStrength_,pH_,pe_,Eh_,&ICResiduals_[0],
        &ICChemicalPotential_[0],&DCMoles_[0],&DCActivityCoeff_[0],
        &GEMPhaseMoles_[0],&GEMPhaseVolume_[0],&GEMPhaseMass_[0],
        &pGEMPhaseStoich_[0],&carrier_[0],&surfaceArea_[0],
        &pSolidStoich_[0]);
  
    /// The results of the thermodynamic calculation are now known, and
    /// the constructor can cast them into appropriate units and set up
    /// the data structure to make correspondences between GEM and microstructure
    ///
    /// Convert all IC and DC molar masses from kg/mol to g/mol
    ///
   
    ICMolarMass_.resize(numICs_,0.0);
    icmolarmass = (node_->pCSD())->ICmm;
    for (i = 0; i < numICs_; i++) {
      // Convert to g per mole
      ICMolarMass_[i] = (1000.0 * (double)(*icmolarmass));
      icmolarmass++;
    }
    DCMolarMass_.resize(numDCs_,0.0);
    dcmolarmass = (node_->pCSD())->DCmm;
    for (i = 0; i < numDCs_; i++) {
      // Convert to g per mole
      DCMolarMass_[i] = (1000.0 * (double)(*dcmolarmass));
      dcmolarmass++;
    }
  
    int maxICnameLength = node_->getMaxICnameLength(); 
    int maxDCnameLength = node_->getMaxDCnameLength(); 
    int maxPHnameLength = node_->getMaxPHnameLength(); 

    string string1;
    for (i = 0; i < numICs_; i++) {
        string1 = node_->xCH_to_IC_name(i);
        if( string1.length() >= maxICnameLength ) string1.resize(maxICnameLength);
        if (verbose_) {
            cout << "IC number " << i << " is " << string1 << endl;
        }
        ICName_.push_back(string1);
        ICIdLookup_.insert(make_pair(string1,i));
    }
    for (i = 0; i < numDCs_; i++) {
        string1 = node_->xCH_to_DC_name(i) ; 
        if( string1.length() >= maxDCnameLength ) string1.resize(maxDCnameLength);
        if (verbose_) {
            cout << "DC number " << i << " is " << string1 << endl;
        }
        DCName_.push_back(string1);
        DCIdLookup_.insert(make_pair(string1,i));
    }
    for (i = 0; i < numGEMPhases_; i++) {
        string1 = node_->xCH_to_Ph_name(i);
        if( string1.length() >= maxPHnameLength ) string1.resize(maxPHnameLength);
        if (verbose_) {
            cout << "PH number " << i << " is " << string1 << endl;
        }
        GEMPhaseName_.push_back(string1);
        GEMPhaseIdLookup_.insert(make_pair(string1,i));
    }

    ///
    /// Set up the stoichiometry matrix for dependent components (DCs) in terms
    /// of independent components (ICs).  This is the GEM CSD A matrix
    ///
 
    vector<double> scplaceholder;
    scplaceholder.clear();
    scplaceholder.resize(numICs_,0);
    DCStoich_.resize(numDCs_,scplaceholder);
    amat = (node_->pCSD())->A;
    for (i = 0; i < numDCs_; i++) {
      for (j = 0; j < numICs_; j++) {
        DCStoich_[i][j] = (double)(*amat);
        amat++;
      } 
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
  
    for (long int i = 0; i < numGEMPhases_; i++) {
        GEMPhaseMass_[i] = (double)(node_->Ph_Mass(i) * 1000.0); // in g, not kg
        prevGEMPhaseMass_[i] = (double)(node_->Ph_Mass(i) * 1000.0); // in g, not kg
    }

    setGEMPhaseVolume();
    setGEMPhaseMolarMass();

    ///
    /// Set up the class codes for ICs, DCs, and phases, based on the type of
    /// component they are.  Refer to the documentation for these individual members
    /// for more detailed information about allowable values of the class codes
    ///

    ICClassCode_.resize(numICs_,' ');
    cc = (node_->pCSD())->ccIC;
    for (i = 0; i < numICs_; i++) {
      ICClassCode_[i] = *cc;
      cc++;
    }

    DCClassCode_.resize(numDCs_,' ');
    cc = (node_->pCSD())->ccDC;
    for (i = 0; i < numDCs_; i++) {
      DCClassCode_[i] = *cc;
      cc++;
    }

    GEMPhaseClassCode_.resize(numGEMPhases_,' ');
    cc = (node_->pCSD())->ccPH;
    for (i = 0; i < numGEMPhases_; i++) {
      GEMPhaseClassCode_[i] = *cc;
      cc++;
    }

    ///
    /// Begin parsing the chemistry input XML file
    ///

    string msg;
    string xmlext = ".xml";
    size_t foundxml = interfaceFileName.find(xmlext);
    try {
      if (foundxml != string::npos) {
        parseDoc(interfaceFileName);
          
        microPhaseVolume_.resize(numMicroPhases_,0.0);
        microPhaseMass_.resize(numMicroPhases_,0.0);
        microPhasePorosity_.resize(numMicroPhases_,0.0);
        if (verbose_) {
            cout << " Setting microPhaseMass size to "
                 << numMicroPhases_ << endl;
            cout.flush();
        }
        microPhaseMassDissolved_.resize(numMicroPhases_,0.0);
      } else {
        msg = "Not an XML file";
        throw FileException("ChemicalSystem","ChemicalSystem",
                            interfaceFileName,msg);
      }
    }
    catch (FileException e) {
      throw e;
    }
    
    ///
    /// Set up the main map that correlates microstructure phases with GEM CSD phases
    ///

    microInitVolume_ = 0.0;
    for (unsigned int i = 0; i < numMicroPhases_; i++) {
      microPhaseToGEMPhase_.insert(make_pair((int)i,microPhaseMembers_[i]));
    }

    ///
    /// Set up the vector of saturation indices for each GEM phase.
    /// This determines the driving force for growth and is also used in calculations
    /// of the crystallization pressure during external sulfate attack

    setSI();
    vector<double> SIforsystem = getSI();

    ///
    /// Set up all the information for the composition of the aqueous solution
    ///

    vector<double> solutionICMoles = getSolution();

    solut_->setICMoles(solutionICMoles);
    try {
        solut_->calculateState(true);
    }
    catch (GEMException gex) {
        gex.printException();
        cout << endl;
    }
    vector<double> solutionSI = solut_->getSI();
}

bool ChemicalSystem::isInputFormatJSON (const char *masterFileName)
{
    ifstream in(masterFileName);
    if (!in) {
        throw FileException("Solution","isInputFormatJSON",
                            masterFileName,"Could not open");
    }

    string filetypeflag;
    if (in.peek() != EOF) {
        in >> filetypeflag;
    } else {
        throw FileException("Solution","isInputFormatJSON",
                            masterFileName,"Bad or corrupt format");
    }
    in.close();
    if (filetypeflag == "-j") {
        return true;
    }
    return false;
}

void ChemicalSystem::getJSONFiles (const char *masterFileName,
                                   string &dchName,
                                   string &ipmName,
                                   string &dbrName)
{
    ifstream in(masterFileName);
    if (!in) {
        throw FileException("Solution","getJSONFiles",
                            masterFileName,"Could not open");
    }

    string filetypeflag;
    if (in.peek() != EOF) {
        in >> filetypeflag;
    } else {
        throw FileException("Solution","getJSONFiles",
                            masterFileName,"Bad or corrupt format");
    }
    if (in.peek() != EOF) {
        in >> dchName;
    } else {
        throw FileException("Solution","getJSONFiles",
                            masterFileName,"Bad or corrupt format");
    }
    if (in.peek() != EOF) {
        in >> ipmName;
    } else {
        throw FileException("Solution","getJSONFiles",
                            masterFileName,"Bad or corrupt format");
    }
    if (in.peek() != EOF) {
        in >> dbrName;
    } else {
        throw FileException("Solution","getJSONFiles",
                            masterFileName,"Bad or corrupt format");
    }

    /// Assuming that all the file names were read correctly, strip their quote marks
    
    dchName.erase(remove(dchName.begin(),dchName.end(),'"'),dchName.end());
    ipmName.erase(remove(ipmName.begin(),ipmName.end(),'"'),ipmName.end());
    dbrName.erase(remove(dbrName.begin(),dbrName.end(),'"'),dbrName.end());
}

vector<double> ChemicalSystem::getSolution (void)
{
  ///
  /// Get IC moles for solution, which fully characterizes the
  /// composition of the aqueous solution
  ///
 
  double waterMass = DCMoles_[getDCId("H2O@")] 
                 * DCMolarMass_[getDCId("H2O@")];

  if (verbose_) {
      cout << "water mass is: "
           << waterMass << endl;
  }
  vector<double> tempICMoles;
  tempICMoles.clear();
  tempICMoles.resize(numICs_,0.0);
  for (unsigned int i = 0; i < numDCs_; i++) {
    char cc = getDCClassCode(i);
    if (cc == 'S' || cc == 'T') {
      double moles = node_->Get_cDC(i) * waterMass * 1.0e-3;
      for(int j = 0; j < (numICs_ - 1); j++) {
        tempICMoles[j] += moles * DCStoich_[i][j];
      }
    }
  }

  // Treat H2O separately

  double waterMoles = DCMoles_[getDCId("H2O@")];
  for (int j = 0; j < numICs_; j++) {
    if (ICName_[j] == "H") tempICMoles[j] += waterMoles * 2;
    if (ICName_[j] == "O") tempICMoles[j] += waterMoles;
  }

  return tempICMoles;  

}

void ChemicalSystem::parseDoc (const string &docName)
{
    // Need to open the docName and scan it somehow for
    // phase names and id numbers

    string msg;
    PhaseData phaseData;
    xmlDocPtr doc;
    xmlChar *key;
    xmlNodePtr cur;
    cout.flush();
    doc = xmlParseFile(docName.c_str());

    /// Check if the xml file is valid
    /// @note This block requires the schema file to be local

    string rxcsd = "chemistry.xsd";
    if(!is_xml_valid(doc,rxcsd.c_str())) {
        cout << "Chemistry xml is NOT valid" <<endl;
        cout.flush();
    } else if (verbose_) {
        cout << "Chemistry xml IS valid" << endl;
        cout.flush();
    }

    if (doc == NULL ) {
        msg = "XML file not parsed successfully";
        throw FileException("ChemicalSystem","parseDoc",
                            docName,msg);
    }

    cur = xmlDocGetRootElement(doc);

    if (cur == NULL) {
        msg = "XML file is empty";
        xmlFreeDoc(doc);
        throw FileException("ChemicalSystem","parseDoc",
                            docName,msg);
    }

    /// XML file is valid and non-empty
    /// Do a first scan through to get the phase names and ids
    ///

    cur = cur->xmlChildrenNode;
    map<string,int> phaseids;
    phaseids.clear();
    if (verbose_) {
        cout << "Parsing phase names" << endl;
        cout.flush();
    }
    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"phase"))) {
            parsePhaseNames (doc, cur, phaseids);
        }
        cur = cur->next;
    }

    ///
    /// Now back up to the beginning and scan properly
    ///
    /// The interface file contains information about each microstructure
    /// phase that is defined, including the list of GEM CSD phases that
    /// are to be associated with that phase, the phase's internal porosity,
    /// dissolved impurities, and visualization properties.
    ///

    if (verbose_) {
        cout << "Back to the top of the xml file to parse properly" << endl;
        cout.flush();
    }
    cur = xmlDocGetRootElement(doc);
    cur = cur->xmlChildrenNode;
    int testnumEntries = 0;
    int satstate = 1;
    isSaturated_ = true;
    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"numentries"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(testnumEntries,st);
            xmlFree(key);
        } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"saturated"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(satstate,st);
            if (satstate == 0) isSaturated_ = false;
        } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"electrolyte"))) {
            parseSolutionComp(doc, cur);
        } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"gas"))) {
            parseGasComp(doc, cur);
        } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"phase"))) {
            try {
                parsePhase(doc, cur, testnumEntries, phaseids, phaseData);
            }
            catch (FileException fex) {
               fex.printException();
               cout << endl;
            }
            catch (GEMException gex) {
               gex.printException();
               cout << endl;
            }
        }
        cur = cur->next;
    }

    xmlFreeDoc(doc);
    return;
}

void ChemicalSystem::parseSolutionComp (xmlDocPtr doc,
                                        xmlNodePtr cur)
{
    // Clear the associative map to initialize it
    
    initialSolutionComposition_.clear();

    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"ICcomp"))) {
            parseICInSolution(doc, cur);
        }
        cur = cur->next;
    }

    return;
}

void ChemicalSystem::parseGasComp (xmlDocPtr doc,
                                        xmlNodePtr cur)
{
    // Clear the associative map to initialize it
    
    xmlChar *key;

    double gassolidratio = 0.0;
    gasComposition_.clear();

    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"gassolidratio"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            from_string(gassolidratio,(char *)key);
            setGasSolidRatio(gassolidratio);
        } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"ICcomp"))) {
            parseICInGas(doc, cur);
        }
        cur = cur->next;
    }

    return;
}

void ChemicalSystem::parseICInSolution (xmlDocPtr doc,
                                        xmlNodePtr cur)
{
    xmlChar *key;
    int ICId = -1;
    string ICName;
    double ICConc = -1.0;

    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"name"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            from_string(ICName,(char *)key);
            ICId = getICId(ICName);
            if (ICConc > 0.0 && ICId > 0) {
                initialSolutionComposition_.insert(make_pair(ICId,ICConc));
                ICId = -1;
                ICConc = -1.0;
                ICName = "Unknown";
            }

            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"conc"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(ICConc,st);
            xmlFree(key);
            if (ICConc > 0.0 && ICId > 0) {
                initialSolutionComposition_.insert(make_pair(ICId,ICConc));
                ICId = -1;
                ICConc = -1.0;
                ICName = "Unknown";
            }
        }

        cur = cur->next;
    }

    return;
}

void ChemicalSystem::parseICInGas (xmlDocPtr doc,
                                   xmlNodePtr cur)
{
    xmlChar *key;
    int ICId = -1;
    string ICName;
    double ICConc = -1.0;

    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"name"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            from_string(ICName,(char *)key);
            ICId = getICId(ICName);
            if (ICConc > 0.0 && ICId > 0) {
                gasComposition_.insert(make_pair(ICId,ICConc));
                ICId = -1;
                ICConc = -1.0;
                ICName = "Unknown";
            }

            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"conc"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(ICConc,st);
            xmlFree(key);
            if (ICConc > 0.0 && ICId > 0) {
                gasComposition_.insert(make_pair(ICId,ICConc));
                ICId = -1;
                ICConc = -1.0;
                ICName = "Unknown";
            }
        }

        cur = cur->next;
    }

    return;
}

void ChemicalSystem::parsePhaseNames (xmlDocPtr doc,
                                      xmlNodePtr cur,
                                      map<string,int> &phaseids)
{

    xmlChar *key;
    cur = cur->xmlChildrenNode;
    int pid;
    string pname;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"id"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(pid,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"thamesname"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            pname = st;
            xmlFree(key);
        }
        cur = cur->next;
    }
    if (verbose_) {
        cout << "    Phase name = " << pname << endl;
        cout.flush();
    }
    phaseids.insert(make_pair(pname,pid));
}


void ChemicalSystem::parsePhase (xmlDocPtr doc,
                                 xmlNodePtr cur,
                                 int numEntries,
                                 map<string,int> phaseids,
                                 PhaseData &phaseData)
{
    xmlChar *key;

    string poreSizeFileName;

    phaseData.growthTemplate.clear();
    phaseData.affinity.clear();

    /// @note The affinity vector is always the same length, one entry for every
    /// microstructure phase, and the default value is zero.  Therefore
    /// the chemistry file does not need to include zero affinity values
    
    phaseData.affinity.resize(numEntries,0);
    phaseData.GEMPhaseId.clear();
    phaseData.DCId.clear();
    phaseData.GEMPhaseName.clear();
    phaseData.microPhaseDCPorosities.clear();
    phaseData.DCName.clear();
    phaseData.stressCalc = 0;
    phaseData.weak = 0;

    cur = cur->xmlChildrenNode;

    /// @note This parsing ignores the kinetic data portion for each
    /// phase.  The kinetic data parsing is handled by the KineticModel class
   
    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"id"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phaseData.id,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"thamesname"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            phaseData.thamesName = st;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"gemphase_data"))) {
            parseGEMPhaseData(doc, cur, phaseData);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"poresizefilename"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(poreSizeFileName,st);
            try {
                parsePoreSizeDistribution(poreSizeFileName, phaseData);
            }
            catch (FileException fex) {
               fex.printException();
               cout << endl;
            }
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"stresscalc"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phaseData.stressCalc,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"weak"))) {
            // Weak means the phase can be damaged by stress
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phaseData.weak,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"display_data"))) {
            parseDisplayData(doc, cur, phaseData);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"impurity_data"))) {
            parseImpurityData(doc, cur, phaseData);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"interface_data"))) {
            parseInterfaceData(doc, cur, phaseids, phaseData);
        }
        cur = cur->next;
    }

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
    for (int i = 0; i < phaseData.microPhaseDCPorosities.size() && !done; ++i) {
        if (phaseData.microPhaseDCPorosities[i] > 0.0 && phaseData.microPhaseDCPorosities[i] < 1.0) {
            porousPhaseName_.push_back(phaseData.thamesName);
            porousPhaseId_.push_back(phaseData.id);
            done = true;
        }
    }

    // The following is a lookup table that will return
    // the ordered list of DC porosities for a given microstructure phase id
  
    microPhaseDCPorosities_.insert(make_pair(phaseData.id,phaseData.microPhaseDCPorosities));

    phaseData.microPhaseDCPorosities.clear();

    microPhaseName_.push_back(phaseData.thamesName);
    microPhaseId_.push_back(phaseData.id);
    microPhaseIdLookup_.insert(make_pair(phaseData.thamesName,phaseData.id));
    randomGrowth_.push_back(phaseData.randomGrowth);
    affinity_.push_back(phaseData.affinity);

    /// Growth template is based on positive affinities only
    
    growthTemplate_.push_back(calcGrowthtemplate(phaseData.affinity));

    poreSizeDistribution_.push_back(phaseData.poreSizeDist);
    if (verbose_) {
        cout << "Pushed pore size distribution data for phase " << phaseData.thamesName << endl;
        cout << "This phase distribution has " << phaseData.poreSizeDist.size() << " entries" << endl;
        cout << "Have now registered " << poreSizeDistribution_.size() << " PSDs" << endl;
        cout.flush();
    }
    phaseData.poreSizeDist.clear();
    grayscale_.push_back(phaseData.gray);
    color_.push_back(phaseData.colors);
    k2o_.push_back(phaseData.k2o);
    na2o_.push_back(phaseData.na2o);
    mgo_.push_back(phaseData.mgo);
    so3_.push_back(phaseData.so3);
    microPhaseMembers_.insert(make_pair(phaseData.id,phaseData.GEMPhaseId));
    microPhaseDCMembers_.insert(make_pair(phaseData.id,phaseData.DCId));

    numMicroPhases_++;

    if (verbose_) {
        cout << "Parsed phase " << phaseData.thamesName << endl;
        cout.flush();
    }

    return;
}

void ChemicalSystem::parsePoreSizeDistribution(string poreSizeFileName,
                                               PhaseData &phaseData)
{
    if (verbose_) {
        cout << "Reading Pore Size Distribution:" << endl;
        cout.flush();
    }

    ifstream in(poreSizeFileName);
    if (!in) {
        throw FileException("ChemicalSystem","parsePoreSizeDistribution",
                            poreSizeFileName,"Could not open");
    }

    // Read the header line
    string headerline;
    getline(in,headerline);

    struct PoreSizeVolume datarow;
    double diam,vfrac;

    phaseData.poreSizeDist.clear();

    // Now read the data row by row
    double sum = 0.0;
    while (!in.eof()) {
        in >> datarow.diam >> datarow.volfrac;
        sum += datarow.volfrac;
        datarow.volume = 0.0;
        phaseData.poreSizeDist.push_back(datarow);
        if (verbose_) {
            cout << "---> " << datarow.diam << "," << datarow.volfrac << endl;
            cout.flush();
        }
        in.peek();
    }

    sum -= datarow.volfrac;
    phaseData.poreSizeDist.erase(phaseData.poreSizeDist.end() - 1);

    if (verbose_) {
        cout << "<---- sum = " << sum << endl;
        cout.flush();
    }
    in.close();

    // Normalize the pore size distribution in case it is not already
    if (sum > 0.0) {
        for (int i = 0; i < phaseData.poreSizeDist.size(); ++i) {
            phaseData.poreSizeDist[i].volfrac *= (1.0/sum);
        }
    }
   
    return;
}

void ChemicalSystem::parseGEMPhaseData (xmlDocPtr doc,
                                        xmlNodePtr cur,
                                        PhaseData &phaseData)
{
    xmlChar *key;
    cur = cur->xmlChildrenNode;

    int GEMPhaseId = 0;
    int dcid = 0;
    string mypstr;
    phaseData.GEMPhaseDCMembers.clear();
    bool scrapeWaterDCs = false;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name,(const xmlChar *)"gemphasename"))) {
            key = xmlNodeListGetString(doc,cur->xmlChildrenNode,1);
            phaseData.GEMPhaseName.push_back((char *)key);
            mypstr = (char *)key;
            if (mypstr == WaterGEMName) {
                scrapeWaterDCs = true;
            } else {
                scrapeWaterDCs = false;
            }
            cout << "GEM Phase name = " << mypstr
                 << ", scrapeWaterDCs = " << scrapeWaterDCs << endl;
            cout.flush();
            // Assign the global microstructure phase name associated with CSH
            if (mypstr == CSHGEMName) {
                CSHMicroName = phaseData.thamesName;
            }
            if (mypstr == MonocarbGEMName) {
                MonocarbMicroName = phaseData.thamesName;
            }
            if (mypstr == HydrotalcGEMName) {
                HydrotalcMicroName = phaseData.thamesName;
            }
            GEMPhaseId = getGEMPhaseId((char *)key);
            phaseData.GEMPhaseId.push_back(GEMPhaseId);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name,(const xmlChar *)"gemdc"))) {
            parseGEMPhaseDCData(doc, cur, phaseData);
        }
        cur = cur->next;
    } 
    GEMPhaseDCMembers_.insert(make_pair(GEMPhaseId,phaseData.GEMPhaseDCMembers));   

    return;

}

void ChemicalSystem::parseGEMPhaseDCData (xmlDocPtr doc,
                                       xmlNodePtr cur,
                                       PhaseData &phaseData)
{
    xmlChar *key;
    string mydcstr;
    cur = cur->xmlChildrenNode;

    int dcid = 0;
    double porosity = 0.0;
    bool scrapeWaterDCs = false;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"gemdcname"))) {
            key = xmlNodeListGetString(doc,cur->xmlChildrenNode,1);
            phaseData.DCName.push_back((char *)key);
            mydcstr = (char *)key;
            cout << "GEM DC name = " << mydcstr
                 << ", scrapeWaterDCs = " << scrapeWaterDCs << endl;
            cout.flush();
            if (mydcstr == AFTDCName) {
                AFTMicroName = phaseData.thamesName;
            }
            if (mydcstr == MonosulfDCName) {
                MonosulfMicroName = phaseData.thamesName;
            }
            if (!scrapeWaterDCs || mydcstr == WaterDCName) {
                // Only for aqueous solution, we keep only the water DC,
                // not all the dissolved components
                dcid = getDCId((char *)key);
                phaseData.DCId.push_back(dcid);
                phaseData.GEMPhaseDCMembers.push_back(dcid);
                cout << "GEM DC id = " << dcid
                     << ", scrapeWaterDCs = " << scrapeWaterDCs << endl;
            }
            // Make certain that there will be a porosity associated
            // with this DC
            if (phaseData.microPhaseDCPorosities.size() < phaseData.DCName.size()) {
                phaseData.microPhaseDCPorosities.push_back(0.0);
            }
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"gemdcporosity"))) {
            key = xmlNodeListGetString(doc,cur->xmlChildrenNode,1);
            string st((char *)key);
            from_string(porosity,st);
            // Make certain that there will be a porosity associated
            // with this DC
            if (phaseData.microPhaseDCPorosities.size() < phaseData.DCName.size()) {
                phaseData.microPhaseDCPorosities.push_back(porosity);
            } else {
                phaseData.microPhaseDCPorosities.at(phaseData.microPhaseDCPorosities.size()-1) = porosity;
            }
            xmlFree(key);
        }

        cur = cur->next;
    }
}

void ChemicalSystem::parseDisplayData (xmlDocPtr doc,
                                       xmlNodePtr cur,
                                       PhaseData &phaseData)
{

    xmlChar *key;
    cur = cur->xmlChildrenNode;
    double red,green,blue;

    red = green = blue = 0.0;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"red"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(red,st);
            // if (verbose_) cout << "        red = " << red << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"green"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(green,st);
            // if (verbose_) cout << "        green = " << green << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"blue"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(blue,st);
            // if (verbose_) cout << "        blue = " << blue << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"gray"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phaseData.gray,st);
            // if (verbose_) cout << "        gray = " << phaseData.gray << endl;
            xmlFree(key);
        }
        cur = cur->next;

    }

    phaseData.colors.clear();
    phaseData.colors.push_back(red);
    phaseData.colors.push_back(green);
    phaseData.colors.push_back(blue);

    return;
}

void ChemicalSystem::parseImpurityData (xmlDocPtr doc,
                                        xmlNodePtr cur,
                                        PhaseData &phaseData)
{
    xmlChar *key;
    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"k2ocoeff"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phaseData.k2o,st);
            // if (verbose_) cout << "        k2ocoeff = " << phaseData.k2o << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"na2ocoeff"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phaseData.na2o,st);
            // if (verbose_) cout << "        na2ocoeff = " << phaseData.na2o << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"mgocoeff"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phaseData.mgo,st);
            // if (verbose_) cout << "        mgocoeff = " << phaseData.mgo << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"so3coeff"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phaseData.so3,st);
            // if (verbose_) cout << "        so3coeff = " << phaseData.so3 << endl;
            xmlFree(key);
        }
        cur = cur->next;

    }
    return;
}

void ChemicalSystem::parseInterfaceData (xmlDocPtr doc,
                                         xmlNodePtr cur,
                                         map<string,int> &phaseids,
                                         PhaseData &phaseData)
{

    xmlChar *key;
    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"randomgrowth"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phaseData.randomGrowth,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"affinity"))) {
            parseAffinityData(doc,cur,phaseids,phaseData);
        }
        cur = cur->next;

    }
    return;
}

void ChemicalSystem::parseAffinityData (xmlDocPtr doc,
                                        xmlNodePtr cur,
                                        map<string,int> &phaseids,
                                        PhaseData &phaseData)
{
    xmlChar *key;
    cur = cur->xmlChildrenNode;
    int testaftyid, testaftyval;
    map<string,int>::iterator it = phaseids.begin();

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"affinityphase"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            map<string,int>::iterator it = phaseids.find(st);
            if (it != phaseids.end()) {
                testaftyid = it->second;
            }
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"affinityvalue"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(testaftyval,st);
            xmlFree(key);
        }
        cur = cur->next;
    }

    // Must be a valid phase name to be added to the list
    if (it != phaseids.end()) {
        phaseData.affinity[testaftyid] = ((int)testaftyval);
    }
    return;
}

ChemicalSystem::ChemicalSystem (const ChemicalSystem &obj)
{

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
    randomGrowth_ = obj.getRandomGrowth();
    ICMoles_ = obj.getICMoles();
    DCMoles_ = obj.getDCMoles();
    ICMolarMass_ = obj.getICMolarMass();
    DCMolarMass_ = obj.getDCMolarMass();
    growthTemplate_ = obj.getGrowthTemplate();
    affinity_ = obj.getAffinity();
    microPhaseMembers_ = obj.getMicroPhaseMembers();
    microPhaseMemberVolumeFraction_ = obj.getMicroPhaseMemberVolumeFraction();
    microPhaseDCMembers_ = obj.getMicroPhaseDCMembers();
    microPhasePorosity_ = obj.getMicroPhasePorosity();
    poreSizeDistribution_ = obj.getPoreSizeDistribution();
    k2o_ = obj.getK2o();
    na2o_ = obj.getNa2o();
    mgo_ = obj.getMgo();
    so3_ = obj.getSo3();
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
    surfaceArea_ = obj.getSurfaceArea();
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
    microInitVolume_ = obj.getMicroInitVolume();
    microPhaseMass_ = obj.getMicroPhaseMass();
    microPhaseMassDissolved_ = obj.getMicroPhaseMassDissolved();
    microVoidVolume_ = obj.getMicroVoidVolume();
    verbose_ = obj.getVerbose();
    warning_ = obj.getWarning();
}

ChemicalSystem::~ChemicalSystem (void)
{
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
    randomGrowth_.clear();
    DCStoich_.clear();
    growthTemplate_.clear();
    affinity_.clear();
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
    grayscale_.clear();
    color_.clear();
    ICClassCode_.clear();
    DCClassCode_.clear();
    GEMPhaseClassCode_.clear();
    GEMPhaseStoich_.clear();

    ///
    /// Free up the dynamically allocated memory
    ///

    delete[]DCLowerLimit_;
    delete[]DCUpperLimit_;
    delete[]surfaceArea_;
    delete[]prevGEMPhaseMass_;
    delete[]prevGEMPhaseVolume_;
    delete[]prevGEMPhaseMoles_;
    delete[]GEMPhaseMass_;
    delete[]GEMPhaseVolume_;
    delete[]carrier_;
    delete[]GEMPhaseMoles_;
    delete[]DCActivityCoeff_;
    delete[]DCMoles_;
    delete[]ICChemicalPotential_;
    delete[]ICResiduals_;
    delete[]ICMoles_;
    delete[]pGEMPhaseStoich_;

    delete node_;
}

void ChemicalSystem::getPGEMPhaseStoich (void)
{
    double *arout = new double[numICs_];
    for (long int i = 0; i < numGEMPhases_; i++) {
        arout = node_->Ph_BC(i,arout);
        for (unsigned int j = 0; j < numICs_; j++) {
            pGEMPhaseStoich_[(i * numICs_) + j] = arout[j];
        }
    }
    delete[] arout;
}

void ChemicalSystem::getGEMPhaseStoich (void)
{
    double minval = 0.0;
    GEMPhaseStoich_.clear();
    vector<double> vplace;
    vplace.clear();
    vplace.resize(numICs_,0.0);
    GEMPhaseStoich_.resize(numGEMPhases_,vplace);
    int indexval,oval;
    for (unsigned int i = 0; i < numGEMPhases_; i++) {
        if (GEMPhaseName_[i] == "aq_gen") {

            ///
            /// Normalize to one mole of oxygen
            ///

            oval = (i * numICs_) + getICId("O");
            if (pGEMPhaseStoich_[oval] > 0.0) {
                for (unsigned int j = 0; j < numICs_; j++) {
                    indexval = (i * numICs_) + j;
                    GEMPhaseStoich_[i][j] =
                        (pGEMPhaseStoich_[indexval]/pGEMPhaseStoich_[oval]);
                }
            }
        } else {
            for (unsigned int j = 0; j < numICs_; j++) {
                indexval = (i * numICs_) + j;
                GEMPhaseStoich_[i][j] =
                        (pGEMPhaseStoich_[indexval]/minval);
            }
        }
    }
}

void ChemicalSystem::writeDb (ostream &stream)
{
    unsigned int i;

    ///
    /// Make the header
    ///

    stream << "--------------------------------------------------------" << endl;
    stream << "CONTENTS OF PHASE DATABASE:" << endl;
    stream <<  endl;

    ///
    /// Format one line at a time
    ///

    for (i = 0; i < numMicroPhases_; i++) {
        writeMember(i,stream);
    }

    stream << endl;
    stream << "--------------------------------------------------------" << endl;
    stream << endl;
}

void ChemicalSystem::writeMember (const unsigned int i,
                                  ostream &stream)
{

    unsigned int idnum;
    if (i >= numMicroPhases_) {
        throw EOBException("ChemicalSystem","writeMember",
                           "microPhaseName_",numMicroPhases_,i);
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

void ChemicalSystem::writeChemSys (void)
{
    unsigned int j;

    ///
    /// First we will list details for the ICs
    ///

    string CSfilename("chemsys.report");
    ofstream out(CSfilename.c_str());
    out << "Report on the Material Database" << endl;
    out << "-------------------------------" << endl << endl;
    out << "List of Independent Components:" << endl << endl;
    for(unsigned int i = 0; i < numICs_; i++) {
        out << i << ")            Name: " << ICName_[i] << endl;
        out << "        classcode: " << ICClassCode_[i] << endl;
        out << "       molar mass: " << ICMolarMass_[i] << endl << endl;
    }

    out << "List of Dependent Components:" << endl << endl;
    for(unsigned int i = 0; i < numDCs_; i++) {
        out << i << ")            Name: " << DCName_[i] << endl;
        out << "        classcode: " << DCClassCode_[i] << endl;
        out << "       molar mass: " << DCMolarMass_[i] << endl << endl;
    }

    out << "List of Phases:" << endl << endl;
    for(unsigned int i = 0; i < numGEMPhases_; i++) {
        out << i << ")            Name: " << GEMPhaseName_[i] << endl;
        out << "        classcode: " << GEMPhaseClassCode_[i] << endl;
    }

    out << "List of Microstructure Phases:" << endl << endl;
    for(unsigned int i = 0; i < numMicroPhases_; i++) {
        out << i << ")       Name: " << microPhaseName_[i] << endl;
        out << "               id: " << microPhaseId_[i] << endl;
        out << "    random growth: " << randomGrowth_[i] << endl;
        for (j = 0; j < affinity_[i].size(); j++) {
            out << "        affinity to " << j
                << ": " << affinity_[i][j] << endl;
        }
        for (j = 0; j < growthTemplate_[i].size(); j++) {
            out << "        growthTemplate: "
                << growthTemplate_[i][j] << endl;
        }
        out << "         porosity: " << microPhasePorosity_[i] << endl;
        out << "              k2o: " << k2o_[i] << endl;
        out << "             na2o: " << na2o_[i] << endl;
        out << "              mgo: " << mgo_[i] << endl;
        out << "              so3: " << so3_[i] << endl;
    }

    return;
}

void ChemicalSystem::setMicroPhaseMass (const unsigned int idx,
                                        const double val)
{
    try {
        microPhaseMass_.at(idx) = val;
    }
    catch (out_of_range &oor) {
        EOBException ex("ChemicalSystem","setMicroPhaseMass",
                           "microPhaseMass_",microPhaseMass_.size(),idx);
        ex.printException();
        exit(1);
    }

    int DCId = 0;
    if (idx == ELECTROLYTEID) {
        DCId = getDCId("H2O@");
    } else if (idx != VOIDID) {
        DCId = getMicroPhaseDCMembers(idx,0);
    }
    if (idx != VOIDID) {
        double v0 = node_->DC_V0(DCId,P_,T_);
        double dcmm = getDCMolarMass(DCId);
        if (verbose_) {
            cout << "    " << microPhaseName_[idx] << ", DC " << getDCName(DCId) << ": v0 = "
                 << v0 << ", dcmm = " << dcmm
                 << ", so volume = ";
            cout.flush();
        }
        if (dcmm < 1.0e-9) {
            FloatException fex("ChemicalSystem","setMicroPhaseMass",
                               "Divide by zero (dcmm)");
            fex.printException();
            exit(1);
        }
        if (verbose_) {
            cout << (val*v0/dcmm) << endl;
            cout << "Setting volume of microphase " << idx << " to "
                 << (val*v0/dcmm) << " (VOIDID = " << VOIDID
                  << ")" << endl;
            cout.flush();
        }
        setMicroPhaseVolume(idx,(val*v0/dcmm));
    }

    return;
}

void ChemicalSystem::calcMicroPhasePorosity (const unsigned int idx)
{

    // No point in doing this if the phase is not present

    double porosity = 0.0;

    // Get the first GEM phase for this microstructure phase
   
    /// @todo Generalize calcMicroPhasePorosity so that we
    /// don't just check the first GEM Phase listed as part
    /// of this microstructure phase
 
    unsigned int gemphaseid = getMicroPhaseToGEMPhase (idx,0);

    // Find all the DC ids for this GEM phase

    vector<int> DClist = getGEMPhaseDCMembers(gemphaseid);

    // Get the porosities for each DC in the microstructure phase

    vector<double> DCporosities = getMicroPhaseDCPorosities(idx);
    
    if (verbose_) {
        cout << "ChemicalSystem::calcMicroPhasePorosity for "
             << getMicroPhaseName(idx) << endl;
        cout << "    This phase's GEM phase id = " << gemphaseid
             << " and volume = " << getMicroPhaseVolume(idx) << " m3" << endl;
        cout << "    The DC members for this GEM phase are:" << endl;
        for (int ii = 0; ii < DClist.size(); ++ii) {
            cout << "        " << DClist[ii] << " (" << getDCName(DClist[ii])
                 << ")" << endl;
            cout.flush();
        }
    }

    /// @todo Do we need to check that DClist and DCporosities are the same size?
    /// @todo Do we need to check that DClist and DCporosities are in the same order?
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
       
    if (DClist.size() == 1) {
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
        for (int i = 0; i < DClist.size(); ++i) {
            DCId = DClist[i];
            conc = getDCConcentration(DCId);
            porosity = DCporosities[i];
            vol = conc * getDCMolarVolume(DCId);
            weightedporosities += (vol * porosity);
            sumvol += vol;
            if (verbose_) {
                cout << "    " << getDCName(DCId) << " (" << DCId << ") concentration = " << conc << endl;
                cout << "    " << getDCName(DCId) << " (" << DCId << ") porosity = " << porosity << endl;
                cout << "    " << getDCName(DCId) << " (" << DCId << ") molarvolume = " << getDCMolarVolume(DCId) << endl;
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
            cout << "    " << getMicroPhaseName(idx) << " subvoxel porosity = " << porosity << endl;
            cout.flush();
        }
    }

    setMicroPhasePorosity(idx,porosity);

    return;
}

int ChemicalSystem::calculateState (double time,
                                    bool isFirst = false)
{
    int status = 0;
    string msg;
 
    // isFirst = true; 
      
    vector<double> oDCMoles;
    oDCMoles.clear();
    oDCMoles.resize(numDCs_,0.0);
  
    for (int i = 0; i < numDCs_; i++) {
      oDCMoles[i] = DCMoles_[i];
    }

    nodeStatus_ = NEED_GEM_AIA;

    if (verbose_) {
        cout << "    Before calculateState, "
             << "printing microPhaseVolumes" << endl;
        vector<double> microPhaseVolumes = getMicroPhaseVolume();
        vector<string> microPhaseNames = getMicroPhaseName();
        for (int i = 0; i < microPhaseVolumes.size(); ++i) {
            cout << "    Phase name " << microPhaseNames[i]
                 << ": volume = " << microPhaseVolumes[i] << endl;
            cout.flush();
        }

        cout << "    Going into ChemicalSystem::calculateState::GEM_from_MT (2)... "
             << endl;
    }

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


    node_->GEM_from_MT(nodeHandle_,nodeStatus_,T_,P_,Vs_,Ms_,
          ICMoles_,DCUpperLimit_,DCLowerLimit_,surfaceArea_,
          DCMoles_,DCActivityCoeff_);

    if (verbose_) {
        cout << "Done!" << endl;
        cout << "    Going into ChemicalSystem::calculateState::GEM_set_MT (2)... ";
        cout.flush();
    }

    /// For passing the current THAMES time and time step into the working instance
    /// of the DBR structure.
    ///
    /// This function returns nothing and appears to be incapable of throwing an exception
    ///
    ///
    /// @todo Check carefully if this function can throw an exception
    ///
    /// @note This function call seems to be optional.  It is not needed for
    /// the GEM data structures or calculations, but mostly for debugging purposes.
    ///
    /// @note This function used to be called GEM_set_MT in older versions of GEMS3K
    ///
    
    node_->GEM_from_MT_time(time,1.0);

    if (verbose_) {
        cout << "Done!" << endl
             << "    Going into "
             << "ChemicalSystem::calculateState::GEM_run() "
             << "(2), with isFirst " << isFirst << endl;
        cout.flush();
        writeICMoles();
    }

    /*
    writeDCMoles();
    */

    ///
    /// Attempt to run GEM with automatic initial approximation (AIA)
    ///
    /// This starts the thermodynamic calculation and returns the results, including
    /// the ionic strength, pH, IC chemical potentials, DC moles, phase moles, phase
    /// volumes, and other results of the calculation.  All of these parameters
    /// are loaded into the THAMES vectors that keep track of these things, since,
    /// they were passed to the GEM calculation by reference.
    /// 
    /// The argument is false if we wamt to use activity coefficients and speciation
    /// from a previous GEM_run, but is true if we want to use the activity coefficients
    /// and speciation stored in a DBR memory structure read from a DBR file
    ///
    /// Possible return values for nodeStatus_:
    ///    0 (NO_GEM_SOLVER): No GEM recalculation needed for node
    ///    1 (NEED_GEM_AIA) : Need GEM calc with LPP (auto initial approx, AIA)
    ///    2 (OK_GEM_AIA)   : OK after GEM calc with LPP AIA
    ///    3 (BAD_GEM_AIA)  : Not fully trusworthy result after calc with LPP AIA
    ///    4 (ERR_GEM_AIA)  : Failure (no result) in GEM calc with LPP AIA
    ///    5 (NEED_GEM_SIA) : Need GEM calc with no-LPP (smart initial approx, SIA)
    ///    6 (OK_GEM_SIA)   : OK after GEM calc with SIA
    ///    7 (BAD_GEM_SIA)  : Not fully trusworthy result after calc with SIA
    ///    8 (ERR_GEM_SIA)  : Failure (no result) in GEM calc with SIA
    ///    9 (T_ERROR_GEM ) : Terminal error (e.g., memory corruption).  Need restart
    ///

    checkICMoles();

    if (isFirst) {
        (node_->pCNode())->NodeStatusCH = NEED_GEM_AIA;
    } else {
        (node_->pCNode())->NodeStatusCH = NEED_GEM_SIA;
    }
    nodeStatus_ = node_->GEM_run(true);

    if (verbose_) {
        cout << "Done!  nodeStatus is " << nodeStatus_ << endl;
        cout.flush();
    }

    if (!(nodeStatus_ == OK_GEM_AIA || nodeStatus_ == OK_GEM_SIA)) {
        bool dothrow = false;
        cerr << "ERROR: Call to GEM_run in "
             << "ChemicalSystem::calculateState had an issue..." << endl;
        cerr << "       nodeStatus_ = " << nodeStatus_;
        switch (nodeStatus_) {
            case NEED_GEM_AIA:
                msg = " Need GEM calc with auto initial approx (AIA)";
                cerr << msg << endl;
                dothrow = false;
                break;
            case BAD_GEM_AIA:
                msg = " Untrustworthy result with auto initial approx (AIA)",
                cerr << msg << endl;
                dothrow = false;
                break;
            case ERR_GEM_AIA:
                msg = " Failed result with auto initial approx (AIA)";
                cerr << msg << ", GEMS failed "
                     << timesGEMFailed_ << " times" << endl;
                node_->GEM_print_ipm("IPM_dump.txt");
                timesGEMFailed_++;
                dothrow = (timesGEMFailed_ > maxGEMFails_) ? true : false;
                break;
            case NEED_GEM_SIA:
                msg =  " Need GEM calc with smart initial approx (SIA)";
                cerr << msg << endl;
                dothrow = false;
                break;
            case BAD_GEM_SIA:
                msg = " Untrustworthy result with smart initial approx (SIA)";
                cerr << msg << endl;
                dothrow = false;
                break;
            case ERR_GEM_SIA:
                msg =  " Failed result with smart initial approx (SIA)";
                cerr << msg << ", GEMS failed "
                     << timesGEMFailed_ << " times" << endl;
                node_->GEM_print_ipm("IPM_dump.txt");
                timesGEMFailed_++;
                dothrow = (timesGEMFailed_ > maxGEMFails_) ? true : false;
                break;
            case T_ERROR_GEM:
                msg = " Terminal GEM error; need restart";
                cerr << msg << endl;
                dothrow = true;
                break;
            case NO_GEM_SOLVER:
                msg =  " No GEM recalculation needed for node";
                cerr << msg << endl;
                dothrow = false;
                break;
        }
        if (dothrow) {
            throw GEMException("ChemicalSystem","calculateState",msg);
        }
    } else {
        timesGEMFailed_ = 0;
    }

    if (timesGEMFailed_ > 0) {
        if (verbose_) {
            cout << "Call to GEM_run has failed "
                 << timesGEMFailed_ << " consecutive times.  "
                 << "Attempt this step again" << endl;
        }
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

    checkICMoles();

    node_->GEM_to_MT(nodeHandle_,nodeStatus_,iterDone_,Vs_,
            Ms_,Gs_,Hs_,ionicStrength_,pH_,pe_,Eh_,&ICResiduals_[0],
            &ICChemicalPotential_[0],&DCMoles_[0],&DCActivityCoeff_[0],
            &solutPhaseMoles_[0],&solutPhaseVolume_[0],
            &solutPhaseMass_[0],&pSolutPhaseStoich_[0],
            &carrier_[0],&surfaceArea_[0],&pSolidStoich_[0]);

    if (verbose_) {
        cout << "Done!" << endl;
        cout << "after GEM_to_MT...Ms_ = " << Ms_ << endl;
        cout.flush();
    }

    /*
    writePhasemoles();
    */

    microVolume_ = 0.0;
    setPGEMPhaseStoich();
    setGEMPhaseStoich();
    setGEMPhaseMass();
    setGEMPhaseVolume();
    setGEMPhaseMolarMass();
  
    if (verbose_) {
        cout << "%%%%%%%%%% Printing GEM Masses and "
             << "Volumes in this Step %%%%%%%" << endl;
        for (long int myid = 0; myid < numGEMPhases_; myid++) {
            cout << "Mass and volume of GEM phase "
                 << node_->pCSD()->PHNL[myid] << " = "
                 << GEMPhaseMass_[myid] << " g and "
                 << GEMPhaseVolume_[myid] << " m3" << endl;
        }
        cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
             << "%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    }

    /// JWB 2023-04-14
    /// Account for subvoxel porosity within each microstructure phase here
    /// We used to do this in the Lattice class but it is cleaner here
   
    double phi;  // local variable to store subvoxel volume fraction of pores of a phase
                 // 0 <= phi <= 1
                 
    for (unsigned int i = 1; i < numMicroPhases_; i++) {
        if (verbose_) {
            cout << "Setting microPhase amounts for "
                 << i << " = " << microPhaseName_[i] << endl;
            cout.flush();
        }
        if (!isKinetic(i)) {
            calcMicroPhasePorosity(i);
            phi = getMicroPhasePorosity(i);
            microPhaseMass_[i] = microPhaseVolume_[i] = 0.0;
            for (unsigned int j = 0; j < microPhaseMembers_[i].size(); j++) {
                microPhaseMass_[i] +=
                     GEMPhaseMass_[microPhaseMembers_[i][j]];

                if (verbose_) {
                    cout << "    Is NOT a KINETIC phase: is composed of "
                         << GEMPhaseName_[microPhaseMembers_[i][j]]
                         << " having mass = "
                         << GEMPhaseMass_[microPhaseMembers_[i][j]]
                         << " and volume = "
                         << GEMPhaseVolume_[microPhaseMembers_[i][j]]
                         << " and porosity = " << phi
                         << endl;
                    cout.flush();
                }

                /// Here is where the subvoxel porosity is included
                /// @todo Need a more disciplined approach to making sure
                /// these are solid phases and not capillary porosity
              
                if (i == ELECTROLYTEID) {
                    microPhaseVolume_[i] += GEMPhaseVolume_[microPhaseMembers_[i][j]];
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
        } else {
            if (verbose_) {
                calcMicroPhasePorosity(i);
                phi = getMicroPhasePorosity(i);

                cout << "    IS a KINETIC phase: is composed of "
                     << GEMPhaseName_[microPhaseMembers_[i][0]]
                     << "  having mass = " << microPhaseMass_[i]
                     << "  and volume = " << microPhaseVolume_[i] << endl;
                cout.flush();
            }

            ///
            /// microPhaseMass and microPhaseVolume for kinetic phases are
            /// already set in KineticModel::calculateKineticStep
            ///

            microVolume_ += (microPhaseVolume_[i] / (1.0 - phi));
        }
    }

    if (isSaturated_) {   // System is saturated

        if (verbose_) {
            cout << "Use water to saturate the "
                 << "porosity." << endl;
        }
  
        double water_molarv, water_molesincr;
        if (microInitVolume_ > microVolume_) {
            double water_molarv, water_molesincr;
            for (int i = 0; i < numMicroPhases_; i++) {
                if (microPhaseName_[i] == "H2O") {
                    water_molarv = node_->DC_V0(getMicroPhaseMembers(i,0), P_, T_);
                    water_molesincr = (microInitVolume_ - microVolume_) / water_molarv;
                    if (verbose_) {
                        cout << "water_molarv = "
                             << water_molarv << endl;
                        cout << "volume increase of water is: "
                             << (microInitVolume_ - microVolume_) << endl;
                        cout << "water_molesincr = "
                             << water_molesincr << endl;
                    }
                }
            }
            for (int i = 0; i < numICs_; i++) {
                if (ICName_[i] == "H") ICMoles_[i] += water_molesincr * 2.0;
                if (ICName_[i] == "O") ICMoles_[i] += water_molesincr;
            }
        }

    }

    if (verbose_) {
        cout << "GEM volume change = "
             << 100.0 * (microVolume_ - microInitVolume_)/(microInitVolume_)
             << " %" << endl;
        cout.flush();
    }

    setGEMPhaseStoich();
   
    ///
    /// Calculate driving force for growth or dissolution
    ///

    // double *soluticmoles;
    // soluticmoles = solut_->getICmoles();
    // for (int i = 0; i < numDCs_; i++) {
    //     char cc;
    //     cc = getDCClassCode(i);
    //     if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M') {
    //         double dissolveDCMoles = oDCMoles[i] - DCMoles_[i];
    //         if (dissolveDCMoles > 0.0) {
    //             for (int j = 0; j < (numICs_ - 1); j++) {
    //                 soluticmoles[j] += dissolveDCMoles * DCStoich_[i][j];
    //             }
    //         }
    //     }
    // }
    // for (int i = 0; i < numICs_; i++) {
    //     solut_->setICMoles(i, soluticmoles[i]);
    // }
    // try {
    //     solut_->calculateState(true);
    // }
    // catch (GEMException gex) {
    //     gex.printException();
    //     cout << endl;
    // }

    ///
    /// Update solution
    ///

    vector<double> solutICMoles = getSolution();
    if (verbose_) cout << "Now update solution IC moles...";
    for (int i = 0; i < numICs_; i++) {
        solut_->setICMoles(i, solutICMoles[i]);
    } 
    if (verbose_) cout << "Done." << endl;

    try {
        // Solution state was calculated for the first
        // time in the constructor.
        solut_->calculateState(false);
    }
    catch (GEMException gex) {
        throw gex;
    }

    if (verbose_) {
        cout << "Leaving ChemicalSystem::calculateState now" << endl;
        cout.flush();
    }

    return timesGEMFailed_;
}
