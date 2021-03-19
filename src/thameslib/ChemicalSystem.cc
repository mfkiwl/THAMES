/**
@file ChemicalSystem.cc
@brief Method definitions for the ChemicalSystem base class
*/

#include "ChemicalSystem.h"

ChemicalSystem::ChemicalSystem (Solution *Solut,
                                const string &GEMfilename,
                                const string &GEMdbrname,
                                const string &Interfacefilename,
                                const bool verbose,
                                const bool warning,
                                const bool debug) : solut_(Solut)
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
    timesGEMfailed_ = 0;
    maxGEMfails_ = 3;
    
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
    debug_ = debug;
    jsonformat_ = false;
    warning_ = warning;
    micphasenum_ = phasenum_ = solutionphasenum_ = 0;
    micimpuritynum_ = 4;
    micphasename_.clear();
    stressphasename_.clear();
    porousphasename_.clear();
    micid_.clear();
    DCname_.clear();
    phasename_.clear();
    ICnum_ = DCnum_ = phasenum_ = 0;
    micphasemembers_.clear();
    micphasemembersvolfrac_.clear();
    micDCmembers_.clear();
    randomgrowth_.clear();
    growthtemplate_.clear();
    affinity_.clear();
    porosity_.clear();
    k2o_.clear();
    na2o_.clear();
    mgo_.clear();
    so3_.clear();
    color_.clear();
    vphasestoich_.clear();
    micidlookup_.clear();
    ICidlookup_.clear();
    DCidlookup_.clear();
    phaseidlookup_.clear();
    mic2phase_.clear();
    mic2DC_.clear();
    kineticphase_.clear();
    thermophase_.clear();
    DCstoich_.clear();
    ICclasscode_.clear();
    DCclasscode_.clear();
    phaseclasscode_.clear();
    Eh_ = 0.0;
    T_ = 298.0;             // Default temperature [K]
    P_ = 101325.0;          // Default pressure in [Pa]
    Vs_ = Ms_ = 1.0;
    Gs_ = Ms_ = 0.0;
    nodestatus_ = nodehandle_ = iterdone_ = 0;
    sattack_time_ = 1.0e10;
    leach_time_ = 1.0e10; 
    ICname_.clear();
    DCname_.clear();
    phasename_.clear();
    ICidlookup_.clear();
    DCidlookup_.clear();
    phaseidlookup_.clear();
    micphasevolfrac_.clear();
    micphasevolume_.clear();
    micphasemass_.clear();
    micphasemassdissolved_.clear();
    SI_.clear();
      
    node_ = new TNode();
   
    ///
    /// Initialize the thermodynamic system for both hydrates and solution 
    /// in order to initialize phasevolume_ 
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
        jsonformat_ = isInputFormatJSON(cGEMfilename);

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

        // if (jsonformat_) {
        //     if (verbose) {
        //         cout << "Detected JSON input file format for ChemicalSystem GEM data files" << endl;
        //         cout.flush();
        //     }
        //     getJSONFiles(cGEMfilename,json_dch,json_ipm,json_dbr);
        //     gemflag = node_->GEM_init(json_dch,json_ipm,json_dbr);
        // } else {
        //     if (verbose) {
        //         cout << "Detected key-value input file format for ChemicalSystem GEM data files" << endl;
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
  
    ICnum_ = (unsigned int)((node_->pCSD())->nIC);
    DCnum_ = (unsigned int)((node_->pCSD())->nDC);
    phasenum_ = (unsigned int)((node_->pCSD())->nPH);
    solutionphasenum_ = (unsigned int)((node_->pCSD())->nPS);
    
    ///
    /// Knowing the dimensions, allocate the memory for all the arrays that
    /// must be created to store thermodynamic calculation results and communicate
    /// them to the microstructure
    ///

    try {
      exmsg = "ICmoles_";
      ICmoles_ = new double [ICnum_];
      exmsg = "ICresiduals_";
      ICresiduals_ = new double [ICnum_];
      exmsg = "ICchempot_";
      ICchempot_ = new double [ICnum_];
      exmsg = "DCmoles_";
      DCmoles_ = new double [DCnum_];
      exmsg = "DCactivitycoeff_";
      DCactivitycoeff_ = new double [DCnum_];
      exmsg = "phasemoles_";
      phasemoles_ = new double [phasenum_];
      exmsg = "solutphasemoles_";
      solutphasemoles_ = new double [phasenum_];
      exmsg = "ophasemoles_";
      ophasemoles_ = new double [phasenum_];
      exmsg = "phasevolume_";
      phasevolume_ = new double [phasenum_];
      exmsg = "solutphasevolume_";
      solutphasevolume_ = new double [phasenum_];
      exmsg = "ophasevolume_";
      ophasevolume_ = new double [phasenum_];
      exmsg = "phasemass_";
      phasemass_ = new double [phasenum_];
      exmsg = "solutphasemass_";
      solutphasemass_ = new double [phasenum_];
      exmsg = "ophasemass_";
      ophasemass_ = new double [phasenum_];
      exmsg = "surfacearea_";
      surfacearea_ = new double [phasenum_];
      exmsg = "carrier_";
      carrier_ = new double [solutionphasenum_];
      exmsg = "DCupperlimit_";
      DCupperlimit_ = new double [DCnum_];
      exmsg = "DClowerlimit_";
      DClowerlimit_ = new double [DCnum_];
      exmsg = "phasestoich_";
      phasestoich_ = new double [phasenum_ * ICnum_];
      exmsg = "solidstoich_";
      solidstoich_ = new double [ICnum_];
      exmsg = "solutphasestoich_";
      solutphasestoich_ = new double [phasenum_ * ICnum_];
      exmsg = "solutsolidstoich_";
      solutsolidstoich_ = new double [ICnum_];
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
    /// Possible return values for nodestatus_:
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

    (node_->pCNode())->NodeStatusCH = NEED_GEM_SIA;
    if (verbose_) {
        cout << "ChemicalSystem::Constructor: Entering GEM_run (1) with node status = "
             << nodestatus_ << endl;
        cout.flush();
    }
    nodestatus_ = node_->GEM_run(true);
    if (verbose_) {
      cout << "Done! nodestatus is " << nodestatus_ << endl;
      cout.flush();
    }
    if (!(nodestatus_ == OK_GEM_AIA || nodestatus_ == OK_GEM_SIA)) {
        bool dothrow = false;
        cerr << "ERROR: Call to GEM_run in ChemicalSystem constructor had an issue..." << endl;
        cerr << "       nodestatus_ = " << nodestatus_;
        switch (nodestatus_) {
            case NEED_GEM_AIA:
                exmsg = " Need GEM calc with auto initial approx (AIA)";
                cerr << exmsg << endl;
                dothrow = true;
                break;
            case BAD_GEM_AIA:
                exmsg = " Untrustworthy result with auto initial approx (AIA)",
                cerr << exmsg << endl;
                dothrow = true;
                break;
            case ERR_GEM_AIA:
                exmsg = " Failed result with auto initial approx (AIA)";
                cerr << exmsg << endl;
                dothrow = true;
                break;
            case NEED_GEM_SIA:
                exmsg =  " Need GEM calc with smart initial approx (SIA)";
                dothrow = true;
                break;
            case BAD_GEM_SIA:
                exmsg = " Untrustworthy result with smart initial approx (SIA)";
                cerr << exmsg << endl;
                dothrow = true;
                break;
            case ERR_GEM_SIA:
                exmsg =  " Failed result with smart initial approx (SIA)";
                cerr << exmsg << endl;
                dothrow = true;
                break;
            case T_ERROR_GEM:
                exmsg = " Terminal GEM error; need restart";
                cerr << exmsg << endl;
                dothrow = true;
                break;
            case NO_GEM_SOLVER:
                exmsg =  " No GEM recalculation needed for node";
                cerr << exmsg << endl;
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
    

    node_->GEM_restore_MT(nodehandle_,nodestatus_,T_,P_,Vs_,Ms_,&ICmoles_[0],
        &DCupperlimit_[0],&DClowerlimit_[0],&surfacearea_[0]);

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

    node_->GEM_to_MT(nodehandle_,nodestatus_,iterdone_,Vs_,
        Ms_,Gs_,Hs_,ionicstrength_,pH_,pe_,Eh_,&ICresiduals_[0],
        &ICchempot_[0],&DCmoles_[0],&DCactivitycoeff_[0],&phasemoles_[0],
        &phasevolume_[0],&phasemass_[0],&phasestoich_[0],
        &carrier_[0],&surfacearea_[0],&solidstoich_[0]);
  
    /// The results of the thermodynamic calculation are now known, and
    /// the constructor can cast them into appropriate units and set up
    /// the data structure to make correspondences between GEM and microstructure
    ///
    /// Convert all IC and DC molar masses from kg/mol to g/mol
    ///
   
    ICmolarmass_.resize(ICnum_,0.0);
    icmolarmass = (node_->pCSD())->ICmm;
    for (i = 0; i < ICnum_; i++) {
      ICmolarmass_[i] = (1000.0 * (double)(*icmolarmass));  // converts to g per mol
      icmolarmass++;
    }
    DCmolarmass_.resize(DCnum_,0.0);
    dcmolarmass = (node_->pCSD())->DCmm;
    for (i = 0; i < DCnum_; i++) {
      DCmolarmass_[i] = (1000.0 * (double)(*dcmolarmass));  // converts to g per mol
      dcmolarmass++;
    }
  
    string string1;
    for (i = 0; i < ICnum_; i++) {
      string1.assign(node_->xCH_to_IC_name(i));
      ICname_.push_back(string1);
      ICidlookup_.insert(make_pair(string1,i));
    }
    for (i = 0; i < DCnum_; i++) {
      // if (verbose_) {
      //    cout << "DC id " << i << " name is " << node_->xCH_to_DC_name(i) << endl;
      //    cout.flush();
      //}
      string1.assign(node_->xCH_to_DC_name(i));
  
      DCname_.push_back(string1);
      DCidlookup_.insert(make_pair(string1,i));
    }
    for (i = 0; i < phasenum_; i++) {
      string1.assign(node_->xCH_to_Ph_name(i));
      phasename_.push_back(string1);
      phaseidlookup_.insert(make_pair(string1,i));
    }

    //if (verbose_) {
    //    cout << "To initialize phasevolume_ and phasemass_, set DCupperlimit to be normal: " 
    //         << endl;
    //    for (int i = 0; i < DCnum_; i++) {
    //      cout << DCname_[i] << ": " << DCupperlimit_[i] << endl;
    //    }
    //}

    ///
    /// Set up the stoichiometry matrix for dependent components (DCs) in terms
    /// of independent components (ICs).  This is the GEM CSD A matrix
    ///
 
    vector<double> scplaceholder;
    scplaceholder.clear();
    scplaceholder.resize(ICnum_,0);
    DCstoich_.resize(DCnum_,scplaceholder);
    amat = (node_->pCSD())->A;
    for (i = 0; i < DCnum_; i++) {
      for (j = 0; j < ICnum_; j++) {
        DCstoich_[i][j] = (double)(*amat);
        amat++;
      } 
    }
  
    ///
    /// Set up the stoichiometry and molar masses of the GEM CSD phases
    ///

    setPhasestoich();
    setVphasestoich();

    ///
    /// Normally we can call setPhasemass() to read in the GEM phase masses
    /// from the GEM CSD, but this first time we have to do it manually
    /// because the units are converted from kg to g and this will mess
    /// up the assignment of the ophasemass_ values, which would still be
    /// in kg if we called the setPhasemass() function.
    ///
  
    for (long int i = 0; i < phasenum_; i++) {
        phasemass_[i] = (double)(node_->Ph_Mass(i) * 1000.0); // in g, not kg
        ophasemass_[i] = (double)(node_->Ph_Mass(i) * 1000.0); // in g, not kg
    }

    setPhasevolume();
    setPhasemolarmass();

    ///
    /// Set up the class codes for ICs, DCs, and phases, based on the type of
    /// component they are.  Refer to the documentation for these individual members
    /// for more detailed information about allowable values of the class codes
    ///

    ICclasscode_.resize(ICnum_,' ');
    cc = (node_->pCSD())->ccIC;
    for (i = 0; i < ICnum_; i++) {
      ICclasscode_[i] = *cc;
      cc++;
    }

    DCclasscode_.resize(DCnum_,' ');
    cc = (node_->pCSD())->ccDC;
    for (i = 0; i < DCnum_; i++) {
      DCclasscode_[i] = *cc;
      cc++;
    }

    phaseclasscode_.resize(phasenum_,' ');
    cc = (node_->pCSD())->ccPH;
    for (i = 0; i < phasenum_; i++) {
      phaseclasscode_[i] = *cc;
      cc++;
    }

    ///
    /// Begin parsing the chemistry input XML file
    ///

    string msg;
    string xmlext = ".xml";
    size_t foundxml = Interfacefilename.find(xmlext);
    try {
      if (foundxml != string::npos) {
        parseDoc(Interfacefilename);
          
        micphasevolfrac_.resize(micphasenum_,0.0);
        micphasevolume_.resize(micphasenum_,0.0);
        micphasemass_.resize(micphasenum_,0.0);
        if (verbose_) {
            cout << " Setting micphasemass size to " << micphasenum_ << endl;
            cout.flush();
        }
        micphasemassdissolved_.resize(micphasenum_,0.0);
      } else {
        msg = "Not an XML file";
        throw FileException("ChemicalSystem","ChemicalSystem",Interfacefilename,msg);
      }
    }
    catch (FileException e) {
      throw e;
    }
    
    ///
    /// Set up the main map that correlates microstructure phases with GEM CSD phases
    ///

    mictotinitvolume_ = 0.0;
    for (unsigned int i = 0; i < micphasenum_; i++) {
      mic2phase_.insert(make_pair((int)i,micphasemembers_[i]));
      mic2DC_.insert(make_pair((int)i,micDCmembers_[i]));
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

    vector<double> solutionICmoles = getSolution();

    solut_->setICmoles(solutionICmoles);
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
        throw FileException("Solution","isInputFormatJSON",masterFileName,"Could not open");
    }

    string filetypeflag;
    if (in.peek() != EOF) {
        in >> filetypeflag;
    } else {
        throw FileException("Solution","isInputFormatJSON",masterFileName,"Bad or corrupt format");
    }
    in.close();
    if (filetypeflag == "-j") {
        return true;
    }
    return false;
}

void ChemicalSystem::getJSONFiles (const char *masterFileName,
                             string &dchname,
                             string &ipmname,
                             string &dbrname)
{
    ifstream in(masterFileName);
    if (!in) {
        throw FileException("Solution","getJSONFiles",masterFileName,"Could not open");
    }

    string filetypeflag;
    if (in.peek() != EOF) {
        in >> filetypeflag;
    } else {
        throw FileException("Solution","getJSONFiles",masterFileName,"Bad or corrupt format");
    }
    if (in.peek() != EOF) {
        in >> dchname;
    } else {
        throw FileException("Solution","getJSONFiles",masterFileName,"Bad or corrupt format");
    }
    if (in.peek() != EOF) {
        in >> ipmname;
    } else {
        throw FileException("Solution","getJSONFiles",masterFileName,"Bad or corrupt format");
    }
    if (in.peek() != EOF) {
        in >> dbrname;
    } else {
        throw FileException("Solution","getJSONFiles",masterFileName,"Bad or corrupt format");
    }

    /// Assuming that all the file names were read correctly, strip their quote marks
    
    dchname.erase(remove(dchname.begin(),dchname.end(),'"'),dchname.end());
    ipmname.erase(remove(ipmname.begin(),ipmname.end(),'"'),ipmname.end());
    dbrname.erase(remove(dbrname.begin(),dbrname.end(),'"'),dbrname.end());
}

vector<double> ChemicalSystem::getSolution ()
{
  ///
  /// Get IC moles for solution, which fully characterizes the
  /// composition of the aqueous solution
  ///
 
  double watermass = DCmoles_[getDCid("H2O@")] 
                 * DCmolarmass_[getDCid("H2O@")];

  if (verbose_) {
      cout << "water mass at the end of hydration is: " << watermass << endl;
  }
  vector<double> tempicmoles;
  tempicmoles.clear();
  tempicmoles.resize(ICnum_,0.0);
  for (unsigned int i = 0; i < DCnum_; i++) {
    char cc = getDCclasscode(i);
    if (cc == 'S' || cc == 'T') {
      double moles = node_->Get_cDC(i) * watermass * 1.0e-3;
      for(int j = 0; j < (ICnum_ - 1); j++) {
        tempicmoles[j] += moles * DCstoich_[i][j];
      }
    }
  }

  // Treat H2O separately

  double watermoles = DCmoles_[getDCid("H2O@")];
  for (int j = 0; j < ICnum_; j++) {
    if (ICname_[j] == "H") tempicmoles[j] += watermoles * 2;
    if (ICname_[j] == "O") tempicmoles[j] += watermoles;
  }

  return tempicmoles;  

}

void ChemicalSystem::parseDoc (const string &docname)
{
    string msg;
    PhaseData phasedata;
    xmlDocPtr doc;
    xmlChar *key;
    xmlNodePtr cur;
    cout.flush();
    doc = xmlParseFile(docname.c_str());

    /// Check if the xml file is valid

    string rxcsd = xsd_files_path;
    rxcsd+="/chemistry.xsd";
    // if (verbose_) {
    //    cout << "Chemistry xsd file is at " << rxcsd << endl;
    // }
    if(!is_xml_valid(doc,rxcsd.c_str())) {
        cout << "Chemistry xml is NOT valid" <<endl;
        cout.flush();
    } else if (verbose_) {
        cout << "Chemistry xml IS valid" << endl;
        cout.flush();
    }

    if (doc == NULL ) {
        msg = "XML file not parsed successfully";
        throw FileException("ChemicalSystem","parseDoc",docname,msg);
    }

    cur = xmlDocGetRootElement(doc);

    if (cur == NULL) {
        msg = "XML file is empty";
        xmlFreeDoc(doc);
        throw FileException("ChemicalSystem","parseDoc",docname,msg);
    }

    ///
    /// Go through the XML file one tag at a time
    ///
    /// The interface file contains information about each microstructure
    /// phase that is defined, including the list of GEM CSD phases that
    /// are to be associated with that phase, the phase's internal porosity,
    /// dissolved impurities, and visualization properties.
    ///

    cur = cur->xmlChildrenNode;
    int testnumentries = 0;
    int satstate = 1;
    saturated_ = true;
    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"numentries"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(testnumentries,st);
            xmlFree(key);
        } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"saturated"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(satstate,st);
            if (satstate == 0) saturated_ = false;
        } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"solution"))) {
            parseSolutionComp(doc, cur);
        } else if ((!xmlStrcmp(cur->name, (const xmlChar *)"phase"))) {
            parsePhase(doc, cur, testnumentries, phasedata);
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
    
    initial_solution_composition_.clear();

    string icname;
    int icid;
    double icconc;

    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"ICcomp"))) {
            parseICinSolution(doc, cur);
        }
        cur = cur->next;
    }

    return;
}

void ChemicalSystem::parseICinSolution (xmlDocPtr doc,
                                        xmlNodePtr cur)
{
    xmlChar *key;
    int icid = -1;
    string icname;
    double icconc = -1.0;

    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"name"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            from_string(icname,(char *)key);
            icid = getICid(icname);
            if (icconc > 0.0 && icid > 0) {
                initial_solution_composition_.insert(make_pair(icid,icconc));
                icid = -1;
                icconc = -1.0;
                icname = "Unknown";
            }

            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"conc"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(icconc,st);
            xmlFree(key);
            if (icconc > 0.0 && icid > 0) {
                initial_solution_composition_.insert(make_pair(icid,icconc));
                icid = -1;
                icconc = -1.0;
                icname = "Unknown";
            }
        }

        cur = cur->next;
    }

    return;
}


void ChemicalSystem::parsePhase (xmlDocPtr doc,
                                 xmlNodePtr cur,
                                 int numentries,
                                 PhaseData &phasedata)
{
    xmlChar *key;
    int proposedgemphaseid,proposedgemdcid;

    phasedata.gtmplt.clear();
    phasedata.atmpvec.clear();

    /// @note The affinity vector is always the same length, one entry for every
    /// microstructure phase, and the default value is zero.  Therefore
    /// the chemistry file does not need to zero affinity values
    
    phasedata.atmpvec.resize(numentries,0);
    phasedata.gemphaseid.clear();
    phasedata.gemdcid.clear();
    phasedata.gemphasename.clear();
    phasedata.gemdcname.clear();
    phasedata.gemphaseid.clear();
    phasedata.gemdcid.clear();
    phasedata.stresscalc = 0;
    phasedata.weak = 0;

    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"id"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.id,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"thamesname"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            phasedata.thamesname = st;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"gemphase_data"))) {
            parseGEMphasedata(doc, cur, phasedata);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"porosity"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.porosity,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"stresscalc"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.stresscalc,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"weak"))) {
            // Weak means the phase can be damaged by stress
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.weak,st);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"display_data"))) {
            parseDisplaydata(doc, cur, phasedata);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"impurity_data"))) {
            parseImpuritydata(doc, cur, phasedata);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"interface_data"))) {
            parseInterfacedata(doc, cur, phasedata);
        }
        cur = cur->next;
    }

    if (phasedata.stresscalc > 0) {
        stressphasename_.push_back(phasedata.thamesname);
        stressphaseid_.push_back(phasedata.id);
    }
    if (phasedata.weak > 0) {
        weakphasename_.push_back(phasedata.thamesname);
        weakphaseid_.push_back(phasedata.id);
    }
    if (phasedata.porosity < 1.0 && phasedata.porosity > 0.0) {
        porousphasename_.push_back(phasedata.thamesname);
        porousphaseid_.push_back(phasedata.id);
    }
    micphasename_.push_back(phasedata.thamesname);
    micid_.push_back(phasedata.id);
    micidlookup_.insert(make_pair(phasedata.thamesname,phasedata.id));
    randomgrowth_.push_back(phasedata.randomgrowth);
    affinity_.push_back(phasedata.atmpvec);

    // Growth template is based on positive affinities only
    
    growthtemplate_.push_back(calcGrowthtemplate(phasedata.atmpvec));

    porosity_.push_back(phasedata.porosity);
    grayscale_.push_back(phasedata.gray);
    color_.push_back(phasedata.colors);
    k2o_.push_back(phasedata.k2o);
    na2o_.push_back(phasedata.na2o);
    mgo_.push_back(phasedata.mgo);
    so3_.push_back(phasedata.so3);
    micphasemembers_.insert(make_pair(phasedata.id,phasedata.gemphaseid));
    micDCmembers_.insert(make_pair(phasedata.id,phasedata.gemdcid));
    micphasenum_++;

    return;
}

void ChemicalSystem::parseGEMphasedata (xmlDocPtr doc,
                                        xmlNodePtr cur,
                                        PhaseData &phasedata)
{
    xmlChar *key;
    cur = cur->xmlChildrenNode;

    int gemphaseid = 0;
    phasedata.gemphasedcmembers.clear();

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name,(const xmlChar *)"gemphasename"))) {
            key = xmlNodeListGetString(doc,cur->xmlChildrenNode,1);
            phasedata.gemphasename.push_back((char *)key);
            gemphaseid = getPhaseid((char *)key);
            phasedata.gemphaseid.push_back(gemphaseid);
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name,(const xmlChar *)"gemdcname"))) {
            key = xmlNodeListGetString(doc,cur->xmlChildrenNode,1);
            phasedata.gemdcname.push_back((char *)key);
            int dcid = getDCid((char *)key);
            phasedata.gemdcid.push_back(dcid);
            phasedata.gemphasedcmembers.push_back(dcid);
            xmlFree(key);
        }
        cur = cur->next;
    } 
    phaseDCmembers_.insert(make_pair(gemphaseid,phasedata.gemphasedcmembers));   
}

void ChemicalSystem::parseDisplaydata (xmlDocPtr doc,
                                       xmlNodePtr cur,
                                       PhaseData &phasedata)
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
            from_string(phasedata.gray,st);
            // if (verbose_) cout << "        gray = " << phasedata.gray << endl;
            xmlFree(key);
        }
        cur = cur->next;

    }

    phasedata.colors.clear();
    phasedata.colors.push_back(red);
    phasedata.colors.push_back(green);
    phasedata.colors.push_back(blue);

    return;
}

void ChemicalSystem::parseImpuritydata (xmlDocPtr doc,
                                        xmlNodePtr cur,
                                        PhaseData &phasedata)
{
    xmlChar *key;
    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"k2ocoeff"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.k2o,st);
            // if (verbose_) cout << "        k2ocoeff = " << phasedata.k2o << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"na2ocoeff"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.na2o,st);
            // if (verbose_) cout << "        na2ocoeff = " << phasedata.na2o << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"mgocoeff"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.mgo,st);
            // if (verbose_) cout << "        mgocoeff = " << phasedata.mgo << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"so3coeff"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.so3,st);
            // if (verbose_) cout << "        so3coeff = " << phasedata.so3 << endl;
            xmlFree(key);
        }
        cur = cur->next;

    }
    return;
}

void ChemicalSystem::parseInterfacedata (xmlDocPtr doc,
                                         xmlNodePtr cur,
                                         PhaseData &phasedata)
{

    xmlChar *key;
    cur = cur->xmlChildrenNode;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"randomgrowth"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(phasedata.randomgrowth,st);
            // if (verbose_) cout << "        randomgrowth = " << phasedata.randomgrowth << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"affinity"))) {
            // if (verbose_) cout << "        Parsing affinity data..." << endl;
            parseAffinitydata(doc,cur,phasedata);
            // if (verbose_) cout << "        Done parsing affinity data." << endl;
        }
        cur = cur->next;

    }
    return;
}

void ChemicalSystem::parseAffinitydata (xmlDocPtr doc,
                                        xmlNodePtr cur,
                                        PhaseData &phasedata)
{
    xmlChar *key;
    cur = cur->xmlChildrenNode;
    int testaftyid,testaftyval;

    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"affinityphaseid"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(testaftyid,st);
            // if (verbose_) cout << "            affinity id = " << testaftyid << endl;
            xmlFree(key);
        }
        if ((!xmlStrcmp(cur->name, (const xmlChar *)"affinityvalue"))) {
            key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            string st((char *)key);
            from_string(testaftyval,st);
            // if (verbose_) cout << "            affinity value = " << testaftyval << endl;
            xmlFree(key);
        }
        cur = cur->next;
    }
    phasedata.atmpvec[testaftyid] = ((int)testaftyval);
    return;
}

ChemicalSystem::ChemicalSystem (const ChemicalSystem &obj)
{

    ///
    /// This is a straightforward copy constructor
    /// Each member and vector is copied into the newly constructed object
    ///

    micphasenum_ = obj.getMicphasenum();
    micimpuritynum_ = obj.getMicimpuritynum();
    ICnum_ = obj.getICnum();
    DCnum_ = obj.getDCnum();
    phasenum_ = obj.getPhasenum();
    solutionphasenum_ = obj.getSolutionphasenum();
    micphasename_ = obj.getMicphasename();
    ICname_ = obj.getICname();
    DCname_ = obj.getDCname();
    phasename_ = obj.getPhasename();
    micid_ = obj.getMicid();
    c3sid_ = obj.getC3sid();
    c2sid_ = obj.getC2sid();
    c3aid_ = obj.getC3aid();
    c4afid_ = obj.getC4afid();
    gypsumid_ = obj.getGypsumid();
    randomgrowth_ = obj.getRandomgrowth();
    ICmoles_ = obj.getICmoles();
    DCmoles_ = obj.getDCmoles();
    ICmolarmass_ = obj.getICmolarmass();
    DCmolarmass_ = obj.getDCmolarmass();
    growthtemplate_ = obj.getGrowthtemplate();
    affinity_ = obj.getAffinity();
    micphasemembers_ = obj.getMicphasemembers();
    micDCmembers_ = obj.getMicDCmembers();
    porosity_ = obj.getPorosity();
    k2o_ = obj.getK2o();
    na2o_ = obj.getNa2o();
    mgo_ = obj.getMgo();
    so3_ = obj.getSo3();
    grayscale_ = obj.getGrayscale();
    color_ = obj.getColor();
    micidlookup_ = obj.getMicidlookup();
    ICidlookup_ = obj.getICidlookup();
    DCidlookup_ = obj.getDCidlookup();
    phaseidlookup_ = obj.getPhaseidlookup();
    mic2phase_ = obj.getMic2phase();
    mic2DC_ = obj.getMic2DC();
    kineticphase_ = obj.getKineticphase();
	thermophase_ = obj.getThermophase();
    ICclasscode_ = obj.getICclasscode();
    DCclasscode_ = obj.getDCclasscode();
    phaseclasscode_ = obj.getPhaseclasscode();
    DCstoich_ = obj.getDCstoich();
    phasestoich_ = obj.getPhasestoich();
    vphasestoich_ = obj.getVphasestoich();
    ICresiduals_ = obj.getICresiduals();
    ICchempot_ = obj.getICchempot();
    DCactivitycoeff_ = obj.getDCactivitycoeff();
    phasemoles_ = obj.getPhasemoles();
    ophasemoles_ = obj.getOphasemoles();
    phasemass_ = obj.getPhasemass();
    ophasemass_ = obj.getOphasemass();
    phasevolume_ = obj.getPhasevolume();
    ophasevolume_ = obj.getOphasevolume();
    carrier_ = obj.getCarrier();
    surfacearea_ = obj.getSurfacearea();
    DClowerlimit_ = obj.getDClowerlimit();
    DCupperlimit_ = obj.getDCupperlimit();
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
    ionicstrength_ = obj.getIonicstrength();
    Gs_ = obj.getGs();
    Hs_ = obj.getHs();
    nodehandle_ = obj.getNodehandle();
    nodestatus_ = obj.getNodestatus();
    iterdone_ = obj.getIterdone();
    micphasevolfrac_ = obj.getMicphasevolfrac();
    micphasevolume_ = obj.getMicphasevolume();
    mictotvolume_ = obj.getMictotvolume();
    mictotinitvolume_ = obj.getMictotinitvolume();
    micphasemass_ = obj.getMicphasemass();
    micphasemassdissolved_ = obj.getMicphasemassdissolved();
    micvoidvolume_ = obj.getMicvoidvolume();
    micvoidvolfrac_ = obj.getMicvoidvolfrac();
    verbose_ = obj.getVerbose();
    debug_ = obj.getDebug();
    warning_ = obj.getWarning();
}

ChemicalSystem::~ChemicalSystem ()
{
    ///
    /// Clear out the maps
    ///

    micidlookup_.clear();
    DCidlookup_.clear();
    ICidlookup_.clear();
    phaseidlookup_.clear();
    mic2phase_.clear();
    mic2DC_.clear();

    ///
    /// Clear out the vectors
    ///

    micphasename_.clear();
    ICname_.clear();
    DCname_.clear();
    phasename_.clear();
    micid_.clear();
    randomgrowth_.clear();
    DCstoich_.clear();
    growthtemplate_.clear();
    affinity_.clear();
    micphasemembers_.clear();
    micphasemembersvolfrac_.clear();
    micphasemass_.clear();
    micphasemassdissolved_.clear();
    micDCmembers_.clear();
    porosity_.clear();
    k2o_.clear();
    na2o_.clear();
    mgo_.clear();
    so3_.clear();
    grayscale_.clear();
    color_.clear();
    ICclasscode_.clear();
    DCclasscode_.clear();
    phaseclasscode_.clear();
    kineticphase_.clear();
	thermophase_.clear();
    vphasestoich_.clear();

    ///
    /// Free up the dynamically allocated memory
    ///

    delete[]DClowerlimit_;
    delete[]DCupperlimit_;
    delete[]surfacearea_;
    delete[]ophasemass_;
    delete[]ophasevolume_;
    delete[]ophasemoles_;
    delete[]phasemass_;
    delete[]phasevolume_;
    delete[]carrier_;
    delete[]phasemoles_;
    delete[]DCactivitycoeff_;
    delete[]DCmoles_;
    delete[]ICchempot_;
    delete[]ICresiduals_;
    delete[]ICmoles_;
    delete[]phasestoich_;

    delete node_;
}

void ChemicalSystem::getGEMPhasestoich ()
{
    double *arout = new double[ICnum_];
    for (long int i = 0; i < phasenum_; i++) {
        arout = node_->Ph_BC(i,arout);
        for (unsigned int j = 0; j < ICnum_; j++) {
            phasestoich_[(i * ICnum_) + j] = arout[j];
        }
    }
    delete[] arout;
}

void ChemicalSystem::getGEMVphasestoich ()
{
    double minval = 0.0;
    vphasestoich_.clear();
    vector<double> vplace;
    vplace.clear();
    vplace.resize(ICnum_,0.0);
    vphasestoich_.resize(phasenum_,vplace);
    int indexval,oval;
    for (unsigned int i = 0; i < phasenum_; i++) {
        if (phasename_[i] == "aq_gen") {

            ///
            /// Normalize to one mole of oxygen
            ///

            oval = (i * ICnum_) + getICid("O");
            if (phasestoich_[oval] > 0.0) {
                for (unsigned int j = 0; j < ICnum_; j++) {
                    indexval = (i * ICnum_) + j;
                    vphasestoich_[i][j] = (phasestoich_[indexval]/phasestoich_[oval]);
                }
            }
        } else {
            for (unsigned int j = 0; j < ICnum_; j++) {
                indexval = (i * ICnum_) + j;
                vphasestoich_[i][j] = (phasestoich_[indexval]/minval);
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

    for (i = 0; i < micphasenum_; i++) {
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
    if (i >= micphasenum_) {
        throw EOBException("ChemicalSystem","writeMember","micphasename_",micphasenum_,i);
    }

    ///
    /// Format the output for one phase
    ///

    stream << "------------------------------------------------------" << endl;
    stream << "DATA FOR MATERIAL " << i << ":" << endl;
    stream << "       Name = " << micphasename_[i] << endl;
    stream << "         Id = " << micid_[i] << endl;
    stream << "   Porosity = " << porosity_[i] << endl;
    stream << "------------------------------------------------------" << endl;
}

void ChemicalSystem::writeChemSys ()
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
    for(unsigned int i = 0; i < ICnum_; i++) {
        out << i << ")            Name: " << ICname_[i] << endl;
        out << "        classcode: " << ICclasscode_[i] << endl;
        out << "       molar mass: " << ICmolarmass_[i] << endl << endl;
    }

    out << "List of Dependent Components:" << endl << endl;
    for(unsigned int i = 0; i < DCnum_; i++) {
        out << i << ")            Name: " << DCname_[i] << endl;
        out << "        classcode: " << DCclasscode_[i] << endl;
        out << "       molar mass: " << DCmolarmass_[i] << endl << endl;
    }

    out << "List of Phases:" << endl << endl;
    for(unsigned int i = 0; i < phasenum_; i++) {
        out << i << ")            Name: " << phasename_[i] << endl;
        out << "        classcode: " << phaseclasscode_[i] << endl;
    }

    out << "List of Microstructure Phases:" << endl << endl;
    for(unsigned int i = 0; i < micphasenum_; i++) {
        out << i << ")       Name: " << micphasename_[i] << endl;
        out << "               id: " << micid_[i] << endl;
        out << "    random growth: " << randomgrowth_[i] << endl;
        for (j = 0; j < affinity_[i].size(); j++) {
            out << "        affinity to " << j << ": " << affinity_[i][j] << endl;
        }
        for (j = 0; j < growthtemplate_[i].size(); j++) {
            out << "        growthtemplate: " << growthtemplate_[i][j] << endl;
        }
        out << "         porosity: " << porosity_[i] << endl;
        out << "              k2o: " << k2o_[i] << endl;
        out << "             na2o: " << na2o_[i] << endl;
        out << "              mgo: " << mgo_[i] << endl;
        out << "              so3: " << so3_[i] << endl;
    }

    out.close();
    return;
}

void ChemicalSystem::writeChemSys (ostream &out)
{
    unsigned int j;

    ///
    /// First we will list details for the ICs
    ///

    out << "Report on the Material Database" << endl;
    out << "-------------------------------" << endl << endl;
    out << "List of Independent Components:" << endl << endl;
    for(unsigned int i = 0; i < ICnum_; i++) {
        out << i << ")            Name: " << ICname_[i] << endl;
        out << "        classcode: " << ICclasscode_[i] << endl;
        out << "       molar mass: " << ICmolarmass_[i] << endl << endl;
    }

    out << "List of Dependent Components:" << endl << endl;
    for(unsigned int i = 0; i < DCnum_; i++) {
        out << i << ")            Name: " << DCname_[i] << endl;
        out << "        classcode: " << DCclasscode_[i] << endl;
        out << "       molar mass: " << DCmolarmass_[i] << endl << endl;
    }

    out << "List of Phases:" << endl << endl;
    for(unsigned int i = 0; i < phasenum_; i++) {
        out << i << ")            Name: " << phasename_[i] << endl;
        out << "        classcode: " << phaseclasscode_[i] << endl;
    }

    out << "List of Microstructure Phases:" << endl << endl;
    for(unsigned int i = 0; i < micphasenum_; i++) {
        out << i << ")       Name: " << micphasename_[i] << endl;
        out << "               id: " << micid_[i] << endl;
        out << "    random growth: " << randomgrowth_[i] << endl;
        for (j = 0; j < affinity_[i].size(); j++) {
            out << "        affinity to " << j << ": " << affinity_[i][j] << endl;
        }
        for (j = 0; j < growthtemplate_[i].size(); j++) {
            out << "        growthtemplate: " << growthtemplate_[i][j] << endl;
        }
        out << "         porosity: " << porosity_[i] << endl;
        out << "              k2o: " << k2o_[i] << endl;
        out << "             na2o: " << na2o_[i] << endl;
        out << "              mgo: " << mgo_[i] << endl;
        out << "              so3: " << so3_[i] << endl;
    }

    return;
}

int ChemicalSystem::calculateState (double time,
                                    bool isfirst = false)
{
    int status = 0;
    string msg;
 
    // isfirst = true; 
      
    vector<double> oDCmoles;
    oDCmoles.clear();
    oDCmoles.resize(DCnum_,0.0);
  
    for (int i = 0; i < DCnum_; i++) {
      oDCmoles[i] = DCmoles_[i];
    }

    nodestatus_ = NEED_GEM_SIA;

    if (verbose_) {
        cout << "    Before calculateState, printing micphasevolumes" << endl;
        vector<double> micvols = getMicphasevolume();
        vector<string> micnames = getMicphasename();
        for (int i = 0; i < micvols.size(); ++i) {
            cout << "    Phase name " << micnames[i] << ": volume = " << micvols[i] << endl;
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


    node_->GEM_from_MT(nodehandle_,nodestatus_,T_,P_,Vs_,Ms_,
          ICmoles_,DCupperlimit_,DClowerlimit_,surfacearea_,
          DCmoles_,DCactivitycoeff_);

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
             << "    Going into ChemicalSystem::calculateState::GEM_run() (2), with isfirst "
             << isfirst << endl;
        cout.flush();
        writeICmoles();
    }

    /*
    writeDCmoles();
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
    /// Possible return values for nodestatus_:
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

    if (isfirst) {
        nodestatus_ = node_->GEM_run(true);
    } else {
        nodestatus_ = node_->GEM_run(true);
    }

    if (verbose_) {
        cout << "Done!  nodestatus is " << nodestatus_ << endl;
        cout.flush();
    }

    if (!(nodestatus_ == OK_GEM_AIA || nodestatus_ == OK_GEM_SIA)) {
        bool dothrow = false;
        cerr << "ERROR: Call to GEM_run in ChemicalSystem::calculateState had an issue..." << endl;
        cerr << "       nodestatus_ = " << nodestatus_;
        switch (nodestatus_) {
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
                cerr << msg << ", GEMS failed " << timesGEMfailed_ << " times" << endl;
                node_->GEM_print_ipm("IPM_dump.txt");
                timesGEMfailed_++;
                dothrow = (timesGEMfailed_ > maxGEMfails_) ? true : false;
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
                cerr << msg << ", GEMS failed " << timesGEMfailed_ << " times" << endl;
                node_->GEM_print_ipm("IPM_dump.txt");
                timesGEMfailed_++;
                dothrow = (timesGEMfailed_ > maxGEMfails_) ? true : false;
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
        timesGEMfailed_ = 0;
    }

    if (timesGEMfailed_ > 0) {
        if (verbose_) {
            cout << "Call to GEM_run has failed " << timesGEMfailed_
                 << " consecutive times.  Attempt this step again" << endl;
        }
        return timesGEMfailed_;
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

    node_->GEM_to_MT(nodehandle_,nodestatus_,iterdone_,Vs_,
            Ms_,Gs_,Hs_,ionicstrength_,pH_,pe_,Eh_,&ICresiduals_[0],
            &ICchempot_[0],&DCmoles_[0],&DCactivitycoeff_[0],&solutphasemoles_[0],
            &solutphasevolume_[0],&solutphasemass_[0],&solutphasestoich_[0],
            &carrier_[0],&surfacearea_[0],&solidstoich_[0]);

    if (verbose_) {
        cout << "Done!" << endl;
        cout << "after GEM_to_MT...Ms_ = " << Ms_ << endl;
        cout.flush();
    }

    /*
    writePhasemoles();
    */

    mictotvolume_ = 0.0;
    setPhasestoich();
    setVphasestoich();
    setPhasemass();
    setPhasevolume();
    setPhasemolarmass();
  
    if (verbose_) {
        cout << "%%%%%%%%%% Printing GEM Masses and Volumes in this Step %%%%%%%" << endl;
        for (long int myid = 0; myid < phasenum_; myid++) {
            cout << "Mass and volume of GEM phase " << node_->pCSD()->PHNL[myid] << " = "
                << phasemass_[myid] << " g and " << phasevolume_[myid] << " m3" << endl;
        }
        cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    }

    /// JWB: 2020 December 22
    /// Testing a change to the way we handle CSH molar volume adjustments
    /// Will let GEMS (thermodynamics) determine volumes of CSH and water
    /// Will adjust volume FRACTIONS of CSH and capillary porosity in the lattice part
    ///
    
    for (unsigned int i = 1; i < micphasenum_; i++) {
        if (verbose_) {
            cout << "Setting micphase amounts for " << i << " = " << micphasename_[i] << endl;
            cout.flush();
        }
        if (!isKineticphase(i)) {
            micphasemass_[i] = micphasevolume_[i] = 0.0;
            for (unsigned int j = 0; j < micphasemembers_[i].size(); j++) {
                if (verbose_) {
                    cout << "    Is NOT a KINETIC phase: is composed of "
                         << phasename_[micphasemembers_[i][j]]
                         << " having mass = " << phasemass_[micphasemembers_[i][j]]
                         << " and volume = " << phasevolume_[micphasemembers_[i][j]] << endl;
                    cout.flush();
                }
                micphasemass_[i] += phasemass_[micphasemembers_[i][j]];
                micphasevolume_[i] += phasevolume_[micphasemembers_[i][j]];
            }
            mictotvolume_ += micphasevolume_[i];

      /// End of test JWB: 2020 December 22
      ///

        } else {
            if (verbose_) {
                cout << "    IS a KINETIC phase: is composed of "
                     << phasename_[micphasemembers_[i][0]]
                     << "  having mass = " << micphasemass_[i]
                     << "  and volume = " << micphasevolume_[i] << endl;
                cout.flush();
            }

            ///
            /// micphasemass and micphasevolume for kinetic phases are
            /// already set in KineticModel::calculateKineticStep
            ///

            mictotvolume_ += micphasevolume_[i];
        }
    }

    if (saturated_) {   // System is saturated

        if (verbose_) cout << "Use water to saturate the porosity." << endl;
  
        if (mictotinitvolume_ > mictotvolume_) {
            double water_molarv, water_molesincr;
            for (int i = 0; i < micphasenum_; i++) {
                if (micphasename_[i] == "H2O") {
                    water_molarv = node_->DC_V0(getMic2DC(i,0), P_, T_);
                    water_molesincr = (mictotinitvolume_ - mictotvolume_) / water_molarv;
                    if (verbose_) {
                        cout << "water_molarv = " << water_molarv << endl;
                        cout << "volume increase of water is: "
                             << (mictotinitvolume_ - mictotvolume_) << endl;
                        cout << "water_molesincr = " << water_molesincr << endl;
                    }
                }
            }
            for (int i = 0; i < ICnum_; i++) {
                if (ICname_[i] == "H") ICmoles_[i] += water_molesincr * 2.0;
                if (ICname_[i] == "O") ICmoles_[i] += water_molesincr;
            }
        }

    }

    if (verbose_) {
        cout << "GEM volume change = "
             << 100.0 * (mictotvolume_ - mictotinitvolume_)/(mictotinitvolume_) << " %" << endl;
        cout.flush();
    }

    ///
    /// Calculate microstructure phase volume fractions
    ///
    
    for (int i = 0; i < micphasenum_; i++) {
        micphasevolfrac_[i] = micphasevolume_[i] / mictotvolume_;
    }

    setPhasestoich();
   
    ///
    /// Calculate driving force for growth or dissolution
    ///

    // double *soluticmoles;
    // soluticmoles = solut_->getICmoles();
    // for (int i = 0; i < DCnum_; i++) {
    //     char cc;
    //     cc = getDCclasscode(i);
    //     if (cc == 'O' || cc == 'I' || cc == 'J' || cc == 'M') {
    //         double dissolveDCmoles = oDCmoles[i] - DCmoles_[i];
    //         if (dissolveDCmoles > 0.0) {
    //             for (int j = 0; j < (ICnum_ - 1); j++) {
    //                 soluticmoles[j] += dissolveDCmoles * DCstoich_[i][j];
    //             }
    //         }
    //     }
    // }
    // for (int i = 0; i < ICnum_; i++) {
    //     solut_->setICmoles(i, soluticmoles[i]);
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

    vector<double> solutICmoles = getSolution();
    if (verbose_) cout << "Now update solution IC moles...";
    for (int i = 0; i < ICnum_; i++) {
        solut_->setICmoles(i, solutICmoles[i]);
    } 
    if (verbose_) cout << "Done." << endl;

    try {
        solut_->calculateState(false);  // Solution state was calculated for the first
                                        // time in the constructor.
    }
    catch (GEMException gex) {
        throw gex;
    }

    if (verbose_) {
        cout << "Leaving ChemicalSystem::calculateState now" << endl;
        cout.flush();
    }

    return timesGEMfailed_;
}
