/**
@file Solution.cc
@brief Define methods for the Solution class

*/
#include "Solution.h"
#include <iostream>

using namespace std;

Solution::Solution (const string &GEMfilename,
                    const string &GEMdbrname,
                    const bool verbose,
                    const bool debug)
{
    pH_ = 7.0;
    pe_ = 1.0;
    Eh_ = 0.0;
    T_ = 0.0;
    P_ = 101325.0;
    Vs_ = Ms_ = 1.0;
    Gs_ = Hs_ = 0.0;
    ionicstrength_ = 0.0;
    nodestatus_ = nodehandle_ = iterdone_ = 0;

    timesGEMfailed_ = 0;
    maxGEMfails_ = 3;

    verbose_ = verbose;
    debug_ = debug;
    jsonformat_ = false;

    ICnum_ = DCnum_ = phasenum_ = solutionphasenum_ = 0;
    ICname_.clear();
    DCname_.clear();
    phasename_.clear();
    SI_.clear();

    char *cGEMfilename = (char*)GEMfilename.c_str();
    char *cGEMdbrname = (char*)GEMdbrname.c_str();
    long int gemflag = 0;
    if (verbose_) {
        cout << "Solution:: Going into GEM_init to read chemical system definition file "
             << cGEMfilename << endl;
    }

    ///
    /// A new GEM3K node is allocated, completely separate from the one
    /// already allocated and used by the simulation's Chemicalsystem object
    ///

    node_ = new TNode();

    ///
    /// Use GEM3K to open and process the chemical system definition (CSD) file
    ///

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
        //         cout << "Detected JSON input file format for Solution GEM data files" << endl;
        //         cout.flush();
        //     }
        //     string json_dch,json_ipm,json_dbr;
        //     getJSONFiles(cGEMfilename,json_dch,json_ipm,json_dbr);
        //     gemflag = node_->GEM_init(json_dch,json_ipm,json_dbr);
        // } else {
        //     if (verbose) {
        //         cout << "Detected key-value input file format for Solution GEM data files" << endl;
        //         cout.flush();
        //     }
            gemflag = node_->GEM_init(cGEMfilename);
        // }
    }
    catch (FileException fex) {
        throw fex;
    }

    if (gemflag == 1) {
        cerr << "Bad return from GEM_init: " << GEMfilename << " missing or corrupt." << endl;
        exit(0);
    } else if (gemflag == -1) {
        cerr << "Bad return from GEM_init: internal memory allocation error." << endl;
        exit(0);
    }

    ///
    /// GEM3K has read the CSD file and now can be queried for class member values
    ///

    ICnum_ = (unsigned int)((node_->pCSD())->nIC);
    if (verbose_) cout << "ICnum_ is: " << ICnum_ << endl;
    DCnum_ = (unsigned int)((node_->pCSD())->nDC);
    if (verbose_) cout << "DCnum_ is: " << DCnum_ << endl;
    phasenum_ = (unsigned int)((node_->pCSD())->nPH);
    if (verbose_) cout << "phasenum_ is: " << phasenum_ << endl;
    solutionphasenum_ = (unsigned int)((node_->pCSD())->nPS);
    if (verbose_) cout << "solutionphasenum_ is: " << solutionphasenum_ << endl;

    ///
    /// Attempt to allocate memory for the various arrays
    ///

    string exmsg;
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
        exmsg = "DCupperlimit_";
        DCupperlimit_ = new double [DCnum_];
        exmsg = "DClowerlimit_";
        DClowerlimit_ = new double [DCnum_];
        exmsg = "phasemoles_";
        phasemoles_ = new double [phasenum_];
        exmsg = "phasemass_";
        phasemass_ = new double [phasenum_];
        exmsg = "phasevolume_";
        phasevolume_ = new double [phasenum_];
        exmsg = "surfacearea_";
        surfacearea_ = new double [phasenum_];
        exmsg = "carrier_";
        carrier_ = new double [solutionphasenum_];
        exmsg = "phasestoich_";
        phasestoich_ = new double [phasenum_ * ICnum_];
        exmsg = "solidstoich_";
        solidstoich_ = new double [ICnum_];
    }
    catch (bad_alloc& ba) {
        cout << endl << "Bad_alloc Exception Thrown:" << endl;
        cout << "    Details:" << endl;
        cout << "    Offending function Solution::Solution" << endl;
        cout << "    Error in allocating memory for array " << exmsg << endl;
        cerr << endl << "Bad_alloc Exception Thrown:" << endl;
        cerr << "    Details:" << endl;
        cerr << "    Offending function Solution::Solution" << endl;
        cerr << "    Error in allocating memory for array " << exmsg << endl;
        exit(0);
    }

    ///
    /// Initialize GEM
    ///

    (node_->pCNode())->NodeStatusCH = NEED_GEM_SIA;
    if (verbose_) {
        cout << "Solution::Constructor: Entering GEM_run (1) with node status = "
             << nodestatus_ << endl;
        cout << "NodeStatusCH is: " << (node_->pCNode())->NodeStatusCH << endl;
        cout.flush();
    }

    ///
    /// Attempt to run GEM
    ///
    /// The argument is false if we want to use activity coefficients and speciation
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
    
    nodestatus_ = node_->GEM_run(false);
    if (verbose_) {
        cout << "Solution::Constructor: Exited GEM_run with node status = "
             << nodestatus_ << endl;
        cout.flush();
    }
    if (!(nodestatus_ == OK_GEM_AIA || nodestatus_ == OK_GEM_SIA)) {
        bool dothrow = false;
        cerr << "ERROR: Call to GEM_run in Solution::constructor had an issue..." << endl;
        cerr << "       nodestatus_ = ";
        switch (nodestatus_) {
            case NEED_GEM_AIA:
                cout <<  " !!!!!Need GEM calc with auto initial approx (AIA)" << endl;
                exmsg = " !!!!!Need GEM calc with auto initial approx (AIA)";
                dothrow = false;
                break;
            case BAD_GEM_AIA:
                cout << " !!!!!Untrustworthy result with auto initial approx (AIA)" << endl;
                exmsg = " !!!!!Untrustworthy result with auto initial approx (AIA)",
                dothrow = false;
                break;
            case ERR_GEM_AIA:
                cout << " !!!!!Failed result with auto initial approx (AIA)" << endl;
                exmsg = " !!!!!Failed result with auto initial approx (AIA)";
                dothrow = false;
                break;
            case NEED_GEM_SIA:
                cout  << " !!!!!Need GEM calc with smart initial approx (SIA)" << endl;
                exmsg =  " !!!!!Need GEM calc with smart initial approx (SIA)";
                dothrow = false;
                break;
            case BAD_GEM_SIA:
                cout << " !!!!!Untrustworthy result with smart initial approx (SIA)" << endl;
                exmsg = " !!!!!Untrustworthy result with smart initial approx (SIA)";
                dothrow = false;
                break;
            case ERR_GEM_SIA:
                cout <<  " !!!!!Failed result with smart initial approx (SIA)" << endl;
                exmsg =  " !!!!!Failed result with smart initial approx (SIA)";
                dothrow = false;
                break;
            case T_ERROR_GEM:
                exmsg = " !!!!!Terminal GEM error; need restart";
                dothrow = true;
                break;
            case NO_GEM_SOLVER:
                 cout <<  " No GEM recalculation needed for node" << endl;
                exmsg =  " No GEM recalculation needed for node";
                dothrow = false;
                break;
        }
        if (dothrow) {
            throw GEMException("Solution","calculateState",exmsg);
        }
    }

    if (verbose_) {
        cout << "Solution::Constructor: Entering GEM_restore_MT (1) ..." << endl;
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
        cout << "T_ is: " << T_ << endl;
        cout << "Solution::Constructor: Entering GEM_to_MT (1)..." << endl;
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

    if (verbose_) {
        cout << "Done!" << endl;
        cout << "Solution::Constructor: Entering GEM_read_dbr (1)..." << endl;
        cout.flush();
    }

    ///
    /// The next function reads the DBR file with input system compsition,
    /// temperature, pressure, etc.  The DBR file must be compatible with
    /// the currently loaded IPM and DCH files.  It will return one of nine
    /// codes (see above for the GEM_init function comment for their meaning).
    ///
    /// @note the function arguments for GEM_read_dbr depend on whether the
    ///  input was given in JSON format or not, so we need two possible different
    ///  ways of calling it.  Seems weird, but that is the way the GEM library
    ///  is written as of ver. 3.8
    ///
   
    GEMS3KGenerator::IOModes type_f = GEMS3KGenerator::f_key_value;
    bool check_dch_compatibility = false;

    if (jsonformat_) {
        if (verbose_) {
            cout << "JSON GEM DBR file name is " << GEMdbrname << endl;
            cout.flush();
        }
        nodestatus_ =  node_->GEM_read_dbr(GEMdbrname,check_dch_compatibility);
    } else {
        if (verbose_) {
            cout << "Key-value GEM DBR file name is " << GEMdbrname << endl;
            cout.flush();
        }
        nodestatus_ =  node_->GEM_read_dbr(cGEMdbrname,type_f);
    }
    switch (nodestatus_) {
        case 1:     // JSON string is empty or corrupt
            cout << "ERROR: Call to GEM_read_dbr failed with nodestatus_ = "
                 << nodestatus_ << ", JSON string empty or corrupt" << endl;
            exit(1);
            break;
        case 2:     // DBR is incompatibile with DCH/IPM
            cout << "ERROR: Call to GEM_read_dbr failed with nodestatus_ = "
                 << nodestatus_ << ", DBR incompatible with DCH/IPM" << endl;
            exit(1);
            break;
        case -1:     // Memory allocation error
            cout << "ERROR: Call to GEM_read_dbr failed with nodestatus_ = "
                 << nodestatus_ << ", memory allocation error" << endl;
            exit(1);
            break;
        default:     // Okay
            break;
    }
    if (verbose_) {
        cout << "Done!" << endl;
        cout.flush();
    }

    ///
    /// Transfer phase names from GEM3K to member variables
    ///

    string string1;
    for (int i = 0; i < ICnum_; i++) {
        string1.assign(node_->xCH_to_IC_name(i));
        ICname_.push_back(string1);
    }
 
    string1.clear();  // This command may be unnecessary
    for (int i = 0; i < DCnum_; i++) {
        string1.assign(node_->xCH_to_DC_name(i));
        DCname_.push_back(string1);
    }
 
    string1.clear();  // This command may be unnecessary
    for (int i = 0; i < phasenum_; i++) {
        string1.assign(node_->xCH_to_Ph_name(i));
        phasename_.push_back(string1);
    }
}

Solution::~Solution ()
{
    ///
    /// Clear out the vectors
    ///

    ICname_.clear();
    DCname_.clear();
    phasename_.clear();
    SI_.clear();

    ///
    /// Delete previously allocated memory
    ///

    delete[]solidstoich_;
    delete[]phasestoich_;
    delete[]carrier_;
    delete[]surfacearea_;
    delete[]phasevolume_;
    delete[]phasemass_;
    delete[]phasemoles_;
    delete[]DClowerlimit_;
    delete[]DCupperlimit_;
    delete[]DCactivitycoeff_;
    delete[]DCmoles_;
    delete[]ICchempot_;
    delete[]ICresiduals_;
    delete[]ICmoles_;

    delete node_;
}

bool Solution::isInputFormatJSON (const char *masterFileName)
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

void Solution::getJSONFiles (const char *masterFileName,
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

void Solution::calculateState (bool isfirst)
{
    int status = 0;
    string msg;

    ///
    /// Load GEM data to the GEM3K library
    ///

    nodestatus_ = NEED_GEM_SIA;
    if (verbose_) {
        cout << "    Going into Solution::calculateState::GEM_from_MT (2) ..." << endl;
        cout.flush();
    }
 
    ///
    /// Next function loads the input data for the THAMES node into the
    /// instance of the DBR structure.  This call precedes the GEM_run call
    ///
    /// This function returns nothing and appears to be incapable of throwing
    /// an exception.
    ///
    /// @todo Check carefully if this function can throw and exception
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
  
        cout << "    Going into Solution::calculateState::GEM_run (2) with isfirst = "
             << isfirst << endl;
        cout.flush();
        cout << "    But first let's print the IC moles:" << endl;
        for (int i = 0; i < ICnum_; i++) {
            cout << ICname_[i] << ": " << ICmoles_[i] << " mol" << endl;
        }
        cout << endl;
        cout.flush();
    }

    ///
    /// Attempt to run a GEM calculation
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
        nodestatus_ = node_->GEM_run(false);
    } else {
        nodestatus_ = node_->GEM_run(true);
    }

    if (verbose_) {
        cout << "Done! nodestatus is " << nodestatus_ << endl;
        cout.flush();
    }
    if (!(nodestatus_ == OK_GEM_AIA || nodestatus_ == OK_GEM_SIA)) {
        bool dothrow = false;
        cerr << "ERROR: Call to GEM_run in Solution::calculateState had an issue..." << endl;
        cerr << "       nodestatus_ = " << nodestatus_;
        switch (nodestatus_) {
            case NEED_GEM_AIA:
                msg = "Need GEM calc with auto initial approx (AIA)";
                cerr << msg << endl;
                dothrow = false;
                break;
            case BAD_GEM_AIA:
                msg = "Untrustworthy result with auto initial approx (AIA)";
                cerr << msg << endl;
                dothrow = false;
                break;
            case ERR_GEM_AIA:
                msg = "Failed result with auto initial approx (AIA)";
                cerr << msg << " and has failed " << timesGEMfailed_ << " in a row" << endl;
                timesGEMfailed_++;
                // dothrow = (timesGEMfailed_ > maxGEMfails_) ? true : false;
                dothrow = false;
                break;
            case NEED_GEM_SIA:
                msg =  "Need GEM calc with smart initial approx (SIA)";
                cerr << msg << endl;
                dothrow = false;
                break;
            case BAD_GEM_SIA:
                msg = "Untrustworthy result with smart initial approx (SIA)";
                cerr << msg << endl;
                dothrow = false;
                break;
            case ERR_GEM_SIA:
                msg =  "Failed result with smart initial approx (SIA)";
                cerr << msg << " and has failed " << timesGEMfailed_ << " in a row" << endl;
                timesGEMfailed_++;
                // dothrow = (timesGEMfailed_ > maxGEMfails_) ? true : false;
                dothrow = false;
                break;
            case T_ERROR_GEM:
                msg = "Terminal GEM error; need restart";
                cerr << msg << endl;
                dothrow = true;
                break;
            case NO_GEM_SOLVER:
                msg =  "No GEM recalculation needed for node";
                cerr << msg << endl;
                dothrow = false;
                break;
        }
        if (dothrow) {
            throw GEMException("Solution","calculateState",msg);
        }
    } else {
        timesGEMfailed_ = 0;
    }

    ///
    /// Get the GEM data back from the GEM3K library, assuming that it ran
    /// without error.  Check for the errors first.
    ///

    if (verbose_) {
        cout << "    Going into Solution::calculateState::GEM_to_MT (2)...";
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

    if (verbose_) {
        cout << "Done!" << endl;
        cout << "after GEM_to_MT...Ms_ = " << Ms_ << endl;
        cout.flush();    
    }

    ///
    /// The GEM calculation did not precipitate or dissolve solid, so it
    /// stores the saturation index of each solid phase.  This now needs to
    /// be assigned to the corresponding array in the Solution object.
    ///

    setSI();

    return;

}

double Solution::calculateCrystrain (double SI,
                                     double porevolfrac,
                                     double Kp,
                                     double Ks)
{
    crystrain_ = 0.0;

    /*
    calculateState(false);
    */
 
    ///
    /// Crystallization pressure only exists if the solution is supersaturated,
    /// which means that \f$\beta > 1\f$.
    ///

    if (SI > 1.0) {

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
  
        if (verbose_) cout << "SI for this phase is: " << SI << endl; 

        ///
        /// Calculate the crystal mean curvature in equilibrium with the
        /// solution with this saturation index (Thompson-Freundlich effect)
        ///

        double kcr = Rg * T_ * log(SI) / (Vc * gamma);

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

        crystrain_ = (1.0 / (3.0 * Kp) - 1.0 / (3.0 * Ks)) * (porevolfrac * pa + pl);

        if (verbose_) cout << "crystrain is: " << crystrain_ << endl;
  
    }

    return crystrain_;
}
