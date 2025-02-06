/**
@file Lattice.cc
@brief Definition of methods for the Lattice class.

*/
#include "Lattice.h"
#include "Interface.h"
#include "RanGen.h"

Lattice::Lattice(ChemicalSystem *cs) : siteneighbors_(NN_NNN), chemSys_(cs) {
  xdim_ = ydim_ = zdim_ = 0;
  time_ = 0.0;
  temperature_ = REFTEMP; // in Kelvin (see global.h)
  oldtemp_ = REFTEMP;     // in Kelvin (see global.h)
  numsites_ = 0;
  resolution_ = REFRES; // in micrometers (see global.h)
  site_.clear();
  deptheffect_ = false;
  masterporevolume_.clear();
#ifdef DEBUG
  verbose_ = true;
#else
  verbose_ = false;
#endif
}

Lattice::Lattice(ChemicalSystem *cs, RanGen *rg, const string &fileName, const bool verbose,
                 const bool warning)
    : siteneighbors_(NN_NNN), chemSys_(cs), rg_(rg) {

  unsigned int i, j, k;
  unsigned int ii;
  string buff;
  int xn, yn, zn;
  unsigned int idn;
  unsigned int pid;
  string msg;

  numMicroPhases_ = chemSys_->getNumMicroPhases();

  ///
  /// Gather data from the periodic table for writing a .cfg file (alternative
  /// to .xyz)
  populateElementData();
  particRadius_ = 0.7;

  xdim_ = ydim_ = zdim_ = 0;
  time_ = 0.0;
  temperature_ = REFTEMP; // in Kelvin (see global.h)
  oldtemp_ = REFTEMP;     // in Kelvin (see global.h)
  numsites_ = 0;
  resolution_ = REFRES; // in micrometers (see global.h)
  site_.clear();
  deptheffect_ = true;
  masterporevolume_.clear();

#ifdef DEBUG
  verbose_ = true;
  warning_ = true;
  cout << "Lattice::Lattice Constructor" << endl;
  cout.flush();
#else
  verbose_ = verbose;
  warning_ = warning;
#endif

  ///
  /// Open the microstructure input file and process it
  ///

  ifstream in(fileName.c_str());
  if (!in) {
    throw FileException("Lattice", "Lattice", fileName, "Could not open");
  }

  in >> buff;
  if (buff == VERSIONSTRING) {
    in >> version_;
    cout << endl << "Lattice::Lattice - .img  input file  generated using "
                    "VCCTL Version : " << version_ << endl;
    in >> buff; // X size string identifier
    in >> xdim_;
    in >> buff; // Y size string identifier
    in >> ydim_;
    in >> buff; // Z size string identifier
    in >> zdim_;
    double testres;
    in >> buff; // Voxel resolution identifier
    in >> testres;
    setResolution(testres);
  } else {

    ///
    /// Image file generated prior to VCCTL Version 3.0.
    /// Allow backward compatibility by defaulting system
    /// size to 100 and resolution to 1.0 micrometers
    ///
    cout << endl << "Lattice::Lattice - .img input file "
         "generated using VCCTL Version prior to VCCTL Version 3.0.0" << endl;
    version_ = "2.0";
    double testres = 1.0;
    setResolution(testres);
    xdim_ = ydim_ = zdim_ = 100;
  }

  ostringstream ostrMAJ, ostrMIN;
  ostrMAJ << VERSION_MAJOR;
  ostrMIN << VERSION_MINOR;
  string majVers(ostrMAJ.str());
  string minVers(ostrMIN.str());
  thamesVersion_ = majVers + "." + minVers + "." + VERSIONBUGFIX;
  cout << endl << "Lattice::Lattice - .img output files generated using "
                  "THAMES Version: " << thamesVersion_ << endl;

  ///
  /// Print out the microstructure size and characteristics
  ///

  if (verbose_) {
    cout << "Lattice::Lattice Read microstructure file header..." << endl;
    cout << "Lattice::Lattice     VCCTL Version = " << version_ << endl;
    cout << "Lattice::Lattice     xdim_ = " << xdim_ << endl;
    cout << "Lattice::Lattice     ydim_ = " << ydim_ << endl;
    cout << "Lattice::Lattice     zdim_ = " << zdim_ << endl;
    cout.flush();
  }
  numsites_ = (unsigned int)(xdim_ * ydim_ * zdim_);
  if (verbose_) {
    cout << "Lattice::Lattice    numsites_ = " << numsites_ << endl;
    cout << "Lattice::Lattice    resolution_ = " << resolution_ << endl;
    cout.flush();
  }

  ///
  /// Allocate a random number generator object and seed it
  ///

  latticeRNGseed_ = -76241; //-142234;
  cout << endl << "Lattice::Lattice   latticeRNGseed_ = " << latticeRNGseed_ << endl;
  //setRNGseed(latticeRNGseed_);
  try {
    rg_->setSeed(latticeRNGseed_);
  } catch (bad_alloc &ba) {
    cout << "Lattice constructor failed when allocating rg_";
    cout.flush();
    exit(1);
  }
  // //cout << endl << "Lattice::Lattice   rg_->getSeed() = " << rg_->getSeed() << endl;
  numRNGcall_0_ = 0;
  numRNGcallLONGMAX_ = 0;
  lastRNG_ = 1.e-16;

  growthInterfaceSize_.clear();
  dissolutionInterfaceSize_.clear();
  growthInterfaceSize_.resize(numMicroPhases_, 0);
  dissolutionInterfaceSize_.resize(numMicroPhases_, 0);
  //growthVector.clear();
  //dissolutionVector.clear();

  //count_.clear();
  count_.resize(numMicroPhases_, 0);

  //SI_.clear();
  SI_.resize(numMicroPhases_, 1.0);

  expansion_.clear();
  expansion_coordin_.clear();

  waterchange_ = 0.0;

  ///
  /// Add the Site objects to the `site_` array, defining the structure
  /// of the lattice.  Not yet assigning a phase to the site because
  /// the file has not yet been read
  ///

  for (k = 0; k < zdim_; k++) {
    for (j = 0; j < ydim_; j++) {
      for (i = 0; i < xdim_; i++) {
        addSite(i, j, k);
      }
    }
  }

  ///
  /// Set up the array of neighbors for each site
  ///

  for (ii = 0; ii < numsites_; ii++) {

    for (j = 0; j < NN_NNN; j++) {
      ///
      /// Store up to 2nd nearest neighbors
      ///

      switch (j) {
      case 0: // West neighbor
        xn = site_[ii].getX() - 1;
        if (xn < 0 && BC != 1)
          xn += xdim_;
        if (xn < 0 && BC == 1)
          xn = 0;
        yn = site_[ii].getY();
        zn = site_[ii].getZ();
        break;
      case 1: // East neighbor
        xn = site_[ii].getX() + 1;
        if (xn >= xdim_ && BC != 1)
          xn -= xdim_;
        if (xn >= xdim_ && BC == 1)
          xn -= 1;
        yn = site_[ii].getY();
        zn = site_[ii].getZ();
        break;
      case 2: // South neighbor
        yn = site_[ii].getY() - 1;
        if (yn < 0 && BC != 2)
          yn += ydim_;
        if (yn < 0 && BC == 2)
          yn = 0;
        xn = site_[ii].getX();
        zn = site_[ii].getZ();
        break;
      case 3: // North neighbor
        yn = site_[ii].getY() + 1;
        if (yn >= ydim_ && BC != 2)
          yn -= ydim_;
        if (yn >= ydim_ && BC == 2)
          yn -= 1;
        xn = site_[ii].getX();
        zn = site_[ii].getZ();
        break;
      case 4: // Down neighbor
        zn = site_[ii].getZ() - 1;
        if (zn < 0 && BC != 3)
          zn += zdim_;
        if (zn < 0 && BC == 3)
          zn = 0;
        xn = site_[ii].getX();
        yn = site_[ii].getY();
        break;
      case 5: // Up neighbor
        zn = site_[ii].getZ() + 1;
        if (zn >= zdim_ && BC != 3)
          zn -= zdim_;
        if (zn >= zdim_ && BC == 3)
          zn -= 1;
        xn = site_[ii].getX();
        yn = site_[ii].getY();
        break;
      case 6: // Southwest neighbor
        xn = site_[ii].getX() - 1;
        if (xn < 0 && BC != 1)
          xn += xdim_;
        if (xn < 0 && BC == 1)
          xn = 0;
        yn = site_[ii].getY() - 1;
        if (yn < 0 && BC != 2)
          yn += ydim_;
        if (yn < 0 && BC == 2)
          yn = 0;
        zn = site_[ii].getZ();
        break;
      case 7: // Northwest neighbor
        xn = site_[ii].getX() - 1;
        if (xn < 0 && BC != 1)
          xn += xdim_;
        if (xn < 0 && BC == 1)
          xn = 0;
        yn = site_[ii].getY() + 1;
        if (yn >= ydim_ && BC != 2)
          yn -= ydim_;
        if (yn >= ydim_ && BC == 2)
          yn -= 1;
        zn = site_[ii].getZ();
        break;
      case 8: // Northeast neighbor
        xn = site_[ii].getX() + 1;
        if (xn >= xdim_ && BC != 1)
          xn -= xdim_;
        if (xn >= xdim_ && BC == 1)
          xn -= 1;
        yn = site_[ii].getY() + 1;
        if (yn >= ydim_ && BC != 2)
          yn -= ydim_;
        if (yn >= ydim_ && BC == 2)
          yn -= 1;
        zn = site_[ii].getZ();
        break;
      case 9: // Southeast neighbor
        xn = site_[ii].getX() + 1;
        if (xn >= xdim_ && BC != 1)
          xn -= xdim_;
        if (xn >= xdim_ && BC == 1)
          xn -= 1;
        yn = site_[ii].getY() - 1;
        if (yn < 0 && BC != 2)
          yn += ydim_;
        if (yn < 0 && BC == 2)
          yn = 0;
        zn = site_[ii].getZ();
        break;
      case 10: // Downwest neighbor
        xn = site_[ii].getX() - 1;
        if (xn < 0 && BC != 1)
          xn += xdim_;
        if (xn < 0 && BC == 1)
          xn = 0;
        yn = site_[ii].getY();
        zn = site_[ii].getZ() - 1;
        if (zn < 0 && BC != 3)
          zn += zdim_;
        if (zn < 0 && BC == 3)
          zn = 0;
        break;
      case 11: // Downnorth neighbor
        xn = site_[ii].getX();
        yn = site_[ii].getY() + 1;
        if (yn >= ydim_ && BC != 2)
          yn -= ydim_;
        if (yn >= ydim_ && BC == 2)
          yn -= 1;
        zn = site_[ii].getZ() - 1;
        if (zn < 0 && BC != 3)
          zn += zdim_;
        if (zn < 0 && BC == 3)
          zn = 0;
        break;
      case 12: // Downeast neighbor
        xn = site_[ii].getX() + 1;
        if (xn >= xdim_ && BC != 1)
          xn -= xdim_;
        if (xn >= xdim_ && BC == 1)
          xn -= 1;
        yn = site_[ii].getY();
        zn = site_[ii].getZ() - 1;
        if (zn < 0 && BC != 3)
          zn += zdim_;
        if (zn < 0 && BC == 3)
          zn = 0;
        break;
      case 13: // Downsouth neighbor
        xn = site_[ii].getX();
        yn = site_[ii].getY() - 1;
        if (yn < 0 && BC != 2)
          yn += ydim_;
        if (yn < 0 && BC == 2)
          yn = 0;
        zn = site_[ii].getZ() - 1;
        if (zn < 0 && BC != 3)
          zn += zdim_;
        if (zn < 0 && BC == 3)
          zn = 0;
        break;
      case 14: // Upwest neighbor
        xn = site_[ii].getX() - 1;
        if (xn < 0 && BC != 1)
          xn += xdim_;
        if (xn < 0 && BC == 1)
          xn = 0;
        yn = site_[ii].getY();
        zn = site_[ii].getZ() + 1;
        if (zn >= zdim_ && BC != 3)
          zn -= zdim_;
        if (zn >= zdim_ && BC == 3)
          zn -= 1;
        break;
      case 15: // Upnorth neighbor
        xn = site_[ii].getX();
        yn = site_[ii].getY() + 1;
        if (yn >= ydim_ && BC != 2)
          yn -= ydim_;
        if (yn >= ydim_ && BC == 2)
          yn -= 1;
        zn = site_[ii].getZ() + 1;
        if (zn >= zdim_ && BC != 3)
          zn -= zdim_;
        if (zn >= zdim_ && BC == 3)
          zn -= 1;
        break;
      case 16: // Upeast neighbor
        xn = site_[ii].getX() + 1;
        if (xn >= xdim_ && BC != 1)
          xn -= xdim_;
        if (xn >= xdim_ && BC == 1)
          xn -= 1;
        yn = site_[ii].getY();
        zn = site_[ii].getZ() + 1;
        if (zn >= zdim_ && BC != 3)
          zn -= zdim_;
        if (zn >= zdim_ && BC == 3)
          zn -= 1;
        break;
      case 17: // Upsouth neighbor
        xn = site_[ii].getX();
        yn = site_[ii].getY() - 1;
        if (yn < 0 && BC != 2)
          yn += ydim_;
        if (yn < 0 && BC == 2)
          yn = 0;
        zn = site_[ii].getZ() + 1;
        if (zn >= zdim_ && BC != 3)
          zn -= zdim_;
        if (zn >= zdim_ && BC == 3)
          zn -= 1;
        break;
      }

      idn = (unsigned int)(xn + (xdim_ * yn) + ((xdim_ * ydim_) * zn));
      site_[ii].setNb(j, &site_[idn]);
    }
  }

  ///
  /// This loop reads the phase for each site and assigns it,
  /// also updating the count of the phase.
  ///

  for (i = 0; i < numsites_; i++) {
    in >> pid;
    site_[i].setMicroPhaseId(pid);
    count_[pid]++;
  }

  ///
  /// Done with the input microstructure file, so close it.
  ///

  in.close();

  ///
  /// With the phase counts known, calculate phase volume fractions
  /// We cannot calculate the mass fractions until we know the phase identities
  /// Phase identities are determined by the KineticModel
  ///

  vector<double> microPhaseMass(numMicroPhases_, 0.0);

  volumefraction_.clear();
  volumefraction_.resize(numMicroPhases_, 0.0);

  initvolumefraction_.clear();
  initvolumefraction_.resize(numMicroPhases_, 0.0);

  double vfrac, mfrac, capfrac, molarMass, molarVolume, density;
  double solidMass = 0.0;
  double cementMass = 0.0;
  int microPhaseId = 0;
  int DCId = 0;
  string myname;
  try {
    if (verbose_) {
      cout << "Lattice::Lattice Calculating volume fractions now ..." << endl;
      for (ii = 0; ii < numMicroPhases_; ii++) {
        cout << "Micro phase "
             << chemSys_->getMicroPhaseName(chemSys_->getMicroPhaseId(ii))
             << ", count = " << count_[ii] << " of " << site_.size() << endl;
      }
      cout.flush();
    }

    if (numsites_ > 0) {
      for (ii = 0; ii < numMicroPhases_; ii++) {
        microPhaseId = chemSys_->getMicroPhaseId(ii);
        myname = chemSys_->getMicroPhaseName(microPhaseId);
        vfrac = ((double)count_[ii]) / numsites_;
        setVolumefraction(microPhaseId, vfrac);
        setInitvolumefraction(microPhaseId, vfrac);
        if (verbose_) {
          cout << "Lattice::Lattice ii = " << ii << ", microPhase = " << myname
               << ", count_[" << ii<< "] = " << count_[ii]
               << ", volume fraction = " << vfrac << endl;
          cout.flush();
        }
        if (microPhaseId == ELECTROLYTEID) {
          DCId = chemSys_->getDCId("H2O@");
        } else if (microPhaseId != VOIDID) {
          DCId = chemSys_->getMicroPhaseDCMembers(microPhaseId, 0);
        }
        if (microPhaseId != VOIDID) {
          molarMass = chemSys_->getDCMolarMass(DCId);     // g/mol
          molarVolume = chemSys_->getDCMolarVolume(DCId); // m3/mol
          density = 0.0;
          if (molarVolume > 1.0e-12) {
            density = molarMass / molarVolume / 1.0e6; // g/cm3
          }
          microPhaseMass[microPhaseId] = vfrac * density;
          if (microPhaseId != ELECTROLYTEID) {
            solidMass += microPhaseMass[microPhaseId];
            if (chemSys_->getCementComponent(microPhaseId))
              cementMass += microPhaseMass[microPhaseId];
          }
          if (verbose_ && vfrac > 0.0) {
            // if (vfrac > 0.0) {
            cout << ii << "\tLattice::Lattice Phase "
                 << chemSys_->getMicroPhaseName(microPhaseId)
                 << "\tmicroPhaseId: " << microPhaseId << endl;
            cout << "Lattice::Lattice     Molar mass = " << molarMass
                 << " g/mol" << endl;
            cout << "Lattice::Lattice     Molar volume = " << molarVolume
                 << " m3/mol" << endl;
            cout << "Lattice::Lattice     Density = " << density << " g/cm3"
                 << endl;
            cout << "Lattice::Lattice     Volume Fraction = " << vfrac << endl;
            cout << "Lattice::Lattice     Mass density = " << (vfrac * density)
                 << " g/cm3 of system" << endl;
            cout.flush();
          }
        }
      }
    } else {
      msg = "Divide by zero error:  site_.size() = 0";
      throw FloatException("Lattice", "Lattice", msg);
    }

    // for (ii = 0; ii < numMicroPhases_; ii++) {
    //   cout << "  ii : " << ii << "\tcementComponent : "
    //        << chemSys_->getCementComponent(ii) << endl;
    // }

    // Set the water-solids mass ratio based on the initial microstructure

    initSolidMass_ = solidMass;

    wsratio_ = microPhaseMass[ELECTROLYTEID] / solidMass;
    wcratio_ = microPhaseMass[ELECTROLYTEID] / cementMass;

    if (verbose_) {
      cout << "Lattice::Lattice Microstructure w/s = " << wsratio_ << endl;
      cout << "Lattice::Lattice Mass of water = "
           << microPhaseMass[ELECTROLYTEID]
           << ", mass of solids = " << solidMass << endl;
      cout.flush();
    }

    // set to zero ALL ICs/DCs in the system (ChemicalSystem)
    int numICs = chemSys_->getNumICs();
    for (int i = 0; i < numICs; i++) {
      chemSys_->setICMoles(i, 0.0);
      // cout << "   " << i << "\t" << chemSys_->getICMoles(i)
      //     << "\t" << chemSys_->getICName(i) << endl;
    }
    int numDCs = chemSys_->getNumDCs();
    for (int i = 0; i < numDCs; i++) {
      chemSys_->setDCMoles(i, 0.0);
      //chemSys_->setDCLowerLimit(i, 0.0);
      // cout << "   " << i << "\t" << chemSys_->getDCMoles(i)
      //      << "\t" << chemSys_->getDCClassCode(i) << "\t" <<
      //      chemSys_->getDCName(i) << endl;
    }

    // Next we set the initial normalized phase masses, microstructure
    // phase volume, and subvoxel porosity.  This is all triggered when we set
    // the mass.

    // calc porosity in
    // normalizePhaseMasses->setMicroPhaseMass->setMicroPhaseVolume->calcMicroPhasePorosity
    normalizePhaseMasses(microPhaseMass, cementMass, solidMass);

    // Set the initial total volume of the microstructure

    double totmicvol = 0.0;
    //int totCount_ = 0;
    // cout << endl << "count_ : " << endl;
    for (int i = 0; i < numMicroPhases_; i++) {
      microPhaseId = chemSys_->getMicroPhaseId(i);
      if (microPhaseId != VOIDID) {
        // cout << i << "\t" << chemSys_->getMicroPhaseVolume(microPhaseId)
        //      << "\t" << chemSys_->getMicroPhaseName(i) << " / " << count_[i]
        //      << endl ;
        totmicvol += chemSys_->getMicroPhaseVolume(microPhaseId);
        // totCount_ += count_[i];
      }
    }
    // cout << "totCount_/totmicvol: " << totCount_ << " / " << totmicvol
    // <<endl;

    chemSys_->setInitMicroVolume(totmicvol);
    initialmicrostructurevolume_ = totmicvol;

    // Initially assume that all free water and void space
    // is capillary volume

    capillaryporevolumefraction_ =
        getVolumefraction(ELECTROLYTEID) + getVolumefraction(VOIDID);

  } catch (FloatException fex) {
    fex.printException();
    exit(1);
  }

  // calc & set wmc
  int phId;
  string nameMicroPhaseTh;
  double rng;
  for (i = 0; i < numsites_; i++) {
    // stId = site_[i].getId();
    phId = site_[i].getMicroPhaseId();
    nameMicroPhaseTh = chemSys_->getMicroPhaseName(phId);
    if (nameMicroPhaseTh == "Electrolyte") {
      site_[i].setWmc0(1.);
    } else if (nameMicroPhaseTh == "CSHQ") {
      rng = callRNG();
      if (rng >= thrPorosityCSH) {
        site_[i].setWmc0(chemSys_->getMicroPhasePorosity(phId));
      } else {
        site_[i].setWmc0(0.);
      }
    } else {
      site_[i].setWmc0(0.);
    }
  }
  double wmcLoc;
  int nbIdLoc;
  for (i = 0; i < numsites_; i++) {
    site_[i].setVisit(0);
    wmcLoc = site_[i].getWmc0();
    for (j = 0; j < NN_NNN; j++) {
      nbIdLoc = (site_[i].nb(j))->getId();
      wmcLoc += site_[nbIdLoc].getWmc0();
    }
    site_[i].setWmc(wmcLoc);
  }

  // voxels without contact with electrolyte i.e. voxels having low probability
  // to dissolve in this step
  // findIsolatedClusters();
}

Lattice::~Lattice() {
  interface_.clear();
  site_.clear();
  masterporevolume_.clear();
  //delete rg_;
}

void Lattice::addSite(const unsigned int x, const unsigned int y,
                      const unsigned int z) {
  string msg;

  try {
    if (x >= xdim_ || x < 0) {
      throw EOBException("Lattice", "addSite", "site_", xdim_, x);
    } else if (y >= ydim_ || y < 0) {
      throw EOBException("Lattice", "addSite", "site_", ydim_, y);
    } else if (z >= zdim_ || z < 0) {
      throw EOBException("Lattice", "addSite", "site_", zdim_, z);
    }
  } catch (EOBException ex) {
    ex.printException();
    exit(1);
  }

  // Assume that the site is WATER by default.  This will break if
  // there is ever a sytem without water.

  site_.push_back(Site(x, y, z, xdim_, ydim_, zdim_, siteneighbors_, chemSys_));
}

void Lattice::normalizePhaseMasses(vector<double> microPhaseMass,
                                   double cementMass, double solidMass) {
  int microPhaseId, DCId;
  double pscaledMass = 0.0;
  double molarMass;
  double totalSolidMass = 0, totalCementMass = 0;

  for (int i = 0; i < numMicroPhases_; i++) {
    microPhaseId = chemSys_->getMicroPhaseId(i);
    if (microPhaseId == ELECTROLYTEID) {
      int waterId = chemSys_->getDCId("H2O@");
      // double waterMolarMass = chemSys_->getDCMolarMass(waterId);
      molarMass = chemSys_->getDCMolarMass(waterId);
      pscaledMass = wsratio_ * 100.0; // Mass of solids scaled to 100 g now

      chemSys_->setDCMoles(waterId, (pscaledMass / molarMass));
      if (verbose_) {
        cout << "Lattice::normalizePhaseMasses Setting initial micphase mass "
                "and volume of "
             << chemSys_->getMicroPhaseName(ELECTROLYTEID) << endl;
        cout.flush();
      }

      chemSys_->setMicroPhaseMass(ELECTROLYTEID, pscaledMass);
      chemSys_->setMicroPhaseMassDissolved(ELECTROLYTEID, 0.0);

    } else if (microPhaseId != VOIDID) {
      pscaledMass = microPhaseMass[microPhaseId] * 100.0 / solidMass;
      DCId = chemSys_->getMicroPhaseDCMembers(microPhaseId, 0);
      chemSys_->setDC_to_MPhID(DCId, microPhaseId);
      molarMass = chemSys_->getDCMolarMass(DCId);
      if (verbose_) {
        cout << "Lattice::normalizePhaseMasses Microstructure scaled mass of "
             << i << "   " << chemSys_->getMicroPhaseName(microPhaseId) << " ("
             << microPhaseId << ") = " << microPhaseMass[microPhaseId]
             << " g out of " << solidMass << " g total" << endl;
        // Setting the phase mass will also automatically calculate the phase
        // volume
        cout.flush();
      }

      totalSolidMass += pscaledMass;
      if (chemSys_->getCementComponent(microPhaseId))
        totalCementMass += pscaledMass;
      chemSys_->setMicroPhaseMass(microPhaseId, pscaledMass);

      chemSys_->setDCMoles(DCId, (pscaledMass / molarMass));

      chemSys_->setMicroPhaseMassDissolved(microPhaseId, 0.0);
    }
  }

  chemSys_->setInitScaledCementMass(cementMass * 100 / solidMass);
  // cout << "normalizePhaseMasses totalSolidMass/totalCementMass = "
  //      << totalSolidMass << " / " << totalCementMass << endl;

  // Up to this point we could not really handle volume of void space, but
  // now we can do so in proportion to electrolyte volume

  double vfv = getVolumefraction(VOIDID);
  double vfe = getVolumefraction(ELECTROLYTEID);
  double ve = chemSys_->getMicroPhaseVolume(ELECTROLYTEID);

  chemSys_->setMicroPhaseVolume(VOIDID, (ve * vfv / vfe));

  // cout << endl << "normalizePhaseMasses DCs:" << endl;
  // for (int i = 0; i < chemSys_->getNumDCs(); i++){
  //   cout << "   " << i << "\t" << chemSys_->getDCMoles(i)
  //        << "\t" << chemSys_->getDCClassCode(i)
  //        << "\t" << chemSys_->getDCName(i) << endl;
  // }
  // cout << endl << "normalizePhaseMasses DCs:" << endl << endl;

  return;
}

void Lattice::findInterfaces(void) {
  unsigned int i, kk;
  unsigned int k;
  vector<Site *> gsite, dsite;
  vector<Site *>::iterator beginLocation, endLocation;
  vector<Isite>::iterator start, end;
  int stId;
  int gvsize = gsite.size();
  int dvsize = dsite.size();

  ///
  /// An interface must have at least one adjacent site that is water or void
  ///
  ///
  cout << endl << "   Lattice::findInterfaces:" << endl;

  interface_.clear();
  for (i = 0; i < numMicroPhases_; i++) {
    if (i != ELECTROLYTEID && i != VOIDID) { // Solid phase of some kind
      gsite.clear();
      dsite.clear();
      for (k = 0; k < numsites_; k++) {
        if (site_[k].getWmc() > 0) { // Is there some water nearby?
          if ((site_[k].getMicroPhaseId() == i)) {
            dsite.push_back(&site_[k]);
            //site_[k].setDissolutionSite(i);
            //if (chemSys_->getAffinity(site_[k].getMicroPhaseId(),site_[k].getMicroPhaseId()) > 10) {//itself affinity!
            for (kk = 0; kk < NN_NNN; kk++) {
              if ((site_[k].nb(kk))->getMicroPhaseId() == ELECTROLYTEID) {
                gsite.push_back(site_[k].nb(kk));
                site_[k].nb(kk)->setGrowthSite(i);
              }
            }
            //}

            /// @note There is no reason to make the phase
            /// itself be a template for itself, nor water be
            /// a growth template for any phase, because we already
            /// tested for those above.

          } else if (chemSys_->isGrowthTemplate(i,
                                                site_[k].getMicroPhaseId())) {
            for (kk = 0; kk < site_[k].nbSize(1); kk++) {
              if ((site_[k].nb(kk))->getMicroPhaseId() == ELECTROLYTEID) {
                gsite.push_back(site_[k].nb(kk));
                site_[k].nb(kk)->setGrowthSite(i);
              }
            }
          }
        }
      }

      // Eliminate duplicate values from the growth site vector
      gvsize = gsite.size();
      if (gvsize > 0) {
        sort(gsite.begin(), gsite.end());
        beginLocation = gsite.begin();
        endLocation = unique(gsite.begin(), gsite.end());
        gsite.erase(endLocation, gsite.end());
      }
      // Eliminate duplicate values from the dissolution site vector
      dvsize = dsite.size();
      if (dvsize > 0) {
        sort(dsite.begin(), dsite.end());
        beginLocation = dsite.begin();
        endLocation = unique(dsite.begin(), dsite.end());
        dsite.erase(endLocation, dsite.end());
      }

      gvsize = gsite.size();
      for (int j = 0; j < gvsize; j++ ){
        stId = gsite[j]->getId();
        site_[stId].setInGrowInterfacePos(i,j);
      }
      dvsize = dsite.size();
      for (int j = 0; j < dvsize; j++ ){
        stId = dsite[j]->getId();
        site_[stId].setInDissInterfacePos(j);
      }

      growthInterfaceSize_[i] = gvsize;
      dissolutionInterfaceSize_[i] = dvsize;

      interface_.push_back(Interface(chemSys_, gsite, dsite, i, verbose_));

    } else { // Not a solid phase (either electrolyte or void)

      interface_.push_back(Interface(verbose_));
    }
  }

  cout
      << endl
      << "                         ***** initial count_ & interface sizes *****"
      << endl;
  cout << "***   numMicroPhases = " << numMicroPhases_ << endl;

  for (int i = 0; i < numMicroPhases_; i++) {
    cout << "  " << setw(3) << i << " : " << setw(15) << left << chemSys_->getMicroPhaseName(i)
         << setw(4) << right << " id:" << setw(3) << chemSys_->getMicroPhaseId(i)
         << "     count_ = " << setw(8) << count_[i] << "     dissolutionInterfaceSize_ =  "
         << setw(8) << dissolutionInterfaceSize_[i]
         << "     growthInterfaceSize_ =  " << setw(8) << growthInterfaceSize_[i]
         << "     porosity : " << chemSys_->getMicroPhasePorosity(i)
         << "     templates : ";
    for (int j = 0; j < numMicroPhases_; j++) {
      if (chemSys_->isGrowthTemplate(i, j)) {
        cout << j << " ";
        //}else{
        //    cout <<"";
      }
    }
    cout << endl;
  }
  cout.flush();
/*
  {
    for (int i = 0; i < numMicroPhases_; i++) {
      int sizeIntLatt = dissolutionInterfaceSize_[i];
      if (sizeIntLatt > 0) {
        int sizeIntInte = interface_[i].getDissolutionSize();
        int intSiteId, intPhId, posSite;
        cout << endl
             << endl
             << "   Lattice::findInterfaces PhId/sizeIntLatt/sizeIntInte : " << i << " / "
             << sizeIntLatt << " / " << sizeIntInte << endl
             << endl;

        for (int jj = 0; jj < sizeIntLatt; jj++) {
          intSiteId = interface_[i].getDissolutionSitesId(jj);
          intPhId = site_[intSiteId].getMicroPhaseId();
          posSite = site_[intSiteId].getInDissInterfacePos();
          if (jj == posSite) {
            cout << "      " << jj << "     intSiteId : " << intSiteId
                 << "     posSite : " << posSite << "     intPhId : " << intPhId << endl;
          } else {
            cout << "      " << jj << "     intSiteId : " << intSiteId
                 << "     posSite : " << posSite << "     intPhId : " << intPhId << "   +++"
                 << endl;
          }
        }
      }
    }
  }
*/

  // cout << "***   findInterfaces exit" << endl; exit(0);

  return;
}

vector<int> Lattice::growPhase(vector<int> growPhaseIDVect,
                               vector<int> numSiteGrowVect,
                               vector<string> growPhNameVect, int &numadded_G,
                               int totalTRC) {

  //*** for controll
  int bcl = 0;
  int static trc_g;
  trc_g++;

  int i, j, jj;

  int growPhaseIDVectSize = growPhaseIDVect.size();
  vector<int> numChange(growPhaseIDVectSize, 0);
  vector<int> dim_isite(growPhaseIDVectSize, 0);
  vector<int> nucleated(growPhaseIDVectSize, 0);
  vector<bool> writeFirst(growPhaseIDVectSize, false);
  vector<int> numLeft = numSiteGrowVect; // numtotake
  vector<int> inGrowInterfacePos;
  vector<Isite> isite;
  vector<unsigned int> plist;
  int posGrPhId, plistSize;
  int phaseInGrowVect;

  bool needUpdate = false;
  bool needNucleation = true;

  double wmcIni, wmcEnd;
  double steWmc, stenbWmc, dwmcval;
  double rng, probRNG, afty;

  unsigned int pid;
  int stId, phId, pos, posGrowVect;
  int isitePos;
  int posVect;
  int siteID, phaseID;

  Site *ste, *stenb;

  string nameMicroPhaseTh;

  //int totInterfaceSize = 0;
  int numChangeTot = 0;
  int numLeftTot = 0;
  for (i = 0; i < growPhaseIDVectSize; i++) {
    numLeftTot += numSiteGrowVect[i];
    //totInterfaceSize += growthInterfaceSize_[growPhaseIDVect[i]];
  }
  if (numLeftTot == 0) {
    cout << endl << "     Lattice::growPhase error => numLeftTot = 0" << endl;
    cout << "     totalTRC/trc_g/bcl :  "
         << "   " << totalTRC << "/" << trc_g << "/" << bcl << endl;
    cout << endl << "     stop program" << endl;
    exit(0);
  }

  // growth probabilities based on affinities
  vector<structGrowVect> growthVector;
  structGrowVect growStruct;

  int posProbVect = 0;
  double affSum = 0;
  for (i = 0; i < growPhaseIDVectSize; i++) {
    phaseID = growPhaseIDVect[i];
    isite = interface_[phaseID].getGrowthSites();
    dim_isite[i] = isite.size();
    for (jj = 0; jj < dim_isite[i]; jj++) {
    //for (jj = 0; jj < 1; jj++) {
      siteID = isite[jj].getId();
      afty = isite[jj].getAffinity();
      affSum += afty;
      growStruct.id = siteID;
      growStruct.posVect = i;
      growStruct.affinity = afty;
      growthVector.push_back(growStruct);
      //growProbStruct(int id_i = 0, int affinity_i = 0, int posVect_i = 0){}
      site_[siteID].setInGrowthVectorPos(phaseID, posProbVect);
      posProbVect++;
    }
  }
  int growthVectorSize = growthVector.size();

  cout << endl
       << "    Lattice::growPhase GROW_INI totalTRC/trc_g/bcl/affSum " << totalTRC
       << "/" << trc_g << "/" << bcl << "/" << affSum << endl;
  cout << "      GROW_INI growPhaseIDVectSize = " << growPhaseIDVectSize
       << "   growthVectorSize = " << growthVectorSize
       << "   numLeftTot = " << numLeftTot
       << "   numChangeTot = " << numChangeTot << endl;
  for (i = 0; i < growPhaseIDVectSize; i++) {
    phaseID = growPhaseIDVect[i];
    cout << "        GROW_INI for i = " << setw(3) << i
         << "  => phaseID phaseName count_ dim_isite numleft numchange  :  "
         << setw(3) << phaseID << "   " << setw(15) << left << growPhNameVect[i] << "   "
         << setw(8) << right << count_[phaseID]
         << "   " << setw(8) << dim_isite[i] << "   " << setw(8) << numLeft[i] << "   "
         << setw(8) << numChange[i] << endl;
  }
  cout.flush();

  while ((numLeftTot > 0) && (growthVectorSize == 0)) {
    cout << endl
       << "    Lattice::growPhase growPhaseIDVectSize = 0" << endl;
    cout << "      => need an initial nucleation for at least one of the growing phases" << endl;
    cout << "      there is(are) " << growPhaseIDVectSize << " growing phase(s):" << endl;
    cout << "      ";
    for (int i = 0; i < growPhaseIDVectSize; i++) {
      cout << "  " << growPhaseIDVect[i];
    }
    int nucPhaseId;
    do {
      nucPhaseId = (int)(rng * growPhaseIDVectSize);
    } while(numLeft[nucPhaseId] == 0);
    cout << endl << "      one of them is chosen randomly for nucleation: " << growPhaseIDVect[nucPhaseId] << endl;

    phaseID = growPhaseIDVect[nucPhaseId];
    dim_isite[nucPhaseId] = growthInterfaceSize_[phaseID];
    if (dim_isite[nucPhaseId] == 0) { // nucleation
      needNucleation = true;
      cout << endl
           << "    *** Lattice::growPhase - need nucleation for phaseID = " << phaseID << endl;
      cout << "      interface dimension dim_site[" << nucPhaseId << "] = " << dim_isite[nucPhaseId]
           << " while numLeft[" << nucPhaseId << "] = " << numLeft[nucPhaseId] << endl;
      cout << "      => for this microPhase (" << growPhNameVect[nucPhaseId] << ") a number of "
           << numLeft[nucPhaseId] << " seed(s)/site(s) will be nucleated" << endl;

      //nucleatePhaseAff(phaseID, numLeft[i]);
      nucleatePhaseRnd(phaseID, numLeft[nucPhaseId]);
      nucleated[nucPhaseId] = numLeft[nucPhaseId];

      numLeftTot = numLeftTot - numLeft[nucPhaseId];
      numChangeTot = numChangeTot + numLeft[nucPhaseId];

      numChange[nucPhaseId] = numChange[nucPhaseId] + numLeft[nucPhaseId];
      numLeft[nucPhaseId] = 0;

      writeFirst[i] = true;
      cout << endl
           << "    Lattice::growPhase GROW_END BY NUCLEATION for i = " << nucPhaseId
           << "   totalTRC/trc_g/bcl " << totalTRC << "/" << trc_g << "/" << bcl << endl;
      cout << "      GROW_END growPhaseIDVectSize = " << growPhaseIDVectSize
           << "   growthVectorSize = " << growthVectorSize << "   numLeftTot = " << numLeftTot
           << "   numChangeTot = " << numChangeTot << endl;
      cout << "        GROW_END phaseid count_ dim_isite numleft numchange  :  "
           << setw(3) << growPhaseIDVect[nucPhaseId] << "   "
           << setw(8) << right << count_[phaseID] << "   "
           << setw(8) << interface_[phaseID].getGrowthSites().size() << "   "
           << setw(8) << numLeft[nucPhaseId] << "   "
           << setw(8) << numChange[nucPhaseId] << endl;
      cout.flush();

      posProbVect = 0;
      affSum = 0;
      for (int j = 0; j < growPhaseIDVectSize; j++) {
        if (numLeft[j] > 0) {
          phaseID = growPhaseIDVect[j];
          isite = interface_[phaseID].getGrowthSites();
          dim_isite[j] = isite.size();

          for (int jj = 0; jj < dim_isite[j]; jj++) {
            afty = isite[jj].getAffinity();
            affSum += afty;
            growStruct.id = isite[jj].getId();
            growStruct.posVect = j;
            growStruct.affinity = afty;
            growthVector.push_back(growStruct);
            site_[isite[jj].getId()].setInGrowthVectorPos(phaseID, posProbVect);
            posProbVect++;
          }
        }
      }
      growthVectorSize = growthVector.size();
    }
  }

  while ((numLeftTot > 0) && (growthVectorSize >= 1)) {
      bcl++;

      // growth probabilities based on affinities
      // affAllPos - works with modified affinities in affinity_ table from
      // ChemicalSystem (all positives!)

      //affSum = 0;
      //for (j = 0; j < growthVectorSize; j++) {
      //  aff = growthVector[j].getAffinity();
      //  affSum += aff;
      //}
      // calc probabilities & choose a site
      rng = callRNG();
      if (affSum > 0) {
        probRNG = growthVector[0].affinity /affSum;
        if (rng <= probRNG) {
          isitePos = 0;
        } else {
          for (isitePos = 1; isitePos < growthVectorSize; isitePos++) {
            probRNG += (growthVector[isitePos].affinity / affSum);
            if (rng <= probRNG)
              break;
          }
        }
      } else {
        isitePos = (int)(rng * growthVectorSize);
      }

      ste = &site_[growthVector[isitePos].id];
      pid = ste->getMicroPhaseId(); // always ELECTROLYTEID !!
      posVect = growthVector[isitePos].posVect;
      phaseID = growPhaseIDVect[posVect]; // phase to grow

      if (pid != ELECTROLYTEID) {
        cout << endl
             << "Lattice::growPhase error: phaseid != ELECTROLYTEID" << endl;
        cout << "  phaseid totalTRC/trc_g/bcl count_ numLeftTot numChangeTot  "
                ":  "
             << phaseID << "   " << totalTRC << "/" << trc_g << "/" << bcl
             << "   " << count_[phaseID] << "   " << numLeftTot << "   "
             << numChangeTot << endl;
        cout << " posVect pid : " << posVect << "   " << pid << endl;
        cout << " ste.id numLeft numChange : " << ste->getId() << "   "
             << numLeft[posVect] << "   " << numChange[posVect] << endl;
        cout << "STOP" << endl;
        exit(0);
      }

      wmcIni = ste->getWmc0();

      plist = ste->getGrowthPhases(); // include phaseID i.e. faza care creste pe acest voxel
      plistSize = plist.size();
      //cout << endl << "    bcl = " << bcl << "   id/pid/posVect/phaseID/plistSize : "
      //<< ste->getId() << " / " << pid << " / " << posVect << " / " << phaseID << " / " << plistSize << endl;
      for (int j = 0; j < plistSize; j++) {
        for (int k = 0; k < growPhaseIDVectSize; k++) {
          if (growPhaseIDVect[k] == plist[j] && numLeft[k] > 0) {
            posGrPhId = ste->getInGrowthVectorPos(plist[j]);
            affSum -= growthVector[posGrPhId].affinity;
            //extractFromgrowthVector.push_back(posGrPhId);
            if (posGrPhId != growthVectorSize - 1) {
              growthVector[posGrPhId] = growthVector[growthVectorSize - 1];
              pos = growthVector[posGrPhId].posVect;
              phaseInGrowVect = growPhaseIDVect[pos];
              site_[growthVector[posGrPhId].id].setInGrowthVectorPos(phaseInGrowVect,posGrPhId);
            }
            growthVector.pop_back();
            growthVectorSize--;
            ste->setInGrowthVectorPos(plist[j], -1);
            break;
          }
        }
        removeGrowthSite_grow(ste, plist[j]);
      }
      ste->clearGrowth();
      setMicroPhaseId(ste, phaseID);

      ///
      /// Weighted mean curvature (wmc) is changed by the difference
      /// between the growing phase's porosity and the template's porosity.
      ///
      /// @todo Determine why the calculation works this way.
      ///

      nameMicroPhaseTh = chemSys_->getMicroPhaseName(phaseID);
      if (nameMicroPhaseTh == "CSHQ") {
        rng = callRNG();
        if (rng >= thrPorosityCSH) {
          wmcEnd = chemSys_->getMicroPhasePorosity(phaseID);
        } else {
          wmcEnd = 0;
        }
      } else {
        wmcEnd = 0;
      }
      ste->setWmc0(wmcEnd);

      // dwmcval = chemSys_->getMicroPhasePorosity(phaseid) -
      //           chemSys_->getMicroPhasePorosity(pid);
      dwmcval = wmcEnd - wmcIni;
      ste->dWmc(dwmcval);

      ///
      /// Now that the site has been added, it is eligible for dissolution
      /// later on, so we add it to the list of dissolution sites.
      ///

      steWmc = ste->getWmc();

      if (steWmc > 0.0) {
        addDissolutionSite(ste, phaseID);
      }

      ///
      /// Update the wmc of the neighboring sites.  This can be done
      /// because the wmc is originally calculated within a box around each
      /// site, so any time the id of a site within that box changes, it
      /// will change the wmc of the site at the box's center.
      ///
      /// NN_NNN = NUM_NEAREST_NEIGHBORS +  NUM_SECONDNEAREST_NEIGHBORS;
      ///

      for (j = 0; j < NN_NNN; j++) {
        stenb = ste->nb(j);
        stenb->dWmc(dwmcval);
        stenbWmc = stenb->getWmc();
        if (stenb->getMicroPhaseId() == ELECTROLYTEID) {
          inGrowInterfacePos = stenb->getInGrowInterfacePosVector();
          for (int k = FIRST_SOLID; k < numMicroPhases_; k++) {
            pos = inGrowInterfacePos[k]; //instenb->getInGrowInterfacePos(k);
            if (pos != -1) {
              afty = chemSys_->getAffinity(k, phaseID);
              interface_[k].updateAffinity(pos, afty);
              posGrowVect = stenb->getInGrowthVectorPos(k);
              if (posGrowVect != -1) {
                affSum += afty;
                afty += growthVector[posGrowVect].affinity;
                growthVector[posGrowVect].affinity = afty;
              }

            } else {

              if (k == phaseID) {
                addGrowthSite(stenb, phaseID);

                if (numLeft[posVect] > 0) {
                  afty = interface_[k].getAffinity(growthInterfaceSize_[k] - 1);
                  affSum += afty;
                  growStruct.id = stenb->getId();
                  growStruct.posVect = posVect;
                  growStruct.affinity = afty;
                  growthVector.push_back(growStruct);
                  //growProbStruct(int id_i = 0, int affinity_i = 0, int posVect_i = 0){}
                  stenb->setInGrowthVectorPos(phaseID, growthVectorSize);
                  growthVectorSize++;
                }
              }
            }
          }
        } else if ((stenbWmc == 0.0) && (stenb->getMicroPhaseId() > ELECTROLYTEID)) {
          removeDissolutionSite(stenb, stenb->getMicroPhaseId());
        }
      }

      for (j = 0; j < NUM_NEAREST_NEIGHBORS; j++) {
        stenb = ste->nb(j);
        if (stenb->getMicroPhaseId() == ELECTROLYTEID) {
          for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases_; phaseTmpl++) {
            if (chemSys_->isGrowthTemplate(phaseTmpl, phaseID)) {
              if (stenb->getInGrowInterfacePos(phaseTmpl) == -1) {
                addGrowthSite(stenb, phaseTmpl);
                for (int kk = 0; kk < growPhaseIDVectSize; kk++){
                  if (growPhaseIDVect[kk] == phaseTmpl && numLeft[kk] > 0) {
                    afty = interface_[phaseTmpl].getAffinity(growthInterfaceSize_[phaseTmpl] - 1);
                    affSum += afty;
                    growStruct.id = stenb->getId();
                    growStruct.posVect = kk;
                    growStruct.affinity = afty;
                    growthVector.push_back(growStruct);
                    //growProbStruct(int id_i = 0, int affinity_i = 0, int posVect_i = 0){}
                    stenb->setInGrowthVectorPos(phaseTmpl,growthVectorSize);
                    growthVectorSize++;
                    break;
                  }
                }
              }
            }
          }
        }
      }

      numLeftTot--;
      numChangeTot++;

      numLeft[posVect]--;
      numChange[posVect]++;

      if (numLeft[posVect] < 0) {
        cout << endl
             << "error_numLeft totalTRC trc_g bcl numLeftTot numChangeTot  :  "
             << totalTRC << "   " << trc_g << "   " << bcl << "   "
             << numLeftTot << "   " << numChangeTot << endl;
        cout << "  phaseid totalTRC/trc_g/bcl count_ numLeftTot numChangeTot  "
                ":  "
             << phaseID << "   " << totalTRC << "/" << trc_g << "/" << bcl
             << "   " << count_[phaseID] << "   " << numLeftTot << "   "
             << numChangeTot << endl;
        cout << " posVect pid : " << posVect << "   " << pid << endl;
        cout << " ste.id numLeft numChange : " << ste->getId() << "   "
             << numLeft[posVect] << "   " << numChange[posVect] << endl;
        cout << "STOP" << endl;
        exit(0);
      }

      // growth probabilities based on affinities
      // test nucleatie - ini
      needUpdate = false;
      needNucleation = true;
      while (needNucleation) {
        needNucleation = false;
        for (i = 0; i < growPhaseIDVectSize; i++) {
          if (numLeft[i] > 0) {
            phaseID = growPhaseIDVect[i];
            dim_isite[i] = growthInterfaceSize_[phaseID];
            if (dim_isite[i] == 0) { // nucleation
              needNucleation = true;
              cout << endl << "    *** Lattice::growPhase - need nucleation for phaseID = " << phaseID << endl;
              cout << "      interface dimension dim_site[" << i << "] = " << dim_isite[i]
                      << " while numLeft[" << i << "] = " << numLeft[i] << endl;
              cout << "      => for this microPhase (" << growPhNameVect[i] << ") a number of "
                   << numLeft[i] << " seed(s)/site(s) will be nucleated" << endl;

              //nucleatePhaseAff(phaseID, numLeft[i]);
              nucleatePhaseRnd(phaseID, numLeft[i]);
              nucleated[i] = numLeft[i];

              numLeftTot = numLeftTot - numLeft[i];
              numChangeTot = numChangeTot + numLeft[i];

              numChange[i] = numChange[i] + numLeft[i];
              numLeft[i] = 0;

              writeFirst[i] = true;
              cout << endl << "    Lattice::growPhase GROW_END BY NUCLEATION for i = " << i
                   << "   totalTRC/trc_g/bcl " << totalTRC << "/" << trc_g << "/" << bcl << endl;
              cout << "      GROW_END growPhaseIDVectSize = " << growPhaseIDVectSize
                   << "   growthVectorSize = " << growthVectorSize << "   numLeftTot = " << numLeftTot
                   << "   numChangeTot = " << numChangeTot << endl;
              cout << "        GROW_END phaseid count_ dim_isite numleft numchange  :  "
                   << setw(3) << growPhaseIDVect[i] << "   "
                   << setw(8) << right << count_[growPhaseIDVect[i]] << "   "
                   << setw(8) << dim_isite[i] << "   "
                   << setw(8) << numLeft[i] << "   "
                   << setw(8) << numChange[i] << endl;
              cout.flush();

              needUpdate = true;
            }
          }
        }
      }
      // test nucleatie - fin

      for (i = 0; i < growPhaseIDVectSize; i++) {
        if (numLeft[i] == 0) {
          if (writeFirst[i]) {
            continue;
          } else {
            writeFirst[i] = true;
            phaseID = growPhaseIDVect[i];
            dim_isite[i] = growthInterfaceSize_[phaseID];
            cout << "    Lattice::growPhase GROW_END for i = " << i << "   totalTRC/trc_g/bcl "
                 << totalTRC << "/" << trc_g << "/" << bcl << endl;
            cout << "      GROW_END growPhaseIDVectSize = " << growPhaseIDVectSize
                 << "   growthVectorSize = " << growthVectorSize << "   numLeftTot = " << numLeftTot
                 << "   numChangeTot = " << numChangeTot << endl;     
            cout << "        GROW_END phaseid count_ dim_isite numleft numchange  :  "
               << setw(3) << growPhaseIDVect[i] << "   "
               << setw(8) << right << count_[growPhaseIDVect[i]] << "   "
               << setw(8) << dim_isite[i] << "   "
               << setw(8) << numLeft[i] << "   "
               << setw(8) << numChange[i] << endl;
            cout.flush();

            needUpdate = true;
          }
        }
      }

      if (needUpdate) {
        for (int i = 0; i < growthVectorSize; i++) {
          stId = growthVector[i].id;
          pid = growPhaseIDVect[growthVector[i].posVect];
          site_[stId].setInGrowthVectorPos(pid, -1);
          //site_[stId].resetInGrowthVectorPos();
        }
        growthVector.clear();
        posProbVect = 0;
        affSum = 0;
        for (int i = 0; i < growPhaseIDVectSize; i++) {
          if (numLeft[i] > 0) {
            phaseID = growPhaseIDVect[i];
            isite = interface_[phaseID].getGrowthSites();
            dim_isite[i] = isite.size();

            for (int jj = 0; jj < dim_isite[i]; jj++) {
              afty = isite[jj].getAffinity();
              affSum += afty;
              growStruct.id = isite[jj].getId();
              growStruct.posVect = i;
              growStruct.affinity = afty;
              growthVector.push_back(growStruct);
              //growProbStruct(int id_i = 0, int affinity_i = 0, int posVect_i = 0){}
              site_[isite[jj].getId()].setInGrowthVectorPos(phaseID, posProbVect);
              posProbVect++;
            }
          }
        }
        growthVectorSize = growthVector.size();
      }
  }

  numadded_G = numChangeTot;

  return (nucleated);

}

void Lattice::nucleatePhaseRnd(int phaseID,int numLeft) {

  int numLeftIni = numLeft;

  vector<int> watersites;
  vector<Site *> localNb;
  int sizeWS;
  double rng;
  int fSiteWS;
  int j, k;
  vector<int> seedID;

  int numSites = 0;
  for (int i = 0; i < numsites_; i++) {
    if (site_[i].getMicroPhaseId() == phaseID) numSites++;
  }
  int numSites0 = numSites;

  string namePhase = chemSys_->getMicroPhaseName(phaseID);
  cout << endl << "      Lattice::nucleatePhaseRnd INI for phaseID = " << phaseID
       << "   namePhase = " << namePhase << endl;
  cout << "        sites to add        : numLeft = " << numLeft << endl;
  cout << "        sites in the system : count_[" << phaseID << "] = "
       << count_[phaseID] << "   &   check numSites = " << numSites << endl;
  cout << "        pore sites          : count_[0] = " << count_[0] << endl;
  cout << "        electrolyte sites   : count_[1] = " << count_[1] << endl;
  cout << "        growthInterfaceSize_ = "
       << growthInterfaceSize_[phaseID] << endl;
  cout << "        dissolutionInterfaceSize_ = "
       << dissolutionInterfaceSize_[phaseID] << endl;

  /*
  cout << "       CHECK ERRORS BEFORE NUCLEATION:";
  if (numSites != 0 ||
      growthInterfaceSize_[phaseID] != 0 ||
      dissolutionInterfaceSize_[phaseID] != 0) {

    vector<Isite> diss = interface_[phaseID].getDissolutionSites();
    int dissSz = diss.size();
    double wmc, wmc0, wmcSum, wmcSum_t;
    string nameNbPhase;
    int nbID, nbPhID;
    bool findSite;
    int siteDiss = 0, siteBulk = 0;
    int numCSH, numELE;
    int error1 = 0, error2 = 0, error3 = 0,
         error4 = 0, error5 = 0;

    int numIsolated = 0;
    for (int i = 0; i < numsites_; i++) {
      if (site_[i].getMicroPhaseId() == phaseID) {
        //cout << endl << "   check possible errors for siteID/phaseID : "
        //     << i << " / " << phaseID << endl;
        numELE = 0;
        for (int k = 0; k < NN_NNN; k++) {
          nbPhID = site_[i].nb(k)->getMicroPhaseId();
          if (nbPhID == ELECTROLYTEID) {
            numELE++;
          }
        }
        if (numELE > 0) {
          //cout << "     error1 i/numELE : "
          //     << i << " / " << numELE << endl;
          error1++;
        }

        wmc = site_[i].getWmc();
        if (wmc > 0) {
          findSite = false;
          for (int j = 0; j < dissSz; j++) {
            if (i == diss[j].getId()) {
              findSite = true;
              if (j != dissSz - 1)
                diss[j] = diss[dissSz - 1];
              dissSz--;
              diss.pop_back();
              siteDiss++;
              break;
            }
          }
          if (!findSite) {
            //cout << "     error2 i/findSite : "
            //     << i << " / findSite = false" << endl;
            error2++;
          }
        }

        wmc0 = site_[i].getWmc0();
        if (namePhase != "CSHQ") {
          if (wmc0 > 0) {
            //cout << "     error3 i/namePhase/wmc0 : " << i
            //     << " / " << namePhase <<  " / " << wmc0 << endl;
            error3++;
          }
        }

        numCSH = 0;
        wmcSum = wmc0;
        for (int k = 0; k < NN_NNN; k++) {
          nbID = site_[i].nb(k)->getId();
          nbPhID = site_[i].nb(k)->getMicroPhaseId();
          if (chemSys_->getMicroPhaseName(nbPhID) == "CSHQ" &&
              site_[nbID].getWmc0() > 0) {
            numCSH++;
            wmcSum += site_[nbID].getWmc0();
          }
        }
        if (abs(wmc - wmcSum) > 1e-5) {
          //cout << "     error4 i/wmcSum/wmc : " << i
          //     << " / " << wmcSum <<  " / " << wmc << endl;
          error4++;
        }

        if (numCSH == 0 && wmc >0) {
          //cout << "     error5 i/namePhase/numCSH/wmc : " << i << " / "
          //     << namePhase << " / " << numCSH
          //     << " / " << wmc << endl;
          error5++;
        }
        if (numCSH == 0) numIsolated++;
      }
    }

    if (error1 > 0 || error2 > 0 || error3 > 0 || error4 > 0 || error5 > 0) {
      cout << endl << "     error: numSites != 0 -> " << numSites
           << " != 0" << endl;
      cout << endl
           << "     error: growthInterfaceSize_ != 0 -> "
           << growthInterfaceSize_[phaseID] << " != 0" << endl;
      cout << endl
           << "     error: dissolutionInterfaceSize_ != 0 -> "
           << dissolutionInterfaceSize_[phaseID] << " != 0"
           << endl;
      cout << endl
           << "     error: numIsolated = " << numIsolated << " (?)" << endl;
      cout << endl
           << "     error1 = " << error1 << endl;
      cout << "     error2 = " << error2 << endl;
      cout << "     error3 = " << error3 << endl;
      cout << "     error4 = " << error4 << endl;
      cout << "     error5 = " << error5 << endl;

      cout << endl << "     exit" << endl;
      bool is_Error = false;
      throw MicrostructureException("Lattice", "nucleatePhaseRnd",
                                    "some errors -> one of error1 .. error5 ",
                                    is_Error);
      //exit(1);
    } else {
      cout << " numIsolated = " << numIsolated << " => NO ERRORS!" << endl;
    }
  } else {
    cout << " NO ERRORS!" << endl;
  }
  */

  seedID.clear();
  if (numLeft <= count_[ELECTROLYTEID]) {
    watersites.clear();
    for (k = 0; k < numsites_; ++k) {
      if (site_[k].getMicroPhaseId() == ELECTROLYTEID) {
        watersites.push_back(k);
      }
    }
    sizeWS = watersites.size();
    cout << "      Lattice::nucleatePhaseRnd -> sizeWS: "
         << sizeWS  << endl;

    while (numLeftIni > 0) {
      rng = callRNG();
      fSiteWS = (int)(rng * sizeWS);

      seedID.push_back(watersites[fSiteWS]);

      numLeftIni--;

      watersites[fSiteWS] = watersites[sizeWS - 1];
      watersites.pop_back();
      sizeWS--;
    }
  } else {
    cout << endl
         << "     Lattice::nucleatePhaseRnd => requested nucleation "
            "numLeft is larger than the electrolyte voxel number!"
         << endl;
    cout << "     numLeft > count_[ELECTROLYTEID] : " << numLeft << " > "
         << count_[ELECTROLYTEID] << endl;
    cout
        << endl
        << "     There is no room to nucleate => normal exit of the program"
        << endl;
    bool is_Error = false;
    throw MicrostructureException("Lattice", "nucleatePhaseRnd",
                                  "no room for nucleation", is_Error);
    // exit(0);
  }

  int sizeSeedID = seedID.size();
  int sID, pid;
  Site *ste, *stenb;
  string nameMicroPhaseTh;
  double dwmcval, steWmc, stenbWmc;
  double wmcIni, wmcEnd;

  for (int i = 0; i < sizeSeedID; i++) {
    sID = seedID[i];
    ste = &site_[sID];
    pid = ste->getMicroPhaseId(); // always ELECTROLYTEID !!

    wmcIni = ste->getWmc0();

    removeGrowthSite_nucleation(ste, phaseID);
    setMicroPhaseId(ste, phaseID);

      ///
      /// Weighted mean curvature (wmc) is changed by the difference
      /// between the growing phase's porosity and the template's porosity.
      ///
      /// @todo Determine why the calculation works this way.
      ///

      nameMicroPhaseTh = chemSys_->getMicroPhaseName(phaseID);
      if (nameMicroPhaseTh == "CSHQ") {
        rng = callRNG();
        if (rng >= thrPorosityCSH) {
          wmcEnd = chemSys_->getMicroPhasePorosity(phaseID);
        } else {
          wmcEnd = 0;
        }
      } else {
        wmcEnd = 0;
      }
      ste->setWmc0(wmcEnd);

      dwmcval = wmcEnd - wmcIni;
      ste->dWmc(dwmcval);

      ///
      /// Now that the site has been added, it is eligible for dissolution
      /// later on, so we add it to the list of dissolution sites.
      ///

      steWmc = ste->getWmc();

      if (steWmc > 0.0) {
        addDissolutionSite(ste, phaseID);
      }

      ///
      /// Update the wmc of the neighboring sites.  This can be done
      /// because the wmc is originally calculated within a box around each
      /// site, so any time the id of a site within that box changes, it
      /// will change the wmc of the site at the box's center.
      ///
      /// NN_NNN = NUM_NEAREST_NEIGHBORS +  NUM_SECONDNEAREST_NEIGHBORS;
      ///

      for (j = 0; j < NN_NNN;j++) {
        stenb = ste->nb(j);
        stenb->dWmc(dwmcval);
        stenbWmc = stenb->getWmc();

        if (stenb->getMicroPhaseId() == ELECTROLYTEID) {
          //addGrowthSite(stenb, phaseID);
          if (stenb->getInGrowInterfacePos(phaseID) == -1) {
            //if (chemSys_->getAffinity(phaseID, phaseID) > 10) {//itself affinity!
              addGrowthSite(stenb, phaseID);
            //}
          } else {
            int pos = stenb->getInGrowInterfacePos(phaseID);
            if(interface_[phaseID].getGrowthSitesId(pos) != stenb->getId()){
              cout << endl
                   << "error Lattice::nucleatePhaseRnd interface_[phaseID]."
                      "getGrowthSitesId(pos) != stenb->getId() : "
                      "ste pid steWmc phaseID stenb"
                      " stenb_pid stenbWmc j : "
                   << ste->getId() << "   " << pid << "   " << ste->getWmc()
                   << "   " << phaseID << "   " << stenb->getId() << "   "
                   << stenb->getMicroPhaseId() << "   " << stenbWmc << "   "
                   << j << endl;
              cout << endl
                   << "pos/stenb->getInGrowInterfacePos(phaseID) : " << pos
                   << " / " << stenb->getInGrowInterfacePos(phaseID) << endl;
              cout.flush();
              cout << endl << "stop program" << endl;
              exit(1);
             }
          }

          if (j < NUM_NEAREST_NEIGHBORS) {
            for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases_;
                 phaseTmpl++) {
              if (chemSys_->isGrowthTemplate(phaseTmpl, phaseID)) {
                // addGrowthSite(stenb, nbgrowthtemp[jj]);
                //addGrowthSite(stenb, phaseTmpl);
                if (stenb->getInGrowInterfacePos(phaseTmpl) == -1) {
                  addGrowthSite(stenb, phaseTmpl);
                } else {
                  int pos = stenb->getInGrowInterfacePos(phaseTmpl);
                  if (interface_[phaseTmpl].getGrowthSitesId(pos) !=
                      stenb->getId()) {
                    cout << endl
                         << "error Lattice:: growPhase interface_[phaseTmpl]."
                            "getGrowthSitesId(pos) != stenb->getId() : "
                            "ste pid steWmc phaseID stenb"
                            " stenb_pid stenbWmc j : "
                         << ste->getId() << "   " << pid << "   "
                         << ste->getWmc() << "   " << phaseID << "   "
                         << "   " << stenb->getId() << "   "
                         << stenb->getMicroPhaseId() << "   " << stenbWmc
                         << "   " << j << endl;
                    cout << endl
                         << "pos/phaseTmpl : " << pos << "   " << phaseID
                         << endl;
                    cout << endl
                         << "stenb->getInGrowInterfacePos(phaseTmpl) : "
                         << stenb->getInGrowInterfacePos(phaseTmpl) << endl;
                    cout.flush();
                    cout << endl << "stop program" << endl;
                    exit(1);
                  }
                }
              }
            }
          }
        } else if ((stenbWmc == 0.0) &&
                   (stenb->getMicroPhaseId() > ELECTROLYTEID)) {
          removeDissolutionSite(stenb, stenb->getMicroPhaseId());
        }
      }
  }
  numSites = 0;
  for (int i = 0; i < numsites_; i++) {
    if (site_[i].getMicroPhaseId() == phaseID) numSites++;
  }
  if (numSites != numSites0 + numLeft) {
    cout << endl << "     error: numSites != numSites0 + numLeft -> "
         << numSites << " != " << numSites0 + numLeft << endl;
    cout << endl << "     error: numSites0 = " << numSites0
         << "  &  numLeft = " << numLeft << endl;
    cout << endl << "     Lattice::nucleatePhaseRnd END for phaseID = "
         << phaseID << endl;

    cout << endl << "     exit" << endl;
    bool is_Error = false;
    throw MicrostructureException("Lattice", "nucleatePhaseRnd",
                                  "not all numLeft have been nucleated ",
                                  is_Error);
    //exit(1);
  }
  cout << "      Lattice::nucleatePhaseRnd END for phaseID = "
       << phaseID << endl;
  cout << "        numNewSites = " << numLeft << endl;
  cout << "        growthInterfaceSize_ = "
       << growthInterfaceSize_[phaseID] << endl;
  cout << "        dissolutionInterfaceSize_ = "
       << dissolutionInterfaceSize_[phaseID] << endl;
}

void Lattice::nucleatePhaseAff(int phaseID,int numLeft) {

  int numLeftIni = numLeft;
  struct localStruct {
    int id;
    double aff;
    double prob;
  };
  localStruct un;
  vector<localStruct> watersites;
  double aff, affSum;
  vector<Site *> localNb;
  int sizeWS;
  double rng;
  int fSiteWS;
  int j, k;
  vector<int> seedID;


  int numSites = 0;
  for (int i = 0; i < numsites_; i++) {
    if (site_[i].getMicroPhaseId() == phaseID) numSites++;
  }
  int numSites0 = numSites;

  string namePhase = chemSys_->getMicroPhaseName(phaseID);
  cout << endl << "      Lattice::nucleatePhaseAff INI for phaseID = " << phaseID
       << "   namePhase = " << namePhase << endl;
  cout << "        sites to add        : numLeft = " << numLeft << endl;
  cout << "        sites in the system : count_[" << phaseID << "] = "
       << count_[phaseID] << "   &   check numSites = " << numSites << endl;
  cout << "        pore sites          : count_[0] = " << count_[0] << endl;
  cout << "        electrolyte sites   : count_[1] = " << count_[1] << endl;
  cout << "        growthInterfaceSize_ = "
       << growthInterfaceSize_[phaseID] << endl;
  cout << "        dissolutionInterfaceSize_ = "
       << dissolutionInterfaceSize_[phaseID] << endl;

/*
  cout << "       CHECK ERRORS BEFORE NUCLEATION:";
  if (numSites != 0 ||
      growthInterfaceSize_[phaseID] != 0 ||
      dissolutionInterfaceSize_[phaseID] != 0) {

    vector<Isite> diss = interface_[phaseID].getDissolutionSites();
    int dissSz = diss.size();
    double wmc, wmc0, wmcSum;
    string namePhase = chemSys_->getMicroPhaseName(phaseID);
    string nameNbPhase;
    int nbID, nbPhID;
    bool findSite;
    int siteDiss = 0, siteBulk = 0;
    int numCSH, numELE;
    int error1 = 0, error2 = 0, error3 = 0,
         error4 = 0, error5 = 0;

    int numIsolated = 0;
    for (int i = 0; i < numsites_; i++) {
      if (site_[i].getMicroPhaseId() == phaseID) {
        //cout << endl << "   check possible errors for siteID/phaseID : "
        //     << i << " / " << phaseID << endl;
        numELE = 0;
        for (int k = 0; k < NN_NNN; k++) {
          nbPhID = site_[i].nb(k)->getMicroPhaseId();
          if (nbPhID == ELECTROLYTEID) {
            numELE++;
          }
        }
        if (numELE > 0) {
          //cout << "     error1 i/numELE : "
          //     << i << " / " << numELE << endl;
          error1++;
        }

        wmc = site_[i].getWmc();
        if (wmc > 0) {
          findSite = false;
          for (int j = 0; j < dissSz; j++) {
            if (i == diss[j].getId()) {
              findSite = true;
              if (j != dissSz - 1)
                diss[j] = diss[dissSz - 1];
              dissSz--;
              diss.pop_back();
              siteDiss++;
              break;
            }
          }
          if (!findSite) {
            //cout << "     error2 i/findSite : "
            //     << i << " / findSite = false" << endl;
            error2++;
          }
        }

        wmc0 = site_[i].getWmc0();

        if (namePhase != "CSHQ") {
          if (wmc0 > 0) {
            //cout << "     error3 i/namePhase/wmc0 : " << i
            //     << " / " << namePhase <<  " / " << wmc0 << endl;
            error3++;
          }
        }

        numCSH = 0;
        wmcSum = wmc0;
        for (int k = 0; k < NN_NNN; k++) {
          nbID = site_[i].nb(k)->getId();
          nbPhID = site_[i].nb(k)->getMicroPhaseId();
          if (chemSys_->getMicroPhaseName(nbPhID) == "CSHQ" &&
              site_[nbID].getWmc0() > 0) {
            numCSH++;
            wmcSum += site_[nbID].getWmc0();
          }
        }
        if (abs(wmc - wmcSum) > 1e-5) {
          //cout << "     error4 i/wmcSum/wmc : " << i
          //     << " / " << wmcSum <<  " / " << wmc << endl;
          error4++;
        }
        if (numCSH == 0 && wmc >0) {
          //cout << "     error5 i/namePhase/numCSH/wmc : " << i << " / "
          //     << namePhase << " / " << numCSH
          //     << " / " << wmc << endl;
          error5++;
        }
        if (numCSH == 0) numIsolated++;
      }
    }

    if (error1 > 0 || error2 > 0 || error3 > 0 || error4 > 0 || error5 > 0) {
      cout << endl << "     error: numSites != 0 -> " << numSites << " != 0"
           << endl;
      cout << endl
           << "     error: growthInterfaceSize_ != 0 -> "
           << growthInterfaceSize_[phaseID]
           << " != 0" << endl;
      cout << endl
           << "     error: dissolutionInterfaceSize_ != 0 -> "
           << dissolutionInterfaceSize_[phaseID] << " != 0" << endl;
      cout << endl
           << "     error: numIsolated = " << numIsolated << " (?)" << endl;
      cout << endl
           << "     error1 = " << error1 << endl;
      cout << "     error2 = " << error2 << endl;
      cout << "     error3 = " << error3 << endl;
      cout << "     error4 = " << error4 << endl;
      cout << "     error5 = " << error5 << endl;

      cout << endl << "     exit" << endl;
      bool is_Error = false;
      throw MicrostructureException("Lattice", "nucleatePhaseAff",
                                  "some errors -> one of error1 .. error5 ",
                                  is_Error);
      //exit(1);
    } else {
      cout << " numIsolated = " << numIsolated << " => NO ERRORS!" << endl;
    }
  } else {
    cout << " NO ERRORS!" << endl;
  }
*/

  seedID.clear();
  if (numLeft <= count_[ELECTROLYTEID]) {
    affSum = 0;
    watersites.clear();
    for (k = 0; k < numsites_; ++k) {
      if (site_[k].getMicroPhaseId() == ELECTROLYTEID) {
        un.id = site_[k].getId();
        aff = 0;
        localNb = site_[k].getNb();
        for (j = 0; j < NN_NNN; j++) {
          aff += chemSys_->getAffinity(phaseID, localNb[j]->getMicroPhaseId());
        }
        un.aff = aff;
        watersites.push_back(un);
        affSum += aff;
      }
    }
    sizeWS = watersites.size();
    cout << "      Lattice::nucleatePhaseAff -> sizeWS/affSum: "
        << sizeWS << " / " << affSum << endl;

    while (numLeftIni > 0) {

      rng = callRNG();

      if (affSum > 0) {
        watersites[0].prob = (watersites[0].aff) / affSum;
        for (k = 1; k < sizeWS; k++) {
          watersites[k].prob = watersites[k - 1].prob + (watersites[k].aff) / affSum;
        }

        for (fSiteWS = 0; fSiteWS < sizeWS; fSiteWS++) {
          if (rng <= watersites[fSiteWS].prob)
            break;
        }
      } else {

        fSiteWS = (int)(rng * sizeWS);
      }

      seedID.push_back(watersites[fSiteWS].id);

      numLeftIni--;

      watersites[fSiteWS] = watersites[sizeWS - 1];
      watersites.pop_back();
      sizeWS--;

      affSum = 0;
      for (k = 0; k < sizeWS; ++k) {
        aff = watersites[k].aff;
        affSum += aff;
      }
    }
  } else {
    cout << endl
         << "      Lattice::nucleatePhaseAff => requested nucleation "
            "numLeft is larger than the electrolyte voxel number!"
         << endl;
    cout << "        numLeft > count_[ELECTROLYTEID] : " << numLeft << " > "
         << count_[ELECTROLYTEID] << endl;
    cout
        << endl
        << "        There is no room to nucleate => normal exit of the program"
        << endl;
    bool is_Error = false;
    throw MicrostructureException("Lattice", "nucleatePhaseAff",
                                  "no room for nucleation", is_Error);
    // exit(0);
  }

  int sizeSeedID = seedID.size();
  int sID, pid;
  Site *ste, *stenb;
  string nameMicroPhaseTh;
  double dwmcval;
  double steWmc, stenbWmc;
  double wmcIni, wmcEnd;

  for (int i = 0; i < sizeSeedID; i++) {
    sID = seedID[i];
    ste = &site_[sID];
    pid = ste->getMicroPhaseId(); // always ELECTROLYTEID !!

    wmcIni = ste->getWmc0();

    removeGrowthSite_nucleation(ste, phaseID);
    setMicroPhaseId(ste, phaseID);

      ///
      /// Weighted mean curvature (wmc) is changed by the difference
      /// between the growing phase's porosity and the template's porosity.
      ///
      /// @todo Determine why the calculation works this way.
      ///

      nameMicroPhaseTh = chemSys_->getMicroPhaseName(phaseID);
      if (nameMicroPhaseTh == "CSHQ") {
        rng = callRNG();
        if (rng >= thrPorosityCSH) {
          wmcEnd = chemSys_->getMicroPhasePorosity(phaseID);
        } else {
          wmcEnd = 0;
        }
      } else {
        wmcEnd = 0;
      }
      ste->setWmc0(wmcEnd);

      dwmcval = wmcEnd - wmcIni;
      ste->dWmc(dwmcval);

      ///
      /// Now that the site has been added, it is eligible for dissolution
      /// later on, so we add it to the list of dissolution sites.
      ///

      steWmc = ste->getWmc();
      if (steWmc > 0.0) {
        addDissolutionSite(ste, phaseID);
      }

      ///
      /// Update the wmc of the neighboring sites.  This can be done
      /// because the wmc is originally calculated within a box around each
      /// site, so any time the id of a site within that box changes, it
      /// will change the wmc of the site at the box's center.
      ///
      /// NN_NNN = NUM_NEAREST_NEIGHBORS +  NUM_SECONDNEAREST_NEIGHBORS;
      ///

      for (j = 0; j < NN_NNN;j++) {
        stenb = ste->nb(j);
        stenb->dWmc(dwmcval);
        stenbWmc = stenb->getWmc();

        if (stenb->getMicroPhaseId() == ELECTROLYTEID) {
          //addGrowthSite(stenb, phaseID);
          if (stenb->getInGrowInterfacePos(phaseID) == -1) {
            addGrowthSite(stenb, phaseID);
          } else {
            int pos = stenb->getInGrowInterfacePos(phaseID);
            if(interface_[phaseID].getGrowthSitesId(pos) != stenb->getId()){
              cout << endl
                   << "error Lattice:: growPhase interface_[phaseID]."
                      "getGrowthSitesId(pos) != stenb->getId() : "
                      "ste pid steWmc phaseID stenb"
                      " stenb_pid stenbWmc j : "
                   << ste->getId() << "   " << pid << "   " << ste->getWmc()
                   << "   " << phaseID << "   " << stenb->getId() << "   "
                   << stenb->getMicroPhaseId() << "   " << stenbWmc << "   "
                   << j << endl;
              cout << endl
                   << "pos/stenb->getInGrowInterfacePos(phaseID) : " << pos
                   << " / " << stenb->getInGrowInterfacePos(phaseID) << endl;
              cout.flush();
              cout << endl << "stop program" << endl;
              exit(1);
             }
          }
          if (j < NUM_NEAREST_NEIGHBORS) {
            for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases_;
                 phaseTmpl++) {
              if (chemSys_->isGrowthTemplate(phaseTmpl, phaseID)) {
                // addGrowthSite(stenb, nbgrowthtemp[jj]);
                //addGrowthSite(stenb, phaseTmpl);
                if (stenb->getInGrowInterfacePos(phaseTmpl) == -1) {
                  addGrowthSite(stenb, phaseTmpl);
                } else {
                  int pos = stenb->getInGrowInterfacePos(phaseTmpl);
                  if (interface_[phaseTmpl].getGrowthSitesId(pos) !=
                      stenb->getId()) {
                    cout << endl
                         << "error Lattice:: growPhase interface_[phaseTmpl]."
                            "getGrowthSitesId(pos) != stenb->getId() : "
                            "ste pid steWmc phaseID stenb"
                            " stenb_pid stenbWmc j : "
                         << ste->getId() << "   " << pid << "   "
                         << ste->getWmc() << "   " << phaseID << "   "
                         << stenb->getId() << "   "
                         << stenb->getMicroPhaseId() << "   " << stenbWmc
                         << "   " << j << endl;
                    cout << endl
                         << "pos/phaseTmpl : " << pos << "   " << phaseID
                         << endl;
                    cout << endl
                         << "stenb->getInGrowInterfacePos(phaseTmpl) : "
                         << stenb->getInGrowInterfacePos(phaseTmpl) << endl;
                    cout.flush();
                    cout << endl << "stop program" << endl;
                    exit(1);
                  }
                }
              }
            }
          }
        } else if ((stenbWmc == 0.0) &&
                   (stenb->getMicroPhaseId() > ELECTROLYTEID)) {
          removeDissolutionSite(stenb, stenb->getMicroPhaseId());
        }
      }
  }

  numSites = 0;
  for (int i = 0; i < numsites_; i++) {
    if (site_[i].getMicroPhaseId() == phaseID) numSites++;
  }
  if (numSites != numSites0 + numLeft) {
    cout << endl << "     error: numSites != numSites0 + numLeft -> "
         << numSites << " != " << numSites0 + numLeft << endl;
    cout << endl << "     error: numSites0 = " << numSites0
         << "  &  numLeft = " << numLeft << endl;
    cout << endl << "     Lattice::nucleatePhaseAff END for phaseID = "
         << phaseID << endl;

    cout << endl << "     exit" << endl;
    bool is_Error = false;
    throw MicrostructureException("Lattice", "nucleatePhaseAff",
                                  "not all numLeft have been nucleated ",
                                  is_Error);
    //exit(1);
  }
  cout << "      Lattice::nucleatePhaseAff END for phaseID = "
       << phaseID << endl;
  cout << "        numNewSites = " << numLeft << endl;
  cout << "        growthInterfaceSize_ = "
       << growthInterfaceSize_[phaseID] << endl;
  cout << "        dissolutionInterfaceSize_ = "
       << dissolutionInterfaceSize_[phaseID] << endl;

}

int Lattice::dissolvePhase(vector<int> dissPhaseIDVect,
                           vector<int> numSiteDissVect,
                           vector<string> dissPhNameVect, int &numadded_D,
                           int totalTRC) {

  //*** controll
  int bcl = 0;
  int static trc_d;
  trc_d++;

  int i, ii, jj;

  Site *ste, *stenb;
  unsigned int pid;
  int isitePos, phaseID;
  int stId, nbid, nbpid;
  int posVect, posnb;
  double rng, probRNG, stWmc;
  double wmcIni, wmcEnd, dwmcval;
  bool phaseid_exist;

  int dissPhaseIDVectSize = dissPhaseIDVect.size();
  vector<int> numChange(dissPhaseIDVectSize, 0);
  vector<int> dim_isite(dissPhaseIDVectSize, 0);
  vector<bool> writeFirst(dissPhaseIDVectSize, false);
  vector<int> numLeft = numSiteDissVect; // numtotake
  vector<Isite> isite;
  vector<unsigned int> growth_local;
  int grLocSize;

  int numChangeTot = 0;
  int numLeftTot = 0;
  for (i = 0; i < dissPhaseIDVectSize; i++) {
    numLeftTot += numSiteDissVect[i];
  }
  if (numLeftTot == 0) {
    cout << "Lattice::dissolvePhase error numLeftTot = 0" << endl;
    cout << "   totalTRC/trc_d/bcl :  "
         << "   " << totalTRC << "/" << trc_d << "/" << bcl << endl;
    cout << "stop program" << endl;
    exit(0);
  }

  // dissolution probabilities based on wmc
  vector<structDissVect> dissolutionVector;
  structDissVect dissStruct;
  int posProbVect = 0;
  double sumWmc = 0;
  for (i = 0; i < dissPhaseIDVectSize; i++) {
    phaseID = dissPhaseIDVect[i];
    isite = interface_[phaseID].getDissolutionSites();
    dim_isite[i] = isite.size();
    for (jj = 0; jj < dim_isite[i]; jj++) {

      stId = isite[jj].getId();
      stWmc =  site_[stId].getWmc();
      sumWmc += stWmc;
      dissStruct.id = stId;
      dissStruct.posVect = i;
      dissStruct.wmc = stWmc;

      dissolutionVector.push_back(dissStruct);
      //dissProbStruct(int id_i = 0, double instab_i = 0, int posVect_i = 0, double prob_i = 0, double prob_01_i = 0) {}

      site_[stId].setInDissolutionVectorPos(posProbVect);
      posProbVect++;
    }
  }
  int dissolutionVectorSize = dissolutionVector.size();

  cout << endl
       << "    Lattice::dissolvePhase DISS_INI totalTRC/trc_d/bcl/sumWmc " << totalTRC
       << "/" << trc_d << "/" << bcl << "/" << sumWmc << endl;
  cout << "      DISS_INI dissPhaseIDVectSize = " << dissPhaseIDVectSize
       << "   dissolutionVectorSize = " << dissolutionVectorSize
       << "   numLeftTot = " << numLeftTot
       << "   numChangeTot = " << numChangeTot << endl;
  for (i = 0; i < dissPhaseIDVectSize; i++) {
    phaseID = dissPhaseIDVect[i];
    // isite = interface_[phaseID].getDissolutionSites();
    // dim_isite = isite.size();
    cout << "        DISS_INI for i = " << setw(3) << i
         << "  => phaseID phaseName count_ dim_isite numleft numchange  :  "
         << setw(3) << phaseID << "   " << setw(15) << left << dissPhNameVect[i] << "   "
         << setw(8) << right << count_[phaseID]
         << "   " << setw(8) << dim_isite[i] << "   " << setw(8) << numLeft[i] << "   "
         << setw(8) << numChange[i] << endl;
  }
  cout.flush();

  //int isitePosError = 0;

  int callGEM = -1;

/*
  if (trc_d == 10) {
    //dissolution
    for (int i = 0; i < numMicroPhases_; i++) {
      int sizeIntLatt = dissolutionInterfaceSize_[i];
      if (sizeIntLatt > 0) {
        int sizeIntInte = interface_[i].getDissolutionSize();
        int intSiteId, intPhId, posSite;
        cout << endl << endl
             << "   Lattice::dissolvePhase PhId/sizeIntLatt/sizeIntInte : " << i << " / "
             << sizeIntLatt << " / " << sizeIntInte << endl
             << endl;

        for (int jj = 0; jj < sizeIntLatt; jj++) {
          intSiteId = interface_[i].getDissolutionSitesId(jj);
          intPhId = site_[intSiteId].getMicroPhaseId();
          posSite = site_[intSiteId].getInDissInterfacePos();
          if (jj == posSite) {
            cout << "      " << jj << "     intSiteId : " << intSiteId
                 << "     posSite : " << posSite << "     intPhId : " << intPhId << endl;
          } else {
            cout << "      " << jj << "     intSiteId : " << intSiteId
                 << "     posSite : " << posSite << "     intPhId : " << intPhId << "   +++"
                 << endl;
          }
        }
      }
    }
    cout << endl << endl << "   *** Lattice::dissolvePhase - ini trc_d " << trc_d << " => exit " << endl;
    exit(0);
  }

  if (trc_d == 10) {
    //growth
    for (int i = 0; i < numMicroPhases_; i++) {
      int sizeIntLatt = growthInterfaceSize_[i];
      if (sizeIntLatt > 0) {
        int sizeIntInte = interface_[i].getGrowthSize();
        int intSiteId, intPhId, posSite;
        cout << endl << endl
             << "   Lattice::dissolvePhase - ini for growth trc_d/PhId/sizeIntLatt/sizeIntInte : " << trc_d <<  " / " << i << " / "
             << sizeIntLatt << " / " << sizeIntInte << endl
             << endl;

        for (int jj = 0; jj < sizeIntLatt; jj++) {
          intSiteId = interface_[i].getGrowthSitesId(jj);
          intPhId = site_[intSiteId].getMicroPhaseId();
          posSite = site_[intSiteId].getInGrowInterfacePos(i);
          if (jj == posSite) {
            cout << "      " << jj << "     intSiteId : " << intSiteId
                 << "     posSite : " << posSite << "     intPhId : " << intPhId << endl;
          } else {
            cout << "      " << jj << "     intSiteId : " << intSiteId
                 << "     posSite : " << posSite << "     intPhId : " << intPhId << "   +++"
                 << endl;
          }
        }
      }
    }
    cout << endl << endl << "   *** Lattice::dissolvePhase - ini for growth trc_d " << trc_d << " => exit " << endl;
    exit(0);
  }
*/
  while ((numLeftTot > 0) && (dissolutionVectorSize >= 1)) {
    try {
      bcl++;

      // dissolution probabilities based on wmc
      rng = callRNG();

      //sumWmc = 0.;
      //for (int i = 0; i < dissolutionVectorSize; i++) {
      //  sumWmc +=  dissolutionVector[i].instab;
      //}
      probRNG = dissolutionVector[0].wmc / sumWmc;
      if (rng <= probRNG) {
        isitePos = 0;
      } else {
        for (isitePos = 1; isitePos < dissolutionVectorSize; isitePos++) {
          probRNG +=  (dissolutionVector[isitePos].wmc / sumWmc);
          if (rng <= probRNG)
            break;
        }
      }

      //if (isitePos == dissolutionVectorSize) {
      //  isitePosError++;
      //  cout << endl << "     *** isitePosError: bcl/isitePos/rng/probRNG = " << bcl
      //       << " / " << isitePos << " / " << rng << " / " << probRNG << endl;
      //  cout << endl << "    exit" << endl; exit(1);
      //}

      ste = &site_[dissolutionVector[isitePos].id];
      pid = ste->getMicroPhaseId(); // intrebare pid diff phaseid ???
      posVect = dissolutionVector[isitePos].posVect;

      wmcIni = ste->getWmc0();
      sumWmc -= ste->getWmc();

      if (ste->getInDissInterfacePos() == -1) {
        cout << endl << "    Lattice::dissolvePhase error: ste->getInDissInterfacePos() = -1" << endl;
        cout << "    Lattice::dissolvePhase error: steId/pid/posVect  " << ste->getId()
               << "/" << pid << "/" << posVect << endl;
        cout << "    Lattice::dissolvePhase error: totalTRC/trc_d/bcl  " << totalTRC
               << "/" << trc_d << "/" << bcl << endl;
        cout << "    Lattice::dissolvePhase error: exit" << endl;
        exit(0);
      }
      removeDissolutionSite(ste, pid);
      setMicroPhaseId(ste, ELECTROLYTEID);

      /// Weighted mean curvature (wmc) is changed by the difference
      /// between the growing phase's porosity and the template's porosity.
      ///
      /// @todo Determine why the calculation works this way.
      ///

      wmcEnd = ELECTROLYTEID;
      ste->setWmc0(wmcEnd);

      dwmcval = wmcEnd - wmcIni;
      ste->dWmc(dwmcval);

      for (i = 0; i < NN_NNN; i++) {
        stenb = ste->nb(i);
        stenb->dWmc(dwmcval);
        nbpid = stenb->getMicroPhaseId();

        if (nbpid > ELECTROLYTEID) {

          ///
          /// Now that the site has been dissolved, it is eligible for growth
          /// later on, so we add it to the list of growth sites for certain
          /// phases.
          ///

          //addDissolutionSite(stenb, nbpid);         
          if (stenb->getInDissInterfacePos() == -1) {
            addDissolutionSite(stenb, nbpid);

            if (stenb->getInDissolutionVectorPos() != -1) {
              //error
              cout << endl << "Lattice::dissolvePhase : stenb->getInDissolutionVectorPos() != -1" << endl;
              cout << endl << "id/pid : " << ste->getId() << " / " << pid << endl;
              cout << endl << "stenb->getInDissolutionVectorPos() = " << stenb->getInDissolutionVectorPos() << endl;
              cout << endl << "dissolutionVectorSize : " << dissolutionVectorSize << endl;
              cout << endl << "i/nbid/nbpid : " << i << " / " << nbid << " / " << nbpid << endl;
              cout << endl
                   << "totalTRC trc_d bcl numLeftTot numChangeTot  :  " << totalTRC
                   << "   " << trc_d << "   " << bcl << "   " << numLeftTot << "   "
                   << numChangeTot << endl;
              cout << endl << "exit" << endl;
              exit(1);
            } else {
              for (int k = 0; k < dissPhaseIDVectSize; k++) {
                if (dissPhaseIDVect[k] == nbpid && numLeft[k] > 0){
                  stWmc = stenb->getWmc();
                  sumWmc += stWmc;
                  stenb->setInDissolutionVectorPos(dissolutionVectorSize);
                  dissStruct.id = stenb->getId();
                  dissStruct.posVect = k;
                  dissStruct.wmc = stWmc;

                  dissolutionVector.push_back(dissStruct);
                  //dissProbStruct(int id_i = 0, double instab_i = 0, int posVect_i = 0, double prob_i = 0, double prob_01_i = 0) {}
                  dissolutionVectorSize++;
                  break;
                }
              }
            }
          } else {
            posnb = stenb->getInDissolutionVectorPos();
            if (posnb != -1) {
              sumWmc += dwmcval;
              stWmc = dwmcval + dissolutionVector[posnb].wmc;
              dissolutionVector[posnb].wmc = stWmc;
            }

            //for (int k = 0; k < dissPhaseIDVectSize; k++) {
            //  if (dissPhaseIDVect[k] == nbpid && numLeft[k] > 0){
            //    posnb = stenb->getInDissolutionVectorPos();
            //    sumWmc += dwmcval;
            //    stWmc = dwmcval + dissolutionVector[posnb].wmc;
            //  }
            //}
          }

          //addGrowthSite(ste, nbpid);
          if (ste->getInGrowInterfacePos(nbpid) == -1) {
            //if (chemSys_->getAffinity(nbpid, nbpid) > 10) {//itself affinity!
              addGrowthSite(ste, nbpid);
            //}
          }
          if (i < NUM_NEAREST_NEIGHBORS) {

            for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases_;
                 phaseTmpl++) {
              if (chemSys_->isGrowthTemplate(phaseTmpl, nbpid)) {
                //addGrowthSite(ste, phaseTmpl);
                if (ste->getInGrowInterfacePos(phaseTmpl) == -1) {
                  addGrowthSite(ste, phaseTmpl);
                }
              }
            }
          }
        } else if (nbpid == ELECTROLYTEID) {
          nbid = stenb->getId();
          growth_local = site_[nbid].getGrowthPhases();
          grLocSize = growth_local.size();
          for (ii = 0; ii < grLocSize; ii++) {
            phaseid_exist = false;
            for (jj = 0; jj < NN_NNN; jj++) {
              nbpid = site_[nbid].nb(jj)->getMicroPhaseId();
              if (growth_local[ii] == nbpid) {
                phaseid_exist = true;
                break;
              }
              if (jj < NUM_NEAREST_NEIGHBORS) {
                if (chemSys_->isGrowthTemplate(growth_local[ii], nbpid)) {
                  phaseid_exist = true;
                  break;
                }
              }
            }
            if (phaseid_exist == false) {
               removeGrowthSite_diss(stenb, growth_local[ii]);
            }
          }
        }
      }

      //update dissolutionVector & involved sites
      site_[dissolutionVector[isitePos].id].setInDissolutionVectorPos(-1);
      if (isitePos != dissolutionVectorSize - 1) {
        dissolutionVector[isitePos] = dissolutionVector[dissolutionVectorSize - 1];
        site_[dissolutionVector[isitePos].id].setInDissolutionVectorPos(isitePos);
      }
      dissolutionVector.pop_back();
      dissolutionVectorSize--;

      numLeftTot--;
      numChangeTot++;

      numLeft[posVect]--;
      numChange[posVect]++;

    } catch (out_of_range &oor) {
      EOBException ex("Lattice", "dissolvePhase", "site_", site_.size(), i);
      ex.printException();

      cout << endl << "Lattice::dissolvePhase error" << endl;
      cout << endl
           << "totalTRC trc_g bcl numLeftTot numChangeTot  :  " << totalTRC
           << "   " << trc_d << "   " << bcl << "   " << numLeftTot << "   "
           << numChangeTot << endl;
      cout << endl
           << "steId pid dissolutionVectorSize :  " << ste->getId() << "   " << pid
           << "   " << dissolutionVectorSize << endl;
      cout.flush();
      for (i = 0; i < dissPhaseIDVectSize; i++) {
        phaseID = dissPhaseIDVect[i];
        isite = interface_[phaseID].getDissolutionSites();
        dim_isite[i] = isite.size();
        cout << "        phaseid count_ dim_isite numleft numchange  :  "
             << phaseID << "   " << "   " << count_[phaseID] << "   "
             << dim_isite[i] << "   " << numLeft[i] << "   " << numChange[i]
             << endl;
      }
      cout << "stop program" << endl;
      exit(1);
    }

    // dissolution probabilities based on wmc
    for (int i = 0; i < dissPhaseIDVectSize; i++) {
      if (numLeft[i] > 0) {
        phaseID = dissPhaseIDVect[i];
        dim_isite[i] = dissolutionInterfaceSize_[phaseID];
        if (dim_isite[i] == 0) {
          cout << endl << "    Lattice::dissolvePhase - anormal end:  dimInterface[" << i
               << "] = 0" << endl;
          cout << "      anormal end: phaseID[" << i << "] = " << phaseID
               << "  &  dim_site[" << i << "] = " << dim_isite[i]
               << " while numLeft[" << i << "] = " << numLeft[i] << endl;
          cout << "      anormal end: totalTRC/trc_d/bcl " << totalTRC << "/"
               << trc_d << "/" << bcl << endl;
          cout << "      need to call GEM after modifying DCLowerLimit_ for "
                  "microPhase "
               << phaseID << endl;
          cout.flush();

          callGEM = i;
          break;
        }

      } else {
        if (writeFirst[i]) {
          continue;
        } else {
          writeFirst[i] = true;

          phaseID = dissPhaseIDVect[i];
          dim_isite[i] = dissolutionInterfaceSize_[phaseID];
          cout << "    Lattice::dissolvePhase DISS_END for i = " << i
               << "   totalTRC/trc_d/bcl " << totalTRC << "/" << trc_d << "/"
               << bcl << endl;
          cout << "      DISS_END dissPhaseIDVectSize = " << dissPhaseIDVectSize
               << "   dissolutionVectorSize = " << dissolutionVectorSize
               << "   numLeftTot = " << numLeftTot
               << "   numChangeTot = " << numChangeTot << endl;
          cout << "        DISS_END phaseid count_ dim_isite numleft numchange  "
                  ":  "
               << setw(3) << dissPhaseIDVect[i] << "   " << setw(8) << right << count_[dissPhaseIDVect[i]]
               << "   " << setw(8) << dim_isite[i] << "   "
               << setw(8) << numLeft[i] << "   " << setw(8) << numChange[i] << endl;
          cout.flush();

          dissolutionVectorSize = dissolutionVector.size();
          for (int i = 0; i < dissolutionVectorSize; i++) {
            site_[dissolutionVector[i].id].setInDissolutionVectorPos(-1);
          }
          dissolutionVector.clear();
          sumWmc = 0.;
          posProbVect = 0;
          for (i = 0; i < dissPhaseIDVectSize; i++) {
            if (numLeft[i] > 0) {
              phaseID = dissPhaseIDVect[i];
              isite = interface_[phaseID].getDissolutionSites();
              dim_isite[i] = isite.size();
              for (jj = 0; jj < dim_isite[i]; jj++) {

                stId = isite[jj].getId();
                stWmc = site_[stId].getWmc();
                sumWmc += stWmc;
                dissStruct.id = stId;
                dissStruct.posVect = i;
                dissStruct.wmc =stWmc;

                dissolutionVector.push_back(dissStruct);
                //dissProbStruct(int id_i = 0, double instab_i = 0, int posVect_i = 0, double prob_i = 0, double prob_01_i = 0) {}

                site_[stId].setInDissolutionVectorPos(posProbVect);
                posProbVect++;
              }
            }
          }
          dissolutionVectorSize = dissolutionVector.size();
        }
      }
    }

    if (callGEM > -1) {
      dissolutionVectorSize = dissolutionVector.size();
      for (int i = 0; i < dissolutionVectorSize; i++) {
        site_[dissolutionVector[i].id].setInDissolutionVectorPos(-1);
      }
      break;
    }

    //dissolutionVectorSize = dissolutionVector.size();
    //cout << "     *** out totalTRC trc_d bcl :  " << totalTRC
    //                 << "   " << trc_d << "   " << bcl << endl;
  }

  /*
  for (int i = 0; i < numsites_; i++) {
    if (site_[i].getInDissolutionVectorPos() != -1) {
      cout << endl << "     *** out error site_[i].getInDissolutionVectorPos() != -1 for i = "
           << i << "   => site_[i].getInDissolutionVectorPos() = " << site_[i].getInDissolutionVectorPos() << endl;
      cout << "     *** out totalTRC trc_d bcl :  " << totalTRC
           << "   " << trc_d << "   " << bcl << endl;
      cout << endl << "     exit" << endl; exit(0);
    }
  }
  */

  numadded_D = numChangeTot;
  return (callGEM);
}

void Lattice::addDissolutionSite(Site *ste, unsigned int pid) {
  //try {
    ste->setInDissInterfacePos(dissolutionInterfaceSize_[pid]);
    interface_[pid].addDissolutionSite(ste);
    dissolutionInterfaceSize_[pid]++;
    ste->clearGrowth();
  //} catch (out_of_range &oor) {
  //  cout << endl << "EOB Lattice::addDissolutionSite pid = " << pid << " => exit" << endl;
  //  exit(1);
  //}
}

void Lattice::addGrowthSite(Site *ste, unsigned int pid) {
  //try {
    interface_[pid].addGrowthSite(ste);
    ste->setInGrowInterfacePos(pid, growthInterfaceSize_[pid]);
    growthInterfaceSize_[pid]++;
    ste->addGrowthPhaseId(pid);
  //} catch (out_of_range &oor) {
  //  cout << endl << "EOB Lattice::addGrowthSite pid = " << pid << " => exit" << endl;
  //  exit(1);
  //}
}

void Lattice::removeDissolutionSite(Site *ste0, unsigned int pid) {
  //try {
    int pos1 = dissolutionInterfaceSize_[pid] - 1;
    int sid = interface_[pid].getDissolutionSitesId(pos1);
    Site *ste1 = &site_[sid];
    int pos0 = ste0->getInDissInterfacePos();
    interface_[pid].removeDissolutionSite(pos0, pos1);
    //if (pos0 != pos1)
    ste1->setInDissInterfacePos(pos0);
    ste0->setInDissInterfacePos(-1);
    dissolutionInterfaceSize_[pid]--;
  //} catch (out_of_range &oor) {
  //  cout << endl << "EOB Lattice::removeDissolutionSite pid = " << pid << " => exit" << endl;
  //  exit(1);
  //}
}

void Lattice::removeGrowthSite_nucleation(Site *ste0, unsigned int pid) {
  //try {
    int pos0, pos1, grPhId;
    Site *ste1;
    int sid;
    vector<unsigned int> plist = ste0->getGrowthPhases();
    int plsize = plist.size();
    for (int j = 0; j < plsize; j++) {
      grPhId = plist[j];
      pos1 = growthInterfaceSize_[grPhId] - 1;
      sid = interface_[grPhId].getGrowthSitesId(pos1);
      ste1 = &site_[sid];
      pos0 = ste0->getInGrowInterfacePos(grPhId);
      interface_[grPhId].removeGrowthSite(pos0, pos1);
      //if (pos0 != pos1)
      ste1->setInGrowInterfacePos(grPhId, pos0);
      ste0->setInGrowInterfacePos(grPhId, -1);

      growthInterfaceSize_[grPhId]--;
    }
    ste0->clearGrowth();
  //} catch (out_of_range &oor) {
  //  cout << endl << "EOB Lattice::removeGrowthSite_nucleation pid = " << pid << " => exit" << endl;
  //  exit(1);
  //}
}

void Lattice::removeGrowthSite_grow(Site *ste0, int grPhId) {
  //try {
    int pos0, pos1, sid;
    Site *ste1;

    pos0 = ste0->getInGrowInterfacePos(grPhId);
    pos1 = growthInterfaceSize_[grPhId] - 1;
    sid = interface_[grPhId].getGrowthSitesId(pos1);

    interface_[grPhId].removeGrowthSite(pos0, pos1);
    //if (pos0 != pos1)
    site_[sid].setInGrowInterfacePos(grPhId, pos0);
    ste0->setInGrowInterfacePos(grPhId, -1);

    growthInterfaceSize_[grPhId]--;
  //} catch (out_of_range &oor) {
  //  cout << endl << "EOB Lattice::removeGrowthSite_grow grPhId = " << grPhId << " => exit" << endl;
  //  exit(1);
  //}

}

void Lattice::removeGrowthSite_diss(Site *ste0, unsigned int pid) {
  //try {
    int pos1 = growthInterfaceSize_[pid] - 1;
    int sid = interface_[pid].getGrowthSitesId(pos1);
    Site *ste1 = &site_[sid];
    int pos0 = ste0->getInGrowInterfacePos(pid);
    interface_[pid].removeGrowthSite(pos0, pos1);
    //if (pos0 != pos1)
    ste1->setInGrowInterfacePos(pid, pos0);
    ste0->setInGrowInterfacePos(pid, -1);

    ste0->removeGrowthSite(pid);

    growthInterfaceSize_[pid]--;
  //} catch (out_of_range &oor) {
  //  cout << endl << "EOB Lattice::removeGrowthSite_diss pid = " << pid << " => exit" << endl;
  //  exit(1);
  //}
}

int Lattice::emptyPorosity(int numsites, int cyc) {
  int maxsearchsize = 3;
  Site *stenb;

  ///
  /// Finding all potential VOID sites.
  ///
  /// @todo Consider removing some of the standard output, or setting a flag for
  /// it.
  ///

  cout << "    Lattice::emptyPorosity - check for cyc = " << cyc
       << " :      numsites = " << numsites;
  cout.flush();

  vector<int> distVect =
      findDomainSizeDistribution(ELECTROLYTEID, numsites, maxsearchsize, cyc, 0);
  int distVectSize = distVect.size();

  // cout << endl << "Lattice::emptyPorosity - check for cyc = " << cyc
  //      << " :      numsites = " << numsites
  cout << "      distVect.size() = " << distVectSize << endl;
  cout.flush();

  ///
  /// We want to empty the sites with the largest pore count
  ///

  if (distVectSize < numsites) {
    cout << endl
         << "    Lattice::emptyPorosity - Ran out of water in the system for "
            "cyc = "
         << cyc << "   distVect.size() < numsites : " << distVectSize
         << " < " << numsites << endl;
    cout << "    Lattice::emptyPorosity -> normal end of the program" << endl;
    //exit(1);
    bool is_Error = false;
    throw MicrostructureException("Lattice", "emptyPorosity",
                                  "ran out of water in the system - normal end", is_Error);
  }

  int numemptied = 0;
  int siteID;
  vector<unsigned int> grVect;
  int dimGrVect;
  Site *ste0, *ste1;
  int sid, pid, pos0, pos1;
  for (int ii = 0; ii < distVectSize; ii++) {
    siteID = distVect[ii];
    ste0 = &site_[siteID];
    setMicroPhaseId(siteID, VOIDID);
    grVect = site_[siteID].getGrowthPhases();
    dimGrVect = grVect.size();
    for (int i = 0; i < dimGrVect; i++) {
      pid = grVect[i];
      pos0 = ste0->getInGrowInterfacePos(pid);
      pos1 = growthInterfaceSize_[pid] - 1;
      sid = interface_[pid].getGrowthSitesId(pos1);
      ste1 = &site_[sid];
      interface_[pid].removeEmptiedSite(pos0, pos1);
      //if (pos0 != pos1)
      ste1->setInGrowInterfacePos(pid, pos0);
      ste0->setInGrowInterfacePos(pid, -1);
      growthInterfaceSize_[pid]--;
    }
    site_[siteID].setWmc0(0);
    site_[siteID].dWmc(-1);
    for (int i = 0; i < NN_NNN; i++) {
      stenb = site_[siteID].nb(i);
      stenb->dWmc(-1);
      if ((stenb->getWmc() == 0) &&
          (stenb->getMicroPhaseId() > ELECTROLYTEID))
        removeDissolutionSite(stenb, stenb->getMicroPhaseId());
    }
    site_[siteID].clearGrowth();
    numemptied++;
  }

  return (numemptied);
}

int Lattice::fillPorosity(int numsites, int cyc) {

  int maxsearchsize = 10;

  ///
  /// Finding all potential VOID sites.
  ///
  /// @todo Consider removing some of the standard output, or setting a flag for
  /// it.
  ///

  if (verbose_) {
    cout << "Lattice::fillPorosity Finding and sorting all "
         << "potential water sites ... ";
    cout.flush();
  }

  //cout << "  Lattice::fillPorosity - check for cyc = " << cyc
  //     << " :      numsites = " << numsites;
  //cout.flush();

  cout << "    Lattice::fillPorosity - check for cyc = " << cyc
       << " :      numsites = " << numsites;
  cout.flush();

  vector<int> distVect =
      findDomainSizeDistribution(VOIDID, numsites, maxsearchsize, cyc, 1);
  int distVectSize = distVect.size();

  cout << "      distVect.size() = " << distVectSize << endl;
  cout.flush();

  ///
  /// We want to fill the sites with the smallest pore count
  ///

  if (distVectSize < numsites) {
    cout << endl
         << "    Lattice::fillPorosity - not enough voids in the system for "
            "cyc = "
         << cyc << "   distVect.size() < numsites : " << distVectSize
         << " < " << numsites << endl;
    //exit(1);
    cout << "    Lattice::fillPorosity -> normal end of the program" << endl;
    bool is_Error = false;
    throw MicrostructureException("Lattice", "fillPorosity",
                                  "not enough voidsr in the system - normal end", is_Error);
  }
  int numfilled = 0;
  int siteID;
  int nbpid, nbid;
  Site *stenb;
  for (int ii = 0; ii < distVectSize; ii++) {
    siteID = distVect[ii];
    setMicroPhaseId(siteID, ELECTROLYTEID);
    site_[siteID].setWmc0(1);
    site_[siteID].dWmc(1);
    site_[siteID].clearGrowth();
    for (int i = 0; i < NN_NNN; i++) {
      stenb = site_[siteID].nb(i);
      // wmcIni = stenb->getWmc();
      stenb->dWmc(1);
      nbpid = stenb->getMicroPhaseId();
      if (nbpid > ELECTROLYTEID) {
        if (stenb->getInDissInterfacePos() == -1) addDissolutionSite(stenb, nbpid);
        if (site_[siteID].getInGrowInterfacePos(nbpid) == -1) addGrowthSite(&site_[siteID], nbpid);
        if (i < NUM_NEAREST_NEIGHBORS) {
          for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases_; phaseTmpl++) {
            if (chemSys_->isGrowthTemplate(phaseTmpl, nbpid)) {
              if (site_[siteID].getInGrowInterfacePos(phaseTmpl) == -1) {
                addGrowthSite(&site_[siteID], phaseTmpl);
              }
            }
          }
        }
      }
    }
    numfilled++;
  }

  return (numfilled);
}

int Lattice::countBox(int boxsize, unsigned int siteid) {
  //string msg;
  int boxhalf = boxsize / 2;
  int nfound = 0;
  int ix, iy, iz, hx, hy, hz;

    int qxlo = site_[siteid].getX() - boxhalf;
    int qxhi = site_[siteid].getX() + boxhalf;
    int qylo = site_[siteid].getY() - boxhalf;
    int qyhi = site_[siteid].getY() + boxhalf;
    int qzlo = site_[siteid].getZ() - boxhalf;
    int qzhi = site_[siteid].getZ() + boxhalf;

    ///
    /// Make sure that the box bounds are all sensible given the
    /// boundary conditions
    ///

    qxlo += checkBC(qxlo, xdim_);
    qxhi += checkBC(qxhi, xdim_);
    qylo += checkBC(qylo, ydim_);
    qyhi += checkBC(qyhi, ydim_);
    qzlo += checkBC(qzlo, zdim_);
    qzhi += checkBC(qzhi, zdim_);

    ///
    /// Count the number of water or void sites in the box
    ///

    for (iz = qzlo; iz <= qzhi; iz++) {
      hz = iz + checkBC(iz, zdim_);
      for (iy = qylo; iy <= qyhi; iy++) {
        hy = iy + checkBC(iy, ydim_);
        for (ix = qxlo; ix <= qxhi; ix++) {
          hx = ix + checkBC(ix, xdim_);
          if (site_[getIndex(hx, hy, hz)].getMicroPhaseId() ==
                  ELECTROLYTEID ||
              site_[getIndex(hx, hy, hz)].getMicroPhaseId() == VOIDID) {
            nfound++;
          }
        }
      }
    }
  return (nfound);
}

void Lattice::setResolution(const double res) {
  ///
  /// Make sure that resolution is a valid value
  ///

    string msg;
    if (res <= 0.001) {
      cout << endl;
      msg = "Lattice resolution <= 0.001";
      throw DataException("Lattice", "setResolution", msg);
    }

    if (verbose_) {
      cout << "Lattice::setResolution Changing lattice resolution from ";
      cout << resolution_ << " to " << res << endl;
      cout.flush();
    }
    resolution_ = res;
}

vector<unsigned int> Lattice::getNeighborhood(const unsigned int sitenum,
                                              const int size) {
  int xp, yp, zp;
  // double dist;

  vector<unsigned int> nh;

  unsigned int xc = site_[sitenum].getX();
  unsigned int yc = site_[sitenum].getY();
  unsigned int zc = site_[sitenum].getZ();

  if (size == 0) {
    nh.push_back(getIndex(xc, yc, zc));
    return nh;
  }

  for (int k = -size; k <= size; k++) {
    zp = zc + k;
    for (int j = -size; j <= size; j++) {
      yp = yc + j;
      for (int i = -size; i <= size; i++) {
        xp = xc + i;
        /*
        dist = static_cast<double>(((xc - xp) * (xc - xp)));
        dist += static_cast<double>(((yc - yp) * (yc - yp)));
        dist += static_cast<double>(((zc - zp) * (zc - zp)));

        if ((sqrt(dist) - 0.5) <= (double)size) {
            nh.push_back(getIndex(xp,yp,zp));
        }
        */
        nh.push_back(getIndex(xp, yp, zp));
      }
    }
  }

  return nh;
}

unsigned int Lattice::getIndex(int ix, int iy, int iz) const {
  if (ix < 0) {
    if (BC != 1) {
      ix += xdim_;
    } else {
      ix = 0;
    }
  } else if (ix >= (int)xdim_) {
    if (BC != 1) {
      ix -= (int)xdim_;
    } else {
      ix = (int)xdim_ - 1;
    }
  }

  if (iy < 0) {
    if (BC != 2) {
      iy += (int)ydim_;
    } else {
      iy = 0;
    }
  } else if (iy >= (int)ydim_) {
    if (BC != 2) {
      iy -= (int)ydim_;
    } else {
      iy = (int)ydim_ - 1;
    }
  }

  if (iz < 0) {
    if (BC != 3) {
      iz += (int)zdim_;
    } else {
      iz = 0;
    }
  } else if (iz >= (int)zdim_) {
    if (BC != 3) {
      iz -= (int)zdim_;
    } else {
      iz = (int)zdim_ - 1;
    }
  }

  return (unsigned int)(ix + (xdim_ * iy) + ((xdim_ * ydim_) * iz));
}

int Lattice::changeMicrostructure(double time, const int simtype,
                                  bool &capWater, int &numDiff, int &phDiff,
                                  string &nameDiff, int repeat, int cyc) {
  unsigned int i, ii;
  int numadded, numadded_actual;
  unsigned int tpid;
  int cursites, newsites, tnetsites;
  int wcursites, wnewsites;
  double td, tvol, tmass;
  int needRecallGEM;
  vector<double> vol_next, vfrac_next;
  vector<int> netsites(numMicroPhases_, 0);
  vector<unsigned int> pid(numMicroPhases_, 0);
  vector<string> phasenames;
  vector<string> phName;
  string phName_t;

  extern string CSHMicroName;
  extern string MonocarbMicroName;
  extern string MonosulfMicroName;
  extern string HydrotalcMicroName;
  extern string AFTMicroName;

  static int totalTRC, normalTRC, totalRepeat;
  totalTRC++;
  if (repeat == 0) {
    normalTRC++;
  } else {
    totalRepeat++;
  }
  cout << endl << "  Lattice::changeMicrostructure -ini- cyc = " << cyc << endl;

  //checkSite(8); cout << endl << " exit Lattice::changeMicrostructure " << endl;

  ///
  /// @todo This function is very large; consider breaking it into small pieces
  /// for ease of maintenance and readability.
  ///
  /// Assign time of this state to the lattice time_ variable
  ///

  time_ = time;

  ///
  /// The next block, if uncommented, reads prior expansion
  /// strain data from a file and loads it into the appropriate
  /// class members.
  ///

  ///
  /// Zero out the amount of water to add to the microstructure.
  ///
  /// @todo Find out why; this class variable is not used again.
  ///

  waterchange_ = 0.0;

  ///
  /// Next, determine how many sites of each phase need to be assigned
  /// The ChemicalSystem object calculates the new volume fractions, so we
  /// just load them up here from the ChemicalSystem object.
  ///

  vol_next = chemSys_->getMicroPhaseVolume();
  phasenames = chemSys_->getMicroPhaseName();
  phName = phasenames;
  int volNextSize = vol_next.size();

  if (verbose_) {
    cout << endl << "Lattice::changeMicrostructure Before adjustMicrostructureVolumes cyc = " << cyc
         << endl;
    for (int iii = 0; iii < volNextSize; ++iii) {
      cout << "Lattice::changeMicrostructure    Volume of " << phasenames[iii]
           << " = " << vol_next[iii] << " m3" << endl;
    }
    cout.flush();
  }

  adjustMicrostructureVolumes(vol_next, volNextSize, cyc);

  if (verbose_) {
    cout << endl << "Lattice::changeMicrostructure After adjustMicrostructureVolumes cyc = " << cyc
         << endl;
    for (int iii = 0; iii < volNextSize; ++iii) {
      cout << "Lattice::changeMicrostructure   Volume of " << phasenames[iii]
           << " = " << vol_next[iii] << " m3" << endl;
    }
  }

  vfrac_next = vol_next;

  adjustMicrostructureVolFracs(phasenames, vol_next, vfrac_next, volNextSize, cyc);

  ///
  /// Calculate number of sites of each phase in next state
  ///

  // Remember we added two extra slots onto the end of vfrac_next
  // to hold the volume fraction of capillary space and subvoxel pore space
  // But those two are not needed here anymore

  if (verbose_) {
    cout << "Lattice::changeMicrostructure Calculating volume "
         << "of each phase to be added... cyc = " << cyc << endl;
      for (i = 0; i < volNextSize; i++) {
        cout << "Lattice::changeMicrostructure ****Volume fraction["
             << phasenames[i]
             << "] in next state should be = " << vfrac_next[i];
        cout << ", or " << (int)((double)(numsites_ * vfrac_next[i]))
             << " sites" << endl;
      }
      cout << "Lattice::changeMicrostructure ****Volume fraction[capillary "
           << "pores] in next state"
           << " should be = " << capillaryporevolumefraction_ << endl;
      cout << "Lattice::changeMicrostructure ****Volume fraction[subvoxel "
              "pores] "
              "in next state"
           << " should be = " << subvoxelporevolumefraction_ << endl;
      cout.flush();
  }

  ///
  /// The next block is executed only if there will eventually be some
  /// sulfate attack during this simulation.
  ///

  if (simtype == SULFATE_ATTACK) {
    /*
    ///
    ///  Normalize to get volume fractions and compute number
    ///  of sites of each phase needed
    ///
    ///  @todo Find out why we need to do all of this just because
    ///  there will eventually be sulfate attack.  Why not wait until
    ///  sulfate attack simulation actually starts?
    ///

    /// Used to be a bunch of sulfate atack stuff here
    */

  } else {
    ///
    /// Sulfate attack will NEVER be done during this simulation.
    ///  Normalize to get volume fractions and compute number
    ///  of sites of each phase needed.
    ///

    cout << "  Lattice::changeMicrostructure  initial values for count_[0]/count_[1] : "
         << count_[0] << " / " << count_[1] << endl;
    for (i = FIRST_SOLID; i < volNextSize; i++) { // from i = 2 !!!
      cursites = count_[i];
      newsites = ceil(numsites_  * vfrac_next[i]);
      //newsites = ceil((numsites_ - count_[0])  * vfrac_next[i]);
      netsites[i] = newsites - cursites;
      pid[i] = i;
      phName[i] = phasenames[i];
      if (netsites[i] != 0) {
        cout << "  Lattice::changeMicrostructure ***for "
             << setw(15) << left << phasenames[i]
             << " => i: " << setw(3) << right << i
             << "   count_[i]: " << setw(8) << right << count_[i]
             << "   vfrac_next[i]: " << right << vfrac_next[i]
             << "   netsites in this state: " << setw(9) << right << netsites[i]
             << "   cursites: " << setw(9) << right << cursites
             << "   newsites: " << setw(9) << right << newsites
             << "   -> phaseID: " << setw(3) << right
             << chemSys_->getMicroPhaseId(phasenames[i]) << endl;
        cout.flush();
      }
    }
  }

  ///
  /// Sort netsites in ascending order, except we will handle
  /// void space differently.

  int netsitesSize_;
  netsitesSize_ = netsites.size();
  vector<int> dissPhaseIDVect, growPhaseIDVect;
  vector<int> numSiteDissVect, numSiteGrowVect;
  vector<string> dissPhNameVect, growPhNameVect;

  int numTotSitesToDissolve = 0, numTotSitesToGrow = 0;

    for (i = FIRST_SOLID; i < netsitesSize_; i++) {
      if (netsites[i] < 0) {
        dissPhaseIDVect.push_back(pid[i]);
        numSiteDissVect.push_back(-netsites[i]);
        dissPhNameVect.push_back(phasenames[pid[i]]);
        numTotSitesToDissolve += -netsites[i];
      } else if (netsites[i] > 0) {
        growPhaseIDVect.push_back(pid[i]);
        numSiteGrowVect.push_back(netsites[i]);
        growPhNameVect.push_back(phasenames[pid[i]]);
        numTotSitesToGrow += netsites[i];
      }
    }

  ///
  /// The next loop starts at FIRST_SOLID because we exclude void and water
  /// phases
  ///
  /// @todo Consider making the starting index more general
  ///

  cout << "  Lattice::changeMicrostructure netsites.size() = " << netsitesSize_
       << " -> normalTRC/repeat (totalTRC/normalTRC/totalRepeat): " << normalTRC
       << " / " << repeat << "   (" << totalTRC << " / " << normalTRC << " / "
       << totalRepeat << ")" << endl;
  cout << "  Lattice::changeMicrostructure number of phases to dissolve for cyc = " << cyc
  << " : " << setw(3) << right << dissPhaseIDVect.size() << endl;
  cout << "  Lattice::changeMicrostructure number of phases to grow for cyc = " << cyc
  << "     : " << setw(3) << right << growPhaseIDVect.size() << endl;

  unsigned int pid_;
  int numadded_D = 0, numadded_G = 0;
  vector<int> nucleatedPhases;
  //nucleatedPhases.clear();
  //nucleatedPhases.resize(growPhaseIDVect.size(), 0);

  try {
    if (dissPhaseIDVect.size() > 0) {
      needRecallGEM = -1;
      needRecallGEM = dissolvePhase(dissPhaseIDVect, numSiteDissVect,
                                    dissPhNameVect, numadded_D, totalTRC);
      if (needRecallGEM > -1) {
        phDiff = dissPhaseIDVect[needRecallGEM];
        nameDiff = dissPhNameVect[needRecallGEM];
        numDiff = count_[phDiff];

        cout << endl << "  Lattice::changeMicrostructure - anormal end" << endl;
        cout << "    phDiff,nameDiff,numDiff,count_[phDiff] : " << phDiff
             << " , " << nameDiff << " , " << numDiff << " , " << count_[phDiff]
             << endl;
        cout << "    numTotSitesToDissolve = " << numTotSitesToDissolve
             << "  while numTotSitesDissolved = " << numadded_D << endl;
        cout << "    => recall GEM after (re)setDCLowerLimit according to the "
                "system configuration (lattice)"
             << endl;

        return 0;
      }
    }

    if (growPhaseIDVect.size() > 0) {
      try {

        nucleatedPhases = growPhase(growPhaseIDVect, numSiteGrowVect,
                                    growPhNameVect, numadded_G, totalTRC);
      } catch (MicrostructureException mex) {
        throw mex;
      } catch (out_of_range &oor) {
        EOBException ex("Lattice", "changeMicrostructure", "after growth",
                        1, 0);
        ex.printException();
        cout << endl << "Lattice::changeMicrostructure -after growth excp- cyc = " << cyc << endl; cout.flush();
        exit(1);
      }
    }

    if ((numadded_D != numTotSitesToDissolve) ||
        (numadded_G != numTotSitesToGrow)) {
      cout << endl << "      Lattice::changeMicrostructure error: " << endl;
      cout << "(numadded_D != numTotSitesToDissolve) || (numadded_G != "
              "numTotSitesToGrow)"
           << endl;
      cout << endl
           << "numadded_D = " << numadded_D
           << "       numTotSitesToDissolve  = " << numTotSitesToDissolve
           << endl;
      cout << endl
           << "numadded_G = " << numadded_G
           << "       numTotSitesToGrow  = " << numTotSitesToGrow << endl;
      cout << endl
           << "   cyc / totalTRC / normalTRC : " << cyc << " / " << totalTRC
           << " / " << normalTRC << endl;
      cout << endl
           << "   repeat / totalRepeat : " << repeat << " / " << totalRepeat
           << endl;
      cout << endl << "   nucleatedPhases :" << endl;
      for (i = 0; i < growPhaseIDVect.size(); i++) {
        cout << " i: " << i << "   phID: " << growPhaseIDVect[i]
             << "   phName: " << growPhNameVect[i]
             << "   toGrow: " << numSiteGrowVect[i]
             << "   grownByNucl: " << nucleatedPhases[i] << endl;
      }
      cout << endl << "stop program" << endl;
      exit(1);
    }
  } catch (out_of_range &oor) {
    throw EOBException("Lattice", "changeMicrostructure", "pid", pid.size(), i);
  }

  cursites = count_[VOIDID];
  newsites = (int)((numsites_ * vfrac_next[VOIDID]) + 0.5);
  wcursites = count_[ELECTROLYTEID];
  wnewsites = wcursites - (newsites - cursites);

  int numEmptyFill = newsites - cursites;
  int numEmpty = 0, numFill = 0;
  // numEmpty = emptyPorosity(newsites - cursites, cyc);
  cout << endl
       << "  Lattice::changeMicrostructure - before emptyPorosity/fillPorosity - check for cyc = "
       << cyc << "   numEmptyFill = " << numEmptyFill
       << "   count_[0]/count_[1] : " << count_[0] << " / " << count_[1] << endl;
  cout.flush();
  if (numEmptyFill > 0) {
    try {
      numEmpty = emptyPorosity(numEmptyFill, cyc);
    } catch (MicrostructureException mex) {
      throw mex;
    }
    cout << "  Lattice::changeMicrostructure - after emptyPorosity - check for cyc = "
         << cyc << "   numEmptyFill = " << numEmptyFill
         << "   numEmpty = " << numEmpty << "   count_[0]/count_[1] : "
         << count_[0] << " / " << count_[1] << endl;cout.flush();
  } else if (numEmptyFill < 0) {
    try {
      numFill = fillPorosity(-numEmptyFill, cyc);
    } catch (MicrostructureException mex) {
      throw mex;
    }
    cout << "  Lattice::changeMicrostructure - after fillPorosity - check for cyc = "
         << cyc << "   numEmptyFill = " << numEmptyFill
         << "   numFill = " << numFill << "   count_[0]/count_[1] : "
         << count_[0] << " / " << count_[1] << endl;cout.flush();
  }

  if (verbose_) {
    cout << "Lattice::changeMicrostructure ***netsites["
         << phasenames[VOIDID]
         << "] in this state = " << (newsites - cursites)
         << "; cursites = " << cursites << " and newsites = " << newsites
         << endl;
    cout << "Lattice::changeMicrostructure ***netsites["
         << phasenames[ELECTROLYTEID]
         << "] in this state = " << (wnewsites - wcursites)
         << "; cursites = " << wcursites << " and newsites = " << wnewsites
         << endl
         << endl;

    // When creating void from water, we should
    // update the target volume fraction of water even though
    // it is not used in any further calculations at this point

    cout << "Lattice::changeMicrostructure Target CAPIILARY WATER "
         << "volume fraction IS " << vfrac_next[ELECTROLYTEID] << endl;

    // vfrac_next[ELECTROLYTEID] -= ((double)(newsites -
    // cursites)/(double)(numsites_)); cout << "But WILL BE " <<
    // vfrac_next[ELECTROLYTEID] << " after creating void space" << endl;

    cout << "Lattice::changeMicrostructure Number CAPILLARY VOXELS "
         << "actually emptied was:  " << numEmptyFill << endl;

    ///
    /// Report on target and actual mass fractions
    ///

    cout << "Lattice::changeMicrostructure "
         << "*******************************" << endl;
    cout.flush();
  }

    int totcount = 0;
    double numSites = numsites_;
    for (i = 0; i < volNextSize; i++) {
      volumefraction_[i] = count_[i] / numSites;
      totcount += count_[i];
      if (verbose_) {
        cout << "  Lattice::changeMicrostructure Phase " << i
             << " Target volume fraction was " << vfrac_next[i]
             << " and actual is " << volumefraction_[i] << ", and "
             << totcount << " of " << site_.size() << " sites claimed so far"
             << endl;
        cout.flush();
      }
    }
    if (totcount != numsites_) {
      cout << endl
           << "Lattice::changeMicrostructure error => totcount != numsites_ : "
           << totcount << " != " << numsites_ << endl;
      cout << "time,simtype,capWater : " << time << " , " << simtype << " , "
           << capWater << endl;
      cout << "repeat  : " << repeat << endl;
      cout << "cyc  : " << cyc << endl;
      cout << "stop program" << endl;
      exit(1);
    }

  ///  This is a local variable and the value is never used.
  ///
  ///  @todo Why not eliminate this line completely?
  ///

  // double surfa = getSurfaceArea(chemSys_->getMicroPhaseId(CSHMicroName));

  if (volumefraction_[ELECTROLYTEID] <= 0.0) {
    capWater = false;
    cout << endl
         << "  Lattice::changeMicrostructure :  capWater = " << capWater << endl;
    cout << "mPhId/vfrac_next_/volumefraction_/count_: " << endl;
    for (int i = 0; i < volNextSize; i++) {
      cout << "   " << i << "  " << vfrac_next[i] << "  " << volumefraction_[i]
           << "  " << count_[i] << endl;
    }
    cout << endl
         << "  time, simtype, capWater : " << time << " , " << simtype << " , "
         << capWater << endl;
    cout << "  repeat  : " << repeat << endl;
    cout << "  cyc  : " << cyc << endl;
    cout << "  no water in the system => normal end of the program" << endl;

    bool is_Error = false;
    throw MicrostructureException("Lattice", "changeMicrostructure",
                                  "no capilary water in the system", is_Error);
    //exit(0);
  }

  cout << "  Lattice::changeMicrostructure -normal end- cyc = " << cyc << endl;

  return 1;
}

void Lattice::adjustMicrostructureVolumes(vector<double> &vol, int volSize,
                                          int cyc) {
  int i = 0;
  // int volSize = vol.size();

#ifdef DEBUG
  cout << "Lattice::adjustMicrostructureVolumes" << endl;
  cout.flush();
#endif

  //try {

    // Find the total system volume according to GEMS, in m3
    // units.  The individual microstructure phase volumes
    // have already been adjusted for subvoxel porosity in the
    // ChemicalSystem::calculateState function

    // The lattice has a fixed volume given by a fixed number
    // of voxels multiplied by the volume per voxel

    // The capillary water volume is now calculated, during
    // which the subvoxel pore volume is calculated as well

    // calcSolidvolumewithpores;
    solidvolumewithpores_ = 0.0;
    for (i = 0; i < volSize; ++i) {
      if (i != ELECTROLYTEID && i != VOIDID)
        solidvolumewithpores_ += vol[i];
    }
    if (solidvolumewithpores_ <= 0.0)
      throw DataException("Lattice", "adjustMicrostructureVolumes",
                          "totvolume is NOT positive");

    // The current microstructure volume as predicted by GEMS
    // The initial microstructure volume is calculated the first
    // time that Lattice::changeMicrostructure is called
    // We currently insist remain the volume of the system, so if
    // the microstructure volume deviates from the system volume,
    // we add or subtract capillary porosity to keep them equal

    // calcSubvoxelporevolume;
    subvoxelporevolume_ = 0.0;
    for (i = 0; i < volSize; ++i) {
      if (i != ELECTROLYTEID && i != VOIDID) {
        subvoxelporevolume_ += (vol[i] * chemSys_->getMicroPhasePorosity(i));
      }
    }
    double solidvolume = solidvolumewithpores_ - subvoxelporevolume_;
    // microstructurevolume_ = chemSys_->getMicroVolume();
    nonsolidvolume_ = initialmicrostructurevolume_ - solidvolume; //
    // nonsolidvolume_ IS capillaryvoidvolume_ + capillarywatervolume_ +
    // subvoxelporevolume_ capillaryporevolume_ = nonsolidvolume_ -
    // subvoxelporevolume_;
    capillaryporevolume_ = initialmicrostructurevolume_ - solidvolumewithpores_;

    // calcCapillarywatervolume
    watervolume_ = vol[ELECTROLYTEID];
    // Up to this point, vol[ELECTROLYTEID] is the total volume
    // of pore solution in the system regardless of whether it is
    // capillary or subvoxel
    // First need to trim off any extra water that lies outside the system

    voidvolume_ = nonsolidvolume_ - watervolume_;

    if (voidvolume_ < 0.0)
      watervolume_ = nonsolidvolume_;

    capillarywatervolume_ = watervolume_ - subvoxelporevolume_;
    if (capillarywatervolume_ < 0.0) {
      subvoxelwatervolume_ = watervolume_;
      capillarywatervolume_ = 0;
    } else {
      subvoxelwatervolume_ = subvoxelporevolume_;
    }
    capillaryvoidvolume_ = capillaryporevolume_ - capillarywatervolume_;
    if (capillaryvoidvolume_ < 0.0)
      capillaryvoidvolume_ = 0.0;

    if (chemSys_->isSaturated()) { // System is saturated
      capillaryvoidvolume_ = 0.0;
      capillarywatervolume_ = watervolume_ - subvoxelporevolume_;
      subvoxelwatervolume_ = subvoxelporevolume_;
    }

    /// From here to the end of this time cycle, vol[ELECTROLYTEID] is the
    /// volume of capillary pore water only, not the total volume of water.
    /// We can always recover the total volume of water by asking ChemicalSystem
    /// for it.

    capillaryporevolume_ = capillaryvoidvolume_ + capillarywatervolume_;
    vol[ELECTROLYTEID] = capillarywatervolume_;
    vol[VOIDID] = capillaryvoidvolume_;

    /// If this has been coded correctly, all we have done is determine how
    /// much of the microstructure water should be assigned to capillary pores
    /// And since only capillary pore water is visible in the microstructure,
    /// we assign the capillary pore water to ELECTROLYTEID and neglect any
    /// subvoxel electrolyte contribution

    if (verbose_) {
      cout << "Lattice::adjustMicrostructureVolumes" << endl;
      cout << "Lattice::adjustMicrostructureVolumesRESULTS:" << endl;
      cout << "Lattice::adjustMicrostructureVolumesAll water volume = "
           << watervolume_ << " m3" << endl;
      cout << "Lattice::adjustMicrostructureVolumesAll void volume = "
           << voidvolume_ << endl;
      cout << "Lattice::adjustMicrostructureVolumesCapillary water volume = "
           << capillarywatervolume_ << " m3" << endl;
      cout << "Lattice::adjustMicrostructureVolumesCapillary void volume = "
           << capillaryvoidvolume_ << " m3" << endl;
      cout << "Lattice::adjustMicrostructureVolumesSubvoxel water volume = "
           << subvoxelwatervolume_ << " m3" << endl;
      cout << "Lattice::adjustMicrostructureVolumesSubvoxel pore volume = "
           << subvoxelporevolume_ << " m3" << endl;
      cout << "Lattice::adjustMicrostructureVolumes" << endl;
      cout.flush();
    }

    //cout << endl
    //     << "   Lattice::adjustMicrostructureVolumes - cyc = " << cyc << endl;
    //cout << "     watervolume_/capPoreVol/capVoidVol/capWaterVol :   "
    //     << watervolume_ << " / " << capillaryporevolume_ << " / "
    //     << capillaryvoidvolume_ << " / " << capillarywatervolume_ << endl;

    ///
    /// End of manual adjustment
    ///

  //}
  //catch (DataException dex) {
  //  throw dex;
  //}
  //catch (out_of_range &oor) {
  //  throw EOBException("Lattice", "adjustMicrostructureVolumes", "volSize",
  //                     volSize, i);
  //}
  //return;
}

void Lattice::adjustMicrostructureVolFracs(vector<string> &names,
                                           const vector<double> vol,
                                           vector<double> &vfrac, int volSize, int cyc) {
  int i = 0;

#ifdef DEBUG
  cout << "Lattice::adjustMicrostructureVolFracs" << endl;
  cout.flush();
#endif

    // Remember there are now two extra slots at the end of vfrac, just
    // like at the end of vol.  These two extra slots hold the capillary
    // pore space and the subvoxel pore space

    // vfrac.clear();
    // vfrac.resize(vol.size(), 0.0);

    // This model tries to reconcile two different volumes:
    //     (1) The total system volume according to GEMS (solid + liquid).
    //         This volume can change during the simulation.  If the
    //         GEMS volume is less than the initial volume of (solid + liquid),
    //         we make up the difference with void space.
    //     (2) The fixed RVE volume which is equal to the number of
    //         voxels multiplied by the volume per voxel
    //
    // We can either let the RVE change volume, which it generally does
    // prior to setting, or we can force it to remain constant, which
    // it approximately does after setting.
    //
    // As of 2023 April 9, we choose to keep the RVE volume fixed, so
    // the only way to reconcile its volume to the GEMS volume is to
    // ensure that the VOLUME FRACTIONS are the same in both.
    //
    // Microstructure volume consists of the sum of (1) the *apparent
    // volumes* of all the solid phases (solid plus subvoxel pores)
    // and (2) the capillary pore volume

    // Find total adjusted microstructure volume, including possibly voids

    // for (i = 0; i < vol.size(); ++i) {
    //     totmicvolume += vol[i];
    //     #ifdef DEBUG
    //         cout << "Lattice::adjustMicrostructureVolFracs Volume("
    //              << names[i] << ") = " << vol[i]
    //              << ", volfrac = " << vfrac[i] << endl;
    //         cout.flush();
    //     #endif
    // }

    /// Above block is to get the current total volume
    /// Replacement below assumes RVE retains constant volume and makes
    /// adjustments based on capillary porosity

    // totmicvolume = chemSys_->getInitMicroVolume();
    // //initialmicrostructurevolume_

#ifdef DEBUvoid
    cout << "Lattice::adjustMicrostructureVolFracsCalculated "
         << "total microstructure volume is " << chemSys_->getInitMicroVolume()
         << endl;
    cout.flush();
#endif

    // Calculate volume fractions based on total GEMS adjusted volume
    for (i = 0; i < volSize; ++i) {
      vfrac[i] = vol[i] / initialmicrostructurevolume_;
    }

    if (verbose_) {
      for (i = 0; i < volSize; ++i) {
        cout << "Lattice::adjustMicrostructureVolFracsVolume "
             << "fraction[" << names[i] << "] should be " << vfrac[i]
             << ", (" << vol[i] << "/"
             << initialmicrostructurevolume_ // totmicvolume
             << ") and volume fraction NOW is "
             << (double)(count_[i]) / (double)(numsites_) << endl;
        cout.flush();
      }
    }

    // Calculate volume fraction of subvoxel porosity and
    // capillary porosity whether saturated or not
    capillaryporevolumefraction_ =
        capillaryporevolume_ / initialmicrostructurevolume_; // totmicvolume;
    subvoxelporevolumefraction_ =
        subvoxelporevolume_ / initialmicrostructurevolume_; // totmicvolume;
}

void Lattice::calcSubvoxelporevolume(vector<double> &vol) {

    // Find the total system volume according to GEMS, in m3
    // units.  The individual microstructure phase volumes
    // have already // been adjusted for subvoxel porosity in the
    // ChemicalSystem::calculateState function

    // The lattice has a fixed volume given by a fixed number
    // of voxels multiplied by the volume per voxel

    // This will hold the subvoxel pore volume (m3)

    subvoxelporevolume_ = 0.0;
    // double phi; // Holds the subvoxel porosity of a microstructurephase
    for (int i = 0; i < vol.size(); ++i) {
      if (i != ELECTROLYTEID && i != VOIDID) {
        subvoxelporevolume_ += (vol[i] * chemSys_->getMicroPhasePorosity(i));
      }
    }

    // The total amount of non-solid space in the microstructure
}

void Lattice::calcSolidvolumewithpores(vector<double> &vol) {

    // Find the total system volume according to GEMS, in m3
    // units.  The individual microstructure phase volumes
    // have already // been adjusted for subvoxel porosity in the
    // ChemicalSystem::calculateState function

    // The lattice has a fixed volume given by a fixed number
    // of voxels multiplied by the volume per voxel

    // This will hold the subvoxel pore volume (m3)

    solidvolumewithpores_ = 0.0;
    for (int i = 0; i < vol.size(); ++i) {
      if (i != ELECTROLYTEID && i != VOIDID) {
        solidvolumewithpores_ += vol[i];
      }
    }
}

void Lattice::calcCapillarywatervolume(vector<double> &vol) {
  calcSubvoxelporevolume(vol);
  capillarywatervolume_ = vol[ELECTROLYTEID] - subvoxelporevolume_;
  if (capillarywatervolume_ < 0.0)
    capillarywatervolume_ = 0.0;
  if (chemSys_->isSaturated()) {
    capillarywatervolume_ = vol[ELECTROLYTEID] - subvoxelporevolume_;
  }
}

void Lattice::calcCapillaryvoidvolume(vector<double> &vol) {
  capillaryporevolume_ = nonsolidvolume_ - subvoxelporevolume_;
  capillaryvoidvolume_ = capillaryporevolume_ - capillarywatervolume_;
}

void Lattice::calculatePoreSizeDistribution(void) {
  // First compose the full pore volume distribution

  vector<double> subpore_volume(volumefraction_.size(), 0.0);

  // Following will hold the subvoxel porosity of a phase
  double phi = 0.0;
  double upper_cutoff_porosity = 0.8;

  // Get the combined pore volume as a function of diameter
  for (int i = 0; i < volumefraction_.size(); ++i) {
    if (i != ELECTROLYTEID && i != VOIDID) {
      phi = chemSys_->getMicroPhasePorosity(i);
      if (phi > upper_cutoff_porosity)
        phi = upper_cutoff_porosity;
      subpore_volume[i] = volumefraction_[i] * phi;
    }
  }

  // At this point subpore_volumes are not normalized
  // but we now know their total volume so we can do that later

  vector<vector<struct PoreSizeVolume>> porevolume;
  porevolume = chemSys_->getPoreSizeDistribution();

  // subpore_volume[i] is the non-normalized pore volume fraction
  // within all of phase i throughout the microstructure.  We
  // now partition that non-normalized volume among the different
  // pore diameters for that phase using the normalized
  // pore volume distribution for that phase
  //
  /// @todo Watch out that pore size distributions are indexed
  /// the same as the THAMES phase id numbers, or else change this
  /// code below.

  for (int i = 0; i < porevolume.size(); ++i) {
    if (i != ELECTROLYTEID && i != VOIDID) {
      for (int j = 0; j < porevolume[i].size(); ++j) {
        porevolume[i][j].volume = porevolume[i][j].volfrac * subpore_volume[i];
      }
    }
  }

  // At this point we have the non-normalized volume of all subvoxel pores
  // as a function of their diameters. Bin them in increments of
  // 1 nanometer

  // Create and initialize the binned volume distribution
  vector<struct PoreSizeVolume> zpsvec;
  vector<vector<struct PoreSizeVolume>> binnedporevolume;
  binnedporevolume.resize(porevolume.size(), zpsvec);

  // Done initializing the binned pore volume distribution

  double maxsize = 1.0;
  double maxmaxsize = 1.0;
  double cumulative_volume = 0.0;
  double cumulative_volfrac = 0.0;
  struct PoreSizeVolume dpsv;

  for (int i = 0; i < porevolume.size(); ++i) {
    if (i != ELECTROLYTEID && i != VOIDID) {
      maxsize = 1.0;
      int j = 0;
      while (j < porevolume[i].size()) {
        while ((j < porevolume[i].size()) && (porevolume[i][j].diam <= maxsize)) {
          cumulative_volume += porevolume[i][j].volume;
          cumulative_volfrac += porevolume[i][j].volfrac;
          j++;
        }
        binnedporevolume[i].push_back(dpsv);
        binnedporevolume[i][binnedporevolume[i].size() - 1].diam = maxsize;
        binnedporevolume[i][binnedporevolume[i].size() - 1].volume =
            cumulative_volume;
        binnedporevolume[i][binnedporevolume[i].size() - 1].volfrac =
            cumulative_volfrac;
        cumulative_volume = cumulative_volfrac = 0.0;
        maxsize += 1.0;
      }
      if (maxsize > maxmaxsize)
        maxmaxsize = maxsize;
    }
  }

  if (verbose_) {
    cout << "Lattice::calculatePoreSizeDistribution  maxmaxsize = "
         << maxmaxsize << endl;
    cout << "Lattice::calculatePoreSizeDistribution  Binned pore distributions"
         << endl;
    for (int i = 0; i < porevolume.size(); ++i) {
      if (i != ELECTROLYTEID && i != VOIDID) {
        cout << "Lattice::calculatePoreSizeDistribution: Distribution " << i
             << " of " << (porevolume.size() - 1) << " now has "
             << porevolume[i].size() << " elements" << endl;
        cout << "Lattice::calculatePoreSizeDistribution  %%%% "
             << chemSys_->getMicroPhaseName(i) << endl;
        for (int j = 0; j < binnedporevolume[i].size(); ++j) {
          cout << "Lattice::calculatePoreSizeDistribution    "
               << "binnedporevolume[" << i << "][" << j
               << "]: diam = " << binnedporevolume[i][j].diam
               << ", volume = " << binnedporevolume[i][j].volume
               << ", volfrac = " << binnedporevolume[i][j].volfrac << endl;
          cout.flush();
          if (binnedporevolume[i][j].volume > 0.0) {
            cout << "Lattice::calculatePoreSizeDistribution: Distribution "
                 << "  %%%%%%% diam = " << binnedporevolume[i][j].diam << " nm,"
                 << " volume = " << binnedporevolume[i][j].volume << ","
                 << " volfrac = " << binnedporevolume[i][j].volfrac << endl;
          }
        }
      }
    }
    cout.flush();
  }

  // Now all subvoxel pores have been binned in 1-nm bins.
  // Combine them into a single distribution for all phases

  double totsubvoxel_volume = 0.0;
  struct PoreSizeVolume dumps;
  dumps.diam = 0.0;
  dumps.volfrac = 0.0;
  dumps.volume = 0.0;
  masterporevolume_.clear();
  masterporevolume_.resize(int(maxmaxsize) - 1, dumps);

  for (int i = 0; i < binnedporevolume.size(); ++i) {
    // cumulative_volume = 0.0;
    // cumulative_volfrac = 0.0;
    for (int j = 0; j < binnedporevolume[i].size(); ++j) {
      masterporevolume_[j].diam = binnedporevolume[i][j].diam;
      masterporevolume_[j].volume += binnedporevolume[i][j].volume;
      totsubvoxel_volume += binnedporevolume[i][j].volume;
      masterporevolume_[j].volfrac = 0.0;
    }
  }

  // Normalize the pore volume distribution relative to the
  // total subvoxel pore volume, stored already in totsubvoxel_volume

  // for (int i = 0; i < masterporevolume_.size(); ++i) {
  //     masterporevolume_[i].volume = masterporevolume_[i].volume
  //                                    / totsubvoxel_volume;
  // }

  if (verbose_) {
    cout << "Lattice::calculatePoreSizeDistribution  Master pore "
         << "size volume fractions" << endl;
    for (int i = 0; i < masterporevolume_.size(); ++i) {
      if (masterporevolume_[i].volume > 0.0) {
        cout << "Lattice::calculatePoreSizeDistribution diam = "
             << masterporevolume_[i].diam << " nm,"
             << " volume = " << masterporevolume_[i].volume << ","
             << " volfrac = " << masterporevolume_[i].volfrac << endl
             << endl;
      }
    }
    cout.flush();
  }

  // At this point we have a complete pore volume distribution
  // for the microstructure.  We next need to determine
  // the saturation state.

  // When we arrive at this function, we know the following
  // information:
  //
  // Volume fraction of saturated capillary pores on a total
  // microstructure volume basis = volumefraction_[ELECTROLYTEID]
  //
  // Volume fraction of unsaturated capillary pores
  // on a total microstructure volume basis
  // = volumefraction_[VOIDID]
  //
  // Volume fraction of capillary pores on a total
  // microstructure volume basis = capillaryporevolumefraction_
  //
  // Volume fraction of subvoxel pores on a total
  // microstructure volume basis = subvoxelporevolumefraction_
  //
  // We do NOT yet know the fraction of subvoxel pores
  // that are saturated

  // water_volume is the absolute volume (m3) of aqueous
  // solution according to GEMS, whether in capillary
  // pores, subvoxel pores, or squeezed outside the
  // system due to lack of porosity to contain it.

  double water_volume = chemSys_->getMicroPhaseVolume(ELECTROLYTEID);

  // microvol is the absolute volume (m3) of all the defined
  // microstructure phases, including subvoxel pores in
  // solid phases

  // double microvol = chemSys_->getMicroVolume();

  // microvol is the absolute initial volume (m3) of all the defined
  // microstructure phases, including subvoxel pores in
  // solid phases

  // double initmicrovol = chemSys_->getInitMicroVolume();

  // This is the volume fraction of liquid water whether
  // in capillary or subvoxel porosity, on a total
  // microstructure volume basis

  // We may have had to push some of the system's
  // capillary water outside the boundary of the
  // microstructure due to lack of space, so maybe
  // water_volume is larger than the actual volume
  // fraction of water available within the microstructure
  // We can check for that now, though.

  double excesswater =
      water_volume - capillaryporevolume_ - subvoxelporevolume_;

  // If excess water > 0 then we have a completely
  // saturated system and it is easy to write out
  // the water partitioning among pore sizes

  // bool isfullysaturated = false;
  // if (excesswater > 0.0)
  //   isfullysaturated = true;

  //double water_volfrac = water_volume / microstructurevolume_;
  double water_volfrac = water_volume / initialmicrostructurevolume_;

  // This is the total porosity including capillary
  // pore volume fraction and subvoxel pore volume
  // fraction, all on a total microstructure volume
  // basis

  // double pore_volfrac =
  //     capillaryporevolumefraction_ + subvoxelporevolumefraction_;

  //  cout << "Lattice::calculatePoreSizeDistribution:" << endl;
  //  cout << "Lattice::calculatePoreSizeDistribution  "
  //       << "==== water_volume = " << water_volume << endl;
  //  cout << "Lattice::calculatePoreSizeDistribution  "
  //       << "======== microstructurevolume = " << microstructurevolume_ <<
  //       endl;
  //  cout << "Lattice::calculatePoreSizeDistribution  "
  //       << "==== initmicrostructurevolume = "
  //        << initialmicrostructurevolume_ << endl;
  //  cout << "Lattice::calculatePoreSizeDistribution  === "
  //       << "water_volfrac = " << water_volfrac << endl;
  //  cout << "Lattice::calculatePoreSizeDistribution  ==== "
  //       << "pore_volfrac = " << pore_volfrac << endl;
  //  cout << "Lattice::calculatePoreSizeDistribution  "
  //       << "====== (cap = " << capillaryporevolumefraction_
  //       << ", subvox = " << subvoxelporevolumefraction_
  //       << ")" << endl;
  //  cout << "Lattice::calculatePoreSizeDistribution  ====== "
  //       << "void fraction = " << volumefraction_[VOIDID] << endl;
  //  cout << "Lattice::calculatePoreSizeDistribution  ====== "
  //       << "is fully saturated? ";
  //  if (isfullysaturated) {
  //      cout << " YES" << endl;
  //  } else {
  //      cout << " NO" << endl;
  //  }
  //  cout.flush();

  // Only worry about unsaturated subvoxel pores if the capillary
  // pores are completely empty of water

  if (volumefraction_[ELECTROLYTEID] >= 0.0) {
    // The next variable represents the volume fraction of subvoxel
    // pores of a given size, on a total microstructure volume basis
    double volfrac_avail = 0.0;
    double volfrac_filled = 0.0;
    // masterporevolume is already a kind of volume fraction
    // because it has been normalized to the total subvoxel pore volume
    //    cout << "Lattice::calculatePoreSizeDistribution  "
    //          << "Master pore size volume filling" << endl;
    //    cout.flush();
    for (int i = 0; i < masterporevolume_.size(); ++i) {
      volfrac_avail = masterporevolume_[i].volume * subvoxelporevolumefraction_;
      volfrac_filled = water_volfrac / volfrac_avail;
      if (volfrac_filled > 1.0)
        volfrac_filled = 1.0;
      masterporevolume_[i].volfrac = volfrac_filled;
      if (!(chemSys_->isSaturated())) { // System is sealed
        water_volfrac -= volfrac_avail;
        if (water_volfrac < 0.0)
          water_volfrac = 0.0;
      }
      //    if (masterporevolume_[i].volume > 0.0) {
      //      cout << "Lattice::calculatePoreSizeDistribution diam = "
      //       << masterporevolume_[i].diam << " nm,"
      //       << " vfrac avail = " << volfrac_avail << ","
      //       << " volume = " << masterporevolume_[i].volume << ","
      //       << " volfilled = " << masterporevolume_[i].volfrac << ","
      //       << " waterfrac left = " << water_volfrac << endl;
      //    cout.flush();
      //    }
    }
    // cout << endl;
    cout.flush();
  }

  return;
}

void Lattice::writePoreSizeDistribution(double curtime, const int simtype,
                                        const string &root) {

  // At this point we have a complete pore volume distribution
  // for the microstructure.  We next need to determine
  // the saturation state.

  // When we arrive at this function, we know the following
  // information:
  //
  // Volume fraction of saturated capillary pores on a total
  // microstructure volume basis = volumefraction_[ELECTROLYTEID]
  //
  // Volume fraction of unsaturated capillary pores
  // on a total microstructure volume basis
  // = volumefraction_[VOIDID]
  //
  // Volume fraction of capillary pores on a total
  // microstructure volume basis = capillaryporevolumefraction_
  //
  // Volume fraction of subvoxel pores on a total
  // microstructure volume basis = subvoxelporevolumefraction_
  //
  // We do NOT yet know the fraction of subvoxel pores
  // that are saturated

  // water_volume is the absolute volume (m3) of aqueous
  // solution according to GEMS, whether in capillary
  // pores, subvoxel pores, or squeezed outside the
  // system due to lack of porosity to contain it.

  if (verbose_) {
    cout << "In Lattice::writePoreSizeDistribution" << endl;
    cout.flush();
  }

  // double water_volume = chemSys_->getMicroPhaseVolume(ELECTROLYTEID);

  // This is the volume fraction of liquid water whether
  // in capillary or subvoxel porosity, on a total
  // microstructure volume basis

  // We may have had to push some of the system's
  // capillary water outside the boundary of the
  // microstructure due to lack of space, so maybe
  // watervolume_ is larger than the actual volume
  // fraction of water available within the microstructure
  // We can check for that now, though.

  double excesswater =
      watervolume_ - capillaryporevolume_ - subvoxelporevolume_;

  // If excess water > 0 then we have a completely
  // saturated system and it is easy to write out
  // the water partitioning among pore sizes

  // bool isfullysaturated = false;
  // if (excesswater > 0.0)
  //   isfullysaturated = true;
  
  //microstructurevolume_ = chemSys_->getMicroVolume();
  //double water_volfrac = watervolume_ / microstructurevolume_;
  double water_volfrac = watervolume_ / initialmicrostructurevolume_;

  // cout << endl << "--> watervolume_         : " << watervolume_ << endl;
  // cout << "--> microstructurevolume_: " << microstructurevolume_ << endl;
  // cout << "--> capillaryporevolume_ : " << capillaryporevolume_ << endl;
  // cout << "--> subvoxelporevolume_  : " << subvoxelporevolume_ << endl;
  // cout << "--> excesswater          : " << excesswater << endl;

  // This is the total porosity including capillary
  // pore volume fraction and subvoxel pore volume
  // fraction, all on a total microstructure volume
  // basis

  double pore_volfrac =
      capillaryporevolumefraction_ + subvoxelporevolumefraction_;

  string ofileName(root);
  ostringstream ostr1, ostr2;
  // Add the time in minutes
  ostr1 << setfill('0') << setw(6) << (int)((curtime * 24.0 * 60.0) + 0.5);
  ostr2 << setprecision(3) << temperature_;
  string timestr(ostr1.str());
  string tempstr(ostr2.str());
  ofileName =
      ofileName + "_PoreSizeDistribution." + timestr + "." + tempstr + ".csv";

  ofstream out(ofileName.c_str());
  try {
    if (!out.is_open()) {
      throw FileException("Lattice", "writePoreSizeDistribution", ofileName,
                          "Could not open");
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  // Write the header

  if (verbose_) {
    cout << "Time = " << (curtime * 24.0) << " h" << endl;
    cout << "Capillary pore volume fraction (> 100 nm) = "
         << capillaryporevolumefraction_ << endl;
    cout << "Capillary void volume fraction = " << volumefraction_[VOIDID]
         << endl;
    cout << "Saturated capillary pore volume fraction = "
         << capillaryporevolumefraction_ - volumefraction_[VOIDID] << endl;
    cout << "Nanopore volume fraction (<= 100 nm) = "
         << subvoxelporevolumefraction_ << endl;
    cout << "Total pore volume fraction = " << pore_volfrac << endl << endl;
    cout << "Total void volume fraction = " << volumefraction_[VOIDID]
         << endl;
    cout << "Pore size saturation data:" << endl;
    cout << "Diameter (nm),Volume Fraction,Fraction Saturated" << endl;
    cout << "Masterporevolume size = " << masterporevolume_.size() << endl;
    cout.flush();
  }

  out << "Time = " << (curtime * 24.0) << " h" << endl;
  out << "Capillary pore volume fraction (> 100 nm) = "
      << capillaryporevolumefraction_ << endl;
  out << "Capillary void volume fraction = " << volumefraction_[VOIDID]
      << endl;
  out << "Saturated capillary pore volume fraction = "
      << capillaryporevolumefraction_ - volumefraction_[VOIDID] << endl;
  out << "Nanopore volume fraction (<= 100 nm) = "
      << subvoxelporevolumefraction_ << endl;
  out << "Total pore volume fraction = " << pore_volfrac << endl;
  out << "Total void volume fraction = " << volumefraction_[VOIDID] << endl;
  out << "Pore size saturation data:" << endl;
  out << "Diameter (nm),Volume Fraction,Fraction Saturated" << endl;
  for (int i = 0; i < masterporevolume_.size(); i++) {
    if (masterporevolume_[i].volume > 0.0) {
      out << masterporevolume_[i].diam << "," << masterporevolume_[i].volume
          << "," << masterporevolume_[i].volfrac << endl;
      out.flush();
    }
  }

  // This is the volume fraction of capillary pore water,
  // on a total microstructure volume basis, already
  // calculated and stored when altering the microstructure

  double capwater_volfrac = water_volfrac;
  double capvoid_volfrac = volumefraction_[VOIDID] +
                           (volumefraction_[ELECTROLYTEID] - water_volfrac);
  double capspace_volfrac = capvoid_volfrac + capwater_volfrac;

  // cout << endl << "--> water_volfrac           : " << water_volfrac << endl;
  // cout << "--> volumefraction_[VOIDID]         : " << volumefraction_[VOIDID] << endl;
  // cout << "--> volumefraction_[ELECTROLYTEID]  : " << volumefraction_[ELECTROLYTEID] << endl;
  // cout << "--> capvoid_volfrac/capspace_volfrac: " << capvoid_volfrac << " / " << capspace_volfrac << endl;

  out << ">" << masterporevolume_[masterporevolume_.size() - 1].diam << ","
      << capspace_volfrac << "," << (1.0 - volumefraction_[VOIDID]) << endl;
  out.flush();

  out.close();

  return;
}

void Lattice::writeMicroColors(const string &root) {
  unsigned int i, j, k;
  string ofileName(root);
  ofileName.append("_Colors.csv");

  ofstream out(ofileName.c_str());
  try {
    if (!out.is_open()) {
      throw FileException("Lattice", "writeLattice", ofileName,
                          "Could not open");
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  int microPhaseId;
  string microPhaseName;
  //int red, green, blue;
  vector<int> colors;
  out << numMicroPhases_ << endl;
  for (int i = 0; i < numMicroPhases_; i++) {
    microPhaseId = chemSys_->getMicroPhaseId(i);
    microPhaseName = chemSys_->getMicroPhaseName(microPhaseId);
    //colors = chemSys_->getColor(microPhaseId);
    colors = chemSys_->getRGB(microPhaseId);
    out << microPhaseId << " " << microPhaseName << " "
        << colors[0] << " " << colors[1] << " " << colors[2] << endl;
  }

  out.flush();
  out.close();
  return;
}

void Lattice::writeLattice(double curtime, const int simtype,
                           const string &root) {
  unsigned int i, j, k;
  string ofileName(root);
  ostringstream ostr1, ostr2;
  ostr1 << setfill('0') << setw(6)
        << (int)((curtime * 24.0 * 60.0) + 0.5); // minutes
  ostr2 << setprecision(3) << temperature_;
  string timestr(ostr1.str());
  string tempstr(ostr2.str());
  ofileName = ofileName + "." + timestr + "m." + tempstr + ".img";
  if (verbose_) {
    cout << "    In Lattice::writeLattice, curtime = " << curtime
         << ", timestr = " << timestr << endl;
    cout.flush();
  }

  ofstream out(ofileName.c_str());
  try {
    if (!out.is_open()) {
      throw FileException("Lattice", "writeLattice", ofileName,
                          "Could not open");
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  // Write image header information first

  out << VERSIONSTRING << " " << thamesVersion_ << endl;
  out << XSIZESTRING << " " << xdim_ << endl;
  out << YSIZESTRING << " " << ydim_ << endl;
  out << ZSIZESTRING << " " << zdim_ << endl;
  out << IMGRESSTRING << " " << resolution_ << endl;

  for (k = 0; k < zdim_; k++) {
    for (j = 0; j < ydim_; j++) {
      for (i = 0; i < xdim_; i++) {
        int index = getIndex(i, j, k);
        out << site_[index].getMicroPhaseId() << endl;
      }
    }
  }

  out.close();

  // The next block is implemented only if we are dealing with sulfate attack
  if (simtype == SULFATE_ATTACK) {

    ofileName = root;
    ofileName = ofileName + "." + timestr + "." + tempstr + ".img.damage";

    ofstream out1(ofileName.c_str());
    try {
      if (!out1.is_open()) {
        throw FileException("Lattice", "writeLattice", ofileName,
                            "Could not open");
      }
    } catch (FileException fex) {
      fex.printException();
      exit(1);
    }

    // Write image header information first

    out1 << VERSIONSTRING << " " << thamesVersion_ << endl;
    out1 << XSIZESTRING << " " << xdim_ << endl;
    out1 << YSIZESTRING << " " << ydim_ << endl;
    out1 << ZSIZESTRING << " " << zdim_ << endl;
    out1 << IMGRESSTRING << " " << resolution_ << endl;

    int DAMAGEID =
        100; // Some large number that cannot represent any other phase
    for (k = 0; k < zdim_; k++) {
      for (j = 0; j < ydim_; j++) {
        for (i = 0; i < xdim_; i++) {
          int index = getIndex(i, j, k);
          if (site_[index].IsDamage()) {
            out1 << DAMAGEID << endl;
          } else {
            out1 << site_[index].getMicroPhaseId() << endl;
          }
        }
      }
    }

    out1.close();
  } // The above block is implemented only if we are dealing with sulfate attack
}

void Lattice::writeLatticeIni(double curtime) {
  unsigned int i, j, k;
  ostringstream ostr1, ostr2;
  ostr1 << setfill('0') << setw(6)
        << (int)((curtime * 24.0 * 60.0) + 0.5); // minutes
  ostr2 << setprecision(3) << temperature_;
  string timestr(ostr1.str());
  string tempstr(ostr2.str());
  if (verbose_) {
    cout << "    In Lattice::writeLatticeIni, curtime = " << curtime
         << ", timestr = " << timestr << endl;
    cout.flush();
  }

  // Write image header information first

  string ofileName1;
  ofileName1 = "imgIni." + timestr + "m." + tempstr + ".img";
  ofstream out1(ofileName1.c_str());
  out1 << VERSIONSTRING << " " << thamesVersion_ << endl;
  out1 << XSIZESTRING << " " << xdim_ << endl;
  out1 << YSIZESTRING << " " << ydim_ << endl;
  out1 << ZSIZESTRING << " " << zdim_ << endl;
  out1 << IMGRESSTRING << " " << resolution_ << endl;
  for (k = 0; k < zdim_; k++) {
    for (j = 0; j < ydim_; j++) {
      for (i = 0; i < xdim_; i++) {
        int index = getIndex(i, j, k);
        out1 << site_[index].getMicroPhaseId() << endl;
      }
    }
  }
  out1.close();
}

void Lattice::writeLatticeXYZ(double curtime, const int simtype,
                              const string &root) {
  unsigned int i, j, k;
  string ofileName(root);
  ostringstream ostr1, ostr2;
  ostr1 << setfill('0') << setw(6)
        << (int)((curtime * 24.0 * 60.0) + 0.5); // minutes
  ostr2 << setprecision(3) << temperature_;
  string timestr(ostr1.str());
  string tempstr(ostr2.str());
  ofileName = ofileName + "allSites." + timestr + "m." + tempstr + ".xyz";
  if (verbose_) {
    cout << "    In Lattice::writeLatticeXYZ, curtime = " << curtime
         << ", timestr = " << timestr << endl;
    cout.flush();
  }

  ofstream out(ofileName.c_str());

  int numvox = xdim_ * ydim_ * zdim_;

  // Write the file headers
  out << numvox << endl; // Number of voxels to visualize
  out << "Lattice=\"" << (float)xdim_ << " 0.0 0.0 0.0 " << (float)ydim_
      << " 0.0 0.0 0.0 " << (float)zdim_ << "\" ";
      out << "Time=" << timestr << " ";
  //out << "Properties=pos:R:3:color:R:3:transparency:R:1 " << endl;
  out << "Properties=phaseID:I:1:element:S:1:pos:R:3:vector_color:R:3:"
         "radius:R:1:transparency:R:1" << endl;

  // Loop over all voxels and write out the solid ones

  float x, y, z;
  int mPhId;
  vector<int> colors;
  string symb;

  for (int i = 0; i < numsites_; i++) {
    x = site_[i].getX();
    y = site_[i].getY();
    z = site_[i].getZ();
    mPhId = site_[i].getMicroPhaseId();
    symb = cfgElem_[mPhId].symb; //mPhIdgetElemSymb(mPhId);
    //mPhName = chemSys_->getMicroPhaseName(mPhId);
    //colors = chemSys_->getColor(mPhId);
    colors = chemSys_->getRGB(mPhId);
    out << mPhId << "\t" << symb << "\t" << x << "\t" << y << "\t" << z
        << "\t" << colors[0] << "\t" << colors[1] << "\t" << colors[2]
        << "\t" << particRadius_ << "\t" << "0.0" << endl;
    //out << x << "\t" << y << "\t" << z
    //    << "\t" << colors[0] << "\t" << colors[1] << "\t" << colors[2]
    //    << "\t" << "0.0" << endl;

  }
  out.close();
}

void Lattice::writeLatticeCFG(double curtime, const int simtype,
                              const string &root) {

  string ofileNameCFG(root);
  string ofileNameUSR(root);
  ostringstream ostr1, ostr2;
  ostr1 << setfill('0') << setw(6)
        << (int)((curtime * 24.0 * 60.0) + 0.5); // minutes
  ostr2 << setprecision(3) << temperature_;
  string timestr(ostr1.str());
  string tempstr(ostr2.str());
  ofileNameCFG = ofileNameCFG + "allSites." + timestr + "m." + tempstr + ".cfg";
  ofileNameUSR = ofileNameUSR + "allSites." + timestr + "m." + tempstr + ".usr";
  if (verbose_) {
    cout << "    In Lattice::writeLatticeCFG, curtime = " << curtime
         << ", timestr = " << timestr << endl;
    cout.flush();
  }

  ofstream outCFG(ofileNameCFG.c_str());
  ofstream outUSR(ofileNameUSR.c_str());

  int ord;
  double x, y, z;
  vector<int> colors;
  int i, mPhId;
  bool mPhNotExist;

  int numPart;
  // numPart = numsites_ - count_[0] - count_[1];
  numPart = numsites_;

  // write cfg header
  outCFG << "Number of particles = " << numPart << endl;

  outCFG << "A = 3.0 Angstrom (basic length-scale)" << endl;

  outCFG << "H0(1,1) = " << xdim_-1 << " A" << endl;
  outCFG << "H0(1,2) = 0 A" << endl;
  outCFG << "H0(1,3) = 0 A" << endl;

  outCFG << "H0(2,1) = 0 A" << endl;
  outCFG << "H0(2,2) = " << ydim_-1 << " A" << endl;
  outCFG << "H0(2,3) = 0 A" << endl;

  outCFG << "H0(3,1) = 0 A" << endl;
  outCFG << "H0(3,2) = 0 A" << endl;
  outCFG << "H0(3,3) = " << zdim_-1 << " A" << endl;
  outCFG << ".NO_VELOCITY." << endl;
  outCFG << "entry_count = 3" << endl;

  // write cfg coordinates & usr file (colors)
  ord = 0;
  for (mPhId = 0; mPhId < numMicroPhases_; mPhId++) {
    // for (mPhId = FIRST_SOLID; mPhId < numMicroPhases_; mPhId++) {
    mPhNotExist = true;
    //colors = chemSys_->getColor(mPhId);
    colors = chemSys_->getRGB(mPhId);
    for (i = 0; i < numsites_; i++) {
      if (site_[i].getMicroPhaseId() == mPhId) {
        if (mPhNotExist) {
          outCFG << "  " << cfgElem_[mPhId].mass << endl;
          outCFG << cfgElem_[mPhId].symb << endl;
          mPhNotExist = false;
        }
        x = (site_[i].getX()) / (double)(xdim_-1);
        y = (site_[i].getY()) / (double)(ydim_-1);
        z = (site_[i].getZ()) / (double)(zdim_-1);
        outCFG << x << "\t" << y << "\t" << z << endl;
        outUSR << ord << "\t" << colors[0] << "\t" << colors[1] << "\t"
               << colors[2] << "\t0.7051" << endl;
        ord++;
      }
    }
  }
  outCFG.close();
  outUSR.close();

  if (ord != numPart) {
    cout << endl
         << "*************** error in writeLatticeCFG! ******************"
         << endl;
    cout << endl << "                    ord != numsites_" << endl;
    cout << "                                   STOP" << endl;
    exit(1);
  }
}

void Lattice::writeDamageLattice(double curtime, const string &root) {
  unsigned int i, j, k;
  string ofileName(root);
  ostringstream ostr1, ostr2;
  ostr1 << setfill('0') << setw(6)
        << (int)((curtime * 24.0 * 60.0) + 0.5); // minutes
  ostr2 << setprecision(3) << temperature_;
  string timestr(ostr1.str());
  string tempstr(ostr2.str());
  ofileName = ofileName + "." + timestr + "m." + tempstr + ".img";

  ofstream out(ofileName.c_str());
  try {
    if (!out.is_open()) {
      throw FileException("Lattice", "writeLattice", ofileName,
                          "Could not open");
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  // Write image header information first

  out << VERSIONSTRING << " " << thamesVersion_ << endl;
  out << XSIZESTRING << " " << xdim_ << endl;
  out << YSIZESTRING << " " << ydim_ << endl;
  out << ZSIZESTRING << " " << zdim_ << endl;
  out << IMGRESSTRING << " " << resolution_ << endl;

  for (k = 0; k < zdim_; k++) {
    for (j = 0; j < ydim_; j++) {
      for (i = 0; i < xdim_; i++) {
        int index = getIndex(i, j, k);
        if (site_[index].IsDamage()) {
          out << "1" << endl;
        } else {
          out << "0" << endl;
        }
      }
    }
  }

  out.close();
}

void Lattice::writeLatticePNG(double curtime, const int simtype,
                              const string &root) {
  unsigned int i, j, k;
  string ofileName(root);
  string ofbasename(root);
  string ofpngname(root);
  string ofpngbasename(root);

  vector<double> dumvec;
  vector<unsigned int> idumvec;
  vector<vector<unsigned int>> image;
  vector<vector<double>> dshade;
  dumvec.resize(ydim_, 0.0);
  idumvec.resize(ydim_, 0);
  dshade.resize(xdim_, dumvec);
  image.resize(xdim_, idumvec);
  bool done;

  ///
  /// Construct the name of the output file
  ///

  ostringstream ostr1, ostr2;
  ostr1 << setfill('0') << setw(6)
        << (int)((curtime * 24.0 * 60.0) + 0.5); // minutes
  ostr2 << setprecision(3) << temperature_;
  string timestr(ostr1.str());
  string tempstr(ostr2.str());
  string buff;
  ofileName = ofileName + "." + timestr + "m." + tempstr + ".ppm";
  ofpngname = ofpngname + "." + timestr + "m." + tempstr + ".png";

  ///
  /// Open the output file
  ///

  ofstream out(ofileName.c_str());
  try {
    if (!out.is_open()) {
      throw FileException("Lattice", "writeLatticePNG", ofileName,
                          "Could not open");
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  ///
  /// Write PPM header for full color image
  ///

  out << "P3" << endl;
  out << xdim_ << " " << ydim_ << endl;
  out << COLORSATVAL << endl;

  unsigned int slice = zdim_ / 2;
  unsigned int nd, izz, valout;
  unsigned int sitenum;
  for (j = 0; j < ydim_; j++) {
    for (i = 0; i < xdim_; i++) {
      if (deptheffect_) {
        done = false;
        nd = 0;
        // sitenum = getIndex(i, j, slice);
        izz = slice;
        do {
          sitenum = getIndex(i, j, izz);
          if (nd == 10 || site_[sitenum].getMicroPhaseId() > 1) {
            done = true;
          } else {
            nd++;
            izz++;
            if (izz >= zdim_)
              izz -= zdim_;
          }
        } while (!done);
        sitenum = getIndex(i, j, izz);
        image[i][j] = site_[sitenum].getMicroPhaseId();
        dshade[i][j] = 0.1 * (10.0 - ((double)nd));
      } else {
        sitenum = getIndex(i, j, slice);
        image[i][j] = site_[sitenum].getMicroPhaseId();
        dshade[i][j] = 1.0;
      }
    }
  }

  int red, green, blue;
  vector<int> colors;
  for (j = 0; j < ydim_; j++) {
    for (i = 0; i < xdim_; i++) {
      //colors = chemSys_->getColor(image[i][j]);
      colors = chemSys_->getRGB(image[i][j]);
      red = dshade[i][j] * colors[0] + 0.5;
      green = dshade[i][j] * colors[1] + 0.5;
      blue = dshade[i][j] * colors[2] + 0.5;
      out << red << " " << green << " " << blue << endl;
    }
  }

  out.close();

  ///
  /// Execute system call to convert PPM to PNG.
  ///
  /// @warning This relies on installation of ImageMagick
  ///

  buff = "convert " + ofileName + " " + ofpngname;
  system(buff.c_str());
  return;
}

void Lattice::writeDamageLatticePNG(double curtime, const string &root) {
  unsigned int i, j, k;
  string ofileName(root);
  string ofbasename(root);
  string ofpngname(root);
  string ofpngbasename(root);

  vector<double> dumvec;
  vector<unsigned int> idumvec;
  vector<vector<unsigned int>> image;
  vector<vector<double>> dshade;
  dumvec.resize(ydim_, 0.0);
  idumvec.resize(ydim_, 0);
  dshade.resize(xdim_, dumvec);
  image.resize(xdim_, idumvec);
  bool done;

  ///
  /// Construct the name of the output file
  ///

  ostringstream ostr1, ostr2;
  ostr1 << (int)((curtime * 100.0) + 0.5); // hundredths of an hour
  ostr2 << setprecision(3) << temperature_;
  string timestr(ostr1.str());
  string tempstr(ostr2.str());
  string buff;
  ofileName = ofileName + "." + timestr + "." + tempstr + ".ppm";
  ofpngname = ofpngname + "." + timestr + "." + tempstr + ".png";

  ///
  /// Open the output file
  ///

  ofstream out(ofileName.c_str());
  try {
    if (!out.is_open()) {
      throw FileException("Lattice", "writeLatticePNG", ofileName,
                          "Could not open");
    }
  } catch (FileException fex) {
    fex.printException();
    exit(1);
  }

  ///
  /// Write PPM header for full color image
  ///

  out << "P3" << endl;
  out << ydim_ << " " << zdim_ << endl;
  out << COLORSATVAL << endl;

  unsigned int slice = zdim_ / 2;
  unsigned int nd, izz, valout;
  unsigned int sitenum;
  for (j = 0; j < ydim_; j++) {
    for (i = 0; i < xdim_; i++) {
      done = false;
      nd = 0;
      // sitenum = getIndex(i, j, slice);
      izz = slice;
      do {
        sitenum = getIndex(i, j, izz);
        if (nd == 10 || site_[sitenum].getMicroPhaseId() > 1) {
          done = true;
        } else {
          nd++;
          izz++;
          if (izz >= zdim_)
            izz -= zdim_;
        }
      } while (!done);
      sitenum = getIndex(i, j, izz);
      if (site_[sitenum].IsDamage()) {
        image[i][j] = 1;
      } else {
        image[i][j] = 5;
      }
      dshade[i][j] = 1.0;
      /*
      dshade[j][k] = 0.1 * (10.0 - ((double)nd));
      */
    }
  }

  int red, green, blue;
  vector<int> colors;
  for (j = 0; j < ydim_; j++) {
    for (i = 0; i < xdim_; i++) {
      //colors = chemSys_->getColor(image[i][j]);
      colors = chemSys_->getRGB(image[i][j]);
      red = dshade[i][j] * colors[0] + 0.5;
      green = dshade[i][j] * colors[1] + 0.5;
      blue = dshade[i][j] * colors[2] + 0.5;
      out << red << " " << green << " " << blue << endl;
    }
  }

  out.close();

  ///
  /// Execute system call to convert PPM to PNG.
  ///
  /// @warning This relies on installation of ImageMagick
  ///

  buff = "convert " + ofileName + " " + ofpngname;
  system(buff.c_str());
  return;
}

void Lattice::makeMovie(const string &root) {
  unsigned int i, j, k;
  string ofileName(root);
  string ofbasename(root);
  string ofgifileName(root);
  string ofgifbasename(root);

  vector<double> dumvec;
  vector<unsigned int> idumvec;
  vector<vector<unsigned int>> image;
  vector<vector<double>> dshade;
  dumvec.resize(ydim_, 0.0);
  idumvec.resize(ydim_, 0);
  dshade.resize(zdim_, dumvec);
  image.resize(zdim_, idumvec);
  bool done;

  ///
  /// Construct the name of the output file
  ///

  string buff;
  ostringstream ostr1, ostr2, ostr3;
  ostr1 << (int)(time_ * 100.0); // hundredths of an hour
  ostr2 << setprecision(3) << temperature_;
  string timestr(ostr1.str());
  string tempstr(ostr2.str());

  ///
  /// Loop over number of slices in the x direction, making one image per slice
  /// and appending it to the end of the master file.
  ///

  for (k = 10; k < zdim_; k++) {

    ///
    /// Open the output file.
    ///

    ostr3.clear();
    ostr3 << (int)(k); // x slice number
    string istr(ostr3.str());
    ofileName =
        ofbasename + "." + timestr + "." + tempstr + "." + istr + ".ppm";
    ofgifileName =
        ofgifbasename + "." + timestr + "." + tempstr + "." + istr + ".gif";

    ofstream out(ofileName.c_str());
    if (!out.is_open()) {
      throw FileException("Lattice", "makeMovie", ofileName, "Could not open");
    }

    ///
    /// Write PPM header for full color image.
    ///

    out << "P3" << endl;
    out << xdim_ << " " << ydim_ << endl;
    out << COLORSATVAL << endl;

    unsigned int slice = k;
    unsigned int nd, izz, valout;
    unsigned int sitenum;
    for (j = 0; j < ydim_; j++) {
      for (i = 0; i < xdim_; i++) {
        done = false;
        nd = 0;
        // sitenum = getIndex(i, j, slice);
        izz = slice;
        do {
          sitenum = getIndex(i, j, izz);
          if (nd == 10 || site_[sitenum].getMicroPhaseId() > 1) {
            done = true;
          } else {
            nd++;
            izz++;
            if (izz >= zdim_)
              izz -= zdim_;
          }
        } while (!done);
        sitenum = getIndex(i, j, izz);
        image[i][j] = site_[sitenum].getMicroPhaseId();
        dshade[i][j] = 0.1 * (10.0 - ((double)nd));
      }
    }

    int red, green, blue;
    vector<int> colors;
    for (j = 0; j < ydim_; j++) {
      for (i = 0; i < xdim_; i++) {
        //colors = chemSys_->getColor(image[i][j]);
        colors = chemSys_->getRGB(image[i][j]);
        red = dshade[i][j] * colors[0] + 0.5;
        green = dshade[i][j] * colors[1] + 0.5;
        blue = dshade[i][j] * colors[2] + 0.5;
        out << red << " " << green << " " << blue << endl;
      }
    }
    out.close();

    ///
    /// Execute system call to convert PPM to GIF.
    ///
    /// @warning This relies on installation of ImageMagick
    ///

    buff = "convert " + ofileName + " " + ofgifileName;
    system(buff.c_str());
  }

  ///
  /// Execute system call to convert GIF frames to animated GIF.
  ///
  /// @warning This relies on installation of gifsicle
  ///

  buff = "gifsicle --delay=10 " + ofgifbasename + "*.gif > " + ofgifbasename +
         ".movie.gif";
  system(buff.c_str());
}

vector<unsigned int> Lattice::writeSubVolume(string fileName, Site *centerste,
                                             int size) {
  ofstream out(fileName.c_str());

  out << "Version: 7.0" << endl;
  out << "X_Size: 3" << endl;
  out << "Y_Size: 3" << endl;
  out << "Z_Size: 3" << endl;
  out << "Image_Resolution: 1" << endl;

  vector<unsigned int> alnb = getNeighborhood(centerste->getId(), size);

  for (int j = 0; j < alnb.size(); j++) {
    int phaseid = site_[alnb[j]].getMicroPhaseId();
    out << phaseid << endl;
  }
  out.close();

  return alnb;
}

void Lattice::applyExp(vector<unsigned int> alnb, double exp) {
  Site *ste;
  if (exp > 0.0) {
    for (int i = 0; i < alnb.size(); i++) {
      ste = &site_[alnb[i]];
      ste->setExpansionStrain(exp);

      map<int, vector<double>>::iterator p = expansion_.find(ste->getId());
      if (p != expansion_.end()) {
        for (int j = 0; j < (p->second).size(); j++) {
          (p->second)[j] = ste->getExpansionStrain();
        }
      } else {
        vector<double> expval;
        expval.clear();
        expval.resize(3, ste->getExpansionStrain());
        expansion_.insert(make_pair(ste->getId(), expval));
      }

      map<int, vector<int>>::iterator pp =
          expansion_coordin_.find(ste->getId());
      if (pp == expansion_coordin_.end()) {
        vector<int> coordin;
        coordin.clear();
        coordin.resize(3, 0);
        coordin[0] = ste->getX();
        coordin[1] = ste->getY();
        coordin[2] = ste->getZ();
        expansion_coordin_.insert(make_pair(ste->getId(), coordin));
      }
    }
  }

  return;
}

double Lattice::getSurfaceArea(int phaseid) {
  double surface1 = 0.0, surface2 = 0.0;
  Site *ste, *stenb;
  vector<Isite> isite = interface_[phaseid].getDissolutionSites();

  ///
  /// Method 1: Surface area is equal to the wmc (prop to volume of surrounding
  /// pores)
  ///
  int isiteSize = isite.size();
  for (int i = 0; i < isiteSize; i++) {
    ste = &site_[isite[i].getId()];
    surface1 += ste->getWmc();
  }

  ///
  /// Method 2: Surface area is related to interior porosity of neighbor sites
  ///

  for (int i = 0; i < site_.size(); i++) {
    ste = &site_[i];
    if (ste->getMicroPhaseId() == phaseid) {
      for (int j = 0; j < ste->nbSize(1); j++) {
        stenb = ste->nb(j);
        surface2 += chemSys_->getMicroPhasePorosity(stenb->getMicroPhaseId());
      }
    }
  }

  if (verbose_) {
    cout << "Lattice::getSurfaceArea Surface area of phase " << phaseid
         << " calculated by method 1 is: " << surface1 << endl;
    cout << "Lattice::getSurfaceArea Surface area of phase " << phaseid
         << " calculated by method 2 is: " << surface2 << endl;
  }

  ///
  /// Use Method 2
  ///

  surfacearea_ = surface2;

  return surfacearea_;
}

int Lattice::findDomainSize(int siteid, int maxsize) {

  int boxhalf = maxsize / 2;
  int nfound = 0;

  int phaseid = (site_[siteid]).getMicroPhaseId();
  int qx = (site_[siteid]).getX();
  int qy = (site_[siteid]).getY();
  int qz = (site_[siteid]).getZ();

  int qxlo = qx - boxhalf;
  int qxhi = qx + boxhalf;
  int qylo = qy - boxhalf;
  int qyhi = qy + boxhalf;
  int qzlo = qz - boxhalf;
  int qzhi = qz + boxhalf;

  /***
   *    Count the number of requisite sites in the
   *    3-D cube box using whatever boundaries are specified
   ***/

  for (int ix = qxlo; ix <= qxhi; ix++) {
    for (int iy = qylo; iy <= qyhi; iy++) {
      for (int iz = qzlo; iz <= qzhi; iz++) {

        /// Count if phase id only

        if (site_[getIndex(ix, iy, iz)].getMicroPhaseId() == phaseid) {
          nfound++;
        }
      }
    }
  }

  return nfound;
}

vector<int> Lattice::findDomainSizeDistribution(int phaseid,
                                                const int numsites,
                                                int maxsize, int cyc,
                                                int sortorder = 0) {
  // if sortorder is 0 => sorting in descending order

  int domainsize = 0;
  int ix, iy, iz;
  int qx, qy, qz;
  int qxlo, qxhi, qylo, qyhi, qzlo, qzhi;
  int sizeCateg, pos, stId;
  double rng;

  int boxhalf = maxsize / 2;
  int dim = maxsize;
  if (maxsize % 2 == 0) dim++;
  int totSize = (int)pow(dim,3) + 1;

  vector<int> init;
  vector<vector<int>> siteDomainSizeDistribution(totSize, init);
  vector<int> selectedSites;

  for (int i = 0; i < numsites_; i++) {
    if (site_[i].getMicroPhaseId() == phaseid) {

      qx = (site_[i]).getX();
      qy = (site_[i]).getY();
      qz = (site_[i]).getZ();

      qxlo = qx - boxhalf;
      qxhi = qx + boxhalf;
      qylo = qy - boxhalf;
      qyhi = qy + boxhalf;
      qzlo = qz - boxhalf;
      qzhi = qz + boxhalf;

      /***
       *    Count the number of requisite sites in the
       *    3-D cube box using whatever boundaries are specified
       ***/

      domainsize = 0;
      for (ix = qxlo; ix <= qxhi; ix++) {
        for (iy = qylo; iy <= qyhi; iy++) {
          for (iz = qzlo; iz <= qzhi; iz++) {

            /// Count if phase id only

            if (site_[getIndex(ix, iy, iz)].getMicroPhaseId() == phaseid) {
              domainsize++;
            }
          }
        }
      }

      siteDomainSizeDistribution[domainsize].push_back(i);
    }
  }

  int j;
  int contor = 0;
  int diff = numsites;
  if (sortorder == 0) { //emptyPorosity
    j = totSize;
    while (diff > 0) {
      sizeCateg = siteDomainSizeDistribution[j - 1].size();
      if (diff < sizeCateg) {
        for (int k = 0; k < diff; k++) {
          rng = callRNG();
          pos = (int)(rng * (sizeCateg - 1));
          stId = siteDomainSizeDistribution[j - 1][pos];
          selectedSites.push_back(stId);
          siteDomainSizeDistribution[j - 1][pos] =
              siteDomainSizeDistribution[j - 1][sizeCateg - 1];

          sizeCateg--;
          contor++;
        }
      } else {
        if (sizeCateg > 0) {
          selectedSites.insert(selectedSites.end(),
                               siteDomainSizeDistribution[j - 1].begin(),
              siteDomainSizeDistribution[j - 1].end());
          contor += sizeCateg;
        }
        j--;
      }
      diff = numsites - contor;
      if (j == 0) break;
    }
  } else { // fillPorosity
    j = 1;
    while (diff > 0) {
      sizeCateg = siteDomainSizeDistribution[j].size();
      if (diff < sizeCateg) {
        for (int k = 0; k < diff; k++) {
          rng = callRNG();
          pos = (int)(rng * (sizeCateg - 1));
          stId = siteDomainSizeDistribution[j][pos];
          selectedSites.push_back(stId);
          siteDomainSizeDistribution[j][pos] =
              siteDomainSizeDistribution[j][sizeCateg - 1];

          sizeCateg--;
          contor++;
        }
      } else {
        if (sizeCateg > 0) {
          selectedSites.insert(selectedSites.end(),
                               siteDomainSizeDistribution[j].begin(),
              siteDomainSizeDistribution[j].end());
          contor += sizeCateg;
        }
        j++;
      }
      diff = numsites - contor;
      if (j > totSize) break;
    }
  }

  return (selectedSites);
}

void Lattice::findIsolatedClusters(void) {
  // voxels without contact with electrolyte i.e. voxels having low probability
  // to dissolve in this step
  // find "isolated" clusters (voxels without contact with electrolyte i.e.
  // voxels having low probability to dissolve in this step) these voxels will
  // not be send to GEMS computing vfrac we use canDissolve vector instead
  // count_ vector
  vector<int> canDissolve;
  int canDissLast;
  bool siteDiss;
  int i, j;

  for (i = 0; i < numMicroPhases_; i++) {
    canDissolve.push_back(0);
    if (i >= FIRST_SOLID) {
      canDissLast = -1;
      while (canDissLast != canDissolve[i]) {
        canDissLast = canDissolve[i];
        for (j = 0; j < numsites_; j++) {
          if (site_[j].getMicroPhaseId() == i) {
            if (site_[j].getVisit() == 0) {
              if (site_[j].getWmc() > 0) {
                canDissolve[i]++;
                site_[j].setVisit(1);
                for (int jj = 0; jj < NN_NNN; jj++) {
                  if (site_[j].nb(jj)->getMicroPhaseId() == i &&
                      site_[j].nb(jj)->getVisit() == 0) {
                    site_[j].nb(jj)->setVisit(1);
                    canDissolve[i]++;
                  }
                }
              } else {
                siteDiss = false;
                for (int jj = 0; jj < NN_NNN; jj++) {
                  if (site_[j].nb(jj)->getMicroPhaseId() == i &&
                      site_[j].nb(jj)->getVisit() == 1) {
                    canDissolve[i]++;
                    site_[j].setVisit(1);
                    siteDiss = true;
                    break;
                  }
                }
                if (siteDiss) {
                  for (int jj = 0; jj < NN_NNN; jj++) {
                    if (site_[j].nb(jj)->getMicroPhaseId() == i &&
                        site_[j].nb(jj)->getVisit() == 0) {
                      site_[j].nb(jj)->setVisit(1);
                      canDissolve[i]++;
                    }
                  }
                }
              }
            }
          }
        } // for
      } // while
    } // if
  } // for

  for (j = 0; j < numsites_; j++) {
    site_[j].setVisit(0);
  }

  cout << endl
       << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
       << endl;
  for (i = 0; i < numMicroPhases_; i++) {
    cout << "     i = " << i << "      canDissolve[" << i
         << "] = " << canDissolve[i] << endl;
  }
  cout << endl
       << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
       << endl;

  return;
}

void Lattice::populateElementData(void) {
  chemElem elem;

  elem.z = 1;
  elem.symb = "H";
  elem.mass = 1.0080;
  cfgElem_.push_back(elem);

  elem.z = 2;
  elem.symb = "He";
  elem.mass = 4.0026;
  cfgElem_.push_back(elem);

  elem.z = 3;
  elem.symb = "Li";
  elem.mass = 6.94;
  cfgElem_.push_back(elem);

  elem.z = 4;
  elem.symb = "Be";
  elem.mass = 9.0122;
  cfgElem_.push_back(elem);

  elem.z = 5;
  elem.symb = "B";
  elem.mass = 10.81;
  cfgElem_.push_back(elem);

  elem.z = 6;
  elem.symb = "C";
  elem.mass = 12.011;
  cfgElem_.push_back(elem);

  elem.z = 7;
  elem.symb = "N";
  elem.mass = 14.007;
  cfgElem_.push_back(elem);

  elem.z = 8;
  elem.symb = "O";
  elem.mass = 15.999;
  cfgElem_.push_back(elem);

  elem.z = 9;
  elem.symb = "F";
  elem.mass = 18.998;
  cfgElem_.push_back(elem);

  elem.z = 10;
  elem.symb = "Ne";
  elem.mass = 20.180;
  cfgElem_.push_back(elem);

  elem.z = 11;
  elem.symb = "Na";
  elem.mass = 22.990;
  cfgElem_.push_back(elem);

  elem.z = 12;
  elem.symb = "Mg";
  elem.mass = 24.305;
  cfgElem_.push_back(elem);

  elem.z = 13;
  elem.symb = "Al";
  elem.mass = 26.982;
  cfgElem_.push_back(elem);

  elem.z = 14;
  elem.symb = "Si";
  elem.mass = 28.085;
  cfgElem_.push_back(elem);

  elem.z = 15;
  elem.symb = "P";
  elem.mass = 30.974;
  cfgElem_.push_back(elem);

  elem.z = 16;
  elem.symb = "S";
  elem.mass = 32.06;
  cfgElem_.push_back(elem);

  elem.z = 17;
  elem.symb = "Cl";
  elem.mass = 35.45;
  cfgElem_.push_back(elem);

  elem.z = 18;
  elem.symb = "Ar";
  elem.mass = 39.95;
  cfgElem_.push_back(elem);

  elem.z = 19;
  elem.symb = "K";
  elem.mass = 39.098;
  cfgElem_.push_back(elem);

  elem.z = 20;
  elem.symb = "Ca";
  elem.mass = 40.078;
  cfgElem_.push_back(elem);

  elem.z = 21;
  elem.symb = "Sc";
  elem.mass = 44.956;
  cfgElem_.push_back(elem);

  elem.z = 22;
  elem.symb = "Ti";
  elem.mass = 47.867;
  cfgElem_.push_back(elem);

  elem.z = 23;
  elem.symb = "V";
  elem.mass = 50.942;
  cfgElem_.push_back(elem);

  elem.z = 24;
  elem.symb = "Cr";
  elem.mass = 51.996;
  cfgElem_.push_back(elem);

  elem.z = 25;
  elem.symb = "Mn";
  elem.mass = 54.938;
  cfgElem_.push_back(elem);

  elem.z = 26;
  elem.symb = "Fe";
  elem.mass = 55.845;
  cfgElem_.push_back(elem);

  elem.z = 27;
  elem.symb = "Co";
  elem.mass = 58.933;
  cfgElem_.push_back(elem);

  elem.z = 28;
  elem.symb = "Ni";
  elem.mass = 58.693;
  cfgElem_.push_back(elem);

  elem.z = 29;
  elem.symb = "Cu";
  elem.mass = 63.546;
  cfgElem_.push_back(elem);

  elem.z = 30;
  elem.symb = "Zn";
  elem.mass = 65.38;
  cfgElem_.push_back(elem);

  elem.z = 31;
  elem.symb = "Ga";
  elem.mass = 69.723;
  cfgElem_.push_back(elem);

  elem.z = 32;
  elem.symb = "Ge";
  elem.mass = 72.630;
  cfgElem_.push_back(elem);

  elem.z = 33;
  elem.symb = "As";
  elem.mass = 74.922;
  cfgElem_.push_back(elem);

  elem.z = 34;
  elem.symb = "Se";
  elem.mass = 78.971;
  cfgElem_.push_back(elem);

  elem.z = 35;
  elem.symb = "Br";
  elem.mass = 79.904;
  cfgElem_.push_back(elem);

  elem.z = 36;
  elem.symb = "Kr";
  elem.mass = 83.798;
  cfgElem_.push_back(elem);

  elem.z = 37;
  elem.symb = "Rb";
  elem.mass = 85.468;
  cfgElem_.push_back(elem);

  elem.z = 38;
  elem.symb = "Sr";
  elem.mass = 87.62;
  cfgElem_.push_back(elem);

  elem.z = 39;
  elem.symb = "Y";
  elem.mass = 88.906;
  cfgElem_.push_back(elem);

  elem.z = 40;
  elem.symb = "Zr";
  elem.mass = 91.224;
  cfgElem_.push_back(elem);

  elem.z = 41;
  elem.symb = "Nb";
  elem.mass = 92.906;
  cfgElem_.push_back(elem);

  elem.z = 42;
  elem.symb = "Mo";
  elem.mass = 95.95;
  cfgElem_.push_back(elem);

  elem.z = 43;
  elem.symb = "Tc";
  elem.mass = 97.0;
  cfgElem_.push_back(elem);

  elem.z = 44;
  elem.symb = "Ru";
  elem.mass = 101.07;
  cfgElem_.push_back(elem);

  elem.z = 45;
  elem.symb = "Rh";
  elem.mass = 102.91;
  cfgElem_.push_back(elem);

  elem.z = 46;
  elem.symb = "Pd";
  elem.mass = 106.42;
  cfgElem_.push_back(elem);

  elem.z = 47;
  elem.symb = "Ag";
  elem.mass = 107.87;
  cfgElem_.push_back(elem);

  elem.z = 48;
  elem.symb = "Cd";
  elem.mass = 112.41;
  cfgElem_.push_back(elem);

  elem.z = 49;
  elem.symb = "In";
  elem.mass = 114.82;
  cfgElem_.push_back(elem);

  elem.z = 50;
  elem.symb = "Sn";
  elem.mass = 118.71;
  cfgElem_.push_back(elem);

  elem.z = 51;
  elem.symb = "Sb";
  elem.mass = 121.76;
  cfgElem_.push_back(elem);

  elem.z = 52;
  elem.symb = "Te";
  elem.mass = 127.60;
  cfgElem_.push_back(elem);

  elem.z = 53;
  elem.symb = "I";
  elem.mass = 126.90;
  cfgElem_.push_back(elem);

  elem.z = 54;
  elem.symb = "Xe";
  elem.mass = 131.29;
  cfgElem_.push_back(elem);

  elem.z = 55;
  elem.symb = "Cs";
  elem.mass = 132.91;
  cfgElem_.push_back(elem);

  elem.z = 56;
  elem.symb = "Ba";
  elem.mass = 137.33;
  cfgElem_.push_back(elem);

  elem.z = 57;
  elem.symb = "La";
  elem.mass = 138.91;
  cfgElem_.push_back(elem);

  elem.z = 58;
  elem.symb = "Ce";
  elem.mass = 140.12;
  cfgElem_.push_back(elem);

  elem.z = 59;
  elem.symb = "Pr";
  elem.mass = 140.91;
  cfgElem_.push_back(elem);

  elem.z = 60;
  elem.symb = "Nd";
  elem.mass = 144.24;
  cfgElem_.push_back(elem);

  elem.z = 61;
  elem.symb = "Pm";
  elem.mass = 145.0;
  cfgElem_.push_back(elem);

  elem.z = 62;
  elem.symb = "Sm";
  elem.mass = 150.36;
  cfgElem_.push_back(elem);

  elem.z = 63;
  elem.symb = "Eu";
  elem.mass = 151.96;
  cfgElem_.push_back(elem);

  elem.z = 64;
  elem.symb = "Gd";
  elem.mass = 157.25;
  cfgElem_.push_back(elem);

  elem.z = 65;
  elem.symb = "Tb";
  elem.mass = 158.93;
  cfgElem_.push_back(elem);

  elem.z = 66;
  elem.symb = "Dy";
  elem.mass = 162.50;
  cfgElem_.push_back(elem);

  elem.z = 67;
  elem.symb = "Ho";
  elem.mass = 164.93;
  cfgElem_.push_back(elem);

  elem.z = 68;
  elem.symb = "Er";
  elem.mass = 167.26;
  cfgElem_.push_back(elem);

  elem.z = 69;
  elem.symb = "Tm";
  elem.mass = 168.93;
  cfgElem_.push_back(elem);

  elem.z = 70;
  elem.symb = "Yb";
  elem.mass = 173.05;
  cfgElem_.push_back(elem);

  elem.z = 71;
  elem.symb = "Lu";
  elem.mass = 174.97;
  cfgElem_.push_back(elem);

  elem.z = 72;
  elem.symb = "Hf";
  elem.mass = 178.49;
  cfgElem_.push_back(elem);

  elem.z = 73;
  elem.symb = "Ta";
  elem.mass = 180.95;
  cfgElem_.push_back(elem);

  elem.z = 74;
  elem.symb = "W";
  elem.mass = 183.84;
  cfgElem_.push_back(elem);

  elem.z = 75;
  elem.symb = "Re";
  elem.mass = 186.21;
  cfgElem_.push_back(elem);

  elem.z = 76;
  elem.symb = "Os";
  elem.mass = 190.23;
  cfgElem_.push_back(elem);

  elem.z = 77;
  elem.symb = "Ir";
  elem.mass = 192.22;
  cfgElem_.push_back(elem);

  elem.z = 78;
  elem.symb = "Pt";
  elem.mass = 195.08;
  cfgElem_.push_back(elem);

  elem.z = 79;
  elem.symb = "Au";
  elem.mass = 196.97;
  cfgElem_.push_back(elem);

  elem.z = 80;
  elem.symb = "Hg";
  elem.mass = 200.59;
  cfgElem_.push_back(elem);

  elem.z = 81;
  elem.symb = "Tl";
  elem.mass = 204.38;
  cfgElem_.push_back(elem);

  elem.z = 82;
  elem.symb = "Pb";
  elem.mass = 207.2;
  cfgElem_.push_back(elem);

  elem.z = 83;
  elem.symb = "Bi";
  elem.mass = 208.98;
  cfgElem_.push_back(elem);

  elem.z = 84;
  elem.symb = "Po";
  elem.mass = 209.0;
  cfgElem_.push_back(elem);

  elem.z = 85;
  elem.symb = "At";
  elem.mass = 210.0;
  cfgElem_.push_back(elem);

  elem.z = 86;
  elem.symb = "Rn";
  elem.mass = 222.0;
  cfgElem_.push_back(elem);

  elem.z = 87;
  elem.symb = "Fr";
  elem.mass = 223.0;
  cfgElem_.push_back(elem);

  elem.z = 88;
  elem.symb = "Ra";
  elem.mass = 226.0;
  cfgElem_.push_back(elem);

  elem.z = 89;
  elem.symb = "Ac";
  elem.mass = 227.0;
  cfgElem_.push_back(elem);

  elem.z = 90;
  elem.symb = "Th";
  elem.mass = 232.04;
  cfgElem_.push_back(elem);

  elem.z = 91;
  elem.symb = "Pa";
  elem.mass = 231.04;
  cfgElem_.push_back(elem);

  elem.z = 92;
  elem.symb = "U";
  elem.mass = 238.03;
  cfgElem_.push_back(elem);

  elem.z = 93;
  elem.symb = "Np";
  elem.mass = 237.0;
  cfgElem_.push_back(elem);

  elem.z = 94;
  elem.symb = "Pu";
  elem.mass = 244.0;
  cfgElem_.push_back(elem);

  elem.z = 95;
  elem.symb = "Am";
  elem.mass = 243.0;
  cfgElem_.push_back(elem);

  elem.z = 96;
  elem.symb = "Cm";
  elem.mass = 247.0;
  cfgElem_.push_back(elem);

  elem.z = 97;
  elem.symb = "Bk";
  elem.mass = 247.0;
  cfgElem_.push_back(elem);

  elem.z = 98;
  elem.symb = "Cf";
  elem.mass = 251.0;
  cfgElem_.push_back(elem);

  elem.z = 99;
  elem.symb = "Es";
  elem.mass = 252.0;
  cfgElem_.push_back(elem);

  elem.z = 100;
  elem.symb = "Fm";
  elem.mass = 257.0;
  cfgElem_.push_back(elem);

  elem.z = 101;
  elem.symb = "Md";
  elem.mass = 258.0;
  cfgElem_.push_back(elem);

  elem.z = 102;
  elem.symb = "No";
  elem.mass = 259.0;
  cfgElem_.push_back(elem);

  elem.z = 103;
  elem.symb = "Lr";
  elem.mass = 266.0;
  cfgElem_.push_back(elem);

  elem.z = 104;
  elem.symb = "Rf";
  elem.mass = 267.0;
  cfgElem_.push_back(elem);

  elem.z = 105;
  elem.symb = "Db";
  elem.mass = 268.0;
  cfgElem_.push_back(elem);

  elem.z = 106;
  elem.symb = "Sg";
  elem.mass = 269.0;
  cfgElem_.push_back(elem);

  elem.z = 107;
  elem.symb = "Bh";
  elem.mass = 270.0;
  cfgElem_.push_back(elem);

  elem.z = 108;
  elem.symb = "Hs";
  elem.mass = 269.0;
  cfgElem_.push_back(elem);

  elem.z = 109;
  elem.symb = "Mt";
  elem.mass = 278.0;
  cfgElem_.push_back(elem);

  elem.z = 110;
  elem.symb = "Ds";
  elem.mass = 281.0;
  cfgElem_.push_back(elem);

  elem.z = 111;
  elem.symb = "Rg";
  elem.mass = 282.0;
  cfgElem_.push_back(elem);

  elem.z = 112;
  elem.symb = "Cn";
  elem.mass = 285.0;
  cfgElem_.push_back(elem);

  elem.z = 113;
  elem.symb = "Nh";
  elem.mass = 286.0;
  cfgElem_.push_back(elem);

  elem.z = 114;
  elem.symb = "Fl";
  elem.mass = 289.0;
  cfgElem_.push_back(elem);

  elem.z = 115;
  elem.symb = "Mc";
  elem.mass = 290.0;
  cfgElem_.push_back(elem);

  elem.z = 116;
  elem.symb = "Lv";
  elem.mass = 293.0;
  cfgElem_.push_back(elem);

  elem.z = 117;
  elem.symb = "Ts";
  elem.mass = 294.0;
  cfgElem_.push_back(elem);

  elem.z = 118;
  elem.symb = "Og";
  elem.mass = 294.0;
  cfgElem_.push_back(elem);

  return;
}
