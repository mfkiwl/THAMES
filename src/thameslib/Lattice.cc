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

Lattice::Lattice(ChemicalSystem *cs, const string &fileName, const bool verbose,
                 const bool warning)
    : siteneighbors_(NN_NNN), chemSys_(cs) {
  unsigned int i, j, k;
  unsigned int ii;
  string buff;
  int xn, yn, zn;
  unsigned int idn;
  unsigned int pid;
  string msg;

  numMicroPhases_ = chemSys_->getNumMicroPhases();

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

    version_ = "2.0";
    double testres = 1.0;
    setResolution(testres);
    xdim_ = ydim_ = zdim_ = 100;
  }

  ///
  /// Print out the microstructure size and characteristics
  ///

  if (verbose_) {
    cout << "Lattice::Lattice Read microstructure file header..." << endl;
    cout << "Lattice::Lattice     Version = " << version_ << endl;
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

  try {
    createRNG();
  } catch (bad_alloc &ba) {
    cout << "Lattice constructor failed when allocating rg_";
    cout.flush();
    exit(1);
  }
  latticeRNGseed_ = -142234;
  numRNGcall_0_ = 0;
  numRNGcallLONGMAX_ = 0;
  setRNGseed(latticeRNGseed_);

  count_.clear();
  count_.resize(numMicroPhases_, 0);

  SI_.clear();
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
    count_.at(pid)++;
  }

  ///
  /// Done with the input microstructure file, so close it.
  ///

  in.close();

  ///
  /// Gather data from the periodic table for writing a .cfg file (alternative
  /// to .xyz)

  populateElementData();

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
    int totCount_ = 0;
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
  delete rg_;
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

      chemSys_->setInitScaledCementMass(cementMass * 100 / solidMass);

      totalSolidMass += pscaledMass;
      if (chemSys_->getCementComponent(microPhaseId))
        totalCementMass += pscaledMass;
      chemSys_->setMicroPhaseMass(microPhaseId, pscaledMass);

      chemSys_->setDCMoles(DCId, (pscaledMass / molarMass));

      chemSys_->setMicroPhaseMassDissolved(microPhaseId, 0.0);
    }
  }
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
            site_[k].setDissolutionSite(i);
            for (kk = 0; kk < NN_NNN; kk++) {
              if ((site_[k].nb(kk))->getMicroPhaseId() == ELECTROLYTEID) {
                gsite.push_back(site_[k].nb(kk));
                site_[k].nb(kk)->setGrowthSite(i);
              }
            }

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
    cout << "  " << i << "   " << chemSys_->getMicroPhaseName(i)
         << "   id: " << chemSys_->getMicroPhaseId(i)
         << "     count_ = " << count_[i] << "     dissInterfaceSize =  "
         << interface_[i].getDissolutionNumSites()
         << "     growInterfaceSize =  " << interface_[i].getGrowthNumSites()
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
  vector<int> numChange;
  numChange.resize(growPhaseIDVectSize, 0);
  int numChangeTot = 0;
  vector<int> numLeft = numSiteGrowVect; // numtotake
  int numLeftTot = 0;
  for (i = 0; i < growPhaseIDVectSize; i++) {
    numLeftTot += numSiteGrowVect[i];
  }
  if (numLeftTot == 0) {
    cout << endl << "     Lattice::growPhase error => numLeftTot = 0" << endl;
    cout << "     totalTRC/trc_g/bcl :  "
         << "   " << totalTRC << "/" << trc_g << "/" << bcl << endl;
    cout << endl << "     stop program" << endl;
    exit(0);
  }
  vector<int> nucleated;
  nucleated.resize(growPhaseIDVectSize, 0);

  Site *ste, *stenb;
  vector<Isite> isite;
  isite.clear();
  int isitePos;
  unsigned int pid;

  int phaseID;
  vector<int> dim_isite;
  dim_isite.resize(growPhaseIDVectSize, 0);

  int aff, affMin, affSum;
  double affSumDbl;
  int valAbs;
  double rng;

  string nameMicroPhaseTh;
  double steWmc, stenbWmc, dwmcval;

  struct growProb_ {
    unsigned int id; /**< The id of the corresponding Site */
    int affinity;    /**< The affinity for growth of a phase at the site */
    double prob; /**< The growth probability of a phase at this site (computed
                    according the affinity) */
    int posVect;
    int interfacePos;
  };

  // growth probabilities based on affinities
  vector<growProb_> growProbVect;
  growProbVect.clear();
  growProb_ growProb;
  int posVect;

  int siteID, sitePhID;
  for (i = 0; i < growPhaseIDVectSize; i++) {
    phaseID = growPhaseIDVect[i];
    isite = interface_[phaseID].getGrowthSites();
    dim_isite[i] = isite.size();
    for (jj = 0; jj < dim_isite[i]; jj++) {
      siteID = isite[jj].getId();
      sitePhID = site_[siteID].getMicroPhaseId();
      if (sitePhID == 0 || sitePhID > 1) {
        cout << endl
             << "     Lattice::growPhase => error for the interface of phase "
                "phaseID = "
             << phaseID << endl;
        cout << "     siteID = " << siteID << "    sitePhID = " << sitePhID
             << endl;
        cout << "     totalTRC/trc_g/bcl :  "
             << "   " << totalTRC << "/" << trc_g << "/" << bcl << endl;
        cout << endl << "     stop program" << endl;
        exit(1);
      }
      growProb.id = isite[jj].getId();
      growProb.affinity = isite[jj].getAffinity();
      growProb.posVect = i;
      growProb.interfacePos = jj;

      growProbVect.push_back(growProb);
    }
  }
  int growProbVectSize = growProbVect.size();

  cout << endl
       << "   Lattice::growPhase GROW_INI totalTRC/trc_g/bcl " << totalTRC
       << "/" << trc_g << "/" << bcl << endl;
  cout << "     GROW_INI growPhaseIDVectSize = " << growPhaseIDVectSize
       << "   growProbVectSize = " << growProbVectSize
       << "   numLeftTot = " << numLeftTot
       << "   numChangeTot = " << numChangeTot << endl;
  for (i = 0; i < growPhaseIDVectSize; i++) {
    phaseID = growPhaseIDVect[i];
    // isite = interface_[phaseID].getDissolutionSites();
    // dim_isite = isite.size();
    cout << "       GROW_INI for i = " << i
         << "  => phaseID phaseName count_ dim_isite numleft numchange  :  "
         << phaseID << "   " << growPhNameVect[i] << "   " << count_[phaseID]
         << "   " << dim_isite[i] << "   " << numLeft[i] << "   "
         << numChange[i] << endl;
  }
  cout.flush();

  int interfacePos;
  vector<bool> writeFirst;
  writeFirst.resize(growPhaseIDVectSize, false);

  // affAllPos = false;
  while ((numLeftTot > 0) && (growProbVectSize >= 1)) {
    try {
      bcl++;
      // growth probabilities based on affinities

      // affAllPos - works with modified affinities in affinity_ table from
      // ChemicalSystem (all positives!)
      affSum = 0;
      for (j = 0; j < growProbVectSize; j++) {
        aff = growProbVect[j].affinity;
        affSum += aff;
      }
      // calc probabilities & choose a site
      rng = callRNG();
      if (affSum != 0) {
        affSumDbl = affSum;
        growProbVect[0].prob = growProbVect[0].affinity / affSumDbl;
        if (rng <= growProbVect[0].prob) {
          isitePos = 0;
        } else {
          for (isitePos = 1; isitePos < growProbVectSize; isitePos++) {
            growProbVect[isitePos].prob =
                growProbVect[isitePos - 1].prob +
                growProbVect[isitePos].affinity / affSumDbl;
            if (rng <= growProbVect[isitePos].prob)
              break;
          }
        }
      } else {
        isitePos = (int)(rng * growProbVectSize);
      }

      ste = &site_[growProbVect[isitePos].id];
      pid = ste->getMicroPhaseId(); // always ELECTROLYTEID !!
      posVect = growProbVect[isitePos].posVect;
      phaseID = growPhaseIDVect[posVect];
      interfacePos = growProbVect[isitePos].interfacePos;

      if (pid != ELECTROLYTEID) {
        cout << endl
             << "Lattice::growPhase error: phaseid != ELECTROLYTEID" << endl;
        cout << "  phaseid totalTRC/trc_g/bcl count_ numLeftTot numChangeTot  "
                ":  "
             << phaseID << "   " << totalTRC << "/" << trc_g << "/" << bcl
             << "   " << count_[phaseID] << "   " << numLeftTot << "   "
             << numChangeTot << endl;
        cout << " posVect interfacePos pid : " << posVect << "   "
             << interfacePos << "   " << pid << endl;
        cout << " ste.id numLeft numChange : " << ste->getId() << "   "
             << numLeft[posVect] << "   " << numChange[posVect] << endl;
        cout << "STOP" << endl;
        exit(0);
      }

      double wmcIni, wmcEnd;
      wmcIni = ste->getWmc0();

      // if (pid == ELECTROLYTEID) {
      // removeGrowthSite(ste, phaseid);
      removeGrowthSite_0(ste, phaseID, interfacePos);
      setMicroPhaseId(ste, phaseID, "grow");

      vector<unsigned int> plist = ste->getGrowthPhases();
      int plsize = plist.size();
      for (j = 0; j < plsize; j++) {
        removeGrowthSite_1(ste, plist[j]);
      }

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
      // if (steWmc < 0.0 || steWmc > 19.0){
      //    cout << endl <<"error Lattice::growPhase for steWmc - ste pid steWmc
      //    phaseid trc trc_g bcl : "
      //          << ste->getId() << "   "
      //          << pid << "   " << steWmc << "   " << phaseid << "   "
      //         << trc << "   " << trc_g << "   " << bcl << endl; cout.flush();
      //     cout << endl << "nn & nnn siteId phaseId siteWmc :" << endl;
      //     for (int i = 0; i < numMicroPhases_; i++){
      //         cout << "   " << i << "   " << (ste->nb(i))->getId() << "   "
      //              << (ste->nb(i))->getMicroPhaseId() << "   " <<
      //              (ste->nb(i))->getWmc() << endl;
      //     }
      //     cout << endl << "stop program" << endl;
      //     exit(1);
      // }
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

      for (j = 0; j < NN_NNN;
           j++) { // for (j = 0; j < NUM_NEAREST_NEIGHBORS; j++) {
        stenb = ste->nb(j);
        stenb->dWmc(dwmcval);
        stenbWmc = stenb->getWmc();
        // if (stenbWmc < 0.0 || stenbWmc > 19.0){
        //     cout << endl <<"error Lattice::growPhase fornumLeft[i] stenbWmc -
        //     ste pid steWmc phaseid trc trc_g bcl stenb stenb_pid stenbWmc j :
        //     "
        //          << ste->getId() << "   "
        //          << pid << "   " << ste->getWmc() << "   " << phaseid << " "
        //          << trc << "   " << trc_g << "   " << bcl << "   " <<
        //          stenb->getId() << "   "
        //          << stenb->getMicroPhaseId() << "   " << stenbWmc << "   " <<
        //          j << endl; cout.flush();
        //     cout << endl << "stop program" << endl;
        //     exit(1);
        // }
        if (stenb->getMicroPhaseId() == ELECTROLYTEID) {
          addGrowthSite(stenb, phaseID);

          if (j < NUM_NEAREST_NEIGHBORS) {
            for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases_;
                 phaseTmpl++) {
              if (chemSys_->isGrowthTemplate(phaseTmpl, phaseID)) {
                // addGrowthSite(stenb, nbgrowthtemp[jj]);
                addGrowthSite(stenb, phaseTmpl);
              }
            }
          }

        } else if ((stenbWmc == 0.0) &&
                   (stenb->getMicroPhaseId() > ELECTROLYTEID)) {

          // removeDissolutionSite(stenb, stenb->getMicroPhaseId());
          removeDissolutionSite_grow(stenb, stenb->getMicroPhaseId());
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
        cout << " posVect interfacePos pid : " << posVect << "   "
             << interfacePos << "   " << pid << endl;
        cout << " ste.id numLeft numChange : " << ste->getId() << "   "
             << numLeft[posVect] << "   " << numChange[posVect] << endl;
        cout << "STOP" << endl;
        exit(0);
      }

    } catch (out_of_range &oor) {
      EOBException ex("Lattice", "growPhase", "site_", site_.size(), j);
      ex.printException();

      cout << endl << "Lattice::growPhase error" << endl;
      cout << endl
           << "totalTRC trc_g bcl numLeftTot numChangeTot  :  " << totalTRC
           << "   " << trc_g << "   " << bcl << "   " << numLeftTot << "   "
           << numChangeTot << endl;
      cout << endl
           << "steId pid growProbVectSize :  " << ste->getId() << "   " << pid
           << "   " << growProbVectSize << endl;
      cout.flush();
      for (i = 0; i < growPhaseIDVectSize; i++) {
        phaseID = growPhaseIDVect[i];
        isite = interface_[phaseID].getGrowthSites();
        dim_isite[i] = isite.size();
        cout << "        phaseid count_ dim_isite numleft numchange  :  "
             << phaseID << "   " << "   " << count_[phaseID] << "   "
             << dim_isite[i] << "   " << numLeft[i] << "   " << numChange[i]
             << endl;
      }
      cout << "stop program" << endl;
      exit(0);
    }

    // growth probabilities based on affinities
    growProbVect.clear();
    for (i = 0; i < growPhaseIDVectSize; i++) {
      if (numLeft[i] > 0) {

        phaseID = growPhaseIDVect[i];
        isite = interface_[phaseID].getGrowthSites();
        dim_isite[i] = isite.size();
        if (dim_isite[i] == 0) { // nucleation
          cout << endl
               << "   Lattice::growPhase - need nucleation for phaseID = "
               << phaseID << endl;
          cout << "     interface dimension dim_site[" << i
               << "] = " << dim_isite[i] << " while numLeft[" << i
               << "] = " << numLeft[i] << endl;
          cout << "     => for this microPhase (" << growPhNameVect[i]
               << ") a new interface having " << numLeft[i]
               << " sites will be created! " << endl;

          try {
            nucleated[i] = createNewGrowthInterface(phaseID, numLeft[i]);
          } catch (MicrostructureException mex) {
            throw mex;
          }

          isite = interface_[phaseID].getGrowthSites();
          dim_isite[i] = isite.size();
        }
        for (jj = 0; jj < dim_isite[i]; jj++) {
          growProb.id = isite[jj].getId();
          growProb.affinity = isite[jj].getAffinity();
          growProb.posVect = i;
          growProb.interfacePos = jj;

          growProbVect.push_back(growProb);
        }
      } else {
        if (writeFirst[i]) {
          continue;
        } else {
          writeFirst[i] = true;
          cout << "   Lattice::growPhase GROW_END for i = " << i
               << "   totalTRC/trc_g/bcl " << totalTRC << "/" << trc_g << "/"
               << bcl << endl;
          cout << "     GROW_END growPhaseIDVectSize = " << growPhaseIDVectSize
               << "   growProbVectSize = " << growProbVectSize
               << "   numLeftTot = " << numLeftTot
               << "   numChangeTot = " << numChangeTot << endl;
          cout << "       GROW_END phaseid count_ dim_isite numleft numchange  "
                  ":  "
               << growPhaseIDVect[i] << "   " << "   "
               << count_[growPhaseIDVect[i]] << "   " << dim_isite[i] << "   "
               << numLeft[i] << "   " << numChange[i] << endl;
          cout.flush();
        }
      }
    }
    growProbVectSize = growProbVect.size();
  }

  // findInterfaces_check();
  numadded_G = numChangeTot;
  return (nucleated);
}

int Lattice::createNewGrowthInterface(int phaseID, int numLeft) {

  struct localStruct {
    int id;
    int aff;
    double prob;
  };
  localStruct un;
  vector<localStruct> watersites;
  int aff, affMin, affSum;
  vector<Site *> localNb;
  int sizeWS;
  int valAbs;
  double rng;
  int fSiteWS;
  double daffsum;
  int j, k;
  int dimNewInterface;

  if (numLeft <= count_[ELECTROLYTEID]) {
    affSum = 0;
    affMin = 1000000;
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
        if (aff < affMin)
          affMin = aff;
      }
    }
    sizeWS = watersites.size();
    cout
        << endl
        << "     Lattice::createNewGrowthInterface -ini- sizeWS/affSum/affMin: "
        << sizeWS << " / " << affSum << " / " << affMin << endl;

    dimNewInterface = 0;
    while (dimNewInterface < numLeft) {

      rng = callRNG();

      if (affSum != 0) {
        if (affMin < 0) {
          valAbs = abs(affMin);
          affSum += (sizeWS * valAbs);
          if (affSum != 0) {
            daffsum = affSum;
            watersites[0].prob = (watersites[0].aff + valAbs) / daffsum;
            for (k = 1; k < sizeWS; k++) {
              watersites[k].prob = watersites[k - 1].prob +
                                   (watersites[k].aff + valAbs) / daffsum;
            }
          } else {
            fSiteWS = (int)(rng * sizeWS);
          }
        } else {
          daffsum = affSum;
          watersites[0].prob = (watersites[0].aff) / daffsum;
          for (k = 1; k < sizeWS; k++) {
            watersites[k].prob =
                watersites[k - 1].prob + (watersites[k].aff) / daffsum;
          }
        }

        for (fSiteWS = 0; fSiteWS < sizeWS; fSiteWS++) {
          if (rng <= watersites[fSiteWS].prob)
            break;
        }

      } else {

        fSiteWS = (int)(rng * sizeWS);
      }

      site_[watersites[fSiteWS].id].setGrowthSite(phaseID);
      interface_[phaseID].addGrowthSite_newInterface(watersites[fSiteWS].id,
                                                     watersites[fSiteWS].aff);
      dimNewInterface++;

      watersites[fSiteWS] = watersites[sizeWS - 1];
      watersites.pop_back();
      sizeWS--;

      affSum = 0;
      affMin = 1000000;
      for (k = 0; k < sizeWS; ++k) {
        aff = watersites[k].aff;
        affSum += aff;
        if (aff < affMin)
          affMin = aff;
      }
    }
  } else {
    cout << endl
         << "     Lattice::createNewGrowthInterface => requested interface "
            "size (numLeft) is larger than the electrolyte voxel number!"
         << endl;
    cout << "     numLeft > count_[ELECTROLYTEID] : " << numLeft << " > "
         << count_[ELECTROLYTEID] << endl;
    cout
        << endl
        << "     There is no room to create a new interface => exit the program"
        << endl;
    bool is_Error = false;
    throw MicrostructureException("Lattice", "createNewGrowthInterface",
                                  "no room for a new interface", is_Error);
    // exit(0);
  }

  int siteID;
  for (int i = 0; i < dimNewInterface; i++) {
    siteID = interface_[phaseID].getGrowthSitesId(i);
    if (site_[siteID].getMicroPhaseId() != ELECTROLYTEID) {
      cout << endl
           << "    Lattice::createNewGrowthInterface error: microPhaseId() for "
              "site_["
           << siteID << "] = " << site_[siteID].getMicroPhaseId() << endl;
      cout << "    Lattice::createNewGrowthInterface i/dimNewInterface: " << i
           << " / " << dimNewInterface << endl;
      cout << "    Lattice::createNewGrowthInterface phaseID/numLeft: "
           << phaseID << " / " << numLeft << endl;
      cout << "    Lattice::createNewGrowthInterface dimNewInterface: "
           << interface_[phaseID].getGrowthSites().size() << endl;
      cout << "    stop" << endl;
      exit(1);
      //} else {
      //    cout << "     i: " << i << "   siteID: " << siteID << "
      //    microPhaseId: " << site_[siteID].getMicroPhaseId() << endl;
    }
  }
  cout << "     Lattice::createNewGrowthInterface -end- dimNewInterface: "
       << dimNewInterface << endl
       << endl;

  return dimNewInterface;
}

int Lattice::dissolvePhase(vector<int> dissPhaseIDVect,
                           vector<int> numSiteDissVect,
                           vector<string> dissPhNameVect, int &numadded_D,
                           int totalTRC) {

  //*** controll
  int bcl = 0;
  int static trc_d;
  trc_d++;

  unsigned int i, ii, jj;

  int dissPhaseIDVectSize = dissPhaseIDVect.size();
  vector<int> numChange;
  numChange.resize(dissPhaseIDVectSize, 0);
  int numChangeTot = 0;
  vector<int> numLeft = numSiteDissVect; // numtotake
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

  double dwmcval;
  Site *ste, *stenb;
  vector<Isite> isite;
  isite.clear();

  int isitePos;
  unsigned int pid;

  vector<unsigned int> growth_local;
  int grLocSize;
  int nbid;
  int nb_id, nb_pid; // for nn & nnn
  bool phaseid_exist;

  struct dissProb {
    int id;
    int posVect;
    int nbElectr;
    double instab;
    double prob;
    double prob_01;
    int interfacePos;
  };
  vector<dissProb> dissProbVect;
  dissProbVect.clear();
  dissProb dissProb_;

  int phaseID;
  vector<int> dim_isite;
  dim_isite.resize(dissPhaseIDVectSize, 0);
  int stId;
  double sumWmc, stWmc;
  sumWmc = 0.0;
  // dissolution probabilities based on wmc
  for (i = 0; i < dissPhaseIDVectSize; i++) {
    phaseID = dissPhaseIDVect[i];
    isite = interface_[phaseID].getDissolutionSites();
    dim_isite[i] = isite.size();
    for (jj = 0; jj < dim_isite[i]; jj++) {
      stId = isite[jj].getId();
      stWmc = site_[stId].getWmc();
      sumWmc += stWmc;

      dissProb_.id = stId;
      dissProb_.posVect = i;
      dissProb_.instab = stWmc;
      dissProb_.interfacePos = jj;

      dissProbVect.push_back(dissProb_);
    }
  }
  int dissProbVectSize = dissProbVect.size();
  int dissProbVectSize_0 = dissProbVectSize;

  double rng;

  cout << endl
       << "   Lattice::dissolvePhase DISS_INI totalTRC/trc_d/bcl " << totalTRC
       << "/" << trc_d << "/" << bcl << endl;
  cout << "     DISS_INI dissPhaseIDVectSize = " << dissPhaseIDVectSize
       << "   dissProbVectSize = " << dissProbVectSize
       << "   numLeftTot = " << numLeftTot
       << "   numChangeTot = " << numChangeTot << endl;
  for (i = 0; i < dissPhaseIDVectSize; i++) {
    phaseID = dissPhaseIDVect[i];
    // isite = interface_[phaseID].getDissolutionSites();
    // dim_isite = isite.size();
    cout << "       DISS_INI for i = " << i
         << "  => phaseID phaseName count_ dim_isite numleft numchange  :  "
         << phaseID << "   " << dissPhNameVect[i] << "   " << count_[phaseID]
         << "   " << dim_isite[i] << "   " << numLeft[i] << "   "
         << numChange[i] << endl;
  }
  cout.flush();

  int posVect, interfacePos;
  bool goodEnd = true;
  int callGEM = -1;
  vector<bool> writeFirst;
  writeFirst.resize(dissPhaseIDVectSize, false);

  while ((numLeftTot > 0) && (dissProbVectSize >= 1)) {
    try {
      bcl++;

      // dissolution probabilities based on wmc
      rng = callRNG();
      dissProbVect[0].prob = (dissProbVect[0].instab) / sumWmc;
      dissProbVect[0].prob_01 = dissProbVect[0].prob;
      if (rng <= dissProbVect[0].prob_01) {
        isitePos = 0;
      } else {
        for (isitePos = 1; isitePos < dissProbVectSize; isitePos++) {
          dissProbVect[isitePos].prob =
              (dissProbVect[isitePos].instab) / sumWmc;
          dissProbVect[isitePos].prob_01 =
              dissProbVect[isitePos - 1].prob_01 + dissProbVect[isitePos].prob;
          if (rng <= dissProbVect[isitePos].prob_01)
            break;
        }
      }

      ste = &site_.at(dissProbVect[isitePos].id);
      pid = ste->getMicroPhaseId(); // intrebare pid diff phaseid ???
      posVect = dissProbVect[isitePos].posVect;
      interfacePos = dissProbVect[isitePos].interfacePos;

      double wmcIni, wmcEnd;
      wmcIni = ste->getWmc0();

      removeDissolutionSite_diss(ste, pid, interfacePos);
      setMicroPhaseId(ste, ELECTROLYTEID, "diss");

      ///
      /// Weighted mean curvature (wmc) is changed by the difference
      /// between the growing phase's porosity and the template's porosity.
      ///
      /// @todo Determine why the calculation works this way.
      ///

      wmcEnd = ELECTROLYTEID;
      ste->setWmc0(wmcEnd);

      dwmcval = wmcEnd - wmcIni;
      ste->dWmc(dwmcval);
      // stWmc = ste->getWmc();
      // if (stWmc < 0.0 || stWmc > 19.0){
      //     cout << endl <<"error Lattice::dissolvePhase for steWmc - ste pid
      //     steWmc phaseid trc trc_d bcl : "
      //          << ste->getId() << "   "
      //          << pid << "   " << stWmc << "   " << phaseid << "   "
      //          << trc << "   " << trc_d << "   " << bcl << endl;
      //          cout.flush();
      //     cout << endl << "nn & nnn siteId phaseId siteWmc :" << endl;
      //     for (int i = 0; i < numMicroPhases_; i++){
      //         cout << "   " << i << "   " << (ste->nb(i))->getId() << "   "
      //              << (ste->nb(i))->getMicroPhaseId() << "   " <<
      //              (ste->nb(i))->getWmc() << endl;
      //     }
      //     cout << endl << "stop program" << endl;
      //     exit(1);
      // }

      // dwmcval = chemSys_->getMicroPhasePorosity(ELECTROLYTEID) -
      //          chemSys_->getMicroPhasePorosity(pid);

      // ste->dWmc(dwmcval);

      for (i = 0; i < NN_NNN; i++) {
        stenb = ste->nb(i);
        stenb->dWmc(dwmcval);
        int nbpid = stenb->getMicroPhaseId();
        if (nbpid > ELECTROLYTEID) {

          ///
          /// Now that the site has been dissolved, it is eligible for growth
          /// later on, so we add it to the list of growth sites for certain
          /// phases.
          ///

          addDissolutionSite(stenb, nbpid);

          addGrowthSite(ste, nbpid);

          if (i < NUM_NEAREST_NEIGHBORS) {

            for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases_;
                 phaseTmpl++) {
              if (chemSys_->isGrowthTemplate(phaseTmpl, nbpid)) {

                addGrowthSite(ste, phaseTmpl);
              }
            }
          }
        } else if (nbpid ==
                   ELECTROLYTEID) { // local update the growthSites_ interfaces
                                    // for all phases in the system
          nbid = stenb->getId();
          growth_local = site_[nbid].getGrowthPhases();
          grLocSize = growth_local.size();

          for (ii = 0; ii < grLocSize; ii++) {
            phaseid_exist = false;
            for (jj = 0; jj < NN_NNN; jj++) {
              nb_pid = site_[nbid].nb(jj)->getMicroPhaseId();
              if (growth_local[ii] == nb_pid) {
                phaseid_exist = true;
                break;
              }
              if (jj < NUM_NEAREST_NEIGHBORS) {
                if (chemSys_->isGrowthTemplate(growth_local[ii], nb_pid)) {
                  phaseid_exist = true;
                  break;
                }
              }
            }
            if (phaseid_exist == false) { // verif removeGrowthSite_1 !!
              removeGrowthSite_1(stenb, growth_local[ii]);
            }
          }
        }
      }

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
           << "steId pid dissProbVectSize :  " << ste->getId() << "   " << pid
           << "   " << dissProbVectSize << endl;
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
    dissProbVect.clear();
    sumWmc = 0.;
    for (i = 0; i < dissPhaseIDVectSize; i++) {
      if (numLeft[i] > 0) {
        phaseID = dissPhaseIDVect[i];
        isite = interface_[phaseID].getDissolutionSites();
        dim_isite[i] = isite.size();
        if (dim_isite[i] == 0) {
          cout << "Lattice::dissolvePhase - anormal end:  dimInterface[" << i
               << "] = 0" << endl;
          cout << "    anormal end: phaseID[" << i << "] = " << phaseID
               << "  &  dim_site[" << i << "] = " << dim_isite[i]
               << " while numLeft[" << i << "] = " << numLeft[i] << endl;
          cout << "    anormal end: totalTRC/trc_d/bcl " << totalTRC << "/"
               << trc_d << "/" << bcl << endl;
          cout << "  need to call GEM after modifying DCLowerLimit_ for "
                  "microPhase "
               << phaseID << endl;
          cout.flush();

          callGEM = i;
          break;
        }
        for (jj = 0; jj < dim_isite[i]; jj++) {
          stId = isite[jj].getId();
          stWmc = site_[stId].getWmc();
          sumWmc += stWmc;

          dissProb_.id = stId;
          dissProb_.posVect = i;
          dissProb_.instab = stWmc;
          dissProb_.interfacePos = jj;

          dissProbVect.push_back(dissProb_);
        }
      } else {
        if (writeFirst[i]) {
          continue;
        } else {
          writeFirst[i] = true;

          cout << "   Lattice::dissolvePhase DISS_END for i = " << i
               << "   totalTRC/trc_d/bcl " << totalTRC << "/" << trc_d << "/"
               << bcl << endl;
          cout << "     DISS_END dissPhaseIDVectSize = " << dissPhaseIDVectSize
               << "   dissProbVectSize = " << dissProbVectSize
               << "   numLeftTot = " << numLeftTot
               << "   numChangeTot = " << numChangeTot << endl;
          cout << "       DISS_END phaseid count_ dim_isite numleft numchange  "
                  ":  "
               << dissPhaseIDVect[i] << "   " << "   "
               << count_[dissPhaseIDVect[i]] << "   " << dim_isite[i] << "   "
               << numLeft[i] << "   " << numChange[i] << endl;
          cout.flush();
        }
      }
    }
    if (callGEM > -1)
      break;

    dissProbVectSize = dissProbVect.size();
  }

  // findInterfaces_check();

  numadded_D = numChangeTot;
  return (callGEM);
}

void Lattice::findInterfaces_check(void) {
  unsigned int i, kk;
  unsigned int k;
  vector<Site *> gsite, dsite;

  struct int_phase {
    unsigned int ds;
    unsigned int gs;
  };

  int_phase locInt;
  vector<int_phase> interfaceSizes;
  interfaceSizes.clear();

  vector<Site *>::iterator beginLocation, endLocation;
  vector<Isite>::iterator start, end;

  ///
  /// An interface must have at least one adjacent site that is water or void
  ///

  cout << endl
       << "***********   from findInterfaces_check   ************" << endl;

  for (i = 0; i < numMicroPhases_; i++) {
    if (i != ELECTROLYTEID && i != VOIDID) { // Solid phase of some kind
      gsite.clear();
      dsite.clear();
      for (k = 0; k < numsites_; k++) {
        if (site_[k].getWmc() > 0) { // Is there some water nearby?
          if ((site_[k].getMicroPhaseId() == i)) {
            dsite.push_back(&site_[k]);
            site_[k].setDissolutionSite(i);
            for (kk = 0; kk < site_[k].nbSize(2); kk++) {
              if ((site_[k].nb(kk))->getMicroPhaseId() == ELECTROLYTEID) {
                gsite.push_back(site_[k].nb(kk));
                site_[k].nb(kk)->setGrowthSite(i);
              }
            }

            /// @note ThinterfaceSizesere is no reason to make the phase
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

      ///
      /// Eliminate duplicate values from the growth site vector
      ///

      if (gsite.size() > 0) {
        sort(gsite.begin(), gsite.end());
        beginLocation = gsite.begin();
        endLocation = unique(gsite.begin(), gsite.end());
        gsite.erase(endLocation, gsite.end());
      }

      ///
      /// Now do the same thing to the dissolution site vector
      ///

      if (dsite.size() > 0) {
        sort(dsite.begin(), dsite.end());
        beginLocation = dsite.begin();
        endLocation = unique(dsite.begin(), dsite.end());
        dsite.erase(endLocation, dsite.end());
      }

      cout << "findInterfaces_check  phaseId/dsite.size()/gsite.size()   " << i
           << "   " << dsite.size() << "   " << gsite.size() << endl;

      locInt.ds = dsite.size();
      locInt.gs = gsite.size();
      interfaceSizes.push_back(locInt);

    } else { // Not a solid phase (either electrolyte or void)

      locInt.ds = 0;
      locInt.gs = 0;
      interfaceSizes.push_back(locInt);
    }
  }
  // return interfaceSizes;

  bool differentInterfaceSize = false;
  for (int i = 2; i < numMicroPhases_; i++) {
    if (i != 6 && i != 13) {
      if (interfaceSizes[i].ds != interface_[i].getDissolutionNumSites()) {
        cout << "       dissolution interfaces are different for phaseId = "
             << i << endl;
        differentInterfaceSize = true;
      }
      if (interfaceSizes[i].gs != interface_[i].getGrowthNumSites()) {
        cout << "            growth interfaces are different for phaseId = "
             << i << endl;
        differentInterfaceSize = true;
      }
    }
  }
  if (differentInterfaceSize) {
    cout << endl << "********** dissolution interfaces **********" << endl;
    for (int i = 0; i < numMicroPhases_; i++) {
      cout << " phaseId" << endl;
      cout << "   " << i << endl;
      cout << "           dissolutionSites_.size()  = "
           << interface_[i].getDissolutionNumSites() << endl;
      cout << "           findInterfaces_check size = " << interfaceSizes[i].ds
           << endl;
    }
    cout << endl << "********** growth interfaces **********" << endl;
    for (int i = 0; i < numMicroPhases_; i++) {
      cout << " phaseId" << endl;
      cout << "   " << i << endl;
      cout << "           growthSites_.Size()         = "
           << interface_[i].getGrowthNumSites() << endl;
      cout << "           findInterfaces_check size   = "
           << interfaceSizes[i].gs << endl;
    }
    cout << endl << "*** stop program ***" << endl;
    exit(1);
  } else {
    cout << "              findInterfaces_check  OK!" << endl << endl;
  }
}

void Lattice::addDissolutionSite(Site *ste, unsigned int pid) {
  try {
    interface_.at(pid).addDissolutionSite(ste);
    // vector<unsigned int> plist = ste->getGrowthPhases();
    // for (unsigned int i = 0; i < plist.size(); i++) {
    //     interface_.at(plist[i]).removeGrowthSite(ste);
    // }
    ste->setDissolutionSite(pid);
  } catch (out_of_range &oor) {
    EOBException ex("Lattice", "addDissolutionSite", "interface_",
                    interface_.size(), pid);
    ex.printException();
    cout << endl << "id pid " << ste->getId() << " " << pid << endl;
    exit(1);
  }
  return;
}

void Lattice::addGrowthSite(Site *ste, unsigned int pid) {
  try {
    interface_.at(pid).addGrowthSite(ste);
    ste->setGrowthSite(pid);
  } catch (out_of_range &oor) {
    EOBException ex("Lattice", "addGrowthSite", "interface_", interface_.size(),
                    pid);
    ex.printException();
    cout << endl << "id pid " << ste->getId() << " " << pid << endl;
    exit(1);
  }
}

/*
// do not delete!
void Lattice::removeDissolutionSite(Site *ste, unsigned int pid) {
  try {
    interface_.at(pid).removeDissolutionSite(ste);
    ste->removeDissolutionSite(pid);
  } catch (out_of_range &oor) {
    EOBException ex("Lattice", "removeDissolutionSite", "interface_",
                    interface_.size(), pid);
    ex.printException();
    exit(1);
  }
}
*/

void Lattice::removeDissolutionSite_diss(Site *ste, unsigned int pid,
                                         int interfacePos) {
  try {
    interface_.at(pid).removeDissolutionSite_diss(ste, interfacePos);
    ste->removeDissolutionSite();
  } catch (out_of_range &oor) {
    EOBException ex("Lattice", "removeDissolutionSite_diss", "interface_",
                    interface_.size(), pid);
    ex.printException();
    cout << endl
         << "id pid interfacePos " << ste->getId() << " " << pid << " "
         << interfacePos << endl;
    exit(1);
  }
}

void Lattice::removeDissolutionSite_grow(Site *ste, unsigned int pid) {
  try {
    interface_.at(pid).removeDissolutionSite_grow(ste);
    ste->removeDissolutionSite();
  } catch (out_of_range &oor) {
    EOBException ex("Lattice", "removeDissolutionSite_grow", "interface_",
                    interface_.size(), pid);
    ex.printException();
    cout << endl << "id pid " << ste->getId() << " " << pid << endl;
    exit(1);
  }
}

/*
// do not delete!
void Lattice::removeGrowthSite(Site *ste, unsigned int pid) {
  try {
    interface_.at(pid).removeGrowthSite(ste);
    ste->removeGrowthSite(pid);
  } catch (out_of_range &oor) {
    EOBException ex("Lattice", "removeGrowthSite", "interface_",
                    interface_.size(), pid);
    ex.printException();
    exit(1);
  }
}
*/

void Lattice::removeGrowthSite_0(Site *ste, unsigned int pid,
                                 int interfacePos) {
  try {
    interface_.at(pid).removeGrowthSite_0(ste, interfacePos);
    ste->removeGrowthSite(pid);
  } catch (out_of_range &oor) {
    EOBException ex("Lattice", "removeGrowthSite_0", "interface_",
                    interface_.size(), pid);
    ex.printException();

    cout << endl
         << "id pid interfacePos " << ste->getId() << " " << pid << " "
         << interfacePos << endl;
    exit(1);
  }
}

void Lattice::removeGrowthSite_1(Site *ste, unsigned int pid) {
  try {
    interface_.at(pid).removeGrowthSite_1(ste);
    ste->removeGrowthSite(pid);
  } catch (out_of_range &oor) {
    EOBException ex("Lattice", "removeGrowthSite_1", "interface_",
                    interface_.size(), pid);
    ex.printException();

    cout << endl << "id pid : " << ste->getId() << " " << pid << endl;
    exit(1);
  }
}

int Lattice::emptyPorosity(int numsites, int cyc) {
  unsigned int i, j;
  int numemptied = 0;
  int maxsearchsize = 3;
  unsigned int cntpore, cntmax;
  bool placed;
  Site *stenb;

  ///
  /// Finding all potential VOID sites.
  ///
  /// @todo Consider removing some of the standard output, or setting a flag for
  /// it.
  ///

  cout << endl
       << "   Lattice::emptyPorosity - check for cyc = " << cyc
       << " :      numsites = " << numsites;
  cout.flush();

  list<Sitesize> distlist =
      findDomainSizeDistribution(ELECTROLYTEID, numsites, maxsearchsize, 0);
  list<Sitesize>::iterator it;

  // cout << endl << "Lattice::emptyPorosity - check for cyc = " << cyc
  //      << " :      numsites = " << numsites
  cout << "      distlist.size() = " << distlist.size() << endl;
  cout.flush();

  if (verbose_) {
    cout << "Lattice::emptyPorosity Found " << distlist.size()
         << " potential void sites." << endl;
    cout.flush();
  }

  ///
  /// We want to empty the sites with the largest pore count
  ///

  if (distlist.size() < numsites) {
    cout << endl
         << "   Lattice::emptyPorosity - Ran out of water in the system for "
            "cyc = "
         << cyc << "   distlist.size(): " << distlist.size()
         << "   numsites: " << numsites << endl;
    string msg = "Ran out of water in the system";
    EOBException ex("Lattice", "emptysite", msg, distlist.size(), numsites);
    ex.printException();
    exit(1);
  }

  numemptied = 0;
  it = distlist.begin();
  int siteID, siteIndex;
  vector<unsigned int> grVect;
  int dimVect;
  int itt = -1;
  double wmcIni, wmcEnd;
  try {

    while (it != distlist.end()) {
      siteIndex = (*it).siteid;
      siteID = site_[siteIndex].getId();
      setMicroPhaseId(siteID, VOIDID);
      grVect = site_[siteID].getGrowthPhases();
      dimVect = grVect.size();
      itt++;
      for (int i = 0; i < dimVect; i++) {
        if (!interface_[grVect[i]].removeEmptiedSite(siteID)) {
          cout << endl
               << "Lattice::emptyPorosity error: siteID = " << siteID
               << "does not belong to the interface phaseID = " << grVect[i]
               << "!!!" << endl;
          cout << "cyc = " << cyc << "   i = " << i
               << "   dimVect = " << dimVect << "   siteIndex = " << siteIndex
               << "   itt = " << itt << "   distlist.size = " << distlist.size()
               << endl;
          cout << "stop program" << endl;
          exit(1);
        }
      }
      site_[siteID].setWmc0(0);
      site_[siteID].dWmc(-1);
      for (int i = 0; i < NN_NNN; i++) {
        stenb = site_[siteID].nb(i);
        stenb->dWmc(-1);
        if ((stenb->getWmc() == 0) &&
            (stenb->getMicroPhaseId() > ELECTROLYTEID))
          removeDissolutionSite_grow(stenb, stenb->getMicroPhaseId());
      }
      site_[siteID].clearGrowthSite();
      numemptied++;
      it++;
    }

  } catch (EOBException eex) {
    eex.printException();
    cout << "Current value of numemptied = " << numemptied;
    cout.flush();
    exit(1);
  }

  return (numemptied);
}

int Lattice::fillPorosity(int numsites, int cyc) {
  unsigned int i, j;
  int numfilled = 0;
  int maxsearchsize = 10;
  unsigned int cntpore, cntmin;
  bool placed;
  Site *stenb;

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

  cout << endl
       << "   Lattice::fillPorosity - check for cyc = " << cyc
       << " :      numsites = " << numsites;
  cout.flush();

  list<Sitesize> distlist =
      findDomainSizeDistribution(VOIDID, numsites, maxsearchsize, 1);
  list<Sitesize>::iterator it;

  // cout << endl << "Lattice::fillPorosity - check for cyc = " << cyc
  //      << " :      numsites = " << numsites
  cout << "      distlist.size() = " << distlist.size() << endl;
  cout.flush();

  ///
  /// We want to fill the sites with the smallest pore count
  ///

  numfilled = 0;
  it = distlist.begin();
  int siteID, siteIndex;
  int nbpid, nbid;
  try {
    while (it != distlist.end()) {
      siteIndex = (*it).siteid;
      siteID = site_[siteIndex].getId();
      setMicroPhaseId(siteID, ELECTROLYTEID);
      site_[siteID].setWmc0(1);
      site_[siteID].dWmc(1);
      for (int i = 0; i < NN_NNN; i++) {
        stenb = site_[siteID].nb(i);
        // wmcIni = stenb->getWmc();
        stenb->dWmc(1);
        nbpid = stenb->getMicroPhaseId();
        if (nbpid > ELECTROLYTEID) {
          addDissolutionSite(stenb, nbpid);
          addGrowthSite(&site_[siteID], nbpid);
          if (i < NUM_NEAREST_NEIGHBORS) {
            for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases_;
                 phaseTmpl++) {
              if (chemSys_->isGrowthTemplate(phaseTmpl, nbpid)) {
                addGrowthSite(&site_[siteID], phaseTmpl);
              }
            }
          }
        }
      }
      numfilled++;
      it++;
    }
  } catch (EOBException eex) {
    eex.printException();
    cout << "Current value of numfilled = " << numfilled;
    cout.flush();
    exit(1);
  }

  return (numfilled);
}

int Lattice::countBox(int boxsize, unsigned int siteid) {
  string msg;
  int boxhalf = boxsize / 2;
  int nfound = 0;
  int ix, iy, iz, hx, hy, hz;

  try {
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
          if (site_.at(getIndex(hx, hy, hz)).getMicroPhaseId() ==
                  ELECTROLYTEID ||
              site_.at(getIndex(hx, hy, hz)).getMicroPhaseId() == VOIDID) {
            nfound++;
          }
        }
      }
    }
  } catch (out_of_range &oor) {
    EOBException ex("Lattice", "countBox", "site_", site_.size(),
                    getIndex(hx, hy, hz));
    ex.printException();
    exit(1);
  }
  return (nfound);
}

void Lattice::setResolution(const double res) {
  ///
  /// Make sure that resolution is a valid value
  ///

  try {
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
  } catch (DataException dex) {
    dex.printException();
    exit(1);
  }
  return;
}

vector<unsigned int> Lattice::getNeighborhood(const unsigned int sitenum,
                                              const int size) {
  int xp, yp, zp;
  // double dist;

  vector<unsigned int> nh;
  nh.clear();

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
  vector<double> vol_next, vfrac_next;
  vector<int> netsites;
  vector<unsigned int> pid;
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
  cout << endl << "Lattice::changeMicrostructure -ini- cyc = " << cyc << endl;
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
    cout << "Lattice::changeMicrostructure Before adjustMicrostructureVolumes:"
         << endl;
    for (int iii = 0; iii < volNextSize; ++iii) {
      cout << "Lattice::changeMicrostructure    Volume of " << phasenames[iii]
           << " = " << vol_next[iii] << " m3" << endl;
    }
    cout.flush();
  }

  try {
    adjustMicrostructureVolumes(vol_next, volNextSize, cyc);

    if (verbose_) {
      cout << "Lattice::changeMicrostructure After adjustMicrostructureVolumes:"
           << endl;
      for (int iii = 0; iii < volNextSize; ++iii) {
        cout << "Lattice::changeMicrostructure   Volume of " << phasenames[iii]
             << " = " << vol_next[iii] << " m3" << endl;
      }
    }

    vfrac_next = vol_next;

    adjustMicrostructureVolFracs(phasenames, vol_next, vfrac_next, volNextSize);
  } catch (DataException dex) {
    throw dex;
  } catch (EOBException ex) {
    throw ex;
  }

  ///
  /// Calculate number of sites of each phase in next state
  ///

  // Remember we added two extra slots onto the end of vfrac_next
  // to hold the volume fraction of capillary space and subvoxel pore space
  // But those two are not needed here anymore

  if (verbose_) {
    cout << "Lattice::changeMicrostructure Calculating volume "
         << "of each phase to be added..." << endl;
    try {
      for (i = 0; i < volNextSize; i++) {
        cout << "Lattice::changeMicrostructure ****Volume fraction["
             << phasenames.at(i)
             << "] in next state should be = " << vfrac_next.at(i);
        cout << ", or " << (int)((double)(numsites_ * vfrac_next.at(i)))
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
    } catch (out_of_range &oor) {
      throw EOBException("Lattice", "changeMicrostructure", "phasenames",
                         phasenames.size(), i);
    }
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

    netsites.clear();
    netsites.resize(numMicroPhases_, 0);
    pid.clear();
    pid.resize(numMicroPhases_, 0);

    try {
      cout << endl;
      for (i = FIRST_SOLID; i < volNextSize; i++) { // from i = 2 !!!
        cursites = (int)(count_.at(i) + 0.5);
        newsites = (int)((numsites_ * vfrac_next.at(i)) + 0.5);
        netsites.at(i) = newsites - cursites;
        pid.at(i) = i;
        phName[i] = phasenames.at(i);
        if (netsites.at(i) != 0) {
          cout << "Lattice::changeMicrostructure ***netsites["
               << phasenames.at(i) << "] in this state = " << netsites.at(i)
               << "; cursites = " << cursites << " and newsites = " << newsites
               << "   -> phaseID = "
               << chemSys_->getMicroPhaseId(phasenames.at(i)) << endl;
          cout.flush();
        }
      }
    } catch (out_of_range &oor) {
      throw EOBException(
          "Lattice", "changeMicrostructure",
          "phasenames or count_ or pid or netsites or vfrac_next",
          phasenames.size(), i);
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
  dissPhaseIDVect.clear();
  growPhaseIDVect.clear();
  numSiteDissVect.clear();
  numSiteGrowVect.clear();
  dissPhNameVect.clear();
  growPhNameVect.clear();
  int numTotSitesToDissolve = 0, numTotSitesToGrow = 0;

  try {
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
  } catch (out_of_range &oor) {
    throw EOBException("Lattice", "changeMicrostructure_0", "pid", pid.size(),
                       ii);
  }

  ///
  /// The next loop starts at FIRST_SOLID because we exclude void and water
  /// phases
  ///
  /// @todo Consider making the starting index more general
  ///

  cout << "Lattice::changeMicrostructure netsites.size() = " << netsitesSize_
       << " -> normalTRC/repeat (totalTRC/normalTRC/totalRepeat): " << normalTRC
       << " / " << repeat << "   (" << totalTRC << " / " << normalTRC << " / "
       << totalRepeat << ")" << endl;

  unsigned int pid_;
  int numadded_D = 0, numadded_G = 0;
  vector<int> nucleatedPhases;
  nucleatedPhases.resize(growPhaseIDVect.size(), 0);

  try {
    if (dissPhaseIDVect.size() > 0) {
      int needRecallGEM = -1;
      needRecallGEM = dissolvePhase(dissPhaseIDVect, numSiteDissVect,
                                    dissPhNameVect, numadded_D, totalTRC);
      if (needRecallGEM > -1) {
        phDiff = dissPhaseIDVect[needRecallGEM];
        nameDiff = dissPhNameVect[needRecallGEM];
        numDiff = count_[phDiff];

        cout << "Lattice::changeMicrostructure - anormal end" << endl;
        cout << "   phDiff,nameDiff,numDiff,count_[phDiff] : " << phDiff
             << " , " << nameDiff << " , " << numDiff << " , " << count_[phDiff]
             << endl;
        cout << "   numTotSitesToDissolve = " << numTotSitesToDissolve
             << "  while numTotSitesDissolved = " << numadded_D << endl;
        cout << "=> recall GEM after (re)setDCLowerLimit according to the "
                "system configuration (lattice)"
             << endl
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

  cursites = count_.at(VOIDID);
  newsites = (int)((numsites_ * vfrac_next.at(VOIDID)) + 0.5);
  wcursites = count_.at(ELECTROLYTEID);
  wnewsites = wcursites - (newsites - cursites);

  int numEmptyFill = newsites - cursites;
  int numEmpty = 0, numFill = 0;
  // numEmpty = emptyPorosity(newsites - cursites, cyc);
  // cout << endl << "Lattice::changeMicrostructure - before empty/fillPorosity
  // - check for cyc = " << cyc << endl;cout.flush();
  if (numEmptyFill > 0) {
    numEmpty = emptyPorosity(numEmptyFill, cyc);
  } else if (numEmptyFill < 0) {
    numFill = fillPorosity(-numEmptyFill, cyc);
  }
  // cout << "Lattice::changeMicrostructure - after empty/fillPorosity - check
  // for cyc = " << cyc << endl;cout.flush();

  if (verbose_) {
    cout << "Lattice::changeMicrostructure ***netsites["
         << phasenames.at(VOIDID)
         << "] in this state = " << (newsites - cursites)
         << "; cursites = " << cursites << " and newsites = " << newsites
         << endl;
    cout << "Lattice::changeMicrostructure ***netsites["
         << phasenames.at(ELECTROLYTEID)
         << "] in this state = " << (wnewsites - wcursites)
         << "; cursites = " << wcursites << " and newsites = " << wnewsites
         << endl
         << endl;

    // When creating void from water, we should
    // update the target volume fraction of water even though
    // it is not used in any further calculations at this point

    cout << "Lattice::changeMicrostructure Target CAPIILARY WATER "
         << "volume fraction IS " << vfrac_next.at(ELECTROLYTEID) << endl;

    // vfrac_next.at(ELECTROLYTEID) -= ((double)(newsites -
    // cursites)/(double)(numsites_)); cout << "But WILL BE " <<
    // vfrac_next.at(ELECTROLYTEID) << " after creating void space" << endl;

    cout << "Lattice::changeMicrostructure Number CAPILLARY VOXELS "
         << "actually emptied was:  " << numEmptyFill << endl;

    ///
    /// Report on target and actual mass fractions
    ///

    cout << "Lattice::changeMicrostructure "
         << "*******************************" << endl;
    cout.flush();
  }

  try {
    int totcount = 0;
    double numSites = numsites_;
    for (i = 0; i < volNextSize; i++) {
      volumefraction_.at(i) = count_.at(i) / numSites;
      totcount += count_.at(i);
      if (verbose_) {
        cout << "Lattice::changeMicrostructure Phase " << i
             << " Target volume fraction was " << vfrac_next[i]
             << " and actual is " << volumefraction_.at(i) << ", and "
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
      cout << "numDiff,phDiff,nameDiff : " << numDiff << " , " << phDiff
           << " , " << nameDiff << endl;
      cout << "repeat  : " << repeat << endl;
      cout << "cyc  : " << cyc << endl;
      cout << "stop program" << endl;
      exit(1);
    }
  } catch (out_of_range &oor) {
    throw EOBException("Lattice", "changeMicrostructure",
                       "volumefraction_ or count_", volumefraction_.size(), i);
  }

  ///  This is a local variable and the value is never used.
  ///
  ///  @todo Why not eliminate this line completely?
  ///

  // double surfa = getSurfaceArea(chemSys_->getMicroPhaseId(CSHMicroName));

  if (volumefraction_[ELECTROLYTEID] <= 0.0) {
    capWater = false;
    cout << endl
         << "Lattice::changeMicrostructure :  capWater = " << capWater << endl;
    cout << "mPhId/vfrac_next_/volumefraction_/count_: " << endl;
    for (int i = 0; i < volNextSize; i++) {
      cout << "   " << i << "  " << vfrac_next[i] << "  " << volumefraction_[i]
           << "  " << count_[i] << endl;
    }
    cout << endl
         << "time, simtype, capWater : " << time << " , " << simtype << " , "
         << capWater << endl;
    cout << "numDiff ,phDiff, nameDiff : " << numDiff << " , " << phDiff
         << " , " << nameDiff << endl;
    cout << "repeat  : " << repeat << endl;
    cout << "cyc  : " << cyc << endl;
    cout << "no water in the system => normal end of the program" << endl;
    exit(0);
  }

  cout << endl
       << "Lattice::changeMicrostructure -normal end- cyc = " << cyc << endl;

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

  try {

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
        solidvolumewithpores_ += vol.at(i);
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
        subvoxelporevolume_ += (vol.at(i) * chemSys_->getMicroPhasePorosity(i));
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
    watervolume_ = vol.at(ELECTROLYTEID);
    // Up to this point, vol.at(ELECTROLYTEID) is the total volume
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

    /// From here to the end of this time cycle, vol.at(ELECTROLYTEID) is the
    /// volume of capillary pore water only, not the total volume of water.
    /// We can always recover the total volume of water by asking ChemicalSystem
    /// for it.

    capillaryporevolume_ = capillaryvoidvolume_ + capillarywatervolume_;
    vol.at(ELECTROLYTEID) = capillarywatervolume_;
    vol.at(VOIDID) = capillaryvoidvolume_;

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

    cout << endl
         << "   Lattice::adjustMicrostructureVolumes - cyc = " << cyc << endl;
    cout << "     watervolume_/capPoreVol/capVoidVol/capWaterVol :   "
         << watervolume_ << " / " << capillaryporevolume_ << " / "
         << capillaryvoidvolume_ << " / " << capillarywatervolume_ << endl;

    ///
    /// End of manual adjustment
    ///

  }

  catch (DataException dex) {
    throw dex;
  }

  catch (out_of_range &oor) {
    throw EOBException("Lattice", "adjustMicrostructureVolumes", "volSize",
                       volSize, i);
  }

  return;
}

void Lattice::adjustMicrostructureVolFracs(vector<string> &names,
                                           const vector<double> vol,
                                           vector<double> &vfrac, int volSize) {
  int i = 0;
  double totmicvolume = 0.0;
  // int volSize = vol.size();

#ifdef DEBUG
  cout << "Lattice::adjustMicrostructureVolFracs" << endl;
  cout.flush();
#endif

  try {

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
    //     totmicvolume += vol.at(i);
    //     #ifdef DEBUG
    //         cout << "Lattice::adjustMicrostructureVolFracs Volume("
    //              << names.at(i) << ") = " << vol.at(i)
    //              << ", volfrac = " << vfrac.at(i) << endl;
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
      vfrac[i] = vol[i] / initialmicrostructurevolume_; // totmicvolume;
    }

    if (verbose_) {
      for (i = 0; i < volSize; ++i) {
        cout << "Lattice::adjustMicrostructureVolFracsVolume "
             << "fraction[" << names.at(i) << "] should be " << vfrac.at(i)
             << ", (" << vol.at(i) << "/"
             << initialmicrostructurevolume_ // totmicvolume
             << ") and volume fraction NOW is "
             << (double)(count_.at(i)) / (double)(numsites_) << endl;
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

  catch (DataException dex) {
    dex.printException();
    exit(1);
  }

  catch (out_of_range &oor) {
    EOBException ex("Lattice", "adjustMicrostructureVolFracs", "volSize",
                    volSize, i);
    ex.printException();
    exit(1);
  }
}

void Lattice::calcSubvoxelporevolume(vector<double> &vol) {
  int i = 0;
  try {

    // Find the total system volume according to GEMS, in m3
    // units.  The individual microstructure phase volumes
    // have already // been adjusted for subvoxel porosity in the
    // ChemicalSystem::calculateState function

    // The lattice has a fixed volume given by a fixed number
    // of voxels multiplied by the volume per voxel

    // This will hold the subvoxel pore volume (m3)

    subvoxelporevolume_ = 0.0;
    // double phi; // Holds the subvoxel porosity of a microstructurephase
    for (i = 0; i < vol.size(); ++i) {
      if (i != ELECTROLYTEID && i != VOIDID) {
        subvoxelporevolume_ += (vol.at(i) * chemSys_->getMicroPhasePorosity(i));
      }
    }

    // The total amount of non-solid space in the microstructure
  }

  catch (out_of_range &oor) {
    throw EOBException("Lattice", "calculateSubvoxelporevolume", "vol",
                       vol.size(), i);
  }
}

void Lattice::calcSolidvolumewithpores(vector<double> &vol) {
  int i = 0;
  try {

    // Find the total system volume according to GEMS, in m3
    // units.  The individual microstructure phase volumes
    // have already // been adjusted for subvoxel porosity in the
    // ChemicalSystem::calculateState function

    // The lattice has a fixed volume given by a fixed number
    // of voxels multiplied by the volume per voxel

    // This will hold the subvoxel pore volume (m3)

    solidvolumewithpores_ = 0.0;
    for (i = 0; i < vol.size(); ++i) {
      if (i != ELECTROLYTEID && i != VOIDID) {
        solidvolumewithpores_ += vol.at(i);
      }
    }
  }

  catch (out_of_range &oor) {
    throw EOBException("Lattice", "calcSolidVolumeWithPores", "vol", vol.size(),
                       i);
  }
}

void Lattice::calcCapillarywatervolume(vector<double> &vol) {
  calcSubvoxelporevolume(vol);
  capillarywatervolume_ = vol.at(ELECTROLYTEID) - subvoxelporevolume_;
  if (capillarywatervolume_ < 0.0)
    capillarywatervolume_ = 0.0;
  if (chemSys_->isSaturated()) {
    capillarywatervolume_ = vol.at(ELECTROLYTEID) - subvoxelporevolume_;
  }
}

void Lattice::calcCapillaryvoidvolume(vector<double> &vol) {
  capillaryporevolume_ = nonsolidvolume_ - subvoxelporevolume_;
  capillaryvoidvolume_ = capillaryporevolume_ - capillarywatervolume_;
}

void Lattice::calculatePoreSizeDistribution(void) {
  // First compose the full pore volume distribution

  vector<double> subpore_volume;
  subpore_volume.resize(volumefraction_.size(), 0.0);

  // Following will hold the subvoxel porosity of a phase
  double phi = 0.0;
  double upper_cutoff_porosity = 0.8;

  // Get the combined pore volume as a function of diameter
  for (int i = 0; i < volumefraction_.size(); ++i) {
    if (i != ELECTROLYTEID && i != VOIDID) {
      phi = chemSys_->getMicroPhasePorosity(i);
      if (phi > upper_cutoff_porosity)
        phi = upper_cutoff_porosity;
      subpore_volume[i] = volumefraction_.at(i) * phi;
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
  vector<vector<struct PoreSizeVolume>> binnedporevolume;
  vector<struct PoreSizeVolume> zpsvec;
  zpsvec.clear();
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
        while ((porevolume[i][j].diam <= maxsize) &&
               (j < porevolume[i].size())) {
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
  // microstructure volume basis = volumefraction_.at(ELECTROLYTEID)
  //
  // Volume fraction of unsaturated capillary pores
  // on a total microstructure volume basis
  // = volumefraction_.at(VOIDID)
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

  double water_volfrac = water_volume / microstructurevolume_;

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
  //       << "void fraction = " << volumefraction_.at(VOIDID) << endl;
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

  if (volumefraction_.at(ELECTROLYTEID) >= 0.0) {
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
  // microstructure volume basis = volumefraction_.at(ELECTROLYTEID)
  //
  // Volume fraction of unsaturated capillary pores
  // on a total microstructure volume basis
  // = volumefraction_.at(VOIDID)
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

  double water_volfrac = watervolume_ / microstructurevolume_;

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
    cout << "Capillary void volume fraction = " << volumefraction_.at(VOIDID)
         << endl;
    cout << "Saturated capillary pore volume fraction = "
         << capillaryporevolumefraction_ - volumefraction_.at(VOIDID) << endl;
    cout << "Nanopore volume fraction (<= 100 nm) = "
         << subvoxelporevolumefraction_ << endl;
    cout << "Total pore volume fraction = " << pore_volfrac << endl << endl;
    cout << "Total void volume fraction = " << volumefraction_.at(VOIDID)
         << endl;
    cout << "Pore size saturation data:" << endl;
    cout << "Diameter (nm),Volume Fraction,Fraction Saturated" << endl;
    cout << "Masterporevolume size = " << masterporevolume_.size() << endl;
    cout.flush();
  }

  out << "Time = " << (curtime * 24.0) << " h" << endl;
  out << "Capillary pore volume fraction (> 100 nm) = "
      << capillaryporevolumefraction_ << endl;
  out << "Capillary void volume fraction = " << volumefraction_.at(VOIDID)
      << endl;
  out << "Saturated capillary pore volume fraction = "
      << capillaryporevolumefraction_ - volumefraction_.at(VOIDID) << endl;
  out << "Nanopore volume fraction (<= 100 nm) = "
      << subvoxelporevolumefraction_ << endl;
  out << "Total pore volume fraction = " << pore_volfrac << endl;
  out << "Total void volume fraction = " << volumefraction_.at(VOIDID) << endl;
  out << "Pore size saturation data:" << endl;
  out << "Diameter (nm),Volume Fraction,Fraction Saturated" << endl;
  for (int i = 0; i < masterporevolume_.size(); ++i) {
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
  double capvoid_volfrac = volumefraction_.at(VOIDID) +
                           (volumefraction_.at(ELECTROLYTEID) - water_volfrac);
  double capspace_volfrac = capvoid_volfrac + capwater_volfrac;
  out << ">" << masterporevolume_[masterporevolume_.size() - 1].diam << ","
      << capspace_volfrac << "," << (1.0 - volumefraction_.at(VOIDID)) << endl;
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
  double red, green, blue;
  vector<double> colors;
  out << numMicroPhases_ << endl;
  for (int i = 0; i < numMicroPhases_; i++) {
    microPhaseId = chemSys_->getMicroPhaseId(i);
    colors = chemSys_->getColor(microPhaseId);
    out << microPhaseId << " " << colors[0] << " " << colors[1] << " "
        << colors[2] << endl;
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

  out << VERSIONSTRING << " " << version_ << endl;
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

    out1 << VERSIONSTRING << " " << version_ << endl;
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
  string ofileName;
  ostringstream ostr1, ostr2;
  ostr1 << setfill('0') << setw(6)
        << (int)((curtime * 24.0 * 60.0) + 0.5); // minutes
  ostr2 << setprecision(3) << temperature_;
  string timestr(ostr1.str());
  string tempstr(ostr2.str());
  ofileName = "imgIni." + timestr + "m." + tempstr + "_TEST.img";
  if (verbose_) {
    cout << "    In Lattice::writeLatticeIni, curtime = " << curtime
         << ", timestr = " << timestr << endl;
    cout.flush();
  }

  ofstream out(ofileName.c_str());

  // Write image header information first

  out << VERSIONSTRING << " " << version_ << endl;
  out << XSIZESTRING << " " << xdim_ << endl;
  out << YSIZESTRING << " " << ydim_ << endl;
  out << ZSIZESTRING << " " << zdim_ << endl;
  out << IMGRESSTRING << " " << resolution_ << endl;

  vector<int> siteTest;
  siteTest.clear();
  siteTest.resize(numsites_, 18);

  for (j = 0; j < numsites_; j++) {
    if (site_[j].getMicroPhaseId() == 5) {
      if (site_[j].getWmc() == 0) {
        siteTest[j] = 5;
        for (i = 0; i < NN_NNN; i++) {
          siteTest[site_[j].nb(i)->getId()] = site_[j].nb(i)->getMicroPhaseId();
        }
      }
    }
  }
  for (j = 0; j < numsites_; j++) {
    out << siteTest[j] << endl;
  }
  out.close();

  string ofileName1;
  ofileName1 = "imgIni." + timestr + "m." + tempstr + ".img";
  ofstream out1(ofileName1.c_str());
  out1 << VERSIONSTRING << " " << version_ << endl;
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
  out << "Properties=pos:R:3:color:R:3:transparency:R:1 ";
  out << "Time=" << timestr << endl;

  // Loop over all voxels and write out the solid ones

  float x, y, z;
  int mPhId;
  vector<double> colors;

  for (int i = 0; i < numsites_; i++) {
    x = site_[i].getX();
    y = site_[i].getY();
    z = site_[i].getZ();
    mPhId = site_[i].getMicroPhaseId();
    colors = chemSys_->getColor(mPhId);
    out << x << "\t" << y << "\t" << z << "\t" << colors[0] << "\t" << colors[1]
        << "\t" << colors[2] << "\t0.0" << endl;
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
  vector<double> colors;
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
    // for (mPhId = 2; mPhId < numMicroPhases_; mPhId++) {
    mPhNotExist = true;
    colors = chemSys_->getColor(mPhId);
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

  out << VERSIONSTRING << " " << version_ << endl;
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

  double red, green, blue;
  vector<double> colors;
  for (j = 0; j < ydim_; j++) {
    for (i = 0; i < xdim_; i++) {
      colors = chemSys_->getColor(image[i][j]);
      red = dshade[i][j] * colors[0] + 0.5;
      green = dshade[i][j] * colors[1] + 0.5;
      blue = dshade[i][j] * colors[2] + 0.5;
      out << (int)(red);
      out << " " << (int)(green);
      out << " " << (int)(blue) << endl;
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

  double red, green, blue;
  vector<double> colors;
  for (j = 0; j < ydim_; j++) {
    for (i = 0; i < xdim_; i++) {
      colors = chemSys_->getColor(image[i][j]);
      red = dshade[i][j] * colors[0] + 0.5;
      green = dshade[i][j] * colors[1] + 0.5;
      blue = dshade[i][j] * colors[2] + 0.5;
      out << (int)(red);
      out << " " << (int)(green);
      out << " " << (int)(blue) << endl;
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

    double red, green, blue;
    vector<double> colors;
    for (j = 0; j < ydim_; j++) {
      for (i = 0; i < xdim_; i++) {
        colors = chemSys_->getColor(image[i][j]);
        red = dshade[i][j] * colors[0] + 0.5;
        green = dshade[i][j] * colors[1] + 0.5;
        blue = dshade[i][j] * colors[2] + 0.5;
        out << (int)(red);
        out << " " << (int)(green);
        out << " " << (int)(blue) << endl;
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

list<Sitesize> Lattice::findDomainSizeDistribution(int phaseid,
                                                   const int numsites,
                                                   int maxsize,
                                                   int sortorder = 0) {
  // if sortorder is 0 => sorting in descending order
  list<Sitesize> distlist;
  list<Sitesize>::iterator it;
  bool found = false;

  Sitesize ss;
  int domainsize = 0;
  // for (int i = 0; i < site_.size(); i++) {
  for (int i = 0; i < numsites_; i++) {
    if (site_[i].getMicroPhaseId() == phaseid) {
      domainsize = findDomainSize(i, maxsize);
      ss.siteid = i;
      ss.nsize = domainsize;
      it = distlist.begin();
      if (distlist.empty()) {
        distlist.push_back(ss);
      } else {
        found = false;
        if (sortorder == 0) {
          it = distlist.begin();
          while (it != distlist.end() && !found) {
            if (ss.nsize >= (*it).nsize) {
              distlist.insert(it, ss);
              found = true;
            } else {
              it++;
            }
          }
        } else {
          it = distlist.begin();
          while (it != distlist.end() && !found) {
            if (ss.nsize <= (*it).nsize) {
              distlist.insert(it, ss);
              found = true;
            } else {
              it++;
            }
          }
        }
        if (!found && (distlist.size() < numsites)) {
          distlist.push_back(ss);
        }
      }
      if (distlist.size() > numsites)
        distlist.pop_back();
    }
  }

  return (distlist);
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
