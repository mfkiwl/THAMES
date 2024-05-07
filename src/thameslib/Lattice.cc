/**
@file Lattice.cc
@brief Definition of methods for the Lattice class.

*/
#include "Lattice.h"
#include "Interface.h"
#include "RanGen.h"
#include "imageutil.h"

Lattice::Lattice(ChemicalSystem *cs, Solution *solut)
    : siteneighbors_(NN_NNN), chemSys_(cs), solut_(solut) {
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

Lattice::Lattice(ChemicalSystem *cs, Solution *solut, const string &fileName,
                 const bool verbose, const bool warning)
    : siteneighbors_(NN_NNN), chemSys_(cs), solut_(solut) {
  unsigned int i, j, k;
  unsigned int ii;
  string buff;
  int xn, yn, zn;
  unsigned int idn;
  unsigned int pid;
  string msg;

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
    double firsttime = true;
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
    rg_ = new RanGen();
  } catch (bad_alloc &ba) {
    cout << "Lattice constructor failed when allocating rg_";
    cout.flush();
    exit(1);
  }

  rg_->setSeed(-142234);
  if (verbose_) {
    cout << "Lattice::Lattice Checking whether I can get a seed "
         << rg_->getSeed() << endl;
    cout << "Lattice::Lattice Checking value of random number "
         << ", " << rg_->Ran3() << endl;
    cout.flush();
  }

  count_.clear();
  count_.resize(chemSys_->getNumMicroPhases(), 0);

  SI_.clear();
  SI_.resize(chemSys_->getNumMicroPhases(), 1.0);

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
  /// With the phase counts known, calculate phase volume fractions
  /// We cannot calculate the mass fractions until we know the phase identities
  /// Phase identities are determined by the KineticModel
  ///

  int numMicroPhases = chemSys_->getNumMicroPhases();

  vector<double> microPhaseMass(numMicroPhases, 0.0);

  volumefraction_.clear();
  volumefraction_.resize(numMicroPhases, 0.0);

  initvolumefraction_.clear();
  initvolumefraction_.resize(numMicroPhases, 0.0);

  double vfrac, mfrac, capfrac, molarMass, molarVolume, density;
  double solidMass = 0.0;
  int microPhaseId = 0;
  int DCId = 0;
  string myname;
  try {
    if (verbose_) {
      cout << "Lattice::Lattice Calculating volume fractions now ..." << endl;
      for (ii = 0; ii < numMicroPhases; ii++) {
        cout << "Micro phase "
             << chemSys_->getMicroPhaseName(chemSys_->getMicroPhaseId(ii))
             << ", count = " << count_[ii] << " of " << site_.size() << endl;
      }
      cout.flush();
    }

    if (site_.size() > 0) {
      for (ii = 0; ii < numMicroPhases; ii++) {
        myname = chemSys_->getMicroPhaseName(chemSys_->getMicroPhaseId(ii));
        microPhaseId = chemSys_->getMicroPhaseId(ii);
        vfrac = ((double)count_[ii]) / ((double)site_.size());
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
          }
          if (verbose_ && vfrac > 0.0) {
            cout << "Lattice::Lattice Phase "
                 << chemSys_->getMicroPhaseName(microPhaseId) << endl;
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

    // Set the water-solids mass ratio based on the initial microstructure

    wsratio_ = microPhaseMass[ELECTROLYTEID] / solidMass;

    if (verbose_) {
      cout << "Lattice::Lattice Microstructure w/s = " << wsratio_ << endl;
      cout << "Lattice::Lattice (Mass of water = "
           << microPhaseMass[ELECTROLYTEID]
           << ", mass of solids = " << solidMass << ")" << endl;
      cout.flush();
    }

    // Next we set the initial normalized phase masses, microstructure
    // phase volume, and subvoxel porosity.  This is all triggered when we set
    // the mass.

    normalizePhaseMasses(microPhaseMass, solidMass);

    // Set the initial total volume of the microstructure

    double totmicvol = 0.0;
    for (int i = 0; i < numMicroPhases; i++) {
      microPhaseId = chemSys_->getMicroPhaseId(i);
      if (microPhaseId != VOIDID) {
        totmicvol += chemSys_->getMicroPhaseVolume(microPhaseId);
      }
    }

    chemSys_->setInitMicroVolume(totmicvol);

    // Initially assume that all free water and void space
    // is capillary volume

    capillaryporevolumefraction_ =
        getVolumefraction(ELECTROLYTEID) + getVolumefraction(VOIDID);

  } catch (FloatException fex) {
    fex.printException();
    exit(1);
  }

  int stId, phId;
  string nameMicroPhaseTh;
  double rng;
  for (ii = 0; ii < site_.size(); ii++) {
    stId = site_[ii].getId();
    phId = site_[ii].getMicroPhaseId();
    nameMicroPhaseTh = chemSys_->getMicroPhaseName(phId);
    if (nameMicroPhaseTh == "Electrolyte") {
      site_[ii].setWmc0(1.);
    } else if (nameMicroPhaseTh == "CSHQ") {
      rng = RanGen::Ran3();
      if (rng >= thrPorosityCSH) {
        site_[ii].setWmc0(chemSys_->getMicroPhasePorosity(phId));
      } else {
        site_[ii].setWmc0(0.);
      }
    } else {
      site_[ii].setWmc0(0.);
    }
  }

  // for (ii = 0; ii < site_.size(); ii++) {
  //     site_[ii].calcWmc();
  // }

  double wmcLoc;
  int nbIdLoc;
  for (ii = 0; ii < site_.size(); ii++) {
    wmcLoc = site_[ii].getWmc0();
    for (int jj = 0; jj < NN_NNN; jj++) {
      nbIdLoc = (site_[ii].nb(jj))->getId();
      wmcLoc += site_[nbIdLoc].getWmc0();
    }
    site_[ii].setWmc(wmcLoc);
  }
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
                                   double solidMass) {
  int microPhaseId;
  double pscaledMass = 0.0;

  int numMicroPhases = chemSys_->getNumMicroPhases();

  for (int i = 0; i < numMicroPhases; i++) {
    microPhaseId = chemSys_->getMicroPhaseId(i);
    if (microPhaseId == ELECTROLYTEID) {
      int waterId = chemSys_->getDCId("H2O@");
      double waterMolarMass = chemSys_->getDCMolarMass(waterId);
      pscaledMass = wsratio_ * 100.0; // Mass of solids scaled to 100 g now
                                      //
      chemSys_->setDCMoles(waterId, (pscaledMass / waterMolarMass));
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
      if (verbose_) {
        cout << "Lattice::normalizePhaseMasses Microstructure scaled mass of "
             << chemSys_->getMicroPhaseName(microPhaseId) << " ("
             << microPhaseId << ") = " << microPhaseMass[microPhaseId]
             << " g out of " << solidMass << " g total" << endl;

        // Setting the phase mass will also automatically calculate the phase
        // volume

        cout << "Lattice::normalizePhaseMasses Setting initial micphase "
             << "mass, volume, and DC moles of "
             << chemSys_->getMicroPhaseName(microPhaseId) << endl;
        cout.flush();
      }

      chemSys_->setMicroPhaseMass(microPhaseId, pscaledMass);

      chemSys_->setMicroPhaseMassDissolved(microPhaseId, 0.0);
    }
  }

  // Up to this point we could not really handle volume of void space, but
  // now we can do so in proportion to electrolyte volume

  double vfv = getVolumefraction(VOIDID);
  double vfe = getVolumefraction(ELECTROLYTEID);
  double ve = chemSys_->getMicroPhaseVolume(ELECTROLYTEID);

  chemSys_->setMicroPhaseVolume(VOIDID, (ve * vfv / vfe));

  return;
}

void Lattice::findInterfaces(void) {
  unsigned int i, kk;
  unsigned int k;
  vector<Site *> gsite, dsite;

  ///
  /// An interface must have at least one adjacent site that is water or void
  ///

  interface_.clear();
  for (i = 0; i < chemSys_->getNumMicroPhases(); i++) {
    if (i != ELECTROLYTEID && i != VOIDID) { // Solid phase of some kind
      gsite.clear();
      dsite.clear();
      for (k = 0; k < site_.size(); k++) {
        if (site_[k].getWmc() > 0) { // Is there some water nearby?
          if ((site_[k].getMicroPhaseId() == i)) {
            dsite.push_back(&site_[k]);
            site_[k].setDissolutionSiteMod(i);
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

      if ((gsite.size() == 0) && (dsite.size() == 0)) {

        cout << "microPhase " << i
             << " had both gsite.size() & dsite.size() equal to 0!" << endl;

        ///
        /// We are dealing with a phase that may need to
        /// nucleate.  Identify the eligible nucleation
        ///  sites for that phase
        ///

        double thresh = (0.5 / pow(resolution_, 3.0));
        double g = 0.0;
        for (k = 0; k < site_.size(); k++) {
          if ((site_[k].getMicroPhaseId() == ELECTROLYTEID)) {
            g = rg_->Ran3();
            if (g < thresh) {
              gsite.push_back(&site_[k]);
              site_[k].setGrowthSite(i);
            }
          }
        }
      }

      interface_.push_back(Interface(chemSys_, rg_, gsite, dsite, i, verbose_));

    } else { // Not a solid phase (either electrolyte or void)

      interface_.push_back(Interface(rg_, verbose_));
    }
  }

  cout << endl
       << "     Lattice::findInterfaces - count_ & interface sizes at the start"
       << endl;
  for (int i; i < chemSys_->getNumMicroPhases(); i++) {
    cout << "  " << i << "   " << chemSys_->getMicroPhaseName(i)
         << "   id: " << chemSys_->getMicroPhaseId(i)
         << "     count_ = " << count_[i] << "     dissInterfaceSize =  "
         << interface_[i].getDissolutionNumSites()
         << "     growInterfaceSize =  " << interface_[i].getGrowthNumSites()
         << "         porosity : " << chemSys_->getMicroPhasePorosity(i)
         << "     templates : ";
    for (int j = 0; j < chemSys_->getNumMicroPhases(); j++) {
      if (chemSys_->isGrowthTemplate(i, j)) {
        cout << j << " ";
        //}else{
        //    cout <<"";
      }
    }
    cout << endl;
  }
  cout.flush();

  // cout << "STOP90" << endl;
  // exit(1);

  /*
  int sizeInt;
  int site_Id, site_PhaseId;
  bool testBool = true, testBool_all = true;
  cout << endl << "     Lattice::findInterfaces - check phaseIds for all
  interfaces" <<endl; cout.flush(); cout << "DISSOLUTION:" << endl; for (int i =
  0; i < chemSys_->getNumMicroPhases(); i++){ sizeInt =
  interface_[i].getDissolutionNumSites(); cout << "   start check for the
  dissolution interface of phaseId = " << i << " having " << sizeInt << " sites"
  << endl; cout.flush(); for (int j = 0; j < sizeInt; j++){ site_Id =
  (interface_[i].getDissolutionSites())[j].getId(); site_PhaseId
  =site_[site_Id].getMicroPhaseId(); if(site_PhaseId != i){ cout << "      " <<
  i << "-dissInterface has on " << j << " position the site " << site_Id << "
  having phaseId = " << site_PhaseId << endl; cout.flush(); testBool = false;
          }
      }
      if (testBool){
          cout << "         OK for phaseId = " << i << endl; cout.flush();
      } else {
          testBool_all = false;
          testBool = true;
      }
  }
  cout << "GROWTH:" << endl;
  testBool = true;
  for (int i = 0; i < chemSys_->getNumMicroPhases(); i++){
      sizeInt = interface_[i].getGrowthNumSites();
      cout << "   start check for the growth interface of phaseId = " << i << "
  having " << sizeInt << " sites" << endl; cout.flush(); for (int j = 0; j <
  sizeInt; j++){ site_Id = (interface_[i].getGrowthSites())[j].getId();
          site_PhaseId =site_[site_Id].getMicroPhaseId();
          if(site_PhaseId != ELECTROLYTEID){
              cout << "      " << i << "-growInterface has on " << j << "
  position the site " << site_Id << " having phaseId = " << site_PhaseId <<
  endl; cout.flush(); testBool = false;
          }
      }
      if (testBool){
          cout << "         OK for phaseId = " << i << endl; cout.flush();
      } else {
          testBool_all = false;
          testBool = true;
      }
  }
  if (testBool_all == false){
      cout << endl << "     ***** there are problems in interface(s) *****" <<
  endl; cout << "stop program" << endl; exit(1); } else { cout << endl << "all
  interfaces are well defined" << endl;
  }
  exit(1);
*/
  return;
}

int Lattice::growPhaseMod(unsigned int phaseid, int numtoadd, int trc) {

  unsigned int j;

  int numMicroPhases = chemSys_->getNumMicroPhases();
  Site *ste, *stenb;
  vector<Isite> isite;
  isite = interface_[phaseid].getGrowthSites();
  int dim_isite = isite.size();
  int numleft = numtoadd, numchange = 0;

  int aff, affMin, affSum;
  double affSumDbl;
  int valAbs;
  double rng;
  int isitePos;
  unsigned int pid;

  double steWmc, stenbWmc, dwmcval;

  //*** for controll
  int static trc_d;
  int bcl;

  try {
    if (numtoadd == 0)
      return 0;
    if (phaseid >= interface_.size()) {
      throw EOBException("Lattice", "growPhaseMod", "interface_",
                         interface_.size(), phaseid);
    }
  } catch (EOBException ex) {
    ex.printException();
    exit(1);
  }

  ///
  /// We need to go through the interface list in
  /// normal order, which is the reverse order that
  /// we use for dissolution.
  ///

  if (verbose_) {
    cout << "Lattice::growPhaseMod -->Phase " << phaseid << " needs to grow at "
         << numtoadd << " sites" << endl;
    cout.flush();
  }

  /*
// TEST DISSOLVE/GROW INTERFACES start
cout<<endl<<endl<<"********************************************"<<endl;
// cout<<endl<<"diss phaseid trc_d/bcl dim_isite numleft numchange  :  " <<
phaseid << "   " << trc_d << "   " << bcl <<
//     "   " << dim_isite << "   " << numleft << "   " << numchange << endl;
// cout.flush();
cout << endl << "     count_ before growth" <<endl;
f  for (i = 0 ; i < numMicroPhases; i++){
    cout << "  " << i << "     count_ = " << count_[i]
         << "     dissInterfaceSize =  " <<
interface_[i].getDissolutionNumSites()
         << "     growInterfaceSize =  " << interface_[i].getGrowthNumSites() <<
"       & is template for : "; for (jj = 0; jj < numMicroPhases; jj++){ if
(chemSys_->isGrowthTemplate(jj,i)){ cout << jj << " ";
        }
    }
    cout << endl;
}
    cout << endl;
}
cout.flush();
// TEST DISSOLVE/GROW INTERFACES end
*/

  trc_d++;
  bcl = 0;

  cout << endl
       << "grow phaseid trc_d/bcl dim_isite numleft numchange  :  " << phaseid
       << "   " << trc_d << "   " << bcl << "   " << dim_isite << "   "
       << numleft << "   " << numchange << endl;
  cout.flush();

  while ((numleft > 0) && (dim_isite >= 1)) {
    try {
      bcl++;
      affSum = 0;
      affMin = 1000000;
      for (j = 0; j < dim_isite; j++) {
        aff = isite[j].getAffinity();
        affSum += aff;
        if (aff < affMin)
          affMin = aff;
      }

      // calc probabilities
      if (affSum != 0) {
        if (affMin < 0) {
          valAbs = abs(affMin);
          affSum += dim_isite * valAbs;
          affSumDbl = affSum;
          if (affSum != 0) {
            isite[0].setProb((isite[0].getAffinity() + valAbs) / affSumDbl);
            for (j = 1; j < dim_isite; j++) {
              isite[j].setProb(isite[j - 1].getProb() +
                               (isite[j].getAffinity() + valAbs) / affSumDbl);
            }
            rng = RanGen::Ran3();
            for (isitePos = 0; isitePos < dim_isite; isitePos++) {
              if (rng <= isite[isitePos].getProb())
                break;
            }
          } else {
            rng = RanGen::Ran3();
            isitePos = (int)(rng * dim_isite);
          }
        } else {
          affSumDbl = affSum;
          isite[0].setProb(isite[0].getAffinity() / affSumDbl);
          for (j = 1; j < dim_isite; j++) {
            isite[j].setProb(isite[j - 1].getProb() +
                             isite[j].getAffinity() / affSumDbl);
          }
          rng = RanGen::Ran3();
          for (isitePos = 0; isitePos < dim_isite; isitePos++) {
            if (rng <= isite[isitePos].getProb())
              break;
          }
        }
      } else {
        rng = RanGen::Ran3();
        isitePos = (int)(rng * dim_isite);
      }

      ste = &site_[isite[isitePos].getId()];
      pid = ste->getMicroPhaseId(); // always ELECTROLYTEID !!

      double wmcIni, wmcEnd;
      wmcIni = ste->getWmc0();

      // if (pid == ELECTROLYTEID) {
      // removeGrowthSite(ste, phaseid);
      removeGrowthSiteMod0_grow(ste, phaseid, isitePos);
      setMicroPhaseIdMod_grow(ste, phaseid);

      vector<unsigned int> plist = ste->getGrowthPhases();
      for (j = 0; j < plist.size(); j++) {
        removeGrowthSiteMod1_grow(ste, plist[j]);
      }

      ///
      /// Weighted mean curvature (wmc) is changed by the difference
      /// between the growing phase's porosity and the template's porosity.
      ///
      /// @todo Determine why the calculation works this way.
      ///

      string nameMicroPhaseTh = chemSys_->getMicroPhaseName(phaseid);
      if (nameMicroPhaseTh == "CSHQ") {
        rng = RanGen::Ran3();
        if (rng >= thrPorosityCSH) {
          wmcEnd = chemSys_->getMicroPhasePorosity(phaseid);
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
      //     cout << endl <<"error removeDissolutionSite for steWmc - ste pid
      //     steWmc phaseid trc trc_d bcl : "
      //          << ste->getId() << "   "
      //          << pid << "   " << steWmc << "   " << phaseid << "   "
      //          << trc << "   " << trc_d << "   " << bcl << endl;
      //          cout.flush();
      //     cout << endl << "nn & nnn siteId phaseId siteWmc :" << endl;
      //     for (int i = 0; i < chemSys_->getNumMicroPhases(); i++){
      //         cout << "   " << i << "   " << (ste->nb(i))->getId() << "   "
      //         << (ste->nb(i))->getMicroPhaseId() << "   " <<
      //         (ste->nb(i))->getWmc() << endl;
      //     }
      //     cout << endl << "stop program" << endl;
      //     exit(1);
      // }
      if (steWmc > 0.0) {
        addDissolutionSiteMod(ste, phaseid);
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
        //     cout << endl <<"error removeDissolutionSite for stenbWmc - ste
        //     pid steWmc phaseid trc trc_d bcl stenb stenb_pid stenbWmc j : "
        //          << ste->getId() << "   "
        //          << pid << "   " << ste->getWmc() << "   " << phaseid << " "
        //          << trc << "   " << trc_d << "   " << bcl << "   " <<
        //          stenb->getId() << "   "
        //          << stenb->getMicroPhaseId() << "   " << stenbWmc << "   " <<
        //          j << endl; cout.flush();
        //     cout << endl << "stop program" << endl;
        //     exit(1);
        // }
        if (stenb->getMicroPhaseId() == ELECTROLYTEID) {
          addGrowthSiteMod(stenb, phaseid);

          if (j < NUM_NEAREST_NEIGHBORS) {
            for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases;
                 phaseTmpl++) {
              if (chemSys_->isGrowthTemplate(phaseTmpl, phaseid)) {
                // addGrowthSite(stenb, nbgrowthtemp[jj]);
                addGrowthSiteMod(stenb, phaseTmpl);
              }
            }
          }

        } else if (stenbWmc == 0.0) {

          // removeDissolutionSite(stenb, stenb->getMicroPhaseId());
          removeDissolutionSiteMod_grow(stenb, stenb->getMicroPhaseId());
        }
      }

      numleft--;
      numchange++;
      //} else {
      //    cout << "error test_removeGrowthSite phaseid numtoadd trc trc_d bcl
      //    steId pid dim_isite numleft numchange  :  "
      //         << phaseid << "   " << numtoadd << "   " << trc << "   " <<
      //         trc_d << "   "
      //         << bcl << "   " << ste->getId() << "   " << pid << "   " <<
      //         dim_isite << "   "
      //         << numleft << "   " << numchange << endl;
      //    cout << endl << "stop program" << endl;
      //    exit(1);
      //    //removeGrowthSite(ste, phaseid);
      //}

    } catch (out_of_range &oor) {
      EOBException ex("Lattice", "dissolvePhase", "site_", site_.size(), j);
      ex.printException();
      cout << endl
           << "error grow phaseid numtoadd trc trc_d bcl steId pid dim_isite "
              "numleft numchange  :  "
           << phaseid << "   " << numtoadd << "   " << trc << "   " << trc_d
           << "   " << bcl << "   " << ste->getId() << "   " << pid << "   "
           << dim_isite << "   " << numleft << "   " << numchange << endl;
      cout.flush();
      exit(1);
    }

    isite = interface_[phaseid].getGrowthSites();
    dim_isite = isite.size();
  }

  // findInterfaces_check();
  return (numchange);
}

int Lattice::dissolvePhaseMod(unsigned int phaseid, int numtotake, int trc) {
  unsigned int i, ii, jj;

  double dwmcval;
  Site *ste, *stenb;
  vector<Isite> isite;
  isite = interface_[phaseid].getDissolutionSites();
  int numleft = numtotake;
  int numchange = 0;
  int isitePos;
  unsigned int pid;

  vector<unsigned int> growth_local;
  int grLocSize;
  int nbid;
  int nb_id, nb_pid; // for nn & nnn
  bool phaseid_exist;
  int numMicroPhases = chemSys_->getNumMicroPhases();

  int dim_isite = isite.size();
  struct dissProb {
    int id;
    int nbElectr;
    double instab;
    double prob;
    double prob_01;
  };
  vector<dissProb> dissProbVect;
  dissProb dissProb_;
  int stId;
  double rng;

  //*** for nn & nnn electrolyte sites
  // int sumEl, sumTot;
  // double sumTot_dbl;
  // double NN_NNN_dbl = NN_NNN; // NUM_NEAREST_NEIGHBORS;

  //*** for wmc
  double sumWmc, stWmc;

  //*** controll
  int bcl = 0;
  int static trc_d;
  trc_d++;

  if (numtotake == 0)
    return 0;
  try {
    if (phaseid >= interface_.size()) {
      throw EOBException("Lattice", "dissolvePhaseMod", "interface_",
                         interface_.size(), phaseid);
    }
  } catch (EOBException ex) {
    ex.printException();
    exit(1);
  }

  try {
    for (i = 0; i < isite.size(); i++) {
      if (site_.at(isite[i].getId()).getMicroPhaseId() != phaseid) {
        if (warning_) {
          cout << "Lattice::dissolvePhaseMod WARNING: Interface " << phaseid
               << " is corrupted with phase "
               << site_.at(isite[i].getId()).getMicroPhaseId() << endl;
          cout << "Lattice::dissolvePhaseMod Offending site is "
               << isite[i].getId() << endl;
          cout.flush();
        }
      }
    }
  } catch (out_of_range &oor) {
    EOBException ex("Lattice", "dissolvePhaseMod", "site_", site_.size(), i);
    ex.printException();
    exit(1);
  }

  if (verbose_) {
    cout << "Lattice::dissolvePhaseMod -->Phase " << phaseid
         << " needs to dissolve at " << numtotake << " sites" << endl;
  }

  // TEST DISSOLVE/GROW INTERFACES start
  cout << endl
       << endl
       << "********************************************" << endl;
  // cout<<endl<<"diss phaseid trc_d/bcl dim_isite numleft numchange  :  " <<
  // phaseid << "   " << trc_d << "   " << bcl <<
  //     "   " << dim_isite << "   " << numleft << "   " << numchange << endl;
  // cout.flush();

  // cout << endl << "     count_ before dissolution" <<endl;
  // for (i = 0 ; i < numMicroPhases; i++){
  //     cout << "  " << i << "     count_ = " << count_[i]
  //          << "     dissInterfaceSize =  " <<
  //          interface_[i].getDissolutionNumSites()
  //          << "     growInterfaceSize =  " <<
  //          interface_[i].getGrowthNumSites() << "       & is template for :
  //          ";
  //     for (jj = 0; jj < numMicroPhases; jj++){
  //         if (chemSys_->isGrowthTemplate(jj,i)){
  //             cout << jj << " ";
  //         }
  //     }
  //     cout << endl;
  // }
  // cout.flush();
  //  TEST DISSOLVE/GROW INTERFACES end

  cout << endl
       << "diss phaseid trc_d/bcl dim_isite numleft numchange  :  " << phaseid
       << "   " << trc_d << "   " << bcl << "   " << dim_isite << "   "
       << numleft << "   " << numchange << endl;
  cout.flush();

  while ((numleft > 0) && (dim_isite >= 1)) {
    try {
      bcl++;
      //    cout<<endl<<"phaseid trc_d/bcl dim_isite numleft numchange  :  " <<
      //    phaseid << "   " << trc_d << "   " << bcl <<
      //        "   " << dim_isite << "   " << numleft << "   " << numchange <<
      //        endl; cout.flush();

      /*
// dissolution probabilities based on the number of electrolyte nn & nnn
dissProbVect.clear();
sumTot = 0;
for (i = 0; i < dim_isite; i++) {
stId=isite[i].getId();
sumEl = 0;
for(jj = 0; jj< NN_NNN; jj++){
    if (site_[stId].nb(jj)->getMicroPhaseId() == 1){
    sumEl++;
  }
}
dissProb_.id = stId;
dissProb_.nbElectr = sumEl;
dissProb_.instab = sumEl/NN_NNN_dbl;
dissProb_.prob = 0.;
dissProb_.prob_01 = 0.;
dissProbVect.push_back(dissProb_);
sumTot += sumEl;
}
sumTot_dbl = sumTot;

//calc electrolyte neighbours probabilities
if (sumTot > 0) { // sumTot at least 1. !
dissProbVect[0].prob = (dissProbVect[0].nbElectr)/sumTot_dbl;
dissProbVect[0].prob_01 = dissProbVect[0].prob;
for (i = 1; i < dim_isite; i++) {
    dissProbVect[i].prob = (dissProbVect[i].nbElectr)/sumTot_dbl;
  dissProbVect[i].prob_01 = dissProbVect[i-1].prob_01 + dissProbVect[i].prob;
}

rng = RanGen::Ran3();
for (isitePos = 0; isitePos < dim_isite; isitePos++){
  if (rng <= dissProbVect[isitePos].prob_01) break;
}
} else {
  rng = RanGen::Ran3();
  isitePos = (int)(rng * dim_isite);
}
*/

      // dissolution probabilities based on wmc
      dissProbVect.clear(); // -> a simpler structure
      sumWmc = 0.;
      for (i = 0; i < dim_isite; i++) {
        stId = isite[i].getId();
        stWmc = site_[stId].getWmc();
        sumWmc += stWmc;
        dissProb_.id = stId;
        // dissProb_.nbElectr = 0;
        dissProb_.instab = stWmc;
        // dissProb_.prob = 0.;
        // dissProb_.prob_01 = 0.;
        dissProbVect.push_back(dissProb_);
      }

      // calc wmc probabilities
      // if (sumWmc > 0) {
      dissProbVect[0].prob = (dissProbVect[0].instab) / sumWmc;
      dissProbVect[0].prob_01 = dissProbVect[0].prob;
      for (i = 1; i < dim_isite; i++) {
        dissProbVect[i].prob = (dissProbVect[i].instab) / sumWmc;
        dissProbVect[i].prob_01 =
            dissProbVect[i - 1].prob_01 + dissProbVect[i].prob;
      }

      rng = RanGen::Ran3();
      for (isitePos = 0; isitePos < dim_isite; isitePos++) {
        if (rng <= dissProbVect[isitePos].prob_01)
          break;
      }
      //} else {
      //    cout << endl << "error in dissolution probabilities based on wmc" <<
      //    endl; cout << endl << "phaseid numtotake : " << phaseid << "   " <<
      //    numtotake << endl; cout << endl << "trc trc_d bcl dim_isite sumWmc :
      //    " << trc << "   " << trc_d << "   " << bcl
      //         << "   " << dim_isite  << "   " << sumWmc << endl << endl;
      //    for (i = 0; i < dim_isite; i++) {
      //        cout << "  dissProbVect[" << i << "] :" << endl;
      //        cout << "      dissProbVect.id       = " << dissProbVect[i].id
      //        << endl;
      //        //cout << "      dissProbVect.nbElectr = " <<
      //        dissProbVect[i].nbElectr << endl; cout << " dissProbVect.wmc = "
      //        << dissProbVect[i].instab << endl; cout << " dissProbVect.prob
      //        = " << dissProbVect[i].prob << endl; cout << "
      //        dissProbVect.prob_01  = " << dissProbVect[i].prob_01 << endl;
      //    }
      //    cout << "program stop" << endl;
      //    exit(1);
      //}

      ste = &site_.at(isite[isitePos].getId());
      pid = ste->getMicroPhaseId(); // intrebare pid diff phaseid ???

      double wmcIni, wmcEnd;
      wmcIni = ste->getWmc0();

      removeDissolutionSiteMod_diss(ste, pid, isitePos);
      setMicroPhaseIdMod_diss(ste, ELECTROLYTEID);

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

      // dwmcval = chemSys_->getMicroPhasePorosity(ELECTROLYTEID) -
      //          chemSys_->getMicroPhasePorosity(pid);

      // ste->dWmc(dwmcval);

      for (i = 0; i < NN_NNN; i++) {
        stenb = ste->nb(i);
        stenb->dWmc(dwmcval);
        int nbpid = stenb->getMicroPhaseId();
        if ((nbpid != ELECTROLYTEID) && (nbpid != VOIDID)) {

          ///
          /// Now that the site has been dissolved, it is eligible for growth
          /// later on, so we add it to the list of growth sites for certain
          /// phases.
          ///

          addDissolutionSiteMod(stenb, nbpid);

          addGrowthSiteMod(ste, nbpid);

          if (i < NUM_NEAREST_NEIGHBORS) {

            numMicroPhases = chemSys_->getNumMicroPhases();
            for (int phaseTmpl = FIRST_SOLID; phaseTmpl < numMicroPhases;
                 phaseTmpl++) {
              if (chemSys_->isGrowthTemplate(phaseTmpl, nbpid)) {

                addGrowthSiteMod(ste, phaseTmpl);
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
            if (phaseid_exist == false) { // verif removeGrowthSiteMod !!
              removeGrowthSiteMod1_grow(stenb, growth_local[ii]);
            }
          }
        }
      }

      numleft--;
      numchange++;

    } catch (out_of_range &oor) {
      EOBException ex("Lattice", "dissolvePhase", "site_", site_.size(), i);
      ex.printException();
      cout << endl
           << "error diss phaseid numtotake trc trc_d bcl steId pid dim_isite "
              "numleft numchange  :  "
           << phaseid << "   " << numtotake << "   " << trc << "   " << trc_d
           << "   " << bcl << "   " << ste->getId() << "   " << pid << "   "
           << dim_isite << "   " << numleft << "   " << numchange << endl;
      cout.flush();
      exit(1);
    }

    isite = interface_[phaseid].getDissolutionSites();
    dim_isite = isite.size();
  }

  // findInterfaces_check();
  // cout << "end first dissolution - stop" << endl; exit(1);
  return (numchange);
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

  for (i = 0; i < chemSys_->getNumMicroPhases(); i++) {
    if (i != ELECTROLYTEID && i != VOIDID) { // Solid phase of some kind
      gsite.clear();
      dsite.clear();
      for (k = 0; k < site_.size(); k++) {
        if (site_[k].getWmc() > 0) { // Is there some water nearby?
          if ((site_[k].getMicroPhaseId() == i)) {
            dsite.push_back(&site_[k]);
            site_[k].setDissolutionSiteMod(i);
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

  cout << endl
       << "***********   from findInterfaces_check   ************" << endl;

  bool differentInterfaceSize = false;
  for (int i = 0; i < chemSys_->getNumMicroPhases(); i++) {
    if (i == 2 || i == 3 || i == 4 || i == 5 || i == 7 || i == 8) {
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
    for (int i = 0; i < chemSys_->getNumMicroPhases(); i++) {
      cout << " phaseId" << endl;
      cout << "   " << i << endl;
      cout << "           dissolutionSites_.size()  = "
           << interface_[i].getDissolutionNumSites() << endl;
      cout << "           findInterfaces_check size = " << interfaceSizes[i].ds
           << endl;
    }
    cout << endl << "********** growth interfaces **********" << endl;
    for (int i = 0; i < chemSys_->getNumMicroPhases(); i++) {
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
    cout << "                      findInterfaces_check  OK!" << endl;
  }
}

void Lattice::addDissolutionSiteMod(Site *ste, unsigned int pid) {
  try {
    interface_.at(pid).addDissolutionSiteMod(ste);
    // vector<unsigned int> plist = ste->getGrowthPhases();
    // for (unsigned int i = 0; i < plist.size(); i++) {
    //     interface_.at(plist[i]).removeGrowthSite(ste);
    // }
    ste->setDissolutionSiteMod(pid);
  } catch (out_of_range &oor) {
    EOBException ex("Lattice", "addDissolutionSiteMod", "interface_",
                    interface_.size(), pid);
    ex.printException();
    cout << endl << "id pid " << ste->getId() << " " << pid << endl;
    exit(1);
  }
  return;
}

void Lattice::addGrowthSiteMod(Site *ste, unsigned int pid) {
  try {
    interface_.at(pid).addGrowthSiteMod(ste);
    ste->setGrowthSite(pid);
  } catch (out_of_range &oor) {
    EOBException ex("Lattice", "addGrowthSiteMod", "interface_",
                    interface_.size(), pid);
    ex.printException();
    cout << endl << "id pid " << ste->getId() << " " << pid << endl;
    exit(1);
  }
}

/*
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

void Lattice::removeDissolutionSiteMod_diss(Site *ste, unsigned int pid,
                                            int interfacePos) {
  try {
    interface_.at(pid).removeDissolutionSiteMod_diss(ste, interfacePos);
    ste->removeDissolutionSiteMod();
  } catch (out_of_range &oor) {
    EOBException ex("Lattice", "removeDissolutionSiteMod_diss", "interface_",
                    interface_.size(), pid);
    ex.printException();
    cout << endl
         << "id pid interfacePos " << ste->getId() << " " << pid << " "
         << interfacePos << endl;
    exit(1);
  }
}

void Lattice::removeDissolutionSiteMod_grow(Site *ste, unsigned int pid) {
  try {
    interface_.at(pid).removeDissolutionSiteMod_grow(ste);
    ste->removeDissolutionSiteMod();
  } catch (out_of_range &oor) {
    EOBException ex("Lattice", "removeDissolutionSiteMod_grow", "interface_",
                    interface_.size(), pid);
    ex.printException();
    cout << endl << "id pid " << ste->getId() << " " << pid << endl;
    exit(1);
  }
}

/*
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

void Lattice::removeGrowthSiteMod0_grow(Site *ste, unsigned int pid,
                                        int interfacePos) {
  try {
    interface_.at(pid).removeGrowthSiteMod0_grow(ste, interfacePos);
    ste->removeGrowthSiteMod_grow(pid);
  } catch (out_of_range &oor) {
    EOBException ex("Lattice", "removeGrowthSiteMod0_grow", "interface_",
                    interface_.size(), pid);
    ex.printException();

    cout << endl
         << "id pid interfacePos " << ste->getId() << " " << pid << " "
         << interfacePos << endl;
    exit(1);
  }
}

void Lattice::removeGrowthSiteMod1_grow(Site *ste, unsigned int pid) {
  try {
    interface_.at(pid).removeGrowthSiteMod1_grow(ste);
    ste->removeGrowthSiteMod_grow(pid);
  } catch (out_of_range &oor) {
    EOBException ex("Lattice", "removeGrowthSiteMod_grow", "interface_",
                    interface_.size(), pid);
    ex.printException();

    cout << endl << "id pid : " << ste->getId() << " " << pid << endl;
    exit(1);
  }
}

int Lattice::emptyPorosity(int numsites) {
  unsigned int i, j;
  int numemptied = 0;
  int maxsearchsize = 3;
  unsigned int cntpore, cntmax;
  bool placed;

  if (numsites == 0)
    return (0);
  if (numsites < 0) {
    numemptied = fillPorosity(-numsites);
    return (-numemptied);
  }

  ///
  /// Finding all potential VOID sites.
  ///
  /// @todo Consider removing some of the standard output, or setting a flag for
  /// it.
  ///

  list<Sitesize> distlist =
      findDomainSizeDistribution(ELECTROLYTEID, numsites, maxsearchsize, 0);
  list<Sitesize>::iterator it;

  if (verbose_) {
    cout << "Lattice::emptyPorosity Found " << distlist.size()
         << " potential void sites." << endl;
    cout.flush();
  }

  ///
  /// We want to empty the sites with the largest pore count
  ///

  if (distlist.size() < numsites) {
    string msg = "Ran out of water in the system";
    EOBException ex("Lattice", "emptysite", msg, distlist.size(), numsites);
    ex.printException();
    exit(1);
  }

  numemptied = 0;
  it = distlist.begin();
  int siteid;
  try {
    while (it != distlist.end()) {
      siteid = (*it).siteid;
      setMicroPhaseId(site_[siteid].getId(), VOIDID);
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

int Lattice::fillPorosity(int numsites) {
  unsigned int i, j;
  int numfilled = 0;
  int maxsearchsize = 10;
  unsigned int cntpore, cntmin;
  bool placed;

  if (numsites == 0)
    return (0);
  if (numsites < 0) {
    numfilled = emptyPorosity(-numsites);
    return (numfilled);
  }

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

  list<Sitesize> distlist =
      findDomainSizeDistribution(VOIDID, numsites, maxsearchsize, 1);
  list<Sitesize>::iterator it;

  ///
  /// We want to fill the sites with the smallest pore count
  ///

  numfilled = 0;
  it = distlist.begin();
  int siteid;
  try {
    while (it != distlist.end()) {
      siteid = (*it).siteid;
      setMicroPhaseId(site_[siteid].getId(), ELECTROLYTEID);
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
  double dist;

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
        dist = static_cast<double>(((xc - xp) * (xc - xp)));
        dist += static_cast<double>(((yc - yp) * (yc - yp)));
        dist += static_cast<double>(((zc - zp) * (zc - zp)));
        /*
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

void Lattice::changeMicrostructureMod(double time, const int simtype,
                                      bool isFirst, bool &capWater) {
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

  extern string CSHMicroName;
  extern string MonocarbMicroName;
  extern string MonosulfMicroName;
  extern string HydrotalcMicroName;
  extern string AFTMicroName;

  if (isFirst) {
    initialmicrostructurevolume_ = chemSys_->getInitMicroVolume();
  }

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

  /*
expansion_.clear();
expansion_coordin_.clear();
string fexpansion = jobroot_ + "_exp.dat";
string fexpansioncoor = jobroot_ + "_exp_coor.dat";
ifstream in(fexpansion.c_str());
if (!in) {
  cout << "can't open file " << fexpansion << ", so exit program." << endl;
  exit(1);
} else {
  while (!in.eof()) {
    int index;
    vector<double> expval;
    expval.clear();
    expval.resize(3,0.0);
    in >> index;
    //cout << "index = " << index << endl;
    in >> expval[0];
    in >> expval[1];
    in >> expval[2];
    //cout << "expval[0]: " << expval[0] << " expval[1]: " << expval[1]
         //<< " expval[2]: " << expval[2] << endl;
    expansion_.insert(make_pair(index,expval));
    site_[index].setExpansionStrain(expval[0]);
  }
  in.close();
}
cout << "open fexpansioncoor file" << endl;
ifstream in1(fexpansioncoor.c_str());
if (!in1) {
  cout << "can't open file " << fexpansioncoor << ", so exit program." <<
endl; exit(1); } else { while (!in1.eof()) { int index; vector<int> coordin;
    coordin.clear();
    coordin.resize(3,0);
    in1 >> index;
    in1 >> coordin[0];
    in1 >> coordin[1];
    in1 >> coordin[2];
    expansion_coordin_.insert(make_pair(index,coordin));
  }
  in1.close();
}

cout << "expansion_.size() = " << expansion_.size() << " at the beginning of "
     << time << endl;
*/

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

  if (verbose_) {
    cout << "Lattice::changeMicrostructure Before adjustMicrostructureVolumes:"
         << endl;
    for (int iii = 0; iii < vol_next.size(); ++iii) {
      cout << "Lattice::changeMicrostructure    Volume of " << phasenames[iii]
           << " = " << vol_next[iii] << " m3" << endl;
    }
    cout.flush();
  }

  try {
    adjustMicrostructureVolumes(phasenames, vol_next);

    if (verbose_) {
      cout << "Lattice::changeMicrostructure After adjustMicrostructureVolumes:"
           << endl;
      for (int iii = 0; iii < vol_next.size(); ++iii) {
        cout << "Lattice::changeMicrostructure   Volume of " << phasenames[iii]
             << " = " << vol_next[iii] << " m3" << endl;
      }
    }

    adjustMicrostructureVolFracs(phasenames, vol_next, vfrac_next);
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
      for (i = 0; i < vfrac_next.size(); i++) {
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

int numMicroPhases = chemSys_->getNumMicroPhases();

netsites.clear();
netsites.resize(numMicroPhases, 0);
pid.clear();
pid.resize(numMicroPhases, 0);

vector<int> growing;
vector<vector<int>> shrinking;
vector<vector<double>> volratios;
growing.clear();
shrinking.clear();
volratios.clear();
vector<int> idummy;
idummy.clear();
vector<double> ddummy;
ddummy.clear();

try {
  // Hard wiring for sulfate attack here
  // @todo Generalize to other phases

  growing.push_back(chemSys_->getMicroPhaseId(AFTMicroName));
  shrinking.resize(growing.size(), idummy);
  volratios.resize(growing.size(), ddummy);
  for (i = 0; i < growing.size(); ++i) {
    if (MonosulfMicroName.length() > 0) {
      shrinking[i].push_back(chemSys_->getMicroPhaseId(MonosulfMicroName));
      volratios[i].push_back(2.288);
    }
    if (MonocarbMicroName.length() > 0) {
      shrinking[i].push_back(chemSys_->getMicroPhaseId(MonocarbMicroName));
      volratios[i].push_back(2.699);
    }
    if (HydrotalcMicroName.length() > 0) {
      shrinking[i].push_back(chemSys_->getMicroPhaseId(HydrotalcMicroName));
      volratios[i].push_back(3.211);
    }
  }

  for (int ii = 0; ii < growing.size(); ++ii) {
    for (i = 0; i < vfrac_next.size(); i++) {
      cursites = (int)(count_.at(i) + 0.5);
      newsites = (int)((numsites_ * vfrac_next.at(i)) + 0.5);
      tnetsites = 0;
      if (i != ELECTROLYTEID && i != VOIDID)
        tnetsites = (newsites - cursites);
      netsites.at(i) = tnetsites;
      pid.at(i) = i;
      if (i == growing[ii] && isFirst) {
        netsites.at(i) = 0;
        count_.at(i) = newsites;
      }
      if (verbose_) {
        if (netsites.at(i) != 0) {
          cout << "Lattice::changeMicrostructure ****netsites["
               << phasenames.at(i) << "] in this state = " << netsites.at(i)
               << endl;
          cout.flush();
        }
      }
    }

    ///
    /// Next block gets executed only if we are now simulating
    /// sulfate attack.
    ///

    if (time_ >= sattack_time_) {
      if (verbose_) {
        cout << "Lattice::changeMicrostructure Crystal-pressure "
             << "transform at time_ = " << time_ << endl;
        cout.flush();
      }

      ///
      /// The relevant stress-free molar volume ratios for sulfate
      /// attack phase transformations.
      ///

      vector<int> numchanged;
      numchanged.clear();
      int growid = growing[ii];
      for (int iii = 0; iii < shrinking[ii].size(); ++iii) {
        numchanged.resize(2, 0);
        int shrinkid = shrinking[ii][iii];
        double volrat = volratios[ii][iii];
        if ((netsites.at(shrinkid) < 0) && (netsites.at(growid) > 0)) {
          numchanged = transform(shrinkid, netsites.at(shrinkid), growid,
                                 netsites.at(growid), volrat);

          netsites.at(shrinkid) += numchanged[0];
          netsites.at(growid) -= numchanged[1];
          if (verbose_) {
            cout << "Lattice::changeMicrostructure netsites.at(" << shrinkid
                 << ") is: " << netsites.at(shrinkid) << endl;
            cout << "Lattice::changeMicrostructure netsites.at(" << growid
                 << ") is: " << netsites.at(growid) << endl;
            cout.flush();
          }
        }
      }
    }
  }
}

catch (out_of_range &oor) {
  throw EOBException(
      "Lattice", "changeMicrostructure",
      "phasenames or count_ or pid or netsites or vfrac_next",
      phasenames.size(), i);
}
*/

  } else {

    ///
    /// Sulfate attack will NEVER be done during this simulation.
    ///  Normalize to get volume fractions and compute number
    ///  of sites of each phase needed.
    ///

    netsites.clear();
    netsites.resize(chemSys_->getNumMicroPhases(), 0);
    pid.clear();
    pid.resize(chemSys_->getNumMicroPhases(), 0);

    try {
      for (i = FIRST_SOLID; i < vfrac_next.size(); i++) {
        cursites = (int)(count_.at(i) + 0.5);
        newsites = (int)((numsites_ * vfrac_next.at(i)) + 0.5);
        netsites.at(i) = newsites - cursites;
        pid.at(i) = i;
        if ((verbose_) && (netsites.at(i) != 0)) {
          cout << "Lattice::changeMicrostructure ***netsites["
               << phasenames.at(i) << "] in this state = " << netsites.at(i)
               << "; cursites = " << cursites << " and newsites = " << newsites
               << endl;
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
  ///
  int netsitesSize_;
  netsitesSize_ = netsites.size();
  try {
    for (i = FIRST_SOLID; i < netsitesSize_; i++) {
      for (ii = i; ii < netsitesSize_; ii++) {
        if (netsites[ii] < netsites[i]) {
          tnetsites = netsites[i];
          netsites[i] = netsites[ii];
          netsites[ii] = tnetsites;
          tpid = pid.at(i);
          pid.at(i) = pid.at(ii);
          pid.at(ii) = tpid;
        }
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

  static int trc;
  trc++;
  cout << "netsites.size() for trc: " << netsites.size() << " " << trc << endl;
  int nuclei;
  unsigned int pid_;
  vector<Isite> gs, ds;

  try {
    for (i = FIRST_SOLID; i < netsitesSize_; i++) {
      pid_ = pid.at(i);
      numadded = 0;
      numadded_actual = 0;

      if (netsites[i] < 0) {
        if (verbose_) {
          cout << "Lattice::changeMicrostructure Going into "
               << "dissolve_phase now... pid = " << pid.at(i) << endl;
          cout.flush();
        }
        // numadded = dissolvePhase(pid.at(i), -netsites[i]);
        numadded = dissolvePhaseMod(pid.at(i), -netsites[i], trc);
        numadded_actual += numadded;
      } else if (netsites[i] > 0) {
        if (verbose_) {
          cout << "Lattice::changeMicrostructure Going into grow_phase "
               << "now... pid = " << pid.at(i) << endl;
          cout.flush();
        }

        numadded = growPhaseMod(pid.at(i), netsites[i], trc);

        numadded_actual += numadded;
        // int diff = netsites[i] - numadded;
        nuclei = netsites[i] - numadded;

        // while (diff > 0) {
        while (nuclei > 0) {
          gs = interface_[pid.at(i)].getGrowthSites();
          ds = interface_[pid.at(i)].getDissolutionSites();

          if (gs.size() == 0) {
            if (nuclei <= count_.at(ELECTROLYTEID)) {
              // create initial new-interface
              struct localStruct {
                int id;
                int aff;
                int count;
                double prob;
              };
              localStruct un;
              vector<localStruct> watersites;
              int aff, affMin, affSum;
              Site *nbId;
              vector<Site *> localNb_;
              int sizeWS_;
              int valAbs;
              int dimNewInterface;
              double rng;
              int fSiteWS;
              int maxFindSameSite;

              affSum = 0;
              affMin = 1000000;
              maxFindSameSite = 3;
              watersites.clear();
              for (int k = 0; k < numsites_; ++k) {
                if (site_[k].getMicroPhaseId() == ELECTROLYTEID) {
                  un.id = site_[k].getId();
                  aff = 0;
                  localNb_ = site_[k].getNb();
                  for (int j = 0; j < NN_NNN; j++) {
                    aff += chemSys_->getAffinity(
                        pid_, localNb_[j]->getMicroPhaseId());
                  }
                  un.aff = aff;
                  un.count = 0;
                  un.prob = 0.0;
                  watersites.push_back(un);
                  affSum += aff;
                  if (aff < affMin)
                    affMin = aff;
                }
              }

              // first probability

              sizeWS_ = watersites.size();
              if (affMin < 0) {
                valAbs = abs(affMin);
                affSum += sizeWS_ * valAbs;
                watersites[0].prob =
                    (watersites[0].aff + valAbs) / (double)affSum;
                for (int k = 1; k < sizeWS_; k++) {
                  watersites[k].prob =
                      watersites[k - 1].prob +
                      (watersites[k].aff + valAbs) / (double)affSum;
                }
              } else {
                watersites[0].prob = (watersites[0].aff) / (double)affSum;
                for (int k = 1; k < sizeWS_; k++) {
                  watersites[k].prob = watersites[k - 1].prob +
                                       (watersites[k].aff) / (double)affSum;
                }
              }

              dimNewInterface = 0;
              if (affSum != 0) {
                // not tested cauta
                while (dimNewInterface < nuclei) {
                  rng = RanGen::Ran3();
                  for (fSiteWS = 0; fSiteWS < sizeWS_; fSiteWS++) {
                    if (rng <= watersites[fSiteWS].prob)
                      break;
                  }
                  if (watersites[fSiteWS].count == 0) {
                    watersites[fSiteWS].count++;
                    site_[watersites[fSiteWS].id].setGrowthSite(pid_);
                    interface_.at(pid_).addGrowthSiteMod_newInterface(
                        watersites[fSiteWS].id, watersites[fSiteWS].aff);
                    dimNewInterface++;
                  } else if (watersites[fSiteWS].count > maxFindSameSite) {
                    if (affMin < 0) {
                      affSum = affSum - valAbs - watersites[fSiteWS].aff;
                      watersites[fSiteWS] = watersites[sizeWS_ - 1];
                      watersites.pop_back();
                      sizeWS_--;
                      watersites[0].prob =
                          (watersites[0].aff + valAbs) / (double)affSum;
                      for (int k = 1; k < sizeWS_; k++) {
                        watersites[k].prob =
                            watersites[k - 1].prob +
                            (watersites[k].aff + valAbs) / (double)affSum;
                      }
                    } else {
                      affSum = affSum - watersites[fSiteWS].aff;
                      watersites[fSiteWS] = watersites[sizeWS_ - 1];
                      watersites.pop_back();
                      sizeWS_--;
                      watersites[0].prob = (watersites[0].aff) / (double)affSum;
                      for (int k = 1; k < sizeWS_; k++) {
                        watersites[k].prob =
                            watersites[k - 1].prob +
                            (watersites[k].aff) / (double)affSum;
                      }
                    }
                  }
                  // add this site to dissolution list and extract it from
                  // growth list
                }
              } else {
                // tested and apparently works  cauta
                while (dimNewInterface < nuclei) {
                  rng = RanGen::Ran3();
                  fSiteWS = (int)(rng * sizeWS_);
                  if (watersites[fSiteWS].count == 0) {
                    watersites[fSiteWS].count++;
                    site_[watersites[fSiteWS].id].setGrowthSite(pid_);
                    interface_.at(pid_).addGrowthSiteMod_newInterface(
                        watersites[fSiteWS].id, watersites[fSiteWS].aff);
                    dimNewInterface++;
                  } else if (watersites[fSiteWS].count > maxFindSameSite) {
                    watersites[fSiteWS] = watersites[sizeWS_ - 1];
                    watersites.pop_back();
                    sizeWS_ = watersites.size();
                  }
                }
              }

              // Now the list is thoroughly shuffled, we just pick
              // the first nuclei elements as new growth sites

              //  for (int k = 0; k < nuclei; ++k) {
              //    addGrowthSite(&site_[watersites[k]], pid.at(i));
              //  }
            } else {
              if (verbose_ || warning_) {
                cout << "Lattice::changeMicrostructure There is no "
                     << "room to grow, so exit the program." << endl;
                cout.flush();
              }
              bool is_Error = false;
              throw MicrostructureException("Lattice",
                                            "changeMicrostructureMod",
                                            "no room to grow phase", is_Error);
            }
          } // gs.size()==0
          // numadded = growPhase(pid.at(i), diff);
          cout << "in second growPhaseMod i pid netsites " << i << " "
               << pid.at(i) << " " << netsites[i] << endl;
          numadded = growPhaseMod(pid.at(i), nuclei, trc);
          cout << "   out2 nuclei numadded count_.at(ELECTROLYTEID) gs.size() "
               << nuclei << " " << numadded << " " << count_.at(ELECTROLYTEID)
               << " " << gs.size() << endl;
          numadded_actual += numadded;
          // diff = diff - numadded;
          nuclei = nuclei - numadded;
        }
      }

      if (numadded_actual * numadded_actual != netsites[i] * netsites[i]) {
        if (warning_) {
          cout << "Lattice::changeMicrostructure WARNING: Needed to switch on "
               << netsites[i] << " of phase " << pid.at(i) << endl;
          cout << "Lattice::changeMicrostructure          But actually did "
               << numadded << " switches" << endl;
          cout.flush();
        }
      }
    }
  }

  catch (out_of_range &oor) {
    throw EOBException("Lattice", "changeMicrostructureMod", "pid", pid.size(),
                       i);
  }

  cursites = count_.at(VOIDID);
  newsites = (int)((numsites_ * vfrac_next.at(VOIDID)) + 0.5);
  wcursites = count_.at(ELECTROLYTEID);
  wnewsites = wcursites - (newsites - cursites);
  int numempty = emptyPorosity(newsites - cursites);
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
         << "actually emptied was:  " << numempty << endl;

    ///
    /// Report on target and actual mass fractions
    ///

    cout << "Lattice::changeMicrostructure "
         << "*******************************" << endl;
    cout.flush();
  }

  try {
    double totcount = 0.0;
    for (i = 0; i < vfrac_next.size(); i++) {
      volumefraction_.at(i) =
          ((double)(count_.at(i))) / ((double)(site_.size()));
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
  }

  catch (out_of_range &oor) {
    throw EOBException("Lattice", "changeMicrostructure",
                       "volumefraction_ or count_", volumefraction_.size(), i);
  }

  //    for (i = FIRST_SOLID; i < interface_.size(); i++) {
  //        interface_[i].sortGrowthSites(site_, i);
  //        interface_[i].sortDissolutionSites(site_, i);
  //    }

  ///  This is a local variable and the value is never used.
  ///
  ///  @todo Why not eliminate this line completely?
  ///

  double surfa = getSurfaceArea(chemSys_->getMicroPhaseId(CSHMicroName));

  if (volumefraction_[ELECTROLYTEID] <= 0.0)
    capWater = false;

  return;
}

void Lattice::adjustMicrostructureVolumes(vector<string> names,
                                          vector<double> &vol) {
  int i = 0;

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

    calcCapillarywatervolume(vol);
    calcCapillaryvoidvolume(vol);
    watervolume_ = vol.at(ELECTROLYTEID);

    calcSolidvolumewithpores(vol);

    if (solidvolumewithpores_ <= 0.0) {
      throw DataException("Lattice", "adjustMicrostructureVolumes",
                          "totvolume is NOT positive");
    }

    // The current microstructure volume as predicted by GEMS
    // The initial microstructure volume is calculated the first
    // time that Lattice::changeMicrostructure is called
    // We currently insist remain the volume of the system, so if
    // the microstructure volume deviates from the system volume,
    // we add or subtract capillary porosity to keep them equal

    microstructurevolume_ = chemSys_->getMicroVolume();

    // The total amount of non-solid space in the microstructure
    // This is based on forcing the total volume to be constant
    // and then calculating how much total material volume there is.

    nonsolidvolume_ =
        microstructurevolume_ - (solidvolumewithpores_ - subvoxelporevolume_);

    /// @note BULLARD trying to assign any change in total microstructure
    /// volume to capillary porosity.  Try to do this with one line here

    nonsolidvolume_ += (initialmicrostructurevolume_ - microstructurevolume_);

    calcCapillaryvoidvolume(vol);

    voidvolume_ = nonsolidvolume_ - watervolume_;

    // Up to this point, vol.at(ELECTROLYTEID) is the total volume
    // of pore solution in the system regardless of whether it is
    // capillary or subvoxel

    // First need to trim off any extra water that lies outside the system

    voidvolume_ = nonsolidvolume_ - watervolume_;
    if (voidvolume_ < 0.0)
      watervolume_ = nonsolidvolume_;

    capillarywatervolume_ = watervolume_ - subvoxelporevolume_;
    capillaryvoidvolume_ = capillaryporevolume_ - capillarywatervolume_;
    subvoxelwatervolume_ = subvoxelporevolume_;

    if (capillaryvoidvolume_ < 0.0)
      capillaryvoidvolume_ = 0.0;
    if (capillarywatervolume_ < 0.0) {
      subvoxelwatervolume_ = watervolume_;
      if (subvoxelwatervolume_ < 0.0)
        subvoxelwatervolume_ = 0.0;
      capillarywatervolume_ = 0.0;
    }

    if (chemSys_->isSaturated()) { // System is saturated

      voidvolume_ = capillaryvoidvolume_ = 0.0;
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

    ///
    /// End of manual adjustment
    ///

  }

  catch (DataException dex) {
    throw dex;
  }

  catch (out_of_range &oor) {
    throw EOBException("Lattice", "adjustMicrostructureVolumes", "names",
                       names.size(), i);
  }

  return;
}

void Lattice::adjustMicrostructureVolFracs(vector<string> &names,
                                           const vector<double> vol,
                                           vector<double> &vfrac) {
  int i = 0;
  double totmicvolume = 0.0;

#ifdef DEBUG
  cout << "Lattice::adjustMicrostructureVolFracs" << endl;
  cout.flush();
#endif

  try {

    // Remember there are now two extra slots at the end of vfrac, just
    // like at the end of vol.  These two extra slots hold the capillary
    // pore space and the subvoxel pore space

    vfrac.clear();
    vfrac.resize(vol.size(), 0.0);

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

    totmicvolume = chemSys_->getInitMicroVolume();

#ifdef DEBUvoid
    cout << "Lattice::adjustMicrostructureVolFracsCalculated "
         << "total microstructure volume is " << totmicvolume << endl;
    cout.flush();
#endif

    // Calculate volume fractions based on total GEMS adjusted volume

    for (i = 0; i < vol.size(); ++i) {
      vfrac.at(i) = vol.at(i) / totmicvolume;
      if (verbose_) {
        cout << "Lattice::adjustMicrostructureVolFracsVolume "
             << "fraction[" << names.at(i) << "] should be " << vfrac.at(i)
             << ", (" << vol.at(i) << "/" << totmicvolume
             << ") and volume fraction NOW is "
             << (double)(count_.at(i)) / (double)(numsites_) << endl;
        cout.flush();
      }
    }

    // Calculate volume fraction of subvoxel porosity and
    // capillary porosity whether saturated or not

    capillaryporevolumefraction_ = capillaryporevolume_ / totmicvolume;
    subvoxelporevolumefraction_ = subvoxelporevolume_ / totmicvolume;

  }

  catch (DataException dex) {
    dex.printException();
    exit(1);
  }

  catch (out_of_range &oor) {
    EOBException ex("Lattice", "adjustMicrostructureVolFracs", "vol",
                    vol.size(), i);
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
    double phi; // Holds the subvoxel porosity of a microstructurephase
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
    cumulative_volume = 0.0;
    cumulative_volfrac = 0.0;
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

  double microvol = chemSys_->getMicroVolume();

  // microvol is the absolute initial volume (m3) of all the defined
  // microstructure phases, including subvoxel pores in
  // solid phases

  double initmicrovol = chemSys_->getInitMicroVolume();

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

  bool isfullysaturated = false;
  if (excesswater > 0.0)
    isfullysaturated = true;

  double water_volfrac = water_volume / microstructurevolume_;

  // This is the total porosity including capillary
  // pore volume fraction and subvoxel pore volume
  // fraction, all on a total microstructure volume
  // basis

  double pore_volfrac =
      capillaryporevolumefraction_ + subvoxelporevolumefraction_;

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
    cout << endl;
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

  double water_volume = chemSys_->getMicroPhaseVolume(ELECTROLYTEID);

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

  bool isfullysaturated = false;
  if (excesswater > 0.0)
    isfullysaturated = true;

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

  int numMicroPhases = chemSys_->getNumMicroPhases();
  int microPhaseId;
  double red, green, blue;
  vector<double> colors;
  out << numMicroPhases << endl;
  for (int i = 0; i < numMicroPhases; i++) {
    microPhaseId = chemSys_->getMicroPhaseId(i);
    colors = chemSys_->getColor(microPhaseId);
    out << chemSys_->getMicroPhaseName(microPhaseId) << "," << colors[0] << ","
        << colors[1] << "," << colors[2] << endl;
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
  string ofbasename(root);
  string ofpngName(root);
  string ofpngbasename(root);
  bitmap_t image;

  image.pixels = NULL;

  vector<double> dumvec;
  vector<unsigned int> idumvec;
  vector<vector<unsigned int>> rawimage;
  vector<vector<double>> dshade;
  dumvec.resize(ydim_, 0.0);
  idumvec.resize(ydim_, 0);
  dshade.resize(xdim_, dumvec);
  rawimage.resize(xdim_, idumvec);
  bool done;

  image.width = xdim_;
  image.height = ydim_;
  image.pixels = pixelvector(image.width * image.height);

  ///
  /// Construct the name of the output file
  ///

  ostringstream ostr1, ostr2;
  ostr1 << setfill('0') << setw(6)
        << static_cast<int>((curtime * 24.0 * 60.0) + 0.5); // minutes
  ostr2 << setprecision(3) << temperature_;
  string timestr(ostr1.str());
  string tempstr(ostr2.str());
  string buff;
  ofpngName = ofpngName + "." + timestr + "m." + tempstr + ".png";

  ///
  /// Write PNG full color image
  ///

  unsigned int slice = zdim_ / 2;
  unsigned int nd, izz, valout;
  unsigned int sitenum;
  for (j = 0; j < ydim_; j++) {
    for (i = 0; i < xdim_; i++) {
      if (deptheffect_) {
        done = false;
        nd = 0;
        sitenum = getIndex(i, j, slice);
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
        rawimage[i][j] = site_[sitenum].getMicroPhaseId();
        dshade[i][j] = 0.1 * (10.0 - (static_cast<double>(nd)));
      } else {
        sitenum = getIndex(i, j, slice);
        rawimage[i][j] = site_[sitenum].getMicroPhaseId();
        dshade[i][j] = 1.0;
      }
    }
  }

  vector<double> colors;
  for (j = 0; j < ydim_; j++) {
    for (i = 0; i < xdim_; i++) {
      colors = chemSys_->getColor(rawimage[i][j]);
      pixel_t *pixel = pixel_at(&image, i, j);
      pixel->red = dshade[i][j] * colors[0] + 0.5;
      pixel->green = dshade[i][j] * colors[1] + 0.5;
      pixel->blue = dshade[i][j] * colors[2] + 0.5;
    }
  }

  save_png_to_file(&image, ofpngName.c_str());

  ///
  /// Execute system call to convert PPM to PNG.
  ///
  /// @warning This relies on installation of ImageMagick
  ///

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
      sitenum = getIndex(i, j, slice);
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
        sitenum = getIndex(i, j, slice);
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

/*
vector<int> Lattice::transform(int shrinkingid, int netsites_shrinkingid,
                               int growingid, int netsites_growingid,
                               double volumeratio) {
  ///
  /// @todo Consider breaking this method into smaller pieces
  ///

  int expindex;

  vector<double> expval;
  expval.clear();
  expval.resize(3, 0.0);

  vector<int> coordin;
  coordin.clear();
  coordin.resize(3, 0);

  vector<Isite> diss;
  diss = interface_[shrinkingid].getDissolutionSites();

  Site *ste;

  ///
  /// Construct the unordered list of sites to dissolve based only
  /// on the phase id and whether or not the total number of sites to dissolve
  /// has been reached.
  ///
  /// @remark Is this task biased by site position?  Has the list of sites been
  /// randomized?
  ///

  if (diss.size() < (-netsites_shrinkingid)) {
    for (int ii = 0; ii < numsites_; ii++) {
      ste = &site_[ii];
      if (ste->getMicroPhaseId() == shrinkingid) {
        addDissolutionSite(ste, shrinkingid);
      }
    }
  }

  diss = interface_[shrinkingid].getDissolutionSites();

  double alreadygrown = 0.0;
  int numtransform = 0;
  int max = (int)volumeratio;

  for (int ii = (diss.size() - 1);
       ii > 0 && alreadygrown < netsites_growingid &&
       numtransform < (-netsites_shrinkingid);
       ii--) {
    ste = &site_[diss[ii].getId()];

    expval.clear();

    vector<Site *> porousneighbor, waterneighbor;
    porousneighbor.clear();
    waterneighbor.clear();
    for (int j = 0; j < ste->nbSize(1); j++) {
      Site *stenb;
      stenb = ste->nb(j);
      if (stenb->getMicroPhaseId() == ELECTROLYTEID) {
        waterneighbor.push_back(stenb);
      } else if (chemSys_->isPorous(stenb->getMicroPhaseId())) {
        porousneighbor.push_back(stenb);
      }
    }

    if ((waterneighbor.size() + 1) <= max) { // count the site itself

      /// Expansion should occur
      ///
      /// 1. Take the subvolume centered on aluminate site that will dissolve
      ///

      string fileName(jobroot_ + "_alsubvol.dat");
      vector<unsigned int> alnb = writeSubVolume(fileName, ste, 1);
      ste->setDamage();
      int numWater, numPorous;
      numWater = numPorous = 0;
      for (int nb = 0; nb < alnb.size(); nb++) {
        Site *alstenb = &site_[alnb[nb]];
        if (alstenb->getMicroPhaseId() == ELECTROLYTEID) {
          numWater++;
        } else if (chemSys_->isPorous(alstenb->getMicroPhaseId())) {
          numPorous++;
          alstenb->setDamage();
        } else if (chemSys_->isWeak(alstenb->getMicroPhaseId())) {
          alstenb->setDamage();
        }
      }

      ///
      /// 2. Calculate the effective bulk modulus of this subvolume
      ///

      double subbulk = FEsolver_->getBulkModulus(fileName);

      subbulk = subbulk * 1.0e3; // convert GPa to MPa
      double subsolidbulk = subbulk;
      subsolidbulk *=
          ((1.0 + (double)numWater / 27.0) / (1.0 - (double)numWater / 27.0));

      ///
      /// @remark This seems like a double conversion.  Hasn't the conversion
      /// been done?
      ///

      subsolidbulk = subsolidbulk * 1.0e3; // convert GPa to MPa

      ///
      /// 3. Calculate crystallization strain in this sub volume;
      ///    porevolfrac is the volume fraction of pore space occupied by
      ///    crystal
      ///
      ///    @todo generalize porous phase porosities in the block below instead
      ///    of 0.25
      ///

      double porevolfrac = 0.0;
      if (numWater != 0 || numPorous != 0) {
        porevolfrac = (double)(volumeratio) / (numWater + (numPorous * 0.25));
      } else {
        porevolfrac = 1.0;
      }

      //  This is hard-wired right now
      //  @todo generalize crystallization pressure to more phases

      double exp = solut_->calculateCrystalStrain(SI_[growingid], porevolfrac,
                                                  subbulk, subsolidbulk);

      ///
      /// 4. Apply expansion strain on each voxel in this sub volume
      ///

      applyExp(alnb, exp);

      setMicroPhaseId(ste, growingid);

      ///
      /// The Al-bearing phase has dissolved, so remove it from the
      /// list of dissolution sites of this phase
      ///

      removeDissolutionSite(ste, shrinkingid);
      numtransform++;
      alreadygrown++;
      count_.at(growingid)++;

      for (int i = 0; i < waterneighbor.size(); i++) {
        setMicroPhaseId(waterneighbor[i], growingid);

        ///
        /// Ettringite has grown here, so remove this site from the
        /// list of growth sites of this phase
        ///

        removeGrowthSite(waterneighbor[i], growingid);
        count_.at(growingid)++;
        alreadygrown++;

        ///
        /// Weighted mean curvature (wmc) is changed by the difference
        /// between the growing phase's porosity and the template's porosity.
        ///
        /// @todo Determine why the calculation works this way.
        ///

        double dwmcval = chemSys_->getMicroPhasePorosity(growingid) -
                         chemSys_->getMicroPhasePorosity(ELECTROLYTEID);
        for (int j = 0; j < waterneighbor[i]->nbSize(2); j++) {
          Site *nb = waterneighbor[i]->nb(j);
          nb->dWmc(dwmcval);
        }
      }

      dWaterchange(volumeratio - (waterneighbor.size() + 1));
      alreadygrown += (volumeratio - (waterneighbor.size() + 1));
      count_.at(growingid) +=
          (int)(volumeratio - (waterneighbor.size() + 1) + 0.5);

    } else {

      ///
      /// Expansion should not occur because there is sufficient
      /// free space for local ettringite growth.
      ///

      setMicroPhaseId(ste, growingid);

      ///
      /// The Al-bearing phase has dissolved, so remove it from the
      /// list of dissolution sites of this phase
      ///

      removeDissolutionSite(ste, shrinkingid);
      numtransform++;
      alreadygrown++;
      count_.at(growingid)++;
      double thresh = 0.0, g = 0.0;
      thresh = volumeratio - max;
      g = rg_->Ran3();

      int upperindex = (g < thresh) ? max : max - 1;
      for (int i = 0; i < upperindex; i++) {
        setMicroPhaseId(waterneighbor[i], growingid);
        removeGrowthSite(waterneighbor[i], growingid);
        alreadygrown++;
        count_.at(growingid)++;

        ///
        /// Weighted mean curvature (wmc) is changed by the difference
        /// between the growing phase's porosity and the template's porosity.
        ///
        /// @todo Determine why the calculation works this way.
        ///

        double dwmcval = chemSys_->getMicroPhasePorosity(growingid) -
                         chemSys_->getMicroPhasePorosity(ELECTROLYTEID);
        for (int j = 0; j < waterneighbor[i]->nbSize(2); j++) {
          Site *nb = waterneighbor[i]->nb(j);
          nb->dWmc(dwmcval);
        }
      }
    }
  } // End of loop over all ettringite sites to form

  //netsites.at(shrinkingid) += (int) numtransform;
  //netsites.at(growingid) -= (int) alreadygrown;

  vector<int> numchanged;
  numchanged.clear();
  numchanged.resize(2, 0);
  numchanged[0] = numtransform;
  numchanged[1] = (int)alreadygrown;

  return numchanged;
}

*/

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

  for (int i = 0; i < isite.size(); i++) {
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
  list<Sitesize> distlist;
  list<Sitesize>::iterator it;
  bool found = false;

  Sitesize ss;
  int domainsize = 0;
  for (int i = 0; i < site_.size(); i++) {
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
