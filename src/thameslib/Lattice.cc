/**
@file Lattice.cc
@brief Definition of methods for the Lattice class.

*/
#include "Lattice.h"
#include "Interface.h"

Lattice::Lattice (ChemicalSystem *cs,
                  Solution *solut)
: siteneighbors_(18),chemSys_(cs),solut_(solut)
{
  xdim_ = ydim_ = zdim_ = 0;
  time_ = 0.0;
  temperature_ = REFTEMP;   // in Kelvin (see global.h)
  oldtemp_ = REFTEMP;       // in Kelvin (see global.h)
  numsites_ = 0;
  resolution_ = REFRES;     // in micrometers (see global.h)
  site_.clear();
  deptheffect_ = false;
  verbose_ = false;
}

Lattice::Lattice (ChemicalSystem *cs,
                  Solution *solut,
                  const string &fileName,
                  const bool verbose,
                  const bool warning,
                  const bool debug)
: siteneighbors_(18),chemSys_(cs),solut_(solut)
{
  unsigned int i,j,k;
  unsigned int ii;
  string buff;
  int xn,yn,zn;
  unsigned int idn;
  unsigned int pid;
  string msg;
    
  xdim_ = ydim_ = zdim_ = 0;
  time_ = 0.0;
  temperature_ = REFTEMP;   // in Kelvin (see global.h)
  oldtemp_ = REFTEMP;       // in Kelvin (see global.h)
  numsites_ = 0;
  resolution_ = REFRES;     // in micrometers (see global.h)
  site_.clear();
  deptheffect_ = true;
    
  verbose_ = verbose;
  debug_ = debug;
  warning_ = warning;

  ///
  /// Open the microstructure input file and process it
  ///

  ifstream in(fileName.c_str());
  if (!in) {
      throw FileException("Lattice","Lattice",fileName,"Could not open");
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
      cout << "Read microstructure file header..." << endl;
      cout.flush();
      cout << "    Version = " << version_ << endl;
      cout.flush();
      cout << "    xdim_ = " << xdim_ << endl;
      cout.flush();
      cout << "    ydim_ = " << ydim_ << endl;
      cout.flush();
      cout << "    zdim_ = " << zdim_ << endl;
      cout.flush();
  }
  numsites_ = (unsigned int)(xdim_ * ydim_ * zdim_);
  if (verbose_) {
      cout << "    numsites_ = " << numsites_ << endl;
      cout.flush();
      cout << "    resolution_ = " << resolution_ << endl;
      cout.flush();  
  }
  
  ///
  /// Allocate a random number generator object and seed it
  ///

  try { rg_ = new RanGen(); }
    catch (bad_alloc &ba) {
        cout << "Lattice constructor failed when allocating rg_";
        cout.flush();
        exit(1);
    }

    rg_->setSeed(-142234);
    if (verbose_) {
        cout << "Checking whether I can get a seed " << rg_->getSeed() << endl;
        cout.flush();
        cout << "Checking value of random number " << ", " << rg_->Ran3() << endl;
        cout.flush();
    }

    count_.clear();
    count_.resize(chemSys_->getNumMicroPhases(),0);

    SI_.clear();
    SI_.resize(chemSys_->getNumMicroPhases(),1.0);

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
                addSite(i,j,k);
            }
        }
    }

    ///
    /// Set up the array of neighbors for each site
    ///

    for (ii = 0; ii < numsites_; ii++) {

      for (j = 0; j < (NUM_NEAREST_NEIGHBORS
                + NUM_SECONDNEAREST_NEIGHBORS); j++) {
        ///
        /// Store up to 2nd nearest neighbors
        ///

        switch (j) {

            case 0: // West neighbor
                xn = (int)site_[ii].getX() - 1;
                if (xn < 0 && BC != 1) xn += (int)xdim_;
                if (xn < 0 && BC == 1) xn = 0;
                yn = (int)site_[ii].getY();
                zn = (int)site_[ii].getZ();
                break;
            case 1: // East neighbor
                xn = (int)site_[ii].getX() + 1;
                if (xn >= (int)xdim_ && BC != 1) xn -= (int)xdim_;
                if (xn >= (int)xdim_ && BC == 1) xn -= 1;
                yn = (int)site_[ii].getY();
                zn = (int)site_[ii].getZ();
                break;
            case 2: // South neighbor
                yn = (int)site_[ii].getY() - 1;
                if (yn < 0 && BC != 2) yn += (int)ydim_;
                if (yn < 0 && BC == 2) yn = 0;
                xn = (int)site_[ii].getX();
                zn = (int)site_[ii].getZ();
                break;
            case 3: // North neighbor
                yn = (int)site_[ii].getY() + 1;
                if (yn >= (int)ydim_ && BC != 2) yn -= (int)ydim_;
                if (yn >= (int)ydim_ && BC == 2) yn -= 1;
                xn = (int)site_[ii].getX();
                zn = (int)site_[ii].getZ();
                break;
            case 4: // Down neighbor
                zn = (int)site_[ii].getZ() - 1;
                if (zn < 0 && BC != 3) zn += (int)zdim_;
                if (zn < 0 && BC == 3) zn = 0;
                xn = (int)site_[ii].getX();
                yn = (int)site_[ii].getY();
                break;
            case 5: // Up neighbor
                zn = (int)site_[ii].getZ() + 1;
                if (zn >= (int)zdim_ && BC != 3) zn -= (int)zdim_;
                if (zn >= (int)zdim_ && BC == 3) zn -= 1;
                xn = (int)site_[ii].getX();
                yn = (int)site_[ii].getY();
                break;
            case 6: // Southwest neighbor
                xn = (int)site_[ii].getX() - 1;
                if (xn < 0 && BC != 1) xn += (int)xdim_;
                if (xn < 0 && BC == 1) xn = 0;
                yn = (int)site_[ii].getY() - 1;
                if (yn < 0 && BC != 2) yn += (int)ydim_;
                if (yn < 0 && BC == 2) yn = 0;
                zn = (int)site_[ii].getZ();
                break;
            case 7: // Northwest neighbor
                xn = (int)site_[ii].getX() - 1;
                if (xn < 0 && BC != 1) xn += (int)xdim_;
                if (xn < 0 && BC == 1) xn = 0;
                yn = (int)site_[ii].getY() + 1;
                if (yn >= (int)ydim_ && BC != 2) yn -= (int)ydim_;
                if (yn >= (int)ydim_ && BC == 2) yn -= 1;
                zn = (int)site_[ii].getZ();
                break;
            case 8: // Northeast neighbor
                xn = (int)site_[ii].getX() + 1;
                if (xn >= (int)xdim_ && BC != 1) xn -= (int)xdim_;
                if (xn >= (int)xdim_ && BC == 1) xn -= 1;
                yn = (int)site_[ii].getY() + 1;
                if (yn >= (int)ydim_ && BC != 2) yn -= (int)ydim_;
                if (yn >= (int)ydim_ && BC == 2) yn -= 1;
                zn = (int)site_[ii].getZ();
                break;
            case 9: // Southeast neighbor
                xn = (int)site_[ii].getX() + 1;
                if (xn >= (int)xdim_ && BC != 1) xn -= (int)xdim_;
                if (xn >= (int)xdim_ && BC == 1) xn -= 1;
                yn = (int)site_[ii].getY() - 1;
                if (yn < 0 && BC != 2) yn += (int)ydim_;
                if (yn < 0 && BC == 2) yn = 0;
                zn = (int)site_[ii].getZ();
                break;
            case 10: // Downwest neighbor
                xn = (int)site_[ii].getX() - 1;
                if (xn < 0 && BC != 1) xn += (int)xdim_;
                if (xn < 0 && BC == 1) xn = 0;
                yn = (int)site_[ii].getY();
                zn = (int)site_[ii].getZ() - 1;
                if (zn < 0 && BC != 3) zn += (int)zdim_;
                if (zn < 0 && BC == 3) zn = 0;
                break;
            case 11: // Downnorth neighbor
                xn = (int)site_[ii].getX();
                yn = (int)site_[ii].getY() + 1;
                if (yn >= (int)ydim_ && BC != 2) yn -= (int)ydim_;
                if (yn >= (int)ydim_ && BC == 2) yn -= 1;
                zn = (int)site_[ii].getZ() - 1;
                if (zn < 0 && BC != 3) zn += (int)zdim_;
                if (zn < 0 && BC == 3) zn = 0;
                break;
            case 12: // Downeast neighbor
                xn = (int)site_[ii].getX() + 1;
                if (xn >= (int)xdim_ && BC != 1) xn -= (int)xdim_;
                if (xn >= (int)xdim_ && BC == 1) xn -= 1;
                yn = (int)site_[ii].getY();
                zn = (int)site_[ii].getZ() - 1;
                if (zn < 0 && BC != 3) zn += (int)zdim_;
                if (zn < 0 && BC == 3) zn = 0;
                break;
            case 13: // Downsouth neighbor
                xn = (int)site_[ii].getX();
                yn = (int)site_[ii].getY() - 1;
                if (yn < 0 && BC != 2) yn += (int)ydim_;
                if (yn < 0 && BC == 2) yn = 0;
                zn = (int)site_[ii].getZ() - 1;
                if (zn < 0 && BC != 3) zn += (int)zdim_;
                if (zn < 0 && BC == 3) zn = 0;
                break;
            case 14: // Upwest neighbor
                xn = (int)site_[ii].getX() - 1;
                if (xn < 0 && BC != 1) xn += (int)xdim_;
                if (xn < 0 && BC == 1) xn = 0;
                yn = (int)site_[ii].getY();
                zn = (int)site_[ii].getZ() + 1;
                if (zn >= (int)zdim_ && BC != 3) zn -= (int)zdim_;
                if (zn >= (int)zdim_ && BC == 3) zn -= 1;
                break;
            case 15: // Upnorth neighbor
                xn = (int)site_[ii].getX();
                yn = (int)site_[ii].getY() + 1;
                if (yn >= (int)ydim_ && BC != 2) yn -= (int)ydim_;
                if (yn >= (int)ydim_ && BC == 2) yn -= 1;
                zn = (int)site_[ii].getZ() + 1;
                if (zn >= (int)zdim_ && BC != 3) zn -= (int)zdim_;
                if (zn >= (int)zdim_ && BC == 3) zn -= 1;
                break;
            case 16: // Upeast neighbor
                xn = (int)site_[ii].getX() + 1;
                if (xn >= (int)xdim_ && BC != 1) xn -= (int)xdim_;
                if (xn >= (int)xdim_ && BC == 1) xn -= 1;
                yn = (int)site_[ii].getY();
                zn = (int)site_[ii].getZ() + 1;
                if (zn >= (int)zdim_ && BC != 3) zn -= (int)zdim_;
                if (zn >= (int)zdim_ && BC == 3) zn -= 1;
                break;
            case 17: // Upsouth neighbor
                xn = (int)site_[ii].getX();
                yn = (int)site_[ii].getY() - 1;
                if (yn < 0 && BC != 2) yn += (int)ydim_;
                if (yn < 0 && BC == 2) yn = 0;
                zn = (int)site_[ii].getZ() + 1;
                if (zn >= (int)zdim_ && BC != 3) zn -= (int)zdim_;
                if (zn >= (int)zdim_ && BC == 3) zn -= 1;
                break;
        }

        idn = (unsigned int)(xn + (xdim_ * yn) +
                ((xdim_ * ydim_) * zn));
        site_[ii].setNb(j,&site_[idn]);
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

    volumefraction_.clear();
    volumefraction_.resize(numMicroPhases,0.0);
    initvolumefraction_.clear();
    initvolumefraction_.resize(numMicroPhases,0.0);

    if (verbose_) cout << "Calculating Volume Fractions now..." << endl;
    double vfrac,mfrac;
    try {
      if (site_.size() > 0) {
        for (ii = 0; ii < numMicroPhases; ii++) {
            vfrac = ((double)count_[ii])/((double)site_.size());
            setVolumefraction(ii,vfrac);
            setInitvolumefraction(ii,vfrac);
            if (getVolumefraction(ii) > 0.0 && verbose_) cout << "***Volume fraction["
                << ii << "] = " << volumefraction_[ii] << endl;
        }
      } else {
        msg = "Divide by zero error:  site_.size() = 0";
        throw FloatException("Lattice","Lattice",msg);
      }
      if (verbose_) cout << "...Done!" << endl;
    }
    catch (FloatException fex) {
      fex.printException();
      exit(1);
    }
    
    if (verbose_) cout << "Calculating weighted mean curvatures now..." << endl;
    for (ii = 0; ii < site_.size(); ii++) {
      site_[ii].calcWmc();
    }
    if (verbose_) cout << "...Done!" << endl;

}

Lattice::~Lattice ()
{
    interface_.clear();
    site_.clear();
    delete rg_;
}

void Lattice::addSite (const unsigned int x,
                       const unsigned int y,
                       const unsigned int z)
{
    string msg;
    
    try {
      if ( x >= xdim_ || x < 0 ) {
        throw EOBException("Lattice","addSite","site_",xdim_,x);
      } else if ( y >= ydim_ || y < 0 ) {
        throw EOBException("Lattice","addSite","site_",ydim_,y);
      } else if ( z >= zdim_ || z < 0 ) {
        throw EOBException("Lattice","addSite","site_",zdim_,z);
      }
    }
    catch (EOBException ex) {
      ex.printException();
      exit(1);
    }
    
    // Assume that the site is WATER by default.  This will break if
    // there is ever a sytem without water.
    
    site_.push_back(Site(x,y,z,xdim_,ydim_,zdim_,siteneighbors_,chemSys_));
}

void Lattice::findInterfaces (void)
{
    unsigned int i,kk;
    unsigned int k;
    vector<Site *> gsite,dsite;
    
    ///
    /// An interface must have at least one adjacent site that is water or void
    ///

    interface_.clear();
    for (i = 0; i < chemSys_->getNumMicroPhases(); i++) { 
      // if (verbose_) {
      //     cout << "  Database item " << i << ": ";
      //     cout.flush();
      // }
      if (i != WATERID && i != VOIDID) {
        // if (verbose_) cout << "Probing for interface... ";
        gsite.clear();
        dsite.clear();
        for (k = 0; k < site_.size(); k++) {
          if (site_[k].getWmc() > 0) {  // Is there some water nearby?
            if ((site_[k].getMicroPhaseId() == i)) {
              dsite.push_back(&site_[k]);
              site_[k].setDissolutionSite(i);
              for (kk = 0; kk < site_[k].nbSize(2); kk++) {
                if ((site_[k].nb(kk))->getMicroPhaseId() == WATERID) {
                  gsite.push_back(site_[k].nb(kk));
                  site_[k].nb(kk)->setGrowthSite(i);
                }
              }
              
              /// @note There is no reason to make the phase
              /// itself be a template for itself, nor water be
              /// a growth template for any phase, because we already
              /// tested for those above.
       
            } else if (chemSys_->isGrowthTemplate(i,site_[k].getMicroPhaseId())) {
              for (kk = 0; kk < site_[k].nbSize(1); kk++) {
                if ((site_[k].nb(kk))->getMicroPhaseId() == WATERID) {
                  gsite.push_back(site_[k].nb(kk));
                  site_[k].nb(kk)->setGrowthSite(i);
                }
              }
            }
          }
        }

        // if (verbose_) cout << " Done! " << dsite.size()
        //                    << " dissolution sites and " << gsite.size()
        //                    << " growth sites" << endl;

        if ((gsite.size() == 0) && (dsite.size() == 0)) {
          //   if (verbose_) cout << "Testing phase " << i << " for nucleation " << endl;

          ///
          /// We are dealing with a phase that may need to
          /// nucleate.  Identify the eligible nucleation
          //  sites for that phase
          ///

          double thresh = (0.5 / pow(resolution_,3.0));
          // if (verbose_) {
          //     cout << "Thresh = " << thresh << endl;
          //     cout.flush();
          // }
          double g = 0.0;
          for (k = 0; k < site_.size(); k++) {
            if ((site_[k].getMicroPhaseId() == WATERID)) {
              g = rg_->Ran3();
              if (g < thresh) {
                gsite.push_back(&site_[k]);
                site_[k].setGrowthSite(i);
              }
            }
          }
        }
        
        // if (verbose_) {
        //     cout << "Trying to add a water interface for phase " << i << "... ";
        //     cout.flush();
        // }
        interface_.push_back(Interface(chemSys_,rg_,gsite,dsite,i,verbose_,debug_));   
        // if (verbose_) {
        //     cout << "Done!" << endl;
        //     cout.flush();
        // }
      } else {
        // cout << "Trying to add a regular interface for phase " << i << "... ";
        // cout.flush();
        interface_.push_back(Interface(rg_,verbose_,debug_));   
      }
    }

    return;
}

int Lattice::growPhase (unsigned int phaseid,
                        int numtoadd)
{
    unsigned int pid;
    unsigned int i,j;
    double dwmcval;
    Site *ste,*stenb;

    try {
        if (numtoadd == 0) return 0;
        if (phaseid >= interface_.size()) {
            throw EOBException("Lattice","growPhase","interface_",
                               interface_.size(),phaseid);
        }
    }
    catch (EOBException ex) {
        ex.printException();
        exit(1);
    }

    vector<Isite> isite = interface_[phaseid].getGrowthSites();
    // if (verbose_) cout << "size of interface_[" << phaseid << "] is "
    //                   << isite.size() << endl; 
    int numleft,numchange = 0;

    ///
    /// We need to go through the interface list in
    /// normal order, which is the reverse order that
    /// we use for dissolution.
    ///

    numleft = numtoadd;
    if (verbose_) {
        cout << "-->Phase " << phaseid
            << " needs to grow at " << numtoadd
            << " sites" << endl;
    }

    while ((numleft > 0) && (isite.size() >= 1)) {
      for (i = 0; (numleft > 0) && (i < isite.size()); i++) {
        ste = &site_[isite[i].getId()];
        pid = ste->getMicroPhaseId();
     
        if (pid == WATERID) {
          removeGrowthSite(ste,phaseid);
          /*
          vector<unsigned int> plist = ste->getGrowthPhases();
          for (int ii = 0; ii < plist.size(); ii++) {
           removeGrowthSite(ste, plist[ii]);
          } 
          */

          ///
          /// Weighted mean curvature (wmc) is changed by the difference
          /// between the growing phase's porosity and the template's porosity.
          ///
          /// @todo Determine why the calculation works this way.
          ///

          dwmcval = chemSys_->getPorosity(phaseid)
                  - chemSys_->getPorosity(pid);
          setMicroPhaseId(ste,phaseid);
          ste->dWmc(dwmcval);

          ///
          /// Now that the site has been added, it is eligible for dissolution
          /// later on, so we add it to the list of dissolution sites.
          ///

          if (ste->getWmc() > 0.0) {
            addDissolutionSite(ste,phaseid);
          }

          ///
          /// Update the wmc of the neighboring sites.  This can be done
          /// because the wmc is originally calculated within a box around each
          /// site, so any time the id of a site within that box changes, it
          /// will change the wmc of the site at the box's center.
          ///

          for (j = 0; j < ste->nbSize(1); j++) {
            stenb = ste->nb(j);
            stenb->dWmc(dwmcval);
            if (stenb->getMicroPhaseId() == WATERID) {
              addGrowthSite(stenb,phaseid);
            } else if (stenb->getWmc() <= 0.0) {
              removeDissolutionSite(stenb,stenb->getMicroPhaseId());
            }
          }

          for (j = ste->nbSize(1); j < ste->nbSize(2); j++) {
            stenb = ste->nb(j);
            stenb->dWmc(dwmcval);
          }

          numleft--;
          numchange++;
        } else {
          removeGrowthSite(ste,phaseid);
        }
      }
  
      isite = interface_[phaseid].getGrowthSites();
    }
    
    return(numchange);
}


int Lattice::dissolvePhase (unsigned int phaseid,
                            int numtotake)
{
    unsigned int pid;
    double dwmcval;
    Site *ste,*stenb;

    if (numtotake == 0) return 0;
    try {
        if (phaseid >= interface_.size()) {
            throw EOBException("Lattice","dissolvePhase",
                           "interface_",interface_.size(),phaseid);
        }
    }
    catch (EOBException ex) {
        ex.printException();
        exit(1);
    }
    
    unsigned int i;
    vector<Isite> isite = interface_[phaseid].getDissolutionSites();
    if (verbose_) cout << "size of interface_[" << phaseid << "] is "
                       << isite.size() << endl; 
    try {
      for (i = 0; i < isite.size(); i++) {
        if (site_.at(isite[i].getId()).getMicroPhaseId() != phaseid) {
          cout << "Uh-oh... interface " << phaseid
               << " is corrupted with phase "
               << site_.at(isite[i].getId()).getMicroPhaseId() << endl;
          cout << "Offending site is " << isite[i].getId() << endl;
        }
      }
    }
    catch (out_of_range &oor) {
      EOBException ex("Lattice","dissolvePhase","site_",site_.size(),i);
      ex.printException();
      exit(1);
    }
    if (verbose_) cout << "-->Phase " << phaseid
                       << " needs to dissolve at " << numtotake
                       << " sites" << endl;
    
    int numleft = numtotake;
    int numchange = 0;
    while ((numleft > 0) && (isite.size() > 1)) {
  
      try {
        for (i = isite.size() - 1; (i > 0) && (numleft > 0); i--) {
          ste = &site_.at(isite[i].getId());
          pid = ste->getMicroPhaseId();
          removeDissolutionSite(ste,pid);
          setMicroPhaseId(ste,WATERID);
      
          ///
          /// Weighted mean curvature (wmc) is changed by the difference
          /// between the growing phase's porosity and the template's porosity.
          ///
          /// @todo Determine why the calculation works this way.
          ///

          dwmcval = chemSys_->getPorosity(WATERID)
                  - chemSys_->getPorosity(pid);
          ste->dWmc(dwmcval);
      
          unsigned int j;
          for (j = 0; j < ste->nbSize(1); j++) {
            stenb = ste->nb(j);
            stenb->dWmc(dwmcval);
            int nbpid = stenb->getMicroPhaseId();
            if ((nbpid != WATERID)
                && (nbpid != VOIDID)) {

              ///
              /// Now that the site has been dissolved, it is eligible for growth
              /// later on, so we add it to the list of growth sites for certain
              /// phases.
              ///

              addDissolutionSite(stenb,nbpid);
              addGrowthSite(ste,nbpid);
    
              vector<int> nbgrowthtemp;
              nbgrowthtemp.clear();
              nbgrowthtemp = chemSys_->getGrowthTemplate(nbpid);
              for (int ii = 0; ii < nbgrowthtemp.size(); ii++) {
                if (nbgrowthtemp[ii] != WATERID) {
                  addGrowthSite(ste,nbgrowthtemp[ii]);
                }
              }
            }
          }

          ///
          /// Update the wmc of the neighboring sites.  This can be done
          /// because the wmc is originally calculated within a box around each
          /// site, so any time the id of a site within that box changes, it
          /// will change the wmc of the site at the box's center.
          ///

          for (j = ste->nbSize(1); j < ste->nbSize(2); j++) {
            stenb = ste->nb(j);
            stenb->dWmc(dwmcval);
          }

          numleft--;
          numchange++;
        }
      }
      catch (out_of_range &oor) {
        EOBException ex("Lattice","dissolvePhase","site_",site_.size(),i);
        ex.printException();
        exit(1);
      }
    }
    
    isite = interface_[phaseid].getDissolutionSites();

    return(numchange);
}


void Lattice::addDissolutionSite (Site *ste,
                                  unsigned int pid)
{
    try {
        interface_.at(pid).addDissolutionSite(ste);
        vector<unsigned int> plist = ste->getGrowthPhases();
        for (unsigned int i = 0; i < plist.size(); i++) {
            interface_.at(plist[i]).removeGrowthSite(ste);
        }
        ste->setDissolutionSite(pid);
    }
    catch(out_of_range &oor) {
        EOBException ex("Lattice","addDissolutionSite","interface_",interface_.size(),pid);
        ex.printException();
        exit(1);
    }
    return;
}

void Lattice::addGrowthSite (Site *ste,
                             unsigned int pid)
{
    try {
        interface_.at(pid).addGrowthSite(ste);
        ste->setGrowthSite(pid);
    }
    catch(out_of_range &oor) {
        EOBException ex("Lattice","addGrowthSite","interface_",interface_.size(),pid);
        ex.printException();
        exit(1);
    }
}

void Lattice::removeDissolutionSite (Site *ste,
                                     unsigned int pid)
{
    try {
        interface_.at(pid).removeDissolutionSite(ste);
        ste->removeDissolutionSite(pid);
    }
    catch(out_of_range &oor) {
        EOBException ex("Lattice","removeDissolutionSite",
                        "interface_",interface_.size(),pid);
        ex.printException();
        exit(1);
    }
}

void Lattice::removeGrowthSite (Site *ste,
                                unsigned int pid)
{
    try {
        interface_.at(pid).removeGrowthSite(ste);
        ste->removeGrowthSite(pid);
    }
    catch(out_of_range &oor) {
        EOBException ex("Lattice","removeGrowthSite",
                        "interface_",interface_.size(),pid);
        ex.printException();
        exit(1);
    }
}

int Lattice::emptyPorosity (int numsites)
{
    unsigned int i,j;
    int numemptied = 0;
    int maxsearchsize = 3;
    unsigned int cntpore,cntmax;
    bool placed;

    if (numsites == 0) return (0);
    if (numsites < 0) {
        if (verbose_) {
            cout << "Going into Lattice::fillPorosity(" << -numsites << ")... " << endl;
            cout.flush();
        }
        numemptied = fillPorosity(-numsites);
        return(-numemptied);
    }
    
    ///
    /// Finding all potential VOID sites.
    ///
    /// @todo Consider removing some of the standard output, or setting a flag for it.
    ///

    if (verbose_) {
        cout << "Finding and sorting all potential void sites ... ";
        cout.flush();
    }
    list<Sitesize> distlist = findDomainSizeDistribution(WATERID,numsites,maxsearchsize,0);
    list<Sitesize>::iterator it;

    if (verbose_) {
        cout << "OK, found " << distlist.size()
             << " potential void sites." << endl;
        cout.flush();
    }
    
    ///
    /// We want to empty the sites with the largest pore count
    ///

    if (distlist.size() < numsites) {
        string msg = "Ran out of water in the system";
        EOBException ex("Lattice","emptysite",msg,distlist.size(),numsites);
        ex.printException();
        exit(1);
    }

    numemptied = 0;
    it = distlist.begin();
    int siteid;
    try {
        while (it != distlist.end()) {
            siteid = (*it).siteid;
            setMicroPhaseId(site_[siteid].getId(),VOIDID);
            numemptied++;
            it++;
        }
    }
    catch (EOBException eex) {
        eex.printException();
        cout << "Current value of numemptied = " << numemptied;
        cout.flush();
        exit(1);
    }

    return(numemptied);
}

int Lattice::fillPorosity (int numsites)
{
    unsigned int i,j;
    int numfilled = 0;
    int maxsearchsize = 10;
    unsigned int cntpore,cntmin;
    bool placed;

    if (verbose_) {
        cout << "In fillPorosity (" << numsites << ")" << endl;
        cout.flush();
    }
    if (numsites == 0) return (0);
    if (numsites < 0) {
        numfilled = emptyPorosity(-numsites);
        return(numfilled);
    }

    ///
    /// Finding all potential VOID sites.
    ///
    /// @todo Consider removing some of the standard output, or setting a flag for it.
    ///

    if (verbose_) cout << "Finding and sorting all potential water sites ... ";
    list<Sitesize> distlist = findDomainSizeDistribution(VOIDID,numsites,maxsearchsize,1);
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
            setMicroPhaseId(site_[siteid].getId(),WATERID);
            numfilled++;
            it++;
        }
    }
    catch (EOBException eex) {
        eex.printException();
        cout << "Current value of numfilled = " << numfilled;
        cout.flush();
        exit(1);
    }

    return(numfilled);
}

int Lattice::countBox (int boxsize,
                       unsigned int siteid)
{
    string msg;
    int boxhalf = boxsize / 2;
    int nfound = 0;
    int ix,iy,iz,hx,hy,hz;

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

        qxlo += checkBC(qxlo,xdim_);
        qxhi += checkBC(qxhi,xdim_);
        qylo += checkBC(qylo,ydim_);
        qyhi += checkBC(qyhi,ydim_);
        qzlo += checkBC(qzlo,zdim_);
        qzhi += checkBC(qzhi,zdim_);

        ///
        /// Count the number of water or void sites in the box
        ///

        for (iz = qzlo; iz <= qzhi; iz++) {
            hz = iz + checkBC(iz,zdim_);
            for (iy = qylo; iy <= qyhi; iy++) {
                hy = iy + checkBC(iy,ydim_);
                for (ix = qxlo; ix <= qxhi; ix++) {
                    hx = ix + checkBC(ix,xdim_);
                    if (site_.at(getIndex(hx,hy,hz)).getMicroPhaseId() == WATERID
                            || site_.at(getIndex(hx,hy,hz)).getMicroPhaseId() == VOIDID) {
                        nfound++;
                    }
                }
            }
        }
    }
    catch (out_of_range &oor) {
        EOBException ex("Lattice","countBox","site_",site_.size(),getIndex(hx,hy,hz));
        ex.printException();
        exit(1);
    }
    return(nfound);
}


void Lattice::setResolution (const double res)
{
    ///
    /// Make sure that resolution is a valid value
    ///

    try {
        string msg;
        if (res <= 0.001) {
            cout << endl;
            msg = "Lattice resolution <= 0.001";
            throw DataException("Lattice","setResolution",msg);
        }

        if (verbose_) {
            cout << "Changing lattice resolution from ";
            cout << resolution_ << " to " << res << endl;
            cout.flush();
        }
        resolution_ = res;
    }
    catch (DataException dex) {
        dex.printException();
        exit(1);
    }
    return;
}

vector<unsigned int> Lattice::getNeighborhood (const unsigned int sitenum,
                                               const int size)
{
    int xp,yp,zp;
    double dist;

    vector<unsigned int> nh;
    nh.clear();

    unsigned int xc = site_[sitenum].getX();
    unsigned int yc = site_[sitenum].getY();
    unsigned int zc = site_[sitenum].getZ();

    if (size == 0) {
        nh.push_back(getIndex(xc,yc,zc));
        return nh;
    }
    
    for (int k = -size; k <= size; k++) {
        zp = zc + k;
        for (int j = -size; j <= size; j++) {
            yp = yc + j;
            for (int i = -size; i <= size; i++) {
                xp = xc + i;
                dist = (double)((xc - xp) * (xc - xp));
                dist += (double)((yc - yp) * (yc - yp));
                dist += (double)((zc - zp) * (zc - zp));
                /*
                if ((sqrt(dist) - 0.5) <= (double)size) {
                    nh.push_back(getIndex(xp,yp,zp));
                }
                */
                nh.push_back(getIndex(xp,yp,zp));
            }
        } 
    }

    return nh;
}

unsigned int Lattice::getIndex (int ix,
                                int iy,
                                int iz) const
{
   if (ix < 0) {
       if (BC != 1) {
           ix += (int)xdim_;
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

   return (unsigned int)(ix + (xdim_ * iy)
           + ((xdim_ * ydim_) * iz)); 
}

void Lattice::changeMicrostructure (double time,
                                    const int simtype,
                                    bool isFirst,
                                    bool &capWater)
{
    unsigned int i,ii;
    int numadded, numadded_actual;
    unsigned int tpid;
    int cursites,newsites,tnetsites;
    int wcursites,wnewsites;
    double td,tvol,tmass;
    vector<double> vol_next,vfrac_next;
    vector<int> netsites;
    vector<unsigned int> pid;
    vector<string> phasenames;

    extern string CSHMicroName;
    extern string MonocarbMicroName;
    extern string MonosulfMicroName;
    extern string HydrotalcMicroName;
    extern string AFTMicroName;

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
      cout << "can't open file " << fexpansioncoor << ", so exit program." << endl;
      exit(1);
    } else {
      while (!in1.eof()) {
        int index;
        vector<int> coordin;
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

    /// JWB: 2020 Dec 22  Do manual adjustment of certain microstructure
    /// phase voxels in a cement paste, accounting for subvoxel porosity, etc.
    /// This is VERY rough right now.
    /// Solid microstructure phase adjustments are always traded with
    /// adjustment to saturated capillary pore volume to keep the total
    /// volume constant.
    //
    /// If this works, then ...
    ///
    /// @todo make a new attribute for microstructure phases which is
    ///       a multiplicative factor on the phase volume fraction.
    ///
   
    try {
        adjustMicrostructureVolumes(phasenames,vol_next);

        adjustMicrostructureVolFracs(phasenames,vol_next,vfrac_next);
    }
    catch (DataException dex) {
        throw dex;
    }
    catch (EOBException ex) {
        throw ex;
    }

    ///
    /// Calculate number of sites of each phase in next state
    ///

    if (verbose_) {
        cout << "Calculating volume of each phase to be added..." << endl;
        try {
          for (i = 0; i < vfrac_next.size(); i++) {
              cout << "****Volume fraction[" << phasenames.at(i)
                   << "] in next state should be = " << vfrac_next.at(i);
              cout << ", or "
                   << (int)((double)(numsites_ * vfrac_next.at(i)))
                   << " sites" << endl;
          }
        }
        catch (out_of_range &oor) {
          throw EOBException("Lattice","changeMicrostructure","phasenames",
                             phasenames.size(),i);
        }
    }

    ///
    /// The next block is executed only if there will eventually be some
    /// sulfate attack during this simulation.
    ///

    if (simtype == SULFATE_ATTACK) {

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
        netsites.resize(numMicroPhases,0);
        pid.clear();
        pid.resize(numMicroPhases,0);

        vector<int> growing;
        vector<vector<int> > shrinking;
        vector<vector<double> > volratios;
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
            shrinking.resize(growing.size(),idummy);
            volratios.resize(growing.size(),ddummy);
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
                    if (i != WATERID && i != VOIDID) tnetsites = (newsites - cursites);
                    netsites.at(i) = tnetsites;
                    pid.at(i) = i;
                    if (i == growing[ii] && isFirst) {
                      netsites.at(i) = 0;
                      count_.at(i) = newsites;
                    }        
                    if (netsites.at(i) != 0 && verbose_) cout << "****netsites["
                                                              << phasenames.at(i)
                                                              << "] in this state = "
                                                              << netsites.at(i)
                                                              << endl;
                  }
            
                  ///
                  /// Next block gets executed only if we are now simulating
                  /// sulfate attack.
                  ///

                  if (time_ >= sattack_time_) {
                      if (verbose_) {
                          cout << "start to crystal-pressure transform at time_ = " 
                               << time_ << endl;
                      }

                      ///
                      /// The relevant stress-free molar volume ratios for sulfate
                      /// attack phase transformations.
                      ///

                      vector<int> numchanged;
                      numchanged.clear();
                      int growid = growing[ii];
                      for (int iii = 0; iii < shrinking[ii].size(); ++iii) {
                          numchanged.resize(2,0);
                          int shrinkid = shrinking[ii][iii];
                          double volrat = volratios[ii][iii];
                          if ((netsites.at(shrinkid) < 0) && (netsites.at(growid) > 0)) {
                              numchanged = transform(shrinkid,
                                                     netsites.at(shrinkid),
                                                     growid,
                                                     netsites.at(growid),
                                                     volrat);

                              netsites.at(shrinkid) += numchanged[0];
                              netsites.at(growid) -= numchanged[1];
                              if (verbose_) {
                                  cout << "netsites.at(" << shrinkid << ") is: "
                                       << netsites.at(shrinkid) << endl;
                                  cout << "netsites.at(" << growid << ") is: "
                                       << netsites.at(growid) << endl;
                              }
                          }

                      }
                  }
            }
        }

        catch (out_of_range &oor) {
            throw EOBException("Lattice","changeMicrostructure",
                    "phasenames or count_ or pid or netsites or vfrac_next",
                    phasenames.size(),i);
        }

    } else {
        
        ///
        /// Sulfate attack will NEVER be done during this simulation.
        ///  Normalize to get volume fractions and compute number
        ///  of sites of each phase needed.
        ///

        netsites.clear();
        netsites.resize(chemSys_->getNumMicroPhases(),0);
        pid.clear();
        pid.resize(chemSys_->getNumMicroPhases(),0);

        try {
            for (i = FIRST_SOLID; i < vfrac_next.size(); i++) {
                cursites = (int)(count_.at(i) + 0.5);
                newsites = (int)((numsites_
                             * vfrac_next.at(i)) + 0.5);
                tnetsites = 0;
                tnetsites = (newsites - cursites);
                netsites.at(i) = tnetsites;
                pid.at(i) = i;
                if (netsites.at(i) != 0 && verbose_) cout << "***netsites["
                                                          << phasenames.at(i)
                                                          << "] in this state = "
                                                          << netsites.at(i)
                                                          << "; cursites = " << cursites
                                                          << " and newsites = " << newsites
                                                          << endl;
            }
        }
        catch (out_of_range &oor) {
            throw EOBException("Lattice","changeMicrostructure",
                    "phasenames or count_ or pid or netsites or vfrac_next",
                     phasenames.size(),i);
        }
    }
    
    if (verbose_) {
        cout << "Sorting non-pore sites... netsites[VOIDID] = netsites["
             << VOIDID << "] = " << netsites[VOIDID] << endl;
        cout.flush();
    }

    ///
    /// Sort netsites in ascending order, except we will handle
    /// void space differently.
    ///

    try {
        for (i = FIRST_SOLID; i < netsites.size(); i++) {
            for (ii = i; ii < netsites.size(); ii++) {
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
    }
    catch (out_of_range &oor) {
        throw EOBException("Lattice","changeMicrostructure","pid", pid.size(),ii);
    }

    Interface ifc;
    vector<Isite> gs,ds;
    if (verbose_) cout << "Getting change vectors for non-pore phases... netsites[VOIDID] = "
                       << netsites[VOIDID] << endl;

    ///
    /// The next loop starts at FIRST_SOLID because we exclude void and water phases
    ///
    /// @todo Consider making the starting index more general
    ///

    for (i = FIRST_SOLID; i < interface_.size(); i++) {
        ifc = interface_[i];
        gs = ifc.getGrowthSites();
        ds = ifc.getDissolutionSites();
    }

    if (verbose_) cout << "Switching phase id values..." << endl;
    try {
        for (i = FIRST_SOLID; i < netsites.size(); i++) {
            numadded = 0;
            numadded_actual = 0;
	  
	        if (netsites[i] < 0) {
                if (verbose_) {
                    cout << "Going into dissolve_phase now... pid = "
                         << pid.at(i) << endl;
                    cout.flush();
                }
                numadded = dissolvePhase(pid.at(i),-netsites[i]);
                numadded_actual += numadded;
                if (verbose_) {
                    cout << "...Done with dissolve_phase for phase "
                         << pid.at(i) << "!" << endl;
                    cout.flush();
                }
            } else if (netsites[i] > 0) {
                if (verbose_) {
                    cout << "Going into grow_phase now... pid = "
                         << pid.at(i) << endl;
            
                    cout.flush();
                }
 
                numadded = growPhase(pid.at(i),netsites[i]);

                numadded_actual += numadded;
                int diff = netsites[i] - numadded;
                while (diff > 0) {
                    gs = interface_[pid.at(i)].getGrowthSites();
                    ds = interface_[pid.at(i)].getDissolutionSites();
                    if (verbose_) {
                        cout << "gs.size() = " << gs.size() << endl;
                        cout << "ds.size() = " << ds.size() << endl;
                    }
                    if (gs.size() == 0) {
                        if (verbose_) cout << "Phase " << pid.at(i)
                                           << " needs to be nucleated. " << endl;

                        int nuclei = diff;
	                    double thresh = (nuclei / (volumefraction_[WATERID] * site_.size()));
	                    if (verbose_) cout << "Thresh = " << thresh << endl;
                        if (thresh < 1) {

                            // Compose a shuffled list of all the water sites
                        
                            vector<int> watersites;
                            watersites.clear();
                            for (int k = 0; k < site_.size(); ++k) {
                                if (site_[k].getMicroPhaseId() == WATERID) watersites.push_back(k);
                            }
                            if (verbose_) {
                                cout << "Found " << watersites.size() << " water sites" << endl;
                                cout.flush();
                            }
                            int numshuffles = watersites.size() / 2;
                            if (verbose_) {
                                cout << "Preparing to perform a shuffle" << endl;
                                cout << "  First ten water sites prior to shuffle:" << endl;
                                for (int k = 0; k < 10; ++k) {
                                    cout << "      " << k << ": " << watersites[k] << endl;
                                    cout.flush();
                                }
                            }

                            rg_->shuffle(watersites,1);

                            if (verbose_) {
                                cout << "Done with shuffle" << endl;
                                cout << "  First ten water sites after to shuffle:" << endl;
                                for (int k = 0; k < 10; ++k) {
                                    cout << "      " << k << ": " << watersites[k] << endl;
                                    cout.flush();
                                }
                            }

                            // Now the list is thoroughly shuffled, we just pick
                            // the first nuclei elements as new growth sites
                        
                            for (int k = 0; k < nuclei; ++k) {
                                addGrowthSite(&site_[watersites[k]],pid.at(i));
    	                    }
                        } else {
                            if (verbose_) {
                                cout << "There is no room to grow, so exit the program."
                                     << endl;
                            }
                            throw MicrostructureException("Lattice","changeMicrostructure","no room to grow phase");
                        }
	                }
                    numadded = growPhase(pid.at(i),diff);
    	            numadded_actual += numadded;
                    diff = diff - numadded;
	                if (verbose_) cout << "diff = " << diff << endl;
	            }
                if (verbose_) {
    	            cout << "...Done with grow_phase for phase "
                         << pid.at(i) << "!" << endl;
                    cout.flush();
                }
            }
        
            if (numadded_actual*numadded_actual != netsites[i]*netsites[i]) {
                if (warning_) {
                    cout << "WARNING: Needed to switch on "
                         << netsites[i] << " of phase " << pid.at(i) << endl;
                    cout << "         But actually did "
                         << numadded << " switches" << endl;
                    cout.flush();
                }
            }
        }
    }

    catch (out_of_range &oor) {
        throw EOBException("Lattice","changeMicrostructure","pid",
                            pid.size(),ii);
    }
    
    if (verbose_) {
        cout << "time_ in the Lattice is: " << time_ << endl;
        cout << "Now emptying the necessary pore volume:  " << endl;
    }

    cursites = count_.at(VOIDID);
    newsites = (int)((numsites_
                 * vfrac_next.at(VOIDID)) + 0.5);
    wcursites = count_.at(WATERID);
    wnewsites = wcursites - (newsites-cursites);
    int numempty = emptyPorosity(newsites - cursites);
    if (verbose_) {
        cout << "***netsites["
             << phasenames.at(VOIDID)
             << "] in this state = "
             << (newsites - cursites)
             << "; cursites = " << cursites
             << " and newsites = " << newsites
             << endl;
        cout << "***netsites["
             << phasenames.at(WATERID)
             << "] in this state = "
             << (wnewsites - wcursites)
             << "; cursites = " << wcursites
             << " and newsites = " << wnewsites
             << endl << endl;

        // When creating void from water, we should
        // update the target volume fraction of water even though
        // it is not used in any further calculations at this point
 
        cout << "Target volume fraction of water IS " << vfrac_next.at(WATERID) << endl;
        // vfrac_next.at(WATERID) -= ((double)(newsites - cursites)/(double)(numsites_));
        // cout << "But WILL BE " << vfrac_next.at(WATERID) << " after creating void space" << endl;

        cout << "Number actually emptied was:  " << numempty << endl;

        /// 
        /// Report on target and actual mass fractions
        ///

        cout << "*******************************" << endl;
    }
    try {
        for (i = 0; i < vfrac_next.size(); i++) {
            volumefraction_.at(i) = ((double)(count_.at(i)))/((double)(site_.size()));
            if (verbose_) {
                cout << "Phase " << i << " Target volume fraction was "
                     << vfrac_next[i] << " and actual is "
                     << volumefraction_.at(i) << endl;
            }
        }
    }

    catch (out_of_range &oor) {
        throw EOBException("Lattice","changeMicrostructure",
                           "volumefraction_ or count_",volumefraction_.size(),i);
    }

    if (verbose_) cout << "*******************************" << endl;
    
    for (i = FIRST_SOLID; i < interface_.size(); i++) {
        interface_[i].sortGrowthSites(site_,i);
        interface_[i].sortDissolutionSites(site_,i);
    }

    /*
    //update expansion and expansion_coordin files
    cout << "in Lattice, expansion_.size() = " << expansion_.size() << endl;
    ofstream out(fexpansion.c_str());
    for (map<int, vector<double> >::iterator it = expansion_.begin();
      it != expansion_.end(); it++) {
      int index = it -> first;
      vector<double> expval = it -> second;
      out << index << " " << expval[0] << " " << expval[1] << " " << expval[2] 
          << endl;
    }
    out.close();

    ofstream out1(fexpansioncoor.c_str());    
    for (map<int, vector<int> >::iterator it = expansion_coordin_.begin();
      it != expansion_coordin_.end(); it++) {
      int index = it -> first;
      vector<int> coor = it -> second;
      out1 << index << " " << coor[0] << " " << coor[1] << " " << coor[2]
           << endl;
    }
    out1.close();
    */

    ///  This is a local variable and the value is never used.
    ///
    ///  @todo Why not eliminate this line completely?
    ///

    double surfa = getSurfaceArea(chemSys_->getMicroPhaseId(CSHMicroName));
    
    if (volumefraction_[WATERID] <= 0.0) capWater = false;

    return;
}

void Lattice::adjustMicrostructureVolumes (vector<string> name,
                                           vector<double> &vol)
{
    int i = 0;
    double cap_watervolume = 0.0;
    double cap_spacevolume = 0.0;
    double voidvolume = 0.0;
    double subvoxel_watervolume = 0.0;
    double subvoxel_volume = 0.0;
    double current_watervolume = 0.0;

    extern string CSHMicroName;

    try {
        
        if (verbose_) {
            cout << endl;
            cout << "@@@@@@@@ ADJUSTING MICROSTRUCTURE PHASE VOLUMES @@@@@@@@" << endl;
            cout << endl;
        }

        // Find the total system volume according to GEMS
        
        double GEM_thinks_totmicvol = chemSys_->getMicroVolume();
        double GEM_thinks_totmicinitvol = chemSys_->getMicroInitVolume();

        // Find the total SOLIDS volume
        
        // vol is the GEMS-based volume not counting
        // subvoxel pores.  The volume fraction of these
        // subvoxel pores has been input by the user as
        // 'porosity'.  We need to alter the apparent volume
        // of the GEMS-based volume to include these pores
        // in the 'apparent volume' observed at the voxel scale
       
        double apparent_solidvol = 0.0;

        // Following is a dummy variable for holding calculations
        double apparent_volume = 0.0;
        // Following will hold the subvoxel porosity of a phase
        double phi = 0.0;
        double upper_cutoff_porosity = 0.8;
        subvoxel_volume = 0.0;

        for (i = 0; i < vol.size(); ++i) {
            if (i != WATERID && i != VOIDID) {
                phi = chemSys_->getPorosity(i);
                if (phi > upper_cutoff_porosity) phi = upper_cutoff_porosity;
                apparent_volume = vol.at(i) * (1.0 + (phi/(1.0 - phi)));
                subvoxel_volume += (vol.at(i) * phi);
                vol.at(i) = apparent_volume;
                apparent_solidvol += apparent_volume;

                if (verbose_) {
                    cout << "    Phase name " << name.at(i) << ": volume = " << vol.at(i) << endl;
                }
            }
        }

        if (verbose_) {
            cout << "Initial microstructure volume was " << GEM_thinks_totmicinitvol << endl;
            cout << "Current microstructure volume is " << GEM_thinks_totmicvol << endl;
        }

        if (apparent_solidvol <= 0.0) {
            throw DataException("Lattice",
                                "adjustMicrostructureVolumes",
                                "totvolume is NOT positive");
        }

        // The capillary SPACE volume is all the volume in the initial microstructure
        // that is not occupied by solids or subvoxel porosity

        cap_spacevolume = (GEM_thinks_totmicinitvol - apparent_solidvol);
        double spaceforwater = cap_spacevolume + subvoxel_volume;

        if (!(chemSys_->isSaturated())) {   // System is sealed

            // Step 1: Total volume of void space in RVE is equal to the initial
            // volume according to GEMS minus the current total volume according
            // to GEMS.  At the moment we do not divide pores into capillary
            // and subvoxel; both of them make up the new_voidvolume
           
            double new_voidvolume = GEM_thinks_totmicinitvol - GEM_thinks_totmicvol;

            // Step 2: Calculate total space occupied by capillary pores (those
            // with a diameter smaller than a single voxel). This capillary SPACE volume
            // is the difference between the *initial* microstructure volume
            // according to GEMS (solids + water) and the *current* apparent volume
            // of solid phases (all solids plus any subvoxel porosity)

            cap_spacevolume = (GEM_thinks_totmicinitvol - apparent_solidvol);

            // Following is the current free liquid water volume according
            // to GEMS, not the volume of capillary pore water.  In
            // reality it will be divided between capillary pores
            // (>= one voxel volume) and subvoxel pores (< one voxel volume)
                                                   
            current_watervolume = vol.at(WATERID);
                                                   
            double current_voidvolume = vol.at(VOIDID);

            // Step 3: Calculate how much of the capillary pore volume (pores
            // with diameter greater or equal to a voxel) is occupied by liquid
            // water.  Free liquid water goes into subvoxel pores first.  If any
            // liquid water is left over after filling the subvoxel pores,
            // then it goes into capillary pores (>= voxel scale)
        
            if (current_watervolume >= subvoxel_volume) {
                cap_watervolume = current_watervolume - subvoxel_volume;
                subvoxel_watervolume = subvoxel_volume;
            } else {
                subvoxel_watervolume = current_watervolume;
                cap_watervolume = 0.0;
            }

            if (verbose_) {
                cout << "System is SEALED:" << endl;
                cout << "Volume change should be made up by voids with volume "
                    << new_voidvolume << endl;
                cout << "Current GEM water volume = " << current_watervolume
                    << endl;
                cout << "   This is divided between 1. " << subvoxel_watervolume
                     << " in subvoxel pores (< one voxel volume)" << endl;
                cout << "                           2. " << cap_watervolume
                                                         << " in capillary pores (>= one voxel volume)"
                                                         << endl;
                cout << "Current GEM void volume = " << current_voidvolume << endl;
            }

            // Step 4: Reassign the volume of void and water so that the water now
            // refers only to water in capillary pores that are larger in diameter
            // than one voxel. We can always recover the total water volume (capillary
            // plus subvoxel) from GEMS if we need to.
           
            if (new_voidvolume > 0) {
                // The volume of VOIDID is still not differentiated between voids in
                // capillary pores and those in subvoxel pores.
                
                vol.at(VOIDID) = new_voidvolume;
                
                // We now change the meaning of vol.at(WATERID).  Up to this point
                // it has been the volume of liquid water in the system as determined
                // by GEMS.  Now, we are instead going to store the volume of
                // liquid water in capillary pores only.
                
                vol.at(WATERID) = cap_watervolume - vol.at(VOIDID);
            } else {
                vol.at(VOIDID) = 0.0;
                vol.at(WATERID) = cap_watervolume;
            }

        } else { // The system is saturated with water, so there will never be voids

            // @note Do we eventually want to check for capillary pore percolation as
            // a requirement for saturated conditions?
           
            vol.at(VOIDID) = 0.0;
            vol.at(WATERID) = spaceforwater;
        }

        if (verbose_) {
            cout << endl;
            cout << "@@@@@@@@ FINISHED ADJUSTING MICROSTRUCTURE PHASE VOLUMES @@@@@@@@" << endl;
            cout << endl;
        }

        ///
        /// End of manual adjustment  
        ///

    }

    catch (DataException dex) {
        throw dex;
    }

    catch (out_of_range &oor) {
      throw EOBException("Lattice","adjustMicrostructureVolumes","name",name.size(),i);
    }

    return;
}

void Lattice::adjustMicrostructureVolFracs (vector<string> &names,
                                            const vector<double> vol,
                                            vector<double> &vfrac)
{
    int i = 0;
    double totmicvolume = 0.0;

    try {

        vfrac.clear();
        vfrac.resize(vol.size(),0.0);

        if (verbose_) {
            cout << endl;
            cout << "@@@@@@@@ ADJUSTING MICROSTRUCTURE PHASE VOLUME FRACTIONS @@@@@@@@" << endl;
            cout << endl;
        }

        // Find total adjusted microstructure volume, including possibly voids
       
        for (i = 0; i < vol.size(); ++i) {
            totmicvolume += vol.at(i);
            if (verbose_) {
                cout << " Volume(" << names.at(i) << ") = " << vol.at(i) << endl;
            }
        }

        if (verbose_) cout << "Calculated total microstructure volume is " << totmicvolume << endl;

        // Calculate volume fractions based on total microstructure volume
        
        for (i = 0; i < vol.size(); ++i) {
            vfrac.at(i) = vol.at(i) / totmicvolume;
            if (verbose_) {
                cout << "Volume fraction[" << names.at(i) << "] should be "
                     << vfrac.at(i) << ", count NOW is "
                     << (double)(count_.at(i))/(double)(numsites_) << endl;
            }
        }
        
        if (verbose_) {
            cout << endl;
            cout << "@@@@@@@@ FINISHED ADJUSTING MICROSTRUCTURE PHASE VOLUME FRACTIONS @@@@@@@@" << endl;
            cout << endl;
        }
    }

    catch (DataException dex) {
        dex.printException();
        exit(1);
    }

    catch (out_of_range &oor) {
      EOBException ex("Lattice","adjustMicrostructureVolFracs","vol",vol.size(),i);
      ex.printException();
      exit(1);
    }
}

void Lattice::writePoreSizeDistribution (double curtime,
                                         const int simtype,
                                         const string &root)
{
    // First get the volumes of each phase
    vector<double> vol = chemSys_->getMicroPhaseVolume();

    // Following is a dummy variable for holding calculations
    vector<double> subpore_volume;
    subpore_volume.resize(vol.size(),0.0);

    // Following will hold the subvoxel porosity of a phase
    double phi = 0.0;
    double upper_cutoff_porosity = 0.8;

    // Get the combined pore volume as a function of diameter 
    for (int i = 0; i < vol.size(); ++i) {
        if (i != WATERID && i != VOIDID) {
            phi = chemSys_->getPorosity(i);
            if (phi > upper_cutoff_porosity) phi = upper_cutoff_porosity;
            subpore_volume[i] = vol.at(i) * phi;
        }
    }

    // Use the normalized pore size distributions to determine the
    // volume at each pore size
    
    vector<vector<struct PoreSizeVolume> > porevolume;
    porevolume = chemSys_->getPoreSizeDistribution();
    if (verbose_) {
        cout << "Number of pore distributions available = " << porevolume.size() << endl;
        for (int i = 0; i < porevolume.size(); ++i) {
            cout << "   Pore distribution " << i << ": number of diameters = "
                 << porevolume[i].size() << endl;
        }
        cout.flush();
    }

    for (int i = 0; i < vol.size(); ++i) {
        if (i != WATERID && i != VOIDID) {
            for (int j = 0; j < porevolume[i].size(); ++j) {
                porevolume[i][j].volume =
                    porevolume[i][j].volfrac * subpore_volume[i];
            }
        }
    }

    // At this point we have the (absolute) volume of all subvoxel pores
    // as a function of their diameters. Bin them in increments of
    // 1 nanometer

    // Create and initialize the binned volume distribution
    vector<vector<struct PoreSizeVolume> > binnedporevolume;
    vector<struct PoreSizeVolume> zpsvec;
    zpsvec.clear();
    binnedporevolume.resize(vol.size(),zpsvec);
    // Done initializing the binned pore volume distribution

    double maxsize = 1.0;
    double maxmaxsize = 1.0;
    double cumulative_volume = 0.0;
    double cumulative_volfrac = 0.0;
    struct PoreSizeVolume dpsv;

    for (int i = 0; i < vol.size(); ++i) {
        if (i != WATERID && i != VOIDID) {
            maxsize = 1.0;
            int j = 0;
            while (j < porevolume[i].size()) {
                while ((porevolume[i][j].diam <= maxsize)
                    && (j < porevolume[i].size())) {
                    if (verbose_) {
                        cout << "Made it 01a, j = " << j;
                        cout << ", max j = " << porevolume[i].size() - 1;
                        cout.flush();
                        cout << ", pore diam " << i << "," << j
                             << " = " << porevolume[i][j].diam;
                        cout.flush();
                        cout << ", vol " << porevolume[i][j].volume;
                        cout.flush();
                        cout << ", volfrac " << porevolume[i][j].volfrac;
                        cout.flush();
                        cout << "maxsize = " << maxsize << endl;
                        cout.flush();
                    }

                    cumulative_volume += porevolume[i][j].volume;
                    cumulative_volfrac += porevolume[i][j].volfrac;
                    j++;
                }
                dpsv.diam = maxsize;
                dpsv.volume = cumulative_volume;
                dpsv.volfrac = cumulative_volfrac;
                binnedporevolume[i].push_back(dpsv);
                dpsv.diam = dpsv.volume = dpsv.volfrac = 0.0;
                cumulative_volume = cumulative_volfrac = 0.0;
                maxsize += 1.0;
            }
            if (maxsize > maxmaxsize) maxmaxsize = maxsize;
        }
    }

    // Now all subvoxel pores have been binned in 1-nm bins.
    // Combine them into a single distribution for all phases
    vector<struct PoreSizeVolume> masterporevolume;
    struct PoreSizeVolume dumps;
    dumps.diam = 0.0;
    dumps.volfrac = 0.0;
    dumps.volume = 0.0;
    masterporevolume.resize(int(maxmaxsize)-1,dumps);

    for (int i = 0; i < binnedporevolume.size(); ++i) {
        cumulative_volume = 0.0;
        cumulative_volfrac = 0.0;
        for (int j = 0; j < binnedporevolume[i].size(); ++j) {
            masterporevolume[j].diam = binnedporevolume[i][j].diam;
            masterporevolume[j].volume += binnedporevolume[i][j].volume;
            masterporevolume[j].volfrac = 0.0;
        }
    }

    // Fill up pores with available liquid water volume,
    // starting with the smallest pores.  This ASSUMES
    // that the water-solid surface energy is the same
    // for all pore surfaces.
    //
    /// @todo Revisit this assumption later

    double watervol = chemSys_->getMicroPhaseVolume(WATERID);
    int k = 0;
    while ((watervol - masterporevolume[k].volume >= 0.0)
            && (k < masterporevolume.size())) {
        watervol -= masterporevolume[k].volume;
        masterporevolume[k].volfrac = 1.0;
        k++;
    }

    // We have either filled up all the sub-voxel pores
    // OR we have reached the size at which the pores are
    // no longer fully saturated

    if (k == masterporevolume.size()) {
        // All sub-voxel pores are filled
    } else {
        masterporevolume[k].volfrac = watervol / masterporevolume[k].volume;
        watervol = 0.0;
    }


    string ofileName(root);
    ostringstream ostr1,ostr2;
    // Add the time in minutes
    ostr1 << setfill('0') << setw(6) << (int)((curtime * 24.0 * 60.0) + 0.5);
    ostr2 << setprecision(3) << temperature_;
    string timestr(ostr1.str());
    string tempstr(ostr2.str());
    ofileName = ofileName + "_PoreSizeDistribution." + timestr + "." + tempstr + ".csv";
    if (verbose_) {
        cout << "    In Lattice::writePoreSizeDistribution, curtime = " << curtime
             << ", timestr = " << timestr << endl;
        cout.flush();
    }

    ofstream out(ofileName.c_str());
    try {
        if (!out.is_open()) {
            throw FileException("Lattice",
                                "writePoreSizeDistribution",
                                ofileName,"Could not open");
        }
    }
    catch (FileException fex) {
        fex.printException();
        exit(1);
    }

    // Write the header
    out << "Diameter (nm),Volume,Fraction Saturated" << endl;
    for (int i = 0; i < masterporevolume.size(); ++i) {
        out << masterporevolume[i].diam << ","
            << masterporevolume[i].volume << ","
            << masterporevolume[i].volfrac << endl;
    }
    out << ">" << masterporevolume[masterporevolume.size()-1].diam << ","
        << watervol << ",??" << endl;

    out.close();

    return;
}

void Lattice::writeLattice (double curtime, const int simtype, const string &root)
{
    unsigned int i,j,k;
    string ofileName(root);
    ostringstream ostr1,ostr2;
    ostr1 << setfill('0') << setw(6) << (int)((curtime * 24.0 * 60.0) + 0.5);	// minutes
    ostr2 << setprecision(3) << temperature_;
    string timestr(ostr1.str());
    string tempstr(ostr2.str());
    ofileName = ofileName + "." + timestr + "." + tempstr + ".img";
    if (verbose_) {
        cout << "    In Lattice::writeLattice, curtime = " << curtime
             << ", timestr = " << timestr << endl;
        cout.flush();
    }

    ofstream out(ofileName.c_str());
    try {
        if (!out.is_open()) {
            throw FileException("Lattice","writeLattice",ofileName,"Could not open");
        }
    }
    catch (FileException fex) {
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
               int index = getIndex(i,j,k);
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
                throw FileException("Lattice","writeLattice",ofileName,"Could not open");
            }
        }
        catch (FileException fex) {
            fex.printException();
            exit(1);
        }

        // Write image header information first

        out1 << VERSIONSTRING << " " << version_ << endl;
        out1 << XSIZESTRING << " " << xdim_ << endl;
        out1 << YSIZESTRING << " " << ydim_ << endl;
        out1 << ZSIZESTRING << " " << zdim_ << endl;
        out1 << IMGRESSTRING << " " << resolution_ << endl;
 
        int DAMAGEID = 100;     // Some large number that cannot represent any other phase
        for (k = 0; k < zdim_; k++) {
            for (j = 0; j < ydim_; j++) {
                for (i = 0; i < xdim_; i++) {
                   int index = getIndex(i,j,k);
                   if (site_[index].IsDamage()) {
                     out1 << DAMAGEID << endl;
                   } else {
                     out1 << site_[index].getMicroPhaseId() << endl;
                   }
                }
            }
        }

        out1.close();
    } //The above block is implemented only if we are dealing with sulfate attack
}

void Lattice::writeDamageLattice (double curtime, const string &root)
{
    unsigned int i,j,k;
    string ofileName(root);
    ostringstream ostr1,ostr2;
    ostr1 << (int)((curtime * 100.0) + 0.5);	// hundredths of a day
    ostr2 << setprecision(3) << temperature_;
    string timestr(ostr1.str());
    string tempstr(ostr2.str());
    ofileName = ofileName + "." + timestr + "." + tempstr + ".img";

    ofstream out(ofileName.c_str());
    try {
        if (!out.is_open()) {
            throw FileException("Lattice","writeLattice",ofileName,"Could not open");
        }
    }
    catch (FileException fex) {
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
               int index = getIndex(i,j,k);
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

void Lattice::writeLatticePNG (double curtime, const int simtype, const string &root)
{
    unsigned int i,j,k;
    string ofileName(root);
    string ofbasename(root);
    string ofpngname(root);
    string ofpngbasename(root);

    vector<double> dumvec;
    vector<unsigned int> idumvec;
    vector<vector<unsigned int> > image;
    vector<vector<double> > dshade;
    dumvec.resize(ydim_,0.0);
    idumvec.resize(ydim_,0);
    dshade.resize(xdim_,dumvec);
    image.resize(xdim_,idumvec);
    bool done;

    ///
    /// Construct the name of the output file
    ///

    ostringstream ostr1,ostr2;
    ostr1 << setfill('0') << setw(6) << (int)((curtime * 24.0 * 60.0) + 0.5);	// minutes
    ostr2 << setprecision(3) << temperature_;
    string timestr(ostr1.str());
    string tempstr(ostr2.str());
    string buff;
    ofileName = ofileName + "." + timestr + "." + tempstr + ".ppm";
    ofpngname = ofpngname + "." + timestr
        + "." + tempstr + ".png";

    ///
    /// Open the output file
    ///

    ofstream out(ofileName.c_str());
    try {
        if (!out.is_open()) {
            throw FileException("Lattice","writeLatticePNG",
                                ofileName,"Could not open");
        }
    }
    catch (FileException fex) {
        fex.printException();
        exit(1);
    }

    ///
    /// Write PPM header for full color image
    ///

    out << "P3" << endl;
    out << xdim_ << " " << ydim_ << endl;
    out << COLORSATVAL << endl;

    unsigned int slice = zdim_/2;
    unsigned int nd,izz,valout;
    unsigned int sitenum;
    for (j = 0; j < ydim_; j++) {
        for (i = 0; i < xdim_; i++) {
           if (deptheffect_) {
               done = false;
               nd = 0;
               sitenum = getIndex(i,j,slice);
               izz = slice;
               do {
                   sitenum = getIndex(i,j,izz);
                   if (nd == 10 || site_[sitenum].getMicroPhaseId() > 1) {
                       done = true;
                   } else {
                       nd++;
                       izz++;
                       if (izz >= zdim_) izz -= zdim_;
                   }
               } while (!done);
               sitenum = getIndex(i,j,izz);
               image[i][j] = site_[sitenum].getMicroPhaseId();
               dshade[i][j] = 0.1 * (10.0 - ((double)nd));
           } else {
               sitenum = getIndex(i,j,slice);
               image[i][j] = site_[sitenum].getMicroPhaseId();
               dshade[i][j] = 1.0;
           }

        }
    }          

    double red,green,blue;
    vector<double> colors;
    for (j = 0; j < ydim_; j++) {
        for (i = 0; i < xdim_; i++) {
           colors = chemSys_->getColor(image[i][j]);
           red = dshade[i][j]*colors[0] + 0.5;
           green = dshade[i][j]*colors[1] + 0.5;
           blue = dshade[i][j]*colors[2] + 0.5;
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

void Lattice::writeDamageLatticePNG (double curtime, const string &root)
{
    unsigned int i,j,k;
    string ofileName(root);
    string ofbasename(root);
    string ofpngname(root);
    string ofpngbasename(root);

    vector<double> dumvec;
    vector<unsigned int> idumvec;
    vector<vector<unsigned int> > image;
    vector<vector<double> > dshade;
    dumvec.resize(ydim_,0.0);
    idumvec.resize(ydim_,0);
    dshade.resize(xdim_,dumvec);
    image.resize(xdim_,idumvec);
    bool done;

    ///
    /// Construct the name of the output file
    ///

    ostringstream ostr1,ostr2;
    ostr1 << (int)((curtime * 100.0) + 0.5);	// hundredths of an hour
    ostr2 << setprecision(3) << temperature_;
    string timestr(ostr1.str());
    string tempstr(ostr2.str());
    string buff;
    ofileName = ofileName + "." + timestr + "." + tempstr + ".ppm";
    ofpngname = ofpngname + "." + timestr
        + "." + tempstr + ".png";

    ///
    /// Open the output file
    ///

    ofstream out(ofileName.c_str());
    try {
        if (!out.is_open()) {
            throw FileException("Lattice","writeLatticePNG",
                                ofileName,"Could not open");
        }
    }
    catch (FileException fex) {
        fex.printException();
        exit(1);
    }

    ///
    /// Write PPM header for full color image
    ///

    out << "P3" << endl;
    out << ydim_ << " " << zdim_ << endl;
    out << COLORSATVAL << endl;

    unsigned int slice = zdim_/2;
    unsigned int nd,izz,valout;
    unsigned int sitenum;
    for (j = 0; j < ydim_; j++) {
        for (i = 0; i < xdim_; i++) {
           done = false;
           nd = 0;
           sitenum = getIndex(i,j,slice);
           izz = slice;
           do {
               sitenum = getIndex(i,j,izz);
               if (nd == 10 || site_[sitenum].getMicroPhaseId() > 1) {
                   done = true;
               } else {
                   nd++;
                   izz++;
                   if (izz >= zdim_) izz -= zdim_;
               }
           } while (!done);
           sitenum = getIndex(i,j,izz);
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

    double red,green,blue;
    vector<double> colors;
    for (j = 0; j < ydim_; j++) {
        for (i = 0; i < xdim_; i++) {
           colors = chemSys_->getColor(image[i][j]);
           red = dshade[i][j]*colors[0] + 0.5;
           green = dshade[i][j]*colors[1] + 0.5;
           blue = dshade[i][j]*colors[2] + 0.5;
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

void Lattice::makeMovie (const string &root)
{
    unsigned int i,j,k;
    string ofileName(root);
    string ofbasename(root);
    string ofgifileName(root);
    string ofgifbasename(root);

    vector<double> dumvec;
    vector<unsigned int> idumvec;
    vector<vector<unsigned int> > image;
    vector<vector<double> > dshade;
    dumvec.resize(ydim_,0.0);
    idumvec.resize(ydim_,0);
    dshade.resize(zdim_,dumvec);
    image.resize(zdim_,idumvec);
    bool done;

    ///
    /// Construct the name of the output file
    ///

    string buff;
    ostringstream ostr1,ostr2,ostr3;
    ostr1 << (int)(time_ * 100.0);	// hundredths of an hour
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
        ostr3 << (int)(k);	// x slice number
        string istr(ostr3.str());
        ofileName = ofbasename + "." + timestr + "."
            + tempstr + "." + istr + ".ppm";
        ofgifileName = ofgifbasename + "." + timestr
            + "." + tempstr + "." + istr + ".gif";

        ofstream out(ofileName.c_str());
        if (!out.is_open()) {
            throw FileException("Lattice","makeMovie",ofileName,"Could not open");
        }

        ///
        /// Write PPM header for full color image.
        ///

        out << "P3" << endl;
        out << xdim_ << " " << ydim_ << endl;
        out << COLORSATVAL << endl;

        unsigned int slice = k;
        unsigned int nd,izz,valout;
        unsigned int sitenum;
        for (j = 0; j < ydim_; j++) {
            for (i = 0; i < xdim_; i++) {
               done = false;
               nd = 0;
               sitenum = getIndex(i,j,slice);
               izz = slice;
               do {
                   sitenum = getIndex(i,j,izz);
                   if (nd == 10 || site_[sitenum].getMicroPhaseId() > 1) {
                       done = true;
                   } else {
                       nd++;
                       izz++;
                       if (izz >= zdim_) izz -= zdim_;
                   }
               } while (!done);
               sitenum = getIndex(i,j,izz);
               image[i][j] = site_[sitenum].getMicroPhaseId();
               dshade[i][j] = 0.1 * (10.0 - ((double)nd));
            }
        }          

        double red,green,blue;
        vector<double> colors;
        for (j = 0; j < ydim_; j++) {
            for (i = 0; i < xdim_; i++) {
               colors = chemSys_->getColor(image[i][j]);
               red = dshade[i][j]*colors[0] + 0.5;
               green = dshade[i][j]*colors[1] + 0.5;
               blue = dshade[i][j]*colors[2] + 0.5;
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

    buff = "gifsicle --delay=10 "
        + ofgifbasename + "*.gif > "
        + ofgifbasename + ".movie.gif";
    system(buff.c_str());
}

vector<int> Lattice::transform (int shrinkingid,
                                int netsites_shrinkingid,
                                int growingid,
                                int netsites_growingid,
                                double volumeratio)
{
    ///
    /// @todo Consider breaking this method into smaller pieces
    ///

    int expindex;

    vector<double> expval;
    expval.clear();
    expval.resize(3,0.0);

    vector<int> coordin;
    coordin.clear();
    coordin.resize(3,0);

    vector<Isite> diss;
    diss = interface_[shrinkingid].getDissolutionSites();
    if (verbose_) cout << "The number of sites of phase "
                       << shrinkingid << " to dissolve is: "
                       << diss.size() << endl;

    Site *ste;

    ///
    /// Construct the unordered list of sites to dissolve based only
    /// on the phase id and whether or not the total number of sites to dissolve
    /// has been reached.
    ///
    /// @remark Is this task biased by site position?  Has the list of sites been randomized?
    ///

    if (diss.size() < (-netsites_shrinkingid)) {
        for (int ii = 0; ii < numsites_; ii++) {
            ste = &site_[ii];
            if (ste->getMicroPhaseId() == shrinkingid) {
                addDissolutionSite(ste,shrinkingid);
            }
        }
    }

    diss = interface_[shrinkingid].getDissolutionSites();
    if (verbose_) cout << "New size of diss of phase " << shrinkingid
                       << " is: " << diss.size() << endl;

    double alreadygrown = 0.0;
    int numtransform = 0;
    int max = (int) volumeratio;
    
    for (int ii = (diss.size() - 1); ii > 0 && alreadygrown < netsites_growingid  
             && numtransform < (-netsites_shrinkingid); ii--) {
        ste = &site_[diss[ii].getId()];

        expval.clear();
        
        vector<Site *> porousneighbor, waterneighbor;
        porousneighbor.clear();
        waterneighbor.clear();
        for (int j = 0; j < ste->nbSize(1); j++) {
            Site *stenb;
            stenb = ste->nb(j);
            if (stenb->getMicroPhaseId() == WATERID) {
              waterneighbor.push_back(stenb);
            } else if (chemSys_->isPorous(stenb->getMicroPhaseId())) {
              porousneighbor.push_back(stenb);
            }
        }
        
        if (verbose_) {
            cout << "Having " << waterneighbor.size() << " water pixels and "
                 << porousneighbor.size() << " porous voxels in the neighborhood."
                 << endl;             
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
              if (alstenb->getMicroPhaseId() == WATERID) {
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
            if (verbose_) cout << "subbulk = " << subbulk << " GPa." << endl;   
            subbulk = subbulk * 1.0e3; // convert GPa to MPa
            double subsolidbulk = subbulk;
            subsolidbulk *= ((1 + (double)numWater / 27.0) / (1 - (double)numWater / 27.0));
            
            ///
            /// @remark This seems like a double conversion.  Hasn't the conversion been done?
            ///

            subsolidbulk = subsolidbulk * 1.0e3; // convert GPa to MPa         

            ///
            /// 3. Calculate crystallization strain in this sub volume;
            ///    porevolfrac is the volume fraction of pore space occupied by crystal
            ///    
            ///    @todo generalize porous phase porosities in the block below instead of 0.25
            ///

            double porevolfrac = 0.0;
            if (numWater != 0 || numPorous != 0) {
              porevolfrac = (double)(volumeratio) / (numWater + (numPorous * 0.25));
            } else {
              porevolfrac = 1.0;
            } 

            //  This is hard-wired right now
            //  @todo generalize crystallization pressure to more phases
           
            double exp = solut_->calculateCrystalStrain(SI_[growingid],
                                                    porevolfrac,
                                                    subbulk,
                                                    subsolidbulk);

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

                double dwmcval = chemSys_->getPorosity(growingid)
                                - chemSys_->getPorosity(WATERID);
                for (int j = 0; j < waterneighbor[i]->nbSize(2); j++) {
                    Site *nb = waterneighbor[i]->nb(j);
                    nb->dWmc(dwmcval);
                }
            }

            dWaterchange(volumeratio - (waterneighbor.size() + 1));
            alreadygrown += (volumeratio - (waterneighbor.size() + 1));
            count_.at(growingid) += (int)(volumeratio - (waterneighbor.size() + 1) + 0.5);

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

                double dwmcval = chemSys_->getPorosity(growingid)
                                 - chemSys_->getPorosity(WATERID);
                for (int j = 0; j < waterneighbor[i]->nbSize(2); j++) {
                    Site *nb = waterneighbor[i]->nb(j);
                    nb->dWmc(dwmcval);
                }
            }

        }
    }       // End of loop over all ettringite sites to form

    /*
    netsites.at(shrinkingid) += (int) numtransform;
    netsites.at(growingid) -= (int) alreadygrown;	
    */

    if (verbose_) {
        cout << "The number of aluminum phase " << shrinkingid
             << " transformed into ETTR is: " << numtransform << endl;
    }

    vector<int> numchanged;
    numchanged.clear();
    numchanged.resize(2,0);
    numchanged[0] = numtransform;
    numchanged[1] = (int)alreadygrown;
    if (verbose_) cout << "The number of alreadygrown is: " << alreadygrown << endl;

    return numchanged;
}

vector<unsigned int> Lattice::writeSubVolume (string fileName,
                                              Site *centerste,
                                              int size)
{
    ofstream out(fileName.c_str());
    
    out << "Version: 7.0" << endl;
    out << "X_Size: 3" << endl;
    out << "Y_Size: 3" << endl;
    out << "Z_Size: 3" << endl;
    out << "Image_Resolution: 1" << endl;

    vector<unsigned int> alnb = getNeighborhood(centerste->getId(),size);

    for (int j = 0; j < alnb.size(); j++) {
      int phaseid = site_[alnb[j]].getMicroPhaseId();
      out << phaseid << endl;
    }
    out.close();

    return alnb;
}

void Lattice::applyExp (vector<unsigned int> alnb,
                        double exp)
{
    Site *ste;
    if (exp > 0.0) {
      for (int i = 0; i < alnb.size(); i++) {
        ste = &site_[alnb[i]];
        ste->setExpansionStrain(exp);

        map<int,vector<double> >::iterator p = expansion_.find(ste->getId());
        if (p != expansion_.end()) {
          for (int j = 0; j < (p->second).size(); j++) {
            (p->second)[j] = ste->getExpansionStrain();
          }
        } else {
          vector<double> expval;
          expval.clear();
          expval.resize(3,ste->getExpansionStrain());
          expansion_.insert(make_pair(ste->getId(),expval));
        }
      
        map<int,vector<int> >::iterator pp = expansion_coordin_.find(ste->getId());
        if (pp == expansion_coordin_.end()) {
          vector<int> coordin;
          coordin.clear();
          coordin.resize(3,0);
          coordin[0] = ste->getX();
          coordin[1] = ste->getY();
          coordin[2] = ste->getZ();
          expansion_coordin_.insert(make_pair(ste->getId(),coordin));
        }
      }
    }
  
    return;
}

double Lattice::getSurfaceArea (int phaseid)
{
    double surface1 = 0.0, surface2 = 0.0;
    Site *ste, *stenb;
    vector<Isite> isite = interface_[phaseid].getDissolutionSites();

    ///
    /// Method 1: Surface area is equal to the wmc (prop to volume of surrounding pores)
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
                surface2 += chemSys_->getPorosity(stenb->getMicroPhaseId());
            }
        }
    }    

    if (verbose_) {
        cout << "surface area of phase " << phaseid << " calculated by method 1 is: "
             << surface1 << endl;
        cout << "surface area of phase " << phaseid << " calculated by method 2 is: "
             << surface2 << endl;
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
                                                   int sortorder = 0)
{
    list<Sitesize> distlist;
    list<Sitesize>::iterator it;
    bool found = false;

    Sitesize ss;
    int domainsize = 0;
    for (int i = 0; i < site_.size(); i++) {
        if (site_[i].getMicroPhaseId() == phaseid) {
            domainsize = findDomainSize(i,maxsize);
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
                            distlist.insert(it,ss);
                            found = true;
                        } else {
                            it++;
                        }
                    }
                } else {
                    it = distlist.begin();
                    while (it != distlist.end() && !found) {
                        if (ss.nsize <= (*it).nsize) {
                            distlist.insert(it,ss);
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
            if (distlist.size() > numsites) distlist.pop_back();
        }

    }

    return(distlist);
}

int Lattice::findDomainSize(int siteid, int maxsize) 
{

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

                if (site_[getIndex(ix,iy,iz)].getMicroPhaseId() == phaseid) {
                    nfound++;
                }
            }
        }
    }

    return nfound;
}   

