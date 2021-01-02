/**
@file Lattice.cc
@brief Definition of methods for the Lattice class.

*/
#include "Lattice.h"
#include "Interface.h"

Lattice::Lattice (ChemicalSystem *cs,
                  Solution *solut)
: siteneighbors_(18),chemsys_(cs),solut_(solut)
{
  xdim_ = ydim_ = zdim_ = 0;
  time_ = 0.0;
  temperature_ = REFTEMP;   // in Kelvin (see global.h)
  oldtemp_ = REFTEMP;       // in Kelvin (see global.h)
  numsites_ = 0;
  resolution_ = REFRES;     // in micrometers (see global.h)
  site_.clear();
  deptheffect_ = false;
}

Lattice::Lattice (ChemicalSystem *cs,
                  Solution *solut,
                  const string &fname)
: siteneighbors_(18),chemsys_(cs),solut_(solut)
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
  deptheffect_ = false;
    
  ///
  /// Open the microstructure input file and process it
  ///

  ifstream in(fname.c_str());
  try {
    if (!in) {
        throw FileException("Lattice","Lattice",fname,"Could not open");
    }
  }
  catch (FileException fex) {
    fex.printException();
    exit(1);
  }
    
  in >> buff;
  if (buff == VERSIONSTRING) {
    in >> version_;
    in >> buff;
    if (buff == XSIZESTRING) {
        in >> xdim_;
        in >> buff;
        in >> ydim_;
        in >> buff;
        in >> zdim_;
    } else if (buff == IMGSIZESTRING) {
        ydim_ = zdim_ = xdim_;
    }
    in >> buff;
    if (buff == IMGRESSTRING) {
        double testres;
        in >> testres;
        setResolution(testres);
    }
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
  numsites_ = (unsigned int)(xdim_ * ydim_ * zdim_);
  cout << "    numsites_ = " << numsites_ << endl;
  cout.flush();
  cout << "    resolution_ = " << resolution_ << endl;
  cout.flush();  
  
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
    cout << "Checking whether I can get a seed " << rg_->getSeed() << endl;
    cout.flush();
    cout << "Checking value of random number " << ", " << rg_->Ran3() << endl;
    cout.flush();

    count_.clear();
    count_.resize(chemsys_->getMicphasenum(),0);

    SI_.clear();
    SI_.resize(chemsys_->getMicphasenum(),1.0);

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
      site_[i].setPhaseId(pid);
      count_[pid] += 1;
    }
    
    ///
    /// Done with the input microstructure file, so close it.
    ///

    in.close();
 
    ///
    /// With the phase counts known, calculate phase volume fractions
    ///

    volumefraction_.clear();
    volumefraction_.resize(chemsys_->getMicphasenum(),0.0);
    cout << "Calculating Volume Fractions now..." << endl;
    try {
      if (site_.size() > 0) {
        for (ii = 0; ii < chemsys_->getMicphasenum(); ii++) {
            volumefraction_[ii] = ((double)count_[ii])/((double)site_.size());
            if (volumefraction_[ii] > 0.0) cout << "***Volume fraction["
                << ii << "] = " << volumefraction_[ii] << endl;
        }
      } else {
        msg = "Divide by zero error:  site_.size() = 0";
        throw FloatException("Lattice","Lattice",msg);
      }
      cout << "...Done!" << endl;
    }
    catch (FloatException fex) {
      fex.printException();
      exit(1);
    }
    
    cout << "Calculating weighted mean curvatures now..." << endl;
    for (ii = 0; ii < site_.size(); ii++) {
      site_[ii].calcWmc();
    }
    cout << "...Done!" << endl;

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
    
    ///
    /// Why is numphase assigned here?  It is a local variable and is not
    /// used.
    ///

    unsigned int numphase = chemsys_->getMicphasenum();

    site_.push_back(Site(x,y,z,xdim_,ydim_,zdim_,siteneighbors_,chemsys_));
}

void Lattice::findInterfaces ()
{
    unsigned int i,kk;
    unsigned int k;
    vector<Site *> gsite,dsite;
    
    ///
    /// An interface must have at least one adjacent site that is water or void
    ///

    interface_.clear();
    for (i = 0; i < chemsys_->getMicphasenum(); i++) { 
      cout << "  Database item " << i << ": ";
      cout.flush();
      if (i != WATERID && i != VOIDID) {
        cout << "Probing for interface... ";
        gsite.clear();
        dsite.clear();
        for (k = 0; k < site_.size(); k++) {
          if (site_[k].getWmc() > 0) {
            if ((site_[k].getPhaseId() == i)) {
              dsite.push_back(&site_[k]);
              site_[k].setDissolutionSite(i);
              for (kk = 0; kk < site_[k].nbSize(2); kk++) {
                if ((site_[k].nb(kk))->getPhaseId() == WATERID) {
                  gsite.push_back(site_[k].nb(kk));
                  site_[k].nb(kk)->setGrowthSite(i);
                }
              }
            } else if (chemsys_->isGrowthtemplate(i,site_[k].getPhaseId())) {
              for (kk = 0; kk < site_[k].nbSize(1); kk++) {
                if ((site_[k].nb(kk))->getPhaseId() == WATERID) {
                  gsite.push_back(site_[k].nb(kk));
                  site_[k].nb(kk)->setGrowthSite(i);
                }
              }
            }
          }
        }

        cout << " Done! " << dsite.size()
            << " dissolution sites and " << gsite.size()
            << " growth sites" << endl;

        if ((gsite.size() == 0) && (dsite.size() == 0)) {
          cout << "Testing phase " << i << " for nucleation ";

          ///
          /// We are dealing with a phase that may need to
          /// nucleate.  Identify the eligible nucleation
          //  sites for that phase
          ///

          double thresh = (0.5 / pow(resolution_,3.0));
          cout << "Thresh = " << thresh << endl;
          cout.flush();
          double g = 0.0;
          for (k = 0; k < site_.size(); k++) {
            if ((site_[k].getPhaseId() == WATERID)) {
              g = rg_->Ran3();
              if (g < thresh) {
                gsite.push_back(&site_[k]);
                site_[k].setGrowthSite(i);
              }
            }
          }
        }
        
        cout << "Trying to add a water interface for phase " << i << "... ";
        cout.flush();
        interface_.push_back(Interface(chemsys_,rg_,gsite,dsite,i));   
        cout << "Done!" << endl;
        cout.flush();
      } else {
        cout << "Trying to add a regular interface for phase " << i << "... ";
        cout.flush();
        interface_.push_back(Interface(rg_));   
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
    cout << "size of interface_[" << phaseid << "] is " << isite.size() << endl; 
    int numleft,numchange = 0;

    ///
    /// We need to go through the interface list in
    /// normal order, which is the reverse order that
    /// we use for dissolution.
    ///

    numleft = numtoadd;
    cout << "-->Phase " << phaseid
        << " needs to grow at " << numtoadd
        << " sites" << endl;

    while ((numleft > 0) && (isite.size() >= 1)) {
      for (i = 0; (numleft > 0) && (i < isite.size()); i++) {
        ste = &site_[isite[i].getId()];
        pid = ste->getPhaseId();
     
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

          dwmcval = chemsys_->getPorosity(phaseid)
                  - chemsys_->getPorosity(pid);
          setPhaseId(ste,phaseid);
          count_[phaseid] += 1;
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
            if (stenb->getPhaseId() == WATERID) {
              addGrowthSite(stenb,phaseid);
            } else if (stenb->getWmc() <= 0.0) {
              removeDissolutionSite(stenb,stenb->getPhaseId());
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
    cout << "size of interface_[" << phaseid << "] is " << isite.size() << endl; 
    try {
      for (i = 0; i < isite.size(); i++) {
        if (site_.at(isite[i].getId()).getPhaseId() != phaseid) {
          cout << "Uh-oh... interface " << phaseid
               << " is corrupted with phase "
               << site_.at(isite[i].getId()).getPhaseId() << endl;
          cout << "Offending site is " << isite[i].getId() << endl;
        }
      }
    }
    catch (out_of_range &oor) {
      EOBException ex("Lattice","dissolvePhase","site_",site_.size(),i);
      ex.printException();
      exit(1);
    }
    cout << "-->Phase " << phaseid
         << " needs to dissolve at " << numtotake
         << " sites" << endl;
    
    int numleft = numtotake;
    int numchange = 0;
    while ((numleft > 0) && (isite.size() > 1)) {
  
      try {
        for (i = isite.size() - 1; (i > 0) && (numleft > 0); i--) {
          ste = &site_.at(isite[i].getId());
          pid = ste->getPhaseId();
          removeDissolutionSite(ste,pid);
      
          setPhaseId(ste,WATERID);
          count_[WATERID] += 1;
      
          ///
          /// Weighted mean curvature (wmc) is changed by the difference
          /// between the growing phase's porosity and the template's porosity.
          ///
          /// @todo Determine why the calculation works this way.
          ///

          dwmcval = chemsys_->getPorosity(WATERID)
                  - chemsys_->getPorosity(pid);
          ste->dWmc(dwmcval);
      
          unsigned int j;
          for (j = 0; j < ste->nbSize(1); j++) {
            stenb = ste->nb(j);
            stenb->dWmc(dwmcval);
            if ((stenb->getPhaseId() != WATERID)
                && (stenb->getPhaseId() != VOIDID)) {
              int nbpid;
              nbpid = stenb->getPhaseId();

              ///
              /// Now that the site has been dissolved, it is eligible for growth
              /// later on, so we add it to the list of growth sites.
              ///

              addDissolutionSite(stenb,nbpid);
              addGrowthSite(ste,nbpid);
    
              vector<int> nbgrowthtemp;
              nbgrowthtemp.clear();
              nbgrowthtemp = chemsys_->getGrowthtemplate(nbpid);
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
        interface_.at(pid).removeDissolutionSite(ste,false);
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
        cout << "Going into Lattice::fillPorosity(" << -numsites << ")... " << endl;
        cout.flush();
        numemptied = fillPorosity(-numsites);
        return(-numemptied);
    }
    
    ///
    /// Finding all potential VOID sites.
    ///
    /// @todo Consider removing some of the standard output, or setting a flag for it.
    ///

    cout << "Finding and sorting all potential void sites ... ";
    cout.flush();
    list<Sitesize> distlist = findDomainSizeDistribution(WATERID,numsites,maxsearchsize,0);
    list<Sitesize>::iterator it;

    cout << "OK, found " << distlist.size()
         << " potential void sites." << endl;
    cout.flush();
    
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
            setPhaseId(site_[siteid].getId(),VOIDID);
	        count_[VOIDID] += 1;
	        count_[WATERID] -= 1;
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

    cout << "In fillPorosity (" << numsites << ")" << endl;
    cout.flush();
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

    cout << "Finding and sorting all potential water sites ... ";
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
            setPhaseId(site_[siteid].getId(),WATERID);
	        count_[WATERID] += 1;
	        count_[VOIDID] -= 1;
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
                    if (site_.at(getIndex(hx,hy,hz)).getPhaseId() == WATERID
                            || site_.at(getIndex(hx,hy,hz)).getPhaseId() == VOIDID) {
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

        cout << "Changing lattice resolution from ";
        cout << resolution_ << " to " << res << endl;
        cout.flush();
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
                                    bool isfirst,
                                    bool &capwater)
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

    vol_next = chemsys_->getMicphasevolume();
    vfrac_next.resize(vol_next.size(),0.0);
    phasenames = chemsys_->getMicphasename();

    /// JWB: 2020 Dec 22  Do manual adjustment of certain microstructure
    /// phase voxels in a cement paste, accounting for gel porosity, etc.
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
   
    adjustMicrostructureVolumeFractions(phasenames,vol_next,vfrac_next);

    ///
    /// Calculate number of sites of each phase in next state
    ///

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
      EOBException ex("Lattice","changeMicrostructure","phasenames",
                    phasenames.size(),i);
      ex.printException();
      exit(1);
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

        netsites.clear();
        netsites.resize(chemsys_->getMicphasenum(),0);
        pid.clear();
        pid.resize(chemsys_->getMicphasenum(),0);

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
            
            growing.push_back(chemsys_->getMicid("ETTR"));
            shrinking.resize(growing.size(),idummy);
            volratios.resize(growing.size(),ddummy);
            for (i = 0; i < growing.size(); ++i) {
                shrinking[i].push_back(chemsys_->getMicid("MONOSULPH"));
                volratios[i].push_back(2.288);
                shrinking[i].push_back(chemsys_->getMicid("AFMC"));
                volratios[i].push_back(2.699);
                shrinking[i].push_back(chemsys_->getMicid("HYDROTALC"));
                volratios[i].push_back(3.211);
            }

            for (int ii = 0; ii < growing.size(); ++ii) {
                for (i = 0; i < vfrac_next.size(); i++) {
                    cursites = (int)(count_.at(i) + 0.5);
                    newsites = (int)((numsites_ * vfrac_next.at(i)) + 0.5);
                    tnetsites = 0;
                    if (i != WATERID && i != VOIDID) tnetsites = (newsites - cursites);
                    netsites.at(i) = tnetsites;
                    pid.at(i) = i;
                    if (i == growing[ii] && isfirst) {
                      netsites.at(i) = 0;
                      count_.at(i) = newsites;
                    }        
                    if (netsites.at(i) != 0) cout << "****netsites["
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
                      cout << "start to crystal-pressure transform at time_ = " 
                       << time_ << endl;

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
                              cout << "netsites.at(" << shrinkid << ") is: "
                                   << netsites.at(shrinkid) << endl;
                              cout << "netsites.at(" << growid << ") is: "
                                   << netsites.at(growid) << endl;
                          }

                      }
                  }
            }
        }

        catch (out_of_range &oor) {
            EOBException ex("Lattice","changeMicrostructure",
                    "phasenames or count_ or pid or netsites or vfrac_next",
                    phasenames.size(),i);
            ex.printException();
            exit(1);
        }

    } else {
        
        ///
        /// Sulfate attack will NEVER be done during this simulation.
        ///  Normalize to get volume fractions and compute number
        ///  of sites of each phase needed.
        ///

        netsites.clear();
        netsites.resize(chemsys_->getMicphasenum(),0);
        pid.clear();
        pid.resize(chemsys_->getMicphasenum(),0);

        try {
            for (i = FIRST_SOLID; i < vfrac_next.size(); i++) {
                cursites = (int)(count_.at(i) + 0.5);
                newsites = (int)((numsites_
                             * vfrac_next.at(i)) + 0.5);
                tnetsites = 0;
                tnetsites = (newsites - cursites);
                netsites.at(i) = tnetsites;
                pid.at(i) = i;
                if (netsites.at(i) != 0) cout << "***netsites["
                                              << phasenames.at(i)
                                              << "] in this state = "
                                              << netsites.at(i)
                                              << "; cursites = " << cursites
                                              << " and newsites = " << newsites
                                              << endl;
            }
        }
        catch (out_of_range &oor) {
            EOBException ex("Lattice","changeMicrostructure",
                    "phasenames or count_ or pid or netsites or vfrac_next",
                     phasenames.size(),i);
            ex.printException();
            exit(1);
        }
    }
    
    cout << "Sorting non-pore sites... netsites[VOIDID] = netsites["
         << VOIDID << "] = " << netsites[VOIDID] << endl;
    cout.flush();

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
        EOBException ex("Lattice","changeMicrostructure","pid",
                     pid.size(),ii);
        ex.printException();
        exit(1);
    }

    Interface ifc;
    vector<Isite> gs,ds;
    cout << "Getting change vectors for non-pore phases... netsites[VOIDID] = "
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

    cout << "Switching phase id values..." << endl;
    try {
        for (i = FIRST_SOLID; i < netsites.size(); i++) {
            numadded = 0;
            numadded_actual = 0;
	  
	        if (netsites[i] < 0) {
                cout << "Going into dissolve_phase now... pid = "
                     << pid.at(i) << endl;
                cout.flush();
                numadded = dissolvePhase(pid.at(i),-netsites[i]);
                numadded_actual += numadded;
                cout << "...Done with dissolve_phase for phase "
                     << pid.at(i) << "!" << endl;
                cout.flush();
            } else if (netsites[i] > 0) {
                cout << "Going into grow_phase now... pid = "
                     << pid.at(i) << endl;
            
                cout.flush();
 
                numadded = growPhase(pid.at(i),netsites[i]);

                numadded_actual += numadded;
                int diff = netsites[i] - numadded;
                while (diff > 0) {
                    gs = interface_[pid.at(i)].getGrowthSites();
                    ds = interface_[pid.at(i)].getDissolutionSites();
                    cout << "gs.size() = " << gs.size() << endl;
                    cout << "ds.size() = " << ds.size() << endl;
                    if (gs.size() == 0) {
                        cout << "Phase " << pid.at(i) << " needs to be nucleated. " << endl;
                        int nuclei = diff;
	                    double thresh = (nuclei / (volumefraction_[WATERID] * site_.size()));
	                    cout << "Thresh = " << thresh << endl;
	                    cout.flush();
                        if (thresh < 1) {
    	                    double g = 0.0;
    	                    for (int k = 0; k < site_.size(); k++) {
    	                        if (site_[k].getPhaseId() == WATERID) {
    	                            g = rg_->Ran3();
    	                            if (g < thresh) addGrowthSite(&site_[k],pid.at(i));
    	                        }
    	                    }
                        } else {
                            cout << "There is no room to grow, so exit the program." << endl;
                            exit(1);
                        }
	                }
                    numadded = growPhase(pid.at(i),diff);
    	            numadded_actual += numadded;
                    diff = diff - numadded;
	                cout << "diff = " << diff << endl;
	            }
    	        cout << "...Done with grow_phase for phase "
                     << pid.at(i) << "!" << endl;
                cout.flush();
            }
        
            if (numadded_actual*numadded_actual != netsites[i]*netsites[i]) {
                cout << "WARNING: Needed to switch on "
                     << netsites[i] << " of phase " << pid.at(i) << endl;
                cout << "         But actually did "
                     << numadded << " switches" << endl;
            }
        }
    }

    catch (out_of_range &oor) {
        EOBException ex("Lattice","changeMicrostructure","pid",
                    pid.size(),ii);
        ex.printException();
        exit(1);
    }
    cout << "time_ in the Lattice is: " << time_ << endl;

    cout << "Now emptying the necessary pore volume:  " << endl;

    cursites = (int)(count_.at(VOIDID) + 0.5);
    newsites = (int)((numsites_
                 * vfrac_next.at(VOIDID)) + 0.5);
    wcursites = (int)(count_.at(WATERID) + 0.5);
    wnewsites = (int)(wcursites - (newsites-cursites));
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
 
    cout << "Target volume fraction of water WAS " << vfrac_next.at(WATERID) << endl;
    vfrac_next.at(WATERID) -= ((double)(newsites - cursites)/(double)(numsites_));
    cout << "But WILL BE " << vfrac_next.at(WATERID) << " after creating void space" << endl;

    int numempty = emptyPorosity(newsites - cursites);
    cout << "Number actually emptied was:  " << numempty << endl;

    /// 
    /// Report on target and actual mass fractions
    ///

    cout << "*******************************" << endl;
    try {
        for (i = 0; i < vfrac_next.size(); i++) {
            volumefraction_.at(i) = ((double)(count_.at(i)))/((double)(site_.size()));
            cout << "Phase " << i << " Target volume fraction was "
                 << vfrac_next[i] << " and actual is "
                 << volumefraction_.at(i) << endl;
        }
    }

    catch (out_of_range &oor) {
        EOBException ex("Lattice","changeMicrostructure",
                    "volumefraction_ or count_",volumefraction_.size(),i);
        ex.printException();
        exit(1);
    }
    cout << "*******************************" << endl;
    
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

    double surfa = getSurfaceArea(chemsys_->getMicid("CSH"));
    
    if (volumefraction_[WATERID] <= 0.0) capwater = false;

    return;
}

void Lattice::adjustMicrostructureVolumeFractions (vector<string> name,
                                                   vector<double> &vol,
                                                   vector<double> &vfrac)
{
    int i = 0;
    double cap_watervolume = 0.0;
    double cap_spacevolume = 0.0;
    double cap_voidvolume = 0.0;
    double gel_watervolume = 0.0;
    double tot_watervolume = 0.0;

    try {
        
        // Initialize the volume fractions
        
        vfrac.clear();
        vfrac.resize(vol.size(),0.0);

        // Find the total system volume according to GEMS
        
        double GEM_thinks_totmicvol = chemsys_->getMictotvolume();
        double GEM_thinks_totmicinitvol = chemsys_->getMictotinitvolume();

        double GEM_micvols_addto = 0.0;

        for (i = 0; i < vol.size(); ++i) {
            GEM_micvols_addto += vol.at(i);
        }

        if (GEM_micvols_addto <= 0.0) {
            throw DataException("Lattice",
                                "adjustMicrostructureVolumeFractions",
                                "totvolume is NOT positive");
        }

        // Adjust the volume of CSH to
        // account for its saturated gel porosity
           
        int cshid = chemsys_->getMicid("CSH");

        tot_watervolume = vol.at(WATERID);
        cap_voidvolume = vol.at(VOIDID);

        gel_watervolume = (chemsys_->getPorosity(cshid)) * vol.at(cshid);
        cap_watervolume = tot_watervolume - gel_watervolume;
        vol.at(cshid) += gel_watervolume;

        // The capillary SPACE volume is the saturated + unsaturated capillary volume
        cap_spacevolume = cap_watervolume + cap_voidvolume;

        double spaceforwater = cap_spacevolume + gel_watervolume;

        if (!(chemsys_->isSaturated())) {   // System is sealed

            // Create void space due to self-desiccation
            cout << "Now deciding how much empty porosity should be made" << endl;
            cout << "%%%%%" << endl;
            cout << "Current microstructure water site fraction = "
                 << ((double)(count_.at(WATERID))) / ((double)(site_.size())) << endl;
            cout << "Current GEM water volume fraction = "
                 << (tot_watervolume / GEM_thinks_totmicvol) << endl;
            cout << "Current GEM water volume = " << tot_watervolume << endl;
            cout << "   This is divided between 1. " << gel_watervolume << " in gel pores" << endl;
            cout << "                           2. " << cap_watervolume
                                                     << " in capillary pores" << endl;
            cout << "                          (3. Plus " << cap_voidvolume <<
                                                       " void space)" << endl;

            double volchange = GEM_thinks_totmicvol - GEM_thinks_totmicinitvol;
            if (volchange < 0) {
                vol.at(VOIDID) = -volchange;
                vol.at(WATERID) = cap_watervolume;
            } else {
                vol.at(VOIDID) = 0.0;
                vol.at(WATERID) = cap_watervolume;
            }
        } else {
            vol.at(VOIDID) = 0.0;
            vol.at(WATERID) = cap_watervolume;
        }

        double totvolume = GEM_thinks_totmicvol;
        cout << "Total volume increases from " << totvolume;
        totvolume += vol.at(VOIDID);
        cout << " to " << totvolume << endl;

        cout << "%%%%%" << endl;

        // Calculate the volume fractions
          
        double sumvolfrac = 0.0;
        for (int i = 0; i < vol.size(); ++i) {
            if (i != WATERID) {
                vfrac.at(i) = (vol.at(i) / totvolume);
                sumvolfrac += vfrac.at(i);
            }
        }
        vfrac.at(WATERID) = 1.0 - sumvolfrac;

        ///
        /// End of manual adjustment  
        ///

        cout << "Total microstructure volume now needs to be: " << totvolume << endl;
        for (i = 0; i < vol.size(); i++) {
            cout << "    " << name.at(i) << "  V needs to be = " << vol.at(i) << ", Vtot = "
                     << totvolume << ", and volume fraction needs to be = " << vfrac.at(i) << endl;
        }
        cout.flush();
        
    }

    catch (DataException dex) {
        dex.printException();
        exit(1);
    }

    catch (out_of_range &oor) {
      EOBException ex("Lattice","adjustMicrostructureVolumeFractions","name",name.size(),i);
      ex.printException();
      exit(1);
    }

    return;
}

void Lattice::writeLattice (double curtime, const int simtype, const string &root)
{
    unsigned int i,j,k;
    string ofname(root);
    ostringstream ostr1,ostr2;
    ostr1 << setfill('0') << setw(6) << (int)((curtime * 24.0 * 60.0) + 0.5);	// minutes
    ostr2 << setprecision(3) << temperature_;
    string timestr(ostr1.str());
    string tempstr(ostr2.str());
    ofname = ofname + "." + timestr + "." + tempstr + ".img";
    cout << "    In Lattice::writeLattice, curtime = " << curtime << ", timestr = " << timestr << endl;
    cout.flush();

    ofstream out(ofname.c_str());
    try {
        if (!out.is_open()) {
            throw FileException("Lattice","writeLattice",ofname,"Could not open");
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
               out << site_[index].getPhaseId() << endl;
            }
        }
    }

    out.close();
    
    // The next block is implemented only if we are dealing with sulfate attack
    if (simtype == SULFATE_ATTACK) {
    
        ofname = root;
        ofname = ofname + "." + timestr + "." + tempstr + ".img.damage";

        ofstream out1(ofname.c_str());
        try {
            if (!out1.is_open()) {
                throw FileException("Lattice","writeLattice",ofname,"Could not open");
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
 
        int DAMAGEID = chemsys_->getMicid("DAMAGE"); 
        for (k = 0; k < zdim_; k++) {
            for (j = 0; j < ydim_; j++) {
                for (i = 0; i < xdim_; i++) {
                   int index = getIndex(i,j,k);
                   if (site_[index].IsDamage()) {
                     out1 << DAMAGEID << endl;
                   } else {
                     out1 << site_[index].getPhaseId() << endl;
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
    string ofname(root);
    ostringstream ostr1,ostr2;
    ostr1 << (int)((curtime * 100.0) + 0.5);	// hundredths of a day
    ostr2 << setprecision(3) << temperature_;
    string timestr(ostr1.str());
    string tempstr(ostr2.str());
    ofname = ofname + "." + timestr + "." + tempstr + ".img";

    ofstream out(ofname.c_str());
    try {
        if (!out.is_open()) {
            throw FileException("Lattice","writeLattice",ofname,"Could not open");
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
    string ofname(root);
    string ofbasename(root);
    string ofpngname(root);
    string ofpngbasename(root);

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

    ostringstream ostr1,ostr2;
    ostr1 << setfill('0') << setw(6) << (int)((curtime * 24.0 * 60.0) + 0.5);	// minutes
    ostr2 << setprecision(3) << temperature_;
    string timestr(ostr1.str());
    string tempstr(ostr2.str());
    string buff;
    ofname = ofname + "." + timestr + "." + tempstr + ".ppm";
    ofpngname = ofpngname + "." + timestr
        + "." + tempstr + ".png";

    ///
    /// Open the output file
    ///

    ofstream out(ofname.c_str());
    try {
        if (!out.is_open()) {
            throw FileException("Lattice","writeLatticePNG",
                                ofname,"Could not open");
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

    unsigned int slice = xdim_/2;
    unsigned int nd,ixx,valout;
    unsigned int sitenum;
    for (k = 0; k < zdim_; k++) {
        for (j = 0; j < ydim_; j++) {
           if (deptheffect_) {
               done = false;
               nd = 0;
               sitenum = getIndex(slice,j,k);
               ixx = slice;
               do {
                   sitenum = getIndex(ixx,j,k);
                   if (nd == 10 || site_[sitenum].getPhaseId() > 1) {
                       done = true;
                   } else {
                       nd++;
                       ixx++;
                       if (ixx >= xdim_) ixx -= xdim_;
                   }
               } while (!done);
               sitenum = getIndex(ixx,j,k);
               image[j][k] = site_[sitenum].getPhaseId();
               dshade[j][k] = 0.1 * (10.0 - ((double)nd));
           } else {
               sitenum = getIndex(slice,j,k);
               image[j][k] = site_[sitenum].getPhaseId();
               dshade[j][k] = 1.0;
           }

        }
    }          

    double red,green,blue;
    vector<double> colors;
    for (k = 0; k < zdim_; k++) {
        for (j = 0; j < ydim_; j++) {
           colors = chemsys_->getColor(image[j][k]);
           red = dshade[j][k]*colors[0] + 0.5;
           green = dshade[j][k]*colors[1] + 0.5;
           blue = dshade[j][k]*colors[2] + 0.5;
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

    buff = "convert " + ofname + " " + ofpngname;
    system(buff.c_str());
    return;
}

void Lattice::writeDamageLatticePNG (double curtime, const string &root)
{
    unsigned int i,j,k;
    string ofname(root);
    string ofbasename(root);
    string ofpngname(root);
    string ofpngbasename(root);

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

    ostringstream ostr1,ostr2;
    ostr1 << (int)((curtime * 100.0) + 0.5);	// hundredths of an hour
    ostr2 << setprecision(3) << temperature_;
    string timestr(ostr1.str());
    string tempstr(ostr2.str());
    string buff;
    ofname = ofname + "." + timestr + "." + tempstr + ".ppm";
    ofpngname = ofpngname + "." + timestr
        + "." + tempstr + ".png";

    ///
    /// Open the output file
    ///

    ofstream out(ofname.c_str());
    try {
        if (!out.is_open()) {
            throw FileException("Lattice","writeLatticePNG",
                                ofname,"Could not open");
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

    unsigned int slice = xdim_/2;
    unsigned int nd,ixx,valout;
    unsigned int sitenum;
    for (k = 0; k < zdim_; k++) {
        for (j = 0; j < ydim_; j++) {
           done = false;
           nd = 0;
           sitenum = getIndex(slice,j,k);
           ixx = slice;
           do {
               sitenum = getIndex(ixx,j,k);
               if (nd == 10 || site_[sitenum].getPhaseId() > 1) {
                   done = true;
               } else {
                   nd++;
                   ixx++;
                   if (ixx >= xdim_) ixx -= xdim_;
               }
           } while (!done);
           sitenum = getIndex(ixx,j,k);
           if (site_[sitenum].IsDamage()) {
               image[j][k] = 1;
           } else {
               image[j][k] = 5;
           }
           dshade[j][k] = 1.0;
           /*
           dshade[j][k] = 0.1 * (10.0 - ((double)nd));
           */
        }
    }          

    double red,green,blue;
    vector<double> colors;
    for (k = 0; k < zdim_; k++) {
        for (j = 0; j < ydim_; j++) {
           colors = chemsys_->getColor(image[j][k]);
           red = dshade[j][k]*colors[0] + 0.5;
           green = dshade[j][k]*colors[1] + 0.5;
           blue = dshade[j][k]*colors[2] + 0.5;
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

    buff = "convert " + ofname + " " + ofpngname;
    system(buff.c_str());
    return;
}

void Lattice::makeMovie (const string &root)
{
    unsigned int i,j,k;
    string ofname(root);
    string ofbasename(root);
    string ofgifname(root);
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

    for (i = 10; i < xdim_; i++) {

        ///
        /// Open the output file.
        ///

        ostr3.clear();
        ostr3 << (int)(i);	// x slice number
        string istr(ostr3.str());
        ofname = ofbasename + "." + timestr + "."
            + tempstr + "." + istr + ".ppm";
        ofgifname = ofgifbasename + "." + timestr
            + "." + tempstr + "." + istr + ".gif";

        ofstream out(ofname.c_str());
        if (!out.is_open()) {
            throw FileException("Lattice","makeMovie",ofname,"Could not open");
        }

        ///
        /// Write PPM header for full color image.
        ///

        out << "P3" << endl;
        out << ydim_ << " " << zdim_ << endl;
        out << COLORSATVAL << endl;

        unsigned int slice = i;
        unsigned int nd,ixx,valout;
        unsigned int sitenum;
        for (k = 0; k < zdim_; k++) {
            for (j = 0; j < ydim_; j++) {
               done = false;
               nd = 0;
               sitenum = getIndex(slice,j,k);
               ixx = slice;
               do {
                   sitenum = getIndex(ixx,j,k);
                   if (nd == 10 || site_[sitenum].getPhaseId() > 1) {
                       done = true;
                   } else {
                       nd++;
                       ixx++;
                       if (ixx >= xdim_) ixx -= xdim_;
                   }
               } while (!done);
               sitenum = getIndex(ixx,j,k);
               image[j][k] = site_[sitenum].getPhaseId();
               dshade[j][k] = 0.1 * (10.0 - ((double)nd));
            }
        }          

        double red,green,blue;
        vector<double> colors;
        for (k = 0; k < zdim_; k++) {
            for (j = 0; j < ydim_; j++) {
               colors = chemsys_->getColor(image[j][k]);
               red = dshade[j][k]*colors[0] + 0.5;
               green = dshade[j][k]*colors[1] + 0.5;
               blue = dshade[j][k]*colors[2] + 0.5;
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

        buff = "convert " + ofname + " " + ofgifname;
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
    cout << "The number of sites of phase " << shrinkingid << " to dissolve is: "
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
            if (ste->getPhaseId() == shrinkingid) {
                addDissolutionSite(ste,shrinkingid);
            }
        }
    }

    diss = interface_[shrinkingid].getDissolutionSites();
    cout << "New size of diss of phase " << shrinkingid << " is: " << diss.size() << endl;

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
            if (stenb->getPhaseId() == WATERID) {
              waterneighbor.push_back(stenb);
            } else if (chemsys_->isPorous(stenb->getPhaseId())) {
              porousneighbor.push_back(stenb);
            }
        }
        cout << "Having " << waterneighbor.size() << " water pixels and "
             << porousneighbor.size() << " porous voxels in the neighborhood."
             << endl;             

        if ((waterneighbor.size() + 1) <= max) { // count the site itself

            /// Expansion should occur
            ///
            /// 1. Take the subvolume centered on aluminate site that will dissolve
            ///

            string fname(jobroot_ + "_alsubvol.dat");
            vector<unsigned int> alnb = writeSubVolume(fname, ste, 1);
            ste->setDamage();  
            int numWater, numPorous;
            numWater = numPorous = 0;
            for (int nb = 0; nb < alnb.size(); nb++) {
              Site *alstenb = &site_[alnb[nb]];
              if (alstenb->getPhaseId() == WATERID) {
                numWater++;
              } else if (chemsys_->isPorous(alstenb->getPhaseId())) {
                numPorous++;
                alstenb->setDamage();
              } else if (chemsys_->isWeak(alstenb->getPhaseId())) {
                alstenb->setDamage();
              }
            }

            ///
            /// 2. Calculate the effective bulk modulus of this subvolume
            ///

            double subbulk = FEsolver_->getBulkModulus(fname);
            cout << "subbulk = " << subbulk << " GPa." << endl;   
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
           
            double exp = solut_->calculateCrystrain(SI_[growingid],
                                                    porevolfrac,
                                                    subbulk,
                                                    subsolidbulk);

            ///
            /// 4. Apply expansion strain on each voxel in this sub volume
            ///
     
            applyExp(alnb, exp);

            setPhaseId(ste, growingid);

            ///
            /// The Al-bearing phase has dissolved, so remove it from the
            /// list of dissolution sites of this phase
            ///

            removeDissolutionSite(ste, shrinkingid);
            numtransform++;
            alreadygrown++;
            count_.at(growingid)++;
            
            for (int i = 0; i < waterneighbor.size(); i++) {
                setPhaseId(waterneighbor[i], growingid);

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

                double dwmcval = chemsys_->getPorosity(growingid)
                                - chemsys_->getPorosity(WATERID);
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

            setPhaseId(ste, growingid);

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
                setPhaseId(waterneighbor[i], growingid);
                removeGrowthSite(waterneighbor[i], growingid);
                alreadygrown++;
                count_.at(growingid)++;

                ///
                /// Weighted mean curvature (wmc) is changed by the difference
                /// between the growing phase's porosity and the template's porosity.
                ///
                /// @todo Determine why the calculation works this way.
                ///

                double dwmcval = chemsys_->getPorosity(growingid)
                                 - chemsys_->getPorosity(WATERID);
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

    cout << "The number of aluminum phase " << shrinkingid
         << " transformed into ETTR is: " << numtransform << endl;

    vector<int> numchanged;
    numchanged.clear();
    numchanged.resize(2,0);
    numchanged[0] = numtransform;
    numchanged[1] = (int)alreadygrown;
    cout << "The number of alreadygrown is: " << alreadygrown << endl;

    return numchanged;
}

vector<unsigned int> Lattice::writeSubVolume (string fname,
                                              Site *centerste,
                                              int size)
{
    ofstream out(fname.c_str());
    
    out << "Version: 7.0" << endl;
    out << "X_Size: 3" << endl;
    out << "Y_Size: 3" << endl;
    out << "Z_Size: 3" << endl;
    out << "Image_Resolution: 1" << endl;

    vector<unsigned int> alnb = getNeighborhood(centerste->getId(),size);

    for (int j = 0; j < alnb.size(); j++) {
      int phaseid = site_[alnb[j]].getPhaseId();
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
        if (ste->getPhaseId() == phaseid) {
            for (int j = 0; j < ste->nbSize(1); j++) {
                stenb = ste->nb(j);
                surface2 += chemsys_->getPorosity(stenb->getPhaseId());
            }
        }
    }    

    cout << "surface area of phase " << phaseid << " calculated by method 1 is: "
         << surface1 << endl;
    cout << "surface area of phase " << phaseid << " calculated by method 2 is: "
         << surface2 << endl;

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
        if (site_[i].getPhaseId() == phaseid) {
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

    int phaseid = (site_[siteid]).getPhaseId();
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

                if (site_[getIndex(ix,iy,iz)].getPhaseId() == phaseid) {
                    nfound++;
                }
            }
        }
    }

    return nfound;
}   

