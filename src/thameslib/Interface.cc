/**
@file Interface.cc
@brief Definition of methods for the Interface class.

*/
#include "Interface.h"
#include "ChemicalSystem.h"
#include "RanGen.h"

bool cmp (const Site *s1,
          const Site *s2)
{
    return s1->getWmc() < s2->getWmc();
}

bool affinitySort (const Isite s1,
                   const Isite s2)
{
    return s1.getAffinity() > s2.getAffinity();
}

Interface::Interface ()
{
    microPhaseId_ = 0;
    growthSites_.clear();
    dissolutionSites_.clear();
}

Interface::Interface (RanGen *rg,
                      const bool verbose)
{
    microPhaseId_ = 0;
    growthSites_.clear();
    dissolutionSites_.clear();
    rg_ = rg;

    #ifdef DEBUG
        verbose_ = true;
    #else
        verbose_ = verbose;
    #endif

}

Interface::Interface (ChemicalSystem *csys,
                      RanGen *rg,
                      vector<Site *> gv,
                      vector<Site *> dv,
                      unsigned int pid,
                      const bool verbose)
{
    unsigned long int j;
    unsigned int i;
    int afty;

    vector<Site *>::iterator beginLocation,endLocation;
    vector<Isite>::iterator start,end;

    #ifdef DEBUG
        verbose_ = true;
    #else
        verbose_ = verbose;
    #endif

    rg_ = rg;
    microPhaseId_ = pid;
    chemSys_ = csys;

    dissolutionSites_.clear();
    growthSites_.clear();

    ///
    /// Eliminate duplicate values from the growth site vector
    ///

    if (gv.size() > 0) {
        sort(gv.begin(),gv.end());
        i = 1;
        beginLocation = gv.begin();
        endLocation = unique(gv.begin(),gv.end());
        gv.erase(endLocation,gv.end());
    }


    ///
    /// Now do the same thing to the dissolution site vector
    ///

    if (dv.size() > 0) {
        sort(dv.begin(),dv.end());
        i = 1;
        beginLocation = dv.begin();
        endLocation = unique(dv.begin(),dv.end());
        dv.erase(endLocation,dv.end());
    }

    ///
    /// Now sort the growth sites according to the affinity
    ///

    for (j = 0; j < gv.size(); j++) {
        afty = 0;
        for (i = 0; i < gv[j]->nbSize(2); i++) {
            afty += chemSys_->getAffinity(pid,gv[j]->nb(i)->getMicroPhaseId());
        }

        ///
        /// Add to the list of Isites.  An Isite is an object consisting
        /// of a pointer to a site and an affinity value
        ///

        growthSites_.push_back(Isite(gv[j]->getId(),afty));
    }

    if (growthSites_.size() > 0) {
        start = growthSites_.begin();
        end = growthSites_.end();

        ///
        /// The sort is built in to the STL for vectors or portions of vectors,
        /// as long as you give it a comparison function, which in this case is the
        /// `affinitySort` function already defined
        ///
        sort(start,end,affinitySort);
    }

    ///
    /// At this point the sites are perfectly sorted in descending order of affinity.
    /// Now shuffle the sorted sites according to the random growth factor.
    /// Each phase has a certain amount of randomness to its growth, which can be
    /// accessed through the `getRandomgrowth` method of the ChemicalSystem object.
    ///

    unsigned long int site1,site2;
    unsigned long int numgsites = growthSites_.size();

    for (j = 0; j < (chemSys_->getRandomGrowth(pid) * numgsites); j++) {

       // Choose two sites at random and switch their places
       site1 = (unsigned long int)(rg_->Ran3() * numgsites);
       site2 = (unsigned long int)(rg_->Ran3() * numgsites);
       swap(growthSites_[site1],growthSites_[site2]);
    }

    ///
    /// Now sort the dissolution sites according to the affinity
    ///

    try {
        for (j = 0; j < dv.size(); j++) {
            afty = 0;
            for (i = 0; i < dv[j]->nbSize(2); i++) {
                afty += chemSys_->getAffinity(pid,dv[j]->nb(i)->getMicroPhaseId());
            }
            dissolutionSites_.push_back(Isite(dv[j]->getId(),afty));
        }
    }
    catch (EOBException e) {
        e.printException();
        exit(0);
    }

    if (dissolutionSites_.size() > 0) {
        start = dissolutionSites_.begin();
        end = dissolutionSites_.end();
        sort(start,end,affinitySort);
    }

    numgsites = dissolutionSites_.size();
    
    ///
    /// The dissolution sites are perfectly ordered by affinity, just like the
    /// growth sites were.  Now add the randomness factor for this phase.
    ///

    try {
        for (j = 0; j < (chemSys_->getRandomGrowth(pid) * numgsites); j++) {

           // Choose two sites at random and switch their places
           site1 = (unsigned long int)(rg_->Ran3() * numgsites);
           site2 = (unsigned long int)(rg_->Ran2() * numgsites);
           swap(dissolutionSites_[site1],dissolutionSites_[site2]);
        }
    }
    catch (EOBException e) {
        e.printException();
        exit(0);
    }
}       // End of constructors

Interface::~Interface ()
{
    growthSites_.clear();
    dissolutionSites_.clear();
}

bool Interface::addGrowthSite (Site *loc)
{
    bool answer = false;
    bool found = false;
    unsigned int i;
    vector<Isite>::iterator p,q,start,end;
    start = growthSites_.begin();
    end = growthSites_.end();

    ///
    /// See if site is already present
    ///

    for (i = 0; (i < growthSites_.size()) && (!found); i++) {
        if (loc->getId() == growthSites_[i].getId()) found = true;
    }

    ///
    /// Add the site only if it was not found already
    ///

    if (!found) {
        int afty = 0;
        for (i = 0; i < loc->nbSize(2); i++) {
            afty += chemSys_->getAffinity(microPhaseId_,loc->nb(i)->getMicroPhaseId());
        }
        Isite tisite(loc->getId(),afty);
        q = lower_bound(start,end,tisite,affinitySort);
        growthSites_.insert(q,tisite);
        answer = true;
    }

    return answer;
}

bool Interface::addDissolutionSite (Site *loc)
{
    bool answer = false;
    bool found = false;
    unsigned int i;
    vector<Isite>::iterator p,q,start,end;
    start = dissolutionSites_.begin();
    end = dissolutionSites_.end();

    ///
    /// See if site is already present
    ///

    for (i = 0; (i < dissolutionSites_.size()) && (!found); i++) {
        if (loc->getId() == dissolutionSites_[i].getId()) found = true;
    }

    ///
    /// Add the site only if it was not found already
    ///

    if (!found) {
        int afty = 0;
        for (i = 0; i < loc->nbSize(2); i++) {
            afty += chemSys_->getAffinity(microPhaseId_,loc->nb(i)->getMicroPhaseId());
        }
        Isite tisite(loc->getId(),afty);
        q = lower_bound(start,end,tisite,affinitySort);
        dissolutionSites_.insert(q,tisite);
        answer = true;
    }

    return answer;
}

bool Interface::sortGrowthSites (vector<Site> &ste,
                                 unsigned int pid)
{
    unsigned int i,j;
    int afty;
    Site gs;

    ///
    /// The list of growth sites already exists or has been constructed, so we
    /// only need to update the affinities for each site
    ///

    for (j = 0; j < growthSites_.size(); j++) {
        afty = 0;
        gs = ste[growthSites_[j].getId()];
        for (i = 0; i < gs.nbSize(2); i++) {
            afty += chemSys_->getAffinity(pid,gs.nb(i)->getMicroPhaseId());
        }
        growthSites_[j].setAffinity(afty);
    }

    ///
    /// Now we sort the list of growth sites in descending order of affinity.
    /// Remember that `affinitySort` is the comparison function, already defined
    /// in this class, that must be passed to the STL sort function.
    ///

    if (growthSites_.size() > 0) {
        vector<Isite>::iterator start,end;
        start = growthSites_.begin();
        end = growthSites_.end();
        sort(start,end,affinitySort);
    }

    ///
    /// At this point the growth sites are perfectly sorted in descending
    /// order of affinity.  Next, we shuffle this sorting somewhat depending
    /// on how much randomization of growth sites is indicated for this particular
    /// phase.  The amount of randomness is obtained from the `getRandomgrowth` method
    /// of the ChemicalSystem object for this simulation
    ///

    unsigned long int site1,site2;
    unsigned long int numgsites = growthSites_.size();
    for (j = 0; j < (chemSys_->getRandomGrowth(pid) * numgsites); j++) {

       ///
       /// Choose two sites at random and switch their places
       ///
 
       site1 = (unsigned long int)(rg_->Ran3() * numgsites);
       site2 = (unsigned long int)(rg_->Ran3() * numgsites);
       swap(growthSites_[site1],growthSites_[site2]);
    }

    return true;  // successful sorting
}

bool Interface::sortDissolutionSites (vector<Site> &ste,
                                      unsigned int pid)
{
    unsigned int i,j;
    int afty;
    Site ds;

    ///
    /// The list of dissolution sites already exists or has been constructed, so we
    /// only need to update the affinities for each site
    ///

    for (j = 0; j < dissolutionSites_.size(); j++) {
       afty = 0;
       ds = ste[dissolutionSites_[j].getId()];
       for (i = 0; i < ds.nbSize(2); i++) {
           afty += chemSys_->getAffinity(pid,ds.nb(i)->getMicroPhaseId());
       }
       dissolutionSites_[j].setAffinity(afty);
    }

    ///
    /// Now we sort the list of dissolution sites in descending order of affinity.
    /// Remember that `affinitySort` is the comparison function, already defined
    /// in this class, that must be passed to the STL sort function.
    ///

    if (dissolutionSites_.size() > 0) {
        vector<Isite>::iterator start,end;
        start = dissolutionSites_.begin();
        end = dissolutionSites_.end();
        sort(start,end,affinitySort);
    }

    ///
    /// At this point the dissolution sites are perfectly sorted in descending
    /// order of affinity.  Next, we shuffle this sorting somewhat depending
    /// on how much randomization of growth sites is indicated for this particular
    /// phase.  The amount of randomness is obtained from the `getRandomgrowth` method
    /// of the ChemicalSystem object for this simulation
    ///

    unsigned long int site1,site2;
    unsigned long int numdsites = dissolutionSites_.size();
    for (j = 0; j < (chemSys_->getRandomGrowth(pid) * numdsites); j++) {

       // Choose two sites at random and switch their places
       site1 = (unsigned long int)(rg_->Ran3() * numdsites);
       site2 = (unsigned long int)(rg_->Ran3() * numdsites);
       swap(dissolutionSites_[site1],dissolutionSites_[site2]);
    }

    return true;
}

bool Interface::removeGrowthSite (Site *loc)
{
    bool found = false;
    vector<Isite>::iterator p = growthSites_.begin();
    while ((p != growthSites_.end()) && (!found)) {
        if (p->getId() == loc->getId()) {
            growthSites_.erase(p);
            found = true;
        }
        p++;
    }

    return found;
}

bool Interface::removeDissolutionSite (Site *loc)
{
     #ifdef DEBUG
        cout << "Interface::removeDissolutionSite " << loc->getId()
             << ", size is " << dissolutionSites_.size()
             << endl;
        cout.flush();
        bool found = false;
        cout << "Interface::removeDissolutionSite Trying to "
             << "declare iterator to Isite vector ...";
        cout.flush();
        vector<Isite>::iterator p;
        cout << "Interface::removeDissolutionSite Trying to "
             << "set it to beginning of dissolutionSites_ ... ";
        cout.flush();
        p = dissolutionSites_.begin();
        for (int i = dissolutionSites_.size() - 1; (i >= 0 && (!found)); i--) {
            if (dissolutionSites_[i].getId() == loc->getId()) {
                p += i;
                dissolutionSites_.erase(p);
                found = true;
            }
        }
    #else
        bool found = false;
        vector<Isite>::iterator p;
        p = dissolutionSites_.begin();
        for (int i = dissolutionSites_.size() - 1; (i >= 0 && (!found)); i--) {
            if (dissolutionSites_[i].getId() == loc->getId()) {
                p += i;
                dissolutionSites_.erase(p);
                found = true;
            }
        }
    #endif
    return found;
}
