/**
@file Interface.cc
@brief Definition of methods for the Interface class.

*/
#include "Interface.h"
#include "ChemicalSystem.h"
#include "RanGen.h"

bool cmp(const Site *s1, const Site *s2) { return s1->getWmc() < s2->getWmc(); }

bool affinitySort(const Isite s1, const Isite s2) {
  return s1.getAffinity() > s2.getAffinity();
}

Interface::Interface() {
  microPhaseId_ = 0;
  growthSites_.clear();
  dissolutionSites_.clear();
}

Interface::Interface(const bool verbose) {
  microPhaseId_ = 0;
  growthSites_.clear();
  dissolutionSites_.clear();


#ifdef DEBUG
  verbose_ = true;
#else
  verbose_ = verbose;
#endif
}

Interface::Interface(ChemicalSystem *csys, vector<Site *> gv,
                     vector<Site *> dv, unsigned int pid, const bool verbose) {
  unsigned int j;
  unsigned int i;
  int afty;

  vector<Site *>::iterator beginLocation, endLocation;
  vector<Isite>::iterator start, end;

#ifdef DEBUG
  verbose_ = true;
#else
  verbose_ = verbose;
#endif


  microPhaseId_ = pid;
  chemSys_ = csys;

  dissolutionSites_.clear();
  growthSites_.clear();

  ///
  /// Eliminate duplicate values from the growth site vector
  ///
  int gvsize = gv.size();
  int dvsize = dv.size();
  if (gvsize > 0) {
    sort(gv.begin(), gv.end());
    // i = 1;
    beginLocation = gv.begin();
    endLocation = unique(gv.begin(), gv.end());
    gv.erase(endLocation, gv.end());
  }

  ///
  /// Now do the same thing to the dissolution site vector
  ///

  if (dvsize > 0) {
    sort(dv.begin(), dv.end());
    // i = 1;
    beginLocation = dv.begin();
    endLocation = unique(dv.begin(), dv.end());
    dv.erase(endLocation, dv.end());
  }

  ///
  /// Now sort the growth sites according to the affinity
  ///

  for (j = 0; j < gv.size(); j++) {
    afty = 0;
    for (i = 0; i < NN_NNN; i++) {
      afty += chemSys_->getAffinity(pid, gv[j]->nb(i)->getMicroPhaseId());
    }

    ///
    /// Add to the list of Isites.  An Isite is an object consisting
    /// of a pointer to a site and an affinity value
    ///

    growthSites_.push_back(Isite(gv[j]->getId(), afty));
  }

  ///
  /// Now sort the dissolution sites according to the affinity
  ///

  try {
    dvsize = dv.size();
    for (j = 0; j < dvsize; j++) {
      afty = 0;
      for (i = 0; i < NN_NNN; i++) {
        afty += chemSys_->getAffinity(pid, dv[j]->nb(i)->getMicroPhaseId());
      }
      dissolutionSites_.push_back(Isite(dv[j]->getId(), afty));
    }
  } catch (EOBException e) {
    e.printException();
    exit(1);
  }

} // End of constructors

Interface::~Interface() {
  growthSites_.clear();
  dissolutionSites_.clear();
}

bool Interface::addGrowthSite_newInterface(int id, int aff) {
  Isite tisite(id, aff);
  growthSites_.push_back(tisite);
  return true;
}

bool Interface::addGrowthSite(Site *loc) {
  bool answer = false;
  bool found = false;
  unsigned int i;
  // vector<Isite>::iterator p, q, start, end;
  // start = growthSites_.begin();
  // end = growthSites_.end();

  /// See if site is already present
  int size = growthSites_.size();
  int id = loc->getId();
  for (i = 0; i < size; i++) {
    if (id == growthSites_[i].getId()) {
      found = true;
      break;
    }
  }

  /// Add the site only if it was not found already
  if (!found) {
    int afty = 0;
    for (i = 0; i < NN_NNN;
         i++) { // NN_NNN = NUM_NEAREST_NEIGHBORS + NUM_SECONDNEAREST_NEIGHBORS
      afty +=
          chemSys_->getAffinity(microPhaseId_, loc->nb(i)->getMicroPhaseId());
    }
    Isite tisite(id, afty);
    // q = lower_bound(start, end, tisite, affinitySort);
    // growthSites_.insert(q, tisite);
    growthSites_.push_back(tisite);
    answer = true;
  }
  return answer;
}

bool Interface::addDissolutionSite(Site *loc) {
  bool answer = false;
  bool found = false;
  unsigned int i;
  unsigned int siteId = loc->getId();

  /// See if site is already present
  int size = dissolutionSites_.size();
  for (i = 0; i < size; i++) {
    if (siteId == dissolutionSites_[i].getId()) {
      found = true;
      break;
    }
  }

  /// Add the site only if it was not found already
  if (!found) {
    Isite tisite(siteId, 0);
    dissolutionSites_.push_back(tisite);
    answer = true;
  }
  return answer;
}

bool Interface::removeGrowthSite_0(Site *loc, int pos) {
  bool found = false;
  // int size = growthSites_.size();
  if (growthSites_[pos].getId() == loc->getId()) {
    growthSites_[pos] = growthSites_[growthSites_.size() - 1];
    growthSites_.pop_back();
    found = true;
  } else {
    cout << endl
         << "error site doesn't belong to this interface or wrong site "
            "position - phaseId/siteId/isitePos/growthSites_.size(): "
         << loc->getMicroPhaseId() << " / " << loc->getId() << " / " << pos
         << " / " << growthSites_.size() << endl;
    cout << "STOP:  Interface::removeGrowthSite_0(Site *loc, int pos)" << endl;
    exit(1);
  }
  return found;
}

bool Interface::removeGrowthSite_1(Site *loc) {
  bool found = false;
  int size = growthSites_.size();
  int siteId = loc->getId();

  for (int i = 0; i < size; i++) {
    if (growthSites_[i].getId() == siteId) {
      growthSites_[i] = growthSites_[size - 1];
      growthSites_.pop_back();
      found = true;
      break;
    }
  }
  if (found == false) {
    cout << endl;
    cout << "error site doesn't belong to this interface - "
            "phaseId/siteId/growthSites_.size(): "
         << loc->getMicroPhaseId() << " / " << loc->getId() << " / "
         << growthSites_.size() << endl;
    cout << "STOP:  Interface::removeGrowthSite_1(Site *loc)" << endl;
    exit(1);
  }
  return found;
}

bool Interface::removeEmptiedSite(int siteID) {

  bool found = false;
  int size = growthSites_.size();

  for (int i = 0; i < size; i++) {
    if (growthSites_[i].getId() == siteID) {
      growthSites_[i] = growthSites_[size - 1];
      growthSites_.pop_back();
      found = true;
      break;
    }
  }
  return found;
}

bool Interface::removeDissolutionSite_diss(Site *loc, int pos) {
  bool found = false;
  // int size = dissolutionSites_.size();
  if (dissolutionSites_[pos].getId() == loc->getId()) {
    dissolutionSites_[pos] = dissolutionSites_[dissolutionSites_.size() - 1];
    dissolutionSites_.pop_back();
    found = true;
  } else {
    cout << endl
         << "error site doesn't belong to this interface or wrong site position"
            "- phaseId/siteId/isitePos/dissolutionSites_.size(): "
         << loc->getMicroPhaseId() << " / " << loc->getId() << " / " << pos
         << " / " << dissolutionSites_.size() << endl;
    cout << "STOP:  Interface::removeDissolutionSite_diss(Site *loc, int pos)"
         << endl;
    exit(1);
  }
  return found;
}

bool Interface::removeDissolutionSite_grow(Site *loc) {
  bool found = false;
  int size_ = dissolutionSites_.size();
  int id = loc->getId();

  for (int i = 0; i < size_; i++) {
    if (id == dissolutionSites_[i].getId()) {
      dissolutionSites_[i] = dissolutionSites_[size_ - 1];
      dissolutionSites_.pop_back();
      found = true;
      break;
    }
  }
  if (found == false) {
    cout << endl
         << "error site doesn't belong to this interface - "
            "phaseId/siteId/growthSites_.size(): "
         << loc->getMicroPhaseId() << " / " << loc->getId() << " / "
         << growthSites_.size() << endl;
    cout << "STOP:  Interface::removeDissolutionSite_grow(Site *loc)" << endl;
    exit(1);
  }
  return found;
}
