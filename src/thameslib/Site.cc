/*
@file Site.cc
@brief Method definitions for the Site class.

*/

#include "Site.h"
#include "ChemicalSystem.h"

Site::Site() {
  x_ = y_ = z_ = 0;
  id_ = 0;
  nb_.clear();
  stressFreeVolume_ = trueVolume_ = 0.0;
  damage_ = false;
  expstrain_ = 0.0;
}

Site::Site(int xp, int yp, int zp, int xs, int ys, int zs, int neigh,
           ChemicalSystem *csys, const bool verbose) {
  x_ = y_ = z_ = 0;
  id_ = 0;
  nb_.clear();
  stressFreeVolume_ = trueVolume_ = 1.0;
  damage_ = false;
  x_ = xp;
  y_ = yp;
  z_ = zp;
  visit_ = 0;

#ifdef DEBUG
  verbose_ = true;
#else
  verbose_ = verbose;
#endif

  growth_.clear();

  id_ = x_ + (xs * y_) + ((xs * ys) * z_);

  if (neigh > 0)
    nb_.resize(neigh, 0);

  chemSys_ = csys;

  expstrain_ = 0.0;

  inGrowInterfacePos_.clear();
  inGrowInterfacePos_.resize(chemSys_->getNumMicroPhases(), -1);
  inDissInterfacePos_ = -1;

  inGrowthVectorPos_.clear();
  inGrowthVectorPos_.resize(chemSys_->getNumMicroPhases(), -1);
  inDissolutionVectorPos_ = -1;
}

void Site::calcWmc(void) {
  wmc_ = chemSys_->getMicroPhasePorosity(getMicroPhaseId());
  for (int i = 0; i < NN_NNN; i++) {
    wmc_ += chemSys_->getMicroPhasePorosity(nb_[i]->getMicroPhaseId());
  }
}
