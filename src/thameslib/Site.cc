/**
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

Site::Site(unsigned int xp, unsigned int yp, unsigned int zp, unsigned int xs,
           unsigned int ys, unsigned int zs, unsigned int neigh,
           ChemicalSystem *csys, const bool verbose) {
  x_ = y_ = z_ = 0;
  id_ = 0;
  nb_.clear();
  stressFreeVolume_ = trueVolume_ = 1.0;
  damage_ = false;
  x_ = xp;
  y_ = yp;
  z_ = zp;

#ifdef DEBUG
  verbose_ = true;
#else
  verbose_ = verbose;
#endif

  dissolution_ = -1;
  growth_.clear();

  id_ = (unsigned int)(x_ + (xs * y_) + ((xs * ys) * z_));

  if (neigh > 0)
    nb_.resize(neigh, 0);

  chemSys_ = csys;

  expstrain_ = 0.0;

  inGrowInterface_.clear();
  inDissInterface_.clear();
  inGrowInterface_.resize(chemSys_->getNumMicroPhases(),-1);
  inDissInterface_.resize(chemSys_->getNumMicroPhases(),-1);
}

void Site::calcWmc(void) {
  wmc_ = chemSys_->getMicroPhasePorosity(getMicroPhaseId());
  for (unsigned int i = 0; i < nb_.size(); i++) {
    wmc_ += chemSys_->getMicroPhasePorosity(nb_[i]->getMicroPhaseId());
  }
}
