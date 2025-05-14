/**
@file Isite.cc
@brief Method definitions for the Isite class.

@remark These could all be placed in the Isite.h file as inline functions.
@todo Place in the Isite.h file as inline functions, and delete this file.
*/
#include "Isite.h"

Isite::Isite() {
  affinity_ = 0;
  id_ = 0;
  verbose_ = false;
  prob_ = 0.;
}

Isite::Isite(int idval, double aftyval, const bool verbose, double prbval) {
  id_ = idval;
  affinity_ = aftyval;
  prob_ = prbval;

#ifdef DEBUG
  verbose_ = true;
#else
  verbose_ = verbose;
#endif
}

Isite::Isite(const Isite &obj) {
  id_ = obj.id_;
  affinity_ = obj.affinity_;
  verbose_ = obj.verbose_;
  prob_ = obj.prob_;
}

Isite& Isite::operator=(const Isite &obj){     // copy assignment operator

    id_=obj.id_;
    affinity_ = obj.affinity_;
    verbose_ = obj.verbose_;
    prob_ = obj.prob_;

    return *this;
}
