/**
@file Isite.cc
@brief Method definitions for the Isite class.

@remark These could all be placed in the Isite.h file as inline functions.
@todo Place in the Isite.h file as inline functions, and delete this file.
*/
#include "Isite.h"

Isite::Isite ()
{
    affinity_ = 0;
    id_ = 0;
    verbose_ = false;
}

Isite::Isite (unsigned long int idval,
              int aftyval,
              const bool verbose)
{
    id_ = idval;
    affinity_ = aftyval;
    verbose_ = verbose;
}

Isite::Isite (const Isite &obj)
{
    id_ = obj.id_;
    affinity_ = obj.affinity_;
    verbose_ = obj.verbose_;
}
