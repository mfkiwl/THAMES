//-------------------------------------------------------------------
// $Id$
/// \file ms_param.cpp
/// Implementation  of default settings for the Interior Points Method
/// (IPM) module for convex programming Gibbs energy minimization
/// as well as I/O functions for IPM files
//
// Copyright (c) 1992,2012 S.Dmytriyeva, D.Kulik, K.Chudnenko
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------
//
#ifdef IPMGEMPLUGIN

#ifdef __unix__
#include <unistd.h>
#endif

#include <math.h>
#include "m_param.h"
#include "num_methods.h"
#include "gdatastream.h"
#include "node.h"

//TProfil* TProfil::pm;

const double R_CONSTANT = 8.31451,
              NA_CONSTANT = 6.0221367e23,
                F_CONSTANT = 96485.309,
                  e_CONSTANT = 1.60217733e-19,
                    k_CONSTANT = 1.380658e-23,
// Conversion factors
                      cal_to_J = 4.184,
                        C_to_K = 273.15,
                          lg_to_ln = 2.302585093,
                            ln_to_lg = 0.434294481,
                             H2O_mol_to_kg = 55.50837344,
                               Min_phys_amount = 1.66e-24;

enum volume_code {  // Codes of volume parameter
    VOL_UNDEF, VOL_CALC, VOL_CONSTR
};

SPP_SETTING pa_ = {
    " GEMS-GUI v.3.1 r.2150 (rc) " " GEMS3K v.3.1 r.690 (rc) ",
    {    // Typical default set (03.04.2012) new PSSC( logSI ) & uDD()
         2,  /* PC */  2,     /* PD */   -5,   /* PRD */
         1,  /* PSM  */ 130,  /* DP */   1,   /* DW */
         0, /* DT */     30000,   /* PLLG */   1,  /* PE */  7000, /* IIM */
         1000., /* DG */   1e-13,  /* DHB */  1e-20,  /* DS */
         1e-6,  /* DK */  0.01,  /* DF */  0.01,  /* DFM */
         1e-5,  /* DFYw */  1e-5,  /* DFYaq */    1e-5,  /* DFYid */
         1e-5,  /* DFYr,*/  1e-5,  /* DFYh,*/   1e-5,  /* DFYc,*/
         1e-6, /* DFYs, */  1e-17,  /* DB */   1.,   /* AG */
         0.,   /* DGC */   1.0,   /* GAR */  1000., /* GAH */
         1e-3, /* GAS */   12.05,  /* DNS */   1e-13,  /* XwMin, */
         1e-13,  /* ScMin, */  1e-33, /* DcMin, */   1e-20, /* PhMin, */
         1e-5,  /* ICmin */   1e-10,  /* EPS */   1e-3,  /* IEPS */
         1e-10,  /* DKIN  */ 0,  /* tprn */
    },
}; // SPP_SETTING



void BASE_PARAM::write(fstream& oss)
{
  short arr[10];

  arr[0] = PC;
  arr[1] = PD;
  arr[2] = PRD;
  arr[3] = PSM;
  arr[4] = DP;
  arr[5] = DW;
  arr[6] = DT;
  arr[7] = PLLG;
  arr[8] = PE;
  arr[9] = IIM;

  oss.write( (char*)arr, 10*sizeof(short) );
  oss.write( (char*)&DG, 28*sizeof(double) );
  oss.write( (char*)&tprn, sizeof(char*) );
}

void SPP_SETTING::write(fstream& oss)
{
    oss.write( ver, TDBVERSION );
    p.write( oss );
}

TProfil::TProfil( TMulti* amulti )
{
    pa= pa_;
    multi = amulti;
    pmp = multi->GetPM();
}

/// GEM IPM calculation of equilibrium state in MULTI
double TProfil::ComputeEquilibriumState( /*long int& RefinLoops_,*/ long int& NumIterFIA_, long int& NumIterIPM_ )
{
   //long int RefineLoops = RefinLoops_;
   //RefineLoops = 0L; // Provisional
   return multi->CalculateEquilibriumState( /*RefineLoops,*/ NumIterFIA_, NumIterIPM_ );
}

#endif

// ------------------ End of ms_param.cpp -----------------------




