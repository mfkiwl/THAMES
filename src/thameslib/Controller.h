/**
@file Controller.h
@brief Declaration of the Controller class.
*/

#ifndef CONTROLLERH
#define CONTROLLERH

#include "KineticController.h"
#include "Lattice.h"
#include "Site.h"
#include "ThermalStrain.h"
#include "global.h"
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;

struct RestoreSite{
    //for each site in site_:
    unsigned int microPhaseId;   // The microstructure phase assignment
    unsigned int dissolution;    // phase that can dissolve at this site or -1
    vector<unsigned int> growth; // Vector of phases that can grow at this site
    double wmc;                  // total porosity ("surface curvature") at this site
    double wmc0;                 // this site internal porosity (its own contribution at wmc_ value)
    int visit;                   // reset to 0
};

struct RestoreInterface{
    //  from Interface
    unsigned int microPhaseId; /**< The phase id of the voxels at this interface */
    vector<Isite> growthSites; /**< The list of all sites eligible foradjacent growth */
    vector<Isite> dissolutionSites; /**< The list of sites eligible for self-dissolution */
    //    for each Isite:
    //      unsigned int id_; /**< The id of the corresponding Site */
    //      int affinity_;    /**< The affinity for growth of a phase at the site */
    //      bool verbose_;    /**< Flag for whether to produce verbose output */
    //      double prob_;     /**< The growth probability of a phase at this site (computed according the affinity) */
    //      double probIni_;
};

struct RestoreSystem {
    //from ChemicalSystem (in fact from KineticController):
    //vector<double> ICMoles;
    vector<double> DCMoles;
    //from Lattice:
    vector<int> count;
    vector<RestoreSite> site;     /**< 1D list of Site objects (site = voxel) */
    //from Interface
    vector<RestoreInterface> interface;
};

//from ChemicalSystem:
//  double *ICMoles_;     /**< List of number of moles of each IC in system */
//  double *DCMoles_;             /**< List of moles of each DC */
//  double *prevGEMPhaseMoles_; /**< List of moles of each phase in the system in the previous time step */
//  double *prevGEMPhaseMass_;  /**< List of mass of each phase in the system in the previous time step */
//  double *prevGEMPhaseVolume_; /**< List of volume of each phase in the system in the previous time step */
//
//from Lattice:
//  vector<Site> site_;     /**< 1D list of Site objects (site = voxel) */
//  for each site in site_:
//    unsigned int microPhaseId_;   // The microstructure phase assignment
//    unsigned int dissolution_;    // phase that can dissolve at this site or -1
//    vector<unsigned int> growth_; // Vector of phases that can grow at this site
//    double wmc_;                  // total porosity ("surface curvature") at this site
//    double wmc0_;                 // this site internal porosity (its own contribution at wmc_ value)
//    >>int visit_;<<               // reset to 0
//vector<Interface> interface_;     //
//  from Interface
//    microPhaseId_; /**< The phase id of the voxels at this interface */
//    vector<Isite> growthSites_; /**< The list of all sites eligible foradjacent growth */
//    vector<Isite> dissolutionSites_; /**< The list of sites eligible for self-dissolution */
//    for each Isite:
//      unsigned int id_; /**< The id of the corresponding Site */
//      int affinity_;    /**< The affinity for growth of a phase at the site */
//      bool verbose_;    /**< Flag for whether to produce verbose output */
//      double prob_;     /**< The growth probability of a phase at this site (computed according the affinity) */
//      double probIni_;
//vector<int> count_;               // recreate or restored


/**
@class Controller
@brief Controls the running of simulation iterations.

The `Controller` class is the hub for THAMES simulations.  It has pointers
to all the other major objects that are instantiated by THAMES, including
the KineticController, the ThermalStrain, the ChemicalSystem, and the Lattice
that are associated with the system.

The `Controller` object is also responsible for running each iteration
of the simulation and deciding which modules to run based on whether
hydration, leaching, or sulfate attack is desired.

In THAMES, the `Lattice` class can be thought of as the system.
By itself, it can identify what materials it contains, the properties of
the materials, the temperature, current age, degree of hydration, etc.
However, the `Lattice` class does not possess a <i>driver</i> to
control the microstructure development.  This latter
functionality is contained in the `Controller` class, which operates
directly on the `Lattice` to determine how to modify the lattice at
each time step specified in the phase input file.

THAMES keeps track of physical and chemical data
about individual material phases, including
specific gravity, internal porosity, composition, molar volume, etc.  All the
phases are stored in a phase database.

The ultimate objective of the `Controller` class is to cycle through
the phase input file data, and to modify the lattice accordingly at each
time step.
*/

class Controller {

protected:
  string jobroot_;   /**< Root name for all output files */
  Lattice *lattice_; /**< Pointer to microstructure lattice object */
  KineticController
      *kineticController_;    /**< Pointer to kinetic controller object */
  ThermalStrain *thermalstr_; /**< Pointer to the finite element model object */

  /**
  @brief Stores the moles of independent components dissolved.

  @warning This member appears to not be used and may be obsolete
  @todo Investigate whether to delete this member
  */
  vector<double> molesdissolved_;

  double imgfreq_;          /**< Frequency to output microstructure image */
  ChemicalSystem *chemSys_; /**< Pointer to `ChemicalSystem` object */
  vector<double> time_;     /**< List of simulation times for each iteration */
  vector<double> output_time_; /**< List of times to output image */
  double statfreq_;            /**< Frequency to output statistics */

  int sim_type_; /**< Hydration, leaching, or sulfate attack for now */

private:
  double sattack_time_; /**< Simulation time at which to begin sulfate attack,
                                in days */
  double leach_time_;   /**< Simulation time at which to begin leaching,
                                in days */
  int damagecount_;     /**< Number of pixels in the lattice that are damaged */

  bool verbose_; /**< Flag for verbose output */
  bool warning_; /**< Flag for warning output */

public:
  /**
  @brief The constructor.

  This is the only Controller constructor provided.  It requires that all the
  auxiliary objects be defined already, including

      - The lattice object
      - The kinetic model object
      - The chemical system object (interface between GEM and THAMES
      - The finite element model for tracking strain and stress

  @param msh is a pointer to the already-instantiated `Lattice` object
  @param kc is a pointer to the already-instantiated `KineticController`
  @param cs is a pointer to the already-instantiated `ChemicalSystem` object
  @param thmstr is a pointer to the already-instantiated `ThermalStrain` object
  @param simtype is the type of simulation to run
  @param parfilename is the name of the input parameter file
  @param jobname is the root name to give to all output files
  @param verbose is true if verbose output should be produced
  @param warning is true if warning output should be produced
  */
  Controller(Lattice *msh, KineticController *kc, ChemicalSystem *cs,
             ThermalStrain *thmstr, const int simtype,
             const string &parfilename, const string &jobname,
             const bool verbose, const bool warning);

  /**
  @brief Run a computational iteration.

  This method launches one computational iteration, which includes

      - Consulting the kinetic controller to determine the number of moles of
          each independent component (IC) to add or subtract from the system
      - Running the GEM thermodynamic calculation
      - Running the finite element code (optionally) to update stress and strain
  states
      - Updating the lattice to reflect the new microstructure

  @param statfilename is the name of the file to store phase statistics
  @param choice is an int flag to specify whether simulating hydration,
  leaching, or sulfate attack
  */
  void doCycle(const string &statfilename, int choice);

  /**
  @brief Calculate the state of the system (called by doCycle).

  This method calculates the change in state of the system during a cycle,
  including

      - Consulting the kinetic model to get IC moles dissolved or added
      - (Optionally) determining IC moles to add from an external sulfate
  solution
      - Launching a thermodynamic calculation
      - (Optionally) determining the AFt saturation index for crystallization
  pressure calculations
      - Updating the microstructure
      - Outputting the microstructure phase volume fractions and other data

  @param time is the simulation time [days]
  @param dt is the change in simulation time used by the kinetic model [days]
  @param isFirst is true if this is the first state calculation
  (initialization)
  */
  //void calculateState(double time, double dt, bool isFirst, int cyc);

  int calculateState(double time, double dt, bool isFirst, int cyc);

  /**
  @brief Parse the input XML file specifying Controller parameters to use.

  Controller parameters that need to be input are

      - Length of time to calculate [days]
      - Frequency to output microstructure images

  @param docname is the name of the XML input file containing the Controller
  parameters
  */
  void parseDoc(const string &docname);

  /**
  @brief Set the simulation time at which to begin sulfate attack simulation.

  @param sattacktime is the simulation time to begin sulfate attack [days]
  */
  void setSattack_time(const double sattacktime) {
    sattack_time_ = sattacktime;
  }

  /**
  @brief Get the simulation time at which to begin sulfate attack simulation.

  @return the simulation time to begin sulfate attack [days]
  */
  double getSattack_time(void) const { return sattack_time_; }

  /**
  @brief Set the simulation type

  @param simtype is the simulation type
  */
  void setSim_type(const double simtype) { sim_type_ = simtype; }

  /**
  @brief Get the simulation type.

  @return the simulation type
  */
  double getSim_type(void) const { return sim_type_; }

  /**
  @brief Set the verbose flag

  @param isverbose is true if verbose output should be produced
  */
  void setVerbose(const bool isverbose) {
    verbose_ = isverbose;
    return;
  }

  /**
  @brief Get the verbose flag

  @return the verbose flag
  */
  bool getVerbose() const { return verbose_; }

  /**
  @brief Set the warning flag

  @param iswarning is true if warning output should be produced
  */
  void setWarning(const bool iswarning) {
    warning_ = iswarning;
    return;
  }

  /**
  @brief Get the warning flag

  @return the warning flag
  */
  bool getWarning() const { return warning_; }

  void writeTxtOutputFiles (double time);
  void writeTxtOutputFiles_onlyICsDCs (double time);

}; // End of Controller class
#endif
