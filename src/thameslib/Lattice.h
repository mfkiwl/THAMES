/**
@file Lattice.h
@brief Declaration of the Lattice class for storing the 3D microstructure

THAMES defines a Lattice class that is instantiated to a Lattice
object at the beginning of the program's execution.  The lattice defines the
three-dimensional environment within which a cement paste microstructure
exists, hydrates, and possibly deteriorates.
*/

#ifndef SRC_THAMESLIB_LATTICE_H_
#define SRC_THAMESLIB_LATTICE_H_

#include <algorithm>
#include <climits>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

#include "../version.h"
#include "AppliedStrain.h"
#include "ChemicalSystem.h"
#include "Interface.h"
#include "Isite.h"
#include "RanGen.h"
#include "Site.h"
#include "global.h"
#include "utils.h"

/**
@struct Sitesize
@brief Structure to catalog site domain sizes
*/
struct Sitesize {
  int siteid; /**< ID of the site in the site_ vector */
  int nsize;  /**< Size of the domain of the phase at that site */
};

struct chemElem {
  int z;
  string symb;
  double mass;
};

struct structGrowVect {
  int id;
  int posVect;
  double affinity;
};

struct structDissVect {
  int id;
  int posVect;
  double wmc;
};

using namespace std;

/**
@class Lattice
@brief Defines and stores the 3D microstructure as a discrete lattice of voxel
sites.

*/
class Lattice {

private:
  string version_; /**< THAMES version for header information */
  string thamesVersion_;
  string jobRoot_; /**< The root name for output files */
  string damageJobRoot_;

  RanGen *rg_; /**< Pointer to random number generator object */
  int latticeRNGseed_;
  long int numRNGcall_0_, numRNGcallLONGMAX_;
  double lastRNG_;

  int xdim_;          /**< Number of sites in the x dimension */
  int ydim_;          /**< Number of sites in the y dimension */
  int zdim_;          /**< Number of sites in the z dimension */
  double resolution_; /**< Voxel edge length [micrometers] */
  vector<Site> site_; /**< 1D list of Site objects (site = voxel) */
  int numSites_;      /**< Total number of sites */
  // unsigned int siteNeighbors_; /**< Number of neighbor sites to a given site
  // */
  ChemicalSystem *chemSys_; /**< Pointer to simulation's ChemicalSystem */
  AppliedStrain *FEsolver_; /**< Pointer to simulation's FE elastic solver */
  vector<Interface> interface_; /**< List of the different interface objects
                                        in the microstructure */

  double areaPerFace_;            /**< Converts a voxel face to m2 units */
  double volumePerVoxel_;         /**< Converts a voxel to its volume in m3
                                   units */
  double wsRatio_;                /**< Water-to-solids mass ratio */
  vector<double> volumeFraction_; /**< Array of volume fractions of each
                                          microstructure phase */
  vector<double> surfaceArea_;    /**< Array of surface areas of each
                                             microstructure phase
                                             (m2 per 100 g of all solid) */
  vector<double>
      specificSurfaceArea_; /**< Array of specific surface areas of each
                               microstructure phase
                               (m2 per kg of that phase) */
  vector<int> count_;       /**< Number of sites of each different type */

  map<int, vector<double>>
      expansion_; /**< Map of expansion strain of each voxel */
  map<int, vector<int>> expansion_coordin_; /**< Map of coordinates of sites
                                               with local expansion strain */
  double waterChange_;          /**< How much water must be added or subtracted
                                        due to hydration or deterioration */
  double microstructureVolume_; /**< Microstructure volume in GEM
volume units */
  double initialMicrostructureVolume_; /**< Initial microstructure volume in GEM
volume units */
  double capillaryPoreVolume_;         /**< Total volume of capillary pores */
  double capillaryPoreVolumeFraction_; /**< Total volume fraction of capillary
                                          pores */
  double subvoxelPoreVolume_;          /**< Total volume of subvoxel pores */
  double nonSolidVolume_;              /**< Total volume not solid */
  double solidVolumeWithPores_;        /** Total solid volume including their
internal pore volume */
  double waterVolume_;                 /** volume of electrolyte in GEM
volume units */
  double voidVolume_;                  /** volume of void in GEM volume
units */
  double capillaryWaterVolume_;        /**< Volume of capillary pore water */
  double capillaryVoidVolume_;         /**< Volume of capillary void space
(no water) */
  double subvoxelWaterVolume_;         /**< Volume of water in subvoxel
pores in GEM units */
  double subvoxelPoreVolumeFraction_;  /**< Total volume fraction of subvoxel
                                          pores */

  vector<struct PoreSizeVolume>
      masterPoreVolume_; /**< Pore size distribution and saturation */

  double time_;              /**< The current simulation time [h] */
  double temperature_;       /**< The current simulation temperature [K] */
  double oldtemp_;           /**< The temperature in the previous
                                     time step [K] */
  double sulfateAttackTime_; /**< Simulation time at which to begin
                                simulation of sulfate attack [h] */
  double leachTime_;         /**< Simulation time at which to begin
                                      simulation of leaching [h] */

  bool depthEffect_; /**< Whether or not PNG images should have
                             depth effect */
  bool verbose_;     /**< Flag to determine verbose output */
  bool warning_;     /**< Flag to determine warning message output */

  vector<chemElem> cfgElem_; /**< Holds periodic table information to output
                                files in cfg format */

  double initSolidMass_;

  double wcRatio_; /**< Water-to-cement mass ratio */

  int numMicroPhases_; /**< Number of microphases */

  double particRadius_; /**< used for graphical representation */

  vector<int>
      growthInterfaceSize_; /**< growth interface size of each microphase */
  vector<int> dissolutionInterfaceSize_; /**< dissolution interface size of each
                                            microphase */

  vector<int> growingVectSA_; /**for SULFATE ATTACK */
  int sizeGrowingVectSA_;
  vector<vector<int>> shrinking_;
  vector<vector<double>> volratios_;

  int waterDCId_; /**< coresp to DCName = "H2O@" */
  double waterMollarMass_;
  double waterMollarVol_;

  int DAMAGEID_;

public:
  /**
  @brief Constructor without input microstructure file name.

  This constructor simply initializes the dimensions and time to zero, sets
  the temperature to the globally defined reference temperature, and
  sets the lattice resolution to the globally defined reference value.

  @note Not currently used in THAMES

  @param cs is a pointer to the ChemicalSystem object for the simulation
  */
  Lattice(ChemicalSystem *cs);

  /**
  @brief Overloaded constructor with input microstructure file name.

  This constructor initializes the dimensions and time to zero, sets
  the temperature to the globally defined reference temperature, and
  sets the lattice resolution to the globally defined reference value.
  Afterward, the input microstructure file is opened and read, so that
  the voxel phase assignments can be made at each site.

  @param cs is a pointer to the ChemicalSystem object for the simulation
  @param rg is a pointer to the random number generator object
  @param seedRNG is the random number seed
  @param fileName is the name of the file containing the microstructure data
  @param verbose is true if extra messages are to be printed
  @param warning is true if warning messages are to be printed
  */
  Lattice(ChemicalSystem *cs, RanGen *rg, int seedRNG, const string &fileName,
          const bool verbose, const bool warning);

  /**
 @brief Destructor.

 This destructor clears out the `interface_` and `site_` vectors, and
 also deletes the allocated memory for the random number generator object,
 since this is the class that allocated the memory for that object.
 */
  ~Lattice();

  void setChemElem(string);
  chemElem getChemElem(int);

  /**
  @brief Set the number of sites in the x dimension.

  @param x is the number of sites in the x dimension
  */
  void setXDim(const int x) {
    xdim_ = x;
    numSites_ = (xdim_ * ydim_ * zdim_);
  }

  /**
  @brief Get the number of sites in the x dimension.

  @return the number of sites in the x dimension
  */
  int getXDim() const { return xdim_; }

  /**
  @brief Set the number of sites in the y dimension.

  @param y is the number of sites in the y dimension
  */
  void setYDim(const int y) {
    ydim_ = y;
    numSites_ = (xdim_ * ydim_ * zdim_);
  }

  /**
  @brief Get the number of sites in the y dimension.

  @return the number of sites in the y dimension
  */
  int getYDim() const { return ydim_; }

  /**
  @brief Set the number of sites in the z dimension.

  @param z is the number of sites in the z dimension
  */
  void setZDim(const int z) {
    zdim_ = z;
    numSites_ = (xdim_ * ydim_ * zdim_);
  }

  /**
  @brief Get the number of sites in the z dimension.

  @return the number of sites in the z dimension
  */
  int getZDim() const { return zdim_; }

  /**
  @brief Get the total number of lattice sites.

  The lattice is rectangular, so the total number of sites is
  `xdim_ * ydim_ * zdim_`, but we store this value as a class member to
  save having to compute it multiple times.

  @return the total number of lattice sites
  */
  int getNumSites() const { return numSites_; }

  /**
  @brief Set the volume fraction of a given microstructure phase.

  @param i is the index of the microstructure phase
  @param vfrac is the volume fraction to assign on a total microstructure basis
  */
  // void setVolumeFraction(const unsigned int i, const double vfrac) {
  //   if (i > -1 && i < volumeFraction_.size()) {
  //     volumeFraction_[i] = vfrac;
  //   } else {
  //     throw EOBException("Lattice", "setVolumeFraction", "volumeFraction_",
  //                        volumeFraction_.size(), i);
  //   }
  // }

  /**
  @brief Set the initial volume fraction of a given microstructure phase.

  @param i is the index of the microstructure phase
  @param vfrac is the volume fraction to assign on a total microstructure basis
  */
  // void setInitVolumeFraction(const int i, const double vfrac) {
  //   if (i > -1 && i < initVolumeFraction_.size()) {
  //     initVolumeFraction_[i] = vfrac;
  //   } else {
  //     throw EOBException("Lattice", "setInitVolumeFraction",
  //                        "initVolumeFraction_", initVolumeFraction_.size(),
  //                        i);
  //   }
  //   initVolumeFraction_[i] = vfrac;
  // }

  /**
  @brief Set the water-solids mass ratio

  @param ws is the water-solids mass ratio
  */
  void setWsRatio(const double ws) {
    wsRatio_ = 0.0;
    if (ws > 0.0) {
      wsRatio_ = ws;
    }
    return;
  }

  /**
  @brief Get the water-solids mass ratio

  @return the water-solids mass ratio
  */
  double getWsRatio(void) const { return wsRatio_; }
  double getWcRatio(void) const { return wcRatio_; }

  /**
  @brief Get the volume fraction of a given microstructure phase.

  This is simply the number of sites with a given phase divided by the
  total number of sites.

  @param i is the index of the microstructure phase
  @return the volume fraction of phase i on a total microstructure basis
  */
  double getVolumeFraction(int i) { return (volumeFraction_[i]); }

  /**
  @brief Calculate the subvoxel pore volume

  @param vol is the array of all microstructure phase volumes
  */
  void calcSubvoxelPoreVolume(vector<double> &vol);

  /**
  @brief Calculate the total volume of solids including
  subvoxel pore volume assigned to solids

  @param vol is the array of all microstructure phase volumes
  it
  */
  void calcSolidVolumeWithPores(vector<double> &vol);

  /**
  @brief Get the total volume of solids including
  subvoxel pore volume assigned to solids

  @return the solid volume including subvoxel pore volume
  */
  double getSolidVolumeWithPores(void) const { return solidVolumeWithPores_; }

  /**
  @brief Calculate the non-solid volume

  @param vol is the array of all microstructure phase volumes
  it
  */
  void calcNonSolidVolume(vector<double> &vol);

  /**
  @brief Get or calculate the non-solid volume

  @return the non-solid volume
  */
  double getNonSolidVolume(void) const { return nonSolidVolume_; }

  /**
  @brief Get the number of neighbor sites each site has.

  This is simply the number of sites with a given phase divided by the
  total number of sites.

  @note NOT USED.

  @return the number of neighbor sites each site has
  */
  // unsigned int getSiteNeighbors() const { return siteNeighbors_; }

  /**
  @brief Set the lattice resolution [meters].

  The lattice resolution is the physical length associated with the edge
  length of a site.

  @param res is the lattice resolution [meters]
  */
  void setResolution(const double res);

  /**
  @brief Get the lattice resolution [meters].

  The lattice resolution is the physical length associated with the edge
  length of a site.

  @note NOT USED.

  @return the lattice resolution [meters]
  */
  double getResolution() const { return resolution_; }

  /**
  @brief Get the areaPerFace_ value [m2].

  @return the areaPerFace_ value [m2]
  */
  double getAreaPerFace() const { return areaPerFace_; }

  /**
  @brief Get the volumePerVoxel_ value [m3].

  @return the volumePerVoxel_ value [m3]
  */
  double getVolumePerVoxel() const { return volumePerVoxel_; }

  /**
  @brief Set the simulation time [hours].

  @note NOT USED.

  @param tval is the simulation time [hours]
  */
  void setTime(const double tval) { time_ = tval; }

  /**
  @brief Get the simulation time [hours].

  @note NOT USED.

  @return the simulation time [hours]
  */
  double getTime() const { return time_; }

  /**
  @brief Get the simulation time at which to start sulfate attack simulation
  [hours].

  @note NOT USED.

  @return the simulation time at which to start sulfate attack [hours]
  */
  double getSulfateAttackTime() const { return sulfateAttackTime_; }

  /**
  @brief Set the simulation time at which to start sulfate attack simulation
  [hours].

  @param sattacktime is the simulation time at which to start sulfate attack
  [hours]
  */
  void setSulfateAttackTime(const double sattacktime) {
    sulfateAttackTime_ = sattacktime;
  }

  /**
  @brief Get the simulation time at which to start leaching simulation [hours].

  @note NOT USED.

  @return the simulation time at which to start leaching [hours]
  */
  double getLeachTime() const { return leachTime_; }

  /**
  @brief Set the simulation time at which to start leaching simulation [hours].

  @param leachtime is the simulation time at which to start leaching [hours]
  */
  void setLeachTime(const double leachtime) { leachTime_ = leachtime; }

  /**
  @brief Set the lattice temperature [K].

  @param tmp is the temperature [K]
  */
  void setTemperature(const double tmp) { temperature_ = tmp; }

  /**
  @brief Get the lattice temperature [K].

  @return the temperature [K]
  */
  double getTemperature() const { return temperature_; }

  /**
  @brief Get the version of THAMES

  @note NOT USED.

  @return the version number as a string
  */
  const string &getVersion() const { return version_; }

  /**
  @brief Set the root name for simulation output files.

  @param jobname is the root name for simulation output files
  */
  void setJobRoot(string jobname) {
    jobRoot_ = jobname;
    // damageJobRoot_ = jobRoot_ + ".damage";
  }

  /**
  @brief Add a site at location (xp,yp,zp) to the lattice.

  The site is checked for valid coordinates.  If valid a new Site object
  is created and pushed back onto the class's `site_` vector.

  @param xp is the x coordinate of the site to add
  @param yp is the y coordinate of the site to add
  @param zp is the z coordinate of the site to add
  */
  void addSite(const int xp, const int yp, const int zp);

  /**
  @brief Get the x coordinate of a site with a given index in the 1D `site_`
  array.

  @param i is the index of the site in the class's `site_` array
  @return the x coordinate
  */
  int getX(const int i) const { return (site_[i].getX()); }

  /**
  @brief Get the y coordinate of a site with a given index in the 1D `site_`
  array.

  @param i is the index of the site in the class's `site_` array
  @return the y coordinate
  */
  int getY(const int i) const { return (site_[i].getY()); }

  /**
  @brief Get the z coordinate of a site with a given index in the 1D `site_`
  array.

  @param i is the index of the site in the class's `site_` array
  @return the x coordinate
  */
  int getZ(const int i) const { return (site_[i].getZ()); }

  /**
  @brief Get a site's index in the 1D `site_` array, given its (x,y,z)
  coordinates.

  @param ix is the x coordinate of the site
  @param iy is the x coordinate of the site
  @param iz is the x coordinate of the site
  @return the index of the site in the `site_` array
  */
  int getIndex(int ix, int iy, int iz) const;

  /**
  @brief Get the collection of site indices neighboring a given site.

  @param sitenum is the index of the site in question
  @param size is the maximum distance defining the neighborhood [sites]
  @return a list of site indices for all neighbors within the maximum distance
  */
  vector<int> getNeighborhood(const int sitenum, const int size);

  /**
  @brief Get a pointer to a Site object at a given index in the `site_` array.

  @param index is the index of the Site object in the `site_` array
  @return a pointer to the Site object in question
  */
  Site *getSite(int index) { return &site_[index]; }

  /**
  @brief Designate a site as damaged.

  The site to be damaged is specified by its index in the `site_` array.

  @note NOT USED.

  @param index is the index of the Site object in the `site_` array
  */
  void setDamage(int index) { site_[index].setDamage(); }

  /**
  @brief Change the wmc (weighted mean curvature) of a site by a prescribed
  amount.

  @param index is the index of the Site object in the `site_` array
  @param dwmcval is the increment to add to the wmc
  */
  void dWmc(int index, double dwmcval) {
    site_[index].setWmc(site_[index].getWmc() + dwmcval);
  }

  void setWmc0(int index, double dwmcval) { site_[index].setWmc0(dwmcval); }

  /**
  @brief Compute normalized initial microstructure phase masses

  Given the initial masses of all phases in 1 cm3 of microstructure,
  this method scales them to 100 grams of solid instead.  In the process,
  this method also sets the initial moles of water in the
  chemical system definition.

  So when this method is finished, microPhaseMass has units of g per 100 g of
  solid

  @param microPhaseMass is a vector of all the microstructure masses
  @param cementMass is the combined mass of all the cement components
  @param solidMass is the combined mass of all the solids
  */
  void normalizePhaseMasses(vector<double> microPhaseMass);

  /**
  @brief Master method to locate the interfaces for each phase in the
  microstructure.

  */
  void findInterfaces(void);

  /**
  @brief Add (grow i.e. switch from electrolyte) the prescribed number of
  sites of each microphase that has to grow

  @param growPhaseIDVect is the vector of microphase IDs that must grow
  @param numSiteGrowVect is a vector containing the number of voxels to add for
  each microphase ID in growPhaseIDVect
  @param growPhNameVect is a vector containing the name of each microphase in
  growPhaseIDVect
  @param numtoadd_G is the number of sites switched by this call
  @param totalTRC is the total call number of the changeMicrostructure method
  @return the actual number of sites that were changed for each microphase ID
  from the input growPhaseIDVect vector
  */
  vector<int> growPhase(vector<int> growPhaseIDVect,
                        vector<int> numSiteGrowVect,
                        vector<string> growPhNameVect, int &numadded_G,
                        int totalTRC);

  /**
  @brief create a new growth interface for a given phase (phaseID) having a
  size of numLeft sites; this is necessary when the "growth" of all requested
  sites for this phase was not possible because the size of the
  corresponding growth interface was zero.

  @param phaseid is the id of the microstructure phase to nucleate
  @param numLeft is the number of sites to nucleate/create for this phase
  @return the actual size of the new interface (must equals numLeft!)
  */
  void nucleatePhaseAff(int phaseID, int numLeft);
  void nucleatePhaseRnd(int phaseID, int numLeft);

  /**
  @brief Remove (dissolve i.e. switch to electrolyte) the prescribed number of
  sites of each microphase that has to dissolve

  @param dissPhaseIDVect is the vector of microphase IDs that must dissolve
  @param numSiteDissVect is a vector containing the number of voxels to dissolve
  for each microphase ID in dissPhaseIDVect
  @param dissPhNameVect is a vector containing the name of each microphase in
  dissPhaseIDVect
  @param numtoadd_D is the number of sites switched by this call
  @param totalTRC is the total call number of the changeMicrostructure method
  @return vector of the number of voxels of each phase that could not dissolve
  */
  vector<int> dissolvePhase(vector<int> dissPhaseIDVect,
                            vector<int> numSiteDissVect,
                            vector<string> dissPhNameVect, int &numadded_D,
                            int totalTRC);

  /**
  @brief Remove the water from a prescribed number of solution-filled sites.

  This method constructs a list of all the <i>potential</i> void sites, based
  on whether there are multiple connected solution-filled sites in a cluster.
  The list is then sorted essentially by the effective pore size.  Only then
  is the list visited and the prescribed number of sites switched to void.

  @param numsites is the number of sites to switch from water to void
  @return the actual number of sites that were changed
  */
  int emptyPorosity(int numsites, int cyc);

  /**
  @brief Add water to a prescribed number of empty pore sites.

  This method constructs a list of all the void sites, based
  on whether there are multiple connected void sites in a cluster.
  The list is then sorted essentially by the effective pore size.  Only then
  is the list visited and the prescribed number of sites switched to water.

  @param numsites is the number of sites to switch from void to water
  @return the actual number of sites that were changed
  */
  int fillPorosity(int numsites, int cyc);

  double fillAllPorosity(int cyc);

  /**
  @brief Count the number of solution sites within a box centered on a given
  site.

  @param boxsize is the linear dimension of the cubic box neighborhood
  @param siteid is the index of the queried site in the `site_` array
  @return the number of solution-filled sites in the box neighborhood
  */
  int countBox(int boxsize, unsigned int siteid);

  /**
  @brief Check whether a linear coordinate is outside the lattice boundaries.

  If a given coordinate is outside the lattice boundaries, then the additive
  adjustment is returned that will locate the equivalent site within the
  lattice, assuming periodic boundary conditions.

  @param pos is the linear coordinate to check
  @param size is the dimension of the lattice in that dimension (number of
  sites)
  @return the additive adjustment to locate the equivalent coordinate within
  the lattice
  */
  int checkBC(int pos, int size) {
    if (pos >= size)
      return (-size);
    if (pos < 0)
      return (size);
    return (0);
  }

  /**
  @brief Get a pointer to the ChemicalSystem object for the simulation.

  @return a pointer to the ChemicalSystem object for the simulation
  */
  ChemicalSystem *getChemSys() const { return chemSys_; }

  /**
  @brief Set the phase id of a given site, specified by a pointer to the Site
  object.

  @param s is a pointer to the Site object
  @param i is the phase index to set at that site
  */
  void setMicroPhaseId(Site *s, const int i) {
    count_[s->getMicroPhaseId()]--;
    s->setMicroPhaseId(i);
    count_[i]++;
  }

  /**
  @brief Set the phase id of a given site, specified by the site's index number.

  @param sitenum is the index of the site in the `site_` array
  @param i is the phase index to set at that site
  */
  void setMicroPhaseId(const int sitenum, const int i) {
    count_[site_[sitenum].getMicroPhaseId()]--;
    site_[sitenum].setMicroPhaseId(i);
    count_[i]++;
  }

  /**
  @brief Get the phase id of a given site, specified by the site's index number.

  @param sitenum is the index of the site in the `site_` array
  @return the microstructure phase id at the site
  */
  int getMicroPhaseId(const int sitenum) {
    return (site_[sitenum].getMicroPhaseId());
  }

  /**
  @brisef Add a site to the list of sites where dissolution of a given phase can
  occur.

  @param loc is a pointer to the Site object to add to the list of potential
  dissolution sites
  @param pid is the microstructure phase id
  */
  void addDissolutionSite(Site *loc, int pid);

  /**
  @brief Add a site to the list of sites where growth of a given phase can
  occur.

  @param loc is a pointer to the Site object to add to the list of potential
  growth sites
  @param pid is the microstructure phase id
  */
  void addGrowthSite(Site *loc, int pid);

  /**
  @brief Remove a site from the list of sites where dissolution of a given phase
  can occur.

  @param loc is a pointer to the Site object to remove from the list of
  potential dissolution sites
  @param pid is the microstructure phase id
  */
  void removeDissolutionSite(Site *loc, int pid);

  /**
  @brief Remove a site from the list of sites where growth of a given phase can
  occur.

  @param loc is a pointer to the Site object to remove from the list of
  potential growth sites
  @param pid is the microstructure phase id
  */
  void removeGrowthSite_diss(Site *loc, int pid);

  void removeGrowthSite_grow(Site *ste0, int pid);

  void removeGrowthSite_nucleation(Site *loc);

  /**
  @brief Master method to update a microstructure during after a given time
  interval.

  Updating the microstructure includes determining how many sites of each phase
  to add and subtract from the lattice, determining which site locations will be
  used to do that, and then actually causing the switches in phase id to happen
  at those sites. The interfaces and lists of dissolution and growth sites are
  updated accordingly, too.

  @note Water is assumed to be chemically reactive only if it is in capillary
  porosity (microstructure id ELECTROLYTEID).  If the capillary water is
  exhausted then some reaction can still happen with water in nanoporosity, but
  for now we assume that the nanopore water is chemically unreactive and cannot
  be removed.

  @todo Generalize to allow water in nanopores to be chemically reactive

  @param time is is the simulation time [hours]
  @param simtype is the type of simulation (hydration, leaching, etc)
  @param capWater is true if there is any capillary pore water in the system.
  @param vectPhNumDiff is the vector of maximum number of voxels belonging to
  each microphase, voxels that can be dissolved according to the system
  configuration (lattice)
  @param vectPhIdDiff is the microphase ID for which a the number of voxels that
  can be dissolved is smaller than the number requested by the corresponding
  kinetic model
  @param vectPhNameDiff is the vector of names of microphases
  @param repeat counts the number of changeMicrostructure calls for a given
  cycle (cyc)
  @param cyc (cycle) is the iteration number in main iteration loop in
  Controller::doCycle - each cycle corresponds to a time step

  @return zero if okay or nonzero if not all requested voxels
  for a certain microphase ID (phDiff) can be dissolved
  */
  int changeMicrostructure(double time, const int simtype, bool &capWater,
                           vector<int> &vectPhNumDiff,
                           vector<int> &vectPhIdDiff,
                           vector<string> &vectPhNameDiff, int repeat, int cyc);

  /**
  @brief Adjust GEMS calculated volumes of microstructure phases

  The volume fractions passed to this function are those coming directly
  from the chemical system.  But the chemical system does not account for
  occluded porosity that may be associated with a solid phase at length
  scales smaller than the lattice spatial resolution.  This method fixes
  those volume fractions, paying special attention to the water distribution.

  @note Water is assumed to be chemically reactive only if it is in capillary
  porosity (microstructure id ELECTROLYTEID).  If the capillary water is
  exhausted then some reaction can still happen with water in nanoporosity, but
  for now we assume that the nanopore water is chemically unreactive and cannot
  be removed.

  @todo Generalize to allow water in nanopores to be chemically reactive

  @param phasenames is a vector of the microstructure phase names
  @param vol is a vector of the pre-adjusted microstructure volumes
  */
  void adjustMicrostructureVolumes(vector<double> &vol, int volSize, int cyc);

  /**
  @brief Calculate microstructure volume fractions

  @param names is a vector of the adjusted microstructure volumes
  @param vol is a vector of the adjusted microstructure volumes
  @param vfrac will hold the microstructure volume fractions
  */
  void adjustMicrostructureVolFracs(vector<string> &names,
                                    const vector<double> vol,
                                    vector<double> &vfrac, int volSize,
                                    int cyc);

  /**
  @brief Calculate the pore size distribution data

  */
  void calculatePoreSizeDistribution(void);

  /**
  @brief Write the pore size distribution data to a file

  @param curtime is the current time in hours
  @param resolvedtime is the current time resolved into y,d,h,m
  @param simtype is the sumulation tyupe
  @param root is the root name of the output file to create
  */
  void writePoreSizeDistribution(const double curtime,
                                 const TimeStruct resolvedTime);

  /**
  @brief Write the microstructure colors to a file

  This is done to save processing the chemistry.json file just to get the colors
  and will make post-processing of images easier.

  @param root is the root name of the output file to create
  */
  void writeMicroColors();

  /**
  @brief Write the 3D microstructure to a file.

  The microstructure output file will indicate the phase id at each site.

  @param curtime is the current time in hours
  @param resolvedtime is the current time resolved into y,d,h,m
  @param root is the root name of the output file to create
  */
  void writeLattice(const double curtime, const TimeStruct resolvedtime);

  void writeLatticeH(const double curtime, const TimeStruct resolvedtime);

  void writeLatticeXYZ(const double curtime, const TimeStruct resolvedtime);

  void appendXYZ(double curtime);

  void writeLatticeCFG(const double curtime, const TimeStruct resolvedtime);

  void writeNewLattice(int newZdim);

  /**
  @brief Write the 3D microstructure to a file.

  The damage output file is binary, each site either being damaged or not.

  @param curtime is the current time in hours
  @param resolvedtime is the current time resolved into y,d,h,m
  @param root is the root name of the output file to create
  */
  void writeDamageLattice(const double curtime, const TimeStruct resolvedtime);

  /**
  @brief Write the 3D microstructure to a png file that can be immediately
  rendered.

  @param curtime is the current time in hours
  @param resolvedtime is the current time resolved into y,d,h,m
  @param root is the root name of the png output file to create
  */
  void writeLatticePNG(const double curtime, const TimeStruct resolvedtime);

  /**
  @brief Write the 3D microstructure to a png file that can be immediately
  rendered.

  The damage output file is binary, each site either being damaged or not.

  @param curtime is the current time in hours
  @param resolvedtime is the current time resolved into y,d,h,m
  @param root is the root name of the png output file to create
  */
  void writeDamageLatticePNG(const double curtime,
                             const TimeStruct resolvedtime);

  /**
  @brief Create files of sequential slices of the microstructure in the x
  direction.

  The slices are individual PPM files of 2D (y,z) microstructure slices,
  written back to back, in the same file.  Once created, the files are each
  converted to GIFs using a system call to Imagemagick, and then the GIFs are
  converted to an animated GIF file using a system call to gifsicle.

  @note NOT USED.

  @warning This method currently depends on system calls
  @warning This method currently depends on having Imagemagick installed
  @warning This method currently depends on having gifsicle installed

  @todo Remove the dependence on system calls, Imagemagick, and gifsicle

  @param root is the root name of the png output file to create
  */
  void makeMovie();

  /**
  @brief Set the expansion strain components of a site specified by its index.

  This function changes the strain components of a site already in the
  list of expansion sites.  If the prescribed site is not already in the
  list of expansion sites, then the site will be added to that list.

  @param index is the index of the site in the `site_` array
  @param val is the vector of expansion strain components to set
  */
  void setExpansion(int index, vector<double> val) {
    map<int, vector<double>>::iterator p = expansion_.find(index);
    if (p != expansion_.end()) {
      p->second = val;
    } else {
      expansion_.insert(make_pair(index, val));
    }
  }

  /**
  @brief Get the expansion strain components of a site specified by its index.

  @param index is the index of the site in the `site_` array
  @return the vector of expansion strain components to set
  */
  vector<double> getExpansion(int index) {
    map<int, vector<double>>::iterator p = expansion_.find(index);
    if (p != expansion_.end()) {
      return p->second;
    } else {
      string msg = "Could not find expansion_ match to index provided";
      throw EOBException("Lattice", "getExpansion", msg, expansion_.size(),
                         index);
    }
  }

  /**
  @brief Get the expansion strain components for all strained sites in the
  lattice.

  @return the map of the strain components, keyed to the site index numbers
  */
  map<int, vector<double>> getExpansion() { return expansion_; }

  /**
  @brief Get the coordinates of local region for calculating expansion stress.

  This gets the coordinates of the center site of a box in the lattice within
  which the expansion strain is calculated in the ThermalStrain model due to
  local crystallization pressure.

  @todo Change the function name to something like getExpansionSiteCoordinates.

  @param index is the index of a site that has crystallization pressure
  @return the (x,y,z) coordinates of the site
  */
  // vector<int> getExpansionCoordin(int index) {
  //   map<int, vector<int>>::iterator p = expansion_coordin_.find(index);
  //   if (p != expansion_coordin_.end()) {
  //     return p->second;
  //   } else {
  //     string msg = "Could not find expansion_coordin_ match to index
  //     provided"; throw EOBException("Lattice", "getExpansionCoordin", msg,
  //                        expansion_coordin_.size(), index);
  //   }
  // }

  /**
  @brief Set the coordinates of local site for calculating expansion stress.

  This gets the coordinates of the center site of a box in the lattice within
  which the expansion strain is calculated in the ThermalStrain model due to
  local crystallization pressure.

  @note NOT USED (commented in Controller)

  @todo Change the function name to something like setExpansionSiteCoordinates

  @param index is the index of a site that has crystallization pressure
  @param coordin is the (x,y,z) triple of the site's coordinates
  @return the (x,y,z) coordinates of the site
  */
  // void setExpansionCoordin(int index, vector<int> coordin) {
  //   map<int, vector<int>>::iterator p = expansion_coordin_.find(index);
  //   if (p == expansion_coordin_.end()) {
  //     expansion_coordin_.insert(make_pair(index, coordin));
  //   }
  // }

  /**
  @brief Get the microstructure volume

  @return the microstructure volume (GEMS volume units)
  */
  double getMicrostructureVolume(void) const {
    return (chemSys_->getMicroVolume());
  }

  /**
  @brief Get the initial microstructure volume

  @return the initial microstructure volume (GEMS volume units)
  */
  double getInitialMicrostructureVolume(void) const {
    return (chemSys_->getInitMicroVolume());
  }

  /**
  @brief Get the total capillary pore volume

  @return the volume of capillary pores (GEMS volume units)
  */
  double getCapillaryPoreVolume(void) const { return capillaryPoreVolume_; }

  /**
  @brief Set the capillary pore volume

  @param capillaryporevolume is the capillary pore volume (GEMS volume units)
  */
  void setCapillaryPoreVolume(double capillaryporevolume) {
    capillaryPoreVolume_ = capillaryporevolume;
  }

  /**
  @brief Get the total capillary pore volume fraction
  This is calculated on a total system volume basis

  @return the volume fraction of capillary pores (microstructure basis)
  */
  double getCapillaryPoreVolumeFraction(void) const {
    return capillaryPoreVolumeFraction_;
  }

  /**
  @brief Set the capillary pore volume fraction
  This is calculated on a total system volume basis

  @param capillaryPoreVolumeFraction is the capillary pore volume
  fraction (microstructure basis)
  */
  void
  setCapillaryPoreVolumeFraction(const double capillaryPoreVolumeFraction) {
    capillaryPoreVolumeFraction_ = capillaryPoreVolumeFraction;
  }

  /**
  @brief Get the total subvoxel pore volume

  @return the volume of subvoxel pores (GEMS volume units)
  */
  double getSubvoxelPoreVolume(void) const { return subvoxelPoreVolume_; }

  /**
  @brief Set the subvoxel pore volume

  @param subvoxelporevolume is the subvoxel pore volume (GEMS volume units)
  */
  void setSubvoxelPoreVolume(const double subvoxelporevolume) {
    subvoxelPoreVolume_ = subvoxelporevolume;
  }

  /**
  @brief Set the subvoxel pore volume

  @param subvoxelporevolume is the subvoxel pore volume (GEMS volume units)
  */
  void setNonSolidVolume(const double nonsolidvolume) {
    nonSolidVolume_ = nonsolidvolume;
  }

  /**
  @brief Get the capillary water volume

  @param vol is the volume of each microstructure phase
  */
  void calcCapillaryWaterVolume(vector<double> &vol);

  /**
  @brief Get the capillary water volume

  @return the capillary water volume
  */
  double getCapillaryWaterVolume(void) const { return capillaryWaterVolume_; }

  /**
  @brief Set the capillary water volume

  @param capillarywatervolume is the capillary water volume (GEMS volume units)
  */
  void setCapillaryWaterVolume(const double capillarywatervolume) {
    capillaryWaterVolume_ = capillarywatervolume;
  }

  /**
  @brief Get the capillary void volume

  @param vol is the volume of each microstructure phase
  @param calc is true only if calculating instead of just returning
  */
  // void calcCapillaryVoidVolume(void);

  /**
  @brief Get the capillary void volume

  @return the capillary void volume
  */
  double getCapillaryVoidVolume(void) const { return capillaryWaterVolume_; }

  /**
  @brief Set the capillary void volume

  @param capillaryvoidvolume is the capillary void volume (GEMS volume units)
  */
  void setCapillaryVoidVolume(const double capillaryvoidvolume) {
    capillaryVoidVolume_ = capillaryvoidvolume;
  }

  /**
  @brief Get the total subvoxel pore volume fraction
  This is calculated on a total system volume basis

  @return the volume fraction of subvoxel pores (microstructure basis)
  */
  double getSubvoxelPoreVolumeFraction(void) const {
    return subvoxelPoreVolumeFraction_;
  }

  /**
  @brief Set the subvoxel pore volume fraction
  This is calculated on a total system volume basis

  @param subvoxelporevolumefraction is the subvoxel pore volume
  fraction (microstructure basis)
  */
  void setSubvoxelPoreVolumeFraction(const double subvoxelporevolumefraction) {
    subvoxelPoreVolumeFraction_ = subvoxelporevolumefraction;
  }

  /**
  @brief Set the master pore volume distribution

  @param masterporevolume is the pore volume distribution
  */
  // void
  // setMasterPoreVolume(const vector<struct PoreSizeVolume> masterporevolume) {
  //   masterPoreVolume_ = masterporevolume;
  //   return;
  // }

  /**
  @brief Set the master pore volume distribution of a particular size

  @param idx is the index to set
  @param diam is the diameter in nm
  @param volume is the volume of pores this size, in nm3
  @param volfrac is the volume fraction of this size filled with electrolyte
  */
  /*
  void setMasterPoreVolume(const int idx, const double diam,
                           const double volume, const double volfrac) {
    try {
      if (idx >= masterPoreVolume_.size()) {
        throw EOBException("Lattice", "setMasterPoreVolume",
                           "masterPoreVolume_", masterPoreVolume_.size(),
                           (int)idx);
      }
      masterPoreVolume_[idx].diam = diam;
      masterPoreVolume_[idx].volume = volume;
      masterPoreVolume_[idx].volfrac = volfrac;
    } catch (EOBException ex) {
      ex.printException();
      exit(1);
    }
    return;
  }
  */

  /**
  @brief Get the master pore volume distribution of a particular size
  @param idx is the index to get
  @return the structure holding the pore size distribution data for that element
  */
  /*
  struct PoreSizeVolume getMasterPoreVolume(const int idx) {
    try {
      if (idx >= masterPoreVolume_.size()) {
        throw EOBException("Lattice", "getMasterPoreVolume",
                           "masterPoreVolume_", masterPoreVolume_.size(),
                           (int)idx);
      }
    } catch (EOBException ex) {
      ex.printException();
      exit(1);
    }
    return (masterPoreVolume_[idx]);
  }
  */

  /**
  @brief Get the diameter of the idx element of the pore volume distribution
  @param idx is the index to get
  @return the diameter of that element in the pore size distribution (nm)
  */
  /*
  double getMasterPoreVolumeDiam(const int idx) {
    try {
      if (idx >= masterPoreVolume_.size()) {
        throw EOBException("Lattice", "getMasterPoreVolumeDiam",
                           "masterPoreVolume_", masterPoreVolume_.size(),
                           (int)idx);
      }
    } catch (EOBException ex) {
      ex.printException();
      exit(1);
    }
    return (masterPoreVolume_[idx].diam);
  }
  */

  /**
  @brief Get the total volume of the idx element of the pore volume distribution
  @param idx is the index to get
  @return the volume of that element in the pore size distribution (nm3)
  */
  /*
  double getMasterPoreVolumeVolume(const int idx) {
    try {
      if (idx >= masterPoreVolume_.size()) {
        throw EOBException("Lattice", "getMasterPoreVolumeVolume",
                           "masterPoreVolume_", masterPoreVolume_.size(),
                           (int)idx);
      }
    } catch (EOBException ex) {
      ex.printException();
      exit(1);
    }
    return (masterPoreVolume_[idx].volume);
  }
  */

  /**
  @brief Get the volume fraction saturated  of the idx element of the pore
  volume distribution
  @param idx is the index to get
  @return the volume fraction saturated of that element in the pore size
  distribution
  */
  /*
  double getMasterPoreVolumeVolfrac(const int idx) {
    try {
      if (idx >= masterPoreVolume_.size()) {
        throw EOBException("Lattice", "getMasterPoreVolumeVolfrac",
                           "masterPoreVolume_", masterPoreVolume_.size(),
                           (int)idx);
      }
    } catch (EOBException ex) {
      ex.printException();
      exit(1);
    }
    return (masterPoreVolume_[idx].volfrac);
  }
  */

  /**
  @brief Get the largest diameter of pores containing electrolyte
  @return the diameter of the largest pore containing electrolyte
  */
  double getLargestSaturatedPore(void) {
    double capsize = 1000.0; // nm of capillary pores
    int size = masterPoreVolume_.size();
    for (int i = 0; i < size; i++) {
      if (masterPoreVolume_[i].volfrac < 1.0) {
        return (masterPoreVolume_[i].diam);
      }
    }
    return (capsize);
  }

  /**
  @brief Get the number of sites of water that must be added after a time step.

  @note Currently only used in sulfate attack simulations.

  @return the amount of water that must be added [site units]
  */
  double getWaterChange(void) const { return waterChange_; }

  /**
  @brief Set the number of sites of water that must be added after a time step.

  @note NOT USED.

  @param the number of sites of water that must be added [site units]
  */
  void setWaterChange(double waterchangeval) { waterChange_ = waterchangeval; }

  /**
  @brief Increment the number of sites of water that must be added after a time
  step.

  @param the extra number of sites of water that must be added [site units]
  */
  void dWaterChange(double dwaterchangeval) { waterChange_ += dwaterchangeval; }

  /**
  @brief Implement conversion of Al-bearing phases to ettringite.

  This method is intended only for simulating sulfate attack.  The method
  locates all the Al-bearing phases that are driven to transform. It also
  calculates the volume of free space adjacent to this site to determine whether
  crystallization pressure will arise.  If so, the method calculates the
  crystallization pressure and crystallization stress-free strain.  It then
  applies the expansion strain so that the new stress field can be calculated by
  the ThermalStrain FE model object.

  @todo Consider breaking this method into smaller pieces for maintenance and
  readability.

  @todo Generalize this for any phase transformation.

  @param alphaseid is the microstructure phase id of the Al-bearing phase to
  dissolve
  @param netsitesAlphaseid is the number of Al-bearing sites to dissolve
  @param ettrid is the microstructure phase id of ettringite
  @param netsitesEttrid is the number of ettringite sites to grow
  @return vector (na,ne) where na is the number of Al-bearing sites actually
  changed, and ne is the number of ettringite sites actually grown
  */
  vector<int> transform(int alphaseid, int netsitesAlphaseid, int ettrid,
                        int netsitesEttrid, double volumeratio);

  /**
  @brief Set a pointer to the AppliedStrain object for the simulation.

  @param elas is a pointer to the AppliedStrain object for the simulation
  */
  void setFEsolver(AppliedStrain *AppliedStrainSolver) {
    FEsolver_ = AppliedStrainSolver;
  }

  /**
  @brief Write a contiguous subvolume of the lattice to a file.

  @param fileName is the file name to which the subvolume will be written
  @param centerste is a pointer to the central site within the subvolume
  @param size is the extent of the subvolume in each direction away from the
  center site
  @return a list of the site indices belonging to the subvolume that was written
  */
  vector<int> writeSubVolume(string fileName, Site *centerste, int size);

  /**
  @brief Assign isotropic expansion strain at a set of prescribed sites.

  This function changes the strain components of a site already in the
  list of expansion sites.  If the prescribed site is not already in the
  list of expansion sites, then the site will be added to that list.

  @todo Consider changing the name of this method to applyExpansion

  @param alnb is the collection of site indices to which strain will be assigned
  @param exp is the isotropic expansion strain to set
  */
  void applyExpansion(vector<int> alnb, double exp);

  /**
  @brief Estimate the surface areas of all solid phases
  with the aqueous solution, in units of m2 per 100 g of
  total solids

  @param solidMass is the mass of all solids in g
  */
  void calcSurfaceAreas(void) {
    surfaceArea_.resize(numMicroPhases_, 0.0);
    specificSurfaceArea_.resize(numMicroPhases_, 0.0);
    for (int i = 0; i < numMicroPhases_; ++i) {
      calcSurfaceArea(i);
    }
  }

  /**
  @brief Estimate the surface area of a phase with the aqueous
  solution, in units of m2 per 100 g of total solids

  @param phaseid is the id of the microstructure phase
  */
  void calcSurfaceArea(int phaseid);

  /**
  @brief Return the current surface area of a phase with the aqueous
  solution.

  @param phaseid is the id of the microstructure phase
  @return the estimated surface area [m2 per 100 g of all solid]
  */
  double getSurfaceArea(int phaseid) {
    // if (phaseid > -1 && phaseid < surfaceArea_.size()) {
    return surfaceArea_[phaseid];
    // }
    // return 0.0;
  }

  vector<double> getSurfaceArea(void) { return surfaceArea_; }

  void resetSurfaceArea(vector<double> vect) { surfaceArea_ = vect; }

  /**
  @brief Return the current specific surface area of a phase with the aqueous
  solution.

  @param phaseid is the id of the microstructure phase
  @return the estimated specific surface area [m2 per g of this phase]
  */
  // double getSpecificSurfaceArea(int phaseid) {
  //   if (phaseid > -1 && phaseid < specificSurfaceArea_.size()) {
  //     return specificSurfaceArea_[phaseid];
  //   }
  //   return 0.0;
  // }

  /**
  @brief Return the combined specific surface area of cementitious
  components

  @return the estimated specific surface area [m2 per kg of all cement]
  */
  double getCementSpecificSurfaceArea(void) {
    double cemmass = 0.0;
    double allsurf = 0.0;
    double cemsurf = 0.0;
    double allsolidmass = 0.0;
    double thismass, thisssa;
    for (int i = 0; i < numMicroPhases_; ++i) {
      if (i != VOIDID && i != ELECTROLYTEID) {
        allsolidmass += chemSys_->getMicroPhaseMass(i);
        allsurf += getSurfaceArea(i);
        if (chemSys_->isCementComponent(i)) {
          thisssa = getSurfaceArea(i); // m2 component /(100g solid)
          thismass =
              chemSys_->getMicroPhaseMass(i); // g component/(100 g solid)
          cemsurf += thisssa;
          cemmass += thismass; // g all cement / (100g solid)
        }
      }
    }
    // At this point ssa is the total surface area of cement per 100 g of all
    // solids, and cemmass is the total g of cement per 100 g of all solids.
    // So ssa/cemmass has units of m2 per g of cement
    // Multiply that by 1000.0 to get units of m2/(kg of cement)
    // if (verbose_) {
    //    cout << "URANIUM all solid mass = " << allsolidmass << " g / (100 g
    //    solid)"
    //          << endl;
    //    cout << "URANIUM all surface = " << allsurf << " m2 / (100 g solid)"
    //            << endl;
    //    cout << "URANIUM cement mass = " << cemmass << " g / (100 g solid)" <<
    //    endl; cout << "URANIUM cement surface = " << cemsurf << " m2 / (100 g
    //    solid)"
    //            << endl;
    //    cout.flush();
    // }
    if (cemmass > 0.0) {
      cemsurf *= (1000.0 / cemmass);
    } else {
      cemsurf = 0.0;
    }
    return cemsurf;
  }

  /**
  @brief Get the sorted distribution of domain sizes

  @param phaseid is the id of the phase to query
  @param numsites is the maximum number of sites to store and sort
  @param maxisze is the maxmimum linear size of interest
  @param sortorder is 0 if sorting in descending order, nonzero otherwise
  @return an STL list of the site ids according to the distribution
  */
  vector<int> findDomainSizeDistribution(int phaseid, const int numsites,
                                         int maxsize, int sortorder);

  /**
  @brief Estimate the <i>linear size</i> of a domain

  @param siteid is the id of the microstructure phase
  @param maxsize is the maxmimum linear size of interest
  @return the edge length of the maximum cube that contains the same phase
  */
  int findDomainSize(int siteid, int maxsize);

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

  vector<int> getCount(void) { return count_; }

  int getCount(int phId) { return count_[phId]; }

  void setCount(vector<int> vect) { count_ = vect; }

  int getInterfaceSize(void) { return interface_.size(); }

  Interface getInterface(int i) { return interface_[i]; }

  void setGrowthSites(int i, vector<Isite> vect) {
    interface_[i].setGrowthSites(vect);
  }

  void setDissolutionSites(int i, vector<Isite> vect) {
    interface_[i].setDissolutionSites(vect);
  }
  void setInterfaceMicroPhaseId(int i, int mPhId) {
    interface_[i].setMicroPhaseId(mPhId);
  }

  double getInitSolidMass(void) { return initSolidMass_; }

  void findIsolatedClusters(void);

  void populateElementData(void);

  string getElemSymb(int index) { return cfgElem_[index].symb; }

  double callRNG(void) {
    numRNGcall_0_++;
    if (numRNGcall_0_ == LONG_MAX) {
      numRNGcallLONGMAX_++;
      numRNGcall_0_ = 0;
    }
    lastRNG_ = rg_->Ran3();
    return lastRNG_;
  }

  long int getNumRNGcall_0(void) { return numRNGcall_0_; }

  long int getNumRNGcallLONGMAX(void) { return numRNGcallLONGMAX_; }

  double getLastRNG(void) { return lastRNG_; }

  void setRNGseed(int seed) { rg_->setSeed(seed); }

  int getRNGseed(void) { return latticeRNGseed_; }

  void resetRNG(long int val_0, long int valLONGMAX, double valRNG) {
    // latticeRNGseed_ = seed;
    rg_->setSeed(latticeRNGseed_);
    numRNGcall_0_ = val_0;
    numRNGcallLONGMAX_ = valLONGMAX;
    // long int count_0 = 0, count_1 = 0;
    long int j0, j1, j11;
    double lastRNGreset = 1.e-16;
    for (j1 = 1; j1 <= numRNGcallLONGMAX_; j1++) {
      for (j11 = 1; j11 <= LONG_MAX; j11++) {
        lastRNGreset = rg_->Ran3();
      }
    }
    for (j0 = 1; j0 <= val_0; j0++) {
      lastRNGreset = rg_->Ran3();
    }
    lastRNG_ = lastRNGreset;

    // cout << endl
    //      << "  Lattice::resetRNG cyc/whileCount/latticeRNGseed_: " << cyc
    //      << " / " << whileCount << " / " << latticeRNGseed_ << endl;
    // cout << "  Lattice::resetRNG "
    //         "numRNGcall_0_/numRNGcallLONGMAX_/lastRNGreset/valRNG: "
    //      << numRNGcall_0_ << " / " << numRNGcallLONGMAX_ << " / "
    //      << lastRNGreset << " / " << valRNG << endl;

    if (abs(lastRNGreset - valRNG) > 1.e-16) {
      cout << endl << "Lattice::resetRNG FAILED => exit" << endl;
      exit(0);
    }
  }

  void increaseLatticeVolume(void);

  void checkSite(int stId) {
    // int phId = site_[stId].getMicroPhaseId();
    cout << endl << " Lattice::checkSite( " << stId << " ):" << endl;
    cout << "    phaseId : " << site_[stId].getMicroPhaseId() << endl;
    cout << "    inDissInterfacePos_ : " << site_[stId].getInDissInterfacePos()
         << endl;
    if (site_[stId].getInDissInterfacePos() != -1) {
      cout << "     in dissInterface on pos inDissInterfacePos_ : "
           << interface_[site_[stId].getMicroPhaseId()].getDissolutionSitesId(
                  site_[stId].getInDissInterfacePos())
           << endl;
    }

    vector<int> growth = site_[stId].getGrowthPhases();
    int size = growth.size();
    int k;
    cout << endl << " growth_.size() : " << size << endl;
    for (k = 0; k < size; k++) {
      cout << "       k = " << k << "   growth_[k] = " << growth[k] << endl;
    }
    cout << endl << " inGrowInterfacePos_ : " << endl;
    for (k = 0; k < numMicroPhases_; k++) {
      cout << "       k = " << k << "   site_[stId].getInGrowInterfacePos(k) = "
           << site_[stId].getInGrowInterfacePos(k) << endl;
    }

    cout << endl
         << "     in growInterfaces on pos inGrowInterfacePos_ :" << endl;
    for (k = 0; k < numMicroPhases_; k++) {
      cout << "       k_ = " << k << endl;
      cout.flush();
      if (site_[stId].getInGrowInterfacePos(k) > -1) {
        size = growthInterfaceSize_[k];
        cout << "       k = " << k
             << "   pos = " << site_[stId].getInGrowInterfacePos(k)
             << "   size = " << size << endl;
        cout.flush();
        cout << "            siteId in grInt = "
             << interface_[k].getGrowthSitesId(
                    site_[stId].getInGrowInterfacePos(k))
             << endl;
        cout.flush();
      }
    }
  }

  vector<int> getGrowthInterfaceSize(void) { return growthInterfaceSize_; }

  void setGrowthInterfaceSize(vector<int> vect) { growthInterfaceSize_ = vect; }

  vector<int> getDissolutionInterfaceSize(void) {
    return dissolutionInterfaceSize_;
  }

  int getDissolutionInterfaceSize(int phId) {
    return dissolutionInterfaceSize_[phId];
  }

  void setDissolutionInterfaceSize(vector<int> vect) {
    dissolutionInterfaceSize_ = vect;
  }

  vector<int> chooseNucleationSitesRND(int phaseID, int numLeft);

  vector<int> chooseNucleationSitesAFF(int phaseID, int numLeft);

  /**
  @brief Convert (switch to electrolyte) the prescribed number of
  sites of each microphase that has to dissolve

  @param dissPhaseIDVect is the vector of microphase IDs that must dissolve
  @param numSiteDissVect is a vector containing the number of voxels to dissolve
  for each microphase ID in dissPhaseIDVect
  @param dissPhNameVect is a vector containing the name of each microphase in
  dissPhaseIDVect
  @param numtoadd_D is the number of sites switched by this call
  @param totalTRC is the total call number of the changeMicrostructure method
  @return -1 if all requested numbers of voxels in numSiteDissVect have been
  switched or the ID of the first microphase (in dissPhaseIDVect) for which this
  was not possible
  */
  vector<int>
  transformPhase(int ettrid, int netsitesEttrid, vector<int> dissPhaseIDVect,
                 vector<int> numSiteDissVect, vector<string> dissPhNameVect,
                 vector<double> volumeRatio, int &numadded_D, int totalTRC);

  void createGrowingVectSA(void);

  void transformGrowPhase(Site *ste, int growPhID, int totalTRC);

  void transformChangePhase(Site *ste, int oldPhId, int newPhId, int totalTRC);

  vector<int> getAllSitesPhId(void) {
    vector<int> allPhId(numSites_, 0);
    for (int i = 0; i < numSites_; i++) {
      allPhId[i] = site_[i].getMicroPhaseId();
    }
    return allPhId;
  }

}; // End of Lattice class
#endif // SRC_THAMESLIB_LATTICE_H_

///
/// The functions below are used to aid in comparison of one site to another,
/// by which means lists of the sites can be sorted.
///

#ifndef CMPFUNCS
#define CMPFUNCS

/**
@brief Compare two sites, returning true is the first site is "less than" the
second.

The comparison is made on the basis of what THAMES loosely calls the
<i>weighted mean curvature</i>, (wmc).  A site with high wmc is a site where
dissolution of a phase is likely to occur, and growth of another phase is
unlikely to occur. Conversely, a site with a low wmc is a site where growth of
a phase is likely to occur but dissolution of a phase is unlikely to occur.

@param s1 is a pointer to the first site in the comparison
@param s2 is a pointer to the second site in the comparison
@return true if the first site has lower wmc than the second, false otherwise
*/
bool cmp(const Site *s1, const Site *s2);

/**
@brief Sort two sites based on their affinity for a given phase.

The comparison is made on the basis of what THAMES loosely calls the
<i>affinity</i>.  A site with high affinity is a site where growth of a phase
is more likely to occur because of an affinity between it and the interface.

@param s1 is the first site in the comparison
@param s2 is the second site in the comparison
@return true if the first site has <i>greater</i> affinity than the second,
false otherwise
*/
bool affinitySort(const Isite s1, const Isite s2);

#endif // CMPFUNCS
