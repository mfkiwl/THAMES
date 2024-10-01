/**
@file Lattice.h
@brief Declaration of the Lattice class for storing the 3D microstructure

THAMES defines a Lattice class that is instantiated to a Lattice
object at the beginning of the program's execution.  The lattice defines the
three-dimensional environment within which a cement paste microstructure
exists, hydrates, and possibly deteriorates.
*/

#ifndef LATTICEH
#define LATTICEH

#include <algorithm>
#include <cmath>
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

#include "AppliedStrain.h"
#include "ChemicalSystem.h"
#include "Interface.h"
#include "Isite.h"
#include "RanGen.h"
#include "Site.h"
#include "global.h"
#include "../version.h"

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

struct structGrowVect
{
  int id;
  int posVect;
  int affinity;
};

struct structDissVect
{
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
  string version_;             /**< THAMES version for header information */
  string thamesVersion_;
  string jobroot_;             /**< The root name for output files */

  RanGen *rg_;                 /**< Pointer to random number generator object */
  int latticeRNGseed_;
  long int numRNGcall_0_, numRNGcallLONGMAX_;
  double lastRNG_;

  unsigned int xdim_;          /**< Number of sites in the x dimension */
  unsigned int ydim_;          /**< Number of sites in the y dimension */
  unsigned int zdim_;          /**< Number of sites in the z dimension */
  double resolution_;          /**< Voxel edge length [micrometers] */
  vector<Site> site_;          /**< 1D list of Site objects (site = voxel) */
  unsigned int numsites_;      /**< Total number of sites */
  unsigned int siteneighbors_; /**< Number of neighbor sites to a given site */
  ChemicalSystem *chemSys_;    /**< Pointer to simulation's ChemicalSystem */
  AppliedStrain *FEsolver_;    /**< Pointer to simulation's FE elastic solver */
  vector<Interface> interface_;   /**< List of the different interface objects
                                          in the microstructure */
  double wsratio_;                /**< Water-to-solids mass ratio */
  vector<double> volumefraction_; /**< Array of volume fractions of each
                                          microstructure phase */
  vector<double> initvolumefraction_; /**< Array of initial volume fractions of
                                         each microstructure phase */
  vector<int> count_; /**< Number of sites of each different type */
  vector<double> SI_; /**< Current saturation indices */

  map<int, vector<double>>
      expansion_; /**< Map of expansion strain of each voxel */
  map<int, vector<int>> expansion_coordin_; /**< Map of coordinates of sites
                                               with local expansion strain */
  double waterchange_;          /**< How much water must be added or subtracted
                                        due to hydration or deterioration */
  double microstructurevolume_; /**< Microstructure volume in GEM
volume units */
  double initialmicrostructurevolume_; /**< Initial microstructure volume in GEM
volume units */
  double capillaryporevolume_;         /**< Total volume of capillary pores */
  double capillaryporevolumefraction_; /**< Total volume fraction of capillary
                                          pores */
  double subvoxelporevolume_;          /**< Total volume of subvoxel pores */
  double nonsolidvolume_;              /**< Total volume not solid */
  double solidvolumewithpores_;        /** Total solid volume including their
internal pore volume */
  double watervolume_;                 /** volume of electrolyte in GEM
volume units */
  double voidvolume_;                  /** volume of void in GEM volume
units */
  double capillarywatervolume_;        /**< Volume of capillary pore water */
  double capillaryvoidvolume_;         /**< Volume of capillary void space
(no water) */
  double subvoxelwatervolume_;         /**< Volume of water in subvoxel
pores in GEM units */
  double subvoxelporevolumefraction_;  /**< Total volume fraction of subvoxel
                                          pores */

  vector<struct PoreSizeVolume>
      masterporevolume_; /**< Pore size distribution and saturation */

  double time_;         /**< The current simulation time [days] */
  double temperature_;  /**< The current simulation temperature [K] */
  double oldtemp_;      /**< The temperature in the previous
                                time step [K] */
  double sattack_time_; /**< Simulation time at which to begin
                                simulation of sulfate attack [days] */
  double leach_time_;   /**< Simulation time at which to begin
                                simulation of leaching [days] */
  double surfacearea_;  /**< Total surface area [m<sup>2</sup>] */

  bool deptheffect_; /**< Whether or not PNG images should have
                             depth effect */
  bool verbose_;     /**< Flag to determine verbose output */
  bool warning_;     /**< Flag to determine warning message output */

  vector<chemElem> cfgElem_; /**< Holds periodic table information to output
                                files in cfg format */

  double initSolidMass_;

  double wcratio_;  /**< Water-to-cement mass ratio */

  int numMicroPhases_;   /**< Number of microphases */

  double particRadius_;  /**< used for graphical representation */

  vector<int> growthInterfaceSize_; /**< growth interface size of each microphase */
  vector<int> dissolutionInterfaceSize_; /**< dissolution interface size of each microphase */

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
  @param fileName is the name of the file containing the microstructure data
  @param verbose is true if extra messages are to be printed
  @param warning is true if warning messages are to be printed
  */
  Lattice(ChemicalSystem *cs, const string &fileName, const bool verbose,
          const bool warning);

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
  void setXDim(const unsigned int x) {
    xdim_ = x;
    numsites_ = (xdim_ * ydim_ * zdim_);
  }

  /**
  @brief Get the number of sites in the x dimension.

  @return the number of sites in the x dimension
  */
  unsigned int getXDim() const { return xdim_; }

  /**
  @brief Set the number of sites in the y dimension.

  @param y is the number of sites in the y dimension
  */
  void setYDim(const unsigned int y) {
    ydim_ = y;
    numsites_ = (xdim_ * ydim_ * zdim_);
  }

  /**
  @brief Get the number of sites in the y dimension.

  @return the number of sites in the y dimension
  */
  unsigned int getYDim() const { return ydim_; }

  /**
  @brief Set the number of sites in the z dimension.

  @param z is the number of sites in the z dimension
  */
  void setZDim(const unsigned int z) {
    zdim_ = z;
    numsites_ = (xdim_ * ydim_ * zdim_);
  }

  /**
  @brief Get the number of sites in the z dimension.

  @return the number of sites in the z dimension
  */
  unsigned int getZDim() const { return zdim_; }

  /**
  @brief Get the total number of lattice sites.

  The lattice is rectangular, so the total number of sites is
  `xdim_ * ydim_ * zdim_`, but we store this value as a class member to
  save having to compute it multiple times.

  @return the total number of lattice sites
  */
  unsigned int getNumsites() const { return numsites_; }

  /**
  @brief Set the volume fraction of a given microstructure phase.

  @param i is the index of the microstructure phase
  @param vfrac is the volume fraction to assign on a total microstructure basis
  */
  void setVolumefraction(unsigned int i, double vfrac) {
    volumefraction_[i] = vfrac;
  }

  /**
  @brief Set the initial volume fraction of a given microstructure phase.

  @param i is the index of the microstructure phase
  @param vfrac is the volume fraction to assign on a total microstructure basis
  */
  void setInitvolumefraction(unsigned int i, double vfrac) {
    initvolumefraction_[i] = vfrac;
  }

  /**
  @brief Set the water-solids mass ratio

  @param ws is the water-solids mass ratio
  */
  void setWsratio(const double ws) {
    wsratio_ = 0.0;
    if (ws > 0.0) {
      wsratio_ = ws;
    }
    return;
  }

  /**
  @brief Get the water-solids mass ratio

  @return the water-solids mass ratio
  */
  double getWsratio(void) const { return wsratio_; }
  double getWcratio(void) const { return wcratio_; }

  /**
  @brief Get the volume fraction of a given microstructure phase.

  This is simply the number of sites with a given phase divided by the
  total number of sites.

  @param i is the index of the microstructure phase
  @return the volume fraction of phase i on a total microstructure basis
  */
  double getVolumefraction(unsigned int i) {
    if (numsites_ == 0) {
      throw FloatException("Lattice", "getVolumefraction",
                           "Divide by zero (numsites_)");
    }
    return (volumefraction_[i]);
  }

  /**
  @brief Get the initial volume fraction of a given microstructure phase.

  This is simply the number of sites with a given phase divided by the
  total number of sites.

  @param i is the index of the microstructure phase
  @return the initial volume fraction of phase i on a total microstructure basis
  */
  double getInitvolumefraction(unsigned int i) {
    if (numsites_ == 0) {
      throw FloatException("Lattice", "getInitialvolumefraction",
                           "Divide by zero (numsites_)");
    }
    return (initvolumefraction_[i]);
  }

  /**
  @brief Calculate the subvoxel pore volume

  @param vol is the array of all microstructure phase volumes
  */
  void calcSubvoxelporevolume(vector<double> &vol);

  /**
  @brief Calculate the total volume of solids including
  subvoxel pore volume assigned to solids

  @param vol is the array of all microstructure phase volumes
  it
  */
  void calcSolidvolumewithpores(vector<double> &vol);

  /**
  @brief Get the total volume of solids including
  subvoxel pore volume assigned to solids

  @return the solid volume including subvoxel pore volume
  */
  double getSolidvolumewithpores(void) const { return solidvolumewithpores_; }

  /**
  @brief Calculate the non-solid volume

  @param vol is the array of all microstructure phase volumes
  it
  */
  void calcNonsolidvolume(vector<double> &vol);

  /**
  @brief Get or calculate the non-solid volume

  @return the non-solid volume
  */
  double getNonsolidvolume(void) const { return nonsolidvolume_; }

  /**
  @brief Get the number of neighbor sites each site has.

  This is simply the number of sites with a given phase divided by the
  total number of sites.

  @note NOT USED.

  @return the number of neighbor sites each site has
  */
  unsigned int getSiteneighbors() const { return siteneighbors_; }

  /**
  @brief Set the lattice resolution [micrometers].

  The lattice resolution is the physical length associated with the edge
  length of a site.

  @param res is the lattice resolution [micrometers]
  */
  void setResolution(const double res);

  /**
  @brief Get the lattice resolution [micrometers].

  The lattice resolution is the physical length associated with the edge
  length of a site.

  @note NOT USED.

  @return the lattice resolution [micrometers]
  */
  double getResolution() const { return resolution_; }

  /**
  @brief Set the simulation time [days].

  @note NOT USED.

  @param tval is the simulation time [days]
  */
  void setTime(const double tval) { time_ = tval; }

  /**
  @brief Get the simulation time [days].

  @note NOT USED.

  @return the simulation time [days]
  */
  double getTime() const { return time_; }

  /**
  @brief Get the simulation time at which to start sulfate attack simulation
  [days].

  @note NOT USED.

  @return the simulation time at which to start sulfate attack [days]
  */
  double getSattack_time() const { return sattack_time_; }

  /**
  @brief Set the simulation time at which to start sulfate attack simulation
  [days].

  @param sattacktime is the simulation time at which to start sulfate attack
  [days]
  */
  void setSattack_time(const double sattacktime) {
    sattack_time_ = sattacktime;
  }

  /**
  @brief Get the simulation time at which to start leaching simulation [days].

  @note NOT USED.

  @return the simulation time at which to start leaching [days]
  */
  double getLeach_time() const { return leach_time_; }

  /**
  @brief Set the simulation time at which to start leaching simulation [days].

  @param leachtime is the simulation time at which to start leaching [days]
  */
  void setLeach_time(const double leachtime) { leach_time_ = leachtime; }

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
  void setJobroot(string jobname) { jobroot_ = jobname; }

  /**
  @brief Add a site at location (xp,yp,zp) to the lattice.

  The site is checked for valid coordinates.  If valid a new Site object
  is created and pushed back onto the class's `site_` vector.

  @param xp is the x coordinate of the site to add
  @param yp is the y coordinate of the site to add
  @param zp is the z coordinate of the site to add
  */
  void addSite(const unsigned int xp, const unsigned int yp,
               const unsigned int zp);

  /**
  @brief Get the x coordinate of a site with a given index in the 1D `site_`
  array.

  @param i is the index of the site in the class's `site_` array
  @return the x coordinate
  */
  unsigned int getX(const unsigned int i) const {
    return (site_[i].getX());
  }

  /**
  @brief Get the y coordinate of a site with a given index in the 1D `site_`
  array.

  @param i is the index of the site in the class's `site_` array
  @return the y coordinate
  */
  unsigned int getY(const unsigned int i) const {
    return (site_[i].getY());
  }

  /**
  @brief Get the z coordinate of a site with a given index in the 1D `site_`
  array.

  @param i is the index of the site in the class's `site_` array
  @return the x coordinate
  */
  unsigned int getZ(const unsigned int i) const {
    return (site_[i].getZ());
  }

  /**
  @brief Get a site's index in the 1D `site_` array, given its (x,y,z)
  coordinates.

  @param ix is the x coordinate of the site
  @param iy is the x coordinate of the site
  @param iz is the x coordinate of the site
  @return the index of the site in the `site_` array
  */
  unsigned int getIndex(int ix, int iy, int iz) const;

  /**
  @brief Get the collection of site indices neighboring a given site.

  @param sitenum is the index of the site in question
  @param size is the maximum distance defining the neighborhood [sites]
  @return a list of site indices for all neighbors within the maximum distance
  */
  vector<unsigned int> getNeighborhood(const unsigned int sitenum,
                                       const int size);

  /**
  @brief Get a pointer to a Site object at a given index in the `site_` array.

  @param index is the index of the Site object in the `site_` array
  @return a pointer to the Site object in question
  */
  Site *getSite(int index) {
    return &site_[index];
  }

  /**
  @brief Designate a site as damaged.

  The site to be damaged is specified by its index in the `site_` array.

  @note NOT USED.

  @param index is the index of the Site object in the `site_` array
  */
  void setDamage(int index) {
    site_[index].setDamage();
  }

  /**
  @brief Change the wmc (weighted mean curvature) of a site by a prescribed
  amount.

  @param index is the index of the Site object in the `site_` array
  @param dwmcval is the increment to add to the wmc
  */
  void dWmc(int index, double dwmcval) {
    site_[index].setWmc(site_[index].getWmc() + dwmcval);
  }

  /**
  @brief Compute normalized initial microstructure phase masses

  Given the initial masses of all phases in the microstructure,
  this method scales them to 100 grams of solid.  In the process,
  this method also sets the initial moles of water in the
  chemical system definition.

  @param microPhaseMass is a vector of all the microstructure masses
  @param solidMass is the combined mass of all the solids
  */
  void normalizePhaseMasses(vector<double> microPhaseMass, double cementMass,
                            double solidMass);

  /**
  @brief Master method to locate the interfaces for each phase in the
  microstructure.

  */
  void findInterfaces(void);

  /**
  @brief Add (grow i.e. switch from electrolyte) the prescribed number of
  sites of each microphase that has to grow

  @param growPhaseIDVect is the vector of microphase IDs that must grow
  @param numSiteGrowVect is a vector containing the number of voxels to add for each
  microphase ID in growPhaseIDVect
  @param growPhNameVect is a vector containing the name of each microphase in growPhaseIDVect
  @param numtoadd_G is the number of sites switched by this call
  @param totalTRC is the total call number of the changeMicrostructure method
  @return the actual number of sites that were changed for each microphase ID from
  the input growPhaseIDVect vector
  */
  vector<int> growPhase(vector<int> growPhaseIDVect, vector<int> numSiteGrowVect,
                        vector<string> growPhNameVect, int &numadded_G, int totalTRC);

  /**
  @brief create a new growth interface for a given phase (phaseID) having a
  size of numLeft sites; this is necessary when the "growth" of all requested
  sites for this phase was not possible because the size of the
  corresponding growth interface was zero.

  @param phaseid is the id of the microstructure phase to nucleate
  @param numLeft is the number of sites to nucleate/create for this phase
  @return the actual size of the new interface (must equals numLeft!)
  */
  void nucleatePhaseAff(int phaseID,int numLeft);
  void nucleatePhaseRnd(int phaseID,int numLeft);

  /**
  @brief Remove (dissolve i.e. switch to electrolyte) the prescribed number of
  sites of each microphase that has to dissolve

  @param dissPhaseIDVect is the vector of microphase IDs that must dissolve
  @param numSiteDissVect is a vector containing the number of voxels to dissolve for each
  microphase ID in dissPhaseIDVect
  @param dissPhNameVect is a vector containing the name of each microphase in dissPhaseIDVect
  @param numtoadd_D is the number of sites switched by this call
  @param totalTRC is the total call number of the changeMicrostructure method
  @return -1 if all requested numbers of voxels in numSiteDissVect have been switched
  or the ID of the first microphase (in dissPhaseIDVect) for which this was not possible
  */
  int dissolvePhase(vector<int> dissPhaseIDVect, vector<int> numSiteDissVect,
                    vector<string> dissPhNameVect, int &numadded_D, int totalTRC);

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
  void setMicroPhaseId(Site *s, const unsigned int i) {
      count_[s->getMicroPhaseId()]--;
      s->setMicroPhaseId(i);
      count_[i]++;
  }

  /**
  @brief Set the phase id of a given site, specified by the site's index number.

  @param sitenum is the index of the site in the `site_` array
  @param i is the phase index to set at that site
  */
  void setMicroPhaseId(const int sitenum, const unsigned int i) {
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
  void addDissolutionSite(Site *loc, unsigned int pid);

  /**
  @brief Add a site to the list of sites where growth of a given phase can
  occur.

  @param loc is a pointer to the Site object to add to the list of potential
  growth sites
  @param pid is the microstructure phase id
  */
  void addGrowthSite(Site *loc, unsigned int pid);

  /**
  @brief Remove a site from the list of sites where dissolution of a given phase
  can occur.

  @param loc is a pointer to the Site object to remove from the list of
  potential dissolution sites
  @param pid is the microstructure phase id
  */
  void removeDissolutionSite(Site *loc, unsigned int pid);

  /**
  @brief Remove a site from the list of sites where growth of a given phase can
  occur.

  @param loc is a pointer to the Site object to remove from the list of
  potential growth sites
  @param pid is the microstructure phase id
  */
  void removeGrowthSite_diss(Site *loc, unsigned int pid);
  void removeGrowthSite_grow(Site *ste0, int pid);
  void removeGrowthSite_nucleation(Site *loc, unsigned int pid);

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

  @param time is is the simulation time [days]
  @param simtype is the type of simulation (hydration, leaching, etc)
  @param capWater is true if there is any capillary pore water in the system.
  @param numDiff is the maximum number of voxels belonging to a given microphase,
  voxels that can be dissolved according to the system configuration (lattice)
  @param phDiff is the microphase ID for which a the number of voxels that can be
  dissolved is smaller than the number requested by the corresponding kinetic model
  @param nameDiff is the name of this microphase
  @param whileCount counts the number of changeMicrostructure calls for a given cycle
  (cyc)
  @param cyc (cycle) is the iteration number in main iteration loop in
  Controller::doCycle - each cycle corresponds to a time step

  @return zero if okay or nonzero if not all requested voxels
  for a certain microphase ID (phDiff) can be dissolved
  */
  int changeMicrostructure(double time, const int simtype, bool &capWater,
                           int &numDiff, int &phDiff, string &nameDiff,
                           int whileCount, int cyc);

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
                                    vector<double> &vfrac, int volSize);

  /**
  @brief Calculate the pore size distribution data

  */
  void calculatePoreSizeDistribution(void);

  /**
  @brief Write the pore size distribution data to a file

  @param curtime is the current time in days
  @param simtype is the sumulation tyupe
  @param root is the root name of the output file to create
  */
  void writePoreSizeDistribution(double curtime, const int simtype,
                                 const string &root);

  /**
  @brief Write the microstructure colors to a file

  This is done to save processing the chemistry.xml file just to get the colors
  and will make post-processing of images easier.

  @param root is the root name of the output file to create
  */
  void writeMicroColors(const string &root);

  /**
  @brief Write the 3D microstructure to a file.

  The microstructure output file will indicate the phase id at each site.

  @param curtime is the current time in days
  @param simtype is the sumulation tyupe
  @param root is the root name of the output file to create
  */
  void writeLattice(double curtime, const int simtype, const string &root);
  void writeLatticeIni(double curtime);
  void writeLatticeXYZ(double curtime, const int simtype, const string &root);
  void writeLatticeCFG(double curtime, const int simtype, const string &root);

  /**
  @brief Write the 3D microstructure to a file.

  The damage output file is binary, each site either being damaged or not.

  @param curtime is the current time in days
  @param root is the root name of the output file to create
  */
  void writeDamageLattice(double curtime, const string &root);

  /**
  @brief Write the 3D microstructure to a png file that can be immediately
  rendered.

  @param curtime is the current time in days
  @param root is the root name of the png output file to create
  */
  void writeLatticePNG(double curtime, const int simtype, const string &root);

  /**
  @brief Write the 3D microstructure to a png file that can be immediately
  rendered.

  The damage output file is binary, each site either being damaged or not.

  @param curtime is the current time in days
  @param root is the root name of the png output file to create
  */
  void writeDamageLatticePNG(double curtime, const string &root);

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
  void makeMovie(const string &root);

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
    string msg;
    map<int, vector<double>>::iterator p = expansion_.find(index);
    if (p != expansion_.end()) {
      return p->second;
    } else {
      msg = "Could not find expansion_ match to index provided";
      throw EOBException("Lattice", "getExpansion", msg, expansion_.size(), 0);
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
  vector<int> getExpansionCoordin(int index) {
    string msg;
    map<int, vector<int>>::iterator p = expansion_coordin_.find(index);
    if (p != expansion_coordin_.end()) {
      return p->second;
    } else {
      msg = "Could not find expansion_coordin_ match to index provided";
      throw EOBException("Lattice", "getExpansionCoordin", msg,
                         expansion_coordin_.size(), 0);
    }
  }

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
  void setExpansionCoordin(int index, vector<int> coordin) {
    string msg;
    map<int, vector<int>>::iterator p = expansion_coordin_.find(index);
    if (p == expansion_coordin_.end()) {
      expansion_coordin_.insert(make_pair(index, coordin));
    }
  }

  /**
  @brief Get the microstructure volume

  @return the microstructure volume (GEMS volume units)
  */
  double getMicrostructurevolume(void) const {
    return (chemSys_->getMicroVolume());
  }

  /**
  @brief Get the initial microstructure volume

  @return the initial microstructure volume (GEMS volume units)
  */
  double getInitialmicrostructurevolume(void) const {
    return (chemSys_->getInitMicroVolume());
  }

  /**
  @brief Get the total capillary pore volume

  @return the volume of capillary pores (GEMS volume units)
  */
  double getCapillaryporevolume(void) const { return capillaryporevolume_; }

  /**
  @brief Set the capillary pore volume

  @param capillaryporevolume is the capillary pore volume (GEMS volume units)
  */
  void setCapillaryporevolume(double capillaryporevolume) {
    capillaryporevolume_ = capillaryporevolume;
  }

  /**
  @brief Get the total capillary pore volume fraction
  This is calculated on a total system volume basis

  @return the volume fraction of capillary pores (microstructure basis)
  */
  double getCapillaryporevolumefraction(void) const {
    return capillaryporevolumefraction_;
  }

  /**
  @brief Set the capillary pore volume fraction
  This is calculated on a total system volume basis

  @param capillaryporevolumefraction is the capillary pore volume
  fraction (microstructure basis)
  */
  void
  setCapillaryporevolumefraction(const double capillaryporevolumefraction) {
    capillaryporevolumefraction_ = capillaryporevolumefraction;
  }

  /**
  @brief Get the total subvoxel pore volume

  @return the volume of subvoxel pores (GEMS volume units)
  */
  double getSubvoxelporevolume(void) const { return subvoxelporevolume_; }

  /**
  @brief Set the subvoxel pore volume

  @param subvoxelporevolume is the subvoxel pore volume (GEMS volume units)
  */
  void setSubvoxelporevolume(const double subvoxelporevolume) {
    subvoxelporevolume_ = subvoxelporevolume;
  }

  /**
  @brief Set the subvoxel pore volume

  @param subvoxelporevolume is the subvoxel pore volume (GEMS volume units)
  */
  void setNonsolidvolume(const double nonsolidvolume) {
    nonsolidvolume_ = nonsolidvolume;
  }

  /**
  @brief Get the capillary water volume

  @param vol is the volume of each microstructure phase
  */
  void calcCapillarywatervolume(vector<double> &vol);

  /**
  @brief Get the capillary water volume

  @return the capillary water volume
  */
  double getCapillarywatervolume(void) const { return capillarywatervolume_; }

  /**
  @brief Set the capillary water volume

  @param capillarywatervolume is the capillary water volume (GEMS volume units)
  */
  void setCapillarywatervolume(const double capillarywatervolume) {
    capillarywatervolume_ = capillarywatervolume;
  }

  /**
  @brief Get the capillary void volume

  @param vol is the volume of each microstructure phase
  @param calc is true only if calculating instead of just returning
  */
  void calcCapillaryvoidvolume(vector<double> &vol);

  /**
  @brief Get the capillary void volume

  @return the capillary void volume
  */
  double getCapillaryvoidvolume(void) const { return capillarywatervolume_; }

  /**
  @brief Set the capillary void volume

  @param capillaryvoidvolume is the capillary void volume (GEMS volume units)
  */
  void setCapillaryvoidvolume(const double capillaryvoidvolume) {
    capillaryvoidvolume_ = capillaryvoidvolume;
  }

  /**
  @brief Get the total subvoxel pore volume fraction
  This is calculated on a total system volume basis

  @return the volume fraction of subvoxel pores (microstructure basis)
  */
  double getSubvoxelporevolumefraction(void) const {
    return subvoxelporevolumefraction_;
  }

  /**
  @brief Set the subvoxel pore volume fraction
  This is calculated on a total system volume basis

  @param subvoxelporevolumefraction is the subvoxel pore volume
  fraction (microstructure basis)
  */
  void setSubvoxelporevolumefraction(const double subvoxelporevolumefraction) {
    subvoxelporevolumefraction_ = subvoxelporevolumefraction;
  }

  /**
  @brief Set the master pore volume distribution

  @param masterporevolume is the pore volume distribution
  */
  void
  setMasterporevolume(const vector<struct PoreSizeVolume> masterporevolume) {
    masterporevolume_ = masterporevolume;
    return;
  }

  /**
  @brief Set the master pore volume distribution of a particular size

  @param idx is the index to set
  @param diam is the diameter in nm
  @param volume is the volume of pores this size, in nm3
  @param volfrac is the volume fraction of this size filled with electrolyte
  */
  void setMasterporevolume(const int idx, const double diam,
                           const double volume, const double volfrac) {
    try {
      if (idx >= masterporevolume_.size()) {
        throw EOBException("Lattice", "setMasterporevolume",
                           "masterporevolume_", masterporevolume_.size(),
                           (int)idx);
      }
      masterporevolume_[idx].diam = diam;
      masterporevolume_[idx].volume = volume;
      masterporevolume_[idx].volfrac = volfrac;
    } catch (EOBException ex) {
      ex.printException();
      exit(1);
    }
    return;
  }

  /**
  @brief Get the master pore volume distribution of a particular size
  @param idx is the index to get
  @return the structure holding the pore size distribution data for that element
  */
  struct PoreSizeVolume getMasterporevolume(const int idx) {
    try {
      if (idx >= masterporevolume_.size()) {
        throw EOBException("Lattice", "getMasterporevolume",
                           "masterporevolume_", masterporevolume_.size(),
                           (int)idx);
      }
    } catch (EOBException ex) {
      ex.printException();
      exit(1);
    }
    return (masterporevolume_[idx]);
  }

  /**
  @brief Get the diameter of the idx element of the pore volume distribution
  @param idx is the index to get
  @return the diameter of that element in the pore size distribution (nm)
  */
  double getMasterporevolumeDiam(const int idx) {
    try {
      if (idx >= masterporevolume_.size()) {
        throw EOBException("Lattice", "getMasterporevolumeDiam",
                           "masterporevolume_", masterporevolume_.size(),
                           (int)idx);
      }
    } catch (EOBException ex) {
      ex.printException();
      exit(1);
    }
    return (masterporevolume_[idx].diam);
  }

  /**
  @brief Get the total volume of the idx element of the pore volume distribution
  @param idx is the index to get
  @return the volume of that element in the pore size distribution (nm3)
  */
  double getMasterporevolumeVolume(const int idx) {
    try {
      if (idx >= masterporevolume_.size()) {
        throw EOBException("Lattice", "getMasterporevolumeVolume",
                           "masterporevolume_", masterporevolume_.size(),
                           (int)idx);
      }
    } catch (EOBException ex) {
      ex.printException();
      exit(1);
    }
    return (masterporevolume_[idx].volume);
  }

  /**
  @brief Get the volume fraction saturated  of the idx element of the pore
  volume distribution
  @param idx is the index to get
  @return the volume fraction saturated of that element in the pore size
  distribution
  */
  double getMasterporevolumeVolfrac(const int idx) {
    try {
      if (idx >= masterporevolume_.size()) {
        throw EOBException("Lattice", "getMasterporevolumeVolfrac",
                           "masterporevolume_", masterporevolume_.size(),
                           (int)idx);
      }
    } catch (EOBException ex) {
      ex.printException();
      exit(1);
    }
    return (masterporevolume_[idx].volfrac);
  }

  /**
  @brief Get the largest diameter of pores containing electrolyte
  @return the diameter of the largest pore containing electrolyte
  */
  double getLargestSaturatedPore(void) {
    double capsize = 1000.0; // nm of capillary pores
    int size = masterporevolume_.size();
    for (int i = 0; i < size; i++) {
      if (masterporevolume_[i].volfrac < 1.0) {
        return (masterporevolume_[i].diam);
      }
    }
    return (capsize);
  }

  /**
  @brief Get the number of sites of water that must be added after a time step.

  @note Currently only used in sulfate attack simulations.

  @return the amount of water that must be added [site units]
  */
  double getWaterchange(void) const { return waterchange_; }

  /**
  @brief Set the number of sites of water that must be added after a time step.

  @note NOT USED.

  @param the number of sites of water that must be added [site units]
  */
  void setWaterchange(double waterchangeval) { waterchange_ = waterchangeval; }

  /**
  @brief Increment the number of sites of water that must be added after a time
  step.

  @param the extra number of sites of water that must be added [site units]
  */
  void dWaterchange(double dwaterchangeval) { waterchange_ += dwaterchangeval; }

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
  vector<unsigned int> writeSubVolume(string fileName, Site *centerste,
                                      int size);

  /**
  @brief Assign isotropic expansion strain at a set of prescribed sites.

  This function changes the strain components of a site already in the
  list of expansion sites.  If the prescribed site is not already in the
  list of expansion sites, then the site will be added to that list.

  @todo Consider changing the name of this method to applyExpansion

  @param alnb is the collection of site indices to which strain will be assigned
  @param exp is the isotropic expansion strain to set
  */
  void applyExp(vector<unsigned int> alnb, double exp);

  /**
  @brief Estimate the <i>internal</i> surface area of a phase with the aqueous
  solution.

  @param phaseid is the id of the microstructure phase
  @return the estimated surface area [site face units]
  */
  double getSurfaceArea(int phaseid);

  /**
  @brief Get the sorted distribution of domain sizes

  @param phaseid is the id of the phase to query
  @param numsites is the maximum number of sites to store and sort
  @param maxisze is the maxmimum linear size of interest
  @param sortorder is 0 if sorting in descending order, nonzero otherwise
  @return an STL list of the site ids according to the distribution
  */
  vector<int> findDomainSizeDistribution(int phaseid, const int numsites,
                                            int maxsize, int cyc, int sortorder);

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
  void setCount(vector<int> vect) {
    int dim = vect.size();
    for (int i = 0; i < dim; i++) {
      count_[i] = vect[i];
    }
  }
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

  void createRNG(void) { rg_ = new RanGen(); }
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



  void resetRNG(long int val_0, long int valLONGMAX, double valRNG, int cyc, int whileCount) {
    //latticeRNGseed_ = seed;
    rg_->setSeed(latticeRNGseed_);
    numRNGcall_0_ = val_0;
    numRNGcallLONGMAX_ = valLONGMAX;
    long int count_0 = 0, count_1 = 0;
    long int j0, j1, j11;
    double lastRNGreset;
    for(j1 = 1; j1 <= numRNGcallLONGMAX_; j1++) {
      for(j11 = 1; j11 <= LONG_MAX; j11++) {
        lastRNGreset = rg_->Ran3();
      }
    }
    for(j0 = 1; j0 <= val_0; j0++) {
      lastRNGreset = rg_->Ran3();
    }
    lastRNG_ = lastRNGreset;

    cout << endl << "Lattice::resetRNG cyc/whileCount/latticeRNGseed_: " << cyc << " / "
         << whileCount << " / " << latticeRNGseed_ << endl;
    cout << "Lattice::resetRNG numRNGcall_0_/numRNGcallLONGMAX_/lastRNGreset/valRNG: "
         << numRNGcall_0_ << " / " << numRNGcallLONGMAX_ << " / " << lastRNGreset
         << " / " << valRNG << endl;
    if (abs(lastRNGreset - valRNG) <= 1.e-16 ) {
      cout << "Lattice::resetRNG OK!" << endl;
    } else {
      cout << endl << "Lattice::resetRNG FAILED => exit" << endl;
      exit(0);
    }
  }

  void shiftAffinityPosVal(void){
    //check affinities values
    int minAff = 100000;
    for (int i = 0; i < numMicroPhases_; i++) {
      for (int j = 0; j < numMicroPhases_; j++) {
        if (chemSys_->getAffinity(i, j) <= minAff) minAff =chemSys_->getAffinity(i, j);
      }
    }
    cout << endl << "   Lattice::shiftAffinityPosVal minAff = " << minAff << endl;
    if (minAff < 0) {
      int newAff;
      int minAff_abs = abs(minAff);
      for (int i = 0; i < numMicroPhases_; i++) {
        for (int j = 0; j < numMicroPhases_; j++) {
          newAff = chemSys_->getAffinity(i, j) + minAff_abs;
          chemSys_->setAffinity(i, j, newAff);
        }
      }
      cout << "   Lattice::shiftAffinityPosVal => all affinities have been shifted with abs(minAff) = " << minAff_abs << endl;
    } else {
      cout << "   Lattice::shiftAffinityPosVal => all affinities are positive!" << endl;
    }
  }

  void increaseLatticeVolume(void);

  void checkSite(int stId){
    int phId = site_[stId].getMicroPhaseId();
    cout << endl << " Lattice::checkSite( " << stId << " ):" << endl;
    cout << "    phaseId : " << site_[stId].getMicroPhaseId() << endl;
    cout << "    inDissInterfacePos_ : " << site_[stId].getInDissInterfacePos() << endl;
    if (site_[stId].getInDissInterfacePos() != -1) {
      cout << "     in dissInterface on pos inDissInterfacePos_ : "
           << interface_[site_[stId].getMicroPhaseId()].getDissolutionSitesId(site_[stId].getInDissInterfacePos()) << endl;
    }

    vector<unsigned int> growth = site_[stId].getGrowthPhases();
    int size = growth.size();
    cout << endl << " growth_.size() : " << size<< endl;
    for (int k =0; k < size; k++) {
      cout << "       k = " << k << "   growth_[k] = " << growth[k] << endl;
    }
    cout << endl << " inGrowInterfacePos_ : " << endl;
    for (int k =0; k < numMicroPhases_; k++) {
      cout << "       k = " << k << "   site_[stId].getInGrowInterfacePos(k) = " << site_[stId].getInGrowInterfacePos(k) << endl;
    }

    cout << endl << "     in growInterfaces on pos inGrowInterfacePos_ :" << endl;
    for (int k = 0; k < numMicroPhases_; k++){
      cout << "       k_ = " << k << endl; cout.flush();
      if (site_[stId].getInGrowInterfacePos(k) > -1) {
        size = growthInterfaceSize_[k];
        cout << "       k = " << k << "   pos = " << site_[stId].getInGrowInterfacePos(k) << "   size = " << size << endl; cout.flush();
        cout << "            siteId in grInt = "
             << interface_[k].getGrowthSitesId(site_[stId].getInGrowInterfacePos(k)) << endl; cout.flush();
      }
    }
  }

  vector<int> getGrowthInterfaceSize(void) {
    return growthInterfaceSize_;
  }
  void setGrowthInterfaceSize(vector<int> vect){
    growthInterfaceSize_ = vect;
  }
  vector<int> getDissolutionInterfaceSize(void) {
    return dissolutionInterfaceSize_;
  }
  void setDissolutionInterfaceSize(vector<int> vect) {
    dissolutionInterfaceSize_ = vect;
  }


vector<int> chooseNucleationSitesRND(int phaseID,int numLeft);
vector<int> chooseNucleationSitesAFF(int phaseID,int numLeft);

}; // End of Lattice class
#endif

///
/// The functions below are used to aid in comparison of one site to another, by
/// which means lists of the sites can be sorted.
///

#ifndef CMPFUNCS
#define CMPFUNCS

/**
@brief Compare two sites, returning true is the first site is "less than" the
second.

The comparison is made on the basis of what THAMES loosely calls the
<i>weighted mean curvature</i>, (wmc).  A site with high wmc is a site where
dissolution of a phase is likely to occur, and growth of another phase is
unlikely to occur. Conversely, a site with a low wmc is a site where growth of a
phase is likely to occur but dissolution of a phase is unlikely to occur.

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

#endif
