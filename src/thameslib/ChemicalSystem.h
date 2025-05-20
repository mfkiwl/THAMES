/**
@file ChemicalSystem.h
@brief Declaration of the ChemicalSystem base class.

A const `ChemicalSystem` object contains the data
relevant to the properties of every phase known to THAMES,
such as composition, growth characteristics, molar volume, specific
gravity, heat capacity, etc.
This class divides matter into three major types:
   -# <b>Independent Components</b> (ICs), which are basically chemical
elements.
   -# <b>Dependent Components</b> (DCs), which are composed of one or more ICs.
       DCs include ion complexes, pure condensed phases, and pure vapors.
   -# <b>Phases</b>, which are composed of one or more DCs.

The operation of THAMES depends heavily upon the design of the
`ChemicalSystem` class because the phase in each cell governs
the types of chemical reactions that can occur.  Therefore, careful attention
must be paid to making this class as useful and, at the same time, easy to use
as possible.
*/

#ifndef CHEMSYSH
#define CHEMSYSH

#include "../GEMS3K-standalone/GEMS3K/node.h"
#include "../Resources/include/nlohmann/json.hpp"
#include "global.h"
#include "utils.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <typeinfo>
#include <vector>

using namespace std;
using json = nlohmann::json;

#ifndef CHEMSYSDATASTRUCT
#define CHEMSYSDATASTRUCT

/**
@struct PoreSizeVolume
@brief Volume fraction of a sub-voxel pore of a giveen effective diameter
*/

struct PoreSizeVolume {
  double diam;
  double volume;
  double volfrac;
};

/**
@struct PhaseData
@brief Stores data about each phase possible in the system for ease of parsing
the input files.

In THAMES, phases are identified either thermodynamically--- in the
GEM data repository--- or microstructurally.  A microstructural phase can
be one, or a combination of more than on, thermodynamically defined phase.

The structure is defined to make it easier to parse the input file for the GEM
chemical system definitions (CSD). It is not used elsewhere in the code.  In
fact, the same members are identified as class variables in the `ChemicalSystem`
class.

Most of the members have self-evident meanings:
    - `id` is the unique integer id of the microstructure phase
    - `randomGrowth` determines how much randomness is associated with the
       growth of the phase, rather than being determined only by mean curvature
    - `stressCalc` determines whether or not crystallization pressure should be
calculated
    - `weak` determines whether or not the phase can be damaged by stress
    - `k2o`, `na2o`, `mgo`, and `so3` are the mass fractions of potassium,
       sodium, magnesium, and sulfur oxides dissolved within the phase.
    - 'RdId' is a vector of all the ICs that can be incorporated into
         growing phases
    - 'RdVal' is a vector of maximum mass fractions of ICs that can be
incorporated into growing phases
    - `porosity` is the volume fraction of internal porosity in the phase,
       (e.g., C-S-H)
    - `red`, `green`, and `blue` are the rgb values for the color that the
       phase will have in simulated color micrographs
    - `gray` is the grayscale index the phase will have in simulated
backscattered electron micrographs
    - `thamesName` is the name the phase will have in the microstructure
    - `GEMPhaseName` is a vector of the GEM (thermodynamic) phases making up the
       THAMES phase, given by name
    - `DCName` is a vector of the GEM DCs (dependent components) making up the
       THAMES phase, given by name
    - `GEMPhaseDCMembers` is a vector of the GEM (thermodynamic) phase DC
       structures making up the THAMES phase, indexed by their GEM DC id number
    - `GEMPhaseId` is a vector of the GEM (thermodynamic) phases making up
       the THAMES phase, given by their GEM phase id number
    - `DCId` is the id of the dependent component in the GEM data base
    - `DCPorosity` is subvoxel porosity of the dependent components
    - `growthTemplate` is a vector of the growth templates for the phase
    - `affinity` the vector of affinities of a phase for another phase
    - `colors` is the rgb colors for displaying the phase in visualizations
*/

struct PhaseData {
  int id;
  // double randomGrowth;
  int stressCalc;
  int weak;
  double k2o;
  double na2o;
  double mgo;
  double so3;
  vector<struct PoreSizeVolume> poreSizeDist;
  int red;
  int green;
  int blue;
  int gray;
  vector<int> colors;
  string thamesName;
  vector<char *> GEMPhaseName, DCName;
  vector<int> GEMPhaseDCMembers;
  vector<int> GEMPhaseId;
  vector<double> microPhaseDCPorosities;
  vector<int> DCId;
  vector<int> growthTemplate;
  // vector<int> affinity;
  vector<double> affinity;
  vector<double> contactAngle;

  //  vector<int>
  //      RdId; /**< Vector of IC ids of the partitioned components in the phase
  //      */
  //  vector<double> RdVal; /**< Vector of Rd values for each IC */
};
#endif

struct elemColor {
  int colorId;
  string altName;
  vector<int> rgb;
  vector<float> rgbf;
  int gray;
  float grayf;
};

struct elMod {
  float K;
  float G;
  float E;
  float n;
  string altName;
};

/**
@class ChemicalSystem
@brief Handles the tracking of phases and communications between GEM and THAMES.

A `ChemicalSystem` object contains the data relevant to the properties of
every phase known to a given run of THAMES, such as composition, growth
characteristics, molar volume, specific gravity, heat capacity, etc.

@section intro_sec Materials
The class divides matter into three major types:
    -# <b>Independent Components</b> (ICs), which are basically chemical
        elements or a unit of charge
    -# <b>Dependent Components</b> (DCs), which are composed of one or
        more ICs.  DCs include ion complexes, pure condensed phases, and
        pure vapors
    -# <b>Phases</b>, which are composed of one or more DCs.

@section Phases
One of the biggest challenges with ChemicalSystem is that it is an interface
between the thermodynamic representation of phases, stored in the associated GEM
Chemical System Definition (CSD), and the phases as they are represented in the
3D microstructure. To accommodate this, ChemicalSystem keeps track of two
entirely different lists of phases, and keeps track of how the members of the
list are associated with each other. For example, the THAMES microstructure
includes a phase called CSHQ, which represents the calcium silicate hydrate
product (C-S-H).  However, the GEM CSD defines several different phases, which
are different compositional and structural end members of a non-ideal solid
solution that makes up C-S-H.  Therefore, THAMES must keep a list of all the GEM
phases that are collectively called CSHQ.  The same thing is true for certain of
the AFm and AFt phases.

@section methods Methods
THAMES communicates frequently with the GEM-IPM library, and so a lot of its
data and methods are devoted to storing and communicating GEM data, and
translating that GEM data into THAMES data that are used to create and change
microstructure.

*/

class ChemicalSystem {

  bool jsonFormat_; /**< True if input GEM data files are JSON format */

  int numMicroPhases_;     /**< Total number of material components that
                                       the microstructure can contain.  The GEM
                                         chemical system definition (CSD) will
                                         list all of the possible phases in the
                                         system, but we may want to combine two or
                                         more phases into a single microstructure
                                         component */
  int numMicroImpurities_; /**< Number of impurities, as oxides, dissolved
                               in the phases that will be tracked
                               Typically the value will be 4 (potassium,
                               sodium, magnesium, and sulfur oxides) */
  int numICs_;             /**< Number of independent components (IC) */
  int numDCs_;             /**< Number of dependent components (DC) */
  int numGEMPhases_;       /**< Number of GEM phases in the CSD */
  int numSolutionPhases_;  /**< Number of GEM solution phases in the CSD;
                                solution phases are non-stoichiometric */
  vector<string> microPhaseName_; /**< Names of phases identified in a THAMES
                                        microstructure */
  vector<string>
      stressPhaseName_; /**< Names of phases that can have crystallization
                                pressure in the microstructure */
  vector<string> weakPhaseName_; /**< Names of solid phases that can be damaged
                                        by stress in the microstructure */
  vector<string> porousPhaseName_; /**< Names of solid phases that have internal
                                          porosity in the microstructure */
  vector<int> stressPhaseId_;   /**< IDs of phases that can have crystallization
                                        pressure in the microstructure */
  vector<int> weakPhaseId_;     /**< IDs of solid phases that can be damaged
                                       by stress in the microstructure */
  vector<int> porousPhaseId_;   /**< IDs of solid phases that have internal
                                       porosity in the microstructure */
  vector<string> ICName_;       /**< Names of ICs in the GEM CSD */
  vector<string> DCName_;       /**< Names of DCs in the GEM CSD */
  vector<string> GEMPhaseName_; /**< Names of phases in the GEM CSD */
  vector<int> microPhaseId_; /**< Unique ids of THAMES microstructure phases */

  // vector<vector<int>>
  //     RdICId_; /**< List of ICs that can be an impurity in each phase */
  // vector<vector<double>> Rd_;  /**< Rd values for each IC in each phase */
  vector<double> ICMolarMass_; /**< One molar mass for each IC [g/mm3ol] */
  vector<double> DCMolarMass_; /**< One molar mass for each DC [g/mol] */
  vector<double>
      GEMPhaseMolarMass_; /**< One molar mass for each GEM phase [g/mol] */
  vector<vector<int>>
      growthTemplate_; /**< A list of the phases on which a given phase
                               is allowed to grow; one list for each phase */
  // vector<vector<int>> affinity_;
  /**< A list of the microstructure phases with
   which a given microstructure phase has an
   affinity to associate when growing */
  vector<vector<double>> affinity_;
  vector<vector<double>> contactAngle_;

  map<int, vector<int>>
      microPhaseMembers_; /**< A list of all the CSD phase ids that are
                                associated with a given microstructure phase */

  map<int, vector<int>>
      microPhaseDCMembers_; /**< A list of all the CSD DC ids that are
                            associated with a given microstructure phase */
  map<int, vector<double>>
      microPhaseDCPorosities_; /**< A list of all the CSD DC porosities that are
                           associated with a given microstructure phase */
  map<int, vector<int>>
      GEMPhaseDCMembers_; /**< A list of all the CSD DC ids that are
                              associated with a given CSD phase */

  vector<bool> isKinetic_;   /**< Whether of not each phase is kinetically
                                   controlled */
  vector<bool> isDCKinetic_; /**< Whether of not each DC is kinetically
                                 controlled */
  /**
  @brief Initial solution composition

  This is a map of key-value pairs.  The key is the integer value of a
  dependent component (DC), and the value is the concentration of that DC in
  molal units [mol/kgw].
  */
  map<int, double> initialSolutionComposition_;

  /**
  @brief Fixed solution composition

  This is a map of key-value pairs.  The key is the integer value of a
  dependent component (DC), and the value is the concentration of that DC in
  molal units [mol/kgw].
  */
  map<int, double> fixedSolutionComposition_;

  /**
  @brief Attack solution composition

  This is a map of key-value pairs.  The key is the integer value of a
  dependent component (DC), and the value is the concentration of that DC in
  molal units [mol/kgw].
  */
  map<int, double> attackSolutionComposition_;

  double gasSolidRatio_; /**< mass ratio of gas to solids */

  /**
  @brief Initial gas composition

  This is a map of key-value pairs.  The key is the integer value of a
  dependent component (DC), and the value is the concentration of that DC in
  molal units [mol/kgw].
  */
  map<int, double> initialGasComposition_;

  /**
  @brief Fixed gas composition

  This is a map of key-value pairs.  The key is the integer value of a
  dependent component (DC), and the value is the concentration of that DC in
  molal units [mol/kgw].
  */
  map<int, double> fixedGasComposition_;

  /**
  @brief Attack gas composition

  This is a map of key-value pairs.  The key is the integer value of a
  dependent component (DC), and the value is the concentration of that DC in
  molal units [mol/kgw].
  */
  map<int, double> attackGasComposition_;

  /**
  @brief Volume fraction of each GEM CSD phase associated with a THAMES phase.

  @warning This variable might not be used
  */
  map<int, vector<double>> microPhaseMemberVolumeFraction_;

  vector<double> microPhasePorosity_; /**< The sub-voxel porosity of a given
                                  phase, such as C-S-H (dimensionless) */

  /**
  @brief Sub-voxel pore size distribution (volume basis) of each phase
  */
  vector<vector<struct PoreSizeVolume>> poreSizeDistribution_;

  vector<double> k2o_;        /**< Mass fraction of K<sub>2</sub>O dissolved in
                                    each phase, in units of
                                    g per 100 g of the phase */
  vector<double> na2o_;       /**< Mass fraction of Na<sub>2</sub>O dissolved in
                                    each phase, in units of
                                    g per 100 g of the phase */
  vector<double> mgo_;        /**< Mass fraction of MgO dissolved in
                                    each phase, in units of
                                    g per 100 g of the phase */
  vector<double> so3_;        /**< Mass fraction of SO<sub>3</sub> dissolved in
                                    each phase, in units of
                                    g per 100 g of the phase */
  vector<int> grayscale_;     /**< A number on [0,255] giving the relative
                                   grayscale brightness of the THAMES
                                   phases in a backscattered electron image */
  vector<vector<int>> color_; /**< A list of <r,g,b> values specifying the
                                    color of the THAMES phases in a false
                                    color micrograph */

  map<string, elemColor> colorN_;

  map<string, int> microPhaseIdLookup_; /**< Map that returns the vector index
                                         of the microstructure phase name */
  map<string, int> ICIdLookup_; /**< Map that returns the vector index of the
                                        IC name */
  map<string, int> DCIdLookup_; /**< Map that returns the vector index of the
                                        DC name */
  map<string, int> GEMPhaseIdLookup_; /**< Map that returns the vector index of
                                         the GEM CSD phase name */

  map<int, vector<int>>
      microPhaseToGEMPhase_; /**< Map that returns the GEM CSD phase for
                              a given microstructure phase */
  vector<vector<double>>
      DCStoich_; /**< List of amount of moles of each IC in a DC */

  vector<double> DCCharge_; /**< Charge associated with a DC */

  double *pGEMPhaseStoich_; /**< List of amount of moles of each IC in a
                                GEM CSD phase (pointer form) */
  /**
  @brief Solid stoichiometry list for communicating with GEM-IPM.

  @warning Only passed back and forth to GEM library
  */
  double *pSolidStoich_;

  /**
  @brief Solution stoichiometry list for communicating with GEM-IPM.

  @warning Only passed back and forth to GEM library
  */
  double *pSolutPhaseStoich_;

  /**
  @brief Solution solid (?) stoichiometry list for communicating with GEM-IPM.

  @warning Only passed back and forth to GEM library
  */
  double *pSolutSolidStoich_;

  vector<vector<double>>
      GEMPhaseStoich_; /**< List of amount of moles of each IC in
                           a given GEM CSD phase (vector form) */

  double *ICMoles_;     /**< List of number of moles of each IC in system */
  double *ICResiduals_; /**< List of errors in IC moles for mass balance */
  double *ICChemicalPotential_; /**< List of chemical potentials of each IC, in
                              the GEM dual solution */
  double *DCMoles_;             /**< List of moles of each DC */
  double *DCActivityCoeff_; /**< List of activity coefficients for each DC */

  double *GEMPhaseMoles_;     /**< List of moles of each phase in the system */
  double *GEMPhaseMass_;      /**< List of mass of each phase in the system */
  double *GEMPhaseVolume_;    /**< List of volume of each phase in the system */
  double *solutPhaseMoles_;   /**< List of moles of each solution phase in
                                      the system */
  double *solutPhaseMass_;    /**< List of mass of each solution phase in
                                       the system */
  double *DCLowerLimitss_;    /**< List of mass of each solution phase in
                                      the system */
  double *solutPhaseVolume_;  /**< List of volume of each phase in the system */
  double *prevGEMPhaseMoles_; /**< List of moles of each phase in the system
                                      in the previous time step */
  double *prevGEMPhaseMass_;  /**< List of mass of each phase in the system
                                in the previous time step */
  double *prevGEMPhaseVolume_;  /**< List of volume of each phase in the system
                                  in the previous time step */
  double *carrier_;             /**< List of moles of carrier (solvent) in
                                        multicomponent asymmetric phases */
  double *specificSurfaceArea_; /**< List of specific surface area of each
                                phase, in m<sup>2</sup>/kg */
  double *DCUpperLimit_;        /**< List of upper bound on moles of each DC */
  double *DCLowerLimit_;        /**< List of lower bound on moles of each DC,
                                        generally non-zero for numerical
                                        stability of the thermodynamic
                                        calculations */

  double *DCH0_; /**< List of molar enthalpies of DCs */

  /**
  @brief List of one-letter class codes for each IC, specifying its kind.

  Each independent component is exactly one of six kinds in GEM database:

      - e = any chemical element except hydrogen
      - h = hydrogen
      - a = component with unspecified stoichiometry, such as Nit
      - i = isotope of a chemical element
      - f = formula unit
      - z = electric charge
  */
  vector<char> ICClassCode_;

  /**
  @brief List of one-letter class codes for each DC, specifying its kind.

  Each dependent component is exactly one of 15 kinds in GEM database:

      - O = a single-component phase
      - I = ideal end member of a solid solution phase (Raoult)
      - J = minor component in a solid solution (Henry)
      - M = majority component in a solid solution (Raoult)
      - S = any aqueous species except for H<sup>+</sup> or electrons
      - T = H<sup>+</sup>
      - E = electron
      - K = surface complex represented as an aqueous species
      - W = water, H<sub>2</sub>O
      - L = other components of a solvent, such as alcohol
      - G = gas
      - V = steam
      - C = carbon dioxide, CO<sub>2</sub>
      - H = hydrogen gas, H<sub>2</sub>
      - N = nitrogen gas, N<sub>2</sub>
  */
  vector<char> DCClassCode_;

  /**
  @brief List of one-letter class codes for each phase defined in GEM CSD,
  specifying its kind.

  Each phase in the GEM CSD is exactly one of seven kinds:

      - a = aqueous electrolyte
      - g = mixture of gases
      - f = supercritical fluid
      - l = non-electrolyte liquid (melt)
      - x = dispersed solid with ion exchange in an aqueous system
      - s = condensed solid solution phase
      - d = dispersed multicomponent solid phase
  */
  vector<char> GEMPhaseClassCode_;

  double T_;             /**< System-wide temperature [K] */
  double P_;             /**< System-wide pressure [Pa] */
  double Vs_;            /**< System total volume [m<sup>3</sup>] */
  double Ms_;            /**< System total mass [kg] */
  double Gs_;            /**< System total Gibbs energy [J] */
  double Hs_;            /**< System total enthalpy [J] */
  double ionicStrength_; /**< Solution ionic strength [mol/kgw] */
  double pH_;            /**< Solution pH */
  double pe_;            /**< Solution pe */
  double Eh_;            /**< Solution Eh [volts] */

  TNode *node_; /**< Pointer to a GEM3K TNode object */

  long int nodeHandle_; /**< integer flag used to identify a node */
  long int nodeStatus_; /**< integer flag used to identify node's status */
  long int iterDone_;   /**< number of iterations performed in the most
                                recent GEM calculation on the node */
  int timesGEMWarned_;  /**< tracks total number of times that the
                                GEM_run calculation warned  */
  int timesGEMFailed_;  /**< tracks number of times in a row that the
                                GEM_run calculation failed */
  int maxGEMFails_;     /**< maximum number of times GEM_run is allowed
                                to fail without throwing an exception */
  /**
  @brief Whether or not the system is saturated with moisture.

  If this variable is nonzero, then the porosity imbibes water from an external
  reservoir as it is consumed by reactions.  If it is set to zero,then the water
  is not replaced and the capillary porosity begins to desiccate.
  */
  bool isSaturated_;

  /**
  @brief Time to begin sulfate solution exposure, in hours.

  THAMES enables one to turn off the hydration simulation (i.e., kinetic
  dissolution of phases according to empirical rate laws in an otherwise closed
  system) and begin a simulation of sulfate attack on that hydrated system. This
  variable designates the time (in hours) at which this switch should happen.
  */
  double sulfateAttackTime_;

  /**
  @brief Time to begin leaching simulation, in hours.

  THAMES enables one to turn of the hydration simulation (i.e., kinetic
  dissolution of phases according to empirical rate laws in an otherwise closed
  system) and begin a simulation of leaching by a low-pH solution.  This
  variable designates the time (in hours) at which this switch should happen.
  */
  double leachTime_;

  vector<double>
      microPhaseMass_; /**< Absolute mass of each microstructure phase */
  vector<double> microPhaseMassDissolved_; /**< Absolute mass dissolved of each
                                              microstructure phase */
  vector<double>
      microPhaseVolume_; /**< Absolute volume of each microstructure phase */

  double microVolume_; /**< Absolute volume of the microstructure */
  // double newMicroVolume_;  /**< Absolute volume of the microstructure */
  double initMicroVolume_; /**< Initial absolute volume of the microstructure */
  double
      microVoidVolume_; /**< Absolute volume of void space in microstrucxture */

  /**
  @brief Saturation index of each phase in the GEM CSD.

  The departure from equilibrium between a given solid phase and an aqueous
  solution is characterized by the saturation index, `SI_`, which is defined as
  the activity product for the dissolution reaction divided by the equilibrium
  value of that activity product.  This variable stores the current SI for each
  solid phase in the GEM CSD.
  */
  vector<double> SI_;
  vector<double> microPhaseSI_;

  bool verbose_; /**< Whether to produce verbose output */
  bool warning_; /**< Whether to produce warning output */

  vector<int>
      DC_to_MPhID_; // microPhaseId for a given DCId: dim(DC_to_MPhID_)=numDCs &
                    // initialVal(DC_to_MPhID_)=-1
  vector<bool> cementComponent_;
  // vector<bool> GEMPhaseBelongsToCement_;
  double initScaledCementMass_;
  double scaledCementMass_;
  vector<int> isParrotKilloh_;
  int isParrotKillohSize_;

  int waterDCId_; /**< coresp to DCName = "H2O@" */
  double waterMollarMass_;
  double waterMollarVol_;
  double beginAttackTime_;
  map<string, elMod> elasticModuli_;
  int aliteDCId_;
  int beliteDCId_;
  int aluminateDCId_;
  int ferriteDCId_;

  double initSolidMass_;

public:
  void setInitSolidMass(double val) { initSolidMass_ = val; }

  /**
  @brief Constructor.

  Only one constructor is provided, which initializes the ChemicalSystem with
  all the information read from the GEM input files.

  @param GEMfilename is the name of the file holding GEM input data
  @param jsonFileName is the name of the file containing information about
  how to relate GEM phases to microstructure phases
  @param verbose is true if producing verbose output
  @param warning is true if producing verbose output
  */
  ChemicalSystem(const string &GEMfilename, const string &jsonFileName,
                 const bool verbose, const bool warning = false);

  /**
  @brief Copy constructor.

  @param obj is the ChemicalSysteme object to copy
  */
  ChemicalSystem(const ChemicalSystem &obj);

  /**
  @brief Destructor (does nothing for now).
  */
  ~ChemicalSystem(void);

  /**
  @brief Determine whether GEM input data files are in JSON format

  @param masterFileName is the lst name containing the others
  @return true if JSON is indicated, false otherwise
  */
  bool isInputFormatJSON(const char *masterFileName);

  /**
  @brief Get the name of the three JSON input files for GEMS data

  @param masterFileName is the lst name containing the others
  @param dchName is the name of the DCH data file
  @param ipmName is the name of the IPM data file
  @param dbrName is the name of the DBR data file
  */
  void getJSONFiles(const char *masterFileName, string &dchName,
                    string &ipmName, string &dbrName);

  /**
  @brief Master function for parsing an input file in JSON format.

  @param docName is the name of the JSON document to parse
  */
  void parseDoc(const string &docName);

  /**
  @brief Parse input about the initial solution composition from a JSON
  document.

  The initial solution composition, if given, is parsed by this function.
  The composition will be held in an associative map of key value pairs:

  * key = integer id of a GEM independent component
  * value = molal concentration of that component in the initial solution
  [mol/kgw]

  @param cdi is an iterator over the JSON data
  */
  void parseSolutionComp(const json::iterator cdi);

  /**
  @brief Parse input about the gas phase composition from a JSON document.

  The gas composition, if given, is parsed by this function.
  The composition will be held in an associative map of key value pairs:

  * key = integer id of a GEM independent component
  * value = molal concentration of that component in the initial solution
  [mol/kgw]

  @param cdi is a JSON iterator pointing to the current location in json data
  object
  */
  void parseGasComp(const json::iterator cdi);

  /**
  @brief Scan a JSON object

  @param cdi is an iterator over the JSON data
  @param phaseids is a map associating phase names with id numbers
  */
  void parseMicroPhaseNames(const json::iterator cdi,
                            map<string, int> &phaseids);

  /**
  @brief Parse input about a microstructure phase from a JSON document.

  @param cdi is an iterator over the JSON data
  @param numEntries is the number of entries in the JSON file
  @param phaseids is a map associating phase names with id numbers
  @param phaseData holds the structure of collected phase data from the document
  */
  void parseMicroPhases(const json::iterator cdi, int numEntries,
                        map<string, int> phaseids, PhaseData &phaseData);

  /**
  @brief Parse input about a GEM CSD phase from a JSON document.

  @param p is an iterator over the JSON data
  @param phaseData holds the structure of collected phase data from the document
  */
  void parseGEMPhaseData(const json::iterator p, PhaseData &phaseData);

  /**
  @brief Parse input about a GEM DC associated with a CSD phase from a JSON
  document.

  @param pp is an iterator over the JSON data
  @param phaseData holds the structure of collected phase data from the document
  */
  void parseGEMPhaseDCData(const json::iterator pp, PhaseData &phaseData);

  /**
  @brief Parse a phase's sub-voxel pore size distribution

  The file must have a header line that will be discarded.  The rest of
  the file must be two column csv with pore diameter in the first column
  (nanometers) and the volume fraction in the second column.  After
  the data are read the volume fraction is normalized.

  @param p is an iterator over the JSON data array
  @param phaseData holds the structure of collected phase data from the document
  */
  void parsePoreSizeDistribution(const json::iterator p, PhaseData &phaseData);

  /**
  @brief Parse the Rd data (impurity partitioning) for one phase in the JSON
  input file.

  @param p is an iterator over the JSON data object
  @param phaseData is a reference to the PhaseData structure for temporarily
  storing the input parameters
  */
  //  void parseRdData(const json::iterator p, struct PhaseData &phaseData);

  /**
  @brief Parse input about how to render a phase in an image.

  @param p is an iterator over the JSON data object
  @param phaseData holds the structure of collected phase data from the document
  */
  void parseDisplayData(const json::iterator p, PhaseData &phaseData);

  /**
  @brief Parse input about dissolved impurities within a phase.

  @param p is an iterator over the JSON data object
  @param phaseData holds the structure of collected phase data from the document
  */
  void parseImpurityData(const json::iterator p, PhaseData &phaseData);

  /**
  @brief Parse input about interfaces associated with a phase.

  @param p is an iterator over the JSON data object
  @param phaseids is a map associating phase names with id numbers
  @param phaseData holds the structure of collected phase data from the document
  */
  void parseInterfaceData(const json::iterator p, map<string, int> &phaseids,
                          PhaseData &phaseData);

  /**
  @brief Parse input about affinity for one phase to grow on another.

  @param pp is an iterator over the JSON data object
  @param phaseids is a map associating phase names with id numbers
  @param phaseData holds the structure of collected phase data from the document
  */
  void parseAffinityData(const json::iterator pp, map<string, int> &phaseids,
                         PhaseData &phaseData);

  /**
  @brief Set the total number of possible microstructure phases in the system.

  @note NOT USED.

  @param val is the total number of possible microstructure phases
  (non-negative)
  */
  void setNumMicroPhases(const int val) { numMicroPhases_ = val; }

  /**
  @brief Get the total number of possible microstructure phases in the system.

  @return the total number of possible microstructure phases (non-negative)
  */
  int getNumMicroPhases(void) const { return numMicroPhases_; }

  /**
  @brief Set the number of possible dissolved impurities in the system.

  @note NOT USED.

  @param val is the total number of possible impurities (non-negative)
  */
  void setNumMicroImpurities(const int val) { numMicroImpurities_ = val; }

  /**
  @brief get the number of possible dissolved impurities in the system.

  @return the total number of possible impurities (non-negative)
  */
  int getNumMicroImpurities(void) const { return numMicroImpurities_; }

  /**
  @brief Set the number of independent components (ICs).

  @note NOT USED.

  @param val is the total number of independent components (non-negative)
  */
  void setNumICs(const int val) { numICs_ = val; }

  /**
  @brief Get the number of independent components (ICs).

  @return the number of independent components (non-negative)
  */
  int getNumICs(void) const { return numICs_; }

  /**
  @brief Set the number of dependent components (DCs).

  @note NOT USED.

  @param val is the total number of dependent components (non-negative)
  */
  void setNumDCs(const int val) { numDCs_ = val; }

  /**
  @brief Get the number of dependent components (DCs).

  @return the number of dependent components (non-negative)
  */
  int getNumDCs(void) const { return numDCs_; }

  /**
  @brief Set the number of phases in the GEM CSD.

  @note NOT USED.

  @param val is the total number of phases recognized in the GEM CSD
  (non-negative)
  */
  void setNumGEMPhases(const int val) { numGEMPhases_ = val; }

  /**
  @brief Get the number of phases in the GEM CSD.

  @return the number of phases in the GEM CSD (non-negative)
  */
  int getNumGEMPhases(void) const { return numGEMPhases_; }

  /**
  @brief Set the number of solution phases in the GEM CSD.

  @note NOT USED.

  @param val is the total number of solution phases recognized in the GEM CSD
  (non-negative)
  */
  void setNumSolutionPhases(const int val) { numSolutionPhases_ = val; }

  /**
  @brief Get the number of solution phases in the GEM CSD.

  @note Used only in this class's copy constructor.

  @return the number of solution phases in the GEM CSD (non-negative)
  */
  int getNumSolutionPhases(void) const { return numSolutionPhases_; }

  /**
  @brief Set the name of a microstructure phase.

  @note NOT USED.

  @param idx is the phase index number (non-negative)
  @param str is the name to assign to the phase
  */
  void setMicroPhaseName(const int idx, const string &str) {
    // try {
    microPhaseName_[idx] = str;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setMicroPhaseName",
    //  "microPhaseName_",
    //                  microPhaseName_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Get the name of a microstructure phase.

  @param idx is the phase index number (non-negative)
  @return the name of the phase
  */
  string &getMicroPhaseName(const unsigned int idx) {
    // try {
    return (string &)microPhaseName_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getMicroPhaseName",
    //  "microPhaseName_",
    //                  microPhaseName_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the list of all microstructure phase names.

  @return the vector of microstructure phase names
  */
  vector<string> getMicroPhaseName(void) const { return microPhaseName_; }

  /**
  @brief Set the name of a microstructure phase that can have crystallization
  pressure.

  @note NOT USED.

  @param idx is the phase index number (non-negative)
  @param str is the name to assign to the phase
  */
  void setStressPhaseName(const unsigned int idx, const string &str) {
    // try {
    stressPhaseName_[idx] = str;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setStressPhaseName",
    //                  "stressPhaseName_", stressPhaseName_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Get the name of a microstructure phase that can have crystallization
  pressure.

  @param idx is the phase index number (non-negative)
  @return the name of the phase
  */
  string &getStressPhaseName(const unsigned int idx) {
    // try {
    return (string &)stressPhaseName_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getStressPhaseName",
    //                  "stressPhaseName_", stressPhaseName_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the list of all microstructure phase names that can have
  crystallization pressure.

  @return the vector of microstructure phase names
  */
  vector<string> getStressPhaseName(void) const { return stressPhaseName_; }

  /**
  @brief Get the list of all microstructure ids that can have crystallization
  pressure

  @return the vector of stress phase ids
  */
  vector<int> getStressPhaseId_(void) const { return stressPhaseId_; }

  /**
  @brief Determine if a given microstructure phase is eligible
  for crystallization pressure.

  The determination is made solely on the basis of whether the phase is
  a member of the stressPhaseName_ vector.

  @param idx is the phase id to check
  @return true if the phase is subject to crystallization pressure
  */
  bool isStress(const int idx) {
    bool istress = false;
    int size = stressPhaseId_.size();
    for (int i = 0; i < size; ++i) {
      if (idx == stressPhaseId_[i]) {
        istress = true;
        break;
      }
    }
    return istress;
  }

  /**
  @brief Determine if a given microstructure phase is eligible
  for crystallization pressure

  The determination is made solely on the basis of whether the phase is
  a member of the stressPhaseName_ vector.

  @param str is the name to check
  @return true if the phase is subject to crystallization pressure
  */
  bool isStress(const string &str) {
    bool istress = false;
    int size = stressPhaseName_.size();
    for (int i = 0; i < size; ++i) {
      if (str == stressPhaseName_[i]) {
        istress = true;
        break;
      }
    }
    return istress;
  }
  /**
  @brief Set the name of a microstructure phase that can be damaged by stress

  @note NOT USED.

  @param idx is the phase index number (non-negative)
  @param str is the name to assign to the phase
  */
  void setWeakPhaseName(const unsigned int idx, const string &str) {
    // try {
    weakPhaseName_[idx] = str;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setWeakPhaseName", "weakPhaseName_",
    //                  weakPhaseName_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Get the name of a microstructure phase that can be damaged by stress

  @param idx is the phase index number (non-negative)
  @return the name of the phase
  */
  string &getWeakPhaseName(const unsigned int idx) {
    // try {
    return (string &)weakPhaseName_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getWeakPhaseName", "weakPhaseName_",
    //                  weakPhaseName_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the list of all microstructure phase names that can be damaged by
  stress

  @return the vector of microstructure phase names
  */
  vector<string> getWeakPhaseName(void) const { return weakPhaseName_; }

  /**
  @brief Get the list of all microstructure ids that can be damaged by stress

  @return the vector of microstructure phase ids
  */
  vector<int> getWeakPhaseId_(void) const { return weakPhaseId_; }

  /**
  @brief Determine if a given microstructure phase is eligible
  for damage by stress

  The determination is made solely on the basis of whether the phase is
  a member of the weakPhaseName_ vector.

  @param idx is the phase id to check
  @return true if the phase is subject to damage
  */
  bool isWeak(const int idx) {
    bool isweak = false;
    int size = weakPhaseId_.size();
    for (int i = 0; i < size; ++i) {
      if (idx == weakPhaseId_[i]) {
        isweak = true;
        break;
      }
    }
    return isweak;
  }

  /**
  @brief Determine if a given microstructure phase is eligible
  for damage by stress

  The determination is made solely on the basis of whether the phase is
  a member of the weakPhaseName_ vector.

  @param str is the name to check
  @return true if the phase is subject to damage
  */
  bool isWeak(const string &str) {
    bool isweak = false;
    int size = weakPhaseName_.size();
    for (int i = 0; i < size; ++i) {
      if (str == weakPhaseName_[i]) {
        isweak = true;
        break;
      }
    }
    return isweak;
  }

  /**
  @brief Set the name of a microstructure phase that has internal porosity.

  @note NOT USED.

  @param idx is the phase index number (non-negative)
  @param str is the name to assign to the phase
  */
  void setPorousPhaseName(const unsigned int idx, const string &str) {
    // try {
    porousPhaseName_[idx] = str;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setPorousPhaseName",
    //                  "porousPhaseName_", porousPhaseName_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Get the name of a microstructure phase that has internal porosity

  @param idx is the phase index number (non-negative)
  @return the name of the phase
  */
  string &getPorousPhaseName(const unsigned int idx) {
    // try {
    return (string &)porousPhaseName_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getPorousPhaseName",
    //                  "porousPhaseName_", porousPhaseName_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the list of all microstructure phase names that have internal
  porosity

  @return the vector of microstructure porous phase names
  */
  vector<string> getPorousPhaseName(void) const { return porousPhaseName_; }

  /**
  @brief Get the list of all microstructure ids with internal porosity

  @return the vector of porous phase ids
  */
  vector<int> getPorousPhaseId_(void) const { return porousPhaseId_; }

  /**
  @brief Determine if a given microstructure phase is porous

  The determination is made solely on the basis of whether the phase is
  a member of the porousPhaseName_ vector.

  @param str is the name to check
  @return true if the phase is porous
  */
  bool isPorous(const string &str) {
    bool isporous = false;
    int size = porousPhaseName_.size();
    for (int i = 0; i < size; ++i) {
      if (str == porousPhaseName_[i]) {
        isporous = true;
        break;
      }
    }
    return isporous;
  }

  /**
  @brief Determine if a given microstructure phase is porous
  The determination is made solely on the basis of whether the phase is
  a member of the porousPhaseName_ vector.

  @param idx is the phase id to check
  @return true if the phase is porous
  */
  bool isPorous(const int idx) {
    bool isporous = false;
    int size = porousPhaseId_.size();
    for (int i = 0; i < size; ++i) {
      if (idx == porousPhaseId_[i]) {
        isporous = true;
        break;
      }
    }
    return isporous;
  }

  /**
  @brief Set the name of an independent component (IC).

  @note NOT USED.

  @param idx is the IC index number (non-negative)
  @param str is the name to assign to the IC
  */
  void setICName(const unsigned int idx, const string &str) {
    // try {
    ICName_[idx] = str;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setICName", "ICName_",
    //  ICName_.size(),
    //                  idx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Get the name of an independent component (IC).

  @param idx is the IC index number (non-negative)
  @return the name of the IC
  */
  string &getICName(const unsigned int idx) {
    // try {
    return (string &)ICName_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getICName", "ICName_",
    //  ICName_.size(),
    //                  idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the list of all independent component (IC) names.

  @return the vector of IC names
  */
  vector<string> getICName(void) const { return ICName_; }

  /**
  @brief Set the name of a dependent component (DC).

  @note NOT USED.

  @param idx is the DC index number (non-negative)
  @param str is the name to assign to the DC
  */
  void setDCName(const unsigned int idx, const string &str) {
    // try {
    DCName_[idx] = str;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setDCName", "DCName_",
    //  DCName_.size(),
    //                  idx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Get the name of a dependent component (DC).

  @param idx is the DC index number (non-negative)
  @return the name of the DC
  */
  string &getDCName(const unsigned int idx) {
    // try {
    return (string &)DCName_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getDCName", "DCName_",
    //  DCName_.size(),
    //                  idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the list of all dependent component (DC) names.

  @return the vector of DC names
  */
  vector<string> getDCName(void) const { return DCName_; }

  /**
  @brief Set the name of a phase in the GEM CSD.

  @note NOT USED.

  @param idx is the GEM phase index number (non-negative)
  @param str is the name to assign to the phase
  */
  void setGEMPhaseName(const unsigned int idx, const string &str) {
    // try {
    GEMPhaseName_[idx] = str;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setGEMPhaseName", "GEMPhaseName_",
    //                  GEMPhaseName_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Get the name of a phase in the GEM CSD.

  @param idx is the GEM phase index number (non-negative)
  @return the name of the GEM phase
  */
  string &getGEMPhaseName(const unsigned int idx) {
    // try {
    return (string &)GEMPhaseName_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getGEMPhaseName", "GEMPhaseName_",
    //                  GEMPhaseName_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the list of all GEM CSD phase names.

  @return the vector of GEM phase names
  */
  vector<string> getGEMPhaseName(void) const { return GEMPhaseName_; }

  /**
  @brief Set the integer id of a microstructure phase.

  @note NOT USED.

  @param idx is the vector element holding the phase id (non-negative)
  @param val is the non-negative integer id to assign to the phase
  */
  void setMicroPhaseId(const unsigned int idx, const int val) {
    // try {
    microPhaseId_[idx] = val;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setMicroPhaseId", "microPhaseId_",
    //                  microPhaseId_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Get the integer id of a microstructure phase by its position in the id
  vector.

  @note NOT USED.

  @param idx is the element of the vector holding the phase's id (non-negative)
  @return the integer id stored at that element
  */
  int getMicroPhaseId(const unsigned int idx) {
    // try {
    return microPhaseId_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getMicroPhaseId", "microPhaseId_",
    //                  microPhaseId_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the integer id of a microstructure phase by its name.

  @param micname is the name of the microstructure phase
  @return the integer id associated with that phase name
  */
  int getMicroPhaseId(const string &micname) {
    map<string, int>::iterator p = microPhaseIdLookup_.find(micname);
    if (p != microPhaseIdLookup_.end()) {
      return p->second;
    } else {
      string msg = "Could not find microPhaseId_ match to " + micname;
      EOBException ex("ChemicalSystem", "getMicroPhaseId", msg,
                      microPhaseId_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  int getMicroPhaseId_SA(const string &micname) {
    map<string, int>::iterator p = microPhaseIdLookup_.find(micname);
    if (p != microPhaseIdLookup_.end()) {
      return p->second;
    } else {
      return -1;
    }
  }

  /**
  @brief Get the list of all phase ids.

  @return the vector holding phase id numbers
  */
  vector<int> getMicroPhaseId(void) const { return microPhaseId_; }

  /**
  @brief Get the list of booleans indicating whether a microstructure
  phase is kinetically controlled

  @note Used only in the copy constructor

  @return the vector holding phase id numbers
  */
  vector<bool> getIsKinetic(void) const { return isKinetic_; }

  /**
  @brief Get the kinetic status of a microstructure phase by its id

  @param idx is a microstructure id number
  @return true if it is kinetically controlled
  */
  bool isKinetic(const unsigned int idx) {
    // try {
    return isKinetic_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "isKinetic", "isKinetic_",
    //                  isKinetic_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the kinetic status of a microstructure phase by its name

  @note NOT USED.

  @param micname is the name of a microstructure phase
  @return true if it is kinetically controlled
  */
  bool isKinetic(const string micname) {
    // try {
    int idx = getMicroPhaseId(micname);
    return isKinetic_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "isKinetic", "isKinetic_",
    //                  isKinetic_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Set the kinetic status of a microstructure phase

  @param idx is the index of a microstructure phase
  */
  void setIsKinetic(vector<bool> &isKinetic) {
    isKinetic_ = isKinetic;
    return;
  }

  /**
  @brief Find out if an IC exists in the system

  @param icname is the name of the proposed IC
  @return true if the IC is recognized
  */
  bool isIC(const string &icname) {
    string msg;
    map<string, int>::iterator p = ICIdLookup_.find(icname);
    if (p != ICIdLookup_.end()) {
      return true;
    }
    return false;
  }

  /**
  @brief Get the integer id of an independent component (IC) by its name.

  @param icname is the name of the IC
  @return the integer id associated with that IC name (non-negative)
  */
  int getICId(const string &icname) {
    map<string, int>::iterator p = ICIdLookup_.find(icname);
    if (p != ICIdLookup_.end()) {
      return p->second;
    } else {
      if (warning_) {
        cout << "WARNING: Could not find ICIdLookup_ match to " << icname
             << endl;
      }
      return (numICs_ + 9999); // nonsense number should be detected
    }
  }

  /**
  @brief Find out if an DC exists in the system

  @param dcname is the name of the proposed DC
  @return true if the DC is recognized
  */
  bool isDC(const string &dcname) {
    map<string, int>::iterator p = DCIdLookup_.find(dcname);
    if (p != DCIdLookup_.end()) {
      return true;
    }
    return false;
  }

  /**
  @brief Find out if named DC is solid

  @param dcname is the name of the DC
  @return true if the DC is solid
  */
  bool isDCSolid(const string &dcname) {
    if (isDC(dcname)) {
      char classcode = getDCClassCode(getDCId(dcname));
      if (classcode == 'I' || classcode == 'J' || classcode == 'M' ||
          classcode == 'O') {
        return true;
      }
    }
    return false;
  }

  /**
  @brief Get the integer id of a dependent component (DC) by its name.

  @param dcname is the name of the DC
  @return the integer id associated with that DC name
  */
  int getDCId(const string &dcname) {
    map<string, int>::iterator p = DCIdLookup_.find(dcname);
    if (p != DCIdLookup_.end()) {
      return p->second;
    } else {
      cout << "WARNING: Could not find DCIdLookup_ match to " << dcname << endl;
      // cout << "WARNING: Here are the ones I know about:" << endl;
      // cout.flush();
      // p = DCIdLookup_.begin();
      // while (p != DCIdLookup_.end()) {
      //     cout << "WARNING:     " << p->first << " ("
      //          << p->second << ")" << endl;
      //     cout.flush();
      //     p++;
      // }
      // cout << "WARNING:" << endl;
      // cout.flush();
      cout << "exit the program" << endl;
      exit(1);
      return (numDCs_ + 9999);
    }
  }

  /**
  @brief Find out if a GEM phase exists in the system

  @param phasename is the name of the proposed phase
  @return true if the GEM phase is recognized
  */
  bool isGEMPhase(const string &phasename) {
    map<string, int>::iterator p = GEMPhaseIdLookup_.find(phasename);
    if (p != GEMPhaseIdLookup_.end()) {
      return true;
    }
    return false;
  }

  /**
  @brief Get the integer id of a GEM CSD phase by its name.

  @param phasename is the name of the GEM phase
  @return the integer id associated with that GEM phase name
  */
  int getGEMPhaseId(const string &phasename) {
    map<string, int>::iterator p = GEMPhaseIdLookup_.find(phasename);
    if (p != GEMPhaseIdLookup_.end()) {
      return p->second;
    } else {
      if (warning_) {
        cout << "Could not find GEMPhaseIdLookup_ match to " << phasename
             << endl;
        cout << "exit the program" << endl;
        exit(1);
      }
      return (numGEMPhases_ + 9999);
    }
  }

  /**
  @brief Get the integer id of a GEM CSD phase by its THAMES id.

  @param thamesid is the index of the phase in the THAMES array
  @return the integer id associated with that GEM phase name
  */
  int getGEMPhaseId(const int thamesid) {
    string gemphasename = GEMPhaseName_[thamesid];
    return (getGEMPhaseId(gemphasename));
  }

  /**
  @brief Get the integer id of a GEM phase associated with a microstructure
  phase id.

  @param i is integer id of the microstructure phase
  @param idx is element position in the vector of associated GEM phases for that
          microstructure phase
  @return the integer id  of the GEM phase stored at position idx in the list
  */
  int getMicroPhaseToGEMPhase(const int i, const int idx) {
    map<int, vector<int>>::iterator p = microPhaseToGEMPhase_.find(i);
    if (p != microPhaseToGEMPhase_.end()) {
      if (idx < (int)((p->second).size())) {
        return (p->second)[idx];
      } else {
        string msg = "microPhaseToGEMPhase_";
        EOBException ex("ChemicalSystem", "getMicroPhaseToGEMPhase(1)", msg,
                        (p->second).size(), idx);
        ex.printException();
        exit(1);
      }
    } else {
      string msg =
          "Could not find microPhaseToGEMPhase_ match to index provided";
      EOBException ex("ChemicalSystem", "getMicroPhaseToGEMPhase(1)", msg,
                      microPhaseToGEMPhase_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the list of integer ids of GEM phases associated with a
  microstructure phase id.

  @param i is integer id of the microstructure phase
  @return the vector of integer ids of the GEM phases for that microstructure
  phase
  */
  vector<int> getMicroPhaseToGEMPhase(const int i) {
    map<int, vector<int>>::iterator p = microPhaseToGEMPhase_.find(i);
    if (p != microPhaseToGEMPhase_.end()) {
      return p->second;
    } else {
      string msg =
          "Could not find microPhaseToGEMPhase_ match to index provided";
      EOBException ex("ChemicalSystem", "getMicroPhaseToGEMPhase(2)", msg,
                      microPhaseToGEMPhase_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Set the integer id of a GEM phase associated with a microstructure
  phase id.

  @note NOT USED.

  @param microPhaseid is integer id of the microstructure phase
  @param GEMphaseid is the integer id of the GEM phase
  */
  void setMicroPhaseToGEMPhase(const int microPhaseid, const int GEMphaseid) {
    map<int, vector<int>>::iterator p =
        microPhaseToGEMPhase_.find(microPhaseid);
    if (p == microPhaseToGEMPhase_.end()) {
      vector<int> GEMphasevector;
      GEMphasevector.clear();
      GEMphasevector.push_back(GEMphaseid);
      microPhaseToGEMPhase_.insert(make_pair(microPhaseid, GEMphasevector));
    } else {
      bool found = false;
      int size = (p->second).size();
      for (int i = 0; i < size; i++) {
        if ((p->second)[i] == GEMphaseid) {
          found = true;
          break;
        }
      }
      if (!found)
        (p->second).push_back(GEMphaseid);
    }
  }

  /**
  @brief Get the list of integer ids of GEM phases associated with a
  microstructure phase name.

  @note NOT USED.

  @param microPhasename is the name of the microstructure phase
  @return the vector of integer ids of the GEM phases for that microstructure
  phase
  */
  vector<int> getMicroPhaseToGEMPhase(const string &microPhasename) {
    int i = (int)(getMicroPhaseId(microPhasename));
    return (getMicroPhaseToGEMPhase(i));
  }

  /**
  @brief Get the map relating every microstructure phase to a list of its
  associated GEM phases.

  @note Used only in this class's copy constructor.

  @return the map relating every microstructure phase to a list of its
  associated GEM phases
  */
  map<int, vector<int>> getMicroPhaseToGEMPhase(void) const {
    return microPhaseToGEMPhase_;
  }

  /**
  @brief Get the random growth tendency parameter for a microstructure phase.

  A given phase, whether hydration product or product of chemical degradation,
  will generally grow in a compact form by make the growth potential highest
  at sites with low mean curvature and lowest at point with high mean curvature.
  The growth sites are ordered from lowest to highest mean curvature when the
  random growth parameter is set to zero.  Higer values of random growth
  parameter cause more severe shuffling of the ordered list of growth sites.

  @param idx is the microstructure phase id
  @return the random growth parameter value for that phase
  */
  // double getRandomGrowth(const unsigned int idx) {
  //   try {
  //     return randomGrowth_[idx];
  //   } catch (out_of_range &oor) {
  //     EOBException ex("ChemicalSystem", "getRandomGrowth", "randomGrowth_",
  //                     randomGrowth_.size(), idx);
  //     ex.printException();
  //     exit(1);
  //   }
  // }

  /**
  @brief Get the whole random growth tendency parameters for all microstructure
  phases.

  A given phase, whether hydration product or product of chemical degradation,
  will generally grow in a compact form by make the growth potential highest
  at sites with low mean curvature and lowest at point with high mean curvature.
  The growth sites are ordered from lowest to highest mean curvature when the
  random growth parameter is set to zero.  Higer values of random growth
  parameter cause more severe shuffling of the ordered list of growth sites.

  @note Used only in this class's copy constructor.

  @return the vector of random growth parameters for all microstructure phases
  */
  // vector<double> getRandomGrowth(void) const { return randomGrowth_; }

  /**
  @brief Construct the growth template based on affinity values

  This function determines the templates of a phase.  A phase is considered
  a template only if the affinity value is greater than zero.

  @param idx is the vector of all affinities
  @return the ordered vector of those phases that have positive affinities
  */
  vector<int> calcGrowthTemplate(vector<double> affty) {
    vector<int> posaffty;
    int size = affty.size();
    // string msg;
    // try {
    for (int i = 0; i < size; ++i) {
      if (affty[i] > 0)
        posaffty.push_back(i);
    }
    return posaffty;
    //} catch (out_of_range &oor) {
    //  msg = "affty index is out of range";
    //  EOBException ex("ChemicalSystem", "calcGrowthTemplate", msg,
    //  affty.size(),
    //                  i);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the list of all microstructure phases on which a given
  microstructure phase can grow.

  Any phase can grow on its own surface, and we also allow phases to grow on the
  surface of other phases, mimicking heterogeneous nucleation and growth.
  This function seeks a requested template for growth of a phase.

  @param idx is the microstructure phase id
  @return the list of all templates for growth of that phase
  */
  vector<int> getGrowthTemplate(const int idx) {
    // string msg;
    // try {
    return growthTemplate_[idx];
    //} catch (out_of_range &oor) {
    //  msg = "growthTemplate_";
    //  EOBException ex("ChemicalSystem", "getGrowthTemplate", msg,
    //                  growthTemplate_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get one of the growth templates for a microstructure phase.

  Any phase can grow on its own surface, and we also allow phases to grow on the
  surface of other phases, mimicking hetergeneous nucleation and growth.
  This function seeks a requested template for growth of a phase, and sets the
  affinity of the growing phase for that template phase.  A high affinity for a
  template mimics a low energy barrier for heterogeneous nucleation.

  @note NOT USED.

  @param idx is the microstructure phase id
  @param jdx is the element to query in the list of templates for phase idx
  @return the integer id of the growth template phase at element jdx
  */
  int getGrowthTemplate(const int idx, const int jdx) {
    string msg;
    if (idx >= (int)growthTemplate_.size()) {
      msg = "growthTemplate_";
      EOBException ex("ChemicalSystem", "getGrowthTemplate", msg,
                      growthTemplate_.size(), idx);
      ex.printException();
      exit(1);
    }
    if (jdx >= (int)growthTemplate_[idx].size()) {
      msg = "growthTemplate_";
      EOBException ex("ChemicalSystem", "getGrowthTemplate", msg,
                      growthTemplate_[idx].size(), jdx);
      ex.printException();
      exit(1);
    }
    return growthTemplate_[idx][jdx];
  }

  /**
  @brief Get all the growth templates of all microstructure phases.

  Any phase can grow on its own surface, and we also allow phases to grow on the
  surface of other phases, mimicking hetergeneous nucleation and growth.
  This function seeks a requested template for growth of a phase, and sets the
  affinity of the growing phase for that template phase.  A high affinity for a
  template mimics a low energy barrier for heterogeneous nucleation.

  @note Used only in this class's copy constructor.

  @return the whole growth template list for all phases
  */
  vector<vector<int>> getGrowthTemplate(void) const { return growthTemplate_; }

  /**
  @brief Determine if a microstructure phase is a template for another phase's
  growth.

  Any phase can grow on its own surface, and we also allow phases to grow on the
  surface of other phases, mimicking hetergeneous nucleation and growth.
  This function seeks a requested template for growth of a phase, and sets the
  affinity of the growing phase for that template phase.  A high affinity for a
  template mimics a low energy barrier for heterogeneous nucleation.

  @param gphaseid is the integer id of the microstructure phase
  @param gtmpid is the id of another phase that may or may not be a template
  @return true if gtmpid is a template for growth of gphaseid, false otherwise
  */
  bool isGrowthTemplate(const int gphaseid, int gtmpid) {
    bool answer = false;
    int size = growthTemplate_[gphaseid].size();
    for (int i = 0; i < size; i++) {
      if (growthTemplate_[gphaseid][i] == gtmpid) {
        answer = true;
        break;
      }
    }
    return answer;
  }

  /**
  @brief Set the list of affinities for growth of a microstructure phase on its
  templates.

  A given phase, whether hydration product or product of chemical degradation,
  will generally grow in a compact form by make the growth potential highest
  at sites with low mean curvature and lowest at point with high mean curvature.
  The growth sites are ordered from lowest to highest mean curvature when the
  random growth parameter is set to zero.  Higher values of random growth
  parameter cause more severe shuffling of the ordered list of growth sites.

  @note NOT USED.

  @param idx is the microstructure phase id
  @param avec is the list of integer affinities for growth of the phase on its
  templates
  */
  void setAffinity(const int idx, vector<double> avec) {
    // string msg;
    // try {
    affinity_[idx] = avec;
    //} catch (out_of_range &oor) {
    //  msg = "affinity_";
    //  EOBException ex("ChemicalSystem", "setAffinity", msg, affinity_.size(),
    //                  idx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Set the affinity for growth of a microstructure phase on one of its
  templates.

  A given phase, whether hydration product or product of chemical degradation,
  will generally grow in a compact form by make the growth potential highest
  at sites with low mean curvature and lowest at point with high mean curvature.
  The growth sites are ordered from lowest to highest mean curvature when the
  random growth parameter is set to zero.  Higer values of random growth
  parameter cause more severe shuffling of the ordered list of growth sites.

  @note NOT USED.

  @param idx is the microstructure phase id
  @param jdx is the element to access in the list of growth affinities
  @param val is the integer affinity to assign to the phase at that element the
  list
  */
  void setAffinity(const int idx, const int jdx, const double val) {
    if (idx >= (int)affinity_.size()) {
      string msg = "affinity_";
      EOBException ex("ChemicalSystem", "setAffinity", msg, affinity_.size(),
                      idx);
      ex.printException();
      exit(1);
    }
    if (jdx >= (int)affinity_[idx].size()) {
      string msg = "affinity_";
      EOBException ex("ChemicalSystem", "setAffinity", msg,
                      affinity_[idx].size(), jdx);
      ex.printException();
      exit(1);
    }
    affinity_[idx][jdx] = val;
    return;
  }

  /**
  @brief Get the list of affinities for growth of a microstructure phase on its
  templates.

  A given phase, whether hydration product or product of chemical degradation,
  will generally grow in a compact form by make the growth potential highest
  at sites with low mean curvature and lowest at point with high mean curvature.
  The growth sites are ordered from lowest to highest mean curvature when the
  random growth parameter is set to zero.  Higer values of random growth
  parameter cause more severe shuffling of the ordered list of growth sites.

  @note NOT USED.

  @param idx is the microstructure phase id of which the affinities are sought
  @return the list of integer affinities for all the templates for growth of
  phase idx
  */
  vector<double> getAffinity(const int idx) {
    // try {
    return affinity_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getAffinity", "affinity_",
    //                  affinity_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the affinitiy for growth of a microstructure phase on one of its
  templates.

  A given phase, whether hydration product or product of chemical degradation,
  will generally grow in a compact form by making the growth potential highest
  at sites with low mean curvature and lowest at point with high mean curvature.
  The growth sites are ordered from lowest to highest mean curvature when the
  random growth parameter is set to zero.  Higher values of random growth
  parameter cause more severe shuffling of the ordered list of growth sites.

  @param idx is the microstructure phase id of which the affinity is sought
  @param jdx is the element in the list of affinities being queried
  @return the affinity for the template phase associated with list element jdx
  */
  double getAffinity(const int idx, const int jdx) {
    if (idx >= (int)affinity_.size()) {
      string msg = "affinity_";
      EOBException ex("ChemicalSystem", "getAffinity", msg, affinity_.size(),
                      idx);
      ex.printException();
      exit(1);
    }
    if (jdx >= (int)affinity_[idx].size()) {
      string msg = "affinity_";
      EOBException ex("ChemicalSystem", "getAffinity", msg,
                      affinity_[idx].size(), jdx);
      ex.printException();
      exit(1);
    }
    return affinity_[idx][jdx];
  }

  /**
  @brief Get the list of affinities for growth of every microstructure phase.

  A given phase, whether hydration product or product of chemical degradation,
  will generally grow in a compact form by make the growth potential highest
  at sites with low mean curvature and lowest at point with high mean curvature.
  The growth sites are ordered from lowest to highest mean curvature when the
  random growth parameter is set to zero.  Higer values of random growth
  parameter cause more severe shuffling of the ordered list of growth sites.

  @note Used only in this class's copy constructor.

  @return the list of all integer affinities for all the templates for growth of
  all phases
  */
  vector<vector<double>> getAffinity(void) const { return affinity_; }

  vector<vector<double>> getContactAngle(void) const { return contactAngle_; }

  /**
  @brief Set the potassium impurity content for a given clinker phase, on an
  oxide basis.

  @note NOT USED.

  @param idx is the id of the microstructure phase
  @param ival is the mass percentage of potassium to assign to this phase on an
  oxide basis, with units of g per 100 g of the clinker phase
  */
  void setK2o(const int idx, double ival) {
    // try {
    k2o_[idx] = ival;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setK2o", "k2o_", k2o_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Set the sodium impurity content for a given phase, on an oxide basis.

  @note NOT USED.

  @param idx is the id of the microstructure phase
  @param ival is the mass percentage of sodium to assign to this phase on an
  oxide basis, with units of g per 100 g of the clinker phase
  */
  void setNa2o(const int idx, double ival) {
    // try {
    na2o_[idx] = ival;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setNa2o", "na2o_", na2o_.size(),
    //  idx); ex.printException(); exit(1);
    //}
    // return;
  }

  /**
  @brief Set the magnesium impurity content for a given phase, on an oxide
  basis.

  @note NOT USED.

  @param idx is the id of the microstructure phase
  @param ival is the mass percentage of magnesium to assign to this phase on an
  oxide basis, with units of g per 100 g of the clinker phase
  */
  void setMgo(const int idx, double ival) {
    // try {
    mgo_[idx] = ival;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setMgo", "mgo_", mgo_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Set the sulfur impurity content for a given phase, on an oxide basis.

  @note NOT USED.

  @param idx is the id of the microstructure phase
  @param ival is the mass percentage of sulfur to assign to this phase on an
  oxide basis, with units of g per 100 g of the clinker phase
  */
  void setSo3(const int idx, double ival) {
    // try {
    so3_[idx] = ival;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setSo3", "so3_", so3_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Get the potassium impurity content for a given clinker phase, on an
  oxide basis.

  @param idx is the id of the microstructure phase
  @return the mass percentage of potassium to assign to this phase on an oxide
  basis, with units of g per 100 g of the clinker phase
  */
  double getK2o(const int idx) {
    // try {
    return k2o_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getK2o", "k2o_", k2o_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the sodium impurity content for a given clinker phase, on an oxide
  basis.

  @param idx is the id of the microstructure phase
  @return the mass percentage of sodium to assign to this phase on an oxide
  basis, with units of g per 100 g of the clinker phase
  */
  double getNa2o(const int idx) {
    // try {
    return na2o_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getNa2o", "na2o_", na2o_.size(),
    //  idx);
    // ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the magnesium impurity content for a given clinker phase, on an
  oxide basis.

  @param idx is the id of the microstructure phase
  @return the mass percentage of magnesium to assign to this phase on an oxide
  basis, with units of g per 100 g of the clinker phase
  */
  double getMgo(const int idx) {
    // try {
    return mgo_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getMgo", "mgo_", mgo_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the sulfur impurity content for a given clinker phase, on an oxide
  basis.

  @param idx is the id of the microstructure phase
  @param ival is the mass percentage of sulfur to assign to this phase on an
  oxide basis, with units of g per 100 g of the clinker phase
  */
  double getSo3(const int idx) {
    // try {
    return so3_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getSo3", "so3_", so3_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the potassium impurity content for all clinker phases, on an oxide
  basis.

  @note Used only in this class's copy constructor.

  @return a list of the mass percentages of potassium in all clinker phases on
  an oxide basis, with units of g per 100 g of each clinker phase
  */
  vector<double> getK2o(void) const { return k2o_; }

  /**
  @brief Get the sodium impurity content for all clinker phases, on an oxide
  basis.

  @note Used only in this class's copy constructor.

  @return a list of the mass percentages of sodium in all clinker phases on an
  oxide basis, with units of g per 100 g of each clinker phase
  */
  vector<double> getNa2o(void) const { return na2o_; }

  /**
  @brief Get the magnesium impurity content for all clinker phases, on an oxide
  basis.

  @note Used only in this class's copy constructor.

  @return a list of the mass percentages of magnesium in all clinker phases on
  an oxide basis, with units of g per 100 g of each clinker phase
  */
  vector<double> getMgo(void) const { return mgo_; }

  /**
  @brief Get the sulfur impurity content for all clinker phases, on an oxide
  basis.

  @note Used only in this class's copy constructor.

  @return a list of the mass percentages of sulfur in all clinker phases on an
  oxide basis, with units of g per 100 g of each clinker phase
  */
  vector<double> getSo3(void) const { return so3_; }

  /**
  @brief Set the internal porosity of a microstructure phase.

  A few phases, mainly C-S-H gel, have finely dispersed porosity that is not
  resolved at the microstructure scale, so these phases are given a property of
  their average internal porosity on a scale of one micrometer.  This has
  important implications when converting mass fractions of a phase to volume
  fractions within a microstructure, because the GEM thermodynamic phase may or
  may not include such porosity in its value of the molar volume and density of
  such phases.

  @warning Each version of the GEMIPM library
  needs to be checked to see what the convention ensures compatibility with that
  library.

  @param idx is the microstructure phase with internal porosity
  @param pval is the value of the volume fraction of subvoxel pores to assign
  */
  void setMicroPhasePorosity(const int idx, double pval) {
    // try {
    microPhasePorosity_[idx] = pval;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setMicroPhasePorosity",
    //                  "microPhasePorosity_", microPhasePorosity_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Get the internal porosity of a microstructure phase.

  This function is used, for two main purposes: (1) Calculating local mean
  curvature by the template method, and (2) adjusting water content an apparent
  molar volume of a microstructure phase.

  Water or aqueous solution has a porosity of 1. In addition,
  a few phases, mainly C-S-H gel, have finely dispersed porosity that is not
  resolved at the microstructure scale, so these phases are given a property of
  their average internal porosity (subvoxel scale).  This has important
  implications when converting mass fractions of a phase to volume fractions
  within a microstructure, because the GEM thermodynamic phase may or may not
  include such porosity in its value of the molar volume and density of such
  phases.

  @warning Each version of the GEMIPM library
  needs to be checked to see what the convention ensures compatibility with that
  library.

  @param idx is the microstructure phase with internal (subvoxel) porosity
  @return the subvoxel pore volume fraction of the phase
  */
  double getMicroPhasePorosity(const int idx) {
    // try {
    return microPhasePorosity_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getMicroPhasePorosity",
    //                  "microPhasePorosity_", microPhasePorosity_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the internal porosity of a microstructure phase.

  @note NOT USED.

  @warning Each version of the GEMIPM library
  needs to be checked to see what the convention ensures compatibility with that
  library.

  @param str is the name of the microstructure phase with internal (subvoxel)
  porosity
  @return the volume fraction of the phase occupied by pores at the scale of one
  micrometer
  */
  double getMicroPhasePorosity(const string &str) {
    int idx = getMicroPhaseId(str);
    // try {
    return microPhasePorosity_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getMicroPhasePorosity",
    //                  "microPhasePorosity_", microPhasePorosity_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the list of internal porosities of all microstructure phases.

  A few phases, mainly C-S-H gel, have finely dispersed porosity that is not
  resolved at the microstructure scale, so these phases are given a property of
  their average internal porosity on a scale of one micrometer.  This has
  important implications when converting mass fractions of a phase to volume
  fractions within a microstructure, because the GEM thermodynamic phase may or
  may not include such porosity in its value of the molar volume and density of
  such phases.

  @warning Each version of the GEMIPM library
  needs to be checked to see what the convention ensures compatibility with that
  library.

  @note Used only in this class's copy constructor.

  @return the list of porosities of all microstructure phases at the scale of
  one micrometer
  */
  vector<double> getMicroPhasePorosity() const { return microPhasePorosity_; }

  /**
  @brief Get the list of sub-voxel pore size distributions.

  A few phases, mainly C-S-H gel, have finely dispersed porosity that is not
  resolved at the microstructure scale, so these phases are given a property of
  their pore size distribution as a probability density function. These data are
  used to provide better refinement to the moisture distribution in a partially
  saturated microstructure.

  @return the list of pore size distributions of all microstructure phases
  at the scale of one voxel
  */
  vector<vector<struct PoreSizeVolume>> getPoreSizeDistribution() const {
    return poreSizeDistribution_;
  }

  /**
  @brief Set the list of all GEM CSD phases that are associated with a given
  microstructure phase.

  @note NOT USED.

  @param idx is the microstructure phase in question
  @param mpvec is the vector of all GEM CSD phase ids associated with the
  microstructure phase
  */
  void setMicroPhaseMembers(const int idx, vector<int> mpvec) {
    map<int, vector<int>>::iterator p = microPhaseMembers_.find(idx);
    if (p != microPhaseMembers_.end()) {
      p->second = mpvec;
    } else {
      microPhaseMembers_.insert(make_pair(idx, mpvec));
    }
  }

  /**
  @brief Set one of the GEM CSD phases that are associated with a given
  microstructure phase.

  @note NOT USED.

  @param idx is the microstructure phase in question
  @param jdx is the element position in the list of all GEM CSD phases for
  microstructure phase idx
  @param val is GEM CSD phase id to assign to this element in the list
  */
  void setMicroPhaseMembers(const int idx, const int jdx, const int val) {
    string msg;
    map<int, vector<int>>::iterator p = microPhaseMembers_.find(idx);
    if (p != microPhaseMembers_.end()) {
      if (jdx < (int)((p->second).size())) {
        (p->second)[jdx] = val;
      } else {
        EOBException ex("ChemicalSystem", "setMicroPhaseMembers",
                        "microPhaseMembers_", microPhaseMembers_.size(), jdx);
        ex.printException();
        exit(1);
      }
    } else {
      msg = "Could not find microPhaseMembers_ match to index provided";
      EOBException ex("ChemicalSystem", "setMicroPhaseMembers", msg,
                      microPhaseMembers_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the list of all GEM CSD phases that are associated with a given
  microstructure phase id number.

  @note NOT USED.

  @param idx is the microstructure phase in question
  @return the vector of all GEM CSD phase ids associated with the microstructure
  phase
  */
  vector<int> getMicroPhaseMembers(const int idx) {
    string msg;
    map<int, vector<int>>::iterator p = microPhaseMembers_.find(idx);
    if (p != microPhaseMembers_.end()) {
      return p->second;
    } else {
      msg = "Could not find microPhaseMembers_ match to index provided";
      EOBException ex("ChemicalSystem", "getMicroPhaseMembers", msg,
                      microPhaseMembers_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the list of all GEM CSD phases that are associated with a given
  microstructure phase name.

  @note NOT USED.

  @param str is the name of the microstructure phase in question
  @return the vector of all GEM CSD phase ids associated with the microstructure
  phase
  */
  vector<int> getMicroPhaseMembers(const string &str) {
    string msg;
    int idx = getMicroPhaseId(str);
    map<int, vector<int>>::iterator p = microPhaseMembers_.find(idx);
    if (p != microPhaseMembers_.end()) {
      return p->second;
    } else {
      msg = "Could not find microPhaseMembers_ match to " + str;
      EOBException ex("ChemicalSystem", "getMicroPhaseMembers", msg,
                      microPhaseMembers_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get one of the GEM CSD phases that are associated with a given
  microstructure phase id.

  @note NOT USED.

  @param idx is the microstructure phase in question
  @param jdx is the element position in the list of all GEM CSD phases for
  microstructure phase idx
  @return the GEM CSD phase id at element jdx in the list
  */
  int getMicroPhaseMembers(const int idx, const int jdx) {
    string msg;
    map<int, vector<int>>::iterator p = microPhaseMembers_.find(idx);
    if (p != microPhaseMembers_.end()) {
      if (jdx < (int)((p->second).size())) {
        return (p->second)[jdx];
      } else {
        EOBException ex("ChemicalSystem", "getMicroPhaseMembers",
                        "microPhaseMembers_", (p->second).size(), jdx);
        ex.printException();
        exit(1);
      }
    } else {
      msg = "Could not find microPhaseMembers_ match to index provided";
      EOBException ex("ChemicalSystem", "getMicroPhaseMembers", msg,
                      microPhaseMembers_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get one of the GEM CSD phases that are associated with a given
  microstructure phase name.

  @note NOT USED.

  @param str is the name of the microstructure phase in question
  @param jdx is the element position in the list of all GEM CSD phases for
  microstructure phase idx
  @return the GEM CSD phase id at element jdx in the list
  */
  int getMicroPhaseMembers(const string &str, const int jdx) {
    string msg;
    int idx = getMicroPhaseId(str);
    map<int, vector<int>>::iterator p = microPhaseMembers_.find(idx);
    if (p != microPhaseMembers_.end()) {
      if (jdx < (int)((p->second).size())) {
        return (p->second)[jdx];
      } else {
        EOBException ex("ChemicalSystem", "getMicroPhaseMembers",
                        "microPhaseMembers_", (p->second).size(), jdx);
        ex.printException();
        exit(1);
      }
    } else {
      msg = "Could not find microPhaseMembers_ match to index provided";
      EOBException ex("ChemicalSystem", "getMicroPhaseMembers", msg,
                      microPhaseMembers_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the entire list of all GEM CSD phase associations for all
  microstructure phases.

  @note Used only in this class's copy constructor.

  @return the map of all GEM CSD associations for all microstructure phases
  */
  map<int, vector<int>> getMicroPhaseMembers(void) const {
    return microPhaseMembers_;
  }

  map<int, vector<double>> getMicroPhaseMemberVolumeFraction(void) const {
    return microPhaseMemberVolumeFraction_;
  }

  /**
  @brief Get the list of volume fractions of all the GEM CSD phases for a given
  microstructure phase.

  @note NOT USED.

  @param idx is the integer id of the microstructure phase
  @return the vector of volume fractions of each GEM CSD phase for this
  microstructure phase
  */
  vector<double> getMicroPhaseMemberVolumeFraction(const int idx) {
    string msg;
    map<int, vector<double>>::iterator p =
        microPhaseMemberVolumeFraction_.find(idx);
    if (p != microPhaseMemberVolumeFraction_.end()) {
      return p->second;
    } else {
      msg = "Could not find microPhaseMemberVolumeFraction_ match to index "
            "provided";
      EOBException ex("ChemicalSystem", "getMicroPhaseMemberVolumeFraction",
                      msg, microPhaseMemberVolumeFraction_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the volume fraction of one of the GEM CSD phases associated with a
  given microstructure phase.

  @note NOT USED.

  @param idx is the microstructure phase in question
  @param jdx is the element position in the list of all GEM CSD phases for
  microstructure phase idx
  @return the volume fraction of that GEM CSD phase id at element jdx in the
  list
  */
  double getMicroPhaseMemberVolumeFraction(const int idx, const int jdx) {
    string msg;
    map<int, vector<double>>::iterator p =
        microPhaseMemberVolumeFraction_.find(idx);
    if (p != microPhaseMemberVolumeFraction_.end()) {
      if (jdx < (int)((p->second).size())) {
        return (p->second)[jdx];
      } else {
        EOBException ex("ChemicalSystem", "getMicroPhaseMemberVolumeFraction",
                        "microPhaseMemberVolumeFraction_", (p->second).size(),
                        jdx);
        ex.printException();
        exit(1);
      }
    } else {
      msg = "Could not find microPhaseMemberVolumeFraction_ match to index "
            "provided";
      EOBException ex("ChemicalSystem", "getMicroPhaseMemberVolumeFraction",
                      msg, microPhaseMemberVolumeFraction_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Set the list of all dependent component (DC) ids for a microstructure
  phase.

  @note NOT USED.

  @param idx is the microstructure phase in question
  @param mpvec is the list of all DCs for that phase
  */
  void setMicroPhaseDCMembers(const int idx, vector<int> mpvec) {
    map<int, vector<int>>::iterator p = microPhaseDCMembers_.find(idx);
    if (p != microPhaseDCMembers_.end()) {
      p->second = mpvec;
    } else {
      microPhaseDCMembers_.insert(make_pair(idx, mpvec));
    }
  }

  /**
  @brief Set one of the DC component ids for a given microstructure phase id.

  @note NOT USED.

  @param idx is the microstructure phase id in question
  @param jdx is the element position in the list of all GEM CSD phases for
  microstructure phase idx
  @param val is the DC component id to set at that position in the list
  */
  void setMicroPhaseDCMembers(const int idx, const int jdx, const int val) {
    string msg;
    map<int, vector<int>>::iterator p = microPhaseDCMembers_.find(idx);
    if (p != microPhaseDCMembers_.end()) {
      if (jdx < (int)((p->second).size())) {
        (p->second)[jdx] = val;
      } else {
        EOBException ex("ChemicalSystem", "setMicroPhaseDCMembers",
                        "microPhaseDCMembers_", (p->second).size(), jdx);
        ex.printException();
        exit(1);
      }
    } else {
      msg = "Could not find microPhaseDCMembers_ match to index provided";
      EOBException ex("ChemicalSystem", "setMicroPhaseDCMembers", msg,
                      microPhaseDCMembers_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the list of DC component ids associated with a given microstructure
  phase id.

  @note NOT USED.

  @param idx is the microstructure phase id in question
  @return the list of all DC component ids for that phase
  */
  vector<int> getMicroPhaseDCMembers(const int idx) {
    map<int, vector<int>>::iterator p = microPhaseDCMembers_.find(idx);
    if (p != microPhaseDCMembers_.end()) {
      return p->second;
    } else {
      string msg =
          "Could not find microPhaseDCMembers_ match to index provided";
      EOBException ex("ChemicalSystem", "getMicroPhaseDCMembers", msg,
                      microPhaseDCMembers_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the list of DC component ids associated with a given microstructure
  phase name.

  @note NOT USED.

  @param str is the name of the microstructure phase in question
  @return the list of all DC component ids for that phase
  */
  vector<int> getMicroPhaseDCMembers(const string &str) {
    string msg;
    int idx = getMicroPhaseId(str);
    map<int, vector<int>>::iterator p = microPhaseDCMembers_.find(idx);
    if (p != microPhaseDCMembers_.end()) {
      return p->second;
    } else {
      msg = "Could not find microPhaseDCMembers_ match to " + str;
      EOBException ex("ChemicalSystem", "getMicroPhaseDCMembers", msg,
                      microPhaseDCMembers_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get one of the DC component ids for a given microstructure phase id.

  @note NOT USED.

  @param idx is the microstructure phase id in question
  @param jdx is the element position in the list of all DCs for microstructure
  phase idx
  @return the DC component id to set at that position in the list
  */
  int getMicroPhaseDCMembers(const int idx, const int jdx) {

    try {
      // return microPhaseDCMembers_[idx][jdx];
      return microPhaseDCMembers_.at(idx).at(jdx);
    } catch (out_of_range &oor) {
      cout << endl
           << "   ChemicalSystem::getMicroPhaseDCMembers error :" << endl;
      cout << "     Could not find microPhaseDCMembers_ match to indexes "
              "provided"
           << endl;
      cout << "     idx = " << idx << " & jdx = " << jdx << endl;
      if (idx > (int)(microPhaseDCMembers_.size()) - 1) {
        cout << "     microPhaseDCMembers_.size() = "
             << microPhaseDCMembers_.size() << endl;
      } else {
        cout << "     microPhaseDCMembers_.size() = "
             << microPhaseDCMembers_.size() << endl;
        cout << "     microPhaseDCMembers_[idx].size() = "
             << microPhaseDCMembers_[idx].size() << endl;
      }
      cout << endl << "   exit" << endl;
      exit(1);
    }
    /*
        try {
          map<int, vector<int>>::iterator p = microPhaseDCMembers_.find(idx);
    #ifdef DEBUG
          cout << "ChemicalSystem::getMicroPhaseDCMembers micro phase id " <<
    idx
               << " looking for dc index " << jdx << endl;
          cout << "ChemicalSystem::getMicroPhaseDCMembers size = "
               << microPhaseDCMembers_.size() << endl;
          cout.flush();
    #endif
          if (p != microPhaseDCMembers_.end()) {
    #ifdef DEBUG
            cout << "micro phase id " << idx << " looking for dc index " << jdx
                 << endl;
            cout.flush();
    #endif
            if (jdx < (p->second).size()) {
              return (p->second)[jdx];
            } else {
              throw EOBException("ChemicalSystem", "getMicroPhaseDCMembers",
                                 "microPhaseDCMembers_", (p->second).size(),
    jdx);
            }
          } else {
            msg = "Could not find microPhaseDCMembers_ match to index provided";
            throw EOBException("ChemicalSystem", "getMicroPhaseDCMembers", msg,
                               microPhaseDCMembers_.size(), 0);
          }
        } catch (EOBException eex) {
          eex.printException();
          cout.flush();
          cerr.flush();
          exit(1);
        }
    */
  }

  /**
  @brief Get one of the DC component ids for a given microstructure phase name.
  @note NOT USED
  @param str is the microstructure phase name in question
  @param jdx is the element position in the list of all DCs for microstructure
  phase idx
  @return the DC component id to set at that position in the list
  */
  int getMicroPhaseDCMembers(const string &str, const int jdx) {
    string msg;
    int idx = getMicroPhaseId(str);
    map<int, vector<int>>::iterator p = microPhaseDCMembers_.find(idx);
    if (p != microPhaseDCMembers_.end()) {
      if (jdx < (int)((p->second).size())) {
        return (p->second)[jdx];
      } else {
        EOBException ex("ChemicalSystem", "getMicroPhaseDCMembers",
                        "microPhaseDCMembers_", (p->second).size(), jdx);
        ex.printException();
        exit(1);
      }
    } else {
      msg = "Could not find microPhaseDCMembers_ match to index provided";
      EOBException ex("ChemicalSystem", "getMicroPhaseDCMembers", msg,
                      microPhaseDCMembers_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get DC component ids for all microstructure phases at once.

  @note Used only in this class's copy constructor.

  @return the map of all DC component ids for every microstructure phase
  */
  map<int, vector<int>> getMicroPhaseDCMembers(void) const {
    return microPhaseDCMembers_;
  }

  /**
  @brief Get the map of DC component ids associated with a given GEM CSD phase
  id.

  @note Used only in this class's copy constructor

  @return the map of all DC member ids for that GEM phase
  */
  map<int, vector<int>> getGEMPhaseDCMembers(void) const {
    return GEMPhaseDCMembers_;
  }

  /**
  @brief Get the map of DC porosities

  @note Used only in this class's copy constructor

  @return the map of all DC member porosities
  */
  map<int, vector<double>> getMicroPhaseDCPorosities(void) const {
    return microPhaseDCPorosities_;
  }

  /**
  @brief Get the map of DC porosities associated with a given microstructure
  phase id.

  @param idx is the index of the microstructure phase
  @return the vector of all DC member porosities for that microstructure phase
  */
  vector<double> getMicroPhaseDCPorosities(const int idx) {
    string msg;
    map<int, vector<double>>::iterator p = microPhaseDCPorosities_.find(idx);
    if (p != microPhaseDCPorosities_.end()) {
      return p->second;
    } else {
      msg = "Could not find microPhaseDCPorosities_ match to index provided";
      EOBException ex("ChemicalSystem", "getMicroPhaseDCPorosities", msg,
                      microPhaseDCPorosities_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the list of DC component ids associated with a given GEM CSD phase
  id.

  @param idx is the GEM phase id in question
  @return the list of all DC component ids for that GEM phase
  */
  vector<int> getGEMPhaseDCMembers(const int idx) {
    string msg;
    map<int, vector<int>>::iterator p = GEMPhaseDCMembers_.find(idx);
    if (p != GEMPhaseDCMembers_.end()) {
      return p->second;
    } else {
      msg = "Could not find GEMPhaseDCMembers_ match to index provided";
      EOBException ex("ChemicalSystem", "getGEMPhaseDCMembers", msg,
                      GEMPhaseDCMembers_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the list of DC component ids associated with a given GEM CSD phase
  name.

  @note NOT USED.

  @param str is the GEM phase name in question
  @return the list of all DC component ids for that GEM phase
  */
  vector<int> getGEMPhaseDCMembers(const string &str) {
    string msg;
    int pidx = getGEMPhaseId(str);
    map<int, vector<int>>::iterator p = GEMPhaseDCMembers_.find(pidx);
    if (p != GEMPhaseDCMembers_.end()) {
      return p->second;
    } else {
      msg = "Could not find GEMPhaseDCMembers_ match to index provided";
      EOBException ex("ChemicalSystem", "getGEMPhaseDCMembers", msg,
                      GEMPhaseDCMembers_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get one of the DC component ids for a given GEM CSD phase name.

  @note NOT USED.

  @param str is the GEM phase name in question
  @param jdx is the element position in the list of all DCs for GEM phase idx
  @return the DC component id to set at that position in the list
  */
  unsigned int getGEMPhaseDCMembers(const string &str, const unsigned int jdx);

  /**
  @brief Set the number of moles of a given independent component (IC) in the
  system.

  @param idx is the GEM IC id to set
  @param val is the number of moles to assign to that IC
  */
  void setICMoles(const int idx, const double val) {
    if (idx < numICs_) {
      ICMoles_[idx] = val;
    } else {
      EOBException ex("ChemicalSystem", "setICMoles", "ICMoles_", numICs_, idx);
      ex.printException();
      exit(1);
    }
    return;
  }

  /**
  @brief Get the number of moles of every independent component (IC) in the
  system.

  @note Used only in this class's copy constructor.

  @return a pointer to the list of number of moles assigned to each IC
  */
  double *getICMoles(void) const { return ICMoles_; }

  /**
  @brief Get the number of moles of an independent component (IC) in the system.

  @param idx is the GEM IC id to query
  @return the number of moles assigned to that IC
  */
  double getICMoles(const int idx) {
    if (idx < numICs_) {
      return ICMoles_[idx];
    } else {
      EOBException ex("ChemicalSystem", "getICMoles", "ICMoles_", numICs_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Make sure that the IC moles are all greater than e-10
  to ensure stability of the IPM calculation.

  The exception is for charge, which should always be set to zero
  for these kinds of simulations.

  Only the ICMoles_ that are less than 10^-9 after the first call
  of calculateKineticStep(...) are set to 10^-9
  */
  void checkICMoles(void) {
    int i, j;
    vector<double> ICMoles;
    ICMoles.resize(numICs_, 0.0);

    for (j = 0; j < numDCs_; j++) {
      for (i = 0; i < numICs_; i++) {
        ICMoles[i] += (DCMoles_[j] - DCLowerLimit_[j]) * getDCStoich(j, i);
      }
    }

    for (i = 0; i < numICs_; i++) {
      if (ICMoles[i] < ICTHRESH)
        ICMoles_[i] = ICTHRESH;
      if (i == numICs_ - 1) // This IC is always charge
        ICMoles_[i] = 0.0;
    }
    return;
  }

  /**
  @brief Output the number of moles of every independent component (IC) in the
  system.

  */
  void writeICMoles(void) {
    cout << endl;
    cout << "Vector of Independent Components:" << endl;
    for (int i = 0; i < numICs_; i++) {
      cout << "    ICId: " << ICIdLookup_[ICName_[i]] << "\t" << ICName_[i]
           << ": " << ICMoles_[i] << " mol" << endl;
    }
    cout << endl;
    cout.flush();
    return;
  }

  /**
  @brief Set the number of moles of a given dependent component (DC) in the
  system.

  @param idx is the GEM DC id to set
  @param val is the number of moles to assign to that DC
  */
  void setDCMoles(const int idx, const double val) {
    if (idx < numDCs_) {
      DCMoles_[idx] = val;
    } else {
      EOBException ex("ChemicalSystem", "setDCMoles", "DCMoles_", numDCs_, idx);
      ex.printException();
      exit(1);
    }
    return;
  }

  /**
  @brief Get the number of moles of every dependent component (DC) in the
  system.

  @note Used only in this class's copy constructor.

  @return a pointer to the list of number of moles assigned to each DC
  */
  double *getDCMoles(void) const { return DCMoles_; }

  /**
  @brief Get the number of moles of an dependent component (DC) in the system.

  @param idx is the GEM DC id to query
  @return the number of moles assigned to that DC
  */
  double getDCMoles(const int idx) {
    if (idx < numDCs_) {
      return DCMoles_[idx];
    } else {
      EOBException ex("ChemicalSystem", "getDCMoles", "DCMoles_", numDCs_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the number of moles of an dependent component (DC) in the system.

  @param str is the GEM DC name to query
  @return the number of moles assigned to that DC
  */
  double getDCMoles(const string &str) {
    int idx = getDCId(str);
    if (idx < numDCs_) {
      return DCMoles_[idx];
    } else {
      EOBException ex("ChemicalSystem", "getDCMoles", "DCMoles_", numDCs_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Output the number of moles of every independent component (IC) in the
  system.

  */
  void writeDCMoles() {
    cout << endl;
    cout << "Vector of Dependent Components:" << endl;
    for (int i = 0; i < numDCs_; i++) {
      cout << "    DCId: " << DCIdLookup_[DCName_[i]] << "\t" << DCName_[i]
           << ": " << DCMoles_[i] << " mol"
           << " \tmolarVolume = " << getDCMolarVolume(i)
           << " \tvolume = " << DCMoles_[i] * getDCMolarVolume(i) << endl;
    }
    cout << endl;
    cout.flush();
  }

  // void writeDCUpperLimit() {
  //   cout << endl;
  //   cout << "Vector of DCUpperLimit:" << endl;
  //   for (int i = 0; i < numDCs_; i++) {
  //     cout << "    DCId: " << DCIdLookup_[DCName_[i]] << "\t" << DCName_[i]
  //          << ": " << DCUpperLimit_[i] << " mol" << endl;
  //   }
  //   cout << endl;
  //   cout.flush();
  // }

  // void writeDCLowerLimit() {
  //   cout << endl;
  //   cout << "Vector of DCLowerLimit:" << endl;
  //   for (int i = 0; i < numDCs_; i++) {
  //     cout << "    DCId: " << DCIdLookup_[DCName_[i]] << "\t" << DCName_[i]
  //          << ": " << DCLowerLimit_[i] << " mol" << endl;
  //   }
  //   cout << endl;
  //   cout.flush();
  // }

  /**
  @brief Set the molar enthalpy of a given dependent component (DC) in the
  system.

  @param idx is the GEM DC id to set
  @param val is the molar enthalpy to assign to that DC (J/mol)
  */
  void setDCH0(const int idx, const double val) {
    if (idx < numDCs_) {
      DCH0_[idx] = val;
    } else {
      EOBException ex("ChemicalSystem", "setDCH0", "DCH0_", numDCs_, idx);
      ex.printException();
      exit(1);
    }
    return;
  }

  /**
  @brief Get the molar enthalpy of every dependent component (DC) in the system.

  @note Used only in this class's copy constructor.

  @return a pointer to the list of molar enthalpies assigned to each DC (J/mol)
  */
  double *getDCH0(void) const { return DCH0_; }

  /**
  @brief Get the molar enthalpy of an dependent component (DC) in the system.

  @param idx is the GEM DC id to query
  @return the molar enthalpy assigned to that DC (J/mol)
  */
  double getDCH0(const int idx) {
    if (idx < numDCs_) {
      return DCH0_[idx];
    } else {
      EOBException ex("ChemicalSystem", "getDCH0", "DCH0_", numDCs_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the number of moles of an dependent component (DC) in the system.

  @param str is the GEM DC name to query
  @return the molar enthalpy assigned to that DC (J/mol)
  */
  double getDCH0(const string &str) {
    int idx = getDCId(str);
    if (idx < numDCs_) {
      return DCH0_[idx];
    } else {
      EOBException ex("ChemicalSystem", "getDCH0", "DCH0_", numDCs_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Output the number of moles of every independent component (IC) in the
  system.

  */
  void writeDCH0() {
    cout << endl;
    cout << "Vector of Dependent Components:" << endl;
    for (int i = 0; i < numDCs_; i++) {
      cout << "    " << DCName_[i] << ": " << DCH0_[i] << " J/mol" << endl;
    }
    cout << endl;
    cout.flush();
  }

  /**
  @brief Get the total enthalpy of an dependent component (DC) in the system.

  @param idx is the GEM DC id to query
  @return the total enthalpy assigned to that DC (J)
  */
  double getDCEnthalpy(const int idx) {
    if (idx < numDCs_) {
      return (DCH0_[idx] * DCMoles_[idx]);
    } else {
      EOBException ex("ChemicalSystem", "getDCH0", "DCH0_", numDCs_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Set the number of moles of a given GEM CSD phase in the system.

  @note NOT USED.

  @param idx is the GEM CSD phase id to set
  @param val is the number of moles to assign to that GEM phase
  */
  void setGEMPhaseMoles(const int idx, const double val) {
    if (idx < numGEMPhases_) {
      GEMPhaseMoles_[idx] = val;
    } else {
      EOBException ex("ChemicalSystem", "setGEMPhaseMoles", "GEMPhaseMoles_",
                      numGEMPhases_, idx);
      ex.printException();
      exit(1);
    }
    return;
  }

  /**
  @brief Get the number of moles of every GEM CSD phase in the system.

  @note Used only in this class's copy constructor.

  @return a pointer to the list of number of moles assigned to each GEM phase
  */
  double *getGEMPhaseMoles(void) const { return GEMPhaseMoles_; }

  /**
  @brief Get the number of moles of a GEM CSD phase in the system.

  @param idx is the GEM phase id to query
  @return the number of moles assigned to that GEM phase
  */
  double getGEMPhaseMoles(const int idx) {
    if (idx < numGEMPhases_) {
      return GEMPhaseMoles_[idx];
    } else {
      EOBException ex("ChemicalSystem", "getGEMPhaseMoles", "GEMPhaseMoles_",
                      numGEMPhases_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the number of moles of a GEM CSD phase in the system by name

  @param name is the GEM phase name to query
  @return the number of moles assigned to that GEM phase
  */
  double getGEMPhaseMoles(const string &name) {
    int idx = getGEMPhaseId(name);
    if (idx < numGEMPhases_) {
      return GEMPhaseMoles_[idx];
    } else {
      EOBException ex("ChemicalSystem", "getGEMPhaseMoles", "GEMPhaseMoles_",
                      numGEMPhases_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Set the number of moles of all GEM CSD phase in the previous time step.

  @note NOT USED.

  This is just a simple copying of the `GEMPhaseMoles_` vector to
  `prevGEMPhaseMoles_`.
  */
  // void setPrevGEMPhaseMoles(void) {
  //   for (int idx = 0; idx < numGEMPhases_; idx++)
  //     prevGEMPhaseMoles_[idx] = GEMPhaseMoles_[idx];
  // }

  /**
  @brief Set the number of moles of a given GEM CSD phase in the previous time
  step.

  @note NOT USED.

  @param idx is the GEM CSD phase id
  @param val is the number of moles to assign to that GEM phase in the previous
  time step
  */
  // void setPrevGEMPhaseMoles(const unsigned int idx, const double val) {
  //   if (idx < numGEMPhases_) {
  //     prevGEMPhaseMoles_[idx] = val;
  //   } else {
  //     EOBException ex("ChemicalSystem", "setPrevGEMPhaseMoles",
  //                     "prevGEMPhaseMoles_", numGEMPhases_, idx);
  //     ex.printException();
  //     exit(1);
  //   }
  //   return;
  // }

  /**
  @brief Get the number of moles of every GEM CSD phase in the previous time
  step.

  @note Used only in this class's copy constructor.

  @return a pointer to the list of number of moles assigned to each GEM phase in
  last time step
  */
  double *getPrevGEMPhaseMoles(void) const { return prevGEMPhaseMoles_; }

  /**
  @brief Get the number of moles of a GEM CSD phase in the previous time step.

  @note NOT USED.

  @param idx is the GEM phase id to query
  @return the number of moles assigned to that GEM phase in the previous time
  step
  */
  double getPrevGEMPhaseMoles(const int idx) {
    if (idx < numGEMPhases_) {
      return prevGEMPhaseMoles_[idx];
    } else {
      EOBException ex("ChemicalSystem", "getPrevGEMPhaseMoles",
                      "prevGEMPhaseMoles_", numGEMPhases_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Output the number of moles of every GEM phase in the system.

  */
  void writePhaseMoles(void) {
    cout << endl;
    cout << "Vector of Phases:" << endl;
    for (int i = 0; i < numGEMPhases_; i++) {
      cout << "    " << GEMPhaseName_[i] << ": " << GEMPhaseMoles_[i] << " mol"
           << endl;
    }
    cout << endl;
    cout.flush();
  }

  /**
  @brief Set the mass of a given GEM CSD phase [g].

  @note NOT USED.

  @param idx is the GEM CSD phase id
  @param val is the mass to assign to that GEM phase [g]
  */
  void setGEMPhaseMass(const int idx, const double val) {
    if (idx < numGEMPhases_) {
      GEMPhaseMass_[idx] = val;
    } else {
      EOBException ex("ChemicalSystem", "setGEMPhaseMass", "GEMPhaseMass_",
                      numGEMPhases_, idx);
      ex.printException();
      exit(1);
    }
    return;
  }

  /**
  @brief Set the mass of all GEM CSD phases [g].

  */
  void setGEMPhaseMass(void) {
    setPrevGEMPhaseMass();
    for (int i = 0; i < numGEMPhases_; i++) {
      GEMPhaseMass_[i] = (double)(node_->Ph_Mass(i) * 1000.0); // in g, not kg
    }
  }

  /**
  @brief Get the mass of every GEM CSD phase in the system [g].

  @note Used only in this class's copy constructor.

  @return a pointer to the list of masses assigned to each GEM phase [g]
  */
  double *getGEMPhaseMass(void) const { return GEMPhaseMass_; }

  /**
  @brief Get the mass of a GEM CSD phase (by id) [g].

  @param idx is the GEM phase id to query
  @return the mass assigned to that GEM phase [g]
  */
  double getGEMPhaseMass(const int idx) {
    if (idx < numGEMPhases_) {
      return GEMPhaseMass_[idx];
    } else {
      EOBException ex("ChemicalSystem", "getGEMPhaseMass", "GEMPhaseMass_",
                      numGEMPhases_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the mass of a GEM CSD phase (by name) [g].

  @param str is the GEM phase name to query
  @return the mass assigned to that GEM phase [g]
  */
  double getGEMPhaseMass(const string &str) {
    string msg;
    int idx = getGEMPhaseId(str);
    if (idx < numGEMPhases_) {
      return GEMPhaseMass_[idx];
    } else {
      msg = "Name " + str + " does not have a valid phase id";
      EOBException ex("ChemicalSystem", "getPhasemass", msg, numGEMPhases_,
                      idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Set the mass of a given GEM CSD phase in the previous time step [g].

  @note NOT USED.

  @param idx is the GEM CSD phase id
  @param val is the mass to assign to that GEM phase in the previous time step
  [g]
  */
  void setPrevGEMPhaseMass(const int idx, const double val) {
    if (idx < numGEMPhases_) {
      prevGEMPhaseMass_[idx] = val;
    } else {
      EOBException ex("ChemicalSystem", "setPrevGEMPhaseMass",
                      "prevGEMPhaseMass_", numGEMPhases_, idx);
      ex.printException();
      exit(1);
    }
    return;
  }

  /**
  @brief Set the mass of all GEM CSD phases in the previous time step [g].

  */
  void setPrevGEMPhaseMass(void) {
    for (int i = 0; i < numGEMPhases_; i++) {
      prevGEMPhaseMass_[i] = GEMPhaseMass_[i];
    }
  }

  /**
  @brief Get the mass of every GEM CSD phase in the system in the previous time
  step [g].

  @note Used only in this class's copy constructor.

  @return a pointer to the list of masses assigned to each GEM phase in the
  previous time step [g]
  */
  double *getPrevGEMPhaseMass(void) const { return prevGEMPhaseMass_; }

  /**
  @brief Get the mass of a GEM CSD phase (by id) in the previous time step [g].

  @param idx is the GEM phase id to query
  @return the mass assigned to that GEM phase in the previous time step [g]
  */
  double getPrevGEMPhaseMass(const int idx) {
    if (idx < numGEMPhases_) {
      return prevGEMPhaseMass_[idx];
    } else {
      EOBException ex("ChemicalSystem", "getPrevGEMPhaseMass",
                      "prevGEMPhaseMass_", numGEMPhases_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the mass of a GEM CSD phase (by name) in the previous time step.

  @note NOT USED.

  @param str is the GEM phase name to query
  @return the mass assigned to that GEM phase in the previous time step
  */
  // double getPrevGEMPhaseMass(const string &str) {
  //   string msg;
  //   unsigned int idx = getGEMPhaseId(str);
  //   if (idx < numGEMPhases_) {
  //     return prevGEMPhaseMass_[idx];
  //   } else {
  //     msg = "Name " + str + " does not have a valid phase id";
  //     EOBException ex("ChemicalSystem", "getPrevGEMPhaseMass", msg,
  //                     numGEMPhases_, 0);
  //     ex.printException();
  //     exit(1);
  //   }
  // }

  /**
  @brief Set the volume of a given GEM CSD phase.

  @note NOT USED.

  @param idx is the GEM CSD phase id
  @param val is the volume to assign to that GEM phase
  */
  void setGEMPhaseVolume(const int idx, const double val) {
    if (idx < numGEMPhases_) {
      GEMPhaseVolume_[idx] = val;
    } else {
      EOBException ex("ChemicalSystem", "setGEMPhaseVolume", "GEMPhaseVolume_",
                      numGEMPhases_, idx);
      ex.printException();
      exit(1);
    }
    return;
  }

  /**
  @brief Set the volume of all GEM CSD phases.

  */
  void setGEMPhaseVolume(void) {
    // setPrevGEMPhaseVolume();
    for (int i = 0; i < numGEMPhases_; i++) {
      GEMPhaseVolume_[i] = (double)(node_->Ph_Volume(i));
      // if (i == 0)
      //   cout << "   i = " << i << "\t\t" << GEMPhaseName_[i] << "\t" <<
      //   GEMPhaseVolume_[i]
      //        << "\t  mole = " << node_->Ph_Mole(i) << "\t  moles = " <<
      //        node_->Ph_Moles(i)
      //        << "\t   volMolar = " << GEMPhaseVolume_[i]/node_->Ph_Mole(i) <<
      //        endl;
    }
  }

  /**
  @brief Get the volume of every GEM CSD phase in the system.

  @note Used only in this class's copy constructor.

  @return a pointer to the list of volumes assigned to each GEM phase
  */
  double *getGEMPhaseVolume() const { return GEMPhaseVolume_; }

  /**
  @brief Get the volume of a GEM CSD phase (by id).

  @note NOT USED.

  @param idx is the GEM phase id to query
  @return the volume assigned to that GEM phase
  */
  double getGEMPhaseVolume(const int idx) {
    if (idx < numGEMPhases_) {
      return GEMPhaseVolume_[idx];
    } else {
      EOBException ex("ChemicalSystem", "getGEMPhaseVolume", "GEMPhaseVolume_",
                      numGEMPhases_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the volume of a GEM CSD phase (by name).

  @param str is the GEM phase name to query
  @return the volume assigned to that GEM phase
  */
  double getGEMPhaseVolume(const string &str) {
    string msg;
    int idx = getGEMPhaseId(str);
    if (idx < numGEMPhases_) {
      return GEMPhaseVolume_[idx];
    } else {
      msg = "Name " + str + " does not have a valid phase id";
      EOBException ex("ChemicalSystem", "getPhasevolume", msg, numGEMPhases_,
                      idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Set the volume of a given GEM CSD phase in the previous time step.

  @note NOT USED.

  @param idx is the GEM CSD phase id
  @param val is the volume to assign to that GEM phase in the previous time step
  */
  void setPrevGEMPhaseVolume(const int idx, const double val) {
    if (idx < numGEMPhases_) {
      prevGEMPhaseVolume_[idx] = val;
    } else {
      EOBException ex("ChemicalSystem", "setPrevGEMPhaseVolume",
                      "prevGEMPhaseVolume_", numGEMPhases_, idx);
      ex.printException();
      exit(1);
    }
    return;
  }

  /**
  @brief Get the volume of every GEM CSD phase in the system in the previous
  time step.

  @return a pointer to the list of volumes assigned to each GEM phase in the
  previous time step
  */
  void setPrevGEMPhaseVolume(void) {
    for (int i = 0; i < numGEMPhases_; i++) {
      prevGEMPhaseVolume_[i] = GEMPhaseVolume_[i];
    }
  }

  /**
  @brief Get the subvoxel porosity of every GEM CSD phase in the system in the
  previous time step.

  @note Does nothing for now, just a placeholder in case it is needed later
  */
  // void setPrevGEMPhasePorosity(void) {
  // for (long int i = 0; i < numGEMPhases_; i++) {
  //     prevGEMPhaseVolume_[i] = GEMPhaseVolume_[i];
  // }
  //  return;
  // }

  /**
  @brief Get the volume of every GEM CSD phase in the system in the previous
  time step.

  @note Used only in this class's copy constructor.

  @return a pointer to the list of volumes assigned to each GEM phase in the
  previous time step
  */
  double *getPrevGEMPhaseVolume(void) const { return prevGEMPhaseVolume_; }

  /**
  @brief Get the volume of a GEM CSD phase (by id) in the previous time step.

  @note NOT USED.

  @param idx is the GEM phase id to query
  @return the volume assigned to that GEM phase in the previous time step
  */
  double getOphasevolume(const int idx) {
    if (idx < numGEMPhases_) {
      return prevGEMPhaseVolume_[idx];
    } else {
      EOBException ex("ChemicalSystem", "getPrevGEMPhaseVolume",
                      "prevGEMPhaseVolume_", numGEMPhases_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the volume of a GEM CSD phase (by name) in the previous time step.

  @note NOT USED.

  @param str is the GEM phase name to query
  @return the volume assigned to that GEM phase in the previous time step
  */
  double getPrevGEMPhaseVolume(const string &str) {
    string msg;
    int idx = getGEMPhaseId(str);
    if (idx < numGEMPhases_) {
      return prevGEMPhaseVolume_[idx];
    } else {
      msg = "Name " + str + " does not have a valid phase id";
      EOBException ex("ChemicalSystem", "getPrevGEMPhaseVolume", msg,
                      numGEMPhases_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Set the moles of a GEM CSD phase identified as a solvent.

  @note NOT USED.

  @param idx is the GEM CSD solvent id
  @param val is the moles to assign to that GEM phase
  */
  void setCarrier(const int idx, const double val) {
    if (idx < numSolutionPhases_) {
      carrier_[idx] = val;
    } else {
      EOBException ex("ChemicalSystem", "setCarrier", "carrier_",
                      numSolutionPhases_, idx);
      ex.printException();
      exit(1);
    }
    return;
  }

  /**
  @brief Get the moles of every GEM CSD phase identified as a solvent.

  @note Used only in this class's copy constructor.

  @return a pointer to the list of moles assigned to each GEM solvent phase
  */
  double *getCarrier(void) const { return carrier_; }

  /**
  @brief Get the moles of a GEM CSD solvent phase (by id).

  @note NOT USED.

  @param idx is the GEM solvent phase name to query
  @return the moles assigned to that GEM solvent phase
  */
  double getCarrier(const int idx) {
    if (idx < numSolutionPhases_) {
      return carrier_[idx];
    } else {
      EOBException ex("ChemicalSystem", "getCarrier", "carrier_",
                      numSolutionPhases_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Set the surface area of a GEM CSD phase (m<sup>2</sup>/kg).

  @note NOT USED.

  @param idx is the GEM CSD phase id
  @param val is the surface area to assign to that GEM phase
  */
  void setSpecificSurfaceArea(const int idx, const double val) {
    if (idx < numGEMPhases_) {
      specificSurfaceArea_[idx] = val;
    } else {
      EOBException ex("ChemicalSystem", "setSpecificSurfaceArea",
                      "specificSurfaceArea_", numGEMPhases_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the specific surface area of all GEM CSD phases (m<sup>2</sup>/kg).

  @note Used only in this class's copy constructor.

  @return a pointer to the list of all phase surface areas
  */
  double *getSpecificSurfaceArea(void) const { return specificSurfaceArea_; }

  /**
  @brief Get the specific surface area of a GEM CSD phase (m<sup>2</sup>/kg).

  @note NOT USED.

  @param idx is the GEM CSD phase id
  @return the specific surface area assigned to that GEM phase (m2/kg)
  */
  double getSpecificSurfaceArea(const int idx) {
    if (idx < numGEMPhases_) {
      return specificSurfaceArea_[idx];
    } else {
      EOBException ex("ChemicalSystem", "getSpecificSurfaceArea",
                      "specificSurfaceArea_", numGEMPhases_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Set the volume of a microstructure phase (by id).

  @param idx is the microstructure phase id
  @param val is the volume to assign to that microstructure phase
  */
  void setMicroPhaseVolume(const int idx, const double val) {
    // try {
    microPhaseVolume_[idx] = val;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setMicroPhaseVolume",
    //                  "microPhaseVolume_", microPhaseVolume_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}

    // Calculate the subvoxel porosity of this phase as well
    if (idx != VOIDID) {
      if (verbose_) {
        cout << "Going into calcMicroPhasePorosity(" << idx << ")" << endl;
        cout.flush();
      }
      calcMicroPhasePorosity(idx);
    }

    // return;
  }

  /**
  @brief Calculate the subvoxel porosity of a microstructure phase (by id).

  @param idx is the microstructure phase id
  */
  void calcMicroPhasePorosity(const unsigned int idx);

  /**
  @brief Get the volumes of all microstructure phases.

  @return a vector of volumes of every microstructure phase
  */
  vector<double> getMicroPhaseVolume(void) const { return microPhaseVolume_; }

  /**
  @brief Get the volume of a microstructure phase (by id).

  @param idx is the microstructure phase id
  @return the volume assigned to that microstructure phase
  */
  double getMicroPhaseVolume(const int idx) {
    // try {
    return microPhaseVolume_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getMicroPhaseVolume",
    //                  "microPhaseVolume_", microPhaseVolume_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Set the mass of a microstructure phase (by id).

  @param idx is the microstructure phase id
  @param val is the mass to assign to that microstructure phase
  */
  void setMicroPhaseMass(const unsigned int idx, const double val);

  /**
  @brief Get the masses of all microstructure phases.

  @note Used only in this class's copy constructor.

  @return a pointer to the list of masses of every microstructure phase
  */
  vector<double> getMicroPhaseMass(void) const { return microPhaseMass_; }

  /**
  @brief Get the mass of a microstructure phase (by id).

  @param idx is the microstructure phase id
  @return the mass assigned to that microstructure phase
  */
  double getMicroPhaseMass(const int idx) {
    // try {
    return microPhaseMass_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getMicroPhaseMass",
    //  "microPhaseMass_",
    //                  microPhaseMass_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Output all the microstructure phases.

  @note NOT USED.
  */
  void writeMicroPhases(void) {
    cout << endl;
    cout << "Microstructure phase quantities:" << endl;
    cout << "Name     Mass (g)     Volume (m3)" << endl;
    cout << "----     --------     -----------" << endl;
    int size = microPhaseName_.size();
    for (int i = 1; i < size; i++) {
      cout << microPhaseName_[i] << "     " << microPhaseMass_[i] << "     "
           << microPhaseVolume_[i] << endl;
    }
    cout << "Void     0.0    " << microVoidVolume_ << endl;
    cout << endl;
    cout.flush();
  }

  /**
  @brief Set the mass dissolved of a microstructure phase (by id).

  @note NOT USED.

  @param idx is the microstructure phase id
  @param val is the dissolved mass to assign to that microstructure phase
  */
  void setMicroPhaseMassDissolved(const int idx, const double val) {
    // try {
    microPhaseMassDissolved_[idx] = val;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setMicroPhaseMassDissolved",
    //                  "microPhaseMassDissolved_",
    //                  microPhaseMassDissolved_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Get the dissolved masses of all microstructure phases.

  @note Used only in this class's copy constructor.

  @return a pointer to the list of dissolved masses of every microstructure
  phase
  */
  vector<double> getMicroPhaseMassDissolved(void) const {
    return microPhaseMassDissolved_;
  }

  /**
  @brief Get the dissolved mass of a microstructure phase (by id).

  @note NOT USED.

  @param idx is the microstructure phase id
  @return the dissolved mass assigned to that microstructure phase
  */
  double getMicroPhaseMassDissolved(const int idx) {
    // try {
    return microPhaseMassDissolved_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getMicroPhaseMassDissolved",
    //                  "microPhaseMassDissolved_",
    //                  microPhaseMassDissolved_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Set the volume of void space in the microstructure.

  @param val is the volume of void space to assign to the microstructure
  */
  void setMicroVoidVolume(const double val) { microVoidVolume_ = val; }

  /**
  @brief Get the volume of void space in the microstructure.

  @note Used only in this class's copy constructor.

  @return the volume of void space in the microstructure
  */
  double getMicroVoidVolume(void) const { return microVoidVolume_; }

  /**
  @brief Set the total microstructure volume.

  @note NOT USED.

  @param val is the total volume to assign to the microstructure
  */
  void setMicroVolume(const double val) { microVolume_ = val; }

  /**
  @brief Get the total microstructure volume.

  @return the total volume of the microstructure
  */
  double getMicroVolume(void) const { return microVolume_; }

  // double getNewMicroVolume(void) const { return newMicroVolume_; }

  /**
  @brief Set the initial total microstructure volume.

  @param val is the initial total volume to assign to the microstructure
  */
  void setInitMicroVolume(const double val) { initMicroVolume_ = val; }

  /**
  @brief Get the initial total microstructure volume.

  @return the initial total volume of the microstructure
  */
  double getInitMicroVolume(void) const { return initMicroVolume_; }

  /**
  @brief Set the upper bound on the moles of a dependent component.

  @param idx is the id of the DC component
  @param val is the upper bound to set (in scaled moles)
  */
  void setDCUpperLimit(const int idx, const double val) {
    if (idx < numDCs_) {
      DCUpperLimit_[idx] = val;
    } else {
      EOBException ex("ChemicalSystem", "setDCUpperLimit", "DCUpperLimit_",
                      numDCs_, idx);
      ex.printException();
      exit(1);
    }
    return;
  }

  /**
  @brief Get the upper bound on the moles of every dependent component (DC).

  @note Used only in this class's copy constructor.

  @return a pointer to the list of DC upper bounds
  */
  double *getDCUpperLimit(void) const { return DCUpperLimit_; }

  /**
  @brief Get the upper bound on the moles of a dependent component (DC).

  @note NOT USED.

  @param idx is the id of the DC component
  @return the upper bound on moles of that DC
  */
  double getDCUpperLimit(const int idx) {
    if (idx < numDCs_) {
      return DCUpperLimit_[idx];
    } else {
      EOBException ex("ChemicalSystem", "getDCUpperLimit", "DCUpperLimit_",
                      numDCs_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the lower bound on the moles of a dependent component (DC).

  @param idx is the id of the DC component
  @param val is the lower bound to set for this DC (in scaled moles)
  */
  void setDCLowerLimit(const int idx, const double val) {
    if (idx < numDCs_) {
      DCLowerLimit_[idx] = val;
    } else {
      EOBException ex("ChemicalSystem", "setDCLowerLimit", "DCLowerLimit_",
                      numDCs_, idx);
      ex.printException();
      exit(1);
    }
    return;
  }

  /**
  @brief Initialize the lower bound on a DC to a prescribed value

  @param val is the lower bound to set for this DC (in scaled moles)
  */
  void initDCLowerLimit(double val) {
    for (int i = 0; i < numDCs_; i++) {
      DCLowerLimit_[i] = val;
    }
  }

  /**
  @brief Get the lower bound on the moles of a dependent component (DC).

  @note Used only in this class's copy constructor.

  @param idx is the id of the DC component
  @return a pointer to the list of DC lower bounds
  */
  double *getDCLowerLimit(void) const { return DCLowerLimit_; }

  /**
  @brief Get the lower bound on the moles of a dependent component (DC).

  @note NOT USED.

  @param idx is the id of the DC component
  @return a pointer to the list of DC lower bounds
  */
  double getDCLowerLimit(const int idx) {
    if (idx < numDCs_) {
      return DCLowerLimit_[idx];
    } else {
      EOBException ex("ChemicalSystem", "getDCLowerLimit", "DCLowerLimit_",
                      numDCs_, idx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Set the molar mass of a particular independent component (IC) [g/mol].

  @note NOT USED.

  @param idx is the id of the IC
  @param val is the molar mass to assign to that IC [g/mol]
  */
  void setICMolarMass(const int idx, const double val) {
    // try {
    ICMolarMass_[idx] = val;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setICMolarMass", "ICMolarMass_",
    //                  ICMolarMass_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Get the molar mass of a particular independent component (IC), by id
  [g/mol].

  @param idx is the id of the IC
  @return the molar mass of that IC, [g/mol]
  */
  double getICMolarMass(const int idx) {
    // try {
    return ICMolarMass_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getICMolarMass", "ICMolarMass_",
    //                  ICMolarMass_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the molar mass of a particular independent component (IC), by name
  [g/mol].

  @param str is the name of the IC
  @return the molar mass of that IC, [g/mol]
  */
  double getICMolarMass(const string &str) {
    // string msg;
    int idx = getICId(str);
    // try {
    return ICMolarMass_[idx];
    //} catch (out_of_range &oor) {
    //  msg = "Name " + str + " does not have a valid phase id";
    //  EOBException ex("ChemicalSystem", "getICMolarMass", msg,
    //                  ICMolarMass_.size(), 0);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the molar masses of all independent components (IC) [g/mol].

  @note Used only in this class's copy constructor.

  @return the vector of IC molar masses [g/mol]
  */
  vector<double> getICMolarMass(void) const { return ICMolarMass_; }

  /**
  @brief Set the molar mass of a particular dependent component (DC) [g/mol].

  @note NOT USED.

  @param idx is the id of the DC
  @param val is the molar mass to assign to that DC [g/mol]
  */
  void setDCMolarMass(const int idx, const double val) {
    // try {
    DCMolarMass_[idx] = val;
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setDCMolarMass", "DCMolarMass_",
    //                  DCMolarMass_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Get the molar mass of a particular dependent component (DC), by id
  [g/mol].

  @param idx is the id of the DC
  @return the molar mass of that DC, [g/mol]
  */
  double getDCMolarMass(const int idx) {
    // try {
    return DCMolarMass_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getDCMolarMass", "DCMolarMass_",
    //                  DCMolarMass_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the molar mass of a particular dependent component (DC), by name
  [g/mol].

  @param str is the name of the DC
  @return the molar mass of that DC, [g/mol]
  */
  double getDCMolarMass(const string &str) {
    // string msg;
    int idx = getDCId(str);
    // try {
    return DCMolarMass_[idx];
    //} catch (out_of_range &oor) {
    //  msg = "Name " + str + " does not have a valid phase id";
    //  EOBException ex("ChemicalSystem", "getDCMolarMass", msg,
    //                  DCMolarMass_.size(), 0);
    // ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the molar masses of all dependent components (DC) [g/mol].

  @note Used only in this class's copy constructor.

  @return the vector of DC molar masses [g/mol]
  */
  vector<double> getDCMolarMass(void) const { return DCMolarMass_; }

  /**
  @brief Get the molar volume of a particular dependent component (DC), by id
  [m3].

  @param dcidx is the id of the DC
  @return the molar volume of that DC, [m3]
  */
  double getDCMolarVolume(const int dcidx) {
    double v0 = node_->DC_V0(dcidx, P_, T_);
    return v0;
  }

  /**
  @brief Get the molar volume of a particular dependent component (DC), by name
  [m3].

  @param str is the name of the DC
  @return the molar volume of that DC, [m3]
  */
  double getDCMolarVolume(const string &str) {
    string msg;
    double V0 = getDCMolarVolume(getDCId(str));
    return V0;
  }

  /**
  @brief Get the molar masses of all GEM CSD phases [g/mol].

  @note NOT USED.

  @return the vector of GEM phase molar masses [g/mol]
  */
  vector<double> getGEMPhaseMolarMass(void) const { return GEMPhaseMolarMass_; }

  /**
  @brief Set the molar masses of all GEM CSD phases [g/mol].

  */
  void setGEMPhaseMolarMass(void) {
    double pmm;
    GEMPhaseMolarMass_.clear();
    GEMPhaseMolarMass_.resize(numGEMPhases_, 0.0);
    int size_pId = GEMPhaseStoich_.size();
    int size_iId;
    for (int pidx = 0; pidx < size_pId; pidx++) {
      pmm = 0.0;
      size_iId = GEMPhaseStoich_[pidx].size();
      for (int icidx = 0; icidx < size_iId; icidx++) {
        pmm += ((GEMPhaseStoich_[pidx][icidx]) * getICMolarMass(icidx));
      }
      GEMPhaseMolarMass_[pidx] = pmm;
    }
  }

  /**
  @brief Get the molar mass of a particular GEM CSD phase, by id [g/mol].

  @note NOT USED.

  @param idx is the id of the GEM phase
  @return the molar mass of that GEM phase [g/mol]
  */
  double getGEMPhaseMolarMass(const int idx) {
    // try {
    return GEMPhaseMolarMass_[idx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getGEMPhaseMolarMass",
    //                  "GEMPhaseMolarMass_", GEMPhaseMolarMass_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the molar mass of a particular GEM CSD phase, by name [g/mol].

  @note NOT USED.

  @param str is the name of the GEM phase
  @return the molar mass of that GEM phase [g/mol]
  */
  double getGEMPhaseMolarMass(const string &str) {
    // string msg;
    int idx = getGEMPhaseId(str);
    // try {
    return GEMPhaseMolarMass_[idx];
    //} catch (out_of_range &oor) {
    //  msg = "Name " + str + " does not have a valid phase id";
    //  EOBException ex("ChemicalSystem", "getGEMPhaseMolarMass", msg,
    //                  GEMPhaseMolarMass_.size(), idx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the number of moles of a particular IC in a particular DC.

  Dependent components (DC) comprise one or (usually) more independent
  components (IC). For example, the DC component CO2(g) has one mole of C and 2
  moles of O per mole of the DC.

  @param dcidx is the id of the DC being queried
  @param icidx is the id of the IC that may be contained within this DC
  @return the number of moles of the ID icidx contained within the DC dcidx
  */
  double getDCStoich(const int dcidx, const int icidx) {
    if (dcidx >= (int)(DCStoich_.size())) {
      EOBException ex("ChemicalSystem", "getDCStoich", "DCStoich_",
                      DCStoich_.size(), dcidx);
      ex.printException();
      exit(1);
    }
    if (icidx >= (int)(DCStoich_[dcidx].size())) {
      EOBException ex("ChemicalSystem", "getDCStoich", "DCStoich_",
                      DCStoich_[dcidx].size(), icidx);
      ex.printException();
      exit(1);
    }
    return DCStoich_[dcidx][icidx];
  }

  /**
  @brief Get the number of moles all ICs in a particular DC.

  Dependent components (DC) comprise one or (usually) more independent
  components (IC). For example, the DC component CO2(g) has one mole of C and 2
  moles of O per mole of the DC.

  @note NOT USED.

  @param dcidx is the id of the DC being queried
  @return the list of moles of each IC in the DC specified by dcidx
  */
  vector<double> getDCStoich(const int dcidx) {
    // try {
    return DCStoich_[dcidx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getDCStoich", "DCStoich_",
    //                  DCStoich_.size(), dcidx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the number of moles all ICs in all DCs.

  Dependent components (DC) comprise one or (usually) more independent
  components (IC). For example, the DC component CO2(g) has one mole of C and 2
  moles of O per mole of the DC.

  @note Used only in this class's copy constructor.

  @return the matrix of moles of each IC in every DCs.
  */
  vector<vector<double>> getDCStoich(void) const { return DCStoich_; }

  /**
  @brief Get the charge of a particular DC.

  Dependent components (DC) comprise one or (usually) more independent
  components (IC). The last component is the electric charge (valence)

  @param dcidx is the id of the DC being queried
  @return the charge of the DC specified by dcidx (valence)
  */
  double getDCCharge(const int dcidx) {
    try {
      return DCCharge_.at(dcidx);
    } catch (out_of_range &oor) {
      EOBException ex("ChemicalSystem", "getDCCharge", "DCCharge_",
                      DCStoich_.size(), dcidx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the charge of a particular DC.

  Dependent components (DC) comprise one or (usually) more independent
  components (IC). The last component is the electric charge (valence)

  @param dcname is the name of the DC being queried
  @return the charge of the DC specified by dcidx (valence)
  */
  double getDCCharge(const string &dcname) {
    double chg = getDCCharge(getDCId(dcname));
    return (chg);
  }

  /**
  @brief Get the charges of all DCs.

  Dependent components (DC) comprise one or (usually) more independent
  components (IC). The last component is the electric charge (valence)

  @note Used only in this class's copy constructor.

  @return the charge of all DCs (valence)
  */
  vector<double> getDCCharge(void) const { return DCCharge_; }

  /**
  @brief Get the number of moles all ICs in all GEM CSD phases.

  GEM phases comprise one or (usually) more independent components (IC).
  For example, one mole of the GEM phase gypsum has one mole of CaSO4 and 1 mole
  of calcium, one mole of sulfur, four moles of hydrogen and six moles of oxygen

  This method makes calls to the GEM node to get the matrix of GEM phase
  stoichiometries, and assigns the results into the 1D `pGEMPhaseStoich_` array
  within the method itself.  For this reason, the function does not return
  anything.

  @todo Consider renaming this method because it does not actually return
  anything to the calling function
  */
  void getPGEMPhaseStoich(void);

  /**
  @brief Set the matrix of moles all ICs in all GEM CSD phases.

  GEM phases comprise one or (usually) more independent components (IC).
  For example, one mole of the GEM phase gypsum has one mole of CaSO4 and 1 mole
  of calcium, one mole of sulfur, four moles of hydrogen and six moles of oxygen

  This particular function just calls the getPGEMPhaseStoich method, which sets
  all the phase stoichiometries into the 1D `pGEMPhaseStoich_` array.

  */
  void setPGEMPhaseStoich(void) {
    getPGEMPhaseStoich();
    return;
  }

  /**
  @brief Set the moles of one IC component of one GEM CSD phase.

  GEM phases comprise one or (usually) more independent components (IC).
  For example, one mole of the GEM phase gypsum has one mole of CaSO4 and 1 mole
  of calcium, one mole of sulfur, four moles of hydrogen and six moles of oxygen

  @note NOT USED.

  @param pidx is the GEM phase id to modify
  @param icidx is the IC component index
  @param val is the number of moles of IC icidx in GEM phase pidx
  */
  void setPGEMPhaseStoich(const int pidx, const int icidx, const double val) {
    int idx = (pidx * numICs_) + icidx;
    if (idx < numGEMPhases_) {
      pGEMPhaseStoich_[idx] = val;
    } else {
      EOBException ex("ChemicalSystem", "setPGEMPhaseStoich",
                      "pGEMPhaseStoich_", numGEMPhases_, idx);
      ex.printException();
      exit(1);
    }
    return;
  }

  /**
  @brief Renormalizes phase stoichiometries per mole of oxygen in the phase.

  This is a poorly worded function because it does not return anything.
  It simply renormalizes the GEM phase stoichiometries by the number of moles
  of oxygen in that phase.

  @todo Find out why it is useful to normalize by moles of oxygen here but not
  for the 1D array version
  @todo Rename this function to something more descriptive
  */
  void getGEMPhaseStoich(void);

  /**
  @brief Caller function for `getGEMPhaseStoich`.

  */
  void setGEMPhaseStoich(void) {
    getGEMPhaseStoich();
    return;
  }

  /**
  @brief Get the number of moles of an IC in a GEM phase.

  This function uses the 1D array pointed to by the `pGEMPhaseStoich_` member
  variable.  This variable lists the IC molar stoichiometry of phase one in the
  first n elements (where n is the number of ICs), and of phase two in the
  second n elements, etc., instead of storing as a 2D matrix

  @note NOT USED.

  @param pidx is the id of the GEM phase to query
  @param icidx is the index of the IC being queried
  @return the molar stoichiometry of the icidx-th IC in phase id pidx
  */
  double getPGEMPhaseStoich(const int pidx, const int icidx) {
    int index = (pidx * numICs_) + icidx;
    if (index < numGEMPhases_) {
      return pGEMPhaseStoich_[index];
    } else {
      EOBException ex("ChemicalSystem", "getPGEMPhaseStoich",
                      "pGEMPhaseStoich_", numGEMPhases_, index);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the list of moles of all ICs in a GEM phase.

  This function uses the 1D array pointed to by the `pGEMPhaseStoich_` member
  variable.  This variable lists the IC molar stoichiometry of phase one in the
  first n elements (where n is the number of ICs), and of phase two in the
  second n elements, etc., instead of storing as a 2D matrix

  @param pidx is the id of the GEM phase to query
  @return a pointer to the list of IC moles in phase id pidx
  */
  double *getPGEMPhaseStoich(const int pidx) {
    if (pidx >= numGEMPhases_) {
      EOBException ex("ChemicalSystem", "getPGEMPhaseStoich",
                      "pGEMPhaseStoich_", numGEMPhases_, pidx);
      ex.printException();
      exit(1);
    }
    double *p = pGEMPhaseStoich_;
    p += (pidx * numICs_);
    return p;
  }

  /**
  @brief Get the number of moles of an IC in a GEM phase by phase id.

  This function uses the 2D matrix form of the molar stoichiometries,
  implemented as a vector of vectors, `GEMPhaseStoich_`.  Each row is a phase
  and each column is an IC.

  @note NOT USED.

  @param pidx is the id of the GEM phase to query
  @param icidx is the index of the IC being queried
  @return the molar stoichiometry of the icidx-th IC in phase id pidx
  */
  double getGEMPhaseStoich(const int pidx, const int icidx) {
    if (pidx >= (int)(GEMPhaseStoich_.size())) {
      EOBException ex("ChemicalSystem", "getGEMPhaseStoich", "GEMPhaseStoich_",
                      GEMPhaseStoich_.size(), pidx);
      ex.printException();
      exit(1);
    }
    if (icidx >= (int)(GEMPhaseStoich_[pidx].size())) {
      EOBException ex("ChemicalSystem", "getGEMPhaseStoich", "GEMPhaseStoich_",
                      GEMPhaseStoich_[pidx].size(), icidx);
      ex.printException();
      exit(1);
    }
    return ((double)(GEMPhaseStoich_[pidx][icidx]));
  }

  /**
  @brief Get the number of moles of an IC in a GEM phase by phase name.

  This function uses the 2D matrix form of the molar stoichiometries,
  implemented as a vector of vectors, `GEMPhaseStoich_`.  Each row is a phase
  and each column is an IC.

  @param str is the name of the GEM phase to query
  @param icidx is the index of the IC being queried
  @return the molar stoichiometry of the icidx-th IC in phase id pidx
  */
  double getGEMPhaseStoich(const string &str, const int icidx) {
    int pidx = getGEMPhaseId(str);
    if (pidx >= (int)(GEMPhaseStoich_.size())) {
      EOBException ex("ChemicalSystem", "getGEMPhaseStoich", "GEMPhaseStoich_",
                      GEMPhaseStoich_.size(), pidx);
      ex.printException();
      exit(1);
    }
    if (icidx >= (int)(GEMPhaseStoich_[pidx].size())) {
      EOBException ex("ChemicalSystem", "getGEMPhaseStoich", "GEMPhaseStoich_",
                      GEMPhaseStoich_[pidx].size(), icidx);
      ex.printException();
      exit(1);
    }
    return ((double)(GEMPhaseStoich_[pidx][icidx]));
  }

  /**
  @brief Get the number of moles of each IC in a GEM phase by phase id.

  This function uses the 2D matrix form of the molar stoichiometries,
  implemented as a vector of vectors, `GEMPhaseStoich_`.  Each row is a phase
  and each column is an IC.

  @note NOT USED.

  @param pidx is the id of the GEM phase to query
  @return the list of moles of each IC in phase id pidx
  */
  vector<double> getGEMPhaseStoich(const int pidx) {
    // try {
    return GEMPhaseStoich_[pidx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getGEMPhaseStoich",
    //  "GEMPhaseStoich_",
    //                  GEMPhaseStoich_.size(), pidx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the number of moles of each IC in a GEM phase by phase name.

  This function uses the 2D matrix form of the molar stoichiometries,
  implemented as a vector of vectors, `GEMPhaseStoich_`.  Each row is a phase
  and each column is an IC.

  @note NOT USED.

  @param str is the name of the GEM phase to query
  @return the list of moles of each IC in phase id pidx
  */
  vector<double> getGEMPhaseStoich(const string &str) {
    int pidx = getGEMPhaseId(str);
    // try {
    return GEMPhaseStoich_[pidx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getGEMPhaseStoich",
    //  "GEMPhaseStoich_",
    //                  GEMPhaseStoich_.size(), pidx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the 1D array of IC stoichiometries of each GEM CSD phase.

  This function uses the 1D array pointed to by the `pGEMPhaseStoich_` member
  variable.  This variable lists the IC molar stoichiometry of phase one in the
  first n elements (where n is the number of ICs), and of phase two in the
  second n elements, etc., instead of storing as a 2D matrix

  @return a pointer to the list of moles of each IC in every GEM phase
  */
  double *getPGEMPhaseStoich(void) const { return pGEMPhaseStoich_; }

  /**
  @brief Get the 2D matrix of IC stoichiometries of each GEM phase.

  This function uses the 2D matrix form of the molar stoichiometries,
  implemented as a vector of vectors, `GEMPhaseStoich_`.  Each row is a phase
  and each column is an IC.

  @note Used only in this class's copy constructor.

  @return the 2D matrix of moles of each IC in each GEM phase
  */
  vector<vector<double>> getGEMPhaseStoich(void) const {
    return GEMPhaseStoich_;
  }

  /**
  @brief Get the chemical potential of an independent component (IC) [mol/mol].

  @note NOT USED.

  @param icidx is the index of the IC being queried
  @return the chemical potential of the IC
  */
  double getICChemicalPotential(const int icidx) {
    if (icidx < numICs_) {
      return ICChemicalPotential_[icidx];
    } else {
      EOBException ex("ChemicalSystem", "getICChemicalPotential",
                      "ICChemicalPotential_", numICs_, icidx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the list of chemical potentials of each independent component (IC)
  [mol/mol].

  @note Used only in this class's copy constructor.

  @return a pointer to the list of IC chemical potentials [mol/mol]
  */
  double *getICChemicalPotential(void) const { return ICChemicalPotential_; }

  /**
  @brief Get the molal activity coefficient of a dependent component (DC).

  @note NOT USED.

  @param dcidx is the index of the DC being queried
  @return the activity coefficient of the DC
  */
  double getDCActivityCoeff(const int dcidx) {
    if (dcidx < numDCs_) {
      return DCActivityCoeff_[dcidx];
    } else {
      EOBException ex("ChemicalSystem", "getDCActivityCoeff",
                      "DCActivityCoeff_", numDCs_, dcidx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the list of molal activity coefficients of each dependent component
  (DC).

  @note Used only in this class's copy constructor.

  @return a pointer to the list of DC molal activity coefficients
  */
  double *getDCActivityCoeff(void) const { return DCActivityCoeff_; }

  /**
  @brief Get the error (residual) in the mass balance for an independent
  component (IC).

  @note NOT USED.

  @param icidx is the index of the IC being queried
  @return the residual for the IC being queried
  */
  double getICResiduals(const int icidx) {
    if (icidx < numICs_) {
      return ICResiduals_[icidx];
    } else {
      EOBException ex("ChemicalSystem", "getICResiduals", "ICResiduals_",
                      numICs_, icidx);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the list of errors (residuals) in the mass balance for each
  independent component (IC).

  @note Used only in this class's copy constructor.

  @return a pointer to the list of IC residuals
  */
  double *getICResiduals(void) const { return ICResiduals_; }

  /**
  @brief Get the chemical potential of a dependent component (DC) [mol/mol or
  J/mol].

  @note NOT USED.

  @param dcidx is the index of the DC being queried
  @param norm is true if units will be in mol/mol, otherwise in J/mol
  @return the chemical potential of the DC
  */
  double getDCChemicalPotential(const long int dcidx, bool norm = false) {
    return (node_->Get_muDC(dcidx, norm));
  }

  /**
  @brief Get the chemical activity of a dependent component (DC).

  @param dcidx is the index of the DC being queried
  @return the chemical activity of the DC
  */
  double getDCActivity(const long int dcidx) { return (node_->Get_aDC(dcidx)); }

  /**
  @brief Get the chemical activity of a dependent component (DC) by name

  @param dcname is the name of the DC being queried
  @return the chemical activity of the DC
  */
  double getDCActivity(const string &dcname) {
    int dcidx = getDCId(dcname);
    if (dcidx < numDCs_) {
      return (node_->Get_aDC(dcidx));
    }
    return (0.0);
  }

  /**
  @brief Get the concentration of a dependent component (DC).

  The units of the returned value depend on the type of DC being queried:

      - Aqueous species [molal]
      - Gas species [partial pressure]
      - Surface complexes [mol/m<sup>2</sup>]
      - Species in other phases [mole fraction]

  @param dcidx is the index of the DC being queried
  @return the concentration of the DC in appropriate units
  */
  double getDCConcentration(const long int dcidx) {
    return (node_->Get_cDC(dcidx));
  }

  /**
  @brief Get the concentration of a dependent component (DC) by name

  The units of the returned value depend on the type of DC being queried:

      - Aqueous species [molal]
      - Gas species [partial pressure]
      - Surface complexes [mol/m<sup>2</sup>]
      - Species in other phases [mole fraction]

  Will return a value of zero if the DC does not exist

  @param dcname is the name of the DC component begin queried
  @return the concentration of the DC in appropriate units
  */
  double getDCConcentration(const string &dcname) {
    int dcidx = getDCId(dcname);
    if (dcidx < numDCs_) {
      return (node_->Get_cDC(dcidx));
    }
    return (0.0);
  }

  /**
  @brief Set the red channel in the rgb triplet for color of a microstructure
  phase.

  @note NOT USED.

  @param mpidx is the index of the microstructure phase
  @param val is the value of the red channel to set
  */
  void setRed(const int mpidx, const int rval) {
    if (mpidx >= (int)(color_.size())) {
      EOBException ex("ChemicalSystem", "setRed", "color_", color_.size(),
                      mpidx);
      ex.printException();
      exit(1);
    }
    color_[mpidx][0] = min(rval, COLORSATVAL);
    return;
  }

  /**
  @brief Get the red channel in the rgb triplet for color of a microstructure
  phase.

  @note NOT USED.

  @param mpidx is the index of the microstructure phase
  @return the value of the red channel for this phase's color
  */
  int getRed(const int mpidx) {
    if (mpidx >= (int)(color_.size())) {
      EOBException ex("ChemicalSystem", "getRed", "color_", color_.size(),
                      mpidx);
      ex.printException();
      exit(1);
    }
    return color_[mpidx][0];
  }

  /**
  @brief Set the green channel in the rgb triplet for color of a microstructure
  phase.

  @note NOT USED.

  @param mpidx is the index of the microstructure phase
  @param val is the value of the green channel to set
  */
  void setGreen(const int mpidx, const int rval) {
    if (mpidx >= (int)(color_.size())) {
      EOBException ex("ChemicalSystem", "setGreen", "color_", color_.size(),
                      mpidx);
      ex.printException();
      exit(1);
    }
    color_[mpidx][1] = min(rval, COLORSATVAL);
    return;
  }
  /**
  @brief Get the green channel in the rgb triplet for color of a microstructure
  phase.

  @note NOT USED.

  @param mpidx is the index of the microstructure phase
  @return the value of the green channel for this phase's color
  */
  int getGreen(const int mpidx) {
    if (mpidx >= (int)(color_.size())) {
      EOBException ex("ChemicalSystem", "getGreen", "color_", color_.size(),
                      mpidx);
      ex.printException();
      exit(1);
    }
    return color_[mpidx][1];
  }

  /**
  @brief Set the blue channel in the rgb triplet for color of a microstructure
  phase.

  @note NOT USED.

  @param mpidx is the index of the microstructure phase
  @param val is the value of the blue channel to set
  */
  void setBlue(const int mpidx, const int rval) {
    if (mpidx >= (int)(color_.size())) {
      EOBException ex("ChemicalSystem", "setBlue", "color_", color_.size(),
                      mpidx);
      ex.printException();
      exit(1);
    }
    color_[mpidx][2] = min(rval, COLORSATVAL);
    return;
  }

  /**
  @brief Get the blue channel in the rgb triplet for color of a microstructure
  phase.

  @note NOT USED.

  @param mpidx is the index of the microstructure phase
  @return the value of the blue channel for this phase's color
  */
  int getBlue(const int mpidx) {
    if (mpidx >= (int)(color_.size())) {
      EOBException ex("ChemicalSystem", "getBlue", "color_", color_.size(),
                      mpidx);
      ex.printException();
      exit(1);
    }
    return color_[mpidx][2];
  }

  /**
  @brief Set the rgb triplet for color of a microstructure phase.

  @note NOT USED.

  @param mpidx is the index of the microstructure phase
  @param cv is the vector of rgb values to set for that phase
  */
  void setColor(const int mpidx, vector<int> cv) {
    if (mpidx >= (int)(color_.size())) {
      EOBException ex("ChemicalSystem", "setColor", "color_", color_.size(),
                      mpidx);
      ex.printException();
      exit(1);
    }
    color_[mpidx] = cv;
    return;
  }

  /**
  @brief Get the rgb triplet for color of a microstructure phase.

  @param mpidx is the index of the microstructure phase
  @return the vector of rgb values defining this phase's color
  */
  vector<int> getColor(const int mpidx) {
    // try {
    return color_[mpidx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getColor", "color_", color_.size(),
    //                  mpidx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Set the rgb triplet for color of every microstructure phase.

  @note NOT USED.

  @param cv is the 2D matrix of rgb values to set each microstructure phase
  */
  void setColor(vector<vector<int>> cv) { color_ = cv; }

  /**
  @brief Get the 2D matrix of rgb triplets for color of every microstructure
  phase.

  @note Used only in this class's copy constructor.

  @return the 2D matrix of rgb values of every microstructure phase
  */
  vector<vector<int>> getColor(void) const { return color_; }

  /**
  @brief Set the grayscale value of a microstructure phase.

  The grayscale value of each microstructure phase is in proportion to its
  relative brightness in a backscattered electron image.

  @note NOT USED.

  @param mpidx is the index of the microstructure phase
  @param rval is the grayscale value to set for that microstructure phase
  */
  void setGrayscale(const int mpidx, const int rval) {
    // try {
    grayscale_[mpidx] = min(rval, COLORSATVAL);
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "setGrayscale", "grayscale_",
    //                  grayscale_.size(), mpidx);
    //  ex.printException();
    //  exit(1);
    //}
    // return;
  }

  /**
  @brief Get the grayscale value of a microstructure phase.

  The grayscale value of each microstructure phase is in proportion to its
  relative brightness in a backscattered electron image.

  @note NOT USED.

  @param mpidx is the index of the microstructure phase
  @return the grayscale value for the microstructure phase
  */
  int getGrayscale(const int mpidx) {
    // try {
    return grayscale_[mpidx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getGrayscale", "grayscale_",
    //                  grayscale_.size(), mpidx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Set the grayscale value of every microstructure phase.

  The grayscale value of each microstructure phase is in proportion to its
  relative brightness in a backscattered electron image.

  @note NOT USED.

  @param gv is the list of grayscale values to set each microstructure phase
  */
  void setGrayscale(vector<int> gv) { grayscale_ = gv; }

  /**
  @brief Get the list of grayscale values of every microstructure phase.

  @note Used only in this class's copy constructor.

  @return the list of grayscale values of every microstructure phase
  */
  vector<int> getGrayscale(void) const { return grayscale_; }

  /**
  @brief Get the map of of the vector index of the microstructure phases by
  name.

  The integer id of each microstructure phase is keyed to its name in this map,
  so one can "look up" a phase id by the name of that phase.

  @note Used only in this class's copy constructor.

  @return the microstructure phase lookup map (look up by name)
  */
  map<string, int> getMicroPhaseIdLookup(void) const {
    return microPhaseIdLookup_;
  }

  /**
  @brief Get the map of of the vector index of the independent components (IC)
  by name.

  The integer id of each independent component is keyed to its name in this map,
  so one can "look up" a IC id by the name of that IC.

  @note Used only in this class's copy constructor.

  @return the IC lookup map (look up by name)
  */
  map<string, int> getICIdLookup(void) const { return ICIdLookup_; }

  /**
  @brief Get the map of of the vector index of the dependent components (DC) by
  name.

  The integer id of each dependent component is keyed to its name in this map,
  so one can "look up" a DC id by the name of that DC.

  @note Used only in this class's copy constructor.

  @return the DC lookup map (look up by name)
  */
  map<string, int> getDCIdLookup(void) const { return DCIdLookup_; }

  /**
  @brief Get the map of of the vector index of the GEM CSD phases by name.

  The integer id of each phase defined in the GEM CSD is keyed to its name in
  this map, so one can "look up" a GEM phase id by the name of that phase.

  @note Used only in this class's copy constructor.

  @return the DC lookup map (look up by name)
  */
  map<string, int> getGEMPhaseIdLookup(void) const { return GEMPhaseIdLookup_; }

  /**
  @brief Get the class code of an independent component (IC).

  @note NOT USED.

  @param icidx is the index of the IC
  @return the class code of the IC
  */
  char getICClassCode(const int icidx) {
    // try {
    return ICClassCode_[icidx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getICClassCode", "ICClassCode_",
    //                  ICClassCode_.size(), icidx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the class codes of all independent components (IC).

  @note Used only in this class's copy constructor.

  @return the list of character class codes of the ICs
  */
  vector<char> getICClassCode(void) const { return ICClassCode_; }

  /**
  @brief Get the class code of a dependent component (DC).

  @param dcidx is the index of the DC
  @return the class code of the DC
  */
  char getDCClassCode(const int dcidx) {
    // try {
    return DCClassCode_[dcidx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getDCClassCode", "DCClassCode_",
    //                  DCClassCode_.size(), dcidx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the class codes of all dependent components (DC).

  @return the list of character class codes of the DCs
  */
  vector<char> getDCClassCode(void) const { return DCClassCode_; }

  /**
  @brief Get the class code of a phase defined in the GEM CSD.

  @param pidx is the index of a GEM phase
  @return the class code of the GEM phase
  */
  char getGEMPhaseClassCode(const unsigned int pidx) {
    // try {
    return GEMPhaseClassCode_[pidx];
    //} catch (out_of_range &oor) {
    //  EOBException ex("ChemicalSystem", "getGEMPhaseClassCode",
    //                  "GEMPhaseClassCode_", GEMPhaseClassCode_.size(), pidx);
    //  ex.printException();
    //  exit(1);
    //}
  }

  /**
  @brief Get the class codes of all phases defined in the GEM CSD.

  @return the list of character class codes of the GEM phases
  */
  vector<char> getGEMPhaseClassCode(void) const { return GEMPhaseClassCode_; }

  /**
  @brief Set the system temperature [K].

  @param val is the absolute temperature [K]
  */
  void setTemperature(const double val) { T_ = val; }

  /**
  @brief Get the system temperature [K].

  @return the absolute temperature [K]
  */
  double getTemperature(void) const { return T_; }

  /**
  @brief Set the system pressure [Pa].

  @note NOT USED.

  @param val is the pressure [Pa]
  */
  void setP(const double val) { P_ = val; }

  /**
  @brief Get the system pressure [Pa].

  @return the pressure [Pa]
  */
  double getP(void) const { return P_; }

  /**
  @brief Set the system volume [m<sup>3</sup>].

  @note NOT USED.

  @param val is the system volume [m<sup>3</sup>]
  */
  void setVs(const double val) { Vs_ = val; }

  /**
  @brief Get the system volume [m<sup>3</sup>].

  @note Used only in this class's copy constructor.

  @return the system volume [m<sup>3</sup>]
  */
  double getVs(void) const { return Vs_; }

  /**
  @brief Set the system mass [kg].

  @note NOT USED.

  @param val is the system mass [kg]
  */
  void setMs(const double val) { Ms_ = val; }

  /**
  @brief Get the system mass [kg].

  @note Used only in this class's copy constructor.

  @return the system mass [kg]
  */
  double getMs(void) const { return Ms_; }

  /**
  @brief Set the system Gibbs energy [J].

  @note NOT USED.

  @param val is the system Gibbs energy [J]
  */
  void setGs(const double val) { Gs_ = val; }

  /**
  @brief Get the system Gibbs energy [J].

  @note Used only in this class's copy constructor.

  @return the system Gibbs energy [J]
  */
  double getGs(void) const { return Gs_; }

  /**
  @brief Set the system enthalpy [J].

  @note NOT USED.

  @param val is the system enthalpy [J]
  */
  void setHs(const double val) { Hs_ = val; }

  /**
  @brief Get the system enthalpy [J].

  @note Used only in this class's copy constructor.

  @return the system enthalpy [J]
  */
  double getHs(void) const { return Hs_; }

  /**
  @brief Set the aqueous solution ionic strength [molal].

  @note NOT USED.

  @param isval is the aqueous solution ionic strength [molal]
  */
  void setIonicStrength(const double isval) { ionicStrength_ = isval; }

  /**
  @brief Get the aqueous solution ionic strength [molal].

  @note Used only in this class's copy constructor.

  @return the aqueous solution ionic strength [molal]
  */
  double getIonicStrength(void) const { return ionicStrength_; }

  /**
  @brief Set the aqueous solution pH.

  @note NOT USED.

  @param val is the aqueous solution pH
  */
  void setPH(const double val) { pH_ = val; }

  /**
  @brief Get the aqueous solution pH.

  @return the aqueous solution pH
  */
  double getPH(void) const { return pH_; }

  /**
  @brief Set the aqueous solution pe.

  @note NOT USED.

  @param val is the aqueous solution pe
  */
  void setPe(const double val) { pe_ = val; }

  /**
  @brief Get the aqueous solution pe.

  @note Used only in this class's copy constructor.

  @return the aqueous solution pe
  */
  double getPe(void) const { return pe_; }

  /**
  @brief Set the aqueous solution Eh [V].

  @note NOT USED.

  @param val is the aqueous solution Eh [V]
  */
  void setEh(const double val) { Eh_ = val; }

  /**
  @brief Get the aqueous solution Eh [V].

  @note Used only in this class's copy constructor.

  @return the aqueous solution Eh [V]
  */
  double getEh(void) const { return Eh_; }

  /**
  @brief Set the handle identification for the GEM-IPM node doing the
  calculations.

  @note NOT USED.

  @param val is the node handle identification number
  */
  void setNodeHandle(const long int val) { nodeHandle_ = val; }

  /**
  @brief Get the handle identification for the GEM-IPM node doing the
  calculations.

  @note Used only in this class's copy constructor.

  @return the node handle identification number
  */
  long int getNodeHandle(void) const { return nodeHandle_; }

  /**
  @brief Set the status of the GEM-IPM node doing the calculations.

  The node status can be:

      - 0 = No GEM re-calculation needed for this node
      - 1 = Need GEM calculation with simplex (automatic) initial approximation
  (IA)
      - 2 = OK after GEM calculation with simplex (automatic) IA
      - 3 = Bad (not fully trustful) result after GEM calculation with simplex
  (automatic) IA
      - 4 = Failure (no result) in GEM calculation with simplex (automatic) IA
      - 5 = Need GEM calculation with no-simplex (smart) IA
              the previous GEM solution (full DATABR lists only)
      - 6 = OK after GEM calculation with no-simplex (smart) IA
      - 7 = Bad (not fully trustful) result after GEM calculation with
  no-simplex (smart) IA
      - 8 = Failure (no result) in GEM calculation with no-simplex (smart) IA
      - 9 = Terminal error has occurred in GEMIPM2K (e.g., memory corruption).
              Restart is required.

  @note NOT USED.

  @param val is the node status number to set
  */
  void setNodeStatus(const long int val) { nodeStatus_ = val; }

  /**
  @brief Get the status of the GEM-IPM node doing the calculations.

  The node status can be:

      - 0 = No GEM re-calculation needed for this node
      - 1 = Need GEM calculation with simplex (automatic) initial approximation
  (IA)
      - 2 = OK after GEM calculation with simplex (automatic) IA
      - 3 = Bad (not fully trustful) result after GEM calculation with simplex
  (automatic) IA
      - 4 = Failure (no result) in GEM calculation with simplex (automatic) IA
      - 5 = Need GEM calculation with no-simplex (smart) IA
              the previous GEM solution (full DATABR lists only)
      - 6 = OK after GEM calculation with no-simplex (smart) IA
      - 7 = Bad (not fully trustful) result after GEM calculation with
  no-simplex (smart) IA
      - 8 = Failure (no result) in GEM calculation with no-simplex (smart) IA
      - 9 = Terminal error has occurred in GEMIPM2K (e.g., memory corruption).
              Restart is required.

  @note Used only in this class's copy constructor.

  @return the node status number
  */
  long int getNodeStatus(void) const { return nodeStatus_; }

  /**
  @brief Set the number of iterations executed by the GEM-IPM calculation.

  @note NOT USED.

  @param val is the number of iterations
  */
  void setIterDone(const long int val) { iterDone_ = val; }

  /**
  @brief Get the number of iterations executed by the GEM-IPM calculation.

  @note Used only in this class's copy constructor.

  @return the number of iterations executed by the GEM-IPM solver
  */
  long int getIterDone(void) const { return iterDone_; }

  /**
  @brief Set the number of times in a row that GEM_run has failed.

  Calculations of equilibrium state sometimes fail to converge.  It is
  possible that some small tweak in the composition could help the
  algorithm converge, so we allow up to `maxGEMFails_' attempts for
  convergence, each time tweaking the composition in some way
  before giving up and throwing an exception.

  @param ntimes is the number of consecutive times GEM_run has failed
  */
  void setTimesGEMFailed(const int ntimes) {
    timesGEMFailed_ = (ntimes >= 0) ? ntimes : 0;
    return;
  }

  /**
  @brief Get the number of consecutive times a call to GEM_run has failed.

  @return the number of consecutive times a call to GEM_run has failed.
  */
  int getTimesGEMFailed(void) const { return timesGEMFailed_; }

  /**
  @brief Set the pointer to the GEM TNode object doing the GEM-IPM calculations.

  @note NOT USED (possibly used in GEM3K library).

  @param np is a pointer to the GEM TNode object doing the GEM-IPM calculations
  */
  void setNode(TNode *np) { node_ = np; }

  /**
  @brief Get the pointer to the GEM TNode object doing the GEM-IPM calculations.

  @return a pointer to the GEM TNode object doing the GEM-IPM calculations
  */
  TNode *getNode(void) { return node_; }

  /**
  @brief Set the flag for the saturation state.

  @param satstate is the flag value
  */
  void setIsSaturated(const bool satstate) { isSaturated_ = satstate; }

  /**
  @brief Get the saturation state flag.

  @return the saturation state flag
  */
  bool isSaturated(void) const { return isSaturated_; }

  /**
  @brief Set the simulated time at which to implement the sulfate attack
  calculations [hours].

  @param sattack_time is the time at which to start sulfate attack [hours]
  */
  void setSulfateAttackTime(const double sattack_time) {
    sulfateAttackTime_ = sattack_time;
  }

  /**
  @brief Get the simulated time at which to start sulfate attack [hours].

  @note NOT USED.

  @return the simulated time at which to begin sulfate attack [hours]
  */
  double getSulfateAttackTime(void) const { return sulfateAttackTime_; }

  /**
  @brief Set the simulated time at which to implement the leaching calculations
  [hours].

  @param leach_time is the time at which to start leaching [hours]
  */
  void setLeachTime(const double leach_time) { leachTime_ = leach_time; }

  /**
  @brief Get the simulated time at which to start leaching [hours].

  @note NOT USED.

  @return the simulated time at which to begin leaching [hours]
  */
  double getLeachTime(void) const { return leachTime_; }

  /**
  @brief Formatted writing of some microstructure phase data to a stream.

  @note NOT USED.

  @param stream is the output stream to which to direct output
  */
  // void writeDb(ostream &stream);

  /**
  @brief Formatted writing of the name, id, and internal porosity of a
  microstructure phase.

  @param i is the index of the microstructure phase
  @param stream is the output stream to which to direct output
  */
  // void writeMember(const unsigned int i, ostream &stream);

  /**
  @brief Formatted writing of ChemicalSystem data to a stream.

  @param out is the output file stream to which to direct output
  */
  void writeChemSys(ofstream &out);

  /**
  @brief Calculates new equilibrium state of the GEM system and relate to
  microstructure.

  One of the main purposes of the [[ChemicalSystem]] class is to calculate
  changes to the state from one time increment to another.  This is accomplished
  by changing the chemical system definition to alter the quantitities
  of one or more independent components (IC), and then call GEM3K
  to calculate the new assemblage of DCs and phases that minimizes the
  Gibbs free energy subject to the IC quantities in the system.

  As a result, we need methods that will allow one to change the quantities
  of ICs and then to call GEM3K functions that will calculate the new
  state and return the values of the quantities of the DCs and phases.

  The method returns the GEM-IPM TNode status, which can be:

      - 0 = No GEM re-calculation needed for this node
      - 1 = Need GEM calculation with simplex (automatic) initial approximation
  (IA)
      - 2 = OK after GEM calculation with simplex (automatic) IA
      - 3 = Bad (not fully trustful) result after GEM calculation with simplex
  (automatic) IA
      - 4 = Failure (no result) in GEM calculation with simplex (automatic) IA
      - 5 = Need GEM calculation with no-simplex (smart) IA
              the previous GEM solution (full DATABR lists only)
      - 6 = OK after GEM calculation with no-simplex (smart) IA
      - 7 = Bad (not fully trustful) result after GEM calculation with
  no-simplex (smart) IA
      - 8 = Failure (no result) in GEM calculation with no-simplex (smart) IA
      - 9 = Terminal error has occurred in GEM (e.g., memory corruption).
              Restart is required.

  @todo Water in CSH gel porosity has a lower chemical potential than in bulk;
  could modify chemical potentials to capture lower reactivity when only gel
  water remains.

  @param time is the simulated time associated with this state [hours]
  @param isFirst is true if this is the first state calculation, false otherwise
  @param cyc is the cycle number for the main controller loop (iteration over
  time)
  @return the node status handle
  */
  int calculateState(double time, bool isFirst, int cyc);

  /**
  @brief Check for chemical composition requirements on
  the electrolyte and set the DC moles if necessary.

  @param isFirst is true if this is the first cycle
  */
  void setElectrolyteComposition(const bool isFirst, bool attack, int cyc);

  /**
  @brief Check for chemical composition requirements on
  the gas and set the DC moles if necessary.

  @param isFirst is true if this is the first cycle
  */
  void setGasComposition(const bool isFirst, bool attack);

  /**
  @brief Update the number of moles of each IC based on changes to a dependent
  component.

  @note NOT USED.

  @param dcid is the index of the DC
  @param moles is the change in number of moles of the DC
  */
  void DCImpact(int dcid, double moles) {
    for (int i = 0; i < numICs_; i++) {
      ICMoles_[i] += DCStoich_[dcid][i] * moles;
    }

    return;
  }

  /**
  @brief Check whether growth of a given dependent component is allowed.

  The change in moles of a particular dependent component (DC) could conceivably
  cause the moles of one of the independent components to decrease to
  effectively zero.  This is not allowed for numerical stability reasons, so
  this function checks whether that would happen.

  @note NOT USED.

  @param dcid is the index of the DC that wants to change in number moles
  @param moles is the proposed change in moles of that DC
  @return true if the change in this DC would not cause an IC mole content to
  drop to near zero
  */
  bool checkDCImpact(int dcid, double moles) {
    bool possible;
    possible = true;
    for (int i = 0; i < numICs_; i++) {
      if ((ICMoles_[i] + DCStoich_[dcid][i] * moles) < 2.0e-17) {
        possible = false;
        cout << "The growth of this phase can cause one or more IC moles to be"
             << " lower than 2.0e-17, so this phase can not grow." << endl;
        break;
      }
    }
    return possible;
  }

  /**
  @brief Set the vector of saturation indices of all GEM CSD phases.
  */
  void setSI(void) {
    SI_.clear();
    double *Falp;
    Falp = (node_->ppmm())->Falp;
    for (int i = 0; i < numGEMPhases_; i++) {
      double si = pow(10, Falp[i]);
      SI_.push_back(si);
      // if (verbose_) {
      //  cout << "logSI for i = " << i << " (" << GEMPhaseName_[i] << ") is:
      //  Falp[i] = "
      //       << Falp[i] << "   &   SI_[i] = " << SI_[i] << endl;
      // }
      if (GEMPhaseName_[i] == "ettr") {
        cout << endl << "ChemicalSystem::setSI : " << endl;
        cout << "\t\t" << i << "\t" << GEMPhaseName_[i] << "\tFalp: " << Falp[i]
             << "\tPh_SatInd: " << node_->Ph_SatInd(i)
             << "\tSI_(Falp): " << SI_[i] << endl;
      }
    }

    // cout << endl << "ChemicalSystem::setSI : " << endl;
    // for (int i = 0; i < numGEMPhases_; i++) {
    //   cout << "\t: " << i << "\tFalp: " << Falp[i]<< "\tSI_: " << SI_[i]
    //        << "\tPh_SatInd: " << node_->Ph_SatInd(i) //pow(10,
    //        node_->Ph_SatInd(i))
    //   //          << "\tFalp: " << Falp[i]
    //        << "\t\t" << GEMPhaseName_[i] << endl;
    // }
    // cout << "ChemicalSystem::setSI : end" << endl;

    return;
  }

  /**
  @brief Get the vector of saturation indices of all GEM CSD phases.

  @return the vector of saturation indices of all GEM CSD phases
  */
  vector<double> getSI(void) { return SI_; }

  /**
  @brief Get the saturation index of a GEM CSD phase, by its id.

  @param phaseid is the index of the GEM phase to query
  @return the saturation index of the GEM phase
  */
  double getSI(int phaseid) {
    try {
      return SI_.at(phaseid);
    } catch (out_of_range &oor) {
      EOBException ex("ChemicalSystem", "getSI", "SI_", SI_.size(), phaseid);
      ex.printException();
      exit(1);
    }
  }

  void initSI(int phaseid, double val) {
    if (SI_.size() == 0)
      SI_.resize(numGEMPhases_, 0.0);
    SI_[phaseid] = val;
  }

  /**
  @brief Get the saturation index of a GEM CSD phase, by its name.

  @param str is the name of the GEM phase to query
  @return the saturation index of the GEM phase
  */
  double getSI(const string &str) { return SI_[getGEMPhaseId(str)]; }

  /**
  @brief Set the vector of saturation indices of all microstructure phases.

  */
  void setMicroPhaseSI(void);
  void setMicroPhaseSI(double t);

  /**
  @brief Set the vector of saturation indices of all microstructure phases.

  @param cyc is the cycle number for time
  */
  // void setMicroPhaseSI(int cyc);

  /**
  @brief Get the vector of saturation indices of all microstructure phases.

  @return the vector of saturation indices of all microstructure phases
  */
  vector<double> getMicroPhaseSI(void) { return microPhaseSI_; }

  /**
  @brief Get the saturation index of a microstructure phase, by its id.

  @param microphaseid is the index of the microstructure phase to query
  @return the saturation index of the microstructure phase
  */
  double getMicroPhaseSI(int microphaseid) {
    // try {
    // cout << "Trying to find microstructure phase id " << microphaseid <<
    // endl; cout.flush();
    return microPhaseSI_[microphaseid];
    //} catch (out_of_range &oor) {
    //    EOBException ex("ChemicalSystem", "getMicroPhaseSI", "microPhaseSI_",
    //                    microPhaseSI_.size(), microphaseid);
    //    ex.printException();
    //    exit(1);
    //}
  }

  /**
  @brief Get the saturation index of a microstructure phase, by its name.

  @param str is the name of the microstructure phase to query
  @return the saturation index of the microstructure phase
  */
  double getMicroPhaseSI(const string &str) {
    return microPhaseSI_[getMicroPhaseId(str)];
  }

  /**
  @brief Get the activity of a GEM DC by its id.

  @param dcstr is the name of the GEM DC to query
  @return the activity of the DC
  */
  double getActivity(const string &dcstr) {
    // try {
    return node_->DC_a(getDCId(dcstr));
    //} catch (out_of_range &oor) {
    //   EOBException ex("ChemicalSystem", "getActivity", dcstr, numDCs_,
    //                   getDCId(dcstr));
    //   ex.printException();
    //   exit(1);
    // }
  }

  /**
  @brief Get the current list of IC moles in the aqueous solution.

  @return the vector of moles of each IC in the aqueous solution
  */
  vector<double> getSolution(void);

  /**
  @brief Get the initial solution composition other than water

  This function returns the map of molal concentrations of each IC in
  the initial solution [mol/kgw].

  @return the initial solute concentration map
  */
  map<int, double> getInitialSolutionComposition(void) {
    return initialSolutionComposition_;
  }

  /**
  @brief Get the fixed solution composition other than water

  This function returns the map of molal concentrations of each DC in
  the fixed solution [mol/kgw].

  @return the fixed solute concentration map
  */
  map<int, double> getFixedSolutionComposition(void) {
    return fixedSolutionComposition_;
  }

  /**
  @brief Get the initial gas composition

  This function returns the map of molal concentrations of each DC in
  the initial solution [mol/kgw].

  @return the initial solute concentration map
  */
  map<int, double> getInitialGasComposition(void) {
    return initialGasComposition_;
  }

  /**
  @brief Get the fixed gas composition

  This function returns the map of molal concentrations of each DC in
  the fixed solution [mol/kgw].

  @return the fixed solute concentration map
  */
  map<int, double> getFixedGasComposition(void) { return fixedGasComposition_; }

  /**
  @brief Set the gas-solid mass ratio

  */
  void setGasSolidRatio(const double gassolidratio) {
    gasSolidRatio_ = (gassolidratio > 0.0) ? gassolidratio : 0.0;
  }

  /**
  @brief Get the gas-solid mass ratio

  @return the mass ratio of gas to initial solids
  */
  double getGasSolidRatio(void) { return gasSolidRatio_; }

  /**
  @brief Set the jsonFormat_ flag

  @param jsonFormat is true if input files are in JSON format
  */
  void setJSONFormat(const bool jsonFormat) {
    jsonFormat_ = jsonFormat;
    return;
  }

  /**
  @brief Get the jsonFormat_ flag

  @return the jsonFormat flag
  */
  bool isJSONFormat(void) const { return jsonFormat_; }

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
  bool getVerbose(void) const { return verbose_; }

  /**
  @brief Set the warning flag

  @param isvsarning is true if warning messages should be produced
  */
  void setWarning(const bool iswarning) {
    warning_ = iswarning;
    return;
  }

  /**
  @brief Get the warning flag

  @return the warning flag
  */
  bool getWarning(void) const { return warning_; }

  //*@*********************************************

  void checkChemSys(void);

  /**
  @brief Get the list of all GEM CSD phases that are associated with a given
  microstructure phase id number.

  @note NOT USED.

  @param idx is the microstructure phase in question
  @return the vector of all GEM CSD phase ids associated with the microstructure
  phase
  */
  vector<int> getMicroPhaseMembers(const unsigned int idx) {
    string msg;
    map<int, vector<int>>::iterator p = microPhaseMembers_.find(idx);
    if (p != microPhaseMembers_.end()) {
      return p->second;
    } else {
      msg = "Could not find microPhaseMembers_ match to index provided";
      EOBException ex("ChemicalSystem", "getMicroPhaseMembers", msg,
                      microPhaseMembers_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  /**
  @brief Get the map of of the vector index of the microstructure phases by
  name.

  The integer id of each microstructure phase is keyed to its name in this map,
  so one can "look up" a phase id by the name of that phase.

  @note Used only in this class's copy constructor.

  @return the microstructure phase lookup map (look up by name)
  */
  int getMicroPhaseIdLookup(string str) {
    string msg;
    map<string, int>::iterator p = microPhaseIdLookup_.find(str);
    if (p != microPhaseIdLookup_.end()) {
      return p->second;
    } else {
      msg = "Could not find microPhaseIdLookup_ match to string provided";
      EOBException ex("ChemicalSystem", "getMicroPhaseIdLookup", msg,
                      microPhaseIdLookup_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  int getICIdLookup(string str) // ICId -> xCH
  {
    string msg;
    map<string, int>::iterator p = ICIdLookup_.find(str);
    if (p != ICIdLookup_.end()) {
      return p->second;
    } else {
      msg = "Could not find ICIdLookup_ match to string provided";
      EOBException ex("ChemicalSystem", "getICIdLookup", msg,
                      ICIdLookup_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  int getDCIdLookup(string str) // DCId -> xCH
  {
    string msg;
    map<string, int>::iterator p = DCIdLookup_.find(str);
    if (p != DCIdLookup_.end()) {
      return p->second;
    } else {
      msg = "Could not find DCIdLookup_ match to string provided";
      EOBException ex("ChemicalSystem", "getDCIdLookup", msg,
                      DCIdLookup_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  int getGEMPhaseIdLookup(string str) {
    string msg;
    map<string, int>::iterator p = GEMPhaseIdLookup_.find(str);
    if (p != GEMPhaseIdLookup_.end()) {
      return p->second;
    } else {
      msg = "Could not find GEMPhaseIdLookup_ match to string provided";
      EOBException ex("ChemicalSystem", "getGEMPhaseIdLookup", msg,
                      GEMPhaseIdLookup_.size(), 0);
      ex.printException();
      exit(1);
    }
  }

  string getMicroPhaseName(int i) { return microPhaseName_[i]; }

  double getDCMolesNode(const unsigned int idx) { return node_->DC_n(idx); }

  bool getIsDCKinetic(int i) { return isDCKinetic_[i]; }

  void setIsDCKinetic(int i, bool val) { isDCKinetic_[i] = val; }

  void setDC_to_MPhID(int i, int val) { DC_to_MPhID_[i] = val; }

  int getDC_to_MPhID(int i) { return DC_to_MPhID_[i]; }

  void addWaterMassAndVolume(double massVal, double volVal) {

    // int wMPhID = getMicroPhaseId("Electrolyte");
    // cout << "wMPhID = " << wMPhID << endl;

    microPhaseMass_[1] += massVal;
    microPhaseVolume_[1] += volVal;

    return;
  }

  bool isCementComponent(int i) { return cementComponent_[i]; }

  void setInitScaledCementMass(double val) {
    initScaledCementMass_ = scaledCementMass_ = val;
  }

  double getInitScaledCementMass(void) { return initScaledCementMass_; }

  double getScaledCementMass(void) { return scaledCementMass_; }

  void setZeroMicroPhaseSI(void) { microPhaseSI_.resize(numMicroPhases_, 0.0); }

  void initColorMap(void);

  vector<int> getRGB(int pid) {
    string mPhName = microPhaseName_[pid];
    map<string, elemColor>::iterator p = colorN_.find(mPhName);
    if (p != colorN_.end()) {
      return colorN_[mPhName].rgb;
    } else {
      cout << endl << "**********************************************" << endl;
      cout << endl
           << "   Microphase " << mPhName
           << " has no associated rgb values by default!" << endl;
      cout << endl << "   => program stops !" << endl;
      cout << endl
           << "Please add in the chemistry.json file before " << mPhName
           << " close phase definition tag (</phase>)," << endl
           << "the following lines replacing VALUE with convenient integer "
              "numbers in [0,255]: "
           << endl;
      cout << endl << "<display_data>" << endl;
      cout << " <red> VALUE </red>" << endl;
      cout << " <green> VALUE </green>" << endl;
      cout << " <blue> VALUE </blue>" << endl;
      cout << " <gray> VALUE </gray>" << endl;
      cout << "</display_data>" << endl;
      cout << endl
           << "The following microphaseses are defined by default in THAMES "
              "3.0.0: "
           << endl;
      int i = 0;
      for (map<string, elemColor>::iterator pp = colorN_.begin();
           pp != colorN_.end(); pp++) {
        cout << "   " << setw(3) << i << " : " << setw(15) << left << pp->first
             << setw(5) << right << "rgb:" << setw(5) << pp->second.rgb[0]
             << setw(5) << pp->second.rgb[1] << setw(5) << pp->second.rgb[2]
             << setw(10) << "gray:" << setw(5) << pp->second.gray << endl;
        i++;
      }
      cout << endl
           << "After modifying and saving the chemistry.json file, please "
              "restart the program."
           << endl
           << endl;
      exit(0);
    }
  }

  vector<float> getRGBf(int pid) {
    string mPhName = microPhaseName_[pid];
    map<string, elemColor>::iterator p = colorN_.find(mPhName);
    if (p != colorN_.end()) {
      return colorN_[mPhName].rgbf;
    } else {
      cout << endl << "**********************************************" << endl;
      cout << endl
           << "   Microphase " << mPhName
           << " has no associated rgb values by default!" << endl;
      cout << endl << "   => program stops !" << endl;
      cout << endl
           << "Please add in the chemistry file before " << mPhName
           << " close phase definition tag (</phase>)," << endl
           << "the following lines replacing VALUE with convenient integer "
              "numbers in [0,255]: "
           << endl;
      cout << endl << "<display_data>" << endl;
      cout << " <red> VALUE </red>" << endl;
      cout << " <green> VALUE </green>" << endl;
      cout << " <blue> VALUE </blue>" << endl;
      cout << " <gray> VALUE </gray>" << endl;
      cout << "</display_data>" << endl;
      cout << endl
           << "The following microphaseses are defined by default in THAMES "
              "3.0.0: "
           << endl;
      int i = 0;
      for (map<string, elemColor>::iterator pp = colorN_.begin();
           pp != colorN_.end(); pp++) {
        cout << "   " << setw(3) << i << " : " << setw(15) << left << pp->first
             << setw(5) << right << "rgb:" << setw(5) << pp->second.rgb[0]
             << setw(5) << pp->second.rgb[1] << setw(5) << pp->second.rgb[2]
             << setw(10) << "gray:" << setw(5) << pp->second.gray << endl;
        i++;
      }
      cout << endl
           << "After modifying and saving the chemistry.json file, please "
              "restart the program."
           << endl
           << endl;
      exit(0);
    }
  }

  void updateMicroPhaseMasses(int idx, double val, int called) {
    int DCId = 0;
    if (idx > ELECTROLYTEID) {
      microPhaseMassDissolved_[idx] = microPhaseMass_[idx] - val;
      microPhaseMass_[idx] = val;

      DCId = getMicroPhaseDCMembers(idx, 0);
      double v0 = node_->DC_V0(DCId, P_, T_);
      double dcmm = getDCMolarMass(DCId);
      if (dcmm < 1.0e-9) {
        FloatException fex("ChemicalSystem", "updateKCMicroPhaseMassess",
                           "Divide by zero (dcmm)");
        fex.printException();
        exit(1);
      }
      // setMicroPhaseVolume(idx, (val * v0 / dcmm));
      microPhaseVolume_[idx] = val * v0 / dcmm;
      if (verbose_) {
        if (called == 0) {
          cout << "    ChemicalSystem::updateMicroPhaseMassess for idx = "
               << setw(3) << right << idx << " : " << setw(15) << left
               << microPhaseName_[idx]
               << " (called = 0) => updated scaledMass = " << val
               << " and volume = " << microPhaseVolume_[idx] << endl;
        } else {
          cout << "    ChemicalSystem::updateMicroPhaseMassess for idx = "
               << setw(3) << right << idx << " : " << setw(15) << left
               << microPhaseName_[idx]
               << " (called = 1) => updated scaledMass = " << val
               << " and volume = " << microPhaseVolume_[idx] << endl;
        }
        cout.flush();
      }
    } else {
      cout << endl
           << "   error in ChemicalSystem::setKCMicroPhaseMasses : idx = "
           << idx << endl;
      cout << endl << "   exit" << endl;
      exit(1);
    }
  }

  /**
  @brief add a component to the isParrotKilloh_ vector
  @param val is the microstructure phase ID of this PK component
  */
  void setIsParrotKilloh(int val) {
    isParrotKilloh_.push_back(val);
    isParrotKillohSize_++;
  }

  /**
  @note NOT USED.

  @return the isParrotKilloh_ vector idx component i.e. a microPhaseId of a
  microphase controlled by Parrot-Killoh model
  */
  // int getIsParrotKilloh(int idx) { return isParrotKilloh_[idx]; }

  /**
  @note NOT USED.

  @return the isParrotKilloh_ vector : all microPhaseIds for microphases
  controlled by Parrot-Killoh model
  */
  vector<int> getIsParrotKilloh(void) { return isParrotKilloh_; }

  /**
  @brief Write out the boundary conditions on the electrolyte
  */
  void writeSatElectrolyteGasConditions(void);

  /**
  @brief Calculate the site strain based on local crystallization pressure.

  The calculation of local crystallization pressure happens in this method,
  based on the saturation index of the phase.  Therefore, the method is
  currently restricted to problems of sulfate attack, and it probably can
  be generalized somewhat.

  Assuming that the saturation index \f$\beta\f$ is known for the phase,
  then the crystallization pressure should be the difference in the Laplace
  pressure between the large pore entrance, \f$r_{cp}\f$, and the size of the
  average gel pore, \f$r_{gp}\f$.  This pressure difference is

  \f[
      p_c = 2 \gamma \left[ \frac{1}{r_{gp} - \delta} - \frac{1}{r_{cp} -
  \delta} \right] \f]

  where we subtract \f$\delta\f$, the liquid film thickness, so that the terms
  in the denominator are the radii of curvature of the actual crystal in these
  two locations.

  In THAMES, we <i>assume</i> that the largest pore entrance to the gel porosity
  is about half the size of a lattice site.  The usual dimension of a lattice
  site in THAMES is 1 micrometer--- <b>although this is not necessary</b>--- so
  we assume that \f$r_{cp} = 500\f$ nm.

  The Thompson-Freundlich relation tells us the size of a crystal that is in
  equilibrium with a supersaturated solution with saturation index, relative to
  that crystal, of \f$\beta\f$.  The condition of equilibrium is that the
  chemical potential of formula units in the crystal be equal to the chemical
  potential of dissolved formula units in the solution:

  \f[
      V_c \gamma \kappa_{cr} = R T \ln \frac{Q}{K} = R T \ln \beta
  \f]

  where \f$V_c\f$ is the stress-free molar volume of the crystal, \f$\gamma\f$
  is the surface energy of the crystal-liquid interface, \f$\kappa_{cr}\f$ is
  the mean curvature of the crystal in equilibrium with the solution, <i>R</i>
  and <i>T</i> are the gas constant and absolute temperature, respectively, and
  <i>Q</i> and <i>K</i> are the activity product and equilibrium constant of the
  crystal dissociation reaction.  Therefore, the mean curvature is

  \f[
      \kappa_{cr} = \frac{R T}{V_c \gamma} \ln \beta
  \f]

  And since \f$\kappa_{cr} \equiv 2/r_{cr}\f$, we have \f$r_{cr} =
  2/\kappa_{cr}\f$ But \f$r_{gp}\f$ is the radius of curvature of the gel pore,
  not the crystal, so we must add the liquid film thickness: \f$r_{gp} = r_{cr}
  + \delta\f$.

  With the crystallization pressure, <i>p</i><sub>c</sub>, calculated, we must
  now estimate the local strain.  Ordinarily, this could be done only with
  detailed knowledge of the local pore structure.  However, to make the
  calculation tractable, we adopt the assumptions of poromechanics, in which
  case the strain is estimated by a volume-averaged equation,

  \f[
      \epsilon_x \approx \left( \frac{1}{3 K_p} - \frac{1}{3 K_s} \right)
      \left[ \phi p_c + pl - p_{atm} \right]
  \f]

  where \f$\epsilon_x\f$ is the linear strain (strain is assumed to be
  isotropic), <i>K</i><sub>p</sub> and <i>K</i><sub>s</sub> are the volume
  averaged bulk and shear moduli of the porous body, \f$\phi\f$ is the local
  porosity, <i>p</i><sub>l</sub> is the hydrostatic pressure in the liquid film,
  and <i>p</i><sub>atm</sub> is atmospheric pressure.

  @todo Consider renaming this function to calcCrystallizationStrain.

  @param SI is the saturation index of the growing crystal
  @param poreVolFrac is the local porosity
  @param Kp is the effective bulk modulus of the porous body
  @param Ks is the effective shear modulus of the porous body
  @return the calculated crystallization strain
  */
  double calculateCrystalStrain(int growPhId, double poreVolFrac, double Kp,
                                double Ks);

  void setIniAttackTime(const double val) { beginAttackTime_ = val; }

  void initElasticModuliMap(void);

  elMod getElasticModuliComp(string str) { return elasticModuli_[str]; }

}; // End of ChemicalSystem class
#endif
