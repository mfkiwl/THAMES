/**
@struct KineticData
@brief Stores data about each phase possible in the system for ease of parsing
the input files.

In THAMES, phases are identified either thermodynamically--- in the
GEM data repository--- or microstructurally.  A microstructural phase can
be one, or a combination of more than on, thermodynamically defined phase.

The structure is defined to make it easier to parse the input file for the
kinetic model. It is not used elsewhere in the code.  In fact,
the same members are identified as class variables in the `KineticModel` class.

Most of the members have self-evident meanings:
    - `name` is the name of the microstructure phase
    - `microPhaseId` is the integer id for the phase in the microstructure
    - `GEMPhaseId` is the integer id of the same phase in the GEM Chemical
System Definition (CSD)
    - `DCId` is the integer id of the GEM dependent component making up the
phase
    - `type` is a string specifying whether the phase is under kinetic control
        or thermodynamic control
    - `scaledMass` is the number of grams of the phase per 100 grams of solid
    - `specificSurfaceArea` is the specific surface area in m2 per kg
    - `refSpecificSurfaceArea` is the reference specific surface area in m2 per
kg
    - `k1` is the Parrot and Killoh <i>K</i><sub>1</sub> parameter for the phase
    - `k2` is the Parrot and Killoh <i>K</i><sub>2</sub> parameter for the phase
    - `k3` is the Parrot and Killoh <i>K</i><sub>3</sub> parameter for the phase
    - `n1` is the Parrot and Killoh <i>N</i><sub>1</sub> parameter for the phase
    - `n3` is the Parrot and Killoh <i>N</i><sub>3</sub> parameter for the phase
    - `dissolutionRateConst` is a generic rate constant for dissolution or
growth (flux units)
    - `diffusionRateConstEarly` is a generic rate constant for early-age
diffusion (flux units)
    - `diffusionRateConstLate` is a generic rate constant for later-age
diffusion (flux units)
    - `dissolvedUnits` is the number of DC units produced by unit dissolution
reaction
    - `siexp` is the exponent on the saturation index in the rate equation
    - `dfexp` is the exponent on the driving force term in the rate equation
    - `dorexp` is the exponent on the degree of reaction in the diffusion rate
equations
    - `ohexp` is the exponent on the hydroxyl ion activity in the rate equation
    - `sio2' is the mole fraction of SiO<sub>2</sub> in the component
    - `al2o3' is the mole fraction of Al<sub>2</sub>O<sub>3</sub> in the
component
    - `cao' is the mole fraction of CaO in the component
    - `activationEnergy` is the activation energy [J/mol/K]
    - `critDOH` is the critical degree of hydration used in the equation for
        calculating the influence of w/c ratio.
*/

#ifndef KINETICDATASTRUCT
#define KINETICDATASTRUCT
struct KineticData {
  string name;                /**< Name of the microstructure phase */
  int microPhaseId;           /**< Integer id of the microstructure phase */
  int GEMPhaseId;             /**< Integer id of the phase in the GEM CSD */
  int DCId;                   /**< Integer id of the DC making up the phase */
  string type;                /**< Specifies kinetic or thermodynamic control */
  double scaledMass;          /**< Mass percent on a total solids basis */
  double specificSurfaceArea; /**< Specific surface area [m2/kg] */
  double refSpecificSurfaceArea; /**< Reference specific surface area [m2/kg] */
  double temperature;            /**< Temperature [K] */
  double reftemperature;         /**< Reference temperature [K] */
  double activationEnergy;       /**< Activation energy [J/mol/ */
  double k1;      /**< Parrot and Killoh <i>K</i><sub>1</sub> parameter */
  double k2;      /**< Parrot and Killoh <i>K</i><sub>2</sub> parameter */
  double k3;      /**< Parrot and Killoh <i>K</i><sub>3</sub> parameter */
  double n1;      /**< Parrot and Killoh <i>N</i><sub>1</sub> parameter */
  double n3;      /**< Parrot and Killoh <i>N</i><sub>3</sub> parameter */
  double critDOH; /**< Critical degree of hydration for w/c effect */
  double dissolutionRateConst;    /**< Generic rate constant [mol/m2/s] */
  double diffusionRateConstEarly; /**< Generic rate constant [mol/m2/s] */
  double diffusionRateConstLate;  /**< Generic rate constant [mol/m2/s] */
  double dissolvedUnits; /**< Number of DC units produced by unit dissolution */
  double siexp;          /**< Exponent on saturation index [unitless] */
  double dfexp;          /**< Exponent on driving force [unitless] */
  double dorexp;         /**< Exponent on degree of reaction [unitless] */
  double ohexp;          /**< Exponent on OH ion activity [unitless] */
  double loi;            /**< Loss on ignition [unitless] */
  double sio2;           /**< Mole fraction SiO2 in material */
  double al2o3;          /**< Mole fraction Al2O3 in material */
  double cao;            /**< Mole fraction CaO in material */
};
#endif
