/**
@file global.h
@brief Declare and assign global variables for THAMES.

*/
#ifndef GLOBALH
#define GLOBALH

// #ifndef DEBUG
// #define DEBUG
// #endif

#include "Exceptions.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace std;

// String indicating THAMES version in the input microstructure file
const string VERSIONSTRING = "#THAMES:Version:";

const string VERSIONBUGFIX = "0";

// String indicating the voxel resolution in micrometers
// for the THAMES input microstructure file
const string IMGRESSTRING("#THAMES:Image_Resolution:");

// String indicating the x dimension in voxels
// for the THAMES input microstructure file
const string XSIZESTRING("#THAMES:X_Size:");

// String indicating the y dimension in voxels
// for the THAMES input microstructure file
const string YSIZESTRING("#THAMES:Y_Size:");

// String indicating the z dimension in voxels
// for the THAMES input microstructure file
const string ZSIZESTRING("#THAMES:Z_Size:");

// String indicating the Alite DC name
const string AliteDCName("C3S");

// String indicating the Belite DC name
const string BeliteDCName("C2S");

// String indicating the Aluminate DC name
const string AluminateDCName("C3A");

// String indicating the Ferrite DC name
const string FerriteDCName("C4AF");

// String indicating the liquid solution GEM Phase name
// @todo Make this general somehow
const string WaterGEMName("aq_gen");

// String indicating the liquid water DC name
// @todo Make this general somehow
const string WaterDCName("H2O@");

// String indicating the CSH GEM Phase name
// @todo Make this general somehow
const string CSHGEMName("CSHQ");

// String indicating the AFt GEM Phase name
// @todo Make this general somehow
// const string AFTDCName("ettringite");
const string AFTGEMName("ettr");

// String indicating the AFt DC name
// @todo Make this general somehow
const string AFTDCName("ettr");

// String indicating the monosulfate GEM Phase name
// @todo Make this general somehow
const string MonosulfGEMName("C4AsH14");

// String indicating the monosulfate DC name
// @todo Make this general somehow
// const string MonosulfGEMName("C4AsH14");
const string MonosulfDCName("monosulf14");

// String indicating the monocarboaluminate GEM Phase name
// @todo Make this general somehow
const string MonocarbGEMName("C4AcH11");

// String indicating the hydrotalcite GEM Phase name
// @todo Make this general somehow
// const string HydrotalcGEMName("OH-hydrotalcite");
const string HydrotalcGEMName("hydrotalc-pyro");

// String indicating a generic kinetic model
const string GenericType("Generic");

// String indicating the Parrot Killoh Model
const string ParrotKillohType("ParrotKilloh");

// String indicating the Standard Dissolution Model
const string StandardType("Standard");

// String indicating the Pozzolanic Reaction Model
const string PozzolanicType("Pozzolanic");

// String for the command to convert between image files (Imagemagick)
const string ConvertCommand("magick");

// Flag to indicate exiting the program
const int QUIT_PROGRAM = 1;

// Flag to indicate simulation of hydration only
const int HYDRATION = 2;

// Flag to indicate simulation of leaching only
const int LEACHING = 3;

// Flag to indicate simulation of external sulfate only
const int SULFATE_ATTACK = 4;

// Flag to indicate kinetic equations handled within code (NOT USED)
const int INTERNAL_KINETICS = 0;

// Flag to indicate kinetics is done a priori outside the code (NOT USED)
const int EXTERNAL_KINETICS = 1;

// Normal return condition flag (NOT USED)
const int RETURN_NORMAL = 0;

// Flag to indicate an element is out of bounds (NOT USED)
const int DB_EOB = 1;

// Flag to indicate call to GEM failed (NOT USED)
const int GEMRUN_ERROR = 1;

// Flag to indicate that a lattice element is out of bounds (NOT USED)
const int MESH_EOB = 1;

// Flag to indicate that an array element is out of bounds (NOT USED)
const int ARRAY_EOB = 1;

// Flag to indicate that a file could not be opened (NOT USED)
const int FILE_OPEN_ERROR = 2;

// Flag to indicate unexpected end of file (NOT USED)
const int PREMATURE_EOF = 3;

// Flag to indicate datum is bad or of wrong kind
const int INVALID_INPUT = 4;

// Special phase ids that are important and must always be the same value
const int VOIDID = 0;
const int ELECTROLYTEID = 1;
const int FIRST_SOLID = ELECTROLYTEID + 1;
const int NUMCLINKERPHASES = 4;
const double thrPorosityCSH = 0.0355255;

// The number of face, edge, and corner neighbors to a cubic lattice site
const int NUM_NEAREST_NEIGHBORS = 6;
const int NUM_SECONDNEAREST_NEIGHBORS = 12;
const int NUM_THIRDNEAREST_NEIGHBORS = 8;

const int NN_NNN = NUM_NEAREST_NEIGHBORS + NUM_SECONDNEAREST_NEIGHBORS;

// Maximum allowed std::string length (NOT USED)
const int MAXSTRING = 128;

// Reference temperature for kinetic calculations [K]
const double REFTEMP = 298.15;

// Reference lattice resolution [micrometers]
const double REFRES = 4.0;

// Ideal gas constant [J/mol/K]
const double GASCONSTANT = 8.314;

// Flags for different kinds of probability distributions (NOT USED)
const int DELTA = +1;
const int UNIFORM = +2;
const int GAUSSIAN = +3;
const int LOGNORMAL = +4;

//
// Next constant specifies boundary condition configuration
// 0 = periodic boundaries everywhere
// 1 = periodic in y and z, greased in x
// 2 = periodic in x and z, greased in y
// 3 = periodic in x and y, greased in z
//

const int BC = 0;

// Growth mode constants
const unsigned int DLA = 0;
const unsigned int WMC = 1;
const unsigned int IWMC = 2;

// Saturation value for colors
const int COLORSATVAL = 255;

// Flags for data formats (NOT USED)
const int GEMSFORMAT = 0;
const int EXTERNALFORMAT = 0;

const double ICTHRESH = 1.0e-9;

const double Pi = 3.14159265359;

const double MIN_PER_S = 0.0166666667;
const double S_PER_MIN = 60.0;
const double H_PER_MIN = 0.0166666667;
const double MIN_PER_H = 60.0;
const double S_PER_H = 3600.0;
const double H_PER_DAY = 24.0;
const double DAY_PER_H = 0.0416666667;
const double H_PER_S = 0.0002777778;
const double DAY_PER_Y = 365.0;
const double Y_PER_DAY = 0.002739726;
const double S_PER_Y = 3.153600000e7;
const double S_PER_DAY = 86400.0;

typedef struct {
  int years;
  int days;
  int hours;
  int minutes;
} TimeStruct;

#endif // GLOBALH
