/***
@brief Header file for viz program

*/
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Preprocessor defines
// VCCTL phase codes
const int POROSITY = 0;            /* 0 */
const int C3S = POROSITY + 1;      /* 1 */
const int C2S = C3S + 1;           /* 2 */
const int C3A = C2S + 1;           /* 3 */
const int C4AF = C3A + 1;          /* 4 */
const int K2SO4 = C4AF + 1;        /* 5 */
const int NA2SO4 = K2SO4 + 1;      /* 6 */
const int GYPSUM = NA2SO4 + 1;     /* 7 */
const int HEMIHYD = GYPSUM + 1;    /* 8 */
const int ANHYDRITE = HEMIHYD + 1; /* 9 */
const int SFUME = ANHYDRITE + 1;   /* 10 */
const int INERT = SFUME + 1;       /* 11 */
const int SLAG = INERT + 1;        /* 12 */
const int INERTAGG = SLAG + 1;     /* 13 */
const int ASG = INERTAGG + 1;      /* 14 */
const int CAS2 = ASG + 1;          /* 15 */
const int AMSIL = CAS2 + 1;        /* 16 */
const int FAC3A = AMSIL + 1;       /* 17 */

const int FLYASH = FAC3A + 1;  /* 18 */
const int CH = FLYASH + 1;     /* 19 */
const int CSH = CH + 1;        /* 20 */
const int C3AH6 = CSH + 1;     /* 21 */
const int ETTR = C3AH6 + 1;    /* 22 */
const int ETTRC4AF = ETTR + 1; /* 23 */

const int AFM = ETTRC4AF + 1;    /* 24 */
const int FH3 = AFM + 1;         /* 25 */
const int POZZCSH = FH3 + 1;     /* 26 */
const int SLAGCSH = POZZCSH + 1; /* 27 */
const int CACL2 = SLAGCSH + 1;   /* 28 */
/* Friedel's salt */
const int FRIEDEL = CACL2 + 1; /* 29 */

/* Stratlingite (C2ASH8) */
const int STRAT = FRIEDEL + 1; /* 30 */

/* Gypsum formed from hemihydrate and anhydrite */
const int GYPSUMS = STRAT + 1;  /* 31 */
const int ABSGYP = GYPSUMS + 1; /* 32 */

const int CACO3 = ABSGYP + 1; /* 33 */
const int AFMC = CACO3 + 1;   /* 34 */

/***
 *	Phases for chloride ingress model and
 *	sulfate attack model
 ***/
const int BRUCITE = AFMC + 1; /* 35 */
const int MS = BRUCITE + 1;   /* 36 */

/* Free lime */
const int FREELIME = MS + 1; /* 37 */

/* Orthorhombic C3A */
const int OC3A = FREELIME + 1; /* 38 */
const int NSPHASES = OC3A;     /* 38 */

/***
 *	Diffusing species
 ***/
const int DIFFCSH = NSPHASES + 1;   /* 39 */
const int DIFFCH = DIFFCSH + 1;     /* 40 */
const int DIFFGYP = DIFFCH + 1;     /* 41 */
const int DIFFC3A = DIFFGYP + 1;    /* 42 */
const int DIFFC4A = DIFFC3A + 1;    /* 43 */
const int DIFFFH3 = DIFFC4A + 1;    /* 44 */
const int DIFFETTR = DIFFFH3 + 1;   /* 45 */
const int DIFFCACO3 = DIFFETTR + 1; /* 46 */
const int DIFFAS = DIFFCACO3 + 1;   /* 47 */
const int DIFFANH = DIFFAS + 1;     /* 48 */
const int DIFFHEM = DIFFANH + 1;    /* 49 */
const int DIFFCAS2 = DIFFHEM + 1;   /* 50 */
const int DIFFCACL2 = DIFFCAS2 + 1; /* 51 */
const int DIFFSO4 = DIFFCACL2 + 1;  /* 52 */

const int NDIFFPHASES = DIFFSO4 + 1; /* 53 */

/***
 *	Special types of porosity
 ***/

const int DRIEDP = NDIFFPHASES; /* 53 */
const int EMPTYDP = DRIEDP + 1; /* 54 */

/***
 *	Empty porosity due to self dessication
 ***/
const int EMPTYP = EMPTYDP + 1; /* 55 */

/***
 *	Crack porosity, defined as the porosity
 *	created when the microstructure is cracked.
 *	Can be saturated or empty, depending on the
 *	application (24 May 2004)
 ***/
const int CRACKP = EMPTYP + 1; /* 56 */

/***
 *	Offset for highlighting potentially
 *	soluble surface pixels in disrealnew
 ***/

const int OFFSET = CRACKP + 1; /* 57 */

/***
 *	Total number of types of pixels, which
 *	INCLUDES diffusing species
 ***/
const int NDIFFUS = OFFSET;

const int SANDINCONCRETE = OFFSET + 3; /* 60 */
const int COARSEAGG01INCONCRETE = SANDINCONCRETE + 1;
const int COARSEAGG02INCONCRETE = COARSEAGG01INCONCRETE + 1;
const int FINEAGG01INCONCRETE = COARSEAGG02INCONCRETE + 1;
const int FINEAGG02INCONCRETE = FINEAGG01INCONCRETE + 1;

const int NPHASES = FINEAGG02INCONCRETE + 1;

const int SAT = 255;

const int R_BROWN = 162;
const int G_BROWN = 117;
const int B_BROWN = 95;

const int R_BLUE = 0;
const int G_BLUE = 0;
const int B_BLUE = SAT;

const int R_CFBLUE = 0;
const int G_CFBLUE = 128;
const int B_CFBLUE = SAT;

const int R_RED = SAT;
const int G_RED = 0;
const int B_RED = 0;

const int R_GREEN = 0;
const int G_GREEN = SAT;
const int B_GREEN = 0;

const int R_WHITE = SAT;
const int G_WHITE = SAT;
const int B_WHITE = SAT;

const int R_BLACK = 0;
const int G_BLACK = 0;
const int B_BLACK = 0;

const int R_AQUA = 0;
const int G_AQUA = SAT;
const int B_AQUA = SAT;

const int R_LTURQUOISE = 174;
const int G_LTURQUOISE = 237;
const int B_LTURQUOISE = 237;

const int R_YELLOW = SAT;
const int G_YELLOW = SAT;
const int B_YELLOW = 0;

const int R_LYELLOW = SAT;
const int G_LYELLOW = SAT;
const int B_LYELLOW = SAT / 2;

const int R_GOLD = SAT;
const int G_GOLD = 215;
const int B_GOLD = 0;

const int R_OLIVE = SAT / 2;
const int G_OLIVE = SAT / 2;
const int B_OLIVE = 0;

const int R_LOLIVE = 150;
const int G_LOLIVE = 150;
const int B_LOLIVE = 0;

const int R_DOLIVE = SAT / 4;
const int G_DOLIVE = SAT / 4;
const int B_DOLIVE = 0;

const int R_DBLUE = 0;
const int G_DBLUE = 0;
const int B_DBLUE = SAT / 2;

const int R_VIOLET = SAT / 2;
const int G_VIOLET = 0;
const int B_VIOLET = SAT;

const int R_LAVENDER = 230;
const int G_LAVENDER = 230;
const int B_LAVENDER = 250;

const int R_PLUM = 238;
const int G_PLUM = 174;
const int B_PLUM = 238;

const int R_FIREBRICK = 178;
const int G_FIREBRICK = 34;
const int B_FIREBRICK = 34;

const int R_MUTEDFIREBRICK = 178;
const int G_MUTEDFIREBRICK = 128;
const int B_MUTEDFIREBRICK = 128;

const int R_SEAGREEN = SAT / 2;
const int G_SEAGREEN = 250;
const int B_SEAGREEN = SAT / 2;

const int R_MAGENTA = SAT;
const int G_MAGENTA = 0;
const int B_MAGENTA = SAT;

const int R_ORANGE = SAT;
const int G_ORANGE = 165;
const int B_ORANGE = 0;

const int R_PEACH = SAT;
const int G_PEACH = 170;
const int B_PEACH = 128;

const int R_WHEAT = 245;
const int G_WHEAT = 222;
const int B_WHEAT = 179;

const int R_TAN = 210;
const int G_TAN = 180;
const int B_TAN = 140;

const int R_DGREEN = 0;
const int G_DGREEN = 100;
const int B_DGREEN = 0;

const int R_LGREEN = SAT / 2;
const int G_LGREEN = SAT;
const int B_LGREEN = SAT / 2;

const int R_LIME = 51;
const int G_LIME = 205;
const int B_LIME = 51;

const int R_DLIME = 26;
const int G_DLIME = 103;
const int B_DLIME = 26;

const int R_LLIME = 128;
const int G_LLIME = 255;
const int B_LLIME = 0;

const int R_GRAY = 178;
const int G_GRAY = 178;
const int B_GRAY = 178;

const int R_DGRAY = SAT / 4;
const int G_DGRAY = SAT / 4;
const int B_DGRAY = SAT / 4;

const int R_CHARCOAL = 50;
const int G_CHARCOAL = 50;
const int B_CHARCOAL = 50;

const int R_LGRAY = 3 * SAT / 4;
const int G_LGRAY = 3 * SAT / 4;
const int B_LGRAY = 3 * SAT / 4;

const int R_DAQUA = SAT / 4;
const int G_DAQUA = SAT / 2;
const int B_DAQUA = SAT / 2;

const int R_SALMON = SAT;
const int G_SALMON = SAT / 2;
const int B_SALMON = SAT / 2;

const int R_SKYBLUE = SAT / 4;
const int G_SKYBLUE = SAT / 2;
const int B_SKYBLUE = SAT;

const int R_DPINK = SAT;
const int G_DPINK = 0;
const int B_DPINK = SAT / 2;

const int R_PINK = SAT;
const int G_PINK = 105;
const int B_PINK = 180;

const int R_ORANGERED = SAT;
const int G_ORANGERED = 69;
const int B_ORANGERED = 0;

// Global variables

using namespace std;

vector<int> Mic;
string RootName;
string TypeName;
string OutName;
int Xsize, Ysize, Zsize;

// Function declarations
int checkargs(int argc, char **argv);
int processImageFiles(vector<string> &names, vector<string> &times);
int getFileNamesAndTimes(vector<string> &names, vector<string> &times);
int writeXYZFile(const string &times);
int countSolid(void);
bool isSolid(int i, int j, int k);
int toIndex(int i, int j, int k);
void getPcolors(vector<float> &red, vector<float> &green, vector<float> &blue);
void getVcolors(vector<float> &red, vector<float> &green, vector<float> &blue);
void getTcolors(vector<float> &red, vector<float> &green, vector<float> &blue);
void printHelp(void);
string getLeftPaddingString(string const &str, int n, char paddedChar);
