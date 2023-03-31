/***
@brief Header file for viz program

*/
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <getopt.h>

// Preprocessor defines
#define VOID   		    0
#define POROSITY 		1
#define C3S				2
#define C3A				3
#define ANHYDRITE		4
#define C2S			    5
#define C4AF			6
#define GYPSUM			7
#define SFUME			8
#define K2SO4           9
#define HEMIHYD			10
#define NA2SO4          11
#define CSH             12
#define CH              13
#define HEMIANH         14
#define CACO3           15
#define ETTR            16
#define MONOSULF        17
#define AFM             18
#define AFMC            19
#define HYDROTALC       20

#define NSPHASES		19

#define NPHASES		    21

#define SAT			255

#define R_BROWN		162
#define G_BROWN		117
#define B_BROWN		95

#define R_BLUE		0
#define G_BLUE		0
#define B_BLUE		SAT

#define R_CFBLUE	0
#define G_CFBLUE	128
#define B_CFBLUE	SAT

#define R_RED		SAT
#define G_RED		0
#define B_RED		0

#define R_GREEN		0
#define G_GREEN		SAT
#define B_GREEN		0

#define R_WHITE		SAT
#define G_WHITE		SAT
#define B_WHITE		SAT

#define R_BLACK		0
#define G_BLACK		0
#define B_BLACK		0

#define R_AQUA		0
#define G_AQUA		SAT
#define B_AQUA		SAT

#define R_LTURQUOISE	174
#define G_LTURQUOISE	237
#define B_LTURQUOISE	237

#define R_YELLOW	SAT
#define G_YELLOW	SAT
#define B_YELLOW	0

#define R_LYELLOW	SAT
#define G_LYELLOW	SAT
#define B_LYELLOW	SAT/2

#define R_GOLD		SAT
#define G_GOLD		215
#define B_GOLD		0

#define R_OLIVE		SAT/2
#define G_OLIVE		SAT/2
#define B_OLIVE		0

#define R_LOLIVE	150
#define G_LOLIVE	150
#define B_LOLIVE	0

#define R_DOLIVE	SAT/4
#define G_DOLIVE	SAT/4
#define B_DOLIVE	0

#define R_DBLUE		0
#define G_DBLUE		0
#define B_DBLUE		SAT/2

#define R_VIOLET	SAT/2
#define G_VIOLET	0
#define B_VIOLET	SAT

#define R_LAVENDER	230
#define G_LAVENDER	230
#define B_LAVENDER	250

#define R_PLUM		238
#define G_PLUM		174
#define B_PLUM		238

#define R_FIREBRICK		178
#define G_FIREBRICK		34
#define B_FIREBRICK		34

#define R_MUTEDFIREBRICK		178
#define G_MUTEDFIREBRICK		128
#define B_MUTEDFIREBRICK		128

#define R_SEAGREEN	SAT/2
#define G_SEAGREEN	250
#define B_SEAGREEN	SAT/2

#define R_MAGENTA	SAT
#define G_MAGENTA	0
#define B_MAGENTA	SAT

#define R_ORANGE	SAT
#define G_ORANGE	165
#define B_ORANGE	0

#define R_PEACH	SAT
#define G_PEACH	170
#define B_PEACH	128

#define R_WHEAT		245
#define G_WHEAT		222
#define B_WHEAT		179

#define R_TAN		210
#define G_TAN		180
#define B_TAN		140

#define R_DGREEN	0
#define G_DGREEN	100
#define B_DGREEN	0

#define R_LGREEN	SAT/2
#define G_LGREEN	SAT
#define B_LGREEN	SAT/2

#define R_LIME		51
#define G_LIME		205
#define B_LIME		51

#define R_DLIME		26
#define G_DLIME		103
#define B_DLIME		26

#define R_LLIME		128
#define G_LLIME		255
#define B_LLIME		0

#define R_LYELLOW	SAT
#define G_LYELLOW	SAT
#define B_LYELLOW	SAT/2

#define R_GRAY		178
#define G_GRAY		178
#define B_GRAY		178

#define R_DGRAY		SAT/4
#define G_DGRAY		SAT/4
#define B_DGRAY		SAT/4

#define R_CHARCOAL	50
#define G_CHARCOAL	50
#define B_CHARCOAL	50

#define R_LGRAY		3*SAT/4
#define G_LGRAY		3*SAT/4
#define B_LGRAY		3*SAT/4

#define R_DAQUA		SAT/4
#define G_DAQUA		SAT/2
#define B_DAQUA		SAT/2

#define R_SALMON	SAT
#define G_SALMON	SAT/2
#define B_SALMON	SAT/2

#define R_SKYBLUE	SAT/4
#define G_SKYBLUE	SAT/2
#define B_SKYBLUE	SAT

#define R_DPINK		SAT
#define G_DPINK		0
#define B_DPINK		SAT/2

#define R_PINK		SAT
#define G_PINK		105
#define B_PINK		180

#define R_ORANGERED		SAT
#define G_ORANGERED		69
#define B_ORANGERED		0

// Global variables

using namespace std;

vector<int> Mic;
string RootName;
string TypeName;
string OutName;
int Xsize,Ysize,Zsize;

// Function declarations
int checkargs(int argc, char **argv);
int processImageFiles(vector<string> &names, vector<string> &times);
int getFileNamesAndTimes(vector<string> &names, vector<string> &times);
int writeXYZFile(const string &times);
int countSolid(void);
bool isSolid(int i, int j, int k);
int toIndex(int i, int j, int k);
void getVcolors(vector<float> &red, vector<float> &green, vector<float> &blue);
void getTcolors(vector<float> &red, vector<float> &green, vector<float> &blue);
void printHelp(void);
string getLeftPaddingString(string const &str, int n, char paddedChar);
