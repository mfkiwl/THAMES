/***
 *	This is the public header file required for
 *	VCCTL programs.  Every VCCTL program that wishes
 *	to operate within the VCCTL system must include
 *	this header file
 *
 *	Programmer:  Jeffrey W. Bullard
 *				 NIST
 *				 100 Bureau Drive Stop 8615
 *				 Gaithersburg, MD  20899-8615
 *
 *				 Phone:	301.975.5725
 *				 Fax:	301.990.6892
 *				 bullard@nist.gov
 *
 *				 15 March 2004
 ***/
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/***
 *	Version number string and version number
 *	for identifying the version under which a
 *	particular file was created
 ***/
#define VERSIONSTRING "#THAMES:Version:"
#define VERSIONNUMBER "2.6"

#define MAXSTRING 500 /* maximum length of strings */

/*******************************************************
 * Variables related to system size and
 * resolution
 *******************************************************/

/* maximum system size in pixels per dimension */
#define MAXSIZE 400
#define DEFAULTSYSTEMSIZE 100

#define LOWRES 1.00
#define MEDLORES 0.75
#define MEDHIRES 0.50
#define HIGHRES 0.25

#define DEFAULTRESOLUTION LOWRES

#define XSIZESTRING "#THAMES:X_Size:"
#define YSIZESTRING "#THAMES:Y_Size:"
#define ZSIZESTRING "#THAMES:Z_Size:"
#define IMGRESSTRING "#THAMES:Image_Resolution:"

/***
 *	Pre-defined strings for info files
 ***/
#define ONEPIXBIASSTRING "One_pixel_bias:"
#define AGGTHICKSTRING "Aggregate_thickness:"
#define PFILESTRING "Particle_file:"
#define NUM1PIXSTRING "One_pixel_particles:"
#define VFALITESTRING "Vol_frac_ALITE:"
#define VFBELITESTRING "Vol_frac_BELITE:"
#define VFC3ASTRING "Vol_frac_C3A:"
#define VFC4AFSTRING "Vol_frac_C4AF:"
#define VFK2SO4STRING "Vol_frac_K2SO4:"
#define VFNA2SO4STRING "Vol_frac_NA2SO4:"
#define SFALITESTRING "Surf_frac_ALITE:"
#define SFBELITESTRING "Surf_frac_BELITE:"
#define SFC3ASTRING "Surf_frac_C3A:"
#define SFC4ASTRING "Surf_frac_C4AF:"
#define SFK2SO4STRING "Surf_frac_K2SO4:"
#define SFNA2SO4STRING "Surf_frac_NA2SO4:"
#define VFGYPSTRING "Vol_frac_Dihydrate:"
#define VFHEMSTRING "Vol_frac_Hemihydrate:"
#define VFANHSTRING "Vol_frac_Anhydrite:"

/******************************************************
 * Define phase identifier for all species
 ******************************************************/

/***
 *	These are phases present in unhydrated
 *	blended cements
 ***/

const int VOID_ID = 0;
const int ELECTROLYTE_ID = 1;
const int ALITE_ID = 2;
const int BELITE_ID = 3;
const int C3A_ID = 4;
const int C4AF_ID = 5;
const int K2SO4_ID = 6;
const int NA2SO4_ID = 7;
const int GYPSUM_ID = 8;
const int BASSANITE_ID = 9;
const int ANHYDRITE_ID = 10;
const int SFUME_ID = 11;
const int GGBS_ID = 12;
const int STEELSLAG_ID = 13;

const int C2AS_ID = 14;
const int CA2S_ID = 15;
const int CAS_ID = 16;
const int CAS2_ID = 17;
const int MULLITE_ID = 18;
const int K6A2S_ID = 19;

const int CLINKER_ID = 100;
const int FLYASH_ID = 200;

const int ClinkerMembers[6] = {ALITE_ID, BELITE_ID, C3A_ID,
                               C4AF_ID,  K2SO4_ID,  NA2SO4_ID};
const int FlyAshMembers[6] = {C2AS_ID, CA2S_ID,    CAS_ID,
                              CAS2_ID, MULLITE_ID, K6A2S_ID};

/*********************************************************
 * Defines RGB colors by name on a SAT scale
 *********************************************************/

#define SAT 255

#define R_BROWN 162
#define G_BROWN 117
#define B_BROWN 95

#define R_BLUE 0
#define G_BLUE 0
#define B_BLUE SAT

#define R_CFBLUE 0
#define G_CFBLUE 128
#define B_CFBLUE SAT

#define R_RED SAT
#define G_RED 0
#define B_RED 0

#define R_GREEN 0
#define G_GREEN SAT
#define B_GREEN 0

#define R_WHITE SAT
#define G_WHITE SAT
#define B_WHITE SAT

#define R_BLACK 0
#define G_BLACK 0
#define B_BLACK 0

#define R_AQUA 0
#define G_AQUA SAT
#define B_AQUA SAT

#define R_LTURQUOISE 174
#define G_LTURQUOISE 237
#define B_LTURQUOISE 237

#define R_YELLOW SAT
#define G_YELLOW SAT
#define B_YELLOW 0

#define R_LYELLOW SAT
#define G_LYELLOW SAT
#define B_LYELLOW SAT / 2

#define R_GOLD SAT
#define G_GOLD 215
#define B_GOLD 0

#define R_OLIVE SAT / 2
#define G_OLIVE SAT / 2
#define B_OLIVE 0

#define R_LOLIVE 150
#define G_LOLIVE 150
#define B_LOLIVE 0

#define R_DOLIVE SAT / 4
#define G_DOLIVE SAT / 4
#define B_DOLIVE 0

#define R_DBLUE 0
#define G_DBLUE 0
#define B_DBLUE SAT / 2

#define R_VIOLET SAT / 2
#define G_VIOLET 0
#define B_VIOLET SAT

#define R_LAVENDER 230
#define G_LAVENDER 230
#define B_LAVENDER 250

#define R_PLUM 238
#define G_PLUM 174
#define B_PLUM 238

#define R_FIREBRICK 178
#define G_FIREBRICK 34
#define B_FIREBRICK 34

#define R_MUTEDFIREBRICK 178
#define G_MUTEDFIREBRICK 128
#define B_MUTEDFIREBRICK 128

#define R_SEAGREEN SAT / 2
#define G_SEAGREEN 250
#define B_SEAGREEN SAT / 2

#define R_MAGENTA SAT
#define G_MAGENTA 0
#define B_MAGENTA SAT

#define R_ORANGE SAT
#define G_ORANGE 165
#define B_ORANGE 0

#define R_PEACH SAT
#define G_PEACH 170
#define B_PEACH 128

#define R_WHEAT 245
#define G_WHEAT 222
#define B_WHEAT 179

#define R_TAN 210
#define G_TAN 180
#define B_TAN 140

#define R_DGREEN 0
#define G_DGREEN 100
#define B_DGREEN 0

#define R_LGREEN SAT / 2
#define G_LGREEN SAT
#define B_LGREEN SAT / 2

#define R_LIME 51
#define G_LIME 205
#define B_LIME 51

#define R_DLIME 26
#define G_DLIME 103
#define B_DLIME 26

#define R_LLIME 128
#define G_LLIME 255
#define B_LLIME 0

#define R_LYELLOW SAT
#define G_LYELLOW SAT
#define B_LYELLOW SAT / 2

#define R_GRAY 178
#define G_GRAY 178
#define B_GRAY 178

#define R_DGRAY SAT / 4
#define G_DGRAY SAT / 4
#define B_DGRAY SAT / 4

#define R_CHARCOAL 50
#define G_CHARCOAL 50
#define B_CHARCOAL 50

#define R_LGRAY 3 * SAT / 4
#define G_LGRAY 3 * SAT / 4
#define B_LGRAY 3 * SAT / 4

#define R_DAQUA SAT / 4
#define G_DAQUA SAT / 2
#define B_DAQUA SAT / 2

#define R_SALMON SAT
#define G_SALMON SAT / 2
#define B_SALMON SAT / 2

#define R_SKYBLUE SAT / 4
#define G_SKYBLUE SAT / 2
#define B_SKYBLUE SAT

#define R_DPINK SAT
#define G_DPINK 0
#define B_DPINK SAT / 2

#define R_PINK SAT
#define G_PINK 105
#define B_PINK 180

#define R_ORANGERED SAT
#define G_ORANGERED 69
#define B_ORANGERED 0

#define HTMLCODE 0
#define PLAINCODE 1

/* typedefs */

/* A coloured pixel. */

typedef struct {
  uint8_t red;
  uint8_t green;
  uint8_t blue;
} pixel_t;

/* Support for image rendering */

/* A picture. */

typedef struct {
  pixel_t *pixels;
  size_t width;
  size_t height;
} bitmap_t;

/* Different types of arrays */

typedef struct {
  int xsize;
  int *val;
} Int1d;

typedef struct {
  int xsize;
  short int *val;
} ShortInt1d;

typedef struct {
  int xsize;
  long int *val;
} LongInt1d;

typedef struct {
  int xsize;
  float *val;
} Float1d;

typedef struct {
  int xsize;
  float *val;
} Double1d;

typedef struct {
  int xsize;
  pixel_t *val;
} Pixel1d;

typedef struct {
  int xsize;
  int ysize;
  int *val;
} Int2d;

typedef struct {
  int xsize;
  int ysize;
  short int *val;
} ShortInt2d;

typedef struct {
  int xsize;
  int ysize;
  long int *val;
} LongInt2d;

typedef struct {
  int xsize;
  int ysize;
  float *val;
} Float2d;

typedef struct {
  int xsize;
  int ysize;
  double *val;
} Double2d;

typedef struct {
  int xsize;
  int ysize;
  char *val;
} Char2d;

typedef struct {
  size_t x;
  size_t y;
  size_t z;
  int *val;
} Int3d;

typedef struct {
  int xsize;
  int ysize;
  int zsize;
  short int *val;
} ShortInt3d;

typedef struct {
  int xsize;
  int ysize;
  int zsize;
  long int *val;
} LongInt3d;

typedef struct {
  int xsize;
  int ysize;
  int zsize;
  float *val;
} Float3d;

typedef struct {
  int xsize;
  int ysize;
  int zsize;
  double *val;
} Double3d;

typedef struct {
  int xsize;
  int ysize;
  int zsize;
  char *val;
} Char3d;

/***
 *	Function declarations needed by createmic
 ***/
#include <auxmiccomplex.h>
#include <auxmicfuncs.h>
