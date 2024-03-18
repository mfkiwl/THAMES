/**
@file viz.cc

*/

#include "viz.h"

/* Make Transparency a global variable */
float transparency;

int main(int argc, char *argv[]) {

  vector<string> names;
  vector<string> times;

  if (checkargs(argc, argv)) {
    exit(1);
  }

  if (processImageFiles(names, times)) {
    exit(1);
  }

  cout << endl << "All done!" << endl << endl;
  exit(0);
}

int checkargs(int argc, char **argv) {
  bool wellformed = false;
  char *iroot, *itype, *tstring;

  if (argc < 3) {
    wellformed = false;
  }

  // Many of the variables here are defined in the getopts.h system header file
  // Can define more options here if we want

  static struct option long_opts[] = {
      {"input", required_argument, 0, 'i'},
      {"format", required_argument, 0, 'f'},
      {"transparency", required_argument, 0, 't'},
      {"help", no_argument, 0, 'h'},
      {NULL, 0, 0, 0}};

  int opt_char;
  int option_index;

  TypeName = "vcctl";

  while ((opt_char = getopt_long(argc, argv, "i:f:t:h", long_opts,
                                 &option_index)) != -1) {
    switch (opt_char) {
    // -i or --input
    case static_cast<int>('i'):
      wellformed = true;
      iroot = optarg;
      RootName = iroot;
      cout << "RootName = " << RootName << endl;
      cout.flush();
      break;
    // -f or --format
    case static_cast<int>('f'):
      itype = optarg;
      TypeName = itype;
      for (int j = 0; j < TypeName.size(); ++j) {
        TypeName[j] = tolower(TypeName[j]);
      }
      cout << "TypeName = " << TypeName << endl;
      cout.flush();
      break;
    // -t or --transparency
    case static_cast<int>('t'):
      tstring = optarg;
      transparency = atof(tstring);
      cout << "Transparency = " << transparency << endl;
      cout.flush();
      break;
    // -h or --help
    case static_cast<int>('h'):
    // Unrecognized option
    case static_cast<int>('?'):
    default:
      wellformed = false;
      break;
    }
  }

  if (!wellformed) {
    printHelp();
    return (1);
  }

  if ((TypeName != "vcctl") && (TypeName != "thames")) {
    printHelp();
    return (1);
  }

  return (0);
}

int processImageFiles(vector<string> &names, vector<string> &times) {
  int valin, iz, iy, ix;
  string buff, buff1;

  if (getFileNamesAndTimes(names, times)) {
    return (1);
  }

  /****
   *    NOTE:  MUST READ
   *            (1) SOFTWARE VERSION OF INPUT FILE (Version)
   *            (1) SYSTEM SIZE (Xsize,Ysize,Zsize)
   *            (2) SYSTEM RESOLUTION (Res)
   *
   *    Then set global variables Syspix, Sizemag, Isizemag
   *    rather than hardwiring them as preprocessor
   *    defines
   ****/

  for (int i = 0; i < names.size(); ++i) {
    ifstream in(names[i].c_str());

    // Read the version stuff
    in >> buff;
    in >> buff1;

    // Read the x,y,z dimensions
    in >> buff;
    in >> Xsize;
    in >> buff;
    in >> Ysize;
    in >> buff;
    in >> Zsize;

    // Read the resolution
    in >> buff >> buff1;

    cout << "Xsize = " << Xsize << "; Ysize = " << Ysize
         << "; Zsize = " << Zsize << endl;

    // Create the Mic vector
    Mic.clear();

    for (iz = 0; iz < Zsize; ++iz) {
      for (iy = 0; iy < Ysize; ++iy) {
        for (ix = 0; ix < Xsize; ++ix) {

          in >> valin;
          Mic.push_back(valin);
        }
      }
    }

    in.close();

    // Now create the xyz file

    cout << "Processing file " << i + 1 << " of " << times.size() << endl;
    cout.flush();
    if (writeXYZFile(times[i])) {
      cout << endl << "ERROR: Problem writing xyz file " << i << endl << endl;
      return (1);
    }
  }

  return (0);
}

/***
 *  Next function uses the .xyz format description
 *
 ***/
int writeXYZFile(const string &times) {

  float ftime = stof(times);

  vector<float> red((int)(NPHASES), 0.0);
  vector<float> blue((int)(NPHASES), 0.0);
  vector<float> green((int)(NPHASES), 0.0);

  if (TypeName == "thames") {
    getTcolors(red, green, blue);
  } else {
    getVcolors(red, green, blue);
  }

  string outname = RootName;
  outname.append("_");
  outname.append(times);
  outname.append(".xyz");

  ofstream out(outname.c_str());

  /*
  int numvox = countSolid();
  */

  int numvox = Xsize * Ysize * Zsize;

  // Write the file headers
  out << numvox << endl; // Number of voxels to visualize
  out << "Time:\t" << ftime << endl;

  // Loop over all voxels and write out the solid ones

  float x, y, z;
  int idx;

  for (int k = 0; k < Zsize; ++k) {
    z = (float)k;
    for (int j = 0; j < Ysize; ++j) {
      y = (float)j;
      for (int i = 0; i < Xsize; ++i) {
        x = (float)i;
        idx = Mic[toIndex(i, j, k)];
        if (isSolid(i, j, k)) {
          /*
          cout << "Phase " << idx << ", " << red[idx] << "\t" << green[idx]
               << "\t" << blue[idx] << "\t" << x << "\t" << y << "\t" << z
               << "\t0.0" << endl;
          */
          out << red[idx] << "\t" << green[idx] << "\t" << blue[idx] << "\t"
              << x << "\t" << y << "\t" << z << "\t0.0" << endl;
        } else {
          /*
          cout << "Phase " << idx << ", " << red[idx] << "\t" << green[idx]
               << "\t" << blue[idx] << "\t" << x << "\t" << y << "\t" << z
               << "\t" << transparency << endl;
          */
          out << red[idx] << "\t" << green[idx] << "\t" << blue[idx] << "\t"
              << x << "\t" << y << "\t" << z << "\t" << transparency << endl;
        }
      }
    }
  }

  out.close();

  return (0);
}

int countSolid(void) {
  int numsol = 0;

  for (int k = 0; k < Zsize; ++k) {
    for (int j = 0; j < Ysize; ++j) {
      for (int i = 0; i < Xsize; ++i) {
        if (isSolid(i, j, k))
          numsol++;
      }
    }
  }
  return (numsol);
}

bool isSolid(int i, int j, int k) {
  bool result = false;

  int id = Mic[toIndex(i, j, k)];

  //  cout << "Mic(" << i << "," << j << "," << k << ") = " << id;
  if (TypeName == "thames") {
    if (id > 1)
      result = true;
  } else {
    if (id != (int)ELECTROLYTE_ID && id <= (int)NSPHASES) {
      result = true;
    }
  }
  //  cout << "; result = " << result << endl;

  return (result);
}

int toIndex(int i, int j, int k) {
  int pixnum = (k * Xsize * Ysize) + (j * Xsize) + i;
  return (pixnum);
}

void getVcolors(vector<float> &red, vector<float> &green, vector<float> &blue) {
  for (int i = 0; i < NPHASES; ++i) {
    switch (i) {
    case ELECTROLYTE_ID:
      red[i] = (float)(R_CHARCOAL) / (float)(SAT);
      green[i] = (float)(G_CHARCOAL) / (float)(SAT);
      blue[i] = (float)(B_CHARCOAL) / (float)(SAT);
      break;
    case EMPTYP:
      red[i] = (float)(R_BLACK) / (float)(SAT);
      green[i] = (float)(G_BLACK) / (float)(SAT);
      blue[i] = (float)(B_BLACK) / (float)(SAT);
      break;
    case EMPTYDP:
      red[i] = (float)(R_BLACK) / (float)(SAT);
      green[i] = (float)(G_BLACK) / (float)(SAT);
      blue[i] = (float)(B_BLACK) / (float)(SAT);
      break;
    case DRIEDP:
      red[i] = (float)(R_BLACK) / (float)(SAT);
      green[i] = (float)(G_BLACK) / (float)(SAT);
      blue[i] = (float)(B_BLACK) / (float)(SAT);
      break;
    case CH:
      red[i] = (float)(R_BLUE) / (float)(SAT);
      green[i] = (float)(G_BLUE) / (float)(SAT);
      blue[i] = (float)(B_BLUE) / (float)(SAT);
      break;
    case CSH:
      red[i] = (float)(R_WHEAT) / (float)(SAT);
      green[i] = (float)(G_WHEAT) / (float)(SAT);
      blue[i] = (float)(B_WHEAT) / (float)(SAT);
      break;
    case C2S:
      red[i] = (float)(R_BROWN) / (float)(SAT);
      green[i] = (float)(G_BROWN) / (float)(SAT);
      blue[i] = (float)(B_BROWN) / (float)(SAT);
      break;
    case C3S:
      red[i] = (float)(R_CFBLUE) / (float)(SAT);
      green[i] = (float)(G_CFBLUE) / (float)(SAT);
      blue[i] = (float)(B_CFBLUE) / (float)(SAT);
      break;
    case C3A:
      red[i] = (float)(R_GRAY) / (float)(SAT);
      green[i] = (float)(G_GRAY) / (float)(SAT);
      blue[i] = (float)(B_GRAY) / (float)(SAT);
      break;
    case C4AF:
      red[i] = (float)(R_WHITE) / (float)(SAT);
      green[i] = (float)(G_WHITE) / (float)(SAT);
      blue[i] = (float)(B_WHITE) / (float)(SAT);
      break;
    case K2SO4:
      red[i] = (float)(R_RED) / (float)(SAT);
      green[i] = (float)(G_RED) / (float)(SAT);
      blue[i] = (float)(B_RED) / (float)(SAT);
      break;
    case NA2SO4:
      red[i] = (float)(R_SALMON) / (float)(SAT);
      green[i] = (float)(G_SALMON) / (float)(SAT);
      blue[i] = (float)(B_SALMON) / (float)(SAT);
      break;
    case GYPSUM:
      red[i] = (float)(R_YELLOW) / (float)(SAT);
      green[i] = (float)(G_YELLOW) / (float)(SAT);
      blue[i] = (float)(B_YELLOW) / (float)(SAT);
      break;
    case ABSGYP:
      red[i] = (float)(R_YELLOW) / (float)(SAT);
      green[i] = (float)(G_YELLOW) / (float)(SAT);
      blue[i] = (float)(B_YELLOW) / (float)(SAT);
      break;
    case GYPSUMS:
      red[i] = (float)(R_YELLOW) / (float)(SAT);
      green[i] = (float)(G_YELLOW) / (float)(SAT);
      blue[i] = (float)(B_YELLOW) / (float)(SAT);
      break;
    case HEMIHYD:
      red[i] = (float)(R_LYELLOW) / (float)(SAT);
      green[i] = (float)(G_LYELLOW) / (float)(SAT);
      blue[i] = (float)(B_LYELLOW) / (float)(SAT);
      break;
    case ANHYDRITE:
      red[i] = (float)(R_GOLD) / (float)(SAT);
      green[i] = (float)(G_GOLD) / (float)(SAT);
      blue[i] = (float)(B_GOLD) / (float)(SAT);
      break;
    case SFUME:
      red[i] = (float)(R_DAQUA) / (float)(SAT);
      green[i] = (float)(G_DAQUA) / (float)(SAT);
      blue[i] = (float)(B_DAQUA) / (float)(SAT);
      break;
    case AMSIL:
      red[i] = (float)(R_AQUA) / (float)(SAT);
      green[i] = (float)(G_AQUA) / (float)(SAT);
      blue[i] = (float)(B_AQUA) / (float)(SAT);
      break;
    case INERT:
      red[i] = (float)(R_PLUM) / (float)(SAT);
      green[i] = (float)(G_PLUM) / (float)(SAT);
      blue[i] = (float)(B_PLUM) / (float)(SAT);
      break;
    case ETTR:
      red[i] = (float)(R_LOLIVE) / (float)(SAT);
      green[i] = (float)(G_LOLIVE) / (float)(SAT);
      blue[i] = (float)(B_LOLIVE) / (float)(SAT);
      break;
    case ETTRC4AF:
      red[i] = (float)(R_LOLIVE) / (float)(SAT);
      green[i] = (float)(G_LOLIVE) / (float)(SAT);
      blue[i] = (float)(B_LOLIVE) / (float)(SAT);
      break;
    case AFM:
      red[i] = (float)(R_OLIVE) / (float)(SAT);
      green[i] = (float)(G_OLIVE) / (float)(SAT);
      blue[i] = (float)(B_OLIVE) / (float)(SAT);
      break;
    case AFMC:
      red[i] = (float)(R_OLIVE) / (float)(SAT);
      green[i] = (float)(G_OLIVE) / (float)(SAT);
      blue[i] = (float)(B_OLIVE) / (float)(SAT);
      break;
    case STRAT:
      red[i] = (float)(R_DOLIVE) / (float)(SAT);
      green[i] = (float)(G_DOLIVE) / (float)(SAT);
      blue[i] = (float)(B_DOLIVE) / (float)(SAT);
      break;
    case CACL2:
      red[i] = (float)(R_PEACH) / (float)(SAT);
      green[i] = (float)(G_PEACH) / (float)(SAT);
      blue[i] = (float)(B_PEACH) / (float)(SAT);
      break;
    case FRIEDEL:
      red[i] = (float)(R_MAGENTA) / (float)(SAT);
      green[i] = (float)(G_MAGENTA) / (float)(SAT);
      blue[i] = (float)(B_MAGENTA) / (float)(SAT);
      break;
    case FH3:
      red[i] = (float)(R_DAQUA) / (float)(SAT);
      green[i] = (float)(G_DAQUA) / (float)(SAT);
      blue[i] = (float)(B_DAQUA) / (float)(SAT);
      break;
    case POZZCSH:
      red[i] = (float)(R_LTURQUOISE) / (float)(SAT);
      green[i] = (float)(G_LTURQUOISE) / (float)(SAT);
      blue[i] = (float)(B_LTURQUOISE) / (float)(SAT);
      break;
    case INERTAGG:
      red[i] = (float)(R_FIREBRICK) / (float)(SAT);
      green[i] = (float)(G_FIREBRICK) / (float)(SAT);
      blue[i] = (float)(B_FIREBRICK) / (float)(SAT);
      break;
    case SANDINCONCRETE:
      red[i] = (float)(R_MUTEDFIREBRICK) / (float)(SAT);
      green[i] = (float)(G_MUTEDFIREBRICK) / (float)(SAT);
      blue[i] = (float)(B_MUTEDFIREBRICK) / (float)(SAT);
      break;
    case COARSEAGG01INCONCRETE:
      red[i] = (float)(R_FIREBRICK) / (float)(SAT);
      green[i] = (float)(G_FIREBRICK) / (float)(SAT);
      blue[i] = (float)(B_FIREBRICK) / (float)(SAT);
      break;
    case COARSEAGG02INCONCRETE:
      red[i] = (float)(R_MAGENTA) / (float)(SAT);
      green[i] = (float)(G_MAGENTA) / (float)(SAT);
      blue[i] = (float)(B_MAGENTA) / (float)(SAT);
      break;
    case FINEAGG01INCONCRETE:
      red[i] = (float)(R_FIREBRICK) / (float)(SAT);
      green[i] = (float)(G_FIREBRICK) / (float)(SAT);
      blue[i] = (float)(B_FIREBRICK) / (float)(SAT);
      break;
    case FINEAGG02INCONCRETE:
      red[i] = (float)(R_MAGENTA) / (float)(SAT);
      green[i] = (float)(G_MAGENTA) / (float)(SAT);
      blue[i] = (float)(B_MAGENTA) / (float)(SAT);
      break;
    case CACO3:
      red[i] = (float)(R_LIME) / (float)(SAT);
      green[i] = (float)(G_LIME) / (float)(SAT);
      blue[i] = (float)(B_LIME) / (float)(SAT);
      break;
    case FREELIME:
      red[i] = (float)(R_LLIME) / (float)(SAT);
      green[i] = (float)(G_LLIME) / (float)(SAT);
      blue[i] = (float)(B_LLIME) / (float)(SAT);
      break;
    case FLYASH:
      red[i] = (float)(R_DGRAY) / (float)(SAT);
      green[i] = (float)(G_DGRAY) / (float)(SAT);
      blue[i] = (float)(B_DGRAY) / (float)(SAT);
      break;
    case FAC3A:
      red[i] = (float)(R_GRAY) / (float)(SAT);
      green[i] = (float)(G_GRAY) / (float)(SAT);
      blue[i] = (float)(B_GRAY) / (float)(SAT);
      break;
    case ASG:
      red[i] = (float)(R_ORANGE) / (float)(SAT);
      green[i] = (float)(G_ORANGE) / (float)(SAT);
      blue[i] = (float)(B_ORANGE) / (float)(SAT);
      break;
    case SLAGCSH:
      red[i] = (float)(R_SEAGREEN) / (float)(SAT);
      green[i] = (float)(G_SEAGREEN) / (float)(SAT);
      blue[i] = (float)(B_SEAGREEN) / (float)(SAT);
      break;
    case SLAG:
      red[i] = (float)(R_DGREEN) / (float)(SAT);
      green[i] = (float)(G_DGREEN) / (float)(SAT);
      blue[i] = (float)(B_DGREEN) / (float)(SAT);
      break;
    case CAS2:
      red[i] = (float)(R_DBLUE) / (float)(SAT);
      green[i] = (float)(G_DBLUE) / (float)(SAT);
      blue[i] = (float)(B_DBLUE) / (float)(SAT);
      break;
    case BRUCITE:
      red[i] = (float)(R_DLIME) / (float)(SAT);
      green[i] = (float)(G_DLIME) / (float)(SAT);
      blue[i] = (float)(B_DLIME) / (float)(SAT);
      break;
    case MS:
      red[i] = (float)(R_ORANGERED) / (float)(SAT);
      green[i] = (float)(G_ORANGERED) / (float)(SAT);
      blue[i] = (float)(B_ORANGERED) / (float)(SAT);
      break;
    default:
      red[i] = (float)(R_LAVENDER) / (float)(SAT);
      green[i] = (float)(G_LAVENDER) / (float)(SAT);
      blue[i] = (float)(B_LAVENDER) / (float)(SAT);
    }
  }

  return;
}

void getTcolors(vector<float> &red, vector<float> &green, vector<float> &blue) {
  for (int i = 0; i < NPHASES; ++i) {
    switch (i) {
    case 0: // THAMES Void
      red[i] = (float)(R_BLACK) / (float)(SAT);
      green[i] = (float)(R_BLACK) / (float)(SAT);
      blue[i] = (float)(R_BLACK) / (float)(SAT);
      break;
    case 1: // THAMES saturated porosity
      red[i] = (float)(R_CHARCOAL) / (float)(SAT);
      green[i] = (float)(R_CHARCOAL) / (float)(SAT);
      blue[i] = (float)(R_CHARCOAL) / (float)(SAT);
      break;
    case 14: // THAMES portlandite
      red[i] = (float)(R_BLUE) / (float)(SAT);
      green[i] = (float)(G_BLUE) / (float)(SAT);
      blue[i] = (float)(B_BLUE) / (float)(SAT);
      break;
    case 11: // THAMES CSH
      red[i] = (float)(R_WHEAT) / (float)(SAT);
      green[i] = (float)(G_WHEAT) / (float)(SAT);
      blue[i] = (float)(B_WHEAT) / (float)(SAT);
      break;
    case 2: // THAMES alite
      red[i] = (float)(R_BROWN) / (float)(SAT);
      green[i] = (float)(G_BROWN) / (float)(SAT);
      blue[i] = (float)(B_BROWN) / (float)(SAT);
      break;
    case 5: // THAMES belite
      red[i] = (float)(R_CFBLUE) / (float)(SAT);
      green[i] = (float)(G_CFBLUE) / (float)(SAT);
      blue[i] = (float)(B_CFBLUE) / (float)(SAT);
      break;
    case 12: // THAMES calcite
      red[i] = (float)(R_CFBLUE) / (float)(SAT);
      green[i] = (float)(G_CFBLUE) / (float)(SAT);
      blue[i] = (float)(B_CFBLUE) / (float)(SAT);
      break;
    case 3: // THAMES C3A
      red[i] = (float)(R_GRAY) / (float)(SAT);
      green[i] = (float)(G_GRAY) / (float)(SAT);
      blue[i] = (float)(B_GRAY) / (float)(SAT);
      break;
    case 6: // THAMES C4AF
      red[i] = (float)(R_WHITE) / (float)(SAT);
      green[i] = (float)(G_WHITE) / (float)(SAT);
      blue[i] = (float)(B_WHITE) / (float)(SAT);
      break;
      //          case 8:  // THAMES arcanite
      //              red[i] = (float)(R_RED)/(float)(SAT);
      //              green[i] = (float)(G_RED)/(float)(SAT);
      //              blue[i] = (float)(B_RED)/(float)(SAT);
      //              break;
    case 10: // THAMES thenardite
      red[i] = (float)(R_SALMON) / (float)(SAT);
      green[i] = (float)(G_SALMON) / (float)(SAT);
      blue[i] = (float)(B_SALMON) / (float)(SAT);
      break;
    case 7: // THAMES gypsum
      red[i] = (float)(R_YELLOW) / (float)(SAT);
      green[i] = (float)(G_YELLOW) / (float)(SAT);
      blue[i] = (float)(B_YELLOW) / (float)(SAT);
      break;
    case 9: // THAMES hemihydrate
      red[i] = (float)(R_LYELLOW) / (float)(SAT);
      green[i] = (float)(G_LYELLOW) / (float)(SAT);
      blue[i] = (float)(B_LYELLOW) / (float)(SAT);
      break;
      //          case 14: // THAMES hemianh
      //              red[i] = (float)(R_LYELLOW)/(float)(SAT);
      //              green[i] = (float)(G_LYELLOW)/(float)(SAT);
      //              blue[i] = (float)(B_LYELLOW)/(float)(SAT);
      //              break;
    case 4: // THAMES anhydrite.
      red[i] = (float)(R_GOLD) / (float)(SAT);
      green[i] = (float)(G_GOLD) / (float)(SAT);
      blue[i] = (float)(B_GOLD) / (float)(SAT);
      break;
    case 8: // THAMES silica fume
      red[i] = (float)(R_AQUA) / (float)(SAT);
      green[i] = (float)(G_AQUA) / (float)(SAT);
      blue[i] = (float)(B_AQUA) / (float)(SAT);
      break;
    case 15: // THAMES AFT
      red[i] = (float)(127.0) / (float)(SAT);
      green[i] = (float)(0.0) / (float)(SAT);
      blue[i] = (float)(255.0) / (float)(SAT);
      break;
    case 16: // THAMES MONOSULF
      red[i] = (float)(238.0) / (float)(SAT);
      green[i] = (float)(174.0) / (float)(SAT);
      blue[i] = (float)(238.0) / (float)(SAT);
      break;
    case 19: // THAMES AFMC
      red[i] = (float)(R_OLIVE) / (float)(SAT);
      green[i] = (float)(G_OLIVE) / (float)(SAT);
      blue[i] = (float)(B_OLIVE) / (float)(SAT);
      break;
    case 18: // THAMES AFM
      red[i] = (float)(R_DOLIVE) / (float)(SAT);
      green[i] = (float)(G_DOLIVE) / (float)(SAT);
      blue[i] = (float)(B_DOLIVE) / (float)(SAT);
      break;
      //          case 19:
      //              red[i] = (float)(26.0)/(float)(SAT);
      //              green[i] = (float)(103.0)/(float)(SAT);
      //              blue[i] = (float)(26.0)/(float)(SAT);
      //              break;
    case 13: // THAMES HYDROTALC
      red[i] = (float)(200.0) / (float)(SAT);
      green[i] = (float)(200.0) / (float)(SAT);
      blue[i] = (float)(208.0) / (float)(SAT);
      break;
    default:
      red[i] = (float)(R_LAVENDER) / (float)(SAT);
      green[i] = (float)(G_LAVENDER) / (float)(SAT);
      blue[i] = (float)(B_LAVENDER) / (float)(SAT);
    }
  }

  return;
}

void printHelp(void) {
  cout << "viz converts one or more VCCTL images to xyz format" << endl;
  cout << "If  more that one file matches the supplied root name," << endl;
  cout << "they will all be processed as time frames, producing" << endl;
  cout << "one xyz file for each img file." << endl << endl;
  cout << "Usage: viz -i,--input rootname -f format -t,--transparency value"
       << endl;
  cout << "       rootname is everything before .img in the file name" << endl
       << endl;
  return;
}

int getFileNamesAndTimes(vector<string> &names, vector<string> &times) {
  string buff, subword;
  string stringtofind = "img";
  int equalOrNot;

  string cmd = "ls ";
  cmd += RootName;
  cmd += ".img* > scratchfilelist.txt";

  cout << "\t" << cmd << endl;
  system(cmd.c_str());

  ifstream in("scratchfilelist.txt");
  string delimdot = ".";
  string delimh = "h";

  names.clear();
  times.clear();

  size_t pos = 0;
  size_t pos1 = 0;
  bool done = false;
  while (!in.eof()) {
    in >> buff;
    if (!in.eof()) {
      cout << buff << endl;
      pos = 0;
      names.push_back(buff);
      times.push_back("0.0");
      done = false;
      while (((pos = buff.find(delimdot)) != string::npos) && !done) {
        subword = buff.substr(0, pos);
        cout << "subword = " << subword << endl;
        cout << "buff pre-erase = " << buff << endl;
        buff.erase(0, pos + delimdot.length());
        cout << "buff post-erase = " << buff << endl;
        equalOrNot = subword.compare(stringtofind);
        if (equalOrNot == 0)
          done = true;
      }
      if (done) {
        done = true;
        pos1 = buff.find(delimh);
        if (pos1 != string::npos) {
          times[times.size() - 1] = buff.substr(0, pos1);
        }
        cout << times[times.size() - 1];
      }
    }
  }

  in.close();

  // Sort both the names and the times in order of time
  // Will need to convert the time strings to floats

  vector<float> timef;
  timef.resize(times.size(), 0.0);
  string tmptime, tmpname;
  float tmpftime;

  for (int i = 0; i < names.size(); ++i) {
    cout << "Presort Name " << i << " = " << names[i] << endl;
  }
  for (int i = 0; i < times.size(); ++i) {
    cout << "Presort Time " << i << " = " << times[i] << endl;
  }

  for (int i = 0; i < times.size(); ++i) {
    timef[i] = stof(times[i]);
  }

  for (int i = 0; i < timef.size() - 1; ++i) {
    for (int j = i; j < timef.size(); ++j) {
      if (timef[j] < timef[i]) {
        tmpftime = timef[i];
        tmptime = times[i];
        tmpname = names[i];
        timef[i] = timef[j];
        times[i] = times[j];
        names[i] = names[j];
        timef[j] = tmpftime;
        times[j] = tmptime;
        names[j] = tmpname;
      }
    }
  }

  // Pad the time strings with zeros to the left
  // for file sorting later

  string padtime;
  for (int i = 0; i < times.size(); ++i) {
    padtime = getLeftPaddingString(times[i], 7, '0');
    times[i] = padtime;
  }

  for (int i = 0; i < names.size(); ++i) {
    cout << "Name " << i << " = " << names[i] << endl;
  }
  for (int i = 0; i < times.size(); ++i) {
    cout << "Time " << i << " = " << times[i] << endl;
  }

  return (0);
}

string getLeftPaddingString(string const &str, int n, char paddedChar = ' ') {
  ostringstream ss;
  ss << right << setfill(paddedChar) << setw(n) << str;
  return ss.str();
}
