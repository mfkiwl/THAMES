#include "vcctl2thames.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  int xsize, ysize, zsize;
  string input_filename;
  float res;

  vector<int> mic;
  map<int, int> idmap;

  checkargs(argc, argv);

  printBanner();

  // The next function reads the VCCTL image
  // and asks the user to make the necessary
  // correspondences between VCCTL and THAMES

  if (readVCCTLImage(input_filename, xsize, ysize, zsize, res, mic, idmap)) {
    exit(1);
  }

  // The next function creates the output file

  if (writeTHAMESImage(input_filename, xsize, ysize, zsize, res, mic, idmap)) {
    exit(1);
  }

  exit(0);
}

/* Function definitions */

void printBanner(void) {
  cout << endl;
  cout << "***********************************" << endl;
  cout << "*  Welcome to vcctl2thames        *" << endl;
  cout << "***********************************" << endl;
  cout << endl;
  return;
}

void checkargs(int argc, char **argv) {

  // Many of the variables here are defined in the getopts.h system header file
  // Can define more options here if we want
  const char *const short_opts = "vdh";
  const option long_opts[] = {{"verbose", no_argument, nullptr, 'v'},
                              {"debug", no_argument, nullptr, 'd'},
                              {"help", no_argument, nullptr, 'h'},
                              {nullptr, no_argument, nullptr, 0}};

  Verbose = false;

  while (true) {

    const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

    if (-1 == opt)
      break; // breaks out of this while loop

    switch (opt) {
    case 'v':
      Verbose = true; // Verbose defined in vcctl2thames.h
      cout << "**Will produce verbose output**" << endl;
      break;
    case 'd':
      Debug = true; // Debug defined in vcctl2thames.h
      cout << "**Will produce debugging output**" << endl;
      break;
    case 'h': // -h or --help
    case '?': // Unrecognized option
    default:
      printHelp();
      exit(0);
      break;
    }
  }
}

void printHelp(void) {
  cout << endl;
  cout << "Usage: \"vcctl2thames [--verbose|-v] [--debug|-d] [--help|-h]\""
       << endl;
  cout << "        --verbose [-v]      Produce verbose output" << endl;
  cout << "        --debug [-d]        Produce debugging output" << endl;
  cout << "        --help [-h]         Print this help message" << endl;
  cout << endl;

  return;
}

int readVCCTLImage(string &input_filename, int &xsize, int &ysize, int &zsize,
                   float &res, vector<int> &mic, map<int, int> &idmap) {
  ifstream fin;
  bool found = false;
  string sval;
  int ival;

  openFile(input_filename, fin);

  if (readVCCTLHeader(fin, xsize, ysize, zsize, res)) {
    return 1;
  }

  // Now we read the voxel ids and store the unique ones
  // in the vector called vcctlid

  int numsize = (xsize * ysize * zsize);
  if (Verbose) {
    cout << "Verbose: numsize = " << numsize << endl;
  }

  mic.clear();
  vector<int> vcctlid;
  vcctlid.clear();

  for (int i = 0; i < numsize; ++i) {

    if (!fin) {
      cout << "ERROR:  File is not open!" << endl;
    }

    fin >> sval;
    if (Verbose)
      cout << "Verbose: Read VCCTL value " << sval << endl;

    ival = stoi(sval);
    mic.push_back(ival);

    // Search to see if this id has already been found

    found = false;
    for (int j = 0; (j < vcctlid.size()) && (!found); ++j) {
      if (ival == vcctlid.at(j))
        found = true;
    }
    if (!found) {
      vcctlid.push_back(ival);
      if (Verbose) {
        cout << "Verbose: Found VCCTL " << ival
             << ". Recognized values so far:" << endl;
        cout << "    ";
        for (int jj = 0; jj < vcctlid.size(); ++jj) {
          cout << vcctlid.at(jj) << " ";
        }
        cout << endl;
      }
    }
  }

  fin.close();

  cout << "*** VCCTL file read successfully, found " << vcctlid.size()
       << " unique ids" << endl
       << endl;

  // Now sort the vector of VCCTL ids in ascending order

  sort(vcctlid.begin(), vcctlid.end());

  // Next, determine the correspondences by user input

  if (getCorrespondences(vcctlid, idmap)) {
    return 1;
  }

  return 0;
}

void openFile(string &fname, ifstream &fin) {
  cout << "Enter the name of the VCCTL image file: ";
  cin >> fname;

  fin.open(fname.c_str());
  if (!fin) {
    cout << endl;
    cout << "ERROR:  Cannot open file " << fname << endl;
    cout << "        Check the file name and path" << endl << endl;
    exit(1);
  }

  return;
}
int readVCCTLHeader(ifstream &fin, int &xsize, int &ysize, int &zsize,
                    float &res) {
  string buff;
  string ver;

  fin >> buff;

  if (buff.compare(VCCTL_Version_string) == 0) {
    fin >> ver;
    if (Verbose)
      cout << "Verbose:" << buff << ver << endl;
    fin >> buff;
    if (Verbose)
      cout << "Verbose:" << buff << endl;
    if (buff.compare(VCCTL_Xsize_string) == 0) {
      fin >> xsize;
      if (Verbose)
        cout << "Verbose:" << xsize << endl;
      fin >> buff >> ysize;
      if (Verbose)
        cout << "Verbose:" << buff << ysize << endl;
      fin >> buff >> zsize;
      if (Verbose)
        cout << "Verbose:" << buff << zsize << endl;
      fin >> buff >> res;
      if (Verbose)
        cout << "Verbose:" << buff << res << endl;
    } else if (buff.compare(VCCTL_ImgSize_string)) {
      fin >> xsize;
      ysize = zsize = xsize;
      res = 1.0;
    }

  } else {
    cout << endl;
    cout << "ERROR: Invalid file format" << endl << endl;
    return 1;
  }

  return 0;
}

int getCorrespondences(vector<int> vcctlid, map<int, int> &corr) {
  corr.clear();
  int i = 0;
  string strinput;
  int idval;
  bool validinput = false;

  cout << "*** Now look at your chemistry.xml file to determine" << endl
       << "*** which THAMES microstructure phase should correspond" << endl
       << "*** to each of the following VCCTL phases:" << endl
       << endl;

  while (i < vcctlid.size()) {
    cout << "%% Found VCCTL phase \"" << Vcctlnames[vcctlid[i]] << "\"" << endl;
    validinput = false;
    do {
      cout << "---> Enter the THAMES id for this phase: ";
      cin >> strinput;
      if (isNaturalNumber(strinput)) {
        idval = stoi(strinput);
        validinput = true;
      }
      if (!validinput) {
        cout << ":( Sorry, " << strinput
             << " is not a natural number. Try again." << endl;
        cout << "    ---> Enter the THAMES id for this phase: ";
      }
    } while (!validinput);

    corr.insert(make_pair(vcctlid.at(i), idval));
    i++;
  }

  cout << "%% All done!  Thanks for your help." << endl << endl;
  return 0;
}

bool isNaturalNumber(string &str) {
  bool isnat = true;

  for (int i = 0; i < str.length(); ++i) {
    if (isdigit(str[i]) == false) {
      return false;
    }
  }
  return true;
}

int writeTHAMESImage(string input_filename, const int xsize, const int ysize,
                     const int zsize, const float res, vector<int> mic,
                     map<int, int> idmap) {
  // Create the output file stream

  string output_filename = input_filename + ".thames";
  cout << "*** Writing the THAMES image, will be called " << output_filename
       << endl;
  ofstream fout(output_filename.c_str());
  if (!fout) {
    cout << endl;
    cout << "ERROR:  Cannot open file " << output_filename << endl;
    cout << "        Check the file name and path" << endl << endl;
    return 1;
  }

  if (writeTHAMESHeader(fout, xsize, ysize, zsize, res)) {
    return 1;
  }

  map<int, int>::iterator it;
  for (int i = 0; i < mic.size(); ++i) {
    it = idmap.find(mic.at(i));
    if (it != idmap.end()) {
      fout << it->second << endl;
    } else {
      cout << endl;
      cout << "ERROR: Could not find VCCTL id mapping" << endl;
      fout.close();
      if (remove(output_filename.c_str()) != 0) {
        cout << "ERROR: Could not remove partial THAMES file" << endl;
      }
      cout << endl;
      return 1;
    }
  }

  fout.close();

  cout << "*** Done!  Exiting now." << endl << endl;

  return 0;
}

int writeTHAMESHeader(ofstream &fout, const int xsize, const int ysize,
                      const int zsize, const float res) {
  fout << THAMES_Version_string << " " << THAMES_Version << endl;
  fout << THAMES_Xsize_string << " " << xsize << endl;
  fout << THAMES_Ysize_string << " " << ysize << endl;
  fout << THAMES_Zsize_string << " " << zsize << endl;
  fout << THAMES_ImgRes_string << " " << fixed << setprecision(4) << res
       << endl;

  return 0;
}
