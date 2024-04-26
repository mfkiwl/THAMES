# INSTALL

## INTRODUCTION

This file describes how to install THAMES.  The instructions in this file
are for the most common use cases, and cover the command line tools.

For further information, or in case of problems, please contact the author,
Jeff Bullard (jeffrey.bullard@nist.gov).

The API documentation of THAMES is available in PDF and HTML formats.  It
is not bundled with the software but can be created during the installation
of the software.

For more information about THAMES, see the API documentation overview, or
the accompanying user guide (coming soon).  More information about THAMES and
its applications can also be found in the following references:

    * Bullard, J.W., Lothenbach, B., Stutzman, P.E., Snyder, K.A., Coupling thermodynamics and digital image models to simulate hydration and microstructure development of portland cement pastes, _Journal of Materials Research_ 26, (2011) 609-622.

    * Feng, P., Miao, C., Bullard, J.W., A model of phase stability, microstructure and properties during leaching of portland cement binders, _Cement and Concrete Composites_ 49, (2014) 9-19.

    * Li, X., Grasley, Z.C., Garboczi, E.J., Bullard, J.W., Modeling the apparent intrinsic viscoelastic relaxation of hydrating cement paste, _Cement and Concrete Composites_ 55, (2014) 322-330.

    * Feng, P., Garboczi, E.J., Miao, C., Bullard, J.W., Microstructural origins of cement paste degradation by external sulfate attack, _Construction and Building Materials_ 96, (2015) 391-403.

    * Li, X., Grasley, Z.C., Bullard, J.W., Garboczi, E.J., Computing the time evolution of the apparent viscoelastic/viscoplastic Poisson's ratio of hydrating cement paste, _Cement and Concrete Composites_, 56 (2015) 121-133.


## PREREQUISITES

* GNU Compiler Collection (gcc/g++ >= 11, REQUIRED)
 * Mac OS no longer uses true gcc/g++, despite the fact that
 Apple provides shortcuts called gcc and g++ to its own Apple Clang
 compiler.  So on a Mac you _must_ install the actual gcc/g++ and
 make it the default compiler, at least long enough to build THAMES.
 Instructions are provided in the section below for building on Mac OS.

* CMake (>= 3.5), the build system used by THAMES
 * Required for building THAMES

* Doxygen (>= 1.8.13), the API documentation software
 * Required for creating the API documentation

* LaTeX 2e, the document preparation system
 * Required only for creating the PDF version of the API documentation

* GEM-Selektor (optional)
 * Needed for producing input data files for thermodynamic calculations.
 * Source code available at https://bitbucket.org/gems4/gems3gui.git
 * Binaries available at https://gems.web.psi.ch

* libxml2

* Imagemagick (optional, for creating PNG image files for the microstructures)

## BUILDING

This section assumes that all the prerequisites are already installed on
your system.  We also assume that your THAMES distribution is installed at /PathToTHAMES/THAMES.

### Mac OS

As already stated, modern Mac OS uses the Apple Clang compiler by default,
but THAMES requires the GNU Compiler Collection (gcc).
Therefore you must install the Gnu compiler suite and then
make it the default C/C++ compiler on your computer. Assuming you
use Homebrew for package management:

* `brew install gcc`
* `cd /opt/homebrew/bin`
* `ln -s gcc-13 gcc`
* `ln -s g++-13 g++`
* Edit your path to ensure that `/opt/homebrew/bin` comes before `/usr/bin`

Proceed to section for Unix below to complete the build.

### Unix

First build and install the GEMS3K library:

1. cd /PathToTHAMES/THAMES/src/GEMS3K-standalone
2. ./install.sh

Next, build and install the z compression library:

1. cd /PathToTHAMES/THAMES/src/zlib
2. `mkdir build`
3. `cd build`
4. `cmake ..`
5. `make`
6. `cp zconf.h ../../Resources/include/.`
7. `cp libz.a ../../Resources/lib/.`
8. `cp libz.dylib ../../Resources/lib/.`

Next, build and install the png library:

1. cd /PathToTHAMES/THAMES/src/libpng
2. `mkdir build`
3. `cmake .. -DZLIB_ROOT=../../Resources`
4. `make`
5. `cp *.dylib ../../Resources/lib/.`
6. `cp *.a ../../Resources/lib/.`
7. `cp *.h ../../Resources/include/.`

Next, build and install THAMES. The recommended way to configure THAMES is to do an out-of-source
build, which means that the original files and directories are left untouched.
Doing this makes the re-compiling and cleaning of the installation files
much simpler.

1. cd /PathToTHAMES/THAMES/build
2. `cmake ..`
3. `make`
4. `make install`
5. `make doc`

This will install the "thames" executable in the /PathToTHAMES/bin directory, and the
static libraries in the PathToTHAMES/lib directory.

### Windows

Windows does not come prebuilt with any kind
of system for compiling C/C++ code.  Therefore
you must first install **MinGW** and **MSYS**. Assuming these are installed, execute the following steps.

1. Open mingw64 shell at C:\msys64\mingw64.exe
2. Build GEMS3K library
	* cd /PathToTHAMES/THAMES/src/GEMS3K-standalone
	* mkdir build
	* cd build
	* cmake .. -G "MinGW Makefile" -DCMAKE\_CSS\_FLAGS=-fPIC -DCMAKE\_BUILD\_TYPE=Release -DCMAKE\_INSTALL\_PREFIX=../Resources
	* mingw32-make.exe
	* mingw32-make.exe install
3. Build and install the z compression library:
    * cd /PathToTHAMES/THAMES/src/zlib
    * `mkdir build`
    * `cd build`
    * `cmake .. -G "MinGW Makefile"`
    * `mingw32-make.exe`
    * `cp zconf.h ../../Resources/include/.`
    * `cp libz.a ../../Resources/lib/.`
    * `cp libz.dylib ../../Resources/lib/.`
4. Build and install the png library:
    * cd /PathToTHAMES/THAMES/src/libpng
    * `mkdir build`
    * `cmake .. -G "MinGW Makefile" -DZLIB_ROOT=../../Resources`
    * `mingw32-make.exe`
    * `cp *.dylib ../../Resources/lib/.`
    * `cp *.a ../../Resources/lib/.`
    * `cp *.h ../../Resources/include/.`
5. Build THAMES
	* cd /PathToTHAMES/THAMES/build
	* cmake .. -G "MinGW Makefile" -DCMAKE\_BUILD\_TYPE=Release -DCMAKE\_INSTALL\_PREFIX=../bin
	* mingw32-make.exe
	* mingw32-make.exe install
	* mingw32-make.exe doc

This will install the "thames" executable in the /PathToTHAMES/bin directory, and the
static libraries in the PathToTHAMES/lib directory.

## UNINSTALLING

To uninstall everything except the original source files, simply delete
recursively everything in the build/ directory.
