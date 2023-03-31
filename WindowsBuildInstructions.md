# Building THAMES using MinGW
1. Install MSYS from https://www.msys2.org/
 * Direct link to the installer is: https://github.com/msys2/msys2-installer/releases/download/2022-06-03/msys2-x86_64-20220603.exe
 * Use default installation folder
 
2. Start MSYS2 MinGW x64 from Start menu
 * $ pacman -Syu
 * window will close automatically
 
3. Start MSYS2 MSYS from Start menu
 * $ pacman -Syu
 * $ pacman -S --needed base-devel mingw-w64-x86_64-toolchain
 * $ pacman -S mingw-w64-x86_64-cmake
 * $ pacman -S mingw-w64-x86_64-make
 * $ pacman -S git
 * close window
 
4. Get THAMES source code
 * Start MSYS2 MinGW x64 from Start menu
 * $ git clone https://github.tamu.edu/jwbullard/THAMES.git
 * Enter userid and password
 
5. Start building GEMS3K-Standalone
 * $ cd THAMES/src/GEMS3K-Standalone
 * $ mkdir build
 * $ cd build
 * $ cmake .. -G "MinGW Makefiles" -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local
 * $ /mingw64/bin/mingw32-make.exe
 * $ /mingw64/bin/mingw32-make.exe install
 
6. Build THAMES
 * Edit Solution.h, line 19 as: #include "../GEMS3K-standalone/GEMS3K/node.h"
 * Edit ChemicalSystem.h line 34 as: #include "../GEMS3K-standalone/GEMS3K/node.h"
 * Edit vcctl2thames.cc
    - Add #include <algorithm>
 * $ cd ../../../build          or back out and go to THAMES/build folder
 * $ cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local
 * $ /mingw64/bin/mingw32-make.exe
 * $ /mingw64/bin/mingw32-make.exe install
 
 7. Edit Path to include C:\msys64\mingw64\bin
  * Start Edit the System Variables form the Start menu
  * Click "Environment Variables" button
  * Under "System Variables", scroll down until you see Path
  * Select that line, and press "Edit"
  * Edit Environment variable window will come up
  * Press "New"
  * Add C:\msys64\mingw64\bin to the bottom
  * OK --> OK --> OK
 
thames.exe should be under C:\msys64\home\<username>\THAMES\bin
 
# Building THAMES under WSL2
1. Install WSL
 * Start Powershell using Start menu - run it as admin
 * PS> wsl --install
 * PS> wsl --install -d Ubuntu-20.04
 * Set up user information
 
2. Start Ubuntu from Start menu
 * $ sudo apt update && sudo apt upgrade
 * You will probably need all the packages/libraries: cmake, libxml2, zlib, etc. If a command is not found, it needs to be installed.
   - You can get them using the following command: sudo apt-get install <package_name>
   - $ sudo apt install cmake
   - $ sudo apt-get install libxml2-dev
   - $ sudo apt-get install doxygen doxygen-doc doxygen-gui graphviz   (optional)
   - $ sudo apt install texlive-latex-extra                            (optional)
 
3. Get THAMES source code
 * $ git clone https://github.tamu.edu/jwbullard/THAMES.git
 * Enter userid and password
 
4. Start building GEMS3K-standalone
 * $ cd THAMES/src/GEMS3K-standalone
 * $ sudo ./install.sh
 
5. Build THAMES
 * Edit vcctl2thames.cc
    - Add #include <algorithm>
 * $ cd ../../build          or back out and go to THAMES/build folder
 * $ cmake ..
 * $ make
 * $ make install
 
 thames.exe should be under THAMES/bin and can be started as:
 $ ./thames           or
 $ ./thames /mnt/c/temp/<input.in> /mnt/c/temp/outout.out    where the paths are the location of the input files
