#!/bin/bash

#./install-dependencies.sh

threads=1
BRANCH_GEMS3K=trunk
BuildType=Release
#BuildType=Debug
InstallPrefix=../../Resources
#InstallPrefix=/home/sveta/tmp
workfolder=${PWD}


mkdir -p $InstallPrefix
mkdir -p build
cd build
cmake .. -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=$BuildType -DCMAKE_INSTALL_PREFIX=$InstallPrefix 

# For new Mac OS, something like below is necessary to specify the GNU compiler,
# so comment the above line and uncomment the one below

#cmake .. -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++ -DCMAKE_CXX_FLAGS=-fPIC -DCMAKE_BUILD_TYPE=$BuildType -DCMAKE_INSTALL_PREFIX=$InstallPrefix 
make -j $threads 
make install

if [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
   sudo ldconfig
fi
