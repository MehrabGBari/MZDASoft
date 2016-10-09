#!/bin/bash
#//////////////////////////////////////////////////////////////////////////////
#//  MZDASoftPPE v.1.0 is under the GPL Version 2 License.
#//  October 2014.
#//  Please refer to the included MZDASoftPPELicense.txt file
#//  for the full text of license.
#//
#//    This program is free software; you can redistribute it and/or modify
#//    it under the terms of the GNU General Public License version 2 
#//    as published by the Free Software Foundation.
#//    http://www.gnu.org/licenses/gpl-2.0.txt
#// 
#//    This program is distributed in the hope that it will be useful,
#//    but WITHOUT ANY WARRANTY; without even the implied warranty of
#//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#//    GNU General Public License for more details.
#//
#//    You should have received a copy of the GNU General Public License Verion 2 
#//    along with this program; if not, write to the Free Software Foundation, Inc.,
#//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#//
#//    Contact Developer(s) at:
#//    email: michelle.zhang@utsa.edu
#//    web: http://compgenomics.utsa.edu/zgroup/MZDASoft
#//    tel: 210-458-6856
#//    mail:
#//    Department of Electrical and Computer Engineering, BSE 1.328
#//    The University of Texas at San Antonio
#//    One UTSA Circle
#//    San Antonio, TX 78249
#//
#//    In collaboration with the RCMI Computational Biology Initiative/
#//    Computational Systems Biology Core at the University of Texas at 
#//    San Antonio. [www.cbi.utsa.edu]
#//
#//    Grant Acknowledgments:
#//    "This work received computational support from 
#//    Computational System Biology Core, funded by the National Institute 
#//    on Minority Health and Health Disparities (G12MD007591) from the 
#//    National Institutes of Health."
#//
#// Copyright (c) 2014, The University of Texas at San Antonio. All rights reserved.
#//////////////////////////////////////////////////////////////////////////////

###############################################################################
# BUILD MZDASOFT ALGORITHMS SOURCE CODE
###############################################################################
cd src
CBILIBPATH="../../cbilib"
MSTOOLKITPATH="../mstoolkit-read-only"
MZPARSERPATH="../mstoolkit-read-only/mzParser"
COMPILERFLAGS="-O3 -fopenmp -fPIC -DMATLABALGORITHMVALIDATIONDEBUG"
RUNTIMELIBRARIES="-lcbi -lmstoolkitlite -lre2"
g++ -c -o ScanGroup.o ScanGroup.cpp  -I$CBILIBPATH/include $COMPILERFLAGS -DOPENMP -DMSDADEBUG -DSCANGROUPDEBUG -DFINDPRELCCANDIDATESDEBUG
g++ -c -o Utilities.o Utilities.cpp -I$CBILIBPATH/include $COMPILERFLAGS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -DOPENMP -DMSDADEBUG -DSCANGROUPDEBUG -DFINDPRELCCANDIDATESDEBUG
g++ -c -o MsDataLoader.o MsDataLoader.cpp -I$CBILIBPATH/include -I$MSTOOLKITPATH -I$MZPARSERPATH $COMPILERFLAGS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -DOPENMP -DMSDADEBUG -DSCANGROUPDEBUG -DFINDPRELCCANDIDATESDEBUG
g++ -c -o Workflow.o Workflow.cpp -I$CBILIBPATH/include -I$MSTOOLKITPATH -I$MZPARSERPATH $COMPILERFLAGS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -DOPENMP -DMSDADEBUG -DSCANGROUPDEBUG -DFINDPRELCCANDIDATESDEBUG
g++ -c -o msda.o msda.cpp -I$CBILIBPATH/include -I$MSTOOLKITPATH -I$MZPARSERPATH $COMPILERFLAGS -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -DOPENMP -DMSDADEBUG -DSCANGROUPDEBUG -DFINDPRELCCANDIDATESDEBUG
g++ -o mzda msda.o ScanGroup.o MsDataLoader.o Utilities.o Workflow.o $RUNTIMELIBRARIES -I$CBILIBPATH/include -I$MSTOOLKITPATH -I$MZPARSERPATH $COMPILERFLAGS -L$CBILIBPATH/lib -L$MSTOOLKITPATH -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -DOPENMP -DMSDADEBUG -DSCANGROUPDEBUG -DFINDPRELCCANDIDATESDEBUG 

###############################################################################
# BUILD SHARED MZDA LIBRARY OBJECTS
###############################################################################
g++ -shared -Wl,-soname,libmzda.so -o libmzda.so Utilities.o

###############################################################################
# INSTALL MZDALIB
###############################################################################
mv libmzda.so ../lib
mv mzda ../bin
cp ./*.h ../include
rm ./*.o
###############################################################################
# BUILD MSDA UNIT TESTS
###############################################################################
cd ../tests
cd Utilities
make;
