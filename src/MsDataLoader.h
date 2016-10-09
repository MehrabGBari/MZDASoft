///////////////////////////////////////////////////////////////////////////////
//  MZDASoftPPE v.1.0 is under the GPL Version 2 License.
//  October 2014.
//  Please refer to the included MZDASoftPPELicense.txt file
//  for the full text of license.
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License version 2 
//    as published by the Free Software Foundation.
//    http://www.gnu.org/licenses/gpl-2.0.txt
// 
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License Verion 2 
//    along with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//
//    Contact Developer(s) at:
//    email: michelle.zhang@utsa.edu
//    web: http://compgenomics.utsa.edu/zgroup/MZDASoft
//    tel: 210-458-6856
//    mail:
//    Department of Electrical and Computer Engineering, BSE 1.328
//    The University of Texas at San Antonio
//    One UTSA Circle
//    San Antonio, TX 78249
//
//    In collaboration with the RCMI Computational Biology Initiative/
//    Computational Systems Biology Core at the University of Texas at 
//    San Antonio. [www.cbi.utsa.edu]
//
//    Grant Acknowledgments:
//    "This work received computational support from 
//    Computational System Biology Core, funded by the National Institute 
//    on Minority Health and Health Disparities (G12MD007591) from the 
//    National Institutes of Health."
//
// Copyright (c) 2014, The University of Texas at San Antonio. All rights reserved.
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/// MsDataLoader.h
///
/// Load the following data types:
///  1) mzXML ( raw data )
///  2) MSDA Level 1 ( Pre-LC Regions of Interest )
///  3) MSDA Level 2 ( LC Candidates )
///  4) MSDA Level 3 ( Peptide Candidates )
///
/// Created on: Aug 3, 2011
///     Author: nelson.ramirez
///////////////////////////////////////////////////////////////////////////////

#ifndef MSDATALOADER_H_
#define MSDATALOADER_H_

#include "ScanGroup.h"
#include "mzParser.h"  // Ramp library for mzXML reading
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

namespace msda {

class MsDataLoader {

public:
	
	///////////////////////////////////////////////////////////////////////////
	///
	/// Constructor.
	/// All work to be done within the specific loader for each type
	/// of Mass Spectrometry data file.
	///
	///////////////////////////////////////////////////////////////////////////
	MsDataLoader();	
	
	///////////////////////////////////////////////////////////////////////////
	///
	/// Destructor
	///
	///////////////////////////////////////////////////////////////////////////
	~MsDataLoader();	

	
	///////////////////////////////////////////////////////////////////////////
	///
	/// mzXML File Loader.
	/// 
	/// Load the ms data
	/// While we can load any mzXML file, the goal of this system is 
	/// to allow loading scan groups, that is sets of scans with some 
	/// overlap between the scans as controlled by a parallel work 
	/// distribution manager that is over the calling of this executable.
	///
	/// We can take any set of user selection criteria regarding the 
	/// subset of the original dataset to process in each workers.
	/// Multiple workers can all be reading from the file in parallel.
	/// Each worker populates a ScanGroup class, which will hold 
	/// all relevant information needed to proceed with additional
	/// analysis.
	///
	///////////////////////////////////////////////////////////////////////////
	int loadMzXmlFile(
			string filename,
			ScanGroup & sg,
			vector<int> & scanList, 
			int level, 
			double minMz, 
			double maxMz);

	///////////////////////////////////////////////////////////////////////////
	///
    /// Level 1 Format Data Loader
	///
    /// This function takes a file in Level 1 format, which
	/// contains the pre-LC regions and loads it into memory.
	///
    /// Level 1 format data can be in 2 formats, either text or binary.
	/// 
	/// 
	/// Inputs:
    /// 1) Level 1 Format version string: "1.0.0","1.0.1",...
	/// 2) Type string: "text" or "binary"
	/// 3) Data selector string: "all","subset"
	/// 4) Start Pre-LC Region ID ( only if Data selector string is set to "subset"
	///    - Note: Only used if the selector variable is set to "all"
	/// 5) End Pre-LC Region ID ( only if Data selector string is set to "subset" )
	///    - Note: Start Pre-LC Region ID and End Pre-LC Region ID are 
	///      only used if the selector variable is set to "all"
    /// 6) Input filename of the Level 1 file ( .csv )
	/// 7) Reference to the ScanGroup Object that will contain all original data
    ///    - Note:  For both reading the Level 1 data from a file and
    ///             loading Level 1 data from memory, there are two
	///             workflows for continuing the processing of the 
    ///             data, from file or from the in memory data structures.
    ///
	///             
	///////////////////////////////////////////////////////////////////////////
    int loadLevel1File(
				string version,
				string type,
				string selector,
				int startPLCRegionId,
				int endPLCRegionId,
				string filename,
				ScanGroup & sg);

private:
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Use the RAMP library(LGPL library) to access mzXML data directly
	// The RAMP library interfacing data structures
	// We must dynamically link the RAMP library.  Only header files
	// and their definitions are allowed to be used within the source code.
	//  
	//
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	// Pointers
	///////////////////////////////////////////////////////////////////////////
	RAMPFILE *mzXML_file;  // File structure for an mzXML file in RAMP library
	ramp_fileoffset_t *pScanIndex;  // Pointer to the scan index data structure
	InstrumentStruct* pInstrStruct; // Information about the instrument
	RAMPREAL *pPeaks;  // contains the data (mz,int)
	
	///////////////////////////////////////////////////////////////////////////
	// Non-Pointers
	///////////////////////////////////////////////////////////////////////////
	ramp_fileoffset_t indexOffset;  // Used access any specific scan in file
	int iLastScan;  // id of the last scan
	struct RunHeaderStruct runHeader;  // Global information about the mzXML file
	struct ScanHeaderStruct scanHeader; // Information about each individual scan
	
	int version;  // Version of the MsDataLoader Class
	
};

} // namespace msda
#endif // MSDATALOADER_H_ 
