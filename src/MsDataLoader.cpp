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
// MsDataLoader.cpp
//
// An MsDataLoader object uses an mzXML reading
// library(RAMP) to a subset of scan into a ScanGroup object.
// For references on using the RAMP library go to:
// The ScanGroup object is composed of a set of ScanLine objects.
// Ideally, each scangroup should overlap a few scan lines, to provide
// for full data analysis coverage.
//
// In order to minimize the contention between workers, we want to 
// release the file as soon as possible.  
//
// http://www.sbeams.org/svn/sbeams/trunk/sbeams/lib/c/Proteomics/getSpectrum/getSpectrumHeader.cpp
// http://sashimi.sourceforge.net/software_glossolalia.html
// http://rss.acs.unt.edu/Rdoc/library/caMassClass/html/mzXML.html
// http://www.sbeams.org/svn/sbeams/trunk/sbeams/lib/c/Proteomics/getSpectrum/ramp/ramp.h
// Visualization:
// Notes: Both gnuplot and matlab 3d visualization are very limited,
// we need an approach that will allow us to visualize every aspect
// of the data.
//														 			Visualization engines
// mzXML --> internal data structures -->  csv files -->R --> mzXML 	-->insilicos viewer
//                                                   					--> matlab
//                                                   					--> gnuplot
// R has a toolkit that can take output in csv format and write mzXML files.
// 
// http://rgm2.lab.nig.ac.jp/RGM2/func.php?rd_id=caMassClass:msc.peaks.IO.mzXML
// http://insilicos.com/products/insilicos-viewer-1
//
//
// Visualization:
// We need to be able to take a scan group and visualize it with insilicos viewer
// Each level of data needs to have the ability to be visualized easily.
//
//
//
//  Created on: Aug 3, 2011
//      Author: nelson.ramirez
///////////////////////////////////////////////////////////////////////////////

#include "MsDataLoader.h"
#include "ScanGroup.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>
#include <sstream>
#include <exception>


using namespace std;


namespace msda {

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
MsDataLoader::MsDataLoader() {
	; // Nothing is to be done at construction time, this is actually 
	  // better than trying to do the reading of files during object 
	  // object construction.  
}


///////////////////////////////////////////////////////////////////////////////
// Destructor
///////////////////////////////////////////////////////////////////////////////
MsDataLoader::~MsDataLoader() {
	// Since we are using vectors for data storage
	// in memory, we need not worry about freeing 
	// memory in this destructor.
	// Dynamic memory is used briefly during the construction
	// phase of this object, with careful handling
	// of potential issues to ensure cleanup of 
	// all dynamic memory alllocated by the RAMP 
	// library.

	// In addition, since this class does NOT claim ownership
	// over any resource such as memory, files, after the 
	// initial construction phase, we want to be absolutely
	// clear that all resources that are temporarily used 
	// during the construction phase, are also released insdie
	// of the constructor.
	;
}


///////////////////////////////////////////////////////////////////////////////
//
// Load an mzXML file into memory ( A ScanGroup ) data structure.
//  
//
//
///////////////////////////////////////////////////////////////////////////////
int MsDataLoader::loadMzXmlFile(string filename,
		ScanGroup & sg, 
		vector<int> & scanList, 
		int level, 
		double minMz, 
		double maxMz){

	///////////////////////////////////////////////////////////////////////////
	// Make sure that we have scans to load
	// throw an exception 
	///////////////////////////////////////////////////////////////////////////
	assert(scanList.size() > 0);   
	if( scanList.size() <= 0 ){ 
		throw "Error in MsDataLoader Constructor: scanList.size() <= 0";
	}
	
	
	///////////////////////////////////////////////////////////////////////////
	// Open the mzXML file for reading
	///////////////////////////////////////////////////////////////////////////
	mzXML_file=rampOpenFile(filename.c_str());
	
	///////////////////////////////////////////////////////////////////////////
	// Make sure the mzXML_file pointer is valid
	///////////////////////////////////////////////////////////////////////////
	if( mzXML_file != 0 ){ 

		try {
			///////////////////////////////////////////////////////////////////
			// Only try to use the mzXML_file pointer if it is valid
			// indexOffset is a 64 bit integer.  This is very important
			// in order to be able to process multi-gigabyte files.
			///////////////////////////////////////////////////////////////////
			indexOffset = getIndexOffset(mzXML_file);  // get the offset into the mzXML index	
			
			///////////////////////////////////////////////////////////////////
			// For an approximately 2GB file, the indexOffset
			// indexOffset=1935769761
			///////////////////////////////////////////////////////////////////
			
			
			///////////////////////////////////////////////////////////////////
			// EXCEPTION TEST
			// It is extremely important to test out what happens in 
			// case of an exception, especially within critical construction
			// phases such as this one.
			///////////////////////////////////////////////////////////////////
			//vector<int> test;
			//test.at(100) =5;
			///////////////////////////////////////////////////////////////////
			// END EXCEPTION TEST
			///////////////////////////////////////////////////////////////////
		} 
		catch( exception &  e){
			///////////////////////////////////////////////////////////////////
			// Perform resource cleanup of all resources we are responsible
			// at this point.
			// The rampCloseFile function performs RAMP library specific
			// cleanup, so this is all we need at this point to cleanup
			// after the exception.
			// For more details refer to ramp.cpp to see exactly what is 
			// getting cleaned up by rampCloseFile.
			///////////////////////////////////////////////////////////////////
			rampCloseFile(mzXML_file);
		
			throw "Error in MsDataLoader Constructor: exception at getIndexOffset(mzXML_file)";
		} // end catch block, handles closing the file

	} // end valid pointer check 


	///////////////////////////////////////////////////////////////////////////
	// Check for a null pointer in the index pointer.
	// The index allows us to reference all other information in the 
	// mzXML file in a random access manner.
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	// Read the scan index into a vector, get LastScan
	// Note: In ramp.cpp, we wee that readIndex returns a pointer which
	// has been allocated with malloc.  Therefore we MUST use free() to 
	// 
	///////////////////////////////////////////////////////////////////////////
	try{
		pScanIndex = readIndex(mzXML_file, indexOffset, &iLastScan);  // must handle null case ( pointer )
		///////////////////////////////////////////////////////////////////
		// EXCEPTION TEST
		// It is extremely important to test out what happens in 
		// case of an exception, especially within critical construction
		// phases such as this one.
		///////////////////////////////////////////////////////////////////
		//vector<int> test;
		//test.at(100) =5;
		///////////////////////////////////////////////////////////////////
		// END EXCEPTION TEST
		///////////////////////////////////////////////////////////////////
	}
	catch( exception & e) {
		///////////////////////////////////////////////////////////////////
		// Perform resource cleanup of all resources we are responsible
		// at this point.
		// The rampCloseFile function performs RAMP library specific
		// cleanup, so this is all we need at this point to cleanup
		// after the exception.
		// For more details refer to ramp.cpp to see exactly what is 
		// getting cleaned up by rampCloseFile.
		///////////////////////////////////////////////////////////////////
		rampCloseFile(mzXML_file);

		throw "Error in MsDataLoader Constructor: exception at readIndex(mzXML_file,indexOffset,&iLastScan)";
	}

	///////////////////////////////////////////////////////////////////////////
	// Read run header information
	// We must first check that that index pointer is valid
	///////////////////////////////////////////////////////////////////////////
	if ( pScanIndex != 0 ) {  // check for the null pointer
		try{
			///////////////////////////////////////////////////////////////////
			// The memory for runHeader has been allocated by the compiler
			// on the stack.
			///////////////////////////////////////////////////////////////////
			//cout<<"runHeader address"<<&(runHeader)<<endl;
			//cout<<"version address"<<&(version)<<endl;
			///////////////////////////////////////////////////////////////////
			// Obtain global information about the run, runHeader structure
			// We only want to try to use the pScanIndex pointer if it is 
			// valid.  
			// If anything goes wrong here, this implies a serious problem
			// and we need to terminate gracefully.
			///////////////////////////////////////////////////////////////////
			readRunHeader(mzXML_file, pScanIndex, &runHeader, iLastScan);
			
			///////////////////////////////////////////////////////////////////
			// Load the data from the run header into temporary ( stack based )
			// data structures 
			///////////////////////////////////////////////////////////////////
			int scanCount = runHeader.scanCount;
			double dEndTime = runHeader.dEndTime;
			double dStartTime = runHeader.dStartTime;
			double endMZ = runHeader.endMZ;
			double highMZ = runHeader.highMZ;
			double lowMZ = runHeader.lowMZ;
			double startMZ = runHeader.startMZ;
			
			///////////////////////////////////////////////////////////////////
			/// DEBUG
			///////////////////////////////////////////////////////////////////
			//cout<<"Run Header Information"<<endl;
			//cout<<"scanCount="<<runHeader.scanCount<<endl;
			//cout<<"dEndTime="<<runHeader.dEndTime<<endl;
			//cout<<"dStartTime="<<runHeader.dStartTime<<endl;
			//cout<<"endMZ="<<runHeader.endMZ<<endl;
			//cout<<"highMZ="<<runHeader.highMZ<<endl;
			//cout<<"lowMZ="<<runHeader.lowMZ<<endl;
			//cout<<"startMZ="<<runHeader.startMZ<<endl;
			
			///////////////////////////////////////////////////////////////////
			// Update the ScanGroup Object
			///////////////////////////////////////////////////////////////////
			sg.addRunHeaderData( scanCount, dEndTime, dStartTime, endMZ, highMZ, lowMZ, startMZ );
			
			///////////////////////////////////////////////////////////////////
			// EXCEPTION TEST
			// It is extremely important to test out what happens in 
			// case of an exception, especially within critical construction
			// phases such as this one.
			///////////////////////////////////////////////////////////////////
			//vector<int> test;
			//test.at(100) =5;
			///////////////////////////////////////////////////////////////////
			// END EXCEPTION TEST
			///////////////////////////////////////////////////////////////////
		}
		catch( exception &  e){


			///////////////////////////////////////////////////////////////////
			// Clean up any resources we are responsible for at
			// this point and let the rest of the program 
			// know that something went wrong via rethrowing 
			// the exception.
			//
			// rampCloseFile calles the operating system's fclose function
			// and then free's the memory for mzXML_file.
			//
			///////////////////////////////////////////////////////////////////
			
			rampCloseFile(mzXML_file); // make sure we close the file if a problem reading the run header
			free(pScanIndex); // free the data allocated for the index data structure
			// The runHeader data structure is allocated on the stack ( no need to use free )
			
			///////////////////////////////////////////////////////
			// Provide some debugging information 
			//
			///////////////////////////////////////////////////////
			cout<<"Exception encountered in MsDataLoader Constructor: Run Header Data loading phase"<<endl;
			cout<<e.what()<<endl;
			///////////////////////////////////////////////////////
			// We want to let the rest of the system know
			// that something wrong, so that appropriate
			// measures can be taken such as closing 
			// the MPI environment if MPI is enabled, etc.
			///////////////////////////////////////////////////////
			throw; // Rethrow the exception after having cleaned up

		} // end catch exception when accessing instrument data structure and updating scan group object

		///////////////////////////////////////////////////////////////////////
		// Successfully used the instrument struct data, clean up the 
		// data structure as soon as possible.
		///////////////////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////////////////
		// We still need to keep the index (pScanIndex) past this point, so we cannot
		// free its memory just yet.
		///////////////////////////////////////////////////////////////////////
		
	} // end processing run header data structure 
	else {
		///////////////////////////////////////////////////////////////////////
		// If the instrument structure pointer was null, there is 
		// something seriously wrong.  We want to throw and exception
		// to allow the error to propagate after first cleaning up 
		// any resources we may have at this point.
		///////////////////////////////////////////////////////////////////////
		throw "Invalid run header data structure exception";
	}


	
	
	///////////////////////////////////////////////////////////////////////////
	// Read the instrument information
	// getInstrumentStruct returns a pointer to memory it allocates 
	// on the heap.
	///////////////////////////////////////////////////////////////////////////
	pInstrStruct =	getInstrumentStruct(mzXML_file);  // must handle null case ( pointer )
	
	
	if ( pInstrStruct != 0 ) {  // check for the null pointer
		try{
			string instrManufacturer( pInstrStruct->manufacturer, INSTRUMENT_LENGTH );  // move to a string data structure
			string instrModel( pInstrStruct->model, INSTRUMENT_LENGTH ); // move to a string data structure 
			string instrIonisation( pInstrStruct->ionisation, INSTRUMENT_LENGTH ); // move to a string data structure
			string instrAnalyzer( pInstrStruct->analyzer,INSTRUMENT_LENGTH); // move to a string data structure
			string instrDetector( pInstrStruct->detector,INSTRUMENT_LENGTH); // move to a string data structure
			///////////////////////////////////////////////////////////
			// Update the ScanGroup Object
			///////////////////////////////////////////////////////////
			sg.addInstrumentData( instrManufacturer, instrModel, instrIonisation, instrAnalyzer, instrDetector );



			///////////////////////////////////////////////////////////////////
			// EXCEPTION TEST
			// It is extremely important to test out what happens in 
			// case of an exception, especially within critical construction
			// phases such as this one.
			///////////////////////////////////////////////////////////////////
			//vector<int> test;
			//test.at(100) =5;
			///////////////////////////////////////////////////////////////////
			// END EXCEPTION TEST
			///////////////////////////////////////////////////////////////////

		}
		catch( exception &  e){

		
			///////////////////////////////////////////////////////
			// Clean up any resources we are responsible for at
			// this point and let the rest of the program 
			// know that something went wrong via rethrowing 
			// the exception.
			///////////////////////////////////////////////////////

			rampCloseFile(mzXML_file); // make sure we close the file if a problem reading the run header
			free(pScanIndex); // free the data allocated for the index data structure ( ramp used malloc to allocate )
			// The runHeader data structure is allocated on the stack ( no need to use free ) 
			free(pInstrStruct); // free the data allocated for the instrument data structure ( ramp used calloc to allocate )
			
			
			///////////////////////////////////////////////////////
			// Provide some debugging information 
			//
			///////////////////////////////////////////////////////
			cout<<"Exception encountered in MsDataLoader Constructor: Instrument Struct data loading phase"<<endl;
			cout<<e.what()<<endl;
			///////////////////////////////////////////////////////
			// We want to let the rest of the system know
			// that something wrong, so that appropriate
			// measures can be taken such as closing 
			// the MPI environment if MPI is enabled, etc.
			///////////////////////////////////////////////////////
			throw; // Rethrow the exception after having cleaned up

		} // end catch exception when accessing instrument data structure and updating scan group object

		
	} // end processing instrument data structure 
	else {
		///////////////////////////////////////////////////////////////////////
		// If the instrument structure pointer was null, there is 
		// something seriously wrong.  We want to throw and exception
		// to allow the error to propagate after first cleaning up 
		// any resources we may have at this point.
		///////////////////////////////////////////////////////////////////////
		throw "Invalid instrument structure exception";
	}
	
	///////////////////////////////////////////////////////////////////////////
	// Test parallel loading of scan lines
	///////////////////////////////////////////////////////////////////////////
	//#ifdef OPENMP
	//#pragma omp parallel 
	//	{
	//		cout<<"Parallel scan line loading test"<<endl;
	//#endif
	//
	//#ifdef OPENMP
	//	}
	//#endif		
	///////////////////////////////////////////////////////////////////////////
	// End Test parallel loading of scan lines
	///////////////////////////////////////////////////////////////////////////
	
	
	///////////////////////////////////////////////////////////////////////////
	// Iterate through each of the items in the scanList as supplied
	// by the input local worker parameter file.
	//
	// In order to minimize memory utilization, since the target of this
	// software is multi-Gigabyte datasets, we need to re-use a temporary
	// buffer where we temporarily hold the data for a single scan line, prior
	// to saving in the scan group class.
	//
	// This leaves the possibility of scanning through the entire file,
	// selecting just a few scans based on user criteria without ever 
	// increasing the total memory utilization of the program.
	//
	///////////////////////////////////////////////////////////////////////////
	vector<int>::size_type i;
	for (  i = 0; i < scanList.size(); i++ ) {
		//for ( int i = 1; i <= runHeader.scanCount; i++ ) {
		
		//cout<<"Looping through all the scan lines:"<<scanList.at(i)<<endl;
		// Read the current scans header information
		// We only want scans that have a level = 1 ( M.Zhang )
		
		
		///////////////////////////////////////////////////////////////////////
		//
		// Try to read the current scan header
		// 
		///////////////////////////////////////////////////////////////////////
		try {
			readHeader(mzXML_file, pScanIndex[scanList.at(i)], &scanHeader);
#ifdef MZDADEBUG
#ifdef MZDATALOADERDEBUG
#ifdef LOADMZXMLFILEDEBUG
			cout<<"*********mzXML Scan Header**********"<<endl;
			cout<<"seqNum="<<scanHeader.seqNum<<endl; // number in sequence observed file (1-based)
			cout<<"acquisitionNum="<<scanHeader.acquisitionNum<<endl;   // scan number as declared in File (may be gaps)
			cout<<"msLevel="<<scanHeader.msLevel<<endl;
			cout<<"peaksCount="<<scanHeader.peaksCount<<endl;
			cout<<"totIonCurrent="<<scanHeader.totIonCurrent<<endl;
			cout<<"retentionTime="<<scanHeader.retentionTime<<endl;        /* in seconds */
			cout<<"basePeakMZ="<<scanHeader.basePeakMZ<<endl;
			cout<<"basePeakIntensity="<<scanHeader.basePeakIntensity<<endl;
			cout<<"collisionEnergy="<<scanHeader.collisionEnergy<<endl;
			cout<<"ionisationEnergy="<<scanHeader.ionisationEnergy<<endl;
			cout<<"lowMZ="<<scanHeader.lowMZ<<endl;
			cout<<"highMZ="<<scanHeader.highMZ<<endl;
#endif
#endif
#endif
			///////////////////////////////////////////////////////////
			// Update the ScanGroup Object with the scan header
			// information
			///////////////////////////////////////////////////////////


			int filePositionInt = scanHeader.filePosition;
			string activationMethodString( scanHeader.activationMethod, SCANTYPE_LENGTH );  // move to a string data structure
			string possibleChargesString( scanHeader.possibleCharges, SCANTYPE_LENGTH ); // move to a string data structure
			string scanTypeString( scanHeader.scanType,SCANTYPE_LENGTH); // move to a string data structure

			///////////////////////////////////////////////////////////////////
			//
			// Save the scan header data structures in the internal
			// Scan Group object.
			// It is very important to be able to go back to the original data
			// in subsequent visualization and data analysis stages.
			//
			// We should immediately save this information to a file,
			// we don't want to store everything in memory, otherwise we 
			// will overwhelm the system.
			///////////////////////////////////////////////////////////////////
			sg.addScanHeaderData( 	
					scanHeader.acquisitionNum,
					scanHeader.mergedScan,
					scanHeader.mergedResultScanNum,
					scanHeader.mergedResultStartScanNum,
					scanHeader.mergedResultEndScanNum,
					scanHeader.msLevel,
					scanHeader.numPossibleCharges,
					scanHeader.peaksCount,
					scanHeader.precursorCharge,
					scanHeader.precursorScanNum,
					scanHeader.scanIndex,
					scanHeader.seqNum,
					filePositionInt,
					scanHeader.basePeakIntensity,
					scanHeader.basePeakMZ,
					scanHeader.collisionEnergy,
					scanHeader.highMZ,
					scanHeader.ionisationEnergy,
					scanHeader.lowMZ,
					scanHeader.precursorIntensity,
					scanHeader.compensationVoltage,
					scanHeader.precursorMZ,
					scanHeader.retentionTime,
					scanHeader.totIonCurrent,
					activationMethodString,
					possibleChargesString,
					scanTypeString );


			///////////////////////////////////////////////////////////////////
			// EXCEPTION TEST
			// It is extremely important to test out what happens in 
			// case of an exception, especially within critical construction
			// phases such as this one.
			///////////////////////////////////////////////////////////////////
			//if ( i == 50) {
			//	vector<int> test;
			//	test.at(100) =5;
			//}
			///////////////////////////////////////////////////////////////////
			// END EXCEPTION TEST
			///////////////////////////////////////////////////////////////////

		}
		catch( exception &  e){


			///////////////////////////////////////////////////////
			// Clean up any resources we are responsible for at
			// this point and let the rest of the program 
			// know that something went wrong via rethrowing 
			// the exception.
			///////////////////////////////////////////////////////

			rampCloseFile(mzXML_file); // make sure we close the file if a problem reading the run header
			free(pScanIndex); // free the data allocated for the index data structure ( ramp used malloc to allocate )
			// The runHeader data structure is allocated on the stack ( no need to use free ) 
			free(pInstrStruct); // free the data allocated for the instrument data structure ( ramp used calloc to allocate )


			///////////////////////////////////////////////////////
			// Provide some debugging information 
			//
			///////////////////////////////////////////////////////
			cout<<"Exception encountered in MsDataLoader Constructor: readHeader# "<<i<<endl;
			cout<<e.what()<<endl;
			///////////////////////////////////////////////////////
			// We want to let the rest of the system know
			// that something wrong, so that appropriate
			// measures can be taken such as closing 
			// the MPI environment if MPI is enabled, etc.
			///////////////////////////////////////////////////////
			throw; // Rethrow the exception after having cleaned up

		} // end catch exception when accessing instrument data structure and updating scan group object

	

		
		
		///////////////////////////////////////////////////////////////////////
		// check the scan header of each non-empty scan
		// only msLevel 1 scans should be used.
		///////////////////////////////////////////////////////////////////////
		if( (scanHeader.msLevel==level) && (scanHeader.peaksCount>0) ) {

			///////////////////////////////////////////////////////////////////
			//cout<<"*********mzXML Scan Header**********"<<endl;
			//cout<<"seqNum="<<scanHeader.seqNum<<endl; // number in sequence observed file (1-based)
			//cout<<"acquisitionNum="<<scanHeader.acquisitionNum<<endl;   // scan number as declared in File (may be gaps)
			//cout<<"msLevel="<<scanHeader.msLevel<<endl;
			//cout<<"peaksCount="<<scanHeader.peaksCount<<endl;
			//cout<<"totIonCurrent="<<scanHeader.totIonCurrent<<endl;
			//cout<<"retentionTime="<<scanHeader.retentionTime<<endl;        /* in seconds */
			//cout<<"basePeakMZ="<<scanHeader.basePeakMZ<<endl;
			//cout<<"basePeakIntensity="<<scanHeader.basePeakIntensity<<endl;
			//cout<<"collisionEnergy="<<scanHeader.collisionEnergy<<endl;
			//cout<<"ionisationEnergy="<<scanHeader.ionisationEnergy<<endl;
			//cout<<"lowMZ="<<scanHeader.lowMZ<<endl;
			//cout<<"highMZ="<<scanHeader.highMZ<<endl;
			///////////////////////////////////////////////////////////////////
			
			///////////////////////////////////////////////////////////////////
			// Read the data for the current scan into a peaks data structure
			// We need to make sure pPeaks is cleaned up in case of an
			// exception.
			//
			// readPeaks returns a pointer allocated using malloc.
			// The caller of the function owns the memory returned 
			// by the readPeaks function.
			// We'll want to add a try catch block here( but only
			// if absolutely necessary, since it implies a performance
			// penalty).
			///////////////////////////////////////////////////////////////////
			pPeaks = readPeaks(mzXML_file, pScanIndex[scanList.at(i)]);
			 
			///////////////////////////////////////////////////////////////////
			// Check for the NULL pointer
			// If the readPeaks functions returns data, then proceed
			// to process it.
			///////////////////////////////////////////////////////////////////
			if ( pPeaks != 0 ) {

				try {
					//////////////////////////////////////////////////////////
					// Read the data for the current scan
					// This data is recycled for each new scan
					//////////////////////////////////////////////////////////
					vector<double> mzdata;
					vector<double> intensitydata;

					///////////////////////////////////////////////////////////
					// totalcount is the total number of data elements,
					// where mz,int,mz,int... is the format
					///////////////////////////////////////////////////////////
					int totalcount = 2*scanHeader.peaksCount;  //total number of data elements

					
					///////////////////////////////////////////////////////////
					// Populate input memory
					///////////////////////////////////////////////////////////
					for ( int k = 0; k < (totalcount-1); k=k+2 ) {
						///////////////////////////////////////////////////////
						// Only load data that meets range requirements
						// High resolution data has many data points that
						// are zero.
						// 08/03/11 meeting with M.Zhang/Z.Wang
						// Discard data points with zero intensity.
						///////////////////////////////////////////////////////
						
						//if( (pPeaks[k] >= minMz) &&   (pPeaks[k] <= maxMz)) {
						if( (pPeaks[k] >= minMz) &&   (pPeaks[k] <= maxMz) &&   (pPeaks[k+1] > 0)) {
							mzdata.push_back(pPeaks[k]);
							intensitydata.push_back(pPeaks[k+1]);
							
							// Data output
							//cout<<scanHeader.retentionTime<<" "<<pPeaks[k]<<" "<<pPeaks[k+1]<<endl;
						}
					}
					
					
					///////////////////////////////////////////////////////////
					// Update the ScanGroup Object
					// Test case:  checked on 11/22/11.
					// Adding an if ( i< 500 ) or some other number 
					// and checking the memory utilization on the system
					// on which this is running.
					// if ( i < 500 ) {  sg.AddScan.. ); }
					///////////////////////////////////////////////////////////
					
					sg.addScan(mzdata, intensitydata, scanHeader.retentionTime, scanHeader.acquisitionNum ); 
					
					
					//cout<<"ScanGroup Object constructed"<<endl;
					//cout<<"sg.size()="<<sg.size()<<endl;
					
					///////////////////////////////////////////////////////////
					// Update the ScanGroup Object
					///////////////////////////////////////////////////////////
					
					
					///////////////////////////////////////////////////////////
					// EXCEPTION TEST
					// It is extremely important to test out what happens in 
					// case of an exception, especially within critical construction
					// phases such as this one.
					///////////////////////////////////////////////////////////
					//vector<int> test;
					//test.at(100) =5;
					///////////////////////////////////////////////////////////
					// END EXCEPTION TEST
					///////////////////////////////////////////////////////////

					
				}
				catch( exception &  e){


					///////////////////////////////////////////////////////
					// Clean up any resources we are responsible for at
					// this point and let the rest of the program 
					// know that something went wrong via rethrowing 
					// the exception.
					///////////////////////////////////////////////////////
					rampCloseFile(mzXML_file); // make sure we close the file if a problem reading the run header
					free(pScanIndex); // free the data allocated for the index data structure ( ramp used malloc to allocate )
					// The runHeader data structure is allocated on the stack ( no need to use free ) 
					free(pInstrStruct); // free the data allocated for the instrument data structure ( ramp used calloc to allocate )
					free(pPeaks);  // free the data allocated by the readPeaks function ( ramp used malloc to allocate )

					///////////////////////////////////////////////////////
					// Provide some debugging information 
					//
					///////////////////////////////////////////////////////
					cout<<"Exception encountered in MsDataLoader Constructor during readPeaks, and sg.addScan block"<<endl;
					cout<<e.what()<<endl;
					///////////////////////////////////////////////////////
					// We want to let the rest of the system know
					// that something wrong, so that appropriate
					// measures can be taken such as closing 
					// the MPI environment if MPI is enabled, etc.
					///////////////////////////////////////////////////////
					throw; // Rethrow the exception after having cleaned up

				}

				///////////////////////////////////////////////////////////////////
				// Read the data for the current scan into a peaks data structure
				// This is the normal program path flow.
				//  
				// To avoid having duplicate memory, we use pPeaks only
				// and exclusively for a single scan.  The data for a scan 
				// is immediately transferred to a properly managed vector
				// object.  At this point, the C style memory allocated by the
				// RAMP library method is de-allocated so that it can 
				// be reused for temporarily interfacing the with C RAMP library.
				// Ramp function readPeaks uses malloc to allocate the data
				//
				//
				// CRITICAL Memory De-allocation
				// A test for this is to search through an entire multi-gigabyte
				// dataset and keep only 1 scan, however, if we do not
				// free the memory every time we read in another scan, we will
				// quickly run out of memory and the system will crash.
				//
				///////////////////////////////////////////////////////////////////
				free(pPeaks);  // Critical de-allocation of memory 



			}
			else {
				///////////////////////////////////////////////////////////////
				// There is no data to read for this scan
				// provide some indication that this scan was invalid
				// in the ScanGroup Data structure.
				// To avoid losing important original file data information,
				// we keep a record of each scan header in the scan group
				// object.
				///////////////////////////////////////////////////////////////
			}
			
			
			
		} // end if check on the scan line properties

	}  // End loop through all scan lines to process
	

	///////////////////////////////////////////////////////////////////////////
	// 
	// Important resource allocation management issue:
	// Since the MsDataLoader class is meant only to interface with the 
	// C RAMP library which uses malloc/free, and all resources allocated
	// in this constructor are merely used as temporary working buffers
	// used to allocate the ScanGroup object, all resources MUST 
	// be given back before exiting this constructor.
	//
	// Note:    Destructors are only called after a complete object
	//          construction, which means that since the internal
	//          data structures in the MsDataLoader class are only
	//          to be used during the constructor call, there is no
	//          reason to try to free resources in the destructor.
	//          This is a very special case, designed very carefully
	//          with the specific intent to handle very lage 
	//          multi-gigabyte files, with full file scanning without
	//          utilizing more than a single scan's memory allocation
	//          at any time.
	//
	// Note:    Before making changes to this approach, the issue of 
	//          processing multi gigabyte datasets must be considered.
	//
	// With this class we handle both low memory utilization of huge file
	// parsing and analysis as well as issues with interfacing to a C
	// parsing library that doesn't manage its memory.
	//
	//
	///////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////
	// Clean up any resources we are responsible for at
	// this point and let the rest of the program 
	// know that something went wrong via rethrowing 
	// the exception.
	///////////////////////////////////////////////////////////////////////////
	rampCloseFile(mzXML_file); // make sure we close the file if a problem reading the run header
	free(pScanIndex); // free the data allocated for the index data structure ( ramp used malloc to allocate )
	// The runHeader data structure is allocated on the stack ( no need to use free ) 
	free(pInstrStruct); // free the data allocated for the instrument data structure ( ramp used calloc to allocate )
	///////////////////////////////////////////////////////////////////////////
	// Important note: pPeaks must NOT be freed here, since it has already 
	// been 
	///////////////////////////////////////////////////////////////////////////
					

	return 0;
	
} // End loadMzXmlFile


///////////////////////////////////////////////////////////////////////////////
///
/// MZDA Level 1 Format Data Loader
///
/// This function takes a file in MZDA Level 1 format, which
/// contains the pre-LC regions and loads it into memory.
///
/// MZDA Level 1 format data can be in 2 formats, either text or binary.
/// 
/// 
/// Inputs:
/// 1) MZDA Level 1 Format version string: "1.0.0","1.0.1",...
/// 2) Type string: "text" or "binary"
/// 3) Data selector string: "all","subset"
/// 4) Start Pre-LC Region ID ( only if Data selector string is set to "subset"
/// 5) End Pre-LC Region ID ( only if Data selector string is set to "subset" )
///    - Note: Start Pre-LC Region ID and End Pre-LC Region ID are 
///      only used if 
/// 6) Input filename of the MZDA Level 1 file ( .csv )
/// 7) Reference to the ScanGroup Object that will contain all original data
///    - Note:  For both reading the MZDA Level 1 data from a file and
///             loading MZDA Level 1 data from memory, there are two
///             workflows for continuing the processing of the 
///             data.
///
/// Outputs:
/// 1) Return code: 
///      0 successful, MZDA Level 1 data has been loaded into
///      -1 invalid version string
///      -2 invalid type string
///      -3 invalid selector string
///      -4 invalid pre-lc start region id
///      -5 invalid pre-lc end region id
///      -6 invalid input filename
///      -7 invalid refernece to ScanGroup object
///
///
///////////////////////////////////////////////////////////////////////////////
int MsDataLoader::loadLevel1File(
		string version,
		string type,
		string selector,
		int startPLCRegionId,
		int endPLCRegionId,
		string filename,
		ScanGroup & sg)  {

	///////////////////////////////////////////////////////////////////////////
	///
    /// Load the MZDA Level 1 File into the MZDA Level 1 Data structures,
	/// if they have not aleady been populated with data from the 
	/// raw mzXML data file.
	///
	///////////////////////////////////////////////////////////////////////////
#ifdef MZDADEBUG
#ifdef MZDATALOADERDEBUG
    cout<<"loadMZDALevel1File Testing"<<endl;
#endif
#endif
	
	///////////////////////////////////////////////////////////////////////////
	///
	/// Perform version specific processing
	///
	///////////////////////////////////////////////////////////////////////////
    if ( version == "0.0.0") {
        cout<<"version= 1.0.0"<<endl;
        ; // add later only if needed

	}

	///////////////////////////////////////////////////////////////////////////
	///
    /// Version "1.0.0"
	///
	///////////////////////////////////////////////////////////////////////////
    else if ( version == "1.0.0") {


        cout<<"Testing mzda level 1 loader version 1.0.0"<<endl;

		///////////////////////////////////////////////////////////////////////
		// 
		// There will be 2 formats a binary and a text format.
		// 
		///////////////////////////////////////////////////////////////////////
		if( type == "text") {  
			if (selector == "all") {
				// Open the file
				// Read the data within the file
				string dataLine;  // buffer to hold the data within a line
				vector<string> data; // buffer to hold all data read in 
				                     // 
				// The memory for this string buffer is 
				// located within the heap.
				try{
					///////////////////////////////////////////////////////////
					//
					// Read data using input file stream
					//
					///////////////////////////////////////////////////////////
					ifstream inputfile;  // input file stream
					inputfile.open(filename.c_str(), ios::in);    // open the streams
					///////////////////////////////////////////////////////////
					// First check if the file is valid 
					///////////////////////////////////////////////////////////
					if (inputfile.is_open())  // make sure the file is valid
					{
						///////////////////////////////////////////////////////
						//
						//  Have an iterative loop going through each
						//  pre-LC region. ( Since pre-LC regions
						//  have a newline separator, we can separate 
						//  them via newlines ), according to the 
                        //  MZDA Level 1 text format specification.
						//
						///////////////////////////////////////////////////////
						while ( inputfile.good() )
						{	
							///////////////////////////////////////////////////
							//
							// read one pre-LC region at a time
							//
							///////////////////////////////////////////////////
							getline(inputfile,dataLine); // get the current line from file
						
							// add only if not empty
							if (!dataLine.empty()){
								///////////////////////////////////////////////
								// 
								//
								// We have some data ready to be placed within
								// the ScanGroup object.  Make a 
								// copy of the data in the dataLine buffer
								// placing it within the ScanGroup object.
								// 
								//
								// Since dataLine is temporary data that will 
								// be lost at soon as the variables
								// dataLine and data go out of scope.
								// 
								// 
								///////////////////////////////////////////////
                                //cout<<"MZDAtaLoader data ready"<<","<<dataLine<<endl;
								///////////////////////////////////////////////
								//
								// Move the data from the temporary heap 
								// buffer dataLine, to the ScanGroup object.
								//
								// Without this double buffering type scheme,
								// reading the data from the file into memroy will take
								// a long time.
								//
								// Since we're using STL vectors, if anything goes
								// wrong such as an exception here,
								// 
								///////////////////////////////////////////////
								data.push_back(dataLine); // Load into temp buffer
								
								
							}  // end processing a non-empty line( a Pre-LC region )
						}  // end while loop
						
						///////////////////////////////////////////////////////
						// Add the data to the ScanGroup object
						///////////////////////////////////////////////////////
                        sg.addLevel1Data(data); // Move buffer to ScanGroup object
						
						inputfile.close();  // close the file after reading data
					} // end if block checking if the file is valid
                    else {

                        // file is not valid
                        return -6;

                    }

				} // end try block 
				catch (exception& e) {
					cout << e.what() << endl;
                    cout <<"loadLevel1File Exception while reading file"<<endl;
				}  // end trycatch block

            } // end read all Level 1 data lines

		} // end text mode processing
		else if ( type == "binary") {  // something to add later on

		} // end binary mode processing

	}  // end version "1.0.1"
	else {
		cerr<<"Invalid version indicator"<<endl;
		return -1; 
	}
	///////////////////////////////////////////////////////////////////////////
	///
	///
	/// Add additional version specific processing here.
	///
	///
	///////////////////////////////////////////////////////////////////////////

    cout<<"Finished loading level 1 file successfully"<<endl;

	return 0; // return successfully
}
	


} // end namespace
