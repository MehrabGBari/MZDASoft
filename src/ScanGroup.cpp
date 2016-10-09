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
/// ScanGroup.cpp
/// Contains all data for a set of overlapping scans
/// A scan group can be processed completely independently of all other
/// scan groups.
/// mzXML files are decomposed into groups of scans with overlap in the
/// mz and retention time dimension, this are called ScanGroups.
///
/// ScanGroups are further processed to generate level 1 and
/// level 2 files.
///
/// Level 1 files contain scangroups and
/// Level 2 files contain LC Candidates
///
///
/// Developer: Nelson.Ramirez@utsa.edu
/// In collaboration with: Dr. Michelle Zhang( Michelle.Zhang@usta.edu )
///
/// Portions of this code are ported algorithms from a Matlab source codes
/// developed by Dr. Michelle Zhang.
///
///////////////////////////////////////////////////////////////////////////////
#include "ScanGroup.h"
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <map>
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace msda {

ScanGroup::ScanGroup() : 
    mzData(0, std::vector<double>(0,0)),
    intensityData(0, std::vector<double>(0,0)),
    retentionTime(0,0),
    scanNum(0,0)	{
    ;  // Nothing to do within this constructor, except initialize variables
}

ScanGroup::~ScanGroup() {

    ///////////////////////////////////////////////////////////////////////////
    // Since all data are STL containers, therefore,
    // since they manage their own memory, we are not
    // responsible to clean their memory.
    ///////////////////////////////////////////////////////////////////////////
    ;
    //std::cout<<"Inside ScanGroup Destructor"<<std::endl;

}

///////////////////////////////////////////////////////////////////////////////
//
// Function to add new scans to the scan group object
//
///////////////////////////////////////////////////////////////////////////////
int ScanGroup::addScan(vector<double> & mz, vector<double> & intensity, double rt, int scan ) {

    ///////////////////////////////////////////////////////////////////////////
    //
    // For each new scan line, we use the vector data structure
    // to automatically allocate memory for each new scanline.
    //
    ///////////////////////////////////////////////////////////////////////////
    mzData.push_back(mz);  // copies the mz data
    intensityData.push_back(intensity); // copies the intensity data
    retentionTime.push_back(rt); // copies the retention time data
    scanNum.push_back(scan); // copies the scan number information

    ///////////////////////////////////////////////////////////////////////////
    // Debug code
    // Since we allow the user to select the level to process, and since
    // we are only processing level 1 information, we need to be clear about
    // the values in the rt and scan variables.
    // The mzData data structure
    ///////////////////////////////////////////////////////////////////////////
    //cout<<"rt="<<rt<<"scan="<<scan<<endl;
    return 0;
}


///////////////////////////////////////////////////////////////////////////////
///
/// Function to add the instrument data to the scan group object
///
///////////////////////////////////////////////////////////////////////////////
int ScanGroup::addInputDataFilename( string & filename ) {

    ///////////////////////////////////////////////////////////////////////////
    // We want to save the original input data filename within the
    // scan group object.
    // We want to provide flexibility to later on handle multiple files.
    ///////////////////////////////////////////////////////////////////////////
    inputDataFilename.push_back(filename);
    return 0; // return successfully
}



///////////////////////////////////////////////////////////////////////////////
///
/// Function to add the instrument data to the scan group object
/// "manufacturer" 	--> string --> instrumentInfoStrings
/// "model" 		--> string --> instrumentInfoStrings
/// "ionisation"   	--> string --> instrumentInfoStrings
/// "analyzer"      --> string --> instrumentInfoStrings
/// "detector"    	--> string --> instrumentInfoStrings
///////////////////////////////////////////////////////////////////////////////
int ScanGroup::addInstrumentData(	
        string & instrManufacturer,
        string & instrModel,
        string & instrIonisation,
        string & instrAnalyzer,
        string & instrDetector) {
    instrumentInfoStrings.insert(std::pair<string, string>("manufacturer", instrManufacturer));
    instrumentInfoStrings.insert(std::pair<string, string>("model", instrModel));
    instrumentInfoStrings.insert(std::pair<string, string>("ionisation", instrIonisation));
    instrumentInfoStrings.insert(std::pair<string, string>("analyzer", instrAnalyzer));
    instrumentInfoStrings.insert(std::pair<string, string>("detector", instrDetector));

    // debug, output the data
    map <string,string>::const_iterator end = instrumentInfoStrings.end();
    for (  map <string,string>::const_iterator it = instrumentInfoStrings.begin(); it != end; ++it) {
        std::cout << "instrument (string) data parameter: " << it->first;
        std::cout << "value: " << it->second << '\n';
    }

    return 0; // return successfully
}



///////////////////////////////////////////////////////////////////////////////
///
/// Function to add the instrument data to the scan group object
/// "scancount" 	--> int --> runHeaderInts
/// "endtime" 		--> double--> runHeaderDoubles
/// "starttime"   	--> double --> runHeaderDoubles
/// "endmz"      	--> double --> runHeaderDoubles
/// "highmz"    	--> double --> runHeaderDoubles
/// "lowmz" 		--> double --> runHeaderDoubles
/// "startmz"		--> double --> runHeaderDoubles
///
///////////////////////////////////////////////////////////////////////////////
int ScanGroup::addRunHeaderData(		
        int & ScanCount,
        double & dEndTime,
        double & dStartTime,
        double & endMZ,
        double & highMZ,
        double & lowMZ,
        double & startMZ) {
    runHeaderInts.insert(std::pair<string,int>("scancount",ScanCount));
    runHeaderDoubles.insert(std::pair<string,double>("endtime",dEndTime));
    runHeaderDoubles.insert(std::pair<string,double>("starttime",dStartTime));
    runHeaderDoubles.insert(std::pair<string,double>("endmz",endMZ));
    runHeaderDoubles.insert(std::pair<string,double>("highmz",highMZ));
    runHeaderDoubles.insert(std::pair<string,double>("lowmz",lowMZ));
    runHeaderDoubles.insert(std::pair<string,double>("startMZ",startMZ));

    // debug, output the integer data
    map <string,int>::const_iterator endInt = runHeaderInts.end();
    for (  map <string,int>::const_iterator itInt = runHeaderInts.begin(); itInt != endInt; ++itInt) {
        std::cout << "run header (int) data parameter: " << itInt->first;
        std::cout << "value: " << itInt->second << '\n';
    }

    // debug, output the double data
    map <string,double>::const_iterator endDouble = runHeaderDoubles.end();
    for (  map <string,double>::const_iterator itDouble = runHeaderDoubles.begin(); itDouble != endDouble; ++itDouble) {
        std::cout << "run header (double) data parameter: " << itDouble->first;
        std::cout << "value: " << itDouble->second << '\n';
    }


    return 0; // return successfully
}



///////////////////////////////////////////////////////////////////////////////
/// 
///
/// Function to add the scan header data to the scan group data structure
/// "acquisitionnum" 			--> int --> originalDataScanHeaderInts ( scan number in file(possible gaps)
/// "mergedscan"        		--> int --> originalDataScanHeaderInts ( only for ms level > 1, MS/MS data )
/// "mergedresultscannum		--> int --> originalDataScanHeaderInts  
/// "mergedresultstartscannum 	--> int --> originalDataScanHeaderInts
/// "mslevel"                   --> int --> originalDataScanHeaderInts ( level, ms = 1, ms/ms = 2 )
/// "numpossiblecharges"  		--> int --> originalDataScanHeaderInts
/// "peakscount"				--> int --> originalDataScanHeaderInts
/// "precursorcharge"			--> int --> originalDataScanHeaderInts ( mslevel >1 )
/// "precursorscannum"			--> int --> originalDataScanHeaderInts ( mslevel >1 )
/// "scanindex"					--> int --> originalDataScanHeaderInts ( A 1 based sequential index for non-sequential scan numbers )
/// "seqnum"					--> int --> originalDataScanHeaderInts  ( number in sequence observed file (1-based) )
/// "fileposition"				--> int --> originalDataScanHeaderInts ( ramp_fileoffset_t, must be 64 bits ) 
///
/// "basepeakintensity" 		--> double--> originalDataScanHeaderDoubles
/// "basepeakmz"   				--> double --> originalDataScanHeaderDoubles
/// "collisionenergy"	      	--> double --> originalDataScanHeaderDoubles
/// "highmz"    				--> double --> originalDataScanHeaderDoubles
/// "ionisationenergy" 			--> double --> originalDataScanHeaderDoubles
/// "lowmz"						--> double --> originalDataScanHeaderDoubles
/// "precursorintensity"		--> double --> originalDataScanHeaderDoubles ( mslevel > 1)
/// "compensationvoltage" 		--> double --> originalDataScanHeaderDoubles ( mslevel > 1)
/// "precursormz" 				--> double --> originalDataScanHeaderDoubles ( mslevel > 1)
/// "retentiontime"				--> double --> originalDataScanHeaderDoubles
/// "totalioncurrent"			--> double --> originalDataScanHeaderDoubles
/// 
/// "activationmethod" 			--> string --> originalDataScanHeaderStrings
/// "possiblecharges"			--> string --> originalDataScanHeaderStrings
/// "scantype"					--> string --> originalDataScanHeaderStrings
///
/// "possiblechargesarray"  	--> bool   --> originalDataScanHeaderBools
///
///////////////////////////////////////////////////////////////////////////////

int ScanGroup::addScanHeaderData(		
        int & acquisitionnum,
        int & mergedscan,
        int & mergedresultscannum,
        int & mergedresultstartscannum,
        int & mergedresultendscannum,
        int & mslevel,
        int & numpossiblecharges,
        int & peakscount,
        int & precursorcharge,
        int & precursorscannum,
        int & scanindex,
        int & seqnum,
        int & fileposition,
        double & basepeakintensity,
        double & basepeakmz,
        double & collisionenergy,
        double & highmz,
        double & ionisationenergy,
        double & lowmz,
        double & precursorintensity,
        double & compensationvoltage,
        double & precursormz,
        double & retentiontime,
        double & totalioncurrent,
        string & activationmethod,
        string & possiblecharges,
        string & scantype ) {

    ///////////////////////////////////////////////////////////////////////////
    // Create temporary data structures to hold the data prior
    // to adding to the vector data structure
    ///////////////////////////////////////////////////////////////////////////

    map<string,string> tempStrStrMap; // string scan header fields
    map<string, int> tempStrIntMap; // integer scan header fields
    map<string, double> tempStrDoubleMap; // double scan header fields
    map<string, bool> tempStrBoolMap; // boolean scan header fields

    ///////////////////////////////////////////////////////////////////////////
    //
    // Insert the scan header fields into the temporary data structures
    //
    ///////////////////////////////////////////////////////////////////////////
    tempStrIntMap.insert(std::pair<string,int>("acquisitionnum",acquisitionnum));
    tempStrIntMap.insert(std::pair<string,int>("mergedscan",mergedscan));
    tempStrIntMap.insert(std::pair<string,int>("mergedresultscannum",mergedresultscannum));
    tempStrIntMap.insert(std::pair<string,int>("mergedresultstartscannum",mergedresultstartscannum));
    tempStrIntMap.insert(std::pair<string,int>("mergedresultendscannum",mergedresultendscannum));
    tempStrIntMap.insert(std::pair<string,int>("mslevel",mslevel));
    tempStrIntMap.insert(std::pair<string,int>("numpossiblecharges",numpossiblecharges));
    tempStrIntMap.insert(std::pair<string,int>("peakscount",peakscount));
    tempStrIntMap.insert(std::pair<string,int>("precursorcharge",precursorcharge));
    tempStrIntMap.insert(std::pair<string,int>("precursorscannum",precursorscannum));
    tempStrIntMap.insert(std::pair<string,int>("scanindex",scanindex));
    tempStrIntMap.insert(std::pair<string,int>("seqnum",seqnum));
    tempStrIntMap.insert(std::pair<string,int>("fileposition",fileposition));
    originalDataScanHeaderInts.push_back(tempStrIntMap);


    ///////////////////////////////////////////////////////////////////////////
    //
    //Move the temporary data structures into the scan header objects
    //
    //
    ///////////////////////////////////////////////////////////////////////////


    return 0; // return successfully
}


///////////////////////////////////////////////////////////////////////////////
// 
//
// Convert PPM to MZ
//
// Mainly used for the following:
//
// For each centroid, we want to find a neighborhood mz range that corresponds
// to a certain parts per million range.
// At a given mz value, i.e. the centroid's mz value, and the specified
// ppm amount, i.e. at a centroid whose peak is at 1000 dalton, if we want
// to use a window of +/- 10 ppm, we would need to search a window
// of +/- 0.01 dalton.
// The bin range of this centroid would become [1000-0.01 to 1000+0.01]
//
//
// Reference:  M. Zhang, meeting on 11/1/11
///////////////////////////////////////////////////////////////////////////////
double ScanGroup::convertPPMtoDaltons(double mzValue, double ppmValue) {
    return (mzValue*ppmValue)/1000000;
}



///////////////////////////////////////////////////////////////////////////////
//
//
// Get PLROI( Pre-LC Region Of Interest Properties )
// Return information about data within the ScanGroup object.
//
// Added on 11/1/12 by Nelson.Ramirez@utsa.edu
///////////////////////////////////////////////////////////////////////////////
int ScanGroup::getNumSignalsPerPlcRoi(vector<int> & numSignalsPerPlcRoi){

    for (int i = 0; i < plcMzData.size(); i++ ) {
        numSignalsPerPlcRoi.push_back(plcMzData.at(i).size());
    }

    return 0;
}





///////////////////////////////////////////////////////////////////////////////
// Generate Pre-LC Regions of Interest(PLCROI)
//
// A Pre-LC Region of Interest can be defined as a region within a certain
// mz range as well as in a certain retention time ( scan range ) that
// meets certain criteria such that the set of mz,intensity values within
// this region can be used as a source data set to apply an lc
// candidate generation algorithm.
//
// A PLCROI represents the lowest level of data organization within the
// MZDASoft framework.
//
//
//
//
// return:
// -1: mzHistBin is empty
///////////////////////////////////////////////////////////////////////////////
int ScanGroup::findPreLCRegionsOfInterest(double ppmDelta) {

    ///////////////////////////////////////////////////////////////////////////
    //
    //  Have the signal processing class available
    //
    ///////////////////////////////////////////////////////////////////////////
    cbi::SignalProcessor spObject;  // signal processing object
    cbi::MathUtils muObject; // math utilities object

    ///////////////////////////////////////////////////////////////////////////
    //
    // We need to generate a dynamic histogram to generate a
    // dynamic set of bins.
    // No centroiding will be used.  All the non-zero intensity data
    // will be used.
    //
    // The bin generation can also be parallelized.
    // The mzBinHistStart and mzBinHistEnd vectors would need to have
    // their data pre-allocated.
    ///////////////////////////////////////////////////////////////////////////

    for (vector< vector<double> >::size_type i = 0; i < intensityData.size(); i++ ){
        for ( vector<double>::size_type j = 0; j < intensityData.at(i).size(); j++ ) {
            ///////////////////////////////////////////////////////////////////
            // Gather the range of data in mz and retention time
            // dimension.
            // The range of each centroid needs to be calculated dynamically
            // based on a ppm (parts per million ) amount.
            //
            //
            ///////////////////////////////////////////////////////////////////
            double mzDelta = convertPPMtoDaltons( mzData.at(i).at(j), ppmDelta);

            ///////////////////////////////////////////////////////////////////
            //
            // Each data items defines a new bin, with dynamically calculated
            // mz ranges based on the user specified PPM value on
            // each side of the centroid.
            //
            ///////////////////////////////////////////////////////////////////
            double mzMinRange = mzData.at(i).at(j)-mzDelta;
            double mzMaxRange = mzData.at(i).at(j)+mzDelta;


            ///////////////////////////////////////////////////////////////////
            // Make absolutely sure that the lower bound is never a negative
            // value and that upper bounds are never exceeded.
            // Checking for the upper bound is not necessary, since
            // having a bin range that is above the maximum mz value in the
            // data is not a problem.
            // However, we do not want to be utilizing negative values
            // for the minimum range of a bin.  Allowing negative
            // bin range values would complicate the range checking
            // process unnecessarily.
            ///////////////////////////////////////////////////////////////////
            if ( mzMinRange < 0 ) {mzMinRange = 0; }


            ///////////////////////////////////////////////////////////////////
            // Each centroid creates a dynamic bin range
            // according to its mz value and the user specified mz
            // values.
            ///////////////////////////////////////////////////////////////////
            mzHistBinStart.push_back(mzMinRange);
            mzHistBinEnd.push_back(mzMaxRange);


        }  // end loop through data in a single scan
    } // end loop through scans

    ///////////////////////////////////////////////////////////////////////////
    // Once the dynamic bins have been created, a consolidation
    // of the data needs to be performed.
    ///////////////////////////////////////////////////////////////////////////
#ifdef MSDADEBUG
#ifdef SCANGROUPDEBUG
#ifdef FINDPRELCCANDIDATESDEBUG
    cout<<"mzHistBinStart.size()"<<mzHistBinStart.size()<<endl;
#endif
#endif
#endif

    ///////////////////////////////////////////////////////////////////////////
    // Make sure that we do not proceed if mzHistBinStart.size() is zero
    ///////////////////////////////////////////////////////////////////////////
    if ( mzHistBinStart.size() < 1 ) {
        return -1; // mzHistBin is empty, no data to process
    }


    ///////////////////////////////////////////////////////////////////////////
    //
    // Sort the centroid bins according to starting mz value
    // This will be critical to enable efficient mapping of centroid
    // to all applicable bin ranges.
    // .. to be implemented........
    //
    // Use CBI library specialized Matlab -like sorting function, that
    // returns the newly sorted indices.
    // These indices are the re-organization of the centroids according
    // to their mz values.  After sorting, we no longer need to perform
    // a linear scan to find the groups.
    //
    // bin 0, [startvalue5:endvalue]
    // bin 1, [startvalue8:endvalue]
    // bin 2, [startvalue1:endvalue]
    // bin 3, [startvalue4:endvalue]
    // bin 4, [startvalue3:endvalue]
    // bin 5, [startvalue9:endvalue]
    //
    //
    //  -->
    //
    // bin 0, [moved from bin 2-->startvalue1:endvalue]
    // bin 1, [moved from bin 4-->startvalue3:endvalue]
    // bin 2, [moved from bin 3-->startvalue4:endvalue]
    // bin 3, [moved from bin 0-->startvalue5:endvalue]
    // bin 4, [startvalue8:endvalue]
    // bin 5, [startvalue9:endvalue]
    //
    // Then take each centroid, and find the closest matching set of bins
    // Then, expand the search on both sides of the closest matching bin until,
    // a bin that is outside of the range is found.
    //
    // The end result of this algorithm should be that all centroids
    // are mapped to ALL related bins.
    //
    // In addition, this will be really helpful for the XIC signal
    // processing, since the XIC signals ARE the list of centroids
    // located in these bins.
    //
    // An XIC signal is therefore a subset of the data contained in
    // the total set of bins.
    // XIC#1 --> bins1,2,3,4
    // XIC#2 --> bins3,4,5,6
    // XIC#3 --> bins....
    // ...
    // XIC#N --> bins....
    //
    //
    //
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // SORT DYNAMIC HISTOGRAM BINS
    // This step is absolutely critical for performance.
    //
    //
    // Generate a linear set of indices
    // 0,1,2,3,4,5,....total number of bins
    ///////////////////////////////////////////////////////////////////////////

    for (vector< double >::size_type i = 0; i< mzHistBinStart.size(); i++ ){
        mzHistBinSortIndex.push_back(i);
    }


    ///////////////////////////////////////////////////////////////////////////
    // Sort the bins, without changing the original data,
    // only returning the sorted order in the
    // mzHistBinSortIndex vector.
    // Sort ascending
    // Since the mzHistBinSortIndex originally contains a
    // list of all the integers from 0 to total number of bins-1,
    // After sortGetIndex function completes, we will have the ascending
    // sorted order of the bins.
    // The sorting is performed according to the mzHistBinStart
    // values.
    //
    // The sorting process can also be parallelized if needed.
    //
    ///////////////////////////////////////////////////////////////////////////
    spObject.sortGetIndex(mzHistBinStart, mzHistBinEnd, mzHistBinSortIndex, true);

    ///////////////////////////////////////////////////////////////////////////
    //
    // Parallel preLC Candidate generation
    // We need to generate cutpoints.  This can be done in parallel
    // with multiple cores.
    //
    // This is a candidate loop for OpenMP parallelism.
    //
    //
    //
    ///////////////////////////////////////////////////////////////////////////
    //int np=0;
    //int iam=0;
    //#pragma omp parallel default(shared) private(iam, np)
    //{
    //np = omp_get_num_threads();
    //iam = omp_get_thread_num();
    //cout<<"np="<<np<<"iam="<<iam<<endl;


    // Using a moving average and a specified ppm amount to determine
    // the cutpoints
    //for ( int k = 0; k < mzHistBinStart.size()-1; k++ ){
    int k = 0;
    preLcCandidateBoundaries.push_back(k);  // initialize first boundary to always be zero.


    for (vector< double >::size_type k = 0; k < mzHistBinStart.size()-1; k++ ){
        int xcurrent=0;
        int ycurrent=0;
        int xnext = 0;
        int ynext = 0;
        //int rCurrent = muObject.mapLinearIndexToCartesianCoordinates(mzHistBinSortIndex.at(k),mzData,xcurrent,ycurrent);
        //int rNext = muObject.mapLinearIndexToCartesianCoordinates(mzHistBinSortIndex.at(k+1),mzData,xnext,ynext);

        muObject.mapLinearIndexToCartesianCoordinates(mzHistBinSortIndex.at(k),mzData,xcurrent,ycurrent);
        muObject.mapLinearIndexToCartesianCoordinates(mzHistBinSortIndex.at(k+1),mzData,xnext,ynext);

        // add visualization for the histogram bin cut regions
        if ( mzData.at(xnext).at(ynext) > mzHistBinEnd.at(mzHistBinSortIndex.at(k))  ) {
            ///////////////////////////////////////////////////////////
            // The first region goes from 0 to 4
            // The second region goes from 5 to 553
            // In general, the boundaries are as follows:
            //
            // Boundary 0 = preLcCandidateBoundary.at(0) to preLcCandidateBoundary.at(1)-1
            // Boundary k = preLcCandidateBoundary.at(k) to preLcCandidateBoundary.at(k+1)-1, for k = 0 to preLcCandidateBoundary.size()-2
            //
            //preLcCandidateBoundary=0
            //preLcCandidateBoundary=5
            //preLcCandidateBoundary=553
            //preLcCandidateBoundary=558
            //preLcCandidateBoundary=566
            //
            ///////////////////////////////////////////////////////////

            preLcCandidateBoundaries.push_back(k+1);
            //cout<<"boundary at k="<<k<<",";
        }
        else {
            //cout<<"noboundary at k="<<k<<",";
        }
        //assert( r == 0 );
        //cout<<"mzHistBinSortIndex.at(k)="<<mzHistBinSortIndex.at(k)<<",next k+1="<<mzHistBinSortIndex.at(k+1)<<",mzdata="<<mzData.at(xcurrent).at(ycurrent)<<",next="<<mzData.at(xnext).at(ynext)<<",scancurrent="<<scanNum.at(xcurrent)<<",intensitycurrent="<<intensityData.at(xcurrent).at(ycurrent)<<endl;
    }
    preLcCandidateBoundaries.push_back(mzHistBinStart.size()-1);  // initialize last boundary to always be the last bin index


    //}  // This is the end of the OpenMP block.


    ///////////////////////////////////////////////////////////////////////////
    //
    //
    //  At this point the preLcCandidateBoundaries vector data structure
    //  contains the list of cut points within the
    //
    //
    //
    ///////////////////////////////////////////////////////////////////////////



    ///////////////////////////////////////////////////////////////////////////
    // Within a mz region, we want to generate the set of pre lc peak
    // candidates.
    // Within an mz region, we must separate the data according to
    // scan number.
    // The data within an mz region must be collapsed into a vector
    // along the scan dimension.
    // Further signal processing can be performed within this vector
    // to separate regions that are a bit too far in the scan dimension.
    //
    ///////////////////////////////////////////////////////////////////////////
    //cout<<"preLcCandidateBoundaries.size()="<<preLcCandidateBoundaries.size()<<endl;

    ///////////////////////////////////////////////////////////////////////////
    //
    //
    //
    // For each mz region, we must generate a vector so that we can view
    // the xic.
    //
    //
    //
    //
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // Generate the list of data items in the first boundary,
    // organizing the data according to scan number.
    // Boundary 0 --> [boundaryStart boundaryEnd]
    // Boundary 1 --> [boundaryStart boundaryEnd]
    // Boundary 2 --> [boundaryStart boundaryEnd]
    // Boundary 3 --> [boundaryStart boundaryEnd]
    // Boundary 4
    // Boundary 5
    // Boundary 6
    // Boundary.. N
    ///////////////////////////////////////////////////////////////////////////
    int boundaryStart = 0;
    int boundaryEnd = 0;
    int preLcCandidateCounter = 0; // Count the total number of preLcCandidates


    ///////////////////////////////////////////////////////////////////////////
    // i.e.
    //preLcCandidateBoundary=0
    //preLcCandidateBoundary=4
    //preLcCandidateBoundary=552
    //preLcCandidateBoundary=557
    //preLcCandidateBoundary=565
    //
    ///////////////////////////////////////////////////////////////////////////
    // The loop below is just for debugging purposes:
    //for ( int i = 0; i< preLcCandidateBoundaries.size(); i++ ) {
    //cout<<"preLcCandidateBoundary="<<preLcCandidateBoundaries.at(i)<<endl;
    //}


    for (vector<int>::size_type i = 1; i< preLcCandidateBoundaries.size()-1; i++ ){
        boundaryStart = preLcCandidateBoundaries.at(i);
        boundaryEnd = preLcCandidateBoundaries.at(i+1)-1;
        //cout<<"preLcCandidateBoundary="<<boundaryStart<<","<<boundaryEnd<<"diff="<<boundaryEnd-boundaryStart+1<<endl;

        int xcurrent=0;
        int ycurrent=0;
        ///////////////////////////////////////////////////////////////////////
        //
        // We need to perform the following on each pre-LC candidate
        //
        // 1. Noise estimation via histogramming (Thresholding )
        // 2. Smoothing of the pre-LC candidate signals
        // 3. Interval cutting
        //
        ///////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////
        // The memory for the following data structures is memory
        // that is recycled for each pre-lc candidate.
        // Since vector, map, multimap are STL containers, they take
        // care of managing their memory on the heap.
        // It is very important that these be NEW before addding data
        // to them, otherwise the aggregate data calculations will include
        // data from multiple pre lc candidates, instead of a single pre
        // lc candidate.
        ///////////////////////////////////////////////////////////////////////
        std::map<int,int> scanVsDataItemCountMap; // We must create a new map for each pre-lc candidate( important )
        std::map<int, double> scanVsIntensitySumMap;  // We must create a new map for each pre-lc candidate( important )
        std::map<int, double> scanVsAverageMzMap; // We must create a new map for each pre-lc candidate( important )
        std::multimap<int,double> scanMzValues; // We must create a new multimap for each pre-lc candidate for the original mz values (important)
        std::multimap<int,double> scanIntensityValues; // We must create a new multimap for each pre-lc candidate for the original intensity values (important)
        std::map<int,double> scanVsRetentionTimeMap; // We must create a new map for each pre-lc candidate ( very important, otherwise memory will be consumed )

        // For the data items inside of a lc candidate region ( that is inside of a boundary )
        for ( int j = boundaryStart; j <= boundaryEnd; j++ ){
            ///////////////////////////////////////////////////////////////////
            //int rCurrent = muObject.mapLinearIndexToCartesianCoordinates(mzHistBinSortIndex.at(j),mzData,xcurrent,ycurrent);
            ///////////////////////////////////////////////////////////////////
            muObject.mapLinearIndexToCartesianCoordinates(mzHistBinSortIndex.at(j),mzData,xcurrent,ycurrent);
            //cout<<"mz at region "<<i<<","<<mzData.at(xcurrent).at(ycurrent)<<","<<intensityData.at(xcurrent).at(ycurrent)<<","<<scanNum.at(xcurrent)<<endl;

            // We must aggregate the data according to scan number
            // Here is the word counting problem in a text document
            scanVsDataItemCountMap[scanNum.at(xcurrent)] = scanVsDataItemCountMap[scanNum.at(xcurrent)]++;
            scanVsIntensitySumMap[scanNum.at(xcurrent)] += intensityData.at(xcurrent).at(ycurrent);
            scanVsAverageMzMap[scanNum.at(xcurrent)] += mzData.at(xcurrent).at(ycurrent);
            scanMzValues.insert(pair<int, double>(scanNum.at(xcurrent), mzData.at(xcurrent).at(ycurrent)));
            scanIntensityValues.insert(pair<int,double>(scanNum.at(xcurrent), intensityData.at(xcurrent).at(ycurrent))  );

            scanVsRetentionTimeMap[scanNum.at(xcurrent)] = retentionTime.at(xcurrent);
        }

        ///////////////////////////////////////////////////////////////////////
        // Copying the local data structures to the global data
        // structures in the ScanGroup Class.
        ///////////////////////////////////////////////////////////////////////
        preLcCandidateScanVsItemCountMap.push_back(scanVsDataItemCountMap);
        preLcCandidateScanVsIntensitySumMap.push_back(scanVsIntensitySumMap);
        preLcCandidateScanVsAverageMzMap.push_back(scanVsAverageMzMap);
        preLcCandidateScanMzValues.push_back(scanMzValues);  // For each pre lc candidate, these are all the mz values at a given scan
        preLcCandidateScanIntensityValues.push_back(scanIntensityValues); // For each pre lc candidate, these are all the intensity values at a given scan
        preLcCandidateScanVsRetentionTimeMap.push_back(scanVsRetentionTimeMap); // For each pre lc candidates, this are a mapping of scan to retention time
        ///////////////////////////////////////////////////////////////////////
        // End Copying the local data structures to the global data
        // structures in the ScanGroup Class.
        ///////////////////////////////////////////////////////////////////////
        preLcCandidateCounter++;   // Increment the pre lc candidate counter


    }  // A new scanVsIntensitySumMap is created each iteration of i



    return 0;
}



///////////////////////////////////////////////////////////////////////////
///
/// MSDA Level 1 Format
/// 
/// Notes:
/// precision = 10 ( default to setprecision(10) )
///
/// Version "1.0.0": Initial version.
///
///////////////////////////////////////////////////////////////////////////
int ScanGroup::createMSDALevel1File( 
        string version,
        string filenamescript,
        string title,
        int scangroupid,
        int idnumwidth,
        int precision,
        int minScan,
        int maxScan,
        double minMz,
        double maxMz,
        double minRetentionTime,
        double maxRetentionTime,
        double mzBinningThreshold,
        int msLevel,
        int totalmsLevel1Scans ){

    ///////////////////////////////////////////////////////////////////////////
    // Generate the output filenames
    // according to the corresponding scangroup id and width as specified
    // by the user.  Multiple scangroups could be processed by a single
    // worker, so the fundamental identifier of a unit of msda work is the
    // scan group id.
    //
    ///////////////////////////////////////////////////////////////////////////
    string valString;
    stringstream scangroupidTemp (stringstream::in | stringstream::out);
    scangroupidTemp << setw(idnumwidth) << setfill('0');
    scangroupidTemp << scangroupid;
    scangroupidTemp >> valString;

    ///////////////////////////////////////////////////////////////////////////
    // Create the MZDA Level 1 file
    ///////////////////////////////////////////////////////////////////////////
    string mzdaLevelOnePathAndFilename = filenamescript;
    mzdaLevelOnePathAndFilename = mzdaLevelOnePathAndFilename.append(valString);
    mzdaLevelOnePathAndFilename = mzdaLevelOnePathAndFilename.append(".csv");

    ///////////////////////////////////////////////////////////////////////////
    //
    // Create the internal data structure information file, each with
    // its scan group identifier.
    //
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    //
    // Create the MZDA Level 1 file
    //
    ///////////////////////////////////////////////////////////////////////////
    ofstream mzdaLevelOneOutputStream;
    mzdaLevelOneOutputStream.open(mzdaLevelOnePathAndFilename.c_str());


    cout<<"ScanGroup::createMSDALevel1File,mzdaLevelOnePathAndFilename="<<mzdaLevelOnePathAndFilename<<endl;

    ///////////////////////////////////////////////////////////////////////////
    // Check the version of the format
    // Process version 1.0.0 MSDA Level 1 format
    ///////////////////////////////////////////////////////////////////////////
    if ( version == "1.0.0") {


        ///////////////////////////////////////////////////////////////////////
        //
        //
        //
        // Loop through the data for each pre-lc regions of interest
        // We need to allow for collaborative work on the processing of each
        // pre-lc candidate region for algorithm development purposes.
        // Dr. Zhang, will be working on an algorithm in Matlab for performing
        // separation of the pre-lc region data.  Once the algorithm is developed
        // it can be moved to C++.
        //
        // Given the differences in low-res to high-res data, the pre-processing
        // algorithm is highly data dependent and so developing a framework
        // to develop new algorithms for new data types is highly beneficial
        // for furthering the proteomics processing software infrastructure.
        // ( Mass Spectrometry Data Analysis Software )
        //
        //
        //
        ///////////////////////////////////////////////////////////////////////

        // Output the header data
        mzdaLevelOneOutputStream<<"mzdaversion="<<version<<","<<"precision="<<precision<<","<<"filenamescript="<<filenamescript<<","<<"datasource="<<title<<","<<"scangroupid="<<scangroupid<<","<<"minScan="<<minScan<<","<<"maxScan="<<maxScan<<","<<"minMz="<<minMz<<","<<"maxMz="<<maxMz<<","<<"minRetentionTime="<<minRetentionTime<<","<<"maxRetentionTime="<<maxRetentionTime<<","<<"mzBinningThresholdPPM="<<mzBinningThreshold<<","<<"msLevel="<<msLevel<<","<<"totalMsLevel1Scans="<<totalmsLevel1Scans<<",totalNumberOfPreLcRegions="<<preLcCandidateScanVsIntensitySumMap.size()<<",totalData_in_each_row=prelcregionid,scan#,retentiontime,intensitysum,averagemz,count,intensities,mzvalues:#"<<endl;

        // Loop through each pre lc candidate component
        for ( vector< map<int,double> >::size_type i = 0; i < preLcCandidateScanVsIntensitySumMap.size(); i++ ){
            // Loop through each pre lc candidate component
            for ( int j = minScan; j <= maxScan; j++ ){
                // only output those pre-lc candidate components that have at least 1 data element
                if ( preLcCandidateScanVsItemCountMap.at(i)[j] > 0 ) {
                    ///////////////////////////////////////////////////////////
                    // output precision set to 10
                    // The variable i contains the index into the preLcCandidateScanVsIntensitySumMap
                    // vector.  The index into these structures correspond to the
                    // prelc candidate region identifier
                    //
                    // The variable j contains the scan we need to obtain data for.
                    //
                    // The preLcCandidateScanVsIntensitySumMap data structure contains
                    // data for the intensity sum at a particular scan number for a given
                    // prelc candidate component.
                    //
                    // preLcCandidateScanVsAverageMzMap.at(i)[j]
                    //
                    ///////////////////////////////////////////////////////////
                    mzdaLevelOneOutputStream<<setprecision(precision)<<i<<","<<j<<","<<preLcCandidateScanVsRetentionTimeMap.at(i)[j]<<","<<preLcCandidateScanVsIntensitySumMap.at(i)[j]<<","<<preLcCandidateScanVsAverageMzMap.at(i)[j]/preLcCandidateScanVsItemCountMap.at(i)[j]<<","<<preLcCandidateScanVsItemCountMap.at(i)[j];
                    ///////////////////////////////////////////////////////////
                    //
                    //
                    // Output the list of mz values and intensity values
                    // corresponding to the current scan number j
                    //
                    //
                    // Per 01/04/12 Meeting with M.Zhang,
                    // We'll want to package C++ software for external building.
                    // -Add the min/max mz, min/max rt range to each
                    // pre-lc region. ( to avoid user of the format
                    // to have to re-scan the
                    // -Begin porting the Matlab code to convert.
                    // -Add version control to the MSDA Level format.
                    // -Create a dedicated MSDA Level v#.#.# reader in
                    // in Matlab, Python, C++.
                    // -Begin an MSDA Software Development Guide
                    //   - Git repo( A new one for msdasoftware )
                    //   - Loading
                    ///////////////////////////////////////////////////////////

                    // The vector of multi-maps facilitates the process of
                    // grouping the data according to scan information
                    pair<multimap<int,double>::iterator, multimap<int,double>::iterator> iiINT;
                    multimap<int,double>::iterator itINT; //Iterator to be used along with iiINT
                    iiINT = preLcCandidateScanIntensityValues.at(i).equal_range(j); //We get the first and last entry in iiINT

                    for(itINT = iiINT.first; itINT != iiINT.second; ++itINT)
                    {
                        mzdaLevelOneOutputStream<<","<<setprecision(precision)<<itINT->second;
                    }

                    // The vector of multi-maps facilitates the process of
                    // grouping the data according to scan information
                    pair<multimap<int,double>::iterator, multimap<int,double>::iterator> iiMZ;
                    multimap<int,double>::iterator itMZ; //Iterator to be used along with iiMZ
                    iiMZ = preLcCandidateScanMzValues.at(i).equal_range(j); //We get the first and last entry in iiMZ

                    for(itMZ = iiMZ.first; itMZ != iiMZ.second; ++itMZ)
                    {
                        mzdaLevelOneOutputStream<<","<<setprecision(precision)<<itMZ->second;
                    }

                    ///////////////////////////////////////////////////////////
                    // Separate data for a single PLCROI
                    ///////////////////////////////////////////////////////////
                    mzdaLevelOneOutputStream<<":"; // row separator
                }
            }
            if ( i < preLcCandidateScanVsIntensitySumMap.size()-1 ){
                mzdaLevelOneOutputStream<<"#\n";
            }
        }

        mzdaLevelOneOutputStream<<"#\n";

    } // End version "1.0.1" MSDA Level 1 format
    ///////////////////////////////////////////////////////////////////////////
    // End version 1.0.1 MSDA Level 1 format
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // Close the file
    ///////////////////////////////////////////////////////////////////////////
    mzdaLevelOneOutputStream.close();


    ///////////////////////////////////////////////////////////////////////////
    // Return successfully
    ///////////////////////////////////////////////////////////////////////////
    return 0;

}



///////////////////////////////////////////////////////////////////////////////
/// Return the number of scans within the ScanGroup
///////////////////////////////////////////////////////////////////////////////
int ScanGroup::size(void) {
    return( mzData.size() );
}


///////////////////////////////////////////////////////////////////////////////
/// Return the minimum retention time in the scan group
///////////////////////////////////////////////////////////////////////////////
double ScanGroup::findMinRetentionTime(void){
    cbi::SignalProcessor spObject;  // signal processing object
    double minRetentionTime=0;
    spObject.findMinValue(retentionTime,minRetentionTime);
    return minRetentionTime;

} //returns the minimum retention time value


///////////////////////////////////////////////////////////////////////////////
/// Return the maximum retention time in the scan group
///////////////////////////////////////////////////////////////////////////////
double ScanGroup::findMaxRetentionTime(void){
    cbi::SignalProcessor spObject;  // signal processing object
    double maxRetentionTime=0;
    spObject.findMaxValue(retentionTime,maxRetentionTime);
    return maxRetentionTime;

} //returns the maximum retention time value



///////////////////////////////////////////////////////////////////////////////
/// Populate buffer data structure when reading data in from an
/// MSDA Level 1 text format file
/// This method ised used from the MsDataLoader class.  The MsDataLoader
/// class takes care of loading the data into heap allocated temporary
/// memory using a vector STL data structure.  However, this buffer
/// needs to be transferred to the ScanGroup class for storage in the 
/// class where subsequent processing is to take place.
/// This way, we separate the code that interfaces with files, and the 
/// code that will process the data contained in those files.
///////////////////////////////////////////////////////////////////////////////
int ScanGroup::addLevel1Data(vector<string>& data){
    level1TextData = data;  // copies the level 1 data from the MsDataLoader
    // buffer to the ScanGroup data structures.
    // The memory alocated during the reading
    // of the file is automatically de-allocated
    // as soon as the buffer in the MsDataLoader
    // class goes out of scope.
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
//
// processLevel1Data
//
// This method looks at the msdaLevel1TextData data structure within the 
// ScanGroup class to check to see if it has been populated.  If it 
// has been populated, we want to process its contents into a 
// format that is easier to work with.
//
// Once the MSDA Level 1 data has been read into the buffer in this
// class.  We need to process the data in the buffer into a form that
// is accessible for further processing.
//
//  This is a port of the msdainterface.m function.
// 
//  It takes an MSDA level 1 format file and generates an in memory 
//  representation of the file for efficient generation of the
//  set of LC Candidates.
//
///////////////////////////////////////////////////////////////////////////////
int ScanGroup::processLevel1Data(void) {

    cout<<"Processing level 1 data"<<endl;

    ///////////////////////////////////////////////////////////////////////
    //
    // Parse the MSDA Level 1, v1.0.1 header data.
    // version ( type = string)
    // filename ( type = string )
    // datasource ( type = string  )
    // scangroupid ( type = int  )
    // minScan ( type = int )
    // maxScan ( type = int  )
    // minMz (
    // maxMz
    // minRetentionTime
    // maxRetentionTime
    // mzBinningThresholdPPM
    // msLevel
    // totalmsLevel1Scans
    // totalNumberOfPreLcRegions
    // row descritor string="totalData_in_each_row=prelcregionid,scan#,retenti
    //               ontime,intensitysum,averagemz,count,intensities,mzvalues:#"
    //   Each row contains the data for a single pre-LC region:
    //     prelcregionid ( type = int )
    //     scan# ( type = int )
    //     retentiontime ( type = double )
    //     intensitysum ( type = double )
    //     averagemz ( type = double )
    //     count ( type = int ), this is the number of
    //     intensities,.,.,.,mzvalues,.,.,
    //     : ( pre-lc region component separator )
    //     # ( pre-lc region separator )
    //
    ///////////////////////////////////////////////////////////////////////
    if( level1TextData.size() <= 1 ) {

        return -1; // insufficient data
    }

    // process header row, ( perhaps more than 1 header row in the future )
    for ( vector<string>::size_type i = 0; i < 1; i++ ){
        cout<<level1TextData.at(i)<<endl;


        // to be added


    }

    // depending on the total number of pre-lc candidates
    int numDataRows = level1TextData.size();



    for ( vector<string>::size_type i = 1; i < numDataRows; i++ ){
        // process each prelc region ofinterest
        string tempString = level1TextData.at(i);
        string tempStringComponent;

        // replace the "," with a blank
        replace(tempString.begin(),tempString.end(),',',' ');
        // replace the "#" with a blank
        //replace(tempString.begin(),tempString.end(),'#',' ');

        // copy the space delimited data to a stringstream
        istringstream tempStringStream(tempString);

        int componentCounter = 0;  // keep track of the components in an lc region

        // get each component of a particular pre-lc candidate

        ///////////////////////////////////////////////////////////////////
        //
        //
        // For each pre-LC Region we need to create a temporary storage
        // location
        //
        ///////////////////////////////////////////////////////////////////
        vector<int> plcTempRegionID;  // the list of plc region IDs for each component within a pre-LC region
        vector<int> plcTempScanNumber; // the list of scan numbers for each component within a pre-LC region
        vector<double> plcTempRetentionTime; // the list of retention times for each component within a pre-LC region
        vector<int> plcTempCount;  // this will tell us how many data points in a component
        vector< vector<double> > plcTempMzData; // this contains the actual data for all the components within a single preLC region
        vector< vector<double> > plcTempIntensityData; // this contains the actual data for all the components within a single preLC region


        ///////////////////////////////////////////////////////////////
        // Loop through all the components of a single pre-LC region of
        // interest.
        ///////////////////////////////////////////////////////////////
        while( getline(tempStringStream,tempStringComponent,':') ) {

            //cout<<"Region="<<i<<"Component="<<componentCounter<<endl;
            //cout<<tempStringComponent<<endl;

            componentCounter++;

            ///////////////////////////////////////////////////////////////
            // Add the component data to a set of temporary
            // variables.  Note that the tempMzData and
            // tempIntensityData vectors get
            ///////////////////////////////////////////////////////////////
            int tempRegionID=0;
            int tempScanNumber=0;
            double tempRetentionTime=0;
            double tempIntensitySum =0;
            double tempAverageMz=0;
            int tempCount=0;  // this will tell us how many data points in a component
            vector<double> tempMzData;
            vector<double> tempIntensityData;
            // skip components starting with '#'.
            // In a newer version of the MSDA Level 1 format,
            // we may want to make sure that this check is not needed.
            if( tempStringComponent[0] == '#' ){
                continue;   // skip this component
            }
            ///////////////////////////////////////////////////////////////
            // Create a string stream so we can parse the data
            // within a single component
            ///////////////////////////////////////////////////////////////
            stringstream tempComponentStringStream(tempStringComponent);


            ///////////////////////////////////////////////////////////////
            // Create a string stream so we can parse the data
            // within a single component.
            // We first need to get the header data from each component
            // The header for each component contains the following
            // information:
            // region identifier, scan number, retention time
            // intensity sum, averagemz, data item count
            ///////////////////////////////////////////////////////////////
            tempComponentStringStream>>tempRegionID>>tempScanNumber>>tempRetentionTime>>tempIntensitySum>>tempAverageMz>>tempCount;

            //cout<<setprecision(10)<<"tempRegionId="<<tempRegionID<<"tempScanNumber="<<tempScanNumber
            //   <<"tempRetentionTime="<<tempRetentionTime<<"tempIntensitySum="
            //  <<tempIntensitySum<<"tempAverageMz="<<tempAverageMz<<"tempCount="<<tempCount<<",";

            ///////////////////////////////////////////////////////////////
            // Load the mz and intensity data according to the
            // count specified in the component's header.
            ///////////////////////////////////////////////////////////////

            ///////////////////////////////////////////////////////////////
            // Get the intensity data
            ///////////////////////////////////////////////////////////////
            double tempIntensityDataBuffer = 0;
            for ( vector<string>::size_type j = 0; j < tempCount; j++ ){
                tempIntensityDataBuffer = 0;
                tempComponentStringStream>>tempIntensityDataBuffer;
                tempIntensityData.push_back(tempIntensityDataBuffer);

                // cout<<setprecision(10)<<tempIntensityDataBuffer<<" ";

            } // end loop through mz,intensity data for a single component


            //cout<<endl;

            ///////////////////////////////////////////////////////////////
            // Get the mz data
            ///////////////////////////////////////////////////////////////
            double tempMzDataBuffer = 0;
            for ( vector<string>::size_type j = 0; j < tempCount; j++ ){
                tempMzDataBuffer = 0;  // clear temporary holding var
                tempComponentStringStream>>tempMzDataBuffer; // use stream as output into buffer
                tempMzData.push_back(tempMzDataBuffer); // push the data onto the

                //cout<<setprecision(10)<<tempMzDataBuffer<<" ";

            } // end loop through mz,intensity data for a single component


            //cout<<endl;



            ///////////////////////////////////////////////////////////////////
            //
            // Add the component's temporary data to the ScanGroup
            // object's MSDA Level 1 in memory storage data structures.
            //
            //
            // The in-memory data structures will be used for all subsequent
            // processing stages, namely LC Candidate detection as well
            // as Peptide Candidate generation.
            //
            //
            //
            ///////////////////////////////////////////////////////////////////
            plcTempRegionID.push_back(tempRegionID);  // the list of plc region IDs for each component within a pre-LC region
            plcTempScanNumber.push_back(tempScanNumber); // the list of scan numbers for each component within a pre-LC region
            plcTempRetentionTime.push_back(tempRetentionTime); // the list of retention times for each component within a pre-LC region
            plcTempCount.push_back(tempCount);  // this will tell us how many data points in a component
            plcTempMzData.push_back(tempMzData); // this contains the actual data for all the components within a single preLC region
            plcTempIntensityData.push_back(tempIntensityData); // this contains the actual data for all the components within a single preLC region



        } // end while loop through components of a single region


        plcRegionID.push_back(plcTempRegionID);  // the list of plc region IDs for each component within a pre-LC region
        plcScanNumber.push_back(plcTempScanNumber); // the list of scan numbers for each component within a pre-LC region
        plcRetentionTime.push_back(plcTempRetentionTime); // the list of retention times for each component within a pre-LC region
        plcCount.push_back(plcTempCount);  // this will tell us how many data points in a component
        plcMzData.push_back(plcTempMzData); // this contains the actual data for all the components within a single preLC region
        plcIntensityData.push_back(plcTempIntensityData); // this contains the actual data for all the components within a single preLC region







    } // end loop through each pre lc region of interest


    return 0;
}  // end processLevel1Data function


///////////////////////////////////////////////////////////////////////////////
//
// This is a peak detection algorithm.  This is a port of the
// SimplestVersionWorking020513.m file located in the docs directory.
//
// This algorithm operates on a PLCROI.   This is the fundamental work unit
//
// Key parameters:
// mzWindowPPM = 5;   ( Used to calculate the window width )
// numSmoothPoint = 3; ( Used as the width of the smoothing filter )
// minLCLength = 3; ( Used by the getInterval and splitInterval methods )
// minMZLength = 3;  ( Used after the intensity trimming and the shape trimming )
// noiseThresholdLevel=1; ( Not Used: later it may be used in getNoiseThreshold and getInterval)
// massResolution=60000; ( Used to calculate peakWidthStd( width around peak), which is further used to generate mzGridIDs )
// R2Threshold = 0.90; ( Not Used )
// LCPeakApexTolerance=1;  ( Used for setting scanIDlow and scanIDhigh )
// rawMZCorrThreshold=0.8; ( This is used by the mz peak shape selection logic )
//
//
///////////////////////////////////////////////////////////////////////////////
int ScanGroup::lccAlgorithm(string sourceFilename,
                            int scanGroupId,
                            int plcRoi,
                            int mzWindowPPM,
                            int numSmoothPoint,
                            int minLCLength,
                            int minMZLength,
                            double noiseThresholdLevel,
                            double massResolution,
                            int LCPeakApexTolerance,
                            double rawMZCorrThreshold,
                            string outputFilename ) {

    return 0;
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//////////////////////Version 2 of the LCC Algorithm///////////////////////////
//////////////////////Feb 26, 2013/////////////////////////////////////////////
//////////////////////Key improvements:////////////////////////////////////////
//////////////////////- Resampling the retention time axis/////////////////////
//////////////////////- No longer using XICs, using resampledRegion////////////
//////////////////////  as the XICs.///////////////////////////////////////////
//////////////////////- Added split interval in generating LCPeakApexMap///////
//////////////////////- Added an estimated elution profile algorithm///////////
//////////////////////  to summarize the data contained within/////////////////
//////////////////////  each peak into a single signal.////////////////////////
//////////////////////- Variable renaming to reflect the use of////////////////
//////////////////////  retention time instead of scan number./////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
// This is a peak detection algorithm.  This is a port of the
// genlccandidatesSimplestVersionWorking022613originalemail.m
// file located in the docs directory.
//
// This algorithm operates on a PLCROI.   This is the fundamental work unit
// for the peak detection algorithms.
//
// Key parameters:
// numSmoothPoint = 3; ( Used as the width of the smoothing filter )
// minLCLength = 3; ( Used by the getInterval and splitInterval methods )
// minMZLength = 3;  ( Used after the intensity trimming and the shape trimming )
// noiseThresholdLevel=20000; ( Used as a key control for the total number of peaks detected )
// massResolution=60000; ( Used to calculate peakWidthStd( width around peak), which is further used to generate mzGridIDs )
// LCPeakApexTolerance=2;  ( Used for setting retention time peak ranges )
// rawMZCorrThreshold=0.8; ( This is used by the mz peak shape selection logic )
//
//
///////////////////////////////////////////////////////////////////////////////
int ScanGroup::lccAlgorithmVersion2(string sourceFilename,
                            int scanGroupId,
                            int plcRoi,
                            int numSmoothPoint,
                            int minLCLength,
                            int minMZLength,
                            double noiseThresholdLevel,
                            double massResolution,
                            int LCPeakApexTolerance,
                            double rawMZCorrThreshold,
                            string outputFilename,
                            int logLevel) {


    cout<<"lccAlgorithm Version2, 030513"<<endl;

    cbi::SignalProcessor sgObject;    // Signal Processing object
    msda::Utilities utilitiesObject;  // Utilities object
 \
    int plcNum = plcRegionID.size();  // get the total number of preLC regions of interest in this scan group


    if ( plcRoi >= plcRegionID.size() ) {
        cout<<"invalid plcRoi requested, exceeds number of available PLCROI"<<endl;
        return -1;
    }

    // Get the number of components within the current PLCROI
    int numScans = plcRegionID.at(plcRoi).size();

    ///////////////////////////////////////////////////////////////////////////
    // originalScanNums maps to plcScanNum.at(plcRoi)
    // rtaxis maps to plcRetentionTime.at(plcRoi)
    //
    // The following is a set of updated Matlab code:
    //  MATLAB SOURCE:
    //  originalScanNums=preLCRegions.scannum{preLCRegionid}; % add by MZ 02/05/13 --> maps to plcScanNum.at(plcRoi)
    //  rtaxis = preLCRegions.retentiontime{preLCRegionid};  % moved here by MZ 02/05/13 --> maps to plcRetentionTime.at(plcRoi)
    //  rtdiff=getDiff(rtaxis,numScans,-1); %add by MZ 02/26/13
    //  rtDelta=min(rtdiff);%add by MZ 02/26/13
    //  rtMin=min(rtaxis);%add by MZ 02/26/13
    //  rtMax=max(rtaxis);%add by MZ 02/26/13
    //  rtGrid=rtMin:rtDelta:rtMax;%add by MZ 02/26/13
    //  rtGridSize=length(rtGrid);%add by MZ 02/26/13
    //
    //
    ///////////////////////////////////////////////////////////////////////////
    vector<double> rtAxis;      //added 03/05/13
    rtAxis = plcRetentionTime.at(plcRoi); //added 03/05/13

    vector<int> originalScanNums;   // added 03/05/13
    originalScanNums = plcScanNumber.at(plcRoi);    // added 03/05/13

    vector<double> rtDiff;      // added 03/05/13
    int rc = utilitiesObject.getDiff(rtAxis, -1, rtDiff); // added 03/05/13

    // Next we need to calculate the minimum retention time delta,
    // this will help us determine how fine a grid we need to use.
    double rtDelta = 0;
    sgObject.findMinValue(rtDiff, rtDelta);


    // find minimum retention time of the PLCROI
    double rtMin=-1.0; // use this as a verification mechanism
    // no retention time value should ever be negative
    sgObject.findMinValue(rtAxis, rtMin );


    // find maximum retention time of the PLCROI
    double rtMax=-1.0; // use this as a verification mechanism
    // no retention time value should ever be negative
    sgObject.findMaxValue(rtAxis, rtMax );


    ///////////////////////////////////////////////////////////////////////////
    //
    // Generate the retention time grid
    // rtGrid=rtMin:rtDelta:rtMax;%add by MZ 02/26/13
    // Note: The following approach is taken to port the Matlab auto-sequencing
    // feature.
    //
    ///////////////////////////////////////////////////////////////////////////
    vector<double> rtGrid; // we need to make sure these are the same length as Matlab variable
    double rtIndex=rtMin;
    for( rtIndex = rtMin; rtIndex <= (rtMax); rtIndex= rtIndex+rtDelta ) {
        rtGrid.push_back(rtIndex);
    }
    vector<double>::size_type rtGridSize = rtGrid.size();  // rt dimension length


    // These scan number values map back to the original file.
    vector<int> scanNumMinAllComponents; // temporary buffer for scan numbers
    vector<int> scanNumMaxAllComponents; // temporary buffer for scan numbers

    vector<double> retentionTimeMinAllComponents; // temporary buffer for retention times
    vector<double> retentionTimeMaxAllComponents; // temporary buffer for retention times

    vector<double> mzMinAllComponents; // temporary buffer for mz values
    vector<double> mzMaxAllComponents; // temporary buffer for mz values

    // we need to keep track of the lengths of each component of a prelc region of interest
    vector<int> mzLength; // temporary buffer for mz length values

    ///////////////////////////////////////////////////////////////////////////
    // Loop through all the components of a single pre lc region of interest
    // Keeping track of minimum and maximum value information for each
    // component, so that we can make decisions based on mz data
    //of all the components within a PLCROI.
    ///////////////////////////////////////////////////////////////////////////
    for ( vector<string>::size_type j = 0; j <  numScans ; j++ ){

        // find minimum scan number of each component
        int minScanNumTemp=-1; // use this as a verification mechanism
        // no scan number value should ever be negative
        sgObject.findMinValue(plcScanNumber.at(plcRoi), minScanNumTemp );
        scanNumMinAllComponents.push_back(minScanNumTemp);

        // find maximum scan number of each component
        int maxScanNumTemp=-1;// use this as a verification mechanism
        // no scan number value should ever be negative
        sgObject.findMaxValue(plcScanNumber.at(plcRoi), maxScanNumTemp );
        scanNumMaxAllComponents.push_back(maxScanNumTemp);

        // find minimum retention time of each component
        double minRetentionTimeTemp=-1.0; // use this as a verification mechanism
        // no retention time value should ever be negative
        sgObject.findMinValue(plcRetentionTime.at(plcRoi), minRetentionTimeTemp );
        retentionTimeMinAllComponents.push_back(minRetentionTimeTemp);

        // find maximum retention time of each component
        double maxRetentionTimeTemp=-1.0; // use this as a verification mechanism
        // no retention time value should ever be negative
        sgObject.findMaxValue(plcRetentionTime.at(plcRoi), maxRetentionTimeTemp );
        retentionTimeMaxAllComponents.push_back(maxRetentionTimeTemp);


        // finding the minimum mz values for this component
        double minMzTemp=-1.0; // use this as a verification mechanism
        // no retention time value should ever be negative
        sgObject.findMinValue(plcMzData.at(plcRoi).at(j), minMzTemp );
        mzMinAllComponents.push_back(minMzTemp);

        // finding the maximum mz values for this component
        double maxMzTemp=-1.0; // use this as a verification mechanism
        // no retention time value should ever be negative
        sgObject.findMaxValue(plcMzData.at(plcRoi).at(j), maxMzTemp );
        mzMaxAllComponents.push_back(maxMzTemp);


        // finding the mzlength values for this component
        int mzLengthTemp=plcMzData.at(plcRoi).at(j).size(); // use this as a verification mechanism
        // no retention time value should ever be negative
        mzLength.push_back(mzLengthTemp);


        // sanity check
        // if these are ever negative, something went wrong, we do not
        // want to continue past this error.
        if ( minScanNumTemp < 0 || maxScanNumTemp < 0 || minRetentionTimeTemp < 0 || maxRetentionTimeTemp < 0  || minMzTemp < 0 || maxMzTemp < 0) {
            cerr<<"Error in finding min/max scannum or retention time value , mz value of each prelc region of interest componet at preLCROI = "<<plcRoi<<", Component= "<<j<<endl;
            return -2;
        }  // end if check

    }  // end loop through all components of a single pre lc region of interest

    ///////////////////////////////////////////////////////////////////////
    // mzGrid is a grid based on the min and max value for all
    // components in a given pre LC region of interest
    // Create an equal spaced box, where the length of each grid
    // column is calculated as follows:
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    // Find the minimum and maximum scan# and retention time
    // values for all components of the current preLC region of interest
    ///////////////////////////////////////////////////////////////////////
    int preLCROIMinScanNumTemp = -1; // use this for global verification
    int preLCROIMaxScanNumTemp = -1; // use this for global verification
    double preLCROIMinRetentionTimeTemp = -1.0; // use this for global verification
    double preLCROIMaxRetentionTimeTemp = -1.0; // use this for global verification
    double preLCROIMinMzTemp = -1.0; // use this for global verification
    double preLCROIMaxMzTemp = -1.0; // use this for global verification


    int preLCROIMinMzLengthTemp = -1.0; // use this for global verification
    int preLCROIMaxMzLengthTemp = -1.0; // use this for global verification
    sgObject.findMinValue(scanNumMinAllComponents, preLCROIMinScanNumTemp );
    sgObject.findMaxValue(scanNumMaxAllComponents, preLCROIMaxScanNumTemp );
    sgObject.findMinValue(retentionTimeMinAllComponents, preLCROIMinRetentionTimeTemp );
    sgObject.findMaxValue(retentionTimeMaxAllComponents, preLCROIMaxRetentionTimeTemp );
    sgObject.findMinValue(mzMinAllComponents, preLCROIMinMzTemp );
    sgObject.findMaxValue(mzMaxAllComponents, preLCROIMaxMzTemp );
    sgObject.findMinValue(mzLength, preLCROIMinMzLengthTemp );
    sgObject.findMaxValue(mzLength, preLCROIMaxMzLengthTemp );

    double mzBegin=preLCROIMinMzTemp;  // Find the component with the minimum mz
    double mzLast=preLCROIMaxMzTemp;    // Find the component with the maximum mz
    ///////////////////////////////////////////////////////////////////////////
    // Key algorithmic update in the generation of the mz grid,
    // By basing the mzGridWidth on the meanMz, we generate
    // a dynamic mzGrid according to the current mz range.  This
    // is especially important since a finer grid is needed as the value
    // of mz increases.
    ///////////////////////////////////////////////////////////////////////////
    double meanMz = (mzBegin+mzLast)/2; // New, algorithmic update by MZ on 04/16/13
    double sqLog4 = 1.177410022515475; // value of sqrt(log(4)) ( avoid dependencies on sqrt,log functions
    double mzTemp1 = sqLog4*massResolution;
    double mzTemp2 = meanMz*3;
    double mzTemp3 = mzTemp2/mzTemp1;
    double mzGridWidth = mzTemp3/10;

    ///////////////////////////////////////////////////////////////////////////
    // The previous mzGridwidth calculation approach has been removed, since
    // it does not change well to different mz values.
    ///////////////////////////////////////////////////////////////////////////
    //int maxMzLength = preLCROIMaxMzLengthTemp; // maximum number of data items in the current mz vector for all components
    //double mzGridWidth=(mzLast-mzBegin)/(maxMzLength)/2;  // fixed bug on 10/16/12(nrz)


    vector<double> mzGrid; // bug, 10/16/12( Matlab, length is 1 more than C++ length )
    double mzIndex=mzBegin;
    ///////////////////////////////////////////////////////////////////////////
    // bug, this is a fix to ensure we match Matlab version of ( mzBegin:mzGridWidth:mzLast)
    // There will likely be cases in which there is a discrepancy
    // We know there are bound to be discrepancies at the right edge of the grid
    // This is an issue with reproducibility of the following
    // Matlab built-in operation:
    // d1:d2:d3 where d1,d2,d3 are doubles
    // Matlab makes minute adjustments to the interval... algorithm they
    // are using is not known.
    // We need to create the mzGrid in a way that avoids the
    // auto-sequence generation in Matlab, otherwise we will
    // have data-dependent discrepancies at a variety of places.
    // Since the mzGrid is a critical component, unless we can ensure the
    // mzGrid is generated in a exactly reproducible manner, downstream
    // processing will diverge between Matlab/C++ versions.
    //
    // Potential solution:
    // Remove the use of the autosequencing feature in the Matlab version
    // for all sequences that are not Integer sequences.
    // Changing
    // mzIndex = mzBegin; mzIndex <= (mzLast); mzIndex= mzIndex+mzGridWidth
    // to
    // mzIndex = mzBegin; mzIndex <= (mzLast+mzGridWidth); mzIndex= mzIndex+mzGridWidth
    //
    // resolves the issue in the particular Level 1, plcroid being
    // used for development.
    //
    ///////////////////////////////////////////////////////////////////////////
    for( mzIndex = mzBegin; mzIndex <= (mzLast); mzIndex= mzIndex+mzGridWidth ) {
        mzGrid.push_back(mzIndex);
    }
    vector<double>::size_type mzGridLength = mzGrid.size();  // mz dimension length

    // validation up to this point ok... 10/16/12@7:20p.m.
    vector< vector<double> > resampledRegion(numScans,vector<double>(mzGridLength,0) );

    // in order to resample along the rt dimension, we need to keep
    // a copy of the transposed resampledRegion so we can
    // directly get access to the retention time signal.
    // added on 03/11/13
    vector< vector<double> > resampledRegionTranspose(mzGridLength, vector<double>(numScans,0)  );

    ///////////////////////////////////////////////////////////////////////////
    //
    // Next we need to populate the resampledRegion data structure:
    //
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////
    // Loop through all the component signal of a pre-LC
    // region of interest.
    ///////////////////////////////////////////////////////////////
    for( vector<double>::size_type scanid = 0; scanid < numScans; scanid++) {
        int rc = 0;
        vector<double> resampledVectorOutput;
        vector<double> resampledVectorOutputDefault(mzGrid.size(),0);

        // Note:
        // mzVector --> maps to plcMzData.at(plcRoi).at(scanid)
        // intensityVector --> maps to plcIntensityData.at(plcRoi).at(scanid)
        rc = utilitiesObject.resampledVector(plcMzData.at(plcRoi).at(scanid),plcIntensityData.at(plcRoi).at(scanid),mzGrid,resampledVectorOutput);
        if ( rc != 0 ){
            // bug found on 07/24/12: when there are too few
            // data items in the mz and intensity vectors, that is
            // less than or equal to 2 items, the resampling
            // process cannot function, and so we
            // should discard this component.
            resampledRegion.at(scanid) = resampledVectorOutputDefault;

            continue;
        }

        assert( mzGrid.size() == resampledVectorOutput.size() );
        ///////////////////////////////////////////////////////////
        // Update the resampledRegion data structure
        // Bug fix: 10/25/12.  We were not saving the
        // resampledVectorOutput into the resampledRegion
        // data structure.
        ///////////////////////////////////////////////////////////
        resampledRegion.at(scanid) = resampledVectorOutput; // added on 10/25/12

    } // end for loop


    ///////////////////////////////////////////////////////////////////////
    //
    // We need to store the data generated by the resampling
    // along the mz dimension within a data structure that
    // will provide direct access to vectors along the retention
    // time dimension, so that we can then resample these vectors
    // along the retention time dimension.
    //
    ///////////////////////////////////////////////////////////////////////
    for( vector<double>::size_type scanid = 0; scanid < numScans; scanid++) {
        for ( vector<double>::size_type mzid = 0 ; mzid < resampledRegion.at(scanid).size(); mzid++ ) {
            double dataItem = resampledRegion.at(scanid).at(mzid);
            resampledRegionTranspose.at(mzid).at(scanid)=dataItem;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // Next we need to generate an XIC data structure, however, this
    // data structure has been changed from being a resampled data
    // structure in the mz dimension to one resampled along the retention
    // time dimension.
    //
    // We need to be careful about the storage format of the updated
    // XICs data structure.
    //
    // We need fast access to all the resampled retention time signal
    // data, not the resampled mz data.  Therefore, we must store
    // complete vectors of resampled retention time signals.
    //
    // So that we can access the XICs data structure as follows:
    // XICs.at(mzgridindex) --> <resampled retention time signal vector>
    //
    ///////////////////////////////////////////////////////////////////////////
    vector< vector<double> > XICs(mzGridLength, vector<double>(rtGridSize,0));
    vector< vector<double> > smoothXICs(mzGridLength, vector<double>(rtGridSize,0));
    vector< vector<double> > resampledRegionXICs(rtGridSize, vector<double>(mzGridLength,0));

    for( vector<double>::size_type mzGridId = 0; mzGridId < XICs.size(); mzGridId++) {
        int rc = 0;
        vector<double> resampledVectorOutput;
        vector<double> resampledVectorOutputDefault(rtGrid.size(),0);
        // retention time axis, intensity data, retention time resampled grid, output resampled retention time vector
        rc = utilitiesObject.resampledVector(rtAxis,resampledRegionTranspose.at(mzGridId),rtGrid,resampledVectorOutput);
        if ( rc != 0 ){
            // bug found on 07/24/12: when there are too few
            // data items in the mz and intensity vectors, that is
            // less than or equal to 2 items, the resampling
            // process cannot function, and so we
            // should discard this component.
            XICs.at(mzGridId)=resampledVectorOutputDefault;
            continue;
        }
        assert( rtGrid.size() == resampledVectorOutput.size() );

        XICs.at(mzGridId) = resampledVectorOutput;
    }

    ///////////////////////////////////////////////////////////////////////////
    //
    // At this point we need to consolidate the PLCROI resampled data
    // into 2 fundamental data structures:
    //   KEY DATA STRUCTURES CONTAINING RESAMPLED
    //   PLCROI DATA:
    //
    //   XICs.at(mzGridId) = resampled retention time signal
    //   resampledRegionXICs.at(rtGridId) = resampled mz signal ( transpose of XICs )
    //
    // Fundamentally, we need fast access to both dimensions of the
    // resampled data.
    //  The XICs data structure provides access to vectors of resampled retention
    //  time data,
    //   and
    //  The resampledRegionXICs data structure provides access to resampled
    //  mz data.
    //
    // Both data structures must contain the same exact data,
    // the key difference is the way the data is organized in each case.
    //
    //
    //
    ///////////////////////////////////////////////////////////////////////////

    // Populate the resampledRegionTranspose data structure
    // This is equal to the transpose of the XICs data structure
    for( vector<double>::size_type idy = 0; idy < XICs.size(); idy++) {
        for ( vector<double>::size_type idx = 0 ; idx < XICs.at(idy).size(); idx++ ) {
            double dataItem = XICs.at(idy).at(idx);
            resampledRegionXICs.at(idx).at(idy) = dataItem;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    //
    // Loop through all the XIC windows
    //
    ///////////////////////////////////////////////////////////////////////////
    rc = 0; // hold return codes
    for ( int i = 0; i < XICs.size(); i++ ) {
        // looping through each window

        vector<double> tempSmoothedXIC(XICs.at(i).size(),0); // temporary storage to store smoothed signal

        ///////////////////////////////////////////////////////////////////////
        // apply smoothing 3 times.
        ///////////////////////////////////////////////////////////////////////
        rc = utilitiesObject.movingAvgSmooth(XICs.at(i), numSmoothPoint, tempSmoothedXIC);
        rc = utilitiesObject.movingAvgSmooth(tempSmoothedXIC, numSmoothPoint, tempSmoothedXIC);
        rc = utilitiesObject.movingAvgSmooth(tempSmoothedXIC, numSmoothPoint, tempSmoothedXIC);


        ///////////////////////////////////////////////////////////////////////
        // temporary output for debugging the smoothing algorithm
        ///////////////////////////////////////////////////////////////////////
        smoothXICs.at(i) = tempSmoothedXIC;

    } // end loop through all the XIC windows

    ///////////////////////////////////////////////////////////////////////////
    //
    // numWindows is the dimension of the mz data
    //
    ///////////////////////////////////////////////////////////////////////////
    int numWindows = XICs.size();  // This maps to numWindows = mzGridLength


    ///////////////////////////////////////////////////////////////////////////
    // create preLCPeakApex and LCPeakApex data structures
    ///////////////////////////////////////////////////////////////////////////

    vector< vector<bool> > preLCPeakApexMap(numWindows, vector<bool>(rtGridSize,0));  // location on XIC: one per retention time peak apex
    vector< vector<double> > preLCPeakApexIntensity(numWindows, vector<double>(rtGridSize,0)); // the intensity at each apex
    vector< vector<int> > preLCPeakApexIntervalStartIndex(numWindows, vector<int>(rtGridSize,0));  // scanid range , start, one per retention time peak apex
    vector< vector<int> > preLCPeakApexIntervalEndIndex(numWindows, vector<int>(rtGridSize,0));  // scanid range , end, one per retention time peak apex
    vector< vector< vector<double> > > preLCPeakApexIntervalProfile(numWindows, vector< vector<double> >(rtGridSize, vector<double>() ));  // the rt intensi

    ///////////////////////////////////////////////////////////////////////////
    // Add debug data structures for preLCPeakApexMap
    ///////////////////////////////////////////////////////////////////////////
    vector<int> plcpeakIntervalRtGridStartId;
    vector<int> plcpeakIntervalRtGridEndId;
    vector< vector<double> > plcpeakIntervalSmoothedSignal;
    vector<double> plcpeakIntervalPeakApexIntensity;
    vector<int> plcpeakIntervalPeakApexMzGridId;
    vector<int> plcpeakIntervalPeakApexRtGridId;
    int plcpeakIntervalCounter = 0;
    ///////////////////////////////////////////////////////////////////////////
    // Add debug data structures for preLCPeakApexMap
    ///////////////////////////////////////////////////////////////////////////

    vector< vector<bool> > LCPeakApexMap(numWindows, vector<bool>(rtGridSize,0));
    vector< vector<double> > LCPeakApexIntensity(numWindows, vector<double>(rtGridSize,0));
    vector< vector<int> > LCPeakApexIntervalStartIndex(numWindows, vector<int>(rtGridSize,0));  // scanid range , start, one per retention time peak apex
    vector< vector<int> > LCPeakApexIntervalEndIndex(numWindows, vector<int>(rtGridSize,0));  // scanid range , end, one per retention time peak apex

    // The result of merging the Retention Time dimension and MZ dimension peaks
    vector< vector<bool> > LCPeakApexMapMerged(numWindows,vector<bool>(rtGridSize,0));

    vector<int> apexWindowId;
    vector<int> apexScanId;
    vector<int> apexWindowIdScanSorted;
    vector<int> apexScanIdScanSorted;

    ///////////////////////////////////////////////////////////////////////////
    // Add the mz value at the peak and the retention time at the peak
    // as information to keep within the apex data structures.
    // Added on 11/27/12
    ///////////////////////////////////////////////////////////////////////////
    vector<double> apexMz;
    vector<int> apexScan;
    vector<double> apexRt;
    vector<double> apexMzSorted;
    vector<int> apexScanSorted;
    vector<double> apexRtSorted;
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    // get the rt intervals ( preLCPeakApexMap )
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    for ( int i = 0; i < smoothXICs.size(); i++ ) {
        vector<int> intervalListStart;
        vector<int> intervalListEnd;
        int rc = 0;

        rc = utilitiesObject.getInterval( 	smoothXICs.at(i),
                                            noiseThresholdLevel,
                                            minLCLength,
                                            intervalListStart,
                                            intervalListEnd );

        // Assert can be disabled by compiling with NDEBUG option
        assert( intervalListStart.size() == intervalListEnd.size() );

        // call the splitInterval function
        // We'll need the diff of the current smoothXIC
        vector<double> diffXICBuffer;
        rc = utilitiesObject.getDiff(smoothXICs.at(i), -1, diffXICBuffer);

        vector<int> outIntervalListStart; // return variable
        vector<int> outIntervalListEnd;  // return variable

        rc = utilitiesObject.splitInterval( smoothXICs.at(i),
                                            diffXICBuffer,
                                            intervalListStart,
                                            intervalListEnd,
                                            minLCLength,
                                            outIntervalListStart,
                                            outIntervalListEnd );

        // Looop through all the intervals, finding the maximum intensity value
        // within each interval and saving the results in the
        // preLCPeakApex data structures.
        for ( int j = 0; j < outIntervalListStart.size(); j++ ) {
            // Populate the preLCPeakApex data structures
            vector<double> intervalIntensityData;     // temporary buffer for the intervals intensity data
            // Bug fix:
            // change < to <= on 10/11/12, size mismatch when comparing
            // elution profiles.
            // Note:  The interval starting and stop indices are inclusive
            for ( int k = outIntervalListStart[j]; k <= outIntervalListEnd[j]; k++ ) {
                intervalIntensityData.push_back(smoothXICs.at(i).at(k));
            } // end loop through the data within an interval ( k variable )

            ///////////////////////////////////////////////////////////////////
            // Find the maximum intensity value within an interval
            // We need the following 3 pieces of information for
            // each peak:
            ///////////////////////////////////////////////////////////////////
            vector<double>::size_type peakIntensityIndex;
            double peakIntensityValue=0;
            rc = sgObject.findMaxValue(intervalIntensityData,peakIntensityIndex,peakIntensityValue);

            // Convert the local index within an interval into the index
            // within the smoothXICs.at(i) data structure
            vector<double>::size_type peakIntensityIndexGlobal = outIntervalListStart[j] + peakIntensityIndex;

            preLCPeakApexMap.at(i).at(peakIntensityIndexGlobal) = 1;  // label this as an rt peak
            preLCPeakApexIntensity.at(i).at(peakIntensityIndexGlobal) = peakIntensityValue; // the intensity at each pre-LC apex
            preLCPeakApexIntervalStartIndex.at(i).at(peakIntensityIndexGlobal) = outIntervalListStart[j]; // one item per pre-lc apex
            preLCPeakApexIntervalEndIndex.at(i).at(peakIntensityIndexGlobal) = outIntervalListEnd[j]; // one item per pre-lc apex
            preLCPeakApexIntervalProfile.at(i).at(peakIntensityIndexGlobal) = intervalIntensityData;  // save the interval's intensity data

            ///////////////////////////////////////////////////////////////////
            // DEBUG
            ///////////////////////////////////////////////////////////////////
            plcpeakIntervalRtGridStartId.push_back(outIntervalListStart[j]);
            plcpeakIntervalRtGridEndId.push_back(outIntervalListEnd[j]);
            plcpeakIntervalSmoothedSignal.push_back(intervalIntensityData);
            plcpeakIntervalPeakApexIntensity.push_back(peakIntensityValue);
            plcpeakIntervalPeakApexMzGridId.push_back(i);
            plcpeakIntervalPeakApexRtGridId.push_back(peakIntensityIndexGlobal);
            plcpeakIntervalCounter = plcpeakIntervalCounter +1;
            ///////////////////////////////////////////////////////////////////
            // DEBUG
            ///////////////////////////////////////////////////////////////////




        } // end loop through all the intervals within this window ( j variable )
        ///////////////////////////////////////////////////////////////////////
    } // end loop through all the windows ( i variable )
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    // Get the mz intervals and generate LCPeakApexMap data
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////


    ///////////////////////////////////////////////////////////////////////////
    // Get the mz intervals and generate LCPeakApexMap data
    // loop through all the scans ( i variable )
    // preLCPeakApexIntensity
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    // Add debug data structures for mz intervals
    ///////////////////////////////////////////////////////////////////////////
    vector<int> mzIntervalMzGridStartId;
    vector<int> mzIntervalMzGridEndId;
    vector< vector<double> > mzIntervalLcIntervalMzSignal;
    vector<double> mzIntervalPeakLocalMaxInt;
    vector<int> mzIntervalPeakLocalMaxId;
    vector<int> mzIntervalPeakMzGridId;
    vector<int> mzIntervalPeakRtGridId;
    int mzIntervalCounter = 0;
    ///////////////////////////////////////////////////////////////////////////
    // Add debug data structures for preLCPeakApexMap
    ///////////////////////////////////////////////////////////////////////////


    for ( int i = 0; i < rtGridSize; i++ ){
        double mzIntervalThreshold = 0; // maybe make an input parameter later on
        int mzIntervalMinLength = 2; // mabe make an input parameter later on
        int scanIndexLow  = sgObject.max(0,i-LCPeakApexTolerance); // make sure we are in bounds(left)
        int scanIndexHigh = sgObject.min(i+LCPeakApexTolerance,rtGridSize-1); // make sure we are in bounds(right)


        // create a vector to hold the sum of each window
        // the length of the sumScansVector is the total number of mz windows
        vector<double> sumScansVector(numWindows,0);

        // generate the summation in the mz dimension
        for ( int k = 0; k < numWindows; k++ ) {
            for ( int j = scanIndexLow; j <= scanIndexHigh; j++ ) {
                sumScansVector.at(k) = sumScansVector.at(k)+preLCPeakApexIntensity.at(k).at(j);
            }  // loop through group of columns
        } // process all mz data for the group of columns

        // Interval detection along the mz dimension
        // now call interval detection to get intervals along the mz dimension
        int rc = 0;
        vector<int> mzIntervalListStart;
        vector<int> mzIntervalListEnd;

        rc = utilitiesObject.getInterval( 	sumScansVector,
                                            mzIntervalThreshold,
                                            mzIntervalMinLength,
                                            mzIntervalListStart,
                                            mzIntervalListEnd );

        ///////////////////////////////////////////////////////////////////////
        //
        // CRITICAL UPDATE:
        // "We didn't split the interval in previous version and this is a
        //  cause for many missing peaks" (MZ 02/26/13)
        //
        // We need to split the interval here.
        //
        //
        ///////////////////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////////////////
        //
        // We need to split the interval along the mz dimension in order
        // to be able to properly find the boundary of the lc candidate.
        //
        ///////////////////////////////////////////////////////////////////////
        vector<double> diffSumScanVector;
        rc = utilitiesObject.getDiff(sumScansVector, -1, diffSumScanVector);

        vector<int> outIntervalListMzStart; // contains split signals along mz dimension
        vector<int> outIntervalListMzEnd;  // contains split signals along mz dimension

        int mzSplitIntervalMinLength = 2; // We need to make this an input parameter

        rc = utilitiesObject.splitInterval( sumScansVector,
                                            diffSumScanVector,
                                            mzIntervalListStart,
                                            mzIntervalListEnd,
                                            mzSplitIntervalMinLength,
                                            outIntervalListMzStart,
                                            outIntervalListMzEnd );

        assert( outIntervalListMzStart.size() == outIntervalListMzEnd.size() );

        // go through each interval found in the mz dimension
        // finding the maximum intensity value within the interval
        // as well as the relative and absolute indices.
        for ( int j = 0; j < outIntervalListMzStart.size(); j++ ) {
            // create a temporary copy of the data within an interval
            vector<double> intervalDataBuffer;
            for (int k = outIntervalListMzStart.at(j); k <= outIntervalListMzEnd.at(j); k++ ){
                intervalDataBuffer.push_back(sumScansVector.at(k));
            }
            // find the maximum value within the interval
            // the findMaxValue function will return
            vector<double>::size_type peakIntensityIndex=0;
            vector<double>::size_type peakIntensityIndexGlobal=0;
            double peakIntensityValue=0;
            int rc = sgObject.findMaxValue(intervalDataBuffer,peakIntensityIndex,peakIntensityValue);
            peakIntensityIndexGlobal =  peakIntensityIndex+outIntervalListMzStart.at(j);

            // Place a marker '1' at each location corresponding to a peak in the mz dimension
            LCPeakApexMap.at(peakIntensityIndexGlobal).at(i) = 1;
            LCPeakApexIntensity.at(peakIntensityIndexGlobal).at(i) = peakIntensityValue;
            LCPeakApexIntervalStartIndex.at(peakIntensityIndexGlobal).at(i) = outIntervalListMzStart.at(j);
            LCPeakApexIntervalEndIndex.at(peakIntensityIndexGlobal).at(i) = outIntervalListMzEnd.at(j);


            ///////////////////////////////////////////////////////////////////
            // DEBUG
            ///////////////////////////////////////////////////////////////////
            mzIntervalMzGridStartId.push_back(outIntervalListMzStart.at(j));  // global index on mz grid of mz interval start
            mzIntervalMzGridEndId.push_back(outIntervalListMzEnd.at(j));  // global index on mz grid of mz interval end
            mzIntervalLcIntervalMzSignal.push_back(intervalDataBuffer); // mz interval signal
            mzIntervalPeakLocalMaxInt.push_back(peakIntensityValue);     // maximum intensity of mz peak local mz axis
            mzIntervalPeakLocalMaxId.push_back(peakIntensityIndex);         // local index of mz peak local mz axis
            mzIntervalPeakMzGridId.push_back(peakIntensityIndexGlobal);  // global index of mz peak on mz grid
            mzIntervalPeakRtGridId.push_back(i);      // global index of mz peak on rt grid
            mzIntervalCounter = mzIntervalCounter+1;
            ///////////////////////////////////////////////////////////////////
            // END DEBUG
            ///////////////////////////////////////////////////////////////////

        } // end loop through mz intervals found
    } // end loop through all the mz dimension signals ( each signal has a single rt value, and varies along the mz dimension )



    ///////////////////////////////////////////////////////////////////////////
    // merged rt and mz interval apex data ( LCPeakApexMap )
    //
    // Note: With the changed dimensions of the PLCROI data structures,
    //
    ///////////////////////////////////////////////////////////////////////////
    for ( int j = 0; j < preLCPeakApexMap.size(); j++ ) {
        for( int k = 0; k < preLCPeakApexMap.at(j).size(); k++ ) {
            bool plcData = preLCPeakApexMap.at(j).at(k);  // convert from int to boolean
            bool lcData = LCPeakApexMap.at(j).at(k);      // convert from int to boolean
            LCPeakApexMapMerged.at(j).at(k) = (plcData && lcData);
            // save the coordinates of the peak apex
            // j --> maps to apexWindowId
            // k --> maps to apexScanId
            if ( LCPeakApexMapMerged.at(j).at(k) == 1) {
                apexWindowId.push_back(j);
                apexScanId.push_back(k);
            } // end check so we can save the peak location ( windowid and scanid )
        } // end loop through all scan#'s( rt dimension )
    } // end loop through all mz windows ( mz dimension )


    ///////////////////////////////////////////////////////////////////////////
    //Bug fix:
    //matlab code compatibility issue:
    //detected: 10/24/12. The find command in Matlab is
    //returning the list of apex peaks sorted in ascending
    //order by the scan dimension.
    //In order to match the Matlab code, we need to
    //obtain the sorted indices of the apexScanId vector
    //in ascending order, and then generated
    //a sorted set of apex index vectors.
    ///////////////////////////////////////////////////////////////////////////
    vector<int> apexScanIndexVector;
    for ( int i = 0; i < apexScanId.size(); i++ ) {
        apexScanIndexVector.push_back(i);
    }

    rc = 0;
    rc = sgObject.sortGetIndex(apexScanId, apexScanId, apexScanIndexVector, true);

    assert(rc == 0);

    // re-organize the apexWindowId and apexScanId lists according
    // to the ordering specified in the apex
    for ( int i = 0; i < apexScanIndexVector.size(); i++ ) {
        apexWindowIdScanSorted.push_back( apexWindowId.at(apexScanIndexVector.at(i)) );
        apexScanIdScanSorted.push_back( apexScanId.at(apexScanIndexVector.at(i)) );
    }


    ///////////////////////////////////////////////////////////////////////////
    // Validation succeeds up to this point as of 10/24/12
    //
    // The following 3 data structures contain the results of all
    // processing up to this point.
    //
    // LCPeakApexMapMerged       --> same size as smoothXIC ( a boolean overlay )
    // apexWindowIdScanSorted    --> index into smoothXIC ( row# )
    // apexScanIdScanSorted      --> index into smoothXIC ( col# )
    //
    // Validated against Matlab code on 10/24/12
    // After Algorithm Update #2 ( 03/20/13 ).  Version with
    // resampled retention time & mz axis.
    // This was a critical algorithmic improvement.
    // Also, performing interval detection on the the retention time and
    // mz axis.
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // ------START----------------
    // LC Candiate Data structures
    // There must be enough information in these data structures
    // such that an LC Candidate can be a completely separate
    // unit ( object ).
    // Important: This is the starting point for the design
    // of an LC Candidate Object.  Once all there general characteristics
    // of what an LC Candidate is are fully worked through, both
    // the PLCROI and LCC data structures can be converted into C++
    // objects.
    ///////////////////////////////////////////////////////////////////////////
    int totalPeakApex = apexWindowIdScanSorted.size();
    vector<int> lccPlcRoi(totalPeakApex,plcRoi);
    vector<double> lccMzGridWidth(totalPeakApex,mzGridWidth);
    vector<double> lccRtDelta(totalPeakApex,rtDelta);
    vector<double> lccCenterMz(totalPeakApex,0);
    vector<double> lccStartMz(totalPeakApex,0);
    vector<double> lccEndMz(totalPeakApex,0);
    vector<double> lccPeakApexElutionTime(totalPeakApex,0);
    vector<double> lccElutionTimeStart(totalPeakApex,0);
    vector<double> lccElutionTimeEnd(totalPeakApex,0);
    vector<int> lccCenterMzGridId(totalPeakApex,0);  // used for estimatedElutionProfile
    vector< vector<double> > lccSmoothElutionProfile(totalPeakApex, vector<double>() );
    vector< vector<double> > lccEstimatedElutionProfile(totalPeakApex, vector<double>() );
    vector< vector<double> > lccMzPeakProfile(totalPeakApex, vector<double>() );
    vector<bool> lccValidCandidate(totalPeakApex,false);  // The index into this array is the lcc's identifier

    //Additional Level 3 Data
    vector<int> lccStartMzGridId(totalPeakApex,0);
    vector<int> lccEndMzGridId(totalPeakApex,0);
    vector<int> lccScanStartID(totalPeakApex,0);
    vector<int> lccScanEndID(totalPeakApex,0);
    vector<int> lccPeakApexRtGridId(totalPeakApex,0);

    vector< vector< vector<double> > > lccRegionData;
    for ( int k = 0; k < totalPeakApex; k++) {
        vector< vector<double> > tempBuffer;
        lccRegionData.push_back(tempBuffer);
    }
    vector< vector<double> > lccPeakMzGrid(totalPeakApex, vector<double>() );  // the piece of mz grid for this lc candidate
    vector< vector<double> > lccPeakRtGrid(totalPeakApex, vector<double>() );  // the piece of rt grid for this lc candidate
    ///////////////////////////////////////////////////////////////////////////
    // ------END-----------------
    // LC Candiate Data structures
    // There must be enough information in these data structures
    // such that an LC Candidate can be a completely separate
    // unit ( object ).
    // Important: This is the starting point for the design
    // of an LC Candidate Object.  Once all there general characteristics
    // of what an LC Candidate is are fully worked through, both
    // the PLCROI and LCC data structures can be converted into C++
    // objects.
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // Use the resampledRegion to estimate the mz boundary of the
    // LC candidate, then populate the LC Candidate data structures.
    // The following 3 data structures contain the results of all
    // processing up to this point.
    //
    // LCPeakApexMapMerged       --> same size as smoothXIC ( a boolean overlay )
    // apexWindowIdScanSorted    --> index into smoothXIC ( row# )
    // apexScanIdScanSorted      --> index into smoothXIC ( col# )
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // DEBUG
    // peakApex debug data structures, one per peakApex
    //
    ///////////////////////////////////////////////////////////////////////////
    vector<int> apexScanRangeStartVector;
    vector<int> apexScanRangeEndVector;
    vector<int> apexWindowRangeStartVector;
    vector<int> apexWindowRangeEndVector;
    vector< vector<double> > apexElutionProfileVector;  // a signal per peakApex
    vector< vector<int> > tempidsVector;  // a signal per peakApex

    // The following are only conditionally set depending on the length of mz
    vector<int> newApexRtGridIDRangeStartVector(totalPeakApex,0);
    vector<int> newApexRtGridIDRangeEndVector(totalPeakApex,0);
    vector<double> peakApexMZVector(totalPeakApex,0);
    vector<int> peakApexmzGridIDVector(totalPeakApex,0);
    vector< vector<double> > apexElutionProfileShortVector(totalPeakApex,vector<double>() );
    vector< vector<double> > smoothApexElutionProfileLongVector(totalPeakApex,vector<double>() );

    // Debugging data structures for finding mz grid limits of lc candidate
    vector<double> sumApexElutionProfileVector(totalPeakApex,0);
    vector<double> threeSTDVector(totalPeakApex,0);
    vector< vector<int> > mzGridIDsVector(totalPeakApex,vector<int>() );
    vector< vector<double> > rcorrMapVector(totalPeakApex,vector<double>() );
    vector< vector<int> > tempidsRcorrMapVector(totalPeakApex,vector<int>() );


    ///////////////////////////////////////////////////////////////////////////
    // DEBUG
    // peakApex debug data structures
    //
    ///////////////////////////////////////////////////////////////////////////

    for ( int i = 0; i < totalPeakApex; i++ ){
        // get the scan range for this apex
        // apexWindowIdScanSortec -> apexRtGridIDRange(Matlab)
        int apexScanRangeStart = preLCPeakApexIntervalStartIndex.at(apexWindowIdScanSorted.at(i)).at(apexScanIdScanSorted.at(i));
        int apexScanRangeEnd = preLCPeakApexIntervalEndIndex.at(apexWindowIdScanSorted.at(i)).at(apexScanIdScanSorted.at(i));



        int apexWindowRangeStart = LCPeakApexIntervalStartIndex.at(apexWindowIdScanSorted.at(i)).at(apexScanIdScanSorted.at(i));
        int apexWindowRangeEnd = LCPeakApexIntervalEndIndex.at(apexWindowIdScanSorted.at(i)).at(apexScanIdScanSorted.at(i));


        // Find the apex elution profile from the XIC data structure
        vector<double> apexElutionProfile = preLCPeakApexIntervalProfile.at(apexWindowIdScanSorted.at(i)).at(apexScanIdScanSorted.at(i));
        vector<double> apexElutionProfileShort;  // holds the retention time cut signals
        vector<double> smoothApexElutionProfileLong; // holds the full smoothed elution profile of the peak apex
        ///////////////////////////////////////////////////////////////////////
        //
        // Trim the elution profile, be removing items that
        // fall below a threshold intensity value.
        //
        ///////////////////////////////////////////////////////////////////////

        // Find the min and max values of the apexElutionProfile
        // Calculate an intensity threshold
        double apexElutionProfileMin=0;
        double apexElutionProfileMax=0;
        double apexElutionProfileDiff = 0;
        double apexElutionProfileThreshold = 0.1;
        double apexElutionProfileThresholdFinal =0;
        rc = sgObject.findMinValue(apexElutionProfile,apexElutionProfileMin);
        rc = sgObject.findMaxValue(apexElutionProfile,apexElutionProfileMax);
        apexElutionProfileDiff = apexElutionProfileMax-apexElutionProfileMin;
        apexElutionProfileThresholdFinal= apexElutionProfileThreshold * apexElutionProfileDiff;


        // intensity threshold indices
        vector<int> tempids;  // maps to tempids

        for ( int j = 0; j < apexElutionProfile.size(); j++ ){
            if(  (apexElutionProfile.at(j) - apexElutionProfileMin) > apexElutionProfileThresholdFinal ) {
                tempids.push_back(j);
            } // end if
        } // end for




        ///////////////////////////////////////////////////////////////////////
        //
        // DEBUG START
        //
        ///////////////////////////////////////////////////////////////////////
        apexScanRangeStartVector.push_back(apexScanRangeStart);
        apexScanRangeEndVector.push_back(apexScanRangeEnd);
        apexWindowRangeStartVector.push_back(apexWindowRangeStart);
        apexWindowRangeEndVector.push_back(apexWindowRangeEnd);
        apexElutionProfileVector.push_back(apexElutionProfile);
        tempidsVector.push_back(tempids);
        ///////////////////////////////////////////////////////////////////////
        //
        // DEBUG END
        //
        ///////////////////////////////////////////////////////////////////////


        ///////////////////////////////////////////////////////////////////////
        // add a check to make sure there are sufficient items in tempids
        // this maps to if length(tempids)>minMZLength%add by MZ 02/05/2013
        // This maps to line 643 of the genlccandidatesSimplestVersionWorking021213.m
        // file, which contains the original logic.
        ///////////////////////////////////////////////////////////////////////
        if ( tempids.size() >= minMZLength ) {
            // minimum mz length criteria met

            ///////////////////////////////////////////////////////////////////////
            //
            // Key improvement:
            // We need to compare equivalent signals.  If we made it to this
            // point, we know that there are enough data points in the mz signal.
            // Instead of using the smoothed elution profile,  we need to
            // use the signal directly on the resampled region.
            //
            ///////////////////////////////////////////////////////////////////////

            // calculate new scanIndexStart and scanIndexEnd
            int scanIndexStartNewLocal = 0;
            int scanIndexEndNewLocal = 0;
            int scanIndexStartNewGlobal = 0;
            int scanIndexEndNewGlobal = 0;

            // Find the min and max values of the intensity profile
            rc = sgObject.findMinValue(tempids,scanIndexStartNewLocal);
            rc = sgObject.findMaxValue(tempids,scanIndexEndNewLocal);

            // Add the global offset back to the new starting and ending indices
            // from the thresholding stage
            scanIndexStartNewGlobal = scanIndexStartNewLocal + apexScanRangeStart;  // add the starting offset
            scanIndexEndNewGlobal = scanIndexEndNewLocal + apexScanRangeStart;  // add the starting offset


            ///////////////////////////////////////////////////////////////////
            // START DEBUG
            ///////////////////////////////////////////////////////////////////
            newApexRtGridIDRangeStartVector.at(i) = scanIndexStartNewGlobal; // The new rt grid start point, based on intensity thresholding
            newApexRtGridIDRangeEndVector.at(i) = scanIndexEndNewGlobal; // The new rt grid end point, based on intensity thresholding
            ///////////////////////////////////////////////////////////////////
            // END DEBUG
            ///////////////////////////////////////////////////////////////////



            // Now we set out to discover the range of mz values on the original mz grid
            assert ( apexWindowIdScanSorted.at(i) <  mzGrid.size() );
            // the apexWindowIdScanSorted is the index within the PLCROI resampled
            // data structures along the mz dimension
            double peakApexMZ  = mzGrid.at( apexWindowIdScanSorted.at(i) );
            int peakApexMzGridId = apexWindowIdScanSorted.at(i);


            ///////////////////////////////////////////////////////////////////
            // START DEBUG
            ///////////////////////////////////////////////////////////////////
            peakApexMZVector.at(i) = peakApexMZ;
            peakApexmzGridIDVector.at(i) = peakApexMzGridId;
            ///////////////////////////////////////////////////////////////////
            // END DEBUG
            ///////////////////////////////////////////////////////////////////



            assert( scanIndexStartNewGlobal < resampledRegionXICs.size() );
            assert( scanIndexEndNewGlobal < resampledRegionXICs.size() );
            // range is inclusive
            for ( int k = scanIndexStartNewGlobal; k <= scanIndexEndNewGlobal; k++) {
                apexElutionProfileShort.push_back( resampledRegionXICs.at(k).at(peakApexMzGridId) );
            }

            ///////////////////////////////////////////////////////////////////
            // START DEBUG
            ///////////////////////////////////////////////////////////////////
            apexElutionProfileShortVector.at(i)= apexElutionProfileShort;
            ///////////////////////////////////////////////////////////////////
            // END DEBUG
            ///////////////////////////////////////////////////////////////////

            // we also need to keep the full smoothed elution profile
            // apexScanRangeStart, apexScanRangeEnd
            assert( apexScanRangeStart < rtGridSize );
            assert( apexScanRangeEnd < rtGridSize );


            for ( int k = apexScanRangeStart; k <= apexScanRangeEnd; k++ ){
                smoothApexElutionProfileLong.push_back(smoothXICs.at(peakApexMzGridId).at(k));
            }


            ///////////////////////////////////////////////////////////////////
            // START DEBUG
            ///////////////////////////////////////////////////////////////////
            smoothApexElutionProfileLongVector.at(i)= smoothApexElutionProfileLong;
            ///////////////////////////////////////////////////////////////////
            // END DEBUG
            ///////////////////////////////////////////////////////////////////



            ///////////////////////////////////////////////////////////////////
            //
            // We need to find the starting and ending mz region
            // using the cut out apexElutionProfileShort signals
            // as the signals whose shape will be compared to determine
            // the mz boundary.
            //
            ///////////////////////////////////////////////////////////////////

            // calculate the sum of the shortened elution profile signal
            // which has been thresholded
            double sumApexElutionProfileShort = 0;
            for ( int k = 0; k < apexElutionProfileShort.size(); k++ ){
                sumApexElutionProfileShort = sumApexElutionProfileShort + apexElutionProfileShort.at(k);
            }



            ///////////////////////////////////////////////////////////////////////
            // Find the range around the peak ( this line maps to:
            // 659 threeSTD=2*peakApexMZ/(massResolution*2)/sqrt(2*log(2));
            ///////////////////////////////////////////////////////////////////////
            double peakWidthStd = 2*peakApexMZ/(massResolution*2)/sqrt(2*log(2)); // maps to threeSTD


            ///////////////////////////////////////////////////////////////////////
            //
            // Use the resampledGrid to determine the mz start and end points
            // of the LC Candidate based on a peak shape metric of neighboring
            // signals.
            //
            // Find the indices of the items in the mzGrid that fall within
            // the range specified by peakWidthStd.
            //
            //    peakApexMZ
            //
            ///////////////////////////////////////////////////////////////////////
            double mzStart = peakApexMZ-peakWidthStd;  // starting mz value from the current apex
            double mzEnd = peakApexMZ+peakWidthStd;    // ending mz value from the current apex
            vector<int> mzGridIds;
            for ( int j = 0; j < mzGrid.size(); j++ ){
                if( (mzGrid.at(j) >= mzStart) && (mzGrid.at(j) <= mzEnd) ){
                    mzGridIds.push_back(j);
                }
            }

            ///////////////////////////////////////////////////////////////////////
            //
            // Create a vector to store the R^2 correlations between signals
            //
            ///////////////////////////////////////////////////////////////////////
            // Find the minimum value
            int minMzGridIds = 0;
            rc = sgObject.findMinValue(mzGridIds,minMzGridIds);

            vector<double> rCorrMap(mzGridIds.size(),0);  // keep track of all correlation values between elution profiles


            for( int j = 0; j < mzGridIds.size(); j++ ){
                vector<double> currentProfile; // temporary buffer( new for each j )

                ///////////////////////////////////////////////////////////////
                // Critical Bug Fix:
                // There are cases where the resampledRegion could
                // be empty, due to there no being sufficient points
                // to perform the resampling process.
                ///////////////////////////////////////////////////////////////

                for( int k = scanIndexStartNewGlobal; k <= scanIndexEndNewGlobal; k++ ){
                    //go through the current signal at a specified mz value
                    if ( resampledRegionXICs.at(k).size() > 0 ) {
                        currentProfile.push_back(resampledRegionXICs.at(k).at(j+minMzGridIds));
                    }  // end check whether signal exists
                }
                rCorrMap.at(j) = 0;
                double rCorrBuffer = 0;

                if ( currentProfile.size() == apexElutionProfileShort.size() ){
                    rc = utilitiesObject.getR2Statistic( currentProfile,
                                                         apexElutionProfileShort,
                                                         rCorrBuffer);
                }

                // populate the rCorrMap vector
                rCorrMap.at(j) = rCorrBuffer;

            } // loop through mzGridIds



            ///////////////////////////////////////////////////////////////////////
            //
            // Obtain a list of the items in the rCorrMap that are above
            // a threshold for signal similarity.
            //
            ///////////////////////////////////////////////////////////////////////
            vector<int> rctempids;
            for ( int j = 0; j < rCorrMap.size(); j++ ){
                if ( rCorrMap.at(j) > rawMZCorrThreshold ) {
                    rctempids.push_back(j);
                }
            }


            ///////////////////////////////////////////////////////////////////
            // DEBUG START
            ///////////////////////////////////////////////////////////////////
            sumApexElutionProfileVector.at(i) = sumApexElutionProfileShort;
            threeSTDVector.at(i)= peakWidthStd;
            mzGridIDsVector.at(i) = mzGridIds;
            rcorrMapVector.at(i) = rCorrMap;
            tempidsRcorrMapVector.at(i) = rctempids;
            ///////////////////////////////////////////////////////////////////
            // DEBUG END
            ///////////////////////////////////////////////////////////////////




            if( rctempids.size() >= minMZLength ) {

                lccValidCandidate.at(i)=true;  // The index into this array is the lcc's identifier
                ///////////////////////////////////////////////////////////////
                // MZ Grid Start and End Indices
                // RT Grid Start and End Indices
                ///////////////////////////////////////////////////////////////
                int minRcTempIds = 0;
                rc = sgObject.findMinValue(rctempids,minRcTempIds);
                lccStartMzGridId.at(i)=minRcTempIds+minMzGridIds;
                int maxRcTempIds = 0;
                rc = sgObject.findMaxValue(rctempids,maxRcTempIds);
                lccEndMzGridId.at(i)=maxRcTempIds+minMzGridIds; //bug fix 03/22/13
                // resolution: We were using minRcTempIds instead of minMzGridIds,
                // and therefore we were not obtaining the correct indices
                // also, we were using minRcTempIds instead of maxRcTempIds...
                //
                lccScanStartID.at(i) = apexScanRangeStart;
                lccScanEndID.at(i) = apexScanRangeEnd;


                ///////////////////////////////////////////////////////////////
                // MZ Grid Start and End Data Values
                // RT Grid Start and End Data Values
                ///////////////////////////////////////////////////////////////
                lccStartMz.at(i)=mzGrid.at(lccStartMzGridId.at(i));
                lccEndMz.at(i)=mzGrid.at(lccEndMzGridId.at(i));
                lccElutionTimeStart.at(i) = rtGrid.at(apexScanRangeStart);
                lccElutionTimeEnd.at(i) = rtGrid.at(apexScanRangeEnd);

                ///////////////////////////////////////////////////////////////
                //Partition the section of the resampledRegion
                //belonging to this lcc.
                ///////////////////////////////////////////////////////////////
                vector< vector<double> > resampledRegionXICsBuffer;
                for (int k = apexScanRangeStart; k <= apexScanRangeEnd; k++ ) {
                    vector<double> mzSectionBuffer;
                    for ( int r = lccStartMzGridId.at(i); r <= lccEndMzGridId.at(i); r++ ){
                        mzSectionBuffer.push_back(resampledRegionXICs.at(k).at(r));
                    }
                    resampledRegionXICsBuffer.push_back(mzSectionBuffer);
                }
                lccRegionData.at(i) = resampledRegionXICsBuffer;


                ///////////////////////////////////////////////////////////////
                // The elution profile data from the resampled region
                ///////////////////////////////////////////////////////////////
                lccSmoothElutionProfile.at(i) = smoothApexElutionProfileLong;


                ///////////////////////////////////////////////////////////////
                // The mz grid corresponding to this LC Candidate
                ///////////////////////////////////////////////////////////////
                vector<double> peakMzGridBuffer;
                for ( int k = lccStartMzGridId.at(i); k <= lccEndMzGridId.at(i); k++){
                    peakMzGridBuffer.push_back(mzGrid.at(k));
                }
                lccPeakMzGrid.at(i)= peakMzGridBuffer;


                ///////////////////////////////////////////////////////////////
                // The rt grid corresponding to this LC Candidate
                ///////////////////////////////////////////////////////////////
                vector<double> peakRtGridBuffer;
                for ( int k = apexScanRangeStart; k <= apexScanRangeEnd; k++){
                    peakRtGridBuffer.push_back(rtGrid.at(k));
                }
                lccPeakRtGrid.at(i)= peakRtGridBuffer;


                ///////////////////////////////////////////////////////////////
                //
                // This is the additive profile collapsing all the rows
                // into a single row vector with the same number of
                // columns.
                ///////////////////////////////////////////////////////////////
                vector<double> mzPeakProfileBuffer( peakMzGridBuffer.size(), 0);
                for(int k = 0; k < resampledRegionXICsBuffer.size(); k++ ){
                    for ( int r = 0; r < resampledRegionXICsBuffer.at(k).size();  r++ ) {
                        mzPeakProfileBuffer.at(r) = mzPeakProfileBuffer.at(r) + resampledRegionXICsBuffer.at(k).at(r);
                    }
                }
                lccMzPeakProfile.at(i) = mzPeakProfileBuffer;



                ///////////////////////////////////////////////////////////////
                // We need to find the maximum value peakMzGridBuffer, both
                // the value and the index location.
                ///////////////////////////////////////////////////////////////




                ///////////////////////////////////////////////////////////////
                // We need to find the maximum value peakMzGridBuffer, both
                // the value and the index location.
                ///////////////////////////////////////////////////////////////
                vector<double>::size_type peakIntensityIndex=-1;
                double peakIntensityValue=0;
                rc = sgObject.findMaxValue(mzPeakProfileBuffer,peakIntensityIndex,peakIntensityValue);

                lccCenterMzGridId.at(i) = peakIntensityIndex + lccStartMzGridId.at(i);

                ///////////////////////////////////////////////////////////////
                // The mz value at the peak
                ///////////////////////////////////////////////////////////////
                lccCenterMz.at(i) = mzGrid.at(lccCenterMzGridId.at(i));

                ///////////////////////////////////////////////////////////////
                // The index into the retention time grid axis at the peak
                ///////////////////////////////////////////////////////////////
                lccPeakApexRtGridId.at(i) = apexScanIdScanSorted.at(i);

                ///////////////////////////////////////////////////////////////
                // The elution time value at the peak
                ///////////////////////////////////////////////////////////////
                lccPeakApexElutionTime.at(i) = rtGrid.at(lccPeakApexRtGridId.at(i));


                ///////////////////////////////////////////////////////////////
                //
                //Key Item:  This is the data that will be used for
                // quantification.
                //
                ///////////////////////////////////////////////////////////////
                vector<double> estimatedElutionProfileBuffer;
                vector<double> peakSignalBuffer;
                for (int k = apexScanRangeStart; k <= apexScanRangeEnd; k++ ) {
                    peakSignalBuffer.push_back(resampledRegionXICs.at(k).at(lccCenterMzGridId.at(i)));
                }

                assert( resampledRegionXICsBuffer.size() == peakSignalBuffer.size() );

                rc = utilitiesObject.estimateElutionProfile(resampledRegionXICsBuffer,peakSignalBuffer,estimatedElutionProfileBuffer);

                if ( rc == 0 ){
                    //cout<<"We were able to estimate the elution profile"<<endl;
                    lccEstimatedElutionProfile.at(i)=estimatedElutionProfileBuffer;
                }
                else {
                    cout<<"elution profile cannot be estimated, rc="<<rc<<endl;
                }



            }
            else {
                lccValidCandidate.at(i) = false;
            }


            ///////////////////////////////////////////////////////////////////
            //For reference only:
            //vector< vector<double> > XICs(mzGridLength, vector<double>(rtGridSize,0));
            //vector< vector<double> > smoothXICs(mzGridLength, vector<double>(rtGridSize,0));
            //vector< vector<double> > resampledRegionXICs(rtGridSize, vector<double>(mzGridLength,0));
            ///////////////////////////////////////////////////////////////////

        } // end if
        else {
            // minimum mz length criteria not met
            lccValidCandidate.at(i) = false;
        } // end else

    } // end loop through all peaks apex


    ///////////////////////////////////////////////////////////////////////////
    //
    // Output the LC Candidate information for this PLCROI
    // to an MZDALevel 2 format file.
    //
    //
    ///////////////////////////////////////////////////////////////////////////

    // include the scan group number and the
    string valStringScanGroup;  // contains the scan group string information
    string valStringPLCROI;     // contains the pre-lc roi information
    string valStringSGPLCROI; // contains both scan group number and prelc roi number
    stringstream scangroupTemp (stringstream::in | stringstream::out);  // convert from int to string
    stringstream plcroiTemp (stringstream::in | stringstream::out); // convert from int to string
    // convert from integer to string
    scangroupTemp<<scanGroupId;
    scangroupTemp>>valStringScanGroup;
    plcroiTemp<<plcRoi;
    plcroiTemp>>valStringPLCROI;
    valStringSGPLCROI = valStringScanGroup.append("_");  // separate the sg and plcroi indicator
    valStringSGPLCROI = valStringSGPLCROI.append(valStringPLCROI);
    string outputFilenamePlcRoi = outputFilename.append(valStringSGPLCROI);
    outputFilenamePlcRoi.append(".m");

    ///////////////////////////////////////////////////////////////////////////
    //
    // Create the MZDA Level 2 file
    //
    ///////////////////////////////////////////////////////////////////////////
    ofstream mzdaLevelTwoOutputStream;
    mzdaLevelTwoOutputStream.open(outputFilenamePlcRoi.c_str(),std::ios_base::out );

    // Set the floating point precision of the output file stream
    mzdaLevelTwoOutputStream<<setprecision(10);





///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// START LOG LEVEL CONTROL BLOCK
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
if ( logLevel == 2 ) {
    // output the original scan number axis for the entire PLCROI
    mzdaLevelTwoOutputStream<<"originalScanNums=[";
    for ( int i = 0; i < originalScanNums.size(); i++ ) {
        mzdaLevelTwoOutputStream<<originalScanNums.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    // output the retention time axis for the entire PLCROI
    mzdaLevelTwoOutputStream<<"rtAxis=[";
    for ( int i = 0; i < rtAxis.size(); i++ ) {
        mzdaLevelTwoOutputStream<<rtAxis.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    // output the rtdiff data for the entire PLCROI
    mzdaLevelTwoOutputStream<<"rtDiff=[";
    for ( int i = 0; i < rtDiff.size(); i++ ) {
        mzdaLevelTwoOutputStream<<rtDiff.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    // output the rtdelta information
    mzdaLevelTwoOutputStream<<"rtDelta=["<<rtDelta<<"];";

    // output the rtmin information
    mzdaLevelTwoOutputStream<<"rtMin=["<<rtMin<<"];";

    // output the rtmin information
    mzdaLevelTwoOutputStream<<"rtMax=["<<rtMax<<"];";

    mzdaLevelTwoOutputStream<<"rtGrid=[";
    for ( int i = 0; i < rtGrid.size(); i++ ) {
        mzdaLevelTwoOutputStream<<rtGrid.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"rtGridSize=["<<rtGrid.size()<<"];";

    mzdaLevelTwoOutputStream<<"mzBegin=["<<mzBegin<<"];";
    mzdaLevelTwoOutputStream<<"mzLast=["<<mzLast<<"];";
    mzdaLevelTwoOutputStream<<"mzGridWidth=["<<mzGridWidth<<"];";

    mzdaLevelTwoOutputStream<<"mzGrid=[";
    for ( int i = 0; i < mzGrid.size(); i++ ) {
        mzdaLevelTwoOutputStream<<mzGrid.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"mzGridLength=["<<mzGridLength<<"];";

    mzdaLevelTwoOutputStream<<"resampledRegionXICs=[";
    for ( int i = 0; i < resampledRegionXICs.size(); i++ ) {
        for ( int j = 0; j < resampledRegionXICs.at(i).size(); j++ ){
            mzdaLevelTwoOutputStream<<resampledRegionXICs.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<";";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"XICs=[";
    for ( int i = 0; i < XICs.size(); i++ ) {
        for ( int j = 0; j < XICs.at(i).size(); j++ ){
            mzdaLevelTwoOutputStream<<XICs.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<";";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"preLCPeakApexMap=[";
    for ( int i = 0; i < preLCPeakApexMap.size(); i++ ) {
        for ( int j = 0; j < preLCPeakApexMap.at(i).size(); j++ ){
            mzdaLevelTwoOutputStream<<preLCPeakApexMap.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<";";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    // output the apex information
    mzdaLevelTwoOutputStream<<"apexWindowIdScanSortedSize=["<<apexWindowIdScanSorted.size()<<"];"<<endl;
    mzdaLevelTwoOutputStream<<"apexScanIdScanSorted=["<<apexScanIdScanSorted.size()<<"];"<<endl;
    mzdaLevelTwoOutputStream<<"apexWindowIdScanSorted=[";
    for ( int i = 0; i < apexWindowIdScanSorted.size(); i++ ) {
        mzdaLevelTwoOutputStream<<apexWindowIdScanSorted.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"apexScanIdScanSorted=[";
    for ( int i = 0; i < apexScanIdScanSorted.size(); i++ ) {
        mzdaLevelTwoOutputStream<<apexScanIdScanSorted.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"lccValidCandidate=[";
    for ( int i = 0; i < lccValidCandidate.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccValidCandidate.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"lccEstimatedElutionProfile={";
    for ( int i = 0; i < lccEstimatedElutionProfile.size(); i++ ) {
        mzdaLevelTwoOutputStream<<"[";
        for ( int j = 0; j < lccEstimatedElutionProfile.at(i).size(); j++ ) {
            mzdaLevelTwoOutputStream<<lccEstimatedElutionProfile.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<"];";
    }
    mzdaLevelTwoOutputStream<<"};"<<endl;
}  // END LOG LEVEL CONTROL BLOCK
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// START LOG LEVEL CONTROL BLOCK
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
if ( logLevel == 2 ) {

    mzdaLevelTwoOutputStream<<"plcpeakIntervalRtGridStartId=[";
    for ( int i = 0; i < plcpeakIntervalRtGridStartId.size(); i++ ) {
        mzdaLevelTwoOutputStream<<plcpeakIntervalRtGridStartId.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;


    mzdaLevelTwoOutputStream<<"plcpeakIntervalRtGridEndId=[";
    for ( int i = 0; i < plcpeakIntervalRtGridEndId.size(); i++ ) {
        mzdaLevelTwoOutputStream<<plcpeakIntervalRtGridEndId.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;


    mzdaLevelTwoOutputStream<<"plcpeakIntervalSmoothedSignal={";
    for ( int i = 0; i < plcpeakIntervalSmoothedSignal.size(); i++ ) {
        mzdaLevelTwoOutputStream<<"[";
        for ( int j = 0; j < plcpeakIntervalSmoothedSignal.at(i).size(); j++ ) {
            mzdaLevelTwoOutputStream<<plcpeakIntervalSmoothedSignal.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<"];";
    }
    mzdaLevelTwoOutputStream<<"};"<<endl;


    mzdaLevelTwoOutputStream<<"plcpeakIntervalPeakApexIntensity=[";
    for ( int i = 0; i < plcpeakIntervalPeakApexIntensity.size(); i++ ) {
        mzdaLevelTwoOutputStream<<plcpeakIntervalPeakApexIntensity.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;


    mzdaLevelTwoOutputStream<<"plcpeakIntervalPeakApexMzGridId=[";
    for ( int i = 0; i < plcpeakIntervalPeakApexMzGridId.size(); i++ ) {
        mzdaLevelTwoOutputStream<<plcpeakIntervalPeakApexMzGridId.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;


    mzdaLevelTwoOutputStream<<"plcpeakIntervalPeakApexRtGridId=[";
    for ( int i = 0; i < plcpeakIntervalPeakApexRtGridId.size(); i++ ) {
        mzdaLevelTwoOutputStream<<plcpeakIntervalPeakApexRtGridId.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;


    mzdaLevelTwoOutputStream<<"plcpeakIntervalCounter=["<<plcpeakIntervalCounter<<"];";


    mzdaLevelTwoOutputStream<<"mzIntervalMzGridStartId=[";
    for ( int i = 0; i < mzIntervalMzGridStartId.size(); i++ ) {
        mzdaLevelTwoOutputStream<<mzIntervalMzGridStartId.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"mzIntervalMzGridEndId=[";
    for ( int i = 0; i < mzIntervalMzGridEndId.size(); i++ ) {
        mzdaLevelTwoOutputStream<<mzIntervalMzGridEndId.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"mzIntervalLcIntervalMzSignal={";
    for ( int i = 0; i < mzIntervalLcIntervalMzSignal.size(); i++ ) {
        mzdaLevelTwoOutputStream<<"[";
        for ( int j = 0; j < mzIntervalLcIntervalMzSignal.at(i).size(); j++ ) {
            mzdaLevelTwoOutputStream<<mzIntervalLcIntervalMzSignal.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<"];";
    }
    mzdaLevelTwoOutputStream<<"};"<<endl;


    mzdaLevelTwoOutputStream<<"mzIntervalPeakLocalMaxInt=[";
    for ( int i = 0; i < mzIntervalPeakLocalMaxInt.size(); i++ ) {
        mzdaLevelTwoOutputStream<<mzIntervalPeakLocalMaxInt.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;


    mzdaLevelTwoOutputStream<<"mzIntervalPeakLocalMaxId=[";
    for ( int i = 0; i < mzIntervalPeakLocalMaxId.size(); i++ ) {
        mzdaLevelTwoOutputStream<<mzIntervalPeakLocalMaxId.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;


    mzdaLevelTwoOutputStream<<"mzIntervalPeakMzGridId=[";
    for ( int i = 0; i < mzIntervalPeakMzGridId.size(); i++ ) {
        mzdaLevelTwoOutputStream<<mzIntervalPeakMzGridId.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;


    mzdaLevelTwoOutputStream<<"mzIntervalPeakRtGridId=[";
    for ( int i = 0; i < mzIntervalPeakRtGridId.size(); i++ ) {
        mzdaLevelTwoOutputStream<<mzIntervalPeakRtGridId.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"mzIntervalCounter=["<<mzIntervalCounter<<"];";
}  // END LOG LEVEL CONTROL BLOCK
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// START LOG LEVEL CONTROL BLOCK
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
if ( logLevel == 2 ) {
    ///////////////////////////////////////////////////////////////////////////
    // LCPeakApexMapMerged       --> same size as smoothXIC ( a boolean overlay )
    // apexWindowIdScanSorted    --> index into smoothXIC ( row# )
    // apexScanIdScanSorted      --> index into smoothXIC ( col# )
    ///////////////////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"LCPeakApexMap=[";
    for ( int i = 0; i < LCPeakApexMap.size(); i++ ) {
        for ( int j = 0; j < LCPeakApexMap.at(i).size(); j++ ){
            mzdaLevelTwoOutputStream<<LCPeakApexMap.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<";";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;


    mzdaLevelTwoOutputStream<<"LCPeakApexMapMerged=[";
    for ( int i = 0; i < LCPeakApexMapMerged.size(); i++ ) {
        for ( int j = 0; j < LCPeakApexMapMerged.at(i).size(); j++ ){
            mzdaLevelTwoOutputStream<<LCPeakApexMapMerged.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<";";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;


    mzdaLevelTwoOutputStream<<"apexWindowIdScanSorted=[";
    for ( int i = 0; i < apexWindowIdScanSorted.size(); i++ ) {
        mzdaLevelTwoOutputStream<<apexWindowIdScanSorted.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;



    mzdaLevelTwoOutputStream<<"apexScanIdScanSorted=[";
    for ( int i = 0; i < apexScanIdScanSorted.size(); i++ ) {
        mzdaLevelTwoOutputStream<<apexScanIdScanSorted.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"totalPeakApex=["<<totalPeakApex<<"];";


    ///////////////////////////////////////////////////////////////////////////
    // PEAK APEX DATA STRUCTURES
    ///////////////////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"apexScanRangeStartVector=[";
    for ( int i = 0; i < apexScanRangeStartVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<apexScanRangeStartVector.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"apexScanRangeEndVector=[";
    for ( int i = 0; i < apexScanRangeEndVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<apexScanRangeEndVector.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"apexWindowRangeStartVector=[";
    for ( int i = 0; i < apexWindowRangeStartVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<apexWindowRangeStartVector.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"apexWindowRangeEndVector=[";
    for ( int i = 0; i < apexWindowRangeEndVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<apexWindowRangeEndVector.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"apexElutionProfileVector={";
    for ( int i = 0; i < apexElutionProfileVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<"[";
        for ( int j = 0; j < apexElutionProfileVector.at(i).size(); j++ ) {
            mzdaLevelTwoOutputStream<<apexElutionProfileVector.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<"];";
    }
    mzdaLevelTwoOutputStream<<"};"<<endl;

    mzdaLevelTwoOutputStream<<"tempidsVector={";
    for ( int i = 0; i < tempidsVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<"[";
        for ( int j = 0; j < tempidsVector.at(i).size(); j++ ) {
            mzdaLevelTwoOutputStream<<tempidsVector.at(i).at(j)+1<<" ";
        }
        mzdaLevelTwoOutputStream<<"];";
    }
    mzdaLevelTwoOutputStream<<"};"<<endl;

}  // END LOG LEVEL CONTROL BLOCK
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// START LOG LEVEL CONTROL BLOCK
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
if ( logLevel == 2 ) {
    ///////////////////////////////////////////////////////////////////////////
    // PEAK APEX DATA - CONDITIONALLY SET DATA STRUCTURES
    ///////////////////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"newApexRtGridIDRangeStartVector=[";
    for ( int i = 0; i < newApexRtGridIDRangeStartVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<newApexRtGridIDRangeStartVector.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"newApexRtGridIDRangeEndVector=[";
    for ( int i = 0; i < newApexRtGridIDRangeEndVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<newApexRtGridIDRangeEndVector.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"peakApexMZVector=[";
    for ( int i = 0; i < peakApexMZVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<peakApexMZVector.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"peakApexmzGridIDVector=[";
    for ( int i = 0; i < peakApexmzGridIDVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<peakApexmzGridIDVector.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"apexElutionProfileShortVector={";
    for ( int i = 0; i < apexElutionProfileShortVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<"[";
        for ( int j = 0; j < apexElutionProfileShortVector.at(i).size(); j++ ) {
            mzdaLevelTwoOutputStream<<apexElutionProfileShortVector.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<"];";
    }
    mzdaLevelTwoOutputStream<<"};"<<endl;

    mzdaLevelTwoOutputStream<<"smoothApexElutionProfileLongVector={";
    for ( int i = 0; i < smoothApexElutionProfileLongVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<"[";
        for ( int j = 0; j < smoothApexElutionProfileLongVector.at(i).size(); j++ ) {
            mzdaLevelTwoOutputStream<<smoothApexElutionProfileLongVector.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<"];";
    }
    mzdaLevelTwoOutputStream<<"};"<<endl;

    mzdaLevelTwoOutputStream<<"sumApexElutionProfileVector=[";
    for ( int i = 0; i < sumApexElutionProfileVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<sumApexElutionProfileVector.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;


    mzdaLevelTwoOutputStream<<"threeSTDVector=[";
    for ( int i = 0; i < threeSTDVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<threeSTDVector.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;


    mzdaLevelTwoOutputStream<<"mzGridIDsVector={";
    for ( int i = 0; i < mzGridIDsVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<"[";
        for ( int j = 0; j < mzGridIDsVector.at(i).size(); j++ ) {
            mzdaLevelTwoOutputStream<<mzGridIDsVector.at(i).at(j)+1<<" ";
        }
        mzdaLevelTwoOutputStream<<"];";
    }
    mzdaLevelTwoOutputStream<<"};"<<endl;


    mzdaLevelTwoOutputStream<<"rcorrMapVector={";
    for ( int i = 0; i < rcorrMapVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<"[";
        for ( int j = 0; j < rcorrMapVector.at(i).size(); j++ ) {
            mzdaLevelTwoOutputStream<<rcorrMapVector.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<"];";
    }
    mzdaLevelTwoOutputStream<<"};"<<endl;


    mzdaLevelTwoOutputStream<<"tempidsRcorrMapVector={";
    for ( int i = 0; i < tempidsRcorrMapVector.size(); i++ ) {
        mzdaLevelTwoOutputStream<<"[";
        for ( int j = 0; j < tempidsRcorrMapVector.at(i).size(); j++ ) {
            mzdaLevelTwoOutputStream<<tempidsRcorrMapVector.at(i).at(j)+1<<" ";
        }
        mzdaLevelTwoOutputStream<<"];";
    }
    mzdaLevelTwoOutputStream<<"};"<<endl;
}  // END LOG LEVEL CONTROL BLOCK
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// START LOG LEVEL CONTROL BLOCK
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
if ( logLevel == 2 ) {
    ///////////////////////////////////////////////////////////////////////////
    // MZDA LEVEL 3 DATA: LCC OUTPUT DATA STRUCTURES
    //
    // - we need to handle some issues here, ...
    //
    ///////////////////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"lccStartMzGridId=[";
    for ( int i = 0; i < lccStartMzGridId.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccStartMzGridId.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    //discrepancy here between Matlab & C++, found on 3/22/13
    //the C++ code appears to be generating local indices instead of global ones
    //resolved on 03/25/13., incorrect indexing
    mzdaLevelTwoOutputStream<<"lccEndMzGridId=[";
    for ( int i = 0; i < lccEndMzGridId.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccEndMzGridId.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    // ISSUE HERE... tbd.... 03/22/13
    mzdaLevelTwoOutputStream<<"lccScanStartID=[";
    for ( int i = 0; i < lccScanStartID.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccScanStartID.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"lccScanEndID=[";
    for ( int i = 0; i < lccScanEndID.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccScanEndID.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    mzdaLevelTwoOutputStream<<"lccStartMz=[";
    for ( int i = 0; i < lccStartMz.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccStartMz.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"lccEndMz=[";
    for ( int i = 0; i < lccEndMz.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccEndMz.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"lccElutionTimeStart=[";
    for ( int i = 0; i < lccElutionTimeStart.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccElutionTimeStart.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"lccElutionTimeEnd=[";
    for ( int i = 0; i < lccElutionTimeEnd.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccElutionTimeEnd.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

}  // END LOG LEVEL CONTROL BLOCK
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// START LOG LEVEL CONTROL BLOCK
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
if ( logLevel == 2 ) {
    ///////////////////////////////////////////////////////////////////////////
    //lccRegionData:
    //A structure of 2-D data
    //
    //test = {[1 2 3; 4 5 6; 7 8 9];[1 2 3; 4 5 6; 7 8 9] }
    ///////////////////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"lccRegionData={";
    for ( int k = 0; k < lccRegionData.size(); k++ ) {
        mzdaLevelTwoOutputStream<<"[";
        for ( int i = 0; i < lccRegionData.at(k).size(); i++ ) {

            for ( int j = 0; j < lccRegionData.at(k).at(i).size(); j++ ) {
                mzdaLevelTwoOutputStream<<lccRegionData.at(k).at(i).at(j)<<" ";
            }  // end  j loop
            mzdaLevelTwoOutputStream<<";";
        } // end i loop
        mzdaLevelTwoOutputStream<<"];"<<endl;
    } // end k loop
    mzdaLevelTwoOutputStream<<"};"<<endl;


    mzdaLevelTwoOutputStream<<"lccSmoothElutionProfile={";
    for ( int i = 0; i < lccSmoothElutionProfile.size(); i++ ) {
        mzdaLevelTwoOutputStream<<"[";
        for ( int j = 0; j < lccSmoothElutionProfile.at(i).size(); j++ ) {
            mzdaLevelTwoOutputStream<<lccSmoothElutionProfile.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<"];";
    }
    mzdaLevelTwoOutputStream<<"};"<<endl;


    mzdaLevelTwoOutputStream<<"lccMzPeakProfile={";
    for ( int i = 0; i < lccMzPeakProfile.size(); i++ ) {
        mzdaLevelTwoOutputStream<<"[";
        for ( int j = 0; j < lccMzPeakProfile.at(i).size(); j++ ) {
            mzdaLevelTwoOutputStream<<lccMzPeakProfile.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<"];";
    }
    mzdaLevelTwoOutputStream<<"};"<<endl;


    mzdaLevelTwoOutputStream<<"lccCenterMzGridId=[";
    for ( int i = 0; i < lccCenterMzGridId.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccCenterMzGridId.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    ///////////////////////////////////////////////////////////////
    // The mz value at the peak
    ///////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"lccCenterMz=[";
    for ( int i = 0; i < lccCenterMz.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccCenterMz.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;



    ///////////////////////////////////////////////////////////////
    // The index into the retention time grid axis at the peak
    ///////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"lccPeakApexRtGridId=[";
    for ( int i = 0; i < lccPeakApexRtGridId.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccPeakApexRtGridId.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    ///////////////////////////////////////////////////////////////
    // The elution time value at the peak
    ///////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"lccPeakApexElutionTime=[";
    for ( int i = 0; i < lccPeakApexElutionTime.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccPeakApexElutionTime.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;
}  // END LOG LEVEL CONTROL BLOCK
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// START LOG LEVEL CONTROL BLOCK
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
if ( logLevel == 1 ) {
    ///////////////////////////////////////////////////////////////////////////
    // Output the input parameters:
    ///////////////////////////////////////////////////////////////////////////
    // output the pre-lc region of interest identifier
    mzdaLevelTwoOutputStream<<"mzdaVersion='1.0.0,LCCAlgV2,Date:03-27-2013';"<<endl;
    mzdaLevelTwoOutputStream<<"precision=10;"<<endl;
    mzdaLevelTwoOutputStream<<"lccAlgorithmVersion=2;"<<endl;
    mzdaLevelTwoOutputStream<<"scanGroupId="<<scanGroupId<<";"<<endl;
    mzdaLevelTwoOutputStream<<"plcRoi="<<plcRoi<<";"<<endl;
    mzdaLevelTwoOutputStream<<"numSmoothPoint="<<numSmoothPoint<<";"<<endl;
    mzdaLevelTwoOutputStream<<"minLCLength="<<minLCLength<<";"<<endl;
    mzdaLevelTwoOutputStream<<"minMZLength="<<minMZLength<<";"<<endl;
    mzdaLevelTwoOutputStream<<"noiseThresholdLevel="<<noiseThresholdLevel<<";"<<endl;
    mzdaLevelTwoOutputStream<<"massResolution="<<massResolution<<";"<<endl;
    mzdaLevelTwoOutputStream<<"LCPeakApexTolerance="<<LCPeakApexTolerance<<";"<<endl;
    mzdaLevelTwoOutputStream<<"rawMZCorrThreshold="<<rawMZCorrThreshold<<";"<<endl;
    // A processing experiment is composed of 1 or more input raw LCMS1 data files( mzXML format )
    mzdaLevelTwoOutputStream<<"sourceFilename='"<<sourceFilename<<"';"<<endl;
    mzdaLevelTwoOutputStream<<"outputFilename='"<<outputFilenamePlcRoi<<"';"<<endl;
    mzdaLevelTwoOutputStream<<"numScans="<<plcRegionID.at(plcRoi).size()<<";"<<endl;
    mzdaLevelTwoOutputStream<<"rtGridSize="<<rtGridSize<<";"<<endl;
    mzdaLevelTwoOutputStream<<"mzGridSize="<<mzGridLength<<";"<<endl;

    ///////////////////////////////////////////////////////////////////////////
    // XICs is the resampled data structure along both the retention
    // time and mz dimensions.
    ///////////////////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"XICs=[";
    for ( int i = 0; i < XICs.size(); i++ ) {
        for ( int j = 0; j < XICs.at(i).size(); j++ ){
            mzdaLevelTwoOutputStream<<XICs.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<";";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    ///////////////////////////////////////////////////////////////////////////
    // The total peak apex count indicates the total number
    // of potential lcc peaks detected.  Some of these are labeled
    // as valid peaks and others are labeled as invalid.
    ///////////////////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"totalPeakApex="<<totalPeakApex<<";"<<endl;


    ///////////////////////////////////////////////////////////////////////////
    // The is a binary mask over the set of peak apex
    // indicating which turned out to be valid ('1') and which
    // did not('0').
    //
    ///////////////////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"lccValidCandidate=[";
    for ( int i = 0; i < lccValidCandidate.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccValidCandidate.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;


    ///////////////////////////////////////////////////////////////////////////
    // These repesent the coordinates of the peak apex within the
    // resampled XICs data structure.
    ///////////////////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"apexWindowIdScanSorted=[";
    for ( int i = 0; i < apexWindowIdScanSorted.size(); i++ ) {
        mzdaLevelTwoOutputStream<<apexWindowIdScanSorted.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"apexScanIdScanSorted=[";
    for ( int i = 0; i < apexScanIdScanSorted.size(); i++ ) {
        mzdaLevelTwoOutputStream<<apexScanIdScanSorted.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    ///////////////////////////////////////////////////////////////////////////
    // These are the resampled retention time and mz axis.
    ///////////////////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"rtGrid=[";
    for ( int i = 0; i < rtGrid.size(); i++ ) {
        mzdaLevelTwoOutputStream<<rtGrid.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"mzGrid=[";
    for ( int i = 0; i < mzGrid.size(); i++ ) {
        mzdaLevelTwoOutputStream<<mzGrid.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    ///////////////////////////////////////////////////////////////////////////
    // These are the indices into the rtGrid of the start and end location
    // of the lcc rt boundary
    ///////////////////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"lccScanStartID=[";
    for ( int i = 0; i < lccScanStartID.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccScanStartID.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"lccScanEndID=[";
    for ( int i = 0; i < lccScanEndID.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccScanEndID.at(i)+1<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    ///////////////////////////////////////////////////////////////////////////
    // These are the values of the starting and ending
    ///////////////////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"lccStartMz=[";
    for ( int i = 0; i < lccStartMz.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccStartMz.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;

    mzdaLevelTwoOutputStream<<"lccEndMz=[";
    for ( int i = 0; i < lccEndMz.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccEndMz.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;


    ///////////////////////////////////////////////////////////////
    // The mz value at the peak
    ///////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"lccCenterMz=[";
    for ( int i = 0; i < lccCenterMz.size(); i++ ) {
        mzdaLevelTwoOutputStream<<lccCenterMz.at(i)<<" ";
    }
    mzdaLevelTwoOutputStream<<"];"<<endl;


    ///////////////////////////////////////////////////////////////
    // The full elution profile
    ///////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"lccSmoothElutionProfile={";
    for ( int i = 0; i < lccSmoothElutionProfile.size(); i++ ) {
        mzdaLevelTwoOutputStream<<"[";
        for ( int j = 0; j < lccSmoothElutionProfile.at(i).size(); j++ ) {
            mzdaLevelTwoOutputStream<<lccSmoothElutionProfile.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<"];";
    }
    mzdaLevelTwoOutputStream<<"};"<<endl;


    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    // output the estimated elution profile
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"lccEstimatedElutionProfile={";
    for ( int i = 0; i < lccEstimatedElutionProfile.size(); i++ ) {
        mzdaLevelTwoOutputStream<<"[";
        for ( int j = 0; j < lccEstimatedElutionProfile.at(i).size(); j++ ) {
            mzdaLevelTwoOutputStream<<lccEstimatedElutionProfile.at(i).at(j)<<" ";
        }
        mzdaLevelTwoOutputStream<<"];";
    }
    mzdaLevelTwoOutputStream<<"};"<<endl;

    ///////////////////////////////////////////////////////////////////////////
    //
    // Output LC Candidate Visualization
    //
    ///////////////////////////////////////////////////////////////////////////
    mzdaLevelTwoOutputStream<<"x=ones(mzGridSize,rtGridSize)*diag(rtGrid);"<<endl;
    mzdaLevelTwoOutputStream<<"y=diag(mzGrid)*ones(mzGridSize,rtGridSize);"<<endl;
    mzdaLevelTwoOutputStream<<"apexScanIdScanSorted = apexScanIdScanSorted';"<<endl;  // transpose required to fix plotting difference between Matlab & C++
    mzdaLevelTwoOutputStream<<"apexWindowIdScanSorted = apexWindowIdScanSorted';"<<endl; // transpose required to fix plotting difference between Matlab & C++
    mzdaLevelTwoOutputStream<<"figure; surf(x,y,XICs);  xlabel('retention time(seconds)'); ylabel('mass/charge ( m/z )'); zlabel('xic intensity');hold on;"<<endl;
    mzdaLevelTwoOutputStream<<"validcounter = 1;"<<endl;
    mzdaLevelTwoOutputStream<<"for peakApexID=1:totalPeakApex"<<endl;
    mzdaLevelTwoOutputStream<<"        if  lccValidCandidate(peakApexID)==1"<<endl;
    mzdaLevelTwoOutputStream<<"        y1= rtGrid(lccScanStartID(peakApexID));"<<endl;
    mzdaLevelTwoOutputStream<<"        y2= rtGrid(lccScanEndID(peakApexID));"<<endl;
    mzdaLevelTwoOutputStream<<"        x1=  lccStartMz(peakApexID);"<<endl;
    mzdaLevelTwoOutputStream<<"        x2=  lccEndMz(peakApexID);"<<endl;
    mzdaLevelTwoOutputStream<<"        plot3(ones(1,2)*y1,[x1 x2],[1 1],'y');"<<endl;
    mzdaLevelTwoOutputStream<<"        plot3(ones(1,2)*y2,[x1 x2],[1 1],'y');"<<endl;
    mzdaLevelTwoOutputStream<<"        plot3([y1 y2],ones(1,2)*x1, [1 1],'y');"<<endl;
    mzdaLevelTwoOutputStream<<"        plot3([y1 y2],ones(1,2)*x2,[1 1],'y');"<<endl;
    mzdaLevelTwoOutputStream<<"        rts=rtGrid(lccScanStartID(peakApexID):lccScanEndID(peakApexID));"<<endl;
    mzdaLevelTwoOutputStream<<"        mzs=ones(lccScanEndID(peakApexID)-lccScanStartID(peakApexID)+1,1)*lccCenterMz(peakApexID);"<<endl;
    mzdaLevelTwoOutputStream<<"        plot3(rts,mzs,lccEstimatedElutionProfile{peakApexID},'r')"<<endl;
    mzdaLevelTwoOutputStream<<"        height= max(lccSmoothElutionProfile{peakApexID})*1.2;"<<endl;
    mzdaLevelTwoOutputStream<<"        stem3(rtGrid(apexScanIdScanSorted(peakApexID)),mzGrid(apexWindowIdScanSorted(peakApexID)),height,'gx');"<<endl;
    mzdaLevelTwoOutputStream<<"        text(rtGrid(apexScanIdScanSorted(peakApexID)),mzGrid(apexWindowIdScanSorted(peakApexID)), height, [num2str(peakApexID)]);"<<endl;
    mzdaLevelTwoOutputStream<<"        validcounter=validcounter+1;"<<endl;
    mzdaLevelTwoOutputStream<<"        end"<<endl;
    mzdaLevelTwoOutputStream<<"end"<<endl;

}  // END LOG LEVEL CONTROL BLOCK
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// START LOG LEVEL CONTROL BLOCK
// MZDA LEVEL 3 DATA
// LOGLEVEL 0 outputs minimal data in level 2 files, minimal
// amount of data needed to create level 3 files.
//
// MZDA LEVEL 3 DATA STRUCTURE
//
//
//
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
if ( logLevel >= 0 ) {

    mzdaLevelTwoOutputStream<<"%Header=centerMz,peakApexElutionTime,mzGridStart,mzGridEnd,mzDelta,rtGridStart,rtGridEnd,rtDelta,mzPeakProfile,estimatedElutionProfile"<<endl;
    for ( int i = 0; i < lccValidCandidate.size(); i++ ) {
        if ( lccValidCandidate.at(i) == true ) {

            ///////////////////////////////////////////////////////////////////
            // MZDA LEVEL 3 DATA PACKET.
            // Per 04/18/13 meeting.
            // 1) lcc_centermz
            // 2) lcc_peakapexelutiontime
            // 3) lcc_mzgridstart
            // 4) lcc_mzgridend
            // 5) lcc_mzgriddelta
            // 6) lcc_rtgridstart
            // 7) lcc_rtgridend
            // 8) lcc_rtgriddelta
            // 9) lcc_mzprofile
            // 10) lcc_estimatedelutionprofile
            ///////////////////////////////////////////////////////////////////
            mzdaLevelTwoOutputStream<<"%LCC"
                                   <<":mC:"<<lccCenterMz.at(i)
                                  <<":rC:"<<rtGrid.at(apexScanIdScanSorted.at(i))
                                 <<":m1:"<<mzGrid.at(lccStartMzGridId.at(i))
                                <<":m2:"<<mzGrid.at(lccEndMzGridId.at(i))
                               <<":"<<mzGridWidth
                              <<":r1:"<<rtGrid.at(lccScanStartID.at(i))
                             <<":r2:"<<rtGrid.at(lccScanEndID.at(i))
                            <<":"<<rtDelta;
            mzdaLevelTwoOutputStream<<":";
            for ( int j = 0; j < lccMzPeakProfile.at(i).size(); j++ ){
                mzdaLevelTwoOutputStream<<lccMzPeakProfile.at(i).at(j)<<",";
            }
            mzdaLevelTwoOutputStream<<":";
            for ( int j = 0; j < lccEstimatedElutionProfile.at(i).size(); j++ ){
                mzdaLevelTwoOutputStream<<lccEstimatedElutionProfile.at(i).at(j)<<",";
            }
            mzdaLevelTwoOutputStream<<":#"<<endl;



            ///////////////////////////////////////////////////////////////////
            //            mzdaLevelTwoOutputStream<<"%LCC:"<<i+1
            //                                   <<":"<<plcRoi
            //                                  <<":mC:"<<lccCenterMz.at(i)
            //                                 <<":rC:"<<rtGrid.at(apexScanIdScanSorted.at(i))
            //                                <<":m1:"<<mzGrid.at(lccStartMzGridId.at(i))
            //                               <<":m2:"<<mzGrid.at(lccEndMzGridId.at(i))
            //                              <<":"<<mzGridWidth
            //                             <<":r1:"<<rtGrid.at(lccScanStartID.at(i))
            //                            <<":r2:"<<rtGrid.at(lccScanEndID.at(i))
            //                           <<":"<<rtDelta;
            ///////////////////////////////////////////////////////////////////
            // The mzGrid needs to be moved to a separate output file
            //            mzdaLevelTwoOutputStream<<":";
            //            for ( int j = 0; j < lccPeakMzGrid.at(i).size(); j++ ){
            //                mzdaLevelTwoOutputStream<<lccPeakMzGrid.at(i).at(j)<<",";
            //            }
            ///////////////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////
            // The rtGrid also need to be moved to a separate output file
            //            mzdaLevelTwoOutputStream<<":";
            //            for ( int j = 0; j < lccPeakRtGrid.at(i).size(); j++ ){
            //                mzdaLevelTwoOutputStream<<lccPeakRtGrid.at(i).at(j)<<",";
            //            }
            ///////////////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////
            //            mzdaLevelTwoOutputStream<<":";
            //            for ( int j = 0; j < lccSmoothElutionProfile.at(i).size(); j++ ){
            //                mzdaLevelTwoOutputStream<<lccSmoothElutionProfile.at(i).at(j)<<",";
            //            }
            ///////////////////////////////////////////////////////////////////
            //            mzdaLevelTwoOutputStream<<":";
            //            for ( int j = 0; j < lccMzPeakProfile.at(i).size(); j++ ){
            //                mzdaLevelTwoOutputStream<<lccMzPeakProfile.at(i).at(j)<<",";
            //            }
            ///////////////////////////////////////////////////////////////////
            //            mzdaLevelTwoOutputStream<<":";
            //            for ( int j = 0; j < lccEstimatedElutionProfile.at(i).size(); j++ ){
            //                mzdaLevelTwoOutputStream<<lccEstimatedElutionProfile.at(i).at(j)<<",";
            //            }
            ///////////////////////////////////////////////////////////////////
            // Outputting the regionData for each inflates the output data too much.
            // We should only output the regionData for debugging&development purposes.
            // Only at LogLevel 2 should we output this information
            ///////////////////////////////////////////////////////////////////
            //            if ( logLevel == 2 ) {
            //                mzdaLevelTwoOutputStream<<":"<<lccRegionData.at(i).size()
            //                                       <<":"<<lccRegionData.at(i).at(0).size()
            //                                      <<":";
            //                for ( int j = 0; j < lccRegionData.at(i).size(); j++ ){
            //                    for ( int k = 0; k < lccRegionData.at(i).at(j).size(); k++) {
            //                        mzdaLevelTwoOutputStream<<lccRegionData.at(i).at(j).at(k)<<",";
            //                    }
            //                }
            //            }
            //            mzdaLevelTwoOutputStream<<":";
            //            mzdaLevelTwoOutputStream<<sourceFilename;
            //            mzdaLevelTwoOutputStream<<":";
            //            mzdaLevelTwoOutputStream<<outputFilenamePlcRoi;
            //            mzdaLevelTwoOutputStream<<":#"<<endl;
            ///////////////////////////////////////////////////////////////////

        }  // end if valid lcc
    } // end loop through all potentially valid peaks



}  // END LOG LEVEL CONTROL BLOCK
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

    mzdaLevelTwoOutputStream.close();



    return 0;
}  // END FUNCTION



} // namespace msda
