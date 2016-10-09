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
/// ScanGroup.h
///
///  Created on: Aug 8, 2011, Last Update: Oct 29, 2012.
///      Author: nelson.ramirez
///  Contains a group of overlapping scan lines
///  These are a subset of the total in an LCMS experiment.
///  A scangroup contains overlapping scanlines to enable
///  parallelism.
///  Each scangroup can be processed independently from all
///  other scangroups.
///  The set of all scangroups contains all the data for
///  an experiment.
///  The scangroup is the fundamental unit of parallelism
///  in this algorithm.
/// Software Developer: Nelson.Ramirez@utsa.edu
///////////////////////////////////////////////////////////////////////////////

#ifndef SCANGROUP_H_
#define SCANGROUP_H_
#include <omp.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include "SignalProcessor.h"
#include "MathUtils.h"
#include "Utilities.h"

using namespace std;

namespace msda {

class ScanGroup {
	
public:
	///////////////////////////////////////////////////////////////////////////
	/// Constructor
	///////////////////////////////////////////////////////////////////////////
	ScanGroup();
	
	///////////////////////////////////////////////////////////////////////////
	/// Destructor
	///////////////////////////////////////////////////////////////////////////
	~ScanGroup();
	
	///////////////////////////////////////////////////////////////////////////
	///
	/// Save the input data filename
	///
	///////////////////////////////////////////////////////////////////////////
	int addInputDataFilename( string & inputDataFilename );

	///////////////////////////////////////////////////////////////////////////
	///
	/// Save the instrument data
	///
	///////////////////////////////////////////////////////////////////////////
	int addInstrumentData(	string & instrManufacturer, 
							string & instrModel, 
							string & instrIonisation, 
							string & instrAnalyzer, 
							string & instrDetector);
	
	///////////////////////////////////////////////////////////////////////////
	///
	/// Save the run header data
	///
	///////////////////////////////////////////////////////////////////////////
	int addRunHeaderData(	int & ScanCount, 
							double & dEndTime, 
							double & dStartTime, 
							double & endMZ,
							double & highMZ,
							double & lowMZ,
							double & startMZ);
	
	///////////////////////////////////////////////////////////////////////////
	///
	/// Save the scan header data
	///
	///////////////////////////////////////////////////////////////////////////
	int addScanHeaderData(		
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
			string & scantype );
	
	///////////////////////////////////////////////////////////////////////////
	/// 
	/// Initially, the scan group object starts off with empty internal
	/// data.  Since all internal data are STL components, this enables
	/// us to dynamically grow the internal scan group data structures
	/// 
	///
	///////////////////////////////////////////////////////////////////////////
	int addScan(vector<double> & mz, vector<double> & intensity, double rt, int scan );
	
	///////////////////////////////////////////////////////////////////////////
	/// Return the total number of scans in the scan group 
	///////////////////////////////////////////////////////////////////////////
	int size(void); //returns the total number of scans in the group
	
	///////////////////////////////////////////////////////////////////////////
	/// Return the minimum retention time in the scan group
	///////////////////////////////////////////////////////////////////////////
	double findMinRetentionTime(void); //returns the minimum retention time value
	
	///////////////////////////////////////////////////////////////////////////
	/// Return the maximum retention time in the scan group
	///////////////////////////////////////////////////////////////////////////
	double findMaxRetentionTime(void); //returns the maximum retention time value

	///////////////////////////////////////////////////////////////////////////
	/// We need complete access into all internal data structures for
	/// visualization purposes, starting from the mzXML global
	/// scan information, the scan header information for each scan,
	/// then, for each scan line that gets taken as input(for example, there
	/// may be subsets of scans that are not read into RAM, depending 
	/// on user selectable criteria such as scan level, etc...
	/// 
	///
	///////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////
	///
	///
	/// Output sorted mz bins with their cutpoints to a file
	///
	///////////////////////////////////////////////////////////////////////////
	int printMzBinningData( string filenamescript,
							string title,
							int scangroupid );
	
	///////////////////////////////////////////////////////////////////////////
	///
	/// Write MSDA Level 1 File
	///
	/// Different versions of the MSDA Level 1 format can be controlled by
	/// provind the version string.
	/// MSDA Level 1 Format
	/// 
	///
    /// Version "1.0.0": Adding a min mz and max mz range to each pre-lc
	///                  region to avoid a repetitive search for 
	/// 				 min mz and max mz ranges.
	///
	///
    ///////////////////////////////////////////////////////////////////////////
    int createMSDALevel1File(
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
            int totalmsLevel1Scans );

	///////////////////////////////////////////////////////////////////////////
	///
	/// Convert parts per million value to its corresponding dalton value(mz)
	///
	///////////////////////////////////////////////////////////////////////////
	double convertPPMtoDaltons(double mzValue, double ppmValue);


    ///////////////////////////////////////////////////////////////////////////
    //
    // Get PLROI( Pre-LC Region Of Interest Properties )
    // Return information about data within the ScanGroup object.
    //
    // Added on 11/1/12 by Nelson.Ramirez@utsa.edu
    ///////////////////////////////////////////////////////////////////////////
    int getNumSignalsPerPlcRoi(vector<int> & numSignalsPerPlcRoi);


	///////////////////////////////////////////////////////////////////////////
	///
    /// Generate a set of PLCROI's( Pre-LC Regions of Interest )
    /// This is  the lowest level of the organization within MZDASoft
    /// framework
    ///
	///////////////////////////////////////////////////////////////////////////
    int findPreLCRegionsOfInterest(double ppmDelta);

	///////////////////////////////////////////////////////////////////////////
	/// Generate 1D signals from each Xic Group.  
	///
	/// Each XIC group consists of a set of centroids meetings a set of 
	/// user specified criteria in both mz and scan(rt) dimensions.
	/// 
	/// In order to process the xic group data as a 1 dimensional signal,
	/// we must merge duplicate items with the same retention time into 
	/// a single data point.
	/// 
	/// The main goal of 
	/// 
	///
	/// 
	///////////////////////////////////////////////////////////////////////////
	int generateXicVectors(void);

	///////////////////////////////////////////////////////////////////////////
	/// 
    /// Level 1 Processing Interfaces
	///
	///
	///////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////
	///
	/// Copy the data buffer from the MsDataLoader class  into the
	/// data buffer in this the ScanGroup class.  This multi-level
	/// data movement strategy from File, to Buffers, to Storage
	/// Management Class ( ScanGroup ) is necessary for speed purposes.
	/// We need to read the data from disk as quickly as possible in large
	/// chucks of data for efficiency, then subsequent processing should
	/// take place within the RAM resident data.
	///
	///
	///////////////////////////////////////////////////////////////////////////
    int addLevel1Data(vector<string>& data);
	
	///////////////////////////////////////////////////////////////////////////
	///
    /// 2 Paths:  Either the Level 1 Data was read in from a file,
	/// or we need to take the data structures that were used to generate
    /// the  Level 1 Data and convert them into a Level 1 in
	/// memory data structure.
    ///
    /// First version of the software will only include option 1.
    /// The level 1 data must be loaded from a file.  Future releases
    /// should include the option to bypass level 1 file generation.
    ///
    ///
	///
	///////////////////////////////////////////////////////////////////////////
    int processLevel1Data(void);

    int lccAlgorithm(string sourceFilename,
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
                     string outputFilename );

    int lccAlgorithmVersion2(   string sourceFilename,
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
                                int logLevel );


	
private:
	
	///////////////////////////////////////////////////////////////////////////
	/// By using the STL data structures and algorithms, we can 
	/// benefit from multi-core research on their parallelization on 
	/// multicore systems.
	///
	///
	/// algo2.iti.kit.edu/singler/mcstl/parallelbulkoperations.pdfSimilar
	/// http://stxxl.sourceforge.net/
	///////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////
	///
	/// RAMP Library Interface data structures
	///
	/// RANDOM ACCESS INTO THE MZXML Data files
	/// The Ramp library has a struct called the RunHeaderStruct
	/// The following is located in mzParser.h
	/// We want to save this information for users to have access to 
	/// all internal data.
	///
	/// Important Goal:  We want to make it easy to expand the data
	/// items that are handled by this software, and so using a data structure
	/// that is dynamic such as a map that takes the name of the data set
	/// feature and generates its value 
	///
	/// It contains the following fields:
	///struct RunHeaderStruct {
	///  int scanCount;
	///	 double	dEndTime;
	///	 double	dStartTime;
	///	 double	endMZ;
	///	 double	highMZ;
	///  double	lowMZ;
	///  double	startMZ;
	///};
    ///
	///
	///typedef struct InstrumentStruct {
	///   char manufacturer[INSTRUMENT_LENGTH];
	///   char model[INSTRUMENT_LENGTH];
	///   char ionisation[INSTRUMENT_LENGTH];
	///   char analyzer[INSTRUMENT_LENGTH];
	///   char detector[INSTRUMENT_LENGTH];
	///} InstrumentStruct;
	///
	///
	/// struct ScanHeaderStruct {
	///   int acquisitionNum; // scan number as declared in File (may be gaps)
	///	  int mergedScan;  /* only if MS level > 1 */
	///   int mergedResultScanNum; /* scan number of the resultant merged scan */
	///   int mergedResultStartScanNum; /* smallest scan number of the scanOrigin for merged scan */
	///   int mergedResultEndScanNum; /* largest scan number of the scanOrigin for merged scan */
	///   int msLevel;
	///	  int numPossibleCharges;
	///   int peaksCount;
	///	  int precursorCharge;  /* only if MS level > 1 */
	///	  int precursorScanNum; /* only if MS level > 1 */
	///	  int scanIndex; //a sequential index for non-sequential scan numbers (1-based)
	///   int seqNum; // number in sequence observed file (1-based)	 
	///   double basePeakIntensity;
	///   double basePeakMZ;
	///   double collisionEnergy;
	///   double highMZ;
	///   double ionisationEnergy;
	///   double lowMZ;
	///	  double precursorIntensity;  /* only if MS level > 1 */
	///	  double compensationVoltage;  /* only if MS level > 1 */
	///   double precursorMZ;  /* only if MS level > 1 */
	///   double retentionTime;        /* in seconds */
	///	  double totIonCurrent;
	///   char activationMethod[SCANTYPE_LENGTH];
	///   char possibleCharges[SCANTYPE_LENGTH];
	///   char scanType[SCANTYPE_LENGTH];
	///   bool possibleChargesArray[CHARGEARRAY_LENGTH]; /* NOTE: does NOT include "precursorCharge" information; only from "possibleCharges" */
	///	 ramp_fileoffset_t	filePosition; /* where in the file is this header? */
	///};
    ///
	///
	///////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////LEVEL 0 DATA////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	vector<string> inputDataFilename;  // Filename of input datasets.  
									   // We use a vector of strings here
									   // because we want to leave the 
									   // possibility open to read data
									   // from multiple files.
	
	///////////////////////////////////////////////////////////////////////////
	/// There is only 1 run header per datasets
	/// Keys:
	/// "scanCount" -->int 		--> runHeaderInts
	/// "dEndTime"  -->double 	--> runHeaderDoubles
	/// "dStartTime"-->double	--> runHeaderDoubles
	/// "endMZ"     -->double	--> runHeaderDoubles
	/// "highMZ"    -->double	--> runHeaderDoubles
	/// "lowMZ"     -->double	--> runHeaderDoubles
	/// "startMZ"   -->double	--> runHeaderDoubles
	///
	/// Add new run header information here
	/// Using a map to store this information allows for flexibility 
	/// when keeping track of the run header information.  We can 
	/// easily add additional header items.
	///
	///////////////////////////////////////////////////////////////////////////
	map< string, string> runHeaderStrings; // string data in run header struct
	map< string, int> runHeaderInts; // integer data in run header struct
	map< string, double> runHeaderDoubles; // double data in run header struct
	
	///////////////////////////////////////////////////////////////////////////
	/// There is only 1 instrument struct per dataset
	/// Keys:
	/// "manufacturer" 	--> string --> instrumentInfoStrings
	/// "model" 		--> string --> instrumentInfoStrings
	/// "ionisation"   	--> string --> instrumentInfoStrings
	/// "analyzer"      --> string --> instrumentInfoStrings
	/// "detector"    	--> string --> instrumentInfoStrings
	///
	/// Add new instrument information here
	///
	///
	///////////////////////////////////////////////////////////////////////////
	map< string, string> instrumentInfoStrings; // string data in instrument struct
	map< string, int> instrumentInfoInts; // integer data in instrument struct
	map< string, double> instrumentInfoDouble; // double data in instrument struct
	
	///////////////////////////////////////////////////////////////////////////
	///
	/// There is 1 scan header struct per scan
	/// For debug and informational purposes, we want to keep all the 
	/// data related to the range of scans the user selected to be 
	/// processed by this job.
	/// Since the user given control over which scans to select for analysis
	/// based on a set of criteria, it is important to know what the data 
	/// in thos scans that were AND were not selected was in order to be able
	/// to know why something was or was not included in the analysis dataset.
	///
	/// While we may not load the data in a scan, we still want to keep 
	/// a record of its scan header so that we can go back to it for 
	/// logging and system verification purposes.
	///
	///////////////////////////////////////////////////////////////////////////
	vector< map< string, string > >originalDataScanHeaderStrings; // string data in scan header structs
	vector< map< string, int > >originalDataScanHeaderInts; // integer data in scan header structs
	vector< map< string, double > >originalDataScanHeaderDoubles; // double data in scan header structs
	vector< map< string, bool> >originalDataScanHeaderBools; // boolean data in scan header structs
	
	///////////////////////////////////////////////////////////////////////////
	///
	/// \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
	///               Scans are selected based on user criteria
	///                        scan level, count > 0
	/// \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
	///
	/// Filtering based on user criteria.
	/// Current criteria specified in the MsDataLoader class.
	/// The MsDataLoader class populates a ScanGroup object.
	/// Within the MsDataLoader class, user criteria is applied as to which
    /// scans are loaded into memory, and what subset of values are 
	/// loaded into memory(i.e., only level 1 scans, which contain 
	/// at least a certain number of data items which have mz and intensity
	/// values in a user specified range. 
	///
	///
	///
	///////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////
	///////////////////////////LEVEL 1 DATA////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	
	
	///////////////////////////////////////////////////////////////////////////
	/// Notes: Level 1 data is a filtered subset of the data in the original
	///        mzXML file.
	///        For example, the user can select to only load the level 1( ms )
	///        scans.  So these would be the only ones that get loaded into
	///        RAM.
	///        Also, 
	///
	/// Selected Scans and Data Items: Raw Data
	/// 
	///
	/// sg.addScan(mzdata, intensitydata, scanHeader.retentionTime, scanHeader.acquisitionNum );
	///
	/// This is all data required for a scan group 
	/// We use a vector of vectors to represent a scan group in memory.
	/// Using a list of vectors was investigated, but we need random
	/// access to any specific scan line.
	///
	/// High frequency data has many 0 intensity data values.
	/// These are not stored in the scangroup to conserve memory, since 
	/// there is no need for keeping data values with intensity of zero.
	/// ( Z.Wang/M.Zhang refernce: 10/3/11 meeting)
	///
	///
	///
	///////////////////////////////////////////////////////////////////////////
	vector< vector<double> > mzData; // the intensity data of a single scan
	vector< vector<double> > intensityData; // the intensity data of a single scan
	vector<double> retentionTime; // maps to scanHeader.retentionTime
	vector<int> scanNum; // maps to scanHeader.acquisitionNumber
	
	///////////////////////////////////////////////////////////////////////////
	// The original data may need to be cleared to conserve memory.
	//
	///////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////
	///
	/// The originalScanGroupIndex allows us to go back to the original
	/// scan header information to be able to see [WHY] we included one 
	/// scan and excluded another scan.
	/// For example, when including only ms level 1 scans, we need to 
	/// be able to easily verify whether a particular scan truly was level 2
	/// data.
	///
	///////////////////////////////////////////////////////////////////////////
	vector<int> originalScanGroupIndex; // An index into the originalDataScanHeader 
										// data structures above 
										// maps the mzData,intensityData,retentionTime,
										// scanNum vectors back to the original 
										// scan headers.
	
	///////////////////////////////////////////////////////////////////////////
	// 
	// 
	// We need a way to create a dynamic set of histogram bins
	// For each centroid we need to create a set of bins
	// centroid # --> [start end]  --> histBinStart[#] histBinEnd[#]
	// centroid # --> [start end]  --> histBinStart[#] histBinEnd[#]
	// centroid # --> [start end]  --> histBinStart[#] histBinEnd[#]
	// 
	// The index into the dynamic histogram is used to create a 
	// map data structure.
	// 
	// 
	///////////////////////////////////////////////////////////////////////////
	vector< double > mzHistBinStart;  //These are unordered bins-> index = data id (mz,int)
	vector< double > mzHistBinEnd;    //These are unordered bins-> index = data id (mz,int)
	vector<int> mzHistBinSortIndex;   //This contains the ordered peak list according to mz value
	
	///////////////////////////////////////////////////////////////////////////
	//
	//  Pre-LC Region Data Structures
	//  Indexed by the xic candidate id
	//  
	//   A Pre-LC Region consists of a set of Pre-LC Region Components.
	//   
	///////////////////////////////////////////////////////////////////////////
	vector<int> preLcCandidateBoundaries; // index into mzHistBinSortIndex
	vector< map<int,int> > preLcCandidateScanVsItemCountMap; // candidate 1--> map[scan#] --> total number of data elements in this scan
	vector< map<int,double> > preLcCandidateScanVsIntensitySumMap; // For each pre lc candidate, this is the sum of the intensities at a given scan  
 	vector< map<int,double> > preLcCandidateScanVsAverageMzMap; // For each pre lc candidate, this is the average of the intensities at a given scan
 	vector< multimap<int,double> > preLcCandidateScanMzValues;  // For each pre lc candidate, these are all the mz values at a given scan
 	vector< multimap<int,double> > preLcCandidateScanIntensityValues; // For each pre lc candidate, these are all the intensity values at a given scan
	vector< map<int,double> > preLcCandidateScanVsRetentionTimeMap; // For each pre lc candidate, this is the retention time at a given scan
	vector< map<int,double> > preLcCandidateScanVsMinMzMap; // For each pre lc candidate, this is the minimum mz value at a given scan
	vector< map<int,double> > preLcCandidateScanVsMaxMzMap; // For each pre lc candidate, this is the maximum mz value at a given scan

	
	///////////////////////////////////////////////////////////////////////////
	//
	//
	// The Pre-LC Region Data structures above are used to generate 
    // a Level 1 file.
	//
	//
	///////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////
    // A Level 1 file has the following structure:
	//
	// The preLCRegion is identified by n, and the component signal
	// of a specific preLCRegion is identified by m.
	// 
	//
	//preLCRegions.header, 
	//preLCRegions.scannum --> preLCRegions.scannum{n} --> [component(m)]
	//preLCRegions.retentiontime --> preLCRegions.retentiontime{n} --> [component(m)]
	//preLCRegions.plcid --> preLCRegions.plcid{n} --> [component(m)]
	//preLCRegions.intensitysum --> preLCRegions.intensitysum{n} --> [component(m)] 
	//preLCRegions.averagemz --> preLCRegions.averagemz{n} --> [component(m)]
	//preLCRegions.intensitydata --> preLCRegions.intensitydata{n} --> [component(m)]
	//preLCRegions.mzdata --> preLCRegions.mzdata{n} --> [component(m)]
	///////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////
	//
    // Level 1 Format internal data structures
	//
	//
	//    There are 2 ways to continue the processing stages.
	//    1) We can continue using the Pre-LC Region Data structures to
	//       generate the LC candidates, OR
    //    2) We can read an MZDALevel1 Data format file from disk, and
	//       populate a set of data structures.
	//
	//  A few important items to note are:
	//  1) We may want to create a PreLcRegions Class to contain 
	//     all data related to a PreLcRegion.
	//  2) We should provide a way to create a PreLcRegions Class both
	//     from the internal data structures located within 
	//     the ScanGroup class and from a file.
	//     
	//
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	// Buffer data structures used to speed up the input of the 
	// data from files.
	///////////////////////////////////////////////////////////////////////////
    vector<string> level1TextData;  // Text data for msda level 1 general structure(text)
    vector<string> level1BinaryData;  // Binary data as read from the Level 1 file(binary)
	///////////////////////////////////////////////////////////////////////////
    // Level 1 File In Memory Representation:
	///////////////////////////////////////////////////////////////////////////
    vector< string > level1HeaderStringData; // All string header data
    vector< int > level1HeaderIntData;  // All integer header data
    vector< double > level1HeaderDoubleData; // All double header data
	

    ///////////////////////////////////////////////////////////////////////////
    //
    // For each pre-LC region, there are a set of components, each pre-LC
    // component consists of the following data items:
    // mz,intensity,[scan#/retentiontime]
    // mz,intensity,[scan#/retentiontime]
    //
    ///////////////////////////////////////////////////////////////////////////




	///////////////////////////////////////////////////////////////////////////
	// 
	// For each pre-LC region, there are a set of components, each pre-LC
	// component consists of the following data items:
	// mz,intensity,[scan#/retentiontime]
	// mz,intensity,[scan#/retentiontime]
	// 
	// Each pre-LC component has a specific scan# / retention time
	// associated with it.
	// 
	// These data structures are populated when reading an MSDA Level 1
	// Parallelism notes:  These lc Candidate data structures must be
	// read-only after they are initially created.  It is particularly 
	// important to keep these data structures are read only during the 
	// pre-lc candidate generation stage, where OpenMP is used 
	// to process each pre LC region of interest in parallel.
	//
    //
	//
	///////////////////////////////////////////////////////////////////////////	
	vector< vector<int> >plcRegionID;  // the list of plc region IDs for each component within a pre-LC region
	vector< vector<int> >plcScanNumber; // the list of scan numbers for each component within a pre-LC region
	vector<	vector<double> > plcRetentionTime; // the list of retention times for each component within a pre-LC region
	vector< vector<int> > plcCount;  // this will tell us how many data points in a component
	vector< vector< vector<double> > >plcMzData; // this contains the actual data for all the components within a single preLC region
	vector< vector< vector<double> > >plcIntensityData; // this contains the actual data for all the components within a single preLC region
	

	
	///////////////////////////////////////////////////////////////////////////
	//
	//
	//
	// The next level of data is the LC Candidate.
	// The first stage is generating a set of Pre-LC Regions with 
	// a data structure that is friendly to signal processing routines
	// such as smoothing, histogramming, peak picking.
	//
	// Each Pre-LC region is made up of a number of component
	// signals.  When all the signals are collapsed and viewed
	// as a summation, they look as depicted below:
	//
	// The issue that arises here is that if there are too 
	// many components, the merged XIC signal becomes too 
	// complex to analyze.
	// 
	// This is where processing within the the pre-LC region
	// comes in.  Each pre-LC region is further cut into
	// LC candidates via peak and interval detection.
	// These LC candidates are used to generate Peptide
	// candidates.  LC candidates and Peptide candidates
	// are the precursors to further MS processing stages.
	// 
	//
	//
	//                   XIC Profile of a PRE-LC Region of Interest
	// This view represents the collapsing of all PRE-LC Region Component
	// signals: 
	//
	//
	//
	// ^intensity
	// |                        
	// |        
	// |         +       +                    +     			      
	// |         ++      +                    + +  +++				
	// |       +++++++  +++++               ++++++++++++           mz delta
	// |    +++++++++++++++++++         +++++++++++++++++++++      
	// |++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++>
 	//														 mz delta
	//
	// The LC Peak candidates are used to generate the set of peptide
	// candidates.
	//
	//
	// Both the LC Peak candidates and Peptide candidates are the foundation
	// for a multitude of additional processing stages, such as alignment,
    // database searching, quantitation.
	//
	///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    //
    // Each Scan Group is composed of a set of PLCROI(PreLC ROI),
    // and a PLCROI contains a set of LC Candidates.
    //
    ///////////////////////////////////////////////////////////////////////////
    // Version 1.0 uses file based I/O.  Future version may use
    // fully internal data structures without outputting multiple files,
    // for improved efficiency.
    ///////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////
	// MSDA Level 2 File Format
	// Data structures to enable the display and further processing of the 
	// data.
    // Level 2 Header
    // Level 1 Header
    // Summary data for the mz,rt ranges within this file(useful to query)
    // LC Candidate Related Data
    //  - XIC
    //  - Resampled Region
    //  - ....
    //  - LC Candidate Specific Information
    //
	///////////////////////////////////////////////////////////////////////////
    // Version 1.0 uses file based I/O.  Future version may use
    // fully internal data structures without outputting multiple files,
    // for improved efficiency.
    ///////////////////////////////////////////////////////////////////////////

	
};

} //namespace msda
#endif // SCANGROUP_H_
