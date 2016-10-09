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
// msda.cpp
// 
// Part of MZDASoft
// 
// A framework for high performance mass spectrometry data analysis.
//
// 
// version 1.0.0
//
// Change History:
// Feature ID, Version, Date, Description, @[ADD,CHG,MOV]DeveloperID
// ft00000, vr1.0.0, 04/11/12, Initial Version, @ADDNRZ
//
//
// Developer's Contact Information:
// Developer ID, Name, Email
// NRZ, Nelson Ramirez, Nelson.Ramirez@utsa.edu
//
// Language: C++
//
// Documentation: Doxygen
// 
// Key Coding Suggested Conventions:
// Refer to CodingGuidelines.txt file
//
// External Libraries Used: 
// 1) mstoolkit, Ramp
//    license = BSD, link=  http://code.google.com/p/mstoolkit/
//
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstring>
#include <ctime>
#include <vector>


///////////////////////////////////////////////////////////////////////////////
//////////////////////Parallelism Libraries////////////////////////////////////
//////////////////////OpenMPI,OpenMP,CUDA,OpenCL///////////////////////////////
//////////////////////Note: For the first version of this software/////////////
//////////////////////we plan to use only the OpenMP libraries/////////////////
//////////////////////Future versions may make use of both MPI, and ///////////
//////////////////////GPU Acceleartion approaches, with////////////////////////
//////////////////////Thrust, Cusp C++ libraries///////////////////////////////
///////////////////////////////////////////////////////////////////////////////
#ifdef OPENMP
#include <omp.h>  // Conditional compiling of omp.h header
#endif

///////////////////////////////////////////////////////////////////////////////
// User classes
///////////////////////////////////////////////////////////////////////////////
#include "msda.h"  // Uses general functionality from cbi hpc library

///////////////////////////////////////////////////////////////////////////////
// CBI Library
///////////////////////////////////////////////////////////////////////////////
#include "cbilib.h"

///////////////////////////////////////////////////////////////////////////////
/////////////////////////Debug Libraries///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
#include <signal.h>
///////////////////////////////////////////////////////////////////////////////
/////
///// Signal handling functionality.  We need to be able to catch
///// signals and then attach to the process dynamically.  This is meant to
///// allow runtime debugging of the program within a debug enabled 
///// environment.
/////
///// Usage example in mainline code: 
///// setup_sigsegv_catcher();
///// When the troublesome code executes:
///// char *a = NULL; *a=0;  --> generates a SIGSEGV signal -->
/////  Operating system calls the previously registered signal handler,
/////  which then places the process into a waiting loop while we 
/////  can attach the debugger to the process.  This will let us
/////  get complete access to the running process and all its internal
/////  structures at the point of failure.
///// 
///////////////////////////////////////////////////////////////////////////////	
#ifdef MSDADEBUG
#ifdef MSDADEBUGWITHSIGNALS
static void sigsegv_handler(int);
static char message[1000];
static int messagelength;
void setup_sigsegv_catcher(void) {
	 cout<<"testing setup_catcher"<<endl;
     sprintf(message, "Interrupt signal: pid %d\n", getpid());
     messagelength = strlen(message);
     signal(SIGSEGV, sigsegv_handler);
 }
 static void sigsegv_handler(int signal) {
     write(1, message, messagelength);
     cout<<"testing sighandler"<<endl;
     //////////////////////////////////////////////////////////////////////////
     // 
     // Add stack trace printing 
     // - Printing the stack trace is very platform specific.
     //////////////////////////////////////////////////////////////////////////
   
     while(1);  // Wait while we can attach a debugger to the process 
     exit(0);
 }
#endif
#endif
///////////////////////////////////////////////////////////////////////////////
/////////////////////////End Debug Libraries///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

using namespace std;

///////////////////////////////////////////////////////////////////////////////
//  START HERE
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {
	cout<<"\n****************** + ** + ************"<<endl;
	cout<<"************Version 1.0 *** + ********"<<endl;
	cout<<"********************* +  + +  +*******"<<endl;
	cout<<"***************************+ + +******"<<endl;
	cout<<"**************************************"<<endl;
        cout<<"************MZDASoftPPE***************"<<endl;
        cout<<"*******Mass Spectrometry Data*********"<<endl;
	cout<<"*********Analysis Software************"<<endl;
	cout<<"**************************************"<<endl;
	cout<<"**************************************"<<endl;
	cout<<"Copyright (c) 2014, The University of"<<endl;
	cout<<"Texas at San Antonio. All rights"<<endl;
	cout<<"reserved."<<endl;
	cout<<"MZDASoftPPE v.1.0 is under the"<<endl;
	cout<<"GPL Version 2 License."<<endl;
	cout<<"Please refer to the"<<endl;
	cout<<"MZDASoftPPELicense.txt file containing"<<endl;
	cout<<"the full license text."<<endl;
	cout<<"**************************************\n"<<endl;

	///////////////////////////////////////////////////////////////////////////
	//
	// Output system information when in debug mode
	//
	///////////////////////////////////////////////////////////////////////////
#ifdef MSDADEBUG
	CBILIB_PRINTLOCATION;
	CBILIB_PRINTDATETIME;
#endif
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Program wide objects
	// The software process is controlled by 2 parameter files
	// The information in these files completely identify all characteristics
	// of the processing.
	// The global parameter file has information that pertains to all workers.
	// The local worker parameter file has information that pertains to a
	//  single worker, this varies depending on the task that  
	//
	//
	///////////////////////////////////////////////////////////////////////////
	cbi::CommandLineProcessor clProcessor; // Command line processor
	cbi::ParameterLoader paraObjGlobal;  // Global parameters
	cbi::ParameterLoader paraObjLocalWorker; // Parameters for each local worker
    cbi::ParameterLoader paraObjPartitions; // Partitions file parameters
    cbi::SignalProcessor spObject;  // Signal processing object
    msda::Utilities utilObject; // Utilities object
    msda::MsDataLoader msdlObject; // Instantiate MsDataLoader Object
    msda::ScanGroup sgObject; // Scan group object
    msda::Workflow wfObject; // Workflow object
	
	
	int rc=0;  // Return code variable.  Contains return codes from all 
	           // function calls within the main program.
	///////////////////////////////////////////////////////////////////////////
	// When running in parallel mode, we need to perform a few runtime validity 
	// checks to ensure that 
	///////////////////////////////////////////////////////////////////////////
	try {
		
	
	
		///////////////////////////////////////////////////////////////////////
		// End parallel runtime checks.
		// When we are running in either mpi and /or openmp mode, we want
		// to perform some sanity checks at the start of the program to make
		// sure that our parallel environment produces correct results for a
		// simple test.
		///////////////////////////////////////////////////////////////////////
		
		///////////////////////////////////////////////////////////////////////
		// Start of code
		// Notes:
		// Comment by Z.Wang.  Allow user to run in multiple modes,
		// - single system local ( A single node with multiple cores )
		// - grid ( Grid mode with Sun Grid Engine )
		///////////////////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////////////////
		//
		// Handle processing the input parameters
		//
		//  - The MSDA/MzDa Software has the following main modes of 
		//    operation.
        //  - Mode 1: Generates the set of scripts and distributed parameter
        //    files necessary to implement the entire workflow end-to-end, in
        //    either the single multi-core node mode using PRLL script
        //    based queueing system or grid multi-node mode using Sun
        //    Grid Engine.
		//
		//
		///////////////////////////////////////////////////////////////////////
		// 
		//  Mode 1) Work distribution management( WorkloadManager Object)
		//          The WorkloadManger object is a framework from the CBI
		// 			HPC infrastucture which provides an infrastructure
		// 			for general distributed work management on one 
		//			or more computational nodes.  The workload manager is 
		// 			a general framework for allowing an application to 
		//			provide information as to the general runtime 
		//			characteristics of the hardware and software via 
		//			configuration files and the WorkloadManager will 
		//			streamline the process of automatically creating the 
		//			set of files needed to run a highly independent set
		//			of workunits on either a single multi-core system,
		//			many dedicated systems, or a grid such as the 
		//			sun grid engine.  The application must be aware
		// 			of both hardware, software, command line parameters,
		//			dependencies, etc, to coordinate a massively 
		//			independent set of work units.  Many times, processing
		//			needs to be done in stages, where each stage has 
		//			different parallelism requirements and capabilites.
		//			These types of applications further demonstrate the 
		//			need for coordination of the work from within the 
		//			application itself via a helper framework such as the 
		//			WorkloadManager framework from the CBI High Performance
		//			Enablement runtime library.
		//			The management of hundreds of independent portions
		// 			of work with multiple stages of work requires a 
		//			workflow management framework. 
		//
		//		 1.0) Load the global parameter file
		// 		 1.1) WorkloadManager: Create a set of local parameter files.  These
		//			  define what each distributed worker needs to be working
		// 			  on.
		//		 1.2) WorkloadManager Create a set of shell scripts that will handle
		// 			  the distributed nature of the data binning and 
		// 			  low level data filtering process.
		//			  Depending on the computational resources specified,
		//   		  the scripts must be tailored to run only a specific
		//            number of processes at any one time.  The approach
		// 			  taken is as follows. 
		//				- Version 1 provides an interface to the PRLL library 
		//				which enables the controlled execution of
		//				batch jobs on a single multi-core system.
		//				- Since the PRLL Library is GPL licensed, we 
        //				must only use it via well defined interfaces. This
        //
		//  	 1.3) Create a set of job scripts that will handle 
        // 			  the distributed nature of the data binning
		//			  and low level data filtering process via a grid 
		//			  engine such as the Sun Grid engine.
		//			  Use job groups to group the work into logical units.
		//			  This will be a next release work item.
		//          
		//  -m 1 ( mode )
        //
        //
		//  -g <globalparameterfilename> ( use absolute path )
		//     ( All user defined criteria is controlled through the 
		// 		 global parameter file. )
		//     This is the only file the user must provide the absolute
		//     path to.
		//     Since the GUI has the global parameter file, it will know 
		//     the full path and filename of the resulting run script.
		//     Once the run script processing is complete, the user/
		//	   gui can then run the dynamically generated script.
		//  
		//
		//
		///////////////////////////////////////////////////////////////////////
		//  Mode 2) Distributed data I/O and data organization into a set
		//          of MSDA Level 1 format files.
		///////////////////////////////////////////////////////////////////////
		//  Mode Flag:
		//  -m 2 : 
		//  
		//  Global Parameter Filename Flag:
		//  -g <filename>: Global paramter file absolute path and filename
		//
		//  Mode 2)
		// 
		//   
        // .. add more details on LC Candidate generation phase...
		//
		//
		//
		//
		//
		//
		//
		///////////////////////////////////////////////////////////////////////
		
#ifdef MSDADEBUG
		cout<<"MZDA COMMANDLINE INTERFACE"<<endl;
#endif
		///////////////////////////////////////////////////////////////////////
		// Load the command line data
		///////////////////////////////////////////////////////////////////////
		rc = clProcessor.loadCommandLineData(argc, argv);
		vector<int> modeParameter; // needed by all run modes
		vector<string> globalParameterFilenameParameter; // needed by all run
														 // modes
		vector<string> localParameterFilenameParameter;  // The local 
														 // parameter file
														 // is tailored to the 
														 // needs of the 
														 // specific work 
														 // unit.
        vector<string> partitionsFilenameParameter; // Contains a list of all partitions files
                                                    // available.  This allows the dynamic
                                                    // processing of any number of input
                                                    // mzXML files.
        vector<string> outputMergedPartitionParameters;  // Contains the results of merging all
                                                // partition information into a single output file.
                                                // this file can be copied into the global parameter
                                                // file for MZDASoft system.
        vector<string> outputClusterJobsFilenameBase; // Contains the base filename that will be used to
                                                      // generate the set of cluster job scripts that
                                                      // will manage a large scale multi-file processing
                                                      // run.

		rc = clProcessor.getParameterValuesAsInt("-m",modeParameter);
		rc = clProcessor.getParameterValuesAsString("-gp", globalParameterFilenameParameter);
		rc = clProcessor.getParameterValuesAsString("-lp", localParameterFilenameParameter);
        rc = clProcessor.getParameterValuesAsString("-partitions",partitionsFilenameParameter);
        rc = clProcessor.getParameterValuesAsString("-outputmergedparms",outputMergedPartitionParameters);
        rc = clProcessor.getParameterValuesAsString("-outputclusterjobsfilenamebase",outputClusterJobsFilenameBase);
		///////////////////////////////////////////////////////////////////////
		// Make sure mode parameter has been specified correctly.
		///////////////////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////////////////
		// Make sure mode parameter has been specified correctly.
		///////////////////////////////////////////////////////////////////////
		if ( modeParameter.size() == 1 ) {
			///////////////////////////////////////////////////////////////////
			// check for modes of operation
			///////////////////////////////////////////////////////////////////
			

            ///////////////////////////////////////////////////////////////////
            //
            // Workload management mode( workflow #0 )
            //
            // -m 0,mode 0 indicates we are in startup mode, either being
            // called directly from the command line, or being called
            // by another software component such as the GUI, this
            // workflow generates all the scripts that will be necessary
            // to perform the end-to-end workflow according to
            // what is specified in the global parameter file.
            //
            // We need to provide a framework that will allow the creation
            // of a variety of different distributed workflows in a way
            // that will be relatively straightforward to extend.
            //
            //
            //
            //
            ///////////////////////////////////////////////////////////////////
            if ( modeParameter.at(0) == 0 ) {
                cout<<"Workload Management Mode(-m 0)"<<endl;

                ///////////////////////////////////////////////////////////////
                // load global parameter file.
                ///////////////////////////////////////////////////////////////
                if( globalParameterFilenameParameter.size() > 1 ) {

                    throw "Error: More than one global parameter file specified";
                }

                ///////////////////////////////////////////////////////////////
                // make sure that a global parameter file was specified
                ///////////////////////////////////////////////////////////////
                if( globalParameterFilenameParameter.size() == 0 ) {

                    throw "Error: No global parameter file specified";
                }

                ///////////////////////////////////////////////////////////////
                // try to load the global parameter file
                ///////////////////////////////////////////////////////////////
                rc = paraObjGlobal.loadParameterDataFromFile(globalParameterFilenameParameter.at(0));  // First argument contains parameter filename
                cout<<"Load of Global Parameter File succeeded, rc = "<<rc<<endl;


                ///////////////////////////////////////////////////////////////
                //
                // Run workflow # 0
                // - After running this workflow, we should have
                //
                ///////////////////////////////////////////////////////////////
                rc = wfObject.workflowZero(paraObjGlobal);
                if( rc == 0 ) {
                    cout<<"Work flow # 0 succeeded"<<endl;
                }
            }  // End workflow #0








			///////////////////////////////////////////////////////////////////
			//
			// Workload management mode( workflow #1 )
			// 
			// -m 1,mode 1 indicates we are in WorkloadManagement mode
			// we need to load the global parameter file and prepare
			// a set of dynamically generated files and scripts that
			// will perform the end-to-end computational workflow
			// from the operating system command line. Depending on the 
			// hardware information available, the scripts are dynamically
			// generated to organize the distributed nature of the 
			// different stages of the computational workload.
			//
			//
			///////////////////////////////////////////////////////////////////
			if ( modeParameter.at(0) == 1 ) {
				cout<<"Workload Management Mode(-m 1)"<<endl;
				
				///////////////////////////////////////////////////////////////
				// load global parameter file.
				///////////////////////////////////////////////////////////////
				if( globalParameterFilenameParameter.size() > 1 ) {

					throw "Error: More than one global parameter file specified";
				}
				
				///////////////////////////////////////////////////////////////
                // make sure that a global parameter file was specified
				///////////////////////////////////////////////////////////////
				if( globalParameterFilenameParameter.size() == 0 ) {

					throw "Error: No global parameter file specified";
				}

				///////////////////////////////////////////////////////////////
				// try to load the global parameter file
				///////////////////////////////////////////////////////////////
				rc = paraObjGlobal.loadParameterDataFromFile(globalParameterFilenameParameter.at(0));  // First argument contains parameter filename
				cout<<"Load of Global Parameter File succeeded, rc = "<<rc<<endl;
				
				
				///////////////////////////////////////////////////////////////
				//
                // Run workflow # 1
				//
				//
				///////////////////////////////////////////////////////////////
				rc = wfObject.workflowOne(paraObjGlobal);
				if( rc == 0 ) {
					cout<<"Work flow # 1 succeeded"<<endl;
				}
            }  // End workflow #1
			
			///////////////////////////////////////////////////////////////////
			//
			// The only mode that needs to be called directly by the user
			// is mode 1.  Mode 1 generates a set of scripts that 
			// automatically manage the rest of the run.  The rest of the 
			// modes should normally be called only from within the 
			// main runManager script.
			//
			//
			///////////////////////////////////////////////////////////////////
			
			
			///////////////////////////////////////////////////////////////////
			//
			// Data Binning & Organization Mode Only( mode #2 )
			// 
			// -m 2,mode 2 indicates we are running as a localized worker
			// that needs to load a portion of a data file.
			// MSDA Level 1 format output file.
			//
			///////////////////////////////////////////////////////////////////
			if ( modeParameter.at(0) == 2) {
				///////////////////////////////////////////////////////////////
				// output msda level 1 format data
				///////////////////////////////////////////////////////////////
				cout<<"Data Binning & Organization Mode Only(-m 2)"<<endl;

				///////////////////////////////////////////////////////////////
				// check global parameter file flag data
				///////////////////////////////////////////////////////////////
				if( globalParameterFilenameParameter.size() > 1 ) {

					throw "Error: More than one global parameter file specified";
				}
				if( globalParameterFilenameParameter.size() == 0 ) {

					throw "Error: No global parameter file specified";
				}

				/////////////////////////////////////////////////////////////// 
				// try to read the global parameter file
				///////////////////////////////////////////////////////////////
				rc = paraObjGlobal.loadParameterDataFromFile(globalParameterFilenameParameter.at(0));  // First argument contains parameter filename
				cout<<"Load of Global Parameter File succeeded, rc = "<<rc<<endl;

				///////////////////////////////////////////////////////////////
				// load local parameter file
				///////////////////////////////////////////////////////////////
				if( localParameterFilenameParameter.size() > 1 ) {

					throw "Error: More than one local parameter file specified";
				}
				if( localParameterFilenameParameter.size() == 0 ) {

					throw "Error: No global local parameter file specified";
				}

				///////////////////////////////////////////////////////////////
				// try to load the local parameter file
				///////////////////////////////////////////////////////////////
                cout<<"localParameterFilename="<<localParameterFilenameParameter.at(0)<<endl;
				rc = paraObjLocalWorker.loadParameterDataFromFile(localParameterFilenameParameter.at(0));  // First argument contains parameter filename
				cout<<"Load of Local Parameter File succeeded, rc = "<<rc<<endl;

				///////////////////////////////////////////////////////////////
				//
				// worflow # 2
				//
				//
				///////////////////////////////////////////////////////////////
                rc = wfObject.workflowTwo(paraObjGlobal,paraObjLocalWorker,sgObject,msdlObject,spObject);


				if( rc == 0 ) {
					cout<<"Work flow # 2 succeeded"<<endl;
				}	
			} // end workflow 2 
			
			
			
			///////////////////////////////////////////////////////////////////
			//
            // LC Candidate Generation ( Creating Level2 output files) ( workflow #3 )
            // - Use Level 1 file as input( LC Regions of Interest )
            // - Produce Level 2 output files ( LC Candidates )
            //
			// -m 3, mode 3 indicates we are running as a localized worker
            // that needs to load a Level 1 file, and run the
            // lc candidate generation algorithm to find the
            // set of lc candidates contained within the regions of interest
            // in the level 1 file.
			//
			//
			//
			///////////////////////////////////////////////////////////////////
			if ( modeParameter.at(0) == 3) {

                cout<<"LC Candidate Generation ( workflow #3 )(-m 3)"<<endl;


                ///////////////////////////////////////////////////////////////
                // check global parameter file flag data
                ///////////////////////////////////////////////////////////////
                if( globalParameterFilenameParameter.size() > 1 ) {

                    throw "Error: More than one global parameter file specified";
                }
                if( globalParameterFilenameParameter.size() == 0 ) {

                    throw "Error: No global parameter file specified";
                }

                ///////////////////////////////////////////////////////////////
                // try to read the global parameter file
                ///////////////////////////////////////////////////////////////
                rc = paraObjGlobal.loadParameterDataFromFile(globalParameterFilenameParameter.at(0));  // First argument contains parameter filename
                cout<<"Load of Global Parameter File succeeded, rc = "<<rc<<endl;

                ///////////////////////////////////////////////////////////////
                // load local parameter file
                ///////////////////////////////////////////////////////////////
                if( localParameterFilenameParameter.size() > 1 ) {

                    throw "Error: More than one local parameter file specified";
                }
                if( localParameterFilenameParameter.size() == 0 ) {

                    throw "Error: No global local parameter file specified";
                }

                ///////////////////////////////////////////////////////////////
                // try to load the local parameter file
                ///////////////////////////////////////////////////////////////
                rc = paraObjLocalWorker.loadParameterDataFromFile(localParameterFilenameParameter.at(0));  // First argument contains parameter filename
                cout<<"Load of Local Parameter File succeeded, rc = "<<rc<<endl;

				///////////////////////////////////////////////////////////////
				//
				//
				// worflow # 3
				//
				//
				///////////////////////////////////////////////////////////////
                rc = wfObject.workflowThree(paraObjGlobal,paraObjLocalWorker,sgObject,msdlObject,spObject);
                cout<<"workflowThree rc ="<<rc<<endl;

                if( rc == 0 ) {
                    cout<<"Work flow # 3 succeeded"<<endl;
                }

            }  // end workflow 3
			

			///////////////////////////////////////////////////////////////////
			//
            // LC Candidate Aggregation Mode
            // - This is a placeholder for a future
			//
			///////////////////////////////////////////////////////////////////
			if ( modeParameter.at(0) == 4) {
                // Aggregate the LC Candidates in MZDA Level 2 files
                cout<<"LC Candidate Data Aggregation  ( workflow #4 )(-m 4)"<<endl;

			}




            ///////////////////////////////////////////////////////////////////
            //
            // This mode dynamically generates a global parameter file
            // based on a set of MZDAPARTITIONS.txt files
            //
            // This is used for automated processing of large
            // numbers of input files.
            //
            //
            ///////////////////////////////////////////////////////////////////
            if( modeParameter.at(0) == 5) {
                 cout<<"Global parameter file generation from a set of partitions files  ( workflow #5 )(-m 5)"<<endl;

                 if ( partitionsFilenameParameter.size() == 0){
                     cout<<"Error: No partitions parameter file specified\n";
                     throw "Error: No partitions parameter file specified";
                 }

                 if( outputMergedPartitionParameters.size() == 0) {
                     cout<<"Error: No output merged partitions parameter file specified\n";
                     throw "Error: No output merged partitions parameter file specified";
                 }
                 rc  = paraObjPartitions.loadParameterDataFromFile(partitionsFilenameParameter.at(0));  // First argument contains parameter filename
                 cout<<"Load of Partitions Parameter File succeeded, rc = "<<rc<<endl;
                 cout<<"partitionsFilenameParameter="<<partitionsFilenameParameter.at(0)<<endl;


                 //////////////////////////////////////////////////////////////
                 // Merged parameter datastructures
                 //////////////////////////////////////////////////////////////
                 //cluster parameters
                 vector<string> clusterinstalldirectoryDefaultMerged;
                 vector<string> clusterworkingdirectoryDefaultMerged;
                 vector<string> clustermzxmlfilenameDefaultMerged;
                 vector<int> sectioncountDefaultMerged;
                 vector<string> mzxmlfilenamebaseDefaultMerged;

                 vector<string> mzXMLFilenameMerged;

                 //////////////////////////////////////////////////////////////
                 // We need to dynamically construct the file name strings
                 // that are expected by each cluster worker, this path
                 // should be a local temporary directory such as /tmp
                 // The localTempDir parameter within the clusterConfig.txt
                 // parameter file is used to generate the full path
                 // expected
                 //////////////////////////////////////////////////////////////
                 vector<string> workermzxmlfilenamebaseDefaultMerged; // worker specific full filename


                 vector<int> mslevelselectorDefaultMerged;
                 vector<int> parglobalminscanDefaultMerged;
                 vector<int>  parglobalmaxscanDefaultMerged;
                 vector<double> parglobalminmzDefaultMerged;
                 vector<double> parglobalmaxmzDefaultMerged;
                 // MZDA Level 1 parameters
                 vector<double> parmzdeltaDefaultMerged;
                 vector<double> parmzdeltaoverlapDefaultMerged;
                 vector<double> mzbinningthresholdppmDefaultMerged;

                 // MZDA Level 2 parameters
                 vector<int> numsmoothpointDefaultMerged;
                 vector<int> minlclengthDefaultMerged;
                 vector<int> minmzlengthDefaultMerged;
                 vector<int> lcpeakapextoleranceDefaultMerged;
                 vector<double> noisethresholdlevelDefaultMerged;
                 vector<double> massresolutionDefaultMerged;
                 vector<double> mzcorrthresholdDefaultMerged;
                 vector<int> minnumprelcxicsignalsDefaultMerged;
                 vector<int> outputlogginglevelDefaultMerged;

                 // Quantification pipeline integration workflow
                 vector<string> MZDAQUANTTPPFILESDefaultMerged;
                 vector<string> MZDAQUANTDBFILESDefaultMerged;
                 vector<int> MZDAQUANTSAMPLEDefaultMerged;
                 vector<int> MZDAQUANTFRACTIONDefaultMerged;
                 vector<int> MZDAQUANTTECHNICALREPLICATEDefaultMerged;



                 // Additional parameters
                 vector<double> parglobalminRetentionTimeDefaultMerged;
                 vector<double> parglobalmaxRetentionTimeDefaultMerged;
                 vector<double> partitionScanCountMerged;


                 //////////////////////////////////////////////////////////////
                 // Open each partition file and generate an aggregate
                 // list of filenames and parameters that will
                 // be used to generate a global parameter file
                 // dynamically.
                 //////////////////////////////////////////////////////////////
                 vector<string> partitionFilenames = paraObjPartitions.getParameterValuesAsString("partitionfiles");


                 vector<string> localTempDirPath = paraObjPartitions.getParameterValuesAsString("localTempDir");
                 //////////////////////////////////////////////////////////////
                 // Sections per input file
                 //////////////////////////////////////////////////////////////
                 vector<int> sectionsPerInputDataFile;

                 // Open each partition file
                 for ( int i = 0; i < partitionFilenames.size(); i++ ){
                    cbi::ParameterLoader paraObjPartitionsItem; // Partitions file parameters
                    cout<<"i ="<<i<<",currentPartitionFilename="<<partitionFilenames.at(i)<<endl;
                    rc = paraObjPartitionsItem.loadParameterDataFromFile(partitionFilenames.at(i));



                    ///////////////////////////////////////////////////////////
                    // Since sectioncount indicates how many sections an input file was partitioned into,
                    // we only need to store the first item
                    ///////////////////////////////////////////////////////////
                    vector<string> mzxmlfilenamebaseDefault = paraObjPartitionsItem.getParameterValuesAsString("mzxmlfilenamebase");
                    mzxmlfilenamebaseDefaultMerged.push_back(mzxmlfilenamebaseDefault.at(0));

                    ///////////////////////////////////////////////////////////
                    // Create the path to the input mzXML files in the
                    // local node's local disk temp space.
                    ///////////////////////////////////////////////////////////
                    vector<string> mzXMLFilename = paraObjPartitionsItem.getParameterValuesAsString("mzxmlfilenamebase");
                    for ( int j = 0; j < mzXMLFilename.size(); j++ ){
                        workermzxmlfilenamebaseDefaultMerged.push_back(localTempDirPath.at(0)+"/"+mzXMLFilename.at(j));
                    }


                    ///////////////////////////////////////////////////////////
                    // Since sectioncount indicates how many sections an input file was partitioned into,
                    // we only need to store the first item
                    ///////////////////////////////////////////////////////////
                    vector<int> sectioncountDefault = paraObjPartitionsItem.getParameterValuesAsInt("sectioncount");
                    sectioncountDefaultMerged.push_back(sectioncountDefault.at(0));







                    //cout<<"  mslevelselectorDefault=";
                    vector<int> mslevelselectorDefault = paraObjPartitionsItem.getParameterValuesAsInt("mslevelselector");
                    for ( int j = 0; j < mslevelselectorDefault.size(); j++ ){
                        //cout<<mslevelselectorDefault.at(j)<<" ";
                        mslevelselectorDefaultMerged.push_back(mslevelselectorDefault.at(j));
                    }
                    //cout<<endl;

                    //cout<<"  parglobalminscanDefault=";
                    vector<int> parglobalminscanDefault = paraObjPartitionsItem.getParameterValuesAsInt("parglobalminscan");
                    for ( int j = 0; j < parglobalminscanDefault.size(); j++ ){
                        //cout<<parglobalminscanDefault.at(j)<<" ";
                        parglobalminscanDefaultMerged.push_back(parglobalminscanDefault.at(j));
                    }
                    //cout<<endl;

                    //cout<<"  parglobalmaxscanDefault=";
                    vector<int> parglobalmaxscanDefault = paraObjPartitionsItem.getParameterValuesAsInt("parglobalmaxscan");
                    for ( int j = 0; j < parglobalmaxscanDefault.size(); j++ ){
                        //cout<<parglobalmaxscanDefault.at(j)<<" ";
                        parglobalmaxscanDefaultMerged.push_back(parglobalmaxscanDefault.at(j));
                    }
                    //cout<<endl;

                    //cout<<"  parglobalminmzDefault=";
                    vector<double> parglobalminmzDefault = paraObjPartitionsItem.getParameterValuesAsDouble("parglobalminmz");
                    for ( int j = 0; j < parglobalminmzDefault.size(); j++ ){
                        //cout<<parglobalminmzDefault.at(j)<<" ";
                        parglobalminmzDefaultMerged.push_back(parglobalminmzDefault.at(j));
                    }
                    //cout<<endl;

                    //cout<<"  parglobalmaxmzDefault=";
                    vector<double> parglobalmaxmzDefault = paraObjPartitionsItem.getParameterValuesAsDouble("parglobalmaxmz");
                    for ( int j = 0; j < parglobalmaxmzDefault.size(); j++ ){
                        //cout<<parglobalmaxmzDefault.at(j)<<" ";
                        parglobalmaxmzDefaultMerged.push_back(parglobalmaxmzDefault.at(j));
                    }
                    //cout<<endl;

                    //cout<<"  parmzdeltaDefault=";
                    vector<double> parmzdeltaDefault = paraObjPartitionsItem.getParameterValuesAsDouble("parmzdelta");
                    for ( int j = 0; j < parmzdeltaDefault.size(); j++ ){
                        //cout<<parmzdeltaDefault.at(j)<<" ";
                        parmzdeltaDefaultMerged.push_back(parmzdeltaDefault.at(j));
                    }
                    //cout<<endl;

                    //cout<<"  parmzdeltaoverlapDefault=";
                    vector<double> parmzdeltaoverlapDefault = paraObjPartitionsItem.getParameterValuesAsDouble("parmzdeltaoverlap");
                    for ( int j = 0; j < parmzdeltaoverlapDefault.size(); j++ ){
                        //cout<<parmzdeltaoverlapDefault.at(j)<<" ";
                        parmzdeltaoverlapDefaultMerged.push_back(parmzdeltaoverlapDefault.at(j));
                    }
                    //cout<<endl;

                    //cout<<"  mzbinningthresholdppmDefault=";
                    vector<double> mzbinningthresholdppmDefault = paraObjPartitionsItem.getParameterValuesAsDouble("mzbinningthresholdppm");
                    for ( int j = 0; j < mzbinningthresholdppmDefault.size(); j++ ){
                        //cout<<mzbinningthresholdppmDefault.at(j)<<" ";
                        mzbinningthresholdppmDefaultMerged.push_back(mzbinningthresholdppmDefault.at(j));
                    }
                    //cout<<endl;

                    ///////////////////////////////////////////////////////////
                    // Adding the MZDA Level 2 parameters
                    ///////////////////////////////////////////////////////////
                    vector<int> numsmoothpointDefault = paraObjPartitionsItem.getParameterValuesAsInt("numsmoothpoint");
                    for ( int j = 0; j < numsmoothpointDefault.size(); j++ ){
                        numsmoothpointDefaultMerged.push_back(numsmoothpointDefault.at(j));
                    }

                    vector<int> minlclengthDefault = paraObjPartitionsItem.getParameterValuesAsInt("minlclength");
                    for ( int j = 0; j < minlclengthDefault.size(); j++ ){

                        minlclengthDefaultMerged.push_back(minlclengthDefault.at(j));
                    }

                    vector<int> minmzlengthDefault = paraObjPartitionsItem.getParameterValuesAsInt("minmzlength");
                    for ( int j = 0; j < minmzlengthDefault.size(); j++ ){
                        minmzlengthDefaultMerged.push_back(minmzlengthDefault.at(j));
                    }

                    vector<int> lcpeakapextoleranceDefault = paraObjPartitionsItem.getParameterValuesAsInt("lcpeakapextolerance");
                    for ( int j = 0; j < lcpeakapextoleranceDefault.size(); j++ ){
                        lcpeakapextoleranceDefaultMerged.push_back(lcpeakapextoleranceDefault.at(j));
                    }

                    vector<double> noisethresholdlevelDefault = paraObjPartitionsItem.getParameterValuesAsDouble("noisethresholdlevel");
                    for ( int j = 0; j < noisethresholdlevelDefault.size(); j++ ){
                        noisethresholdlevelDefaultMerged.push_back(noisethresholdlevelDefault.at(j));
                    }

                    vector<double> massresolutionDefault = paraObjPartitionsItem.getParameterValuesAsDouble("massresolution");
                    for ( int j = 0; j < massresolutionDefault.size(); j++ ){
                        massresolutionDefaultMerged.push_back(massresolutionDefault.at(j));
                    }

                    vector<double> mzcorrthresholdDefault = paraObjPartitionsItem.getParameterValuesAsDouble("mzcorrthreshold");
                    for ( int j = 0; j < mzcorrthresholdDefault.size(); j++ ){
                        mzcorrthresholdDefaultMerged.push_back(mzcorrthresholdDefault.at(j));
                    }

                    vector<int> minnumprelcxicsignalsDefault = paraObjPartitionsItem.getParameterValuesAsInt("minnumprelcxicsignals");
                    for ( int j = 0; j < minnumprelcxicsignalsDefault.size(); j++ ){
                        minnumprelcxicsignalsDefaultMerged.push_back(minnumprelcxicsignalsDefault.at(j));
                    }

                    vector<int> outputlogginglevelDefault = paraObjPartitionsItem.getParameterValuesAsInt("outputlogginglevel");
                    for ( int j = 0; j < outputlogginglevelDefault.size(); j++ ){
                        outputlogginglevelDefaultMerged.push_back(outputlogginglevelDefault.at(j));
                    }

                    ////////////////////////////////////////////////////////////
                    // START: Quantification Pipeline Automation Worklow Parameters
                    ////////////////////////////////////////////////////////////
                    vector<string> MZDAQUANTTPPFILESDefault = paraObjPartitionsItem.getParameterValuesAsString("MZDAQUANTTPPFILES");
                    MZDAQUANTTPPFILESDefaultMerged.push_back(MZDAQUANTTPPFILESDefault.at(0));  // Only use the first item in list

                    vector<string> MZDAQUANTDBFILESDefault = paraObjPartitionsItem.getParameterValuesAsString("MZDAQUANTDBFILES");
                    MZDAQUANTDBFILESDefaultMerged.push_back(MZDAQUANTDBFILESDefault.at(0));   // Only use the first item in list

                    vector<int> MZDAQUANTSAMPLEDefault = paraObjPartitionsItem.getParameterValuesAsInt("MZDAQUANTSAMPLE");
                    MZDAQUANTSAMPLEDefaultMerged.push_back(MZDAQUANTSAMPLEDefault.at(0));     // Only use the first item in list

                    vector<int> MZDAQUANTFRACTIONDefault = paraObjPartitionsItem.getParameterValuesAsInt("MZDAQUANTFRACTION");
                    MZDAQUANTFRACTIONDefaultMerged.push_back(MZDAQUANTFRACTIONDefault.at(0)); // Only use the first item in list

                    vector<int> MZDAQUANTTECHNICALREPLICATEDefault = paraObjPartitionsItem.getParameterValuesAsInt("MZDAQUANTTECHNICALREPLICATE");
                    MZDAQUANTTECHNICALREPLICATEDefaultMerged.push_back(MZDAQUANTTECHNICALREPLICATEDefault.at(0)); // Only use the first item in list

                    ////////////////////////////////////////////////////////////
                    // END: Quantification Pipeline Automation Worklow Parameters
                    ////////////////////////////////////////////////////////////

                    ///////////////////////////////////////////////////////////
                    // Additional parameters
                    ///////////////////////////////////////////////////////////
                    //cout<<"  parglobalminRetentionTimeDefault=";
                    vector<double> parglobalminRetentionTimeDefault = paraObjPartitionsItem.getParameterValuesAsDouble("parglobalminretentiontime");
                    for ( int j = 0; j < parglobalminRetentionTimeDefault.size(); j++ ){
                        //cout<<parglobalminRetentionTimeDefault.at(j)<<" ";
                        parglobalminRetentionTimeDefaultMerged.push_back(parglobalminRetentionTimeDefault.at(j));
                    }
                    //cout<<endl;

                    //cout<<"  parglobalmaxRetentionTimeDefault=";
                    vector<double> parglobalmaxRetentionTimeDefault = paraObjPartitionsItem.getParameterValuesAsDouble("parglobalmaxretentiontime");
                    for ( int j = 0; j < parglobalmaxRetentionTimeDefault.size(); j++ ){
                        //cout<<parglobalmaxRetentionTimeDefault.at(j)<<" ";
                        parglobalmaxRetentionTimeDefaultMerged.push_back(parglobalmaxRetentionTimeDefault.at(j));
                    }
                    //cout<<endl;

                    //cout<<"  partitionScanCount=";
                    vector<int> partitionScanCount = paraObjPartitionsItem.getParameterValuesAsInt("partitionscancount");
                    for ( int j = 0; j < partitionScanCount.size(); j++ ){
                        //cout<<partitionScanCount.at(j)<<" ";
                        partitionScanCountMerged.push_back(partitionScanCount.at(j));
                    }
                    //cout<<endl;

                    ///////////////////////////////////////////////////////////
                    // Set the number of sections in this file
                    ///////////////////////////////////////////////////////////
                    sectionsPerInputDataFile.push_back(mzXMLFilename.size());



                 }  // end loop through all partition files (i)

                 //////////////////////////////////////////////////////////////
                 // Workflow:
                 // Call the partition file processing workflow
                 //
                 //////////////////////////////////////////////////////////////
                 int totalInputFiles = partitionFilenames.size(); // total number of mzXML input files
                 int totalSections = mslevelselectorDefaultMerged.size();  // total number of sections


                 cout<<"totalInputFiles="<<totalInputFiles<<endl;
                 cout<<"totalSections="<<totalSections<<endl;
                 cout<<"mzxmlfile,mslevel,minscan,maxscan,minmz,maxmz,mzdelta,mzoverlap,mzbinningthreshold,minrt,maxrt,scancount"<<endl;
                 for ( int i = 0; i < mzXMLFilenameMerged.size(); i++ ){
                     cout<<"section="<<i;

                     cout<<","<<clusterinstalldirectoryDefaultMerged.at(i);
                     cout<<","<<clusterworkingdirectoryDefaultMerged.at(i);
                     cout<<","<<clustermzxmlfilenameDefaultMerged.at(i);
                     cout<<","<<mzXMLFilenameMerged.at(i);
                     cout<<","<<mslevelselectorDefaultMerged.at(i);
                     cout<<","<<parglobalminscanDefaultMerged.at(i);
                     cout<<","<<parglobalmaxscanDefaultMerged.at(i);
                     cout<<","<<parglobalminmzDefaultMerged.at(i);
                     cout<<","<<parglobalmaxmzDefaultMerged.at(i);
                     cout<<","<<parmzdeltaDefaultMerged.at(i);
                     cout<<","<<parmzdeltaoverlapDefaultMerged.at(i);
                     cout<<","<<mzbinningthresholdppmDefaultMerged.at(i);

                     // MZDA Level 2 parameters
                     cout<<","<<numsmoothpointDefaultMerged.at(i);
                     cout<<","<<minlclengthDefaultMerged.at(i);
                     cout<<","<<minmzlengthDefaultMerged.at(i);
                     cout<<","<<lcpeakapextoleranceDefaultMerged.at(i);
                     cout<<","<<noisethresholdlevelDefaultMerged.at(i);
                     cout<<","<<massresolutionDefaultMerged.at(i);
                     cout<<","<<mzcorrthresholdDefaultMerged.at(i);
                     cout<<","<<minnumprelcxicsignalsDefaultMerged.at(i);
                     cout<<","<<outputlogginglevelDefaultMerged.at(i);




                     cout<<","<<parglobalminRetentionTimeDefaultMerged.at(i);
                     cout<<","<<parglobalmaxRetentionTimeDefaultMerged.at(i);
                     cout<<","<<partitionScanCountMerged.at(i);
                     cout<<endl;
                 }


                 //////////////////////////////////////////////////////////////
                 // The Quantification Pipeline integration parameter structures
                 // should be the same length as the number of input
                 // mzXML files, as there is 1 partition file per input mzXML
                 // file. @nramirez 10/27/14
                 //////////////////////////////////////////////////////////////
                 assert( MZDAQUANTTPPFILESDefaultMerged.size() == totalInputFiles );
                 assert( MZDAQUANTDBFILESDefaultMerged.size() == totalInputFiles );
                 assert( MZDAQUANTSAMPLEDefaultMerged.size() == totalInputFiles );
                 assert( MZDAQUANTFRACTIONDefaultMerged.size() == totalInputFiles );
                 assert( MZDAQUANTTECHNICALREPLICATEDefaultMerged.size() == totalInputFiles );

                 //////////////////////////////////////////////////////////////
                 // Generate the output merged parameter file
                 //////////////////////////////////////////////////////////////
                 std::ofstream outputFileStream;
                 outputFileStream.open(outputMergedPartitionParameters.at(0).c_str(), std::ios::out);
                 try{
                     if( outputFileStream.is_open()) {
                         outputFileStream<<"totalInputFiles="<<totalInputFiles<<endl;
                         outputFileStream<<"totalSections="<<totalSections<<endl;

                         outputFileStream<<"sectioncount=";
                         for ( int i = 0; i < sectioncountDefaultMerged.size(); i++ ){
                             outputFileStream<<sectioncountDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;


                         outputFileStream<<"mzxmlfilenamebase=";
                         for ( int i = 0; i < mzxmlfilenamebaseDefaultMerged.size(); i++ ){
                             outputFileStream<<mzxmlfilenamebaseDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;



                         //////////////////////////////////////////////////////
                         // Cluster ( a.k.a.) mzXML filenames
                         //   - For example, this should be the paths
                         //     that the compute nodes will be searching
                         //     when looking for their respective mzMXL
                         //     files.
                         //   - This paths should normally be the node's
                         //     local /tmp directory.
                         //
                         //
                         //////////////////////////////////////////////////////  
                         // we'll want to remove the workermzxmlfilename field
                         // and only have a single field for input filenames.
                         outputFileStream<<"workermzxmlfilename=";
                         for ( int i = 0; i < workermzxmlfilenamebaseDefaultMerged.size(); i++ ){
                             outputFileStream<<workermzxmlfilenamebaseDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;

                         // added on 07/24/13 to enable cluster workflow
                         outputFileStream<<"mzxmlfilename=";
                         for ( int i = 0; i < workermzxmlfilenamebaseDefaultMerged.size(); i++ ){
                             outputFileStream<<workermzxmlfilenamebaseDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;

                         outputFileStream<<"mslevelselector=";
                         for ( int i = 0; i < mslevelselectorDefaultMerged.size(); i++ ){
                             outputFileStream<<mslevelselectorDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;

                         outputFileStream<<"parglobalminscan=";
                         for ( int i = 0; i < parglobalminscanDefaultMerged.size(); i++ ){
                             outputFileStream<<parglobalminscanDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;

                         outputFileStream<<"parglobalmaxscan=";
                         for ( int i = 0; i < parglobalmaxscanDefaultMerged.size(); i++ ){
                             outputFileStream<<parglobalmaxscanDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;


                         outputFileStream<<"parglobalminmz=";
                         for ( int i = 0; i < parglobalminmzDefaultMerged.size(); i++ ){
                             outputFileStream<<parglobalminmzDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;


                         outputFileStream<<"parglobalmaxmz=";
                         for ( int i = 0; i < parglobalmaxmzDefaultMerged.size(); i++ ){
                             outputFileStream<<parglobalmaxmzDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;


                         outputFileStream<<"parmzdelta=";
                         for ( int i = 0; i < parmzdeltaDefaultMerged.size(); i++ ){
                             outputFileStream<<parmzdeltaDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;


                         outputFileStream<<"parmzdeltaoverlap=";
                         for ( int i = 0; i < parmzdeltaoverlapDefaultMerged.size(); i++ ){
                             outputFileStream<<parmzdeltaoverlapDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;


                         outputFileStream<<"mzbinningthresholdppm=";
                         for ( int i = 0; i < mzbinningthresholdppmDefaultMerged.size(); i++ ){
                             outputFileStream<<mzbinningthresholdppmDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;


                         // MZDA Level 2 parameters
                         outputFileStream<<"numsmoothpoint=";
                         for ( int i = 0; i < numsmoothpointDefaultMerged.size(); i++ ){
                             outputFileStream<<numsmoothpointDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;
                         outputFileStream<<"minlclength=";
                         for ( int i = 0; i < minlclengthDefaultMerged.size(); i++ ){
                             outputFileStream<<minlclengthDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;
                         outputFileStream<<"minmzlength=";
                         for ( int i = 0; i < minmzlengthDefaultMerged.size(); i++ ){
                             outputFileStream<<minmzlengthDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;
                         outputFileStream<<"noisethresholdlevel=";
                         for ( int i = 0; i < noisethresholdlevelDefaultMerged.size(); i++ ){
                             outputFileStream<<noisethresholdlevelDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;
                         outputFileStream<<"massresolution=";
                         for ( int i = 0; i < massresolutionDefaultMerged.size(); i++ ){
                             outputFileStream<<massresolutionDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;

                         outputFileStream<<"lcpeakapextolerance=";
                         for ( int i = 0; i < lcpeakapextoleranceDefaultMerged.size(); i++ ){
                             outputFileStream<<lcpeakapextoleranceDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;

                         outputFileStream<<"mzcorrthreshold=";
                         for ( int i = 0; i < mzcorrthresholdDefaultMerged.size(); i++ ){
                             outputFileStream<<mzcorrthresholdDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;

                         outputFileStream<<"minnumprelcxicsignals=";
                         for ( int i = 0; i < minnumprelcxicsignalsDefaultMerged.size(); i++ ){
                             outputFileStream<<minnumprelcxicsignalsDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;

                         outputFileStream<<"outputlogginglevel=";
                         for ( int i = 0; i < outputlogginglevelDefaultMerged.size(); i++ ){
                             outputFileStream<<outputlogginglevelDefaultMerged.at(i)<<" ";
                         }
                         outputFileStream<<std::endl;

                         //////////////////////////////////////////////////////
                         // START: QUANTIFICATION PIPELINE INTEGRATION
                         // We only want to use the first item in the
                         // list from the partition files for the
                         // Quantification integration, as the Quantification
                         // modules
                         //////////////////////////////////////////////////////
                         outputFileStream<<"MZDAQUANTTPPFILES=";
                         for ( int i = 0; i < MZDAQUANTTPPFILESDefaultMerged.size()-1; i++ ){
                             outputFileStream<<MZDAQUANTTPPFILESDefaultMerged.at(i)<<",";  // Quant Module uses comma separator
                         }
                         outputFileStream<<MZDAQUANTTPPFILESDefaultMerged.at(MZDAQUANTTPPFILESDefaultMerged.size()-1)<<" ";  // Remove last comma
                         outputFileStream<<std::endl;

                         outputFileStream<<"MZDAQUANTDBFILES=";
                         for ( int i = 0; i < MZDAQUANTDBFILESDefaultMerged.size()-1; i++ ){
                             outputFileStream<<MZDAQUANTDBFILESDefaultMerged.at(i)<<",";  // Quant Module uses comma separator
                         }
                         outputFileStream<<MZDAQUANTDBFILESDefaultMerged.at(MZDAQUANTDBFILESDefaultMerged.size()-1)<<" ";  // Remove last comma
                         outputFileStream<<std::endl;

                         outputFileStream<<"MZDAQUANTSAMPLE=";
                         for ( int i = 0; i < MZDAQUANTSAMPLEDefaultMerged.size()-1; i++ ){
                             outputFileStream<<MZDAQUANTSAMPLEDefaultMerged.at(i)<<",";  // Quant Module uses comma separator
                         }
                         outputFileStream<<MZDAQUANTSAMPLEDefaultMerged.at(MZDAQUANTSAMPLEDefaultMerged.size()-1)<<" ";  // Remove last comma
                         outputFileStream<<std::endl;

                         outputFileStream<<"MZDAQUANTFRACTION=";
                         for ( int i = 0; i < MZDAQUANTFRACTIONDefaultMerged.size()-1; i++ ){
                             outputFileStream<<MZDAQUANTFRACTIONDefaultMerged.at(i)<<",";  // Quant Module uses comma separator
                         }
                         outputFileStream<<MZDAQUANTFRACTIONDefaultMerged.at(MZDAQUANTFRACTIONDefaultMerged.size()-1)<<" ";  // Remove last comma
                         outputFileStream<<std::endl;


                         outputFileStream<<"MZDAQUANTTECHNICALREPLICATE=";
                         for ( int i = 0; i < MZDAQUANTTECHNICALREPLICATEDefaultMerged.size()-1; i++ ){
                             outputFileStream<<MZDAQUANTTECHNICALREPLICATEDefaultMerged.at(i)<<",";  // Quant Module uses comma separator
                         }
                         outputFileStream<<MZDAQUANTTECHNICALREPLICATEDefaultMerged.at(MZDAQUANTTECHNICALREPLICATEDefaultMerged.size()-1)<<" ";  // Remove last comma
                         outputFileStream<<std::endl;

                         //////////////////////////////////////////////////////
                         // START: QUANTIFICATION PIPELINE INTEGRATION
                         //////////////////////////////////////////////////////


                     }
                     else{
                         // output file is not open, no need to try to close it
                         return -2;
                     }
                 }
                 catch (std::exception& e) {
                     std::cout << e.what() <<std::endl;
                     std::cout <<"mzda: exception while trying to open output merged parameter file"<<std::endl;
                     // make sure to close the file stream
                     outputFileStream.close();
                 }  // end trycatch block

                 // make sure to close the output file
                 outputFileStream.close();
                 //////////////////////////////////////////////////////////////
                 //////////////////////////////////////////////////////////////


                 //////////////////////////////////////////////////////////////
                 //////////////////////////////////////////////////////////////
                 //// GENERATE A SET OF SHELL SCRIPTS AND
                 //// THEIR CORRESPONDING JOB SCRIPTS
                 //// FOR USE ON A SHARED HPC CLUSTER.
                 ////
                 //// One job script per input mzXML file will be
                 //// generated.
                 ////
                 //// An mzXML file is partitioned into a set of
                 //// Sections in the retention time dimension.
                 //// Each section is further partitioned into a set
                 //// of ScanGroups.
                 //// Each ScanGroup is further partitioned into a
                 //// set of PLCROI's ( Pre-LC Regions of Interest.
                 ////
                 //// The LC Candidate identification algorithm is performed
                 //// on these PLCROI's.
                 ////
                 //// For cluster enablement, the processing for
                 //// each mzXML file is managed by 1 job script.
                 //// Each of these manages a set of lower level
                 //// scripts that manage the MZDA Level1,2,3 processing
                 //// stages while automatically cleaning up
                 //// the working directory before and after each
                 //// processing stage to make the processing
                 //// scalable across any number of files, using a
                 //// controlled amount of disk, cpu, ram on each
                 //// compute core.
                 //// Without this control, a processing run would
                 //// eventually cause < any > cluster to run out
                 //// of resource eventually as the number of
                 //// input files exceeds a certain number.
                 //// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 //// Our approach scales to ANY number of input
                 //// files using the SAME amount of user configurable
                 //// compute/disk/ram resources.
                 //// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 ////
                 //////////////////////////////////////////////////////////////
                 //////////////////////////////////////////////////////////////
                 vector<string> sharedInputDataDir = paraObjPartitions.getParameterValuesAsString("sharedInputDataDir");
                 vector<string> sharedOutputDataDir = paraObjPartitions.getParameterValuesAsString("sharedOutputDataDir");
                 vector<string> sharedInstallDir = paraObjPartitions.getParameterValuesAsString("sharedInstallDir");
                 vector<string> localTempDir = paraObjPartitions.getParameterValuesAsString("localTempDir");

                 //////////////////////////////////////////////////////////////
                 // Quantification Module Pipeline Integration
                 //////////////////////////////////////////////////////////////
                 vector<string> mzdaQuantInstallDir = paraObjPartitions.getParameterValuesAsString("mzdaQuantInstallDir");
                 vector<int> MZDAQUANTMAXFRACTIONS = paraObjPartitions.getParameterValuesAsInt("MZDAQUANTMAXFRACTIONS");

                 //////////////////////////////////////////////////////////////
                 // Below are the job script specific variables for a general HPC Cluster,
                 // These values are contained within the clusterConfig.txt
                 // file which is dynamically generated by the
                 // mzdaController.sh script contained within the product
                 // directory.
                 //////////////////////////////////////////////////////////////
                 vector<string> mcrPath = paraObjPartitions.getParameterValuesAsString("mcrPath");
                 vector<string> nodesPerJob = paraObjPartitions.getParameterValuesAsString("nodesPerJob");
                 vector<string> coresPerJob = paraObjPartitions.getParameterValuesAsString("coresPerJob");
                 vector<string> queueName = paraObjPartitions.getParameterValuesAsString("queueName");
                 vector<string> maxRuntimePerJob = paraObjPartitions.getParameterValuesAsString("maxRuntimePerJob");
                 vector<string> submitCommand = paraObjPartitions.getParameterValuesAsString("submitCommand");
                 vector<string> submitDelayTime= paraObjPartitions.getParameterValuesAsString("submitDelayTime");
                 vector<string> parallelEnvironment = paraObjPartitions.getParameterValuesAsString("parallelEnvironment");
                 vector<string> queueNamePartTwo = paraObjPartitions.getParameterValuesAsString("queueNamePartTwo");
                 vector<string> gridManager = paraObjPartitions.getParameterValuesAsString("gridManager");

                 //////////////////////////////////////////////////////////////
                 // TPPXTandem Input Parameters
                 //////////////////////////////////////////////////////////////
                 vector<string> tppXTandemParameterFile = paraObjPartitions.getParameterValuesAsString("tppXTandemParameterFile");
                 vector<string> tppXTandemDatabaseFastaFile = paraObjPartitions.getParameterValuesAsString("tppXTandemDatabaseFastaFile");
                 vector<string> tppXTandemBinDir = paraObjPartitions.getParameterValuesAsString("tppXTandemBinDir");
                 vector<string> tppXTandemLocalTempDir = paraObjPartitions.getParameterValuesAsString("tppXTandemLocalTempDir"); // We need to be able to control output director of tpp xtandem

                 //////////////////////////////////////////////////////////////
                 // clusterConfig.txt parameter validation
                 //////////////////////////////////////////////////////////////
                 assert(sharedInputDataDir.size()>0);
                 assert(sharedOutputDataDir.size()>0);
                 assert(sharedInstallDir.size()>0);
                 assert(localTempDir.size()>0);
                 assert(mzdaQuantInstallDir.size()>0);
                 assert(mcrPath.size()>0);
                 assert(nodesPerJob.size()>0);
                 assert(coresPerJob.size()>0);
                 assert(queueName.size()>0);
                 assert(maxRuntimePerJob.size()>0);
                 assert(submitCommand.size()>0);
                 assert(submitDelayTime.size()>0);
                 assert(parallelEnvironment.size()>0);
                 assert(queueNamePartTwo.size()>0);
                 assert(gridManager.size()>0);


                 assert(tppXTandemParameterFile.size() >0);
                 assert(tppXTandemDatabaseFastaFile.size() >0);
                 assert(tppXTandemBinDir.size()>0);
                 assert(tppXTandemLocalTempDir.size()>0);


                 assert(MZDAQUANTMAXFRACTIONS.size()>0);   // make sure this parameter is specified
                 assert(MZDAQUANTMAXFRACTIONS.size()==1);  // make sure only 1 number specified
                 assert(MZDAQUANTMAXFRACTIONS.at(0)<500);  // a precautionary limit of max number fractions

                 //////////////////////////////////////////////////////////////
                 //////////////////////////////////////////////////////////////
                 //////////////////////////////////////////////////////////////

                 int currentSectionStart = 1;
                 int currentSectionEnd = 1;
                 for ( int i = 0; i < totalInputFiles; i++ ) {


                 std::ofstream outputFileStream;
                 stringstream ss;
                 ss<<(i+1); // integer to string conversion
                 string fileNumber = ss.str(); // converting stringstream to string object

                 string outputClusterJobFilename= outputClusterJobsFilenameBase.at(0) + fileNumber + ".sh";

                 outputFileStream.open(outputClusterJobFilename.c_str(), std::ios::out);
                 try{
                     if( outputFileStream.is_open()) {

                         // The start and end points
                         // update the files section end point
                         currentSectionEnd = currentSectionStart + sectioncountDefaultMerged.at(i) -1;

                         //////////////////////////////////////////////////////
                         //////////////////////////////////////////////////////
                         outputFileStream<<"#!/bin/bash"<<endl;
                         outputFileStream<<"JOBNUMBER="<<fileNumber<<endl;
                         outputFileStream<<"SECTIONNUMBERSTART="<<"\""<<currentSectionStart<<"\""<<endl;
                         outputFileStream<<"SECTIONNUMBEREND="<<"\""<<currentSectionEnd<<"\""<<endl;
                         outputFileStream<<"FILENAME="<<mzxmlfilenamebaseDefaultMerged.at(i)<<endl;
                         outputFileStream<<"INSTALLDIR="<<sharedInstallDir.at(0)<<endl; // shared location where mzdasoft is installed
                         outputFileStream<<"WORKDIR="<<localTempDir.at(0)<<endl; // node local disk ( usually /tmp/... )
                         outputFileStream<<"MCRPATH="<<mcrPath.at(0)<<endl; // path to MCR installed HPC Center
                         //////////////////////////////////////////////////////
                         // Current version requires at MOST 1 job per system at any given time.
                         // more elaborate schemes can be easily implemented within this dynamically
                         // generated script.
                         //////////////////////////////////////////////////////
                         //{ $WORK/mzdataccwork/mzdaWorkflow2_1.sh && sleep 10 && echo "mzdasoft level1 complete" && sleep 10 && $WORK/mzdataccwork/mzdaWorkflow3_1.sh && sleep 10 && echo "mzdasoft level 2 complete" && sleep 10 && $WORK/mzdataccwork/mzdaWorkflow4_1.sh && sleep 10 && echo "mzdasoft level 3 complete" && wc $WORK/mzdataccwork/mzdaLevelThree1.dat; } &&
                         //{ $WORK/mzdataccwork/mzdaWorkflow2_2.sh && sleep 10 && echo "mzdasoft level1 complete" && sleep 10 && $WORK/mzdataccwork/mzdaWorkflow3_2.sh && sleep 10 && echo "mzdasoft level 2 complete" && sleep 10 && $WORK/mzdataccwork/mzdaWorkflow4_2.sh && sleep 10 && echo "mzdasoft level 3 complete" && wc $WORK/mzdataccwork/mzdaLevelThree2.dat; } &&
                         //cat < /work/01759/nramirez/mzdataccwork/mzdaLevelThree1.dat >> /work/01759/nramirez/mzdataccwork/mzdaLevelThreeFile1.dat &&
                         //cat < /work/01759/nramirez/mzdataccwork/mzdaLevelThree2.dat >> /work/01759/nramirez/mzdataccwork/mzdaLevelThreeFile1.dat &&
                         //echo "File 1 Processing Complete" &&
                         //wc /work/01759/nramirez/mzdataccwork/mzdaLevelThreeFile1.dat &&
                         //rm $WORK/mzdataccwork/mzdaLevelThree1.dat &&
                         //rm $WORK/mzdataccwork/mzdaLevelThree2.dat


                         // Run MZDA Level 1 and Level 2 processing
                         for ( int k = currentSectionStart; k <= currentSectionEnd; k++ ){
                             outputFileStream<<"{ $WORKDIR/mzdaWorkflow2_"<<k<<".sh"<<" && sleep 10 && echo \"mzdasoft level1 complete\" && sleep 10 && $WORKDIR/mzdaWorkflow3_"<<k<<".sh"<<" && sleep 10 && echo \"mzdasoft level 2 complete\" && sleep 10 && $WORKDIR/mzdaWorkflow4_"<<k<<".sh"<<" && sleep 10 && echo \"mzdasoft level 3 complete\" && wc $WORKDIR/mzdaLevelThree"<<k<<".dat; } &&"<<endl;
                         }

                         // Concatenate results info the mzdaLevelThreeFile#.dat file

                         //////////////////////////////////////////////////////
                         // Per 08/16/13 meeting with Michelle Zhang,
                         // We need to change the filename from using
                         // the file number in the list to also use
                         // the source mzXML filename.
                         //
                         //////////////////////////////////////////////////////
                         string outputFilenameNameExtension = mzxmlfilenamebaseDefaultMerged.at(i);
                         assert( outputFilenameNameExtension.size() > 0);

                         for ( int k = currentSectionStart; k <= currentSectionEnd; k++ ){
                             outputFileStream<<"cat < $WORKDIR/mzdaLevelThree"<<k<<".dat >> $WORKDIR/mzdaLevelThreeFile"<<fileNumber<<"-"<<outputFilenameNameExtension<<".dat &&"<<endl;
                         }

                         // Provide a completion indicator
                         outputFileStream<<"echo \"File "<<fileNumber<<" Processing Complete\" &&"<<endl;

                         // Count the total number of LCC's found within the file being processed.
                         outputFileStream<<"wc $WORKDIR/mzdaLevelThreeFile"<<fileNumber<<"-"<<outputFilenameNameExtension<<".dat &&"<<endl;


                         for ( int k = currentSectionStart; k <= currentSectionEnd; k++ ){
                             outputFileStream<<"rm $WORKDIR/mzdaLevelThree"<<k<<".dat"<<endl;
                         }

                         //////////////////////////////////////////////////////
                         // Post-processing, MZDA Level 3(.dat) to Matlab(.mat) format conversion
                         //////////////////////////////////////////////////////
                         outputFileStream<<"cd $INSTALLDIR/interfaces/matlabconversiontool"<<endl;
                         outputFileStream<<"./run.sh $MCRPATH $WORKDIR/mzdaLevelThreeFile"<<fileNumber<<"-"<<outputFilenameNameExtension<<".dat "<<"$WORKDIR/mzdaLevelThreeFile"<<fileNumber<<"-"<<outputFilenameNameExtension<<".mat "<<endl;



                         //////////////////////////////////////////////////////
                         // DATABASE CONVERSION: Convert MZDA Level 3 output in .dat
                         //  format to SQLite database format(.db) with special indices
                         //  for high speed range queries on m/z and retention time
                         //  ranges.
                         //
                         // The MZDA Level 3 to SQLITE database format converter.
                         // [nelson@puma matlabconversiontool]$ time ./mzdalevelthreequery ~/mzdaLevelThreeFile8-Yeast-SILAC-1to1repeat-1.mzXML.dat ~/Desktop/querytest111213.db
                         //
                         //
                         ///////////////////////////////////////////////////////
                         outputFileStream<<"cd $INSTALLDIR/interfaces/matlabconversiontool"<<endl;
                         outputFileStream<<"./mzdalevelthreequery $WORKDIR/mzdaLevelThreeFile"<<fileNumber<<"-"<<outputFilenameNameExtension<<".dat "<<"$WORKDIR/mzdaLevelThreeFile"<<fileNumber<<"-"<<outputFilenameNameExtension<<".db "<<endl;


                         //////////////////////////////////////////////////////
                         // DATABASE VALIDATION PROCESS
                         /////////////////////////////////////////////////////
                         outputFileStream<<"echo \"Allowing 60 seconds for post .dat to .db query process\""<<endl;
                         outputFileStream<<"sleep 60"<<endl;
                         outputFileStream<<"cd $INSTALLDIR/product"<<endl;
                         string outputFilenameNameExtensionNoFileType = outputFilenameNameExtension;
                         size_t location = outputFilenameNameExtensionNoFileType.find(".mzXML");
                         outputFilenameNameExtensionNoFileType.replace(location,std::string(".mzXML").length(),"");
                         outputFileStream<<"./validateDatabaseFile.sh $INSTALLDIR $WORKDIR "<<"mzdaLevelThreeFile"<<fileNumber<<"-"<<outputFilenameNameExtensionNoFileType<<endl;
                         outputFileStream<<"sleep 60"<<endl;
                         outputFileStream<<"echo \"MZDASoft Parallel Peak Extractor Complete\""<<endl;
                         //////////////////////////////////////////////////////
                         // MZDAQuantify Module:
                         // - This is a Matlab Compiled executable module
                         //   that will run the Quantification workflow.
                         // - The Quantification module uses the Matlab
                         //   interface to SQLITE to access the SQLITE
                         //   database file.
                         // - Quantification is performed on a per mzXML file
                         //   basis, and the results are subsequently
                         //   merged and or analyzed by post processing
                         //   data mining tools.
                         // - Once the Quantification module has been
                         //   compiled.
                         //////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// To be removed, testing a variety of distributed workflows
///////////////////////////////////////////////////////////////////////////////
//                         outputFileStream<<"#rm -rf $LOCALTEMPDIR"<<endl; // make sure that the local temp directory is clear prior to starting
//                         outputFileStream<<"mkdir $LOCALTEMPDIR"<<endl; // make sure that the local temp directory is recreated for each job
//                         outputFileStream<<"cd $SHAREDINSTALLDIR"<<endl;
//                         outputFileStream<<"$SHAREDINSTALLDIR/product/start.sh"<<endl;
//                         outputFileStream<<"cp $SHAREDINPUTDATADIR/$FILENAME $LOCALTEMPDIR/$FILENAME"<<endl;
//                         outputFileStream<<"for SECTIONNUMBER in `seq $SECTIONNUMBERSTART $SECTIONNUMBEREND`"<<endl;
//                         outputFileStream<<"do"<<endl;
//                         outputFileStream<<"echo $SECTIONNUMBER"<<endl;
//                         outputFileStream<<"# Run MZDASoft Level 1,2,3 processing"<<endl;
//                         outputFileStream<<"$LOCALTEMPDIR/mzdaWorkflowCluster_$SECTIONNUMBER.sh"<<endl;
//                         outputFileStream<<"# Copy results back to shared directory"<<endl;
//                         outputFileStream<<"cp $LOCALTEMPDIR/mzdaLevelThree.$SECTIONNUMBER.*.dat $SHAREDOUTPUTDATADIR"<<endl;
//                         outputFileStream<<"#Cleanup intermediate files in local temp space: A critical feature for scalability acroos any amount of data"<<endl;
//                         outputFileStream<<"#rm $LOCALTEMPDIR/mzdaLevelThree.$SECTIONNUMBER.*.dat"<<endl;
//                         outputFileStream<<"#rm -rf $LOCALTEMPDIR/$SECTIONNUMBER"<<endl;
//                         outputFileStream<<"done"<<endl;
//                         outputFileStream<<"# remove the mzXML file from the local temp directory"<<endl;
//                         outputFileStream<<"#rm $LOCALTEMPDIR/$FILENAME"<<endl;
//                         outputFileStream<<"# remove the temp working directory"<<endl;
//                         outputFileStream<<"#rm -rf $LOCALTEMPDIR"<<endl;
///////////////////////////////////////////////////////////////////////////////
                         //////////////////////////////////////////////////////
                         //////////////////////////////////////////////////////

                         assert(currentSectionEnd-currentSectionStart+1 == sectioncountDefaultMerged.at(i));
                         // debug output:
                         //cout<<"currentdiff="<<currentSectionEnd-currentSectionStart+1<<"sectioncount="<<sectioncountDefaultMerged.at(i)<<endl;

                         // update the end point
                         currentSectionStart = currentSectionStart + sectioncountDefaultMerged.at(i);

                     }
                     else{
                         // output file is not open, no need to try to close it
                         return -2;
                     }
                 }
                 catch (std::exception& e) {
                     std::cout << e.what() <<std::endl;
                     std::cout <<"mzda: exception while trying to open output merged parameter file"<<std::endl;
                     // make sure to close the file stream
                     outputFileStream.close();
                 }  // end trycatch block

                 // make sure to close the output file
                 outputFileStream.close();






                 //////////////////////////////////////////////////////////////
                 // Create the job files for a cluster
                 // - Current code generates job files for the
                 //   TACC Stampede Cluster.
                 //
                 //
                 //#!/bin/bash
                 //#SBATCH -J mzdasoftfile1
                 //#SBATCH -o mzdasoft.o%j
                 //#SBATCH -N 1
                 //#SBATCH -n 16
                 //#SBATCH -p development
                 //#SBATCH -t 01:00:00
                 //$WORK/mzdataccwork/taccJob1.sh
                 //////////////////////////////////////////////////////////////
                 string outputClusterJobScriptFilename= outputClusterJobsFilenameBase.at(0) + fileNumber + ".job";

                 outputFileStream.open(outputClusterJobScriptFilename.c_str(), std::ios::out);
                 try{
                     if( outputFileStream.is_open()) {

                         // The start and end points
                         // update the files section end point
                         currentSectionEnd = currentSectionStart + sectioncountDefaultMerged.at(i) -1;

                         //////////////////////////////////////////////////////
                         // These HPC cluster specific items
                         //////////////////////////////////////////////////////
                         vector<string> nodesPerJob = paraObjPartitions.getParameterValuesAsString("nodesPerJob");
                         vector<string> coresPerJob = paraObjPartitions.getParameterValuesAsString("coresPerJob");
                         vector<string> queueName = paraObjPartitions.getParameterValuesAsString("queueName");
                         vector<string> queueNamePartTwo = paraObjPartitions.getParameterValuesAsString("queueNamePartTwo");
                         vector<string> maxRuntimePerJob = paraObjPartitions.getParameterValuesAsString("maxRuntimePerJob");
                         vector<string> parallelEnvironment = paraObjPartitions.getParameterValuesAsString("parallelEnvironment");
                         vector<string> gridManager = paraObjPartitions.getParameterValuesAsString("gridManager");
                         vector<string> projectAllocation = paraObjPartitions.getParameterValuesAsString("projectAllocation");

                         string nodesPerJobDefault = "1";      // 1 system is default
                         string coresPerJobDefault = "16";     // 16 cores default
                         string queueNameDefault = "normal";   // normal queue
                         string queueNamePartTwoDefault = "ib_only"; // sge cluster only
                         string maxRuntimePerJobDefault = "02:00:00"; // 2 hours
                         string parallelEnvironmentDefault = "orte2";  // sge cluster only
                         string gridManagerDefault = "slurm";
                         string projectAllocationDefault = "";
                         //////////////////////////////////////////////////////
                         //
                         // Extract the parameters from the cluster
                         // configuration file.
                         //
                         // To run on a different or new cluster, more
                         // parameters can be added here....
                         //
                         // This provides flexibility/ adaptibility to
                         // run on novel clusters.
                         //
                         //////////////////////////////////////////////////////
                         if ( nodesPerJob.size() == 1 &&
                              coresPerJob.size() == 1 &&
                              queueName.size() == 1 &&
                              queueNamePartTwo.size() == 1 &&
                              maxRuntimePerJob.size() == 1 &&
                              parallelEnvironment.size() == 1 &&
                              gridManager.size() == 1 &&
                              projectAllocation.size() == 1 ) {

                             //////////////////////////////////////////////////
                             //  If user has specified parameters
                             //  within the Cluster
                             //////////////////////////////////////////////////
                             nodesPerJobDefault = nodesPerJob.at(0);
                             coresPerJobDefault = coresPerJob.at(0);
                             queueNameDefault = queueName.at(0);
                             queueNamePartTwoDefault = queueNamePartTwo.at(0);
                             maxRuntimePerJobDefault = maxRuntimePerJob.at(0);
                             parallelEnvironmentDefault = parallelEnvironment.at(0);
                             gridManagerDefault = gridManager.at(0);
                             projectAllocationDefault = projectAllocation.at(0);
                         }


                         //////////////////////////////////////////////////////
                         // SLURM Job Manager
                         //////////////////////////////////////////////////////
                         if( gridManagerDefault == "slurm" ) {
                             outputFileStream<<"#!/bin/bash"<<endl;
                             outputFileStream<<"#SBATCH -J mzda"<<fileNumber<<endl;
                             outputFileStream<<"#SBATCH -o mzda"<<fileNumber<<".o%j"<<endl;
                             outputFileStream<<"#SBATCH -N "<<nodesPerJobDefault<<endl; // Number of nodes with 1 job
                             outputFileStream<<"#SBATCH -n "<<coresPerJobDefault<<endl; // Number of cores on a cluster node
                             outputFileStream<<"#SBATCH -p "<<queueNameDefault<<endl;  // queue: normal, development
                             outputFileStream<<"#SBATCH -t "<<maxRuntimePerJobDefault<<endl; // max runtime: hh
                             outputFileStream<<"#SBATCH -A "<<projectAllocationDefault<<endl; // project allocation information required if multiple projects
                             outputFileStream<<outputClusterJobsFilenameBase.at(0) + fileNumber + ".sh"<<endl;
                         }

                         //////////////////////////////////////////////////////
                         // Grid Engine Job Manager
                         //////////////////////////////////////////////////////
                         if( gridManagerDefault == "sge" ) {
                             outputFileStream<<"#!/bin/bash"<<endl;
                             outputFileStream<<"#$ -N mzda"<<fileNumber<<endl;
                             outputFileStream<<"#$ -j y"<<endl;
                             outputFileStream<<"#$ -cwd"<<endl;
                             outputFileStream<<"#$ -pe "<<parallelEnvironmentDefault<<" "<<coresPerJobDefault<<endl; // Keep 8 cores as the default.
                             outputFileStream<<"#$ -q "<<queueNameDefault<<endl;
                             outputFileStream<<"#$ -l "<<queueNamePartTwoDefault<<endl;
                             outputFileStream<<outputClusterJobsFilenameBase.at(0) + fileNumber + ".sh"<<endl;
                             outputFileStream<<"exit 0"<<endl;
                         }



                         //////////////////////////////////////////////////////
                         // Add new job managers here.
                         //////////////////////////////////////////////////////





                     }
                     else{
                         // output file is not open, no need to try to close it
                         return -2;
                     }
                 }
                 catch (std::exception& e) {
                     std::cout << e.what() <<std::endl;
                     std::cout <<"mzda: exception while trying to open cluster job script file"<<std::endl;
                     // make sure to close the file stream
                     outputFileStream.close();
                 }  // end trycatch block

                 // make sure to close the output file
                 outputFileStream.close();
                 //////////////////////////////////////////////////////////////
                 // End generating the TACC Stampede Cluster
                 // specific job scripts.
                 //////////////////////////////////////////////////////////////



                 //////////////////////////////////////////////////////////////
                 // START
                 // Add new HPC job workflows here:  04/24/14
                 // The following block of code adds the TPP XTandem
                 // Workflow to the MZDASoft HPC Pipeline.
                 //
                 //
                 //
                 //////////////////////////////////////////////////////////////


                 //////////////////////////////////////////////////////////////
                 // TPPXtandem Workflow: Control Files for the XTandem
                 // Workflow.
                 // Input XML Files:
                 //
                 // input_searchGUI_MZXMLFILENAMEBASE_TPPXTandem#.xml
                 //
                 //                 <?xml version="1.0"?>
                 //                 <bioml>
                 //                     <note type="input" label="list path, default parameters">parameters_searchGUI.xml</note>  --> We need to make this an input parameter
                 //                     <note type="input" label="list path, taxonomy information">taxonomy_searchGUI_MZXMLFILENAMEBASE.xml</note>
                 //                     <note type="input" label="protein, taxon">all</note>
                 //                     <note type="input" label="spectrum, path">FULLPATHTOMZXMLFILE.mzXML</note>  --> already an input parameter
                 //                     <note type="input" label="output, path">FULLPATHTOMZXMLFILE.xml</note>
                 //                 </bioml>
                 //
                 // taxonomy_searchGUI_MZXMLFILENAMEBASE_TPPXTandem#.xml
                 //                 <?xml version="1.0"?>
                 //                 <bioml label="x! taxon-to-file matching list">
                 //                     <taxon label="all">
                 //                         <file format="peptide" URL="FULLPATHTOFASTADATABASE.fasta" />  --> This needs to be an input parameter( 1 per processing run )
                 //                     </taxon>
                 //                 </bioml>
                 //
                 //////////////////////////////////////////////////////////////



                 //////////////////////////////////////////////////////////////
                 // Generate the taxonomy xml file for the TPP XTandem
                 // Workflow.  The taxonomy xml file contains the
                 // reference to the search database fasta file to be used.
                 //
                 //////////////////////////////////////////////////////////////
                 cout<<"Logging Point 1, 042814"<<endl;
                 std::ofstream ofsTPPXTandemTaxonomyXML;
                 stringstream ssTempTPPXTandemTaxonomyXML;
                 ssTempTPPXTandemTaxonomyXML<<(i+1);
                 string fileNumberTPPXTandemTaxonomyXML = ssTempTPPXTandemTaxonomyXML.str();
                 string outputClusterTPPXTandemTaxonomyXMLFile = outputClusterJobsFilenameBase.at(0) + fileNumberTPPXTandemTaxonomyXML + "TaxonomyTPPXTandem"  + ".xml";
                 ofsTPPXTandemTaxonomyXML.open(outputClusterTPPXTandemTaxonomyXMLFile.c_str(), std::ios::out);
                 try{
                     if( ofsTPPXTandemTaxonomyXML.is_open() ) {
                         ofsTPPXTandemTaxonomyXML<<"<?xml version=\"1.0\"?>"<<endl;
                         ofsTPPXTandemTaxonomyXML<<" <bioml label=\"x! taxon-to-file matching list\">"<<endl;
                         ofsTPPXTandemTaxonomyXML<<"  <taxon label=\"all\">"<<endl;
                         ofsTPPXTandemTaxonomyXML<<"   <file format=\"peptide\" URL=\""+ tppXTandemDatabaseFastaFile.at(0) + "\" />"<<endl;
                         ofsTPPXTandemTaxonomyXML<<"  </taxon>"<<endl;
                         ofsTPPXTandemTaxonomyXML<<" </bioml>"<<endl;

                     }
                     else{
                         // output file is not open, no need to try to close it
                         return -2;
                     }
                 }
                 catch (std::exception& e) {
                     std::cout << e.what() <<std::endl;
                     std::cout <<"mzda: exception while trying to open output merged parameter file"<<std::endl;
                     // make sure to close the file stream
                     ofsTPPXTandemTaxonomyXML.close();
                 }  // end trycatch block
                 ofsTPPXTandemTaxonomyXML.close(); // Make sure we close the file in normal case


                 //////////////////////////////////////////////////////////////
                 //
                 // Generate the input parameter file for the TPP XTandem
                 // Workflow.
                 //
                 //////////////////////////////////////////////////////////////
                 cout<<"Logging Point 2, 042814"<<endl;
                 std::ofstream ofsTPPXTandemInputXML;
                 stringstream ssTempTPPXTandemInputXML;
                 ssTempTPPXTandemInputXML<<(i+1);
                 string fileNumberTPPXTandemInputXML = ssTempTPPXTandemInputXML.str();
                 string outputClusterTPPXTandemInputXMLFile =  outputClusterJobsFilenameBase.at(0) + fileNumberTPPXTandemInputXML +"InputTPPXTandem" + ".xml";
                 ofsTPPXTandemInputXML.open(outputClusterTPPXTandemInputXMLFile.c_str(), std::ios::out);
                 try{
                     if( ofsTPPXTandemInputXML.is_open() ) {
                         ofsTPPXTandemInputXML<<"<?xml version=\"1.0\"?>"<<endl;
                         ofsTPPXTandemInputXML<<" <bioml>"<<endl;
                         ofsTPPXTandemInputXML<<"  <note type=\"input\" label=\"list path, default parameters\">" + tppXTandemParameterFile.at(0) + "</note>"<<endl;  // We need to make this an input parameter
                         ofsTPPXTandemInputXML<<"  <note type=\"input\" label=\"list path, taxonomy information\">"+ outputClusterTPPXTandemTaxonomyXMLFile+ "</note>"<<endl;
                         ofsTPPXTandemInputXML<<"  <note type=\"input\" label=\"protein, taxon\">all</note>"<<endl;
                         ofsTPPXTandemInputXML<<"  <note type=\"input\" label=\"spectrum, path\">"+ localTempDir.at(0)+ "/" + mzxmlfilenamebaseDefaultMerged.at(i) +"</note>"<<endl;


                         //////////////////////////////////////////////////////
                         // here we need to change the extension for the output of the TPPXTandem result, we simply
                         // need to change the extension from .mzXML to xml
                         //////////////////////////////////////////////////////

                         string temp = mzxmlfilenamebaseDefaultMerged.at(i);
                         string temp2 = "mzXML";
                         temp.replace(temp.find("mzXML"),temp2.length(),"xml");
                         ofsTPPXTandemInputXML<<"  <note type=\"input\" label=\"output, path\">" + tppXTandemLocalTempDir.at(0) + "/" + temp +"</note>"<<endl;  // output to tppxtandem temp dir instead of global mzdasoft tempdir
                         ofsTPPXTandemInputXML<<" </bioml>"<<endl;
                     }
                     else{
                         // output file is not open, no need to try to close it
                         return -2;
                     }
                 }
                 catch (std::exception& e) {
                     std::cout << e.what() <<std::endl;
                     std::cout <<"mzda: exception while trying to open output merged parameter file"<<std::endl;
                     // make sure to close the file stream
                     ofsTPPXTandemInputXML.close();
                 }  // end trycatch block
                 ofsTPPXTandemInputXML.close(); // Make sure we close the file in normal case


                 //////////////////////////////////////////////////////////////
                 // Shell Script Run Manager: This script controls all the
                 // work that will be done by the Tandem MS search engine
                 // workflow for TPPXTandem.  As new search engines are
                 // added, a new workflow simply needs to be created
                 // as shown below.
                 // Note: Some new workflows may need new parameters to
                 // be available on the
                 //////////////////////////////////////////////////////////////
                 cout<<"Logging Point 3, 042814"<<endl;
                 std::ofstream outputFileStreamTPPXTandemScript;
                 stringstream ssTPPXTandemScript;
                 ssTPPXTandemScript<<(i+1); // integer to string conversion
                 string fileNumberTPPXTandemScript = ssTPPXTandemScript.str(); // converting stringstream to string object

                 string outputClusterJobFilenameTPPXTandemScript= outputClusterJobsFilenameBase.at(0) + fileNumberTPPXTandemScript + "TPPXTandem" +  ".sh";

                 outputFileStreamTPPXTandemScript.open(outputClusterJobFilenameTPPXTandemScript.c_str(), std::ios::out);
                 try{
                     if( outputFileStreamTPPXTandemScript.is_open()) {
                         //////////////////////////////////////////////////////
                         // Generate the data contained in the TPPXTandem Workflow
                         //#!/bin/bash
                         //FASTADATABASE=/scratch/01759/nramirez/work/updated.ipi.HUMAN.v3.83.fasta
                         //TPPBINDIR=/work/01759/nramirez/mzdatpp/TPP-4.6.3/trans_proteomic_pipeline/build/CentOS-x86_64
                         //TANDEMEXE=$TPPBINDIR/tandem
                         //TANDEM2XML=$TPPBINDIR/Tandem2XML
                         //INTERACTPARSER=$TPPBINDIR/InteractParser
                         //DATABASEPARSER=$TPPBINDIR/DatabaseParser
                         //REFRESHPARSER=$TPPBINDIR/RefreshParser
                         //PEPTIDEPROPHETPARSER=$TPPBINDIR/PeptideProphetParser
                         //PROPHETMODELS=$TPPBINDIR/ProphetModels.pl
                         //#Next we need to update the path to the spectra file, should be changed to point to .mzXML file: input_searchGUI.xml
                         //$TANDEMEXE /scratch/01759/nramirez/work/input_searchGUI_$1.xml
                         //mv /scratch/01759/nramirez/work/$1.*.t.xml  /scratch/01759/nramirez/work/$1.tandem.xml
                         //# Following the Tandem Search we run through the TPP Pipeline process up to PeptideProphet.
                         //#$TANDEM2XML /scratch/01759/nramirez/work/$1.tandem.xml /scratch/01759/nramirez/work/$1.tandem.pepXML
                         //# InteractParser is adding the [] modification fields.
                         //#$INTERACTPARSER /scratch/01759/nramirez/work/$1.tandem.interact.pepXML /scratch/01759/nramirez/work/$1.tandem.pepXML
                         //# The DatabaseParser utility prints the full path and filename of database references in a pep.xml file
                         //#$DATABASEPARSER /scratch/01759/nramirez/work/$1.tandem.interact.pepXML
                         //# The RefreshParser add Protein mappings to each Peptide using the .fasta database file
                         //#$REFRESHPARSER /scratch/01759/nramirez/work/$1.tandem.interact.pepXML $FASTADATABASE
                         //# Prepares to run PeptideProphet
                         //#$PEPTIDEPROPHETPARSER /scratch/01759/nramirez/work/$1.tandem.interact.pepXML
                         //# Runs PeptideProphet
                         //#$PROPHETMODELS -i /scratch/01759/nramirez/work/$1.tandem.interact.pepXML
                         ///////////////////////////////////////////////////////
                         outputFileStreamTPPXTandemScript<<"#!/bin/bash"<<endl;
                         outputFileStreamTPPXTandemScript<<"echo \"File "<<fileNumber<<" : Start Processing TPPXTandem Workflow\""<<endl;
                         outputFileStreamTPPXTandemScript<<"FASTADATABASE="<<tppXTandemDatabaseFastaFile.at(0)<<endl;
                         outputFileStreamTPPXTandemScript<<"TPPBINDIR="<<tppXTandemBinDir.at(0)<<endl;
                         outputFileStreamTPPXTandemScript<<"TANDEMEXE="<<"$TPPBINDIR/tandem"<<endl;
                         outputFileStreamTPPXTandemScript<<"TANDEM2XML="<<"$TPPBINDIR/Tandem2XML"<<endl;
                         outputFileStreamTPPXTandemScript<<"INTERACTPARSER="<<"$TPPBINDIR/InteractParser"<<endl;
                         outputFileStreamTPPXTandemScript<<"DATABASEPARSER="<<"$TPPBINDIR/DatabaseParser"<<endl;
                         outputFileStreamTPPXTandemScript<<"REFRESHPARSER="<<"$TPPBINDIR/RefreshParser"<<endl;
                         outputFileStreamTPPXTandemScript<<"PEPTIDEPROPHETPARSER="<<"$TPPBINDIR/PeptideProphetParser"<<endl;
                         outputFileStreamTPPXTandemScript<<"PROPHETMODELS="<<"$TPPBINDIR/ProphetModels.pl"<<endl;


                         //////////////////////////////////////////////////////
                         //Setup a new directory structure that is
                         //completely separate from the Parallel Peak
                         //Extractor.
                         //
                         //////////////////////////////////////////////////////


                         //////////////////////////////////////////////////////
                         // Running the XTandem executable
                         //
                         //////////////////////////////////////////////////////
                         outputFileStreamTPPXTandemScript<<"echo \"File "<<fileNumber<<" : Running XTandem Search Engine\""<<endl;
                         outputFileStreamTPPXTandemScript<<"$TANDEMEXE "<<outputClusterTPPXTandemInputXMLFile<<endl;  // The input parameter file

                         //////////////////////////////////////////////////////
                         // Dynamically generate the rename command
                         // to rename the output of the XTandem search engine.
                         // The XTandem search engine automatically adds
                         // the current date to the filename.
                         // We need to perform a series of transformations
                         // to have a consistent naming scheme for the
                         // output of TPPXTandem.
                         // mv /scratch/01759/nramirez/work/$1.*.t.xml  /scratch/01759/nramirez/work/$1.tandem.xml
                         //////////////////////////////////////////////////////
                         string temp1 = mzxmlfilenamebaseDefaultMerged.at(i);
                         string temp2 = "mzXML";
                         temp1.replace(temp1.find("mzXML"),temp2.length(),"*.t.xml");
                         string sourceFilePattern = temp1;  // $1.*.t.xml

                         string temp3 = mzxmlfilenamebaseDefaultMerged.at(i);
                         string temp4 = "mzXML";
                         temp3.replace(temp3.find("mzXML"),temp4.length(),"tandem.xml");
                         string destFilePattern = temp3;  // $1.tandem.xml
                         // generate a file that no longer has the date time stamp within the filename
                         outputFileStreamTPPXTandemScript<<"echo \"File "<<fileNumber<<" : Running XTandem results file renaming process\""<<endl;
                         outputFileStreamTPPXTandemScript<<"mv "<<tppXTandemLocalTempDir.at(0)+"/"+sourceFilePattern<<" "<<tppXTandemLocalTempDir.at(0)<<"/"<<destFilePattern<<endl;

                         //////////////////////////////////////////////////////
                         // Dynamically generate the command to convert
                         // the output of the XTandem Search Engine
                         // to pepXML format
                         // //#$TANDEM2XML /scratch/01759/nramirez/work/$1.tandem.xml /scratch/01759/nramirez/work/$1.tandem.pepXML
                         //////////////////////////////////////////////////////
                         temp1 = mzxmlfilenamebaseDefaultMerged.at(i);
                         temp2 = "mzXML";
                         temp1.replace(temp1.find("mzXML"),temp2.length(),"tandem.xml");
                         sourceFilePattern = temp1;  //$1.tandem.xml

                         temp3 = mzxmlfilenamebaseDefaultMerged.at(i);
                         temp4 = "mzXML";
                         temp3.replace(temp3.find("mzXML"),temp4.length(),"tandem.pepXML");
                         destFilePattern = temp3;  // $1.tandem.pepXML
                         // generate the pep XML file
                         outputFileStreamTPPXTandemScript<<"echo \"File "<<fileNumber<<" : Running Tandem2XML\""<<endl;
                         outputFileStreamTPPXTandemScript<<"$TANDEM2XML "<<tppXTandemLocalTempDir.at(0)+"/"+sourceFilePattern<<" "<<tppXTandemLocalTempDir.at(0)<<"/"<<destFilePattern<<endl;

                         //////////////////////////////////////////////////////
                         //# InteractParser is adding the [] modification fields.
                         //
                         //
                         //////////////////////////////////////////////////////
                         //#$INTERACTPARSER /scratch/01759/nramirez/work/$1.tandem.interact.pepXML /scratch/01759/nramirez/work/$1.tandem.pepXML
                         temp1 = mzxmlfilenamebaseDefaultMerged.at(i);
                         temp2 = "mzXML";
                         temp1.replace(temp1.find("mzXML"),temp2.length(),"tandem.interact.pepXML");
                         sourceFilePattern = temp1;  //tandem.interact.pepXML

                         temp3 = mzxmlfilenamebaseDefaultMerged.at(i);
                         temp4 = "mzXML";
                         temp3.replace(temp3.find("mzXML"),temp4.length(),"tandem.pepXML");
                         destFilePattern = temp3;  // $1.tandem.pepXML
                         // generate the pep XML file
                         outputFileStreamTPPXTandemScript<<"echo \"File "<<fileNumber<<" : Running InteractParser\""<<endl;
                         outputFileStreamTPPXTandemScript<<"$INTERACTPARSER "<<tppXTandemLocalTempDir.at(0)+"/"+sourceFilePattern<<" "<<tppXTandemLocalTempDir.at(0)<<"/"<<destFilePattern<<endl;


                         //////////////////////////////////////////////////////
                         //# DatabaseParser
                         //
                         //
                         //////////////////////////////////////////////////////
                         //#$DATABASEPARSER /scratch/01759/nramirez/work/$1.tandem.interact.pepXML
                         temp1 = mzxmlfilenamebaseDefaultMerged.at(i);
                         temp2 = "mzXML";
                         temp1.replace(temp1.find("mzXML"),temp2.length(),"tandem.interact.pepXML");
                         sourceFilePattern = temp1;  //$tandem.interact.pepXML
                         outputFileStreamTPPXTandemScript<<"echo \"File "<<fileNumber<<" : Running DatabaseParser\""<<endl;
                         outputFileStreamTPPXTandemScript<<"$DATABASEPARSER "<<tppXTandemLocalTempDir.at(0)+"/"+sourceFilePattern<<endl;

                         //////////////////////////////////////////////////////
                         //#$REFRESHPARSER /scratch/01759/nramirez/work/$1.tandem.interact.pepXML $FASTADATABASE
                         //
                         //
                         //////////////////////////////////////////////////////
                         temp1 = mzxmlfilenamebaseDefaultMerged.at(i);
                         temp2 = "mzXML";
                         temp1.replace(temp1.find("mzXML"),temp2.length(),"tandem.interact.pepXML");
                         sourceFilePattern = temp1;  //$tandem.interact.pepXML
                         outputFileStreamTPPXTandemScript<<"echo \"File "<<fileNumber<<" : Running RefreshParser\""<<endl;
                         outputFileStreamTPPXTandemScript<<"$REFRESHPARSER "<<tppXTandemLocalTempDir.at(0)+"/"+sourceFilePattern<<" $FASTADATABASE"<<endl;

                         //////////////////////////////////////////////////////
                         //#$PEPTIDEPROPHETPARSER /scratch/01759/nramirez/work/$1.tandem.interact.pepXML
                         //#$PROPHETMODELS -i /scratch/01759/nramirez/work/$1.tandem.interact.pepXML
                         //
                         //
                         //////////////////////////////////////////////////////
                         temp1 = mzxmlfilenamebaseDefaultMerged.at(i);
                         temp2 = "mzXML";
                         temp1.replace(temp1.find("mzXML"),temp2.length(),"tandem.interact.pepXML");
                         sourceFilePattern = temp1;  //$tandem.interact.pepXML
                         outputFileStreamTPPXTandemScript<<"echo \"File "<<fileNumber<<" : Running PeptideProphetParser\""<<endl;
                         outputFileStreamTPPXTandemScript<<"$PEPTIDEPROPHETPARSER "<<tppXTandemLocalTempDir.at(0)+"/"+sourceFilePattern<<endl;

                         temp1 = mzxmlfilenamebaseDefaultMerged.at(i);
                         temp2 = "mzXML";
                         temp1.replace(temp1.find("mzXML"),temp2.length(),"tandem.interact.pepXML");
                         sourceFilePattern = temp1;  //$tandem.interact.pepXML
                         outputFileStreamTPPXTandemScript<<"echo \"File "<<fileNumber<<" : Running ProphetModels\""<<endl;
                         outputFileStreamTPPXTandemScript<<"$PROPHETMODELS -i "<<tppXTandemLocalTempDir.at(0)+"/"+sourceFilePattern<<endl;




                     }
                     else{
                         // output file is not open, no need to try to close it
                         return -2;
                     }
                 }
                 catch (std::exception& e) {
                     std::cout << e.what() <<std::endl;
                     std::cout <<"mzda: exception while trying to open output merged parameter file"<<std::endl;
                     // make sure to close the file stream
                     outputFileStreamTPPXTandemScript.close();
                 }  // end trycatch block

                 // make sure to close the output file
                 outputFileStreamTPPXTandemScript.close();
                 //////////////////////////////////////////////////////////////
                 // Job scipt file generation.  For each run script we need
                 // to dynamically create the job scheduler output for
                 // each type of job scheduler supported.
                 //////////////////////////////////////////////////////////////
                 std::ofstream outputFileStreamTPPXTandemJobFile;
                 stringstream ssTPPXTandemJobFile;
                 ssTPPXTandemJobFile<<(i+1); // integer to string conversion
                 string fileNumberTPPXTandemJobFile = ssTPPXTandemJobFile.str(); // converting stringstream to string object


                 //////////////////////////////////////////////////////////////
                 // HPC Job Manager Job File
                 //////////////////////////////////////////////////////////////
                 //////////////////////////////////////////////////////////////
                 // These are tacc specific items
                 //////////////////////////////////////////////////////////////
                 vector<string> nodesPerJob = paraObjPartitions.getParameterValuesAsString("nodesPerJob");
                 vector<string> coresPerJob = paraObjPartitions.getParameterValuesAsString("coresPerJob");
                 vector<string> queueName = paraObjPartitions.getParameterValuesAsString("queueName");
                 vector<string> queueNamePartTwo = paraObjPartitions.getParameterValuesAsString("queueNamePartTwo");
                 vector<string> maxRuntimePerJob = paraObjPartitions.getParameterValuesAsString("maxRuntimePerJob");
                 vector<string> parallelEnvironment = paraObjPartitions.getParameterValuesAsString("parallelEnvironment");
                 vector<string> gridManager = paraObjPartitions.getParameterValuesAsString("gridManager");
                 vector<string> projectAllocation = paraObjPartitions.getParameterValuesAsString("projectAllocation");


                 string nodesPerJobDefault = "1";      // 1 system is default
                 string coresPerJobDefault = "16";     // 16 cores default
                 string queueNameDefault = "normal";   // normal queue
                 string queueNamePartTwoDefault = "ib_only"; // sge cluster only
                 string maxRuntimePerJobDefault = "02:00:00"; // 2 hours
                 string parallelEnvironmentDefault = "orte2";  // sge cluster only
                 string gridManagerDefault = "slurm";
                 string projectAllocationDefault = "";
                 //////////////////////////////////////////////////////////////
                 //
                 // Extract the parameters from the cluster
                 // configuration file.
                 //
                 // To run on a different or new cluster, more
                 // parameters can be added here....
                 //
                 // This provides flexibility/ adaptibility to
                 // run on novel clusters.
                 //
                 //////////////////////////////////////////////////////////////
                 if ( nodesPerJob.size() == 1 &&
                      coresPerJob.size() == 1 &&
                      queueName.size() == 1 &&
                      queueNamePartTwo.size() == 1 &&
                      maxRuntimePerJob.size() == 1 &&
                      parallelEnvironment.size() == 1 &&
                      gridManager.size() == 1 &&
                      projectAllocation.size() == 1 ) {

                     //////////////////////////////////////////////////
                     //  If user has specified parameters
                     //  within the Cluster
                     //////////////////////////////////////////////////
                     nodesPerJobDefault = nodesPerJob.at(0);
                     coresPerJobDefault = coresPerJob.at(0);
                     queueNameDefault = queueName.at(0);
                     queueNamePartTwoDefault = queueNamePartTwo.at(0);
                     maxRuntimePerJobDefault = maxRuntimePerJob.at(0);
                     parallelEnvironmentDefault = parallelEnvironment.at(0);
                     gridManagerDefault = gridManager.at(0);
                     projectAllocationDefault = projectAllocation.at(0);
                 }

                 string outputClusterJobFilenameTPPXTandemJobFile= outputClusterJobsFilenameBase.at(0) +  fileNumberTPPXTandemJobFile + "TPPXTandem"  + ".job";
                 outputFileStreamTPPXTandemJobFile.open(outputClusterJobFilenameTPPXTandemJobFile.c_str(), std::ios::out);
                 try{
                     if( outputFileStreamTPPXTandemJobFile.is_open()) {
                         //////////////////////////////////////////////////////
                         // Generate the job file for the TPPXTandem Workflow
                         //////////////////////////////////////////////////////

                         //////////////////////////////////////////////////////
                         // SLURM Job manager
                         //////////////////////////////////////////////////////
                         if( gridManagerDefault == "slurm" ) {
                             outputFileStreamTPPXTandemJobFile<<"#!/bin/bash"<<endl;
                             outputFileStreamTPPXTandemJobFile<<"#SBATCH -J mzdaTPP"<<fileNumberTPPXTandemJobFile<<endl;
                             outputFileStreamTPPXTandemJobFile<<"#SBATCH -o mzdaTPP"<<fileNumberTPPXTandemJobFile<<".o%j"<<endl;
                             outputFileStreamTPPXTandemJobFile<<"#SBATCH -N "<<nodesPerJobDefault<<endl; // Number of nodes with 1 job
                             outputFileStreamTPPXTandemJobFile<<"#SBATCH -n "<<coresPerJobDefault<<endl; // Number of cores on a cluster node
                             outputFileStreamTPPXTandemJobFile<<"#SBATCH -p "<<queueNameDefault<<endl;  // queue: normal, development
                             outputFileStreamTPPXTandemJobFile<<"#SBATCH -t "<<maxRuntimePerJobDefault<<endl; // max runtime: hh
                             outputFileStreamTPPXTandemJobFile<<"#SBATCH -A "<<projectAllocationDefault<<endl; // project allocation information required if multiple projects
                             outputFileStreamTPPXTandemJobFile<<outputClusterJobFilenameTPPXTandemScript<<endl;  // the name of the runtime script
                         }

                         //////////////////////////////////////////////////////
                         // Grid Engine Job Manager
                         //////////////////////////////////////////////////////
                         if( gridManagerDefault == "sge" ) {
                             outputFileStreamTPPXTandemJobFile<<"#!/bin/bash"<<endl;
                             outputFileStreamTPPXTandemJobFile<<"#$ -N mzdaTPP"<<fileNumberTPPXTandemJobFile<<endl;
                             outputFileStreamTPPXTandemJobFile<<"#$ -j y"<<endl;
                             outputFileStreamTPPXTandemJobFile<<"#$ -cwd"<<endl;
                             outputFileStreamTPPXTandemJobFile<<"#$ -pe "<<parallelEnvironmentDefault<<" "<<coresPerJobDefault<<endl; // Keep 8 cores as the default.
                             outputFileStreamTPPXTandemJobFile<<"#$ -q "<<queueNameDefault<<endl;
                             outputFileStreamTPPXTandemJobFile<<"#$ -l "<<queueNamePartTwoDefault<<endl;
                             outputFileStreamTPPXTandemJobFile<<outputClusterJobFilenameTPPXTandemScript<<endl;  // the name of the runtime script
                             outputFileStreamTPPXTandemJobFile<<"exit 0"<<endl;
                         }

                         //////////////////////////////////////////////////////
                         // Add new job managers here.
                         //////////////////////////////////////////////////////


                     }
                     else{
                         // output file is not open, no need to try to close it
                         return -2;
                     }
                 }
                 catch (std::exception& e) {
                     std::cout << e.what() <<std::endl;
                     std::cout <<"mzda: exception while trying to open output merged parameter file"<<std::endl;
                     // make sure to close the file stream
                     outputFileStreamTPPXTandemJobFile.close();
                 }  // end trycatch block
                 // make sure to close the output file
                 outputFileStreamTPPXTandemJobFile.close();



                 //////////////////////////////////////////////////////////////
                 // Quantification Workflow.
                 //
                 // The Quantification workflow uses a Matlab Compiled
                 // module, that requires all the MZDASoft Parallel Peak
                 // Extractor results to be present as well as all the
                 // Tandem results in .pepXML format.
                 // - In a future release we can create a dependency
                 //   workflow where the user would not need to
                 //   submit the Quantification Workflow after the
                 //   main workflow.  For version 1.0.0 we have
                 //   this as an independent workflow.
                 //
                 // Only 1 job script is needed for this workflow
                 // A single job will carry out the processing of the
                 // Quantification workflow.  This job can be submitted
                 // after both the Parallel Peak Extractor and the
                 // and the Tandem MS searches have completed.
                 //
                 // In a future release we can parallelize the Quantification
                 // code by partitioning the set of Peptides detected
                 // by the Tandem Database Searches.  Multiple jobs
                 // can be dynamically created so that the set of
                 // Quantifcations can be partitioned and can
                 // proceed in parallel by creating jobs of
                 // disjoint groupings of peptides.
                 //
                 // For this release however, we run an integrated single
                 // job that performs the Quantification workflow.
                 //
                 // For the quantification workflow, we have all the
                 // required information to proceed with the global
                 // parameter file.  Using this global information,
                 // at this point we can dynamically generate everything
                 // that will be needed by the Quantification workflow.
                 //
                 // We will likely need a shell file script and a .job file
                 // We can consider adding job dependencies in a future release.
                 //
                 //
                 //
                 //////////////////////////////////////////////////////////////


                ///////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////
                } // end loop through the total number of input files




                 //////////////////////////////////////////////////////////////
                 //MZDAQuant Pipeline Integration
                 //START BLOCK
                 //////////////////////////////////////////////////////////////
                 //////////////////////////////////////////////////////////////
                 //TPP Merge Process:
                 //mzdaClusterJobTPPMerge.job
                 //////////////////////////////////////////////////////////////
                 string outputClusterJobFilenameTPPMerge=outputClusterJobsFilenameBase.at(0)+"TPPMerge.job";
                 std::ofstream outputFileStreamTPPMerge;
                 outputFileStreamTPPMerge.open(outputClusterJobFilenameTPPMerge.c_str(), std::ios::out);
                 try{
                     if( outputFileStreamTPPMerge.is_open()) {
                         //////////////////////////////////////////////////////
                         // These HPC cluster specific items
                         //////////////////////////////////////////////////////
                         vector<string> nodesPerJob = paraObjPartitions.getParameterValuesAsString("nodesPerJob");
                         vector<string> coresPerJob = paraObjPartitions.getParameterValuesAsString("coresPerJob");
                         vector<string> queueName = paraObjPartitions.getParameterValuesAsString("queueName");
                         vector<string> queueNamePartTwo = paraObjPartitions.getParameterValuesAsString("queueNamePartTwo");
                         vector<string> maxRuntimePerJob = paraObjPartitions.getParameterValuesAsString("maxRuntimePerJob");
                         vector<string> parallelEnvironment = paraObjPartitions.getParameterValuesAsString("parallelEnvironment");
                         vector<string> gridManager = paraObjPartitions.getParameterValuesAsString("gridManager");
                         vector<string> projectAllocation = paraObjPartitions.getParameterValuesAsString("projectAllocation");
                         string nodesPerJobDefault = "1";      // 1 system is default
                         string coresPerJobDefault = "16";     // 16 cores default
                         string queueNameDefault = "normal";   // normal queue
                         string queueNamePartTwoDefault = "ib_only"; // sge cluster only
                         string maxRuntimePerJobDefault = "02:00:00"; // 2 hours
                         string parallelEnvironmentDefault = "orte2";  // sge cluster only
                         string gridManagerDefault = "slurm";
                         string projectAllocationDefault = "";
                         //////////////////////////////////////////////////////
                         //
                         // Extract the parameters from the cluster
                         // configuration file.
                         //
                         // To run on a different or new cluster, more
                         // parameters can be added here....
                         //
                         // This provides flexibility/ adaptibility to
                         // run on novel clusters.
                         //
                         //////////////////////////////////////////////////////
                         if ( nodesPerJob.size() == 1 &&
                              coresPerJob.size() == 1 &&
                              queueName.size() == 1 &&
                              queueNamePartTwo.size() == 1 &&
                              maxRuntimePerJob.size() == 1 &&
                              parallelEnvironment.size() == 1 &&
                              gridManager.size() == 1 &&
                              projectAllocation.size() == 1 ) {
                             //////////////////////////////////////////////////
                             //  If user has specified parameters
                             //  within the Cluster
                             //////////////////////////////////////////////////
                             nodesPerJobDefault = nodesPerJob.at(0);
                             coresPerJobDefault = coresPerJob.at(0);
                             queueNameDefault = queueName.at(0);
                             queueNamePartTwoDefault = queueNamePartTwo.at(0);
                             maxRuntimePerJobDefault = maxRuntimePerJob.at(0);
                             parallelEnvironmentDefault = parallelEnvironment.at(0);
                             gridManagerDefault = gridManager.at(0);
                             projectAllocationDefault = projectAllocation.at(0);
                         }

                         //////////////////////////////////////////////////////
                         // Variables from the global parameter file
                         // needed to dynamically
                         //////////////////////////////////////////////////////
                         string mcrPathTemp = mcrPath.at(0);
                         string sharedInstallDirTemp = sharedInstallDir.at(0);
                         string sharedOutputDataDirTemp = sharedOutputDataDir.at(0);
                         string mzdaQuantInstallDirTemp = mzdaQuantInstallDir.at(0);
                         // todo: add a mzdasoftquantinstallpath input parameter
                         //////////////////////////////////////////////////////
                         // SLURM Job Manager
                         //////////////////////////////////////////////////////
                         if( gridManagerDefault == "slurm" ) {
                             outputFileStreamTPPMerge<<"#!/bin/bash"<<endl;
                             outputFileStreamTPPMerge<<"#SBATCH -J mzdaTPPMerge"<<endl;
                             outputFileStreamTPPMerge<<"#SBATCH -o mzdaTPPMerge"<<".o%j"<<endl;
                             outputFileStreamTPPMerge<<"#SBATCH -N "<<nodesPerJobDefault<<endl; // Number of nodes with 1 job
                             outputFileStreamTPPMerge<<"#SBATCH -n "<<coresPerJobDefault<<endl; // Number of cores on a cluster node
                             outputFileStreamTPPMerge<<"#SBATCH -p "<<queueNameDefault<<endl;  // queue: normal, development
                             outputFileStreamTPPMerge<<"#SBATCH -t "<<maxRuntimePerJobDefault<<endl; // max runtime: hh
                             outputFileStreamTPPMerge<<"#SBATCH -A "<<projectAllocationDefault<<endl; // project allocation information required if multiple projects
                             outputFileStreamTPPMerge<<mzdaQuantInstallDirTemp<<"/run_mzdaQuantMergeTPP.sh "<<mcrPathTemp<<" "<<sharedOutputDataDirTemp<<"/parameters.txt"<<endl;
                         }
                         //////////////////////////////////////////////////////
                         // Grid Engine Job Manager
                         //////////////////////////////////////////////////////
                         if( gridManagerDefault == "sge" ) {
                             outputFileStreamTPPMerge<<"#!/bin/bash"<<endl;
                             outputFileStreamTPPMerge<<"#$ -N mzdaTPPMerge"<<endl;
                             outputFileStreamTPPMerge<<"#$ -j y"<<endl;
                             outputFileStreamTPPMerge<<"#$ -cwd"<<endl;
                             outputFileStreamTPPMerge<<"#$ -pe "<<parallelEnvironmentDefault<<" "<<coresPerJobDefault<<endl; // Keep 8 cores as the default.
                             outputFileStreamTPPMerge<<"#$ -q "<<queueNameDefault<<endl;
                             outputFileStreamTPPMerge<<"#$ -l "<<queueNamePartTwoDefault<<endl;
                             outputFileStreamTPPMerge<<mzdaQuantInstallDirTemp<<"/run_mzdaQuantMergeTPP.sh "<<mcrPathTemp<<" "<<sharedOutputDataDirTemp<<"/parameters.txt"<<endl;
                             outputFileStreamTPPMerge<<"exit 0"<<endl;
                         }
                         //////////////////////////////////////////////////////
                         // Add new job managers here.
                         //////////////////////////////////////////////////////
                     }
                     else{
                         // output file is not open, no need to try to close it
                         return -2;
                     }
                 }
                 catch (std::exception& e) {
                     std::cout << e.what() <<std::endl;
                     std::cout <<"mzda: exception while trying to open cluster job script file"<<std::endl;
                     // make sure to close the file stream
                     outputFileStreamTPPMerge.close();
                 }  // end trycatch block
                 // make sure to close the output file
                 outputFileStreamTPPMerge.close();



                 //////////////////////////////////////////////////////////////
                 //Fraction Level Peptide Quantification:
                 //mzdaClusterJobPeptideQuant#.job
                 //Need to loop through all the fractions
                 //////////////////////////////////////////////////////////////
                 // Here, we need to create 1 job per fraction
                 // Each fraction can be processed independently, in parallel
                 //
                 //////////////////////////////////////////////////////////////

                 int maxFractions = MZDAQUANTMAXFRACTIONS.at(0);  // need to read this variable in from parameter file
                 int fractionCounter=0;

                 for ( fractionCounter = 0; fractionCounter < maxFractions; fractionCounter++ ) {

                 stringstream ssPeptideQuantJobFile;
                 ssPeptideQuantJobFile<<(fractionCounter+1); // integer to string conversion
                 string fileNumberPeptideQuantJobFile = ssPeptideQuantJobFile.str(); // converting stringstream to string object

                 string outputClusterJobFilenamePeptideQuant=outputClusterJobsFilenameBase.at(0)+fileNumberPeptideQuantJobFile+"PeptideQuant.job";
                 std::ofstream outputFileStreamPeptideQuant;
                 outputFileStreamPeptideQuant.open(outputClusterJobFilenamePeptideQuant.c_str(), std::ios::out);
                 try{
                     if( outputFileStreamPeptideQuant.is_open()) {
                         //////////////////////////////////////////////////////
                         // These HPC cluster specific items
                         //////////////////////////////////////////////////////
                         vector<string> nodesPerJob = paraObjPartitions.getParameterValuesAsString("nodesPerJob");
                         vector<string> coresPerJob = paraObjPartitions.getParameterValuesAsString("coresPerJob");
                         vector<string> queueName = paraObjPartitions.getParameterValuesAsString("queueName");
                         vector<string> queueNamePartTwo = paraObjPartitions.getParameterValuesAsString("queueNamePartTwo");
                         vector<string> maxRuntimePerJob = paraObjPartitions.getParameterValuesAsString("maxRuntimePerJob");
                         vector<string> parallelEnvironment = paraObjPartitions.getParameterValuesAsString("parallelEnvironment");
                         vector<string> gridManager = paraObjPartitions.getParameterValuesAsString("gridManager");
                         vector<string> projectAllocation = paraObjPartitions.getParameterValuesAsString("projectAllocation");
                         string nodesPerJobDefault = "1";      // 1 system is default
                         string coresPerJobDefault = "16";     // 16 cores default
                         string queueNameDefault = "normal";   // normal queue
                         string queueNamePartTwoDefault = "ib_only"; // sge cluster only
                         string maxRuntimePerJobDefault = "02:00:00"; // 2 hours
                         string parallelEnvironmentDefault = "orte2";  // sge cluster only
                         string gridManagerDefault = "slurm";
                         string projectAllocationDefault = "";
                         //////////////////////////////////////////////////////
                         //
                         // Extract the parameters from the cluster
                         // configuration file.
                         //
                         // To run on a different or new cluster, more
                         // parameters can be added here....
                         //
                         // This provides flexibility/ adaptibility to
                         // run on novel clusters.
                         //
                         //////////////////////////////////////////////////////
                         if ( nodesPerJob.size() == 1 &&
                              coresPerJob.size() == 1 &&
                              queueName.size() == 1 &&
                              queueNamePartTwo.size() == 1 &&
                              maxRuntimePerJob.size() == 1 &&
                              parallelEnvironment.size() == 1 &&
                              gridManager.size() == 1 &&
                              projectAllocation.size() == 1 ) {
                             //////////////////////////////////////////////////
                             //  If user has specified parameters
                             //  within the Cluster
                             //////////////////////////////////////////////////
                             nodesPerJobDefault = nodesPerJob.at(0);
                             coresPerJobDefault = coresPerJob.at(0);
                             queueNameDefault = queueName.at(0);
                             queueNamePartTwoDefault = queueNamePartTwo.at(0);
                             maxRuntimePerJobDefault = maxRuntimePerJob.at(0);
                             parallelEnvironmentDefault = parallelEnvironment.at(0);
                             gridManagerDefault = gridManager.at(0);
                             projectAllocationDefault = projectAllocation.at(0);
                         }

                         //////////////////////////////////////////////////////
                         // Variables from the global parameter file
                         // needed to dynamically
                         //////////////////////////////////////////////////////
                         string mcrPathTemp = mcrPath.at(0);
                         string sharedInstallDirTemp = sharedInstallDir.at(0);
                         string sharedOutputDataDirTemp = sharedOutputDataDir.at(0);
                         string mzdaQuantInstallDirTemp = mzdaQuantInstallDir.at(0);
                         //////////////////////////////////////////////////////
                         // SLURM Job Manager
                         //////////////////////////////////////////////////////
                         if( gridManagerDefault == "slurm" ) {
                             outputFileStreamPeptideQuant<<"#!/bin/bash"<<endl;
                             outputFileStreamPeptideQuant<<"#SBATCH -J mzdaPepQuant"<<fileNumberPeptideQuantJobFile<<endl;
                             outputFileStreamPeptideQuant<<"#SBATCH -o mzdaPepQuant"<<fileNumberPeptideQuantJobFile<<".o%j"<<endl;
                             outputFileStreamPeptideQuant<<"#SBATCH -N "<<nodesPerJobDefault<<endl; // Number of nodes with 1 job
                             outputFileStreamPeptideQuant<<"#SBATCH -n "<<coresPerJobDefault<<endl; // Number of cores on a cluster node
                             outputFileStreamPeptideQuant<<"#SBATCH -p "<<queueNameDefault<<endl;  // queue: normal, development
                             outputFileStreamPeptideQuant<<"#SBATCH -t "<<maxRuntimePerJobDefault<<endl; // max runtime: hh
                             outputFileStreamPeptideQuant<<"#SBATCH -A "<<projectAllocationDefault<<endl; // project allocation information required if multiple projects
                             outputFileStreamPeptideQuant<<mzdaQuantInstallDirTemp<<"/run_mzdaPeptideQuant.sh "<<mcrPathTemp<<" "<<sharedOutputDataDirTemp<<"/parameters.txt "<<fileNumberPeptideQuantJobFile<<" TPPXTandem"<<endl;
                         }
                         //////////////////////////////////////////////////////
                         // Grid Engine Job Manager
                         //////////////////////////////////////////////////////
                         if( gridManagerDefault == "sge" ) {
                             outputFileStreamPeptideQuant<<"#!/bin/bash"<<endl;
                             outputFileStreamPeptideQuant<<"#$ -N mzdaPepQuant"<<endl;
                             outputFileStreamPeptideQuant<<"#$ -j y"<<endl;
                             outputFileStreamPeptideQuant<<"#$ -cwd"<<endl;
                             outputFileStreamPeptideQuant<<"#$ -pe "<<parallelEnvironmentDefault<<" "<<coresPerJobDefault<<endl; // Keep 8 cores as the default.
                             outputFileStreamPeptideQuant<<"#$ -q "<<queueNameDefault<<endl;
                             outputFileStreamPeptideQuant<<"#$ -l "<<queueNamePartTwoDefault<<endl;
                             outputFileStreamPeptideQuant<<mzdaQuantInstallDirTemp<<"/run_mzdaPeptideQuant.sh "<<mcrPathTemp<<" "<<sharedOutputDataDirTemp<<"/parameters.txt "<<fileNumberPeptideQuantJobFile<<" TPPXTandem"<<endl;
                             outputFileStreamPeptideQuant<<"exit 0"<<endl;
                         }
                         //////////////////////////////////////////////////////
                         // Add new job managers here.
                         //////////////////////////////////////////////////////
                     }
                     else{
                         // output file is not open, no need to try to close it
                         return -2;
                     }
                 }
                 catch (std::exception& e) {
                     std::cout << e.what() <<std::endl;
                     std::cout <<"mzda: exception while trying to open cluster job script file (MZDASoft PeptideQuant)"<<std::endl;
                     // make sure to close the file stream
                     outputFileStreamPeptideQuant.close();
                 }  // end trycatch block
                 // make sure to close the output file
                 outputFileStreamPeptideQuant.close();

                 } // end loop through total number of fractions


                 //////////////////////////////////////////////////////////////
                 //Protein Quantification: Merge of all fraction level quantification data
                 //mzdaClusterJobProteinQuant.job
                 //////////////////////////////////////////////////////////////
                 string outputClusterJobFilenameProteinQuant=outputClusterJobsFilenameBase.at(0)+"ProteinQuant.job";
                 std::ofstream outputFileStreamProteinQuant;
                 outputFileStreamProteinQuant.open(outputClusterJobFilenameProteinQuant.c_str(), std::ios::out);
                 try{
                     if( outputFileStreamProteinQuant.is_open()) {
                         //////////////////////////////////////////////////////
                         // These HPC cluster specific items
                         //////////////////////////////////////////////////////
                         vector<string> nodesPerJob = paraObjPartitions.getParameterValuesAsString("nodesPerJob");
                         vector<string> coresPerJob = paraObjPartitions.getParameterValuesAsString("coresPerJob");
                         vector<string> queueName = paraObjPartitions.getParameterValuesAsString("queueName");
                         vector<string> queueNamePartTwo = paraObjPartitions.getParameterValuesAsString("queueNamePartTwo");
                         vector<string> maxRuntimePerJob = paraObjPartitions.getParameterValuesAsString("maxRuntimePerJob");
                         vector<string> parallelEnvironment = paraObjPartitions.getParameterValuesAsString("parallelEnvironment");
                         vector<string> gridManager = paraObjPartitions.getParameterValuesAsString("gridManager");
                         vector<string> projectAllocation = paraObjPartitions.getParameterValuesAsString("projectAllocation");
                         string nodesPerJobDefault = "1";      // 1 system is default
                         string coresPerJobDefault = "16";     // 16 cores default
                         string queueNameDefault = "normal";   // normal queue
                         string queueNamePartTwoDefault = "ib_only"; // sge cluster only
                         string maxRuntimePerJobDefault = "02:00:00"; // 2 hours
                         string parallelEnvironmentDefault = "orte2";  // sge cluster only
                         string gridManagerDefault = "slurm";
                         string projectAllocationDefault = "";
                         //////////////////////////////////////////////////////
                         //
                         // Extract the parameters from the cluster
                         // configuration file.
                         //
                         // To run on a different or new cluster, more
                         // parameters can be added here....
                         //
                         // This provides flexibility/ adaptibility to
                         // run on novel clusters.
                         //
                         //////////////////////////////////////////////////////
                         if ( nodesPerJob.size() == 1 &&
                              coresPerJob.size() == 1 &&
                              queueName.size() == 1 &&
                              queueNamePartTwo.size() == 1 &&
                              maxRuntimePerJob.size() == 1 &&
                              parallelEnvironment.size() == 1 &&
                              gridManager.size() == 1 &&
                              projectAllocation.size() == 1 ) {
                             //////////////////////////////////////////////////
                             //  If user has specified parameters
                             //  within the Cluster
                             //////////////////////////////////////////////////
                             nodesPerJobDefault = nodesPerJob.at(0);
                             coresPerJobDefault = coresPerJob.at(0);
                             queueNameDefault = queueName.at(0);
                             queueNamePartTwoDefault = queueNamePartTwo.at(0);
                             maxRuntimePerJobDefault = maxRuntimePerJob.at(0);
                             parallelEnvironmentDefault = parallelEnvironment.at(0);
                             gridManagerDefault = gridManager.at(0);
                             projectAllocationDefault = projectAllocation.at(0);
                         }

                         //////////////////////////////////////////////////////
                         // Variables from the global parameter file
                         // needed to dynamically
                         //////////////////////////////////////////////////////
                         string mcrPathTemp = mcrPath.at(0);
                         string sharedInstallDirTemp = sharedInstallDir.at(0);
                         string sharedOutputDataDirTemp = sharedOutputDataDir.at(0);
                         string mzdaQuantInstallDirTemp = mzdaQuantInstallDir.at(0);
                         //////////////////////////////////////////////////////
                         // SLURM Job Manager
                         //////////////////////////////////////////////////////
                         if( gridManagerDefault == "slurm" ) {
                             outputFileStreamProteinQuant<<"#!/bin/bash"<<endl;
                             outputFileStreamProteinQuant<<"#SBATCH -J mzdaProtQuant"<<endl;
                             outputFileStreamProteinQuant<<"#SBATCH -o mzdaProtQuant"<<".o%j"<<endl;
                             outputFileStreamProteinQuant<<"#SBATCH -N "<<nodesPerJobDefault<<endl; // Number of nodes with 1 job
                             outputFileStreamProteinQuant<<"#SBATCH -n "<<coresPerJobDefault<<endl; // Number of cores on a cluster node
                             outputFileStreamProteinQuant<<"#SBATCH -p "<<queueNameDefault<<endl;  // queue: normal, development
                             outputFileStreamProteinQuant<<"#SBATCH -t "<<maxRuntimePerJobDefault<<endl; // max runtime: hh
                             outputFileStreamProteinQuant<<"#SBATCH -A "<<projectAllocationDefault<<endl; // project allocation information required if multiple projects
                             outputFileStreamProteinQuant<<mzdaQuantInstallDirTemp<<"/run_mzdaProteinQuant.sh "<<mcrPathTemp<<" "<<sharedOutputDataDirTemp<<"/parameters.txt"<<endl;
                         }
                         //////////////////////////////////////////////////////
                         // Grid Engine Job Manager
                         //////////////////////////////////////////////////////
                         if( gridManagerDefault == "sge" ) {
                             outputFileStreamProteinQuant<<"#!/bin/bash"<<endl;
                             outputFileStreamProteinQuant<<"#$ -N mzdaProtQuant"<<endl;
                             outputFileStreamProteinQuant<<"#$ -j y"<<endl;
                             outputFileStreamProteinQuant<<"#$ -cwd"<<endl;
                             outputFileStreamProteinQuant<<"#$ -pe "<<parallelEnvironmentDefault<<" "<<coresPerJobDefault<<endl; // Keep 8 cores as the default.
                             outputFileStreamProteinQuant<<"#$ -q "<<queueNameDefault<<endl;
                             outputFileStreamProteinQuant<<"#$ -l "<<queueNamePartTwoDefault<<endl;
                             outputFileStreamProteinQuant<<mzdaQuantInstallDirTemp<<"/run_mzdaProteinQuant.sh "<<mcrPathTemp<<" "<<sharedOutputDataDirTemp<<"/parameters.txt"<<endl;
                             outputFileStreamProteinQuant<<"exit 0"<<endl;
                         }
                         //////////////////////////////////////////////////////
                         // Add new job managers here.
                         //////////////////////////////////////////////////////
                     }
                     else{
                         // output file is not open, no need to try to close it
                         return -2;
                     }
                 }
                 catch (std::exception& e) {
                     std::cout << e.what() <<std::endl;
                     std::cout <<"mzda: exception while trying to open cluster job script file (MZDASoft PeptideQuant)"<<std::endl;
                     // make sure to close the file stream
                     outputFileStreamProteinQuant.close();
                 }  // end trycatch block
                 // make sure to close the output file
                 outputFileStreamProteinQuant.close();


                 //////////////////////////////////////////////////////////////
                 //MZDAQuant Pipeline Integration
                 //END BLOCK
                 //////////////////////////////////////////////////////////////


                 //////////////////////////////////////////////////////////////
                 //START-DATABASE FILE VALIDATION JOBS
                 //////////////////////////////////////////////////////////////


                 int totalDBFiles = totalInputFiles;  // need to read this variable in from parameter file
                 int totalDBFilesCounter=0;

                 for ( totalDBFilesCounter = 0; totalDBFilesCounter < totalDBFiles; totalDBFilesCounter++ ) {

                 string outputFilenameNameExtension = mzxmlfilenamebaseDefaultMerged.at(totalDBFilesCounter);
                 assert( outputFilenameNameExtension.size() > 0);

                 stringstream ssDBValidateJobFile;
                 ssDBValidateJobFile<<(totalDBFilesCounter+1); // integer to string conversion
                 string fileNumberDBValidateJobFile = ssDBValidateJobFile.str(); // converting stringstream to string object

                 string outputClusterJobFilenameDBValidateJobFile=outputClusterJobsFilenameBase.at(0)+fileNumberDBValidateJobFile+"DBValidate.job";
                 std::ofstream outputFileStreamDBValidateJobFile;
                 outputFileStreamDBValidateJobFile.open(outputClusterJobFilenameDBValidateJobFile.c_str(), std::ios::out);
                 try{
                     if( outputFileStreamDBValidateJobFile.is_open()) {
                         //////////////////////////////////////////////////////
                         // These HPC cluster specific items
                         //////////////////////////////////////////////////////
                         vector<string> nodesPerJob = paraObjPartitions.getParameterValuesAsString("nodesPerJob");
                         vector<string> coresPerJob = paraObjPartitions.getParameterValuesAsString("coresPerJob");
                         vector<string> queueName = paraObjPartitions.getParameterValuesAsString("queueName");
                         vector<string> queueNamePartTwo = paraObjPartitions.getParameterValuesAsString("queueNamePartTwo");
                         vector<string> maxRuntimePerJob = paraObjPartitions.getParameterValuesAsString("maxRuntimePerJob");
                         vector<string> parallelEnvironment = paraObjPartitions.getParameterValuesAsString("parallelEnvironment");
                         vector<string> gridManager = paraObjPartitions.getParameterValuesAsString("gridManager");
                         vector<string> projectAllocation = paraObjPartitions.getParameterValuesAsString("projectAllocation");
                         string nodesPerJobDefault = "1";      // 1 system is default
                         string coresPerJobDefault = "16";     // 16 cores default
                         string queueNameDefault = "normal";   // normal queue
                         string queueNamePartTwoDefault = "ib_only"; // sge cluster only
                         string maxRuntimePerJobDefault = "02:00:00"; // 2 hours
                         string parallelEnvironmentDefault = "orte2";  // sge cluster only
                         string gridManagerDefault = "slurm";
                         string projectAllocationDefault = "";
                         //////////////////////////////////////////////////////
                         //
                         // Extract the parameters from the cluster
                         // configuration file.
                         //
                         // To run on a different or new cluster, more
                         // parameters can be added here....
                         //
                         // This provides flexibility/ adaptibility to
                         // run on novel clusters.
                         //
                         //////////////////////////////////////////////////////
                         if ( nodesPerJob.size() == 1 &&
                              coresPerJob.size() == 1 &&
                              queueName.size() == 1 &&
                              queueNamePartTwo.size() == 1 &&
                              maxRuntimePerJob.size() == 1 &&
                              parallelEnvironment.size() == 1 &&
                              gridManager.size() == 1 &&
                              projectAllocation.size() == 1 ) {
                             //////////////////////////////////////////////////
                             //  If user has specified parameters
                             //  within the Cluster
                             //////////////////////////////////////////////////
                             nodesPerJobDefault = nodesPerJob.at(0);
                             coresPerJobDefault = coresPerJob.at(0);
                             queueNameDefault = queueName.at(0);
                             queueNamePartTwoDefault = queueNamePartTwo.at(0);
                             maxRuntimePerJobDefault = maxRuntimePerJob.at(0);
                             parallelEnvironmentDefault = parallelEnvironment.at(0);
                             gridManagerDefault = gridManager.at(0);
                             projectAllocationDefault = projectAllocation.at(0);
                         }

                         //////////////////////////////////////////////////////
                         // Variables from the global parameter file
                         // needed to dynamically
                         //////////////////////////////////////////////////////
                         string mcrPathTemp = mcrPath.at(0);
                         string sharedInstallDirTemp = sharedInstallDir.at(0);
                         string sharedOutputDataDirTemp = sharedOutputDataDir.at(0);
                         string mzdaQuantInstallDirTemp = mzdaQuantInstallDir.at(0);

                         string outputFilenameNameExtensionNoFileType = outputFilenameNameExtension;
                         size_t location = outputFilenameNameExtensionNoFileType.find(".mzXML");
                         outputFilenameNameExtensionNoFileType.replace(location,std::string(".mzXML").length(),"");
                         //////////////////////////////////////////////////////
                         // SLURM Job Manager
                         //////////////////////////////////////////////////////
                         if( gridManagerDefault == "slurm" ) {
                             outputFileStreamDBValidateJobFile<<"#!/bin/bash"<<endl;
                             outputFileStreamDBValidateJobFile<<"#SBATCH -J mzdaDBValidate"<<fileNumberDBValidateJobFile<<endl;
                             outputFileStreamDBValidateJobFile<<"#SBATCH -o mzdaDBValidate"<<fileNumberDBValidateJobFile<<".o%j"<<endl;
                             outputFileStreamDBValidateJobFile<<"#SBATCH -N "<<nodesPerJobDefault<<endl; // Number of nodes with 1 job
                             outputFileStreamDBValidateJobFile<<"#SBATCH -n "<<coresPerJobDefault<<endl; // Number of cores on a cluster node
                             outputFileStreamDBValidateJobFile<<"#SBATCH -p "<<queueNameDefault<<endl;  // queue: normal, development
                             outputFileStreamDBValidateJobFile<<"#SBATCH -t "<<maxRuntimePerJobDefault<<endl; // max runtime: hh
                             outputFileStreamDBValidateJobFile<<"#SBATCH -A "<<projectAllocationDefault<<endl; // project allocation information required if multiple projects
                             outputFileStreamDBValidateJobFile<<sharedInstallDirTemp<<"/product/validateDatabaseFile.sh "<<sharedInstallDirTemp<<" "<<sharedOutputDataDirTemp<<" mzdaLevelThreeFile"<<fileNumberDBValidateJobFile<<"-"<<outputFilenameNameExtensionNoFileType<<endl;



                         }
                         //////////////////////////////////////////////////////
                         // Grid Engine Job Manager
                         //////////////////////////////////////////////////////
                         if( gridManagerDefault == "sge" ) {
                             outputFileStreamDBValidateJobFile<<"#!/bin/bash"<<endl;
                             outputFileStreamDBValidateJobFile<<"#$ -N mzdaDBValidate"<<endl;
                             outputFileStreamDBValidateJobFile<<"#$ -j y"<<endl;
                             outputFileStreamDBValidateJobFile<<"#$ -cwd"<<endl;
                             outputFileStreamDBValidateJobFile<<"#$ -pe "<<parallelEnvironmentDefault<<" "<<coresPerJobDefault<<endl; // Keep 8 cores as the default.
                             outputFileStreamDBValidateJobFile<<"#$ -q "<<queueNameDefault<<endl;
                             outputFileStreamDBValidateJobFile<<"#$ -l "<<queueNamePartTwoDefault<<endl;
                             outputFileStreamDBValidateJobFile<<sharedInstallDirTemp<<"/product/validateDatabaseFile.sh "<<sharedInstallDirTemp<<" "<<sharedOutputDataDirTemp<<" mzdaLevelThreeFile"<<fileNumberDBValidateJobFile<<"-"<<outputFilenameNameExtensionNoFileType<<endl;


                             outputFileStreamDBValidateJobFile<<"exit 0"<<endl;
                         }
                         //////////////////////////////////////////////////////
                         // Add new job managers here.
                         //////////////////////////////////////////////////////
                     }
                     else{
                         // output file is not open, no need to try to close it
                         return -2;
                     }
                 }
                 catch (std::exception& e) {
                     std::cout << e.what() <<std::endl;
                     std::cout <<"mzda: exception while trying to open cluster job script file (MZDASoft PeptideQuant)"<<std::endl;
                     // make sure to close the file stream
                     outputFileStreamDBValidateJobFile.close();
                 }  // end trycatch block
                 // make sure to close the output file
                 outputFileStreamDBValidateJobFile.close();

                 } // end loop through total number of fractions


                 //////////////////////////////////////////////////////////////
                 //END-DATABASE FILE VALIDATION JOBS
                 //////////////////////////////////////////////////////////////



            ///////////////////////////////////////////////////////////////////
            ///////////////////////////////////////////////////////////////////
            }   // end mode 5
			
		}
		else {
			
			// return a parameter error
			throw "invalid -m parameter values, multiple items for the -m parameter are not allowed.";
		}
		

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////	
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/////////////////////////////End of user code//////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
	} //end try block
	catch (char const* str) {
		//
		// This needs be changed to a derived class from std::exception.
		// We need to design a general error class for library.
		//
		cout<<"Exception produced: "<<str<<endl;
    } // end catch block
	///////////////////////////////////////////////////////////////////////////
	// msda Catch All Exception Handler
	///////////////////////////////////////////////////////////////////////////
	catch( exception &  e){
		cout<<"Exception encountered at Catch-All handler"<<endl;
		cout<<e.what()<<endl;
	}
	///////////////////////////////////////////////////////////////////////////
	// Return successfully
	///////////////////////////////////////////////////////////////////////////

	return 0;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////	
///////////////////////////////////////////////////////////////////////////////	
}  // end main

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

