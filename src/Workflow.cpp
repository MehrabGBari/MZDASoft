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
// Workflow.cpp
// 
// Part of MZDASoft.
// Manages the runtime workflow of the system.
// 
// mzdaWorkflow1.sh
//   - Prepackaged script produced as a top level wrapper for mzda executable.
//     It sets up the LD_LIBRARY_PATH environment variable and allows for
//     running the mzda executable for workflow #1.
//     workflow #1 generates all necessary parameter files to distribute
//     the processing of one or more source data files as well as the
//     the scripts that will be needed to perform subsequent stages of
//     processing.
//
// mzdaWorkflow2.sh
// mzdaWorkflow3.sh
// mzdaWorkflow4.sh
// 
// version 1.0.0
//
// Change History:
// Feature ID, Version, Date, Description, @[ADD,CHG,MOV]DeveloperID
// ft00000, vr1.0.0, 04/23/12-09/25/12 Initial Version, @ADDNRZ
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
//
///////////////////////////////////////////////////////////////////////////////
#include "Workflow.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <exception>
#include <cassert>


///////////////////////////////////////////////////////////////////////////////
//std namespace 
///////////////////////////////////////////////////////////////////////////////
using namespace std;


///////////////////////////////////////////////////////////////////////////////
//Start msda namespace
///////////////////////////////////////////////////////////////////////////////
namespace msda{


///////////////////////////////////////////////////////////////////////////////
//Constructor
///////////////////////////////////////////////////////////////////////////////
Workflow::Workflow(){
;
}  // end Workflow()

///////////////////////////////////////////////////////////////////////////////
//Destructor
///////////////////////////////////////////////////////////////////////////////
Workflow::~Workflow(){
;
} // end ~Workflow()




///////////////////////////////////////////////////////////////////////////////
// Name:  workflowZero
//  The purpose of method is to implement workflow #0.
//  Workflow# 0.  This is the startup workflow.  It uses only the
//  global parameter file and so
//
//
//
// Algorithm:
//
//
// Input:
//
//  1. ParameterLoader object reference, containing global parameters
//
// Output:
//  1. [ a set of output files are generated, according to global parameters ]
//
// Return code: 0 successful,
//
/// Nelson.Ramirez@utsa.edu, 05/31/12
///////////////////////////////////////////////////////////////////////////////
int Workflow::workflowZero( cbi::ParameterLoader & paraObjGlobal  ) {

    ///////////////////////////////////////////////////////////////////////////
    //
    // workflow 0
    //
    ///////////////////////////////////////////////////////////////////////////
    int rc = 0; // used to keep track of return codes

    ///////////////////////////////////////////////////////////////////////////
    // Get the runtime parameters.
    // This is information that is used from the global parameter file to
    // generate localized information for each worker.
    // We need to minimize the number of code changes, when a parameter
    // name needs to be updated, therefore, the workflow# methods
    // should be the only place where changes should be needed when
    // changing the name / value of a parameter.
    ///////////////////////////////////////////////////////////////////////////

    // Get the working directory
    vector<string> workingdirectory = paraObjGlobal.getParameterValuesAsString("workingdirectory");
    cout<<"workingdirectory.size()="<<workingdirectory.size()<<endl;
    cout<<"workingdirectory.at(0)="<<workingdirectory.at(0)<<endl;

    // The install directory is the location of the mzda executable
    // The run scripts will be tailored to the working directory environment
    // and so the run scripts will be placed
    vector<string> installdirectory = paraObjGlobal.getParameterValuesAsString("installdirectory");
    cout<<"installdirectory.size()="<<installdirectory.size()<<endl;
    cout<<"installdirectory.at(0)="<<installdirectory.at(0)<<endl;

    // The cbilib install directory is the location os the
    vector<string> cbilibinstalldirectory = paraObjGlobal.getParameterValuesAsString("cbilibinstalldirectory");
    cout<<"cbilibinstalldirectory.size()="<<cbilibinstalldirectory.size()<<endl;
    cout<<"cbilibinstalldirectory.at(0)="<<cbilibinstalldirectory.at(0)<<endl;


    cout<<"workflowZero, generating the script file for workflow one"<<endl;


    // todo: extend to create a new script for each input file
    rc = generateWorkflowOneScript(workingdirectory.at(0), installdirectory.at(0), cbilibinstalldirectory.at(0) );  // this needs to be extended to create a script for each input file

    return 0;
}  // end workflowZero





///////////////////////////////////////////////////////////////////////////////
// Name:  workflowOne
//  The purpose of method is to implement workflow #1.
//  Workflow# 1 contains all the code necessary to generate a set of 
//  scripts that will implement the functionality end-to-end, either
//  within a single multi-core system via shell scripts or in a grid
//  environment with grid job scripts.  This is the only mode
//  that a user should need to call directly.  
//
//  All subsequent processing will be done via this interface layer.
//  This will provide the flexibility needed to support a variety of 
//  run modes in the future.  
//
//
//
//
// Algorithm:
//
//   
// Input:
//
//  1. ParameterLoader object reference, containing global parameters
// 
// Output: 
//  1. [ a set of output files are generated, according to global parameters ]
//
// Return code: 0 successful, 
// 
/// Nelson.Ramirez@utsa.edu, 04/16/12-09/25/12
///////////////////////////////////////////////////////////////////////////////
int Workflow::workflowOne( cbi::ParameterLoader & paraObjGlobal  ) {

	///////////////////////////////////////////////////////////////////////////
	//
	// workflow 1
	//
	///////////////////////////////////////////////////////////////////////////
	int rc = 0; // used to keep track of return codes

    // Get the working directory
    vector<string> workingdirectory = paraObjGlobal.getParameterValuesAsString("workingdirectory");
    cout<<"workingdirectory.size()="<<workingdirectory.size()<<endl;
    cout<<"workingdirectory.at(0)="<<workingdirectory.at(0)<<endl;

    // The install directory is the location of the mzda executable
    // The run scripts will be tailored to the working directory environment
    // and so the run scripts will be placed
    vector<string> installdirectory = paraObjGlobal.getParameterValuesAsString("installdirectory");
    cout<<"installdirectory.size()="<<installdirectory.size()<<endl;
    cout<<"installdirectory.at(0)="<<installdirectory.at(0)<<endl;

    // The cbilib install directory is the location where the CBILIB is installed.
    vector<string> cbilibinstalldirectory = paraObjGlobal.getParameterValuesAsString("cbilibinstalldirectory");
    cout<<"cbilibinstalldirectory.size()="<<cbilibinstalldirectory.size()<<endl;
    cout<<"cbilibinstalldirectory.at(0)="<<cbilibinstalldirectory.at(0)<<endl;

    // The prllpath is the location of the installation of the PRLL system
    vector<string> prllpath= paraObjGlobal.getParameterValuesAsString("prllpath");
    cout<<"prllpath.size()="<<prllpath.size()<<endl;
    cout<<"prllpath.at(0)="<<prllpath.at(0)<<endl;

    // The prllsourcepath is the location of the installation of the PRLL system shell script
    vector<string> prllsourcepath= paraObjGlobal.getParameterValuesAsString("prllsourcepath");
    cout<<"prllsourcepath.size()="<<prllsourcepath.size()<<endl;
    cout<<"prllsourcepath.at(0)="<<prllsourcepath.at(0)<<endl;


    // mzxmlfilename contains one or more input filenames in mzXML format.
    vector<string> mzxmlfilename = paraObjGlobal.getParameterValuesAsString("mzxmlfilename");
    cout<<"mzxmlfilename.size()="<<mzxmlfilename.size()<<endl;
    for ( int i = 0; i < mzxmlfilename.size(); i++ ) {
        cout<<"mzxmlfilename.at("<<i<<")="<<mzxmlfilename.at(i)<<endl;

        // create the set of local parameters for each input file
        stringstream ss;
        ss<<(i+1); // integer to string conversion
        string folderNumber = ss.str(); // converting stringstream to string object
        string parmFilenameString = workingdirectory.at(0)+folderNumber+ "/"+ "localparameters";
        // separate the extension, since each parameter file will be numbered
        string parmFilenameStringExtension = ".txt";

        cout<<"DEBUG 1"<<endl;
        vector<double> parglobalminscan = paraObjGlobal.getParameterValuesAsDouble("parglobalminscan");
        vector<double> parglobalmaxscan = paraObjGlobal.getParameterValuesAsDouble("parglobalmaxscan");
        vector<double> parglobalminmz = paraObjGlobal.getParameterValuesAsDouble("parglobalminmz");
        vector<double> parglobalmaxmz = paraObjGlobal.getParameterValuesAsDouble("parglobalmaxmz");
        vector<double> parmzdelta = paraObjGlobal.getParameterValuesAsDouble("parmzdelta");
        vector<double> parmzdeltaoverlap = paraObjGlobal.getParameterValuesAsDouble("parmzdeltaoverlap");
        vector<double> mzbinningthresholdppm = paraObjGlobal.getParameterValuesAsDouble("mzbinningthresholdppm");

        vector<int> corespermachine = paraObjGlobal.getParameterValuesAsInt("corespermachine");

        // we could later on extend this to have different parameters for each input file
        double minScan = parglobalminscan.at(i);
        double maxScan = parglobalmaxscan.at(i);
        double startMZValue = parglobalminmz.at(i);
        double endMZValue = parglobalmaxmz.at(i);
        double mzDelta = parmzdelta.at(i);
        double mzOverlap = parmzdeltaoverlap.at(i);
        double mzBinningThreshold = mzbinningthresholdppm.at(i);

        cout<<"DEBUG 2, i="<<i<<endl;

        // number of cores
        int numCores = corespermachine.at(0);
        ///////////////////////////////////////////////////////////////////////////
        //
        // Create the job scripts that will manage the processing of the
        // local parameters files
        //
        ///////////////////////////////////////////////////////////////////////////
        rc = generateParameterFiles(	parmFilenameString,
                                        parmFilenameStringExtension,
                                        i+1,
                                        minScan,
                                        maxScan,
                                        startMZValue,
                                        endMZValue,
                                        mzDelta,
                                        mzOverlap,
                                        mzBinningThreshold);
        if ( rc != 0 ) {
            throw "Error in workflow one, generateParameterFiles method";
        }

         cout<<"DEBUG 3, i="<<i<<endl;
        ///////////////////////////////////////////////////////////////////////
        //
        // Create the level 2 scripts that will run the binning algorithms
        //
        //
        ///////////////////////////////////////////////////////////////////////
        cout<<"generateWorkflowTwoScript for file#"<<folderNumber<<endl;
        rc = generateWorkflowTwoScript( workingdirectory.at(0),
                                        installdirectory.at(0),
                                        cbilibinstalldirectory.at(0),
                                        prllpath.at(0),
                                        prllsourcepath.at(0),
                                        i+1,
                                        numCores);
        if ( rc != 0 ) {
            throw "Error in workflow two, generateWorkflowTwoScript method";
        }

         cout<<"DEBUG 4, i="<<i<<endl;
        ///////////////////////////////////////////////////////////////////////
        //
        // Create the level 3 scripts that will run the binning algorithms
        //
        //
        ///////////////////////////////////////////////////////////////////////
        cout<<"generateWorkflowThreeScript for file#"<<folderNumber<<endl;
        rc = generateWorkflowThreeScript(workingdirectory.at(0),
                                         installdirectory.at(0),
                                         cbilibinstalldirectory.at(0),
                                         prllpath.at(0),
                                         prllsourcepath.at(0),
                                         i+1,
                                         numCores);

        if ( rc != 0 ) {
            throw "Error in workflow three, generateWorkflowThreeScript method";
        }

         cout<<"DEBUG 5, i="<<i<<endl;

        ///////////////////////////////////////////////////////////////////////
        //
        // Here we create a cluster workflow for hpc clusters and/or
        // compute cloud infrastructure.
        //
        ///////////////////////////////////////////////////////////////////////
//        cout<<"generateWorkflowClusterScript for file#"<<folderNumber<<endl;
//        rc = generateWorkflowClusterScript(workingdirectory.at(0),
//                                         installdirectory.at(0),
//                                         cbilibinstalldirectory.at(0),
//                                         prllpath.at(0),
//                                         prllsourcepath.at(0),
//                                         i+1,
//                                         numCores);

//        if ( rc != 0 ) {
//            throw "Error in workflow three, generateWorkflowClusterScript method";
//        }

//         cout<<"DEBUG 6, i="<<i<<endl;


	///////////////////////////////////////////////////////////////////////////
	//
	// Create the job scripts that will manage the processing of the
	// local parameters files
    //
    // For each input file to be processed, create a directory to store
    // its level 1, level 2 results.
    //
    // The following directory structure is to be followed:
    //    workingdirectory###
    //      ./1
    //          ./1/level1
    //          ./1/level2
    //          ./1/localparameters#####.txt
    //              - each stage may append more information to the local
    //                parameter files
    //
    //      ./2
    //          ./2/level1
    //          ./2/level2
    //          ./2/localparameters#####.txt
    //
    //      ./n
    //          ./n/level1
    //          ./n/level2
    //          ./n/localparameters#####.txt
    //
    //      ./level3
    //
    //
    //
    //      mzdaWorkflow1.sh
    //      mzdaWorkflow2_file#.sh, file numbers map to directories 1,2,...,n
    //      mzda
    //      parameters.txt   ( global parameter file )
    //
    //      Note: input mzXML files can reside anywhere
    //
    //      This organizational structure will allow separating
    //      the processing of any number of input files to generate
    //      all the level 1 format data within each input files
    //      corresponding directory.
    //
    //
    //
	// 
	///////////////////////////////////////////////////////////////////////////


    ///////////////////////////////////////////////////////////////////////
    // Workflow 4 script:
    // Here we aggregate all the MZDALevel 3 items ( LCC's ) contained
    // within a section's directory
    //
    ///////////////////////////////////////////////////////////////////////
         cout<<"generateWorkflowFourScript for file#"<<folderNumber<<endl;

         rc = generateWorkflowFourScript(workingdirectory.at(0),
                                         installdirectory.at(0),
                                         cbilibinstalldirectory.at(0),
                                         prllpath.at(0),
                                         prllsourcepath.at(0),
                                         i+1,
                                         numCores,
                                         true,
                                         true);

         if ( rc != 0 ) {
             throw "Error in workflow four, generateWorkflowFourScript method";
         }


    }





    return 0;
}  // end workflowOne





///////////////////////////////////////////////////////////////////////////////
// Name:  workflowTwo
//  The purpose of method is to implement workflow #2.
//  Workflow# 2 reads in a section of an mzXML file and proceeds to perform
//  the first level binning algorithm.  This is a critical first stage
//  of processing, generating an organized set of data within a certain
//  retention time(scan) and mz range.
//
//  This workflow is handled by a multitude of separate operating system
//  processes, each reading in a separate portion of the input dataset.
//
//  Each process generates an output file in MZDA Level 1 format.
//  MZDA Level 1 format is a text based format containing the
//  organized regions of interest within the source dataset.
//
//  MZDA Level 1 format data can subsequently be used for a variety
//  of further analysis.
//
// Algorithm:
//  Binning approach
//
// Input:
//
//  1. ParameterLoader object reference, containing global parameters
//  2. ParameterLoader object reference, containing local parameters
//  3. ScanGroup object reference, a container for the region of interest data subset
//  4. MsDataLoader object reference, an object in charge of loading the mzXML data subset into a ScanGroup container
//  5. SignalProcessor object reference, an object in charge of a number of "Matlab -like" signal processing functions
//
// Output:
//  1. An MZDA Level 1 format .csv text file containing only that data which
//     the worker running this workflow is responsible for.
//     As a separate workflow, we could also merge the binning workflow with the
//     lc candidate generation workflow for performance efficiency.
//     However, both options should exist, that of going through the
//     generation of an MZDA Level 1 format csv file, as well as
//     going directly to the set of LC candidates.
//
//
//
// Return code: 0 successful,
//
/// Nelson.Ramirez@utsa.edu, 05/30/12
///////////////////////////////////////////////////////////////////////////////
int Workflow::workflowTwo( cbi::ParameterLoader & paraObjGlobal ,
                           cbi::ParameterLoader & paraObjLocal,
                           msda::ScanGroup & sgObject ,
                           msda::MsDataLoader & msdlObject ,
                           cbi::SignalProcessor & spObject) {

    ///////////////////////////////////////////////////////////////////////////
    //
    // workflow 2
    //
    ///////////////////////////////////////////////////////////////////////////

    int rc = 0; // used to keep track of return codes

    ///////////////////////////////////////////////////////////////////////////
    // Finding the folder number
    // We need to first find the filenum identifier to know the directory
    // this process is working on. Note how each directory corresponds
    // to a single input file.
    ///////////////////////////////////////////////////////////////////////////
    vector<int> filenum = paraObjLocal.getParameterValuesAsInt("filenum");
    vector<string> filenumStr = paraObjLocal.getParameterValuesAsString("filenum");
    cout<<"filenum.size()="<<filenum.size()<<endl;
    cout<<"filenum.at(0)="<<filenum.at(0)<<endl;
    int folderNumber = filenum.at(0);
    int folderNumberIndex = folderNumber -1;


    ///////////////////////////////////////////////////////////////////////////
    // Get the level of the ms scans to process.
    // MSDA handles level 1 scans only.
    // Level 1 scan are
    ///////////////////////////////////////////////////////////////////////////
    vector<int> mslevel = paraObjGlobal.getParameterValuesAsInt("mslevelselector");
    cout<<"mslevel="<<mslevel.at(folderNumberIndex)<<endl;

    ///////////////////////////////////////////////////////////////////////////
    // Load the local parameter file data
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    // Get the mzXMLFilename(s)
    // For a future release, we will want to be able to process
    // multiple files in this same type of distributed manner.
    // For example, a use case can be as follows:
    // For the following 60 mzXML files, get data within the
    // following retention time & mz region and perform some analysis
    // on this data...
    // For this purpose, the parameter object framework provides
    // for a very powerful and flexible mechanism for distributed processing.
    ///////////////////////////////////////////////////////////////////////////
    vector<string> mzxmlfilename = paraObjGlobal.getParameterValuesAsString("mzxmlfilename");
    string inputfilename = mzxmlfilename.at(folderNumberIndex);  // vectors start at 0
    cout<<"inputfilename="<<inputfilename<<endl;

    ///////////////////////////////////////////////////////////////////////////
    // Get the list of scan numbers that this executable instance
    // is responsible for.  This expands the {} syntax and allows for
    // straightforward specification of multiple sets of scan ranges
    // within the local parameter file.
    ///////////////////////////////////////////////////////////////////////////
    vector<int> scanList = paraObjLocal.getParameterValuesAsIntRange("localworkerscansrange");

    ///////////////////////////////////////////////////////////////////////////
    // Get the identifier of the scan group this executable is working
    // on.
    ///////////////////////////////////////////////////////////////////////////
    vector<int> scanGroupId = paraObjLocal.getParameterValuesAsInt("scangroupid");
    int sgId = scanGroupId.at(0);

    ///////////////////////////////////////////////////////////////////////////
    // Get the local worker minimum MZ value we will be processing
    ///////////////////////////////////////////////////////////////////////////
    vector<double> minMzLocal = paraObjLocal.getParameterValuesAsDouble("localminmz");
    double minMz = minMzLocal.at(0);

    ///////////////////////////////////////////////////////////////////////////
    // Get the local worker maximum MZ value we will be processing
    ///////////////////////////////////////////////////////////////////////////
    vector<double> maxMzLocal = paraObjLocal.getParameterValuesAsDouble("localmaxmz");
    double maxMz = maxMzLocal.at(0);

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    /////////////////KEY CONTROL PARAMETER: ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    // Get the mzBinning threshold in PPM.
    // This threshold is used to cut the mz bins
    // This control value is the first critical user controlled
    // parameter.
    ///////////////////////////////////////////////////////////////////////////
    vector<double> mzBinningThresholdInPPMLocal = paraObjLocal.getParameterValuesAsDouble("mzbinningthresholdinppm");
    double mzBinningThresholdInPPM = mzBinningThresholdInPPMLocal.at(0);

    ///////////////////////////////////////////////////////////////////////////
    // Handle a scan list parameter error
    // Throw an exception.  ( Later change to a derived exception object )
    // There should be at least one scan to process, and the
    // scan number must be a positive integer as required by the
    // getParameterValuesAsIntRange api specification.
    ///////////////////////////////////////////////////////////////////////////
    if ( scanList.size() == 1 ) {
        if ( scanList.at(0) < 0 ) {
            throw "Invalid scan list range exception";
        }
    }


    ///////////////////////////////////////////////////////////////////////////
    //
    // Load the raw LC/MS data from an mzXML file.
    // This data must be raw data, not peak picked data, although
    // later on we could add the option to process already peak picked
    // data in a format such as Mascot's MGF file format.
    //
    ///////////////////////////////////////////////////////////////////////////
    rc = 0;// clear return code
    rc = msdlObject.loadMzXmlFile(inputfilename,
                                  sgObject,
                                  scanList,
                                  mslevel.at(folderNumberIndex),
                                  minMz,
                                  maxMz);
    if ( rc != 0 ){
        throw "Error during call to msdlObject.loadMzXmlFile";
    }

    ///////////////////////////////////////////////////////////////////////////
    // Output Scan Group Information to Standard Output
    ///////////////////////////////////////////////////////////////////////////
    int minScan=0; // always initialize local variables( stack variables )
    int maxScan=0; // always initialize local variables ( stack variables )
    double minRetentionTime=0;
    double maxRetentionTime=0;
    spObject.findMinValue(scanList,minScan); // Find minimum scan number
    spObject.findMaxValue(scanList,maxScan); // Find maximum scan number
    minRetentionTime = sgObject.findMinRetentionTime(); // Find minimum retention time
    maxRetentionTime = sgObject.findMaxRetentionTime(); // Find maximum retention time

    ///////////////////////////////////////////////////////////////////////////
    // Generate an MZDA Level 1 format file.
    // This represents the first level of processing
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // Find the pre LC candidate regions of interest
    // These represent the starting point for all further processing
    //
    // Interesting effect noticed:
    //
    // As the number of scans increases, for the same mzBinning
    // threshold, the number of bins are reduced.
    //
    //
    ///////////////////////////////////////////////////////////////////////////
    rc = 0; // clear return code
    rc = sgObject.findPreLCRegionsOfInterest(mzBinningThresholdInPPM);

    if ( rc != 0 ){
        if ( rc == -1 ){
            throw "Not enough data matches mz and retention time criteria to proceed";
        }
        else {
            throw "Error during call to sgObject.findPreLcCandidates";
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // MZDA Level 1 File Output
    // This format summarizes the data within an mzXML file into a
    // region in mz and retention time, sorted into bins.
    //
    // The convention used for the MZDA Level 1 format files
    // will be mzdaOne######.csv
    //
    // This text format has up to 4.5x compression ratio when compressed
    // as a tar.bz2 format file.
    ///////////////////////////////////////////////////////////////////////////
    rc = 0; // clear return code


    // Get the working directory
    vector<string> workingdirectory = paraObjGlobal.getParameterValuesAsString("workingdirectory");
    cout<<"workingdirectory.size()="<<workingdirectory.size()<<endl;
    cout<<"workingdirectory.at(0)="<<workingdirectory.at(0)<<endl;

    // construct the output filename pre-fix for level 1 output files
    string levelOneFilename1 = workingdirectory.at(0);
    string levelOneFilename2 =  filenumStr.at(0)+"/level1/mzdaOne";
    string outputFilenameString = levelOneFilename1 + levelOneFilename2;
    cout<<"outputFilenameString="<<outputFilenameString<<endl;

    rc = sgObject.createMSDALevel1File("1.0.0",
                                       outputFilenameString,
                                       inputfilename,
                                       sgId,
                                       6,
                                       10,
                                       minScan,
                                       maxScan,
                                       minMz,
                                       maxMz,
                                       minRetentionTime,
                                       maxRetentionTime,
                                       mzBinningThresholdInPPM,
                                       mslevel.at(folderNumberIndex),
                                       sgObject.size() );
    if ( rc != 0 ){
        throw "Error during call to sgObject.createMSDALevel1File";
    }

    ///////////////////////////////////////////////////////////////////////////
    // At this point we should have an MZDA Level 1 file created.
    // All subsequent processing stages can proceed from these set of files.
    ///////////////////////////////////////////////////////////////////////////

    return 0;
}  // end workflowTwo



///////////////////////////////////////////////////////////////////////////////
// Name:  workflowThree
//  The purpose of method is to implement workflow #3.
//  Workflow# 3 reads in an Level 1 file and proceeds to perform
//  the LC Candidate generation.
//
//  Since the LC Candidate generation phase is a completely
//  independent process, we generate a set of LC Candidates within
//  an Level 2 file.  Each row of this file will contain
//  the data for a single LC Candidate.
//
//  The output of each worker processing an Level 1 file is an
//  Level 2 file containing a set of LC Candidates.
//
//  The format of an Level 2 file is as follows:
//  header row containing information to allow us to get back to the raw data
//    - A mapping from scan # to retention time would be very useful in the header of each file.
//
//
//  lccandidateid, lccandidate data
//  0, elution profile,maxelutionpeakscannumber,mzpeak(mz,intensity,mz,intensity,mz,intensity),mzstart,mzend,scanstart,scanend,rtstart,rtend, resampledregionid,centermz
//  1, elution profile,maxelutionpeakscannumber,mzpeak(mz,intensity,mz,intensity,mz,intensity),mzstart,mzend,scanstart,scanend,rtstart,rtend, resampledregionid,centermz
//
//  elution profile: contains the collapsed view of the <mz,int> data items
//  along the retention time dimension.
//
//  maxelutionpeakscannumber:  this is the scan number that corresponds to the
//   peak location of the elution profile.
//
//  mzpeak:  this is the data of the mz peak as mz, intensity value pairs
//
//
// Algorithm:
//  - A "Box" organization approach
//
// Input:
//
//  1. ParameterLoader object reference, containing global parameters
//  2. ParameterLoader object reference, containing local parameters
//  3. ScanGroup object reference, a container for the region of interest data subset
//  4. MsDataLoader object reference, an object in charge of loading the mzXML data subset into a ScanGroup container
//  5. SignalProcessor object reference, an object in charge of a number of "Matlab -like" signal processing functions
//
//
//
// Output:
//  1. A Level 2 format .csv text file containing the
//     set of LC Candidates generated from a given work unit.
//
// Return code: 0 successful,
//
/// Nelson.Ramirez@utsa.edu, 07/09/12
///////////////////////////////////////////////////////////////////////////////
int Workflow::workflowThree( cbi::ParameterLoader & paraObjGlobal ,
                           cbi::ParameterLoader & paraObjLocal,
                           msda::ScanGroup & sgObject ,
                           msda::MsDataLoader & msdlObject ,
                           cbi::SignalProcessor & spObject ) {


    ///////////////////////////////////////////////////////////////////////////
    //
    // workflow 3: LC Candidate Generation
    //
    ///////////////////////////////////////////////////////////////////////////

    int rc = 0; // used to keep track of return codes


    ///////////////////////////////////////////////////////////////////////////
    // Finding the folder number
    // We need to first find the filenum identifier to know the directory
    // this process is working on. Note how each directory corresponds
    // to a single input file.
    ///////////////////////////////////////////////////////////////////////////
    vector<int> filenum = paraObjLocal.getParameterValuesAsInt("filenum");
    vector<string> filenumStr = paraObjLocal.getParameterValuesAsString("filenum");
    cout<<"filenum.size()="<<filenum.size()<<endl;
    cout<<"filenum.at(0)="<<filenum.at(0)<<endl;
    int folderNumber = filenum.at(0);
    int folderNumberIndex = folderNumber -1;

    ///////////////////////////////////////////////////////////////////////////
    // Get the level of the ms scans to process.
    // MSDA handles level 1 scans only.
    // Level 1 scan are
    ///////////////////////////////////////////////////////////////////////////
    vector<int> mslevel = paraObjGlobal.getParameterValuesAsInt("mslevelselector");


    ///////////////////////////////////////////////////////////////////////////
    // Load the local parameter file data
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    // Get the mzXMLFilename(s)
    // For a future release, we will want to be able to process
    // multiple files in this same type of distributed manner.
    // For example, a use case can be as follows:
    // For the following 60 mzXML files, get data within the
    // following retention time & mz region and perform some analysis
    // on this data...
    // For this purpose, the parameter object framework provides
    // for a very powerful and flexible mechanism for distributed processing.
    ///////////////////////////////////////////////////////////////////////////
    vector<string> mzxmlfilename = paraObjGlobal.getParameterValuesAsString("mzxmlfilename");
    string inputfilename = mzxmlfilename.at(folderNumberIndex);  // vectors start at 0


    ///////////////////////////////////////////////////////////////////////////
    // Get the list of scan numbers that this executable instance
    // is responsible for.  This expands the {} syntax and allows for
    // straightforward specification of multiple sets of scan ranges
    // within the local parameter file.
    ///////////////////////////////////////////////////////////////////////////
    vector<int> scanList = paraObjLocal.getParameterValuesAsIntRange("localworkerscansrange");

    ///////////////////////////////////////////////////////////////////////////
    // Get the identifier of the scan group this executable is working
    // on.
    ///////////////////////////////////////////////////////////////////////////
    vector<int> scanGroupId = paraObjLocal.getParameterValuesAsInt("scangroupid");
    int sgId = scanGroupId.at(0);

    ///////////////////////////////////////////////////////////////////////////
    // Get the local worker minimum MZ value we will be processing
    ///////////////////////////////////////////////////////////////////////////
    vector<double> minMzLocal = paraObjLocal.getParameterValuesAsDouble("localminmz");
    double minMz = minMzLocal.at(0);

    ///////////////////////////////////////////////////////////////////////////
    // Get the local worker maximum MZ value we will be processing
    ///////////////////////////////////////////////////////////////////////////
    vector<double> maxMzLocal = paraObjLocal.getParameterValuesAsDouble("localmaxmz");
    double maxMz = maxMzLocal.at(0);


    ///////////////////////////////////////////////////////////////////////////
    // Load the necessary parameters for LC Candidate generation from the
    // global parameter file.
    ///////////////////////////////////////////////////////////////////////////

    vector<int> numSmoothPoint = paraObjGlobal.getParameterValuesAsInt("numsmoothpoint");
    int numSmoothPointFolderIndexed = numSmoothPoint.at(folderNumberIndex);

    vector<int> minLCLength = paraObjGlobal.getParameterValuesAsInt("minlclength");
    int minLCLengthFolderIndexed = minLCLength.at(folderNumberIndex);

    vector<double> noiseThresholdLevel = paraObjGlobal.getParameterValuesAsDouble("noisethresholdlevel");
    double noiseThresholdLevelFolderIndexed =  noiseThresholdLevel.at(folderNumberIndex);

    vector<double> massResolution = paraObjGlobal.getParameterValuesAsDouble("massresolution");
    double massResolutionFolderIndexed = massResolution.at(folderNumberIndex);

    vector<int> minNumPreLcXicSignals = paraObjGlobal.getParameterValuesAsInt("minnumprelcxicsignals");
    int minNumPreLcXicSignalsIndexed = minNumPreLcXicSignals.at(folderNumberIndex);


    // new parameters added on 02/14/13
    // default minMZLength = 3
    vector<int> minMZLength = paraObjGlobal.getParameterValuesAsInt("minmzlength");
    int minMZLengthFolderIndexed = minMZLength.at(folderNumberIndex);

    // default lcPeakApexTolerance = 1
    vector<int> lcPeakApexTolerance = paraObjGlobal.getParameterValuesAsInt("lcpeakapextolerance");
    int lcPeakApexToleranceFolderIndexed = lcPeakApexTolerance.at(folderNumberIndex);

    // default mzCorrThreshold = 0.8
    vector<double> mzCorrThreshold = paraObjGlobal.getParameterValuesAsDouble("mzcorrthreshold");
    double mzCorrThresholdFolderIndexed = mzCorrThreshold.at(folderNumberIndex);


    ///////////////////////////////////////////////////////////////////////////
    //
    //
    // Load the Level 1 file containing the set of LC Regions of Interest,
    // run the lc candidate generation algorithm, and
    //
    //
    //
    ///////////////////////////////////////////////////////////////////////////

    // Get the working directory
    vector<string> workingdirectory = paraObjGlobal.getParameterValuesAsString("workingdirectory");


    // construct the filename pre-fix for level 1 input files
    string levelOneFilename1 = workingdirectory.at(0);
    string levelOneFilename2 =  filenumStr.at(0)+"/level1/mzdaOne";
    string inputFilenameString= levelOneFilename1 + levelOneFilename2;


    // construct the filename pre-fix for level 2 output files
    string levelTwoFilename1 = workingdirectory.at(0);
    string levelTwoFilename2 =  filenumStr.at(0)+"/level2/mzdaTwo";
    string outputFilenameString= levelTwoFilename1 + levelTwoFilename2;


    ///////////////////////////////////////////////////////////////////////////
    //
    // Generate the input & output filenames
    // according to the corresponding scangroup id and width as specified
    // by the user.  Multiple scangroups could be processed by a single
    // worker, so the fundamental identifier of a unit of work is the
    // scan group id.
    //
    ///////////////////////////////////////////////////////////////////////////

    string valString;
    stringstream scangroupidTemp (stringstream::in | stringstream::out);
    scangroupidTemp << setw(6) << setfill('0');  // ??? This needs to be made more general
    scangroupidTemp <<  sgId;
    scangroupidTemp >> valString;

    ///////////////////////////////////////////////////////////////////////////
    // Create the Level 1 input file string
    ///////////////////////////////////////////////////////////////////////////
    string levelOnePathAndFilename = inputFilenameString;
    levelOnePathAndFilename = levelOnePathAndFilename.append(valString);
    levelOnePathAndFilename = levelOnePathAndFilename.append(".csv");

    ///////////////////////////////////////////////////////////////////////////
    // Create the Level 2 output file string
    ///////////////////////////////////////////////////////////////////////////
    string levelTwoPathAndFilename = outputFilenameString;

    ///////////////////////////////////////////////////////////////////////////
    //
    // Load a level 1 file into memory to a scan group object.
    //
    //
    // This loads the data into a block of memory in as a blob of data
    // for I/O efficiency.
    //
    //
    ///////////////////////////////////////////////////////////////////////////
    rc = msdlObject.loadLevel1File("1.0.0","text","all",0,0,levelOnePathAndFilename,sgObject);
    if ( rc != 0 ){
     throw "Error during call to msdlObject.loadLevel1File";
    }


    ///////////////////////////////////////////////////////////////////////////
    //
    // Apply the LC Candidate generation algorithm
    // to the data within a level 1 format file.
    // Take the in memory buffer initialized by the loadLevel1File
    // and generate internal data structures to be used for LC
    // Candidate generation.
    //
    ///////////////////////////////////////////////////////////////////////////
    rc = sgObject.processLevel1Data();
    if ( rc != 0 ) {
        throw "Error during call the sgObject.processLevel1Data() within Workflow.cpp";
    }


    ///////////////////////////////////////////////////////////////////////////
    // Apply the LC Candidate Generation Algorithm
    // on each PLCROI.
    // Add a loop here after debugging initial algorithm
    //
    // Parallelism note: At this level each call to lccAlgorithm is
    // completely independent, so it can be parallelized across
    // as many cores as are available.
    //
    ///////////////////////////////////////////////////////////////////////////

    // Loop through all PLCROI's: Pre-LC Regions of Interest
    vector<int> numSignalsPerPlcRoi;
    rc = sgObject.getNumSignalsPerPlcRoi(numSignalsPerPlcRoi);

    cout<<"numSignalsPerPlcRoi.size()="<<numSignalsPerPlcRoi.size()<<endl;

    for ( int plcRoiId = 0; plcRoiId < numSignalsPerPlcRoi.size(); plcRoiId++ ){

        // this check is needed to handle PLCROI's that
        // have too few data items.  We need to look into this
        // a bit further.
        if ( numSignalsPerPlcRoi.at(plcRoiId) > minNumPreLcXicSignalsIndexed ){

            ///////////////////////////////////////////////////////////////////
            // Prior version of the algorithm:
            // Note:  Adding new algorithms is simply a matter of
            // creating a new function as shown below.
            // In the future we can make this process more plug-in enabled,
            //
            ///////////////////////////////////////////////////////////////////
            //        rc = sgObject.lccAlgorithm(levelOnePathAndFilename,
            //                                   sgId,
            //                                   plcRoiId,
            //                                   mzWindowPPMFolderIndexed,
            //                                   numSmoothPointFolderIndexed,
            //                                   minLCLengthFolderIndexed,
            //                                   minMZLengthFolderIndexed,
            //                                   noiseThresholdLevelFolderIndexed,
            //                                   massResolutionFolderIndexed,
            //                                   lcPeakApexToleranceFolderIndexed,
            //                                   mzCorrThresholdFolderIndexed,
            //                                   levelTwoPathAndFilename);
            ///////////////////////////////////////////////////////////////////

            ///////////////////////////////////////////////////////////////////
            // Update to new algorithm, 03/05/13
            // Update on 0/27/13, loglevel updates
            ///////////////////////////////////////////////////////////////////
            int logLevel = 1;
            rc = sgObject.lccAlgorithmVersion2(levelOnePathAndFilename,
                                       sgId,
                                       plcRoiId,
                                       numSmoothPointFolderIndexed,
                                       minLCLengthFolderIndexed,
                                       minMZLengthFolderIndexed,
                                       noiseThresholdLevelFolderIndexed,
                                       massResolutionFolderIndexed,
                                       lcPeakApexToleranceFolderIndexed,
                                       mzCorrThresholdFolderIndexed,
                                       levelTwoPathAndFilename,
                                       logLevel );

        } // end check to see if there are enough data within this prelc region of interest( PLCROI )
    }  // end loop through all prelc regions of interest within the current scan group


    return 0;
}  // end workflowThree




///////////////////////////////////////////////////////////////////////////////
// Name:  workflowFour
//  The purpose of method is to implement workflow #4.
//  Workflow# 4 aggregates all the LC Candidates in the set of
//  directories generated by Workflow# 3.
//
//  This is a placeholder for a future non-grep data aggregation solution.
//  For the first version of MZDASoft, we will use grep to
//  aggregate all the LC Candidate Information from the MZDA Level 2 files.
//
// Return code: 0 successful.
//
/// Nelson.Ramirez@utsa.edu, 11/14/12
///////////////////////////////////////////////////////////////////////////////
int Workflow::workflowFour( cbi::ParameterLoader & paraObjGlobal ,
                           cbi::ParameterLoader & paraObjLocal ) {


    ///////////////////////////////////////////////////////////////////////////
    // Placeholder method
    //
    // workflow 4: LC Candidate Data Aggregation
    //
    // grep "MZDASoft found an LC Candidate" ./*/level2/*.m > lcc.txt
    //
    ///////////////////////////////////////////////////////////////////////////
    cout<<"Workflow 4: LC Candidate Data Aggregation"<<endl;


    return 0;
}  // end workflowFour

///////////////////////////////////////////////////////////////////////////////
//
// Private Helper methods
//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//
// Generate the local parameter files:
// i.e. parameterfilename00080.txt
// localworkerscansrange={2000:1:6000}
// localMinMz=2000
// localMaxMz=2020
// mzBinningThresholdInPPM=10
// scangroupid=80
//
// Also generate the set of PRLL job scheduling parameter files.
//
//
//
///////////////////////////////////////////////////////////////////////////////
int Workflow::generateParameterFiles(	string parmFilenameString, 
										string parmFilenameStringExtension,
                                        int fileNum,
										double minScan, 
										double maxScan,
										double startMZValue, 
										double endMZValue, 
										double mzDelta, 
										double mzOverlap,
                                        double mzBinningThresholdPPM){

    //cout<<"debug point: generateParameterFles: 1, mzDelta="<<mzDelta<<endl;
	///////////////////////////////////////////////////////////////////////////
	//
	//
	// Generate the output filenames
	// according to the corresponding scangroup id and width as specified
	// by the user.  Multiple scangroups could be processed by a single
	// worker, so the fundamental identifier of a unit of msda work is the 
	// scan group id.  
	// 
	//
	///////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////
	//
	// Generate a list of scan groups according to the parGlobalMinScan.at(0)
	// and parGlobalMaxScan.at(0) as well as the parGlobalMinMz.at(0) and
	// parGlobalMaxMz.at(0) and  parGlobalMinMz.at(0), as well as the 
	// parMzDelta and parMzDeltaOverlap.
	// 
	///////////////////////////////////////////////////////////////////////////

	
	///////////////////////////////////////////////////////////////////////////
	// Generate the sequence of start and end mz values for each 
	// distributed job.
	///////////////////////////////////////////////////////////////////////////
	vector<double> startMZ;
	vector<double> endMZ;
	double counter = 0;

    assert(mzDelta > 0);
    assert(startMZValue <= endMZValue );

	for ( counter = startMZValue; counter <= endMZValue; counter = counter + mzDelta ) {
		///////////////////////////////////////////////////////////////////////
		//
		// For each combination, we need to create a new parameter file 
		// for this combination.
		//
		///////////////////////////////////////////////////////////////////////
		startMZ.push_back(counter);
		endMZ.push_back(counter+mzDelta);
        //cout<<"debug point: generateParameterFles: 1a, counter="<<counter<<endl;
	}
    //cout<<"debug point: generateParameterFles: 2"<<endl;
	///////////////////////////////////////////////////////////////////////////
	// End mz range generator
	///////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////
	// 
	// Add the overlap calculation later... 
	// A second pass over the mz intervals 
	//
	///////////////////////////////////////////////////////////////////////////


    ///////////////////////////////////////////////////////////////////////////
    //
    // Generate a PRLL local parameter file-list:
    //  - When running via a PRLL interface, we need to generate a file
    //    containing the list of local parameter files.  The PRLL system
    //    will take care of managing the queueing of multiple processes
    //    within the same system so as to run a new process with a different
    //    item from the list.
    //  - Convention, we'll append the word list to one of the files
    //    in the same parameters directory.
    //
    //
    ///////////////////////////////////////////////////////////////////////////
    string pListFilename = parmFilenameString;
    pListFilename = pListFilename.append("list.txt");
    // cout<<pListFilename<<endl;
    //cout<<"debug point: generateParameterFles: 3"<<endl;
    ///////////////////////////////////////////////////////////////////////////
    //
    // Create the parameter file list
    //
    ///////////////////////////////////////////////////////////////////////////
    ofstream pListFile;  // the ofstream object owns the file.
    // it will handle cleanup if anything
    // happens after the file structure
    // has been allocated.
    pListFile.open(pListFilename.c_str());

     //cout<<"debug point: generateParameterFles: 4"<<endl;
	for ( vector<double>::size_type i = 0 ; i < startMZ.size(); i++ ){
		///////////////////////////////////////////////////////////////////////
		// Generate the output filenames
		// according to the corresponding scangroup id and width as specified
		// by the user.  Multiple scangroups could be processed by a single
		// worker, so the fundamental identifier of a unit of msda work is the 
		// scan group id.  
		// 
		///////////////////////////////////////////////////////////////////////
		int scangroupid= i;  // must convert from an integer to a string
		string valString; // string to hold the integer with leading zeros.
		///////////////////////////////////////////////////////////////////////
		// Use a string stream to hold the converted integer to the 
		// leading zero based string.
		///////////////////////////////////////////////////////////////////////
		stringstream scangroupidTemp (stringstream::in | stringstream::out);
		///////////////////////////////////////////////////////////////////////
		//
		// Set the fill option so that files get named with leading 
		// zeros.
		//
		///////////////////////////////////////////////////////////////////////
		scangroupidTemp << std::setw(5) << std::setfill('0');
		///////////////////////////////////////////////////////////////////////
		//
		// Send the current scangroupid information to the string stream 
		// so that we can generate the identification string.
		//
		///////////////////////////////////////////////////////////////////////
		scangroupidTemp << scangroupid;
		///////////////////////////////////////////////////////////////////////
		//
		// Move the string stream over to a holding string.
		//
		///////////////////////////////////////////////////////////////////////
		scangroupidTemp >> valString;

		///////////////////////////////////////////////////////////////////////
		//
		// Create the pre lc candidate filename
		//
		///////////////////////////////////////////////////////////////////////
		string pFilename = parmFilenameString;
		pFilename = pFilename.append(valString);
		pFilename = pFilename.append(parmFilenameStringExtension);
        //cout<<pFilename<<endl;

        ///////////////////////////////////////////////////////////////////////
        //
        // output the filename to the list file.  The PRLL queueing
        // system will automatically manage the processing of all
        // items within a list on a single multi-core system.
        ///////////////////////////////////////////////////////////////////////
        //
        // Add the section# and scangroup#, this is needed to be
        // able to dynamically remove level1 and level2 files are they
        // are no longer needed in the cluster running case.
        //
        ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////
        pListFile<<pFilename<<" "<<fileNum<<" "<<scangroupid<<endl;


		///////////////////////////////////////////////////////////////////////
		//
		// Create the internal data structure information file, each with
		// its scan group identifier.
		//
		///////////////////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////////////////
		//
		// Create the parameter file
		// 
		///////////////////////////////////////////////////////////////////////
		ofstream pFile;  // the ofstream object owns the file.
		// it will handle cleanup if anything 
		// happens after the file structure
		// has been allocated.
		pFile.open(pFilename.c_str());

		///////////////////////////////////////////////////////////////////////
		// Write out the contents of the file
		///////////////////////////////////////////////////////////////////////
		pFile<<"localworkerscansrange={"<<minScan<<":"<<1<<":"<<maxScan<<"}"<<endl;
        pFile<<"localminmz="<<startMZ.at(i)<<endl;
        pFile<<"localmaxmz="<<endMZ.at(i)<<endl;
        pFile<<"mzbinningthresholdinppm="<<mzBinningThresholdPPM<<endl;
		pFile<<"scangroupid="<<scangroupid<<endl;
        stringstream ss;
        ss<<(fileNum); // integer to string conversion
        string folderNumber = ss.str(); // converting stringstream to string object
        pFile<<"filenum="<<folderNumber<<endl; // internal identifier of the file we are processing
		///////////////////////////////////////////////////////////////////////
        // Close the parameter file
		///////////////////////////////////////////////////////////////////////
		pFile.close();
        //cout<<"debug point: generateParameterFles: 5,i="<<i<<endl;

	}  // end parameter file 


    ///////////////////////////////////////////////////////////////////////////
    // Close the parameter list file
    ///////////////////////////////////////////////////////////////////////////
    pListFile.close();
    //cout<<"debug point: generateParameterFles: 6"<<endl;
	return 0;
}  // end method




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////
//
// Generate mzdaWorkflow1 wrapper script
// This is the main script, managing the entire run.  It calls other
// dynamically generated scripts as appropriate.
//
//
//
// For example:
// # MZDASoft WorkFlow #1
// # The LD_LIBRARY_PATH must be set relative to the current
// # working directory
// cd ../
// export LD_LIBRARY_PATH=$(pwd)/product/lib:$(pwd)/cbilib/external/trng-4.13/lib:$(pwd)/cbilib/external/lib:$(pwd)/cbilib/external/gsl-1.5/lib:$LD_LIBRARY_PATH
// cd ./product
// ./mzda -m 1 -gp /home/nelson.ramirez/Desktop/softwareprojects/mzdasoftResearch/mzda/product/parameters.txt
// chmod +x ./*.sh
// ./mzdaWorkflow2.sh
//     --> Uses the PRLL shell based single system queue to utilize
//         as many cores as are available on the system to get through a
//         set of many large datasets.
//     --> The output of this stage is a set of level 1 data files.
//
//
//
//
//
//
//
///////////////////////////////////////////////////////////////////////////////
int Workflow::generateWorkflowOneScript(string wDir, string iDir, string cbiDir) {

    string outputScriptFilename1 = wDir;
    string outputScriptFilename2 = "mzdaWorkflow1.sh";
    string outputScriptFilename = outputScriptFilename1 + outputScriptFilename2;

    string header = "#MZDASoft Workflow #1";
    string ldLibraryPathBase1 = "export LD_LIBRARY_PATH=";
    string ldLibraryPathBase2 = iDir;
    string ldLibraryPathBase3 = "lib/:";
    string ldLibraryPathBase4 = cbiDir;
    string ldLibraryPathBase5 = "lib/:";
    string ldLibraryPathBase6 = "$LD_LIBRARY_PATH";
    string ldLibraryPath = ldLibraryPathBase1+ldLibraryPathBase2+ldLibraryPathBase3+ldLibraryPathBase4+ldLibraryPathBase5+ldLibraryPathBase6;

    string runCommandA1 = iDir;
    string runCommandA2 = "bin/mzda";
    string runCommandA3 = " -m 1";
    string runCommandA4 = " -gp ";
    string runCommandA5 = wDir;
    string runCommandA6 = "parameters.txt";

    string runCommandB1 = "chmod +x ";
    string runCommandB2 = wDir;
    string runCommandB3 = "*.sh";


    // generate calls to all workflow2 scripts, since there is one workflow2 script
    // for each data file being processed
    string runCommandC1 = wDir;
    string runCommandC2 = " ./mzdaWorkflow2.sh";

    string runCommand1 = runCommandA1 + runCommandA2 + runCommandA3 + runCommandA4 + runCommandA5 + runCommandA6;
    string runCommand2 = runCommandB1 + runCommandB2 + runCommandB3;
    string runCommand3 = runCommandC1 + runCommandC2;

    ///////////////////////////////////////////////////////////////////////////
    //
    // Create the run script
    //
    ///////////////////////////////////////////////////////////////////////////
    ofstream pFile;  // the ofstream object owns the file.
    // it will handle cleanup if anything
    // happens after the file structure
    // has been allocated.
    pFile.open(outputScriptFilename.c_str());
    cout<<"Making workflow one script"<<endl;
    pFile<<header<<endl;
    pFile<<ldLibraryPath<<endl;
    pFile<<runCommand1<<endl;
    pFile<<runCommand2<<endl;
    pFile.close();

    return 0;
}



///////////////////////////////////////////////////////////////////////////////
//
// Generate mzdaWorkflow2 wrapper script( A PRLL runscript )
//
//
//
//
// export LD_LIBRARY_PATH=/home/nelson.ramirez/testgit/mzda/product/lib:/home/nelson.ramirez/testgit/mzda/cbilib/external/trng-4.13/lib:/home/nelson.ramirez/testgit/mzda/cbilib/external/lib:/home/nelson.ramirez/testgit/mzda/cbilib/external/gsl-1.5/lib:/home/nelson.ramirez/testgit/mzda/cbilib/external/prll-0.6.2:$LD_LIBRARY_PATH
// export PATH=/home/nelson.ramirez/testgit/mzda/cbilib/external/prll-0.6.2:$PATH
// source /home/nelson.ramirez/testgit/mzda/cbilib/external/prll-0.6.2/prll.sh
// cd /home/nelson.ramirez/testgit/mzda/product
// prll -c 4 -s '/home/nelson.ramirez/testgit/mzda/product/mzda -m 2 -gp /home/nelson.ramirez/testgit/mzda/product/parameters.txt -lp "$1"' -p < /home/nelson.ramirez/testgit/mzda/parameterlist.txt
//
///////////////////////////////////////////////////////////////////////////////
int Workflow::generateWorkflowTwoScript( string wDir, string iDir, string cbiDir, string prllpath, string prllsource, int fileNum, int numCores ) {


    ///////////////////////////////////////////////////////////////////////////
    //
    // We need to make these completely general, only depending upon
    // a "working processing directory for the current run"
    //
    //
    //
    ///////////////////////////////////////////////////////////////////////////
    //we need to generate a workflow script for each file being processed
    //the allows separating the generation of level 1 files across more than 1 system in addition
    //to making use of the multi-core capabilities of the each individual system
    stringstream ss;
    ss<<(fileNum); // integer to string conversion
    string folderNumber = ss.str(); // converting stringstream to string object

    //generating the output
    string outputScriptFilename1 = wDir;
    string outputScriptFilename2 = "mzdaWorkflow2_"+folderNumber+".sh";
    string outputScriptFilename = outputScriptFilename1 + outputScriptFilename2;
    //create the LD_LIBRARY_PATH
    string ldLibraryPathBase1 = "export LD_LIBRARY_PATH=";
    string ldLibraryPathBase2 = iDir;
    string ldLibraryPathBase3 = "lib/:";
    string ldLibraryPathBase4 = cbiDir;
    string ldLibraryPathBase5 = "lib/:";
    string ldLibraryPathBase6 = "$LD_LIBRARY_PATH";
    string ldLibraryPath = ldLibraryPathBase1+ldLibraryPathBase2+ldLibraryPathBase3+ldLibraryPathBase4+ldLibraryPathBase5+ldLibraryPathBase6;


    string path1 = "export PATH=";
    string path2 = prllpath;
    string path3 = ":$PATH";
    string path = path1+path2+path3;

    string prllsourcepath1 = "source ";
    string prllsourcepath2 = prllsource;
    string prllsourcepath = prllsourcepath1+prllsourcepath2;

    string runDirectory1 = "cd ";
    string runDirectory2 = wDir;
    string runDirectory3 = folderNumber+"/";
    string runDirectory = runDirectory1 + runDirectory2 + runDirectory3;

    // build the prll path file
    //
    // prll path
    // number of cores
    // executable filename
    // mode
    // global parameter file
    // local parameter list file
    stringstream ssNumCores;
    ssNumCores<<numCores; // integer to string conversion
    string numCoresStr = ssNumCores.str(); // converting stringstream to string object


    string prllCommand1 = "prll";
    string prllCommand2 = " -c " + numCoresStr;
    string prllCommand3 = " -s '" + iDir + "bin/mzda";
    string prllCommand4 = " -m 2"; // the mode
    string prllCommand5 = " -gp " + wDir + "parameters.txt";
    string prllCommand6 = " -lp \"$1\"'";
    string prllCommand7 = " -p < ";
    string prllCommand8 = wDir+folderNumber+"/localparameterslist.txt";

    string prllCommand = prllCommand1+prllCommand2+prllCommand3+prllCommand4 + prllCommand5 + prllCommand6 + prllCommand7 + prllCommand8;

    ///////////////////////////////////////////////////////////////////////////
    //
    // Create the run script
    //
    ///////////////////////////////////////////////////////////////////////////
    ofstream pFile;  // the ofstream object owns the file.
    // it will handle cleanup if anything
    // happens after the file structure
    // has been allocated.
    pFile.open(outputScriptFilename.c_str());
    pFile<<ldLibraryPath<<endl;
    pFile<<path<<endl;
    pFile<<prllsourcepath<<endl;
    pFile<<runDirectory<<endl;
    pFile<<prllCommand<<endl;

    pFile.close();

    return 0;
}



///////////////////////////////////////////////////////////////////////////////
//
// Generate mzdaWorkflow3 wrapper script( A PRLL runscript )
//
//
//
//
// export LD_LIBRARY_PATH=/home/nelson.ramirez/testgit/mzda/product/lib:/home/nelson.ramirez/testgit/mzda/cbilib/external/trng-4.13/lib:/home/nelson.ramirez/testgit/mzda/cbilib/external/lib:/home/nelson.ramirez/testgit/mzda/cbilib/external/gsl-1.5/lib:/home/nelson.ramirez/testgit/mzda/cbilib/external/prll-0.6.2:$LD_LIBRARY_PATH
// export PATH=/home/nelson.ramirez/testgit/mzda/cbilib/external/prll-0.6.2:$PATH
// source /home/nelson.ramirez/testgit/mzda/cbilib/external/prll-0.6.2/prll.sh
// cd /home/nelson.ramirez/testgit/mzda/product
// prll -c 4 -s '/home/nelson.ramirez/testgit/mzda/product/msda -m 3 -gp /home/nelson.ramirez/testgit/mzda/product/parameters.txt -lp "$1"' -p < /home/nelson.ramirez/testgit/mzda/parameterlist.txt
//
///////////////////////////////////////////////////////////////////////////////
int Workflow::generateWorkflowThreeScript(string wDir, string iDir, string cbiDir, string prllpath, string prllsource, int fileNum, int numCores ) {

    ///////////////////////////////////////////////////////////////////////////
    //
    //
    // We need to make these completely general, only depending upon
    // a "working processing directory for the current run"
    //
    //
    //
    ///////////////////////////////////////////////////////////////////////////

    //we need to generate a workflow script for each file being processed
    //the allows separating the generation of level 1 files across more than 1 system in addition
    //to making use of the multi-core capabilities of the each individual system
    stringstream ss;
    ss<<(fileNum); // integer to string conversion
    string folderNumber = ss.str(); // converting stringstream to string object

    //generating the output
    string outputScriptFilename1 = wDir;
    string outputScriptFilename2 = "mzdaWorkflow3_"+folderNumber+".sh";
    string outputScriptFilename = outputScriptFilename1 + outputScriptFilename2;


    //create the LD_LIBRARY_PATH
    string ldLibraryPathBase1 = "export LD_LIBRARY_PATH=";
    string ldLibraryPathBase2 = iDir;
    string ldLibraryPathBase3 = "lib/:";
    string ldLibraryPathBase4 = cbiDir;
    string ldLibraryPathBase5 = "lib/:";
    string ldLibraryPathBase6 = "$LD_LIBRARY_PATH";
    string ldLibraryPath = ldLibraryPathBase1+ldLibraryPathBase2+ldLibraryPathBase3+ldLibraryPathBase4+ldLibraryPathBase5+ldLibraryPathBase6;


    string path1 = "export PATH=";
    string path2 = prllpath;
    string path3 = ":$PATH";
    string path = path1+path2+path3;

    string prllsourcepath1 = "source ";
    string prllsourcepath2 = prllsource;
    string prllsourcepath = prllsourcepath1+prllsourcepath2;

    string runDirectory1 = "cd ";
    string runDirectory2 = wDir;
    string runDirectory3 = folderNumber+"/";
    string runDirectory = runDirectory1 + runDirectory2 + runDirectory3;
    //
    // build the prll path file
    //
    // prll path
    // number of cores
    // executable filename
    // mode
    // global parameter file
    // local parameter list file
    //
    stringstream ssNumCores;
    ssNumCores<<numCores; // integer to string conversion
    string numCoresStr = ssNumCores.str(); // converting stringstream to string object


    string prllCommand1 = "prll";
    string prllCommand2 = " -c " + numCoresStr;
    string prllCommand3 = " -s '" + iDir + "bin/mzda";
    string prllCommand4 = " -m 3"; // the mode
    string prllCommand5 = " -gp " + wDir + "parameters.txt";
    string prllCommand6 = " -lp \"$1\"'";
    string prllCommand7 = " -p < ";
    string prllCommand8 = wDir+folderNumber+"/localparameterslist.txt";

    string prllCommand = prllCommand1+prllCommand2+prllCommand3+prllCommand4 + prllCommand5 + prllCommand6 + prllCommand7 + prllCommand8;

    ///////////////////////////////////////////////////////////////////////////
    //
    // Create the run script
    //
    ///////////////////////////////////////////////////////////////////////////
    ofstream pFile;  // the ofstream object owns the file.
    // it will handle cleanup if anything
    // happens after the file structure
    // has been allocated.
    pFile.open(outputScriptFilename.c_str());
    pFile<<ldLibraryPath<<endl;
    pFile<<path<<endl;
    pFile<<prllsourcepath<<endl;
    pFile<<runDirectory<<endl;
    pFile<<prllCommand<<endl;

    pFile.close();

    return 0;

}






///////////////////////////////////////////////////////////////////////////////
//
// Generate mzdaWorkflow4 wrapper script
// Workflow 4 will aggregate all mzda level 2 information onto
// a single file.
//
//
// grep -h "%LCC:"  /work/01759/nramirez/mzdataccwork/1/level2/*.m > /work/01759/nramirez/mzdataccwork/mzdaLevelThree1.dat
// wc /work/01759/nramirez/mzdataccwork/mzdaLevelThree1.dat > /work/01759/nramirez/mzdataccwork/job1summary.txt
// rm /work/01759/nramirez/mzdataccwork/1/level2/*.m
// rm /work/01759/nramirez/mzdataccwork/1/level1/*.csv
//
//
///////////////////////////////////////////////////////////////////////////////
int Workflow::generateWorkflowFourScript(string wDir, string iDir, string cbiDir, string prllpath, string prllsource, int fileNum, int numCores, bool levelOneClear, bool levelTwoClear ) {

    ///////////////////////////////////////////////////////////////////////////
    //
    //
    // We need to make these completely general, only depending upon
    // a "working processing directory for the current run"
    //
    //
    //
    ///////////////////////////////////////////////////////////////////////////

    //we need to generate a workflow script for each file being processed
    //the allows separating the generation of level 1 files across more than 1 system in addition
    //to making use of the multi-core capabilities of the each individual system
    stringstream ss;
    ss<<(fileNum); // integer to string conversion
    string folderNumber = ss.str(); // converting stringstream to string object

    //generating the output
    string outputScriptFilename1 = wDir;
    string outputScriptFilename2 = "mzdaWorkflow4_"+folderNumber+".sh";
    string outputScriptFilename = outputScriptFilename1 + outputScriptFilename2;


    //create the LD_LIBRARY_PATH
    string ldLibraryPathBase1 = "export LD_LIBRARY_PATH=";
    string ldLibraryPathBase2 = iDir;
    string ldLibraryPathBase3 = "lib/:";
    string ldLibraryPathBase4 = cbiDir;
    string ldLibraryPathBase5 = "lib/:";
    string ldLibraryPathBase6 = "$LD_LIBRARY_PATH";
    string ldLibraryPath = ldLibraryPathBase1+ldLibraryPathBase2+ldLibraryPathBase3+ldLibraryPathBase4+ldLibraryPathBase5+ldLibraryPathBase6;


    string path1 = "export PATH=";
    string path2 = prllpath;
    string path3 = ":$PATH";
    string path = path1+path2+path3;

    string prllsourcepath1 = "source ";
    string prllsourcepath2 = prllsource;
    string prllsourcepath = prllsourcepath1+prllsourcepath2;

    ///////////////////////////////////////////////////////////////////////////
    //
    // build the prll path file
    //
    // prll path
    // number of cores
    // executable filename
    // mode
    // global parameter file
    // local parameter list file
    ///////////////////////////////////////////////////////////////////////////
    stringstream ssNumCores;
    ssNumCores<<numCores; // integer to string conversion
    string numCoresStr = ssNumCores.str(); // converting stringstream to string object


    ///////////////////////////////////////////////////////////////////////////
    //
    // Create the run script
    //
    ///////////////////////////////////////////////////////////////////////////
    ofstream pFile;  // the ofstream object owns the file.
    // it will handle cleanup if anything
    // happens after the file structure
    // has been allocated.
    pFile.open(outputScriptFilename.c_str());
    pFile<<ldLibraryPath<<endl;
    pFile<<path<<endl;
    pFile<<prllsourcepath<<endl;



    ///////////////////////////////////////////////////////////////////////////
    //
    //We need to first extract the LCC's from the MZDA Level 2 files.
    //grep -h "%LCC:"  /work/01759/nramirez/mzdataccwork/1/level2/*.m > /work/01759/nramirez/mzdataccwork/mzdaLevelThree1.dat
    //
    ///////////////////////////////////////////////////////////////////////////
    string grepCommand1 = "grep -h ";
    string grepCommand2 = " \"%LCC:\" ";
    string grepCommand3 = wDir+folderNumber+"/level2/*.m > ";
    string grepCommand4 = wDir+"mzdaLevelThree"+ folderNumber + ".dat";
    string grepCommandResult = grepCommand1+grepCommand2+grepCommand3+grepCommand4;
    pFile<<grepCommandResult<<endl;

    ///////////////////////////////////////////////////////////////////////////
    //
    //We count the number of LCC's found for validation / test purposes.
    //wc /work/01759/nramirez/mzdataccwork/mzdaLevelThree1.dat > /work/01759/nramirez/mzdataccwork/job1summary.txt
    //
    //
    ///////////////////////////////////////////////////////////////////////////
    string wcCommand1 = "wc ";
    string wcCommand2 = wDir+"mzdaLevelThree"+folderNumber+".dat > ";
    string wcCommand3 = wDir+"sectionsummary"+folderNumber+".txt";
    string wcCommandResult = wcCommand1 + wcCommand2 + wcCommand3;
    pFile<<wcCommandResult<<endl;


    ///////////////////////////////////////////////////////////////////////////
    //We need to remove the temporary files in level1 and level2 file
    //folders.
    //In a future release we can make this an option. ( e.g. clearIntermediateData )
    //
    //We should at some point provide a user controllable option to allow
    //users to control whether to clear the intermediate MZDA Level 1 data
    //as well as the intermediate MZDA Level 2 data.
    //
    //The MZDA Level 1 and MZDA Level 2 data allow the visualization of the
    //algorithm.
    //
    //It is critical to provide the user with an option
    //
    //rm /work/01759/nramirez/mzdataccwork/1/level1/*.csv
    //rm /work/01759/nramirez/mzdataccwork/1/level2/*.m
    //
    ///////////////////////////////////////////////////////////////////////////
    // Remove intermediate MZDASoft Level 1 files
    if ( levelOneClear == true ){

        string checkCommand1 = "levelonefiles=$(ls ";
        string checkCommand2 = wDir+folderNumber+"/level1/mzda*.csv 2> /dev/null | wc -l)";
        string checkCommand3 = "if [ \"$levelonefiles\" != \"0\" ]";
        string checkCommand4 = "    then";
        string checkCommand5 = "    echo \"Some level 1 files to remove\"";
        string checkCommand6 = "fi";
        string checkCommand7 = "if [ \"$levelonefiles\" == \"0\" ]";
        string checkCommand8 = "    then";
        string checkCommand9 = "    echo \"No level 1 files generated, nothing to remove\"";
        string checkCommand10 = "fi";


        string rmCommand1 = "rm ";
        string rmCommand2 = wDir+folderNumber+"/level1/mzda*.csv";
        string rmCommandResult = rmCommand1 + rmCommand2;


        pFile<<checkCommand1<<checkCommand2<<endl<<checkCommand3<<endl<<checkCommand4<<endl<<checkCommand5<<endl<<"     "<<rmCommandResult<<endl<<checkCommand6<<endl;
        pFile<<checkCommand7<<endl<<checkCommand8<<endl<<checkCommand9<<endl<<checkCommand10<<endl;


        //string rmCommand1 = "rm ";
        //string rmCommand2 = wDir+folderNumber+"/level1/*.csv";
        //string rmCommandResult = rmCommand1 + rmCommand2;
        //pFile<<rmCommandResult<<endl;
    }

    // Remove intermediate MZDASoft Level 2 files
    if ( levelTwoClear == true ){

        ///////////////////////////////////////////////////////////////////////
        // This fix resolves an issue with the script ending
        // prematurely due to trying to remove files that do not exists.
        // we need to add a check to only remove files only if they exist.
        // We also add the same checks to MZDA Level 1 intermediate file
        // removal.
        ///////////////////////////////////////////////////////////////////////
        // For example:
        //        leveltwofiles=$(ls /work/01759/nramirez/mzdataccwork/54/level2/*.m 2> /dev/null | wc -l)
        //        if [ "$leveltwofiles" != "0" ]
        //         then
        //                echo "Some level 2 files to remove"
        //                rm /work/01759/nramirez/mzdataccwork/54/level2/*.m
        //        fi
        //        if [ "$leveltwofiles" == "0" ]
        //         then
        //                echo "No level 2 files generated, nothing to remove"
        //        fi
        ///////////////////////////////////////////////////////////////////////
        string checkCommand1 = "leveltwofiles=$(ls ";
        string checkCommand2 = wDir+folderNumber+"/level2/mzda*.m 2> /dev/null | wc -l)";
        string checkCommand3 = "if [ \"$leveltwofiles\" != \"0\" ]";
        string checkCommand4 = "    then";
        string checkCommand5 = "    echo \"Some level 2 files to remove\"";
        string checkCommand6 = "fi";
        string checkCommand7 = "if [ \"$leveltwofiles\" == \"0\" ]";
        string checkCommand8 = "    then";
        string checkCommand9 = "    echo \"No level 2 files generated, nothing to remove\"";
        string checkCommand10 = "fi";


        string rmCommand1 = "rm ";
        string rmCommand2 = wDir+folderNumber+"/level2/mzda*.m";
        string rmCommandResult = rmCommand1 + rmCommand2;


        pFile<<checkCommand1<<checkCommand2<<endl<<checkCommand3<<endl<<checkCommand4<<endl<<checkCommand5<<endl<<"     "<<rmCommandResult<<endl<<checkCommand6<<endl;
        pFile<<checkCommand7<<endl<<checkCommand8<<endl<<checkCommand9<<endl<<checkCommand10<<endl;


    }









    pFile.close();

    return 0;

}






//int Workflow::generateWorkflowFourScript(string wDir, string iDir, string cbiDir, int totalSections, int fileNum) {



//    stringstream ss;
//    ss<<(fileNum); // integer to string conversion
//    string folderNumber = ss.str(); // converting stringstream to string object

//    //generating the output
//    string outputScriptFilename1 = wDir;
//    string outputScriptFilename2 = "mzdaWorkflow4_"+folderNumber+".sh";
//    string outputScriptFilename = outputScriptFilename1 + outputScriptFilename2;


//    //create the LD_LIBRARY_PATH
//    string ldLibraryPathBase1 = "export LD_LIBRARY_PATH=";
//    string ldLibraryPathBase2 = iDir;
//    string ldLibraryPathBase3 = "lib/:";
//    string ldLibraryPathBase4 = cbiDir;
//    string ldLibraryPathBase5 = "lib/:";
//    string ldLibraryPathBase6 = "$LD_LIBRARY_PATH";
//    string ldLibraryPath = ldLibraryPathBase1+ldLibraryPathBase2+ldLibraryPathBase3+ldLibraryPathBase4+ldLibraryPathBase5+ldLibraryPathBase6;


//    string path1 = "export PATH=";
//    string path2 = prllpath;
//    string path3 = ":$PATH";
//    string path = path1+path2+path3;

//    string prllsourcepath1 = "source ";
//    string prllsourcepath2 = prllsource;
//    string prllsourcepath = prllsourcepath1+prllsourcepath2;

//    string runDirectory1 = "cd ";
//    string runDirectory2 = wDir;
//    string runDirectory3 = folderNumber+"/";
//    string runDirectory = runDirectory1 + runDirectory2 + runDirectory3;





//    //generating the output
//    string outputScriptFilename1 = wDir;
//    string outputScriptFilename2 = "mzdaWorkflow4.sh";
//    string outputScriptFilename = outputScriptFilename1 + outputScriptFilename2;


//    ///////////////////////////////////////////////////////////////////////////
//    // Generate an MZDA Level 3 format file: A summary of LC Candidates
//    // An MZDA Level 3 format file
//    ///////////////////////////////////////////////////////////////////////////
//    //grep "MZDASoft found an LC Candidate" ./*/level2/*.m > lcc.txt

//    ///////////////////////////////////////////////////////////////////////////
//    // Generate complete aggregated data file
//    ///////////////////////////////////////////////////////////////////////////
//    string prllCommand1 = "grep -h ";
//    string prllCommand2 = "\"%LCC:\" ";
//    string prllCommand3 = " " + wDir + "*/level2/*.m > ";
//    string prllCommand4 = wDir + "mzdaLevelThree.dat";
//    string prllCommand = prllCommand1+prllCommand2+prllCommand3+prllCommand4;


//    ///////////////////////////////////////////////////////////////////////////
//    // Generate each level 3 file separately, one file per
//    // section.  As the data increases in size,it is important to
//    // be able to generate the MZDA Level 3 files separately.
//    ///////////////////////////////////////////////////////////////////////////
//    vector<string> prllCommandResultVec;


//    for ( int i = 1; i <= totalSections; i++ ) {
//        stringstream sectionId;
//        sectionId<<i; // integer to string conversion
//        string sectionIdStr = sectionId.str(); // converting stringstream to string object

//        string prllCommand1VecTemp = "grep -h ";;
//        string prllCommand2VecTemp = "\"%LCC:\" ";
//        string prllCommand3VecTemp = " " + wDir + sectionIdStr + "/level2/*.m > ";
//        string prllCommand4VecTemp = wDir + "mzdaLevelThree"+sectionIdStr+".dat";
//        string prllCommandResultVecTemp = prllCommand1VecTemp + prllCommand2VecTemp + prllCommand3VecTemp + prllCommand4VecTemp;

//        prllCommandResultVec.push_back(prllCommandResultVecTemp);

//    }

//    ///////////////////////////////////////////////////////////////////////////
//    //
//    // Create the run script
//    //
//    ///////////////////////////////////////////////////////////////////////////
//    ofstream pFile;  // the ofstream object owns the file.
//    // it will handle cleanup if anything
//    // happens after the file structure
//    // has been allocated.
//    pFile.open(outputScriptFilename.c_str());

//    pFile<<prllCommand<<endl;

//    for ( int i = 0; i < prllCommandResultVec.size(); i++ ){
//        pFile<<prllCommandResultVec.at(i)<<endl;
//    }

//    pFile.close();

//    return 0;
//}






///////////////////////////////////////////////////////////////////////////////
//
// Generate cluster worker workflow script
//
//
///////////////////////////////////////////////////////////////////////////////
//int Workflow::generateWorkflowClusterScript(string wDir, string iDir, string cbiDir, string prllpath, string prllsource, int fileNum, int numCores ) {

//    ///////////////////////////////////////////////////////////////////////////
//    //
//    //
//    // We need to make these completely general, only depending upon
//    // a "working processing directory for the current run"
//    //
//    //
//    //
//    ///////////////////////////////////////////////////////////////////////////

//    //we need to generate a workflow script for each file being processed
//    //the allows separating the generation of level 1 files across more than 1 system in addition
//    //to making use of the multi-core capabilities of the each individual system
//    stringstream ss;
//    ss<<(fileNum); // integer to string conversion
//    string folderNumber = ss.str(); // converting stringstream to string object

//    //generating the output
//    string outputScriptFilename1 = wDir;
//    string outputScriptFilename2 = "mzdaWorkflowCluster_"+folderNumber+".sh";
//    string outputScriptFilename = outputScriptFilename1 + outputScriptFilename2;


//    //create the LD_LIBRARY_PATH
//    string ldLibraryPathBase1 = "export LD_LIBRARY_PATH=";
//    string ldLibraryPathBase2 = iDir;
//    string ldLibraryPathBase3 = "lib/:";
//    string ldLibraryPathBase4 = cbiDir;
//    string ldLibraryPathBase5 = "lib/:";
//    string ldLibraryPathBase6 = "$LD_LIBRARY_PATH";
//    string ldLibraryPath = ldLibraryPathBase1+ldLibraryPathBase2+ldLibraryPathBase3+ldLibraryPathBase4+ldLibraryPathBase5+ldLibraryPathBase6;


//    string path1 = "export PATH=";
//    string path2 = prllpath;
//    string path3 = ":$PATH";
//    string path = path1+path2+path3;

//    string prllsourcepath1 = "source ";
//    string prllsourcepath2 = prllsource;
//    string prllsourcepath = prllsourcepath1+prllsourcepath2;

//    string runDirectory1 = "cd ";
//    string runDirectory2 = wDir;
//    string runDirectory3 = folderNumber+"/";
//    string runDirectory = runDirectory1 + runDirectory2 + runDirectory3;
//    //
//    // build the prll path file
//    //
//    // prll path
//    // number of cores
//    // executable filename
//    // mode
//    // global parameter file
//    // local parameter list file
//    //
//    stringstream ssNumCores;
//    ssNumCores<<numCores; // integer to string conversion
//    string numCoresStr = ssNumCores.str(); // converting stringstream to string object


//    ///////////////////////////////////////////////////////////////////////////
//    // Create the level 1 cluster workflow
//    ///////////////////////////////////////////////////////////////////////////
//    string prllCommand3_m2 = iDir + "bin/mzda";
//    string prllCommand4_m2 = " -m 2"; // the mode
//    string prllCommand5_m2 = " -gp " + wDir + "parameters.txt";
//    string prllCommand6_m2 = " -lp $prll_arg_1;";
//    string prllCommand_m2 = prllCommand3_m2+prllCommand4_m2 + prllCommand5_m2 + prllCommand6_m2;

//    ///////////////////////////////////////////////////////////////////////////
//    // Create the level 2 cluster workflow
//    ///////////////////////////////////////////////////////////////////////////
//    string prllCommand3_m3 = iDir + "bin/mzda";
//    string prllCommand4_m3 = " -m 3"; // the mode
//    string prllCommand5_m3 = " -gp " + wDir + "parameters.txt";
//    string prllCommand6_m3 = " -lp $prll_arg_1;";
//    string prllCommand_m3 = prllCommand3_m3+prllCommand4_m3 + prllCommand5_m3 + prllCommand6_m3;

//    ///////////////////////////////////////////////////////////////////////////
//    // Create the level 3 cluster workflow
//    ///////////////////////////////////////////////////////////////////////////
//    string grepCommand1 = "grep -h";
//    string grepCommand2 = " \"%LCC:\" " + wDir + "$prll_arg_2/level2/mzdaTwo$prll_arg_3*.m > ";
//    string grepCommand3 = wDir + "mzdaLevelThree.$prll_arg_2.$prll_arg_3.dat";
//    string grepCommand = grepCommand1+grepCommand2+grepCommand3;

//    ///////////////////////////////////////////////////////////////////////////
//    //
//    // Level 1 file find and remove.  This is needed to free up
//    // disk space in large runs.
//    //
//    ///////////////////////////////////////////////////////////////////////////
//    string findLevelOneCommand1 = "find ";
//    string findLevelOneCommand2 = wDir + "$prll_arg_2/level1/mzdaOne*0$prll_arg_3.csv -exec rm {} \\;";
//    string findLevelOneCommand = findLevelOneCommand1 + findLevelOneCommand2;

//    ///////////////////////////////////////////////////////////////////////////
//    //
//    // Level 2 file find and remove.  This is needed to free up
//    // disk space in large runs.
//    //
//    ///////////////////////////////////////////////////////////////////////////
//    string findLevelTwoCommand1 = "find ";
//    string findLevelTwoCommand2 = wDir + "$prll_arg_2/level2/mzdaTwo$prll_arg_3*.m -exec rm {} \\;";
//    string findLevelTwoCommand = findLevelTwoCommand1  + findLevelTwoCommand2;

//    ///////////////////////////////////////////////////////////////////////////
//    //
//    // Build the PRLL Command
//    //
//    ///////////////////////////////////////////////////////////////////////////
//    string prllCommand1 = "prll";
//    string prllCommand2 = " -c " + numCoresStr;
//    string prllCommand3 = " -s 'mzdaclusterworker \"$1\";'";
//    string prllCommand4 = " -p < ";
//    string prllCommand5 = wDir+folderNumber+"/localparameterslist.txt";
//    string prllCommand = prllCommand1+prllCommand2+prllCommand3+prllCommand4 + prllCommand5;

//    ///////////////////////////////////////////////////////////////////////////
//    //
//    // Create the run script
//    //
//    ///////////////////////////////////////////////////////////////////////////
//    ofstream pFile;  // the ofstream object owns the file.
//    // it will handle cleanup if anything
//    // happens after the file structure
//    // has been allocated.
//    pFile.open(outputScriptFilename.c_str());
//    pFile<<ldLibraryPath<<endl;
//    pFile<<path<<endl;
//    pFile<<prllsourcepath<<endl;
//    pFile<<runDirectory<<endl;

//    ///////////////////////////////////////////////////////////////////////////
//    //
//    // Create the cluster node workflow
//    //
//    ///////////////////////////////////////////////////////////////////////////
//    pFile<<endl<<"mzdaclusterworker() {"<<endl;
//    pFile<<"prll_splitarg"<<endl;
//    pFile<<"echo $prll_arg_1 # path to the local parameter file list"<<endl;
//    pFile<<"echo $prll_arg_2 # the section number, same as the folder name"<<endl;
//    pFile<<"echo $prll_arg_3 # the scan group identifier"<<endl;
//    pFile<<"# run mzda level 1 workflow"<<endl;
//    pFile<<prllCommand_m2<<endl;
//    pFile<<"# run mzda level 2 workflow"<<endl;
//    pFile<<prllCommand_m3<<endl;
//    pFile<<"# run mzda level 3 workflow"<<endl;
//    pFile<<"#"<<grepCommand<<endl;
//    pFile<<"#Find level 1 files that we can safely remove"<<endl;
//    pFile<<"#"<<findLevelOneCommand<<endl;
//    pFile<<"#Find level 2 files that we can safely remove"<<endl;
//    pFile<<"#"<<findLevelTwoCommand<<endl;
//    pFile<<"}"<<endl;
//    pFile<<prllCommand<<endl;




//    pFile.close();

//    return 0;

//}








///////////////////////////////////////////////////////////////////////////////
//
// Add new workflows here.
//
///////////////////////////////////////////////////////////////////////////////






} // End namespace msda

