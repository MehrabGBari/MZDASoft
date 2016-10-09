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
///@file Workflow.h
/// 
///@brief Part of MSDA/MzDa Software Infrastructure 
/// A Workflow object manages the various runtimes modes of the 
/// software.  The Workflow object is where everything comes
/// together. 
/// 
///@version 1.0.0
///
/// Change History:
/// Feature ID, Version, Date, Description, ADD,CHG,MOV DeveloperID
/// ft00000, vr1.0.0, 04/11/12, Initial Version, @ADDNRZ
///
///
/// Developer Contact Information:
/// Developer ID, Name, Email
/// NRZ, Nelson Ramirez, Nelson.Ramirez@utsa.edu
///
/// Language: C++
/// 
/// Documentation: Doxygen
///////////////////////////////////////////////////////////////////////////////
#ifndef WORKFLOW_H_
#define WORKFLOW_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "ParameterLoader.h"
#include "MsDataLoader.h"
#include "ScanGroup.h"
#include "SignalProcessor.h"

using namespace std;


///////////////////////////////////////////////////////////////////////////////
///
/// @namespace msda
/// @brief msda
///
///
///////////////////////////////////////////////////////////////////////////////
namespace msda{


///////////////////////////////////////////////////////////////////////////////
/// @class Workflow
/// @brief This class manages the runtime environment
///
///
///
///////////////////////////////////////////////////////////////////////////////
class Workflow {
public:
	///////////////////////////////////////////////////////////////////////////
	/// @fn void Workflow()
	/// @brief Constructor
	/// @param Void
	/// @exception None
	/// @return Void
	///////////////////////////////////////////////////////////////////////////
	Workflow();
	
	///////////////////////////////////////////////////////////////////////////
	/// @fn void ~Workflow()
	/// @brief Destructor
	/// @param Void
	/// @exception None
	/// @return Void
	///////////////////////////////////////////////////////////////////////////
	~Workflow();
	

	///////////////////////////////////////////////////////////////////////////
	//
    // Workflow Generators
	//
	// Nelson.Ramirez@utsa.edu, 04/16/12
	///////////////////////////////////////////////////////////////////////////
    int workflowZero(cbi::ParameterLoader & globalParameters);
	int workflowOne(cbi::ParameterLoader & globalParameters);
    int workflowTwo( cbi::ParameterLoader & paraObjGlobal,
                     cbi::ParameterLoader & paraObjLocal,
                     msda::ScanGroup & sgObject,
                     msda::MsDataLoader & msdlObject,
                     cbi::SignalProcessor & spObject);

    int workflowThree( cbi::ParameterLoader & paraObjGlobal,
                     cbi::ParameterLoader & paraObjLocal,
                     msda::ScanGroup & sgObject,
                     msda::MsDataLoader & msdlObject,
                     cbi::SignalProcessor & spObject);
	
    // aggregate all lc candidates from MZDA Level 2 files
    // this is a placeholder for a non-grep solution.
    int workflowFour( cbi::ParameterLoader & paraObjGlobal ,
                      cbi::ParameterLoader & paraObjLocal);






private:	

	///////////////////////////////////////////////////////////////////////////
	/// Private methods
	///////////////////////////////////////////////////////////////////////////		

	int generateParameterFiles(	string parmFilenameString, 
								string parmFilenameStringExtension,
                                int fileNum,
								double minScan, 
								double maxScan,
								double startMZValue, 
								double endMZValue, 
								double mzDelta, 
								double mzOverlap,
                                double mzBinningThresholdPPM);


    int generateWorkflowOneScript(string wDir, string iDir, string cbiDir ); // Sets up work distribution for top level domain decomposition into scangroups
    int generateWorkflowTwoScript(string wDir, string iDir, string cbiDir, string prllpath, string prllsourcepath,int fileNum, int numCores); // Uses PRLL shell queueing system for binning
    int generateWorkflowThreeScript(string wDir, string iDir, string cbiDir, string prllpath, string prllsource, int fileNum, int numCores); // Uses PRLL shell queueing system for lc candidate generation
    int generateWorkflowFourScript(string wDir, string iDir, string cbiDir, string prllpath, string prllsource, int fileNum, int numCores, bool levelOneClear, bool levelTwoClear ); // Uses grep tool for aggregating all level 2 data ( lc candidates ) into a single file

    ///////////////////////////////////////////////////////////////////////////
    // Call the workflowCluster function to create the workflow files
    // for a cluster run with large numbers of input files.
    ///////////////////////////////////////////////////////////////////////////
    //int generateWorkflowClusterScript(string wDir, string iDir, string cbiDir, string prllpath, string prllsource, int fileNum, int numCores); // Uses PRLL shell queueing system for lc candidate generation

    ///////////////////////////////////////////////////////////////////////////
    // Create a workflow for starting web services running
    // over the query indices. ( potential interfacing idea )
    ///////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////
	/// Private variables
	///////////////////////////////////////////////////////////////////////////							

};  /// End Class Declaration

} /// End Namespace msda


#endif /// WORKFLOW_H_ 
