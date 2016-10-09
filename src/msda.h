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
///@file CommandLineProcessor.h
/// 
///@brief Part of CBILIB HPC software library
/// A command line processor object makes it easier to manage many 
/// command line parameters.  It organizes the parameters by user flag,
/// and within each flag, it parses and expands the contents of each
/// parameter set associated with a given flag using regular expressions.
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
///
/// External Libraries Used: 
/// 1) mstoolkit, Ramp
///    license = BSD, link=  http://code.google.com/p/mstoolkit/
///
///////////////////////////////////////////////////////////////////////////////
#ifndef MSDA_H_
#define MSDA_H_


///////////////////////////////////////////////////////////////////////////////
//
// MSDA Serial and Parallel Workflow
//
//  A domain decomposition based parallelization approach
//
//
//
//    GUI -->
//
//    Global Parameter File -->
//
//    Parallel Manager[msdaParallelManager.exe] --> 
//			Generate a set of overlapping scan groups that can
//          be processed completely in parallel.
// 
//			{jobscript####.job, localworkerparameters####.txt}  -->
//		
//
//    	GlobalParameters.txt,LocalParameters.txt -->    [ msda.exe ] -->
//
//		--> LocalResultsFile####.txt
//
//
//	  msda Results Merger[msdaResultsMerge.exe]  --> 
//      --> Takes all LocalResultsFiles and produces a final
//          merged results file.
//      --> For different data sets, this will be different.
//          mass#1..#N for each system
//
//
//
//
///////////////////////////////////////////////////////////////////////////////
#ifdef OPENMP
#include <omp.h>
#endif
///////////////////////////////////////////////////////////////////////////////
// CBI Library Includes(CBI Namespace)
// The following headers are from the CBI Library
// Each header file has specific copyright information
///////////////////////////////////////////////////////////////////////////////
#include "CommandLineProcessor.h"  // Process Command Line Parameters
#include "ParameterLoader.h"  // Definiitions for a single parameter
#include "Parameter.h"  // A parameter handling class 
#include "MathUtils.h"  // Utilities aimed at generally useful math algorithms
#include "SignalProcessor.h" // Signal processing utilities such as filters
#include "Histogram.h"  // Histogram processing utilities

///////////////////////////////////////////////////////////////////////////////
//
// External libraries( A variety of Namespaces )
// 
//
///////////////////////////////////////////////////////////////////////////////
#include "mzParser.h"  // Ramp library for mzXML reading
#include <Eigen/Dense> // Matrix library ( dense matrices )

///////////////////////////////////////////////////////////////////////////////
//
// mzdasoft Specific Includes(msda Namespace)
//
///////////////////////////////////////////////////////////////////////////////
#include "Workflow.h"
#include "MsDataLoader.h"    // Handles loading a subset of the total 

///////////////////////////////////////////////////////////////////////////////
//
// MSDA Logging
//
// Logging is needed for 
// Debugging, runtime stats, records, traceability.
//
// 
///////////////////////////////////////////////////////////////////////////////
#ifdef __MSDALOGGING__
#define MSDALOGINGDIRECTORY "./"
#define MSDALOGFILE "./msdalogfile.log"   
#endif



///////////////////////////////////////////////////////////////////////////////
//
// GUI Design:
//
// 
///////////////////////////////////////////////////////////////////////////////
//The goal of the GUI will be to generate a set of job scripts that can be 
// submitted to a one or more compute systems.
//This design therefore provides for the following key design points:
//1) Loosely coupled modules
//2) Data Domain decomposition parallelism is controlled and validated by 
//   the user.  In the case of MSDA, the user
//   controls how one or more mzXML LCMS data files are to be processed.  
//   The GUI provides the user with a convenient way to generate the 
//   specifications about how the work is to be parallelized.
//   As well as provides the user with control on a variety of 
//   domain and data specific parameters.
//3) User interface platform independence.  The GUI should be completely 
//   separate from the compute work.  The GUI should therefore be able to be 
//   run from any operating system ( Linux, Windows, Mac ) and "connect" to the 
//   computational part via a loosely couple interface.   The interface to be 
//   used is simply a set of job scripts( files ), that will be given to the 
//   computational engine to process.  These jobs scripts contain all the 
//   information needed to perform all of the computation according to 
//   the user specified parameters.
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
// Research References:
// [1] http://en.wikipedia.org/wiki/Mass_%28mass_spectrometry%29
// [2] arxiv.org/pdf/1101.1154 ( LCMS Overview )
// [3] http://omics.pnl.gov/software/Decon2LS.php
// [4] http://www.biomedcentral.com/1471-2105/12/74  ( MRCQuant )
// [5] http://www.chm.bris.ac.uk/ms/faq.html
///////////////////////////////////////////////////////////////////////////////

#endif // msda_H_ 
