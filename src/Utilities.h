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
///
/// Utilities.h
///
///  Created on: Aug 8, 2011
///      Author: nelson.ramirez
///
///  Contains utilities that are general to 
///  any part of the MSDA processing pipeline.
///
///  This will signal processing routines 
///  that are required for the Mass Spectrometry
///  data analysis software.
///  
///  In addition, ports of signal processing
///  algorithms developed by M.Zhang. will
///  be placed in this class.
///
/// Software Developer: Nelson.Ramirez@utsa.edu
///
///////////////////////////////////////////////////////////////////////////////
#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <vector>
#include <fstream>
#include "MathUtils.h"
#include "SignalProcessor.h"
#include "ParameterLoader.h"



using namespace std;

namespace msda {

class Utilities {
public:
	Utilities();
	~Utilities();






    ///////////////////////////////////////////////////////////////////////////
    // Interval clustering algorithm
    //
    // Custom algorithm developed to create a cluster of the interval
    // peaks.
    //
    // - In the future, this function can be updated to call any clustering
    //   algorithm. This is a generic interface to processing
    //
    // Step 1: mz*rt
    // Step 2: sort (mz*rt) ascending
    // Step 3: break the bins according to a data resolution dependent parameter
    // Step 4: each group represents an LC Candidate.
    //
    // Input: An intervalGraph data structure, consisting of the following data
    //   1.  label for each data point( maps to the global interval data structure)
    //   2.  mz data ( corresponds to the window mz value ) of the interval peak
    //   3.  rt data ( corresponds to the retention time information) of the interval peak
    //
    // Output:
    //   1. groups  ( Note:  This must be an empty vector )
    //      groups[0] --> label#, label#, label#
    //      groups[1] --> label#, label#
    //      groups[2] --> label#, label#, label#, label#, label#
    //
    //  The label# allows us to go back to the global data structure and get
    //  all the detailed information about the interval to which this
    //
    //
    // Return codes:
    //
    //
    //
    // Developed by:
    // Michelle.Zhang@utsa.edu, Nelson.Ramirez@utsa.edu 10/3/2012.
    // Copyright CBI,UTSA, All rights reserved.  2012.
    //
    ///////////////////////////////////////////////////////////////////////////
    int intervalCluster(     vector<int> & peakIntervalId,
                             vector<double> & peakIntervalMz,
                             vector<double> & peakIntervalRt,
                             vector<double> & peakIntervalIntensity,
                             vector< vector<int> >  & groupsIdtoIntervalIdVectorsResult,
                             vector< vector<double> > & groupsIdtoIntervalMzVectorsResult,
                             vector< vector<double> > & groupsIdtoIntervalRtVectorsResult,
                             vector< vector<double> > & groupsIdtoIntervalIntensityVectorsResult,
                             vector<int> &  intervalIdtoGroupIdMappingResult);

    ///////////////////////////////////////////////////////////////////////////
    // Version 2 of a peak cluster algorithm,
    // Based on matlab code from M.Zhang.10/08/12
    //
    //  The input to this routine are a set of interval peaks
    //  generated after the following key stages:
    //
    //  1) Interval Detection in the scan ( retention time ) dimension
    //  2) Interval Splitting in the scan ( retention time ) dimension
    //  3) Maximum intensity value finding within each interval in retention time
    //     dimension.
    //
    //  The algorithm proceeds as follows:
    //   - Use the preLCPeakApexMap and the preLCPeakApexIntensity
    //     matrices to generate a runnig sum of groups of columns.
    //     The number of neighbor columns is set '1' by default.
    //
    //     This means that we're going to collapse groups of 3
    //     columns( except at the first and last column which will
    //     only consist of 2 columns.
    //
    //   - Interval detection will then be applied to this collapsed
    //     signal in the mz dimension.
    //
    //   - For each interval in the mz dimension, the peak in the
    //     mz dimension interval is determined.
    //
    //   - This will allow us to determine the 2 dimensional location
    //     of the mz/rt peak.
    //
    //
    //
    //  Initial default values:
    //  int LCPeakApexTolerance = 1
    //  double mzIntervalThreshold = 0
    //  double mzIntervalMinLength = 2
    //
    //
    //
    //
    //
    ///////////////////////////////////////////////////////////////////////////
    int intervalClusterVersion2(      vector<int> & peakIntervalId,
                                      vector<double> & peakIntervalMz,
                                      vector<double> & peakIntervalRt,
                                      vector<double> & peakIntervalIntensity,
                                      vector< vector<int>  > & preLCPeakApexMap,
                                      vector< vector<double> > & preLCPeakApexIntensity,
                                      int LCPeakApexTolerance,
                                      double mzIntervalThreshold,
                                      double mzIntervalMinLength,
                                      vector< vector<int>  > & LCPeakApexMap,
                                      vector< vector<double> > & LCPeakApexIntensity,
                                      vector<int> & apexWindowId,
                                      vector<int> & apexScanId,
                                      vector<double> & apexMz,
                                      vector<double> & apexRt,
                                      vector<int> & apexScanNum,
                                      vector<double> & apexIntensity,
                                      vector< vector< vector<int> > > & preLCPeakApexIntervalStartIndex,
                                      vector< vector< vector<int> > > & preLCPeakApexIntervalEndIndex,
                                      vector< vector< vector<double> > > & preLCPeakApexIntervalProfile,
                                      vector< vector< vector<int> > > & LCPeakApexIntervalStartIndex,
                                      vector< vector< vector<int> > > & LCPeakApexIntervalEndIndex );


    ///////////////////////////////////////////////////////////////////////////
    // Name:  logKLDistance
    // The purpose of this method is to calculate the KL Distance of
    // 2 vectors.
    // The result is the log KL Value.
    // This function is used during the Peptide Candidate generation
    // process.
    //
    // Parameters:
    //
    // Input:
    //    1) inputA: reference to a const vector of doubles
    //    2) inputB: reference to a const vector of doubles
    //
    // Output:log of the KL distance
    //     When the input vectors contain the same data:
    //     3) result = -10000000; ( instead of returning -Inf...., since log(0) = -Inf)
    //     4) equal = true;
    //
    //     When the input vector are not the same:
    //     3) result = a double value
    //     4) equal = false;
    //
    //
    // Return code: 0 successful
    //              -1(inputA size is 0)
    //              -2(inputB size is 0)
    //              -3(input sizes do not match)
    //              -4(inputA contains a negative value)
    //              -5(inputB contains a negative value)
    //              -6(divide by zero)
    //              -7(divide by zero)
    //
    // Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 06/15/12
    // Original matlab function: KL_calculate(Vector01,Vector02) by M.Zhang
    ///////////////////////////////////////////////////////////////////////////
    int logKLDistance(   const vector<double> & inputA,
                         const vector<double> & inputB,
                         double & result,
                         bool & equal);



    ///////////////////////////////////////////////////////////////////////////
    // Name:  calcIsoPatternAvg
    // The purpose of this method is to generate the average
    // isotope pattern
    //
    // This function is used during the Peptide Candidate generation
    // process.
    //
    // For a given mass & maximum number of isotopes, there is an
    // average theoretical average isotope pattern.
    //
    // Parameters:
    //
    // Input:
    //    1) mass: reference to a const double
    //    2) maxIso: reference to a const int. Maximum isotope number to
    //               consider.
    //
    // Output:
    //    1) averageIsotopePattern: reference to a vector of doubles
    //
    //
    // Return code: 0 successful
    //              -1
    //
    // Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 06/19/12
    // Original matlab function:
    //   isopatternavg=calculateisopatternavg(mass,maxiso) by M.Zhang
    ///////////////////////////////////////////////////////////////////////////
    int calcIsoPatternAvg(  const double & mass, const int & maxiso, vector<double> & averageIsotopePattern);



    ///////////////////////////////////////////////////////////////////////////
    // Name:  checkModelfitness
    //
    // The purpose of this method is to determine the quality of the
    // quantification process.
    //
    // Key cases:
    // 1) Labeled SILAC
    //
    //
    // This function is used after the quantification process.
    //
    // Parameters:
    //
    // Input:
    //    1) const int labelingId
    //    2) const vector< vector<double> > & adjustedElutionProfiles
    //
    // Output:
    //    1)
    //
    //
    // Return code:
    //
    // Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 06/19/12
    // Original matlab function:
    //   [peptideFitnessScore]=checkModelfitness(peptideCandidateListItem, labelingMethod,isoPatternNatural)
    //    by M.Zhang
    ///////////////////////////////////////////////////////////////////////////
    int checkModelFitness( const int labelingId,
                           const vector< vector<double> > & adjustedElutionProfiles,
                           const vector<double> & massTemplate,
                           const vector<double> & labelingMethodMassShiftVector,
                           const vector<double> & isoPatternNatural);

	///////////////////////////////////////////////////////////////////////////
	//
	//
	// Vector differential calculation algorithm:
	// Original algorithm developed in Matlab by: Jianqiu Zhang
	//
	// Parameters:
	// Input: xic, a reference to a vector of doubles
	//        direction, an integer indicating whether left(-1) or right(1) differencing
	// Output: 
	// Return code: 0 successful, -1 (xic length error), -2(direction error)
	//
	// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 01/25/12
	//
	///////////////////////////////////////////////////////////////////////////
	int getDiff(vector<double> & xic, int direction, vector<double> & diffxic);

	///////////////////////////////////////////////////////////////////////////
	//
	//
	// Vector differential calculation algorithm:
	// Original algorithm developed in Matlab by: Jianqiu Zhang
	//
	// Parameters:
	// Input: xic, a reference to a vector of ints
	//        direction, an integer indicating whether left(-1) or right(1) differencing
	// Output: 
	// Return code: 0 successful, -1 (xic length error), -2(direction error)
	//
	// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 02/20/12
	//
	///////////////////////////////////////////////////////////////////////////
	int getDiff(vector<int> & xic, int direction, vector<int> & diffxic);

	///////////////////////////////////////////////////////////////////////////
	//
	//
	// getR2 correlation statistic algorithm:
	// Original algorithm developed in Matlab by: Jianqiu Zhang
	// Parameters:
	// Input: datavector1: a reference to a vector of doubles 
	//        datavector2: a reference to a vector of doubles
	// Output: 
	//        rcorr: a reference to a double where the result will be placed.
	//
	//
	// Return code: 0 successful, -1 (xic length error), -2(direction error)
	//
	// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 01/27/12
	//
	///////////////////////////////////////////////////////////////////////////
	int getR2Statistic(vector<double> & datavector1, vector<double> & datavector2, double & rcorr);

	///////////////////////////////////////////////////////////////////////////
	//
	// getNoiseThreshold statistic algorithm:
	//
	// Original algorithm developed in Matlab by: Jianqiu Zhang
	//  function noiseThreshold=getNoiseThreshold(smoothXICs,XICs,noiseThresholdLevel);
	//            noiseVector=smoothXICs-XICs;
	//            noiseMean=mean(noiseVector(noiseVector>0));  (  use only those values that are greater than 0 )
	//            noiseStd=std(noiseVector(noiseVector>0)); ( use only those values that are greater than 0 )
	//            noiseThreshold=noiseMean+noiseStd*noiseThresholdLevel;  ( apply the threshold value )
	//  end   
	//
	// Parameters:
	// Input: datavector1, a reference to a vector of doubles, containing smoothed signal data(smoothXIs)
	//        originalData: a reference to a vector of doubles, containing original signal data
	//		  noiseThresholdLevel: the threshold level, the number of standard deviations ( pass by value )
	//
	// Output: 
	//         noiseThreshold: the resulting noise threshold value
	//
	//
	//
	// Return code: 0 successful, -1 ( vectors are not same size ), -2( datavector1 sum == 0), -3( datavector2 sum == 0 ), -4( Stotal == 0 )
	//
	// Refer to: /msda/matlab/interfacing/NelsonCode012412/src/getNoiseThreshold.m
	//
	// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 01/30/12
	//
	///////////////////////////////////////////////////////////////////////////
	int getNoiseThreshold(vector<double> & datavector1, vector<double> & datavector2, double noiseThresholdLevel, double & noiseThreshold);

	///////////////////////////////////////////////////////////////////////////
	//
	// moving average smoothing algorithm:
	//
	// This algorithm takes a composite data structure, the XIC, which 
	// contains a set of signals, that are within a narrow mz range and a 
	// specfic scan number
	//
	// Original algorithm developed in Matlab by: Jianqiu Zhang
	// function smoothXIC=movingAvgSmooth(XIC,numScans,numSmoothPoint)
	// smoothXIC=zeros(1,numScans);
	// for scanid=1:numScans
	//     if scanid<=numSmoothPoint;  % smoothing boundary
	//        smoothBegin=1;
	//     else
	//        smoothBegin=scanid-numSmoothPoint;
	//     end
	//     % 
	//     if scanid>numScans-numSmoothPoint;
	//        smoothEnd=numScans;
	//     else
	//        smoothEnd=scanid+numSmoothPoint; 
	//     end
	//     tempSum=0;
	//     for tempid=smoothBegin:smoothEnd  
	//        tempSum=XIC(tempid)+tempSum;
	//     end   
	//    
	//     % automatic handling of boundaries.
	//     smoothXIC(scanid)=tempSum/(smoothEnd-smoothBegin+1);
	// end
	//
	//
	// Parameters:
	// Input: datavector, a reference to a vector of doubles, containing the data to be smoothed
	//        numScans, a double value containing the number of scans to include in the smoothing process
	//		  numSmoothPoint, a double value containing the width of the smoothing kernel
	//
	// Output: 
	//         outputdatavector: the resulting  smoothed vector
	//
	//
	//
	// Return code: 0 successful, -1 ( vectors are not same size ), -2( datavector1 sum == 0), -3( datavector2 sum == 0 ), -4( Stotal == 0 )
	//
	// Refer to: /msda/matlab/interfacing/NelsonCode012412/src/getNoiseThreshold.m
	//
	// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 01/31/12
	//
	///////////////////////////////////////////////////////////////////////////
	int movingAvgSmooth(vector<double> & inputvector, int numSmoothPoints, vector<double> & outputvector);
	
	
	///////////////////////////////////////////////////////////////////////////
	// Name:
	// addNewInterval algorithm:
	//
	// Description:
	// Original algorithm developed in Matlab by: Jianqiu Zhang
	//
	// Parameters:
	// Input: 
	//
	//
	//
	//
	//
	// Output: 
	//         
	//
	// Return code: 0 successful
	//
	// Refer to: /msda/matlab/interfacing/NelsonCode02072012/src/getNoiseThreshold.m
	//
	// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 02/08/12
	//
	///////////////////////////////////////////////////////////////////////////
	int addNewIntervals(			vector<int> & intervalListStart, 
									vector<int> & intervalListEnd,
									vector<double> & XIC,
									double & windowMZ,
									vector<int> & plcGroupID,
									vector<double> & plcMz,
									vector<int> & plcElutionStart,
									vector<int> & plcElutionEnd,
									vector < vector<double> > & plcElutionPeak,
									double & massResolution,
									double & overlapPercentage,
									vector<int> & groupMap,
									double & R2Threshold,
									vector<int> & outputGroupMap);
	
	///////////////////////////////////////////////////////////////////////////
	// Name:  getLCPeakInfo
	// 
	// Calculate statistics about an LC Peak.  
	//
	// Parameters:
	// Input:
	//  1) Data Matrix 
	//     i.e. a 7x7 matrix containing the following data: 
	//     37220.6479772296	37648.0677223798	38075.4874675455	42540.5661674301	70170.2340017539	87391.8931233253	90561.3326282671
	//     59739.1694604254	60274.6156762713	60810.0618921366	69273.3910507242	123220.044366547	160538.145947509	166604.338174930
	//	   33445.5387220645	34076.6389345445	47261.1794838858	83088.0273710790	128355.223433288	160687.906049677	164637.495401330
	//     75068.5056606047	73414.5363775556	85600.2227289044	124734.572700347	185673.015251091	249607.787304643	285226.427586895
	//     34135.1273332247	36327.9120208148	49186.1905321339	83660.3950779063	139803.012647750	184758.766459419	202754.277571757
	//     45458.6558691267	43823.3706442475	42188.0854193090	46315.0287194532	83500.6473008280	111944.527257098	120304.605698007
	//     41851.0904228137	36365.6409656633	47489.7893436464	88636.5820301244	142643.208766860	185240.313170782	198777.642080997
    //
	//  2) mzGridStartID.  This is an index into the mzGrid vector.
	//     i.e. a double value containing the number 65.
	// 
	//  3) mzGrid.  This is a vector containing a linear grid of mz values.
	//     i.e. a 1x79 vector of doubles
	//     750.467311782051	 750.470448833333	750.473585884615	750.476722935897	750.479859987180	750.482997038462	
	//
	// 
	//
	// Output: 
	// Return code: 0 successful
	//
	// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 02/15/12
	//
	///////////////////////////////////////////////////////////////////////////
	int getLCPeakInfo(	vector< vector<double > > & dataMatrix, 
						vector<double>::size_type mzGridStartID, 
						vector<double> & mzGrid,
						vector<double> & elutionProfile,
						double & centerMass,
						bool & maxOnLeft,
						bool & maxOnRight);
	
	///////////////////////////////////////////////////////////////////////////
	// Name:  getInterval
	// 
	// Used to generate a set of intervals of the xic signal sorted by the
	// height of each peak.
	// Each interval should always only have 1 peak.
	//
    // The intervalListStart and intervalListEnd data structures
    // are used to keep track of both the  starting and ending
    //
    //
    //
	///////////////////////////////////////////////////////////////////////////
	int getInterval( 	vector <double> & xic, 
						double threshold, 
						int mininterval, 
						vector<int> & intervalListStart, 
						vector<int> & intervalListEnd );
	
	
	
	
	///////////////////////////////////////////////////////////////////////////
	// splitInterval.m port
	///////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////
	// Name:  splitInterval
	// The purpose of this method is to perform additional cutting of an
	// xic interval.
	// Each interval is split to at most two intervals.
	//
	//
	//
	// Output: 
	// Return code: 0 successful, ..., tbd
	//
	// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 02/20/12
	///////////////////////////////////////////////////////////////////////////
	int splitInterval(		vector <double> & inSmoothXic, 
							vector <double> & inDiffXic,
							vector<int> & inIntervalListStart,
							vector<int> & inIntervalListEnd,
							int & inMinLCLength,
							vector<int> & outIntervalListStart,
							vector<int> & outIntervalListEnd );

	///////////////////////////////////////////////////////////////////////////
	// Name:  getResampledVector
	// The purpose of this method is to perform a linear resampling of 
	// a vector.
	//
	// Input:
	//  1. mzVector:  A 1d vector of doubles
	//  2. intensity: A 1d vector of doubles
	//  3. mzGrid: A 1d vector of doubles
	// 
	//
	// Output:
	//  1. resampledVector: A 1d vector of doubles. ( intensity values ) 
	// Return code: 0 successful, ..., tbd
	//
	// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 02/27/12
	///////////////////////////////////////////////////////////////////////////
	int resampledVector( 	vector<double> & mzVector,
							vector<double> & intensity,
							vector<double> & mzGrid,
							vector<double> & resampledVector);

	///////////////////////////////////////////////////////////////////////////
	//
	// Name:  removeShortIntervals
	// The purpose of this method is to remove intervals that are below a
	// certain minLCLength threshold.
	//
	//
	//
	// Input:
	//  1. 
	//  2. 
	//  3. 
	// Output:
	//  1. 
	//  2. 
	// 	3. 
	// Return Codes:
	//
	// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 03/01/12
	///////////////////////////////////////////////////////////////////////////
	int removeShortIntervals (	vector<int> & inIntervalListStart,
								vector<int> & inIntervalListEnd,
								int & inMinLCLength,
								vector<int> & outIntervalListStart,
								vector<int> & outIntervalListEnd );
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Name:  sortIntervals
	// The purpose of this method is to sort intervals that are below a
	// certain minLCLength threshold.
	//
	//
	//
	// Input:
	//  1. 
	//  2. 
	//  3. 
	// Output:
	//  1. 
	//  2. 
	// Return Codes:
	//
	// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 03/01/12
	///////////////////////////////////////////////////////////////////////////
	int sortIntervals(	vector<int> & inIntervalListStart,
						vector<int> & inIntervalListEnd,
						int & inMinLCLength,
						vector<int> & outIntervalListStart,
						vector<int> & outIntervalListEnd );

	///////////////////////////////////////////////////////////////////////////
	//
	// Name: mergeIntervals
	// The purpose of this method is to remove intervals that meet certain
	// requirements.  Mainly proximity to each other, together with
	// other factors.
	// Full criteria is still to be determined.
	//
	// Input:
	//  1. 
	//  2. 
	//  3. 
	// Output:
	//  1. 
	//  2. 
	// Return Codes:
	//
	// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 03/01/12
	///////////////////////////////////////////////////////////////////////////
	int mergeIntervals(		vector<int> & inIntervalListStart,
							vector<int> & inIntervalListEnd,
							vector<int> & outIntervalListStart,
							vector<int> & outIntervalListEnd );
	

    ///////////////////////////////////////////////////////////////////////////
    //
    // Name: estimateElutionProfile
    //
    // Input:
    //  1. regionData, the piece of the resampled region that corresponds
    //       to this Pre-LC Region of Interest.
    //  2. peakSignal, This signal will be used as the
    //       reference signal, it must have the same number of
    //
    //
    // Output:
    //  1. output: contains the estimatedElutionProfile
    //
    // Return Codes:
    //
    // Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 02/22/13
    ///////////////////////////////////////////////////////////////////////////
    int estimateElutionProfile(	vector< vector<double> > & regionData,
                                vector<double> & peakSignal,
                                vector<double> & output);

    ///////////////////////////////////////////////////////////////////////////
    //
    // Name:
    // The purpose of this method is to generate the set of parallel parameter
    // files together with the set of job scripts that will be used to
    // process the large input data sets according to a work partitioning
    // design.
    //
    // Input:
    //  1.
    //  2.
    //  3.
    // Output:
    //  1.
    //  2.
    // Return Codes:
    //
    // Nelson.Ramirez@utsa.edu, 03/08/12
    ///////////////////////////////////////////////////////////////////////////
    int generateParallelJobFiles(cbi::ParameterLoader & globalParameters );

	
	

    ///////////////////////////////////////////////////////////////////////////
    //
    // Name:
    //
    //
    // Input:
    //  1.
    //  2.
    //  3.
    // Output:
    //  1.
    //  2.
    // Return Codes:
    //
    // Nelson.Ramirez@utsa.edu, 04/16/12
    ///////////////////////////////////////////////////////////////////////////
    int workflowOne(cbi::ParameterLoader & globalParameters );
    int workflowTwo( cbi::ParameterLoader & paraObjGlobal, cbi::ParameterLoader & paraObjLocal);
    int workflowThree( cbi::ParameterLoader & paraObjGlobal, cbi::ParameterLoader & paraObjLocal );




private:
	
	int version;  // a simple counter to count major updates to this class.
	
};

} // namespace msda 
#endif // UTILITIES_H_ 
