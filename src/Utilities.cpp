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
/// Utilities.cpp
/// General utility functions specific to msda software
///  Created on: Aug 8, 2011
///      Author: nelson.ramirez
///
/// Software Developer: Nelson.Ramirez@utsa.edu
///
///////////////////////////////////////////////////////////////////////////////
#include "Utilities.h"
#include <vector>
#include <cmath>
#include <cassert>
#include <string>
#include <sstream>
#include <iomanip>
using namespace std;

namespace msda {

///////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////
Utilities::Utilities() {
	version = 1;
}

///////////////////////////////////////////////////////////////////////////////
// Destructor
// As of version 1, this class does NOT own any resources.
// If any future addition to this class adds a managed resources,
// i.e. the use of new operator, malloc, fopen, .. versus
//      vector, fstream, ... objects that own their resources...
// Then we must implement both a destructor as well as a copy and assignment
// operator. The default shallow copy is almost never correct for classes
// that manage their own memory.
//
//
///////////////////////////////////////////////////////////////////////////////
Utilities::~Utilities() {
	;
}



///////////////////////////////////////////////////////////////////////////
// Version 2 of a peak cluster algorithm,
// Based on matlab code from M.Zhang.10/08/12
//
///////////////////////////////////////////////////////////////////////////
int Utilities::intervalClusterVersion2(      vector<int> & peakIntervalId,
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
                                             vector< vector< vector<int> > > & LCPeakApexIntervalEndIndex) {

    ///////////////////////////////////////////////////////////////////////////
    // create a cbi::SignalProcessor object, needed for the sorting
    // routine.
    ///////////////////////////////////////////////////////////////////////////
    cbi::SignalProcessor spObject;
//    cout<<"Inside Utilities::intervalClusterVersion2"<<endl;

    ///////////////////////////////////////////////////////////////////////////
    //
    // Perform a running summation
    //
    ///////////////////////////////////////////////////////////////////////////
    int numWindows = preLCPeakApexMap.size();
    int numScans = preLCPeakApexMap.at(0).size();

//    cout<<"numWindows="<<numWindows<<endl;
//    cout<<"numScans="<<numScans<<endl;

    for ( int i = 0; i < numScans; i++ ) {
        int scanIndexLow  = spObject.max(0,i-LCPeakApexTolerance); // make sure we are in bounds(left)
        int scanIndexHigh = spObject.min(i+LCPeakApexTolerance,numScans-1); // make sure we are in bounds(right)

        // create a vector to hold the sum of each window
        // the length of the sumScansVector is the total number of mz windows
        vector<double> sumScansVector(numWindows,0);

        // generate the summation in the mz dimension
        for ( int k = 0; k < numWindows; k++ ) {
            for ( int j = scanIndexLow; j <= scanIndexHigh; j++ ) {
                sumScansVector.at(k) = sumScansVector.at(k)+preLCPeakApexIntensity.at(k).at(j);
            }
        }

        // now call interval detection to get intervals along the mz dimension
        int rc = 0;
        vector<int> mzIntervalListStart;
        vector<int> mzIntervalListEnd;

        rc = getInterval( 	sumScansVector,
                            mzIntervalThreshold,
                            mzIntervalMinLength,
                            mzIntervalListStart,
                            mzIntervalListEnd );

        // make sure we have at least 1 interval in the mz dimension
        if( mzIntervalListStart.size() > 0) {
            // cout<<"sumScansVector.size()"<<sumScansVector.size()<<endl;
            // cout<<"mzIntervalListStart.size()="<<mzIntervalListStart.size()<<endl;
            // cout<<"mzIntervalListEnd.size()="<<mzIntervalListEnd.size()<<endl;

            // the number of starting and ending values for intervals must always be the same
            assert(mzIntervalListStart.size() == mzIntervalListEnd.size() );
            // when we find an interval we next want to calculate the summation
            // of the sumScansVector


            // we need a place to store the maximum value of each interval: lcIntervalMzSumMax
            // one storage location is reserved to keep track of the peaks in the
            // mz
            vector<double> lcIntervalMzSumMax(mzIntervalListStart.size(),0);    // Maps to intensity
            vector<int> lcIntervalMzSumMaxIndex(mzIntervalListStart.size(),0);  // Mz location ( maps to mz window# of XIC data structure)
            vector<int> lcIntervalScanLocation(mzIntervalListStart.size(),0);   // Scan location ( maps to scan dimension of XIC data structure)
            vector<int> lcIntervalMzIntervalStart(mzIntervalListStart.size(),0);   //Mz Interval Start( maps to mz window# of XIC data structure)
            vector<int> lcIntervalMzIntervalEnd(mzIntervalListStart.size(),0);     //Mz Interval End( maps to mz window# of XIC data structure)


            // go through each interval found in the mz dimension
            for ( int j = 0; j < mzIntervalListStart.size(); j++ ) {

                // create a temporary copy of the data within an interval
                vector<double> intervalDataBuffer;
                for (int k = mzIntervalListStart.at(j); k <= mzIntervalListEnd.at(j); k++ ){
                    intervalDataBuffer.push_back(sumScansVector.at(k));
                }

                // find the maximum value within the interval
                // the findMaxValue function will return
                vector<double>::size_type peakIntensityIndex;
                double peakIntensityValue;

                int rc = spObject.findMaxValue(intervalDataBuffer,peakIntensityIndex,peakIntensityValue);

                //cout<<"peakIntensityValue="<<peakIntensityValue<<endl;

                // save the peak intensity value as well as its index
                lcIntervalMzSumMax.at(j) = peakIntensityValue; // save the summed intensity value at the peak
                lcIntervalMzSumMaxIndex.at(j) = peakIntensityIndex+mzIntervalListStart.at(j); // save the mz index of peak
                lcIntervalScanLocation.at(j) = i; // save the scan number
                lcIntervalMzIntervalStart.at(j) = mzIntervalListStart.at(j); // save the mz value of the start of the interval along mz dimension
                lcIntervalMzIntervalEnd.at(j) = mzIntervalListEnd.at(j); // save the mz value of the end of the interval along mz dimension


                // save the mz interval start and end range
                LCPeakApexIntervalStartIndex.at(lcIntervalMzSumMaxIndex.at(j)).at(i) = vector<int>(1,lcIntervalMzIntervalStart.at(j));
                LCPeakApexIntervalEndIndex.at(lcIntervalMzSumMaxIndex.at(j)).at(i) = vector<int>(1,lcIntervalMzIntervalEnd.at(j));


            } // end loop through mz intervals found

            ///////////////////////////////////////////////////////////////////////
            //
            // Output the mz interval peak information
            //
            ///////////////////////////////////////////////////////////////////////
//            cout<<"i="<<i<<"lcIntervalMzSumMax.size()="<<lcIntervalMzSumMax.size()<<endl;
//            cout<<"i="<<i<<"lcIntervalMzSumMaxIndex.size()="<<lcIntervalMzSumMaxIndex.size()<<endl;
//            cout<<"i="<<i<<"lcIntervalScanLocation.size()="<<lcIntervalScanLocation.size()<<endl;
//            cout<<"i="<<i<<"lcIntervalMzIntervalStart.size()="<<lcIntervalMzIntervalStart.size()<<endl;
//            cout<<"i="<<i<<"lcIntervalMzIntervalEnd.size()="<<lcIntervalMzIntervalEnd.size()<<endl;

            // save result to the global data structure
            // go through each interval found in the mz dimension
            for ( int j = 0; j < lcIntervalMzSumMax.size(); j++ ) {
                LCPeakApexMap.at(lcIntervalMzSumMaxIndex.at(j)).at(lcIntervalScanLocation.at(j)) = 1;
                LCPeakApexIntensity.at(lcIntervalMzSumMaxIndex.at(j)).at(lcIntervalScanLocation.at(j)) = lcIntervalMzSumMax.at(j);
                // we also want to save the mz start and end ranges...
            }


        } // end if check: make sure there is at least one interval in the mz dimension found
    } // end loop through all the mz dimension signals ( each signal has a single rt value, and varies along the mz dimension )



    // create a temporary buffer to hold the AND operation results
    vector< vector<int> > LCPeakApexMerged(numWindows,vector<int>(numScans,0));

    // Perform a logical AND operation using the data items within the
    // preLCPeakApexMap(Retention time(scan) dimension peaks) and the LCPeakApexMap(mz dimension peaks)
    assert(preLCPeakApexMap.size() == LCPeakApexMap.size() );

    for ( int j = 0; j < preLCPeakApexMap.size(); j++ ) {
        for( int k = 0; k < preLCPeakApexMap.at(j).size(); k++ ) {
            bool plcData = preLCPeakApexMap.at(j).at(k);  // convert from int to boolean
            bool lcData = LCPeakApexMap.at(j).at(k);      // convert from int to boolean
            LCPeakApexMerged.at(j).at(k) = (plcData && lcData);

            // save the coordinates of the peak apex
            // j --> maps to apexWindowId
            // k --> maps to apexScanId
            if ( LCPeakApexMerged.at(j).at(k) == 1) {
                apexWindowId.push_back(j);
                apexScanId.push_back(k);
            } // end check so we can save the peak location ( windowid and scanid )
        } // end loop through all scan#'s( rt dimension )
    } // end loop through all mz windows ( mz dimension )


    // Output merged mz and rt peaks
    LCPeakApexMap = LCPeakApexMerged;


    return 0;
}



///////////////////////////////////////////////////////////////////////////////
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
// Input:
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
//  -1,-2,-3 ( input data size mismatch )
//  -4 ( error in the sorting process )
//
// Developed by:
// Michelle.Zhang@utsa.edu, Nelson.Ramirez@utsa.edu 10/3/2012.
// Copyright CBI,UTSA, All rights reserved.  2012.
//
///////////////////////////////////////////////////////////////////////////////
int Utilities::intervalCluster(     vector<int> & peakIntervalId,
                                    vector<double> & peakIntervalMz,
                                    vector<double> & peakIntervalRt,
                                    vector<double> & peakIntervalIntensity,
                                    vector< vector<int> >  & groupsIdtoIntervalIdVectorsResult,
                                    vector< vector<double> > & groupsIdtoIntervalMzVectorsResult,
                                    vector< vector<double> > & groupsIdtoIntervalRtVectorsResult,
                                    vector< vector<double> > & groupsIdtoIntervalIntensityVectorsResult,
                                    vector<int> &  intervalIdtoGroupIdMappingResult ) {


    ///////////////////////////////////////////////////////////////////////////
    // create a cbi::SignalProcessor object, needed for the sorting
    // routine.
    ///////////////////////////////////////////////////////////////////////////
    cbi::SignalProcessor spObject;

    // store return codes from function calls
    int rc = 0;

    // validity checking on the input vector sizes
    if ( peakIntervalMz.size() != peakIntervalRt.size() ) {
        return -1;
    }

    // validity checking on the input vector sizes
    if ( peakIntervalMz.size() != peakIntervalId.size() ) {
        return -2;
    }

    // validity checking on the input vector sizes
    if ( peakIntervalMz.size() != peakIntervalIntensity.size() ) {
        return -3;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Step 1. mz*rt ( reduce the dimensionality of the data )
    ///////////////////////////////////////////////////////////////////////////
    vector<double> mzRtMultiplyBuffer;

    vector<int> peakIntervalIdTemp;
    vector<double> peakIntervalMzTemp;
    vector<double> peakIntervalRtTemp;
    vector<double> peakIntervalIntensityTemp;
    vector<int> sortIndex;  // buffer with linear indices

    int count = 0;
    ///////////////////////////////////////////////////////////////////////////
    // Unusable data removal
    ///////////////////////////////////////////////////////////////////////////
    for( vector< double >::size_type i = 0; i <  peakIntervalMz.size(); i++ ) {
        //cout<<setprecision(10)<<"index="<<peakIntervalId.at(i)<<",mz="<<peakIntervalMz.at(i)<<",rt="<<peakIntervalRt.at(i)<<"intensity="<<peakIntervalIntensity.at(i)<<endl;
        if ( (peakIntervalMz.at(i) > 0) && (peakIntervalRt.at(i) > 0) ) {
            mzRtMultiplyBuffer.push_back(peakIntervalMz.at(i)* peakIntervalRt.at(i));

            peakIntervalIdTemp.push_back(peakIntervalId.at(i));
            peakIntervalMzTemp.push_back(peakIntervalMz.at(i));
            peakIntervalRtTemp.push_back(peakIntervalRt.at(i));
            peakIntervalIntensityTemp.push_back(peakIntervalIntensity.at(i));


            sortIndex.push_back(count);
            count++;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // Step 2. sort the data in ascending order
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // Sort without changing the input data
    // Order: ascending
    ///////////////////////////////////////////////////////////////////////////
    rc = spObject.sortGetIndex(mzRtMultiplyBuffer, mzRtMultiplyBuffer, sortIndex, true);

    if ( rc != 0 ){
        return -4; // error in the sortGetIndex function
    }

    ///////////////////////////////////////////////////////////////////////////
    // Step 3. cut the sorted data into bins
    ///////////////////////////////////////////////////////////////////////////
    vector<int> groupCutPointsLevel1(mzRtMultiplyBuffer.size(),0); // contains the cutpoints for each level 1 partitioning
    vector<int> groupCutPointsLevel2(mzRtMultiplyBuffer.size(),0); // contains the cutpoints for each level 2 partitioning

    double runningAverageMzRt = 0;
    int runningAverageMzRtCount = 1;

    ///////////////////////////////////////////////////////////////////////////
    //
    // Multi-stage lc candidate refinement algorithm:
    //   The index gets us back to the smoothed XIC data of the
    //   intervals.
    //
    ///////////////////////////////////////////////////////////////////////////
    int groupIdStage1 = 1;  // First level rough binning ( by mz*rt )
    int groupIdStage2 = 1;  // Second level binning by ( rt )


    // Go through the interval data by stepping through the mz*rt column,
    // which is sorted in ascending order.
    for (vector< double >::size_type k = 0; k < mzRtMultiplyBuffer.size()-1; k++ ){
        // Update mz*rt running average
        runningAverageMzRt = (runningAverageMzRt + mzRtMultiplyBuffer.at(sortIndex.at(k)))/runningAverageMzRtCount;
        runningAverageMzRtCount++;

        // Apply cut criteria
        if ( mzRtMultiplyBuffer.at(sortIndex.at(k+1)) > (mzRtMultiplyBuffer.at(sortIndex.at(k))+0.05*runningAverageMzRt)  ) {
            runningAverageMzRtCount=1;
            groupIdStage1++;
        }
        groupCutPointsLevel1.at(k+1) = groupIdStage1;


        ///////////////////////////////////////////////////////////////////////////
        // Step 4. perform a second pass using rt information
        //         to refine the classification.
        ///////////////////////////////////////////////////////////////////////////
        // Apply cut criteria for rt based segmentation
        if (  abs(peakIntervalRtTemp.at(sortIndex.at(k+1)) - peakIntervalRtTemp.at(sortIndex.at(k))) > 0.001  ) {
            groupIdStage2++;
        }
        groupCutPointsLevel2.at(k+1) = groupIdStage2;

    }

    ///////////////////////////////////////////////////////////////////////////
    //
    // Debug output
    //
    ///////////////////////////////////////////////////////////////////////////
//    for( vector< double >::size_type i = 0; i <  mzRtMultiplyBuffer.size(); i++ ) {
//        cout<<setprecision(10)<<"groupIdStage1="<<groupCutPointsLevel1.at(i)<<"groupIdStage2="<<groupCutPointsLevel2.at(i)<<",intervalIndex="<<peakIntervalIdTemp.at(sortIndex.at(i))<<",mz="<<peakIntervalMzTemp.at(sortIndex.at(i))<<",rt="<<peakIntervalRtTemp.at(sortIndex.at(i))<<",intensity="<<peakIntervalIntensityTemp.at(sortIndex.at(i))<<endl;
//    }


    ///////////////////////////////////////////////////////////////////////////
    //
    // At this point the peak of each interval have been organized
    // into
    // groupIdStage1=34groupIdStage2=430,intervalIndex=168,mz=1988.931708,rt=2061.919922intensity=3537.733123
    // groupIdStage1=34groupIdStage2=430,intervalIndex=172,mz=1988.951644,rt=2061.919922intensity=4084.022902
    // groupIdStage1=34groupIdStage2=430,intervalIndex=176,mz=1988.971581,rt=2061.919922intensity=6160.655486
    //
    ///////////////////////////////////////////////////////////////////////////
    vector< vector<int> >  groupsIdtoIntervalIdVectors;  // group # --> interval 1,2,3,5 ( groupsId to interval Id )
    vector< vector<double> >  groupsIdtoIntervalMzVectors;   // group # --> mz1,mz2,mz3,mz5 ( groupsId to interval mz value at peak )
    vector< vector<double> >  groupsIdtoIntervalRtVectors;   // group # --> rt1,rt2,rt3,rt5 ( groupsId to interval rt value at peak )
    vector< vector<double> >  groupsIdtoIntervalIntensityVectors; // group # --> intensity1, intensity2, intensity3, intensity5 ( groupsId to interval intensity value at peak )
    // We should be able to go back from interval# to the group# associated with this interval
    vector< int >  intervalIdtoGroupIdMapping;  // interval# --> group# mapping ( mappint interval ids to group ids )


    if ( groupCutPointsLevel2.size() == 0 )  {
        return -5; // no cutpoints found...
    }

    ///////////////////////////////////////////////////////////////////////////
    //
    //
    // Intermediate vectors to store the data for each group id
    //
    ///////////////////////////////////////////////////////////////////////////
    vector<int> groupsIdtoIntervalId;
    vector<double> groupsIdtoIntervalMz;
    vector<double> groupsIdtoIntervalRt;
    vector<double> groupsIdtoIntervalIntensity;

    // initialize with first data point
    groupsIdtoIntervalId.push_back(peakIntervalIdTemp.at(sortIndex.at(0)));
    groupsIdtoIntervalMz.push_back(peakIntervalMzTemp.at(sortIndex.at(0)));
    groupsIdtoIntervalRt.push_back(peakIntervalRtTemp.at(sortIndex.at(0)));
    groupsIdtoIntervalIntensity.push_back(peakIntervalIntensityTemp.at(sortIndex.at(0)));

    ///////////////////////////////////////////////////////////////////////
    // obtain a map from the interval id to the group#
    ///////////////////////////////////////////////////////////////////////
    intervalIdtoGroupIdMapping.push_back(peakIntervalIdTemp.at(sortIndex.at(0)));


    for( vector< double >::size_type i = 1; i <  mzRtMultiplyBuffer.size(); i++ ) {
        //cout<<setprecision(10)<<"groupIdStage1="<<groupCutPointsLevel1.at(i)<<"groupIdStage2="<<groupCutPointsLevel2.at(i)<<",intervalIndex="<<peakIntervalIdTemp.at(sortIndex.at(i))<<",mz="<<peakIntervalMzTemp.at(sortIndex.at(i))<<",rt="<<peakIntervalRtTemp.at(sortIndex.at(i))<<",intensity="<<peakIntervalIntensityTemp.at(sortIndex.at(i))<<endl;

        ///////////////////////////////////////////////////////////////////////
        // obtain a map from the interval id to the group#
        ///////////////////////////////////////////////////////////////////////
       intervalIdtoGroupIdMapping.push_back(peakIntervalIdTemp.at(sortIndex.at(i)));

        if (  groupCutPointsLevel2.at(i) == groupCutPointsLevel2.at(i-1) ) {
            // append to current set of vectors
            groupsIdtoIntervalId.push_back(peakIntervalIdTemp.at(sortIndex.at(i)));
            groupsIdtoIntervalMz.push_back(peakIntervalMzTemp.at(sortIndex.at(i)));
            groupsIdtoIntervalRt.push_back(peakIntervalRtTemp.at(sortIndex.at(i)));
            groupsIdtoIntervalIntensity.push_back(peakIntervalIntensityTemp.at(sortIndex.at(i)));
        }
        else {

            // update outer data structure with vector containing
            groupsIdtoIntervalIdVectors.push_back(groupsIdtoIntervalId);
            groupsIdtoIntervalMzVectors.push_back(groupsIdtoIntervalMz);
            groupsIdtoIntervalRtVectors.push_back(groupsIdtoIntervalRt);
            groupsIdtoIntervalIntensityVectors.push_back(groupsIdtoIntervalIntensity);

            // clear the current intermediate vectors
            groupsIdtoIntervalId.clear();
            groupsIdtoIntervalMz.clear();
            groupsIdtoIntervalRt.clear();
            groupsIdtoIntervalIntensity.clear();

            // re-initialize with the first data item
            groupsIdtoIntervalId.push_back(peakIntervalIdTemp.at(sortIndex.at(i)));
            groupsIdtoIntervalMz.push_back(peakIntervalMzTemp.at(sortIndex.at(i)));
            groupsIdtoIntervalRt.push_back(peakIntervalRtTemp.at(sortIndex.at(i)));
            groupsIdtoIntervalIntensity.push_back(peakIntervalIntensityTemp.at(sortIndex.at(i)));

        }  // end else

    } // end loop through all interval peaks

//    cout<<"groupsIdtoIntervalIdVectors.size()="<<groupsIdtoIntervalIdVectors.size()<<endl;  // group # --> interval 1,2,3,5 ( groupsId to interval Id )
//    cout<<"groupsIdtoIntervalMzVectors.size()="<<groupsIdtoIntervalMzVectors.size()<<endl;   // group # --> mz1,mz2,mz3,mz5 ( groupsId to interval mz value at peak )
//    cout<<"groupsIdtoIntervalRtVectors.size()="<<groupsIdtoIntervalRtVectors.size()<<endl;   // group # --> rt1,rt2,rt3,rt5 ( groupsId to interval rt value at peak )
//    cout<<"groupsIdtoIntervalIntensityVectors.size()="<<groupsIdtoIntervalIntensityVectors.size()<<endl; // group # --> intensity1, intensity2, intensity3, intensity5 ( groupsId to interval intensity value at peak )
//    // We should be able to go back from interval# to the group# associated with this interval
//    cout<<"intervalIdtoGroupIdMapping.size()="<<intervalIdtoGroupIdMapping.size()<<endl;  // interval# --> group# mapping ( mappint interval ids to group ids )



    groupsIdtoIntervalIdVectorsResult = groupsIdtoIntervalIdVectors;
    groupsIdtoIntervalMzVectorsResult = groupsIdtoIntervalMzVectors;
    groupsIdtoIntervalRtVectorsResult = groupsIdtoIntervalRtVectors;
    groupsIdtoIntervalIntensityVectorsResult = groupsIdtoIntervalIntensityVectors;
    intervalIdtoGroupIdMappingResult = intervalIdtoGroupIdMapping;




    return 0;
} // end clustering function



///////////////////////////////////////////////////////////////////////////////
//
//
// KL distance algorithm:
// Original algorithm developed in Matlab by: Jianqiu Zhang
// function  log_KL_Value=KL_calculate(Vector01,Vector02)
//%%%%%%%%%%%% calculate KL based on add small values to zeros position
//ID_zeros01=find(Vector01==0);  % ??? This find routine needs to be implemented
//ID_zeros02=find(Vector02==0);
//Vector01v1=Vector01;
//Vector02v1=Vector02;
//if ~isempty(ID_zeros01)
//    Vector01v1(ID_zeros01)=max(Vector01)/10000000;
//end
//if ~isempty(ID_zeros02)
//    Vector02v1(ID_zeros02)=max(Vector02)/10000000;
//end
//Vector01v2=Vector01v1./sum(Vector01v1);
//Vector02v2=Vector02v1./sum(Vector02v1);
//
//
//log_KL_Value=log(sum(Vector01v2.*log(Vector01v2./Vector02v2)));
//
//end
//
//
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
//
///////////////////////////////////////////////////////////////////////////////
int Utilities::logKLDistance(   const vector<double> & inputA,
                                const vector<double> & inputB,
                                double & result,
                                bool & equal) {



    // create a cbi::SignalProcessor object
    cbi::SignalProcessor spObject;

    // store return codes from function calls
    int rc = 0;

    //
    //
    //
    // key validity items
    // log(0) =-inf, when both inputA and inputB are the same, this condition
    //               needs to be handled.
    // If any of the data items in inputA or inputB are zero,
    // they need to be converted to a small value
    //
    // if the length of inputA is zero
    //
    // if the length of inputB is zero
    //
    // input vectors must have the same length
    //
    // if the inputs contain negative values, we must
    // return an error code.  The KL distance
    // is not valid if the input data contains
    // negative values.
    //
    //
    //

    if ( inputA.size() < 1) {
       return -1;
    }

    if ( inputB.size() < 1 ){
       return -2;
    }

    if ( inputA.size() != inputB.size() ){
        return -3;
    }


    // we need to make absolutely sure we do not modify the inputs
    // make temporary variables
    vector<double> tempA = inputA;
    vector<double> tempB = inputB;


    // check to see if the inputs contain negative data
    for ( vector<double>::iterator ii = tempA.begin(); ii != tempA.end(); ii++ ){
        if ( *ii < 0 ) {
            return -4;  // indicate inputA contains negative data items
        }
    }

    // fix zero items by mapping to a very small but non-zero value
    for ( vector<double>::iterator ii = tempB.begin(); ii != tempB.end(); ii++ ){
        if ( *ii < 0 ) {
            return -5; // indicate inputB contains negative data items
        }
    }

    // check to see if the vectors contain the same data
    bool inputEqualityCheck = false; // initialize
    // check if the two input vectors are equal
    spObject.isEqual(tempA,tempB,0.0000001,inputEqualityCheck);

    //
    // set equal output variable
    //
    if ( inputEqualityCheck == true){
        result = -10000000;
        equal = true;
        return 0;
    }

    // variables to hold maximum values of input vectors
    double maxA=0;
    double maxB=0;
    double divisor = 10000000;
    double sumA=0; // contains the sum of the inputA vector
    double sumB=0; // contains the sum of the inputB vector

    spObject.findMaxValue(tempA, maxA );

    spObject.findMaxValue(tempB, maxB );

    // fix zero items by mapping to a very small but non-zero value
    // since we need to perform division by the sum of both of these
    // vectors
    for ( vector<double>::iterator ii = tempA.begin(); ii != tempA.end(); ii++ ){
        if ( *ii == 0 ) {
            *ii = maxA/divisor;
        }
    }

    // fix zero items by mapping to a very small but non-zero value
    for ( vector<double>::iterator ii = tempB.begin(); ii != tempB.end(); ii++ ){
        if ( *ii == 0 ) {
            *ii = maxB/divisor;
        }
    }

    //
    // Find the sum of the input vectors
    //
    rc = spObject.sum(tempA,sumA);
    rc = spObject.sum(tempB,sumB);

    // make sure sumA is never 0
    // make sure sumB is never 0

    if (sumA <= 0 ) {
        // prevent a potential divide-by-zero exception
        return -6;  // this should never happen at this point
    }

    if ( sumB <= 0) {
        // prevent a potential divide-by-zero exception
        return -7; // this should never happen at this point
    }

    //
    // perform an element by element division
    //
    for ( vector<double>::iterator ii = tempA.begin(); ii != tempA.end(); ii++ ){
       *ii = *ii / sumA;
    }

    //
    // perform an element by element division
    //
    for ( vector<double>::iterator ii = tempB.begin(); ii != tempB.end(); ii++ ){
        *ii = *ii / sumB;
    }

    //
    //
    // calculate the log(KL) value
    //
    //
    vector<double> tempResult;

    //
    //
    // Perform an element by element division operation
    // of the items in tempA ./ tempB
    //
    //
    // get iterators and point them to the start of each input
    vector<double>::iterator ia = tempA.begin();
    vector<double>::iterator ib = tempB.begin();

    // if the vectors are the same size, then we need to check
    // their contents.
    // Note: We've already checked to make sure that both
    // vectors have the same size.
    //
    for ( ia, ib; ia != tempA.end(); ia++,ib++) {
        double tempValue = 0;
        tempValue = (*ia) / (*ib);
        tempValue = log(tempValue);
        tempValue = (*ia) * tempValue;
        tempResult.push_back(tempValue);
    }

    //
    // calculate the running sum of the tempResult
    //
    double sumTempValue = 0;
    rc = spObject.sum(tempResult,sumTempValue);

    //
    // Generate final results
    //
    result = log(sumTempValue); // return result variable
    equal = false;  // indicate vectors were not equal
    return 0; // return successful

}





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
int Utilities::calcIsoPatternAvg(  const double & mass, const int & maxiso, vector<double> & averageIsotopePattern){

    // create a cbi::SignalProcessor object
    cbi::MathUtils muObject;

    cout<<"Testing calcIsoPatternAvg"<<endl;

    //carbon13ratio=1.1078;
    //oxy17ratio=0.0372;
    //oxy18ratio=0.2004;
    //nitrogen15ratio=0.3663;
    //hydrogen2ratio=0.015574;
    //iso1ratio=[1.1078 0.015574 0.3663 0.0372 0.750]/100;
    //%iso2ratio=[0       0        0     0.2004 4.215];
    //%iso3ratio=[0       0        0     0       0.017];
    //atomcount=[4.9384 7.7583 1.3577 1.4773 0.0417];
    //averageamino=111.1254;


    vector<double> iso1ratio(5,0);
    iso1ratio.at(0) = 1.1078/100;
    iso1ratio.at(1) = 0.015574/100;
    iso1ratio.at(2) = 0.3663/100;
    iso1ratio.at(3) = 0.0372/100;
    iso1ratio.at(4) = 0.750/100;

    vector<double> atomcount(5,0);
    atomcount.at(0) = 4.9384;
    atomcount.at(1) = 7.7583;
    atomcount.at(2) = 1.3577;
    atomcount.at(3) = 1.4773;
    atomcount.at(4) = 0.0417;

    double averageAmino = 111.1254; // Average amino acid mass

    double atomnumber = round((mass/averageAmino)* atomcount.at(0));


    double carbonPoissonAvg= atomnumber * iso1ratio.at(0);  // lambda parameter for poisspdf


    // output vector
    vector<double> isoPatternAverageTemp;

    // call the poisspdf function

    // make sure to check that maxiso-1 is within the limits
    // of the poisspdf function
    int rc = 0;
    int min = 0;
    int max = maxiso -1;
    double lambda = carbonPoissonAvg;
    rc = muObject.poissonpdf(min, max, lambda, isoPatternAverageTemp);

    cout<<"rc = "<<endl;

    // copy internal buffer to output
    averageIsotopePattern = isoPatternAverageTemp;
    return 0; // return successful
}



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
//    1) labelingID, integer
//       - identifies the type of labeling method used
//
//    2) peptideCandidateAdjustedElutionProfiles
//       - In Matlab, each column contains the adjusted elution profile for one of the
//         lc candidates that are part of a peptide candidate.
//            % We got area factors which should be applied to each elution profile.
//            % We want to estimate the monoistopic mass
//            % The taller the peak, the higher the confidence
//
//       - Note how for each peptideCandidateListItem, there are a set of
//         LC Candidates associated with each peptideCandidateListItem.
//
//       - Each column of the adjustedElutionProfiles and elutionProfiles
//         data structure contains the summarized elution profile for one
//         of the LC Candidates associated with this peptide candidate.
//
//
//            numOfCandidates=length(peptideCandidateListItem.AssociatedLCCandidateIDs);
//            adjustedElutionProfiles=elutionProfiles;
//            for nid=1:numOfCandidates
//                  adjustedElutionProfiles(:,nid)=elutionProfiles(:,nid)*AreaFactor(nid);
//            end
//
//      - In C++, we use the vector container to hold generalized
//        vector data.
//
//          - We use a vector of vectors to hold the data
//            for this data structure.
//
//      - Therefore, the peptideCandidateAdjustedElutionProfile
//        contains the following:
//
//
//
//          peptideCandidateAdjustedElutionProfile.at(0) --> Summarized elution profile of the first LC Candidate
//              associated with this peptide candidate
//              - Each LC Candidate has a set of elution signals associated with it.
//
//                For example, LC Candidate # 0
//                     lcc=0,signal=0 = [scan#=400],[(mz,intensity),(mz,intensity),...(mz,intensity)]
//                     lcc=0,signal=1 = [scan#=401],[(mz,intensity),(mz,intensity),...(mz,intensity)]
//                     lcc=0,signal=2 = [scan#=402],[(mz,intensity),(mz,intensity),...(mz,intensity)]
//                     lcc=0,signal=3 = [scan#=403],[(mz,intensity),(mz,intensity),...(mz,intensity)]
//                     lcc=0,signal=4 = [scan#=404],[(mz,intensity),(mz,intensity),...(mz,intensity)]
//                     lcc=0,signal=5 = [scan#=405],[(mz,intensity),(mz,intensity),...(mz,intensity)]
//
//                You can view this as a little box of data, the intensity is the z axis.
//                The mz is the x axis, and the scan# is the y axis.
//
//                For an LC Candidate, the mz values should all be very close to each other, and the
//                intensities should form a 3d peak.
//
//                Over this fundamental data structure, we can define a set of summarized properties
//                that can be used to work with 1d signals instead of using the entire 3d data structure.
//
//                The first signal that can be defined is the elution profile of an LC Candidate
//
//                This is a view of this data by projecting the data such that the horizontal axis
//                is the scan dimension, and each y axis data value is the sum of all the
//                intensity values for that signal.
//
//                Therefore, the elution profile is a 2d signal.
//
//
//               Cumulative
//               Intensity
//                ^
//                |            +
//                |        +      +
//                |     +            +
//                |   +                 +
//                | +                      +
//                |                          +
//                -------------------------------> Scan#
//
//                Since the mz dimension information is not displayed in this view, we can associate
//                a statistic derived from the mz values of each of the signals belonging to this
//                lc candidate, such as the median mz, mean mz, ...
//
//
//
//                For example, LC Candidate # 1
//                     lcc=1,signal=0 = [scan#],[(mz,intensity),(mz,intensity),...(mz,intensity)]
//                     lcc=1,signal=1 = [scan#],[(mz,intensity),(mz,intensity),...(mz,intensity)]
//                     lcc=1,signal=2 = [scan#],[(mz,intensity),(mz,intensity),...(mz,intensity)]
//
//
//
//
//
//
// Output:
//    1)
//
//
// Return code:
//
// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 06/22/12
// Original Matlab function:
//   [peptideFitnessScore]=checkModelfitness(peptideCandidateListItem, labelingMethod,isoPatternNatural)
//    by M.Zhang
///////////////////////////////////////////////////////////////////////////
int Utilities::checkModelFitness( const int labelingId,
                                  const vector< vector<double> > & peptideCandidateAdjustElutionProfiles,
                                  const vector<double> & massTemplate,
                                  const vector<double> & labelingMethodMassShiftVector,
                                  const vector<double> & isoPatternNatural){

    // create a cbi::SignalProcessor object
    cbi::SignalProcessor spObject;


    // each labeling method has an assocated mass template




    // the sum of the intensities of all the adjusted elution profiles in a peptide candidate
    // list item
    double sumIntensities  = 0;

    //
    //
    //
    //
    // loop through each adjusted elution profile.
    // Clearly define what is an adjusted elution profile:
    //
    //
    //
    //
    //

    int sumRows(vector<vector<double> > & data , vector<double> &result); // 2d summation( row-wise, dimension = 2 )


    for (vector<double>::size_type i = 0; i <= peptideCandidateAdjustElutionProfiles.size(); i++ ) {



    }

    //
    // to get the length of the mass shift vector, just use
    // labelingMethodMassShiftVector.length()
    //



    return 0;
}





///////////////////////////////////////////////////////////////////////////////
//
//
// Vector differential calculation algorithm:
// Original algorithm developed in Matlab by: Jianqiu Zhang
// "% Created by Jianqiu Zhang for calculating the differential of a vector XIC
//% with length ll.
//function diff=getDiff(XIC,ll,direction)
//diff=zeros(ll,1);
//if direction==-1
//    for i=1:ll
//        if i==1
//           diff(i)=XIC(i);
//        else
//           diff(i)=XIC(i)-XIC(i-1);
//        end   
//    end 
//else   
//    for i=1:ll
//        if i==ll
//           diff(i)=XIC(i);
//        else
//           diff(i)=XIC(i+1)-XIC(i);
//        end   
//    end 
//end"  
//
// Parameters:
// Input: xic, a reference to a vector of doubles
//        direction, an integer indicating whether left(-1) or right(1) differencing
// Output: 
// A vector of doubles that is 1 element longer than the input vector.( bug fix, 10/17/12, avoid leaving out information )
//
// Return code: 0 successful, -1 (xic length error), -2(direction error)
//
// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 01/25/12
//
///////////////////////////////////////////////////////////////////////////////
int Utilities::getDiff(vector<double> & xic, int direction, vector<double> & diffxic) {
	
	///////////////////////////////////////////////////////////////////////////
	// check for invalid length
	///////////////////////////////////////////////////////////////////////////
	if ( xic.size() < 2) {
		return -1;
	}
	
	///////////////////////////////////////////////////////////////////////////
	// check for invalid direction
	///////////////////////////////////////////////////////////////////////////
	if( direction != -1 && direction != 1) {
		return -2;
	}
	
	///////////////////////////////////////////////////////////////////////////
	// Allocate a working buffer, of the same size as the input
	// xic vector, and initialize it to zero.
	// Since tempdiff is a vector, if any exception were to take place
	// after it is fully constructed, we can be sure that its 
	// memory will be de-allocated by its destructor.
    ///////////////////////////////////////////////////////////////////////////

    vector<double> tempdiff(xic.size(),0); // initialize a temporary buffer to hold working results
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Evaluate the left difference
	//
	///////////////////////////////////////////////////////////////////////////
	if ( direction == -1 ) {
        for ( vector< double>::size_type i = 0; i < xic.size(); i++ ) {
            if ( i == 0) {  // handle first item
				tempdiff.at(i) = xic.at(i);
			}
			else {
				tempdiff.at(i) = xic.at(i) - xic.at(i-1);
			}
		}
	}
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Evaluate the right difference
	//
	///////////////////////////////////////////////////////////////////////////
	else if ( direction == 1) {
        for ( vector< double>::size_type i = 0; i < xic.size(); i++ ) {
			if ( i == xic.size()-1) {
				tempdiff.at(i) = xic.at(i);
			}
			else {
				tempdiff.at(i) = xic.at(i+1) - xic.at(i);
			}
		}
	}
	else {
		return -2;
	}

	///////////////////////////////////////////////////////////////////////////
	//
	//  Create a deep copy of the tempdiff object
	//  This will copy the contents of the tempdiff vector into the memory of
	//  diffxic vector.  Then, as soon as this function goes out of scope,
	//  the destructor of the tempdiff object will get called and all
	//  the data that was allocated for the tempdiff vector will be freed
	//  from the heap.  However, since we copied its contents into another
	//  vector that lives outside of this function we made sure to keep the
	//  results.
	//  The copy constructor of diffxic will get called, allocating all 
	//  required memory for performing the data copy automatically.
	//
	///////////////////////////////////////////////////////////////////////////
	diffxic = tempdiff;   // copy contents of tempdiff into diffxic.
	
	return 0; // return successfully
	
}  // end getdiff method


///////////////////////////////////////////////////////////////////////////////
//
//
// Vector differential calculation algorithm:
// Original algorithm developed in Matlab by: Jianqiu Zhang
// "% Created by Jianqiu Zhang for calculating the differential of a vector XIC
//% with length ll.
//function diff=getDiff(XIC,ll,direction)
//diff=zeros(ll,1);
//if direction==-1
//    for i=1:ll
//        if i==1
//           diff(i)=XIC(i);
//        else
//           diff(i)=XIC(i)-XIC(i-1);
//        end   
//    end 
//else   
//    for i=1:ll
//        if i==ll
//           diff(i)=XIC(i);
//        else
//           diff(i)=XIC(i+1)-XIC(i);
//        end   
//    end 
//end"  
//
//
//
//
//
// Parameters:
// Input: xic, a reference to a vector of ints
//        direction, an integer indicating whether left(-1) or right(1) differencing
// Output: 
//        diffxic, the output vector, where the first item is a copy of the
//          source vector.  This first item contains whatever the
//          input vector contains as its first data item.
//          Care must be taken not to utilize this first item in
//
//
// Return code: 0 successful, -1 (xic length error), -2(direction error)
//
// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 01/25/12
//
///////////////////////////////////////////////////////////////////////////////
int Utilities::getDiff(vector<int> & xic, int direction, vector<int> & diffxic) {
	
	///////////////////////////////////////////////////////////////////////////
	// check for invalid length
	///////////////////////////////////////////////////////////////////////////
	if ( xic.size() < 2) {
		return -1;
	}
	
	///////////////////////////////////////////////////////////////////////////
	// check for invalid direction
	///////////////////////////////////////////////////////////////////////////
	if( direction != -1 && direction != 1) {
		return -2;
	}
	
	///////////////////////////////////////////////////////////////////////////
	// Allocate a working buffer, of the same size as the input
	// xic vector, and initialize it to zero.
	// Since tempdiff is a vector, if any exception were to take place
	// after it is fully constructed, we can be sure that its 
	// memory will be de-allocated by its destructor.
	///////////////////////////////////////////////////////////////////////////
	vector<int> tempdiff(xic.size(),0); // initialize a temporary buffer to hold working results
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Evaluate the left difference
	//
	///////////////////////////////////////////////////////////////////////////
	if ( direction == -1 ) {
		for ( vector<int>::size_type i = 0; i < xic.size(); i++ ) {
			if ( i == 0) {
				tempdiff.at(i) = xic.at(i);
			}
			else {
				tempdiff.at(i) = xic.at(i) - xic.at(i-1);
			}
		}
	}
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Evaluate the right difference
	//
	///////////////////////////////////////////////////////////////////////////
	else if ( direction == 1) {
		for ( vector<int>::size_type i = 0; i < xic.size(); i++ ) {
			if ( i == xic.size()-1) {
				tempdiff.at(i) = xic.at(i);
			}
			else {
				tempdiff.at(i) = xic.at(i+1) - xic.at(i);
			}
		}
	}
	else {
		return -2;
	}

	///////////////////////////////////////////////////////////////////////////
	//
	//  Create a deep copy of the tempdiff object
	//  This will copy the contents of the tempdiff vector into the memory of
	//  diffxic vector.  Then, as soon as this function goes out of scope,
	//  the destructor of the tempdiff object will get called and all
	//  the data that was allocated for the tempdiff vector will be freed
	//  from the heap.  However, since we copied its contents into another
	//  vector that lives outside of this function we made sure to keep the
	//  results.
	//  The copy constructor of diffxic will get called, allocating all 
	//  required memory for performing the data copy automatically.
	//
	///////////////////////////////////////////////////////////////////////////
	diffxic = tempdiff;   // copy contents of tempdiff into diffxic.
	
	return 0; // return successfully
	
}  // end getdiff method













///////////////////////////////////////////////////////////////////////////////
//
// getR2 correlation statistic algorithm:
// Original algorithm developed in Matlab by: Jianqiu Zhang
//function rcorr=getR2Statistics(profile1,profile2,plength);
//if plength==0
//    rcorr=0;
//    return;
//    % add checks for sum(profile1,and profile2) > 0, and Stotal > 0
//else
//    Serr=0; Stotal=0;
//    normalizedprofile1=profile1/sum(profile1);
//    normalizedprofile2=profile2/sum(profile2);
//    for index=1:plength
//         Serr=Serr+(normalizedprofile1(index)-normalizedprofile2(index))^2; % Error squared between signals
//         Stotal=Stotal+(normalizedprofile1(index)+normalizedprofile2(index))^2/4; % average of 2 heights squared
//    end
//    rcorr=1-(Serr/Stotal);  % correlation( r^2 statistic )
//end    
//
// Parameters:
// Input: datavector1: a reference to a vector of doubles 
//        datavector2: a reference to a vector of doubles
// Output: 
//        rcorr:  a referece to a double, ( data will only be written to the result variable if return code == 0 )
// Return code: 0 successful, -1 ( vectors are not same size ), -2( datavector1 sum == 0), -3( datavector2 sum == 0 ), -4( Stotal == 0 )
//
// Refer to: /msda/matlab/interfacing/NelsonCode012412/src/getR2Statistics.m
//
// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 01/27/12
//
//
///////////////////////////////////////////////////////////////////////////////
int Utilities::getR2Statistic(vector<double> & datavector1, vector<double> & datavector2, double & rcorr) {
	

	///////////////////////////////////////////////////////////////////////////
	// 
	// Calculate the sum of each input vector
	// 
	///////////////////////////////////////////////////////////////////////////
	vector<double>::iterator itr;
	vector<double>::iterator itr1;
	vector<double>::iterator itr2;
	double datavector1Sum = 0; // sum of the contents of datavector1 ( always initialize stack variables )
	double datavector2Sum = 0; // sum of the contents of datavector2 ( always initialize stack variables )
	double Serr = 0;
	double Stotal = 0;
	
	///////////////////////////////////////////////////////////////////////////
	// 
	// 
	// The input data vectors must be of the same size.
	// 
	//  
	///////////////////////////////////////////////////////////////////////////
	if ( datavector1.size() != datavector2.size() ) {
		return -1;
	}
	
	//////////////////////////////////////////////////////////////////////////
	// Handle zero length vectors
	// The correlation of 2 zero length vectors will be returned as zero.
	//////////////////////////////////////////////////////////////////////////
	if( datavector1.size() == 0 ) {
		return -1;
	}
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Calculate the sum of the contents of data vector 1
	//
	///////////////////////////////////////////////////////////////////////////
	for ( itr = datavector1.begin(); itr < datavector1.end(); itr++ ) {
		datavector1Sum = *itr + datavector1Sum;  // keep running total
	}
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Calculate the sum of the contents of data vector 2
	//
	///////////////////////////////////////////////////////////////////////////
	for ( itr = datavector2.begin(); itr < datavector2.end(); itr++ ) {
			datavector2Sum = *itr + datavector2Sum;
	}
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Make sure to check if the sum of the input vectors is zero, since
	// we need to divide by them, this is critically important, since
	// division by zero is cause for an exception, and an unhandled excetion
	// is the cause of a crash.
	//
	///////////////////////////////////////////////////////////////////////////
	if ( datavector1Sum == 0 ) {
		return -2;	// return error code condition ( first input vector's contents add up to zero )
	}
	
	if ( datavector2Sum == 0 ) {
		return -3; // return error code condition ( second input vector's contents add up to zero )
	}
	
	///////////////////////////////////////////////////////////////////////////
	//  
	//  Data normalization:  Divide each vector data point by the 
	//  sum of the entire vector.
	//
	//  Perform an element by element division of each vector by the 
	//  sum of said vector.
	//  
	///////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////
	// Normalize input 1 
	///////////////////////////////////////////////////////////////////////////
	for ( itr = datavector1.begin(); itr < datavector1.end(); itr++ ) {
		*itr = *itr / datavector1Sum;  // normalize the first input vector by its sum
	}

	///////////////////////////////////////////////////////////////////////////
	// Normalize input 2
	///////////////////////////////////////////////////////////////////////////
	for ( itr = datavector2.begin(); itr < datavector2.end(); itr++ ) {
		*itr = *itr / datavector2Sum;  // normalize the second input vector by its sum
	}

	///////////////////////////////////////////////////////////////////////////
	// Evaluate error metric and correlation totals
	// At this point we know that both input vectors are the same size,
	// so we can use the size of either to loop over the data in order to 
	// evaluate 
	///////////////////////////////////////////////////////////////////////////
	for ( (itr1 = datavector1.begin(), itr2 = datavector2.begin()) ; itr1 < datavector1.end(); (itr1++,itr2++) ) {
		Serr = Serr + ((*itr1) - (*itr2)) * ((*itr1) - (*itr2));  // difference squared running sum
		Stotal = Stotal + (((*itr1) + (*itr2)) * ((*itr1) + (*itr2)))/4;  // average of 2 heights squared running sum
	}
	
	///////////////////////////////////////////////////////////////////////////
	//
	//  Important, we need to make sure that we do not divide by zero.
	//
	///////////////////////////////////////////////////////////////////////////
	if( Stotal == 0 ) {
		return -4;
	}
	
	
	///////////////////////////////////////////////////////////////////////////
	// Calculate the final correlation metric
	// The metric can be negative, 0, or positive double data item
	//
	//
	//  The result variable is ONLY modified on successful completion
	///////////////////////////////////////////////////////////////////////////
	rcorr = 1- (Serr/Stotal);  // copy result to output data

	return 0; // return successfully
	
}  // end getdiff method


///////////////////////////////////////////////////////////////////////////////
//
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
// Input: datavector1, a reference to a vector of doubles, containing smoothed signal data
//        datavector2: a reference to a vector of doubles, containing original signal data
//		  noiseThresholdLevel:  threshold level ( pass by value ), number of standard deviations( >= 0 ) 
//
// Output: 
//         noiseThreshold: the resulting noise threshold value
//
//
// Return code: 0 successful, -1 ( vectors are not same size ), -2( zero length input data vectors)
//				-3( Invalid noise threshold level )
//
// Refer to: /msda/matlab/interfacing/NelsonCode012412/src/getNoiseThreshold.m
//
// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 01/30/12
//
///////////////////////////////////////////////////////////////////////////////
int Utilities::getNoiseThreshold(vector<double> & datavector1, vector<double> & datavector2, double noiseThresholdLevel, double & noiseThreshold) {
	
	///////////////////////////////////////////////////////////////////////////
	// 
	// Calculate the sum of each input vector
	// 
	///////////////////////////////////////////////////////////////////////////
	vector<double>::iterator itr;
	vector<double>::iterator itr1;
	vector<double>::iterator itr2;
	vector<double>::iterator itr3;
	///////////////////////////////////////////////////////////////////////////
	// Instantiate temporary vector to hold the difference between
	// the two input vectors.
	///////////////////////////////////////////////////////////////////////////
    vector<double> noiseVector;

	double dataStd=0;
	double dataMean=0;
	///////////////////////////////////////////////////////////////////////////
	// 
	// 
	// The input data vectors must be of the same size.
	// 
	//  
	///////////////////////////////////////////////////////////////////////////
	if ( datavector1.size() != datavector2.size() ) {
		return -1;
	}
	
	///////////////////////////////////////////////////////////////////////////
	// Handle zero length vectors
	// The correlation of 2 zero length vectors will be returned as zero.
	///////////////////////////////////////////////////////////////////////////
	if( datavector1.size() == 0 ) {
		return -2;
	}
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Check the noiseThresholdLevel validity
	//
	///////////////////////////////////////////////////////////////////////////
	if ( noiseThresholdLevel < 0 ) {
		return -3; 
	}
	
	///////////////////////////////////////////////////////////////////////////
	// Calculate the difference between the two vectors, placing the 
	// results in a temporary vector.
	// Loop through the vector and evaluate the running difference.
	///////////////////////////////////////////////////////////////////////////
    for ( vector<double>::size_type i = 0; i < datavector1.size(); i++ ) {
        double temp  = datavector1.at(i) - datavector2.at(i);
        if ( temp > 0 ) {
            noiseVector.push_back( temp) ;
        }
    }
	


	
	///////////////////////////////////////////////////////////////////////////
	//  
	//  Calculate the mean , variance, and standard deviation
	//
	///////////////////////////////////////////////////////////////////////////
	cbi::MathUtils muObject;  // math utilities object
	
	int rc = 0;
	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	cout<<"noiseVector.size()="<<noiseVector.size()<<"noiseVector.at(0)"<<noiseVector.at(0)<<endl;
#endif
#endif
	
	//rc = muObject.stdAndMean(noiseVector,0,dataStd,dataMean);
	//rc = muObject.mean(noiseVector,dataMean);
	
	///////////////////////////////////////////////////////////////////////////
	// variance test
	// 8.341666666666667e+04
	///////////////////////////////////////////////////////////////////////////
	//rc = muObject.var(noiseVector,dataVar);
	//cout<<"Within Utilities::getNoiseThreshold"<<",rc="<<rc<<","<<dataVar<<endl;
	
	///////////////////////////////////////////////////////////////////////////
	// conditional variance test  
	// Only those items greater than 0
	//  temp = data1-data2;
	//  >> r = var(temp(temp>0))
	// 8.308350000000000e+04
	///////////////////////////////////////////////////////////////////////////
	//rc = muObject.var(noiseVector,0,dataVar);
	//cout<<"Within Utilities::getNoiseThreshold"<<",rc="<<rc<<","<<dataVar<<endl;
		
	///////////////////////////////////////////////////////////////////////////
	// conditional variance test  
	// Only those items greater than 0
	//  temp = data1-data2;
	//  >> r = var(temp(temp>0))
	// 8.308350000000000e+04
	///////////////////////////////////////////////////////////////////////////
	//rc = muObject.var(noiseVector,0,dataVar);
	//cout<<"Within Utilities::getNoiseThreshold"<<",rc="<<rc<<","<<dataVar<<endl;
		
	///////////////////////////////////////////////////////////////////////////
	// conditional variance test  
	// Only those items greater than 0
	//  temp = data1-data2;
	//  >> r = var(temp(temp>0))
	// 8.308350000000000e+04
	///////////////////////////////////////////////////////////////////////////
	rc = muObject.stdAndMean(noiseVector,0,dataStd,dataMean);

	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	cout<<"Within Utilities::getNoiseThreshold,stdAndMean"<<",rc="<<rc<<","<<dataStd<<","<<dataMean<<endl;
#endif
#endif
	
	///////////////////////////////////////////////////////////////////////////
	// Return the noise threshold
	///////////////////////////////////////////////////////////////////////////
	noiseThreshold = dataMean+dataStd*noiseThresholdLevel;
	
	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	cout<<"noiseThreshold="<<noiseThreshold<<endl;
#endif
#endif
	
	
	
	return 0; // return successfully
	
}  // end getNoiseThreshold method


///////////////////////////////////////////////////////////////////////////////
//
//
// movingAvgSmooth algorithm
//
//Original algorithm developed in Matlab by: Jianqiu Zhang
//function smoothXIC=movingAvgSmooth(XIC,numScans,numSmoothPoint)
//smoothXIC=zeros(1,numScans);
//for scanid=1:numScans
//    if scanid<=numSmoothPoint;
//       smoothBegin=1;
//    else
//       smoothBegin=scanid-numSmoothPoint;
//    end
//
//    if scanid>numScans-numSmoothPoint;
//       smoothEnd=numScans;
//    else
//       smoothEnd=scanid+numSmoothPoint;
//    end
//    tempSum=0;
//    for tempid=smoothBegin:smoothEnd
//       tempSum=XIC(tempid)+tempSum;
//    end    
//    smoothXIC(scanid)=tempSum/(smoothEnd-smoothBegin+1);
//
//end
//
//
// Parameters:
// Input: inputvector, a reference to a vector of doubles, original signal to be smoothed
//		  numSmoothPoints, the width of the smoothing kernel
//
// Output: 
//         outputvector, the smoothed result
//
//
// Return code: 0 successful, -1 ( vectors are not same size ), -2( zero length input data vectors),
//             -3 division by 0.
//
// Refer to: /msda/matlab/interfacing/NelsonCode012412/src/movingAvgSmooth.m
//
// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 02/07/12
//
///////////////////////////////////////////////////////////////////////////////
int Utilities::movingAvgSmooth(vector<double> & inputvector, int numSmoothPoints, vector<double> & outputvector) {
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Create a temporary vector to keep intermediate results
	//
	///////////////////////////////////////////////////////////////////////////
	int smoothBegin = 0;
	int smoothEnd = 0;
	double tempSum= 0;
	double denominator = 1;
	vector<double> tempResult(inputvector.size(),0);
	///////////////////////////////////////////////////////////////////////////
	// 
	// 
	// The input data vectors must be of the same size.
	// 
	//  
	///////////////////////////////////////////////////////////////////////////
	if ( inputvector.size() != outputvector.size() ) {
		return -1;
	}
	
	///////////////////////////////////////////////////////////////////////////
	// Handle zero length vectors
	// The correlation of 2 zero length vectors will be returned as zero.
	///////////////////////////////////////////////////////////////////////////
	if( inputvector.size() == 0 ) {
		return -2;
	}
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Perform the smoothing operation
	//
	///////////////////////////////////////////////////////////////////////////
	for ( vector<double>::size_type i = 0; i < inputvector.size(); i++ ) {
		///////////////////////////////////////////////////////////////////////
		//
		//  Check the start of the vector.  Since Matlab vectors start at
		//  index 1 and C++ vectors start at index 0, we need to compensate
		//  for the difference by subtracting 1 from numSmoothPoints.
		//  
		///////////////////////////////////////////////////////////////////////
		if(i <= (unsigned int)(numSmoothPoints-1) ){
			smoothBegin = 0;  // C++ vectors start at 0
		}
		else {
			smoothBegin= i-numSmoothPoints;  // Offset by numSmoothPoints
		}

		///////////////////////////////////////////////////////////////////////
		//
		// Check the end of the vector.  We need to compensate for difference
		// with Matlab's indices starting at 1 and C++ starting at 0. 
		//
		///////////////////////////////////////////////////////////////////////
		if( i >= (unsigned int)(inputvector.size()-numSmoothPoints)) {
			smoothEnd = inputvector.size()-1;  // set to max index
		}
		else {
			smoothEnd = i+numSmoothPoints;  // set right boundary
		}

		///////////////////////////////////////////////////////////////////////
		//
		// Initially the variable that will store the temporary
		// kernel sum.
		//
		///////////////////////////////////////////////////////////////////////
		tempSum = 0;

		///////////////////////////////////////////////////////////////////////
		//
		// Go through the data calculating an average kernel.
		// [smoothBegin smoothEnd] inclusive
		//
		///////////////////////////////////////////////////////////////////////
		for ( int j = smoothBegin; j <= smoothEnd; j++ ) {
			tempSum = inputvector.at(j) + tempSum;
		}
		
		///////////////////////////////////////////////////////////////////////
		// The denominator contains the normalization factor
		///////////////////////////////////////////////////////////////////////
		denominator = smoothEnd-smoothBegin+1;
		
		///////////////////////////////////////////////////////////////////////
		//
		//  Evaluate the new value to place at location i.
		//  Note:  We must make sure that we check any divisions for 
		//  a divide-by-zero exception.
		//
		///////////////////////////////////////////////////////////////////////
		if (denominator > 0) {
			tempResult.at(i) = tempSum/denominator;
		}
		else {
			// divide by zero case
			return -3;
		}

	}  // end for loop

	///////////////////////////////////////////////////////////////////////////
	// Return result 
	//  We only want to change the 
	///////////////////////////////////////////////////////////////////////////
	outputvector = tempResult;
	return 0; // return successfully
}
	
	


///////////////////////////////////////////////////////////////////////////////
//
// This function adds a new interval.
// Port from Matlab code:
// msda/matlab/interfacing/NelsonCode02072012/src/addNewIntervals.m
//
// Inputs:
// 1) intervalList:   1x2 double ( one row per interval )
//                 333 347
//                 next interval start and end value
//  - Note:  Implemented as intervalListStart and intervalListEnd
// 	  An intervalList is a 2 column table, with a row for each interval
//    i.e.  row1 484, 508
//          row2 432, 445
//          row3 446, 461
//    intervalListStart contains the start value of an interval
//    intervalListEnd contains the end value of an interval
//    i.e. intervalListStart.at(0) --> should yield 484 
//    i.e. intervalListEnd.at(0)   --> should yield 500
//    i.e. intervalListStart.at(1) --> should yield 432
//    i.e. intervalListEnd.at(1)   --> should yield 445
//  
//    intervalListStart.at(#), intervalListEnd.at(#) are used to 
//    index the XIC vector.  The intervals specify the sections of the 
//    XIC vector
//
//
//
//  intervalListStart    intervalListEnd		XIC vector
//  1						5			--> XIC.at(1),..XIC.at(5) 
//  20						30			--> XIC.at(20),..XIC.at(30)
//	50						60			--> XIC.at(50),..XIC.at(60)
//	75						100  		--> XIC.at(75),..XIC.at(100)
//
// 2) XIC:  1xN vector of doubles.
//    This is the current XIC signal that is being analyzed to see 
//    how each section of the signal as specified by the intervalList
//    data structure maps to a set of "boxes" identified by groupID.
//
//    The XIC signal has already been cut into intervals, where each
//    interval ideally has only one peak within it.  This generated
//    by the interval detection algorithm, performed prior to 
//    calling this method.
//    
// 3) windowMZ:  This is the center MZ value of the XIC signal.
//    Since an XIC signal is a retention time axis view of the mz versus
//    intensity data, an XIC signal is the sum of the intensities within 
//    an mz window.
//
//
//
//
//
// 
//
//
//
//
// - Note: intervalCount:  Not needed, implemented as intervalListStart.size()
// - Note: For each invocation, there will be a new interval list.
// - Note: By using STL vectors we can leave memeory cleanup to the vector 
//         object.
//
//
// XIC:  1x541 double
//
//
//
//
// windowMZ:  1x1 double ( i.e. 750.2703 )
// preLC: 
// 
// 
// 
///////////////////////////////////////////////////////////////////////////////
int Utilities::addNewIntervals(		vector<int> & intervalListStart, 
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
									vector<int> & outputGroupMap ) {
	///////////////////////////////////////////////////////////////////////////
	//% intervalList --> newly cut regions from current XIC signal to be appended
	//% to preLC structure.
	//    Each XIC signal maps to the data from a single color in the 
	//    MSDA Level 1 viewer.  
	//     interval # --> start, end
	//     interval # --> start, end
	//     interval # --> start, end
	//     interval # --> start, end
	//  --> I need to more closely diagram the meaning of the intervalList
	//
	//
	//% intervalCount --> # of new intervals
	//  --> How are the number of newintervals determined
	//% XIC --> signal
	//
	//
	//% windowMZ --> center MZ value of the XIC
	//
	//
	//% preLC --> data structure ( contains the list of elements and their box
	//
	//
	//% number(aka group ID )
	//
	//
	//% massResolution --> instrument specific determines peak width
	//
	//
	//% preLCLength--> total # of elements already contained in the preLC structure
	//
	//
	//% GroupMap--> Box ids
	//
	//
	//% R2Threshold --> correlation threshold
	//
	//
	///////////////////////////////////////////////////////////////////////////
	
	
	//iteration =0 --> length(groupmap)=1,preLCLength=1,intervalcount=1  --> newGroupMap length = preLCLength+intervalcount
	//iteration =1 --> length(groupmap)=2,preLCLength=2,intervalcount=1
	//iteration =2 --> length(groupmap)=3,preLCLength=3,intervalcount=2
	//iteration =3 --> length(groupmap)=5,preLCLength=5,intervalcount=1
	//iteration =4 --> length(groupmap)=6,preLCLength=6,intervalcount=2

	
	
	
	
	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	// We need to be able to see the input parameter properties for 
	// debug purposes.
	cout<<"Utilities::addNewIntervals"<<endl;
	cout<<"intervalListStart.size()="<<intervalListStart.size()<<endl;
	cout<<"intervalListEnd.size()="<<intervalListEnd.size()<<endl;
	cout<<"XIC.size()="<<XIC.size()<<endl;
	cout<<"plcGroupID.size()="<<plcGroupID.size()<<endl;
	cout<<"plcMz.size()="<<plcMz.size()<<endl;
	cout<<"plcElutionStart.size()="<<plcElutionStart.size()<<endl;
	cout<<"plcElutionEnd.size()="<<plcElutionEnd.size()<<endl;
	cout<<"plcElutionPeak.size()="<<plcElutionPeak.size()<<endl;
	cout<<"massResolution="<<massResolution<<endl;
	cout<<"overlapPercentage="<<overlapPercentage<<endl;
	cout<<"groupMap.size()="<<groupMap.size()<<endl;
	cout<<"R2Threshold="<<R2Threshold<<endl<<endl;
	cout<<"outputGroupMap.size()"<<outputGroupMap.size()<<endl;
#endif
#endif
	
	
	///////////////////////////////////////////////////////////////////////////
	//  Input Validity tests
	///////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////
	//  Case:1  
	//  intervalListStart and intervalListEnd must map into valid
	//  ranges of the XIC input vector.
	//
	//  intervalListStart    intervalListEnd
	//  1						5			--> XIC.at(1),..XIC.at(5) 
	//  20						30			--> XIC.at(20),..XIC.at(30)
	//	50						60			--> XIC.at(50),..XIC.at(60)
	//	75						100  		--> XIC.at(75),..XIC.at(100)
	//
	//
	///////////////////////////////////////////////////////////////////////////
	
	assert(intervalListStart.size() == intervalListEnd.size() );
	
	if( intervalListStart.size() != intervalListEnd.size() ) {
		return -1;  // interval length mismatch
	}
	
	int maxIntervalListStart = 0;
	int maxIntervalListEnd = 0;
	cbi::SignalProcessor sgProcessor;
	sgProcessor.findMaxValue(intervalListStart,maxIntervalListStart );  // find max value of intervalListStart
	sgProcessor.findMaxValue(intervalListEnd, maxIntervalListEnd); // find max value of intervalListEnd
	
	if (maxIntervalListStart < 0 ){
		return -2; // invalid maxIntervalListStart
	}
	if(maxIntervalListEnd < 0) {
		return -3; // invalid maxIntervalListEnd
	}
	
	

	
	
	
	
	
	///////////////////////////////////////////////////////////////////////////
	// START Local variables
	///////////////////////////////////////////////////////////////////////////
	double mzPeakWidth = 0;
	
	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	cout<<"Utilities::addNewIntervals, plcGroupID.size()="<<plcGroupID.size()<<endl;
	cout<<"Utilities::addNewIntervals, intervalListStart.size()="<<intervalListStart.size()<<endl;
#endif
#endif
	
	///////////////////////////////////////////////////////////////////////////
	// Pre-Allocate memory for the output group map
	// We know exactly how many data elements we will need in the
	// output group map.
	///////////////////////////////////////////////////////////////////////////
	//% First decide preLC that has matching MZ range with windowMZ. If the
	//% difference is within mzPeakWidth, then return the group ID. 
	
	vector<int> tempOutputGroupMap;
	vector<int> matchingPreLCIDs;
	
	
	int preLCIndex = 0;
	int preLCLength = 0;
	int maxGroupID = 0;  // Will contain the maximum value in the GroupMap vector
	int rc = 0; // return code holding variable
	
	///////////////////////////////////////////////////////////////////////////
	// Decision data 
	//  Each of these boolean variables is used to keep the result
	//  of one of the decisions within the logic of this method.
	///////////////////////////////////////////////////////////////////////////
	bool check1 = false;
	bool check2 = false;
	bool check3 = false;
	bool check4 = false;
	bool check5 = false;
	bool check6 = false;
	bool check7 = false;
	bool check8 = false;
	
	///////////////////////////////////////////////////////////////////////////
	//  We need to store the original size of the plc data structures so that 
	//  we can keep a record of the original amount, since we add new data
	//  items into the plc data structures.
	///////////////////////////////////////////////////////////////////////////
	preLCLength= plcGroupID.size();  // save the original size of the plc data structure
	
	///////////////////////////////////////////////////////////////////////////
	// END Local variables
	///////////////////////////////////////////////////////////////////////////
	if (massResolution == 0 ){
		
		
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
		cout<<"Utilities::addNewIntervals, massResolution  == 0"<<endl;
#endif
#endif
		return -2;
	}
	
	///////////////////////////////////////////////////////////////////////////
	// Calculate the mzPeakWidth 
	///////////////////////////////////////////////////////////////////////////
	mzPeakWidth = windowMZ/massResolution/(2*sqrt(log(4)));  // checked massResoltion divide by zero case

	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
		// Matlab = 0.007965261443519
		// C++ = 0.00796526
		// Basic test passes.
		cout<<"Utilities::addNewIntervals,mzPeakWidth="<<mzPeakWidth<<endl;
		
		cout<<"Utilities::tempOutputGroupMap.size()="<<tempOutputGroupMap.size()<<endl;
#endif
#endif	
	
	///////////////////////////////////////////////////////////////////////////
	// 
	// Create a copy of the first plcGroupID.size() elements in groupMap
	// and place them into tempOutputGroupMap.
	// 
	// 
	// 
	// K>> intervalList
	//
	//	intervalList =
	//
	//	   240   278
	//	   362   394
	//	    69   107
	//	   496   512
	//	   513   526
	//	   479   485
	//  
	// 
	///////////////////////////////////////////////////////////////////////////
	
	// The vectors containing the start and end interval regions must have the same
	// size.

	assert( intervalListStart.size() == intervalListEnd.size() );	
	if( intervalListStart.size() != intervalListEnd.size() ) {
		return -2;  // interval list start and end vectors are not the same size
	}
	
	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	cout<<"Utilities::intervalListStart.size()"<<intervalListStart.size()<<endl;
#endif
#endif

	
	///////////////////////////////////////////////////////////////////////////
	// Copy the old group map into a new larger group map
	// This loop is suspect....  Is the size of tempOutputGroupMap always
	// greater than or equal to the size of groupMap
	// newGroupMap=zeros(preLCLength+intervalCount,1);
	// newGroupMap(1:preLCLength)=GroupMap;
	//
	// That is, is preLCLength+intervalCount >= length(GroupMap)
	///////////////////////////////////////////////////////////////////////////
	for( vector<int>::size_type i = 0; i < groupMap.size(); i++ ) {
		tempOutputGroupMap.push_back( groupMap.at(i) ); 
	}
	
	///////////////////////////////////////////////////////////////////////////
	// M.Zhang Matlab comment:"
	//	% First decide preLC that has matching MZ range with windowMZ. If the
	//  % difference is within mzPeakWidth, then return the group ID. "
	// Note:
	// vector<int> matchingPreLCIDs(plcGroupID.size(),0);
	//
	// Maps to:
	// matchingPreLCIDs=zeros(preLCLength,1);
	//
	//
	//
	///////////////////////////////////////////////////////////////////////////
	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	cout<<"preLCIndex="<<preLCIndex<<endl;
	cout<<"preLCLength="<<preLCLength<<endl;
#endif
#endif
	
	///////////////////////////////////////////////////////////////////////////
	//% Go through all the preLCs, 
	//% We check the membership in each box...  it has to be matched
	//% in mz first( pre-filter ).  Since the whole xic only has a single mz
	//% value, we only need to go through this ones.
	//% Each new xic, needs to be checked with the list of all previously
	//% determined boxes, to see which boxes it falls on in mz.
	//% Idea: Keep a stack of the latest boxes,  the results are always in the list
	//% of latest boxes. ( this will minimize the number of total comparisons).
	///////////////////////////////////////////////////////////////////////////
	for( preLCIndex = 0; preLCIndex < preLCLength; preLCIndex++ ) {
		///////////////////////////////////////////////////////////////////////
		//
		// Go through each pre-LC candidate, and check whether it 
		// the mz value of the plc candidate falls within the current 
		// mz window.  
		// 
		///////////////////////////////////////////////////////////////////////
		// we need to check the membership in each box
		if (abs(plcMz.at(preLCIndex)-windowMZ ) < mzPeakWidth) {   // find the distance in mz
			
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
			cout<<"MatchingPreLCIDs.size()"<<matchingPreLCIDs.size()<<endl;
			cout<<"abs(plcMz.at(preLCIndex)-windowMZ ) < mzPeakWidth"<<endl;
#endif
#endif
			///////////////////////////////////////////////////////////////////
			//  We do not need to keep an counts variable
			//  since we can 
			///////////////////////////////////////////////////////////////////
			
			// save the index of the match
			matchingPreLCIDs.push_back(preLCIndex);
		} // end finding the distance in mz
	}  // end looping through all the preLCs, finding the membership in each box
	// The whole xic only has a single mz value assocated with it at this point
	
	///////////////////////////////////////////////////////////////////////////
	// Find the maximum value contained within the GroupMap
	//
	//K>> maxGroupID
	//maxGroupID =
	//     2
	//K>> GroupMap
	//GroupMap =
	//     0
	//     1
	//     1
	//     2
	//     1
	//     2
	//
	//
	///////////////////////////////////////////////////////////////////////////
	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	vector<int>::iterator itr;
	cout<<"groupMap="<<endl;
	for ( itr = groupMap.begin(); itr < groupMap.end(); itr++ ) {
		
		cout<<*itr<<endl;
	}
#endif
#endif
	
	cbi::SignalProcessor spObject;  // Signal Processor Object
	rc = spObject.findMaxValue(groupMap, maxGroupID );  // Find the maximum box id #

	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	cout<<"rc= "<<maxGroupID<<endl;
#endif
#endif
	
	
	///////////////////////////////////////////////////////////////////////////
	//
	// The Matlab code pre-allocates the matchingPreLCIDs data structure to 
	// be the same length as the # of preLC items.
	// Since not all of the preLC items will match the current box, the 
	// Matlab code resize the matchingPreLCIDs data structure to be
	// the same size as the value in matching PreLCIDsCounts.
	//
	//
	///////////////////////////////////////////////////////////////////////////
	
		
	///////////////////////////////////////////////////////////////////////////
	//
	// %Now for each interval, search if it belongs to one group, one interval
	// %can belong to multiple groups by checking overlapping retention time
	// %range.
	//
	// % Check for scan range... the be considered together, they must also
	// % overlap in scan#( rt dimension)
	//
	///////////////////////////////////////////////////////////////////////////
	
	// Loop through all intervals
	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	cout<<"intervalListStart.size()"<<intervalListStart.size()<<endl;
	cout<<"matchingPreLCIDs.size()"<<matchingPreLCIDs.size()<<endl;
#endif
#endif
	
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	//  End loop through all intervals
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	// i maps to intervalIndex
	for( vector<int>::size_type i = 0; i < intervalListStart.size(); i++ ) {
		///////////////////////////////////////////////////////////////////////
		// Local variable that will hold the list of of intervals
		// that belong to multiple groups.
		///////////////////////////////////////////////////////////////////////
		vector<int> doubleMatchedIDs; // this is a loop local variable
		// This variable will get deleted and recreated each time 
		// through this loop.
		///////////////////////////////////////////////////////////////////////
		//
		// Loop through each matchingPreLCID
		//
		//  % finding the double matches, gets us to the final box...
        //  % checks rt dimension overlap for those who passed the mz dimension
        //  % filter ( matchingPreLCIDs contain the list of intervals that
        //  % match in mz values )... Go through those preLC elements......
        //  % If the end is less than the start of new interval( case 1) -->
        //  % no overlap
        //  % If the start is greater than the end of the new interval( case 2) --> no overlap
        //  % If these 2 conditions are not met, then we have an overlap
		//
		//
		//
		///////////////////////////////////////////////////////////////////////
		// j maps to index
		for( vector<int>::size_type j =0; j < matchingPreLCIDs.size(); j++ ) {
			
			
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
			cout<<"i="<<i<<"j="<<j<<endl;
#endif
#endif
			
			///////////////////////////////////////////////////////////////////
			//
			// Calculate overlap and make decision based on the overlap
			// For speed improvement, we can have 2 versions, one for
			// debugging and one for deployed version.
			// debugg version can use the checked at() method, and the 
			// deployed version can use the unchecked [] operator. 
			//
			//  % finding the double matches, gets us to the final box...
			//  % checks rt dimension overlap for those who passed the mz dimension
			//  % filter ( matchingPreLCIDs contain the list of intervals that
			//  % match in mz values )... Go through those preLC elements......
			//  % If the end is less than the start of new interval( case 1) -->
			//  % no overlap
			//  % If the start is greater than the end of the new interval( case 2) --> no overlap
			//  % If these 2 conditions are not met, then we have an overlap
			//
			///////////////////////////////////////////////////////////////////

			check1 = plcElutionEnd.at(matchingPreLCIDs.at(j)) < intervalListStart.at(i);
			check2 = plcElutionStart.at(matchingPreLCIDs.at(j)) > intervalListEnd.at(i);
			check3 = ~(check1 || check2);
			if ( check3 ) {
				
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
				cout<<"within checkFinal if block"<<endl;
#endif
#endif
				///////////////////////////////////////////////////////////////
				// Calculate the length of the interval:
				// intervalListEnd.at(i) - intervalListStart.at(i) + 1
				// 
				///////////////////////////////////////////////////////////////
				int intervall = 0;
				intervall =  intervalListEnd.at(i) - intervalListStart.at(i) + 1;
				doubleMatchedIDs.push_back(matchingPreLCIDs.at(j));


				///////////////////////////////////////////////////////////////
				//  Check overlap in the retention time dimension
				///////////////////////////////////////////////////////////////
				//   % if there is overlap, it should be significant.......
				//           %  Make the 0.3 a parameter.... [overlappercentage].
				//            %   0.3 --> 30% overlap in boxes
				//            %
				//            %  When there is a significant overlap,
				///////////////////////////////////////////////////////////////
				check4 = abs(plcElutionEnd.at(matchingPreLCIDs.at(j)) - intervalListEnd.at(i) ) < (overlapPercentage*intervall);
				check5 = abs(plcElutionStart.at(matchingPreLCIDs.at(j)) - intervalListStart.at(i) ) < (overlapPercentage*intervall);
				check6 = check4 && check5;
				
				if( check6 ) {
					///////////////////////////////////////////////////////////
					//   % then we declare a double match
					///////////////////////////////////////////////////////////
					doubleMatchedIDs.push_back(matchingPreLCIDs.at(j));
				}  // end retention time dimension check

			} // End overlap test

		}  // End loop through all matchingPreLCIDs 


		///////////////////////////////////////////////////////////////////////
		//
		// Process double-matched regions
		// 
		//   % For those double matched prelc elements,  we check for the correlations
		//   % Those peaks that got into the same box.. we want to find their
		//   % correlation.
		//
		///////////////////////////////////////////////////////////////////////
		if ( doubleMatchedIDs.size() == 0 ) {
			///////////////////////////////////////////////////////////////////
			//
			//  preLC(preLCLength+intervalIndex).GroupID=maxGroupID+1;
	        //  maxGroupID=maxGroupID+1;
			//
			///////////////////////////////////////////////////////////////////
			// % For this prelc element, we didn't find any box to put it in,
			// % therefore, we need to make a new box....
			///////////////////////////////////////////////////////////////////
			
			///////////////////////////////////////////////////////////////////
			// Append an item to the preLC data structure
			// Since we didn't find a box to put this xic into,  we need to 
			// create a new box.
			// By pushing maxGroupID + 1 onto the plcGroupID data 
			// structure, we are creating a new box, ( a new Pre-LC
			// candidate.
			///////////////////////////////////////////////////////////////////
			
			plcGroupID.push_back( maxGroupID + 1);
			
			///////////////////////////////////////////////////////////////////
			// 
			//  We do not need to call the findMax function again,
			//  to keep track of the maximum value of the group identifier.
			//  All we have to do is to incremenet it by 1 each time
			//  we create a new box( a new plcGroupID ) element.
			// 
			// 
			///////////////////////////////////////////////////////////////////
			maxGroupID++; // increment by 1, a new box was created.
			
		}
		else {
			///////////////////////////////////////////////////////////////////
			//
			//% Now we find a group of matched pre-LC peaks. Calculate the peak shape
			//% correlation and decide if the pre-LC peak should be assigned to an
			//% existing group or should be assigned to a new group.
			//
			// The current XIC matches one or more boxes.  It's match was
			// made based on MZ and Retention time distance.  However,
			// this is not enough to add a new XIC to a box.  We must also
			// make sure that the peak shapes that are already within
			// the boxes we are considering are similar in shape to our current
			// peak.,  
			// For example, even if our current xic signal falls inside the 
			// mz/scan# range of one or more boxes, but if the shape is 
			// completely different from the peak shapes already belonging
			// to those boxes, we don't want to be adding our new xic signal
			// into these boxes.  It would be better to create a new box
			// to add this xic into.
			// In this way, we are in essense clustering by mz, then scan#(retention time),
			// and finally by peak shape.
			//  
			// This helps to reduce the amount of computation required for evaluating
			// the lc candidates.  There is no need in calculating a correlation
			// metric between 2 xic signals that are aren't close in mz, or retention time.
			//
			//
			//
			///////////////////////////////////////////////////////////////////
			
			///////////////////////////////////////////////////////////////////
			//
			//  % We found at least 1 box that we can put this preLC signal in.
			//  Before attaching this preLC signal to this box, we must 
			//  first make sure that the peak shape is similar to those
			//  preLC signals that are already inside of that box
			//
			//
			///////////////////////////////////////////////////////////////////
			
			vector<double> rcorr;  	// a holder for the correlation results,
									// one result per double match.
			vector<int> groupIDs;  	// keep track of the group number of the double matches
			
			///////////////////////////////////////////////////////////////////
			//
			// Variables to keep track of the overlap positions
			// 
			// These 4 variables are used to cut the source data into
			// profile1 and profile2 vectors.  The correlation is performed
			// with the profile1 and profile2 vectors.
			///////////////////////////////////////////////////////////////////
			int preLCProfileStart = 0;   // Used to index 
			int currentProfileStart = 0;
			int currentProfileEnd = 0;
			int preLCProfileEndDiff = 0;
			int pLength = 0; // Profile vector length
			
			///////////////////////////////////////////////////////////////////
			//
			// Loop through the double matches:
			//
			///////////////////////////////////////////////////////////////////
			for ( vector<int>::size_type k = 0; k< doubleMatchedIDs.size(); k++ ) {
					///////////////////////////////////////////////////////////
					//
					//
					// % We have 2 preLC xic's, but they have different start and end
					// % scans, so we need to take out the common part to be able to
					// % perform a valid comparison.  
					// % Find overlap...
					// % We have a preLC element, then we have a current profile...  the
					// % start of the preLC element is greater than the start of the
					// % new interval.  
					//
					//
					///////////////////////////////////////////////////////////
				
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
					cout<<"doubleMatchedIDs,k="<<k<<endl;
#endif
#endif
					///////////////////////////////////////////////////////////
					// % If the current element starts early, then use the starting
					// % point of the preLC element.
					///////////////////////////////////////////////////////////
					check7 = plcElutionStart.at(doubleMatchedIDs.at(k)) >= intervalListStart.at(i);
					
					if ( check7 ) {
						preLCProfileStart = 1;  // Need to check where this is used.
						currentProfileStart = plcElutionStart.at(doubleMatchedIDs.at(k));
					}   // end if block
					else {
						///////////////////////////////////////////////////////
						// If the current element starts later, then use the 
						// starting value of the current element.
						///////////////////////////////////////////////////////
						preLCProfileStart = intervalListStart.at(i) - plcElutionStart.at(doubleMatchedIDs.at(k)) +1 ;  // Need to check where this is used.
						currentProfileStart = intervalListStart.at(i);					
					}  // end else block
					
							
					///////////////////////////////////////////////////////////
					// 
					// 
					//   % Handle the end case
					//   % take care of the end.
			        //   % If the preLC element has an that is earlier than the current
			        //   % interval, then the end of the preLC element is used and we must
			        //   % cut the current element to match the preLC element's end.
					// 
					// 
					///////////////////////////////////////////////////////////
					check8 = plcElutionEnd.at(doubleMatchedIDs.at(k)) <= intervalListEnd.at(i);
					
					if( check8 ) {
						///////////////////////////////////////////////////////
						//
						// preLC element is earlier than the current interval,
						// cut the current element to match the preLC element's end
						//
						///////////////////////////////////////////////////////
						preLCProfileEndDiff = 0;
						currentProfileEnd = plcElutionEnd.at(doubleMatchedIDs.at(k));
					}
					else {
						///////////////////////////////////////////////////////
						// % If the preLC element has an end that is later than the
						// % current interval, the the end of the current element is used
						// % and we must do nothing to the current element.
						///////////////////////////////////////////////////////
						preLCProfileEndDiff = plcElutionEnd.at(doubleMatchedIDs.at(k)) - intervalListEnd.at(i);
						currentProfileEnd = intervalListEnd.at(i);
						
					}
					
					
					///////////////////////////////////////////////////////////
					//  
					//  Correlation calculation.  
					//  By this point the two signals have been trimmed so that
					//  they contain the same number of data elements.
					//  Now, we can call the getR2Statistics function to obtain
					// 
					// 
					///////////////////////////////////////////////////////////
					
					///////////////////////////////////////////////////////////
					//   % Profile 1 is from the previous LC.  Cut the preLC element, as
				    //   % profile 1
					///////////////////////////////////////////////////////////
					
					///////////////////////////////////////////////////////////
					// We need to make a temporary copy of the data in 
					// profile 1.
					///////////////////////////////////////////////////////////
					vector<double> profile1;  // Will be populate from the previous LC
					vector<double> profile2;  // Will be populated from the current LC
					
					///////////////////////////////////////////////////////////
					// 
					//
					// Copy the data from the elutionPeak data structure
					// of the preLC data structure according to the cut 
					// points determined in the previous overlap analysis stage.
					// 
					// 
					///////////////////////////////////////////////////////////
					for ( vector<int>::size_type n = preLCProfileStart; n <= plcElutionPeak.at(doubleMatchedIDs.at(k)).size() - preLCProfileEndDiff; n++ ) {
						profile1.push_back(plcElutionPeak.at(doubleMatchedIDs.at(k)).at(n));
					}
					
					
					///////////////////////////////////////////////////////////
					// 
					//
					// Copy the data from the XIC 
					//
					// 
					///////////////////////////////////////////////////////////
					for( vector<int>::size_type n = (unsigned)currentProfileStart; n<= (unsigned)currentProfileEnd; n++ ) {
						profile2.push_back(XIC.at(n));
					}
					
					///////////////////////////////////////////////////////////
					//
					// Calculate the profile length
					// 
					///////////////////////////////////////////////////////////
					pLength  = currentProfileEnd-currentProfileStart +1;
					
					///////////////////////////////////////////////////////////
					// The length of the two profiles should match
					///////////////////////////////////////////////////////////
					assert(profile1.size() == profile2.size());
					
					
					///////////////////////////////////////////////////////////
					//
					// I'm pretty sure the following is a bug:
					// index should be index2.....
				    // rcorr(index)=getR2Statistics(profile1,profile2,plength);
				    // groupIDs(index)=preLC(doubleMatchedIDs(index2)).GroupID;
					// rCorrValues(i) 
					// At this point the loop that uses index has already completed
					//
					///////////////////////////////////////////////////////////
					
					
					///////////////////////////////////////////////////////////
					// Evaluate the correlation between the 2 signals
					///////////////////////////////////////////////////////////
					double tempR2Corr = 0;
					rc = Utilities::getR2Statistic(profile1,profile2,tempR2Corr);
					rcorr.push_back(tempR2Corr);
					
					///////////////////////////////////////////////////////////
					// Record the group ids for later reference.  This
					// tells us which boxes this correlation belongs to.
					///////////////////////////////////////////////////////////
					groupIDs.push_back(plcGroupID.at(doubleMatchedIDs.at(k)));
					
					
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
					cout<<"rc="<<rc<<endl;
					
#endif
#endif
					
					///////////////////////////////////////////////////////////
					//  At this point both rcorr and groupIDs must have 
					//  the same length
					///////////////////////////////////////////////////////////
					assert(rcorr.size() == groupIDs.size());
					
					
			} // end loop through the double matches
			
			
			///////////////////////////////////////////////////////////////////
			//  Find the maximum correlation, we also need to know 
			//  the index of the maximum correlation.
			///////////////////////////////////////////////////////////////////
			cbi::SignalProcessor sgObject;  // signal processing utilities object
			vector<int>::size_type maxIndex = 0;
			double maxCorr = 0;
			
			///////////////////////////////////////////////////////////////////
			//  Find the maximum value 
			//  Note, this maps to 
			// [maxcorr,maxID]=max(rcorr);  % Find the maximum correlation
			///////////////////////////////////////////////////////////////////
			rc = sgObject.findMaxValue(rcorr,maxIndex,maxCorr);
			
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
			cout<<"rc="<<rc<<endl;
#endif
#endif
		
			///////////////////////////////////////////////////////////////////
			//
			// Check threshold criteria:
			// 
			// Design question:
			// The following check only handles the case with the highest
			// correlation,   what about cases where there are more 
			// than one correlation value that meets the threshold criteria?
			// 
			// 
			// 
			///////////////////////////////////////////////////////////////////
			if( maxCorr > R2Threshold ) {
				
				///////////////////////////////////////////////////////////////
				//
				// Update the preLC data structure
				// preLC(preLCLength+intervalIndex).GroupID=groupIDs(maxID);
				// For example,
				// groupIDs = 0 0 0 0 0 5  ( This was a bug in the
				//  first version of the Matlab code. ).  Instead of 
				//  using all the data items, only the last one was
				//  ever getting used, the other ones were
				//  being left at 0.
				// maxID = 6
				// groupIDs(maxID) = 5
				// preLCLength = 22 ( the original data structure length )
				// preLC(27).GroupID = groupIDs(maxID) = 5
				// 
				// If we pass the threshold test, then we need only update the 
				// groupID information in the plc data structure.
				// 
				// If we don't pass the threshold test, this means that
				// the membership criteria has not been met to include
				// this xic within any of the boxes, therefore:
				//
				// We need to make a new box:  This is handled during the
				// else clause.
				//
				// Save the box# of the box that our signal had the highest
				// correlation to.
				//  
				// An analogy, we're trying to find the "box" where our
				// signal should go.
				//
				//
				//
				//
				///////////////////////////////////////////////////////////////
				
				plcGroupID.at( preLCLength +  i)= groupIDs.at(maxIndex); // update group membership ( a.k.a change box # )
				
				
			}  // End processing 
			else {
				///////////////////////////////////////////////////////////////
				// If the highest "box" membership criteria has not been 
				// met by the highest correlation case, then neither will
				// the others, and so we need a new "box".  However,
				// this presents an interesting design point, on how to 
				// handle more than 1 box meeting the membership elegibility
				// criteria.
				//  ( 02/14/12 Design meeting with M.Zhang: 
				//  - Add additional criteria, such that multiple 
				//    boxes that have above threshold correlations are 
				//    then passed through an MZ distance filter.
				//    1) The one with the closest value in MZ is selected, if 
				//       there are multiple boxes with a similar above
				//       threshold correlation.
				///////////////////////////////////////////////////////////////
				
				///////////////////////////////////////////////////////////////
				//
				// Update the preLC data structure
				// increment the maxGroupID index
				//
				///////////////////////////////////////////////////////////////
				plcGroupID.at( preLCLength +  i)= maxGroupID+1;
				
				///////////////////////////////////////////////////////////////
				// 
				// Since we created a new box, we need to increment our
				// box count.
				// 
				///////////////////////////////////////////////////////////////
				maxGroupID++; // increment by 1, a new box was created.
				
			}  // end else block, make new "box" case.
			
			
			
		}  // end else block
		
		

		///////////////////////////////////////////////////////////////////////
		//
		//  Finish filling in preLC data structure
		//
		//   % We've determined which boxes our current interval belongs to, then,
		//   % populate the rest of the preLC data structure fields. 
		//   % Populating the other fields of the LC element.
		// 
		//
		//
		///////////////////////////////////////////////////////////////////////
		
		tempOutputGroupMap.push_back(plcGroupID.at(preLCLength + i));
		
		///////////////////////////////////////////////////////////////////////
		// Update the PLC data structure
		///////////////////////////////////////////////////////////////////////
		plcMz.push_back(windowMZ);
		plcElutionStart.push_back(intervalListStart.at(i));
		plcElutionEnd.push_back(intervalListEnd.at(i));
	
		///////////////////////////////////////////////////////////////////////
		// Update the PLC elution profile data
		///////////////////////////////////////////////////////////////////////
		vector<double> tempPlcElutionPeak;
		//for( vector<int>::size_type n = (unsigned)intervalListStart.at(i); n<= (unsigned)intervalListEnd.at(i); n++ ) {
		//	tempPlcElutionPeak.push_back(XIC.at(n));
		//}
		// Temporary,  we need to make sure that the intervalListStart and intervalListEnd are 
		// valid intervals, that can be used to index into the XIC vector.
		plcElutionPeak.push_back(tempPlcElutionPeak);
		
		///////////////////////////////////////////////////////////////////////
		// At this point the following needs to be valid:
		// All PLC data structure vectors need to be the same length.
		// If the PLC data structures are not the same length, then 
		// there is a problem and we need to exit with an error code.
		///////////////////////////////////////////////////////////////////////
		
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
		cout<<"plc sizes()"<<plcGroupID.size()<<","<<plcMz.size()<<","<<plcElutionStart.size()<<","<< plcElutionEnd.size()<<","<<plcElutionPeak.size()<<endl;
		//assert(plcGroupID.size() == plcMz.size() == plcElutionStart.size() == plcElutionEnd.size() == plcElutionPeak.size() );
#endif
#endif
		
	
		

		
	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	} // End loop through all intervals
	
	
	///////////////////////////////////////////////////////////////////////////
	// The original Matlab code has the following line here:
	// preLCLength=preLCLength+intervalCount;  
	// We do not need to return the specific preLCLength in a separate
	// variable since we are using C++ STL vector data structures, which
	// keep track of their own length.
	///////////////////////////////////////////////////////////////////////////
	
	
	///////////////////////////////////////////////////////////////////////////
	// For debugging purposes, output the test decision vector.
	// On a complex decision path, we want to be able to have visibility
	// into all decision points from a single location.
	///////////////////////////////////////////////////////////////////////////
	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	cout<<"decision vector="<<check1<<check2<<check3<<check4<<check5<<check6<<check7<<check8<<endl;
#endif
#endif
	
	
	
	
	return 0;  // return successfully
}



//////////////////////////////////////////////////////////////////////////////
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
//  4) elutionProfile.  This is a vector containing the peak data.
// 
//
// Output:
//  1) centerMass.  Reference to a double. Contains result value
//  2) maxOnLeft.  Reference to a bool. Contains indicator of whether peak
//                 is located on the left boundary.
//  3) maxOnRight.  Reference to a bool.  Contains indicator of whether peak
//                 is located on the right boundary.
// Return code: 	0 successful, 
// 					-1 dataMatrix.size() is zero,
//					-2 invalid mzGrid size
// 					-3 invalid mzGridStartID ( too large )
//                  -4 invalid mzGridStartID ( too small )
//                  -5 invalid dataMatrix row element
// 				    -6 invalid weightSum calculation
//   				-10  error at sgProc.matrixDiagMultiply call
//               	-20 error at mC.size <=0 
// 					-30 error at rc = sumColumns(mC,weightScan)
//  				-40 error at rc = sgProc.findMaxValue(weightScan,maxID, maxheight )
// 					-50 invalid mzGridCenterMassIndex
//
//
//
// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 02/15/12
//
///////////////////////////////////////////////////////////////////////////////
int Utilities::getLCPeakInfo(	vector< vector<double > > & dataMatrix, 
						vector<double>::size_type mzGridStartID, 
						vector<double> & mzGrid,
						vector<double> & elutionProfile,
						double & centerMass,
						bool & maxOnLeft,
						bool & maxOnRight) {
	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	cout<<"Inside Utiliites::getLCPeakInfo()"<<endl;
	cout<<"dataMatrix.size()="<<dataMatrix.size()<<endl;
	cout<<"mzGridStartID="<<mzGridStartID<<endl;
	cout<<"mzGrid.size()="<<mzGrid.size()<<endl;
	cout<<"elutionProfile.size()="<<elutionProfile.size()<<endl;
	cout<<"centerMass="<<centerMass<<endl;
	cout<<"maxOnLeft="<<maxOnLeft<<endl;
	cout<<"maxOnRight="<<maxOnRight<<endl;
#endif
#endif
	///////////////////////////////////////////////////////////////////////////
	// Local stack variables
	///////////////////////////////////////////////////////////////////////////
	vector<double> tempElutionProfile; // used to save a copy of the weight vector prior to normalizing it.
	vector<double> weight; // normalized row-wise vector
	double weightSum = 0; // row-wise summation of dataMatrix temporary  variable
	int rc = 0; // return code for function calls
	cbi::SignalProcessor sgProc;  // Signal processing routines
	
	///////////////////////////////////////////////////////////////////////////
	// check input parms
	///////////////////////////////////////////////////////////////////////////
	if(dataMatrix.size() == 0 ) {
		return -1; // invalid dataMatrix size
	}
	
	if(mzGrid.size() == 0 ) {
		return -2; // invalid mzGrid size
	}
	
	if (  mzGridStartID >= mzGrid.size() ) {
		return -3; // invalid mzGridStartID, out of range 
	}
	
	if( mzGridStartID < 0 ) {
		return -4;  // invalid mzGridStartID, out of range
	}
	
 	///////////////////////////////////////////////////////////////////////////
	//
	//  First make sure that we get a square dataset.
	//  Each row of data within the dataMatrix must have the same
	//  size as dataMatrix.size()
	//
	///////////////////////////////////////////////////////////////////////////
	for( vector<double>::size_type i = 0; i< dataMatrix.size(); i++ )  {
		if ( dataMatrix.at(i).size() != dataMatrix.size() ) {
			return -5; // invalid dataMatrix row element
		}
	}
	

	///////////////////////////////////////////////////////////////////////////
	// 
	// Calculate the sum of the dataMatrix.  The weight vector is a row-wise
	// summation of the dataMatrix.
	//
	///////////////////////////////////////////////////////////////////////////
	rc = sgProc.sumRows(dataMatrix,weight);
	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	cout<<"rc=sgProc.sum(dataMatrix,weight)"<<rc<<endl;
	cout<<"weight.size()="<<weight.size()<<endl;
	cout<<"weight.at(0)="<<weight.at(0)<<endl;
	cout<<"weight.at(1)="<<weight.at(1)<<endl;
#endif
#endif
	///////////////////////////////////////////////////////////////////////////
	//
	// Must keep a copy of the weight vector prior to making modifications
	// to it.  The weight vector is returned as the elutionProfile.
    // In general, we only ever want to modify output variables until
    // we are confident that the function completed successfully.
    //
	///////////////////////////////////////////////////////////////////////////
	tempElutionProfile = weight;  // Perform a deep copy, vector assignment 
	
	
	///////////////////////////////////////////////////////////////////////////
	// 
	// Calculate the sum of the weight vector.  The weight vector is a row-wise
	// summation of the dataMatrix.
	//
	///////////////////////////////////////////////////////////////////////////
	rc = sgProc.sum(weight,weightSum);
	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	cout<<"rc=sgProc.sum(weight,weightSum)"<<rc<<endl;
	cout<<"weightSum="<<weightSum<<endl;
#endif
#endif
	///////////////////////////////////////////////////////////////////////////
	//
	// 
	// Find the centerMass
	// 
	//
	///////////////////////////////////////////////////////////////////////////
	if( weightSum == 0 ){
		return -6;  // indicate an error condition to the caller
	}
	
	
	///////////////////////////////////////////////////////////////////////////
	//
	// 
	// Normalize the weight vector by performing an element by element 
	// division by the weightSum.
	// 
	//
	///////////////////////////////////////////////////////////////////////////
	for( vector<double>::size_type i = 0; i< weight.size(); i++ )  {
		weight.at(i) = weight.at(i) / weightSum;
	}
	
	///////////////////////////////////////////////////////////////////////////
	// End matrix library development testing
	///////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////
	//
	// 
	// Next, we need to find the diagonal of the dataMatrix.
	// The diagonal is a 1d signal.  
	// We then need to perform matrix multiplication 
	//   dataMatrx = rand(7,4)
	//		0.2511    0.9172    0.0540    0.0119
    //		0.6160    0.2858    0.5308    0.3371
    //		0.4733    0.7572    0.7792    0.1622
    //		0.3517    0.7537    0.9340    0.7943
    //		0.8308    0.3804    0.1299    0.3112
    //		0.5853    0.5678    0.5688    0.5285
    //		0.5497    0.0759    0.4694    0.1656
	// 
	//  weight = sum(dataMatrix,2) 
	//      1.2341
    //      1.7698
    //      2.1718
    //      2.8337
    //		1.6524
    //		2.2504
    //		1.2606
	//
	//
	//  diag(weight),  expands weight vector into a 2d matrix
	//      1.2341    0         0         0         0         0         0
    //		0    1.7698         0         0         0         0         0
    //		0         0    2.1718         0         0         0         0
    //		0         0         0    2.8337         0         0         0
    //		0         0         0         0    1.6524         0         0
    //		0         0         0         0         0    2.2504         0
    //		0         0         0         0         0         0    1.2606
	// 
	// diag(weight)*dataMatrix, performs matrix multiplication of the 
	//    7x7 diag(weight) and the 7x4 dataMatrix
	//
	// The result is a 7x4 matrix containing 
	//     0.3099    1.1319    0.0666    0.0147
    //     1.0903    0.5059    0.9394    0.5966
    //     1.0279    1.6445    1.6922    0.3522
    //     0.9965    2.1358    2.6467    2.2508
    //     1.3729    0.6286    0.2147    0.5143
    //     1.3171    1.2778    1.2801    1.1894
    //     0.6930    0.0956    0.5917    0.2088
	//
	//
	//
	//
	//  7x7  by 7x4 = 7x4 matrix
	//
	//
	//
	//
	///////////////////////////////////////////////////////////////////////////
	//vector< double > vA;  // vector A
	//vector< vector<double> > mA(7,vector<double>(7,0)); // matrix A
	//vector< vector<double> > mB(7,vector<double>(4,0)); ; // matrix B
	//vector< vector<double> > mC(7,vector<double>(4,0)); // matrix C
	
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Setup the test matrices
	//
	///////////////////////////////////////////////////////////////////////////
	//double count= 0;
	//for( vector<double>::size_type i = 0; i< 7; i++) {
	//	for( vector<double>::size_type j = 0; j < 7; j++ ) {
	//		mA.at(i).at(j) = count; 
	//		count++;
	//	}
	//}
	//
	//count= 0;
	//for( vector<double>::size_type i = 0; i< 7; i++) {
	//	for( vector<double>::size_type j = 0; j < 4; j++ ) {
	//		mB.at(i).at(j) = count; 
	//		count++;
	//	}
	//}
	//
	//for ( vector<double>::size_type i = 0; i < 7; i++ ) {
	//	vA.push_back(i+1);
	//}
	
	///////////////////////////////////////////////////////////////////////////
	//   Perform a diagonalization and matrix multiplication
	//   The notation used 
	///////////////////////////////////////////////////////////////////////////
	//rc = sgProc.matrixDiagMultiply(vA,mB,mC,0); 
	//rc = sgProc.matrixMultiply(mA, mB, mC, 0);
	//cout<<"sgProc.matrixMultiply, rc = "<<rc<<endl;
	//
	//for( vector<double>::size_type i = 0; i< 7; i++) {
	//	for( vector<double>::size_type j = 0; j < 4; j++ ) {
	//		cout<<mC.at(i).at(j)<<",";
	//	}
	//	cout<<endl;
	//}
	//
	//rc = sgProc.matrixMultiply(mA, mB, mC, 1);
	//cout<<"sgProc.matrixMultiply, rc = "<<rc<<endl;
	//
	//for( vector<double>::size_type i = 0; i< 7; i++) {
	//	for( vector<double>::size_type j = 0; j < 4; j++ ) {
	//		cout<<mC.at(i).at(j)<<",";
	//	}
	//	cout<<endl;
	//}
	//
	//rc = sgProc.matrixDiagMultiply(vA, mB, mC, 0);
	//cout<<"sgProc.matrixDiagMultiply, rc = "<<rc<<endl;
	//cout<<"mC.size()= "<<mC.size()<<endl;
	//cout<<"mC.at(0).size()="<<mC.at(0).size()<<endl;
	//
	//for( vector<double>::size_type i = 0; i< 7; i++) {
	//	for( vector<double>::size_type j = 0; j < 4; j++ ) {
	//		cout<<mC.at(i).at(j)<<",";
	//	}
	//	cout<<endl;
	//}
	//
	//rc = sgProc.matrixDiagMultiply(vA, mB, mC, 1);
	//cout<<"sgProc.matrixDiagMultiply, rc = "<<rc<<endl;
	//cout<<"mC.size()= "<<mC.size()<<endl;
	//cout<<"mC.at(0).size()="<<mC.at(0).size()<<endl;
	//
	//for( vector<double>::size_type i = 0; i< 7; i++) {
	//	for( vector<double>::size_type j = 0; j < 4; j++ ) {
	//		cout<<mC.at(i).at(j)<<",";
	//	}
	//	cout<<endl;
	//}	
	///////////////////////////////////////////////////////////////////////////
	// End matrix library development testing
	///////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Create a temporary data buffer to hold the result matrix
	//
	///////////////////////////////////////////////////////////////////////////
	vector< vector<double> > mC(weight.size(),vector<double>(dataMatrix.at(0).size(),0)); // matrix C
	
	///////////////////////////////////////////////////////////////////////////
	//
	//  Perform matrix multiplication with a auto-expanded diagonal
	//  vector.
	// 
	// 
	//  We can change to different matrix multiplcation libraries, by
	//  using the CBI interface, and simply changing the library selection
	//  indicator ( 4th parameter )
	//
	//
	///////////////////////////////////////////////////////////////////////////
	rc = sgProc.matrixDiagMultiply(weight, dataMatrix, mC, 0);
	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	cout<<"Utilities::getLCPeakInfo, sgProc.matrixDiagMultiply, rc = "<<rc<<endl;
#endif
#endif
	
	if ( rc != 0 ) {
		return -10; // indicate error at sgProc.matrixDiagMultiply call
	}
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Add up the data in the result of the matrix multiplication above
	// along the first dimension.  
	// This would be adding up all the data items along each column.
	// Therefore, the resultant, data would be a vector with size equal to 
	// the number of columns in mC.
	// 
	// 
	///////////////////////////////////////////////////////////////////////////
	if( mC.size() <= 0 ) {
		return -20; // indicate error at mC.size <=0 
	}
	
	///////////////////////////////////////////////////////////////////////////
	//
	//
	// Loop through all columns of the mC
	//  1 2 3 4                 0.02      0.02  0    0
	//  5 6 7 8 <------   diag( 2.05 )  = 0     2.05 0  x dataMatrix
	//  9 0 1 2                 4.45      0     4.45 0 
	// ==========              weight        
	// 15 8 11 14  weightScan
	//
	//
	// Performance notes,  since C/C++ stores data in row-major format,
	// we need to perform this column-wise summation one row at a time,
	// otherwise we would be having to touch memory locations that are potentially
	// very far away, within the same loop iteration, which would cause
	// cache line misses and subsequent pipeline efficiency problems, with 
	// the processor pipelines having to get flushed, and restarted.
	// 
	//
	///////////////////////////////////////////////////////////////////////////
	vector<double> weightScan;  // create a temp vector to hold the result of the column
	rc = sgProc.sumColumns(mC,weightScan); 
	
	if ( rc != 0 ) {
		return -30;  // indicate error at rc = sumColumns(mC,weightScan)
	}
	
	///////////////////////////////////////////////////////////////////////////
	// Each row must have the same number of data elements.
	///////////////////////////////////////////////////////////////////////////
	if ( weightScan.size() != mC.at(0).size() ) {
		return -40; // indicate invalid sumColumns result, difference # of columns 
	}
	
	///////////////////////////////////////////////////////////////////////////
	//
	//  Find the maximum of the weightScan vector, returning the 
	//  index of the maximum value.
	//      [maxheight,maxID]=max(weightScan);
	// 
	///////////////////////////////////////////////////////////////////////////
	vector<double>::size_type maxID;
	double maxHeight;
	
	rc = sgProc.findMaxValue(weightScan,maxID, maxHeight );
	if ( rc != 0 ) {
		return -40;  // indicate error at rc = sgProc.findMaxValue(weightScan,maxID, maxheight )
	}
	
	///////////////////////////////////////////////////////////////////////////
	// Find the center mass of the peak
	///////////////////////////////////////////////////////////////////////////
	vector<double>::size_type mzGridCenterMassIndex = mzGridStartID + maxID - 1;
	
	///////////////////////////////////////////////////////////////////////////
	// Make sure that the index of the center mass is valid
	///////////////////////////////////////////////////////////////////////////
	if ( (mzGridCenterMassIndex >= mzGrid.size()) || (mzGridCenterMassIndex < 0) ) {
		return -50; // invalid mzGridCenterMassIndex
	}
	
	
	///////////////////////////////////////////////////////////////////////////
	//  Copy results to output data structures
	//  Note:  Since we are porting this algorithm from Matlab, the 
	//  indexing needs to be carefully tested, to ensure that the 
	//  Matlab starting index format of 1 matches the C/C++ indexing
	//  scheme starting from 
	///////////////////////////////////////////////////////////////////////////
	centerMass = mzGrid.at(mzGridCenterMassIndex); 
	
	///////////////////////////////////////////////////////////////////////////
	//
	//  Provide an indication as to whether the maximum value of the 
	//  weights is located on the boundary.
	//  This is the purpose of the maxOnLeft and maxOnRight variables.
	//  If they are returned as 1, this indicates that the maximum 
	//  values of the weightScan vector was at either the left or right
	//  boundary.
	//
	///////////////////////////////////////////////////////////////////////////
	if ( maxID == 0 ) {
		maxOnLeft = true;
	}
	///////////////////////////////////////////////////////////////////////////
	// Since dataMatrix must be a square matrix, we only need to know
	// the number of rows to be able to figure out the size of the square.
	///////////////////////////////////////////////////////////////////////////
	if ( maxID == dataMatrix.size() ) {
		maxOnRight = true;
	}
	
	
	return 0; // return successfully
}


///////////////////////////////////////////////////////////////////////////////
// Name:  getInterval
// 
// Used to generate a set of intervals of the xic signal sorted by the
// height of each peak.
// Each interval should always only have 1 peak.
//
// A valid threshold value must always be provided to this function.  The 
// default value of 20 should be used, however it is the caller of this
// functions responsibility to set a valid threshold value.
// 
//
// Output: 
// Return code: 0 successful, ..., tbd
//
// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 02/15/12
///////////////////////////////////////////////////////////////////////////////
int Utilities::getInterval( vector <double> & xic, 
							double threshold, 
							int mininterval, 
							vector<int> & intervalListStart, 
                            vector<int> & intervalListEnd ) {
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Local variable
	//
	///////////////////////////////////////////////////////////////////////////
	int rc = 0;  // return code holder
	
	///////////////////////////////////////////////////////////////////////////
	//
	//
	// Create an interval map the same length as the input xic vector
	//
	//
	///////////////////////////////////////////////////////////////////////////
	if ( xic.size() <= 0 ) {
		return -1;  // invlaid xic input size
	}
	///////////////////////////////////////////////////////////////////////////
	// Check for invalid threshold, set to a default value
	///////////////////////////////////////////////////////////////////////////
	if( threshold < 0 ) {
		return -2; // invalid threshold value
	}
	
    vector<int> intervalMap(xic.size()+1,0);  // create local vector to store the cutpoints( bug fix, 10/17/12, add 1 to the size )
    vector<int> intervalMapDiff(xic.size()+1,0); // create local vector contain the diff of the intervalMap, bug fix, 10/17/12, add 1 to the size )
	
	///////////////////////////////////////////////////////////////////////////
	//
	//  Go through the xic vector, finding all those data items that have
	//  a value above the threshold.
	//
	//  For example:
	//  This would say that the positions where the values are set to 1
	//  exceeded the threshold value.
	//  00000000000000111111111100000000000011111111111111000000000000000000
	//    < threshold  >threshold <threshold    >threshold    < threshold
	// 
	///////////////////////////////////////////////////////////////////////////
	for( vector<double>::size_type i = 0; i < xic.size(); i++ ) {
		if ( xic.at(i) > threshold ) {
			intervalMap.at(i) = 1;
		}
	}
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Calculate the sum of the interval map
	//
	///////////////////////////////////////////////////////////////////////////
	cbi::SignalProcessor sgProc;  // Signal processing routines
	
    int sumIntervalMap = 0; // temporary variable to hold the summation of the threshold resultant vector
	rc = sgProc.sum(intervalMap,sumIntervalMap);  // evaluate the sum
	
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	cout<<"rc = sgProc.sum(intervalMap,sumIntervalMap),rc=,sumIntervalMap= "<<rc<<","<<sumIntervalMap<<endl;
#endif
#endif
	///////////////////////////////////////////////////////////////////////////
	//  If the sum of the intervalMap is zero, we know that there
	//  are no regions of the XIC that meet the threshold criteria and so we 
	//  want to return an empty 
	///////////////////////////////////////////////////////////////////////////
	vector<int> tempIntervalListStart;  // create empty interval list vector ( starting values )
	vector<int> tempIntervalListEnd;    // create empty interval list vector ( ending values )
	vector<int> tempIntervalListStartEmpty;  // create empty interval list vector ( starting values )
	vector<int> tempIntervalListEndEmpty;    // create empty interval list vector ( ending values )
    vector<int> finalIntervalListStart;  // create empty interval list vector ( starting values )
    vector<int> finalIntervalListEnd;    // create empty interval list vector ( ending values )
	
	///////////////////////////////////////////////////////////////////////////
	// 
	// No data to process
	//
	///////////////////////////////////////////////////////////////////////////
	if ( sumIntervalMap == 0 ) {
		intervalListStart = tempIntervalListStartEmpty; // deep copy of empty vector
		intervalListEnd =   tempIntervalListEndEmpty; // deep copy of empty vector
		///////////////////////////////////////////////////////////////////////
		//  Make sure that we return empty vectors for the interval list
		//  so that the caller of this method can check their sizes to find
		//  if there were any valid intervals found.
		///////////////////////////////////////////////////////////////////////
		return 0;
	}
	///////////////////////////////////////////////////////////////////////////
	// 
	// End - No data to process
	//
	///////////////////////////////////////////////////////////////////////////
	
	
	///////////////////////////////////////////////////////////////////////////
	//
	// If we get here, we know there is at least some data that meets the
    // threshold criteria.
	// Next we want to find the difference in the map intervalmap(i)-
    // intervalmap(i-1)
	//
	///////////////////////////////////////////////////////////////////////////
    rc = getDiff(intervalMap, -1, intervalMapDiff);  // evaluate the i - (i-1) difference  ( length changed by default to ll+1 from Matlab code )
	
	if ( rc < 0 ) {
		return -3; // getDiff error
    }
	

	///////////////////////////////////////////////////////////////////////////
	//
	//  Loop through all the items in the intervalMap, using the 
	//  information in the intervalMapDiff vector to make the decisions
	//  as to when to add a new interval to the intervalListStart and 
	//  intervalListEnd data structures.
    //
    //  Use a state machine to manage the detection of the interval
    //  found pattern: This is a much for concise and efficient
    //  approach than using many conditionals.
    //
    //  A valid interval consists of a 1 0000 -1 with the intermediate
    //  zeros.
    //
    //  This state machine allows for the detection of only valid
    //  intervals.
    //
    //
    //  State machine design by:nrz
	///////////////////////////////////////////////////////////////////////////

    typedef enum{ S0 = 0, S1 = 1, S2 = 2, S3= 3 } state;

    state currentState=S0;
    state nextState=S0;


    //set the current state;
    currentState = S0;

    int intervalStartIndex = 0;
    int intervalEndIndex = 0;

     int intervalCount = 0;
     for ( vector<int>::size_type i = 0; i < intervalMapDiff.size(); i++ )  {
           int dataValue = intervalMapDiff[i];
           if (intervalMapDiff[i] == 1){
               intervalCount++;
               tempIntervalListStart.push_back(i);
           }
           if( intervalMapDiff[i] == -1) {
               tempIntervalListEnd.push_back(i-1);
           }
     }

    //state transition loop
    //note: intervalMapDiff data starts at index 1 since
    //the diff operation is a forward difference.
//    for ( vector<int>::size_type i = 0; i < intervalMapDiff.size(); i++ )  {  // change i = 1 to i = 0, 10/17/12,bug fix to account for Matlab starting at 1 instead of 0.
//        // analyze the intervalMapDiff vector
//        int dataValue = intervalMapDiff[i];
//        switch(currentState) {
//        case S0:
//            switch(dataValue) {
//            case 0: nextState = S0; break;
//            case 1: nextState = S1; break;
//            case -1: nextState = S0; break;
//            default:  nextState = S0; break;// should never happen
//            }
//            break;
//        case S1:
//            // keep a record of the index at this state ( start of the interval )
//            intervalStartIndex  = i-1;  // bug fix, change from i to i-1 to match Matlab version( 10/17/12 )
//            //cout<<"Interval start found:"<<intervalStartIndex<<" ";
//            switch(dataValue) {
//            case 0: nextState = S2; break;
//            case 1: nextState = S1; break;
//            case -1: nextState = S3; break;  // accept state
//            default:  nextState = S0; break;// should never happen
//            }
//            break;
//        case S2:
//            switch(dataValue) {
//            case 0: nextState = S2; break;
//            case 1: nextState = S1; break;
//            case -1: nextState = S3; break;  // accept state
//            default:  nextState = S0; break;// should never happen
//            }
//            break;
//        case S3:
//            // interval was found state
//            intervalEndIndex = i-2; // bug fix, change from i-1 to i-2 to match Matlab version( 10/17/12 )
//            //cout<<"Interval found:"<<intervalStartIndex<<","<<intervalEndIndex<<endl;
//            // only when we find a valid interval do we
//            tempIntervalListStart.push_back(intervalStartIndex);
//            tempIntervalListEnd.push_back(intervalEndIndex);

//            switch(dataValue) {
//            case 0: nextState = S0; break;
//            case 1: nextState = S1; break;
//            case -1: nextState = S0; break;
//            default:  nextState = S0; break;// should never happen
//            }
//            break;
//        default: return -4;  // should never happen
//        }
//        // transition the state machine
//        //cout<<intervalMapDiff[i]<<" ";
//        currentState = nextState;

//    }
    //cout<<endl;
    ///////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////
	//
	//  Handle the case where there is at least 1 interval found,
	//  we need to make sure to handle whether or not there is an ending 
	//  value found, that is, we need each and every starting value to 
	//  have a corresponding end value.
	// 
	//  This is a very important check, that both of these match, 
	//  if they do not, this implies an error condition.
	//
	///////////////////////////////////////////////////////////////////////////
    if ( tempIntervalListStart.size() != tempIntervalListEnd.size() ) {
        return -4; // error, non matching interval start and end sizes
    }
		
	///////////////////////////////////////////////////////////////////////////
	//
	// Filter the decision,  we need to check each XIC interval, and make
	// a decision as to whether or not we want to include it in the final
	// list.
	// Loop through all intervals
	//
	///////////////////////////////////////////////////////////////////////////
    int intervalLength = 0; // length of each interval
    for ( vector<int>::size_type i = 0; i < tempIntervalListEnd.size(); i++ ) {
        intervalLength = tempIntervalListEnd.at(i) - tempIntervalListStart.at(i) +1;  // calculate interval length
        ///////////////////////////////////////////////////////////////////////
        // If the interval meets length criteria, keep in final
        // interval list.
        ///////////////////////////////////////////////////////////////////////
        if( intervalLength > mininterval ) {
            //cout<<"mininterval="<<mininterval<<endl;
            finalIntervalListStart.push_back(tempIntervalListStart.at(i));
            finalIntervalListEnd.push_back(tempIntervalListEnd.at(i));
        }
    }


	///////////////////////////////////////////////////////////////////////////
	// 
	// Generate final set of intervals
	// 
	///////////////////////////////////////////////////////////////////////////
    if ( finalIntervalListStart.size() > 0 ) {
        intervalListStart = finalIntervalListStart;  // deep copy of final interval start list
        intervalListEnd = finalIntervalListEnd;  // deep copy of final interval end list
    }
    else {
        intervalListStart = tempIntervalListStartEmpty; // deep copy of empty vector
        intervalListEnd =   tempIntervalListEndEmpty; // deep copy of empty vector
    }

	///////////////////////////////////////////////////////////////////////////
	// 
	// Return successfully
	// 
	///////////////////////////////////////////////////////////////////////////
	return 0;
}



///////////////////////////////////////////////////////////////////////////////
// Name:  splitInterval
// The purpose of this method is to perform additional cutting of an
// xic interval.
// Each interval is split to at most two intervals.
//
// Algorithm:
// 1) Find zero crossings using the diff signal
//	 2)  Find the valleys 
// 			- valleys have neighbors 
//    	3) If there is more than one valley, then there are places we 
//   	   need to split it.
//    	
//
// Output: 
// Return code: 0 successful, 
//
// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 02/20/12
///////////////////////////////////////////////////////////////////////////////
int Utilities::splitInterval(	vector <double> & inSmoothXic, 
                                vector <double> & inDiffXic,
                                vector<int> & inIntervalListStart,
                                vector<int> & inIntervalListEnd,
                                int & inMinLCLength,
                                vector<int> & outIntervalListStart,
                                vector<int> & outIntervalListEnd ) {

#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
    cout<<"Inside Utilities::splitInterval method"<<endl;
    cout<<"inSmoothXic.size()="<<inSmoothXic.size()<<endl;
    cout<<"inDiffXic.size()="<<inDiffXic.size()<<endl;
    cout<<"inIntervalListStart.size()="<<inIntervalListStart.size()<<",inIntervalListStart.at(0)"<<inIntervalListStart.at(0)<<endl;
    cout<<"inIntervalListEnd.size()="<<inIntervalListEnd.size()<<",inIntervalListEnd.at(0)"<<inIntervalListEnd.at(0)<<endl;
    cout<<"inMinLCLength="<<inMinLCLength<<endl;
#endif
#endif
#endif
    ///////////////////////////////////////////////////////////////////////////
    // Instantiate a math utilities object
    ///////////////////////////////////////////////////////////////////////////
    cbi::MathUtils muObject;
    cbi::SignalProcessor spObject;

    ///////////////////////////////////////////////////////////////////////////
    // Make sure the sizes of the interval list start and end values
    // match.
    ///////////////////////////////////////////////////////////////////////////
    if( inIntervalListStart.size() != inIntervalListEnd.size()) {
        cout<<"intervalstart and end size mismatch"<<endl;
        return -1; // mismatching input interval list start and end sizes
    }
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
    cout<<"splitinterval cp1"<<endl;
#endif
#endif
#endif
    ///////////////////////////////////////////////////////////////////////////
    //
    // Matlab newIntervalIndex variable maps to i variable in this program.
    //
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // Decision variables
    //  We need to be able to easily keep track of what decisions
    //  have been made and where in any function.  Specially
    //  for complex decision points with many composite checks, we need
    //  to be able to output a vector at the end of the function
    ///////////////////////////////////////////////////////////////////////////
    bool check1 = false;
    bool check2 = false;
    bool check3 = false;

    ///////////////////////////////////////////////////////////////////////////
    //
    // Create a copy of the input intervals
    //
    ///////////////////////////////////////////////////////////////////////////
    vector<int> emptyIntervalListStart;  // empty vector for output
    vector<int> emptyIntervalListEnd;    // empty vector for output
    vector<int> inIntervalListStartTemp		= 	inIntervalListStart; // temporary buffer output interval list starting points
    vector<int> inIntervalListEndTemp		=	inIntervalListEnd; // temporary buffer output interval list ending points
    vector<int> newIntervalListStart;
    vector<int> newIntervalListEnd;


    ///////////////////////////////////////////////////////////////////////////
    // BLOCK 0
    // Local variable
    // Loop through all the intervals
    //
    //  intervalListStart 		intervalListEnd
    //   scan#						scan#
    //   scan#						scan#
    // 	 scan#						scan#
    //
    // For each interval, we need to have a local set of data structures
    // per interval.
    //
    // Note:  Everything indexed by zCPid needs to be recreated for
    //        each iteration processed.
    //  1) zeroCrossingPoints
    //  2) zCPheight
    //  3) valleyPoints  ( a subset of the zeroCrossingPoints )
    //
    ///////////////////////////////////////////////////////////////////////////

#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
    cout<<"inIntervalListStart.size()="<<inIntervalListStart.size()<<endl;
    cout<<"splitinterval cp2"<<endl;
#endif
#endif
#endif



    for ( vector<int>::size_type i = 0; i < inIntervalListStart.size(); i++ ) {

        ///////////////////////////////////////////////////////////////////////
        // LOOP LOCAL VARIABLES
        // Make temporary local variables to be used as working buffers
        // for the data to be outputted.
        // For each interval we're processing,
        ///////////////////////////////////////////////////////////////////////
        vector<int> zeroCrossingPoints;  // contains the set of indices for each zero crossing point
        vector<double> zCPHeight;// contains the set of signal heights for each zero crossing point
        vector<int> valleyPoints; // contains the set of indices for each valley point
        int valleyPointCount=0;
        double heightv=0;  // contains the maximum value of the signal in each output interval
        // This is used to sort the intervals by maximum height
        int intervalLength = 0; // re-initialize the intervalLength variable (important to re-initialize every iteration )

        ///////////////////////////////////////////////////////////////////////
        //
        // END LOOP LOCAL VARIABLES
        //
        ///////////////////////////////////////////////////////////////////////


        ///////////////////////////////////////////////////////////////////////
        //
        // BLOCK 1
        //
        // Calculating interval length
        //
        //  Matlab arrays:  1 2 3 4 5 6 7 8 9   10
        //  C++ arrays:     0 1 2 3 4 5 6 7 8 9 10
        //
        //  Matlab interval: 5:9 : length = 9-5+1 = 4+1 = 5
        //  C++ interval   : 4:8 : length = 8-4+1 = 4+1 = 5
        //
        //  This is an important note, that for a given input interval
        //  List, the Matlab version assumes a starting point of 1, and the
        //  C++ version must assume a starting point of 0.
        //
        //  The intervalLength must be the same for both the Matlab
        //  version and the C++ version.
        //
        ///////////////////////////////////////////////////////////////////////
        intervalLength = inIntervalListEnd.at(i) - inIntervalListStart.at(i) +1;  // calculate interval length

#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
        cout<<"intervalLength="<<intervalLength<<endl;
#endif
#endif
#endif


        ///////////////////////////////////////////////////////////////////////
        // Loop local variables
        // Important, these vectors will be reset each time through the
        // loop.  It is therefore important to make a copy of the important
        // data in these temporary loop scope vectors prior to leaving the
        // current iteration.
        // zeroCrossingPoints and zCPheight have been declared above.
        //
        // Note:  The zero crossing point index was set to zCPid = 1
        // in the Matlab code, however we must set it to 0 in the C++ code,
        // since C++ arrays start at 0.
        //
        ///////////////////////////////////////////////////////////////////////
        int zCPid = 0; // zero crossing point index.

        ///////////////////////////////////////////////////////////////////////
        //
        // BLOCK 2
        //   We want to place the starting value of the current interval,
        //   into the zeroCrossingPoints data structure.
        //   From the diff process, +1 = a zero crossing point, and
        //   -1 = a zero crossing point, but one is at a place of a
        //   positive signal slope(+1) and the other is at a place of
        //   negative(-1) signal slope.
        //
        ///////////////////////////////////////////////////////////////////////


        ///////////////////////////////////////////////////////////////////////
        // Start with the smoothed XIC data, starting at the current
        // intervals starting index.  Place this intensity
        // data within the zCPheight vector.
        ///////////////////////////////////////////////////////////////////////

#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
        cout<<"splitinterval cp3"<<endl;
#endif
#endif
#endif


        ///////////////////////////////////////////////////////////////////////
        //
        //  BLOCK 3.  We want to count the total number of zero
        //            crossing points.
        //
        //  Loop through the indices within each interval
        //  For example:
        //  intervalListStart 		intervalListEnd
        //   scan#						scan#  -->
        //   scan#						scan#  -->
        //   scan#   					scan#  -->
        //
        // We need to be very careful when calculating the starting and
        // ending indices.  In Matlab, arrays start at 1, however in C++,
        // since arrays start at 0, we need to modifiy the range
        // calculations accordingly.
        //
        // Matlab source:
        //  for tempid=intervalList(intervalIndex,1)+1:intervalList(intervalIndex,2)-1
        ///////////////////////////////////////////////////////////////////////
        vector<int>::size_type startIndex=0;  // use for zero crossing loop counting
        vector<int>::size_type endIndex=0;    // use for zero crossing loop counting
        startIndex = inIntervalListStart.at(i);  // intervals input are 0 based arrays
        endIndex = inIntervalListEnd.at(i);  // intervals input are 0 based arrays

        ///////////////////////////////////////////////////////////////////////
        //  Zero Crossing Point Checking:
        //  Count the total number of zero crossing points.
        //  This loop must be inclusive [startIndex endIndex]
        ///////////////////////////////////////////////////////////////////////

#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
        cout<<"splitinterval cp4"<<endl;
#endif
#endif
#endif

        ///////////////////////////////////////////////////////////////////////
        //         2 cases:  There was a problem in the port,  issue handling
        //                   these 2 cases.  Going through the Matlab code
        //                   again revealed there are really only 2 cases:
        //                   ( This logic is the current interpretation of the
        //                     Matlab code )
        //          Case 1: There are no zero crossing points.
        //            We want to return the startIndex and endIndex as the as the
        //            set of zero crossing points.
        //                   p0 = startIndex
        //                   p1 = endIndex
        //          Case 2: There is at least one zero crossing point.
        //            We want to append the endIndex to the list of zero crossing
        //            points
        //                   p0 = first zero crossing point detected by loop ( not startIndex )
        //                        Note: startIndex is getting overwritten in Matlab version
        //
        //
        ///////////////////////////////////////////////////////////////////////

#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
        cout<<"tempid="<<startIndex+1<<" to "<<endIndex-1<<" incluse at the ends"<<endl;
#endif
#endif
#endif

        vector<int> zeroCrossingPointsTemp;
        vector<double> zCPHeightTemp;


        for ( vector<int>::size_type j = startIndex+1; j <= endIndex-1; j++ ) {  // bug fix 10/18/12
            check1 = ((muObject.sign(inDiffXic.at(j)) + muObject.sign(inDiffXic.at(j-1))) == 0 ); // cases 1,8,9
            check2 = (inDiffXic.at(j) == 0); // cases 6,7
            check3 = (check1 || check2);   // cases 1,6,7,8,9

#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
            cout<<"splitinterval cp5"<<endl;
            cout<<"check1="<<check1<<" check2="<<check2<<" check3="<<check3<<endl;
#endif
#endif
#endif


            if ( check3 ) {
                zeroCrossingPointsTemp.push_back(j); // keep the zero crossing points index
                zCPHeightTemp.push_back(inSmoothXic.at(j));// keep the height of the zero crossing point, bug fix 10/18/12
                zCPid++; // increase the zero crossing point counter

#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
                cout<<"splitinterval cp6"<<endl;
#endif
#endif
#endif

            }  // end if check3
        } // end zero crossing point loop

        // check to see if there were any zero crossing points(excluding the endpoints)
        if ( zCPid > 0 ){
            zeroCrossingPoints = zeroCrossingPointsTemp;
            zCPHeight = zCPHeightTemp;

            // There was at least 1 zero-crossing point found

            // To match the Matlab code, we will want to
            // 1) Append the ending index of the interval
            zeroCrossingPoints.push_back(endIndex); // keep the zero crossing points index
            zCPHeight.push_back(inSmoothXic.at(endIndex));// keep the height of the zero crossing point, bug fix 10/18/12
        }



#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
        cout<<"zeroCrossingPointsTemp.size()="<<zeroCrossingPoints.size()<<endl;
        for ( int idx = 0; idx < zeroCrossingPoints.size(); idx++  ){
            cout<<"zeroCrossingPoints.at(idx)="<<zeroCrossingPoints.at(idx)<<endl;
            cout<<"zCPHeight.at(idx)="<<zCPHeight.at(idx)<<endl;
        }
#endif
#endif
#endif


        // make sure there are at least 2 zero crossing points
        if( zeroCrossingPoints.size() > 1 ){
            // Find those zero crossing points that are valleys
            // The first one cannot be a valley

#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
            cout<<"zeroCrossingPoints.size()="<<zeroCrossingPoints.size()<<endl;
#endif
#endif
#endif


            for ( int idx = 1; idx < zeroCrossingPoints.size()-1; idx++ ){

#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
                cout<<"idx= "<<idx<<endl;
#endif
#endif
#endif

                bool check1 = (zCPHeight.at(idx) < zCPHeight.at(idx-1));
                bool check2 = (zCPHeight.at(idx) < zCPHeight.at(idx+1));
                if(  check1  && check2  ) {
                    valleyPointCount++;
                    valleyPoints.push_back(zeroCrossingPoints.at(idx));
                }
            }

            // if we found any valleys
            if ( valleyPoints.size() > 0 ) {

                vector<int> intervalListStartBuffer = inIntervalListStartTemp;
                vector<int> intervalListEndBuffer = inIntervalListEndTemp;
                vector<int> newIntervalListStartBuffer;
                vector<int> newIntervalListEndBuffer;

#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
                cout<<"Found one or more valleys"<<" valleyPoints.size()="<<valleyPoints.size()<<endl;
                for ( int idx2 = 0; idx2 < valleyPoints.size(); idx2++ ){
                    cout<<"valleyPoints.at(idx2)="<<valleyPoints.at(idx2)<<endl;
                }
#endif
#endif
#endif


                // handle the first valley point
                intervalListEndBuffer.at(i) = valleyPoints.at(0);
                newIntervalListStartBuffer.push_back(valleyPoints.at(0)+1);

                // handle the internal valley points( make sure there are at Least 2 points )
                if ( valleyPoints.size() > 1 ){

                    for ( int vpindex = 1; vpindex < valleyPoints.size(); vpindex++ ){
                        newIntervalListEndBuffer.push_back(valleyPoints.at(vpindex));
                        newIntervalListStartBuffer.push_back(valleyPoints.at(vpindex)+1);
                    }
                }

                // handle the last valley point
                newIntervalListEndBuffer.push_back(endIndex);

                // save the temporary buffer of intervals

                assert( newIntervalListStartBuffer.size() == newIntervalListEndBuffer.size() );

                // save the temporary buffers
                inIntervalListStartTemp = intervalListStartBuffer;
                inIntervalListEndTemp = intervalListEndBuffer; // the endpoint of the current interval has been modified

                ///////////////////////////////////////////////////////////////
                // Bug fix for missing intervals:
                // We are only keeping the first part of a cut interval
                // Refer to Utilites_test_10.txt unit test.
                // - Instead of overwriting the newIntervalListStart
                //   and newIntervalListEnd, we need to save the
                //   new intervals.
                ///////////////////////////////////////////////////////////////
                // Remove:               newIntervalListStart = newIntervalListStartBuffer;
                // Remove:               newIntervalListEnd = newIntervalListEndBuffer;
                // The following loop fixes the issue of missing intervals:
                // Test 9 continues to pass.
                // Start 10/22/12 @ 7:16 p.m. bug fix
                for ( int idx3 = 0; idx3 <  newIntervalListStartBuffer.size(); idx3++ ){
                    newIntervalListStart.push_back(newIntervalListStartBuffer.at(idx3));
                    newIntervalListEnd.push_back(newIntervalListEndBuffer.at(idx3));
                }
                // End 10/22/12 @ 7:16 p.m. bug fix


#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
                cout<<"inIntervalListStartTemp.size()="<<inIntervalListStartTemp.size()<<endl;
                cout<<"inIntervalListEndTemp.size()="<<inIntervalListEndTemp.size()<<endl;
                for ( int idx3 = 0; idx3 < inIntervalListStartTemp.size(); idx3++ ){
                    cout<<"inIntervalListStartTemp.at(idx3)="<<inIntervalListStartTemp.at(idx3)<<endl;
                    cout<<"inIntervalListEndTemp.at(idx3)="<<inIntervalListEndTemp.at(idx3)<<endl;
                }

                cout<<"newIntervalListStart.size()="<<newIntervalListStart.size()<<endl;
                cout<<"newIntervalListEnd.size()="<<newIntervalListEnd.size()<<endl;
                for ( int idx3 = 0; idx3 < newIntervalListStart.size(); idx3++ ){
                    cout<<"newIntervalListStart.at(idx3)="<<newIntervalListStart.at(idx3)<<endl;
                    cout<<"newIntervalListEnd.at(idx3)="<<newIntervalListEnd.at(idx3)<<endl;
                }
#endif
#endif
#endif


            }  // if check to make sure there are enough valley points

        } // if check to make sure that there is at least 2 zero crossing points

    }// end loop through intervals

#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
    cout<<"newIntervalListStart.size()="<<newIntervalListStart.size()<<endl;
    cout<<"newIntervalListEnd.size()="<<newIntervalListEnd.size()<<endl;
    for ( int idx3 = 0; idx3 < newIntervalListStart.size(); idx3++ ){
        cout<<"newIntervalListStart.at(idx3)="<<newIntervalListStart.at(idx3)<<endl;
        cout<<"newIntervalListEnd.at(idx3)="<<newIntervalListEnd.at(idx3)<<endl;
    }
#endif
#endif
#endif



    // Block 14.
    // concatenate the two lists
    for ( int idx3 = 0; idx3 < newIntervalListStart.size(); idx3++ ){
        inIntervalListStartTemp.push_back(newIntervalListStart.at(idx3));
        inIntervalListEndTemp.push_back(newIntervalListEnd.at(idx3));
    }

    // remove intervals that are too short
    vector<int> intervalLengthCheckBufferStart;
    vector<int> intervalLengthCheckBufferEnd;

    for ( int idx3 = 0; idx3 < inIntervalListStartTemp.size(); idx3++ ){
        int intervalLength = inIntervalListEndTemp.at(idx3)-inIntervalListStartTemp.at(idx3)+1;
        if ( intervalLength > inMinLCLength ) {
            intervalLengthCheckBufferStart.push_back(inIntervalListStartTemp.at(idx3));
            intervalLengthCheckBufferEnd.push_back(inIntervalListEndTemp.at(idx3));
        }
    }

    // update the intervals
    inIntervalListStartTemp = intervalLengthCheckBufferStart;
    inIntervalListEndTemp = intervalLengthCheckBufferEnd;


#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
    cout<<"updated interval list="<<endl;
    cout<<"inIntervalListStartTemp.size()="<<inIntervalListStartTemp.size()<<endl;
    cout<<"inIntervalListEndTemp.size()="<<inIntervalListEndTemp.size()<<endl;
    for ( int idx3 = 0; idx3 < inIntervalListStartTemp.size(); idx3++ ){
        cout<<"inIntervalListStartTemp.at(idx3)="<<inIntervalListStartTemp.at(idx3)<<endl;
        cout<<"inIntervalListEndTemp.at(idx3)="<<inIntervalListEndTemp.at(idx3)<<endl;
    }
#endif
#endif
#endif


    // if no valid intervals, return empty vectors
    if ( inIntervalListStartTemp.size() == 0) {
        outIntervalListStart = emptyIntervalListStart;
        outIntervalListEnd = emptyIntervalListEnd;
        return 0;
    }

#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
    // At this point there is at least 1 valid interval
    cout<<"There is at least 1 valid interval"<<endl;
#endif
#endif
#endif


    // Block 20
    // sort the intervals based on the intensity at the peak of the interval
    // find the maximum intensity of the smoothed XIC signal
    // for each interval, find the peak intensity back in the smoothed XIC
    vector<double> heightv;
    for ( int i = 0; i < inIntervalListStartTemp.size(); i++ ) {

        double tempMaxHeight = 0;

        vector<double> intervalIntensitySignalBuffer;

        // out indices must be within array bounds
        assert( inIntervalListStartTemp.at(i) < inSmoothXic.size() );
        assert( inIntervalListEndTemp.at(i) < inSmoothXic.size() );

        for ( int j = inIntervalListStartTemp.at(i); j <= inIntervalListEndTemp.at(i); j++ ){
            intervalIntensitySignalBuffer.push_back(inSmoothXic.at(j));
        }

        spObject.findMaxValue( intervalIntensitySignalBuffer, tempMaxHeight);
        // keep a record of the maximum intensity value at each interval
        heightv.push_back(tempMaxHeight);
    }

#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
    cout<<"Maximum height finding process completed"<<endl;
    cout<<"heightv.size()="<<heightv.size()<<endl;
    for ( int idx3 = 0; idx3 < heightv.size(); idx3++ ){
        cout<<"heightv.at(idx3)="<<heightv.at(idx3)<<endl;
    }
#endif
#endif
#endif


    // Sort the intervals based on their intensity heights at their peaks
    vector<int> indexVector;
    for ( int i = 0; i < heightv.size(); i++ ){
        indexVector.push_back(i);
    }
    // use a descending sorter
    int rc = 0;
    rc = spObject.sortGetIndex(heightv, heightv, indexVector, false);


#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
    cout<<"sortGetIndex, rc ="<<rc<<endl;
#endif
#endif
#endif


    // the indexVector should now contain the indices of the
    // data items sorted by maximum intensity height
    // e.g. indexVector = 0 1 2 3 4 5
    //      heightv     = 100.5 29.1 51.12 500.16 98.1 85.3
    //
    // The resulting indexVector should be
    //      indexVector = 3 0 4 5 2 1
    //
    // Therefore, we can populate the output lists according
    // to the sorted indices.
    vector<int> sortedIntervalListStartBuffer;
    vector<int> sortedIntervalListEndBuffer;
    for ( int i = 0; i < indexVector.size(); i++ ){

#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
        cout<<"indexVector.at(i)="<<indexVector.at(i)<<endl;
#endif
#endif
#endif

        sortedIntervalListStartBuffer.push_back(inIntervalListStartTemp.at(indexVector.at(i)));
        sortedIntervalListEndBuffer.push_back(inIntervalListEndTemp.at(indexVector.at(i)));
    }
    // return results
    outIntervalListStart = sortedIntervalListStartBuffer;
    outIntervalListEnd = sortedIntervalListEndBuffer;


#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
#ifdef SPLITINTERVALDEBUG
    for ( int idx3 = 0; idx3 < outIntervalListStart.size(); idx3++ ){
        cout<<"outIntervalListStart.at(idx3)="<<outIntervalListStart.at(idx3)<<endl;
        cout<<"outIntervalListEnd.at(idx3)="<<outIntervalListEnd.at(idx3)<<endl;
    }
#endif
#endif
#endif

    return 0;

}  // splitInterval





///////////////////////////////////////////////////////////////////////////////
// Name:  resampledVector
//  The purpose of this algorithm is to resample a set of mz and 
//  intensity value pairs 
//
// Algorithm:
//   Uses linear interpolation.   It may distort the peaks.
//   The algorithms maps a list of <mz,intensity> pairs to 
//   a uniform mzGrid.
// 	 
//   
// Input:
//
//  1. mzVector:  A 1d vector of doubles ( mz units )
//  2. intensity: A 1d vector of doubles ( intensity units )
//  3. mzGrid: A 1d vector of doubles ( mz units )
// 
//
// Output: 
//  1. resampledVector: A 1d vector of doubles. ( intensity units ) 
// Return code: 0 successful, 
// 				-1: mismatch input vector size
//				-2: invalid input size
//				-3: out of bounds array access attempt
//				-4: invalid currentMZVectorID( currentMZVectorID > mzVector.size() )
// 
// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 02/20/12
// Original Matlab code by: Michelle.Zhang@utsa.edu
///////////////////////////////////////////////////////////////////////////////
int Utilities::resampledVector( 	vector<double> & mzVector,
									vector<double> & intensity,
									vector<double> & mzGrid,
									vector<double> & resampledVectorResult) {
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
	cout<<"Inside resampledVector method"<<endl;
#endif
#endif
	
	///////////////////////////////////////////////////////////////////////////
	//  We need to snap all ofthe input mz and intensity pairs to the 
	//  new grid.
	//  The resampledVectorResultTemp vector object needs to be the 
	//  same size as the mzGrid vector,
	///////////////////////////////////////////////////////////////////////////
	vector<double> resampledVectorResultTemp(mzGrid.size(),0); // temporary vector initialized to zero
	vector<double>::size_type currentMZVectorID = 0; // running count of the current mz vector id

	///////////////////////////////////////////////////////////////////////////
	// Input validity check
	///////////////////////////////////////////////////////////////////////////
	if ( mzVector.size() != intensity.size()  ) {
		return -1; // mismatch input vector size
	}

	///////////////////////////////////////////////////////////////////////////
	//
	// Input validity check
	// We need to make sure that our input contains at least 2 elements,
	// since the algorithm is a "look ahead by 1" algorithm.
	// 
	// It is better to return an error code than continue.
	// 
	//
	// 
	///////////////////////////////////////////////////////////////////////////
	if ( mzVector.size() <= 2 ) {
		return -2; // invalid input size
	}

	///////////////////////////////////////////////////////////////////////////
	// Loop through the mzGrid.  Avoiding the use of the find function
	// for efficiency. 
	//
	// For each element in the grid.
	// 
	///////////////////////////////////////////////////////////////////////////
	for ( vector<double>::size_type i = 0; i< mzGrid.size(); i++ ) {	
		
		
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
		cout<<"i="<<i<<endl;
#endif
#endif
		///////////////////////////////////////////////////////////////////////
		//  For each iteration, we need to ask the following questions:
		//  1) Is the current mz vector in range
		//  2) As long as the grid value is between the current one and the 
		//     next one.  
		///////////////////////////////////////////////////////////////////////
		
		///////////////////////////////////////////////////////////////////////
		//
		//  The question asked is: A safeguard.  If the current mz vector is in
		//  range.  Question 2:  As long as the grid value is between the
		//  current one and the next one.  Then
		//    If mzVectorLength = 10
		//      in Matlab:  1 2 3 4 5 6 7 8 9 10 
		//         currentMZVectorID < mzVectorLength = [1-->9]
		//      in C++   :  0 1 2 3 4 5 6 7 8 9
		//         currentMZVectorID < mzVectorLength --> needs to be updated
		//         to
		//         currentMZVectorID < mzVectorLength-1 = [0-->8]
		//       This is needed to match the same Matlab logic in C++ version.
		//       Since we reference the currentMZVectorID+1, we need 
		//       to make sure we do not extend past the limits of the 
		//       vector.
		//        This makes sure that we can always make a reference
		//        to the last element in the vector.
		//
		//
		//
		//  At this point we can be sure that the result from the following
		//  line will be valid mzVector.size() -1, will not generate 
		//  negative values.
		//     
		///////////////////////////////////////////////////////////////////////       
		      
		if ( currentMZVectorID < (mzVector.size() - 1) ) {  // we've checked to make sure that mzVector.size() > 0
			///////////////////////////////////////////////////////////////////
			//
			// Look through the mzVector until we find a crossing point of the 
			// grid and the vector.  Loop through the grid until the mzVector
			// hits the mzGrid from the left side.
			//  
			// Go through the current mz vector ( that is , the vector
			// that is being fitted to the grid ).
			//
			// Stop when the current mz vector exceeds the value of the 
			// current grid point at location i.
			//
			//	mzVector.at(currentMZVectorID)  +1  
 			//								 |   |    ( currentMZVectorID)--> moving -->
			//                               |   |    [[[Stop when d3 becomes equal to or 
			//                               \/  \/      greater than v3]]]
			//  mz vector             =  d1  d2  d3  d4  d5  d6
			//  currentMZVectorID     =  0   1   2   3   4   5
			//                       
			//  mzGrid.at(i)           = v1  v2  v3  v4  v5  v6
			//  Grid                 i = 0   1   2   3   4   5 ... 
			//                                   ^
			//                                   |   ( i is fixed during while )
			//                                   |     loop
			//                              mzGrid.at(2) = v3
			//                        
			//  
			//  As soon as the stopping condition is met, we want to 
			//  update the currentMZVectorID by incrementing it by 1, 
			//  that is, the current index into the mz vector will be pointing
			//  to the value that met the stopping criteria, in the example
			//  above it would be indexing the value d3.
			//     
			// 
			// 
			///////////////////////////////////////////////////////////////////
			
			///////////////////////////////////////////////////////////////////
			//
			// The goal of this loop is to find the next value of
			// currentMZVectorID.
			// 
			///////////////////////////////////////////////////////////////////
			
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
			cout<<"i="<<i<<"mzGrid.at("<<i<<")="<<mzGrid.at(i)<<","<<"mzVector.at("<<currentMZVectorID+1<<")="<<mzVector.at(currentMZVectorID+1)<<endl;
#endif
#endif

			while (  mzGrid.at(i)  >= mzVector.at(currentMZVectorID+1)  )  {
				
				
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
				cout<<"Inside while loop, i="<<i<<"mzGrid.at("<<i<<")="<<mzGrid.at(i)<<endl;
#endif
#endif
				///////////////////////////////////////////////////////////////
				// As soon as we hit the boundary, we want to annotate 
				// where the boundary is located.
				///////////////////////////////////////////////////////////////
				currentMZVectorID = currentMZVectorID + 1;
				
				///////////////////////////////////////////////////////////////
				//  Right boundary check
				// 
				//  This is a stopping criteria check.
				//
				//
				//
				//  Since the algorithm is a "look ahead by 1" algorithm, 
				//  we must be able to know when we've reached "one before
				//  the end". 
				// 
				//  Otherwise, we will have an array out of bounds exception,
				//  for checked vectors .at(), or a memory corruption
				//  problem / crash if [] operator is used.
				// 
				//  Remember: for STL vectors, [] is not ranged checked, 
				//  while .at() is range checked.  
				//
				//  By default, .at() should be used except when 
				//  absolutely necessary in certain performance sensitive
				//  loops.  However, this must be done with a great 
				//  deal of caution, since without range checks, the 
				//  program will very likely crash if an index goes past
				//  the end of the array.
				// 
				//  We need to use mzVector.size() -2 because we will
				//  be checking currentMZVectorID+1 in the next iteration
				//  of the while loop.  If we use -1, then this will
				//  generate an out-of-bounds exception, since C++ 
				//  arrays are indexed from:
				//  [0 .... .size()-1]
				//
				///////////////////////////////////////////////////////////////
				if( currentMZVectorID > (mzVector.size()-2) ) {
					///////////////////////////////////////////////////////////
					// 
					//
					//
					// Boundary condition
					//
					//
					// If we get in this loop we need to copy the contents
					// of our temporary resampled vector to the output 
					// vector and return with a successful return code.
					// 
					// This if check simply tells us when to stop.
					// 
					//
					//
					//
					//
					///////////////////////////////////////////////////////////
					resampledVectorResult = resampledVectorResultTemp;
					return 0;  // out of bounds array access intent
				}  // end right boundary check
			}  // end while loop going from the currentMZVectorID index to one before the end of the mzVector
			
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
			cout<<"testpoint2"<<endl;
#endif
#endif
			///////////////////////////////////////////////////////////////////
			//
			//  The previous while loop has performed the action of finding
			//  the next value of currentMZVectorID.  This means we are 
			//  now at the boundary of where the input mz vector hits
			//  the boundary of the mzGrid.
			// 
			//
			//
			//  If the grid value is greater than the left and smaller than 
			//  right value
			//
			///////////////////////////////////////////////////////////////////
			
			
			///////////////////////////////////////////////////////////////////
			//
			// Look through the mzVector until we find a crossing point of the 
			// grid and the vector.  Loop through the grid until the mzVector
			// hits the mzGrid from the left side.
			//  
			// Go through the current mz vector ( that is , the vector
			// that is being fitted to the grid ).
			//
			// Stop when the current mz vector exceeds the value of the 
			// current grid point at location i.
			//
			//	mzVector.at(currentMZVectorID)  +1  
			//								 |   |    ( currentMZVectorID)--> moving -->
			//                               |   |    [[[Stop when d3 becomes equal to or 
			//                               \/  \/      greater than v3]]]
			//  mz vector             =  d1  d2  d3  d4  d5  d6
			//  currentMZVectorID     =  0   1   2   3   4   5
			//                       
			//  mzGrid.at(i)           = v1  v2  v3  v4  v5  v6
			//  Grid                 i = 0   1   2   3   4   5 ... 
			//                                   ^
			//                                   |   ( i is fixed during while )
			//                                   |     loop
			//                              mzGrid.at(2) = v3
			//                        
			//  
			//  As soon as the stopping condition is met, we want to 
			//  update the currentMZVectorID by incrementing it by 1, 
			//  that is, the current index into the mz vector will be pointing
			//  to the value that met the stopping criteria, in the example
			//  above it would be indexing the value d3.
			//     
			//
			//  Afterwards, we need to perform a linear approximation to 
			//  enable the estimation of the intensity value that "should"
			//  be at the mzGrid point if estimated linearly, since the
			//  intensity that is tied to the mz value d3 will not be 
			//  the same as the intensity mapped to value v3 is linearly
			//  interpolated.
			//
			//  check1 = v3 >= d2
			//  check2 = v3 < d3
			//  check3 = check1 and check2
			//
			//  We need to make sure that v3 is between  [ d2 and d3 ) 
			//   d2  <=   v3  <   d3
			//  Inclusive on the left and non-inclusive on the right
			// 
			//  
			//  
			//     
			//  
			///////////////////////////////////////////////////////////////////
			bool check1 = mzGrid.at(i) >= mzVector.at(currentMZVectorID);
			bool check2 = mzGrid.at(i) <  mzVector.at(currentMZVectorID+1); 
			bool check3 = check1 && check2;
			///////////////////////////////////////////////////////////////////
			//
			// Boundary condition safeguard
			//  Making sure that the grid value is greater than the 
			//  left and smaller than the right value.
			//
			//
			///////////////////////////////////////////////////////////////////
			
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
			cout<<"testpoint3, check3="<<check3<<endl;
#endif
#endif
			if( check3 ) { 
				///////////////////////////////////////////////////////////////
				//  We need to make sure that this check has the same
				//  meaning in a zero based array language as in a 
				//  1 based language:
				// 
				//
				// 
				///////////////////////////////////////////////////////////////
				
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
				cout<<"testpoint4, currentMZVectorID="<<currentMZVectorID<<endl;
#endif
#endif
				
				if( currentMZVectorID < (mzVector.size()-1) ) {
					
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
					cout<<"testpoint5,mzVector.size()-1="<<mzVector.size()-1<<endl;
#endif
#endif
					///////////////////////////////////////////////////////////
					// Calculate the y value of the line segment connecting
					// the two endpoints.  The ratio of the height/width of 
					// the two triangles should match. 
					// The large triange and the triangle in the center of the
					// current region.
					// The mzGrid value's height is determined using triangle
					// congruency property:
					//						 ^ (y axis = intensity units) 
					//				tempVar2 |            +
					//						 |       (?)+ |	
					//		dy =        ? -->|        +   |
					//	tempVar3			 |      + |   |
					//				tempVar1 |    +--------
					//						 |        ^
					//                       |________|_________________-> ( x-axis = mz units ( currentMZVectorID ))
					//					          1   2   3  ( x axis label)
					//               tempVar1 = intensity.at(1)
					//               tempVar2 = intensity.at(3)
					//						1 = mzVector.at(currentMZVectorID)
					//						2 = mzVector.at(currentMZVectorID+1)
					//
					//				
					// The Grid=              _________________________-> ( x-axis = grid's mz units )
					//  						      i
					//
					// check1 makes sure that the mzGrid.at(i) >= mzVector.at(currentMZVectorID)
					//   - It makes sure that the i point in the grid axis is within the two mzVector points.
					//
					// From the prior diagrams above, we see that the mzGrid point
					// becomes situated between the two mzVector 
					//
					//  1 --> d2
					//  2 --> v3
					//  3 --> d3
					//
					//  
					// GOAL:
					//
					// We are trying to estimate the height of the intensity 
					// at the grid point (v3) with x axis label(2).
					// This is the entire reason behind requiring that it be between
					// [d2 and d3).
					//
					// Since we know that intensity values corresponding to 
					// d2 and d3, we can use the triangle congruence property to
					// estimate the intensity value mapped to the grid point (v3)
					// with axis label (2) above.
					// 
					// The height ( that is the intensity value ) at the grid 
					// point is then a simple calculation by the triangle congruence:
					//
					//  The larger triangle is formed by the mzVector mz and intensity values.
					//  
					//  The smaller triangle is formed by the first mzVector value, and 
					//  the grid mz value which should have been at the 
					//  center of the two mzVector values.
					//
					//  Large Triangle:
					//  tempVar2-tempVar1 ( delta intensity value )  = height1
					//  3-1 ( delta mz value of mzVector ) = width1
					//  
					//  Small Internal Triangle whose right vertex is the intensity
					//  value we are trying to estimate:
					// 
					//  x-tempVar1 (delta intensity value ) = height2,   x is the unknown intensity of the mzGrid point
					//  2-1 ( delta mz value of mzGrid-mzVector) = width2
					//
					//  width2       height2
					// --------  =   -------- 
					//  width1       height1
					//
					// height2 = (width2*height1)/width1
					//  
					//  height2 contains the variable x, that we do not know and want to
					//   estimate.  x is the estimated intensity at the mzGrid point that 
					//   was between the two mzVector mz values.
					//
					//  height2 = (height1*width2) / width1
					//  x-tempVar1 = (height1*width2) / width1 
					//  x = tempVar1 + ( (height1*width2)/width1   )
					//  x = intensity.at(1) +  ( (intensity.at(3)-intensity.at(1))) * (mzGrid.at(2)-mzVector.at(1)) )/ (mzVector(3)-mzVector(1))
					//
					//
					//
					// Triangle Congruence:
					//
					// Matlab Code:
					//
					//  resampledVector(mzGridID)=intensity(currentMZVectorID)+ (intensity(currentMZVectorID+1)-intensity(currentMZVectorID))* ...
                    //       ( mzGrid(mzGridID)-mzVector(currentMZVectorID))/(mzVector(currentMZVectorID+1)-mzVector(currentMZVectorID));
					//
					//
					//
					// 
					///////////////////////////////////////////////////////////
					
					double tempVar1 = intensity.at(currentMZVectorID);    // part of height1
					double tempVar2 = intensity.at(currentMZVectorID+1);  // part of height1
					double tempVar3 = tempVar2-tempVar1;     // part of height1
					double tempVar4 = mzGrid.at(i)-mzVector.at(currentMZVectorID);  // part of width2
					double tempVar5 = mzVector.at(currentMZVectorID+1) - mzVector.at(currentMZVectorID); // part of width1
					double tempVar6 = tempVar4 / tempVar5; // ratio of the widths ( the triangle bases, width2 /width1 )
					double tempVar7 = tempVar3*tempVar6;   // height1 * ( width2/width1)
					double tempVar8 = tempVar1 + tempVar7; // calculate linearly approximated small triangle height. 
					
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
					cout<<setprecision(12)<<"testpoint5, tempVar8 = "<<tempVar8<<endl;
#endif
#endif
					
					///////////////////////////////////////////////////////////
					// Update temporary results vector
					///////////////////////////////////////////////////////////
					resampledVectorResultTemp.at(i) = tempVar8;

				} // end if case
				///////////////////////////////////////////////////////////////
				// We need to subtract one from mzVector.size() in order
				// to avoid overrunning the memory boundaries ( an
				// issue with C++ versus Matlab arrays. )
				// Matlab starts indexing at 1, C++ starts at 0.
				///////////////////////////////////////////////////////////////
				else if ( currentMZVectorID == (mzVector.size()-1) ) {
					///////////////////////////////////////////////////////////
					//
					// If the its the last element, output result
					// and do not continue loop
					//
					///////////////////////////////////////////////////////////
					
					resampledVectorResultTemp.at(i) = intensity.at(currentMZVectorID);
					resampledVectorResult = resampledVectorResultTemp;
					return 0;  // don't continue
					
				}  // end else if case
				else {
					return -3; // invalid currentMZVectorID( currentMZVectorID > mzVector.size() )
				}  // end else case, invalid currentMZVectorID
				///////////////////////////////////////////////////////////////
				//
				// End - Boundary condition checking
				//
				///////////////////////////////////////////////////////////////
				
			} // end check3 
		}  // end if check
	}  // end for loop
	///////////////////////////////////////////////////////////////////////////
	//
	// Return successful
	// Only copy internal buffer to output if processing succeeded.
	//
	///////////////////////////////////////////////////////////////////////////
	resampledVectorResult = resampledVectorResultTemp; // copy to output
	return 0; // return successfully
}  // end resampledVector method


///////////////////////////////////////////////////////////////////////////////
// Name:  removeShortIntervals
//  The purpose of this algorithm is to remove intervals ranges that are
//  too short.
//
//
// Algorithm:
//
//   
// Input:
//	1.
//  2.
//  3.
// 
//
// Output: 
//  1. 
// 
//  Error codes:
// 				-1: 
//				-2: 
//				-3: 
//				-4: 
// 
// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 03/01/12
// Original Matlab code by: Michelle.Zhang@utsa.edu
///////////////////////////////////////////////////////////////////////////////
int Utilities::removeShortIntervals (		vector<int> & inIntervalListStart,
								vector<int> & inIntervalListEnd,
								int & inMinLCLength,
								vector<int> & outIntervalListStart,
								vector<int> & outIntervalListEnd ) {
	
	return 0; // return successfully
	
}
	


///////////////////////////////////////////////////////////////////////////////
// Name:  sortIntervals
//  The purpose of this algorithm is to sort the set of intervals
//
//
// Algorithm:
//
//   
// Input:
//	1.
//  2.
//  3.
// 
//
// Output: 
//  1. 
// 
//  Error codes:
// 				-1: 
//				-2: 
//				-3: 
//				-4: 
// 
// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 03/01/12
// Original Matlab code by: Michelle.Zhang@utsa.edu
///////////////////////////////////////////////////////////////////////////////
int Utilities::sortIntervals(	vector<int> & inIntervalListStart,
					vector<int> & inIntervalListEnd,
					int & inMinLCLength,
					vector<int> & outIntervalListStart,
					vector<int> & outIntervalListEnd ) {
	
	return 0; // return successfully
	
}
	





///////////////////////////////////////////////////////////////////////////////
// Name:  mergeIntervals
//  The purpose of this algorithm is to merge intervals that meet a 
//  certain criteria.
//  The criteria is yet to be specified.
//
// Algorithm:
//
//   
// Input:
//	1.
//  2.
//  3.
// 
//
// Output: 
//  1. 
// 
//  Error codes:
// 				-1: 
//				-2: 
//				-3: 
//				-4: 
// 
// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 03/02/12
// Original Matlab code by: Michelle.Zhang@utsa.edu
///////////////////////////////////////////////////////////////////////////////
int Utilities::mergeIntervals(	vector<int> & inIntervalListStart,
					vector<int> & inIntervalListEnd,
					vector<int> & outIntervalListStart,
					vector<int> & outIntervalListEnd ) {
	
	return 0; // return successfully
	
}



///////////////////////////////////////////////////////////////////////////////
// Name:  estimateElutionProfile
//
// Algorithm:
//  Provided the section of resampledRegion section correscponding to the LC
//  signal along with signal at the peak apex location.
//
//  The algorithm is as follows:
//  The data input to this algorithm
//
//  regionData:
//
//  -----------------------------------> mz
//  |
//  |  1 2 3 4 20 6 7 8 9
//  |  1 2 3 4 20 6 7 8 9
//  |  1 2 3 4 20 6 7 8 9
//  |  1 2 3 4 20 6 7 8 9
//  |
//  |
//  |
//  |
//  |
//  V rt
//
//  In this case, the signal at location
//
// Input:
//	1. regionData
//  2. peakSignal
//
// Output:
//  1. output: contains the estimatedElutionProfile
//
//  Error codes:
// 				-1: rows = 0
//				-2: mismatch in number of colums
//				-3: invalid peakSignal length
//
// Port from Matlab to C++ by: Nelson.Ramirez@utsa.edu, 02/21/13
// Original Matlab code by: Michelle.Zhang@utsa.edu/Nelson.Ramirez@utsa.edu
///////////////////////////////////////////////////////////////////////////////
int Utilities::estimateElutionProfile(	vector< vector<double> > & regionData,
                                        vector<double> & peakSignal,
                                        vector<double> & output) {

    // Find the number of rows
    // Find the number of columns
    int rows = regionData.size();
    assert( rows > 0 );
    if ( rows <= 0 ) { return -1; }

    int cols = regionData.at(0).size();
    // note: all columns must always be the same length

    // debug loop
    for ( int i = 0; i < rows; i++){
        if ( regionData.at(i).size() != cols ){
            return -2;
        }
    }
    // make sure peak mz index has valid length,
    // it must be equal to the number of rows
    if ( peakSignal.size() != rows ) {
        return -3;
    }

    // At this point: rows,cols,peakColumnIndex have been verified

    //cout<<"estimatedElutionProfile,rows="<<rows<<"cols="<<cols<<endl;
    if ( cols == 0 ){
        return -4;
    }
    // Note: the output refers to a vector object, the data contained
    // within the vector object will be updated by this function

    // Find the sum of the peak's elution profile
    double sumPeakElutionProfile = 0;
    for ( int j = 0; j < rows; j++ ){
        sumPeakElutionProfile = sumPeakElutionProfile + peakSignal.at(j);
    }

    // Create a column-wise weighted running sum of each elution profile in the regionData
    vector<double> estimatedElutionProfileBuffer(rows,0); // temporary buffer
    vector<double> estimatedElutionProfile(rows,0);     // final result, to be copied to output

    for( int j = 0; j < rows; j++ ){
        for ( int i = 0; i < cols; i++ ){
            estimatedElutionProfileBuffer.at(j) = estimatedElutionProfileBuffer.at(j) +  regionData.at(j).at(i);
        }
    }

    // The estimatedElutionProfileBuffer is used to calculate the quotient value
    double quotientValue = 0;
    for( int j = 0; j < rows; j++ ){
        quotientValue = quotientValue + estimatedElutionProfileBuffer.at(j);
    }

    if( quotientValue == 0 ){
        return -5;
    }
    // At this point:
    // var1 == lcResampledRegionSection
    // var2 == peakSignal
    // var3 == estimatedElutionProfileBuffer
    // var4 == quotientValue
    // var5 == sumPeakElutionProfile
    double scalingFactor = sumPeakElutionProfile/quotientValue;

    //cout<<"scalingFactor="<<scalingFactor<<endl;

    // Evaluate the final value of the estimatedElutionProfile
    // Loop through all the rows of the estimated elutionprofile
    for( int j = 0; j < rows; j++ ){
            estimatedElutionProfileBuffer.at(j) = estimatedElutionProfileBuffer.at(j) * scalingFactor;
//            cout<<"estimatedElutionProfile.at(j)="<<estimatedElutionProfileBuffer.at(j)<<endl;
    }


    // Perform a deep copy of the estimatedElutionProfile vector
    output = estimatedElutionProfileBuffer;

    return 0; // return successfully
}



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
int Utilities::generateParallelJobFiles(	cbi::ParameterLoader & paraObjGlobal ){
	cout<<"Inside Utilities::generateParallelJobFiles method"<<endl;


	///////////////////////////////////////////////////////////////////////
	// Get the mzXMLFilename
	// This is information that is used from the global parameter file to
	// generate localized information for each worker.
	///////////////////////////////////////////////////////////////////////
	vector<string> globalparameterFilename = paraObjGlobal.getParameterValuesAsString("globalparameterFilename"); // required parm
	vector<string> runDirectory = paraObjGlobal.getParameterValuesAsString("runDirectory"); // required parm
	vector<string> exportCommand = paraObjGlobal.getParameterValuesAsString("exportcommand"); // required parm
	vector<int> parGlobalMinScan = paraObjGlobal.getParameterValuesAsInt("parglobalminscan"); // required parm
	vector<int> parGlobalMaxScan = paraObjGlobal.getParameterValuesAsInt("parglobalmaxscan"); // required parm
	vector<double> parGlobalMinMz = paraObjGlobal.getParameterValuesAsDouble("parglobalminmz");		 // required parm
	vector<double> parGlobalMaxMz = paraObjGlobal.getParameterValuesAsDouble("parglobalmaxmz"); // required parm
	vector<double> parMzDelta = paraObjGlobal.getParameterValuesAsDouble("parmzdelta"); // required parm
	vector<double> parMzDeltaOverlap = paraObjGlobal.getParameterValuesAsDouble("parmzdeltaoverlap"); // required parm
	vector<double> mzBinningThresholdPPM = paraObjGlobal.getParameterValuesAsDouble("mzBinningThresholdPPM"); // required parm

	vector<int> numberOfMachines = paraObjGlobal.getParameterValuesAsInt("numberOfMachines"); // required parm
	vector<int> coresPerMachine = paraObjGlobal.getParameterValuesAsInt("coresPerMachine"); // required parm
	vector<string> jobscriptfilenamebase = paraObjGlobal.getParameterValuesAsString("jobscriptfilenamebase"); // required parm
	vector<string> jobscriptfilenameext = paraObjGlobal.getParameterValuesAsString("jobscriptfilenameext"); // required parm
	vector<string> localparameterfilenamebase = paraObjGlobal.getParameterValuesAsString("localparameterfilenamebase"); // required parm
	vector<string> localparameterfilenameext = paraObjGlobal.getParameterValuesAsString("localparameterfilenameext"); // required parm
	vector<string> shell = paraObjGlobal.getParameterValuesAsString("shell"); // required parm
	vector<string> queuesubmitcommand = paraObjGlobal.getParameterValuesAsString("queuesubmitcommand"); // optional parm
	vector<string> envvarcommand = paraObjGlobal.getParameterValuesAsString("envvarcommand"); // required parm
	vector<string> queueType = paraObjGlobal.getParameterValuesAsString("queueType"); //  optional parm
	vector<string> jobscriptoutputdirectory = paraObjGlobal.getParameterValuesAsString("jobscriptoutputdirectory"); // required parm
	vector<string> localparameterfileoutputdirectory = paraObjGlobal.getParameterValuesAsString("localparameterfileoutputdirectory"); // required parm
	vector<string> msdaexepath = paraObjGlobal.getParameterValuesAsString("msdaexepath"); // required parm
	vector<string> msdaexe = paraObjGlobal.getParameterValuesAsString("msdaexe"); // required parm
	vector<string> exportldlibrarypathlinux = paraObjGlobal.getParameterValuesAsString("exportldlibrarypathlinux"); // required parm
	vector<int> mpienabled = paraObjGlobal.getParameterValuesAsInt("mpienabled"); // required parm
	vector<string> mpicommand = paraObjGlobal.getParameterValuesAsString("mpicommand"); // optional parm
	vector<string> mpiprocessesflag = paraObjGlobal.getParameterValuesAsString("mpiprocessesflag"); // optional parm
	vector<int> nummpiprocesses = paraObjGlobal.getParameterValuesAsInt("nummpiprocesses"); // optional parm
	vector<string> mpinodesflag = paraObjGlobal.getParameterValuesAsString("mpinodesflag");  // optional parm
	vector<int> nummpinodes = paraObjGlobal.getParameterValuesAsInt("nummpinodes"); // optional parm

	///////////////////////////////////////////////////////////////////////////
	// Validate the input parms.
	// Additional validation should be added passed the "empty" check
	///////////////////////////////////////////////////////////////////////////
	if ( 	globalparameterFilename.size() == 0 ||
			runDirectory.size() == 0 ||
			exportCommand.size() == 0 ||
			parGlobalMinScan.size() == 0 ||
			parGlobalMaxScan.size() == 0 || 
			parGlobalMinMz.size()   == 0 ||
			parGlobalMaxMz.size()   == 0 ||
			parMzDelta.size()    	== 0 ||
			parMzDeltaOverlap.size()== 0 ||
			mzBinningThresholdPPM.size() == 0 ||
			numberOfMachines.size() == 0 ||
			coresPerMachine.size() == 0 ||
			jobscriptfilenamebase.size() == 0 ||
			jobscriptfilenameext.size() == 0 ||
			localparameterfilenamebase.size() == 0 ||
			localparameterfilenameext.size() == 0 ||
			shell.size() == 0 ||
			queuesubmitcommand.size() == 0 ||
			envvarcommand.size() == 0 ||
			queueType.size() == 0 ||
			jobscriptoutputdirectory.size() == 0 ||
			localparameterfileoutputdirectory.size() == 0 ||
			msdaexepath.size() == 0 ||
			msdaexe.size() == 0 ||
			exportldlibrarypathlinux.size() == 0 ||
			mpienabled.size() == 0
			) {
		return -1;  // required parameters are empty
	} 
	


#ifdef MSDADEBUG
	cout<<"runDirectory.at(0)="<<runDirectory.at(0)<<endl;
	cout<<"exportCommand.at(0)="<<exportCommand.at(0)<<endl;
	cout<<"parGlobalMinScan.at(0)="<<parGlobalMinScan.at(0)<<endl;
	cout<<"parGlobalMaxScan.at(0)="<<parGlobalMaxScan.at(0)<<endl;
	cout<<"parGlobalMinMz.at(0)="<<parGlobalMinMz.at(0)<<endl;
	cout<<"parGlobalMaxMz.at(0)="<<parGlobalMaxMz.at(0)<<endl;
	cout<<"parMzDelta.at(0)="<<parMzDelta.at(0)<<endl;
	cout<<"parMzDeltaOverlap.at(0)="<<parMzDeltaOverlap.at(0)<<endl;
	cout<<"mzBinningThresholdPPM.at(0)="<<mzBinningThresholdPPM.at(0)<<endl;
	cout<<"numberOfMachines.at(0)="<<numberOfMachines.at(0)<<endl;
	cout<<"coresPerMachine.at(0)=";
	for ( vector<int>::size_type i = 0; i < coresPerMachine.size(); i++ ) {
		cout<<coresPerMachine.at(i)<<" ";
	}
	cout<<endl;
	cout<<"jobscriptfilenamebase.at(0)="<<jobscriptfilenamebase.at(0)<<endl;
	cout<<"jobscriptfilenameext.at(0)="<<jobscriptfilenameext.at(0)<<endl;
	cout<<"localparameterfilenamebase.at(0)="<<localparameterfilenamebase.at(0)<<endl;
	cout<<"localparameterfilenameext.at(0)="<<localparameterfilenameext.at(0)<<endl;
	cout<<"shell.at(0)="<<shell.at(0)<<endl;
	cout<<"queuesubmitcommand.at(0)="<<queuesubmitcommand.at(0)<<endl;
	cout<<"envvarcommand.at(0)="<<envvarcommand.at(0)<<endl;
	cout<<"queueType.at(0)="<<queueType.at(0)<<endl;
	cout<<"jobscriptoutputdirectory.at(0)="<<jobscriptoutputdirectory.at(0)<<endl;
	cout<<"localparameterfileoutputdirectory.at(0)="<<localparameterfileoutputdirectory.at(0)<<endl;
	cout<<"msdaexepath.at(0)="<<msdaexepath.at(0)<<endl;
	cout<<"msdaexe.at(0)="<<msdaexe.at(0)<<endl;
	cout<<"exportldlibrarypathlinux.at(0)="<<exportldlibrarypathlinux.at(0)<<endl;
#endif

	
	///////////////////////////////////////////////////////////////////////////
	// Check if MPI is enabled, if it is, this makes it required
	// to have a few additional parameters in the global
	// parameters file.
	///////////////////////////////////////////////////////////////////////////
	if( mpienabled.at(0) > 0 ) { // this implies mpi is enabled
		if( 	mpicommand.size() == 0 || 
				mpiprocessesflag.size() == 0 ||
				nummpiprocesses.size() == 0 ||
				mpinodesflag.size() == 0 ||
				nummpinodes.size() == 0  ) {
			return -2; // required mpi parameters are empty
		}
#ifdef MSDADEBUG
		cout<<"mpicommand.at(0)="<<mpicommand.at(0)<<endl;
		cout<<"mpiprocessesflag.at(0)="<<mpiprocessesflag.at(0)<<endl;
		cout<<"nummpiprocesses.at(0)="<<nummpiprocesses.at(0)<<endl;
		cout<<"mpinodesflag.at(0)="<<mpinodesflag.at(0)<<endl;
		cout<<"nummpinodes.at(0)="<<nummpinodes.at(0)<<endl;
#endif
	}

	
	
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
	
	string parmFilenameString = localparameterfileoutputdirectory.at(0).append(localparameterfilenamebase.at(0)); // parameter files
	double startMZValue = parGlobalMinMz.at(0);
	double endMZValue = parGlobalMaxMz.at(0);
	double mzDelta = parMzDelta.at(0);
	double mzOverlap = parMzDeltaOverlap.at(0);
	///////////////////////////////////////////////////////////////////////////
	// Generate the sequence of start and end mz values for each 
	// distributed job.
	///////////////////////////////////////////////////////////////////////////
	vector<double> startMZ;
	vector<double> endMZ;
	double counter = 0;
	for ( counter = startMZValue; counter <= endMZValue; counter = counter + parMzDelta.at(0) ) {
		///////////////////////////////////////////////////////////////////////
		//
		// For each combination, we need to create a new parameter file 
		// for this combination.
		//
		///////////////////////////////////////////////////////////////////////
		startMZ.push_back(counter);
		endMZ.push_back(counter+parMzDelta.at(0));
	}

	///////////////////////////////////////////////////////////////////////////
	// 
	///////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////
	// 
	// Add the overlap calculation later... 
	// A second pass over the mz intervals 
	//
	///////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////
	// 
	// Generate the parameter files followed by 
	// 
	///////////////////////////////////////////////////////////////////////////
		
	if ( startMZ.size() > 200 ) {
		return -1; // to many files will be created, we want to abort the process...
	}
	
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
		scangroupidTemp << setw(5) << setfill('0');
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
		pFilename = pFilename.append(localparameterfilenameext.at(0));
		cout<<pFilename<<endl;
		
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
		pFile<<"localworkerscansrange={"<<parGlobalMinScan.at(0)<<":"<<1<<":"<<parGlobalMaxScan.at(0)<<"}"<<endl;
		pFile<<"localMinMz="<<startMZ.at(i)<<endl;
		pFile<<"localMaxMz="<<endMZ.at(i)<<endl;
		pFile<<"mzBinningThresholdInPPM="<<mzBinningThresholdPPM.at(0)<<endl;
		pFile<<"scangroupid="<<scangroupid<<endl;
		///////////////////////////////////////////////////////////////////////
		// Close the file
		///////////////////////////////////////////////////////////////////////
		pFile.close();




		
	}  // end parameter file 
	

	///////////////////////////////////////////////////////////////////////////
	//
	// Create the job scripts that will manage the processing of the
	// local parameters files
	// 
	///////////////////////////////////////////////////////////////////////////
	string jobScriptFilenameString = jobscriptoutputdirectory.at(0).append(jobscriptfilenamebase.at(0)); // job script files
	
	
	for ( vector<double>::size_type i = 0 ; i < numberOfMachines.at(0); i++ ) {
#ifdef MSDADEBUG
#ifdef UTILITIESDEBUG
		cout<<"numberOfMachines.at(0)="<<numberOfMachines.at(0)<<endl;
		cout<<"jobScriptFilenameString="<<jobScriptFilenameString<<endl;
#endif
#endif

		///////////////////////////////////////////////////////////////////////
		//
		// Generate the output filenames
		// according to the corresponding scangroup id and width as specified
		// by the user.  Multiple scangroups could be processed by a single
		// worker, so the fundamental identifier of a unit of msda work is the 
		// scan group id.  
		// 
		///////////////////////////////////////////////////////////////////////
		int jobScriptId= i;  // must convert from an integer to a string
		string valString; // string to hold the integer with leading zeros.
		///////////////////////////////////////////////////////////////////////
		// Use a string stream to hold the converted integer to the 
		// leading zero based string.
		///////////////////////////////////////////////////////////////////////
		stringstream jobScriptIdTemp (stringstream::in | stringstream::out);
		///////////////////////////////////////////////////////////////////////
		//
		// Set the fill option so that files get named with leading 
		// zeros.
		//
		///////////////////////////////////////////////////////////////////////
		jobScriptIdTemp << setw(5) << setfill('0');
		///////////////////////////////////////////////////////////////////////
		//
		// Send the current scangroupid information to the string stream 
		// so that we can generate the identification string.
		//
		///////////////////////////////////////////////////////////////////////
		jobScriptIdTemp << jobScriptId;
		///////////////////////////////////////////////////////////////////////
		//
		// Move the string stream over to a holding string.
		//
		///////////////////////////////////////////////////////////////////////
		jobScriptIdTemp >> valString;

		///////////////////////////////////////////////////////////////////////
		//
		// Create the pre lc candidate filename
		//
		///////////////////////////////////////////////////////////////////////
		string jsFilename = jobScriptFilenameString;
		jsFilename = jsFilename.append(valString);
		jsFilename = jsFilename.append(jobscriptfilenameext.at(0));
		cout<<jsFilename<<endl;

		///////////////////////////////////////////////////////////////////////
		//
		// Create the internal data structure information file, each with
		// its scan group identifier.
		//
		///////////////////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////////////////
		//
		// Create the pre LC Candidate file
		// 
		///////////////////////////////////////////////////////////////////////
		ofstream jsFile; 	// the ofstream object owns the file.
		// it will handle cleanup if anything 
		// happens after the file structure
		// has been allocated.
		
		jsFile.open(jsFilename.c_str()); // open the file.  Notem pFile owns the file,
										// the ofstream object pFile is responsible
										// for managing the file resource, 
										// in case if anything happens 
		
		///////////////////////////////////////////////////////////////////////
		//
		// Write out the contents of the file
		// 
		///////////////////////////////////////////////////////////////////////
		jsFile<<exportCommand.at(0)<<" "<<exportldlibrarypathlinux.at(0)<<endl;
		jsFile<<"cd "<<runDirectory.at(0)<<endl;
		
		///////////////////////////////////////////////////////////////////////
		// Depending on the number of compute nodes and the number of 
		// cores on each compute node, generate calls to the msda 
		// executable accordingly.
		///////////////////////////////////////////////////////////////////////
		
		for ( vector<double>::size_type j = 0; j < startMZ.size(); j++ ){
			
			///////////////////////////////////////////////////////////////////
			// Generate the local parameter file ids 
			///////////////////////////////////////////////////////////////////
			int parmFileId = j; 
			string parmFileIdString; // contains final string representation
			stringstream parmFileIdStringStream (stringstream::in | stringstream::out);
			parmFileIdStringStream <<setw(5)<<setfill('0');
			parmFileIdStringStream << parmFileId;
			parmFileIdStringStream >> parmFileIdString;
			
			
			///////////////////////////////////////////////////////////////////
			//
			//
			//
			///////////////////////////////////////////////////////////////////
			jsFile<<mpicommand.at(0)<<" "<<mpiprocessesflag.at(0)<<" "<<nummpiprocesses.at(0)<<" "<<msdaexe.at(0)<<" "<<globalparameterFilename.at(0)<<" "<<parmFilenameString<<parmFileIdString<<localparameterfilenameext.at(0)<<endl;
			///////////////////////////////////////////////////////////////////
			//
			//
			//
			///////////////////////////////////////////////////////////////////
			
		}
		
		///////////////////////////////////////////////////////////////////////
		// Close the file
		///////////////////////////////////////////////////////////////////////
		jsFile.close();

	}  // end for loop generating job scripts
		

	return 0;
}

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
/// Nelson.Ramirez@utsa.edu, 04/16/12
///////////////////////////////////////////////////////////////////////////////
int Utilities::workflowOne( cbi::ParameterLoader & paraObjGlobal ) {

	///////////////////////////////////////////////////////////////////////////
	//
	// workflow 1
	//
	///////////////////////////////////////////////////////////////////////////
	cout<<"Inside Utilities::workflowOne"<<endl;


		///////////////////////////////////////////////////////////////////////
		// Get the mzXMLFilename
		// This is information that is used from the global parameter file to
		// generate localized information for each worker.
		// We need to minimize the number of code changes, when a parameter
		// name needs to be updated, therefore, the workflow# methods
		// should be the only place where changes should be needed when 
		// changing the name / value of a parameter.
		///////////////////////////////////////////////////////////////////////
		vector<string> globalparameterFilename = paraObjGlobal.getParameterValuesAsString("globalparameterFilename"); // required parm
		vector<string> runDirectory = paraObjGlobal.getParameterValuesAsString("runDirectory"); // required parm
		vector<string> exportCommand = paraObjGlobal.getParameterValuesAsString("exportcommand"); // required parm
		vector<int> parGlobalMinScan = paraObjGlobal.getParameterValuesAsInt("parglobalminscan"); // required parm
		vector<int> parGlobalMaxScan = paraObjGlobal.getParameterValuesAsInt("parglobalmaxscan"); // required parm
		vector<double> parGlobalMinMz = paraObjGlobal.getParameterValuesAsDouble("parglobalminmz");		 // required parm
		vector<double> parGlobalMaxMz = paraObjGlobal.getParameterValuesAsDouble("parglobalmaxmz"); // required parm
		vector<double> parMzDelta = paraObjGlobal.getParameterValuesAsDouble("parmzdelta"); // required parm
		vector<double> parMzDeltaOverlap = paraObjGlobal.getParameterValuesAsDouble("parmzdeltaoverlap"); // required parm
		vector<double> mzBinningThresholdPPM = paraObjGlobal.getParameterValuesAsDouble("mzBinningThresholdPPM"); // required parm

		vector<int> numberOfMachines = paraObjGlobal.getParameterValuesAsInt("numberOfMachines"); // required parm
		vector<int> coresPerMachine = paraObjGlobal.getParameterValuesAsInt("coresPerMachine"); // required parm
		vector<string> jobscriptfilenamebase = paraObjGlobal.getParameterValuesAsString("jobscriptfilenamebase"); // required parm
		vector<string> jobscriptfilenameext = paraObjGlobal.getParameterValuesAsString("jobscriptfilenameext"); // required parm
		vector<string> localparameterfilenamebase = paraObjGlobal.getParameterValuesAsString("localparameterfilenamebase"); // required parm
		vector<string> localparameterfilenameext = paraObjGlobal.getParameterValuesAsString("localparameterfilenameext"); // required parm
		vector<string> shell = paraObjGlobal.getParameterValuesAsString("shell"); // required parm
		vector<string> queuesubmitcommand = paraObjGlobal.getParameterValuesAsString("queuesubmitcommand"); // optional parm
		vector<string> envvarcommand = paraObjGlobal.getParameterValuesAsString("envvarcommand"); // required parm
		vector<string> queueType = paraObjGlobal.getParameterValuesAsString("queueType"); //  optional parm
		vector<string> jobscriptoutputdirectory = paraObjGlobal.getParameterValuesAsString("jobscriptoutputdirectory"); // required parm
		vector<string> localparameterfileoutputdirectory = paraObjGlobal.getParameterValuesAsString("localparameterfileoutputdirectory"); // required parm
		vector<string> msdaexepath = paraObjGlobal.getParameterValuesAsString("msdaexepath"); // required parm
		vector<string> msdaexe = paraObjGlobal.getParameterValuesAsString("msdaexe"); // required parm
		vector<string> exportldlibrarypathlinux = paraObjGlobal.getParameterValuesAsString("exportldlibrarypathlinux"); // required parm
		vector<int> mpienabled = paraObjGlobal.getParameterValuesAsInt("mpienabled"); // required parm
		vector<string> mpicommand = paraObjGlobal.getParameterValuesAsString("mpicommand"); // optional parm
		vector<string> mpiprocessesflag = paraObjGlobal.getParameterValuesAsString("mpiprocessesflag"); // optional parm
		vector<int> nummpiprocesses = paraObjGlobal.getParameterValuesAsInt("nummpiprocesses"); // optional parm
		vector<string> mpinodesflag = paraObjGlobal.getParameterValuesAsString("mpinodesflag");  // optional parm
		vector<int> nummpinodes = paraObjGlobal.getParameterValuesAsInt("nummpinodes"); // optional parm

		///////////////////////////////////////////////////////////////////////////
		// Validate the input parms.
		// Additional validation should be added passed the "empty" check
		///////////////////////////////////////////////////////////////////////////
		if ( 	globalparameterFilename.size() == 0 ||
				runDirectory.size() == 0 ||
				exportCommand.size() == 0 ||
				parGlobalMinScan.size() == 0 ||
				parGlobalMaxScan.size() == 0 || 
				parGlobalMinMz.size()   == 0 ||
				parGlobalMaxMz.size()   == 0 ||
				parMzDelta.size()    	== 0 ||
				parMzDeltaOverlap.size()== 0 ||
				mzBinningThresholdPPM.size() == 0 ||
				numberOfMachines.size() == 0 ||
				coresPerMachine.size() == 0 ||
				jobscriptfilenamebase.size() == 0 ||
				jobscriptfilenameext.size() == 0 ||
				localparameterfilenamebase.size() == 0 ||
				localparameterfilenameext.size() == 0 ||
				shell.size() == 0 ||
				queuesubmitcommand.size() == 0 ||
				envvarcommand.size() == 0 ||
				queueType.size() == 0 ||
				jobscriptoutputdirectory.size() == 0 ||
				localparameterfileoutputdirectory.size() == 0 ||
				msdaexepath.size() == 0 ||
				msdaexe.size() == 0 ||
				exportldlibrarypathlinux.size() == 0 ||
				mpienabled.size() == 0
				) {
			return -1;  // required parameters are empty
		} 
		


	#ifdef MSDADEBUG
		cout<<"runDirectory.at(0)="<<runDirectory.at(0)<<endl;
		cout<<"exportCommand.at(0)="<<exportCommand.at(0)<<endl;
		cout<<"parGlobalMinScan.at(0)="<<parGlobalMinScan.at(0)<<endl;
		cout<<"parGlobalMaxScan.at(0)="<<parGlobalMaxScan.at(0)<<endl;
		cout<<"parGlobalMinMz.at(0)="<<parGlobalMinMz.at(0)<<endl;
		cout<<"parGlobalMaxMz.at(0)="<<parGlobalMaxMz.at(0)<<endl;
		cout<<"parMzDelta.at(0)="<<parMzDelta.at(0)<<endl;
		cout<<"parMzDeltaOverlap.at(0)="<<parMzDeltaOverlap.at(0)<<endl;
		cout<<"mzBinningThresholdPPM.at(0)="<<mzBinningThresholdPPM.at(0)<<endl;
		cout<<"numberOfMachines.at(0)="<<numberOfMachines.at(0)<<endl;
		cout<<"coresPerMachine.at(0)=";
		for ( vector<int>::size_type i = 0; i < coresPerMachine.size(); i++ ) {
			cout<<coresPerMachine.at(i)<<" ";
		}
		cout<<endl;
		cout<<"jobscriptfilenamebase.at(0)="<<jobscriptfilenamebase.at(0)<<endl;
		cout<<"jobscriptfilenameext.at(0)="<<jobscriptfilenameext.at(0)<<endl;
		cout<<"localparameterfilenamebase.at(0)="<<localparameterfilenamebase.at(0)<<endl;
		cout<<"localparameterfilenameext.at(0)="<<localparameterfilenameext.at(0)<<endl;
		cout<<"shell.at(0)="<<shell.at(0)<<endl;
		cout<<"queuesubmitcommand.at(0)="<<queuesubmitcommand.at(0)<<endl;
		cout<<"envvarcommand.at(0)="<<envvarcommand.at(0)<<endl;
		cout<<"queueType.at(0)="<<queueType.at(0)<<endl;
		cout<<"jobscriptoutputdirectory.at(0)="<<jobscriptoutputdirectory.at(0)<<endl;
		cout<<"localparameterfileoutputdirectory.at(0)="<<localparameterfileoutputdirectory.at(0)<<endl;
		cout<<"msdaexepath.at(0)="<<msdaexepath.at(0)<<endl;
		cout<<"msdaexe.at(0)="<<msdaexe.at(0)<<endl;
		cout<<"exportldlibrarypathlinux.at(0)="<<exportldlibrarypathlinux.at(0)<<endl;
	#endif

	
		
		
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
		
		string parmFilenameString = localparameterfileoutputdirectory.at(0).append(localparameterfilenamebase.at(0)); // parameter files
		double startMZValue = parGlobalMinMz.at(0);
		double endMZValue = parGlobalMaxMz.at(0);
		double mzDelta = parMzDelta.at(0);
		double mzOverlap = parMzDeltaOverlap.at(0);
		///////////////////////////////////////////////////////////////////////////
		// Generate the sequence of start and end mz values for each 
		// distributed job.
		///////////////////////////////////////////////////////////////////////////
		vector<double> startMZ;
		vector<double> endMZ;
		double counter = 0;
		for ( counter = startMZValue; counter <= endMZValue; counter = counter + parMzDelta.at(0) ) {
			///////////////////////////////////////////////////////////////////////
			//
			// For each combination, we need to create a new parameter file 
			// for this combination.
			//
			///////////////////////////////////////////////////////////////////////
			startMZ.push_back(counter);
			endMZ.push_back(counter+parMzDelta.at(0));
		}

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
		// Generate the parameter files.  We need to have a limit
		// to the total number of files generated.  This should 
		// be a global paramter.
		//
		///////////////////////////////////////////////////////////////////////////
			
		if ( startMZ.size() > 200 ) {
			return -1; // to many files will be created, we want to abort the process...
		}
		
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
			scangroupidTemp << setw(5) << setfill('0');
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
			pFilename = pFilename.append(localparameterfilenameext.at(0));
			cout<<pFilename<<endl;
			
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
			pFile<<"localworkerscansrange={"<<parGlobalMinScan.at(0)<<":"<<1<<":"<<parGlobalMaxScan.at(0)<<"}"<<endl;
			pFile<<"localMinMz="<<startMZ.at(i)<<endl;
			pFile<<"localMaxMz="<<endMZ.at(i)<<endl;
			pFile<<"mzBinningThresholdInPPM="<<mzBinningThresholdPPM.at(0)<<endl;
			pFile<<"scangroupid="<<scangroupid<<endl;
			///////////////////////////////////////////////////////////////////////
			// Close the file
			///////////////////////////////////////////////////////////////////////
			pFile.close();




			
		}  // end parameter file 
		

		///////////////////////////////////////////////////////////////////////////
		//
		// Create the job scripts that will manage the processing of the
		// local parameters files
		// 
		///////////////////////////////////////////////////////////////////////////
		string jobScriptFilenameString = jobscriptoutputdirectory.at(0).append(jobscriptfilenamebase.at(0)); // job script files
		
		
		for ( vector<double>::size_type i = 0 ; i < numberOfMachines.at(0); i++ ) {
	#ifdef MSDADEBUG
	#ifdef UTILITIESDEBUG
			cout<<"numberOfMachines.at(0)="<<numberOfMachines.at(0)<<endl;
			cout<<"jobScriptFilenameString="<<jobScriptFilenameString<<endl;
	#endif
	#endif

			///////////////////////////////////////////////////////////////////////
			//
			// Generate the output filenames
			// according to the corresponding scangroup id and width as specified
			// by the user.  Multiple scangroups could be processed by a single
			// worker, so the fundamental identifier of a unit of msda work is the 
			// scan group id.  
			// 
			///////////////////////////////////////////////////////////////////////
			int jobScriptId= i;  // must convert from an integer to a string
			string valString; // string to hold the integer with leading zeros.
			///////////////////////////////////////////////////////////////////////
			// Use a string stream to hold the converted integer to the 
			// leading zero based string.
			///////////////////////////////////////////////////////////////////////
			stringstream jobScriptIdTemp (stringstream::in | stringstream::out);
			///////////////////////////////////////////////////////////////////////
			//
			// Set the fill option so that files get named with leading 
			// zeros.
			//
			///////////////////////////////////////////////////////////////////////
			jobScriptIdTemp << setw(5) << setfill('0');
			///////////////////////////////////////////////////////////////////////
			//
			// Send the current scangroupid information to the string stream 
			// so that we can generate the identification string.
			//
			///////////////////////////////////////////////////////////////////////
			jobScriptIdTemp << jobScriptId;
			///////////////////////////////////////////////////////////////////////
			//
			// Move the string stream over to a holding string.
			//
			///////////////////////////////////////////////////////////////////////
			jobScriptIdTemp >> valString;

			///////////////////////////////////////////////////////////////////////
			//
			// Create the pre lc candidate filename
			//
			///////////////////////////////////////////////////////////////////////
			string jsFilename = jobScriptFilenameString;
			jsFilename = jsFilename.append(valString);
			jsFilename = jsFilename.append(jobscriptfilenameext.at(0));
			cout<<jsFilename<<endl;

			///////////////////////////////////////////////////////////////////////
			//
			// Create the internal data structure information file, each with
			// its scan group identifier.
			//
			///////////////////////////////////////////////////////////////////////

			///////////////////////////////////////////////////////////////////////
			//
			// Create the pre LC Candidate file
			// 
			///////////////////////////////////////////////////////////////////////
			ofstream jsFile; 	// the ofstream object owns the file.
			// it will handle cleanup if anything 
			// happens after the file structure
			// has been allocated.
			
			jsFile.open(jsFilename.c_str()); // open the file.  Notem pFile owns the file,
											// the ofstream object pFile is responsible
											// for managing the file resource, 
											// in case if anything happens 
			
			///////////////////////////////////////////////////////////////////////
			//
			// Write out the contents of the file
			// 
			///////////////////////////////////////////////////////////////////////
			jsFile<<exportCommand.at(0)<<" "<<exportldlibrarypathlinux.at(0)<<endl;
			jsFile<<"cd "<<runDirectory.at(0)<<endl;
			
			///////////////////////////////////////////////////////////////////////
			// Depending on the number of compute nodes and the number of 
			// cores on each compute node, generate calls to the msda 
			// executable accordingly.
			///////////////////////////////////////////////////////////////////////
			
			for ( vector<double>::size_type j = 0; j < startMZ.size(); j++ ){
				
				///////////////////////////////////////////////////////////////////
				// Generate the local parameter file ids 
				///////////////////////////////////////////////////////////////////
				int parmFileId = j; 
				string parmFileIdString; // contains final string representation
				stringstream parmFileIdStringStream (stringstream::in | stringstream::out);
				parmFileIdStringStream <<setw(5)<<setfill('0');
				parmFileIdStringStream << parmFileId;
				parmFileIdStringStream >> parmFileIdString;
				
				
				///////////////////////////////////////////////////////////////////
				//
				//  Output the command string
				//
				///////////////////////////////////////////////////////////////////
				jsFile<<msdaexe.at(0)<<" "<<globalparameterFilename.at(0)<<" "<<parmFilenameString<<parmFileIdString<<localparameterfilenameext.at(0)<<endl;
				///////////////////////////////////////////////////////////////////
				//
				// End - Output the command string
				//
				///////////////////////////////////////////////////////////////////
				
			}
			
			///////////////////////////////////////////////////////////////////////
			// Close the file
			///////////////////////////////////////////////////////////////////////
			jsFile.close();

		}  // end for loop generating job scripts
			

	return 0;
}  // end workflowOne


///////////////////////////////////////////////////////////////////////////////
// Name:  workflowTwo
//  The purpose of method is to implement run mode #2.
//  Run mode # 2 contains all the code necessary to implement the data input
//  and data organization and binning process.
//  A dataset is loaded into memory via random file i/o access using 
//  the mzXML file format. 
//
//  This mode should normally be called by an autogenerated script.
//  At this point both the global parameter file and the local parameter
//  file must be available.
//
//  This mode generates MSDA Level 1 Format Data,the first interface
//  layer to the MSDA HPC proteomics data processing architecture.
//
//  MSDA Level 1 format takes a section of raw data and processes it
//  to generate a ScanGroup object.  The data within a scan group is
//  a set of pre-lc regions of interest. 
//
//  In the next processing stage, the LC Candidate generation algorithm
//  will be applied at the Pre-LC region of interest level. 
//
//
//  From each Pre-LC region of interest, a set of LC Candidates will
//  be generated.  
//
//
//
// Algorithm:
//   MSDA Data Organization Algorithm 
//   
// Input:
//
//  1. ParameterLoader object reference, containing global parameters
//  2. ParameterLoader object reference, containing local parameters
// 
// Output: 
//  1. A single output file in MSDA Level 1 data format.
//
//
// Return code: 0 successful, 
// 
// Nelson.Ramirez@utsa.edu, 04/16/12
///////////////////////////////////////////////////////////////////////////////

int Utilities::workflowTwo( cbi::ParameterLoader & paraObjGlobal, cbi::ParameterLoader & paraObjLocalWorker) {
	


	///////////////////////////////////////////////////////////////////////////
	//
	// workflow 2
	//
	///////////////////////////////////////////////////////////////////////////
	cout<<"Running workflow 2"<<endl;
	
	
	
	///////////////////////////////////////////////////////////////////////
	// Get the mzXMLFilename 
	///////////////////////////////////////////////////////////////////////
	vector<string> mzxmlfilename = paraObjGlobal.getParameterValuesAsString("mzXMLFilename");

	///////////////////////////////////////////////////////////////////////
	// Get the list of scan numbers that this executable instance
	// is responsible for.
	///////////////////////////////////////////////////////////////////////
	vector<int> scanList = paraObjLocalWorker.getParameterValuesAsIntRange("localworkerscansrange");
	
	///////////////////////////////////////////////////////////////////////
	// Get the identifier of the scan group this executable is working
	// on.
	///////////////////////////////////////////////////////////////////////
	vector<int> scanGroupId = paraObjLocalWorker.getParameterValuesAsInt("scangroupid");
	
	///////////////////////////////////////////////////////////////////////
	// Handle a scan list parameter error
	// Throw an exception.  ( Later change to a derived exception object )
	// There should be at least one scan to process, and the 
	// scan number must be a positive integer as required by the 
	// getParameterValuesAsIntRange api specification.
	///////////////////////////////////////////////////////////////////////
	if ( scanList.size() == 1 ) {  
		if ( scanList.at(0) < 0 ) {
			throw "Invalid scan list range exception";
		}
	}
	
	///////////////////////////////////////////////////////////////////////
	// Get the minimum MZ value we will be processing
	///////////////////////////////////////////////////////////////////////
	vector<double> minMz = paraObjGlobal.getParameterValuesAsDouble("minMz");
	
	///////////////////////////////////////////////////////////////////////
	// Get the maximum MZ value we will be processing
	///////////////////////////////////////////////////////////////////////
	vector<double> maxMz = paraObjGlobal.getParameterValuesAsDouble("maxMz");
	
	

	///////////////////////////////////////////////////////////////////////
	// Get the level of the ms scans to process.
	// MSDA handles level 1 scans only.
	// Level 1 scan are 
	///////////////////////////////////////////////////////////////////////
	vector<int> mslevel = paraObjGlobal.getParameterValuesAsInt("mslevelselector");
	
	///////////////////////////////////////////////////////////////////////
	// Get the basename for visualization output files
	// This is the full path and base of the file name
	// ../../filenamebase-######.plt
	// ../../filenamebase-######.dat
	// ../../filenamebase-######.csv ( MSDA Level 1 Format )
	///////////////////////////////////////////////////////////////////////
	vector<string> visfile = paraObjGlobal.getParameterValuesAsString("visfilenamebase");
			
	
	///////////////////////////////////////////////////////////////////////
	// Get the local worker minimum MZ value we will be processing
	///////////////////////////////////////////////////////////////////////
	vector<double> minMzLocal = paraObjLocalWorker.getParameterValuesAsDouble("localMinMz");
		
	
	///////////////////////////////////////////////////////////////////////
	// Get the local worker maximum MZ value we will be processing
	///////////////////////////////////////////////////////////////////////
	vector<double> maxMzLocal = paraObjLocalWorker.getParameterValuesAsDouble("localMaxMz");
		
	
	///////////////////////////////////////////////////////////////////////
	// Get the mzBinning threshold in PPM.
	// This threshold is used to cut the mz bins
	// This control value is the first critical user controlled
	// parameter.
	///////////////////////////////////////////////////////////////////////
	vector<double> mzBinningThresholdInPPMLocal = paraObjLocalWorker.getParameterValuesAsDouble("mzBinningThresholdInPPM");
			

	///////////////////////////////////////////////////////////////////////
	//
	// Load parameters related to LC Candidate Generation
	//
	///////////////////////////////////////////////////////////////////////
	vector<int> genlcParmNumSmoothPoint = paraObjGlobal.getParameterValuesAsInt("numSmoothPoint");
	vector<int> genlcParmMinLCLength= paraObjGlobal.getParameterValuesAsInt("minLCLength");
	vector<double> genlcParmNoiseThresholdLevel= paraObjGlobal.getParameterValuesAsDouble("noiseThresholdLevel");
	vector<double> genlcParmMassResolution= paraObjGlobal.getParameterValuesAsDouble("massResolution");
	vector<double> genlcParmR2Threshold= paraObjGlobal.getParameterValuesAsDouble("R2Threshold");
	

	///////////////////////////////////////////////////////////////////////
	// Use asserts to check validity of input parameters
	///////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////
	// Make sure that there is a valid filename
	///////////////////////////////////////////////////////////////////////
	assert(mzxmlfilename.size() > 0);
	
	///////////////////////////////////////////////////////////////////////
	// Make sure that there is a valid filename
	///////////////////////////////////////////////////////////////////////
	assert(scanList.size() > 0);
	
	
	///////////////////////////////////////////////////////////////////////
	// Make sure minMz and maxMz parameters are valid
	///////////////////////////////////////////////////////////////////////
	assert( minMz.at(0) <= maxMz.at(0) );
	
	
	///////////////////////////////////////////////////////////////////////
	// Make sure mzBinning Threshold parameter is valid
	///////////////////////////////////////////////////////////////////////
	assert( mzBinningThresholdInPPMLocal.size() > 0); // Make sure that this parameter was specified
	assert( mzBinningThresholdInPPMLocal.at(0) > 0 ); // Make sure that the value of the parameter is valid
	///////////////////////////////////////////////////////////////////////
	// Parameter error handling
	///////////////////////////////////////////////////////////////////////
	if ( minMz.at(0) > maxMz.at(0) ) {  
		throw "minMz > maxMz exception";
	}
	if ( minMz.at(0) <= 0 ) {  
		throw "minMz <= 0 exception";
	}
	if ( minMzLocal.at(0) > maxMzLocal.at(0) ) {  
		throw "minMzLocal > maxMzLocal exception";
	}
	if ( minMzLocal.at(0) <= 0 ) {  
		throw "minMzLocal <= 0 exception";
	}
	if ( mslevel.size() < 1 ) {  
		throw "mslevel.size() <= 1 exception";
	}
	if ( scanGroupId.size() < 1){
		throw "scanGroupId.size() < 1 exception";
	}
	if ( visfile.size() < 1){
		throw "visfile.size() < 1 exception";
	}
	if ( mzBinningThresholdInPPMLocal.size() < 1){
		throw "mzBinningThresholdInPPMLocal.size() < 1 exception";
	}
	

	if( genlcParmNumSmoothPoint.size() < 1 ) {
		throw "genlcParmNumSmoothPoint.size() < 1 exception";
	}

	if( genlcParmMinLCLength.size() < 1 ) {
		throw "genlcParmMinLCLength.size() < 1 exception";
	}
	
	if( genlcParmNoiseThresholdLevel.size() < 1 ) {
		throw "genlcParmNoiseThresholdLevel.size() < 1 exception";
	}
	
	if( genlcParmMassResolution.size() < 1 ) {
		throw "genlcParmMassResolution.size() < 1 exception";
	}

	if( genlcParmR2Threshold.size() < 1 ) {
		throw "genlcParmR2Threshold.size() < 1 exception";
	}

		
	

#ifdef MSDADEBUG
	cout<<"LCCandidateGenerationParameter, numSmoothPoint="<<genlcParmNumSmoothPoint.at(0)<<endl;
	cout<<"LCCandidateGenerationParameter, minLCLength="<<genlcParmMinLCLength.at(0)<<endl;
	cout<<"LCCandidateGenerationParameter, noiseThresholdLevel="<<genlcParmNoiseThresholdLevel.at(0)<<endl;
	cout<<"LCCandidateGenerationParameter, massResolution="<<genlcParmMassResolution.at(0)<<endl;
	cout<<"LCCandidateGenerationParameter, R2Threshold="<<genlcParmR2Threshold.at(0)<<endl;
#endif
	

	
	return 0;
}


///////////////////////////////////////////////////////////////////////////////
// Name:  runModeThree
//  The purpose of method is to implement run mode #2.
//  Run mode # 2 contains all the code necessary to implement the data input
//  and data organization and binning process.
//  A dataset is loaded into memory via random file i/o access using 
//  the mzXML file format. 
//
//  This mode should normally be called by an autogenerated script.
//  At this point both the global parameter file and the local parameter
//  file must be available.
//
//  This mode generates MSDA Level 1 Format Data.  The first interface
//  layer to the HPC proteomics data processing architecture.
//
//  MSDA Level 1 format takes a section of raw data and processes it
//  to generate a ScanGroup object.  The data within a scan group is
//  a set of pre-lc regions of interest. 
//
//  In the next processing stage, the LC Candidate generation algorithm
//  will be applied at the Pre-LC region of interest level. 
//  From each Pre-LC region of interest, a set of LC Candidates will
//  be generated.  
//
//
//
// Algorithm:
//   MSDA Data Organization Algorithm 
//   
// Input:
//
//  1. ParameterLoader object reference, containing global parameters
//  2. ParameterLoader object reference, containing local parameters
// 
// Output: 
//  1. A single output file in MSDA Level 1 data format.
//
//
// Return code: 0 successful, 
// 
// Nelson.Ramirez@utsa.edu, 04/16/12
///////////////////////////////////////////////////////////////////////////////
int Utilities::workflowThree( cbi::ParameterLoader & paraObjGlobal, cbi::ParameterLoader & paraObjLocal ) {
	
	///////////////////////////////////////////////////////////////////////////
	//
	// Workflow 3
	//
	///////////////////////////////////////////////////////////////////////////
	
	return 0;
}





///////////////////////////////////////////////////////////////////////////////
//
// End MSDA namespace
//
///////////////////////////////////////////////////////////////////////////////
} // namespace msda


// General utilities required for msda algorithm
