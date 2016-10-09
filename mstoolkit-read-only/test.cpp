// This is test code to understand how to work with mzXML data format
//  Most source from:
//  msalign_uTOF - automatic alignment of 2 LC-MS datasets using feature extraction and GA 


#include <iostream>
#include "mzParser.h"


#define HPLUS_MASS 1.00727646688

#define MAX_DATAPOINTS 50000
#define N_CANDIDATES 300
#define MAX_BREAKPOINTS 24      // original value was 12
#define N_GENERATIONS 10000    // original value was 10000 (before even 5000)
#define FRACTION_KEPT 0.5 
#define MAX_FEATURES 100000

typedef struct {
  int size;
  double * xval;
  double * yval;
} spectStrct;

void freeSpectStrct(spectStrct spectrum)
{
  free(spectrum.xval);
  free(spectrum.yval);
  return;
}

ramp_fileoffset_t *pScanIndex;
int iLastScan;
struct ScanHeaderStruct scanHeader;
struct RunHeaderStruct runHeader;
ramp_fileoffset_t indexOffset;

spectStrct *spects;
double *wghs;
int spectNum;
long i;
RAMPREAL *pPeaks;
int n, n_MS_spectra=0;
RAMPFILE *pFI;

int isLittleEndian() {
long int testInt = 0x12345678;
            char *pMem;

pMem = (char *) testInt;
if (pMem[0] == 0x78)
return(1);
else
return(0);
}


int getMsSpect(spectStrct *msSpect, RAMPFILE *pFI, int scanNum[2])
{
  void getCmbSpect(spectStrct *cmbSpect, int spectNum, spectStrct *spects, double *wghs);

  msSpect->size = -1;

  // initialize
  scanNum[0] = scanNum[0] > 1 ? scanNum[0] : 1;
  scanNum[1] = scanNum[1] < iLastScan ? scanNum[1] : iLastScan;
  spectNum = scanNum[1] - scanNum[0] + 1;
  if(spectNum < 1){
    printf("invalid scan number: %d-%d (full scan range: 1-%d)\n",scanNum[0], scanNum[1], iLastScan);
    fflush(stdout);
    free (pScanIndex);
    return -1;    
  }

  spects = (spectStrct *) calloc(spectNum, sizeof(spectStrct));
  spectNum = 0;
  for (i = scanNum[0]; i <= scanNum[1]; ++i) 
    {
      if((scanHeader.msLevel==1)&&(scanHeader.peaksCount>0)) /* MS ? */
	{                         
	  spects[spectNum].size = scanHeader.peaksCount;
	  spects[spectNum].xval = (double *) calloc(spects[spectNum].size, sizeof(double));
	  spects[spectNum].yval = (double *) calloc(spects[spectNum].size, sizeof(double));
	  
	  pPeaks = readPeaks (pFI, pScanIndex[i]);
	  
	  spects[spectNum].size = 0;
	  n = 0;
	  while (pPeaks[n] != -1)
	    {
	      spects[spectNum].xval[spects[spectNum].size] = pPeaks[n];
	      n++;
	      spects[spectNum].yval[spects[spectNum].size] = pPeaks[n];
	      n++;
	      ++(spects[spectNum].size);
	    }
	  free (pPeaks);
	  n_MS_spectra++;
	  if(spects[spectNum].size > 0) ++spectNum; 
	  else freeSpectStrct(spects[spectNum]);
	}
      
    } 
  
  if(spectNum > 0) {
    wghs = (double *) calloc(spectNum, sizeof(double));
    for (i = 0; i < spectNum; ++i)
      wghs[i] = 1.;
    getCmbSpect(msSpect, spectNum, spects, wghs);
    free(wghs);
  }
  else 
    {
      printf("cannot find an MS spectrum\n"); fflush(stdout);
      for (i = 0; i < spectNum; ++i) freeSpectStrct(spects[i]);
      free(spects);
      return -1;
    }
  
  for (i = 0; i < spectNum; ++i) freeSpectStrct(spects[i]);
  free(spects);
  
  return 0;
}


int getCidSpect(double *mz, double *et, spectStrct *cidSpect, RAMPFILE *pFI, int scanNum[2])
{
  void getCmbSpect(spectStrct *cmbSpect,int spectNum, spectStrct *spects, double *wghs);

  cidSpect->size = -1;

  scanNum[0] = scanNum[0] > 1 ? scanNum[0] : 1;
  scanNum[1] = scanNum[1] < iLastScan ? scanNum[1] : iLastScan;
  spectNum = scanNum[1] - scanNum[0] + 1;
  if(spectNum < 1){
    printf("invalid scan number: %d-%d (full scan range: 1-%d)\n",scanNum[0], scanNum[1], iLastScan);
    fflush(stdout);
    free (pScanIndex);
    return -1;    
  }
  spects = (spectStrct *) calloc(spectNum, sizeof(spectStrct));
  *mz = 0.;
  *et = 0.;
  spectNum = 0;
  for (i = scanNum[0]; i <= scanNum[1]; ++i) 
    {
      if((scanHeader.msLevel==2)&&(scanHeader.peaksCount>0)) 
	{                         
	  *mz += scanHeader.precursorMZ;
	  *et += scanHeader.retentionTime/60;
	  
	  spects[spectNum].size = scanHeader.peaksCount;
	  spects[spectNum].xval = (double *) calloc(spects[spectNum].size, sizeof(double));
	  spects[spectNum].yval = (double *) calloc(spects[spectNum].size, sizeof(double));
	  
	  pPeaks = readPeaks (pFI, pScanIndex[i]);
	  
	  spects[spectNum].size = 0;
	  n = 0;
	  while (pPeaks[n] != -1)
	    {
	      spects[spectNum].xval[spects[spectNum].size] = pPeaks[n];
	      n++;
	      spects[spectNum].yval[spects[spectNum].size] = pPeaks[n];
	      n++;
	      ++(spects[spectNum].size);
	    }
	  free (pPeaks);
	  
	  if(spects[spectNum].size > 0) ++spectNum;
	  else freeSpectStrct(spects[spectNum]);
	} 
      
    } 
  
  if(spectNum > 0) {
    *mz /= spectNum;
    *et /= spectNum;
    wghs = (double *) calloc(spectNum, sizeof(double));
    for (i = 0; i < spectNum; ++i)
      wghs[i] = 1.;
    getCmbSpect(cidSpect, spectNum, spects, wghs);
    free(wghs);
  }
  else 
    {
      printf("cannot find an MS/MS spectrum\n"); fflush(stdout); 
      for (i = 0; i < spectNum; ++i) freeSpectStrct(spects[i]);
      free(spects);
      return -1;
    }
  
  for (i = 0; i < spectNum; ++i) freeSpectStrct(spects[i]);
  free(spects);
  
  return 0;
}


void getCmbSpect(spectStrct *cmbSpect, int spectNum, spectStrct *spects, double *wghs)
{
  void copySpectStrct(spectStrct * tgtSpect, spectStrct srcSpect);
  spectStrct tmpSpect[2];
  int indx, indx1, indx2;
  double tmpWghs[2] = {1., 1.};
  int i;
  
  if (spectNum < 1)
    return;
  
  // single spectrum
  if(spectNum == 1) {
    copySpectStrct(cmbSpect, spects[0]);
    if(wghs[0] != 1.) {
      for(i = 0; i < cmbSpect->size; ++i)
	cmbSpect->yval[i] *= wghs[0];
    }
    return;
  } 
  
  // 2 spectra
  if(spectNum == 2) {
    tmpSpect[0].size = spects[0].size + spects[1].size;
    tmpSpect[0].xval = (double *) calloc(tmpSpect[0].size, sizeof(double));
    tmpSpect[0].yval = (double *) calloc(tmpSpect[0].size, sizeof(double));

    indx1 = 0;
    indx2 = 0;
    indx = 0;
    while(indx1 < spects[0].size || indx2 < spects[1].size) {
      
      if(indx1 >= spects[0].size){
	tmpSpect[0].xval[indx] = spects[1].xval[indx2];
	tmpSpect[0].yval[indx] = spects[1].yval[indx2]*wghs[1];      
	++indx2;
	++indx;
      }
      else if (indx2 >= spects[1].size) {
	tmpSpect[0].xval[indx] = spects[0].xval[indx1];
	tmpSpect[0].yval[indx] = spects[0].yval[indx1]*wghs[0];      
	++indx1;
	++indx;
      }
      else if(spects[0].xval[indx1] == spects[1].xval[indx2]) {
	tmpSpect[0].xval[indx] = spects[0].xval[indx1];
	tmpSpect[0].yval[indx] = spects[0].yval[indx1]*wghs[0] 
	  + spects[1].yval[indx2]*wghs[1];      
	++indx1;
	++indx2;
	++indx;
      }
      else if(spects[0].xval[indx1] < spects[1].xval[indx2]) {
	tmpSpect[0].xval[indx] = spects[0].xval[indx1];
	tmpSpect[0].yval[indx] = spects[0].yval[indx1]*wghs[0];      
	++indx1;
	++indx;
      }
      else {
	tmpSpect[0].xval[indx] = spects[1].xval[indx2];
	tmpSpect[0].yval[indx] = spects[1].yval[indx2]*wghs[1];      
	++indx2;
	++indx;
      }
    } 
    tmpSpect[0].size = indx;
    
    copySpectStrct(cmbSpect, tmpSpect[0]);
    freeSpectStrct(tmpSpect[0]);

    return;
  }

  // at least three spectra
  indx1 = spectNum/2;
  indx2 = spectNum - spectNum/2;
  getCmbSpect(&(tmpSpect[0]), indx1, spects, wghs);
  getCmbSpect(&(tmpSpect[1]), indx2, spects+indx1, wghs+indx1);
  getCmbSpect(cmbSpect, 2, tmpSpect, tmpWghs);
  
  freeSpectStrct(tmpSpect[0]);
  freeSpectStrct(tmpSpect[1]);

  return;
}

void copySpectStrct(spectStrct * tgtSpect, spectStrct srcSpect)
{
  int i;

  tgtSpect->size = srcSpect.size;
  tgtSpect->xval = (double *) calloc(tgtSpect->size, sizeof(double));
  tgtSpect->yval = (double *) calloc(tgtSpect->size, sizeof(double));

  for (i = 0; i < tgtSpect->size; ++i) {
    tgtSpect->xval[i] = srcSpect.xval[i];
    tgtSpect->yval[i] = srcSpect.yval[i];
  }

  return;
}






int main(int argc, char** argv) {

FILE *inp, *outp;

  typedef struct {
    double B[MAX_BREAKPOINTS][2];
    int nB;
  } PLF_type; /* piecewise linear function */
 
  PLF_type best_C;
  spectStrct mzXML_spectrum;
  RAMPFILE *mzXML_file;
  char pepXML_filename[100], mzXML_filename_1[100], mzXML_filename_2[100], output_filename[100], scan_range[100], line[4000], *p, temp[30], new_mass;
  ramp_fileoffset_t offset, *scan_index;
  struct ScanHeaderStruct scan_header;
  struct RunHeaderStruct run_header;
  int nographics, MS_start_scan, MS_end_scan, SetFeatures, CountFeatures, assumed_charge, n_masses, n_unique_masses, n_unique_masses1, n_unique_masses2, *SIC_1_argmax, *SIC_2_argmax, *SIC_2a_argmax, range[2], seen_before, LC_sigma_set, Xmax_set, Ymax_set, non_matched;
  long i, j, k, n;
  double max_mass_measurement_error, start_scan, end_scan, intensity_sum, *SIC_1_max, *SIC_2_max, *SIC_2a_max, A[MAX_DATAPOINTS][2], scan_diff[MAX_DATAPOINTS], scan_diff_median, intensity_max, scan_intensity_max, *SIC_1_max_mz, *SIC_2_max_mz, *SIC_2a_max_mz, error[MAX_DATAPOINTS], swap, sum_squared_deviation, LC_sigma, Xmax, Xmax_temp, Ymax, background, oldbackground1, background1, background2, oldbackground2, costs, argFitness;
  int go, r, SumFeatures, SetBackground, feature_index, CountAttempts; 
  
  
 	cout<<"Checking mzXML file"<<endl;
  	mzXML_file=rampOpenFile(argv[1]);
  	
  	// We first need to find where in the mzXML file the index is located
  	// the index contains the references to the scan data for each scan
  	
    /* Read the offset of the index */
    indexOffset = getIndexOffset (mzXML_file);
	cout<<"indexOffset="<<indexOffset<<endl;
  
    /* Read the scan index into a vector, get LastScan */
    pScanIndex = readIndex(mzXML_file, indexOffset, &iLastScan);
  
    //cout<<iLastScan<<endl;

    readRunHeader(mzXML_file, pScanIndex, &runHeader, iLastScan);
	cout<<runHeader.scanCount<<endl;
	//cout<<runHeader.dStartTime<<endl;
	//cout<<runHeader.dEndTime<<endl;
    //cout<<runHeader.endMZ<<endl;
    //cout<<runHeader.highMZ<<endl;
    //cout<<runHeader.lowMZ<<endl;
    //cout<<runHeader.startMZ<<endl;
    
    // Read the data in the scan index for scans 1 through 100
    // the index will contain a -1 if the scan number does not exist
    
    // Loop through all scans in an mzXML file
    for ( int i = 1; i <= runHeader.scanCount; i++ ) {
    
    	
    	readHeader(mzXML_file, pScanIndex[i], &scanHeader);
    	
    	// check the scan header of each non-empty scan
    	if((scanHeader.msLevel==1) && (scanHeader.peaksCount>0)) {
    		cout<<"peaksCount="<<scanHeader.peaksCount<<endl;
    	//	cout<<"acquisitionNum="<<scanHeader.acquisitionNum<<endl;
    	//	cout<<"mergedScan="<<scanHeader.mergedScan<<endl;
    	//	cout<<"scanIndex="<<scanHeader.scanIndex<<endl;
    	//	cout<<"seqNum="<<scanHeader.seqNum<<endl;
    		cout<<"basePeakIntensity="<<scanHeader.basePeakIntensity<<endl;
    		cout<<"retentionTime="<<scanHeader.retentionTime<<endl;
    		cout<<"lowMZ="<<scanHeader.lowMZ<<endl;
    		cout<<"highMZ"<<scanHeader.highMZ<<endl;
    	}
    	
  
   
    	 pPeaks = readPeaks(mzXML_file, pScanIndex[i]);
    	 cout<<pPeaks[0]<<endl;  
    	 cout<<pPeaks[1]<<endl;
    	 //cout<<endl<<endl;
    	 //vector<double> mzdata;
    	 //vector<double> intensitydata;
    	 //int count = 0;
    	 //for ( int k = 0; k < 2*scanHeader.peaksCount; k++ ) {
    	  	//cout<<pPeaks[k]<<" ";
    	  	//mzdata.push_back(pPeaks[count]);
    	  	//intensitydata.push_back(pPeaks[count+1]);
    	  	//count++;
    	 //}
    	 //cout<<endl<<endl;
    	 
    	 
    	 //for ( int i = 0; i < mzdata.size(); i++ ) {
    	 //	cout<<hex<<mzdata.at(i)<<endl;
    	 //}
    	
    	
    	//cout<<"size of mzdata = "<<mzdata.size()<<endl;
    	//cout<<"size of intensity = "<<intensitydata.size()<<endl;
    	
    	 //cout<<"isLittleEndian"<<isLittleEndian()<<endl;
    	//if(getMsSpect(&mzXML_spectrum,mzXML_file,range)<0) {
    	//	cout<<"skipping spectrum"<<endl;
    	//	continue; /* skip spectrum */
    	//}
    	//free(pPeaks);
    	
    }  // end loop through all scans
    
/*	runHeader->scanCount=0;
	runHeader->dEndTime=0.0;
	runHeader->dStartTime=0.0;
	runHeader->endMZ=0.0;
	runHeader->highMZ=0.0;
  runHeader->lowMZ=0.0;
  runHeader->startMZ=0.0;*/



	//	MSReadeLite msr;
	//BasicSpectrum * bs = new BasicSpectrum();
	//MzParser mzp(bs);	
	//mzp.load(argv[1]);

	//int r = checkFileType(argv[1]);
	//cout<<r<<endl;

	//RAMPFILE *rf = rampOpenFile(argv[1]);
	//ramp_fileoffset_t rfot =  getIndexOffset(rf);
	//cout<<"fileoffset = "<<rfot<<endl;

	//
	//RunHeaderStruct rhs;
    //    readMSRun(rf, &rhs);
	//cout<<"scan count="<<rhs.scanCount<<endl;


	//
	//ScanHeaderStruct sh;
	//int start;
	//int end;
	//getScanSpanRange(&sh,&start,&end);
	//cout<<"scanstart= "<<start<<" scanend= "<<end<<endl;
	//cout<<"peaksCount= "<<sh.peaksCount<<endl;
        //readHeader(rf, rfot,&sh);

	



/*
ramp_fileoffset_t getIndexOffset(RAMPFILE *pFI);
InstrumentStruct* getInstrumentStruct(RAMPFILE *pFI);
void getScanSpanRange(const struct ScanHeaderStruct *scanHeader, int *startScanNum, int *endScanNum);
void rampCloseFile(RAMPFILE *pFI);
string rampConstructInputFileName(const string &basename);
char* rampConstructInputFileName(char *buf,int buflen,const char *basename);
char* rampConstructInputPath(char *buf, int inbuflen, const char *dir_in, const char *basename);
const char** rampListSupportedFileTypes();
RAMPFILE* rampOpenFile(const char *filename);
char* rampValidFileType(const char *buf);
void readHeader(RAMPFILE *pFI, ramp_fileoffset_t lScanIndex, struct ScanHeaderStruct *scanHeader);
ramp_fileoffset_t* readIndex(RAMPFILE *pFI, ramp_fileoffset_t indexOffset, int *iLastScan);
int readMsLevel(RAMPFILE *pFI, ramp_fileoffset_t lScanIndex);
void readMSRun(RAMPFILE *pFI, struct RunHeaderStruct *runHeader);
RAMPREAL* readPeaks(RAMPFILE *pFI, ramp_fileoffset_t lScanIndex);
int readPeaksCount(RAMPFILE *pFI, ramp_fileoffset_t lScanIndex);
void readRunHeader(RAMPFILE *pFI, ramp_fileoffset_t *pScanIndex, struct RunHeaderStruct *runHeader, int iLastScan);	
*/




	return 0;
}

