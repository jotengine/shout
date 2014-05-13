/***************************************************************************
 *   This file is part of the 'Shout LVCS Recognition toolkit'.            *
 ***************************************************************************
 *   Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010 by Marijn Huijbregts *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; version 2 of the License.               *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "standard.h"
#include "featureextraction.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "FFTReal.h"
#include "gaussian.h"
#include "vector.h"
#include "shout-misc.h"

using namespace StringFunctions;
using namespace WriteFileLittleBigEndian;

#define STRING_SIZE      (5000)


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

FeatureExtraction::FeatureExtraction(const char *audioIn, const char *audioOut, int mfccS,
                                     int useE, int useZC, int useD, double vtln, bool useSqrt10, bool offs,
                                     bool onlyPrepare, bool onlyM, bool pCTS) :
                                                                           Gaussian(mfccS + useE + useZC*2 + useD*(mfccS+useE+useZC*2))
{
  assert((useE == 0) || (useE == 1) || (useE == 3));
  assert((useD == 0) || (useD == 1) || (useD == 2));
  audioFileOut  = NULL;
  gaussianization = NULL;
  onlyMEL       = onlyM;
  doCMN         = true;
  doCVN         = false;
  performCTS    = pCTS;
  doSqrt10      = useSqrt10;
  normType      = FEATURENORM_SEGMENT;
  normData      = NULL;
  audioFile     = NULL;
  offset        = offs;
  useEnergy     = useE;
  useDeltas     = useD;
  useZeroCross  = useZC;
  featureVector = NULL;
  nrFrames      = 0;
  mfccSize      = mfccS;
  inverseHistogram = NULL;
  nonDeltaSize  = mfccS + useEnergy + useZeroCross*2;
  vectorSize    = nonDeltaSize + useDeltas*nonDeltaSize;

  MEL_BANKLENGTH = DEFAULT_MEL_BANKLENGTH;
  CEP_LIFTER     = DEFAULT_CEP_LIFTER;
  if(mfccS > SWITCH_MELANDCEPLIFTER)
  {
    MEL_BANKLENGTH = LARGE_MEL_BANKLENGTH;
    CEP_LIFTER     = LARGE_CEP_LIFTER;
  }
  
  assert(!onlyMEL || (vectorSize >= MEL_BANKLENGTH));
  
  // Just printing vectors:
  if((audioIn != NULL) && (strcmp(audioIn,"NULL") == 0))
  {
    FILE *feaFile   = fopen(audioOut,"rb");

    int check;
    freadEndianSafe(&check,1,sizeof(check),feaFile);
    if(check != vectorSize)
    {
      PRINTF2("You said that the vector size is %d components long\n",vectorSize);
      USER_ERROR2_NUM("The total vector size of this file doesn't match your input. Vector size of file = ",check);
    }

    float gemm[nonDeltaSize];
    for(int i=0;i<nonDeltaSize;i++)
    {
      gemm[i] = 0.0;
    }
    int len=0;
    int time = 0;
    while(!feof(feaFile))
    {
      printf("Fea%04d: ",time);
      time++;
      for(int i=0;i<vectorSize;i++)
      {
        float f;
        freadEndianSafe(&f,1,sizeof(float),feaFile);
        if(!feof(feaFile))
        {
          if(i == nonDeltaSize)
          {
            printf("\n  Delta: ");
          }
          if(i == nonDeltaSize*2)
          {
            printf("\n    D-D: ");
          }
          if(i<nonDeltaSize)
          {
            gemm[i] += f;
          }
          printf("%10g ",f);
        }
      }
      if(!feof(feaFile))
      {
        len++;
        printf("\n");
      }
    }
    fclose(feaFile);
    printf("LENGTH %d\n",len);
    for(int i=0;i<mfccSize+1;i++)
    {
      printf("Mean %d: %g\n",i,gemm[i]/((float)len));
    }
    exit(0);
  }
  
  // Calculating some helper variable that will remain constant over the entire process:
  fea_delta_noemer = 0.0;
  for(int i=1;i<=FEA_DELTA;i++)
  {
    fea_delta_noemer += (i*i);
  }
  for(int i=0;i<FEA_WINDOWSIZE;i++)
  {
    hammingTemplate[i] = (0.54-0.46*cos(TWO_PI*i / (FEA_WINDOWSIZE-1)));
  }
  fftProc = new FFTReal(FEA_WINDOWSIZE+FEA_ADDEDZEROS);

  // Preparing the MEL filter bank:
  setVTLN(vtln);

  // Prepare CEP_LIFTER:
  cepLifterFact = new Vector(vectorSize);
  cepLifterFact->setAllValues(1.0);
  if(CEP_LIFTER > 0)
  {
    double cepLdiv2  = ((double)CEP_LIFTER) / 2.0;
    double cepPidivL = PI / ((double)CEP_LIFTER);
    for(int i=0;i<mfccSize;i++)
    {
      cepLifterFact->setValue(i,(1.0 + cepLdiv2*sin(cepPidivL*((double)i))));
    }
  }

  if(!(audioIn == NULL && audioOut == NULL))
  {
    doFeatureExtraction(audioIn,audioOut, onlyPrepare);
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeatureExtraction::setOnlyMEL()
{
  onlyMEL = true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeatureExtraction::doFeatureExtraction(const char *audioIn, const char *audioOut, bool onlyPrepare)
{
  FILE *configFile = NULL;
  char str[STRING_SIZE];
  char str2[STRING_SIZE];
  // Parsing an entire list of audio files:
  if(strcmp(audioIn,"LIST") == 0)
  {
    assert(!onlyPrepare);
    configFile = fopen(audioOut,"rt");
    if(configFile == NULL)
    {
      USER_ERROR("There is no such config-list file, or I am not allowed to read it!");
    }
    audioIn    = NULL;
    char *res = fgets(str,STRING_SIZE,configFile);
    if((res != NULL) && !feof(configFile))
    {
      splitList(str,str2);
      audioIn  = str;
      audioOut = str2;
    }
  }

  // Let's go to work:
  while(audioIn != NULL)
  {
    if(featureVector != NULL)
    {
      for(int i=0;i<nrFrames+GMMLOOKAHEAD;i++)
      {
        if(featureVector[i] != NULL)
        {
          delete featureVector[i];
        }
      }
      delete[] featureVector;
      featureVector = NULL;
      nrFrames      = 0;
    }
    // Initialising our buffers:
    assert(nrFrames == 0);
    for(int i=0;i<FEA_REMEMBERWINSIZE;i++)
    {
      rememberBuffer[i] = 0;
    }

    if(audioFile != NULL)
    {
      fclose(audioFile);
    }

    if(strcmp(audioIn,"STDIN") == 0)
    {
      audioFile = stdin;
    }
    else
    {
      audioFile = fopen(audioIn, "rb");
    }

    if(audioFile == NULL)
    {
      USER_ERROR2("The audio file does not exist or I am not allowed to read it:",audioIn);
    }

    // Determine the length of the audio file:
    short int audioWindow[FEA_WINDOWSIZE];
/*    if(offset)
    {
      freadEndianSafe(&(audioWindow[0]),FEA_WINDOWSIZE,FEA_BYTES_PER_SAMPLE/2,audioFile);
    }*/
    if(audioFile == stdin)
    {
      nrFrames = FEA_ONLINE_BUFFER;
    }
    else
    {
      while(!feof(audioFile))
      {
        if(readAudioWindow(audioWindow))
        {
          nrFrames++;
        }
      }
      fseek(audioFile,0,SEEK_SET);
    }
/*    if(offset)
    {
      freadEndianSafe(&(audioWindow[0]),FEA_WINDOWSIZE,FEA_BYTES_PER_SAMPLE/2,audioFile);
    }*/
    for(int i=0;i<FEA_REMEMBERWINSIZE;i++)
    {
      rememberBuffer[i] = 0;
    }
    featureVector = new Vector*[nrFrames+GMMLOOKAHEAD];
    for(int i=0;i<nrFrames;i++)
    {
      featureVector[i] = NULL;
    }
    for(int i=nrFrames;i<nrFrames+GMMLOOKAHEAD;i++)
    {
      featureVector[i] = new Vector(vectorSize);
    }

    FILE *feaFile = NULL;
    if(audioOut != NULL)
    {
      feaFile = fopen(audioOut,"wb");
      fwrite(&vectorSize,1,sizeof(vectorSize),feaFile);
    }

    numberOfFeaturesProcessed = 0;
    if(audioFile == stdin)
    {
      createFeaturesUntilFrame(FEA_ONLINE_BUFFER-1);

      trainFinish(false,1.0e-14);                   // Determine mean...
      varianceVector->sqrtElements(true);  // Interested in the real variance...

      for(int i=0;i<FEA_ONLINE_BUFFER;i++)
      {
        assert(featureVector[i] != NULL);
        if(doCMN)
        {
          featureVector[i]->substractElements(true, meanVector);
        }
        if(doCVN)
        {
          featureVector[i]->divideElements(true,varianceVector);
        }
      }
      if(useDeltas != 0)
      {
        // Calculate delta's:
        for(int i=FEA_DELTA;i<FEA_ONLINE_BUFFER-FEA_DELTA;i++)
        {
          calculateDelta(i,0);
        }
        // Fill the first and final FEA_DELTA delta spots:
        for(int i=0;i<FEA_DELTA;i++)
        {
          for(int i2=nonDeltaSize;i2<nonDeltaSize*2;i2++)
          {
            featureVector[i]->setValue(i2, featureVector[FEA_DELTA]->getValue(i2));
            featureVector[FEA_ONLINE_BUFFER-1-i]->setValue(i2, featureVector[nrFrames-1-FEA_DELTA]->getValue(i2));
          }
        }
        if(useDeltas > 1)
        {
          // Calculate the delta-delta:
          for(int i=FEA_DELTA;i<FEA_ONLINE_BUFFER-FEA_DELTA;i++)
          {
            calculateDelta(i,1);
          }
          // Fill the first and final FEA_DELTA delta-delta spots:
          for(int i=0;i<FEA_DELTA;i++)
          {
            for(int i2=nonDeltaSize*2;i2<nonDeltaSize*3;i2++)
            {
              featureVector[i]->setValue(i2, featureVector[FEA_DELTA]->getValue(i2));
              featureVector[FEA_ONLINE_BUFFER-1-i]->setValue(i2, featureVector[nrFrames-1-FEA_DELTA]->getValue(i2));
            }
          }
        }
      }
    }
    else if(!onlyPrepare)
    {
      createFeaturesUntilFrame(nrFrames-1);
      finishExtraction();
      if(configFile != NULL)
      {
        char *res = fgets(str,STRING_SIZE,configFile);
        if((res != NULL) && !feof(configFile))
        {
          splitList(str,str2);
          audioIn  = str;
          audioOut = str2;
        }
      }
      // Write the vectors to disc:
      if(audioOut != NULL)
      {
        for(int i=0;i<nrFrames;i++)
        {
          featureVector[i]->storeFloatData(feaFile);
        }
        fclose(feaFile);
      }
    }
    audioIn = NULL;
    if(configFile != NULL)
    {
      fclose(configFile);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////
void FeatureExtraction::setAudioFileOut(const char *fileName)
{
  if(audioFileOut != NULL)
  {
    fclose(audioFileOut);
    audioFileOut = NULL;
  }
  if(fileName != NULL)
  {
    audioFileOut = fopen(fileName,"wb");
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

Vector *FeatureExtraction::getOnlineVector(double *silThreshold, bool isSil)
{
  assert(audioFile == stdin);
  int vectNr = FEA_ONLINE_BUFFER - 1;
  delete featureVector[vectNr - 3*FEA_DELTA];
  for(int i= vectNr-3*FEA_DELTA; i < vectNr;i++)
  {
    featureVector[i] = featureVector[i+1];
  }
  featureVector[vectNr]     = NULL;
  numberOfFeaturesProcessed = vectNr;
  int normID = -1;
  if(isSil)
  {
    normID = -2;
  }
  if(!createFeaturesUntilFrame(vectNr,normID))
  {
    return NULL; // End of stream..
  }

  if(featureVector[vectNr]->getValue(ASR_DEFAULT_MFCC_NUMBER) < *silThreshold)
  {
    *silThreshold = -1.0;
  }
  else
  {
    *silThreshold = 1.0;
  }

  trainFinish(false,1.0e-14);                   // Determine mean...
  varianceVector->sqrtElements(true);  // Interested in the real variance...
  assert(featureVector[vectNr] != NULL);
  if(doCMN)
  {
    featureVector[vectNr]->substractElements(true, meanVector);
  }
  if(doCVN)
  {
    featureVector[vectNr]->divideElements(true,varianceVector);
  }
  if(useDeltas != 0)
  {
    calculateDelta(vectNr-FEA_DELTA,0);
    if(useDeltas > 1)
    {
      calculateDelta(vectNr-2*FEA_DELTA,1);
    }
  }
  return featureVector[vectNr-2*FEA_DELTA];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool FeatureExtraction::createFeaturesUntilFrame(int lastFrame, int normID)
{
  int spectSub = 0;
  if((normType == FEATURENORM_CLUSTER) && (normData == NULL))
  {
    normType = FEATURENORM_SEGMENT;
  }
  // Let's calculate all static parameters:
  for(int i=numberOfFeaturesProcessed;i<=lastFrame;i++)
  {
    featureVector[i] = createMfccFrame(normID, spectSub);
    spectSub++;
    if(featureVector[i] == NULL)
    {
      return false;
    }
    if((normType == FEATURENORM_CLUSTER) && (normID >= 0))
    {
      assert(normData != NULL);
      assert(normData[normID].cmvnCluster != NULL);
      normData[normID].cmvnCluster->train(1.0,featureVector[i]);
    }
    else if(normID != -2) // for online feature extraction: do not train silence.
    {
      train(1.0,featureVector[i]);
    }
  }
  if(normType == FEATURENORM_SEGMENT)
  {
    trainFinish(true,1.0e-14);                   // Determine mean...
    varianceVector->sqrtElements(true);  // Interested in the real variance...
    if(useZeroCross > 0)
    {
      for(int aa=useEnergy;aa<useEnergy+useZeroCross+1;aa++)
      {
        meanVector->setValue    (mfccSize+aa,  0.0);
        varianceVector->setValue(mfccSize+aa,  1.0);
      }
    }
    // Normalize the mean of the mfcc's:
    for(int i=numberOfFeaturesProcessed;i<=lastFrame;i++)
    {
      assert(featureVector[i] != NULL);
      if(doCMN)
      {
        featureVector[i]->substractElements(true, meanVector);
      }
      if(doCVN)
      {
        featureVector[i]->divideElements(true,varianceVector);
      }
    }
  }
  numberOfFeaturesProcessed = lastFrame+1;
  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeatureExtraction::setVTLN(double vtln)
{
  // Preparing the MEL filter bank:
  double melMin = MEL_FREQ(MIN_MEL_WINDOW);
  double melMax = MEL_FREQ(MAX_MEL_WINDOW);
  if(performCTS)
  {
    melMin      = MEL_FREQ(MIN_CTS_MEL_WINDOW);
    melMax      = MEL_FREQ(MAX_CTS_MEL_WINDOW);
  }
  double melLength   = (melMax - melMin) / (MEL_BANKLENGTH+1.0);
  double curMelFreq  = 0.0;
  int counter = 0;
  double melCTSFreq = MEL_FREQ(MAX_CTS_MEL_WINDOW);
  melIndexAboveCTS = 0;
  for(int i=0;i<FEA_FREQWINDOWSIZE;i++)
  {
    curMelFreq = MEL_FREQ((FREQ(i)*vtln));
    if((melIndexAboveCTS == 0) && (curMelFreq > melCTSFreq))
    {
      melIndexAboveCTS = counter+1;
    }
    if(curMelFreq >= (melLength*counter+melMin))
    {
      melTemplate[counter] = i;
      if(counter > MEL_BANKLENGTH)
      {
        i = FEA_FREQWINDOWSIZE;
      }
      counter++;
    }
  }
  assert(counter > MEL_BANKLENGTH);
  if(counter == (MEL_BANKLENGTH+1))
  {
    melTemplate[counter] = FEA_FREQWINDOWSIZE-1;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeatureExtraction::finishExtraction()
{
  // Should we normalize for Mean and variance?
  if(normType == FEATURENORM_ALL)
  {
    trainFinish(true,1.0e-14);                   // Determine mean...
    varianceVector->sqrtElements(true);  // Interested in the real variance...
    if(useZeroCross > 0)
    {
      for(int aa=useEnergy;aa<useEnergy+useZeroCross+1;aa++)
      {
        meanVector->setValue    (mfccSize+aa,  0.0);
        varianceVector->setValue(mfccSize+aa,  1.0);
      }
    }
    // Normalize the mean of the mfcc's:
    for(int i=0;i<numberOfFeaturesProcessed;i++)
    {
      assert(featureVector[i] != NULL);
      if(doCMN)
      {
        featureVector[i]->substractElements(true, meanVector);
      }
      if(doCVN)
      {
        featureVector[i]->divideElements(true,varianceVector);
      }
    }
  }

  // If we want to use delta's:
  if(useDeltas != 0)
  {
    // Calculate delta's:
    for(int i=FEA_DELTA;i<nrFrames-FEA_DELTA;i++)
    {
      calculateDelta(i,0);
    }
    // Fill the first and final FEA_DELTA delta spots:
    for(int i=0;i<FEA_DELTA;i++)
    {
      for(int i2=nonDeltaSize;i2<nonDeltaSize*2;i2++)
      {
        featureVector[i]->setValue(i2, featureVector[FEA_DELTA]->getValue(i2));
        featureVector[nrFrames-1-i]->setValue(i2, featureVector[nrFrames-1-FEA_DELTA]->getValue(i2));
      }
    }
    if(useDeltas > 1)
    {
      // Calculate the delta-delta:
      for(int i=FEA_DELTA;i<nrFrames-FEA_DELTA;i++)
      {
        calculateDelta(i,1);
      }
  
      // Fill the first and final FEA_DELTA delta-delta spots:
      for(int i=0;i<FEA_DELTA;i++)
      {
        for(int i2=nonDeltaSize*2;i2<nonDeltaSize*3;i2++)
        {
          featureVector[i]->setValue(i2, featureVector[FEA_DELTA]->getValue(i2));
          featureVector[nrFrames-1-i]->setValue(i2, featureVector[nrFrames-1-FEA_DELTA]->getValue(i2));
        }
      }
    }
  }
  if(audioFile != NULL)
  {
    fclose(audioFile);
    audioFile = NULL;
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

FeatureExtraction::~FeatureExtraction()
{
  if(audioFile != NULL)
  {
    fclose(audioFile);
  }
  if(inverseHistogram != NULL)
  {
    for(int i=0;i<FEATURE_NORM_INVERSE_HISTOGRAM_BINS;i++)
    {
      delete inverseHistogram[i];
    }
    delete[] inverseHistogram;
  }
  if(normData != NULL)
  {
    for(int i=0;i<nrNormClusters;i++)
    {
      if(normData[i].gaussianizationSource != NULL)
      {
        delete[] normData[i].gaussianizationSource;
        delete[] normData[i].gaussianizationDest;
      }
      if(normData[i].cmvnCluster != NULL)
      {
        delete normData[i].cmvnCluster;
      }
      if(normData[i].histogram != NULL)
      {
        for(int i2=0;i2<FEATURE_NORM_HISTOGRAM_BINS;i2++)
        {
          if(normData[i].histogram[i2] != NULL)
          {
            delete normData[i].histogram[i2];
          }
        }
        delete[] normData[i].histogram;
      }
    }
    delete[] normData;
  }
  if(featureVector != NULL)
  {
    for(int i=0;i<nrFrames+GMMLOOKAHEAD;i++)
    {
      if(featureVector[i] != NULL)
      {
        delete featureVector[i];
      }
    }
    delete[] featureVector;
    featureVector = NULL;
    nrFrames      = 0;
  }
  delete fftProc;
  delete cepLifterFact;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool FeatureExtraction::readAudioWindow(short int *audioWindow)
{
  memcpy(audioWindow,rememberBuffer,FEA_REMEMBERWINSIZE*FEA_BYTES_PER_SAMPLE);
  int length = freadEndianSafe(&(audioWindow[FEA_REMEMBERWINSIZE]),FEA_ADDWINSIZE,FEA_BYTES_PER_SAMPLE,audioFile);
  if(length < FEA_ADDWINSIZE*FEA_BYTES_PER_SAMPLE)
  {
    return false;
  }
  if(audioFileOut != NULL)
  {
    fwrite(&(audioWindow[FEA_REMEMBERWINSIZE]),FEA_ADDWINSIZE,FEA_BYTES_PER_SAMPLE,audioFileOut);
  }
  memcpy(rememberBuffer,&(audioWindow[FEA_ADDWINSIZE]),FEA_REMEMBERWINSIZE*FEA_BYTES_PER_SAMPLE);
  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///\todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeatureExtraction::storeFeatureVectors(FILE *file)
{
  for(int i=0;i<nrFrames+GMMLOOKAHEAD;i++)
  {
    assert(featureVector[i] != NULL);
    if(onlyMEL)
    {
      featureVector[i]->storeFloatData(file,MEL_BANKLENGTH);
    }
    else
    {
      featureVector[i]->storeFloatData(file);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeatureExtraction::performPCA(Vector **pca, int len)
{
  for(int i=0;i<nrFrames;i++)
  {
    if(featureVector[i] != NULL)
    {
      Vector *pcaVector = new Vector(len);
      for(int a=0;a<len;a++)
      {
        pcaVector->setValue(a, featureVector[i]->multiplyVector(pca[a]));
      }
      delete featureVector[i];
      featureVector[i] = pcaVector;
    }
  }
  vectorSize = len;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

Vector *FeatureExtraction::createMfccFrame(int normID, int spectralSubtraction)
{
  short int audioWindow    [FEA_WINDOWSIZE];
  float     hammingWindow  [FEA_WINDOWSIZE+FEA_ADDEDZEROS];
  float     fftWindow      [FEA_WINDOWSIZE+FEA_ADDEDZEROS];
  float     melVector      [MEL_BANKLENGTH];

  Vector *mfccVector = new Vector(vectorSize);

  if(!readAudioWindow(audioWindow))
  {
    return NULL;
  }

  // abuse hammingWindow for a while (safe memory):
  // Normalize for the DC offset
  int mean = 0;
  for(int i=0;i<FEA_WINDOWSIZE;i++)
  {
    mean += audioWindow[i];
  }

  float flmean = ((float)mean)/((float)FEA_WINDOWSIZE);
  for(int i=0;i<FEA_WINDOWSIZE;i++)
  {
    hammingWindow[i] = ((float)audioWindow[i]) - flmean;
  }

  if(useZeroCross > 0)
  {
    int zeroCount = 0;
    for(int i=1;i<FEA_WINDOWSIZE;i++)
    {
      if(hammingWindow[i-1] < 0.0 && hammingWindow[i] >= 0)
      {
        zeroCount++;
      }
    }
    mfccVector->setValue(mfccSize+useEnergy,((double)zeroCount));
  }

  // Make added zero...
  for(int i=FEA_WINDOWSIZE;i<FEA_WINDOWSIZE+FEA_ADDEDZEROS;i++)
  {
    hammingWindow[i] = 0.0;
  }

  // Pre-emphesis
  for(int i=FEA_WINDOWSIZE-1;i>0;i--)
  {
    hammingWindow[i] = hammingWindow[i] - 0.97*hammingWindow[i-1];
  }

  // Hamming window and energy calculation:
  double energy     = 0.0;
  for(int i=0;i<FEA_WINDOWSIZE;i++)
  {
    hammingWindow[i] = hammingTemplate[i]*hammingWindow[i];
    energy += (hammingWindow[i]*hammingWindow[i]);
  }
  if(useEnergy > 0)
  {
    mfccVector->setValue(mfccSize,log10(1.0 + energy)); // Energy...
  }

  // Convert to frequency domain:
  fftProc->do_fft(fftWindow,hammingWindow);
  fftProc->rescale(fftWindow);

  // finish frequency calc and add VTLN:
  double maxFreq  = 0.0;
  for(int i=0;i<FEA_FREQWINDOWSIZE;i++)
  {
    int i2 = FEA_FREQWINDOWSIZE+i;
#ifdef DO_SPECTRAL_SUBTRACTION
    if((spectralSubtraction >= 0) && (vectorSize == ASR_DEFAULT_VECTORSIZE))
    {
      fftWindow[i] = (fftWindow[i]*fftWindow[i]+fftWindow[i2]*fftWindow[i2]);
    }
    else
    {
      fftWindow[i] = sqrt(fftWindow[i]*fftWindow[i]+fftWindow[i2]*fftWindow[i2]);
    }
#else
      fftWindow[i] = sqrt(fftWindow[i]*fftWindow[i]+fftWindow[i2]*fftWindow[i2]);
#endif
    if((useZeroCross > 0) && (fftWindow[i] > maxFreq))
    {
      maxFreq  = fftWindow[i];
    }
  }

#ifdef DO_SPECTRAL_SUBTRACTION
  if((spectralSubtraction >= 0) && (vectorSize == ASR_DEFAULT_VECTORSIZE))
  {
    maxFreq = sqrt(maxFreq);
  }
#endif

#ifdef DO_SPECTRAL_SUBTRACTION
  if((spectralSubtraction >= 0) && (vectorSize == ASR_DEFAULT_VECTORSIZE))
  {
    if(spectralSubtraction == 0)
    {
      for(int i=0;i<FEA_FREQWINDOWSIZE;i++)
      {
        fftWindowSpectralSubtract[i] = 0.0;
      }
    }
    if(spectralSubtraction < 5)
    {
      for(int i=0;i<FEA_FREQWINDOWSIZE;i++)
      {
        fftWindowSpectralSubtract[i] += 0.001*(fftWindow[i] / 5.0);
      }
    }
    else
    {
      for(int i=0;i<FEA_FREQWINDOWSIZE;i++)
      {
        float band = (fftWindow[i] - fftWindowSpectralSubtract[i]) / fftWindow[i];
        if(band < 0.001)
        {
          band = 0.001;
        }
        fftWindow[i] = sqrt(band * fftWindow[i]);
      }
    }
  }
#endif

  if(useZeroCross > 0)
  {
    int freqCount = 0;
    maxFreq      /= 2.0;
    for(int i=0;i<FEA_FREQWINDOWSIZE;i++)
    {
      if(fftWindow[i] > maxFreq)
      {
        freqCount++;
      }
    }
    mfccVector->setValue(mfccSize+useEnergy+1,((double)freqCount));
  }

  // MEL bank filtering:
  double energyLow  = 0.0;
  double energyHigh = 0.0;
  for(int a=0;a<MEL_BANKLENGTH;a++)
  {
    int max = melTemplate[a+1]-melTemplate[a];
    melVector[a] = 0.0;
    for(int i=0;i<max;i++)
    {
      melVector[a] += fftWindow[melTemplate[a]+i]*((float)i)/((float)max);
    }
    max = melTemplate[a+2]-melTemplate[a+1];
    for(int i=0;i<=max;i++)
    {
      melVector[a] += fftWindow[melTemplate[a+1]+i]*(((float)max)-((float)i))/((float)max);
    }
    if(!onlyMEL)
    {
      if(!doSqrt10) // default
      {
        melVector[a] = log10(1.0 + melVector[a]);
      }
      else
      {
        melVector[a] = pow(melVector[a],0.1);    // Will this help? See "Quantile based histogram equalization for noise robust large vocabulary speech recognition
      }
      if(a < melIndexAboveCTS)
      {
        energyLow += melVector[a];
      }
      else
      {
        energyHigh += melVector[a];
      }
      if(useEnergy == 3)
      {
        mfccVector->setValue(mfccSize+1,energyLow);
        mfccVector->setValue(mfccSize+2,energyHigh);
      }
    }
  }
  if(onlyMEL)
  {
    assert(MEL_BANKLENGTH <= vectorSize);
    for(int j=0;j<MEL_BANKLENGTH;j++)
    {
      mfccVector->setValue(j,melVector[j]);
    }
  }
  else
  {
#ifdef HISTOGRAM_NORM_SPECTRUM
    // Should we perform histogram normalization?
    if(normID >= 0 && inverseHistogram != NULL && normData[normID].histogram != NULL)
    {

#ifdef PRINT_HISTNORM
      for(int TEMPID = 0;TEMPID < MEL_BANKLENGTH;TEMPID++)
      {
        char str[80];
        sprintf(str,"/home/parlevink/huijbreg/gnuplot/mapping.%02d.txt",TEMPID);
        FILE *file = fopen(str,"wt");
        for(int i2=0;i2<FEATURE_NORM_HISTOGRAM_BINS-2;i2++)
        {
          float res  = normData[normID].histogram[FEATURE_NORM_HISTOGRAM_BINS-2]->getValue(TEMPID) +
                      ((float)i2)*normData[normID].histogram[FEATURE_NORM_HISTOGRAM_BINS-1]->getValue(TEMPID);
          int histI = (int)((res - normData[normID].histogram[FEATURE_NORM_HISTOGRAM_BINS-2]->getValue(TEMPID)) /
                                  normData[normID].histogram[FEATURE_NORM_HISTOGRAM_BINS-1]->getValue(TEMPID)     );
          if(histI < 0)
          {
            histI = 0;
          }
          if(histI > FEATURE_NORM_HISTOGRAM_BINS - 3)
          {
            histI = FEATURE_NORM_HISTOGRAM_BINS - 3;
          }
          res = inverseHistogram[(int)(normData[normID].histogram[histI]->getValue(TEMPID)*(FEATURE_NORM_INVERSE_HISTOGRAM_BINS-1))]->getValue(TEMPID);
          fprintf(file,"%f %f\n",res,normData[normID].histogram[i2]->getValue(TEMPID));
        }
        fclose(file);
      }
exit(0);
#endif


      for(int j=0;j<MEL_BANKLENGTH;j++)
      {
        int histI = (int)((melVector[j] - normData[normID].histogram[FEATURE_NORM_HISTOGRAM_BINS-2]->getValue(j)) /
                                          normData[normID].histogram[FEATURE_NORM_HISTOGRAM_BINS-1]->getValue(j)     );
        if(histI < 0)
        {
          histI = 0;
        }
        if(histI > FEATURE_NORM_HISTOGRAM_BINS - 3)
        {
          histI = FEATURE_NORM_HISTOGRAM_BINS - 3;
        }
assert(false); // Should not be here!
        melVector[j] = inverseHistogram[(int)(normData[normID].histogram[histI]->getValue(j)*(FEATURE_NORM_INVERSE_HISTOGRAM_BINS-1))]->getValue(j);
      }
    }
#endif

    for(int i=1;i<mfccSize+1;i++)
    {
      float sum = 0.0;
      for(int j=0;j<MEL_BANKLENGTH;j++)
      {
        sum += melVector[j]*cos( ((float)PI*i) / ((float)MEL_BANKLENGTH) * (((float)j)+0.5) );
      }
      mfccVector->setValue(i-1,sum);
    }
    // Cepstrum liftering:
    mfccVector->multiplyElements(true,cepLifterFact);
  }
  return mfccVector;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeatureExtraction::calculateDelta(int time, int offset)
{
  for(int i2=nonDeltaSize*offset;i2<nonDeltaSize*(1+offset);i2++)
  {
    float d = 0;
    for(int i=1;i<=FEA_DELTA;i++)
    {
      d += ((float)i)*(featureVector[time+i]->getValue(i2) - featureVector[time-i]->getValue(i2));
    }
    featureVector[time]->setValue(i2+nonDeltaSize,d / fea_delta_noemer);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeatureExtraction::setNormalization(FEATURENORM_TYPE normT, int nrClusters, bool doMean, bool doVar, Vector **hist)
{
  normType  = normT;
  doCMN     = doMean;
  doCVN     = doVar;
  if(hist != NULL)
  {
    int histVect = hist[0]->len();
    if(inverseHistogram != NULL)
    {
      for(int i=0;i<FEATURE_NORM_INVERSE_HISTOGRAM_BINS;i++)
      {
        delete inverseHistogram[i];
      }
      delete[] inverseHistogram;
    }
    inverseHistogram = new Vector*[FEATURE_NORM_INVERSE_HISTOGRAM_BINS];
    for(int i=0;i<FEATURE_NORM_INVERSE_HISTOGRAM_BINS;i++)
    {
      inverseHistogram[i] = new Vector(histVect);
    }
    for(int i2=0;i2<histVect;i2++)
    {
      int index = 0;
      float min   = hist[FEATURE_NORM_HISTOGRAM_BINS-2]->getValue(i2);
      float delta = hist[FEATURE_NORM_HISTOGRAM_BINS-1]->getValue(i2);
      float value = min;
      inverseHistogram[0]->setValue(i2,value);
      for(int i=1;i<FEATURE_NORM_HISTOGRAM_BINS-2;i++)
      {
        int index2  = (int)(hist[i]->getValue(i2)*(FEATURE_NORM_INVERSE_HISTOGRAM_BINS-1));
        int steps   = index2-index;
        if(steps > 0)
        {
          float d2 = delta / (float)(steps);
          float a  = 0.0;
          while(index < index2)
          {
            inverseHistogram[index]->setValue(i2,value+d2*a);
            a++;
            index++;
          }
        }
        value+=delta;
      }
      int index2 = FEATURE_NORM_INVERSE_HISTOGRAM_BINS;
      int steps   = index2-index;
      if(steps > 0)
      {
        float d2 = delta / (float)(steps);
        float a  = 0.0;
        while(index < index2)
        {
          inverseHistogram[index]->setValue(i2,value+d2*a);
          a++;
          index++;
        }
      }
      value+=delta;
    }
  }

  if(gaussianization != NULL)
  {
    for(int i=0;i<nrFrames+GMMLOOKAHEAD;i++)
    {
      if(gaussianization[i] != NULL);
      {
        delete gaussianization[i];
      }
    }
    delete[] gaussianization;
  }

  if(normData != NULL)
  {
    for(int i=0;i<nrNormClusters;i++)
    {
      if(normData[i].gaussianizationSource != NULL)
      {
        delete[] normData[i].gaussianizationSource;
        delete[] normData[i].gaussianizationDest;
      }
      if(normData[i].cmvnCluster != NULL)
      {
        delete normData[i].cmvnCluster;
      }
      if(normData[i].histogram != NULL)
      {
        for(int i2=0;i2<FEATURE_NORM_HISTOGRAM_BINS;i2++)
        {
          if(normData[i].histogram[i2] != NULL)
          {
            delete normData[i].histogram[i2];
          }
        }
        delete[] normData[i].histogram;
      }
    }
    delete[] normData;
    normData = NULL;
  }
  nrNormClusters = nrClusters;
  if(nrNormClusters > 0 && normType == FEATURENORM_CLUSTER)
  {
    normData    = new NormData[nrNormClusters];
    for(int i=0;i<nrNormClusters;i++)
    {
      normData[i].cmvnCluster     = new Gaussian(this->dim());
      normData[i].gaussianizationSource = NULL;
      normData[i].gaussianizationDest   = NULL;
      normData[i].histogram       = NULL;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeatureExtraction::finishClusterCMVN()
{
  assert(normType == FEATURENORM_CLUSTER);
  for(int i=0;i<nrNormClusters;i++)
  {
    normData[i].cmvnCluster->trainFinish(true,1.0e-14);                   // Determine mean...
    normData[i].cmvnCluster->getVariance()->sqrtElements(true);  // Interested in the real variance...
    if(useZeroCross > 0)
    {
      for(int aa=useEnergy;aa<useEnergy+useZeroCross+1;aa++)
      {
        normData[i].cmvnCluster->getMean()->setValue    (mfccSize+aa,  0.0);
        normData[i].cmvnCluster->getVariance()->setValue(mfccSize+aa,  1.0);
      }
    }
  }
  numberOfFeaturesProcessed = 0; // Reset for performClusterCMVNUntilFrame ...
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeatureExtraction::performClusterCMVNUntilFrame(int lastFrame, int normID)
{
  // Normalize the mean of the mfcc's:
  if(normID >= 0)
  {
    assert(normData != NULL);
    assert(normData[normID].cmvnCluster != NULL);
    for(int i=numberOfFeaturesProcessed;i<=lastFrame;i++)
    {
      assert(featureVector[i] != NULL);
      if(doCMN)
      {
        featureVector[i]->substractElements(true, normData[normID].cmvnCluster->getMean());
      }
      if(doCVN)
      {
        featureVector[i]->divideElements(true,    normData[normID].cmvnCluster->getVariance());
      }
    }
  }
  numberOfFeaturesProcessed = lastFrame+1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeatureExtraction::setHistNormModel(int normID, Vector **h)
{
  assert(normData != NULL);
  if(normData[normID].histogram != NULL)
  {
    for(int i2=0;i2<FEATURE_NORM_HISTOGRAM_BINS;i2++)
    {
      if(normData[normID].histogram[i2] != NULL)
      {
        delete normData[normID].histogram[i2];
      }
    }
    delete[] normData[normID].histogram;
  }
  normData[normID].histogram = h;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeatureExtraction::initializeHistNorm()
{
  numberOfFeaturesProcessed = 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeatureExtraction::performHistNormUntilFrame(int lastFrame, int normID)
{
  // Normalize the mean of the mfcc's:
  if(normID >= 0 && inverseHistogram != NULL && normData[normID].histogram != NULL)
  {
    int histSize = normData[normID].histogram[FEATURE_NORM_HISTOGRAM_BINS-2]->len() - 1; // Forget about energy...
    for(int j=0;j<histSize;j++)
    {
      float mean = normData[normID].histogram[FEATURE_NORM_HISTOGRAM_BINS-2]->getValue(j);
      float step = normData[normID].histogram[FEATURE_NORM_HISTOGRAM_BINS-1]->getValue(j);
      for(int i=numberOfFeaturesProcessed;i<=lastFrame;i++)
      {
        int histI = (int)((featureVector[i]->getValue(j) - mean) / step);
        if(histI < 0)
        {
          histI = 0;
        }
        if(histI > FEATURE_NORM_HISTOGRAM_BINS - 3)
        {
          histI = FEATURE_NORM_HISTOGRAM_BINS - 3;
        }
        featureVector[i]->setValue(j,
            inverseHistogram[(int)(normData[normID].histogram[histI]->getValue(j)*(FEATURE_NORM_INVERSE_HISTOGRAM_BINS-1))]->getValue(j));
      }
    }
  }
  numberOfFeaturesProcessed = lastFrame+1;
}

