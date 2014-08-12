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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "standard.h"
#include "featurepool.h"
#include "vector.h"
#include "gaussian.h"
#include "phonemodel.h"
#include "trainphonemodel.h"
#include "shout-misc.h"


using namespace StringFunctions;
using namespace WriteFileLittleBigEndian;


#define NR_SUBSET_BINS  (50)

#define PCA_SIZE (20)

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor will populate the vector pool by reading an input file.
/////////////////////////////////////////////////////////////////////////////////////////////////////

FeaturePool::FeaturePool(char *featureFileName, FILE_TYPE fType, bool noBuffering, double vtln, bool doOffset, bool onlyPrepare)
{
  onlyMEL     = false;
  multiStream = NULL;
  textFeatureDistribution = NULL;
  infoBlock.numberOfChannels = 1;
  for(int i=0;i<FEATUREPOOL_MAX_NR_STREAMS;i++)
  {
    infoBlock.weight[i]        = 0.0;
    infoBlock.weightBackup[i]  = 0.0;
    infoBlock.distrG[i]        = new Gaussian(1);
  }
  infoBlock.weight[0] = 1.0;
  infoBlock.weightBackup[0] = 1.0;

  fileType      = fType;
  for(int i=0;i<2;i++)
  {
    feaEnergy[i].feaExtract = NULL;
    feaEnergy[i].buffer     = NULL;
  }
  featureFile   = NULL;
  featureNorm   = FEATURENORM_SEGMENT;

  createNewFeatureExtraction(NULL,vtln,doOffset,onlyPrepare);

  filter             = NULL;
  filterEnabled      = false;
  buffer             = NULL;
  for(int i=0;i<MAX_NR_SEGMENTATIONS;i++)
  {
    freeID[i]          = true;
    segList[i]         = NULL;
    fillPool_curSeg[i] = NULL;
  }
  createNewPool(featureFileName, NULL, NULL, noBuffering,onlyPrepare,-1,vtln);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor will populate the vector pool by reading an input file.
/////////////////////////////////////////////////////////////////////////////////////////////////////

FeaturePool::FeaturePool(int nrStreams, FeaturePool_Multistream_file *files, double *initWeight, bool onlyPrepare)
{
  onlyMEL = false;
  assert(nrStreams <= FEATUREPOOL_MAX_NR_STREAMS);
  multiStream = files;
  textFeatureDistribution    = NULL;
  for(int i=0;i<FEATUREPOOL_MAX_NR_STREAMS;i++)
  {
    infoBlock.weight[i]        = 0.0;
    infoBlock.weightBackup[i]  = 0.0;
    infoBlock.distrG[i]        = new Gaussian(1);
  }
  infoBlock.weight[0]        = 1.0;
  infoBlock.weightBackup[0]  = 1.0;

  for(int i=0;i<nrStreams;i++)
  {
    PRINTF2("Loading info for channel %d...\n",i);fflush(stdout);
    files[i].featurePool = new FeaturePool(files[i].fileName,files[i].fileType,false,1.0,false,onlyPrepare);
    int p = files[i].featurePool->getPoolLength();
    if(i == 0 || p<poolLength)
    {
      poolLength = p;
    }
  }

  featureFile                = NULL;
  for(int i=0;i<2;i++)
  {
    feaEnergy[i].feaExtract  = NULL;
    feaEnergy[i].buffer      = NULL;
  }
  filter                     = NULL;
  filterEnabled              = false;
  featureNorm                = FEATURENORM_CLUSTER;
  fileOffset                 = 0;
  fileType                   = FILE_TYPE_SHOUT; // The buffer is ours..
  buffer                     = new Vector*[poolLength];
  infoBlock.numberOfChannels = nrStreams;
  infoBlock.lastComponent[0] = files[0].featurePool->getVectorSize() - 1;
  for(int i=1;i<nrStreams;i++)
  {
    infoBlock.lastComponent[i] = infoBlock.lastComponent[i-1] + files[i].featurePool->getVectorSize();
  }
  setStreamWeights(initWeight, true);
  vectorSize                 = infoBlock.lastComponent[nrStreams-1] + 1;
  poolFilled                 = poolLength;
  restrictStart              = 0;
  restrictEnd                = poolLength - 1;
  for(int i=0;i<MAX_NR_SEGMENTATIONS;i++)
  {
    freeID[i]                = true;
    segList[i]               = NULL;
    fillPool_curSeg[i]       = NULL;
  }
  if(!onlyPrepare)
  {
    createNewMultiStreamPool(files, -1);
  }
  PRINTF2("Created a multi-stream feature pool: %d. Setting weights at: ", infoBlock.numberOfChannels);
  for(int i=0;i<infoBlock.numberOfChannels;i++)
  {
    PRINTF3("%f (%d) ",initWeight[i], files[i].featurePool->getVectorSize());
  }
  PRINTF("\n");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::createNewMultiStreamPool(FeaturePool_Multistream_file *files, int segID)
{
  int nrChannels = infoBlock.numberOfChannels;
  for(int i=0;i<nrChannels;i++)
  {
    files[i].featurePool->createNewPool(files[i].fileName, NULL, NULL,false,false,segID);
  }
  /// \todo generalize this code to N channels!
  assert(nrChannels < 4);
  for(int i=0;i<poolLength;i++)
  {
    if(nrChannels == 2)
    {
      buffer[i] = new Vector(files[0].featurePool->getVector(NULL,i),files[1].featurePool->getVector(NULL,i),NULL);
    }
    else
    {
      buffer[i] = new Vector(files[0].featurePool->getVector(NULL,i),files[1].featurePool->getVector(NULL,i),files[2].featurePool->getVector(NULL,i));
    }
  }
  for(int i=0;i<nrChannels;i++)
  {
    delete files[i].featurePool;
    files[i].featurePool = NULL;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor initializes the pool, but it does not fill it with data. This can be done
/// later with the addFeatures() method. This way of populating the pool is typically done when
/// the size of the pool is not known beforhand (for example when collecting adaptation data).
/////////////////////////////////////////////////////////////////////////////////////////////////////

FeaturePool::FeaturePool(int poolS, int vectorS)
{
  onlyMEL     = false;
  multiStream = NULL;
  textFeatureDistribution = NULL;
  for(int i=0;i<MAX_NR_SEGMENTATIONS;i++)
  {
    freeID[i]          = true;
    segList[i]         = NULL;
    fillPool_curSeg[i] = NULL;
  }
  fileOffset         = 0;
  poolFilled         = 0;
  filter             = NULL;
  filterEnabled      = false;
  featureFile        = NULL;
  poolLength         = poolS;
  vectorSize         = vectorS;
  buffer             = new Vector*[poolLength];
  feaExtract         = NULL;
  for(int i=0;i<2;i++)
  {
    feaEnergy[i].feaExtract = NULL;
    feaEnergy[i].buffer     = NULL;
  }
  for(int i=0;i<poolLength;i++)
  {
    buffer[i] = NULL;
  }
  restrictStart = 0;
  restrictEnd   = poolLength - 1;

  for(int i=0;i<FEATUREPOOL_MAX_NR_STREAMS;i++)
  {
    infoBlock.weight[i]        = 0.0;
    infoBlock.weightBackup[i]  = 0.0;
    infoBlock.distrG[i]        = new Gaussian(1);
  }
  infoBlock.weight[0] = 1.0;
  infoBlock.weightBackup[0] = 1.0;
  infoBlock.numberOfChannels = 1;
  infoBlock.lastComponent[0] = vectorSize-1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor cleans up administration
/////////////////////////////////////////////////////////////////////////////////////////////////////

FeaturePool::~FeaturePool()
{
  if(textFeatureDistribution != NULL)
  {
    delete textFeatureDistribution;
    textFeatureDistribution = NULL;
  }
  for(int i=0;i<FEATUREPOOL_MAX_NR_STREAMS;i++)
  {
    delete infoBlock.distrG[i];
  }
  if(buffer != NULL && ((fileType == FILE_TYPE_SHOUT) || (fileType == FILE_TYPE_TEXT)))  // Else the buffer isn't ours...
  {
    for(int i=0;i<poolLength;i++)
    {
      if(buffer[i] != NULL)
      {
        delete buffer[i];
      }
    }
    delete[] buffer;
  }
  if(feaExtract != NULL)
  {
    delete feaExtract;
    feaExtract = NULL;
  }
  if(featureFile != NULL)
  {
    fclose(featureFile);
  }
  for(int i=0;i<MAX_NR_SEGMENTATIONS;i++)
  {
    while(segList[i] != NULL)
    {
      SegmentationList *t = segList[i]->next;
      if(segList[i]->infoString != NULL)
      {
        delete[] segList[i]->infoString;
      }
      delete segList[i];
      segList[i] = t;
    }
  }
  if(filter != NULL)
  {
    delete[] filter;
  }
  for(int i=0;i<2;i++)
  {
    if(feaEnergy[i].feaExtract != NULL)
    {
      delete feaEnergy[i].feaExtract;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///\todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::setStreamWeights(double *w, bool forEver)
{
  for(int i=0;i<infoBlock.numberOfChannels;i++)
  {
    infoBlock.weight[i] = w[i];
    if(forEver)
    {
      infoBlock.weightBackup[i] = w[i];
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///\todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::resetStreamWeights()
{
  for(int i=0;i<infoBlock.numberOfChannels;i++)
  {
    infoBlock.weight[i] = infoBlock.weightBackup[i];
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///\todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::createNewFeatureExtraction(char *featureFileName, double vtln, bool doOffset, bool onlyPrepare)
{
  switch(fileType)
  {
    case FILE_TYPE_RAW_DIARIZATION:
      featureNorm= DIARIZATION_DEFAULT_FEATURENORM;
      feaExtract = new FeatureExtraction(featureFileName,NULL,DIARIZATION_DEFAULT_MFCC_NUMBER,
                                                          DIARIZATION_DEFAULT_MFCC_ENERGY,
                                                          DIARIZATION_DEFAULT_MFCC_ZEROCROSS,
                                                          DIARIZATION_DEFAULT_MFCC_DELTA,
                                                          vtln,
                                                          DIARIZATION_DEFAULT_SQRT10,
                                                          doOffset,onlyPrepare);
      break;
    case FILE_TYPE_RAW_HISTOGRAM_NORM:
      featureNorm= VTLN_DEFAULT_FEATURENORM;
#ifdef HISTOGRAM_NORM_SPECTRUM
      feaExtract = new FeatureExtraction(featureFileName,NULL,DEFAULT_MEL_BANKLENGTH,0,0,0,1.0,ASR_DEFAULT_SQRT10,doOffset,onlyPrepare,true);
#else
      feaExtract = new FeatureExtraction(featureFileName,NULL,VTLN_DEFAULT_MFCC_NUMBER,
                                                          VTLN_DEFAULT_MFCC_ENERGY,
                                                          VTLN_DEFAULT_MFCC_ZEROCROSS,
                                                          VTLN_DEFAULT_MFCC_DELTA,
                                                          vtln,
                                                          VTLN_DEFAULT_SQRT10,
                                                          doOffset,onlyPrepare);
#endif
      break;
    case FILE_TYPE_RAW_VTLN:
      featureNorm= VTLN_DEFAULT_FEATURENORM;
      feaExtract = new FeatureExtraction(featureFileName,NULL,VTLN_DEFAULT_MFCC_NUMBER,
                                                          VTLN_DEFAULT_MFCC_ENERGY,
                                                          VTLN_DEFAULT_MFCC_ZEROCROSS,
                                                          VTLN_DEFAULT_MFCC_DELTA,
                                                          vtln,
                                                          VTLN_DEFAULT_SQRT10,
                                                          doOffset,onlyPrepare);
      break;
    case FILE_TYPE_RAW_SAD:
      featureNorm= SAD_DEFAULT_FEATURENORM;
      feaExtract = new FeatureExtraction(featureFileName,NULL,SAD_DEFAULT_MFCC_NUMBER,
                                                          SAD_DEFAULT_MFCC_ENERGY,
                                                          SAD_DEFAULT_MFCC_ZEROCROSS,
                                                          SAD_DEFAULT_MFCC_DELTA,
                                                          vtln,
                                                          SAD_DEFAULT_SQRT10,
                                                          doOffset,onlyPrepare);
      break;
    case FILE_TYPE_RAW_DECODING:
      featureNorm= ASR_DEFAULT_FEATURENORM;
      feaExtract = new FeatureExtraction(featureFileName,NULL,ASR_DEFAULT_MFCC_NUMBER,
                                                          ASR_DEFAULT_MFCC_ENERGY,
                                                          ASR_DEFAULT_MFCC_ZEROCROSS,
                                                          ASR_DEFAULT_MFCC_DELTA,
                                                          vtln,
                                                          ASR_DEFAULT_SQRT10,
                                                          doOffset,onlyPrepare);
      break;
    case FILE_TYPE_RAW_SPKREC:
      featureNorm= SPKREC_DEFAULT_FEATURENORM;
      feaExtract = new FeatureExtraction(featureFileName,NULL,SPKREC_DEFAULT_MFCC_NUMBER,
                                                          SPKREC_DEFAULT_MFCC_ENERGY,
                                                          SPKREC_DEFAULT_MFCC_ZEROCROSS,
                                                          SPKREC_DEFAULT_MFCC_DELTA,
                                                          vtln,
                                                          SPKREC_DEFAULT_SQRT10,
                                                          doOffset,onlyPrepare,false,true);
      break;
    case FILE_TYPE_TEXT:
    case FILE_TYPE_SHOUT:
    default:
      feaExtract = NULL;
    break;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///\todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::setFeatureNormalization(FEATURENORM_TYPE m)
{
  featureNorm = m;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///\todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::setOnlyMEL()
{
  onlyMEL = true;
  feaExtract->setOnlyMEL();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///\todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::storeFeatureVectors(FILE *file)
{
  feaExtract->storeFeatureVectors(file);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///\todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::createNewPool(char *featureFileName, char *histNormName, char *histNormSpkName, bool noBuffering, bool onlyPrepare, int vtlnSegID, double globalVtln)
{
  SegmentationAdmin readPool_curSeg;
  fileOffset         = 0;
  poolFilled         = 0;

  if(buffer != NULL && ((fileType == FILE_TYPE_SHOUT) || (fileType == FILE_TYPE_TEXT)))  // it is our buffer..
  {
    if(!onlyPrepare && vtlnSegID >= 0)
    {
      // There was a 'onlyPrepare' call that filled the buffer. Let's first use it to determine the distribution of the data for normalization!
      if(textFeatureDistribution != NULL)
      {
        delete textFeatureDistribution;
      }
      textFeatureDistribution = new Gaussian(vectorSize);
      readPool_curSeg.curSeg    = segList[vtlnSegID];
      readPool_curSeg.timeStamp = 0;
      while(readPool_curSeg.curSeg != NULL)
      {
        int timeFrom = readPool_curSeg.curSeg->firstFrame;
        int timeTo = readPool_curSeg.curSeg->lastFrame;
        for(int i=timeFrom;i<=timeTo;i++)
        {
          textFeatureDistribution->train(1.0,buffer[i]);
        }
        readPool_curSeg.curSeg = readPool_curSeg.curSeg->next;
      }
      textFeatureDistribution->trainFinish(true,1.0e-14);
      textFeatureDistribution->getVariance()->sqrtElements(true);
      PRINTF("The variance for delay-normalization is: ");
      textFeatureDistribution->getVariance()->printVector();
    }
    // Now delete the previous buffer.
    for(int i=0;i<poolLength;i++)
    {
      if(buffer[i] != NULL)
      {
        delete buffer[i];
      }
    }
    delete[] buffer;
  }

  vectorSize         = 0;
  poolLength         = 0;
  if(fileType == FILE_TYPE_TEXT)
  {
    featureFile        = fopen(featureFileName,"rt");
  }
  else
  {
    featureFile        = fopen(featureFileName,"rb");
  }
  if(featureFile == NULL)
  {
    USER_ERROR2("Couldn't open the feature file:",featureFileName);
  }


  if(globalVtln >= 0.0 && (fileType != FILE_TYPE_TEXT) && (fileType != FILE_TYPE_SHOUT))
  {
    feaExtract->setVTLN(globalVtln);
  }

  bool doMeanNorm = true;
  bool doVarNorm  = true;
  switch(fileType)
  {
    case FILE_TYPE_RAW_DIARIZATION:
      doMeanNorm = DIARIZATION_DEFAULT_CMN;
      doVarNorm  = DIARIZATION_DEFAULT_CVN;
      break;
    case FILE_TYPE_RAW_HISTOGRAM_NORM:
      doMeanNorm = false;
      doVarNorm  = false;
      break;
    case FILE_TYPE_RAW_VTLN:
      doMeanNorm = VTLN_DEFAULT_CMN;
      doVarNorm  = VTLN_DEFAULT_CVN;
      break;
    case FILE_TYPE_RAW_SAD:
      doMeanNorm = SAD_DEFAULT_CMN;
      doVarNorm  = SAD_DEFAULT_CVN;
      break;
    case FILE_TYPE_RAW_DECODING:
      doMeanNorm = ASR_DEFAULT_CMN;
      doVarNorm  = ASR_DEFAULT_CVN;
      break;
    case FILE_TYPE_RAW_SPKREC:
      doMeanNorm = SPKREC_DEFAULT_CMN;
      doVarNorm  = SPKREC_DEFAULT_CVN;
      break;
    case FILE_TYPE_TEXT:
    case FILE_TYPE_SHOUT:
    default:
      doMeanNorm = false;
      doVarNorm  = false;
    break;
  }

  Vector **histogram = NULL;
  PhoneModel *trainHistModel = NULL;
  int histVectorSize = 0;
  // First prepare for histogram normalization:
  if(!onlyPrepare && histNormName != NULL && strlen(histNormName) > 0)
  {
    FILE *histFile = fopen(histNormName,"rb");
    if(histFile == NULL)
    {
      USER_ERROR("Could not open the histogram feature normalization file.");
    }
    int histNumberOfModels;
    int histNumberOfClusters;
    freadEndianSafe(&histVectorSize,1,sizeof(histVectorSize),histFile);
    freadEndianSafe(&histNumberOfModels,1,sizeof(histNumberOfModels),histFile);
    freadEndianSafe(&histNumberOfClusters,1,sizeof(histNumberOfClusters),histFile);
    if(histNumberOfModels != 1 || histNumberOfClusters != 1)
    {
      USER_ERROR("The histogram feature normalization file is corrupt.");
    }
    histogram = new Vector*[FEATURE_NORM_HISTOGRAM_BINS];
    for(int i=0;i<FEATURE_NORM_HISTOGRAM_BINS;i++)
    {
      histogram[i] = new Vector(histVectorSize);
    }
    trainHistModel = new PhoneModel(histFile,histVectorSize);
    fclose(histFile);
    trainHistModel->getSilSumHist(histogram);


#ifdef PRINT_HISTNORM
    for(int TEMPID = 0;TEMPID < histVectorSize;TEMPID++)
    {
      char str[80];
      sprintf(str,"/home/parlevink/huijbreg/gnuplot/global.%02d.txt",TEMPID);
      FILE *file = fopen(str,"wt");
      float x  = histogram[FEATURE_NORM_HISTOGRAM_BINS-2]->getValue(TEMPID);
      float dx = histogram[FEATURE_NORM_HISTOGRAM_BINS-1]->getValue(TEMPID);
      for(int i2=0;i2<FEATURE_NORM_HISTOGRAM_BINS-2;i2++)
      {
        float res = histogram[i2]->getValue(TEMPID);
        fprintf(file,"%f %f\n",x,res);
        x+=dx;
      }
      fclose(file);
    }
#endif


  }


  switch(fileType)
  {
    case FILE_TYPE_RAW_HISTOGRAM_NORM:
    case FILE_TYPE_RAW_DIARIZATION:
    case FILE_TYPE_RAW_SPKREC:
    case FILE_TYPE_RAW_DECODING:
    case FILE_TYPE_RAW_VTLN:
    case FILE_TYPE_RAW_SAD:
      feaExtract->setNormalization(featureNorm,0,doMeanNorm,doVarNorm,histogram);
    // Reading header information (extracting vector size and file length):
      if(vtlnSegID >= 0)
      {
        feaExtract->setNormalization(featureNorm,numberOfSpkIDs[vtlnSegID],doMeanNorm,doVarNorm,histogram);
        int timeSil = 0;
        int time    = 0;
        feaExtract->doFeatureExtraction(featureFileName, NULL, true);
        poolLength = feaExtract->getPoolSize();
        vectorSize = feaExtract->getVectorSize();
        readPool_curSeg.curSeg    = segList[vtlnSegID];
        readPool_curSeg.timeStamp = 0;
        timeSil = 0;
        time    = 0;

        while(readPool_curSeg.curSeg != NULL)
        {
          timeSil = readPool_curSeg.curSeg->firstFrame-1;
          if(timeSil > time)
          {
            feaExtract->createFeaturesUntilFrame(timeSil,-1);
          }
          time = readPool_curSeg.curSeg->lastFrame;
          if(globalVtln < 0.0)
          {
            double vtln = VTLN_START + VTLN_STEPSIZE*((double)getSegmentID(TRAIN_CONTEXTRIGHT,&readPool_curSeg));
            feaExtract->setVTLN(vtln);
          }
          int normID = 0;
          while(normID < numberOfSpkIDs[vtlnSegID] && spkID[vtlnSegID][normID] != getSegmentID(0,&readPool_curSeg))
          {
            normID++;
          }
          assert(normID < numberOfSpkIDs[vtlnSegID]);
          feaExtract->createFeaturesUntilFrame(time,normID);
          readPool_curSeg.curSeg = readPool_curSeg.curSeg->next;
        }
        if(time < poolLength)
        {
          feaExtract->createFeaturesUntilFrame(poolLength-1,-1);
        }

        if(!onlyMEL)
        {
          if(featureNorm == FEATURENORM_CLUSTER && (doMeanNorm || doVarNorm))
          {
            // Now repeat for Mean and Variance normalization:
            feaExtract->finishClusterCMVN();
            readPool_curSeg.curSeg = segList[vtlnSegID];
            timeSil = 0;
            time    = 0;
            while(readPool_curSeg.curSeg != NULL)
            {
              timeSil = readPool_curSeg.curSeg->firstFrame-1;
              if(timeSil > time)
              {
                feaExtract->performClusterCMVNUntilFrame(timeSil,0);
              }

              time = readPool_curSeg.curSeg->lastFrame;
              int normID = 0;
              while(normID < numberOfSpkIDs[vtlnSegID] && spkID[vtlnSegID][normID] != getSegmentID(0,&readPool_curSeg))
              {
                normID++;
              }
              feaExtract->performClusterCMVNUntilFrame(time,normID);
              readPool_curSeg.curSeg = readPool_curSeg.curSeg->next;
            }
            if(time < poolLength)
            {
              feaExtract->performClusterCMVNUntilFrame(poolLength-1,0);
            }
          }
          // Done, let's finish up:
          feaExtract->finishExtraction();
        }

        if(false)
        {
/*
          readPool_curSeg.curSeg = segList[vtlnSegID];
          timeSil = 0;
          time    = 0;
          feaExtract->initGaussianization();
          while(readPool_curSeg.curSeg != NULL)
          {
            timeSil = readPool_curSeg.curSeg->firstFrame-1;
            if(timeSil > time)
            {
              feaExtract->performGaussianizationUntilFrame(timeSil,0);
            }
            time = readPool_curSeg.curSeg->lastFrame;
            int normID = 0;
            while(normID < numberOfSpkIDs[vtlnSegID] && spkID[vtlnSegID][normID] != getSegmentID(0,&readPool_curSeg))
            {
              normID++;
            }
            feaExtract->performGaussianizationUntilFrame(time,normID);
            readPool_curSeg.curSeg = readPool_curSeg.curSeg->next;
          }
          if(time < poolLength)
          {
            feaExtract->performGaussianizationUntilFrame(poolLength-1,0);
          }
          feaExtract->finishGaussianization();
*/

          FILE *pca = fopen("/RAID/t0/people/huijbreg/projects/evalita/ubm/pca.male.dat","rb");
          assert(pca != NULL);
          Vector *pcaVector[vectorSize];
          for(int i=0;i<vectorSize;i++)
          {
            pcaVector[i] = new Vector(pca,vectorSize);
          }

          feaExtract->performPCA(pcaVector,PCA_SIZE);

          for(int i=0;i<vectorSize;i++)
          {
            delete pcaVector[i];
          }

          vectorSize = PCA_SIZE;
          fclose(pca);
        }

      }
      else
      {
        feaExtract->doFeatureExtraction(featureFileName, NULL, onlyPrepare);
        poolLength = feaExtract->getPoolSize();
        vectorSize = feaExtract->getVectorSize();
      }
    break;

    case FILE_TYPE_TEXT:
      poolLength = readHeader_TextFormat();
    break;
    case FILE_TYPE_SHOUT:  // this is the default:
    default:
      poolLength = readHeader_ShoutFormat();
    break;
  }

  // Filling buffer (if wanted):
  if(noBuffering)
  {
    poolLength = 1;
    buffer     = new Vector*[1];
    buffer[0]  = NULL;
  }
  else
  {
    // Be aware: the following method MUST fill all entries OR make unfilled entries NULL.
    switch(fileType)
    {
      case FILE_TYPE_RAW_HISTOGRAM_NORM:
      case FILE_TYPE_RAW_DIARIZATION:
      case FILE_TYPE_RAW_SPKREC:
      case FILE_TYPE_RAW_DECODING:
      case FILE_TYPE_RAW_VTLN:
      case FILE_TYPE_RAW_SAD:
        buffer = feaExtract->getPool();
      break;
      case FILE_TYPE_TEXT:
        buffer = new Vector*[poolLength];
        fillBuffer_TextFormat(poolLength);
      break;
      case FILE_TYPE_SHOUT:  // this is the default:
      default:
        buffer = new Vector*[poolLength];
        fillBuffer_ShoutFormat(poolLength);
      break;
    }
    poolFilled = poolLength;
    if(featureFile != NULL)
    {
      fclose(featureFile);
      featureFile = NULL;
    }
  }


  restrictStart = 0;
  restrictEnd   = poolLength - 1;

  infoBlock.numberOfChannels = 1;
  infoBlock.lastComponent[0] = vectorSize-1;
  infoBlock.weight[0]        = 1.0;
  for(int i=1;i<FEATUREPOOL_MAX_NR_STREAMS;i++)
  {
    infoBlock.weight[i] = 0.0;
  }

  if(textFeatureDistribution != NULL)
  {
    delete textFeatureDistribution;
    textFeatureDistribution = NULL;
  }
  if(trainHistModel != NULL)
  {
    delete trainHistModel;
    trainHistModel = NULL;
  }
  if(histogram != NULL)
  {
    for(int i=0;i<FEATURE_NORM_HISTOGRAM_BINS;i++)
    {
      if(histogram[i] != NULL)
      {
        delete histogram[i];
      }
    }
    delete[] histogram;
    histogram = NULL;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::retrieveEnergy(char *featureFileName, int channelNr)
{
  assert(channelNr > 0 && channelNr < 3);
  channelNr--;
  if(feaEnergy[channelNr].feaExtract != NULL)
  {
    delete feaEnergy[channelNr].feaExtract;
    feaEnergy[channelNr].feaExtract = NULL;
  }
  feaEnergy[channelNr].feaExtract = new FeatureExtraction(NULL,NULL,0,ZEROCROSS_COMP,1,0,1.0);
  feaEnergy[channelNr].feaExtract->doFeatureExtraction(featureFileName, NULL);
  feaEnergy[channelNr].buffer = feaEnergy[channelNr].feaExtract->getPool();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::maxSegmentLength(int segNr, int labelNr, int size, int minRemainder)
{
  assert(segNr >= 0 && segNr < MAX_NR_SEGMENTATIONS);
  assert(labelNr >= 0);

  SegmentationList *seg = segList[segNr];
  while(seg != NULL)
  {
    if(seg->ID[0] == labelNr)
    {
      if(seg->lastFrame - seg->firstFrame + 1 > size+minRemainder)
      {
        SegmentationList *segNew = new SegmentationList;
        segNew->infoString = NULL;
        segNew->next = seg->next;
        segNew->firstFrame = seg->firstFrame + size;
        segNew->lastFrame  = seg->lastFrame;
        for(int i=0;i<4;i++)
        {
          segNew->ID[i] = seg->ID[i];
        }
        seg->lastFrame = segNew->firstFrame - 1;
        seg->next      = segNew;
      }
    }
    seg = seg->next;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::splitOnEnergy(int segNr, int sourceLabel, int silLabel, int soundLabel, int nrSamplesSil, int nrSamplesSound, int comp)
{
  assert(feaEnergy[0].feaExtract != NULL);
  PRINTF("Now determine the NON-SPEECH split!\n");
  assert(segNr >= 0 && segNr < MAX_NR_SEGMENTATIONS);
  assert(sourceLabel >= 0);

  double min     =  9.0e300;
  double max     = -9.0e300;
  double binSize = 0.0;
  int    bins[NR_SUBSET_BINS*2];

  for(int i=0;i<NR_SUBSET_BINS*2;i++)
  {
    bins[i]   = 0;
  }

  int    tot   = 0;
  SegmentationList *seg = segList[segNr];
  while(seg != NULL)
  {
    if(seg->ID[0] == sourceLabel)
    {
      tot++;
    }
    seg = seg->next;
  }

  double meanSet[tot];
  int    sizeSet[tot];

  int curSet = 0;
  seg = segList[segNr];
  while(seg != NULL)
  {
    double score   = 0.0;
    int    nrS     = 0;
    if(seg->ID[0] == sourceLabel)
    {
      for(int i=seg->firstFrame;i<=seg->lastFrame;i++)
      {
        score += feaEnergy[0].buffer[i]->getValue(comp);
        nrS++;
      }
      sizeSet[curSet]   = nrS;
      meanSet[curSet]   = score/((double)nrS);
      if(meanSet[curSet] > max)
      {
        max = meanSet[curSet];
      }
      if(meanSet[curSet] < min)
      {
        min = meanSet[curSet];
      }
      curSet++;
    }
    seg = seg->next;
  }
  binSize = (max-min) / ((double)NR_SUBSET_BINS*2);

  for(int i=0;i<tot;i++)
  {
    int index = ((int)((meanSet[i] - min) / binSize));
    if(index < 0)
    {
      index = 0;
    }
    if(index >= NR_SUBSET_BINS*2)
    {
      index = NR_SUBSET_BINS*2-1;
    }
    bins[index]+=sizeSet[i];
  }

  int index = NR_SUBSET_BINS*2-1;
  int totS  = 0;
  while((totS < nrSamplesSound) && index > NR_SUBSET_BINS)
  {
    totS += bins[index];
    index--;
  }
  double minProbSound = min + index*binSize;
  index = 0;
  totS  = 0;
  while((totS < nrSamplesSil) && index < NR_SUBSET_BINS)
  {
    totS += bins[index];
    index++;
  }
  double maxProbSil = min + index*binSize;


  curSet = 0;
  seg = segList[segNr];
  while(seg != NULL)
  {
    if(seg->ID[0] == sourceLabel)
    {
      if(meanSet[curSet] >= minProbSound)
      {
        seg->ID[0] = soundLabel;
      }
      else if(meanSet[curSet] <= maxProbSil)
      {
        seg->ID[0] = silLabel;
      }
      else
      {
        seg->ID[0] = 525;
      }
      curSet++;
    }
    seg = seg->next;
  }











/*
  assert(feaEnergy != NULL);
  PRINTF("Now determine the NON-SPEECH split!\n");
  assert(segNr >= 0 && segNr < MAX_NR_SEGMENTATIONS);
  assert(sourceLabel >= 0);
  assert(destLabel >= 0);


  SegmentationList *seg = segList[segNr];
  int nrSeg   = 0;
  while(seg != NULL)
  {
    if(seg->ID[0] == sourceLabel)
    {
      nrSeg++;
    }
    seg = seg->next;
  }
  seg = segList[segNr];

  Vector *data[nrSeg];
  int curSeg = 0;
  while(seg != NULL)
  {
    Vector mean(bufferEnergy[0]->len());
    mean.setAllValues(0.0);
    int    tot  = 0;
    if(seg->ID[0] == sourceLabel)
    {
      for(int i=seg->firstFrame;i<=seg->lastFrame;i++)
      {
        mean.addElements(true,bufferEnergy[i]);
        tot++;
      }
      mean.multiplyElements(true,1.0/((double)tot));
      data[curSeg] = mean.copyVector();
      curSeg++;
    }
    seg = seg->next;
  }

  if(threshold < -10.0)
  {

    Gaussian *g1 = new Gaussian(bufferEnergy[0]->len());

    for(int i=0;i<nrSeg;i++)
    {
      g1->train(1.0,data[i]);
    }
    g1->trainFinish();


    printf("Variance, mean: %f %f\n",g1->getVariance()->getValue(0),g1->getMean()->getValue(0));

    Gaussian *g2 = new Gaussian(g1);
    Gaussian *g4 = new Gaussian(g1);
    Gaussian *g3 = new Gaussian(g1,false);
    g3->shiftMean();
    g1->shiftMean();

    for(int it=0;it<10;it++)
    {
      for(int i=0;i<nrSeg;i++)
      {
        double s1  = g1->getP(data[i]);
        double s2  = g2->getP(data[i]);
        double s3  = g3->getP(data[i]);
        double s4  = g4->getP(data[i]);
        if(s3>s4)
        {
          g3->train(s3,data[i]);
        }
        else
        {
          g4->train(s4,data[i]);
        }
        g1->train(s1,data[i]);
        g2->train(s2,data[i]);
      }
      printf("%f  %f  %f\n",g1->getTrainingDenominator(),g2->getTrainingDenominator(),g3->getTrainingDenominator());
      g1->trainFinish(true,0.0001);
      g2->trainFinish(true,0.0001);
      g3->trainFinish(true,0.0001);
      g4->trainFinish(true,0.0001);
    }
    printf("Variance, mean: %f %f\n",g1->getVariance()->getValue(0),g1->getMean()->getValue(0));
    printf("Variance, mean: %f %f\n",g2->getVariance()->getValue(0),g2->getMean()->getValue(0));
    printf("Variance, mean: %f %f\n",g3->getVariance()->getValue(0),g3->getMean()->getValue(0));
    printf("Variance, mean: %f %f\n",g4->getVariance()->getValue(0),g4->getMean()->getValue(0));


    printf("SHOULD WE??: %f  %f  %f\n",g1->getKLDistance(g2), g2->getKLDistance(g1) ,g1->getKLDistance(g2) - g2->getKLDistance(g1));

    if(true) //g1->getKLDistance(g2) - g2->getKLDistance(g1) > 0.5)
    {
      seg = segList[segNr];
      curSeg = 0;
      while(seg != NULL)
      {
        if(seg->ID[0] == sourceLabel)
        {
          if(g3->getP(data[curSeg]) < g4->getP(data[curSeg]))
          {
            seg->ID[0] = destLabel;
          }
          curSeg++;
        }
        seg = seg->next;
      }
    }
    delete g1;
    delete g2;
    delete g3;
    delete g4;
  }
  else
  {
    seg = segList[segNr];
    curSeg = 0;
    while(seg != NULL)
    {
      if(seg->ID[0] == sourceLabel)
      {
        printf("%f\n",data[curSeg]->getValue(0));
        if(data[curSeg]->getValue(0) > threshold)
        {
          seg->ID[0] = destLabel;
          printf("   BINGO!\n");
        }
        curSeg++;
      }
      seg = seg->next;
    }
  }
  for(int i=0;i<nrSeg;i++)
  {
    delete data[i];
  }
  */
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::removeLowEnergy(int segID, int channelNr, bool doShort)
{
  int minSpeech    = CTS_MINSPEECH;
  int minSilence   = CTS_MINSILENCE;
  int silColor     = CTS_SILCOLOR;
  double silFactor = CTS_SILFACTOR;

  if(doShort)
  {
    minSpeech    = CTS_SHORT_MINSPEECH;
    minSilence   = CTS_SHORT_MINSILENCE;
    silColor     = CTS_SHORT_SILCOLOR;
    silFactor    = CTS_SHORT_SILFACTOR;
  }

  int tempID         = claimSegmentationID();
  double threshold   = getAverageEnergy(segID,channelNr,silFactor);

  PRINTF3("The minimum energy threshold for channel %d is %f\n",channelNr,threshold);
  bool isLast        = false;
  bool inFilter      = true;
  SegmentationAdmin readPool_curSeg;
  getFirstVectorFirstSegment(&readPool_curSeg,segID,-1,true, &isLast, &inFilter);
  int time           = getCurSegmentStart(&readPool_curSeg);
  Vector *energyV    = feaEnergy[channelNr-1].buffer[time];
  int startTime      = time;
  int silenceCounter = 0;
  while(energyV != NULL)
  {
    if(energyV->getValue(0) > threshold)
    {
      // It is not possible that the counter is higher than CTS_MINSILENCE...
      silenceCounter = 0;
    }
    else
    {
      silenceCounter = 1;
    }
    while(!isLast)
    {
      getNextVector(&readPool_curSeg,&isLast, &inFilter);
      time++;
      energyV = feaEnergy[channelNr-1].buffer[time];
      if(energyV->getValue(0) > threshold)
      {
        if(silenceCounter > minSilence)
        {
          // Remove this silence region
          int endTime = time - silenceCounter + silColor;
          if(endTime - startTime > minSpeech)
          {
            addSegment(tempID,startTime,endTime,channelNr);
          }
          startTime = time - silColor;
        }
        silenceCounter = 0;
      }
      else
      {
        silenceCounter++;
      }
      if(isLast)
      {
        if(silenceCounter > minSilence)
        {
          int endTime = time - silenceCounter + silColor;
          if(endTime - startTime > minSpeech)
          {
            addSegment(tempID,startTime,endTime,channelNr);
          }
          startTime = time - silColor;
        }
        else
        {
          if(time - startTime > minSpeech)
          {
            addSegment(tempID,startTime,time,channelNr);
          }
          startTime = time;
        }
      }
    }
    energyV = getFirstVectorNextSegment(&readPool_curSeg,true, &isLast, &inFilter);
    if(energyV != NULL)
    {
      time         = getCurSegmentStart(&readPool_curSeg);
      energyV      = feaEnergy[channelNr-1].buffer[time];
    }
    startTime      = time;
    silenceCounter = 0;
  }




  resetSegmentation(segID);
  segList[segID]          = segList[tempID];
  segList[tempID]         = NULL;
  fillPool_curSeg[segID]  = NULL;
  fillPool_curSeg[tempID] = NULL;
  freeID[tempID]          = true;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

double FeaturePool::getAverageEnergy(int segNr, int channelNr, double varianceFactor)
{
  assert(channelNr > 0 && channelNr < 3);
  channelNr--;
  assert(feaEnergy[channelNr].feaExtract != NULL);
  assert(segNr >= 0 && segNr < MAX_NR_SEGMENTATIONS);

  // Assuming that energy has been added with addEnergyToFeatures if not there by default:
  int dimension = 0;


  Gaussian meanG(feaEnergy[channelNr].buffer[0]->len());
  SegmentationList *seg = segList[segNr];

  while(seg != NULL)
  {
    for(int i=seg->firstFrame;i<=seg->lastFrame;i++)
    {
      meanG.train(1.0,feaEnergy[channelNr].buffer[i]);
    }
    seg = seg->next;
  }
  meanG.trainFinish(true,1.0e-14);

  double threshold = meanG.getMean()->getValue(dimension);

  // Now start over with all samples > threshold...
  seg = segList[segNr];

  while(seg != NULL)
  {
    for(int i=seg->firstFrame;i<=seg->lastFrame;i++)
    {
      if(feaEnergy[channelNr].buffer[i]->getValue(dimension) > threshold)
      {
        meanG.train(1.0,feaEnergy[channelNr].buffer[i]);
      }
    }
    seg = seg->next;
  }
  meanG.trainFinish(true,1.0e-14);

  PRINTF3("Average energy: %f. Variation: %f\n",meanG.getMean()->getValue(dimension),meanG.getVariance()->getValue(dimension));
  return meanG.getMean()->getValue(dimension) - varianceFactor * sqrt(meanG.getVariance()->getValue(dimension));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Add a filter to the pool. All frames that are filtered out are never offered to the user.
/// with enableFilter() it is possible to enable or disable the filter.
/// Only shout-format feature files are processed. If no file is available (NULL), the pool itself is used.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::setFilter(char *featureFileName, int dimension, double threshold)
{
  FILE *fFile = fopen(featureFileName,"rb");
  if(filter != NULL)
  {
    delete[] filter;
  }
  filter = new bool[poolLength];
  double *value = new double[poolLength];
  for(int i=0;i<poolLength;i++)
  {
    filter[i] = false;
    value[i]  = 0.0;
  }
  int vectSize = vectorSize;
  if(fFile != NULL)
  {
    freadEndianSafe(&vectSize,1,sizeof(vectorSize),fFile);
    assert(vectSize > dimension);
  }
  int i = 0;
  Gaussian meanG(vectSize);



  Vector *v = buffer[0];
  while((i < poolLength) && (v != NULL) && (fFile == NULL || !feof(fFile)))
  {
    if(fFile != NULL)
    {
      v = new Vector(fFile,vectSize,true);
    }
    value[i] = v->getValue(dimension);
    meanG.train(1.0,v);
    if(fFile != NULL)
    {
      delete v;
    }
    i++;
    if(i<poolLength)
    {
      v = buffer[i];
    }
  }
  meanG.trainFinish(true,1.0e-14);
  double minValue = meanG.getMean()->getValue(dimension) - meanG.getVariance()->getValue(dimension)*threshold;
  for(int a=0;a<i;a++)
  {
    filter[a] = (value[a] > minValue);
  }
  delete[] value;
  filterEnabled = false; //true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Add a new feature file to the pool. The fileOffset will be updated so that new read RTTM
/// files will point to the correct position in time...
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::addFeatureFile(char *featureFileName, FILE_TYPE fType, int ignoreFrames, int trainFrames, double vtln)
{

  fileType = fType;

  fileOffset = poolFilled;
  // Reading header information (extracting vector size and file length):
  int length = 0;
  switch(fileType)
  {
    case FILE_TYPE_RAW_HISTOGRAM_NORM:
    case FILE_TYPE_RAW_DIARIZATION:
    case FILE_TYPE_RAW_SPKREC:
    case FILE_TYPE_RAW_DECODING:
    case FILE_TYPE_RAW_VTLN:
    case FILE_TYPE_RAW_SAD:
      if(feaExtract != NULL)
      {
        delete feaExtract;
        feaExtract = NULL;
      }
      createNewFeatureExtraction(featureFileName,vtln,false,false);
      length = feaExtract->getPoolSize();
      assert(vectorSize == feaExtract->getVectorSize());
    break;
    case FILE_TYPE_TEXT:
      featureFile = fopen(featureFileName,"rt");
      if(featureFile == NULL)
      {
        USER_ERROR2("Couldn't open the feature file:",featureFileName);
      }
      length = readHeader_TextFormat();
    case FILE_TYPE_SHOUT:  // this is the default:
    default:
      featureFile = fopen(featureFileName,"rb");
      if(featureFile == NULL)
      {
        USER_ERROR2("Couldn't open the feature file:",featureFileName);
      }
      length = readHeader_ShoutFormat();
    break;
  }
  if(poolLength < (length + poolFilled))
  {
    USER_ERROR("The total pool size is too small!");
  }
  if(length < (ignoreFrames + trainFrames))
  {
    trainFrames  = length-ignoreFrames;
    if(trainFrames <= 5)
    {
      trainFrames  = length;
      ignoreFrames = 0;
    }
  }


  // Be aware: the following method MUST fill all entries OR make unfilled entries NULL.
  switch(fileType)
  {
    case FILE_TYPE_RAW_HISTOGRAM_NORM:
    case FILE_TYPE_RAW_DIARIZATION:
    case FILE_TYPE_RAW_SPKREC:
    case FILE_TYPE_RAW_DECODING:
    case FILE_TYPE_RAW_VTLN:
    case FILE_TYPE_RAW_SAD:
      fillBuffer_fromFeatureExtraction(trainFrames,ignoreFrames);
    break;
    case FILE_TYPE_TEXT:
      fillBuffer_TextFormat(trainFrames,ignoreFrames);
      fclose(featureFile);
    break;
    case FILE_TYPE_SHOUT:  // this is the default:
    default:
      fillBuffer_ShoutFormat(trainFrames,ignoreFrames);
      fclose(featureFile);
    break;
  }
  poolFilled += trainFrames;
  featureFile = NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Enables or disables the filter
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::enableFilter(bool enable)
{
  filterEnabled = false; //enable;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// All filtered frames are SPK01, the remaining frames SPK02..
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::filterToSegmentation(int segNr, FeaturePool *feaPool)
{
  if(feaPool == NULL)
  {
    feaPool = this;
  }
  assert(feaPool->filter != NULL);
  resetSegmentation(segNr);
  bool active = feaPool->filter[0];
  int oldTime = 0;
  for(int i=0;i<feaPool->poolLength;i++)
  {
    if(active != feaPool->filter[i])
    {
      // Check if this will be long enough (for silence):
      bool switching = true;
      if(active)
      {
        int len = 0;
        while((len < 10) && (len+i<feaPool->poolLength) && (feaPool->filter[len+i] != active))
        {
          len++;
        }
        if((len < 10) && (len+i<feaPool->poolLength) && (active))
        {
          i+=len;
          switching = false;
        }
      }
      if(switching)
      {
        if(active)
        {
          addSegment(segNr,oldTime,i-1,2);
        }
        else
        {
          addSegment(segNr,oldTime,i-1,1);
        }
        active = feaPool->filter[i];
        oldTime = i;
      }
    }
  }
  if(active)
  {
    addSegment(segNr,oldTime,feaPool->poolLength-1,2);
  }
  else
  {
    addSegment(segNr,oldTime,feaPool->poolLength-1,1);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::segmentationSubset(int sourceSegNr, int destSegNr, int useLabel, PhoneModel *model, int nrSamples)
{
  assert(useLabel >= 0);
  resetSegmentation(destSegNr);
  bool isLast   = false;
  bool inFilter = true;
  Vector *trainVector = NULL;
  double mean    = 0.0;
  double max     = -9.0e300;
  double binSize = 0.0;
  int    bins[NR_SUBSET_BINS];
  for(int i=0;i<NR_SUBSET_BINS;i++)
  {
    bins[i] = 0;
  }

  SegmentationAdmin readPool_curSeg;
  trainVector = getFirstVectorFirstSegment(&readPool_curSeg, sourceSegNr,useLabel,useLabel < 0, &isLast, &inFilter);
  int    tot   = 0;
  while(trainVector != NULL)
  {
    tot++;
    trainVector = getFirstVectorNextSegment(&readPool_curSeg,useLabel < 0, &isLast, &inFilter);
  }

  double meanSet[tot];
  int    sizeSet[tot];
  int curSet = 0;
  trainVector = getFirstVectorFirstSegment(&readPool_curSeg,sourceSegNr,useLabel,useLabel < 0, &isLast, &inFilter);
  while(trainVector != NULL)
  {
    double score = 0.0;
    int    nrS   = 0;
    int contextKey = getSegmentID(TRAIN_CONTEXTKEY,&readPool_curSeg);
    if(inFilter)
    {
      score += model->getLogPDFProbability(contextKey,trainVector);
      nrS++;
    }
    while(!isLast)
    {
      trainVector = getNextVector(&readPool_curSeg,&isLast, &inFilter);
      if(inFilter)
      {
        score += model->getLogPDFProbability(contextKey,trainVector);
        nrS++;
      }
    }
    sizeSet[curSet] = nrS;
    meanSet[curSet] = score/((double)nrS);
    if(meanSet[curSet] > max)
    {
      max = meanSet[curSet];
    }
    mean += meanSet[curSet] / ((double)tot);
    curSet++;
    trainVector = getFirstVectorNextSegment(&readPool_curSeg,useLabel < 0, &isLast, &inFilter);
  }
  binSize = (max-mean) / ((double)NR_SUBSET_BINS);

  for(int i=0;i<tot;i++)
  {
    int index = ((int)((meanSet[i] - mean) / binSize));
    if(index < 0)
    {
      index = 0;
    }
    if(index >= NR_SUBSET_BINS)
    {
      index = NR_SUBSET_BINS-1;
    }
    bins[index]+=sizeSet[i];
  }

  int index = NR_SUBSET_BINS-1;
  int totS  = bins[index];
  while((totS < nrSamples) && index > 0)
  {
    index--;
    totS += bins[index];
  }
  double minProb = mean + index*binSize;

  trainVector = getFirstVectorFirstSegment(&readPool_curSeg,sourceSegNr,useLabel,useLabel<0, &isLast, &inFilter);
  while(trainVector != NULL)
  {
    double score = 0.0;
    int nrS   = 0;
    int contextKey = getSegmentID(TRAIN_CONTEXTKEY,&readPool_curSeg);
    if(inFilter)
    {
      score += model->getLogPDFProbability(contextKey,trainVector);
      nrS++;
    }
    while(!isLast)
    {
      trainVector = getNextVector(&readPool_curSeg,&isLast, &inFilter);
      if(inFilter)
      {
        score += model->getLogPDFProbability(contextKey,trainVector);
        nrS++;
      }
    }

    if(score/((double)nrS) > minProb)
    {
      addSegment(destSegNr, readPool_curSeg.curSeg->firstFrame, readPool_curSeg.curSeg->lastFrame, readPool_curSeg.curSeg->ID[0],
                readPool_curSeg.curSeg->ID[TRAIN_CONTEXTKEY],readPool_curSeg.curSeg->ID[TRAIN_CONTEXTLEFT],readPool_curSeg.curSeg->ID[TRAIN_CONTEXTRIGHT]);
    }
    trainVector = getFirstVectorNextSegment(&readPool_curSeg,useLabel<0, &isLast, &inFilter);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// In order to make it possible to share a feature pool between objects, a user may claim a free
/// segmentation ID. This ID will be freed using freeSegmentationID().
/////////////////////////////////////////////////////////////////////////////////////////////////////

int FeaturePool::claimSegmentationID()
{
  int id = 0;
  while(id < MAX_NR_SEGMENTATIONS && !freeID[id])
  {
    id++;
  }
  assert(id < MAX_NR_SEGMENTATIONS);
  resetSegmentation(id);
  freeID[id] = false;
  return id;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will free an ID claimed by claimSegmentationID().
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::freeSegmentationID(int segNr)
{
  assert(segNr >= 0 && segNr < MAX_NR_SEGMENTATIONS);
  resetSegmentation(segNr);
  freeID[segNr] = true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The pool can have multiple segmentations. In order to reset one of those segmentations,
/// this method needs to be called.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::resetSegmentation(int segNr)
{
  assert(segNr >= 0 && segNr < MAX_NR_SEGMENTATIONS);
  while(segList[segNr] != NULL)
  {
    SegmentationList *t = segList[segNr]->next;
    if(segList[segNr]->infoString != NULL)
    {
      delete[] segList[segNr]->infoString;
    }
    delete segList[segNr];
    segList[segNr] = t;
  }
  fillPool_curSeg[segNr] = NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Adds a segment in the current segmentation linked-list.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::addSegment(int segNr, int firstFrame, int lastFrame, int ID, int cKey, int cLeft, int cRight)
{
  assert(firstFrame >= 0 && lastFrame < poolLength);
  assert(lastFrame-firstFrame >= 0);

  if(fillPool_curSeg[segNr] == NULL)
  {
    fillPool_curSeg[segNr] = new SegmentationList;
    segList[segNr]         = fillPool_curSeg[segNr];
  }
  else
  {
    fillPool_curSeg[segNr]->next = new SegmentationList;
    fillPool_curSeg[segNr]       = fillPool_curSeg[segNr]->next;
  }
  fillPool_curSeg[segNr]->infoString             = NULL;
  fillPool_curSeg[segNr]->ID[0]                  = ID;
  fillPool_curSeg[segNr]->ID[TRAIN_CONTEXTKEY]   = cKey;
  fillPool_curSeg[segNr]->ID[TRAIN_CONTEXTLEFT]  = cLeft;
  fillPool_curSeg[segNr]->ID[TRAIN_CONTEXTRIGHT] = cRight;
  fillPool_curSeg[segNr]->firstFrame             = firstFrame;
  fillPool_curSeg[segNr]->lastFrame              = lastFrame;
  fillPool_curSeg[segNr]->next                   = NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Adds a segment in the current segmentation linked-list at the very first slot.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::addSegmentToFront(int segNr, int firstFrame, int lastFrame, int ID, int cKey, int cLeft, int cRight)
{
  assert(firstFrame >= 0 && lastFrame < poolLength);
  assert(lastFrame-firstFrame >= 0);

  SegmentationList *newSeg       = new SegmentationList;
  newSeg->infoString             = NULL;
  newSeg->ID[0]                  = ID;
  newSeg->ID[TRAIN_CONTEXTKEY]   = cKey;
  newSeg->ID[TRAIN_CONTEXTLEFT]  = cLeft;
  newSeg->ID[TRAIN_CONTEXTRIGHT] = cRight;
  newSeg->firstFrame             = firstFrame;
  newSeg->lastFrame              = lastFrame;
  newSeg->next                   = segList[segNr];
  segList[segNr]                 = newSeg;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Once a segment is initialized, it can be filled from file using this method.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::setCurSegmentVectors(int segNr, FILE *inFile, bool inputIsFloat)
{
  assert(fillPool_curSeg[segNr] != NULL);
  assert(inFile != NULL);
  assert(buffer != NULL);

  for(int i=fillPool_curSeg[segNr]->firstFrame;i<=fillPool_curSeg[segNr]->lastFrame;i++)
  {
    if(buffer[i] != NULL)
    {
      delete buffer[i];
    }
    buffer[i] = new Vector(inFile,vectorSize,inputIsFloat);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Some users might want to fill the pool one Vector at a time. That's possible in combination
/// with the FeaturePool(int poolS, int vectorS) constructor. Segments can be added with
/// addSegment. Once a segment is added, this method can fill this segment with Vector objects.
///
/// The offset is the number of frames from the start of the segment.
///
/// This method will COPY the vector. The user still needs to delete its own Vector object.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::setCurSegmentVector(int segNr, int offset, Vector *vector)
{
  assert(fillPool_curSeg[segNr] != NULL);
  int time = offset + fillPool_curSeg[segNr]->firstFrame;
  assert(offset >= 0 && time < fillPool_curSeg[segNr]->lastFrame);
  assert(vector != NULL);
  assert(buffer != NULL);

  if(buffer[time] == NULL)
  {
    buffer[time] = new Vector(vectorSize);
  }
  buffer[time]->setAllValues(vector);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Reads the header information of a shout feature vector file.
/// Actually, this file type doesn't have the length of the file in its header, so the length is
/// calculated by just reading trough the file...
/////////////////////////////////////////////////////////////////////////////////////////////////////

int FeaturePool::readHeader_ShoutFormat()
{
  int length = 0;
  assert(featureFile != NULL);
  freadEndianSafe(&vectorSize,1,sizeof(vectorSize),featureFile);
  Vector dummy(vectorSize);
  int readBytes = dummy.readFloatData(featureFile);
  while(readBytes > 0)
  {
    length++;
    readBytes = dummy.readFloatData(featureFile);
  }
  fseek(featureFile,sizeof(vectorSize),SEEK_SET);  // sizeof(vectorSize) to get past the header..
  return length;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

int FeaturePool::readHeader_TextFormat()
{
  int length = 0;
  char str[MAXSTR];
  char str2[MAXSTR];
  assert(featureFile != NULL);
  strcpy(str,"");
  char *resChar = fgets(str,MAXSTR,featureFile);
  assert(resChar != NULL);
  splitList(str,str2);

  ///todo BUGFIX THIS!!!!
  if(strlen(str) <= 0)
  {
    strcpy(str,str2);
    splitList(str,str2);
  }

  vectorSize = 1;
  while(strlen(str2)>0)
  {
    vectorSize++;
    strcpy(str,str2);
    splitList(str,str2);
  }
  while(strlen(str)>0)
  {
    length++;
    strcpy(str,"");
    char *resChar = fgets(str,MAXSTR,featureFile);
  }
  fseek(featureFile,0,SEEK_SET);
  return length;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Creates a vector buffer and fills it with the vectors of the text-formatted file.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::fillBuffer_TextFormat(int length, int ignoreFrames)
{
  char str[MAXSTR];
  char str2[MAXSTR];
  assert(featureFile != NULL);
  for(int i=0;i<ignoreFrames;i++)
  {
    int length = 0;
    strcpy(str,"");
    char *resChar = fgets(str,MAXSTR,featureFile);
    assert(resChar != NULL);
  }
  for(int i=0;i<length;i++)
  {
    buffer[i+poolFilled] = new Vector(vectorSize);
    strcpy(str,"");
    char *resChar = fgets(str,MAXSTR,featureFile);
    assert(resChar != NULL);
    splitList(str,str2);

    ///todo BUGFIX THIS!!!!
    if(strlen(str) <= 0)
    {
      strcpy(str,str2);
      splitList(str,str2);
    }

    int dim = 0;
    while(strlen(str2)>0)
    {
      double val = atof(str);
      buffer[i+poolFilled]->setValue(dim,val);
      dim++;
      strcpy(str,str2);
      splitList(str,str2);
    }
    double val = atof(str);


    buffer[i+poolFilled]->setValue(dim,val);
    if(textFeatureDistribution != NULL)
    {
      buffer[i+poolFilled]->substractElements(true, textFeatureDistribution->getMean());
      buffer[i+poolFilled]->divideElements(true, textFeatureDistribution->getVariance());
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Creates a vector buffer and fills it with the vectors of the shout-formatted file.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::fillBuffer_ShoutFormat(int length, int ignoreFrames)
{
  for(int i=0;i<ignoreFrames;i++)
  {
    Vector *v = new Vector(featureFile,vectorSize,true);
    delete v;
  }
  for(int i=0;i<length;i++)
  {
    buffer[i+poolFilled] = new Vector(featureFile,vectorSize,true);
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::fillBuffer_fromFeatureExtraction(int length, int ignoreFrames)
{
  Vector **feaBuffer = feaExtract->getPool();
  for(int i=0;i<length;i++)
  {
    buffer[i+poolFilled] = feaBuffer[i+ignoreFrames]->copyVector();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the total length of the pool (the number of vectors)
/////////////////////////////////////////////////////////////////////////////////////////////////////

int FeaturePool::getPoolLength()
{
  return poolFilled;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the length of all labels of segmentation segNr (the number of vectors)
/////////////////////////////////////////////////////////////////////////////////////////////////////

int FeaturePool::getClusterLength(int segNr,int label)
{
  assert(poolLength > 1);
  assert(segNr >= 0 && segNr < MAX_NR_SEGMENTATIONS);
  assert(restrictStart == 0);
  assert(restrictEnd >= poolLength-1);

  int res = 0;
  SegmentationList *seg = segList[segNr];
  while(seg != NULL)
  {
    if(seg->ID[0] == label)
    {
      int start = seg->firstFrame;
      if(start < restrictStart)
      {
        start = restrictStart;
      }
      int end  = seg->lastFrame;
      if(end > restrictEnd)
      {
        end = restrictEnd;
      }
      int length = (end - start + 1);
      if(length < 0)
      {
        length = 0;
      }
      res += length;
    }
    seg = seg->next;
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns a pointer to the vector stored at timestamp 'time'.It will assert when time is bigger
/// than the poolLength. If 'noBuffering' was true the poolLength is one. getVector(0) will then
/// return the pointer to the last queried vector.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Vector **FeaturePool::getVectorList(SegmentationAdmin *admin, int time)
{
  assert(buffer != NULL);
  assert(time >= 0 && time < poolLength);
  if(admin != NULL)
  {
    admin->timeStamp = time;
  }
  return &buffer[time];
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns a pointer to the vector that follows the last queried vector. If no one follows,
/// it returns NULL.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Vector **FeaturePool::getNextVectorList(SegmentationAdmin *admin, bool *isLastFromSegment, bool *inFilter)
{

  if(isLastFromSegment != NULL)
  {
    *isLastFromSegment = true;
  }
  // restrictEnd == poolLength-1;
  if((featureFile == NULL) && (admin->timeStamp >= restrictEnd))
  {
    return NULL;
  }
  assert((featureFile != NULL) || (admin->timeStamp+1 < poolLength));
  assert(buffer != NULL);
  if(inFilter != NULL)
  {
    *inFilter = true;
  }
  // If the file is still open, there will be no buffer, so just return the next vector from file:
  if(featureFile != NULL)
  {
    if(buffer[0] != NULL)
    {
      delete buffer[0];
    }
    buffer[0] = new Vector(featureFile,vectorSize,true);
    return buffer;
  }
  else if(admin->timeStamp+1 <= restrictEnd)
  {
    admin->timeStamp++;
    if(isLastFromSegment != NULL && admin->curSeg != NULL)
    {
      *isLastFromSegment = ((admin->timeStamp == admin->curSeg->lastFrame) || (admin->timeStamp == restrictEnd));
    }
    if(inFilter != NULL && filterEnabled)
    {
      *inFilter = filter[admin->timeStamp];
    }
    return &buffer[admin->timeStamp];
  }
  return NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Searches for the first segment with ID[0] = ID and if it finds it, it will return the first
/// Vector of that segment.
/// If no Vector is found, NULL is returned.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Vector **FeaturePool::getFirstVectorFirstSegmentList(SegmentationAdmin *admin,int segmentationNr, int ID, bool dontCareID, bool *isLastFromSegment, bool *inFilter)
{
  assert(poolLength > 1);
  assert(segmentationNr >= 0 && segmentationNr < MAX_NR_SEGMENTATIONS);
  SegmentationList *readPool_curSeg = admin->curSeg;

  if(isLastFromSegment != NULL)
  {
    *isLastFromSegment = true;
  }
  if(inFilter != NULL)
  {
    *inFilter = true;
  }

  readPool_curSeg = segList[segmentationNr];
  admin->prevSeg  = NULL;
  if(restrictStart > 0)
  {
    while(readPool_curSeg != NULL && readPool_curSeg->lastFrame < restrictStart)
    {
      admin->prevSeg  = readPool_curSeg;
      readPool_curSeg = readPool_curSeg->next;
    }
  }
  if(!dontCareID)
  {
    while(readPool_curSeg != NULL && readPool_curSeg->ID[0] != ID)
    {
      admin->prevSeg  = readPool_curSeg;
      readPool_curSeg = readPool_curSeg->next;
    }
  }
  if(readPool_curSeg != NULL)
  {
    if(isLastFromSegment != NULL)
    {
      *isLastFromSegment = ((readPool_curSeg->firstFrame == readPool_curSeg->lastFrame) ||
                            (readPool_curSeg->firstFrame == restrictEnd));
    }
    if(inFilter != NULL && filterEnabled)
    {
      *inFilter = filter[readPool_curSeg->firstFrame];
    }
    if(restrictStart > readPool_curSeg->firstFrame)
    {
      admin->curSeg = readPool_curSeg;
      return getVectorList(admin,restrictStart);
    }
    admin->curSeg = readPool_curSeg;
    return getVectorList(admin, readPool_curSeg->firstFrame);
  }
  admin->curSeg = readPool_curSeg;
  return NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Checks if there are any segments with ID in segmentation segNr
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool FeaturePool::doSegmentsExist(int segmentationNr, int ID)
{
  assert(poolLength > 1);
  assert(segmentationNr >= 0 && segmentationNr < MAX_NR_SEGMENTATIONS);

  SegmentationList *s = segList[segmentationNr];
  while(s != NULL && s->ID[0] != ID)
  {
    s = s->next;
  }
  return (s != NULL);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The next segment with the same ID[0] as the current segment will be set active and the first
/// Vector of this segment will be returned.
/// If no Vector is found, NULL is returned.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Vector **FeaturePool::getFirstVectorNextSegmentList(SegmentationAdmin *admin, bool dontCareID, bool *isLastFromSegment, bool *inFilter)
{
  assert(poolLength > 1);
  SegmentationList *readPool_curSeg = admin->curSeg;
  assert(readPool_curSeg != NULL);
  if(isLastFromSegment != NULL)
  {
    *isLastFromSegment = true;
  }
  if(inFilter != NULL)
  {
    *inFilter = true;
  }

  int ID = readPool_curSeg->ID[0];
  admin->prevSeg  = readPool_curSeg;
  readPool_curSeg = readPool_curSeg->next;
  if(!dontCareID)
  {
    while(readPool_curSeg != NULL && readPool_curSeg->ID[0] != ID)
    {
      admin->prevSeg  = readPool_curSeg;
      readPool_curSeg = readPool_curSeg->next;
    }
  }
  if((readPool_curSeg != NULL) && (readPool_curSeg->firstFrame <= restrictEnd))
  {
    if(isLastFromSegment != NULL)
    {
      *isLastFromSegment = ((readPool_curSeg->firstFrame == readPool_curSeg->lastFrame) || (readPool_curSeg->firstFrame == restrictEnd));
    }
    if(inFilter != NULL && filterEnabled)
    {
      *inFilter = filter[readPool_curSeg->firstFrame];
    }
    admin->curSeg = readPool_curSeg;
    return getVectorList(admin,readPool_curSeg->firstFrame);
  }
  admin->curSeg = readPool_curSeg;
  return NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Within the current active segment, it returns the vector with offset 'offset'.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Vector **FeaturePool::getCurSegmentVectorList(SegmentationAdmin *admin, int offset)
{
  assert(buffer != NULL);
  assert(admin->curSeg != NULL);
  admin->timeStamp = offset + admin->curSeg->firstFrame;
  assert(offset >= 0 && admin->timeStamp <= admin->curSeg->lastFrame);
  return &buffer[admin->timeStamp];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

double FeaturePool::getCurSegmentEnergy(SegmentationAdmin *admin,int offset)
{
  assert(feaEnergy[0].buffer != NULL);
  assert(admin->curSeg != NULL);
  admin->timeStamp = offset + admin->curSeg->firstFrame;
  assert(offset >= 0 && admin->timeStamp <= admin->curSeg->lastFrame);
  return feaEnergy[0].buffer[admin->timeStamp]->getValue(0);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Fill gaps in a segmentation with labels labelNr.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::segmentationFillGaps(int segNr, int labelNr)
{
  assert(segNr >= 0 && segNr < MAX_NR_SEGMENTATIONS);

  SegmentationList *source = segList[segNr];

  while(source != NULL)
  {
    if((source->next != NULL) && (source->lastFrame +1 < source->next->firstFrame))
    {
       SegmentationList *add = new SegmentationList;
       add->infoString = NULL;
       add->ID[0]      = labelNr;
       add->firstFrame = source->lastFrame+1;
       add->lastFrame  = source->next->firstFrame-1;
       add->next       = source->next;
       source->next    = add;
    }
    source = source->next;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Copy a segmentation to another segmentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::segmentationCopy(int sourceSegNr, int destSegNr)
{
  assert(sourceSegNr >= 0 && sourceSegNr < MAX_NR_SEGMENTATIONS);
  assert(destSegNr >= 0 && destSegNr < MAX_NR_SEGMENTATIONS);

  resetSegmentation(destSegNr);

  SegmentationList *source = segList[sourceSegNr];

  while(source != NULL)
  {
    addSegment(destSegNr, source->firstFrame, source->lastFrame, source->ID[0],
                            source->ID[TRAIN_CONTEXTKEY],source->ID[TRAIN_CONTEXTLEFT],source->ID[TRAIN_CONTEXTRIGHT]);
    source = source->next;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// All labels in the sourceSegNr are split in two. The first half is assigned to destSegNr1 and
/// the second to destSegNr2..
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::segmentationSplit(int sourceSegNr, int destSegNr1, int destSegNr2)
{
  assert(sourceSegNr >= 0 && sourceSegNr < MAX_NR_SEGMENTATIONS);
  assert(destSegNr1 >= 0 && destSegNr1 < MAX_NR_SEGMENTATIONS);
  assert(destSegNr2   >= 0 && destSegNr2   < MAX_NR_SEGMENTATIONS);

  resetSegmentation(destSegNr1);
  resetSegmentation(destSegNr2);

  SegmentationList *source = segList[sourceSegNr];

  while(source != NULL)
  {
    int firstFrame1 = source->firstFrame;
    int lastFrame2  = source->lastFrame;
    int lastFrame1  = firstFrame1 + (lastFrame2 - firstFrame1 + 1) / 2 - 1;
    int firstFrame2 = lastFrame1 + 1;

    // If the two new segments are too short, just leave it one segment:
    if((lastFrame1 >= lastFrame2) || (firstFrame2 <= firstFrame1))
    {
      addSegment(destSegNr1, firstFrame1, lastFrame2, source->ID[0],
                            source->ID[TRAIN_CONTEXTKEY],source->ID[TRAIN_CONTEXTLEFT],source->ID[TRAIN_CONTEXTRIGHT]);
    }
    else
    {
      addSegment(destSegNr1, firstFrame1, lastFrame1, source->ID[0],
                            source->ID[TRAIN_CONTEXTKEY],source->ID[TRAIN_CONTEXTLEFT],source->ID[TRAIN_CONTEXTRIGHT]);
      addSegment(destSegNr2, firstFrame2, lastFrame2, source->ID[0],
                            source->ID[TRAIN_CONTEXTKEY],source->ID[TRAIN_CONTEXTLEFT],source->ID[TRAIN_CONTEXTRIGHT]);
    }
    source = source->next;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// It will copy the source segmentation into the destination segmentation, but only the parts
/// that are also in the filter segmentation with the filter labelNr (-1 means all labels).
/// This algorithm only works when the segments are added sequencially!
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::segmentationIntersection(int sourceSegNr, int filterSegNr, int filterLabel, int destSegNr)
{
  assert(sourceSegNr >= 0 && sourceSegNr < MAX_NR_SEGMENTATIONS);
  assert(filterSegNr >= 0 && filterSegNr < MAX_NR_SEGMENTATIONS);
  assert(destSegNr   >= 0 && destSegNr   < MAX_NR_SEGMENTATIONS);

  resetSegmentation(destSegNr);

  SegmentationList *source = segList[sourceSegNr];
  SegmentationList *filter = segList[filterSegNr];

  while(source != NULL && filter != NULL)
  {
    int firstFrame = source->firstFrame;
    int lastFrame  = source->lastFrame;
    if(firstFrame < filter->firstFrame)
    {
      firstFrame = filter->firstFrame;
    }
    if(lastFrame > filter->lastFrame)
    {
      lastFrame = filter->lastFrame;
    }
    if((firstFrame < lastFrame) && (filterLabel < 0 || filter->ID[0] == filterLabel))
    {
      addSegment(destSegNr, firstFrame, lastFrame, source->ID[0],
                            source->ID[TRAIN_CONTEXTKEY],source->ID[TRAIN_CONTEXTLEFT],source->ID[TRAIN_CONTEXTRIGHT]);
    }
    if(source->lastFrame < filter->lastFrame)
    {
      source = source->next;
    }
    else
    {
      filter = filter->next;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::segmentationCreateOverlap(int segNr, int label)
{
  SegmentationList *seg = segList[segNr];
  int lastID = -1;
  while(seg != NULL)
  {
    if(seg->ID[0] == label)
    {
      if(lastID >= 0 && lastID != label)
      {
        seg->ID[0] = lastID;
      }
      if(lastID < 0 && seg->next != NULL)
      {
        seg->ID[0] = seg->next->ID[0];
      }
      else if(seg->next != NULL && seg->next->ID[0] != lastID && seg->next->ID[0] != label)
      {
        SegmentationList *segOverlap = new SegmentationList;
        segOverlap->ID[0] = seg->next->ID[0];
        segOverlap->ID[1] = seg->next->ID[1];
        segOverlap->ID[2] = seg->next->ID[2];
        segOverlap->ID[3] = seg->next->ID[3];
        segOverlap->firstFrame = seg->firstFrame;
        segOverlap->lastFrame  = seg->lastFrame;
        segOverlap->infoString = NULL;
        segOverlap->next = seg->next;
        seg->next = segOverlap;
        seg = seg->next;
      }
    }
    lastID = seg->ID[0];
    seg = seg->next;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::segmentationSubstract(int sourceSegNr, int filterSegNr, int destSegNr, int filterLabel)
{
  assert(sourceSegNr >= 0 && sourceSegNr < MAX_NR_SEGMENTATIONS);
  assert(filterSegNr >= 0 && filterSegNr < MAX_NR_SEGMENTATIONS);
  assert(destSegNr   >= 0 && destSegNr   < MAX_NR_SEGMENTATIONS);

  resetSegmentation(destSegNr);

  SegmentationList *source = segList[sourceSegNr];
  SegmentationList *filter = segList[filterSegNr];

  bool newSet[poolLength];

  for(int i=0;i<poolLength;i++)
  {
    newSet[i] = false;
  }
  while(source != NULL)
  {
    for(int i=source->firstFrame;i<=source->lastFrame;i++)
    {
      newSet[i] = true;
    }
    source = source->next;
  }
  while(filter != NULL)
  {
    if(filterLabel < 0 || filter->ID[0] == filterLabel)
    {
      for(int i=filter->firstFrame;i<=filter->lastFrame;i++)
      {
        newSet[i] = false;
      }
    }
    filter = filter->next;
  }
  int i=0;
  int firstFrame = -1;
  source = segList[sourceSegNr];
  while(i<poolLength)
  {
    if(!newSet[i] && firstFrame >= 0)
    {
      addSegment(destSegNr, firstFrame, i, source->ID[0],
                              source->ID[TRAIN_CONTEXTKEY],source->ID[TRAIN_CONTEXTLEFT],source->ID[TRAIN_CONTEXTRIGHT]);
      firstFrame = -1;
    }
    while((source != NULL) && (i > source->lastFrame))
    {
      source = source->next;
    }
    if(newSet[i] && firstFrame < 0)
    {
      firstFrame = i;
    }
    i++;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Merge the labels label1 and label2 in the segmentation segNr (to label1).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::mergeLabels(int segNr, int label1, int label2)
{
  assert(segNr >= 0 && segNr < MAX_NR_SEGMENTATIONS);
  SegmentationList *seg = segList[segNr];

  while(seg != NULL)
  {
    if(seg->ID[0] == label2)
    {
      seg->ID[0] = label1;
    }
    seg = seg->next;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Prints the current segmentation in RTTM format to the screen (outFile = NULL), or to a file.
/// The label is placed as label in each line, the prefix is added to the segment ID.
/// If nrPr is > 0, the prefix string contains all speaker labels, it is supposed to
/// be nrSpeakers*nrPr long..
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::writeRTTM(int segmentationNr, const char *label, const char *prefix, FILE *outFile, int nrPr, int overwriteVTLN)
{
  assert(segmentationNr >= 0 && segmentationNr < MAX_NR_SEGMENTATIONS);
  UniqueIDList *idList = NULL;

  SegmentationList *readPool_curSeg = segList[segmentationNr];

  char speaker[100];

  while(readPool_curSeg != NULL)
  {
    // First check if this is a segment with a new ID:
    UniqueIDList *idItem = idList;
    int segmentID        = readPool_curSeg->ID[0];
    if(nrPr <= 0)
    {
      if(segmentID == 0)
      {
        strcpy(speaker,"SOUND");
      }
      else
      {
        sprintf(speaker,"%s%02d",prefix,segmentID);
      }
    }
    else
    {
      strncpy(speaker,&(prefix[nrPr*segmentID]),nrPr-1);
    }
    while(idItem != NULL && idItem->ID != segmentID)
    {
      idItem = idItem->next;
    }
    if(idItem == NULL)
    {
      // Found a new segment..
      idItem       = new UniqueIDList;
      idItem->ID   = segmentID;
      idItem->next = idList;
      idList       = idItem;
      if(outFile)
      {
        fprintf(outFile,"SPKR-INFO %s 1 <NA> <NA> <NA> unknown %s <NA>\n",label, speaker);
      }
      else
      {
        printf("SPKR-INFO %s 1 <NA> <NA> <NA> unknown %s <NA>\n",label, speaker);
      }
    }

    double tS = ((double)(readPool_curSeg->firstFrame))/100.0;
    double tL = ((double)(readPool_curSeg->lastFrame - readPool_curSeg->firstFrame + 1))/100.0;
    int vtlnFactor = overwriteVTLN;
    if(vtlnFactor < 0)
    {
      vtlnFactor = readPool_curSeg->ID[TRAIN_CONTEXTRIGHT];
    }
    if(outFile)
    {
      fprintf(outFile,"SPEAKER %s %d %.3f %.3f <NA> <NA> %s <NA>",label,vtlnFactor,tS,tL,speaker);
    }
    else
    {
      printf("SPEAKER %s %d %.3f %.3f <NA> <NA> %s <NA>",label,vtlnFactor,tS,tL,speaker);
    }
    if(readPool_curSeg->infoString != NULL)
    {
      if(outFile)
      {
        fprintf(outFile," %s\n",readPool_curSeg->infoString);
      }
      else
      {
        printf(" %s\n",readPool_curSeg->infoString);
      }
    }
    else
    {
      if(outFile)
      {
        fprintf(outFile,"\n");
      }
      else
      {
        printf("\n");
      }
    }
    readPool_curSeg = readPool_curSeg->next;
  }

  // Clean up the unique ID list:
  while(idList != NULL)
  {
    UniqueIDList *idItem = idList->next;
    delete idList;
    idList = idItem;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Makes all labels this segmentation a bit wider (nrFrames) except for soundNr.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::segmentationWiden(int segNr, int soundNr, int maxNrFrames, int onlyWidenLabel)
{
  SegmentationList *seg = segList[segNr];
  // First segment:
  if(seg != NULL)
  {
    seg->firstFrame -= maxNrFrames;
    if(seg->firstFrame < 0)
    {
      seg->firstFrame = 0;
    }
  }

  while(seg != NULL)
  {
    if(seg->next != NULL)
    {
      int length = seg->next->firstFrame - seg->lastFrame - 1;
      int add    = length / 2;
      if(onlyWidenLabel < 0)
      {
        if((seg->next->ID[0] == soundNr) && (seg->ID[0] == soundNr))
        {
          add = 0;
        }
        else if((seg->next->ID[0] == soundNr) || (seg->ID[0] == soundNr))
        {
          add = length;
        }
      }
      else
      {
        if(seg->ID[0] != onlyWidenLabel)
        {
          add = 0;
        }
        else if(seg->next->ID[0] != onlyWidenLabel)
        {
          add = length;
        }
      }
      if(add > maxNrFrames)
      {
        add = maxNrFrames;
      }


      if(seg->ID[0] != soundNr && ((onlyWidenLabel < 0) || (seg->ID[0] == onlyWidenLabel)))
      {
        seg->lastFrame += add;
      }
      if(seg->next->ID[0] != soundNr && ((onlyWidenLabel < 0) || (seg->next->ID[0] == onlyWidenLabel)))
      {
        seg->next->firstFrame -= add;
      }
    }
    seg = seg->next;
  }
  // The last segment is out of luck. It will not get widened.
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Filter a label-seq in a segmentation so that short segments of maximum duration maxSelf
/// neighbouring (identical) long segments of minimum duration minNeighbours are merged into the
/// neighbouring label.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::filterShortSegments(int segID, int labelNr, int labelSelf, int minNeighbours)
{
  SegmentationList *seg  = segList[segID];
  SegmentationList *prev = NULL;

  int length;
  int prevLength = minNeighbours+1;

  while(seg != NULL)
  {
    length = seg->lastFrame - seg->firstFrame + 1;
    if((seg->ID[0] == labelSelf) && (prev != NULL) && (seg->next != NULL) && (prev->ID[0] == labelNr) &&
       (prev->ID[0] == seg->next->ID[0]) && (prevLength > minNeighbours) && (seg->next->lastFrame - seg->next->firstFrame+1> minNeighbours))
    {
      length += prevLength;
      prev->next      = seg->next->next;
      prev->lastFrame = seg->next->lastFrame;
      delete seg->next;
      delete seg;
      seg = prev;
    }
    prevLength = length;
    prev       = seg;
    seg        = seg->next;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
////////////////////////////////////////////////////////////////////////////////////////////////////

int FeaturePool::filterLongSegments(int segID, int labelNr, int newLabel, int minLength)
{
  int nrChanges = 0;
  SegmentationList *seg  = segList[segID];

  while(seg != NULL)
  {
    int length = seg->lastFrame - seg->firstFrame + 1;
    if((seg->ID[0] == labelNr) && (length > minLength))
    {
      nrChanges+=length;
      seg->ID[0] = newLabel;
    }
    seg        = seg->next;
  }
  return nrChanges;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Reads a segmentation from file (RTTM format) and stores the segmentation in segment number
/// segmentationNr. ID[0] will be determined using the mapCluster. This variable contains an array
/// of strings. When one of the strings matches, the ID will be the index of the string + 1. If
/// no string matches, the ID will be zero. (the number of strings is mapNr, the length of each
/// string is strLen).
/// When prefixLength > 0, the first prefixLength characters of the RTTM 'speaker ID' will not be
/// used in the match. When mapNr < 0, instead of matching strings, the remaining string
/// will be converted to a number.
/// If fillerID > 0, all gaps in the segmentation are filled with segments of fillerID ID.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::readRTTM(int segmentationNr, char *mapCluster, int mapNr, int strLen, int prefixLength, FILE *inSegFile,
                           int fillerID, int ignoreFrames, int trainFrames, char *labelName)
{
  if(multiStream != NULL)
  {
    for(int i=0;i<infoBlock.numberOfChannels;i++)
    {
      multiStream[i].featurePool->readRTTM(segmentationNr, mapCluster, mapNr, strLen, prefixLength, inSegFile, fillerID, ignoreFrames, trainFrames, labelName);
      fseek(inSegFile,0,SEEK_SET);
    }
  }
  assert(inSegFile != NULL);
  assert(segmentationNr >= 0 && segmentationNr < MAX_NR_SEGMENTATIONS);
  SegmentationList *list = NULL;
  numberOfSpkIDs[segmentationNr] = 0;
  if(trainFrames <= 0)
  {
    trainFrames = poolLength;
  }
  if(fileOffset == 0)
  {
    resetSegmentation(segmentationNr);
  }
  else
  {
    list = segList[segmentationNr];
    if(list != NULL)
    {
      while(list->next != NULL)
      {
        list = list->next;
      }
    }
  }

  char line[MAXSTR];
  char subline[MAXSTR];
  char *label;

  int  lastFrame         = fileOffset;
  bool passedTrainFrames = false;
  while(!feof(inSegFile) && !passedTrainFrames)
  {
    if(fgets(line,MAXSTR,inSegFile) != NULL)
    {
      splitList(line,subline);
      if(strcmp(line,"SPEAKER")==0)
      {
        /// \todo add checking file format...
        // Now presuming the file is well-formed:
        SegmentationList *seg = new SegmentationList;
        seg->infoString = NULL;
        seg->next = NULL;

        // Initialise all id's:
        for(int i=0;i<4;i++)
        {
          seg->ID[i] = 0;
        }
        splitList(subline,line);  // label is in subline.
        if((labelName != NULL) && (list == NULL))
        {
          strcpy(labelName,subline);
        }
        splitList(line,subline);  // 'vtln' is in line.
        seg->ID[TRAIN_CONTEXTRIGHT] = atoi(line);
        splitList(subline,line);  // start time is in subline
        seg->firstFrame = ((int)(atof(subline)*100.0+0.5));
        splitList(line,subline);  // length is in line
        seg->lastFrame = ((int)(atof(line)*100.0+0.5)) + seg->firstFrame - 1 - 1; /// \todo This last -1 is a bug-fix: segments should not overlap...

        if(seg->lastFrame < ignoreFrames)
        {
          delete seg;
        }
        else
        {
          seg->firstFrame -= ignoreFrames;
          seg->lastFrame  -= ignoreFrames;
          if(seg->firstFrame < 0)
          {
            seg->firstFrame = 0;
          }
          if(seg->lastFrame >= trainFrames)
          {
            seg->lastFrame = trainFrames-1;
            passedTrainFrames = true;
          }
          if(seg->firstFrame >= trainFrames)
          {
            delete seg;
            passedTrainFrames = true;
          }
          else
          {
            seg->firstFrame += fileOffset;
            seg->lastFrame  += fileOffset;

            splitList(subline,line);  // <NA> is in subline
            splitList(line,subline);  // <NA> is in line
            splitList(subline,line);  // label is in subline

            // Now determine the ID:
            label = &(subline[0]);
            if(prefixLength > 0)
            {
              label = &(subline[prefixLength]);
            }
            if(mapNr > 0)
            {
              int i      = 0;
              int offset = 0;
              while(i<mapNr && strcmp(label,&(mapCluster[offset])) != 0)
              {
                i++;
                offset+=strLen;
              }
              if(i>=mapNr)
              {
                USER_ERROR2("Unknown label in RTTM file:",label);
              }
              else
              {
                seg->ID[0] = i;
              }
            }
            else
            {
              seg->ID[0] = atoi(label);
            }
            // Determine if this is a new speaker and, if so, add it to the list:
            int  i = 0;
            while((i < numberOfSpkIDs[segmentationNr]) && (spkID[segmentationNr][i] != seg->ID[0]))
            {
              i++;
            }
            if(i>=numberOfSpkIDs[segmentationNr])
            {
              if(numberOfSpkIDs[segmentationNr] >= MAX_NR_SPKIDS)
              {
                USER_ERROR("Sorry, because of a stupid implementation choice, it is not possible to have more than 10000 speakers in your RTTM.");
              }
              spkID[segmentationNr][numberOfSpkIDs[segmentationNr]] = seg->ID[0];
              numberOfSpkIDs[segmentationNr]++;
            }

            splitList(line,subline);  // <NA> is in line
            seg->infoString = new char[strlen(subline)+1];
            strcpy(seg->infoString,subline); // The info string (alignment string) is in subline.

            if((fillerID > 0) && (lastFrame + 5 < seg->firstFrame))
            {
              SegmentationList *segFiller = new SegmentationList;
              for(int i=0;i<4;i++)
              {
                segFiller->ID[i] = 0;
              }
              segFiller->infoString  = NULL;
              segFiller->ID[0]       = fillerID;
              segFiller->firstFrame  = lastFrame+1;
              segFiller->lastFrame   = seg->firstFrame - 1;
              segFiller->next        = seg;
              seg = segFiller;
            }
            if(list == NULL)
            {
              segList[segmentationNr] = seg;
            }
            else
            {
              list->next = seg;
              if(list->lastFrame >= seg->firstFrame)
              {
                PRINTF3("Error at time stamps: %d and %d (last frame, first frame)\n",list->lastFrame,seg->firstFrame);
                USER_ERROR("The RTTM file is not sorted properly. The segments should be sorted on begin time and should not overlap each other!");
              }
            }
            while(seg != NULL)
            {
              list      = seg;
              lastFrame = seg->lastFrame;
              seg       = seg->next;
            }
          }
        }
      }
    }
  }
  /// \todo This is a bug-fix so we're able to use ICSI's feature files in combination with SRI's SAD input file...
  if(list != NULL)
  {
    if(list->lastFrame >= poolLength)
    {
      list->lastFrame = poolLength - 1;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::segmentationAddLabel(int sourceSegNr, int destSegNr, int sourceLabel, int destLabel)
{
  assert(sourceSegNr >= 0 && sourceSegNr < MAX_NR_SEGMENTATIONS);
  assert(destSegNr >= 0 && destSegNr < MAX_NR_SEGMENTATIONS);

  SegmentationList *source = segList[sourceSegNr];
  SegmentationList *dest   = segList[destSegNr];

  SegmentationList *destPrev = NULL;
  // As long as dest > source:
  bool proceed = (source != NULL && dest != NULL);
  while(proceed)
  {
    while(source != NULL && source->ID[0] != sourceLabel)
    {
      source = source->next;
    }
    proceed = ((source != NULL) && (dest != NULL) && (source->lastFrame < dest->firstFrame));
    if(proceed)
    {
      SegmentationList *add = new SegmentationList;
      for(int i=1;i<4;i++)
      {
        add->ID[i] = source->ID[i];
      }
      add->infoString = NULL;
      add->ID[0]      = destLabel;
      add->firstFrame = source->firstFrame;
      add->lastFrame  = source->lastFrame;
      add->next       = dest;
      if(destPrev == NULL)
      {
        segList[destSegNr] = add;
      }
      else
      {
        destPrev->next = add;
      }
      destPrev = add;
      source = source->next;
    }
  }

  // Now the first frame in
  while(source != NULL && dest != NULL)
  {
    while(source != NULL && source->ID[0] != sourceLabel)
    {
      source = source->next;
    }
    if(source != NULL)
    {
      while(dest->next != NULL && dest->next->firstFrame < source->firstFrame)
      {
        dest = dest->next;
      }
      SegmentationList *add = new SegmentationList;
      for(int i=1;i<4;i++)
      {
        add->ID[i] = source->ID[i];
      }
      add->infoString = NULL;
      add->ID[0]      = destLabel;
      add->firstFrame = source->firstFrame;
      add->lastFrame  = source->lastFrame;
      add->next       = dest->next;
      dest->next      = add;

      dest   = dest->next;
      source = source->next;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

int FeaturePool::mergeSegmentations(int *seg, int nrSeg, int outSeg)
{
  int INITIAL_NUMBER_OF_CLUSTERS = 17;
  int mergeKey[INITIAL_NUMBER_OF_CLUSTERS][INITIAL_NUMBER_OF_CLUSTERS][INITIAL_NUMBER_OF_CLUSTERS][INITIAL_NUMBER_OF_CLUSTERS];
  int mergeFilter[poolLength];
  int nrSpkrs = 0;
  for(int i=0;i<INITIAL_NUMBER_OF_CLUSTERS;i++)
  {
    for(int i2=0;i2<INITIAL_NUMBER_OF_CLUSTERS;i2++)
    {
      for(int i3=0;i3<INITIAL_NUMBER_OF_CLUSTERS;i3++)
      {
        for(int i4=0;i4<INITIAL_NUMBER_OF_CLUSTERS;i4++)
        {
          mergeKey[i][i2][i3][i4] = -1;
        }
      }
    }
  }
  SegmentationList *mergeList[nrSeg];
  for(int i=0;i<nrSeg;i++)
  {
    mergeList[i] = segList[seg[i]];
  }
  for(int i=0;i<poolLength;i++)
  {
    int key[4];
    for(int i2=0;i2<4;i2++)
    {
      key[i2] = 0;
    }
    bool allKeyZero = true;
    for(int i2=0;i2<nrSeg;i2++)
    {
      while((mergeList[i2] != NULL) && (mergeList[i2]->lastFrame < i))
      {
        mergeList[i2] = mergeList[i2]->next;
      }
      if((mergeList[i2] == NULL) || (mergeList[i2]->firstFrame > i))
      {
        key[i2] = 0;
      }
      else
      {
        allKeyZero = false;
        key[i2]    = mergeList[i2]->ID[0];
      }
    }
    mergeFilter[i] = mergeKey[key[0]][key[1]][key[2]][key[3]];
    if((!allKeyZero) && (mergeFilter[i] == -1))
    {
      mergeKey[key[0]][key[1]][key[2]][key[3]] = nrSpkrs;
      mergeFilter[i] = nrSpkrs;
      nrSpkrs++;
    }
  }
  int totals[nrSpkrs];
  for(int i=0;i<nrSpkrs;i++)
  {
    totals[i] = 0;
  }
  for(int i=0;i<poolLength;i++)
  {
    totals[mergeFilter[i]]++;
  }


  resetSegmentation(outSeg);
  int mapSpkr[nrSpkrs];
  int lowSpkr = 1;
  for(int i=0;i<nrSpkrs;i++)
  {
    if(totals[i] > 1000)
    {
      mapSpkr[i] = lowSpkr;
      lowSpkr++;
    }
    else if(totals[i] < 100)
    {
      mapSpkr[i] = 0;
    }
    else
    {
      mapSpkr[i] = -1;
    }
  }
  int active  = mergeFilter[0];
  int oldTime = 0;
  for(int i=0;i<poolLength;i++)
  {
    if(active != mergeFilter[i])
    {
      // Check if this will be long enough (for silence):
      bool switching = true;
      if(switching)
      {
        if((active != -1) && (mapSpkr[active]>=0))
        {
          addSegment(outSeg,oldTime,i-1,mapSpkr[active]);
        }
        active = mergeFilter[i];
        oldTime = i;
      }
    }
  }
  if((active != -1) && (mapSpkr[active] >= 0))
  {
    addSegment(outSeg,oldTime,poolLength-1,mapSpkr[active]);
  }

  for(int i=0;i<nrSpkrs;i++)
  {
    PRINTF4("Nr frames for speaker SPK%02d: %12d (%d)\n",i,totals[i],mapSpkr[i]);
  }
  return lowSpkr;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::concatSegmentation(int destSeg, int copySeg)
{
  SegmentationList *seg1 = segList[destSeg];
  SegmentationList *seg2 = segList[copySeg];
  if(seg1 == NULL)
  {
    segList[destSeg] = segList[copySeg];
  }
  else if(seg2 != NULL)
  {
    while(seg1->next != NULL)
    {
      seg1 = seg1->next;
    }
    if(seg2->lastFrame < seg1->lastFrame)
    {
      PRINTF5("Last and first: (%d,%d) : (%d,%d)\n",seg1->firstFrame, seg1->lastFrame, seg2->firstFrame, seg2->lastFrame);
      USER_ERROR("Can not concatinate these segments: the two segmentations do not overlap properly");
    }
    if(seg1->ID[0] == seg2->ID[0])
    {
      seg1->lastFrame = seg2->lastFrame;
      seg1->next      = seg2->next;
      delete seg2;
    }
    else
    {
      seg1->next = seg2;
    }
  }
  fillPool_curSeg[destSeg] = fillPool_curSeg[copySeg];
  fillPool_curSeg[copySeg] = NULL;
  segList[copySeg]         = NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FeaturePool::restrictPool(int start, int end)
{
  if(start < 0 || end < 0)
  {
    restrictStart = 0;
    restrictEnd   = poolLength - 1;
  }
  else
  {
    restrictStart = start;
    restrictEnd   = end;
    if(restrictEnd > poolLength - 1)
    {
      restrictEnd = poolLength - 1;
    }
    PRINTF3("Start and End: %d %d\n",restrictStart,restrictEnd);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool FeaturePool::createUEMSegmentation(int segID, int channelNumber, char *label, char *UEM)
{
  resetSegmentation(segID);
  bool res = false;
  char filter[MAXSTR];
  sprintf(filter,"%s %%s %%s %%s",label);
  char str[MAXSTR];
  FILE *uemFile = fopen(UEM,"rt");
  if(uemFile == NULL)
  {
    USER_ERROR("Could not open the UEM file.");
  }
  char beginTime[100];
  char endTime[100];
  char channel[100];
  int  lastEndFrame = -1;
  while(!feof(uemFile))
  {
    if(fscanf(uemFile,filter,channel,beginTime,endTime) != 0 && !feof(uemFile))
    {
      // Assuming 10ms windows:
      int channelNr  = atoi(channel);
      int beginFrame = ((int)(atof(beginTime) * 100.0));
      int endFrame   = ((int)(atof(endTime)   * 100.0)) - 1;
      if((channelNumber < 0) || (channelNumber == channelNr))
      {
        if(beginFrame <= lastEndFrame)
        {
          PRINTF5("Beginframe, lastEndFrame, channelNr (channelNumber) = %d,%d,%d (%d)\n",beginFrame,lastEndFrame,channelNr,channelNumber);
          USER_ERROR("The UEM is corrupt: please make sure the segments are sorted assending in time");
        }
        res = true;
        if(channelNumber < 0)
        {
          addSegment(segID,beginFrame,endFrame,0);
        }
        else
        {
          addSegment(segID,beginFrame,endFrame,channelNumber);
        }
        lastEndFrame = endFrame;
      }
    }
    else
    {
      char *resChar = fgets(str,MAXSTR,uemFile);
    }
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

int FeaturePool::getNumberOfSpeakers(int segNr)
{
  assert(segNr < MAX_NR_SEGMENTATIONS);
  return numberOfSpkIDs[segNr];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

int FeaturePool::getNumberOfSegments(int segNr, int ID)
{
  assert(segNr < MAX_NR_SEGMENTATIONS);
  int count = 0;

  SegmentationAdmin readPool;
  Vector *v = getFirstVectorFirstSegment(&readPool, segNr, ID);
  while(v != NULL)
  {
    count++;
    v = getFirstVectorNextSegment(&readPool);
  }
  return count;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

int FeaturePool::getSpeakerID(int segNr, int spkrNr)
{
  assert(segNr < MAX_NR_SEGMENTATIONS);
  assert(spkrNr < MAX_NR_SPKIDS);
  return spkID[segNr][spkrNr];
}



