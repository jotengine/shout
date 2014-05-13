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
#include <sys/time.h>

#include "standard.h"
#include "shout_cluster.h"
#include "shout-misc.h"
#include "shoutconfig.h"

using namespace StringFunctions;

#include <math.h>

extern int SEGMENT_INITIAL_DATA_SEGMENTS;
// extern int MINPASSES;
extern int MERGE_METHOD;
extern int STOP_METHOD;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This is the main function for the shout_cluster application.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_CLUSTER, argc, argv);

  char *feature2 = NULL;
  if(myConfig.parIsDefined(ARGUMENT_SECOND_STREAM))
  {
    feature2 = myConfig.getStringValue(ARGUMENT_SECOND_STREAM);
  }
  char *feature3 = NULL;
  if(myConfig.parIsDefined(ARGUMENT_THIRD_STREAM))
  {
    feature3 = myConfig.getStringValue(ARGUMENT_THIRD_STREAM);
  }
  int mergeAtaTime = 1;
  if(myConfig.parIsDefined(ARGUMENT_FAST_MERGE))
  {
    mergeAtaTime = (int)(myConfig.getFloatValue(ARGUMENT_FAST_MERGE));
    if(mergeAtaTime <= 0)
    {
      USER_ERROR("The 'fast merge' number must be 1 or higher.");
    }
  }
  char fileSAD2[MAXSTR];
  strcpy(fileSAD2,"");

  if(myConfig.parIsDefined(ARGUMENT_META_IN))
  {
    strcpy(fileSAD2,myConfig.getStringValue(ARGUMENT_META_IN));
    strcat(fileSAD2,".train");
  }
  int maxNumberOfClusters = -1;
  if(myConfig.parIsDefined(ARGUMENT_MAX_CLUSTERS))
  {
    maxNumberOfClusters = (int)(myConfig.getFloatValue(ARGUMENT_MAX_CLUSTERS));
  }
  Shout_Cluster Shout_Cluster(myConfig.getStringValue(ARGUMENT_AUDIOFILE),
                              myConfig.getStringValue(ARGUMENT_LABEL),
                              fileSAD2,
                              myConfig.getStringValue(ARGUMENT_META_IN),
                              myConfig.getStringValue(ARGUMENT_META_OUT),
                              feature2,
                              0.9,
                              feature3,
                              mergeAtaTime,
                              maxNumberOfClusters,
                              myConfig.parIsDefined(ARGUMENT_WIDEN_SEGMENTS));
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor reads a feature file in order to segment and cluster it.
///
/// Energy based Speaker Activity Detection is used for a first rough speech/non-speech clustering.
/// After that, the method in Train_Speaker_Segmenter is used to obtain a speaker segmentation.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Shout_Cluster::Shout_Cluster(char *featureFile, char *label, char *inputRTTM_train, char *inputRTTM_decode, 
                             char *outputRTTM, char *featureFile2, double initWeight, char *featureFile3,
                             int fastMerge, int maxClusters, bool widen)
{
  FeaturePool_Multistream_file files[3];

  char refCheat[MAXSTR];
  refCheat[0] = 0;

  strcpy(refCheat, inputRTTM_decode);
  strcat(refCheat,".ref");

  FILE *test = fopen(inputRTTM_train,"rt");
  if(test == NULL)
  {
    inputRTTM_train = inputRTTM_decode;
  }
  else
  {
    fclose(test);
  }

  files[0].fileName = featureFile;
  files[0].fileType = FILE_TYPE_RAW_DIARIZATION;
  PRINTF2("Feature file: %s\n",featureFile);
  FILE *testFile = fopen(featureFile2,"rb");
  bool useSecondStream = false;
  bool useThirdStream  = false;
  if(testFile != NULL)
  {
    PRINTF2("Second feature file: %s\n",featureFile2);
    useSecondStream = true;
    fclose(testFile);
    files[1].fileName = featureFile2;
    files[1].fileType = FILE_TYPE_TEXT;

    testFile = fopen(featureFile3,"rb");
    if(testFile != NULL)
    {
      PRINTF2("Third feature file: %s\n",featureFile3);
      useThirdStream = true;
      fclose(testFile);
      files[2].fileName = featureFile3;
      files[2].fileType = FILE_TYPE_TEXT;
    }
    else
    {
      PRINTF("Not using a third feature file.\n");
    }
  }
  else
  {
    PRINTF("Not using a second feature file.\n");
  }

  PRINTF2("Label: %s\n",label);
  PRINTF2("RTTM for training: %s\n",inputRTTM_train);
  PRINTF2("RTTM for decoding: %s\n",inputRTTM_decode);
  PRINTF2("RTTM for cheating: %s\n",refCheat);
  PRINTF2("Output RTTM: %s\n",outputRTTM);

  PRINTF("Loading feature vector file...\n");
  int len = strlen(featureFile);
  FeaturePool *feaPool = NULL;
  FILE_TYPE fileType = FILE_TYPE_SHOUT;
  if((featureFile[len-3] == 'r' && featureFile[len-2] == 'a' && featureFile[len-1] == 'w') ||
     (featureFile[len-3] == 'R' && featureFile[len-2] == 'A' && featureFile[len-1] == 'W')    )
  {
    PRINTF("Raw file detected...\n");
    fileType = FILE_TYPE_RAW_DIARIZATION;
  }
  if(useThirdStream)
  {
    double weight[3] = {0.75,0.25,0.05};
    feaPool = new FeaturePool(3,files,weight,true);
  }
  else if(useSecondStream)
  {
    double weight[2] = {initWeight,1.0-initWeight};
    feaPool = new FeaturePool(2,files,weight,true);
  }
  else
  {
    feaPool = new FeaturePool(featureFile,fileType,false,1.0,false,true);
  }

  char filterName[2000];
  strcpy(filterName,featureFile);
  strcat(filterName,".asr");  
  
//  feaPool->setFilter(filterName,FEATURE_SPEAKER_SIZE,0.8);  
  feaPool->enableFilter(false);
  
  int sadID_train  = feaPool->claimSegmentationID();
  int sadID_decode = feaPool->claimSegmentationID();
  int segID_cheat  = -1;


//  segID_cheat = maxClusters;
//  maxClusters = 40;


  PRINTF2("Done, total number of frames: %d\n",feaPool->getPoolLength());
  
  // Use the input RTTM file to create segments:
  PRINTF("Reading SAD RTTM file for training...\n");
  FILE *inSegFile = fopen(inputRTTM_train,"rt");
  if(inSegFile == NULL)
  {
    USER_ERROR("ERROR: Could not read the input RTTM file for training!\n");
  }
  char segmentString[60];
  strcpy(segmentString,"SPEECH");
  strcpy(&(segmentString[20]),"SIL");
  strcpy(&(segmentString[40]),"SOUND");
  feaPool->readRTTM(sadID_train,segmentString,3,20,0,inSegFile);
  fclose(inSegFile);


  PRINTF("Reading SAD RTTM file for decoding...\n");
  inSegFile = fopen(inputRTTM_decode,"rt");
  if(inSegFile == NULL)
  {
    USER_ERROR("ERROR: Could not read the input RTTM file for decoding!\n");
  }
  feaPool->readRTTM(sadID_decode,segmentString,3,20,0,inSegFile);
  fclose(inSegFile);
  
  PRINTF("Reading reference RTTM file for cheating...\n");
  inSegFile = NULL;
  if(refCheat != NULL && strlen(refCheat) > 0)
  {
    inSegFile = fopen(refCheat,"rt");
  }
  if(inSegFile == NULL)
  {
    PRINTF("Could not open reference file. I guess we're not going to cheat today...\n");
  }
  else
  {
    segID_cheat = feaPool->claimSegmentationID();
    feaPool->readRTTM(segID_cheat,NULL,-1,-1,3,inSegFile);
    fclose(inSegFile);
  }


  PRINTF("Creating features...\n"); fflush(stdout);
  if(useSecondStream)
  {
    feaPool->createNewMultiStreamPool(files, sadID_train);
  }
  else
  {
    feaPool->createNewPool(featureFile, NULL, NULL, false, false, sadID_train);
  }


  // Now train!
  PRINTF("Start training..\n");

//  segmenter = new Train_Segmenter(feaPool,sadID_train,sadID_decode,label, fastMerge, 14, 18, widen);
  segmenter = new Train_Segmenter(feaPool,sadID_train,sadID_decode,label, fastMerge, 16, maxClusters, widen, NULL);
  segmenter->train(0, label, outputRTTM, NULL, segID_cheat);
  delete segmenter;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor is empty.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Shout_Cluster::~Shout_Cluster()
{
}


