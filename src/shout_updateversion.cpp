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
#include "shout_updateversion.h"
#include "segmenter.h"
#include "featurepool.h"
#include "shout-misc.h"
using namespace WriteFileLittleBigEndian;


#define NR_CHOISES (3)

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Handles the command parameters for FeatureExtraction. It is the main function of the 
/// feature extraction application. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  if(argc != 6)
  {
    printf("Usage: %s amName audio rttm clusterAM resultFile\n",(char*)argv[0]);
  }
  Shout_UpdateVersion Shout_UpdateVersion((char*)argv[1],(char*)argv[2],(char*)argv[3],(char*)argv[4],(char*)argv[5]);
  return 0;
}


Shout_UpdateVersion::Shout_UpdateVersion(char *amName, char *audio, char *rttm, char *amClusterName, char *resFile)
{
  // ++++++++++++++++++++++++ CHECK INPUT ++++++++++++++++++++++++++++++++
/*
  FILE *modelFile = fopen(amName, "rb");
  if(modelFile==NULL)
  {
    USER_ERROR("No valid AM file passed!");
  }
  FILE *rttmFile   = fopen(rttm, "rt");
  if(rttmFile==NULL)
  {
    USER_ERROR("No valid META file passed!");
  }
  FILE *modelFileCluster = fopen(amClusterName, "rb");
  if(modelFileCluster==NULL)
  {
    USER_ERROR("No valid AM file passed!");
  }

  // +++++++++++++++++++++++++++++ LOAD FEATURES +++++++++++++++++++++++++

  FeaturePool *featurePool = new FeaturePool(audio,FILE_TYPE_RAW_DECODING,false,1.0,false,true);  // Only prepare...
  int segID = featurePool->claimSegmentationID();
  featurePool->readRTTM(segID,NULL,-1,-1,3,rttmFile);
  fclose(rttmFile);
  featurePool->createNewPool(audio, NULL, NULL, false, false, segID);

  // ++++++++++++++++++++++++ CREATE SEGMENTER +++++++++++++++++++++++++++

  Segmenter *segmenter = new Segmenter(modelFile);
  segmenter->createLexicalTree(3);
  segmenter->setFeaturePool(featurePool);

  int outSeg = featurePool->claimSegmentationID();

// ++++++++++++++++++++++++ INITIALIZE MAPPINGS ++++++++++++++++++++++++++


int nrOfClusters = featurePool->getNumberOfSpeakers(segID);
int mapping[nrOfClusters][NR_CHOISES];
int phoneMap[segmenter->getNumberOfModels()+1];

for(int i=0;i<segmenter->getNumberOfModels()+1;i++)
{
  phoneMap[i] = 0;
}
for(int i=0;i<nrOfClusters;i++)
{
  for(int i2=0;i2<NR_CHOISES;i2++)
  {
    mapping[i][i2] = -1;
  }
}


// ++++++++++++++++++++++++ LET'S MAP! +++++++++++++++++++++++++++++++++++

  segmenter->overwritePrunePars(true,true,true, 20.0, 1.0, 1.0, false, 1, 500);
  for(int i=0;i<nrOfClusters;i++)
  {
    int id = featurePool->getSpeakerID(segID,i);
    if(featurePool->getClusterLength(segID,id) > 10)
    {
      segmenter->mapSegmentation(segID,id,outSeg);
      printf("MODEL-MAP %03d:", id);
      int score[segmenter->getNumberOfModels()+1];
      for(int i2=0;i2<segmenter->getNumberOfModels()+1;i2++)
      {
        score[i2] = featurePool->getNumberOfSegments(outSeg,i2);
      }
      for(int aha=0;aha<NR_CHOISES;aha++)
      {
        int max   = 0;
        int res   = -1;
        int total = 0;
        for(int i2=0;i2<segmenter->getNumberOfModels()+1;i2++)
        {
          int nr = score[i2];
          if(nr > max)
          {
            max = nr;
            res = i2;
          }
          total += nr;
        }
        printf(" %02d",res-1);
        mapping[i][aha] = res - 1;
        phoneMap[res-1] = -1;
        score[res] = -1;
      }
      printf("\n"); fflush(stdout);
    }
    else
    {
      printf("MODEL-MAP %03d: none\n", id);
    }
  }


  segmenter->overwritePrunePars(true,true,true, 2000.0, 1.0, 1.0, false, 1, 5);
  for(int i2=2;i2<segmenter->getNumberOfModels()+1;i2++)
  {
    if(phoneMap[i2-1] == 0)
    {
      int winner      = -1;
      int winnerModel = -1;
      float max  = -9.0e31;
      for(int i=0;i<nrOfClusters;i++)
      {
          int id = featurePool->getSpeakerID(segID,i);
          float res = segmenter->modelScore(segID,id,i2);
          if(max < res)
          {
            max         = res;
            winner      = id;
            winnerModel = i;
          }
      }
      printf("PHONE-MAP %02d: %d (%d) %f\n",i2-1,winner,winnerModel,max); fflush(stdout);
      phoneMap[i2-1] = winnerModel;
    }
  }

  // ++++++++++++++++++++++ STORE THE EXTERNAL AM +++++++++++++++++++++++++

printf("Storing external AM scores...\n"); fflush(stdout);

  delete segmenter;

  int vectorSize = 0;
  int nrClusters = 0;
  int notUsed    = 0;
  freadEndianSafe(&vectorSize,1,sizeof(vectorSize),modelFileCluster);
  freadEndianSafe(&nrClusters,1,sizeof(nrClusters),modelFileCluster);
  freadEndianSafe(&notUsed,1,sizeof(notUsed),modelFileCluster);
  PhoneModel **models = new PhoneModel*[nrClusters];
  models[0] = NULL;
  for(int i=1;i<nrClusters;i++)
  {
    models[i] = new PhoneModel(modelFileCluster,vectorSize);
  }

  modelFile = fopen(amName, "rb");
  if(modelFile==NULL)
  {
    USER_ERROR("No valid AM file passed!");
  }
  int nrOfModels;
  freadEndianSafe(&notUsed,1,sizeof(notUsed),modelFile);
  freadEndianSafe(&nrOfModels,1,sizeof(nrOfModels),modelFile);
  freadEndianSafe(&notUsed,1,sizeof(notUsed),modelFile);
  PhoneModel *silModel = new PhoneModel(modelFile,vectorSize);
  fclose(modelFile);

  FILE *externalAM = fopen(resFile, "wb");
  int poolLength = featurePool->getPoolLength();
  SegmentationAdmin readPool_curSeg;
  Vector **vector = featurePool->getVectorList(&readPool_curSeg,0);
  while(vector != NULL)
  {
    float res = silModel->getLogPDFProbability(0,vector[0]);
    fwriteEndianSafe(&res,1,sizeof(float),externalAM);
    vector = featurePool->getNextVectorList(&readPool_curSeg);
  }

  int nrMappedModels = 0;
  for(int i=1; i<nrOfModels;i++)
  {
    if(phoneMap[i] < 0)
    {
      nrMappedModels++;
      printf("Writing data for model %d (combined)..\n",i);
    }
    else
    {
      printf("Writing data for model %d (fixed)..\n",i);
    }
    vector = featurePool->getVectorList(&readPool_curSeg,0);
    int time = 0;
    while(vector != NULL)
    {
      float res = 0.0;
      float nrM = 0.0;
      if(phoneMap[i] > 0)
      {
        res = models[phoneMap[i]]->getLogPDFProbability(0,vector[0]);
      }
      else
      {
        for(int a=1;a<nrOfClusters;a++)
        {
          for(int a2=0;a2<NR_CHOISES;a2++)
          {
            if(mapping[a][a2] == i)
            {
              nrM += (((float)(NR_CHOISES - a2)) / ((float)NR_CHOISES));
              res += models[a]->getPDFProbability(0,vector[0],0,time);
            }
          }
        }
        res = log(res / nrM);
        if(isnan(res) || res < MINIMUM_LOGP)
        {
          res = MINIMUM_LOGP;
        }
      }
      fwriteEndianSafe(&res,1,sizeof(float),externalAM);
      vector = featurePool->getNextVectorList(&readPool_curSeg);
      time++;
    }
  }
  fclose(externalAM);
  printf("Done, there are %d combined models!\n",nrMappedModels);

  // +++++++++++++++++++++++++++++ DONE ++++++++++++++++++++++++++++++++++
  for(int i=1;i<nrClusters;i++)
  {
    delete models[i];
  }
  delete[] models;
  delete silModel;
  delete featurePool;
  printf("Bye bye!\n");
*/
}


Shout_UpdateVersion::~Shout_UpdateVersion()
{
}


