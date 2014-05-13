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
#include "shout_vtln.h"
#include "lexicaltree.h"
#include "featurepool.h"
#include "phonemodel.h"
#include "shout-misc.h"
#include "featureextraction.h"

using namespace StringFunctions;
using namespace WriteFileLittleBigEndian;

#define MAX_STRING  (5000)


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

Shout_VTLN::Shout_VTLN(char *method, char *AMFile, char *audioRttmIn, char *audioName, char *audioRttmOut, double *totals)
{
  lastSpeakerVTLN = VTLN_STEPS / 2;

  SegmentationAdmin segmentList;
  int vectorSize       = VTLN_DEFAULT_MFCC_NUMBER+VTLN_DEFAULT_MFCC_ENERGY + VTLN_DEFAULT_MFCC_ZEROCROSS +
                         2*VTLN_DEFAULT_MFCC_DELTA*(VTLN_DEFAULT_MFCC_NUMBER+VTLN_DEFAULT_MFCC_ENERGY+VTLN_DEFAULT_MFCC_ZEROCROSS);
  int numberOfModels   = 0;
  int numberOfClusters = 0;
  int useModel         = -1;
  PhoneModel **models  = NULL;

  if(audioRttmOut != NULL)
  {
    PRINTF("# Reading acoustic models...\n"); fflush(stdout);
  }
  // Creating the acoustic models:  
  FILE *modelFile = fopen(AMFile, "rb");
  if(modelFile==NULL)
  {
    USER_ERROR("No valid AM file passed\n");
  }
  int vectSize;
  freadEndianSafe(&vectSize,1,sizeof(vectSize),modelFile);
  if(vectorSize != vectSize)
  {
    USER_ERROR("This AM file contains models with incompatible vector size!");
  }
  freadEndianSafe(&numberOfModels,1,sizeof(numberOfModels),modelFile);
  freadEndianSafe(&numberOfClusters,1,sizeof(numberOfClusters),modelFile);
  if(numberOfClusters > 1)
  {
    USER_WARNING("This AM contains clusters. Only the first cluster will be used!\n");
  }
  if(numberOfModels <= 0 || numberOfClusters <= 0)
  {
    USER_ERROR("This AM file does not contain any acoustic models!\n");
  }
  else
  {
    models = new PhoneModel*[numberOfModels];
    for(int i=0;i<numberOfModels;i++)
    {
      models[i] = new PhoneModel(modelFile,vectorSize,true);
      if(audioRttmOut != NULL)
      {
        PRINTF3("# model %02d: %s\n",i,models[i]->getStatistics()->name); fflush(stdout);
      }
      if(strncmp(models[i]->getStatistics()->name,"SPEECH    ",6)==0)
      {
        useModel = i;
      }
    }
    if(audioRttmOut != NULL)
    {
      PRINTF2("# Using model %02d\n",useModel); fflush(stdout);
    }
    if(useModel < 0)
    {
      USER_ERROR("Sorry, no 'SPEECH' AM model in this file.");
    }
    fclose(modelFile);
    modelFile = NULL;
  }

  FeaturePool *featurePool[VTLN_STEPS];
  double vtln = VTLN_START;
  char label[MAX_STRING];
  FILE *listFile;
  int segID   = 0;
  for(int vtcount=0;vtcount<VTLN_STEPS;vtcount++)
  {
    if(audioRttmOut != NULL)
    {
      PRINTF2("Creating feature pool %d..\n",vtcount); fflush(stdout);
    }
    featurePool[vtcount] = new FeaturePool(audioName,FILE_TYPE_RAW_VTLN,false,vtln,false,true);  // Only prepare...
    listFile = fopen(audioRttmIn,"rt");
    int s = featurePool[vtcount]->claimSegmentationID();
    if(vtcount == 0)
    {
      segID = s;
    }
    featurePool[vtcount]->readRTTM(s,NULL,-1,-1,3,listFile,-1,0,0,label);
    fclose(listFile);
    featurePool[vtcount]->createNewPool(audioName, NULL, NULL, false, false, s, vtln);  // Segments are used to find SAD...
    vtln+=VTLN_STEPSIZE;
  }

  if(strcmp(method,"SEGMENT") == 0)
  {
    Vector *vector = featurePool[0]->getFirstVectorFirstSegment(&segmentList,segID,-1,true);
    while(vector != NULL)
    {
      int bestVtln = VTLN_STEPS/2 + 1;
      int time    = featurePool[0]->getCurSegmentStart(&segmentList);
      int length  = featurePool[0]->getCurSegmentLen(&segmentList);
      double bestScore = -9.0e300;
      for(int vtcount=0;vtcount<VTLN_STEPS;vtcount++)
      {
        double like = 0.0;
        for(int i=time;i<time+length;i++)
        {
          like += models[useModel]->getLogPDFProbability(0,featurePool[vtcount]->getVector(NULL,i));
        }
        if(like > bestScore)
        {
          bestScore = like;
          bestVtln  = vtcount;
        }
      }
      featurePool[0]->setSegmentID(TRAIN_CONTEXTRIGHT,bestVtln,&segmentList);
      vector = featurePool[0]->getFirstVectorNextSegment(&segmentList,true);
    }

    if(audioRttmOut != NULL)
    {
      listFile = fopen(audioRttmOut,"wt");
      featurePool[0]->writeRTTM(segID, label, "SPK", listFile, 0);
      fclose(listFile);
    }
  }
  else if(strcmp(method,"SPEAKER") == 0)
  {
    for(int spk=0;spk<featurePool[0]->getNumberOfSpeakers(segID);spk++)
    {
      int spkID = featurePool[0]->getSpeakerID(segID,spk);
      int bestVtln = VTLN_STEPS/2 + 1;
      double bestScore = -9.0e300;
      for(int vtcount=0;vtcount<VTLN_STEPS;vtcount++)
      {
        double like = 0.0;
        if(totals != NULL)
        {
          like = totals[vtcount];
        }
        Vector *vector = featurePool[0]->getFirstVectorFirstSegment(&segmentList,segID,spkID,false);
        while(vector != NULL)
        {
          int time    = featurePool[0]->getCurSegmentStart(&segmentList);
          int length  = featurePool[0]->getCurSegmentLen(&segmentList);
          for(int i=time;i<time+length;i++)
          {
            like += models[useModel]->getLogPDFProbability(0,featurePool[vtcount]->getVector(NULL,i));
          }
          vector = featurePool[0]->getFirstVectorNextSegment(&segmentList,false);
        }
        if(audioRttmOut != NULL)
        {
          PRINTF3("Score using VTLN factor number %02d\t%f\n",vtcount,like);
        }
        if(totals != NULL)
        {
          totals[vtcount] = like;
        }
        if(like > bestScore)
        {
          bestScore = like;
          bestVtln  = vtcount;
        }
      }
      lastSpeakerVTLN = bestVtln;
      // Set the ID for each segment:
      Vector *vector = featurePool[0]->getFirstVectorFirstSegment(&segmentList,segID,spkID,false);
      while(vector != NULL)
      {
        featurePool[0]->setSegmentID(TRAIN_CONTEXTRIGHT,bestVtln,&segmentList);
        vector = featurePool[0]->getFirstVectorNextSegment(&segmentList,false);
      }
    }
    if(audioRttmOut != NULL)
    {
      listFile = fopen(audioRttmOut,"wt");
      featurePool[0]->writeRTTM(segID, label, "SPK", listFile, 0);
      fclose(listFile);
    }
  }
  else
  {
    USER_ERROR("Sorry, this VTLN method is not yet implemented... Use 'SEGMENT' for now.\n");
  }

/*
  for(int i=0;i<numberOfModels;i++)
  {
    if(models[i] != NULL)
    {
      delete models[i];
    }
    delete[] models;
  }
*/
//   for(int vtcount=0;vtcount<VTLN_STEPS;vtcount++)
//   {
//     if(featurePool[vtcount] != NULL)
//     {
//       delete featurePool[vtcount];
//       featurePool[vtcount] = NULL;
//     }
//   }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

int Shout_VTLN::getVtlnFactor()
{
  return lastSpeakerVTLN;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

Shout_VTLN::~Shout_VTLN()
{
}

