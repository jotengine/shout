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
#include "segmenter.h"
#include "phonemodel.h"
#include "shout-misc.h"

using namespace WriteFileLittleBigEndian;


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor initialises the cluster parameters and sets the LexicalTree parameters for 
/// speaker clustering (defined in standard.h).
/////////////////////////////////////////////////////////////////////////////////////////////////////

Segmenter::Segmenter() : LexicalTree(stdout)
{
  discrTrain       = NULL;
  discrTrainMask   = NULL;
  numberOfClusters = 0;
  overwritePrunePars(CLUSTER_PERFORM_PRUNING_HISTOGRAM,CLUSTER_PERFORM_PRUNING_STATE,CLUSTER_PERFORM_PRUNING_ENDSTATE,
                     CLUSTER_BEAM, CLUSTER_STATE_BEAM, CLUSTER_END_STATE_BEAM, CLUSTER_LMLA,
                     CLUSTER_HISTOGRAM_STATE_PRUNING, CLUSTER_HISTOGRAM_PRUNING);
  overwriteWeightPars(CLUSTER_LM_SCALE,CLUSTER_TRANS_PENALTY,CLUSTER_SIL_PENALTY);
  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor initialises the cluster parameters and sets the LexicalTree parameters for 
/// speaker clustering (defined in standard.h). It also loads an existing clustering from inFile.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Segmenter::Segmenter(FILE *inFile)  : LexicalTree(stdout)
{
  discrTrain       = NULL;
  discrTrainMask   = NULL;
  numberOfClusters = 0;

  if(inFile != NULL)
  {
    int vectorSize = 0;
    freadEndianSafe(&vectorSize,1,sizeof(vectorSize),inFile);
    freadEndianSafe(&numberOfClusters,1,sizeof(numberOfClusters),inFile);
    int notUsed = 0;
    freadEndianSafe(&notUsed,1,sizeof(notUsed),inFile);
    assert(notUsed == 1);
    numberOfPhones = numberOfClusters+1; // admin... Model 0 has special meaning...
    phoneModels = new PhoneModel*[numberOfPhones];
    phoneModels[0] = NULL;
    for(int i=1;i<numberOfPhones;i++)
    {
      phoneModels[i] = new PhoneModel(inFile,vectorSize);
    }
    fclose(inFile);
  }
  overwritePrunePars(CLUSTER_PERFORM_PRUNING_HISTOGRAM,CLUSTER_PERFORM_PRUNING_STATE,CLUSTER_PERFORM_PRUNING_ENDSTATE,
                     CLUSTER_BEAM, CLUSTER_STATE_BEAM, CLUSTER_END_STATE_BEAM, CLUSTER_LMLA,
                     CLUSTER_HISTOGRAM_STATE_PRUNING, CLUSTER_HISTOGRAM_PRUNING);
  overwriteWeightPars(CLUSTER_LM_SCALE,CLUSTER_TRANS_PENALTY,CLUSTER_SIL_PENALTY);  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Given the current number of models, create the lexical tree that will be used to align the data.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Segmenter::createLexicalTree(int minNumberOfFrames, int *forEachFrames1, double *forEachFrames2)
{
  assert(minNumberOfFrames > 2);

  // First delete the old tree:
  initialiseSystem();
  deleteTree(grammarStart,false);

  grammarStart = NULL;
  // Now create the new tree!
  currentlyAligning = true;

  for(int i=1;i<numberOfPhones;i++)
  {
    if(forEachFrames1 != NULL)
    {
      minNumberOfFrames = forEachFrames1[i];
    }

    if(phoneModels[i] != NULL)
    {
      phoneModels[i]->resetPhoneAdmin();

      LexicalNode *n = new LexicalNode;
      n->compressedTreeIndex          = 0;
      if(grammarStart == NULL)
      {
        grammarStart                  = n;
        grammarStart->parallel        = NULL;
        grammarStart->inputToken      = new (TokenType*);
        *grammarStart->inputToken     = NULL;
        startOfSentenceWord           = i;
        endOfSentenceWord             = i;
      }
      else
      {
        n->inputToken                 = grammarStart->inputToken;
        n->parallel                   = grammarStart->parallel;
        grammarStart->parallel        = n;
      }
      n->tokenSeqLength               = minNumberOfFrames;
      n->tokenSeq                     = new TokenType*[minNumberOfFrames];
      for(int a=0;a<minNumberOfFrames;a++)
      {
        n->tokenSeq[a]                = NULL;
      }
      if(forEachFrames2 != NULL)
      {
        phoneModels[i]->setSilTrans(forEachFrames2[i]);
      }
/*      else
      {
        phoneModels[i]->setSilTrans(log(0.1 / (double)(numberOfPhones-1)));
      }
*/
      n->modelID                      = i;
      n->contextPrev                  = 0;
      n->contextNext                  = 0;
      n->contextKey                   = 0;
      n->depth                        = -1;
      n->next                         = NULL;
      n->nextTree                     = grammarStart;
      n->wordID                       = i;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor clears any used memory.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Segmenter::~Segmenter()
{
  if(phoneModels != NULL)
  {
    for(int i=0;i<numberOfPhones;i++)
    {
      if(phoneModels[i] != NULL)
      {
        delete phoneModels[i];
      }
    }
    delete[] phoneModels; 
  }
  if(discrTrain != NULL)
  {
    for(int i=0;i<discrTrainLen;i++)
    {
      delete discrTrain[i];
    }
    delete[] discrTrain;
  }
  if(discrTrainMask != NULL)
  {
    for(int i=0;i<numberOfPhones;i++)
    {
      delete discrTrainMask[i];
    }
    delete[] discrTrainMask;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Given the current feature pool and input segmentation and label to use, a new output segmentation
/// is created. If the inputSegID is smaller than zero, the entire featurePool is used as one segment.
/// If inputLabelID is smaller than zero, all labels in the input segmentation are used.
/////////////////////////////////////////////////////////////////////////////////////////////////////

float Segmenter::segmentFeaturePool(int inputSegID, int inputLabelID, int outputSegID)
{
  assert(featurePool != NULL);
  SegmentationAdmin readPool_curSeg;
  initialiseTree(featurePool->startTimeStamp());

  float res  = 0.0;
  if(inputSegID >= 0)
  {
    int tempID   = featurePool->claimSegmentationID();
    bool isLast = false;
    Vector **vector = featurePool->getFirstVectorFirstSegmentList(&readPool_curSeg, inputSegID,inputLabelID,inputLabelID < 0, &isLast);

    int nrSamples  = 0;
    while(vector != NULL)
    {
      int time       = featurePool->getCurSegmentStart(&readPool_curSeg);
      processVector(vector,time);

      
      while(!isLast)
      {
        vector = featurePool->getNextVectorList(&readPool_curSeg,&isLast);        
        if(vector != NULL)
        {
          time++;

          processVector(vector,time);
          nrSamples++;
        }
      }
      vector = featurePool->getFirstVectorNextSegmentList(&readPool_curSeg,inputLabelID < 0, &isLast);
    }
    featurePool->resetSegmentation(tempID);
    res = createSegments(tempID,&readPool_curSeg);

    featurePool->segmentationIntersection(tempID,inputSegID,inputLabelID,outputSegID);
    featurePool->freeSegmentationID(tempID);
  }
  else
  {
    int time = featurePool->startTimeStamp();
    Vector **vector = featurePool->getVectorList(&readPool_curSeg,time);
    while(vector != NULL)
    {
      processVector(vector,time);
      time++;
      vector = featurePool->getNextVectorList(&readPool_curSeg);
    }
    res = createSegments(outputSegID,&readPool_curSeg);
  }
  initialiseTree(featurePool->startTimeStamp());
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Uses decode result to create some segments, (part of) the pool should have been processed.
/////////////////////////////////////////////////////////////////////////////////////////////////////

float Segmenter::createSegments(int segID, SegmentationAdmin *admin)
{
  float totScore = -9.0e300;
  TokenType *useToken = findBestToken(false,false,false);
  if(useToken != NULL)
  {
    totScore = useToken->likelihood;
PRINTF2("TOTSCORE: %g\n",totScore);

    WLRType *path = useToken->path;
    while(path != NULL && path->previous != NULL)
    {
      if(path->timeStamp >= 0)
      {
        int start  = path->previous->timeStamp+1;
        if(start < 0)
        {
          start = 0;
        }
        int winner = 0;
        if(!path->isSil)
        {
          winner = path->lmHistory[0];
        }
        featurePool->addSegmentToFront(segID,start,path->timeStamp,winner);
      }
      path = path->previous;
    }
  }
  else
  {
    printf("WARNING: too short input stream! I will give the entire stream to one model.\n");
    // This fragment was too small to begin with... Let's find the one with the highest score and give it to that cluster:
    // Each first token will contain a value...
    totScore = -9.0e300;
    int winner     = -1;
    LexicalNode *l = grammarStart;
    while(l != NULL)
    {
      if(l->tokenSeq[0] != NULL && l->tokenSeq[0]->likelihood > totScore)
      {
        totScore = l->tokenSeq[0]->likelihood;
        winner   = l->modelID;
      }
      l = l->parallel;
    }
    if(winner >= 0)
    {
      featurePool->addSegment(segID,0,featurePool->getCurSegmentLen(admin)-1,winner);    
    }
    else
    {
      printf("  WARNING: Could not find a winner for this segment...\n");
    }    
  }
  return totScore;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Sets the feature pool.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Segmenter::setFeaturePool(FeaturePool *fp)
{
  featurePool = fp;
}
