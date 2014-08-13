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
#include "shoutsegment.h"
#include "featurepool.h"
#include "adapt_am_treenode.h"
#include "trainphonemodel.h"
#include "shoutconfig.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern int MINPASSES;

#define USE_GARBAGE (true)
#define FAST_MERGE  (false)

#define MIN_DATA    (100)


#ifdef MAX_NR_THREADS
#include <pthread.h>

struct Thread_Train_Data
{
  TrainPhoneModel *model;
  int id;
  int nr;
};

Thread_Train_Data  trainThread[MAX_NR_THREADS];
bool               shoutsegment_threadRunning[MAX_NR_THREADS];
bool               shoutsegment_threadMayStart[MAX_NR_THREADS];
pthread_mutex_t    condition_shoutSegment_mutexStart[MAX_NR_THREADS];
pthread_cond_t     condition_shoutSegment_threadStart[MAX_NR_THREADS];
pthread_t          shoutsegment_thread[MAX_NR_THREADS];
pthread_mutex_t    condition_shoutsegment_mutexDone;
pthread_cond_t     condition_shoutsegment_threadDone;



void *thread_train( void *ptr )
{
  Thread_Train_Data *data = ((Thread_Train_Data*)ptr);
  while(true)
  {
    pthread_mutex_lock( &condition_shoutSegment_mutexStart[data->id]);
    while(!shoutsegment_threadMayStart[data->id])
    {
      pthread_cond_wait( &condition_shoutSegment_threadStart[data->id], &condition_shoutSegment_mutexStart[data->id]);
    }
    shoutsegment_threadMayStart[data->id] = false;
    pthread_mutex_unlock( &condition_shoutSegment_mutexStart[data->id] );


    PRINTF3("Thread %d is training (%p)...\n",data->id,data->model); fflush(stdout);
    data->model->train(data->nr,true);
    pthread_mutex_lock( &condition_shoutsegment_mutexDone);
    shoutsegment_threadRunning[data->id] = false;
    pthread_cond_signal( &condition_shoutsegment_threadDone);
    pthread_mutex_unlock( &condition_shoutsegment_mutexDone);
  }
}
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Helper function for threads, so that we don't need compile switches all over my code.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void ShoutSegment::trainModel(TrainPhoneModel *model, int nr, bool waitForMe)
{
PRINTF3("Request for training a model %p (waiting=%d)\n",model,waitForMe); fflush(stdout);
#ifdef MAX_NR_THREADS
  int freeID = -1;
  // First find a free thread:
  pthread_mutex_lock( &condition_shoutsegment_mutexDone );
  while(freeID < 0)
  {
    int i = 0;
    while(shoutsegment_threadRunning[i] && i < MAX_NR_THREADS)
    {
      i++;
    }
    if(i< MAX_NR_THREADS)
    {
      freeID = i;
    }
    if(freeID < 0)
    {
      pthread_cond_wait( &condition_shoutsegment_threadDone, &condition_shoutsegment_mutexDone );
    }
  }
  assert(freeID >= 0 && freeID < MAX_NR_THREADS);
  shoutsegment_threadRunning[freeID] = true;
  PRINTF2("  Assigned ID %d.\n",freeID);fflush(stdout);
  pthread_mutex_unlock( &condition_shoutsegment_mutexDone );

  trainThread[freeID].model = model;
  trainThread[freeID].nr    = nr;

  pthread_mutex_lock( &condition_shoutSegment_mutexStart[freeID] );
  shoutsegment_threadMayStart[freeID] = true;
  pthread_cond_signal( &condition_shoutSegment_threadStart[freeID] );
  pthread_mutex_unlock( &condition_shoutSegment_mutexStart[freeID] );


  // If we must wait, let's wait until all threads are done.
  pthread_mutex_lock( &condition_shoutsegment_mutexDone );
  while(waitForMe)
  {
    int i = 0;
    while(!shoutsegment_threadRunning[i] && i < MAX_NR_THREADS)
    {
      i++;
    }
    waitForMe = (i < MAX_NR_THREADS);
    if(waitForMe)
    {
      pthread_cond_wait( &condition_shoutsegment_threadDone, &condition_shoutsegment_mutexDone );
    }
  }
  pthread_mutex_unlock( &condition_shoutsegment_mutexDone );
PRINTF("Done!\n");fflush(stdout);
#else
  model->train(nr,true);
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Handles the command parameters for ShoutSegment. It is the main function of the
/// shout_segment application.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_SEGMENT, argc, argv);

  FILE *modelFile = fopen(myConfig.getStringValue(ARGUMENT_AM_SEGMENT),"rb");

  if(( myConfig.parIsDefined(ARGUMENT_CTS_IN) && !myConfig.parIsDefined(ARGUMENT_CTS_OUT)) ||
     (!myConfig.parIsDefined(ARGUMENT_CTS_IN) &&  myConfig.parIsDefined(ARGUMENT_CTS_OUT)))
  {
    USER_ERROR("For CTS segmentation, define the second input AND second output file.");
  }

  if(modelFile == NULL && !myConfig.parIsDefined(ARGUMENT_CTS_IN))
  {
    USER_ERROR("Could not open the acoustic model file");
  }
  char *uemFile = NULL;
  if(myConfig.parIsDefined(ARGUMENT_UEM))
  {
    uemFile = myConfig.getStringValue(ARGUMENT_UEM);
  }
#ifdef MAX_NR_THREADS
  int max = 3;
  PRINTF("Max threads defined\n");
  if(MAX_NR_THREADS < 3)
  {
    max = MAX_NR_THREADS;
  }
  for(int i=0;i<max;i++)
  {
    PRINTF2("Inside thread loop %i\n",i);
    shoutsegment_threadRunning[i]   = false;
    shoutsegment_threadMayStart[i]  = false;
    trainThread[i].id  = i;
    pthread_create( &shoutsegment_thread[i], NULL, thread_train, (void*) &trainThread[i]);
    PRINTF("Outside of thread creation...\n");
    fflush(stdout);
  }
#endif
  ShoutSegment ShoutSegment(myConfig.getStringValue(ARGUMENT_AUDIOFILE), modelFile, myConfig.getStringValue(ARGUMENT_META_OUT),
                            myConfig.getStringValue(ARGUMENT_LABEL), uemFile,
                            myConfig.parIsDefined(ARGUMENT_WIDEN_SEGMENTS),
                            myConfig.getStringValue(ARGUMENT_CTS_IN),
                            myConfig.getStringValue(ARGUMENT_CTS_OUT),
                            myConfig.parIsDefined(ARGUMENT_ENERGY_BASED));
  PRINTF("ShoutSegment created\n...");
  fflush(stdout);
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

#define FIRST_RUN          (5)
#define SECOND_RUN         (5)

#define MIN_GARBAGE         (75)
#define MIN_SPEECH          (75)
#define FINAL_MIN_GARBAGE  (200)
#define FINAL_MIN_SPEECH   (200)
#define SEGMENT_WINDOWSIZE (60000)

#define MIN_GARBAGE_DATA   (1600)

#define SAMPLE_PER_GAUSSIAN (2000)

ShoutSegment::ShoutSegment(char *featureFile, FILE *amFile, char *rttmOut, char *label, char *UEM, bool widenSegments,
                           char *ctsInName, char *ctsOutName, bool energyBased) : Segmenter(amFile)
{
  if(widenSegments)
  {
    USER_ERROR("Sorry, widening segments after SAD is not implemented yet!");
  }

  PRINTF("Creating ShoutSegment...\n");
  fflush(stdout);
  SegmentationAdmin segmentList;
  bool performCTS = (ctsInName != NULL && strlen(ctsInName) > 0);

  char classNames[80];
  strncpy(  classNames,     "SIL"   ,20);
  strncpy(&(classNames[20]),"SIL"   ,20);
  strncpy(&(classNames[40]),"SPEECH",20);
  strncpy(&(classNames[60]),"SOUND" ,20);
  MINPASSES = 10;
  int forEachModel1[4]     = {1,30,MIN_SPEECH,MIN_GARBAGE};
  PRINTF("Reading input file..\n");
  fflush(stdout);
  featurePool = new FeaturePool(featureFile,FILE_TYPE_RAW_SAD);
  if(!energyBased && (phoneModels[1]->dim() != featurePool->getVectorSize()))
  {
    USER_ERROR("Phonemodel dimensionality and feature file do not match!");
  }

  featurePool->retrieveEnergy(featureFile,1);

  if(energyBased || featurePool->getPoolLength() < 1000)
  {
    PRINTF("Starting energy-based segmentation.\n");
    int   id;
    bool  content = false;
    id = featurePool->claimSegmentationID();
    featurePool->addSegment(id,0,featurePool->getPoolLength()-1,1,1,1,1);
    featurePool->removeLowEnergy(id,1,true);

    // int id2 = featurePool->claimSegmentationID();
    // featurePool->addSegment(id2,0,featurePool->getPoolLength()-1,0);
    // int id3 = featurePool->claimSegmentationID();
    // featurePool->segmentationSubstract(id2,id,id3);

    // Done. Let's write to file...
    FILE *rttm = fopen(rttmOut,"wt");
    if(rttm == NULL)
    {
      USER_ERROR("Could not write to the output RTTM file for the left channel.");
    }
    featurePool->writeRTTM(id,label,"SPEECH",rttm,0,VTLN_NEUTRAL);
    fclose(rttm);
  }
  else if(performCTS)
  {
    PRINTF("Starting segmentation for CTS!\n");
    featurePool->retrieveEnergy(ctsInName,2);
    int   ctsID[2];
    bool  ctsContent[2] = {false,false};
    if(UEM == NULL)
    {
      USER_ERROR("For CTS segmentation, a UEM is needed");
    }
    for(int i=0;i<2;i++)
    {
      ctsID[i] = featurePool->claimSegmentationID();
      ctsContent[i] = featurePool->createUEMSegmentation(ctsID[i],i+1,label,UEM);
      if(ctsContent[i])
      {
        featurePool->removeLowEnergy(ctsID[i],i+1);
      }
    }

    if(ctsContent[0])
    {
      // Done. Let's write to file...
      FILE *rttm = fopen(rttmOut,"wt");
      if(rttm == NULL)
      {
        USER_ERROR("Could not write to the output RTTM file for the left channel.");
      }
      featurePool->writeRTTM(ctsID[0],label,"SPK",rttm,0,VTLN_NEUTRAL);
      fclose(rttm);
    }
    if(ctsContent[1])
    {
      FILE *rttm = fopen(ctsOutName,"wt");
      if(rttm == NULL)
      {
        USER_ERROR("Could not write to the output RTTM file for the right channel.");
      }
      featurePool->writeRTTM(ctsID[1],label,"SPK",rttm,0,VTLN_NEUTRAL);
      fclose(rttm);
    }
  }
  else
  {
    PRINTF("Done reading!\n"); fflush(stdout);
    uemID = featurePool->claimSegmentationID();
    if(UEM != NULL)
    {
      featurePool->createUEMSegmentation(uemID,-1,label,UEM);
    }
    else
    {
      featurePool->resetSegmentation(uemID);
      featurePool->addSegment(uemID,0,featurePool->getPoolLength()-1,0);
    }

    int finalID    = featurePool->claimSegmentationID();
    int segID      = featurePool->claimSegmentationID();
    int silID      = featurePool->claimSegmentationID();
    int speechID   = featurePool->claimSegmentationID();
    int garID      = featurePool->claimSegmentationID();
    int tempID     = featurePool->claimSegmentationID();

    featurePool->resetSegmentation(finalID);

    if(numberOfPhones != 3)
    {
      USER_ERROR("The acoustic model file must contain two models: SIL and SPEECH!");
    }

    for(int i=1;i<numberOfPhones;i++)
    {
      PRINTF4("Model %d (%d gaussians): %s\n",i,phoneModels[i]->getNumberOfGaussians(),phoneModels[i]->getStatistics()->name);
    }

    PhoneModel *orgSil     = phoneModels[1];
    PhoneModel *orgSpeech  = phoneModels[2];

    delete[] phoneModels;
    numberOfClusters = 3;
    numberOfPhones   = 4;
    phoneModels      = new PhoneModel*[4];
    phoneModels[0]   = NULL;
    phoneModels[1]   = orgSil;
    phoneModels[2]   = orgSpeech;
    phoneModels[3]   = NULL;
    createLexicalTree(100, forEachModel1, NULL);


    // Splitting in SEGMENT_WINDOWSIZE samples (approx)...
    int totalLength = 0;
    int chunkLength = 0;
    int nrChunks    = 0;
    totalLength = featurePool->getClusterLength(uemID,0);
    double part     = ((double)totalLength) / ((double)SEGMENT_WINDOWSIZE);
    nrChunks    = ((int)part);
    if(nrChunks <= 0)
    {
      nrChunks = 1;
    }
    else if(part - ((double)nrChunks) > 0.5)
    {
      nrChunks++;
    }
    chunkLength = totalLength / nrChunks + 1;

    // Now find the other positions:

    int i     = 0;
    int count = 0;
    featurePool->resetSegmentation(tempID);
    Vector *vector = featurePool->getFirstVectorFirstSegment(&segmentList,uemID,0,false);
    while(vector != NULL)
    {
      int length  = featurePool->getCurSegmentLen(&segmentList);
      int tot     = featurePool->getCurSegmentStart(&segmentList);
      int ind     = 0;
      int indPrev = 0;
      while(ind+1 < length)
      {
        if((length - ind) < (chunkLength - count))
        {
          // Easy does it, use the entire sentence:
          int lframe = tot+length-ind-1;
          if(lframe - tot > 0)
          {
            featurePool->addSegment(tempID,tot,lframe,i);
          }
          tot   += (length - ind);
          count += (length - ind);
          ind    = length;
        }
        else
        {
          int lframe = tot+(chunkLength - count)-1;
          if(lframe - tot > 0)
          {
            featurePool->addSegment(tempID,tot,lframe,i);
          }
          tot  += (chunkLength - count);
          ind  += (chunkLength - count);
          count = chunkLength;
        }
        indPrev = ind;
        if(count >= chunkLength)
        {
          count = 0;
          if(i+1 < nrChunks)
          {
            i++;
          }
        }
      }
      vector = featurePool->getFirstVectorNextSegment(&segmentList,false);
    }

    featurePool->resetSegmentation(uemID);
    featurePool->segmentationCopy(tempID,uemID);

    // Just some logging:
    vector = featurePool->getFirstVectorFirstSegment(&segmentList, uemID,-1,true);
    while(vector != NULL)
    {
      int length  = featurePool->getCurSegmentLen(&segmentList);
      int tot     = featurePool->getCurSegmentStart(&segmentList);
      int id      = featurePool->getSegmentID(0,&segmentList);
      PRINTF4("UEM: Start: %d End: %d ID: %d\n",tot,tot+length-1,id);
      vector = featurePool->getFirstVectorNextSegment(&segmentList,true);
    }
    PRINTF2("Number of chunks: %d\n",nrChunks);
    PRINTF2("Chunk size      : %d\n",chunkLength);
    PRINTF2("Pool length     : %d\n",featurePool->getPoolLength());


    for(int chunks = 0;chunks < nrChunks;chunks++)
    {
      PRINTF2("Handling chunk number %d\n",chunks); fflush(stdout);
      createLexicalTree(100, forEachModel1, NULL);
      PRINTF("Segmenting pool: first run.\n"); fflush(stdout);
      segmentFeaturePool(uemID,chunks,segID);
      featurePool->segmentationCopy(segID,speechID);

      TrainPhoneModel *speechModel  = NULL;
      TrainPhoneModel *silModel     = NULL;
      TrainPhoneModel *garbageModel = NULL;

      if(USE_GARBAGE)
      {
        int silG    = 2;
        int garG    = 2;
        int speechG = 6;

        if(silModel != NULL)
        {
          delete silModel;
        }
        if(garbageModel != NULL)
        {
          delete garbageModel;
        }
        PRINTF("Training models...\n"); fflush(stdout);
        silModel     = new TrainPhoneModel("SIL",1,1,true,featurePool->getVectorSize());
        garbageModel = new TrainPhoneModel("SOUND",1,1,true,featurePool->getVectorSize());;

        for(int i=0;i<FIRST_RUN;i++)
        {
          featurePool->segmentationCopy(segID,silID);
          featurePool->resetSegmentation(segID);
          featurePool->segmentationIntersection(silID,speechID,1,segID);

          PRINTF5("First run iteration %d.. %d - %d - %d\n",i,silG,garG,speechG);
          int silSize = 0;
          int garSize = 0;

          if(i<3)
          {
            featurePool->mergeLabels(segID,1,3);
            featurePool->maxSegmentLength(segID,1,100,50);
            featurePool->splitOnEnergy(segID,1,1,3,(i+1)*SAMPLE_PER_GAUSSIAN,5*SAMPLE_PER_GAUSSIAN,ENERGY_COMP);
            featurePool->splitOnEnergy(segID,3,0,3,SAMPLE_PER_GAUSSIAN,(i+1)*SAMPLE_PER_GAUSSIAN,ZEROCROSS_COMP);
          }
          featurePool->segmentationCopy(segID,garID);
          featurePool->segmentationCopy(segID,silID);
         // featurePool->segmentationSubset(segID,silID,1,phoneModels[1],(i+1)*SAMPLE_PER_GAUSSIAN);

          garSize = featurePool->getClusterLength(garID,3);
          silSize = featurePool->getClusterLength(silID,1);
          PRINTF2("Size of models: SIL = %d\n",silSize);
          PRINTF2("                GAR = %d\n",garSize);

          if(garSize <= MIN_DATA)
          {
            delete garbageModel;
            garbageModel   = NULL;
          }
          if(silSize < MIN_DATA)
          {
            delete silModel;
            silModel       = NULL;
          }

          if(silModel != NULL)
          {
            silModel->setTrainingData(featurePool,silID,1);
            trainModel(silModel,silG,garbageModel == NULL);
          }
          if(garbageModel != NULL)
          {
            garbageModel->setTrainingData(featurePool,garID,3);
            trainModel(garbageModel,garG,true);
          }

          phoneModels[1] = ((PhoneModel*)silModel);
          phoneModels[3] = ((PhoneModel*)garbageModel);

          createLexicalTree(100, forEachModel1, NULL);
          featurePool->resetSegmentation(segID);
          PRINTF2("Segmenting pool: run %d.\n",i);
          segmentFeaturePool(uemID,chunks,segID);
          if(i<3)
          {
            garG+=2;
          }
        }

        featurePool->segmentationCopy(segID,silID);
        featurePool->resetSegmentation(garID);
        featurePool->segmentationIntersection(silID,speechID,2,garID);

        if(featurePool->getClusterLength(garID,2) > MIN_DATA)
        {
          featurePool->resetSegmentation(segID);
          featurePool->segmentationCopy(garID,segID);
          featurePool->resetSegmentation(garID);

          if(speechModel != NULL)
          {
            delete speechModel;
          }
          speechModel  = new TrainPhoneModel("SPEECH",1,1,true,featurePool->getVectorSize());
          speechModel->setTrainingData(featurePool,segID,2);
          trainModel(speechModel,speechG,true);
          phoneModels[2] = ((PhoneModel*)speechModel);
          createLexicalTree(100, forEachModel1, NULL);

          featurePool->resetSegmentation(segID);
          segmentFeaturePool(uemID,chunks,segID);

          for(int i=0;i<SECOND_RUN;i++)
          {

            speechG+=2;
            garG   +=2;
            silG   ++;

            PRINTF5("Second run iteration %d.. %d - %d - %d\n",i,silG,garG,speechG);

            if(garbageModel != NULL)
            {
              if(featurePool->getClusterLength(segID,3) < MIN_DATA)
              {
                PRINTF("No data for SOUND model, deleting...\n");
                delete garbageModel;
                garbageModel   = NULL;
                phoneModels[3] = NULL;
                createLexicalTree(100, forEachModel1, NULL);
              }
              else
              {
                garbageModel->setTrainingData(featurePool,segID,3);
                trainModel(garbageModel,garG,(silModel==NULL && speechModel == NULL));
              }
            }
            if(silModel != NULL)
            {
              silModel->setTrainingData(featurePool,segID,1);
              trainModel(silModel,silG,(speechModel == NULL));
            }
            if(speechModel != NULL)
            {
              speechModel->setTrainingData(featurePool,segID,2);
              trainModel(speechModel,speechG,true);
            }
            featurePool->resetSegmentation(segID);
            segmentFeaturePool(uemID,chunks,segID);
          }
          // Okay, now look if we should keep this garbage model:
          if(garbageModel != NULL)
          {
            garbageModel->setTrainingData(featurePool,segID,3);
            speechModel->setTrainingData(featurePool,segID,2);
            TrainPhoneModel *merge = new TrainPhoneModel(garbageModel,speechModel);
            merge->setTrainingData(featurePool,segID,2,3);
            double mergeScore = merge->train(-1,true,true) - garbageModel->getTrainSilP(3) - speechModel->getTrainSilP(2);
            delete merge;
            PRINTF2("===============>>>>>>>  check: %f\n",mergeScore);
            if(mergeScore > 0.0)
            {
              if(FAST_MERGE)
              {
                PRINTF("Using 'fast merge' to merge the sound model and the speech model..\n");
                featurePool->mergeLabels(segID,2,3);
              }
              else
              {
                PRINTF("Using 'extended merge' to merge the sound model and the speech model..\n");
                // Sorry, have to start over again....
                phoneModels[1] = orgSil;
                phoneModels[2] = orgSpeech;
                phoneModels[3] = NULL;
                createLexicalTree(100, forEachModel1, NULL);
                featurePool->resetSegmentation(segID);
                segmentFeaturePool(uemID,chunks,segID);
                doSAD_noGarbage(featurePool,segID,silID,speechID,chunks);
              }
            }
          }
        }
        else
        {
          // We are assuming that there MUST be speech.. So the garbage model sucked all the speech out of the pool.
          // Let's start over without the garbage model. The real garbage will go into the silence model.
          phoneModels[1] = orgSil;
          phoneModels[2] = orgSpeech;
          phoneModels[3] = NULL;
          createLexicalTree(100, forEachModel1, NULL);
          featurePool->resetSegmentation(segID);
          segmentFeaturePool(uemID,chunks,segID);
          doSAD_noGarbage(featurePool,segID,silID,speechID,chunks);
        }
      }
      else
      {
        doSAD_noGarbage(featurePool,segID,silID,speechID,chunks);
      }

      PRINTF("Starting post-processing...\n");
     // createLexicalTree(100, forEachModelFinal, NULL);
     // featurePool->resetSegmentation(segID);
     // segmentFeaturePool(uemID,chunks,segID);
     // featurePool->filterShortSegments(segID,2,1,MIN_SPEECH-1);
      PRINTF("  filter short SIL segments.\n");
      featurePool->filterShortSegments(segID,3,1,MIN_GARBAGE-1);

      PRINTF("  delete models.\n");
      assert(silModel != orgSil);
      assert(speechModel != orgSpeech);
      if(garbageModel != NULL)
      {
        delete garbageModel;
        garbageModel   = NULL;
      }
      if(silModel != NULL)
      {
        delete silModel;
        silModel = NULL;
      }
      if(speechModel != NULL)
      {
        delete speechModel;
        speechModel = NULL;
      }
      PRINTF("  Concatinate segmentations.\n");
      featurePool->concatSegmentation(finalID,segID);
      if(rttmOut == NULL || (strlen(rttmOut) == 0))
      {
        PRINTF("  write result to STDOUT.\n");
        featurePool->writeRTTM(finalID,label,classNames,NULL,20,VTLN_NEUTRAL);
      }

      PRINTF("  reset models to original...\n");
      phoneModels[1] = orgSil;
      phoneModels[2] = orgSpeech;
      phoneModels[3] = NULL;
    }


    if(rttmOut != NULL && (strlen(rttmOut) > 0))
    {
      PRINTF("  write result to file.\n"); fflush(stdout);
      FILE *rttm = fopen(rttmOut,"wt");
      if(rttm == NULL)
      {
        USER_ERROR("Could not write to the output RTTM file.");
      }
      featurePool->writeRTTM(finalID,label,classNames,rttm,20,VTLN_NEUTRAL);
      fclose(rttm);
    }

    PRINTF("  delete original (SIL/SPEECH) models.\n");

    if(orgSil != NULL)
    {
      delete orgSil;
    }
    if(orgSpeech != NULL)
    {
      delete orgSpeech;
    }
    phoneModels[0] = NULL;
    phoneModels[1] = NULL;
    phoneModels[2] = NULL;
    phoneModels[3] = NULL;
    createLexicalTree(100, forEachModel1, NULL);
  }

  PRINTF("Deleting pool\n");  fflush(stdout);
  delete featurePool;
  PRINTF("Done.\n");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

ShoutSegment::~ShoutSegment()
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

void ShoutSegment::doSAD_noGarbage(FeaturePool *featurePool, int segID, int silID, int speechID, int chunks)
{
  int silG    = 2;
  int speechG = 2;
  MINPASSES   = 8;

  TrainPhoneModel *silModel     = new TrainPhoneModel("SIL",   1,1,true,featurePool->getVectorSize());
  TrainPhoneModel *speechModel  = new TrainPhoneModel("SPEECH",1,1,true,featurePool->getVectorSize());
  featurePool->segmentationSubset(segID,silID,   1,phoneModels[1],200*SAMPLE_PER_GAUSSIAN);
  featurePool->segmentationCopy(segID,speechID);

  silModel->setTrainingData(featurePool,silID,1);
  speechModel->setTrainingData(featurePool,speechID,2);
  int silSize    = featurePool->getClusterLength(silID,1);
  int speechSize = featurePool->getClusterLength(speechID,2);
  PRINTF3("SIL AND SPEECH: %d  %d\n",silSize,speechSize);
  trainModel(silModel,silG,false);
  trainModel(speechModel,speechG,true);
  phoneModels[1] = ((PhoneModel*)silModel);
  phoneModels[2] = ((PhoneModel*)speechModel);

  featurePool->resetSegmentation(segID);
  segmentFeaturePool(uemID,chunks,segID);
  if(featurePool->getClusterLength(segID,2) > MIN_DATA)
  {
    featurePool->segmentationSubset(segID,silID,   1,phoneModels[1],  SAMPLE_PER_GAUSSIAN);
    featurePool->segmentationSubset(segID,speechID,2,phoneModels[2],  SAMPLE_PER_GAUSSIAN);
    delete silModel;
    delete speechModel;
    silModel     = new TrainPhoneModel("SIL",   1,1,true,featurePool->getVectorSize());
    speechModel  = new TrainPhoneModel("SPEECH",1,1,true,featurePool->getVectorSize());
    silModel->setTrainingData(featurePool,silID,1);
    speechModel->setTrainingData(featurePool,speechID,2);
    trainModel(silModel,silG,false);
    trainModel(speechModel,speechG,true);
    phoneModels[1] = ((PhoneModel*)silModel);
    phoneModels[2] = ((PhoneModel*)speechModel);
    featurePool->resetSegmentation(segID);
    segmentFeaturePool(uemID,chunks,segID);


    for(int i2=0;i2<7;i2++)
    {
      silModel->setTrainingData(featurePool,segID,1);
      speechModel->setTrainingData(featurePool,segID,2);
      trainModel(silModel,silG,false);
      trainModel(speechModel,speechG,true);

      featurePool->resetSegmentation(segID);
      segmentFeaturePool(uemID,chunks,segID);

      if(i2<5)
      {
        speechG+=2;
      }
      if(i2<3)
      {
        silG++;
      }
    }
  }
  MINPASSES = 10;
}
