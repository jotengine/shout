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

#include <time.h>
#include <sys/time.h>
#include <math.h>

#include "standard.h"
#include "train_segmenter.h"
#include "languagemodel_segmenter.h"
#include "shout-misc.h"

using namespace WriteFileLittleBigEndian;

#define absolute(x)      (((x)<(0.0))?(-1.0*x):(x))

#define LOGGING_ON

//#define ICASSP2010GMM
//int CLUSTER_MIN_NUMBER_FRAMES = 4;

// Perfect topology:
// #define DIARIZATION_CHEAT_EXP1
// Speech activity detection (only change *.rttm (decode))
// #define DIARIZATION_CHEAT_EXP1
// Merging algorithm (also change *.rttm.train):
//#define DIARIZATION_CHEAT_EXP3
// Model initialization
//#define DIARIZATION_CHEAT_EXP4
// Merge candidate selection:
//#define DIARIZATION_CHEAT_EXP5
// Stop criterion (normal system...)
// none..
// Perfect initialization, but further normal...
//#define DIARIZATION_CHEAT_EXP7
// Perfect initialization, no re-alignment...
//#define DIARIZATION_CHEAT_EXP8
// Do not train with false alarms...
// #define DIARIZATION_CHEAT_NOFA

// #define GLOBAL_STOP

#define MINPASSES_ALIGN  (3)

// #define STARTWITHLAST

int SEGMENT_INITIAL_DATA_SEGMENTS  = 2;

extern int MINPASSES;
extern int nrOfStateHistPruning;

#if (defined DIARIZATION_CHEAT_EXP3) || (defined DIARIZATION_CHEAT_EXP4) || (defined DIARIZATION_CHEAT_EXP5) || (defined DIARIZATION_CHEAT_EXP7) || (defined GLOBAL_STOP)
int globalIt     = 1;
#endif

int MERGE_METHOD = 0;
int STOP_METHOD  = 0;

//#define SAD_TEST



#ifdef MAX_NR_THREADS
#include <pthread.h>

struct Thread_Train_Data_Cluster
{
  TrainPhoneModel *model;
  double          *trainPRes;
  int segid;
  int id;
  int nr;
};

Thread_Train_Data_Cluster  trainThreadCluster[MAX_NR_THREADS];
bool               cluster_threadRunning[MAX_NR_THREADS];
bool               cluster_threadMayStart[MAX_NR_THREADS];
pthread_mutex_t    condition_cluster_mutexStart[MAX_NR_THREADS];
pthread_cond_t     condition_cluster_threadStart[MAX_NR_THREADS];
pthread_t          cluster_thread[MAX_NR_THREADS];
pthread_mutex_t    condition_cluster_mutexDone;
pthread_cond_t     condition_cluster_threadDone;



void *thread_train_cluster( void *ptr )
{
  Thread_Train_Data_Cluster *data = ((Thread_Train_Data_Cluster*)ptr);
  while(true)
  {
    pthread_mutex_lock( &condition_cluster_mutexStart[data->id]);
    while(!cluster_threadMayStart[data->id])
    {
      pthread_cond_wait( &condition_cluster_threadStart[data->id], &condition_cluster_mutexStart[data->id]);
    }
    cluster_threadMayStart[data->id] = false;
    pthread_mutex_unlock( &condition_cluster_mutexStart[data->id] );

    data->model->train(data->nr,(CLUSTER_MIN_NUMBER_FRAMES != 3),true,NULL,NULL);
    *data->trainPRes = data->model->getTrainSilP(data->segid);

    pthread_mutex_lock( &condition_cluster_mutexDone);
    cluster_threadRunning[data->id] = false;
    pthread_cond_signal( &condition_cluster_threadDone);
    pthread_mutex_unlock( &condition_cluster_mutexDone);
  }
}
#endif


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor will initialise a LexicalTree topology with (hopefully) too many clusters.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Train_Segmenter::Train_Segmenter(FeaturePool *fp, int sID_train, int sID_decode, const char *name, int nrMerges, int minClusters, int maxClusters, bool widen, const char *nonMergeableModels) : Segmenter()
{
  maxDataPoints  = 0;

  for(int i=0;i<10;i++)
  {
    bic_meanoffset[i] = 0.0;
  }

  prepareForASR = widen;
  if(nrMerges < 1)
  {
    nrMerges = 1;
  }
  fastMerge = nrMerges;

  srand((unsigned)time(0));  
  int numberOfVectors = fp->getClusterLength(sID_decode,0);
  PRINTF2("Total number of speech frames: %d\n",numberOfVectors);
  int initialNrClusters = (int) floor(0.5 + ((double) numberOfVectors) / ((double) (INITIAL_DATA_PER_GAUSSIAN * MIN_NUMBER_OF_GAUSSIANS))) + 1;

  if(initialNrClusters < minClusters+1)
  {
    initialNrClusters = minClusters+1;
  }
  if(initialNrClusters > maxClusters+1)
  {
    initialNrClusters = maxClusters+1;
  }
//  initialNrClusters = 17;
  PRINTF2("Number of initial clusters: %d\n",initialNrClusters);

  assert(fp != NULL);
  setFeaturePool(fp);

  nrNonMergeModels = 0;
  trcl = NULL;
  if(nonMergeableModels != NULL)
  {
    FILE *noMergeFile = fopen(nonMergeableModels,"rb");
    if(noMergeFile != NULL)
    {
      bool proceed = true;
      int temp;
      freadEndianSafe(&temp,1,sizeof(temp),noMergeFile);
      PRINTF2("-> %d\n",temp);
      freadEndianSafe(&temp,1,sizeof(temp),noMergeFile);
      PRINTF2("-> %d\n",temp);
      freadEndianSafe(&temp,1,sizeof(temp),noMergeFile);
      PRINTF2("-> %d\n",temp);
      while(proceed)
      {
        char str[11];
        strcpy(str,"FIXED");
        trcl = new TrainPhoneModel(str,1,1,true, featurePool->getVectorSize(),featurePool->getStreamInfo());
        proceed = trcl->readModel(noMergeFile);
        if(proceed)
        {
          nrNonMergeModels++;
        }
        delete trcl;
        trcl = NULL;
      }
      fclose(noMergeFile);
    }
    PRINTF2("NUMBER OF EXTRA MODELS: %d\n",nrNonMergeModels);
  }

  mergeable = new bool[initialNrClusters+nrNonMergeModels];
  for(int i=0;i<initialNrClusters;i++)
  {
    mergeable[i] = true;
  }
  for(int i=initialNrClusters;i<initialNrClusters+nrNonMergeModels;i++)
  {
    mergeable[i] = false;
  }
  initialNrClusters += nrNonMergeModels;

  distanceVect   = NULL;
  name = name; // Not used..
  numberOfMerges = 0;
  assert(initialNrClusters >= 0);
  sadID_train      = sID_train;
  sadID_decode     = sID_decode;
  detectOverlap    = false; //(sID_train == sID_decode);
  trainCluster     = new TrainPhoneModel*[initialNrClusters+1];  // overlapModel
  clusterScore     = new double[initialNrClusters];
  clusterSize      = new int[initialNrClusters];
  phoneModels      = ((PhoneModel**)trainCluster);
  numberOfWords    = initialNrClusters;
  numberOfClusters = initialNrClusters;
  vocabulary       = new char*[numberOfWords];
  languageModel    = new LanguageModel_Segmenter(initialNrClusters+1); // overlapModel
  numberOfPhones   = initialNrClusters;
  nrMergeIterations= 0;

  for(int i=0;i<initialNrClusters-nrNonMergeModels;i++)
  {
    clusterScore[i] = 0.0;
    char str[11];
    sprintf(str,"SPK%02d",i);
    vocabulary[i] = new char[11];
    strcpy(vocabulary[i],str);  
    vocabulary[i][strlen(str)] = 0;
    trainCluster[i] = new TrainPhoneModel(str,1,1,true, featurePool->getVectorSize(),featurePool->getStreamInfo()); 
    trainCluster[i]->doNotuseBordersForTraining(false); //detectOverlap);
  }
  FILE *noMergeFile = fopen(nonMergeableModels,"rb");

  if(noMergeFile != NULL)
  {
    int temp;
    freadEndianSafe(&temp,1,sizeof(temp),noMergeFile);
    freadEndianSafe(&temp,1,sizeof(temp),noMergeFile);
    freadEndianSafe(&temp,1,sizeof(temp),noMergeFile);

    for(int i=2;i<initialNrClusters-nrNonMergeModels;i++)
/*    for(int i=initialNrClusters-nrNonMergeModels;i<initialNrClusters;i++)*/
    {
      PRINTF2("Reading model %02d\n",i);
      char str[11];
      sprintf(str,"FIX%02d",i-initialNrClusters-nrNonMergeModels+1);
      if(trainCluster[i] != NULL)
      {
        delete trainCluster[i];
      }
      trainCluster[i] = new TrainPhoneModel(str,1,1,true, featurePool->getVectorSize(),featurePool->getStreamInfo());
      assert(trainCluster[i]->readModel(noMergeFile));
      trainCluster[i]->doNotuseBordersForTraining(false); //detectOverlap);
      trainCluster[i]->setTrainingData(featurePool,segID,i);
      clusterScore[i] = 0.0;
      strcpy(trainCluster[i]->getStatistics()->name,str);
      vocabulary[i] = new char[11];
      strcpy(vocabulary[i],str);
      vocabulary[i][strlen(str)] = 0;
    }
    fclose(noMergeFile);
  }
  trainCluster[initialNrClusters] = NULL;

  compareClusters = new int[initialNrClusters*initialNrClusters];

  createOverlapTree(CLUSTER_MIN_NUMBER_FRAMES);
  segID        = featurePool->claimSegmentationID();
  segID_decode = segID;
//  segID_decode = featurePool->claimSegmentationID();

//  overwriteWeightPars(10.0,CLUSTER_TRANS_PENALTY,CLUSTER_SIL_PENALTY);
  overwriteWeightPars(0.0,CLUSTER_TRANS_PENALTY,CLUSTER_SIL_PENALTY);

#ifdef MAX_NR_THREADS
  for(int i=0;i<MAX_NR_THREADS;i++)
  {
    cluster_threadRunning[i]   = false;
    cluster_threadMayStart[i]  = false;
    trainThreadCluster[i].id  = i;
    pthread_create( &cluster_thread[i], NULL, thread_train_cluster, (void*) &trainThreadCluster[i]);
  }
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool loadedOverlapAM = false;

void Train_Segmenter::createOverlapTree(int min)
{
  if(trainCluster[numberOfWords] == NULL)
  {
    createLexicalTree(min);
  }
  else
  {
    numberOfPhones++; // for overlap
    int forEachM[numberOfPhones];
    for(int i=0;i<numberOfWords;i++)
    {
      forEachM[i] = min;
    }
    forEachM[numberOfWords] = 10;
    createLexicalTree(min,forEachM);
    numberOfPhones--; // real nr, without overlap
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor will delete all claimed memory (the clusters)
/////////////////////////////////////////////////////////////////////////////////////////////////////

Train_Segmenter::~Train_Segmenter()
{
  featurePool->freeSegmentationID(segID);
//  featurePool->freeSegmentationID(segID_decode);

  for(int i=0;i<numberOfPhones;i++)
  {
    if(trainCluster[i] != NULL)
    {
      delete trainCluster[i];
    }
  }
  delete[] compareClusters;
  delete[] trainCluster;
  phoneModels = NULL; // pointed to trainCluster..
  delete[] clusterScore;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// We have done training, let's start merging!
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Train_Segmenter::createInitialModels(int segAmount, char *outputName)
{
  segAmount  = segAmount;   // Not used for this version...
  outputName = outputName;  // Not used for this version...

  int nrSpeech = featurePool->getClusterLength(sadID_train,0); // - featurePool->getClusterLength(segID,numberOfWords);
  double score = trainClusters(MIN_NUMBER_OF_GAUSSIANS);

  score /= ((double)nrSpeech);
#ifdef TWOSTREAMTHING
  double score2 = 1.0;
  if(score < -23.0)
  {
    score2 = 1.0;
  }
  else
  {
    score2 = 1.0 - 1.0 / ((23.0 + score)*7.0);
  }
#else
  double score2 = 10.0 + 5.0*(-24.2 - score);
  if(score2 < 5.0)
  {
    score2 = 5.0;
  }
  if(score2 > 15.0)
  {
    score2 = 15.0;
  }
#endif
//  PRINTF2("Setting Prior scaling factor (LM-scale) to %f (%f)\n",score2,score);
//  overwriteWeightPars(score2,CLUSTER_TRANS_PENALTY,CLUSTER_SIL_PENALTY);

  score = trainIteration(2);

  nrSpeech = featurePool->getClusterLength(sadID_train,0) - featurePool->getClusterLength(segID,numberOfWords);
  score /= ((double)nrSpeech);
#ifdef TWOSTREAMTHING
  if(score < -22.3)
  {
    score2 = 1.0;
  }
  else
  {
    score2 = 1.0 - 1.0 / ((22.3 + score)*7.0);
  }
#else
  score2 = 10.0 + 5.0*(-23.4 - score);
  if(score2 < 5.0)
  {
    score2 = 5.0;
  }
  if(score2 > 15.0)
  {
    score2 = 15.0;
  }
#endif
//  printf("Setting Prior scaling factor (LM-scale) to %f (%f)\n",score2,score);
//  overwriteWeightPars(score2,CLUSTER_TRANS_PENALTY,CLUSTER_SIL_PENALTY);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// We have decided on two merge candidates. Now determine if we should actually do it, or just stop.
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool Train_Segmenter::proceedMerge(int model1, int model2, int method)
{
  return false; 
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will store the trained clusters to outFile.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Train_Segmenter::storeClusters(FILE *outFile, char *outputName)
{
  assert(featurePool != NULL);
  char tempName[MAXSTR];

  int tempSeg = featurePool->claimSegmentationID();
  FILE *rttmFile = fopen(outputName,"rt");
  featurePool->readRTTM(tempSeg,NULL,-1,-1,3,rttmFile);
  fclose(rttmFile);


  numberOfClusters = featurePool->getNumberOfSpeakers(tempSeg);
  int vectorSize = featurePool->getVectorSize();
  fwrite(&vectorSize,1,sizeof(vectorSize),outFile);
  fwrite(&numberOfClusters,1,sizeof(numberOfClusters),outFile);
  int temp = 1;
  fwrite(&temp,1,sizeof(temp),outFile);

  int numberOfVectors = featurePool->getClusterLength(sadID_decode,0);

  int totL = 0;
  int num = 0;

  for(int i=0;i<featurePool->getNumberOfSpeakers(tempSeg);i++)
  {
    int id = featurePool->getSpeakerID(tempSeg,i);
    PRINTF2("Storing model %d...\n",id);
    trainCluster[id]->writeModel(outFile);
    int l = featurePool->getClusterLength(tempSeg,i);
/*
      sprintf(tempName,"%s.posteriors.SPK%02d", outputName, id);
      FILE *file = fopen(tempName,"wb");
      SegmentationAdmin readPool_curSeg;
      Vector **vector = featurePool->getVectorList(&readPool_curSeg,0);
      while(vector != NULL)
      {
        float res = trainCluster[id]->getPDFProbability(0,vector[0]);
        fwrite(&res,1,sizeof(float),file);
        vector = featurePool->getNextVectorList(&readPool_curSeg);
      }
      fclose(file);
*/
    num++;
  }
  PRINTF3("Total number of models: %d (%d)\n",num, numberOfClusters);

//  USED FOR GMM training (ICASSP 2010 paper)
//   for(int i=0;i<114;i++)
//   {
//     int id = 1;
//     printf("Storing model %d...\n",id);
//     trainCluster[id]->writeModel(outFile);
//     sprintf(tempName,"%s.posteriors.SPK%02d", outputName, i);
//       FILE *file = fopen(tempName,"wb");
//       SegmentationAdmin readPool_curSeg;
//       Vector **vector = featurePool->getVectorList(&readPool_curSeg,0);
//       while(vector != NULL)
//       {
//         float res = trainCluster[id]->getPDFProbability(0,vector[0],i);
//         fwrite(&res,1,sizeof(float),file);
//         vector = featurePool->getNextVectorList(&readPool_curSeg);
//       }
//       fclose(file);
//     num++;
//   }
//   printf("Total number of models: %d (%d)\n",num, numberOfClusters);



//   for(int i=0;i<57;i++)
//   {
//     int id = i+1;
//     printf("Storing model %d...\n",id);
//     sprintf(tempName,"%s.posteriors.SPK%02d", outputName, id);
//     FILE *file = fopen(tempName,"wb");
//     SegmentationAdmin readPool_curSeg;
//     Vector **vector = featurePool->getVectorList(&readPool_curSeg,0);
//     while(vector != NULL)
//     {
//       float res = vector[0]->getValue(i);
//       fwrite(&res,1,sizeof(float),file);
//       vector = featurePool->getNextVectorList(&readPool_curSeg);
//     }
//     fclose(file);
//   }
//   printf("Total number of models: %d (%d)\n",num, numberOfClusters);

  fclose(outFile);
  outFile = NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will load the trained clusters from inFile.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Train_Segmenter::loadClusters(FILE *inFile)
{
  for(int i=0;i<numberOfWords;i++)
  {
    if(trainCluster[i] != NULL)
    {
      delete trainCluster[i];
      trainCluster[i] = NULL;
    }
  }
  assert(featurePool != NULL);
  int vectorSize;
  freadEndianSafe(&vectorSize,1,sizeof(vectorSize),inFile);
  assert(featurePool->getVectorSize() == vectorSize);
  freadEndianSafe(&numberOfClusters,1,sizeof(numberOfClusters),inFile);
  for(int i=1;i<numberOfClusters;i++)
  {
    char str[11];
    sprintf(str,"SPK%02d",i);
    trainCluster[i] = new TrainPhoneModel(str,1,1,true, featurePool->getVectorSize(),featurePool->getStreamInfo());
    trainCluster[i]->readModel(inFile);
    trainCluster[i]->setTrainingData(featurePool,segID,i);
    trainCluster[i]->doNotuseBordersForTraining(false); //detectOverlap);
  }
  fclose(inFile);
  inFile = NULL;
  createOverlapTree(CLUSTER_MIN_NUMBER_FRAMES);
#ifndef DIARIZATION_CHEAT_EXP8
  segmentFeaturePool(sadID_train,0,segID);
#endif
  #ifdef DIARIZATION_CHEAT_NOFA
  int tempSeg = featurePool->claimSegmentationID();
  featurePool->segmentationIntersection(segID,segID_cheat,-1,tempSeg);
  featurePool->freeSegmentationID(segID);
  segID        = tempSeg;
  segID_decode = tempSeg;
  for(int i=1;i<numberOfClusters;i++)
  {
    if(trainCluster[i] != NULL)
    {
      trainCluster[i]->setTrainingData(featurePool,segID,i);
    }
  }
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Merges two single models and calculates the total score..
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Train_Segmenter::getMergeModelScore(int model1, int model2, double *mergeScore)
{
  TrainPhoneModel *merge = new TrainPhoneModel(trainCluster[model1],trainCluster[model2]);
  merge->setTrainingData(featurePool,segID,model1,model2);
  merge->train(-1,(CLUSTER_MIN_NUMBER_FRAMES != 3),true);

  int nrChannels = 1;
  if(featurePool->getStreamInfo() != NULL)
  {
    nrChannels = featurePool->getStreamInfo()->numberOfChannels; 
  }
  double weights[nrChannels];
  PRINTF3("Merge of %02d and %02d = ",model1,model2);
  for(int i=0;i<nrChannels;i++)
  {
    weights[i] = 0.0;
  }
  for(int c=0;c<nrChannels;c++)
  {
    weights[c] = 1.0;
    featurePool->setStreamWeights(weights,false);


/*
    int    n    = 0;
    double f[2] = {0.0, 0.0};
    double s[2] = {0.0, 0.0};
    bool isLast    = false;
    SegmentationAdmin readPool_curSeg;
    Vector *vector = featurePool->getFirstVectorFirstSegment(&readPool_curSeg,segID,model1,false,&isLast);
    while(vector != NULL)
    {
      n++;
      double ubm = merge->getLogPDFProbability(0,vector);
      double s1  = (trainCluster[model1]->getLogPDFProbability(0,vector) - ubm);
      double s2  = (trainCluster[model2]->getLogPDFProbability(0,vector) - ubm);
      f[0]+= s1;
      f[1]+= s2;
      s[0]+= (s1*s1);
      s[1]+= (s2*s2);
      while(!isLast)
      {
        vector = featurePool->getNextVector(&readPool_curSeg,&isLast);
        n++;
        double ubm = merge->getLogPDFProbability(0,vector);
        double s1  = (trainCluster[model1]->getLogPDFProbability(0,vector) - ubm);
        double s2  = (trainCluster[model2]->getLogPDFProbability(0,vector) - ubm);
        f[0]+= s1;
        f[1]+= s2;
        s[0]+= (s1*s1);
        s[1]+= (s2*s2);
      }
      vector = featurePool->getFirstVectorNextSegment(&readPool_curSeg,false,&isLast);
    }

    isLast    = false;
    vector = featurePool->getFirstVectorFirstSegment(&readPool_curSeg,segID,model2,false,&isLast);
    while(vector != NULL)
    {
      n++;
      double ubm = merge->getLogPDFProbability(0,vector);
      double s1  = (trainCluster[model2]->getLogPDFProbability(0,vector) - ubm);
      double s2  = (trainCluster[model1]->getLogPDFProbability(0,vector) - ubm);
      f[0]+= s1;
      f[1]+= s2;
      s[0]+= (s1*s1);
      s[1]+= (s2*s2);
      while(!isLast)
      {
        vector = featurePool->getNextVector(&readPool_curSeg,&isLast);
        n++;
        double ubm = merge->getLogPDFProbability(0,vector);
        double s1  = (trainCluster[model2]->getLogPDFProbability(0,vector) - ubm);
        double s2  = (trainCluster[model1]->getLogPDFProbability(0,vector) - ubm);
        f[0]+= s1;
        f[1]+= s2;
        s[0]+= (s1*s1);
        s[1]+= (s2*s2);
      }
      vector = featurePool->getFirstVectorNextSegment(&readPool_curSeg,false,&isLast);
    }
    double m[2] = {0.0, 0.0};
    double v[2] = {0.0, 0.0};
    m[0] = f[0] / n;
    m[1] = f[1] / n;
    v[0] = s[0] / n - m[0]*m[0];
    v[1] = s[1] / n - m[1]*m[1];
    double mean = m[0] - m[1];
    if(mean < 0.0)
    {
      mean = m[1] - m[0];
    }
    mergeScore[c] = 100.0 * sqrt(v[0]+v[1]) / (mean * sqrt(n));
*/

    mergeScore[c] = merge->getTrainSilP(model1) + merge->getTrainSilP(model2) - 
                    trainCluster[model1]->getTrainSilP(model1) - trainCluster[model2]->getTrainSilP(model2);
    PRINTF2("%f\t",mergeScore[c]);
    weights[c] = 0.0;
  }
  PRINTF("\n");
  featurePool->resetStreamWeights();

  delete merge;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Merges two single models (over loaded by Adapt_Segmenter)
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Train_Segmenter::mergeModels(int model1, int model2)
{
#ifdef DIARIZATION_CHEAT_EXP5
  char unWidened[MAXSTR];
  sprintf(unWidened,"%s.%d.mergeChoise",outputFileName,globalIt-1);
  FILE *textClFile = fopen(unWidened,"wt");
  fprintf(textClFile,"SPK%02d\tSPK%02d\n",model1,model2);
  fclose(textClFile);
  fflush(stdout);
#if !(defined DIARIZATION_CHEAT_EXP3) && !(defined DIARIZATION_CHEAT_EXP4)
  char fileName[2000];
  sprintf(fileName,"/home/huijbreg/nist/diarytools/cheatPickMerge.sh %s.%d %s",outputFileName,globalIt-1,label);
  int resNotUsed = system(fileName);
#endif
#endif

  PRINTF3("Merging SPK%02d and SPK%02d..\n",model1,model2);
  TrainPhoneModel *merge = new TrainPhoneModel(trainCluster[model1],trainCluster[model2]);
  merge->setTrainingData(featurePool,segID,model1,model2);

  merge->train(-1,(CLUSTER_MIN_NUMBER_FRAMES != 3),true);  // -1 for maxgaussians: do not split any more...

  delete trainCluster[model1];
  delete trainCluster[model2];
  trainCluster[model1] = merge;
  trainCluster[model2] = NULL;
#ifdef DIARIZATION_CHEAT_EXP8
 featurePool->mergeLabels(segID,model1,model2);
#endif
#ifdef DIARIZATION_CHEAT_NOFA
  int tempSeg = featurePool->claimSegmentationID();
  featurePool->segmentationIntersection(segID,segID_cheat,-1,tempSeg);
  featurePool->freeSegmentationID(segID);
  segID        = tempSeg;
  segID_decode = tempSeg;
  for(int i=1;i<numberOfClusters;i++)
  {
    if(trainCluster[i] != NULL)
    {
      trainCluster[i]->setTrainingData(featurePool,segID,i);
    }
  }

#endif

  PRINTF("Done merging!\n");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Trains a single model with nrG gaussians (over loaded by Adapt_Segmenter)
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Train_Segmenter::trainModel(int model , int nrG, double *trainPRes)
{
  *trainPRes = -9.0e300;
  if((trainCluster[model] != NULL) && mergeable[model])
  {
//  USED FOR GMM training (ICASSP 2010 paper)
/*     trainCluster[model]->setTrainingData(featurePool,sadID_train,0);
     if(!featurePool->doSegmentsExist(sadID_train,0))*/
    trainCluster[model]->setTrainingData(featurePool,segID,model);
    if(!featurePool->doSegmentsExist(segID,model))
    {
      PRINTF2("Wow, this model (%d) does not contain any data! Let's delete..\n", model);
      delete trainCluster[model];
      trainCluster[model] = NULL;
      createOverlapTree(CLUSTER_MIN_NUMBER_FRAMES);
    }
    else
    {
#ifndef MAX_NR_THREADS
      if(discrTrainMask != NULL)
      {
        trainCluster[model]->train(nrG,(CLUSTER_MIN_NUMBER_FRAMES != 3),true,discrTrain,discrTrainMask[model]);
      }
      else
      {
//  USED FOR GMM training (ICASSP 2010 paper)
//        trainCluster[model]->train(114,(CLUSTER_MIN_NUMBER_FRAMES != 3),true,NULL,NULL);
        trainCluster[model]->train(nrG,(CLUSTER_MIN_NUMBER_FRAMES != 3),true,NULL,NULL);
      }
      *trainPRes = trainCluster[model]->getTrainSilP(model);

#else
      int freeID = -1;
      // First find a free thread:
      pthread_mutex_lock( &condition_cluster_mutexDone );
      while(freeID < 0)
      {
        int i = 0;
        while(cluster_threadRunning[i] && i < MAX_NR_THREADS)
        {
          i++;
        }
        if(i< MAX_NR_THREADS)
        {
          freeID = i;
        }
        if(freeID < 0)
        {
          pthread_cond_wait( &condition_cluster_threadDone, &condition_cluster_mutexDone );
        }
      }
      assert(freeID >= 0 && freeID < MAX_NR_THREADS);
      cluster_threadRunning[freeID] = true;
      pthread_mutex_unlock( &condition_cluster_mutexDone );

      trainThreadCluster[freeID].model = trainCluster[model];
      trainThreadCluster[freeID].nr    = nrG;
      trainThreadCluster[freeID].segid = model;
      trainThreadCluster[freeID].trainPRes = trainPRes;

      pthread_mutex_lock( &condition_cluster_mutexStart[freeID] );
      cluster_threadMayStart[freeID] = true;
      pthread_cond_signal( &condition_cluster_threadStart[freeID] );
      pthread_mutex_unlock( &condition_cluster_mutexStart[freeID] );
#endif
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The training iterations are performed by this method.
/////////////////////////////////////////////////////////////////////////////////////////////////////

double Train_Segmenter::trainClusters(int nrG)
{
  double res = 0.0;
  for(int i=1;i<numberOfWords;i++)
  {
    if((trainCluster[i] != NULL) && mergeable[i])
    {
      int n = nrG;
      if(n <= 0)
      {
        n = trainCluster[i]->maxNrOfGaussians();
      }
      trainModel(i,n,&clusterScore[i]);
    }
  }
#ifdef MAX_NR_THREADS
  // Wait until all threads are done.
  pthread_mutex_lock( &condition_cluster_mutexDone );
  bool waitForMe = true;
  while(waitForMe)
  {
    int i = 0;
    while(!cluster_threadRunning[i] && i < MAX_NR_THREADS)
    {
      i++;
    }
    waitForMe = (i < MAX_NR_THREADS);
    if(waitForMe)
    {
      pthread_cond_wait( &condition_cluster_threadDone, &condition_cluster_mutexDone );
    }
  }
  pthread_mutex_unlock( &condition_cluster_mutexDone );
#endif

  for(int i=1;i<numberOfWords;i++)
  {
    if((trainCluster[i] != NULL) && mergeable[i])
    {
      if(clusterScore[i] > -9.0e100)
      {
        res += clusterScore[i];
      }
    }
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The main method of a train iteration:
/// - align the data using Segmenter::segmentFeaturePool()
/// - train the clusters (states) with trainClusters()
/////////////////////////////////////////////////////////////////////////////////////////////////////

double Train_Segmenter::trainIteration(int nrIterations)
{
  double score = 0.0;
  double oldScore = -9.0e100;
  int count = 0;

  int countIterations = 0;
// || ((score - oldScore) / absolute(oldScore) > MINIMUM_TRAIN_IMPROVEMENT_CLUSTER))
  while((count < nrIterations))
  {
    countIterations++;
    oldScore = score;
/*
    if(count >= nrIterations)
    {
      count = 0;
      oldScore = score;
    }
*/    
    // Make a nice alignment using segmentation 0, which is speech/non-speech (segment 0 which is SPEECH!):
    // Store it at segmentID 1...
/*    
#ifdef SAD_TEST        
  int forEachModel1[4] = {20,50,50,500};
  double forEachModel2[4] = {-1.0,-1.0,-50.0,-500.0};
  createLexicalTree(50,forEachModel1, forEachModel2);
#endif        
  */  
    
    PRINTF("Re-segment...\n"); fflush(stdout);
    createOverlapTree(CLUSTER_MIN_NUMBER_FRAMES);
#ifndef DIARIZATION_CHEAT_EXP8
    PRINTF("Segment now...\n"); fflush(stdout);
    segmentFeaturePool(sadID_train,0,segID);
#endif
#ifdef DIARIZATION_CHEAT_NOFA
  int tempSeg = featurePool->claimSegmentationID();
  featurePool->segmentationIntersection(segID,segID_cheat,-1,tempSeg);
  featurePool->freeSegmentationID(segID);
  segID        = tempSeg;
  segID_decode = tempSeg;
  for(int i=1;i<numberOfClusters;i++)
  {
    if(trainCluster[i] != NULL)
    {
      trainCluster[i]->setTrainingData(featurePool,segID,i);
    }
  }

#endif

    LanguageModel_Segmenter *lms = (LanguageModel_Segmenter *)languageModel;
    SegmentationAdmin readPool_curSeg;
    Vector *vector = featurePool->getFirstVectorFirstSegment(&readPool_curSeg,segID,-1,true);
    while(vector != NULL)
    {
      int word = featurePool->getSegmentID(0,&readPool_curSeg);
      lms->transferTo(word);
      vector = featurePool->getFirstVectorNextSegment(&readPool_curSeg,true);
    }
    lms->finishModel();
    lms->addPenalty(numberOfWords,0.3);

    PRINTF("Done segmenting, let's train\n");fflush(stdout);
    score = trainClusters();
    count++;
  }
  PRINTF2("Resegmented %d times.\n",count);
  return score;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The main method of this class. It will 'train' new clusters with maximal maxClusters. If 
/// maxClusters is equal to zero, regular speaker clustering is performed. If not, the following
/// steps are added:
/// - Merge small clusters (clusters with not enough data used during training) to their closest 
/// clusters until all clusters contain are big enough.
/// - Merge the remaining clusters until no more than maxClusters clusters are left.
/////////////////////////////////////////////////////////////////////////////////////////////////////


FeaturePool *pass2Pool   = NULL;
int pass2PoolindexWeight = 0;
int pass2PoolindexShift  = 0;

void Train_Segmenter::train(int maxClusters, char *l, char *outputName, char *feaPosteriors, int segID_cht)
{
  segID_cheat = segID_cht;
  SegmentationAdmin readPool_curSeg;
  double score = -9.0e300;   
  int milliSec;
  struct timeval time1;
  struct timeval time2;
  gettimeofday( &time1, NULL );
  
  
PRINTF2("Number of CHANNELS: %d\n",featurePool->getStreamInfo()->numberOfChannels);
for(int i=0;i<featurePool->getStreamInfo()->numberOfChannels;i++)
{
  PRINTF2("Weights: %f\n",featurePool->getStreamInfo()->weight[i]);
}


  strcpy(label,l);
  
  strcpy(outputFileName,outputName);
  
  if(maxClusters > 0)
  {
    maxClusters++; // for administration reasons...
  }

//  FILE *of = fopen("/kletsmajoor/projects/rt09s/overlap/alloverlap.am","rb");

//   char filename[1000];
//   strcpy(filename,outputName);
//   strcat(filename,"overlap.am.in");
//   FILE *of = fopen(filename,"rb");
//   if(of != NULL)
//   {
//     loadedOverlapAM             = true;
//     trainCluster[numberOfWords] = new TrainPhoneModel(of,featurePool->getVectorSize(),featurePool->getStreamInfo());
//     fclose(of);
//   }
//   strcpy(filename,outputName);
//   strcat(filename,".pass2");
//   pass2Pool = new FeaturePool(filename,FILE_TYPE_TEXT);


  featurePool->resetSegmentation(segID);
  int numberOfVectors = featurePool->getClusterLength(sadID_train,0);
  int segAmount = numberOfVectors / (numberOfClusters-1) / SEGMENT_INITIAL_DATA_SEGMENTS;
  PRINTF2("Total number of speech frames: %d\n",numberOfVectors);


  PRINTF2("Number of GAUSSIANS: %d\n",MIN_NUMBER_OF_GAUSSIANS);


/////////////////////////////////////
// INITIALIZATION STEP.
/////////////////////////////////////

#ifndef STARTWITHLAST

#if !(defined DIARIZATION_CHEAT_EXP1) && !(defined DIARIZATION_CHEAT_EXP3) && !(defined DIARIZATION_CHEAT_EXP7) && !(defined DIARIZATION_CHEAT_EXP8)

    /////////////////// Initially divide data:
    initialiseSystem();

#ifndef ICASSP2010GMM

    int i = 1;
    int count = 0;
    Vector *vector = featurePool->getFirstVectorFirstSegment(&readPool_curSeg,sadID_train,0,false);
    while(vector != NULL)
    {
      int length = featurePool->getCurSegmentLen(&readPool_curSeg);
      int tot    = featurePool->getCurSegmentStart(&readPool_curSeg);
      int ind = 0;
      int indPrev = 0;
      while(ind+1 < length)
      {
        if((length - ind) < (segAmount - count))
        {
          // Easy does it, use the entire sentence:
          int lframe = tot+length-ind-1;
          if(lframe - tot > 0)
          {
            featurePool->addSegment(segID,tot,lframe,i);
          }
          tot   += (length - ind);
          count += (length - ind);
          ind    = length;
        } 
        else
        {
          int lframe = tot+(segAmount - count)-1;
          if(lframe - tot > 0)
          {
            featurePool->addSegment(segID,tot,lframe,i);
          }
          tot  += (segAmount - count);
          ind  += (segAmount - count);
          count = segAmount;
        }
        indPrev = ind;
        if(count >= segAmount)
        {
          count = 0;
          if(i+1 < numberOfClusters)
          {
            i++;
          }
          else
          {
            i = 1;
            if((segAmount*SEGMENT_INITIAL_DATA_SEGMENTS)  < (numberOfVectors / (numberOfClusters-1)))
            {
              segAmount++;
            }
          }
        }
      }
      vector = featurePool->getFirstVectorNextSegment(&readPool_curSeg,false);
    }

    int    nrChannels = 1;
    if(featurePool->getStreamInfo() != NULL)
    {
      nrChannels = featurePool->getStreamInfo()->numberOfChannels;
    }
    createInitialModels(numberOfVectors,outputName);
#endif
#ifdef ICASSP2010GMM
/*  for(int i=0;i<numberOfClusters-nrNonMergeModels;i++)
  {
    if(trainCluster[i] != NULL && trcl != NULL && i>0)
    {
      printf("Setting gaussian for model %d...\n",i);
      trainCluster[i]->setGaussian(0,0,trcl->getGaussian(0,i-1));
    }
  }*/
  featurePool->resetSegmentation(segID);
  createLexicalTree(CLUSTER_MIN_NUMBER_FRAMES);
  segmentFeaturePool(sadID_train,0,segID);
#endif

#endif
#ifdef DIARIZATION_CHEAT_EXP1
//    Determine amount of data per speaker:
    int spkrData[numberOfWords];
    int totData = 0;
    for(int i=1;i<numberOfWords;i++)
    {
      Vector *vector = featurePool->getFirstVectorFirstSegment(&readPool_curSeg,segID_cheat,i,false);
      spkrData[i] = 0;
      while(vector != NULL)
      {
        spkrData[i] += featurePool->getCurSegmentLen(&readPool_curSeg);
        vector = featurePool->getFirstVectorNextSegment(&readPool_curSeg,false);
      }
      totData += spkrData[i];
      PRINTF3("Data for speaker %02d: %d\n",i,spkrData[i]);
    }    
    PRINTF2("Data total: %d\n",totData);
//    numberOfClusters = 17;
    double gaussianPerVector = ((double)((numberOfClusters-1)*MIN_NUMBER_OF_GAUSSIANS)) / ((double)totData);
    for(int i=1;i<numberOfWords;i++)
    {
      int nrG = (int)(((double)(spkrData[i]))*gaussianPerVector) + 1;
      PRINTF3("Gaussians for speaker %02d: %d\n",i,nrG);
    }

    int backupSeg = segID;
    segID         = segID_cheat;
    initialiseSystem();
    for(int i=1;i<numberOfWords;i++)
    {
      if(spkrData[i] > 0)
      {
        int nrG = (int)(((double)(spkrData[i]))*gaussianPerVector) + 1;
        if(nrG < MIN_NUMBER_OF_GAUSSIANS)
        {
          nrG = MIN_NUMBER_OF_GAUSSIANS;
        }

        double res = -9.0e30;
        trainModel(i,nrG,&res);
      }
      else
      {
        delete trainCluster[i];
        trainCluster[i] = NULL;
      }
    }
#ifdef MAX_NR_THREADS
    // Wait until all threads are done.
    pthread_mutex_lock( &condition_cluster_mutexDone );
    bool waitForMe = true;
    while(waitForMe)
    {
      int i = 0;
      while(!cluster_threadRunning[i] && i < MAX_NR_THREADS)
      {
        i++;
      }
      waitForMe = (i < MAX_NR_THREADS);
      if(waitForMe)
      {
        pthread_cond_wait( &condition_cluster_threadDone, &condition_cluster_mutexDone );
      }
    }
    pthread_mutex_unlock( &condition_cluster_mutexDone );
#endif

    segID         = backupSeg;
    createLexicalTree(CLUSTER_MIN_NUMBER_FRAMES);
//    trainIteration(1);

#endif
#if (defined DIARIZATION_CHEAT_EXP3) || (defined DIARIZATION_CHEAT_EXP7) || (defined DIARIZATION_CHEAT_EXP8)

    // Determine number of speakers in the reference file:
    int nrSpeakers = 0;
    Vector *vector = featurePool->getFirstVectorFirstSegment(&readPool_curSeg,segID_cheat,nrSpeakers+1,false);
    while(vector != NULL)
    {
      nrSpeakers++;
      vector = featurePool->getFirstVectorFirstSegment(&readPool_curSeg,segID_cheat,nrSpeakers+1,false);
    }
    // Determine amount of data per speaker:
    int spkrData[nrSpeakers];
    int totData = 0;
    for(int i=0;i<nrSpeakers;i++)
    {
      vector = featurePool->getFirstVectorFirstSegment(&readPool_curSeg,segID_cheat,i+1,false);
      spkrData[i] = 0;
      while(vector != NULL)
      {
        spkrData[i] += featurePool->getCurSegmentLen(&readPool_curSeg);       
        vector = featurePool->getFirstVectorNextSegment(&readPool_curSeg,false);
      }
      totData += spkrData[i];
      PRINTF3("Data for speaker %02d: %d\n",i+1,spkrData[i]);
    }    
    PRINTF2("Data total: %d\n",totData);
    int datPerSpeaker = totData / (numberOfClusters-1);
    int spkrPerRef[nrSpeakers];
    int tot = 0;
    PRINTF2("Data per speaker: %d\n",datPerSpeaker);
    for(int i=0;i<nrSpeakers;i++)
    {
      spkrPerRef[i] = spkrData[i] / datPerSpeaker;
      if(spkrPerRef[i] < 1)
      {
        spkrPerRef[i] = 1;
      }
      tot+=spkrPerRef[i];
    }
    if((numberOfClusters-1)-tot > 0) // Missing models
    {
      for(int a=0;a<(numberOfClusters-1)-tot;a++)
      {
        double minError = 9.0e300;
        int    minI     = 0;
        for(int i=0;i<nrSpeakers;i++)
        {
          double orgError = pow(datPerSpeaker - (spkrData[i]/spkrPerRef[i]),2)*spkrPerRef[i];
          double newError = pow(datPerSpeaker - (spkrData[i]/(spkrPerRef[i]+1)),2)*spkrPerRef[i];
          if(newError - orgError < minError)
          {
            minError = newError - orgError;
            minI     = i;
          }
        }
        spkrPerRef[minI]++;          
      }    
    }
    else if((numberOfClusters-1)-tot < 0) // too many models
    {
      for(int a=0;a<tot-(numberOfClusters-1);a++)
      {
        double minError = 9.0e300;
        int    minI     = 0;
        for(int i=0;i<nrSpeakers;i++)
        {
          if(spkrPerRef[i] > 1)
          {
            double orgError = pow(datPerSpeaker - (spkrData[i]/spkrPerRef[i]),2)*spkrPerRef[i];
            double newError = pow(datPerSpeaker - (spkrData[i]/(spkrPerRef[i]-1)),2)*spkrPerRef[i];
            if(newError - orgError < minError)
            {
              minError = newError - orgError;
              minI     = i;
            }
          }
        }
        spkrPerRef[minI]--;          
      }        
    }
    for(int i=0;i<nrSpeakers;i++)
    {
      PRINTF3("Number of speakers for %02d: %d\n",i+1,spkrPerRef[i]);
    }
    int spkrID     = 1;    
    for(int i=0;i<nrSpeakers;i++)
    {
      int size = 0;
      int sizePerSpkr = spkrData[i]/spkrPerRef[i];
      int nextSpkrID = spkrID + spkrPerRef[i] - 1;
      vector = featurePool->getFirstVectorFirstSegment(&readPool_curSeg,segID_cheat,i+1,false);
      while(vector != NULL)
      {
        int len = featurePool->getCurSegmentLen(&readPool_curSeg);
        if(size + len < sizePerSpkr)
        {
          size+=len;
          featurePool->addSegment(segID,featurePool->getCurSegmentStart(&readPool_curSeg),featurePool->getCurSegmentStart(&readPool_curSeg)+len-1,spkrID);
        }        
        else
        {
          int part = sizePerSpkr - size;
          featurePool->addSegment(segID,featurePool->getCurSegmentStart(&readPool_curSeg),featurePool->getCurSegmentStart(&readPool_curSeg)+part-1,spkrID);
          size = len - part;
          if(spkrID < nextSpkrID)
          {
            spkrID++;
            assert(size >= 0);
            if(size > 0)
            {
              featurePool->addSegment(segID,featurePool->getCurSegmentStart(&readPool_curSeg)+part,featurePool->getCurSegmentStart(&readPool_curSeg)+len-1,spkrID);          
            }          
          }
        }
        vector = featurePool->getFirstVectorNextSegment(&readPool_curSeg,false);
      }
      spkrID++;
    }    
    PRINTF2("SpeakerID: %d\n",spkrID-1);
    for(int i=0;i<spkrID-1;i++)
    {
      vector = featurePool->getFirstVectorFirstSegment(&readPool_curSeg,segID,i+1,false);
      int dat = 0;
      while(vector != NULL)
      {
        dat += featurePool->getCurSegmentLen(&readPool_curSeg);       
        vector = featurePool->getFirstVectorNextSegment(&readPool_curSeg,false);
      }
      PRINTF3("Data for speaker %02d: %d\n",i+1,dat);
    }
    createInitialModels(numberOfVectors,outputName);
#endif

/////////////////////////////////////
// TRAINING STEP.
/////////////////////////////////////


#ifndef DIARIZATION_CHEAT_EXP1

  PRINTF2("Number of clusters (real): %d\n",numberOfClusters-1);  
// Initialise models:
  if(feaPosteriors != NULL)
  {
    PRINTF("Going to train models for creation of posterior-features.\n");
    trainIteration(MINPASSES_ALIGN);
    writePosteriors(feaPosteriors);
    return;
  }
  trainIteration(MINPASSES_ALIGN);

/*
  printf("MARIJN: TRAIN TRIPPLE GMM\n");fflush(stdout);

  CLUSTER_MIN_NUMBER_FRAMES = 3;
  for(int i=0;i<numberOfWords;i++)
  {
    if(trainCluster[i] != NULL)
    {
      delete trainCluster[i];
      clusterScore[i] = 0.0;
      char str[11];
      sprintf(str,"TRI%02d",i);
      trainCluster[i] = new TrainPhoneModel(str,1,1,false, featurePool->getVectorSize(),featurePool->getStreamInfo());
      trainCluster[i]->doNotuseBordersForTraining(false); //detectOverlap);
    }
  }
printf("MARIJN: TRAIN 1\n");fflush(stdout);
  createLexicalTree(CLUSTER_MIN_NUMBER_FRAMES);
printf("MARIJN: TRAIN 2\n");fflush(stdout);
  createInitialModels(numberOfVectors,outputName);
printf("MARIJN: TRAIN 3\n");fflush(stdout);
  trainIteration(MINPASSES_ALIGN);
*/



#endif // notdef DIARIZATION_CHEAT_EXP1


#if (defined DIARIZATION_CHEAT_EXP3) || (defined DIARIZATION_CHEAT_EXP4) || (defined DIARIZATION_CHEAT_EXP5) || (defined DIARIZATION_CHEAT_EXP7) || (defined GLOBAL_STOP)
  char unWidened[MAXSTR];

  int tempSeg = featurePool->claimSegmentationID();
  featurePool->resetSegmentation(tempSeg);
  createLexicalTree(CLUSTER_MIN_NUMBER_FRAMES);
//  createLexicalTree(100);
  segmentFeaturePool(sadID_decode,0,tempSeg);
  sprintf(unWidened,"%s.%d",outputName,globalIt);
  FILE *textClFile = fopen(unWidened,"wt");
  featurePool->writeRTTM(tempSeg,label,"SPK",textClFile,0,1);
  fclose(textClFile);
  PRINTF("Writing rttm to file:\n");  fflush(stdout);
  PRINTF2("%s\n",unWidened); fflush(stdout);
  globalIt++;
  fflush(stdout);
  featurePool->freeSegmentationID(tempSeg);
  createLexicalTree(CLUSTER_MIN_NUMBER_FRAMES);
#endif

/*  char str[MAXSTR];
  gettimeofday( &time1, NULL );
  sprintf(str,"%s.models.%d",outputName,globalIt);
  printf("Writing model to file: %s\n",str); fflush(stdout);
  FILE *modelFile = fopen(str,"wb");
  storeClusters(modelFile);*/
#else
  char str[MAXSTR];
  gettimeofday( &time1, NULL );
  sprintf(str,"%s.models.%d",outputName,segID_cheat);
  PRINTF2("Starting from iteration model-file: %s\n",str); fflush(stdout);
  FILE *modelFile = fopen(str,"rb");
  loadClusters(modelFile);
  PRINTF("Finished loading clusters...\n"); fflush(stdout);

  mergeClusters(0,false,false);

#endif //STARTWITHLAST



/////////////////////////////////////
// MERGING STEP.
/////////////////////////////////////



#ifndef STARTWITHLAST
#ifndef DIARIZATION_CHEAT_EXP1

//  getOverlap(NULL);

  // Regular clustering:
//  USED FOR GMM training (ICASSP 2010 paper)
//  while(false)
  while(mergeClusters(nrNonMergeModels+maxClusters, false,true))
  {
    score = trainIteration(MINPASSES_ALIGN);
    PRINTF2("**************** SCORE: %g\n", score);

#ifdef GLOBAL_STOP
    char unWidened[MAXSTR];
    sprintf(unWidened,"%s.%d.score",outputFileName,globalIt);
    FILE *tf = fopen(unWidened,"wt");
    fprintf(tf,"%d\t%f\n",globalIt,score);
    fclose(tf);
    fflush(stdout);
#endif

#if (defined DIARIZATION_CHEAT_EXP3) || (defined DIARIZATION_CHEAT_EXP4) || (defined DIARIZATION_CHEAT_EXP5) || (defined DIARIZATION_CHEAT_EXP7) || (defined GLOBAL_STOP)
    int tempSeg = featurePool->claimSegmentationID();
    featurePool->resetSegmentation(tempSeg);
    createLexicalTree(CLUSTER_MIN_NUMBER_FRAMES);
//    createLexicalTree(100);
    segmentFeaturePool(sadID_decode,0,tempSeg);
    sprintf(unWidened,"%s.%d",outputName,globalIt);
    FILE *textClFile = fopen(unWidened,"wt");
    featurePool->writeRTTM(tempSeg,label,"SPK",textClFile,0,1);
    fclose(textClFile);
    PRINTF("Writing rttm to file:\n");  fflush(stdout);
    PRINTF2("%s\n",unWidened); fflush(stdout);
    globalIt++;
    fflush(stdout);
    featurePool->freeSegmentationID(tempSeg);
    createLexicalTree(CLUSTER_MIN_NUMBER_FRAMES);
#endif


    gettimeofday( &time2, NULL );
    milliSec = (time2.tv_sec*1000 + time2.tv_usec/1000) - (time1.tv_sec*1000 + time1.tv_usec/1000);
    PRINTF2("** USED TIME          : %g\n",((float)milliSec)/1000);
    fflush(stdout);
    gettimeofday( &time1, NULL );
  }

#endif //EXP1
#endif //STARTWITHLAST

  if(detectOverlap)
  {
    char overlapName[MAXSTR];
    sprintf(overlapName,"%s.un-widened.no-overlap",outputName);
    FILE *textFile = fopen(overlapName,"wt");
    getOverlap(textFile);
    fclose(textFile);
  }
  else if(loadedOverlapAM)
  {
    char overlapName[MAXSTR];
    sprintf(overlapName,"%s.un-widened.no-overlap",outputName);
    FILE *textFile = fopen(overlapName,"wt");
    featurePool->resetSegmentation(segID);
    createLexicalTree(100);
    segmentFeaturePool(sadID_decode,0,segID);
    featurePool->writeRTTM(segID,label,"SPK",textFile,0,1);
    fclose(textFile);

    featurePool->resetSegmentation(segID);
    createOverlapTree(100);
    PRINTF("Segmenting feature pool (with overlap)..\n");fflush(stdout);
    segmentFeaturePool(sadID_decode,0,segID);
    PRINTF2("Done segmenting with overlap model (%d)...\n", numberOfWords);fflush(stdout);
    featurePool->segmentationCreateOverlap(segID,numberOfWords);
  }
  else
  {
    featurePool->resetSegmentation(segID);
    createLexicalTree(CLUSTER_MIN_NUMBER_FRAMES);
    segmentFeaturePool(sadID_decode,0,segID_decode);
  }

  if(prepareForASR)
  {
    if(outputName != NULL)
    {
      char unWidened[MAXSTR];
      sprintf(unWidened,"%s.un-widened",outputName);
      FILE *textClFile = fopen(unWidened,"wt");
      featurePool->writeRTTM(segID_decode,label,"SPK",textClFile,0,1);
      fclose(textClFile);
    }
    featurePool->segmentationWiden(segID,0,50);
  }

  if(outputName != NULL)
  {
//     double streamScore[nrChannels*numberOfWords*numberOfWords];
//     for(int i=0;i<numberOfWords;i++)
//     {
//       for(int i2=0;i2<numberOfWords;i2++)
//       {
//         for(int a=0;a<nrChannels;a++)
//         {
//           streamScore[(i*numberOfWords+i2)*nrChannels+a] = -42.42e30;
//         }
//       }
//     }
// 
//     printf("=================================================================\n");
//     for(int i=1;i<numberOfWords;i++)
//     {
//       if((trainCluster[i] != NULL))
//       {
//         for(int i2=i+1;i2<numberOfWords;i2++)
//         {
//           if((trainCluster[i2] != NULL))
//           {
//             getMergeModelScore(i,i2,&(streamScore[(i*numberOfWords+i2)*nrChannels]));
//             printf("FINAL COMPARE:\tSPK%02d\tSPK%02d\t%f\n",i,i2,streamScore[(i*numberOfWords+i2)*nrChannels]);
//           }
//         }
//       }
//     }

    FILE *textClFile = fopen(outputName,"wt");
    featurePool->writeRTTM(segID_decode,label,"SPK",textClFile,0,1);
    fclose(textClFile);

/*
    char unWidened[MAXSTR];
    sprintf(unWidened,"%s.AM",outputName);
    FILE *AMFile = fopen(unWidened,"wb");
    storeClusters(AMFile, outputName);
*/
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Train_Segmenter::writePosteriors(char *fileName)
{
    
  FILE *outFile = fopen(fileName,"wt");
  int len = featurePool->getPoolLength();
  double likelihood[numberOfWords];
  double sum;
  for(int i=0;i<len;i++)
  {
    Vector *vector = featurePool->getVector(NULL,i);
    sum = 0.0;
    for(int i2=1;i2<numberOfWords;i2++)
    {
      if(trainCluster[i2] != NULL)
      { 
        likelihood[i2] = trainCluster[i2]->getPDFProbability(0,vector,0,i)/30.0;
        sum+=likelihood[i2];
      }
    }
    for(int i2=1;i2<numberOfWords;i2++)
    {
      if(trainCluster[i2] != NULL)
      { 
        likelihood[i2] /= sum;
        fprintf(outFile,"%f ",likelihood[i2]);
      }
    }
    fprintf(outFile,"\n");
  }
  fclose(outFile);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Train_Segmenter::getOverlap(FILE *file)
{
  int length = 50;
  int overlapModelID  = numberOfWords;
  SegmentationAdmin readPool_curSeg;
  bool vectorNotNull  = true;
  bool isLast         = false;
  Vector *trainVector = featurePool->getFirstVectorFirstSegment(&readPool_curSeg,segID,-1,true, &isLast);
  vectorNotNull = (trainVector != NULL);
  for(int i=1;i<numberOfWords;i++)
  {
    if(trainCluster[i] != NULL)
    {
      trainCluster[i]->doNotuseBordersForTraining(false);
      clusterScore[i] = trainCluster[i]->getTrainSilP(i);
      clusterSize[i]  = featurePool->getClusterLength(segID,i);
    }
  }

  char str[11];
  sprintf(str,"OVERLAP");
  vocabulary[overlapModelID] = new char[11];
  strcpy(vocabulary[overlapModelID],str);
  vocabulary[overlapModelID][strlen(str)] = 0;
  trainCluster[overlapModelID] = new TrainPhoneModel(str,1,1,true, featurePool->getVectorSize(),featurePool->getStreamInfo());

  int overlapSeg = featurePool->claimSegmentationID();

  while(vectorNotNull)
  {
    int observationLength = featurePool->getCurSegmentLen(&readPool_curSeg);
    int time              = featurePool->getCurSegmentStart(&readPool_curSeg);
    int contextKey        = featurePool->getSegmentID(TRAIN_CONTEXTKEY,&readPool_curSeg);
    int stateNr           = 0;
    double likelihood[2]  = {0.0, 0.0};
    int    count[2]       = {0, 0};

    int cnt   = 0;
    int mxIn  = length;
    int mxOut = length;
    if(readPool_curSeg.prevSeg != NULL && (time - readPool_curSeg.prevSeg->lastFrame > length || readPool_curSeg.prevSeg->ID[0] == readPool_curSeg.curSeg->ID[0]))
    {
      mxIn = 0;
    }
    if(readPool_curSeg.curSeg->next != NULL && 
       (readPool_curSeg.curSeg->next->firstFrame - (time + observationLength) > length || readPool_curSeg.curSeg->next->ID[0] == readPool_curSeg.curSeg->ID[0]))
    {
      mxOut = 0;
    }

    if(observationLength > mxIn + mxOut)
    {

      if(cnt < mxIn)
      {
        likelihood[0] += trainCluster[readPool_curSeg.curSeg->ID[0]]->getLogPDFProbability(0,trainVector);
        count[0]++;
      }
      else if(observationLength - cnt < mxOut)
      {
        likelihood[1] += trainCluster[readPool_curSeg.curSeg->ID[0]]->getLogPDFProbability(0,trainVector);
        count[1]++;
      }
      cnt++;

      while(!isLast)
      {
        trainVector = featurePool->getNextVector(&readPool_curSeg,&isLast);
        time++;
        if(cnt < mxIn)
        {
          likelihood[0] += trainCluster[readPool_curSeg.curSeg->ID[0]]->getLogPDFProbability(0,trainVector);
          count[0]++;
        }
        else if(observationLength - cnt < mxOut)
        {
          likelihood[1] += trainCluster[readPool_curSeg.curSeg->ID[0]]->getLogPDFProbability(0,trainVector);
          count[1]++;
        }
        cnt++;
      }
      double res0 = (likelihood[0] / ((double)count[0])) - (clusterScore[readPool_curSeg.curSeg->ID[0]] / ((double)clusterSize[readPool_curSeg.curSeg->ID[0]]));
      double res1 = (likelihood[1] / ((double)count[1])) - (clusterScore[readPool_curSeg.curSeg->ID[0]] / ((double)clusterSize[readPool_curSeg.curSeg->ID[0]]));
      double tS0 = ((double)(readPool_curSeg.curSeg->firstFrame))/100.0;
      double tS1 = ((double)(readPool_curSeg.curSeg->lastFrame - mxOut))/100.0;
      double tS2 = tS0;
      double tL0 = ((double)(mxIn))/100.0;
      double tL1 = ((double)(mxOut))/100.0;
      double tL2 = ((double)(readPool_curSeg.curSeg->lastFrame - readPool_curSeg.curSeg->firstFrame + 1))/100.0;
      PRINTF6("Result: (%f / %d) - (%f / %d)  = %f\t",likelihood[0],count[0],clusterScore[readPool_curSeg.curSeg->ID[0]]
              ,clusterSize[readPool_curSeg.curSeg->ID[0]],res0);
      PRINTF6("Result: (%f / %d) - (%f / %d)  = %f\n",likelihood[1],count[1],clusterScore[readPool_curSeg.curSeg->ID[0]]
              ,clusterSize[readPool_curSeg.curSeg->ID[0]],res1);
      if(count[0] > 0)
      {
        if(res0 < 0.0)
        {
          featurePool->addSegment(overlapSeg,readPool_curSeg.curSeg->firstFrame, readPool_curSeg.curSeg->firstFrame + mxIn - 1, overlapModelID);
//          fprintf(file,"SPEAKER %s 1 %.3f %.3f <NA> <NA> SPK01 <NA>\n",label,tS0,tL0);
        }
        else
        {
//          fprintf(file,"SPEAKER %s 1 %.3f %.3f <NA> <NA> SPK10 <NA>\n",label,tS0,tL0);
        }
      }
      if(count[1] > 0)
      {
        if(res1 < 0.0)
        {
          featurePool->addSegment(overlapSeg,readPool_curSeg.curSeg->lastFrame - mxOut, readPool_curSeg.curSeg->lastFrame - 1, overlapModelID);
//          fprintf(file,"SPEAKER %s 1 %.3f %.3f <NA> <NA> SPK02 <NA>\n",label,tS1,tL1);
        }
        else
        {
//          fprintf(file,"SPEAKER %s 1 %.3f %.3f <NA> <NA> SPK20 <NA>\n",label,tS1,tL1);
        }
      }
    }
/*    else
    {
      double tS = ((double)(readPool_curSeg.curSeg->firstFrame))/100.0;
      double tL = ((double)(readPool_curSeg.curSeg->lastFrame - readPool_curSeg.curSeg->firstFrame + 1))/100.0;
      fprintf(file,"SPEAKER %s 1 %.3f %.3f <NA> <NA> SPK03 <NA>\n",label,tS,tL);
    }*/
    trainVector = featurePool->getFirstVectorNextSegment(&readPool_curSeg,true, &isLast);
    vectorNotNull = (trainVector != NULL);
  }
  trainCluster[overlapModelID]->setTrainingData(featurePool,overlapSeg, overlapModelID);
  trainCluster[overlapModelID]->train(5,(CLUSTER_MIN_NUMBER_FRAMES != 3),true);
  PRINTF("Done training the overlap model..\n");fflush(stdout);
  featurePool->freeSegmentationID(overlapSeg);

  trainCluster[overlapModelID]->setTrainingData(featurePool,segID,overlapModelID);


  PRINTF("Returning training data to pool..\n");fflush(stdout);
  createOverlapTree(CLUSTER_MIN_NUMBER_FRAMES);

  for(int i=0;i<3;i++)
  {
    segmentFeaturePool(sadID_train,0,segID);
    double score = trainClusters();
    trainCluster[overlapModelID]->train(5,(CLUSTER_MIN_NUMBER_FRAMES != 3),true);
    PRINTF2("Overall score with overlap model: %f\n",score); fflush(stdout);
  }

  char filename[1000];
  strcpy(filename,outputFileName);
  strcat(filename,"overlap.am");
  FILE *outFile = fopen(filename,"wb");
  trainCluster[overlapModelID]->writeModel(outFile);
  fclose(outFile);

  if(mergeClusters(0,false,false))
  {
    PRINTF("THE OVERLAP MODEL WANTS TO BE MERGED!\n");
  }

  createLexicalTree(100);
  PRINTF("Segmenting feature pool (without overlap)..\n");fflush(stdout);
  segmentFeaturePool(sadID_decode,0,segID);

  if(file != NULL)
  {
    featurePool->writeRTTM(segID,label,"SPK",file,0,1);
  }

  createOverlapTree(100);
  PRINTF("Segmenting feature pool (with overlap)..\n");fflush(stdout);
  segmentFeaturePool(sadID_decode,0,segID);
  PRINTF2("Done segmenting with overlap model (%d)...\n", overlapModelID);fflush(stdout);

//   int cleanSadSeg = featurePool->claimSegmentationID();
//   featurePool->segmentationSubstract(sadID_train,segID,cleanSadSeg,overlapModelID);
//   char segmentString[60];
//   strcpy(segmentString,"SPEECH");
//   strcpy(&(segmentString[20]),"SIL");
//   strcpy(&(segmentString[40]),"SOUND");
//   featurePool->writeRTTM(cleanSadSeg,label,segmentString,file,20);
//   featurePool->freeSegmentationID(cleanSadSeg);

  featurePool->segmentationCreateOverlap(segID,overlapModelID);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Before we check which model to merge, we do some preparations in this method.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Train_Segmenter::startMergeIteration()
{
  for(int i=1;i<numberOfWords;i++)
  {
    for(int i2=1;i2<numberOfWords;i2++)
    {
      compareClusters[i+i2*numberOfWords] = 0;//(((warpIndex[i]+3) >= warpIndex[i2]) && ((warpIndex[i]-3) <= warpIndex[i2]));
    }
  }
  return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method reduces the network with a single state (cluster).
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool Train_Segmenter::mergeClusters(int maxClusters, bool smallestWins, bool actuallyMerge)
{
//  bool   merged      = false;

  if(numberOfClusters <= 2)  // +1 for administration...
  {
    return false;
  }

  nrMergeIterations++;

  int nrOfModels = 0;
  for(int i=1;i<numberOfWords;i++)
  {
    if((trainCluster[i] != NULL) && mergeable[i])
    {
      nrOfModels++;
    }
  }
  if(nrOfModels <= 1)
  {
    PRINTF("Done, we only have one model left!\n");
    return false;
  }
  PRINTF("Start merging\n");


  int    nrChannels = 1;
  if(featurePool->getStreamInfo() != NULL)
  {
    nrChannels = featurePool->getStreamInfo()->numberOfChannels;
  }

  startMergeIteration();

  double streamScore[nrChannels*numberOfWords*numberOfWords];
  double mergeScore[numberOfWords][numberOfWords];
  for(int i=0;i<numberOfWords;i++)
  {
    for(int i2=0;i2<numberOfWords;i2++)
    {
      mergeScore[i][i2] = -42.42e30;
      for(int a=0;a<nrChannels;a++)
      {
        streamScore[(i*numberOfWords+i2)*nrChannels+a] = -42.42e30;
      }
    }
  }

  for(int i=1;i<numberOfWords;i++)
  {
    if((trainCluster[i] != NULL) && mergeable[i])
    {
      for(int i2=i+1;i2<numberOfWords;i2++)
      {
        if((trainCluster[i2] != NULL) && mergeable[i2])
        {
          getMergeModelScore(i,i2,&(streamScore[(i*numberOfWords+i2)*nrChannels]));
        }
      }
    }
  }

  if(nrChannels > 1)
  {
    // Now determine the other part of the new weights: by BIC comparison:
    double streamVariance[nrChannels];
    int    nrDataPoints[nrChannels];
    Gaussian *posG[nrChannels];
    for(int i=0;i<nrChannels;i++)
    {
      streamVariance[i] = 1.0e300;
      posG[i]           = new Gaussian(1);
      nrDataPoints[i]   = 0;
    }
    for(int cd1 = 1; cd1 < numberOfWords; cd1++)
    {
      if((trainCluster[cd1] != NULL) && mergeable[cd1])
      {
        for(int cd2 = cd1+1; cd2< numberOfWords; cd2++)
        {
          if((trainCluster[cd2] != NULL) && mergeable[cd2])
          {
            for(int a=0;a<nrChannels;a++)
            {
              assert(streamScore[a + (cd1*numberOfWords+cd2)*nrChannels] != -42.42e30);
              streamScore[a + (cd1*numberOfWords+cd2)*nrChannels] -= bic_meanoffset[a];
              double score = streamScore[a + (cd1*numberOfWords+cd2)*nrChannels];
              score = sqrt(score*score);
              Vector *v = new Vector(1);
              v->setValue(0,score);
              posG[a]->train(1.0,v);
              nrDataPoints[a]++;
              delete v;
            }
          }
        }
      }
    }
    double m[nrChannels];
    for(int a=0;a<nrChannels;a++)
    {
      posG[a]->trainFinish(true);
      double var = posG[a]->getVariance()->getValue(0);
      m[a] = posG[a]->getMean()->getValue(0);
      streamVariance[a] = sqrt(m[a]*m[a] + var);
      delete posG[a];
      posG[a] = NULL;
    }

    double total = 0.0;
    int    minDataPoints = nrDataPoints[0];
    for(int a=0;a<nrChannels;a++)
    {
      PRINTF3("Stream variance: %d = %f\n",a,streamVariance[a]);
      streamVariance[a] = 1.0 / streamVariance[a];
      total += streamVariance[a];
      if(nrDataPoints[a] < minDataPoints)
      {
        minDataPoints = nrDataPoints[a]; 
      }
    }
    double weight[nrChannels];
    PRINTF("Quality Weight : ");
    for(int a=0;a<nrChannels;a++)
    {
      weight[a] = streamVariance[a] / total;
      PRINTF2("%f  ",weight[a]);
    }
    PRINTF("\n");

    if(nrMergeIterations < 3)
    {
      bic_meanoffset[0] = 0.0;
      for(int a=1;a<nrChannels;a++)
      {
        bic_meanoffset[a] = (m[a] - m[0]);

//        bic_meanoffset[a] = pass2Pool->getVector(NULL,pass2PoolindexShift)->getValue(a+1);

      }
//      pass2PoolindexShift++;

      PRINTF2("OFFSET CHANNEL 1: %f\n",bic_meanoffset[1]);


//       for(int a=0;a<nrChannels;a++)
//       {
//         weight[a] = pass2Pool->getVector(NULL,pass2PoolindexWeight)->getValue(a);
//       }
//       pass2PoolindexWeight++;

      featurePool->setStreamWeights(weight,true);
      PRINTF("Setting final Weights: ");
      for(int a=0;a<nrChannels;a++)
      {
        PRINTF2("%f  ",weight[a]);
      }
      PRINTF("\n"); 
    }
  }

  for(int i=1;i<numberOfWords;i++)
  {
    if((trainCluster[i] != NULL) && mergeable[i])
    {
      for(int i2=i+1;i2<numberOfWords;i2++)
      {
        if((trainCluster[i2] != NULL) && mergeable[i2])
        {
          mergeScore[i][i2] = 0.0;
          for(int c=0;c<nrChannels;c++)
          {
            assert(streamScore[c+(i*numberOfWords+i2)*nrChannels] != -42.42e30);
            mergeScore[i][i2] += featurePool->getStreamWeight(c)*streamScore[c+(i*numberOfWords+i2)*nrChannels];
          }
          if(smallestWins && featurePool->getClusterLength(segID,i) > 6000)
          {
            mergeScore[i][i2] = -9e300;
          }
          mergeScore[i2][i] = mergeScore[i][i2];
          PRINTF4("Total merge of %02d and %02d = %f\n",i,i2,mergeScore[i][i2]);
        }
      }
    }
  }

  int merge1[numberOfWords];
  int merge2[numberOfWords];
  int mapSpk[numberOfWords];
  int nrLeft = 0;
  for(int i=0;i<numberOfWords;i++)
  {
    if((trainCluster[i] != NULL) && mergeable[i])
    {
      mapSpk[i] = i;
      nrLeft++;
    }
    else
    {
      mapSpk[i] = -1;
    }
  }

  int pos = 0;
  double bestScore = 1.0;
  while(bestScore > 0.0)
  {
    bestScore = -9.0e200;
    for(int i=1;i<numberOfWords;i++)
    {
      if(((trainCluster[i] != NULL) && mergeable[i]) && (mapSpk[i] == i))
      {
        for(int i2=i+1;i2<numberOfWords;i2++)
        {
          if(((trainCluster[i2] != NULL) && mergeable[i2]) && (mapSpk[i2] == i2))
          {
            assert(mergeScore[i][i2] != -42.42e30);
            if(mergeScore[i][i2] > bestScore)
            {
              bestScore   = mergeScore[i][i2];
              merge1[pos] = i;
              merge2[pos] = i2;
            }
          }
        }
      } 
    }
#if (defined DIARIZATION_CHEAT_EXP3) || (defined DIARIZATION_CHEAT_EXP4) || (defined DIARIZATION_CHEAT_EXP5) 
//|| (defined GLOBAL_STOP)
    bestScore = 1.0;
#endif
    if(bestScore > 0.0 || (smallestWins && bestScore > -9.0e200) || ((nrOfModels - pos > maxClusters) && (bestScore > -9.0e200)))
    {
      bool doIt = true;
      if(!smallestWins)
      {
        for(int i=1;i<numberOfWords;i++)
        {
          if(mapSpk[i] == mapSpk[merge1[pos]])
          {
            for(int i2=1;i2<numberOfWords;i2++)
            {
              if(mapSpk[i2] == mapSpk[merge2[pos]])
              {
                assert(mergeScore[i][i2] != -42.42e30);
                int p1 = i;
                int p2 = i2;
                if(i2 < i)
                {
                  p1 = i2;
                  p2 = i;
                }
                if(mergeScore[p1][p2] <= 0.0)
                {
#if !(defined DIARIZATION_CHEAT_EXP3) && !(defined DIARIZATION_CHEAT_EXP4) && !(defined DIARIZATION_CHEAT_EXP5) 
//&& !(defined GLOBAL_STOP)
                  if(maxClusters <=0)
                  {
                    doIt = false;
                  }
#else
                  if(pos > 0)
                  {
                    doIt = false;
                  }
#endif
                }
              }
            }
          }
        }
      }
      mergeScore[merge1[pos]][merge2[pos]] = -9.0e200;
      mergeScore[merge2[pos]][merge1[pos]] = -9.0e200;
      if(doIt)
      {
        for(int i=1;i<numberOfWords;i++)
        {
          if(mapSpk[i] == mapSpk[merge2[pos]])
          {
            mapSpk[i] = mapSpk[merge1[pos]];
          }
        }
        PRINTF7("Merge %02d: SPK%02d - SPK%02d (%d,%d) %f\n",pos,merge1[pos],merge2[pos],mapSpk[merge1[pos]],mapSpk[merge2[pos]],bestScore);
        pos++;
        nrLeft--;
#if (defined DIARIZATION_CHEAT_EXP3) || (defined DIARIZATION_CHEAT_EXP4)
        fastMerge = 1;
#endif
        if(smallestWins)
        {
          fastMerge = 1;
        }
        if((maxClusters > 0) && (nrOfModels - pos > maxClusters))
        {
          bestScore = 1.0;
        }
        if(((nrMergeIterations > 2) && (fastMerge > 1)) || (fastMerge > 10))
        {
          if(nrLeft > 20)
          {
            if(pos>=fastMerge)
            {
              bestScore = -100.0;
            }
          }
          else if(nrLeft >= 16)
          {
            if(pos>=fastMerge/2)
            {
              bestScore = -100.0;
            }
          }
          else if(pos>=1)
          {
            bestScore = -100.0;
          }
        }
        else
        {
          bestScore = -100.0;
        }
      }
    }
  }
  if(pos <= 0)
  {
    return false;
  }



#if (defined DIARIZATION_CHEAT_EXP3) || (defined DIARIZATION_CHEAT_EXP4)
  // Cheating!
  char unWidened[MAXSTR];
  sprintf(unWidened,"%s.%d.mergeChoise",outputFileName,globalIt-1);
  FILE *textClFile = fopen(unWidened,"wt");
  fprintf(textClFile,"SPK%02d\tSPK%02d\n",merge1[0],merge2[0]);
  fclose(textClFile);
  fflush(stdout);

  char fileName[2000];
  cheatSpk1 = -1;
  cheatSpk2 = -1;
  sprintf(fileName,"/home/huijbreg/nist/diarytools/cheatPickMerge.sh %s.%d %s",outputFileName,globalIt-1,label);
  int resNotUsed = system(fileName);
  sprintf(fileName,"%s.%d.mergemax",outputFileName,globalIt-1);
  FILE *cheatFile = fopen(fileName,"rt");
  if(cheatFile != NULL)
  {
    char str[200];
    char *str2;
    if(!feof(cheatFile))
    {
      char *resChar = fgets(str,200,cheatFile);
      assert(resChar != NULL);
      str[2] = 0;
      str2 = &(str[3]);
      cheatSpk1 = atoi(str);
      cheatSpk2 = atoi(str2);
      PRINTF3("CHEATING-MERGE: %d %d\n",cheatSpk1,cheatSpk2);
      if(cheatSpk1 == -1 || cheatSpk2 == -1)
      {
        return false;
      }
      else
      {
        mergeModels(cheatSpk1,cheatSpk2);
      }
    }
    else
    {
      USER_ERROR("Could not cheat!");
    }    
    fclose(cheatFile);
  }  
  else
  {
    USER_WARNING("Could not cheat 2!");
    return false;
  }
#else
  PRINTF2("WE HAVE %d MERGES TO PERFORM...\n",pos);
  if(actuallyMerge)
  {
    for(int i=0;i<pos;i++)
    {
      numberOfMerges++;
      mergeModels(merge1[i],merge2[i]);
    }
  }
#endif
  return true;
}









