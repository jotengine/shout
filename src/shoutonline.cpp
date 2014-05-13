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

#include "shoutonline.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "standard.h"
#include "shoutconfig.h"
#include "featureextraction.h"
#include "shout-misc.h"
#include "shout_vtln.h"

#include <time.h>
#include <math.h>
#include <sys/time.h>

using namespace WriteFileLittleBigEndian;
using namespace StringFunctions;



#define MAX_ADAPT_LENGTH                  (24)
//#define MAX_ADAPT_LENGTH                   (8)
#define VTLN_ADAPT_MINSPEECHFRAMES       (250)


#ifdef MAX_NR_THREADS
#include <pthread.h>

struct Thread_online_Data
{
  int threadID;
  int quit;

  char configDir[MAXSTR];
  int  number;

  int result;
};

#define NR_ONLINE_THREADS   (1)
#define ONLINE_VTLN_THREAD  (0)


bool threadMayStartOnline[NR_ONLINE_THREADS];
bool vtlnThreadFinished;
pthread_t                   threadOnline[NR_ONLINE_THREADS];
Thread_online_Data          dataThreadOnline[NR_ONLINE_THREADS];
pthread_mutex_t             condition_mutexStartOnline[NR_ONLINE_THREADS];
pthread_cond_t              condition_threadStartOnline[NR_ONLINE_THREADS];
pthread_mutex_t             condition_mutexDoneOnline;
pthread_cond_t              condition_threadDoneOnline;


void *thread_online( void *ptr )
{
  int threadID = ((Thread_online_Data*)ptr)->threadID;
  double totalScores[VTLN_STEPS];
  for(int i=0;i<VTLN_STEPS;i++)
  {
    totalScores[i] = 0.0;
  }

  while(true)
  {
    pthread_mutex_lock( &condition_mutexStartOnline[threadID] );
    while(!threadMayStartOnline[threadID])
    {
      pthread_cond_wait( &condition_threadStartOnline[threadID],
                         &condition_mutexStartOnline[threadID] );
    }
    pthread_mutex_lock( &condition_mutexDoneOnline);
    vtlnThreadFinished = false;
    pthread_mutex_unlock( &condition_mutexDoneOnline);

    threadMayStartOnline[threadID] = false;
    pthread_mutex_unlock( &condition_mutexStartOnline[threadID] );

    Thread_online_Data *data = (Thread_online_Data*)ptr;

    if(data->quit != 0)
    {
      return NULL;
    }

    char audioIn[MAXSTR];
    char rttmIn[MAXSTR];
    char alignFile[MAXSTR];

    sprintf(audioIn,"%s/audio-%04d.raw",data->configDir,data->number);
    sprintf(rttmIn,"%s/hypothesis-%04d.rttm",data->configDir,data->number);
    sprintf(alignFile,"%s/hypothesis-%04d.alignment",data->configDir,data->number);

//     Shout_VTLN *doVtln = new Shout_VTLN("SPEAKER", "models/shout.vtln", rttmIn, audioIn, NULL, totalScores);
//     data->result = doVtln->getVtlnFactor();
//     delete doVtln;

     data->result = 5;

    // Acoustic model adaptation...
    // =========== First: forced alignment =====================================
//     if(data->number > 2)
//     {
// printf("There we go!!\n");fflush(stdout);
//     char shoutcomm[MAXSTR];
//     sprintf(shoutcomm,"~/shout-toolkit/optimized/src/shout --am-phone %s/adapted.am --align --meta-in %s --audio %s --dct  models/9904.65K.release-003.dct.bin --phone-based -hyp %s",data->configDir,rttmIn,audioIn,alignFile );
//     system(shoutcomm);
//     sprintf(shoutcomm,"~/shout-toolkit/optimized/src/shout_maketrainset -mi %s -tr %s/traindir -hyp %s -amt PHONE -pl models/ut-nl2006.phoneset -a %s > %s/temp.txt", data->configDir,rttmIn,alignFile,audioIn,data->configDir);
//     system(shoutcomm);
//     sprintf(shoutcomm,"~/shout-toolkit/optimized/src/shout_adapt_am -tr %s/traindir -amp %s/adapted.am -amo %s/adapted.am.out > %s/temp.txt; mv %s/adapted.am.out %s/adapted.am", data->configDir,data->configDir,data->configDir,data->configDir,data->configDir,data->configDir);
//     system(shoutcomm);
//     }

    pthread_mutex_lock( &condition_mutexDoneOnline);
    vtlnThreadFinished = true;
    pthread_cond_signal( &condition_threadDoneOnline);
    pthread_mutex_unlock( &condition_mutexDoneOnline);
  }
}

#endif

#define SILTHRESHOLD (20)

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Main function of the shout application.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_SHOUTONLINE, argc, argv);
  ShoutOnline ShoutOnline(myConfig.getStringValue(ARGUMENT_AM_PHONE), myConfig.getStringValue(ARGUMENT_DCT_BIN),
                          myConfig.getStringValue(ARGUMENT_LM_BIN),myConfig.getStringValue(ARGUMENT_LATTICE),myConfig.getStringValue(ARGUMENT_CONFIGDIR));
  return 0;
}


ShoutOnline::ShoutOnline(char *amName, char *lexName, char *lmName, char *latName, char *configDir)
{
  bool performAdaptation = ((configDir != NULL) && (strlen(configDir)>0));
  int  adaptCounter      = 1;
  int  adaptTotalSpeech  = 0;
  FILE *hypFile          = NULL;
  bool vtlnRunning       = false;
  int  vtlnFactor        = 6;

#ifdef MAX_NR_THREADS
  for(int i=0;i<NR_ONLINE_THREADS;i++)
  {
    threadMayStartOnline[i]         = false;
    dataThreadOnline[i].quit        = 0;
    dataThreadOnline[i].threadID    = i;
    strcpy(dataThreadOnline[i].configDir,configDir);
    pthread_create( &threadOnline[i], NULL, thread_online, (void*) &dataThreadOnline[i]);
  }
#endif

  bool doLattice = (latName != NULL && (strlen(latName) > 0));

  if((latName != NULL) && (lmName != NULL) && (strlen(latName) > 0) && (strlen(lmName) > 0))
  {
    USER_ERROR("You specified a language model and a lattice. Please only define one of them.");
  }

  if(performAdaptation)
  {
    char name[MAXSTR];
    sprintf(name,"%s/adapted.am",configDir);
    FILE *testFile = fopen(name,"rb");
    if(testFile != NULL)
    {
      fclose(testFile);
      strcpy(amName,name);
    }
  }

  FILE *modelFile = fopen(amName, "rb");
  if(modelFile==NULL)
  {
    USER_ERROR("No valid AM file passed!");
  }
  freadEndianSafe(&vectorSize,1,sizeof(vectorSize),modelFile);
  freadEndianSafe(&numberOfModels,1,sizeof(numberOfModels),modelFile);
  freadEndianSafe(&numberOfClusters,1,sizeof(numberOfClusters),modelFile);
  assert(numberOfClusters == 1);
  if(numberOfModels <= 0 || numberOfClusters <= 0)
  {
    USER_ERROR("This AM file does not contain any acoustic models!");
  }
  else
  {
    models = new PhoneModel*[numberOfModels];
    for(int i=0;i<numberOfModels;i++)
    {
      models[i] = new PhoneModel(modelFile,vectorSize,true);
    }
  }
  fclose(modelFile);

  FILE *lexFile = fopen(lexName, "rb");
  if(lexFile==NULL)
  {
    USER_ERROR("No valid dictionary (lexical tree) file passed!");
  }
  // Load the lexical tree:
  myTree = new LexicalTree(stdout, lexFile);
  fclose(lexFile);

  myTree->setAMs(models);

  if(doLattice)
  {
    FILE *latFile = fopen(latName, "rt");
    if(latFile==NULL)
    {
      USER_ERROR("Could not open the lattice file...");
    }
    myTree->setLattice(latFile);
    fclose(latFile);
    myTree->overwritePrunePars(PERFORM_PRUNING_HISTOGRAM,PERFORM_PRUNING_STATE,PERFORM_PRUNING_ENDSTATE,
                              200.0, 1.0, 1.0, false, 1, 30000);
    myTree->overwriteWeightPars(1.0,0.0,0.0);
  }
  else
  {
    myLM = new LanguageModel(lmName);
    myTree->setLM(myLM);
    myTree->overwritePrunePars(PERFORM_PRUNING_HISTOGRAM,PERFORM_PRUNING_STATE,PERFORM_PRUNING_ENDSTATE,
                              150.0, 40.0, 30.0, true, 160, 35000);
    myTree->overwriteWeightPars(25.0,0.0,0.0);
  }


  int  currentSlide = 0;
  int  nrSlides     = 13;
  int slide[nrSlides];

  slide[0] = myTree->getWordID("SHoUT");
  slide[1] = myTree->getWordID("nederlands");
  slide[2] = myTree->getWordID("onderzoek");
  slide[3] = myTree->getWordID("verdediging");
  slide[4] = myTree->getWordID("realtime");
  slide[5] = myTree->getWordID("online");
  slide[6] = myTree->getWordID("waterpolo");
  slide[7] = myTree->getWordID("Gonnagles");
  slide[8] = myTree->getWordID("todo-list");
  slide[9] = myTree->getWordID("afspraken");
  slide[10] = myTree->getWordID("todo-lijst");
  slide[11] = myTree->getWordID("stabiel");
  slide[12] = myTree->getWordID("welkom");

  char comm[MAXSTR];
  sprintf(comm,"gconftool -t string -s /desktop/gnome/background/picture_filename /home/huijbreg/slides/slide-%02d.png &",currentSlide+1);
  int resNotUsed = system(comm);


  myTree->printInitialSettings(amName, lexName, "", lmName, false);
  printf("I'm ready...\n\n\n\n\n");

  short int sample;
  FeatureExtraction *feaExtract = new FeatureExtraction(NULL,NULL,ASR_DEFAULT_MFCC_NUMBER,
                                                                  ASR_DEFAULT_MFCC_ENERGY, 
                                                                  ASR_DEFAULT_MFCC_ZEROCROSS, 
                                                                  ASR_DEFAULT_MFCC_DELTA,
                                                                  VTLN_START + VTLN_STEPSIZE*((double)vtlnFactor),
                                                                  ASR_DEFAULT_SQRT10, false, true);
  feaExtract->setNormalization(FEATURENORM_ALL,1,true,true,NULL);
  feaExtract->doFeatureExtraction("STDIN",NULL,true);

  if(performAdaptation)
  {
    char name[MAXSTR];
    sprintf(name,"%s/config.dat",configDir);
    hypFile = fopen(name,"rb");
    if(hypFile != NULL)
    {
      freadEndianSafe(&adaptCounter,1,sizeof(adaptCounter),hypFile);
      freadEndianSafe(&vtlnFactor,1,sizeof(vtlnFactor),hypFile);
      fclose(hypFile);
      hypFile = NULL;
    }
adaptCounter = 1000;
    vtlnFactor = 5;
    feaExtract->setVTLN(VTLN_START + VTLN_STEPSIZE*((double)vtlnFactor));
    if(adaptCounter < MAX_ADAPT_LENGTH)
    {
      sprintf(name,"%s/audio-%04d.raw",configDir,adaptCounter);
      feaExtract->setAudioFileOut(name);
      sprintf(name,"%s/hypothesis-%04d.rttm",configDir,adaptCounter);
      hypFile = fopen(name,"wt");
    }
  }

  double silenceThreshold = feaExtract->getMean()->getValue(ASR_DEFAULT_MFCC_NUMBER) +
                            feaExtract->getVariance()->getValue(ASR_DEFAULT_MFCC_NUMBER);

struct timeval time1;
struct timeval time2;
gettimeofday( &time1, NULL );


  silCounter = SILTHRESHOLD+1;
  int timeShow  = 0;
  int time      = 0;
  int timeStart = 0;
  int idList[100];
  int lineLength = 0;
  bool isSil[GMMLOOKAHEAD];
  double testSil[GMMLOOKAHEAD];
  int action[GMMLOOKAHEAD];
  myTree->initialiseTree(time);
  time++;
  Vector **v = new Vector*[GMMLOOKAHEAD];
  for(int la=0;la<GMMLOOKAHEAD;la++)
  {
    isSil[la] = (la==0 || isSil[la-1]);
    testSil[la] = silenceThreshold;
    v[la] = feaExtract->getOnlineVector(&testSil[la],isSil[la])->copyVector();
    action[la] = determineAction(testSil[la],&isSil[la]);
  }

  while(v[0] != NULL)
  {
    if(action[0] == 1)
    {
//      printf("\n<RECOGNIZED %g>\t",myTree->getBestRecognitionScore());
      while(lineLength > 0)
      {
        printf("\b");
        lineLength--;
      }
      int nr = myTree->getBestIDSequence(idList,100,true,false);
      if(hypFile != NULL)
      {
        fprintf(hypFile,"SPEAKER hypothesis-%04d %d %7.3f %7.3f <NA> <NA> SPEECH <NA> <s> <s> ", 
                adaptCounter,vtlnFactor, ((double)timeStart)/100.0,((double)time-timeStart)/100.0);
        adaptTotalSpeech += (time-timeStart); 
      }
      for(int i=nr-1;i>=0;i--)
      {
        if(idList[i] != 0)
        {
          printf("%s ",myTree->getWord(idList[i]));
          if(hypFile != NULL)
          {
            fprintf(hypFile,"%s ",myTree->getWord(idList[i]));
          }
        }
      }
      printf("\n");
      if(hypFile != NULL)
      {
        fprintf(hypFile,"\n");
      }

      if((nr > 1) && strcmp(myTree->getWord(idList[1]),"oktober")==0 && strcmp(myTree->getWord(idList[0]),"<s>") == 0)
      {
        gettimeofday( &time2, NULL );
        int milliSec = (time2.tv_sec*1000 + time2.tv_usec/1000) - (time1.tv_sec*1000 + time1.tv_usec/1000);
        printf("%d ms!!!!\n",milliSec);
        exit(0);
      }
    }
    if(action[0] == 2)
    {
      timeStart = time;
      myTree->initialiseTree(time);

//      printf("<LISTENING> \t"); fflush(stdout);
      lineLength = 0;
    }
    if(action[0] == 3)
    {
      // Speaker is silent... Let's check if we should adapt the VTLN factor:
      if(performAdaptation && !vtlnRunning && (adaptCounter <= MAX_ADAPT_LENGTH) && 
         (adaptTotalSpeech > VTLN_ADAPT_MINSPEECHFRAMES*adaptCounter))
      {
        assert(hypFile != NULL);
        char name[MAXSTR];
        adaptCounter++;
        fclose(hypFile);
        hypFile = NULL;
        feaExtract->setAudioFileOut(NULL);
        if(adaptCounter < MAX_ADAPT_LENGTH)
        {
          sprintf(name,"%s/audio-%04d.raw",configDir,adaptCounter);
          feaExtract->setAudioFileOut(name);
          sprintf(name,"%s/hypothesis-%04d.rttm",configDir,adaptCounter);
          hypFile = fopen(name,"wt");
        }
        adaptTotalSpeech = 0;
        vtlnRunning      = true;
        time             = 0;

#ifdef MAX_NR_THREADS
        pthread_mutex_lock( &condition_mutexStartOnline[ONLINE_VTLN_THREAD] );
        dataThreadOnline[ONLINE_VTLN_THREAD].number = adaptCounter-1;
        threadMayStartOnline[ONLINE_VTLN_THREAD] = true;
        pthread_cond_signal( &condition_threadStartOnline[ONLINE_VTLN_THREAD] );
        pthread_mutex_unlock( &condition_mutexStartOnline[ONLINE_VTLN_THREAD] );
#endif
      }
    }
    if(!isSil[0])
    {
      myTree->processVector(v,time);
      timeShow++;
      if(timeShow > 10)
      {
        int nr = myTree->getBestIDSequence(idList,15,false,true);
        while(lineLength > 0)
        {
          printf("\b");
          lineLength--;
        }
        if(nr == 15)
        {
          printf("... ");
          lineLength += 4;
        }
        for(int i=nr-1;i>=0;i--)
        {
          lineLength += strlen(myTree->getWord(idList[i]))+1;
          printf("%s ",myTree->getWord(idList[i])); fflush(stdout);
        }
        timeShow = 0;

        for(int i=0; i<nr;i++)
        {
          for(int i2=0;i2<nrSlides;i2++)
          {
            if(slide[i2] == idList[i])
            {
              i=nr;
              if(currentSlide != i2)
              {
                currentSlide = i2;
                char comm[MAXSTR];
                sprintf(comm,"gconftool -t string -s /desktop/gnome/background/picture_filename /home/huijbreg/slides/slide-%02d.png",i2+1);
                int resNotUsed = system(comm);
              }
            }
          }
        }
      }
    }
    // Check if vtln is finished:
    if(performAdaptation && vtlnRunning)
    {

#ifdef MAX_NR_THREADS
      pthread_mutex_lock( &condition_mutexDoneOnline );
      if(vtlnThreadFinished)
      {
        vtlnRunning = false;
        vtlnFactor = dataThreadOnline->result;
        vtlnThreadFinished = false;
      }
      pthread_mutex_unlock( &condition_mutexDoneOnline );
#endif
      if(!vtlnRunning)
      {
        feaExtract->setVTLN(VTLN_START + VTLN_STEPSIZE*((double)vtlnFactor));
        char name[MAXSTR];
        sprintf(name,"%s/config.dat",configDir);
        FILE *configFile = fopen(name,"wb");
        if(configFile != NULL)
        {
          fwriteEndianSafe(&adaptCounter,1,sizeof(adaptCounter),configFile);
          fwriteEndianSafe(&vtlnFactor,1,sizeof(vtlnFactor),configFile);
          fclose(configFile);
        }
      }
    }

    time++;

    delete v[0];
    for(int i=1;i<GMMLOOKAHEAD;i++)
    {
      v[i-1]       = v[i];
      isSil[i-1]   = isSil[i];
      testSil[i-1] = testSil[i];
      action[i-1]  = action[i];
    }
    testSil[GMMLOOKAHEAD-1] = silenceThreshold;
    isSil[GMMLOOKAHEAD-1]   = isSil[GMMLOOKAHEAD-2];
    v[GMMLOOKAHEAD-1] = feaExtract->getOnlineVector(&testSil[GMMLOOKAHEAD-1],isSil[GMMLOOKAHEAD-1])->copyVector();
    action[GMMLOOKAHEAD-1]  = determineAction(testSil[GMMLOOKAHEAD-1],&isSil[GMMLOOKAHEAD-1]);
  }
  for(int i=0;i<GMMLOOKAHEAD;i++)
  {
    delete v[i];
  }
  delete[] v;
}


ShoutOnline::~ShoutOnline()
{
#ifdef MAX_NR_THREADS
  for(int i=0;i<NR_ONLINE_THREADS;i++)
  {
    dataThreadOnline[i].quit = 1;
    pthread_mutex_lock( &condition_mutexStartOnline[i] );
    threadMayStartOnline[i] = true;
    pthread_cond_signal( &condition_threadStartOnline[i] );
    pthread_mutex_unlock( &condition_mutexStartOnline[i] );
    pthread_join( threadOnline[i], NULL);
  }
#endif

  /* We are done. Let's clean up: */
  if(models!=NULL)
  {
    for(int i=0;i<numberOfModels;i++)
    {
      delete models[i];
    }
    delete[] models;
  }
  if(myLM!=NULL)
  {
    delete myLM;
  }
  if(myTree != NULL)
  {
    delete myTree;
  }
}


int ShoutOnline::determineAction(double testSil, bool *isSil)
{
  int action = 0;

  if(!*isSil)
  {
    if(testSil < 0.0)
    {
      silCounter++;
    }
    else
    {
      silCounter = 0;
    }
    if(silCounter > SILTHRESHOLD)
    {
      action = 1; // Recognized, now it's silence.
      *isSil = true;
    }
  }
  else
  {
    if(testSil > 0.0)
    {
      silCounter = 0;
      *isSil = false;
      action = 2; // Start listening..
    }
    else
    {
      action = 3; // silence (check adaptation)..
    }
  }
  return action;
}


