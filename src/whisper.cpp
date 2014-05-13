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
#include "whisper.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include "shout_maketrainset.h"
#include "adapt_am_treenode.h"
#include "shout-misc.h"
#include "featurepool.h"
#include "articulatorystream.h"

extern double weights[37][40];

using namespace StringFunctions;
using namespace WriteFileLittleBigEndian;

extern bool doPhoneAlignment;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Main function of the shout application.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_SHOUT, argc, argv);

  Whisper myWhisper(&myConfig);
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// All code is run in this constructor. It loads the LanguageModel, the PhoneModel and the LexicalTree.
/// The decoding is done in the LexicalTree.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Whisper::Whisper(ShoutConfig *myConfig)
{
  // Decoder settings from configuration file:
  double myLM_SCALE                 = (double) myConfig->getFloatValue(ARGUMENT_LM_SCALE);
  double myTRANS_PENALTY            = (double) myConfig->getFloatValue(ARGUMENT_TRANS_PENALTY);
  double mySIL_PENALTY              = (double) myConfig->getFloatValue(ARGUMENT_SIL_PENALTY);
  double myBEAM                     = (double) myConfig->getFloatValue(ARGUMENT_BEAM);
  double mySTATE_BEAM               = (double) myConfig->getFloatValue(ARGUMENT_STATE_BEAM);
  double myEND_STATE_BEAM           = (double) myConfig->getFloatValue(ARGUMENT_END_STATE_BEAM);
  int    myHISTOGRAM_STATE_PRUNING  = (int)    myConfig->getFloatValue(ARGUMENT_HISTOGRAM_STATE_PRUNING);
  int    myHISTOGRAM_PRUNING        = (int)    myConfig->getFloatValue(ARGUMENT_HISTOGRAM_PRUNING);
  bool   myLMLA                     = myConfig->parIsDefined(ARGUMENT_LMLA);
  bool   doForcedAlignment          = myConfig->parIsDefined(ARGUMENT_ALIGN);
  bool   doLatticeRescoring         = myConfig->parIsDefined(ARGUMENT_LATTICERESCORE);
  bool   doAnalyse                  = myConfig->parIsDefined(ARGUMENT_ANALYSE);
  bool   doConfidence               = myConfig->parIsDefined(ARGUMENT_PHONECONF);
  bool   outputXML                  = myConfig->parIsDefined(ARGUMENT_XML);
         doPhoneAlignment           = myConfig->parIsDefined(ARGUMENT_PHONEBASED);
  char *latticeTraining             = myConfig->getStringValue(ARGUMENT_LATTICETRAINING);
  char *amName                      = myConfig->getStringValue(ARGUMENT_AM_PHONE);
  char *rawName                     = myConfig->getStringValue(ARGUMENT_AUDIOFILE);
  char *backName                    = myConfig->getStringValue(ARGUMENT_BACKDCT_BIN);
  char *lexName                     = myConfig->getStringValue(ARGUMENT_DCT_BIN);
  char *lmName                      = myConfig->getStringValue(ARGUMENT_LM_BIN);
  char *rttmName                    = myConfig->getStringValue(ARGUMENT_META_IN);
  char *nbestDir                    = myConfig->getStringValue(ARGUMENT_NBEST);
  char *latticeDir                  = myConfig->getStringValue(ARGUMENT_LATTICE);
  char *tokenDistDir                = myConfig->getStringValue(ARGUMENT_TOKENDIST);
  char *hypName                     = myConfig->getStringValue(ARGUMENT_HYPFILE);
  char *method                      = myConfig->getStringValue(ARGUMENT_PERLINE);
  char *histNormName                = myConfig->getStringValue(ARGUMENT_AM_HISTNORM);

  if(strlen(latticeTraining) > 0)
  {
    doLatticeRescoring  = true;
  }
  bool addEndOfSentence = false;

  if(doAnalyse && !(doPhoneAlignment || (strcmp(latticeDir,"") != 0)))
  {
    USER_ERROR("For analysis, both phone alignment and lattice generation needs to be set");
  }

  if(doForcedAlignment && doLatticeRescoring)
  {
    USER_ERROR("Forced alignement and lattice rescoring cannot be done simultaneously");
  }

  if(doLatticeRescoring && (strlen(latticeDir) == 0 || strcmp(latticeDir,"NULL") == 0))
  {
    USER_ERROR("For lattice rescoring, the lattice directory must be defined");
  }

  // Other variables:
  char tokenDistFile[MAXSTR];
  char latticeFile[MAXSTR];

  // Initialise:
  myBackgroundTree                  = NULL;
  models                            = NULL;
  myTree                            = NULL;
  myPhoneLoop                       = NULL;
  myLM                              = NULL;
  clusterModels                     = NULL;
  int clusterID                     = 0;
  FILE *tdFile                      = NULL; // Token Distribution File...
  FILE *outFile                     = NULL;


  // Load the language model:
  myLM = new LanguageModel(lmName);
  if(myLM->getNumberOfWords() <= 0)
  {
    printf("#\n# WARNING!!!!  No valid Language Model passed. Going to work without a LM!\n#\n");
    delete myLM;
    myLM = NULL;
  }

  FILE *modelFile = fopen(amName, "rb");
  if(modelFile==NULL)
  {
    USER_ERROR("No valid AM file passed!");
  }
  FILE *backFile = NULL;
  if(strlen(backName) > 0)
  {
    backFile = fopen(backName, "rb");
    if(backFile==NULL)
    {
      USER_ERROR("No valid background DCT (lexical tree) file passed!");
    }
  }
  FILE *rttmFile   = fopen(rttmName, "rb");
  if(rttmFile==NULL)
  {
    USER_ERROR("No valid META file passed!");
  }
  FILE *lexFile = fopen(lexName, "rb");
  if(lexFile==NULL)
  {
    USER_ERROR("No valid dictionary (lexical tree) file passed!");
  }

  if(strlen(hypName) != 0)
  {
    outFile = fopen(hypName,"wt");
  }
  if(outFile == NULL)
  {
    myTree = new LexicalTree(stdout, lexFile);
  }
  else
  {
    myTree = new LexicalTree(outFile, lexFile);
  }
  fclose(lexFile);

  if(doLatticeRescoring && myLMLA)
  {
    myLMLA            = false;
  }

  if(doForcedAlignment)
  {
    myLMLA = false;
    myTree->overwritePrunePars(PERFORM_PRUNING_HISTOGRAM,PERFORM_PRUNING_STATE,PERFORM_PRUNING_ENDSTATE,
                              myBEAM*5.0, 1.0, 1.0, false, 1, 30000);
    myTree->overwriteWeightPars(1.0,0.0,0.0);
  }
  else
  {
    myTree->overwritePrunePars(PERFORM_PRUNING_HISTOGRAM,PERFORM_PRUNING_STATE,PERFORM_PRUNING_ENDSTATE,
                              myBEAM, mySTATE_BEAM, myEND_STATE_BEAM, myLMLA, myHISTOGRAM_STATE_PRUNING, myHISTOGRAM_PRUNING);
    myTree->overwriteWeightPars(myLM_SCALE,myTRANS_PENALTY,mySIL_PENALTY);
  }

  if(backFile != NULL)
  {
    if(outFile == NULL)
    {
      myBackgroundTree = new LexicalTree(stdout, backFile);
    }
    else
    {
      myBackgroundTree = new LexicalTree(outFile, backFile);
    }
    fclose(backFile);
  }

  // Load the lexical tree:
  freadEndianSafe(&vectorSize,1,sizeof(vectorSize),modelFile);
  freadEndianSafe(&numberOfModels,1,sizeof(numberOfModels),modelFile);
  freadEndianSafe(&numberOfClusters,1,sizeof(numberOfClusters),modelFile);
  numberOfClusteredModels = numberOfModels*numberOfClusters;
  if(numberOfModels <= 0 || numberOfClusters <= 0)
  {
    USER_ERROR("This AM file does not contain any acoustic models!");
  }
  else
  {
    models = new PhoneModel*[numberOfClusteredModels];
    for(int i=0;i<numberOfClusteredModels;i++)
    {
      models[i] = new PhoneModel(modelFile,vectorSize,!(latticeTraining != NULL && strlen(latticeTraining) > 0));
    }
    if(numberOfClusters > 1)
    {
      clusterModels = new PhoneModel*[numberOfClusters];
      for(int i=0;i<numberOfClusters;i++)
      {
        clusterModels[i] = new PhoneModel(modelFile,vectorSize);
      }
    }
    else
    {
      fclose(modelFile);
    }
  }

  if(latticeTraining != NULL && strlen(latticeTraining) > 0)
  {
    for(int i=0;i<numberOfClusteredModels;i++)
    {
      models[i]->adapt_clear();
    }
  }


  int clusterIDTrain                  = -1;


  myTree->printInitialSettings(amName, lexName, backName, lmName, outputXML);


  for(int i=0;i<numberOfClusters;i++)
  {
    myTree->checkAMs(numberOfModels,&(models[i*numberOfModels]), outputXML);
  }

  myTree->setAMs(&(models[clusterID*numberOfModels]));
  if(myLM != NULL)
  {
    myTree->setLM(myLM);
  }


  myTree->setLatticeGeneration((strlen(latticeDir) > 0) || (strlen(nbestDir) > 0));


  if(doConfidence)
  {
    if(outFile == NULL)
    {
      myPhoneLoop = new LexicalTree(stdout);
    }
    else
    {
      myPhoneLoop = new LexicalTree(outFile);
    }
    myPhoneLoop->setPhoneLoop(numberOfModels,models);
    myPhoneLoop->overwritePrunePars(true,true,true,30.0,5.0,5.0,false,1,2000);
    myPhoneLoop->overwriteWeightPars(1.0,0.0,0.0);
  }

  char label[MAXSTR];
  char beginTime[MAXSTR];
  char endTime[MAXSTR];

  int totTime     = 0;
  int totMilliSec = 0;
  clusterIDTrain  = 0;

  FeaturePool *featurePool = new FeaturePool(rawName,FILE_TYPE_RAW_DECODING,false,1.0,false,true);  // Only prepare...
  int segID = featurePool->claimSegmentationID();
  featurePool->readRTTM(segID,NULL,-1,-1,3,rttmFile);
  fclose(rttmFile);
//  featurePool->setFeatureNormalization(CMVN_CLUSTER);
  featurePool->createNewPool(rawName, histNormName, NULL, false, false, segID);


  int sentenceID = 1;
  bool isLast   = false;
  Vector *vector = featurePool->getFirstVectorFirstSegment(&readPool_curSeg, segID,-1,true, &isLast);
  while(vector != NULL)
  {
    if(strcmp(latticeDir,"") != 0)
    {
      sprintf(label,"SPK00-%03d",sentenceID);
    }
    else
    {
      sprintf(label,"SPK%02d-%03d",featurePool->getSegmentID(0,&readPool_curSeg),sentenceID);
    }
    sentenceID++;
    struct timeval time1;
    struct timeval time2;
    gettimeofday( &time1, NULL );
    int time      = featurePool->getCurSegmentStart(&readPool_curSeg);
    int length    = featurePool->getCurSegmentLen(&readPool_curSeg);

    sprintf(beginTime,"%0.2f",((double)time)/100.0);
    sprintf(endTime,  "%0.2f",((double)time+length)/100.0);

    bool notOOV = true;
    if(doLatticeRescoring)
    {
      strcpy(latticeFile,latticeDir);
      strcat(latticeFile,"/");
      strcat(latticeFile,label);
      strcat(latticeFile,".shout-lattice");
      FILE *latFile = fopen(latticeFile,"rt");
      if(latFile == NULL)
      {
        USER_ERROR2("Could not read from:",latticeFile);
      }
      myTree->setLattice(latFile);
      fclose(latFile);
    }
    else if(doForcedAlignment || doAnalyse)
    {
      char alignString[MAXSTR];
      char alignRestString[MAXSTR];
      char word[MAXSTR];
      strcpy(alignString,featurePool->getCurInfoString(&readPool_curSeg));
      if(strlen(alignString) > 0)
      {
        myTree->setForcedAlign(NULL);   // This will initialise the tree!
        splitList(alignString,alignRestString);
        strcpy(word,alignString);
        strcpy(alignString,alignRestString);
        splitList(alignString,alignRestString);
        myTree->setInitialLMHistory(word,alignString);
        strcpy(alignString,alignRestString);
        while(notOOV && (strlen(alignString) > 0))
        {
          splitList(alignString,alignRestString);
          notOOV = myTree->setForcedAlign(alignString);
          if(!notOOV)
          {
            if(myBackgroundTree != NULL)
            {
              notOOV = myTree->addForcedAlignOOV(myBackgroundTree->borrowPhoneString(alignString));
              if(!notOOV)
              {
                USER_WARNING2("No alignment due to OOV:",alignString);
              }
            }
            else
            {
              notOOV = false;
              USER_WARNING2("No alignment due to OOV:",alignString);
            }
          }
          strcpy(alignString,alignRestString);
        }
      }
      else
      {
        USER_ERROR("Transcription needed for forced alignment is missing.");
      }
    }



    if(notOOV)
    {
      if(doAnalyse)
      {
        myLMLA = false;
//        doPhoneAlignment = true;
        myTree->overwritePrunePars(PERFORM_PRUNING_HISTOGRAM,PERFORM_PRUNING_STATE,PERFORM_PRUNING_ENDSTATE,
                                  myBEAM*5.0, 1.0, 1.0, false, 1, 30000);
        myTree->overwriteWeightPars(1.0,0.0,0.0);
        // First align:
        myTree->initialiseTree(time);
        for(int i=time;i<time+length;i++)
        {
          Vector **feaVector = featurePool->getVectorList(&readPool_curSeg,i);
          myTree->processVector(feaVector,i);
        }
        myLMLA = true;
//        doPhoneAlignment = false;
        myTree->overwritePrunePars(PERFORM_PRUNING_HISTOGRAM,PERFORM_PRUNING_STATE,PERFORM_PRUNING_ENDSTATE,
                                  myBEAM, mySTATE_BEAM, myEND_STATE_BEAM, myLMLA, myHISTOGRAM_STATE_PRUNING, myHISTOGRAM_PRUNING);
        myTree->overwriteWeightPars(myLM_SCALE,myTRANS_PENALTY,mySIL_PENALTY);
        // Now reset the tree and do a regular search:
        // Done the alignment. Safe the best path:
        if(!myTree->safeBestRecognition(addEndOfSentence))
        {
          printf("## ERROR: No alignment. It will not be possible to check the pruning beams for this sentence!\n");
        }
        myTree->setForcedAlign(NULL);   // This will initialise the tree!
      }

      if(latticeTraining != NULL && strlen(latticeTraining) > 0)
      {
        if(outFile != NULL)
        {
          fprintf(outFile, "Calculating accumulators for %s... ",label); fflush(outFile);
        }
        else
        {
          printf("Calculating accumulators for %s... ",label); fflush(stdout);
        }
        myTree->initialiseTree(time);
        int numberOfNodes  = myTree->latticeBaumWelch_numberNodes(NULL,0,true);
        int numberOfStates = numberOfNodes*3;
        int unmark = myTree->latticeBaumWelch_numberNodes(NULL,0,false);

        // Create and initialize administration:
        double *Q = new double[length];
        double *latticeBaumWelchAlfa = new double[length*numberOfStates];
        double *latticeLikelihood    = new double[length*numberOfStates];
        double *latticeBaumWelchBeta = new double[(length+1)*numberOfStates];
        for(int i=0;i<length;i++)
        {
          myTree->latticeBaumWelch_setLikelihoods(NULL,featurePool->getVector(&readPool_curSeg,time+i),i,numberOfStates,latticeLikelihood);
          unmark = myTree->latticeBaumWelch_numberNodes(NULL,0,true);
          if(unmark != numberOfNodes)
          {
            printf("Only unmarked %d nodes (out of %d)!\n",unmark,numberOfNodes);
          }
          unmark = myTree->latticeBaumWelch_numberNodes(NULL,0,false);

          for(int i2=0;i2<numberOfStates;i2++)
          {
            latticeBaumWelchAlfa[i*numberOfStates+i2] = 0.0;
            latticeBaumWelchBeta[i*numberOfStates+i2] = 0.0;
          }
          Q[i] = 0.0;
        }
        for(int i2=0;i2<numberOfStates;i2++)
        {
          latticeBaumWelchBeta[length*numberOfStates+i2] = 0.0;
        }
        // Forward pass:
        myTree->latticeBaumWelch_initForward((double*)latticeBaumWelchAlfa);
        for(int i2=0;i2<numberOfStates;i2++)
        {
          Q[0] += latticeBaumWelchAlfa[i2];
        }
        for(int i2=0;i2<numberOfStates;i2++)
        {
          latticeBaumWelchAlfa[i2] /= Q[0];
        }

        for(int i=1;i<length;i++)
        {
          myTree->latticeBaumWelch_forward(latticeLikelihood,NULL,0.0,latticeBaumWelchAlfa,i,numberOfStates);
          unmark = myTree->latticeBaumWelch_numberNodes(NULL,0,true);
          if(unmark != numberOfNodes)
          {
            printf("1. Only unmarked %d nodes (out of %d)!\n",unmark,numberOfNodes);
          }
          unmark = myTree->latticeBaumWelch_numberNodes(NULL,0,false);

          // Normalize
          int offset = i*numberOfStates;
          for(int i2=0;i2<numberOfStates;i2++)
          {
            Q[i] += latticeBaumWelchAlfa[offset+i2];
          }

          for(int i2=0;i2<numberOfStates;i2++)
          {
            latticeBaumWelchAlfa[offset+i2] /= Q[i];
          }
        }
 
       // Backward pass:
        myTree->latticeBaumWelch_initBackward((double*)latticeBaumWelchBeta,length*numberOfStates);
        double resArray[numberOfNodes];
        for(int i=0;i<numberOfNodes;i++)
        {
          resArray[i] = -1.0e30;
        }
        myTree->latticeBaumWelch_backward(latticeLikelihood,NULL,latticeBaumWelchBeta,length-1,numberOfStates, Q[length-1],resArray);
        for(int i=length-2;i>=0;i--)
        {
          for(int i2=0;i2<numberOfNodes;i2++)
          {
            resArray[i2] = -1.0e30;
          }
          myTree->latticeBaumWelch_backward(latticeLikelihood,NULL,latticeBaumWelchBeta,i,numberOfStates, Q[i],resArray);
        }
        // Calculate posteriors:
        for(int i=1;i<length;i++)
        {
          double *posteriors = new double[numberOfStates];
          for(int i2=0;i2<numberOfStates;i2++)
          {
            posteriors[i2] = 0.0;
          }
          myTree->latticeBaumWelch_calculatePosteriors(latticeLikelihood,NULL,1.0,latticeBaumWelchAlfa,latticeBaumWelchBeta,
                                                       posteriors,i,numberOfStates);
          int unmark = myTree->latticeBaumWelch_numberNodes(NULL,0,true);
          if(unmark != numberOfNodes)
          {
            printf("2. Posteriors init: Only unmarked %d nodes (out of %d)!\n",unmark,numberOfNodes);
          }
          unmark = myTree->latticeBaumWelch_numberNodes(NULL,0,false);

          // Normalize the posteriors:
          double normFact = 0.0;
          for(int i2=0;i2<numberOfStates;i2++)
          {
            posteriors[i2] = exp(log(posteriors[i2])/(myLM_SCALE*2.0));
            normFact += posteriors[i2];
          }
          if(normFact > TRAIN_WEIGHT_MINIMUM)
          {
            for(int i2=0;i2<numberOfStates;i2++)
            {
              posteriors[i2] /= normFact;
            }
            myTree->latticeBaumWelch_mmi_accumulatorsPosteriors(NULL,posteriors,numberOfStates,featurePool->getVector(&readPool_curSeg,time+i));
//            myTree->latticeBaumWelch_printPosteriors(NULL,posteriors,i,numberOfStates,time);
            int unmark = myTree->latticeBaumWelch_numberNodes(NULL,0,true);
            if(unmark != numberOfNodes)
            {
              printf("3. Posteriors: Only unmarked %d nodes (out of %d)!\n",unmark,numberOfNodes);
            }
            myTree->latticeBaumWelch_numberNodes(NULL,0,false);

          }
          delete[] posteriors;
        }
        delete[]Q;
        delete[]latticeLikelihood;
        delete[]latticeBaumWelchAlfa;
        delete[]latticeBaumWelchBeta;

        gettimeofday( &time2, NULL );
        int milliSec = (time2.tv_sec*1000 + time2.tv_usec/1000) - (time1.tv_sec*1000 + time1.tv_usec/1000);
        if(outFile != NULL)
        {
          fprintf(outFile, "In %f times realtime (%d ms, %d frames)\n",
            ((double)milliSec/10.0)/((double)length),milliSec,length); fflush(outFile);
        }
        else
        {
          printf("In %f times realtime (%d ms, %d frames)\n",
            ((double)milliSec/10.0)/((double)length),milliSec,length); fflush(stdout);
        }
      }
      else
      {
        if(!(doForcedAlignment && myTree->alignmentIsEmpty()))
        {
          myTree->initialiseTree(time);
          for(int i=time;i<time+length;i++)
          {
            Vector **feaVector = featurePool->getVectorList(&readPool_curSeg,i);
            myTree->processVector(feaVector,i);
          }
        }

        float *myPhoneLoopProb = NULL;
        if(doConfidence)
        {
          USER_ERROR("Sorry, not implemented right now!");
/*          myPhoneLoopProb = new float[length];
          myTree->setPhoneLoopConfidence(myPhoneLoopProb, time);        // For printing
          myPhoneLoop->setPhoneLoopConfidence(myPhoneLoopProb, time);   // For filling
  
          myPhoneLoop->initialiseTree();
          for(int i=0;i<length;i++)
          {
            Vector *feaVector = featurePool->getVector(time+i);
            myPhoneLoop->processVector(feaVector,i);
            myPhoneLoop->storePLConfidence(i);
          }
          myPhoneLoop->initialiseTree();*/
        }
  
        gettimeofday( &time2, NULL );
        int milliSec = (time2.tv_sec*1000 + time2.tv_usec/1000) - (time1.tv_sec*1000 + time1.tv_usec/1000);
        totMilliSec += milliSec;
        totTime     += length;

        if(outputXML)
        {
          myTree->getBestRecognition(addEndOfSentence,true,true,label,milliSec,length,beginTime,endTime);
        }
        else
        {
          myTree->getBestRecognition(addEndOfSentence,false,true,label,milliSec,length,beginTime,endTime);
          myTree->getBestRecognition(addEndOfSentence,false,false,label);
        }
        if(doConfidence)
        {
          delete[] myPhoneLoopProb;
        }
        if(!doLatticeRescoring && strlen(latticeDir) > 0 && strcmp(latticeDir,"NULL") != 0)
        {
          strcpy(latticeFile,latticeDir);
          strcat(latticeFile,"/");
          strcat(latticeFile,label);
          strcat(latticeFile,".shout-lattice");
          FILE *latFile = fopen(latticeFile,"wt");
          if(latFile == NULL)
          {
            USER_ERROR2("Could not write to:",latticeFile);
          }
          myTree->printLattice(latFile, label,time+length);
          fclose(latFile);
        }
        if(strlen(nbestDir) > 0)
        {
          strcpy(latticeFile,nbestDir);
          strcat(latticeFile,"/");
          strcat(latticeFile,label);
          strcat(latticeFile,".txt");
          FILE *latFile = fopen(latticeFile,"wt");
          if(latFile == NULL)
          {
            USER_ERROR2("Could not write to:",latticeFile);
          }
          myTree->printNBestList(latFile);
          fclose(latFile);
        }
      }
    }
    else
    {
      if(outFile != NULL)
      {
        fprintf(outFile, "## No alignment due to OOV of reference...\n");
        fprintf(outFile,"  (%s)\n", label);
        fflush(outFile);
      }
      else
      {
        printf("## No alignment due to OOV of reference...\n");
        printf("  (%s)\n", label);
        fflush(stdout);
      }
    }
    fflush(stdout);
    vector = featurePool->getFirstVectorNextSegment(&readPool_curSeg,true, &isLast);
  }
  myTree->initialiseTree();


  myTree->printFinalSettings(outputXML, (int)totMilliSec,(int)totTime);

  if(outFile != NULL)
  {
    fclose(outFile);
  }

  if(latticeTraining != NULL && strlen(latticeTraining) > 0)
  {
    strcpy(latticeFile,latticeTraining);
/*    FILE *latFile = fopen(latticeFile,"rb");
    if(latFile != NULL)
    {
      for(int i=0;i<numberOfClusteredModels;i++)
      {
        models[i]->addAccumulators(latFile);
      }
      fclose(latFile);
    }*/
    FILE *latFile = fopen(latticeFile,"wb");
    for(int i=0;i<numberOfClusteredModels;i++)
    {
      models[i]->writeAccumulators(latFile);
    }
    fclose(latFile);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor cleans up all models.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Whisper::~Whisper()
{
  /* We are done. Let's clean up: */
  if(clusterModels!=NULL)
  {
    for(int i=0;i<numberOfClusters;i++)
    {
      delete clusterModels[i];
    }
    delete[] clusterModels;
  }
  if(models!=NULL)
  {
    for(int i=0;i<numberOfClusteredModels;i++)
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
  if(myBackgroundTree != NULL)
  {
    delete myBackgroundTree;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Determine the cluster ID
/////////////////////////////////////////////////////////////////////////////////////////////////////

int Whisper::getClusterID(FeaturePool *fPool)
{
  assert(fPool != NULL);
  int bestID = 0;
  if(clusterModels != NULL)
  {
    double like[numberOfClusters];
    for(int i=0;i<numberOfClusters;i++)
    {
      like[i] = 0.0;
    }
    Vector *feaVector = fPool->getVector(&readPool_curSeg,0);
    while(feaVector != NULL)
    {
      for(int i=0;i<numberOfClusters;i++)
      {
        like[i] += clusterModels[i]->getLogPDFProbability(0,feaVector);
      }
      feaVector = fPool->getNextVector(&readPool_curSeg);
    }
    double bestScore = -9.0e300;
    for(int i=0;i<numberOfClusters;i++)
    {
      if(like[i] > bestScore)
      {
        bestScore = like[i];
        bestID    = i;
      }
    }
  }
//  printf("# debug: clusterID = %02d\n",bestID);
  return bestID;
}


