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
#include "shoutconfig.h"
#include "trainphonemodel.h"
#include "shouttrainmodel.h"
#include "shout-misc.h"

using namespace WriteFileLittleBigEndian;

int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_TRAIN_MODEL, argc, argv);

  ShoutTrainModel myModel(myConfig.getStringValue(ARGUMENT_TRAINDIR),
                          myConfig.getStringValue(ARGUMENT_DECISIONTREE),
                          myConfig.getStringValue(ARGUMENT_AM_CREATE),
                          (int)(myConfig.getFloatValue(ARGUMENT_TRAINMODELNR)),
                          (int)(myConfig.getFloatValue(ARGUMENT_MAXGAUSSIANS)),
                          (int)(myConfig.getFloatValue(ARGUMENT_TRAINSILP)*100.0),
                          (int)(myConfig.getFloatValue(ARGUMENT_TRAINSILMAX)),
                          myConfig.getStringValue(ARGUMENT_DOFAST));
  return 0;
}

ShoutTrainModel::ShoutTrainModel(char *dataDir, char *decisionRulesName, char *existingAM, int modelNumber, int maxGaussians, int trainSilP,
                                 int trainSilMax, char *doFast)
{
  bool             doFastTraining    = (doFast != NULL && strlen(doFast) > 0);
  bool             createClusters    = false;
  int              numberOfModels    = 1;
  int              numberOfSil       = 1;
  int              numberOfClusters  = 1;
  int              vectorSize        = 0;
  DataStats       *data              = NULL;
  FeaturePool     *featurePool       = NULL;
  PhoneModel      *satModel          = NULL;
  bool             dataHasNoHeader   = false;  // Only for SRE models.
  char str[MAXSTR];
  TrainPhoneModel *model[modelNumber+1];
  for(int i=0;i<modelNumber+1;i++)
  {
    model[i] = NULL;
  }

  trainSilMax = -1;
  trainSilP   =  100;

  strcpy(str,dataDir);
  strcat(str,"/");
  strcat(str,PHONEFILE);

  FILE *checkData     = NULL;
  FILE *trainDataFile = NULL;
  FILE *phoneFile     = fopen64(str, "rb");
  if(phoneFile==NULL)
  {
    USER_WARNING("Could not read phones.lst, you'd better be training SRE models!");
    vectorSize     = 26;
    numberOfModels =  1;
    numberOfSil    =  1;
    data = new DataStats[numberOfModels];
    strcpy(data[0].name,"UBM");
    data[0].nrOfSamples[0] = 1;
    dataHasNoHeader        = true;
  }
  else
  {
    freadEndianSafe(&vectorSize,1,sizeof(vectorSize),phoneFile);
    freadEndianSafe(&numberOfModels,1,sizeof(numberOfModels),phoneFile);
    freadEndianSafe(&numberOfSil,1,sizeof(numberOfSil),phoneFile);
    data = new DataStats[numberOfModels];
  
    if(modelNumber >= numberOfModels || modelNumber < 0)
    {
      USER_ERROR("The model number is invalid");
    }
  
    // We will check all models, although we are only going to train one. Don't want to find out something is wrong at the final model..
    for(int i=0;i<numberOfModels;i++)
    {
      char name[11];
      freadEndianSafe(name,10,1,phoneFile);
      name[10] = 0;
      strcpy(data[i].name,name);
      freadEndianSafe(data[i].nrOfSamples,MAX_CLUSTERS,sizeof(int),phoneFile);
      strcpy(str,dataDir);
      strcat(str,"/");
      strncat(str,data[i].name,10);
  
      if(i == modelNumber)
      {
        printf("The model I'm going to train: %s (%d of %d) has %d training samples\n",name,i,numberOfModels,data[i].nrOfSamples[0]);
      }
      if(data[i].nrOfSamples[0] > 0)
      {
        checkData = fopen64(str,"rb");
        if(checkData == NULL)
        {
          USER_ERROR2("Data samples for one or more phones are missing:",data[i].name);
        }
        int cdmID                 = 0;
        int numberOfObservations  = 0;
        int startOfState2         = 0;
        int startOfState3         = 0;
        int contextLeft           = 0;
        int contextRight          = 0;
  
        freadEndianSafe(&cdmID,1,sizeof(cdmID),checkData);
        freadEndianSafe(&numberOfObservations,1,sizeof(numberOfObservations),checkData);
        freadEndianSafe(&startOfState2,1,sizeof(startOfState2),checkData);
        freadEndianSafe(&startOfState3,1,sizeof(startOfState3),checkData);
        freadEndianSafe(&contextLeft,1,sizeof(contextLeft),checkData);
        freadEndianSafe(&contextRight,1,sizeof(contextRight),checkData);
        Vector *v = new Vector(checkData);
        if(v == NULL)
        {
          USER_ERROR2("Could not read data for model",str);
        }
        delete v;
        fclose(checkData);
        checkData = NULL;
      }
    }
    fclose(phoneFile);
    phoneFile = NULL;
  }

  char rulestring[MAXSTR];
  int highestRule = -1;
  int decisionMatrix[MAX_NUMBER_OF_DECISION_RULES*numberOfModels];

  // Opening the decision tree rules file..
  FILE *rulesFile = fopen(decisionRulesName, "rt");
  if(rulesFile==NULL)
  {
    USER_WARNING("The decision rule file is invalid or does not exist, I am going to work without clustering! (this is only okay if you are training ART models)");
  }
  else
  {
    while(!feof(rulesFile))
    {
      if(fscanf(rulesFile,"%s",rulestring)!=0)
      {
        if(rulestring[0] == '$')
        {
          highestRule++;
          for(int i=0;i<numberOfModels;i++)
          {
            decisionMatrix[highestRule*numberOfModels+i] = 0;
          }
        }
        else
        {
          if(highestRule < 0 || highestRule >= MAX_NUMBER_OF_DECISION_RULES)
          {
            USER_ERROR("Either too little or too much rules in the decision rule file!\n");
          }
          int i=0;
          strcat(rulestring,"          ");
          while((strncmp(rulestring,data[i].name,10)!=0) && (i<numberOfModels))
          {
            i++;
          }
          if(i>=numberOfModels)
          {
            USER_ERROR2("An unknown phone was found in the decision rule file: ",rulestring);
          }
          decisionMatrix[highestRule*numberOfModels+i] = 1;
        }
      }
    }
    fclose(rulesFile);
    rulesFile = NULL;
  }

  // Opening existing model or creating a new one:

  FILE *modelFile = fopen(existingAM, "rb");
  if(modelFile!=NULL)
  {
    printf("Using existing model in the AM file as initial model...\n");
    int nrM   = 0;
    int vectS = 0;
    freadEndianSafe(&vectS,1,sizeof(vectS),modelFile);
    if(vectS != vectorSize)
    {
      USER_ERROR("The vector size of this model is not compatible. Please use the standard ASR vector size...");
    }
    freadEndianSafe(&nrM,1,sizeof(nrM),modelFile);
    freadEndianSafe(&numberOfClusters,1,sizeof(numberOfClusters),modelFile);
    if(numberOfClusters != 1)
    {
      printf("WARNING: The acoustic model file contains clustered-models (%d clusters).\n",numberOfClusters);
      printf("         All data will be used to train 1 set of models!\n");
    }
    numberOfClusters = 1;
    if(nrM != numberOfModels)
    {
      USER_ERROR("The number of models in this file is not the same\nas the number of models in the training data.\nSorry, stopping here...");
    }
    else
    {
      for(int i=0;i<=modelNumber;i++)
      {
        model[i] = new TrainPhoneModel(modelFile,vectorSize);  // Reading the phone model
      }
    }
    fclose(modelFile);
    modelFile = NULL;
  }
  else  // Create model
  {
    model[modelNumber] = new TrainPhoneModel(data[modelNumber].name,numberOfModels,numberOfModels,(modelNumber<numberOfSil),vectorSize);
  }

  assert(model[modelNumber] != NULL);

  
/*
  if(doSat != NULL && strlen(doSat) > 0)
  {
    FILE *modelFile = fopen(doSat, "rb");
    if(modelFile!=NULL)
    {
      int nrM   = 0;
      int vectS = 0;
      freadEndianSafe(&vectS,1,sizeof(vectS),modelFile);
      if(vectS != vectorSize)
      {
        USER_ERROR("The vector size of the SAT-org model is not compatible. Please use the standard ASR vector size...");
      }
      freadEndianSafe(&nrM,1,sizeof(nrM),modelFile);
      freadEndianSafe(&numberOfClusters,1,sizeof(numberOfClusters),modelFile);
      if(numberOfClusters != 1)
      {
        printf("WARNING: The acoustic model file contains clustered-models (%d clusters).\n",numberOfClusters);
        printf("         All data will be used to train 1 set of models!\n");
      }
      numberOfClusters = 1;
      if(nrM != numberOfModels)
      {
        USER_ERROR("The number of models in this file is not the same\nas the number of models in the training data.\nSorry, stopping here...");
      }
      else
      {
        for(int i=0;i<=modelNumber;i++)
        {
          if(satModel != NULL)
          {
            delete satModel;
            satModel = NULL;
          }
          satModel = new PhoneModel(modelFile,vectorSize);  // Reading the phone model
        }
      }
      fclose(modelFile);
      modelFile = NULL;
    }
  }
*/

  createClusters    = (model[modelNumber]->getStatistics()->nrOfContexts <= 3 && model[modelNumber]->getStatistics()->nrOfGaussians <= 1);
  int numberOfRules = highestRule+1;

  if(model[modelNumber]->getStatistics()->nrOfGaussians > maxGaussians)
  {
    printf("I'm done: the existing model contains more gaussians than I'm allowed to train...\n");
  }
  else   // Yes: Let's train!!!
  {
    if(featurePool != NULL)
    {
      delete featurePool;
    }
    featurePool = new FeaturePool(TRAIN_MAX_FEATURE_POOL_SIZE, vectorSize);
    featurePool->resetSegmentation(TRAIN_ID1);  // Train samples
    featurePool->resetSegmentation(TRAIN_CLUSTER_ID1);
    featurePool->resetSegmentation(TRAIN_CLUSTER_ID2);
    featurePool->resetSegmentation(TRAIN_CLUSTER_ID3);
    model[modelNumber]->setTrainingData(featurePool,TRAIN_ID1,-1,-1,trainSilP, trainSilMax);

    if(createClusters)
    {
      printf("Using the decision rules for clustering...\n");
      printf("Number of rules: %d, number of phones: %d\n",numberOfRules,numberOfModels);
      model[modelNumber]->setDecisionMatrix(numberOfModels,numberOfRules,decisionMatrix);
    }

    strcpy(str,dataDir);
    strcat(str,"/");
    strncat(str,data[modelNumber].name,10);
    trainDataFile = fopen64(str,"rb");
    if(trainDataFile == NULL)
    {
      USER_ERROR("The training directory is invalid: model file is missing.");
    }
    else
    {
      int numTot    = 0;
      int numTot2   = 0;
      int timeStamp = 0;
      bool stopLoop = false;
      while(!stopLoop && (!feof(trainDataFile) && timeStamp < TRAIN_MAX_FEATURE_POOL_SIZE))
      {
        int cdmID  = 0;
        int num    = 0;
        int startOfState2;
        int startOfState3;
        int contextLeft;
        int contextRight;
        if(dataHasNoHeader)
        {
          stopLoop = true;
          cdmID    = 0;
          num      = 0;
          fseek(trainDataFile,0,SEEK_END);
          num = ftell(trainDataFile) / (sizeof(float)*vectorSize);
          fseek(trainDataFile,0,SEEK_SET);
          startOfState2 = -1;
          startOfState3 = -1;
          contextLeft   =  0;
          contextRight  =  0;
        }
        else
        {
          freadEndianSafe(&cdmID,1,sizeof(cdmID),trainDataFile);  // Reading the cluster-dependend-model ID: cdmID
          cdmID = 0;                          // For now, the cdm-ID is not used for training, one model is created!!
          freadEndianSafe(&num,1,sizeof(int),trainDataFile);
          freadEndianSafe(&startOfState2,1,sizeof(startOfState2),trainDataFile);
          freadEndianSafe(&startOfState3,1,sizeof(startOfState3),trainDataFile);
          freadEndianSafe(&contextLeft,1,sizeof(contextLeft),trainDataFile);
          freadEndianSafe(&contextRight,1,sizeof(contextRight),trainDataFile);
        }
        if(!feof(trainDataFile))
        {
          if(timeStamp+num-1 < TRAIN_MAX_FEATURE_POOL_SIZE)
          {
            numTot+=num;
            numTot2++;
            featurePool->addSegment(TRAIN_ID1,timeStamp,timeStamp+num-1,0,contextLeft*numberOfModels+contextRight);
            featurePool->setCurSegmentVectors(TRAIN_ID1,trainDataFile, dataHasNoHeader);
            if(createClusters)
            {
              if((modelNumber<numberOfSil)) // SIL
              {
                featurePool->addSegment(TRAIN_CLUSTER_ID1,timeStamp,timeStamp+num-1,0,
                                              contextLeft*numberOfModels+contextRight,contextLeft,contextRight);
              }
              else
              {
                featurePool->addSegment(TRAIN_CLUSTER_ID1,timeStamp,               timeStamp+startOfState2-1,0,
                                              contextLeft*numberOfModels+contextRight,contextLeft,contextRight);
                featurePool->addSegment(TRAIN_CLUSTER_ID2,timeStamp+startOfState2, timeStamp+startOfState3-1,0,
                                              contextLeft*numberOfModels+contextRight,contextLeft,contextRight);
                featurePool->addSegment(TRAIN_CLUSTER_ID3,timeStamp+startOfState3, timeStamp+num-1,0,
                                              contextLeft*numberOfModels+contextRight,contextLeft,contextRight);
              }
            }
          }
          else
          {
            printf("WARNING: Hit the roof of the training pool after %d samples\n", numTot2);
            for(int nn=0;nn<num;nn++)
            {
              Vector *v = new Vector(trainDataFile,vectorSize);
              delete v;
            }
          }
          timeStamp += num;
        }
      }
      printf("Total data used: %d  %d\n",numTot,numTot2);
      fclose(trainDataFile);
      trainDataFile = NULL;

//       if(true)
//       {
//         SegmentationAdmin readPool_curSeg;
//         Gaussian *trainGaussian = new Gaussian(featurePool->getVector(&readPool_curSeg,0)->len());
//         bool isLast            = false;
//         Vector *vector = featurePool->getFirstVectorFirstSegment(&readPool_curSeg, TRAIN_ID1,-1,true, &isLast);
//         while(vector != NULL)
//         {
//           int time      = featurePool->getCurSegmentStart(&readPool_curSeg);
//           int length    = featurePool->getCurSegmentLen(&readPool_curSeg);
//           for(int aa=time;aa<time+length;aa++)
//           {
//             Vector *feaVector = featurePool->getVector(&readPool_curSeg,aa);
//             trainGaussian->train(1.0,feaVector);
//             int cntList = 0;
//           }
//           vector = featurePool->getFirstVectorNextSegment(&readPool_curSeg,true, &isLast);
//         }
//         FILE *fN = fopen("total.N.stats","wt");
//         FILE *fF = fopen("total.F.stats","wt");
//         FILE *fS = fopen("total.S.stats","wt");
//         trainGaussian->writeAccumulators(fN,fF,fS);
//         fprintf(fN,"\n");
//         fprintf(fF,"\n");
//         fprintf(fS,"\n");
//         fclose(fN);
//         fclose(fF);
//         fclose(fS);
//         delete trainGaussian;
//         exit(0);
//       }

      printf("Starting the training!\n");
      ModelStats *stats = model[modelNumber]->getStatistics();
      printf("Number of training samples: %d\n",stats->nrOfTrainOcc);
      printf("Going to train until max Gaussians: %d\n",maxGaussians);

      model[modelNumber]->train(maxGaussians,(modelNumber<numberOfSil), false, NULL, NULL, satModel, doFastTraining);

      printf("Done! Writing model to file...\n");
      strcpy(str,dataDir);
      strcat(str,"/");
      strncat(str,data[modelNumber].name,10);
      if(satModel != NULL)
      {
        strcat(str,".sat");
      }
      else
      {
        strcat(str,".am");
      }
      FILE *outFile = fopen(str,"wb");
      if(outFile == NULL)
      {
        USER_ERROR("Could not write in the training directory.");
      }
      if(satModel != NULL)
      {
        model[modelNumber]->writeSAT(outFile);
      }
      else
      {
        model[modelNumber]->writeModel(outFile);
      }
      fclose(outFile);
    }
  }
  if(satModel != NULL)
  {
    delete satModel;
    satModel = NULL;
  }
  for(int i=0;i<=modelNumber;i++)
  {
    if(model[i] != NULL)
    {
      delete model[i];
      model[i] = NULL;
    }
  }
  if(data != NULL)
  {
    delete[] data;
  }
  if(featurePool != NULL)
  {
    delete featurePool;
  }
}


ShoutTrainModel::~ShoutTrainModel()
{
}


