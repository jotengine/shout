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
#include <time.h>
#include <sys/time.h>

#include "standard.h"
#include "shoutprepareadapt.h"
#include "shout-misc.h"
#include "trainphonemodel.h"
#include "featurepool.h"
#include "shoutconfig.h"

using namespace StringFunctions;
using namespace WriteFileLittleBigEndian;

#define INITIAL_CLUSTERS         (40)
#define MINIMUM_ADAPT_DATA      (500)

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Handles the command parameters for ShoutPrepareAdapt. It is the main function of the 
/// shout_prepare_adapt application. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_PREPARE_ADAPT, argc, argv);
  exit(0);

  if(argc != 4)
  {
    printf("%s%s%s\n",WELCOME_HEADER1,CURRENT_VERSION,WELCOME_HEADER2);    
    printf("Usage: shout_prepare_adapt [Cluster-AM file] [segment list (raw vtln clusterID)] [output segment list]\n");
    printf("If the cluster-AM file does not exist, an initial clustering will be done using\n");
    printf("a rough clustering method. The actual cluster-am file can be created with shout_maketrainset\n");
    printf("and shout_train_master/helper. After this a new iteration with shout_prepare_adapt can be started.\n");
    printf("Output is a new segment list (path to the audio file (raw), vtln warping factor and cluster ID)\n");
    exit(0);
  }  
//  int clusters = atoi((char*)argv[3]);
  FILE *am = fopen((char*)argv[1],"rb");
  ShoutPrepareAdapt ShoutPrepareAdapt(am,(char*)argv[2], 4, (char*)argv[3]);
  return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

ShoutPrepareAdapt::ShoutPrepareAdapt(FILE *inFile, char *dblName, int nrClusters, char *outDblName)
{
  int numberOfClusters     = 0;
  PhoneModel **phoneModels = NULL;
  
  char str[MAXSTR];
  char rawFileName[MAXSTR];
  char speakerID[MAXSTR];
  char vtlnString[MAXSTR];
  FILE *dblFile = fopen(dblName,"rt");
  if(dblFile == NULL)
  {
    USER_ERROR("Could not read the DBL file!");
  }
  
  if(inFile == NULL)
  {
    printf("# going to create initial clustering (Cluster-AM file does not exist)\n");
    srand((unsigned)time(0)); 
    FeaturePool *featurePool = new FeaturePool(TRAIN_MAX_FEATURE_POOL_SIZE, ASR_DEFAULT_VECTORSIZE);
    
    FeaturePool *adaptPool[INITIAL_CLUSTERS];
    for(int i=0;i<INITIAL_CLUSTERS;i++)
    {
      adaptPool[i] = new FeaturePool(TRAIN_MAX_FEATURE_POOL_SIZE/INITIAL_CLUSTERS, ASR_DEFAULT_VECTORSIZE);
    }
    
    int nrLines = 0;
    while(!feof(dblFile))
    {
      char *resChar = fgets(str,MAXSTR,dblFile);
      nrLines++;
    }
    fseek(dblFile,0,SEEK_SET);
    // We will use at least 100 lines:
    int skipLines = nrLines/100 - 1;
    while(!feof(dblFile))
    {
      for(int i=0;i<skipLines;i++)
      {
        char *resChar = fgets(str,MAXSTR,dblFile);
        assert(resChar != NULL);
      }
      if(fscanf(dblFile,"%s %s %s %s\n",rawFileName,speakerID,vtlnString,str) != 0)
      {
        double vtln   = atof(vtlnString);
        featurePool->addFeatureFile(rawFileName,FILE_TYPE_RAW_DECODING,50,1000,vtln);
      }
      else
      {
        USER_ERROR("The DBL file is corrupt!");
      }
    }    
    fseek(dblFile,0,SEEK_SET);
    
    TrainPhoneModel *mainModel = new TrainPhoneModel("CLUSTER",1,1,true,ASR_DEFAULT_VECTORSIZE);
    int trainID = featurePool->claimSegmentationID();
    featurePool->resetSegmentation(trainID);
    featurePool->addSegment(trainID,0,featurePool->getPoolLength()-1,0);
    mainModel->setTrainingData(featurePool,trainID,0);
    mainModel->train(4,true);
    
    // Now start dividing data over the INITIAL_CLUSTERS clusters:
    skipLines = nrLines/1000 - 1;
    while(!feof(dblFile))
    {
      for(int i=0;i<skipLines;i++)
      {
        char *resChar = fgets(str,MAXSTR,dblFile);
        assert(resChar != NULL);
      }
      if(fscanf(dblFile,"%s %s %s %s\n",rawFileName,speakerID,vtlnString,str) != 0)
      {
        double vtln   = atof(vtlnString);
        int randval   = int(double(INITIAL_CLUSTERS)*rand()/(RAND_MAX));
        adaptPool[randval]->addFeatureFile(rawFileName,FILE_TYPE_RAW_DECODING,50,1000,vtln);
      }
      else
      {
        USER_ERROR("The DBL file is corrupt!");
      }
    }    
    fseek(dblFile,0,SEEK_SET);
    
    TrainPhoneModel *adaptModel[INITIAL_CLUSTERS];
    double           distance[INITIAL_CLUSTERS][INITIAL_CLUSTERS];

    int model1      = 0;
    int model2      = 0;
    double biggestD = 0.0;        
    for(int i=0;i<INITIAL_CLUSTERS;i++)
    {
      adaptModel[i] = NULL;
      printf("Model: %02d has %d data\n",i,adaptPool[i]->getPoolLength());
      if(adaptPool[i]->getPoolLength() >= MINIMUM_ADAPT_DATA)
      {
        int adaptID = adaptPool[i]->claimSegmentationID();
        adaptPool[i]->resetSegmentation(adaptID);
        adaptPool[i]->addSegment(adaptID,0,adaptPool[i]->getPoolLength()-1,0);
        
        adaptModel[i] = new TrainPhoneModel(mainModel,0);
        for(int a=0;a<4;a++)
        {
          adaptModel[i]->adapt_clear();
          adaptModel[i]->adapt_setAcTrain(0,adaptID,adaptPool[i]);
          adaptingModel(adaptModel[i]);
        }
        for(int i2=0;i2<i;i2++)
        {
          if(adaptModel[i2] != NULL)
          {
            distance[i][i2] = adaptModel[i]->getKLDistance(adaptModel[i2]);
          }
          else
          {
            distance[i][i2] = 0.0;
          }
          distance[i2][i] = distance[i][i2];
          if(biggestD < distance[i][i2])
          {
            biggestD = distance[i][i2];
            model1   = i;
            model2   = i2;
          }
          printf("Distance %02d %02d = %f\n",i,i2,distance[i][i2]);
        }
      }
      delete adaptPool[i];
      adaptPool[i] = NULL;
    }

    int model3       = 0;
    int model4       = 0;
    double biggestD3 = 0.0;        
    double biggestD4 = 0.0;        
            
    for(int i=0;i<INITIAL_CLUSTERS;i++)
    {
      if(i != model1 && i != model2)
      {
        if(distance[i][model1] < distance[i][model2])
        {
          if(biggestD3 < distance[i][model1])
          {
            biggestD3 = distance[i][model1];
            model3    = i;
          }
        }
        else
        {
          if(biggestD4 < distance[i][model2])
          {
            biggestD4 = distance[i][model2];
            model4    = i;
          }
        }
      }
    }

    printf("Winners %02d - %02d = %f\n",model1,model2,distance[model1][model2]);
    printf("Winners %02d - %02d = %f\n",model1,model3,distance[model1][model3]);
    printf("Winners %02d - %02d = %f\n",model1,model4,distance[model1][model4]);
    printf("Winners %02d - %02d = %f\n",model2,model3,distance[model2][model3]);
    printf("Winners %02d - %02d = %f\n",model2,model4,distance[model2][model4]);
    printf("Winners %02d - %02d = %f\n",model3,model4,distance[model3][model4]);

    // Preparing models:
    numberOfClusters = 4;
    phoneModels      = new PhoneModel*[numberOfClusters];
    phoneModels[0]   = (PhoneModel*)adaptModel[model1]; 
    phoneModels[1]   = (PhoneModel*)adaptModel[model2]; 
    phoneModels[2]   = (PhoneModel*)adaptModel[model3]; 
    phoneModels[3]   = (PhoneModel*)adaptModel[model4]; 
    adaptModel[model1] = NULL;    
    adaptModel[model2] = NULL;    
    adaptModel[model3] = NULL;    
    adaptModel[model4] = NULL;    
            
    for(int i=0;i<INITIAL_CLUSTERS;i++)
    {
      if(adaptModel[i] != NULL)
      {
        delete adaptModel[i];
      }
    }
    delete mainModel;
  }
  else
  {
    int vectorSize = 0;
    freadEndianSafe(&vectorSize,1,sizeof(vectorSize),inFile);
    freadEndianSafe(&numberOfClusters,1,sizeof(numberOfClusters),inFile);
    int notUsed = 0;
    freadEndianSafe(&notUsed,1,sizeof(notUsed),inFile);
    assert(notUsed == 1);
    phoneModels = new PhoneModel*[numberOfClusters];
    for(int i=0;i<numberOfClusters;i++)
    {
      phoneModels[i] = new PhoneModel(inFile,vectorSize);
    }
    fclose(inFile);
  
    if(numberOfClusters != nrClusters)
    {
      USER_ERROR("The Cluster-AM file is not compatible to the number of clusters requested!");
    }
  }

  FILE *outFile = fopen(outDblName,"wt");
  if(outFile == NULL)
  {
    USER_ERROR("Could not write to the output segment-list file!");
  }
  
  // Now we can assume that we have our phone models ready and rocking!
  printf("Determining new cluster distribution...\n");
  while(!feof(dblFile))
  {
    if(fscanf(dblFile,"%s %s %s %s\n",rawFileName,speakerID,vtlnString,str) != 0)
    {
      double vtln   = atof(vtlnString);
      FeaturePool *fPool = new FeaturePool(rawFileName,FILE_TYPE_RAW_DECODING,false,vtln);
      double like[numberOfClusters];
      for(int i=0;i<numberOfClusters;i++)
      {
        like[i] = 0.0;
      }
      for(int t=0;t<fPool->getPoolLength();t++)
      {
        for(int i=0;i<numberOfClusters;i++)
        {
          like[i] += phoneModels[i]->getLogPDFProbability(0,fPool->getVector(NULL,t));
        }
      }
      double max = -9.0e300;
      int    win = 0;
      for(int i=0;i<numberOfClusters;i++)
      {
        if(like[i] > max)
        {
          max = like[i];
          win = i;
        }
      }
      fprintf(outFile,"%s %s %s %d\n",rawFileName,speakerID,vtlnString,win);
      delete fPool;      
    }
    else
    {
      USER_ERROR("The DBL file is corrupt!");
    }
  }    
  fclose(outFile);
  fclose(dblFile);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor cleans up administration.
/////////////////////////////////////////////////////////////////////////////////////////////////////

ShoutPrepareAdapt::~ShoutPrepareAdapt()
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

void ShoutPrepareAdapt::adaptingModel(PhoneModel *model)
{
  int dim = model->dim();
  Adapt_AM_TreeNode *adaptationTree = new Adapt_AM_TreeNode(dim);
  adaptationTree->meanGaussian = new Gaussian(dim);  
  model->adapt_setInitialNode(adaptationTree); 
  adaptationTree->setHelperMatrices_1(((PhoneModel**)&model), 1,model->isSilModel()); 
  adaptationTree->setHelperMatrices_2();
  model->adapt_adapt();
  model->adapt_clear();
  delete adaptationTree;
}



