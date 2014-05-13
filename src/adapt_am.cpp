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
#include "adapt_am.h"
#include "phonemodel.h"
#include "vector.h"
#include "train_segmenter.h"
#include "featurepool.h"
#include "shoutconfig.h"

#define PHONEFILE                "/phones.lst"

#include "shout-misc.h"
using namespace WriteFileLittleBigEndian;


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The constructor will check the directory that contains the adaptation data, load the data
/// and pass it to Adapt_AM_TreeNode. The resulting model will be stored to disc.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Adapt_AM::Adapt_AM(char *adapt_dir, char *amName, char *clusterAM, char *newAmName)
{
  char                str[MAXSTR];
  char                dataDir[MAXSTR];
  int                 numberOfModels;
  int                 numberOfSil;
  DataStats          *data   = NULL;
  models = NULL;
  FILE               *trainDataFile;
  int numberOfClusters = 1;
  int vectorSize       = ASR_DEFAULT_VECTORSIZE;
  
  TrainPhoneModel   **clusterModels = NULL;  

  
  FILE *modelFile = fopen(clusterAM,"rb");
  if(modelFile == NULL)
  {
    USER_WARNING("The clusterAM could not be openend. I'm assuming we're not using clustered adaptation...\n");
  }
  else
  {
    int vect,nrC;
    freadEndianSafe(&vect,1,sizeof(vect),modelFile);
    freadEndianSafe(&numberOfClusters,1,sizeof(numberOfClusters),modelFile);
    freadEndianSafe(&nrC,1,sizeof(nrC),modelFile);
    if(vect != vectorSize)
    {
      USER_ERROR("The clusterAM is not compatible! Wrong vector size...");
    }
    if(nrC != 1)
    {
      USER_ERROR("The clusterAM is not compatible! Number of 'clusters' in file should be one (no clusters within clusters)...");
    }
    clusterModels = new TrainPhoneModel*[numberOfClusters];
    for(int i=0;i<numberOfClusters;i++)
    {
      clusterModels[i] = new TrainPhoneModel(modelFile,vect);
    }            
    fclose(modelFile);
  }  
  
  strcpy(dataDir,adapt_dir);
  
  strcpy(str,dataDir);
  strcat(str,"/");
  strcat(str,PHONEFILE);
  
  FILE *phoneFile     = fopen(str, "rb");
  FILE *checkData     = NULL;
  if(phoneFile!=NULL)
  {    
    freadEndianSafe(&vectorSize,1,sizeof(vectorSize),phoneFile);
    freadEndianSafe(&numberOfModels,1,sizeof(numberOfModels),phoneFile);
    freadEndianSafe(&numberOfSil,1,sizeof(numberOfSil),phoneFile);
    PRINTF3("Number of models, number of SIL-models: %d %d\n",numberOfModels,numberOfSil);
    data = new DataStats[numberOfModels];
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
      if(data[i].nrOfSamples > 0)
      {
        checkData = fopen(str,"rb");
        if(checkData)
        {
          fclose(checkData);
        }
        else
        {
          fclose(phoneFile);
          delete[] data;
          data           = NULL;
          numberOfModels = 0;                  
          USER_ERROR("Data samples for one or more phones are missing!");
        }
      }
    }    
    fclose(phoneFile);
  }  

  PRINTF3("Checked the data directory. Seems okay. Number of models: %d with %d clusters.\n",numberOfModels,numberOfClusters);
    
  int numberOfClusteredModels = numberOfModels*numberOfClusters;  
  models = new TrainPhoneModel*[numberOfClusteredModels];

  FeaturePool *pools[numberOfClusteredModels];
  int          poolsTime[numberOfClusteredModels];
  int          maxPoolSize = TRAIN_MAX_FEATURE_POOL_SIZE/100;
  if(numberOfModels >= 5)
  {
    maxPoolSize = TRAIN_MAX_FEATURE_POOL_SIZE/(numberOfModels/5);
  }
  PRINTF("Initializing the adaptation data pools..\n");
  for(int i=0;i<numberOfClusteredModels;i++)
  {
    
    pools[i]     = new FeaturePool(maxPoolSize, vectorSize);
    poolsTime[i] = 0;
  }
  
  // Opening existing Acoustic Model:  
  modelFile = fopen(amName, "rb");
  if(modelFile!=NULL)
  {
    PRINTF("Loading AM file...\n");
    int nrM = 0;   
    int nrC = 0;
    int vectS = 0;
    freadEndianSafe(&vectS,1,sizeof(vectS),modelFile);
    assert(vectS == vectorSize);
    freadEndianSafe(&nrM,1,sizeof(nrM),modelFile);
    freadEndianSafe(&nrC,1,sizeof(nrC),modelFile);    
    if(nrC != 1)
    {
      USER_WARNING("The initial AM file contains clustered-models. Only the first model will be used for re-training!\n");
    }
    if(nrM != numberOfModels)
    {
      fclose(modelFile);
      USER_ERROR("The number of models in this file is not the same\nas the number of models in the training data!");
    }
    else
    {
      for(int i=0;i<numberOfModels;i++)
      {
        models[i] = new TrainPhoneModel(modelFile,vectorSize);
        models[i]->setTrainingData(pools[i],TRAIN_ID1,-1);
      }
      for(int j=1;j<numberOfClusters;j++)
      {
        for(int i=0;i<numberOfModels;i++)
        {
          models[i+j*numberOfModels] = new TrainPhoneModel(models[i],vectorSize);
          models[i+j*numberOfModels]->setTrainingData(pools[i+j*numberOfModels],TRAIN_ID1,-1);
        }
      }      
    }
    fclose(modelFile);
  } 
  
  PRINTF("Loading adaptation data...\n");
  for(int i=0;i<numberOfModels;i++)
  {
    strcpy(str,dataDir);
    strcat(str,"/");
    strncat(str,data[i].name,10);
    PRINTF2("  ...Opening %s\n",str);
    trainDataFile = fopen(str,"rb");
    if(trainDataFile != NULL)
    {
      int num;
      int notUsed;
      int contextLeft;
      int contextRight;
      int cdmID;
      while(!feof(trainDataFile))
      {   
        freadEndianSafe(&cdmID,1,sizeof(cdmID),trainDataFile);  // Reading the cluster-dependend-model ID: cdmID
        freadEndianSafe(&num,1,sizeof(int),trainDataFile);        
        freadEndianSafe(&notUsed,1,sizeof(notUsed),trainDataFile);             // Start of state 2
        freadEndianSafe(&notUsed,1,sizeof(notUsed),trainDataFile);             // Start of state 3
        freadEndianSafe(&contextLeft,1,sizeof(contextLeft),trainDataFile);        
        freadEndianSafe(&contextRight,1,sizeof(contextRight),trainDataFile);        
        if(!feof(trainDataFile))
        {
          int ID = i+(numberOfModels*cdmID);
          pools[ID]->addSegment(TRAIN_ID1, poolsTime[ID], poolsTime[ID]+num-1,0,contextLeft*numberOfModels+contextRight);
          pools[ID]->setCurSegmentVectors(TRAIN_ID1,trainDataFile);
          poolsTime[ID] += num;
        }
      }
      fclose(trainDataFile);
    }
    else
    {
      USER_ERROR("The training directory is empty or does not exist!\n");
    }
  }

  
  PRINTF("Calculating accumulators for the Gaussian components (training)...\n");
  for(int i=0;i<numberOfModels;i++)
  {
    PRINTF2("Training model %s Cluster: ",models[i]->getStatistics()->name);
    fflush(stdout);
    for(int j=0;j<numberOfClusters;j++)
    {
      PRINTF2("%02d  ",j);
      fflush(stdout);
      models[i+(numberOfModels*j)]->adapt_clear();
      models[i+(numberOfModels*j)]->adapt_setAcTrain();
    }
    PRINTF("\n");
  }

  if((numberOfModels == 1) && (numberOfSil == 1) && (strncmp(data[0].name,"UBM",3) == 0)) // this is a speaker model...
  {
    PRINTF("This is the UBM for SPEAKER RECOGNITION. I'm going to update the Gaussian means using MAP adaptation.\n"); 
    models[0]->mapAdaptMeans();
  }
  else
  {

    PRINTF("Creating structural tree for SMAPLR adaptation...\n");
    Adapt_AM_TreeNode **adaptationSilTree = new Adapt_AM_TreeNode*[numberOfClusters];
    Adapt_AM_TreeNode **adaptationTree    = new Adapt_AM_TreeNode*[numberOfClusters];
    for(int i=0;i<numberOfClusters;i++)
    {
      adaptationTree[i]                  = new Adapt_AM_TreeNode(vectorSize);
      adaptationTree[i]->meanGaussian    = new Gaussian(vectorSize);
      adaptationSilTree[i]               = new Adapt_AM_TreeNode(vectorSize);
      adaptationSilTree[i]->meanGaussian = new Gaussian(vectorSize);
    }    
    for(int j=0;j<numberOfClusters;j++)
    {
      for(int i=0;i<numberOfModels;i++)
      {
        if(models[i+j*numberOfModels]->isSilModel())
        {
          models[i+j*numberOfModels]->adapt_setInitialNode(adaptationSilTree[j]);
        }
        else
        {
          models[i+j*numberOfModels]->adapt_setInitialNode(adaptationTree[j]);
        }
      }
    }
  
    for(int j=0;j<numberOfClusters;j++)
    {
      int totNodes = 1;  
      int splitRes = adaptationTree[j]->split(1);
      while(splitRes > 0)
      {
        totNodes += splitRes;
        for(int i=0;i<numberOfModels;i++)
        {
          models[i+j*numberOfModels]->adapt_setNode();
        }
        splitRes = adaptationTree[j]->split(1);  
      }
      PRINTF3("Phone Tree %d created (%d nodes)\n",j,totNodes);
  
  
      PRINTF("Filling the helper matrices and distribute them through the tree\n");
      // Determine the first summations of Z and W... (going up)
      adaptationTree[j]->setHelperMatrices_1((PhoneModel**)(&models[j*numberOfModels]), numberOfModels);
      PRINTF("Calculating W for each node\n");
      adaptationTree[j]->setHelperMatrices_2();   // Determine the adaptation matrices... (going down)
  
      PRINTF("Adapting gaussians...\n");
      for(int i=0;i<numberOfModels;i++)
      {
        models[i+j*numberOfModels]->adapt_adapt();
      }
  
  
  
      PRINTF("Filling helper matrix for variance adaptation..\n");
      for(int i=0;i<numberOfModels;i++)
      {
        models[i+j*numberOfModels]->adapt_setVarTrans();
      }
      PRINTF("Adapting variances..\n");
      for(int i=0;i<numberOfModels;i++)
      {
        models[i+j*numberOfModels]->adapt_adaptVar();
      }
  
    }

    PRINTF("Done. Cleaning up adaptation matrices...\n"); fflush(stdout);
    for(int j=0;j<numberOfClusters;j++)
    {
      if(adaptationTree != NULL && adaptationTree[j] != NULL)
      {
        delete adaptationTree[j];
      }
      if(adaptationSilTree != NULL && adaptationSilTree[j] != NULL)
      {
        delete adaptationSilTree[j];
      }
    }
    if(adaptationTree != NULL)
    {
      delete[] adaptationTree;
    }
    if(adaptationSilTree != NULL)
    {
      delete[] adaptationSilTree;
    }

  }

  modelFile = fopen(newAmName, "wb");
  if(modelFile!=NULL)
  { 
    PRINTF("Writing adapted AM file to disc...\n");
    fwrite(&vectorSize,1,sizeof(vectorSize),modelFile);
    fwrite(&numberOfModels,1,sizeof(numberOfModels),modelFile);
    fwrite(&numberOfClusters,1,sizeof(numberOfClusters),modelFile);
    for(int i=0;i<numberOfClusteredModels;i++)
    {
      PRINTF2("  ... Writing model %02d ...\n",i);
      models[i]->writeModel(modelFile);
    }
    /*
    if(clusterModels != NULL)
    {
      for(int i=0;i<numberOfClusters;i++)
      {
        clusterModels[i]->writeModel(modelFile);
      }
    }
    */
    fclose(modelFile);
  } 
  PRINTF("Done. Have a nice day.\n");    
  for(int j=0;j<numberOfClusters;j++)
  {
    if(clusterModels != NULL && clusterModels[j] != NULL)
    {
      delete clusterModels[j];
    }
  }
  if(clusterModels != NULL)
  {
    delete[] clusterModels;
  }

  for(int i=0;i<numberOfClusteredModels;i++)
  {
    delete pools[i];
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor is empty.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Adapt_AM::~Adapt_AM()
{
}


