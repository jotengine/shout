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
#include "shouttrainfinishsat.h"

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
  ShoutConfig myConfig(APPID_TRAIN_FINISH_SAT, argc, argv);

  ShoutTrainFinishSAT myModel(myConfig.getStringValue(ARGUMENT_TRAINDIR),
                              myConfig.getStringValue(ARGUMENT_SATLIST),
                              myConfig.getStringValue(ARGUMENT_AM_CREATE),
                              (int)(myConfig.getFloatValue(ARGUMENT_TRAINMODELNR)));
  return 0;
}

ShoutTrainFinishSAT::ShoutTrainFinishSAT(char *traindir, char *satlist, char *am, int modelNumber)
{
  int              numberOfModels    = 1;
  int              numberOfSil       = 1;
  int              numberOfClusters  = 1;
  int              vectorSize        = 0;
  DataStats       *data              = NULL;
  TrainPhoneModel *model             = NULL;
  char str[MAXSTR];

  strcpy(str,traindir);
  strcat(str,"/");
  strcat(str,PHONEFILE);

  FILE *checkData     = NULL;
  FILE *trainDataFile = NULL;
  FILE *phoneFile     = fopen64(str, "rb");
  if(phoneFile==NULL)
  {
    USER_ERROR("The training directory is not valid! (no phones.lst)");
  }
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
    strcpy(str,traindir);
    strcat(str,"/");
    strncat(str,data[i].name,10);

    if(i == modelNumber)
    {
      printf("Found the model to train: %s\n",name);
    }
  }
  fclose(phoneFile);
  phoneFile = NULL;

  // Opening existing model:

  FILE *modelFile = fopen(am, "rb");
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
        if(model != NULL)
        {
          delete model;
          model = NULL;
        }
        model = new TrainPhoneModel(modelFile,vectorSize);  // Reading the phone model
      }
    }
    fclose(modelFile);
    modelFile = NULL;
  }
  else  // Create model
  {
    USER_ERROR("Could not read the input acoustic model!");
  }

  assert(model != NULL);

  if(model->isSilModel())
  {
    printf("This is a SIL model. Just copying it from the original AM...\n");
  }
  else
  {
    FILE *listFile = fopen(satlist, "rt");
    while(listFile!=NULL && !feof(listFile))
    {
      strcpy(str,"");
      int resNotUsed = fscanf(listFile,"%s\n",str);
      printf("Appending: %s\n",str);
      strcat(str,"/");
      strncat(str,data[modelNumber].name,10);
      strcat(str,".sat");
      FILE *satFile = fopen(str,"rb");
      if(satFile == NULL)
      {
        USER_WARNING2("Skipping file",str);
      }
      else
      {
        model->appendSAT(satFile);
        fclose(satFile);
      }
    }
    fclose(listFile);
  
    printf("Finishing training... Total training data count: %f\n",model->finishSAT());
  }

  printf("Done! Writing model to file...\n");
  strcpy(str,traindir);
  strcat(str,"/");
  strncat(str,data[modelNumber].name,10);
  strcat(str,".am");
  FILE *outFile = fopen(str,"wb");
  if(outFile == NULL)
  {
    USER_ERROR("Could not write in the training directory.");
  }

  model->writeModel(outFile);
  fclose(outFile);


  if(model != NULL)
  {
    delete model;
  }
  if(data != NULL)
  {
    delete[] data;
  }
}


ShoutTrainFinishSAT::~ShoutTrainFinishSAT()
{
}


