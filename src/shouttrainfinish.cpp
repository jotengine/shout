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
#include "shouttrainfinish.h"
#include "shout-misc.h"

using namespace WriteFileLittleBigEndian;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Handles the command parameters for ShoutTrainFinish. It is the main function of the 
/// shout_train_finish application. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_TRAIN_FINISH, argc, argv);

  ShoutTrainFinish myFinish(myConfig.getStringValue(ARGUMENT_TRAINDIR), myConfig.getStringValue(ARGUMENT_AM_CREATE));
  return 0;
}


ShoutTrainFinish::ShoutTrainFinish(char *trainDir, char *createAM)
{
  int              numberOfModels    = 1;
  int              numberOfSil       = 1;
  int              numberOfClusters  = 1;
  int              vectorSize        = 0;
  DataStats       *data              = NULL;
  TrainPhoneModel *model             = NULL;
  char str[MAXSTR];

  strcpy(str,trainDir);
  strcat(str,"/");
  strcat(str,PHONEFILE);

  FILE *phoneFile     = fopen64(str, "rb");
  if(phoneFile==NULL)
  {
    USER_WARNING("The training directory is not complete! (no phones.lst). I assume you are training an SRE model!");
    vectorSize     = 26;
    numberOfModels =  1;
    numberOfSil    =  1;
    data = new DataStats[numberOfModels];
    strcpy(data[0].name,"UBM");
    data[0].nrOfSamples[0] = 1;
  }
  else
  {
    freadEndianSafe(&vectorSize,1,sizeof(vectorSize),phoneFile);
    freadEndianSafe(&numberOfModels,1,sizeof(numberOfModels),phoneFile);
    freadEndianSafe(&numberOfSil,1,sizeof(numberOfSil),phoneFile);
    data = new DataStats[numberOfModels];
  
    // We will check all models, although we are only going to train one. Don't want to find out something is wrong at the final model..
    for(int i=0;i<numberOfModels;i++)
    {
      char name[11];
      freadEndianSafe(name,10,1,phoneFile);
      name[10] = 0;
      strcpy(data[i].name,name);
      freadEndianSafe(data[i].nrOfSamples,MAX_CLUSTERS,sizeof(int),phoneFile);
    }
    fclose(phoneFile);
    phoneFile = NULL;
  }

  // Opening our data files:
  FILE *modelFile = fopen(createAM, "wb");
  if(modelFile==NULL)
  {
    USER_ERROR("Could not write to, or create the AM");
  }

  fwrite(&vectorSize,1,sizeof(vectorSize),modelFile);
  fwrite(&numberOfModels,1,sizeof(numberOfModels),modelFile);
  fwrite(&numberOfClusters,1,sizeof(numberOfClusters),modelFile);
//  int numberOfClusteredModels = numberOfModels*numberOfClusters;
  for(int i=0;i<numberOfModels;i++)
  {
    strcpy(str,trainDir);
    strcat(str,"/");
    strncat(str,data[i].name,10);
    strcat(str,".am");
    FILE *outFile = fopen(str,"rb");
    if(outFile == NULL)
    {
      USER_ERROR("Could not read one of the acoustic models.");
    }
    if(model != NULL)
    {
      delete model;
      model = NULL;
    }
    model = new TrainPhoneModel(outFile,vectorSize);  // Reading the phone model
    fclose(outFile);
    model->writeModel(modelFile);
  }
  // Check if there is a decisiontree file. Append if so..
  strcpy(str,trainDir);
  strcat(str,ARTICULATORY_FILENAME);
  FILE *dataFile = fopen(str, "rb");
  if(dataFile != NULL)
  {
    int decisionMatrix[MAX_NUMBER_OF_DECISION_RULES*numberOfModels];
    int artModelNumber[MAX_NUMBER_OF_DECISION_RULES];
    int nrOfModels;
    freadEndianSafe(&nrOfModels,1,sizeof(nrOfModels),dataFile);
    fwrite(&nrOfModels,1,sizeof(nrOfModels),modelFile);
    for(int i=0;i<nrOfModels;i++)
    {
      freadEndianSafe(&(decisionMatrix[i*MAX_NUMBER_OF_DECISION_RULES]),MAX_NUMBER_OF_DECISION_RULES,sizeof(int),dataFile);
      fwrite(&(decisionMatrix[i*MAX_NUMBER_OF_DECISION_RULES]),MAX_NUMBER_OF_DECISION_RULES,sizeof(int),modelFile);
    }
    freadEndianSafe(&(artModelNumber[0]),MAX_NUMBER_OF_DECISION_RULES,sizeof(int),dataFile);
    fwrite(&(artModelNumber[0]),MAX_NUMBER_OF_DECISION_RULES,sizeof(int),modelFile);
    fclose(dataFile);
  }
  fclose(modelFile);
}


ShoutTrainFinish::~ShoutTrainFinish()
{
}


