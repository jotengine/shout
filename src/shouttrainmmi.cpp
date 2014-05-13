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
#include "shout-misc.h"
#include "shouttrainmmi.h"

using namespace WriteFileLittleBigEndian;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_TRAIN_MMI, argc, argv);

  ShoutTrainMMI   myModel(myConfig.getStringValue(ARGUMENT_AM_PHONE),
                          myConfig.getStringValue(ARGUMENT_AM_OUT),
                          myConfig.getStringValue(ARGUMENT_LATTICETRAINING_ENUM),
                          myConfig.getStringValue(ARGUMENT_LATTICETRAINING_DENOM));
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

ShoutTrainMMI::ShoutTrainMMI(char *amIn, char *amOut, char *enumFile, char *denomFile)
{
  int numberOfModels    = 1;
  int numberOfSil       = 1;
  int numberOfClusters  = 1;
  int vectorSize        = 0;
  FILE *modelFile = fopen(amIn, "rb");
  if(modelFile==NULL)
  {
    USER_ERROR("Could not read the input acoustic model");
  }

  freadEndianSafe(&vectorSize,1,sizeof(vectorSize),modelFile);
  freadEndianSafe(&numberOfModels,1,sizeof(numberOfModels),modelFile);
  freadEndianSafe(&numberOfClusters,1,sizeof(numberOfClusters),modelFile);
  if(numberOfClusters != 1)
  {
    printf("WARNING: The acoustic model file contains clustered-models (%d clusters).\n",numberOfClusters);
    printf("         All data will be used to train 1 set of models!\n");
  }
  numberOfClusters = 1;

  TrainPhoneModel *model[numberOfModels];
  for(int i=0;i<numberOfModels;i++)
  {
    model[i] = new TrainPhoneModel(modelFile,vectorSize);  // Reading the phone model
  }
  fclose(modelFile);
  modelFile = NULL;


  printf("Reading statistics files...\n");
  FILE *mmiTrainFileEnum  = fopen(enumFile,"rb");
  if(mmiTrainFileEnum == NULL)
  {
    USER_ERROR("Could not read the enumurator file...");
  }
  FILE *mmiTrainFileDenom = fopen(denomFile,"rb");
  if(mmiTrainFileDenom == NULL)
  {
    USER_ERROR("Could not read the denomonator file...");
  }
  for(int i=0;i<numberOfModels;i++)
  {
    printf("Training phone %d...\n",i); fflush(stdout);
    model[i]->trainMMI(mmiTrainFileEnum,mmiTrainFileDenom);
  }
  printf("Finished training. Let's write the output AM to file...\n"); fflush(stdout);
  fclose(mmiTrainFileEnum);
  fclose(mmiTrainFileDenom);




  // Done, now write to file..
  modelFile = fopen(amOut, "wb");
  if(modelFile==NULL)
  {
    USER_ERROR("Could not write to, or create the AM");
  }
  fwrite(&vectorSize,1,sizeof(vectorSize),modelFile);
  fwrite(&numberOfModels,1,sizeof(numberOfModels),modelFile);
  fwrite(&numberOfClusters,1,sizeof(numberOfClusters),modelFile);
  for(int i=0;i<numberOfModels;i++)
  {
    model[i]->writeModel(modelFile);
    delete model[i];
  }
  fclose(modelFile);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

ShoutTrainMMI::~ShoutTrainMMI()
{
}


