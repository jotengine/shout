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
#include <math.h>
#include "standard.h"
#include "vector.h"
#include "shout_maketrainset.h"
#include "phonefilereader.h"
#include "shout-misc.h"
#include "shout_preprocess.h"
#include "shoutconfig.h"

using namespace StringFunctions;
using namespace WriteFileLittleBigEndian;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Handles the command parameters for Shout_Preprocess. It is the main function of the 
/// shout_preprocess application. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_PREPROCESS, argc, argv);

  Shout_Preprocess Shout_Preprocess(myConfig.getStringValue(ARGUMENT_PHONE_LIST),
                                    myConfig.getStringValue(ARGUMENT_TRAINDIR),
                                    myConfig.getStringValue(ARGUMENT_AUDIO_LIST),
                                    myConfig.getStringValue(ARGUMENT_AM_PHONE)
                                   );
  return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

Shout_Preprocess::Shout_Preprocess(char *phoneList, char *traindir, char *audioList, char *amPhone)
{
  models          = NULL;

  doStateLabeling = (phoneList != NULL && strlen(phoneList) > 0 && amPhone != NULL && strlen(amPhone) > 0);

  char makeTrDir[MAXSTR];
  sprintf(makeTrDir,"mkdir %s",traindir);
  int resNotUsed = system(makeTrDir);

  if(doStateLabeling)
  {
    printf("=== Going to perform state labeling ===\n");
    printf("Reading phone list...\n");
    PhoneFileReader readPhones(phoneList, &data, &nrOfModels, &nrOfSil);
    printf("Reading acoustic model...\n");

    FILE *modelFile = fopen(amPhone, "rb");
    if(modelFile==NULL)
    {
      USER_ERROR("No valid AM file passed!");
    }
    int numberOfClusters = 0;
    int numberOfModels   = 0;
    freadEndianSafe(&vectorSize,1,sizeof(vectorSize),modelFile);
    freadEndianSafe(&numberOfModels,1,sizeof(numberOfModels),modelFile);
    freadEndianSafe(&numberOfClusters,1,sizeof(numberOfClusters),modelFile);
    assert(numberOfClusters == 1);
    if(numberOfModels != nrOfModels)
    {
      USER_ERROR("The number of phones in the phone-list do not match the number of phones in the AM!");
    }
    models = new PhoneModel*[nrOfModels];
    for(int i=0;i<nrOfModels;i++)
    {
      models[i] = new PhoneModel(modelFile);
    }
    fclose(modelFile);
  }
  else
  {
    printf("==!! NOT Going to perform state labeling !!==\n");
  }

  FeaturePool *featurePool = NULL;
  bool proceed             = true;
  FILE *listFile           = fopen(audioList,"rt");
  if(listFile == NULL)
  {
    USER_ERROR("Unable to read the task list");
  }
  FILE *inFile   = NULL;
  char fileID[MAXSTR];
  char rawName[MAXSTR];
  char rttmFileName[MAXSTR];
  char datFile[MAXSTR];
  char outFileName[MAXSTR];
  int  fileNr = 0;
  if(fscanf(listFile,"%s %s %s",rawName,rttmFileName,datFile) == 0)
  {
    USER_ERROR("Invalid task-list file!");
  }
  while(proceed)
  {
    int i=strlen(rawName)-1;
    while((i>=0) && (rawName[i] != '/'))
    {
      i--;
    }
    i++;
    strcpy(fileID,&rawName[i]);
    i=strlen(fileID)-1;
    while((i>=0) && (fileID[i] != '.'))
    {
      i--;
    }
    fileID[i] = 0;

    printf("Pre-processing '%s'\n",fileID); fflush(stdout);


    FILE *listF = fopen(rttmFileName, "rt");
    if(listF==NULL)
    {
      USER_ERROR("No valid RTTM file passed!");
    }

    featurePool = new FeaturePool(rawName,FILE_TYPE_RAW_DECODING,false,1.0,false,true);  // Only prepare...
    int segID = featurePool->claimSegmentationID();
    featurePool->readRTTM(segID,NULL,-1,-1,3,listF);
    fclose(listF);
    featurePool->setOnlyMEL();
    featurePool->createNewPool(rawName, NULL, NULL, false, false, segID);
    sprintf(outFileName,"%s/%s.melscale",traindir,fileID);
    FILE *outFile = fopen(outFileName,"wb");
    if(outFile == NULL)
    {
      USER_ERROR("Could not write the feature file in the training directory.");
    }
    featurePool->storeFeatureVectors(outFile);
    fclose(outFile);
    int poolLength = featurePool->getPoolLength();
    delete featurePool;
    featurePool = NULL;

    if(doStateLabeling)
    {
      inFile       = fopen(datFile, "rt");
      if(inFile == NULL)
      {
        USER_ERROR("Invalid hypothesis file");
      }
      int label[poolLength];
      for(int i=0;i<poolLength;i++)
      {
        label[i] = -1;
      }
      int a,b,c,d,e,f,i;
      int state           = 0;
      int contextLeft     = 0;
      int contextRight    = 0;
      char strPhone[200];
      char strPhone2[200];
      char strPhone3[200];  
      char strID[MAXSTR];
      char str[MAXSTR];
      while(!feof(inFile))
      {
        // First look what we have: label-definition, phone, or other crap:
        if(fscanf(inFile,"###################################################-| %s |-#################################",strID) != 0)
        {
        }
        else if(fscanf(inFile,"          #Phone# %s %s %s | %d %d | %d %d | %d %d |",strPhone,strPhone2,strPhone3,&a,&b,&c,&d,&e,&f) != 0)
        {
          contextLeft=0;
          i=0;
          contextRight=0;
          strcat(strPhone,"          ");
          while((strncmp(strPhone,data[contextLeft].name,10)!=0) && (contextLeft<nrOfModels))
          {
            contextLeft++;
          }
          if(contextLeft>=nrOfModels)
          {
            USER_ERROR("A phone (left context) was found in the HYP file, but does not appear in the phone list!");
          }

          strcat(strPhone2,"          ");
          while((strncmp(strPhone2,data[i].name,10)!=0) && (i<nrOfModels))
          {
            i++;
          }
          if(i>=nrOfModels)
          {
            USER_ERROR("A phone (middle) was found in the HYP file, but does not appear in the phone list!");
          }

          strcat(strPhone3,"          ");
          while((strncmp(strPhone3,data[contextRight].name,10)!=0) && (contextRight<nrOfModels))
          {
            contextRight++;
          }
          if(contextRight>=nrOfModels)
          {
            USER_ERROR("A phone (right context) was found in the HYP file, but does not appear in the phone list!");
          }

          assert((i < nrOfModels) && (i >= 0));
          int contextKey = contextLeft*nrOfModels + contextRight;

          if(a>=0)      // LEFT STATE
          {
            int stateNr = i*MAX_NUMBER_OF_MIXTURES_PER_PHONE + models[i]->getStateNr(contextKey,0);
            for(int jj=a;jj<=b;jj++)
            {
              label[jj] = stateNr;
            }
          }
          if(c>=0)      // MIDDLE STATE
          {
            int stateNr = i*MAX_NUMBER_OF_MIXTURES_PER_PHONE + models[i]->getStateNr(contextKey,1);
            for(int jj=c;jj<=d;jj++)
            {
              label[jj] = stateNr;
            }
          }
          if(e>=0)      // RIGHT STATE
          {
            int stateNr = i*MAX_NUMBER_OF_MIXTURES_PER_PHONE + models[i]->getStateNr(contextKey,2);
            for(int jj=e;jj<=f;jj++)
            {
              label[jj] = stateNr;
            }
          }
        }
        else
        {
          char *resChar = fgets(str,MAXSTR,inFile);
          assert(resChar != NULL);
        }
      }
      fclose(inFile);
      inFile = NULL;

      sprintf(outFileName,"%s/%s.labels",traindir,fileID);
      outFile = fopen(outFileName,"wt");
      if(outFile == NULL)
      {
        USER_ERROR("Could not write the label file in the training directory.");
      }
      for(int ii=0;ii<poolLength;ii++)
      {
        fprintf(outFile,"%d\t%d\n",ii,label[ii]);
      }
      fclose(outFile);
    }
    proceed = !(listFile == NULL || feof(listFile));
    strcpy(rawName,"");
    if(proceed)
    {
      if(fscanf(listFile,"%s %s %s",rawName,rttmFileName,datFile) == 0)
      {
        USER_ERROR("Invalid train-list file!");
      }
      proceed = (strlen(rawName) > 0);
    }
  }
  if(listFile != NULL)
  {
    fclose(listFile);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

Shout_Preprocess::~Shout_Preprocess()
{
  if(models != NULL)
  {
    for(int i=0;i<nrOfModels;i++)
    {
      delete models[i];
    }
    delete[] models;
    models = NULL;
  }
}


