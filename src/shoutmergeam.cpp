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

#include <iostream>
#include "standard.h"
#include "shoutmergeam.h"
#include "trainphonemodel.h"
#include "shoutconfig.h"

#include <stdlib.h>

#include "shout-misc.h"

using namespace WriteFileLittleBigEndian;


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Handles the command parameters for ShoutMergeAm. It is the main function of the 
/// shout_merge_am application. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_MERGE_AM, argc, argv);
  ShoutMergeAm ShoutMergeAm(myConfig.getStringValue(ARGUMENT_AM_MERGE1),myConfig.getStringValue(ARGUMENT_AM_MERGE2),myConfig.getStringValue(ARGUMENT_AM_OUT));
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add docs.
/////////////////////////////////////////////////////////////////////////////////////////////////////

ShoutMergeAm::ShoutMergeAm(char *amName1, char *amName2, char *amOutName)
{
  models1 = NULL;
  models2 = NULL;
  FILE *modelFile = fopen(amName1, "rb");
  bool classificationModel = false;
  int  vectorSize          = ASR_DEFAULT_VECTORSIZE;
  if(modelFile!=NULL)
  { 
    int nrM = 0;
    int numberOfClusters = 0;
    freadEndianSafe(&vectorSize,1,sizeof(vectorSize),modelFile);
    freadEndianSafe(&nrM,1,sizeof(nrM),modelFile);
    freadEndianSafe(&numberOfClusters,1,sizeof(numberOfClusters),modelFile);
    if(numberOfClusters < 0) // This is a classification model! (no delta's and delta-delta's)
    {
      classificationModel = true;
    }
    numberOfModels1  = nrM;
    models1          = new TrainPhoneModel*[numberOfModels1];
    for(int i=0;i<numberOfModels1;i++)
    {      
      models1[i]           = new TrainPhoneModel(modelFile,vectorSize);
    }
    fclose(modelFile);
  }
  else
  {
    USER_ERROR("The first Acoustic Model file does not exist!");
  }
  modelFile = fopen(amName2, "rb");
  if(modelFile!=NULL)
  { 
    int nrM = 0;
    int numberOfClusters = 0;
    int vectS = ASR_DEFAULT_VECTORSIZE;
    freadEndianSafe(&vectS,1,sizeof(vectS),modelFile);
    if(vectorSize != vectS)
    {
      USER_ERROR("The two files contain vectors of different sizes!");
    }
    freadEndianSafe(&nrM,1,sizeof(nrM),modelFile);
    freadEndianSafe(&numberOfClusters,1,sizeof(numberOfClusters),modelFile);
    if((classificationModel && (numberOfClusters >= 0)) || (!classificationModel && (numberOfClusters < 0)))
    {
      USER_ERROR("The models are not compatible! Pick either classification models or ASR models!");
    }
    numberOfModels2  = nrM;
    models2          = new TrainPhoneModel*[numberOfModels2];
    for(int i=0;i<numberOfModels2;i++)
    {
      models2[i]           = new TrainPhoneModel(modelFile,vectorSize);
    }
    fclose(modelFile);
  }
  else
  {
    USER_ERROR("The second Acoustic Model file does not exist!");
  } 
  
  int numberOfModelsSmallest = numberOfModels1;
  if(numberOfModelsSmallest > numberOfModels2)
  {
    numberOfModelsSmallest = numberOfModels2;
  }
  
  printf("Merging the two AM files can be done using one or more of the following three methods:\n");
  printf("1) For the first %d models, pick the corresponding model from one of the two files\n", numberOfModelsSmallest);
  printf("2) Pick one or more models from the first file to append to the resulting file.\n");
  printf("3) Pick one or more models from the second file to append to the resulting file.\n");
  printf("\n");
  printf("File A contains %d phones.\n",numberOfModels1);  
  printf("File B contains %d phones.\n",numberOfModels2);
  printf("How many phones should the resulting file contain? ");
  fflush(stdout);
  
  char temp[100];
  std::cin.getline(temp, 100);
  numberOfModelsTot = strtol(temp,0,10);  
  
  int numberOfClusters = 1;
  if(classificationModel)
  {
    numberOfClusters = -1;
  }
  // Writing models:  
  modelFile = fopen(amOutName, "wb");
  if(modelFile==NULL)
  { 
    USER_ERROR("Could not write to output file!\n");
  }

  fwrite(&vectorSize,1,sizeof(vectorSize),modelFile);
  fwrite(&numberOfModelsTot,1,sizeof(numberOfModelsTot),modelFile);
  fwrite(&numberOfClusters,1,sizeof(numberOfClusters),modelFile);

  int modelsUsed = 0;
  
  printf("Do you want to use merging method 1 (pick from file A or B)? (y/n) ");
  fflush(stdout);
  int res = std::cin.get();
  while(res != 'y' && res != 'Y' && res != 'n' && res != 'N')
  { 
    res = std::cin.get();
  }
  if(res == 'y' || res == 'Y')
  {
    int i=0;
    while(i<numberOfModelsSmallest && modelsUsed < numberOfModelsTot)
    {
      printf("Do you want to use model A or model B for phone (%s - %s): (A/B)",models1[i]->getStatistics()->name
                                                                               ,models2[i]->getStatistics()->name);
      fflush(stdout);
      int res = std::cin.get();
      while(res != 'a' && res != 'A' && res != 'b' && res != 'B')
      { 
        res = std::cin.get();
      }
      if(res == 'a' || res == 'A')
      {
        models1[i]->writeModel(modelFile);
      }
      else
      {
        models2[i]->writeModel(modelFile);
      }
      i++;
      modelsUsed++;
    }
  } 
  
  
  if(modelsUsed < numberOfModelsTot)
  {
    printf("Do you want to use merging method 2 (add from file A)? (y/n) ");
    fflush(stdout);
    int res = std::cin.get();
    while(res != 'y' && res != 'Y' && res != 'n' && res != 'N')
    { 
      res = std::cin.get();
    }
    if(res == 'y' || res == 'Y')
    {
      int i=0;
      while(i<numberOfModels1 && modelsUsed < numberOfModelsTot)
      {
        printf("Do you want to use phone %s from file A: (y/n)",models1[i]->getStatistics()->name);
        fflush(stdout);
        int res = std::cin.get();
        while(res != 'y' && res != 'Y' && res != 'n' && res != 'N')
        { 
          res = std::cin.get();
        }
        if(res == 'y' || res == 'Y')
        {
          models1[i]->writeModel(modelFile);
          modelsUsed++;
        }
        i++;
      }
    }
  }
  
  if(modelsUsed < numberOfModelsTot)
  {
    printf("Do you want to use merging method 3 (add from file B)? (y/n) ");
    fflush(stdout);
    int res = std::cin.get();
    while(res != 'y' && res != 'Y' && res != 'n' && res != 'N')
    { 
      res = std::cin.get();
    }
    if(res == 'y' || res == 'Y')
    {
      int i=0;
      while(i<numberOfModels2 && modelsUsed < numberOfModelsTot)
      {
        printf("Do you want to use phone %s from file B: (y/n)",models2[i]->getStatistics()->name);
        fflush(stdout);
        int res = std::cin.get();
        while(res != 'y' && res != 'Y' && res != 'n' && res != 'N')
        { 
          res = std::cin.get();
        }
        if(res == 'y' || res == 'Y')
        {
          models2[i]->writeModel(modelFile);
          modelsUsed++;
        }
        i++;
      }
    }
  }
     
  fclose(modelFile);
  if(modelsUsed != numberOfModelsTot)
  {
    USER_ERROR("THE RESULT FILE IS INVALID: THE NUMBER OF USED MODELS DOES NOT MATCH THE INITIAL GIVEN NUMBER!\n");
  }
  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Add docs.
/////////////////////////////////////////////////////////////////////////////////////////////////////

ShoutMergeAm::~ShoutMergeAm()
{
  if(models1 != NULL)
  {
    for(int i=0;i<numberOfModels1;i++)
    {
      delete models1[i];
    }
    delete[] models1;
  }
  if(models2 != NULL)
  {
    for(int i=0;i<numberOfModels2;i++)
    {
      delete models2[i];
    }
    delete[] models2;
  }
}


