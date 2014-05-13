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
#include "mergetrainset.h"
#include "shoutconfig.h"
#include "shout-misc.h"

using namespace WriteFileLittleBigEndian;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_MERGETRAINSET, argc, argv);

  MergeTrainSet MergeTrainSet(myConfig.getStringValue(ARGUMENT_AM_TYPE),
                              myConfig.getStringValue(ARGUMENT_TRAINDIR),
                              myConfig.getStringValue(ARGUMENT_TRAINDIR2),
                              myConfig.getStringValue(ARGUMENT_PHONE_LIST));
  return 0;
}


MergeTrainSet::MergeTrainSet(char *task, char *traindir1, char *traindir2, char *phoneList) : Shout_MakeTrainSet(task,phoneList,traindir1)
{
  char str[MAXSTR];
  strcpy(str,traindir2);
  strcat(str,"/");
  strcat(str,PHONEFILE);

  FILE *checkData     = NULL;
  FILE *phoneInfo     = fopen64(str, "rb");
  if(phoneInfo==NULL)
  {
    USER_ERROR("The second training directory is not valid! (no phones.lst)");
  }
  int vS, nM, nS;
  freadEndianSafe(&vS,1,sizeof(vS),phoneInfo);
  freadEndianSafe(&nM,1,sizeof(nM),phoneInfo);
  freadEndianSafe(&nS,1,sizeof(nS),phoneInfo);

  printf("Second training directory: (Vector size, number of models, number of SIL-models) = (%d,%d,%d)\n",vS,nM,nS);

  if(vS != vectorSize)
  {
    USER_ERROR("The vector size of the features in the first and second training directory do not match.");
  }
  if(strcmp(task,"HIST") != 0 && nM != nrOfModels)
  {
    USER_ERROR("The number of models in the first and second training directory do not match.");
  }
  if(nS != nrOfSil)
  {
    USER_ERROR("The number of SIL models in the first and second training directory do not match.");
  }

  for(int i=0;i<nrOfModels;i++)
  {
    DataStats d;
    char name[11];
    freadEndianSafe(name,10,1,phoneInfo);
    name[10] = 0;
    if(strcmp(data[i].name,name) != 0)
    {
      printf("Traindir1: %s\nTraindir2: %s\n",data[i].name,name);
      USER_ERROR("The names of the phones of the two training directories do not match.");
    }
    freadEndianSafe(d.nrOfSamples,MAX_CLUSTERS,sizeof(int),phoneInfo);
    strcpy(str,traindir2);
    strcat(str,"/");
    strncat(str,data[i].name,10);

    printf("Processing: (%d - %d) %s\n",data[i].nrOfSamples[0],d.nrOfSamples[0], data[i].name);
    if(d.nrOfSamples[0] > 0)
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

      while(!feof(checkData))
      {
        numberOfObservations = 0;
        freadEndianSafe(&cdmID,1,sizeof(cdmID),checkData);
        freadEndianSafe(&numberOfObservations,1,sizeof(numberOfObservations),checkData);
        freadEndianSafe(&startOfState2,1,sizeof(startOfState2),checkData);
        freadEndianSafe(&startOfState3,1,sizeof(startOfState3),checkData);
        freadEndianSafe(&contextLeft,1,sizeof(contextLeft),checkData);
        freadEndianSafe(&contextRight,1,sizeof(contextRight),checkData);

        if(!feof(checkData) && numberOfObservations > 0)
        {
          totalPhonesStored++;
          fwrite(&cdmID,1,sizeof(cdmID),phoneFile[i]);
          fwrite(&numberOfObservations,1,sizeof(numberOfObservations),phoneFile[i]);
          fwrite(&startOfState2,1,sizeof(startOfState2),phoneFile[i]);
          fwrite(&startOfState3,1,sizeof(startOfState3),phoneFile[i]);
          fwrite(&contextLeft,1,sizeof(contextLeft),phoneFile[i]);
          fwrite(&contextRight,1,sizeof(contextRight),phoneFile[i]);
          if(performWhat != TASK_TYPE_ARTICULATORY)
          {
            data[i].nrOfSamples[cdmID]++;
          }
          for(int b=0;b<numberOfObservations;b++)
          {
            Vector *v = new Vector(checkData,vectorSize);
            if(v==NULL)
            {
              USER_ERROR("Not enough vectors in the feature file..");
            }
            if(isnan(v->getValue(1)))
            {
              USER_ERROR("One of the observation values is NotANumber!!");
            }
            v->storeData(phoneFile[i]);
            delete v;
          }
          fflush(phoneFile[i]);
        }
      }
      fclose(checkData);
      checkData = NULL;
    }
  }
  fclose(phoneInfo);
  phoneInfo = NULL;
}


MergeTrainSet::~MergeTrainSet()
{
}


