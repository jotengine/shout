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
#include "speakerrecognition.h"
#include "shoutconfig.h"
#include "trainphonemodel.h"
#include "shout-misc.h"
#include "featurepool.h"

using namespace WriteFileLittleBigEndian;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Handles the command parameters for SpeakerRecognition. It is the main function of the 
/// shout_spkrec application. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_SPKREC, argc, argv);

  SpeakerRecognition SpeakerRecognition(myConfig.getStringValue(ARGUMENT_SPKREC_UBM_MALE),
                                        myConfig.getStringValue(ARGUMENT_SPKREC_UBM_FEMALE),
                                        myConfig.getStringValue(ARGUMENT_SPKREC_TASKLIST),
                                        myConfig.getStringValue(ARGUMENT_SPKREC_MODELDIR),
                                        myConfig.getStringValue(ARGUMENT_SPKREC_OUTPUT_PREFIX),
                                        myConfig.getStringValue(ARGUMENT_SPKREC_ZNORM_LIST),
                                        myConfig.parIsDefined(ARGUMENT_SPKREC_OUTPUT_STATS));
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

SpeakerRecognition::SpeakerRecognition(char *ubmMaleFile, char *ubmFemaleFile, char *tasklist, 
                                       char *modelDir, char *outputPrefix, char *znormList, bool outputS)
{
  outputStats = outputS;
  for(int i=0;i<MAX_UBMSIZE;i++)
  {
    ubmMale[i]   = NULL;
    ubmFemale[i] = NULL;
  }
  vectorSize     = 0;
  zNormAvailable = false;
  zNormMean      = 0.0;
  zNormVariance  = 1.0;
  nrUbmMale      = 0;
  nrUbmFemale    = 0;

  readUbmModel(ubmMaleFile,true);
  readUbmModel(ubmFemaleFile,false);

  if(outputStats)
  {
    FILE *outMeans     = fopen("m.ubm_means","wt");
    FILE *outVariances = fopen("m.ubm_variances","wt");
    FILE *outWeights   = fopen("m.ubm_weights","wt");
    ubmMale[0]->printModel(outMeans,outVariances,outWeights);
    fclose(outMeans);
    fclose(outVariances);
    fclose(outWeights);
    outMeans           = fopen("f.ubm_means","wt");
    outVariances       = fopen("f.ubm_variances","wt");
    outWeights         = fopen("f.ubm_weights","wt");
    ubmFemale[0]->printModel(outMeans,outVariances,outWeights);
    fclose(outMeans);
    fclose(outVariances);
    fclose(outWeights);
  }

  if((znormList != NULL) && (strlen(znormList)>0))  // Performing Z-Normalization
  {
    FILE *zNormFile = fopen(znormList,"rt");
    char zNormSpeaker[MAXSTR];
    char gender[MAXSTR];
    if(zNormFile != NULL)
    {
      while(!feof(zNormFile))
      {
        strcpy(zNormSpeaker,"");
        int resNotUsed = fscanf(zNormFile,"%s %s",zNormSpeaker, gender);
        if((strlen(zNormSpeaker) > 0))
        {
          fprintf(stderr,"========== Starting run for speaker %s ==============\n",zNormSpeaker);

          char speakerModelName[MAXSTR];
          strcpy(speakerModelName, modelDir);
          if(strcmp(gender,"m") == 0)
          {
            strcat(speakerModelName, "male-");
          }
          else if(strcmp(gender,"f") == 0)
          {
            strcat(speakerModelName, "female-");
          }
          else
          {
            USER_ERROR("Could not determine gender in the z-norm list file.");
          }
          strcat(speakerModelName,zNormSpeaker);
          strcat(speakerModelName,".am");
          PhoneModel *speakerModel = readModel(speakerModelName,false);

          Gaussian *znorm = runTrials(tasklist, modelDir, outputPrefix, speakerModel);
          assert(znorm != NULL);

          zNormAvailable = true;
          zNormMean      = znorm->getMean()->getValue(0);
          zNormVariance  = sqrt(znorm->getVariance()->getValue(0));
          fprintf(stderr,"Distribution for speaker %s: mean %f, variance %f\n",zNormSpeaker, zNormMean, zNormVariance);
          writeModel(speakerModelName, speakerModel);
          delete znorm;
          delete speakerModel;
          zNormAvailable = false;
          zNormMean      = 0.0;
          zNormVariance  = 1.0;
        }
      }
      fclose(zNormFile);
    }
    else
    {
      USER_ERROR("Could not read the Z-Norm list file");
    }
  }
  else
  {
    Gaussian *znorm = runTrials(tasklist, modelDir, outputPrefix);
    if(znorm != NULL) // in case we are running without z-normation, a new one will be calculated automatically. Remove it...
    {
      delete znorm;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

SpeakerRecognition::~SpeakerRecognition()
{
 for(int i=0;i<MAX_UBMSIZE;i++)
  {
    if(ubmMale[i] != NULL)
    {
      delete ubmMale[i];
      ubmMale[i] = NULL;
    }
    if(ubmFemale[i] != NULL)
    {
      delete ubmFemale[i];
      ubmFemale[i] = NULL;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

PhoneModel *SpeakerRecognition::readModel(char *modelName, bool applyZNorm)
{
  PhoneModel *res = NULL;
  FILE *modelFile = fopen(modelName, "rb");
  if(modelFile!=NULL)
  {
    int nrM, nrC, vectS;
    freadEndianSafe(&vectS,1,sizeof(vectS),modelFile);
    if(vectorSize != vectS)
    {
      USER_ERROR("This model is not compatible!");
    }
    freadEndianSafe(&nrM,1,sizeof(nrM),modelFile);
    freadEndianSafe(&nrC,1,sizeof(nrC),modelFile);
    if(nrM != 1)
    {
      fclose(modelFile);
      USER_ERROR("The file is not a model file.");
    }
    res = new PhoneModel(modelFile,vectorSize);
    if(applyZNorm)
    {
      if(freadEndianSafe(&zNormMean,1,sizeof(zNormMean),modelFile) == sizeof(zNormMean))
      {
        if(freadEndianSafe(&zNormVariance,1,sizeof(zNormVariance),modelFile) == sizeof(zNormVariance))
        {
//          zNormVariance = 1.0;
          zNormAvailable = true;
        }
      }
    }
    else
    {
      zNormAvailable = false;
      zNormMean      = 0.0;
      zNormVariance  = 1.0;
    }
    fclose(modelFile);
  }
  else
  {
    USER_ERROR("Could not read the model file.");
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void SpeakerRecognition::readUbmModel(char *modelName, bool male)
{
  int mf = 0;
  if(male)
  {
    mf = 1;
  }
  PhoneModel *res = NULL;
  FILE *modelFile = fopen(modelName, "rb");
  int nr = 0;
  PhoneModel **ubm = ubmFemale;
  if(male)
  {
    ubm = ubmMale;
  }
  if(modelFile == NULL)
  {
    USER_ERROR("Could not read the model file.");
  }

  while(!feof(modelFile))
  {
    int nrM, nrC, vectS;
    int r = freadEndianSafe(&vectS,1,sizeof(vectS),modelFile);
    if(r == sizeof(vectS))
    {
      if(nr == 0)
      {
        vectorSize = vectS;
      }
      if(vectorSize != vectS)
      {
        USER_ERROR("This model is not compatible!");
      }
      freadEndianSafe(&nrM,1,sizeof(nrM),modelFile);
      freadEndianSafe(&nrC,1,sizeof(nrC),modelFile);
      if(nrM != 1)
      {
        fclose(modelFile);
        USER_ERROR("The file is not a model file.");
      }
      ubm[nr] = new PhoneModel(modelFile,vectorSize);
      if(nr != 0)
      {
        ubmZNormAvailable[mf] = false;
        if(freadEndianSafe(&ubmZNormMean[mf][nr],1,sizeof(ubmZNormMean[mf][nr]),modelFile) == sizeof(ubmZNormMean[mf][nr]))
        {
          if(freadEndianSafe(&ubmZNormVariance[mf][nr],1,sizeof(ubmZNormVariance[mf][nr]),modelFile) == sizeof(ubmZNormVariance[mf][nr]))
          {
            ubmZNormAvailable[mf] = true;
          }
        }
      }
      else
      {
//        zNormAvailable = false;
//        zNormMean      = 0.0;
//        zNormVariance  = 1.0;
      }
//assert(ubmZNormAvailable[mf]);
      nr++;
    }
  }
  fclose(modelFile);
  if(male)
  {
    nrUbmMale = nr;
  }
  else
  {
    nrUbmFemale = nr;
  }
  printf("Number of models in this UBM: %d\n",nr);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void SpeakerRecognition::writeModel(char *modelName, PhoneModel *model)
{
  FILE *modelFile = fopen(modelName, "wb");
  if(modelFile!=NULL)
  {
    int nrM = 1, nrC = 1, vectS = vectorSize;
    fwriteEndianSafe(&vectS,1,sizeof(vectS),modelFile);
    fwriteEndianSafe(&nrM,1,sizeof(nrM),modelFile);
    fwriteEndianSafe(&nrC,1,sizeof(nrC),modelFile);
    model->writeModel(modelFile);
    if(zNormAvailable)
    {
      fwriteEndianSafe(&zNormMean,1,sizeof(zNormMean),modelFile);
      fwriteEndianSafe(&zNormVariance,1,sizeof(zNormVariance),modelFile);
    }
    fclose(modelFile);
  }
  else
  {
    USER_ERROR("Could not read the model file.");
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

Gaussian *SpeakerRecognition::runTrials(char *tasklist, char *modelDir, char *outputPrefix, PhoneModel *zNormSpeakerModel)
{
  Gaussian *resGaussian = NULL;
  Vector   *scoreVector = NULL;

  if(!zNormAvailable)
  {
    resGaussian = new Gaussian(1);
    scoreVector = new Vector(1);
  }

  FILE *taskFile = fopen(tasklist,"rt");
  char spkid[MAXSTR];
  char gender[MAXSTR];
  char audioLabel[MAXSTR];
  char audioName[MAXSTR];
  char sadName[MAXSTR];
  char channel[MAXSTR];
  char trueFalse[2];
  trueFalse[0] = 'F';
  trueFalse[1] = 'T';
  if(taskFile != NULL)
  {
    while(!feof(taskFile))
    {
      strcpy(spkid,"");
      strcpy(gender,"");
      strcpy(audioLabel,"");
      strcpy(audioName,"");
      strcpy(sadName,"");
      strcpy(channel,"");
      int resNotUsed = fscanf(taskFile,"%s %s %s %s %s %s",spkid,gender,audioLabel, channel, audioName, sadName);
      if(strlen(spkid) > 0)
      {
        ////////////////// -------------------------------------------------------------------------------------------
        // Now do the actual run:
        ////////////////// -------------------------------------------------------------------------------------------


        int          decision     = 0;
        double       score        = -9.0e10;
        FeaturePool *featurePool  = NULL;
        PhoneModel  *speakerModel = NULL;

        featurePool = new FeaturePool(audioName,FILE_TYPE_RAW_SPKREC,false,1.0,false,true);  // Only prepare...
        int segID = featurePool->claimSegmentationID();
        FILE *listF = fopen(sadName, "rt");
        if(listF == NULL)
        {
          USER_ERROR("Could not read SAD RTTM file.");
        }
        featurePool->readRTTM(segID,NULL,-1,-1,3,listF);
        fclose(listF);
        featurePool->createNewPool(audioName, NULL, NULL, false, false, segID);
        assert(vectorSize == featurePool->getVectorSize());

        PhoneModel **ubm = NULL;
        int mf    = 0;
        int nrUbm = 0;
        char speakerModelName[MAXSTR];
        strcpy(speakerModelName, modelDir);
        if(strcmp(gender,"m") == 0)
        {
          mf = 0;
          strcat(speakerModelName, "male-");
          ubm = ubmMale;
          nrUbm = nrUbmMale;
        }
        else if(strcmp(gender,"f") == 0)
        {
          mf = 1;
          strcat(speakerModelName, "female-");
          ubm = ubmFemale;
          nrUbm = nrUbmFemale;
        }
        else
        {
          USER_ERROR("Could not determine gender.");
        }

        if(outputStats)
        {
          char outFileName[MAXSTR];
          sprintf(outFileName,"%s.%s.N.stats",gender,spkid);
          FILE *outFileN = fopen(outFileName,"ab");
          sprintf(outFileName,"%s.%s.F.stats",gender,spkid);
          FILE *outFileF = fopen(outFileName,"ab");
/*
          sprintf(outFileName,"%s.%s.features",gender,spkid);
          FILE *outFileFEA = fopen(outFileName,"at");
          bool isLast = false;
          Vector *vector = featurePool->getFirstVectorFirstSegment(&readPool_curSeg, segID,-1,true, &isLast);
          while(vector != NULL)
          {
            int time      = featurePool->getCurSegmentStart(&readPool_curSeg);
            int length    = featurePool->getCurSegmentLen(&readPool_curSeg);
            for(int aa=time;aa<time+length;aa++)
            {
              Vector *feaVector = featurePool->getVector(&readPool_curSeg,aa);
              for(int ii=0;ii<feaVector->len();ii++)
              {
                fprintf(outFileFEA,"%g ",feaVector->getValue(ii));
              }
              fprintf(outFileFEA,"\n");
            }
            vector = featurePool->getFirstVectorNextSegment(&readPool_curSeg,true, &isLast);
          }
          fclose(outFileFEA);
*/
          assert(outFileN != NULL);
          assert(outFileF != NULL);
          ubm[0]->adapt_clear();
          ubm[0]->adapt_setAcumulators(-1, segID, featurePool);
          ubm[0]->writeAccumulators(outFileN, outFileF);
          fclose(outFileN);
          fclose(outFileF);
        }
        else
        {
          strcat(speakerModelName,spkid);
          strcat(speakerModelName,".am");
  
          if(zNormSpeakerModel == NULL)
          {
            speakerModel = readModel(speakerModelName,true);
          }
          else
          {
            speakerModel = zNormSpeakerModel;
          }
  
          int  observationLength = 0;
          int  lookAheadSize     = vectorSize*HALF_GMMLOOKAHEAD;
          double speakerScore    = 0.0;
          double ubmScore[nrUbm];
          for(int aa=0;aa<nrUbm;aa++)
          {
            ubmScore[aa] = 0.0;
          }
          bool isLast            = false;
          Vector *vector = featurePool->getFirstVectorFirstSegment(&readPool_curSeg, segID,-1,true, &isLast);
          while(vector != NULL)
          {
            int time      = featurePool->getCurSegmentStart(&readPool_curSeg);
            int length    = featurePool->getCurSegmentLen(&readPool_curSeg);
            observationLength += length;
            for(int aa=time;aa<time+length;aa++)
            {
              Vector **feaVector = featurePool->getVectorList(&readPool_curSeg,aa);
              double *vList = new double[lookAheadSize];
              int cntList = 0;
              for(int i=0;i<vectorSize;i++)
              {
                for(int i2=0;i2<HALF_GMMLOOKAHEAD;i2++)
                {
                  vList[cntList] = feaVector[i2]->getValue(i);
                  cntList++;
                }
              }
              speakerScore += speakerModel->getLookaheadLogP(vList,aa,false);
              for(int ubmi=0;ubmi<nrUbm;ubmi++)
              {
                ubmScore[ubmi] += ubm[ubmi]->getLookaheadLogP(vList,aa,false);
              }
              delete[] vList;
            }
            vector = featurePool->getFirstVectorNextSegment(&readPool_curSeg,true, &isLast);
          }
          if(nrUbm == 1)
          {
            score = (speakerScore - ubmScore[0]) / ((double)observationLength);
          }
          else
          {
            Gaussian *nNormG = new Gaussian(1);
            Vector *nNormV = new Vector(1);
            double tRes = 0.0;
            for(int ubmi=1;ubmi<nrUbm;ubmi++)
            {
              tRes = (ubmScore[ubmi] - ubmScore[0])  / ((double)observationLength);
              tRes = (tRes - ubmZNormMean[mf][ubmi]) / ubmZNormVariance[mf][ubmi];
              nNormV->setValue(0,tRes);
              nNormG->train(1.0,nNormV);
            }
            nNormG->trainFinish();
            double nm = nNormG->getMean()->getValue(0);
            double nv = sqrt(nNormG->getVariance()->getValue(0));
            delete nNormG;
            delete nNormV;
            score = (((speakerScore - ubmScore[0]) / ((double)observationLength)) - zNormMean) / zNormVariance;
            score = (score - nm) / nv;
          }
  /*        if(zNormAvailable)
          {
            score = (score - zNormMean) / zNormVariance;
          }
          else
          {
            scoreVector->setValue(0,score);
            resGaussian->train(1.0,scoreVector);
          }
  */
          if(zNormSpeakerModel == NULL)
          {
            delete speakerModel;
          }
        }
        delete featurePool;

        ////////////////// -------------------------------------------------------------------------------------------
        // Print it and start the next run...
        ////////////////// -------------------------------------------------------------------------------------------
        if((outputPrefix != NULL) && strlen(outputPrefix)>0)
        {
          printf("%s ",outputPrefix);
        }
        printf("%s %s %s X %c %f\n",gender, spkid, audioLabel,trueFalse[decision],score);
        fflush(stdout);
      }
    }
    fclose(taskFile);
  }
  else
  {
    USER_ERROR("Could not read the tasklist");
  }

  if(!zNormAvailable)
  {
    delete scoreVector;
    resGaussian->trainFinish();
  }
  return resGaussian;
}
