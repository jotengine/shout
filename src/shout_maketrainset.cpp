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

using namespace StringFunctions;
using namespace WriteFileLittleBigEndian;

Shout_MakeTrainSet::Shout_MakeTrainSet(char *task, char *phoneSetFileName, char *tD)
{
  skippedLines = 0;
  vectorSize   = ASR_DEFAULT_VECTORSIZE;
  totalPhonesStored   = 0;
  trainDirName        = new char[strlen(tD)+1];
  strcpy(trainDirName,tD);
  performCluster   = false;

  if(strcmp(task,"PHONE") == 0)
  {
    performWhat = TASK_TYPE_PHONE;
    vectorSize   = ASR_DEFAULT_VECTORSIZE;
    PhoneFileReader readPhones(phoneSetFileName, &data, &nrOfModels, &nrOfSil);
  }
  else if(strcmp(task,"DIAR") == 0)
  {
    performWhat = TASK_TYPE_DIAR;
    vectorSize   = DIARIZATION_DEFAULT_MFCC_NUMBER;
    PhoneFileReader readPhones(phoneSetFileName, &data, &nrOfModels, &nrOfSil);
  }
  else if(strcmp(task,"HIST") == 0)
  {
    performWhat = TASK_TYPE_HIST_NORM;
#ifdef HISTOGRAM_NORM_SPECTRUM
    vectorSize  = DEFAULT_MEL_BANKLENGTH;
#else
    vectorSize  = (VTLN_DEFAULT_MFCC_NUMBER + VTLN_DEFAULT_MFCC_ENERGY + VTLN_DEFAULT_MFCC_ZEROCROSS);
#endif
    if(VTLN_DEFAULT_MFCC_DELTA != 0)
    {
      vectorSize *= 3;
    }
    nrOfModels  = 1;
    nrOfSil     = 1;
    data        = new DataStats[1];
    strcpy(data[0].name,"HIST-NORM ");
    for(int cli=0;cli<MAX_CLUSTERS;cli++)
    {
      data[0].nrOfSamples[cli] = 0;
    }
  }
  else
  {
    USER_ERROR("Sorry, this application is only implemented for PHONE, HIST or DIAR traindirs...");
  }
  checkDataDir(false);
  openPhoneFiles();
  printf("First training directory: (Vector size, number of models, number of SIL-models) = (%d,%d,%d)\n",vectorSize,nrOfModels,nrOfSil);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor will read a phone definition file and a hypothesis file.
/// If both are consistent, it will write the features of the phones from the MLF to the
/// output training directory. ShoutTrainMaster reads the data from this file during training.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Shout_MakeTrainSet::Shout_MakeTrainSet(char *task, char *phoneSetFileName, char *datFile, char *rawName, char *rttmFileName,
                                       char *tD, char *decisionTree, char *method, char *trainList, char *speakerFilter, float useP)
{
  char  str[MAXSTR];

  char *histNormName    = NULL; // not used at the moment. Want to use it again in the future...
  char *histNormSpkName = NULL; // same.

  useProbability = useP;


  /// \todo I'm in a hurry because our harddrive is full with accumulator data! Make sure to clean up this bit of code later!
  if(false) //accMergeName != NULL && strlen(accMergeName)>0)
  {
    char accMergeName[MAXSTR];
    strcpy(accMergeName,"fast-fix-merge");
    if(trainList == NULL || strlen(trainList) == 0)
    {
      USER_ERROR("Please provide a list of accumulators to merge (trainList)!");
    }
    FILE *acc = fopen(trainList,"rt");
    int nrContexts = 0;

    Vector **totMean = NULL;
    Vector **totVar  = NULL;
    double *totDenom = NULL;
    while(!feof(acc))
    {
      strcpy(str,"");
      int resNotUsed = fscanf(acc,"%s",str); fflush(stdout);
      if(nrContexts == 0)
      {
        FILE *accBin = fopen(str,"rb");
        assert(accBin != NULL);
        while(!feof(accBin))
        {
          // assuming standard vector size!
          nrContexts++;
          Vector *numeratorMean     = new Vector(accBin,ASR_DEFAULT_VECTORSIZE,false);
          Vector *numeratorVariance = new Vector(accBin,ASR_DEFAULT_VECTORSIZE,false);
          delete numeratorMean;
          delete numeratorVariance;
          double trainingDenominator = 0;
          freadEndianSafe(&trainingDenominator,1,sizeof(trainingDenominator),accBin);
        }
        fclose(accBin);
        printf("Total number of Gaussians is %d..\n",nrContexts);fflush(stdout);
        totMean = new Vector*[nrContexts];
        totVar  = new Vector*[nrContexts];
        totDenom = new double[nrContexts];
        for(int i=0;i<nrContexts;i++)
        {
          totDenom[i] = 0.0;
          totMean[i] = new Vector(ASR_DEFAULT_VECTORSIZE);
          totVar[i]  = new Vector(ASR_DEFAULT_VECTORSIZE);
          totMean[i]->setAllValues(0.0);
          totVar[i]->setAllValues(0.0);
        }
      }
      printf("Processing %s...\n",str); fflush(stdout);
      FILE *accBin = fopen(str,"rb");
      int num = 0;
      if(accBin != NULL)
      {
        while(!feof(accBin))
        {
          // assuming standard vector size!
          Vector *numeratorMean     = new Vector(accBin,ASR_DEFAULT_VECTORSIZE,false);
          Vector *numeratorVariance = new Vector(accBin,ASR_DEFAULT_VECTORSIZE,false);
          double trainingDenominator = 0;
          freadEndianSafe(&trainingDenominator,1,sizeof(trainingDenominator),accBin);
          totMean[num]->addElements(true,numeratorMean);
          totVar[num]->addElements(true,numeratorVariance);
          totDenom[num] += trainingDenominator;
if(num==18)
{
  printf("Denom, tot: %f  %f\n",trainingDenominator,totDenom[num]);
}
          delete numeratorMean;
          delete numeratorVariance;
          assert(num < nrContexts);
          num++;
        }
      fclose(accBin);
      }
    }
    FILE *outFile = fopen(accMergeName,"wb");
    for(int i=0;i<nrContexts-1;i++)
    {
      totMean[i]->storeData(outFile);
      totVar[i]->storeData(outFile);
if(i==18)
{
  printf("nrContexts, denom: %d  %f\n",nrContexts, totDenom[i]);
}
      fwrite(&(totDenom[i]),1,sizeof(totDenom[i]),outFile);
    }
    fclose(outFile);   
    fclose(acc);
    exit(0);
  }






  FILE_TYPE  fileType    = FILE_TYPE_RAW_DECODING;
  int clusterID          = 0;
  totalPhonesStored      = 0;

  performCluster   = false;
  performWhat = TASK_TYPE_UNDEF;
  if(strcmp(task,"PHONE") == 0)
  {
    performWhat = TASK_TYPE_PHONE;
    fileType    = FILE_TYPE_RAW_DECODING;
    vectorSize = ASR_DEFAULT_VECTORSIZE;
    if(VTLN_DEFAULT_MFCC_DELTA != 0)
    {
      vectorSize *= 3;
    }
  }
  else if(strcmp(task,"ART") == 0)
  {
    performWhat = TASK_TYPE_ARTICULATORY;
    fileType    = FILE_TYPE_RAW_DECODING;
  }
  else if(strcmp(task,"VTLN") == 0)
  {
    performWhat = TASK_TYPE_VTLN;
    fileType    = FILE_TYPE_RAW_VTLN;
  }
  else if(strcmp(task,"SAD") == 0)
  {
    performWhat = TASK_TYPE_SAD;
    fileType    = FILE_TYPE_RAW_SAD;
  }
  else if(strcmp(task,"HIST") == 0)
  {
    performWhat = TASK_TYPE_HIST_NORM;
    fileType    = FILE_TYPE_RAW_HISTOGRAM_NORM;
  }
  else if(strcmp(task,"CLUSTER") == 0)
  {
    performWhat    = TASK_TYPE_PHONE;
    fileType       = FILE_TYPE_RAW_DECODING;
    performCluster = true;
  }
  else if(strcmp(task,"DIAR") == 0)
  {
    performWhat    = TASK_TYPE_DIAR;
    fileType       = FILE_TYPE_RAW_DIARIZATION;
  }
  else if(strcmp(task,"SPKREC") == 0)
  {
    performWhat    = TASK_TYPE_PHONE;
    fileType       = FILE_TYPE_RAW_SPKREC;
  }


  if(performWhat == TASK_TYPE_UNDEF)
  {
    USER_ERROR("No valid training type passed!");
  }

  trainDirName   = new char[strlen(tD)+1];
  strcpy(trainDirName,tD);  
  nrOfModels     = 0;  
  nrOfSil        = 0;

  PhoneFileReader readPhones(phoneSetFileName, &data, &nrOfModels, &nrOfSil);


  // Let's open the decision file:
  int decisionMatrix[MAX_NUMBER_OF_DECISION_RULES*nrOfModels];
  int artModelNumber[MAX_NUMBER_OF_DECISION_RULES];
  nrOfDecisionRules = 0;
  if(performWhat == TASK_TYPE_ARTICULATORY)
  {
    FILE *rulesFile = fopen(decisionTree, "rt");
    char rulestring[MAXSTR];
    int highestRule = -1;
    if(rulesFile!=NULL)
    {
      while(!feof(rulesFile))
      {
        if(fscanf(rulesFile,"%s",rulestring)!=0)
        {
          if(rulestring[0] == '$')
          {
            highestRule++;
            for(int i=0;i<nrOfModels;i++)
            {
              decisionMatrix[highestRule*nrOfModels+i] = 0;
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
            while((strncmp(rulestring,data[i].name,10)!=0) && (i<nrOfModels))
            {
              i++;
            }
            if(i>=nrOfModels)
            {
              USER_ERROR2("An unknown phone was found in the decision rule file: ",rulestring);
            }
            else
            {
              decisionMatrix[highestRule*nrOfModels+i] = 1;
            }
          }
        }
      }
    }
    else
    {
      USER_ERROR("The decision rule file is invalid or does not exist. I need one when I'm processing articulatory models\n");
    }
    for(int i=0;i<highestRule;i++)
    {
      int nr = 0;
      for(int i2=0;i2<nrOfModels;i2++)
      {
        if(decisionMatrix[i*nrOfModels+i2] == 1)
        {
          nr++;
        }
      }
      if(nr > 1 && nr < (nrOfModels - 3))
      {
        artModelNumber[nrOfDecisionRules] = i;
        nrOfDecisionRules++;
      }
    }

    char strName[MAXSTR];
    strcpy(strName,trainDirName);
    strcat(strName,ARTICULATORY_FILENAME);
    FILE *dataFile = fopen(strName, "wb");
    if(dataFile == NULL)
    {
      sprintf(strName,"mkdir %s",trainDirName);
      int resNotUsed = system(strName);
      strcpy(strName,trainDirName);
      strcat(strName,ARTICULATORY_FILENAME);
      dataFile = fopen(strName, "wb");
      if(dataFile == NULL)
      {
        USER_ERROR("I don't have disk access!");
      }
    }
    fwrite(&nrOfModels,1,sizeof(nrOfModels),dataFile);
    for(int i=0;i<MAX_NUMBER_OF_DECISION_RULES;i++)
    {
      fwrite(&(decisionMatrix[i*nrOfModels]),sizeof(int),nrOfModels,dataFile);
    }
    fwrite(&(artModelNumber[0]),sizeof(int),MAX_NUMBER_OF_DECISION_RULES,dataFile);
    fclose(dataFile);
  }

  printf("Initialising data directory...\n");

  checkDataDir(false);
  printf("\nInitial statistics of the data directory:\n");
  printf("===========================================\n");    
  if(performWhat != TASK_TYPE_ARTICULATORY)
  {
    printf("Phone set contains %d phones:\n",nrOfModels);
    printf("===========================================\n");    
    for(int aa=0;aa<nrOfModels;aa++)
    {
      printf("%s\t",data[aa].name);
      for(int cli=0;cli<MAX_CLUSTERS;cli++)
      {
        printf("%05d ", data[aa].nrOfSamples[cli]);
      }
      printf("\n");
    }  
  }
  else
  {
    printf("Data set contains %d models:\n",nrOfDecisionRules*2);
    printf("===========================================\n");    
    for(int aa=0;aa<nrOfDecisionRules;aa++)
    {
      printf("art%02d-y\t",aa);
      printf("%05d ", nrOfDecisionRuleSamples[aa][0]);
      printf("\n");
    }
    for(int aa=nrOfDecisionRules;aa<nrOfDecisionRules*2;aa++)
    {
      printf("art%02d-n\t",aa-nrOfDecisionRules);
      printf("%05d ", nrOfDecisionRuleSamples[aa][0]);
      printf("\n");
    }
  }

  int state = 0;
  int i = 0;
  int progress = 0;
  int progressCounter = 0;
  int contextLeft     = 0;
  int contextRight    = 0;
  skippedLines        = 0;

  openPhoneFiles();

  
  
  bool proceed   = true;
  FILE *listFile = fopen(trainList,"rt");
  FILE *inFile   = NULL;
  char rawStr[MAXSTR];
  char metaStr[MAXSTR];
  char hypStr[MAXSTR];
  int  fileNr = 0;
  if(listFile != NULL && fscanf(listFile,"%s %s %s",rawStr,metaStr,hypStr) == 0)
  {
    USER_ERROR("Invalid train-list file!");
  }
  while(proceed)
  {
    if(listFile != NULL)
    {
      rawName      = rawStr;
      datFile      = hypStr;
      rttmFileName = metaStr;
      featurePool  = NULL;
    }

    featurePool = NULL;
    printf("Creating features for '%s'\n",rttmFileName);
    inFile       = fopen(datFile, "rt");
    FILE *listF = fopen(rttmFileName, "rt");
    if(listF==NULL)
    {
      USER_ERROR("No valid RTTM file passed!");
    }
    featurePool = new FeaturePool(rawName,fileType,false,1.0,false,true);  // Only prepare...
    int segID = featurePool->claimSegmentationID();
    featurePool->readRTTM(segID,NULL,-1,-1,3,listF);
    fclose(listF);
    if(histNormSpkName == NULL || strlen(histNormSpkName) == 0)
    {
      featurePool->createNewPool(rawName, histNormName, NULL, false, false, segID);
    }
    else
    {
      featurePool->createNewPool(rawName, histNormName, histNormSpkName, false, false, segID);
    }

    vectorSize = featurePool->getVectorSize();


    if(inFile == NULL)
    {
      USER_ERROR("The input file (mlf or hyp) does not exist!");
    }
    fileNr++;
    printf("%03d: Storing training data %s...\n",fileNr, rawName);
    fflush(stdout);

    int a,b,c,d,e,f;
    char strPhone[200];
    char strPhone2[200];
    char strPhone3[200];  
    char strID[2000];
    char kl[10];
    int dataAmount[4] = {0,0,0,0};
    bool useSegment = true;
    char speakerF[MAXSTR];
    sprintf(speakerF,"%s-",speakerFilter);
    int  speakerFLen = strlen(speakerF);
    int processed    = 0;
    int notProcessed = 0;

    {
      while(!feof(inFile))
      {
        // First look what we have: label-definition, phone, or other crap:
       if(fscanf(inFile,"###################################################-| %s |-#################################",strID) != 0)
       {
          if((fileType == FILE_TYPE_RAW_SPKREC) || (speakerFLen == 1 || strncmp(strID,speakerF,speakerFLen)==0))
          {
              useSegment = true;
              processed++;
//            clusterID = 0; // Rttm based hyp that does not contain sentence labels.
          }
          else
          {
             useSegment = false;
             notProcessed++;
//            clusterID = 0;
  /*
            assert(maxSpeaker < 499);
            clusterID = -1;
            unsigned int startSpeaker = 0;
            unsigned int endSpeaker   = 0;
            while((startSpeaker < strlen(strID)) && (strID[startSpeaker] != 's'))
            {
              startSpeaker++;
            }
            startSpeaker++;
            if(startSpeaker >= strlen(strID))
            {
              clusterID = 0;
            }
            else
            {
              endSpeaker = startSpeaker+1;
              while((endSpeaker < strlen(strID)) && (strID[endSpeaker] != '-'))
              {
                endSpeaker++;
              }
              endSpeaker = endSpeaker - startSpeaker;          
              for(int i=0;i<maxSpeaker;i++)
              {
                if(strncmp(&(strID[startSpeaker]),speakers[i],endSpeaker) == 0)
                {
                  clusterID = i;
                }
              }
              if(clusterID < 0)
              {
                clusterID = maxSpeaker;
                speakers[clusterID] = new char[endSpeaker+1];
                strncpy(speakers[clusterID],&(strID[startSpeaker]),endSpeaker);
                maxSpeaker++;
              }
            }
  */
          }
        }
        else if(useSegment && fscanf(inFile,"          #Phone# %s %s %s | %d %d | %d %d | %d %d |",strPhone,strPhone2,strPhone3,&a,&b,&c,&d,&e,&f) != 0)
        {
          contextLeft=0;
          i=0;
          contextRight=0;
          strcat(strPhone,"          ");
          while((strncmp(strPhone,data[contextLeft].name,10)!=0) && (contextLeft<nrOfModels))
          {
            contextLeft++;
          }
          if(contextLeft>=nrOfModels && (performWhat == TASK_TYPE_PHONE))
          {
            fclose(inFile);
            USER_ERROR("A phone (left context) was found in the HYP file, but does not appear in the phone list!");
          }
  
          strcat(strPhone2,"          ");
          while((strncmp(strPhone2,data[i].name,10)!=0) && (i<nrOfModels))
          {
            i++;
          }
          if(i>=nrOfModels && (performWhat == TASK_TYPE_PHONE))
          {
            fclose(inFile);
            USER_ERROR("A phone (middle) was found in the HYP file, but does not appear in the phone list!");
          }
  
          strcat(strPhone3,"          ");
          while((strncmp(strPhone3,data[contextRight].name,10)!=0) && (contextRight<nrOfModels))
          {
            contextRight++;
          }
          if(contextRight>=nrOfModels && (performWhat == TASK_TYPE_PHONE))
          {
            fclose(inFile);
            USER_ERROR("A phone (right context) was found in the HYP file, but does not appear in the phone list!");
          }
  
          int numberOfObservations = 0;
          int startOfState2        = 0;
          int startOfState3        = 0;
          int cdmID                = 0;
          if(f >= 0)  // Not sil
          {
            numberOfObservations = f - a + 1;    
            startOfState2        = c - a;
            startOfState3        = e - a;
            if(performWhat != TASK_TYPE_PHONE && performWhat != TASK_TYPE_ARTICULATORY)
            {
              if(i < nrOfModels && contextLeft < nrOfModels && contextRight < nrOfModels && contextLeft != 0 && contextRight != 0)
              {
                startOfState2        = -1;
                startOfState3        = -1;            
                contextLeft  = 0;
                contextRight = 0;
                i            = 1;
                if(performWhat == TASK_TYPE_HIST_NORM)
                {
                  i = 0;
                }
              }
              else
              {
                i = nrOfModels;              
              }
            }        
          }
          else
          {
            numberOfObservations = b - a + 1;    
            startOfState2        = -1;
            startOfState3        = -1;            
            if(performWhat != TASK_TYPE_PHONE && performWhat != TASK_TYPE_ARTICULATORY)
            {
              if((performWhat != TASK_TYPE_HIST_NORM) && (i < nrOfModels && contextLeft < nrOfModels && contextRight < nrOfModels))
              {
                contextLeft  = 0;
                contextRight = 0;
                i            = 0;
              }
              else
              {
                i = nrOfModels;
              }
            }
          }
          if((startOfState2*useProbability > TRAIN_MAX_STATE_LENGTH) || (startOfState3*useProbability > TRAIN_MAX_STATE_LENGTH) || (numberOfObservations*useProbability > TRAIN_MAX_TOTAL_LENGTH))
          {
            skippedLines++;
          }
          else if(numberOfObservations > TRAIN_MIN_TOTAL_LENGTH)
          {
            if(i < nrOfModels)
            {
  /*          
              if(performWhat == TASK_TYPE_CLUSTER)
              {
                i = clusterID;
              }
  */            
              if(performWhat != TASK_TYPE_PHONE && performWhat != TASK_TYPE_ARTICULATORY)
              {
                dataAmount[i] += numberOfObservations;
              }
              if(performCluster)
              {
                cdmID = clusterID;
              }
              if(i >= 0)
              {
                if(performWhat == TASK_TYPE_ARTICULATORY)
                {
                  if(startOfState2 >= 0)
                  {
                    for(int aa=0;aa<nrOfDecisionRules;aa++)
                    {
                      if(decisionMatrix[artModelNumber[aa]*nrOfModels+i] == 1)
                      {
                        storePhone(cdmID,0,aa,0,a,-1,-1,numberOfObservations);
                        nrOfDecisionRuleSamples[aa][0]++;
                      }
                      else
                      {
                        storePhone(cdmID,0,aa+nrOfDecisionRules,0,a,-1,-1,numberOfObservations);
                        nrOfDecisionRuleSamples[aa+nrOfDecisionRules][0]++;
                      }
                    }
                  }
                }
                else
                {
                  storePhone(cdmID,contextLeft,i,contextRight,a,startOfState2,startOfState3,numberOfObservations);
                }
              }
            }
          }
        }
        else
        {
          char *resChar = fgets(str,MAXSTR,inFile);
          assert(resChar != NULL);
        }
      }
    }

    printf("\nNumber of sentences processed: %d\nNumber of sentences skipped: %d\n",processed,notProcessed);

    fclose(inFile);
    inFile = NULL;
    proceed = !(listFile == NULL || feof(listFile));
    strcpy(rawStr,"");
    if(proceed)
    {
      if(listFile != NULL && fscanf(listFile,"%s %s %s",rawStr,metaStr,hypStr) == 0)
      {
        USER_ERROR("Invalid train-list file!");
      }
      proceed = (strlen(rawStr) > 0);
    }
    printf("\n"); fflush(stdout);

    if(featurePool != NULL)
    {
      delete featurePool;
      featurePool = NULL;
    }
  }
  if(listFile != NULL)
  {
    fclose(listFile);
  }
} 

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor gets rid of some administration...
/////////////////////////////////////////////////////////////////////////////////////////////////////

Shout_MakeTrainSet::~Shout_MakeTrainSet()
{
  printf("Statistics of the training set:\n");
  printf("===========================================\n");  
  /* Now data is updated, lets safe the phone data definition: */
  closePhoneFiles(true);
  printf("===========================================\n");  
  printf("NUMBER OF SENTENCES CORRUPT (AND SKIPPED FROM THE CORRUPT POINT): %d\n",skippedLines);
  printf("===========================================\n");  

  if(trainDirName != NULL)
  {
    delete trainDirName;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Creates a file handle to each phone data file.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Shout_MakeTrainSet::openPhoneFiles()
{
  char str2[MAXSTR];
  if(performWhat == TASK_TYPE_SAD || performWhat == TASK_TYPE_VTLN)
  {
    phoneFile = new FILE*[2];  
    strcpy(str2,trainDirName);
    strcat(str2,"/");
    strncat(str2,"SIL       ",10);
    phoneFile[0] = fopen64(str2,"ab");
    if(phoneFile[0] == NULL)
    {
      USER_ERROR("Could not write in the data directory!");
    }
    strcpy(str2,trainDirName);
    strcat(str2,"/");
    strncat(str2,"SPEECH    ",10);
    phoneFile[1] = fopen64(str2,"ab");
    if(phoneFile[1] == NULL)
    {
      USER_ERROR("Could not write in the data directory!");
    }
  }
  else if(performWhat == TASK_TYPE_HIST_NORM)
  {
    phoneFile = new FILE*[1];
    strcpy(str2,trainDirName);
    strcat(str2,"/");
    strncat(str2,"HIST-NORM ",10);
    phoneFile[0] = fopen64(str2,"ab");
    if(phoneFile[0] == NULL)
    {
      USER_ERROR("Could not write in the data directory!");
    }
  }
  else if(performWhat == TASK_TYPE_ARTICULATORY)
  {
    phoneFile = new FILE*[nrOfDecisionRules*2];
    for(int i=0;i<nrOfDecisionRules;i++)
    {
      sprintf(str2,"%s/art%02d-y   ",trainDirName,i);
      phoneFile[i] = fopen64(str2,"ab");
      if(phoneFile[i] == NULL)
      {
        USER_ERROR("Could not write in the data directory!");
      }
    }
    for(int i=nrOfDecisionRules;i<nrOfDecisionRules*2;i++)
    {
      sprintf(str2,"%s/art%02d-n   ",trainDirName,i-nrOfDecisionRules);
      phoneFile[i] = fopen64(str2,"ab");
      if(phoneFile[i] == NULL)
      {
        USER_ERROR("Could not write in the data directory!");
      }
    }
  }
/*  
  else if(performWhat == TASK_TYPE_CLUSTER)
  {
    phoneFile = new FILE*[4];  
    for(int i=0;i<4;i++)
    {
      sprintf(str2,"%s/CLUST%02d   ",trainDirName,i);
      phoneFile[i] = fopen64(str2,"ab");
      if(phoneFile[i] == NULL)
      {            
        USER_ERROR("Could not write in the data directory!");
      }
    }
  }
*/  
  else
  {
    phoneFile = new FILE*[nrOfModels]; 
    for(int i=0;i<nrOfModels;i++)
    {
      strcpy(str2,trainDirName);
      strcat(str2,"/");
      strncat(str2,data[i].name,10);
      phoneFile[i] = fopen64(str2,"ab");
      if(phoneFile[i] == NULL)
      {            
        USER_ERROR("Could not write in the data directory!");
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Writes final statistics are written to disc and if printOutput is true, some statistics are
/// written to standard output.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Shout_MakeTrainSet::closePhoneFiles(bool printOutput)
{  
  char str[MAXSTR];
  int totalPhones = 0;
  if(printOutput)
  {
    printf("\n###################################################\n");
    printf("## Stored all recognised phones. The new statistics:\n");
    printf("###################################################\n");
  }
  if(phoneFile != NULL)
  {
    strcpy(str,trainDirName);
    strcat(str,PHONEFILE);
    FILE *dataFile = fopen64(str, "wb");
    if(performWhat == TASK_TYPE_SAD || performWhat == TASK_TYPE_VTLN)
    {
      int mod = 2; // Two models, two SIL...
      fwrite(&vectorSize,1,sizeof(vectorSize),dataFile);
      fwrite(&(mod),1,sizeof(mod),dataFile);
      fwrite(&(mod),1,sizeof(mod),dataFile);
      strcpy(data[0].name,"SIL       ");
      strcpy(data[1].name,"SPEECH    ");
      for(int i=0;i<2;i++)
      {
        if(printOutput)
        {
          printf("# %s\t",data[i].name);
        }
        for(int cli=0;cli<MAX_CLUSTERS;cli++)
        {
          if(printOutput)
          {
            printf("%05d ", data[i].nrOfSamples[cli]);
          }
          totalPhones+=data[i].nrOfSamples[cli];
        }
        if(printOutput)
        {
          printf("\n");
        }
        fwrite(&(data[i].name),1,10,dataFile);
        fwrite(data[i].nrOfSamples,MAX_CLUSTERS,sizeof(int),dataFile);
        fclose(phoneFile[i]);
      }
    }
    else if(performWhat == TASK_TYPE_HIST_NORM)
    {
      int mod = 1; // One model.
      fwrite(&vectorSize,1,sizeof(vectorSize),dataFile);
      fwrite(&(mod),1,sizeof(mod),dataFile);
      fwrite(&(mod),1,sizeof(mod),dataFile);
      strcpy(data[0].name,"HIST-NORM ");
      for(int i=0;i<1;i++)
      {
        if(printOutput)
        {
          printf("# %s\t",data[i].name);
        }
        for(int cli=0;cli<MAX_CLUSTERS;cli++)
        {
          if(printOutput)
          {
            printf("%05d ", data[i].nrOfSamples[cli]);
          }
          totalPhones+=data[i].nrOfSamples[cli];
        }
        if(printOutput)
        {
          printf("\n");
        }
        fwrite(&(data[i].name),1,10,dataFile);
        fwrite(data[i].nrOfSamples,MAX_CLUSTERS,sizeof(int),dataFile);
        fclose(phoneFile[i]);
      }
    }
    else if(performWhat == TASK_TYPE_ARTICULATORY)
    {
      fwrite(&vectorSize,1,sizeof(vectorSize),dataFile);
      int nrD = nrOfDecisionRules*2;
      fwrite(&(nrD),1,sizeof(nrD),dataFile);
      fwrite(&(nrD),1,sizeof(nrD),dataFile);
      for(int aa=0;aa<nrOfDecisionRules;aa++)
      {
        char name[MAXSTR];
        sprintf(name,"art%02d-y     ",aa);
        printf("%s\t",name);
        for(int cli=0;cli<MAX_CLUSTERS;cli++)
        {
          printf("%05d ", nrOfDecisionRuleSamples[aa][cli]);
        }
        printf("\n");
        fwrite(&(name),1,10,dataFile);
        fwrite(nrOfDecisionRuleSamples[aa],MAX_CLUSTERS,sizeof(int),dataFile);
        fclose(phoneFile[aa]);
      }
      for(int aa=nrOfDecisionRules;aa<nrOfDecisionRules*2;aa++)
      {
        char name[MAXSTR];
        sprintf(name,"art%02d-n     ",aa-nrOfDecisionRules);
        printf("%s\t",name);
        for(int cli=0;cli<MAX_CLUSTERS;cli++)
        {
          printf("%05d ", nrOfDecisionRuleSamples[aa][cli]);
        }
        printf("\n");
        fwrite(&(name),1,10,dataFile);
        fwrite(nrOfDecisionRuleSamples[aa],MAX_CLUSTERS,sizeof(int),dataFile);
        fclose(phoneFile[aa]);
      }
    }
/*
    else if(performWhat == TASK_TYPE_CLUSTER)
    {
      int mod = 4;
      fwrite(&vectorSize,1,sizeof(vectorSize),dataFile);
      fwrite(&(mod),1,sizeof(mod),dataFile);
      fwrite(&(mod),1,sizeof(mod),dataFile);
      for(int i=0;i<mod;i++)
      {
        sprintf(data[i].name,"CLUST%02d   ",i);
        if(printOutput)
        {
          printf("# %s\t",data[i].name);
        }
        for(int cli=0;cli<MAX_CLUSTERS;cli++)
        {
          if(printOutput)
          {
            printf("%05d ", data[i].nrOfSamples[cli]);
          }
          totalPhones+=data[i].nrOfSamples[cli];
        }
        if(printOutput)
        {
          printf("\n");
        }
        fwrite(&(data[i].name),1,10,dataFile);
        fwrite(data[i].nrOfSamples,MAX_CLUSTERS,sizeof(int),dataFile);
        fclose(phoneFile[i]);
      }
    }
*/
    else
    {    
      fwrite(&vectorSize,1,sizeof(vectorSize),dataFile);
      fwrite(&(nrOfModels),1,sizeof(nrOfModels),dataFile);
      fwrite(&(nrOfSil),1,sizeof(nrOfSil),dataFile);
      for(int i=0;i<nrOfModels;i++)
      {
        if(printOutput)
        {
          printf("# %s\t",data[i].name);
        }
        for(int cli=0;cli<MAX_CLUSTERS;cli++)
        {
          if(printOutput)
          {
            printf("%05d ", data[i].nrOfSamples[cli]);
          }
          totalPhones+=data[i].nrOfSamples[cli];
        }
        if(printOutput)
        {
          printf("\n");
        }
        fwrite(&(data[i].name),1,10,dataFile);
        fwrite(data[i].nrOfSamples,MAX_CLUSTERS,sizeof(int),dataFile);
        fclose(phoneFile[i]);
      }
    }
    fclose(dataFile);    
    delete[] phoneFile;
    if(printOutput)
    {
      printf("# TOTAL NUMBER OF PHONES, NUMBER OF PHONES ADDED: %d AND %d\n",(int)totalPhones,(int)totalPhonesStored);
    }
  }
  if(printOutput)
  {
    printf("###################################################\n");
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Checks if the acoustic data in the training directory is valid. If allowCDM is false,
/// a warning is printed when the data contains Cluster Dependent Model information.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Shout_MakeTrainSet::checkDataDir(bool allowCDM)
{
  char  str[MAXSTR];
  strcpy(str,trainDirName);
  strcat(str,PHONEFILE);
  FILE *dataFile = fopen64(str, "rb");
  if(dataFile == NULL)
  {
    sprintf(str,"mkdir %s",trainDirName);
    int resNotUsed = system(str);
    strcpy(str,"");  
    for(int i=0;i<nrOfModels;i++)
    {
      for(int cli=0;cli<MAX_CLUSTERS;cli++)
      {
        data[i].nrOfSamples[cli] = 0;
      }
    }
    for(int i=0;i<nrOfDecisionRules*2;i++)
    {
      for(int cli=0;cli<MAX_CLUSTERS;cli++)
      {
        nrOfDecisionRuleSamples[i][cli] = 0;
      }
    }
  }
  else
  {
    int v, nr, nrS;
    freadEndianSafe(&(v),1,sizeof(v),dataFile);     // This is the vector size
    freadEndianSafe(&(nr),1,sizeof(nr),dataFile);   // This is the number of phones
    freadEndianSafe(&(nrS),1,sizeof(nrS),dataFile); // This is the number of SIL phones
    if(vectorSize != v)
    {
      printf("(Vector, phones, sil): Existing traindir=(%d, %d, %d); New data=(%d, %d, %d)\n",vectorSize,nrOfModels,nrOfSil,v,nr,nrS);
      USER_ERROR("The existing AM data directory is not compatible with the new data! (the vector size does not match (SAD,VTLN,ART,PHONE mismatch?))");
    }
    if(performWhat == TASK_TYPE_PHONE)
    {
      if(nr != nrOfModels)
      {
        USER_ERROR("The existing AM data directory is not compatible with the new data! (number of phones do not match)");
      }
      if(nrS != nrOfSil)
      {
        USER_ERROR("The existing AM data directory is not compatible with the new data! (number of SIL phones do not match)");
      }
      for(int i=0;i<nrOfModels;i++)
      {
        DataStats d;
        freadEndianSafe(d.name,10,1,dataFile);
        d.name[10] = 0;
        freadEndianSafe(d.nrOfSamples,MAX_CLUSTERS,sizeof(int),dataFile);
        if(strcmp(d.name,data[i].name) != 0)
        {
          USER_ERROR("The existing AM data directory is not compatible with the new data!");
        }
        else
        {
          for(int cli=0;cli<MAX_CLUSTERS;cli++)
          {
            data[i].nrOfSamples[cli] = d.nrOfSamples[cli];
          }
        }
      }
    }
    else if(performWhat == TASK_TYPE_ARTICULATORY)
    {
      printf("Checking for articulatory features (nr of phones=%d, decision rules=%d)...\n",nr,nrOfDecisionRules);
      if(nr != nrOfDecisionRules*2)
      {
        USER_ERROR("The existing AM data directory is not compatible with the new data! (number of models do not match)");
      }
      if(nrS != nrOfDecisionRules*2)
      {
        USER_ERROR("The existing AM data directory is not compatible with the new data! (number of models do not match)");
      }
      for(int i=0;i<nrOfDecisionRules*2;i++)
      {
        assert(!feof(dataFile));
        DataStats d;
        freadEndianSafe(d.name,10,1,dataFile);
        d.name[10] = 0;
        freadEndianSafe(d.nrOfSamples,MAX_CLUSTERS,sizeof(int),dataFile);
        for(int cli=0;cli<MAX_CLUSTERS;cli++)
        {
          nrOfDecisionRuleSamples[i][cli] = d.nrOfSamples[cli];
        }
      }
    }
    else
    {
      if(performWhat == TASK_TYPE_HIST_NORM && (nr != 1 || nrS != 1))
      {
        USER_ERROR("The existing AM data directory (for HIST) is not compatible with the new data!");
      }
      if(performWhat != TASK_TYPE_HIST_NORM && performWhat != TASK_TYPE_DIAR && (nr != 2 || nrS != 2))
      {
        USER_ERROR("The existing AM data directory (for SAD or VTLN) is not compatible with the new data!");
      }
      for(int i=0;i<nr;i++)
      {
        DataStats d;
        freadEndianSafe(d.name,10,1,dataFile);
        freadEndianSafe(d.nrOfSamples,MAX_CLUSTERS,sizeof(int),dataFile);
        for(int cli=0;cli<MAX_CLUSTERS;cli++)
        {
          data[i].nrOfSamples[cli] = d.nrOfSamples[cli];
        }
      }
    }
    fclose(dataFile);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Stores a phone to the data directory.
/// - cdmID determines for which speaker cluster this data is used.
/// - contextLeft and contextRight determine which phones are at the left- and right of this phone.
/// (not always used)
/// - phoneID is the ID of the phone itself.
/// - bl, b2, b3 and e3 provide state alignment information.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Shout_MakeTrainSet::storePhone(int cdmID, int contextLeft, int phoneID, int contextRight, int b1, int b2, int b3, int numberOfObservations)
{
  assert(featurePool != NULL);

  int onceEvery = (int)(1.0 / useProbability + 0.01);
  int usingObservations = (int)(numberOfObservations / onceEvery);
  if(numberOfObservations != usingObservations * onceEvery)
  {
    usingObservations++;
  }

  if(cdmID < 0)
  {
    printf("WARNING: Cluster ID is invalid! Adding phone to cluster-0.\n");
    cdmID = 0;
  }  
  int startOfState2        = b2 / onceEvery;
  int startOfState3        = b3 / onceEvery;
  if(usingObservations > TRAIN_MAX_TOTAL_LENGTH)
  {
    USER_ERROR("Number of observations out of bounds!");
  }
  if(usingObservations>TRAIN_MIN_TOTAL_LENGTH || phoneID < nrOfSil)
  {
    totalPhonesStored++;
    fwrite(&cdmID,1,sizeof(cdmID),phoneFile[phoneID]);
    fwrite(&usingObservations,1,sizeof(usingObservations),phoneFile[phoneID]);
    fwrite(&startOfState2,1,sizeof(startOfState2),phoneFile[phoneID]);
    fwrite(&startOfState3,1,sizeof(startOfState3),phoneFile[phoneID]);                
    fwrite(&contextLeft,1,sizeof(contextLeft),phoneFile[phoneID]);
    fwrite(&contextRight,1,sizeof(contextRight),phoneFile[phoneID]);  
    if(performWhat != TASK_TYPE_ARTICULATORY)
    {
      data[phoneID].nrOfSamples[cdmID]++;
    }
    int cnt = 0;
    for(int b=b1;b<b1+numberOfObservations;b+=onceEvery)
    {
      cnt++;
      Vector *v = NULL;
      v = featurePool->getVector(NULL,b);
      if(v==NULL)
      {
        USER_ERROR("Not enough vectors in the feature file..");
      }
      if(isnan(v->getValue(1)))
      {
        USER_ERROR("One of the observation values is NotANumber!!");
      }
      v->storeData(phoneFile[phoneID]);
    }
    assert(cnt == usingObservations);
    fflush(phoneFile[phoneID]);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Sets a new feature pool for storing phones..
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Shout_MakeTrainSet::setFeaturePool(FeaturePool *fp)
{
  featurePool = fp;
}


