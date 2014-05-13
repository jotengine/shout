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

#include "standard.h"
#include "lexicaltree.h"
#include "phonemodel.h"
#include "articulatorystream.h"

int deltas[37][40];
double scores[37][40];
double weights[37][40];

#include "shout-misc.h"
using namespace WriteFileLittleBigEndian;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

ArticulatoryStream::ArticulatoryStream(char *fileName)
{
  assert(fileName != NULL);
  FILE *modelFile = fopen(fileName,"rb");
  if(modelFile == NULL)
  {
    USER_ERROR("Could not open the articulatory-AM file");
  }
  int numberOfClusters;
  freadEndianSafe(&vectorSize,1,sizeof(vectorSize),modelFile);
  freadEndianSafe(&numberOfModels,1,sizeof(numberOfModels),modelFile);
  if((numberOfModels / 2)*2 != numberOfModels)
  {
    USER_ERROR("The articulatory-AM is invalid!");
  }
  numberOfModels /= 2;
  freadEndianSafe(&numberOfClusters,1,sizeof(numberOfClusters),modelFile);
  assert(numberOfClusters = 1);
  if(numberOfModels <= 0 || numberOfClusters <= 0)
  {
    USER_ERROR("This articulatory-AM file does not contain any acoustic models!");
  }
  else
  {
    modelsYes = new PhoneModel*[numberOfModels];
    modelsNo  = new PhoneModel*[numberOfModels];
    for(int i=0;i<numberOfModels;i++)
    {
      modelsYes[i] = new PhoneModel(modelFile,vectorSize);
      if(modelsYes[i] == NULL)
      {
        USER_ERROR("This articulatory-AM file is corrupt! (1)");
      }
    }
    for(int i=0;i<numberOfModels;i++)
    {
      modelsNo[i]  = new PhoneModel(modelFile,vectorSize);
      if(modelsNo[i] == NULL)
      {
        USER_ERROR("This articulatory-AM file is corrupt! (2)");
      }
    }
    if(feof(modelFile))
    {
      USER_ERROR("This articulatory-AM file is invalid!");
    }
    freadEndianSafe(&numberOfPhoneModels,1,sizeof(numberOfPhoneModels),modelFile);
    probs     = new double[numberOfPhoneModels];
#ifdef TEST_ARTICULATORY
    testRT     = new double[numberOfModels];
    testRF     = new double[numberOfModels];
    testYes    = new double[numberOfPhoneModels];
    testNo     = new double[numberOfPhoneModels];
    testYesTot = new double[numberOfPhoneModels];
    testNoTot  = new double[numberOfPhoneModels];
    testOK     = new int[numberOfPhoneModels];
    testFP     = new int[numberOfPhoneModels];
    testMISS   = new int[numberOfPhoneModels];
    testAllOK     = new int[numberOfPhoneModels];
    testAllFP     = new int[numberOfPhoneModels];
    testAllMISS   = new int[numberOfPhoneModels];

    for(int i=0;i<numberOfModels;i++)
    {
      testRF[i]   = 0.0;
      testRT[i]   = 0.0;
      for(int i2=0;i2<numberOfPhoneModels;i2++)
      {
        testRulesTT[i][i2] = 0;
        testRulesFT[i][i2] = 0;
        testRulesTF[i][i2] = 0;
        testRulesFF[i][i2] = 0;
      }
    }
    for(int i=0;i<numberOfPhoneModels;i++)
    {
      for(int i2=0;i2<numberOfPhoneModels;i2++)
      {
        confusion[i][i2] = 0;
      }

      testYes[i]  = 0.0;
      testNo[i]   = 0.0;
      testYesTot[i]  = 0.0;
      testNoTot[i]   = 0.0;
      testOK[i]   = 0;
      testFP[i]   = 0;
      testMISS[i] = 0;
      testAllOK[i]   = 0;
      testAllFP[i]   = 0;
      testAllMISS[i] = 0;
    }



    classifier = new LexicalNode[numberOfPhoneModels];
    for(int i=1;i<numberOfPhoneModels;i++)
    {
      classifier[i].inputToken = new TokenType*;
      (*classifier[i].inputToken) = NULL;
      classifier[i].tokenSeq = new TokenType*[NONSIL_NUMBER_OF_STATES];
      classifier[i].tokenSeqLength = NONSIL_NUMBER_OF_STATES;
      for(int ii=0;ii<NONSIL_NUMBER_OF_STATES;ii++)
      {
        classifier[i].tokenSeq[ii] = NULL;
      }
      classifier[i].modelID = i;
      classifier[i].contextPrev = 0;
      classifier[i].contextNext = 0;
      classifier[i].contextKey  = 0;
    }

#endif
    decisionMatrix = new int[MAX_NUMBER_OF_DECISION_RULES*numberOfPhoneModels];
    for(int i=0;i<MAX_NUMBER_OF_DECISION_RULES;i++)
    {
      freadEndianSafe(&(decisionMatrix[i*numberOfPhoneModels]),numberOfPhoneModels,sizeof(int),modelFile);
    }
    if(feof(modelFile))
    {
      USER_ERROR("This articulatory-AM file is invalid!");
    }
    freadEndianSafe(&(artModelNumber[0]),MAX_NUMBER_OF_DECISION_RULES,sizeof(int),modelFile);
    fclose(modelFile);
  }


  for(int i2=0;i2<numberOfPhoneModels;i2++)
  {
    double fact = 0.0 / ((double)numberOfModels);
    for(int i=0;i<numberOfModels+1;i++)
    {
      scores[i][i2]  = 0.0;
      deltas[i][i2]  = 0;
      weights[i][i2] = fact;
    }
    weights[numberOfModels][i2] = 1.0;
  }
  FILE *wFile = fopen("/home/parlevink/huijbreg/articulatory/weights.dat","rb");
  if(wFile != NULL)
  {
    for(int i2=0;i2<numberOfPhoneModels;i2++)
    {
      for(int i=0;i<numberOfModels+1;i++)
      {
        freadEndianSafe(&(weights[i][i2]),1,sizeof(double),wFile);
      }
    }
    fclose(wFile);
  }

  PRINTF("{\n");
  for(int i2=0;i2<numberOfModels+1;i2++)
  {
    for(int i=0;i<numberOfPhoneModels;i++)
    {
      PRINTF2("%lf, ",weights[i2][i]);
    }
  }
  PRINTF("\n};\n\n\n");

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

ArticulatoryStream::~ArticulatoryStream()
{
  if(decisionMatrix != NULL)
  {
    delete[] decisionMatrix;
  }
  if(modelsYes != NULL)
  {
    for(int i=0;i<numberOfModels;i++)
    {
      delete modelsYes[i];
    }
    delete[] modelsYes;
  }
  if(modelsNo != NULL)
  {
    for(int i=0;i<numberOfModels;i++)
    {
      delete modelsNo[i];
    }
    delete[] modelsNo;
  }
  if(probs != NULL)
  {
    delete[] probs;
#ifdef TEST_ARTICULATORY
    delete[] testRT;
    delete[] testRF;

    delete[] testOK;
    delete[] testFP;
    delete[] testMISS;
    delete[] testAllOK;
    delete[] testAllFP;
    delete[] testAllMISS;
    delete[] testYes;
    delete[] testNo;
    delete[] testYesTot;
    delete[] testNoTot;
#endif
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void ArticulatoryStream::processVector(Vector *v)
{
  double probsYes[numberOfModels];
  double probsNo[numberOfModels];
  for(int i=0;i<numberOfModels;i++)
  {
    probsYes[i] = (modelsYes[i]->getLogPDFProbability(0,v));
    probsNo[i]  = (modelsNo[i]->getLogPDFProbability(0,v));

#ifdef TEST_ARTICULATORY
    testRT[i] += (probsYes[i]);
    testRF[i] += (probsNo[i]);
#endif
  }

  for(int i=0;i<numberOfPhoneModels;i++)
  {
    probs[i] = 0.0;
#ifdef TEST_ARTICULATORY
    double scoreYes = 0.0;
    double scoreNo  = 0.0;
#endif
    for(int d=0;d<numberOfModels;d++)
    {
      double weight = weights[d][i];
      if(decisionMatrix[artModelNumber[d]*numberOfPhoneModels+i] == 1)
      {
        probs[i] += probsYes[d]*weight;
        scores[d][i] += probsYes[d];
#ifdef TEST_ARTICULATORY
        scoreYes += probsYes[d]*weight;
        scoreNo  += probsNo[d]*weight;
#endif
      }
      else
      {
        probs[i] += probsNo[d]*weight;
        scores[d][i] += probsNo[d];
#ifdef TEST_ARTICULATORY
        scoreYes += probsNo[d]*weight;
        scoreNo  += probsYes[d]*weight;
#endif
      }
    }

#ifdef TEST_ARTICULATORY
    testYes[i] += (scoreYes);
    testNo[i]  += (scoreNo);
#endif
    if((probs[i] < MINIMUM_LOGP) || (isnan(probs[i])))
    {
      probs[i] = MINIMUM_LOGP;
    }
    else
    {
//      printf("%f\n",probs[i]);
    }
  }
}

#ifdef TEST_ARTICULATORY
void ArticulatoryStream::setPool(FeaturePool *p, PhoneModel **pM)
{
  pool = p;
  phoneModels = pM;

  wlrStart                 = new WLRType;
  wlrStart->nBest          = NULL;
  wlrStart->lattice        = NULL;
  wlrStart->isSil          = 1;
  wlrStart->previous       = NULL;
  wlrStart->adminNext      = NULL;
  wlrStart->timeStamp      = -1;
  wlrStart->COMBlikelihood =  0.0;
  wlrStart->LMlikelihood   =  0.0;
  wlrStart->usedAt         =    0;
  wlrStart->phoneAlignment = NULL;

}

void ArticulatoryStream::testStream(int phoneID, int firstFrame, int lastFrame, DecoderSettings  *settings, int contextKey)
{
  for(int a=1;a<numberOfPhoneModels;a++)
  {
    PhoneModel::initialiseToken(classifier[a].inputToken);
    for(int ai=0;ai<classifier[a].tokenSeqLength;ai++)
    {
      PhoneModel::initialiseToken(&(classifier[a].tokenSeq[ai]));
    }
    *classifier[a].inputToken = new TokenType;
    (*classifier[a].inputToken)->phonePath   = new PLRType;
    (*classifier[a].inputToken)->phonePath->phoneID    = -1;
    (*classifier[a].inputToken)->phonePath->timeStamp  =  firstFrame;
    (*classifier[a].inputToken)->phonePath->likelihood =  0.0;
    (*classifier[a].inputToken)->phonePath->stateOffset[0] = 0;  
    (*classifier[a].inputToken)->phonePath->stateOffset[1] = 0;  
    (*classifier[a].inputToken)->phonePath->previous   =  NULL;
    (*classifier[a].inputToken)->path        = wlrStart;
    (*classifier[a].inputToken)->lmLookAhead = NULL;
    (*classifier[a].inputToken)->lookAheadV  = 0.0;
    (*classifier[a].inputToken)->likelihood  = 0.0;
    (*classifier[a].inputToken)->next        = NULL;
  }

  for(int i=0;i<numberOfPhoneModels;i++)
  {
    for(int i2=0;i2<numberOfModels+1;i2++)
    {
      scores[i2][i] = 0.0;
    }
    testYes[i] = 0.0;
    testNo[i]  = 0.0;
  }

  for(int i=0;i<numberOfModels;i++)
  {
    testRT[i] = 0.0;
    testRF[i] = 0.0;
  }

  for(int i=firstFrame;i<=lastFrame;i++)
  {
    Vector *v = pool->getVector(i);
    processVector(v);
    for(int a=1;a<numberOfPhoneModels;a++)
    {
      phoneModels[classifier[a].modelID]->processVector(contextKey,v,i,-1,classifier[a].tokenSeq,classifier[a].tokenSeqLength,
                                                        (*classifier[a].inputToken),settings);

      *classifier[a].inputToken = NULL;

      TokenType *resToken = NULL;
      phoneModels[classifier[a].modelID]->getOutput(contextKey, classifier[a].tokenSeq, classifier[a].tokenSeqLength, &resToken, settings);
      if(resToken != NULL)
      {
        scores[numberOfModels][a] = resToken->likelihood;
      }
      PhoneModel::initialiseToken(&resToken);
    }
  }

  scores[numberOfModels][0] = -9.0e300;
  for(int i=0;i<numberOfPhoneModels;i++)
  {
    testYes[i] += scores[numberOfModels][i]*weights[numberOfModels][i];
  }

  for(int i=0;i<numberOfModels;i++)
  {
    if(decisionMatrix[artModelNumber[i]*numberOfPhoneModels+phoneID] == 1)
    {
      if(testRT[i] > testRF[i])
      {
        testRulesTT[i][phoneID]++;
      }
      else
      {
        testRulesFT[i][phoneID]++;
      }
    }
    else
    {
      if(testRT[i] <= testRF[i])
      {
        testRulesTF[i][phoneID]++;
      }
      else
      {
        testRulesFF[i][phoneID]++;
      }
    }
  }



  double winnerScore = testYes[0];
  int    winner      = 0;
  int    second      = 0;
  for(int i=1;i<numberOfPhoneModels;i++)
  {
    if(scores[numberOfModels][i]*weights[numberOfModels][i] > winnerScore)
    {
      winnerScore = scores[numberOfModels][i]*weights[numberOfModels][i];
      winner      = i;
    }
  }
  winnerScore = testYes[0];
  for(int i=1;i<numberOfPhoneModels;i++)
  {
    if(i != winner && testYes[i] > winnerScore)
    {
      winnerScore = testYes[i];
      second      = i;
    }
  }
  confusion[phoneID][winner]++;


/*
  winnerScore = scores[numberOfModels][0]*weights[numberOfModels][0];
  winner      = 0;
  second      = 0;
  for(int i=1;i<numberOfPhoneModels;i++)
  {
    if(scores[numberOfModels][i]*weights[numberOfModels][i] > winnerScore)
    {
      winnerScore = scores[numberOfModels][i]*weights[numberOfModels][i];
      winner      = i;
    }
  }
  winnerScore = scores[numberOfModels][0]*weights[numberOfModels][0];
  for(int i=1;i<numberOfPhoneModels;i++)
  {
    if(i != winner && scores[numberOfModels][i]*weights[numberOfModels][i] > winnerScore)
    {
      winnerScore = scores[numberOfModels][i]*weights[numberOfModels][i];
      second      = i;
    }
  }
*/

  if(phoneID != winner)
  {
    second = winner;
    for(int i=0;i<numberOfModels+1;i++)
    {
      if(scores[i][phoneID] > scores[i][second])
      {
        deltas[i][phoneID]++;
      }
    }
  }






  if(winner == phoneID)
  {
    testAllOK[winner]++;
  }
  else
  {
    testAllMISS[phoneID]++;
    testAllFP[winner]++;
  }
  for(int i=0;i<numberOfPhoneModels;i++)
  {
    double len = ((double)lastFrame - firstFrame);
    if(phoneID == i)
    {
      testOK[i]++;
      testYesTot[i]+=testYes[i]/len;
    }
    else
    {
      testFP[i]++;
      testNoTot[i]+=testYes[i]/len;
    }
  }
}

void ArticulatoryStream::testPrint(FILE *outputFile)
{
  for(int i=0;i<numberOfPhoneModels;i++)
  {
    PRINTF2("Phone%02d:\t",i);
    for(int i2=0;i2<numberOfModels+1;i2++)
    {
      PRINTF2("%d ",deltas[i2][i]);
    }
    PRINTF("\n");
  }
  PRINTF("\n\n\n");


  for(int i=0;i<numberOfModels;i++)
  {
  int a=0,b=0,c=0,d=0;
    for(int i2=0;i2<numberOfPhoneModels;i2++)
    {
      a+=testRulesTT[i][i2];
      b+=testRulesFT[i][i2];
      c+=testRulesTF[i][i2];
      d+=testRulesFF[i][i2];
    }
       fprintf(outputFile,"#%1 Rule %02d: %d %d %d %d\n",i,a,b,c,d);
  }

  for(int i=0;i<numberOfPhoneModels;i++)
  {
    fprintf(outputFile,"#%3 Phone%02d:\t",i);
    for(int i2=0;i2<numberOfPhoneModels;i2++)
    {
      fprintf(outputFile,"%d ",confusion[i][i2]);
    }
    fprintf(outputFile,"\n");
  }

}

void ArticulatoryStream::testPrint2(FILE *outputFile)
{
  int tnf = 0;
  int tnc = 0;
  for(int i=0;i<numberOfPhoneModels;i++)
  {
    int nf = 0;
    int nc = 0;
    for(int i2=0;i2<numberOfPhoneModels;i2++)
    {
      if(i==i2)
      {
        nc = confusion[i][i2];
      }
      else
      {
        nf += confusion[i][i2];
      }
    }
    tnf+=nf;
    tnc+=nc;
  }
  PRINTF4("Total score: %d  %d  %f\n",tnc,tnf,((double)tnc)/((double)(tnc+tnf)));
  PRINTF("\n");

  for(int i=0;i<numberOfPhoneModels;i++)
  {
    PRINTF2("Phone%02d:\t",i);
    for(int i2=0;i2<numberOfModels+1;i2++)
    {
      PRINTF2("%d ",deltas[i2][i]);
    }
    PRINTF("\n");
  }
  PRINTF("\n");


  for(int i=0;i<numberOfPhoneModels;i++)
  {
    int tot = 0;
    for(int i2=0;i2<numberOfModels+1;i2++)
    {
      tot += deltas[i2][i];
    }
    if(tot > 0)
    {
      double fact = 0.05 / ((double)tot);
      for(int i2=0;i2<numberOfModels+1;i2++)
      {
        weights[i2][i] += fact*((double)deltas[i2][i]);
      }
    }
  }

  for(int i=0;i<numberOfPhoneModels;i++)
  {
    for(int i2=0;i2<numberOfModels+1;i2++)
    {
      if(weights[i2][i] < 0.0)
      {
        weights[i2][i] = 0.0;
      }
    }
  }
  for(int i=0;i<numberOfPhoneModels;i++)
  {
    double totW = 0.0;
    for(int i2=0;i2<numberOfModels+1;i2++)
    {
      totW += weights[i2][i];
    }
    for(int i2=0;i2<numberOfModels+1;i2++)
    {
      weights[i2][i] /= totW;
    }
  }

  PRINTF("{\n");
  for(int i2=0;i2<numberOfModels+1;i2++)
  {
    for(int i=0;i<numberOfPhoneModels;i++)
    {
      PRINTF2("%lf, ",weights[i2][i]);
    }
  }
  PRINTF("\n};\n\n\n");


  FILE *wFile = fopen("/home/parlevink/huijbreg/articulatory/weights.dat","wb");
  if(wFile != NULL)
  {
    for(int i2=0;i2<numberOfPhoneModels;i2++)
    {
      for(int i=0;i<numberOfModels+1;i++)
      {
        fwrite(&(weights[i][i2]),1,sizeof(double),wFile);
      }
    }
    fclose(wFile);
  }


}

#endif
