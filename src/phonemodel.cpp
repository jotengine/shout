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
#include "phonemodel.h"
#include "lexicaltree.h"

#include "shout-misc.h"
using namespace WriteFileLittleBigEndian;

#ifdef SEARCH_STATISTICS_ON
  int nrOfStateHistPruning = 0;
#endif

int plrCount   = 0;
int wlrCount   = 0;
int tokCount   = 0;
int trackCount = 0;

void PhoneModel::setSilTrans(double tr)
{
 mixtureSetData[0].transitionP_toNext = tr; //isSil = true;
}

extern TokenType *bestTokenInSystem;
int uniqueNumber = 291076;

//using namespace FastMath;
#define INITLOG

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This is an empty constructor
/////////////////////////////////////////////////////////////////////////////////////////////////////

PhoneModel::PhoneModel(int dim)
{
  INITLOG;

  weightRinglastPos = NULL;
  silRinglastPos  = 1;
  dimensions      = dim;
  mixtureSetData  = NULL;
  timeStamp       = NULL;
  stateMix_1      = NULL;
  stateMix_2      = NULL;
  stateMix_3      = NULL;
  isSil           = false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method deletes an entire token list. The token 'token' and all tokens that are next in the
/// list are deleted. When phone alignment is activated, also the phone history paths are deleted.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::initialiseToken(TokenType **token)
{
  TokenType *t  = (*token);
  TokenType *tn;
  if(t != NULL)
  {
    tn = t->next;
    if(doPhoneAlignment)
    {
      initialisePhonePath(t->phonePath); 
    }
    delete t;
    COUNTER_MINUS(tokCount);
    *token = NULL;
    t = tn;
    while(t != NULL)
    {
      if(doPhoneAlignment)
      {
        initialisePhonePath(t->phonePath); 
      }
      tn = t->next;
      delete t;
      COUNTER_MINUS(tokCount);
      t = tn;
    }
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor is used to load an existing phone model from file. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

PhoneModel::PhoneModel(FILE *inFile, int dim, bool onlyUseFastP)
{
  INITLOG;

  weightRinglastPos = NULL;
  silRinglastPos = 1;
  char notUsed;
  dimensions = dim;
  freadEndianSafe(statistics.name,10,1,inFile);
  freadEndianSafe(&notUsed,1,sizeof(notUsed),inFile);
  freadEndianSafe(&notUsed,1,sizeof(notUsed),inFile);
  freadEndianSafe(&statistics.nrOfGaussians,1,sizeof(statistics.nrOfGaussians),inFile);
  freadEndianSafe(&statistics.nrOfContexts,1,sizeof(statistics.nrOfContexts),inFile);
  freadEndianSafe(&statistics.maxNrOfContexts,1,sizeof(statistics.maxNrOfContexts),inFile);
  freadEndianSafe(&statistics.nrOfTrainOcc,1,sizeof(statistics.nrOfTrainOcc),inFile);
  freadEndianSafe(&statistics.likelihood,1,sizeof(statistics.likelihood),inFile);
  freadEndianSafe(&statistics.frameMeanLikelihood,1,sizeof(statistics.frameMeanLikelihood),inFile);
  freadEndianSafe(&statistics.isSil,1,sizeof(statistics.isSil),inFile);

  mixtureSetData = new MixtureSet[statistics.nrOfContexts];
  timeStamp      = new int[statistics.maxNrOfContexts]; 
  stateMix_1     = new int[statistics.maxNrOfContexts];
  stateMix_2     = new int[statistics.maxNrOfContexts];
  stateMix_3     = new int[statistics.maxNrOfContexts];
  for(int i=0;i<statistics.maxNrOfContexts;i++)
  {
    timeStamp[i] = 9272312;
  }  
  freadEndianSafe(stateMix_1,statistics.maxNrOfContexts,sizeof(int),inFile);
  freadEndianSafe(stateMix_2,statistics.maxNrOfContexts,sizeof(int),inFile);
  freadEndianSafe(stateMix_3,statistics.maxNrOfContexts,sizeof(int),inFile);

  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    freadEndianSafe(&(mixtureSetData[i].transitionP_toSelf),1,sizeof((mixtureSetData[i].transitionP_toSelf)),inFile);  
    freadEndianSafe(&(mixtureSetData[i].transitionP_toNext),1,sizeof((mixtureSetData[i].transitionP_toNext)),inFile);  
    mixtureSetData[i].state = new MixGaussian(inFile,dimensions, onlyUseFastP);
  }
  isSil = (statistics.isSil == 1.0);
  if(isSil)
  {
    mixtureSetData[0].transitionP_toSelf = -0.1249;
    mixtureSetData[0].transitionP_toNext = -0.6021;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

int PhoneModel::printInfo(Vector *v)
{
  int tot = 0;
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    tot += mixtureSetData[i].state->printInfo(v);
  }
  return tot;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::printModel(FILE *fileMean, FILE *fileVariance, FILE *fileWeight)
{
  mixtureSetData[0].state->printModel(fileMean,fileVariance,fileWeight);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::getSilSumHist(Vector **histogram)
{
  int stateNr = stateMix_1[0];  // Only possible to use the first (and only) one.
  mixtureSetData[stateNr].state->getSumHist(histogram);
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::resetPhoneAdmin()
{
  silRinglastPos = 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

PhoneModel::PhoneModel(MixGaussian *gmm, double trans)
{
  INITLOG;

  weightRinglastPos = NULL;
  silRinglastPos = 1;
  dimensions = gmm->dim();
  statistics.nrOfContexts     = 1;
  statistics.maxNrOfContexts  = 1;
  statistics.isSil            = 1.0;
  isSil = true;
  strcpy(statistics.name,"-");
  mixtureSetData = new MixtureSet[statistics.nrOfContexts];
  timeStamp      = new int[statistics.maxNrOfContexts];
  stateMix_1     = new int[statistics.maxNrOfContexts];
  stateMix_2     = NULL;
  stateMix_3     = NULL;
  timeStamp[0]   = 9272311;
  stateMix_1[0]  = 0;
  mixtureSetData[0].transitionP_toSelf = log(1.0-trans);
  mixtureSetData[0].transitionP_toNext = log(trans);
  mixtureSetData[0].state              = new MixGaussian(gmm);
  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// During destruction, the context-information (stateMix_x, mixtureSetData and timestamp) is deleted.
/////////////////////////////////////////////////////////////////////////////////////////////////////

PhoneModel::~PhoneModel()
{
  if(weightRinglastPos != NULL)
  {
    delete[] weightRinglastPos;
    weightRinglastPos = NULL;
  }
  if(stateMix_1 != NULL)
  {
    delete[] stateMix_1;
    stateMix_1 = NULL;
  }
  if(stateMix_2 != NULL)
  {
    delete[] stateMix_2;
    stateMix_2 = NULL;
  }
  if(stateMix_3 != NULL)
  {
    delete[] stateMix_3;
    stateMix_3 = NULL;
  }
  
  if(mixtureSetData != NULL)
  {
    for(int i=0;i<statistics.nrOfContexts;i++)
    {
      delete mixtureSetData[i].state;
    }
    delete[] mixtureSetData;
    mixtureSetData = NULL;
  }
  
  if(timeStamp != NULL)
  {
    delete[] timeStamp;
    timeStamp = NULL;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The copyGaussians() method copies at most maxNmbr gaussians of each HMM state to the HMM of
/// destMixGaussian. The gaussians that represent the HMM state the most will be used.
///
/// We use the first SIL model as context...
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::copyGaussians(MixGaussian *destMixGaussian, int maxNmbr)
{
  int max = statistics.nrOfContexts;
  if(max > 10)
  {
    max = 10;
  }
  for(int i=0;i<max;i++)
  {
    mixtureSetData[i].state->copyGaussians(destMixGaussian,maxNmbr);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Provides the user with some statistics about the model.
/////////////////////////////////////////////////////////////////////////////////////////////////////

ModelStats *PhoneModel::getStatistics(void)
{
  return &statistics;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// copyPhonePath() is used for phone alignment. It copies the time-history path of a phone.
/////////////////////////////////////////////////////////////////////////////////////////////////////

PLRType *PhoneModel::copyPhonePath(PLRType *pP)
{
  PLRType *phonePath = NULL;
  PLRType *p         = NULL;
  PLRType *pOld      = NULL;
  while(pP != NULL)
  {
    p = new PLRType;
    COUNTER_PLUS(plrCount);
    p->stateOffset[0] = pP->stateOffset[0];
    p->stateOffset[1] = pP->stateOffset[1];
    p->likelihood     = pP->likelihood;
    p->phoneID        = pP->phoneID;
    p->contextKey     = pP->contextKey;
    p->timeStamp      = pP->timeStamp;
    p->previous       = NULL;
    if(pOld != NULL)
    {
      pOld->previous  = p;
    }
    else
    {
      phonePath       = p;
    }
    pOld = p; 
    pP   = pP->previous;
  }  
  return phonePath;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the total number of gaussians in this model (sum of all clusters)
/////////////////////////////////////////////////////////////////////////////////////////////////////

int PhoneModel::getNumberOfGaussians()
{
  int sum = 0;
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    sum += mixtureSetData[i].state->getNumberOfGaussians();
  }
  return sum;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method deletes an entire phone list. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::initialisePhonePath(PLRType *t)
{
  while(t != NULL)
  {
    PLRType *tt = t->previous;
    delete t;
    COUNTER_MINUS(plrCount);
    t = tt;
  }
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
/// With this method, the user can obtain the output token list of this phone model at this moment in
/// time. 
///
/// With processVector(), new input can be added to the model. The 'index' parameter is used when
/// language model look-ahead is used. It will re-arrange the sorted output list according to the LMLA
/// tables output.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::getOutput(int contextKey, TokenType **token, int tLength, TokenType **outToken, DecoderSettings *settings, float *bestL)
{
  assert(contextKey >= 0);
  TokenType *t;
  if(!isSil)
  {
    // This is not SIL!
    t = token[2];
    if(t != NULL)
    {
      float passLike = mixtureSetData[stateMix_3[contextKey]].transitionP_toNext;
      if(isnan(passLike))
      {
        USER_ERROR2("(1). Found one or more tokens with 'not a number' (badly trained model) value in phone:",statistics.name);  
      }    
      else
      {
        addChain(outToken,passLike,t,0,0,settings,bestL,false);
      }
    }      
  }
  else
  {
    // This is SIL!
    float transP = mixtureSetData[stateMix_1[contextKey]].transitionP_toNext;
    if(tLength > 3)
    {
      assert(weightRinglastPos != NULL);
      t = token[silRinglastPos];
      if(t != NULL)
      {
        t->likelihood = weightRinglastPos[silRinglastPos];
      }
    }
    else
    {
      t = token[0];
    }
    if(t != NULL)
    {
      if(isnan(transP))
      {
printf("%d  %d  %g\n",contextKey,stateMix_1[contextKey],mixtureSetData[stateMix_1[contextKey]].transitionP_toNext);
        USER_ERROR2("(2). Found one or more tokens with 'not a number' (badly trained model) value in phone:",statistics.name);  
      }    
      else
      {
        addChain(outToken,transP,t,0,0,settings,bestL,false);      
      }
    }  
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the PDF probability of the first state of this phone.
/////////////////////////////////////////////////////////////////////////////////////////////////////

double PhoneModel::getLogPDFProbability(int contextKey, Vector *v)
{ 
  return mixtureSetData[stateMix_1[contextKey]].state->getLogP(v);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the state number of this phone.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int PhoneModel::getStateNr(int contextKey, int state)
{
  assert(state >= 0 && state < 3);
  int res        = -2;
  if(state == 0)
  {
    res = stateMix_1[contextKey];
  }
  if(state == 1)
  {
    res = stateMix_2[contextKey];
  }
  if(state == 2)
  {
    res = stateMix_3[contextKey];
  }
  return res;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the PDF probability of the first state of this phone.
/////////////////////////////////////////////////////////////////////////////////////////////////////

double PhoneModel::getLookaheadLogP(double *vectorList, int timeStamp, bool doSecondHalf)
{ 
  return mixtureSetData[0].state->getLookaheadLogP(vectorList, timeStamp, doSecondHalf);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the PDF probability of one of the states of this phone.
/////////////////////////////////////////////////////////////////////////////////////////////////////

double PhoneModel::getPDFProbability(int contextKey, Vector *v, int stateNr, int time)
{
  assert(stateNr >= 0 && stateNr < 3);
  int state = 0;
  if(stateNr == 0)
  {
   state = stateMix_1[contextKey];
  }
  else if(stateNr == 1)
  {
   state = stateMix_2[contextKey];
  }
  else
  {
   state = stateMix_3[contextKey];
  }

  if((time != timeStamp[state]))
  {
    timeStamp[state] = time;
    MixGaussian *model = mixtureSetData[state].state;
    mixtureSetData[state].currentVectorP = model->getP(v);
  }
  return mixtureSetData[state].currentVectorP;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

double PhoneModel::getTransition (int contextKey, bool toSelf, int stateNr)
{
  assert(stateNr >= 0 && stateNr < 3);
  if(toSelf)
  {
    if(stateNr == 0)
    {
      return exp(mixtureSetData[stateMix_1[contextKey]].transitionP_toSelf);
    }
    else if(stateNr == 1)
    {
      return exp(mixtureSetData[stateMix_2[contextKey]].transitionP_toSelf);
    }
    else
    {
      return exp(mixtureSetData[stateMix_3[contextKey]].transitionP_toSelf);
    }
  }
  else
  {
    if(stateNr == 0)
    {
      return exp(mixtureSetData[stateMix_1[contextKey]].transitionP_toNext);
    }
    else if(stateNr == 1)
    {
      return exp(mixtureSetData[stateMix_2[contextKey]].transitionP_toNext);
    }
    else
    {
      return exp(mixtureSetData[stateMix_3[contextKey]].transitionP_toNext);
    }
  }
  return 0.0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

int PhoneModel::touchPDF(int contextKey, int time, MixGaussian **updateThese, double **resultHere)
{
  int nrUpdate = 0;
  int state[3];
  state[0] = stateMix_1[contextKey];
  if((time != timeStamp[state[0]]))
  {
    timeStamp[state[0]]   = time;
    updateThese[nrUpdate] = mixtureSetData[state[0]].state;
    resultHere[nrUpdate]  = &(mixtureSetData[state[0]].currentVectorP);
    nrUpdate++;
  }
  if(!isSil)
  {
    state[1] = stateMix_2[contextKey];
    state[2] = stateMix_3[contextKey];
    if((time != timeStamp[state[1]]))
    {
      timeStamp[state[1]]   = time;
      updateThese[nrUpdate] = mixtureSetData[state[1]].state;
      resultHere[nrUpdate]  = &(mixtureSetData[state[1]].currentVectorP);
      nrUpdate++;
    }
    if((time != timeStamp[state[2]]))
    {
      timeStamp[state[2]]   = time;
      updateThese[nrUpdate] = mixtureSetData[state[2]].state;
      resultHere[nrUpdate]  = &(mixtureSetData[state[2]].currentVectorP);
      nrUpdate++;
    }
  }
  return nrUpdate;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// For each observation (Vector) of a sequence, processVector() needs to be called to update the
/// internal structure of the phone model. The user can request for the output probability token list
/// at any time by calling getOutput(). 
///
/// The inToken variable may be a single input token, but it is also allowed to be an entire list of
/// tokens, sorted on likelihood. The most likely token needs to be the first in the list.
///
/// The token pointer shall point to an array of three token lists; each of the lists form the tokens 
/// for one of the HMM states. The inToken list will be inserted in these lists according to the 
/// HMM likelihoods. The token variable is updated and needs to be provided by the user when the next
/// observation Vector is processed.
///
/// The time variable is used to add timing history and to determine if 
/// preprocessing is needed (by preProcessVector()).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::processVector(int contextKey,Vector *v, int time, int index, 
                                     TokenType **token, int tLength, TokenType *inToken, DecoderSettings *settings, float *bestL)
{
if(contextKey < 0)
{
  printf("\n%d\n",contextKey);
}
  assert(contextKey >= 0);

  TokenType *tokenNew = NULL;
  TokenType *t;

  int state[3];
  state[0] = stateMix_1[contextKey];
  if(!isSil)
  {
    state[1] = stateMix_2[contextKey];
    state[2] = stateMix_3[contextKey];
  }

  if((time != timeStamp[state[0]]))
  {
    timeStamp[state[0]] = time;
    MixGaussian *model = mixtureSetData[state[0]].state;
    mixtureSetData[state[0]].currentVectorP = model->getLogP(v);
  }

  if(!isSil)
  {
    if((time != timeStamp[state[1]]))
    {
      timeStamp[state[1]] = time;
      MixGaussian *model = mixtureSetData[state[1]].state;
      mixtureSetData[state[1]].currentVectorP = model->getLogP(v);
    }
    if((time != timeStamp[state[2]]))
    {
      timeStamp[state[2]] = time;
      MixGaussian *model = mixtureSetData[state[2]].state;
      mixtureSetData[state[2]].currentVectorP = model->getLogP(v);
    }
//    printf("VECTOR: Not SIL\n");fflush(stdout);
    for(int i=2;i>0;i--)
    {
      uniqueNumber++;
//     printf("VECTOR: Start %d\n",i);fflush(stdout);
      int j = i-1;
      tokenNew = NULL;
      t = token[j];
      if(t != NULL)
      {
        addChain(&tokenNew,(mixtureSetData[state[j]].transitionP_toNext+mixtureSetData[state[i]].currentVectorP),t,i,time,settings,bestL,false);
      }
//      printf("VECTOR: Continue\n");fflush(stdout);
      t = token[i];
      if(t != NULL)
      {
        addChain(&tokenNew,(mixtureSetData[state[i]].transitionP_toSelf+mixtureSetData[state[i]].currentVectorP),t,i,time,settings,bestL,true);
      }
//      printf("VECTOR: Init and replace\n");fflush(stdout);
      
      initialiseToken(&(token[i]));
      token[i] = tokenNew;
    }
//    printf("Input vector..\n");fflush(stdout);

    uniqueNumber++;
    t = token[0];
    TokenType *pr = NULL;
    float addScore = (mixtureSetData[state[0]].transitionP_toSelf+mixtureSetData[state[0]].currentVectorP);
    if(!isnan(addScore))
    {
      float lowestLikelihood = *bestL - settings->prune_Beam;
      if(t != NULL)
      {
        t->likelihood += addScore;
        if(t->likelihood - settings->prune_StateBeam > lowestLikelihood)
        {
          lowestLikelihood = t->likelihood - settings->prune_StateBeam;
        }
        if(t->likelihood > *bestL)
        {
          *bestL = t->likelihood;
          bestTokenInSystem = t;
        }
        if(t->lmLookAhead != NULL)
        {
          t->lmLookAhead->collissionTime = uniqueNumber;
          t->lmLookAhead->collissionNode = t;
        }
        pr = t;
        t  = t->next;
      }
      while((t != NULL))
      {
        t->likelihood += addScore;
        if(t->likelihood < lowestLikelihood)
        {
          pr->next = t->next;
          if(doPhoneAlignment)
          {
            PhoneModel::initialisePhonePath(t->phonePath);
          }
          delete t;
          COUNTER_MINUS(tokCount);
          t = pr->next;
        }
        else
        {
          if(t->lmLookAhead != NULL)
          {
            t->lmLookAhead->collissionTime = uniqueNumber;
            t->lmLookAhead->collissionNode = t;
          }

          pr = t;
          t  = t->next;
        }
      }
    }
    else
    {
      USER_ERROR2("(3). Found one or more tokens with 'not a number' (badly trained model) value in phone:",statistics.name);  
    }    

//    printf("addChain\n");fflush(stdout);
    t = inToken;
    if(t != NULL)
    {
      if(index != -1)
      {
        addChain_lmla(&(token[0]),(mixtureSetData[state[0]].currentVectorP),t,index,settings,bestL);
      }
      else
      {
        addChain(&(token[0]),(mixtureSetData[state[0]].currentVectorP),t,0,0,settings,bestL,true);
      }
    }
//    printf("Done vector\n");fflush(stdout);



  }
  else
  {
//    printf("VECTOR: THIS IS SIL!  %d\n",tLength);fflush(stdout);
    // This is SIL:
    uniqueNumber++;
    float addScore = mixtureSetData[state[0]].currentVectorP;
    if(tLength > 3)
    {
      if(weightRinglastPos == NULL)
      {
        weightRinglastPos = new float[tLength];
        for(int a=0;a<tLength;a++)
        {
          weightRinglastPos[a] = 0.0;
        }
      }

      // delete last token in the ring buffer. Now this is the first token:
      initialiseToken(&(token[silRinglastPos]));

      TokenType *bestToken = token[0];
      TokenType *t         = token[0];
      while(t != NULL)
      {
        if(t->likelihood > bestToken->likelihood)
        {
          bestToken = t;
        }
        t = t->next;
      }

      token[silRinglastPos] = NULL;

      if(bestToken != NULL)
      {
        token[silRinglastPos] = new TokenType;
        COUNTER_PLUS(tokCount);
        token[silRinglastPos]->likelihood  = bestToken->likelihood;
        token[silRinglastPos]->lmLookAhead = bestToken->lmLookAhead;
        token[silRinglastPos]->lookAheadV  = bestToken->lookAheadV;
        token[silRinglastPos]->path        = bestToken->path;
        if(doPhoneAlignment)
        {
          token[silRinglastPos]->phonePath = PhoneModel::copyPhonePath(bestToken->phonePath);
        }
        token[silRinglastPos]->next        = NULL;
        weightRinglastPos[silRinglastPos]  = token[silRinglastPos]->likelihood;
      }
      // Now add the score.. (we KNOW that there is only one token in the list!)
      for(int i=1;i<tLength;i++)
      {
        weightRinglastPos[i] += addScore;
      }
      silRinglastPos = silRinglastPos+1;
      if(silRinglastPos == tLength)
      {
        silRinglastPos = 1;
      }
    }
    else
    {
//    printf("SIL 1...\n");fflush(stdout);
      addScore += mixtureSetData[state[0]].transitionP_toSelf; // We don't use transP if we use minimum duration...
//    printf("SIL 2...\n");fflush(stdout);
    }
    tokenNew = NULL;

    t = token[0];
    TokenType *pr = NULL;
    if(!isnan(addScore))
    {
//    printf("SIL 3...\n");fflush(stdout);
      float lowestLikelihood = *bestL - settings->prune_Beam;
      if(t != NULL)
      {
        t->likelihood += addScore;
        if(t->likelihood - settings->prune_StateBeam > lowestLikelihood)
        {
          lowestLikelihood = t->likelihood - settings->prune_StateBeam;
        }
        if(t->likelihood > *bestL)
        {
          *bestL = t->likelihood;
          bestTokenInSystem = t;
        }
        if(t->lmLookAhead != NULL)
        {
          t->lmLookAhead->collissionTime = uniqueNumber;
          t->lmLookAhead->collissionNode = t;
        }
        pr = t;
        t  = t->next;
      }
//    printf("SIL 4...\n");fflush(stdout);

      while((t != NULL))
      {
        t->likelihood += addScore;
        if(t->likelihood < lowestLikelihood)
        {
          pr->next = t->next;
          if(doPhoneAlignment)
          {
            PhoneModel::initialisePhonePath(t->phonePath);
          }
          delete t;
          COUNTER_MINUS(tokCount);
          t = pr->next;
        }
        else
        {
          if(t->lmLookAhead != NULL)
          {
            t->lmLookAhead->collissionTime = uniqueNumber;
            t->lmLookAhead->collissionNode = t;
          }

          pr = t;
          t  = t->next;
        }
      }
//    printf("SIL 5...\n");fflush(stdout);

    }
    else
    {
      USER_ERROR2("(4). Found one or more tokens with 'not a number' (badly trained model) value in phone:",statistics.name);  
    }

//    printf("SIL 6...\n");fflush(stdout);

    t = inToken;
    if(t != NULL)
    {
      if(index != -1)
      {
        addChain_lmla(&(token[0]),(mixtureSetData[state[0]].currentVectorP-settings->weights_SilPenalty),t,index,settings,bestL);
      }
      else
      {
        addChain(&(token[0]),(mixtureSetData[state[0]].currentVectorP-settings->weights_SilPenalty),t,0,0,settings,bestL,true);
      }
    }
//    printf("SIL 7...\n");fflush(stdout);

  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This is a helper function for addChain(). It checks if the language model history of the new
/// token is identical to the history of a token from the existing list. 
/// If so, it returns true, otherwise it returns false. When the token history is identical, the
/// function checks if the new token likelihood is better than the old. If this is the case, the 
/// old token is replaced by the new token. 
///
/// Note that this method is responsible for the fact that internally, only the first-best 
/// recognition is stored.
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool PhoneModel::replaceTokenLM(TokenType *nt,TokenType *token, float like, float *bestL)
{
  int i=0;
  while(i<LM_NGRAM_DEPTH && (nt->path->lmHistory[i] == token->path->lmHistory[i]))
  {
    i++;
  }
  if(i==LM_NGRAM_DEPTH)
  {
    if(nt->likelihood < like)
    {
      nt->likelihood  = like;
      nt->path        = token->path;
      if(like > *bestL)
      {
        *bestL            = like;
        bestTokenInSystem = nt;
      }
      if(doPhoneAlignment)
      {
        initialisePhonePath(nt->phonePath); 
        nt->phonePath = PhoneModel::copyPhonePath(token->phonePath);
      }
    }
    return true;
  }
  return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Because of the internal structure of Shout, it has to be possible to store multiple tokens in
/// one state. Therefore each state has its own token list.
/// Each token contains language model history information. Only tokens with an identical
/// history are merged together.
///
/// The addChain() method adds a token list (an input token list or a list from another state)
/// to an existing token list. A fixed likelihood (input parameter) is added to the likelihoods 
/// of the input token list. This is done, because the state transition and observation likelihood 
/// that are calculated (by processVector()) and need to be added to these likelihoods in order to get
/// the new correct value are the same for all tokens in the list.
///
/// The result token list is directly sorted to speed up pruning in a later stage. A static beam value 
/// is defined for the biggest difference between best likelihood and worst likelihood inside one state.
/// If this threshold is reached, the token is not added to the list to begin with.
/////////////////////////////////////////////////////////////////////////////////////////////////////


void PhoneModel::addChain_ordered(TokenType **tokenNew, float likelihood, TokenType *token, 
                          int stateNr, int curTime, DecoderSettings *settings, float *bestL)
{
  TokenType *nt  = *tokenNew;
  float lowestLikelihood = *bestL - settings->prune_Beam;
  if((token != NULL) && (token->likelihood+likelihood-settings->prune_StateBeam > lowestLikelihood))
  {
    lowestLikelihood = token->likelihood+likelihood-settings->prune_StateBeam;
  }
  if((nt != NULL) && (nt->likelihood-settings->prune_StateBeam > lowestLikelihood))
  {
    lowestLikelihood = nt->likelihood-settings->prune_StateBeam;
  }
  TokenType *nt2 = NULL;
  // First sort the new tokens into the existing list:  
  if(nt != NULL)
  {
    float like   = token->likelihood + likelihood;
    while(token != NULL && nt != NULL && like > lowestLikelihood)
    {
      bool replaced = false;

      while((nt != NULL) && !(replaced=replaceTokenLM(nt,token,like,bestL)) && (nt->likelihood > like))
      {
        nt2 = nt;        
        nt  = nt->next;      
      }
      if(replaced)
      {
        if(doPhoneAlignment)
        {
          if(stateNr > 0 && nt->phonePath != NULL && nt->phonePath->stateOffset[stateNr-1] == 0)
          {
            nt->phonePath->stateOffset[stateNr-1] = (unsigned char) (curTime - nt->phonePath->timeStamp);
          }
        }
        nt2 = nt;        
        nt  = nt->next;      
      }
      else
      {
        if(nt2 == NULL) // We are still at the first token!
        {
          *tokenNew = new TokenType;
          COUNTER_PLUS(tokCount);
          (*tokenNew)->likelihood  = like;
          (*tokenNew)->lmLookAhead = token->lmLookAhead;
          (*tokenNew)->lookAheadV  = token->lookAheadV;
          if(like > *bestL)
          {
            *bestL            = like;
            bestTokenInSystem = *tokenNew;
            if(*bestL - settings->prune_Beam > lowestLikelihood)
            {
              lowestLikelihood  = *bestL - settings->prune_Beam;
            }
          }
          (*tokenNew)->path        = token->path;
          if(doPhoneAlignment)
          {
            (*tokenNew)->phonePath   = PhoneModel::copyPhonePath(token->phonePath);
            if(stateNr > 0 && (*tokenNew)->phonePath != NULL && (*tokenNew)->phonePath->stateOffset[stateNr-1] == 0)
            {
              (*tokenNew)->phonePath->stateOffset[stateNr-1] = (unsigned char) (curTime - (*tokenNew)->phonePath->timeStamp);
            }
          }
          (*tokenNew)->next        = nt;
          nt2 = *tokenNew;
        }
        else
        {
          TokenType *tt            = new TokenType;
          COUNTER_PLUS(tokCount);
          tt->likelihood           = like;
          tt->lmLookAhead          = token->lmLookAhead;
          if(like > *bestL)
          {
            *bestL            = like;
            bestTokenInSystem = tt;
            if(*bestL - settings->prune_Beam > lowestLikelihood)
            {
              lowestLikelihood  = *bestL - settings->prune_Beam;
            }
          }
          tt->lookAheadV           = token->lookAheadV;
          tt->path                 = token->path;   
          if(doPhoneAlignment)
          {
            tt->phonePath            = PhoneModel::copyPhonePath(token->phonePath);   
            if(stateNr > 0 && tt->phonePath != NULL && tt->phonePath->stateOffset[stateNr-1] == 0)
            {
              tt->phonePath->stateOffset[stateNr-1] = (unsigned char) (curTime - tt->phonePath->timeStamp);
            }
          }
          tt->next                 = nt;
          nt2->next                = tt;
          nt2 = tt;
        }
      }
      token = token->next;
      if(token!=NULL)
      {
        like  = token->likelihood + likelihood;
      }
      
      /// \todo Solve this strange bug!
      if(nt2 == NULL)
      {
        printf("Major ERROR in PhoneModel. Please contact the developer (0x32)\n");
      }
    }
    
  }
  else
  {
    if(token != NULL)
    {
      *tokenNew = new TokenType;
      COUNTER_PLUS(tokCount);
      float like = token->likelihood + likelihood;
      (*tokenNew)->likelihood  = like;
      if(like > *bestL)
      {
        *bestL            = like;
        bestTokenInSystem = *tokenNew;
        if(*bestL - settings->prune_Beam > lowestLikelihood)
        {
          lowestLikelihood  = *bestL - settings->prune_Beam;
        }
      }
      (*tokenNew)->lmLookAhead = token->lmLookAhead;
      (*tokenNew)->lookAheadV  = token->lookAheadV;
      (*tokenNew)->path        = token->path;   
      if(doPhoneAlignment)
      {
        (*tokenNew)->phonePath   = PhoneModel::copyPhonePath(token->phonePath);   
        if(stateNr > 0 && (*tokenNew)->phonePath != NULL && (*tokenNew)->phonePath->stateOffset[stateNr-1] == 0)
        {
          (*tokenNew)->phonePath->stateOffset[stateNr-1] = (unsigned char) (curTime - (*tokenNew)->phonePath->timeStamp);
        }
      }
      (*tokenNew)->next        = NULL;              
      token = token->next;
      nt  = NULL;
      nt2 = *tokenNew;
    }
  }
  // Now extent the list with a copy of the new tokens:


  while(token != NULL && (token->likelihood+likelihood) > lowestLikelihood)
  {
    if(nt2 == NULL)
    {
      printf ("NT2 is NULL!! 1\n");
    }
    if(tokenNew == NULL || *tokenNew == NULL)
    {
      printf("tokenNew is NULL!\n");
    }
    
    float like               = token->likelihood + likelihood;
    TokenType *tt            = new TokenType;
    COUNTER_PLUS(tokCount);
    tt->likelihood           = like;
    tt->lmLookAhead          = token->lmLookAhead;
    if(like > *bestL)
    {
      *bestL            = like;
      bestTokenInSystem = tt;
      if(*bestL - settings->prune_Beam > lowestLikelihood)
      {
        lowestLikelihood  = *bestL - settings->prune_Beam;
      }
    }
    tt->lookAheadV           = token->lookAheadV;
    tt->path                 = token->path;   
    if(doPhoneAlignment)
    {
      tt->phonePath            = PhoneModel::copyPhonePath(token->phonePath);   
      if(stateNr > 0 && tt->phonePath != NULL && tt->phonePath->stateOffset[stateNr-1] == 0)
      {
        tt->phonePath->stateOffset[stateNr-1] = (unsigned char) (curTime - tt->phonePath->timeStamp);
      }
    }
    if(nt2 == NULL)
    {
      printf ("NT2 is NULL!! 2\n");
    }
    tt->next                 = NULL;
    nt2->next                = tt;
    nt2   = nt2->next;    
    token = token->next;
  }
}

void PhoneModel::addChain(TokenType **tokenNew, float likelihood, TokenType *token, int stateNr, int curTime, DecoderSettings *settings, float *bestL, bool checkCollission)
{
  assert(tokenNew != NULL);
  float lowestLikelihood = *bestL - settings->prune_Beam;
  if((*tokenNew != NULL) && ((*tokenNew)->likelihood - settings->prune_StateBeam > lowestLikelihood))
  {
    lowestLikelihood = (*tokenNew)->likelihood - settings->prune_StateBeam;
  }

  while(token != NULL)
  {
    float like = token->likelihood + likelihood;

    if(like > lowestLikelihood)
    {
      TokenType *addT  = new TokenType;
      COUNTER_PLUS(tokCount);

      addT->path       = token->path;

      if(doPhoneAlignment)
      {
        addT->phonePath   = PhoneModel::copyPhonePath(token->phonePath);
        if(stateNr > 0 && addT->phonePath != NULL && addT->phonePath->stateOffset[stateNr-1] == 0)
        {
          addT->phonePath->stateOffset[stateNr-1] = (unsigned char) (curTime - addT->phonePath->timeStamp);
        }
      }
      addT->lmLookAhead = token->lmLookAhead;
      addT->likelihood  = like;
      addT->lookAheadV  = token->lookAheadV;

      if(settings->prune_Lmla)
      {
        if(checkCollission && (addT->lmLookAhead->collissionTime == uniqueNumber))
        {
          if(addT->lmLookAhead->collissionNode->likelihood < like)
          {
            addT->lmLookAhead->collissionNode->likelihood = -9.0e300;
            addT->lmLookAhead->collissionNode = addT;
          }
          else
          {
            addT->likelihood = -9.0e300;
            like             = -9.0e300;
          }
        }
        else
        {
          addT->lmLookAhead->collissionTime = uniqueNumber;
          addT->lmLookAhead->collissionNode = addT;
        }
      }

      if((*tokenNew == NULL) || (like > (*tokenNew)->likelihood)) 
      {
        if(like > *bestL)
        {
          *bestL = like;
          bestTokenInSystem = *tokenNew;
        }
        if(like - settings->prune_StateBeam > lowestLikelihood)
        {
          lowestLikelihood = like - settings->prune_StateBeam;
        }
        addT->next = *tokenNew;
        *tokenNew  = addT;
      }
      else
      {
        addT->next      = (*tokenNew)->next;
        (*tokenNew)->next = addT;
      }
    }
    token = token->next;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// addChain_lmla() does exactly the same as addChain(), only the tokens are sorted with the 
/// lookahead value of the node with 'index'. This will speed up LMLA for output tokens.
///
/// Because addChain_lmla() is only used for output tokens, which are always coming from one single 
/// node and the tokens coming from one node are always sorted, it is not needed to sort the list
/// except for lookahead!
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::addChain_lmla(TokenType **tokenNew, float likelihood, TokenType *token, int index, DecoderSettings *settings, float *bestL)
{
    TokenType addT;
    addT.next       = NULL;
    while(token != NULL)
    {
      addT.path       = token->path;
      if(doPhoneAlignment)
      {
        addT.phonePath= token->phonePath;
      }
      addT.lmLookAhead= token->lmLookAhead;
      addT.likelihood = token->likelihood - token->lookAheadV;  // if LMLA is not used, lookAheadV will be zero...
      if(settings->prune_Lmla)
      {
        addT.lookAheadV = addT.lmLookAhead->lookAhead[index];
        addT.likelihood += addT.lookAheadV;
      }
      addChain(tokenNew,likelihood,&addT,0,0,settings,bestL,true);
      token = token->next;
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The SMAPLR adaptation starts with creating one single cluster (node) at which all gaussians in
/// the system will need to register. The method adapt_setInitialNode() will call the 
/// MixGaussian::adapt_setInitialNode() of all Gaussian objects that are in the mix (of this MixGaussian object) 
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::adapt_setInitialNode(Adapt_AM_TreeNode *node)
{
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    mixtureSetData[i].state->adapt_setInitialNode(node);
  }  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The PhoneModel has registered at least once to an adaptation node when this method is called.
/// The first time adapt_setInitialNode() shall be used. This method will call the Mix
/// Gaussian::adapt_setNode() method of the MixGaussian object of all contexts of the model
/// (of this PhoneModel object).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::adapt_setNode()
{
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    mixtureSetData[i].state->adapt_setNode();
  }  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The LexicalTree will align acoustic data and call this method with the state- and 
/// phone context of the observed feature Vector. 
/// The Vector will be used for training the MixGaussian object that is used in the model 
/// at this state and for this context.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::adapt_addAcumulatorData(int state, int contextKey, Vector *observation, double probability)
{
  int mix = -1;
  if(state == 0)
  {
    mix = stateMix_1[contextKey];
  }
  else if(state == 1)
  {
    mix = stateMix_2[contextKey];  
  }
  else
  {
    mix = stateMix_3[contextKey];  
  }
  if(mix != -1)
  {
    mixtureSetData[mix].state->train(false,probability,observation);
  } 
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::mapAdaptMeans()
{
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    mixtureSetData[i].state->mapAdaptMeans();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Performs the adaptation training-run. This run will set all acumulators for SIL models only!
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::adapt_setAcumulators(int useLabel, int useSegmentation, FeaturePool *usePool)
{
  assert(usePool != NULL);
  assert(useSegmentation >= 0);
  assert(isSil);
  bool isLast   = false;
  bool inFilter = true;
  Vector *trainVector = NULL;
  
  SegmentationAdmin curSeg;
  trainVector = usePool->getFirstVectorFirstSegment(&curSeg,useSegmentation,useLabel,useLabel < 0, &isLast, &inFilter);
  while(trainVector != NULL)
  {
    int contextKey        = usePool->getSegmentID(TRAIN_CONTEXTKEY,&curSeg);
    int stateNr  = stateMix_1[contextKey];
        
    if(inFilter)
    {      
      mixtureSetData[stateNr].state->train(false,1.0,trainVector);    
    }
    while(!isLast)
    {
      trainVector = usePool->getNextVector(&curSeg,&isLast, &inFilter);
      if(inFilter)
      {
        mixtureSetData[stateNr].state->train(false,1.0,trainVector);    
      }
    }
    trainVector = usePool->getFirstVectorNextSegment(&curSeg,useLabel < 0, &isLast, &inFilter);
  }          
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will call the MixGaussian::adapt_setHelperMatrices()
/// method of the MixGaussian object of all contexts in this model.  
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::adapt_setHelperMatrices()
{
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    mixtureSetData[i].state->adapt_setHelperMatrices();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// After all adaptation matrices are calculated, the adapt_adapt() method will adapt the PhoneModel.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::adapt_adapt()
{
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    mixtureSetData[i].state->adapt_adapt();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::adapt_adaptVar()
{
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    mixtureSetData[i].state->adapt_adaptVar();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::adapt_setVarTrans()
{
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    mixtureSetData[i].state->adapt_setVarTrans();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// After adapting the PhoneModel, the old mean vectors of all gaussians can be restored by calling 
/// adapt_unAdapt(). This can not be done recursivly (only one old Gaussian is stored, calling 
/// adapt_adapt twice will delete the original vectors). 
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::adapt_unAdapt()
{
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    mixtureSetData[i].state->adapt_unAdapt();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will clear all adaptation data and make the PhoneModel ready for a new adaptation run.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::adapt_clear()
{
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    mixtureSetData[i].state->adapt_clear();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will write training accumulators to file.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::writeAccumulators(FILE *file, FILE *fileF, FILE *fileST, bool doBinary)
{
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    mixtureSetData[i].state->writeAccumulators(file,fileF,fileST,doBinary);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will read training accumulators from file.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::addAccumulators(FILE *file)
{
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    mixtureSetData[i].state->addAccumulators(file);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// WriteModel writes this model to disc. It directly stores statistics and the state transition 
/// values. It uses MixGaussian::storeData() to store the state likelihoods.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void PhoneModel::writeModel(FILE *outFile)
{
  char notUsed;
  fwrite(statistics.name,1,10,outFile);
  fwrite(&notUsed,1,sizeof(notUsed),outFile);
  fwrite(&notUsed,1,sizeof(notUsed),outFile);
  fwrite(&statistics.nrOfGaussians,1,sizeof(statistics.nrOfGaussians),outFile);
  fwrite(&statistics.nrOfContexts,1,sizeof(statistics.nrOfContexts),outFile);
  fwrite(&statistics.maxNrOfContexts,1,sizeof(statistics.maxNrOfContexts),outFile);
  fwrite(&statistics.nrOfTrainOcc,1,sizeof(statistics.nrOfTrainOcc),outFile);
  fwrite(&statistics.likelihood,1,sizeof(statistics.likelihood),outFile);
  fwrite(&statistics.frameMeanLikelihood,1,sizeof(statistics.frameMeanLikelihood),outFile);
  fwrite(&statistics.isSil,1,sizeof(statistics.isSil),outFile);

  fwrite(stateMix_1,statistics.maxNrOfContexts,sizeof(int),outFile);  
  fwrite(stateMix_2,statistics.maxNrOfContexts,sizeof(int),outFile);  
  fwrite(stateMix_3,statistics.maxNrOfContexts,sizeof(int),outFile);  
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    fwrite(&(mixtureSetData[i].transitionP_toSelf),1,sizeof((mixtureSetData[i].transitionP_toSelf)),outFile);  
    fwrite(&(mixtureSetData[i].transitionP_toNext),1,sizeof((mixtureSetData[i].transitionP_toNext)),outFile);  
    mixtureSetData[i].state->storeData(outFile);  
  }

  printf("# Writing model. Gaussians: %d, totGaussians: %d, Contexts: %d\n",
          statistics.nrOfGaussians, getNumberOfGaussians(), statistics.nrOfContexts);

}


