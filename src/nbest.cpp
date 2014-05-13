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
#include "nbest.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The constructor will create the N-Best list structure. The original WLR strings are flatten.
/// For each sentence occurence, the one with the best AM score is retained, duplicates are removed.
/////////////////////////////////////////////////////////////////////////////////////////////////////

NBest::NBest(WLRType *path, char **voc, int stS, LanguageModel *lm, double orgLm, double orgTrP)
{
  // Parameter initialisation:
  admin               = NULL;
  languageModel       = lm;
  vocabulary          = voc;
  startOfSentenceWord = stS;
  org_LmScale         = orgLm;
  org_transPenalty    = orgTrP;
  int nrOfPaths       = 1;  
  // Get the size of the initial NBest matrix:
  int length          = getNBestArraySize(path,&nrOfPaths);
  length++;   // for the extra parameter at the end of each sentence: path to parent and index at parent... 
  // Create the first (temporary) matrix and initialise:
  NBestType nbestArray[length*nrOfPaths];
  for(int a = 0;a<nrOfPaths;a++)
  {
    int key = a*length;
    for(int b = 0;b<length;b++)
    {
      nbestArray[key+b].wordID = 0;
      nbestArray[key+b].AM     = 0.0;
      nbestArray[key+b].LM     = 0.0;      
    }
  }
  
  // Fill the matrix (not flatten yet and LM is not valid for all paths!):
  int nr      = 0;
  int emptyNr = length;
  fillNBestArray(path,nbestArray,&nr,&emptyNr,length);

  // Now flatten the matrix (no cross-paths and correct LM scores):
  for(int a = 0;a<nrOfPaths;a++)
  {
    int keySource = a*length;    
    int keyDest   = keySource;
    int hist[2]   = {startOfSentenceWord,-1};
    while(nbestArray[keySource].wordID != -5555)
    {
      int r = nbestArray[keySource].wordID; 
      while(r != -5555 && r < 0)
      {
        keySource = -1*r;
        r = nbestArray[keySource].wordID; 
      }
      if(r != -5555)
      {
        nbestArray[keyDest].wordID = nbestArray[keySource].wordID;
        nbestArray[keyDest].AM     = nbestArray[keySource].AM;
        
        if(nbestArray[keyDest].wordID == startOfSentenceWord)
        {
          nbestArray[keyDest].LM = 0.0;
        }
        else
        {
          nbestArray[keyDest].LM = languageModel->getP(nbestArray[keyDest].wordID,hist,hist);
        }
        
        keyDest++;
        keySource++;
      }
    }
    nbestArray[keyDest].wordID = -5555;
  }

  adminLength = 0;
  /// \todo Do we need to check which duplicate is the best?
  // Now remove duplicates:
  for(int a=0;a<nrOfPaths;a++)
  {
    int keySource = a*length; 
    if(nbestArray[keySource].wordID != -5555)
    {
      adminLength++;
      for(int b=a+1;b<nrOfPaths;b++)
      {
        int keyDest = b*length;
        // Are they the same?
        int c = keySource;
        int d = keyDest;
        while(nbestArray[c].wordID == nbestArray[d].wordID && nbestArray[c].wordID != -5555)
        {
          c++;
          d++;
        }
        if(nbestArray[c].wordID == -5555 && nbestArray[d].wordID == -5555)  // Bingo, these are identical!
        {
          nbestArray[keyDest].wordID = -5555; // Removed from the system...
        }
      }        
    }
  }
  // Let's create the final N-Best list:
  admin = new NBestList[adminLength];
  int ad = 0;
  for(int a = 0;a<nrOfPaths;a++)
  {
    int key = a*length;
    if(nbestArray[key].wordID != -5555)
    {
      admin[ad].totAM = 0.0;
      admin[ad].totLM = 0.0;
      int b = 0;
      int s = 0;
      while(nbestArray[key+b].wordID != -5555)
      {       
        admin[ad].totAM += nbestArray[key+b].AM;
        admin[ad].totLM += nbestArray[key+b].LM;
        if(nbestArray[key+b].wordID != startOfSentenceWord)
        {
          s++;
        }
        b++;
      }
      admin[ad].nrWords    = b;
      admin[ad].noSilWords = s;
      admin[ad].wordList   = new NBestType[b];
      b = 0;
      while(nbestArray[key+b].wordID != -5555)
      {       
        admin[ad].wordList[b].AM     = nbestArray[key+b].AM;
        admin[ad].wordList[b].LM     = nbestArray[key+b].LM;
        admin[ad].wordList[b].wordID = nbestArray[key+b].wordID;
        b++;
      }
      ad++;
    }
  }
 // if(adminLength < NBEST_MAXNR_SCORE)  
  {
    printf("# warning: length of N-best list is shorter than NBEST_MAXNR_SCORE (%d and %d)\n",adminLength, NBEST_MAXNR_SCORE);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Destructor will delete the administration
/////////////////////////////////////////////////////////////////////////////////////////////////////

NBest::~NBest()
{
  if(admin != NULL)
  {
    for(int i=0;i<adminLength;i++)
    {
      delete[] admin[i].wordList;
    }
    delete[] admin;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Fills the array of the N-Best list.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int NBest::fillNBestArray(WLRType *w, NBestType *array, int *nr, int *emptyNr, int length)
{
  int res = 0;
  if(w != NULL)
  {
    w->usedAt = -4;
    res = fillNBestArray(w->previous,array,nr, emptyNr, length); 
    
    // Now for myself:    
    if(w->timeStamp >= 0)
    {
      if(w->isSil)
      {
        array[*nr+res].wordID = startOfSentenceWord;
        array[*nr+res].LM     = 0;
        array[*nr+res].AM     = w->COMBlikelihood - w->previous->COMBlikelihood;        
      }
      else
      {
        int dummy[LM_NGRAM_DEPTH];
        array[*nr+res].wordID = languageModel->getLastWordID(w->lmHistory);      
        array[*nr+res].LM     = languageModel->getP(array[*nr+res].wordID,w->previous->lmHistory,dummy)*org_LmScale
                                 - org_transPenalty;
        array[*nr+res].AM     = (w->COMBlikelihood - w->previous->COMBlikelihood) - array[*nr+res].LM;
      }           
      res++;
      array[*nr+res].wordID = -5555;
      // Now the alternatives:
      int parentNr = *nr;  
      for(int i=0;i<NBEST_DEPTH;i++)
      {
        if(w->nBest[i] != NULL && w->nBest[i]->usedAt != -4)
        {          
          *nr       = *emptyNr;
          *emptyNr += length;
          int a = fillNBestArray(w->nBest[i],array,nr,emptyNr,length);
          array[*nr+a].wordID = -1*(parentNr+res);
        }
      }
      *nr = parentNr;
    }       
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will determine how many words the longest path in the N-Best list contains and 
/// how big N is (nrOfPaths, unfiltered)
/////////////////////////////////////////////////////////////////////////////////////////////////////

int NBest::getNBestArraySize(WLRType *w, int *nrOfPaths)
{ 
  int res = 0;
  if(w != NULL)
  {
    res = getNBestArraySize(w->previous,nrOfPaths);
    
    if(w->timeStamp >= 0)
    {
      res++;
      w->usedAt = -8;    
      for(int i=0;i<NBEST_DEPTH;i++)
      {
        if(w->nBest[i] != NULL && w->nBest[i]->usedAt != -8)
        {
          (*nrOfPaths)++;
          int r2 = getNBestArraySize(w->nBest[i],nrOfPaths);
          if(r2 > res)
          {
            res = r2;
          }      
        }
      }
    }  
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Sets a reference hypothese. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

void NBest::setReference(WLRType *w)
{
  reference.nrWords    = 0;
  reference.noSilWords = 0;
  reference.totAM      = 0.0;  
  reference.totLM      = 0.0;
  
  while(w != NULL)
  {
    if(w->timeStamp >= 0)
    {
      if(w->isSil)
      {
        reference.totAM += (w->COMBlikelihood - w->previous->COMBlikelihood);
      }
      else
      {
        int dummy[LM_NGRAM_DEPTH];
        float lm = languageModel->getP(languageModel->getLastWordID(w->lmHistory),w->previous->lmHistory,dummy)*org_LmScale
                    - org_transPenalty;
        reference.totLM += lm;
        reference.totAM += (w->COMBlikelihood - w->previous->COMBlikelihood - lm);
        reference.noSilWords++;
      }           
      reference.nrWords++;
    }
    w = w->previous;  
  }
  reference.totLM = (reference.totLM+org_transPenalty*reference.noSilWords)/org_LmScale;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Calculates the score of the reference with LM-scale lmScale..
/////////////////////////////////////////////////////////////////////////////////////////////////////

int NBest::getScore(float lmScale, float transPenalty)
{
  for(int i=0;i<adminLength;i++)
  {
    admin[i].processed = false;
  }  
  float refTot = reference.totAM + reference.totLM*lmScale - transPenalty*reference.noSilWords;
  int count = 0;
  while(count < adminLength && bigger(lmScale,transPenalty,refTot))
  {
    count++;
  }
  return ((adminLength - count) * NBEST_MAXNR_SCORE / adminLength);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will return true if there is a hypothesis bigger than target that is 
/// not yet 'processed'
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool NBest::bigger(float lmScale, float transPenalty, float target)
{
  int bestI       = -1;
  int i           =  0;
  while(i<adminLength && bestI < 0)
  {
    float score = admin[i].totAM + admin[i].totLM * lmScale - transPenalty*admin[i].noSilWords;
    if(score > target && !admin[i].processed)
    {
      bestI     = i;
      admin[i].processed = true;
    }
    i++;
  }
  return (bestI != -1);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Print the N-Best hypothesis (N == max)
/////////////////////////////////////////////////////////////////////////////////////////////////////

void NBest::print(int max)
{
  int count = 0;
  for(int i=0;i<adminLength;i++)
  {
    admin[i].processed = false;
  }  
  while(count < max && printNextBest())
  {
    count++;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Print the best hypothesis which is not just as good as 'upperLimit'...
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool NBest::printNextBest()
{
  float bestScore = -9.0e300;
  int bestI       = -1;
  for(int i=0;i<adminLength;i++)
  {
    float score = admin[i].totAM + admin[i].totLM * org_LmScale - org_transPenalty*admin[i].noSilWords;
    if(score > bestScore && !admin[i].processed)
    {
      bestScore = score;
      bestI     = i;
    }
  }
  if(bestI >= 0)
  {
    admin[bestI].processed = true;
    for(int a=0;a<admin[bestI].nrWords;a++)
    {
      printf("%s ", vocabulary[admin[bestI].wordList[a].wordID]);
    }  
    printf("\n");  
  }
  return (bestI != -1);
}


