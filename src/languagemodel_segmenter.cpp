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
#include "languagemodel_segmenter.h"

#include <math.h>

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

LanguageModel_Segmenter::LanguageModel_Segmenter(int uniNr) : LanguageModel(uniNr)
{
  transfer_uni = new int[uni_tableLength];
  transfer_bi  = new int[uni_tableLength*uni_tableLength];
  probs_bi     = new float[uni_tableLength*uni_tableLength];

  for(int i=0;i<uni_tableLength;i++)
  {
    transfer_uni[i] = 1;
  }
  for(int i=0;i<uni_tableLength*uni_tableLength;i++)
  {
    transfer_bi[i] = 0;
    probs_bi[i]    = 0.0;
  }
  currentWord   = 0;
  totalTransfer = 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

LanguageModel_Segmenter::~LanguageModel_Segmenter()
{
  if(transfer_bi != NULL)
  {
    delete[] transfer_bi;
  }
  if(transfer_uni != NULL)
  {
    delete[] transfer_uni;
  }
  if(probs_bi != NULL)
  {
    delete[] probs_bi;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LanguageModel_Segmenter::addPenalty(int word, double penalty)
{
  for(int i=0;i<uni_tableLength;i++)
  {
    if(i != word)
    {
      probs_bi[i+uni_tableLength*word] -= penalty;
    }
  }
  uni_lmData[word].p -= penalty;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LanguageModel_Segmenter::transferTo(int word)
{
  totalTransfer++;
  transfer_uni[word]++;
  transfer_bi[currentWord+word*uni_tableLength]++;
  currentWord = word;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LanguageModel_Segmenter::finishModel()
{
  transferTo(0);

  for(int i=0;i<uni_tableLength;i++)
  {
    uni_lmData[i].p = log(((float)transfer_uni[i]) / ((float)totalTransfer));
  }

  for(int i2=0;i2<uni_tableLength;i2++)
  {
    for(int i=0;i<uni_tableLength;i++)
    {
      if(i != i2 && transfer_bi[i+uni_tableLength*i2] > 0)
      {
        probs_bi[i+uni_tableLength*i2] = log(((float)transfer_bi[i+uni_tableLength*i2]) / ((float)transfer_uni[i]));
      }
      else
      {
        probs_bi[i+uni_tableLength*i2] = -10;
      }
    }
  }

  for(int i=0;i<uni_tableLength;i++)
  {
    transfer_uni[i] = 1;
  }
  for(int i=0;i<uni_tableLength*uni_tableLength;i++)
  {
    transfer_bi[i] = 1;
  }
  currentWord   = 0;
  totalTransfer = 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

float LanguageModel_Segmenter::getP(int newWordID,int *wordHistory, int *wordHistoryOut)
{
  float res = -100;
  if(wordHistory[0] > 0)
  {
    res = probs_bi[wordHistory[0] + uni_tableLength*newWordID];
  }
  else
  {
    res = uni_lmData[newWordID].p;
  }
  wordHistoryOut[0] = newWordID;
  wordHistoryOut[1] = -1;
  wordHistoryOut[2] = -1;
  return res;
}
