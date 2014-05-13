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
#include "hash.h"
#include "languagemodel.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The standard constructor only initialises some internal variables.
/////////////////////////////////////////////////////////////////////////////////////////////////////

LanguageModel::LanguageModel()
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// If needed, the destructor releases memory that is claimed for the LM arrays and the Hash tables.
/////////////////////////////////////////////////////////////////////////////////////////////////////

LanguageModel::~LanguageModel()
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor (normally used in the shout application) loads its LM data from the given file handle.
/////////////////////////////////////////////////////////////////////////////////////////////////////

LanguageModel::LanguageModel(FILE *binlmFile)
{
  binlmFile = binlmFile;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Given the current word history, the last wordID is returned.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int LanguageModel::getLastWordID(int *wordHistory)
{
  return wordHistory[0];
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The unigram probability is retrieved from the unigram array using index wordID. The validity of 
/// wordID is not checked by this method! 
/////////////////////////////////////////////////////////////////////////////////////////////////////

float LanguageModel::getUnigram(int wordID)
{
  wordID = wordID;
  return 0.0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The LM probability for the new word 'newWordID' given the LM history 'wordHistory' is returned by
/// this method. The probability is returned in the log domain.
///
/// The new LM history is stored in 'wordHistoryOut'. This last parameter may be ignored and the old 
/// history may be maintained. This makes it possible to use this method for Language Model Look-Ahead.
/// 
/// This method will check which probability (unigram, bigram or trigram) needs to be calculated and if
/// needed multiply it with a backoff value (actualy, because all probabilities are stored in the log
/// domain, the backoff is added). 
/////////////////////////////////////////////////////////////////////////////////////////////////////

float LanguageModel::getP(int newWordID,int *wordHistory, int *wordHistoryOut)
{
  newWordID      = newWordID;
  wordHistory    = wordHistory;
  wordHistoryOut = NULL;
  return 0.0;  
}

