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
#include "languagemodel.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include "hash.h"
#include "shout-misc.h"

using namespace WriteFileLittleBigEndian;


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The standard constructor only initialises some internal variables.
/////////////////////////////////////////////////////////////////////////////////////////////////////

LanguageModel::LanguageModel()
{
  uni_lmData      = NULL;
  bi_lmData       = NULL;
  tri_lmData      = NULL;
  four_lmData     = NULL;
  uni_tableLength = 0;
  bi_tableLength  = 0;
  tri_tableLength = 0;
  four_tableLength= 0;
  bi_Hash         = NULL;
  tri_Hash        = NULL;
  four_Hash       = NULL;

  tri_HashListSearch  = NULL;
  bi_lmDataList       = NULL;
  tri_lmDataList      = NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Creates a unigram LM
/////////////////////////////////////////////////////////////////////////////////////////////////////

LanguageModel::LanguageModel(int nrUni)
{
  uni_tableLength = nrUni;
  uni_lmData      = new LMEntryType_1[uni_tableLength];
  bi_lmData       = NULL;
  tri_lmData      = NULL;
  four_lmData     = NULL;
  bi_tableLength  = 0;
  tri_tableLength = 0;
  four_tableLength= 0;
  bi_Hash         = NULL;
  tri_Hash        = NULL;
  four_Hash       = NULL;
  for(int i=0;i<uni_tableLength;i++)
  {
    uni_lmData[i].p       = 0.0;
    uni_lmData[i].backoff = 0.0;
  }

  tri_HashListSearch  = NULL;
  bi_lmDataList       = NULL;
  tri_lmDataList      = NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// If needed, the destructor releases memory that is claimed for the LM arrays and the Hash tables.
/////////////////////////////////////////////////////////////////////////////////////////////////////

LanguageModel::~LanguageModel()
{
  if(!DO_MEMORY_MAPPED_FILE_LM)
  {
    if(uni_lmData != NULL)
      delete[] uni_lmData;
    if(bi_lmData != NULL)
      delete[] bi_lmData;
    if(tri_lmData != NULL)
      delete[] tri_lmData;
    if(four_lmData != NULL)
      delete[] four_lmData;
    if(bi_lmDataList != NULL)
      delete[] bi_lmDataList;
    if(tri_lmDataList != NULL)
      delete[] tri_lmDataList;
  }
  if(bi_Hash != NULL)
    delete bi_Hash;
  if(tri_Hash != NULL)
    delete tri_Hash;
  if(four_Hash != NULL)
    delete four_Hash;
  if(tri_HashListSearch != NULL)
    delete tri_HashListSearch;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor (normally used in the shout application) loads its LM data from the given file handle.
/////////////////////////////////////////////////////////////////////////////////////////////////////

LanguageModel::LanguageModel(const char *binlmFileName)
{
  uni_lmData      = NULL;
  bi_lmData       = NULL;
  tri_lmData      = NULL;
  four_lmData     = NULL;
  uni_tableLength = 0;
  bi_tableLength  = 0;
  tri_tableLength = 0;
  four_tableLength= 0;
  bi_Hash         = NULL;
  tri_Hash        = NULL;
  four_Hash       = NULL;

  long int fileOffsets[4]  = {0,0,0,0};
  long int offset          = 0;

  tri_HashListSearch  = NULL;
  bi_lmDataList       = NULL;
  tri_lmDataList      = NULL;


  char *memoryPool = NULL;
  int   memoryFileDescriptor;

  FILE *binlmFile = fopen(binlmFileName, "rb");
  if(binlmFile == NULL)
  {
    return;
  }

  if(DO_MEMORY_MAPPED_FILE_LM)
  {
    struct stat statbuf;
    memoryFileDescriptor = open(binlmFileName, O_RDONLY);
    if(memoryFileDescriptor < 0)
    {
      USER_ERROR("Could not read the language model");
    }
    fstat(memoryFileDescriptor,&statbuf);
    memoryPool = (char*) mmap(0, statbuf.st_size, PROT_READ, MAP_SHARED, memoryFileDescriptor, 0);
    assert(memoryPool != (caddr_t) - 1);
  }

  int depth = 0;
  freadEndianSafe(&(depth),1,sizeof(depth),binlmFile);
  assert(depth == LM_NGRAM_DEPTH);

  freadEndianSafe(&(uni_tableLength),1,sizeof(uni_tableLength),binlmFile);

  if(DO_MEMORY_MAPPED_FILE_LM)
  {
    fileOffsets[0] = ftell(binlmFile);
    offset         = fileOffsets[0] + uni_tableLength * (sizeof(float)*2);
//    fseek(binlmFile,uni_tableLength * (sizeof(float)*2), SEEK_CUR);
  }
  else
  {
    uni_lmData = new LMEntryType_1[uni_tableLength];
    for(int ientry=0;ientry<uni_tableLength;ientry++)
    {
      freadEndianSafe(&(uni_lmData[ientry].p),      1,sizeof(uni_lmData[ientry].p),binlmFile);
      freadEndianSafe(&(uni_lmData[ientry].backoff),1,sizeof(uni_lmData[ientry].backoff),binlmFile);
    }
  }
  if(DO_MEMORY_MAPPED_FILE_LM)
  {
    bi_tableLength = *((int*)(memoryPool+offset));
    offset += sizeof(int);
  }
  else
  { 
    freadEndianSafe(&(bi_tableLength),1,sizeof(bi_tableLength),binlmFile);
  }
  if(bi_tableLength > 0)
  {
    if(DO_MEMORY_MAPPED_FILE_LM)
    {
      fileOffsets[1] = offset;
      offset += bi_tableLength * (sizeof(float)*2 + sizeof(int)*2);
//      fseek(binlmFile,bi_tableLength * (sizeof(float)*2 + sizeof(int)*2), SEEK_CUR);
    }
    else
    {
      bi_lmData  = new LMEntryType_2[bi_tableLength];
      for(int ientry=0;ientry<bi_tableLength;ientry++)
      {
        freadEndianSafe(bi_lmData[ientry].words,       2,sizeof(bi_lmData[ientry].words[0]),binlmFile);
        freadEndianSafe(&(bi_lmData[ientry].p),        1,sizeof(bi_lmData[ientry].p),binlmFile);
        #if LM_NGRAM_DEPTH > 1
        freadEndianSafe(&(bi_lmData[ientry].backoff),  1,sizeof(bi_lmData[ientry].backoff),binlmFile);
        #endif
      }
    }
/*
    int temp;
    freadEndianSafe(&(temp),1,sizeof(temp),binlmFile);
    assert(temp == uni_tableLength);
    if(DO_MEMORY_MAPPED_FILE_LM)
    {
      bi_lmDataList = &((LMListSearch_2*)(memoryPool+ftell(binlmFile)))[0];
      fseek(binlmFile,uni_tableLength * (sizeof(int)*2), SEEK_CUR);
    }
    else
    {
      bi_lmDataList  = new LMListSearch_2[uni_tableLength];
      for(int ientry=0;ientry<uni_tableLength;ientry++)
      {
        freadEndianSafe(&(bi_lmDataList[ientry].lowestIndex), 1,sizeof(bi_lmDataList[ientry].lowestIndex),binlmFile);
        freadEndianSafe(&(bi_lmDataList[ientry].length),      1,sizeof(bi_lmDataList[ientry].length),binlmFile);
      }
    }
*/
    bi_Hash  = new Hash(binlmFile, memoryPool,&offset);
  }
  #if LM_NGRAM_DEPTH > 1

  if(DO_MEMORY_MAPPED_FILE_LM)
  {
    tri_tableLength = *((int*)(memoryPool+offset));
    offset += sizeof(int);
  }
  else
  {
    freadEndianSafe(&(tri_tableLength),1,sizeof(tri_tableLength),binlmFile);
  }

    if(tri_tableLength > 0)
    {

      if(DO_MEMORY_MAPPED_FILE_LM)
      {
        fileOffsets[2] = offset;
        offset += tri_tableLength * (sizeof(float) + sizeof(int)*3);
//        fseek(binlmFile,tri_tableLength * (sizeof(float) + sizeof(int)*3), SEEK_CUR);
      }
      else
      {
        tri_lmData = new LMEntryType_3[tri_tableLength];
        for(int ientry=0;ientry<tri_tableLength;ientry++)
        {
          freadEndianSafe(tri_lmData[ientry].words,       3,sizeof(tri_lmData[ientry].words[0]),binlmFile);
          freadEndianSafe(&(tri_lmData[ientry].p),        1,sizeof(tri_lmData[ientry].p),binlmFile);
          #if LM_NGRAM_DEPTH > 2
          freadEndianSafe(&(tri_lmData[ientry].backoff),  1,sizeof(tri_lmData[ientry].backoff),binlmFile);
          #endif
        }
      }
      tri_Hash = new Hash(binlmFile, memoryPool, &offset);

/*
      freadEndianSafe(&(tri_tableLengthList),1,sizeof(tri_tableLengthList),binlmFile);

      if(DO_MEMORY_MAPPED_FILE_LM)
      {
        tri_lmDataList = &((LMListSearch_3*)(memoryPool+ftell(binlmFile)))[0];
        fseek(binlmFile,tri_tableLengthList * (sizeof(int)*4), SEEK_CUR);
      }
      else
      {
        tri_lmDataList  = new LMListSearch_3[tri_tableLengthList];
        for(int ientry=0;ientry<tri_tableLengthList;ientry++)
        {
          freadEndianSafe(tri_lmDataList[ientry].words,          2,sizeof(tri_lmDataList[ientry].words[0]),binlmFile);
          freadEndianSafe(&(tri_lmDataList[ientry].lowestIndex), 1,sizeof(tri_lmDataList[ientry].lowestIndex),binlmFile);
          freadEndianSafe(&(tri_lmDataList[ientry].length),      1,sizeof(tri_lmDataList[ientry].length),binlmFile);
        }
      }
      tri_HashListSearch = new Hash(binlmFile, memoryPool);
*/
    }
  #endif
/*
  #if LM_NGRAM_DEPTH > 2  

  four_tableLength = *((int*)(memoryPool+offset));
  offset += sizeof(int);
//    freadEndianSafe(&(four_tableLength),1,sizeof(four_tableLength),binlmFile);
    if(four_tableLength > 0)
    {
      four_lmData = new LMEntryType_4[four_tableLength];
      for(int ientry=0;ientry<four_tableLength;ientry++)
      {
        freadEndianSafe(four_lmData[ientry].words,       4,sizeof(four_lmData[ientry].words[0]),binlmFile);
        freadEndianSafe(&(four_lmData[ientry].p),        1,sizeof(four_lmData[ientry].p),binlmFile);
        #if LM_NGRAM_DEPTH > 3
        freadEndianSafe(&(four_lmData[ientry].backoff),  1,sizeof(four_lmData[ientry].backoff),binlmFile);
        #endif
      }
      four_Hash = new Hash(binlmFile);  
    }
  #endif
*/

  fclose(binlmFile);

  if(DO_MEMORY_MAPPED_FILE_LM)
  {
    uni_lmData      = &((LMEntryType_1*)(memoryPool+fileOffsets[0]))[0];
    bi_lmData       = &((LMEntryType_2*)(memoryPool+fileOffsets[1]))[0];
    tri_lmData      = &((LMEntryType_3*)(memoryPool+fileOffsets[2]))[0];
    //four_lmData     = ((LMEntryType_4*)(memoryPool+fileOffsets[3]))[0];
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Given the current word history, the last wordID is returned.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int LanguageModel::getLastWordID(int *wordHistory)
{
  int wordID = -1;
  if(wordHistory != NULL)
  {
    #if LM_NGRAM_DEPTH > 2
      if(wordHistory[2] >= 0)
      {
        wordID = tri_lmData[wordHistory[2]].words[2];
      }
    #endif
    #if LM_NGRAM_DEPTH > 1
      if(wordID < 0 && wordHistory[1] >= 0)
      {
        wordID = bi_lmData[wordHistory[1]].words[1];
      }
    #endif
    if(wordID < 0)
    {
      wordID = wordHistory[0];
    }
  }
  return wordID;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LanguageModel::getAllP(int *lmHistory, float *allP)
{
  if(lmHistory == NULL)
  {
    for(int i=0;i<uni_tableLength;i++)
    {
      allP[i] = uni_lmData[i].p;
    }
  }
  else
  {
    int w1 = -1;
    int w2 = -1;
    if(lmHistory[1] >= 0)
    {
      w1 = lmHistory[0];
      w2 = bi_lmData[lmHistory[1]].words[1];
    }
    else if(bi_tableLength > 0)
    {
      w2 = lmHistory[0];
    }

    bool pSet[uni_tableLength];
    for(int i=0;i<uni_tableLength;i++)
    {
      pSet[i] = false;
      allP[i] = 0.0;
    }
    if(w1 >=0 && w2 >=0)
    {
      int wordHistoryNew[2] = {w1,w2};
      int wIndex = tri_HashListSearch->getIndex(wordHistoryNew);
      if(tri_lmDataList[wIndex].words[0] == wordHistoryNew[0] && tri_lmDataList[wIndex].words[1] == wordHistoryNew[1])
      {
        int start = tri_lmDataList[wIndex].lowestIndex;
        int end   = start + tri_lmDataList[wIndex].length;
        for(int i=start;i<end;i++)
        {
          pSet[tri_lmData[i].words[2]] = true;
          allP[tri_lmData[i].words[2]] = tri_lmData[i].p;
        }
      }
      int index = getBiGramIndex(wordHistoryNew);
      float backoff = bi_lmData[index].backoff;
      for(int i=0;i<uni_tableLength;i++)
      {
        if(!pSet[i])
        {
          allP[i] = backoff;
        }
      }
      w1 = -1;
    }
    if(w2>=0)
    {
      int start = bi_lmDataList[w2].lowestIndex;
      int end   = start + bi_lmDataList[w2].length;
      for(int i=start;i<end;i++)
      {
        int w = bi_lmData[i].words[1];
        if(!pSet[w])
        {
          pSet[w]  = true;
          allP[w] += bi_lmData[i].p;
        }
      }
      float backoff = uni_lmData[w2].backoff;
      for(int i=0;i<uni_tableLength;i++)
      {
        if(!pSet[i])
        {
          allP[i] += backoff;
        }
      }
      w1 = -1;
    }
    for(int i=0;i<uni_tableLength;i++)
    {
      if(!pSet[i])
      {
        allP[i] += uni_lmData[i].p;
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LanguageModel::setUnigram(int wordID, float p)
{
  assert(wordID < uni_tableLength);
  uni_lmData[wordID].p = p;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The unigram probability is retrieved from the unigram array using index wordID. The validity of 
/// wordID is not checked by this method! 
/////////////////////////////////////////////////////////////////////////////////////////////////////

float LanguageModel::getUnigram(int wordID)
{
  return uni_lmData[wordID].p;
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
  float res = SMALLEST_NEGATIVE;
  if(wordHistory[0] == -1)       // We do not have a history yet!
  {
    res = getUniGramP(newWordID,wordHistoryOut);
    #if LM_NGRAM_DEPTH > 1
    wordHistoryOut[1] = -1;
    #endif
    #if LM_NGRAM_DEPTH > 2
    wordHistoryOut[2] = -1;
    #endif
  }
  #if LM_NGRAM_DEPTH > 1
  else if(tri_tableLength > 0 && wordHistory[1] == -1)  // We only have a uni-gram history
  {
    res = getBiGramP(newWordID,wordHistory,wordHistoryOut);
    #if LM_NGRAM_DEPTH > 2
    wordHistoryOut[2] = -1;
    #endif      
  }
  #endif
  #if LM_NGRAM_DEPTH > 2
  else if(four_tableLength > 0 && wordHistory[2] == -1)  // We only have a bi-gram history
  {
    res = getTriGramP(newWordID,wordHistory,wordHistoryOut);
  }
  #endif
  else  // Full house..
  {
    int w[LM_NGRAM_DEPTH+1];
    if(tri_tableLength == 0)
    {
      res = getBiGramP(newWordID,wordHistory,w);
      if(w[1] != -1)
      {
        wordHistoryOut[0] = bi_lmData[w[1]].words[1];
      }
      else
      {
        wordHistoryOut[0] = w[0];
      }
      #if LM_NGRAM_DEPTH > 1
      wordHistoryOut[1] = -1;
      #endif
      #if LM_NGRAM_DEPTH > 2
      wordHistoryOut[2] = -1;
      #endif      
    }
    else if(four_tableLength == 0)
    {
      #if LM_NGRAM_DEPTH > 1
      res = getTriGramP(newWordID,wordHistory,w);
      if(w[2] != -1)
      {
        int w2[2] = {tri_lmData[w[2]].words[1],tri_lmData[w[2]].words[2]};
        int r    = getBiGramIndex(w2);
        if(r == -1)
        {
          wordHistoryOut[0] = tri_lmData[w[2]].words[2];
          wordHistoryOut[1] = -1;
        }
        else
        {
          wordHistoryOut[0] = tri_lmData[w[2]].words[1];
          wordHistoryOut[1] = r;
        }
      }
      else
      {
        wordHistoryOut[0] = w[0];
        wordHistoryOut[1] = w[1];
      }
      #endif
      #if LM_NGRAM_DEPTH > 2
      wordHistoryOut[2] = -1;
      #endif
    }
    else
    {
      #if LM_NGRAM_DEPTH > 2
      res = get4GramP(newWordID,wordHistory,w);
      if(w[3] != -1)
      {
        int w2[3] = {four_lmData[w[3]].words[1],four_lmData[w[3]].words[2],four_lmData[w[3]].words[3]};
        int r    = getTriGramIndex(w2);
        int r2   = getBiGramIndex(w2);
        if(r == -1)
        {
          if(r2 == -1)
          {
            wordHistoryOut[0] = four_lmData[w[3]].words[3];
            wordHistoryOut[1] = -1;
            wordHistoryOut[2] = -1;
          }
          else
          {
            wordHistoryOut[0] = four_lmData[w[3]].words[2];
            wordHistoryOut[1] = r2;
            wordHistoryOut[2] = -1;
          }
        }
        else
        {
          wordHistoryOut[0] = four_lmData[w[3]].words[1];
          wordHistoryOut[0] = r2;
          wordHistoryOut[2] = r;
        }
      }
      else
      {
        wordHistoryOut[0] = w[0];
        wordHistoryOut[1] = w[1];
        wordHistoryOut[2] = w[2];
      }
      #endif
    }
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method returns the unigram probability and updates the new LM history. 
/// (this method is a helper method of getP() )
/////////////////////////////////////////////////////////////////////////////////////////////////////

float LanguageModel::getUniGramP(int newWordID, int *wordHistoryNew)
{
  wordHistoryNew[0] = newWordID;
  return uni_lmData[newWordID].p;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method returns the bigram probability and updates the new LM history. 
/// If needed it will backoff to getUniGramP(). (this method is a helper method of getP() )
/////////////////////////////////////////////////////////////////////////////////////////////////////

float LanguageModel::getBiGramP(int newWordID,int *wordHistory, int *wordHistoryNew)
{
  wordHistoryNew[0] = wordHistory[0];
  wordHistoryNew[1] = newWordID;
  wordHistoryNew[1] = getBiGramIndex(wordHistoryNew);
  if(wordHistoryNew[1] < 0)
  {
    if(bi_tableLength > 0)
    {
      return (uni_lmData[wordHistory[0]].backoff + getUniGramP(newWordID,wordHistoryNew));
    }
    return getUniGramP(newWordID,wordHistoryNew);
  }
  else
  {
    return bi_lmData[wordHistoryNew[1]].p;
  }
  return SMALLEST_NEGATIVE;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method returns the trigram probability and updates the new LM history. 
/// If needed it will backoff to getBiGramP(). (this method is a helper method of getP() )
/////////////////////////////////////////////////////////////////////////////////////////////////////

float LanguageModel::getTriGramP(int newWordID,int *wordHistory, int *wordHistoryNew)
{
  #if LM_NGRAM_DEPTH > 1
  /* returns a triple history. This triple is cut back to 2 at getP... */
  wordHistoryNew[0] = wordHistory[0];
  wordHistoryNew[1] = bi_lmData[wordHistory[1]].words[1];
  wordHistoryNew[2] = newWordID;

  wordHistoryNew[2] = getTriGramIndex(wordHistoryNew);  
  if(wordHistoryNew[2] < 0)
  {
    wordHistoryNew[0] = wordHistoryNew[1];
    wordHistoryNew[1] = -1;
    wordHistoryNew[2] = -1;
    if(tri_tableLength > 0)
    {
      return (bi_lmData[wordHistory[1]].backoff + getBiGramP(newWordID,wordHistoryNew,wordHistoryNew));
    }
    return getBiGramP(newWordID,wordHistoryNew,wordHistoryNew);
  }
  else
  {
    wordHistoryNew[1] = wordHistory[1];

    return tri_lmData[wordHistoryNew[2]].p;
  }
  #else
    // Not used in this case:
    newWordID      = newWordID;
    wordHistory    = wordHistory;
    wordHistoryNew = wordHistoryNew;
  #endif
  return SMALLEST_NEGATIVE;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method returns the 4-gram probability and updates the new LM history. 
/// If needed it will backoff to getTriGramP(). (this method is a helper method of getP() )
/////////////////////////////////////////////////////////////////////////////////////////////////////

float LanguageModel::get4GramP(int newWordID,int *wordHistory, int *wordHistoryNew)
{
  #if LM_NGRAM_DEPTH > 2
  /* returns a 'quatro' ;-) history. This quatro is cut back to 3 at getP... */
  wordHistoryNew[0] = wordHistory[0];
  wordHistoryNew[1] = tri_lmData[wordHistory[2]].words[1];
  wordHistoryNew[2] = tri_lmData[wordHistory[2]].words[2];
  wordHistoryNew[3] = newWordID;
  wordHistoryNew[3] = get4GramIndex(wordHistoryNew);  
  if(wordHistoryNew[3] < 0)
  {
    wordHistoryNew[0] = wordHistoryNew[1];
    int w[2] = {wordHistoryNew[1],wordHistoryNew[2]};
    wordHistoryNew[1] = getBiGramIndex(w);
    assert(wordHistoryNew[1] >= 0);
    wordHistoryNew[2] = -1;
    wordHistoryNew[3] = -1;
    if(four_tableLength > 0)
    {
      return (tri_lmData[wordHistory[2]].backoff + getTriGramP(newWordID,wordHistoryNew,wordHistoryNew));
    }
    return getTriGramP(newWordID,wordHistoryNew,wordHistoryNew);
  }
  else
  {
    wordHistoryNew[1] = wordHistory[1];
    wordHistoryNew[2] = wordHistory[2];    
    return four_lmData[wordHistoryNew[3]].p;  
  }
  #else
    // Not used in this case:
    newWordID      = newWordID;
    wordHistory    = wordHistory;
    wordHistoryNew = wordHistoryNew;
  #endif
  return SMALLEST_NEGATIVE;
}
  
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method returns the bigram index for the bigram array for the word pair w[0] and w[1].
/// If this bigram does not exist, it will return -1.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int LanguageModel::getBiGramIndex(int *w)
{
  if(bi_Hash != NULL)
  {
    int res = bi_Hash->getIndex(w);
    if(res>=0 && compareKeys(bi_lmData[res].words,w,2))
    {
      return res;
    }
  }
  return -1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method returns the trigram index for the trigram array for the word pair w[0], w[1] and w[2].
/// If this trigram does not exist, it will return -1.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int LanguageModel::getTriGramIndex(int *w)
{
  int res = tri_Hash->getIndex(w);
  if(res>=0 && compareKeys(tri_lmData[res].words,w,3))
  {
    return res;
  }
  return -1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method returns the 4-gram index for the 4-gram array for the word pair w[0], w[1], w[2] and w[3].
/// If this trigram does not exist, it will return -1.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int LanguageModel::get4GramIndex(int *w)
{
  int res = four_Hash->getIndex(w);
  if(res>=0 && compareKeys(four_lmData[res].words,w,4))
  {
    return res;
  }
  return -1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method is only used for debugging purposes. It checks if the language model is consistant with
/// the hash tables.
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool LanguageModel::checkValidity()
{
  // First check the bigram table:  
  for(int i=0;i<bi_tableLength;i++)
  {
    int w[2] = {bi_lmData[i].words[0],bi_lmData[i].words[1]};
    int res = bi_Hash->getIndex(w);
    if(res < 0)
    {
      PRINTF2("No return from %d\n",i);
    }
    else if(res != i)
    {  
      PRINTF5("ERROR IN: %d\t%s\t%s\t%f\n",i,(char*) bi_lmData[i].words[0],
                                          (char*) bi_lmData[i].words[1],
                                          bi_lmData[i].p);
      return false;
    }
  }
  // Same for the trigram:
  for(int i=0;i<tri_tableLength;i++)
  {
    int w[3] = {tri_lmData[i].words[0],tri_lmData[i].words[1],tri_lmData[i].words[2]};
    int res = tri_Hash->getIndex(w);
    if(res < 0)
    {
      PRINTF2("No return from %d\n",i);
    }
    else if(res != i)
    {  
      PRINTF6("ERROR IN: %d\t%s\t%s\t%s\t%f\n",i,(char*) tri_lmData[i].words[0],
                                                (char*) tri_lmData[i].words[1],
                                                (char*) tri_lmData[i].words[2],
                                                tri_lmData[i].p);
      return false;
    }
  }

  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method prints some information about the language model to standard output.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LanguageModel::printInfo()
{
  printf("1-gram table: %d\n", uni_tableLength);
  printf("2-gram table: %d\n", bi_tableLength);
  printf("3-gram table: %d\n", tri_tableLength);
  printf("4-gram table: %d\n", four_tableLength);
}





