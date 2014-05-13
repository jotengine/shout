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
#ifndef LANGUAGEMODEL_H
#define LANGUAGEMODEL_H

#include "standard.h"
#include <stdio.h>
#include "hash.h"

struct LMEntryType_1;
struct LMEntryType_2;
struct LMEntryType_3;
struct LMEntryType_4;


////////////////////////////////////////////////////////////////////////////
// Lattice:
// Op de arcs:
// -------------------------------------
// woord
// nodes niet samenvoegen! (LM uniek)
// lm score (niet scalen)
// am score
// begin/eind tijd per woord
//   (offset phone tijden + score)
////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////
/// \brief This structure-type is part of the LanguageModel data structure.
///
/// The LanguageModel object contains an array with elements of this type.
/// It represents all LM unigram probabilities.
/// This array is mapped one-to-one with all wordIDs.
////////////////////////////////////////////////////////////////////////////

struct LMEntryType_1
{
  float               p;
  float               backoff;
};

////////////////////////////////////////////////////////////////////////////
/// \brief This structure-type is part of the LanguageModel data structure.
///
/// The LanguageModel object contains an array with elements of this type.
/// It represents all LM bigram probabilities.
/// The indices of this array are mapped to wordID pairs using a (perfect 
/// minimal) Hash table.
////////////////////////////////////////////////////////////////////////////

struct LMEntryType_2
{
  int                 words[2];
  float               p;
  #if LM_NGRAM_DEPTH > 1
  float               backoff;
  #endif
};

struct LMListSearch_2
{
  int                 lowestIndex;
  int                 length;
};


////////////////////////////////////////////////////////////////////////////
/// \brief This structure-type is part of the LanguageModel data structure.
///
/// The LanguageModel object contains an array with elements of this type.
/// It represents all LM trigram probabilities.
/// The indices of this array are mapped to wordID triples using a (perfect 
/// minimal) Hash table.
////////////////////////////////////////////////////////////////////////////

struct LMEntryType_3
{
  int                 words[3];
  float               p;
  #if LM_NGRAM_DEPTH > 2
  float               backoff;
  #endif
};

struct LMListSearch_3
{
  int                 words[2];
  int                 lowestIndex;
  int                 length;
};

////////////////////////////////////////////////////////////////////////////
/// \brief This structure-type is part of the LanguageModel data structure.
///
/// The LanguageModel object contains an array with elements of this type.
/// It represents all LM 4-gram probabilities.
/// The indices of this array are mapped to wordID triples using a (perfect 
/// minimal) Hash table.
////////////////////////////////////////////////////////////////////////////

struct LMEntryType_4
{
  int                 words[4];
  float               p;
  #if LM_NGRAM_DEPTH > 3
//  float               backoff;
  #endif
};

////////////////////////////////////////////////////////////////////////////
/// \brief This class contains all language model functionality.
///
/// Each LanguageModel object stores its data in three tables (uni-, bi- and tri-gram
/// tables). The bigram and trigram tables are queried using the two (perfect
/// minimal) Hash tables. The data is stored using the class Shout_lm2bin with
/// the application shout_lm2bin.
///
/// With the method getUnigram() the unigram probability of a word can be retrieved.
/// The method getP() is used to retrieve the LM probability given a certain 
/// word history. getP() will handle backoff if needed.
////////////////////////////////////////////////////////////////////////////

class LanguageModel
{
  protected:
    int  uni_tableLength;            ///< The length of the uni-gram (direct) lookup-table.
    int  bi_tableLength;             ///< The length of the bi-gram (hash) lookup-table.
    int  tri_tableLength;            ///< The length of the tri-gram (hash) lookup-table.
    int  four_tableLength;           ///< The length of the 4-gram (hash) lookup-table.
    int  uni_startData;
    int  bi_startData;
    int  tri_startData;


    Hash *bi_Hash;                   ///< The Hash-function of the bi-gram table.
    Hash *tri_Hash;                  ///< The Hash-function of the tri-gram table.
    Hash *four_Hash;                 ///< The Hash-function of the tri-gram table.

    LMEntryType_1 *uni_lmData;       ///< The uni-gram direct lookup-table.
    LMEntryType_2 *bi_lmData;        ///< The bi-gram hash table.
    LMEntryType_3 *tri_lmData;       ///< The tri-gram hash table.
    LMEntryType_4 *four_lmData;      ///< The tri-gram hash table.

    Hash *tri_HashListSearch;
    int  tri_tableLengthList;
    LMListSearch_2 *bi_lmDataList;
    LMListSearch_3 *tri_lmDataList;

  public:
    LanguageModel();
    LanguageModel(int numberOfWords);
    LanguageModel(const char *lmbinFileName);
   ~LanguageModel();
   
  protected:
    float      getUniGramP        (int newWordID, int *wordHistoryNew);
    float      getBiGramP         (int newWordID, int *wordHistory, int *wordHistoryNew);
    float      getTriGramP        (int newWordID, int *wordHistory, int *wordHistoryNew);
    float      get4GramP          (int newWordID, int *wordHistory, int *wordHistoryNew);
    int        getBiGramIndex     (int *w);
    int        getTriGramIndex    (int *w);        
    int        get4GramIndex      (int *w);        
  
  public:
    void       printInfo          ();
    void       setUnigram         (int wordID, float p);
    float      getUnigram         (int wordID);
    int        getLastWordID      (int *wordHistory);
    void       getAllP            (int *lmHistory, float *allP);
    virtual float getP            (int newWordID,int *wordHistoryIn, int *wordHistoryOut);
    bool       checkValidity      ();
    inline int getNumberOfWords   ()          const {return uni_tableLength; };
    inline bool compareKeys(int *key1, int *key2, int length) const {for(int i=0;i<length;i++){if(key1[i] != key2[i]) return false;}; return true;};
};


#endif
