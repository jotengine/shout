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
#include "shout_lm2bin.h"
#include "lexicaltree.h"
#include "trainhash.h"
#include "shout-misc.h"
#include "shoutconfig.h"

using namespace StringFunctions;
using namespace WriteFileLittleBigEndian;


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#define WARNING(a)             printf("Warning: %s\n",a); return;

#define MAPbi_lmData(i)    ((LMEntryType_2*)memMap->map(i))
#define MAPtri_lmData(i)   ((LMEntryType_3*)memMap->map(i))
#define MAPfour_lmData(i)  ((LMEntryType_4*)memMap->map(i))

int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_LM2BIN, argc, argv);
  Shout_lm2bin myLm2Bin(myConfig.getStringValue(ARGUMENT_DCT_BIN),myConfig.getStringValue(ARGUMENT_LM_ARPA),myConfig.getStringValue(ARGUMENT_LM_BIN));
  return 0;
}

void Shout_lm2bin::sortWords(int numberOfItems, int *wordList, int *sortList)
{
  int nrWords[numberOfWords];
  for(int i=0;i<numberOfWords;i++)
  {
    nrWords[i] = 0;
  }
  for(int i=0;i<numberOfItems;i++)
  {
    nrWords[wordList[i]]++;
  }
  int comm = 0;
  for(int i=0;i<numberOfWords;i++)
  {
    int oldComm = comm;
    comm += nrWords[i];
    nrWords[i]   = oldComm;
  }
  for(int i=0;i<numberOfItems;i++)
  {
    sortList[nrWords[wordList[i]]] = i;
    nrWords[wordList[i]]++;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The main function of the shout_lm2bin application will call this constructor. All work is done
/// by this constructor.
///
/// The lexical tree is needed to obtain the wordID for each word in the ARPA model. Looking up the 
/// IDs is done by just calling the getWordID() method from the LexicalTree object.
///
/// After loading the lexical tree, the ARPA model is loaded and the uni-, bi- and tri-gram probabilities
/// and backoff values are loaded and translated into LanguageModel format.
///
/// Finally, the binary language model is stored to file in the format so that it can be loaded by 
/// LanguageModel.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Shout_lm2bin::Shout_lm2bin(char *lexName, char *lmName, char *binLmName): LanguageModel()
{
  memMap            = NULL;
  tri_tableLength   = 0;
  bi_tableLength    = 0;
  #if LM_NGRAM_DEPTH > 2
  four_tableLength  = 0;
  #endif

  struct timeval time1;
  struct timeval time2;
  gettimeofday( &time1, NULL );

  FILE *lexFile = fopen(lexName,"rb");
  if(lexFile == NULL)
  {
    USER_ERROR("The lexical tree file could not be opened!");
  }
  freadEndianSafe(&numberOfWords,1,sizeof(numberOfWords),lexFile);
  freadEndianSafe(&numberOfWords,1,sizeof(numberOfWords),lexFile);
  vocabulary = new char*[numberOfWords];
  for(int i=0;i<numberOfWords;i++)
  {
    int len;
    freadEndianSafe(&len,1,sizeof(len),lexFile);
    vocabulary[i] = new char[len+1];
    freadEndianSafe(vocabulary[i],len,1,lexFile);
    vocabulary[i][len] = 0;
  }
  fclose(lexFile);
  wordTree = new StringLookup(numberOfWords,vocabulary);
  
  FILE *lmFile = fopen64(lmName,"rt");
  if(lmFile == NULL)
  {
    USER_ERROR("The language model file could not be opened!");
  }

  // Now store the stuff!
  printf("\nCreating output file...\n");
  FILE *binlmFile = fopen(binLmName,"wb");
  int depth = LM_NGRAM_DEPTH;
  fwriteEndianSafe(&(depth),1,sizeof(depth),binlmFile); 


  char   str[MAXSTR];
  char   str2[MAXSTR];
  char   str3[MAXSTR];
  char   str4[MAXSTR];
  char   str5[MAXSTR];
  float  p, backoff;
  int    wordID, wordID2, wordID3, wordID4;
  int state = 0;
  int counter = 0;
  int displayCounter = 0;

  while(state < 5 && !feof(lmFile))
  {
    char *resChar = fgets(str,MAXSTR,lmFile);
    assert(resChar != NULL);
    switch(state)
    {
      case 0:    // Expecting the start of the uni-gram list:
        if(strncmp("ngram 2=",str,8)==0)
        {
            strcpy(str2,&str[8]);
            bi_tableLength = atoi(str2);
//            fwriteEndianSafe(&(bi_tableLength),1,sizeof(bi_tableLength),binlmFile);
            bi_lmData      = NULL;
        }
        else if(strncmp("ngram 3=",str,8)==0)
        {
          #if LM_NGRAM_DEPTH > 1
            strcpy(str2,&str[8]);
            tri_tableLength = atoi(str2);
//            fwriteEndianSafe(&(tri_tableLength),1,sizeof(tri_tableLength),binlmFile);
            tri_lmData      = NULL;
          #endif
        }
        else if(strncmp("ngram 4=",str,8)==0)
        {
          #if LM_NGRAM_DEPTH > 2
            strcpy(str2,&str[8]);
            four_tableLength = atoi(str2);
//            fwriteEndianSafe(&(four_tableLength),1,sizeof(four_tableLength),binlmFile);
            four_lmData      = NULL;
          #endif
        }
        else if(strncmp("\\1-grams:",str,9)==0)
        {
          state++;
          uni_tableLength = numberOfWords;
          fwriteEndianSafe(&(uni_tableLength),1,sizeof(uni_tableLength),binlmFile);
          uni_lmData      = new LMEntryType_1[uni_tableLength];
          for(int i=0;i<uni_tableLength;i++)
          {
            uni_lmData[i].p         = SMALLEST_NEGATIVE;
            uni_lmData[i].backoff   = SMALLEST_NEGATIVE;
          }
        }
        break;
      case 1:    // Reading unigram stats:
        if(strncmp("\\2-grams:",str,9)==0)
        {
          for(int i=0;i<uni_tableLength;i++)
          {
            if(uni_lmData[i].p == SMALLEST_NEGATIVE)
            {
              printf("Warning: word \"%s\" has no unigram statistics! This word will be ignored...\n",vocabulary[i]);
            }
          }
          printf("------------------------------------------------------------ \n");
          printf("Number of unigrams : %10d\n",uni_tableLength);
          counter = 0;
          printf("Number of bigrams  : %10d",counter);
          fflush(stdout);
          state++;

          for(int ientry=0;ientry<uni_tableLength;ientry++)
          {
            fwriteEndianSafe(&(uni_lmData[ientry].p),      1,sizeof(uni_lmData[ientry].p),binlmFile);
            fwriteEndianSafe(&(uni_lmData[ientry].backoff),1,sizeof(uni_lmData[ientry].backoff),binlmFile);
          }

          if(memMap != NULL)
          {
            delete memMap;
          }
          memMap = new MemMappedFile(sizeof(LMEntryType_2),bi_tableLength,sizeof(LMEntryType_2)*bi_tableLength, "tmp.lm");
          break;
        }
        splitList(str,str2);
        p = atof(str);
        splitList(str2,str);
        if(strlen(str)>0)
        {
          backoff = atof(str);
        }
        else
        {
          backoff = 0.0;        
        }
        if(p == -99.0)
        {
          p = -9e30; // Small, but not to be ignored!
        }
        wordID = wordTree->getWordID(str2);
        if(wordID >= 0)
        {
          if(wordID >= uni_tableLength)
          {
            USER_ERROR2("The lexical tree contains at least one word that does not exist in the language model file:",vocabulary[wordID]);
          }
          uni_lmData[wordID].p         = p;
          uni_lmData[wordID].backoff   = backoff;
        }
        else
        { 
          if(strlen(str2)>0)
          {
            printf("Warning: word \"%s\" does not exist in the lexical tree. This word will be ignored...\n",str2);
          }
        } 
        break;
      case 2:      // Reading bigram stats:
        if((strncmp("\\3-grams:",str,9)==0) || (strncmp("\\end\\",str,5)==0))
        {
          printf("\b\b\b\b\b\b\b\b\b\b%10d\n",counter);
          fflush(stdout);
          bi_tableLength = counter;  // We may have won some entries...

/*
          // Sorting:
          printf("Sorting the bigrams...\n"); fflush(stdout);
          int *sortList  = new int[bi_tableLength];
          assert(sortList  != NULL);
          int *sortList2 = new int[bi_tableLength];
          assert(sortList2 != NULL);
          int *sortList3 = new int[bi_tableLength];
          assert(sortList3 != NULL);
          int *wordList  = new int[bi_tableLength];
          assert(wordList != NULL);

          bi_lmDataList  = new LMListSearch_2[numberOfWords];
          assert(bi_lmDataList != NULL);

          for(int i=0;i<numberOfWords;i++)
          {
            bi_lmDataList[i].lowestIndex = 0;
            bi_lmDataList[i].length      = 0;
          }
          for(int i=0;i<bi_tableLength;i++)
          {
            sortList[i]     = 0;
            sortList2[i]    = 0;
            sortList3[i]    = 0;
            wordList[i]     = bi_lmData[i].words[0];
          }
          sortWords(bi_tableLength,wordList,sortList);
          int index = 0;
          int sum   = 0;
          bi_lmDataList[0].lowestIndex = 0;
          for(int i=1;i<bi_tableLength;i++)
          {
            sum++;
            int diff = wordList[sortList[i]] - wordList[sortList[i-1]];
            while(diff > 0)
            {
              bi_lmDataList[index].length      = sum - bi_lmDataList[index].lowestIndex;
              index++;
              bi_lmDataList[index].lowestIndex = sum;
              diff--;
            }
          }
          bi_lmDataList[numberOfWords-1].length = bi_tableLength - bi_lmDataList[numberOfWords-1].lowestIndex;
          for(int i=0;i<numberOfWords;i++)
          {
            int start  = bi_lmDataList[i].lowestIndex;
            int length = bi_lmDataList[i].length;
            if(length > 0)
            {
              for(int i2=0;i2<length;i2++)
              {
                sortList2[i2]=-1;
                sortList3[i2]=-1;
              }
              for(int i2=0;i2<length;i2++)
              {
                wordList[i2] = bi_lmData[sortList[i2+start]].words[1];
              }
              sortWords(length,wordList,sortList2);
              for(int i2=0;i2<length;i2++)
              {
                sortList3[i2] = sortList[sortList2[i2]+start];
              }
              for(int i2=0;i2<length;i2++)
              {
                sortList[i2+start] = sortList3[i2];
              }
            }
          }

          LMEntryType_2 *reOrderBi = new LMEntryType_2[bi_tableLength];
          for(int i=0;i<bi_tableLength;i++)
          {
            reOrderBi[i].words[0] = bi_lmData[sortList[i]].words[0];
            reOrderBi[i].words[1] = bi_lmData[sortList[i]].words[1];
            reOrderBi[i].p        = bi_lmData[sortList[i]].p;
            #if LM_NGRAM_DEPTH > 1
            reOrderBi[i].backoff  = bi_lmData[sortList[i]].backoff;
            #endif
          }
          delete[] sortList;
          delete[] sortList2;
          delete[] sortList3;
          delete[] wordList;
          delete[] bi_lmData;
          bi_lmData = reOrderBi;
          reOrderBi = NULL;
*/

          // Now create the bigram hash table!  
          TrainHash *hash = NULL;
          bool done = false;
          printf("Bigram hash    (f) : %10d",0);
          fflush(stdout);
          while(!done)
          {
            if(hash != NULL)
            {
              delete hash;
            }
            hash    = new TrainHash(bi_tableLength,2,numberOfWords);

            memMap->setWritePermission(false);
            int factorCounter = 0;
            while(!done && factorCounter < 999)
            {
              printf("\b\b\b%03d",factorCounter);
              fflush(stdout);
              factorCounter++;
              hash->initialiseMapping();
              int i = 0;
              done = true;
              while(i<bi_tableLength && done)
              {
                done = hash->fillMapping(MAPbi_lmData(i)->words);
                i++;
              }
            }
          }
          hash->finalizeHash();
          for(int i=0;i<bi_tableLength;i++)
          {
            int hashResult = hash->getIndex((int*)MAPbi_lmData(i));
            if(hashResult != i)
            {
              printf("\nError in hash: (%d,%d) - %d  %d\n",i,hashResult,MAPbi_lmData(i)->words[0],MAPbi_lmData(i)->words[1]);
              exit(1);
            }
          }
          // Done with the bigram hash table.
          fwriteEndianSafe(&(bi_tableLength),1,sizeof(bi_tableLength),binlmFile);
          for(int ientry=0;ientry<bi_tableLength;ientry++)
          {
            fwriteEndianSafe(MAPbi_lmData(ientry)->words,       2,sizeof(MAPbi_lmData(ientry)->words[0]),binlmFile);
            fwriteEndianSafe(&(MAPbi_lmData(ientry)->p),        1,sizeof(MAPbi_lmData(ientry)->p),binlmFile);
            #if LM_NGRAM_DEPTH > 1
            fwriteEndianSafe(&(MAPbi_lmData(ientry)->backoff),  1,sizeof(MAPbi_lmData(ientry)->backoff),binlmFile);
            #endif
          }

/*
          fwriteEndianSafe(&(numberOfWords),1,sizeof(numberOfWords),binlmFile);
          for(int ientry=0;ientry<numberOfWords;ientry++)
          {
            fwriteEndianSafe(&(bi_lmDataList[ientry].lowestIndex), 1,sizeof(bi_lmDataList[ientry].lowestIndex),binlmFile);
            fwriteEndianSafe(&(bi_lmDataList[ientry].length),      1,sizeof(bi_lmDataList[ientry].length),binlmFile);
          }
*/

          hash->storeHash(binlmFile);
          delete hash;
          bi_Hash        = NULL;
          hash           = NULL;

          if(memMap != NULL)
          {
            delete memMap;
          }
          memMap = new MemMappedFile(sizeof(LMEntryType_3),tri_tableLength,4096*1000, "tmp.lm");

          counter        = 0;
          displayCounter = 0;
          if(strncmp("\\end\\",str,5)==0)
          {
            state=100; // Done.
          }
          else
          {
            #if LM_NGRAM_DEPTH > 1
              printf("\nNumber of trigrams : %10d",counter);
              fflush(stdout);
            #else
              printf("\n");
              fflush(stdout);
              state = 100;  // Done..
            #endif          
            state++;
          }
          break;
        }

        splitList(str,str2);
        p = atof(str);
        if(p == -99.0)
        {
          p = -9e30; // Small, but not to be ignored!
        }
        splitList(str2,str);
        splitList(str,str3);
        RightTrim(str3);
        #if LM_NGRAM_DEPTH > 1
        if(strlen(str3)>0)
        {
          backoff = atof(str3);
        }
        else
        {
          backoff = 0.0;        
        }
        #endif
        wordID  = wordTree->getWordID(str2);
        wordID2 = wordTree->getWordID(str);
        if(wordID >=0 && wordID2 >= 0)
        {
          if(counter >= bi_tableLength)
          {
            USER_ERROR("The arpa LM is invalid: more bi-grams than noted at the top of the file!");
          }
          MAPbi_lmData(counter)->words[0] = wordID;
          MAPbi_lmData(counter)->words[1] = wordID2;
          MAPbi_lmData(counter)->p       = p;
          #if LM_NGRAM_DEPTH > 1
          MAPbi_lmData(counter)->backoff = backoff;
          #endif
          counter++;
          displayCounter++;
          if(displayCounter > 9999)
          {
            displayCounter = 0;
            printf("\b\b\b\b\b\b\b\b\b\b%10d",counter);
            fflush(stdout);
          }
        }        
        break;
      case 3:      // Reading trigram stats:
        if((strncmp("\\4-grams:",str,9)==0) || (strncmp("\\end\\",str,5)==0))
        {
          printf("\b\b\b\b\b\b\b\b\b\b%10d",counter);
          fflush(stdout);
          tri_tableLength = counter;  // We may have won some entries...


/*
          // Sorting:
          printf("\nSorting the trigrams...\n"); fflush(stdout);
          int *sortList  = new int[tri_tableLength];
          int *sortList2 = new int[tri_tableLength];
          int *sortList3 = new int[tri_tableLength];
          int *wordList  = new int[tri_tableLength];
          int sortHelp[numberOfWords+1];
          for(int i=0;i<tri_tableLength;i++)
          {
            sortList[i]     = -1;
            sortList2[i]    =  0;
            sortList3[i]    =  0;
            wordList[i]     = tri_lmData[i].words[0];
          }
          sortWords(tri_tableLength,wordList,sortList);
          for(int i=0;i<tri_tableLength;i++)
          {
            assert(sortList[i] >= 0);
          }
          int index = 0;
          int sum   = 0;
          sortHelp[0] = 0;
          for(int i=1;i<tri_tableLength;i++)
          {
            sum++;
            int diff = wordList[sortList[i]] - wordList[sortList[i-1]];
            while(diff > 0)
            {
              index++;
              sortHelp[index] = sum;
              diff--;
            }
          }

          for(int i=index+1;i<numberOfWords+1;i++)
          {
            sortHelp[i] = tri_tableLength;
          }
          for(int i=0;i<numberOfWords;i++)
          {
            int start  = sortHelp[i];
            int length = sortHelp[i+1] - start;
            if(length > 0)
            {
              for(int i2=0;i2<length;i2++)
              {
                wordList[i2] = tri_lmData[sortList[i2+start]].words[1];
              }
              sortWords(length,wordList,sortList2);
              for(int i2=0;i2<length;i2++)
              {
                sortList3[i2] = sortList[sortList2[i2]+start];
              }
              for(int i2=0;i2<length;i2++)
              {
                sortList[i2+start] = sortList3[i2];
              }
            }
          }


          LMEntryType_3 *reOrderTri = new LMEntryType_3[tri_tableLength];
          for(int i=0;i<tri_tableLength;i++)
          {
            reOrderTri[i].words[0] = tri_lmData[sortList[i]].words[0];
            reOrderTri[i].words[1] = tri_lmData[sortList[i]].words[1];
            reOrderTri[i].words[2] = tri_lmData[sortList[i]].words[2];
            reOrderTri[i].p        = tri_lmData[sortList[i]].p;
            #if LM_NGRAM_DEPTH > 2
            reOrderTri[i].backoff  = tri_lmData[sortList[i]].backoff;
            #endif
          }
          delete[] tri_lmData;
          tri_lmData = reOrderTri;
          reOrderTri = NULL;


          tri_tableLengthList = 1;
          int prevID1         = tri_lmData[0].words[0];
          int prevID2         = tri_lmData[0].words[1];
          for(int li=1;li<tri_tableLength;li++)
          {
            if(prevID1 != tri_lmData[li].words[0] || prevID2 != tri_lmData[li].words[1])
            {
              tri_tableLengthList++;
              prevID1 = tri_lmData[li].words[0];
              prevID2 = tri_lmData[li].words[1];
            }
          }
          printf("Number of trigram list-search items: %d\n",tri_tableLengthList); fflush(stdout);
          tri_lmDataList = new LMListSearch_3[tri_tableLengthList];

          tri_lmDataList[0].words[0]    = tri_lmData[0].words[0];
          tri_lmDataList[0].words[1]    = tri_lmData[0].words[1];
          tri_lmDataList[0].lowestIndex = 0;
          tri_lmDataList[0].length      = 1;
          prevID1                       = tri_lmData[0].words[0];
          prevID2                       = tri_lmData[0].words[1];
          int listL                     = 0;
          for(int li=0;li<tri_tableLength;li++)
          {
            if(prevID1 != tri_lmData[li].words[0] || prevID2 != tri_lmData[li].words[1])
            {
              listL++;
              tri_lmDataList[listL].words[0]    = tri_lmData[li].words[0];
              tri_lmDataList[listL].words[1]    = tri_lmData[li].words[1];
              tri_lmDataList[listL].lowestIndex = li;
              tri_lmDataList[listL].length      = 1;
              prevID1 = tri_lmData[li].words[0];
              prevID2 = tri_lmData[li].words[1];
            }
            else
            {
              tri_lmDataList[listL].length++;
            }
          }
          listL=0;
          for(int li=0;li<tri_tableLengthList;li++)
          {
            listL += tri_lmDataList[li].length;
          }
          printf("Number of trigrams in sorted lists: %d\n",listL);
          delete[] sortList;
          delete[] sortList2;
          delete[] sortList3;
          delete[] wordList;
*/

          // Create the 3-gram hash:
          #if LM_NGRAM_DEPTH > 1
            TrainHash *hash = NULL;
            bool done = false;
            printf("\nTrigram hash   (f) : %10d",0);
            fflush(stdout);
            while(!done)
            {
              if(hash != NULL)
              {
                delete hash;
              }
              hash    = new TrainHash(tri_tableLength,3,numberOfWords);

              memMap->setWritePermission(false);
              int factorCounter = 0;
              while(!done && factorCounter < 999)
              {
                printf("\b\b\b%03d",factorCounter);
                fflush(stdout);
                factorCounter++;
                hash->initialiseMapping(); 
                int i = 0;
                done = true;
                while(i<tri_tableLength && done)
                {
                  done = hash->fillMapping(MAPtri_lmData(i)->words);
                  i++;
                }
              }
            }
            hash->finalizeHash();
            for(int i=0;i<tri_tableLength;i++)
            {
              int hashResult = hash->getIndex((int*)MAPtri_lmData(i));
              if(hashResult != i)
              {
                printf("Error in hash: (%d,%d)- %d %d %d\n",
                                    i,hashResult,MAPtri_lmData(i)->words[0],MAPtri_lmData(i)->words[1],MAPtri_lmData(i)->words[2]);
                exit(1);
              }
            }

            fwriteEndianSafe(&(tri_tableLength),1,sizeof(tri_tableLength),binlmFile);

            if(tri_tableLength > 0)
            {
              for(int ientry=0;ientry<tri_tableLength;ientry++)
              {
                fwriteEndianSafe(MAPtri_lmData(ientry)->words,       3,sizeof(MAPtri_lmData(ientry)->words[0]),binlmFile);
                fwriteEndianSafe(&(MAPtri_lmData(ientry)->p),        1,sizeof(MAPtri_lmData(ientry)->p),binlmFile);
                #if LM_NGRAM_DEPTH > 2
                fwriteEndianSafe(&(MAPtri_lmData(ientry)->backoff),  1,sizeof(MAPtri_lmData(ientry)->backoff),binlmFile);
                #endif
              }
              hash->storeHash(binlmFile);
            }
            delete hash;
            tri_Hash = NULL;
            hash = NULL;

            if(tri_tableLength > 0)
            {
/*
              // Now create the trigram-list hash table!  
              hash = NULL;
              done = false;
              printf("Trigram-list hash    (f) : %10d",0);
              fflush(stdout);
              while(!done)
              {
                if(hash != NULL)
                {
                  delete hash;
                }
                hash    = new TrainHash(tri_tableLengthList,2,numberOfWords);
                tri_HashListSearch = hash;
    
                int factorCounter = 0;
                while(!done && factorCounter < 999)
                {
                  printf("\b\b\b%03d",factorCounter);
                  fflush(stdout);
                  factorCounter++;
                  hash->initialiseMapping();
                  int i = 0;
                  done = true;
                  while(i<tri_tableLengthList && done)
                  {
                    done = hash->fillMapping(tri_lmDataList[i].words);
                    i++;
                  }
                }
              }
              hash->finalizeHash();
              for(int i=0;i<tri_tableLengthList;i++)
              {
                int hashResult = hash->getIndex((int*)&tri_lmDataList[i]);
                if(hashResult != i)
                {
                  printf("\nError in hash: (%d,%d) - %d\n",i,hashResult,tri_lmDataList[i].words[0]);
                  exit(1);
                }
              }
              // Sorting information:
              fwriteEndianSafe(&(tri_tableLengthList),1,sizeof(tri_tableLengthList),binlmFile);
              for(int ientry=0;ientry<tri_tableLengthList;ientry++)
              {
                fwriteEndianSafe(tri_lmDataList[ientry].words,          2,sizeof(tri_lmDataList[ientry].words[0]),binlmFile);
                fwriteEndianSafe(&(tri_lmDataList[ientry].lowestIndex), 1,sizeof(tri_lmDataList[ientry].lowestIndex),binlmFile);
                fwriteEndianSafe(&(tri_lmDataList[ientry].length),      1,sizeof(tri_lmDataList[ientry].length),binlmFile);
              }
              // Done with the trigram hash table.
              hash->storeHash(binlmFile);
              delete hash;
              tri_HashListSearch    = NULL;
              hash                  = NULL;
            */

            }

            if(memMap != NULL)
            {
              delete memMap;
            }
            memMap = new MemMappedFile(sizeof(LMEntryType_4),four_tableLength,4096*1000, "tmp.lm");

          #endif




          // Done with the 3-gram hash.

          counter         = 0;
          displayCounter  = 0;

          if(strncmp("\\end\\",str,5)==0)
          {
            state=100; // Done.
          }
          else
          {
            #if LM_NGRAM_DEPTH > 2
              printf("\nNumber of 4-grams  : %10d",counter);
              fflush(stdout);
            #else
              printf("\n");
              fflush(stdout);
              state = 100;  // Done..
            #endif
            state++;
          }
          break;
        }

        splitList(str,str2);
        p = atof(str);
        if(p == -99.0)
        {
          p = -9e30; // Small, but not to be ignored!
        }
        splitList(str2,str);
        splitList(str,str3);
        splitList(str3,str4);
        RightTrim(str4);
        #if LM_NGRAM_DEPTH > 2        
        if(strlen(str4)>0)
        {
          backoff = atof(str4);
        }
        else
        {
          backoff = 0.0;        
        }
        #endif
        wordID  = wordTree->getWordID(str2);
        wordID2 = wordTree->getWordID(str);
        wordID3 = wordTree->getWordID(str3);

        if(wordID >=0 && wordID2 >= 0 && wordID3 >=0)
        {
          if(counter >= tri_tableLength)
          {
            USER_ERROR("The arpa LM is invalid: more tri-grams than noted at the top of the file!");
          }
          MAPtri_lmData(counter)->words[0] = wordID;
          MAPtri_lmData(counter)->words[1] = wordID2;
          MAPtri_lmData(counter)->words[2] = wordID3;
          MAPtri_lmData(counter)->p        = p;
          #if LM_NGRAM_DEPTH > 2
          MAPtri_lmData(counter)->backoff  = backoff;
          #endif
          counter++;
          displayCounter++;
          if(displayCounter > 9999)
          {
            displayCounter = 0;
            printf("\b\b\b\b\b\b\b\b\b\b%10d",counter);
            fflush(stdout);
          }
        }
        break;
      case 4:      // Reading 4-gram stats:
        if(strncmp("\\end\\",str,5)==0)
        {
          printf("\b\b\b\b\b\b\b\b\b\b%10d\n",counter);
          fflush(stdout);
          four_tableLength = counter;  // We may have won some entries...

          // Create the 4-gram hash:
          #if LM_NGRAM_DEPTH > 2
            TrainHash *hash = NULL;
            bool done = false;
            printf("\n4-gram hash    (f) : %10d",0);
            fflush(stdout);
            while(!done)
            {
              if(hash != NULL)
              {
                delete hash;
              }
              hash    = new TrainHash(four_tableLength,4,numberOfWords);

              memMap->setWritePermission(false);
              int factorCounter = 0;
              while(!done && factorCounter < 999)
              {
                printf("\b\b\b%03d",factorCounter);
                fflush(stdout);
                factorCounter++;
                hash->initialiseMapping(); 
                int i = 0;
                done = true;
                while(i<four_tableLength && done)
                {
                  done = hash->fillMapping(i,MAPfour_lmData(i)->words);
                  i++;
                }       
              }
            }
            hash->finalizeHash();
            for(int i=0;i<four_tableLength;i++)
            {
              int hashResult = hash->getIndex((int*)MAPfour_lmData(i));
              if(hashResult != i)
              {
                printf("Error in hash: (%d,%d)- %d %d %d %d\n", i,hashResult,
                                    MAPfour_lmData(i)->words[0],MAPfour_lmData(i)->words[1],
                                    MAPfour_lmData(i)->words[2],MAPfour_lmData(i)->words[3]);
                exit(1);
              }
            }
            fwriteEndianSafe(&(four_tableLength),1,sizeof(four_tableLength),binlmFile);
            if(four_tableLength > 0)
            {
              for(int ientry=0;ientry<four_tableLength;ientry++)
              {
                fwriteEndianSafe(MAPfour_lmData(ientry)->words,       4,sizeof(MAPfour_lmData(ientry)->words[0]),binlmFile);
                fwriteEndianSafe(&(MAPfour_lmData(ientry)->p),        1,sizeof(MAPfour_lmData(ientry)->p),binlmFile);
                #if LM_NGRAM_DEPTH > 3
                fwriteEndianSafe(&(MAPfour_lmData(ientry)->backoff),  1,sizeof(MAPfour_lmData(ientry)->backoff),binlmFile);
                #endif
              }
              hash->storeHash(binlmFile);
            }
            delete hash;
            four_Hash = NULL;
            hash = NULL;
          #endif

          // Done.

          state++;
          break;
        }
        splitList(str,str2);
        p = atof(str);
        if(p == -99.0)
        {
          p = -9e30; // Small, but not to be ignored!
        }
        splitList(str2,str);
        splitList(str,str3);
        splitList(str3,str4);
        splitList(str4,str5);
        RightTrim(str5);
        #if LM_NGRAM_DEPTH > 3
/*        if(strlen(str5)>0)
        {
          backoff = atof(str5);
        }
        else
        {
          backoff = 0.0;        
        }*/
        #endif
        wordID  = wordTree->getWordID(str2);
        wordID2 = wordTree->getWordID(str);
        wordID3 = wordTree->getWordID(str3);
        wordID4 = wordTree->getWordID(str4);
        
        if(wordID >=0 && wordID2 >= 0 && wordID3 >=0 && wordID4 >= 0)
        {         
          if(counter >= four_tableLength)
          {
            USER_ERROR("The arpa LM is invalid: more 4-grams than noted at the top of the file!");
          }
          MAPfour_lmData(counter)->words[0] = wordID;
          MAPfour_lmData(counter)->words[1] = wordID2;
          MAPfour_lmData(counter)->words[2] = wordID3;
          MAPfour_lmData(counter)->words[3] = wordID4;
          MAPfour_lmData(counter)->p       = p;
          #if LM_NGRAM_DEPTH > 3
//          MAPfour_lmData(counter)->backoff = backoff;
          #endif
          counter++;
          displayCounter++;
          if(displayCounter > 9999)
          {
            displayCounter = 0;
            printf("\b\b\b\b\b\b\b\b\b\b%10d",counter);
            fflush(stdout);
          }
        }                
        break;        
      default:
        break;
    }  
  }

  if(tri_tableLength == 0)
  {
    fwriteEndianSafe(&(tri_tableLength),1,sizeof(tri_tableLength),binlmFile);
  }

  fclose(binlmFile);

  gettimeofday( &time2, NULL );
  int sec = (time2.tv_sec) - (time1.tv_sec);
  int min = ((int)sec/60);
  sec = sec - min*60;
  printf("\n----------------------------------------------------\n");
  printf("-- Done. It took %d minutes and %d seconds.\n", min,sec);
  printf("----------------------------------------------------\n");    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor is empty.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Shout_lm2bin::~Shout_lm2bin()
{
}



