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
#include <sys/time.h>
#include "standard.h"
#include "shoutconfig.h"
#include "shout_dct2lextree.h"
#include "phonefilereader.h"


#define WARNING(a)             printf("Warning: %s\n",a); return;


int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_DCT2LEXTREE, argc, argv);

  Shout_dct2lextree myLex2Dct(myConfig.getStringValue(ARGUMENT_PHONE_LIST),myConfig.getStringValue(ARGUMENT_DCT_TEXT),myConfig.getStringValue(ARGUMENT_DCT_BIN));
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor is called by the main function of shout_dct2lextree. It loads a phone definition
/// file and a dictionary. After creating a PPT, it is stored to file in a binary format so that it
/// can be read by LexicalTree.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Shout_dct2lextree::Shout_dct2lextree(char *phoneListFile, char *dctFileName, char *lexFileName): LexicalTree(stdout)
{
  DataStats *data = NULL;
  char str[MAXSTR];
  int i;
  int numberOfSil;
  
  PhoneFileReader readPhones(phoneListFile, &data, &numberOfPhones, &numberOfSil);
  
  printf("Reading dictionary...\n");

  FILE *inFile = fopen(dctFileName,"rt");
  // Count the number of words: 
  numberOfWords = 0;
  while(!feof(inFile))
  {
    if(fgets(str,MAXSTR,inFile)>0)
    {
      numberOfWords++;
    }
  }
  vocabulary = new char*[numberOfWords];

  printf("Total number of words (including pronunciation variants): %d\n",numberOfWords);

  printf("Filling lexical tree...\n");

  // Fill the lexical tree! 
  fseek(inFile,0,SEEK_SET);
  char phoneList[MAXSTR];
  numberOfWords = 0;
  LexicalNode *previousNode = NULL;  
  LexicalNode *currentNode  = NULL;  
  bool adding = false;
  int totNumOfPhones = 0;
  int totNumOfNodes  = 0;

  // Initialise tree:
  treeStart = new LexicalNode*[numberOfPhones];
  for(int a=0;a<numberOfPhones;a++)
  {
    treeStart[a] = NULL;
  }

  while(!feof(inFile))
  {
    if(fgets(str,MAXSTR,inFile)>0)
    {
      splitList(str,phoneList);
      // Check if the word is already in the list:
      int wordID = 0;
      while((wordID < numberOfWords) && (strcmp(vocabulary[wordID],str) != 0))
      {
        wordID++;
      }

      if(wordID == numberOfWords)
      {
        vocabulary[numberOfWords] = new char[strlen(str)+1];
        strncpy(vocabulary[numberOfWords],str,strlen(str));
        vocabulary[numberOfWords][strlen(str)] = 0;
        numberOfWords++;
      }
      strcpy(str,phoneList);
      adding           = false;
      int phoneIDHist  =  0;
      int phoneID1     = -1;
      int phoneID2     =  0;

      if(strlen(str)>0)
      {
        splitList(str,phoneList);     
        for(int ai=strlen(str);ai<10;ai++)
        {
          str[ai] = ' ';
        }
        str[10] = 0;
        phoneID1 = 0;
        while((phoneID1<numberOfPhones) && (strcmp(data[phoneID1].name,str)!=0))
        {
          phoneID1++;
        }
        if(phoneID1==numberOfPhones)
        {
          USER_ERROR2("Encountered an unknown phone:",str);
        }
        strcpy(str,phoneList);
        currentNode      = treeStart[phoneID1];
      }
      else
      {
        USER_ERROR2("The following word is lacking pronunciation:",vocabulary[wordID]);
      }

      while(phoneID2 >= 0)
      {
        if(strlen(str)>0)
        {
          splitList(str,phoneList);
          for(int ai=strlen(str);ai<10;ai++)
          {
            str[ai] = ' ';
          }
          str[10] = 0;
          phoneID2 = 0;
          while((phoneID2<numberOfPhones) && (strcmp(data[phoneID2].name,str)!=0))
          {
            phoneID2++;
          }
          if(phoneID2==numberOfPhones)
          {
            USER_ERROR2("Encountered an unknown phone:",str);
          }
          strcpy(str,phoneList);
        }
        else
        {
          phoneID2 = -1;
        }

        totNumOfPhones++;
        if(currentNode == NULL)
        {
          // This is the first node in the network!
          adding      = true;
          currentNode = new LexicalNode;
          currentNode->modelID      = phoneID1;
          currentNode->contextPrev  = phoneIDHist;
          currentNode->contextNext  = phoneID2;
          currentNode->next         = NULL;
          currentNode->parallel     = NULL;
          currentNode->wordID       = -1;
          treeStart[phoneID1] = currentNode;
          totNumOfNodes++;
        }
        else if(adding)
        {
          // We are currently adding new phones:
          LexicalNode *n    = new LexicalNode;
          n->modelID        = phoneID1;
          n->contextPrev    = phoneIDHist;
          n->contextNext    = phoneID2;
          n->next           = NULL;
          n->parallel       = NULL;
          n->wordID         = -1;
          totNumOfNodes++;
          currentNode->next = n;
          currentNode = n;
        }
        else
        {
          // We are browsing an existing node:
          if(!(currentNode->modelID == phoneID1 && currentNode->contextNext == phoneID2))
          {
            // Hey! Not the same phone. First check if there are parallel paths:
            while((currentNode->parallel != NULL) && 
                  !(currentNode->parallel->modelID == phoneID1 && currentNode->parallel->contextNext == phoneID2))
            {
              currentNode = currentNode->parallel;
            }
            if(currentNode->parallel == NULL)
            {
              // The correct phone was not found. Let's make our own!
              adding = true;
              LexicalNode *n = new LexicalNode;
              n->modelID            = phoneID1;
              n->contextPrev        = phoneIDHist;
              n->contextNext        = phoneID2;              
              n->next               = NULL;
              n->parallel           = NULL;
              n->wordID             = -1;
              currentNode->parallel = n;              
              currentNode = n;              
              totNumOfNodes++;
            }
            else
            {
              currentNode = currentNode->parallel;
              assert(currentNode->contextPrev == phoneIDHist);
            }                      
          }
          if(!adding)  // Now we know that the current node has modelID i! (either because we went to a parallel path or directly)
          {
            // Still the same phone sequence, let's see what the new 'currentNode' should be:
            if(currentNode->next == NULL)
            {
              // Bingo, lets extent this path, but make this node an end-node:
              adding = true;
            }
            else
            {
              previousNode = currentNode;
              currentNode  = currentNode->next;
            }
          }
        }

        if(phoneID2 != -1)
        {
          phoneIDHist = phoneID1;
          phoneID1    = phoneID2;
        }
      }
      assert(currentNode != NULL);
      if(adding && currentNode->wordID < 0)
      {
        currentNode->wordID = wordID;    
      }
      else if(phoneID1 != -1)
      {
        // Wow: this word is a sub-word of some other word... 
        // Make an extra node, even though it will contain the same phone as the current node:
        if(adding)
        {
          previousNode = currentNode;
        }
        LexicalNode *n         = new LexicalNode;
        n->modelID             = phoneID1;
        n->contextPrev         = phoneIDHist;
        n->contextNext         = phoneID2;                      
        n->next                = NULL;
        n->parallel            = previousNode->parallel;
        n->wordID              = wordID;
        previousNode->parallel = n;
        totNumOfNodes++;
      }
    }
  }
    
  printf("------------------------ Done! -------------------------------\n");

  fclose(inFile);
  printf("Number of phones: %d\n",numberOfPhones);
  printf("Number of unique words : %d\n",numberOfWords);
  printf("Total number of phones in the dictionary  : %d\n",totNumOfPhones);
  printf("Total number of nodes in the lexical tree : %d\n",totNumOfNodes);

  // Store the lexicon tree 
  FILE *outFile = fopen(lexFileName,"wb");
  fwrite(&numberOfPhones,1,sizeof(numberOfPhones),outFile);
  fwrite(&numberOfWords,1,sizeof(numberOfWords),outFile);
  printf("Writing to file:\n");
  printf("-- Writing words...\n");
  for(i=0;i<numberOfWords;i++)
  {
    int len = strlen(vocabulary[i]);
    fwrite(&len,1,sizeof(len),outFile);
    fwrite(vocabulary[i],1,len,outFile);
  }
  printf("-- Writing lexical tree...\n");
  for(int a=0;a<numberOfPhones;a++)
  {
    char treeExists = 0;
    if(treeStart[a] == NULL)
    {
      treeExists = 0;
      fwrite(&treeExists,1,sizeof(treeExists),outFile);    
    }
    else
    {
      treeExists = 1;
      fwrite(&treeExists,1,sizeof(treeExists),outFile);
      writeTree(treeStart[a], outFile);
    }    
  }
  fclose(outFile);
  printf("Finished writing to file.\n");    
    
  
  // Clean up local data 
  delete[] data; 
}





/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This destructor is empty.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Shout_dct2lextree::~Shout_dct2lextree() 
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will write the tree starting with LexicalNode 'node' to a binary file so that it can
/// be read by LexicalTree.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Shout_dct2lextree::writeTree(LexicalNode *node, FILE *outFile)
{
  if(node != NULL)
  {
    char state;
    if(node->next != NULL)
    {
      state = 0;
      fwrite(&state,1,sizeof(state),outFile);      
      writeTree(node->next,outFile);
    }
    else
    {
      state = 9;
      fwrite(&state,1,sizeof(state),outFile);      
    } 
    if(node->parallel != NULL)
    {
      state = 1;
      fwrite(&state,1,sizeof(state),outFile);      
      writeTree(node->parallel,outFile);
    }
    else
    {
      state = 9;
      fwrite(&state,1,sizeof(state),outFile);      
    }     
    state = 2;
    fwrite(&state,1,sizeof(state),outFile);          
    fwrite(&(node->modelID),1,sizeof(node->modelID),outFile);
    fwrite(&(node->contextPrev),1,sizeof(node->contextPrev),outFile);
    fwrite(&(node->contextNext),1,sizeof(node->contextNext),outFile);   
    fwrite(&(node->wordID),1,sizeof(node->wordID),outFile);   
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This is a helper function which will split a string 'str' into two strings ('str' and 'str2') at
/// the first space or tab character of 'str'.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Shout_dct2lextree::splitList(char *str, char *str2)
{
  int k;
  k = 0;
  int n = strlen(str);
  while (k<n && !(str[k]==' ' || str[k]=='\t'))
  {
    k++;
  }
  int kstop = k;
  while (k<n && (str[k]==' ' || str[k]=='\t'))
  {
    k++;
  }
  strcpy(str2,&(str[k]));  
  k = strlen(str2) - 1;  
  while(k >= 0)
  {
    if(str2[k] == '\n')
    {
      str2[k] = 0;
    }
    k--;
  } 
  str[kstop] = 0;
}
