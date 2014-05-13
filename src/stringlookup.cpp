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
#include "stringlookup.h"
#include "shout-misc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace StringFunctions;


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////////

StringLookup::~StringLookup()
{
  if(wordTree != NULL)
  {
    deleteWordTree(wordTree);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Creates a tree-structure of all words. This will make getWordID a lot faster, but it costs memory!
/////////////////////////////////////////////////////////////////////////////////////////////////////////

StringLookup::StringLookup(int numberOfWords, char **vocabulary)
{
  wordTree = NULL;
  if(wordTree == NULL)
  {
    wordTree = new WordStringNode;
    wordTree->wordID = -1;
    wordTree->c    = 0;
    wordTree->next = NULL;
    wordTree->par  = NULL;
    for(int i=0;i<numberOfWords;i++)
    {
      WordStringNode *w = wordTree;
      assert(vocabulary[i] != NULL);
      int maxChar = strlen(vocabulary[i]);
      for(int i2=0;i2<maxChar;i2++)
      {
        char letter = vocabulary[i][i2];
        if(w->next == NULL)
        {
          w->next = new WordStringNode;
          w = w->next;
          w->c    = letter;
          w->next = NULL;
          w->par  = NULL;
          w->wordID = -1;
        }
        else
        {
          WordStringNode *wpar = w->next;
          while(wpar->par != NULL && wpar->c != letter)
          {
            wpar = wpar->par;
          }
          if(wpar->c == letter)
          {
            w = wpar; 
          }
          else
          {
            wpar->par = new WordStringNode;
            w = wpar->par;
            w->c    = letter;
            w->next = NULL;
            w->par  = NULL;            
            w->wordID = -1;
          }
        }
      }
      w->wordID = i;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the ID of the word 'word'. If the word is not in the vocabulary (OOV), -1 is returned.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int StringLookup::getWordID(char *word)
{  
  if(strlen(word) == 0)
  {
    return -1;
  } 
  if(wordTree != NULL)
  {
    int i=0;
    WordStringNode *w = wordTree;
    int wordLength = strlen(word);
    while((i < wordLength) && w != NULL)
    {
      char letter = word[i];
      WordStringNode *wpar = w->next;
      if(wpar == NULL)
      {
        return -1;
      }
      while(wpar->par != NULL && wpar->c != letter)
      {
        wpar = wpar->par;
      }
      if(wpar->c == letter)
      {
        w = wpar; 
      }
      else
      {
        return -1;
      }
      i++;
    }
    if(w!=NULL)
    {
      return w->wordID;
    }  
  }
  return -1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Deletes the word tree.
/////////////////////////////////////////////////////////////////////////////////////////////////////////


void StringLookup::deleteWordTree(WordStringNode *w)
{
  if(w != NULL)
  {
    deleteWordTree(w->par);
    deleteWordTree(w->next);
    delete w;
  }
}
