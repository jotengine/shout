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
#ifndef STRINGLOOKUP_H
#define STRINGLOOKUP_H


////////////////////////////////////////////////////////////////////////////
/// \brief A WordStringNode is a node from the tree-structure that contains 
/// the text versions of all words in a vocabulary. Organising the 
/// words in a tree makes look-up faster. This is mainly important when
/// converting ARPA language models to Shout language models, but also
/// when parsing application arguments.
////////////////////////////////////////////////////////////////////////////

struct WordStringNode
{
  int wordID;
  WordStringNode *next;
  WordStringNode *par;
  char c;
};


////////////////////////////////////////////////////////////////////////////
/// \todo Docs
////////////////////////////////////////////////////////////////////////////

class StringLookup
{
protected:
  WordStringNode *wordTree;

public:
  StringLookup(int numberOfWords, char **vocabulary);
 ~StringLookup();
  
  int getWordID(char *word);

protected:
  void deleteWordTree(WordStringNode *w);
  
};

#endif
