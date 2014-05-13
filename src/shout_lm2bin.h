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
#ifndef SHOUT_LM2BIN_H
#define SHOUT_LM2BIN_H

#include "languagemodel.h"
#include "lexicaltree.h"
#include "stringlookup.h"
#include "memmappedfile.h"
#include "trainhash.h"

#include <stdio.h>

////////////////////////////////////////////////////////////////////////////
/// \brief With this class (main class of the application shout_lm2bin) an
/// ARPA language model can be stored in a LanguageModel object and written
/// to disk.
////////////////////////////////////////////////////////////////////////////

class Shout_lm2bin  : public LanguageModel
{
  protected:
    MemMappedFile  *memMap;
    int             numberOfWords;
    char          **vocabulary;
    StringLookup   *wordTree;
  protected:
    int  getNumberOfValidEntries(char *word1, char *word2, FILE *file);
    void sortWords(int numberOfItems, int *wordList, int *sortList);
  
  public:
    Shout_lm2bin(char *lexName, char *lmName, char *binLmName);
   ~Shout_lm2bin();
};

#endif
