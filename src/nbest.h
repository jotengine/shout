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
#ifndef NBEST_H
#define NBEST_H

#include "phonemodel.h"
#include "languagemodel.h"


////////////////////////////////////////////////////////////////////////////
/// \brief Used for N-Best calculations:
////////////////////////////////////////////////////////////////////////////

struct NBestType
{
  float AM;
  float LM;
  int   wordID;
};

////////////////////////////////////////////////////////////////////////////
/// \brief Used for N-Best calculations:
////////////////////////////////////////////////////////////////////////////

struct NBestList
{
  int        nrWords;
  int        noSilWords;
  NBestType *wordList;
  float      totAM;
  float      totLM;  
  bool       processed;
};

////////////////////////////////////////////////////////////////////////////
/// \brief This class will translate recognition administration (a bunch
/// of WLRType objects) into an N-Best list.
////////////////////////////////////////////////////////////////////////////

class NBest
{
  protected:
    NBestList             reference;
    char                **vocabulary;
    LanguageModel        *languageModel;
    int                   startOfSentenceWord;
    double                org_LmScale;
    double                org_transPenalty;
    int                   adminLength;
    NBestList            *admin;
  public:
    NBest(WLRType *path, char **voc, int stS, LanguageModel *lm, double orgLm, double orgTrP);
   ~NBest();
   
  protected:
    int                  getNBestArraySize      (WLRType *w, int *nrOfPaths);
    int                  fillNBestArray         (WLRType *w, NBestType *array, int *nr, int *emptyNr, int length);
    bool                 printNextBest          ();
    bool                 bigger                 (float lmScale, float transPenalty, float target);
    
  public:
    void                 print                  (int max);
    void                 setReference           (WLRType *w);
    int                  getScore               (float lmScale, float transPenalty);
};

#endif
