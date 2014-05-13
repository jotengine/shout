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

#ifndef LANGUAGEMODEL_SEGMENTER_H
#define LANGUAGEMODEL_SEGMENTER_H

#include <languagemodel.h>

class LanguageModel_Segmenter : public LanguageModel
{
  protected:
    int *transfer_uni;
    int *transfer_bi;
    float  *probs_bi;
    int  currentWord;
    int  totalTransfer;
  public:
    LanguageModel_Segmenter (int uniNr);
   ~LanguageModel_Segmenter ();

    void transferTo         (int word);
    void finishModel        ();
    void addPenalty         (int word, double penalty);
    virtual float getP      (int newWordID,int *wordHistoryIn, int *wordHistoryOut);
};

#endif
