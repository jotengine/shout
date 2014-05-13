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
#ifndef WHISPER_H
#define WHISPER_H

#include "standard.h"
#include "vector.h"
#include "phonemodel.h"
#include "lexicaltree.h"
#include "languagemodel.h"
#include "shoutconfig.h"


////////////////////////////////////////////////////////////////////////////
/// \brief Whisper is the top-level class of the decoder.
///
/// The constructor Whisper::Whisper() handles everything. It requires
/// file names for the acoustic models (PhoneModel), the lexical tree (LexicalTree)
/// the language model (LanguageModel) and a list of files to decode.
/// 
/// Next to LM decoding, it is also possible to do forced allignment with Whisper.
/// In this case the file-list file must contain the correct sentences after the
/// feature labels. Each sentence must start with the language model history.
/// In most cases this is (<s> <s>).
////////////////////////////////////////////////////////////////////////////

class Whisper
{
  protected:
    int                       numberOfModels;
    int                       numberOfClusters;
    int                       numberOfClusteredModels;
    int                       vectorSize;
    PhoneModel              **clusterModels;
    PhoneModel              **models;
    LexicalTree              *myTree;
    LexicalTree              *myBackgroundTree;
    LexicalTree              *myPhoneLoop;
    LanguageModel            *myLM;
    FeaturePool              *vtlnStream[VTLN_STEPS];
    SegmentationAdmin        readPool_curSeg;

    
    int                       getClusterID(FeaturePool *fPool);
  public:
    Whisper(ShoutConfig *myConfig);
   ~Whisper();

};

#endif
