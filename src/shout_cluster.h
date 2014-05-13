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
#ifndef SHOUT_CLUSTER_H
#define SHOUT_CLUSTER_H


#include "train_segmenter.h"
#include "adapt_segmenter.h"

////////////////////////////////////////////////////////////////////////////
/// \brief This class is the base class for the shout_cluster application.
/// It segments an audio file based on the speaker clustering procedure
/// from the Train_Speaker_Segmenter class.
////////////////////////////////////////////////////////////////////////////

class Shout_Cluster
{
  protected:
#ifdef ADAPT_DIARIZATION
    Adapt_Segmenter *segmenter;
#else
    Train_Segmenter *segmenter;
#endif
    
  public:
    Shout_Cluster(char *featureFile, char *label, char *inputRTTM_train, char *inputRTTM_decode, 
                  char *outputRTTM, char *fea2, double initWeight, char *fea3,
                  int fastMerge, int maxClusters, bool widen);
   ~Shout_Cluster();

};

#endif
