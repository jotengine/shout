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
  
#ifndef SPEAKER_SEGMENTER_H
#define SPEAKER_SEGMENTER_H

#include "standard.h"
#include "mixgaussian.h"
#include "phonemodel.h"
#include "vector.h"
#include "lexicaltree.h"

class  LexicalTree;


////////////////////////////////////////////////////////////////////////////
/// \brief This class can determine which cluster a certain stream of
/// feature vectors is closest to. The class Train_Segmenter
/// can train these clusters.
////////////////////////////////////////////////////////////////////////////

class Segmenter : public LexicalTree
{
  protected:
    int               numberOfClusters;
    FeaturePool      *featurePool;
    Vector          **discrTrain;
    Vector          **discrTrainMask;
    int               discrTrainLen;

  protected:
    void createLexicalTree                (int minNumberOfFrames, int *forEachFrames1 = NULL, double *forEachFrames2 = NULL);
    float createSegments                  (int segID, SegmentationAdmin *admin); 
  public:
    // Constructor and destructor:
    Segmenter                             ();
    Segmenter                             (FILE *inFile);
   ~Segmenter                             ();

    // methods:
   float segmentFeaturePool               (int inputSegID, int inputLabelID, int outputSegID);
   void setFeaturePool                    (FeaturePool *fp);

   // INLINE FUNCTIONS FOR CODE OPTIMALIZATION:
   inline int getNumberOfModels           ()  const {return numberOfClusters;}

};

#endif // SPEAKER_SEGMENTER_H
