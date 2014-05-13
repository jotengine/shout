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
#ifndef ADAPT_SPEAKER_SEGMENTER_H
#define ADAPT_SPEAKER_SEGMENTER_H

#include "train_segmenter.h"
#include "trainphonemodel.h"

class Adapt_Segmenter : public Train_Segmenter
{
  protected:
    int                winnerM1;
    int                winnerM2;
    double             mergeWinThreshold;
    TrainPhoneModel  **compareModels;
    int                nrCompareModels;
  
    int                ubm_word_amount;
    int                fixedSpeakers;
    int                systemID;
    int                segAmount;
    double             maxDistance;
    double            *pTable1;    
    double            *pTable2;
    int                pTlen;
    char               label[100];
    int               *modelMapping;
    int                segInitID;
    TrainPhoneModel  **masterAdapt;
    TrainPhoneModel  **adaptedModel;
    TrainPhoneModel  **mergeTestMod;
    int               *helpID;
        
  public:
    Adapt_Segmenter               (FeaturePool *fp, int sID_train, int sID_decode, char *name, int sysID);
   ~Adapt_Segmenter               ();

    void   createFixedModels      (int nrSpkrs, int fixedSeg);

  protected:
    
    void   startMergeIteration    ();
    double getMergeModelScore     (int model1, int model2, int method);
    bool   proceedMerge           (int model1, int model2, int method);
    void   mergeModels            (int model1, int model2);
    
#ifdef NOT_USED_NOW    
    double getMatrixScore         (int a, int model1, int model2, int data1, int data2, int seg);
    double trainModel             (int model , int nrG);
#endif    
    void   adaptModel             (TrainPhoneModel *model);
};

#endif
