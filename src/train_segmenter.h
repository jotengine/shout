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
#ifndef TRAIN_SPEAKER_SEGMENTER_H
#define TRAIN_SPEAKER_SEGMENTER_H

#include "standard.h"
#include "segmenter.h"
#include "trainphonemodel.h"

////////////////////////////////////////////////////////////////////////////
/// \brief This class can 'train' new speaker clusters. 
/// One could call determining the clusters training because the phone (SIL)
/// training procedures are used. But in fact, training is done on the
/// target data (not on a training set), so it is not really training but
/// processing. See our Spring 2006 NIST Rich Transcription Speaker Diarization
/// paper for more information.
/// This is ment to be the most simple form of the HMM-based merging
/// diarization method. Overloaded classes may be created for testing
/// new algorithms.
////////////////////////////////////////////////////////////////////////////

class Train_Segmenter : public Segmenter
{
  protected:
    char              outputFileName[2000];
    char              label[200];
    int               sadID_train;
    int               sadID_decode;
    int               segID;
    int               segID_decode;
    int               segID_cheat;
    TrainPhoneModel **trainCluster;
    TrainPhoneModel **scoreCluster;
    int               nrScoreClusters;
    Vector          **distanceVect;
    bool             *mergeable;
    int               nrNonMergeModels;
    double           *clusterScore;
    int              *clusterSize;
    int               cheatSpk1;
    int               cheatSpk2;
    int               numberOfMerges;
    int              *compareClusters;
    PhoneModel       *vtlnModel;
    int               fastMerge;
    bool              prepareForASR;
    int               maxDataPoints;
    bool              detectOverlap;
    int               nrMergeIterations;
    double            bic_meanoffset[10];

    TrainPhoneModel  *trcl;

  public:
    Train_Segmenter                       (FeaturePool *fp, int sID_train, int sID_decode, const char *name, 
                                           int nrMerges, int minClusters, int maxClusters, bool widen, const char *nonMergeableModels);
    virtual ~Train_Segmenter              (); 
  
  protected:
    double trainIteration                 (int nrIterations);
    double trainClusters                  (int nrG = -1);
    bool   mergeClusters                  (int maxClusters, bool smallestWins, bool actuallyMerge);

            void   createInitialModels    (int segAmount, char *outputName = NULL);
            void   createOverlapTree      (int min);
    virtual void   startMergeIteration    ();
    virtual void   getMergeModelScore     (int model1, int model2, double *mergeScore);
    virtual bool   proceedMerge           (int model1, int model2, int method);
    virtual void   mergeModels            (int model1, int model2);
    virtual void   trainModel             (int model , int nrG, double *trainPRes);

    void   getOverlap                     (FILE *file);
    
    void   writePosteriors                (char *fileName);
    
  public:
    void   train                          (int maxClusters, char *label, char *tempStr, char *feaPosteriors, int segID_cheat = -1);
    inline int getSegmentation            ()  const {return segID_decode;}
    void   storeClusters                  (FILE *outFile, char *outputName);
    void   loadClusters                   (FILE *inFile);
};

#endif
