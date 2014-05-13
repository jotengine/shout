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
#ifndef TRAINPHONEMODEL_H
#define TRAINPHONEMODEL_H

#include "phonemodel.h"
#include "mixgaussian.h"

#include "featurepool.h"

class Vector;
class ShoutTrainMaster;


typedef bool (*TrainCallback)();

////////////////////////////////////////////////////////////////////////////
/// \brief With help of this class, the acoustic models are trained.
///
/// Before training the HMM, training samples need to be added to the 
/// system, using the addTrainingSample() method. Once all samples are added
/// the HMM may be trained with the train() method. It is possible to either
/// receive and send the model over a socket connection (also see the 
/// Socket_Server and Socket_Client classes) with receiveModel() and 
/// sendModel() or write the model to disk with writeModel().
/// The viterbi() method is used to determine the best path for one training
/// sample. During use of the models (see the PhoneModel class), the viterbi
/// token passing, a special form of viterbi, is used.
////////////////////////////////////////////////////////////////////////////

class TrainPhoneModel : public PhoneModel
{
  // Attributes:
  protected:
    int              trainSilP;
    int              trainSilMax;
    int             *decision_Matrix;
    int              decision_numberOfModels;
    int              decision_numberOfRules;
    int              totalLength;    
    FeaturePool     *trainingPool;
    int              trainingSegment;
    int              trainingLabel;
    int              guestTrainingLabel;
    FeaturePoolInfo *channelInfoBlock;
    bool             trainWithoutBorders;
    
  public:
    TrainPhoneModel(const char *n,int contextLeft,int contextRight, bool isSil, int dim, FeaturePoolInfo *infoBlock = NULL);
    TrainPhoneModel(MixGaussian *gmm, double trans, const char *name);
    TrainPhoneModel(TrainPhoneModel *model1, TrainPhoneModel *model2, int maxGaussians = -1);
    TrainPhoneModel(TrainPhoneModel *model1, TrainPhoneModel *model2, double rate);
    TrainPhoneModel(TrainPhoneModel *orgModel, int shiftLeftRight = 0);
    TrainPhoneModel(FILE *inFile, int dim,FeaturePoolInfo *infoBlock = NULL);
   ~TrainPhoneModel();   
      
   // Methods:
  protected:
    double      baumWelch          (int trainWhat, PhoneModel *doSat = NULL);
    double      viterbi            (int trainWhat);
    double      getSilP            (int useLabel, int useSegmentation, FeaturePool *usePool);

  public:
    void        trainMMI           (FILE *fileEnum, FILE *fileDenom);
    void        doNotuseBordersForTraining(bool useBordersNot);
    int         maxNrOfGaussians   ();
    double      getTrainSilP        (int useLabel = -1, int useSegmentation = -1, FeaturePool *usePool = NULL);    
    void        adapt_setAcTrain   (int useLabel = -1, int useSegmentation = -1, FeaturePool *usePool = NULL);
//    void        writeModel         (FILE *outFile);  // moved to PhoneModel
    bool        readModel          (FILE *inFile);
    void        writeSAT           (FILE *outFile);
    void        appendSAT          (FILE *outFile);
    double      finishSAT          ();
    void        setDecisionMatrix  (int numberOfModels, int numberOfRules, int *dMatrix);
    
    void        setTrainingData    (FeaturePool *fp, int segmentationID, int labelID, int guestID = -1, int tSilP = 100, int tSilMax = -1);
    
    double      train              (int maxGaussians,bool isSil,bool neverPrune = false, Vector **trainDiscr = NULL, 
                                    Vector *trainDiscrMask = NULL, PhoneModel *doSAT = NULL, bool doFastTraining = false);
    
    double      getCoSim           (TrainPhoneModel *t1, TrainPhoneModel *t2);
    double      getKLDistance      (TrainPhoneModel *t2);
    double      getNormDistance    ();
    
    void        startCount         ();
    void        stopCount          ();
    void        count              (Vector *observation);
    int         getDominantGaussian(); 
    void        addCountedGaussians(TrainPhoneModel *source, int nmbr);
    
    void        moveModelGaussians (TrainPhoneModel *model, double factor);
    
    void        addGaussian        (Vector *v);
    
    void        normalize          ();
    
    void        setMaxGaussians    (int maxGaussians);

    #ifdef SPLIT_CLUSTERING
    double      getClusterP        (Vector *observation);
    void        fillDistanceArray  (int *distA);
    #endif
};

#endif
