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

#include "standard.h"
#include "trainphonemodel.h"
#include "multimixgaussian.h"
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <stdlib.h>
#include "shout-misc.h"

using namespace WriteFileLittleBigEndian;


#define absolute(x)    (((x)<(0.0))?(-1.0*x):(x))
#define ROOFVALUE(x,a)   if(x > TRAIN_BIGGEST_P){x = TRAIN_BIGGEST_P; roofNumber++;}

int MINPASSES = -10;  // add one for the number of iterations...

void *testP;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor is used to create an empty training model. The name (input parameter) of the
/// phone is stored in the statistics structure. The number of gaussians per state is set to one.
/////////////////////////////////////////////////////////////////////////////////////////////////////

TrainPhoneModel::TrainPhoneModel(const char *n, int contextLeft, int contextRight, bool isSil, int dim, FeaturePoolInfo *infoBlock)
{
  trainWithoutBorders    = false;
  trainSilMax            = -1;
  trainSilP              = 100;
  dimensions             = dim;
  trainingPool           = NULL;
  trainingSegment        = 0;
  trainingLabel          = -1;  // meaning don't care..
  guestTrainingLabel     = -1;  // Not there...
  channelInfoBlock       = infoBlock;

  strncpy(statistics.name,n,10);
  statistics.nrOfTrainOcc  = 0;
  statistics.nrOfGaussians = 1;
  statistics.likelihood          = -9e300;
  statistics.frameMeanLikelihood =    0.0;

  decision_Matrix        = NULL;
  statistics.maxNrOfContexts = contextLeft*contextRight;

  timeStamp  = new int[statistics.maxNrOfContexts];
  stateMix_1 = new int[statistics.maxNrOfContexts];
  stateMix_2 = new int[statistics.maxNrOfContexts];
  stateMix_3 = new int[statistics.maxNrOfContexts];

  for(int i=0;i<statistics.maxNrOfContexts;i++)
  {
    timeStamp[i] = 9272312;
  }
  if(isSil)
  {
    statistics.isSil        = 1.0;
    statistics.nrOfContexts = 1;
  }
  else
  {
    statistics.isSil        = 0.0;
    statistics.nrOfContexts = 3;
  }

  mixtureSetData = new MixtureSet[statistics.nrOfContexts];
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    if(channelInfoBlock == NULL || channelInfoBlock->numberOfChannels == 1)
    {
      mixtureSetData[i].state = new MixGaussian(dim);
    }
    else
    {
      mixtureSetData[i].state = new MultiMixGaussian(channelInfoBlock, dim);
    }
    mixtureSetData[i].transitionP_toSelf = 0.0;
    mixtureSetData[i].transitionP_toNext = 0.0;
  }

  for(int i=0;i<statistics.maxNrOfContexts;i++)
  {
    stateMix_1[i] = 0;
    stateMix_2[i] = MAX_NUMBER_OF_MIXTURES_PER_PHONE/3;
    stateMix_3[i] = MAX_NUMBER_OF_MIXTURES_PER_PHONE/3*2;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor is used to create a new PhoneModel object with only one state (SIL),
/// outgoing transition probability trans and GMM gmm. (used by Train_Speaker_Segmenter)
/////////////////////////////////////////////////////////////////////////////////////////////////////

TrainPhoneModel::TrainPhoneModel(MixGaussian *gmm, double trans, const char *name)
{
  trainWithoutBorders    = false;
  trainSilMax            = -1;
  trainSilP              = 100;
  channelInfoBlock       = NULL;
  dimensions     = gmm->dim();
  statistics.nrOfContexts    = 1;
  statistics.maxNrOfContexts = 1;
  strcpy(statistics.name,name);
  trainingPool           = NULL;
  trainingSegment        = 0;
  trainingLabel          = -1;  // meaning don't care..
  guestTrainingLabel     = -1;  // Not there...
  decision_Matrix        = NULL;
  mixtureSetData         = new MixtureSet[statistics.nrOfContexts];
  timeStamp              = new int[statistics.maxNrOfContexts];
  stateMix_1             = new int[statistics.maxNrOfContexts];
  stateMix_2             = new int[statistics.maxNrOfContexts];
  stateMix_3             = new int[statistics.maxNrOfContexts];
  timeStamp[0]           = 9272311;
  stateMix_1[0]          = 0;
  stateMix_2[0]          = 0;
  stateMix_3[0]          = 0;
  mixtureSetData[0].transitionP_toSelf = log(1.0-trans);
  mixtureSetData[0].transitionP_toNext = log(trans);
  mixtureSetData[0].state              = new MixGaussian(gmm);
  statistics.nrOfGaussians = mixtureSetData[0].state->getNumberOfGaussians();
  statistics.nrOfTrainOcc = 1;
  statistics.likelihood = 1.0;
  statistics.frameMeanLikelihood = 1.0;
  statistics.isSil = 1;
  isSil          = true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor will merge two models (it doesn't train them, use train for that after setting
/// the training data setTrainingData() with a guestID!)
/// THIS IMPLEMENTATION DOES ONLY MERGE THE FIRST CONTEXT OF SIL PHONES!!!
/////////////////////////////////////////////////////////////////////////////////////////////////////

TrainPhoneModel::TrainPhoneModel(TrainPhoneModel *model1, TrainPhoneModel *model2, int maxGaussians)
{
  assert(model1 != NULL);
  assert(model2 != NULL);
  assert(model1->trainingPool != NULL);
  assert(model2->trainingPool != NULL);
  assert(model1->trainingPool == model2->trainingPool);
  assert(model1->trainingSegment == model2->trainingSegment);
  assert(model1->dimensions == model2->dimensions);

  trainWithoutBorders= model1->trainWithoutBorders;
  channelInfoBlock   = model1->channelInfoBlock;
  trainingPool       = model1->trainingPool;
  trainingSegment    = model1->trainingSegment;
  trainingLabel      = model1->trainingLabel;
  guestTrainingLabel = model2->trainingLabel;
  dimensions         = model1->dimensions;

  strcpy(statistics.name,"COMB");

  trainSilMax            = -1;
  trainSilP              = 100;
  statistics.frameMeanLikelihood = 0.0;
  statistics.likelihood      = -9e300;
  statistics.isSil           = 1.0;
  statistics.nrOfTrainOcc    = model1->statistics.nrOfTrainOcc  + model2->statistics.nrOfTrainOcc;
  statistics.nrOfGaussians   = model1->statistics.nrOfGaussians + model2->statistics.nrOfGaussians;
  statistics.maxNrOfContexts = 1;
  decision_Matrix            = NULL;
  if(model1->isSil)
  {
    statistics.nrOfContexts    = 1;
    isSil                      = true;
    assert(model2->isSil);
    assert(model1->statistics.nrOfContexts == 1);
    assert(model2->statistics.nrOfContexts == 1);
  }
  else
  {
    assert(!model2->isSil);
    assert(model1->statistics.nrOfContexts == 3);
    assert(model2->statistics.nrOfContexts == 3);
    statistics.nrOfContexts    = 3;
    isSil                      = false;
  }
  timeStamp                  = new int[statistics.maxNrOfContexts];
  stateMix_1                 = new int[statistics.nrOfContexts];
  stateMix_2                 = new int[statistics.nrOfContexts];
  stateMix_3                 = new int[statistics.nrOfContexts];
  mixtureSetData             = new MixtureSet[statistics.nrOfContexts];

  for(int i=0;i<statistics.maxNrOfContexts;i++)
  {
    timeStamp[i] = 9272312;
  }
  if(model1->isSil)
  {
    stateMix_1[0] = 0;
    mixtureSetData[0].transitionP_toSelf = 0.0;
    mixtureSetData[0].transitionP_toNext = 0.25;
    if(channelInfoBlock == NULL || channelInfoBlock->numberOfChannels == 1)
    {
      mixtureSetData[0].state = new MixGaussian(model1->mixtureSetData[0].state, model2->mixtureSetData[0].state);
    }
    else
    {
      mixtureSetData[0].state = new MultiMixGaussian(model1->mixtureSetData[0].state, model2->mixtureSetData[0].state);
    }
  }
  else
  {
    stateMix_1[0] = 0;
    stateMix_2[0] = 1;
    stateMix_3[0] = 2;
    for(int i = 0;i < 3; i++)
    {
/*      mixtureSetData[i].transitionP_toSelf = model1->mixtureSetData[i].transitionP_toSelf;
      mixtureSetData[i].transitionP_toNext = model1->mixtureSetData[i].transitionP_toNext;*/
      mixtureSetData[i].transitionP_toSelf = log( (exp(model1->mixtureSetData[i].transitionP_toSelf) +
                                                  exp(model1->mixtureSetData[i].transitionP_toSelf)   ) / 2.0 );
      mixtureSetData[i].transitionP_toNext = log( (exp(model1->mixtureSetData[i].transitionP_toNext) +
                                                  exp(model1->mixtureSetData[i].transitionP_toNext)   ) / 2.0 );
      if(channelInfoBlock == NULL || channelInfoBlock->numberOfChannels == 1)
      {
        mixtureSetData[i].state = new MixGaussian(model1->mixtureSetData[i].state, model2->mixtureSetData[i].state);
      }
      else
      {
        mixtureSetData[i].state = new MultiMixGaussian(model1->mixtureSetData[i].state, model2->mixtureSetData[i].state);
      }
    }
  }

  if((maxGaussians > 0) && (maxNrOfGaussians() > maxGaussians))
  {
    mixtureSetData[0].state->normalizeWeights(maxGaussians);
    if(!model1->isSil)
    {
      mixtureSetData[1].state->normalizeWeights(maxGaussians);
      mixtureSetData[2].state->normalizeWeights(maxGaussians);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor will merge two models (it doesn't train them, use train for that after setting
/// the training data setTrainingData() with a guestID!)
/// THIS IMPLEMENTATION DOES ONLY MERGE THE FIRST CONTEXT OF SIL PHONES!!!
/// The number of gaussians in the models must be the same! In fact, it only works for models
/// that are adapted from the same UBM!!!
/////////////////////////////////////////////////////////////////////////////////////////////////////

TrainPhoneModel::TrainPhoneModel(TrainPhoneModel *model1, TrainPhoneModel *model2, double rate)
{
  assert(model1 != NULL);
  assert(model2 != NULL);
  assert(model1->isSil);
  assert(model2->isSil);
  assert(model1->trainingPool != NULL);
  assert(model2->trainingPool != NULL);
  assert(model1->trainingPool == model2->trainingPool);
  assert(model1->trainingSegment == model2->trainingSegment);
  assert(model1->statistics.nrOfGaussians == model2->statistics.nrOfGaussians);
  assert(model1->dimensions == model2->dimensions);
  assert(model1->statistics.nrOfContexts == 1);
  assert(model2->statistics.nrOfContexts == 1);

  trainWithoutBorders= model1->trainWithoutBorders;
  channelInfoBlock   = model1->channelInfoBlock;
  trainingPool       = model1->trainingPool;
  trainingSegment    = model1->trainingSegment;
  trainingLabel      = model1->trainingLabel;
  guestTrainingLabel = model2->trainingLabel;
  dimensions         = model1->dimensions;

  strcpy(statistics.name,model1->statistics.name);

  trainSilMax            = -1;
  trainSilP              = 100;

  statistics.frameMeanLikelihood = 0.0;
  statistics.likelihood      = -9e300;
  statistics.isSil           = 1.0;
  statistics.nrOfTrainOcc    = model1->statistics.nrOfTrainOcc  + model2->statistics.nrOfTrainOcc;
  statistics.nrOfGaussians   = model1->statistics.nrOfGaussians;
  statistics.maxNrOfContexts = 1;
  statistics.nrOfContexts    = 1;
  isSil                      = true;
  decision_Matrix            = NULL;
  timeStamp                  = new int[1];
  stateMix_1                 = new int[1];
  stateMix_2                 = new int[1];
  stateMix_3                 = new int[1];
  mixtureSetData             = new MixtureSet[1];

  for(int i=0;i<statistics.maxNrOfContexts;i++)
  {
    timeStamp[i] = 9272312;
  }
  stateMix_1[0] = 0;
  mixtureSetData[0].transitionP_toSelf = 0.0;
  mixtureSetData[0].transitionP_toNext = 0.25;
  if(channelInfoBlock == NULL || channelInfoBlock->numberOfChannels == 1)
  {
    mixtureSetData[0].state  = new MixGaussian(model1->mixtureSetData[0].state, model2->mixtureSetData[0].state,rate);
  }
  else
  {
    mixtureSetData[0].state  = new MultiMixGaussian(model1->mixtureSetData[0].state, model2->mixtureSetData[0].state,rate);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The constructor only initialises some variables. Not very interesting. The parameters are
/// only added to make inherentece possible with the PhoneModel class.
/////////////////////////////////////////////////////////////////////////////////////////////////////

TrainPhoneModel::TrainPhoneModel(FILE *inFile, int dim, FeaturePoolInfo *infoBlock)
{
  trainWithoutBorders    = false;
  silRinglastPos         = 1;
  trainSilMax            = -1;
  trainSilP              = 100;
  channelInfoBlock       = infoBlock;
  trainingPool           = NULL;
  trainingSegment        = 0;
  trainingLabel          = -1;  // meaning don't care..
  guestTrainingLabel     = -1;  // not There
  dimensions             = dim;
  char notUsed;
  decision_Matrix        = NULL;

  freadEndianSafe(statistics.name,10,1,inFile);
  freadEndianSafe(&notUsed,1,sizeof(notUsed),inFile);
  freadEndianSafe(&notUsed,1,sizeof(notUsed),inFile);
  freadEndianSafe(&statistics.nrOfGaussians,1,sizeof(statistics.nrOfGaussians),inFile);
  freadEndianSafe(&statistics.nrOfContexts,1,sizeof(statistics.nrOfContexts),inFile);
  freadEndianSafe(&statistics.maxNrOfContexts,1,sizeof(statistics.maxNrOfContexts),inFile);
  freadEndianSafe(&statistics.nrOfTrainOcc,1,sizeof(statistics.nrOfTrainOcc),inFile);
  freadEndianSafe(&statistics.likelihood,1,sizeof(statistics.likelihood),inFile);
  freadEndianSafe(&statistics.frameMeanLikelihood,1,sizeof(statistics.frameMeanLikelihood),inFile);
  freadEndianSafe(&statistics.isSil,1,sizeof(statistics.isSil),inFile);

  PRINTF4("# Reading model %s (Gaussians:%3d,Contexts:%3d)\n",
                  statistics.name,statistics.nrOfGaussians,statistics.nrOfContexts);

  mixtureSetData = new MixtureSet[statistics.nrOfContexts];

  timeStamp      = new int[statistics.maxNrOfContexts];
  stateMix_1     = new int[statistics.maxNrOfContexts];
  stateMix_2     = new int[statistics.maxNrOfContexts];
  stateMix_3     = new int[statistics.maxNrOfContexts];

  for(int i=0;i<statistics.maxNrOfContexts;i++)
  {
    timeStamp[i] = 9272312;
  }
  freadEndianSafe(stateMix_1,statistics.maxNrOfContexts,sizeof(int),inFile);
  freadEndianSafe(stateMix_2,statistics.maxNrOfContexts,sizeof(int),inFile);
  freadEndianSafe(stateMix_3,statistics.maxNrOfContexts,sizeof(int),inFile);

  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    freadEndianSafe(&(mixtureSetData[i].transitionP_toSelf),1,sizeof((mixtureSetData[i].transitionP_toSelf)),inFile);
    freadEndianSafe(&(mixtureSetData[i].transitionP_toNext),1,sizeof((mixtureSetData[i].transitionP_toNext)),inFile);
    if(channelInfoBlock == NULL || channelInfoBlock->numberOfChannels == 1)
    {
      mixtureSetData[i].state = new MixGaussian(inFile,dim);
    }
    else
    {
      mixtureSetData[i].state = new MultiMixGaussian(inFile, channelInfoBlock, dim);
    }
  }
  isSil = (statistics.isSil == 1.0);
  if(isSil)
  {
    mixtureSetData[0].transitionP_toSelf = -0.1249;
    mixtureSetData[0].transitionP_toNext = -0.6021;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The constructor copies the settings of the original TrainPhoneModel.
/// This constructor is used for creating clustered-based acoustic models.
/////////////////////////////////////////////////////////////////////////////////////////////////////

TrainPhoneModel::TrainPhoneModel(TrainPhoneModel *org, int shiftLeftRight)
{
  assert(org != NULL);

  trainWithoutBorders    = org->trainWithoutBorders;
  trainSilMax            = -1;
  trainSilP              = 100;
  channelInfoBlock       = org->channelInfoBlock;
  trainingPool           = NULL;
  trainingSegment        = 0;
  trainingLabel          = -1;  // meaning don't care..
  guestTrainingLabel     = -1;
  dimensions             = org->dimensions;
  decision_Matrix        = NULL;

  isSil                  = org->isSil;
  totalLength            = org->totalLength;

  strncpy(statistics.name,(org->statistics).name,10);
  statistics.nrOfGaussians       = org->statistics.nrOfGaussians;
  statistics.nrOfContexts        = org->statistics.nrOfContexts;
  statistics.maxNrOfContexts     = org->statistics.maxNrOfContexts;
  statistics.nrOfTrainOcc        = org->statistics.nrOfTrainOcc;
  statistics.likelihood          = org->statistics.likelihood;
  statistics.frameMeanLikelihood = org->statistics.frameMeanLikelihood;
  statistics.isSil               = org->statistics.isSil;

  mixtureSetData = new MixtureSet[statistics.nrOfContexts];

  timeStamp      = new int[statistics.maxNrOfContexts];
  stateMix_1     = new int[statistics.maxNrOfContexts];
  stateMix_2     = new int[statistics.maxNrOfContexts];
  stateMix_3     = new int[statistics.maxNrOfContexts];

  for(int i=0;i<statistics.maxNrOfContexts;i++)
  {
    timeStamp[i] = 9272312;
  }

  memcpy(stateMix_1, org->stateMix_1, statistics.maxNrOfContexts*sizeof(int));
  memcpy(stateMix_2, org->stateMix_2, statistics.maxNrOfContexts*sizeof(int));
  memcpy(stateMix_3, org->stateMix_3, statistics.maxNrOfContexts*sizeof(int));

  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    memcpy(&(mixtureSetData[i].transitionP_toSelf), &(org->mixtureSetData[i].transitionP_toSelf),
           sizeof((mixtureSetData[i].transitionP_toSelf)));
    memcpy(&(mixtureSetData[i].transitionP_toNext), &(org->mixtureSetData[i].transitionP_toNext),
           sizeof((mixtureSetData[i].transitionP_toNext)));
    if(channelInfoBlock == NULL || channelInfoBlock->numberOfChannels == 1)
    {
      mixtureSetData[i].state = new MixGaussian(org->mixtureSetData[i].state);
    }
    else
    {
      mixtureSetData[i].state = new MultiMixGaussian(org->mixtureSetData[i].state);
    }
  }
  if(shiftLeftRight != 0)
  {
    for(int i=0;i<statistics.nrOfContexts;i++)
    {
      mixtureSetData[i].state->shiftBestGaussian(shiftLeftRight);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor is responsible for deleting the entire training pool.
/////////////////////////////////////////////////////////////////////////////////////////////////////

TrainPhoneModel::~TrainPhoneModel()
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Load training parameters of a model from disc.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainPhoneModel::appendSAT(FILE *outFile)
{
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    mixtureSetData[i].state->appendSAT(outFile);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Finish the training for SAT
/////////////////////////////////////////////////////////////////////////////////////////////////////

double TrainPhoneModel::finishSAT()
{
  double res = 0.0;
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    res += mixtureSetData[i].state->trainFinish(true);
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// WriteSAT writes the training parameters of a model to disc. These data will later be used
/// by the SAT training application.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainPhoneModel::writeSAT(FILE *outFile)
{
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    mixtureSetData[i].state->storeSAT(outFile);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// readModel reads a model from disc. Especially needed for MultiMixGaussian PDFs
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool TrainPhoneModel::readModel(FILE *inFile)
{
  char notUsed;
  if(freadEndianSafe(statistics.name,1,10,inFile) > 0)
  {
    freadEndianSafe(&notUsed,1,sizeof(notUsed),inFile);
    freadEndianSafe(&notUsed,1,sizeof(notUsed),inFile);
    freadEndianSafe(&statistics.nrOfGaussians,1,sizeof(statistics.nrOfGaussians),inFile);
    freadEndianSafe(&statistics.nrOfContexts,1,sizeof(statistics.nrOfContexts),inFile);
    int max = statistics.maxNrOfContexts;
    freadEndianSafe(&statistics.maxNrOfContexts,1,sizeof(statistics.maxNrOfContexts),inFile);
    assert(max == statistics.maxNrOfContexts);
    freadEndianSafe(&statistics.nrOfTrainOcc,1,sizeof(statistics.nrOfTrainOcc),inFile);
    freadEndianSafe(&statistics.likelihood,1,sizeof(statistics.likelihood),inFile);
    freadEndianSafe(&statistics.frameMeanLikelihood,1,sizeof(statistics.frameMeanLikelihood),inFile);
    freadEndianSafe(&statistics.isSil,1,sizeof(statistics.isSil),inFile);

    freadEndianSafe(stateMix_1,statistics.maxNrOfContexts,sizeof(int),inFile);
    freadEndianSafe(stateMix_2,statistics.maxNrOfContexts,sizeof(int),inFile);
    freadEndianSafe(stateMix_3,statistics.maxNrOfContexts,sizeof(int),inFile);
    for(int i=0;i<statistics.nrOfContexts;i++)
    {
      freadEndianSafe(&(mixtureSetData[i].transitionP_toSelf),1,sizeof((mixtureSetData[i].transitionP_toSelf)),inFile);
      freadEndianSafe(&(mixtureSetData[i].transitionP_toNext),1,sizeof((mixtureSetData[i].transitionP_toNext)),inFile);
      delete mixtureSetData[i].state;
      if(channelInfoBlock == NULL || channelInfoBlock->numberOfChannels == 1)
      {
        mixtureSetData[i].state = new MixGaussian(inFile,dimensions);
      }
      else
      {
        mixtureSetData[i].state = new MultiMixGaussian(inFile, channelInfoBlock, dimensions);
      }
    }
    isSil = statistics.isSil;
    return true;
  }
  return false;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will set a new decision-rule matrix. Each row in this matrix contains all phones that
/// are part of a particular tree-based clustering rule.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainPhoneModel::setDecisionMatrix(int numberOfModels, int numberOfRules, int *dMatrix)
{
  decision_Matrix = dMatrix;
  decision_numberOfModels = numberOfModels;
  decision_numberOfRules  = numberOfRules;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The FeaturePool contains all training data!
/// Deleting the training pool after it has no more use, is the responsability of the user, because
/// one pool maybe used by multiple models!
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainPhoneModel::setTrainingData(FeaturePool *fp, int trainID, int labelID, int guestID, int tSilP, int tSilMax)
{
  trainSilP          = tSilP;
  trainSilMax        = tSilMax;
  trainingPool       = fp;
  trainingSegment    = trainID;
  trainingLabel      = labelID;
  guestTrainingLabel = guestID;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method retrieves the number of gausians of each state and returns the maximum value.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int TrainPhoneModel::maxNrOfGaussians()
{
  int res = -1;
  // This is not a sil phone. Those get trained elsewhere...
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    int s = mixtureSetData[i].state->getNumberOfGaussians();
    if(res < s)
    {
      res = s;
    }
  }
  return res;
}

#ifdef SPLIT_CLUSTERING
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Fills the distance array: an array of size numberOfGaussians (mixtureSetData[0]). The two
/// most distant gaussians are in the first and last spot. The others are sorted in between.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainPhoneModel::fillDistanceArray(int *distA)
{
  mixtureSetData[0].state->fillDistanceArray(distA);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the probability of the mixtureSet 0 on the input vector.
/////////////////////////////////////////////////////////////////////////////////////////////////////

double TrainPhoneModel::getClusterP(Vector *observation)
{
  return mixtureSetData[0].state->getP(observation);
}

#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Set the maximum number of gaussians...
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainPhoneModel::setMaxGaussians(int maxGaussians)
{
  if(isSil && maxNrOfGaussians() > maxGaussians)
  {
    mixtureSetData[0].state->normalizeWeights(maxGaussians);
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method has the same function as the method viterbi(), to train the acoustic models.
/// The algorithm used in this case is Baum-Welch. the method viterbi() uses Viterbi.
///
/// This method uses the Baum-Welch algorithm on all training samples from the pool. The total likelihood
/// (the product of all likelihoods) is returned. The heigher this value, the better the training samples
/// match the HMM. When the train parameter is set to true, the number of transitions into each state
/// are stored and the MixGaussian::train() methods are called. The transition statistics are used
/// to determine new state transition values and MixGaussian::trainFinish() is called to finish
/// the state likelihood training.
///
/// The viterbi algorithm used in this method is as follows:
/// For every observation (in one training sample),
/// the highest probable state transition into each state is calculated.
/// Basically, this means comparing between two transitions (for example, state 2 may get input
/// from state 1 or from itself). The transition choosen is stored in the variable transPath.
/// After all observations are handled, the best path is calculated by starting at the final state
/// at time observationLength and 'walking' back in the transPath matrix.
/// The transition weights are calculated directly from the 'best path' transitions by counting the
/// number of transitions into a state and the number of transitions into that state using a specific
/// route. When all samples are handled, those two numbers are devided, and the transition chance
/// is the result. The state likelihoods are trained with those observations from the best path that are
/// mapped onto a state->
/// For the procedure to calculate the state likelihoods out of these observation sequences, see
/// the TrainMixGaussian::train() method.
///
/// Literature: Frederick Jelinek, "Statistical Methods for Speech Recognition": page 22 for the viterbi
/// training procedure and page 30 for transition weight calculation.
/////////////////////////////////////////////////////////////////////////////////////////////////////

double TrainPhoneModel::viterbi(int trainWhat)
{
  assert(trainingPool != NULL);

/* BE WARE: ALL IS IN LOG!!!!!! */

  int state[3];

  int        transitionCount_Self[statistics.nrOfContexts];
  int        transitionCount_Next[statistics.nrOfContexts];
  int        t,i,j;
  Vector    *trainVector = NULL;
  double     totalLogP    = 0.0;
  double     currentP[3];
  double     p[2];

  if(trainWhat == 4) // fast training
  {
    trainingPool->resetSegmentation(TRAIN_CLUSTER_ID1);
    trainingPool->resetSegmentation(TRAIN_CLUSTER_ID2);
    trainingPool->resetSegmentation(TRAIN_CLUSTER_ID3);
  }

  for(i=0;i<statistics.nrOfContexts;i++)
  {
    transitionCount_Self[i] = 0;
    transitionCount_Next[i] = 0;
  }

  statistics.frameMeanLikelihood = 0.0;
  SegmentationAdmin readPool_curSeg;
  trainVector = trainingPool->getFirstVectorFirstSegment(&readPool_curSeg, trainingSegment,trainingLabel,trainingLabel < 0);
  statistics.nrOfTrainOcc = 0;
  while(trainVector != NULL)
  {
    statistics.nrOfTrainOcc++;
    int observationLength = trainingPool->getCurSegmentLen(&readPool_curSeg);
    int contextKey        = trainingPool->getSegmentID(TRAIN_CONTEXTKEY,&readPool_curSeg);

    // This is not a sil phone. Those get trained elsewhere...
    // So let's determine the state mix for this observation!
    state[0] = stateMix_1[contextKey];
    state[1] = stateMix_2[contextKey];
    state[2] = stateMix_3[contextKey];

    int transPath[observationLength][3];

    for(t=observationLength-1;t>=0;t--)
    {
      transPath[t][0] = -4;
      transPath[t][1] = -4;
      transPath[t][2] = -4;
    }

    // Observation 0:
    currentP[0] = mixtureSetData[state[0]].state->getLogP(trainVector);
    currentP[1] = SMALLEST_NEGATIVE;
    currentP[2] = SMALLEST_NEGATIVE;

    // Viterbi:

    // Observation 1 to observationLength:
    for(t=1;t<observationLength;t++)
    {
      trainVector = trainingPool->getCurSegmentVector(&readPool_curSeg,t);
      // STATE 0:
      transPath[t][0] = 0;
      for(i=2;i>0;i--)
      {
        j = i-1;
        p[0] = mixtureSetData[state[j]].transitionP_toNext+currentP[j];
        p[1] = mixtureSetData[state[i]].transitionP_toSelf+currentP[i];
        if(p[0] > p[1])
        {
          currentP[i]     = p[0]+mixtureSetData[state[i]].state->getLogP(trainVector);
          transPath[t][i] = j;
        }
        else
        {
          currentP[i]     = p[1]+mixtureSetData[state[i]].state->getLogP(trainVector);
          transPath[t][i] = i;
        }
      }
      currentP[0] = mixtureSetData[state[0]].transitionP_toSelf +currentP[0]+ mixtureSetData[state[0]].state->getLogP(trainVector);
    }

    if(trainWhat == 1)   // Transitions:
    {
      i = 2;
      for(t=observationLength-1;t>0;t--)
      {
        j = transPath[t][i];
        if(i==j)
        {
          transitionCount_Self[state[j]]++;
        }
        else
        {
          transitionCount_Next[state[j]]++;
        }
        i = j;
      }
      transitionCount_Next[state[2]]++;                // Special: outgoing transition
    }
    else if(trainWhat == 2 || trainWhat == 3) // Gaussians:
    {
      i = 2;
      for(t=observationLength-1;t>0;t--)
      {
        trainVector = trainingPool->getCurSegmentVector(&readPool_curSeg,t);
        mixtureSetData[state[i]].state->train(false,1.0,trainVector);
        j = transPath[t][i];
        i = j;
      }
      trainVector = trainingPool->getCurSegmentVector(&readPool_curSeg,0);
      mixtureSetData[state[0]].state->train(false,1.0,trainVector);
    }
    else if(trainWhat == 4) // fast training preparation
    {
      i = 2;
      int contextLeft  = trainingPool->getSegmentID(TRAIN_CONTEXTLEFT,&readPool_curSeg);
      int contextRight = trainingPool->getSegmentID(TRAIN_CONTEXTRIGHT,&readPool_curSeg);
      int timeS[4];
      timeS[0] = trainingPool->getCurSegmentStart(&readPool_curSeg);
      timeS[3] = timeS[0] + observationLength;
      timeS[1] = timeS[3];
      timeS[2] = timeS[3];
      for(t=observationLength-1;t>0;t--)
      {
        if(i>0)
        {
          timeS[i] = t + timeS[0];
        }
        j = transPath[t][i];
        i = j;
      }
      trainingPool->addSegment(TRAIN_CLUSTER_ID1,timeS[0], timeS[1]-1,0,
                                    contextKey,contextLeft,contextRight);
      trainingPool->addSegment(TRAIN_CLUSTER_ID2,timeS[1], timeS[2]-1,0,
                                    contextKey,contextLeft,contextRight);
      trainingPool->addSegment(TRAIN_CLUSTER_ID3,timeS[2], timeS[3]-1,0,
                                    contextKey,contextLeft,contextRight);
    }

    totalLogP += (currentP[2]+mixtureSetData[state[2]].transitionP_toNext);

    trainVector = trainingPool->getFirstVectorNextSegment(&readPool_curSeg,trainingLabel < 0);
  }


  if(trainWhat == 1)  // Transitions
  {
    for(i=0;i<statistics.nrOfContexts;i++)
    {
      double totCount  = (double)(transitionCount_Self[i] + transitionCount_Next[i]);
      double transSelf = ((double)transitionCount_Self[i]) / totCount;
      if(transSelf < TRANSITION_MINIMUM)
      {
        PRINTF4("Problems in transition probabilities: %d:  %d - %d\n",i, transitionCount_Self[i], transitionCount_Next[i]);
        mixtureSetData[i].transitionP_toSelf = log(TRANSITION_MINIMUM);
        mixtureSetData[i].transitionP_toNext = log(1.0 - TRANSITION_MINIMUM);
      }
      else
      {
        mixtureSetData[i].transitionP_toSelf = log(transSelf);
        mixtureSetData[i].transitionP_toNext = log(((double)transitionCount_Next[i]) / totCount);
      }
    }
  }
  else if(trainWhat == 2)  // Gaussians:
  {
    for(i=0;i<statistics.nrOfContexts;i++)
    {
      mixtureSetData[i].state->trainFinish();
    }
  }

  return totalLogP;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method has the same function as the method baumWelch(), to train the acoustic models.
/// The algorithm used in this case is Viterbi. the method baumWelch() uses Baum-Welch.
///
/// This method uses the viterbi algorithm on all training samples from the pool. The total likelihood
/// (the product of all likelihoods) is returned. The heigher this value, the better the training samples
/// match the HMM. When the train parameter is set to true, the number of transitions into each state
/// are stored and the MixGaussian::train() methods are called. The transition statistics are used
/// to determine new state transition values and MixGaussian::trainFinish() is called to finish
/// the state likelihood training.
///
/// The Baum-Welch algorithm used in this method is as follows:
/// For every observation (in one training sample),
/// all state transition into each state are calculated. This is done in a Forward- and Backward pass.
/// The probability of each path is calculated. The transition weights are calculated according
/// the these probabilities by counting the
/// number of transitions into a state and the number of transitions into that state using a specific
/// route (multiplied with the path probability).
/// When all samples are handled, those two numbers are devided, and the transition chance
/// is the result. The state likelihoods are trained with the observations from a specific path that are
/// mapped onto a state with the weight of the path probability.
/// For the procedure to calculate the state likelihoods out of these observation sequences, see
/// the TrainMixGaussian::train() method.
///
/// Literature: the HTK-BOOK, page 130-132.
 /////////////////////////////////////////////////////////////////////////////////////////////////////

double TrainPhoneModel::baumWelch(int trainWhat, PhoneModel *doSat)
{
  assert(trainingPool != NULL);

  int        state[3];
  double     transitionP_Self[3];
  double     transitionP_Next[3];
  double     transitionCount_Self_Teller[statistics.nrOfContexts];
  double     transitionCount_Next_Teller[statistics.nrOfContexts];
  double     transitionCount_Noemer[statistics.nrOfContexts];
  Vector    *trainVector;
  double     totalLogP = 0.0;
  int        outLayer[statistics.nrOfContexts];
  int        nrObserv[statistics.nrOfContexts];
  int        roofNumber = 0;

  statistics.frameMeanLikelihood = 0.0;
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    outLayer[i]                    = 0;
    nrObserv[i]                    = 0;
    transitionCount_Self_Teller[i] = 0.0;
    transitionCount_Next_Teller[i] = 0.0;
    transitionCount_Noemer[i]      = 0.0;
  }
  SegmentationAdmin readPool_curSeg;
  trainVector = trainingPool->getFirstVectorFirstSegment(&readPool_curSeg,trainingSegment,trainingLabel,trainingLabel < 0);
  while(trainVector != NULL)
  {
    int observationLength = trainingPool->getCurSegmentLen(&readPool_curSeg);
    if(observationLength < 3)
    {
      PRINTF2("# Skipping observation with length %d..\n",observationLength);
    }
    else
    {
      int contextKey        = trainingPool->getSegmentID(TRAIN_CONTEXTKEY,&readPool_curSeg);

      // This is not a sil phone. Those get trained elsewhere...
      // So let's determine the state mix for this observation!
      state[0]            = stateMix_1[contextKey];
      state[1]            = stateMix_2[contextKey];
      state[2]            = stateMix_3[contextKey];
      transitionP_Self[0] = exp(mixtureSetData[state[0]].transitionP_toSelf);
      transitionP_Self[1] = exp(mixtureSetData[state[1]].transitionP_toSelf);
      transitionP_Self[2] = exp(mixtureSetData[state[2]].transitionP_toSelf);
      transitionP_Next[0] = exp(mixtureSetData[state[0]].transitionP_toNext);
      transitionP_Next[1] = exp(mixtureSetData[state[1]].transitionP_toNext);
      transitionP_Next[2] = exp(mixtureSetData[state[2]].transitionP_toNext);


      // Baum-Welch:
      double alfa[observationLength][3];
      double Q[observationLength];
      double logQ;
      double beta[observationLength][3];

      // Forward pass:
      alfa[0][0] = 1.0;
      Q[0]       = mixtureSetData[state[0]].state->getP(trainVector);
      logQ       = log(Q[0]);
      alfa[0][1] = 0.0;
      alfa[0][2] = 0.0;

      for(int t=1;t<observationLength;t++)
      {
        trainVector = trainingPool->getCurSegmentVector(&readPool_curSeg,t);
        double b[3];
        b[0] = (mixtureSetData[state[0]].state->getP(trainVector));
        b[1] = (mixtureSetData[state[1]].state->getP(trainVector));
        b[2] = (mixtureSetData[state[2]].state->getP(trainVector));
        alfa[t][0] = alfa[t-1][0]*transitionP_Self[0]*b[0];
        alfa[t][1] = alfa[t-1][0]*transitionP_Next[0]*b[1]+alfa[t-1][1]*transitionP_Self[1]*b[1];
        alfa[t][2] = alfa[t-1][1]*transitionP_Next[1]*b[2]+alfa[t-1][2]*transitionP_Self[2]*b[2];
        Q[t]       = alfa[t][0]+alfa[t][1]+alfa[t][2];
        logQ      += log(Q[t]);
        alfa[t][0] = alfa[t][0] / Q[t];
        alfa[t][1] = alfa[t][1] / Q[t];
        alfa[t][2] = alfa[t][2] / Q[t];
        ROOFVALUE(alfa[t][0],"alfa[t][0]");
        ROOFVALUE(alfa[t][1],"alfa[t][1]");
        ROOFVALUE(alfa[t][2],"alfa[t][2]");
      }

      // Backward pass:
      beta[observationLength-1][0] = 0.0;
      beta[observationLength-1][1] = 0.0;
      beta[observationLength-1][2] = transitionP_Next[2] / Q[observationLength-1];
      for(int t=observationLength-2;t>=0;t--)
      {
        trainVector = trainingPool->getCurSegmentVector(&readPool_curSeg,t+1);
        double b[3];
        b[0] = (mixtureSetData[state[0]].state->getP(trainVector));
        b[1] = (mixtureSetData[state[1]].state->getP(trainVector));
        b[2] = (mixtureSetData[state[2]].state->getP(trainVector));
        beta[t][0] = (beta[t+1][0]*transitionP_Self[0]*b[0] + beta[t+1][1]*transitionP_Next[0]*b[1]) / Q[t];
        beta[t][1] = (beta[t+1][1]*transitionP_Self[1]*b[1] + beta[t+1][2]*transitionP_Next[1]*b[2]) / Q[t];
        beta[t][2] = (beta[t+1][2]*transitionP_Self[2]*b[2])                                         / Q[t];
        ROOFVALUE(beta[t][0],"beta[t][0]");
        ROOFVALUE(beta[t][1],"beta[t][1]");
        ROOFVALUE(beta[t][2],"beta[t][2]");
      }
      totalLogP += log(alfa[observationLength-1][2]*transitionP_Next[2])+logQ;
      ROOFVALUE(totalLogP,"totalLogP");
      nrObserv[state[0]]++;
      nrObserv[state[1]]++;
      nrObserv[state[2]]++;

      if(trainWhat == 1)   // Transitions:
      {
          double teller[6] = {0.0,0.0,0.0,0.0,0.0,1.0};
          for(int t=0;t<observationLength-1;t++)
          {
            trainVector = trainingPool->getCurSegmentVector(&readPool_curSeg,t+1);
            double b[3];
            b[0] = mixtureSetData[state[0]].state->getP(trainVector);
            b[1] = mixtureSetData[state[1]].state->getP(trainVector);
            b[2] = mixtureSetData[state[2]].state->getP(trainVector);
            // Transition from 0 to 0:
            teller[0] += (alfa[t][0]*beta[t+1][0]*b[0]*transitionP_Self[0]);
            // Transition from 1 to 1:
            teller[1] += (alfa[t][1]*beta[t+1][1]*b[1]*transitionP_Self[1]);
            // Transition from 2 to 2:
            teller[2] += (alfa[t][2]*beta[t+1][2]*b[2]*transitionP_Self[2]);
            // Transition from 0 to 1:
            teller[3] += (alfa[t][0]*beta[t+1][1]*b[1]*transitionP_Next[0]);
            // Transition from 1 to 2:
            teller[4] += (alfa[t][1]*beta[t+1][2]*b[2]*transitionP_Next[1]);
          }
          teller[5]  = alfa[observationLength-1][2]*transitionP_Next[2];
          ROOFVALUE(teller[0],"teller[0]");
          ROOFVALUE(teller[1],"teller[0]");
          ROOFVALUE(teller[2],"teller[0]");
          ROOFVALUE(teller[3],"teller[0]");
          ROOFVALUE(teller[4],"teller[0]");
          ROOFVALUE(teller[5],"teller[0]");
          transitionCount_Self_Teller[state[0]] += teller[0];
          transitionCount_Self_Teller[state[1]] += teller[1];
          transitionCount_Self_Teller[state[2]] += teller[2];

          transitionCount_Next_Teller[state[0]] += teller[3];
          transitionCount_Next_Teller[state[1]] += teller[4];
          transitionCount_Next_Teller[state[2]] += teller[5];

          transitionCount_Noemer[state[0]]      += teller[0]+teller[3];
          transitionCount_Noemer[state[1]]      += teller[1]+teller[4];
          transitionCount_Noemer[state[2]]      += teller[2]+teller[5];
      }
      else if(trainWhat == 2 || trainWhat == 3) // Gaussians or mean-adaptation:
      {
          // First observation:
          trainVector = trainingPool->getCurSegmentVector(&readPool_curSeg,0);
          mixtureSetData[state[0]].state->train(false,1.0,trainVector);
          // Following observations:
          for(int t=1;t<observationLength;t++)
          {
            trainVector = trainingPool->getCurSegmentVector(&readPool_curSeg,t);
            double a[3];
            double b[3];
            b[0] = mixtureSetData[state[0]].state->getP(trainVector);
            b[1] = mixtureSetData[state[1]].state->getP(trainVector);
            b[2] = mixtureSetData[state[2]].state->getP(trainVector);
            a[0] = (alfa[t-1][0]*beta[t][0]*b[0]*transitionP_Self[0]);
            a[1] = (alfa[t-1][1]*beta[t][1]*b[1]*transitionP_Self[1])+(alfa[t-1][0]*beta[t][1]*b[1]*transitionP_Next[0]);
            a[2] = (alfa[t-1][2]*beta[t][2]*b[2]*transitionP_Self[2])+(alfa[t-1][1]*beta[t][2]*b[2]*transitionP_Next[1]);
            double sum = a[0]+a[1]+a[2];
            a[0] /= sum;
            a[1] /= sum;
            a[2] /= sum;
            ROOFVALUE(a[0],"a[0]");
            ROOFVALUE(a[1],"a[1]");
            ROOFVALUE(a[2],"a[2]");
            // Probability of being in state 0 at time t:
            if(a[0] > TRAIN_WEIGHT_MINIMUM)
            {
              if(doSat == NULL)
              {
                mixtureSetData[state[0]].state->train(false,a[0],trainVector);
              }
              else
              {
                mixtureSetData[state[0]].state->train(false,a[0],trainVector,doSat->mixtureSetData[state[0]].state);
              }
            }
            // Probability of being in state 1 at time t:
            if(a[1] > TRAIN_WEIGHT_MINIMUM)
            {
              if(doSat == NULL)
              {
                mixtureSetData[state[1]].state->train(false,a[1],trainVector);
              }
              else
              {
                mixtureSetData[state[1]].state->train(false,a[1],trainVector,doSat->mixtureSetData[state[1]].state);
              }
            }
            // Probability of being in state 2 at time t:
            if(a[2] > TRAIN_WEIGHT_MINIMUM)
            {
              if(doSat == NULL)
              {
                mixtureSetData[state[2]].state->train(false,a[2],trainVector);
              }
              else
              {
                mixtureSetData[state[2]].state->train(false,a[2],trainVector,doSat->mixtureSetData[state[2]].state);
              }
            }
        }
      }
    }
    trainVector = trainingPool->getFirstVectorNextSegment(&readPool_curSeg,trainingLabel < 0);
  }
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    if(outLayer[i] > 100)
    {
      PRINTF4("# WARNING: Number of outlayers[%d] is %d from %d!\n",i,outLayer[i],nrObserv[i]);
    }
    if(roofNumber > 100)
    {
      PRINTF2("# WARNING: Number of calculations maximized to the ROOF value is %d!\n",roofNumber);
    }
  }

  if(trainWhat == 1)  // Transitions
  {
    for(int i=0;i<statistics.nrOfContexts;i++)
    {
       double curTrans = transitionCount_Self_Teller[i] / transitionCount_Noemer[i];
       if(curTrans < TRANSITION_MINIMUM)
       {
         PRINTF4("# WARNING: Problems in transition probabilities! %d: %g - %g\n",
                 i,transitionCount_Self_Teller[i],transitionCount_Noemer[i]);
         mixtureSetData[i].transitionP_toSelf = log(TRANSITION_MINIMUM);
         mixtureSetData[i].transitionP_toNext = log(1.0 - TRANSITION_MINIMUM);
       }
       else
       {
         mixtureSetData[i].transitionP_toSelf = log(curTrans);
         mixtureSetData[i].transitionP_toNext = log(transitionCount_Next_Teller[i] / transitionCount_Noemer[i]);
       }
     }
  }
  else if(trainWhat == 2)  // Gaussians:
  {
    for(int i=0;i<statistics.nrOfContexts;i++)
    {
      mixtureSetData[i].state->trainFinish();
    }
  }
  else if(trainWhat == 3)     // Adaptation (or SAT)
  {

  }

  return totalLogP;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Used by the train() method to calculate the current score if the phone is a SIL phone.
/// Effectively the score is the MixGaussian probability (MixGaussian::getP()) of all training
/// samples.
/// It is possible to use an alternative training pool to calculate the score.
/////////////////////////////////////////////////////////////////////////////////////////////////////

double TrainPhoneModel::getTrainSilP(int useLabel, int useSegmentation, FeaturePool *usePool)
{
  assert(trainingPool != NULL || usePool != NULL);
  assert(usePool == NULL || useSegmentation >= 0);

  statistics.frameMeanLikelihood = 0.0;

  if(usePool == NULL)
  {
    usePool = trainingPool;
  }
  if(useLabel < 0)
  {
    useLabel      = trainingLabel;
  }
  if(useSegmentation < 0)
  {
    useSegmentation = trainingSegment;
  }

  return getSilP(useLabel,useSegmentation,usePool);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

double TrainPhoneModel::getSilP(int useLabel, int useSegmentation, FeaturePool *usePool)
{
  double res = 0.0;
  bool isLast   = false;
  bool inFilter = true;
  SegmentationAdmin curSeg;
  Vector *trainVector = usePool->getFirstVectorFirstSegment(&curSeg,useSegmentation,useLabel,useLabel < 0, &isLast, &inFilter);
  while(trainVector != NULL)
  {
    int stateNr = stateMix_1[usePool->getSegmentID(TRAIN_CONTEXTKEY,&curSeg)];
    testP = (void*)mixtureSetData[stateNr].state;
    if(inFilter)
    {
      res += mixtureSetData[stateNr].state->getLogP(trainVector);
    }
    while(!isLast)
    {
      trainVector = usePool->getNextVector(&curSeg,&isLast, &inFilter);
      if(inFilter)
      {
        res += mixtureSetData[stateNr].state->getLogP(trainVector);
      }
    }
    trainVector = usePool->getFirstVectorNextSegment(&curSeg,useLabel < 0, &isLast, &inFilter);
  }
  return res;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will read training accumulators from file and perform MMI training
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainPhoneModel::trainMMI(FILE *fileEnum, FILE *fileDenom)
{
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    mixtureSetData[i].state->trainMMI(fileEnum,fileDenom);
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method is responsible for the training iteration. First it checks if the HMM is a SIL HMM or a regular one.
/// SIL HMMs only have one state and only this one state is trained. Apart from this difference, both
/// types of HMMs are trained as follows:
///
/// First it is checked if an existing HMM is already loaded. If not, an empty HMM is created.
/// Then, in a loop, the system calculates the total likelihood on all training samples by
/// calling viterbi() with its training parameter set to true.
/// When the result is more than MINIMUM_TRAIN_IMPROVEMENT better than the last result, another viterbi iteration
/// is performed. Otherwise, it is checked if it is allowed to add more gaussians to each state-> This
/// is allowed as long as during the previous iterations, the number of gaussians has not been reduced
/// (it is possible that a gaussian is pruned, when not enough training samples represent it) and the
/// maximum number of gaussians has not yet been reached. The MixGaussian::splitBestGaussian() method
/// is used to split the gaussians. After splitting, the entire procedure is repeated. If splitting
/// was not permitted, the training run has finished.
///
/// When the training run is finished, the training sample pool is deleted. This means that it is
/// not possible to train the system again, without providing new samples (but who would want to?)
/////////////////////////////////////////////////////////////////////////////////////////////////////

double TrainPhoneModel::train(int maxGaussians, bool isS, bool neverPrune, Vector **trainDiscr, Vector *trainDiscrMask, PhoneModel *doSat, bool doFastTraining)
{
  SegmentationAdmin readPool_curSeg;


  bool       proceedDoubling = true;
  double     likelihood_old    = SMALLEST_NEGATIVE;
  double     likelihood_new    = SMALLEST_NEGATIVE;
  int        lastNrOfGaussians = -1;

  int        nrOfG = statistics.nrOfGaussians;

  isSil = isS;
  if(isSil)
  {
    statistics.isSil = 1.0;
  }
  else
  {
    statistics.isSil = 0.0;
  }

  if(isSil)
  {
    PRINTF("# Going to train a SILENCE model.\n");
    /* This is a silence model. Lets keep it easy and only use the first state: */
    proceedDoubling = true;
    while(proceedDoubling)
    {
      if(!doFastTraining)
      {
        statistics.likelihood = getTrainSilP();
        if(guestTrainingLabel >= 0)
        {
          statistics.likelihood += getTrainSilP(guestTrainingLabel);
        }
      }
      else
      {
        statistics.likelihood = 0.0;
      }
      likelihood_old = -9.0e100;
      int    pass    = 0;
      int    pass2   = 0;
      int    count   = 0;
      while((pass < MINIMUM_TRAIN_ITERATIONS) || (((pass < MAXIMUM_TRAIN_ITERATIONS) && !(doFastTraining && (pass >= 5))) &&
                                                   ((statistics.likelihood - likelihood_old) / absolute(likelihood_old) > MINIMUM_TRAIN_IMPROVEMENT)))
      {
        pass++;
        pass2++;
        if(MINPASSES > 0)
        {
          if(pass2 >= MINPASSES)
          {
            pass = MAXIMUM_TRAIN_ITERATIONS + 10;   // No more iterations...
          }
          else
          {
           pass = 0;  // force iterations...
          }
        }

        if(!doFastTraining)
        {
          likelihood_old = statistics.likelihood;
        }
        int nrOfSamplesPerCluster[statistics.nrOfContexts];
        for(int i=0;i<statistics.nrOfContexts;i++)
        {
          nrOfSamplesPerCluster[i] = 0;
          mixtureSetData[i].transitionP_toNext = 0.0;
        }


        int useLabel = trainingLabel;
        bool proceed = true;
        count    = 0;
        while(proceed)  // Extra loop for training merged models (with a valid guestTrainingLabel)
        {
          /* Train */
          bool isLast   = false;
          bool inFilter = true;

          bool vectorNotNull = true;
          Vector *trainVector = trainingPool->getFirstVectorFirstSegment(&readPool_curSeg,trainingSegment,useLabel,useLabel < 0, &isLast, &inFilter);
          vectorNotNull = (trainVector != NULL);

          statistics.nrOfTrainOcc = 0;
          while(vectorNotNull)
          {
            if((trainSilMax < 0) || ((statistics.nrOfTrainOcc < trainSilMax) && (int(100.0*rand()/(RAND_MAX)) > trainSilP)))
            {
              statistics.nrOfTrainOcc++;
              int observationLength = trainingPool->getCurSegmentLen(&readPool_curSeg);
              int time              = trainingPool->getCurSegmentStart(&readPool_curSeg);
              int contextKey        = trainingPool->getSegmentID(TRAIN_CONTEXTKEY,&readPool_curSeg);
              int stateNr           = stateMix_1[contextKey];

              mixtureSetData[stateNr].transitionP_toNext += observationLength;
              nrOfSamplesPerCluster[stateNr] ++;
              int cnt   = 0;
              int mxIn  = -1;
              int mxOut = -1;
              if(trainWithoutBorders)
              {
                int length = 50;
                mxIn  = length;
                mxOut = length;

                if(readPool_curSeg.prevSeg == NULL || (time - readPool_curSeg.prevSeg->lastFrame > length || readPool_curSeg.prevSeg->ID[0] == readPool_curSeg.curSeg->ID[0]))
                {
                  mxIn = 0;
                }
                if(readPool_curSeg.curSeg->next != NULL &&
                  (readPool_curSeg.curSeg->next->firstFrame - (time + observationLength) > length || readPool_curSeg.curSeg->next->ID[0] == readPool_curSeg.curSeg->ID[0]))
                {
                  mxOut = 0;
                }
              }

              if(!trainWithoutBorders || observationLength > mxIn + mxOut)
              {
                if(inFilter)
                {
                  count++;
                  double weight = 1.0;
                  if(trainWithoutBorders && ((cnt < mxIn) || (observationLength - cnt < mxOut)))
                  {
                    weight = 0.3;
                  }
                  cnt++;
                  if(trainDiscr != NULL)
                  {
                    weight = trainDiscr[time]->multiplyVector(trainDiscrMask) / trainDiscr[time]->getValue(0);
                  }
                  if(weight > 0.0)
                  {
                    mixtureSetData[stateNr].state->train(false,weight,trainVector);
                  }
                }
                while(!isLast)
                {
                  trainVector = trainingPool->getNextVector(&readPool_curSeg,&isLast, &inFilter);
                  time++;
                  if(inFilter)
                  {
                    count++;
                    double weight = 1.0;
                    if(trainWithoutBorders && ((cnt < mxIn) || (observationLength - cnt < mxOut)))
                    {
                      weight = 0.3;
                    }
                    cnt++;
                    if(trainDiscr != NULL)
                    {
                      weight = trainDiscr[time]->multiplyVector(trainDiscrMask) / trainDiscr[time]->getValue(0);
                    }
                    if(weight > 0.0)
                    {
                      mixtureSetData[stateNr].state->train(false,weight,trainVector);
                    }
                  }
                }
              }
            }
            trainVector = trainingPool->getFirstVectorNextSegment(&readPool_curSeg,useLabel < 0, &isLast, &inFilter);
            vectorNotNull = (trainVector != NULL);
          }
          if((useLabel != guestTrainingLabel) && (guestTrainingLabel >= 0))
          {
            useLabel = guestTrainingLabel;
            proceed  = true;
          }
          else
          {
            proceed  = false;
          }
        }

        for(int a=0;a<statistics.nrOfContexts;a++)
        {
          mixtureSetData[a].transitionP_toNext /= ((double)nrOfSamplesPerCluster[a]);
          double transSelf = 1.0 - 1.0/mixtureSetData[a].transitionP_toNext;
          if(transSelf < TRANSITION_MINIMUM)
          {
            PRINTF("# WARNING: INVALID TRANSITION VALUE!\n");
            mixtureSetData[a].transitionP_toSelf = log(TRANSITION_MINIMUM);
            mixtureSetData[a].transitionP_toNext = log(1.0 - TRANSITION_MINIMUM);
          }
          else
          {
            mixtureSetData[a].transitionP_toSelf = log(transSelf);
            mixtureSetData[a].transitionP_toNext = log(1.0 / mixtureSetData[a].transitionP_toNext);
          }
          mixtureSetData[a].state->trainFinish(neverPrune);
        }
        if(!doFastTraining)
        {
          statistics.likelihood = getTrainSilP();
          if(guestTrainingLabel >= 0)
          {
            statistics.likelihood += getTrainSilP(guestTrainingLabel);
          }
        }
      }

      int r = -1;
      int r2;
      for(int a=0;a<statistics.nrOfContexts;a++)
      {
        r2 = mixtureSetData[a].state->getNumberOfGaussians();
        if(r2>r)
        {
          r = r2;
        }
      }
      statistics.nrOfGaussians = r;
      nrOfG = statistics.nrOfGaussians;

      PRINTF6("# %s (SINGLE-STATE,%d) Training in progress: %d - %12f (%d)\n",statistics.name,pass2, nrOfG,statistics.likelihood,count);

//      writeModel(modelFile);

      // If still allowed, double the gaussians for each state and train again!
      // (Because gaussians may be deleted as well, check if the number of gaussians is bigger than last time..
      if(lastNrOfGaussians < nrOfG && nrOfG < maxGaussians)
      {
        PRINTF("Start splitting procedure...\n");
        lastNrOfGaussians = nrOfG;
        int r = -1;
        int r2;
        for(int a=0;a<statistics.nrOfContexts;a++)
        {
          if(doFastTraining)
          {
            mixtureSetData[a].state->splitAllGaussians(true);
          }
          else
          {
            mixtureSetData[a].state->splitBestGaussian();
          }
          r2 = mixtureSetData[a].state->getNumberOfGaussians();
          if(r2>r)
          {
            r = r2;
          }
        }
        statistics.nrOfGaussians = r;
        nrOfG = statistics.nrOfGaussians;
      }
      else
      {
        proceedDoubling = false;
      }
    }
  }
  else    ////////// ----------> This is NOT a SIL model:
  {
//    printf("Going to train a TRI-PHONE model.. (%d) \n", nrOfG);
    // Do we need to initialise the model?
    if(nrOfG < 2)
    {
      if(decision_Matrix != NULL)
      {
        double stateMixL[statistics.maxNrOfContexts];
        int    totNrS[statistics.maxNrOfContexts];
        for(int i=0;i<statistics.maxNrOfContexts;i++)
        {
          stateMixL[i] = 0.0;
          totNrS[i]    = 0;
        }
        int *stateMix        = NULL;
        int clusterNr        = 0;
        int numberOfClusters = 0;
        int clusterID        = 1;
        for(int i=0;i<3;i++)
        {
          int clusterNrState = 0;
          PRINTF4("# Calculating phone cluster %03d  (%d,%d)...\n",i+1,numberOfClusters,clusterNr);
          if(i==0)
          {
            stateMix  = stateMix_1;
            clusterID = TRAIN_CLUSTER_ID1;
          }
          else if(i==1)
          {
            stateMix = stateMix_2;
            clusterID = TRAIN_CLUSTER_ID2;
          }
          else
          {
            stateMix = stateMix_3;
            clusterID = TRAIN_CLUSTER_ID3;
          }
          numberOfClusters++;
          for(int a=0;a<decision_numberOfModels;a++)
          {
            for(int b=0;b<decision_numberOfModels;b++)
            {
              stateMix[a*decision_numberOfModels+b] = clusterNr;
            }
          }

          // Okay, let's populate our 'context gaussian' list:
          Gaussian gClust[statistics.maxNrOfContexts];

          bool isLast = false;
          Vector *trainVector = trainingPool->getFirstVectorFirstSegment(&readPool_curSeg,clusterID,trainingLabel,trainingLabel < 0, &isLast);
          while(trainVector != NULL)
          {
            int contextKey = trainingPool->getSegmentID(TRAIN_CONTEXTKEY,&readPool_curSeg);
            gClust[contextKey].train(1.0,trainVector);
            while(!isLast)
            {
              trainVector = trainingPool->getNextVector(&readPool_curSeg,&isLast);
              gClust[contextKey].train(1.0,trainVector);
            }
            trainVector = trainingPool->getFirstVectorNextSegment(&readPool_curSeg,trainingLabel < 0, &isLast);
          }

          // Calculate the overall score first:
          double totL = 0.0;
          Gaussian gTot;

          for(int aa=0;aa<statistics.maxNrOfContexts;aa++)
          {
            if(stateMix[aa] == clusterNr)
            {
              gTot.train(&(gClust[aa]));
            }
          }
          gTot.trainFinish();

          isLast = false;
          trainVector = trainingPool->getFirstVectorFirstSegment(&readPool_curSeg,clusterID,trainingLabel,trainingLabel < 0, &isLast);
          while(trainVector != NULL)
          {
            int contextKey = trainingPool->getSegmentID(TRAIN_CONTEXTKEY,&readPool_curSeg);

            if(stateMix[contextKey] == clusterNr)
            {
              totNrS[clusterNr]++;
              totL += gTot.getLogP(trainVector);
              while(!isLast)
              {
                trainVector = trainingPool->getNextVector(&readPool_curSeg,&isLast);
                totL += gTot.getLogP(trainVector);
              }
            }
            trainVector = trainingPool->getFirstVectorNextSegment(&readPool_curSeg,trainingLabel < 0, &isLast);
          }
          stateMixL[clusterNr] = totL;

          // Start splitting the first cluster.
          // keep on it until the final cluster is reached...
          while((clusterNr < numberOfClusters) && (clusterNrState < MAXIMUM_NUMBER_OF_STATE_CLUSTERS))
          {
            // First calculate the total likelihood of this cluster:
            double bestL      = stateMixL[clusterNr];
            double bestL1     = 0.0;
            double bestL2     = 0.0;
            int    bestRule   = -1;
            bool   bestIsLeft = true;
            int    bestNrOfSamples1 = 0;
            int    bestNrOfSamples2 = 0;

            for(int ruleL=0;ruleL<decision_numberOfRules;ruleL++)
            {
              Gaussian g1;
              Gaussian g2;
              int nrOfSamples1 = 0;
              int nrOfSamples2 = 0;

              for(int aL=0;aL<decision_numberOfModels;aL++)
              {
                int ai = aL*decision_numberOfModels;
                for(int aR=0;aR<decision_numberOfModels;aR++)
                {
                  if(stateMix[ai] == clusterNr)
                  {
                    if(decision_Matrix[ruleL*decision_numberOfModels+aL] == 1)
                    {
                      g1.train(&(gClust[ai]));
                    }
                    else
                    {
                      g2.train(&(gClust[ai]));
                    }
                  }
                  ai++;
                }
              }
              g1.trainFinish();
              g2.trainFinish();

              Vector *trainVector = trainingPool->getFirstVectorFirstSegment(&readPool_curSeg,clusterID,trainingLabel,trainingLabel < 0);
              while(trainVector != NULL)
              {
                int contextKey = trainingPool->getSegmentID(TRAIN_CONTEXTKEY,&readPool_curSeg);
                if(stateMix[contextKey] == clusterNr)
                {
                  int contextLeft = trainingPool->getSegmentID(TRAIN_CONTEXTLEFT,&readPool_curSeg);
                  if(decision_Matrix[ruleL*decision_numberOfModels+contextLeft] == 1)
                  {
                    nrOfSamples1++;
                  }
                  else
                  {
                    nrOfSamples2++;
                  }
                }
                trainVector = trainingPool->getFirstVectorNextSegment(&readPool_curSeg,trainingLabel < 0);
              }

              if(nrOfSamples1 > MIN_SAMPLES_PER_CLUSTER && nrOfSamples2 > MIN_SAMPLES_PER_CLUSTER)
              {
                double currentL1 = 0.0;
                double currentL2 = 0.0;
                bool isLast      = false;
                Vector *trainVector = trainingPool->getFirstVectorFirstSegment(&readPool_curSeg,clusterID,trainingLabel,trainingLabel < 0,&isLast);
                while(trainVector != NULL)
                {
                  int contextKey = trainingPool->getSegmentID(TRAIN_CONTEXTKEY,&readPool_curSeg);
                  if(stateMix[contextKey] == clusterNr)
                  {
                    int contextLeft = trainingPool->getSegmentID(TRAIN_CONTEXTLEFT,&readPool_curSeg);
                    if(decision_Matrix[ruleL*decision_numberOfModels+contextLeft] == 1)
                    {
                      currentL1 += g1.getLogP(trainVector);
                      while(!isLast)
                      {
                        trainVector = trainingPool->getNextVector(&readPool_curSeg,&isLast);
                        currentL1 += g1.getLogP(trainVector);
                      }
                    }
                    else
                    {
                      currentL2 += g2.getLogP(trainVector);
                      while(!isLast)
                      {
                        trainVector = trainingPool->getNextVector(&readPool_curSeg,&isLast);
                        currentL2 += g2.getLogP(trainVector);
                      }
                    }
                  }
                  trainVector = trainingPool->getFirstVectorNextSegment(&readPool_curSeg,trainingLabel < 0, &isLast);
                }
                if(currentL1 + currentL2 > bestL)
                {
                  bestL            = currentL1 + currentL2;
                  bestL1           = currentL1;
                  bestL2           = currentL2;
                  bestRule         = ruleL;
                  bestIsLeft       = true;
                  bestNrOfSamples1 = nrOfSamples1;
                  bestNrOfSamples2 = nrOfSamples2;
                }
              }
            }

            for(int ruleR=0;ruleR<decision_numberOfRules;ruleR++)
            {
              Gaussian g1;
              Gaussian g2;
              int nrOfSamples1 = 0;
              int nrOfSamples2 = 0;

              for(int aL=0;aL<decision_numberOfModels;aL++)
              {
                int ai = aL*decision_numberOfModels;
                for(int aR=0;aR<decision_numberOfModels;aR++)
                {
                  if(stateMix[ai] == clusterNr)
                  {
                    if(decision_Matrix[ruleR*decision_numberOfModels+aR] == 1)
                    {
                      g1.train(&(gClust[ai]));
                    }
                    else
                    {
                      g2.train(&(gClust[ai]));
                    }
                  }
                  ai++;
                }
              }
              g1.trainFinish();
              g2.trainFinish();

              Vector *trainVector = trainingPool->getFirstVectorFirstSegment(&readPool_curSeg,clusterID,trainingLabel,trainingLabel < 0);
              while(trainVector != NULL)
              {
                int contextKey = trainingPool->getSegmentID(TRAIN_CONTEXTKEY,&readPool_curSeg);
                if(stateMix[contextKey] == clusterNr)
                {
                  int contextRight = trainingPool->getSegmentID(TRAIN_CONTEXTRIGHT,&readPool_curSeg);
                  if(decision_Matrix[ruleR*decision_numberOfModels+contextRight] == 1)
                  {
                    nrOfSamples1++;
                  }
                  else
                  {
                    nrOfSamples2++;
                  }
                }
                trainVector = trainingPool->getFirstVectorNextSegment(&readPool_curSeg,trainingLabel < 0);
              }

              if(nrOfSamples1 > MIN_SAMPLES_PER_CLUSTER && nrOfSamples2 > MIN_SAMPLES_PER_CLUSTER)
              {
                double currentL1 = 0.0;
                double currentL2 = 0.0;
                bool isLast      = false;
                Vector *trainVector = trainingPool->getFirstVectorFirstSegment(&readPool_curSeg,clusterID,trainingLabel,trainingLabel < 0,&isLast);
                while(trainVector != NULL)
                {
                  int contextKey = trainingPool->getSegmentID(TRAIN_CONTEXTKEY,&readPool_curSeg);
                  if(stateMix[contextKey] == clusterNr)
                  {
                    int contextRight = trainingPool->getSegmentID(TRAIN_CONTEXTRIGHT,&readPool_curSeg);
                    if(decision_Matrix[ruleR*decision_numberOfModels+contextRight] == 1)
                    {
                      currentL1 += g1.getLogP(trainVector);
                      while(!isLast)
                      {
                        trainVector = trainingPool->getNextVector(&readPool_curSeg,&isLast);
                        currentL1 += g1.getLogP(trainVector);
                      }
                    }
                    else
                    {
                      currentL2 += g2.getLogP(trainVector);
                      while(!isLast)
                      {
                        trainVector = trainingPool->getNextVector(&readPool_curSeg,&isLast);
                        currentL2 += g2.getLogP(trainVector);
                      }
                    }
                  }
                  trainVector = trainingPool->getFirstVectorNextSegment(&readPool_curSeg,trainingLabel < 0, &isLast);
                }
                if(currentL1 + currentL2 > bestL)
                {
                  bestL            = currentL1 + currentL2;
                  bestL1           = currentL1;
                  bestL2           = currentL2;
                  bestRule         = ruleR;
                  bestIsLeft       = false;
                  bestNrOfSamples1 = nrOfSamples1;
                  bestNrOfSamples2 = nrOfSamples2;
                }
              }
            }
            // Update the clusters!
            double perc = (bestL - stateMixL[clusterNr])/((double)totNrS[clusterNr]);
            if(bestRule >= 0 && (perc > MINIMUM_CLUSTER_IMPROVEMENT))
            {
              clusterNrState++;
              PRINTF2("New splitting rule for cluster %03d: Rule ",clusterNr);
              if(bestIsLeft)
              {
                PRINTF("L-");
              }
              else
              {
                PRINTF("R-");
              }
              PRINTF5("%03d\t(%07.3f) (%d,%d)\n", bestRule, perc,bestNrOfSamples1,bestNrOfSamples2);

              stateMixL[clusterNr]        = bestL2;
              stateMixL[numberOfClusters] = bestL1;
              totNrS[clusterNr]           = bestNrOfSamples2;
              totNrS[numberOfClusters]    = bestNrOfSamples1;

              for(int a=0;a<decision_numberOfModels;a++)
              {
                for(int b=0;b<decision_numberOfModels;b++)
                {
                  int index = a*decision_numberOfModels+b;
                  if(stateMix[index] == clusterNr)
                  {
                    if(( bestIsLeft && decision_Matrix[bestRule*decision_numberOfModels+a]==1) ||
                       (!bestIsLeft && decision_Matrix[bestRule*decision_numberOfModels+b]==1)    )
                    {
                      stateMix[index] = numberOfClusters;
                    }
                  }
                }
              }
              numberOfClusters++;
            }
            else            // Choose new cluster number:
            {
              clusterNr++;
            }
          }
          clusterNr = numberOfClusters;
        }
        statistics.nrOfContexts = numberOfClusters;
        PRINTF2("------- DONE CLUSTERING: Number of clusters = %d --------\n", numberOfClusters);
        delete[] mixtureSetData;
        mixtureSetData = new MixtureSet[statistics.nrOfContexts];
        for(int ai=0;ai<statistics.nrOfContexts;ai++)
        {
          mixtureSetData[ai].state = new MixGaussian();
        }
      }


      int initTransitionsLength[statistics.nrOfContexts];
      int nrOfSamplesPerCluster[statistics.nrOfContexts];
      for(int a=0;a<statistics.nrOfContexts;a++)
      {
        initTransitionsLength[a] = 0;
        nrOfSamplesPerCluster[a] = 0;
      }


      if(decision_Matrix != NULL)  // Initialise with pre aligned phones:
      {
        bool isUsed[statistics.nrOfContexts];
        for(int i=0;i<statistics.nrOfContexts;i++)
        {
          isUsed[i] = false;
        }
        int clusterID = 1;
        PRINTF("Initialising clusters using the decision matrix...\n");
        for(int i=0;i<3;i++)
        {
          int *stateMix;
          if(i==0)
          {
            stateMix = stateMix_1;
            clusterID = TRAIN_CLUSTER_ID1;
          }
          else if(i==1)
          {
            stateMix = stateMix_2;
            clusterID = TRAIN_CLUSTER_ID2;
          }
          else
          {
            stateMix = stateMix_3;
            clusterID = TRAIN_CLUSTER_ID3;
          }
          bool isLast = false;
          Vector *trainVector = trainingPool->getFirstVectorFirstSegment(&readPool_curSeg,clusterID,trainingLabel,trainingLabel < 0,&isLast);
          while(trainVector != NULL)
          {
            int contextKey        = trainingPool->getSegmentID(TRAIN_CONTEXTKEY,&readPool_curSeg);
            int observationLength = trainingPool->getCurSegmentLen(&readPool_curSeg);

            int stateNr = stateMix[contextKey];
            isUsed[stateNr] = true;
            initTransitionsLength[stateNr] += observationLength;
            nrOfSamplesPerCluster[stateNr]++;

            mixtureSetData[stateNr].state->train(true,1.0,trainVector);
            while(!isLast)
            {
              trainVector = trainingPool->getNextVector(&readPool_curSeg,&isLast);
              mixtureSetData[stateNr].state->train(true,1.0,trainVector);
            }
            trainVector = trainingPool->getFirstVectorNextSegment(&readPool_curSeg,trainingLabel < 0, &isLast);
          }
        }
      }
      else                        // Initialise with the training set:
      {
        PRINTF("Initialising using the entire training set...\n");

        stateMix_1[0] = 0;
        stateMix_2[0] = 1;
        stateMix_3[0] = 2;

        bool isLast = false;
        Vector *trainVector = trainingPool->getFirstVectorFirstSegment(&readPool_curSeg,trainingSegment,trainingLabel,trainingLabel < 0,&isLast);
        while(trainVector != NULL)
        {
          int observationLength = trainingPool->getCurSegmentLen(&readPool_curSeg);
          int contextKey        = trainingPool->getSegmentID(TRAIN_CONTEXTKEY,&readPool_curSeg);
          int stateNr1 = stateMix_1[contextKey];
          int stateNr2 = stateMix_2[contextKey];
          int stateNr3 = stateMix_3[contextKey];
          initTransitionsLength[stateNr1] += (observationLength/3);
          initTransitionsLength[stateNr2] += (observationLength/3);
          initTransitionsLength[stateNr3] += (observationLength/3);
          nrOfSamplesPerCluster[stateNr1]++;
          nrOfSamplesPerCluster[stateNr2]++;
          nrOfSamplesPerCluster[stateNr3]++;

          mixtureSetData[stateNr1].state->train(true,1.0,trainVector);
          mixtureSetData[stateNr2].state->train(true,1.0,trainVector);
          mixtureSetData[stateNr3].state->train(true,1.0,trainVector);
          while(!isLast)
          {
            trainVector = trainingPool->getNextVector(&readPool_curSeg,&isLast);
            mixtureSetData[stateNr1].state->train(true,1.0,trainVector);
            mixtureSetData[stateNr2].state->train(true,1.0,trainVector);
            mixtureSetData[stateNr3].state->train(true,1.0,trainVector);
          }
          trainVector = trainingPool->getFirstVectorNextSegment(&readPool_curSeg,trainingLabel < 0, &isLast);
        }
      }

      for(int a=0;a<statistics.nrOfContexts;a++)
      {
        mixtureSetData[a].state->trainFinish();
        double trans = ((double)initTransitionsLength[a]) / ((double)nrOfSamplesPerCluster[a]);
        mixtureSetData[a].transitionP_toSelf = log( ((double)trans-1.0)/((double)trans) );
        mixtureSetData[a].transitionP_toNext = log(   1.0              /((double)trans) );
      }

      PRINTF("Initialised...\n");
    }

    //////////// Now training the model:

    // Loop for doubling the gaussians
    if(doFastTraining)
    {
      viterbi(4);
      bool isLast       = false;
      bool inFilter     = false;
      lastNrOfGaussians = -1;
      proceedDoubling   = true;
      while(proceedDoubling)
      {
        for(int pass=0;pass<3;pass++)
        {
          PRINTF2("Pass %d...\n",pass);
          for(int tr=0;tr<3;tr++)
          {
            int trID = TRAIN_CLUSTER_ID1;
            if(tr==1)
            {
              trID = TRAIN_CLUSTER_ID2;
            }
            else if(tr == 2)
            {
              trID = TRAIN_CLUSTER_ID3;
            }
            Vector *trainVector = trainingPool->getFirstVectorFirstSegment(&readPool_curSeg,trID,trainingLabel,trainingLabel < 0, &isLast, &inFilter);
            statistics.nrOfTrainOcc = 0;
            while(trainVector != NULL)
            {
              statistics.nrOfTrainOcc++;
              int contextKey = trainingPool->getSegmentID(TRAIN_CONTEXTKEY,&readPool_curSeg);
              int state[3];
              state[0] = stateMix_1[contextKey];
              state[1] = stateMix_2[contextKey];
              state[2] = stateMix_3[contextKey];
              mixtureSetData[state[tr]].state->train(false,1.0,trainVector);
              while(!isLast)
              {
                trainVector = trainingPool->getNextVector(&readPool_curSeg,&isLast, &inFilter);
                mixtureSetData[state[tr]].state->train(false,1.0,trainVector);
              }
              trainVector = trainingPool->getFirstVectorNextSegment(&readPool_curSeg,trainingLabel < 0, &isLast, &inFilter);
            }
          }
          for(int a=0;a<statistics.nrOfContexts;a++)
          {
            mixtureSetData[a].state->trainFinish();
          }
        }
        statistics.likelihood          = likelihood_new;
        statistics.nrOfGaussians       = maxNrOfGaussians();
        nrOfG                          = statistics.nrOfGaussians;
//        printf("Should we split? lastNrOfGaussians = %d, nrOfG = %d, maxGaussians = %d\n",lastNrOfGaussians,nrOfG,maxGaussians);
        if(lastNrOfGaussians < nrOfG && nrOfG < maxGaussians)
        {
          lastNrOfGaussians = nrOfG;
          for(int a=0;a<statistics.nrOfContexts;a++)
          {
            mixtureSetData[a].state->splitBestGaussian();
          }
          statistics.nrOfGaussians = maxNrOfGaussians();
          nrOfG = statistics.nrOfGaussians;
        }
        else
        {
          proceedDoubling = false;
        }
      }
      PRINTF("Starting final Baum-Welch training..\n");
      for(int i=0;i<3;i++)
      {
        likelihood_new   = baumWelch(2);
        PRINTF2(" ... score = %f\n",likelihood_new);
        likelihood_new   = baumWelch(1);
        PRINTF2(" ... score = %f\n",likelihood_new);
        likelihood_new   = baumWelch(2);
        PRINTF2(" ... score = %f\n",likelihood_new);
        likelihood_new   = baumWelch(2);
      }
      statistics.likelihood          = likelihood_new;
      statistics.nrOfGaussians       = maxNrOfGaussians();
      nrOfG                          = statistics.nrOfGaussians;

    }
    else if(!doSat)
    {
      proceedDoubling = true;
      lastNrOfGaussians = -1;
      int baumIteration = 0;
      while(proceedDoubling)
      {
        // Ready to start the training until the likelihood does not improve enough
        likelihood_new   = viterbi(2);
        PRINTF3("Viterbi training value: %d - %f\n",nrOfG,likelihood_new);
        fflush(stdout);

        bool pass = true;
        while(pass || ((likelihood_new - likelihood_old) / absolute(likelihood_old) > MINIMUM_TRAIN_IMPROVEMENT))
        {
          pass                           = false;
          likelihood_old                 = likelihood_new;
          likelihood_new                 = viterbi(2);
          likelihood_new                 = viterbi(1);
          likelihood_new                 = viterbi(2);
          statistics.likelihood          = likelihood_new;
          statistics.nrOfGaussians       = maxNrOfGaussians();
          nrOfG                          = statistics.nrOfGaussians;
          PRINTF3("Viterbi training in progress: %d - %f\n",nrOfG,likelihood_new);
          fflush(stdout);
        }

        baumIteration++;
        if(baumIteration > TRAIN_BAUM_WELCH_FREQUENCY)
        {
          baumIteration = 0;
          // Refine with Baum-Welch:
          likelihood_new   = baumWelch(2);
          PRINTF3("Baum-Welch initial training value: %d - %f\n",nrOfG,likelihood_new);
          fflush(stdout);
          statistics.likelihood = likelihood_new;
          bool pass = true;
          while(pass || ((likelihood_new - likelihood_old) / absolute(likelihood_old) > MINIMUM_TRAIN_IMPROVEMENT))
          {
            pass             = false;
            likelihood_old   = likelihood_new;
            likelihood_new   = baumWelch(2);
            likelihood_new   = baumWelch(1);
            likelihood_new   = baumWelch(2);
            statistics.likelihood          = likelihood_new;
            statistics.nrOfGaussians       = maxNrOfGaussians();
            nrOfG                          = statistics.nrOfGaussians;
            PRINTF3("Baum-Welch training in progress: %d - %f\n",nrOfG,likelihood_new);
            fflush(stdout);
          }
        }
        // If still allowed, increase the gaussians for each state and train again!
        // (Because gaussians may be deleted as well, check if the number of gaussians is bigger than last time..
//        printf("Should we split? lastNrOfGaussians = %d, nrOfG = %d, maxGaussians = %d\n",lastNrOfGaussians,nrOfG,maxGaussians);
        if(lastNrOfGaussians < nrOfG && nrOfG < maxGaussians)
        {
          lastNrOfGaussians = nrOfG;
          for(int a=0;a<statistics.nrOfContexts;a++)
          {
            mixtureSetData[a].state->splitBestGaussian();
          }
          statistics.nrOfGaussians = maxNrOfGaussians();
          nrOfG = statistics.nrOfGaussians;
        }
        else
        {
          proceedDoubling = false;
        }
      }
      for(int i=0;i<6;i++)
      {
        // Final Refine with Baum-Welch:
        likelihood_new   = baumWelch(1);
        likelihood_new   = baumWelch(2);
        likelihood_new   = baumWelch(2);
        PRINTF3("Baum-Welch final training initial value: %d - %f\n",nrOfG,likelihood_new);
        fflush(stdout);
        statistics.likelihood = likelihood_new;
        bool pass = true;
        while(pass || ((likelihood_new - likelihood_old) / absolute(likelihood_old) > MINIMUM_TRAIN_IMPROVEMENT_FINAL))
        {
          pass             = false;
          likelihood_old   = likelihood_new;
          likelihood_new   = baumWelch(1);
          likelihood_new   = baumWelch(2);
          likelihood_new   = baumWelch(2);
          statistics.likelihood          = likelihood_new;
          statistics.nrOfGaussians       = maxNrOfGaussians();
          nrOfG                          = statistics.nrOfGaussians;
          PRINTF3("Final Baum-Welch training in progress: %d - %f\n",nrOfG,likelihood_new);
          fflush(stdout);
        }
      }
    }
    else if(doSat)
    {
      baumWelch(3,doSat);
    }
    PRINTF2("Training done: %f\n",likelihood_new);
    fflush(stdout);
  }
  return statistics.likelihood;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainPhoneModel::doNotuseBordersForTraining(bool useBordersNot)
{
  trainWithoutBorders = useBordersNot;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainPhoneModel::moveModelGaussians(TrainPhoneModel *model, double factor)
{
  // Only look at the first context!
  mixtureSetData[0].state->moveModelGaussians(model->mixtureSetData[0].state,factor);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Performs the adaptation training-run. This run will set all acumulators. You can only
/// use other segments/labels/pools for sil models!
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainPhoneModel::adapt_setAcTrain(int useLabel, int useSegmentation, FeaturePool *usePool)
{
  if(usePool == NULL)
  {
    usePool = trainingPool;
  }
  if(useLabel < 0)
  {
    useLabel      = trainingLabel;
  }
  if(useSegmentation < 0)
  {
    useSegmentation = trainingSegment;
  }

  if(isSil)  // This is SIL!
  {
    adapt_setAcumulators(useLabel,useSegmentation,usePool);
  }
  else  // Not SIL..
  {
    assert(usePool == trainingPool);
    assert(useLabel == trainingLabel);
    assert(useSegmentation == trainingSegment);
    baumWelch(3);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

double TrainPhoneModel::getCoSim(TrainPhoneModel *t1, TrainPhoneModel *t2)
{
  // Only look at the first context!
  return mixtureSetData[0].state->getCoSim(t1->mixtureSetData[0].state,t2->mixtureSetData[0].state);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainPhoneModel::addGaussian(Vector *v)
{
  mixtureSetData[0].state->addGaussian(v);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

double TrainPhoneModel::getKLDistance(TrainPhoneModel *t2)
{
  // Only look at the first context!
  return mixtureSetData[0].state->getKLDistance(t2->mixtureSetData[0].state);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

double TrainPhoneModel::getNormDistance()
{
  // Only look at the first context!
  return mixtureSetData[0].state->getNormDistance();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Add the best nmbr gaussians from the model source to this model. The best gaussians are
/// determined during earlier counting by count().
/// (mixtureset 0, we expect this model to be SIL without context (used by Train_Speaker_Segmenter).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainPhoneModel::addCountedGaussians(TrainPhoneModel *source, int nmbr)
{
  mixtureSetData[0].state->addCountedGaussians(source->mixtureSetData[0].state,nmbr);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Starts a counting run for determining the most important gaussians of this model
/// (mixtureset 0, we expect this model to be SIL without context (used by Train_Speaker_Segmenter).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainPhoneModel::startCount()
{
  mixtureSetData[0].state->startCount();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Stops the counting run for determining the most important gaussians of this model
/// (mixtureset 0, we expect this model to be SIL without context (used by Train_Speaker_Segmenter).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainPhoneModel::stopCount()
{
  mixtureSetData[0].state->stopCount();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Counts the importance of each gaussian in the model given the feature Vector observation.
/// (mixtureset 0, we expect this model to be SIL without context (used by Train_Speaker_Segmenter).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainPhoneModel::count(Vector *observation)
{
  mixtureSetData[0].state->count(observation);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Determines the dominant Gaussian in the model...
/////////////////////////////////////////////////////////////////////////////////////////////////////

int TrainPhoneModel::getDominantGaussian()
{
  return mixtureSetData[0].state->getBestCount();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainPhoneModel::normalize()
{
  Gaussian g(1);
  Vector v(1);
  v.setValue(0,-5.0);
  g.train(1.0,&v);
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    PRINTF2("%f\n",mixtureSetData[i].transitionP_toNext);
    v.setValue(0,mixtureSetData[i].transitionP_toNext);
    g.train(1.0,&v);
  }
  g.trainFinish();
  double minValue = g.getMean()->getValue(0) - g.getVariance()->getValue(0)*0.5;
  PRINTF3("mean, min: %f  %f\n",g.getMean()->getValue(0),minValue);
  for(int i=0;i<statistics.nrOfContexts;i++)
  {
    if(mixtureSetData[i].transitionP_toNext < minValue)
    {
      mixtureSetData[i].transitionP_toNext = minValue;
    }
  }
}

