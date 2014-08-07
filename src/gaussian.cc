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
#include "gaussian.h"
#include <math.h>
#include "shout-misc.h"

using namespace WriteFileLittleBigEndian;

#define GAUSSIAN_SHIFT_FACTOR  (0.3)

#define MAP_ADAPT_GAUSSIAN_RELEVANCE_FACTOR  (10.0)
#define MMI_LAMBDA  (0.05)

/// \todo This was ment to be 'code efficient' but I need to get rid of it...
int AUGMENTED_MATRIX_SIZE    =     (63960);
int AUG_X                    =      (1640);
int AUG_Y                    =        (41);
int AUG_Y_TIMES_VECTOR       =      (1599);
int AUG_P_PLUS_ONE           =      (1639);
int AUG_P                    =      (1638);


bool LOG_ON = false;


void Gaussian::printInfo(Vector *v)
{
  v->addElements(true,meanVector);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::printModel(FILE *fileMean, FILE *fileVariance)
{
  for(int i=0;i<meanVector->len();i++)
  {
    fprintf(fileMean,"%g ",meanVector->getValue(i));
    fprintf(fileVariance,"%g ",sqrt(varianceVector->getValue(i)));
  }
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The constructor initialises the gaussian. The variance is set to 100 (very big), and the
/// mean is set to all zero.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Gaussian::Gaussian(int dimensions)
{
  mmi_trainPars     = NULL;
  trainPars         = NULL;
  varianceVector    = new Vector(dimensions);
  meanVector        = new Vector(dimensions);
  varianceVector->setAllValues(100.0);
  meanVector->setAllValues(0.0);
  calculateCP();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor creates an identical copy from another Gaussian.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Gaussian::Gaussian(Gaussian *copyFrom, bool shiftGaussian)
{
  mmi_trainPars     = NULL;
  trainPars         = NULL;
  varianceVector    = copyFrom->varianceVector->copyVector();
  meanVector        = copyFrom->meanVector->copyVector();
  if(shiftGaussian)
  {
//    Vector *v = varianceVector->multiplyElements(false,GAUSSIAN_SHIFT_FACTOR);
    Vector *v = varianceVector->sqrtElements(false);
    v->getMaxElement(true);
    v->multiplyElements(true,GAUSSIAN_SHIFT_FACTOR);
    meanVector->addElements(true,v);
    delete v;
  }
  calculateCP();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor creates a Gaussian with mean inbetween the mean of model1 and model2.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Gaussian::Gaussian(Gaussian *model1, Gaussian *model2, double rate)
{
  mmi_trainPars     = NULL;
  trainPars         = NULL;
  varianceVector    = model1->varianceVector->copyVector();
  meanVector        = model1->meanVector->copyVector();
  meanVector->multiplyElements(true,rate);
  rate = 1.0-rate;
  Vector *m2        = model2->meanVector->copyVector();
  m2->multiplyElements(true,rate);
  meanVector->addElements(true,m2);
  delete m2;
  calculateCP();
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Reads the Gaussian parameters (mean and variance Vectors) from file.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Gaussian::Gaussian(FILE *inFile, int dim)
{
  mmi_trainPars     = NULL;
  trainPars         = NULL;
  meanVector        = new Vector(inFile,dim);
  varianceVector    = new Vector(inFile,dim);
  calculateCP();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor deletes the variance Vector and the mean Vector.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Gaussian::~Gaussian()
{
  delete varianceVector;
  delete meanVector;
  if(trainPars != NULL)
  {
    if(trainPars->unAdaptedMeanVector != NULL)
    {
      delete trainPars->unAdaptedMeanVector;
    }
    delete trainPars->numeratorMean;
    delete trainPars->numeratorVariance;
    delete trainPars;
  }
  if(mmi_trainPars != NULL)
  {
    if(mmi_trainPars->unAdaptedMeanVector != NULL)
    {
      delete mmi_trainPars->unAdaptedMeanVector;
    }
    delete mmi_trainPars->numeratorMean;
    delete mmi_trainPars->numeratorVariance;
    delete mmi_trainPars;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

double Gaussian::getCoSim(Gaussian *g1, Gaussian *g2)
{
  Vector *v1 = meanVector->substractElements(false,g1->meanVector);
  Vector *v2 = meanVector->substractElements(false,g2->meanVector);
  Vector *v0 = new Vector(meanVector->len());
  v0->setAllValues(0.0);
  double res = v1->multiplyVector(v2) / (sqrt(v0->distance(v1))*sqrt(v0->distance(v2)));
  delete v0;
  delete v1;
  delete v2;
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

double Gaussian::getKLDistance(Gaussian *g2)
{
  double res  = 0.0;
  double dist = 0.0;
  for(int i=0;i<meanVector->len();i++)
  {
    dist = (meanVector->getValue(i) - g2->meanVector->getValue(i));
    res += (dist*dist / varianceVector->getValue(i));
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

double Gaussian::getNormDistance()
{
  double res  = 0.0;
  for(int i=0;i<meanVector->len();i++)
  {
    res += 1.0/varianceVector->getValue(i);
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// In order to safe memory, Gaussian training parameters are grouped in a TrainGaussian structure.
/// Each Gaussian object contains a NULL pointer (trainPars) of type TrainGaussian. When a training
/// method is called and the trainPars variable is NULL, the method enterTrainingMode() is automatically
/// called. This method will create the TrainGaussian structure and initialise its variables.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::enterTrainingMode(bool mmiTraining)
{
  int dim = meanVector->len();
  if(mmiTraining)
  {
    mmi_trainPars = new TrainGaussian;
    mmi_trainPars->adaptNode = NULL;
    mmi_trainPars->numeratorMean      = new Vector(dim);
    mmi_trainPars->numeratorVariance  = new Vector(dim);
    mmi_trainPars->numeratorMean->setAllValues(0.0);
    mmi_trainPars->numeratorVariance->setAllValues(0.0);
    mmi_trainPars->lastTrainingWeight  = 0.0;
    mmi_trainPars->trainingDenominator = 0.0;
    mmi_trainPars->trainingSamples     = 0;
    mmi_trainPars->unAdaptedMeanVector = NULL;
  }
  else
  {
    trainPars = new TrainGaussian;
    trainPars->adaptNode = NULL;
    trainPars->numeratorMean      = new Vector(dim);
    trainPars->numeratorVariance  = new Vector(dim);
    trainPars->numeratorMean->setAllValues(0.0);
    trainPars->numeratorVariance->setAllValues(0.0);
    trainPars->lastTrainingWeight  = 0.0;
    trainPars->trainingDenominator = 0.0;
    trainPars->trainingSamples     = 0;
    trainPars->unAdaptedMeanVector = NULL;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Load training data from file and append it to the training parameters.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::appendSAT(FILE *inFile)
{
  if(trainPars == NULL)
  {
    enterTrainingMode();
  }
  Vector *m = new Vector(inFile,trainPars->numeratorMean->len(),false);
  trainPars->numeratorMean->addElements(true,m);
  delete m;
  m = new Vector(inFile,trainPars->numeratorVariance->len(),false);
  trainPars->numeratorVariance->addElements(true,m);
  delete m;
  double d = 0.0;
  freadEndianSafe(&d,1,sizeof(d),inFile);
  trainPars->trainingDenominator += d;
  trainPars->lastTrainingWeight   = trainPars->trainingDenominator;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Writes training data to file.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::storeSAT(FILE *outFile)
{
  if(trainPars == NULL)
  {
    enterTrainingMode();
  }
  trainPars->numeratorMean->storeData(outFile);
  trainPars->numeratorVariance->storeData(outFile);
  fwrite(&(trainPars->trainingDenominator),1,sizeof(trainPars->trainingDenominator),outFile);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Writes the Gaussian parameters (mean and variance Vectors) to file.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::storeData(FILE *outFile)
{
  meanVector->storeData(outFile);
  varianceVector->storeData(outFile);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Sets the variance of the Gaussian
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::setVariance(const Vector *v)
{
  // the call to setAllValues() below will check to make sure
  // the the sizes of the vectors are the same, so don't check here...
  varianceVector->setAllValues(v);
  calculateCP();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Sets the mean of the Gaussian
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::setMean(Vector *v)
{
  // the call to setAllValues() below will check to make sure
  // the the sizes of the vectors are the same, so don't check here...
  meanVector->setAllValues(v);
  calculateCP();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Calculates the CP and logCP. The CP is 1 devided by the square root from the determinant of
/// the variance vector:
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::calculateCP(void)
{
  int len = varianceVector->len();
  CP = 0.0;
  for(int i=0;i<len;i++)
  {
    CP += log(TWO_PI*varianceVector->getValue(i));
  }
  CP*=-0.5;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the chance that an observation is produced by this Gaussian.
/////////////////////////////////////////////////////////////////////////////////////////////////////

double Gaussian::getP(Vector *observation)
{
  double res = observation->calculateEPart(meanVector,varianceVector)+CP;
  if(res < MINIMUM_LOGP)
  {
    res = MINIMUM_LOGP;
  }
  return exp(res);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the chance that an observation is produced by this Gaussian (in log)
/////////////////////////////////////////////////////////////////////////////////////////////////////

double Gaussian::getLogP(Vector *observation)
{
  double res = observation->calculateEPart(meanVector,varianceVector)+CP;
  if(res < MINIMUM_LOGP)
  {
    res = MINIMUM_LOGP;
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the training denominator. This is basically the number of training observations used
/// during the training phase.
/////////////////////////////////////////////////////////////////////////////////////////////////////

double Gaussian::getTrainingDenominator(void)
{
  if(trainPars == NULL)
  {
    enterTrainingMode();
  }
  return trainPars->lastTrainingWeight;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Shifts the mean Vector by -0.2*shiftFactor of the variance Vector. shiftFactor may be negative
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::shiftMean(int shiftFactor)
{
//  Vector *v = varianceVector->multiplyElements(false,-0.2*((double)shiftFactor));
  Vector *v = varianceVector->sqrtElements(false);
  v->getMaxElement(true);
  v->multiplyElements(true,-GAUSSIAN_SHIFT_FACTOR*((double)shiftFactor));
  meanVector->addElements(true,v);
  delete(v);
  calculateCP();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::moveModel(Gaussian *model, double factor)
{
  Vector *temp = model->meanVector->multiplyElements(false,factor);
  meanVector->addElements(true,temp);
  delete temp;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Train the Gaussian as follows:
///
/// The training data calculated for another Gaussian is added to the training data of this Gaussian.
/// This is used for quickly determine an optimal clustering of single Gaussian observation models.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::train(Gaussian *trainG)
{
  if(trainPars == NULL)
  {
    enterTrainingMode();
  }
  if(trainG->trainPars != NULL)
  {
    trainPars->trainingDenominator+=trainG->trainPars->trainingDenominator;
    trainPars->trainingSamples+=trainG->trainPars->trainingSamples;
    trainPars->lastTrainingWeight  =trainPars->trainingDenominator;
    trainPars->numeratorMean->addElements(true,trainG->trainPars->numeratorMean);
    trainPars->numeratorVariance->addElements(true,trainG->trainPars->numeratorVariance);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Train the Gaussian as follows:
///
/// The training factor (input parameter) is added to the denominator. When factor is one, this means
/// that an entire observation is added to the training pool. If 0.0 < factor < 1.0, the observation
/// is seen as being partly mapped on this gaussian.
///
/// Literature: Frederick Jelinek, "Statistical Methods for Speech Recognition": page 30 for
/// parameter (mean and variance) estimation.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::train(double factor, Vector *observation, Gaussian *doSat)
{
  if(trainPars == NULL)
  {
    enterTrainingMode();
  }
  if(factor > 0.0)
  {
    trainPars->trainingSamples++;
    trainPars->trainingDenominator+=factor;
    trainPars->lastTrainingWeight  =trainPars->trainingDenominator;
    // SAT training?
    // Mean:
    if(doSat != NULL)
    {
      Vector *o = observation->copyVector();
      Vector *sat = meanVector->substractElements(false,doSat->meanVector);
      o->substractElements(true,sat);
      delete sat;
      Vector *v = o->multiplyElements(false,factor);
      trainPars->numeratorMean->addElements(true,v);
      // Variance:
      v->multiplyElements(true,o);
      trainPars->numeratorVariance->addElements(true,v);
      delete v;
      delete o;
    }
    else
    {
      Vector *v = observation->multiplyElements(false,factor);
      trainPars->numeratorMean->addElements(true,v);
      // Variance:
      v->multiplyElements(true,observation);
      trainPars->numeratorVariance->addElements(true,v);
      delete v;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Used to finish the training procedure.
///
/// Literature: Frederick Jelinek, "Statistical Methods for Speech Recognition": page 30 for
/// parameter (mean and variance) estimation.
/////////////////////////////////////////////////////////////////////////////////////////////////////

double Gaussian::trainFinish(bool cleanUp, double minVar)
{
  double denom = 0.0;
  if(trainPars == NULL)
  {
    enterTrainingMode();
  }
  denom = trainPars->trainingDenominator;
  if(denom > 0.0)
  {
    /* First the new mean: */
    delete(meanVector);
    meanVector = trainPars->numeratorMean->multiplyElements(false,(double)(1.0 / denom));

    /* Now the new variance: */
    Vector *v = meanVector->multiplyElements(false,meanVector);
    delete(varianceVector);
    varianceVector = trainPars->numeratorVariance->multiplyElements(false,(double)(1.0 / denom));
    varianceVector->substractElements(true,v);
    varianceVector->minimizeElements(minVar);
    delete(v);
    /* Calculate new CP and trainingP: */
    trainPars->trainP = log(varianceVector->getDeterminant()*pow(TWO_PI,varianceVector->len()))*denom;
    calculateCP();
    /* Clean up: */
    if(cleanUp)
    {
      trainPars->numeratorVariance->setAllValues(0.0);
      trainPars->numeratorMean->setAllValues(0.0);
      trainPars->trainingDenominator = 0.0;
      trainPars->trainingSamples = 0;
    }
  }
  return denom;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::mapAdaptMean()
{
  double denom = 0.0;
  assert(trainPars != NULL);
  denom = trainPars->trainingDenominator;
  if(denom > 0.0)
  {
    double alfa = denom / (denom + (double)(MAP_ADAPT_GAUSSIAN_RELEVANCE_FACTOR));
    // Formula for standard MAP adaptation of the mean:
    // mean-new = alfa * numeratorMean / denom + (1 - alfa) * mean-old
    Vector *v  = trainPars->numeratorMean->multiplyElements(false,(alfa / (double)(denom)));
    Vector *v2 = meanVector->multiplyElements(false,(1.0 - alfa));
    delete(meanVector);
    meanVector = v->addElements(false,v2);
    delete v;
    delete v2;

    /* Calculate new CP and trainingP: */
    trainPars->trainP = log(varianceVector->getDeterminant()*pow(TWO_PI,varianceVector->len()))*denom;
    calculateCP();
    /* Clean up: */
    trainPars->numeratorVariance->setAllValues(0.0);
    trainPars->numeratorMean->setAllValues(0.0);
    trainPars->trainingDenominator = 0.0;
    trainPars->trainingSamples = 0;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The SMAPLR adaptation starts with creating one single cluster (node) at which all gaussians in
/// the system will need to register. The method adapt_setInitialNode() will set the Gaussian in
/// training mode (enterTrainingMode()) and add its mean vector to the training pool of the
/// adaptation cluster (Adapt_AM_TreeNode).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::adapt_setInitialNode(Adapt_AM_TreeNode *node)
{
  if(trainPars == NULL)
  {
    enterTrainingMode();
  }
  trainPars->adaptNode = node;
  // Using the adaptation data weight as weight for clustering!
  trainPars->adaptNode->meanGaussian->train(trainPars->trainingDenominator,meanVector);
  //trainPars->adaptNode->meanGaussian->train(1.0,meanVector);
  trainPars->adaptNode->denominator += trainPars->trainingDenominator;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The Gaussian has registered at least once to an adaptation node when this method is called.
/// The first time adapt_setInitialNode() shall be used. This method will calculate which one of
/// the two Adapt_AM_TreeNode children has a mean vector closest to the mean vector of the Gaussian
/// itself. It will register to the closest child adaptation node and add its mean vector to the
/// adaptation training pool.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::adapt_setNode()
{
  if(trainPars == NULL)
  {
    enterTrainingMode();
  }
  if(trainPars->adaptNode != NULL && trainPars->adaptNode->childLeft != NULL && trainPars->adaptNode->childRight != NULL)
  {

    double resL = trainPars->adaptNode->childLeft->meanGaussian->getP(meanVector);
    double resR = trainPars->adaptNode->childRight->meanGaussian->getP(meanVector);
    if(resL > resR)
    {
      trainPars->adaptNode = trainPars->adaptNode->childLeft;
    }
    else
    {
      trainPars->adaptNode = trainPars->adaptNode->childRight;
    }
    // Using the adaptation data weight as weight for clustering!
    trainPars->adaptNode->meanGaussian->train(trainPars->trainingDenominator,meanVector);
//    trainPars->adaptNode->meanGaussian->train(1.0,meanVector);
    trainPars->adaptNode->denominator += trainPars->trainingDenominator;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method adds its mean vector statistics to the helper matrix of the adaptation node it is
/// registered to.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::adapt_setHelperMatrices()
{
  int vectorSize = varianceVector->len();
  if(trainPars == NULL)
  {
    enterTrainingMode();
  }
  Vector *numeratorMean     = trainPars->numeratorMean->divideElements(false,varianceVector);;
  Vector *numeratorVariance = new Vector(numeratorMean->len());
  numeratorVariance->setAllValues(trainPars->trainingDenominator);
  numeratorVariance->divideElements(true,varianceVector);
  if(trainPars->adaptNode != NULL)
  {
    Adapt_AM_TreeNode *trainNode = trainPars->adaptNode;
    double *augmentMatrix = trainNode->AUGMENTED_ALL;
    int     offset_k_x    = 0;
    int     offset_k_y    = 0;
    for(int k=0;k<vectorSize;k++)   // In this case: if i != k, then the result is zero!
    {
      double numerator_meanK = numeratorMean->getValue(k);   // Code optimalization
      double mean_meanK      = meanVector->getValue(k);                // Code optimalization
      int offset_j = offset_k_x + (vectorSize+1);
      for(int j=0;j<vectorSize;j++)
      {
        // This is Z(i,j)
        double mean_meanJ = meanVector->getValue(j);                   // Code optimalization
        double meanJKProd = mean_meanK * mean_meanJ;                   // Code optimalization
        augmentMatrix[offset_j] += numerator_meanK * mean_meanJ;
        int offset_i = offset_k_y + j;
        for(int i=0;i<vectorSize;i++)
        {
          augmentMatrix[offset_i] += (numeratorVariance->getValue(i) * meanJKProd);
          offset_i += AUG_X;
        }
        offset_j += AUG_Y;
      }
      // Now also do p+1 (j+1):
      augmentMatrix[offset_k_x + AUG_P_PLUS_ONE] += numerator_meanK;
      int offset_1 = offset_k_y + vectorSize;
      int offset_2 = AUG_Y_TIMES_VECTOR + k;
      int offset_i_x = 0;
      for(int i=0;i<vectorSize;i++)
      {
        double t = (numeratorVariance->getValue(i)*mean_meanK);
        augmentMatrix[offset_i_x + offset_1] += t;
        augmentMatrix[offset_i_x + offset_2] += t;
        offset_i_x += AUG_X;
      }
      augmentMatrix[offset_k_x + AUG_P] += (numeratorVariance->getValue(k));
      offset_k_x += AUG_X;
      offset_k_y += AUG_Y;
    }
  }
  delete numeratorVariance;
  delete numeratorMean;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::adapt_setVarTrans()
{
  assert(trainPars != NULL);
  if(trainPars->adaptNode != NULL)
  {
    trainPars->adaptNode->varDenom += trainPars->trainingDenominator;
    Vector *resV = trainPars->numeratorVariance->copyVector();
    Vector *m    = trainPars->numeratorMean->multiplyElements(false,meanVector);
    m->multiplyElements(true,2.0);
    resV->substractElements(true,m);
    delete m;
    m = meanVector->multiplyElements(false,meanVector);
    m->multiplyElements(true,trainPars->trainingDenominator);
    resV->addElements(true,m);
    delete m;
    resV->divideElements(true,varianceVector);
    trainPars->adaptNode->varTrans->addElements(true,resV);
    delete resV;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The Gaussian object will adapt its mean Vector with use of the adaptation matrix that is
/// calculated by the Adapt_AM_TreeNode adaptation node that the Gaussian is registered to. The
/// old mean Vector is kept in memory and can be restored by calling adapt_unAdapt().
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::adapt_adapt()
{
  if(trainPars == NULL)
  {
    enterTrainingMode();
  }
  if(trainPars->adaptNode != NULL)
  {
    int vectorSize = varianceVector->len();
    Vector newMean(vectorSize);
    if(trainPars->unAdaptedMeanVector != NULL)
    {
      delete trainPars->unAdaptedMeanVector;
      trainPars->unAdaptedMeanVector = NULL;
    }
    // Now adapt this gaussian using the W matrix of the adaptation node we belong to:
    for(int row=0;row<vectorSize;row++)
    {
      double res = 0.0;
      for(int col=0;col<vectorSize;col++)
      {
        res += meanVector->getValue(col)*trainPars->adaptNode->W[row][col];
      }
      newMean.setValue(row,res+trainPars->adaptNode->W[row][vectorSize]);
    }
//    trainPars->unAdaptedMeanVector = meanVector;
    delete meanVector;
    meanVector = newMean.copyVector();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::adapt_adaptVar()
{
  Vector *trans = trainPars->adaptNode->varTrans->multiplyElements(false,1.0 / trainPars->adaptNode->varDenom);
  varianceVector->multiplyElements(true,trans);
  delete trans;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The Gaussian object will adapt its mean Vector with use of the adaptation matrix that is
/// calculated by the Adapt_AM_TreeNode adaptation node that the Gaussian is registered to when the
/// method adapt_adapt() is called. The method adapt_unAdapt() will restore the old
/// mean Vector that is kept in memory.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::adapt_unAdapt()
{
  assert(false); // Unadapt is switched of for this release because of memory issues.
  if(trainPars == NULL)
  {
    enterTrainingMode();
  }
  if(trainPars->unAdaptedMeanVector != NULL)
  {
    delete meanVector;
    meanVector = trainPars->unAdaptedMeanVector;
    trainPars->unAdaptedMeanVector = NULL;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will clear all adaptation settings, so that a new adaptation iteration can be started.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::adapt_clear()
{
  if(trainPars == NULL)
  {
    enterTrainingMode();
  }
  trainPars->numeratorVariance->setAllValues(0.0);
  trainPars->numeratorMean->setAllValues(0.0);
  trainPars->trainingDenominator = 0.0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will write accumulators to file
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::writeAccumulators(FILE *outFile, FILE *outFileF, FILE *outFileST, bool doBinary)
{
  if(trainPars == NULL)
  {
    enterTrainingMode();
  }
  if(outFileF == NULL)
  {
    trainPars->numeratorMean->storeData(outFile);
    trainPars->numeratorVariance->storeData(outFile);
    fwrite(&(trainPars->trainingDenominator),1,sizeof(trainPars->trainingDenominator),outFile);
  }
  else
  {
    if(doBinary)
    {
      fwrite(&(trainPars->trainingDenominator),1,sizeof(trainPars->trainingDenominator),outFile);
      for(int i=0;i<trainPars->numeratorMean->len();i++)
      {
        double m = trainPars->numeratorMean->getValue(i);
        fwrite(&m,1,sizeof(m), outFileF);
      }
      if(outFileST != NULL)
      {
        for(int i=0;i<trainPars->numeratorVariance->len();i++)
        {
          double m = trainPars->numeratorVariance->getValue(i);
          fwrite(&m,1,sizeof(m), outFileST);
        }
      }
    }
    else
    {
      fprintf(outFile,"%g ",trainPars->trainingDenominator);
      for(int i=0;i<trainPars->numeratorMean->len();i++)
      {
        fprintf(outFileF,"%g ",trainPars->numeratorMean->getValue(i));
      }
      if(outFileST != NULL)
      {
        for(int i=0;i<trainPars->numeratorVariance->len();i++)
        {
          fprintf(outFileST,"%g ",trainPars->numeratorVariance->getValue(i));
        }
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will read accumulators from file
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::addAccumulators(FILE *inFile)
{
  if(trainPars == NULL)
  {
    enterTrainingMode();
  }
  double mmiDenom = 0.0;
  Vector *mmiMean = new Vector(inFile,trainPars->numeratorMean->len(),false);
  Vector *mmiVar  = new Vector(inFile,trainPars->numeratorMean->len(),false);
  freadEndianSafe(&mmiDenom,1,sizeof(mmiDenom),inFile);
  trainPars->numeratorMean->addElements(true,mmiMean);
  trainPars->numeratorVariance->addElements(true,mmiVar);
  trainPars->trainingDenominator += mmiDenom;
  delete mmiMean;
  delete mmiVar;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will read accumulators from file and perform MMI
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::trainMMI(FILE *fileEnum, FILE *fileDenom)
{
  if(trainPars == NULL)
  {
    enterTrainingMode();
  }
  if(mmi_trainPars == NULL)
  {
    enterTrainingMode(true);
  }
  int len = trainPars->numeratorMean->len();
  delete trainPars->numeratorMean;
  delete trainPars->numeratorVariance;
  delete mmi_trainPars->numeratorMean;
  delete mmi_trainPars->numeratorVariance;

  trainPars->trainingDenominator = 0.0;
  trainPars->numeratorMean       = new Vector(fileEnum,len,false);
  trainPars->numeratorVariance   = new Vector(fileEnum,len,false);
  freadEndianSafe(&(trainPars->trainingDenominator),1,sizeof(trainPars->trainingDenominator),fileEnum);

  mmi_trainPars->trainingDenominator = 0.0;
  mmi_trainPars->numeratorMean       = new Vector(fileDenom,len,false);
  mmi_trainPars->numeratorVariance   = new Vector(fileDenom,len,false);
  freadEndianSafe(&(mmi_trainPars->trainingDenominator),1,sizeof(mmi_trainPars->trainingDenominator),fileDenom);

  if((trainPars->trainingDenominator > 0.0) && (mmi_trainPars->trainingDenominator > 0.0))
  {
    double D = 2.0*mmi_trainPars->trainingDenominator;
    bool tooSmall      = true;
    double denom    = 1.0 / trainPars->trainingDenominator;
    Vector *meanNew = NULL;
    Vector *varNew  = NULL;
    for(int iii=0;iii<2;iii++)
    {
      int nrIt = 0;
      tooSmall = true;
      while(tooSmall && (nrIt<100))
      {
        nrIt++;
        PRINTF2("D = %f... ",D);fflush(stdout);
        denom  = 1.0 / (trainPars->trainingDenominator - mmi_trainPars->trainingDenominator + D);
        Vector *nMean = meanVector->multiplyElements(false,D);
        nMean->addElements(true, trainPars->numeratorMean);
        nMean->substractElements(true,mmi_trainPars->numeratorMean);
        Vector *nVar  = meanVector->multiplyElements(false,meanVector);
        nVar->addElements(true,varianceVector);
        nVar->multiplyElements(true,D);
        nVar->addElements(true,trainPars->numeratorVariance);
        nVar->substractElements(true,mmi_trainPars->numeratorVariance);

        if(meanNew != NULL)
        {
          delete meanNew;
          delete varNew;
        }
        meanNew = nMean->multiplyElements(false,denom);
        /* Now the new variance: */
        Vector *v = meanNew->multiplyElements(false,meanNew);
        varNew = nVar->multiplyElements(false,denom);
        varNew->substractElements(true,v);
        tooSmall = varNew->minimizeElements(GAUSSIAN_VARIANCE_MINIMUM);
        delete nMean;
        delete nVar;
        delete v;
        D *= 2.0;
        PRINTF2("%d\n",tooSmall);fflush(stdout);
      }
    }
    assert(meanNew != NULL);
    assert(varNew != NULL);
    if(!tooSmall)
    {
      varNew->minimizeElements(GAUSSIAN_VARIANCE_MINIMUM);
      delete meanVector;
      delete varianceVector;
      meanVector     = meanNew;
      varianceVector = varNew;
    }
    else
    {
      delete meanNew;
      delete varNew;
    }
    /* Calculate new CP and trainingP: */
    trainPars->trainP = log(varianceVector->getDeterminant()*pow(TWO_PI,varianceVector->len())) / denom;
    calculateCP();
    /* Clean up: */
    trainPars->numeratorVariance->setAllValues(0.0);
    trainPars->numeratorMean->setAllValues(0.0);
    trainPars->trainingDenominator = 0.0;
    trainPars->trainingSamples = 0;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Prints the gaussian to the standard output for debugging purposes.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Gaussian::printGaussian(void)
{
#ifdef DEBUG_ON
  printf("Means:\n");
  meanVector->printVector();
  printf("Variance:\n");
  varianceVector->printVector();
#endif
}


