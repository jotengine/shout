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

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "standard.h"
#include "mixgaussian.h"
#include "shout-misc.h"

using namespace WriteFileLittleBigEndian;


int     maxGaussianTrainThreads = -1;



int MixGaussian::printInfo(Vector *v)
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->printInfo(v);
  }
  return numberOfGaussians;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::printModel(FILE *fileMean, FILE *fileVariance, FILE *fileWeight)
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->printModel(fileMean,fileVariance);
    fprintf(fileWeight,"%g ",myGaussianWeights[i]);
  }
  fprintf(fileWeight,"\n");
  fprintf(fileMean,"\n");
  fprintf(fileVariance,"\n");
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The constructor initialises the object. It creates a single gaussian (mixture of one) and makes
/// a weight table of size one with weight one.
/////////////////////////////////////////////////////////////////////////////////////////////////////

MixGaussian::MixGaussian(int dim)
{
  fastP_CPandWeight_meanAndVar = NULL;
  fastP_res                    = NULL;

  myGaussiansCounts       = NULL;
  myGaussiansSorted       = NULL;
  numberOfGaussians       = 1;
  myGaussians             = new Gaussian*[numberOfGaussians];
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i] = new Gaussian(dim);
  }
  myGaussianWeights       = new double[numberOfGaussians];
  myGaussianWeights[0]    = 1.0;
  noSplitting             = false;
  vectorDimension         = dim;
  initializeFastCalc(false);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor will merge the two input gaussians.
/////////////////////////////////////////////////////////////////////////////////////////////////////

MixGaussian::MixGaussian(MixGaussian *model1, MixGaussian *model2)
{
  fastP_CPandWeight_meanAndVar = NULL;
  fastP_res                    = NULL;

  myGaussiansCounts = NULL;
  myGaussiansSorted = NULL;
  myGaussians       = NULL;
  myGaussianWeights = NULL;
  noSplitting       = false;  
  mergeMixGaussians(model1,model2);
  vectorDimension   = myGaussians[0]->dim();
  initializeFastCalc(false);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor will merge the two input gaussians.
/////////////////////////////////////////////////////////////////////////////////////////////////////

MixGaussian::MixGaussian(MixGaussian *model1, MixGaussian *model2, double rate)
{ 
  fastP_CPandWeight_meanAndVar = NULL;
  fastP_res                    = NULL;

  assert(model1 != NULL);
  assert(model2 != NULL);
  assert(model1->dim() == model2->dim());
  assert(model1->numberOfGaussians == model2->numberOfGaussians);
  myGaussiansCounts = NULL;
  myGaussiansSorted = NULL;
  myGaussians       = NULL;
  myGaussianWeights = NULL;
  
  numberOfGaussians = model1->numberOfGaussians;
  myGaussians       = new Gaussian*[numberOfGaussians];
  myGaussianWeights = new double[numberOfGaussians];
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]       = new Gaussian(model1->myGaussians[i],model2->myGaussians[i],rate);
    myGaussianWeights[i] = model1->myGaussianWeights[i];
  }
  vectorDimension  = myGaussians[0]->dim();
  noSplitting      = false;
  initializeFastCalc(false);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::initializeFastCalc(bool onlyOnce)
{
  assert(myGaussians != NULL);
  lookaheadPointer = GMMLOOKAHEAD*5;
  lastObservation  = NULL;
  if(fastP_res != NULL)
  {
    delete[] fastP_res;
  }
  if(fastP_CPandWeight_meanAndVar != NULL)
  {
    delete[] fastP_CPandWeight_meanAndVar;
  }


  int len           = dim();
  fastP_res         = new double[numberOfGaussians*HALF_GMMLOOKAHEAD];

#ifdef CELL_PROCESSOR
  if(posix_memalign((void**)(&fastP_CPandWeight_meanAndVar), 128, (numberOfGaussians*2*len + numberOfGaussians) * sizeof(float)) != 0)
  {
    assert(false);
  }
#else

  fastP_CPandWeight_meanAndVar  = new double[numberOfGaussians*2*len + numberOfGaussians];
#endif

  int cnt = 0;
  double *fastP_CPandWeight = fastP_CPandWeight_meanAndVar;
  double *fastP_meanAndVar  = &fastP_CPandWeight_meanAndVar[numberOfGaussians];
  for(int i2=0;i2<len;i2++)
  {
    for(int i=0;i<numberOfGaussians;i++)
    {
      fastP_CPandWeight[i]    = myGaussians[i]->getCP() + log(myGaussianWeights[i]);
      fastP_meanAndVar[cnt]   = myGaussians[i]->getMean()->getValue(i2);
      fastP_meanAndVar[cnt+1] = -0.5 / myGaussians[i]->getVariance()->getValue(i2);
      cnt += 2;
    }
  }
  if(onlyOnce)
  {
    if(myGaussians != NULL)
    {
      for(int i=0;i<numberOfGaussians;i++)
      {  
        delete myGaussians[i];
      }
      delete[] myGaussians;
      myGaussians = NULL;
    }
    if(myGaussianWeights != NULL)
    {
      delete[] myGaussianWeights;  
      myGaussianWeights = NULL;
    }
    if(myGaussiansCounts != NULL)
    {
      delete[] myGaussiansCounts;
      myGaussiansCounts = NULL;
    }
    if(myGaussiansSorted != NULL)
    {
      delete[] myGaussiansSorted;
      myGaussiansSorted = NULL;
    }
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Copy the settings from the org MixGaussian.
/////////////////////////////////////////////////////////////////////////////////////////////////////

MixGaussian::MixGaussian(MixGaussian *org)
{
  fastP_CPandWeight_meanAndVar  = NULL;
  fastP_res                     = NULL;

  myGaussiansCounts       = NULL;
  myGaussiansSorted       = NULL;
  numberOfGaussians       = org->numberOfGaussians;
  myGaussians             = new Gaussian*[numberOfGaussians];
  myGaussianWeights       = new   double[numberOfGaussians];
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i] = new Gaussian(org->myGaussians[i],false);
    myGaussianWeights[i] = org->myGaussianWeights[i];
  }
  noSplitting             = org->noSplitting;
  vectorDimension         = org->myGaussians[0]->dim();
  initializeFastCalc(false);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor reads all its data from file. (number of gaussians, the weight table and the 
/// gaussian parameters)
/////////////////////////////////////////////////////////////////////////////////////////////////////

MixGaussian::MixGaussian(FILE *inFile, int dim, bool noReTraining)
{
  fastP_CPandWeight_meanAndVar = NULL;
  fastP_res                    = NULL;

  myGaussiansCounts       = NULL;
  myGaussiansSorted       = NULL;
  noSplitting             = false;

  // Read number of gaussians:
  freadEndianSafe(&numberOfGaussians,1,sizeof(numberOfGaussians),inFile);
  // Initialize with new number of gaussians:
  myGaussians             = new Gaussian*[numberOfGaussians];
  myGaussianWeights       = new double[numberOfGaussians];

  // Read the gaussian values:  
  for(int i=0;i<numberOfGaussians;i++)
  {
    freadEndianSafe(&(myGaussianWeights[i]),1,sizeof(double),inFile);
    myGaussians[i] = new Gaussian(inFile,dim);
  }
  vectorDimension = dim;
  initializeFastCalc(noReTraining);
}
  
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor deletes all gaussians from the mix and also deletes the weights table.
/////////////////////////////////////////////////////////////////////////////////////////////////////

MixGaussian::~MixGaussian()
{
  if(fastP_CPandWeight_meanAndVar != NULL)
  {
    delete[] fastP_res;
#ifdef CELL_PROCESSOR
    free(fastP_CPandWeight_meanAndVar);
#else
    delete[] fastP_CPandWeight_meanAndVar;
#endif
  }

  if(myGaussians != NULL)
  {
    for(int i=0;i<numberOfGaussians;i++)
    {  
      delete myGaussians[i];
    }
    delete[] myGaussians;
  }
  if(myGaussianWeights != NULL)
  {
    delete[] myGaussianWeights;  
  }
  if(myGaussiansCounts != NULL)
  {
    delete[] myGaussiansCounts;
  }
  if(myGaussiansSorted != NULL)
  {
    delete[] myGaussiansSorted;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::getSumHist(Vector **histogram)
{
  int vectorLength = histogram[0]->len();
  Vector minVector(vectorLength);
  Vector maxVector(vectorLength);
  Vector stepVector(vectorLength);
  Vector v(vectorLength);
  Vector minv(vectorLength);
  Vector maxv(vectorLength);
  for (int i=0; i<numberOfGaussians;i++)
  {
    v.setAllValues(myGaussians[i]->getVariance());
    v.sqrtElements(true);
    v.multiplyElements(true,(double)(FEATURE_NORM_TIMESVAR));
    minv.setAllValues(myGaussians[i]->getMean());
    maxv.setAllValues(myGaussians[i]->getMean());
    minv.substractElements(true,&v);
    maxv.addElements(true,&v);
    if(i == 0)
    {
      minVector.setAllValues(&minv);
      maxVector.setAllValues(&maxv);
    }
    else
    {
      for(int i2=0;i2<vectorLength;i2++)
      {
        if(minVector.getValue(i2) > minv.getValue(i2))
        {
          minVector.setValue(i2,minv.getValue(i2));
        }
        if(maxVector.getValue(i2) < maxv.getValue(i2))
        {
          maxVector.setValue(i2,maxv.getValue(i2));
        }
      }
    }
  }  
  stepVector.setAllValues(&maxVector);
  stepVector.substractElements(true,&minVector);
  stepVector.multiplyElements(true,1.0 / ((double)FEATURE_NORM_HISTOGRAM_BINS-3)); // (-1 and -2 are reserved for step size and min)

  histogram[FEATURE_NORM_HISTOGRAM_BINS-1]->setAllValues(&stepVector);
  histogram[FEATURE_NORM_HISTOGRAM_BINS-2]->setAllValues(&minVector);
  histogram[0]->setAllValues(0.0);

  Gaussian **my1DGaussians = new Gaussian*[numberOfGaussians];
  for(int i=0;i<numberOfGaussians;i++)
  {
    my1DGaussians[i] = new Gaussian(1);
  }
  Vector var1D(1);
  Vector mean1D(1);
  Vector observation(1);

  for(int i2=0;i2<vectorLength;i2++)
  {
    for(int i=0;i<numberOfGaussians;i++)
    {
      mean1D.setValue(0,myGaussians[i]->getMean()->getValue(i2));
      var1D.setValue(0,myGaussians[i]->getVariance()->getValue(i2));
      my1DGaussians[i]->setMean(&mean1D);
      my1DGaussians[i]->setVariance(&var1D);
    }    
    observation.setValue(0,minVector.getValue(i2));
    for(int i=0;i<FEATURE_NORM_HISTOGRAM_BINS-2;i++)
    {
      if(i != 0)
      {
        histogram[i]->setValue(i2,histogram[i-1]->getValue(i2));
      }
      double res = 0.0; 
      for(int i3=0;i3<numberOfGaussians;i3++)
      {
        res += myGaussianWeights[i3] * my1DGaussians[i3]->getP(&observation);
      }
      histogram[i]->addValue(i2,res);
      observation.addValue(0,stepVector.getValue(i2));
    }
  }
  
  for(int i2=0;i2<numberOfGaussians;i2++)
  {
    delete my1DGaussians[i2];
  }
  delete[] my1DGaussians;
  for(int i=0;i<FEATURE_NORM_HISTOGRAM_BINS-2;i++)
  {
    histogram[i]->divideElements(true,histogram[FEATURE_NORM_HISTOGRAM_BINS-3]);
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Set the variance for all of the gaussians in this mixture
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::setVariance(const Vector *v)
{
  assert(v->len() == dim());
  for (int i=0; i<numberOfGaussians;i++)
    myGaussians[i]->setVariance(v);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Merge the two input gaussians.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::mergeMixGaussians(MixGaussian *model1, MixGaussian *model2)
{
  assert(model1 != NULL);
  assert(model2 != NULL);
  assert(model1->dim() == model2->dim());

  if(myGaussians != NULL)
  {
    for(int i=0;i<numberOfGaussians;i++)
    {  
      delete myGaussians[i];
    }
    delete[] myGaussians;
  }
  if(myGaussianWeights != NULL)
  {
    delete[] myGaussianWeights;  
  }
  
  numberOfGaussians = model1->numberOfGaussians + model2->numberOfGaussians;
  myGaussians       = new Gaussian*[numberOfGaussians];
  myGaussianWeights = new double[numberOfGaussians];
  int model1G       = model1->numberOfGaussians;
  double weight     = ((double)(model1G)/((double)numberOfGaussians));
  for(int i=0;i<model1G;i++)
  {
    myGaussians[i]       = new Gaussian(model1->myGaussians[i],false);
    myGaussianWeights[i] = model1->myGaussianWeights[i]*weight;
  }
  weight = 1.0 - weight;
  for(int i=0;i<model2->numberOfGaussians;i++)
  {
    myGaussians[i+model1G]       = new Gaussian(model2->myGaussians[i],false);
    myGaussianWeights[i+model1G] = model2->myGaussianWeights[i]*weight;
  }
  noSplitting      = false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the number of gaussians in the mix.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int MixGaussian::getNumberOfGaussians(void)
{
  return numberOfGaussians;
}

#ifdef SPLIT_CLUSTERING    
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Fills the distance array: an array of size numberOfGaussians. The two
/// most distant gaussians are in the first and last spot. The others are sorted in between.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::fillDistanceArray(int *distA)
{
  double biggestDist = -1.0;
  int    biggestDi1  = -1;
  int    biggestDi2  = -1;
  for(int i=0;i<numberOfGaussians-1;i++)
  {
    if(myGaussianWeights[i] > 0.0)
    {
      for(int i2=i+1;i2<numberOfGaussians;i2++)
      {
        if(myGaussianWeights[i2] > 0.0)
        {
          double dist = myGaussians[i]->meanVector->distance(myGaussians[2]->meanVector);
          if(dist > biggestDist)
          {
            biggestDist = dist;
            biggestDi1  = i;
            biggestDi2  = i2;
          }
        }
      }
    }
  }
  double distances[numberOfGaussians];
  for(int i=0;i<numberOfGaussians;i++)
  {
    distances[i] = -1.0;
    distA[i]     = -1;
  }  
  distA[0]                   = biggestDi1;
  distA[numberOfGaussians-1] = biggestDi2;
  Vector *v = myGaussians[biggestDi1]->meanVector;
  for(int i=0;i<numberOfGaussians;i++)
  {
    if((i != biggestDi1) && (i != biggestDi2) && (myGaussianWeights[i] > 0.0))
    {
      distances[i] = v->distance(myGaussians[i]->meanVector);
    }  
  }
  for(int i=1;i<numberOfGaussians-1;i++)
  {
    double smallest = biggestDist + 1.0;
    for(int i2=0;i2<numberOfGaussians;i2++)
    {
      if(distances[i2] >= 0.0 && distances[i2] < smallest)
      {
        smallest = distances[i2];
        distA[i] = i2;
      }
    }
    distances[distA[i]] = -1.0;
  }    
}
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// GetLogP() returns the chance of an observation beeing produced by this mixture of gaussians (in log)
/////////////////////////////////////////////////////////////////////////////////////////////////////


double MixGaussian::getLookaheadLogP(double *vectorList, int timeStamp, bool doSecondHalf)
{
  int pointer = timeStamp - lookaheadPointer;

  if(doSecondHalf && (pointer >= GMMLOOKAHEAD))
  {
    return -9.0e30;
  }

  if(pointer == 0 && doSecondHalf)
  {
    timeStamp = HALF_GMMLOOKAHEAD;
  }
  else if((pointer >= 0) && (pointer < GMMLOOKAHEAD))
  {
    return lookaheadRes[pointer];
  }
  else
  {
    lookaheadPointer = timeStamp;
    timeStamp = 0;
  }

  assert(fastP_CPandWeight_meanAndVar != NULL);
  assert(fastP_res != NULL);

  double *runThroughFastP = &fastP_CPandWeight_meanAndVar[numberOfGaussians];
  double tmp[HALF_GMMLOOKAHEAD];
  double threshold[HALF_GMMLOOKAHEAD];
  double *runThroughFastP_res[HALF_GMMLOOKAHEAD];
  for(int i=0;i<HALF_GMMLOOKAHEAD;i++)
  {
    runThroughFastP_res[i] = &fastP_res[i*numberOfGaussians];
    memcpy(runThroughFastP_res[i],fastP_CPandWeight_meanAndVar,sizeof(double)*numberOfGaussians);
    lookaheadRes[i+timeStamp] = 0.0;
    threshold[i] = -9.0e100;
  }
  int len = dim();
  for(int i2=0;i2<len-1;i2++)
  {
    for(int i=0;i<numberOfGaussians;i++)
    {
      double *vl = vectorList;
      for(int la=0;la<HALF_GMMLOOKAHEAD;la++)
      {
        tmp[la] = pow(*vl-*runThroughFastP,2);
        vl++;
      }
      runThroughFastP++;
      for(int la=0;la<HALF_GMMLOOKAHEAD;la++)
      {
        *runThroughFastP_res[la] += tmp[la]*(*runThroughFastP);
         runThroughFastP_res[la]++;
      }
      runThroughFastP++;
    }
    vectorList+=HALF_GMMLOOKAHEAD;
    for(int i=0;i<HALF_GMMLOOKAHEAD;i++)
    {
      runThroughFastP_res[i] = &fastP_res[i*numberOfGaussians];
    }
  }
  for(int i=0;i<numberOfGaussians;i++)
  {
    double *vl = vectorList;
    for(int la=0;la<HALF_GMMLOOKAHEAD;la++)
    {
      tmp[la] = pow(*vl-*runThroughFastP,2);
      vl++;
    }

    runThroughFastP++;
    for(int la=0;la<HALF_GMMLOOKAHEAD;la++)
    {
      *runThroughFastP_res[la] += tmp[la]*(*runThroughFastP);
      if(*runThroughFastP_res[la] > threshold[la])
      {
        threshold[la] = *runThroughFastP_res[la];
      }
      runThroughFastP_res[la]++;
    }
    runThroughFastP++;
  }
  for(int i=0;i<HALF_GMMLOOKAHEAD;i++)
  {
    runThroughFastP_res[i] = &fastP_res[i*numberOfGaussians];
    threshold[i] -= 7.0;
  }

  for(int i=0;i<numberOfGaussians;i++)
  {
    for(int la=0;la<HALF_GMMLOOKAHEAD;la++)
    {
      if((*runThroughFastP_res[la] > threshold[la]))
      {
        *runThroughFastP_res[la] = exp(*runThroughFastP_res[la]);
        lookaheadRes[la+timeStamp] += *runThroughFastP_res[la];
      }
      else
      {
        *runThroughFastP_res[la] = 0.0;
      }
      runThroughFastP_res[la]++;
    }
  }
  for(int i=0;i<HALF_GMMLOOKAHEAD;i++)
  {
    lookaheadRes[i+timeStamp] = log(lookaheadRes[i+timeStamp]);
    if(isnan(lookaheadRes[i+timeStamp]) || lookaheadRes[i+timeStamp] < MINIMUM_LOGP)
    {
      lookaheadRes[i+timeStamp] = MINIMUM_LOGP;
    }
  }
  return lookaheadRes[0];
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// GetLogP() returns the chance of an observation beeing produced by this mixture of gaussians (in log)
/////////////////////////////////////////////////////////////////////////////////////////////////////

double MixGaussian::getLogP(Vector *observation)
{
  lastObservation = observation;
  assert(fastP_CPandWeight_meanAndVar != NULL);
  assert(fastP_res != NULL);

  memcpy(fastP_res,fastP_CPandWeight_meanAndVar,sizeof(double)*numberOfGaussians);

  double *runThroughFastP = &fastP_CPandWeight_meanAndVar[numberOfGaussians];
  double tmp;
  double threshold = -9.0e100;
  int len = dim();
  double *vArray = observation->getVectorArray();
  for(int i2=0;i2<len-1;i2++)
  {
    for(int i=0;i<numberOfGaussians;i++)
    {
      tmp = pow(*vArray-*runThroughFastP,2);
      runThroughFastP++;
      fastP_res[i] += tmp*(*runThroughFastP);
      runThroughFastP++;
    }
    vArray++;
  }
  double res = 0.0;
  for(int i=0;i<numberOfGaussians;i++)
  {
    tmp = pow(*vArray-*runThroughFastP,2);
    runThroughFastP++;
    fastP_res[i] += tmp*(*runThroughFastP);
    if(fastP_res[i] > threshold)
    {
      threshold = fastP_res[i];
    }
    runThroughFastP++;
  }

  threshold -= 7.0;
  for(int i=0;i<numberOfGaussians;i++)
  {
    if((fastP_res[i] > threshold))
    {
      fastP_res[i] = exp(fastP_res[i]);
      res += fastP_res[i];
    }
    else
    {
      fastP_res[i] = 0.0;
    }
  }
  res = log(res);
  if(isnan(res) || res < MINIMUM_LOGP)
  {
    res = MINIMUM_LOGP;
  }
  return res;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// GetP() returns the chance of an observation beeing produced by this mixture of gaussians.
/////////////////////////////////////////////////////////////////////////////////////////////////////

double MixGaussian::getP(Vector *observation)
{
  lastObservation = observation->getKey();
  assert(fastP_CPandWeight_meanAndVar != NULL);
  assert(fastP_res != NULL);

  memcpy(fastP_res,fastP_CPandWeight_meanAndVar,sizeof(double)*numberOfGaussians);

  double *runThroughFastP = &fastP_CPandWeight_meanAndVar[numberOfGaussians];
  double tmp,v;
  double threshold = -9.0e100;
  int len = dim();
  double *vArray = observation->getVectorArray();
  for(int i2=0;i2<len-1;i2++)
  {
    for(int i=0;i<numberOfGaussians;i++)
    {
      tmp = pow(*vArray-*runThroughFastP,2);
      runThroughFastP++;
      fastP_res[i] += tmp*(*runThroughFastP);
      runThroughFastP++;
    }
    vArray++;
  }
  double res = 0.0;
  for(int i=0;i<numberOfGaussians;i++)
  {
    tmp = pow(*vArray-*runThroughFastP,2);
    runThroughFastP++;
    fastP_res[i] += tmp*(*runThroughFastP);
    if(fastP_res[i] > threshold)
    {
      threshold = fastP_res[i];
    }
    runThroughFastP++;
  }

  threshold -= 7.0;
  for(int i=0;i<numberOfGaussians;i++)
  {
    if((fastP_res[i] > threshold))
    {
      fastP_res[i] = exp(fastP_res[i]);
      res += fastP_res[i];
    }
    else
    {
      fastP_res[i] = 0.0;
    }
  }
  if(isnan(res) || res < SMALLEST_NEGATIVE)
  {
    res = SMALLEST_NEGATIVE;
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method searches the "best" gaussian, and shifts its mean vector.
///
/// The best gaussian is the gaussian with the largest weight. The weight is chosen, because indirectly,
/// the weight represents the amount of training samples defining the gaussian. The mean of the new
/// gaussian is 0.2*shiftFactor times shifted (shiftFactor may be negative).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::shiftBestGaussian(int shiftFactor)
{
  double biggestW    = 0.0;
  int winner         = -1;
  for(int i=0;i<numberOfGaussians;i++)
  {
    if(biggestW < myGaussianWeights[i])
    {
      biggestW =  myGaussianWeights[i];
      winner   = i;
    }
  } 
  if(winner >= 0)
  {
    myGaussians[winner]->shiftMean(shiftFactor);  
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method searches the "best" gaussian, and splits it into two new gaussians as follows:
///
/// The best gaussian is the gaussian with the largest weight. The weight is chosen, because indirectly,
/// the weight represents the amount of training samples defining the gaussian. The mean of the new
/// gaussian is 0.2 times the variance bigger than the old mean. The mean of the existing gaussian is 
/// made 0.2 times the variance smaller. The variance is not changed. Note that the mean is shifted
/// in all dimensions. One could argue if this is a good initialisation step. The weights of the two
/// new gaussians is set half the weight of the original gaussian. Also, it could be argued if this
/// is good practice.
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool MixGaussian::splitBestGaussian(bool alwaysSplit)
{
  if(noSplitting && !alwaysSplit)
  {
    printf("No splitting for this PDF... (%d,%d)\n",noSplitting,alwaysSplit);
    return false;
  }
  int i;
  double biggestW    = 0.0;
  int winner         = -1;
  for(i=0;i<numberOfGaussians;i++)
  {
    if(biggestW < myGaussianWeights[i])
    {
      biggestW =  myGaussianWeights[i];
      winner   = i;
    }
  }

  double    *tw      = new    double[numberOfGaussians+1];
  Gaussian **tg      = new Gaussian*[numberOfGaussians+1];
  for(i=0;i<numberOfGaussians;i++)
  {
    tw[i] = myGaussianWeights[i];
    tg[i] = myGaussians[i];
  }
  if(winner >=0 && (myGaussians[winner]->getTrainingDenominator() > MIN_TRAIN_SAMPLES_PER_GAUSSIAN*3.0 || alwaysSplit))
  {
    tg[numberOfGaussians] = new Gaussian(myGaussians[winner]);
    tg[winner]->shiftMean();
    tw[numberOfGaussians]     = myGaussianWeights[winner] / 2;
    myGaussianWeights[winner] = tw[numberOfGaussians];
    numberOfGaussians++;

    delete[] myGaussianWeights;
    myGaussianWeights = tw;
    delete[] myGaussians;
    myGaussians = tg;  
    initializeFastCalc(false);
    return true;
  }
  else
  {
    printf("       Winner: %d\n",winner);
    if(winner>=0)
    {
      printf("       Denom: %f\n",myGaussians[winner]->getTrainingDenominator());
    }
    PRINTF("   Splitting stopped during splitting..\n");
    noSplitting = true;
  }
  return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool MixGaussian::splitAllGaussians(bool alwaysSplit)
{
  if(noSplitting && !alwaysSplit)
  {
    printf("No splitting for this PDF... (%d,%d)\n",noSplitting,alwaysSplit);
    return false;
  }

  double    *tw      = new    double[numberOfGaussians*2];
  Gaussian **tg      = new Gaussian*[numberOfGaussians*2];
  for(int i=0;i<numberOfGaussians;i++)
  {
    tw[i] = myGaussianWeights[i];
    tg[i] = myGaussians[i];
  }
  for(int i=0;i<numberOfGaussians;i++)
  {
    tg[i+numberOfGaussians] = new Gaussian(myGaussians[i]);
    tg[i+numberOfGaussians]->shiftMean();
    tw[i+numberOfGaussians] = myGaussianWeights[i] / 2;
    myGaussianWeights[i]    = tw[i+numberOfGaussians];
  }
  numberOfGaussians*=2;

  delete[] myGaussianWeights;
  myGaussianWeights = tw;
  delete[] myGaussians;
  myGaussians = tg;
  initializeFastCalc(false);
  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::addGaussian(Vector *v)
{
  int i;
  double    *tw        = new    double[numberOfGaussians+1];
  Gaussian **tg        = new Gaussian*[numberOfGaussians+1];  
  double weight        = 1.0/(numberOfGaussians+1);
  for(i=0;i<numberOfGaussians;i++)
  {
    tw[i] = weight;
    tg[i] = myGaussians[i];
  }
  tw[numberOfGaussians] = weight;
  tg[numberOfGaussians] = new Gaussian(myGaussians[0],false);
  tg[numberOfGaussians]->setMean(v);
  if(myGaussianWeights != NULL)
  {
    delete[] myGaussianWeights;
  }
  myGaussianWeights = tw;
  if(myGaussians != NULL)
  {    
    delete[] myGaussians;
  }
  myGaussians = tg;  
  numberOfGaussians++;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This is a helper method for copyGaussians(). It adds a gaussian directly to the mix. 
/// The gaussian itself and the weight are provided as input parameters.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::addGaussian(Gaussian *gaussian, double weight)
{
  int i;
  double    *tw        = new    double[numberOfGaussians+1];
  Gaussian **tg        = new Gaussian*[numberOfGaussians+1];  
  for(i=0;i<numberOfGaussians;i++)
  {
    tw[i] = myGaussianWeights[i];
    tg[i] = myGaussians[i];
  }
  tw[numberOfGaussians] = weight;
  tg[numberOfGaussians] = new Gaussian(gaussian,false);
  if(myGaussianWeights != NULL)
  {
    delete[] myGaussianWeights;
  }
  myGaussianWeights = tw;
  if(myGaussians != NULL)
  {    
    delete[] myGaussians;
  }
  myGaussians = tg;  
  numberOfGaussians++;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// With this method, the sum of the weights from the gaussian mix is forced to be one and if
/// the mix contains more than maxGaussians, the gaussians with the smallest weight are pruned.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::normalizeWeights(int maxGaussians)
{
  // First remove the gaussians:
  if(maxGaussians > 0 && maxGaussians < numberOfGaussians)
  {
    double    *tw = new    double[maxGaussians];
    Gaussian **tg = new Gaussian*[maxGaussians];  
    int curG = 0;
    while(curG < maxGaussians)
    {
      double maxW = 0.0;
      int bestI = -1;
      for(int i=0;i<numberOfGaussians;i++)
      {
        if(myGaussianWeights[i] > maxW)
        {
          bestI = i;
          maxW  = myGaussianWeights[i];
        }
      }
      tw[curG] = maxW;
      tg[curG] = myGaussians[bestI];
      myGaussianWeights[bestI] = -900.0;
      curG++;
    }
    for(int i=0;i<numberOfGaussians;i++)
    {
      if(myGaussianWeights[i] != -900.0)
      {
        delete myGaussians[i];
      }
    }    
    delete[] myGaussianWeights;
    myGaussianWeights = tw;
    delete[] myGaussians;
    myGaussians = tg;  
    numberOfGaussians = maxGaussians;    
  }

  double totWeight = 0.0;
  for(int i=0;i<numberOfGaussians;i++)
  {
    totWeight += myGaussianWeights[i];
  }
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussianWeights[i] = myGaussianWeights[i]/totWeight;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method copies the maxNmbr gaussians with the biggest weights into the destination gaussian mix. 
/// This method is used by the SpeechDetection class to make speech/non-speech models.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::copyGaussians(MixGaussian *destMixGaussian, int maxNmbr)
{
  bool beenThere[numberOfGaussians];
  if((maxNmbr < 0) || (maxNmbr > numberOfGaussians))
  {
    maxNmbr = numberOfGaussians;
  }
  
  for(int j=0;j<numberOfGaussians;j++)
  {
    beenThere[j] = false;
  }
  for(int i=0;i<maxNmbr;i++)
  {
    double maxWeight = 0;
    int    luckyG    = -1;
    for(int j=0;j<numberOfGaussians;j++)
    {
      if(!beenThere[j] && myGaussianWeights[j] > maxWeight)
      {
        maxWeight = myGaussianWeights[j];
        luckyG    = j;
      }
    }
    if(luckyG != -1)
    {
      beenThere[luckyG] = true;
      destMixGaussian->addGaussian(myGaussians[luckyG],myGaussianWeights[luckyG]);
    }
    else
    {
      /* Done, there are no more gaussians! */
      return;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method loads training data from file.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::appendSAT(FILE *outFile)
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->appendSAT(outFile);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method writes training data to file.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::storeSAT(FILE *outFile)
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->storeSAT(outFile);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method writes all its data to file. (number of gaussians, the weight table and the 
/// gaussian parameters)
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::storeData(FILE *outFile)
{
  fwrite(&numberOfGaussians,1,sizeof(numberOfGaussians),outFile);
  for(int i=0;i<numberOfGaussians;i++)
  {
    fwrite(&(myGaussianWeights[i]),1,sizeof(double),outFile);
    myGaussians[i]->storeData(outFile);
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Add all acumulator data to the data of this MixGaussian. The number of gaussians of both
/// models should be the same!
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::train(MixGaussian *trainMG)
{
  if(numberOfGaussians == trainMG->numberOfGaussians)
  {
    for(int i=0;i<numberOfGaussians;i++)
    {
      myGaussians[i]->train(trainMG->myGaussians[i]);
    }    
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Given a single sample-observation (only one Vector), the mixed gaussian is trained. If K-means
/// is used, the observation is used to train only the best fitting Gaussian. K-means is (currently)
/// only used during initialisation of a TrainPhoneModel. If K-means is not used, all gaussians in 
/// the mix are updated with the observation, using the probability that the observation would map 
/// on the perticular gaussian. See Gaussian::train() for more information on training single
/// gaussians. The weight is passed to each of the gaussians (used by Baum-Welch). 
///
/// After feeding all observations to this method, trainFinish() shall be called to finalize the
/// training procedure.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::train(bool useKMeans, double weight, Vector *observation, MixGaussian *doSat)
{
  if(weight > TRAIN_OUTLAYER_P)
  {
    if(numberOfGaussians <= 1)
    {
      if(doSat == NULL)
      {
        myGaussians[0]->train(1.0,observation);
      }
      else
      {
        myGaussians[0]->train(1.0,observation, doSat->myGaussians[0]);
      }
    }
    else if(useKMeans)
    {
      if(observation->getKey() != lastObservation)
      {
        getP(observation);  // Now we have all the gaussian info in fastP_res!!!
      }
      double maxP = 0.0;
      int k = 0;
      for(int i=0;i<numberOfGaussians;i++)
      {
        double P = fastP_res[i];
        if(P > maxP)
        {
          maxP = P;
          k    = i; 
        }
      }
      if(doSat == NULL)
      {
        myGaussians[k]->train(weight,observation);
      }
      else
      {
        myGaussians[k]->train(weight,observation,doSat->myGaussians[k]);
      }
    }
    else
    {
      double factor = 0.0;
      if(observation->getKey() != lastObservation)
      {
        factor = getP(observation);  // Now we have all the gaussian info in fastP_res!!!
      }
      else
      {
        for(int i=0;i<numberOfGaussians;i++)
        {
          factor += fastP_res[i];
        }
        if(isnan(factor) || factor < SMALLEST_NEGATIVE)
        {
          factor = SMALLEST_NEGATIVE;
        }
      }
      if(factor > TRAIN_OUTLAYER_P)
      {
        factor = weight/factor;

        for(int i=0;i<numberOfGaussians;i++)
        {
          double gw = fastP_res[i]*factor;
          if(gw > TRAIN_BIGGEST_P)
          {
            gw = TRAIN_BIGGEST_P;
            printf("  Warning: The weight of a single Gaussian of a Gaussian mixture has been roofed!\n");
          }
          if(gw > TRAIN_OUTLAYER_P)
          {
            if(doSat == NULL)
            {
              myGaussians[i]->train(gw,observation);
            }
            else
            {
              myGaussians[i]->train(gw,observation,doSat->myGaussians[i]);
            }
          }
        }
      }
      else
      {
        printf("  Warning: Total probability of this observation is too small! (skipped, %d, %g)\n",numberOfGaussians,factor);
      }      
    }  
  }
  else
  {
    printf("  Warning: The weight for this mixture is too small! (skipped)\n");
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method checks if all gaussians from the mix are updated by enough (MIN_TRAIN_SAMPLES_PER_GAUSSIAN)
/// training observations. If not, the gaussian is pruned. After this, the mixture weights are updated
/// according to the number of observations used for each gaussian, and the single gaussian training 
/// procedures are finalized by calling Gaussian::trainFinish().
/// This method will return the amount of used training data.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int MixGaussian::trainFinish(bool neverSplit, double minVar)
{   
  double sumD = 0.0;
  double denom;
  int numPruned = 0;
  for(int i=0;i<numberOfGaussians;i++)
  {
    denom = myGaussians[i]->getTrainingDenominator();
    if(!neverSplit && denom <= MIN_TRAIN_SAMPLES_PER_GAUSSIAN)
    {
      PRINTF4("  Pruning Gaussians %d with denom: %g (%g)\n",i,denom,myGaussianWeights[i]);
      delete myGaussians[i];
      numberOfGaussians--;
      myGaussians[i] = myGaussians[numberOfGaussians];
      myGaussians[numberOfGaussians] = NULL;
      myGaussianWeights[i] = myGaussianWeights[numberOfGaussians];
      myGaussianWeights[numberOfGaussians] = 0.0;
      numPruned++;
    }
    sumD += denom;    
  }
  if(numPruned != 0)
  {
    PRINTF3("Number of Gaussians pruned: %d (remaining %d)\n",numPruned, numberOfGaussians);
  }
  if(sumD > 0.0)
  {
    for(int i=0;i<numberOfGaussians;i++)
    {
      myGaussianWeights[i] = myGaussians[i]->getTrainingDenominator() / sumD;
      myGaussians[i]->trainFinish(true,minVar);
    }  
  }  
  else
  {
    return 0;
  }
  initializeFastCalc(false);
  return ((int)sumD);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

double MixGaussian::finishSAT(double minVar)
{
  double res = 0.0;
  for(int i=0;i<numberOfGaussians;i++)
  {
    if(myGaussians[i]->getTrainingDenominator() > MIN_TRAIN_SAMPLES_PER_GAUSSIAN)
    {
      res += myGaussians[i]->trainFinish(true,minVar);
    }
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::mapAdaptMeans()
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->mapAdaptMean();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The SMAPLR adaptation starts with creating one single cluster (node) at which all gaussians in
/// the system will need to register. The method adapt_setInitialNode() will call the Gaussian::adapt_setInitialNode()
/// of all Gaussian objects that are in the mix (of this MixGaussian object) 
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::adapt_setInitialNode(Adapt_AM_TreeNode *node)
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->adapt_setInitialNode(node);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The MixGaussian has registered at least once to an adaptation node when this method is called.
/// The first time adapt_setInitialNode() shall be used. This method will call the Gaussian::adapt_setNode()
/// method of all Gaussian objects in the mix (of this MixGaussian object).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::adapt_setNode()
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->adapt_setNode();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will clear all adaptation settings, so that a new adaptation iteration can be started.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::adapt_clear()
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->adapt_clear();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will write accumulator data to file
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::writeAccumulators(FILE *file, FILE *fileF, FILE *fileST, bool doBinary)
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->writeAccumulators(file,fileF,fileST,doBinary);
  }
  if(fileF != NULL)
  {
    fprintf(file,"\n");
    fprintf(fileF,"\n");
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will read accumulator data from file
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::addAccumulators(FILE *file)
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->addAccumulators(file);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will read accumulator data from file and perform MMI training
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::trainMMI(FILE *fileEnum, FILE *fileDenom)
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->trainMMI(fileEnum,fileDenom);
  }
  initializeFastCalc(false);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::adapt_adaptVar()
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->adapt_adaptVar();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::adapt_setVarTrans()
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->adapt_setVarTrans();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will call the Gaussian::adapt_setHelperMatrices()
/// method of all Gaussian objects in the mix (of this MixGaussian object).  
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::adapt_setHelperMatrices()
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->adapt_setHelperMatrices();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Once all adaptation matrices are calculated, this method will adapt all Gaussians in the mix
/// according to those matrices.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::adapt_adapt()
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->adapt_adapt();
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The mean Vectors of the Gaussian objects in the mix can be restored to the values before 
/// adaptation by calling adapt_unAdapt(). 
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::adapt_unAdapt()
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->adapt_unAdapt();
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// 
/////////////////////////////////////////////////////////////////////////////////////////////////////

double MixGaussian::getCoSim(MixGaussian *mg1, MixGaussian *mg2)
{
  double res = 0.0;
  for(int i=0;i<numberOfGaussians;i++)
  {
    res += myGaussianWeights[i]*myGaussians[i]->getCoSim(mg1->myGaussians[i],mg2->myGaussians[i]);
  }  
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// 
/////////////////////////////////////////////////////////////////////////////////////////////////////

double MixGaussian::getKLDistance(MixGaussian *mg2)
{
  double res = 0.0;
  for(int i=0;i<numberOfGaussians;i++)
  {
    res += myGaussianWeights[i]*myGaussians[i]->getKLDistance(mg2->myGaussians[i]);
  }  
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// 
/////////////////////////////////////////////////////////////////////////////////////////////////////

double MixGaussian::getNormDistance()
{
  double res = 0.0;
  for(int i=0;i<numberOfGaussians;i++)
  {
    res += myGaussianWeights[i]*myGaussians[i]->getNormDistance();
  }  
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// 
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::moveModelGaussians(MixGaussian *model, double factor)
{
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussians[i]->moveModel(model->myGaussians[i],factor);
  }  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method is added for the Train_Speaker_Clustering class.
/// The amount of data assigned to each Gaussian of the mix is counted by count(). Using this 
/// method, the nmbr best gaussians from source are added to the mix. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::addCountedGaussians(MixGaussian *source, int nmbr)
{
  if(nmbr > source->getNumberOfGaussians())
  {
    nmbr = source->getNumberOfGaussians();
  }
  for(int i=0;i<nmbr;i++)
  {
    addGaussian(source->myGaussians[source->myGaussiansSorted[i]],source->myGaussianWeights[source->myGaussiansSorted[i]]);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method is added for the Train_Speaker_Clustering class.
/// The amount of data assigned to each Gaussian of the mix is counted by count(). 
/// Before counting is started, startCount should be called. This method will delete earlier 
/// allocated count-data and initialises the counting parameters.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::startCount()
{
  if(myGaussiansCounts != NULL)
  {
    delete[] myGaussiansCounts;
  }
  if(myGaussiansSorted != NULL)
  {
    delete[] myGaussiansSorted;
  }
  myGaussiansCounts = new double[numberOfGaussians];
  myGaussiansSorted = new    int[numberOfGaussians];
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussiansCounts[i] = 0.0;
    myGaussiansSorted[i] = i;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method is added for the Train_Speaker_Clustering class.
/// The amount of data assigned to each Gaussian of the mix is counted by count(). When done,
/// stopCount() will sort the counted gaussians so that the best X can be chosen.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::stopCount()
{
  // Bubble sort:
  bool swapped = true;
  while(swapped)
  {
    swapped = false;
    for(int i=1;i<numberOfGaussians;i++)
    {
      if(myGaussiansCounts[i] > myGaussiansCounts[i-1])
      {
        double help = myGaussiansCounts[i];
        myGaussiansCounts[i]   = myGaussiansCounts[i-1];
        myGaussiansCounts[i-1] = help;
        swapped = true;
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method is added for the Train_Speaker_Clustering class.
/// The amount of data assigned to each Gaussian of the mix is counted by count(). When done,
/// this method will return the best Gaussian...
/////////////////////////////////////////////////////////////////////////////////////////////////////

int MixGaussian::getBestCount()
{
  int winner = 0;
  double winScore = SMALLEST_NEGATIVE;
  for(int i=0;i<numberOfGaussians;i++)
  {
    if(myGaussiansCounts[i] > winScore)
    {
      winScore = myGaussiansCounts[i];
      winner   = i;
    }
  }
  return winner;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Given the input parameter observation, this method counts the weights of the single Gaussian 
/// objects. These counts are later used by the Train_Speaker_Cluster class.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::count(Vector *observation)
{  
  for(int i=0;i<numberOfGaussians;i++)
  {
    myGaussiansCounts[i] += myGaussianWeights[i] * myGaussians[i]->getP(observation);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Prints the mixture of gaussians to the standard output for debugging purposes.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MixGaussian::printMixGaussian(void)
{
#ifdef DEBUG_ON
  printf("Number of gaussians: %d\n",numberOfGaussians);
  for(int i=0;i<numberOfGaussians;i++)
  {
    printf("(%f)      ",myGaussianWeights[i]);
    myGaussians[i]->printGaussian();
  }  
#endif
}
