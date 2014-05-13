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
#include "multimixgaussian.h"
#include <math.h>
#include "shout-misc.h"

#define MINP      (5.0e-43)

using namespace WriteFileLittleBigEndian;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

MultiMixGaussian::MultiMixGaussian(FeaturePoolInfo *infoBlock, int dim) : MixGaussian(dim)
{
  theSubVector = new Vector(0); // This is an empty vector that will borrow its value from others...
  channelInfo = infoBlock;
  channels    = new MixGaussian*[channelInfo->numberOfChannels];
  int length  = channelInfo->lastComponent[0]+1;
  for(int i=0;i<channelInfo->numberOfChannels;i++)
  {
    if(i>0)
    {
      length = channelInfo->lastComponent[i] - channelInfo->lastComponent[i-1];
    }
    channels[i] = new MixGaussian(length);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

MultiMixGaussian::MultiMixGaussian(FILE *inFile, FeaturePoolInfo *infoBlock, int dim)
{
  theSubVector = new Vector(0); // This is an empty vector that will borrow its value from others...
  channelInfo = infoBlock;
  channels    = new MixGaussian*[channelInfo->numberOfChannels];
  int length  = channelInfo->lastComponent[0]+1;
  for(int i=0;i<channelInfo->numberOfChannels;i++)
  {
    if(i>0)
    {
      length = channelInfo->lastComponent[i] - channelInfo->lastComponent[i-1];
    }
    freadEndianSafe(&(channelInfo->weight[i]), 1,sizeof(channelInfo->weight[i]),inFile);
    channelInfo->weightBackup[i] = channelInfo->weight[i]; 
    channels[i] = new MixGaussian(inFile,length);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MultiMixGaussian::storeData(FILE *outFile)
{

  for(int i=0;i<channelInfo->numberOfChannels;i++)
  {
    fwrite(&(channelInfo->weight[i]), 1,sizeof(channelInfo->weight[i]),outFile);
    channels[i]->storeData(outFile);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

MultiMixGaussian::MultiMixGaussian(MixGaussian *org) : MixGaussian(org->dim())
{
  theSubVector = new Vector(0); // This is an empty vector that will borrow its value from others...
  MultiMixGaussian *realOrg = ((MultiMixGaussian*)org);
  channelInfo   = realOrg->channelInfo;
  channels      = new MixGaussian*[channelInfo->numberOfChannels];

  for(int i=0;i<channelInfo->numberOfChannels;i++)
  {
    channels[i] = new MixGaussian(realOrg->channels[i]);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

MultiMixGaussian::MultiMixGaussian(MixGaussian *model1, MixGaussian *model2) : MixGaussian(model1->dim())
{
  theSubVector = new Vector(0); // This is an empty vector that will borrow its value from others...
  MultiMixGaussian *realModel1 = ((MultiMixGaussian*)model1);
  MultiMixGaussian *realModel2 = ((MultiMixGaussian*)model2);
  channelInfo = realModel1->channelInfo;
  channels      = new MixGaussian*[channelInfo->numberOfChannels];
  for(int i=0;i<channelInfo->numberOfChannels;i++)
  {
    channels[i] = new MixGaussian(realModel1->channels[i],realModel2->channels[i]);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

MultiMixGaussian::MultiMixGaussian(MixGaussian *model1, MixGaussian *model2, double rate) : MixGaussian(model1->dim())
{
  theSubVector = new Vector(0); // This is an empty vector that will borrow its value from others...
  if(fastP_CPandWeight_meanAndVar != NULL)
  {
    delete[] fastP_res;
    delete[] fastP_CPandWeight_meanAndVar;
    fastP_res                    = NULL;
    fastP_CPandWeight_meanAndVar = NULL;
  }

  MultiMixGaussian *realModel1 = ((MultiMixGaussian*)model1);
  MultiMixGaussian *realModel2 = ((MultiMixGaussian*)model2);
  channelInfo = realModel1->channelInfo;
  channels      = new MixGaussian*[channelInfo->numberOfChannels];

  for(int i=0;i<channelInfo->numberOfChannels;i++)
  {
    channels[i] = new MixGaussian(realModel1->channels[i],realModel2->channels[i],rate);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////


MultiMixGaussian::~MultiMixGaussian()
{
  if(channels != NULL)
  {
    for(int i=0;i<channelInfo->numberOfChannels;i++)
    {
      if(channels[i] != NULL)
      {
        delete channels[i];
      }
    }
    delete[] channels;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MultiMixGaussian::mergeMixGaussians(MixGaussian *model1, MixGaussian *model2)
{
  MultiMixGaussian *realModel1 = ((MultiMixGaussian*)model1);
  MultiMixGaussian *realModel2 = ((MultiMixGaussian*)model2);
  for(int i=0;i<channelInfo->numberOfChannels;i++)
  {
    channels[i]->mergeMixGaussians(realModel1->channels[i],realModel2->channels[i]);
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

double MultiMixGaussian::getLogP(Vector *observation)
{
  double res = 0.0;
  for(int i=0;i<channelInfo->numberOfChannels;i++)
  {
    if(channelInfo->weight[i] > 0.0)
    {
      res += channels[i]->getLogP(createVector(observation,i))*channelInfo->weight[i];
    }
  }
  if(isnan(res) || res < MINIMUM_LOGP)
  {
    res = MINIMUM_LOGP;
  }
  return ((double)res);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

double MultiMixGaussian::getP(Vector *observation)
{
  double res = exp(getLogP(observation));
  if(isnan(res) || res < MINP)
  {
    res = MINP;
  }
  return ((double)res);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MultiMixGaussian::shiftBestGaussian(int shiftFactor)
{
  channels[0]->shiftBestGaussian(shiftFactor);
  for(int i=1;i<channelInfo->numberOfChannels;i++)
  {
    channels[i]->shiftBestGaussian(shiftFactor);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool MultiMixGaussian::splitBestGaussian(bool alwaysSplit)
{
  bool res = channels[0]->splitBestGaussian(alwaysSplit);
/*  for(int i=1;i<channelInfo->numberOfChannels;i++)
  {
    channels[i]->splitBestGaussian(alwaysSplit);
  }*/
  return res;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MultiMixGaussian::normalizeWeights(int maxGaussians)
{
  for(int i=0;i<channelInfo->numberOfChannels;i++)
  {
    channels[i]->normalizeWeights(maxGaussians);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MultiMixGaussian::copyGaussians(MixGaussian *destMixGaussian, int maxNmbr)
{
  MultiMixGaussian *realdestMixGaussian = ((MultiMixGaussian*)destMixGaussian);

  for(int i=0;i<channelInfo->numberOfChannels;i++)
  {
    channels[i]->copyGaussians(realdestMixGaussian->channels[i], maxNmbr);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MultiMixGaussian::train(MixGaussian *trainMG)
{
  MultiMixGaussian *realTrainMG = ((MultiMixGaussian*)trainMG);

  for(int i=0;i<channelInfo->numberOfChannels;i++)
  {
    channels[i]->train(realTrainMG->channels[i]);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MultiMixGaussian::train(bool useKMeans, double weight, Vector *observation, MixGaussian *doSat)
{
  for(int i=0;i<channelInfo->numberOfChannels;i++)
  {
    channels[i]->train(useKMeans,weight,createVector(observation,i),doSat);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

int MultiMixGaussian::getNumberOfGaussians(void)
{
  return channels[0]->getNumberOfGaussians();
}
    
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MultiMixGaussian::addGaussian(Vector *vector)
{
  for(int i=0;i<channelInfo->numberOfChannels;i++)
  {
    channels[i]->addGaussian(createVector(vector,i));
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

int MultiMixGaussian::trainFinish(bool neverSplit, double minVar)
{
  int res = channels[0]->trainFinish(neverSplit, minVar);
  for(int i=1;i<channelInfo->numberOfChannels;i++)
  {
    channels[i]->trainFinish(neverSplit,minVar);
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

Vector *MultiMixGaussian::createVector(Vector *totVect, int nr)
{   
  int first   = 0;
  int last    = channelInfo->lastComponent[nr];
  if(nr > 0)
  {
    first = channelInfo->lastComponent[nr-1]+1;
  }
  theSubVector->borrowVector(totVect,first,last);
  return theSubVector;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// GetLogP() returns the chance of an observation beeing produced by this mixture of gaussians (in log)
/////////////////////////////////////////////////////////////////////////////////////////////////////


double MultiMixGaussian::getLookaheadLogP(double *vectorList, int timeStamp, bool doSecondHalf)
{
  double res = 0.0;
  double *vectorPoint = vectorList;
  for(int i=0;i<channelInfo->numberOfChannels;i++)
  {
    if(channelInfo->weight[i] > 0.0)
    {
      res += channels[i]->getLookaheadLogP(vectorPoint,timeStamp,doSecondHalf)*channelInfo->weight[i];
    }
    vectorPoint = vectorList + HALF_GMMLOOKAHEAD*(channelInfo->lastComponent[i]+1);
  }
  if(isnan(res) || res < MINIMUM_LOGP)
  {
    res = MINIMUM_LOGP;
  }
  return ((double)res);
}



