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
#ifndef MULTIMIXGAUSSIAN_H
#define MULTIMIXGAUSSIAN_H

#include "mixgaussian.h"
#include "featurepool.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

class MultiMixGaussian : public MixGaussian
{
  protected:
    FeaturePoolInfo  *channelInfo;
    MixGaussian     **channels;
    Vector           *theSubVector;

  public:
   ~MultiMixGaussian                   ();
    MultiMixGaussian           (FeaturePoolInfo *infoBlock, int dim = ASR_DEFAULT_VECTORSIZE);
    MultiMixGaussian           (FILE *inFile, int dim);
    MultiMixGaussian           (FILE *inFile, FeaturePoolInfo *infoBlock, int dim);
    MultiMixGaussian           (MixGaussian *org);
    MultiMixGaussian           (MixGaussian *model1, MixGaussian *model2);
    MultiMixGaussian           (MixGaussian *model1, MixGaussian *model2, double rate);

  protected:
    Vector *createVector(Vector *totVect, int nr);
  public:
    virtual void             storeData          (FILE *outFile);
    
    virtual void     mergeMixGaussians          (MixGaussian *model1, MixGaussian *model2);
    virtual double getP                       (Vector *observation);
    virtual double getLogP                    (Vector *observation);
    
    virtual int      getNumberOfGaussians       (void);
    virtual void     addGaussian                (Vector *vector);

    virtual void     train                      (MixGaussian *trainMG);
    virtual void     train                      (bool useKMeans, double weight, Vector *observation, MixGaussian *doSat = NULL);
    
    virtual int      trainFinish                (bool neverSplit = false, double minVar = GAUSSIAN_VARIANCE_MINIMUM);
    virtual bool     splitBestGaussian          (bool alwaysSplit = false);
    virtual void     shiftBestGaussian          (int shiftFactor);
    virtual void     copyGaussians              (MixGaussian *destMixGaussian, int maxNmbr);
    virtual void     normalizeWeights           (int maxGaussians = -1);
    virtual double getLookaheadLogP           (double *vectorList, int timeStamp, bool doSecondHalf);
};

#endif
