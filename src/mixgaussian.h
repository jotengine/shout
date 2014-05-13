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
 
#ifndef MIXGAUSSIAN_H
#define MIXGAUSSIAN_H

#include "standard.h"
#include "gaussian.h"
#include "adapt_am_treenode.h"

class Gaussian;
class Adapt_AM_TreeNode;

////////////////////////////////////////////////////////////////////////////
/// \brief The MixGaussian class handles mixture sets of gaussians.
///
/// Gaussian objects offer functionality for training single gaussian parameter
/// values and determining chances given specific observations. The MixGaussian
/// class makes it possible to train and use mixtures of gaussians. 
/// The user may increase the number of gaussians in the mixture using splitBestGaussian().
/// This method splits the Gaussian that has been trained with the most observation
/// samples.
////////////////////////////////////////////////////////////////////////////

class MixGaussian
{
  protected:
    // Attributes:       
    Gaussian **myGaussians;
    double    *myGaussianWeights;
    double     trainP;
    int        numberOfGaussians;
    bool       noSplitting;
    
    double    *myGaussiansCounts;
    int       *myGaussiansSorted;

    double   *fastP_res;
    void     *lastObservation;
    int       vectorDimension;
    double    lookaheadRes[GMMLOOKAHEAD];
    int       lookaheadPointer;

  public:
    double *fastP_CPandWeight_meanAndVar;
    // Methods:   
  protected:
  
  public:  
    // Constructor and destructor:
             MixGaussian                (int dim = ASR_DEFAULT_VECTORSIZE);
             MixGaussian                (FILE *inFile, int dim, bool noReTrainig = false);
             MixGaussian                (MixGaussian *org);
             MixGaussian                (MixGaussian *model1, MixGaussian *model2);
             MixGaussian                (MixGaussian *model1, MixGaussian *model2, double rate);
   virtual ~MixGaussian                ();
    
    // INLINE FUNCTIONS FOR CODE OPTIMALIZATION:
    inline int dim                      () const {return vectorDimension;}
    // Methods:
    virtual void     mergeMixGaussians          (MixGaussian *model1, MixGaussian *model2);
    virtual int      getNumberOfGaussians       (void);
            void     setVariance                (const Vector *v);
    virtual double getP                       (Vector *observation);
    virtual double getLogP                    (Vector *observation);

    virtual double getLookaheadLogP           (double *vectorList, int time, bool doSecondHalf);


            void     getSumHist                 (Vector **histogram);

    virtual void     train                      (MixGaussian *trainMG);
    virtual void     train                      (bool useKMeans, double weight, Vector *observation, MixGaussian *doSat = NULL);
    void     mapAdaptMeans              ();
    void     adapt_setInitialNode       (Adapt_AM_TreeNode *node);
    void     adapt_setHelperMatrices    ();
    void     adapt_adapt                ();
    void     adapt_unAdapt              ();
    void     initializeFastCalc         (bool onlyOnce);

    virtual int      trainFinish                (bool neverSplit = false, double minVar = GAUSSIAN_VARIANCE_MINIMUM);
    virtual double   finishSAT                  (double minVar = GAUSSIAN_VARIANCE_MINIMUM);
    virtual bool     splitAllGaussians          (bool alwaysSplit = false);
    virtual bool     splitBestGaussian          (bool alwaysSplit = false);
    virtual void     shiftBestGaussian          (int shiftFactor);
    virtual void     copyGaussians              (MixGaussian *destMixGaussian, int maxNmbr);
    virtual void     addGaussian                (Gaussian *gaussian, double weight);
    virtual void     addGaussian                (Vector *vector);
    virtual void     normalizeWeights           (int maxGaussians = -1);
    virtual void     storeData                  (FILE *outFile);
    void     storeSAT                   (FILE *outFile);
    void     appendSAT                  (FILE *outFile);
    void     adapt_setNode              ();
    void     adapt_clear                ();
    void     writeAccumulators          (FILE *file, FILE *fileF = NULL, FILE *fileST = NULL, bool doBinary = false);
    void     addAccumulators            (FILE *file);
    void     trainMMI                   (FILE *fileEnum, FILE *fileDenom);

    void     adapt_setVarTrans          ();
    void     adapt_adaptVar             ();
    
    double   getCoSim                   (MixGaussian *mg1, MixGaussian *mg2);
    double   getKLDistance              (MixGaussian *mg2);
    double   getNormDistance            ();
        
    void     moveModelGaussians         (MixGaussian *model, double factor);
    
    void startCount                     ();
    void count                          (Vector *observation);
    void stopCount                      ();
    int  getBestCount                   ();
    void addCountedGaussians            (MixGaussian *source, int nmbr);
    #ifdef SPLIT_CLUSTERING
    void fillDistanceArray              (int *distA);
    #endif

    void printModel(FILE *fileMean, FILE *fileVariance, FILE *fileWeight);
    int  printInfo (Vector *v);
    void     printMixGaussian           (void);
};

#endif // MIXGAUSSIAN_H
