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
  
#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include "standard.h"
#include "vector.h"
#include "adapt_am_treenode.h"

class Adapt_AM_TreeNode;

struct TrainGaussian
{
  Vector            *numeratorMean;
  Vector            *numeratorVariance;
  double             trainingDenominator;
  double             lastTrainingWeight;
  double             trainP;
  Vector            *unAdaptedMeanVector;    
  Adapt_AM_TreeNode *adaptNode;  
  int                trainingSamples;
};

////////////////////////////////////////////////////////////////////////////
/// \brief Using Gaussian objects, single gaussians can be trained and used. 
///
/// In order to form a gaussian, the user must provide Gaussian objects with
/// training samples. The object will set its mean Vector and variance Vector
/// according to these samples.
/// With getP() the user can request the likelihood that the gaussian would
/// produce a particular observation.
////////////////////////////////////////////////////////////////////////////

class Gaussian
{
  protected:
    double  CP;                     ///< CP is the part of the normal distribution calculation that can be done before-hand.
    Vector *varianceVector;
    Vector *meanVector;
    TrainGaussian *trainPars;
    TrainGaussian *mmi_trainPars;
    double  helper_meanJKProd;
    // Methods:   
    void     calculateCP                ();
    void     enterTrainingMode          (bool mmiTraining = false);
    
  public:  

    
    // Constructor and destructor:
             Gaussian                   (int dim = ASR_DEFAULT_VECTORSIZE);
             Gaussian                   (Gaussian *copyFrom, bool shiftGaussian = true);
             Gaussian                   (Gaussian *model1, Gaussian *model2, double rate);
             Gaussian                   (FILE *inFile, int dim);
            ~Gaussian                   ();
    
    // INLINE FUNCTIONS FOR CODE OPTIMALIZATION:
    inline int getNrTrainingSamples     ()  const {assert(trainPars != NULL); return trainPars->trainingSamples;}
    inline int dim                      ()  const {assert(varianceVector != NULL); return varianceVector->len();}
    inline Vector *getMean              ()  const {assert(meanVector != NULL);     return meanVector;}
    inline double  getCP                ()  const {return CP;}
    inline Vector *getVariance          ()  const {assert(varianceVector != NULL); return varianceVector;}
    
    
    // Methods:
    virtual double getP                 (Vector *observation);
    virtual double getLogP              (Vector *observation);
    void     setVariance                (const Vector *v);
    void     setMean                    (Vector *v);
    
    void     storeData                  (FILE *outFile);
    void     storeSAT                   (FILE *outFile);
    void     appendSAT                  (FILE *outFile);

    void     train                      (Gaussian *trainG);
    void     train                      (double factor, Vector *observation, Gaussian *doSat = NULL);
    double   trainFinish                (bool cleanUp = true, double minVar = GAUSSIAN_VARIANCE_MINIMUM);
    void     shiftMean                  (int shiftFactor = 1);    
    double   getTrainingDenominator     (void);
  
    double   getCoSim                   (Gaussian *g1, Gaussian *g2);
    double   getKLDistance              (Gaussian *g2);
    double   getNormDistance            ();
    
    void     moveModel                  (Gaussian *model, double factor);

    void     mapAdaptMean               ();
    void     adapt_setInitialNode       (Adapt_AM_TreeNode *node);
    void     adapt_setNode              ();
    void     adapt_setHelperMatrices    ();
    void     adapt_setVarTrans          ();
    void     adapt_adaptVar             ();
    void     adapt_adapt                ();
    void     adapt_unAdapt              ();
    void     adapt_clear                ();
    void     writeAccumulators          (FILE *file, FILE *fileF = NULL, FILE *fileST = NULL, bool doBinary = false);
    void     trainMMI                   (FILE *fileEnum, FILE *fileDenom);
    void     addAccumulators            (FILE *file);

    void     printModel                 (FILE *fileMean, FILE *fileVariance);
    void     printInfo                  (Vector *v);
    void     printGaussian              (void);

};

#endif // GAUSSIAN_H
