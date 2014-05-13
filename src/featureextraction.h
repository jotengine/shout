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
#ifndef FEATUREEXTRACTION_H
#define FEATUREEXTRACTION_H

// Stuff that is not meant to be changed (therefore not in standard.h)
#define FEA_ADDEDZEROS           (0)
#define FEA_WINDOWSIZE         (512)
#define FEA_REMEMBERWINSIZE    (352)
#define FEA_ADDWINSIZE         (160)
#define FEA_FREQWINDOWSIZE     (256)
#define FEA_BYTES_PER_SAMPLE     (2)
#define FEA_DELTA                (2)

#define DEFAULT_MEL_BANKLENGTH  (24)
#define LARGE_MEL_BANKLENGTH    (24)

#define SWITCH_MELANDCEPLIFTER  (14)
#define DEFAULT_CEP_LIFTER       (0)
#define LARGE_CEP_LIFTER         (0)
#define MEL_FREQ(a)             (2595.0*log10(1.0+a/700.0))
#define FREQ(a)             (8000.0*(a)/((double)FEA_FREQWINDOWSIZE))

#define MIN_MEL_WINDOW        (94.0)
#define MAX_MEL_WINDOW      (6438.0)

#define MIN_CTS_MEL_WINDOW    (94.0)
#define MAX_CTS_MEL_WINDOW  (3800.0)

#include <stdio.h>
#include "standard.h"
#include "FFTReal.h"
#include "gaussian.h"
#include "mixgaussian.h"

typedef enum
{
  FEATURENORM_SEGMENT  = 0,
  FEATURENORM_CLUSTER  = 1,
  FEATURENORM_ALL      = 2
} FEATURENORM_TYPE;

struct WaveHeaderType
{
  char         chunk1[20];
  char         audioFormat[2];
  char         channels[2];
  char         sampleRate[4];
  char         byteRate[4];
  char         blockAlign[2];
  char         bitsPerSample[2];
  char         chunk2[8];
};

struct NormData
{
  Gaussian     *cmvnCluster;
  Vector      **gaussianizationSource;
  Vector     ***gaussianizationDest;
  int           gaussianizationPointer;
  Vector      **histogram;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief This class handles all feature extraction issues.
///
/// The following feature extraction steps are taken:
/// - Create overlapping windows: default: 32 ms windows, every 10 ms.
/// - DC offset removal. (is this needed? Does the pre-emphesis solve this?)
/// - Pre-emphesis
/// - Apply Hamming window
/// - Calculate energy
/// - Fast Fourier Transformation (magnitude)
/// - Vocal Tract Length Normalization
/// - Mel bank filtering (default 20 banks, logaritmic)
/// - DCT transformation creating MFCC coefficients
/// - Cepstrum liftering
/// - Normalize the coefficients: mean substraction
/// - Calculate delta's and delta-delta's
/// 
/// This class can do feature extraction for a single audio file (16K16, raw PCM audio) or for a 
/// batch file.
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

#include "vector.h"

class FeatureExtraction : public Gaussian
{
  protected:
    int       melIndexAboveCTS;
    FILE     *audioFile;
    FILE     *audioFileOut;
    bool      offset;
    bool      onlyMEL;
    int       CEP_LIFTER;
    int       MEL_BANKLENGTH;
    Vector   *cepLifterFact;
    float     fea_delta_noemer;
    FFTReal  *fftProc;
    float     hammingTemplate[FEA_WINDOWSIZE];
    short int rememberBuffer [FEA_REMEMBERWINSIZE];
    int       melTemplate    [LARGE_MEL_BANKLENGTH+2];      
    float     fftWindowSpectralSubtract[FEA_FREQWINDOWSIZE];

    Vector  **featureVector;
    Vector  **gaussianization;
    int       nrFrames;
    int       mfccSize;
    int       vectorSize;
    int       useDeltas;
    int       useZeroCross;
    int       useEnergy;
    int       nonDeltaSize;
    int       numberOfFeaturesProcessed;
    int       nrNormClusters;
    NormData *normData;
    bool      doCMN;
    bool      doCVN;
    bool      doSqrt10;
    bool      performCTS;
    FEATURENORM_TYPE normType;

  public:
    Vector  **inverseHistogram;

    FeatureExtraction             (const char *audioIn, const char *audioOut, int mfccSize, int useEnergy, 
                                   int useZeroCross, int useDeltas, double vtln, bool useSqrt10 = false,
                                   bool offs = false, bool onlyPrepare = false, bool onlyMEL = false, bool pCTS = false);
   ~FeatureExtraction             ();

    // INLINE FUNCTIONS FOR CODE OPTIMALIZATION:
    inline int getVectorSize      ()          const {return vectorSize; };
    inline int getPoolSize        ()          const {return nrFrames; };
    inline Vector **getPool       ()          const {return featureVector; };

    Vector *getOnlineVector       (double *silThreshold, bool isSil);
    void    setVTLN               (double vtln);
    void    setBackgroundHistNorm (Vector **hist, int numBins);
    bool    createFeaturesUntilFrame(int lastFrame, int normID = 0);
    void    finishExtraction      ();
    void    doFeatureExtraction   (const char *audioIn, const char *audioOut, bool onlyPrepare = false);    
    void    setNormalization      (FEATURENORM_TYPE normT, int nrClusters, bool doMean, bool doVar, Vector **hist);
    void    finishClusterCMVN     ();
    void    performClusterCMVNUntilFrame(int lastFrame, int cmvnID);
    void    performHistNormUntilFrame(int lastFrame, int histID);
    void    initializeHistNorm    ();
    void    setHistNormModel      (int normID, Vector **h);
    void    setAudioFileOut       (const char *fileName);
    void    performPCA            (Vector **pca, int len);
    void    setOnlyMEL            ();
    void    storeFeatureVectors   (FILE *file);
  protected:
    bool    readAudioWindow       (short int *audioWindow);
    void    calculateDelta        (int time, int offset);
    Vector *createMfccFrame       (int normID, int spectralSubtraction = -1);
};

#endif
