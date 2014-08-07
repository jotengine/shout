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

#ifndef STANDARD_H
#define STANDARD_H

#include <stdio.h>
#include <cassert>

// #define NDEBUG

#ifdef __APPLE__
# define fopen64 fopen
# define open64 open
#endif

#define SMALLEST_NEGATIVE            (-9.0e60)

#define MAXSTR                       (1000)

#define USER_ERROR(a)                 printf("\nERROR: %s\n\n",a); exit(1);
#define USER_ERROR2(a,b)              printf("\nERROR: %s %s\n\n",a,b); exit(1);
#define USER_ERROR2_NUM(a,b)          printf("\nERROR: %s %d\n\n",a,b); exit(1);
#define USER_WARNING(a)               printf("WARNING: %s\n",a); fflush(stdout);
#define USER_WARNING2(a,b)            printf("WARNING: %s %s\n",a,b); fflush(stdout);

#define CURRENT_VERSION "2010-RELEASE, 01-12-2010\n"

#define WELCOME_HEADER1 "############################################################################\n\
### This application is part of the 'SHoUT LVCS Recognition toolkit'.    ###\n\
### SHoUT is developed by Marijn Huijbregts,                             ###\n\
### Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010                     ###\n\
### For more info visit: http://sourceforge.net/projects/shout-toolkit/  ###\n\
############################################################################\n\
### This program is free software; you can redistribute it and/or modify ###\n\
### it under the terms of the GNU General Public License as published by ###\n\
### the Free Software Foundation; version 2 of the License.              ###\n\
###                                                                      ###\n\
### This program is distributed in the hope that it will be useful,      ###\n\
### but WITHOUT ANY WARRANTY; without even the implied warranty of       ###\n\
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        ###\n\
### GNU General Public License for more details.                         ###\n\
###                                                                      ###\n\
### You should have received a copy of the GNU General Public License    ###\n\
### along with this program; if not, write to the                        ###\n\
### Free Software Foundation, Inc.,                                      ###\n\
### 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.            ###\n\
############################################################################\n\
### Version: "

#define WELCOME_HEADER2 "############################################################################\n"

// #define DEBUG_ON

#ifdef DEBUG_ON
  #define PRINTF(a)           printf(a)
  #define PRINTF2(a,b)        printf(a,b)
  #define PRINTF3(a,b,c)      printf(a,b,c)
  #define PRINTF4(a,b,c,d)    printf(a,b,c,d)
  #define PRINTF5(a,b,c,d,e)  printf(a,b,c,d,e)
  #define PRINTF6(a,b,c,d,e,f)  printf(a,b,c,d,e,f)
  #define PRINTF7(a,b,c,d,e,f,g)  printf(a,b,c,d,e,f,g)
#else
  #define PRINTF(a)
  #define PRINTF2(a,b)
  #define PRINTF3(a,b,c)
  #define PRINTF4(a,b,c,d)
  #define PRINTF5(a,b,c,d,e)
  #define PRINTF6(a,b,c,d,e,f)
  #define PRINTF7(a,b,c,d,e,f,g)
#endif


// #define HISTOGRAM_NORM_SPECTRUM
// #define PRINT_HISTNORM (true)


//////////////////////////  DEBUGGING AND ANALYSIS OPTIONS: //////////////////////////

// #define DEBUG_MEMORY
//#define SEARCH_STATISTICS_ON

#ifdef DEBUG_MEMORY
  extern int plrCount;
  extern int wlrCount;
  extern int tokCount;
  extern int trackCount;
  extern int gaussianCount;
  extern int lmlaCount;
  extern int pdfCount;
  #define COUNTER_PRINT()   printf("# TOKEN COUNTERS:\n# wlrCount: %d\n# tokCount: %d\n# plrCount: %d\n# trackCount: %d\n# gaussianCount: %d\n# LMLA: %d\n# pdfCount: %d\n",wlrCount,tokCount,plrCount, trackCount,gaussianCount,lmlaCount,pdfCount);
  #define COUNTER_PLUS(a)   a++;
  #define COUNTER_MINUS(a)  a--;
#else
  #define COUNTER_PRINT()   ;
  #define COUNTER_PLUS(a)   ;
  #define COUNTER_MINUS(a)  ;
#endif


//////////////////////////  STANDARD (MATH) CONSTANTS: ////////////////

#define TWO_PI                          (6.283185307179586)
#define PI                              (3.141592653589793)

//////////////////////////  FEATURE VECTOR OPTIONS: ///////////////////

#define ADAPT_DEFAULT_VECTORSIZE                       (50)

#define ASR_DEFAULT_VECTORSIZE                         (39)

#define ASR_DEFAULT_MFCC_NUMBER                        (12)
#define ASR_DEFAULT_MFCC_ENERGY                         (1)
#define ASR_DEFAULT_MFCC_ZEROCROSS                      (0)
#define ASR_DEFAULT_MFCC_DELTA                          (2)
#define ASR_DEFAULT_CMN                              (true)
#define ASR_DEFAULT_CVN                              (true)
#define ASR_DEFAULT_FEATURENORM       (FEATURENORM_CLUSTER)
#define ASR_DEFAULT_SQRT10                          (false)

#define SPKREC_DEFAULT_MFCC_NUMBER                     (12)
#define SPKREC_DEFAULT_MFCC_ENERGY                      (1)
#define SPKREC_DEFAULT_MFCC_ZEROCROSS                   (0)
#define SPKREC_DEFAULT_MFCC_DELTA                       (1)
#define SPKREC_DEFAULT_CMN                           (true)
#define SPKREC_DEFAULT_CVN                           (true)
#define SPKREC_DEFAULT_FEATURENORM    (FEATURENORM_CLUSTER)
#define SPKREC_DEFAULT_SQRT10                       (false)

#define DIARIZATION_DEFAULT_MFCC_NUMBER                (19)
#define DIARIZATION_DEFAULT_MFCC_ENERGY                 (0)
#define DIARIZATION_DEFAULT_MFCC_ZEROCROSS              (0)
#define DIARIZATION_DEFAULT_MFCC_DELTA                  (0)
#define DIARIZATION_DEFAULT_CMN                      (true)
#define DIARIZATION_DEFAULT_CVN                      (true)
#define DIARIZATION_DEFAULT_FEATURENORM   (FEATURENORM_CLUSTER)
#define DIARIZATION_DEFAULT_SQRT10                  (false)

#define VTLN_DEFAULT_MFCC_NUMBER                       (12)
#define VTLN_DEFAULT_MFCC_ENERGY                        (1)
#define VTLN_DEFAULT_MFCC_ZEROCROSS                     (0)
#define VTLN_DEFAULT_MFCC_DELTA                         (0)
#define VTLN_DEFAULT_CMN                             (true)
#define VTLN_DEFAULT_CVN                            (false)
#define VTLN_DEFAULT_FEATURENORM      (FEATURENORM_CLUSTER)
#define VTLN_DEFAULT_SQRT10                         (false)

#define SAD_DEFAULT_MFCC_NUMBER                        (12)
#define SAD_DEFAULT_MFCC_ENERGY                         (0)
#define SAD_DEFAULT_MFCC_ZEROCROSS                      (1)
#define SAD_DEFAULT_MFCC_DELTA                          (2)
#define SAD_DEFAULT_CMN                              (true)
#define SAD_DEFAULT_CVN                             (false)
#define SAD_DEFAULT_FEATURENORM           (FEATURENORM_ALL)
#define SAD_DEFAULT_SQRT10                          (false)

// The max-vector size is including energy and delta, d-d.
#define FEA_MAXVECTORSIZE                             (100)
#define FEA_ONLINE_BUFFER                             (100)

#define FEATURE_NORM_INVERSE_HISTOGRAM_BINS         (50000)
#define FEATURE_NORM_HISTOGRAM_BINS                 (50002)
#define FEATURE_NORM_TIMESVAR                           (4)

//////////////////////////  DECODER OPTIONS: //////////////////////////

struct DecoderSettings
{
  bool   doBeam;
  bool   doEndState;
  bool   doHist;
  double prune_Beam;
  double prune_StateBeam;
  double prune_EndStateBeam;
  bool   prune_Lmla;
  int    prune_HistState;
  int    prune_Hist;
  double weights_LmScale;
  double weights_TransPenalty;
  double weights_SilPenalty;
};

// Can be changed:
#define LM_SCALE                                     (27.5)
#define TRANS_PENALTY                                 (0.0)
#define SIL_PENALTY                                   (0.0)
#define BEAM                                        (200.0)
#define STATE_BEAM                                   (75.0)
#define END_STATE_BEAM                               (50.0)
#define LMLA                                           "ON"
#define HISTOGRAM_STATE_PRUNING                       (160)
#define HISTOGRAM_PRUNING                           (35000)

#define PERFORM_PRUNING_HISTOGRAM                    (true)
#define PERFORM_PRUNING_STATE                        (true)
#define PERFORM_PRUNING_ENDSTATE                     (true)

// Can not be changed:

#define LMLA_UNIGRAM_ONLY

#define SHORT_WORD_PENALTY                            (0.0)
#define SHORT_WORD_LENGTH                               (3)

// DO NOT CHANGE THIS VALUE!! The architecture is not developed to do so...
#define NONSIL_NUMBER_OF_STATES                         (3)

// Can not be changed:
#define LMLA_CACHE                                    (200)
#define LMLA_INTERVAL                                  (50)
#define NR_OF_HISTOGRAM_BINS                          (100)

// Can not be changed:
#define NBEST_DEPTH                                     (9)
#define NBEST_MAXNR_SCORE                            (1000)

//////////////////////  LANGUAGE MODEL OPTIONS: ///////////////////////

// LM_NGRAM_DEPTH must be bigger than 0 and (for now) smaller than 4...
// (trigram = 2)
#define LM_NGRAM_DEPTH                                  (2)
#define DO_MEMORY_MAPPED_FILE_LM                    (false)

//////////////////////////  VTLN OPTIONS: //////////////////////////

#define VTLN_STEPS                                     (11)
#define VTLN_NEUTRAL                                    (5)
#define VTLN_START                                   (0.80)
#define VTLN_STEPSIZE                                (0.04)
#define VTLN_NR_TRAIN_SENTENCES                        (30)

//////////////////////////  TRAINER OPTIONS: //////////////////////////

#define PHONEFILE                             "/phones.lst"

#define TRAIN_MAX_STATE_LENGTH                        (100)
#define TRAIN_MAX_TOTAL_LENGTH                       (1000)
#define TRAIN_MIN_TOTAL_LENGTH                          (2)

#define TRAIN_MAX_FEATURE_POOL_SIZE              (14000000)
#define TRAIN_ID1                                       (0)
#define TRAIN_CLUSTER_ID1                               (1)
#define TRAIN_CLUSTER_ID2                               (2)
#define TRAIN_CLUSTER_ID3                               (3)
#define TRAIN_CONTEXTKEY                                (1)
#define TRAIN_CONTEXTLEFT                               (2)
#define TRAIN_CONTEXTRIGHT                              (3)

#define MINIMUM_TRAIN_IMPROVEMENT                  (0.0005)
#define MINIMUM_TRAIN_IMPROVEMENT_FINAL             (0.003)
#define MINIMUM_TRAIN_ITERATIONS                        (3)
#define MAXIMUM_TRAIN_ITERATIONS                       (30)
#define MIN_TRAIN_SAMPLES_PER_GAUSSIAN               (20.0)


#define TRAIN_BAUM_WELCH_FREQUENCY                      (3)
#define TRAIN_BIGGEST_P                           (1.0e300)
#define TRAIN_OUTLAYER_P                         (1.0e-300)
#define TRAIN_WEIGHT_MINIMUM                     (1.0e-100)

#define GAUSSIAN_VARIANCE_MINIMUM                    (0.01)

#define MINIMUM_LOGP                                (-1000)

#define TRANSITION_MINIMUM                        (0.00001)

#define MIN_SAMPLES_PER_CLUSTER                       (250)
#define MAXIMUM_NUMBER_OF_STATE_CLUSTERS              (250)
#define MINIMUM_CLUSTER_IMPROVEMENT                  (0.10)
#define MAX_NUMBER_OF_MIXTURES_PER_PHONE             (1000)
#define MAX_NUMBER_OF_DECISION_RULES                  (100)

//////////////////////////  ADAPTATION OPTIONS:  ////////////////////////

#define ADAPT_MIN_GAUSSIANS_PER_NODE               (2000.0)
#define ADAPT_MAX_TREE_DEPTH                            (7)

//////////////////////////  SPEAKER CLUSTERING:  ////////////////////////

#define FEATURE_SPEAKER_SIZE                           (12)


#define MINIMUM_TRAIN_IMPROVEMENT_CLUSTER           (0.001)

#define CLUSTER_MIN_NUMBER_FRAMES                     (250)
#define CLUSTER_LM_SCALE                              (0.0)
#define CLUSTER_TRANS_PENALTY                         (0.0)
#define CLUSTER_SIL_PENALTY                           (0.0)

#define DIARIZATION_FAST_MERGE                          (2)
#define DIARIZATION_MAX_CLUSTERS                       (40)

#define CLUSTER_BEAM                               (9.0e30)
#define CLUSTER_STATE_BEAM                            (1.0)
#define CLUSTER_END_STATE_BEAM                        (1.0)
#define CLUSTER_LMLA                                (false)
#define CLUSTER_HISTOGRAM_STATE_PRUNING                 (1)
#define CLUSTER_HISTOGRAM_PRUNING                   (40000)

#define CLUSTER_PERFORM_PRUNING_HISTOGRAM            (true)
#define CLUSTER_PERFORM_PRUNING_STATE                (true)
#define CLUSTER_PERFORM_PRUNING_ENDSTATE             (true)


#define INITIAL_DATA_PER_GAUSSIAN                     (800)
#define MIN_NUMBER_OF_GAUSSIANS                         (5)

// FOR THE SPLIT VARIANT:
#define DATA_PER_MERGE                              (18000)

#define MAX_DATA_PER_CLUSTER                        (20000)
#define MIN_SPLIT_GAUSSIANS                             (4)
#define MAX_SPLIT_GAUSSIANS                             (8)
#define MAX_DATA_PER_GAUSSIAN                         (300)

//////////////////////////  CONFIDENCE MEASURES:  ////////////////////////

#define PHONE_CONFIDENCE_THRESHOLD                  (-50.0)

//////////////////////////  CLUSTERED AMs:  ////////////////////////

#define CLUSTEREDAM_MINIMUM_DATA_SIZE                (8000)

#define PHONEFILE                             "/phones.lst"
#define ARTICULATORY_FILENAME  "/articulatory.decisiontree"


////////////////////// MONITORING: ////////////////////////////////////

#define TOKENDISTRIBUTION_TREEDEPTH                     (3)

////////////////////// Fast P /////////////////////////////////

#define MAX_NR_THREADS                                  (4)
#define GMMLOOKAHEAD                                   (10)
#define HALF_GMMLOOKAHEAD                               (5)

////////////////////// MISC ////////////////////////////////////

#define ABSOLUTE(a)          if(a<0.0){a=-1.0*a;}
#define ABSOLUTE_INT(a)      if(a<0)  {a=-1*a;}
#define MODULO(a,b,c)        a=b%c; if(a<0){a+=c;}

#endif // STANDARD_H

