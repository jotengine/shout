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
#ifndef SHOUTCONFIG_H
#define SHOUTCONFIG_H

#include <stdio.h>
#include "stringlookup.h"
#include "standard.h"

/////////////////////////////////////// global:
#define ARGUMENT_HELP           (0)
#define ARGUMENT_CONFIG         (1)
#define ARGUMENT_CREATE_CONFIG  (2)
/////////////////////////////////////// shout_segment
#define ARGUMENT_AUDIOFILE      (3)
#define ARGUMENT_AM_SEGMENT     (4)
#define ARGUMENT_META_OUT       (5)
#define ARGUMENT_LABEL          (6)
#define ARGUMENT_UEM            (7)
#define ARGUMENT_CTS_IN         (8)
#define ARGUMENT_CTS_OUT        (9)
#define ARGUMENT_WIDEN_SEGMENTS (10)
#define ARGUMENT_ENERGY_BASED   (11)

/////////////////////////////////////// shout_cluster
// #define ARGUMENT_AUDIOFILE
// #define ARGUMENT_LABEL
// #define ARGUMENT_META_OUT
// #define ARGUMENT_WIDEN_SEGMENTS
#define ARGUMENT_META_IN        (12)
#define ARGUMENT_SECOND_STREAM  (13)
#define ARGUMENT_FAST_MERGE     (14)
#define ARGUMENT_MAX_CLUSTERS   (15)
#define ARGUMENT_THIRD_STREAM   (16)

/////////////////////////////////////// shout_vtln
// #define ARGUMENT_AUDIOFILE
// #define ARGUMENT_META_IN
// #define ARGUMENT_META_OUT
#define ARGUMENT_PERLINE        (17)
#define ARGUMENT_AM_VTLN        (18)
#define ARGUMENT_AM_HISTNORM    (19)

/////////////////////////////////////// Shout:
// #define ARGUMENT_AUDIOFILE
// #define ARGUMENT_META_IN
#define ARGUMENT_LM_BIN         (20)
#define ARGUMENT_DCT_BIN        (21)
#define ARGUMENT_PHONEBASED     (22)
#define ARGUMENT_LATTICE        (23)
#define ARGUMENT_NBEST          (24)
#define ARGUMENT_ALIGN          (25)
#define ARGUMENT_LATTICERESCORE (26)
#define ARGUMENT_ANALYSE        (27)
#define ARGUMENT_TOKENDIST      (28)
#define ARGUMENT_PHONECONF      (29)
#define ARGUMENT_BACKDCT_BIN    (30)
#define ARGUMENT_TRAINDIR       (31)
#define ARGUMENT_XML            (32)
#define ARGUMENT_LM_SCALE       (33)
#define ARGUMENT_TRANS_PENALTY  (34)
#define ARGUMENT_SIL_PENALTY    (35)
#define ARGUMENT_BEAM           (36)
#define ARGUMENT_STATE_BEAM     (37)
#define ARGUMENT_END_STATE_BEAM (38)
#define ARGUMENT_LMLA           (39)
#define ARGUMENT_HISTOGRAM_STATE_PRUNING (40)
#define ARGUMENT_HISTOGRAM_PRUNING       (41)
#define ARGUMENT_HYPFILE                 (42)
#define ARGUMENT_AM_PHONE                (43)
#define ARGUMENT_LATTICETRAINING         (44)
/////////////////////////////////////// shout_maketrainset
// #define ARGUMENT_AUDIOFILE
// #define ARGUMENT_HYPFILE
// #define ARGUMENT_TRAINDIR
// #define ARGUMENT_META_IN
#define ARGUMENT_PHONE_LIST     (45)
#define ARGUMENT_AM_TYPE        (46)
#define ARGUMENT_AUDIO_LIST     (47)
#define ARGUMENT_SPEAKER_FILTER (48)

/////////////////////////////////////// shout_adapt_am
// #define ARGUMENT_TRAINDIR
// #define ARGUMENT_AM_PHONE
#define ARGUMENT_AM_OUT         (50)
#define ARGUMENT_AM_CLUSTER     (51)


/////////////////////////////////////// shout_dct2lextree
// #define ARGUMENT_DCT_BIN
#define ARGUMENT_DCT_TEXT       (52)

/////////////////////////////////////// shout_lm2bin
// #define ARGUMENT_DCT_BIN
// #define ARGUMENT_LM_BIN
// #define ARGUMENT_PHONE_LIST
#define ARGUMENT_LM_ARPA        (53)

/////////////////////////////////////// shout_train_model
// #define ARGUMENT_TRAINDIR
#define ARGUMENT_MAXGAUSSIANS   (54)
#define ARGUMENT_TRAINSILP      (55)
#define ARGUMENT_TRAINSILMAX    (56)
#define ARGUMENT_TRAINMODELNR   (57)
#define ARGUMENT_AM_CREATE      (58)
#define ARGUMENT_DECISIONTREE   (59)
#define ARGUMENT_DOFAST         (60)
// ARGUMENT_LATTICETRAINING

/////////////////////////////////////// shout_train_finish
// #define ARGUMENT_TRAINDIR
// #define ARGUMENT_AM_CREATE

/////////////////////////////////////// shout_prepare_adapt
// For now not in use...

/////////////////////////////////////// shout_merge_am
// #define ARGUMENT_AM_OUT
#define ARGUMENT_AM_MERGE1      (61)
#define ARGUMENT_AM_MERGE2      (62)

/////////////////////////////////////// shout_mergetrainset
#define ARGUMENT_TRAINDIR2      (63)

/////////////////////////////////////// shout_train_finish_sat
//ARGUMENT_TRAINDIR
//ARGUMENT_AM_CREATE
//ARGUMENT_TRAINMODELNR
#define ARGUMENT_SATLIST        (64)


/////////////////////////////////////// shout_expand_vocab
//#define ARGUMENT_DCT_BIN
//#define ARGUMENT_BACKDCT_BIN
// #define ARGUMENT_HYPFILE

/////////////////////////////////////// shout_train_mmi
#define ARGUMENT_LATTICETRAINING_ENUM    (65)
#define ARGUMENT_LATTICETRAINING_DENOM   (66)

/////////////////////////////////////// shout_online
#define ARGUMENT_CONFIGDIR               (67)

////////////////////////////////////// shout_preprocess

//#define ARGUMENT_PHONE_LIST
//#define ARGUMENT_AUDIO_LIST
//#define ARGUMENT_AM_PHONE
//#define ARGUMENT_TRAINDIR
#define ARGUMENT_NOT_USED                (49)


/////////////////////////////////////// shout_spkrec
#define ARGUMENT_SPKREC_UBM_MALE         (68)
#define ARGUMENT_SPKREC_UBM_FEMALE       (69)
#define ARGUMENT_SPKREC_TASKLIST         (70)
#define ARGUMENT_SPKREC_MODELDIR         (71)
#define ARGUMENT_SPKREC_OUTPUT_PREFIX    (72)
#define ARGUMENT_SPKREC_ZNORM_LIST       (73)
#define ARGUMENT_SPKREC_OUTPUT_STATS     (74)

#define ARGUMENTS_NUMBER                 (75)


/////////////////////////////////////////////////////////////////////////////////////////
/// \brief This structure contains all applications that may register to ShoutConfig
/////////////////////////////////////////////////////////////////////////////////////////

// Currently, only one file type is supported:
typedef enum
{
  APPID_SEGMENT,
  APPID_CLUSTER,
  APPID_VTLN,
  APPID_SHOUT,
  APPID_MAKETRAINSET,
  APPID_ADAPT_AM,
  APPID_DCT2LEXTREE,
  APPID_LM2BIN,
  APPID_TRAIN_MMI,
  APPID_TRAIN_MODEL,
  APPID_TRAIN_FINISH,
  APPID_MERGE_AM,
  APPID_PREPROCESS,
  APPID_NORMALIZE_AM,
  APPID_PREPARE_ADAPT,
  APPID_MERGETRAINSET,
  APPID_TRAIN_FINISH_SAT,
  APPID_AVG_ENERGY,
  APPID_EXPAND_VOCAB,
  APPID_UPDATE_VERSION,
  APPID_SPKREC,
  APPID_SPKRECSTATS,
  APPID_SHOUTONLINE
} APPLICATION_ID;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct
{
  char value[MAXSTR];
  char description[MAXSTR];
  char descrVal[MAXSTR];
} CommandParameterType;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

class ShoutConfig
{
  protected:
    char **parShort;
    char **parLong;
    CommandParameterType *parameter;
    int                   nrPars;
    APPLICATION_ID        applicationID;

  public:
    ShoutConfig(APPLICATION_ID appID,int argc,char *argv[]);
   ~ShoutConfig();

  public:
    bool  parIsDefined(int parameter);
    char *getStringValue(int parameter);
    float getFloatValue(int parameter);

  protected:
    void initializeParameter(int id, const char *shortDiscr, const char *longDiscr, const char *text, const char *valDiscr, double val);
    void initializeParameter(int id, const char *shortDiscr, const char *longDiscr, const char *text, const char *valDiscr, const char *val);
    int  parInApp(int i);
    void readConfigFile(FILE *config, StringLookup *lookupLong);
    bool fillParameter(char *line, StringLookup *lookupLong);

};

#endif
