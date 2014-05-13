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
#ifndef SHOUT_MAKETRAINSET_H
#define SHOUT_MAKETRAINSET_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vector.h"
#include "phonemodel.h"
#include "featurepool.h"

typedef enum
{
  TASK_TYPE_PHONE             = 0,
  TASK_TYPE_SAD               = 1,
  TASK_TYPE_HIST_NORM         = 2,
  TASK_TYPE_VTLN              = 3,
  TASK_TYPE_ARTICULATORY      = 4,
  TASK_TYPE_CLUSTER           = 5,
  TASK_TYPE_DIAR              = 6,
  TASK_TYPE_UNDEF             = 7
} TASK_TYPE;

////////////////////////////////////////////////////////////////////////////
/// \brief This class is the base class for the shout_maketrainset application.
/// It creates a new AM trainset given a Master Label File and an output directory.
////////////////////////////////////////////////////////////////////////////

class Shout_MakeTrainSet
{
  protected:
    bool                performCluster;
    TASK_TYPE           performWhat;
    FILE              **phoneFile;
    DataStats          *data;
    int                 nrOfModels;
    int                 nrOfDecisionRules;
    int                 nrOfDecisionRuleSamples[MAX_NUMBER_OF_DECISION_RULES][MAX_CLUSTERS];
    int                 nrOfSil;
    char               *trainDirName;
    FeaturePool        *featurePool;
    int                 totalPhonesStored;
    int                 vectorSize;
    int                 skippedLines;
    float               useProbability;

  protected:
    void openPhoneFiles();
    void checkDataDir(bool allowCDM);
  public:
    Shout_MakeTrainSet(char *task, char *phoneSetFileName, char *trainDirName);
    Shout_MakeTrainSet(char *task, char *phoneSetFile, char *datFile, char *rawName, char *rttmFileName, char *trainDirName, 
                       char *decisionTree, char *method, char *trainList, char *speakerFilter, float useP);

   ~Shout_MakeTrainSet();

    void closePhoneFiles(bool printOutput);
    void setFeaturePool(FeaturePool *fp);
    void storePhone(int cdmID, int contextLeft, int phoneID, int contextRight, int b1, int b2, int b3, int numberOfObservations);
};

#endif
