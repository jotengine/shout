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
#ifndef SPEAKERRECOGNITION_H
#define SPEAKERRECOGNITION_H

#define MAX_UBMSIZE (200)

#include "phonemodel.h"
#include "gaussian.h"

class SpeakerRecognition
{
  private:
    int                      vectorSize;
    PhoneModel              *ubmMale[MAX_UBMSIZE];
    PhoneModel              *ubmFemale[MAX_UBMSIZE];
    double                   ubmZNormMean[2][MAX_UBMSIZE];
    double                   ubmZNormVariance[2][MAX_UBMSIZE];
    bool                     ubmZNormAvailable[2];
    SegmentationAdmin        readPool_curSeg;
    bool                     zNormAvailable;
    double                   zNormMean;
    double                   zNormVariance;
    int                      nrUbmMale;
    int                      nrUbmFemale;
    bool                     outputStats;

  public:
    SpeakerRecognition(char *ubmMaleFile, char *ubmFemaleFile, char *tasklist, 
                       char *modelDir, char *outputPrefix, char *znormList, bool outputStats);
   ~SpeakerRecognition();

  private:
    void        readUbmModel(char *modelName, bool male);
    PhoneModel *readModel(char *modelName, bool applyZNorm);
    void        writeModel(char *modelName, PhoneModel *model);
    Gaussian *runTrials(char *tasklist, char *modelDir, char *outputPrefix, PhoneModel *zNormSpeakerModel = NULL);
};

#endif
