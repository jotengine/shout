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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "standard.h"
#include "spkrecstats.h"
#include "shoutconfig.h"
#include "trainphonemodel.h"
#include "shout-misc.h"
#include "featurepool.h"

using namespace WriteFileLittleBigEndian;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Handles the command parameters for this application.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_SPKRECSTATS, argc, argv);

  SpkrecStats SpkrecStats(myConfig.getStringValue(ARGUMENT_SPKREC_UBM_MALE),
                          myConfig.getStringValue(ARGUMENT_AUDIOFILE),
                          myConfig.getStringValue(ARGUMENT_META_IN),
                          myConfig.parIsDefined(ARGUMENT_SPKREC_OUTPUT_STATS));
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

SpkrecStats::SpkrecStats(char *ubmFile, char *audioFile, char *metaInFile, bool outputStats)
{
  ubm = NULL;
  if(strlen(ubmFile)>0)
  {
    readUbmModel(ubmFile);
  
    if(outputStats)
    {
      FILE *outMeans     = fopen("ubm_means","wt");
      FILE *outVariances = fopen("ubm_variances","wt");
      FILE *outWeights   = fopen("ubm_weights","wt");
      ubm->printModel(outMeans,outVariances,outWeights);
      fclose(outMeans);
      fclose(outVariances);
      fclose(outWeights);
    }
    else
    {
      FeaturePool *featurePool  = NULL;
      printf("Initializing feature pool...\n");
      featurePool = new FeaturePool(audioFile,FILE_TYPE_RAW_SPKREC,false,1.0,false,true);  // Only prepare
      int segID = featurePool->claimSegmentationID();
      printf("Reading RTTM meta input file...\n");
      FILE *metaIn = fopen(metaInFile, "rt");
      if(metaIn == NULL)
      {
        USER_ERROR("Could not read the input RTTM file.");
      }
      featurePool->readRTTM(segID,NULL,-1,-1,3,metaIn);
      fclose(metaIn);
      printf("Generating features...\n");
      featurePool->createNewPool(audioFile, NULL, NULL, false, false, segID, 1.0);  // NO VTLN!
      assert(vectorSize == featurePool->getVectorSize());
  
      printf("Calculating statistics...\n");
  
      for(int spk=0;spk<featurePool->getNumberOfSpeakers(segID);spk++)
      {
        int spkID = featurePool->getSpeakerID(segID,spk);
        printf("  ...for speaker SPK%02d.\n",spkID);
        char outFileName[MAXSTR];
        sprintf(outFileName,"SPK%02d.N.stats",spkID);
        FILE *outFileN = fopen(outFileName,"wb");
        sprintf(outFileName,"SPK%02d.F.stats",spkID);
        FILE *outFileF = fopen(outFileName,"wb");
        assert(outFileN != NULL);
        assert(outFileF != NULL);
        ubm->adapt_clear();
        ubm->adapt_setAcumulators(spkID, segID, featurePool);
        ubm->writeAccumulators(outFileN, outFileF,NULL,true);
        fclose(outFileN);
        fclose(outFileF);
      }
      printf("Done. Bye-bye!\n");
    }
  }
  else  // Just output features.
  {
    FeaturePool *featurePool  = NULL;
    printf("Initializing feature pool...\n");
    featurePool = new FeaturePool(audioFile,FILE_TYPE_RAW_SPKREC,false,1.0,false,false);
    featurePool->storeFeatureVectors(stdout);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

SpkrecStats::~SpkrecStats()
{
  if(ubm != NULL)
  {
    delete ubm;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void SpkrecStats::readUbmModel(char *modelName)
{
  PhoneModel *res = NULL;
  FILE *modelFile = fopen(modelName, "rb");

  if(modelFile == NULL)
  {
    USER_ERROR("Could not read the model file.");
  }

  int nrM, nrC;
  int r = freadEndianSafe(&vectorSize,1,sizeof(vectorSize),modelFile);
  freadEndianSafe(&nrM,1,sizeof(nrM),modelFile);
  freadEndianSafe(&nrC,1,sizeof(nrC),modelFile);
  if(nrM != 1)
  {
    fclose(modelFile);
    USER_ERROR("The file is not a model file.");
  }
  ubm = new PhoneModel(modelFile,vectorSize);
  fclose(modelFile);
}

