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
#include "normalizeam.h"
#include "trainphonemodel.h"
#include "shoutconfig.h"

#include <stdlib.h>
#include <string.h>
#include "featurepool.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Handles the command parameters for Shout_MakeTrainSet. It is the main function of the 
/// shout_maketrainset application. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_NORMALIZE_AM, argc, argv);

  NormalizeAM NormalizeAM(myConfig.getStringValue(ARGUMENT_AUDIOFILE), myConfig.getStringValue(ARGUMENT_META_IN));
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Doc 
/////////////////////////////////////////////////////////////////////////////////////////////////////

NormalizeAM::NormalizeAM(char *inFile, char *rttmName)
{
  FeaturePool *featurePool = new FeaturePool(inFile,FILE_TYPE_RAW_DECODING,false,1.0,false,true);  // Only prepare...
  int segID = featurePool->claimSegmentationID();


  char segmentString[60];
  strcpy(segmentString,"SPEECH");
  strcpy(&(segmentString[20]),"SIL");
  strcpy(&(segmentString[40]),"SOUND");
  FILE *inSegFile = fopen(rttmName,"rt");
  if(inSegFile == NULL)
  {
    USER_ERROR("ERROR: Could not read the input RTTM file for decoding!\n");
  }
  featurePool->readRTTM(segID,segmentString,3,20,0,inSegFile);
  fclose(inSegFile);

  featurePool->createNewPool(inFile, NULL, NULL, false, false, segID);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Doc 
/////////////////////////////////////////////////////////////////////////////////////////////////////

NormalizeAM::~NormalizeAM()
{
}


