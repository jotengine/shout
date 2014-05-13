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
#include "io_speaker_segmenter.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor will read speaker segmentation information from inFile.
/////////////////////////////////////////////////////////////////////////////////////////////////////

IO_Speaker_Segmenter::IO_Speaker_Segmenter(FILE *inFile)
{
  numberOfClusters = 0;
  for(int i=0;i<INITIAL_NUMBER_OF_CLUSTERS;i++)
  {
    trainCluster[i] = NULL;
  }
  if(inFile != NULL)
  {
    fread(&numberOfClusters,1,sizeof(numberOfClusters),inFile);
    for(int i=0;i<numberOfClusters;i++)
    {
      trainCluster[i] = new TrainPhoneModel(inFile,FEATURE_SPEAKER_SIZE);
    }  
    fclose(inFile);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor clears allocated data.
/////////////////////////////////////////////////////////////////////////////////////////////////////


IO_Speaker_Segmenter::~IO_Speaker_Segmenter()
{
  for(int i=0;i<INITIAL_NUMBER_OF_CLUSTERS;i++)
  {
    if(trainCluster[i] != NULL)
    {
      delete trainCluster[i];
    }
  }    

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method stores the speaker clustering information to disc (loaded by the constructor).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void IO_Speaker_Segmenter::storeClusters(FILE *outFile)
{
  fwrite(&numberOfClusters,1,sizeof(numberOfClusters),outFile);
  for(int i=0;i<INITIAL_NUMBER_OF_CLUSTERS;i++)
  {
    if(trainCluster[i] != NULL)
    { 
      trainCluster[i]->writeModel(outFile);
    }
  }  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method returns the number of clusters. If no clusters are available, one is returned.
/// (all data is assigned the same clusterID)
/////////////////////////////////////////////////////////////////////////////////////////////////////

int IO_Speaker_Segmenter::getNumberOfClusters()
{
  if(numberOfClusters == 0)
  {
    return 1;
  }
  return numberOfClusters;  
}
