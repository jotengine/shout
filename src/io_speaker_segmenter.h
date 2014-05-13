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
#ifndef IO_SPEAKER_SEGMENTER_H
#define IO_SPEAKER_SEGMENTER_H

#include "trainphonemodel.h"

////////////////////////////////////////////////////////////////////////////
/// \brief The AM adaptation application does not need any speaker clustering
/// information. The prepare adaptation application handles clustering.
/// Unfortunately the adaptation application needs to read and write the 
/// speaker clustering information from and to disc. For this purpose this 
/// small class is created with only read and write methods.  
////////////////////////////////////////////////////////////////////////////

class IO_Speaker_Segmenter 
{
  
  protected:  
    int               numberOfClusters;
    TrainPhoneModel  *trainCluster[INITIAL_NUMBER_OF_CLUSTERS];
    
  public:
    IO_Speaker_Segmenter                   (FILE *inFile);
   ~IO_Speaker_Segmenter                   ();

    int    getNumberOfClusters             ();
    void   storeClusters                   (FILE *outFile);
   
};

#endif
