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
#ifndef SHOUTPREPAREADAPT_H
#define SHOUTPREPAREADAPT_H

#include "standard.h"
#include "phonemodel.h"

////////////////////////////////////////////////////////////////////////////
/// \brief This class is the main class of shout_prepare_adapt.
////////////////////////////////////////////////////////////////////////////

class ShoutPrepareAdapt
{
  protected:
    void adaptingModel(PhoneModel *model);
  public:
    ShoutPrepareAdapt(FILE *amFile, char *dblFile, int maxClusters, char *outDblName);
   ~ShoutPrepareAdapt();
};

#endif
