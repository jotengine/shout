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
#include "adapt_am.h"
#include "shoutconfig.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Handles the command parameters for AdaptAM. It is the main function of the 
/// shout_adapt_am application. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  ShoutConfig myConfig(APPID_ADAPT_AM, argc, argv);

  Adapt_AM Adapt_AM(myConfig.getStringValue(ARGUMENT_TRAINDIR),myConfig.getStringValue(ARGUMENT_AM_PHONE),
                    myConfig.getStringValue(ARGUMENT_AM_CLUSTER),myConfig.getStringValue(ARGUMENT_AM_OUT));
  return 0;
}

