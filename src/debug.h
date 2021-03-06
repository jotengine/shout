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

#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>
#include <iomanip>

#define DEBUG_OFF

#ifdef DEBUG_ON
  #define DEBUG_PRINT(a) std::cout << a << std::endl
  #define DEBUG_DIE(a) std::cout << a << std::endl; exit(1);
  #define DEBUG_DIE_IF(a,b) if(a){ std::cout << b << std::endl; exit(1); }
#else
  #define DEBUG_PRINT(a) ;
  #define DEBUG_DIE(a) ;
  #define DEBUG_DIE_IF(a,b) ;
#endif


#endif // DEBUG_H
