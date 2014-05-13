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

#ifndef SHOUT_MISC_H
#define SHOUT_MISC_H

#define FASTLOG_2_POW_20  (1048576)
#define FASTLOG_N         (20)
#define FASTLOG_23_MIN_N  (3)

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Some methods for reading and writing to file, big/little-endian safe.
/////////////////////////////////////////////////////////////////////////////////////////////////////

namespace WriteFileLittleBigEndian
{
  int fwriteEndianSafe(void *var, int nrTimes, int nrBytes, FILE* file);
  int freadEndianSafe(void *var,  int nrTimes, int nrBytes, FILE* file);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Some methods that operate on strings
/////////////////////////////////////////////////////////////////////////////////////////////////////

namespace StringFunctions
{
  void RightTrim(char *str);
  void splitList(char *str, char *str2);
  void splitRef(char *str, char *str2);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Fast implementation of log by Laurent de Soras, adjusted to table-lookup by people at ICSI
/////////////////////////////////////////////////////////////////////////////////////////////////////

namespace FastMath
{

  static float logLookupTable[FASTLOG_2_POW_20];

  void initLog();

  inline float log(float val)
  {
    int * const  exp_ptr = reinterpret_cast <int *> (&val);
    int          x = *exp_ptr;
    const int    log_2 = ((x >> 23) & 255) - 128;
    x &= ~(255 << 23);
    x += 127 << 23;
    *exp_ptr = x;

    return (val + log_2);

/*    int *const     exp_ptr = ((int*)&val);
    int            x = *exp_ptr;
    const int      log_2 = ((x >> 23) & 255) - 127;
    x &= 0x7FFFFF;
    x = x >> (FASTLOG_23_MIN_N);
    val = logLookupTable[x];
    return ((val + log_2)* 0.69314718);*/
  }
}


#endif
