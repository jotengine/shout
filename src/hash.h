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
#ifndef HASH_H
#define HASH_H

#include "standard.h"
#include <stdio.h>


#define HASH_FACTOR_VERSION  (3.0)


////////////////////////////////////////////////////////////////////////////
/// \brief The Hash class handles minimal perfect hashing.
///
/// The theory from: ["An Optimal Algorithm for Generating Minimal Perfect Hash Functions", 
/// Z.J. Czech, G. Havas and B.S. Majewski]. is used to construct and use the hash function.
/// This class only loads, stores and uses the hash function. The class TrainHash is able
/// to create new hash functions. 
///
/// The implementation of Hash is fine tuned for use in the shout software package. 
/// It is possible to query for data slots with as key either a bigram, a trigram 
/// or 4-gram. (getIndex with *w different length, depending on how it is trained).
////////////////////////////////////////////////////////////////////////////

class Hash
{
  public:
    int    randomTableLength;
    int    nrOfSamples;
    int    gLength;
    int   *random1, *random2;
    int   *gTable;
    int    depth;
    bool   doMemoryMapping;

  public:
    Hash(FILE *file, char *memoryPool, long int *offset);
   ~Hash();

  public:   
    int  getIndex                 (int *w);   
};

#endif
