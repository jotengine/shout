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
#include <stdlib.h>
#include "hash.h"
#include "shout-misc.h"

using namespace WriteFileLittleBigEndian;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor will load the hash function.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Hash::Hash(FILE *file, char *memoryPool, long int *offset)
{
  random1 = NULL;
  random2 = NULL;

  doMemoryMapping = (DO_MEMORY_MAPPED_FILE_LM && (memoryPool != NULL));
  if(doMemoryMapping && offset != NULL)
  {
    randomTableLength = *((int*)(memoryPool+*offset));
    *offset += sizeof(int);
    nrOfSamples = *((int*)(memoryPool+*offset));
    *offset += sizeof(int);
    gLength = *((int*)(memoryPool+*offset));
    *offset += sizeof(int);
    if(int(HASH_FACTOR_VERSION*nrOfSamples) != gLength)
    {
      PRINTF4("%f %d %d\n",HASH_FACTOR_VERSION,nrOfSamples,gLength);
      USER_ERROR("The language model binary file is invalid. You might be using an old version. Please rebuild the LM with shout_lm2bin.");
    }

    depth = *((int*)(memoryPool+*offset));
    *offset += sizeof(int);

    random1 = ((int*)(memoryPool+*offset));
    *offset += depth*randomTableLength*sizeof(int);
    random2 = ((int*)(memoryPool+*offset));
    *offset += depth*randomTableLength*sizeof(int);

    gTable = ((int*)(memoryPool+*offset));
    *offset += gLength*sizeof(int);
  }
  else if(file!=NULL)
  {
    freadEndianSafe(&randomTableLength,1,sizeof(randomTableLength),file);
    freadEndianSafe(&nrOfSamples,1,sizeof(nrOfSamples),file);
    freadEndianSafe(&gLength,1,sizeof(gLength),file);
// A check on LM version to make sure the version improvement isn't going to backfire on us:
    if(int(HASH_FACTOR_VERSION*nrOfSamples) != gLength)
    {
      USER_ERROR("The language model binary file is invalid. You might be using an old version. Please rebuild the LM with shout_lm2bin.");
    }
    freadEndianSafe(&depth,1,sizeof(depth),file);
    random1 = new int[depth*randomTableLength];
    random2 = new int[depth*randomTableLength];
    int pointer = 0;
    for(int i=0;i<depth;i++)
    {
      freadEndianSafe(&random1[pointer],randomTableLength,sizeof(int),file);
      pointer += randomTableLength;
    }
    pointer = 0;
    for(int i=0;i<depth;i++)
    {
      freadEndianSafe(&random2[pointer],randomTableLength,sizeof(int),file);
      pointer += randomTableLength;
    }
    gTable = new int[gLength]; 
    freadEndianSafe(gTable,gLength,sizeof(int),file);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor deletes the data needed for the hash function.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Hash::~Hash()
{
  if(!doMemoryMapping)
  {
    if(gTable != NULL)
    {
      delete[] gTable;
    }
    if(random1 != NULL)
    {
      delete[] random1;
      delete[] random2;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Query for a data slot index given a key of depth 'depth' (w is a pointer to int[depth]). 
/// The user is responsible for checking if the returned index is valid. (Because this is a minimal 
/// perfect hash, all slots will contain data. 
/// If a non-existing key is queried, it will be randomly mapped to an existing slot.)
/////////////////////////////////////////////////////////////////////////////////////////////////////

int Hash::getIndex(int *w)
{
  int result;
  int u = 0;
  int v = 0;
  int pointer = 0;
  for(int i=0;i<depth;i++)
  { 
    int ww = w[i];
    u += random1[pointer+ww];
    v += random2[pointer+ww];
    pointer += randomTableLength;
  }
  MODULO(u, u, gLength);
  MODULO(v, v, gLength);  
  MODULO(result,(gTable[u]+gTable[v]),nrOfSamples);
  return result;
}


