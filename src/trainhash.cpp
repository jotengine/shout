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

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "standard.h"
#include "trainhash.h"
#include "shout-misc.h"

using namespace WriteFileLittleBigEndian;



/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The constructor initialises the training parameters.
/////////////////////////////////////////////////////////////////////////////////////////////////////

TrainHash::TrainHash(int nrOfS, int d, int nrOfW) : Hash(NULL, NULL,NULL)
{
printf("\nInitialize\n");
  nrOfEdges         = 0;
  randomTableLength = nrOfW;
  nrOfSamples       = nrOfS;
  gLength           = int(HASH_FACTOR_VERSION*nrOfSamples);
printf("Hash: %d * %f = %d\n",nrOfSamples,HASH_FACTOR_VERSION,gLength);

  edges             = new Edge[nrOfSamples];
  gTable            = new int[gLength];
  vertexSet         = new VertexSet[gLength];

  for(int i=0;i<gLength;i++)
  {
    gTable[i]           = 0;
    vertexSet[i].setNr  = i;
    vertexSet[i].next   = -1;
  }
  depth   = d;
  random1 = new int[depth*randomTableLength];
  random2 = new int[depth*randomTableLength];
  srand((unsigned)time(0));
printf("Go!\n");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor deletes all claimed memory.
/////////////////////////////////////////////////////////////////////////////////////////////////////

TrainHash::~TrainHash()
{
  delete[] edges;
  delete[] vertexSet;
  // gTable and random are deleted in ~Hash.
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method deletes all existing data and will generate new random numbers for a new mapping.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainHash::initialiseMapping()
{
printf("another init\n");
  for(int i=0;i<randomTableLength;i++)
  {
    for(int i2 = 0;i2<depth;i2++)
    {
      random1[i2*randomTableLength+i] = int(double(gLength)*rand()/(RAND_MAX));
      random2[i2*randomTableLength+i] = int(double(gLength)*rand()/(RAND_MAX));
    }
  }
  nrOfEdges    = 0;
  for(int i=0;i<gLength;i++)
  {
    gTable[i]           = 0;
    vertexSet[i].setNr  = i;
    vertexSet[i].next   = -1;
  }
printf("Go2!\n");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Given the by initialiseMapping() random generated numbers, the mapping is filled. If the mapping
/// is not a-cyclic (and therefor it is useless), this method will return false. If a succesfull
/// mapping was possible, it will return true.
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool TrainHash::fillMapping(int *w)
{
  assert(nrOfEdges >= 0 && nrOfEdges < nrOfSamples);
  int u = 0;
  int v = 0;
  for(int i=0;i<depth;i++)
  {
    u += random1[i*randomTableLength+w[i]];
    v += random2[i*randomTableLength+w[i]];
  }
  MODULO(u, u, gLength);
  MODULO(v, v, gLength);

  gTable[u]      = 0;
  gTable[v]      = 0;

  int theNumber = vertexSet[u].setNr;
  if(vertexSet[v].setNr < theNumber)
  {
    int aa = v;
    v = u;
    u = aa;
    theNumber = vertexSet[u].setNr;
  }

  if(theNumber == vertexSet[v].setNr)
  {
    return false;
  }
  else
  {
    int nr = u;
    while(vertexSet[nr].next >= 0)
    {
      nr = vertexSet[nr].next;
    }
    vertexSet[nr].next = vertexSet[v].setNr;
    while(vertexSet[nr].next >= 0)
    {
      nr = vertexSet[nr].next;
      vertexSet[nr].setNr = theNumber;
    }
    // Update data structure:
    edges[nrOfEdges].vertexID_A = v;
    edges[nrOfEdges].vertexID_B = u;
    nrOfEdges++;
  }
  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method may be used when fillMapping() has returned true for all samples. It will fill
/// the gTable table with the hash appropriate hash indices.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainHash::finalizeHash()
{
printf("Finalize\n");

  // If the graph is acyclic (is tested earlier), do the traverse step:
  for(int i=0;i<gLength;i++)
  {
    if(vertexSet[i].setNr == i && vertexSet[i].next >= 0)
    {
      vertexSet[i].next = 1;
    }
    else
    {
      vertexSet[i].next = 0;
    }
    gTable[i]      = 0;
  }
printf("Traverse\n");
  while(traverse())
  {
  }
printf("Done!\n");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This is a helper function for finalizeHash().
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool TrainHash::traverse()
{
  bool didSomething = false;
  for(int i=0;i<nrOfEdges;i++)
  {
    int vertexID_A = edges[i].vertexID_A;
    if(vertexID_A >= 0)
    {
      int vertexID_B    = edges[i].vertexID_B;
      bool prVertexID_A = (vertexSet[vertexID_A].next == 1);
      bool prVertexID_B = (vertexSet[vertexID_B].next == 1);
      if(prVertexID_A && !prVertexID_B)
      {
        MODULO(gTable[vertexID_B], (i - gTable[vertexID_A]), nrOfSamples);
        didSomething = true;
        vertexSet[vertexID_B].next = 1;
        edges[i].vertexID_A = -1;
      }
      if(prVertexID_B && !prVertexID_A)
      {
        MODULO(gTable[vertexID_A], (i - gTable[vertexID_B]), nrOfSamples);
        didSomething = true;
        vertexSet[vertexID_A].next = 1;
        edges[i].vertexID_A = -1;
      }
    }
  }
  return didSomething;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void TrainHash::storeHash(FILE *file)
{
  fwriteEndianSafe(&randomTableLength,1,sizeof(randomTableLength),file);
  fwriteEndianSafe(&nrOfSamples,1,sizeof(nrOfSamples),file);
  fwriteEndianSafe(&gLength,1,sizeof(gLength),file);
  fwriteEndianSafe(&depth,1,sizeof(depth),file);
  for(int i=0;i<depth;i++)
  {
    fwriteEndianSafe(&random1[i*randomTableLength],randomTableLength,sizeof(int),file);
  }
  for(int i=0;i<depth;i++)
  {
    fwriteEndianSafe(&random2[i*randomTableLength],randomTableLength,sizeof(int),file);
  }
  for(int i=0;i<gLength;i++)
  {
    fwriteEndianSafe(&gTable[i],1,sizeof(int),file);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

int TrainHash::getIndex(int *w)
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


