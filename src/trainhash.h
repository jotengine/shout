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
#ifndef HASHTRAINHASH_H
#define HASHTRAINHASH_H

#include "hash.h"


struct VertexSet;
struct Edge;


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief This structure can store a set of vertices. 
///
/// A vertex is a node in a graph. Each vertex from the
/// VertexSet contains a list of edges stored in Edge and an ID of the set the vertex is in. This ID 
/// makes it possible to check very quickly if the graph is still a-cyclic.
///
/// VertexSets are connected by next and previous pointers.
/////////////////////////////////////////////////////////////////////////////////////////////////////

struct VertexSet
{
  int        setNr;
  int        next;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief This structure stores the edges (connections) of the graph.
/////////////////////////////////////////////////////////////////////////////////////////////////////

struct Edge
{
  int   vertexID_A;
  int   vertexID_B;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief This class handles creation of new Hash functions.
///
/// The theory from: ["An Optimal Algorithm for Generating Minimal Perfect Hash Functions", 
/// Z.J. Czech, G. Havas and B.S. Majewski]. is used to construct and use the hash function.
/// This class only loads, stores and uses the hash function. The class TrainHash is able
/// to create new hash functions. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

class TrainHash : public Hash
{
  protected:
    Edge          *edges;
    VertexSet     *vertexSet;
    int            nrOfEdges;

  public:
    TrainHash                      (int nrOfS, int d, int nrOfW);
   ~TrainHash                      ();

   protected:
     bool traverse                 ();
   public:
     void storeHash                (FILE *file);
     void initialiseMapping        ();
     bool fillMapping              (int *w);
     void finalizeHash             ();
     int  getIndex                 (int *w);
};

#endif
