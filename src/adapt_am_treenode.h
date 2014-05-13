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
#ifndef ADAPT_AM_TREENODE_H
#define ADAPT_AM_TREENODE_H

#include "gaussian.h"
#include "vector.h"


class PhoneModel;
class Gaussian;


//////////////////////////////////////////////////////////////////////////////////////////////
/// \brief This class represents a single cluster in the SMAPLR adaptation method.
/// method. Each cluster is structured in a node of a tree, hence the name Adapt_AM_TreeNode. 
/// Each node has access to a number of Gaussians of the acoustic models. The class creates
/// an adaptation matrix using MAPLR adaptation and passes its matrix as M-matrix to its children.
///
/// The description of the general order of method-calls can be found in the description of the
/// Adapt_AM class. 
/////////////////////////////////////////////////////////////////////////////////////////////

class Adapt_AM_TreeNode
{
  protected:
    double C_factor;
    int    dimensions;
  public:
    double W[ADAPT_DEFAULT_VECTORSIZE][ADAPT_DEFAULT_VECTORSIZE+1];                              ///< The adaptation matrix
    double M[ADAPT_DEFAULT_VECTORSIZE][ADAPT_DEFAULT_VECTORSIZE+1];                              ///< Prior matrix (I at top, otherwise the W matrix of parent)
    double Z[ADAPT_DEFAULT_VECTORSIZE][ADAPT_DEFAULT_VECTORSIZE+1];                              ///< Helper matrix for calculating W.
    double AUGMENTED_ALL[ADAPT_DEFAULT_VECTORSIZE*(ADAPT_DEFAULT_VECTORSIZE+1)*(ADAPT_DEFAULT_VECTORSIZE+2)];  ///< All data from the member-gaussians.
    Gaussian *meanGaussian;                                          ///< The mean Gaussian is used for splitting the cluster.
    Adapt_AM_TreeNode *parent;                                       ///< The parent will deliver the M matrix
    Adapt_AM_TreeNode *childLeft;       ///< The gaussians of the children are also member of this cluster. 
    Adapt_AM_TreeNode *childRight;      ///< The gaussians of the children are also member of this cluster. 
    Vector *varTrans;
    double  varDenom;
    double  denominator;
  protected:
    void    iterate_setHelperMatrices_1    ();
  public:
    Adapt_AM_TreeNode(int dim);
   ~Adapt_AM_TreeNode();
      
    int     split                          (int depth);  
    void    setHelperMatrices_1            (PhoneModel **models, int numberOfModels, bool doSil);
    void    setHelperMatrices_1            (PhoneModel **models, int numberOfModels);
    void    setHelperMatrices_2            ();
};

#endif
