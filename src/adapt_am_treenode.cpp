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

#include "standard.h"
#include "adapt_am_treenode.h"
#include "phonemodel.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Code efficiency:
extern int AUGMENTED_MATRIX_SIZE;
extern int AUG_X;
extern int AUG_Y;
extern int AUG_Y_TIMES_VECTOR;
extern int AUG_P_PLUS_ONE;
extern int AUG_P;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The constructor checks the corectness of some pre-compiler settings (defines) and initialises
/// the augmentation matrix for this tree node. The child-nodes are NULL.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Adapt_AM_TreeNode::Adapt_AM_TreeNode(int dim)
{  
  assert(dim <= ADAPT_DEFAULT_VECTORSIZE);
  
  dimensions               = dim;
  AUGMENTED_MATRIX_SIZE    = dimensions*(dimensions+1)*(dimensions+2);
  AUG_X                    = (dimensions+1)*(dimensions+2);
  AUG_Y                    = (dimensions+2);
  AUG_Y_TIMES_VECTOR       = dimensions*AUG_Y;
  AUG_P_PLUS_ONE           = dimensions*AUG_Y+(dimensions+1);
  AUG_P                    = dimensions*AUG_Y+dimensions;

  denominator  = 0.0;
  varDenom     = 0.0;
  varTrans     = new Vector(dim);
  varTrans->setAllValues(0.0);
  meanGaussian = NULL;
  parent       = NULL;
  childLeft    = NULL;
  childRight   = NULL;
  C_factor     = 0.01;
  
  
  for(int k=0;k<AUGMENTED_MATRIX_SIZE;k++)
  {
    AUGMENTED_ALL[k] = 0.0;
  }
  for(int k=0;k<dimensions;k++)
  {
    for(int l=0;l<dimensions+1;l++)
    {
      W[k][l] = 0.0;
      M[k][l] = 0.0;      
    }  
    M[k][k] = 1.0;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor deletes its children and its mean-Gaussian vector.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Adapt_AM_TreeNode::~Adapt_AM_TreeNode()
{
  if(varTrans != NULL)
  {
    delete varTrans;
    varTrans = NULL;
  }
  if(meanGaussian != NULL)
  {
    delete meanGaussian;
    meanGaussian = NULL;
  }
  if(childLeft != NULL)
  {
    delete childLeft;
    childLeft = NULL;
  }
  if(childRight != NULL)
  {
    delete childRight;
    childRight = NULL;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The first helper matrix contains the statistics of all acoustic models that are a member
/// of this tree node (cluster). This method calls the PhoneModel::adapt_setHelperMatrices() method
/// of all acoustic models. This method will call MixGaussian::MixGaussian() for each mixed gaussian
/// (triphone model) and this method will call Gaussian::adapt_setHelperMatrices() for each single
/// Gaussian in the model. Each Gaussian will then add its statistics to its own Adapt_AM_TreeNode
/// object. Therefore, PhoneModel::adapt_setNode() must be called at least once for each model.
///
/// After adding the model statistics, iterate_setHelperMatrices_1() is called to iterate the 
/// statistics from the leaf-nodes into the helper matrices of the parent nodes.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Adapt_AM_TreeNode::setHelperMatrices_1(PhoneModel **models, int numberOfModels, bool doSil)
{
  for(int i=0;i<numberOfModels;i++)
  {
    if(models[i]->isSilModel() == doSil)
    {
      models[i]->adapt_setHelperMatrices();
    }
  }
  iterate_setHelperMatrices_1();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The first helper matrix contains the statistics of all acoustic models that are a member
/// of this tree node (cluster). This method calls the PhoneModel::adapt_setHelperMatrices() method
/// of all acoustic models. This method will call MixGaussian::MixGaussian() for each mixed gaussian
/// (triphone model) and this method will call Gaussian::adapt_setHelperMatrices() for each single
/// Gaussian in the model. Each Gaussian will then add its statistics to its own Adapt_AM_TreeNode
/// object. Therefore, PhoneModel::adapt_setNode() must be called at least once for each model.
///
/// After adding the model statistics, iterate_setHelperMatrices_1() is called to iterate the 
/// statistics from the leaf-nodes into the helper matrices of the parent nodes.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Adapt_AM_TreeNode::setHelperMatrices_1(PhoneModel **models, int numberOfModels)
{
  for(int i=0;i<numberOfModels;i++)
  {
    models[i]->adapt_setHelperMatrices();
  }
  iterate_setHelperMatrices_1();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will go iterativly into the child nodes and then add the statistics of these nodes
/// to the statistics of the parent node (left and right child statistics are added). The top node
/// contains the statistics from all gaussians of all models.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Adapt_AM_TreeNode::iterate_setHelperMatrices_1()
{
  if(childLeft != NULL && childRight != NULL)
  {
    PRINTF3("%p  %p\n",childLeft,childRight);
    childLeft->iterate_setHelperMatrices_1();
    childRight->iterate_setHelperMatrices_1();
  }
  else if(parent != NULL)
  {
    for(int i=0;i<AUGMENTED_MATRIX_SIZE;i++)
    {
      parent->AUGMENTED_ALL[i] += AUGMENTED_ALL[i];
    }    
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will train a Gaussian-model with the statistics of all acoustic models. These statistics
/// (for now the mean vector of each Gaussian) are obtained with setHelperMatrices_1().
/// The resulting Gaussian is split in two and its variance is multiplied with 1.2. The new created 
/// child-nodes will use these two Gaussians as meanGaussian. The samples are devided betweeen the two
/// child-nodes (the sample will be added to the node which mean vector is closest to the sample).
/// split() will return the amount of added child nodes. It will not split a node when the number
/// of samples is less than ADAPT_MIN_GAUSSIANS_PER_NODE.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int Adapt_AM_TreeNode::split(int depth)
{
  int res = 0;
  if(childLeft != NULL && childRight != NULL)
  {
    res += childLeft->split(depth+1);
    res += childRight->split(depth+1);
  }
  else 
  {
    meanGaussian->trainFinish();
    if(depth >= ADAPT_MAX_TREE_DEPTH)
    {
      C_factor = 1.0;
    }
    else
    {
      C_factor = (1.0 - 0.01) * (depth-1) / (ADAPT_MAX_TREE_DEPTH - 1) + 0.01;
    }
    PRINTF4("Should we split at depth %d: %g training samples, C_factor=%g\n",depth,denominator,C_factor);
    if(denominator > ADAPT_MIN_GAUSSIANS_PER_NODE && depth < 2.0*ADAPT_MAX_TREE_DEPTH)
    {
      PRINTF("  Yep...\n");
      childLeft                = new Adapt_AM_TreeNode(dimensions);
      childRight               = new Adapt_AM_TreeNode(dimensions);
      childLeft->meanGaussian  = new Gaussian(meanGaussian);
      meanGaussian->shiftMean();
      childRight->meanGaussian = meanGaussian;
      meanGaussian             = NULL;        
      childLeft->parent        = this;
      childRight->parent       = this;    
      childLeft->C_factor      = C_factor;  
      childRight->C_factor     = C_factor;        
      res = 2;
    }
  }
  return res;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// See the SMAPLR documentation for more information on the mathimatics used in this method.
/// This method calulates the adaptation matrix (W) and copies this matrix to the M matrix of 
/// its child-nodes. During the calculation of W, the helper matrices (AUGMENTED_ALL) need to be
/// filled by the setHelperMatrices_1() method.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Adapt_AM_TreeNode::setHelperMatrices_2()
{
  // Finish the Z matrix:
  for(int i=0;i<dimensions;i++)
  {
    AUGMENTED_ALL[i*AUG_X+i*AUG_Y+(dimensions+1)] += M[i][i]*C_factor;
  }
  // Finish the W_MULT matrix:
  for(int i=0;i<dimensions;i++)
  {
    for(int k=0;k<dimensions+1;k++)
    {
      AUGMENTED_ALL[i*AUG_X+k*AUG_Y+k] += C_factor;
    }
  }
  
  // Now iterate over i: calculate W row by row:
  for(int i=0;i<dimensions;i++)
  {
    // Use the Gauss-Jordan method to solve the augmented matrix AUGMENTED_ALL[i]...               
    // Now normalise all rows:
    for(int allRows=0;allRows<dimensions+1;allRows++) 
    {            
      if(AUGMENTED_ALL[i*AUG_X+allRows*AUG_Y] == 0.0)
      {
        // Move this row down as far as possible:
        int row = dimensions;
        while(row>allRows && AUGMENTED_ALL[i*AUG_X+row*AUG_Y] == 0.0)
        {
          row--;
        }
        // Now switch row and allRows:
        double switchRow;
        for(int switchR=0;switchR<dimensions+2;switchR++)
        {          
          switchRow = AUGMENTED_ALL[i*AUG_X+allRows*AUG_Y+switchR];
          AUGMENTED_ALL[i*AUG_X+allRows*AUG_Y+switchR] =
                             AUGMENTED_ALL[i*AUG_X+row*AUG_Y+switchR];
          AUGMENTED_ALL[i*AUG_X+row*AUG_Y+switchR] = switchRow;
        }
      }
      if(AUGMENTED_ALL[i*AUG_X+allRows*AUG_Y] != 0.0)
      {
        double factor = AUGMENTED_ALL[i*AUG_X+allRows*AUG_Y];
        for(int k=0;k<dimensions+2;k++)
        {
          AUGMENTED_ALL[i*AUG_X+allRows*AUG_Y+k] /= factor;
        }
      }
    }
    
    // Now iterate over the other rows.. (the last one is not needed, hence j<VECTORSIZE)
    for(int j=0;j<dimensions;j++) 
    {
      bool zeroFound = false;
      // Make all x-th columns zero
      for(int l=j+1;l<dimensions+1;l++)
      {
        if(AUGMENTED_ALL[i*AUG_X+l*AUG_Y+j] != 0.0)
        {
          assert(!zeroFound);
          for(int k=j;k<dimensions+2;k++)
          {
            AUGMENTED_ALL[i*AUG_X+l*AUG_Y+k] -=                           
                   AUGMENTED_ALL[i*AUG_X+j*AUG_Y+k];
          }        
        }
        else
        {
          zeroFound = true;
        }
      }
    
      // Now normalise all remaining rows:
      for(int allRows=j+1;allRows<dimensions+1;allRows++) 
      {
        if(AUGMENTED_ALL[i*AUG_X+allRows*AUG_Y+(j+1)] == 0.0)
        {
          // Move this row down as far as possible:
          int row = dimensions;
          while(row>allRows && AUGMENTED_ALL[i*AUG_X+row*AUG_Y+(j+1)] == 0.0)
          {
            row--;
          }
          // Now switch row and allRows:
          double switchRow;
          for(int switchR=j+1;switchR<dimensions+2;switchR++)
          {
            switchRow = AUGMENTED_ALL[i*AUG_X+allRows*AUG_Y+switchR];
            AUGMENTED_ALL[i*AUG_X+allRows*AUG_Y+switchR] = 
                   AUGMENTED_ALL[i*AUG_X+row*AUG_Y+switchR];
            AUGMENTED_ALL[i*AUG_X+row*AUG_Y+switchR] = switchRow;
          }
        }
        if(AUGMENTED_ALL[i*AUG_X+allRows*AUG_Y+(j+1)] != 0.0)
        {
          double factor = AUGMENTED_ALL[i*AUG_X+allRows*AUG_Y+(j+1)];
          for(int k=j+1;k<dimensions+2;k++)
          {        
            AUGMENTED_ALL[i*AUG_X+allRows*AUG_Y+k] /= factor;
          }
        }
      }                 
    }        
    // We have finished normalising matrix AUGMENTED_ALL[i]. Now calculate the values for the row w[i,-]:
    // Counting backwards:
    
    for(int j = dimensions; j >= 0; j--)
    {
      W[i][j] = AUGMENTED_ALL[i*AUG_X+j*AUG_Y+(dimensions+1)];
      for(int k = j-1; k >= 0; k--)
      {
        AUGMENTED_ALL[i*AUG_X+k*AUG_Y+(dimensions+1)] -=
           (W[i][j]*AUGMENTED_ALL[i*AUG_X+k*AUG_Y+j]);
      }      
    }
  }

  // Okay, Now M of our children is W of ourself. Let's copy the matrix and call our children...
  if(childLeft != NULL && childRight != NULL)
  {
    for(int i=0;i<dimensions;i++)
    {
      for(int j=0;j<dimensions+1;j++)
      {
        childLeft->M[i][j] = W[i][j];
        childRight->M[i][j] = W[i][j];
      }
    }  
    childLeft->setHelperMatrices_2();
    childRight->setHelperMatrices_2();    
  }
}

