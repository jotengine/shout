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
#ifndef ADAPT_AM_H
#define ADAPT_AM_H

#include "standard.h"
#include "trainphonemodel.h"

////////////////////////////////////////////////////////////////////////////
/// \brief The Adapt_AM class is the main class of the shout_adapt_am application. 
///
/// The constructor will be called by the main() function. All work will be
/// done in the constructor. The acoustic data will be loaded and passed to
/// the structured tree class Adapt_AM_TreeNode.
///
/// The general call order of the methods used is:
/// - TrainPhoneModel::adapt_setAcumulators()
/// - Create a new Adapt_AM_TreeNode object.
/// - Create the meanGaussian (Gaussian) object of the Adapt_AM_TreeNode object
/// - TrainPhoneModel::adapt_setInitialNode()
/// - Call Adapt_AM_TreeNode::split() and re-set the gaussians with TrainPhoneModel::adapt_setNode() (recursively)
/// - Call Adapt_AM_TreeNode::setHelperMatrices_1() and Adapt_AM_TreeNode::setHelperMatrices_2()
/// - The adaptation matrices are now calculated and can be used by Gaussian::adapt_adapt()
///   and (un-used) by Gaussian::adapt_unAdapt()
////////////////////////////////////////////////////////////////////////////

class Adapt_AM
{
  public:
    TrainPhoneModel **models;
  public:
    Adapt_AM(char *adapt_dir, char *amName, char *clusterAM, char *newAmName);
   ~Adapt_AM();

};

#endif
