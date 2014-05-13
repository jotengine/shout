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
#ifndef SHOUT_DCT2LEXTREE_H
#define SHOUT_DCT2LEXTREE_H

#include "vector.h"
#include "phonemodel.h"
#include "lexicaltree.h"
#include "languagemodel.h"



////////////////////////////////////////////////////////////////////////////
/// \brief This class will translate a plain text pronunciation dictionary into 
/// a PPT stored in LexicalTree format. This class is used by the application
/// shout_dct2lextree.
////////////////////////////////////////////////////////////////////////////

class Shout_dct2lextree : public LexicalTree
{
  public:
    Shout_dct2lextree(char *phoneListFile, char *dctFileName, char *lexFileName);
   ~Shout_dct2lextree();
   
  protected:
    void writeTree(LexicalNode *node, FILE *outFile);
    void splitList(char *str, char *str2);
};

#endif
