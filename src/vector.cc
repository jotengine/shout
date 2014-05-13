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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "standard.h"
#include "vector.h"
#include "shout-misc.h"

using namespace WriteFileLittleBigEndian;


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor will initialize the Vector with value 0.0
/////////////////////////////////////////////////////////////////////////////////////////////////////

Vector::Vector(int length)
{
  vectorSize = length;
  if(length > 0)
  {
    theVector = new double[vectorSize];
    for(int i=0;i<vectorSize;i++)
    {    
      setValue(i, 0.0);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor will initialize the Vector with an alternative length and will copy the values
/// of th as values for the vector.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Vector::Vector(int length,double *th)
{
  vectorSize = length;
  theVector = new double[vectorSize];
  assert(th != NULL);
  memcpy(theVector,th,sizeof(double)*vectorSize);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor will set the vector size and will read the Vector from file.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Vector::Vector(FILE *inFile, int length, bool inputIsFloat)
{  
  vectorSize = length;
  theVector  = new double[vectorSize];
  if(inputIsFloat)
  {
    readFloatData(inFile);
  }
  else
  {
    readData(inFile);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor will concatinate the two vectors in a new vector.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Vector::Vector(Vector *vect1, Vector *vect2, Vector *vect3)
{
  int len1   = vect1->vectorSize;
  int len2   = vect2->vectorSize;
/*  int len3   = 0;
  if(vect3 != NULL)
  {
    len3 = vect3->vectorSize;
  }*/
  vectorSize = len1 + len2; // + len3;
  theVector  = new double[vectorSize];
  for(int i=0;i<len1;i++)
  {
    theVector[i] = vect1->theVector[i];
  }
  for(int i=0;i<len2;i++)
  {
    theVector[i+len1] = vect2->theVector[i];
  }
/*  for(int i=0;i<len3;i++)
  {
    theVector[i+len1+len2] = vect3->theVector[i];
  }*/
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor is empty. No vector data is stored dynamically and therefore no data needs to be
/// deleted.
/////////////////////////////////////////////////////////////////////////////////////////////////////

Vector::~Vector()
{
  if(theVector != NULL)
  {
    delete[] theVector;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// With StoreData(), the vector data is stored to an open file handle. This function is used to 
/// store Vector data embedded in the acoustic model (mean vector of Gaussian etc)
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Vector::storeData(FILE *outFile)
{
  fwrite(theVector,vectorSize,sizeof(double),outFile);  
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// With StoreFloatData(), the vector data is stored to an open file handle. All components are
/// written in 'float' format. This is usefull for storing feature vectors.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Vector::storeFloatData(FILE *outFile, int maxIndex)
{
  float t;
  if((maxIndex < 0) || (maxIndex > vectorSize))
  {
    maxIndex = vectorSize;
  }
  for(int i=0;i<maxIndex;i++)
  {
    t = ((float)theVector[i]);
    fwrite(&t,1,sizeof(float),outFile);  
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// ReadData() will read Vector data from an open file handler. This function is used to 
/// read Vector data embedded in the acoustic model (mean vector of Gaussian etc)
/// It returns the number of bytes actually read.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int Vector::readData(FILE *inFile)
{  
  return freadEndianSafe(theVector,vectorSize,sizeof(double),inFile);  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// ReadFloatData() is used to read Vector data from a feature file. In feature files, the vector 
/// data is stored as a floating number. Internaly, Vector uses double numbers.
/// It will return the actual number of bytes read.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int Vector::readFloatData(FILE *inFile)
{
  float t;
  int read = 0;
  for(int i=0;i<vectorSize;i++)
  {
    read += freadEndianSafe(&t,1,sizeof(float),inFile);
    theVector[i] = (double)(t);
  }
  return read;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

Vector *Vector::getMaxElement(bool toSelf)
{
  Vector *res = this;
  if(!toSelf)
  {
    res = new Vector(vectorSize);  
  }
  double max = theVector[0];
  int winner = 0;
  for(int i=1;i<vectorSize;i++)
  {
    if(theVector[i] > max)
    {
      max = theVector[i];
      winner = i;
    }
  }
  res->setAllValues(0.0);
  res->setValue(winner, max);
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Adapts a vector
/////////////////////////////////////////////////////////////////////////////////////////////////////

Vector *Vector::adapt()
{
  double tmp = 0.0;
  Vector *v = new Vector();
  for(int i=0;i<vectorSize;i++)
  {
    tmp = 0.0;
    for(int j=0;j<vectorSize+1;j++)
    {
      tmp += theVector[i*(vectorSize+1)+j];
    }
    v->theVector[i] = tmp;
  }
  return v;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Prints the Vector. Only used for debugging purposes.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Vector::printVector()
{
#ifdef DEBUG_ON
  int telMin = 0;
  if(false) //vectorSize == ASR_DEFAULT_VECTORSIZE)
  {
    int i = 0;
    printf("Feature: ");
    while(i < vectorSize/3)
    {
      printf("%10g ",getValue(i));
      i++;
    } 
    printf("\n");  
/*    printf("Delta  : ");
    while(i < vectorSize/3*2)
    {
      printf("%10g ",getValue(i));
      i++;
    } 
    printf("\n");  
    printf("Delta-D: ");
    while(i < vectorSize)
    {
      printf("%10g ",getValue(i));
      i++;
    } 
    printf("\n");  
*/
  
  }
  else
  {
    for(int i=0;i<vectorSize;i++)
    {    
      if(getValue(i) <= GAUSSIAN_VARIANCE_MINIMUM)
      {
        telMin++;
      }
      printf("%g ",getValue(i));
    } 
    printf("\n");
  //  printf("Number of minimum variance values in the variance vector: %d\n",telMin);
  }
#endif
}

