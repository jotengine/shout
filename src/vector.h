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
 
#ifndef VECTOR_H
#define VECTOR_H

#include "standard.h"
#include <math.h>

/*
extern "C"
{
#include </n/shokuji/da/marijn/ATLAS/include/cblas.h>
}
*/

////////////////////////////////////////////////////////////////////////////
/// \brief This class handles basic vector calculations that are needed for
/// the system. 
///
/// Audio feature data is stored as a Vector, but also the Gaussian (and MixGaussian)
/// parameters mean and variance are stored in Vector format.
/// All methods are optimized for use in Shout. It is possible to store and
/// load data from disc or send and receive them over a socket connection.
////////////////////////////////////////////////////////////////////////////

class Vector
{
  protected:
    // Attributes:
    double     *theVector;
    int         vectorSize;
    // Methods:   
    
  public:  
    // Constructor and destructor:
             Vector                     (Vector *vect1, Vector *vect2, Vector *vect3);
             Vector                     (int length = ASR_DEFAULT_VECTORSIZE);
             Vector                     (int length, double *th);
             Vector                     (FILE *inFile, int length = ASR_DEFAULT_VECTORSIZE, bool inputIsFloat = false);
            ~Vector                     ();
    
    // INLINE FUNCTIONS FOR CODE OPTIMALIZATION:
    inline double getValue(int element) const 
    {
      assert(element >=0 && element < vectorSize); 
      return theVector[element];
    }

    inline double* getVectorArray() const
    {
      return theVector;
    }

    inline int len() const 
    {
      return vectorSize;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Part of the used formula is: sum of all elements(e) { epow((vector[e]-mean[e])^2 / variance[e] ) }
    /// This method calculates this part of the formula for the Gaussian object.
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline double calculateEPart(Vector *mean, Vector *variance)
    {
      double res = 0.0;
      double tmp = 0.0;
      double *mean_theVector     = mean->theVector;      // code optimalisation
      double *variance_theVector = variance->theVector;  // code optimalisation
      for(int i=0;i<vectorSize;i++)
      {
        tmp = (theVector[i] - mean_theVector[i]);
        res += (tmp/variance_theVector[i]*tmp);
      }    
      return (-0.5*res);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \todo docs
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    inline void *getKey()
    {
      return ((void*)theVector);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Adds 'val' to vector element 'element'. The elements from the Vector are in the range of
    /// zero to vectorSize-1.
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline void addValue(int element, double val)
    {
      assert(element >= 0 && element < vectorSize);
      theVector[element] += val;  
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Sets vector element 'element' with value 'val'. The elements from the Vector are in the range of
    /// zero to vectorSize-1.
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline void setValue(int element, double val)
    {
      assert(element >= 0 && element < vectorSize);
      theVector[element] = val;  
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Sets all vector elements with value 'val'. The elements from the Vector are in the range of
    /// zero to vectorSize-1.
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline void setAllValues(double val)
    {
      for(int i=0;i<vectorSize;i++)
      {    
        theVector[i] = val;  
      }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \todo Docs
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline void isSmaller(Vector *compare, Vector *mask)
    {
      double *compareVector = compare->theVector;
      double *maskVector    = mask->theVector;
      for(int i=0;i<vectorSize;i++)
      {
        if(theVector[i] < compareVector[i])
        {
          maskVector[i]++;
        }
      }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Borrows the vector of another object..
    /// Be sure to set it to NULL before you kill this object!
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline void borrowVector(Vector *totVect, int first, int last)
    {
      vectorSize = last - first + 1;
      if(totVect != NULL)
      {
        theVector = &(totVect->theVector[first]);
      }
      else
      {
        theVector = NULL;
      }
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Sets all vector elements with the same value as in the input vector.
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline void setAllValues(const Vector *inVector)
    {
      double *resVector = inVector->theVector;
      assert(inVector->len() == vectorSize);
      for(int i=0;i<vectorSize;i++)
      {    
        theVector[i] = resVector[i];  
      }
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Calculates the distance between the input vector and this vector, using: Sum( (v1[i]-v2[i])^2 )
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline double distance(Vector *inVector)
    {
      double res = 0.0;
      double tmp = 0.0;
      assert(vectorSize == inVector->len());
      for(int i=0;i<vectorSize;i++)
      {
        tmp = (theVector[i] - inVector->theVector[i]);
        res += tmp*tmp;
      }
      return res;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// This method returns the multiplication of the Vector object and the input parameter Vector m. 
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline double multiplyVector(Vector *m)
    {
      double res = 0.0;
      assert(m->len() == vectorSize);
//      cblas_dgemv (CblasRowMajor, CblasNoTrans, 3, 3, 1.0, theVector, 3, m->theVector, 1, 0.0, m->theVector, 1);      
      for(int i=0;i<vectorSize;i++)
      {    
        res += theVector[i] * m->theVector[i];
      }
      return res;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// All elements of input Vector m are multiplied with the values of this Vector object. The results
    /// are either stored in a new Vector object, of in the current object, depending on the value of the
    /// boolean input variable 'toSelf'.
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline Vector* multiplyElements(bool toSelf, Vector *m)
    {
      if(toSelf)
      {
        for(int i=0;i<vectorSize;i++)
        {    
          theVector[i] *= m->theVector[i];
        }
        return this;
      }
      else
      {
        assert(m->len() == vectorSize);
        Vector *res = new Vector(vectorSize);  
        for(int i=0;i<vectorSize;i++)
        {    
          res->theVector[i] = theVector[i] * m->theVector[i];
        }
        return res;    
      }
      return NULL;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// All elements of this Vector object are divided by the elements of input Vector m. The results
    /// are either stored in a new Vector object, of in the current object, depending on the value of the
    /// boolean input variable 'toSelf'.
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline Vector* divideElements(bool toSelf, Vector *m)
    {
      Vector *res = this;
      if(!toSelf)
      {
        assert(m->len() == vectorSize);
        res = new Vector(vectorSize);  
      }
      for(int i=0;i<vectorSize;i++)
      {    
        res->theVector[i] = theVector[i] / m->theVector[i];
      }
      return res;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// All elements of input Vector m are added to the values of this Vector object. The results
    /// are either stored in a new Vector object, of in the current object, depending on the value of the
    /// boolean input variable 'toSelf'.
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline Vector* addElements(bool toSelf, Vector *m)
    {
      Vector *res = this;
      if(!toSelf)
      {
        assert(m->len() == vectorSize);
        res = new Vector(vectorSize);  
      }
      for(int i=0;i<vectorSize;i++)
      {    
        res->theVector[i] = theVector[i] + m->theVector[i];
      }
      return res;
    }
        
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// All elements of of this Vector object are multiplied with the input variable m (double). The results
    /// are either stored in a new Vector object, of in the current object, depending on the value of the
    /// boolean input variable 'toSelf'.
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline Vector* multiplyElements(bool toSelf, double m)
    {
      if(toSelf)
      {
        for(int i=0;i<vectorSize;i++)
        {    
          theVector[i] *= m;
        }
        return this;
      }
      else
      {
        Vector *res = new Vector(vectorSize);  
        for(int i=0;i<vectorSize;i++)
        {    
          res->theVector[i] = theVector[i] * m;
        }
        return res;    
      }
      return NULL;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// All elements of input Vector m are substracted from the values of this Vector object. The results
    /// are either stored in a new Vector object, of in the current object, depending on the value of the
    /// boolean input variable 'toSelf'.
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline Vector* substractElements(bool toSelf, Vector *m)
    {
      Vector *res = this;
      if(!toSelf)
      {
        assert(m->len() == vectorSize);
        res = new Vector(vectorSize);  
      }
      for(int i=0;i<vectorSize;i++)
      {    
        res->theVector[i] = theVector[i] - m->theVector[i];
      }
      return res;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// The input parameter m (double) is added to all of the values of this Vector object. The results
    /// are either stored in a new Vector object, of in the current object, depending on the value of the
    /// boolean input variable 'toSelf'.
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline Vector* addElements(bool toSelf, double m)
    {
      
      Vector *res = this;
      if(!toSelf)
      {
        res = new Vector(vectorSize);  
      }
      for(int i=0;i<vectorSize;i++)
      {    
        res->theVector[i] = theVector[i] + m;
      }
      return res;
    }

    inline Vector* sqrtElements(bool toSelf)
    {
      
      Vector *res = this;
      if(!toSelf)
      {
        res = new Vector(vectorSize);  
      }
      for(int i=0;i<vectorSize;i++)
      {    
        res->theVector[i] = sqrt(theVector[i]);
      }
      return res;
    }

    inline Vector* pow2Elements(bool toSelf)
    {
      
      Vector *res = this;
      if(!toSelf)
      {
        res = new Vector(vectorSize);  
      }
      for(int i=0;i<vectorSize;i++)
      {    
        res->theVector[i] = theVector[i]*theVector[i];
      }
      return res;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Makes a copy of the current Vector object and returns it.
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline Vector* copyVector(void)
    {
      Vector *cpy = new Vector(vectorSize,theVector);
      return cpy;
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Returns the determinant of the current Vector object.
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline double getDeterminant(void)
    {
      double res = 1.0; 
      for(int i=0;i<vectorSize;i++)
      {    
        res *= theVector[i];
      }
      return res; 
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Checks if one of the elements of this Vector is negative. If so, the method returns true, otherwise
    /// it returns false.
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline bool checkNegativeElements(void)
    {
      for(int i=0;i<vectorSize;i++)
      {
        if(theVector[i]<0)
        {
          return true; 
        }
      }
      return false;  
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Used in order to make sure that all elements of the Vector are bigger than 'min'. If one of the 
    /// elements is too small, calculating the determinant will return zero. And dividing by the determinant
    /// would mean crashing.. Therefore, elements from the variance Vector may not be too small.
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    inline bool minimizeElements(double min)
    {
      bool tooSmall = false;
      double res;
      for(int i=0;i<vectorSize;i++)
      {
        res = theVector[i];
        if((res < min) || (!(res <= 0) && !(res >= 0)))
        {
          tooSmall     = true;
          theVector[i] = min;
        }
      }
      return tooSmall;
    }
    

    
        
    // Methods:
/*    
    double   distance                   (Vector *inVector);
    void     setValue                   (int element, double value);
    void     addValue                   (int element, double value);
    void     setAllValues               (double value);
    void     setAllValues               (const Vector *inVector);
    double   multiplyVector             (Vector *m);
    Vector*  multiplyElements           (bool toSelf, Vector *m);
    Vector*  multiplyElements           (bool toSelf, double m);    
    Vector*  divideElements             (bool toSelf, Vector *m);
    Vector*  substractElements          (bool toSelf, Vector *m);
    Vector*  addElements                (bool toSelf, Vector *m);
    Vector*  addElements                (bool toSelf, double m);
    void     minimizeElements           (double min);
    double   getDeterminant             (void);      
    Vector*  copyVector                 (void);
    bool     checkNegativeElements      (void);    
*/    
    
    Vector  *adapt                      ();
    void     storeData                  (FILE *outFile);
    void     storeFloatData             (FILE *outFile, int maxIndex = -1);
    int      readData                   (FILE *inFile);
    int      readFloatData              (FILE *inFile);
    Vector  *getMaxElement              (bool toSelf);

    // Helper Methods:
    void     printVector                ();
};

#endif // VECTOR_H
