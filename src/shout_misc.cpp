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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "standard.h"
#include "shout-misc.h"
#include <inttypes.h>

#ifdef CELL_PROCESSOR
  #define little2host_short(x)     ((x>>8)&0xff)+((x << 8)&0xff00)
  #define little2host_long(x)      ((x & 0xff)<<24)+((x & 0xff00)<<8)+((x & 0xff0000)>>8)+((x >> 24)&0xff)
  #define little2host_longlong(x)  ((little2host_long(((x >> 32) & 0xffffffff)) | ((uint64_t)(little2host_long((x & 0xffffffff))) << 32)))
#else
  #define little2host_short(x)     (x)
  #define little2host_long(x)      (x)
  #define little2host_longlong(x)  (x)
#endif


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

int WriteFileLittleBigEndian::fwriteEndianSafe(void *var, int nrTimes, int nrBytes, FILE* file)
{
  int res = 0;
  while(nrTimes > 0)
  {
    switch(nrBytes)
    {
      case 1:
      {
        fwrite((char*)var,1,1,file);
      }
      break;
      case 2:
      {
        uint16_t wr;
        wr = little2host_short(*((uint16_t*) var)); 
        fwrite(&wr,1,2,file);
      }
      break;
      case 4:
      {
        uint32_t wr;
        wr = little2host_long(*((uint32_t*) var)); 
        fwrite(&wr,1,4,file);
      }
      break;
      case 8:
      {
        uint64_t wr;
        wr = little2host_longlong(*((uint64_t*) var)); 
        fwrite(&wr,1,8,file);
      }
      break;
      default:
      {
        assert(false);
      }
      break;
    }
    res += nrBytes;
    nrTimes--;
    var = (void*)&(((char *)var)[nrBytes]);
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

int WriteFileLittleBigEndian::freadEndianSafe(void *var, int nrTimes, int nrBytes, FILE* file)
{
  int res = 0;
  while(nrTimes > 0)
  {
    switch(nrBytes)
    {
      case 1:
      {
        if(fread((char*)var,1,1,file) != (unsigned int)nrBytes) return res;
      }
      break;
      case 2:
      {
        uint16_t wr;
        if(fread(&wr,1,2,file) != (unsigned int)nrBytes) return res;
        *((uint16_t*)var) = little2host_short(wr); 
      }
      break;
      case 4:
      {
        uint32_t wr;
        if(fread(&wr,1,4,file) != (unsigned int)nrBytes) return res;
        *((uint32_t*)var) = little2host_long(wr); 
      }
      break;
      case 8:
      {
        uint64_t wr;
        if(fread(&wr,1,8,file) != (unsigned int)nrBytes) return res;
        *((uint64_t*)var) = little2host_longlong(wr); 
      }
      break;
      default:
      {
        assert(false);
      }
      break;
    }
    res += nrBytes;
    nrTimes--;
    var = (void*)&(((char *)var)[nrBytes]);
  }
  return res;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This global function strips all spaces and tabs at the end of a string.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void StringFunctions::RightTrim(char *str)
{
  int k,n;
  n = k = strlen(str);
  if(k)
  {
    n--;
    while (n>=0 && (str[n]==' ' || str[n]=='\t')) n--;
    str[n+1]='\0';
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This function splits the first string in two parts at the place of the first space or tab 
/// character from the left. The second part is placed in str2, the first in str.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void StringFunctions::splitList(char *str, char *str2)
{
  strcpy(str2,"");
  int k;
  k = 0;
  int n = strlen(str);
  while (k<n && !(str[k]==' ' || str[k]=='\t'))
  {
    k++;
  }
  int kstop = k;
  while (k<n && (str[k]==' ' || str[k]=='\t'))
  {
    k++;
  }
  if(k < n)
  {
    strcpy(str2,&(str[k]));  
  }
  k = strlen(str2) - 1;  
  while(k >= 0)
  {
    if(str2[k] == '\n')
    {
      str2[k] = 0;
    }
    k--;
  } 
  k = strlen(str) - 1;  
  while(k >= 0)
  {
    if(str[k] == '\n')
    {
      str[k] = 0;
    }
    k--;
  } 
  str[kstop] = 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This function removes the label from the reference string. The label is stored in str2.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void StringFunctions::splitRef(char *str, char *str2)
{
  int k;
  k = 0;
  int n = strlen(str);
  while (k<n && str[k]!='(')
  {
    k++;
  }
  int kstop = k+1;
  while (k<n && str[k]!=')')
  {
    k++;
  }
  if(k==n)
  {
    printf("Error in REF file: %s\n",str);
    exit(1);
  }
  str[k] = 0;
  strcpy(str2,&(str[kstop]));  
  str[kstop-1] = 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Fast implementation of log by Laurent de Soras, adjusted to table-lookup by people at ICSI
/////////////////////////////////////////////////////////////////////////////////////////////////////

void FastMath::initLog()
{
  int i;
  float numlog;
  int *const exp_ptr = ((int*)&numlog);
  int x = *exp_ptr;
  x = 0x3F800000;
  *exp_ptr = x;
  int incr = 1 << (FASTLOG_23_MIN_N);
  for(i=0;i<FASTLOG_2_POW_20;i++)
  {
    logLookupTable[i] = log2(numlog);
    x+=incr;
    *exp_ptr = x;
  }
}

  
