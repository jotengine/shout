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
#include <time.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/resource.h>

#include "standard.h"
#include "memmappedfile.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

MemMappedFile::MemMappedFile(int objSize, int nrObj, int memSize, const char *fileN, bool isWr)
{
  isWritable   = isWr;
  pageSize     = sysconf(_SC_PAGE_SIZE);
  objectSize   = objSize;
  totNrObjects = nrObj;
  firstObject  = 0;
  memBlockSize = (memSize/pageSize + 1) * pageSize;
  nrObjects    = memBlockSize / objectSize - 1;

  assert(nrObjects > 1);

  FILE *memFile = fopen64(fileN,  "wb");
  char tmp      = 1;
  int written   = 0;
  int writeSize = pageSize*100;
  while(written < objectSize * totNrObjects + memBlockSize)
  {
    fseek(memFile, writeSize-1, SEEK_CUR);
    fwrite(&tmp,1,1,memFile);
    written += writeSize;
  }
  fclose(memFile);

  fileDescriptor = open64(fileN, O_RDWR);
  assert(fileDescriptor >= 0);
  memoryPool = (char*) mmap(0, memBlockSize, PROT_READ | PROT_WRITE, MAP_SHARED, fileDescriptor, 0);
  assert(memoryPool != (caddr_t) - 1);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

MemMappedFile::~MemMappedFile()
{
  if(isWritable)
  {
    msync(memoryPool,memBlockSize,MS_SYNC);
  }
  munmap(memoryPool,memBlockSize);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

char *MemMappedFile::map(int i)
{
  int a = i - firstObject;
  if(a >= 0 && a < nrObjects)
  {
    return &memoryPool[a*objectSize];
  }
  else
  {
    if(isWritable)
    {
      msync(memoryPool,memBlockSize,MS_SYNC);
    }
    munmap(memoryPool,memBlockSize);
    int page = objectSize*i / pageSize;
    firstObject = page*pageSize / objectSize;

    if(firstObject*objectSize < page*pageSize)
    {
      firstObject++;
    }
    int byteDiff = firstObject*objectSize - page*pageSize;

    memoryPool = (char*) mmap(0, memBlockSize, PROT_READ | PROT_WRITE, MAP_SHARED, fileDescriptor, page*pageSize);
    assert(memoryPool != (caddr_t) - 1);
    a = i - firstObject;
    return &memoryPool[a*objectSize];
  }
  return NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void MemMappedFile::setWritePermission(bool isWr)
{
  isWritable = isWr;
}
