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
#include "phonefilereader.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

PhoneFileReader::PhoneFileReader(char *phoneListFile, DataStats **data, int *numberOfPhones, int *numberOfSil)
{
  int  i = 0;
  int  state = 0;
  char str[MAXSTR];
  int  phoneNum = 0;
  
  printf("Reading phones...\n");
  FILE *inFile = fopen(phoneListFile,"rt");
  while(!feof(inFile))
  {          
          char *resChar = fgets(str,MAXSTR,inFile);
          assert(resChar != NULL);
          switch(state)
          {
            case 0:        
              if(strncmp(str,"Number of phones:",17) == 0)
              {
                *numberOfPhones = atoi(&(str[18]));
                if(*numberOfPhones >= 1)
                {                
                  (*data) = new DataStats[*numberOfPhones];
                } 
                else
                {
                  USER_ERROR("The phonelist file is corrupt. Number of phones is smaller than one!");
                }            
                state++;
              }
              else
              {
                USER_ERROR("The phonelist file is corrupt! First line shout be: 'Number of phones: [number]'");
              }
              break;
            case 1:
              if(strncmp(str,"Number of SIL's: ",17) == 0)
              {
                *numberOfSil = atoi(&(str[18]));
                if(*numberOfSil < 1)
                {                
                  USER_ERROR("The phonelist file is corrupt. Number of SIL phones is smaller than one!");
                }            
                state++;
              }
              else
              {
                USER_ERROR("The phonelist file is corrupt! Second line shout be: 'Number of SIL's: [number]'");
              }
              break;
            case 2:
              if(phoneNum < *numberOfPhones)
              {            
                if(strlen(str) > 9)
                {
                  USER_ERROR("Sorry, the length of a phone name may not exceed 9.");
                }
                i=0;
                while(i<10)
                {
                  if(str[i] == '\r' || str[i] == ' ' || str[i] == '\t' || str[i] == '\n' || str[i] == 0)
                  {
                    while(i<10)
                    {
                      str[i] = ' '; 
                      i++;                     
                    }
                  } 
                  i++;
                }       
                str[10] = 0;
                strncpy((*data)[phoneNum].name,str,11);
                phoneNum++;
              }
              break;
            default:
              USER_ERROR("The phonelist file is corrupt!");
              break;
          }
  }
  fclose(inFile);    
}


PhoneFileReader::~PhoneFileReader()
{
}


