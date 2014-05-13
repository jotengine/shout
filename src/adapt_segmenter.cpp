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
#include "adapt_segmenter.h"

#include <math.h>

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Initialize...
/////////////////////////////////////////////////////////////////////////////////////////////////////


Adapt_Segmenter::Adapt_Segmenter(FeaturePool *fp, int sID_train, int sID_decode, char *name, int sysID)
                :Train_Segmenter(fp, sID_train, sID_decode, name,1,8,50,false, NULL)
{ 
  compareModels   = NULL;
  for(int i=1;i<numberOfPhones;i++)
  {
    compareClusters[i] = 0;
  }  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Clean up..
/////////////////////////////////////////////////////////////////////////////////////////////////////

Adapt_Segmenter::~Adapt_Segmenter()
{
/*
  if(compareModels != NULL)
  {
    for(int a=0;a<nrCompareModels;a++)
    {
      if(compareModels[a] != NULL)
      {
        delete compareModels[a];
      }
    }
    delete[] compareModels;
  }


  if(mergeTestMod != NULL)
  {
    for(int i=0;i<numberOfWords;i++)
    {
      if(mergeTestMod[i] != NULL)
      {
        delete mergeTestMod[i];
      }
    }
    delete[] mergeTestMod;
  }
  if(adaptedModel != NULL)
  {
    for(int i=0;i<numberOfWords;i++)
    {
      for(int i2=0;i2<ubm_word_amount;i2++)
      {
        if(adaptedModel[i+i2*numberOfWords] != NULL)
        {
          if(trainCluster[i] == adaptedModel[i+i2*numberOfWords])
          {
            trainCluster[i] = NULL;
          }
          delete adaptedModel[i+i2*numberOfWords];
        }
      }
    }
    delete[] adaptedModel;
  }
  if(masterAdapt != NULL)
  {
    for(int i=0;i<UBM_AMOUNT;i++)
    {
      if(masterAdapt[i] != NULL)
      {
        delete masterAdapt[i];
      }
    }
    delete[] masterAdapt;
  }
  if(helpID != NULL)
  {
    delete[] helpID;
  }
  if(pTable1 != NULL)
  {
    delete[] pTable1;
  }
  if(pTable2 != NULL)
  {
    delete[] pTable2;
  }
  */
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo add documentation
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Adapt_Segmenter::adaptModel(TrainPhoneModel *model)
{  
  int dim = model->dim();
  Adapt_AM_TreeNode *adaptationTree = new Adapt_AM_TreeNode(dim);
  adaptationTree->meanGaussian = new Gaussian(dim);  
  model->adapt_setInitialNode(adaptationTree); 
  adaptationTree->setHelperMatrices_1(((PhoneModel**)&model), 1,model->isSilModel()); 
  adaptationTree->setHelperMatrices_2();
  model->adapt_adapt();
  model->adapt_clear();
  delete adaptationTree;
}


  


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Adapt_Segmenter::startMergeIteration()
{

}



/////////////////////////////////////////////////////////////////////////////////////////////////////
/// We have decided on two merge candidates. Now determine if we should actually do it, or just stop.
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool Adapt_Segmenter::proceedMerge(int model1, int model2, int method)
{
  PRINTF("Deleting old backup-models..\n");
  if(compareModels != NULL)
  {
    for(int a=1;a<numberOfPhones;a++)
    {
      if(compareModels[a] != NULL)
      {
        delete compareModels[a];
      }
    }
    delete[] compareModels;
  }
  compareModels = new TrainPhoneModel*[numberOfPhones];

  compareModels[0] = NULL;
  
  int segTryID = featurePool->claimSegmentationID();
  
  featurePool->segmentationCopy(segID,segTryID);
  
  int    winner      = 0;
  double score[numberOfPhones];
  int    data1[numberOfPhones];
  double totScore    = 0;
  PRINTF("Start the party..\n");
  for(int i=1;i<numberOfPhones;i++)
  {
    if(trainCluster[i] != NULL)
    {
      trainCluster[i]->setTrainingData(featurePool,segTryID,i);
      data1[i] = featurePool->getClusterLength(segTryID,i);        
      double s = trainCluster[i]->getTrainSilP();
      score[i] = s / ((double)data1[i]);
      totScore += s;
    }
  }
  
  bool doAgain = true;
  while(doAgain)
  {
    featurePool->segmentationCopy(segID,segTryID);
    
    // Picking candidate:
    double lowestScore = 9.0e300;
    winner = 0;
    for(int i=1;i<numberOfPhones;i++)
    {
      if(trainCluster[i] == NULL)
      {
        compareModels[i] = NULL;
      }
      else
      {
        compareModels[i] = new TrainPhoneModel(trainCluster[i]);
      }
    }  
    if(winner == 0)
    {
      for(int i=1;i<numberOfPhones;i++)
      {
        if(trainCluster[i] != NULL)
        {
          if((score[i] < lowestScore) && (compareClusters[i] == 0))
          {
            lowestScore = score[i];
            winner = i;
          }
        }
      }
    }
    if(winner > 0)
    {
      winnerM1 =  1;
      winnerM2 =  1;
  
      int nrGaussians = trainCluster[winner]->getNumberOfGaussians();
      // Now actually merge...
      delete trainCluster[winner];
      trainCluster[winner] = NULL;
assert(false);
//      createLexicalTree(CLUSTER_MIN_NUMBER_FRAMES);
      segmentFeaturePool(sadID_train,0,segTryID);
      
      int data2[numberOfPhones];
      int addG[numberOfPhones];
      int remainG[numberOfPhones];
      int difData = 0;
      for(int i=1;i<numberOfPhones;i++)
      {
        if(trainCluster[i] != NULL)
        {
          data2[i] = featurePool->getClusterLength(segTryID,i);
          if(data2[i]>data1[i])
          {
            difData += (data2[i]-data1[i]);
          }
        }
      }
      PRINTF3("Dividing %d gaussians over %d data samples..\n",nrGaussians,difData);
      double gaussianPerData = ((double)nrGaussians) / ((double)difData);
      int added = 0;
      for(int i=1;i<numberOfPhones;i++)
      {
        if(trainCluster[i] != NULL)
        {
          if(data2[i]>data1[i])
          {
            double a = ((double)(data2[i]-data1[i]))*gaussianPerData;
            addG[i]    = ((int)a);
            remainG[i] = ((int)((a-((double)addG[i]))*100.0));
            added += addG[i];
          }
          else
          {
            addG[i] = 0;
          }
        }
      }
      for(int i=1;i<numberOfPhones;i++)
      {
        if(trainCluster[i] != NULL)
        {
          PRINTF4("SPK%02d: %2d %4d\n",i,addG[i],remainG[i]);
        }
      }
            
      
      PRINTF3("Added, still have to add: %d - %d\n",added,nrGaussians-added);
      
      int i=1;
      while((added < nrGaussians) && (i<numberOfPhones))
      {
        int highest = 0;
        i = numberOfPhones;
        for(int i2=1;i2<numberOfPhones;i2++)
        {
          if((trainCluster[i2] != NULL) && (data2[i2]>data1[i2]) && (remainG[i2] > highest))
          {
            highest = remainG[i2];
            i = i2;
          }
        }
        if(i<numberOfPhones)
        {
          remainG[i] = 0;
          addG[i]++;
          added++;
          PRINTF2("One gaussian added to SPK%02d\n",i);
        }
      }
      
      
      
      for(int i=1;i<numberOfPhones;i++)
      {
        if(trainCluster[i] != NULL)
        {
          trainCluster[i]->setTrainingData(featurePool,segTryID,i);
          trainCluster[i]->train(trainCluster[i]->getNumberOfGaussians()+addG[i],true,true);
        }
      }
      if(segmentFeaturePool(sadID_train,0,segTryID) < totScore)
      {
        PRINTF2("Sorry, no luck with this one (SPK%02d).. Let's start over...\n",winner);
        compareClusters[winner] = 1;
        
        for(int i=1;i<numberOfPhones;i++)
        {
          if(compareModels[i] != NULL)
          {
            if(trainCluster[i] != NULL)
            {
              delete trainCluster[i];
            }
            trainCluster[i]  = compareModels[i];
            compareModels[i] = NULL;
          }
        }
      }
      else
      {
        PRINTF2("YEAH, WE GOT RID OF SPK%02d!!!!\n",winner);
        doAgain = false;
      }
    }
    else
    {
      winnerM1 = -1;
      winnerM2 = -1;
      doAgain  = false;
    }
  }
  
  if(compareModels != NULL)
  {
    for(int i=1;i<numberOfPhones;i++)
    {
      compareModels[i] = NULL;
    }
  }
  
  for(int i=1;i<numberOfPhones;i++)
  {
    if(compareClusters[i] > 0)
    {
      compareClusters[i]--;
    }
  }
  PRINTF("Done merging...\n");
  return (winnerM1>0);
}

 
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Merges two single models (over loaded by Adapt_Segmenter)
/////////////////////////////////////////////////////////////////////////////////////////////////////

void Adapt_Segmenter::mergeModels(int model1, int model2)
{
/*
  // Let's cut back on gaussians:
  for(int i=1;i<numberOfPhones;i++)
  {
    if((trainCluster[i] != NULL) && (trainCluster[i]->getNumberOfGaussians() > 50))
    {
      trainCluster[i]->setMaxGaussians(25);
      trainCluster[i]->train(-1,true,true);
    }
  }
  */
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Merges two single models and determine the score
/////////////////////////////////////////////////////////////////////////////////////////////////////


double Adapt_Segmenter::getMergeModelScore(int model1, int model2, int useMethod)
{
  return 100.0;
}



