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
#include "lexicaltree.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "standard.h"
#include "shout-misc.h"

using namespace WriteFileLittleBigEndian;

#define START_OF_SENTENCE       "<s>"
#define END_OF_SENTENCE        "</s>"

#define MIN_LIKELIHOOD         (-9e300)
#define absolute(x)            (((x)<(0.0))?(-1.0*x):(x))


#ifdef SEARCH_STATISTICS_ON
  extern int nrOfStateHistPruning;
#endif    

#ifdef DEBUG_MEMORY
  int gaussianCount = 0;
  int lmlaCount     = 0;
  int pdfCount      = 0;
#endif

#define MAX_ANALYSIS_NBEST_PATHS  (50000)
#define MAX_NDEPTH                (50000)
extern bool doPhoneAlignment;

TokenType *bestTokenInSystem;


#ifdef MAX_NR_THREADS
#include <pthread.h>

#ifdef CELL_PROCESSOR
  #include <libspe2.h>
  extern spe_program_handle_t dma_getfastp_spu;
  extern spe_program_handle_t dma_createlmla_spu;
#endif

struct Thread_LMLACalculation_Data
{
  LMLAGlobalListType *lmlaGlobal;
  float *allLMP;
  int numberOfCompressedNodes;
  FastCompressedTree *fastCompressedTree;
};

bool lmla_threadMayStart      = false;
bool lmla_threadFinished      = true;
pthread_t                   lmla_thread;
Thread_LMLACalculation_Data lmla_dataThread;
pthread_mutex_t             lmla_condition_mutexStart;
pthread_cond_t              lmla_condition_threadStart;
pthread_mutex_t             lmla_condition_mutexDone;
pthread_cond_t              lmla_condition_threadDone;


void *thread_lmlaCalculation( void *ptr )
{
  FastCompressedTree *fastCompressedTree = lmla_dataThread.fastCompressedTree;
  int numberOfCompressedNodes            = lmla_dataThread.numberOfCompressedNodes;
// Maybe... maybe we can use the SPEs! Yeeaj!
#ifdef CELL_PROCESSOR
    float                      *result_mem;
    control_block_createlmla_t *control_block;

    int rc1 = posix_memalign((void**)(&result_mem),      128, numberOfCompressedNodes * sizeof(float));
    int rc2 = posix_memalign((void**)(&control_block),   128, sizeof(control_block_createlmla_t));
    if (rc1 != 0 || rc2 != 0)
    {
      USER_ERROR("PS3: Failed to allocate memory for SPE communication");
    }


    union 
    {
      unsigned long long ull;
      unsigned int ui[2];
    } long_addr;
    /* Create SPE threads to execute 'dma_example_spu'. */
    spe_context_ptr_t spe_context = spe_context_create(0, NULL);
    if (spe_context == NULL)
    {
      USER_ERROR("PS3: Failed creating SPE context");
    }
    /* Load SPE program into context */
    if (spe_program_load(spe_context, &dma_createlmla_spu))
    {
      USER_ERROR("PS3: Failed loading SPE program");
    }
    long_addr.ui[0] = 0;
    long_addr.ui[1] = (unsigned int)result_mem;
    control_block->out_addr = long_addr.ull;
#endif


  while(true)
  {
    pthread_mutex_lock( &lmla_condition_mutexStart );
    while(!lmla_threadMayStart)
    {
      pthread_cond_wait( &lmla_condition_threadStart, &lmla_condition_mutexStart );
    }
    lmla_threadMayStart = false;
    pthread_mutex_unlock( &lmla_condition_mutexStart );


    LMLAGlobalListType *lmlaGlobal = lmla_dataThread.lmlaGlobal;
    float *allLMP                  = lmla_dataThread.allLMP;
    float fastCompressedTreeBackValue[numberOfCompressedNodes];
    float *bv              = fastCompressedTreeBackValue;
    float prefRes = 0.0;
    for(int i=0;i<numberOfCompressedNodes;i++)
    {
      float res = SMALLEST_NEGATIVE;
      if(fastCompressedTree[i].wordID == -10)
      {
        res = 0.0;
      }
      else if(fastCompressedTree[i].wordID >= 0)
      {
        res = allLMP[fastCompressedTree[i].wordID];
      }
      else if(fastCompressedTree[i].nextIndex >= 0)
      {
        res = fastCompressedTreeBackValue[fastCompressedTree[i].nextIndex];
      }
      lmlaGlobal->lookAhead[i] = res;
  
      if(fastCompressedTree[i].hasPar && prefRes > res)
      {
        res = prefRes;
      }
      prefRes    = res;
      *bv        = res;
      bv++;
    }
    delete[] lmla_dataThread.allLMP;
    lmla_dataThread.allLMP = NULL;
    COUNTER_MINUS(lmlaCount);

    pthread_mutex_lock( &lmla_condition_mutexDone);
    lmla_threadFinished = true;
    pthread_cond_signal( &lmla_condition_threadDone);
    pthread_mutex_unlock( &lmla_condition_mutexDone);
  }
}






struct Thread_pdfCalculation_Data
{
  int threadID;
  int length;
  MixGaussian **mixGaussian;
  double      **result;
  double       *vList;
  bool          doSecond;
  int           time;
};

bool threadMayStart[MAX_NR_THREADS];
int  nrDataThreadsFinished = 0;
pthread_t                   thread[MAX_NR_THREADS];
Thread_pdfCalculation_Data  dataThread[MAX_NR_THREADS];
pthread_mutex_t             condition_mutexStart[MAX_NR_THREADS];
pthread_cond_t              condition_threadStart[MAX_NR_THREADS];
pthread_mutex_t             condition_mutexDone;
pthread_cond_t              condition_threadDone;


void *thread_pdfCalculation( void *ptr )
{
  int threadID = ((Thread_pdfCalculation_Data*)ptr)->threadID;


// Maybe... maybe we can use the SPEs! Yeeaj!
#ifdef CELL_PROCESSOR
    float                    *result_mem;
    control_block_getfastp_t *control_block;

    int rc1 = posix_memalign((void**)(&result_mem),      128, MAX_PDF_CALC_SIZE * sizeof(float));
    int rc2 = posix_memalign((void**)(&control_block),   128, sizeof(control_block_getfastp_t));
    if (rc1 != 0 || rc2 != 0)
    {
      USER_ERROR("PS3: Failed to allocate memory for SPE communication");
    }


    union 
    {
      unsigned long long ull;
      unsigned int ui[2];
    } long_addr;
    /* Create SPE threads to execute 'dma_example_spu'. */
    spe_context_ptr_t spe_context = spe_context_create(0, NULL);
    if (spe_context == NULL)
    {
      USER_ERROR("PS3: Failed creating SPE context");
    }
    /* Load SPE program into context */
    if (spe_program_load(spe_context, &dma_getfastp_spu))
    {
      USER_ERROR("PS3: Failed loading SPE program");
    }
    long_addr.ui[0] = 0;
    long_addr.ui[1] = (unsigned int)result_mem;
    control_block->out_addr = long_addr.ull;
#endif


  while(true)
  {
    pthread_mutex_lock( &condition_mutexStart[threadID] );
    while(!threadMayStart[threadID])
    {
      pthread_cond_wait( &condition_threadStart[threadID], &condition_mutexStart[threadID] );
    }
    threadMayStart[threadID] = false;
    pthread_mutex_unlock( &condition_mutexStart[threadID] );


    Thread_pdfCalculation_Data *data = (Thread_pdfCalculation_Data*)ptr;
    if(data->length <= 0)
    {
#ifdef CELL_PROCESSOR
      free(result_mem);
      free(control_block);
#endif
      return NULL;
    }


#ifdef CELL_PROCESSOR
  unsigned int entry = SPE_DEFAULT_ENTRY;
  double resultSPE[MAX_PDF_CALC_SIZE];

  control_block->nrPDFs   = data->length;
  for(int i=0;i<39;i++)
  {
    control_block->vector[i] = (float)data->vector->getValue(i);
  }
  for(int i=0;i<data->length;i++)
  {
    long_addr.ui[0] = 0;
    long_addr.ui[1] = (unsigned int)data->mixGaussian[i]->fastP_CPandWeight_meanAndVar;
    control_block->in_addr[i]   = long_addr.ull;
    control_block->in_length[i] = data->mixGaussian[i]->getNumberOfGaussians();
  }
  /* start the SPE thread. We want to pass the pointer to the control block to the SPE */
  if (spe_context_run(spe_context, &entry, 0, control_block, NULL, NULL) < 0) 
  {
    USER_ERROR("PS3: Failed to run application on SPE");
  }
  MixGaussian **pdfUpdateList = data->mixGaussian;
  double      **resultList    = data->result;
  for(int i=0;i<data->length;i++)
  {
    *resultList[i] = (double)(result_mem[i]);
  }
#else
    MixGaussian **pdfUpdateList = data->mixGaussian;
    double      **resultList    = data->result;
    bool          doSecond      = data->doSecond;
    double       *vList         = data->vList;
    int           time          = data->time;
    if(doSecond)
    {
      for(int i=0;i<data->length;i++)
      {
        pdfUpdateList[i]->getLookaheadLogP(vList,time,true);
      }
    }
    else
    {
      for(int i=0;i<data->length;i++)
      {
        *resultList[i] = pdfUpdateList[i]->getLookaheadLogP(vList,time,false);
      }
    }
#endif
    pthread_mutex_lock( &condition_mutexDone);
    nrDataThreadsFinished++;
    pthread_cond_signal( &condition_threadDone);
    pthread_mutex_unlock( &condition_mutexDone);
  }
}

#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor only initialises some variables.
/////////////////////////////////////////////////////////////////////////////////////////////////////

LexicalTree::LexicalTree(FILE *outFile)
{
  lookAheadUnigram = NULL;
  vList          = NULL;
  borrowWordList = NULL;
  fastCompressedTree = NULL;
  nodeArray      = NULL;
  nodeListLength = 0;
  outputFile     = outFile;

  overwritePrunePars(PERFORM_PRUNING_HISTOGRAM,PERFORM_PRUNING_STATE,PERFORM_PRUNING_ENDSTATE,
                     BEAM, STATE_BEAM, END_STATE_BEAM, LMLA, HISTOGRAM_STATE_PRUNING, HISTOGRAM_PRUNING);
  overwriteWeightPars(LM_SCALE,TRANS_PENALTY,SIL_PENALTY);


  for(int i=0;i<3;i++)
  {
    totCategories[i]  = 0;
  }
  for(int i=0;i<NUMBER_OF_BLAME_CLUSTERS;i++)
  {
    blameAssingment[i]       = 0;
    analyse_deletions[i]     = 0;
    analyse_insertions[i]    = 0;
    analyse_substitutions[i] = 0;
  }  
  analyse_totRefWords = 0;

  startTime           = 0;
  latticeGeneration   = false;
  phoneLoopConfidence = NULL;
  phoneLoopConfidenceOffset = 0;
  grammarStartContext = 0;
  currentlyAligning   = false;
  startOfSentenceWord = 0;
  endOfSentenceWord   = 0;
  endNode             = NULL;
  languageModel       = NULL;
  nodeList            = NULL;
  timeStamp           = 0;
  bestL               = SMALLEST_NEGATIVE;
  numberOfPhones      = -1;
  numberOfWords       = -1;
  transitionPenalty   = NULL;
  treeStart           = NULL;
  treeEnd             = NULL;
  grammarStart        = NULL;
  wlrStart            = NULL;
  wlrNBest            = NULL;
  latticeAdmin        = NULL;
  vocabulary          = NULL;
  wordLength          = NULL;
  analysisSettings    = NULL;
  tdFile              = NULL;
  tokenDepthAdmin     = NULL;
  latticeWordList     = NULL;
  initialLMHist[0]    = startOfSentenceWord;
  for(int i=1;i<LM_NGRAM_DEPTH;i++)
  {
    initialLMHist[i]  = -1;
  }
  
  for(int i=0;i<LMLA_CACHE;i++)
  {
    lmla_list[i] = NULL;
  }
#ifdef MAX_NR_THREADS
  for(int i=0;i<MAX_NR_THREADS;i++)
  {
    threadMayStart[i]         = false;
    dataThread[i].threadID    = i;
    pthread_create( &thread[i], NULL, thread_pdfCalculation, (void*) &dataThread[i]);
  }
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This constructor reads a binary PPT file created with Shout_dct2lextree. It creates the lexical 
/// nodes and initialises all variables.
/////////////////////////////////////////////////////////////////////////////////////////////////////

LexicalTree::LexicalTree(FILE *outFile, FILE *treeFile, bool useT)
{
  lookAheadUnigram = NULL;
  vList          = NULL;
  borrowWordList = NULL;
  fastCompressedTree          = NULL;
  nodeListLength = 0;
  nodeArray      = NULL;
  outputFile     = outFile;

  overwritePrunePars(PERFORM_PRUNING_HISTOGRAM,PERFORM_PRUNING_STATE,PERFORM_PRUNING_ENDSTATE,
                     BEAM, STATE_BEAM, END_STATE_BEAM, LMLA, HISTOGRAM_STATE_PRUNING, HISTOGRAM_PRUNING);
  overwriteWeightPars(LM_SCALE,TRANS_PENALTY,SIL_PENALTY);


  for(int i=0;i<3;i++)
  {
    totCategories[i]  = 0;     
  }
  for(int i=0;i<NUMBER_OF_BLAME_CLUSTERS;i++)
  {
    blameAssingment[i]       = 0;
    analyse_deletions[i]     = 0;
    analyse_insertions[i]    = 0;
    analyse_substitutions[i] = 0;
  }
  analyse_totRefWords = 0;

  latticeGeneration    = false;
  startTime            = 0;
  phoneLoopConfidence  = NULL;
  phoneLoopConfidenceOffset = 0;
  grammarStartContext  = 0;
  currentlyAligning    = false;
  endNode              = NULL;
  languageModel        = NULL;
  nodeList             = NULL;
  timeStamp            = 0;
  bestL                = SMALLEST_NEGATIVE;
  treeStart            = NULL;
  grammarStart         = NULL;
  wlrStart             = NULL;
  wlrNBest             = NULL;
  latticeAdmin         = NULL;
  analysisSettings     = NULL;
  tdFile               = NULL;
  tokenDepthAdmin      = NULL;
  latticeWordList      = NULL;
  startOfSentenceWord  = -1;
  endOfSentenceWord    = -1;
  freadEndianSafe(&numberOfPhones,1,sizeof(numberOfPhones),treeFile);
  freadEndianSafe(&numberOfWords,1,sizeof(numberOfWords),treeFile);

  vocabulary = new              char*[numberOfWords];
  wordLength = new                int[numberOfWords];
  for(int i=0;i<numberOfWords;i++)
  {
    int len;
    freadEndianSafe(&len,1,sizeof(len),treeFile);
    vocabulary[i] = new char[len+1];
    wordLength[i] = 0;
    freadEndianSafe(vocabulary[i],len,1,treeFile);
    vocabulary[i][len] = 0;
    if(strcmp(vocabulary[i],START_OF_SENTENCE)==0)
    {
      startOfSentenceWord = i;
    }
    else if(strcmp(vocabulary[i],END_OF_SENTENCE)==0)
    {
      endOfSentenceWord = i;
    }
  }
  if(startOfSentenceWord == -1)
  {
    USER_WARNING("The dictionary does not contain a start of sentence word (\"<s>\"). Mapping it to the first word..");
    startOfSentenceWord = 0;
  }
  if(endOfSentenceWord == -1)
  {
    USER_WARNING("The dictionary does not contain an end of sentence word (\"</s>\"). Mapping it to the first word..");
    endOfSentenceWord = 0;
  }

  treeStart = new LexicalNode*[numberOfPhones*numberOfPhones];
  treeEnd   = new LexicalNode*[numberOfPhones*numberOfPhones];
  for(int a=0;a<numberOfPhones*numberOfPhones;a++)
  {
    treeStart[a] = NULL;
    treeEnd[a]   = NULL;
  }

  // Read the model from file:
  for(int a=0;a<numberOfPhones;a++)
  {
    char treeExists;
    freadEndianSafe(&treeExists,1,sizeof(treeExists),treeFile);
    if(treeExists == 1)
    {
      treeStart[a] = new LexicalNode;
      treeStart[a]->inputToken          = new TokenType*;
      (*treeStart[a]->inputToken)       = NULL;
      treeStart[a]->tokenSeqLength      = NONSIL_NUMBER_OF_STATES;
      treeStart[a]->tokenSeq            = new TokenType*[NONSIL_NUMBER_OF_STATES];
      for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
      {
        treeStart[a]->tokenSeq[ai]      = NULL;
      }
      treeStart[a]->nextTree            = NULL;
      treeStart[a]->nodeIsActive        = false;
      treeStart[a]->toBeDeletedFromList = false;
      readTree(treeStart[a], treeFile, 1);
    }
  }

  // The structure is set, handle all 
  // context-dependent stuff etc: 
  setTreeStartEndMatrix();

  // Pre-calculate some stuff:
  transitionPenalty = new float[numberOfWords];
  for(int a=1;a<numberOfWords;a++)
  {
    transitionPenalty[a] = -settings.weights_TransPenalty;
    if(wordLength[a] <= SHORT_WORD_LENGTH)
    {
      transitionPenalty[a] -= SHORT_WORD_PENALTY;
    }
  }
#ifdef MAX_NR_THREADS
  for(int i=0;i<MAX_NR_THREADS;i++)
  {
    threadMayStart[i]         = false;
    dataThread[i].threadID    = i;
    pthread_create( &thread[i], NULL, thread_pdfCalculation, (void*) &dataThread[i]);
  }
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Determines if we should create administration for lattices...
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::setLatticeGeneration(bool setting)
{
  latticeGeneration = setting;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Sets the weight decoder settings for this tree.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::overwriteWeightPars(double lmScale, double transPenalty, double silPenalty)
{
  settings.weights_LmScale      = lmScale;
  settings.weights_TransPenalty = transPenalty;
  settings.weights_SilPenalty   = silPenalty;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Sets the pruning decoder settings for this tree.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::overwritePrunePars(bool doHist, bool doBeam, bool doEndState,
                                     double beam, double state_beam, double endstate_beam,bool lmla,int histState, int hist)
{
  assert((doHist && doBeam) || (!doHist));
  
  settings.doHist             = doBeam;
  settings.doEndState         = doEndState;
  settings.doBeam             = doBeam;
  settings.prune_Beam         = beam;
  settings.prune_StateBeam    = state_beam;
  settings.prune_EndStateBeam = endstate_beam;
  settings.prune_Lmla         = lmla;
  settings.prune_HistState    = histState;
  settings.prune_Hist         = hist;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method is added to check if all phones used in this lexical tree are robust: if a good amount
/// of training occurences were used during creation. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::checkTreeRobustness(LexicalNode *node)
{
  if(node->next != NULL)
  {
    checkTreeRobustness(node->next);
  }
  if(node->parallel != NULL)
  {
    checkTreeRobustness(node->parallel);
  }
  assert(node->modelID < numberOfPhones);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method recursively reads a binary PPT file.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::readTree(LexicalNode *node, FILE *treeFile, int length)
{
  char state;
  freadEndianSafe(&state,1,sizeof(state),treeFile);
  if(state == 0)
  {
    LexicalNode *n = new LexicalNode;
    n->inputToken          = new (TokenType*);
    (*n->inputToken)       = NULL;
    n->tokenSeqLength      = NONSIL_NUMBER_OF_STATES;
    n->tokenSeq            = new TokenType*[NONSIL_NUMBER_OF_STATES];
    for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
    {
      n->tokenSeq[ai]      = NULL;
    }
    n->nextTree            = NULL;
    n->nodeIsActive        = false;
    n->toBeDeletedFromList = false;
    node->next = n;
    readTree(n,treeFile,length+1);
  }
  else
  {
    node->next = NULL;
  }
  freadEndianSafe(&state,1,sizeof(state),treeFile);
  if(state == 1)
  {
    LexicalNode *n = new LexicalNode;
    n->inputToken          = node->inputToken;
    n->tokenSeqLength      = NONSIL_NUMBER_OF_STATES;
    n->tokenSeq            = new TokenType*[NONSIL_NUMBER_OF_STATES];
    for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
    {
      n->tokenSeq[ai]      = NULL;
    }
    n->nextTree            = NULL;
    n->nodeIsActive        = false;
    n->toBeDeletedFromList = false;
    node->parallel = n;
    readTree(n,treeFile,length);
  }
  else
  {
    node->parallel = NULL;
  }
  freadEndianSafe(&state,1,sizeof(state),treeFile);
  
  assert(state == 2);
  freadEndianSafe(&(node->modelID),1,sizeof(node->modelID),treeFile);
  freadEndianSafe(&(node->contextPrev),1,sizeof(node->contextPrev),treeFile);
  freadEndianSafe(&(node->contextNext),1,sizeof(node->contextNext),treeFile);
  freadEndianSafe(&(node->wordID),1,sizeof(node->wordID),treeFile);
  if(node->wordID >= 0)
  {
    wordLength[node->wordID] = length;
    node->nextTree           = NULL;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will change the depth of all nodes to depth.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::setDepthLevel(LexicalNode *node, int phone, int depth)
{
  node->depth = phone + depth*numberOfPhones;
    
  if(node->parallel != NULL)
  {
    setDepthLevel(node->parallel,phone,depth);
  }
  if(node->next != NULL)
  {
    if(depth<TOKENDISTRIBUTION_TREEDEPTH)
    {
      depth++;
    }
    setDepthLevel(node->next,phone,depth);
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::prepareLMLACreation(LexicalNode *node)
{
  if(node->next != NULL)
  {
    prepareLMLACreation(node->next);
  }
  if(node->parallel != NULL)
  {
    prepareLMLACreation(node->parallel);
  }

  if(node->compressedTreeIndex >= 0)
  {
    fastCompressedTree[node->compressedTreeIndex].wordID = node->wordID;
    if((node->wordID >= 0) && (node->modelID == 0))
    {
      fastCompressedTree[node->compressedTreeIndex].wordID = -10;
    }
    LexicalNode *l = node->next;
    while(l != NULL && l->compressedTreeIndex < 0)
    {
      l=l->next;
    }
    if(l != NULL)
    {
      fastCompressedTree[node->compressedTreeIndex].nextIndex = l->compressedTreeIndex;
    }
    else
    {
      fastCompressedTree[node->compressedTreeIndex].nextIndex = -1;
    }
    fastCompressedTree[node->compressedTreeIndex].hasPar = (node->parallel != NULL);
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Each node in the tree contains a number of parameters that give information on the location of the
/// node in the tree. These parameters are set with this method.
///
/// The most important parameter is compressedTreeIndex. This parameter defines if a node is also
/// part of the 'compressed tree'.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::setNodeLocationPars(LexicalNode *node, bool fromParallel)
{
  // Okay, first set the context key:
  node->contextKey = node->contextPrev*numberOfPhones + node->contextNext;

  if(node->next != NULL)
  {
    setNodeLocationPars(node->next,false);
  }

  if(node->parallel != NULL)
  {
    setNodeLocationPars(node->parallel,true);
  }

  if(fromParallel || node->wordID >= 0 || node->parallel != NULL)
  {
    node->compressedTreeIndex = numberOfCompressedNodes;
    numberOfCompressedNodes++;
  }
  else
  {
    node->compressedTreeIndex = -1;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method deletes all LM lookahead data (the global lookahead list).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::deleteLookAheadList()
{
  if(lookAheadUnigram != NULL)
  {
    delete[] lookAheadUnigram->lookAhead;
    delete lookAheadUnigram;
    lookAheadUnigram = NULL;
  }
 
  for(int i=0;i<LMLA_CACHE;i++)
  {
    LMLAGlobalListType *l = lmla_list[i];
    LMLAGlobalListType *l2;
    while(l != NULL)
    {
      l2 = l->next;
#ifndef LMLA_UNIGRAM_ONLY
      if(l->lookAhead != NULL)
      {
        delete[] l->lookAhead;
      }
#endif
      delete l;
      l = l2;
    }
    lmla_list[i] = NULL;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method deletes the lexical tree (PPT) and all its nodes.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::deleteTree(LexicalNode *node, bool isPar)
{
  if(node != NULL)
  {
    if(node->depth >= 0)
    {
      deleteTree(node->next,false);
      deleteTree(node->nextTree,false);
    }
    deleteTree(node->parallel,true);
    if(!isPar)
    {
      delete node->inputToken;
      COUNTER_MINUS(tokCount);
    }
    delete[] node->tokenSeq;  // Tokens itself are already deleted..
    delete node;
    COUNTER_MINUS(gaussianCount);    
    node = NULL;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor will delete the tree, remaining tokens and LMLA tables.
/////////////////////////////////////////////////////////////////////////////////////////////////////

LexicalTree::~LexicalTree()
{
  if(borrowWordList != NULL)
  {
    for(int i=0;i<numberOfWords;i++)
    {
      if(borrowWordList[i] != NULL)
      {
        deleteTree(borrowWordList[i],false);
      }
    }
    delete[] borrowWordList;
  }

#ifdef MAX_NR_THREADS
  for(int i=0;i<MAX_NR_THREADS;i++)
  {
    dataThread[i].length = 0;
    pthread_mutex_lock( &condition_mutexStart[i] );
    threadMayStart[i] = true;
    pthread_cond_signal( &condition_threadStart[i] );
    pthread_mutex_unlock( &condition_mutexStart[i] );
    pthread_join( thread[i], NULL);
  }
#endif

  if(fastCompressedTree != NULL)
  {
    delete[] fastCompressedTree;
    fastCompressedTree = NULL;
  }
  if(nodeArray != NULL)
  {
    delete[] nodeArray;
    nodeArray      = NULL;
    nodeListLength = 0;
  }
  while(nodeList != NULL)
  {
    LexicalNodeList *l = nodeList;
    nodeList = nodeList->next;
    delete l;
  }
  nodeList = NULL;

  if(transitionPenalty != NULL)
  {
    delete[] transitionPenalty;
  }
  deleteLookAheadList();
  if(grammarStart != NULL)
  {
    if(latticeWordList != NULL)
    {
      initialiseSystem();
      LexicalNode *n2 = grammarStart->parallel;
      LexicalNode *n3 = grammarStart;
      while(n2 != NULL)
      {
        if(n2->depth != 2 && n2->depth != -2)
        {
          n2->inputToken    = new (TokenType*);
          *n2->inputToken   = NULL;
          n3->parallel = NULL;
        }
        else
        {
          n2->inputToken = n3->inputToken;
        }
        n3 = n2;
        n2 = n2->parallel;
      }

      grammarStart->parallel = NULL;
      for(int i=0;i<latticeWordListLength;i++)
      {
        deleteTree(latticeWordList[i],false);
      }

      delete[] latticeWordList;
      latticeWordList = NULL;
    }
    else
    {
      deleteTree(grammarStart,false);
      grammarStart      = NULL;
      currentlyAligning = true;
      initialiseSystem();
    }
    grammarStart          = NULL;
    endNode               = NULL;
    latticeWordListLength = 0;
    currentlyAligning     = false;
  }

  if(treeStart != NULL)
  {
    initialiseSystem();
    for(int a=0;a<numberOfPhones;a++)
    {
      if(treeStart[a] != NULL)
      {
        deleteTree(treeStart[a],false);
        treeStart[a] = NULL;
      }
    }

    for(int a=0;a<numberOfPhones;a++)
    {
      for(int a2=0;a2<numberOfPhones;a2++)
      {
        int index  = a*numberOfPhones+a2;
        int index2 = a2*numberOfPhones+a;
        LexicalNode *n = treeStart[index];
        if(n!=NULL)
        {
          delete n->inputToken;
          n->inputToken = NULL;
          COUNTER_MINUS(tokCount);
        }
        while(n!=NULL)
        {
          LexicalNode *n2 = n->parallel;
          delete[] n->tokenSeq;  // Tokens itself are already deleted..
          delete n;
          n = n2;
        }
        
        n = treeEnd[index2];
        if(n!=NULL)
        {
          delete n->inputToken;
          n->inputToken = NULL;
          COUNTER_MINUS(tokCount);
        }
        while(n!=NULL)
        {
          LexicalNode *n2 = n->parallel;
          delete[] n->tokenSeq;  // Tokens itself are already deleted..
          delete n;
          n = n2;
        }      
      }    
    }
    delete[] treeStart;
    delete[] treeEnd;
  }

  for(int i=0;i<numberOfWords;i++)
  {
    if(vocabulary[i] != NULL)
    {
      delete(vocabulary[i]);
    }
  }
  delete[] vocabulary;  
  if(analysisSettings != NULL)
  {
    if(analysisSettings->prunePathAnalyses != NULL)
    {
      for(int i=0;i<analysisSettings->prunePathAnalysesLength;i++)
      {
        if(doPhoneAlignment)
        {
          PhoneModel::initialisePhonePath(analysisSettings->prunePathAnalyses[i]->w.phoneAlignment);
        }
        delete analysisSettings->prunePathAnalyses[i];
        COUNTER_MINUS(trackCount);
      }
      delete[] analysisSettings->prunePathAnalyses;
    }
    if(analysisSettings->binDistr != NULL)
    {
      delete[] analysisSettings->binDistr;
      delete[] analysisSettings->bestScore;
      delete[] analysisSettings->correctScore;
      delete[] analysisSettings->bestStateScore;
      delete[] analysisSettings->stateRanking;
    }
    delete analysisSettings;
    analysisSettings = NULL;
  }  

  if(tokenDepthAdmin != NULL)
  {
    delete[] tokenDepthAdmin;
    tokenDepthAdmin = NULL;
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The application can add a LanguageModel to the lexical tree system. If a model is added, when 
/// tokens are re-propagated to the root of the tree, the language model probability is calculated and 
/// added to the likelihood of the token. Also if the compiler switch LMLA is set, this LM will be used
/// to calculate LMLA tables.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::setLM(LanguageModel *lm)
{
  if(lm == NULL)
  {
    lm = new LanguageModel(numberOfWords);
  }
  languageModel = lm;
  if(languageModel->getNumberOfWords() < numberOfWords)
  {
    USER_ERROR("The language model has less words than the dictionary! Something is wrong!");
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Uses the AM to make a phone-loop lexical tree...
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::setPhoneLoop(int nrP, PhoneModel **models)
{
  numberOfPhones = nrP;
  if(phoneModels == NULL && treeStart == NULL)
  {
    grammarStartContext = 0;
    currentlyAligning   = false;
    endNode             = NULL;
    languageModel       = NULL;
    nodeList            = NULL;
    timeStamp           = 0;
    bestL               = SMALLEST_NEGATIVE;
    treeStart           = NULL;
    grammarStart        = NULL;
    wlrStart            = NULL;
    wlrNBest            = NULL;
    latticeAdmin        = NULL;
    numberOfWords = numberOfPhones;
    vocabulary = new     char*[numberOfWords];
    wordLength = new       int[numberOfWords];
    startOfSentenceWord = 0;
    endOfSentenceWord   = 0;    
    for(int i=0;i<numberOfWords;i++)
    {
      vocabulary[i] = new char[11];
      wordLength[i] = 1;
    }    
    treeStart = new LexicalNode*[numberOfPhones*numberOfPhones];
    treeEnd   = new LexicalNode*[numberOfPhones*numberOfPhones];
    for(int a=0;a<numberOfPhones*numberOfPhones;a++)
    {
      treeStart[a] = NULL;
      treeEnd[a]   = NULL;
    }
  
    // Set the phone loop:
    for(int a=0;a<numberOfPhones;a++)
    {
      treeStart[a] = new LexicalNode;
      treeStart[a]->inputToken          = new (TokenType*);
      (*treeStart[a]->inputToken)       = NULL;
      treeStart[a]->tokenSeqLength      = NONSIL_NUMBER_OF_STATES;
      treeStart[a]->tokenSeq            = new TokenType*[NONSIL_NUMBER_OF_STATES];
      for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
      {
        treeStart[a]->tokenSeq[ai]      = NULL;
      }
      treeStart[a]->nextTree            = NULL;
      treeStart[a]->nodeIsActive        = false;
      treeStart[a]->toBeDeletedFromList = false;
      treeStart[a]->next                = NULL;
      treeStart[a]->nextTree            = NULL;
      treeStart[a]->parallel            = NULL;
      treeStart[a]->modelID             = a;
      treeStart[a]->wordID              = a;
      treeStart[a]->contextPrev         = 0;
      treeStart[a]->contextNext         =-1;
    }

    // The structure is set, handle all 
    // context-dependent stuff etc:
    setTreeStartEndMatrix();

    // Set the models so that we can 
    checkAMs(numberOfPhones,models,false);
    for(int i=0;i<numberOfWords;i++)
    {
      strcpy(vocabulary[i],phoneModels[i]->getStatistics()->name);
    }
  }
  else
  {
    USER_ERROR("Can't create a phone-loop lexical tree if a tree or AMs are already loaded!\n");
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The tree skeleton has been constructed. Now finish the start- and end-matrices...
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::setTreeStartEndMatrix()
{
    biggestNodeDepth = numberOfPhones*(TOKENDISTRIBUTION_TREEDEPTH+1);
    numberOfCompressedNodes = 0;
    for(int a=0;a<numberOfPhones;a++)
    {
      if(treeStart[a] != NULL)
      {
        setNodeLocationPars(treeStart[a],false);
        LexicalNode *n = treeStart[a];
        while(n != NULL)
        {
          n->depth = a;
          if(n->next != NULL)
          {
            setDepthLevel(n->next,a,1);
          }
          n = n->parallel;
        }
      }
    }
    assert(fastCompressedTree == NULL);
    fastCompressedTree           = new FastCompressedTree[numberOfCompressedNodes];
    for(int a=0;a<numberOfPhones;a++)
    {
      if(treeStart[a] != NULL)
      {
        prepareLMLACreation(treeStart[a]);
      }
    }

#ifdef MAX_NR_THREADS
  lmla_dataThread.numberOfCompressedNodes        = numberOfCompressedNodes;
  lmla_dataThread.fastCompressedTree             = fastCompressedTree;
  pthread_create( &lmla_thread, NULL, thread_lmlaCalculation, (void*) &lmla_dataThread);
#endif


    // Create the start nodes with help of the template:
    for(int a=0;a<numberOfPhones;a++)        // the model itself
    {
      for(int a2=1;a2<numberOfPhones;a2++)   // contextPrev
      {
        LexicalNode *p = treeStart[a];
        int index = a2*numberOfPhones+a;
        while(p != NULL)
        {      
          LexicalNode *n = new LexicalNode;
          n->compressedTreeIndex = p->compressedTreeIndex;
          n->depth         = p->depth;
          n->toBeDeletedFromList = false;
          n->nodeIsActive  = false;
          n->inputToken    = NULL;
          n->tokenSeqLength=NONSIL_NUMBER_OF_STATES;
          n->tokenSeq      = new TokenType*[NONSIL_NUMBER_OF_STATES];
          for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
          {
            n->tokenSeq[ai]      = NULL;
          }
          n->nextTree      = NULL;
          n->next          = p->next;
          n->contextPrev   = a2;
          n->modelID       = a;        
          n->contextNext   = p->contextNext;
          n->contextKey    = n->contextPrev*numberOfPhones + n->contextNext;
          n->wordID        = p->wordID;
          n->parallel      = treeStart[index];
          treeStart[index] = n;
          p = p->parallel;
        }
        if(treeStart[index] != NULL)
        {
          treeStart[index]->inputToken    = new (TokenType*);
          (*treeStart[index]->inputToken) = NULL;
          p = treeStart[index]->parallel;
          while(p != NULL)
          {
            p->inputToken = treeStart[index]->inputToken;
            p = p->parallel;
          }
        }            
      }    
    }  
      
    for(int a=0;a<numberOfPhones;a++)        // contextPrev
    {
      for(int a2=0;a2<numberOfPhones;a2++)   // The model itself
      {
        int index = a2*numberOfPhones+a;
        for(int a3=0;a3<numberOfPhones;a3++) // contextNext
        {
          if(treeStart[a3+numberOfPhones*a2] != NULL)
          {
            LexicalNode *n = new LexicalNode;
            n->compressedTreeIndex = -1;
            n->toBeDeletedFromList = false;          
            n->nodeIsActive  = false;
            n->inputToken    = NULL;
            n->tokenSeqLength=NONSIL_NUMBER_OF_STATES;
            n->tokenSeq      = new TokenType*[NONSIL_NUMBER_OF_STATES];
            for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
            {
              n->tokenSeq[ai]= NULL;
            }
            n->depth         = biggestNodeDepth + a2;
            n->next          = NULL;
            n->nextTree      = treeStart[a3+numberOfPhones*a2];
            n->contextPrev   = a;
            n->modelID       = a2;        
            n->contextNext   = a3;
            n->contextKey    = n->contextPrev*numberOfPhones + n->contextNext;
            n->wordID        = -1;
            n->parallel      = treeEnd[index];
            treeEnd[index]   = n;
          }
        }   
        if(treeEnd[index] != NULL)
        {
          treeEnd[index]->inputToken    = new (TokenType*);
          (*treeEnd[index]->inputToken) = NULL;
          LexicalNode *p = treeEnd[index]->parallel;
          while(p != NULL)
          {
            p->inputToken = treeEnd[index]->inputToken;
            p->depth      = biggestNodeDepth + a2;
            p = p->parallel;
          }
        }
      }
    }
     
    biggestNodeDepth += numberOfPhones; 
    tokenDepthAdmin = new int[biggestNodeDepth];
    for(int i=0;i<biggestNodeDepth;i++)
    {
      tokenDepthAdmin[i] = 0;
    }   
       
    for(int i=0;i<LMLA_CACHE;i++)
    {
      lmla_list[i] = NULL;
    }
        
    initialLMHist[0] = startOfSentenceWord;
    for(int i=1;i<LM_NGRAM_DEPTH;i++)
    {
      initialLMHist[i]  = -1;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Checks if the AMs are consistent with the lexical tree structure. This method also sets
/// the models as new LexicalTree acoustic models.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::checkAMs(int nrM, PhoneModel **models, bool outXML)
{
  outXML = true;  // Don't print phone info...
  if(nrM != numberOfPhones)
  {
    USER_ERROR("The number of acoustic models does not match the number of models in the lexical tree!\n");
  }
  phoneModels = models;
  for(int a=0;a<numberOfPhones;a++)
  {
    if(treeStart[a] != NULL)
    {
      checkTreeRobustness(treeStart[a]);
    }
  }
  int nrClusters  = 0;
  int nrGaussians = 0;
  for(int a=0;a<numberOfPhones;a++)
  {
    if(!outXML)
    {
      fprintf(outputFile, "## phone %s: Number of clusters: %d, number of gaussians: %d, g/c: %f\n",
        phoneModels[a]->getStatistics()->name,
        phoneModels[a]->getStatistics()->nrOfContexts,
        phoneModels[a]->getNumberOfGaussians(),
       ((float)phoneModels[a]->getNumberOfGaussians()) / ((float)phoneModels[a]->getStatistics()->nrOfContexts));
      fprintf(outputFile, "## ");
    }
    nrClusters  += phoneModels[a]->getStatistics()->nrOfContexts;
    nrGaussians += phoneModels[a]->getNumberOfGaussians();
  }

  if(!outXML)
  {
    fprintf(outputFile, "## TOTAL: Number of clusters: %d, number of gaussians: %d\n",nrClusters,nrGaussians);
  
    for(int a=0;a<numberOfPhones;a++)
    {
      if(phoneModels[a]->getStatistics()->nrOfTrainOcc < MIN_TRAIN_SAMPLES_PER_GAUSSIAN*2)  
      {
        fprintf(outputFile, "## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ##\n");
        fprintf(outputFile, "## WARNING: The AM is not robust! It uses phone '%s' (%d) trained with only %d number of samples!\n",
                        phoneModels[a]->getStatistics()->name,
                        a,
                        phoneModels[a]->getStatistics()->nrOfTrainOcc );
        fprintf(outputFile, "## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ##\n");
      }
    }    
    fprintf(outputFile, "# logP mean per frame of all phones:\n");
  }
  for(int i=0;i<numberOfPhones;i++)
  {
    phoneModels[i]->getStatistics()->frameMeanLikelihood = phoneModels[i]->getStatistics()->likelihood /  
                                                           phoneModels[i]->getStatistics()->nrOfTrainOcc;
    if(!outXML)
    {
      fprintf(outputFile, "# %s (%d) -  %f\n",phoneModels[i]->getStatistics()->name,
                                phoneModels[i]->getStatistics()->nrOfTrainOcc,
                                phoneModels[i]->getStatistics()->frameMeanLikelihood);
    }
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The nodes of the lexical tree each contain modelIDs. These IDs point to acoustical models. 
/// With this method a set of models (in the 'models' array) is added to the system. The index of this
/// model-array corresponds to the modelIDs of the tree.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::setAMs(PhoneModel **models)
{
  phoneModels = models;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will initialise a node and its parallel nodes, following nodes and nodes of following 
/// trees (if any). Initialising a node means deleting all input and internal token lists.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::initialiseNode(LexicalNode *node)
{
  if(node != NULL)
  {
    if(latticeWordList == NULL || node->compressedTreeIndex >= 0)
    {
      if(latticeWordList != NULL)
      {
        node->compressedTreeIndex *= -1;
      }
      node->nodeIsActive         = false;
      node->toBeDeletedFromList  = false;
      initialiseNode(node->parallel);
      if(node->depth >= 0 || latticeWordList != NULL)
      {
        initialiseNode(node->next);
        initialiseNode(node->nextTree);
      }
      PhoneModel::initialiseToken(node->inputToken);
      for(int ai=0;ai<node->tokenSeqLength;ai++)
      {
        PhoneModel::initialiseToken(&(node->tokenSeq[ai]));
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The first token that is fed to the lexical tree root will need an initial language model history.
/// This method changes the initial history to {word1,word2,-1}. This is needed for forced alignment
/// tasks where the initial history is not equal to <sil>. (Used for blame assignment).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::setInitialLMHistory(const char *word1, const char *word2)
{
  initialLMHist[0] = startOfSentenceWord;
  for(int i=1;i<LM_NGRAM_DEPTH;i++)
  {
    initialLMHist[i] = -1;
  }
  int w1 = getWordID(word1);
  int w2 = getWordID(word2);
  if(w1 >= 0)
  {
    initialLMHist[0] = w1;
  }
  if(w2>=0 && languageModel!=NULL)
  {
    languageModel->getP(w2,initialLMHist,initialLMHist);
  }
  else
  {
    w2 = w1;
  }
  grammarStartContext = -1;
  int a = 0;  
  while(a<numberOfPhones && grammarStartContext < 0)
  {
    grammarStartContext = getLastModelForContext(treeStart[a],w2);
    a++;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method initialises the entire tree. The application needs to call this method before 
/// recognition can start. All LMLA tables and token lists are deleted and a single new token is created
/// that is fed to the root of the tree. Also, statistics are reset.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::initialiseTree(int sTime)
{
  startTime = sTime;
  COUNTER_PRINT();
  initialiseSystem();
  COUNTER_PRINT();

  wlrStart                 = new WLRType;
  COUNTER_PLUS(wlrCount);
  wlrStart->nBest          = NULL;
  wlrStart->lattice        = NULL;
  copyLMHistory(wlrStart->lmHistory,initialLMHist);
  wlrStart->isSil          = 1;
  wlrStart->previous       = NULL;
  wlrStart->adminNext      = NULL;
  wlrStart->timeStamp      = startTime-1;
  wlrStart->COMBlikelihood =  0.0;
  wlrStart->LMlikelihood   =  0.0;
  wlrStart->usedAt         =    0;
  wlrStart->phoneAlignment = NULL;
  if(grammarStart == NULL)           // Regular recognition
  {
    assert(grammarStartContext >= 0);
    nodeListLength = 0;
    for(int a=0;a<numberOfPhones;a++)
    {
      if(treeStart[a+grammarStartContext*numberOfPhones] != NULL)
      {
        nodeListLength++;
      }
    }
    assert(nodeArray == NULL);
    nodeArray = new LexicalNode*[nodeListLength];
    nodeListLength = 0;
    for(int a=0;a<numberOfPhones;a++)
    {
      LexicalNode *lexTree = treeStart[a+grammarStartContext*numberOfPhones];
      if(lexTree != NULL)
      {
        nodeArray[nodeListLength]                        = lexTree;
        nodeArray[nodeListLength]->nodeIsActive          = true;
        nodeArray[nodeListLength++]->toBeDeletedFromList = false;
        *lexTree->inputToken = new TokenType;
        COUNTER_PLUS(tokCount);
        if(doPhoneAlignment)
        {
          (*lexTree->inputToken)->phonePath   = new PLRType;
          COUNTER_PLUS(plrCount);
          (*lexTree->inputToken)->phonePath->phoneID    = -1;
          (*lexTree->inputToken)->phonePath->timeStamp  =  startTime;
          (*lexTree->inputToken)->phonePath->likelihood =  0.0;
          (*lexTree->inputToken)->phonePath->stateOffset[0] = 0;  
          (*lexTree->inputToken)->phonePath->stateOffset[1] = 0;  
          (*lexTree->inputToken)->phonePath->previous   =  NULL;    
        }
        (*lexTree->inputToken)->path        = wlrStart;
        (*lexTree->inputToken)->lmLookAhead = NULL;
        (*lexTree->inputToken)->lookAheadV  = 0.0;
        (*lexTree->inputToken)->likelihood  = 0.0;
        (*lexTree->inputToken)->next        = NULL;
      }
    }  
  }
  else                  // Forced alignment
  {
    assert(nodeArray == NULL);
    nodeListLength = 1;
    nodeArray = new LexicalNode*[nodeListLength];
    nodeArray[0]                      = grammarStart;
    nodeArray[0]->nodeIsActive        = true;
    nodeArray[0]->toBeDeletedFromList = false;
    *grammarStart->inputToken = new TokenType;
    COUNTER_PLUS(tokCount);
    if(doPhoneAlignment)
    {
      (*grammarStart->inputToken)->phonePath   = new PLRType;
      COUNTER_PLUS(plrCount);
      (*grammarStart->inputToken)->phonePath->phoneID    = -1;
      (*grammarStart->inputToken)->phonePath->timeStamp  =  startTime;
      (*grammarStart->inputToken)->phonePath->likelihood =  0.0;  
      (*grammarStart->inputToken)->phonePath->stateOffset[0] = 0;  
      (*grammarStart->inputToken)->phonePath->stateOffset[1] = 0;      
      (*grammarStart->inputToken)->phonePath->previous   =  NULL;    
    }
    (*grammarStart->inputToken)->path        = wlrStart;
    (*grammarStart->inputToken)->lookAheadV  = 0.0;
    (*grammarStart->inputToken)->likelihood  = 0.0;
    (*grammarStart->inputToken)->next        = NULL;
    (*grammarStart->inputToken)->lmLookAhead = NULL;

    if(latticeWordList != NULL)   // lattice rescoring
    {
      LexicalNode *ln = grammarStart;
      while(ln != NULL)
      {
        LexicalNode *targetNode = ln->nextTree;
        *targetNode->inputToken = new TokenType;
        COUNTER_PLUS(tokCount);
        if(doPhoneAlignment)
        {
          (*targetNode->inputToken)->phonePath   = new PLRType;
          COUNTER_PLUS(plrCount);
          (*targetNode->inputToken)->phonePath->phoneID    = -1;
          (*targetNode->inputToken)->phonePath->timeStamp  =  startTime;
          (*targetNode->inputToken)->phonePath->likelihood =  0.0;  
          (*targetNode->inputToken)->phonePath->stateOffset[0] = 0;  
          (*targetNode->inputToken)->phonePath->stateOffset[1] = 0;      
          (*targetNode->inputToken)->phonePath->previous   =  NULL;    
        }
        (*targetNode->inputToken)->path        = wlrStart;
        (*targetNode->inputToken)->lmLookAhead = NULL;
        (*targetNode->inputToken)->lookAheadV  = 0.0;
        (*targetNode->inputToken)->likelihood  = 0.0;
        (*targetNode->inputToken)->next        = NULL;
        (*targetNode->inputToken)->lmLookAhead = NULL;

        ln = ln->parallel;
      }
    }
  }
  timeStamp = startTime;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// All pending data will be cleaned up and all parameters are reset to their initial values,
/// so that a new recognition run can start over cleanly.
/////////////////////////////////////////////////////////////////////////////////////////////////////


void LexicalTree::initialiseSystem()
{
  if(treeStart == NULL && grammarStart == NULL)
  {
    // Not aligning or decoding, so all is already initialised!
    return;
  }
  if(nodeArray != NULL)
  {
    delete[] nodeArray;
    nodeListLength = 0;
    nodeArray      = NULL;
  }
  while(nodeList != NULL)
  {
    LexicalNodeList *l = nodeList;
    nodeList = nodeList->next;
    delete l;
  }
  nodeList = NULL;

  bestRecPath = NULL;
#ifdef SEARCH_STATISTICS_ON
  nrOfStateHistPruning               = 0;
  sentenceStats.nrOfHistPruning      = 0;
  sentenceStats.nrOfCleanUps         = 0;
  sentenceStats.nrOfLMLATables       = 0;
  sentenceStats.nrOfLMLATables_Temp  = 0;
  sentenceStats.nrOfActiveNodes      = 0;
  sentenceStats.nrOfActiveTokens     = 0;
  sentenceStats.nrOfLMLATables_Prev  = 0;
  sentenceStats.nrOfActiveNodes_Prev = 0;
  sentenceStats.nrOfActiveTokens_Prev= 0;  
  sentenceStats.nrOfLMLATables_Max   = 0;
  sentenceStats.nrOfActiveNodes_Max  = 0;
  sentenceStats.nrOfActiveTokens_Max = 0;  
#endif  

  WLRTypeList *ww = wlrNBest;
  while(ww != NULL)
  {
    wlrNBest = ww->next;
    delete ww;
    ww       = wlrNBest;
  }
  deleteLatticeAdmin();

  if(grammarStart == NULL)            // Normal recognition
  {
    for(int a=0;a<numberOfPhones;a++)
    {
      if(treeStart[a]!=NULL)
      {
        initialiseNode(treeStart[a]);
      }
    }
    for(int a=0;a<numberOfPhones;a++)
    {
      for(int a2=0;a2<numberOfPhones;a2++)
      {
        int index = a*numberOfPhones+a2;
        LexicalNode *n = treeStart[index];
        while(n!=NULL)
        {
          n->toBeDeletedFromList = false;
          n->nodeIsActive        = false;
          PhoneModel::initialiseToken(n->inputToken);
          for(int ai=0;ai<n->tokenSeqLength;ai++)
          {
            PhoneModel::initialiseToken(&(n->tokenSeq[ai]));
          }
          n = n->parallel;
        }
        n = treeEnd[index];
        while(n!=NULL)
        {
          n->toBeDeletedFromList = false;
          n->nodeIsActive        = false;
          PhoneModel::initialiseToken(n->inputToken);
          for(int ai=0;ai<n->tokenSeqLength;ai++)
          {
            PhoneModel::initialiseToken(&(n->tokenSeq[ai]));
          }
          n = n->parallel;
        }      
      }    
    }
  }
  else                    // Forced alignment
  {
//    latticeBaumWelch_numberNodes(grammarStart);
//    latticeBaumWelch_unmarkNodes(endNode);
    initialiseNode(grammarStart);
    initialiseNode(endNode);
  }
  deleteLookAheadList();

#ifdef SEARCH_STATISTICS_ON
  intervalTimer = 0;  
#else
  if(settings.weights_LmScale)
  {
    intervalTimer = 0;
  }
#endif  

  WLRType *w = wlrStart;
  WLRType *w2;
  while(w != NULL)
  {
    w2 = w;
    if(doPhoneAlignment)
    {
      PhoneModel::initialisePhonePath(w->phoneAlignment);
      w->phoneAlignment = NULL;
    }
    w = w->adminNext;
    delete w2;
    COUNTER_MINUS(wlrCount);    
  }
  wlrStart    = NULL;

  deleteNodes();

  nrOfTokens   = 0;
  bestL        = SMALLEST_NEGATIVE;
  timeStamp    = 0;
}  

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method searches for a word with word ID 'wordID' and returns the last phone model ID.
/// This method is used to find the correct left-context for the first phone in forced alignment.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int LexicalTree::getLastModelForContext(LexicalNode *node, int wordID)
{
  int res = -1;
  if(node != NULL && node->wordID == wordID)
  {
    return node->modelID;
  }
  if(node != NULL && node->next != NULL)
  {
    res = getLastModelForContext(node->next,wordID);
  }
  if(res == -1 && node != NULL && node->parallel != NULL)
  {
    res = getLastModelForContext(node->parallel,wordID);  
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method returns a string of nodes with the phones in it belonging to the word that matches
/// the input string.
/// borrowPhoneString() is used for grammars and forced alignment tasks. It is typically used by
/// a background lexical tree and the result is passed to the tree that does not have this word in 
/// its system. The word is then named 'OOV'.
/////////////////////////////////////////////////////////////////////////////////////////////////////

LexicalNode *LexicalTree::borrowPhoneString(const char *string)
{
  LexicalNode *res = NULL;
  if(borrowWordList == NULL)
  {
    borrowWordList = new LexicalNode*[numberOfWords];
    for(int w=0;w<numberOfWords;w++)
    {
      borrowWordList[w] = NULL;
    }
  }

  int w = getWordID(string);
  if(w >= 0)
  {
    if(borrowWordList[w] == NULL)
    {
      int a = 0;
      while(a<numberOfPhones && borrowWordList[w] == NULL)
      {
        if(treeStart[a] != NULL)
        {
          borrowWordList[w] = getPhoneString(treeStart[a],w);
        }
        a++;
      }
    }
    res = getPhoneString(borrowWordList[w],w);
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method returns a string of nodes with the phones in it belonging to the word with ID 'wordID'.
/// getPhoneString() is used for grammars and forced alignment tasks. It will give the most left and
/// most right phones a context of -1.
/////////////////////////////////////////////////////////////////////////////////////////////////////

LexicalNode *LexicalTree::getPhoneString(LexicalNode *node, int wordID)
{
  LexicalNode *res = NULL;
  LexicalNode *n   = NULL; 
  if(node->wordID == wordID)
  {
    n                = new LexicalNode;
    n->inputToken    = new (TokenType*);
    *n->inputToken   = NULL;
    n->tokenSeqLength= NONSIL_NUMBER_OF_STATES;
    n->tokenSeq      = new TokenType*[NONSIL_NUMBER_OF_STATES];
    for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
    {
      n->tokenSeq[ai]= NULL;
    }
    n->nextTree      = NULL;
    n->modelID       = node->modelID;
    n->contextPrev   = -1; //node->contextPrev;
    n->contextNext   = -1; //node->contextNext;
    n->contextKey    = -1;
    n->next          = NULL;
    n->parallel      = NULL;
    n->wordID        = wordID;
    n->nextTree            = NULL;
    n->depth               = 0;
    n->compressedTreeIndex = 0;
    n->nodeIsActive        = false;
    n->toBeDeletedFromList = false;
  }
  else
  {
    if(node->parallel != NULL)
    {
      res = getPhoneString(node->parallel,wordID);
    }   
    if(res!=NULL)
    {
      return res;
    }    
    if(node->next != NULL)
    {
      res = getPhoneString(node->next,wordID);
    }
    if(res!=NULL)
    {
      n                = new LexicalNode;
      n->inputToken    = new (TokenType*);
      *n->inputToken   = NULL;
      n->tokenSeqLength= NONSIL_NUMBER_OF_STATES;
      n->tokenSeq      = new TokenType*[NONSIL_NUMBER_OF_STATES];
      for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
      {
        n->tokenSeq[ai]= NULL;
      }
      n->modelID       = node->modelID;
      n->contextPrev   = -1; //node->contextPrev;
      n->contextNext   = -1; //node->contextNext;
      n->contextKey    = -1;
      n->next          = res;
      n->parallel      = NULL;
      n->wordID        = -1;
      n->nextTree            = NULL;
      n->depth               = 0;
      n->compressedTreeIndex = 0;
      n->nodeIsActive        = false;
      n->toBeDeletedFromList = false;
    }
  }
  return n;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// FOR FUTURE USE! This message is called when a grammar needs parallel words. Currently this 
/// method is not used. It is only possible to do forced alignment without alternative paths.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::setAlignParallel()
{
  alignParallel = true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// When the lexical tree is in alignment mode and it needs to align an OOV word, this method can
/// be used to pass it a word that is borrowed from another lexical tree (background lexicon) using
/// the method borrowPhoneString(). The LM probabilities will be set to 1.
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool LexicalTree::addForcedAlignOOV(LexicalNode *oovNode)
{
  if(oovNode == NULL || !currentlyAligning)
  {
    return false;
  }
  else
  {
    // Now go to the end of the string and change the wordID to -1,
    // so that we do not match this word to an incorrect ID:
    LexicalNode *n;
    if(endNode == NULL)
    {
      n                = new LexicalNode;
      n->inputToken    = new (TokenType*);
      *n->inputToken   = NULL;
      n->tokenSeqLength= NONSIL_NUMBER_OF_STATES;
      n->tokenSeq      = new TokenType*[NONSIL_NUMBER_OF_STATES];
      for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
      {
        n->tokenSeq[ai]= NULL;
      }
      n->modelID       = 0;
      n->contextPrev   = -1;
      n->contextNext   = -1;
      n->next          = NULL;
      n->nextTree      = NULL;
      n->parallel      = NULL;
      n->wordID        = endOfSentenceWord;    
      n->depth         = MAX_NDEPTH;
      endNode          = n;      
    }
    n = oovNode;
    while(n != NULL)
    {
      n->wordID = -1;
      if(n->next == NULL)
      {
        n->nextTree = endNode;
      }
      n = n->next;
    }    
  }
  return addWordStringToAlignment(oovNode);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool LexicalTree::alignmentIsEmpty()
{
  return (grammarStart == NULL);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The application may call this method repeadingly to make a forced alignment setting. The 'string'
/// variable must be a single word. 
///
/// If it this word is not in the vocabulary, this method will return 
/// false so that the application knows that the sentence (or blame region) has an OOV in it. 
///
/// The nodes for each word will be added to the end of the alignment tree structure. If no alignment
/// tree has been made yet, a new one will be formed. 
///
/// Starting recognition will result in performing a forced alignment, until a NULL 'string' is send.
/// In that case, the alingment tree structure is deleted and the normal PPT is used for recognition.
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool LexicalTree::setForcedAlign(const char *string)
{
  bool res = true;
  if(string == NULL)
  {
    initialiseSystem();
    if(grammarStart != NULL)
    {
      if(latticeWordList != NULL)
      {
        for(int i=0;i<latticeWordListLength;i++)
        {
          deleteTree(latticeWordList[i],false);
        }
        delete[] latticeWordList;
        latticeWordList = NULL;
      }
      else
      {
        deleteTree(grammarStart,false);
      }
    }
    endNode           = NULL;
    grammarStart      = NULL;
    currentlyAligning = false;
    return true;
  }
  else
  {
    currentlyAligning = true;    
    LexicalNode *word = borrowPhoneString(string);


    if((word != NULL) && (word->wordID == startOfSentenceWord))
    {
      while(word != NULL)
      {
        LexicalNode *w2 = word->next;
        delete word;
        word = w2;
      }
      return true;
    }

    if(word != NULL)
    {
      // Add the endNode!
      LexicalNode *n;
      if(endNode == NULL)
      {
        n                = new LexicalNode;
        n->inputToken    = new (TokenType*);
        *n->inputToken   = NULL;
        n->tokenSeqLength= NONSIL_NUMBER_OF_STATES;
        n->tokenSeq      = new TokenType*[NONSIL_NUMBER_OF_STATES];
        for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
        {
          n->tokenSeq[ai]= NULL;
        }
        n->modelID       = 0;
        n->contextPrev   = -1;
        n->contextNext   = -1;
        n->next          = NULL;
        n->nextTree      = NULL;
        n->parallel      = NULL;
        n->wordID        = endOfSentenceWord;
        n->depth         = MAX_NDEPTH;
        endNode          = n;
      }
      n = word;
      while(n->next != NULL)
      {
        n = n->next;
      }
      n->nextTree = endNode;
    }    
    res = addWordStringToAlignment(word);
  }
  return res;
}    


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Calculates the last wordID given a WLRType object.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int LexicalTree::getWordFromWLR(WLRType *wlr)
{
  int wordID = -1;
  if(wlr->isSil == 1)
  {
    wordID = startOfSentenceWord;
  }
  else
  {
    if(languageModel != NULL)
    {
      wordID = languageModel->getLastWordID(wlr->lmHistory);
    }
    else
    {
      wordID = wlr->lmHistory[0];
    }
  } 
  return wordID;
}


 
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This helper-method actually adds a phone-string to the back of the alignment-tree. It is
/// called by setForcedAlign() and addForcedAlignOOV().
/////////////////////////////////////////////////////////////////////////////////////////////////////
   
bool LexicalTree::addWordStringToAlignment(LexicalNode *word)
{
  if(word != NULL)
  {
    if(grammarStart == NULL)
    {
      grammarStart = word;
      LexicalNode *dn = grammarStart;
      while(dn != NULL)
      {
        dn->depth = 0;
        dn = dn->next;
      }
      if(grammarStart == NULL)
      {
        USER_ERROR("No grammar loaded!\n");
      }
      LexicalNode *n   = new LexicalNode;
      n->inputToken    = grammarStart->inputToken;
      n->tokenSeqLength= NONSIL_NUMBER_OF_STATES;
      n->tokenSeq      = new TokenType*[NONSIL_NUMBER_OF_STATES];
      for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
      {
        n->tokenSeq[ai]= NULL;
      }
      n->nextTree      = NULL;
      n->modelID       = 0;    // This is SIL!
      n->next          = NULL;
      n->parallel      = NULL;
      n->contextPrev   = -1;
      n->contextNext   = -1;
      n->contextKey    = -1;
      n->wordID        = endOfSentenceWord;
      n->depth         = 0;
      grammarStart->parallel = n;

      n                = new LexicalNode;
      n->inputToken    = new (TokenType*);
      *n->inputToken   = NULL;
      n->tokenSeqLength=NONSIL_NUMBER_OF_STATES;
      n->tokenSeq      = new TokenType*[NONSIL_NUMBER_OF_STATES];
      for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
      {
        n->tokenSeq[ai]= NULL;
      }
      n->contextPrev   = -2;
      n->contextNext   = -2;
      n->contextKey    = -2;
      n->nextTree      = grammarStart->nextTree;
      n->next          = grammarStart->next;
      n->modelID       = grammarStart->modelID;
      n->parallel      = NULL;
      n->wordID        = grammarStart->wordID;
      n->depth         = -1;
      grammarStart->parallel->nextTree = n;
    }
    else
    {
      LexicalNode *n = grammarStart;
      int counter = 1;
      while(n->next != NULL || (n->nextTree != NULL && n->nextTree != endNode))
      {
        if(n->next != NULL)
        {
          n = n->next;          
        }
        else
        {
          counter++;
          n = n->nextTree; 
        }
      }

      n->nextTree = word;         
      
      LexicalNode *dn = n->nextTree;
      while(dn != NULL)
      {
        dn->depth = counter;
        dn = dn->next;
      }
      n->nextTree->parallel = NULL;
      
      assert(n->nextTree != NULL);
      
      if(n != grammarStart) // \todo: THIS IS A BUG-FIX. Find out why we need it!
      {
      LexicalNode *n2   = new LexicalNode;
      n2->inputToken    = n->inputToken;
      n2->tokenSeqLength=NONSIL_NUMBER_OF_STATES;
      n2->tokenSeq      = new TokenType*[NONSIL_NUMBER_OF_STATES];
      for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
      {
        n2->tokenSeq[ai]= NULL;
      }
      n2->nextTree      = NULL;
      n2->modelID       = n->modelID;
      n2->next          = NULL;
      n2->parallel      = n->parallel;
      n2->wordID        = n->wordID;
      n2->depth         = n->depth;
      n->parallel       = n2;


      n2                = new LexicalNode;
      n2->inputToken    = new (TokenType*);
      *n2->inputToken   = NULL;
      n2->tokenSeqLength=NONSIL_NUMBER_OF_STATES;
      n2->tokenSeq      = new TokenType*[NONSIL_NUMBER_OF_STATES];
      for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
      {
        n2->tokenSeq[ai]= NULL;
      }
      n2->nextTree      = NULL;
      n2->modelID       = 0;    // This is SIL!
      n2->next          = NULL;
      n2->parallel      = NULL;
      n->contextPrev   = -3;
      n->contextNext   = -3;
      n->contextKey    = -3;
      n2->wordID        = endOfSentenceWord;
      n2->depth         = counter;
      n->parallel->nextTree = n2;

      n2                = new LexicalNode;
      n2->inputToken    = new (TokenType*);
      *n2->inputToken   = NULL;
      n2->tokenSeqLength=NONSIL_NUMBER_OF_STATES;
      n2->tokenSeq      = new TokenType*[NONSIL_NUMBER_OF_STATES];
      for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
      {
        n2->tokenSeq[ai]= NULL;
      }
      n2->contextPrev   = -4;
      n2->contextNext   = -4;
      n2->contextKey    = -4;
      n2->nextTree      = NULL;
      n2->next          = n->nextTree->next;
      n2->modelID       = n->nextTree->modelID;
      n2->parallel      = NULL;
      n2->wordID        = n->nextTree->wordID;
      n2->depth         = -1;
      n->parallel->nextTree->nextTree = n2;
      }
    }
  }
  else
  {
    return false;
  }  
  alignParallel = false;
  // setting the parallel nextTree (can't do that the first time. We don't know the next word yet...:
  // This is always NULL except for words of one phone!
  LexicalNode *pr = grammarStart;
  while(pr!=NULL)
  {
    LexicalNode *par = pr->parallel;
    if(par != NULL)
    {
      if(par->nextTree->modelID == 0 && par->nextTree->nextTree != NULL)
      {
        par = par->nextTree->nextTree;
        par->nextTree = pr->nextTree->nextTree;
      }
    }
    if(pr->next != NULL)
    {
      pr = pr->next;
    }
    else
    {
      pr = pr->nextTree;
    }
  }


  // Setting the context:
  assert(grammarStartContext >= 0);
  
  int contextPrev = grammarStartContext;
  if(grammarStart != NULL)
  {
    LexicalNode *pr = grammarStart;
    while(pr!=NULL)
    {

      pr->contextPrev = contextPrev;

      LexicalNode *par = pr->parallel;
      if(par != NULL)
      {
        if(grammarStart != pr)
        {
          par->contextPrev = contextPrev;
          par->contextNext = 0;
          par->contextKey  = contextPrev*numberOfPhones+par->contextNext;
          assert(par->nextTree != NULL);
          par->nextTree->contextPrev = par->modelID;
          par = par->nextTree;             // Now we are SIL!

          par->contextNext = par->nextTree->modelID;
          par->contextKey  = par->contextPrev*numberOfPhones+par->contextNext;
          assert(par->nextTree != NULL);
          par = par->nextTree;          // Now we are the phone after SIL!
          par->contextPrev = 0;
          if(par->next != NULL)
          {
            par->contextNext = par->next->modelID;
          }
          else if(par->nextTree != NULL && par->nextTree->modelID >= 0)
          {
            par->contextNext = par->nextTree->modelID;
          }
          else
          {
            par->contextNext = 0;
          }
          par->contextKey  = par->contextPrev*numberOfPhones+par->contextNext;
        }
        else
        {
          assert(par->nextTree != NULL);
          par->contextPrev = contextPrev;
          par->contextKey  = par->modelID;
          par->contextNext = par->nextTree->modelID;
          par->contextKey  = par->contextPrev*numberOfPhones+par->contextNext;
          par->nextTree->contextPrev = par->modelID;
          par = par->nextTree;
          if(par->next != NULL)
          {
            par->contextNext = par->next->modelID;
          }
          else if(par->nextTree != NULL && par->nextTree->modelID >= 0)
          {
            par->contextNext = par->nextTree->modelID;
          }
          else
          {
            par->contextNext = 0;
          }
          par->contextKey  = par->contextPrev*numberOfPhones+par->contextNext;
        }
      }

      contextPrev = pr->modelID;
      if(pr->nextTree != NULL && pr->nextTree->modelID != -1)
      {
        pr->contextNext = pr->nextTree->modelID;
        pr->contextKey  = pr->contextPrev*numberOfPhones+pr->contextNext;
        pr = pr->nextTree;
      }
      else if(pr->next != NULL)
      {
        pr->contextNext = pr->next->modelID;
        pr->contextKey  = pr->contextPrev*numberOfPhones+pr->contextNext;
        pr = pr->next;
      }
      else
      {
        pr->contextNext = 0;
        pr->contextKey  = pr->contextPrev*numberOfPhones+pr->contextNext;
        pr = NULL;              
      }
    }
  }
  return true;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
// Keeps track of the NBEST WLRS of this timestamp.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::setWlrNBest(WLRType *wlr)
{
  if(wlr->timeStamp >= 0)
  {
    int w2 = getWordFromWLR(wlr);
    int t = 0;
    for(int i=0;i<NBEST_DEPTH;i++)
    {
      if(wlrNBest->nBest[i] != NULL)
      {
        int w1 = getWordFromWLR(wlrNBest->nBest[i]);
        if(w1 == w2)        
        {
          t = 1;
          if(wlrNBest->nBest[i]->COMBlikelihood < wlr->COMBlikelihood)
          {
            wlrNBest->nBest[i] = wlr;
          }
          i = NBEST_DEPTH;
        }
      }
    }
    if(t == 0)
    {
      for(int i=0;i<NBEST_DEPTH;i++)
      {
        if(wlrNBest->nBest[i] == NULL)
        {
          t = i;
          i = NBEST_DEPTH;
        }
        else if(wlrNBest->nBest[i]->COMBlikelihood < wlrNBest->nBest[t]->COMBlikelihood)
        {
          t = i;
        }
      }
      if(wlrNBest->nBest[t] == NULL || wlr->COMBlikelihood > wlrNBest->nBest[t]->COMBlikelihood)
      {
        wlrNBest->nBest[t] = wlr;
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method is called when a token has reached a leaf node. If available, the LM probability will
/// be added. The word history of the token is updated (using a WLRType: Word Link Record) and the 
/// token is passed to the root of the tree.
///
/// When forced alignment is being performed, the token will be passed to the next tree, containing the
/// next word to align.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::processWord(int wordID, TokenType *token, char isSil, LexicalNode *resultNode)
{
//printf("r1 ");fflush(stdout);

  int newHistory[LM_NGRAM_DEPTH];
  for(int i=0;i<LM_NGRAM_DEPTH;i++)
  {
    newHistory[i] = -1;
  }
  float likelihood   = 0.0;
  float lmLikelihood = 0.0;

//printf("r2 ");fflush(stdout);
  if(isSil == 0)
  {
    likelihood   = token->likelihood;
    if(!currentlyAligning)
    {
      likelihood -= token->lookAheadV;  // If LMLA is switched off, the lookAheadV will be zero always..
    }
    lmLikelihood -= settings.weights_TransPenalty;
    if(wordID >= 0)
    {
      if(wordLength != NULL && wordLength[wordID] <= SHORT_WORD_LENGTH)
      {
        lmLikelihood -= SHORT_WORD_PENALTY;
      }
      // 3. Calculate real LM likelihood:
      if(languageModel!=NULL)
      {
        lmLikelihood   += languageModel->getP(wordID,token->path->lmHistory,newHistory)*settings.weights_LmScale;
      }
      else
      {
        newHistory[0] = wordID;
      }
    }
  }
  else
  {
    likelihood   =  token->likelihood;
    if(!currentlyAligning)
    {
      likelihood -= token->lookAheadV; // If LMLA is switched off, the lookAheadV will be zero always..
    }
  }

//printf("r3 ");fflush(stdout);


  WLRType   *wlr = new WLRType;    
  COUNTER_PLUS(wlrCount);    
  wlr->phoneAlignment = NULL;
  wlr->nBest          = NULL;
  wlr->lattice        = NULL;
  wlr->adminNext      = wlrStart;
  wlr->previous       = token->path;
  wlrStart            = wlr;

  if(currentlyAligning)
  {
    wlr->timeStamp    = timeStamp;
    wlr->usedAt       = timeStamp;
  }
  else
  {
    wlr->timeStamp = -1;
    wlr->usedAt    = -1;    
  }

//printf("r4 ");fflush(stdout);

  wlr->COMBlikelihood = (float)(likelihood+lmLikelihood);
  wlr->LMlikelihood   = wlr->previous->LMlikelihood + lmLikelihood;

  TokenType out;
  if(isSil == 0)
  {
    copyLMHistory(wlr->lmHistory, newHistory);
    wlr->isSil        = 0;
  }
  else
  {
    copyLMHistory(wlr->lmHistory, token->path->lmHistory);
    wlr->isSil        = isSil;
  }

//printf("r5 ");fflush(stdout);

  if(!currentlyAligning && latticeGeneration)
  {
    setWlrNBest(wlr);
  }

  if(currentlyAligning)
  {
    out.lmLookAhead = token->lmLookAhead;
  }
  else
  {
    out.lmLookAhead = NULL;
  }
  out.lookAheadV  = 0.0; 
  out.likelihood  = (likelihood+lmLikelihood);
  out.next        = NULL;
  out.path        = wlr;

//printf("r6 ");fflush(stdout);

  if(doPhoneAlignment)
  {
    out.phonePath   = PhoneModel::copyPhonePath(token->phonePath);
    if(!currentlyAligning)
    {
      out.phonePath->likelihood = -lmLikelihood;
    }
  }

  assert(resultNode != NULL);
  PhoneModel::addChain_ordered(resultNode->inputToken,0.0,&out,0,0,&settings,&bestL);

//printf("r7 ");fflush(stdout);

  if(doPhoneAlignment)
  {
    PhoneModel::initialisePhonePath(out.phonePath);
    out.phonePath = NULL;
  }

//printf("r8 ");fflush(stdout);

  if(token->next != NULL)
  {
    processWord(wordID,token->next,isSil,resultNode);
  }
//printf("r? ");fflush(stdout);

}




/////////////////////////////////////////////////////////////////////////////////////////////////////
/// In order to know if a Word Link Record is still being used by any token (the tokens using a particular
/// WLR may be pruned away), this method will "touch" all WLRs that are in the token list 'token'.
///
/// Touching the WLRs every few iterations is faster than checking if all available WLRs are still valid
/// all the time.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::touchWLRs(TokenType *token)
{
  touchWLRpath(token->path);
  if(token->next != NULL)
  {
    touchWLRs(token->next);
  }
}

void LexicalTree::touchWLRpath(WLRType *w)
{
  if((w!=NULL) && (w->usedAt != timeStamp))
  {
    w->usedAt = timeStamp;
    if(w->nBest != NULL)
    {
      w->nBest[NBEST_DEPTH] = ((WLRType*)timeStamp);
      for(int i=0;i<NBEST_DEPTH;i++)
      {
        touchWLRpath(w->nBest[i]);
      }
    }
    touchWLRpath(w->previous);
  }
}

void LexicalTree::createLatticeNodeGroups(WLRType *w, double lmScore2, int wordID2)
{
  if(w!=NULL)
  {
    int wordID = languageModel->getLastWordID(w->lmHistory);
    WLRType *w3 = w->previous;
    if((w3 != NULL) && !(wordID == startOfSentenceWord && w3->previous == NULL))
    {
      LatticeNode *toNode = findLatticeNodes(w);
      w->lattice  = toNode;
      WLRList *wl = new WLRList;
      wl->amScore  = (w->COMBlikelihood - w->LMlikelihood) - (w3->COMBlikelihood - w3->LMlikelihood);
      wl->lmScore  = ((w->LMlikelihood) - (w3->LMlikelihood)) / settings.weights_LmScale;
      wl->totScore = wl->amScore / settings.weights_LmScale + wl->lmScore;
      wl->wlr     = w3;
      wl->next    = toNode->inArcs;
      toNode->inArcs = wl;
    }
    else
    {
      w->lattice = latticeAdmin;
    }
    createLatticeNodeGroups(w3,w->LMlikelihood,wordID);
    if((w3 != NULL) &&  (w3->nBest != NULL))
    {
      for(int i=0;i<NBEST_DEPTH;i++)
      {
        if((w3->nBest[i] != NULL) && (w3->nBest[i] != w3))
        {
          double amLike      = (w->COMBlikelihood - w->LMlikelihood) - (w3->COMBlikelihood - w3->LMlikelihood);
          int dummy[LM_NGRAM_DEPTH];
          double lmLike1     = 0.0;
          double lmLike2     = 0.0;
          if(wordID != startOfSentenceWord)
          {
            lmLike1 = languageModel->getP(wordID,w3->nBest[i]->lmHistory,dummy);
          }
          if(wordID2 != startOfSentenceWord && wordID2 >= 0)
          {
            lmLike2 = languageModel->getP(wordID2,dummy,dummy);
          }
          WLRList *wl = new WLRList;
          wl->amScore  = amLike;
          if(wordID2 < 0 || wordID2 == startOfSentenceWord)
          {
            wl->lmScore  = lmLike1;
          }
          else
          {
            wl->lmScore  = lmLike1; // + lmLike2 - (lmScore2 - w->LMlikelihood)/settings.weights_LmScale;
          }
          wl->totScore = wl->amScore / settings.weights_LmScale + wl->lmScore;
          wl->wlr     = w3->nBest[i];
          wl->next    = w->lattice->inArcs;
          w->lattice->inArcs = wl;

          createLatticeNodeGroups(w3->nBest[i],w->LMlikelihood,wordID);
          w3->nBest[i]       = NULL;
        }
      }
    }
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will add a node to the list of active nodes. This is done when a node has active tokens
/// that need to be processed in later time frames.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::addNodeToList(LexicalNode *node)
{
  if(node != NULL)
  {
    if((node->toBeDeletedFromList))
    {
      // Already in list, so reset the delete trigger:
      node->toBeDeletedFromList = false;
      node->nodeIsActive        = true;
    }
    else if(!node->nodeIsActive)
    {
      // Add to list:
      LexicalNodeList *l        = new LexicalNodeList;
      l->next                   = nodeList;
      l->node                   = node;

      nodeList                  = l;
      node->nodeIsActive        = true;
      node->toBeDeletedFromList = false;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method handles all administration needed when a tokenlist enters the final node of a word
/// (a leaf-node). The Word-Link-Record is updated and the tokens are passed to the next tree.
/// This tree is the phase-in matrix of the lexical tree.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::processNodeOutput(LexicalNode *node)
{
  if((node != NULL) && (!node->toBeDeletedFromList))
  {
    bool killResToken = true;
    TokenType *resToken = NULL;
    phoneModels[node->modelID]->getOutput(node->contextKey, node->tokenSeq, node->tokenSeqLength, &resToken, &settings, &bestL);

    if(doPhoneAlignment)
    {
      if(resToken != NULL && resToken->phonePath != NULL)
      {
        TokenType *rt = resToken;
        while (rt!=NULL) 
        {
          if(rt->phonePath != NULL)
          {
            if(rt->phonePath->likelihood == 0.0)  // Regular phone, not the last one! (or aligning!)
            {
              rt->phonePath->likelihood = (rt->likelihood - rt->lookAheadV); // If LMLA is switched off: lookAheadV = 0...
              if(rt->path != NULL)
              {
                rt->phonePath->likelihood -= rt->path->COMBlikelihood;
              }
            }
            else     // This is the last phone of the word:
            {
              rt->phonePath->likelihood += rt->likelihood;  // no look-ahead
              // A new WLR is already added, so remove the COMBlikelihood of the previous:
              if(rt->path != NULL && rt->path->previous != NULL)
              {
                rt->phonePath->likelihood -= rt->path->previous->COMBlikelihood;
              }
            }
            PLRType *m                = rt->phonePath->previous;
            while(m != NULL)
            {
              rt->phonePath->likelihood -= m->likelihood;
              m = m->previous;
            }
            rt->phonePath->phoneID     = node->modelID;
            rt->phonePath->contextKey  = node->contextKey;
            PLRType *pt         = new PLRType;
            COUNTER_PLUS(plrCount);
            pt->stateOffset[0]  = 0;  
            pt->stateOffset[1]  = 0;      
            pt->likelihood      = 0.0;
            pt->phoneID         = -1;
            pt->previous        = rt->phonePath;
            pt->timeStamp       = timeStamp+1;
            rt->phonePath       = pt;
          }
          rt = rt->next;
        }
      }
    }

    if(node->next != NULL)
    {
      if(resToken != NULL)
      {
        if((*node->next->inputToken) == NULL) // I'm the first one (from an internal or start token)
        {
          (*node->next->inputToken) = resToken;
          killResToken = false;
          addNodeToList(node->next);
        }
        else  // I'm not the first. This must be a start token!
        {
          PhoneModel::addChain(node->next->inputToken,0.0,resToken,0,0,&settings,&bestL,false);
        }

        touchWLRs(resToken); 
      }  
    }
    else if(node->nextTree != NULL && resToken != NULL)
    {
      // This is it: Update the WLR and add the token to the start of the next tree!
      if(currentlyAligning)
      {
        char sil = 0;
        if(node->modelID == 0)
        {
          sil = 1;
        }
        processWord(node->wordID,resToken,sil,node->nextTree);

        if(doPhoneAlignment)
        {
          TokenType *t  = *node->nextTree->inputToken;
          while(t!=NULL)
          {
            if(t->path->phoneAlignment == NULL) // Could be a token from another node (or self-loop)
            {
              t->path->phoneAlignment = t->phonePath;
              PLRType *phonePath = new PLRType;
              COUNTER_PLUS(plrCount);
              phonePath->stateOffset[0] = 0;
              phonePath->stateOffset[1] = 0;
              phonePath->phoneID    = -1;
              phonePath->timeStamp  =  timeStamp+1;
              phonePath->likelihood =  0.0; // I think this was an error, it was: phonePath->likelihood = t->likelihood;  
              phonePath->previous   = NULL;    
              t->phonePath      = phonePath;
            }
            t = t->next;
          }
        }
        addNodeToList(node->nextTree);
      }
      else
      {
        TokenType *t = resToken;
        while(t!=NULL)
        {
          WLRType *w = new WLRType;
          COUNTER_PLUS(wlrCount);
          if(latticeGeneration)
          {
            w->nBest        = wlrNBest->nBest;
          }
          else
          {
            w->nBest        = NULL;
          }
          w->lattice        = NULL;
          w->adminNext      = wlrStart;
          wlrStart          = w;
          w->COMBlikelihood = t->likelihood; 
          w->LMlikelihood   = t->path->LMlikelihood;
          copyLMHistory(w->lmHistory, t->path->lmHistory);
          w->isSil          = t->path->isSil;
          if(latticeGeneration)
          {
            setWlrNBest(w);          
          }
          w->previous       = t->path->previous;
          if(doPhoneAlignment)
          {
            w->phoneAlignment = t->phonePath;
            PLRType *phonePath = new PLRType;
            COUNTER_PLUS(plrCount);
            phonePath->stateOffset[0] = 0;
            phonePath->stateOffset[1] = 0;
            phonePath->phoneID    = -1;
            phonePath->timeStamp  =  timeStamp+1;
            phonePath->likelihood =  0.0;
            phonePath->previous   = NULL;    
            t->phonePath      = phonePath;
          }
          w->timeStamp      = timeStamp;
          t->path           = w;

          t = t->next;
        }
        touchWLRs(resToken);
        if((*node->nextTree->inputToken) == NULL)  // I'm the first...
        {
          (*node->nextTree->inputToken) = resToken;
          killResToken = false;    
          addNodeToList(node->nextTree);
        }
        else  // I'm not the first.
        {
          PhoneModel::addChain(node->nextTree->inputToken,0.0,resToken,0,0,&settings,&bestL,false);
        }
      }
    }


    if(killResToken == true)
    {
      PhoneModel::initialiseToken(&resToken);
    }


    if(!currentlyAligning)
    {
      //// Now touch the WLRs of the internal token-lists: ////
      for(int ai=0;ai<node->tokenSeqLength;ai++)
      {
        if(node->tokenSeq[ai] != NULL)
        {        
          touchWLRs(node->tokenSeq[ai]);
        }
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// All nodes that are in the global active node list will be passed to this method every time frame.
/// The tokens in these nodes will be updated using the acoustic models (and observation v)
/// language model and lookahead values. See paper [XXXX] for more information.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::processNode(LexicalNode *node, Vector *v)
{
  if(node != NULL)
  {
    if(tdFile != NULL && !currentlyAligning)
    {
      tokenDepthAdmin[node->depth]++;
    }

    int nodeIndex = -1;
    if(!currentlyAligning && settings.prune_Lmla)
    {
      nodeIndex = node->compressedTreeIndex;
    }
    phoneModels[node->modelID]->processVector(node->contextKey,v,timeStamp,nodeIndex,node->tokenSeq,node->tokenSeqLength, (*node->inputToken),&settings,&bestL);

    // This node is currently in the list. (otherwise this method will never be called..)
    // We will prune it from the list at the end of out processVector method if needed.
    // We don't want unnecesary swapping...
    if(settings.doBeam)
    {
      int i = 0;
      while((i<node->tokenSeqLength) && (node->tokenSeq[i] == NULL))
      {
        i++;
      }    
      node->toBeDeletedFromList = (i>=node->tokenSeqLength);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// All nodes in the lexical node list 'node' are initialised (active tokens are deleted) and after that
/// the entries in the lexical node list are all deleted (not the nodes itself!)
/// \todo The names of the variables are unclear. Change them!
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::deleteNodes()
{
  for(int i=0;i<nodeListLength;i++)
  {
    LexicalNode *node = nodeArray[i];
    if(node != NULL)
    {
      PhoneModel::initialiseToken(node->inputToken);
      for(int ai=0;ai<node->tokenSeqLength;ai++)
      {
        PhoneModel::initialiseToken(&(node->tokenSeq[ai]));
      }
    }
  }
  delete[] nodeArray;
  nodeArray      = NULL;
  nodeListLength = 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Every time frame, some administrational data needs to be cleaned up.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::processVector_administrationCleanup()
{
///////////////////////////////////////////////// Prune unused WLRs: //////////////////////////  
  // Clean up the NBest WLRs:

  WLRTypeList *ww     = wlrNBest;
  WLRTypeList *wwPrev = NULL;
  while(ww != NULL)
  {
    if(ww->nBest[NBEST_DEPTH] != ((WLRType*)timeStamp))
    {
      if(ww == wlrNBest)
      {
        wlrNBest = ww->next;
        delete ww;
        ww       = wlrNBest;
      }
      else
      {
        wwPrev->next = ww->next;
        delete ww;
        ww           = wwPrev->next;
      }
    }
    else
    {
      wwPrev = ww;
      ww     = ww->next;
    }
  }


  if(!currentlyAligning)
  {
    // Clean up the WLRs:  
    WLRType *w  = wlrStart;
    WLRType *w2 = NULL;
    while(w != NULL)
    {
      if(w->usedAt != timeStamp)
      {
        if(w2 == NULL)
        {
          wlrStart = w->adminNext;
          if(doPhoneAlignment)
          {
            PhoneModel::initialisePhonePath(w->phoneAlignment);
            w->phoneAlignment = NULL;
          }
          delete w;
          COUNTER_MINUS(wlrCount);
          w = wlrStart;
        }
        else
        {
          w2->adminNext = w->adminNext;
          if(doPhoneAlignment)
          {
            PhoneModel::initialisePhonePath(w->phoneAlignment);
            w->phoneAlignment = NULL;
          }

          delete w;
          COUNTER_MINUS(wlrCount);
          w = w2->adminNext;
        }
      }
      else
      {
        w2 = w;
        w = w->adminNext;
      }
    }
  }


/////////////////////// Prune lookahead tables if the interval has passed:   ////////
#ifdef SEARCH_STATISTICS_ON
  intervalTimer++;
  if(intervalTimer >= LMLA_INTERVAL)
  {
    intervalTimer = 0;
    updateStats();
    if(settings.prune_Lmla)
    {
      pruneLMLA();
    }
  }
#else
  if(settings.prune_Lmla)
  {
    intervalTimer++;
    if(intervalTimer >= LMLA_INTERVAL)
    {
      intervalTimer = 0;
      pruneLMLA();
    }  
  }
#endif  

/////////////////////// Update the node list: Delete all 'unused nodes'   ////////
//   printf("Update node list\n");

  if(settings.doBeam)
  {
    for(int i=0;i<nodeListLength;i++)
    {
      LexicalNode *node = nodeArray[i];
      if(node->toBeDeletedFromList)
      {
        node->toBeDeletedFromList = false;
        node->nodeIsActive        = false;
        nodeArray[i]              = NULL;
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// For forced alignment, which uses basically a simple grammar instead of a lexical tree,
/// the last node in the network, the endNode, contains the result of the decoding run.
/// All tokens in this node should be contained. Therefore, the 'usedAt' flags need all to be
/// kept up to date.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::processVector_grammar()
{
  if(endNode != NULL)
  {
    if((*endNode->inputToken) != NULL)
    {
      WLRType *w = (*endNode->inputToken)->path;
      while(w!=NULL)
      {
        w->usedAt = timeStamp;
        w = w->previous;
      }  
    }
    if(endNode->tokenSeq[0] != NULL)
    {
      WLRType *w = endNode->tokenSeq[0]->path;
      while(w!=NULL)
      {
        w->usedAt = timeStamp;
        w = w->previous;
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// All tokens that are in the phase-in matrix of the lexical tree at the point that this method
/// is called, have just been assigned a new language model probability. At this point, a sharper 
/// beam-prune threshold can be used. This method applies this beam. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::processVector_pruneLM()
{
  if(settings.doEndState)
  {
    float bestLMend   = SMALLEST_NEGATIVE;
    for(int a=0;a<numberOfPhones;a++)
    {
      for(int a2=0;a2<numberOfPhones;a2++)
      {
        int index = a+a2*numberOfPhones;
        if(treeStart[index] != NULL && (*treeStart[index]->inputToken) != NULL)
        {
          if((*treeStart[index]->inputToken)->likelihood > bestLMend)
          {
            bestLMend = (*treeStart[index]->inputToken)->likelihood;
          }
        }
      }
    }
    float minimumLMPruneEnd   = bestLMend   - settings.prune_EndStateBeam;
    for(int a=0;a<numberOfPhones;a++)
    {
      for(int a2=0;a2<numberOfPhones;a2++)
      {
        int index = a+a2*numberOfPhones;
        if(treeStart[index] != NULL && (*treeStart[index]->inputToken) != NULL)
        {
          pruneToken(treeStart[index]->inputToken,minimumLMPruneEnd);
          if((*treeStart[index]->inputToken) != NULL)
          {
            touchWLRs((*treeStart[index]->inputToken));
          }
        }
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// First pruning is performed:
/// - global beam pruning
/// - global histogram pruning
/// - state beam pruning (see paper)
/// - state histogram pruning
/// Language Model Beam pruning is performed by processVector_pruneLM()
/// Then the method determines which tokens should be passed to processNodeOutput().
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::processVector_prune_processNodesOutput()
{
  if(!currentlyAligning)
  {
//printf("a");fflush(stdout);
    // Beam pruning will be performed before histogram pruning, so we know that all tokens will be 
    // between bestL and (bestL-BEAM) for histogram pruning!
    // Let's make a couple of bins out of this range and start counting:
    int bins[NR_OF_HISTOGRAM_BINS];
    float binSize = settings.prune_Beam / ((float)(NR_OF_HISTOGRAM_BINS));
    for(int i=0;i<NR_OF_HISTOGRAM_BINS;i++)
    {
      bins[i] = 0;
    }
// printf("b");fflush(stdout);
    for(int loopList=0; loopList < nodeListLength; loopList++)
    {
      LexicalNode *node = nodeArray[loopList];
      // Prune all tokens, either with the global beam or with the state beam:
      if(node->wordID == -1)
      {
        float minLikelihood_0 = bestL - settings.prune_Beam;
        for(int tokenNmbr = 0; tokenNmbr < node->tokenSeqLength; tokenNmbr++)
        {
          if(node->tokenSeq[tokenNmbr] != NULL)
          {
            TokenType *findMin = node->tokenSeq[tokenNmbr];
            float minLikelihood_1 = -9.0e300;
            while(findMin != NULL)
            {
              if(findMin->likelihood > minLikelihood_1)
              {
                minLikelihood_1 = findMin->likelihood;
              }
              findMin = findMin->next;
            }
            minLikelihood_1 -= settings.prune_StateBeam;
            if(minLikelihood_0 >= minLikelihood_1)
            {
              pruneToken(&(node->tokenSeq[tokenNmbr]),minLikelihood_0,binSize, bins);
            }
            else
            {
              pruneToken(&(node->tokenSeq[tokenNmbr]),minLikelihood_1,binSize, bins);
            }
          }
        }
      }
    }
// printf("c");fflush(stdout);

  #ifdef SEARCH_STATISTICS_ON  
      sentenceStats.nrOfHistPruning++;  
  #endif
    int binTotal = 0;
    int maxBin   = -1;
    bool onePlus = false;
    while(maxBin < NR_OF_HISTOGRAM_BINS-1 && binTotal < settings.prune_Hist)
    {
      onePlus = (binTotal > 0);
      maxBin++;
      binTotal+=bins[maxBin];
    }
    maxBin--;
    assert(maxBin > 0 && maxBin <= NR_OF_HISTOGRAM_BINS);
  
  
    bool doHistPrune = onePlus && (binTotal >= settings.prune_Hist);
    float minL       = bestL - (maxBin+1)*binSize;

// printf("d");fflush(stdout);

    for(int loopList=0; loopList < nodeListLength; loopList++)
    {
      LexicalNode *node = nodeArray[loopList];
      if(doHistPrune)
      {
        for(int i=0;i<node->tokenSeqLength;i++)
        {
          if(node->tokenSeq[i] != NULL)
          {
            pruneToken(&(node->tokenSeq[i]), minL);
          }
        }
      }
      if(node->wordID == -1 || currentlyAligning)
      {
// printf("z");fflush(stdout);
        processNodeOutput(node);
// printf("!");fflush(stdout);

      }
    }
// printf("e");fflush(stdout);

  }
  else
  {
    for(int loopList=0; loopList < nodeListLength; loopList++)
    {
      LexicalNode *node = nodeArray[loopList];
      if(node->wordID == -1 || currentlyAligning)
      {
        processNodeOutput(node);
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Helper method to check if two LM histories are equal (all depths) or unequal (all depths).
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool LexicalTree::compareLMHistory(int *lmHistory1, int *lmHistory2)
{
  for(int i=0;i<LM_NGRAM_DEPTH;i++)
  {
    if(lmHistory1[i] != lmHistory2[i])
    {
      return false;
    }
  }  
  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Helper method to copy one LM history administration node to another
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::copyLMHistory(int *lmHistory1, int *lmHistory2)
{
  for(int i=0;i<LM_NGRAM_DEPTH;i++)
  {
    lmHistory1[i] = lmHistory2[i];
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// All tokens from the token list 'token' that have a likelihood less than 'minLikelihood' will
/// be pruned by this method.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::pruneToken(TokenType **token, float minLikelihood, float binSize, int *bins)
{
  TokenType *t = *token;
  TokenType *p = NULL;
  while(t != NULL)
  {
    if(t->likelihood <= minLikelihood)
    {
      // delete node...
      if(doPhoneAlignment)
      {
        PhoneModel::initialisePhonePath(t->phonePath);
      }
      TokenType *del = t;
      t = t->next;
      if(p == NULL)
      {
        *token = t;
      }
      else
      {
        p->next = t;
      }
      delete del;
      COUNTER_MINUS(tokCount);
    }
    else
    {
      if(bins != NULL)
      {
        int binIndex = (int)((bestL - t->likelihood) / binSize);
        if(binIndex >= NR_OF_HISTOGRAM_BINS)
        {
          binIndex = NR_OF_HISTOGRAM_BINS - 1;
        }
        if(binIndex < 0)
        {
          printf("( %f  -  %f ) /  %f = %d\n",bestL,t->likelihood,binSize,binIndex);
        }
        assert(binIndex >= 0);
        bins[binIndex]++;
      }
      p = t;
      t = t->next;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::createActiveNodesList(void)
{
  // Now process the updated list:
  LexicalNodeList *l = nodeList;
  int newListLength = 0;
  while(l != NULL)
  {
    if(l->node->nodeIsActive)
    {
      newListLength++;
    }
    l = l->next;
  }
  for(int i=0;i<nodeListLength;i++)
  {
    if(nodeArray[i] != NULL)
    {
      newListLength++;
    }
  }

  LexicalNode **newList = new LexicalNode*[newListLength];
  newListLength = 0;
  for(int i=0;i<nodeListLength;i++)
  {
    if(nodeArray[i] != NULL)
    {
      newList[newListLength++] = nodeArray[i];
    }
  }
  l = nodeList;
  LexicalNodeList *l2 = NULL;
  while(l != NULL)
  {
    if(l->node->nodeIsActive)
    {
      newList[newListLength++] = l->node;
    }
    l2 = l;
    l = l->next;
    delete l2;
  }
  nodeList = NULL;

  if(nodeArray != NULL)
  {
    delete[] nodeArray;
    nodeArray = NULL;
  }
  nodeArray      = newList;
  nodeListLength = newListLength;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// processVector_processNodes() is part of the inner-loop methods. It provides the observation 
/// Vector v to all nodes that are active. 
///
/// The method consits of multiple steps:
/// - The active nodes list is extended with all nodes that are parralel nodes of a node that 
/// contains input tokens. 
/// - The 'fan out' nodes, or end-nodes that contain input tokens are added to the list.
/// - All nodes are processed except for the nodes that contain wordIDs. The tokens for these
/// nodes are already passed to the endNodes.
/// - All inputToken lists are deleted (they will be replaced by the new inputToken lists by
/// the method processVector_prune_processNodesOutput() ).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::processVector_processNodes(Vector **v)
{
  // First add all parallel nodes to the list when the input is not NULL, and copy 
  // the input list to the end nodes if needed:

  for(int i=0;i<nodeListLength;i++)
  {
    LexicalNode *node = nodeArray[i];
    if(node != NULL && (*node->inputToken) != NULL)
    {
      if((node->wordID >= 0) && !currentlyAligning)
      {
        char sil = 0;
        if(node->modelID == 0)
        {
          sil = 1;
        }
        processWord(node->wordID,(*node->inputToken),sil,
                    treeEnd[node->modelID*numberOfPhones+node->contextPrev]);
      }
      LexicalNode *ln = node->parallel;
      while(ln != NULL && !ln->nodeIsActive)
      {
        addNodeToList(ln);
        if((ln->wordID >= 0) && !currentlyAligning)
        {
          char sil = 0;
          if(node->modelID == 0)
          {
            sil = 1;
          }
          processWord(ln->wordID,(*ln->inputToken),sil,
                      treeEnd[ln->modelID*numberOfPhones+ln->contextPrev]);
        }
        ln = ln->parallel;
      }
    }
  }

//  printf("2!\n");fflush(stdout);

  if(!currentlyAligning && settings.doEndState)
  {
    float bestLMend   = SMALLEST_NEGATIVE;
    for(int a=0;a<numberOfPhones;a++)
    {
      for(int a2=0;a2<numberOfPhones;a2++)
      {
        int index = a+a2*numberOfPhones;
        if(treeEnd[index] != NULL && (*treeEnd[index]->inputToken) != NULL)
        {
          if((*treeEnd[index]->inputToken)->likelihood > bestLMend)
          {
            bestLMend = (*treeEnd[index]->inputToken)->likelihood;
          }
        }
      }
    }
    float minimumLMPruneEnd   = bestLMend   - settings.prune_StateBeam;

    for(int a=0;a<numberOfPhones;a++)
    {
      for(int a2=0;a2<numberOfPhones;a2++)
      {
        int index = a+a2*numberOfPhones;
        if(treeEnd[index] != NULL && (*treeEnd[index]->inputToken) != NULL)
        {
          pruneToken(treeEnd[index]->inputToken,minimumLMPruneEnd);
          if((*treeEnd[index]->inputToken) != NULL)
          {
            touchWLRs((*treeEnd[index]->inputToken));
          }
        }
      }
    }
  }

//  printf("3!\n");fflush(stdout);

  if(settings.prune_Lmla)
  {
    processVector_LMLAReordering_prepare();
  }

//  printf("4!\n");fflush(stdout);

  if(grammarStart == NULL)
  {
    // Add all end nodes to the list that have input:
    for(int a=0;a<numberOfPhones;a++)
    {
      for(int a2=0;a2<numberOfPhones;a2++)
      {
        int index = a+numberOfPhones*a2;
        if(treeEnd[index] != NULL && (*treeEnd[index]->inputToken) != NULL)
        {
          LexicalNode *ln = treeEnd[index];
          while(ln != NULL)
          {
            addNodeToList(ln);
            ln = ln->parallel;
          }
        }
      }
    }
  }

//  printf("5!\n");fflush(stdout);

  createActiveNodesList();

//  printf("6!\n");fflush(stdout);

#ifdef MAX_NR_THREADS
  if(vList != NULL)
  {
    pthread_mutex_lock( &condition_mutexDone );
    while(nrDataThreadsFinished < threadsRunning)
    {
      pthread_cond_wait( &condition_threadDone, &condition_mutexDone );
    }
    nrDataThreadsFinished = 0;
    pthread_mutex_unlock( &condition_mutexDone );

    threadsRunning = 0;
    delete[] vList;
    vList = NULL;
    delete[] pdfUpdateList;
    pdfUpdateList = NULL;
    COUNTER_MINUS(pdfCount);
  }

  // First calculate all gaussian mixtures:
  pdfUpdateList = new MixGaussian*[nodeListLength*3];
  COUNTER_PLUS(pdfCount);

  double      *resultList[nodeListLength*3];
  int cnt = 0;
  for(int i=0;i<nodeListLength;i++)
  {
    if(nodeArray[i] != NULL && (currentlyAligning || nodeArray[i]->wordID == -1))
    {
      LexicalNode *node = nodeArray[i];
      assert(node->modelID >= 0);
      cnt += phoneModels[node->modelID]->touchPDF(node->contextKey,timeStamp,&pdfUpdateList[cnt],&resultList[cnt]);
    }
  }

  if(cnt > 0)
  {
    vList = new double[v[0]->len()*HALF_GMMLOOKAHEAD];
    int cntList = 0;
    for(int i=0;i<v[0]->len();i++)
    {
      for(int i2=0;i2<HALF_GMMLOOKAHEAD;i2++)
      {
        vList[cntList] = v[i2]->getValue(i);
        cntList++;
      }
    }
  
    int totalLength = 0;
    while(totalLength < cnt)
    {
      int dataPerThread = (cnt-totalLength) / MAX_NR_THREADS;
      if(dataPerThread > MAX_PDF_CALC_SIZE)
      {
        dataPerThread = MAX_PDF_CALC_SIZE;
      }
      for(int i=0;i<MAX_NR_THREADS;i++)
      {
        dataThread[i].length      = dataPerThread;
        dataThread[i].mixGaussian = &pdfUpdateList[totalLength];
        dataThread[i].result      = &resultList[totalLength];
  
        dataThread[i].vList       = vList;
        dataThread[i].doSecond    = false;
  
        dataThread[i].time        = timeStamp;
        totalLength              += dataPerThread;
      }
      if((totalLength < cnt) && (dataPerThread + cnt - totalLength  < MAX_PDF_CALC_SIZE))
      {
        dataThread[MAX_NR_THREADS-1].length += cnt - totalLength;
        totalLength = cnt;
      }
  
      // start threads
      threadsRunning = 0;
      for(int i=0;i<MAX_NR_THREADS;i++)
      {
        if(dataThread[i].length > 0)
        {
          threadsRunning++;
          pthread_mutex_lock( &condition_mutexStart[i] );
          threadMayStart[i] = true;
          pthread_cond_signal( &condition_threadStart[i] );
          pthread_mutex_unlock( &condition_mutexStart[i] );
        }
      }
      pthread_mutex_lock( &condition_mutexDone );
      while(nrDataThreadsFinished < threadsRunning)
      {
        pthread_cond_wait( &condition_threadDone, &condition_mutexDone );
      }
      nrDataThreadsFinished = 0;
      pthread_mutex_unlock( &condition_mutexDone );
    }
    delete[] vList;
    vList = NULL;
  
    vList = new double[v[0]->len()*HALF_GMMLOOKAHEAD];
    cntList = 0;
    for(int i=0;i<v[0]->len();i++)
    {
      for(int i2=HALF_GMMLOOKAHEAD;i2<GMMLOOKAHEAD;i2++)
      {
        vList[cntList] = v[i2]->getValue(i);
        cntList++;
      }
    }
  
    totalLength = 0;
    while(totalLength < cnt)
    {
      int dataPerThread = (cnt-totalLength) / MAX_NR_THREADS;
      if(dataPerThread > MAX_PDF_CALC_SIZE)
      {
        dataPerThread = MAX_PDF_CALC_SIZE;
      }
      for(int i=0;i<MAX_NR_THREADS;i++)
      {
        dataThread[i].length      = dataPerThread;
        dataThread[i].mixGaussian = &pdfUpdateList[totalLength];
        dataThread[i].result      = NULL;
  
        dataThread[i].vList       = vList;
        dataThread[i].doSecond    = true;
  
        dataThread[i].time        = timeStamp;
        totalLength              += dataPerThread;
      }
      if((totalLength < cnt) && (dataPerThread + cnt - totalLength  < MAX_PDF_CALC_SIZE))
      {
        dataThread[MAX_NR_THREADS-1].length += cnt - totalLength;
        totalLength = cnt;
      }
  
      // start threads
      int threadsRunning = 0;
      for(int i=0;i<MAX_NR_THREADS;i++)
      {
        if(dataThread[i].length > 0)
        {
          threadsRunning++;
          pthread_mutex_lock( &condition_mutexStart[i] );
          threadMayStart[i] = true;
          pthread_cond_signal( &condition_threadStart[i] );
          pthread_mutex_unlock( &condition_mutexStart[i] );
        }
      }
    }
  }
#endif

//   printf("7!\n");fflush(stdout);

  for(int i=0;i<nodeListLength;i++)
  {
    if(currentlyAligning || nodeArray[i]->wordID == -1)
    {
      processNode(nodeArray[i],v[0]);
    }
    else
    {
      nodeArray[i]->toBeDeletedFromList = true;
    }
  }

//   printf("8!\n");fflush(stdout);

  // Finally, delete all input token lists:
  for(int i=0;i<nodeListLength;i++)
  {
    if((*nodeArray[i]->inputToken) != NULL)
    {
      PhoneModel::initialiseToken(nodeArray[i]->inputToken);
    }
  }

//   printf("89!\n");fflush(stdout);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will provide each input node (fan-in node) from treeStart with an LMLA table pointer.
/// This pointer is either looked-up or created by getLMLATable().
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::processVector_LMLAReordering()
{
  assert(treeStart != NULL);
  // If this method is called, LMLA is switched on!
  for(int a=0;a<numberOfPhones;a++)
  {
    for(int a2=0;a2<numberOfPhones;a2++)
    {
      int index = a+a2*numberOfPhones;
      if(treeStart[index] != NULL && *(treeStart[index]->inputToken) != NULL)
      {
        TokenType *t = *(treeStart[index]->inputToken);
        while(t != NULL)
        {
          t->lookAheadV  = 0.0;
          t->lmLookAhead = getLMLATable(t->path->lmHistory,false);
          t = t->next;
        }
      }
    }
  }
#ifdef MAX_NR_THREADS
    pthread_mutex_lock( &lmla_condition_mutexDone );
    while(!lmla_threadFinished)
    {
      pthread_cond_wait( &lmla_condition_threadDone, &lmla_condition_mutexDone );
    }
    pthread_mutex_unlock( &lmla_condition_mutexDone );
#endif
}


void LexicalTree::processVector_LMLAReordering_prepare()
{
  assert(treeStart != NULL);
  // If this method is called, LMLA is switched on!
  for(int a=0;a<numberOfPhones;a++)
  {
    for(int a2=0;a2<numberOfPhones;a2++)
    {
      int index = a+a2*numberOfPhones;
      if(treeEnd[index] != NULL && *(treeEnd[index]->inputToken) != NULL)
      {
        TokenType *t = *(treeEnd[index]->inputToken);
        while(t != NULL)
        {
          t->lookAheadV  = 0.0;
          t->lmLookAhead = getLMLATable(t->path->lmHistory,true);
          t = t->next;
        }
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Use the last recognition to adapt the Acoustic Models (fill the numerators)
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::adaptAMs(Vector *v, int time)
{ 
  if(doPhoneAlignment)
  {
    if(bestRecPath == NULL)
    {
      TokenType *useToken = findBestToken(false, true, false);
      if(useToken != NULL && useToken->path!=NULL)
      {   
        bestRecPath = useToken->path;
      }    
    }
    if(bestRecPath != NULL)
    {
      WLRType *w = bestRecPath;
      while(w != NULL && w->previous != NULL && w->previous->timeStamp >= time)
      {
        w = w->previous;
      }
      if(w != NULL && w->previous != NULL)   // Meaning that w is now the correct word.
      {
        // Now find the correct phone:      
        PLRType *p = w->phoneAlignment;
        while(p != NULL && p->timeStamp > time)
        {
          p = p->previous;
        }
        if(p != NULL) // Found the correct phone!
        {
          unsigned char offset = (unsigned char)(time - p->timeStamp);
          if(p->stateOffset[0] == 0)  // SIL...
          {
            phoneModels[p->phoneID]->adapt_addAcumulatorData(0,p->contextKey,v);
          }
          else
          {
            if(offset < p->stateOffset[0])
            {
              phoneModels[p->phoneID]->adapt_addAcumulatorData(0,p->contextKey,v);
            }
            else if(offset < p->stateOffset[1])
            {
              phoneModels[p->phoneID]->adapt_addAcumulatorData(1,p->contextKey,v);
            }
            else
            {
              phoneModels[p->phoneID]->adapt_addAcumulatorData(2,p->contextKey,v);
            }
          }
        }    
      }  
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// For the phone-loop: stores the best likelihood minus the number of phones times the transition penalty..
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::storePLConfidence(int time)
{
  assert(phoneLoopConfidence != NULL);
  /*
  if(bestToken != NULL)
  {
    // First find the number of phones:  
    int nrP = 0;
    WLRType *w = bestToken->path;
    while(w != NULL)
    {
      if(w->timeStamp > 0)
      {
        nrP++;
      }
      w = w->previous;
    }
    phoneLoopConfidence[time] = bestL + settings.weights_TransPenalty*nrP;
  }
  else
  {
    phoneLoopConfidence[time] = 0.0;  
  }
*/
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Inner loop method. Every time-frame this method is called with the observation-vector.
/// The tasks in this method are split up to make it possible to analyse efficiency with gprof (profiler).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::processVector(Vector **v, int time)
{
  bestL                  = SMALLEST_NEGATIVE;

  bestTokenInSystem      = NULL;
  timeStamp              = time;
  nrOfTokens             = 0;
  if(currentlyAligning || latticeWordList != NULL)
  {
    processVector_processNodes(v);
    processVector_prune_processNodesOutput();
    processVector_administrationCleanup();
  }
  else
  {
//printf("P");fflush(stdout);
    if(latticeGeneration)
    {
      // NBEST/LATTICE:
      WLRTypeList *ww     = new WLRTypeList;
      ww->nBest           = new WLRType*[NBEST_DEPTH + 1];
      ww->next            = wlrNBest;
      wlrNBest            = ww;
      for(int i=0;i<NBEST_DEPTH+1;i++)
      {
        ww->nBest[i] = NULL;
      }
    }
//printf("O");fflush(stdout);
    if(settings.prune_Lmla)
    {
      processVector_LMLAReordering();
    }

//printf("X");fflush(stdout);
    processVector_processNodes(v);
// printf("Y");fflush(stdout);
    if(tdFile != NULL)
    {
      printTokenDistribution();
    }
// printf("Z");fflush(stdout);
    processVector_prune_processNodesOutput();
// printf("A");fflush(stdout);
    processVector_pruneLM();
// printf("B");fflush(stdout);
    processVector_administrationCleanup();
// printf("C\n");fflush(stdout);
  }
//printf("%d\t%f\n",time,bestL); fflush(stdout);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Prints the token distribution for this time frame..
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::printTokenDistribution()
{
  fprintf(tdFile,"<frame ");
  for(int i=0;i<biggestNodeDepth;i++)
  {
    fprintf(tdFile,"n%04d=\"%d\" ", i,tokenDepthAdmin[i]);
    tokenDepthAdmin[i] = 0;
  }
  fprintf(tdFile, "/>\n");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method prints a welcoming message and initial statistics (such as tree depth etc)
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::printInitialSettings(const char *amName, const char *dctName, const char *backName, const char*lmName, bool outXML)
{
  if(outXML)
  {
    fprintf(outputFile, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(outputFile, "<!--\n");
    fprintf(outputFile, "  ##################################################################################################\n");
    fprintf(outputFile, "  ###              Shout is the decoder of the 'SHoUT LVCS Recognition toolkit'.                 ###\n");
    fprintf(outputFile, "  ###                            This toolkit is developed by:                                   ###\n");
    fprintf(outputFile, "  ###                   Marijn Huijbregts, HMI, University of Twente.                            ###\n");
    fprintf(outputFile, "  ###                      http://wwwhome.cs.utwente.nl/~huijbreg/                               ###\n");
    fprintf(outputFile, "  ###                           marijn.huijbregts@utwente.nl                                     ###\n");
    fprintf(outputFile, "  ##################################################################################################\n");
    fprintf(outputFile, "-->\n");
    fprintf(outputFile, "<shout_metadata>\n");
    fprintf(outputFile, "  <model_info>\n");
    if(amName != NULL && strlen(amName) > 0)
    {
      fprintf(outputFile, "    <AM>%s</AM>\n",amName);
    }
    if(dctName != NULL && strlen(dctName) > 0)
    {
      fprintf(outputFile, "    <DCT>%s</DCT>\n",dctName);
    }
    if(backName != NULL && strlen(backName) > 0)
    {
      fprintf(outputFile, "    <BACK-DCT>%s</BACK-DCT>\n",backName);
    }
    if(lmName != NULL && strlen(lmName) > 0)
    {
      fprintf(outputFile, "    <LM>%s</LM>\n",lmName);
    }
    fprintf(outputFile, "  </model_info>\n");
    fprintf(outputFile, "  <decoding_settings>\n");
    fprintf(outputFile, "    <LM_SCALE>%.1f</LM_SCALE>\n",  settings.weights_LmScale);
    fprintf(outputFile, "    <TRANSITION_PENALTY>%.1f</TRANSITION_PENALTY>\n",  settings.weights_TransPenalty);
    fprintf(outputFile, "    <SHORT_WORD_PENALTY>%.1f</SHORT_WORD_PENALTY>\n",  SHORT_WORD_PENALTY);
    fprintf(outputFile, "    <SHORT_WORD_LENGTH>%d</SHORT_WORD_LENGTH>\n",    SHORT_WORD_LENGTH);  
    fprintf(outputFile, "    <SIL_PENALTY>%.1f</SIL_PENALTY>\n",  settings.weights_SilPenalty);
    fprintf(outputFile, "    <GLOBAL_BEAM>%.1f</GLOBAL_BEAM>\n",  settings.prune_Beam);
    fprintf(outputFile, "    <NODE_BEAM>%.1f</NODE_BEAM>\n",  settings.prune_StateBeam);
    fprintf(outputFile, "    <ENDWORD_BEAM>%.1f</ENDWORD_BEAM>\n",  settings.prune_EndStateBeam);
    if(settings.prune_Lmla)
    {
      fprintf(outputFile, "    <LMLA>ON</LMLA>\n");
      fprintf(outputFile, "    <LMLA_CACHESIZE>%d</LMLA_CACHESIZE>\n",    LMLA_CACHE);
      fprintf(outputFile, "    <LMLA_CLEANUPINTERVAL>%d</LMLA_CLEANUPINTERVAL>\n",    LMLA_INTERVAL);
    }
    else
    {
      fprintf(outputFile, "    <LMLA>OFF</LMLA>\n");
      fprintf(outputFile, "    <LMLA_CLEANUPINTERVAL>%d</LMLA_CLEANUPINTERVAL>\n",    LMLA_INTERVAL);
    }
    fprintf(outputFile, "    <MAX_TOKENS_PER_STATE>%d</MAX_TOKENS_PER_STATE>\n",    settings.prune_HistState);
    fprintf(outputFile, "    <MAX_TOKENS_TOTAL>%d</MAX_TOKENS_TOTAL>\n",    settings.prune_Hist);          
    fprintf(outputFile, "  </decoding_settings>\n");
    fprintf(outputFile, "  <segments>\n");    
  }
  else
  {  
    fprintf(outputFile, "##################################################################################################\n");
    fprintf(outputFile, "###              Shout is the decoder of the 'Shout LVCS Recognition toolkit'.                 ###\n");
    fprintf(outputFile, "###                             Shout is developed by:                                         ###\n");
    fprintf(outputFile, "###                   Marijn Huijbregts, HMI, University of Twente.                            ###\n");
    fprintf(outputFile, "###                      http://wwwhome.cs.utwente.nl/~huijbreg/                               ###\n");
    fprintf(outputFile, "###                           marijn.huijbregts@utwente.nl                                     ###\n");
    fprintf(outputFile, "##################################################################################################\n");
    fprintf(outputFile, "### \n");
    fprintf(outputFile, "### ------------------------------------------------------------------------\n");
    fprintf(outputFile, "###  Models used:                                      \n");
    fprintf(outputFile, "### ------------------------------------------------------------------------\n");
    if(amName != NULL && strlen(amName) > 0)
    {
      fprintf(outputFile, "###  AM:             %s\n",amName);
    }
    if(dctName != NULL && strlen(dctName) > 0)
    {
      fprintf(outputFile, "###  DCT:            %s\n",dctName);
    }
    if(backName != NULL && strlen(backName) > 0)
    {
      fprintf(outputFile, "###  Background DCT: %s\n",backName);
    }
    if(lmName != NULL && strlen(lmName) > 0)
    {
      fprintf(outputFile, "###  LM:             %s\n",lmName);
    }
    fprintf(outputFile, "### ------------------------------------------------------------------------\n");
    fprintf(outputFile, "###  Decoder settings:\n");
    fprintf(outputFile, "### ------------------------------------------------------------------------\n");
    fprintf(outputFile, "###  Language Model scale:     %8.1f\n",  settings.weights_LmScale);
    fprintf(outputFile, "###  Transition penalty:       %8.1f\n",  settings.weights_TransPenalty);
    fprintf(outputFile, "###  Short word penalty:       %8.1f\n",  SHORT_WORD_PENALTY);
    fprintf(outputFile, "###  Short word length:        %8d\n",    SHORT_WORD_LENGTH);  
    fprintf(outputFile, "###  Silence penalty:          %8.1f\n",  settings.weights_SilPenalty);
    fprintf(outputFile, "###  Global beam:              %8.1f\n",  settings.prune_Beam);
    fprintf(outputFile, "###  Node beam:                %8.1f\n",  settings.prune_StateBeam);
    fprintf(outputFile, "###  End words beam:           %8.1f\n",  settings.prune_EndStateBeam);
    if(settings.prune_Lmla)
    {
      fprintf(outputFile, "###  LMLA:                           ON\n");
      fprintf(outputFile, "###  LMLA cache size:          %8d\n",    LMLA_CACHE);
      fprintf(outputFile, "###  LMLA cleanup interval:    %8d\n",    LMLA_INTERVAL);
    }
    else
    {
      fprintf(outputFile, "###  LMLA:                          OFF\n");
      fprintf(outputFile, "###  LMLA cleanup interval:    %8d\n",    LMLA_INTERVAL);
    }
    fprintf(outputFile, "###  Max tokens per state:     %8d\n",    settings.prune_HistState);
    fprintf(outputFile, "###  Max tokens:               %8d\n",    settings.prune_Hist);
    fprintf(outputFile, "### \n");
    fprintf(outputFile, "##################################################################################################\n");
    fprintf(outputFile, "# \n");
  }
  
#ifdef SEARCH_STATISTICS_ON  
  globalStats.nrOfCleanUps          = 0;
  globalStats.nrOfLMLATables        = 0;
  globalStats.nrOfLMLATables_Temp   = 0;
  globalStats.nrOfActiveNodes       = 0;
  globalStats.nrOfActiveTokens      = 0;
  globalStats.nrOfLMLATables_Prev   = 0;
  globalStats.nrOfActiveNodes_Prev  = 0;
  globalStats.nrOfActiveTokens_Prev = 0;
  globalStats.nrOfLMLATables_Max    = 0;
  globalStats.nrOfActiveNodes_Max   = 0;
  globalStats.nrOfActiveTokens_Max  = 0;
#endif  
  fflush(outputFile);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method prints some final statistics such as average RealTime factor.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::printFinalSettings(bool outXML, int totMilliSec, int totTime)
{
  fprintf(outputFile, "\n");
  
  if(outXML)
  {
    fprintf(outputFile, "  </segments>\n");
    fprintf(outputFile, "  <statistics>\n");
  }
  else
  {  
    fprintf(outputFile, "##################################################################################################\n");
  }
   
#ifdef SEARCH_STATISTICS_ON  
  if(outXML)
  {
    float nr = ((float)(globalStats.nrOfCleanUps));
    fprintf(outputFile, "    <search_statistics>\n");    
    fprintf(outputFile, "      <LMLA_TABLES AVG=\"%.2f\" MAX=\"%d\"/>\n",
                  globalStats.nrOfLMLATables/nr,globalStats.nrOfLMLATables_Max);
    fprintf(outputFile, "      <ACTIVE_NODES AVG=\"%.2f\" MAX=\"%d\"/>\n",
                  globalStats.nrOfActiveNodes/nr,globalStats.nrOfActiveNodes_Max);
    fprintf(outputFile, "      <ACTIVE_TOKENS AVG=\"%.2f\" MAX=\"%d\"/>\n",
                  globalStats.nrOfActiveTokens/nr,globalStats.nrOfActiveTokens_Max);
    fprintf(outputFile, "    </search_statistics>\n");    
  }
  else
  {
    fprintf(outputFile, "### \n");
    fprintf(outputFile, "### ---------------------------------------------------\n");
    fprintf(outputFile, "###  Session statistics (averages per %d frames): \n", LMLA_INTERVAL);
    fprintf(outputFile, "### ---------------------------------------------------\n");
    float nr = ((float)(globalStats.nrOfCleanUps));
    fprintf(outputFile, "###  Number of LMLA tables   (avg,max) : %10.2f %6d\n",
                  globalStats.nrOfLMLATables/nr,globalStats.nrOfLMLATables_Max);
    fprintf(outputFile, "###  Number of active nodes  (avg,max) : %10.2f %6d\n",
                  globalStats.nrOfActiveNodes/nr,globalStats.nrOfActiveNodes_Max);
    fprintf(outputFile, "###  Number of active tokens (avg,max) : %10.2f %6d\n", 
                  globalStats.nrOfActiveTokens/nr,globalStats.nrOfActiveTokens_Max);
    fprintf(outputFile, "##################################################################################################\n");
  }
#endif
  if(analysisSettings != NULL)
  {
    int totSub = 0;
    int totDel = 0;
    int totIns = 0;
    for(int i=0;i<NUMBER_OF_BLAME_CLUSTERS;i++)
    {
      totSub += analyse_substitutions[i];
      totIns += analyse_insertions[i];
      totDel += analyse_deletions[i];
    }
    
    if(outXML)
    {
      fprintf(outputFile, "    <error_statistics>\n");
      fprintf(outputFile, "      <NR_WORDS_IN_REF>%d</NR_WORDS_IN_REF>\n",analyse_totRefWords);
      fprintf(outputFile, "      <word_error_rate INS=\"%d\" SUB=\"%d\" DEL=\"%d\" WER=\"%f\"/>\n",
                    totIns,totSub,totDel,100.0*(double)(totSub+totIns+totDel)/(double)(analyse_totRefWords));
      fprintf(outputFile, "      <blame_assignment>\n");
    }
    else
    {    
      fprintf(outputFile, "### \n");
      fprintf(outputFile, "### Number of UNKNOWN error regions: %d\n",blameAssingment[UNKNOWN]);
      fprintf(outputFile, "### Number of OOV     error regions: %d\n",blameAssingment[OOV]);
      fprintf(outputFile, "### Number of SEARCH  error regions: %d\n",blameAssingment[SEARCH]);
      fprintf(outputFile, "### Number of AM/LM   error regions: %d\n",blameAssingment[AMLM]);
      fprintf(outputFile, "### Number of LM      error regions: %d\n",blameAssingment[LM]);
      fprintf(outputFile, "### Number of AM      error regions: %d\n",blameAssingment[AM]);
      fprintf(outputFile, "### \n"); 
  
      fprintf(outputFile, "### Total: Number of words in REF : %d\n",analyse_totRefWords);
      fprintf(outputFile, "### Total: Number of insertions   : %d\n",totIns);
      fprintf(outputFile, "### Total: Number of substitutions: %d\n",totSub);
      fprintf(outputFile, "### Total: Number of deletions    : %d\n",totDel);
      fprintf(outputFile, "### Total: Estimated WER: %f%%\n",100.0*(double)(totSub+totIns+totDel)/(double)(analyse_totRefWords));        
      fprintf(outputFile, "### \n"); 
    }
    char prefix[20];
    for(int i=0;i<NUMBER_OF_BLAME_CLUSTERS;i++)
    {
      switch(i)
      {
        case 0:
          if(outXML)
          {
            fprintf(outputFile, "      <blame_class NAME=\"Unknown\" ");          
          }
          else
          {    
            strcpy(prefix, "### Unknown:");
          }
          break;
        case 1: 
          if(outXML)
          {
            fprintf(outputFile, "      <blame_class NAME=\"OOV\" ");
          }
          else
          {    
            strcpy(prefix, "### OOV:");
          }
          break;
        case 2: 
          if(outXML)
          {
            fprintf(outputFile, "      <blame_class NAME=\"Search\" ");
          }
          else
          {    
            strcpy(prefix, "### Search:");
          }
          break;
        case 3: 
          if(outXML)
          {
            fprintf(outputFile, "      <blame_class NAME=\"AM/LM\" ");
          }
          else
          {    
            strcpy(prefix, "### AM/LM:");
          }
          break;
        case 4: 
          if(outXML)
          {
            fprintf(outputFile, "      <blame_class NAME=\"LM\" ");
          }
          else
          {    
            strcpy(prefix, "### LM:");
          }
          break;
        case 5: 
          if(outXML)
          {
            fprintf(outputFile, "      <blame_class NAME=\"AM\" ");
          }
          else
          {    
            strcpy(prefix, "### AM:");
          }
          break;
      }
      if(outXML)
      {
        fprintf(outputFile, ">\n      <OCCURENCES>%d</OCCURENCES>\n",blameAssingment[i]);
        fprintf(outputFile, "      <word_error_rate INS=\"%d\" SUB=\"%d\" DEL=\"%d\" WER=\"%.2f\"/>\n",
                    analyse_insertions[i],
                    analyse_substitutions[i],
                    analyse_deletions[i],
                  100.0*(double)(analyse_substitutions[i]+analyse_insertions[i]+analyse_deletions[i]) /
                  (double)(analyse_totRefWords));        
        fprintf(outputFile, "      </blame_class>\n");
      }
      else
      {    
        fprintf(outputFile, "%s Number of insertions   : %d\n", prefix, analyse_insertions[i]);
        fprintf(outputFile, "%s Number of substitutions: %d\n", prefix, analyse_substitutions[i]);
        fprintf(outputFile, "%s Number of deletions    : %d\n", prefix, analyse_deletions[i]);
        fprintf(outputFile, "%s Estimated WER: %f%%\n", prefix,
                  100.0*(double)(analyse_substitutions[i]+analyse_insertions[i]+analyse_deletions[i]) /
                  (double)(analyse_totRefWords));        
        fprintf(outputFile, "### \n");
      }
    } 
    if(outXML)
    {
      fprintf(outputFile, "      </blame_assignment>\n");
      fprintf(outputFile, "    </error_statistics>\n");    
    }
    else
    {    
      fprintf(outputFile, "##################################################################################################\n"); 
    }
  }
  if(outXML)
  {
    fprintf(outputFile, "    <real_time milliseconds=\"%d\" ",totMilliSec);
    fprintf(outputFile, "frames=\"%d\" ",totTime);
    fprintf(outputFile, "RTF=\"%.4f\"/>\n",((double)totMilliSec/10.0)/((double)totTime));    
    fprintf(outputFile, "  </statistics>\n");
    fprintf(outputFile, "</shout_metadata>\n");  
  }
  else
  {    
    fprintf(outputFile, "#\n"); 
    fprintf(outputFile, "###################################################\n");
    fprintf(outputFile, "## Total processing time: %d milliseconds / %d frames = %2.4f times realtime.\n",
                    totMilliSec,totTime,((double)totMilliSec/10.0)/((double)totTime));      
    fprintf(outputFile, "###################################################\n");    
  }  
  fflush(outputFile);
}

      
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method determines the error regions and performs blame assignment.
/// The depth should be set to zero. The method uses this variable for recursive calls.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::errorAnalysis(WLRType *wlr, int depth, bool outputXML)  
{
  if(analysisSettings != NULL)
  {
    if(depth == 0)    // First time the function is called: initialise..
    {
      analysisSettings->errorRegionActive = 0;
      analysisSettings->containsOOV       = false;
      analysisSettings->prunePathAnalyses[0]->errorRegionID = -1;
      analysisSettings->prunePathAnalyses[0]->errorCategory = UNKNOWN;
      for(int i=1;i<analysisSettings->prunePathAnalysesLength;i++)
      {
        analysisSettings->prunePathAnalyses[i]->errorRegionID = 0;
        analysisSettings->prunePathAnalyses[i]->errorCategory = UNKNOWN;
      }
    }
    if((wlr != NULL) && (wlr->timeStamp >= 0))
    {
      // First do earlier words:
      errorAnalysis(wlr->previous, 1, outputXML);
      int wordID = getWordFromWLR(wlr);
      int match = 0;
      int i = analysisSettings->errorRegionActive+1;
      while(match == 0 && i < analysisSettings->prunePathAnalysesLength)
      {
        WLRType *w = &analysisSettings->prunePathAnalyses[i]->w;
        if(w->timeStamp == wlr->timeStamp && analysisSettings->prunePathAnalyses[i]->wordID == wordID)
        {
          if(compareLMHistory(w->lmHistory, wlr->lmHistory))
          {
            match = 1;
            analysisSettings->prunePathAnalyses[i]->linkWord = wlr;
          }
          w = w->previous;
          if(compareLMHistory(w->lmHistory, wlr->previous->lmHistory))
          {
            match = 2;
            analysisSettings->prunePathAnalyses[i]->linkWord = wlr;
          }
        }
        i++;
      }
      if(match == 2)  // CORRECT!
      {
        analysisSettings->errorRegionActive = i-1;
        analysisSettings->prunePathAnalyses[i-1]->errorRegionID = -1;
      }
      else if(match == 1)  // Next will be in new region (no interference of this one!)
      {
        analysisSettings->errorRegionActive = i-1;
        analysisSettings->prunePathAnalyses[i-1]->errorRegionID = -2;
      }

      if(depth == 0)
      {
        if(outputXML)
        {
          fprintf(outputFile, "      <error_regions>\n");
        }

        int  id          = 0;
        int  firstWrong  = -1;
        bool oov         = false;
        int  nrWordRef   = 0;
        int totalDepth   = 0;
        int nrTotWrong   = 0;
        analysisSettings->prunePathAnalyses[analysisSettings->prunePathAnalysesLength-1]->linkWord = wlr;
        for(int i=1;i<analysisSettings->prunePathAnalysesLength;i++)
        {
          int wordID = analysisSettings->prunePathAnalyses[i]->wordID;
          if(wordID != startOfSentenceWord)
          {
            totalDepth++;
          } 
          if(analysisSettings->prunePathAnalyses[i]->errorRegionID >= 0)
          {
            if(wordID != startOfSentenceWord)
            {
              nrWordRef++;
            }
            analysisSettings->prunePathAnalyses[i]->errorRegionID = id;
            if(firstWrong < 0)
            {
              firstWrong = i;
              oov = false;
            }
            if(analysisSettings->prunePathAnalyses[i]->wordID < 0)
            {
              oov = true;
            } 
          }
          if((analysisSettings->prunePathAnalyses[i]->errorRegionID >= 0 && i == analysisSettings->prunePathAnalysesLength-1) ||
              analysisSettings->prunePathAnalyses[i]->errorRegionID == -2) // last in region
          {
            if(wordID != startOfSentenceWord)
            {
              nrWordRef++;
            }
            if(firstWrong < 0)
            {
              firstWrong = i;
              oov = false;
            }
            if(analysisSettings->prunePathAnalyses[i]->wordID < 0)
            {
              oov = true;
            } 
            // Now do blame assignment:
            BlameCluster errorCategory = UNKNOWN;
            double totScore[2]= {0.0, 0.0};
            double amScore[2] = {0.0, 0.0};
            float  lmScore[2] = {0.0, 0.0};
            WLRType *lastWrong   = analysisSettings->prunePathAnalyses[i]->linkWord;
            WLRType *lastCorrect = NULL;
            totScore[0] = analysisSettings->prunePathAnalyses[i]->w.COMBlikelihood;
            lmScore[0]  = analysisSettings->prunePathAnalyses[i]->w.LMlikelihood;
            totScore[1] = lastWrong->COMBlikelihood;
            if(firstWrong > 1)
            {
              lastCorrect = analysisSettings->prunePathAnalyses[firstWrong-1]->linkWord;
              totScore[0]-= analysisSettings->prunePathAnalyses[firstWrong-1]->w.COMBlikelihood;
              lmScore[0] -= analysisSettings->prunePathAnalyses[firstWrong-1]->w.LMlikelihood;
              totScore[1]-= lastCorrect->COMBlikelihood;              
            }
            else
            {
              firstWrong = 1;
            }
            int nrWordHyp     = 0;
            int deletions     = 0;
            int insertions    = 0;
            int substitutions = 0;
            int nrCorrect     = 0;
            calcErrorRegionStats(lastWrong,lastCorrect, firstWrong, i, &nrWordHyp, &lmScore[1]);            
            for(int i2=firstWrong;i2<=i;i2++)
            {
              if(analysisSettings->prunePathAnalyses[i2]->linkWord != NULL)
              {
                nrCorrect++;
              }
            }
            // Test if the last word in the string is wrong! If so, decrease nrCorrect, because
            // it will have a linkWord...
            if(analysisSettings->prunePathAnalyses[i]->errorRegionID != -2) 
            {
              nrCorrect--;
            }
            
            analysisSettings->prunePathAnalyses[i]->errorRegionID = id;
            
            if(nrWordHyp > nrWordRef)
            {
              insertions = nrWordHyp - nrWordRef;
            }
            else if(nrWordHyp < nrWordRef)
            {
              deletions = nrWordRef - nrWordHyp;
            }
            substitutions = nrWordRef - nrCorrect - deletions;
            

            amScore[0] = totScore[0] - lmScore[0];
            lmScore[0] *= settings.weights_LmScale;
            totScore[0] = amScore[0] + lmScore[0];
            amScore[1] = totScore[1] - lmScore[1];
            
            // 0 = ref, 1 = hyp:
            if(oov)
            {
              errorCategory = OOV;
              analysisSettings->containsOOV = true;
            }
            else if(totScore[0] > totScore[1])
            {
              errorCategory = SEARCH;
            }
            else if(amScore[1] >= amScore[0] && lmScore[1] >= lmScore[0])
            {
              errorCategory = AMLM;
            }
            else if(amScore[1] < amScore[0] && lmScore[1] >= lmScore[0])
            {
              errorCategory = LM;
            }            
            else if(amScore[1] >= amScore[0] && lmScore[1] < lmScore[0])
            {
              errorCategory = AM;
            }            

            blameAssingment[errorCategory]++;
            analyse_deletions[errorCategory]     += deletions;
            analyse_insertions[errorCategory]    += insertions;
            analyse_substitutions[errorCategory] += substitutions;
            nrTotWrong                           += deletions + substitutions;
            
            for(int a=firstWrong;a<=i;a++)
            {
              analysisSettings->prunePathAnalyses[a]->errorCategory = errorCategory;
            }

            int amsim     = 0;
            if(false) //doPhoneAlignment)
            {
              // Okay, finally check for AM simularity:
              WLRType *wRef = &(analysisSettings->prunePathAnalyses[i]->w);
              WLRType *wHyp =   analysisSettings->prunePathAnalyses[i]->linkWord;
              if(wRef != NULL && wHyp != NULL)
              {
                int nrWords      = i - firstWrong + 1; 
                int correctPhone = 0;
                int nrPhone      = 0;
                int timeRef      = 0;
                int timeHyp      = 0;
                PLRType *pHyp = wHyp->phoneAlignment;
                timeHyp       = wHyp->timeStamp;
                if(pHyp->phoneID < 0)
                {
                  pHyp        = pHyp->previous;
                }

                while(nrWords > 0 && wRef != NULL && pHyp != NULL)
                {
                  PLRType *pRef = wRef->phoneAlignment;
                  timeRef       = wRef->timeStamp;
                  if(pRef->phoneID < 0)
                  {
                    pRef        = pRef->previous;
                  }
                  while(pRef != NULL && pHyp != NULL)
                  {
                    if((pRef->timeStamp-1) <= pHyp->timeStamp && (pRef->timeStamp+1) >= pHyp->timeStamp && 
                      (timeRef-1) <= timeHyp && (timeRef+1) >= timeHyp &&
                        pRef->phoneID == pHyp->phoneID)
                    {
                      correctPhone++;
                    }          
                    timeRef = pRef->timeStamp;
                    timeHyp = pHyp->timeStamp;
                    nrPhone++;      
                    pRef = pRef->previous;
                    pHyp = pHyp->previous;
                    if(pHyp == NULL)
                    {
                      wHyp = wHyp->previous;
                      if(wHyp != NULL && wHyp->timeStamp <= 0)
                      {
                        wHyp = NULL;
                      }
                      if(wHyp != NULL)
                      {
                        pHyp    = wHyp->phoneAlignment;
                        timeHyp = wHyp->timeStamp;
                        if(pHyp->phoneID < 0)
                        {
                          pHyp  = pHyp->previous;
                        }
                      }
                    }
                  }
                  nrWords--;
                  wRef = wRef->previous;
                  if(wRef->timeStamp <= 0)
                  {
                    wRef = NULL;
                  }
                }       
                amsim = 100*correctPhone / nrPhone;
              }
            }

            if(outputXML)
            {
              fprintf(outputFile, "        <error_region ID=\"%d\" cat=\"(",id);
              printErrorString(errorCategory);
              fprintf(outputFile, "\" am-sim=\"%d\" worderrors=\"%d\"/>\n",amsim, deletions + substitutions + insertions);
            }
            else
            {            
              fprintf(outputFile, "## Error region %02d: (",id);
              printErrorString(errorCategory);
              fprintf(outputFile, " AM simularity: %3d%% (%d errors)\n", amsim, deletions + substitutions + insertions);
            }
            
            
            id++;
            firstWrong = -1;
            oov        = false;
            nrWordRef  = 0;
          }
        }

        if(outputXML)
        {
          fprintf(outputFile, "      </error_regions>\n");
          fprintf(outputFile, "      <WORDS_CORRECT_NON_SIL>%d</WORDS_CORRECT_NON_SIL>\n",analysisSettings->refStats.noSilWords - nrTotWrong);
          fprintf(outputFile, "      <REF-ALIGNMENT>\n");
          fprintf(outputFile, "        <wordsequence>\n");
        }
        else
        {                
          fprintf(outputFile, "## Number of words correct (non-SIL)               : %d\n", analysisSettings->refStats.noSilWords - nrTotWrong);
          fprintf(outputFile, "## -----------------------------------------------------------\n");
        }          

        for(int i=1;i<analysisSettings->prunePathAnalysesLength;i++)
        {          
          WLRType *wlr = &analysisSettings->prunePathAnalyses[i]->w;
          int wordID = analysisSettings->prunePathAnalyses[i]->wordID;
          
          
          if(outputXML)
          {
            fprintf(outputFile, "      <word confidence=\"1\" wordID=");
            if(wordID >= 0)
            {
              char tmp[80];
              strcpy(tmp,vocabulary[wordID]);
              for(int i=0;i<(int)strlen(tmp);i++)
              {
                if(tmp[i] == '<')
                {
                  tmp[i] = '[';
                }
                if(tmp[i] == '>')
                {
                  tmp[i] = ']';
                }
              }            
              fprintf(outputFile, "\"%s\" ",tmp);
            }
            else
            {
              fprintf(outputFile, "\"[OOV]\" ");
            }
            
            fprintf(outputFile, "beginTime=\"%.3f\" ",((double)wlr->previous->timeStamp+1.0)/100.0);
            fprintf(outputFile, "endTime=\"%.3f\">\n",((double)wlr->timeStamp+1.0)/100.0);
            
            fprintf(outputFile, "          <error_region ID=\"%d\" cat=\"(",analysisSettings->prunePathAnalyses[i]->errorRegionID);
            printErrorString(analysisSettings->prunePathAnalyses[i]->errorCategory);
            fprintf(outputFile, "\"/>\n");
            
            fprintf(outputFile, "        <score AM=\"%.5f\" ",wlr->COMBlikelihood - wlr->LMlikelihood);
            fprintf(outputFile, "LM=\"%.5f\" ",wlr->LMlikelihood);
            fprintf(outputFile, "COMBINED=\"%.5f\"/>\n",wlr->COMBlikelihood);
            if(doPhoneAlignment)
            {
              fprintf(outputFile, "        <phonesequence>\n");
              getPhoneAlignment(NULL, wlr->phoneAlignment, wlr->timeStamp,true);          // NULL = XML...
              fprintf(outputFile, "        </phonesequence>\n");
            }
            fprintf(outputFile, "      </word>\n");
          }
          else
          {          
            if(doPhoneAlignment)
            {
              getPhoneAlignment((char*)"#ANALYSE#", wlr->phoneAlignment, wlr->timeStamp,true);
            }

            char tmpStr[80];
            fprintf(outputFile, "#ANALYSE# ");
            int len = 40;

            if(analysisSettings->prunePathAnalyses[i]->errorRegionID >= 0)
            {
              sprintf(tmpStr,"(error_%02d:",analysisSettings->prunePathAnalyses[i]->errorRegionID);
              len -= strlen(tmpStr);
              fprintf(outputFile, "%s",tmpStr);
              len -= printErrorString(analysisSettings->prunePathAnalyses[i]->errorCategory);
            }
                      
            if(wordID >= 0)
            {
              fprintf(outputFile, "%s",vocabulary[wordID]);
              len -= strlen(vocabulary[wordID]);  
            }
            else
            {
              fprintf(outputFile, "<OOV>");
              len -= 5;
            }                  
              
            while(len>0)
            {
              fprintf(outputFile, " ");
              len--;
            }
            fprintf(outputFile, "%6d      %12.5f    %12.5f    %12.5f\n",
                                    (int)wlr->timeStamp,
                                    wlr->COMBlikelihood - wlr->LMlikelihood,
                                    wlr->LMlikelihood,
                                    wlr->COMBlikelihood);
          }                             
        }
        if(outputXML)
        {
          fprintf(outputFile, "        </wordsequence>\n");
          fprintf(outputFile, "      </REF-ALIGNMENT>\n");
        }

        analyse_totRefWords += totalDepth;
      }    
    }
  }
}

int LexicalTree::printErrorString(int errorID)
{
  int len = 0;
  switch(errorID)
  {
    case OOV:
      fprintf(outputFile, "OOV) ");
      len = 5;
      break;
    case SEARCH:
      fprintf(outputFile, "SEARCH) ");
      len = 7;
      break;
    case AMLM:
      fprintf(outputFile, "AM/LM) ");
      len = 7;
      break;
    case LM:
      fprintf(outputFile, "LM) ");
      len = 4;
      break;
    case AM:
      fprintf(outputFile, "AM) ");
      len = 4;
      break;
    default:
      fprintf(outputFile, "UNKNOWN) ");
      len = 9;
      break;
  }
  return len;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Calculates the LM score of an error region.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::calcErrorRegionStats(WLRType *wlr, WLRType *lastCorrect, int firstRef, int lastRef, 
                                       int *nrWordHyp, float *scoreLM)
{  
  *scoreLM   = 0.0;
  *nrWordHyp = 0;
  if((wlr != NULL) && (wlr->timeStamp >= 0))
  {
    if(wlr->previous != lastCorrect)
    {
      calcErrorRegionStats(wlr->previous,lastCorrect, firstRef, lastRef, nrWordHyp, scoreLM);
    }

    int wordID = -1;
    if(wlr->isSil == 1)
    {
      wordID = startOfSentenceWord;
    }
    else
    {
      (*nrWordHyp)++;  // SIL is not counted!
      if(languageModel != NULL)
      {
        wordID = languageModel->getLastWordID(wlr->lmHistory);
        if(wordID >= 0)
        {
          int dummy[LM_NGRAM_DEPTH];
          (*scoreLM) += languageModel->getP(wordID,wlr->previous->lmHistory,dummy)*settings.weights_LmScale;
        }
      }
      else
      {
        wordID = wlr->lmHistory[0];
      }
      int i = firstRef;
      while(i <= lastRef && !(analysisSettings->prunePathAnalyses[i]->linkWord == NULL &&
                              analysisSettings->prunePathAnalyses[i]->wordID == wordID)    )
      {
        i++;
      }
      if(i <= lastRef && analysisSettings->prunePathAnalyses[i]->linkWord == NULL)
      {
        analysisSettings->prunePathAnalyses[i]->linkWord = wlr;
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method recursively goes through the Word Link Records of a (normaly the best) recognition and
/// prints to recognized words to the string 'string'.
///
/// If 'complete' is true, also alignment and likelihood data is print.
///
/// If phone alignment is switched on (configuration switch) also the phone alignment data
/// is printed to the string.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::getBestPath(WLRType *wlr, bool outputXML, bool complete)
{  
  if((wlr != NULL) && (wlr->timeStamp >= startTime))
  {
    getBestPath(wlr->previous,outputXML,complete);
    int wordID = -1;
    if(wlr->isSil == 1)
    {
      wordID = startOfSentenceWord;
    }
    else
    {
      if(languageModel != NULL)
      {
        wordID = languageModel->getLastWordID(wlr->lmHistory);
      }
      else
      {
        wordID = wlr->lmHistory[0];
      }
    } 
    
    double confident = 1.0;   
    if(phoneLoopConfidence != NULL)
    {
      int totTime = wlr->timeStamp - wlr->previous->timeStamp;
      if(totTime < 1)
      {
        totTime = 1;
      }
      double amLike   = (wlr->COMBlikelihood - wlr->LMlikelihood) - (wlr->previous->COMBlikelihood - wlr->previous->LMlikelihood);
      double confLike = phoneLoopConfidence[wlr->timeStamp - phoneLoopConfidenceOffset] - 
                        phoneLoopConfidence[wlr->previous->timeStamp - phoneLoopConfidenceOffset];
      confident = exp((amLike - confLike))/((double)totTime);
    }
    
    if(complete)
    {  
      if(outputXML)
      {
        if(phoneLoopConfidence != NULL)
        {
          fprintf(outputFile, "      <word confidence=\"%6g\" wordID=",confident);
        }
        else
        {
          fprintf(outputFile, "      <word wordID=");
        }
        if(wordID >= 0)
        {
          char tmp[80];
          strcpy(tmp,vocabulary[wordID]);
          for(int i=0;i<(int)strlen(tmp);i++)
          {
            if(tmp[i] == '<')
            {
              tmp[i] = '[';
            }
            if(tmp[i] == '>')
            {
              tmp[i] = ']';
            }
          }        
          fprintf(outputFile, "\"%s\" ",tmp);
        }
        else
        {
          fprintf(outputFile, "\"[OOV]\" ");
        }
        fprintf(outputFile, "beginTime=\"%.3f\" ",((double)wlr->previous->timeStamp+1.0)/100.0);
        fprintf(outputFile, "endTime=\"%.3f\">\n",((double)wlr->timeStamp+1.0)/100.0);
        fprintf(outputFile, "        <score AM=\"%.5f\" ",wlr->COMBlikelihood - wlr->LMlikelihood);
        fprintf(outputFile, "LM=\"%.5f\" ",wlr->LMlikelihood);
        fprintf(outputFile, "COMBINED=\"%.5f\"/>\n",wlr->COMBlikelihood);
        if(doPhoneAlignment)
        {
          fprintf(outputFile, "        <phonesequence>\n");
          getPhoneAlignment(NULL, wlr->phoneAlignment, wlr->timeStamp,confident > PHONE_CONFIDENCE_THRESHOLD); // <- NULL means XML!
          fprintf(outputFile, "        </phonesequence>\n");
        }
        fprintf(outputFile, "      </word>\n");
      }
      else
      { 
        if(doPhoneAlignment)
        {
          getPhoneAlignment((char*)"#", wlr->phoneAlignment, wlr->timeStamp,confident > PHONE_CONFIDENCE_THRESHOLD);
        }
        fprintf(outputFile, "# ");
        int len = 35;
        if(wordID >= 0)
        {
          fprintf(outputFile, "%s",vocabulary[wordID]);
          len = 40 - strlen(vocabulary[wordID]);  
        }
        else
        {
          fprintf(outputFile, "<OOV>");
        }

        while(len>0)
        {
          fprintf(outputFile, " ");
          len--;
        }
        fprintf(outputFile, "%6d      %12.5f    %12.5f    %12.5f\n",
                                (int)wlr->timeStamp,
                                wlr->COMBlikelihood - wlr->LMlikelihood,
                                wlr->LMlikelihood,
                                wlr->COMBlikelihood);
      }
    }
    else
    {
      if(!confident)
      {
        fprintf(outputFile, "*");
      }
      if(wordID >= 0)
      {
        fprintf(outputFile, "%s", vocabulary[wordID]);
      }
      else
      {
        fprintf(outputFile, "<OOV>");
      }
      fprintf(outputFile, " ");
    }
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will safe the word-path of the best recognition. This golden recognition is used 
/// for statistical analysis (blame assignment).
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool LexicalTree::safeBestRecognition(bool addEndOfSentence)
{
  // First delete old safings:
  if(analysisSettings != NULL)
  {
    if(analysisSettings->prunePathAnalyses != NULL)
    {
      for(int i=0;i<analysisSettings->prunePathAnalysesLength;i++)
      {
        if(doPhoneAlignment)
        {
          PhoneModel::initialisePhonePath(analysisSettings->prunePathAnalyses[i]->w.phoneAlignment);
        }
        delete analysisSettings->prunePathAnalyses[i];
        COUNTER_MINUS(trackCount);
      }
      delete[] analysisSettings->prunePathAnalyses;
    }
    if(analysisSettings->binDistr != NULL)
    {
      delete[] (analysisSettings->binDistr);
      delete[] (analysisSettings->bestScore);
      delete[] (analysisSettings->correctScore);
      delete[] (analysisSettings->bestStateScore);
      delete[] (analysisSettings->stateRanking);
    }    
    delete analysisSettings;
    analysisSettings = NULL;
  }

  // Create a new path:

  analysisSettings = new AnalysisSettings;
  analysisSettings->prunePathAnalyses        = NULL;
  analysisSettings->binDistr                 = NULL;
  analysisSettings->bestScore                = NULL;
  analysisSettings->correctScore             = NULL;
  analysisSettings->bestStateScore           = NULL;
  analysisSettings->stateRanking             = NULL;
  analysisSettings->endStateActive           = NULL;
  analysisSettings->refStats.noSilWords      = 0;
  analysisSettings->refStats.nrWords         = 0;
  analysisSettings->refStats.totAM           = 0.0;
  analysisSettings->refStats.totLM           = 0.0;
  for(int i=0;i<3;i++)
  {
    analysisSettings->beamCategories[i]      = 0;
  }
  analysisSettings->prunePathAnalysesLength  = 0;
  analysisSettings->containsOOV              = false;
  analysisSettings->errorRegionActive        = 0;
  analysisSettings->maxBeam                  = -10.0;
  analysisSettings->maxStateBeam             = -10.0;
  analysisSettings->maxEndStateBeam          = -10.0;
  analysisSettings->maxHistogram             = -10;
  analysisSettings->maxStateHistogram        = -10;
  for(int i=0;i<LM_NGRAM_DEPTH;i++)
  {
    analysisSettings->currentLMHist[i]       = startOfSentenceWord;
  }

  TokenType *useToken = findBestToken(addEndOfSentence, true, false);
  if(useToken != NULL && useToken->path!=NULL)
  { 
    WLRType *w  = useToken->path;
    analysisSettings->prunePathAnalysesLength = 0;
    while(w != NULL)
    {
      if(w->timeStamp >= 0)
      {
        analysisSettings->prunePathAnalysesLength++;
      }
      w = w->previous;
    }

    analysisSettings->refStats.nrWords = analysisSettings->prunePathAnalysesLength-1;
    if(analysisSettings->prunePathAnalysesLength > 0)  //It always is!
    {
      w  = useToken->path;
      if(w->timeStamp < 0)
      {
        w = w->previous;
      }
      analysisSettings->refStats.totAM = (w->COMBlikelihood - w->LMlikelihood);

      analysisSettings->prunePathAnalyses = new WLRTracker*[analysisSettings->prunePathAnalysesLength];
      for(int i=0;i<analysisSettings->prunePathAnalysesLength;i++)
      {
        analysisSettings->prunePathAnalyses[i] = NULL;
      }

      int nr = analysisSettings->prunePathAnalysesLength-1;
      analysisSettings->prunePathAnalyses[nr] = new WLRTracker;
      COUNTER_PLUS(trackCount);

      if(w->isSil == 1)
      {
        analysisSettings->prunePathAnalyses[nr]->wordID = startOfSentenceWord;
      }
      else if(languageModel != NULL)
      {
        int dummy[LM_NGRAM_DEPTH];
        int wordID = languageModel->getLastWordID(w->lmHistory);
        analysisSettings->refStats.noSilWords++;
        if(w->previous != NULL)
        {
          analysisSettings->refStats.totLM += languageModel->getP(wordID,w->previous->lmHistory,dummy);
        }
        analysisSettings->prunePathAnalyses[nr]->wordID = wordID;
      }
      copyLMHistory(analysisSettings->prunePathAnalyses[nr]->w.lmHistory, w->lmHistory);
      analysisSettings->prunePathAnalyses[nr]->w.isSil         = w->isSil;
      analysisSettings->prunePathAnalyses[nr]->w.timeStamp     = w->timeStamp;
      analysisSettings->prunePathAnalyses[nr]->w.COMBlikelihood= w->COMBlikelihood;
      analysisSettings->prunePathAnalyses[nr]->w.LMlikelihood  = w->LMlikelihood;
      analysisSettings->prunePathAnalyses[nr]->contextNext     = -1;
      analysisSettings->prunePathAnalyses[nr]->linkWord        = NULL;
      if(doPhoneAlignment)
      {
        analysisSettings->prunePathAnalyses[nr]->w.phoneAlignment = PhoneModel::copyPhonePath(w->phoneAlignment);
      }
      w = w->previous;
      nr--;
      while(nr >= 0)
      {
        analysisSettings->prunePathAnalyses[nr] = new WLRTracker;
        COUNTER_PLUS(trackCount);
        if(w->isSil == 1)
        {
          analysisSettings->prunePathAnalyses[nr]->wordID = startOfSentenceWord;
        }
        else if(languageModel != NULL)
        {
          int dummy[LM_NGRAM_DEPTH];
          int wordID = languageModel->getLastWordID(w->lmHistory);
          analysisSettings->refStats.noSilWords++;
          if(w->previous != NULL)
          {
            analysisSettings->refStats.totLM += languageModel->getP(wordID,w->previous->lmHistory,dummy);
          }
          analysisSettings->prunePathAnalyses[nr]->wordID = wordID;
        }
        copyLMHistory(analysisSettings->prunePathAnalyses[nr]->w.lmHistory, w->lmHistory);
        analysisSettings->prunePathAnalyses[nr]->w.isSil         = w->isSil;
        analysisSettings->prunePathAnalyses[nr]->w.timeStamp     = w->timeStamp;
        analysisSettings->prunePathAnalyses[nr]->w.COMBlikelihood= w->COMBlikelihood;  
        analysisSettings->prunePathAnalyses[nr]->w.LMlikelihood  = w->LMlikelihood;
        analysisSettings->prunePathAnalyses[nr]->contextNext     = -1;
        analysisSettings->prunePathAnalyses[nr]->linkWord        = NULL;
        if(doPhoneAlignment)
        {
          if(w->phoneAlignment != NULL)
          {
            int prev = (int)(w->phoneAlignment->contextKey/numberOfPhones);
            analysisSettings->prunePathAnalyses[nr]->contextNext = (int)(w->phoneAlignment->contextKey - prev*numberOfPhones);
          }
          analysisSettings->prunePathAnalyses[nr]->w.phoneAlignment = PhoneModel::copyPhonePath(w->phoneAlignment);
        }
        nr--;
        w = w->previous;
      }
      analysisSettings->prunePathAnalyses[0]->w.previous = NULL;
      for(int i=1;i<analysisSettings->prunePathAnalysesLength;i++)
      {
        analysisSettings->prunePathAnalyses[i]->w.previous = &(analysisSettings->prunePathAnalyses[i-1]->w);
      }
    }
    // Okay, now prepare our PRUNING analysis:
    int totLength = analysisSettings->prunePathAnalyses[analysisSettings->prunePathAnalysesLength-1]->w.timeStamp + 1;
    analysisSettings->binDistr           = new int[NR_OF_HISTOGRAM_BINS*totLength];
    analysisSettings->bestScore          = new float[totLength];
    analysisSettings->correctScore       = new float[totLength];
    analysisSettings->bestStateScore     = new float[totLength];
    analysisSettings->stateRanking       = new int[totLength];
    analysisSettings->maxBeam            = -10.0;
    analysisSettings->maxStateBeam       = -10.0;
    analysisSettings->maxEndStateBeam    = -10.0;
    analysisSettings->maxHistogram       = -10;
    analysisSettings->maxStateHistogram  = -10;
  }
  return (analysisSettings != NULL);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// The best grammar end-token is returned.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::getBestGrammarEndToken(TokenType **bestT, bool notFinishedAllowed)
{
  float res = SMALLEST_NEGATIVE;
  *bestT    = NULL;
  if((*endNode->inputToken) != NULL)
  {
    TokenType *t = (*endNode->inputToken);
    *bestT   = (*endNode->inputToken);
    res     = (*endNode->inputToken)->likelihood;
    while(t != NULL)
    {
      if(t->likelihood > res)
      {
        *bestT   = t;
        res      = t->likelihood;
      }
      t = t->next;
    }
  }
  if(*bestT == NULL && notFinishedAllowed && bestTokenInSystem != NULL)
  {
    *bestT = bestTokenInSystem;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will look in the lexical tree for the best token. If no token is found, it will
/// return NULL. If a pointer to a string is provided, logging will be stored in this string.
/////////////////////////////////////////////////////////////////////////////////////////////////////

TokenType *LexicalTree::findBestToken(bool addEndOfSentence, bool complete, bool outputXML, bool notFinishedAllowed)
{
  TokenType *useToken = NULL;
  float useAddScore   = 0.0;
  if(endNode != NULL)
  {
    getBestGrammarEndToken(&useToken,notFinishedAllowed);
  }
  else if(currentlyAligning)
  {
    useToken = *(grammarStart->inputToken);
  }
  if(useToken == NULL && treeStart != NULL)
  {
    if((latticeWordList == NULL) && (!currentlyAligning || complete))
    {
      for(int a=0;a<numberOfPhones;a++)
      {
        for(int a2=0;a2<numberOfPhones;a2++)
        {
          int index = a+numberOfPhones*a2;
          if(treeStart[index] != NULL && (*treeStart[index]->inputToken) != NULL)
              // phoneModels[treeStart[index]->modelID]->isSilModel())  : if you want to force the sentence to end with <s>!
          {
            if(addEndOfSentence)  // We'll have to sort on endOfSentence...
            {
              TokenType *searchToken = (*treeStart[index]->inputToken);
              while(searchToken != NULL)
              {
                int dummy[LM_NGRAM_DEPTH];
                float lmLikelihood = languageModel->getP(endOfSentenceWord,searchToken->path->lmHistory,dummy)*settings.weights_LmScale;
                if((useToken == NULL || useToken->likelihood + useAddScore < (*treeStart[index]->inputToken)->likelihood+lmLikelihood))
                {
                  useToken    = searchToken;
                  useAddScore = lmLikelihood;
                }
                searchToken = searchToken->next;
              }

            }
            else
            {
              TokenType *t = *treeStart[index]->inputToken;
              while(t != NULL)
              {
                if(useToken == NULL || t->likelihood > useToken->likelihood)
                {
                  useToken = t;
                }
                t = t->next;
              }
            }
          }
        }
      }
    }
    else if(complete)
    {
      fprintf(outputFile, "\n## WARNING: THE DECODER HAS NOT REACHED THE FINAL STATE, WORDS WILL BE MISSING FROM THE ALIGNMENT: %s\n",
               getAlignmentString());
    }
  }
  return useToken;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Cleans up the lattice administration.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::deleteLatticeAdmin()
{
  LatticeNode *ln = latticeAdmin;
  while(ln != NULL)
  {
    WLRList *l = ln->inArcs;
    while(l != NULL)
    {
      ln->inArcs = l->next;
      delete l;
      l = ln->inArcs;
    }
    l = ln->outArcs;
    while(l != NULL)
    {
      ln->outArcs = l->next;
      delete l;
      l = ln->outArcs;
    }    
    latticeAdmin = ln->adminNext;
    delete ln;
    ln = latticeAdmin;
  }
}



/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This lattice will search for existing lattice-nodes that are needed for wlr. If either of them
/// does not exists, it will be created.
/////////////////////////////////////////////////////////////////////////////////////////////////////

LatticeNode *LexicalTree::findLatticeNodes(WLRType *wlr)
{
  int timeFrom = wlr->previous->timeStamp;
  int timeTo   = wlr->timeStamp;
  int wordID   = getWordFromWLR(wlr);

  LatticeNode *l = latticeAdmin->adminNext;
  float beginDiff = l->timeBegin - timeFrom;
  float endDiff   = l->timeEnd - timeTo;
  beginDiff = sqrt(beginDiff*beginDiff);
  endDiff   = sqrt(endDiff*endDiff);
  while((l != NULL) && !((l->wordID == wordID) && (beginDiff < 3.0) && (endDiff < 3.0) && (l->timeBegin < timeTo) && (l->timeEnd > timeFrom)))
  {
    l = l->adminNext;
    if(l != NULL)
    {
      beginDiff = l->timeBegin - timeFrom;
      endDiff   = l->timeEnd - timeTo;
      beginDiff = sqrt(beginDiff*beginDiff);
      endDiff   = sqrt(endDiff*endDiff);
    }
  }
  if(l == NULL) // Node does not exist yet. Let's make it:
  {
    l = new LatticeNode;
    l->exampleWord = wlr;
    l->inArcs    = NULL;
    l->outArcs   = NULL;
    l->timeBegin = timeFrom;
    l->timeEnd   = timeTo;
    l->wordID    = wordID;
    l->adminNext = latticeAdmin->adminNext;
    latticeAdmin->adminNext = l;
  }
  return l;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::lattice_copyNonSilArcs(LatticeNode *source, LatticeNode *dest, int loopID)
{
  dest->nodeNr   = loopID;
  source->nodeNr = loopID;
  if(source->wordID == startOfSentenceWord && source->exampleWord != NULL)
  {
    double amScore  = 0.0;
    double lmScore  = 0.0;
    double totScore = 0.0;
    // First mark the original arc source-dest for deletion:
    WLRList *wl     = dest->inArcs;
    while(wl != NULL)
    {
      if((wl->wlr->lattice == source) && !(wl->totScore==2910.0 && wl->amScore==-2910.0 && wl->lmScore == 0.0))
      {
        amScore = wl->amScore;
        lmScore = wl->lmScore;
        totScore= wl->totScore;
        wl->totScore = 2910.0;
        wl->amScore  = -2910.0;
        wl->lmScore  = 0.0;
      }
      wl = wl->next;
    }
    // Now copy the inArcs from source to dest:
    wl = source->inArcs;
    while(wl != NULL)
    {
      if(wl->wlr->lattice->nodeNr != loopID)
      {
        if(wl->wlr->lattice->wordID == startOfSentenceWord && wl->wlr->lattice->exampleWord != NULL)
        {
          lattice_copyNonSilArcs(wl->wlr->lattice,dest,loopID);
        }
        else
        {
          WLRList *copyW = new WLRList;
          copyW->amScore   = wl->amScore + amScore;
          copyW->lmScore   = wl->lmScore + lmScore;
          copyW->totScore  = wl->totScore + totScore;
          copyW->wlr       = wl->wlr;
          copyW->next      = dest->inArcs;
          dest->inArcs     = copyW;
        }
      }
      wl = wl->next;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will create the lattice structure, given a set with WLR objects (after decoding).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::createLattice()
{
  // First clean up the administration:
  // Find best token:
  latticeN = 0;
  latticeL = 0;
  TokenType *useToken = findBestToken(false, false, false);
  if(useToken == NULL)
  {
    return;
  }
  // Delete last-minute, bad paths:
  if(useToken->next != NULL)
  {
    TokenType *t = useToken->next;
    PhoneModel::initialiseToken(&t);
    useToken->next = NULL;
  }

  // Start cleaning:
  createActiveNodesList();
  timeStamp++;
  touchWLRs(useToken);
  processVector_administrationCleanup();
  timeStamp--;
  deleteLatticeAdmin();

  latticeN = 0;
  latticeL = 0;

  // Now create the extra WLR's that are in the nBest list:
  // So that there is a path from each Nbest node to the best.
  // Note that the LM history is not valid without creating at least N extra nodes!
  assert(useToken->path != NULL);

  // First node (time: zero, wordID: startOfSentence)
  latticeAdmin = new LatticeNode;
  latticeAdmin->adminNext = NULL;
  latticeAdmin->inArcs    = NULL;
  latticeAdmin->outArcs   = NULL;
  latticeAdmin->timeBegin =  startTime;
  latticeAdmin->timeEnd   =  startTime;
  latticeAdmin->exampleWord = NULL;
  latticeAdmin->wordID    = startOfSentenceWord;

  LatticeNode *lastNode   = new LatticeNode;
  latticeAdmin->adminNext = lastNode;
  lastNode->adminNext = NULL;
  lastNode->inArcs    = NULL;
  lastNode->outArcs   = NULL;
  lastNode->exampleWord = NULL;
  lastNode->timeBegin =  -1;
  lastNode->timeEnd   =  -1;
  lastNode->wordID    = -1;


  createLatticeNodeGroups(useToken->path,0.0,-1);
  WLRList *lastW = new WLRList;
  lastW->amScore   = 0.0;
  lastW->lmScore   = 0.0;
  lastW->totScore  = 0.0;
  lastW->wlr       = useToken->path;
  lastW->next      = lastNode->inArcs;
  lastNode->inArcs = lastW;

  if(useToken->path->nBest != NULL)
  {
    for(int i=0;i<NBEST_DEPTH;i++)
    {
      if((useToken->path->nBest[i] != NULL) && (useToken->path->nBest[i] != useToken->path))
      {
        createLatticeNodeGroups(useToken->path->nBest[i],0.0,-1);
        WLRList *lastW = new WLRList;
        lastW->amScore   = 0.0;
        lastW->lmScore   = 0.0;
        lastW->totScore  = 0.0;
        lastW->wlr       = useToken->path->nBest[i];
        lastW->next      = lastNode->inArcs;
        lastNode->inArcs = lastW;
      }
    }
  }

  LatticeNode *l = NULL;
  lattice_removeDoubleArcs();

//   // Remove all SIL words... First set all nodeIDs to zero..
//   l = latticeAdmin;
//   while(l != NULL)
//   {
//     l->nodeNr = 0;
//     l = l->adminNext;
//   }
//   // now copy inarcs from silence words.
//   l = latticeAdmin;
//   while(l != NULL)
//   {
//     int lID = 2910; // my birthday..
//     if(l->wordID != startOfSentenceWord)
//     {
//       WLRList *wl = l->inArcs;
//       while(wl != NULL)
//       {
//         lID++;
//         lattice_copyNonSilArcs(wl->wlr->lattice,l,lID);
//         wl = wl->next;
//       }
//       // Now delete inArcs that are not pointing to anything any more:
//       wl = l->inArcs;
//       WLRList *wlPrev = NULL;
//       while(wl != NULL)
//       {
//         if(wl->totScore==2910.0 && wl->amScore==-2910.0 && wl->lmScore == 0.0) // marked for deletion..
//         {
//           if(wlPrev == NULL)
//           {
//             l->inArcs = wl->next;
//           }
//           else
//           {
//             wlPrev->next = wl->next;
//           }
//           WLRList *delwl = wl;
//           wl = wl->next;
//           delete delwl;
//         }
//         else
//         {
//           wlPrev = wl;
//           wl = wl->next;
//         }
//       }
//     }
//     l = l->adminNext;
//   }
//   // Now remove the silence words:
//   LatticeNode *lPrev = NULL;
//   l = latticeAdmin;
//   while(l != NULL)
//   {
//     if(l->wordID == startOfSentenceWord && l->exampleWord != NULL)
//     {
//       LatticeNode *removeL = l;
//       l = l->adminNext;
//       if(l != NULL)
//       {
//         if(lPrev == NULL)
//         {
//           latticeAdmin = l;
//         }
//         else
//         {
//           lPrev->adminNext = l;
//         }
//         WLRList *wl = removeL->inArcs;
//         while(wl != NULL)
//         {
//           WLRList *removeWl = wl;
//           wl = wl->next;
//           delete removeWl;
//         }
//         delete removeL;
//       }
//       else
//       {
//         lPrev = removeL;  // doing nothing..
//       }
//     }
//     else
//     {
//       lPrev = l;
//       l = l->adminNext;
//     }
//   }
// 
// 
//   lattice_removeDoubleArcs();


  // numbering nodes and arcs:
  l = latticeAdmin;
  while(l != NULL)
  {
    l->nodeNr = latticeN;
    latticeN++;
    WLRList *wl = l->inArcs;
    int t = 0;
    if(wl != NULL)
    {
      assert(wl->wlr != NULL);
      t = wl->wlr->timeStamp;
    }
    while(wl != NULL)
    {
      latticeL++;
      wl = wl->next;
    }

    l = l->adminNext;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo
/////////////////////////////////////////////////////////////////////////////////////////////////////


void LexicalTree::lattice_removeDoubleArcs()
{
  // Remove all duplicate arcs (same from/to node) and
  // determine all outgoing arcs:
  LatticeNode *l = latticeAdmin;
  while(l != NULL)
  {
    WLRList *wl = l->inArcs;
    while(wl != NULL)
    {
      WLRType *w = wl->wlr;
      
      LatticeNode *l2 = w->lattice;
      assert(l2 != NULL);
      
      // Now first remove arcs in inArcs, after this arc that have the same from-node:
      WLRList *wlPrev = wl;
      WLRList *wl2    = wl->next;
      while(wl2 != NULL)
      {
        if(wl2->wlr->lattice == l2) // Same!
        {
          if(wl->totScore < wl2->totScore) // The original is not as good... Switch arcs before deleting:
          {
            wl->amScore  = wl2->amScore;
            wl->lmScore  = wl2->lmScore;
            wl->totScore = wl2->totScore;
            wl->wlr      = wl2->wlr;
          }
          wlPrev->next = wl2->next;
          delete wl2;
          wl2          = wlPrev;
        }
        wlPrev = wl2;
        wl2 = wl2->next;
      }
      wl = wl->next;
    }
    l = l->adminNext;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will create a lattice, rescore the LM-parameters in various ways and print the 
/// resulting statistics.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::printLMParStats(bool outputXML)
{
  int SCALESTEPS   = 21;
  int TRANSSTEPS   = 1;
  
  float SCALESTART = 0.0;
  float SCALESIZE  = 2.5;
  
  float TRANSSTART = 0.0;
  float TRANSSIZE  = 5.0;
  
  
  if(analysisSettings == NULL)
  {
    USER_ERROR("Sorry, it is not possible to print LM-parameter statistics when ANALYSE is not switched on!");
  }

  if(!analysisSettings->containsOOV)
  {
    if(latticeAdmin == NULL)
    {
      createLattice();
    }

    sortLatticePaths(latticeAdmin);
    int nrPaths = countLatticePaths(latticeAdmin);
    if(outputXML)
    {
      fprintf(outputFile, "      <model_scaling NBEST_LENGTH=\"%d\">\n",nrPaths);
    }
    else
    {
      fprintf(outputFile, "# Number of Paths in the N-Best list: %d\n",nrPaths);
    }
    if(nrPaths > MAX_ANALYSIS_NBEST_PATHS)
    {
      nrPaths = MAX_ANALYSIS_NBEST_PATHS;
    }
    NBestList scoreList[nrPaths];
    int hist[LM_NGRAM_DEPTH];
    copyLMHistory(hist,initialLMHist);
    int usedPaths =  createLatticeLMRescoring(latticeAdmin, scoreList, 0.0, 0.0, 0, 0, hist);
    if(!outputXML)
    {
      fprintf(outputFile, "# Number of Paths in the N-Best list actually used: %d\n",usedPaths);
    }


    if(!outputXML)
    {    
      // Now start rescoring and printing:
      fprintf(outputFile, "#PAR#        ");
      int   scaleStep    = 0;
      float scaleSize    = SCALESTART;    
      while(scaleStep < SCALESTEPS)
      {
        fprintf(outputFile, "%5g ",scaleSize);
        scaleSize += SCALESIZE;
        scaleStep++;
      }
      fprintf(outputFile, "\n");
    }
        
    int   transStep    = 0;
    float transPenalty = TRANSSTART;
    while(transStep < TRANSSTEPS)
    {
      if(outputXML)
      {
        fprintf(outputFile, "      <TR ");
      }
      else
      {
        fprintf(outputFile, "#PAR# %5g: ", transPenalty);
      }
      int   scaleStep    = 0;
      float scaleSize    = SCALESTART;
      while(scaleStep < SCALESTEPS)
      {
        int   nrWorse = 0;
        float refScore = analysisSettings->refStats.totAM + 
                         analysisSettings->refStats.totLM*scaleSize
                         - analysisSettings->refStats.noSilWords*transPenalty;
        for(int p=0;p<nrPaths;p++)
        {
          float score = scoreList[p].totAM + scoreList[p].totLM*scaleSize
                         - scoreList[p].noSilWords*transPenalty;
          if(score <= refScore)
          {
            nrWorse++;
          }
        }
        int res = nrWorse*100/nrPaths;
        if(outputXML)
        {
          fprintf(outputFile, "LM%02d=\"%d\" ",scaleStep,res);
        }
        else
        {
          fprintf(outputFile, "%5d ",res);
        }
        
        scaleSize += SCALESIZE;
        scaleStep++;
      }
      if(outputXML)
      {
        fprintf(outputFile, "/>\n");      
      }
      else
      {
        fprintf(outputFile, "\n");
      }
      transPenalty += TRANSSIZE;
      transStep++;
    }
    if(outputXML)
    {
      fprintf(outputFile, "      </model_scaling>\n");
    }    
  }  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// We need an array with all possible paths through the lattice. Once we have this array (filled
/// with AM and LM scores) we can rescore very quickly!
/////////////////////////////////////////////////////////////////////////////////////////////////////

int LexicalTree::createLatticeLMRescoring(LatticeNode *l,NBestList *scoreList,float amTot,float lmTot,
                                          int nrWords,int noSilWords, int *lmHist, int sentenceID)
{
  if(l != NULL)
  {
    if(l->wordID >= 0)
    {
      nrWords++;      
      if(l->wordID != startOfSentenceWord)
      {
        noSilWords++;
        if(languageModel != NULL)
        {
          lmTot += languageModel->getP(l->wordID,lmHist,lmHist);
        }
      }
    }  
    WLRList *wl = l->outArcs;
    if(wl == NULL)
    {
      scoreList[sentenceID].nrWords     = nrWords;
      scoreList[sentenceID].noSilWords  = noSilWords;
      scoreList[sentenceID].totAM       = amTot;
      scoreList[sentenceID].totLM       = lmTot;
      sentenceID++;
    }
    while(wl != NULL && sentenceID < MAX_ANALYSIS_NBEST_PATHS)
    {
      if(wl->wlr->lattice != NULL)
      {
        float am = amTot + wl->amScore;
        sentenceID = createLatticeLMRescoring(wl->wlr->lattice,scoreList,am,lmTot,nrWords,noSilWords, lmHist, sentenceID);
      }
      wl = wl->next;
    }
  }
  return sentenceID;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Determines how many paths there are through the lattice.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int LexicalTree::countLatticePaths(LatticeNode *l, int sentenceID)
{
  if(l != NULL)
  {
    WLRList *wl = l->outArcs;
    if(wl == NULL)
    {
      sentenceID++;
    }
    while(wl != NULL)
    {
      if(wl->wlr->lattice != NULL)
      {
        sentenceID = countLatticePaths(wl->wlr->lattice, sentenceID);
      }
      wl = wl->next;
    }
  }
  return sentenceID;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo DOCS
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::sortLatticePaths(LatticeNode *l)
{
  if(l != NULL)
  {
    // bubble sort:
    bool change = true;
    while(change)
    {
      change = false;
      WLRList *w0 = NULL;
      WLRList *w1 = l->outArcs;
      if(w1 != NULL)
      {
        WLRList *w2 = w1->next;
        while(w2 != NULL)
        {
          if(w2->wlr->COMBlikelihood > w1->wlr->COMBlikelihood)
          {
            if(w0 == NULL)
            {
              l->outArcs = w2;
            }
            else
            {
              w0->next = w2;
            }
            w1->next = w2->next;
            w2->next = w1;
            w0 = w2;
            w1 = w1;
            w2 = w1->next;
            change = true;
          }
          else
          {
            w0 = w1;
            w1 = w2;
            w2 = w2->next;
          }
        }
      }
    }
  }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will print the lattice output. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::printLattice(FILE *latFile, const char *label, int timeEnd)
{
  if(latticeAdmin == NULL)
  {
    createLattice();
  }
  if(latticeN == 0)
  {
    fprintf(latFile,"name %s\n", label);
    fprintf(latFile,"nodes 3 <s> <s> NULL\n");
    fprintf(latFile,"initial 0\n");
    fprintf(latFile,"final 2\n");
    fprintf(latFile,"transitions 2\n");
    fprintf(latFile,"0 1 lm=0.0,am=0.0,tot=0.0\n");
    fprintf(latFile,"1 2 lm=0.0,am=0.0,tot=0.0\n");
  }
  else
  {
    fprintf(latFile,"name %s\n", label);
    fprintf(latFile,"nodes %d <s> ",latticeN);
    LatticeNode *l = latticeAdmin->adminNext;
    while(l != NULL)
    {
      if(l->wordID < 0 || l->wordID == startOfSentenceWord)
      {
        fprintf(latFile,"NULL ");
      }
      else
      {
        fprintf(latFile,"%s ",vocabulary[l->wordID]);
      }
      l = l->adminNext;
    }
    fprintf(latFile,"\ninitial 0\n");
    fprintf(latFile,"final %d\n", latticeN-1);
    fprintf(latFile,"transitions %d\n",latticeL);
  
    
    l = latticeAdmin;  
    while(l != NULL)
    {
      WLRList *wl = l->inArcs;
      while(wl != NULL)
      {
        if(l->wordID < 0)
        {
          if(l->nodeNr == latticeN-1) // final node
          {
            fprintf(latFile,"%d %d %g %g %g\n",wl->wlr->lattice->nodeNr,l->nodeNr,wl->amScore,
                                               ((float)wl->wlr->timeStamp)/100.0,((float)wl->wlr->timeStamp)/100.0);
/*            fprintf(latFile,"%d %d lm=%g,am=%g,tot=%g NULL start-time=%g,end-time=%g\n",wl->wlr->lattice->nodeNr,l->nodeNr,
                  wl->lmScore, wl->amScore, wl->totScore,
                  ((float)wl->wlr->timeStamp)/100.0,((float)wl->wlr->timeStamp)/100.0);*/
          }
          else
          {
            fprintf(latFile,"%d %d %g %g %g\n",wl->wlr->lattice->nodeNr,l->nodeNr, wl->amScore,
                                               ((float)wl->wlr->timeStamp)/100.0,((float)l->timeEnd)/100.0);
/*            fprintf(latFile,"%d %d lm=%g,am=%g,tot=%g NULL start-time=%g,end-time=%g\n",wl->wlr->lattice->nodeNr,l->nodeNr,
                  wl->lmScore, wl->amScore, wl->totScore,
                  ((float)wl->wlr->timeStamp)/100.0,((float)l->timeEnd)/100.0);*/
          }
        }
        else
        {
          fprintf(latFile,"%d %d %g %g %g\n",wl->wlr->lattice->nodeNr,l->nodeNr, wl->amScore,
                                             ((float)wl->wlr->timeStamp)/100.0,((float)l->timeEnd)/100.0);
/*          fprintf(latFile,"%d %d lm=%g,am=%g,tot=%g %s start-time=%g,end-time=%g\n",wl->wlr->lattice->nodeNr,l->nodeNr,
                  wl->lmScore, wl->amScore, wl->totScore,
                  vocabulary[l->wordID], ((float)wl->wlr->timeStamp)/100.0,((float)l->timeEnd)/100.0);*/
        }
        wl = wl->next;
      }
      l = l->adminNext;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::latticeBaumWelch_initBackward(double *beta, int offset)
{
  LexicalNode *n = endNode;
  while(n != NULL)
  {
    beta[offset+n->compressedTreeIndex*3] = 1.0;
    n = n->parallel;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::latticeBaumWelch_setLikelihoods(LexicalNode *node, Vector *t, int time, int numberOfStates, double *latticeLikelihood)
{
  if(node == NULL)
  {
    node = grammarStart;
  }
  if(node != NULL)
  {
    int contextKey    = node->contextKey;
    int postIndex     = node->compressedTreeIndex*3;
    if(postIndex < 0)
    {
      postIndex *= -1;
    }
    int index         = time*numberOfStates     + postIndex;
    PhoneModel *m     = phoneModels[node->modelID];
    if(node->compressedTreeIndex >= 0)
    {
      if(node->parallel != NULL)
      {
        latticeBaumWelch_setLikelihoods(node->parallel,t,time,numberOfStates,latticeLikelihood);
      }
      int isSilOffset   = 0;
      latticeLikelihood[index] = m->getPDFProbability(contextKey,t,0,time);
      if(m->getStatistics()->isSil != 1)
      {
        latticeLikelihood[index+1] = m->getPDFProbability(contextKey,t,1,time);
        latticeLikelihood[index+2] = m->getPDFProbability(contextKey,t,2,time);
      }
      node->compressedTreeIndex *= -1;  // Mark this node so that the self-transition score is only added once!
      if(node->next != NULL)
      {
        latticeBaumWelch_setLikelihoods(node->next,t,time,numberOfStates,latticeLikelihood);
      }
      if(node->nextTree != NULL)
      {
        latticeBaumWelch_setLikelihoods(node->nextTree,t,time,numberOfStates,latticeLikelihood);
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

double LexicalTree::latticeBaumWelch_backward(double *latticeLikelihood, LexicalNode *node, double *beta, int time, int numberOfStates, double normFactor, double *resArray)
{
  if(node == NULL)
  {
    node = grammarStart;
  }
  if(node != NULL)
  {
    int postIndex     = node->compressedTreeIndex*3;
    if(node == endNode)
    {
      beta[time*numberOfStates + postIndex] = 0.0;
      return beta[(time+1)*numberOfStates + postIndex];
    }
    else
    {
      assert(node->compressedTreeIndex >= 0);
      if(resArray[node->compressedTreeIndex] < -1.0e20)
      {
        resArray[node->compressedTreeIndex] = 0.0;

        if(node->parallel != NULL)
        {
          latticeBaumWelch_backward(latticeLikelihood,node->parallel,beta,time,numberOfStates,normFactor,resArray);
        }
        PhoneModel *m     = phoneModels[node->modelID];
        int contextKey    = node->contextKey;
        double nextScore = 0.0;
        int index         = time*numberOfStates     + postIndex;
        int nextIndex     = (time+1)*numberOfStates + postIndex;
        int isSilOffset   = 0;
        if(m->getStatistics()->isSil != 1)
        {
          isSilOffset = 2;
        }
        if(node->next != NULL)
        {
          nextScore = latticeBaumWelch_backward(latticeLikelihood,node->next,beta,time,numberOfStates,normFactor,resArray);
        }
        if(node->nextTree != NULL)
        {
          nextScore = latticeBaumWelch_backward(latticeLikelihood,node->nextTree,beta,time,numberOfStates,normFactor,resArray);
        }

        double toSelf[3];
        double toNext[3];
        toSelf[0] = m->getTransition(contextKey,true, 0);
        toNext[0] = m->getTransition(contextKey,false,0);
        if(isSilOffset == 0)
        {
          beta[index] = (nextScore * toNext[0] + beta[nextIndex] * latticeLikelihood[index] * toSelf[0]) / normFactor;
        }
        else
        {
          toSelf[1] = m->getTransition(contextKey,true, 1);
          toNext[1] = m->getTransition(contextKey,false,1);
          toSelf[2] = m->getTransition(contextKey,true, 2);
          double lmProb = 1.0;
/*          if(node->wordID >= 0)
          {
            lmProb = exp(languageModel->getUnigram(node->wordID)/30.0);
          }*/
          toNext[2] = m->getTransition(contextKey,false,2) * lmProb;
          beta[index  ] = (beta[nextIndex+1] * latticeLikelihood[index+1] * toNext[0] + beta[nextIndex  ] * latticeLikelihood[index] * toSelf[0]) / normFactor;
          beta[index+1] = (beta[nextIndex+2] * latticeLikelihood[index+2] * toNext[1] + beta[nextIndex+1] * latticeLikelihood[index+1] * toSelf[1]) / normFactor;
          beta[index+2] = (nextScore * toNext[2] + beta[nextIndex+2] * latticeLikelihood[index+2] * toSelf[2]) / normFactor;
        }
        resArray[node->compressedTreeIndex] = beta[nextIndex] * latticeLikelihood[index];
        if(resArray[node->compressedTreeIndex] < -1.0e20)
        {
          resArray[node->compressedTreeIndex] = -1.0e19;
        }
      }
    }
    int number = 0;
    double res = 0.0;
    LexicalNode *par = node;
    while(par != NULL)
    {
      if(!(par->next == NULL && par->nextTree == NULL))
      {
        number++;
        res += resArray[par->compressedTreeIndex];
      }
      par = par->parallel;
    }
    if(number == 0)
    {
      return resArray[node->compressedTreeIndex];
    }
    res /= ((double)number);
    return res;
  }
  return -9.0e30;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::latticeBaumWelch_initForward(double *latticeBaumWelchAlfa)
{
  LexicalNode *n = grammarStart;
  while(n != NULL)
  {
    latticeBaumWelchAlfa[n->compressedTreeIndex*3] = 1.0;
    n = n->parallel;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::latticeBaumWelch_forward(double *latticeLikelihood, LexicalNode *node, double incomingScore,
                                           double *alfa, int time, int numberOfStates)
{
  if(node == NULL)
  {
    node = grammarStart;
  }
  if(node != NULL)
  {
    int contextKey    = node->contextKey;
    int postIndex     = node->compressedTreeIndex*3;
    if(postIndex < 0)
    {
      postIndex *= -1;
    }
    int index         = time*numberOfStates     + postIndex;
    PhoneModel *m     = phoneModels[node->modelID];
    alfa[index]   += incomingScore * latticeLikelihood[index];

    if(node->compressedTreeIndex >= 0)
    {
      if(node->parallel != NULL)
      {
        latticeBaumWelch_forward(latticeLikelihood,node->parallel,incomingScore,alfa,time,numberOfStates);
      }
      int previousIndex = (time-1)*numberOfStates + postIndex;
      int isSilOffset   = 0;
      if(m->getStatistics()->isSil != 1)
      {
        isSilOffset = 2;
      }

      double toSelf[3];
      double toNext[3];
      toSelf[0] = m->getTransition(contextKey,true, 0);
      toNext[0] = m->getTransition(contextKey,false,0);
      alfa[index]   += alfa[previousIndex] * latticeLikelihood[index] * toSelf[0];
      if(isSilOffset > 0)
      {
        toSelf[1] = m->getTransition(contextKey,true, 1);
        toNext[1] = m->getTransition(contextKey,false,1);
        toSelf[2] = m->getTransition(contextKey,true, 2);
        double lmProb = 1.0;
/*        if(node->wordID >= 0)
        {
          lmProb = exp(languageModel->getUnigram(node->wordID)/30.0);
        }*/
        toNext[2] = m->getTransition(contextKey,false,2) * lmProb;
        alfa[index+1] = alfa[previousIndex+1] * latticeLikelihood[index+1] * toSelf[1] + alfa[previousIndex]   * latticeLikelihood[index+1] * toNext[0];
        alfa[index+2] = alfa[previousIndex+2] * latticeLikelihood[index+2] * toSelf[2] + alfa[previousIndex+1] * latticeLikelihood[index+2] * toNext[1];
      }
      node->compressedTreeIndex *= -1;  // Mark this node so that the self-transition score is only added once!
      if(node->next != NULL)
      {
        latticeBaumWelch_forward(latticeLikelihood,node->next,alfa[previousIndex+isSilOffset] * toNext[isSilOffset],alfa,time,numberOfStates);
      }
      if(node->nextTree != NULL)
      {
        latticeBaumWelch_forward(latticeLikelihood,node->nextTree,alfa[previousIndex+isSilOffset] * toNext[isSilOffset],alfa,time,numberOfStates);
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::latticeBaumWelch_mmi_accumulatorsPosteriors(LexicalNode *node, double *posteriors, int numberOfStates, Vector *observation)
{
  if(node == NULL)
  {
    node = grammarStart;
  }
  if(node != NULL)
  {
    if(node->compressedTreeIndex >= 0)
    {
      int index = node->compressedTreeIndex*3;
      int contextKey    = node->contextKey;
      node->compressedTreeIndex *= -1;  // Mark this node so that the self-transition score is only added once!
      if(node->parallel != NULL)
      {
        latticeBaumWelch_mmi_accumulatorsPosteriors(node->parallel,posteriors,numberOfStates,observation);
      }
      for(int i=0;i<3;i++)
      {
        if(posteriors[index+i]>TRAIN_WEIGHT_MINIMUM && posteriors[index+i] <= 1.0)
        {
          phoneModels[node->modelID]->adapt_addAcumulatorData(i,node->contextKey,observation,posteriors[index+i]);
        }
      }
      if(node->next != NULL)
      {
        latticeBaumWelch_mmi_accumulatorsPosteriors(node->next,posteriors,numberOfStates,observation);
      }
      if(node->nextTree != NULL)
      {
        latticeBaumWelch_mmi_accumulatorsPosteriors(node->nextTree,posteriors,numberOfStates,observation);
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::latticeBaumWelch_printPosteriors (LexicalNode *node, double *posteriors, int time, int numberOfStates, int timeOffset)
{
  if(node == NULL)
  {
    node = grammarStart;
  }
  if(node != NULL)
  {
    if(node->compressedTreeIndex >= 0)
    {
      if(node->parallel != NULL)
      {
        latticeBaumWelch_printPosteriors(node->parallel,posteriors,time,numberOfStates, timeOffset);
      }
      int index = node->compressedTreeIndex*3;
      int contextKey    = node->contextKey;
      node->compressedTreeIndex *= -1;  // Mark this node so that the self-transition score is only added once!
      for(int i=0;i<3;i++)
      {
        if(posteriors[index+i]>TRAIN_WEIGHT_MINIMUM && posteriors[index+i] <= 1.0)
        {
          fprintf(outputFile, "%d\t%d\t%d\t%d\t%d\t%g\n",time+timeOffset,node->compressedTreeIndex*-1,i,node->modelID,contextKey,posteriors[index+i]);
        }
      }
      if(node->next != NULL)
      {
        latticeBaumWelch_printPosteriors(node->next,posteriors,time,numberOfStates, timeOffset);
      }
      if(node->nextTree != NULL)
      {
        latticeBaumWelch_printPosteriors(node->nextTree,posteriors,time,numberOfStates, timeOffset);
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////


void LexicalTree::latticeBaumWelch_calculatePosteriors (double *latticeLikelihood, LexicalNode *node, double incomingScore,
                                                        double *alfa, double *beta, double *posteriors, int time, int numberOfStates)
{
  if(node == NULL)
  {
    node = grammarStart;
  }
  if(node != NULL)
  {
    int postIndex     = node->compressedTreeIndex*3;
    if(postIndex < 0)
    {
      postIndex *= -1;
    }
    int index         = time*numberOfStates     + postIndex;
    PhoneModel *m     = phoneModels[node->modelID];
    int contextKey    = node->contextKey;

    posteriors[postIndex] += incomingScore * beta[index] * latticeLikelihood[index];


    if(node->compressedTreeIndex >= 0)
    {
      if(node->parallel != NULL)
      {
        latticeBaumWelch_calculatePosteriors(latticeLikelihood,node->parallel,incomingScore,alfa,beta,posteriors,time,numberOfStates);
      }
      int previousIndex = (time-1)*numberOfStates + postIndex;
      int isSilOffset   = 0;
      if(m->getStatistics()->isSil != 1)
      {
        isSilOffset = 2;
      }
      assert(postIndex >=0);
      assert(index >=0);
      assert(previousIndex >=0);

      double toSelf[3];
      double toNext[3];
      toSelf[0] = m->getTransition(contextKey,true, 0);
      toNext[0] = m->getTransition(contextKey,false,0);
      posteriors[postIndex] += alfa[previousIndex]*beta[index] * latticeLikelihood[index] * toSelf[0];
      if(isSilOffset > 0)
      {
        toSelf[1] = m->getTransition(contextKey,true, 1);
        toNext[1] = m->getTransition(contextKey,false,1);
        toSelf[2] = m->getTransition(contextKey,true, 2);
        double lmProb = 1.0;
/*        if(node->wordID >= 0)
        {
          lmProb = exp(languageModel->getUnigram(node->wordID)/30.0);
        }*/
        toNext[2] = m->getTransition(contextKey,false,2)*lmProb;
        posteriors[postIndex+1] += (alfa[previousIndex+1] * toSelf[1] + alfa[previousIndex] * toNext[0]) * beta[index+1] * latticeLikelihood[index+1];
        posteriors[postIndex+2] += (alfa[previousIndex+2] * toSelf[2] + alfa[previousIndex+1] * toNext[1]) * beta[index+2] * latticeLikelihood[index+2];
      }
      node->compressedTreeIndex *= -1;  // Mark this node so that the self-transition score is only added once!
      if(node->next != NULL)
      {
        latticeBaumWelch_calculatePosteriors(latticeLikelihood,node->next,alfa[previousIndex+isSilOffset]*toNext[isSilOffset],alfa,beta,posteriors,time,numberOfStates);
      }
      if(node->nextTree != NULL)
      {
        latticeBaumWelch_calculatePosteriors(latticeLikelihood,node->nextTree,alfa[previousIndex+isSilOffset]*toNext[isSilOffset],alfa,beta,posteriors,time,numberOfStates);
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

int LexicalTree::latticeBaumWelch_numberNodes(LexicalNode *node, int number, bool clear)
{
  int constant = 29101976;
  if(node == NULL)
  {
    node = grammarStart;
  }
  assert(node != NULL);

  if((!clear && ((node->compressedTreeIndex == constant) || (node->compressedTreeIndex == -1*constant))) ||
     (clear && (node->compressedTreeIndex != constant))                                                        )
  {
    if(clear)
    {
      node->compressedTreeIndex = constant;
    }
    else
    {
      node->compressedTreeIndex = number;  // using compressedTreeIndex because LMLA is not used during lattice rescoring anyway.
    }
    number++;
    if(node->parallel != NULL)
    {
      number = latticeBaumWelch_numberNodes(node->parallel,number,clear);
    }
    if(node->next != NULL)
    {
      number = latticeBaumWelch_numberNodes(node->next,number,clear);
    }
    if(node->nextTree != NULL)
    {
      number = latticeBaumWelch_numberNodes(node->nextTree,number,clear);
    }
  }
  assert(number < constant);
  return number;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will create a lattice structure for rescoring.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::setLattice(FILE *latFile)
{
  initialiseSystem();
  if(grammarStart != NULL)
  {
    if(latticeWordList != NULL)
    {
      LexicalNode *n2 = grammarStart->parallel;
      LexicalNode *n3 = grammarStart;
      while(n2 != NULL)
      {
        if(n2->depth != 2 && n2->depth != -2)
        {
          n2->inputToken    = new (TokenType*);
          *n2->inputToken   = NULL;
          n3->parallel = NULL;
        }
        else
        {
          n2->inputToken = n3->inputToken;
        }
        n3 = n2;
        n2 = n2->parallel;
      }
      grammarStart->parallel = NULL;
      for(int i=0;i<latticeWordListLength;i++)
      {
        deleteTree(latticeWordList[i],false);
      }
      delete[] latticeWordList;
      latticeWordList = NULL;
    }
    else
    {
      deleteTree(grammarStart,false);
    }
  }
  endNode           = NULL;
  grammarStart      = NULL;
  currentlyAligning = true;

  char label[MAXSTR];
  if(fscanf(latFile,"name %s\n", label) == 0)
  {
    USER_ERROR("Could not read lattice file. The first line should be: 'name [name]'");
  }
  latticeWordListLength = 0;
  if(fscanf(latFile,"nodes %d\n",&latticeWordListLength) == 0)
  {
    USER_ERROR("Could not read lattice file. The second line should be: 'nodes [number] [list of words]'");
  }

  latticeWordList = new LexicalNode*[latticeWordListLength];
  for(int i=0;i<latticeWordListLength;i++)
  {
    if(fscanf(latFile,"%s",label)==0)
    {
      USER_ERROR("Could not read all words in the lattice file (second line)");
    }
    if(i != 0 && strcmp(label,"NULL") == 0)
    {
      LexicalNode *n   = new LexicalNode;
      n->inputToken    = new (TokenType*);
      *n->inputToken   = NULL;
      n->tokenSeqLength=NONSIL_NUMBER_OF_STATES;
      n->tokenSeq      = new TokenType*[NONSIL_NUMBER_OF_STATES];
      for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
      {
        n->tokenSeq[ai]= NULL;
      }
      n->modelID       = 0;
      n->contextPrev   = 0;
      n->contextNext   = 0;
      n->contextKey    = 0;
      n->next          = NULL;
      n->nextTree      = NULL;
      n->parallel      = NULL;
      n->wordID        = endOfSentenceWord;
      n->depth         = -1*MAX_NDEPTH;
      n->compressedTreeIndex = -1;
      latticeWordList[i] = n;
    }
    else
    {
      if(i==0)
      {
        latticeWordList[i] = borrowPhoneString("<s>");
      }
      else
      {
        latticeWordList[i] = borrowPhoneString(label);
      }
      if(latticeWordList[i] == NULL)
      {
        USER_ERROR2("The following lattice word does not occur in the dictionary:",label);
      }
      LexicalNode *n = latticeWordList[i];
      LexicalNode *s = new LexicalNode;
      while(n->next != NULL)
      {
        n->depth = 0;
        n->compressedTreeIndex = -1;
        n = n->next;
      }
      n->depth = -1;
    }
  }
  char *resChar = fgets(label,MAXSTR,latFile);
  assert(resChar != NULL);

  int firstWord;
  if(fscanf(latFile,"initial %d\n", &firstWord) == 0)
  {
    USER_ERROR("Could not read lattice file. The third line should be: 'initial [number]'");
  }
  if(firstWord < 0 || firstWord >= latticeWordListLength)
  {
    USER_ERROR("The ID of the initial word is invalid");
  }
  if(latticeWordList[firstWord]->depth == -1*MAX_NDEPTH)
  {
    USER_ERROR("The first word cannot be a NULL word (except if it is word number zero).");
  }

  grammarStart = latticeWordList[firstWord];
  int finalWord;
  if(fscanf(latFile,"final %d\n", &finalWord) == 0)
  {
    USER_ERROR("Could not read lattice file. The fourth line should be: 'final [number]'");
  }
  if(finalWord < 0 || finalWord >= latticeWordListLength)
  {
    USER_ERROR("The ID of the final word is invalid");
  }
  if(latticeWordList[finalWord]->depth != -1*MAX_NDEPTH)
  {
    USER_ERROR("The final word in the lattice is not a 'NULL' node. I'm sorry, but for lattice recording this is a requirement...");
  }
  endNode = latticeWordList[finalWord];

  int nrTransitions;
  if(fscanf(latFile,"transitions %d\n", &nrTransitions) == 0)
  {
    USER_ERROR("Could not read lattice file. The fifth line should be: 'transitions [number]'");
  }

  for(int i=0;i<nrTransitions;i++)
  {
    int fromW = 0;
    int toW   = 0;
    if(fscanf(latFile,"%d %d",&fromW, &toW) == 0)
    {
      USER_ERROR("Could not read all transitions from the lattice file");
    }
    if(fromW < 0 || fromW >= latticeWordListLength || toW < 0 || toW >= latticeWordListLength)
    {
      printf("%04d:\t%d - %d\n",i,fromW,toW);
      USER_ERROR("One of the transitions is invalid");
    }
    char *resChar = fgets(label,MAXSTR,latFile);
    assert(resChar != NULL);
    if(fromW == finalWord)
    {
      // Ignore transitions from the final word back into the lattice!
    }
    else
    {
      // Now make a connection between two words!
      LexicalNode *fromNode = latticeWordList[fromW];
      while(fromNode->next != NULL)
      {
        fromNode = fromNode->next;
      }
      LexicalNode *toNode = latticeWordList[toW];
      if(fromNode->nextTree == NULL && toNode->contextPrev < 0)
      {
        fromNode->nextTree  = toNode;
        toNode->contextPrev = fromNode->modelID;
      }
      else
      {
        LexicalNode *n2   = new LexicalNode;
        n2->inputToken    = fromNode->inputToken;
        n2->compressedTreeIndex = -1;
        n2->tokenSeqLength=NONSIL_NUMBER_OF_STATES;
        n2->tokenSeq      = new TokenType*[NONSIL_NUMBER_OF_STATES];
        for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
        {
          n2->tokenSeq[ai]= NULL;
        }
        n2->nextTree      = NULL;  // will be overwritten below.
        n2->modelID       = fromNode->modelID;
        n2->next          = NULL;
        n2->parallel      = NULL;
        n2->wordID        = fromNode->wordID;
        n2->depth         = 2;
        n2->parallel      = fromNode->parallel;
        n2->contextPrev   = -1; //node->contextPrev;
        n2->contextNext   = -1; //node->contextNext;
        n2->contextKey    = -1;
        n2->nodeIsActive        = false;
        n2->toBeDeletedFromList = false;
        fromNode->parallel = n2;
  
        if(toW == finalWord)
        {
          n2->nextTree = toNode;
          n2->depth    = -2;
        }
        else
        {
          n2                = new LexicalNode;
          n2->compressedTreeIndex = -1;
          n2->inputToken    = new (TokenType*);
          *n2->inputToken   = NULL;
          n2->tokenSeqLength=NONSIL_NUMBER_OF_STATES;
          n2->tokenSeq      = new TokenType*[NONSIL_NUMBER_OF_STATES];
          for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
          {
            n2->tokenSeq[ai]= NULL;
          }
          n2->nextTree      = NULL;
          n2->next          = toNode;  // On purpose!! This will be changed to toNode->next later.
          n2->modelID       = toNode->modelID;
          n2->parallel      = NULL;
          n2->wordID        = toNode->wordID;
          n2->depth         = MAX_NDEPTH+1;
          n2->contextPrev   = -1; //node->contextPrev;
          n2->contextNext   = -1; //node->contextNext;
          n2->contextKey    = -1;
          n2->nodeIsActive        = false;
          n2->toBeDeletedFromList = false;
          fromNode->parallel->nextTree = n2;
        }
      }
    }
  }

  // Now set the correct values for the parallel pronunciation nodes. 
  // (do it now, because now we know the nextTree values: we might need to set multiple nextTree's)
  for(int i=0;i<latticeWordListLength;i++)
  {
    LexicalNode *pr = latticeWordList[i];
    if(pr!=NULL)
    {
      while(pr->next != NULL)
      {
        pr = pr->next;
      }
      LexicalNode *par = pr->parallel;
      while(par != NULL)
      {
        if(par->nextTree != NULL && par->nextTree->depth == MAX_NDEPTH+1)
        {
          par->nextTree->depth = -1;

          LexicalNode *nPar = par->nextTree;
          LexicalNode *n2   = nPar->next;
          nPar->next        = n2->next;     // First the next and nextTree pointers are set to the correct node.
          nPar->nextTree    = n2->nextTree;
          // Now check if there are parallel nodes that we need to copy:
          n2 = n2->parallel;
          while(n2 != NULL)
          {
            LexicalNode *n3   = new LexicalNode;
            n3->compressedTreeIndex = -1;
            n3->inputToken    = nPar->inputToken;
            n3->tokenSeqLength=NONSIL_NUMBER_OF_STATES;
            n3->tokenSeq      = new TokenType*[NONSIL_NUMBER_OF_STATES];
            for(int ai=0;ai<NONSIL_NUMBER_OF_STATES;ai++)
            {
              n3->tokenSeq[ai]= NULL;
            }
            n3->nextTree      = n2->nextTree;
            n3->modelID       = n2->modelID;
            n3->next          = NULL;
            n3->parallel      = nPar->parallel;
            nPar->parallel    = n3;
            n3->wordID        = n2->wordID;
            n3->depth         = -1;
            n3->contextPrev   = -1; //node->contextPrev;
            n3->contextNext   = -1; //node->contextNext;
            n3->contextKey    = -1;
            n3->nodeIsActive        = false;
            n3->toBeDeletedFromList = false;

            n2 = n2->parallel;
          }
          nPar->depth    = -1;
        }
        par = par->parallel;
      }
    }
  }

  // Setting the context:
  assert(grammarStartContext >= 0);
  setNodeContext(grammarStart,grammarStartContext);

  assert(grammarStart != NULL);
  assert(endNode != NULL);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::setNodeContext(LexicalNode *node,int leftContext)
{
  if(node != NULL && node->contextKey < 0)
  {
    setNodeContext(node->parallel,leftContext);
    node->contextPrev = leftContext;
    if(node->next != NULL)
    {
      node->contextNext = node->next->modelID;
    }
    else if(node->nextTree != NULL)
    {
      node->contextNext = node->nextTree->modelID;
    }
    else
    {
      node->contextNext = 0;
    }
    node->contextKey = node->contextPrev*numberOfPhones+node->contextNext;
    if(node->next != NULL)
    {
      setNodeContext(node->next,node->modelID);
    }
    if(node->nextTree != NULL)
    {
      setNodeContext(node->nextTree,node->modelID);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will print the N-Best list output.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::printNBestList(FILE *nbestFile, LatticeNode *l)
{
  if(latticeAdmin == NULL)
  {
    createLattice();
  }
 // The first in the list is the first in the lattice:
  if(l == NULL)
  {
    l = latticeAdmin;
  }  
  if(l != NULL)
  {
    if(l->wordID >= 0)
    {
      if(nbestFile == NULL)
      {
        fprintf(outputFile, "%s ",vocabulary[l->wordID]);
      }
      else
      {
        fprintf(nbestFile,"%s ",vocabulary[l->wordID]);
      }
    }
    WLRList *wl = l->outArcs;
    if(wl == NULL)
    {
      if(nbestFile == NULL)
      {
        fprintf(outputFile, "\n");
      }
      else
      {
        fprintf(nbestFile,"\n");
      }
    }
    while(wl != NULL)
    {
      if(wl->wlr->lattice != NULL)
      {
        printNBestList(nbestFile,wl->wlr->lattice);
      }
      wl = wl->next;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will return the total score of the best recognition (used for VTLN)
/////////////////////////////////////////////////////////////////////////////////////////////////////

double LexicalTree::getBestRecognitionScore()
{
  TokenType *useToken = findBestToken(false, false, false,false);
  if(useToken != NULL && useToken->path!=NULL)
  {
    return ((double) useToken->path->COMBlikelihood);
  }
  return -9.0e300;
}

#ifdef TEST_ARTICULATORY
void LexicalTree::testArticulatory(ArticulatoryStream *s)
{
  myArtStream = s;
}
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo Docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

int LexicalTree::getBestIDSequence(int *idList, int maxLength, bool showSil, bool notFinishedAllowed)
{
  int realLength = 0;
  TokenType *useToken = findBestToken(false, false, false, notFinishedAllowed);
  if(useToken != NULL && useToken->path!=NULL)
  {
    WLRType *wlr = useToken->path;
    while((wlr != NULL) && (wlr->timeStamp >= startTime) && (realLength < maxLength))
    {
      int wordID = -1;
      if(showSil || wlr->isSil != 1)
      {
        if(wlr->isSil)
        {
          wordID = startOfSentenceWord;
        }
        else if(languageModel != NULL)
        {
          wordID = languageModel->getLastWordID(wlr->lmHistory);
        }
        else
        {
          wordID = wlr->lmHistory[0];
        }
        idList[realLength] = wordID;
        realLength++;
      }
      wlr = wlr->previous;
    }
  }
  return realLength;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will find the best found recognition and print it to standard output.
/// If 'complete' is true, the string will also contain timing and likelihood data.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::getBestRecognition(bool addEndOfSentence, bool outputXML, bool complete, const char *label,
                                     int milliSec, int totLength, const char *beginTime, const char *endTime)
{
  if(complete)
  {
    if(outputXML)
    {
      fprintf(outputFile, "  <speech label=\"%s\" ",label);
      if(beginTime != NULL && endTime != NULL)
      {
        fprintf(outputFile, "begintime=\"%s\" ",beginTime);
        fprintf(outputFile, "endtime=\"%s\" ",endTime);        
      }
      fprintf(outputFile, ">\n");
      if(milliSec != 0 && totLength != 0)
      {
        fprintf(outputFile, "    <real_time milliseconds=\"%d\" ",milliSec);
        fprintf(outputFile, "frames=\"%d\" ",totLength);
        fprintf(outputFile, "RTF=\"%.4f\"/>\n",((double)milliSec/10.0)/((double)totLength));    
      }
    }
    else
    {
      fprintf(outputFile, "###################################################-| %s |-#################################\n",label);    
      if(beginTime != NULL && endTime != NULL)
      {
        fprintf(outputFile, "## Begin time: %s\n## End time:   %s\n",beginTime,endTime);
      }
      if(milliSec != 0 && totLength != 0)
      {
        fprintf(outputFile, "## Sentence processed in: %d milliseconds / %d frames = %2.4f times realtime.\n",
                          milliSec,(int)totLength, ((float)milliSec/10.0)/((float)totLength));      
      }
    }
  }
  TokenType *useToken = findBestToken(addEndOfSentence, complete, outputXML);
  if(useToken != NULL && useToken->path!=NULL)
  {   
    bestRecPath = useToken->path;
    if(complete)
    {    
       // If LMLA is switched off, the lookAheadV will be zero always..
      if(outputXML)
      {
        fprintf(outputFile, "    <score COMBINED=\"%.5f\"/>\n",useToken->likelihood - useToken->lookAheadV);
      }
      else
      {
        fprintf(outputFile, "## Sentence score     : %10.5f\n",useToken->likelihood - useToken->lookAheadV);
      }

#ifdef SEARCH_STATISTICS_ON      
      float nr = (float)(sentenceStats.nrOfCleanUps);
      if(outputXML)
      {
        fprintf(outputFile, "    <search_statistics>\n");    
        fprintf(outputFile, "      <LMLA_TABLES AVG=\"%.2f\" MAX=\"%d\"/>\n",
                      sentenceStats.nrOfLMLATables/nr,sentenceStats.nrOfLMLATables_Max);
        fprintf(outputFile, "      <ACTIVE_NODES AVG=\"%.2f\" MAX=\"%d\"/>\n",
                      sentenceStats.nrOfActiveNodes/nr,sentenceStats.nrOfActiveNodes_Max);
        fprintf(outputFile, "      <ACTIVE_TOKENS AVG=\"%.2f\" MAX=\"%d\"/>\n",
                      sentenceStats.nrOfActiveTokens/nr,sentenceStats.nrOfActiveTokens_Max);
        fprintf(outputFile, "    </search_statistics>\n");    
      }
      else
      {      
        fprintf(outputFile, "## LMLA tables   (avg,max) : %10.2f %6d\n",sentenceStats.nrOfLMLATables/nr,sentenceStats.nrOfLMLATables_Max);
        fprintf(outputFile, "## Active nodes  (avg,max) : %10.2f %6d\n",sentenceStats.nrOfActiveNodes/nr,sentenceStats.nrOfActiveNodes_Max);
        fprintf(outputFile, "## Active tokens (avg,max) : %10.2f %6d\n",sentenceStats.nrOfActiveTokens/nr,sentenceStats.nrOfActiveTokens_Max);
        
        if(settings.prune_Hist > 0)
        {
          fprintf(outputFile, "## Histogram prunes :       %6d\n", sentenceStats.nrOfHistPruning);
        }
        fprintf(outputFile, "## Histogram state prunes : %6d\n", nrOfStateHistPruning);        
        fprintf(outputFile, "## ----------------- ANALYSE ---------------------------------\n");
      }
      
#endif


      if(analysisSettings != NULL)
      {
        if(outputXML)
        {
          fprintf(outputFile, "      <analyse ");          
          fprintf(outputFile, "MIN_NEEDED_HISTOGRAM_PRUNING=\"%d\" ", analysisSettings->maxHistogram);
          fprintf(outputFile, "MIN_NEEDED_STATE_HISTOGRAM_PRUNING=\"%d\" ", analysisSettings->maxStateHistogram);
          fprintf(outputFile, "MIN_NEEDED_END_STATE_BEAM=\"%g\" ", analysisSettings->maxEndStateBeam);      
          fprintf(outputFile, "MIN_NEEDED_BEAM=\"%g\" ", analysisSettings->maxBeam);
          fprintf(outputFile, "MIN_NEEDED_STATE_BEAM=\"%g\" ", analysisSettings->maxStateBeam);      
          fprintf(outputFile, "NUMBER_WORDS_IN_REF_NON_SIL=\"%d\">\n", analysisSettings->refStats.noSilWords);
//          printLMParStats(outputXML);
          errorAnalysis(useToken->path,0, outputXML);
          fprintf(outputFile, "      </analyse>\n");
        }
        else
        {              
          fprintf(outputFile, "## Maximum needed HISTOGRAM PRUNING VALUE          : %d\n", analysisSettings->maxHistogram);
          fprintf(outputFile, "## Maximum needed STATE HISTOGRAM PRUNING VALUE    : %d\n", analysisSettings->maxStateHistogram);
          fprintf(outputFile, "## Maximum needed END STATE BEAM VALUE             : %g\n", analysisSettings->maxEndStateBeam);      
          fprintf(outputFile, "## Maximum needed BEAM VALUE                       : %g\n", analysisSettings->maxBeam);
          fprintf(outputFile, "## Maximum needed STATE BEAM VALUE                 : %g\n", analysisSettings->maxStateBeam);      
          fprintf(outputFile, "## Number of words in reference (non-SIL)          : %d\n", analysisSettings->refStats.noSilWords);
//          printLMParStats(outputXML);
          errorAnalysis(useToken->path,0, outputXML);
          fprintf(outputFile, "## ----------------- HYPOTHESIS ------------------------------\n");
        }
      }
      if(outputXML)
      {
        fprintf(outputFile, "    <wordsequence>\n");            
      }
      else
      {
        fprintf(outputFile, "#                                          Timestamp        TotAM           TotLM         TotComb\n");      
      }
    }    
    getBestPath(useToken->path,outputXML,complete);    
    if(outputXML)
    {
      fprintf(outputFile, "    </wordsequence>\n");    
    }
    else if(!complete)
    {
      fprintf(outputFile, " (%s)\n",label);
    }
#ifdef TEST_ARTICULATORY
    if(complete)
    {
      myArtStream->testPrint(outputFile);
    }
#endif

  }
  else
  {
//    fprintf(outputFile, "NO RECOGNITION!!!");
  }
  if(outputXML)
  {
    fprintf(outputFile, "  </speech>\n");
  }
  fflush(outputFile);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the ID of the word 'word'. If the word is not in the vocabulary (OOV), -1 is returned.
/////////////////////////////////////////////////////////////////////////////////////////////////////

int LexicalTree::getWordID(const char *word)
{
  if(strlen(word) == 0)
  {
    return -1;
  } 
  int res = 0;
  while((res<numberOfWords) && (strcmp(vocabulary[res],word)!=0))
  {
    res++;
  }
  if(res<numberOfWords)  
  {
    return res;
  }
  else
  {
    return -1;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Given a wordID, the corresponding word string is returned. If the wordID is not valid, NULL is returned.
/////////////////////////////////////////////////////////////////////////////////////////////////////

char *LexicalTree::getWord(int wordID)
{
  if(wordID < numberOfWords && wordID >= 0)
  {
    return vocabulary[wordID];
  }
  return NULL;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the number of words in the current vocabulary (PPT)
/////////////////////////////////////////////////////////////////////////////////////////////////////

int LexicalTree::getNumberOfWords()
{
  return numberOfWords;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Print all words followed with their pronunciation that are reachable from node 'node'. If the 
/// root node is used, all words from the PPT are printed.
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool LexicalTree::printWordPronunciation(LexicalNode *node, int wordID)
{
  bool printLater = false;
  bool printNow = false;
  if(node->wordID == wordID)
  {
    fprintf(outputFile, "%s\t", getWord(node->wordID));
    fflush(outputFile);
    printNow = true;
  }  
  if(node->parallel)
  {
    printLater = printWordPronunciation(node->parallel,wordID);
  }
  if(node->next)
  {
    printNow = printNow || printWordPronunciation(node->next,wordID);
  }
  if(printNow)
  {
    fprintf(outputFile, "%s ",phoneModels[node->modelID]->getStatistics()->name);
    fflush(outputFile);
  }
  return (printNow || printLater);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns a string of all phones in the alignment task. 
/////////////////////////////////////////////////////////////////////////////////////////////////////

char *LexicalTree::getAlignmentString()
{
  static char res[5000];
  char tmp[20];
  strcpy(res,"");
  LexicalNode *node = grammarStart;
  while(node != NULL)
  {
    if(node->modelID >= 0)
    {      
      strcpy(tmp,phoneModels[node->modelID]->getStatistics()->name);
      int i=0;
      while(i<19 && tmp[i] != ' ')
      {
        i++;
      }
      tmp[i] = 0;
      strcat(res,tmp);
      strcat(res," ");
    }
    if(node->next == NULL)
    {
      node = node->nextTree;
    }
    else
    {
      node = node->next;
    }
  }
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method will prune away all LMLA tables that are currently not being used. An LMLA is not used
/// when there is no token in the PPT with the same LM history as the history of the LMLA table.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::pruneLMLA()
{
  // If this method is called, LMLA is switched on!
    int deleted = 0;
    int tel     = 0;
    for(int i=0;i<LMLA_CACHE;i++)
    {
      LMLAGlobalListType *l = lmla_list[i];
      LMLAGlobalListType *l2 = NULL;
      while(l!=NULL)
      {
        if(l->collissionTime == 0)
        {
          if(l2 == NULL)
          {
            l2 = l;
            l = l->next;
#ifndef LMLA_UNIGRAM_ONLY
            if(l2->lookAhead != NULL)
            {
              delete[] l2->lookAhead;
            }
#endif
            delete l2;
#ifdef SEARCH_STATISTICS_ON            
            sentenceStats.nrOfLMLATables_Temp--;
#endif            
            deleted++;
            l2 = NULL;
            lmla_list[i] = l;
          }
          else
          {
            LMLAGlobalListType *l3 = l;
            l2->next = l->next;
            l = l->next;
#ifndef LMLA_UNIGRAM_ONLY
            if(l3->lookAhead != NULL)
            {
              delete[] l3->lookAhead;
            }
#endif
            delete l3;
#ifdef SEARCH_STATISTICS_ON            
            sentenceStats.nrOfLMLATables_Temp--;
#endif            
            deleted++;
          }
        }
        else
        {
          tel++;
          l2 = l;
          l->collissionTime = 0;
          l  = l->next;
        }
      }
    }
//    fprintf(outputFile, "Deleted,still here: %d  %d\n",deleted,tel);      
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method calculates a hash key for an LMLA table.
///
/// A simple modulo function is used: modulo-LMLA_CACHE(sum(lmHistory[i]))
/////////////////////////////////////////////////////////////////////////////////////////////////////

int LexicalTree::getLMLAHashKey(int *lmHistory)
{
  // If this method is called, LMLA is switched on!
  int res = 0;
  for(int i=0;i<LM_NGRAM_DEPTH;i++)
  {
    res += lmHistory[i];
  }
  int MODULO(key,res,LMLA_CACHE);
  return key;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method creates a new LMLA table and adds it to the global LMLA list.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::createLMLAs(LMLAGlobalListType *lmlaGlobal, float *allLMP)
{
#ifdef MAX_NR_THREADS
  pthread_mutex_lock( &lmla_condition_mutexDone );
  while(!lmla_threadFinished)
  {
    pthread_cond_wait( &lmla_condition_threadDone, &lmla_condition_mutexDone );
  }
  lmla_threadFinished = false;
  pthread_mutex_unlock( &lmla_condition_mutexDone );

  lmla_dataThread.lmlaGlobal              = lmlaGlobal;
  lmla_dataThread.allLMP                  = allLMP;
  // start thread
  pthread_mutex_lock( &lmla_condition_mutexStart );
  lmla_threadMayStart = true;
  pthread_cond_signal( &lmla_condition_threadStart );
  pthread_mutex_unlock( &lmla_condition_mutexStart );

#else
  float fastCompressedTreeBackValue[numberOfCompressedNodes];
  float *bv              = fastCompressedTreeBackValue;
  float prefRes = 0.0;
  for(int i=0;i<numberOfCompressedNodes;i++)
  {
    float res = SMALLEST_NEGATIVE;
    if(fastCompressedTree[i].wordID == -10)
    {
      res = 0.0;
    }
    else if(fastCompressedTree[i].wordID >= 0)
    {
      res = allLMP[fastCompressedTree[i].wordID];
    }
    else if(fastCompressedTree[i].nextIndex >= 0)
    {
      res = fastCompressedTreeBackValue[fastCompressedTree[i].nextIndex];
    }
    lmlaGlobal->lookAhead[i] = res;

    if(fastCompressedTree[i].hasPar && prefRes > res)
    {
      res = prefRes;
    }
    prefRes    = res;
    *bv        = res;
    bv++;
  }
  delete[] allLMP;
  COUNTER_MINUS(lmlaCount);
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Sets the token-distribution file handler for the current utterence. 
/// Each time frame of a recognition (not of a forced alignment), the distribution file will be 
/// updated.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::setTokenDistributionFile(FILE *tdF)
{
  if(tdFile != NULL)
  {
    fprintf(tdFile,"</token-distribution>\n");  
    fclose(tdFile);
  }
  tdFile = tdF;
  if(tdFile != NULL)
  {
    fprintf(tdFile,"<token-distribution numberOfPhones=\"%d\" numberOfColumns=\"%d\">\n",numberOfPhones,
                                                                                         TOKENDISTRIBUTION_TREEDEPTH);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Sets the confidence array of the phone loop. This is used for confidence measuring.
/// Currently, the confidence is only printed on the screen. Analysis is done externally...
/////////////////////////////////////////////////////////////////////////////////////////////////////


void LexicalTree::setPhoneLoopConfidence(float *phoneConf, int offset)
{
  phoneLoopConfidence       = phoneConf;
  phoneLoopConfidenceOffset = offset;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Given the LM history 'lmHistory' a pointer to the corresponding LMLA table is returned. If the
/// table does not exist, it is created.
/////////////////////////////////////////////////////////////////////////////////////////////////////

LMLAGlobalListType *LexicalTree::getLMLATable(int *lmHistory, bool onlyPrepare)
{
  // If this method is called, LMLA is switched on!
  int hashKey = getLMLAHashKey(lmHistory);
  LMLAGlobalListType *lglob = lmla_list[hashKey];
  while(lglob != NULL && !(compareLMHistory(lglob->wordHistory, lmHistory)))
  {
    lglob = lglob->next;
  }
  if(lglob == NULL || lglob->lookAhead == NULL) // Does not exist yet! Let's make a new one...
  {
    if(lglob == NULL)
    {
  #ifdef SEARCH_STATISTICS_ON
      sentenceStats.nrOfLMLATables_Temp++;
  #endif
      lglob                 = new LMLAGlobalListType;
      lglob->next           = lmla_list[hashKey];
      copyLMHistory(lglob->wordHistory, lmHistory);
      lmla_list[hashKey]    = lglob;
    }
    if(onlyPrepare)
    {
      lglob->lookAhead    = NULL;
    }
    else
    {
#ifdef LMLA_UNIGRAM_ONLY
      if(lookAheadUnigram == NULL)
      {
        lookAheadUnigram = new LMLAGlobalListType;
        lookAheadUnigram->lookAhead    = new float[numberOfCompressedNodes];
        float *allLMP = new float[numberOfWords];
        COUNTER_PLUS(lmlaCount);
        languageModel->getAllP(NULL,allLMP);
        for(int a=0;a<numberOfWords;a++)
        {
          allLMP[a] = allLMP[a]*settings.weights_LmScale + transitionPenalty[a];
        }
        createLMLAs(lookAheadUnigram,allLMP);
      }
      lglob->lookAhead = lookAheadUnigram->lookAhead;
#else
      lglob->lookAhead    = new float[numberOfCompressedNodes];
      // Now do the lookahead for all nodes:
      float *allLMP = new float[numberOfWords];
      COUNTER_PLUS(lmlaCount);
      languageModel->getAllP(lmHistory,allLMP);
      for(int a=0;a<numberOfWords;a++)
      {
        allLMP[a] = allLMP[a]*settings.weights_LmScale + transitionPenalty[a];
      }
      createLMLAs(lglob,allLMP);
      // allLMP will be deleted by createLMLAs...
#endif
    }
  }
  lglob->collissionTime = 1;
  return lglob;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method updates the global lexical tree statistics. Global statistiscs are statistics on
/// all recognition tasks in one session (like average amount of nodes needed for all utterences).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::updateGlobalStats()
{
#ifdef SEARCH_STATISTICS_ON
  globalStats.nrOfCleanUps          += sentenceStats.nrOfCleanUps;
  globalStats.nrOfLMLATables        += sentenceStats.nrOfLMLATables;
  globalStats.nrOfActiveNodes       += sentenceStats.nrOfActiveNodes;
  globalStats.nrOfActiveTokens      += sentenceStats.nrOfActiveTokens;
  
  if(sentenceStats.nrOfLMLATables_Max > globalStats.nrOfLMLATables_Max)
  {
    globalStats.nrOfLMLATables_Max = sentenceStats.nrOfLMLATables_Max;
  }
  if(sentenceStats.nrOfActiveNodes_Max > globalStats.nrOfActiveNodes_Max)
  {
    globalStats.nrOfActiveNodes_Max = sentenceStats.nrOfActiveNodes_Max;
  }
  if(sentenceStats.nrOfActiveTokens_Max > globalStats.nrOfActiveTokens_Max)
  {
    globalStats.nrOfActiveTokens_Max = sentenceStats.nrOfActiveTokens_Max;
  }
  
  globalStats.nrOfLMLATables_Prev   = globalStats.nrOfLMLATables;
  globalStats.nrOfActiveNodes_Prev  = globalStats.nrOfActiveNodes;
  globalStats.nrOfActiveTokens_Prev = globalStats.nrOfActiveTokens;  
#endif  
}

#ifdef SEARCH_STATISTICS_ON

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// This method updates the internal lexical tree statistics. Internal statistiscs are statistics on
/// one recognition task of a session (like average amount of nodes needed for the current utterence).
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::updateStats()
{
  LexicalNodeList *l = nodeList;   
  while(l != NULL)
  {
    TokenType *t = (*l->node->inputToken);
    sentenceStats.nrOfActiveNodes++;
    while(t != NULL)
    {
      t = t->next;
      sentenceStats.nrOfActiveTokens++;
    }   
    for(int i=0;i<l->node->tokenSeqLength;i++)
    { 
      t = l->node->tokenSeq[i];
      while(t != NULL)
      {
        t = t->next;
        sentenceStats.nrOfActiveTokens++;
      }   
    }
    l = l->next;
  }

  sentenceStats.nrOfCleanUps++;
  
  if(sentenceStats.nrOfLMLATables_Temp > sentenceStats.nrOfLMLATables_Max)
  {
    sentenceStats.nrOfLMLATables_Max = sentenceStats.nrOfLMLATables_Temp;
  }
  sentenceStats.nrOfLMLATables += sentenceStats.nrOfLMLATables_Temp;
  int nr = sentenceStats.nrOfActiveNodes - sentenceStats.nrOfActiveNodes_Prev;
  if(nr > sentenceStats.nrOfActiveNodes_Max)
  {
    sentenceStats.nrOfActiveNodes_Max = nr;
  }
  nr = sentenceStats.nrOfActiveTokens - sentenceStats.nrOfActiveTokens_Prev;
  if(nr > sentenceStats.nrOfActiveTokens_Max)
  {
    sentenceStats.nrOfActiveTokens_Max = nr;
  }
  
  sentenceStats.nrOfLMLATables_Prev   = sentenceStats.nrOfLMLATables;
  sentenceStats.nrOfActiveNodes_Prev  = sentenceStats.nrOfActiveNodes;
  sentenceStats.nrOfActiveTokens_Prev = sentenceStats.nrOfActiveTokens;  
}

#endif    


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the node in the lexical tree that correstponds to the PLRType string *pt.
/// wordID should only be > 0 if the phone pt is the last phone of the word.
/////////////////////////////////////////////////////////////////////////////////////////////////////

LexicalNode *LexicalTree::findCorrectNode(PLRType *pt, int wordID)
{
  LexicalNode *res = NULL;
  if(pt != NULL)
  {
    int prev = (int)(pt->contextKey/numberOfPhones);
    
    if(wordID >= 0) // We are the last phone:
    {
      res = treeEnd[pt->phoneID*numberOfPhones + prev];
      while(res != NULL && res->contextKey != pt->contextKey)
      {
        res = res->parallel;
      }     
    } 
    else
    {  
      if(pt->previous == NULL)  // First node: let's find 'm!
      {
        res = treeStart[pt->phoneID + prev*numberOfPhones];
        while(res != NULL && res->contextKey != pt->contextKey)
        {
          res = res->parallel;
        }
        assert(res != NULL);
      }
      else                      // Not the first. Let's find the correct one first...
      {
        res = findCorrectNode(pt->previous, -1);
        assert(res != NULL);
        res = res->next;      
        while(res != NULL && (res->contextKey != pt->contextKey))
        {
          res = res->parallel;
        }
      }    
    }
  }  
  return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// Returns the phone alignment. Called by getBestPath() when the configuration switch "phone-based"
/// is defined.
/////////////////////////////////////////////////////////////////////////////////////////////////////

void LexicalTree::getPhoneAlignment(const char *prefix, PLRType *pt, int lastFrame, bool confident)
{
  if(pt != NULL)
  {
    char string[10];
    getPhoneAlignment(prefix, pt->previous,pt->timeStamp - 1, confident);
    if(pt->phoneID != -1)
    {
      if(prefix!=NULL)
      {
        if(confident)
        {
          fprintf(outputFile, "%s          #Phone#  ",prefix);
        }
        else
        {
          fprintf(outputFile, "%s*         #Phone#  ",prefix);
        }
      }
      else
      {
        if(confident)
        {
          fprintf(outputFile, "          <phone confident=\"1\" phoneID=\"%d\" phoneStr=\"",pt->phoneID);
        }
        else
        {
          fprintf(outputFile, "          <phone confident=\"0\" phoneID=\"%d\" phoneStr=\"",pt->phoneID);
        }
      }
      int prev = (int)(pt->contextKey/numberOfPhones);
      int next = (int)(pt->contextKey - prev*numberOfPhones);      
      if(prev >= 0)
      {
        strncpy(string,phoneModels[prev]->getStatistics()->name,3);
        if(prefix==NULL)
        {
          for(int i=0;i<3;i++)
          {
            if(string[i]==' ')
            {
              string[i] = 0;
            }
          }
        }
        string[3]=0;
        fprintf(outputFile, "%s ",string);
      }
      else
      {
        fprintf(outputFile, "!!! ");
      }
      strncpy(string,phoneModels[pt->phoneID]->getStatistics()->name,3);
      if(prefix==NULL)
      {
        for(int i=0;i<3;i++)
        {
          if(string[i]==' ')
          {
            string[i] = 0;
          }
        }
      }
      string[3]=0;      
      fprintf(outputFile, "%s ",string);
      if(next >= 0)
      {
        strncpy(string,phoneModels[next]->getStatistics()->name,3);
        if(prefix==NULL)
        {
          for(int i=0;i<3;i++)
          {
            if(string[i]==' ')
            {
              string[i] = 0;
            }
          }
        }
        string[3]=0;        
        fprintf(outputFile, "%s ",string);
      }
      else
      {
        fprintf(outputFile, "!!! ");
      }
      if(prefix!=NULL)
      {
        fprintf(outputFile, "           ");            
      }
      else
      {
        fprintf(outputFile, "\" ");
      }
      float phoneConf = 0.0;
      if(phoneLoopConfidence != NULL)
      {
        phoneConf = phoneLoopConfidence[lastFrame - phoneLoopConfidenceOffset];
        if(pt->timeStamp > 0)
        {
          phoneConf -= phoneLoopConfidence[((int)pt->timeStamp-1) - phoneLoopConfidenceOffset];
        }
      }
      if(phoneModels[pt->phoneID]->getStatistics()->isSil == 1.0)
      {
        if(prefix!=NULL)
        {
          fprintf(outputFile, "| %5d %5d | %5d %5d | %5d %5d |   %12.5f  %12.5f\n",
                                              (int)pt->timeStamp,lastFrame,
                                              -1,-1,
                                              -1,-1,
                                               pt->likelihood, phoneConf);      
        }
        else
        {
          fprintf(outputFile, "beginFr=\"%d\" ",(int)pt->timeStamp);
          fprintf(outputFile, "endFr=\"%d\" ", lastFrame);
          fprintf(outputFile, "AML=\"%.4f\" ", pt->likelihood);
          fprintf(outputFile, "confidence=\"%.4f\"/>\n", phoneConf);
        }                                               
      }
      else
      {
        if(prefix!=NULL)
        {
          fprintf(outputFile, "| %5d %5d | %5d %5d | %5d %5d |   %12.5f  %12.5f\n",
                                              (int)pt->timeStamp,(int)(pt->timeStamp+pt->stateOffset[0]-1),
                                              (int)(pt->timeStamp+pt->stateOffset[0]),(int)(pt->timeStamp+pt->stateOffset[1]-1),
                                              (int)(pt->timeStamp+pt->stateOffset[1]),lastFrame,
                                                    pt->likelihood,phoneConf);
        }
        else
        {
          fprintf(outputFile, "beginFr=\"%d\" ",(int)pt->timeStamp);
          fprintf(outputFile, "beginFrMid=\"%d\" ",(int)(pt->timeStamp+pt->stateOffset[0]));
          fprintf(outputFile, "endFrMid=\"%d\" ", (int)(pt->timeStamp+pt->stateOffset[1]-1));
          fprintf(outputFile, "endFr=\"%d\" ", lastFrame);
          fprintf(outputFile, "AML=\"%.4f\" ", pt->likelihood);
          fprintf(outputFile, "confidence=\"%.4f\"/>\n", phoneConf);
        }
      }
#ifdef TEST_ARTICULATORY
      if(pt->timeStamp >= 0 && phoneModels[pt->phoneID]->getStatistics()->isSil != 1.0)
      {
        myArtStream->testStream((int)pt->phoneID,(int)pt->timeStamp,(int)lastFrame,&settings, pt->contextKey);
//        myArtStream->testStream((int)pt->phoneID,(int)(pt->timeStamp+pt->stateOffset[0]),(int)(pt->timeStamp+pt->stateOffset[1]-1));
      }
#endif

    }
  }
}

