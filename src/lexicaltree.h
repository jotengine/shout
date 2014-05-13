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
#ifndef LEXICALTREE_H
#define LEXICALTREE_H

#include "standard.h"
#include "phonemodel.h"
#include "languagemodel.h"
#include "shout_maketrainset.h"

#ifdef TEST_ARTICULATORY
  #include "articulatorystream.h"
#endif

struct LMLAType;

#define NUMBER_OF_BLAME_CLUSTERS  (6)

#define MAX_PDF_CALC_SIZE (1024000)

#ifdef CELL_PROCESSOR
typedef struct _control_block
{
  unsigned long long out_addr;  //beginning address of the output array
  unsigned long long in_addr[MAX_PDF_CALC_SIZE];
  unsigned char      in_length[MAX_PDF_CALC_SIZE];
  float              vector[39];
  int                nrPDFs;
  unsigned char      pad[24];
} control_block_getfastp_t;

typedef struct _control_block_lmla
{
  unsigned long long out_addr;  //beginning address of the output array
  unsigned long long in_addr;
  int                numberOfCompressedNodes;
  unsigned char      pad[8];
} control_block_createlmla_t;
#endif

////////////////////////////////////////////////////////////////////////////
/// \brief Used for N-Best/lattice calculations:
////////////////////////////////////////////////////////////////////////////

struct NBestList
{
  int        nrWords;
  int        noSilWords;
  float      totAM;
  float      totLM;  
};

////////////////////////////////////////////////////////////////////////////
/// \brief This type is used to keep track of analysis administration.
////////////////////////////////////////////////////////////////////////////

struct AnalysisSettings
{  
   WLRTracker          **prunePathAnalyses;
   int                   prunePathAnalysesLength;
   int                   beamCategories[3];
   bool                  containsOOV;
   int                   errorRegionActive;    // In this record so that we do not have to pass it in the recursive error-function.
   
   int                  *binDistr;
   float                *bestScore;
   float                *correctScore;
   float                *bestStateScore;
   int                  *stateRanking;
   
   float                 maxBeam;
   float                 maxStateBeam;
   float                 maxEndStateBeam;
   int                   maxHistogram;
   int                   maxStateHistogram;
   LexicalNode          *endStateActive;
   int                   currentLMHist[LM_NGRAM_DEPTH];
   
   NBestList             refStats;
};

////////////////////////////////////////////////////////////////////////////
/// \brief This type is used to make a global cache table for Lanuage Model
/// Look-Ahead (LMLA).
/// 
/// The 'used' parameter is added for making pruning faster. The wordHistory
/// defines the LM history for which this record is created. 'lookAhead' is 
/// an array (of length 'number of compressed noded') and next is used to
/// link the structures together.
////////////////////////////////////////////////////////////////////////////

struct LMLAGlobalListType
{
  int                 wordHistory[LM_NGRAM_DEPTH];
  int                 collissionTime;
  TokenType          *collissionNode;
  LMLAGlobalListType *next;
  float              *lookAhead;
};

////////////////////////////////////////////////////////////////////////////
/// \brief A lexical node is a single node of the lexical tree or Pronounciation
/// Prefix Tree (PPT). It contains token, AM and word information.
///
/// Because we use a static PPT (no copies are made), for each HMM state in
/// a node, it is possible that there are multiple valid tokens linked to it.
/// 
/// The modelID points to the PhoneModel object (AM model). If wordID is not 
/// -1, this node is a word-end node.
///
/// Each node can have a following node (next) in the tree. It may also have
/// a parallel node. Finaly it is possible that a node points to another 
/// lexical tree structure. This typically happens with alignment- or grammar-tasks.
////////////////////////////////////////////////////////////////////////////

struct LexicalNode
{
  TokenType     **tokenSeq;        ///< Normally of length NONSIL_NUMBER_OF_STATES (3)
  LexicalNode    *nextTree;        ///< If wordID != -1, this points to the start node of the next tree (for LVCSR back to the top)
  LexicalNode    *next;            ///< Next node. Multiple 'next' nodes are linked together (in the next node) with 'parallel'.
  LexicalNode    *parallel;        ///< Parallel node. 
  TokenType     **inputToken;
  int             modelID; 
  int             contextPrev;     ///< Previous modelID (left context)
  int             contextNext;     ///< Next modelID (right context)
  int             contextKey;      ///< For speed-up reasons pre-calculated: contextPrev * numberOfPhones + contextNext
  int             wordID;          ///< -1 if not an end node. Otherwise the ID of the word.
  int             depth;
  int             tokenSeqLength;
  int             compressedTreeIndex;
  bool            nodeIsActive;    ///< One or more states of this node are being used. 
  bool            toBeDeletedFromList;  ///< When this switch is TRUE. The node will be deleted from the 'active list'.   
};

////////////////////////////////////////////////////////////////////////////
/// \brief The LexicalNodeList is needed for administration of all currently
/// active nodes.
////////////////////////////////////////////////////////////////////////////

struct LexicalNodeList
{
  LexicalNodeList *next;
  LexicalNode     *node;
};


#ifdef SEARCH_STATISTICS_ON

////////////////////////////////////////////////////////////////////////////
/// \brief This structure keeps record of different kind of statistics during
/// decoding. It may be switched on and off by the compiler switch SEARCH_STATISTICS_ON
////////////////////////////////////////////////////////////////////////////

struct SearchStatistics
{
  int nrOfHistPruning;
  int nrOfCleanUps;
  int nrOfLMLATables_Temp;
  int nrOfLMLATables;
  int nrOfActiveNodes;
  int nrOfActiveTokens;
  int nrOfLMLATables_Prev;
  int nrOfActiveNodes_Prev;
  int nrOfActiveTokens_Prev;  
  int nrOfLMLATables_Max;
  int nrOfActiveNodes_Max;
  int nrOfActiveTokens_Max;  
};

#endif


struct FastCompressedTree
{
  int  wordID;
  int  nextIndex;
  bool hasPar;
};

////////////////////////////////////////////////////////////////////////////
/// \brief This class contains all lexical tree functionality and also almost
/// all decoder steps.
/// 
/// The lexical tree (or Pronunciation Prefix Tree (PPT)) is a special form
/// of a pronunciation vocabulary. All words that start with the same phones 
/// are put together. As soon as the pronunciation differs, the tree splits into
/// two sub-trees. (for example: "word" and "worst" will share the first two 
/// phones.)
/// 
/// This class handles the structure of the tree using LexicalNode, but it 
/// also handles tokens being propagated through the tree (including pruning), 
/// LM probabilities being added to tokens at tree leaves, LMLA, keeping track
/// of word- and phone-histories, picking the final recognition and storing
/// statistics.  
////////////////////////////////////////////////////////////////////////////

class LexicalTree
{
  protected:
   LMLAGlobalListType *lookAheadUnigram;
#ifdef TEST_ARTICULATORY
   ArticulatoryStream   *myArtStream;
#endif
   FastCompressedTree   *fastCompressedTree;
   int                   numberOfCompressedNodes;

   int                   nodeListLength;
   LexicalNode         **nodeArray;
   bool                  latticeGeneration;
   LexicalNode         **latticeWordList;
   int                   latticeWordListLength;
   int                  *tokenDepthAdmin;
   int                   biggestNodeDepth;
   int                   startTime;
   
   FILE                 *outputFile;             ///< Output hyp file.
   FILE                 *tdFile;                 ///< Token Distribution File.
   AnalysisSettings     *analysisSettings;       ///< Used to keep track of analysis administration (blame assignment).
   int                   blameAssingment[NUMBER_OF_BLAME_CLUSTERS];
   int                   analyse_deletions[NUMBER_OF_BLAME_CLUSTERS];
   int                   analyse_insertions[NUMBER_OF_BLAME_CLUSTERS];
   int                   analyse_substitutions[NUMBER_OF_BLAME_CLUSTERS];
   int                   analyse_totRefWords;
   int                   totCategories[3];   
   DecoderSettings       settings;
   WLRType              *bestRecPath;            ///< Used for adaptation, it would be inconvenient to find the best token each frame.

   int                   initialLMHist[LM_NGRAM_DEPTH]; ///< The history words at the start of a new utterence.
   int                   nrOfTokens;             ///< The number of tokens in the tree.
   int                   timeStamp;
   int                   numberOfPhones;         ///< The number of acoustic models types (not the number in the tree).
   int                   numberOfWords;          ///< The number of words in the DCT.
   LexicalNode         **treeEnd;                ///< The matrix of end-of-word nodes (so called phase-out).
   LexicalNode         **treeStart;              ///< The matrix of start-of-word nodes (phase-in)
   LexicalNode          *grammarStart;           ///< The start of the tree or (forced-align) grammar
   LexicalNode          *endNode;                ///< The last node. Only used for forced-align grammars.
   int                   grammarStartContext;    ///< Initial left context when starting the system (for forced-align).
   WLRType              *wlrStart;               ///< The first Word-Link-Record (linked together in WLRType).
   WLRTypeList          *wlrNBest;               ///< Each time frame the N-Best wlrs are marked as being one of the best!
   LatticeNode          *latticeAdmin;           ///< A list of all lattice nodes.
   int                   latticeN;
   int                   latticeL;
   char                **vocabulary;             ///< The vocabulary.
   PhoneModel          **phoneModels;            ///< An array of all distinct acoustic models.
   int                  *wordLength;
   TokenType            *bestToken;              ///< A pointer to the best token in the tree.
   float                 bestL;                  ///< The likelihood of token 'bestToken'.
   LanguageModel        *languageModel;          ///< The used language model.
   LexicalNodeList      *nodeList;               ///< A list of all active nodes for the following time frame (so that the tree does not have to be searched)
   bool                  alignParallel;          ///< Currently not used!
   LMLAGlobalListType   *lmla_list[LMLA_CACHE];  ///< The LMLA data structure. (see paper)
   float                *transitionPenalty;      ///< The transition penalty and short-word-penalty can be calculated once in the constructor.

   int                   intervalTimer;   

   float                *phoneLoopConfidence;
   int                   phoneLoopConfidenceOffset;
   LexicalNode         **borrowWordList;

   double *vList;
   int threadsRunning;
   MixGaussian **pdfUpdateList;

#ifdef SEARCH_STATISTICS_ON
   SearchStatistics      sentenceStats;
   SearchStatistics      globalStats;
#endif  
   int                   startOfSentenceWord;    ///< The index of the \<s\> word.
   int                   endOfSentenceWord;      ///< The index of the \</s\> word.
   
   bool                  currentlyAligning;      ///< If TRUE, the system will be aligning. Otherwise it does LVCSR.
   
  public:
                         LexicalTree         (FILE *outFile);
                         LexicalTree         (FILE *outFile, FILE *treeFile, bool useT = true);
                        ~LexicalTree         ();

  public: // for multi-threading
    void                 processNode         (LexicalNode *node, Vector *v);

  protected:

float oldcreateLMLAs(LMLAGlobalListType *lmlaGlobal, LexicalNode *node, float *allLMP);

    void setNodeContext  (LexicalNode *node,int leftContext);
    void createLatticeNodeGroups(WLRType *w, double lmScore, int wordID);

    void                 copyLMHistory(int *lmHistory1, int *lmHistory2);
    bool                 compareLMHistory(int *lmHistory1, int *lmHistory2);
    
    void                 printTokenDistribution();
    int                  createLatticeLMRescoring(LatticeNode *l, NBestList *scoreList, float amTot, float lmTot,
                                                  int nrWords, int noSilWords, int *lmHist, int sentenceID = 0);
    int                  countLatticePaths(LatticeNode *l, int sentenceID = 0);
    void                 sortLatticePaths(LatticeNode *l);
    void                 setWlrNBest(WLRType *wlr);  
    int                  getWordFromWLR(WLRType *wlr);
    LatticeNode         *findLatticeNodes(WLRType *wlr);
    void                 printLMParStats     (bool outputXML);
    void                 createLattice();
    void                 lattice_copyNonSilArcs(LatticeNode *source, LatticeNode *dest, int loopID);
    void                 lattice_removeDoubleArcs();
    void                 deleteLatticeAdmin();
    bool                 addWordStringToAlignment(LexicalNode *word);
    int                  getLastModelForContext(LexicalNode *node, int wordID);
    
    void                 calcErrorRegionStats(WLRType *wlr, WLRType *lastCorrect, int firstRef, int lastRef, 
                                              int *nrWordHyp, float *scoreLM);
      
    int                  printErrorString(int errorID);
    void                 processVector_LMLAReordering();  
    void                 processVector_LMLAReordering_prepare();  
    void                 processVector_processNodes(Vector **v);  
    void                 createActiveNodesList();

    void                 processVector_prune_processNodesOutput();
    void                 processVector_pruneLM();
    void                 processVector_grammar();
    void                 processVector_administrationCleanup();      

    void                 getBestGrammarEndToken(TokenType **bestT, bool notFinishedAllowed);
  
    void                 setNodeLocationPars (LexicalNode *node, bool fromParallel);
    void                 setDepthLevel       (LexicalNode *node, int phone, int depth);
    void                 readTree            (LexicalNode *node, FILE *treeFile, int length);
    void                 deleteLookAheadList ();
    void                 deleteTree          (LexicalNode *node, bool isPar);
    void                 initialiseNode      (LexicalNode *node);
    void                 pruneWithMinBeam    (LexicalNode *node, float minLikelihood_0);
    void                 deleteNodes         ();

    virtual void         processNodeOutput   (LexicalNode *node);
    void                 processWord         (int wordID, TokenType *token, char isSil, LexicalNode *resultNode);
    void                 getBestPath         (WLRType *wlr, bool outputXML, bool complete);
    void                 errorAnalysis       (WLRType *wlr, int depth, bool outputXML);
    void                 addNodeToList       (LexicalNode *node);

void pruneToken(TokenType **token, float minLikelihood, float binSize = 0.0, int *bins = NULL);

    void                 touchWLRpath        (WLRType *w);
    void                 touchWLRs           (TokenType *token);
    bool                 printWordPronunciation(LexicalNode *node, int wordID);
    void                 checkTreeRobustness (LexicalNode *node);
    void                 initialiseSystem    ();
    void                 setTreeStartEndMatrix ();
    TokenType           *findBestToken       (bool addEndOfSentence, bool complete, bool outputXML, bool notFinishedAllowed = false);
    LMLAGlobalListType  *getLMLATable        (int *lmHistory, bool onlyPrepare);
    int                  getLMLAHashKey      (int *lmHistory);
    void                 createLMLAs         (LMLAGlobalListType *lmlaGlobal, float *allLMP);
    void                 prepareLMLACreation (LexicalNode *node);
    void                 pruneLMLA           ();
    LexicalNode         *getPhoneString      (LexicalNode *node, int wordID);

    void                 getPhoneAlignment   (const char *prefix, PLRType *pt, int lastFrame, bool confident);
    LexicalNode         *findCorrectNode     (PLRType *pt, int wordID);

#ifdef SEARCH_STATISTICS_ON
    void                 updateStats         ();
#endif    
  public:
    bool                 alignmentIsEmpty    ();
    int                  latticeBaumWelch_numberNodes (LexicalNode *node, int number, bool clear);
    void                 latticeBaumWelch_setLikelihoods      (LexicalNode *node, Vector *t, int time, int numberOfStates, double *latticeLikelihood);
    void                 latticeBaumWelch_calculatePosteriors (double *latticeLikelihood, LexicalNode *node, double incomingScore, double *latticeBaumWelchAlfa,
                                                               double *latticeBaumWelchBeta, double *posteriors, int time, int numberOfStates);

    void                 latticeBaumWelch_mmi_accumulatorsPosteriors(LexicalNode *node, double *posteriors, int numberOfStates, Vector *observation);
    void                 latticeBaumWelch_printPosteriors (LexicalNode *node, double *posteriors, int time, int numberOfStates, int timeOffset);
    void                 latticeBaumWelch_forward     (double *latticeLikelihood, LexicalNode *node, double incomingScore,
                                                       double *latticeBaumWelchAlfa, int time, int numberOfStates);
    void                 latticeBaumWelch_initForward (double *latticeBaumWelchAlfa);
    double               latticeBaumWelch_backward    (double *latticeLikelihood, LexicalNode *node,
                                                       double *latticeBaumWelchBeta, int time, int numberOfStates, double normFactor, double *resArray);
    void                 latticeBaumWelch_initBackward(double *latticeBaumWelchBeta, int offset);
    void                 setLattice          (FILE *latFile);
    char                *getAlignmentString  ();
    void                 printInitialSettings(const char *amName, const char *dctName, const char *backName, const char *lmName, bool outXML);
    void                 printFinalSettings  (bool outXML, int totMilliSec, int totTime);
    void                 updateGlobalStats   ();
    void                 checkAMs            (int nrM, PhoneModel **models, bool outXML);
    void                 setAMs              (PhoneModel **models);
    void                 setLM               (LanguageModel *lm);
    int                  getWordID           (const char *word);
    char*                getWord             (int wordID);
    int                  getNumberOfWords    ();
    bool                 setForcedAlign      (const char *string);
    void                 setAlignParallel    ();
        
    virtual void         initialiseTree      (int startTime = 0);
    void                 setInitialLMHistory (const char *word1, const char *word2);

    void                 processVector       (Vector **v, int time);
    void                 adaptAMs            (Vector *v, int time);
    void                 getLogging          (const char *string);
#ifdef TEST_ARTICULATORY
    void                 testArticulatory    (ArticulatoryStream *s);
#endif
    int                  getBestIDSequence   (int *idList, int maxLength, bool showSil, bool notFinishedAllowed);
    void                 getBestRecognition  (bool addEndOfSentence, bool outputXML, bool complete, const char *label = NULL,
                                              int milliSec = 0, int totLength = 0, 
                                              const char *beginTime = NULL, const char *endTime = NULL);    
    double               getBestRecognitionScore();                                              
    void                 storePLConfidence   (int time);
    void                 setPhoneLoop        (int nrP, PhoneModel **models);
    void                 overwritePrunePars  (bool doHist, bool doBeam, bool doEndBeam,
                                              double beam, double state_beam, double endstate_beam,bool lmla,int histState, int hist);
    void                 overwriteWeightPars (double lmScale, double transPenalty, double silPenalty);
    bool                 safeBestRecognition (bool addEndOfSentence);
    LexicalNode         *borrowPhoneString   (const char *string);
    bool                 addForcedAlignOOV   (LexicalNode *oovNode);
    void                 createWordTree      ();
    
    void                 setPhoneLoopConfidence (float *phoneConf, int offset = 0);

    void                 printLattice        (FILE *latFile, const char *label, int timeEnd);
    void                 printNBestList      (FILE *nbestFile = NULL, LatticeNode *node = NULL);
    
    void                 setTokenDistributionFile (FILE *tdFile);
    
    void                 setLatticeGeneration(bool setting);
};

#endif
