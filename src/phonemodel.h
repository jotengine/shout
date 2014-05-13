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
#ifndef MODEL_H
#define MODEL_H

#include "standard.h"
#include "mixgaussian.h"
#include "vector.h"
#include "featurepool.h"

struct WLRType;
struct WLRTracker;
struct PLRType;
struct TokenType;
struct LatticeNode;
struct WLRList;
struct LMLAGlobalListType;
struct LexicalNode;

// So that we can make initialiseToken static
extern bool doPhoneAlignment;


typedef enum
{
  UNKNOWN = 0,
  OOV,
  SEARCH,
  AMLM,
  LM,
  AM
} BlameCluster; 

#define MAX_CLUSTERS  (10)

////////////////////////////////////////////////////////////////////////////
/// \brief This structure contains the name of an acoustic model and the 
/// number of samples that are used during training (per cluster).
////////////////////////////////////////////////////////////////////////////

struct DataStats
{
  char name[11];
  int nrOfSamples[MAX_CLUSTERS];
};

////////////////////////////////////////////////////////////////////////////
/// \brief The PLRType, the Phone Link Record Type is the structure that contains 
/// the phone history information for a single word stored in a WLR (of struct WLRType).
/// 
/// This structure is ONLY used when phone alignment is enabled!
///
/// It contains the phoneID that is beeing pronounced, the time at which the 
/// sound started (timeStamp) and the likelihood that the particular phone is 
/// pronounced given the history defined by the previous variable.
////////////////////////////////////////////////////////////////////////////

struct PLRType  // Phone Link Record
{
  int      phoneID;
  int      contextKey;
  unsigned char stateOffset[2];
  int      timeStamp;
  float    likelihood;  
  PLRType *previous;  
};

////////////////////////////////////////////////////////////////////////////
/// \brief The WLRType, the Word Link Record is the structure that contains 
/// word history information for tokens (of struct TokenType).
///
/// The wordID variable defines the word that is recognized. 
/// The timeStamp variable contains the time at which the word is started. 
/// COMBlikelihood contains the combined AM/LM likelihood.
/// (likelihoods are in the log scale). The variable usedAt is added for 
/// administrative purposes. This variable is updated each time frame as long as 
/// there is at least a single token history containing this WLR. If the usedAt
/// variable is not updated, the WLR is deleted. The variable adminNext is used
/// to make create a long list of WLRs. The previous variable points to the 
/// previous word pronounced. This variable is used to track the total word history.
///
/// The phoneAlignment variable points to the phone alignment history of this 
/// word. This variable is ONLY used when phone alignment is enabled!
/// 
////////////////////////////////////////////////////////////////////////////

struct WLRType  // Word Link Record
{
  int          lmHistory[LM_NGRAM_DEPTH];
  char         isSil;
  int          timeStamp;
  LatticeNode *lattice;
  WLRType    **nBest;
  float        COMBlikelihood;
  float        LMlikelihood;
  int          usedAt;
  WLRType     *adminNext;
  WLRType     *previous; 
  PLRType     *phoneAlignment;
};

////////////////////////////////////////////////////////////////////////////
/// The NBest list administration is done using this data type:
////////////////////////////////////////////////////////////////////////////

struct WLRTypeList
{
  WLRType    **nBest;   
  WLRTypeList *next;  
};

////////////////////////////////////////////////////////////////////////////
/// The lattice administration is done using this data type:
////////////////////////////////////////////////////////////////////////////

struct LatticeNode
{
  int            nodeNr;
  int            wordID;
  int            timeBegin;
  int            timeEnd;
  WLRType       *exampleWord;
  WLRList       *inArcs;   
  WLRList       *outArcs;
  LatticeNode   *adminNext;
};

/////////////////////////////////////////////////////////////////////////////////
/// The WLRList is a list of WLR's :-) (used for arc-administration in lattices)
/////////////////////////////////////////////////////////////////////////////////

struct WLRList
{
  WLRType *wlr;
  WLRList *next;
  float    amScore; // Also available in wlr, but this safes calculations.
  float    lmScore; // Also available in wlr, but this safes calculations.
  float    totScore; // Also available in wlr, but this safes calculations.
};


////////////////////////////////////////////////////////////////////////////
/// \brief The WLRTracker type, is created to track the single correct recognition
/// obtained during forced alignment. It can be used to check when the correct
/// path was left because of (for example) pruning.
///
/// The contextNext parameter can only be used when phone recognition is enabled.
////////////////////////////////////////////////////////////////////////////

struct WLRTracker
{
  WLRType  w;
  WLRType *linkWord;
  int wordID;
  int errorRegionID;
  BlameCluster errorCategory;
  int contextNext;
};

////////////////////////////////////////////////////////////////////////////
/// \brief The PhoneModel class defines the phone models, the LexicalTree class the tree structure and the TokenType
/// struct is the glue between them. LexicalTree is responsible for the token flow between phones and for 
/// language model lookahead, while PhoneModels will fill the likelihood variables and decises if a
/// token may be passed within the phone model.
///
/// The TokenType structure is the data type of each token. Because a token 
/// will always be in a list and it shall be possible to search forwards and backwards
/// in that list, the next and previous variables are added to the type.
///
/// The likelihood variable and lookAheadV(alue) contain the current likelihood of the token
/// and the part of this likelihood that is added due to language model lookahead.
/// lookAheadV is stored in order to quickly being able to substract the lookahead from
/// the real likelihood and replace it with a new value.
///
/// The lookahead value is calculated with help of an entire LM lookahead tree. This tree contains all
/// possible lookahead values and is shared between all tokens with the same LM history. The variable 
/// lmLookAhead contains a pointer to this lookahead tree.
///
/// In order to make backtracking word history possible a Word Link Record list is stored in the path variable.
/// If phone alignment is enabled, the phonePath variable contains the phone history information.
///
////////////////////////////////////////////////////////////////////////////

struct TokenType
{
  TokenType           *next;
  float                likelihood;
  WLRType             *path;
  PLRType             *phonePath;
  float                lookAheadV;
  LMLAGlobalListType  *lmLookAhead;
};

////////////////////////////////////////////////////////////////////////////
/// \brief This structure contains all data for one single gaussian mixture
/// state set. This includes the PDF and state transition probabilities.
////////////////////////////////////////////////////////////////////////////

struct MixtureSet
{
  MixGaussian *state;
  double       transitionP_toSelf;
  double       transitionP_toNext;
  double       currentVectorP;
};

////////////////////////////////////////////////////////////////////////////
/// \brief This structure contains a summary of the data in an acoustic 
/// model stored as a PhoneModel object.
/// 
/// The parameters 'frameMeanLikelihood' and 'unused' are added for future use.
/// Please ignore them (especialy frameMeanLikelihood, which is filled falsely).
////////////////////////////////////////////////////////////////////////////

/// \todo Change this data type! (there is an 'unUsed' in read/write and isSil is double).
struct ModelStats
{
  char   name[10];
  int    nrOfGaussians;
  int    nrOfContexts;
  int    maxNrOfContexts;
  int    nrOfTrainOcc;
  double likelihood;
  double frameMeanLikelihood;
  double isSil;              // Now a double, because there was already a double variable here "unused"...
};

////////////////////////////////////////////////////////////////////////////
/// \brief The PhoneModel class handles likelihood calculation of phones given an observation sequence.
///
/// For this implementation, the phone models are represented by Hidden Markov Models (HMM).  
/// Training the HMM is done by the TrainPhoneModel class. TrainPhoneModel stores the HMM paramaters
/// in a binary file (all models together form the binary acoustic model file) and PhoneModel
/// can load the model in memory at startup.
///
/// PhoneModel objects are able to determine the likelihood that a phone is pronounced,
/// given an observation sequence AND a TokenType string. PhoneModel does not store
/// state tokens itself, it only handles the (static) parameters needed to calculate HMM likelihoods.
/// Because the user handles the token administration (in LexicalTree), it is possible
/// to use a phone model in more than one node, without copying its parameters.
///
/// The acoustic models are context-dependent models. During training (TrainPhoneModel) it is decided
/// which contexts share one or more states and transition probabilities. Each model contains a pool
/// of states, stored in the variable mixtureSetData. The arrays stateMix_1, stateMix_2 and stateMix_3 
/// contain for each context the index of mixtureSetData. stateMix_1 is used to determine which state
/// of the pool-of-states should be used for the first state of the HMM. stateMix_2 is used for the 
/// second and stateMix_3 for the third. The index for stateMix_x is calculated as followes: 
///
/// contextKey = leftContext * numberOfPhones + rightContext
///
/// Therefore, the first state of the context 'A' - 'l' - 's' (where this model is 'l' and the left/right
/// context is 'A'/'s' and 'A' is the 2nd phone and 's' the 8th) can be found by:
///
/// First HMM-state = mixtureSetData[stateMix_1[2*54+8]]
///
////////////////////////////////////////////////////////////////////////////

class PhoneModel
{
  public:                         /// \todo This is a quick fix to check if this can help SAT training. Will change it back to protected later!
    MixtureSet  *mixtureSetData;  ///< The entire set of states and transition parameters used by this model.

  protected:
    // Attributes:
    ModelStats   statistics;      ///< Statistical information about this acoustic model.
    int         *timeStamp;       ///< The model only calculates some context info once for every 'timestamp' frame. 
    int         *stateMix_1;      ///< An array with the index i for mixtureSetData[i] for every context (state 1).
    int         *stateMix_2;      ///< An array with the index i for mixtureSetData[i] for every context (state 2).
    int         *stateMix_3;      ///< An array with the index i for mixtureSetData[i] for every context (state 3).
    bool         isSil;
    int          dimensions;
    int          silRinglastPos;
    float       *weightRinglastPos;

  public:
    // Constructor and destructor:
                PhoneModel       (int dim = ASR_DEFAULT_VECTORSIZE);
                PhoneModel       (FILE *inFile, int dim = ASR_DEFAULT_VECTORSIZE, bool onlyUseFastP = false);
                PhoneModel       (MixGaussian *mix, double toNext);
               ~PhoneModel       ();
 
    // Methods:    
  protected:
    static bool replaceTokenLM   (TokenType *nt,TokenType *token, float  like, float  *bestL);
    static void addChain_lmla    (TokenType **tokenNew, float  likelihood, TokenType *token, int index, DecoderSettings *settings, float  *bestL);    

  public:
    
    // INLINE FUNCTIONS FOR CODE OPTIMALIZATION:
    inline int dim               ()  const {return dimensions;}
    inline bool isSilModel       ()  const {return isSil;}

    static void addChain         (TokenType **tokenNew, float  likelihood, TokenType *token,
                                  int stateNr, int curTime, DecoderSettings *settings, float  *bestL, bool checkCollission);    
    static void addChain_ordered(TokenType **tokenNew, float  likelihood, TokenType *token,
                                  int stateNr, int curTime, DecoderSettings *settings, float  *bestL);    


    static void initialiseToken(TokenType **token);

    void        setSilTrans      (double tr);
    

    int         getStateNr       (int contextKey, int state);
        
    ModelStats *getStatistics    (void);
    int         getNumberOfGaussians();
    
    virtual int         touchPDF         (int contextKey, int t, MixGaussian **updateThese, double **resultHere);
    virtual void        processVector    (int contextKey, Vector *v, int t, int index,
                                  TokenType **token, int tLength, TokenType  *inToken, DecoderSettings *settings, float  *bestL);
    virtual void        getOutput        (int contextKey, TokenType **token, int tLength, TokenType **outToken, DecoderSettings *settings, float  *bestL);

    virtual double      getLookaheadLogP  (double *vectorList, int timeStamp, bool doSecondHalf);


    void        copyGaussians    (MixGaussian *destMixGaussian, int maxNmbr);

    void        resetPhoneAdmin  ();

    void        mapAdaptMeans           ();
    void        adapt_setInitialNode    (Adapt_AM_TreeNode *node);
    void        adapt_setNode           ();
    void        adapt_addAcumulatorData (int state, int contextKey, Vector *observation, double probability = 1.0);
    void        adapt_setHelperMatrices ();
    void        adapt_clear             ();
    void        adapt_adapt             ();
    void        adapt_setVarTrans       ();
    void        adapt_adaptVar          ();
    void        adapt_unAdapt           ();
    void        adapt_setAcumulators    (int useLabel, int useSegmentation, FeaturePool *usePool);
    void        writeAccumulators       (FILE *file, FILE *fileF = NULL, FILE *fileST = NULL, bool doBinary = false);
    void        addAccumulators         (FILE *file);

    void        writeModel              (FILE *outFile);

    double      getLogPDFProbability           (int contextKey, Vector *v);
    double getPDFProbability            (int contextKey, Vector *v, int stateNr, int time);
    double      getTransition                  (int contextKey, bool toSelf, int stateNr);


    void        getSilSumHist           (Vector **histogram);

    static PLRType *copyPhonePath       (PLRType *pP);
    static void     initialisePhonePath (PLRType *t);

    void printModel(FILE *fileMean, FILE *fileVariance, FILE *fileWeight);
    int  printInfo(Vector *v);
};

#endif
