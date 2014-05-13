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
#ifndef FEATUREPOOL_H
#define FEATUREPOOL_H

#include "vector.h"
#include "featureextraction.h"

#define MAX_NR_SEGMENTATIONS   (25)
#define MAX_NR_SPKIDS          (10000)

#define HISTNORM_NR_GAUSSIANS   (2)

#define CTS_MINSPEECH          (30)
#define CTS_MINSILENCE         (30)
#define CTS_SILCOLOR            (5)
#define CTS_SILFACTOR         (2.0)

#define CTS_SHORT_MINSPEECH    (5)
#define CTS_SHORT_MINSILENCE   (5)
#define CTS_SHORT_SILCOLOR      (0)
#define CTS_SHORT_SILFACTOR   (1.0)

#define ENERGY_COMP             (0)
#define ZEROCROSS_COMP          (1)


#define FEATUREPOOL_MAX_NR_STREAMS (3)

/////////////////////////////////////////////////////////////////////////////////////////
/// \brief This structure contains all feature file formats that are supported by shout.
/////////////////////////////////////////////////////////////////////////////////////////

class FeaturePool;


// Currently, only one file type is supported:
typedef enum
{
  FILE_TYPE_SHOUT              = 0,
  FILE_TYPE_RAW_DECODING       = 1,
  FILE_TYPE_RAW_HISTOGRAM_NORM = 2,
  FILE_TYPE_RAW_DIARIZATION    = 3,
  FILE_TYPE_RAW_VTLN           = 4,
  FILE_TYPE_RAW_SAD            = 5,
  FILE_TYPE_RAW_SPKREC         = 6,
  FILE_TYPE_TEXT               = 7
} FILE_TYPE;


/////////////////////////////////////////////////////////////////////////////////////////
/// \brief Linked list that defines a segmentation of the feature vector pool.
///
/// The ID variable can be used for different purposes. During training, the IDs can be
/// used to store the left- and right contexts of the phones. If needed, pointers to
/// bigger data structures can be stored in ID. ID[0] is used as 'cluster ID'.
/////////////////////////////////////////////////////////////////////////////////////////

struct SegmentationList
{
  int               firstFrame;
  int               lastFrame;
  int               ID[4];
  char             *infoString;
  SegmentationList *next;
};

struct SegmentationAdmin
{
  SegmentationList *curSeg;
  SegmentationList *prevSeg;
  int               timeStamp;
};

/////////////////////////////////////////////////////////////////////////////////////////
/// \brief Linked list that stores a list of unique IDs.
///
/// This structure is only used by FeaturePool::writeRTTM() in order to write the
/// SPKR-INFO lines.
/////////////////////////////////////////////////////////////////////////////////////////

struct UniqueIDList
{
  int           ID;
  UniqueIDList *next;
};

struct FeaturePoolInfo
{
  int     numberOfChannels;
  int     lastComponent[FEATUREPOOL_MAX_NR_STREAMS];
  double  weight[FEATUREPOOL_MAX_NR_STREAMS];
  double  weightBackup[FEATUREPOOL_MAX_NR_STREAMS];
  Gaussian *distrG[FEATUREPOOL_MAX_NR_STREAMS];
};

struct FeaturePool_Multistream_file
{
  char        *fileName;
  FILE_TYPE    fileType;
  FeaturePool *featurePool;
};

struct FeatureEnergy
{
  FeatureExtraction  *feaExtract;
  Vector            **buffer;
};

////////////////////////////////////////////////////////////////////////////
/// \brief FeaturePool objects contain a feature vector file in memory.
///
/// It is also possible to read out one vector at a time so that it is not
/// needed to keep the file in memory (for decoding large files). If you want
/// this, set 'noBuffering' in the constructor to true.
///
/// For phone training or adaptation, it is often not known how big
/// the entire pool with phone data is going to be. For this case it is
/// possible to initialise the pool with a maximum length and
/// add phone data while making the segmentation.
////////////////////////////////////////////////////////////////////////////

class FeaturePool
{
  protected:
    // Variables:
    bool                onlyMEL;
    FILE               *featureFile;
    FeatureExtraction  *feaExtract;
    FeatureEnergy       feaEnergy[2];  // Left and right channel...
    FILE_TYPE           fileType;

    Vector            **buffer;
    int                 poolLength;

    int                 fileOffset;
    int                 poolFilled;

    bool               *filter;
    bool                filterEnabled;

    int                 restrictStart;
    int                 restrictEnd;

    int                 vectorSize;

    int                 numberOfSpkIDs[MAX_NR_SEGMENTATIONS];
    int                 spkID[MAX_NR_SEGMENTATIONS][MAX_NR_SPKIDS];
    bool                freeID[MAX_NR_SEGMENTATIONS];
    FEATURENORM_TYPE    featureNorm;
    SegmentationList   *segList[MAX_NR_SEGMENTATIONS];
    SegmentationList   *fillPool_curSeg[MAX_NR_SEGMENTATIONS];
    FeaturePoolInfo     infoBlock;
    Gaussian           *textFeatureDistribution;
    FeaturePool_Multistream_file *multiStream;

  public:
    // Constructors:
            FeaturePool                  (int nrStreams, FeaturePool_Multistream_file *files, double *iw, bool onlyPrepare);
            FeaturePool                  (char *featureFileName, FILE_TYPE fileType = FILE_TYPE_SHOUT, 
                                          bool noBuffering = false, double vtln = 1.0, bool doOffset = false, bool onlyPrepare = false);
            FeaturePool                  (int poolS, int vectorS);
           ~FeaturePool                  ();

    void    setStreamWeights             (double *w, bool forEver);
    void    resetStreamWeights           ();

    // Public methods:
    void    storeFeatureVectors          (FILE *file);
    void    setOnlyMEL                   ();
    void    createNewMultiStreamPool     (FeaturePool_Multistream_file *files, int segID);
    int     mergeSegmentations           (int *seg, int nrSeg, int outSeg);
    void    retrieveEnergy               (char *featureFileName, int channelNr);
    void    removeLowEnergy              (int segID, int channelNr, bool doShort = false);
    double  getAverageEnergy             (int segID, int channelNr, double varianceFactor);
    double  getCurSegmentEnergy          (SegmentationAdmin *admin,int t);
    void    splitOnEnergy                (int segNr, int sourceLabel, int silLabel, int soundLabel, 
                                          int nrSamplesSil, int nrSamplesSound, int comp);
    
    void    createNewPool                (char *featureFileName, char *histNormName, char *histNormSpkName,
                                          bool noBuffering = false, bool onlyPrepare = false, int vtlnSegID = -1, double globalVtln = -1.0);
    void    addFeatureFile               (char *featureFileName, FILE_TYPE fileType, int ignoreFrames, int trainFrames, double vtln = 1.0);
    void    setFilter                    (char *featureFileName, int dimension, double threshold);
    void    enableFilter                 (bool enable);
    void    filterToSegmentation         (int segNr, FeaturePool *feaPool = NULL);
    void    setFeatureNormalization      (FEATURENORM_TYPE m);
    int     getNumberOfSpeakers          (int segNr);
    int     getNumberOfSegments          (int segNr, int spkrNr);
    int     getSpeakerID                 (int segNr, int spkrNr);
    
    
    bool    createUEMSegmentation        (int segID, int channelNumber, char *label, char *UEM);
    
    int     claimSegmentationID          ();
    void    freeSegmentationID           (int segNr);
    
    void    resetSegmentation            (int segNr);
    void    addSegment                   (int segNr, int firstFrame, int lastFrame, int ID, int cKey = 0, int cLeft = 0, int cRight = 0);
    void    addSegmentToFront            (int segNr, int firstFrame, int lastFrame, int ID, int cKey = 0, int cLeft = 0, int cRight = 0);
    void    setCurSegmentVectors         (int segNr, FILE *inFile, bool inputIsFloat = false);
    void    setCurSegmentVector          (int segNr, int offset, Vector *vector);
    int     getPoolLength                ();
    int     getClusterLength             (int segNr,int label);
    bool    doSegmentsExist              (int segNr, int ID);


    Vector **getVectorList               (SegmentationAdmin *admin, int time);
    Vector **getFirstVectorFirstSegmentList(SegmentationAdmin *admin, int segNr, int ID, bool dontCareID = false, 
                                          bool *isLastFromSegment = NULL, bool *inFilter = NULL);
    Vector **getFirstVectorNextSegmentList(SegmentationAdmin *admin, bool dontCareID = false, bool *isLastFromSegment = NULL, 
                                          bool *inFilter = NULL);
    Vector **getCurSegmentVectorList      (SegmentationAdmin *admin, int offset);
    Vector **getNextVectorList            (SegmentationAdmin *admin, bool *isLastFromSegment = NULL, bool *inFilter = NULL);

    
    void    segmentationCreateOverlap    (int segNr, int label);
    void    segmentationAddLabel         (int sourceSegNr, int destSegNr, int sourceLabel, int destLabel);
    void    segmentationSubset           (int sourceSegNr, int destSegNr, int label, PhoneModel *model, int nrSamples);
    void    segmentationCopy             (int sourceSegNr, int destSegNr);
    void    segmentationIntersection     (int sourceSegNr, int filterSegNr, int filterLabel, int destSegNr);
    void    segmentationFillGaps         (int segNr, int fillLabel);
    void    segmentationSubstract        (int sourceSegNr, int filterSegNr, int destSegNr, int filterLabel = -1);
    void    segmentationSplit            (int sourceSegNr, int destSegNr1, int destSegNr2);
    void    segmentationWiden            (int segNr, int soundNr, int maxNrFrames, int onlyWidenLabel = -1);
    int     filterLongSegments           (int segID, int labelNr, int newLabel, int minLength);
    void    filterShortSegments          (int segID, int labelNr, int labelSelf, int minNeighbours);
    
    void    maxSegmentLength             (int segNr, int labelNr, int size, int minRemainder);
    
    void    concatSegmentation           (int destSeg, int copySeg);
    
    void    restrictPool                 (int start, int end);
    
    void    mergeLabels                  (int segNr, int label1, int label2);

    void    readRTTM                     (int segNr, char *mapCluster, int mapNr, int strLen, int prefixLength, FILE *inFile, int fillerID=-1,
                                          int ignoreFrames = 0, int trainFrames = 0, char *labelName = NULL);
    void    writeRTTM                    (int segNr, const char *label, const char *prefix, FILE *outFile = NULL, int nrPr = 0, int overwriteVTLN = -1);

    // INLINE FUNCTIONS FOR CODE OPTIMALIZATION:

inline Vector *getFirstVectorFirstSegment(SegmentationAdmin *admin, int segNr, int ID, bool dontCareID = false, 
                                          bool *isLastFromSegment = NULL, bool *inFilter = NULL)
{
 Vector **res = getFirstVectorFirstSegmentList(admin,segNr, ID, dontCareID, isLastFromSegment, inFilter);
 if(res != NULL)
 {
   return res[0];
 }
 return NULL;
};

inline Vector *getVector(SegmentationAdmin *admin, int time)
{
 Vector **res = getVectorList(admin,time);
 if(res != NULL)
 {
   return res[0];
 }
 return NULL;
};

inline Vector *getFirstVectorNextSegment(SegmentationAdmin *admin, bool dontCareID = false, bool *isLastFromSegment = NULL, bool *inFilter = NULL)
{
  Vector **res = getFirstVectorNextSegmentList(admin,dontCareID,isLastFromSegment,inFilter);
  if(res != NULL)
  {
    return res[0];
  }
  return NULL;
};

inline Vector *getCurSegmentVector(SegmentationAdmin *admin, int offset)
{
  Vector **res = getCurSegmentVectorList(admin,offset);
  if(res != NULL)
  {
    return res[0];
  }
  return NULL;
};
inline Vector *getNextVector(SegmentationAdmin *admin, bool *isLastFromSegment = NULL, bool *inFilter = NULL)
{
  Vector **res = getNextVectorList(admin,isLastFromSegment,inFilter);
  if(res != NULL)
  {
    return res[0];
  }
  return NULL;
};


    inline double getStreamWeight(int i)    const {return infoBlock.weight[i];};
    inline int getVectorSize    ()          const {return vectorSize; };
    inline int getSegmentID     (int idNr,SegmentationAdmin *admin)  const { 
                                                    assert(idNr >= 0 && idNr < 4 && admin->curSeg != NULL);
                                                    return admin->curSeg->ID[idNr]; };
    inline void setSegmentID    (int idNr, int val,SegmentationAdmin *admin)  const {
                                                              assert(idNr >= 0 && idNr < 4 && admin->curSeg != NULL);
                                                              admin->curSeg->ID[idNr] = val; };
    inline int startTimeStamp   ()          const { return restrictStart; };
    inline int getCurSegmentLen (SegmentationAdmin *admin)          const {
                                                   assert(admin->curSeg != NULL);
                                                   int start = admin->curSeg->firstFrame;
                                                   int end   = admin->curSeg->lastFrame;
                                                   if(start < restrictStart) { start = restrictStart; }
                                                   if(end > restrictEnd) { end = restrictEnd; }
                                                   int length = end-start+1;
                                                   if(length < 0){length = 0;}
                                                   return length; };
    inline int getCurSegmentStart (SegmentationAdmin *admin)        const {
                                                    assert(admin->curSeg != NULL); 
                                                   if(admin->curSeg->firstFrame < restrictStart){ return restrictStart; }
                                                   return admin->curSeg->firstFrame; };
    inline char *getCurInfoString (SegmentationAdmin *admin)        const {
                                                    assert(admin->curSeg != NULL); 
                                                   return admin->curSeg->infoString; };
    inline FeaturePoolInfo *getStreamInfo()       {return &infoBlock;};
  protected:
    // Protected methods:
    void    createNewFeatureExtraction   (char *featureFileName, double vtln, bool doOffset, bool onlyPrepare);
    int     readHeader_ShoutFormat       ();
    void    fillBuffer_ShoutFormat       (int length, int ignoreFrames = 0);
    int     readHeader_TextFormat        ();
    void    fillBuffer_TextFormat        (int length, int ignoreFrames = 0);
    void    fillBuffer_fromFeatureExtraction(int length, int ignoreFrames = 0);
};

#endif
