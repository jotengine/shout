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

#include <stdlib.h>
#include <string.h>
#include "standard.h"
#include "shoutconfig.h"
#include "shout-misc.h"
#include "vector.h"

using namespace StringFunctions;


bool doPhoneAlignment = false;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

ShoutConfig::ShoutConfig(APPLICATION_ID appID,int argc,char *argv[])
{
  parShort      = new char*[ARGUMENTS_NUMBER];
  parLong       = new char*[ARGUMENTS_NUMBER];
  applicationID = appID;
  parameter     = new CommandParameterType[ARGUMENTS_NUMBER];

/////////////////////////////////////// global
  initializeParameter(ARGUMENT_HELP,           "h", "help",          "Shows this message. When 'ALL' is specified, all toolkit parameters are shown", "'ALL'", "");
  initializeParameter(ARGUMENT_CONFIG,         "c", "config",        "All parameters may be defined in this config file using the long parameter names", "filename", "");
  initializeParameter(ARGUMENT_CREATE_CONFIG, "cc", "create-config", "Create a config file for this application, or a global file when 'ALL' is passed",   "'ALL'", "");
/////////////////////////////////////// shout_segment
  initializeParameter(ARGUMENT_AUDIOFILE,      "a", "audio",         "The audio file that will be processed (must be raw, 16KhZ, 16 bit, little endian)",
                                                                     "filename", "");
  initializeParameter(ARGUMENT_ENERGY_BASED,  "eb", "energy",        "Segmentation will be only based on mean energy.","define","");
  initializeParameter(ARGUMENT_AM_SEGMENT,     "ams", "am-segment", "Acoustic model for SAD", "filename", "");
  initializeParameter(ARGUMENT_META_OUT    ,   "mo", "meta-out", "The resulting meta data file (see manual for the format)", "filename", "");
  initializeParameter(ARGUMENT_LABEL,          "l", "label", "Label (ID) of audio file", "label", "SpeechNonSpeech");
  initializeParameter(ARGUMENT_UEM,            "u", "uem", "UEM file to restrict regions of interest", "filename", "");
  initializeParameter(ARGUMENT_CTS_IN,         "ctsi", "cts-in", "Used as second channel input (CTS)", "second input file name", "");
  initializeParameter(ARGUMENT_CTS_OUT,        "ctso", "cts-out", "Used as second output info file (CTS)", "second output file name", "");
  initializeParameter(ARGUMENT_WIDEN_SEGMENTS, "w", "widen", "After processing, the final segments will be widened", "define", "");
/////////////////////////////////////// shout_cluster
//   initializeParameter(ARGUMENT_AUDIOFILE
//   initializeParameter(ARGUMENT_LABEL
//   initializeParameter(ARGUMENT_META_OUT
//   initializeParameter(ARGUMENT_WIDEN_SEGMENTS
  initializeParameter(ARGUMENT_META_IN,        "mi", "meta-in", "The input meta data file (see manual for the format)", "filename", "");
  initializeParameter(ARGUMENT_SECOND_STREAM,  "s", "second-stream", "A second stream (floating-point vectors) will be used", "filename", "");
  initializeParameter(ARGUMENT_THIRD_STREAM,  "t", "third-stream", "A third stream (posteriors) will be used", "filename", "");
  initializeParameter(ARGUMENT_FAST_MERGE,     "fm", "fast-merge", "For large files: faster merging (maximum number of merges per run)", "number", DIARIZATION_FAST_MERGE);
  initializeParameter(ARGUMENT_MAX_CLUSTERS,   "mc", "max-clusters", "Maximum number of clusters at initialization (not defined=unlimited)", "number", DIARIZATION_MAX_CLUSTERS);



/////////////////////////////////////// shout_vtln
//   initializeParameter(ARGUMENT_AUDIOFILE
//   initializeParameter(ARGUMENT_META_IN
//   initializeParameter(ARGUMENT_META_OUT
  initializeParameter(ARGUMENT_PERLINE,        "fnrm", "feat-norm",   "Method used for normalizing features (VTLN,CMN,CVN)", "SEGMENT|ALL|SPEAKER", "SPEAKER");
  initializeParameter(ARGUMENT_AM_VTLN,        "amv", "am-vtln", "The VTLN acoustic model", "filename", "");
  initializeParameter(ARGUMENT_AM_HISTNORM,    "amh", "am-hist", "The global histogram feature normalization model", "filename", "");
//  initializeParameter(ARGUMENT_AM_HISTNORMSPK, "amhs", "am-hist-speaker", "The speaker histogram feature normalization model", "filename", "");

/////////////////////////////////////// Shout:
//   initializeParameter(ARGUMENT_AUDIOFILE
//   initializeParameter(ARGUMENT_META_IN
//   initializeParameter(ARGUMENT_META_OUT
  initializeParameter(ARGUMENT_HYPFILE,        "hyp", "hyp", "The hypothesis file (speech recognition result)", "filename", "");
  initializeParameter(ARGUMENT_AM_PHONE,       "amp", "am-phone", "This file must be in shout format", "acoustic models file", "");
  initializeParameter(ARGUMENT_LM_BIN,         "lm", "lm", "This file must be in shout format", "language model file", "");
  initializeParameter(ARGUMENT_DCT_BIN,        "dct", "dct", "This file must be in shout format", "dictionary file", "");
  initializeParameter(ARGUMENT_PHONEBASED,     "pb", "phone-based", "The recognition will be aligned on phone level (when defined)", "define", "");
  initializeParameter(ARGUMENT_LATTICE,        "lat", "lattice", "will write lattice files to the output directory in PSFG format. If lattice-rescore is defined, it will use this directory for reading. For shout_online, this lattice will be used for decoding (provide full path, not dir)", "input/output dir", "");
  initializeParameter(ARGUMENT_NBEST,          "nb", "nbest", "will write N-best files to the output directory in plain text format", "output dir", "");
  initializeParameter(ARGUMENT_ALIGN,          "al", "align", "will use the transcription in the DBL file for forced alignment", "define", "");
  initializeParameter(ARGUMENT_LATTICERESCORE, "lr", "lattice-rescore", "Enable lattice rescoring (read from lattice directory)", "define", "");

  initializeParameter(ARGUMENT_ANALYSE,        "ana", "analyse", "will use the transcription in the DBL file to analyse the recognition", "define", "");
  initializeParameter(ARGUMENT_TOKENDIST,      "td", "token-distance", "will output token/active-node distribution statistics to file", "output dir", "");
  initializeParameter(ARGUMENT_PHONECONF,      "pc", "phone-confidence", "will calculate the maximum phone-likelihood at each timeframe", "define", "");
  initializeParameter(ARGUMENT_BACKDCT_BIN,    "bdct", "background-dct", "used to make alignment possible despite OOV words (must be in shout format).", "dictionary file", "");
  initializeParameter(ARGUMENT_TRAINDIR,       "tr", "traindir", "recognised or aligned phones will be placed in this directory for (re-)training", "traindir", "");
  initializeParameter(ARGUMENT_TRAINDIR2,      "tr2", "traindir2", "Training directory that needs to be added to 'traindir'", "traindir", "");
  initializeParameter(ARGUMENT_XML,            "xml", "xml", "output will be in XML format. The DTD is stored in the package data directory", "define", "");
  initializeParameter(ARGUMENT_LM_SCALE,       "lms", "lm-scale", "Language model scaling factor", "floating point number", LM_SCALE);
  initializeParameter(ARGUMENT_TRANS_PENALTY,  "trp", "transition-penalty", "Transition penalty", "floating point number",  TRANS_PENALTY);
  initializeParameter(ARGUMENT_SIL_PENALTY,    "silp", "silence-penalty", "Silence penalty", "floating point number",       SIL_PENALTY);
  initializeParameter(ARGUMENT_BEAM,           "b", "beam", "Maximum global beam", "floating point number",                 BEAM);
  initializeParameter(ARGUMENT_STATE_BEAM,     "stb", "state-beam", "Maximum beam per state", "floating point number",      STATE_BEAM);
  initializeParameter(ARGUMENT_END_STATE_BEAM, "endb", "end-state-beam", "Maximum beam at the end state", "floating point number", END_STATE_BEAM);
  initializeParameter(ARGUMENT_LMLA,           "lmla", "lmla", "Determines if language model lookahead will be used", "define", LMLA);
  initializeParameter(ARGUMENT_HISTOGRAM_STATE_PRUNING,  "hsp", "histogram-state-prune", "Maximum number of tokens allowed in each state", "floating point number", HISTOGRAM_STATE_PRUNING);
  initializeParameter(ARGUMENT_HISTOGRAM_PRUNING,        "hp", "histogram-prune", "Maximum number of tokens allowed in the network", "floating point number", HISTOGRAM_PRUNING);

  initializeParameter(ARGUMENT_CONFIGDIR, "cd", "configdir", "Directory for configuring the current speaker", "directory name", "");


/////////////////////////////////////// shout_train_mmi
  initializeParameter(ARGUMENT_LATTICETRAINING, "ltr", "lattice-training", "file with (enum or denom) lattice accumulators for MMI training", "file name", "");
  initializeParameter(ARGUMENT_LATTICETRAINING_ENUM, "lte", "lattice-enum", "file with the enumurator lattice accumulators for MMI training", "file name", "");
  initializeParameter(ARGUMENT_LATTICETRAINING_DENOM, "ltd", "lattice-denom", "file with the denomurator lattice accumulators for MMI training", "file name", "");



/////////////////////////////////////// shout_maketrainset
//   initializeParameter(ARGUMENT_AUDIOFILE
//   initializeParameter(ARGUMENT_HYPFILE
//   initializeParameter(ARGUMENT_TRAINDIR
//   initializeParameter(ARGUMENT_META_IN
//   initializeParameter(ARGUMENT_META_OUT
  initializeParameter(ARGUMENT_PHONE_LIST,     "pl", "phone-list", "List of phones in Shout format", "filename", "");
  initializeParameter(ARGUMENT_AM_TYPE,        "amt", "am-type", "Type of model to create", "SAD|VTLN|DIAR|HIST|PHONE", "");
  initializeParameter(ARGUMENT_AUDIO_LIST,     "trl", "train-list", "A list of audio-, meta-, hyp-files.", "filename", "");
  initializeParameter(ARGUMENT_SPEAKER_FILTER, "spk", "speaker", "Only use speech from this speaker", "label", "");



/////////////////////////////////////// shout_adapt_am
//   initializeParameter(ARGUMENT_TRAINDIR
//   initializeParameter(ARGUMENT_AM_PHONE
  initializeParameter(ARGUMENT_AM_CLUSTER,     "amc", "am-cluster", "Acoustic model used for clustering", "filename", "");
  initializeParameter(ARGUMENT_AM_OUT,         "amo", "am-out", "Adapted acoustic model", "filename", "");

/////////////////////////////////////// shout_dct2lextree
//   initializeParameter(ARGUMENT_DCT_BIN
//   initializeParameter(ARGUMENT_PHONE_LIST
  initializeParameter(ARGUMENT_DCT_TEXT,       "d", "dct-text", "Dictionary in text format (word phone-list)", "filename", "");

/////////////////////////////////////// shout_lm2bin
//   initializeParameter(ARGUMENT_DCT_BIN
//   initializeParameter(ARGUMENT_LM_BIN
  initializeParameter(ARGUMENT_LM_ARPA,        "arpa", "lm-arpa", "ARPA language model", "filename", "");

/////////////////////////////////////// shout_train_model
//   initializeParameter(ARGUMENT_TRAINDIR

  initializeParameter(ARGUMENT_DOFAST,         "fst", "fast", "Do fast, less accurate, training", "define", "");
  initializeParameter(ARGUMENT_DECISIONTREE,   "dtr", "decision-tree", "Decision tree in Shout (text) format", "filename", "");
  initializeParameter(ARGUMENT_MAXGAUSSIANS,   "mg", "max-gaussians", "Maximum number of gaussians per model", "number", 32.0);
  initializeParameter(ARGUMENT_TRAINSILP,      "tsp", "train-sil-p",   "Probability that a SIL training sample will be used for training", "number", 1.0);
  initializeParameter(ARGUMENT_TRAINSILMAX,    "tsm", "train-sil-max", "Maximum number of SIL training samples that will be used for training", "number", -1.0);
  initializeParameter(ARGUMENT_TRAINMODELNR,   "nr", "number", "The number of the model to train (starting from zero)", "number", "");
  initializeParameter(ARGUMENT_AM_CREATE,      "am", "am", "The acoustic model that will be created", "filename", "");
  initializeParameter(ARGUMENT_SATLIST,        "sl", "SAT-list", "A list with training directories containing SAT files", "filename", "");

/////////////////////////////////////// shout_train_finish
// #define ARGUMENT_TRAINDIR
// #define ARGUMENT_AM_CREATE

/////////////////////////////////////// shout_prepare_adapt
// For now not in use...

/////////////////////////////////////// shout_merge_am
//  initializeParameter(ARGUMENT_AM_OUT,
  initializeParameter(ARGUMENT_AM_MERGE1,      "am1", "am-merge-1", "First AM to merge", "filename", "");
  initializeParameter(ARGUMENT_AM_MERGE2,      "am2", "am-merge-2", "Second AM to merge", "filename", "");

/////////////////////////////////////// shout_spkrec
  initializeParameter(ARGUMENT_SPKREC_UBM_MALE,    "ubm", "ubm-male",   "Male UBM model for speaker recognition", "filename", "");
  initializeParameter(ARGUMENT_SPKREC_UBM_FEMALE,  "ubf", "ubm-female", "Female UBM model for speaker recognition", "filename", "");
  initializeParameter(ARGUMENT_SPKREC_TASKLIST,    "tl", "tasklist", "List of runs to score (ID m/f audio-name channel)", "filename", "");
  initializeParameter(ARGUMENT_SPKREC_MODELDIR,    "md", "model-dir", "directory (including task prefix) with speaker models.", "filename", "");
  initializeParameter(ARGUMENT_SPKREC_OUTPUT_PREFIX,"op", "output-prefix", "I will prefix this string to each line of output", "string", "");
  initializeParameter(ARGUMENT_SPKREC_ZNORM_LIST   ,"zn", "z-norm", "A list of speakers to perform Z-norm on", "filename", "");
  initializeParameter(ARGUMENT_SPKREC_OUTPUT_STATS ,"st", "stats", "Output the statistics for each speaker in the tasklist", "define", "");


  initializeParameter(ARGUMENT_NOT_USED ,"nt", "not", "NOT USED", "define", "");


  StringLookup *lookupShort = new StringLookup(ARGUMENTS_NUMBER,parShort);
  StringLookup *lookupLong  = new StringLookup(ARGUMENTS_NUMBER,parLong);

  FILE *globalConfig = fopen(getenv("SHOUT_CONFIG"),"rt");
  if(globalConfig != NULL)
  {
    readConfigFile(globalConfig, lookupLong);
    fclose(globalConfig);
  }


  // Process the argument line, but only to find the config file (rough, but how can we do it otherwise?)
  for(int i=1;i<argc;i++)
  {
    if(strlen(argv[i]) > 1)
    {
      char str[MAXSTR];
      strcpy(str,argv[i]);
      if(i+1 < argc && argv[i+1][0] != '-')
      {
        strcat(str," ");
        strcat(str,argv[i+1]);
        i++;
      }
      if(str[0] == '-')
      {
        if(str[1] == '-')
        {
          if(!fillParameter(&(str[2]),lookupLong))
          {
            USER_ERROR2("Unknown input argument:",str);
          }
        }
        else if(!fillParameter(&(str[1]),lookupShort))
        {
          USER_ERROR2("Unknown input argument:",str);
        }
      }
      else
      {
        USER_ERROR2("Unknown input argument:",str);
      }
    }
  }

  if(parIsDefined(ARGUMENT_CONFIG))
  {
    FILE *configFile = fopen(parameter[ARGUMENT_CONFIG].value,"rt");
    if(configFile != NULL)
    {
      readConfigFile(configFile, lookupLong);
      fclose(configFile);
    }
    else
    {
      printf("%s%s%s\n",WELCOME_HEADER1,CURRENT_VERSION,WELCOME_HEADER2);
      printf("ERROR: The configuration file does not exist or is not readable (%s)\n\n",parameter[ARGUMENT_CONFIG].value);
      exit(0);
    }
  }

// Now the command line arguments are processed:

  for(int i=1;i<argc;i++)
  {
    if(strlen(argv[i]) > 1)
    {
      char str[MAXSTR];
      strcpy(str,argv[i]);
      if(i+1 < argc && argv[i+1][0] != '-')
      {
        strcat(str," ");
        strcat(str,argv[i+1]);
        i++;
      }
      if(str[0] == '-')
      {
        if(str[1] == '-')
        {
          if(!fillParameter(&(str[2]),lookupLong))
          {
            USER_ERROR2("Unknown input argument:",str);
          }
        }
        else if(!fillParameter(&(str[1]),lookupShort))
        {
          USER_ERROR2("Unknown input argument:",str);
        }
      }
      else
      {
        USER_ERROR2("Unknown input argument:",str);
      }
    }
  }


/////////////////////////////// HELP ///////////////////////////////////////////////////
/////////////////////////////// HELP ///////////////////////////////////////////////////
/////////////////////////////// HELP ///////////////////////////////////////////////////

  if(parIsDefined(ARGUMENT_HELP))
  {
  #ifdef PERFORM_CTS_FEATURE_EXTRACTION
    printf("\nWATCH OUT!!!! THIS VERSION IS COMPILED FOR CONVERSATIONAL TELEPHONE SPEECH!!!!! (features go upto 3.8KHz)\n");
    printf("WATCH OUT!!!! THIS VERSION IS COMPILED FOR CONVERSATIONAL TELEPHONE SPEECH!!!!! (features go upto 3.8KHz)\n");
    printf("WATCH OUT!!!! THIS VERSION IS COMPILED FOR CONVERSATIONAL TELEPHONE SPEECH!!!!! (features go upto 3.8KHz)\n");
    printf("WATCH OUT!!!! THIS VERSION IS COMPILED FOR CONVERSATIONAL TELEPHONE SPEECH!!!!! (features go upto 3.8KHz)\n");
    printf("WATCH OUT!!!! THIS VERSION IS COMPILED FOR CONVERSATIONAL TELEPHONE SPEECH!!!!! (features go upto 3.8KHz)\n\n");
  #endif

    printf("%s%s%s\n",WELCOME_HEADER1,CURRENT_VERSION,WELCOME_HEADER2);

    switch(applicationID)
    {
      case APPID_MERGETRAINSET:
        printf("This application will add the training data of traindir2 to the data of the first traindir.\n");
        break;
      case APPID_MAKETRAINSET:
        printf("This application will create a directory with training data for the various acoustic models needed by Shout.\n");
        printf("If the directory already exists, the new data will be added to the available data.\n");
        printf("Currently, four different types of models are supported: SAD, VTLN, DIAR(ization) and PHONE models.\n");
        printf("The META file is needed to determine the warping factors for VTLN. if VTLN is not used, the warping factor column\n");
        printf("should be set to the average (%d).\n",VTLN_NEUTRAL);
        printf("If PHONE models are trained, all phone labels found in the HYP file should be listed in the phone file. Otherwise,\n");
        printf("only the phones available in the phone file will be used for training (the others will be ignored).\n");
        printf("If the train list (audio/meta/hyp) is passed, the audio file and meta file are ignored and the list is used instead\n");
        printf("If 'speaker' is defined, only speech with the label 'speaker' will be processed.\n");
        printf("If 'mel' is defined, the training directory will be filled with binary MELscale feature files and phone-state files\n");
      break;
      case APPID_SEGMENT:
        printf("If no CTS segmentation is required, no AM file is needed.\n");
        printf("Otherwise, the AM file will be used to segment the feature file in the first pass only.\n");
        printf("The UEM file can be used if not all data should be processed, but only the regions defined by the label are used.\n");
      break;
      case APPID_DCT2LEXTREE:
       printf("This application will create a lexical tree (binary dictionary) using a text based dictionary.\n");
       printf("The text based dictionary contains one word per line: the word itself followed by a list of phones.\n");
       printf("WARNING: Pronunciation variations should be placed on sepparate lines, directly after each other!\n");
      break;
      case APPID_PREPROCESS:
        printf("WARNING: This preprocess application is NOT needed for regular use of this toolkit. It's an extra\n");
        printf("         application written for some interesting research.\n");
        printf("Use the audio-list to provide the application with a list of files to process: audio-file meta-file hyp-file\n");
        printf("If you don't provide the phone list or the phone model, only MELscale feature vectors will be generated.\n");
        printf("Otherwise, also a file with state numbers is generated for each feature file.\n");
        printf("Results will be written in the traindir.\n");
        printf("Enjoy ;-)\n");
      break;
      case APPID_PREPARE_ADAPT:
        printf("If the cluster-AM file does not exist, an initial clustering will be done using\n");
        printf("a rough clustering method. The actual cluster-am file can be created with shout_maketrainset\n");
        printf("and shout_train_master/helper. After this a new iteration with shout_prepare_adapt can be started.\n");
        printf("Output is a new segment list (path to the audio file (raw), vtln warping factor and cluster ID)\n");
        break;
      case APPID_UPDATE_VERSION:
        printf("This application will update feature vector files, acoustic models, lexical trees (dictionaries)\n");
        printf("or language models from an older version of shout to the current version. Unfortunately, this may\n");
        printf("not always be possible though. Here is what this version of shout_update_version will do:\n");
        printf("\n");
        printf("The first 4 bytes of each feature vector file or acoustic model needs\n");
        printf("to contain the vector length (dimension) of the vectors in the file. This application will add the\n");
        printf("vector length given as the first argument.\n");
        printf("\n");
        printf("If you provide a random number, the string 'LM' and a language model, some information about the model is written.\n");
      break;
      case APPID_TRAIN_MMI:
        printf("For training acoustic models with the Maximum Mutual Information method\n");
      break;
      case APPID_TRAIN_MODEL:
        printf("Used to train acoustic models:\nThis application will do the actual training.\n");
        printf("It will train the specified model and write the result in the training directory.\n");
        printf("The application shout_train_finish can be used to create the final combined AM\n"); 
        printf("For silence phones, the training time can be reduced by setting train-sil-p and train-sil-max. If set,\n");
        printf("each sample will be used with probability train-sil-p until train-sil-max samples are processed.\n");
        printf("If train-sil-max < 0, the entire set will be processed.\n");
        printf("If an existing acoustic model is defined, this model will be used as a starting point for training the new model.\n");
      break;
      case APPID_TRAIN_FINISH_SAT:
        printf("Finish SAT training... (docs TODO!)\n");
      break;
      case APPID_TRAIN_FINISH:
        printf("Used to train acoustic models:\nThis application will combine all models trained with shout_train_model into one AM\n");
        printf("The training directory is created with shout_maketrainset\n");
      break;
      case APPID_VTLN:
        printf("This application determines the VTLN warping factor for each segment in the META file.\n");
        printf("It can do this either by processing each segment (SEGMENT) sepperately, or by combining all\n");
        printf("segments of each speaker (SPEAKER) or by creating one warping factor out of all lines (ALL).\n");
      break;
      case APPID_ADAPT_AM:
        printf("This application will adapt the input acoustic models on the data in the training directory using SMAPLR\n");
        printf("The resulting AM will be stored to disc. Except when '-sat' is defined. In that case the output AM is\n");
        printf("assumed to be another training directory. This directory will then be prepared for SAT training.\n");
      break;
      case APPID_CLUSTER:
        printf("This application performs speaker clustering based on a speaker diarization algorithm that is still under development.\n");
      break;
      case APPID_AVG_ENERGY:
        printf("This application determines the average energy in the audio file, only measuring the meta-in segments.\n");
      break;
      case APPID_EXPAND_VOCAB:
        printf("This is an experimental application...\n");
      break;
      case APPID_SPKREC:
        printf("This application will score a list of identification runs (tasklist) given a UBM and all speaker models.\n");
      break;
      case APPID_SPKRECSTATS:
        printf("This application will output the training statistics for JFA speaker recognition.\n");
      break;
      case APPID_SHOUT:
      case APPID_SHOUTONLINE:
      case APPID_LM2BIN:
      case APPID_MERGE_AM:
      case APPID_NORMALIZE_AM:
      default:
      break;
    }

    if(strcmp(parameter[ARGUMENT_HELP].value,"ALL") == 0)
    {
      printf("\nThe following is a list with ALL parameters used in the Shout toolkit.\n");
    }
    else
    {
      printf("\nThe following is a list with all parameters that influence the behaviour of this application.\n");
    }
    printf("<parameter value> - the parameter may not be omitted\n");
    printf("[parameter value] - the parameter is optional\n");
    printf("[define]          - 'define' means that any value will enable this parameter\n");
    printf("The parameters may be defined in the global configuration file ($SHOUT_CONFIG),\nin the config file (-c [config file]) or on the argument line.\n\n");


    for(int i=0;i<ARGUMENTS_NUMBER;i++)
    {
      int needPar = parInApp(i);
      if(strcmp(parameter[ARGUMENT_HELP].value,"ALL") == 0 || (needPar>0))
      {
        printf("--%s (-%s) ",parLong[i],parShort[i]);
        for(int i2=7+strlen(parLong[i])+strlen(parShort[i]);i2<30;i2++)
        {
          printf(" ");
        }
        if(needPar == 1)
        {
          printf("[%s] ",parameter[i].descrVal);
        }
        else
        {
          printf("<%s> ",parameter[i].descrVal);
        }
        for(int i2=strlen(parameter[i].descrVal);i2<25;i2++)
        {
          printf(" ");
        }
        printf("%s\n",parameter[i].description);
      }
    }
    printf("\n");
  #ifdef PERFORM_CTS_FEATURE_EXTRACTION
    printf("WATCH OUT!!!! THIS VERSION IS COMPILED FOR CONVERSATIONAL TELEPHONE SPEECH!!!!! (features go upto 3.8KHz)\n");
    printf("WATCH OUT!!!! THIS VERSION IS COMPILED FOR CONVERSATIONAL TELEPHONE SPEECH!!!!! (features go upto 3.8KHz)\n");
    printf("WATCH OUT!!!! THIS VERSION IS COMPILED FOR CONVERSATIONAL TELEPHONE SPEECH!!!!! (features go upto 3.8KHz)\n");
    printf("WATCH OUT!!!! THIS VERSION IS COMPILED FOR CONVERSATIONAL TELEPHONE SPEECH!!!!! (features go upto 3.8KHz)\n");
    printf("WATCH OUT!!!! THIS VERSION IS COMPILED FOR CONVERSATIONAL TELEPHONE SPEECH!!!!! (features go upto 3.8KHz)\n\n");
  #endif
    exit(0);
  }
/////////////////////////////// CONFIG ///////////////////////////////////////////////////
/////////////////////////////// CONFIG ///////////////////////////////////////////////////
/////////////////////////////// CONFIG ///////////////////////////////////////////////////
  else if(parIsDefined(ARGUMENT_CREATE_CONFIG))
  {
    printf("%s%s%s\n",WELCOME_HEADER1,CURRENT_VERSION,WELCOME_HEADER2);
    for(int i=0;i<ARGUMENTS_NUMBER;i++)
    {
      if(i != ARGUMENT_CREATE_CONFIG && i != ARGUMENT_HELP && i != ARGUMENT_CONFIG)
      {
        if(strcmp(parameter[ARGUMENT_CREATE_CONFIG].value,"ALL") == 0 || (parInApp(i)>0))
        {
          if(strcmp(parameter[i].value,"") == 0)
          {
            printf("# (not enabled) ");
          }
          printf("%s",parLong[i]);
          for(int i2=strlen(parLong[i]);i2<50;i2++)
          {
            printf(" ");
          }
          printf("%s\n",parameter[i].value);
        }
      }
    }
    exit(0);
  }
  else // Checking if all needed variables are available...
  {
    bool allParsAvailable = true;
    for(int i=0;i<ARGUMENTS_NUMBER;i++)
    {
      if(i != ARGUMENT_CREATE_CONFIG && i != ARGUMENT_HELP && i != ARGUMENT_CONFIG)
      {
        if((parInApp(i)>1) && (strcmp(parameter[i].value,"") == 0))
        {
          if(allParsAvailable)
          {
            printf("%s%s%s\n",WELCOME_HEADER1,CURRENT_VERSION,WELCOME_HEADER2);
            allParsAvailable = false;
          }
          printf("ERROR: No value defined for the parameter '%s'\n",parLong[i]);
        }
      }
    }
    if(!allParsAvailable)
    {
      printf("\nFor more information use the '--help' parameter. Bye, bye...\n\n");
      exit(0);
    }
  }

  delete lookupShort;
  delete lookupLong;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

int ShoutConfig::parInApp(int i)  // Return 0 if not, 1 if optional, 2 if needed.
{
  switch(applicationID)
  {
    case APPID_ADAPT_AM:
    {
      int incl[4] = {ARGUMENT_TRAINDIR,ARGUMENT_AM_PHONE,ARGUMENT_AM_CLUSTER,ARGUMENT_AM_OUT};
      int need[4] = {2,2,1,2};
      int a = 3; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_CLUSTER:
    {
      int incl[9] = {ARGUMENT_AUDIOFILE,ARGUMENT_LABEL,ARGUMENT_META_IN,ARGUMENT_META_OUT,ARGUMENT_WIDEN_SEGMENTS,ARGUMENT_SECOND_STREAM,ARGUMENT_THIRD_STREAM,ARGUMENT_FAST_MERGE,ARGUMENT_MAX_CLUSTERS};
      int need[9] = {2,2,2,2,1,1,1,1,1};
      int a = 8; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_DCT2LEXTREE:
    {
      int incl[3] = {ARGUMENT_DCT_BIN,ARGUMENT_PHONE_LIST,ARGUMENT_DCT_TEXT};
      int need[3] = {2,2,2};
      int a = 2; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_LM2BIN:
    {
      int incl[3] = {ARGUMENT_DCT_BIN,ARGUMENT_LM_BIN,ARGUMENT_LM_ARPA};
      int need[3] = {2,2,2};
      int a = 2; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_MAKETRAINSET:
    {
      int incl[11] = {ARGUMENT_AUDIO_LIST,ARGUMENT_AUDIOFILE,ARGUMENT_PHONE_LIST,ARGUMENT_HYPFILE,ARGUMENT_META_IN,ARGUMENT_TRAINDIR,ARGUMENT_AM_TYPE,ARGUMENT_DECISIONTREE,ARGUMENT_PERLINE,ARGUMENT_SPEAKER_FILTER,ARGUMENT_TRAINSILP};
      int need[11] = {1,1,2,1,1,2,2,1,1,1,1};
      int a = 10; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_MERGE_AM:
    {
      int incl[3] = {ARGUMENT_AM_OUT,ARGUMENT_AM_MERGE1,ARGUMENT_AM_MERGE2};
      int need[3] = {2,2,2};
      int a = 2; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_MERGETRAINSET:
    {
      int incl[4] = {ARGUMENT_AM_TYPE,ARGUMENT_PHONE_LIST,ARGUMENT_TRAINDIR,ARGUMENT_TRAINDIR2};
      int need[4] = {2,2,2,2};
      int a = 3; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_AVG_ENERGY:
    {
      int incl[2] = {ARGUMENT_AUDIOFILE,ARGUMENT_META_IN};
      int need[2] = {2,2};
      int a = 1; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_SHOUT:
    {
      int incl[29] = {ARGUMENT_AUDIOFILE,ARGUMENT_META_IN,ARGUMENT_HYPFILE,ARGUMENT_AM_PHONE,ARGUMENT_LM_BIN,ARGUMENT_DCT_BIN,ARGUMENT_PHONEBASED,ARGUMENT_LATTICE,ARGUMENT_NBEST,ARGUMENT_ALIGN,ARGUMENT_LATTICERESCORE,ARGUMENT_ANALYSE,ARGUMENT_TOKENDIST,ARGUMENT_PHONECONF,ARGUMENT_BACKDCT_BIN,ARGUMENT_TRAINDIR,ARGUMENT_XML,ARGUMENT_LM_SCALE,ARGUMENT_TRANS_PENALTY,ARGUMENT_SIL_PENALTY,ARGUMENT_BEAM,ARGUMENT_STATE_BEAM,ARGUMENT_END_STATE_BEAM,ARGUMENT_LMLA,ARGUMENT_HISTOGRAM_STATE_PRUNING,ARGUMENT_HISTOGRAM_PRUNING,ARGUMENT_PERLINE,ARGUMENT_AM_HISTNORM,ARGUMENT_LATTICETRAINING};
      int need[29] = {2,2,1,2,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
      int a = 28; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_SPKREC:
    {
      int incl[7] = {ARGUMENT_SPKREC_UBM_MALE,ARGUMENT_SPKREC_UBM_FEMALE,ARGUMENT_SPKREC_TASKLIST,
                     ARGUMENT_SPKREC_MODELDIR,ARGUMENT_SPKREC_OUTPUT_PREFIX,ARGUMENT_SPKREC_ZNORM_LIST,ARGUMENT_SPKREC_OUTPUT_STATS};
      int need[7] = {2,2,2,1,1,1,1};
      int a = 6; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_SPKRECSTATS:
    {
      int incl[4] = {ARGUMENT_SPKREC_UBM_MALE,ARGUMENT_AUDIOFILE,ARGUMENT_META_IN,ARGUMENT_SPKREC_OUTPUT_STATS};
      int need[4] = {1,2,1,1};
      int a = 3; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_SHOUTONLINE:
    {
      int incl[5] = {ARGUMENT_AM_PHONE,ARGUMENT_LM_BIN,ARGUMENT_DCT_BIN,ARGUMENT_LATTICE,ARGUMENT_CONFIGDIR};
      int need[5] = {2,1,2,1,1};
      int a = 4; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_TRAIN_MMI:
    {
      int incl[4] = {ARGUMENT_AM_PHONE,ARGUMENT_AM_OUT,ARGUMENT_LATTICETRAINING_ENUM,ARGUMENT_LATTICETRAINING_DENOM};
      int need[4] = {2,2,2,2};
      int a = 3; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    case APPID_TRAIN_MODEL:
    {
      int incl[8] = {ARGUMENT_TRAINDIR,ARGUMENT_DECISIONTREE,ARGUMENT_MAXGAUSSIANS,ARGUMENT_TRAINSILP,ARGUMENT_TRAINSILMAX,ARGUMENT_TRAINMODELNR,
                     ARGUMENT_AM_CREATE,ARGUMENT_DOFAST};
      int need[8] = {2,1,2,1,1,2,1,1};
      int a = 7; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_SEGMENT:
    {
      int incl[9] = {ARGUMENT_AUDIOFILE,ARGUMENT_AM_SEGMENT,ARGUMENT_META_OUT,ARGUMENT_LABEL,ARGUMENT_UEM,ARGUMENT_WIDEN_SEGMENTS,ARGUMENT_CTS_IN,ARGUMENT_CTS_OUT,ARGUMENT_ENERGY_BASED};
      int need[9] = {2,2,2,1,1,1,1,1,1};
      int a = 8; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_TRAIN_FINISH_SAT:
    {
      int incl[4] = {ARGUMENT_TRAINDIR,ARGUMENT_AM_CREATE,ARGUMENT_TRAINMODELNR,ARGUMENT_SATLIST};
      int need[4] = {2,2,2,2};
      int a = 3; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_TRAIN_FINISH:
    {
      int incl[2] = {ARGUMENT_TRAINDIR,ARGUMENT_AM_CREATE};
      int need[2] = {2,2};
      int a = 1; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_VTLN:
    {
      int incl[5] = {ARGUMENT_AUDIOFILE,ARGUMENT_META_IN,ARGUMENT_META_OUT,ARGUMENT_AM_VTLN,ARGUMENT_PERLINE};
      int need[5] = {2,2,2,2,1};
      int a = 4; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_EXPAND_VOCAB:
    {
      int incl[3] = {ARGUMENT_DCT_BIN, ARGUMENT_BACKDCT_BIN, ARGUMENT_HYPFILE};
      int need[3] = {2,2,2};
      int a = 2; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_NORMALIZE_AM:
    {
      int incl[2] = {ARGUMENT_AUDIOFILE, ARGUMENT_META_IN};
      int need[2] = {2,2};
      int a = 1; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_PREPROCESS:
    {
      int incl[4] = {ARGUMENT_PHONE_LIST, ARGUMENT_TRAINDIR, ARGUMENT_AUDIO_LIST, ARGUMENT_AM_PHONE};
      int need[4] = {1,2,2,1};
      int a = 3; while(a>=0 && incl[a] != i){ a--;} if(a>=0){return need[a];} else {return 0;}
    }
    break;
    case APPID_UPDATE_VERSION:
    case APPID_PREPARE_ADAPT:
    default:
    {
      return 0;
    }
    break;
  }
  return false;
}








/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

ShoutConfig::~ShoutConfig()
{
  for(int i=0;i<ARGUMENTS_NUMBER;i++)
  {
    delete[] parShort[i];
    delete[] parLong[i];
  }
  delete[] parShort;
  delete[] parLong;
  delete[] parameter;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void ShoutConfig::readConfigFile(FILE *config, StringLookup *lookupLong)
{
  char str[MAXSTR];
  int line = 0;

  while(!feof(config))
  {
    str[0]  = 0;
    line++;
    char *resChar = fgets(str,MAXSTR,config);
    assert(resChar != NULL);
    if(strlen(str) > 0 && str[0] != '#' && str[0] != '\n')
    {
      if(!fillParameter(str,lookupLong))
      {
        USER_ERROR2_NUM("The configuration file is invalid at line",line);
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool ShoutConfig::fillParameter(char *line, StringLookup *lookupLong)
{
  char str2[MAXSTR];
  str2[0] = 0;
  RightTrim(line);
  if(strlen(line) > 0)
  {
    int parID = lookupLong->getWordID(line);
    if(parID < 0)
    {
      splitList(line,str2);
      parID = lookupLong->getWordID(line);
    }
    else
    {
      strcpy(str2,"YES");
    }
    if(parID < 0)
    {
      return false;
    }
    if(strlen(str2) > MAXSTR)
    {
      USER_ERROR2("Sorry, the following parameter value is too long:",str2);
    }
    strcpy(parameter[parID].value,str2);
  }
  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void ShoutConfig::initializeParameter(int id, const char *shortDiscr,const char *longDiscr,const char *text,const char *valDescr,const char *val)
{
  parShort[id] = new char[5];  strncpy(parShort[id],shortDiscr,5);
  parLong[id]  = new char[50]; strncpy(parLong[id], longDiscr,50);
  strcpy(parameter[id].description, text);
  strcpy(parameter[id].descrVal,    valDescr);
  strcpy(parameter[id].value,       val);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

void ShoutConfig::initializeParameter(int id, const char *shortDiscr,const char *longDiscr,const char *text,const char *valDescr, double val)
{
  char tmp[100];
  sprintf(tmp,"%10.2f",val);
  initializeParameter(id,shortDiscr,longDiscr,text,valDescr,tmp);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

char *ShoutConfig::getStringValue(int par)
{
  assert(parInApp(par));
  return parameter[par].value;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

float ShoutConfig::getFloatValue(int par)
{
  assert(parInApp(par));
  return atof(parameter[par].value);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \todo docs
/////////////////////////////////////////////////////////////////////////////////////////////////////

bool ShoutConfig::parIsDefined(int par)
{
  return ((strlen(parameter[par].value) > 0) && (strcmp(parameter[par].value,"NO") != 0) && (strcmp(parameter[par].value,"OFF") != 0) && (strcmp(parameter[par].value,"FALSE") != 0));
}

