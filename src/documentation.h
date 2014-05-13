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

// gmake -f Makefile.cvs
// mkdir optimized
// cd optimized
// CXXFLAGS="-O3 -funroll-loops -march=pentium4 -malign-double -mfpmath=sse -msse -msse2" ../configure
// gmake -j1

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \mainpage Large Vocabulary Continuous Speech Recognition
///
/// During my PhD research I have developed a large vocabulary continuous speech recognition toolkit that I named SHoUT.
/// SHoUT is a Dutch acronym for: 'Speech Recognition Research at the University of Twente' (I could have called the toolkit 'SRRatUoT', but 'SHoUT' just sounds better).
/// I think (read: 'hope') that the toolkit is now mature enough to be used by other researchers.
/// You can download the latest release from the \ref download "download page". All comments and suggestions are welcome!
/// For information on how to use the toolkit, have a look at the \ref user_manual "manual" (for users) and
/// at the \ref programmers_manual "API reference" (for programmers). For background information (architecture, algorithms etc.) please read my <a href="img/thesis_Marijn_Huijbregts.pdf">thesis</a>.
/// <img src="img/comic-small.png" alt="" class="center" />
/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \page contact Contact
///
/// You can reach me for questions about SHoUT at: marijn.huijbregts@gmail.com
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \page download Download
///
/// Shout is a research tool that I have mainly written for myself. It is *not* a commercial tool that you can use out of the box!
/// If you are not a researcher on ASR, this tool might not be very useful to you. Having said that... Have fun!
///
/// <div class="meta"></div>
///
/// \section version_0_3 Version 0.3 released at 01/12/2010
/// <A HREF="http://wwwhome.cs.utwente.nl/~huijbreg/shout/download/shout_0_2_complete.tar.gz">Source and Linux binaries tar-ball</A>
///
/// <div class="meta"></div>
///
/// \section version_0_2 Version 0.2 released at 01/12/2008
/// <A HREF="http://wwwhome.cs.utwente.nl/~huijbreg/shout/download/shout_0_2_complete.tar.gz">Source and Linux binaries tar-ball</A>
///
/// <div class="meta"></div>
///
/// \section version_0_1 Version 0.1 released at 06/11/2007
/// <A HREF="http://wwwhome.cs.utwente.nl/~huijbreg/shout/download/shout_0_1_complete.tar.gz">Source and Linux binaries tar-ball</A>
/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \page programmers_manual Reference guide
/// This is your starting point for programming in or with the toolkit. Class reference can be found in the menu above.
///
/// I use the great application "kdevelop" to edit and build my source code, but if
/// you want to build the source code manually, this is how kdevelop does it. Start in the main shout directory and type:
/// - gmake -f Makefile.cvs
/// - mkdir optimized
/// - cd optimized
/// - CXXFLAGS="-O3 -funroll-loops -march=pentium4 -malign-double -mfpmath=sse -msse -msse2" ../configure
/// - gmake -j1
///
/// <div class="meta"></div>
///
/// \section building_blocks The building blocks of Shout
/// This software package contains multiple applications. A short description of each application is given below. If you just want
/// to use the toolkit and you do not have the intention to develop yourself, it is better to read the \ref user_manual "user manual".
///
/// <div class="meta"></div>
///
/// \section app_shout shout
/// The hart of the software package, the decoder, is called shout. This is where all the models come together... During the early
/// days of the development the decoder was called 'whisper'. Unfortunately another decoder with the same name already existed.
/// I have changed the decoder name, but the main class that is handling the top-level recognition is still called Whisper.
/// Whisper will load all needed models. After that, most work is done by the LexicalTree class.
///
/// <div class="meta"></div>
///
/// \section app_shout_adapt_am shout_adapt_am
/// This application reads an acoustic model file and a training/adapting phone directory and creates a new acoustic model file
/// that is adapted to the training data using the Structured Maximum a Posteriori Linear Regression (SMAPLR) method. The main class
/// of shout_adapt_am is Adapt_AM. 
///
/// <div class="meta"></div>
///
/// \section app_shout_cluster shout_cluster
/// This is the speaker diarization application. See Shout_Cluster.
///
/// <div class="meta"></div>
///
/// \section app_shout_dct2lextree shout_dct2lextree
/// This application translates a dictionary file (text file containing rows of: word - tabs or spaces - phone pronunciation) into a 
/// lexical Prefix Tree, also called a Pronunciation Prefix Tree (PPT), suitable for the decoder to read. The main class of 
/// shout_dct2lextree is Shout_dct2lextree.
///
/// <div class="meta"></div>
///
/// \section app_shout_lm2bin shout_lm2bin
/// Shout can handle uni-, bi-, tri- and four-gram ARPA language models (depending on how the distribution is compiled).
/// The application shout_lm2bin will read an ARPA LM and translate it to a binary format suitable
/// for the decoder. The main class of shout_lm2bin is Shout_lm2bin.
///
/// <div class="meta"></div>
///
/// \section app_shout_maketrainset shout_maketrainset
/// Shout_maketrainset will read an hypothesis file (the output of Shout in native format) or a Master Label File (MLF) and will
/// store all phones in the training directory. This data directory
/// is later used by the training application (shout_train_master). The main class of shout_maketrainset is Shout_MakeTrainSet.
/// It is possible to create a training directory for either ASR models, SAD models, diarization models or VTLN models.
///
/// <div class="meta"></div>
///
/// \section app_shout_merge_am shout_merge_am
/// In case you don't want to re-train all phones, with this application you can choose phone models from two AM files to store
/// in a new AM file. See ShoutMergeAm.
///
/// <div class="meta"></div>
///
/// \section app_shout_segment shout_segment
/// This is the speech/non-speech detector. The main class is ShoutSegment.
///
/// <div class="meta"></div>
///
/// \section app_shout_vtln shout_vtln
/// Determines the VTLN warping factor. The main class of shout_vtln is Shout_VTLN.
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \page user_manual User manual
///
/// SHoUT is written in C++ on a Linux platform. The code should run on other platforms (Windows, Mac) without adjustments,
/// but I do not support the porting. You can build the toolkit by running the script configure-make.sh. The binaries will appear in release/src.
/// Installation (by root) of the binaries is not needed.
///
/// The various applications of the Shout toolkit will be discussed in the following sections. Each section
/// describes a typical use case. For the precise syntax of all application parameters, just run the
/// applications with the parameter --help (or -h). A short help text will appear.
///
/// <div class="meta"></div>
///
/// \section input_output Input/output
/// The shout toolkit does not support different audio and video file types. You will need to convert your multimedia files into raw audio files
/// yourself. Input should be headerless raw audio: mono, 16 KHz, 16 bits little endian.
///
/// Each application will use a 'meta-data' file for input and/or output. On each line of this file, an audio segment is specified. What is done with these segments depends on the application that uses the file. The format of each line of a meta-data file is as follows:
/// - SPEAKER [label] [VTLN factor] [begin time] [length] \<NA\> \<NA\> [SPK ID] \<NA\> [\<s\> \<s\>  followed by a word based transcription]
///
/// - label       - is the identifier of the file.
/// - VTLN factor - is the factor calculated by shout_vtln used for VTLN.
/// - begin time  - begin time of the audio segment
/// - length      - length of the segment
/// - \<NA\>      - Not applicable (for compatibility with NIST RTTM files).
/// - SPK ID      - label of the segment (can be SPEECH/SIL or SPK ID).
///
/// I use these meta-data files so that it is not needed to cut-up the audio files in actual segments. Each application simply reads the audio file and uses the bits of it that are defined in the meta-data file.
///
/// <div class="meta"></div>
///
/// \section use_cases Typical use-cases
/// For my work on spoken document retrieval I run the system as depicted in the figure below. First, I run
/// Speech Activity Detection (SAD) to find the segments in the audio that contain speech. Next, I run diarization
/// to determine: 'who spoke when?'. Third, I run Vocal Tract Length Normalization (VTLN) to normalize the speech of
/// each speaker for his/her vocal tract length. After this I perform the actual decoding. But before decoding is possible,
/// first acoustic models, a dictionary and a language model need to be created. After the first decoding iteration it is
/// possible to perform acoustic model adaptation and run a second decoding iteration. All these steps are explained in the use cases below.
///
/// <img src="img/system-overview-small.png" alt="" class="center" />
/// - \ref use_case_sad
/// - \ref use_case_diarization
/// - \ref use_case_vtln
/// - \ref use_case_dct_lm
/// - \ref use_case_training
/// - \ref use_case_recognition
/// - \ref use_case_adapt
/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \page use_case_dct_lm Preparing your dictionary and language model for Shout
///
/// For decoding, three binary files are needed: a lexical tree file, a language model file and an 
/// acoustic model file. The acoustic models need to be trained by Shout. The language model
/// and lexical tree are created using the applications shout_dct2lextree and shout_lm2bin.
///
/// \section preparing_dictionary Preparing a dictionary
/// The application shout_dct2lextree needs two input files and will output a binary 'lexical tree'
/// file. The input phone list file consists of the entire list of phones. The first two lines of this format define
/// the total number of phone and non-speech models. It uses the following syntax:
/// - "Number of phones:" [number of phone models]
/// - "Number of SIL's:" [number of non-speech models]
/// - One non-speech model name per line (times the specified number of non-speech models)
/// - One phone model name per line (times the specified number of phone models)
///
/// The pronunciation dictionary contains one word per line followed by a string of phone or non-speech model names 
/// separated by one or more spaces. Make sure that the first two lines in your DCT are as follows:
/// - \<s\>   SIL
/// - \</s\>  SIL
///
/// <div class="meta"></div>
///
/// \section preparing_lm Preparing a language model
/// The application 'shout_lm2bin' needs a lexical tree file (the output file of shout_dct2lextree)
/// and an ARPA language model. It will create a binary language model file fit for shout.
/// Currently, unigram, bigram, trigram and four-gram language models are supported.
/// During compile time the system can be optimized by setting a maximum depth (from bi- to four-grams) with the
/// parameter LM_NGRAM_DEPTH in standard.h The standard setting of LM_NGRAM_DEPTH is for trigrams.
/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \page use_case_training Setting up the system for training your acoustic models
///
/// Training new acoustic files is done in three repetitive steps:
/// - Create an hypothesis file in native shout format. Alignments in other formats need to be re-formatted.
/// - Make a training set using the alignment file.
/// - Run shout_train_am for each phone and run shout_train_finish to combine all phones in a single file. 
/// The resulting AM can be used in the first step to create a better alignment.
///
/// <div class="meta"></div>
///
/// \section create_alignment Create an alignment of the training audio
/// In order to train new models, a set of training phones are needed. These training sets can only be
/// created when an alignment of the audio and the speech transcription on phone level is available.
/// This alignment must be in native shout format.
///
/// A shout alignment can be created by calling shout with a meta-data file in the following format
/// - SPEAKER [label] [VTLN factor] [begin time] [length] \<NA\> \<NA\> [SPK ID] \<NA\> \<s\> \<s\> [word based transcription]
/// Note that this file format is basically the RTTM format from NIST.
///
/// The \<s\> symbol represents silence. The first two symbols in the transcription are not used to align the audio,
/// but only to create a language model history. When two \<s\> symbols are added, the language model history is reset to 'start of sentence'.
///
/// <div class="meta"></div>
///
/// \section create_trainset Make a training set
/// Once the alignment file is created, a training directory can be made by the application shout_maketrainset.
/// This application uses a phone-set file, it needs to know the path to the audio files (all files
/// need to be in the same directory), the alignment file (HYP file, printed to stdout by shout) and the
/// path to the training set directory that will be created.
///
/// Alternatively, shout_maketrainset can use a phone-set file and a file with a list of audio-metadata-hyp files.
///
/// In order to create acoustic phone models, make sure to set the type to PHONE in shout_maketrainset.
///
/// <div class="meta"></div>
///
/// \section train_phones Training the phones
/// You can train each individual phone with shout_train_am. The training procedure is split-up like this so that you
/// can run the training process in parallel. See shout_train_am -h for more info.
///
/// Once all phones are trained, you can use shout_train_finish to combine the models into one single binary file.
/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \page use_case_recognition Decoding
///
/// The one and only reason why all other applications are developed is being able to decode!
/// That's why the decoder, the hart of the toolkit, is just called... Shout!
///
/// Run ./shout with the output meta-data file of shout_vtln (or of shout_cluster if no VTLN is needed).
/// Next is a short description of the most important parameters. Please run ./shout -h for more help.
///
/// \section decoder_settings_model Model settings
/// The decoder needs a language model file (lm), acoustic model file (amp) and a lexical tree file (dct). 
/// All files should be binary files created by the shout toolkit.
///
/// <div class="meta"></div>
///
/// \section decoder_settings_search Search settings
/// The search space of the decoder is restricted using five parameters. If these paramters are not assigned a value, the default values (shown when shout is started with -cc) will be used.
///
/// The five search restriction parameters:
/// - BEAM                     (floating point number)
/// - STATE_BEAM               (floating point number)
/// - END_STATE_BEAM           (floating point number)
/// - HISTOGRAM_STATE_PRUNING  (positive number)
/// - HISTOGRAM_PRUNING        (positive number)
///
/// <div class="meta"></div>
///
/// \section decoder_settings_amlm AM and LM scaling settings
/// The most likely paths in the jungle of feature vectors are calculated using a language model and acoustic models.
/// The scaling between the two types of models influences the outcome of the trip through this jungle. This scaling is 
/// set using three parameters in the formula: 
///
/// Score(LM_SCALE,TRANS_PENALTY,SIL_PENALTY) = ln(AMSCORE) + LM_SCALE*lm(LMSCORE) + TRANS_PENALTY*NR_WORDS + SIL_PENALTY*NR_SIL
///
/// Shout has implemented an efficient method of incorporating the LM score in the search. This method, Language Model Look-Ahead, 
/// is switched on by default, but it can be toggled on or off in the configuration file.
/// - LM_SCALE                 (floating point number)
/// - TRANS_PENALTY            (floating point number)
/// - SIL_PENALTY              (floating point number)
/// - LMLA                     (1=on, 0=off)
///
/// <div class="meta"></div>
///
/// \section decoder_settings_alignment Alignment
/// You can specify a special background dictionary if you want to perform alignment with OOV marking. For performing alignment instead
/// of ASR, simply set the forced-alignment parameter (see ./shout -h). Make sure to add the utterance to align in the meta-data file, starting
/// with \<s\> \<s\>. See \ref use_case_training "the training use-case" for more information.
///
/// <div class="meta"></div>
///
/// \section decoder_settings_output Type of output (text or XML)
/// - XML       (output will be in XML format)
///
/// <div class="meta"></div>
///
/// \section decoder_settings_lattice Lattice output
/// Shout can generate lattices and output them in PSFG format. If this is wanted (decoding will be a bit slower)
/// the lattice parameter should be used, with a path to the directory where the lattice
/// files need to be written to. One lattice will be created for each line in the meta-data file. The lattice files will be
/// named after the label of each audio file (defined in the meta-data file). 
///
/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \page use_case_sad Speech/non-speech segmentation
///
/// The decoder can only handle audio containing solely speech. Some silence models may be trained
/// (like silence, lip smack, etc) but when audio contains other sources like for example  
/// music, jingles or formula-one cars, it is a good idea to segment the audio before feeding it to
/// the decoder. Shout provides this possibility with the shout_segment application. This application can
/// cluster audio in three categories: "speech", "silence" and "audible non-speech".
///
/// Refer to my thesis for detailed information on how the application works (see menu). In short
/// this is what happens:
/// - Using a speech/silence AM, an initial segmentation is created.
/// - Iteratively, using the high confidence fragments of the initial segmentation, three new GMMs are trained: "speech", "silence" and "audible non-speech"
/// - Using BIC, the application checks if the "speech" and "audible non-speech" models are actually different from each other. If they are identical, all three models
/// are discarded and a new set of GMMs, only containing "speech" and "silence", is trained iteratively.
/// - A final Viterbi run determines the final alignment and clustering.
///
/// <div class="meta"></div>
///
/// \section train_sad Training the GMM for the initial training
/// In a future implementation I'm planning to change the first step where models are needed for creating the initial segmentation to a method that does not require models at all.
/// But for now, the speech and non-speech models needs to be trained. Training these GMMs is done the same as \ref use_case_training "training phones", except that 
/// it is not needed to perform multiple training iterations (especially when the initial alignment is the final alignment of the training step for phones).
///
/// Training new speech/non-speech models is done in three steps:
/// - Create an hypothesis file in native shout format. Alignments in other formats need to be re-formatted.
/// - Make a training set using the alignment.
/// - Run shout_train_am and shout_train_finish. The resulting AM can be used in the first step to create
/// a better alignment.
///
/// <div class="meta"></div>
///
/// \section sad_create_alignment Create an alignment of the training audio
/// For creating the alignment see the section on \ref create_alignment "training acoustic models".
///
/// <div class="meta"></div>
///
/// \section sad_create_trainset shout_maketrainset
/// You can create a training set with the application make_trainset. Be sure to choose SAD for the type of trainset you want to create.
/// Run shout_maketrainset with the -h option for usage of this application.
///
/// Note 1: The feature vectors used for creating these models are different from the phone feature vectors. That's why it is not possible to use the training directory for phones.
///
/// Note 2: Not all phone occurrences are used for training. Only the once that are not directly next to a silence (non-speech) phone.
///
/// <div class="meta"></div>
///
/// \section sad_train_am shout_train_am and shout_train_finish
/// Once you have a training set for your SAD-AM, you can run shout_train_am to train your SIL and SPEECH models (run the application twice) and after that you can run shout_train_finish to generate a binary AM file from the two models.
/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \page use_case_diarization Speaker Diarization
///
/// Speaker diarization is the task of: 'Who spoke when?'. If you have a speech/non-speech meta-data file
/// from shout_segment it is very easy to perform diarization using the application shout_cluster.
/// Type ./shout_cluster -h for more information.
/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \page use_case_adapt Adapting your acoustic models (SMAPLR)
///
/// \todo Under development. Come back later :-)
/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// \page use_case_vtln Vocal tract length normalization
///
/// Vocal Tract Length Normalization (VTLN) is a normalization step that improves ASR considerably.
/// Simply run shout_vtln on the output of the diarization step (shout_cluster) to get the normalization
/// factors for each speaker (incorporated in the output meta-data file).
/////////////////////////////////////////////////////////////////////////////////////////////////////

