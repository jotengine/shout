INCLUDES = $(all_includes)

# these are the headers for your project
noinst_HEADERS = adapt_am.h adapt_am_treenode.h articulatorystream.h \
	featureextraction.h featurepool.h io_speaker_segmenter.h memmappedfile.h mergetrainset.h \
	multimixgaussian.h nbest.h normalizeam.h phonefilereader.h shout_cluster.h \
	shout_dct2lextree.h shout_lm2bin.h shout_maketrainset.h shout_updateversion.h shout_vtln.h \
	shoutconfig.h shoutmergeam.h shoutonline.h shoutprepareadapt.h shouttrainfinish.h \
	shouttrainfinishsat.h shouttrainmmi.h shouttrainmodel.h speakerrecognition.h stringlookup.h \
	train_segmenter.h trainphonemodel.h whisper.h


#########################################################################
# APPLICATION SECTION
#########################################################################
# this is the program that gets installed.  it's name is used for all
# of the other Makefile.am variables
bin_PROGRAMS = shout shout_adapt_am shout_avgEnergy shout_cluster \
	shout_dct2lextree shout_lm2bin shout_maketrainset shout_merge_am shout_merge_trainset \
	shout_normalize_am shout_online shout_prepare_adapt shout_preprocess shout_segment \
	shout_spkrec shout_spkrec_stats shout_train_finish shout_train_finish_sat \
	shout_train_mmi shout_train_model shout_update_version shout_vtln

# the application source, library search path, and link libraries

shout_dct2lextree_LDFLAGS = $(all_libraries)
shout_dct2lextree_SOURCES = socket_standard.h standard.h shout_dct2lextree.cpp shout_maketrainset.cpp shout_maketrainset.h
shout_lm2bin_LDFLAGS = $(all_libraries)
shout_lm2bin_SOURCES = memmappedfile.cpp shout_lm2bin.cpp \
	shout_maketrainset.cpp shout_maketrainset.h standard.h trainhash.cpp trainhash.h

shout_maketrainset_LDFLAGS = $(all_libraries) -pg
shout_maketrainset_SOURCES = shout_maketrainset.cpp shout_maketrainset_main.cpp

shout_LDFLAGS = $(all_libraries)
shout_SOURCES = standard.h whisper.cpp whisper.h shout_maketrainset.cpp \
		shout_maketrainset.h adapt_am_treenode.cpp adapt_am_treenode.h segmenter.cpp segmenter.h \
	articulatorystream.cpp

shout_adapt_am_LDFLAGS = $(all_libraries)
shout_adapt_am_SOURCES = adapt_am.cpp adapt_am_main.cpp adapt_am_treenode.cpp \
	languagemodel_segmenter.cpp languagemodel_segmenter.h segmenter.cpp segmenter.h standard.h \
	train_segmenter.cpp train_segmenter.h
shout_prepare_adapt_LDFLAGS = $(all_libraries)
shout_prepare_adapt_SOURCES = adapt_am_treenode.cpp adapt_am_treenode.h debug.h standard.h shoutprepareadapt.cpp
shout_cluster_LDFLAGS = $(all_libraries)
shout_cluster_SOURCES = standard.h shout_cluster.cpp adapt_am_treenode.cpp \
	adapt_am_treenode.h adapt_segmenter.cpp adapt_segmenter.h segmenter.cpp segmenter.h \
	train_segmenter.cpp train_segmenter.h languagemodel_segmenter.cpp languagemodel_segmenter.h
lib_LIBRARIES =  libshout_basics.a
libshout_basics_a_SOURCES = standard.h gaussian.cc gaussian.h languagemodel.cpp 	languagemodel.h lexicaltree.cpp lexicaltree.h mixgaussian.cc mixgaussian.h 	phonemodel.cpp phonemodel.h vector.cc vector.h hash.cpp hash.h phonefilereader.cpp 	featurepool.cpp featureextraction.cpp featureextraction.h FFTReal.cpp FFTReal.h 	multimixgaussian.cpp shout_misc.cpp shout-misc.h shoutconfig.cpp shoutconfig.h 	stringlookup.cpp trainphonemodel.cpp trainphonemodel.h
shout_cluster_LDADD = $(top_builddir)/src/libshout_basics.a
shout_prepare_adapt_LDADD = $(top_builddir)/src/libshout_basics.a
shout_adapt_am_LDADD = $(top_builddir)/src/libshout_basics.a
shout_LDADD = $(top_builddir)/src/libshout_basics.a
shout_lm2bin_LDADD = $(top_builddir)/src/libshout_basics.a
shout_dct2lextree_LDADD = $(top_builddir)/src/libshout_basics.a
shout_maketrainset_LDADD = $(top_builddir)/src/libshout_basics.a

shout_merge_am_SOURCES = shoutmergeam.cpp
shout_merge_am_LDADD = $(top_builddir)/src/libshout_basics.a

shout_update_version_LDADD = $(top_builddir)/src/libshout_basics.a
shout_update_version_SOURCES = standard.h shout_updateversion.cpp

shout_segment_LDADD = $(top_builddir)/src/libshout_basics.a
shout_segment_SOURCES = shoutsegment.cpp shoutsegment.h segmenter.cpp segmenter.h
shout_normalize_am_SOURCES = normalizeam.cpp
shout_normalize_am_LDADD = $(top_builddir)/src/libshout_basics.a
shout_vtln_LDADD = $(top_builddir)/src/libshout_basics.a
shout_vtln_SOURCES = shout_vtln.cpp shout_vtln_main.cpp standard.h





shout_train_model_LDADD = $(top_builddir)/src/libshout_basics.a
shout_train_model_SOURCES = shouttrainmodel.cpp
shout_train_finish_SOURCES = shouttrainfinish.cpp
shout_train_finish_LDADD = $(top_builddir)/src/libshout_basics.a
shout_merge_trainset_LDADD = $(top_builddir)/src/libshout_basics.a
shout_merge_trainset_SOURCES = mergetrainset.cpp shout_maketrainset.cpp \
	shout_maketrainset.h
shout_train_finish_sat_SOURCES = shouttrainfinishsat.cpp
shout_train_finish_sat_LDADD = $(top_builddir)/src/libshout_basics.a
shout_online_SOURCES = shout_vtln.cpp shoutonline.cpp
shout_online_LDADD = $(top_builddir)/src/libshout_basics.a
shout_train_mmi_LDADD = $(top_builddir)/src/libshout_basics.a
shout_train_mmi_SOURCES = shouttrainmmi.cpp
shout_spkrec_LDADD = $(top_builddir)/src/libshout_basics.a
shout_spkrec_SOURCES = speakerrecognition.cpp
shout_preprocess_SOURCES = shout_preprocess.cpp shout_preprocess.h
shout_preprocess_LDADD = $(top_builddir)/src/libshout_basics.a
shout_spkrec_stats_LDADD = $(top_builddir)/src/libshout_basics.a
shout_spkrec_stats_SOURCES = spkrecstats.cpp spkrecstats.h
shout_avgEnergy_LDADD = $(top_builddir)/src/libshout_basics.a
shout_avgEnergy_SOURCES = avgenergy.cpp avgenergy.h
