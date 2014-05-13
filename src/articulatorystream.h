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
#ifndef ARTICULATORYSTREAM_H
#define ARTICULATORYSTREAM_H

#include "phonemodel.h"
#include "vector.h"
#include "featurepool.h"

struct DecoderSettings;
struct WLRType;
struct WLRTracker;
struct PLRType;
struct TokenType;
struct LatticeNode;
struct WLRList;
struct LMLAGlobalListType;
struct LexicalNode;

class ArticulatoryStream
{
  public:
    double                   *probs;

  protected:
    PhoneModel              **modelsYes;
    PhoneModel              **modelsNo;
    int                       vectorSize;
    int                       numberOfModels;
    int                       numberOfPhoneModels;
    int                      *decisionMatrix;
    int                       artModelNumber[MAX_NUMBER_OF_DECISION_RULES];

#ifdef TEST_ARTICULATORY
    LexicalNode         *classifier;
    WLRType             *wlrStart;

    FeaturePool *pool;
    double                   *testRT;
    double                   *testRF;
    int                      confusion[40][40];
    int                      testRulesTT[36][40];
    int                      testRulesFT[36][40];
    int                      testRulesTF[36][40];
    int                      testRulesFF[36][40];

    double                   *testYes;
    double                   *testNo;
    double                   *testYesTot;
    double                   *testNoTot;
    int                      *testOK;
    int                      *testFP;
    int                      *testMISS;
    int                      *testAllOK;
    int                      *testAllFP;
    int                      *testAllMISS;

    PhoneModel              **phoneModels;
#endif

  public:
    ArticulatoryStream        (char *fileName);
   ~ArticulatoryStream        ();

    void processVector        (Vector *v);

#ifdef TEST_ARTICULATORY
    void setPool (FeaturePool *p, PhoneModel **pM);
    void testStream(int phoneID, int firstFrame, int lastFrame, DecoderSettings  *settings, int contextKey);
    void testPrint(FILE *outputFile);
    void testPrint2(FILE *outputFile);
#endif
};

#endif
