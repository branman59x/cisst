/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-    */
/* ex: set filetype=cpp softtabstop=4 shiftwidth=4 tabstop=4 cindent expandtab: */

/*
  Author(s):  Anton Deguet, Brandon Cohen
  Created on: 2015-07-14
  
  (C) Copyright 2015 Johns Hopkins University (JHU), All Rights Reserved.

--- begin cisst license - do not edit ---

This software is provided "as is" under an open source license, with
no warranty.  The complete license can be found in license.txt and
http://www.cisst.org/cisst/license.txt.

--- end cisst license ---
*/


#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "cisstRobot/robLSPB.h"

class robLSPBTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(robLSPBTest);
    {
        CPPUNIT_TEST(PositiveViZero);
        CPPUNIT_TEST(PositiveViZeroPlateau);
        CPPUNIT_TEST(PositiveViPositive);
        CPPUNIT_TEST(PositiveViPositivePlateau);
        CPPUNIT_TEST(PositiveViAbove);
        CPPUNIT_TEST(PositiveViAbovePlateau);
        CPPUNIT_TEST(PositiveViPositiveOvershot);
        CPPUNIT_TEST(PositiveViPositiveOvershotPlateau);
        CPPUNIT_TEST(PositiveViNegativeOvershot);
        CPPUNIT_TEST(PositiveViNegativeOvershotPlateau);

        CPPUNIT_TEST(NegativeViZero);
        CPPUNIT_TEST(NegativeViZeroPlateau);
        CPPUNIT_TEST(NegativeViNegative);
        CPPUNIT_TEST(NegativeViNegativePlateau);
        CPPUNIT_TEST(NegativeViBelow);
        CPPUNIT_TEST(NegativeViBelowPlateau);
        CPPUNIT_TEST(NegativeViNegativeOvershot);
        CPPUNIT_TEST(NegativeViNegativeOvershotPlateau);
        CPPUNIT_TEST(NegativeViPositiveOvershot);
        CPPUNIT_TEST(NegativeViPositiveOvershotPlateau);

        CPPUNIT_TEST(MultipleJoints);
    }
    CPPUNIT_TEST_SUITE_END();
    
 public:
    void setUp(void) {
    }
    
    void tearDown(void) {
    }
    
    void PositiveViZero(void);
    void PositiveViZeroPlateau(void);
    void PositiveViPositive(void);
    void PositiveViPositivePlateau(void);
    void PositiveViAbove(void);
    void PositiveViAbovePlateau(void);
    void PositiveViPositiveOvershot(void);
    void PositiveViPositiveOvershotPlateau(void);
    void PositiveViNegativeOvershot(void);
    void PositiveViNegativeOvershotPlateau(void);

    void NegativeViZero(void);
    void NegativeViZeroPlateau(void);
    void NegativeViNegative(void);
    void NegativeViNegativePlateau(void);
    void NegativeViBelow(void);
    void NegativeViBelowPlateau(void);
    void NegativeViNegativeOvershot(void);
    void NegativeViNegativeOvershotPlateau(void);
    void NegativeViPositiveOvershot(void);
    void NegativeViPositiveOvershotPlateau(void);

    void MultipleJoints(void);

 protected:

    /*! Set dimension of all data members */
    void SetDimension(const size_t mDimension);

    /*! Use the trajectory to create a csv file with position, velocity and acceleration. */
    void Log(const robLSPB & trajectory,
             const std::string & name);

    /*! Basic tests of continuity and bounds on position, velocity and acceleration */
    void TestContinuity(const robLSPB & trajectory);

    size_t mDimension;
    vctDoubleVec
        mStart,
        mFinish,
        mMaxVelocity,
        mMaxAcceleration,
        mPosition,
        mVelocity,
        mAcceleration,
        mInitialVelocity;
};

CPPUNIT_TEST_SUITE_REGISTRATION(robLSPBTest);
