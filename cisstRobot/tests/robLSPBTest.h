/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-    */
/* ex: set filetype=cpp softtabstop=4 shiftwidth=4 tabstop=4 cindent expandtab: */

/*
  Author(s):  Anton Deguet
  Created on: 2014-07-14
  
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
        CPPUNIT_TEST(Test1);
        CPPUNIT_TEST(Test2);
        CPPUNIT_TEST(Test3);
        CPPUNIT_TEST(Test4);
        CPPUNIT_TEST(Test5);
        CPPUNIT_TEST(Test6);
    }
    CPPUNIT_TEST_SUITE_END();
    
 public:
    void setUp(void) {
    }
    
    void tearDown(void) {
    }
    
    void Test1(void);
    void Test2(void);
    void Test3(void);
    void Test4(void);
    void Test5(void);
    void Test6(void);

 protected:

    void SetDimension(const size_t mDimension);
    void LogAndTestContinuity(robLSPB & trajectory,
                              const std::string & name);

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

