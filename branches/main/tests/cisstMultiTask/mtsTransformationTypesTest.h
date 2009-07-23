/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-    */
/* ex: set filetype=cpp softtabstop=4 shiftwidth=4 tabstop=4 cindent expandtab: */

/*
  $Id: mtsTransformationTypesTest.h 456 2009-06-13 03:11:44Z adeguet1 $
  
  Author(s):  Anton Deguet
  Created on: 2009-04-29
  
  (C) Copyright 2009 Johns Hopkins University (JHU), All Rights
  Reserved.

--- begin cisst license - do not edit ---

This software is provided "as is" under an open source license, with
no warranty.  The complete license can be found in license.txt and
http://www.cisst.org/cisst/license.txt.

--- end cisst license ---
*/

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include <cisstMultiTask/mtsMatrix.h>
#include <cisstMultiTask/mtsTransformationTypes.h>

class mtsTransformationTypesTest: public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(mtsTransformationTypesTest);

    CPPUNIT_TEST(TestSerializeRawDoubleQuat);
    CPPUNIT_TEST(TestSerializeRawDoubleQuatRot3);
    CPPUNIT_TEST(TestSerializeRawDoubleAxAnRot3);
    CPPUNIT_TEST(TestSerializeRawDoubleRodRot3);
    CPPUNIT_TEST(TestSerializeRawDoubleMatRot3);
    CPPUNIT_TEST(TestSerializeRawDoubleQuatFrm3);
    CPPUNIT_TEST(TestSerializeRawDoubleMatFrm3);

    CPPUNIT_TEST(TestSerializeRawDoubleMat);

    CPPUNIT_TEST_SUITE_END();
    
public:
    void setUp(void) {}
    
    void tearDown(void) {}
    
    /*! Test the SerializeRaw method */
    template <class _elementType> void TestSerializeRaw(_elementType & initial);
    void TestSerializeRawDoubleQuat(void);
    void TestSerializeRawDoubleQuatRot3(void);
    void TestSerializeRawDoubleAxAnRot3(void);
    void TestSerializeRawDoubleRodRot3(void);
    void TestSerializeRawDoubleMatRot3(void);
    void TestSerializeRawDoubleQuatFrm3(void);
    void TestSerializeRawDoubleMatFrm3(void);

    void TestSerializeRawDoubleMat(void);
};


CPPUNIT_TEST_SUITE_REGISTRATION(mtsTransformationTypesTest);

