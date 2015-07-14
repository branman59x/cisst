/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-    */
/* ex: set filetype=cpp softtabstop=4 shiftwidth=4 tabstop=4 cindent expandtab: */

/*
  Author(s):  Anton Deguet
  Created on: 2015-07-14
  
  (C) Copyright 2015 Johns Hopkins University (JHU), All Rights Reserved.

--- begin cisst license - do not edit ---

This software is provided "as is" under an open source license, with
no warranty.  The complete license can be found in license.txt and
http://www.cisst.org/cisst/license.txt.

--- end cisst license ---
*/


#include "robLSPBTest.h"


void robLSPBTest::Test1(void) {
    const size_t dimension = 1;
    vctDoubleVec
        start,
        finish,
        maxVelocity,
        maxAcceleration,
        position,
        velocity,
        acceleration,
        initialVelocity;
    start.SetSize(dimension);
    finish.SetSize(dimension);
    maxVelocity.SetSize(dimension);
    maxAcceleration.SetSize(dimension);
    position.SetSize(dimension);
    velocity.SetSize(dimension);
    acceleration.SetSize(dimension);
    initialVelocity.SetSize(dimension);

    start[0] = 1.0;
    finish[0] = 30.0;
    maxVelocity[0] = 10.0;
    maxAcceleration[0] = 50.0;
    initialVelocity[0] = 0.0;

    const double startTime = 2.0;

    robLSPB trajectory;
    trajectory.Set(start, finish,
                   maxVelocity, maxAcceleration,initialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    const double duration = trajectory.Duration();
    const double extraPlotTime = 2.0;
    const double plotTime = extraPlotTime + duration + extraPlotTime;
    const size_t nbSteps = 2000;
    const double step = plotTime / nbSteps;

    std::cout << "duration: " << duration << std::endl;

    std::ofstream log, logHeader;
    const char * logName = "robLSPBTest1.txt";
    const char * logHeaderName = "robLSPBTest1-header.txt";
    log.open(logName);
    logHeader.open(logHeaderName);
    // header for logs
    logHeader << cmnData<double>::SerializeDescription(duration, ',', "time")
              << ','
              << cmnData<vctDoubleVec>::SerializeDescription(position, ',', "position")
              << ','
              << cmnData<vctDoubleVec>::SerializeDescription(velocity, ',', "velocity")
              << ','
              << cmnData<vctDoubleVec>::SerializeDescription(acceleration, ',', "acceleration")
              << std::endl;

    vctDoubleVec previousPosition(start);
    vctDoubleVec previousVelocity(initialVelocity);
    double previousTime = startTime - extraPlotTime - step;

    for (size_t i = 0; i < nbSteps; ++i) {
        double now = (startTime - extraPlotTime) + i * step;
        trajectory.Evaluate(now , position, velocity, acceleration);

        // tests
        CPPUNIT_ASSERT(velocity.LesserOrEqual(maxVelocity));
        CPPUNIT_ASSERT(velocity.GreaterOrEqual(-maxVelocity));
        CPPUNIT_ASSERT(acceleration.LesserOrEqual(maxAcceleration));
        CPPUNIT_ASSERT(acceleration.GreaterOrEqual(-maxAcceleration));

        // compare previous and current to find derivatives
        double deltaTime = now - previousTime;

        vctDoubleVec deltaPosition(dimension);
        deltaPosition.DifferenceOf(position, previousPosition);
        deltaPosition.Divide(deltaTime);
        std::cerr << deltaPosition[0] - 10.0 << std::endl;
        CPPUNIT_ASSERT(deltaPosition.LesserOrEqual(maxVelocity * 1.001));
        CPPUNIT_ASSERT(deltaPosition.GreaterOrEqual(-maxVelocity * 1.001));

        vctDoubleVec deltaVelocity(dimension);
        deltaPosition.DifferenceOf(velocity, previousVelocity);
        deltaPosition.Divide(deltaTime);
        CPPUNIT_ASSERT(deltaVelocity.LesserOrEqual(maxAcceleration * 1.001));
        CPPUNIT_ASSERT(deltaVelocity.GreaterOrEqual(-maxAcceleration * 1.001));

        previousTime = now;
        previousPosition.Assign(position);
        previousVelocity.Assign(velocity);

        // log to csv file
        cmnData<double>::SerializeText(now, log);
        log << ',';
        cmnData<vctDoubleVec>::SerializeText(position, log);
        log << ',';
        cmnData<vctDoubleVec>::SerializeText(velocity, log);
        log << ',';
        cmnData<vctDoubleVec>::SerializeText(acceleration, log);
        log << std::endl;
    }

    log.close();
    logHeader.close();
}

void robLSPBTest::Test2(void) {
    const size_t dimension = 1;
    vctDoubleVec
        start,
        finish,
        maxVelocity,
        maxAcceleration,
        position,
        velocity,
        acceleration,
        initialVelocity;
    start.SetSize(dimension);
    finish.SetSize(dimension);
    maxVelocity.SetSize(dimension);
    maxAcceleration.SetSize(dimension);
    position.SetSize(dimension);
    velocity.SetSize(dimension);
    acceleration.SetSize(dimension);
    initialVelocity.SetSize(dimension);

    start[0] = 1.0;
    finish[0] = 30.0;
    maxVelocity[0] = 10.0;
    maxAcceleration[0] = 50.0;
    initialVelocity[0] = 29.0;

    const double startTime = 2.0;

    robLSPB trajectory;
    trajectory.Set(start, finish,
                   maxVelocity, maxAcceleration,initialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    const double duration = trajectory.Duration();
    const double extraPlotTime = 2.0;
    const double plotTime = extraPlotTime + duration + extraPlotTime;
    const size_t nbSteps = 2000;
    const double step = plotTime / nbSteps;

    std::cout << "duration: " << duration << std::endl;

    std::ofstream log, logHeader;
    const char * logName = "robLSPBTest2.txt";
    const char * logHeaderName = "robLSPBTest2-header.txt";
    log.open(logName);
    logHeader.open(logHeaderName);
    // header for logs
    logHeader << cmnData<double>::SerializeDescription(duration, ',', "time")
              << ','
              << cmnData<vctDoubleVec>::SerializeDescription(position, ',', "position")
              << ','
              << cmnData<vctDoubleVec>::SerializeDescription(velocity, ',', "velocity")
              << ','
              << cmnData<vctDoubleVec>::SerializeDescription(acceleration, ',', "acceleration")
              << std::endl;

    vctDoubleVec previousPosition(start);
    vctDoubleVec previousVelocity(initialVelocity);
    double previousTime = startTime - extraPlotTime - step;

    for (size_t i = 0; i < nbSteps; ++i) {
        double now = (startTime - extraPlotTime) + i * step;
        trajectory.Evaluate(now , position, velocity, acceleration);

        // tests
        //CPPUNIT_ASSERT(velocity.LesserOrEqual(maxVelocity));
        //CPPUNIT_ASSERT(velocity.GreaterOrEqual(-maxVelocity));
        CPPUNIT_ASSERT(acceleration.LesserOrEqual(maxAcceleration));
        CPPUNIT_ASSERT(acceleration.GreaterOrEqual(-maxAcceleration));

        // compare previous and current to find derivatives
        double deltaTime = now - previousTime;

        vctDoubleVec deltaPosition(dimension);
        deltaPosition.DifferenceOf(position, previousPosition);
        deltaPosition.Divide(deltaTime);
        std::cerr << deltaPosition[0] << std::endl;
        vctDoubleVec totallyMaxVelocity(maxVelocity);
        totallyMaxVelocity.ElementwiseMax(initialVelocity);
        CPPUNIT_ASSERT(deltaPosition.LesserOrEqual(totallyMaxVelocity * 1.001));
        CPPUNIT_ASSERT(deltaPosition.GreaterOrEqual(-maxVelocity * 1.001));

        vctDoubleVec deltaVelocity(dimension);
        deltaPosition.DifferenceOf(velocity, previousVelocity);
        deltaPosition.Divide(deltaTime);
        CPPUNIT_ASSERT(deltaVelocity.LesserOrEqual(maxAcceleration * 1.001));
        CPPUNIT_ASSERT(deltaVelocity.GreaterOrEqual(-maxAcceleration * 1.001));
        previousTime = now;
        previousPosition.Assign(position);
        previousVelocity.Assign(velocity);

        // log to csv file
        cmnData<double>::SerializeText(now, log);
        log << ',';
        cmnData<vctDoubleVec>::SerializeText(position, log);
        log << ',';
        cmnData<vctDoubleVec>::SerializeText(velocity, log);
        log << ',';
        cmnData<vctDoubleVec>::SerializeText(acceleration, log);
        log << std::endl;
    }

    log.close();
    logHeader.close();
}
