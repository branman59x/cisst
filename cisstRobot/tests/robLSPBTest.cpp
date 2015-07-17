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

void robLSPBTest::SetDimension(const size_t dim)
{
    mDimension = dim;
    mStart.SetSize(mDimension);
    mFinish.SetSize(mDimension);
    mMaxVelocity.SetSize(mDimension);
    mMaxAcceleration.SetSize(mDimension);
    mPosition.SetSize(mDimension);
    mVelocity.SetSize(mDimension);
    mAcceleration.SetSize(mDimension);
    mInitialVelocity.SetSize(mDimension);
}

void robLSPBTest::LogAndTestContinuity(robLSPB & trajectory,
                                       const std::string & name)
{
    const double startTime = trajectory.StartTime();
    const double duration = trajectory.Duration();
    const double extraPlotTime = 2.0;
    const double plotTime = extraPlotTime + duration + extraPlotTime;
    const size_t nbSteps = 200000;
    const double step = plotTime / nbSteps;

    // First log all data without tests
    std::ofstream log, logHeader;
    const std::string logName = "robLSPB" + name + ".txt";
    const std::string logHeaderName = "robLSPB" + name + "-header.txt";
    log.open(logName.c_str());
    logHeader.open(logHeaderName.c_str());

    // header for logs
    logHeader << cmnData<double>::SerializeDescription(duration, ',', "time")
              << ','
              << cmnData<vctDoubleVec>::SerializeDescription(mPosition, ',', "position")
              << ','
              << cmnData<vctDoubleVec>::SerializeDescription(mVelocity, ',', "velocity")
              << ','
              << cmnData<vctDoubleVec>::SerializeDescription(mAcceleration, ',', "acceleration")
              << std::endl;

    for (size_t i = 0; i < nbSteps; ++i) {
        double now = (startTime - extraPlotTime) + i * step;
        trajectory.Evaluate(now , mPosition, mVelocity, mAcceleration);
        // log to csv file
        cmnData<double>::SerializeText(now, log);
        log << ',';
        cmnData<vctDoubleVec>::SerializeText(mPosition, log);
        log << ',';
        cmnData<vctDoubleVec>::SerializeText(mVelocity, log);
        log << ',';
        cmnData<vctDoubleVec>::SerializeText(mAcceleration, log);
        log << std::endl;
    }

    log.close();
    logHeader.close();
    std::cout<<"NOW TESTING CONTINUITY"<< std::endl;
    // Then test continuity
    vctDoubleVec previousPosition(mStart);
    vctDoubleVec previousVelocity(mInitialVelocity);
    vctDoubleVec previousTime(mDimension); //= startTime - extraPlotTime - step;
    previousTime.SetAll(startTime - extraPlotTime - step);
    vctDoubleVec maxVelocity(mDimension);
    maxVelocity.AbsOf(mInitialVelocity);
    maxVelocity.ElementwiseMax(mMaxVelocity);

    std::cerr << "max velocities " << maxVelocity << std::endl;

    for (size_t j = 0; j < nbSteps; j++) {
        double now = (startTime - extraPlotTime) + j * step;
        trajectory.Evaluate(now , mPosition, mVelocity, mAcceleration);
        std::cerr << now << std::endl;

        for (size_t i = 0; i < mDimension;++i){
            // basic tests
            CPPUNIT_ASSERT(mVelocity[i] <= maxVelocity[i]);
            CPPUNIT_ASSERT(mVelocity[i] >= -maxVelocity[i]);
            CPPUNIT_ASSERT(mAcceleration[i] <= mMaxAcceleration[i]);
            CPPUNIT_ASSERT(mAcceleration[i] >= -mMaxAcceleration[i]);

            // compare previous and current to find derivatives
            vctDoubleVec deltaTime(mDimension);
            deltaTime[i] = now - previousTime[i];
            vctDoubleVec deltaPosition(mDimension);
            deltaPosition[i] = mPosition[i] - previousPosition[i];
            std::cerr << "previous " << previousPosition[i] << std::endl
                      << "current  " << mPosition[i] << std::endl
                      << "delta    " << deltaPosition[i] << std::endl
                      << "maxVelocity " << maxVelocity[i] << std::endl;
            deltaPosition[i] /= deltaTime[i];

            CPPUNIT_ASSERT(deltaPosition[i] <= (maxVelocity[i] * 1.1));
            CPPUNIT_ASSERT(deltaPosition[i] >= (-maxVelocity[i] * 1.1));

            vctDoubleVec deltaVelocity(mDimension);
            deltaVelocity[i] = mVelocity[i] - previousVelocity[i];
            deltaVelocity[i] /= deltaTime[i];
            CPPUNIT_ASSERT(deltaVelocity[i] <= (mMaxAcceleration[i] * 1.1));
            CPPUNIT_ASSERT(deltaVelocity[i] >= (-mMaxAcceleration[i] * 1.1));

            previousTime[i] = now;
            previousPosition[i] = mPosition[i];
            previousVelocity[i] = mVelocity[i];
        }
    }
}

void robLSPBTest::Test1(void)
{
    SetDimension(1);
    mStart[0] = 0.0;
    mFinish[0] = 10.0;
    mMaxVelocity[0] = 10.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = 5.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    LogAndTestContinuity(trajectory, "Test1");
}

void robLSPBTest::Test2(void)
{
    SetDimension(1);
//    mStart[0] = 1.0;
//    mFinish[0] = 30.0;
//    mMaxVelocity[0] = 10.0;
//    mMaxAcceleration[0] = 50.0;
//    mInitialVelocity[0] = 29.0;
    mStart[0] = 0;
    mFinish[0] = 10;
    mMaxVelocity[0] = 4;
    mMaxAcceleration[0] = 1;
    mInitialVelocity[0] = 0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    LogAndTestContinuity(trajectory, "Test2");
}


void robLSPBTest::Test3(void)
{
    SetDimension(6);
    mStart.Assign(          0.0, 2.0, 3.0, 0.0, 2.0, 3.0);
    mFinish.Assign(        10.0, 12.0, 5.0, 10.0, 12.0, 15.0);
    mMaxVelocity.Assign(    4.0, 2.0, 5.0, 4.0, 2.0, 5.0);
    mMaxAcceleration.Assign(1.0, 1.0, 2.0, 1.0, 1.0, 2.0);
    mInitialVelocity.Assign(0.0, 1.0, 4.0, 5.0, -1.0, -6.0);

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    LogAndTestContinuity(trajectory, "Test3");
}

void robLSPBTest::Test4(void)
{
    SetDimension(2);
    mStart.Assign(0,2);
    mFinish.Assign(10,12);
    mMaxVelocity.Assign(4,2);
    mMaxAcceleration.Assign(1,1);
    mInitialVelocity.Assign(0,1);

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    LogAndTestContinuity(trajectory, "Test4");
}

void robLSPBTest::Test5(void)
{
    SetDimension(1);
    mStart[0] = 0;
    mFinish[0] = 10;
    mMaxVelocity[0] = 4;
    mMaxAcceleration[0] = 1;
    mInitialVelocity[0] = 0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    LogAndTestContinuity(trajectory, "Test5");
}

void robLSPBTest::Test6(void)
{
    SetDimension(1);
    mStart[0] = 2;
    mFinish[0] = 12;
    mMaxVelocity[0] = 2;
    mMaxAcceleration[0] = 1;
    mInitialVelocity[0] = 1;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    LogAndTestContinuity(trajectory, "Test6");
}
