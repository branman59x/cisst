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

void robLSPBTest::Log(const robLSPB & trajectory,
                      const std::string & name)
{
    const double startTime = trajectory.StartTime();
    const double duration = trajectory.Duration();
    const double extraPlotTime = 2.0;
    const double plotTime = extraPlotTime + duration + extraPlotTime;
    const size_t nbSteps = 2000;
    const double timeStep = plotTime / nbSteps;

    vctDoubleVec previousPosition(mDimension);
    vctDoubleVec previousVelocity(mDimension);

    // estimate velocities based on positions reported by LSPB
    vctDoubleVec deltaPosition(mDimension);
    deltaPosition.SetAll(0.0);
    // estimate accelerations based on velocities reported by LSPB
    vctDoubleVec deltaVelocity(mDimension);
    deltaVelocity.SetAll(0.0);

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
              << ','
              << cmnData<vctDoubleVec>::SerializeDescription(deltaPosition, ',', "deltaPosition")
              << ','
              << cmnData<vctDoubleVec>::SerializeDescription(deltaVelocity, ',', "deltaVelocity")
              << std::endl;

    for (size_t i = 0; i < nbSteps; ++i) {
        double now = (startTime - extraPlotTime) + i * timeStep;

        trajectory.Evaluate(now , mPosition, mVelocity, mAcceleration);

        if (i > 0) {
            deltaPosition.DifferenceOf(mPosition, previousPosition);
            deltaPosition.Divide(timeStep);
            deltaVelocity.DifferenceOf(mVelocity, previousVelocity);
            deltaVelocity.Divide(timeStep);
            // log to csv file
            cmnData<double>::SerializeText(now, log);
            log << ',';
            cmnData<vctDoubleVec>::SerializeText(mPosition, log);
            log << ',';
            cmnData<vctDoubleVec>::SerializeText(mVelocity, log);
            log << ',';
            cmnData<vctDoubleVec>::SerializeText(mAcceleration, log);
            log << ',';
            cmnData<vctDoubleVec>::SerializeText(deltaPosition, log);
            log << ',';
            cmnData<vctDoubleVec>::SerializeText(deltaVelocity, log);
            log << std::endl;
        }
        previousPosition.Assign(mPosition);
        previousVelocity.Assign(mVelocity);
    }

    log.close();
    logHeader.close();
}


void robLSPBTest::TestContinuity(const robLSPB & trajectory)
{
    const double startTime = trajectory.StartTime();
    const double duration = trajectory.Duration();
    const size_t nbSteps = 2000;
    const double timeStep = duration / nbSteps;

    vctDoubleVec previousPosition(mStart);
    vctDoubleVec previousVelocity(mInitialVelocity);
    vctDoubleVec previousAcceleration(mDimension);

    // max velocity might be higher because of initial velocity, take max
    vctDoubleVec maxVelocity(mDimension);
    maxVelocity.AbsOf(mInitialVelocity);
    maxVelocity.ElementwiseMax(mMaxVelocity);

    // std::cerr << "max velocities " << maxVelocity << std::endl;
    vctDoubleVec deltaPosition(mDimension);
    vctDoubleVec avgVelocity(mDimension);
    vctDoubleVec deltaVelocity(mDimension);

    for (size_t j = 0; j < nbSteps; j++) {
        double now = startTime + j * timeStep;

        trajectory.Evaluate(now , mPosition, mVelocity, mAcceleration);


        for (size_t i = 0; i < mDimension;++i) {
            // basic tests
            CPPUNIT_ASSERT(mVelocity[i] <= maxVelocity[i]);
            CPPUNIT_ASSERT(mVelocity[i] >= -maxVelocity[i]);
            CPPUNIT_ASSERT(mAcceleration[i] <= mMaxAcceleration[i]);
            CPPUNIT_ASSERT(mAcceleration[i] >= -mMaxAcceleration[i]);

            // compare previous and current to find derivatives
            deltaPosition[i] = (mPosition[i] - previousPosition[i]) / timeStep;
#if 0
            std::cerr << "previous " << previousPosition[i] << std::endl
                      << "current  " << mPosition[i] << std::endl
                      << "delta    " << deltaPosition[i] << std::endl
                      << "maxVelocity " << maxVelocity[i] << std::endl;
#endif
            CPPUNIT_ASSERT(deltaPosition[i] <= (maxVelocity[i] * 1.01));
            CPPUNIT_ASSERT(deltaPosition[i] >= (-maxVelocity[i] * 1.01));
            avgVelocity[i] = (mVelocity[i] + previousVelocity[i]) / 2.0;
            deltaVelocity[i] = (mVelocity[i] - previousVelocity[i]) / timeStep;

            // some tests can't be performed on first iteration because previous state is not well defined
            if (now > (startTime + timeStep)) {
#if 0
                std::cout << "previous " << previousVelocity[i] << std::endl
                          << "current  " << mVelocity[i] << std::endl
                          << "average  " << avgVelocity[i] << std::endl /* this is a C comment */
                          << "prevPos  " << previousPosition[i] << std::endl
                          << "currPos  " << mPosition[i] << std::endl
                          << "deltaPos " << deltaPosition[i] << std::endl
                          << "now      " << now << std::endl
                          << "deltaVel " << deltaVelocity[i]<<std::endl
                          << "accel    " << mAcceleration[i] << std::endl
                          << "index    " << i << std::endl;
#endif
                CPPUNIT_ASSERT((deltaPosition[i] - avgVelocity[i]) <= 0.01);
                CPPUNIT_ASSERT((deltaPosition[i] - avgVelocity[i]) >= -0.01);

                // only test when acceleration is constant
                if (previousAcceleration[i] == mAcceleration[i]) {
                    CPPUNIT_ASSERT((deltaVelocity[i] - mAcceleration[i]) <= (0.01 * mMaxAcceleration[i]));
                    CPPUNIT_ASSERT((deltaVelocity[i] - mAcceleration[i]) >= -(0.01 * mMaxAcceleration[i]));
                }

                CPPUNIT_ASSERT(deltaVelocity[i] <= (mMaxAcceleration[i] * 1.01));
                CPPUNIT_ASSERT(deltaVelocity[i] >= (-mMaxAcceleration[i] * 1.01));
            }

            previousPosition[i] = mPosition[i];
            previousVelocity[i] = mVelocity[i];
            previousAcceleration[i] = mAcceleration[i];
        }
    }
}

void robLSPBTest::PositiveViZero(void)
{
    SetDimension(1);
    mStart[0] = 2;
    mFinish[0] = 10;
    mMaxVelocity[0] = 200.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = 0.0;
    const double startTime = 2.0;

    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "PositiveViZero");
    TestContinuity(trajectory);
}

void robLSPBTest::PositiveViZeroPlateau(void)
{
    SetDimension(1);
    mStart[0] = 2;
    mFinish[0] = 10;
    mMaxVelocity[0] = 2.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = 0.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "PositiveViZeroPlateau");
    TestContinuity(trajectory);
}

void robLSPBTest::PositiveViPositive(void)
{
    SetDimension(1);
    mStart[0] = 2;
    mFinish[0] = 10;
    mMaxVelocity[0] = 200.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = 1.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "PositiveViPositive");
    TestContinuity(trajectory);
}

void robLSPBTest::PositiveViPositivePlateau(void)
{
    SetDimension(1);
    mStart[0] = 2;
    mFinish[0] = 10;
    mMaxVelocity[0] = 2.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = 1.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "PositiveViPositivePlateau");
    TestContinuity(trajectory);
}

void robLSPBTest::PositiveViAbove(void)
{
    SetDimension(1);
    mStart[0] = 0.0;
    mFinish[0] = 4.0;
    mMaxVelocity[0] = 2.0;
    mMaxAcceleration[0] = 1.0;
    mInitialVelocity[0] = 2.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "PositiveViAbove");
    TestContinuity(trajectory);
}

void robLSPBTest::PositiveViAbovePlateau(void)
{
    SetDimension(1);
    mStart[0] = 2;
    mFinish[0] = 10;
    mMaxVelocity[0] = 3.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = 4.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "PositiveViAbovePlateau");
    TestContinuity(trajectory);
}

void robLSPBTest::PositiveViPositiveOvershot(void)
{
    SetDimension(1);
    mStart[0] = 3.0;
    mFinish[0] = 4.0;
    mMaxVelocity[0] = 8.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = 5.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "PositiveViPositiveOvershot");
    TestContinuity(trajectory);
}

void robLSPBTest::PositiveViPositiveOvershotPlateau(void)
{
    SetDimension(1);
    mStart[0] = 3.0;
    mFinish[0] = 5.0;
    mMaxVelocity[0] = 2.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = 8.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "PositiveViPositiveOvershotPlateau");
    TestContinuity(trajectory);
}




void robLSPBTest::PositiveViNegativeOvershot(void)
{
    SetDimension(1);
    mStart[0] = 2.0;
    mFinish[0] = 10.0;
    mMaxVelocity[0] = 200.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = -1.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "PositiveViNegativeOvershot");
    TestContinuity(trajectory);
}

void robLSPBTest::PositiveViNegativeOvershotPlateau(void)
{
    SetDimension(1);
    mStart[0] = 2.0;
    mFinish[0] = 10.0;
    mMaxVelocity[0] = 2.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = -1.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "PositiveViNegativeOvershotPlateau");
    TestContinuity(trajectory);
}

void robLSPBTest::NegativeViZero(void)
{
    SetDimension(1);
    mStart[0] = 10.0;
    mFinish[0] = 2.0;
    mMaxVelocity[0] = 200.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = 0.0;
    const double startTime = 2.0;

    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "NegativeViZero");
    TestContinuity(trajectory);
}

void robLSPBTest::NegativeViZeroPlateau(void)
{
    SetDimension(1);
    mStart[0] = 10.0;
    mFinish[0] = 2.0;
    mMaxVelocity[0] = 2.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = 0.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "NegativeViZeroPlateau");
    TestContinuity(trajectory);
}

void robLSPBTest::NegativeViNegative(void)
{
    SetDimension(1);
    mStart[0] = 10.0;
    mFinish[0] = 2.0;
    mMaxVelocity[0] = 200.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = -1.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "NegativeViNegative");
    TestContinuity(trajectory);
}

void robLSPBTest::NegativeViNegativePlateau(void)
{
    SetDimension(1);
    mStart[0] = 10.0;
    mFinish[0] = 2.0;
    mMaxVelocity[0] = 2.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = -1.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "NegativeViNegativePlateau");
    TestContinuity(trajectory);
}

void robLSPBTest::NegativeViBelow(void)
{
    SetDimension(1);
    mStart[0] = 4.0;
    mFinish[0] = 0.0;
    mMaxVelocity[0] = 2.0;
    mMaxAcceleration[0] = 1.0;
    mInitialVelocity[0] = -2.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "NegativeViBelow");
    TestContinuity(trajectory);
}

void robLSPBTest::NegativeViBelowPlateau(void)
{
    SetDimension(1);
    mStart[0] = 10;
    mFinish[0] = 2;
    mMaxVelocity[0] = 3.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = -4.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "NegativeViBelowPlateau");
    TestContinuity(trajectory);
}

void robLSPBTest::NegativeViNegativeOvershot(void)
{
    SetDimension(1);
    mStart[0] = 4.0;
    mFinish[0] = 3.0;
    mMaxVelocity[0] = 8.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = -5.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "NegativeViNegativeOvershot");
    TestContinuity(trajectory);
}

void robLSPBTest::NegativeViNegativeOvershotPlateau(void)
{
    SetDimension(1);
    mStart[0] = 5.0;
    mFinish[0] = 3.0;
    mMaxVelocity[0] = 2.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = -8.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "NegativeViNegativeOvershotPlateau");
    TestContinuity(trajectory);
}

void robLSPBTest::NegativeViPositiveOvershot(void)
{
    SetDimension(1);
    mStart[0] = 10.0;
    mFinish[0] = 2.0;
    mMaxVelocity[0] = 200.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = 1.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "NegativeViPositiveOvershot");
    TestContinuity(trajectory);
}

void robLSPBTest::NegativeViPositiveOvershotPlateau(void)
{
    SetDimension(1);
    mStart[0] = 10.0;
    mFinish[0] = 2.0;
    mMaxVelocity[0] = 2.0;
    mMaxAcceleration[0] = 2.0;
    mInitialVelocity[0] = 1.0;

    const double startTime = 2.0;
    robLSPB trajectory;
    trajectory.Set(mStart, mFinish,
                   mMaxVelocity, mMaxAcceleration,mInitialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    Log(trajectory, "NegativeViPositiveOvershotPlateau");
    TestContinuity(trajectory);
}

void robLSPBTest::MultipleJoints(void)
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

    Log(trajectory, "MultipleJoints");
    TestContinuity(trajectory);
}
