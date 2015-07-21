/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-    */
/* ex: set filetype=cpp softtabstop=4 shiftwidth=4 tabstop=4 cindent expandtab: */

/*
  Author(s):  Anton Deguet, Brandon Cohen
  Created on: 2014-10-27

(C) Copyright 2014 Johns Hopkins University (JHU), All Rights Reserved.

--- begin cisst license - do not edit ---
This software is provided "as is" under an open source license, with
no warranty.  The complete license can be found in license.txt and
http://www.cisst.org/cisst/license.txt.
--- end cisst license ---

*/

#include <cisstRobot/robLSPB.h>

robLSPB::robLSPB(void):
    mIsSet(false)
{}

robLSPB::robLSPB(const vctDoubleVec & start,
                 const vctDoubleVec & finish,
                 const vctDoubleVec & velocity,
                 const vctDoubleVec & acceleration,
                 const vctDoubleVec & initialVelocity,
                 const double startTime,
                 const CoordinationType coordination)
{
    Set(start, finish, velocity, acceleration,initialVelocity, startTime, coordination);
}

void robLSPB::Set(const vctDoubleVec & start,
                  const vctDoubleVec & finish,
                  const vctDoubleVec & velocity,
                  const vctDoubleVec & acceleration,
                  const vctDoubleVec & initialVelocity,
                  const double startTime,
                  const CoordinationType coordination)
{
    mIsSet = false;

    // sanity checks
    mDimension = start.size();
    if (finish.size() != mDimension) {
        cmnThrow("robLSPB::Set: finish point doesn't match start point dimension");
    }
    if (velocity.size() != mDimension) {
        cmnThrow("robLSPB::Set: velocity doesn't match start point dimension");
    }
    if (acceleration.size() != mDimension) {
        cmnThrow("robLSPB::Set: acceleration doesn't match start point dimension");
    }
    // store information and resize data members
    mCoordination = coordination;
    mStartTime = startTime;
    mStart.ForceAssign(start);
    mFinish.ForceAssign(finish);
    mVelocity.ForceAssign(velocity);
    mAcceleration.ForceAssign(acceleration);
    mAccelerationTime.SetSize(mDimension);
    mDecelerationTime.SetSize(mDimension);
    mFinishTime.SetSize(mDimension);
    mInitialVelocity.ForceAssign(initialVelocity);
    mInitialDirection.SetSize(mDimension);
    mInitialDirection.SetAll(1.0);
    mSecondAccelDistance.SetSize(mDimension);//Distance of Acceleration of a trajectory after an overshoot
    mSecondDistance.SetSize(mDimension);//Total distance of a trajectory after an overshoot
    mOvershotTime.SetSize(mDimension);//Time of Overshoot
    mOvershotTime.Zeros();
    mOvershotDistance.SetSize(mDimension);
    mAccelerationDistance.SetSize(mDimension);
    mDecelerationDistance.SetSize(mDimension);
    mAccelerationDistance.Zeros();
    mDecelerationDistance.Zeros();
    mDirection.SetSize(mDimension);
    mOvershotInitialVelocity.SetSize(mDimension);
    mOvershotDirection.SetSize(mDimension);
    mOvershotInitialVelocity.SetSize(mDimension);
    mOvershotStart.SetSize(mDimension);

    // compute trajectory parameters
    for (size_t i = 0;
         i < mDimension;
         ++i) {
        // compute direction
        double displacement = mFinish[i] - mStart[i];
        if (displacement < 0.0) {
            // sets variables so that the graph will flip as if the displacement was positive
            mDirection[i] = -1.0;
        } else {
            mDirection[i] = 1.0;
            if (displacement == 0 && mInitialVelocity[i] == 0) {
                mVelocity[i] = 0.0;
                mAcceleration[i] = 0.0;
            }
        }
        // compute time if distance != 0
        if (displacement != 0 || mInitialVelocity[i] != 0) {
            if (mVelocity[i] == 0.0) {
                cmnThrow("robLSPB::Set: velocity must be greater than zero");
            }
            if (mAcceleration[i] == 0.0) {
                cmnThrow("robLSPB::Set: acceleration must be greater than zero");
            }
            // look for overshoot, first case we are going in wrong direction
            if (mInitialVelocity[i] < 0.0) {
                mInitialDirection[i] = -1.0;
            }
            // tests if displacement is in the opposite direction of the velocity
            if (mInitialVelocity[i] * mDirection[i] < 0.0) {
                mOvershotTime[i] = -(mInitialVelocity[i] * mDirection[i] / mAcceleration[i]);
                mOvershotDistance[i] =
                        0.5 * mAcceleration[i] * mOvershotTime[i] * mOvershotTime[i] * mInitialDirection[i];
                mOvershotInitialVelocity[i] = 0.0;
                mOvershotDirection[i] = mDirection[i];
            } else {
                mOvershotTime[i] = fabs(mInitialVelocity[i] / mAcceleration[i]);
                mOvershotDistance[i] =
                        0.5 * mAcceleration[i] * mOvershotTime[i] * mOvershotTime[i] * mInitialDirection[i];
                // test if this is really an overshot
                if (fabs(mOvershotDistance[i]) <= fabs(mFinish[i] - mStart[i])) {
                    // test if we start with an initial velocity too high and we need to decelerate
                    if (fabs(mInitialVelocity[i]) > mVelocity[i]) {
                        mOvershotTime[i] = (fabs(mInitialVelocity[i]) - mVelocity[i]) / mAcceleration[i];
                        mOvershotDistance[i] = mInitialVelocity[i] * mOvershotTime[i]
                                - 0.5 * mAcceleration[i] * mOvershotTime[i] * mOvershotTime[i] * mInitialDirection[i];
                        mOvershotInitialVelocity[i] = mVelocity[i];
                        mOvershotDirection[i] = mDirection[i];
                    } else {
                        // really not an overshot
                        mOvershotTime[i] = 0.0;
                        mOvershotDistance[i] = 0.0;
                        mOvershotInitialVelocity[i] = mInitialVelocity[i];
                        mOvershotDirection[i] = mDirection[i];
                    }
                } else {
                    // we had an overshot, we stopped
                    mOvershotInitialVelocity[i] = 0.0;
                    // and then go backward
                    mOvershotDirection[i] = -mDirection[i];
                }
            }
            mOvershotStart[i] = mStart[i] + mOvershotDistance[i];
            displacement = mFinish[i] - (mStart[i] + mOvershotDistance[i]);
            // plan trajectory from the overshot point
            // calculates the max velocity that the velocity will hit (may be lower than predetermined)
            const double peakVelocity = sqrt((2.0 * mAcceleration[i] * fabs(displacement)
                                              + mOvershotInitialVelocity[i] * mOvershotInitialVelocity[i])
                                             / 2.0);
            // detect case where we don't have time to reach max velocity
            if (peakVelocity < mVelocity[i]) {
                mVelocity[i] = peakVelocity;
            }

            // calculates default times and distances
            mAccelerationTime[i] = (mVelocity[i] - fabs(mOvershotInitialVelocity[i])) / mAcceleration[i];
            mDecelerationTime[i] = mVelocity[i] / mAcceleration[i];
            mAccelerationDistance[i] =
                    fabs(mOvershotInitialVelocity[i]) * mAccelerationTime[i]
                    + 0.5 * mAcceleration[i] * mAccelerationTime[i] * mAccelerationTime[i];
            mDecelerationDistance[i] =
                    0.5 * mAcceleration[i] * mDecelerationTime[i] * mDecelerationTime[i];
            // calculates the finish time for a movement that reaches a constant velocity
            mFinishTime[i] =
                    mOvershotTime[i]
                    + (fabs(displacement) - fabs(mAccelerationDistance[i]) - mDecelerationDistance[i]) / mVelocity[i] // Constant velocity phase
                    + mAccelerationTime[i]
                    + mDecelerationTime[i];
        }
    }

    // compute max time
    mDuration = mFinishTime.MaxElement();
    // scale time to all arrive at same time
    if (mCoordination == LSPB_DURATION) {
        mTimeScale.SetSize(mDimension);
        if (mDuration > 0) {
            mTimeScale.RatioOf(mFinishTime, mDuration);
        } else {
            mTimeScale.SetAll(1.0);
        }
    }
    mIsSet = true;
}

void robLSPB::Evaluate(const double absoluteTime,
                       vctDoubleVec & position,
                       vctDoubleVec & velocity,
                       vctDoubleVec & acceleration) const
{

    // sanity checks
    if (!mIsSet) {
        cmnThrow("robLSPB::Evaluate trajectory parameters are not set yet");
    }
    if (position.size() != mDimension) {
        cmnThrow("robLSPB::Evaluate: position doesn't match dimension");
    }
    if (velocity.size() != mDimension) {
        cmnThrow("robLSPB::Evaluate: velocity doesn't match dimension");
    }
    if (acceleration.size() != mDimension) {
        cmnThrow("robLSPB::Evaluate: acceleration doesn't match dimension");
    }
    for (size_t i = 0;
         i < mDimension;
         ++i) {
        const double time = absoluteTime - mStartTime;
        double dimTime;
        if (mCoordination == LSPB_DURATION) {
            dimTime = time * mTimeScale[i];
        } else {
            dimTime = time;
        }
        const double time2 = dimTime * dimTime;

        // before trajectory
        if (time <= 0) {
            velocity[i] = mInitialVelocity[i];
            position[i] = mStart[i] + mInitialVelocity[i] * dimTime;
            acceleration[i] = 0;
        }

        // after trajectory
        else if (dimTime >= mFinishTime[i]) {
            position[i] = mFinish[i];
            velocity[i] = 0.0;
            acceleration[i] = 0.0;
        }

        // immediate deceleration phase to overshoot the desired position
        else if (dimTime <= mOvershotTime[i]) {
            const double overshotTime = time * mTimeScale[i];
            position[i] =
                    mStart[i]
                    + mInitialVelocity[i] * overshotTime
                    + mInitialDirection[i]
                    * (-0.5 * mAcceleration[i] * overshotTime * overshotTime);
            velocity[i] =
                    mInitialVelocity[i]
                    - mInitialDirection[i]* mAcceleration[i] * overshotTime;
            acceleration[i] = - mInitialDirection[i] * mAcceleration[i];
        }

        // acceleration phase
        else if (dimTime <= mAccelerationTime[i] + mOvershotTime[i]) {
            const double accelerationTime = (time - mOvershotTime[i]) * mTimeScale[i];
            position[i] =
                    mOvershotStart[i]
                    + mOvershotInitialVelocity[i] * accelerationTime
                    + mOvershotDirection[i]
                    * (0.5 * mAcceleration[i] * accelerationTime * accelerationTime);
            velocity[i] =
                    mOvershotInitialVelocity[i]
                    + mOvershotDirection[i] * (mAcceleration[i] * accelerationTime);
            acceleration[i] = mOvershotDirection[i] * mAcceleration[i];
        }

        // deceleration phase
        else if (dimTime >= (mFinishTime[i] - mDecelerationTime[i])) {
            position[i] =
                    mFinish[i]
                    + mOvershotDirection[i]
                    * (- 0.5 * mAcceleration[i] * mFinishTime[i] * mFinishTime[i]
                       + mAcceleration[i] * mFinishTime[i] * dimTime
                       - 0.5 * mAcceleration[i] * time2);
            velocity[i] =
                    mOvershotDirection[i]
                    * (mAcceleration[i] * mFinishTime[i]
                       - mAcceleration[i] * dimTime);
            acceleration[i] = - mOvershotDirection[i] * mAcceleration[i];
        }

        // constant velocity
        else {
            const double constantTime = (time - (mAccelerationTime[i] + mOvershotTime[i])) * mTimeScale[i];
            position[i] =
                    mOvershotStart[i]
                    + mOvershotDirection[i]
                    * (mAccelerationDistance[i]
                       + mVelocity[i] * constantTime);
            velocity[i] = mOvershotDirection[i] * mVelocity[i];
            acceleration[i] = 0.0;
        }

        if (mCoordination == LSPB_DURATION) {
            velocity[i] *= mTimeScale[i];
            acceleration[i] *= (mTimeScale[i] * mTimeScale[i]);
        }
    }
}

void robLSPB::Evaluate(const double absoluteTime,
                       vctDoubleVec & position) const
{
    vctDoubleVec temp(mDimension);
    Evaluate(absoluteTime, position, temp, temp);
}

const double & robLSPB::StartTime(void) const {
    return mStartTime;
}

const double & robLSPB::Duration(void) const {
    return mDuration;
}
