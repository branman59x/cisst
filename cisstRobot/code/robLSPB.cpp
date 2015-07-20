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
            mInitialVelocity[i] *= -1.0;
            displacement *= -1.0;
            mStart[i] *= -1.0;
            mFinish[i] *= -1.0;
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
                mOvershotTime[i] = -(mInitialVelocity[i] / mAcceleration[i]);
                mOvershotDistance[i] =
                        - 0.5 * mAcceleration[i] * mOvershotTime[i] * mOvershotTime[i];
                mOvershotInitialVelocity[i] = 0.0;
                mOvershotDirection[i] = 1.0;
            } else {
                mOvershotTime[i] = mInitialVelocity[i] / mAcceleration[i];
                mOvershotDistance[i] =
                        0.5 * mAcceleration[i] * mOvershotTime[i] * mOvershotTime[i];
                // test if this is really an overshot
                if (mOvershotDistance[i] <= mFinish[i]) {
                    // test if we start with an initial velocity too high and we need to deccelerate
                    if (mInitialVelocity[i] > mVelocity[i]) {
                        mOvershotTime[i] = (mInitialVelocity[i] - mVelocity[i]) / mAcceleration[i];
                        mOvershotDistance[i] = mInitialVelocity[i] * mOvershotTime[i]
                                 - 0.5 * mAcceleration[i] * mOvershotTime[i] * mOvershotTime[i];
                        std::cout<<mOvershotDistance[i]<<" OvershotDist "<<mOvershotTime[i]<<" OvershotTime\n";
                        mOvershotInitialVelocity[i] = mVelocity[i];
                        mOvershotDirection[i] = 1.0;
                    } else {
                        // really not an overshot
                        mOvershotTime[i] = 0;
                        mOvershotDistance[i] = 0;
                        mOvershotInitialVelocity[i] = mInitialVelocity[i];
                        mOvershotDirection[i] = 1.0;
                    }
                } else {
                    // we had an overshot, we stopped
                    mOvershotInitialVelocity[i] = 0.0;
                    // and then go backward
                    mOvershotDirection[i] = -1.0;
                }
            }
            std::cout<<"OD "<<mOvershotDistance[i]<<"\n";
            mOvershotStart[i] = mStart[i] + mOvershotDistance[i];
            displacement = mFinish[i] - (mStart[i] + mOvershotDistance[i]);
            // plan trajectory from the overshot point
            // calculates the max velocity that the velocity will hit (may be lower than predetermined)
            const double peakVelocity = sqrt((2.0 * mAcceleration[i] * displacement * mOvershotDirection[i]
                                              + mOvershotInitialVelocity[i] * mOvershotInitialVelocity[i])
                                             / 2.0);
            // detect case where we don't have time to reach max velocity
            if (peakVelocity < mVelocity[i]) {
                mVelocity[i] = peakVelocity;
            }

            // calculates default times and distances if we reached constant velocity
            mAccelerationTime[i] = (mVelocity[i] - mOvershotInitialVelocity[i]) / mAcceleration[i];
            mDecelerationTime[i] = mVelocity[i] / mAcceleration[i];
            mAccelerationDistance[i] =
                    mOvershotInitialVelocity[i] * mAccelerationTime[i]
                    + 0.5 * mAcceleration[i] * mAccelerationTime[i] * mAccelerationTime[i];
            mDecelerationDistance[i] =
                    0.5 * mAcceleration[i] * mDecelerationTime[i] * mDecelerationTime[i];
            std::cout<<"AT,DT,AD,DD "<<mAccelerationTime[i]<<" "<<mDecelerationTime[i]<<" "<<mAccelerationDistance[i]<<" "<<mDecelerationDistance[i]<<"\n";
            // calculates the finish time for a movement that reaches a constant velocity
            std::cout<<"Displacement "<<displacement<<"\n";
            std::cout<<"Const Phase "<<(displacement - mAccelerationDistance[i] - mDecelerationDistance[i] - mOvershotDistance[i]) / mVelocity[i]<<"\n";
            mFinishTime[i] =
                    mOvershotTime[i]
                    + (displacement - mAccelerationDistance[i] - mDecelerationDistance[i]) / mVelocity[i] // Constant velocity phase
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
        if (time <= 0) {
            velocity[i] = mInitialVelocity[i];
            position[i] = mStart[i] + mInitialVelocity[i] * time;
            acceleration[i] = 0;
            return;
        }
        if (dimTime >= mFinishTime[i]) {
            position[i] = mFinish[i];
            velocity[i] = 0.0;
            acceleration[i] = 0.0;}
        else if(dimTime <= mOvershotTime[i]){
            // immediate deceleration phase to overshoot the desired position
            std::cout<<"dimTime "<<dimTime<<" Pos:"<<position[i]<<" mOvershotTime "<<mOvershotTime[i]<<" Overshooting\n";
            position[i] =
                    mStart[i] + mInitialVelocity[i]*dimTime - 0.5*mAcceleration[i]*time2;
            velocity[i] =
                    mInitialVelocity[i] - mAcceleration[i] * dimTime;
            acceleration[i] = -fabs(mAcceleration[i]);
        }
        else if (dimTime <= mAccelerationTime[i]){
            // acceleration phase
            std::cout<<"ACCELERATING dimTime "<<dimTime<<" Pos:"<<position[i]<<"\n";
            position[i] = mOvershotStart[i] + mOvershotInitialVelocity[i]*dimTime + 0.5*mAcceleration[i]*time2;
            velocity[i] = mAcceleration[i] * dimTime + mOvershotInitialVelocity[i];
            acceleration[i] = fabs(mAcceleration[i]);
        } else if (dimTime >= (mFinishTime[i] - mDecelerationTime[i])) {
            // deceleration phase
            std::cout<<"DECELERATING dimTime "<<dimTime<<" Pos:"<<position[i]<<"\n";
            position[i] =
                    mFinish[i]
                    - 0.5 * mAcceleration[i] * mFinishTime[i] * mFinishTime[i]
                    + mAcceleration[i] * mFinishTime[i] * dimTime
                    - 0.5 * mAcceleration[i] * time2;
            velocity[i] =
                    mAcceleration[i] * mFinishTime[i]
                    - mAcceleration[i] * dimTime;
            acceleration[i] = -fabs(mAcceleration[i]);
        } else {
            std::cout<<"CONSTANT| dimTime "<<dimTime<<" Pos:"<<position[i]<<"\n";
            // constant velocity phase
            std::cout<<"in here\n";
            position[i] = mAccelerationDistance[i] + mOvershotStart[i] + mVelocity[i] * (dimTime - mAccelerationTime[i] - mOvershotTime[i]);
            velocity[i] = mVelocity[i];
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
