/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-    */
/* ex: set filetype=cpp softtabstop=4 shiftwidth=4 tabstop=4 cindent expandtab: */

/*
  Author(s):  Anton Deguet
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
    mSecondAccelTime.SetSize(mDimension);
    mSecondAccelDistance.SetSize(mDimension);
    mSecondDistance.SetSize(mDimension);
    mOvershotTime.SetSize(mDimension);
    mOvershotTime.Zeros();
    mOvershoot = false;
    mTotalTime.SetSize(mDimension);
    //mDecelToMaxTime.SetSize(mDimension);
    //mDecelToMax
    // compute trajectory parameters
    for (size_t i = 0;
         i < mDimension;
         ++i) {
        // compute direction
        const double distance = finish[i] - start[i];
        if (distance < 0.0) {
            mVelocity[i] *= -1.0;;
            mAcceleration[i] *= -1.0;
        } else if (distance == 0) {
            mVelocity[i] = 0.0;
            mAcceleration[i] = 0.0;
        }

        // compute time if distance != 0
        if (distance != 0) {
            if (mVelocity[i] == 0.0) {
                cmnThrow("robLSPB::Set: velocity must be greater than zero");
            }
            if (mAcceleration[i] == 0.0) {
                cmnThrow("robLSPB::Set: acceleration must be greater than zero");
            }
            mAccelerationTime[i] = (mVelocity[i]-mInitialVelocity[i])/mAcceleration[i];
            mDecelerationTime[i] = mVelocity[i] / mAcceleration[i];
            // check distance over max accel and max decel at end of
            // acceleration phase to see if we're past mid-point.
            mAccelerationDistance = mInitialVelocity[i] * mAccelerationTime[i] + 0.5 * mAcceleration[i]
                    * mAccelerationTime[i] * mAccelerationTime[i];
            mDecelerationDistance = 0.5 * mAcceleration[i]
                    * mDecelerationTime[i] * mDecelerationTime[i];
            if (fabs(mAccelerationDistance + mDecelerationDistance) >= fabs(distance)) {
                //calculates the times with an initial velocity of 0
                if(mInitialVelocity[i]==0){
                    mFinishTime[i] = 2 * sqrt(distance / mAcceleration[i]);
                    mAccelerationTime[i] = mFinishTime[i] * 0.5;
                    mDecelerationTime[i] = mFinishTime[i] * 0.5;
                }
                else{
                    //calculates the times with an initial velocity that is not 0
                    std::cout<<"Hello\n";
                    double tempMaxVelocity = mVelocity[i];
                    mVelocity[i] = sqrt((2*mAcceleration[i]*distance + mInitialVelocity[i]*mInitialVelocity[i])/2);
                    if(mVelocity[i] > tempMaxVelocity)
                    {
                        mVelocity[i] = tempMaxVelocity;
                    }
                    mAccelerationTime[i] = (mVelocity[i]-mInitialVelocity[i])/mAcceleration[i];
                    mDecelerationTime[i] = mVelocity[i]/mAcceleration[i];
                    if(mAccelerationTime[i] <= 0)
                        mDecelerationTime[i] += fabs(mAccelerationTime[i]);

                    std::cout<<"FIRST DECELERATION TIME "<<(mDecelerationTime[i])<<"\n";
                    mFinishTime[i] = mAccelerationTime[i] + mDecelerationTime[i];
                }
                std::cout<<"Setting Distances Triangle\n";
                mAccelerationDistance = mInitialVelocity[i] * mAccelerationTime[i] + 0.5 * mAcceleration[i]
                        * mAccelerationTime[i] * mAccelerationTime[i];
                mDecelerationDistance = 0.5 * mAcceleration[i]
                        * mDecelerationTime[i] * mDecelerationTime[i];
            } else {
                //calculates the finish time for a movement that reaches a constant velocity
                std::cout<<"Constant V\n";
                mFinishTime[i] = (distance-mAccelerationDistance-mDecelerationDistance) / mVelocity[i] + mAccelerationTime[i] + mDecelerationTime[i];
            }
            mTotalTime[i] = mFinishTime[i];
//            if(mVelocity[i] < mInitialVelocity[i])
//            {
//                double dif = mInitialVelocity[i] - mVelocity[i];
//                double decelToMaxTime = fabs(dif/mAcceleration[i]);
//                std::cout<<"Acceleration Time "<<mAccelerationTime[i]<<"\n";
//                mTotalTime[i] = mTotalTime[i] + decelToMaxTime - mAccelerationTime[i];
//                mFinishTime[i] = mFinishTime[i] + decelToMaxTime - mAccelerationTime[i];
//                mAccelerationTime[i] = 0;
//                std::cout<<decelToMaxTime<<" DTMT "<<mTotalTime<<" Total Time " <<mFinishTime[i]<< " FinishTime " <<"\n";
//            }
            if(mDecelerationDistance > distance)
            {
                std::cout<<"Overshot\n";
                mOvershoot = true;
                mSecondDistance[i] = abs(mFinish[i] - mDecelerationDistance) + mStart[i];
                std::cout<<mDecelerationDistance<<": Deceleration Distance\n";
                double tempMaxVelocity = mVelocity[i];
                mVelocity[i] = sqrt(fabs(mAcceleration[i]*mSecondDistance[i]));
                if(mVelocity[i] > tempMaxVelocity)
                {
                    mVelocity[i] = tempMaxVelocity;
                }
                mSecondAccelTime[i] = mVelocity[i]/mAcceleration[i];
                std::cout<<mSecondAccelTime[i]<< " Second Accel Time\n";
                std::cout<<"mVelocity "<<mVelocity[i]<<"\n";
                mSecondAccelDistance[i] = 0.5*mAcceleration[i]*mSecondAccelTime[i]*mSecondAccelTime[i];
                mFinishTime[i] = mSecondAccelTime[i]*2 +
                        (fabs(mSecondDistance[i])-mSecondAccelDistance[i]*2)/mVelocity[i];
                mTotalTime[i] = mFinishTime[i] + mDecelerationTime[i];
                std::cout<<"Second Distance: "<<mSecondDistance[i]<<" Second Accel Distance "<<mSecondAccelDistance[i]<<"\n";
            }
        } else {
            std::cout<<"Zeros\n";
            mAccelerationTime[i] = 0.0;
            mDecelerationTime[i] = 0.0;
            mFinishTime[i] = 0.0;
        }
    }
    // compute max time
    std::cout<<mFinishTime<< " : Finish Time\n";
    mDuration = mTotalTime.MaxElement();
    // scale time to all arrive at same time
    if (mCoordination == LSPB_DURATION) {
        mTimeScale.SetSize(mDimension);
        if (mDuration > 0) {
            mTimeScale.RatioOf(mTotalTime, mDuration);
        } else {
            mTimeScale.SetAll(1.0);
        }
    }
    mIsSet = true;
}

void robLSPB::Evaluate(const double absoluteTime,
                       vctDoubleVec & position,
                       vctDoubleVec & velocity,
                       vctDoubleVec & acceleration)
{
  // std::cout<<"AccelTime "<<mAccelerationTime<<" DecelTime"<<mDecelerationTime<<" FinishTime"<<mFinishTime<<" Initial V"<<mInitialVelocity<<"\n";
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

    if(mOvershoot)
    {
        std::cout<<"Im Overshooting\n";
        for (size_t i = 0;
             i < mDimension;
             ++i) {
            const double time = absoluteTime - mStartTime;
            double dimTime;
            if (mCoordination == LSPB_DURATION) {
                std::cout<<"st: "<<mStartTime<<" , AT: "<<absoluteTime<<"\n";
                dimTime = time * mTimeScale[i];
            } else {
                std::cout<<"thing2\n";
                dimTime = time;
            }
            const double time2 = dimTime * dimTime;
            if (time <= 0) {
                position.ForceAssign(mStart);
                velocity.Zeros();
                velocity = mInitialVelocity;
                acceleration.Zeros();
                return;
            }
        // deceleration phase
            std::cout<<mInitialVelocity[i]<<":IV\n";
        position[i] =
                mStart[i] + mInitialVelocity[i]*dimTime + -0.5*fabs(mAcceleration[i])*time2;
        velocity[i] =
                mInitialVelocity[i] + acceleration[i] * dimTime;
        acceleration[i] = -mAcceleration[i];
        std::cout<<velocity[i]<<":Velocity||Position: "<<position[i]<<"\n";
        std::cout<<dimTime<<":dimTime\n";
        if(dimTime >= mDecelerationTime[i])
        {
            mStart[i] += mDecelerationDistance;
            mOvershoot = false;
            mDecelerationDistance = mSecondAccelDistance[i];
            mAccelerationDistance = mSecondAccelDistance[i];
            std::cout<<mDecelerationDistance<<"DecelDist\n";
            mAccelerationTime[i] = mSecondAccelTime[i];
            mOvershotTime[i] = mDecelerationTime[i];
            std::cout<<mOvershotTime[i]<<" Overshot Time\n";
            mDecelerationTime[i] = mSecondAccelTime[i];
            mInitialVelocity[i] = velocity[i];
            mVelocity[i] *= -1;
            mAcceleration[i] *= -1;
            std::cout<<mFinishTime<<mDecelerationTime<<mAccelerationTime<<" FT,DT,AT\n";
        }
      }
    }
    else
    {
        const double time = absoluteTime - mStartTime;
        if (time <= 0) {
            position.ForceAssign(mStart);
            velocity.Zeros();
            velocity = mInitialVelocity;
            acceleration.Zeros();
            return;
        }

        for (size_t i = 0;
             i < mDimension;
             ++i) {
            double dimTime;
            double conditionTime;
            if (mCoordination == LSPB_DURATION) {
                dimTime = time * mTimeScale[i];
                conditionTime = time*mTimeScale[i] - mOvershotTime[i];
            } else {
                dimTime = time;
                conditionTime = time - mOvershotTime[i];
            }
            //const double time2 = dimTime * dimTime;
            if (conditionTime >= mFinishTime[i]) {
                position[i] = mFinish[i];
                velocity[i] = 0.0;
                acceleration[i] = 0.0;
            } else {
                // acceleration phase
                if (conditionTime <= mAccelerationTime[i]) {
                    std::cout<<"InitialVelocity: "<<mInitialVelocity<<"\n";
                    std::cout<<"Accelerating: dimTime:" <<dimTime<<" conditionTime "<<conditionTime<<" Velocity: "<<velocity[i]<<" Position: "<<position[i]<< "Acceleration: "<<acceleration[i]<<"\n";
                    position[i] =
                        mStart[i]
                            + 0.5*mAcceleration[i] * conditionTime*conditionTime + mInitialVelocity[i]*conditionTime;
                    velocity[i] = mAcceleration[i] * conditionTime + mInitialVelocity[i];
                    acceleration[i] = mAcceleration[i];
                } else if (conditionTime >= (mFinishTime[i] - mDecelerationTime[i])|| -conditionTime > mAccelerationTime[i]) {
                    // deceleration phase
                    std::cout<<"Decelerating: dimTime:" <<dimTime<<" conditionTime "<<conditionTime<<" Velocity: "<<velocity[i]<<" Position: "<<position[i]<< "Acceleration: "<<acceleration[i]<<" DecelTime "<<mDecelerationTime[i]<<"\n";
                    if(-conditionTime > mAccelerationTime[i])
                    {
                        position[i] =
                                mInitialVelocity[i]*conditionTime + -0.5*fabs(mAcceleration[i])*conditionTime*conditionTime;
                        velocity[i] =
                                mInitialVelocity[i] + acceleration[i] * conditionTime;
                    }
                    else
                    {
                        position[i] =
                            mFinish[i]
                            - 0.5 * mAcceleration[i] * mFinishTime[i] * mFinishTime[i]
                            + mAcceleration[i] * mFinishTime[i] * conditionTime
                            - 0.5 * mAcceleration[i] * conditionTime*conditionTime;
                        velocity[i] =
                            mAcceleration[i] * mFinishTime[i]
                            - mAcceleration[i] * conditionTime;
                    }
                    acceleration[i] = -mAcceleration[i];
                } else {
                    std::cout<<"Constant\n";
                    // constant velocity
                    position[i] = 0.5 * mAcceleration[i] * mAccelerationTime[i] * mAccelerationTime[i] + mStart[i]
                            + mInitialVelocity[i]*mAccelerationTime[i] + mVelocity[i]*(conditionTime - mAccelerationTime[i]);
                    velocity[i] = mVelocity[i];
                    acceleration[i] = 0.0;
                }
            }
        }
    }
}

void robLSPB::Evaluate(const double absoluteTime,
                       vctDoubleVec & position)
{
    mTemp.SetSize(mDimension);
    Evaluate(absoluteTime, position, mTemp, mTemp);
}

double & robLSPB::StartTime(void) {
    return mStartTime;
}

double robLSPB::Duration(void) const {
    return mDuration;
}
