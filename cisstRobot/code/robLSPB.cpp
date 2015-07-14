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
    mSecondAccelTime.SetSize(mDimension);//Time of Acceleration of a trajectory after an overshoot
    mSecondAccelDistance.SetSize(mDimension);//Distance of Acceleration of a trajectory after an overshoot
    mSecondDistance.SetSize(mDimension);//Total distance of a trajectory after an overshoot
    mOvershotTime.SetSize(mDimension);//Time of Overshoot
    mOvershotTime.Zeros();
    mAccelerationDistance.SetSize(mDimension);
    mDecelerationDistance.SetSize(mDimension);
    mAccelerationDistance.Zeros();
    mDecelerationDistance.Zeros();
    //mOvershoot = false;//Overshot Test
    mTotalTime.SetSize(mDimension);//Time that includes overshot for time scale
    mOvershoot.SetSize(mDimension);
    mOvershoot.Zeros();
    //bool mOvershoot [mDimension];
    // compute trajectory parameters
    for (size_t i = 0;
         i < mDimension;
         ++i) {
        // compute direction
        const double distance = finish[i] - start[i];
        if (distance < 0.0) {
            mVelocity[i] *= -1.0;;
            mAcceleration[i] *= -1.0;
        } else if (distance == 0 && mInitialVelocity[i] == 0) {
            mVelocity[i] = 0.0;
            mAcceleration[i] = 0.0;
        }
        // compute time if distance != 0
        if (distance != 0 || mInitialVelocity[i] != 0) {
            if (mVelocity[i] == 0.0) {
                cmnThrow("robLSPB::Set:  be greater than zero");
            }
            if (mAcceleration[i] == 0.0) {
                cmnThrow("robLSPB::Set: acceleration must be greater than zero");
            }
            mAccelerationTime[i] = (mVelocity[i]-mInitialVelocity[i])/mAcceleration[i];
            mDecelerationTime[i] = mVelocity[i] / mAcceleration[i];
            mAccelerationDistance[i] = mInitialVelocity[i] * mAccelerationTime[i] + 0.5 * mAcceleration[i]
                    * mAccelerationTime[i] * mAccelerationTime[i];
            mDecelerationDistance[i] = 0.5 * mAcceleration[i]
                    * mDecelerationTime[i] * mDecelerationTime[i];
            // check distance over max accel and max decel at end of
            // acceleration phase to see if we're past mid-point.
            if (fabs(mAccelerationDistance[i] + mDecelerationDistance[i]) >= fabs(distance)) {
                //calculates the times with an initial velocity of 0
                if(mInitialVelocity[i]==0){
                    mFinishTime[i] = 2 * sqrt(distance / mAcceleration[i]);
                    mAccelerationTime[i] = mFinishTime[i] * 0.5;
                    mDecelerationTime[i] = mFinishTime[i] * 0.5;
                }
                else{
                    //calculates the times with an initial velocity that is not 0
                    double tempMaxVelocity = mVelocity[i];
                    std::cout<<tempMaxVelocity<<"|Temp\n";
                    //calculates the max velocity that the velocity will hit (may be lower than predetermined)
                    mVelocity[i] = sqrt((2*mAcceleration[i]*distance + mInitialVelocity[i]*mInitialVelocity[i])/2);
                    std::cout<<mVelocity[i]<<"|mV\n";
                    if(fabs(mVelocity[i]) > fabs(tempMaxVelocity))
                    {
                        mVelocity[i] = tempMaxVelocity;
                    }
                    mAccelerationTime[i] = (mVelocity[i]-fabs(mInitialVelocity[i]))/fabs(mAcceleration[i]);
                    std::cout<<mVelocity[i]<<" : mVelocity "<<mAccelerationTime[i]<<" AccelTime\n";
                    mDecelerationTime[i] = mVelocity[i]/fabs(mAcceleration[i]);
                    if(mAccelerationTime[i] <= 0)
                        mDecelerationTime[i] += fabs(mAccelerationTime[i]);
                    mFinishTime[i] = mAccelerationTime[i] + mDecelerationTime[i];
                }
                //calculates the distances for a movement that does not reach a constant velocity
                //and has an initial velocity that does not equal 0
                mAccelerationDistance[i] = mInitialVelocity[i] * mAccelerationTime[i] + 0.5 * mAcceleration[i]
                        * mAccelerationTime[i] * mAccelerationTime[i];
                mDecelerationDistance[i] = 0.5 * mAcceleration[i]
                        * mDecelerationTime[i] * mDecelerationTime[i];
            } else {
                //calculates the finish time for a movement that reaches a constant velocity
                std::cout<<"Constant V\n";
                //if(mAccelerationTime[i] <= 0)
                //    mDecelerationTime[i] += fabs(mAccelerationTime[i]);
                std::cout<<"AccelDistance "<<mAccelerationDistance[i]<<"\n";
                std::cout<<"Constant Distance "<<distance - fabs(mAccelerationDistance[i]) - mDecelerationDistance[i]<<"\n";
                mFinishTime[i] = (distance-fabs(mAccelerationDistance[i])-mDecelerationDistance[i]) / mVelocity[i] + fabs(mAccelerationTime[i]) + mDecelerationTime[i];
                std::cout<<mFinishTime[i]<< " FinishTime\n";
            }
            mTotalTime[i] = mFinishTime[i];
            //checks for overshoot
            std::cout<<mDecelerationDistance[i]<<":decelDistance "<< distance << " distance "<<mDecelerationTime[i]<<":decelTIME "<< mAccelerationTime[i]<<":acelTime\n";
            if(fabs(mDecelerationDistance[i]) > fabs(distance))
            {
                std::cout<<"joint "<<i<<"\n";
                mOvershoot[i] = true;
                //prepares variables for an immediate deceleration
                mAcceleration[i] = fabs(mAcceleration[i]);
                mDecelerationTime[i] = mInitialVelocity[i]/mAcceleration[i];
                mDecelerationDistance[i] = 0.5 * mAcceleration[i]
                        * mDecelerationTime[i] * mDecelerationTime[i];
                //calculates variables for a new trajectory after deceleration
                mSecondDistance[i] = fabs(mStart[i] + mDecelerationDistance[i]) - mFinish[i];
                double tempMaxVelocity = mVelocity[i];
                mVelocity[i] = sqrt(fabs(mAcceleration[i]*mSecondDistance[i]));
                if(fabs(mVelocity[i]) > fabs(tempMaxVelocity))
                {
                    mVelocity[i] = fabs(tempMaxVelocity);
                }
                mSecondAccelTime[i] = fabs(mVelocity[i]/mAcceleration[i]);
                mSecondAccelDistance[i] = 0.5*mAcceleration[i]*mSecondAccelTime[i]*mSecondAccelTime[i];
                mFinishTime[i] = mSecondAccelTime[i]*2 +
                        (fabs(mSecondDistance[i])-mSecondAccelDistance[i]*2)/mVelocity[i];
                mTotalTime[i] = mFinishTime[i] + mDecelerationTime[i];
            }
            else
                mOvershoot[i] = false;
        } else{
            mAccelerationTime[i] = 0.0;
            mDecelerationTime[i] = 0.0;
            mFinishTime[i] = 0.0;
        }
    }
    std::cout<<mOvershoot[0]<<" "<<mOvershoot[1]<<" "<<mOvershoot[2]<<"\n";
    // compute max time
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
        if(mOvershoot[i])
        {
            std::cout<<mOvershoot[0]<<" "<<mOvershoot[1]<<" "<<mOvershoot[2]<<"\n";
            const double time = absoluteTime - mStartTime;
            double dimTime;
            if (mCoordination == LSPB_DURATION) {
                dimTime = time * mTimeScale[i];
            } else {
                dimTime = time;
            }
            std::cout<<"OVERSHOOT decelTime: "<<i<<","<<mDecelerationTime[i]<<"\n";
            std::cout<<"dimTime: "<<dimTime<<"\n";
            const double time2 = dimTime * dimTime;
            if (time <= 0) {
                position.ForceAssign(mStart);
                velocity.Zeros();
                velocity.ForceAssign(mInitialVelocity);
                acceleration.Zeros();
                return;
            }
            // immediate deceleration phase to overshoot the desired position
            std::cout<<"Decelerating in the Overshoot|Velocity "<<velocity[i]<<" |Position "<<position[i]<<" |Acceleration "<<mAcceleration[i]<<"\n";
            position[i] =
                    mStart[i] + mInitialVelocity[i]*dimTime + -0.5*fabs(mAcceleration[i])*time2;
            velocity[i] =
                    mInitialVelocity[i] - mAcceleration[i] * dimTime;
            acceleration[i] = -mAcceleration[i];
            //checks if deceleration has finished
            if(dimTime >= mDecelerationTime[i])
            {
                std::cout<<"dimTime > decelTime: " << mDecelerationTime[i] << "\n";
                mOvershoot[i] = false;
                //sets variables for the trajectory to get to the finish point
                mStart[i] += mDecelerationDistance[i];
                mDecelerationDistance[i] = mSecondAccelDistance[i];
                mAccelerationDistance[i] = mSecondAccelDistance[i];
                mAccelerationTime[i] = mSecondAccelTime[i];
                mOvershotTime[i] = mDecelerationTime[i];
                mDecelerationTime[i] = mSecondAccelTime[i];
                mInitialVelocity[i] = velocity[i];
                mVelocity[i] *= -1;
                mAcceleration[i] *= -1;
                std::cout<<"AT|"<<mAccelerationTime[i]<<" |DT| "<<mDecelerationTime[i]<<"\n";
            }

        }
        if(!mOvershoot[i])
        {
            const double time = absoluteTime - mStartTime;
            if (time <= 0) {
                position.ForceAssign(mStart);
                velocity.Zeros();
                velocity.ForceAssign(mInitialVelocity);
                acceleration.Zeros();
                return;
            }

            //        for (size_t i = 0;
            //             i < mDimension;
            //             ++i) {
            double conditionTime;
            if (mCoordination == LSPB_DURATION) {
                conditionTime = time*mTimeScale[i] - mOvershotTime[i];
            } else {
                conditionTime = time - mOvershotTime[i];
            }
            if (conditionTime >= mFinishTime[i]) {
                position[i] = mFinish[i];
                velocity[i] = 0.0;
                acceleration[i] = 0.0;
            } else {
                // acceleration phase
                if (conditionTime <= mAccelerationTime[i]) {
                    std::cout<<"Accelerating|ConditionTime "<<conditionTime<<" |Velocity "<<velocity[i]<<" |Position "<<position[i]<<" |Acceleration "<<mAcceleration[i]<<"\n";
                    position[i] =
                            mStart[i]
                            + 0.5*mAcceleration[i] * conditionTime*conditionTime + mInitialVelocity[i]*conditionTime;
                    velocity[i] = mAcceleration[i] * conditionTime + mInitialVelocity[i];
                    acceleration[i] = mAcceleration[i];
                } else if (conditionTime >= (mFinishTime[i] - mDecelerationTime[i])|| -conditionTime > mAccelerationTime[i]) {
                    // deceleration phase
                    std::cout<<"Decelerating|ConditionTime "<<conditionTime<<" |Velocity "<<velocity[i]<<" |Position "<<position[i]<<" |Acceleration "<<mAcceleration[i]<<"\n";
                    // deceleration when the max velocity is lower than the initial velocity
                    if(-conditionTime > mAccelerationTime[i])
                    {
                        std::cout<<"decel\n";
                        position[i] =
                                mStart[i] + mInitialVelocity[i]*conditionTime + -0.5*fabs(mAcceleration[i])*conditionTime*conditionTime;
                        velocity[i] =
                                mInitialVelocity[i] + acceleration[i] * conditionTime;
                    }
                    // regular deceleration phase
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
                    std::cout<<"Constant|ConditionTime "<<conditionTime<<" |Velocity "<<velocity[i]<<" |Position "<<position[i]<<" |Acceleration "<<mAcceleration[i]<<"\n";
                    // constant velocity phase
                    position[i] = 0.5 * mAcceleration[i] * mAccelerationTime[i] * mAccelerationTime[i] + mStart[i]
                            + mInitialVelocity[i]*fabs(mAccelerationTime[i]) + mVelocity[i]*(conditionTime - fabs(mAccelerationTime[i]));
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
