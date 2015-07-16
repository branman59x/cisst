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
    mTotalTime.SetSize(mDimension);//Time that includes overshot for time scale
    mOvershoot.SetSize(mDimension);
    mOvershoot.Zeros();
    mFlipper.SetSize(mDimension);// flips the graph if there is a negative dispalcement
    mFlipper.Zeros();

    // compute trajectory parameters
    for (size_t i = 0;
         i < mDimension;
         ++i) {
        // make sure that the max velocity and max acceleration are positive
        mVelocity[i] = fabs(mVelocity[i]);
        mAcceleration[i] = fabs(mAcceleration[i]);
        // compute direction
        double distance = finish[i] - start[i];
        if (distance < 0.0) {
            //sets variables so that the graph will flip as if the displacement was positive
            mInitialVelocity[i] *= -1;
            distance *= -1;
            double temp = mStart[i];
            mStart[i] = mFinish[i];
            mFinish[i] = temp;
            mFlipper[i] = true;
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
            //calculates default times and distances
            mAccelerationTime[i] = (mVelocity[i]-mInitialVelocity[i])/mAcceleration[i];
            mDecelerationTime[i] = mVelocity[i] / mAcceleration[i];
            mAccelerationDistance[i] = mInitialVelocity[i] * mAccelerationTime[i] + 0.5 * mAcceleration[i]
                    * mAccelerationTime[i] * mAccelerationTime[i];
            mDecelerationDistance[i] = 0.5 * mAcceleration[i]
                    * mDecelerationTime[i] * mDecelerationTime[i];
            // check distance over max accel and max decel at end of
            // acceleration phase to see if we're past mid-point.
            if ((mAccelerationDistance[i] + mDecelerationDistance[i]) >= distance) {
                //calculates the times with an initial velocity of 0
                if(mInitialVelocity[i]==0){
                    mFinishTime[i] = 2 * sqrt(distance / mAcceleration[i]);
                    mAccelerationTime[i] = mFinishTime[i] * 0.5;
                    mDecelerationTime[i] = mFinishTime[i] * 0.5;
                }
                else{
                    //calculates the max velocity that the velocity will hit (may be lower than predetermined)
                    double tempMaxVelocity = mVelocity[i]; 
                    mVelocity[i] = sqrt((2*mAcceleration[i]*distance + mInitialVelocity[i]*mInitialVelocity[i])/2);
                    if(mVelocity[i] > tempMaxVelocity)
                    {
                        mVelocity[i] = tempMaxVelocity;
                    }
                    //calculates the times with an initial velocity that is not 0
                    mAccelerationTime[i] = (mVelocity[i]-mInitialVelocity[i])/(mAcceleration[i]);
                    mDecelerationTime[i] = mVelocity[i]/mAcceleration[i];
                    std::cout<<"AccelTime" <<mAccelerationTime[i] <<" DecelTime "<<mDecelerationTime[i]<<"\n";
                    if(mAccelerationTime[i] <= 0)
                        mDecelerationTime[i] += mAccelerationTime[i];
                    mFinishTime[i] = mAccelerationTime[i] + mDecelerationTime[i];
                    std::cout<<"AccelTime" <<mAccelerationTime[i] <<" DecelTime "<<mDecelerationTime[i]<<"\n";
                    std::cout<<mFinishTime[i]<<" FT\n";
                }
                //calculates the distances for a movement that does not reach a constant velocity
                //and has an initial velocity that does not equal 0
                mAccelerationDistance[i] = mInitialVelocity[i] * mAccelerationTime[i] + 0.5 * mAcceleration[i]
                        * mAccelerationTime[i] * mAccelerationTime[i];
                mDecelerationDistance[i] = 0.5 * mAcceleration[i]
                        * mDecelerationTime[i] * mDecelerationTime[i];
            } else {
                //calculates the finish time for a movement that reaches a constant velocity
                mFinishTime[i] = (distance-fabs(mAccelerationDistance[i])-mDecelerationDistance[i]) / mVelocity[i] + fabs(mAccelerationTime[i]) + mDecelerationTime[i];
            }
            mTotalTime[i] = mFinishTime[i];

            if(mAccelerationDistance[i] < 0)
                mDecelerationDistance[i] -= mAccelerationDistance[i];
            //checks for overshoot or if displacement is opposite initial velocity
            if(mDecelerationDistance[i] > distance || (mInitialVelocity[i] != 0 && mInitialVelocity[i]/fabs(mInitialVelocity[i]) != distance/fabs(distance)))
            {
                mOvershoot[i] = true;
                //prepares variables for an immediate deceleration
                if((mInitialVelocity[i] != 0 && mInitialVelocity[i]/fabs(mInitialVelocity[i]) != distance/fabs(distance)))
                    mAcceleration[i] = -mAcceleration[i];
                mDecelerationTime[i] = fabs(mInitialVelocity[i])/fabs(mAcceleration[i]);
                mDecelerationDistance[i] = 0.5 * mAcceleration[i]
                        * mDecelerationTime[i] * mDecelerationTime[i];

                //calculates variables for a new trajectory after deceleration
                mSecondDistance[i] = mFinish[i] - (mStart[i] + mDecelerationDistance[i]);
                double tempMaxVelocity = mVelocity[i];
                mVelocity[i] = sqrt(fabs(mAcceleration[i]*mSecondDistance[i]));
                if(mVelocity[i] > tempMaxVelocity)
                {
                    mVelocity[i] = (tempMaxVelocity);
                }
                mSecondAccelTime[i] = mVelocity[i]/fabs(mAcceleration[i]);
                mSecondAccelDistance[i] = 0.5*mAcceleration[i]*mSecondAccelTime[i]*mSecondAccelTime[i];
                mFinishTime[i] = mSecondAccelTime[i]*2 +
                        (fabs(mSecondDistance[i])-fabs(mSecondAccelDistance[i]*2))/mVelocity[i];
                mTotalTime[i] = mFinishTime[i] + mDecelerationTime[i];
            }
            else
                mOvershoot[i] = false;
        } else{
            mAccelerationTime[i] = 0.0;
            mDecelerationTime[i] = 0.0;
            mFinishTime[i] = 0.0;
        }
        // changes variables back to their initial states in order to flip the graph in the evaluation
        if(mFlipper[i]){
            double temp = mStart[i];
            mStart[i] = mFinish[i];
            mFinish[i] = temp;
            mAcceleration[i] *= -1;
            if(!mOvershoot[i])
                mVelocity[i] *= -1;
            mInitialVelocity[i] *= -1;
        }
    }
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
            const double time = absoluteTime - mStartTime;
            double dimTime;
            if (mCoordination == LSPB_DURATION) {
                dimTime = time * mTimeScale[i];
            } else {
                dimTime = time;
            }
            const double time2 = dimTime * dimTime;
            if (time <= 0) {
                position.ForceAssign(mStart);
                velocity.Zeros();
                velocity.ForceAssign(mInitialVelocity);
                acceleration.Zeros();
                return;
            }
            // immediate deceleration phase to overshoot the desired position
            position[i] =
                    mStart[i] + mInitialVelocity[i]*dimTime - 0.5*mAcceleration[i]*time2;
            velocity[i] =
                    mInitialVelocity[i] - mAcceleration[i] * dimTime;
            acceleration[i] = -fabs(mAcceleration[i]);
            //checks if deceleration has finished
            if(dimTime >= mDecelerationTime[i])
            {
                mOvershoot[i] = false;
                //sets variables for the trajectory to get to the finish point
                if(mFlipper[i])
                    mStart[i] -= mDecelerationDistance[i];
                else
                    mStart[i] += mDecelerationDistance[i];
                mDecelerationDistance[i] = mSecondAccelDistance[i];
                mAccelerationDistance[i] = mSecondAccelDistance[i];
                mAccelerationTime[i] = mSecondAccelTime[i];
                mOvershotTime[i] = mDecelerationTime[i];
                mDecelerationTime[i] = mSecondAccelTime[i];
                mAcceleration[i] *= -1;
                if(mInitialVelocity[i] >= 0)
                    mVelocity[i]*=-1;
                mInitialVelocity[i] = velocity[i];
            }

        }
        if (!mOvershoot[i])
        {
            const double time = absoluteTime - mStartTime;
            if (time <= 0) {
                position.ForceAssign(mStart);
                velocity.Zeros();
                velocity.ForceAssign(mInitialVelocity);
                acceleration.Zeros();
                return;
            }
            double dimTime;
            if (mCoordination == LSPB_DURATION) {
                dimTime = time*mTimeScale[i] - mOvershotTime[i];
            } else {
                dimTime = time - mOvershotTime[i];
            }
            double time2 = dimTime * dimTime;
            if (dimTime >= mFinishTime[i]) {
                position[i] = mFinish[i];
                velocity[i] = 0.0;
                acceleration[i] = 0.0;
            } else {
                // acceleration phase
                if (dimTime <= mAccelerationTime[i]) {
                    position[i] = mStart[i] + mInitialVelocity[i]*dimTime + 0.5*mAcceleration[i]*time2;
                    velocity[i] = mAcceleration[i] * dimTime + mInitialVelocity[i];
                    acceleration[i] = fabs(mAcceleration[i]);                    
                } else if (dimTime >= (mFinishTime[i] - mDecelerationTime[i])|| -dimTime > mAccelerationTime[i]) {
                    // deceleration phase

                    // deceleration when the max velocity is lower than the initial velocity
                    if(-dimTime > mAccelerationTime[i])
                    {
                        position[i] =
                                mStart[i] + mInitialVelocity[i]*dimTime - 0.5*mAcceleration[i]*time2;
                        velocity[i] =
                                mInitialVelocity[i] - mAcceleration[i] * dimTime;
                    }
                    // regular deceleration phase
                    else
                    {
                        position[i] =
                                mFinish[i]
                                - 0.5 * mAcceleration[i] * mFinishTime[i] * mFinishTime[i]
                                + mAcceleration[i] * mFinishTime[i] * dimTime
                                - 0.5 * mAcceleration[i] * time2;
                        velocity[i] =
                                mAcceleration[i] * mFinishTime[i]
                                - mAcceleration[i] * dimTime;
                    }
                    acceleration[i] = -fabs(mAcceleration[i]);
                } else {
                    // constant velocity phase
                    // checks if initial velocity is greater than the max velocity and has not overshot
                    if((mAcceleration[i] > 0 && mInitialVelocity[i] > 0 && mInitialVelocity[i] > mVelocity[i]) || (mAcceleration[i] < 0 && mInitialVelocity[i] < 0 && fabs(mInitialVelocity[i]) > fabs(mVelocity[i])))
                        position[i] = -0.5 * mAcceleration[i] * mAccelerationTime[i] * mAccelerationTime[i] + mStart[i]
                                + mInitialVelocity[i]*fabs(mAccelerationTime[i]) + mVelocity[i]*(dimTime - fabs(mAccelerationTime[i]));
                    else
                        position[i] = 0.5 * mAcceleration[i] * mAccelerationTime[i] * mAccelerationTime[i] + mStart[i]
                                + mInitialVelocity[i]*fabs(mAccelerationTime[i]) + mVelocity[i]*(dimTime - fabs(mAccelerationTime[i]));

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
