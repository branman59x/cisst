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

#include <iostream>
#include <fstream>
#include <cisstCommon/cmnConstants.h>
#include <cisstVector/vctDynamicVectorTypes.h>
#include <cisstRobot/robLSPB.h>

int main(int CMN_UNUSED(argc), char ** CMN_UNUSED(argv))
{
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

    // set parameters
    //0 vel,  non0 vel, overshoot reg, max less than init, opposite direction, opposite direction and max less than init
//    start.Assign(          0.0, 2.0, 3.0);
//    finish.Assign(        10.0,12.0, 15.0);
//    maxVelocity.Assign(    4.0, 2.0, 5.0);
//    maxAcceleration.Assign(1.0, 1.0, 2.0);
//    initialVelocity.Assign(5.0, -1.0, -6.0);
    std::cout<<"I DID IT\n";

    start[0] = 3;
    finish[0] = 15;
    maxVelocity[0] = 5;
    maxAcceleration[0] = 2;
    initialVelocity[0] = -6;

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
    const char * logName = "robLSPB.txt";
    const char * logHeaderName = "robLSPB-header.txt";
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

    for (size_t i = 0; i < nbSteps; ++i) {
        double now = (startTime - extraPlotTime) + i * step;
        trajectory.Evaluate(now , position, velocity, acceleration);

        // csv file
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

    std::cout << "trajectory saved in " << logName << " (format description: " << logHeaderName << ")" << std::endl;

    return 0;
}
