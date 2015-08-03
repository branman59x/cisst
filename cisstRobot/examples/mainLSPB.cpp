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
    const size_t dimension = 6;
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

    start.Assign(          0.0, 2.0, 3.0, 0.0, 2.0, 3.0);
    finish.Assign(        10.0, 12.0, 5.0, 10.0, 12.0, 15.0);
    maxVelocity.Assign(    4.0, 2.0, 5.0, 4.0, 2.0, 5.0);
    maxAcceleration.Assign(1.0, 1.0, 2.0, 1.0, 1.0, 2.0);
    initialVelocity.Assign(0.0, 1.0, 4.0, 5.0, -1.0, -6.0);

//    start[0] = 10;
//    finish[0] = 10;
//    maxVelocity[0] = 3;
//    maxAcceleration[0] = 2;
//    initialVelocity[0] = 0;

//    start.Assign(2.0, 0.0);
//    finish.Assign(10.0, 10.0);
//    maxVelocity.Assign(2.0, 4.0);
//    maxAcceleration.Assign(2.0, 1.0);
//    initialVelocity.Assign(1.0, 5.0);

//    start.Assign(0.00000, 0.486378, 0.0820529, 0.00000, 4.55546e-07, 0.00000, 0.00000);
//    finish.Assign(0.00000, 0.261799, 0.0462000, 0.00000, 4.55546e-07, 0.00000, 0.00000);
//    maxVelocity.Assign(1.57080, 1.57080, 0.200000, 1.57080, 1.57080, 1.57080, 1.57080);
//    maxAcceleration.Assign(1.57080, 1.57080, 0.200000, 1.57080, 1.57080, 1.57080, 1.57080);
//    initialVelocity.Assign(0.00230181, 0.486223, 0.114260, -0.000875546, -0.00111248, 0.000160613, -0.00222943);

    const double startTime = 2.0;

    robLSPB trajectory;
    trajectory.Set(start, finish,
                   maxVelocity, maxAcceleration, initialVelocity,
                   startTime, robLSPB::LSPB_DURATION); // default is LSPB_NONE

    const double duration = trajectory.Duration();
    const double extraPlotTime = 2.0;
    const double plotTime = extraPlotTime + duration + extraPlotTime;
    const size_t nbSteps = 20000;
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
