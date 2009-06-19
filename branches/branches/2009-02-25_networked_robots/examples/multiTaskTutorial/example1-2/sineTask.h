/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-    */
/* ex: set filetype=cpp softtabstop=4 shiftwidth=4 tabstop=4 cindent expandtab: */
/* $Id: sineTask.h 75 2009-02-24 16:47:20Z adeguet1 $ */

#ifndef _sineTask_h
#define _sineTask_h

// include for the whole cisstMultiTask library
#include <cisstMultiTask.h>

class sineTask: public mtsTaskPeriodic {
    // used to control the log level, 5 by default
    CMN_DECLARE_SERVICES(CMN_NO_DYNAMIC_CREATION, CMN_LOG_LOD_RUN_ERROR);
 protected:
    // data generated by the sine wave generator
    mtsDouble SineData;
    // this is the amplitude of the sine wave
    double SineAmplitude;

    void SetAmplitude(const mtsDouble &amp) { SineAmplitude = amp.Data; }
	
 public:
    // provide a name for the task and define the frequency (time
    // interval between calls to the periodic Run).  Also used to
    // populate the interface(s)
    sineTask(const std::string & taskName, double period);
    ~sineTask() {};
    // all four methods are pure virtual in mtsTask
    void Configure(const std::string & CMN_UNUSED(filename)) {};
    void Startup(void);    // set some initial values
    void Run(void);        // performed periodically
    void Cleanup(void) {}; // user defined cleanup

    void CommandVoidTest(void);
};

CMN_DECLARE_SERVICES_INSTANTIATION(sineTask);

#endif // _sineTask_h

/*
  Author(s):  Ankur Kapoor, Peter Kazanzides, Anton Deguet, Min Yang Jung
  Created on: 2004-04-30

  (C) Copyright 2004-2009 Johns Hopkins University (JHU), All Rights
  Reserved.

--- begin cisst license - do not edit ---

This software is provided "as is" under an open source license, with
no warranty.  The complete license can be found in license.txt and
http://www.cisst.org/cisst/license.txt.

--- end cisst license ---

*/
