/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-    */
/* ex: set filetype=cpp softtabstop=4 shiftwidth=4 tabstop=4 cindent expandtab: */
// $Id$

#ifndef _sineTask_h
#define _sineTask_h

#include <cisstMultiTask.h>

class sineTask: public mtsTaskPeriodic {
    CMN_DECLARE_SERVICES(CMN_NO_DYNAMIC_CREATION, CMN_LOG_LOD_RUN_ERROR);
    
 protected:
    // data generated by the sine wave generator
    mtsDouble SineData;
    // this is the amplitude of the sine wave
    double SineAmplitude;
    // trigger value, i.e. when the data generated reaches this value
    // an event is generated.
    double TriggerValue;
    // methods to set SineAmplitude and TriggerValue
    void SetAmplitude(const mtsDouble & amp) { SineAmplitude = amp.Data; }
    void SetTrigger(const mtsDouble & val) { TriggerValue = val.Data; }
    // method used to reset the trigger, we will add a command to the
    // interface for this method
    void ResetTrigger(void);
    // function bound to the command used to send the event, one could
    // use a command instead but this is somewhat more convenient
    mtsFunctionWrite TriggerEvent;
    // internal flag
    bool TriggerEnabled;
	
 public:
    // as in previous example
    sineTask(const std::string & taskName, double period);
    ~sineTask() {};
    void Configure(const std::string & CMN_UNUSED(filename) = "") {};
    void Startup(void);
    void Run(void);
    void Cleanup(void) {};
};

CMN_DECLARE_SERVICES_INSTANTIATION(sineTask);

#endif // _sineTask_h

/*
  Author(s):  Ankur Kapoor, Peter Kazanzides, Anton Deguet
  Created on: 2004-04-30

  (C) Copyright 2004-2008 Johns Hopkins University (JHU), All Rights
  Reserved.

--- begin cisst license - do not edit ---

This software is provided "as is" under an open source license, with
no warranty.  The complete license can be found in license.txt and
http://www.cisst.org/cisst/license.txt.

--- end cisst license ---

*/
