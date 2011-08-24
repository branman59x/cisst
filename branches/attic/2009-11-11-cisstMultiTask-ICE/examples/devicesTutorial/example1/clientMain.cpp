/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-    */
/* ex: set filetype=cpp softtabstop=4 shiftwidth=4 tabstop=4 cindent expandtab: */
/* $Id: clientMain.cpp 671 2009-08-13 02:41:31Z adeguet1 $ */

#include <cisstCommon.h>
#include <cisstOSAbstraction.h>
#include <cisstMultiTask.h>

#include "displayTask.h"


int main(int argc, char * argv[])
{

    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " GlobalManagerIP ServerTaskIP" << std::endl;
        exit(-1);
    }

    std::string globalTaskManagerIP(argv[1]);
    std::string serverTaskIP(argv[2]);

    // log configuration
    cmnLogger::SetLoD(CMN_LOG_LOD_VERY_VERBOSE);
    cmnLogger::GetMultiplexer()->AddChannel(std::cout, CMN_LOG_LOD_VERY_VERBOSE);
    // add a log per thread
    osaThreadedLogFile threadedLog("example9Client");
    cmnLogger::GetMultiplexer()->AddChannel(threadedLog, CMN_LOG_LOD_VERY_VERBOSE);
    // specify a higher, more verbose log level for these classes
    cmnClassRegister::SetLoD("mtsTaskInterface", CMN_LOG_LOD_VERY_VERBOSE);
    cmnClassRegister::SetLoD("mtsTaskManager", CMN_LOG_LOD_VERY_VERBOSE);
    cmnClassRegister::SetLoD("displayTask", CMN_LOG_LOD_VERY_VERBOSE);

    // create our server task
    const double PeriodClient = 10 * cmn_ms; // in milliseconds
    displayTask * displayTaskObject = new displayTask("DISP", PeriodClient);

    // Get the TaskManager instance and set operation mode
    mtsTaskManager * taskManager = mtsTaskManager::GetInstance();
    taskManager->AddTask(displayTaskObject);        
    taskManager->SetGlobalTaskManagerIP(globalTaskManagerIP);
    taskManager->SetServerTaskIP(serverTaskIP);
    
    // Set the type of task manager either as a server or as a client.
    // mtsTaskManager::SetTaskManagerType() should be called before
    // mtsTaskManager::Connect()
    taskManager->SetTaskManagerType(mtsTaskManager::TASK_MANAGER_CLIENT);

    //
    // TODO: Hide this waiting routine inside mtsTaskManager using events or other things.
    //
    osaSleep(0.5 * cmn_s);

    // Connect the tasks across networks
    taskManager->Connect("DISP", "Robot", "Omni", "Omni1");
    taskManager->Connect("DISP", "Button1", "Omni", "Omni1Button1");
    taskManager->Connect("DISP", "Button2", "Omni", "Omni1Button2");

    // create the tasks, i.e. find the commands
    taskManager->CreateAll();
    // start the periodic Run
    taskManager->StartAll();

    while (1) {
        osaSleep(10 * cmn_ms);
    }
    
    // cleanup
    taskManager->KillAll();

    return 0;
}

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