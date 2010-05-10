/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-    */
/* ex: set filetype=cpp softtabstop=4 shiftwidth=4 tabstop=4 cindent expandtab: */

/*
  $Id: vctPlot2DBase.cpp 1238 2010-02-27 03:16:01Z auneri1 $

  Author(s):  Anton Deguet
  Created on: 2010-05-05

  (C) Copyright 2010 Johns Hopkins University (JHU), All Rights Reserved.

--- begin cisst license - do not edit ---

This software is provided "as is" under an open source license, with
no warranty.  The complete license can be found in license.txt and
http://www.cisst.org/cisst/license.txt.

--- end cisst license ---
*/

#include "vctPlot2DBase.h"


vctPlot2DBase::Trace::Trace(const std::string & name, size_t numberOfPoints):
    Active(true),
    Name(name),
    IndexFirst(0),
    IndexLast(0),
    Color(1.0, 1.0, 1.0),
    LineWidth(1.0)
{
    this->Data.SetSize(numberOfPoints);
}


void vctPlot2DBase::Trace::AddPoint(const vctDouble2 & point)
{
    // look where to store this point
    IndexLast = (IndexLast + 1) % this->Data.size();
    if (IndexFirst == IndexLast ) {
        IndexFirst = (IndexLast + 1) % this->Data.size();
    }
    Data.Element(IndexLast).Assign(point); 
}


void vctPlot2DBase::Trace::UpdateMinAndMax(vctDouble2 & min, vctDouble2 & max)
{
    size_t pointIndex;
    size_t numberOfPoints;
    double value;
    size_t nb = 0;
    numberOfPoints = this->Data.size();
    for (pointIndex = this->IndexFirst;
         pointIndex != this->IndexLast;
         pointIndex = (pointIndex + 1) % numberOfPoints) {
        nb++;
        value = this->Data.Element(pointIndex).X(); // X
        if (value < min.X()) {
            min.X() = value;
        } else if (value > max.X()) {
            max.X() = value;
        }
        value = this->Data.Element(pointIndex).Y(); // Y
        if (value < min.Y()) {
            min.Y() = value;
        } else if (value > max.Y()) {
            max.Y() = value;
        }
    }
}


vctPlot2DBase::vctPlot2DBase(void):
    NumberOfPoints(200),
    BackgroundColor(0.1, 0.1, 0.1)
{
}


bool vctPlot2DBase::AddTrace(const std::string & name, size_t & traceId)
{
    // check if the name already exists
    TracesIdType::const_iterator found = this->TracesId.find(name);
    const TracesIdType::const_iterator end = this->TracesId.end();
    if (found == end) {
        // new name
        traceId = Traces.size();
        Trace * newTrace = new Trace(name, this->NumberOfPoints);
        Traces.push_back(newTrace);
        TracesId[name] = traceId;
        return true;
    }
    return false;
}


void vctPlot2DBase::SetNumberOfPoints(size_t numberOfPoints)
{
    this->NumberOfPoints = numberOfPoints;
    if (this->NumberOfPoints < 3) {
        this->NumberOfPoints = 3;
    }
}


void vctPlot2DBase::AddPoint(size_t traceId, const vctDouble2 & point)
{
    this->Traces[traceId]->AddPoint(point);
} 


void vctPlot2DBase::SetBackgroundColor(const vctDouble3 & color)
{
    // should add check for [0,1]
    this->BackgroundColor = color;
}
