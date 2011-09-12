/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-    */
/* ex: set filetype=cpp softtabstop=4 shiftwidth=4 tabstop=4 cindent expandtab: */

/*
  Author(s):	Rajesh Kumar, Anton Deguet
  Created on:	2008-03-12

  (C) Copyright 2008 Johns Hopkins University (JHU), All Rights
  Reserved.

--- begin cisst license - do not edit ---

This software is provided "as is" under an open source license, with
no warranty.  The complete license can be found in license.txt and
http://www.cisst.org/cisst/license.txt.

--- end cisst license ---
*/


#ifndef _prmMotionBase_h
#define _prmMotionBase_h

// basic includes
#include <cisstCommon/cmnClassRegister.h>

#include <cisstParameterTypes/prmTypes.h>
#include <cisstParameterTypes/prmMacros.h>

// Always include last
#include <cisstParameterTypes/prmExport.h>

/*!
  base class for motion arguments 
*/
class CISST_EXPORT prmMotionBase: public cmnGenericObject
{
 public:
    /* default constructor*/
    inline prmMotionBase():
        BlockingFlagMember(NO_WAIT),
        BlendingFactorMember(0),
        IsPreemptableMember(true),
        IsCoordinatedMember(true),
        IsGoalOnlyMember(false)
    {}

    /*! Set all parameters
      \param 
    */
    inline prmMotionBase(const prmBlocking & blockingFlag, 
                         const prmBlending & blendingFactor,
                         const prmTime & timeLimit,
                         const bool & isPreemptable,
                         const bool & isCoordinated,
                         const bool & isGoalOnly
                         ):
        BlockingFlagMember(blockingFlag), 
        BlendingFactorMember(blendingFactor),
        TimeLimitMember(timeLimit),
        IsPreemptableMember(isPreemptable),
        IsCoordinatedMember(isCoordinated),
        IsGoalOnlyMember(isGoalOnly)
    {}
    
    /*! destructor
     */
    virtual ~prmMotionBase();
    
    /*! Set all parameters from the base class
      \param 
    */
    inline void SetBaseParameter(const prmBlocking & blockingFlag, 
                                 const prmBlending & blendingFactor, 
                                 const prmTime & timeLimit,
                                 const bool & isPreemptable,
                                 const bool & isCoordinated,
                                 const bool & isGoalOnly)
    { 
        this->BlockingFlagMember = blockingFlag; 
        this->BlendingFactorMember = blendingFactor; 
        this->TimeLimitMember = timeLimit; 
        this->IsPreemptableMember = isPreemptable; 
        this->IsCoordinatedMember = isCoordinated; 
        this->IsGoalOnlyMember = isGoalOnly; 
    }
   
    /*! Set and get methods for blocking flag. */
    //@{
    PRM_DECLARE_MEMBER_AND_ACCESSORS(prmBlocking, BlockingFlag);
    //@}

    /*! Set and Get method for the blending parameter.  Whether or not
      to blend motions, and the tolerance (to be implemented later)
    */
    //@{
    PRM_DECLARE_MEMBER_AND_ACCESSORS(prmBlending, BlendingFactor);
    //@}

    /*! Set and Get methods for TimeLimit.  Describes time allowed or
      requested for a motion.
    */
    //@{
    PRM_DECLARE_MEMBER_AND_ACCESSORS(prmTime, TimeLimit);
    //@}

    /*! Set and Get methods for Preemptable parameter.  Default is
      true.*/
    //@{
    PRM_DECLARE_MEMBER_AND_ACCESSORS(bool, IsPreemptable);
    //@}

    /*! Set and Get methods for Coordinated parameter.  Force to sync
      motions between joints or rotation/translation to start and end
      at the same time by lowering speed on some joints. */
    //@{
    PRM_DECLARE_MEMBER_AND_ACCESSORS(bool, IsCoordinated);
    //@}

    /*! Set and Get methods for IsGoalOnly.
        if True Use current device settings for all except for goal .
        Eg. Accelerations,Velcities are chosen by the robot, 
        goal and motion base should be set by the user.
        False by default.
    */
    //@{
    PRM_DECLARE_MEMBER_AND_ACCESSORS(bool, IsGoalOnly);
    //@}

};


#endif  // _prmMotionBase_h