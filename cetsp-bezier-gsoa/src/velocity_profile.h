/*
 * Filename: velocity_profile.h
 * Date:     2017/06/22
 * Author:   Petr Vana and Jan Faigl
 */

#ifndef __VELOCITY_PROFILE_H__
#define __VELOCITY_PROFILE_H__

#include "coords.h"

namespace bezier {

   struct InterpolationPoint {
      Coords coords;
      double curvature;
      double maxAccel;
      double speed;
      double maxSpeed;
      double time;
   };

   typedef std::vector<InterpolationPoint> InterpolationPointVector;

   class VelocityProfile {
      public:
         double horizontalSpeed;
         double horizontalAccel;
         double verticalSpeed;
         double verticalAccel;

         // WARNING! - this is open loop function
         // to create close loop copy first state at the end of the list
         // create velocity profile for given list of states
         InterpolationPointVector interpolateOpenLoop(const StateVector &states, const double endSpeed = 0) const;

         // interpolate given velocity profile for uniform time step 
         InterpolationPointVector interpolateUniformTime(double timeStep, const InterpolationPointVector &interpolation) const;

      private:
         void limitSpeed(InterpolationPoint &from, InterpolationPoint &to) const;
   };

} // end namespace bezier

std::ostream &operator<<(std::ostream &os, const bezier::InterpolationPoint &pt);

#endif

/* end of velocity_profile.h */
