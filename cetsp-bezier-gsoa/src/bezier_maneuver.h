/*
 * Filename: bezier_maneuver.h
 * Date:     2017/06/17 15:13
 * Author:   Petr Vana and Jan Faigl
 */

#ifndef __BEZIER_MANEUVER_H_
#define __BEZIER_MANEUVER_H_

#include "coords.h"

namespace bezier {

   static int BEZIER_SAMPLES = 200;
   static int BEZIER_SAMPLES_VISUAL = 20;

   /// ----------------------------------------------------------------------------
   /// Class Maneuver
   /// ----------------------------------------------------------------------------
   class Maneuver {
      public:
         State start, end;
         Coords p1, p2, p3, p4;
         Maneuver();
         Maneuver(State s1, State s2);

         double getLength() const;

         // t from [0, 1]
         State getStateRelative(double t) const;

         // t from [0, 1]
         Coords getPointRelative(double t) const;

         State getClosestState(const Coords &p) const;

         StateAtDistance getClosestStateAndDistance(const Coords &p) const;

      private:
         double calculateLength() const;
   };

} // end namespace bezier

std::ostream &operator<<(std::ostream &os, const bezier::Maneuver &d);

#endif

/* end of bezier_maneuver.h */
