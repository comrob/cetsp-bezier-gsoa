/*
 * File name: bezier_maneuver.cc
 * Date:      2017/06/17 15:23
 * Author:    Petr Vana and Jan Faigl
 */

#include <crl/logging.h>

#include "bezier_maneuver.h"

using crl::logger;

using namespace bezier; 

/// - constructor --------------------------------------------------------------
Maneuver::Maneuver() {}

/// - constructor --------------------------------------------------------------
Maneuver::Maneuver(State s1, State s2) : start(s1), end(s2)
{
   p1 = s1.coords;
   p2 = s1.coords + s1.direction.vector.normalize() * (s1.direction.a);
   p3 = s2.coords - s2.direction.vector.normalize() * (s2.direction.b);
   p4 = s2.coords;
}

/// - public method ------------------------------------------------------------
double Maneuver::getLength() const { return calculateLength(); }

/// - public method ------------------------------------------------------------
State Maneuver::getStateRelative(double t) const
{ // t from [0, 1]
   State state;
   state.coords = getPointRelative(t);
   const Coords next = getPointRelative(t + 1e-5);
   state.direction.vector = (next - state.coords);
   return state;
}

/// - public method ------------------------------------------------------------
Coords Maneuver::getPointRelative(double t) const
{ // t from [0, 1]
   double t2 = 1 - t;

   Coords m1 = p1 * t2 + p2 * t;
   Coords m2 = p2 * t2 + p3 * t;
   Coords m3 = p3 * t2 + p4 * t;

   Coords n1 = m1 * t2 + m2 * t;
   Coords n2 = m2 * t2 + m3 * t;

   Coords ret = n1 * t2 + n2 * t;

   return ret;
}

/// - public method ------------------------------------------------------------
State Maneuver::getClosestState(const Coords &p) const
{
   return getClosestStateAndDistance(p).state;
}

/// - public method ------------------------------------------------------------
StateAtDistance Maneuver::getClosestStateAndDistance(const Coords &p) const
{
   StateAtDistance closest;

   double minDist = p.distance(start.coords);
   double minT = 0;

   const int C = BEZIER_SAMPLES;
   for (int i = 0; i <= C; i++) {
      double t = ((double)i) / C;
      const Coords act = getPointRelative(t);
      double actDst = act.distance(p);
      if (actDst < minDist) {
         minDist = actDst;
         minT = t;
      }
   }

   closest.state = getStateRelative(minT);
   closest.distance = minT;

   return closest;
}

/// - private method -----------------------------------------------------------
double Maneuver::calculateLength() const
{
   double length = 0;
   Coords last = getPointRelative(0);
   Coords llast = last;
   const int C = BEZIER_SAMPLES;
   for (int i = 0; i <= C; i++) {
      double t = ((double)i) / C;
      const Coords act = getPointRelative(t);
      length += act.distance(last);
      llast = last;
      last = act;
   }
   return length;
}

/* end of bezier_maneuver.cc */
