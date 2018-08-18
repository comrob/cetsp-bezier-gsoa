/*
 * Filename: velocity_profile.cc
 * Date:     2017/06/22
 * Author:   Petr Vana
 */

#include <cmath>

#include <crl/logging.h>

#include "bezier_maneuver.h"
#include "velocity_profile.h"

using namespace bezier;
using crl::logger;

/// - public method ------------------------------------------------------------
InterpolationPointVector VelocityProfile::interpolateOpenLoop(const StateVector &states, const double endSpeed) const
{
   InterpolationPointVector interpolation;
   if (states.size() < 2) {
      WARN("No input trajectory - less than 2 states!");
   }
   if (states.size() < 1) {
      ERROR("No input trajectory - no state!");
   }
   // 1 - interpolate points
   const int count = states.size();
   for (int i = 0; i < count - 1; i++) {
      const State &s1 = states[i];
      const State &s2 = states[i + 1];
      const Maneuver maneuver = bezier::Maneuver(s1, s2);
      for (int j = 0; j <= BEZIER_SAMPLES; j++) {
         double t = ((double)j) / BEZIER_SAMPLES;
         InterpolationPoint p;
         p.coords = maneuver.getPointRelative(t);
         // p.accel = 0;
         p.maxAccel = horizontalAccel;
         p.speed = horizontalSpeed;
         p.curvature = 0;
         p.time = NAN;
         p.maxSpeed = horizontalSpeed;
         interpolation.push_back(p);
      }
   }

   // set endSpeed velocity
   interpolation.front().speed = endSpeed;
   interpolation.front().time = 0;
   interpolation.back().speed = endSpeed;

   // 2 - remove duplicity points
   InterpolationPointVector copy = interpolation;
   interpolation.clear();
   Coords last;
   bool first = true;
   for (const InterpolationPoint &p : copy) {
      if (!first) {
         if (p.coords.distance(last) > horizontalSpeed * 1e-5) {
            interpolation.push_back(p);
         }
      } else {
         interpolation.push_back(p);
      }
      last = p.coords;
      first = false;
   }

   // 3 - compute curvature and maximum speed
   int ints = interpolation.size();
   for (int i = 1; i < ints - 1; i++) {
      const InterpolationPoint &prev = interpolation[i - 1];
      const InterpolationPoint &next = interpolation[i + 1];
      InterpolationPoint &cur = interpolation[i];

      Coords d1 = cur.coords - prev.coords;
      Coords d2 = next.coords - cur.coords;

      double cosinus = d1.normalize().dotProduct(d2.normalize());
      double alpha = std::acos(cosinus);

      double radius = (d1.length() + d2.length()) / 2 / alpha;

      cur.curvature = 1 / radius;
      cur.maxSpeed =
         std::fmin(cur.speed, std::sqrt(horizontalAccel / cur.curvature));
   }

   // 4 - compute speed forward limitation
   for (int i = 0; i < ints - 1; i++) {
      InterpolationPoint &act = interpolation[i];
      InterpolationPoint &next = interpolation[i + 1];
      limitSpeed(act, next);
   }

   // 5 - compute speed backward limitation
   for (int i = ints - 1; i; i--) {
      InterpolationPoint &act = interpolation[i];
      InterpolationPoint &next = interpolation[i - 1];
      limitSpeed(act, next);
   }

   // 6 - compute time
   for (int i = 1; i < ints; i++) {
      InterpolationPoint &act = interpolation[i];
      InterpolationPoint &prev = interpolation[i - 1];
      double dst = act.coords.distance(prev.coords);
      act.time = prev.time + dst / (act.speed + prev.speed) * 2;
   }
   return interpolation;
}

/// - public method ------------------------------------------------------------
InterpolationPointVector VelocityProfile::interpolateUniformTime(double timeStep, const InterpolationPointVector &interpolation) const
{
   InterpolationPointVector result;
   result.push_back(interpolation.front());
   bool more = true;
   InterpolationPointVector::const_iterator prev = interpolation.begin();
   InterpolationPointVector::const_iterator next = prev++;
   for (double time = timeStep; more; time += timeStep) {
      while (next.base()->time < time) {
         prev = next;
         if (++next == interpolation.end()) {
            more = false;
            break;
         }
      }
      if (!more)
         break;

      double t1 = prev->time;
      double t2 = next->time;
      double alpha = (t2 - time) / (t2 - t1);
      InterpolationPoint point = *(next.base());
      point.time = time;
      point.coords = prev->coords * alpha + next->coords * (1 - alpha);
      point.speed = prev->speed * alpha + next->speed * (1 - alpha);
      point.curvature = prev->curvature * alpha + next->curvature * (1 - alpha);
      point.maxAccel = prev->maxAccel * alpha + next->maxAccel * (1 - alpha);

      result.push_back(point);
   }

   return result;
}

/// - private method -----------------------------------------------------------
void VelocityProfile::limitSpeed(InterpolationPoint &from, InterpolationPoint &to) const
{
   double dst = from.coords.distance(to.coords);
   // centrifugal acceleration
   double cent_accel = from.speed * from.speed * from.curvature;
   double max_tangent_accel_squared = horizontalAccel * horizontalAccel - cent_accel * cent_accel;
   double max_tangent_accel = max_tangent_accel_squared > 0 ? std::sqrt(max_tangent_accel_squared) : 0;
   // limit vertical acceleration
   double verAccel = verticalAccel;
   const Coords dir = to.coords - from.coords;
   double dirxy = std::sqrt(dir.x * dir.x + dir.y * dir.y);
   double dirz = std::fabs(dir.z);
   double dirxyz = dir.length();
   // compute maximum acceleration (verAccel) to meet vertical limitation
   if (dirz > 1e-5 * dirxy) {
      verAccel = verticalAccel / dirz * dirxyz;
   }

   // choose smaller acceleration
   double accel = std::min(verAccel, max_tangent_accel);

   double D = from.speed * from.speed + 2 * accel * dst;
   // compute maximum speed at the end of this segment
   double next_speed = from.speed + (std::sqrt(D) - from.speed);

   // limit vertical speed
   double ver_speed = std::max(verticalSpeed, horizontalSpeed);
   if (dirz > 1e-5 * dirxy) {
      ver_speed = verticalSpeed / dirz * dirxyz; // TODO - fix, compute by triangle
   }
   // limit the speed at the end of the segment by the given limits
   to.speed = std::min(std::min(to.speed, next_speed), std::min(to.maxSpeed, ver_speed));
}

std::ostream &operator<<(std::ostream &os, const bezier::InterpolationPoint &pt)
{
   os << "Point " << pt.coords << " time " << pt.time << " speed " << pt.speed;
   return os;
}

/* end of velocity_profile.cc */
