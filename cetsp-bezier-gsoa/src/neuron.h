/*
 * File name: neuron.h
 * Date:      2016/02/05 10:21
 * Author:    Jan Faigl
 */

#ifndef __NEURON_H__
#define __NEURON_H__

#include "bezier_maneuver.h"

namespace gsoa {

   struct SNeuron {
      SNeuron *prev;
      SNeuron *next;

      Coords alternateGoal;

      int targetOnTourStep;
      int targetOnTour;
      int activityStep;

      State state;
      bezier::Maneuver maneuver;

      bool depot;

      SNeuron()
      {
         state.coords.x = state.coords.y = state.coords.z = state.direction.vector.y = state.direction.vector.z = NAN;
         state.direction.vector.x = NAN;
         clear();
         depot = false;
      }

      ~SNeuron() {}

      void setPoint(const Coords &pt)
      {
         state.coords = pt;
         clear();
      }

      void clear(void)
      {
         targetOnTour = -1;
         targetOnTourStep = -1;
         activityStep = -1;
      }

      void adapt(const Coords &pt, const double beta)
      {
         state.coords = state.coords + (pt - state.coords) * beta;
      }

      inline bool isInhibited(int step) { return targetOnTourStep >= step; }
   };

} // end namespace gsoa

#endif

/* end of neuron.h */
