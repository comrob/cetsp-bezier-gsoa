/*
 * File name: target.h
 * Date:      2016/02/05 10:20
 * Author:    Jan Faigl
 */

#ifndef __TARGET_H__
#define __TARGET_H__

#include <vector>

#include "coords.h"

namespace gsoa {

   struct STarget {
      const int label;
      const Coords coords;
      bool depot;

      int stepWinnerSelected;
      double ringDistance;

      STarget(const int label, const Coords &pt)
         : label(label), coords(pt), depot(false), stepWinnerSelected(-1)
      {
      } // point()
      ~STarget() {}

      inline double squared_distance(const Coords &c) const
      {
         return coords.squared_distance(c);
      }
   };

   typedef std::vector<STarget *> TargetPtrVector;

} // end namespace gsoa

#endif

/* end of target.h */
