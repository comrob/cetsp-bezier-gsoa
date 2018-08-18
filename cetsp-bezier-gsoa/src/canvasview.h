/*
 * File name: canvasview_coords.h
 * Date:      2013/10/13 09:26
 * Author:    Jan Faigl
 */

#ifndef __CANVASVIEW_H__
#define __CANVASVIEW_H__

#include <crl/gui/canvas.h>
#include <crl/gui/colors.h>
#include <crl/gui/gui.h>
#include <crl/gui/renderer.h>

#include "bezier_maneuver.h"
#include "coords.h"

/// ----------------------------------------------------------------------------
inline crl::gui::CCanvasBase &operator<<(crl::gui::CCanvasBase &canvas,
      const bezier::Maneuver &maneuver)
{
   Coords last = maneuver.getPointRelative(0);
   const int count = bezier::BEZIER_SAMPLES_VISUAL;
   for (int i = 0; i <= count; i++) {
      const Coords act = maneuver.getPointRelative(((double)i) / count);
      canvas << crl::gui::canvas::LINE;
      double t = std::fmax(std::fmin((act.z - 5) / 10, 1.0), 0);
      crl::gui::SColor color;
      color.set(t, 1 - t, 0, 1.0);
      canvas << color << last.x << last.y << act.x << act.y;
      last = act;
   }

   canvas << crl::gui::canvas::END;
   return canvas;
}

/// ----------------------------------------------------------------------------
inline crl::gui::CCanvasBase & operator<<(crl::gui::CCanvasBase &canvas, const std::vector<bezier::Maneuver> &path)
{
   for (const bezier::Maneuver &maneuver : path) {
      canvas << maneuver;
   }
   return canvas;
}

/// ----------------------------------------------------------------------------
inline crl::gui::CCanvasBase &operator<<(crl::gui::CCanvasBase &canvas,
      const std::vector<State> &states)
{
   for (const State &state : states) {
      canvas << state.coords.x << state.coords.y;
   }
   return canvas;
}

/// ----------------------------------------------------------------------------
inline crl::gui::CCanvasBase &operator<<(crl::gui::CCanvasBase &canvas,
      const Coords &coords)
{
   canvas << coords.x << coords.y;
   return canvas;
}

/// ----------------------------------------------------------------------------
inline crl::gui::CCanvasBase &operator<<(crl::gui::CCanvasBase &canvas,
      const std::vector<Coords> &points)
{
   for (const Coords &pt : points) {
      canvas << pt;
   }
   return canvas;
}

#endif

/* end of canvasview_coords.h */
