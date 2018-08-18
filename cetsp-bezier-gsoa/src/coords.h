/*
 * File name: coords.h
 * Date:      2013/10/13 09:23
 * Author:    Jan Faigl, Petr Vana
 */

#ifndef __COORDS_H__
#define __COORDS_H__

#include <cmath>
#include <iostream>
#include <vector>

/// ----------------------------------------------------------------------------
/// @brief Coords
/// ----------------------------------------------------------------------------
struct Coords {
   double x;
   double y;
   double z;

   Coords(Coords &c) : x(c.x), y(c.y), z(c.z) {}

   Coords(const Coords &c) : x(c.x), y(c.y), z(c.z) {}

   Coords() : x(0), y(0), z(0) {}

   Coords(double x, double y, double z) : x(x), y(y), z(z) {}

   Coords &operator=(const Coords &c)
   {
      if (this != &c) {
         x = c.x;
         y = c.y;
         z = c.z;
      }
      return *this;
   }

   inline Coords operator-(const Coords &r) const
   {
      return Coords(x - r.x, y - r.y, z - r.z);
   }

   inline Coords operator+(const Coords &right) const
   {
      return Coords(x + right.x, y + right.y, z + right.z);
   }

   inline Coords &operator+=(const Coords &rhs)
   {
      this->x += rhs.x;
      this->y += rhs.y;
      this->z += rhs.z;
      return *this;
   }

   inline Coords &operator-=(const Coords &rhs)
   {
      this->x -= rhs.x;
      this->y -= rhs.y;
      this->z -= rhs.z;
      return *this;
   }

   inline Coords operator*(const double b) const
   {
      return Coords(x * b, y * b, z * b);
   }

   inline Coords operator/(const double div) const
   {
      if (div != 0)
         return Coords(x / div, y / div, z / div);
      else
         return Coords();
   }

   inline Coords interpolate(const Coords &p, double &alpha) const
   {
      double beta = 1 - alpha;
      return Coords(beta * x + alpha * p.x, beta * y + alpha * p.y,
            beta * z + alpha * p.z);
   }

   inline double dotProduct(const Coords &b) const
   {
      return x * b.x + y * b.y + z * b.z;
   }

   inline bool operator==(const Coords &b) const
   {
      return (x == b.x) && (y == b.y) && (z == b.z);
   }

   inline double squared_distance(const Coords &c) const
   {
      return squared_distance(*this, c);
   }

   inline static double squared_distance(const Coords &c1, const Coords &c2)
   {
      double dx = c1.x - c2.x;
      double dy = c1.y - c2.y;
      double dz = c1.z - c2.z;
      return dx * dx + dy * dy + dz * dz;
   }

   inline double distance(const Coords &to) const
   {
      return std::sqrt(squared_distance(to));
   }

   inline Coords normalize() const
   {
      return *this / distance(Coords(0, 0, 0));
   }

   inline double length() const
   {
      return std::sqrt(squared_distance(Coords(0, 0, 0)));
   }

   inline Coords negate() const { return Coords(-x, -y, -z); }
};

typedef std::vector<Coords> CoordsVector;
typedef std::vector<CoordsVector> CoordsVectorVector;

struct Direction {
   Coords vector;
   // lengths of tangent for Bezier
   // control points are:
   // coords + a * vector
   // coords - b * vector
   // to the next
   double a;
   // from the previous
   double b;
};

class State {
   public:
      Coords coords;
      Direction direction;

      State() { direction.a = direction.b = NAN; };

      State(Coords coords, Coords vector, double a, double b) : coords(coords)
   {
      direction.vector = vector;
      direction.a = a;
      direction.b = b;
   };

      void init(Coords c, double angle, double a, double b)
      {
         coords = c;
         direction.vector.x = std::cos(angle);
         direction.vector.y = std::sin(angle);
         direction.vector.z = 0;
         direction.a = a;
         direction.b = b;
      }
};

typedef std::vector<State> StateVector;
typedef std::vector<StateVector> StateVectorVector;

typedef std::vector<double> DoubleVector;

struct StateAtDistance {
   State state;
   double distance;
};

typedef std::vector<StateAtDistance> StateAtDistanceVector;

inline std::istream &operator>>(std::istream &is, Coords &pt)
{
   return is >> pt.x >> pt.y >> pt.z;
}

inline std::ostream &operator<<(std::ostream &os, const Coords &pt)
{
   os << "(" << pt.x << "," << pt.y << "," << pt.z << ")";
   return os;
}

#endif

/* end of coords.h */
