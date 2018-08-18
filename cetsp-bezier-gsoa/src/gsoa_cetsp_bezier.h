/*
 * File name: gsoa_cetsp_bezier.h
 * Date:      2016/02/05 08:59
 * Author:    Jan Faigl
 */

#ifndef __GSOA_CETSP_BEZIER_H__
#define __GSOA_CETSP_BEZIER_H__

#include <memory>
#include <thread>
#include <vector>

#include <crl/alg/algorithm.h>
#include <crl/config.h>

#include "coords.h"
#include "neuron.h"
#include "target.h"
#include "ring.h"

#include "velocity_profile.h"

namespace gsoa {

   class CGSOACETSPBezier : public crl::CAlgorithm
   {
      typedef crl::CAlgorithm Base;
      typedef std::vector<int> IntVector;
      typedef std::vector<SNeuron *> NeuronPtrVector;
      typedef enum { SOMHOM_PARAM, ICANN_PARAM } TParamChange;
      typedef enum { RING_LENGTH_EUCLID, RING_LENGTH_TIME } TRingLength;

      struct SWinner {
         SNeuron *winner;
         SNeuron *prev;
         SNeuron *next;
         SRing *ring;
         int ringIDX;
      };

      typedef std::vector<SRing *> RingPtrVector;

      public:
      static crl::CConfig &getConfig(crl::CConfig &config);

      static std::string getName(void) { return "cetsp-bezier-gsoa"; }

      static TParamChange getParamChange(const std::string &str);
      static std::string getParamChange(TParamChange param);
      static TRingLength getRingLength(const std::string &str);
      static std::string getRingLength(TRingLength ringLength);

      CGSOACETSPBezier(crl::CConfig &config, const std::string &problemFile);
      ~CGSOACETSPBezier();

      std::string getRevision(void);
      std::string getMethod(void) { return getName(); }
      void setTimers(const crl::CTimerN &load, const crl::CTimerN &init);
      std::string getVersion(void);

      void solve(void);

      protected:
      //# methods from CAlgorithm
      void load(void);
      void initialize(void);
      void after_init(void);
      void iterate(int iter);
      void save(void);
      void release(void);
      void visualize(void);
      void defineResultLog(void);
      void fillResultRecord(int iter, double length, double time, int steps);

      void allocate_neurons(void);   // allocate neurons
      void deallocate_neurons(void); // release neurons
      void initialize_neurons(void); // place neurons around the first target

      double refine(int step, double errorMax);

      private:
      SNeuron *selectWinner(int step, STarget *target, bool allrings, int ringIDX, double &errorToTarget, SWinner &winner);
      void adapt(SRing *ring, SNeuron *neuron, const SWinner &winner, const int step);

      void drawRing(int step);

      void drawPath(void);

      void savePic(int step, bool detail = false, const std::string &dir_suffix = "");
      void releaseTargets(void);

      void ring_regenerate(int step);

      void load_targets_coords(const std::string &filename);

      bool isSolutionAdmissibleCoverageSimple(int step, DoubleVector &lengths, double &solutionError, StateVectorVector &solutions);

      CoordsVector &get_ring_path(int step, const SRing *const ring, CoordsVector &path) const;

      // determine heading of the center point
      void determineCenterHeading(Coords p1, State &state, Coords p3);

      double calculateTime(const State &s1, const State &s2, const State &s3);
      void optimizeDirection(const State &p1, State &state, const State &p3);

      void updateHeadings(int step, SRing *ring, bool fixWinners);

      bool getWinnersDTPSolution(int step, DoubleVector &solutionLength, StateVectorVector &solutions);

      double getRingWinnerPenalty(const int ringIDX);
      double getRingLength(const SRing *const ring) const;
      double getRingLengthEuclid(const SRing *const ring) const;
      double getRingLengthTime(const SRing *const ring) const;

      /// ----------------------------------------------------------------------------
      /// Variables
      /// ----------------------------------------------------------------------------
      protected:
      IntVector permutation;

      double mi;

      double G; // variable
      int n;    // number targets

      // solution representation
      SRing **rings;

      crl::CTimerN loadTimer;
      crl::CTimerN initTimer;
      bool savePicEnabled;
      const TParamChange PARAM_CHANGE;
      const TRingLength RING_LENGTH;
      double RADIUS;
      const double GAIN_DECREASING_RATE;
      const double NEIGHBORHOOD_FACTOR;
      const int NEURON_NUMBER_HEADINGS;
      const bool ALTERNATE_TARGET;
      const double ICANN_GAIN_DECREASING_RATE;
      const bool RING_REGENERATE_MIN_WINNER_DISTANCE;
      const bool RING_REGENERATE_PRESERVE_ACTIVE_NEURONS;
      const bool FORCE_ZERO_ERROR_QUIT;
      const double COVER_DISTANCE;
      const double COVER_DISTANCE_SQ;

      const double BEZIER_OPTIMIZATION_STEP;
      const int BEZIER_OPTIMIZATION_MAX;
      const int BEZIER_OPTIMIZATION_CYCLES;
      const double MAX_HORIZONTAL_SPEED;
      const double MAX_HORIZONTAL_ACCEL;
      const double MAX_VERTICAL_SPEED;
      const double MAX_VERTICAL_ACCEL;
      const double SAMPLING_STEP;

      const bool SAVE_RESULTS;
      const bool SAVE_INFO;
      const bool SAVE_SETTINGS;

      const int NRINGS;

      TargetPtrVector targets;
      const double BORDER;
      bool drawNeurons;
      bool drawAlternateGoals;
      bool drawWinners;
      bool drawTourRepresentedByRing;

      crl::gui::CShape shapeRing;
      crl::gui::CShape shapeNeuron;
      crl::gui::CShape shapeWinner;
      crl::gui::CShape shapeAlternateGoal;
      crl::gui::CShape shapeTourRepresentedByRing;
      bool drawRingIter;
      bool drawRingEnabled;

      DoubleVector ringLengths;
      double avgLength;
      StateVectorVector finalPaths;

      IntVector depots;

      // velocity profile with speed settings
      bezier::VelocityProfile velocityProfile;
   };

} // end namespace gsoa

#endif

/* end of gsoa_cetsp_bezier.h */
