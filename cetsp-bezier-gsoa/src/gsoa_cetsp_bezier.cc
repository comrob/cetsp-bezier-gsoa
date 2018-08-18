/*
 * File name: gsoa_cetsp_bezier.cc
 * Date:      2016/02/05 08:59
 * Author:    Jan Faigl, Petr Vana
 */

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <thread>

#include <crl/assert.h>
#include <crl/file_utils.h>
#include <crl/logging.h>
#include <crl/perf_timer.h>
#include <crl/stringconversions.h>
#include <crl/timerN.h>

#include <crl/gui/colormap.h>
#include <crl/gui/shapes.h>

#include <crl/alg/text_result_log.h>

#include "canvasview.h"

#include "gsoa_cetsp_bezier.h"
#include "velocity_profile.h"

using namespace gsoa;
using namespace crl;
using namespace crl::gui;

#define foreach(x, y) for (x : y)

static const double NEURON_COORDS_IDENTITY2 = 1e-5 * 1e-5;

typedef std::vector<int> IntVector;

/// ----------------------------------------------------------------------------
template <class T> inline double getMax(const T &array)
{
   double ret = array.empty() ? -1.0 : array[0];
   for (int i = 1; i < array.size(); ++i) {
      if (ret < array[i]) {
         ret = array[i];
      }
   }
   return ret;
}

/// ----------------------------------------------------------------------------
std::string getString(const DoubleVector &a, const int n)
{
   std::stringstream ss;
   if (a.size() == n) {
      for (int r = 0; r < n; ++r) {
         if (r != 0) {
            ss << ", ";
         }
         ss << a[r];
      }
   } else {
      ss << "-";
   }
   return ss.str();
}

/// ----------------------------------------------------------------------------
class CompareIntersection
{
   public:
      bool operator()(const StateAtDistance &a, const StateAtDistance &b) const
      {
         return a.distance < b.distance;
      }
};

/// ----------------------------------------------------------------------------
void createPermutation(int number, IntVector &permutation)
{
   permutation.clear();
   for (int i = 0; i < number; ++i) {
      permutation.push_back(i);
   }
}

/// ----------------------------------------------------------------------------
void permute(IntVector &permutation)
{
   int k, tmp;
   crl::CRandom::randomize();
   for (int i = permutation.size(); i > 0; --i) {
      k = crl::CRandom::random() % i;
      tmp = permutation[i - 1];
      permutation[i - 1] = permutation[k];
      permutation[k] = tmp;
   }
}

void checkRing(SRing *ring)
{
   // check values
   SNeuron *cur = ring->neurons;
   for (int i = 0; i < ring->nbrNeurons; ++i) {
      if (cur->next && cur->targetOnTourStep > 2) {
         if (cur->state.direction.vector.length() < 0.5) {
            ERROR("Some values are not initialized yet! (7)");
            // exit(30);
         }
         if (std::isnan(cur->state.direction.a)) {
            ERROR("Some values are NAN! (11)");
            // exit(30);
         }
         cur = cur->next;
      }
   }
}

/// ----------------------------------------------------------------------------
/// Class CGSOACETSPBezier
/// ----------------------------------------------------------------------------

/// - static method ------------------------------------------------------------
crl::CConfig &CGSOACETSPBezier::getConfig(crl::CConfig &config)
{
   // basic config
   config.add<std::string>("output", "output directory to store particular results and outputs", "./");
   config.add<std::string>("results", "result log file, it will be placed in output directory", "results.log"); // TODO .db for sqlite, .txt text logs
   config.add<std::string>("info","information file, it will be placed in particular experiment directory", "info.txt");
   config.add<std::string>("settings", "store configurations in boost::program_options config file format", "settings.txt");
   config.add<std::string>("result-path", "file name for the final found path (ring) as sequence of points", "path");
   config.add<std::string>("result-ext", "File name extension for the final found path (ring) as sequence of points", ".txt");
   config.add<std::string>("result-canvas-output", "file name for final canvas output (eps,png,pdf,svg) are supported");  
   config.add<std::string>("result-canvas-suffixes", "comman separated list of exptensis e.g. png,pdf. If specified several result images are created, with particular suffix to the resultCanvasOutput");
   config.add<std::string>("name", "name used in result log as user identification if not set a default values (cities) is used"); 
   config.add<int>("iteration", "set particular interation, otherwise all interations (batch) are performed", -1);
   config.add<int>("batch", "number of iterations from 0 to batch (not included) ", -1);
   config.add<bool>("continue", "in batch mode, partial results are loaded and checked, only missing iterations are computed ", false);
   config.add<bool>("save-results", "disable/enable save results,configs and so on", true);
   config.add<bool>("save-info", "disable/enable save info", true);
   config.add<bool>("save-settings", "disable/enable save settings", true);
   config.add<bool>("save-visual", "disable/enable save visualization results as the canvas", true);
   config.add<bool>("verbose-result-log", "disable/enable printing results log into logger", false);
   // end basic config
   //
   config.add<std::string>("param-change", "somhom|icann", "somhom");
   config.add<std::string>("ring-length", "ring_euclid|ring_time", "ring_eulid");
   config.add<double>("learning-rate", "neuron adaptation parametr in activation function mi*exp(d)", 0.6); // 0.6
   config.add<double>("number-neurons-multiplication", "multiplication parametr to get number of neurons as multiplication of number of cities", 2.5);
   config.add<double>("gain-decreasing-rate", "Decrasing rate, higher value means faster convergence", 1e-4);
   config.add<double>("neighborhood-factor", "factor to compute number of neighborhood neurons as division of number of neurons, neighborhood are considered on both sides (left, right) of winning neuron so d*2 neurons are moved ", 5);
   config.add<int>("termination-max-steps", "Stop adaptation if number of steps is large than given", 180);
   config.add<double>("termination-error", "Stop adaptation when the current error is less than this value", 0.001);
   config.add<bool>("best-path", "Enable/disable considering best path found during the evaluation", false);
   config.add<int>("neuron-number-headings", "Number of alternative headings per neuron for each direction", 10);
   config.add<bool>("alternate-target", "Enable/disable alternate target, the winner is not set to the presented target", true);
   config.add<int>("vehicles", "Number of vehicles", 1);
   config.add<std::string>("start-targets", "A list of target labels (starting from 0) that denote the starting locations (separated by ','). Must match with number of vehicles", "0");
   config.add<double>("icann-gain-decreasing-rate", "Decrasing reate for the icann param change, higher value means faster convergence", 1e-4);
   config.add<bool>("ring-regenerate-min-distance", "If enabled only neurons that are a bit far from each other are preserved", true);
   config.add<bool>("ring-regenerate-preserve-active", "If enabled also active neurons are preserved in the ring regeneration", true);
   config.add<bool>("force-zero-error-quit", "If enabled, the adaptation is termined when winners fit the targets and the solution is not improving", true);
   config.add<double>("target-cover-distance", "If > 0 and dubins_coverage select winner, the target is considered covered if the current ring is within the given distance from the target", -1.0);
   //
   config.add<std::string>("pic-dir", "relative directory in result directory to store pictures from each iteration");
   config.add<std::string>("pic-ext", "extension of pic, eps, png, pdf, svg (supported by particular gui renderer rendered", "png");
   config.add<bool>("save-pic", "enable/disable saving pictures (after each refine)");
   //
   config.add<bool>("draw-targets", "Enable/Disable drawing targets", true);
   config.add<bool>("draw-ring", "Enable/Disable drawing ring in the final shoot", true);
   config.add<bool>("draw-path", "Enable/Disable drawing ring in the final shoot", true);
   config.add<bool>("draw-alternate-goals", "Enable/Disable drawing alternate goals in drawing ring", false);
   config.add<double>("canvas-border", "Free space around the canvas", 10);
   config.add<bool>("draw-cover-circle", "Enable/Disable drawing cover circle around target", true); 
   config.add<std::string>("draw-shape-target", "Shape of the target", Shape::CITY());
   config.add<std::string>("draw-shape-neuron", "Shape of the neurons", Shape::NEURON());
   config.add<std::string>("draw-shape-winner", "Shape of the winners", Shape::NEURON());
   config.add<std::string>("draw-shape-depot", "Shape of the depot", Shape::DEPOT());
   config.add<std::string>("draw-shape-path", "Shape of the path", Shape::RED_LINE());
   config.add<std::string>("draw-shape-ring", "Shape of the ring", Shape::GREEN_BOLD_LINE());
   config.add<std::string>("draw-shape-tour-represented-by-ring", Shape::BLACK_BOLD_LINE());
   config.add<std::string>("draw-shape-alternate-goal", Shape::MAP_VERTEX());
   config.add<std::string>("draw-shape-cover-circle", Shape::REGION_GOLD());

   config.add<bool>("draw-neurons", "enable/disable drawing neurons", false);
   config.add<bool>("draw-winners", "enable/disable drawing winner using a different shape", false);
   config.add<bool>("draw-depots", "enable/disable drawing depots using a different shape", false);
   config.add<bool>("draw-ring-iter", "enable/disable drawing ring at each iteration", false);
   config.add<bool>("draw-path-vertices", "enable/disable drawing path vertices(nodes)", true);
   config.add<bool>("draw-tour-represented-by-ring", "enable/disable drawing tour represented by ring", false);
   //
   config.add<std::string>("input-path", "If this input file is given, a solution is computed from it considering the given radius", "");
   //
   config.add<double>("bezier-optimization-step", "Step size for local optimization of the Bezier curve", 0.01);
   config.add<int>("bezier-optimization-max", "Maximum count of optimization iteration of Bezier curve", 100);
   config.add<int>("bezier-optimization-cycles", "Maximum number of cycles of optimization iteration of Bezier curve", 3);
   //
   config.add<int>("bezier-samples", "Numver of samples used for numberic solution of the Bezier curve", 200);
   config.add<int>("bezier-samples-visual", "Number of samples used for the visualization of the Bezier curve", 20);
   //
   config.add<double>("max-horizontal-speed", "Maximum horizontal speed of the vehicle [m/s]", 5);
   config.add<double>("max-horizontal-accel", "Maximum horizontal acceleration of the vehicle [m/s^2]", 2);
   //
   config.add<double>("max-vertical-speed", "Maximum vertical speed of the vehicle [m/s]", 2);
   config.add<double>("max-vertical-accel", "Maximum vertical acceleration of the vehicle [m/s^2]", 1);
   //
   config.add<double>("sampling-step", "Time step for trajectory sampling [s]", 0.2);
   return config;
}

/// - static method ------------------------------------------------------------
CGSOACETSPBezier::TParamChange CGSOACETSPBezier::getParamChange(const std::string &str)
{
   if (str == "somhom") {
      return SOMHOM_PARAM;
   }
   if (str == "icann") {
      return ICANN_PARAM;
   }
   ASSERT_ARGUMENT(false, "Unknown param change");
}

/// - static method ------------------------------------------------------------
std::string CGSOACETSPBezier::getParamChange(TParamChange param)
{
   switch (param) {
      case SOMHOM_PARAM:
         return "somhom";
      case ICANN_PARAM:
         return "icann";
   }
   ASSERT_ARGUMENT(false, "Unknown change parameter");
}

/// - static method ------------------------------------------------------------
CGSOACETSPBezier::TRingLength CGSOACETSPBezier::getRingLength(const std::string &str)
{
   if (str == "ring_euclid") {
      return RING_LENGTH_EUCLID;
   }
   if (str == "ring_time") {
      return RING_LENGTH_TIME;
   }
   ASSERT_ARGUMENT(false, "Unknown ring length");
}

/// - static method ------------------------------------------------------------
std::string CGSOACETSPBezier::getRingLength(TRingLength ringLength)
{
   switch (ringLength) {
      case RING_LENGTH_EUCLID:
         return "ring_euclid";
         break;
      case RING_LENGTH_TIME:
         return "ring_time";
         break;
   }
}

/// - constructor --------------------------------------------------------------
CGSOACETSPBezier::CGSOACETSPBezier(crl::CConfig &config, const std::string &problemFile)
   : Base(config),
   PARAM_CHANGE(getParamChange(config.get<std::string>("param-change"))),
   RING_LENGTH(getRingLength(config.get<std::string>("ring-length"))),
   GAIN_DECREASING_RATE(config.get<double>("gain-decreasing-rate")),
   NEIGHBORHOOD_FACTOR(config.get<double>("neighborhood-factor")),
   NEURON_NUMBER_HEADINGS(config.get<int>("neuron-number-headings")),
   ALTERNATE_TARGET(config.get<bool>("alternate-target")),
   ICANN_GAIN_DECREASING_RATE(config.get<double>("icann-gain-decreasing-rate")),
   RING_REGENERATE_MIN_WINNER_DISTANCE(config.get<bool>("ring-regenerate-min-distance")),
   RING_REGENERATE_PRESERVE_ACTIVE_NEURONS(config.get<bool>("ring-regenerate-preserve-active")),
   FORCE_ZERO_ERROR_QUIT(config.get<bool>("force-zero-error-quit")),
   COVER_DISTANCE(config.get<double>("target-cover-distance")),
   COVER_DISTANCE_SQ(COVER_DISTANCE * COVER_DISTANCE),
   SAVE_RESULTS(config.get<bool>("save-results")),
   SAVE_SETTINGS(config.get<bool>("save-settings")),
   BORDER(config.get<double>("canvas-border")),
   SAVE_INFO(config.get<bool>("save-info")),
   NRINGS(config.get<int>("vehicles")),
   BEZIER_OPTIMIZATION_STEP(config.get<double>("bezier-optimization-step")),
   BEZIER_OPTIMIZATION_MAX(config.get<int>("bezier-optimization-max")),
   BEZIER_OPTIMIZATION_CYCLES(config.get<int>("bezier-optimization-cycles")),
   MAX_HORIZONTAL_SPEED(config.get<double>("max-horizontal-speed")),
   MAX_HORIZONTAL_ACCEL(config.get<double>("max-horizontal-accel")),
   MAX_VERTICAL_SPEED(config.get<double>("max-vertical-speed")),
   MAX_VERTICAL_ACCEL(config.get<double>("max-vertical-accel")),
   SAMPLING_STEP(config.get<double>("sampling-step"))
{
   // the prefered radius is calculated based on the speed and acceleration
   // that is the reason why this constant is not const
   RADIUS =
      (MAX_HORIZONTAL_SPEED * MAX_HORIZONTAL_SPEED / MAX_HORIZONTAL_ACCEL);

   bezier::BEZIER_SAMPLES = config.get<int>("bezier-samples");
   bezier::BEZIER_SAMPLES_VISUAL = config.get<int>("bezier-samples-visual");

   savePicEnabled = config.get<bool>("save-pic");

   drawNeurons = config.get<bool>("draw-neurons");
   drawAlternateGoals = config.get<bool>("draw-alternate-goals");
   drawWinners = config.get<bool>("draw-winners");
   shapeRing.setShape(config.get<std::string>("draw-shape-ring"));
   shapeNeuron.setShape(config.get<std::string>("draw-shape-neuron"));
   shapeWinner.setShape(config.get<std::string>("draw-shape-winner"));
   drawRingIter = config.get<bool>("draw-ring-iter");
   drawRingEnabled = config.get<bool>("draw-ring");
   drawTourRepresentedByRing =
      config.get<bool>("draw-tour-represented-by-ring");
   shapeAlternateGoal.setShape(
         config.get<std::string>("draw-shape-alternate-goal"));
   shapeTourRepresentedByRing.setShape(
         config.get<std::string>("draw-shape-tour-represented-by-ring"));

   velocityProfile.horizontalSpeed = MAX_HORIZONTAL_SPEED;
   velocityProfile.horizontalAccel = MAX_HORIZONTAL_ACCEL;
   velocityProfile.verticalSpeed = MAX_VERTICAL_SPEED;
   velocityProfile.verticalAccel = MAX_VERTICAL_ACCEL;

   load_targets_coords(problemFile);
   ASSERT_ARGUMENT(!targets.empty(), "Empty targets loaded");

   rings = 0;
   ringLengths.resize(NRINGS, 0.0);

   ASSERT_ARGUMENT(NRINGS > 0, "At least one vehicle must be set");
   // parse depots
   crl::parse_list(config.get<std::string>("start-targets"), ",", depots);
   ASSERT_ARGUMENT(depots.size() == NRINGS, "No. of lables in the "
         "start-targets does not much with "
         "the number of vehicles");

   foreach (int depot, depots) {
      ASSERT_ARGUMENT(depot >= 0 and depot < targets.size(),
            "Depot label out of range of the current targets");
      targets[depot]->depot = true;
   }
}

/// - destructor --------------------------------------------------------------
CGSOACETSPBezier::~CGSOACETSPBezier()
{
   deallocate_neurons();
   releaseTargets();
}

/// ----------------------------------------------------------------------------
std::string CGSOACETSPBezier::getRevision(void) { return "$Id$"; }

/// ----------------------------------------------------------------------------
void CGSOACETSPBezier::setTimers(const crl::CTimerN &load, const crl::CTimerN &init)
{
   loadTimer = load;
   initTimer = init;
}

/// ----------------------------------------------------------------------------
std::string CGSOACETSPBezier::getVersion(void) { return "CGSOACETSPBezier Bezier 0.1"; }

/// ----------------------------------------------------------------------------
void CGSOACETSPBezier::solve(void)
{
   crl::CRandom::randomize();
   Base::solve();
}

/// - protected method ---------------------------------------------------------
void CGSOACETSPBezier::load(void)
{
   // nothing to load, structures are passed to the constructor
   n = targets.size();
   if (canvas) {
      *canvas << canvas::AREA;
      foreach(const STarget *target, targets) {
         *canvas << (target->coords + Coords(1, 1, 0) * BORDER);
         *canvas << (target->coords - Coords(1, 1, 0) * BORDER);
      }
      *canvas << canvas::END;
      if (config.get<bool>("draw-targets")) {
         CShape shapeTarget(config.get<std::string>("draw-shape-target"));
         *canvas << "targets" << canvas::POINT;
         if (config.get<bool>("draw-depots")) {
            CShape shapeDepot(config.get<std::string>("draw-shape-depot"));

            foreach (const STarget *target, targets) {
               *canvas << (target->depot ? shapeDepot : shapeTarget)
                  << target->coords << canvas::END;
            }
         } else {
            *canvas << shapeTarget;

            foreach (const STarget *target, targets) {
               *canvas << target->coords << canvas::END;
            }
         }
      }

      if (COVER_DISTANCE > 0.0 and config.get<bool>("draw-cover-circle")) {
         *canvas << canvas::ARC
            << CShape(
                  config.get<std::string>("draw-shape-cover-circle"));
         foreach (const STarget *target, targets) {
            *canvas << target->coords << COVER_DISTANCE << 0.0 << 2 * M_PI;
         }
      }

      canvas->redraw();
   } // end canvas
}

/// - protected method ---------------------------------------------------------
void CGSOACETSPBezier::initialize(void)
{
   createPermutation(targets.size(), permutation); // create permutation
}

/// - protected method ---------------------------------------------------------
void CGSOACETSPBezier::iterate(int iter)
{
   allocate_neurons();

   switch (PARAM_CHANGE) {
      case SOMHOM_PARAM:
         G = 12.41 * n + 0.06;
         break;
      case ICANN_PARAM:
         G = 10;
         break;
         //     G = n / 3; //Stepan Dvorak
   }
   mi = config.get<double>("learning-rate");
   // thresholds for the termination conditions
   const double MAX_ERROR = config.get<double>("termination-error");
   const int MAX_STEPS = config.get<int>("termination-max-steps");
   if (canvas) {
      *canvas << canvas::CLEAR << "path"
         << "path";
   }
   double error = 2 * MAX_ERROR;
   int step = 0;
   // this initialization must be done in each iteration
   initialize_neurons(); // place neurons around first target
   int bestPathStep = -1;
   DoubleVector path_time(NRINGS);
   DoubleVector bestPathTime(NRINGS, std::numeric_limits<double>::max());
   double bestSolutionError;
   double e;
   StateVectorVector bestPaths;
   int stepWithZeroError = 0;
   StateVectorVector solutions;
   while (!((error < MAX_ERROR)) &&
         (step < MAX_STEPS)) { // perform adaptation step
      error = refine(step, error);
      ring_regenerate(step);
      INFO("error: " << error << " COVER_DISTANCE: " << COVER_DISTANCE);
      const bool admissibleSolution =
         error < COVER_DISTANCE and
         isSolutionAdmissibleCoverageSimple(step, path_time, e, solutions);
      if (canvas and drawRingIter) {
         drawRing(step);
         if (savePicEnabled) {
            savePic(step, true);
         }
      }
      if (admissibleSolution) {
         INFO("Step "
               << step
               << " current solution is admissible with the max length: "
               << getMax(path_time) << " current best solution: "
               << getMax(bestPathTime) << " e: " << e);
         if (getMax(path_time) < getMax(bestPathTime)) { // keep the best solution
            bestPathTime = path_time;
            bestPathStep = step;
            bestSolutionError = e;
            bestPaths = solutions;
         }
         if (FORCE_ZERO_ERROR_QUIT and e < 1e-6) { // almost zero
            stepWithZeroError += 1;
            if (getMax(path_time) > getMax(bestPathTime) and stepWithZeroError > 10) {
               DEBUG("Force zero error quit: " << stepWithZeroError);
               DEBUG("MAX_ERROR: " << MAX_ERROR << " e: " << e << " stepWithZeroError: " << stepWithZeroError);
               error = e; // force quit
            }
         }
      } else {
         INFO("Solution is not admissible\n");
         for (int r = 0; r < NRINGS; ++r) {
            if (path_time[r] > 0.0) {
               DEBUG("Step: " << step << " ring " << r << " current solution ( " << rings[r]->nbrNeurons << " neurons) is not admissible " << "length: " << path_time[r] << " error: " << e);
            } else {
               DEBUG("Step: " << step << " ring " << r << " current solution ( " << rings[r]->nbrNeurons << " neurons) is not admissible");
            }
         } // end all rings
      }
      switch (PARAM_CHANGE) {
         case SOMHOM_PARAM:
            G = G * (1 - GAIN_DECREASING_RATE);
            break;
         case ICANN_PARAM:
            G = G * (1 - ICANN_GAIN_DECREASING_RATE * (step + 1));
            break;
      }
      step += 1;
      for (int r = 0; r < NRINGS; ++r) {
         DEBUG("No. of neurons ring[" << r << "]: " << rings[r]->nbrNeurons);
      }
   } // end step loop
   INFO("Adaptation terminated error: " << error << " step: " << step);
   tSolve.stop();
   INFO("Best path with the length: " << getMax(bestPathTime) << " found in: "
         << bestPathStep << " steps");
   DoubleVector pathLength;
   if (config.get<bool>("best-path")) {
      pathLength = bestPathTime;
      finalPaths = bestPaths;
   } else {
      bestPathStep = -1;
   }
   double maxPathLength = getMax(pathLength);
   DoubleVector times;

   double maxTime = 0.0;
   for (int r = 0; r < NRINGS; ++r) {
      StateVector tour = finalPaths[r];
      tour.push_back(tour.front());
      times.push_back(velocityProfile.interpolateUniformTime(SAMPLING_STEP, velocityProfile.interpolateOpenLoop(tour)).back().time);
      if (maxTime < times.back()) {
         maxTime = times.back();
      }
   }
   fillResultRecord(iter, maxPathLength, maxTime, step); // default log results (iter, length, steps)
   resultLog << (bestPathStep == -1 ? "0" : "1") << bestPathStep << bestSolutionError << targets.size();
   resultLog << getString(pathLength, NRINGS) << getString(times, NRINGS) << crl::result::endrec;
   DEBUG("number targets: " << targets.size() << " maxpath: " << maxPathLength);
}

/// - protected method ---------------------------------------------------------

void CGSOACETSPBezier::save(void)
{
   std::string dir;
   updateResultRecordTimes(); // update timers as load and initilization is
   // outside class
   DEBUG("LOAD_TIME_CPU: " << tLoad.cpuTime());
   DEBUG("INIT_TIME_CPU: " << tInit.cpuTime());
   DEBUG("SAVE_TIME_CPU: " << tSave.cpuTime());
   DEBUG("SOLVE_TIME_CPU: " << tSolve.cpuTime());
   if (SAVE_SETTINGS) {
      saveSettings(getOutputIterPath(config.get<std::string>("settings"), dir));
   }
   if (SAVE_INFO) {
      saveInfo(getOutputIterPath(config.get<std::string>("info"), dir));
   }
   if (SAVE_RESULTS) {
      std::string file = getOutputIterPath(config.get<std::string>("result-path"), dir);
      crl::assert_io(createDirectory(dir), "Can not create file in path'" + file + "'");
      if (finalPaths.size() == NRINGS) {
         for (int r = 0; r < NRINGS; ++r) {
            const StateVector &finalPath = finalPaths[r];
            std::stringstream ss; // Files with Bezier maneuvers
            ss << file << std::setw(3) << std::setfill('0') << r << config.get<std::string>("result-ext");
            std::ofstream ofs(ss.str());
            crl::assert_io(not ofs.fail(), "Cannot create path '" + ss.str() + "'");
            ofs << std::setprecision(14);
            foreach (const State &pt, finalPath) {
               ofs << pt.coords << " " << pt.direction.vector << std::endl;
            }
            assert_io(not ofs.fail(), "Error occur during path saving");
            ofs.close();
            // Files with interpolated points
            std::stringstream ss2;
            ss2 
               << file << std::setw(3) << std::setfill('0') << r
               << "-samples" << config.get<std::string>("result-ext");
            std::string filename2 = ss2.str();
            std::ofstream ofs2(filename2.c_str());
            crl::assert_io(not ofs2.fail(),
                  "Cannot create path '" + filename2 + "'");
            ofs2 << std::setprecision(13);
            // create close loop
            StateVector tour = finalPath;
            tour.push_back(finalPath.front());
            bezier::InterpolationPointVector interpolation = velocityProfile.interpolateOpenLoop(tour);
            foreach(const bezier::InterpolationPoint &p, interpolation) {
               ofs2 
                  << p.coords.x << " " << p.coords.y << " " << p.coords.z
                  << " " << p.speed << " " << p.maxSpeed << " " << p.time
                  << " " << p.curvature << " " << p.maxAccel
                  << std::endl;
            }
            assert_io(not ofs2.fail(), "Error occur during path saving");
            ofs2.close();
            // Files with interpolated points
            std::stringstream ss3;
            ss3 << file << std::setw(3) << std::setfill('0') << r << "-samples-time" << config.get<std::string>("result-ext");
            std::string filename3 = ss3.str();
            std::ofstream ofs3(filename3.c_str());
            crl::assert_io(not ofs3.fail(), "Cannot create path '" + filename3 + "'");
            ofs3 << std::setprecision(13);
            const bezier::InterpolationPointVector interpolationUniform =  velocityProfile.interpolateUniformTime(SAMPLING_STEP, interpolation);
            foreach(const bezier::InterpolationPoint &p, interpolationUniform) {
               ofs3 
                  << p.coords.x << " " << p.coords.y << " " << p.coords.z
                  << " " << p.speed << " " << p.maxSpeed << " " << p.time
                  << " " << p.curvature << " " << p.maxAccel
                  << std::endl;
            }
            assert_io(not ofs3.fail(), "Error occur during path saving");
            ofs3.close();
         } // end all rings
      } else {
         DEBUG("Path not found, thus it cannot be saved");
      }
   }
   if (canvas) { // map must be set
      *canvas << canvas::CLEAR << "ring";
      if (config.get<bool>("draw-path")) {
         drawPath();
      } else if (drawRingEnabled) {
         drawRing(-1);
      }
      saveCanvas();
   }
}

/// - protected method ---------------------------------------------------------
void CGSOACETSPBezier::release(void) {}

/// - protected method ---------------------------------------------------------
void CGSOACETSPBezier::visualize(void) {}

/// - protected method ---------------------------------------------------------
void CGSOACETSPBezier::defineResultLog(void)
{
   static bool resultLogInitialized = false;
   if (!resultLogInitialized) {
      resultLog << result::newcol << "NAME";
      resultLog << result::newcol << "METHOD";
      resultLog << result::newcol << "TSP";
      resultLog << result::newcol << "ITER";
      resultLog << result::newcol << "RTIME";
      resultLog << result::newcol << "CTIME";
      resultLog << result::newcol << "UTIME";
      resultLog << result::newcol << "MAX_LENGTH";
      resultLog << result::newcol << "MAX_TIME";
      resultLog << result::newcol << "STEPS";
      resultLog << result::newcol << "ADMISSIBLE";
      resultLog << result::newcol << "SOLUTION_STEP";
      resultLog << CResultLog::AddColumn("ERROR", "x");
      resultLog << CResultLog::AddColumn("NBR_TARGETS", "x");
      resultLog << result::newcol << "LENGTHS";
      resultLog << result::newcol << "TIMES";
      resultLogInitialized = true;
   }
}

/// - protected method ---------------------------------------------------------
void CGSOACETSPBezier::fillResultRecord(int iter, double length, double time,
      int steps)
{
   resultLog << result::newrec << name << getMethod() << NRINGS << iter;
   long t[3] = {0l, 0l, 0l};
   tLoad.addTime(t);
   tInit.addTime(t);
   tSolve.addTime(t);
   tSave.addTime(t);
   resultLog << t[0] << t[1] << t[2] << length << time << steps;
}

/// - protected method ---------------------------------------------------------
void CGSOACETSPBezier::allocate_neurons(void)
{ // allocate neurons
   const double neuronMul =
      config.get<double>("number-neurons-multiplication");
   const int m = (int)round((targets.size() * neuronMul) / NRINGS);
   deallocate_neurons();
   rings = new SRing *[NRINGS];
   for (int r = 0; r < NRINGS; ++r) {
      rings[r] = new SRing();
      rings[r]->allocate_neurons(m);
   }
}

/// - protected method ---------------------------------------------------------
void CGSOACETSPBezier::deallocate_neurons(void)
{ // release neurons
   if (rings) {
      for (int i = 0; i < NRINGS; ++i) {
         delete rings[i];
      }
      delete[] rings;
      rings = 0;
   }
}

/// - private method -----------------------------------------------------------
void CGSOACETSPBezier::initialize_neurons(void)
{
   if (!targets.empty()) { // place neurons in circle around the problem center
      crl::CBBox bbox;
      double minz = std::numeric_limits<double>::max();
      double maxz = std::numeric_limits<double>::min();
      foreach (const STarget *t, targets) {
         bbox.add(t->coords.x, t->coords.y);
         minz = std::fmin(minz, t->coords.z);
         maxz = std::fmax(maxz, t->coords.z);
      }
      const Coords center(bbox.minx() + (bbox.maxx() - bbox.minx()) / 2.0, bbox.miny() + (bbox.maxy() - bbox.miny()) / 2.0, (minz + maxz) / 2.0);
      const double INITIAL_RADIUS = (bbox.maxx() + bbox.maxy() - bbox.minx() - bbox.miny()) / 2;
      avgLength = 0.0;
      for (int r = 0; r < NRINGS; ++r) {
         SRing *ring = rings[r];
         SNeuron *n = ring->neurons;
         const double pom = 2 * M_PI / ring->nbrNeurons;
         const double ppom = 2 * M_PI / NRINGS;
         for (int i = 0; i < ring->nbrNeurons; ++i) {
            double a = crl::CRandom::random(-1, 1);
            double b = crl::CRandom::random(-1, 1);
            Coords position(center.x + INITIAL_RADIUS * (a / 2), center.y + INITIAL_RADIUS * (b / 2), center.z);
            n->state.init(position, 0, RADIUS, RADIUS);
            n->targetOnTourStep = -1;
            n = n->next;
         } // end initial number of neurons
         ringLengths[r] = getRingLength(ring);
         DEBUG("Ring " << r << " ring length: " << ringLengths[r]);
         avgLength += ringLengths[r];
      } // end all rings
      avgLength /= NRINGS;
      DEBUG("Avg ring length: " << avgLength);
   }
}

/// - protected method ---------------------------------------------------------
void CGSOACETSPBezier::after_init(void)
{
   tLoad.append(loadTimer);
   tInit.append(initTimer);
}

/// - protected method ---------------------------------------------------------
double CGSOACETSPBezier::refine(int step, double errorMax)
{

   double errorToTarget = errorMax;
   double error = 0.0;
   SWinner winner;

   // 1st adapt individual rings to
   for (int r = 0; r < NRINGS; ++r) {
      STarget *target = targets[depots[r]];
      SNeuron *neuron = selectWinner(step, target, false, r, errorToTarget, winner);
      if (neuron) { // winner neuron has been selected
         target->stepWinnerSelected = step;
         adapt(winner.ring, neuron, winner, step);
         { // recompute avgLength
            avgLength = 0.0;
            ringLengths[winner.ringIDX] = getRingLength(winner.ring);
            for (int r = 0; r < NRINGS; ++r) {
               avgLength += ringLengths[r];
            }
            avgLength /= NRINGS;
         }
         neuron->activityStep = step; // force to preserve first neuron
      }                                // end adaptation
   }                                    // end adapt to depots
   permute(permutation);
   // 2nd select winner neuron to each target
   foreach(int i, permutation) {
      STarget *target = targets[i];
      if (target->depot) {
         continue;
      } // skip depots
      // select winner from all rings starting from the first ring with the
      // index 0
      if (SNeuron *neuron =
            selectWinner(step, target, true, 0, errorToTarget,
               winner)) { // winner neuron has been selected
         target->stepWinnerSelected = step;
         adapt(winner.ring, neuron, winner, step);
         { // recompute avgLength
            avgLength = 0.0;
            ringLengths[winner.ringIDX] = getRingLength(winner.ring);
            for (int r = 0; r < NRINGS; ++r) {
               avgLength += ringLengths[r];
            }
            avgLength /= NRINGS;
         }
         if (error < errorToTarget) {
            error = errorToTarget; // update error
         }
      }         // end adaptation
   }             // end permutation of all targets
   return error; // return largest error to city
}

/// - private method -----------------------------------------------------------
SNeuron *CGSOACETSPBezier::selectWinner(int step, STarget *target, bool ALL_RINGS,
      int RING_IDX, double &errorToTarget,
      SWinner &winner)
{
   winner.winner = 0;
   winner.prev = 0;
   winner.next = 0;
   winner.ring = 0;
   winner.ringIDX = -1;
   Coords alternateTarget;
   bool hasAlteranteTarget = false;
   bool hasWinner = false;
   double bestCost = std::numeric_limits<double>::max();
   double bestError = bestCost;

   // Closest neuron based on its euclidean distance is selected
   SNeuron *winnerNeuron = 0;
   State newState;
   for (int r = RING_IDX; r < NRINGS; ++r) {
      const double RING_WEIGHT = getRingWinnerPenalty(r);
      SRing *ring = rings[r];
      // 1. determine the closet point to the target
      SNeuron *cur = ring->neurons;
      for (int i = 0; i < ring->nbrNeurons; ++i) {
         cur->maneuver = bezier::Maneuver(cur->state, cur->next->state);
         // find the closes state on the maneuver  (euclidean distance)
         State st = cur->maneuver.getClosestState(target->coords);
         const double dist = st.coords.squared_distance(target->coords);
         const double d = RING_WEIGHT * sqrt(dist);
         if (d < bestCost) {
            bestCost = d;
            bestError = dist;
            winnerNeuron = cur;
            winner.ring = ring;
            winner.ringIDX = r;
            newState = st;
            hasWinner = true;
         }
         cur = cur->next;
      } // end all neurons in the ring = rings[r]
      if (not ALL_RINGS) {
         break;
      }
   } // end all rings
   bool covered = false;
   if ((hasWinner and bestError < COVER_DISTANCE_SQ) || !ALL_RINGS) {
      if (!winnerNeuron->isInhibited(step)) {
         winner.winner = winnerNeuron;
         covered = false;
      }
   }
   if (not covered) { // 2. determine the neighboring nodes
      // hasWinner not been found
      ASSERT_ARGUMENT(winner.ring, "Ring must be determined");
      if (winnerNeuron) {
         SNeuron *newNeuron = new SNeuron();
         winner.ring->insertNeuron(winnerNeuron, newNeuron);
         newNeuron->state = newState;
         determineCenterHeading(newNeuron->state.coords, newNeuron->next->state, newNeuron->next->next->state.coords);
         determineCenterHeading(newNeuron->prev->state.coords, newNeuron->state, newNeuron->next->state.coords);
         determineCenterHeading(newNeuron->prev->prev->state.coords, newNeuron->prev->state, newNeuron->state.coords);
         if (ALTERNATE_TARGET) {
            const Coords &pt = newNeuron->state.coords;
            const double d = sqrt(pt.squared_distance(target->coords));
            const double d2 = d - COVER_DISTANCE * 0.9;
            if (d2 > 0) {
               alternateTarget = pt + (target->coords - pt) * (d2 / d);
               hasAlteranteTarget = true;
            }
         }
         winner.winner = newNeuron;
         winner.prev = newNeuron->prev;
         winner.next = newNeuron->next;
      }
   }                                  // end determine the neighboring nodes
   if (hasWinner and winner.winner) { // winner has been selected
      if (ALTERNATE_TARGET and hasAlteranteTarget) {
         winner.winner->alternateGoal = alternateTarget;
      } else {
         winner.winner->alternateGoal = target->coords;
      }
      winner.winner->targetOnTour = target->label;
      winner.winner->targetOnTourStep = step;
      errorToTarget = sqrt(bestError);
      if (target->depot) {
         winner.winner->depot = true;
      }
   }
   return winner.winner;
}

/// - private method -----------------------------------------------------------
void CGSOACETSPBezier::adapt(SRing *ring, SNeuron *neuron, const SWinner &winner,
      const int step)
{
   SNeuron *neuronBackward = neuron->prev;
   SNeuron *neuronForward = neuron->next;
   double dd = 1.0;
   const int d = ring->nbrNeurons / NEIGHBORHOOD_FACTOR;
   const int TO = d < (ring->nbrNeurons / 2) ? d : ring->nbrNeurons / 2;
   const Coords target = neuron->alternateGoal;
   // Bezier adapt is simple euclidean adapt (tangents for Bezier curves are
   // not updated)
   if (winner.winner) {
      winner.winner->adapt(target, mi);
      bool prevReached = winner.prev ? false : true; // only if it has been determined
      bool nextReached = winner.next ? false : true; // only if it has been determined
      NeuronPtrVector prevNeurons;
      NeuronPtrVector nextNeurons;
      for (int i = 0; i < TO; ++i) {
         // adapt neighborhood according to the ETSP but adjust headings
         // for neurons between winner.prev and winner.next
         const double b = mi * exp((-dd * dd) / (G * G));
         if (neuronBackward and neuronBackward != neuron) {
            neuronBackward->adapt(target, b);
            if (not prevReached) {
               if (neuronBackward == winner.prev) {
                  prevReached = true;
               } else {
                  prevNeurons.push_back(neuronBackward);
               }
            }
            neuronBackward = neuronBackward->prev;
         }
         if (neuronForward and neuronForward != neuron) {
            neuronForward->adapt(target, b);
            if (not nextReached) {
               if (neuronForward == winner.next) {
                  nextReached = true;
               } else {
                  nextNeurons.push_back(neuronForward);
               }
            }
            neuronForward = neuronForward->next;
         }
         dd += 1.0;
      } // end update neighboring nodes
   }
}

/// - private method -----------------------------------------------------------
void CGSOACETSPBezier::drawRing(int step)
{
   if (canvas and rings) {
      *canvas << canvas::CLEAR << "ring"
         << "ring" << canvas::LINESTRING << shapeRing;
      if (drawNeurons) {
         *canvas << canvas::POINT << shapeNeuron;
         for (int r = 0; r < NRINGS; ++r) {
            const SRing *const ring = rings[r];
            const SNeuron *neuron = ring->neurons;
            for (int i = 0; i < ring->nbrNeurons; ++i) {
               if (drawWinners and neuron->targetOnTourStep == step) {
                  *canvas << shapeWinner << neuron->state.coords << shapeNeuron;
               } else {
                  *canvas << neuron->state.coords;
               }
               neuron = neuron->next;
            } // end all neurons
         }     // end rings
      } else if (drawWinners) {
         *canvas << canvas::POINT << shapeWinner;
         for (int r = 0; r < NRINGS; ++r) {
            const SRing *const ring = rings[r];
            const SNeuron *neuron = ring->neurons;
            for (int i = 0; i < ring->nbrNeurons; ++i) {
               if (drawWinners and neuron->targetOnTourStep == step) {
                  *canvas << shapeWinner << neuron->state.coords;
               } // end winner
               neuron = neuron->next;
            } // end all neurons
         }     // end rings
      }
      if (drawTourRepresentedByRing) {
         for (int r = 0; r < NRINGS; ++r) {
            const SRing *const ring = rings[r];
            CoordsVector path;
            get_ring_path(step, ring, path);
            if (path.size() > 1) {
               *canvas << shapeTourRepresentedByRing << canvas::LINESTRING << path;
               *canvas << path.front();
               *canvas << canvas::END;
            }
         } // end rings
      }
      for (int r = 0; r < NRINGS; ++r) {
         const SRing *const ring = rings[r];
         SNeuron *neuron = ring->neurons;
         // DEBUG("Draw ring as a sequence of Dubins maneuvers");
         *canvas << CShape(config.get<std::string>("draw-shape-ring"));
         for (int i = 0; i < ring->nbrNeurons; ++i) {
            *canvas << neuron->maneuver;
            neuron = neuron->next;
         }
      } // end rings
   }
}

/// - private method -----------------------------------------------------------
void CGSOACETSPBezier::drawPath(void)
{
   if (canvas) {
      *canvas << canvas::CLEAR << "path"
         << "path" << CShape(config.get<std::string>("draw-shape-path"))
         << Fill(false);
      foreach (const StateVector &finalPath, finalPaths) {
         for (int i = 0; i < finalPath.size(); ++i) {
            const State &st1 = finalPath[i];
            const State &st2 = finalPath[(i + 1) % finalPath.size()];
            *canvas << bezier::Maneuver(st1, st2);
         }
      } // end all vehicles
      if (config.get<bool>("draw-path-vertices")) {
         foreach (const StateVector &finalPath, finalPaths) {
            *canvas << canvas::POINT << Shape::MAP_VERTEX << finalPath;
         } // end all vehicles
      }     // end draw-path-vertices
   }         // end if canvas
}

/// - private method -----------------------------------------------------------
void CGSOACETSPBezier::savePic(int step, bool detail, const std::string &dir_suffix)
{
   static int lastStep = step;
   static int i = 0;
   if (lastStep != step) {
      i = 0;
   }
   if (canvas) {
      canvas->redraw();
      std::string dir;
      std::string file = getOutputIterPath(config.get<std::string>("pic-dir") + dir_suffix, dir);
      crl::assert_io(createDirectory(file), "Can not create file in path '" + file + "'");
      std::stringstream ss;
      ss << file << "/" << "iter-" << std::setw(3) << std::setfill('0') << step;
      ss << "-" << std::setw(4) << std::setfill('0') << i;
      std::string suffixes(config.get<std::string>("pic-ext"));
      if (!suffixes.empty()) {
         std::string::size_type cur = 0;
         std::string::size_type next;
         do {
            next = suffixes.find(',', cur);
            const std::string &ext = suffixes.substr(cur, next - cur);
            if (!ext.empty()) {
               crl::assert_io(canvas->save(ss.str() + "." + ext), "Can create output canvas file '" + file + "'");
            }
            cur = next + 1;
         } while (next != std::string::npos);
      } else {
         ss << "." << config.get<std::string>("pic-ext");
         crl::assert_io(canvas->save(ss.str()), "Can create output canvas file '" + ss.str() + "'");
      }
   }
   lastStep = step;
   i += 1;
}

/// - private method -----------------------------------------------------------
void CGSOACETSPBezier::releaseTargets(void)
{
   foreach (STarget *st, targets) {
      delete st;
   }
   targets.clear();
}

/// - private method -----------------------------------------------------------
void CGSOACETSPBezier::ring_regenerate(int step)
{
   for (int r = 0; r < NRINGS; ++r) {
      SRing *ring = rings[r];
      NeuronPtrVector del;
      NeuronPtrVector preserve; // TODO gui only
      SNeuron *cur = ring->neurons;
      for (int i = 0; i < ring->nbrNeurons; ++i) {
         if (cur->isInhibited(step) or
               (RING_REGENERATE_PRESERVE_ACTIVE_NEURONS and cur->activityStep >= step)) {
            if (not RING_REGENERATE_MIN_WINNER_DISTANCE or preserve.empty()) {
               preserve.push_back(cur);
            } else {
               if (preserve.back()->state.coords.squared_distance(cur->state.coords) > 1e-1 * 1e-1) {
                  preserve.push_back(cur);
               } else {
                  del.push_back(cur);
               }
            }
         } else {
            del.push_back(cur);
         }
         cur = cur->next;
      }
      foreach (SNeuron *n, del) { // delete all non winners and active neurons
         if (ring->nbrNeurons > 3) { // dont remove if there are only 3 neurons left
            ring->removeNeuron(n);
         }
      }
      updateHeadings(step, ring, true); // optimize the trajectory
      cur = ring->neurons;              // generate maneuvers
      for (int i = 0; i < ring->nbrNeurons; ++i) {
         if (cur->next) {
            cur->maneuver = bezier::Maneuver(cur->state, cur->next->state);
            cur = cur->next;
         }
      }
   } // end all rings
}

/// - private method -----------------------------------------------------------
void CGSOACETSPBezier::load_targets_coords(const std::string &filename)
{
   releaseTargets();
   Coords pt;
   std::ifstream in(filename.c_str());
   while (in >> pt) {
      STarget *t = new STarget(targets.size(), pt);
      t->ringDistance = COVER_DISTANCE;
      targets.push_back(t);
   }
   in.close();
}

/// - private method -----------------------------------------------------------
bool CGSOACETSPBezier::isSolutionAdmissibleCoverageSimple(
      int step, DoubleVector &times, double &solutionError,
      StateVectorVector &solutions)
{
   double maxMinError = 0.0;
   solutions.clear();
   times.clear();
   times.resize(NRINGS, 0.0);
   for (int r = 0; r < NRINGS; ++r) {
      const SRing *const ring = rings[r];
      SNeuron *cur = ring->neurons;
      for (int i = 0; i < ring->nbrNeurons; ++i) { // precompute all maneuvers
         cur->maneuver = bezier::Maneuver(cur->state, cur->next->state);
         cur = cur->next;
      }
   }
   // for each dubins maneuver get all targets that can be covered
   // order the states according to the distance than
   foreach (STarget *target, targets) {
      target->ringDistance = std::numeric_limits<double>::max();
   }
   std::vector<bool> covered(targets.size(), false);
   for (int r = 0; r < NRINGS; ++r) {
      const SRing *const ring = rings[r];
      SNeuron *cur = ring->neurons;
      for (int i = 0; i < ring->nbrNeurons; ++i) { //
         if (cur->isInhibited(step)) {
            foreach (const STarget *target, targets) {
               if (cur->state.coords.squared_distance(target->coords) < COVER_DISTANCE_SQ) {
                  covered[target->label] = true;
               }
            }
         } // end winner
         cur = cur->next;
      }
   } // end all rings
   // Dermine new winners for the targets not covered by winners by from
   // manuvers  Force solution to be admissible
   bool allCovered = true;
   foreach (const STarget *target, targets) {
      if (not covered[target->label]) {
         DEBUG("Target: " << target->label << " not covered distance to ring: " << target->ringDistance);
         allCovered = false;
         break;
      }
   } // end check all targets covered
   for (int r = 0; r < NRINGS; ++r) {
      const SRing *const ring = rings[r];
      solutions.push_back(StateVector());
      StateVector &solution = solutions.back();
      SNeuron *cur = ring->neurons;
      while (!cur->depot) {
         cur = cur->next;
      }
      for (int i = 0; i < ring->nbrNeurons; ++i) { // all winners must be evaluated
         if (cur->isInhibited(step)) {
            solution.push_back(cur->state);
         }
         cur = cur->next;
      } // end all neurons
   }    // end all vehicles (rings)
   times.clear();
   foreach (const StateVector &solution, solutions) {
      times.push_back(0.0);
      double &time = times.back();
      StateVector tour = solution;
      tour.push_back(solution.front());
      time = velocityProfile.interpolateOpenLoop(tour).back().time;
   } // end all solutions
   if (allCovered) {
      solutionError = sqrt(maxMinError);
   } else {
      solutions.clear();
   }
   DEBUG("isSolutionAdmissibleCoverage -- neurons:target all covered: " << allCovered);
   DEBUG("Step: " << step << " all covered: " << allCovered << " isSolutionAdmissibleCoverage max min error : " << maxMinError);
   return allCovered;
}

/// - private method -----------------------------------------------------------
CoordsVector &CGSOACETSPBezier::get_ring_path(int step, const SRing *const ring, CoordsVector &path) const
{
   path.clear();
   SNeuron *cur = ring->neurons;
   for (int i = 0; i < ring->nbrNeurons; ++i) {
      if (cur->isInhibited(step) and cur->targetOnTourStep >= 0) {
         path.push_back(targets[cur->targetOnTour]->coords);
      }
      cur = cur->next;
   }
   return path;
}

/// - private method -----------------------------------------------------------
void CGSOACETSPBezier::determineCenterHeading(Coords p1, State &state, Coords p3)
{
   const Coords &p2 = state.coords;
   const Coords d1 = (p2 - p1);
   const Coords d2 = (p3 - p2);
   double alpha = d1.length() / (d1.length() + d2.length());
   Coords d = d1 * (1 - alpha) + d2 * alpha;

   state.direction.vector = d.normalize();
   state.direction.b = std::min(RADIUS, d1.length() / 2.5);
   state.direction.a = std::min(RADIUS, d2.length() / 2.5);
}

/// - private method -----------------------------------------------------------
double CGSOACETSPBezier::calculateTime(const State &s1, const State &s2,
      const State &s3)
{
   StateVector path;
   path.push_back(s1);
   path.push_back(s2);
   path.push_back(s3);
   return velocityProfile.interpolateOpenLoop(path, MAX_HORIZONTAL_SPEED).back().time;
}

/// - private method -----------------------------------------------------------
void CGSOACETSPBezier::optimizeDirection(const State &s1, State &s2, const State &s3)
{
   determineCenterHeading(s1.coords, s2, s3.coords);
   const Coords &p1 = s1.coords;
   const Coords &p2 = s2.coords;
   const Coords &p3 = s3.coords;
   State bestState = s2;
   double bestTime = std::numeric_limits<double>::max();
   StateVector possibleStates;
   const double step = BEZIER_OPTIMIZATION_STEP;
   int i;
   for (i = 0; i < BEZIER_OPTIMIZATION_MAX; i++) {
      possibleStates.clear();

      State candidate = bestState;
      candidate.direction.a *= 1 - step;
      possibleStates.push_back(candidate);

      candidate = bestState;
      candidate.direction.a *= 1 + step;
      possibleStates.push_back(candidate);

      candidate = bestState;
      candidate.direction.a *= 1 - step;
      possibleStates.push_back(candidate);

      candidate = bestState;
      candidate.direction.a *= 1 + step;
      possibleStates.push_back(candidate);

      candidate = bestState;
      candidate.direction.vector.z += step;
      candidate.direction.vector = candidate.direction.vector.normalize();
      possibleStates.push_back(candidate);

      candidate = bestState;
      candidate.direction.vector.z -= step;
      candidate.direction.vector = candidate.direction.vector.normalize();
      possibleStates.push_back(candidate);

      candidate = bestState;
      candidate.direction.vector.z -= step;
      candidate.direction.vector = candidate.direction.vector.normalize();
      possibleStates.push_back(candidate);

      const Coords orig = bestState.direction.vector;
      double vsin = std::sin(step);
      double vcos = std::cos(step);

      candidate = bestState;
      candidate.direction.vector.x = vcos * orig.x - vsin * orig.y;
      candidate.direction.vector.y = vsin * orig.x + vcos * orig.y;
      possibleStates.push_back(candidate);

      candidate = bestState;
      candidate.direction.vector.x = vcos * orig.x + vsin * orig.y;
      candidate.direction.vector.y = -vsin * orig.x + vcos * orig.y;
      possibleStates.push_back(candidate);
      bool change = false;
      foreach (const State &st, possibleStates) {
         double time = calculateTime(candidate, st, s3);
         if (time < bestTime) {
            bestTime = time;
            bestState = st;
            change = true;
         }
      }
      if (!change) {
         break;
      }
   }
   if (i == BEZIER_OPTIMIZATION_MAX) {
      WARN("Maximum number of "
            << i << " reached during Bezier curve optimization");
   }
   s2 = bestState;
}

/// - private method -----------------------------------------------------------
void CGSOACETSPBezier::updateHeadings(int step, SRing *ring, bool fixWinners)
{
   SNeuron *cur = ring->neurons;
   for (int i = 0; i < ring->nbrNeurons * BEZIER_OPTIMIZATION_CYCLES; ++i) {
      SNeuron *prev = cur->prev;
      SNeuron *next = cur->next;
      if (prev && next) {
         const State &prev_state = prev->state;
         State &cur_state = cur->state;
         const State &next_state = next->state;
         optimizeDirection(prev_state, cur_state, next_state);
      } else {
         ERROR("Only close loop rings are supported");
         exit(42);
      }
      cur = cur->next; // iterate
   }
}

/// - private method -----------------------------------------------------------
double CGSOACETSPBezier::getRingWinnerPenalty(const int ringIDX)
{
   return (1 + (ringLengths[ringIDX] - avgLength) / avgLength);
}

/// - private method -----------------------------------------------------------
double CGSOACETSPBezier::getRingLength(const SRing *const ring) const
{
   switch (RING_LENGTH) {
      case RING_LENGTH_EUCLID:
         return getRingLengthEuclid(ring);
      case RING_LENGTH_TIME:
         return getRingLengthTime(ring);
   }
   ASSERT_ARGUMENT(false, "Unknown ring length method");
}

/// - private method -----------------------------------------------------------
double CGSOACETSPBezier::getRingLengthEuclid(const SRing *const ring) const
{
   double length = 0.0;
   SNeuron *cur = ring->neurons;
   for (int i = 0; i < ring->nbrNeurons; ++i) {
      length += bezier::Maneuver(cur->state, cur->next->state).getLength();
      cur = cur->next;
   }
   return length;
}

/// - private method -----------------------------------------------------------
double CGSOACETSPBezier::getRingLengthTime(const SRing *const ring) const
{
   StateVector tour;
   SNeuron *cur = ring->neurons;
   for (int i = 0; i <= ring->nbrNeurons; ++i) { // <= to create close loop
      tour.push_back(cur->state);
      cur = cur->next;
   }
   return velocityProfile.interpolateOpenLoop(tour).back().time;
}

/* end of gsoa_cetsp_bezier.cc */
