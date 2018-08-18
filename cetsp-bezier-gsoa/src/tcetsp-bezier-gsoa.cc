/*
 * File name: tbezier-cetsp-gsoa.cc
 * Date:      2016/02/05 08:54
 * Author:    Jan Faigl
 */

#include <fstream>
#include <limits.h>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>

#include <crl/boost_args_config.h>
#include <crl/config.h>

#include <crl/exceptions.h>
#include <crl/gui/guifactory.h>
#include <crl/gui/shapes.h>
#include <crl/gui/win_adjust_size.h>
#include <crl/logging.h>
#include <crl/perf_timer.h>
#include <crl/timerN.h>

#include "canvasview.h"

#include "gsoa_cetsp_bezier.h"

#define foreach(x, y) for (x : y)

using crl::logger;

namespace po = boost::program_options;
namespace fs = boost::filesystem;

const std::string CETSP_BEZIER_GSOA_VERSION = "1.0";

using namespace crl;
using namespace crl::gui;

typedef crl::gui::CCanvasBase Canvas;
typedef gsoa::CGSOACETSPBezier Solver;

typedef std::vector<Coords> CoordsVector;
typedef std::vector<CoordsVector> CoordsVectorVector;

/// ----------------------------------------------------------------------------
/// Program options variables
/// ----------------------------------------------------------------------------
std::string guiType = "none";

crl::CConfig guiConfig;
crl::CConfig solverConfig;
std::string canvasOutput = "";

/// ----------------------------------------------------------------------------
/// Global variable
/// ----------------------------------------------------------------------------
crl::gui::CGui *g = 0;
#define GUI(x)                                                                 \
   if (gui) {                                                                 \
      x;                                                                     \
   }

/// ----------------------------------------------------------------------------

bool parseArgs(int argc, char *argv[])
{
   bool ret = true;
   std::string configFile;
   std::string guiConfigFile;
   std::string loggerCfg = "";

   po::options_description desc("General options");
   desc.add_options()
      ("help,h", "produce help message")
      ("config,c", po::value<std::string>(&configFile)->default_value(std::string(argv[0]) + ".cfg"), "configuration file")
      ("logger-config,l", po::value<std::string>(&loggerCfg)->default_value(loggerCfg), "logger configuration file")
      ("config-gui", po::value<std::string>(&guiConfigFile)->default_value(std::string(argv[0]) + "-gui.cfg"), "dedicated gui configuration file");
   try {
      po::options_description guiOptions("Gui options");
      crl::gui::CGuiFactory::getConfig(guiConfig);
      crl::gui::CWinAdjustSize::getConfig(guiConfig);
      guiConfig.add<double>("gui-add-x", "add the given value to the loaded goals x coord to determine the canvas size and transformation", 0);
      guiConfig.add<double>("gui-add-y", "add the given value to the loaded goals y coord to determine the canvas size and transformation", 0);
      boost_args_add_options(guiConfig, "", guiOptions);
      guiOptions.add_options()("canvas-output", po::value<std::string>(&canvasOutput), "result canvas outputfile");

      po::options_description solverOptions("Solver options");
      solverConfig.add<std::string>("problem", "File with the targets to be visited by the vehicle", ""); 
      boost_args_add_options(Solver::getConfig(solverConfig), "", solverOptions);

      po::options_description cmdline_options;
      cmdline_options.add(desc).add(guiOptions).add(solverOptions);

      po::variables_map vm;
      po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
      po::notify(vm);

      std::ifstream ifs(configFile.c_str());
      store(parse_config_file(ifs, cmdline_options), vm);
      po::notify(vm);
      ifs.close();
      ifs.open(guiConfigFile.c_str());
      store(parse_config_file(ifs, cmdline_options), vm);
      po::notify(vm);
      ifs.close();

      if (vm.count("help")) {
         std::cerr << std::endl;
         std::cerr << "CETSP with Bezier curves solver based on GSOA ver. " << CETSP_BEZIER_GSOA_VERSION
            << std::endl;
         std::cerr << cmdline_options << std::endl;
         ret = false;
      }
      if (ret && loggerCfg != "" && fs::exists(fs::path(loggerCfg))) {
         crl::initLogger("bezier", loggerCfg.c_str());
      } else {
         crl::initLogger("bezier");
      }
      const std::string &problemFile = solverConfig.get<std::string>("problem");
      if (!fs::exists(fs::path(problemFile))) {
         ERROR("Problem file '" + problemFile + "' does not exists");
         ret = false;
      }
   } catch (std::exception &e) {
      std::cerr << std::endl;
      std::cerr << "Error in parsing arguments: " << e.what() << std::endl;
      ret = false;
   }
   return ret;
}

/// - local helper function ----------------------------------------------------
CoordsVector &load_goals_coords(const std::string &filename, CoordsVector &pts)
{
   Coords pt;
   std::ifstream in(filename.c_str());
   while (in >> pt.x >> pt.y) {
      pts.push_back(pt);
   }
   return pts;
}

/// - Main program -------------------------------------------------------------
int main(int argc, char *argv[])
{
   Canvas *canvas = 0;
   int ret = -1;
   if (parseArgs(argc, argv)) {
      INFO("Start Logging");
      try {
         crl::CPerfTimer t("Load problem time real:");
         CoordsVector pts;
         crl::CTimerN tLoad;
         tLoad.start();
         const std::string &problemFile = solverConfig.get<std::string>("problem");
         load_goals_coords(problemFile, pts);
         crl::gui::CWinAdjustSize::adjust(pts, guiConfig);
         if ((g = gui::CGuiFactory::createGui(guiConfig)) != 0) {
            INFO("Start gui " + guiConfig.get<std::string>("gui"));
            canvas = new Canvas(*g);
         }
         tLoad.stop();
         crl::CTimerN tInit;

         Solver solver(solverConfig, problemFile);
         solver.setTimers(tLoad, tInit);
         solver.setCanvas(canvas);
         t.start("TSP solve time: ");
         solver.solve();
         t.stop();
         INFO("End Logging");
         if (canvas) {
            if (canvasOutput.size()) {
               canvas->save(canvasOutput);
            }
            if (!guiConfig.get<bool>("nowait")) {
               INFO("click to exit");
               canvas->click();
            }
            delete canvas;
            delete g;
         }
      } catch (crl::exception &e) {
         ERROR("Exception " << e.what() << "!");
      } catch (std::exception &e) {
         ERROR("Runtime error " << e.what() << "!");
      }
      ret = EXIT_SUCCESS;
   }
   crl::shutdownLogger();
   return ret;
}

/* end of tbezier-cetsp-gsoa.cc */
